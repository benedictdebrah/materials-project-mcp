import os
import logging
from typing import Optional, List, Union
from mcp.server.fastmcp import FastMCP
from pydantic import Field
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.plotter import PhononBSPlotter
from pymatgen.analysis.wulff import WulffShape
from emmet.core.electronic_structure import BSPathType
from typing import Literal
import base64
import io 
from dotenv import load_dotenv

# Materials Project client
from mp_api.client import MPRester

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("materials_project_mcp")
load_dotenv()
API_KEY = os.environ.get("MP_API_KEY")

# Create the MCP server instance
mcp = FastMCP(
    name="Materials Project",
    version="0.0.1",
    description=(
        "A Model Context Protocol (MCP) server that exposes query tools "
        "for the Materials Project database using the mp_api client."
    ),
)


def _get_mp_rester() -> MPRester:
    """
    Initialize and return a MPRester session with the user's API key.
    
    Returns:
        MPRester: An authenticated MPRester instance for querying the Materials Project API.
        
    Note:
        If no API key is found in environment variables, attempts to initialize without key.
    """
    if not API_KEY:
        logger.warning(
            "No MP_API_KEY found in environment. Attempting MPRester() without key."
        )
        return MPRester()
    return MPRester(API_KEY)


@mcp.tool()
async def search_materials(
    elements: Optional[List[str]] = Field(
        default=None,
        description="List of element symbols to filter by (e.g. ['Si', 'O']). If None, searches across all elements.",
    ),
    band_gap_min: float = Field(
        default=0.0, 
        description="Lower bound for band gap filtering in eV. Materials with band gaps below this value will be excluded.",
    ),
    band_gap_max: float = Field(
        default=10.0, 
        description="Upper bound for band gap filtering in eV. Materials with band gaps above this value will be excluded.",
    ),
    is_stable: bool = Field(
        default=False,
        description="If True, only returns materials that are thermodynamically stable (energy above hull = 0). If False, returns all materials.",
    ),
    max_results: int = Field(
        default=50, 
        ge=1, 
        le=200, 
        description="Maximum number of results to return. Must be between 1 and 200.",
    ),
) -> str:
    """
    Search for materials in the Materials Project database using various filters.
    
    This function allows searching for materials based on their elemental composition,
    band gap range, and thermodynamic stability. Results are returned in a formatted
    markdown string containing material IDs, formulas, band gaps, and energy above hull values.
    
    Args:
        elements: Optional list of element symbols to filter by (e.g. ['Si', 'O'])
        band_gap_min: Minimum band gap in eV (default: 0.0)
        band_gap_max: Maximum band gap in eV (default: 10.0)
        is_stable: Whether to only return stable materials (default: False)
        max_results: Maximum number of results to return (default: 50, max: 200)
        
    Returns:
        str: A formatted markdown string containing the search results
        
    Example:
        >>> search_materials(elements=['Si', 'O'], band_gap_min=1.0, band_gap_max=5.0)
        Returns materials containing Si and O with band gaps between 1 and 5 eV
    """
    logger.info("Starting search_materials query...")
    with _get_mp_rester() as mpr:
        docs = mpr.materials.summary.search(
            elements=elements,
            band_gap=(band_gap_min, band_gap_max),
            is_stable=is_stable,
            fields=["material_id", "formula_pretty", "band_gap", "energy_above_hull"],
        )

    # Truncate results to max_results
    docs = list(docs)[:max_results]

    if not docs:
        return "No materials found matching your criteria."

    results_md = (
        f"## Materials Search Results\n\n"
        f"- **Elements**: {elements or 'Any'}\n"
        f"- **Band gap range**: {band_gap_min} eV to {band_gap_max} eV\n"
        f"- **Stable only**: {is_stable}\n\n"
        f"**Showing up to {max_results} matches**\n\n"
    )
    for i, mat in enumerate(docs, 1):
        results_md += (
            f"**{i}.** ID: `{mat.material_id}` | Formula: **{mat.formula_pretty}** | "
            f"Band gap: {mat.band_gap:.3f} eV | E above hull: {mat.energy_above_hull:.3f} eV\n"
        )
    return results_md


@mcp.tool()
async def get_structure_by_id(
    material_id: str = Field(
        ..., 
        description="Materials Project ID (e.g. 'mp-149'). Must be a valid MP ID."
    )
) -> str:
    """
    Retrieve and format the crystal structure for a given material from the Materials Project.
    
    This function fetches the final computed structure for a material and returns a
    formatted summary including the lattice parameters, number of sites, and chemical formula.
    
    Args:
        material_id: The Materials Project ID of the material (e.g. 'mp-149')
        
    Returns:
        str: A formatted markdown string containing the structure information
        
    Example:
        >>> get_structure_by_id('mp-149')
        Returns the crystal structure information for silicon (mp-149)
    """
    logger.info(f"Fetching structure for {material_id}...")
    with _get_mp_rester() as mpr:
        structure = mpr.get_structure_by_material_id(material_id)

    if not structure:
        return f"No structure found for {material_id}."

    formula = structure.composition.reduced_formula
    lattice = structure.lattice
    sites_count = len(structure)
    text_summary = (
        f"## Structure for {material_id}\n\n"
        f"- **Formula**: {formula}\n"
        f"- **Lattice**:\n"
        f"   a = {lattice.a:.3f} Å, b = {lattice.b:.3f} Å, c = {lattice.c:.3f} Å\n"
        f"   α = {lattice.alpha:.2f}°, β = {lattice.beta:.2f}°, γ = {lattice.gamma:.2f}°\n"
        f"- **Number of sites**: {sites_count}\n"
        f"- **Reduced formula**: {structure.composition.reduced_formula}\n"
    )
    return text_summary


@mcp.tool()
async def get_electronic_bandstructure(
    material_id: str = Field(
        ..., 
        description="Materials Project ID (e.g. 'mp-149'). Must be a valid MP ID."
    ),
    path_type: Literal["setyawan_curtarolo", "hinuma", "latimer_munro", "uniform"] = Field(
        default="setyawan_curtarolo",
        description="Type of k-point path to use for the band structure plot. Options are:\n"
                   "- setyawan_curtarolo: Standard path for cubic systems\n"
                   "- hinuma: Standard path for hexagonal systems\n"
                   "- latimer_munro: Alternative path for cubic systems\n"
                   "- uniform: Uniform k-point sampling (not recommended for plotting)"
    ),
) -> str:
    """
    Generate and return a electronic band structure plot for a given material.
    
    This function fetches the band structure data from the Materials Project and creates
    a plot showing the electronic band structure along high-symmetry k-points. The plot
    is returned as a base64-encoded PNG image embedded in a markdown string.
    
    Args:
        material_id: The Materials Project ID of the material (e.g. 'mp-149')
        path_type: The type of k-point path to use for the band structure plot
        
    Returns:
        A plot of the electronic band structure 
        
    Example:
        >>> get_electronic_bandstructure('mp-149', path_type='setyawan_curtarolo')
        Returns a band structure plot for silicon using the standard cubic path
    """
    logger.info(f"Plotting band structure for {material_id} with path_type: {path_type}")

    with _get_mp_rester() as mpr:
        if path_type == "uniform":
            bs = mpr.get_bandstructure_by_material_id(material_id, line_mode=False)
        else:
            bs = mpr.get_bandstructure_by_material_id(
                material_id, path_type=BSPathType(path_type)
            )

    if not isinstance(bs, BandStructureSymmLine):
        return f"Cannot plot `{path_type}` band structure. Only line-mode paths are plottable."

    plotter = BSPlotter(bs)
    ax = plotter.get_plot()
    fig = ax.get_figure()  
    fig.show()  


@mcp.tool()
async def get_electronic_dos_by_id(
    material_id: str = Field(
        ..., 
        description="Materials Project ID (e.g. 'mp-149'). Must be a valid MP ID."
    ),
) -> str:
    """
    Retrieve the electronic density of states (DOS) data for a given material.
    
    This function fetches the electronic density of states data from the Materials Project
    for the specified material. The DOS data includes information about the
    electronic states available to electrons in the material.
    
    Args:
        material_id: The Materials Project ID of the material (e.g. 'mp-149')
        
    Returns:
        str: A string containing the density of states information
        
    Example:
        >>> get_electronic_dos_by_id('mp-149')
        Returns the electronic density of states data for silicon
    """   
    logger.info(f"Fetching electronic density of states for {material_id}...")
    with _get_mp_rester() as mpr:
        dos = mpr.get_dos_by_material_id(material_id)

    if not dos:
        return f"No density of states found for {material_id}."
    
    return f"Electronic density of states for {material_id}: {dos}"

#phonons
@mcp.tool()
async def get_phonon_bandstructure(
    material_id: str = Field(
        ..., 
        description="Materials Project ID (e.g. 'mp-149'). Must be a valid MP ID."
    ),
) -> str:
    """
    Retrieve the phonon band structure for a given material.
    
    This function fetches the phonon band structure data from the Materials Project
    for the specified material. The phonon band structure includes information about
    the vibrational modes and frequencies of the material.
    
    Args:
        material_id: The Materials Project ID of the material (e.g. 'mp-149')
        
    Returns:
        A plot of the phonon band structure 
    """

    logger.info(f"Fetching phonon band structure for {material_id}...")
    with _get_mp_rester() as mpr:
        bs = mpr.get_phonon_bandstructure_by_material_id(material_id)

    if not isinstance(bs, PhononBandStructureSymmLine):
        return "Cannot plot phonon band structure. Only line-mode paths are plottable."    

    plotter = PhononBSPlotter(bs)
    fig = plotter.get_plot()
    
    plt.title(f"Phonon Band Structure for {material_id}")
    plt.ylabel("Frequency (THz)")
    plt.tight_layout()
    
    return fig  
 



@mcp.tool()
async def get_phonon_dos_by_id(
    material_id: str = Field(
        ..., 
        description="Materials Project ID (e.g. 'mp-149'). Must be a valid MP ID."
    ),
) -> str:
    """
    Retrieve the phonon density of states (DOS) data for a given material.
    
    This function fetches the phonon density of states data from the Materials Project
    for the specified material. The DOS data includes information about the
    vibrational modes and frequencies of the material.
    """
    logger.info(f"Fetching phonon density of states for {material_id}...")
    with _get_mp_rester() as mpr:
        dos = mpr.get_phonon_dos_by_material_id(material_id)

    if not dos:
        return f"No density of states found for {material_id}."

    return f"Phonon density of states for {material_id}: {dos}"
    

@mcp.tool()
async def get_ion_reference_data_for_chemsys(
    chemsys: Optional[Union[List, str]] = Field(
        ..., 
        description="Chemical system string comprising element symbols separated by dashes, e.g., 'Li-Fe-O' or List of element symbols, e.g., ['Li', 'Fe', 'O']"
    )
) -> str: 
    """
    Downloads aqueouse  ion reference data used in the contruction Pourbaix 
    The data returned from this method can be passed to get_ion_entries(). 

    Args:
        chemsys (str | list):  Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].

    Returns:
            str: markdown format of the reference data for ions 
    """

    logger.info("Fetch reference data for ion by Chemsys")
    mpr_rester = _get_mp_rester()

    with mpr_rester as mpr: 
        ion_reference_data = mpr.get_ion_reference_data_for_chemsys(chemsys=chemsys)

    if not ion_reference_data: 
        logger.info(f"data not found for {chemsys}")
        return f"No ion reference data for {chemsys}"


    ion_data = f"Ion Reference Data for Chemical System: {chemsys}\n\n"

    for idx, ion in enumerate(ion_reference_data, 1): 
        identifier = ion.get("identifier", "Unknown")
        formula = ion.get("formula", "Unknown")
        data = ion.get("data", {})

        # get the properties for idx 
        charge_info = data.get("charge", {})
        charge_value = charge_info.get('value', 0)
        charge_display = charge_info.get('display', str(charge_value))
            
        delta_gf_info = data.get('ΔGᶠ', {})
        delta_gf_value = delta_gf_info.get('value', 'N/A')
        delta_gf_display = delta_gf_info.get('display', f'{delta_gf_value} kJ/mol' if delta_gf_value != 'N/A' else 'N/A')

        maj_elements = data.get('MajElements', 'Unknown')
        ref_solid = data.get('RefSolid', 'Unknown')

        ref_solid_info = data.get('ΔGᶠRefSolid', {})
        ref_solid_value = ref_solid_info.get('value', 'N/A')
        ref_solid_display = ref_solid_info.get('display', f'{ref_solid_value} kJ/mol' if ref_solid_value != 'N/A' else 'N/A')

        reference = data.get('reference', 'No reference provided')

        ion_data +=  f"""## {idx}. {identifier}

            | Property | Value |
            |----------|--------|
            | *Formula* | {formula} |
            | *Charge* | {charge_display} |
            | *Formation Energy (ΔGᶠ)* | {delta_gf_display} |
            | *Major Elements* | {maj_elements} |
            | *Reference Solid* | {ref_solid} |
            | *Ref. Solid ΔGᶠ* | {ref_solid_display} |

            *Reference:* {reference}
        """

    return ion_data


@mcp.tool()
async def get_cohesive_energy(
        material_ids: List[str] = Field(
            ...,
            description="List of Material IDs to compute cohesive energies"
        ),
        normalization: str = Field(
            default="atom",
            description="The normalization to use, whether to normalize cohesive energy by number of atoms (deflaut) or by number of formula units  "
        )
):
    """
    Obtain the cohesive energy of the structure(s) corresponding to single or multiple material IDs

    Args:
        material_ids: List to material IDs to compute their cohesive energy
        normalization: Whether to normalize cohesive energy using number of atoms or number of formula

    Returns:
        str: The Markdown of  cohesive energies (in eV/atom or eV/formula unit) for
            each material, indexed by Material IDs .

    """
    logger.info("Getting cohesive energy for material IDs")

    with _get_mp_rester() as mpr:
        cohesive_energies = mpr.get_cohesive_energy(material_ids=material_ids, normalization=normalization)

    if not cohesive_energies:
        logger.info(f"No cohesive energy was retrived for {material_ids}")
        return f"No cohesive energies found for these Material IDs: {material_ids}"


    energies = f"## Cohesive Energies \n"
    for identifier, energy in cohesive_energies.items():
        unit = "eV/atom" if normalization == "atom" else "eV/formula unit"
        energies += f"-- **{identifier}** : {energy} {unit}\n"

    return energies


@mcp.tool()
async def get_atom_reference_data(
        funcs: tuple[str, ...] = Field(
            default=("PBE",),
            description="list of functionals to retrieve data for "
        )
) -> str:
    """
    Retrieve reference energies of isolated neutral atoms. this energies can be used to calculate formations energies of compounds,
    Write the meaning of these funcs eg thier full names
    Args:
        funcs ([str] or None ) : list of functionals to retrieve data for.
    Returns:
        str : Markdown containing isolated atom energies 
    """
    logger.info("Getting Atom Reference Data")
    with _get_mp_rester() as mpr:
        atom_data = mpr.get_atom_reference_data(funcs=funcs)

    if not atom_data:
        return f"No atom data retrieved for functionals {funcs}"

    atom_references = "| Element | Reference Energy (eV/atom) |\n"

    for element, energy in atom_data.items():
        atom_references += f"| **{element}** | {energy} | \n"

    return atom_references



@mcp.tool()
async def get_magnetic_data_by_id(
        material_ids: list[str] = Field(
            ...,
            description="Material ID of the material"
        ),
) -> str: 
    """
    Get magnetic data using material ID. The materials api provides computed
    magnetic propertics from Density Functional Theory (DFT) calculations. This includes
    1. Magnetic ordering
    2. Total Magnetization
    3. Site-projected Magnetic Moments
    4. Spin-polarized electronic structures

    Args:
        material_id: Material ID of the material e.g., mp-20664, which is Mn2Sb
    
    Returns:
        (str): returns a markdown string containing the magnetic data for the material.
    """
    logger.info(f"Getting magnetic data for material{material_ids}")
    with _get_mp_rester() as mpr:
        magnetic_data = mpr.magnetism.search(material_ids=material_ids)


    if not magnetic_data:
        logger.info(f"Not data collected for {material_ids}")
        return f"No magnetic data found for material {material_ids}"
    
    data_md  = f"|##      Magnetic Data for Material IDs    |\n\n"
    for idx, model in enumerate(magnetic_data): 
        data_md += f"idx : {idx}"
        data = model.model_dump()
        for key, value in data.items(): 
            data_md += f"| **{key}       :         {value}   |\n\n"
            
    
    return data_md
    


@mcp.tool()
async def get_charge_density_by_id(
        material_ids: str = Field(
            ...,
            description="Material ID of the material"
        )
):
    """
    Get charge density data for a given materials project ID
    Args:
        material_id: Material Project ID

    Returns:
        str :

    """
    logging.info(f"Getting charge density of material {material_ids}")
    with _get_mp_rester() as mpr:
        charge_density = mpr.get_charge_density_from_material_id(material_id=material_ids)
        logger.info(f"Charge density data retrieved for {material_ids}")
        

    if not charge_density:
        return f"No data found for material {charge_density}"

    density_data = f"""
            ## Material ID: {material_ids}
        
            ### Structure Summary:
            {charge_density.structure}
        
            ### Charge Density (Total):
            {charge_density.data["total"]}
            
            ### Is charge Density Polarized : 
            {charge_density.is_spin_polarized}
            
        """

    return density_data


@mcp.tool()
async def get_dielectric_data_by_id(
        material_id: str = Field(
            ..., 
            description="Material ID of the material"
    )
) -> str: 
    """
    Gets the dielectric data for a given material. Dielectric is a 
    material the can be polarized by an applied electric field. 
    The mathematical description of the dielectric effect is a tensor 
    constant of proportionality that relates an externally applied electric 
    field to the field within the material
    
    Args: 
        material_id (str): Material ID for the material

    Returns: 
        str: markdown of the dielectric data


    """
    logger.info(f"Getting Dielectric data for material: {material_id}")
    with _get_mp_rester() as mpr: 
        dielectric_data = mpr.materials.dielectric.search(material_id)
    
    if not dielectric_data: 
        logger.info(f"No data found for material {material_id}")
        return f"No data for the material: {material_id}"

    data_md  = f"|##    Dielectric Data  for Material IDs    |\n\n"
    for idx, model in enumerate(dielectric_data): 
        data_md += f"idx : {idx}"
        data = model.model_dump()
        for key, value in data.items(): 
            data_md += f"| **{key}    :    {value}   |\n\n"
            
    return data_md
    


@mcp.tool()
async def get_diffraction_patterns(
    material_id : str = Field(
        ..., 
        description="Material ID of the material "
    )
) -> str: 
    """
    Gets diffraction patterns of a material given its ID. 
    Diffraction occurs when waves (electrons, x-rays, neutrons)
    scattering from obstructions act as a secondary sources of propagations

    Args: 
        material id (str): the material id of the material to get the diffracton pattern 


    Return: 
        str: markdown of the patterns 
    
    """
    logger.info(f"Getting the Diffraction Pattern of element: {material_id}")
  
    with _get_mp_rester() as mpr: 
        # first retrieve the relevant structure 
        structure = mpr.get_structure_by_material_id(material_id)
    try: 
        sga = SpacegroupAnalyzer(structure=structure)
        conventional_structure  = sga.get_conventional_standard_structure()
        calculator = XRDCalculator(wavelength="CuKa")
        pattern = calculator.get_pattern(conventional_structure)
        return str(pattern)
    except: 
        logging.error("Error occurred when function get_diffraction_patterns ")
        return f"No diffraction pattern retrieved for material : {material_id}"
    


    
@mcp.tool()
async def get_xRay_absorption_spectra(
    material_ids: List[str] = Field(
        ..., 
        description="Material ID of the material"
    )
) -> str:
    """
    Obtain X-ray Absorption Spectra using single or multiple IDs,
    following the methodology as discussed by Mathew et al and Chen et al.

    Args: 
        material_ids (List[str]) : material_ids of the elements
    
    Return: 
        str: 
    
    """
    logging.info("")
    with _get_mp_rester() as mpr: 
        xas_doc = mpr.materials.xas.search(material_ids=material_ids)

    if not xas_doc: 
        logging.info(f"No data retrieve for material(s) : {material_ids}")
        return f"No data retrieve for material(s) : {material_ids}"

    data_md = f"|##  X-ray absorption spectra for Material IDs    |\n\n"
    for idx, model in enumerate(xas_doc):
        data_md += f"idx : {idx}"
        data = model.model_dump()
        for key, value in data.items():
            data_md += f"| **{key}  :  {value}   |\n\n"

    return data_md


@mcp.tool()
async def get_elastic_constants(
    material_ids: List[str] = Field(
        ..., 
        description="Material ID of the material"
    )
):
    """
    Obtain Elastic constants given material IDs.
    Elasticity describes a material's ability to resist deformations
    (i.e. size and shape) when subjected to external forces.

    :param material_ids :   material ID(s) of the elements

    :return:
        str: markdown of the elastic constants

    """
    logging.info(f"Getting Elastic Constant for material(s): {material_ids}")
    with _get_mp_rester() as mpr:
        elasticity_doc = mpr.materials.elasticity.search(material_ids=material_ids)

    if not elasticity_doc:
        return f"No Elasticity data retrieved for material: {material_ids}"

    data_md = f"|##     Elastic Constants   |\n\n"
    for idx, model in enumerate(elasticity_doc):
        data_md += f"idx : {idx}"
        data = model.model_dump()
        for key, value in data.items():
            data_md += f"| **{key}  :  {value}   |\n\n"

    return data_md



 

@mcp.tool()
async def get_suggested_substrates(
    material_id: str = Field(
        ..., 
        description="Material ID of the material"
    )
) -> str: 
    """
    Obtains Suggested substrates for a film material. 
    It helps to find suitable substrate materials for thin films 
    
    Args: 
        material_id (str): material ID of the material 
    
    Returns: 
        str: markdown of the data 
    
    """
    logging.info(f"Getting suggested substrates for the material : {material_id}")
    with _get_mp_rester() as mpr: 
        substrates_doc = mpr.materials.substrates.search(film_id=material_id)

    if not substrates_doc: 
        return f"No substrates gotten for material: {material_id}"

    sub_md = f""
    for idx, data in enumerate(substrates_doc):
        sub_md += f"## Substrate {idx + 1}\n\n"
        # Create a detailed view for each substrate
        sub_md += f"- **Index**: {idx}\n"
        sub_md += f"- **Substrate Formula**: {getattr(data, 'sub_form', 'N/A')}\n"
        sub_md += f"- **Substrate ID**: {getattr(data, 'sub_id', 'N/A')}\n"
        sub_md += f"- **Film Orientation**: {getattr(data, 'film_orient', 'N/A')}\n"
        sub_md += f"- **Area**: {getattr(data, 'area', 'N/A')}\n"
        sub_md += f"- **Energy**: {getattr(data, 'energy', 'N/A')}\n"
        sub_md += f"- **Film ID**: {getattr(data, 'film_id', 'N/A')}\n"
        sub_md += f"- **Orientation**: {getattr(data, 'orient', 'N/A')}\n\n"
        sub_md += "---\n\n"
    
    return sub_md



@mcp.tool()
async def get_thermo_stability(
    material_ids: List[str] = Field(
        ..., 
        description="Materials IDs of the material"
    ), 
    thermo_types: List[str] = Field(
        default=["GGA_GGA+U_R2SCAN"], 
        description=""
    )
) -> str: 
    """
    Obtains thermodynamic stability data for a material

    Args: 
        material_ids (List[str]) : A list of the material ID(s) eg. ["mp-861883"]
        thermo_types (List[str]) : 

    Returns: 
        str: Markdown of the thermodynamic stability data 
    
    """
    logging.info(f"Getting thermodynamic stability for material(s):{material_ids}")
    with _get_mp_rester() as mpr: 
        thermo_docs = mpr.materials.thermo.search(
            material_ids=material_ids, 
            thermo_types=thermo_types
        )

    if not thermo_docs: 
        logging.info("No thermodynamic stability data retrieved for ")
        return f"No thermodynamic stability data retrieved for materials: {material_ids}"
    
    thermo_md = f"Thermodynamic Stability for: {material_ids}"
    for idx, data in enumerate(thermo_docs):
        thermo_md += f"\n--- Material {idx + 1} ---\n"
       
        energy_above_hull = getattr(data, "energy_above_hull", "Not available")
        thermo_md += f"| ** | energy_above_hull : {energy_above_hull} | \n"
      
        formation_energy = getattr(data, "formation_energy_per_atom", "Not available")
        thermo_md += f"| ** | formation_energy_per_atom : {formation_energy} | \n"
        
        thermo_type = getattr(data, "thermo_type", "Not available")
        thermo_md += f"| ** | thermo_type : {thermo_type} | \n"
        
        
        is_stable = getattr(data, "is_stable", False)
        thermo_md += f"| ** | is_stable : {is_stable} | \n"
    
        formula_pretty = getattr(data, 'formula_pretty', 'Not available')
        thermo_md += f"| ** | formula : {formula_pretty} | \n"

    return thermo_md 





if __name__ == "__main__":
    logger.info("Starting Materials Project MCP server...")
    mcp.run()