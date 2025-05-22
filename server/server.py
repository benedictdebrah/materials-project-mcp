import os
import logging
from typing import Optional, List, Union
from mcp.server.fastmcp import FastMCP
from pydantic import Field
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.plotter import PhononBSPlotter
from emmet.core.electronic_structure import BSPathType
from typing import Literal

# Materials Project client
from mp_api.client import MPRester

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("materials_project_mcp")

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


    return ion_reference_data






if __name__ == "__main__":
    logger.info("Starting Materials Project MCP server...")
    mcp.run()