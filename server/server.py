import os
import logging
from typing import Optional, List
from mcp.server.fastmcp import FastMCP
from pydantic import Field

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
    Helper function to initialize a MPRester session with the user's API key.
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
        description="List of element symbols to filter by (e.g. ['Si', 'O']).",
    ),
    band_gap_min: float = Field(
        default=0.0, description="Lower bound for band gap filtering in eV."
    ),
    band_gap_max: float = Field(
        default=10.0, description="Upper bound for band gap filtering in eV."
    ),
    is_stable: bool = Field(
        default=False,
        description="Whether to only retrieve stable materials (True) or all (False).",
    ),
    max_results: int = Field(
        default=20, ge=1, le=200, description="Maximum number of results to return."
    ),
) -> str:
    """
    Search for materials in the Materials Project database using basic filters:
    - elements (list of elements to include)
    - band_gap (range in eV)
    - is_stable
    Returns a formatted list of matches with their material_id, formula, and band gap.
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
    material_id: str = Field(..., description="Materials Project ID (e.g. 'mp-149')")
) -> str:
    """
    Retrieve the final computed structure for a given material_id from the Materials Project.
    Returns a plain text summary of the structure (lattice, sites, formula).
    """
    logger.info(f"Fetching structure for {material_id}...")
    with _get_mp_rester() as mpr:
        # Shortcut method to get just the final structure
        structure = mpr.get_structure_by_material_id(material_id)

    if not structure:
        return f"No structure found for {material_id}."

    # Summarize the structure
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


if __name__ == "__main__":
    logger.info("Starting Materials Project MCP server...")
    mcp.run()