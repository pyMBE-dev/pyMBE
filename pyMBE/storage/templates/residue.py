from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class ResidueTemplate(PMBBaseModel):
    """
    Template defining a residue in a pyMBE simulation.

    Attributes:
        pmb_type (str): Fixed type identifier. Always "residue".
        name (str): Unique name of the residue template.
        central_bead (str): Name of the central bead representing the residue.
        side_chains (List[str]): List of side-chain names attached to the central bead.
            Defaults to an empty list if no side chains are present.
    """
    pmb_type: str = Field(default="residue", frozen=True)
    name: str
    central_bead: str 
    side_chains: list[str] = []
                