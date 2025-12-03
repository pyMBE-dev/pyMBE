from pyMBE.storage.base_type import PMBBaseModel
from pydantic import Field

class ProteinTemplate(PMBBaseModel):
    """
    Template defining a protein in a pyMBE simulation.

    Attributes:
        pmb_type (str): Fixed type identifier. Always "protein".
        name (str): Unique name of the protein template.
        model (str): Name or type of the model used for this protein.
        residue_list (List[str]): Ordered list of residue names that compose the protein.
        sequence (List[str]): Ordered sequence of residues representing the protein's structure.
    """
    pmb_type: str = Field(default="protein", frozen=True)
    name: str
    model: str
    residue_list: list[str] 
    sequence: list[str] 