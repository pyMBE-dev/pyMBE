from typing import Dict, Literal
from ..base_type import PMBBaseModel
from ..pint_quantity import PintQuantity


class BondTemplate(PMBBaseModel):
    pmb_type: Literal["bond"] = "bond"
    name: str                           # e.g. "HARMONIC_default"
    bond_type: str                      # "HARMONIC", "FENE"
    parameters: Dict[str, PintQuantity] # k, r0, d_r_max...
    l0: PintQuantity                    # initial bond length
