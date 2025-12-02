# pyMBE/storage/quantity_field.py
from dataclasses import dataclass
from typing import Any
from pint import UnitRegistry, Quantity

# dimension -> representative unit used to check dimensionality
_DIMENSION_REPRESENTATIVE = {
    "length": "nm",
    "energy": "meV",
    "energy/length**2": "meV/nm**2",
    "dimensionless": "dimensionless",
    # extend as needed
}


@dataclass
class PintQuantity:
    """
    Internal, SI-based stored representation of a Pint quantity.
    Stores magnitude and unit string using base/SI units.
    """

    magnitude: float
    units: str        # string representation of base units (e.g. "meter", "joule")
    dimension: str    # logical dimension: "length", "energy", ...

    @classmethod
    def from_quantity(cls, q: Quantity, expected_dimension: str, ureg: UnitRegistry):
        """
        Validate `q` has the expected dimension using the provided ureg,
        convert to base units (SI-like) and store magnitude + units.
        """
        if not isinstance(q, Quantity):
            raise TypeError("from_quantity expects a pint.Quantity")

        # Build a representative unit for the dimension using the provided registry
        if expected_dimension not in _DIMENSION_REPRESENTATIVE:
            raise ValueError(f"Unknown expected_dimension '{expected_dimension}'")

        rep_unit = ureg(_DIMENSION_REPRESENTATIVE[expected_dimension])

        # Use pint's dimensionality check
        try:
            if not q.check(rep_unit):
                raise ValueError(f"Quantity {q} does not have expected dimension '{expected_dimension}'")
        except Exception as e:
            # If check fails because registries differ, try converting via string (best-effort)
            raise

        # Use the dimension representative unit
        rep_unit_str = _DIMENSION_REPRESENTATIVE[expected_dimension]
        rep_unit = ureg(rep_unit_str)

        # Validate dimensionality
        if not q.check(rep_unit):
            raise ValueError(f"Quantity {q} does not match expected dimension '{expected_dimension}'")

        # Convert to the representative SI unit
        q_base = q.to(rep_unit)

        # Store magnitude and unit name
        mag = float(q_base.magnitude)
        unit_str = rep_unit_str
        return cls(magnitude=mag, units=unit_str, dimension=expected_dimension)

    def to_quantity(self, ureg: UnitRegistry) -> Quantity:
        """
        Reconstruct a pint.Quantity using the provided UnitRegistry.
        The units string should be parseable by ureg.
        """
        return self.magnitude * ureg(self.units)

    def to_dict(self) -> dict:
        return {"magnitude": self.magnitude, "units": self.units, "dimension": self.dimension}

    @classmethod
    def from_dict(cls, d: dict):
        return cls(magnitude=d["magnitude"], units=d["units"], dimension=d["dimension"])

