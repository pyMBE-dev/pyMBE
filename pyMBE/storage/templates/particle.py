#
# Copyright (C) 2025 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from typing import Dict, Literal, Optional 
from pydantic import Field, field_validator

from ..base_type import PMBBaseModel
from ..pint_quantity import PintQuantity

class ParticleState(PMBBaseModel):
    """
    Represents a single state of a particle in pyMBE.

    Attributes:
        pmb_type (Literal["particle_state"]): Fixed type identifier. Always "particle_state".
        name (str): Name of the particle state, e.g., "HA", "A-", "H+".
        z (int): Charge of the particle in this state.
        es_type (float): Identifier for the state used in Espresso simulations.
    """
    pmb_type: Literal["particle_state"] = "particle_state"
    name: str                      # e.g. "HA", "A-", "H+"
    z: int
    es_type: int                  # label in espresso


class ParticleTemplate(PMBBaseModel):
    """
    Template describing a particle type, including interaction parameters and allowed states.

    Attributes:
        pmb_type (str): Fixed type identifier. Always "particle".
        sigma (PintQuantity): Particle diameter or size parameter.
        epsilon (PintQuantity): Depth of the LJ potential well (interaction strength).
        cutoff (PintQuantity): Cutoff distance for the LJ potential.
        offset (PintQuantity): Offset distance for the LJ potential.
        states (Dict[str, ParticleState]): Dictionary of allowed particle states.
            Keys are state names, values are ParticleState instances.
        initial_state (Optional[str]): Name of the default particle state.
            If not provided explicitly, the first added state becomes the initial state.
    """

    pmb_type: str = Field(default="particle", frozen=True)
    name : str
    sigma: PintQuantity
    cutoff: PintQuantity
    offset: PintQuantity
    epsilon: PintQuantity
    states: Dict[str, ParticleState] = {}
    initial_state: Optional[str] = None

    def add_state(self, state):
        """
        Add a new state to the particle template.

        This method registers a new `ParticleState` in the template's `states` dictionary.
        If a state with the same name already exists, a `ValueError` is raised.

        Args:
            state (ParticleState): The particle state to add.

        Raises:
            ValueError: If a state with the same name already exists in the template.
        """
        if state.name in self.states:
            raise ValueError(f"State {state.name} already exists in template {self.name}")
        self.states[state.name] = state

         # Automatically assign initial state if this is the first state
        if self.initial_state is None:
            self.initial_state = state.name

    def get_lj_parameters(self, ureg):
        """
        Retrieve the Lennard-Jones interaction parameters for the particle template.

        Args:
            ureg (pint.UnitRegistry) : Pint unit registry used to reconstruct physical quantities from storage.
        
        Returns:
            Dict[str, pint.Quantity]:
                A dictionary containing the following LJ parameters: sigma, epsilon, cutoff, offset.

        Example:
            >>> tpl = ParticleTemplate(...)
            >>> params = tpl.get_lj_parameters()
            >>> params["sigma"]
            <Quantity(1.0, 'nanometer')>
        """
        return {"sigma": self.sigma.to_quantity(ureg), 
                "epsilon": self.epsilon.to_quantity(ureg),
                "cutoff": self.cutoff.to_quantity(ureg),
                "offset": self.offset.to_quantity(ureg)}

