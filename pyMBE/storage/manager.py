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

import pandas as pd
import json
import re
import numpy as np
import logging
import warnings

from typing import Dict, Any
from pyMBE.storage.templates.particle import ParticleTemplate
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.reactions.reaction import Reaction
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.templates.hydrogel import HydrogelTemplate
from pyMBE.storage.instances.hydrogel import HydrogelInstance
from pyMBE.storage.templates.lj import LJInteractionTemplate

TemplateType = Any  # union of template classes (ParticleTemplate, ResidueTemplate, ...)
InstanceType = Any  # union of instance classes (ParticleInstance, ResidueInstance, ...)

class Manager:
    """
    The canonical database manager for pyMBE.

    This class stores all templates, instances, and reactions in structured,
    explicit dictionaries.

    All I/O operations (CSV/JSON save/load) operate through DFManager.

    Attributes
    ----------
    ureg : UnitRegistry
        Pint unit registry used to reconstruct physical quantities from storage.

    templates : dict[str, dict[str, TemplateType]]
        Templates indexed by type and name.
        Example: templates["particle"]["A"] → ParticleTemplate

    instances : dict[str, dict[int, InstanceType]]
        Instances indexed by type and id.
        Example: instances["particle"][5] → ParticleInstance

    reactions : dict[str, Reaction]
        Chemical reactions keyed by reaction name.
    """

    def __init__(self,units):
        """
        Initialize an empty structured database.

        Args:
            ureg (UnitRegistry): Pint unit registry used to rebuild quantities.
        """
        self._units = units
        self._templates: Dict[str, Dict[str, TemplateType]] = {}
        self._instances: Dict[str, Dict[int, InstanceType]] = {}
        self._reactions: Dict[str, Reaction] = {}
        
        self._molecule_like_types = ["molecule",
                                     "peptide",
                                     "protein"]
        self._assembly_like_types = ["hydrogel"]
        self._pmb_types =  ["particle", "residue"] + self._molecule_like_types + self._assembly_like_types

    def _delete_bonds_of_particle(self, pid):
        """
        Delete all bond instances involving a given particle instance.

        Args:
            pid (int): The particle ID whose associated bonds should be deleted.

        Notes:
            - If no `"bond"` instances are present in the database, the method
            exits immediately.
            - This method does not raise errors if no bonds involve the particle.
            - It is intended for internal use by cascade-deletion routines.
        """
        if "bond" not in self._instances:
            return
        bonds_to_delete = [
            b_id for b_id, b in list(self._instances["bond"].items())
            if b.particle_id1 == pid or b.particle_id2 == pid
        ]
        for b_id in bonds_to_delete:
            del self._instances["bond"][b_id]
        if "bond" in self._instances and not self._instances["bond"]:
            del self._instances["bond"]

            if "bond" not in self._instances:
                return
            bonds_to_delete = [
                b_id for b_id, b in list(self._instances["bond"].items())
                if b.particle_id1 == pid or b.particle_id2 == pid
            ]
            for b_id in bonds_to_delete:
                del self._instances["bond"][b_id]
            if "bond" in self._instances and not self._instances["bond"]:
                del self._instances["bond"]

    def _find_instance_ids_by_attribute(self, pmb_type, attribute, value):
        """
        Return a list of instance IDs for a given pmb_type where a given attribute
        matches the requested value.

        Args:
            pmb_type (str): The pyMBE type to search within.
            attribute (str): The attribute name to match on (e.g. "residue_id", "molecule_id").
            value: The attribute value to match.

        Returns:
            List[int]: IDs of matching instances.
        """
        if pmb_type not in self._instances:
            return []
        results = []
        for inst_id, inst in self._instances[pmb_type].items():
            if hasattr(inst, attribute) and getattr(inst, attribute) == value:
                results.append(inst_id)
        return results

    def _find_instance_ids_by_name(self, pmb_type, name):
        """
        Return the IDs of all instances of a given pyMBE type that use a
        specific template name.

        This method inspects the instance registry stored under
        ``self._instances[pmb_type]`` and collects all instance identifiers
        whose ``instance.name`` matches the provided template name.

        Args:
            pmb_type (str):
                The instance category to search within.

            name (str):
                The template name associated with the instances of interest.

        Returns:
            list[int]:
                A list of instance IDs whose underlying template name matches
                ``name``. The list is empty if no such instances exist.

        Raises:
            KeyError:
                If ``pmb_type`` is not a recognized instance category.

        Examples:
            >>> db._find_instance_ids_by_name("particle", "A")
            [0, 3, 7][]

        Notes:
            - Only exact name matches are considered.
            - This method does not validate whether the corresponding template
            actually exists; it only inspects registered *instances*.
        """
        if pmb_type not in self._instances:
            return []

        result = []
        for iid, inst in self._instances[pmb_type].items():
            if hasattr(inst, "name") and inst.name == name:
                result.append(iid)

        return result

    def _find_template_types(self, name):
        """
        Return all pyMBE template categories that contain a template
        with a given name.

        Args:
            name (str):
                The template name to search for.

        Returns:
            list[str]:
                A list of PMB types (e.g., ``["particle", "residue"]``) in
                which a template named ``name`` exists. The list is empty if
                no such template is found.

        Examples:
            >>> db._find_template_types("A")
            ["particle"]

            >>> db._find_template_types("nonexistent")
            []
        """
        found = []
        for pmb_type, group in self._templates.items():
            if name in group:
                found.append(pmb_type)
        return found


    def _get_instances_df(self, pmb_type):
        """
        Returns a DataFrame containing all instance objects of a given pyMBE type.

        Args:
            pmb_type (str):
                The instance type to query. Must be a key in
                `self._instances`, such as `"particle"` or `"residue"`.

        Returns:
            pandas.DataFrame:
                A DataFrame where each row corresponds to one registered
                instance of the specified PMB type. If no instances exist,
                an empty DataFrame is returned.

        Notes:
            - Missing integer identifiers (e.g., `residue_id`) are stored as
            `pandas.NA` to ensure proper nullable integer handling.
            - Particle and residue instances receive custom row structures;
            all other instance types use direct model dumps.
        """
        rows = []
        if pmb_type not in self._instances:
            return pd.DataFrame(rows)
        for inst in self._instances[pmb_type].values():
            if pmb_type == "particle":
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "particle_id": inst.particle_id,
                    "initial_state": inst.initial_state,
                    "residue_id": int(inst.residue_id) if inst.residue_id is not None else pd.NA,
                    "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else pd.NA,
                    "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else pd.NA 
                })
            elif pmb_type == "residue":
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "residue_id": inst.residue_id,
                    "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else pd.NA,
                    "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else pd.NA
                })
            elif pmb_type in ["molecule","peptide","protein"]:
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "molecule_id": inst.molecule_id,
                    "assembly_id": int(inst.assembly_id) if inst.assembly_id is not None else pd.NA
                })

            else:
                # Generic representation for other types
                rows.append(inst.model_dump())
        return pd.DataFrame(rows)

    def _get_reactions_df(self):
        """
        Returns a DataFrame summarizing all registered chemical reactions.

        Returns:
            pandas.DataFrame:
                A DataFrame where each row corresponds to one reaction.

        Notes:
            - Participant objects are expected to expose ``state_name`` and
            ``coefficient`` attributes.
            - Stoichiometry is stored as a single dictionary per row to allow
            flexible downstream manipulation.
        """  
        rows = []
        for r in self._reactions.values():
            stoich = {
                f"{p.state_name}": p.coefficient
                for p in r.participants
            }
            rows.append({
                "reaction": r.name,
                "stoichiometry": stoich,
                "pK": r.pK,
                "reaction_type": r.reaction_type,
                "metadata": r.metadata,
            })
        return pd.DataFrame(rows)

    def _get_templates_df(self, pmb_type):
        """
        Returns a DataFrame containing all template definitions of a PMB type.

        Args:
            pmb_type (str):
                The template type to query, e.g. `"particle"`, `"residue"`,
                `"molecule"`.

        Returns:
            pandas.DataFrame:
                A DataFrame representing all templates of the given type.
                Particle templates expand to multiple rows, one per state.
                Empty DataFrame if no templates for that type exist.

        Notes:
            - Unit-bearing fields are converted to plain quantities through
            ``to_quantity(self._units)`` to maintain consistent I/O.
        """
        rows = []
        if pmb_type not in self._templates:
            return pd.DataFrame(rows)
        for tpl in self._templates[pmb_type].values():
            if pmb_type == "particle":
                for sname, st in tpl.states.items():
                    rows.append({
                    "pmb_type": tpl.pmb_type,
                    "name": tpl.name,
                    "sigma": tpl.sigma.to_quantity(self._units),
                    "epsilon": tpl.epsilon.to_quantity(self._units),
                    "cutoff": tpl.cutoff.to_quantity(self._units),
                    "offset": tpl.offset.to_quantity(self._units),
                    "initial_state": tpl.initial_state,
                    "state": sname,
                    "z": st.z,
                    "es_type": st.es_type
                })  
            elif pmb_type == "lj":
                shift = tpl.shift
                if isinstance(shift, dict) and {"magnitude", "units", "dimension"}.issubset(shift.keys()):
                    shift = tpl.shift.to_quantity(self._units)
                rows.append({
                    "pmb_type": tpl.pmb_type,
                    "name": tpl.name,
                    "state1": tpl.state1,
                    "state2": tpl.state2,
                    "sigma": tpl.sigma.to_quantity(self._units),
                    "epsilon": tpl.epsilon.to_quantity(self._units),
                    "cutoff": tpl.cutoff.to_quantity(self._units),
                    "offset": tpl.offset.to_quantity(self._units),
                    "shift": shift
                })

            elif pmb_type == "bond":
                parameters = {}
                for key in tpl.parameters.keys():
                    parameters[key] = tpl.parameters[key].to_quantity(self._units)
                rows.append({
                    "pmb_type": tpl.pmb_type,
                    "name": tpl.name,
                    "bond_type": tpl.bond_type,
                    "particle_name1": tpl.particle_name1,
                    "particle_name2": tpl.particle_name2,
                    "parameters": parameters,
                })

            else:
                # Generic representation for other types
                rows.append(tpl.model_dump())
        return pd.DataFrame(rows)

    def _has_instance(self, pmb_type, instance_id):
        """
        Check whether an instance with a given ID exists under a specific pyMBE type.

        Args:
            pmb_type (str):
                The instance category to search in. 

            instance_id (int):
                The unique identifier of the instance.

        Returns:
            bool:
                ``True`` if the instance exists in the given category,
                ``False`` otherwise.

        Raises:
            KeyError:
                If ``pmb_type`` is not a known instance category in the database.

        Examples:
            >>> db._has_instance("particle", 3)
            True

            >>> db._has_instance("nonexistent_type", 5)
            KeyError
        """
        if pmb_type not in self._instances:
            raise KeyError(f"Instance type '{pmb_type}' not found in the database.")

        return instance_id in self._instances[pmb_type]

    def _has_template(self, pmb_type, name):
        """
        Check whether a template with a given name exists within a specific pyMBE type.

        Args:
            pmb_type (str):
                The template category to search in (e.g. ``"particle"``,
                ``"bond"``, ``"molecule"``, ``"lj"``, etc.).
            name (str):
                The template name to check for.

        Returns:
            bool:
                ``True`` if a template named ``name`` exists under ``pmb_type``;
                ``False`` otherwise.

        Raises:
            KeyError:
                If ``pmb_type`` is not a recognized template category in the database.

        Examples:
            >>> db.has_template("particle", "A")
            True
        """
        if pmb_type not in self._templates:
            raise KeyError(f"Template type '{pmb_type}' not found in the database.")
        template_in_db = name in self._templates.get(pmb_type, {})
        return template_in_db

    def _register_instance(self, instance):
        """
        Register an instance of a pyMBE object.

        Args:
            instance: Any instance conforming to the pyMBE instance models.

        Raises:
            ValueError: If the id or instance model does not exist or is duplicated.
        """
        # infer pmb_type from instance class
        if isinstance(instance, ParticleInstance):
            pmb_type = "particle"
            iid = instance.particle_id
        elif isinstance(instance, ResidueInstance):
            pmb_type = "residue"
            iid = instance.residue_id
        elif isinstance(instance, MoleculeInstance):
            pmb_type = "molecule"
            iid = instance.molecule_id
        elif isinstance(instance, PeptideInstance):
            pmb_type = "peptide"
            iid = instance.molecule_id
        elif isinstance(instance, ProteinInstance):
            pmb_type = "protein"
            iid = instance.molecule_id
        elif isinstance(instance, BondInstance):
            pmb_type = "bond"
            iid = instance.bond_id
        elif isinstance(instance, HydrogelInstance):
            pmb_type = "hydrogel"
            iid = instance.assembly_id
        else:
            raise TypeError("Unsupported instance type")

        self._instances.setdefault(pmb_type, {})

        if iid in self._instances[pmb_type]:
            raise ValueError(f"Instance id {iid} already exists in type '{pmb_type}'")

        # validate template exists
        if instance.name not in self._templates.get(pmb_type, {}):
            raise ValueError(f"Template '{instance.name}' not found for type '{pmb_type}'")

        # validate state for particle instances
        if pmb_type == "particle":
            tpl: ParticleTemplate = self._templates[pmb_type][instance.name]
            if instance.initial_state not in tpl.states:
                raise ValueError(f"State '{instance.initial_state}' not defined in template '{instance.name}'")

        self._instances[pmb_type][iid] = instance

    def _register_reaction(self, reaction):
        """
        Register a chemical or physical reaction.

        Args:
            reaction (Reaction): Reaction object.

        Raises:
            ValueError: If reaction name already exists.
        """
        if reaction.name in self._reactions:
            raise ValueError(f"Reaction '{reaction.name}' already exists.")

        self._reactions[reaction.name] = reaction

    def _register_template(self, template):
        """
        Register a template.

        Args:
            template: Any template object conforming to the pyMBE template models.

        Raises:
            ValueError: If a template with the same name already exists.
        """
        pmb_type = getattr(template, "pmb_type", None)
        if pmb_type is None:
            # infer from class
            if isinstance(template, ParticleTemplate):
                pmb_type = "particle"
            elif isinstance(template, ResidueTemplate):
                pmb_type = "residue"
            elif isinstance(template, MoleculeTemplate):
                pmb_type = "molecule"
            elif isinstance(template, PeptideTemplate):
                pmb_type = "peptide"
            elif isinstance(template, ProteinTemplate):
                pmb_type = "protein"
            elif isinstance(template, HydrogelTemplate):
                pmb_type = "hydrogel"
            elif isinstance(template, BondTemplate):
                pmb_type = "bond"
            elif isinstance(template, LJInteractionTemplate):
                pmb_type = "lj"
            else:
                raise TypeError("Unknown template type; set attribute pmb_type or use supported templates")

        self._templates.setdefault(pmb_type, {})

        if template.name in self._templates[pmb_type]:
            raise ValueError(f"Template '{template.name}' exists in '{pmb_type}'")

        # particle templates must define at least one state
        if pmb_type == "particle":
            if not hasattr(template, "states") or len(template.states) == 0:
                raise ValueError("ParticleTemplate must define at least one state.")
            # ensure default_state valid if set
            if getattr(template, "default_state", None) is not None and template.default_state not in template.states:
                raise ValueError("default_state not in template states")

        self._templates[pmb_type][template.name] = template

    def _update_instance(self, instance_id, pmb_type, attribute, value):
        """
        Updates a single attribute of a registered instance.
        Only a restricted set of attributes is allowed for each PMB type,
        ensuring database consistency. 

        Args:
            instance_id (Hashable):
                Unique identifier of the instance to update. 
            pmb_type (str):
                Instance category, such as ``"particle"`` or ``"residue"``.
            attribute (str):
                Name of the field to update. 
            value (Any):
                New value to assign to the specified attribute.

        Raises:
            KeyError:
                If the provided ``instance_id`` does not exist for the given
                ``pmb_type``.
            ValueError:
                If attempting to modify an attribute that is not permitted
                for the instance's PMB type.

        Notes:
            - Allowed updates:
                * ``particle``: ``initial_state``, ``residue_id``, ``molecule_id``, ``assembly_id`` 
                * ``residue``: ``molecule_id``, ``assembly_id``
                * ``molecule``: ``assembly_id``
                * All other types: no attribute updates allowed.
            - The method replaces the instance with a new Pydantic model
            using ``model_copy(update=...)`` to maintain immutability and
            avoid partial mutations of internal state.
        """
        
        if instance_id not in self._instances[pmb_type]:
            raise KeyError(f"Instance '{instance_id}' not found for type '{pmb_type}' in the pyMBE database.")
                                
        if pmb_type == "particle":
            allowed = ["initial_state", "residue_id", "molecule_id", "assembly_id"]
        elif pmb_type == "residue":
            allowed = ["molecule_id", "assembly_id"]
        elif pmb_type == "molecule":
            allowed = ["assembly_id"]
        else:
            allowed = [None]  # No attributes allowed for other types        
                        
        if attribute not in allowed:
            raise ValueError(f"Attribute '{attribute}' not allowed for {pmb_type}. Allowed attributes: {allowed}")
        
        self._instances[pmb_type][instance_id] = self._instances[pmb_type][instance_id].model_copy(update={attribute: value})

    def _propagate_id(self, root_type, root_id, attribute, value):
        """
        Recursively updates an attribute (e.g., molecule_id, assembly_id)
        on an instance and all of its hierarchical descendants.

        Supported relationships:
            assembly → molecules → residues → particles
            molecule → residues → particles
            residue  → particles
            particle → (nothing)

        Args:
            root_type (str):
                One of {"assembly", "molecule", "residue", "particle"}.
            root_id (int):
                Instance ID of the root object to update.
            attribute (str):
                The attribute to update (e.g., "molecule_id", "assembly_id").
            value:
                The new value to assign.

        Returns:
            list[int]:
                A flat list of all instance IDs updated (including root).

        Raises:
            KeyError:
                If the root instance does not exist.
            ValueError:
                If an unsupported type or attribute is given.
        """
        updated = []
        # Map each type to its own identity attribute
        self_id_attribute = {
            "hydrogel": "assembly_id",
            "molecule": "molecule_id",
            "peptide": "molecule_id",
            "protein": "molecule_id",
            "residue": "residue_id",
            "particle": "particle_id",
        }
        assembly_types = ["hydrogel"]
        molecule_types = ["molecule", "peptide", "protein"]
        # 1) Update ROOT (unless attribute corresponds to its own ID)
        if attribute != self_id_attribute.get(root_type):
            self._update_instance(instance_id=root_id,
                                  pmb_type=root_type,
                                  attribute=attribute,
                                  value=value,)
            updated.append((root_type, root_id))
        # 2) Descendants: assembly → molecules
        if root_type in assembly_types:
            for mtype in molecule_types:
                molecule_ids = self._find_instance_ids_by_attribute(pmb_type=mtype,
                                                                    attribute="assembly_id",
                                                                    value=root_id)
                for mid in molecule_ids:
                    updated += self._propagate_id(root_type=mtype,   
                                                  root_id=mid,
                                                  attribute=attribute,
                                                  value=value)
        # 3) Descendants: molecule → residues
        if root_type in molecule_types:
            residue_ids = self._find_instance_ids_by_attribute(pmb_type="residue",
                                                               attribute="molecule_id",
                                                               value=root_id)
            for rid in residue_ids:
                updated += self._propagate_id(root_type="residue",
                                              root_id=rid,
                                              attribute=attribute,
                                              value=value)
        # 4) Descendants: residue → particles
        if root_type == "residue":
            particle_ids = self._find_instance_ids_by_attribute(pmb_type="particle",
                                                                attribute="residue_id",
                                                                value=root_id,)
            for pid in particle_ids:
                self._update_instance(instance_id=pid,
                                      pmb_type="particle",
                                      attribute=attribute,
                                      value=value,)
                updated.append(("particle", pid))
        return updated


    def _update_reaction_participant(self, reaction_name, particle_name, state_name, coefficient):
        """
        Append a new participant to an existing reaction in the database.
    
        Args:
            reaction_name (str):
                Name of the reaction to be updated.
            particle_name (str):
                Name of the particle template participating in the reaction.
            state_name (str):
                Specific state of the particle (e.g., protonation or charge state).
            coefficient (int):
                Stoichiometric coefficient for the new participant:
                - ``coefficient < 0`` → reactant  
                - ``coefficient > 0`` → product  
                Zero is not allowed.

        """
        if reaction_name not in self._reactions:
            raise KeyError(f"Reaction '{reaction_name}' not found in the pyMBE database.")

        rxn = self._reactions[reaction_name].add_participant(particle_name=particle_name,
                                                             state_name=state_name,
                                                             coefficient=coefficient)
        self._register_reaction(rxn)
        self._reactions.pop(reaction_name)
    
    def _propose_instance_id(self, pmb_type):
        """
        Propose the next available id for a new TypeInstance.

        If no instances of the given pmb_type exist, the proposed
        identifier is ``0``. Otherwise, the next available integer after the
        current maximum is returned.

        Returns:
            int: A non-negative integer that is not already used in the pyMBE database.

        Notes:
            - The method does not fill gaps; it always returns ``max + 1``.
        """
        if pmb_type in self._molecule_like_types:
            used_ids = []
            for t in self._molecule_like_types:
                if t in self._instances:
                    used_ids.extend(self._instances[t].keys())
            if not used_ids:
                return 0
        else:
            if pmb_type not in self._instances or len(self._instances[pmb_type]) == 0:
                return 0
            used_ids = list(self._instances[pmb_type].keys())
        return max(used_ids) + 1

    def delete_template(self, pmb_type, name):
        """
        Delete a template from the pyMBE database.

        This method removes a template identified by its pyMBE type and name.
        Before deletion, it checks whether any instance in the database uses
        this template. If any instance depends on it, a ``ValueError`` is raised
        to prevent breaking database integrity.

        Args:
            pmb_type (str):
                The template category.
            name (str):
                The name of the template to delete.

        Raises:
            KeyError:
                If the template type or name does not exist.
            ValueError:
                If one or more instances reference the template.
        """
        # Check template exists
        if pmb_type not in self._templates:
            raise KeyError(f"Template type '{pmb_type}' not found.")
        if name not in self._templates[pmb_type]:
            raise KeyError(f"Template '{name}' not found in type '{pmb_type}'.")

        # Check if any instance depends on this template
        if pmb_type in self._instances:
            for inst in self._instances[pmb_type].values():
                if getattr(inst, "name", None) == name:
                    raise ValueError(
                        f"Cannot delete template '{name}' from '{pmb_type}': "
                        f"Instance with ID {getattr(inst, pmb_type + '_id')} depends on it."
                    )

        # Delete
        del self._templates[pmb_type][name]

        # Delete empty groups
        if not self._templates[pmb_type]:
            del self._templates[pmb_type]

    def delete_instance(self, pmb_type, instance_id, cascade=False):
        """
        Delete an instance from the pyMBE database.

        Supports cascade deletion through the hierarchy:
            assembly → molecules → residues → particles → bonds
            molecule → residues → particles → bonds
            residue  → particles → bonds
            particle → bonds
            bond     → nothing

        Args:
            pmb_type (str):
                Category of the instance (particle, residue, molecule, peptide,
                protein, hydrogel, bond).
            instance_id (int):
                Unique identifier of the instance.
            cascade (bool):
                If True, automatically delete dependent objects.

        Raises:
            KeyError:
                If the instance does not exist.
            ValueError:
                If cascade is False but dependencies exist.
        """

        # ---- Basic checks ----
        if pmb_type not in self._instances:
            raise KeyError(f"Instance type '{pmb_type}' not found.")
        if instance_id not in self._instances[pmb_type]:
            raise KeyError(f"Instance ID '{instance_id}' not found in '{pmb_type}'.")
        inst = self._instances[pmb_type][instance_id]
        # ===============  CASCADE DELETION  =========================
        if cascade:
            # --- Delete children of ASSEMBLY-like objects ---
            if pmb_type in self._assembly_like_types:
                for mtype in self._molecule_like_types:
                    mids = self._find_instance_ids_by_attribute(pmb_type=mtype,
                                                                attribute="assembly_id",
                                                                value=instance_id,)
                    for mid in mids:
                        self.delete_instance(pmb_type=mtype, 
                                             instance_id=mid, 
                                             cascade=True)
                # delete particles inside the assembly *even if they have no residue/molecule* (e.g. nodes)
                pids = self._find_instance_ids_by_attribute(pmb_type="particle", 
                                                            attribute="assembly_id",
                                                            value=instance_id)
                for pid in pids:
                    self.delete_instance(pmb_type="particle", 
                                         instance_id=pid, 
                                         cascade=True)
            # --- Delete children of MOLECULE-like objects ---
            if pmb_type in self._molecule_like_types:
                residues = self._find_instance_ids_by_attribute(pmb_type="residue",
                                                                attribute="molecule_id",
                                                                value=instance_id,)
                for rid in residues:
                    self.delete_instance(pmb_type="residue", 
                                         instance_id=rid, 
                                         cascade=True)

            # --- Delete children of RESIDUE ---
            if pmb_type == "residue":
                particles = self._find_instance_ids_by_attribute(pmb_type="particle",
                                                                 attribute="residue_id",
                                                                 value=instance_id,)
                for pid in particles:
                    self.delete_instance(pmb_type="particle", 
                                         instance_id=pid, 
                                         cascade=True)

            # --- Delete children of PARTICLE (only bonds) ---
            if pmb_type == "particle":
                self._delete_bonds_of_particle(instance_id)

        # ===============  NON-CASCADE (SAFE DELETE)  ================
        else:
            # ---- ASSEMBLY-like: forbid deletion if molecules belong to it ----
            if pmb_type in self._assembly_like_types:
                for mtype in self._molecule_like_types:
                    mids = self._find_instance_ids_by_attribute(pmb_type=mtype,
                                                                attribute="assembly_id",
                                                                value=instance_id,)
                    if mids:
                        raise ValueError(f"{pmb_type} {instance_id} contains {mtype} instances {mids}. Use cascade=True to delete.")
            # ---- MOLECULE-like: check residues ----
            if pmb_type in self._molecule_like_types:
                residues = self._find_instance_ids_by_attribute(pmb_type="residue",
                                                                attribute="molecule_id",
                                                                value=instance_id,)
                if residues:
                    raise ValueError(f"{pmb_type} {instance_id} has residues {residues}. Use cascade=True to delete.")
            # ---- RESIDUE: check particles ----
            if pmb_type == "residue":
                particles = self._find_instance_ids_by_attribute(pmb_type="particle",
                                                                 attribute="residue_id",
                                                                 value=instance_id)
                if particles:
                    raise ValueError(f"Residue {instance_id} contains particles {particles}. Use cascade=True.")
            # ---- PARTICLE: check bonds and belonging ----
            if pmb_type == "particle":
                if inst.residue_id is not None:
                    raise ValueError(f"Particle {instance_id} belongs to residue {inst.residue_id}. Use cascade=True.")
                if inst.molecule_id is not None:
                    raise ValueError(f"Particle {instance_id} belongs to molecule {inst.molecule_id}. Use cascade=True.")
                if inst.assembly_id is not None:
                    raise ValueError(f"Particle {instance_id} belongs to assembly {inst.assembly_id}. "f"Use cascade=True.")
                bonds = [b_id for b_id, b in self._instances.get("bond", {}).items() if b.particle_id1 == instance_id or b.particle_id2 == instance_id]
                if bonds:
                    raise ValueError(f"Particle {instance_id} participates in bonds {bonds}. Use cascade=True.")
        # ===============  FINAL DELETION STEP  ======================
        del self._instances[pmb_type][instance_id]
        if not self._instances[pmb_type]:
            del self._instances[pmb_type]

    def get_instance(self, pmb_type, instance_id):
        """
        Retrieve a stored instance by type and instance_id.

        Looks up an instance within the internal instance registry
        (`self._instances`) using its pyMBE type (e.g., "particle", "residue",
        "bond", ...) and its unique id. If the instance does not exist,
        a `KeyError` is raised.

        Args:
            pmb_type (str): The instance pyMBE category.
            name (str): The unique name of the template to retrieve.

        Returns:
            InstanceType: The stored InstanceTemplate instance corresponding to the
            provided type and name.

        Raises:
            KeyError: If no template with the given type and name exists in
            the internal registry.
        """
        if instance_id not in self._instances[pmb_type]:
            raise KeyError(f"InstanceTemplate with id = '{instance_id}' not found in type '{pmb_type}'.")
        else:
            return self._instances[pmb_type][instance_id]

    def get_template(self, pmb_type, name):
        """
        Retrieve a stored template by type and name.

        Looks up a template within the internal template registry
        (`self._templates`) using its pyMBE type (e.g., "particle", "residue",
        "bond", ...) and its unique name. If the template does not exist,
        a `KeyError` is raised.

        Args:
            pmb_type (str): The template pyMBE category.
            name (str): The unique id of the template to retrieve.

        Returns:
            TemplateType: The stored template instance corresponding to the
            provided type and name.

        Raises:
            ValueError: If no template with the given type and name exists in
            the internal registry.
        """
        if name not in self._templates[pmb_type]:
            raise ValueError(f"Template '{name}' not found in type '{pmb_type}'.")
        else:
            return self._templates[pmb_type][name]

    def get_es_types_map(self):
        """
        Iterates over all particle templates and extracts the ESPResSo type (`es_type`)
        defined for each state. 

        Returns:
            dict[str, int]:
                A dictionary mapping each particle state to its corresponding ESPResSo type.

        """
        if "particle" not in self._templates:
            return {}          
        result = {}
        for _, tpl in self._templates["particle"].items():
            for state_name, state in tpl.states.items():
                result[state_name] = state.es_type
        return result
    
    def get_particle_id_map(self, object_name):
        """
        Collect all particle IDs associated with an object of given name in the
        pyMBE database. Works for particles, residues, molecules, proteins,
        peptides, and assemblies.

        Args:
            object_name (str): Name of the object.

        Returns:
            dict: {"all": [particle_ids],
                   "residue_map": {residue_id: [particle_ids]},
                   "molecule_map": {molecule_id: [particle_ids]},
                   "assembly_map": {assembly_id: [particle_ids]},}
        """
        # --- Determine object type by searching in the DB ------------------------
        object_type = None
        object_ids = []
        for pmb_type in self._pmb_types:
            if pmb_type in self._instances:
                for inst_id, inst in self._instances[pmb_type].items():
                    if getattr(inst, "name", None) == object_name:
                        object_type = pmb_type
                        object_ids.append(inst_id)

        if object_type is None:
            raise KeyError(f"No object named '{object_name}' found in database.")
        # Maps to return
        id_list = []
        residue_map = {}
        molecule_map = {}
        assembly_map = {}
        # Shortcut access to all particle instances
        particles = self._instances.get("particle", {})
        # Helper: group particle IDs by attribute (molecule_id, residue_id, assembly_id)
        def add_to_map(target_map, key, pid):
            if key is None:
                return
            target_map.setdefault(key, []).append(pid)
        # Case 1: object is a molecule-like type (molecule, protein, peptide)
        if object_type in self._molecule_like_types:
            for mol_id in object_ids:
                molecule_map[mol_id] = []
                for pid, p in particles.items():
                    if p.molecule_id == mol_id:
                        id_list.append(pid)
                        molecule_map[mol_id].append(pid)
                        add_to_map(residue_map, p.residue_id, pid)
                        add_to_map(assembly_map, p.assembly_id, pid)
        # Case 2: object is a residue
        elif object_type == "residue":
            for res_id in object_ids:
                residue_map[res_id] = []
                for pid, p in particles.items():
                    if p.residue_id == res_id:
                        id_list.append(pid)
                        residue_map[res_id].append(pid)
                        add_to_map(molecule_map, p.molecule_id, pid)
                        add_to_map(assembly_map, p.assembly_id, pid)
        # Case 3: object is a particle
        elif object_type == "particle":
            id_list.extend(object_ids)
            for pid in object_ids:
                p = particles[pid]
                add_to_map(residue_map, p.residue_id, pid)
                add_to_map(molecule_map, p.molecule_id, pid)
                add_to_map(assembly_map, p.assembly_id, pid)
        # Case 4: object is an assembly
        elif object_type == "assembly":
            for assembly_id in object_ids:
                assembly_map[assembly_id] = []
                for pid, p in particles.items():
                    if p.assembly_id == assembly_id:
                        id_list.append(pid)
                        assembly_map[assembly_id].append(pid)
                        add_to_map(molecule_map, p.molecule_id, pid)
                        add_to_map(residue_map, p.residue_id, pid)
        # Deduplicate + sort IDs
        id_list = sorted(set(id_list))
        return {"all": id_list, "molecule_map": molecule_map,  "residue_map": residue_map, "assembly_map": assembly_map,}