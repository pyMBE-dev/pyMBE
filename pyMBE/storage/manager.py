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
        self.units = units
        self.templates: Dict[str, Dict[str, TemplateType]] = {}
        self.instances: Dict[str, Dict[int, InstanceType]] = {}
        self.reactions: Dict[str, Reaction] = {}

    def _get_instances_df(self, pmb_type):
        """
        Returns a DataFrame containing all instance objects of a given pyMBE type.

        Args:
            pmb_type (str):
                The instance type to query. Must be a key in
                `self.instances`, such as `"particle"` or `"residue"`.

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
        if pmb_type not in self.instances:
            return pd.DataFrame(rows)
        for inst in self.instances[pmb_type].values():
            if pmb_type == "particle":
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "particle_id": inst.particle_id,
                    "initial_state": inst.initial_state,
                    "residue_id": int(inst.residue_id) if inst.residue_id is not None else pd.NA,
                    "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else pd.NA,
                })
            elif pmb_type == "residue":
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "residue_id": inst.residue_id,
                    "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else pd.NA,
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
        for r in self.reactions.values():
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
            ``to_quantity(self.units)`` to maintain consistent I/O.
        """
        rows = []
        if pmb_type not in self.templates:
            return pd.DataFrame(rows)
        for tpl in self.templates[pmb_type].values():
            if pmb_type == "particle":
                for sname, st in tpl.states.items():
                    rows.append({
                    "particle": tpl.name,
                    "sigma": tpl.sigma.to_quantity(self.units),
                    "epsilon": tpl.epsilon.to_quantity(self.units),
                    "cutoff": tpl.cutoff.to_quantity(self.units),
                    "offset": tpl.offset.to_quantity(self.units),
                    "state": sname,
                    "z": st.z,
                    "es_type": st.es_type
                })    
            else:
                # Generic representation for other types
                rows.append(tpl.model_dump())
        return pd.DataFrame(rows)

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
            iid = instance.hydrogel_id
        else:
            raise TypeError("Unsupported instance type")

        self.instances.setdefault(pmb_type, {})

        if iid in self.instances[pmb_type]:
            raise ValueError(f"Instance id {iid} already exists in type '{pmb_type}'")

        # validate template exists
        if instance.name not in self.templates.get(pmb_type, {}):
            raise ValueError(f"Template '{instance.name}' not found for type '{pmb_type}'")

        # validate state for particle instances
        if pmb_type == "particle":
            tpl: ParticleTemplate = self.templates[pmb_type][instance.name]
            if instance.initial_state not in tpl.states:
                raise ValueError(f"State '{instance.initial_state}' not defined in template '{instance.name}'")

        self.instances[pmb_type][iid] = instance

    def _register_reaction(self, reaction):
        """
        Register a chemical or physical reaction.

        Args:
            reaction (Reaction): Reaction object.

        Raises:
            ValueError: If reaction name already exists.
        """
        if reaction.name in self.reactions:
            raise ValueError(f"Reaction '{reaction.name}' already exists.")

        self.reactions[reaction.name] = reaction

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
            elif isinstance(template, BondTemplate):
                pmb_type = "bond"
            elif isinstance(template, PeptideTemplate):
                pmb_type = "peptide"
            elif isinstance(template, ProteinTemplate):
                pmb_type = "protein"
            elif isinstance(template, HydrogelTemplate):
                pmb_type = "hydrogel"
            else:
                raise TypeError("Unknown template type; set attribute pmb_type or use supported templates")

        self.templates.setdefault(pmb_type, {})

        if template.name in self.templates[pmb_type]:
            raise ValueError(f"Template '{template.name}' exists in '{pmb_type}'")

        # particle templates must define at least one state
        if pmb_type == "particle":
            if not hasattr(template, "states") or len(template.states) == 0:
                raise ValueError("ParticleTemplate must define at least one state.")
            # ensure default_state valid if set
            if getattr(template, "default_state", None) is not None and template.default_state not in template.states:
                raise ValueError("default_state not in template states")

        self.templates[pmb_type][template.name] = template

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
                * ``particle``: ``initial_state``, ``residue_id``, ``molecule_id``
                * ``residue``: ``molecule_id``
                * All other types: no attribute updates allowed.
            - The method replaces the instance with a new Pydantic model
            using ``model_copy(update=...)`` to maintain immutability and
            avoid partial mutations of internal state.
        """
        
        if instance_id not in self.instances[pmb_type]:
            raise KeyError(f"Instance '{instance_id}' not found in type '{pmb_type}'.")
                                
        if pmb_type == "particle":
            allowed = ["initial_state", "residue_id", "molecule_id"]
        elif pmb_type == "residue":
            allowed = ["molecule_id"]
        else:
            allowed = [None]  # No attributes allowed for other types        
                        
        if attribute not in allowed:
            raise ValueError(f"Attribute '{attribute}' not allowed for {pmb_type}. Allowed attributes: {allowed}")
        
        self.instances[pmb_type][instance_id] = self.instances[pmb_type][instance_id].model_copy(update={attribute: value})


    class _NumpyEncoder(json.JSONEncoder):
        """
        Custom JSON encoder that converts NumPy arrays to Python lists
        and NumPy scalars to Python scalars.
        """
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.generic):
                return obj.item()
            return super().default(obj)

    @classmethod
    def _add_bond_in_df(cls, df, particle_id1, particle_id2, use_default_bond=False):
        """
        Adds a bond entry on the `pymbe.df` storing the particle_ids of the two bonded particles.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            particle_id1(`int`): particle_id of the type of the first particle type of the bonded particles
            particle_id2(`int`): particle_id of the type of the second particle type of the bonded particles
            use_default_bond(`bool`, optional): Controls if a bond of type `default` is used to bond particle whose bond types are not defined in `pmb.df`. Defaults to False.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with the new bond added.
            index(`int`): Row index where the bond information has been added in pmb.df.
        """
        particle_name1 = df.loc[(df['particle_id']==particle_id1) & (df['pmb_type']=="particle")].name.values[0]
        particle_name2 = df.loc[(df['particle_id']==particle_id2) & (df['pmb_type']=="particle")].name.values[0]

        bond_key = cls._find_bond_key(df = df, 
                                      particle_name1 = particle_name1,
                                      particle_name2 = particle_name2,
                                      use_default_bond = use_default_bond)
        if not bond_key:
            return None
        df = cls._copy_df_entry(df = df,
                                name = bond_key,
                                column_name = 'particle_id2',
                                number_of_copies = 1)
        indexs = np.where(df['name'] == bond_key)
        index_list = list(indexs[0])
        used_bond_df = df.loc[df['particle_id2'].notnull()]
        #without this drop the program crashes when dropping duplicates because the 'bond' column is a dict
        used_bond_df = used_bond_df.drop([('bond_object','')],axis =1 )
        used_bond_index = used_bond_df.index.to_list()
        if not index_list:
            return None
        for index in index_list:
            if index not in used_bond_index:
                cls._clean_df_row(df = df,
                                  index = int(index))
                df.at[index,'particle_id'] = particle_id1
                df.at[index,'particle_id2'] = particle_id2
                break
        return df, index

    @classmethod
    def _add_value_to_df(cls, df, index,key,new_value, non_standard_value=False, overwrite=False):
        """
        Adds a value to a cell in the `pmb.df` DataFrame.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            index(`int`): index of the row to add the value to.
            key(`str`): the column label to add the value to.
            non_standard_value(`bool`, optional): Switch to enable insertion of non-standard values, such as `dict` objects. Defaults to False.
            overwrite(`bool`, optional): Switch to enable overwriting of already existing values in pmb.df. Defaults to False.
        """

        token = "#protected:"

        def protect(obj):
            if non_standard_value:
                return token + json.dumps(obj, cls=cls._NumpyEncoder)
            return obj

        def deprotect(obj):
            if non_standard_value and isinstance(obj, str) and obj.startswith(token):
                return json.loads(obj.removeprefix(token))
            return obj

        # Make sure index is a scalar integer value
        index = int(index)
        assert isinstance(index, int), '`index` should be a scalar integer value.'
        idx = pd.IndexSlice
        if cls._check_if_df_cell_has_a_value(df=df, index=index, key=key):
            old_value = df.loc[index,idx[key]]
            if not pd.Series([protect(old_value)]).equals(pd.Series([protect(new_value)])):
                name= df.loc[index,('name','')]
                pmb_type= df.loc[index,('pmb_type','')]
                logging.debug(f"You are attempting to redefine the properties of {name} of pmb_type {pmb_type}")    
                if overwrite:
                    logging.info(f'Overwritting the value of the entry `{key}`: old_value = {old_value} new_value = {new_value}')
                if not overwrite:
                    logging.debug(f"pyMBE has preserved of the entry `{key}`: old_value = {old_value}. If you want to overwrite it with new_value = {new_value}, activate the switch overwrite = True ")
                    return

        df.loc[index,idx[key]] = protect(new_value)
        if non_standard_value:
            df[key] = df[key].apply(deprotect)
        return
    
    @classmethod
    def _assign_molecule_id(cls, df, molecule_index):
        """
        Assigns the `molecule_id` of the pmb object given by `pmb_type`
        
        Args:
            molecule_index(`int`): index of the current `pmb_object_type` to assign the `molecule_id`
        Returns:
            molecule_id(`int`): Id of the molecule
        """
        cls._clean_df_row(df = df,
                          index = int(molecule_index))

        if df['molecule_id'].isnull().values.all():
            molecule_id = 0        
        else:
            molecule_id = df['molecule_id'].max() +1
        cls._add_value_to_df(df = df,
                             key = ('molecule_id',''),
                             index = int(molecule_index),
                             new_value = molecule_id)
        return molecule_id

    @staticmethod
    def _check_if_df_cell_has_a_value(df, index, key):
        """
        Checks if a cell in the `pmb.df` at the specified index and column has a value.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            index(`int`): Index of the row to check.
            key(`str`): Column label to check.

        Returns:
            `bool`: `True` if the cell has a value, `False` otherwise.
        """
        idx = pd.IndexSlice
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return not pd.isna(df.loc[index, idx[key]])

    @staticmethod
    def _check_if_name_is_defined_in_df(name, df):
        """
        Checks if `name` is defined in `pmb.df`.

        Args:
            name(`str`): label to check if defined in `pmb.df`.
            df(`DataFrame`): dataframe with pyMBE information.

        Returns:
            `bool`: `True` for success, `False` otherwise.
        """
        return name in df['name'].unique()

    @staticmethod
    def _check_if_multiple_pmb_types_for_name(name, pmb_type_to_be_defined, df):
        """
        Checks if `name` is defined in `pmb.df` with multiple pmb_types.

        Args:
            name(`str`): label to check if defined in `pmb.df`.
            pmb_type_to_be_defined(`str`): pmb object type corresponding to `name`.
            df(`DataFrame`): dataframe with pyMBE information.

        Returns:
            `bool`: `True` for success, `False` otherwise.
        """
        if name in df['name'].unique():
            current_object_type = df[df['name']==name].pmb_type.values[0]
            if current_object_type != pmb_type_to_be_defined:
                raise ValueError (f"The name {name} is already defined in the df with a pmb_type = {current_object_type}, pyMBE does not support objects with the same name but different pmb_types")

    @classmethod
    def _clean_df_row(cls, df, index, columns_keys_to_clean=("particle_id", "particle_id2", "residue_id", "molecule_id")):
        """
        Cleans the columns of `pmb.df` in `columns_keys_to_clean` of the row with index `index` by assigning them a pd.NA value.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            index(`int`): Index of the row to clean.
            columns_keys_to_clean(`list` of `str`, optional): List with the column keys to be cleaned. Defaults to [`particle_id`, `particle_id2`, `residue_id`, `molecule_id`].
        """   
        for column_key in columns_keys_to_clean:
            cls._add_value_to_df(df = df,
                                 key = (column_key,''),
                                 index = index,
                                 new_value = pd.NA)
        df.fillna(pd.NA, 
                  inplace = True)

    @staticmethod
    def _clean_ids_in_df_row(df, row):
        """
        Cleans particle, residue and molecules ids in `row`.
        If there are other repeated entries for the same name, drops the row.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            row(pd.DataFrame): A row from the DataFrame to clean.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with cleaned ids in `row`
        """
        columns_to_clean = ['particle_id',
                            'particle_id2', 
                            'residue_id', 
                            'molecule_id']
        if len(df.loc[df['name'] == row['name'].values[0]]) > 1:
            df = df.drop(row.index).reset_index(drop=True)
            
        else:
            for column_name in columns_to_clean:
                df.loc[row.index, column_name] = pd.NA
        return df

    @staticmethod
    def _copy_df_entry(df, name, column_name, number_of_copies):
        '''
        Creates 'number_of_copies' of a given 'name' in `pymbe.df`.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            name(`str`): Label of the particle/residue/molecule type to be created. `name` must be defined in `pmb.df`
            column_name(`str`): Column name to use as a filter. 
            number_of_copies(`int`): number of copies of `name` to be created.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with the new copies of `name` added.

        Note:
            - Currently, column_name only supports "particle_id", "particle_id2", "residue_id" and "molecule_id" 
        '''
        valid_column_names=["particle_id", "residue_id", "molecule_id", "particle_id2" ]
        if column_name not in valid_column_names:
            raise ValueError(f"{column_name} is not a valid column_name, currently only the following are supported: {valid_column_names}")
        df_by_name = df.loc[df.name == name]
        if number_of_copies != 1:                           
            df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies-1), ignore_index=True)
            # Concatenate the new particle rows to  `df`
            df = pd.concat ([df,df_by_name_repeated], ignore_index=True)
        else:
            if not df_by_name[column_name].isnull().values.any():     
                df_by_name = df_by_name[df_by_name.index == df_by_name.index.min()] 
                df_by_name_repeated = pd.concat ([df_by_name]*(number_of_copies), ignore_index=True)
                df_by_name_repeated[column_name] = pd.NA
                df = pd.concat ([df,df_by_name_repeated], ignore_index=True)
        return df

    @staticmethod
    def _create_variable_with_units(variable, units_registry):
        """
        Returns a pint object with the value and units defined in `variable`.

        Args:
            variable(`dict` or `str`): {'value': value, 'units': units}
            units_registry(`pint.UnitRegistry`): pyMBE UnitRegistry object.

        Returns:
            variable_with_units(`obj`): variable with units using the pyMBE UnitRegistry.
        """        
        if isinstance(variable, dict):
            value=variable.pop('value')
            units=variable.pop('units')
        elif isinstance(variable, str):
            value = float(re.split(r'\s+', variable)[0])
            units = re.split(r'\s+', variable)[1]
        variable_with_units = value * units_registry(units)
        return variable_with_units

    @classmethod
    def _convert_columns_to_original_format(cls,df,units_registry):
        """
        Converts the columns of the Dataframe to the original format in pyMBE.
        
        Args:
            df(`DataFrame`): dataframe with pyMBE information as a string
            units_registry(`pint.UnitRegistry`): pyMBE UnitRegistry object.  
        
        """

        columns_dtype_int = ['particle_id','particle_id2', 'residue_id','molecule_id', ('state_one','es_type'),('state_two','es_type'),('state_one','z'),('state_two','z') ]  

        columns_with_units = ['sigma', 'epsilon', 'cutoff', 'offset']

        columns_with_list_or_dict = ['residue_list','side_chains', 'parameters_of_the_potential','sequence', 'chain_map', 'node_map']

        for column_name in columns_dtype_int:
            df[column_name] = df[column_name].astype(pd.Int64Dtype())
            
        for column_name in columns_with_list_or_dict:
            if df[column_name].isnull().all():
                df[column_name] = df[column_name].astype(object)
            else:
                df[column_name] = df[column_name].apply(lambda x: json.loads(x) if pd.notnull(x) else x)

        for column_name in columns_with_units:
            df[column_name] = df[column_name].apply(lambda x: cls._create_variable_with_units(x, units_registry) if pd.notnull(x) else x)

        df['bond_object'] = df['bond_object'].apply(lambda x: cls._convert_str_to_bond_object(x) if pd.notnull(x) else x)
        df["l0"] = df["l0"].astype(object)
        df["pka"] = df["pka"].astype(object)

    @staticmethod
    def _convert_str_to_bond_object(bond_str):
        """
        Convert a row read as a `str` to the corresponding ESPResSo bond object. 

        Args:
            bond_str(`str`): string with the information of a bond object.

        Returns:
            bond_object(`obj`): ESPResSo bond object.

        Note:
            Currently supported bonds are: HarmonicBond and FeneBond
        """
        import espressomd.interactions

        supported_bonds = ['HarmonicBond', 'FeneBond']
        m = re.search(r'^([A-Za-z0-9_]+)\((\{.+\})\)$', bond_str)
        if m is None:
            raise ValueError(f'Cannot parse bond "{bond_str}"')
        bond = m.group(1)
        if bond not in supported_bonds:
            raise NotImplementedError(f"Bond type '{bond}' currently not implemented in pyMBE, accepted types are {supported_bonds}")
        params = json.loads(m.group(2))
        bond_id = params.pop("bond_id")
        bond_object = getattr(espressomd.interactions, bond)(**params)
        bond_object._bond_id = bond_id
        return bond_object

    @staticmethod
    def _delete_entries_in_df(df, entry_name):
        """
        Deletes entries with name `entry_name` from the DataFrame if it exists.

        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            entry_name (`str`): The name of the entry in the dataframe to delete.

        Returns:
            df(`DataFrame`): dataframe with pyMBE information with the entry deleted.
        """
        if entry_name in df["name"].values:
            df = df[df["name"] != entry_name].reset_index(drop=True)
        return df

    @staticmethod
    def _find_bond_key(df, particle_name1, particle_name2, use_default_bond=False):
        """
        Searches for the `name` of the bond between `particle_name1` and `particle_name2` in `pymbe.df` and returns it.
        
        Args:
            df(`DataFrame`): dataframe with pyMBE information.
            particle_name1(`str`): label of the type of the first particle type of the bonded particles.
            particle_name2(`str`): label of the type of the second particle type of the bonded particles.
            use_default_bond(`bool`, optional): If it is activated, the "default" bond is returned if no bond is found between `particle_name1` and `particle_name2`. Defaults to 'False'. 

        Returns:
            bond_key (str): `name` of the bond between `particle_name1` and `particle_name2` if a matching bond exists

        Note:
            - If `use_default_bond`=`True`, it returns "default" if no key is found.
        """
        bond_keys = [f'{particle_name1}-{particle_name2}', f'{particle_name2}-{particle_name1}']
        bond_defined=False
        for bond_key in bond_keys:
            if bond_key in df["name"].values:
                bond_defined=True
                correct_key=bond_key
                break
        if bond_defined:
            return correct_key
        elif use_default_bond:
            return 'default'
        else:
            return None

    @staticmethod
    def _setup_df():
        """
        Sets up the pyMBE's dataframe `pymbe.df`.

        Returns:
            columns_names(`obj`): pandas multiindex object with the column names of the pyMBE's dataframe
        """
        
        columns_dtypes = {
            'name': {
                '': str},
            'pmb_type': {
                '': str},
            'particle_id': {
                '': pd.Int64Dtype()},
            'particle_id2':  {
                '': pd.Int64Dtype()},
            'residue_id':  {
                '': pd.Int64Dtype()},
            'molecule_id':  {
                '': pd.Int64Dtype()},
            'acidity':  {
                '': str},
            'pka':  {
                '': object},
            'central_bead':  {
                '': object},
            'side_chains': {
                '': object},
            'residue_list': {
                '': object},
            'model': {
                '': str},
            'sigma': {
                '': object},
            'cutoff': {
                '': object},
            'offset': {
                '': object},
            'epsilon': {
                '': object},
            'state_one': {
                'label': str,
                'es_type': pd.Int64Dtype(),
                'z': pd.Int64Dtype()},
            'state_two': {
                'label': str,
                'es_type': pd.Int64Dtype(),
                'z': pd.Int64Dtype()},
            'sequence': {
                '': object},
            'bond_object': {
                '': object},
            'parameters_of_the_potential':{
                '': object},
            'l0': {
                '': float},
            'node_map':{
                '':object},
            'chain_map':{
                '':object}}
        
        df = pd.DataFrame(columns=pd.MultiIndex.from_tuples([(col_main, col_sub) for col_main, sub_cols in columns_dtypes.items() for col_sub in sub_cols.keys()]))
        
        for level1, sub_dtypes in columns_dtypes.items():
            for level2, dtype in sub_dtypes.items():
                df[level1, level2] = df[level1, level2].astype(dtype)

        columns_names = pd.MultiIndex.from_frame(df)
        columns_names = columns_names.names
                
        return df