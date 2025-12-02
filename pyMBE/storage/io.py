# pyMBE/storage/io.py
import os
import json
from pathlib import Path
from typing import Any, Dict

import pandas as pd
from pint import UnitRegistry

from pyMBE.storage.df_management import _DFManagement
from pyMBE.storage.pint_quantity import PintQuantity
from pyMBE.storage.templates.particle import ParticleTemplate, ParticleState
from pyMBE.storage.templates.residue import ResidueTemplate
from pyMBE.storage.templates.molecule import MoleculeTemplate
from pyMBE.storage.templates.bond import BondTemplate
from pyMBE.storage.instances.particle import ParticleInstance
from pyMBE.storage.instances.residue import ResidueInstance
from pyMBE.storage.instances.molecule import MoleculeInstance
from pyMBE.storage.instances.bond import BondInstance
from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
from pyMBE.storage.templates.peptide import PeptideTemplate
from pyMBE.storage.instances.peptide import PeptideInstance
from pyMBE.storage.templates.protein import ProteinTemplate
from pyMBE.storage.instances.protein import ProteinInstance
from pyMBE.storage.templates.hydrogel import HydrogelTemplate, HydrogelNode, HydrogelChain
from pyMBE.storage.instances.hydrogel import HydrogelInstance

# ----------------------------------------------------------------------
# Helpers for JSON encode/decode
# ----------------------------------------------------------------------
def _encode(obj: Any) -> str:
    """
    Convert Python object -> JSON string.
    - If obj is a PintQuantity (PintQuantity dataclass), convert via to_dict().
    - If obj is a Pydantic model, prefer its model_dump() result.
    - Otherwise, json.dumps(obj).
    """
    if obj is None:
        return ""
    # PintQuantity dataclass (has to_dict)
    if isinstance(obj, PintQuantity):
        return json.dumps(obj.to_dict(), separators=(",", ":"), ensure_ascii=False)

    # If it's already a dict/list/scalar, json-dump it
    try:
        return json.dumps(obj, separators=(",", ":"), ensure_ascii=False)
    except TypeError:
        # Last resort: convert to string
        return json.dumps(str(obj), separators=(",", ":"), ensure_ascii=False)


def _decode(s: Any) -> Any:
    """
    Robust JSON decoder:
      - Returns None for None/NaN/empty/non-string values
      - If s is a JSON string, returns parsed Python object
      - If s is already a dict/list/number, returns it unchanged
    """
    # None / pandas NA / nan handling
    if s is None:
        return None
    # pandas often gives float('nan') or numpy.nan
    if isinstance(s, float):
        # NaN -> None
        if pd.isna(s):
            return None
        return s
    # If already native Python container
    if isinstance(s, (dict, list, int, bool)):
        return s
    # Must be a string to parse JSON
    if not isinstance(s, str):
        return None
    s_str = s.strip()
    if s_str == "" or s_str.lower() == "nan":
        return None
    try:
        return json.loads(s_str)
    except Exception:
        # If it fails, try to interpret as plain string
        return s_str

# ----------------------------------------------------------------------
# SAVE
# ----------------------------------------------------------------------
def save_database_csv(db: _DFManagement, folder: str):
    """
    Save all database content to CSV files in the specified folder.
    Files produced:
      - templates_particle.csv
      - templates_residue.csv
      - templates_molecule.csv
      - templates_bond.csv
      - instances_particle.csv
      - instances_residue.csv
      - instances_molecule.csv
      - instances_bond.csv
      - reactions.csv
    """
    os.makedirs(folder, exist_ok=True)

    # -----------------------------
    # TEMPLATES
    # -----------------------------
    for pmb_type, tpl_dict in db.templates.items():
        rows = []
        for tpl in tpl_dict.values():
            # PARTICLE TEMPLATE: explicit custom encoding
            if pmb_type == "particle" and isinstance(tpl, ParticleTemplate):
                rows.append({
                    "name": tpl.name,
                    "sigma": _encode(tpl.sigma.to_dict()),
                    "epsilon": _encode(tpl.epsilon.to_dict()),
                    "cutoff": _encode(tpl.cutoff.to_dict()),
                    "offset": _encode(tpl.offset.to_dict()),
                    # states: dict state_name -> ParticleState.model_dump()
                    "states": _encode({sname: st.model_dump() for sname, st in tpl.states.items()}),
                })

            # RESIDUE TEMPLATE
            elif pmb_type == "residue" and isinstance(tpl, ResidueTemplate):
                rows.append({
                    "name": tpl.name,
                    "central_bead": tpl.central_bead,
                    "side_chains": _encode(tpl.side_chains),
                })

            # MOLECULE TEMPLATE
            elif pmb_type == "molecule" and isinstance(tpl, MoleculeTemplate):
                rows.append({
                    "name": tpl.name,
                    "residue_list": _encode(tpl.residue_list),
                })
            
            elif pmb_type == "peptide" and isinstance(tpl, PeptideTemplate):
                rows.append({
                    "name": tpl.name,
                    "model": tpl.model,
                    "residue_list": _encode(tpl.residue_list),
                    "sequence": _encode(tpl.sequence),
                })
            elif pmb_type == "protein" and isinstance(tpl, ProteinTemplate):
                rows.append({
                    "name": tpl.name,
                    "model": tpl.model,
                    "residue_list": _encode(tpl.residue_list),
                    "sequence": _encode(tpl.sequence),
                })
            # BOND TEMPLATE
            elif pmb_type == "bond" and isinstance(tpl, BondTemplate):
                # parameters: dict[str, scalar or PintQuantity]
                params_serial = {}
                for k, v in tpl.parameters.items():
                    if isinstance(v, PintQuantity):
                        params_serial[k] = v.to_dict()
                    else:
                        # assume scalar serializable
                        params_serial[k] = v
                rows.append({
                    "name": tpl.name,
                    "bond_type": tpl.bond_type,
                    "parameters": _encode(params_serial),
                    "l0": _encode(tpl.l0.to_dict()),
                })
            # HYDROGEL TEMPLATE
            elif pmb_type == "hydrogel" and isinstance(tpl, HydrogelTemplate):
                rows.append({
                    "name": tpl.name,
                    "node_map": _encode([node.model_dump() for node in tpl.node_map]),
                    "chain_map": _encode([chain.model_dump() for chain in tpl.chain_map]),
                })
            else:
                # Generic fallback: try model_dump()
                try:
                    rows.append(tpl.model_dump())
                except Exception:
                    rows.append({"name": getattr(tpl, "name", None)})

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(folder, f"templates_{pmb_type}.csv"), index=False)

    # -----------------------------
    # INSTANCES
    # -----------------------------
    for pmb_type, inst_dict in db.instances.items():
        rows = []
        for inst in inst_dict.values():
            if pmb_type == "particle" and isinstance(inst, ParticleInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "particle_id": int(inst.particle_id),
                    "initial_state": inst.initial_state,
                    "residue_id": int(inst.residue_id) if inst.residue_id is not None else "",
                    "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else "",
                })
            elif pmb_type == "residue" and isinstance(inst, ResidueInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "residue_id": int(inst.residue_id),
                    "molecule_id": int(inst.molecule_id) if inst.molecule_id is not None else "",
                })
            elif pmb_type == "molecule" and isinstance(inst, MoleculeInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "molecule_id": int(inst.molecule_id),
                })
            elif pmb_type == "peptide" and isinstance(inst, PeptideInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "molecule_id": int(inst.molecule_id),
                })
            elif pmb_type == "protein" and isinstance(inst, ProteinInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "molecule_id": int(inst.molecule_id),
                })
            elif pmb_type == "bond" and isinstance(inst, BondInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "bond_id": int(inst.bond_id),
                    "particle_id1": int(inst.particle_id1),
                    "particle_id2": int(inst.particle_id2),
                })
            elif pmb_type == "hydrogel" and isinstance(inst, HydrogelInstance):
                rows.append({
                    "pmb_type": pmb_type,
                    "name": inst.name,
                    "hydrogel_id": int(inst.hydrogel_id),
                    "molecule_ids": _encode(inst.molecule_ids),
                })
            else:
                # fallback to model_dump
                try:
                    rows.append(inst.model_dump())
                except Exception:
                    rows.append({"name": getattr(inst, "name", None)})

        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(folder, f"instances_{pmb_type}.csv"), index=False)

    # -----------------------------
    # REACTIONS
    # -----------------------------
    rows = []
    for rx in db.reactions.values():
        rows.append({
            "name": rx.name,
            "participants": _encode([p.model_dump() for p in rx.participants]),
            "pK": rx.pK if hasattr(rx, "pK") else None,
            "reaction_type": rx.reaction_type,
            "metadata": _encode(rx.metadata) if getattr(rx, "metadata", None) is not None else "",
        })
    pd.DataFrame(rows).to_csv(os.path.join(folder, "reactions.csv"), index=False)


# ----------------------------------------------------------------------
# LOAD
# ----------------------------------------------------------------------
def load_database_csv(db: _DFManagement, folder: str):
    """
    Load CSV files from folder into the provided _DFManagement instance.
    This mutates db.templates, db.instances and db.reactions in place.

    Important:
      - The PintQuantity.from_dict(...) returns a stored PintQuantity dataclass
        (no registry required). Reconstruction to a pint.Quantity occurs when
        users call .to_quantity(db.ureg) on the PintQuantity.
    """
    folder = Path(folder)
    if not folder.exists():
        raise FileNotFoundError(f"Folder '{folder}' does not exist.")

    # target pmb types we support
    pyMBE_types = ["particle", 
                   "residue", 
                   "molecule", 
                   "bond", 
                   "peptide",
                   "protein",
                   "hydrogel"]

    # -----------------------------
    # TEMPLATES
    # -----------------------------
    for pmb_type in pyMBE_types:
        csv_file = folder / f"templates_{pmb_type}.csv"
        if not csv_file.exists():
            continue
        df = pd.read_csv(csv_file, dtype=str).fillna("")

        templates: Dict[str, Any] = {}

        for _, row in df.iterrows():
            # row values are strings (or empty string)
            if pmb_type == "particle":
                sigma_d = _decode(row["sigma"])
                epsilon_d = _decode(row["epsilon"])
                cutoff_d = _decode(row["cutoff"])
                offset_d = _decode(row["offset"])
                states_d = _decode(row["states"])

                sigma = PintQuantity.from_dict(sigma_d) if sigma_d is not None else None
                epsilon = PintQuantity.from_dict(epsilon_d) if epsilon_d is not None else None
                cutoff = PintQuantity.from_dict(cutoff_d) if cutoff_d is not None else None
                offset = PintQuantity.from_dict(offset_d) if offset_d is not None else None

                states: Dict[str, ParticleState] = {}
                if isinstance(states_d, dict):
                    for sname, sdata in states_d.items():
                        # sdata expected to be a dict matching ParticleState fields
                        states[sname] = ParticleState(**sdata)

                tpl = ParticleTemplate(
                    name=row["name"],
                    sigma=sigma,
                    epsilon=epsilon,
                    cutoff=cutoff,
                    offset=offset,
                    states=states,
                )
                templates[tpl.name] = tpl

            elif pmb_type == "residue":
                sc = _decode(row.get("side_chains", "")) or []
                if not isinstance(sc, list):
                    sc = list(sc)
                tpl = ResidueTemplate(
                    name=row["name"],
                    central_bead=row.get("central_bead", ""),
                    side_chains=sc
                )
                templates[tpl.name] = tpl

            elif pmb_type == "molecule":
                rl = _decode(row.get("residue_list", "")) or []
                if not isinstance(rl, list):
                    rl = list(rl)
                tpl = MoleculeTemplate(
                    name=row["name"],
                    residue_list=rl
                )
                templates[tpl.name] = tpl
            elif pmb_type == "peptide":
                rl = _decode(row.get("residue_list", "")) or []
                if not isinstance(rl, list):
                    rl = list(rl)
                seq = _decode(row.get("sequence", "")) or []
                if not isinstance(seq, list):
                    seq = list(seq)
                tpl = PeptideTemplate(
                    name=row["name"],
                    model=row.get("model", ""),
                    residue_list=rl,
                    sequence=seq
                )
                templates[tpl.name] = tpl
            elif pmb_type == "protein":
                rl = _decode(row.get("residue_list", "")) or []
                if not isinstance(rl, list):
                    rl = list(rl)
                seq = _decode(row.get("sequence", "")) or []
                if not isinstance(seq, list):
                    seq = list(seq)
                tpl = ProteinTemplate(
                    name=row["name"],
                    model=row.get("model", ""),
                    residue_list=rl,
                    sequence=seq
                )
                templates[tpl.name] = tpl
            elif pmb_type == "bond":
                params_raw = _decode(row.get("parameters", "")) or {}
                parameters: Dict[str, Any] = {}
                for k, v in params_raw.items():
                    # if v is a dict, assume PintQuantity dict
                    if isinstance(v, dict) and {"magnitude", "units", "dimension"}.issubset(v.keys()):
                        parameters[k] = PintQuantity.from_dict(v)
                    else:
                        parameters[k] = v
                tpl = BondTemplate(
                    name=row["name"],
                    bond_type=row.get("bond_type", ""),
                    parameters=parameters,
                    l0=PintQuantity.from_dict(_decode(row["l0"])))
                templates[tpl.name] = tpl
            elif pmb_type == "hydrogel":
                node_map_raw = _decode(row.get("node_map", "")) or []
                chain_map_raw = _decode(row.get("chain_map", "")) or []
                
                node_map = [HydrogelNode(**n) for n in node_map_raw if isinstance(n, dict)]
                chain_map = [HydrogelChain(**c) for c in chain_map_raw if isinstance(c, dict)]
                tpl = HydrogelTemplate(
                    name=row["name"],
                    node_map=node_map,
                    chain_map=chain_map
                )
                templates[tpl.name] = tpl
        db.templates[pmb_type] = templates

    # -----------------------------
    # INSTANCES
    # -----------------------------
    for pmb_type in pyMBE_types:
        csv_file = folder / f"instances_{pmb_type}.csv"
        if not csv_file.exists():
            continue
        df = pd.read_csv(csv_file, dtype=str).fillna("")

        instances: Dict[Any, Any] = {}

        for _, row in df.iterrows():
            if pmb_type == "particle":
                # some fields might be empty strings -> map to None
                residue_val = row.get("residue_id", "") or ""
                molecule_val = row.get("molecule_id", "") or ""
                inst = ParticleInstance(
                    name=row["name"],
                    particle_id=int(row["particle_id"]),
                    initial_state=row["initial_state"],
                    residue_id=None if residue_val == "" else int(residue_val),
                    molecule_id=None if molecule_val == "" else int(molecule_val),
                )
                instances[inst.particle_id] = inst

            elif pmb_type == "residue":
                mol_val = row.get("molecule_id", "") or ""
                inst = ResidueInstance(
                    name=row["name"],
                    residue_id=int(row["residue_id"]),
                    molecule_id=None if mol_val == "" else int(mol_val),
                )
                instances[inst.residue_id] = inst

            elif pmb_type == "molecule":
                inst = MoleculeInstance(
                    name=row["name"],
                    molecule_id=int(row["molecule_id"]),
                )
                instances[inst.molecule_id] = inst
            elif pmb_type == "peptide":
                inst = PeptideInstance(
                    name=row["name"],
                    molecule_id=int(row["molecule_id"]),
                )
                instances[inst.molecule_id] = inst
            elif pmb_type == "protein":
                inst = ProteinInstance(
                    name=row["name"],
                    molecule_id=int(row["molecule_id"]),
                )
                instances[inst.molecule_id] = inst
            elif pmb_type == "bond":
                inst = BondInstance(
                    name=row["name"],
                    bond_id=int(row["bond_id"]),
                    particle_id1=int(row["particle_id1"]),
                    particle_id2=int(row["particle_id2"]),
                )
                instances[inst.bond_id] = inst
            elif pmb_type == "hydrogel":
                mol_ids = _decode(row.get("molecule_ids", "")) or []
                if not isinstance(mol_ids, list):
                    mol_ids = list(mol_ids)
                inst = HydrogelInstance(
                    name=row["name"],
                    hydrogel_id=int(row["hydrogel_id"]),
                    molecule_ids=mol_ids
                )
                instances[inst.hydrogel_id] = inst
        db.instances[pmb_type] = instances

    # -----------------------------
    # REACTIONS
    # -----------------------------
    rx_file = folder / "reactions.csv"
    reactions: Dict[str, Reaction] = {}
    if rx_file.exists():
        df = pd.read_csv(rx_file, dtype=str).fillna("")
        for _, row in df.iterrows():
            participants_raw = _decode(row.get("participants", "")) or []
            participants = [ReactionParticipant(**p) for p in participants_raw]
            metadata = _decode(row.get("metadata", "")) or None
            rx = Reaction(
                name=row["name"],
                participants=participants,
                pK=float(row["pK"]) if (row.get("pK") not in (None, "", "nan")) else None,
                reaction_type=row.get("reaction_type", None),
                metadata=metadata
            )
            reactions[rx.name] = rx
    db.reactions = reactions


