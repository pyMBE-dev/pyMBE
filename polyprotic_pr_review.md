# Polyprotic Acid/Base Support — PR Review Document

All changes are in `pyMBE/pyMBE.py` (modified) and `testsuite/polyprotic_acidbase_tests.py` (new file).
The `define_particle` function is **unchanged** — monoprotic behavior is fully preserved.

Base commit: `3147a11` (pre-polyprotic, `ganesh` branch)

---

## Table of Contents

1. [Modified: `_check_pka_set`](#1-modified-_check_pka_set) — line 174
2. [Modified: `calculate_HH`](#2-modified-calculate_hh) — line 503
3. [New: `define_polyprotic_acidbase_reactions`](#3-new-define_polyprotic_acidbase_reactions) — line 1567
4. [New: `define_polyprotic_particle_states`](#4-new-define_polyprotic_particle_states) — line 1638
5. [New: `define_polyprotic_particle`](#5-new-define_polyprotic_particle) — line 1746
6. [Modified: `get_pka_set`](#6-modified-get_pka_set) — line 2345
7. [Modified: `load_pka_set`](#7-modified-load_pka_set) — line 2504
8. [Modified: Reaction filters](#8-modified-reaction-filters) — lines 2786, 3102, 3293
9. [New: `testsuite/polyprotic_acidbase_tests.py`](#9-new-test-file)

---

## 1. Modified: `_check_pka_set`

**File:** `pyMBE/pyMBE.py:174`
**What changed:** Extended to accept either `"pka_value"` (monoprotic) or `"pka_values"` (polyprotic), rejecting both or neither.

```diff
     def _check_pka_set(self, pka_set):
         """
         Checks that 'pka_set' has the formatting expected by pyMBE.
-
+
         Args:
-            pka_set ('dict'):
-                {"name" : {"pka_value": pka, "acidity": acidity}}
-        """
-        required_keys=['pka_value','acidity']
-        for required_key in required_keys:
-            for pka_name, pka_entry in pka_set.items():
-                if required_key not in pka_entry:
-                    raise ValueError(f'missing a required key "{required_key}" in entry "{pka_name}" of pka_set ("{pka_entry}")')
+            pka_set ('dict'):
+                Monoprotic: {"name" : {"pka_value": pka, "acidity": acidity}}
+                Polyprotic: {"name" : {"pka_values": [pka1, pka2, ...], "acidity": acidity}}
+        """
+        for pka_name, pka_entry in pka_set.items():
+            if "acidity" not in pka_entry:
+                raise ValueError(f'missing required key "acidity" in entry "{pka_name}" of pka_set ("{pka_entry}")')
+            has_mono = "pka_value" in pka_entry
+            has_poly = "pka_values" in pka_entry
+            if not has_mono and not has_poly:
+                raise ValueError(f'missing "pka_value" or "pka_values" in entry "{pka_name}" of pka_set ("{pka_entry}")')
+            if has_mono and has_poly:
+                raise ValueError(f'entry "{pka_name}" has both "pka_value" and "pka_values", use only one')
         return
```

**Why:** `_check_pka_set` is called by `calculate_HH` and `load_pka_set`. It needs to validate both formats since polyprotic particles use `"pka_values": [list]` while monoprotic uses `"pka_value": float`.

---

## 2. Modified: `calculate_HH`

**File:** `pyMBE/pyMBE.py:503`
**What changed:** Added polyprotic branch using the coupled Henderson-Hasselbalch formula with cumulative pKa sums. Monoprotic formula untouched.

```diff
     def calculate_HH(self, template_name, pH_list=None, pka_set=None):
         """
         Calculates the charge in the template object according to the ideal  Henderson–Hasselbalch titration curve.

         Args:
             template_name ('str'):
                 Name of the template.

             pH_list ('list[float]', optional):
                 pH values at which the charge is evaluated.
                 Defaults to 50 values between 2 and 12.

             pka_set ('dict', optional):
-                Mapping: {particle_name: {"pka_value": 'float', "acidity": "acidic"|"basic"}}
+                Mapping per particle. Monoprotic entries use
+                {"pka_value": float, "acidity": ...}, polyprotic entries use
+                {"pka_values": [float, ...], "acidity": ...}.

         Returns:
             'list[float]':
                 Net molecular charge at each pH value.
         """
         if pH_list is None:
             pH_list = np.linspace(2, 12, 50)
         if pka_set is None:
             pka_set = self.get_pka_set()
         self._check_pka_set(pka_set=pka_set)
         particle_counts = self.db.get_particle_templates_under(template_name=template_name,
                                                                return_counts=True)
         if not particle_counts:
             return [None] * len(pH_list)
         charge_number_map = self.get_charge_number_map()
         def formal_charge(particle_name):
-            tpl = self.db.get_template(name=particle_name,
+            tpl = self.db.get_template(name=particle_name,
                                        pmb_type="particle")
             state = self.db.get_template(name=tpl.initial_state,
                                          pmb_type="particle_state")
             return charge_number_map[state.es_type]
         Z_HH = []
         for pH in pH_list:
             Z = 0.0
             for particle, multiplicity in particle_counts.items():
                 if particle in pka_set:
-                    pka = pka_set[particle]["pka_value"]
-                    acidity = pka_set[particle]["acidity"]
+                    entry = pka_set[particle]
+                    acidity = entry["acidity"]
                     if acidity == "acidic":
                         psi = -1
                     elif acidity == "basic":
                         psi = +1
                     else:
                         raise ValueError(f"Unknown acidity '{acidity}' for particle '{particle}'")
-                    charge = psi / (1.0 + 10.0 ** (psi * (pH - pka)))
+                    if "pka_values" in entry:
+                        pka_list = entry["pka_values"]
+                        n = len(pka_list)
+                        cumsum_pK = np.cumsum(pka_list)
+                        terms = [10.0 ** (-cumsum_pK[i] + (i + 1) * pH) for i in range(n)]
+                        denominator = 1.0 + sum(terms)
+                        numerator = sum((i + 1) * terms[i] for i in range(n))
+                        if acidity == "acidic":
+                            charge = -numerator / denominator
+                        else:
+                            charge = n - numerator / denominator
+                    else:
+                        pka = entry["pka_value"]
+                        charge = psi / (1.0 + 10.0 ** (psi * (pH - pka)))
                     Z += multiplicity * charge
                 else:
                     Z += multiplicity * formal_charge(particle)
             Z_HH.append(Z)
-        return Z_HH
+        return Z_HH
```

**Why:** A naive sum-of-independent-terms approach (treating each pKa independently) is incorrect for polyprotic systems. The correct coupled formula uses cumulative pKa sums and a shared partition function denominator:

- **Acidic:** `Z = -sum(i * 10^(-cumsum_pK[i] + i*pH)) / (1 + sum(10^(-cumsum_pK[i] + i*pH)))`
- **Basic:** `Z = n - sum(i * 10^(-cumsum_pK[i] + i*pH)) / (1 + sum(10^(-cumsum_pK[i] + i*pH)))`

This formula comes from `~/ntnu/pbs_simulations/scripts/post_processing.py:calculate_Z` and has been validated against both the analytical reference and ideal constant-pH simulation data.

**Key design decision:** The monoprotic formula `psi / (1 + 10^(psi*(pH-pka)))` is kept as-is for `"pka_value"` entries. Keeping them separate avoids any risk of breaking existing behavior.

---

## 3. New: `define_polyprotic_acidbase_reactions`

**File:** `pyMBE/pyMBE.py:1567`
**Placed after:** `define_monoprototic_acidbase_reaction` (line 1520)

```diff
+    def define_polyprotic_acidbase_reactions(self, particle_name, state_names, pka_list, acidity, metadata=None):
+        """
+        Defines stepwise acid-base reactions for a polyprotic particle.
+
+        Creates N reactions from N+1 ordered states, one per adjacent pair.
+
+        Args:
+            particle_name ('str'):
+                Unique label that identifies the particle template.
+
+            state_names ('list[str]'):
+                Ordered list of state names from most protonated to least,
+                as returned by define_polyprotic_particle_states.
+
+            pka_list ('list[float]'):
+                pKa values for each deprotonation step. Must have len(state_names) - 1 entries.
+
+            acidity ('str'):
+                Identifies whether the particle is 'acidic' or 'basic'.
+
+            metadata ('dict', optional):
+                Additional information to be stored in each reaction. Defaults to None.
+        """
+        if len(pka_list) != len(state_names) - 1:
+            raise ValueError(f"Expected {len(state_names) - 1} pKa values for {len(state_names)} states, got {len(pka_list)}.")
+        supported_acidities = ["acidic", "basic"]
+        if acidity not in supported_acidities:
+            raise ValueError(f"Unsupported acidity '{acidity}' for particle '{particle_name}'. Supported acidities are {supported_acidities}.")
+        if acidity == "basic":
+            reaction_type = "polyprotic_base"
+        else:
+            reaction_type = "polyprotic_acid"
+        for i, pka in enumerate(pka_list):
+            reactant_state = state_names[i]
+            product_state = state_names[i + 1]
+            reaction = Reaction(participants=[ReactionParticipant(particle_name=particle_name,
+                                                                  state_name=reactant_state,
+                                                                  coefficient=-1),
+                                              ReactionParticipant(particle_name=particle_name,
+                                                                  state_name=product_state,
+                                                                  coefficient=1)],
+                                reaction_type=reaction_type,
+                                pK=pka,
+                                metadata=metadata)
+            self.db._register_reaction(reaction)
```

**Why:** Each deprotonation step is a separate reaction with the same 2-participant structure as monoprotic. Uses `reaction_type="polyprotic_acid"` / `"polyprotic_base"` to distinguish from monoprotic — this is critical for `get_pka_set` which raises errors on duplicate monoprotic entries but accumulates polyprotic ones into a list.

**Example:** For triprotic acid "A" with pKa = [2.15, 7.20, 12.35]:
- Reaction 1: H3A -> H2A (pK=2.15)
- Reaction 2: H2A -> HA  (pK=7.20)
- Reaction 3: HA  -> A   (pK=12.35)

---

## 4. New: `define_polyprotic_particle_states`

**File:** `pyMBE/pyMBE.py:1638`
**Placed after:** `define_monoprototic_particle_states` (line 1613)

```diff
+    def define_polyprotic_particle_states(self, particle_name, n, acidity):
+        """
+        Defines particle states for a polyprotic particle template.
+
+        Generates N+1 states with auto-naming convention:
+            - Acidic with n=3, particle "A": H3A (z=0), H2A (z=-1), HA (z=-2), A (z=-3)
+            - Basic with n=2, particle "B":  H2B (z=+2), HB (z=+1), B (z=0)
+
+        Args:
+            particle_name ('str'):
+                Unique label that identifies the particle template.
+
+            n ('int'):
+                Number of protons that can dissociate (e.g. 2 for diprotic, 3 for triprotic).
+                Must be >= 2. For n=1, use define_monoprototic_particle_states instead.
+
+            acidity ('str'):
+                Identifies whether the particle is 'acidic' or 'basic'.
+
+        Returns:
+            ('list[str]'):
+                Ordered list of generated state names, from most protonated to least.
+        """
+        if n < 2:
+            raise ValueError(f"n={n} is not polyprotic. Use define_monoprototic_particle_states for n=1.")
+        acidity_valid_keys = ['acidic', 'basic']
+        if acidity not in acidity_valid_keys:
+            raise ValueError(f"Acidity '{acidity}' for particle '{particle_name}' is not supported. Valid keys are: {acidity_valid_keys}")
+        states = []
+        for k in range(n, -1, -1):
+            # Build name: H3A, H2A, HA, A
+            if k == 0:
+                name = particle_name
+            elif k == 1:
+                name = f"H{particle_name}"
+            else:
+                name = f"H{k}{particle_name}"
+            if acidity == "acidic":
+                z = -(n - k)
+            else:
+                z = k
+            states.append({"name": name, "z": z})
+        self.define_particle_states(particle_name=particle_name, states=states)
+        return [s["name"] for s in states]
```

**Auto-naming convention table:**

| k (protons) | Name format    | Example (n=3, "A") | Charge (acidic) | Charge (basic) |
|:-----------:|:--------------:|:-------------------:|:----------------:|:--------------:|
| n           | `H{n}{name}`   | `H3A`              | 0                | +3             |
| n-1         | `H{n-1}{name}` | `H2A`              | -1               | +2             |
| 1           | `H{name}`      | `HA`               | -(n-1)           | +1             |
| 0           | `{name}`       | `A`                | -n               | 0              |

---

## 5. New: `define_polyprotic_particle`

**File:** `pyMBE/pyMBE.py:1746`
**Placed after:** `define_particle` (line 1683)

```diff
+    def define_polyprotic_particle(self, name, sigma, epsilon, n, acidity, pka_list, cutoff=pd.NA, offset=pd.NA):
+        """
+        Defines a polyprotic particle template in the pyMBE database.
+
+        Creates N+1 states with auto-naming and N stepwise acid-base reactions.
+
+        Args:
+            name('str'):
+                 Unique label that identifies this particle type.
+
+            sigma('pint.Quantity'):
+                Sigma parameter used to set up Lennard-Jones interactions for this particle type.
+
+            epsilon('pint.Quantity'):
+                Epsilon parameter used to setup Lennard-Jones interactions for this particle tipe.
+
+            n('int'):
+                Number of dissociable protons (e.g. 2 for diprotic, 3 for triprotic).
+                Must be >= 2.
+
+            acidity('str'):
+                Identifies whether the particle is 'acidic' or 'basic'.
+
+            pka_list('list[float]'):
+                pKa values for each deprotonation step, ordered from first to last.
+                Must have exactly n entries.
+
+            cutoff('pint.Quantity', optional):
+                Cutoff parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
+
+            offset('pint.Quantity', optional):
+                Offset parameter used to set up Lennard-Jones interactions for this particle type. Defaults to pd.NA.
+
+        Notes:
+            - Auto-naming convention: H3A, H2A, HA, A for triprotic acid 'A'.
+            - 'sigma', 'cutoff' and 'offset' must have a dimensitonality of '[length]' and should be defined using pmb.units.
+            - 'epsilon' must have a dimensitonality of '[energy]' and should be defined using pmb.units.
+            - 'cutoff' defaults to '2**(1./6.) reduced_length'.
+            - 'offset' defaults to 0.
+        """
+        if len(pka_list) != n:
+            raise ValueError(f"Expected {n} pKa values for {n}-protic particle, got {len(pka_list)}.")
+        if pd.isna(cutoff):
+            cutoff=self.units.Quantity(2**(1./6.), "reduced_length")
+        if pd.isna(offset):
+            offset=self.units.Quantity(0, "reduced_length")
+        state_names = self.define_polyprotic_particle_states(particle_name=name,
+                                                             n=n,
+                                                             acidity=acidity)
+        initial_state = state_names[0]
+        self.define_polyprotic_acidbase_reactions(particle_name=name,
+                                                  state_names=state_names,
+                                                  pka_list=pka_list,
+                                                  acidity=acidity)
+        tpl = ParticleTemplate(name=name,
+                               sigma=PintQuantity.from_quantity(q=sigma, expected_dimension="length", ureg=self.units),
+                               epsilon=PintQuantity.from_quantity(q=epsilon, expected_dimension="energy", ureg=self.units),
+                               cutoff=PintQuantity.from_quantity(q=cutoff, expected_dimension="length", ureg=self.units),
+                               offset=PintQuantity.from_quantity(q=offset, expected_dimension="length", ureg=self.units),
+                               initial_state=initial_state)
+        self.db._register_template(tpl)
```

**Why a separate function instead of extending `define_particle`:** Explicitly requested to avoid any risk of breaking existing monoprotic `define_particle` behavior. `define_polyprotic_particle` has explicit `n` and `pka_list` parameters — no ambiguous routing based on whether `pka` is a scalar or list.

**What it does in sequence:**
1. Validates `len(pka_list) == n`
2. Calls `define_polyprotic_particle_states` -> creates N+1 states, returns ordered state names
3. Calls `define_polyprotic_acidbase_reactions` -> creates N stepwise reactions
4. Registers the `ParticleTemplate` with `initial_state` = most protonated state

---

## 6. Modified: `get_pka_set`

**File:** `pyMBE/pyMBE.py:2345`

```diff
     def get_pka_set(self):
         """
         Retrieve the pKa set for all titratable particles in the pyMBE database.

         Returns:
-            ('dict'):
+            ('dict'):
                 Dictionary of the form:
-                {"particle_name": {"pka_value": float,
-                                   "acidity": "acidic" | "basic"}}
-        Notes:
-            - If a particle participates in multiple acid/base reactions, an error is raised.
+                For monoprotic particles:
+                    {"particle_name": {"pka_value": float,
+                                       "acidity": "acidic" | "basic"}}
+                For polyprotic particles:
+                    {"particle_name": {"pka_values": [float, ...],
+                                       "acidity": "acidic" | "basic"}}
         """
         pka_set = {}
-        supported_reactions = ["monoprotic_acid",
-                               "monoprotic_base"]
+        supported_reactions = ["monoprotic_acid", "monoprotic_base",
+                               "polyprotic_acid", "polyprotic_base"]
         for reaction in self.db._reactions.values():
             if reaction.reaction_type not in supported_reactions:
                 continue
-            # Identify involved particle(s)
             particle_names = {participant.particle_name for participant in reaction.participants}
             particle_name = particle_names.pop()
-            if particle_name in pka_set:
-                raise ValueError(f"Multiple acid/base reactions found for particle '{particle_name}'.")
-            pka_set[particle_name] = {"pka_value": reaction.pK}
-            if reaction.reaction_type == "monoprotic_acid":
-                acidity = "acidic"
-            elif reaction.reaction_type == "monoprotic_base":
-                acidity = "basic"
-            pka_set[particle_name]["acidity"] = acidity
+            if reaction.reaction_type in ["monoprotic_acid", "monoprotic_base"]:
+                if particle_name in pka_set:
+                    raise ValueError(f"Multiple acid/base reactions found for particle '{particle_name}'.")
+                acidity = "acidic" if reaction.reaction_type == "monoprotic_acid" else "basic"
+                pka_set[particle_name] = {"pka_value": reaction.pK, "acidity": acidity}
+            else:
+                acidity = "acidic" if reaction.reaction_type == "polyprotic_acid" else "basic"
+                if particle_name not in pka_set:
+                    pka_set[particle_name] = {"pka_values": [], "acidity": acidity}
+                pka_set[particle_name]["pka_values"].append(reaction.pK)
         return pka_set
```

**Why:** Monoprotic particles have exactly 1 reaction, so duplicates indicate an error. Polyprotic particles have N reactions for the same particle — the `reaction_type` field (`"polyprotic_acid"` vs `"monoprotic_acid"`) is what lets `get_pka_set` know to accumulate instead of raising.

**Return format difference:**
- Monoprotic: `{"A": {"pka_value": 4.0, "acidity": "acidic"}}`
- Polyprotic: `{"PO4": {"pka_values": [2.15, 7.20, 12.35], "acidity": "acidic"}}`

This format is consumed by `calculate_HH` (line 503) and `_check_pka_set` (line 174).

---

## 7. Modified: `load_pka_set`

**File:** `pyMBE/pyMBE.py:2504`

```diff
     def load_pka_set(self, filename):
         """
         Load a pKa set and attach chemical states and acid–base reactions
         to existing particle templates.

         Args:
-            filename ('str'):
+            filename ('str'):
                 Path to a JSON file containing the pKa set. Expected format:
-                {"metadata": {...},
-                  "data": {"A": {"acidity": "acidic", "pka_value": 4.5},
-                           "B": {"acidity": "basic",  "pka_value": 9.8}}}
+                Monoprotic:
+                    {"metadata": {...},
+                      "data": {"A": {"acidity": "acidic", "pka_value": 4.5},
+                               "B": {"acidity": "basic",  "pka_value": 9.8}}}
+                Polyprotic:
+                    {"metadata": {...},
+                      "data": {"PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]}}}

         Returns:
-            ('dict'):
+            ('dict'):
                 Dictionary with bibliographic metadata about the original work were the pKa set was determined.
-
-        Notes:
-            - This method is designed for monoprotic acids and bases only.
         """
         with open(filename, "r") as f:
             pka_data = json.load(f)
         pka_set = pka_data["data"]
         metadata = pka_data.get("metadata", {})
         self._check_pka_set(pka_set)
         for particle_name, entry in pka_set.items():
             acidity = entry["acidity"]
-            pka = entry["pka_value"]
-            self.define_monoprototic_acidbase_reaction(particle_name=particle_name,
-                                                       pka=pka,
-                                                       acidity=acidity,
-                                                       metadata=metadata)
+            if "pka_values" in entry:
+                pka_list = entry["pka_values"]
+                n = len(pka_list)
+                # Build state names using the same H{k}{name} convention
+                state_names = []
+                for k in range(n, -1, -1):
+                    if k == 0:
+                        state_names.append(particle_name)
+                    elif k == 1:
+                        state_names.append(f"H{particle_name}")
+                    else:
+                        state_names.append(f"H{k}{particle_name}")
+                self.define_polyprotic_acidbase_reactions(particle_name=particle_name,
+                                                         state_names=state_names,
+                                                         pka_list=pka_list,
+                                                         acidity=acidity,
+                                                         metadata=metadata)
+            else:
+                pka = entry["pka_value"]
+                self.define_monoprototic_acidbase_reaction(particle_name=particle_name,
+                                                           pka=pka,
+                                                           acidity=acidity,
+                                                           metadata=metadata)
         return metadata
```

**Why:** `load_pka_set` creates **reactions only** (not states), matching the existing monoprotic pattern. States are created separately by the user via `define_polyprotic_particle_states` or `define_polyprotic_particle`. The state names are reconstructed inline using the same `H{k}{name}` convention so the reaction participants reference the correct state names.

**New JSON format for polyprotic:**
```json
{
  "metadata": {"source": "CRC Handbook"},
  "data": {
    "PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]},
    "A":   {"acidity": "acidic", "pka_value": 4.5}
  }
}
```

---

## 8. Modified: Reaction filters

Three methods filter reactions to find acid-base ones for setting up ESPResSo reaction methods. Each had the same change: expanding the filter to include polyprotic reaction types.

### `setup_cpH` — `pyMBE/pyMBE.py:2786`

```diff
+        acidbase_reaction_types = ["monoprotic_acid", "monoprotic_base",
+                                    "polyprotic_acid", "polyprotic_base"]
         for reaction in self.db.get_reactions():
-            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
+            if reaction.reaction_type not in acidbase_reaction_types:
                 continue
```

### `setup_grxmc_reactions` — `pyMBE/pyMBE.py:3102`

```diff
+        acidbase_reaction_types = ["monoprotic_acid", "monoprotic_base",
+                                    "polyprotic_acid", "polyprotic_base"]
         for reaction in self.db.get_reactions():
-            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
+            if reaction.reaction_type not in acidbase_reaction_types:
                 continue
```

### `setup_grxmc_unified` — `pyMBE/pyMBE.py:3293`

```diff
+        acidbase_reaction_types = ["monoprotic_acid", "monoprotic_base",
+                                    "polyprotic_acid", "polyprotic_base"]
         for reaction in self.db.get_reactions():
-            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
+            if reaction.reaction_type not in acidbase_reaction_types:
                 continue
```

**Why:** Without this, polyprotic reactions would be silently skipped when setting up the simulation reaction ensemble, making the polyprotic particles inert in actual simulations.

---

## 9. New: Test file

**File:** `testsuite/polyprotic_acidbase_tests.py` (345 lines, 17 tests)

```diff
+#
+# Copyright (C) 2024-2026 pyMBE-dev team
+#
+# This file is part of pyMBE.
+#
+# pyMBE is free software: you can redistribute it and/or modify
+# it under the terms of the GNU General Public License as published by
+# the Free Software Foundation, either version 3 of the License, or
+# (at your option) any later version.
+#
+# pyMBE is distributed in the hope that it will be useful,
+# but WITHOUT ANY WARRANTY; without even the implied warranty of
+# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
+# GNU General Public License for more details.
+#
+# You should have received a copy of the GNU General Public License
+# along with this program.  If not, see <http://www.gnu.org/licenses/>.
+
+import json
+import os
+import tempfile
+import pyMBE
+import unittest as ut
+from pyMBE.storage.reactions.reaction import Reaction, ReactionParticipant
+
+
+class TestPolyprotic(ut.TestCase):
+
+    def test_triprotic_acid_states(self):
+        """
+        Test that a triprotic acid creates 4 states with correct names and charges.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="PO4",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.15, 7.20, 12.35])
+        states = pmb.db.get_particle_states_templates("PO4")
+        self.assertEqual(len(states), 4)
+        expected = {"H3PO4": 0, "H2PO4": -1, "HPO4": -2, "PO4": -3}
+        for name, z in expected.items():
+            self.assertIn(name, states)
+            self.assertEqual(states[name].z, z)
+        # All es_types must be unique
+        es_types = [s.es_type for s in states.values()]
+        self.assertEqual(len(es_types), len(set(es_types)))
+
+    def test_triprotic_acid_initial_state(self):
+        """
+        Test that the initial state is the most protonated one.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="PO4",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.15, 7.20, 12.35])
+        tpl = pmb.db.get_template("particle", "PO4")
+        self.assertEqual(tpl.initial_state, "H3PO4")
+
+    def test_triprotic_acid_reactions(self):
+        """
+        Test that a triprotic acid creates 3 stepwise reactions with correct pKa values.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="PO4",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.15, 7.20, 12.35])
+        reactions = list(pmb.db.get_reactions())
+        self.assertEqual(len(reactions), 3)
+        for r in reactions:
+            self.assertEqual(r.reaction_type, "polyprotic_acid")
+        pks = sorted([r.pK for r in reactions])
+        self.assertAlmostEqual(pks[0], 2.15)
+        self.assertAlmostEqual(pks[1], 7.20)
+        self.assertAlmostEqual(pks[2], 12.35)
+
+    def test_triprotic_acid_reaction_participants(self):
+        """
+        Test that each reaction has the correct reactant/product pair.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="A",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.0, 7.0, 12.0])
+        reactions = sorted(pmb.db.get_reactions(), key=lambda r: r.pK)
+        expected_pairs = [("H3A", "H2A"), ("H2A", "HA"), ("HA", "A")]
+        for r, (reactant, product) in zip(reactions, expected_pairs):
+            reactants = [p for p in r.participants if p.coefficient < 0]
+            products = [p for p in r.participants if p.coefficient > 0]
+            self.assertEqual(len(reactants), 1)
+            self.assertEqual(len(products), 1)
+            self.assertEqual(reactants[0].state_name, reactant)
+            self.assertEqual(products[0].state_name, product)
+
+    def test_diprotic_base_states(self):
+        """
+        Test that a diprotic base creates 3 states with correct charges.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="B",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=2,
+                                       acidity="basic",
+                                       pka_list=[6.0, 10.0])
+        states = pmb.db.get_particle_states_templates("B")
+        self.assertEqual(len(states), 3)
+        expected = {"H2B": 2, "HB": 1, "B": 0}
+        for name, z in expected.items():
+            self.assertIn(name, states)
+            self.assertEqual(states[name].z, z)
+
+    def test_diprotic_base_reactions(self):
+        """
+        Test that a diprotic base creates 2 reactions with type polyprotic_base.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="B",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=2,
+                                       acidity="basic",
+                                       pka_list=[6.0, 10.0])
+        reactions = list(pmb.db.get_reactions())
+        self.assertEqual(len(reactions), 2)
+        for r in reactions:
+            self.assertEqual(r.reaction_type, "polyprotic_base")
+
+    def test_monoprotic_backward_compat(self):
+        """
+        Test that a single pka float still works (monoprotic path).
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_particle(name="A",
+                            sigma=0.35*pmb.units.nm,
+                            epsilon=1*pmb.units("reduced_energy"),
+                            acidity="acidic",
+                            pka=4.0)
+        states = pmb.db.get_particle_states_templates("A")
+        self.assertEqual(len(states), 2)
+        self.assertIn("AH", states)
+        self.assertIn("A", states)
+        reactions = list(pmb.db.get_reactions())
+        self.assertEqual(len(reactions), 1)
+        self.assertEqual(reactions[0].reaction_type, "monoprotic_acid")
+
+    def test_get_pka_set_polyprotic(self):
+        """
+        Test that get_pka_set returns pka_values list for polyprotic particles.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="PO4",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.15, 7.20, 12.35])
+        pka_set = pmb.get_pka_set()
+        self.assertIn("PO4", pka_set)
+        self.assertIn("pka_values", pka_set["PO4"])
+        self.assertEqual(pka_set["PO4"]["acidity"], "acidic")
+        self.assertEqual(pka_set["PO4"]["pka_values"], [2.15, 7.20, 12.35])
+
+    def test_get_pka_set_mixed(self):
+        """
+        Test get_pka_set with both monoprotic and polyprotic particles.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_particle(name="A",
+                            sigma=0.35*pmb.units.nm,
+                            epsilon=1*pmb.units("reduced_energy"),
+                            acidity="acidic",
+                            pka=4.0)
+        pmb.define_polyprotic_particle(name="PO4",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.15, 7.20, 12.35])
+        pka_set = pmb.get_pka_set()
+        self.assertIn("pka_value", pka_set["A"])
+        self.assertIn("pka_values", pka_set["PO4"])
+
+    def test_load_pka_set_polyprotic(self):
+        """
+        Test that load_pka_set creates reactions only (not states), matching monoprotic behavior.
+        User creates states separately via define_polyprotic_particle_states.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pka_data = {
+            "metadata": {"source": "test"},
+            "data": {
+                "PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]}
+            }
+        }
+        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
+            json.dump(pka_data, f)
+            tmpfile = f.name
+        try:
+            metadata = pmb.load_pka_set(tmpfile)
+            self.assertEqual(metadata, {"source": "test"})
+            # Only reactions created, no states yet
+            reactions = list(pmb.db.get_reactions())
+            self.assertEqual(len(reactions), 3)
+            for r in reactions:
+                self.assertEqual(r.reaction_type, "polyprotic_acid")
+            # Now create states separately (as user would do after load_pka_set)
+            pmb.define_polyprotic_particle_states("PO4", n=3, acidity="acidic")
+            states = pmb.db.get_particle_states_templates("PO4")
+            self.assertEqual(len(states), 4)
+        finally:
+            os.unlink(tmpfile)
+
+    def test_load_pka_set_mixed(self):
+        """
+        Test that load_pka_set handles a mix of monoprotic and polyprotic in one file.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pka_data = {
+            "metadata": {},
+            "data": {
+                "A": {"acidity": "acidic", "pka_value": 4.5},
+                "PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]}
+            }
+        }
+        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
+            json.dump(pka_data, f)
+            tmpfile = f.name
+        try:
+            pmb.load_pka_set(tmpfile)
+            reactions = list(pmb.db.get_reactions())
+            mono_rxns = [r for r in reactions if r.reaction_type == "monoprotic_acid"]
+            poly_rxns = [r for r in reactions if r.reaction_type == "polyprotic_acid"]
+            self.assertEqual(len(mono_rxns), 1)
+            self.assertEqual(len(poly_rxns), 3)
+        finally:
+            os.unlink(tmpfile)
+
+    def test_check_pka_set_validation(self):
+        """
+        Test that _check_pka_set rejects invalid formats.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        # Missing both pka_value and pka_values
+        with self.assertRaises(ValueError):
+            pmb._check_pka_set({"X": {"acidity": "acidic"}})
+        # Missing acidity
+        with self.assertRaises(ValueError):
+            pmb._check_pka_set({"X": {"pka_value": 4.0}})
+        # Both pka_value and pka_values
+        with self.assertRaises(ValueError):
+            pmb._check_pka_set({"X": {"acidity": "acidic", "pka_value": 4.0, "pka_values": [1, 2]}})
+
+    def test_define_polyprotic_particle_states_errors(self):
+        """
+        Test error handling in define_polyprotic_particle_states.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        # n < 2
+        with self.assertRaises(ValueError):
+            pmb.define_polyprotic_particle_states("X", n=1, acidity="acidic")
+        # Invalid acidity
+        with self.assertRaises(ValueError):
+            pmb.define_polyprotic_particle_states("X", n=2, acidity="neutral")
+
+    def test_define_polyprotic_reactions_pka_count_mismatch(self):
+        """
+        Test that mismatched pka_list length raises ValueError.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        state_names = ["H3A", "H2A", "HA", "A"]
+        # 4 states need 3 pKa values, not 2
+        with self.assertRaises(ValueError):
+            pmb.define_polyprotic_acidbase_reactions(
+                particle_name="A",
+                state_names=state_names,
+                pka_list=[2.0, 7.0],
+                acidity="acidic")
+
+    def test_auto_naming_convention(self):
+        """
+        Test the H{k}{name} naming convention for various n values.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        # n=2 (diprotic): H2X, HX, X
+        names = pmb.define_polyprotic_particle_states("X", n=2, acidity="acidic")
+        self.assertEqual(names, ["H2X", "HX", "X"])
+        # n=4 (tetraprotic): H4Y, H3Y, H2Y, HY, Y
+        names = pmb.define_polyprotic_particle_states("Y", n=4, acidity="acidic")
+        self.assertEqual(names, ["H4Y", "H3Y", "H2Y", "HY", "Y"])
+
+
+    def test_calculate_HH_triprotic_acid(self):
+        """
+        Test Henderson-Hasselbalch charge for a triprotic acid at extreme pH values.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="PO4",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=3,
+                                       acidity="acidic",
+                                       pka_list=[2.15, 7.20, 12.35])
+        pmb.define_molecule(name="mol", residue_list=[])
+        pmb.define_residue(name="res", central_bead="PO4", side_chains=[])
+        pmb.define_molecule(name="test_mol", residue_list=["res"])
+        Z = pmb.calculate_HH(template_name="test_mol", pH_list=[0, 14])
+        # At pH 0: fully protonated, charge ~ 0
+        self.assertAlmostEqual(Z[0], 0.0, places=1)
+        # At pH 14: fully deprotonated, charge ~ -3
+        self.assertAlmostEqual(Z[1], -3.0, places=1)
+
+    def test_calculate_HH_diprotic_base(self):
+        """
+        Test Henderson-Hasselbalch charge for a diprotic base at extreme pH values.
+        """
+        pmb = pyMBE.pymbe_library(seed=42)
+        pmb.define_polyprotic_particle(name="B",
+                                       sigma=0.35*pmb.units.nm,
+                                       epsilon=1*pmb.units("reduced_energy"),
+                                       n=2,
+                                       acidity="basic",
+                                       pka_list=[6.0, 10.0])
+        pmb.define_residue(name="res", central_bead="B", side_chains=[])
+        pmb.define_molecule(name="test_mol", residue_list=["res"])
+        Z = pmb.calculate_HH(template_name="test_mol", pH_list=[0, 14])
+        # At pH 0: fully protonated, charge ~ +2
+        self.assertAlmostEqual(Z[0], 2.0, places=2)
+        # At pH 14: fully deprotonated, charge ~ 0
+        self.assertAlmostEqual(Z[1], 0.0, places=2)
+
+
+if __name__ == "__main__":
+    ut.main()
```

### Test summary

| # | Test name | What it verifies |
|---|-----------|-----------------|
| 1 | `test_triprotic_acid_states` | 4 states with correct names (H3PO4, H2PO4, HPO4, PO4) and charges (0, -1, -2, -3). Unique es_types. |
| 2 | `test_triprotic_acid_initial_state` | Initial state is most protonated (H3PO4). |
| 3 | `test_triprotic_acid_reactions` | 3 reactions with type `polyprotic_acid` and correct pKa values. |
| 4 | `test_triprotic_acid_reaction_participants` | Each reaction pairs the correct adjacent states (H3A->H2A, H2A->HA, HA->A). |
| 5 | `test_diprotic_base_states` | 3 states with correct charges (H2B=+2, HB=+1, B=0). |
| 6 | `test_diprotic_base_reactions` | 2 reactions with type `polyprotic_base`. |
| 7 | `test_monoprotic_backward_compat` | `define_particle(pka=float)` still works: 2 states, 1 reaction, type `monoprotic_acid`. |
| 8 | `test_get_pka_set_polyprotic` | Returns `"pka_values": [list]` for polyprotic particle. |
| 9 | `test_get_pka_set_mixed` | Returns both `"pka_value"` and `"pka_values"` when mixing mono/polyprotic. |
| 10 | `test_load_pka_set_polyprotic` | JSON with `"pka_values"` creates 3 polyprotic reactions. States created separately. |
| 11 | `test_load_pka_set_mixed` | Mixed JSON creates correct mono (1) and polyprotic (3) reactions. |
| 12 | `test_check_pka_set_validation` | Rejects missing both/neither `pka_value`/`pka_values`, missing `acidity`, and having both keys. |
| 13 | `test_define_polyprotic_particle_states_errors` | n<2 raises ValueError. Invalid acidity raises ValueError. |
| 14 | `test_define_polyprotic_reactions_pka_count_mismatch` | Mismatched pka_list length raises ValueError. |
| 15 | `test_auto_naming_convention` | n=2 -> [H2X, HX, X]. n=4 -> [H4Y, H3Y, H2Y, HY, Y]. |
| 16 | `test_calculate_HH_triprotic_acid` | At pH=0: charge ~ 0. At pH=14: charge ~ -3. |
| 17 | `test_calculate_HH_diprotic_base` | At pH=0: charge ~ +2. At pH=14: charge ~ 0. |

---

## Usage example

```python
import pyMBE

pmb = pyMBE.pymbe_library(seed=42)

# Triprotic acid (e.g. phosphoric acid)
pmb.define_polyprotic_particle(
    name="PO4",
    sigma=0.35 * pmb.units.nm,
    epsilon=1 * pmb.units("reduced_energy"),
    n=3,
    acidity="acidic",
    pka_list=[2.15, 7.20, 12.35]
)
# Creates:
#   States:    H3PO4 (z=0), H2PO4 (z=-1), HPO4 (z=-2), PO4 (z=-3)
#   Reactions: H3PO4 -> H2PO4 (pK=2.15)
#              H2PO4 -> HPO4  (pK=7.20)
#              HPO4  -> PO4   (pK=12.35)

# Diprotic base
pmb.define_polyprotic_particle(
    name="B",
    sigma=0.35 * pmb.units.nm,
    epsilon=1 * pmb.units("reduced_energy"),
    n=2,
    acidity="basic",
    pka_list=[6.0, 10.0]
)
# Creates:
#   States:    H2B (z=+2), HB (z=+1), B (z=0)
#   Reactions: H2B -> HB (pK=6.0)
#              HB  -> B  (pK=10.0)

# Monoprotic still uses define_particle (unchanged)
pmb.define_particle(
    name="A",
    sigma=0.35 * pmb.units.nm,
    epsilon=1 * pmb.units("reduced_energy"),
    acidity="acidic",
    pka=4.0
)
```

---

## Files not part of PR (development/validation only)

- `polyprotic_changelog.md` — development log
- `polyprotic_pr_review.md` — this file
- `tests/compare_polyprotic_HH.py` — validation vs analytical + simulation data (acidic)
- `tests/compare_polyprotic_HH_base.py` — validation vs analytical (basic)
- `pyMBE/.pyMBE.py.swp` — vim swap file
- `tutorials/peptide.png` — unrelated
