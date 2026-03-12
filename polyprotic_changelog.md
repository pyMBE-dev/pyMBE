# Polyprotic Acid/Base Support — Change Log

## Stage 1: Storage Layer
**Status:** No changes needed.

`define_particle_states` already accepts arbitrary lists of `{"name", "z"}` dicts.
`Reaction` and `ReactionParticipant` are fully generic. The database can hold N states and N-1 reactions as-is.

---

## Stage 2: `define_polyprotic_particle_states`
**Status:** Done.
**File:** `pyMBE/pyMBE.py` (inserted after `define_monoprototic_particle_states`)

New method generates N+1 states with auto-naming:
- Acidic n=3, particle "A": `H3A` (z=0), `H2A` (z=-1), `HA` (z=-2), `A` (z=-3)
- Basic n=2, particle "B": `H2B` (z=+2), `HB` (z=+1), `B` (z=0)

Returns list of state names for use in reaction definition.

```diff
+    def define_polyprotic_particle_states(self, particle_name, n, acidity):
+        ...
+        states = []
+        for k in range(n, -1, -1):
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

---

## Stage 3: `define_polyprotic_acidbase_reactions`
**Status:** Done.
**File:** `pyMBE/pyMBE.py` (inserted after `define_monoprototic_acidbase_reaction`)

Creates N stepwise reactions from N+1 ordered states. Each reaction pairs adjacent states
(e.g. `H3A → H2A`, `H2A → HA`, `HA → A`) with its corresponding pKa. Uses
`reaction_type="polyprotic_acid"` or `"polyprotic_base"`.

Takes `state_names` (from Stage 2) and `pka_list` as input, validates
`len(pka_list) == len(state_names) - 1`.

```diff
+    def define_polyprotic_acidbase_reactions(self, particle_name, state_names, pka_list, acidity, metadata=None):
+        ...
+        for i, pka in enumerate(pka_list):
+            reactant_state = state_names[i]
+            product_state = state_names[i + 1]
+            reaction = Reaction(participants=[
+                ReactionParticipant(..., state_name=reactant_state, coefficient=-1),
+                ReactionParticipant(..., state_name=product_state,  coefficient=1)],
+                reaction_type=reaction_type, pK=pka, metadata=metadata)
+            self.db._register_reaction(reaction)
```

---

## Stage 4: New `define_polyprotic_particle`
**Status:** Done.
**File:** `pyMBE/pyMBE.py` (new method, inserted after `define_particle`)

Separate function with explicit `n` and `pka_list` parameters.
`define_particle` is **untouched** — monoprotic behavior fully preserved.

```diff
+    def define_polyprotic_particle(self, name, sigma, epsilon, n, acidity, pka_list, cutoff=pd.NA, offset=pd.NA):
+        if len(pka_list) != n:
+            raise ValueError(...)
+        state_names = self.define_polyprotic_particle_states(particle_name=name, n=n, acidity=acidity)
+        initial_state = state_names[0]
+        self.define_polyprotic_acidbase_reactions(particle_name=name, state_names=state_names,
+                                                  pka_list=pka_list, acidity=acidity)
+        tpl = ParticleTemplate(name=name, sigma=..., epsilon=..., cutoff=..., offset=...,
+                               initial_state=initial_state)
+        self.db._register_template(tpl)
```

---

## Stage 5: Extend `setup_cpH` filter
**Status:** Done.
**File:** `pyMBE/pyMBE.py` line ~2701

```diff
+        acidbase_reaction_types = ["monoprotic_acid", "monoprotic_base",
+                                    "polyprotic_acid", "polyprotic_base"]
         for reaction in self.db.get_reactions():
-            if reaction.reaction_type not in ["monoprotic_acid", "monoprotic_base"]:
+            if reaction.reaction_type not in acidbase_reaction_types:
```

---

## Stage 6: Extend `setup_grxmc_reactions` and `setup_grxmc_unified` filters
**Status:** Done.
**File:** `pyMBE/pyMBE.py` lines ~3017, ~3208

Same change as Stage 5 applied to both methods.

---

## Stage 7: Rework `get_pka_set`
**Status:** Done.
**File:** `pyMBE/pyMBE.py` (modified `get_pka_set`)

- Now returns `"pka_values": [list]` for polyprotic (vs `"pka_value": float` for monoprotic)
- Polyprotic particles accumulate pKa values across multiple reactions
- Monoprotic still raises on duplicates

```diff
+            if reaction.reaction_type in ["monoprotic_acid", "monoprotic_base"]:
+                # single pka_value (unchanged behavior)
+            else:
+                # accumulate into pka_values list
```

---

## Stage 8: Extend `load_pka_set` and `_check_pka_set`
**Status:** Done.
**File:** `pyMBE/pyMBE.py`

`_check_pka_set` now accepts either `"pka_value"` (float) or `"pka_values"` (list), rejects both/neither.

`load_pka_set` routes by key:
- `"pka_values"` present → calls `define_polyprotic_particle_states` + `define_polyprotic_acidbase_reactions`
- `"pka_value"` present → calls `define_monoprototic_acidbase_reaction` (unchanged)

New JSON format for polyprotic:
```json
{"data": {"PO4": {"acidity": "acidic", "pka_values": [2.15, 7.20, 12.35]}}}
```

---

## Tests
**File:** `testsuite/polyprotic_acidbase_tests.py` (new, 15 tests)

Covers: triprotic acid states/reactions/initial_state, diprotic base, monoprotic backward compat,
reaction participant pairing, auto-naming convention, get_pka_set (polyprotic/mixed),
load_pka_set (polyprotic/mixed), _check_pka_set validation, error handling.

All polyprotic tests use `define_polyprotic_particle(n=..., pka_list=[...])`.
Monoprotic backward compat test uses unchanged `define_particle(pka=float)`.

All 24 tests pass (15 new + 9 existing acidity).

---

## Tutorial
**File:** `tutorials/pyMBE_tutorial.ipynb`

Added new section "How to define polyprotic particles" between the polyampholyte exercise
and the peptides section. Contains:
- Explanation of polyprotic vs monoprotic
- Triprotic phosphoric acid example using `define_polyprotic_particle`
- Side-by-side comparison with monoprotic `define_particle`
- Auto-naming convention table
