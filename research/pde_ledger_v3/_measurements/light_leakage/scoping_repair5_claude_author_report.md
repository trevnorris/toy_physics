# Light-leakage inventory, repair 5: author report (fresh Claude author, 2026-10-07)

**Artifact:** `light_leakage_scoping.md`, edited in place from v4 (`214d1fe6…`). The frozen copy
`light_leakage_scoping_v4.md` is untouched and still hashes to `214d1fe6…`. No lines were added or removed
(161 lines). Three lines changed: 26, 27 and 41. A cell-level comparison against v4 (split on unescaped `|`)
shows one changed cell on each line: row 3 mapping (line 26), row 4 mapping (line 27), and the taxonomy
note (line 41). Every table row keeps its cell count (8 pieces).

## The finding holds, so it was applied

Mechanical lookups (literal output):

```
$ sed -n '125p' research/pde_ledger_v3/directives/S11c_decisions.md
inhomogeneity `η`**: the transverse↔thickness coupling is `O(εη)`; the leakage probability/rate `O(ε²η²)`.
$ sed -n '194,198p' research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
Wave perturbations `u`, `ζ_c`, `δW`, `θ`, and the first variations of the face and bulk fields carry the
independent amplitude bookkeeper `ε`. Every computed object must be multigraded by `(ε,η,σ_W)` from its
actual data dependency. No term is removed merely because it contains both a wave and a background
bookkeeper; the requested truncation is first order in wave amplitude and first shape order in each
background bookkeeper.
$ sed -n '366,368p' research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
⛔ Do **not** take the first variation of the perturbation (wave) energy and then set the perturbations to
zero: that energy is bilinear in the perturbations, so its first variation is linear in them and vanishes
identically at `𝔅⁰`, producing a vacuous operand with no data dependence on the background it is meant to
$ sed -n '48p' research/pde_ledger_v3/steps/S11c_SCOPE.md
- ⭐ **In S11c** (non-uniform, linear): the full **variable-coefficient** slab spectrum, actual **leakage
```

`δW` (thickness) is a wave perturbation (SPa:194). DEC:125 places the transverse↔thickness coupling at
`O(εη)`, which is first order in `ε`. So v4's clauses on lines 26 and 41 ("a coupling of transverse light to
another wave perturbation lies beyond the records' 'first order in wave amplitude' truncation"; "enters as a
product of wave amplitudes, beyond the first-order truncation") contradict the records. They also contradict
the inventory's unchanged E1 route (line 57, "nonuniform transverse↔`{θ,e_W,u_L}` operator, including
gradient-thickness leg").

## Changes

| Line | Cell | Change | Record support |
|---|---|---|---|
| 26 | Row 3, mapping (col 3) | Deleted ", and, on this inventory's reading, a coupling of transverse light to another wave perturbation lies beyond the records' “first order in wave amplitude” truncation (SPa:194–198; taxonomy note below)". The pointer to the note moves into the remaining citation: "(S11:157–159; SC:50–51; taxonomy note below)". The C3 pointer, the AB “linear” gloss, the three gaps and "Not O7" are unchanged. | Removal only. The remaining basis is v4's existing text, S11:157–159 and SC:50–51. |
| 41 | Rows 3–4 taxonomy note | Deleted ": on held profiles, a coupling of transverse light to another wave perturbation enters as a product of wave amplitudes, beyond the first-order truncation (this inventory's reading of SPa:194–198)". Added in its place: "This pointer does not claim that a coupling between wave branches lies, as such, beyond the “first order in wave amplitude” truncation: S11c is “non-uniform, linear” (SC:48), DEC:125 gives “the transverse↔thickness coupling is `O(εη)`”, and E1 (§2.1) carries that coupling as an equation-level route." The rest of the note is byte-identical, including the C3 pointer "whatever the initiation … or WB gain class", the two gaps, the WB thermal/vacuum caveat, "Neither row points to O7" and the closing no-verdict sentence. | SC:48 and DEC:125, both quoted verbatim. E1's gradient-to-`e_W` leg is the inventory's existing statement (line 54: b:54–57; SC:16–17). |

## Other cells changed for consistency (line 27, row 4 mapping, col 3)

- **Colon → period:** "…or WB gain class: the records carry…" now reads "…or WB gain class. The records carry…".
  In v4 the colon made the reading "a light–acoustic-type coupling is a coupling among wave perturbations,
  not a release of the background hold" the stated reason for routing to C3. The note's "beyond the
  first-order truncation" step is now gone. Without that step, the colon would make "being a coupling among
  wave perturbations" the reason for C3 on its own. That is the removed claim in implicit form, and it would
  contradict E1, which is also a coupling among wave perturbations but is an equation-level route, not C3.
  After the change, the reading supports only "not a release of the background hold", which is the not-O7
  point (note: "Neither row points to O7…"). The wording of the reading is otherwise unchanged.
- **"quadratic cross-block is zero" → "the homogeneous quadratic cross-block is zero":** the cited record limits
  the zero to the homogeneous slab: "the homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches
  have zero linear cross-block under the selected action. This does not establish decoupling at nonlinear
  order, on a nonuniform slab, …" (S11:74–76, already cited in the cell as S11:68–76). Without the qualifier, the
  sentence leading to the C3 pointer could be read as a zero first-order transverse↔longitudinal coupling in
  general. That reading would contradict E1's nonuniform transverse↔`u_L` operator and the note's new DEC:125
  sentence. This text dates from v1.

No count changes: the mapping accounting (line 39, NOT ADDRESSED 6) and the route accounting (line 75, 15
items) are unaffected. No status, vertex or route was added.

## Noticed but left alone

- **The basis for routing the spontaneous forms to C3.** With the over-broad clause gone, no record line
  checked here places *spontaneous* Raman or Brillouin in the nonlinear program. The records name C3 as
  "nonlinear intensity coupling" and the "DC/harmonic/sideband radiation audit" (S11:158; SC:50–51), and none
  names Raman or Brillouin. AB §4.2 calls spontaneous Raman the "linear polarisation" term. The C3 pointer for
  spontaneous initiation therefore rests on the r3 F1 routing decision, which the user accepted, not on a
  record statement. I added no justification, because the brief says the routing is not in question and bars
  new readings. The bracket author may want to know this.
- **Line 27's "on this inventory's reading"** (wave perturbations on held profiles, so not a release of the
  hold): this is v4 text, and both r4 legs found the not-O7 routing faithful. Kept. It is now attached only to
  the not-O7 point.
- **Line 41's "request a truncation “first order in wave amplitude” (SPa:194–198)":** an accurate record
  quote. It is kept. It no longer supports any exclusion, and the new sentence sits beside it.
- **S11c_b SHARED_PHYSICS:366–368** (wave energy "bilinear in the perturbations, so its first variation is
  linear in them"): the finding cites it, but I did not add it. Its context is the background-order
  admissibility operand evaluated at `𝔅⁰`, so using it as a general energy-versus-equation order rule would be
  my own reading. DEC:125 and SC:48 state the needed point directly.
- **DEC:125's second clause** ("the leakage probability/rate `O(ε²η²)`"): not quoted, because it is not needed
  for the repair.
- **Line 92 (WB oracle row, "No three-wave vertex is recorded")** and **line 124 (C3 matrix, "Beyond linear
  S11c; quadratic zero insufficient")**: consistent with the change. Byte-identical.
- **Checksum file:** `scoping_v4.sha256` now reports FAILED for the edited file. This is expected; the frozen
  copy still matches `214d1fe6…`.
