# DIRECTIVE — S11b-A Mathematica audit, REPAIR PASS

**Target (edit in place):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bA_interface_response_mathematica_audit.wl`

Run with `math -script <path>`. Iterate until it completes with no errors and no unevaluated output. Then
stop and exit — ⛔ do not write a report or a second script.

⛔ **The read-bar list of the shared specification still applies in full**, including that no `.py` and no
other engine's work may be read. This script is still blind.

---

## ⛔⛔ WHAT THIS REPAIR IS, IN ONE SENTENCE

Two independent reviews ran **27 ablations** against this script. **Eight of its checks fired correctly.**
The rest of its conclusions are emitted inside association objects that **no check ever reads**, so
`VERDICT: PASS` survives corruption of things it appears to certify.

⇒ ⭐ **This is ONE defect with many instances, ⛔ not many defects.** Repair the class.

## ⭐⭐ THE STRUCTURAL REQUIREMENT

> **Every emitted tag whose value is a physical conclusion must be gated by a check that FAILS when that
> value is corrupted.**

For each such tag:
1. Add a check that reads **the same expression that is emitted** — ⛔ not a parallel re-derivation of it,
   and ⛔ not the CAS's own algebra. ⚠ A check that re-derives a value independently and then compares
   against a *separately written* literal certifies the literal, not the emitted value.
2. Register it in the script's `allChecks` list so it can reach `VERDICT`.
3. **Ablate it**: corrupt the emitted value, confirm `VERDICT` becomes `INTERNAL_CONTRADICTION`, restore.
4. ⭐ **Print an ablation-coverage line naming every tag now gated**, so coverage is machine-visible rather
   than a claim in prose.

⛔⛔ **DO NOT adjust a check so that an existing value passes.** ⭐ **If gating a tag causes its value to
change, that is a RESULT — report it prominently and keep the corrected value.** ⚠ Rigging a check to
preserve a prior output is the single worst outcome of this pass.

### Tags currently ungated — the work list

`S11BA_PROJECTION_FINITE` · `S11BA_PROJECTION_INFINITE` · `S11BA_PARITY_EVEN_JW` ·
`S11BA_PARITY_ODD_JW` · `S11BA_DYNAMIC_WINDOW_EXTRA_TERMS` · `S11BA_Z_BY_REGIME` ·
`S11BA_Z_BY_PARITY` · `S11BA_GRAZING_BEHAVIOUR` · `S11BA_RELATIVE_FLUX_CHANNELS` ·
`S11BA_CONTROL_WINDOW_PARITY` · `S11BA_CONTROL_INTERVAL_SYMMETRY` ·
`S11BA_DISCARDED_CONVECTIVE_ORDER` · `S11BA_ADDED_MASS` (its sign sub-field) ·
`S11BA_DEGENERATE_LOCI` (its undriven sub-fields)

⚠ **Specific patterns to eliminate wherever they appear**, all found by ablation:
- a status gated on a predicate containing **no model quantity** (e.g. `Reduce[0 == 0, …]`);
- a conclusion written as a **literal** beside a computation that does not produce it;
- a check that is a **`Solve` round-trip** or other invertible relabelling, which is invariant under
  swapping the very things it appears to distinguish;
- a dimensional check of the form `(X − 2Y) + 2Y == X`, where the symbol under test is *defined* by the
  equation being asserted.

## ⭐⭐ THE ONE GENUINE GAP — `S11BA_Z_BY_PARITY`

⛔ **This one is not a weak check; the problem was never posed.** The reported per-face response is
obtained from a ratio that returns its own input for **any** face motion whatsoever — including motions
that are not thickness or centre-shift modes at all. ⇒ the tag is an identity, and it certifies nothing.

⭐ **Repair by actually solving it:**
1. Write the bulk potential in **each** half-space with its own unknown amplitude.
2. Impose the face condition at **each** face.
3. Solve the resulting system for the amplitudes.
4. Derive the per-face response **from that solution**, then evaluate it on the thickness and centre-shift
   combinations.
5. ⭐ **Report whether the per-face response depends on the combination**, as a computed statement.
6. ⭐ **Ablate with a motion that is neither combination** and confirm the reported quantities change.

⚠ The specification warns that boundedness alone leaves four amplitudes against two face conditions.
⛔ That warning must be **answered** by the solve above, not sidestepped. State explicitly how many
amplitudes remain after the radiation condition and how many conditions fix them.

## ⭐ TWO TAGS THAT CLAIM MORE THAN THE COMPUTATION EARNS

**`S11BA_DISCARDED_CONVECTIVE_ORDER`** reports a **relative** order in `v₀/c_s0`, but the computation
performed yields only the **absolute** power of `v₀`. ⇒ either perform the comparison that makes it
relative — which requires relating a spatial derivative to a time derivative through the sound speed — or
⭐ emit `NOT_ESTABLISHED` and say which step is missing. ⛔ Do not restate the current claim.

**`S11BA_DIM_ROUTE_KIND_*`** labels are drawn from a hand-maintained list rather than determined by the
route. ⇒ **derive the label**: a route whose asserted equation *defines* the symbol under test is
`DEFINITIONAL`; only a route that reaches the symbol through an independent physical relation is
`INDEPENDENT`. ⚠ Two labels currently read `INDEPENDENT` for routes structurally identical to one labelled
`DEFINITIONAL`. ⛔ Do not resolve that by relabelling the third one; determine all of them by rule.

## ⛔ WHAT NOT TO DO

- ⛔ Do not add a checking framework, a harness, or checks that verify other checks. One check per gated
  tag, in the existing style.
- ⛔ Do not change the script's structure beyond what the above requires.
- ⛔ Do not weaken or delete any check that currently fires. **Eight do.**
- ⛔ Do not add new tags beyond the ablation-coverage line.

## Script conventions

Unchanged: standalone, print-only, no arguments, no exports, no external file reads; `WL_` prefix at emit;
strip `ConditionalExpression[0, …]`; test poles with `1/expr == 0`; total runtime under **10 minutes**.
