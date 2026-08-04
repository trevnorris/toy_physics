# Independent physics review — a RETRACTION, and whether it is itself correct

## Artifact

Commit `3e41b463` in `/var/projects/toy_physics` (branch `ledger-v2-rebuild`). Read it with
`git show 3e41b463`. Files touched:

- `research/pde_ledger_v3/steps/S11bB_interface_assembly.md`
- `research/pde_ledger_v3/steps/S11b_HANDOFF.md`
- `research/pde_ledger_v3/paper/steps/S11b_interface_coupling_law.tex`
- `research/pde_ledger_v3/V3_STEP_PLAN.md`

## Background — what happened

Two blind engines assembled a linear brane–bulk interface law:

- `research/pde_ledger_v3/mathematica/S11bB_interface_assembly_mathematica_audit.wl`
- `research/pde_ledger_v3/scripts/S11bB_interface_assembly_sympy_audit.py`
- shared spec: `research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md` (§2b is the relevant part)

They computed a parameter region from non-negative interfacial power, plus a **conditional**
Onsager–Casimir relation. The step record then wrote this up as **"thermodynamically FORBIDDEN"**, killing
a result from the previous sub-step (`steps/S11bA_interface_response.md`). Commit `3e41b463` **retracts**
that, on the grounds that (a) the engines emitted a classification not a gate, (b) Onsager needs
microscopic time-reversibility the model never postulates, and (c) the reference state is driven by a
steady background flow `v₀`, so it is not an equilibrium.

## What to check — in priority order

**1. ⭐⭐ IS THE RETRACTION'S CENTRAL PHYSICS CLAIM TRUE? Attack this hardest.**
The retraction asserts the reference state is a **non-equilibrium steady state with a reservoir the
interfacial fluctuations can draw on**, because of the steady background normal flow `v₀`.
⛔ **A uniform background flow is removable by a Galilean boost, and about a boosted equilibrium both
passivity and Onsager reciprocity are restored.** So the claim survives only if `v₀` is a genuine
throughput that no boost removes — e.g. a steady mass flux *across* the interface, which is a preferred
surface. **Determine which it is, from the spec and the scripts, and say so.** If `v₀` as actually
implemented is a boost-removable uniform advection, ⛔ **the retraction's second reason collapses and this
is a blocking finding.**

**2. ⭐ Did the retraction over-correct?** Identify anything the computed region legitimately *does*
establish that the new text now under-claims or muddies. In particular: `Λ_A⁰ ≥ 0` — is that also merely
conditional, or does it follow from something weaker and safer than the second law?

**3. ⭐ Is "a non-passive coupling is admissible only with a NAMED reservoir and a STATED power budget"
operational, or a slogan?** State what a step would concretely have to produce to satisfy it. If it cannot
be discharged in practice, say so.

**4. Fidelity.** Every statement in the `.tex` card traceable to the step records; numbers, signs, symbols
and conditions checked term by term.

**5. ⭐ MACRO HAZARD.** Read `research/pde_ledger_v3/paper/macros.tex`. ⚠ Some stage-field-style macros are
**SUPPRESSED in the default build**. ⛔ Anything a reader must see that sits in a suppressed field is
invisible in the PDF and is a real finding. Verify against the built PDF, ⛔ not the `.tex` alone.

**6. Does a reader who stops halfway** come away believing either (a) the velocity leak is a live
prediction, or (b) that it was proven impossible? ⛔ Both are wrong; it is a conditional departure carrying
a debt.

## Required method

Derive independently before reading the retraction's reasoning. Ablate any load-bearing check and report
its literal output. ⛔ Do not accept the retraction's framing as the frame of your review.

## Physics filter

Report a finding only if it catches a way the **physics** could be wrong. ⛔ Do not report "the script would
be wrong on a different input."

## Do not read

`research/pde_ledger_v3/_scratch/`, any other reviewer's output, anything named `PREREGISTERED`/`PREREG`,
`research/pde_audit/`.

## Output — under 35 lines

1. ⭐ One-line verdict: **RETRACTION SOUND** or **RETRACTION DEFECTIVE**.
2. Per numbered check: pass/fail with the specific evidence.
3. ⭐ Your independent answer to check 1, with the evidence from spec or script that settles it.

⭐ A clean result is acceptable; ⛔ do not manufacture a finding.
