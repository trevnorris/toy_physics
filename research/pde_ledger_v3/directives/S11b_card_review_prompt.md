# FIDELITY REVIEW — the S11b ledger card

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S11b_interface_coupling_law.tex`

## Sources of truth — ⭐ the card must not exceed these

- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S11bA_interface_response.md`
- `/var/projects/toy_physics/research/pde_ledger_v3/steps/S11bB_interface_assembly.md`

⚠ **The second overturns a headline result of the first.** ⛔ Read both in full before judging the card.

## ⛔⛔ THE SPECIFIC RISK THIS CARD CARRIES

S11b-A produced an attractive, quotable result: *the velocity coupling dissipates iff `ωτ ≠ 0`* — "slow
disturbances see the leak, fast ones do not." ⭐⭐ **S11b-B established that thermodynamic admissibility
forbids that channel: a nonzero velocity coupling is admissible only at zero relaxation time, which is
exactly where it does no work.**

⇒ ⛔⛔ **THE FAILURE MODE IS A CARD THAT LEAVES THAT DEAD RESULT READING AS LIVE PHYSICS** — in the summary,
in a figure caption, in an abstract line, in a macro field, or simply by stating it early and retracting it
late. ⭐ **It belongs ONLY in departures: derived, then ruled out, with the reason.**

⭐ **Check specifically:** does a reader who stops halfway through this card come away believing the leak is
a live prediction of the model? ⛔ If yes, that is a blocking finding.

## What to check

**1. ⭐⭐ THE DEAD RESULT.** Is it anywhere except departures? Is its retirement stated with its reason?

**2. ⭐⭐ MACRO HAZARD.** Read `/var/projects/toy_physics/research/pde_ledger_v3/paper/macros.tex`. ⚠ Some
stage-field-style macros are **SUPPRESSED in the default build**. ⛔ **Anything a reader must see that sits
in a suppressed field is invisible in the PDF and is a real finding.** ⭐ Check specifically that the
**departures**, the **known limits**, and the **does-not-settle-confinement** caveat all actually render.
⚠ Verify against the built PDF or the rendered output, ⛔ not by reading the `.tex` alone.

**3. ⭐ OVERCLAIM AND FIDELITY.** Every statement traceable to the step records. Check numbers, signs,
symbols and conditions term by term — in particular the admissibility condition, the reciprocity relation,
the stability boundary `K₀`, the added inertia (`ρ_m/(2α)` per face on the thickness coordinate), the
affinity's normalization, and the count of independent quadratic invariants.

**4. ⭐ THE KNOWN LIMITS — are they present, visible, and unsoftened?** ⚠ The load-bearing one: **the
supplied objects cannot be falsified by the build**, because both engines transcribed them; their
verification comes from independent reviewer derivations. ⭐ Also required: the **slice** scope of the root
result, the two-port check's blind spots, and that the unequal-relaxation-time face response has no
independent standard.

**5. ⛔ Is the engines' `VERDICT: PASS` presented as a verdict on the physics?** It means only that each
script's internal checks did not contradict each other.

**6. ⛔ Does the card narrate process** — revisions, review rounds, defect counts, repair passes? Those
belong in the step records, ⛔ not in the ledger.

**7. ⭐ The transverse null.** Is it stated with the explicit caveat that it does **not** settle whether
light's confinement is unconditional, and that its coefficient's dimension is **undetermined** because the
coupling vanishes identically?

**8. ⭐ The growing root.** Is it stated as **not** an energy-conservation violation — the stored energy
having no minimum in that direction, with energy accounting closing exactly?

**9. Structure.** Is it written as **one** derivation of the interface law? ⛔ It must not be structured as
"part A then part B", must not narrate the sub-step split, and must not read as a result followed by its
retraction.

**10. Build.** Confirm the paper builds clean and report the page count.

## Do not read

`research/pde_ledger_v3/_scratch/`, any other reviewer's output, anything named `PREREGISTERED`/`PREREG`,
`research/pde_audit/`.

## Output — under 35 lines

1. ⭐ **A one-line verdict: SHIP or DO NOT SHIP.**
2. Per numbered check: pass/fail plus the specific evidence.
3. ⭐ Every fidelity discrepancy in full — the card's value and the step record's value.
4. ⭐ Anything a reader must see that does not render.

⭐ A clean result is acceptable; ⛔ do not manufacture a finding.
