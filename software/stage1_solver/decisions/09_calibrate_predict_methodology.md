# Decision 09 — Path-A test methodology: calibrate-and-predict (supersedes the blunt no-fine-tuning rule)

**Date:** 2026-06-17
**Driver:** USER methodology decision (the user is the authority on how the toy analog is validated), grounded
in Claude+Codex discussion. The user corrected an over-purist orchestrator framing. See
`[[feedback-calibrate-predict-methodology]]`. Methodology-flaw review:
`software/stage1_solver/_scratch/pathA_methodology_review_codex.{md,log}`.

## The decision

The Path-A `R_norm` test uses **calibrate-and-predict** (the Standard-Model methodology), NOT
blind-freeze-and-hope. This matches the project's own bootstrap from the start (`research/1pn_orbital_dynamics/`:
freeze `n=5` → `β=3` by subsystem consistency → `κ_PV` → …) and how the SM is validated (calibrate ~19 params
on a data subset → predict far more). Credibility = the **predictive surplus** = (independent held-out matches)
− (free parameters tuned). String theory is the cautionary tale (too many knobs → no surplus). The
discriminator is the surplus, not whether any parameter was fit.

**This SUPERSEDES decision-08's blunt "may NOT tune `S_Σ` to force `D0→0`" non-rescue rule** — that rule
wrongly forbade legitimate calibration. decision-08's PHYSICS stands unchanged (Fork B narrow: reciprocity ⇒ no
new numerator knob; Path A's only `R_norm` lever is the self-consistent background/`D0`). Only the *test design*
around it is reframed.

## The reframed test (now in `docs/pathA_preregistration.md`)

- **Calibration anchor (stated openly):** `R_norm = 0` (the GR quadrupole, the only externally-benchmarked
  Stage-1 constant). `S_Σ` parameters may be calibrated to it. `R_norm`/`chi_Q` are therefore the anchor, NOT
  held-out predictions.
- **Central question (reframed):** can a **narrow, physically-principled, predeclared** `S_Σ` family reach the
  GR near-pole `D0→0` while staying physically admissible (stable, `D0>0`, passivity)? Reaching it within a
  narrow family is a *mild* positive; the scientific content is the held-out surplus.
- **Held-out surplus = the test:** observables NOT used in calibration, judged by surplus-over-DOF.
- **TWO-FREEZE protocol:** (1) freeze model class/forms target-blind + the DOF list; (2) calibrate declared
  params only to `R_norm`; (3) freeze calibrated values + hash; (4) evaluate held-out, no residual refit.

## Codex methodology-flaw review — verdict + the accounting fixes (all folded into the prereg)

**Verdict: SOUND as adopted; NO fatal methodological flaw** (Codex did not relitigate the decision — it
accepted calibrate-and-predict as defensible for a toy analog). The accounting fixes it found, now in the
prereg:
- **Feasibility is degenerate under a BROAD `S_Σ`** (you can prescribe `U_{Σ,R}`/`U_{Σ,RR}` to hit any `D0`) →
  feasibility is informative ONLY relative to a narrow predeclared family; feasibility ≠ surplus (§L). *(This
  corrected an orchestrator overstatement that "reaching `D0→0` is a win" — it is only with a narrow family.)*
- **`S_Σ` is infinite-dim unless GATE-A pins finite forms** → GATE-A must declare finite forms + count ALL DOF
  (coefficients, knots, weights, branch switches, normalization, post-hoc family selection) (§E/§K/§M).
- **Held-out independence is limited** → after calibrating `D0`, `R_pole/P2/P4` share it (partially
  contaminated); `chi_Q` not independent at all. Only the non-`D0` bundle parts are genuinely held out (§H map).
- **Stage-1 surplus is thin** → report Stage 1 as "anchor calibration + internal-consistency held-outs"; real
  external surplus is downstream (WP3 response, Stage 2, EM side) (§H/§L).
- **Near-pole numerics can fake surplus** (`P2∝D0⁻²`, `P4∝D0⁻³`) → propagate full calibration covariance into
  every held-out residual; margin-to-Schur error bars mandatory (§J item 8).

## Status

`docs/pathA_preregistration.md` revised to this methodology (DRAFT for the conceptual GO; the narrow `S_Σ`
family + DOF count is the required pre-freeze GATE-A deliverable, post-GO). Conceptual GO (the `S_Σ` promotion)
remains the USER gate.
