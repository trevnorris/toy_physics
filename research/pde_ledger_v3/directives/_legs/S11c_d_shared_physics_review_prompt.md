# Independent physics review — S11c-d SHARED PHYSICS spec (v6)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md`

This is an **orchestrator-written physics spec** — the physics authority that TWO blind CAS engines (one SymPy, one
Wolfram that imports nothing and re-derives) will read to independently construct the same objects, plus their
comparator. An error here makes **both** engines agree on the same wrong thing, so it is reviewed until clear.

⚠ This is **spec v6** (a Codex-authored base, folded across several two-leg gates; the latest fold surgically
re-did the §1c Fourier-carrier treatment and §8: the c2 carrier convention is no longer supplied — each engine must
EMIT the 3-D→1-D reduction of its OWN closed-kernel Fourier symbols with both operands, and the comparator joins those;
`f̂_red` is the only supplied 1-D convention). Give it a **full, fresh independent review** — ⛔ do not assume the fold
is correct or complete, and ⛔ do not merely check a prior finding list; derive your own view from the sources and
report anything that changes what an engine computes or what the spec may claim. In particular, confirm the §1c
Fourier treatment is now both **executable** on the real consume-set (the c2 kernel `s11cc2ClosedCouplingKernel` — note
c2 strips its momentum-deltas, which live on c1's `dtn_kernel`) and **leak-free / blindness-safe** (no supplied `(2π)`
map, no pointer into any engine's construction script), and that §8 lists `f̂_red` SUPPLIED but the reduction factor
COMPUTED.

## Your role
Form your **own** view of what S11c-d must be from the requirements and the real c2 exports **first**, then read the
spec and report every way it could be **wrong** — physics wrong, an engine misled, an answer leaked, a cross-engine
closure over-claimed, a varying quantity frozen, or an order/label misstated. A finding is something that catches a way
the result could be wrong; ⛔ do not report "it would be wrong on a different input," and ⛔ do not manufacture a
finding to fill a quota — "nothing survives the filter" is a valid outcome if it is true.

## Read these FIRST (sources of truth), and form your own view before opening the spec
- `research/pde_ledger_v3/directives/S11c_decisions.md` — the governing requirements. In particular **N5** (:83),
  **N6** (:94), **N7** (:106), **N10** (:169), **N11** (:179), **N12** (:123), **N13** (:130), **N14** (:135),
  **N15** (:143), and the S11c-d row (:52).
- `research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md` — what c2 actually established vs owes (read
  "Established … vs owed", "Method notes", "Carry-forward", "What's next").
- `research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` — the immediate house-format precedent and the exact
  objects/status S11c-d consumes.
- `research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_disposition.md` — the N6 / Path-B disposition and the
  operand DEBT.
- The three-way profile-class consult that fixed §1c–§1d:
  `research/pde_ledger_v3/directives/_legs/S11c_d_profile_class_consult.md` (context only — assess the physics on its
  own merits, ⛔ not by whether it matches the consult).

## What to scrutinize (each is a QUESTION — derive/decide it yourself; ⛔ I have not told you the answer)
1. **Profile class / regime (§1c–§1d, N5/N12/N14).** Is naming a *localized interface* (asymptotically constant,
   `∫W₀′ ≠ 0`) a correct, non-global object, and is the *interface-vs-bump* distinction (`∫W₀′ ≠ 0` vs `= 0`) stated
   correctly? Are the three grades (`η` contrast, `σ_W` sharpness, `kL_W` kinematic) genuinely independent, and is the
   ⛔"weak-gradient is forbidden" framing right — i.e. would expanding in `σ_W` / `L_W→∞` actually damage the object?
   Is the ⛔"do not set `η→O(1)` in the first-shape-order operators" claim correct?
2. **Order bookkeeping (§3c, N12).** Verify the orders independently: converted amplitude `O(εη)`, absolute flux
   `O(ε²η²)`, incident flux `O(ε²)`, and the **fractional** conversion `C = J_conv/J_in = O(η²)` (the `ε²` cancels in a
   **linear** theory). Is emitting BOTH the `O(ε²η²)` and `O(η²)` labels correct, and does the spec avoid letting a
   mis-ordered term read as the excluded nonlinear program (N10/N12)? A short symbolic check is welcome.
3. **Strong-edge bridge (§3d, N7).** Is it honest that S11c-d gives only the weak matching coefficient `F′(0)` and the
   lab bounds `F(1)`? Verify the counterexample `C = sin²(ηG) → 0` at finite `ηG` and the claim that a nonzero Born
   coefficient gives **no** lower bound on strong-edge conversion. Is the order-unity edge correctly OUT of scope and
   named as a NEW construction (not a reduction)? Is anything about the strong edge over- or under-claimed?
4. **Honest supply of the c2 import (§1b, rule 6/16, M2).** Does the spec correctly mark: the c2 closed
   operator/kernel VALUES as **per-engine SOUND only**; the cross-engine content as the **N6 covariance thread only**
   (Reading B, the matched zeros are `(0)−(0)`, ⛔ not operand agreement); the operand DEBT (carrier 40 / source 76 /
   Φ 18) as UNADJUDICATED and **material to this consumer**; F and G as **WITHDRAWN**; and does it preserve BOTH
   `R_N6 = 18/288` and `R_cov` no-nonzero? Flag any place that quietly upgrades any of these to "closed," or that lets
   "representational-difference-unadjudicated" become "known to be just thickness."
5. **N6 control (§5a, N6/decisions:94-104).** Is the control the **independent shape/coordinate route + one-sided
   corruption**, and NOT the vacuous uniform limit? Are the tilt (N3) and advection (N4) probes correct, is
   `∇W₀→0` correctly **rejected** as a corruption, and is the `RHO4_CONSTANT` structural-absence handling
   (⛔ no `A−A`) right? Is corrupting an *anchoring* correctly excluded?
6. **Two photon-kill channels (§3b, N13).** Is "confinement = survival of the transverse channel" applied correctly,
   and are the continuum-conversion and **bound-mode capture** channels genuinely **distinct emitted objects**? Is the
   bound pole (a weak 1D well binds) real and correctly *not* a Bloch band? Is the confinement question emitted as a
   computed object, ⛔ not asserted?
7. **Answer/recipe discipline (M2/M3).** Does the spec **name the object** or does it manufacture a derivation-path
   question or leak an expected value/sign/order beyond the required `(ε,η,σ_W)` power-counting contract? Is the
   falsification numeric bound / `O(1)` reductio correctly WITHHELD (orchestrator-side)? Is any varying quantity frozen
   (M3/N14)?
8. **Completeness / internal consistency.** Any missing premise an engine needs; any internal contradiction; any
   obligation from the c2 house template (chain output / comparator / blind-WL / supplied-vs-computed) dropped or
   weakened; any `N11` carry-in silently lost.

## Method (this is a DOCUMENT, not a script)
Read the sources of truth first and form your own view; only then read the spec; quote **both** sides for every
finding (the spec's text and the source it contradicts or the derivation that refutes it). Where a claim is a
physics/order question (items 2, 3), a prose "I checked it" is weak — show the derivation or a short symbolic check and
its result. ⛔ There is no do-not-read list; ⛔ do not modify the working tree.

## Output
For each finding: **the spec location (quote it)**, **why it is wrong (with the source quote or your derivation)**, and
**the minimal fix**. End with an explicit verdict: **SOUND** (nothing outstanding changes what is computed or may be
claimed) or **NOT-SOUND** (list the findings). Separate must-fix (changes what an engine computes or the spec may
claim) from nits.
