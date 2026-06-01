---
batch: V.3
range: 188-200
total_stages: 13
verified: 13
findings_count: 13
findings_resolved: 13
findings_blocked_legitimate: 0
material_change_count: 0
clean_stages: [190, 192, 194, 196, 197, 198, 199]
status_only: []
dirty_stages: [188, 189, 191, 193, 195, 200]
checkpoints: [200]
consult: redteam/codex_reviews/_consult_V3.md
audit_date: 2026-06-01
verify_date: 2026-06-01
status: closed
---

# Red-team batch V.3 — Branch observables, isotropic target, home stretch

## Summary

13-stage audit unit for V.3 (`Part V.3 — Branch observables, isotropic target,
home stretch`), the **second batch of the resumed first pass**. **200 is a checkpoint**
(higher verify bar); no status-only units. All 13 reached `verified`; **7 clean**
(190, 192, 194, 196, 197, 198, 199), **6 dirty** (188, 189, 191, 193, 195, 200).
13 findings, **all resolved, 0 blocked, 0 stop-cold**. **Material change: 0** —
every change is a strengthening/route change (a new second engine, a
de-tautologized check, a demoted-to-definition identity, or a banner relabel);
no derived value, constant, identity target, or paper number moved, so no
`upstream_stale` propagation.

### Headline — the dual-engine rule correction (user-clarified 2026-06-01)

The defining event of this batch is a **standing-drift correction to the
dual-engine rule**, not a single finding. Before V.3, stages 188–199 were
**SymPy-only** (no `.wl`). During V.3 the user clarified that a Mathematica `.wl`
is **REQUIRED on every stage where Mathematica CAN independently verify the
result** — the test for omitting a `.wl` is "**is it possible**," NOT "is it
necessary." Single-engine acceptance had accreted on the late SymPy-only frontier
via orchestrator/Claude+Codex discretion without user sign-off; the auditor
prompt's line 118 had in fact required it all along.

Consequence for V.3:

- **12 NEW independent-route Mathematica `.wl` scripts were created (188–199)** by
  Codex — every one verified as a **GENUINE INDEPENDENT route** (0 transliterations
  accepted), and
- the checkpoint **200's pre-existing `.wl`** (which the audit caught as a
  transliteration) was **de-transliterated**.

So all 13 V.3 stages are now genuine dual-engine. The labor split was strictly
enforced: **Claude reviews** (audit + verify); **Codex writes ALL script code**
(designs and writes the new `.wl`); the directives stated only the requirement +
acceptance criteria, never script code.

## Per-stage findings tally

| Stage | Status | Findings | Notes |
|-------|--------|----------|-------|
| 188 | dirty | 2 + `.wl` | F1 vacuous "zero-set equivalence" (`M*0==0`) → generic-`obs` bijection round-trips + nonzero-determinant checks; F2 tautological `Cstar*Ctr−1` (`Cstar:=1/Ctr`) → `Cstar` from the independent appendix closed form, cross-checked `Cstar − 1/Ctr`. New `.wl`: compilers as Jacobians via `Solve` of the log-drift relations |
| 189 | dirty | 2 + `.wl` | **F1** tautological self-inverted selected-branch identity (see iteration-2); **F2** tautological transfer identities re-asserting their own definitions → derived from INDEPENDENT observable defect symbols (`Theta1`, `Xi1`, `Sigma_eta`). Banner relabel "STAGE 172"→"STAGE 189". New `.wl`. **One iteration-2 rework** (Section II) |
| 190 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: logarithmic Euler operator for the slippage laws |
| 191 | dirty | 1 + `.wl` | Banner relabel "STAGE 174"→"STAGE 191" (string-only, no math). New `.wl`: `Coefficient` + `D[]` Taylor cross-check, metric-built projectors |
| 192 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: constrained `LinearSolve` on the augmented 8×8 + `Solve` equivalence |
| 193 | dirty | 1 + `.wl` | F1 insufficient §IV "scalar/geometry firewall" — the hand-written quadratic `Deff` → an explicit block operator with LINEAR `chi` coupling, taking the genuine Schur complement (proves linear→quadratic). New `.wl`: `ArrayFlatten` block Schur via `Inverse`, `Solve[deltaPole==0]` |
| 194 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: native `SphericalHankelH1`, `Solve` for the pole/σ from coefficient matches |
| 195 | dirty | 2 + `.wl` | F1 X−X self-echo of the headline `m̂₀²χ_Q N_Q=1` → derived from the observable odd condition `(m̂₀²Γ̄₅−Γ̄₅^target).subs(P0, N_Q*P0_target)`; F2 deleted two definitional-echo checks. New `.wl`: derives `χ_Q` from the Hankel operator |
| 196 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: the tail wired in as a free symbol `l7`, `D[coeff,l7]=0` through z⁵ |
| 197 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: `χ_Q` two ways + `Reduce`/`Equivalent` |
| 198 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: `M_*` reconstructed as a Jacobian `Table[D[...]]`, `LinearSolve` orbit |
| 199 | **clean** | 0 + `.wl` | clean-verify confirmed. New `.wl`: `M_*` from additive exponent-weight vectors, `LinearSolve` projectors |
| 200 | dirty (ckpt) | 2, `.wl` de-translit | F1 transliteration (Section I hand-collapsed ratios → a `ratioSubs` helper-monomial-quotient route feeding the `Mderived` Jacobian); F2 HIGH tautological (Section III `Log[a^b]` collapse → full `ctrMonomial[...TActual...]/CtrTarget` exercising the exponents). The `.wl` was de-transliterated, not newly created |

**Totals:** 13 findings, 6 dirty, 7 clean, 0 status-only, 0 blocked. **12 new
independent-route `.wl` (188–199)** + 1 de-transliterated checkpoint `.wl` (200).

## Mathematica mirror policy — the dual-engine retrofit

The policy itself changed this batch: **dual-engine is required wherever an
independent Mathematica verification is POSSIBLE** (not merely where it is
"necessary"). **0 sanctioned mirrors were accepted in V.3.** All 12 new `.wl`
(188–199) are GENUINE INDEPENDENT routes (representative routes per stage are
listed in the tally table above), and the checkpoint 200's pre-existing `.wl` was
**de-transliterated** (Section I `ratioSubs` helper-monomial quotients feeding the
`Mderived` Jacobian; Section III `ctrMonomial[...]/CtrTarget` exercising the
exponents instead of a `Log[a^b]` collapse). All 12 of 188–199 are added to the
Independent-Mirror Set in `MATHEMATICA_MIRROR_POLICY.md`. Transliteration share
for V.3 was effectively the entire range as a coverage gap, all closed with
independent routes.

## SymPy-side de-tautologizations

Concentrated in the "passes-for-the-wrong-reason" findings:

- **188 F1/F2** — vacuous `M*0==0` zero-set equivalence → generic-`obs` bijection
  round-trips + nonzero-determinant checks; tautological `Cstar*Ctr−1` (with
  `Cstar:=1/Ctr`) → `Cstar` from an independent appendix closed form,
  cross-checked `Cstar − 1/Ctr`.
- **189 F1/F2** — tautological self-inverted selected-branch identity (resolved in
  iteration-2, below); tautological transfer identities re-asserting their own
  definitions → derived from INDEPENDENT observable defect symbols (`Theta1`,
  `Xi1`, `Sigma_eta`).
- **193 F1** — insufficient §IV "scalar/geometry firewall": the hand-written
  quadratic `Deff` was replaced by an explicit block operator with a LINEAR `chi`
  coupling, taking the genuine Schur complement (proves linear→quadratic).
- **195 F1/F2** — X−X self-echo of the headline `m̂₀²χ_Q N_Q=1` → derived from the
  observable odd condition `(m̂₀²Γ̄₅−Γ̄₅^target).subs(P0, N_Q*P0_target)`; deleted
  two definitional-echo checks.

## Stale stage-number label relabels (paper_misalignment, NOT escalated)

Resolved via the SETTLED canonical-stage-number convention = direction (a),
Codex-CONCUR — **not** escalated to the user. **189** banners "STAGE 172"→"STAGE
189"; **191** banners "STAGE 174"→"STAGE 191". String-only, no math.

The notes-PROSE legacy numbers are OUT of red-team scope (notes are paper files):
189's notes call it "Stage 240" / source 239; 191's notes call it "Stage 242" /
`..._stage242_...`. These are logged to `PAPER_CLEANUP_TRACKER.md` (P4-51) for a
later paper-side pass, alongside the non-blocking cosmetic banner-drift cluster on
the stages whose banners were NOT relabeled (188, 190, 194, 196, 197, 198, 199, 200).

## Iteration-2 rework (orchestrator-review catch) — 189 Section II

The iter-1 F1 fix left a **RESIDUAL tautology**: `Rtarget_oneport =
Lambda0*(1-epseta)/T2_direct` is a back-definition, so
`Rtarget_oneport*T2_direct − Lambda0*(1-epseta)` ≡ 0 for any input. The verify
agent flagged the "structurally identity by construction" smell but **passed** it;
the orchestrator **REJECTED** it on review (no-rubber-stamp, same class as the V.2
checkpoint-185 catch).

A read-only Claude+Codex consult (`redteam/codex_reviews/_consult_V3.md`,
codex_session `019e843e`, **CONCUR full**) confirmed:

- the selected-branch identity `R_target·T²=Λ₀(1−ε_η)` is **DEFINITIONAL** (rank-2
  compatibility) and must be **DEMOTED to a printed definition**;
- the genuine content is the **direct-slope bridge** `δln T_A² = ε·λ_A·Ξ₁`,
  perturbing the concrete continuum `T2_coh`, with `Ξ₁` supplied **INDEPENDENTLY**
  (not back-derived).

Iter-2 applied this in BOTH engines; re-verified — `Xi1_closed` confirmed
independent (built from the input perturbation amplitudes before the perturbation
path exists; touches neither `direct_slope` nor the T2 perturbation).

## Claude+Codex consult (`redteam/codex_reviews/_consult_V3.md`, session 019e843e)

ONE read-only consult, covering **189 Section II only**, CONCUR (full). Nothing
conceptual; nothing escalated to the user.

## Orchestrator catches

- **Post-verify catch (189 F1):** rejected the verifier's `verified` verdict on
  iter-1 because the re-pointed selected-branch identity was still a
  back-definition (≡ 0 for any input). Routed to the consult for the genuinely
  load-bearing anchor (the direct-slope bridge with an independent `Ξ₁`). This is
  the checkpoint higher bar + no-rubber-stamp rule working as designed (precedent:
  the V.2 checkpoint-185 F1 catch).

## Verification

All 13 verification files under `redteam/verifications/stage_188.md` …
`stage_200.md`. Final verdicts:
- `verified` (13): 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200.
- `needs_rework` → reworked → re-`verified`: 189 (1 loop, Section II).
- `blocked_unfixable` (0).

Material change: **0** (`material_change: false` on all 13 — new-second-engine /
de-tautologized / demoted-to-definition / banner-relabel only; no derived value,
constant, identity target, or paper number moved).

## Cumulative

Range 001-200 paper-aligned at v2 depth. **200/253 stages red-team verified
(79.1%)** (was 187 in the MANIFEST after V.2's range close; V.3 adds 13 across
188–200). **Second batch of the resumed first pass**; zero stop-cold, zero
material change, and the standing dual-engine drift on the 188–199 SymPy-only
frontier now closed (12 new independent `.wl` + the de-transliterated 200 `.wl`).

Next batch (sequential-audit-chunks rule, awaits explicit user authorization):
**V.4 onward = stages 201–253** (all currently `pending`). The planned full
end-to-end **second pass** remains a later cross-check, only after the first pass
reaches stage 253.
