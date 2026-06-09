---
batch: V.3
pass: 2
range: 188-200
total_stages: 13
verified: 13
findings_count: 1   # one load-bearing finding (200 transliteration); the rest are stale_output (refresh-only) + a paper-side card-text-lag cluster (user-deferred)
material_change_count: 0
clean_stages: [189, 191, 198]
status_only: []
dirty_stages: [200]   # the only stage with a source-code change (the .wl re-author); 188/190/192-199 had refresh-only stale_output + the deferred card-lag
checkpoints: [200]
value_recon_total: 171
value_recon_misaligned: 0
audit_date: 2026-06-09
verify_date: 2026-06-09
status: closed
---

# Red-team PASS-2 batch V.3 — Branch observables, isotropic target, home stretch

## Summary

Pass-2 re-audit of V.3 (`Part V.3`), stages 188–200 (13), isolated under
`redteam/pass2/`. **All 13 reached `verified`; `material_change: false` on all 13**
(the only source change is the checkpoint 200 `.wl` re-author — a method/route change;
no derived value, constant, identity target, or paper number moved). 0 stop-cold,
0 blocked, 0 Codex deviations; 200's codex-invoke exit 0 on iter-1 (no iter-2).

Value reconciliation (pass-2 augmentation): **171 deliverable values checked
batch-wide, 0 misaligned** (188=16, 189=21, 190=19, 191=11, 192=15, 193=7, 194=14,
195=15, 196=13, 197=6, 198=10, 199=16, 200=8). No genuine `paper_misalignment` →
**ZERO substantive paper/notes edits.**

## Headline — checkpoint 200 `.wl` re-author (first-pass de-transliteration was insufficient)

V.3 was the first pass's **dual-engine RETROFIT** batch (it created the 188–199 `.wl`
and "de-transliterated" 200's). The pass-2 heads-up flagged the standing risk that a
first-pass retrofit can be insufficient (V.1-173, IV.5-139, V.2's 8-of-12 cluster) —
and that risk was realized **on the checkpoint**.

The orchestrator ground-truth `.wl`-vs-`.py` read found 200's de-transliteration
**insufficient**: both engines still produced the load-bearing compiler matrix `M_*`
by the SAME autodiff-Jacobian of the SAME log-ratios (SymPy `q_pair.jacobian(Dvec)`
vs Mathematica `Table[D[qPair, Dvec]]`), §III posited the same hand-written orbit
closed forms, and §V used the same `Series`. The per-stage audit agent had called it
"borderline → accept" on a "each CAS runs its own simplifier" defense — the same
under-call framing the orchestrator overturned at 100/129/139/158. A dedicated
calibration agent (run against the 200 standard) confirmed the OTHER M_*-family
stages (188/192/198/199) are genuinely independent, isolating 200 as the lone port.

**Re-author-vs-accept is USER-LEVEL → the user authorized the re-author.** Per the
labor split (Claude reviews; Codex designs + writes all script code) the directive
stated requirement + acceptance only; Codex re-authored the `.wl`:
- **§I `M_*`** assembled from primitive monomial exponent-weight vectors
  (`chiCoreWeights`/`thermalWeights`/…); `qPair` survives only as the cross-check
  `qPair − M_*·Δx == 0` (no longer the derivation). [matches sibling 199's route]
- **§III orbit** SOLVED via `Coefficient` → `orbitCoeffMatrix` → `LinearSolve` of the
  log-linear residual system, then exponentiate (not posited). [matches sibling 198]
- **§V Packet-A** linearized by a base-point derivative `D[chiPerturbed,eps]/.eps->0`
  (not `Series`). [the IV.6-158 / V.1-173 distinction]

Verify-confirmed GENUINELY INDEPENDENT; **all checkpoint deliverables preserved** (the
carried Stage-192 `Mexpected`, §II witness-invariance zeros, §III mismatch chart
`q=((1+chi0_*)ln m_T, ln m_mu−ln m_K−F_*ln m_T, −ln m_K)`, §IV cocycle, §V Packet-A
coeff `eps(5 eps_beta+dSigma0/(3S)+9 dSigma5/S)`); **no checkpoint constant changed**;
the SymPy `.py` was UNTOUCHED. 200's Independent-Mirror-Set row updated to the pass-2
route (the first-pass `ratioSubs`+`Mderived` Jacobian row is superseded).

## The other 12 (188–199) — independence re-confirmed, 0 ports

Orchestrator ground-truth read + the calibration agent confirmed all 12 genuinely
independent. Discriminator vs 200: their `.py` POSITS/types the load-bearing objects
(literal matrices / closed forms / matrix-inverse) while the `.wl` DERIVES by a
different operation — 188 Series→`Coefficient` log-drift + `Solve` + `Outer[D]`;
192 constrained `LinearSolve` (augmented 8×8); 198 `Coefficient`+`LinearSolve` orbit
+ `D[]` Jacobian; 199 weight-vector `M_*` + native `Solve` Φ + `LinearSolve`
projectors; 189/190/191/193/194/195/196/197 logEuler operator / native
`SphericalHankelH1` vs explicit `j2+i·y2` / `Solve` / `D[]`-Taylor cross-routes.
All 13 V.3 stages are genuine dual-engine; **0 sanctioned mirrors.** Spot-confirmations
that held under fresh attack: 189's iter-2 demotion (`Ξ₁` independent), 195's X−X
de-taut, 193's linear→quadratic firewall, 196's non-vacuous `l7` free-symbol derivative.

## Stale outputs (INFRA, refresh-only, no Codex, not material)

11 committed SymPy outputs (188,190,192–200) carried pre-renumber `STAGE 171–183`
banners (the III.2/III.3 stale-committed-output class; the script source banners were
already canonical); 200's MMA output likewise (`STAGE 183`). The orchestrator
reliability-gate re-run (26 exec, all exit 0 — incl. the re-authored 200 running clean
independently) refreshed them to canonical `STAGE NNN`. 189/191 SymPy + the 188–199 MMA
outputs were already fresh (re-run byte-identical). Arbiter grep on the refreshed
outputs CLEAN (no stale self-epoch 171–183 banner, no `168π²`/`168%` class, no FAIL).

## Card-text lag (paper-side, USER-DEFERRED)

~9 V.3 cards' `\stagefield{Verification}` still say "Mathematica audit: none yet"
despite the retrofit `.wl` now existing and passing — a stale STATUS annotation, NOT a
value/identity mismatch. The audit agents applied an inconsistent threshold (some filed
it as `paper_misalignment`, some noted it); the orchestrator adjudicated batch-wide and
the **user decided to DEFER** the class to `PAPER_CLEANUP_TRACKER` (appended to the
existing first-pass P4-51 entry, kept `open`; the later paper pass is Codex-applied +
Claude-reviewed). Non-blocking; V.3 red-team (scripts) is complete.

## Numbering / hygiene

Cards CLEAN of the `+17 \stagefield{Purpose}` self-label class (205–217 absent). No
wrong-root path-typo introduced. Seat policy held: 200 = 1 `.wl`-touching Codex
session solo; orchestrator exec sequential after Codex done (no seat overlap, no
MANIFEST race). Pass-1 `MANIFEST.yaml` untouched (isolation held).

## Verification

All 13 verification files `redteam/pass2/verifications/stage_188.md` …
`stage_200.md`, final verdict `verified`. 6 trackers synced (PAPER_CLEANUP **P5-18**
= ZERO substantive paper/notes edits; MIRROR_POLICY 200 row updated; CHECKPOINT_TRUST
200 re-author logged; STAGE_VERIFICATION_COVERAGE pass-2 now I.1…V.3).

## Cumulative

Pass-2 progress: **I.1 … V.3 verified (200/253 at pass-2 depth)**. Next pass-2 batch
(sequential-audit-chunks rule, awaits explicit user authorization): **VI.1 = stages
201–218**.
