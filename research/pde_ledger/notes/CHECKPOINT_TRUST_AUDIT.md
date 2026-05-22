# Checkpoint Trust Audit

This document classifies the `25` checkpoint stages in
`CITATION_SUPPORT_SET.md` by current audit strength.

It does not replace:

- `STAGE_PROVENANCE_INDEX.md`, which keeps exact file paths,
- `STAGE_VERIFICATION_COVERAGE.md`, which tracks repo-wide coverage classes,
- `CHECKPOINT_CONSTANT_PROVENANCE.md`, which tracks the no-magic-numbers log.

Snapshot date: `2026-05-22` (batch III.1 close — no checkpoint stages in
range 037-048, snapshot tier unchanged)

## Scope

- This is a stage-level trust baseline for the checkpoint support set.
- It uses the canonical stage cards, the existing review notes, the
  constant-provenance log, and representative spot-checks of the matching
  Mathematica mirrors.
- It is not a claim that every Mathematica mirror is fully independent.
  Independence is now governed separately by `MATHEMATICA_MIRROR_POLICY.md`.
- The checkpoint tier follows the weakest exposed issue in the current audit
  stack, not the strongest-looking subsection of a script.

## Assessment Rule

- `strong`
  exact symbolic re-derivation, genuinely nontrivial symbolic theorem check, or
  a deliberately narrow status-boundary checkpoint whose carried inputs are
  source-anchored and whose audit exactly matches the stated claim, with no
  open review findings that narrow that claim
- `moderate`
  useful audit, but narrower than the checkpoint claim, still dependent on a
  convention/assumption caveat, or framed as a status-consistency check whose
  carried inputs or summary scope are still too weak to stand alone
- `weak`
  tautological, self-substituting, or assumption-fragile enough that the
  checkpoint should not yet be trusted as citation support

## Disposition Rule

- `strong` -> `ready for citation support`
- `moderate` -> `internally useful but not citation-grade`
- `weak` -> `needs hardening`

## Snapshot

- `25` checkpoints are currently `strong`
- `0` checkpoints are currently `moderate`
- `0` checkpoints are currently `weak`
- Benchmark-only or probe-only numeric slices appear at `221`, `242`, `248`,
  and `253`. They are not currently trust defects because the theorem path is
  symbolic and the literals are explicitly labeled.
- Dedicated numerical-stress coverage now exists at `003`, `185`, `248`, and
  `253`.

## Current Trust Tiers

### Strong

Disposition: `ready for citation support`

| Stage(s) | Why currently strong | Immediate next action |
|---|---|---|
| `051` | Exact symbolic product law, endpoint limits, threshold rewrite, and closed root solve; no open review findings. | None urgent. Add numerical stress only if this threshold is later used numerically. |
| `069` | Exact reduced threshold-window and profile-penalty algebra; no open review findings. | None urgent. |
| `096` | Re-derives the isotropic `l=0 <-> l=2` decoupling and then evaluates the carried Stage 092 obstruction formula on the isotropic branch to recover the `3/4 + 1/4` conservative module. | None urgent. |
| `003` | Exact Schur-complement replay, exact one-mode pole split, grouped real `P_2` isotropy/anomaly bookkeeping, independent Mathematica mirror (red-team batch I.1 patched a multi-line `lRed = ...` continuation defect that had captured only kinetic terms -- downstream results unaffected, flowed through `mMat/kMat/cMat/oMat`; red-team batch I.2 caught the same continuation-defect class at non-checkpoint stage 021 and patched it in iter 2), and now-runnable shared numerical stress. | None urgent. |
| `022` | Exact grouped-`P_2` normalization bridge, explicit Stage-021 dictionary round-trip for `N0/N2/N4`, and independent Mathematica replay of the normalization-product solve; red-team batch I.2 (2026-05-21) reverified end-to-end with no material change after rewriting the Mathematica mirror Sections I/II/IV/V to `LinearSolve` + `SphericalHankelH1` (F1) and removing a tautological round-trip block (F2). | None urgent. |
| `023` | Exact weighted-projector calculus, representative one-port reconstruction of `Z_n/N_n`, full grouped-bundle assembly, and isotropic prefactor laws. Matching Mathematica execution coverage exists, but the trust grade rests on the symbolic theorem path rather than mirror independence. Red-team batch I.2 (2026-05-21) reverified end-to-end with no material change after replacing a tautological additivity check on `grouped_parts` (F1), swapping two tautological substitutions for closed-form comparisons (F2/F3), and adding a numerical-substitution route plus direct small-z Bessel expansion to the Mathematica mirror (F4). | None urgent. Add numerical stress later only if the grouped-bundle checkpoint starts carrying quantitative downstream claims. |
| `024` | Exact STF harmonic Gram/source-map closure, unequal-lane witness checks, exact `Y_20` triple-overlap matrix, and the grouped `(1,1/2,-1)` splitting law. Matching Mathematica execution coverage exists, but the trust grade rests on the symbolic theorem path rather than mirror independence. Red-team batch II.1 (2026-05-22) reverified end-to-end with no material change after rewriting the Mathematica mirror to use direct `Integrate[..., {theta, 0, Pi}, {phi, 0, 2 Pi}]` over the 2-sphere and removing the SymPy-named `pairings`/`i4`/`i6` helpers and `deltaPair/sPair/qPair/hPair/pPair` shorthands, replacing pre-substituted lane forms with an `xLane[lam_]` parameterizer (F1 transliteration). | None urgent. Add numerical stress only if future papers start using quantitative overlap-sensitivity claims rather than the symbolic angular law. |
| `036` | Exact support-feasibility function, exact loading split, onset boundary, and near-onset / endpoint structure under the explicit `0 <= xi < 1` boundary. Matching Mathematica execution coverage exists, but the trust grade rests on the symbolic theorem path rather than mirror independence. Red-team batch II.1 (2026-05-22) reverified end-to-end with no material change after adding a substantive symbolic kappa-based `F`-`R_target` identity, deleting tautological `R_target - X = X - X` and `gMaxTarget - alphaCrit` checks, re-deriving `dGTarget`/`gMaxTarget`/`gSeriesTarget` natively via polynomial form / `Limit` form / coefficient extraction, and adding a `disc + 72*delta^2 == 0` discriminant check (6 findings closed: tautologies, transliteration, insufficient_verification). | None urgent. Add numerical stress only if future papers start using this frontier quantitatively rather than as a symbolic admissibility gate. |
| `089` | Exact minimal-isotropic demand versus explicit Family-1 window/ceiling comparisons in both CAS layers; a closed arithmetic theorem once the upstream minimal module is accepted. | None urgent. |
| `090` | Narrow but explicit theorem-status boundary: the carried minimal module and carried Family-1 thresholds are all source-anchored, and both CAS layers replay the exact status verdict claimed in the note. | None urgent. |
| `163` | Exact off-family normal-coordinate transport packet in both CAS layers; the Family-1 numbers are explanatory readbacks only and do not enter the symbolic theorem checks. | None urgent. Add numerical stress later only if downstream work starts relying on the canonical-point coefficient readbacks quantitatively. |
| `105` | Exact outgoing-DtN fixing of `chi_Q`; symbolic theorem path is clean. | None urgent. |
| `112` | Exact hybrid / Robin branch solve with no open review issues. | None urgent. |
| `001--002` | The foundational harmonic bookkeeping, confinement sign, densitized-versus-weighted convention, monopole bridge, `4\pi` overlap factor, and conservative `(a,L)` / grouped-`P_2` reductions are now explicit in the human-facing stage material and checked in both CAS layers; red-team batch I.1 (2026-05-21) reverified end-to-end with no material change. | None urgent. Add numerical stress only if later downstream papers need quantitative sensitivity tests at this foundation layer. |
| `185` | Primitive microscopic ratios, tracking/nontracking/dressing monomial compilers, observable complement law, and zero-defect solve now all check symbolically in both CAS layers, with the earlier tautology concern removed. | None urgent. Keep the existing `185--187` numerical stress as secondary coverage. |
| `200` | Exact reference-free packet identities, cocycle law, mismatch compiler, and linearized four-scalar compiler in both CAS layers. | None urgent. |
| `203` | Exact graph-kernel theorem, inverse compiler, and repair-vector algebra in both CAS layers. | None urgent. |
| `218` | Exact combinatorial closure, exhaustive splice theorem, and sourced budget arithmetic; the theorem is symbolic / combinatorial rather than numerical. | None urgent. |
| `221` | Exact line-shape tradeoff laws and survival-window formulas; the numeric slice is explicitly probe-only. | Optional later numerical stress if this gate becomes a recurring quantitative citation. |
| `239` | Exact rigid-mouth physical chart, dependent-plane compiler, support-blindness, and orbit-lock theorem in both CAS layers. | None urgent. |
| `242` | Exact support-placement / orbit-lock compiler; rational sample is explicitly probe-only. | None urgent. |
| `243` | Exact leakage/work lane, non-rigid solve, recovery slice, and short-range limit firewall in both CAS layers. | None urgent. |
| `248` | Exact symbolic theorem path for energy conservation, threshold speeds, Coulomb reference formulas, and near-top action; benchmark numerics remain labeled, and dedicated dual-CAS stress now exercises multiple admissible corridors plus direct action quadrature. | None urgent. Widen only if future papers need denser dynamic sweeps or a concrete potential family. |
| `253` | Exact symbolic theorem path for lattice-turnover, calibration recovery, stiffness map, temperature ceiling, and screening ratios; benchmark numerics remain labeled, and dedicated dual-CAS stress now exercises micro, legacy-boundary, and targeted failure slices. | None urgent. Widen only if future work starts screening real candidate-host datasets. |

### Moderate

Disposition: `internally useful but not citation-grade`

No checkpoint in the current support set is classified `moderate`.

### Weak

Disposition: `needs hardening`

No checkpoint in the current support set is classified `weak`.

## Strongest And Weakest Current Checkpoints

### Strongest Current Shortlist

- `096`
  direct geometric decoupling plus derived conservative module
- `185`
  direct microscopic monomial compiler with repaired symbolic independence
- `200`
  exact reference-free packet and linearized compiler
- `203`
  exact graph-kernel / inverse / repair compiler
- `239`
  exact physical-chart and orbit-lock compiler
- `242`
  exact support-placement and coherent orbit-lock compiler
- `243`
  exact relaxed-branch lift with concrete leakage/work and recovery firewall

### Weakest Current Frontier

- No symbolic checkpoint currently sits below `strong`.
- No checkpoint currently has an immediate numerical-hardening defect inside the
  support set.
- Optional next wideners are `200`, `221`, and `242` if those checkpoints begin
  carrying quantitative downstream claims rather than symbolic citation support.

## Immediate Hardening Queue

1. Decide whether to widen numerical-stress hardening beyond the current four
   stress-backed checkpoints: `003`, `185`, `248`, and `253`.
2. If quantitative downstream use grows, prioritize `200`, `221`, and `242`
   before widening to the rest of the strong support set.

## Operational Result

The checkpoint support set now has a conservative trust baseline:

- `strong` checkpoints can support downstream citation routing now,
- `moderate` checkpoints should be cited with caution or hardened first,
- there is currently no checkpoint blocked on a known weak symbolic audit.
