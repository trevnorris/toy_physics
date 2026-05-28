# Checkpoint Trust Audit

This document classifies the `25` checkpoint stages in
`CITATION_SUPPORT_SET.md` by current audit strength.

It does not replace:

- `STAGE_PROVENANCE_INDEX.md`, which keeps exact file paths,
- `STAGE_VERIFICATION_COVERAGE.md`, which tracks repo-wide coverage classes,
- `CHECKPOINT_CONSTANT_PROVENANCE.md`, which tracks the no-magic-numbers log.

Snapshot date: `2026-05-28` (batch V.1 close — first-pass paper-grounded audit on stages 164-175. **No checkpoints in V.1 range.** The 25 checkpoint trust tiers are unaffected by V.1 (no checkpoint stage in 164-175); cumulative checkpoint trust unchanged. Previous IV.6 entry retained below.)

Snapshot date prior: `2026-05-28` (batch IV.6 close — first-pass paper-grounded
audit on stages 151-163. **No checkpoints in IV.6 range.** Cumulative checkpoint
status (105 from IV.2, 096 from IV.1, 089/090 from III.5, 069 from III.3, 051
from III.2) unchanged. IV.6 closed range 001-163 paper-aligned at v2 depth.
Renumbering pattern (banners −17 across 9 stages, notes H1 −85 for 159/160,
body citations −51/−85/−102), plus stage 158 paper-card Checks downgrade —
**first forward (downstream) carry-forward in v2 batches** (items 2-3 cite
`\ref{stage:159}` / `\ref{stage:162}` / `\ref{stage:163}`). Twelfth consecutive
zero-redirection batch. One material_change (stage 151 SymPy rewrite to mpmath
after `sp.integrate` with free `Pi_star` hung; downstream impact zero — script
verifies algebraic identities not numeric carry-forwards). Previous IV.5 entry
retained below.

Snapshot date prior: `2026-05-27` (batch IV.5 close — first-pass paper-grounded
audit on stages 139-150. **No checkpoints in IV.5 range.** Cumulative checkpoint
status (105 from IV.2, 096 from IV.1, 089/090 from III.5, 069 from III.3, 051
from III.2) unchanged. IV.5 closed range 001-150 paper-aligned at v2 depth.
Same renumbering pattern as IV.4 (banners −17, notes H1 +102, body citations
−51/−102), plus stage 144 paper-card Checks downgrade mirroring stage 134.
Eleventh consecutive zero-redirection batch. Previous IV.4 entry retained below.

Snapshot date prior: `2026-05-27` (batch IV.4 close — first-pass paper-grounded
audit on stages 127-138. **No checkpoints in IV.4 range.** Cumulative checkpoint
status (105 from IV.2, 096 from IV.1, 089/090 from III.5, 069 from III.3, 051
from III.2) unchanged. IV.4 closed range 001-138 paper-aligned at v2 depth.
Previous IV.3 entry retained below.

batch IV.3 close — first-pass paper-grounded
audit on stages 115-126. **No checkpoints in IV.3 range.** Cumulative checkpoint
status (105 from IV.2, 096 from IV.1, 089/090 from III.5, 069 from III.3, 051
from III.2) unchanged. IV.3 closed range 001-126 paper-aligned at v2 depth.
Previous IV.2 entry retained below.

batch IV.2 close — first-pass paper-grounded
audit on stages 103-114. **Checkpoint stage `105` (chi_Q fix from outgoing
DtN) returned `verified` after first-pass v2 audit + apply + verifier
cycle at the higher-bar standard.** 105 cleared 2 findings: F1
(mathematica_transliteration) — full `.wl` re-author with a structurally
distinct independent path (unfactored ratio + `Apart` round-trip
verification of the canonical 3/4 + 1/(4D) partial-fraction form +
`SeriesCoefficient` operator-form coefficient extraction + `Reduce`-over-
reals for chi_Q + polynomial inversion of `Λ_def·Y_def = -3` for the
deformed branch). All variable names changed (no `yRet`/`lamDef`/`yDef`/
`omegaQ`/`sigmaCan`); the forbidden substring `(1/4)/(1 - omega^2/omegaQ^2`
is absent from the rewritten script. F2 (paper_misalignment, Cluster C)
— banner/docstring labels updated from stale "Stage 88/STAGE 088" to
"Stage 105/STAGE 105" matching `\label{stage:105}`. Both engines exit 0
with 10 PASS lines on Mathematica (sigma_Q^can identity + Apart-form
check + 3 series-coefficient checks + chi_Q-1 + 4 deformed-branch
coefficient checks) and 13 explicit `expect_zero` results on SymPy. The
substantive `chi_Q = 1` derivation passes via two independent paths
(SymPy's `sp.solve` on coeff(omega^5)/I = a^5/(27 c_s^5) vs Mathematica's
`Reduce[..., chiQ, Reals]`); the deformed-branch closure also passes via
two genuinely independent paths (SymPy: series of `-3/Lambda_def`;
Mathematica: polynomial inversion via `Solve[Lambda_def · y_ansatz = -3]`
matching coefficients). This is the third paper-grounded v2 checkpoint
verified at the higher bar (after 069 III.3 v2, 089/090 III.5, 096 IV.1).

Snapshot date prior: `2026-05-27` (batch IV.1 close — first-pass paper-grounded
audit on stages 091-102. **Checkpoint stage `096` (geometry-lane check verdict)
returned `verified` after first-pass v2 audit + apply + verifier cycle.**
096 cleared 3 findings: F1 (tautological_check) — deleted the trivially-true
`eps_2`/`eps_4`/`zeta_req - c_pole/c_geom` asserts (these were forced-zero by
the line-above definitions), replaced with `Print`/`Print` of values so the
transcript still shows the carried zeros; F2 (insufficient_verification) —
appended a printed-only "HYPOTHESIS CARRIED" block annotating the Part III
minimal-isotropic-module + grouped real `P_2` carrier conditioning per paper
card Check (iii); F3 (banner-label correction) — STAGE 079 → STAGE 096 in
both engines, plus a docstring "Stage 75" → "Stage 092" fix for the
upstream-obstruction reference. The 5 substantive final-ledger asserts
(`c_pole - 1/4`, `c_geom - 3/4`, `Yhat_Q^cons - […]`, `rho_alpha - 4/3`,
`zeta_req - 1/3`) remain in both engines as documented arithmetic
confirmations of paper Check (i). The 15-mode orthogonality block (5 Y2A
labels × 3 checks each) at the top of the script is the substantive
geometric-decoupling content that justifies the verdict; PASS-line count of
20 (15 orthogonality + 5 final ledger) confirms full coverage.

Snapshot date prior: `2026-05-27` (batch III.5 close — first-pass paper-grounded
audit on stages 085-090. Checkpoint stages `089` (Family-1 minimal isotropic
verdict) and `090` (updated reduced status) fall in III.5's range and **both
returned `verified` after first-pass v2 audit + apply + verifier cycle.**
089 cleared 4 findings including the paper_misalignment chain-closure gap
(Omega(Pe→0)=1 limit + zeta_F1(0)=A_F1 link + explicit Pe_req=0 assertion,
locking the paper-side boxed Output that scripts previously did not verify);
090 cleared 3 findings including a Mathematica-side definitional tautology
(rho_alpha and zeta_req now derived from c_contact=3/4, c_pole=1/4 rather
than hardcoded) plus a both-engines Pe_req=0 carry-forward proxy.

Checkpoint stage `069` falls in III.3's range (061-072)
and returned **clean (0 findings)** under v2 — the three-zone classification
combining Stage 066's matched window with Stage 068's resonance penalty
verified at v2 depth via parameterized `W_match` generator + monotonicity
check (SymPy) and `Cres2Prim` primitive + `Pres = 1/Cres2` derivation +
`PresGap` via `Solve` (Mathematica). Engine independence genuine: SymPy
uses `Pres_gap` primitive with multiplication-by-Pres; Mathematica uses
`Cres2Prim` primitive with division-by-Cres2 and Solve-based extraction.
The batch closed with all 12 stages verified, **4 paper_misalignment
(2 substantive + 2 banner relabels) + 9 script-side findings**;
**zero `material_change: true` flags**. Stage 068 (v1 material_change:true)
returned clean at v2 — Solve-derived `Wfail_res`/`Wfail_match` preserved
on resonance-corrected premises (Mathematica upgraded from `Solve` to
`Reduce` without value change). **One orchestrator hot-fix on stage 064
Mathematica**: `Integrate[]` does not factor constant multipliers when
the integrand has unspecified symbolic functions — verify integrand
equality first (NEW pitfall #9 candidate). With III.3 v2 closed, entire
range 001-072 is now paper-aligned at v2 depth.

Previous (III.2 v2) snapshot: Checkpoint stage `051` falls in III.2's range
(049-060) and returned **clean (0 findings)** under v2 — load-bearing boxed
criterion `Pi_tr <= 2 C_mix = 16 Lambda(1-eps)/pi^2` plus all six
notes-side deliverables (closed `Pi_tr`, endpoint limits, `zeta_req`
threshold, `Lambda_twin,req`, `M_mix^(twin,req)=G_tr/2`,
`Z_W^(twin,req)`, `xi_(2x)` quadratic root) verified at v2 depth across
8 SymPy and 9 Mathematica assertions. The batch closed with all 12
stages verified, **2 paper_misalignment + 14 script-side findings**;
**zero `material_change: true` flags**. Stage 050's paper card was
extended with a fifth boxed equation `S_n^{twin}(x;ε) < S_n^{max}(ε)
:= 1 + (1-ε)/((2n+1)^2 - ε)` (label `eq:app-stage050-Sn-max`) per Q1=(a)
— scripts unchanged. Stage 057 gained a local Pe-monotonicity numerical
sweep per Q2=(a) anchoring `partial_Pe zeta > 0` after destination
check confirmed Stage 056's scripts only verify the covariance identity
`dOmega_Pe/dPe = Cov(chi_0,s)/I_W`, not the sign — making (b) carry-forward
unsound. **One orchestrator hot-fix on stage 058**: Codex iter2's full
`sp.dsolve` / `DSolve` symbolic BVP solve for F2 hung sympy 7+ hours at
100% CPU; orchestrator replaced both engines' dsolve blocks with the
equivalent kernel-integral identity (sympy uses a numerical
`Delta = integral(K * Sigma_Pe)` sweep over 4 `(α,η,Pe)` tuples;
Mathematica relies on its pre-existing `delta independent integral
matches combination form` check at L84). **NEW pitfall #8 candidate**:
heavy BVP `dsolve` checks are not worth the symbolic cost; Green-function
identity verifies the same content in seconds. Stage 060 (v1
`material_change: true`) returned clean — v1 gain definition is sound
at v2 depth. With III.2 v2 closed, the entire range 001-060 is now
paper-aligned at v2 depth. Previous batch III.1 v2 text now superseded
by this entry; the I.1/I.2/II.1 v2 close text remains accurate for
those ranges.

Previous (II.1 v2) snapshot: Checkpoints `024` and `036` fall in II.1's range.
Stage 024 v2 had F1 mathematica_transliteration (Sections III/V port —
remediated to a 2x2 matrix-inverse-based independent derivation with
`-R` off-diagonal derived from `paper/parts/part01_parent_geometry.tex:956`
`+R_l A_l W_l` Lagrangian), F2 insufficient_verification (SymPy
`g.T * Mpair.inv() * g` anchor added with same `-R` sign), F3
tautological_check (C_alpha self-equality + equal-lane substitutions
replaced with `x20/x21/x22 reassembled` arbitrary-lane checks), and F4
insufficient_verification (explicit O(3)-collapse `D_{20,n}=D_{21,n}=D_{22,n}`
witness plus lane-breaking witness with non-zero linear-in-delta
coefficient on both engines). Section IV symbol-context reset +
`i4`/`i6` sphere-integral memoization fixed a separate performance hang
(>18min CPU → 25.05s total runtime). Stage 036 v2 had F1 tautological_check
(M_mix admissible/inadmissible witnesses replaced with parameter-derived
`Mmix_expr` evaluations) and F2 tautological_check labelling
(definitional self-consistency comment blocks added). Both checkpoints
verify with `material_change: false`. Three other II.1 stages had
paper_misalignment items (029 F1 docstring relabel, 029 F4 `alpha_crit`
trim with destination stage 031, 035 F1 paper polynomial coefficient
fix `206→189`, `138→121`) — none required user redirection (Codex
first-pass recommendations all held up; orchestrator cross-verified
Q2's destination claim). Zero `material_change: true` flags this batch
overall — all script-side edits were additions or removals of
tautological self-checks; no downstream-visible closed forms changed.
With II.1 v2 closed, the entire range 001-036 is now paper-aligned at
v2 depth. I.1 v2 close text remains accurate for batches III.1-III.4
and for the I.1 stage range.)

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
| `069` | Exact reduced threshold-window and profile-penalty algebra; no open review findings. Red-team batch III.3 (2026-05-22) reverified end-to-end with no tier shift after replacing the `Cres2`/`Wfail_res`/`delta_fail` definitional identities with a parameterized `W_match` generator + monotonicity check (SymPy) and a `Cres2Prim` primitive + `Pres = 1/Cres2` derivation + `PresGap` via `Solve` (Mathematica); upstream `Pres`/`Wfail_match` carry-forward now carries an explicit provenance comment. | None urgent. |
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
