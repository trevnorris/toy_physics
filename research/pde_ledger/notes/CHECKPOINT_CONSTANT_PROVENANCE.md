# Checkpoint Constant Provenance Audit

This document records constant-provenance findings for the checkpoint stages in
`CITATION_SUPPORT_SET.md`.

The goal is narrow: make sure the checkpoint audits do not hide unexplained
literals behind apparently passing CAS scripts.

Snapshot date: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 3** — re-verification of IV.3-range stages 117/118/119/122 whose first-pass IV.3 fixes had been tainted-applied (committed in the IV.3 pass `b4e02d8`) while Codex was bypassed; all REMEDIATION-`verified` on 2026-05-29 with **`material_change: false` on every stage**. **NO checkpoints in the batch-3 range {117–122}** — IV.3 has no checkpoints, and the only nearby checkpoint (105) was batch 2 — so there is no checkpoint-constant provenance to log this batch. **Cumulative checkpoint-constant provenance is unchanged from the IV.2 close (105 retained at the higher-bar standard);** no derived or carried constant moved on any of 117/118/119/122 (e.g. 118's λ-sign resolution to MINUS aligns section V with the script's own already-correct section IV and the notes' boxed forms — no value changed). The 117/118 resolution directions came from a Claude+Codex consult (Codex base `019e74f7`, all CONCUR), recorded at `redteam/codex_reviews/_consult_batch3.md`. Previous batch-2 entry retained below.)

Snapshot date prior: `2026-05-29` (**IV.x/V.1 orchestrator-direct integrity remediation, batch 2** — re-verification of stages 105/106/109/112 whose first-pass IV.2 fixes had been applied orchestrator-direct while Codex was bypassed. **Checkpoint stage `105` is the only checkpoint in this batch.** 105 was REMEDIATION-`verified` (2026-05-29) with **`material_change: false`** — the `.wl` chi_Q derivation was rewritten along an independent residue/`Reduce`-witness route, but no derived or carried constant changed, so **105's constant provenance is unchanged from the IV.2 close** (`sigma_Q^can = 4 a^5/(27 c_s^5)` still pole-scale-derived in-script; `chi_Q = 1` still derived non-tautologically — now via the residue/`Reduce` witness on the Mathematica side; deformed-branch coefficients still imported with explicit provenance; zero unexplained literals). The IV.2 105 provenance entry below remains accurate. The other three batch-2 stages (106, 109, 112) are NOT checkpoints, so no further checkpoint-constant provenance to log; cumulative unchanged. Previous V.1 entry retained below.)

Snapshot date prior: `2026-05-28` (batch V.1 close — first-pass paper-grounded audit on stages 164-175. **No checkpoints in V.1 range**, so no checkpoint-constant provenance to log; cumulative unchanged. Previous IV.6 entry retained below.)

Snapshot date prior: `2026-05-28` (batch IV.6 close — first-pass paper-grounded
audit on stages 151-163. **No checkpoints in IV.6 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.5 entry retained below.

Snapshot date prior: `2026-05-27` (batch IV.5 close — first-pass paper-grounded
audit on stages 139-150. **No checkpoints in IV.5 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.4 entry retained below.

Snapshot date prior: `2026-05-27` (batch IV.4 close — first-pass paper-grounded
audit on stages 127-138. **No checkpoints in IV.4 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.3 entry retained below.

batch IV.3 close — first-pass paper-grounded
audit on stages 115-126. **No checkpoints in IV.3 range.** Cumulative
checkpoint-constant provenance unchanged from IV.2 close (105 retained at
higher-bar standard). Previous IV.2 entry below.

batch IV.2 close — first-pass paper-grounded
audit on stages 103-114. **Checkpoint stage `105` (chi_Q fix from outgoing
DtN) verified after first-pass cycle at the higher-bar standard.**
Constant provenance assessment for 105: the load-bearing constant
`sigma_Q^can = 4 a^5/(27 c_s^5)` is *derived* in-script from the pole-scale
identity `sigma_Q^can = (9/8)/Omega_Q^5` with `Omega_Q = 3 c_s/(2 a)`,
checked via `expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", ...)` (rather than
literal-asserted). The canonical odd-coefficient identification
`chi_Q = 1` is derived non-tautologically via two independent paths: SymPy
solves `sp.Eq(Yret.coeff(omega, 5)/I, a^5/(27 c_s^5))` directly for
`chi_Q`; Mathematica uses `Reduce[c5/I == a^5/(27 c_s^5), chiQ, Reals]` on
an `Apart`-decomposed retarded-module form. The deformed-branch
coefficients `(1, 1/9, 4/81, xi_Q/27)` for `Y_def` are derived in SymPy
via `sp.series(-3/Lambda_def)` and in Mathematica via a structurally
distinct polynomial inversion `Solve[Lambda_def · y_ansatz = -3]`, with
neither engine importing the other's RHS. **Zero unexplained literals in
checkpoint 105's scripts**; every constant is either pole-scale-derived,
imported with explicit provenance from stages 074/088 (now 104) on the
SymPy docstring, or pinned via a paper-quoted carry-in (`Lambda_out`
fingerprint coefficients from stage 104).

Snapshot date prior: `2026-05-27` (batch IV.1 close — first-pass paper-grounded
audit on stages 091-102. **Checkpoint stage `096` (geometry-lane check verdict)
verified after first-pass cycle.** Constant provenance assessment for 096:
the four cardinal constants (`c_pole = 1/4`, `c_geom = 3/4`, `rho_alpha = 4/3`,
`zeta_req = 1/3`) are *derived* in-script from the static-limit hypotheses
`eps_2 = eps_4 = 0` (themselves established by stage 094's orthogonality
output `K_{g,2} = K_{g,4} = 0`) via the `c_pole = (1 + eps_4)/(4(1 + eps_2)^2)`
obstruction formula carried from stage 092. The `Yhat_Q^cons(omega) = 3/4 +
(1/4)/(1 - omega^2/Omega_Q^2)` partial-fraction form is built from these
derived constants, not literal-asserted. The 15-mode orthogonality block at
the top derives `K_{g,2} = K_{g,4} = 0` from explicit angular integrals
(5 Y2A labels × 3 checks each: overlap, Laplace eigenvalue, gradient cross).
The `l = 2` Laplace eigenvalue `6 = ell(ell+1)` is documented in the
constant-provenance docstring. **Zero unexplained literals in checkpoint
096's scripts**; every load-bearing constant is either derived in-script
from orthogonality / obstruction-formula inputs or carry-forward with explicit
source anchor.

Snapshot date prior: `2026-05-27` (batch III.5 close — first-pass paper-grounded
audit on stages 085-090. Checkpoint stages `089` and `090` (both in III.5)
verified after first-pass cycle. **Constant provenance assessment for III.5
checkpoints**: 089's `Pe_suff_chi = 96.5285247264386` and `Pe_fail_chi =
11220.5441626259` are now anchored by an explicit provenance comment
naming `scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt` as the
upstream source (SymPy side; per pitfall #10 SymPy nsolve was not used).
On the Mathematica side, the same Pe values are *rederived* via
`FindRoot[zetaF1[pe] == zetaTarget, {pe, …}]` from notes-quoted
`rho_target - 1`, giving a second-engine independent path. 090's
`c_contact = 3/4` and `c_pole = 1/4` are paper-quoted minimal-isotropic
module coefficients; both engines derive `rho_alpha = 1/c_contact` and
`zeta_req = c_pole/c_contact` from them (no longer hardcoded on the
Mathematica side). Both checkpoints add a `Pe_req = 0` carry-forward
proxy with source-anchored comments. Every load-bearing constant in
III.5 checkpoint scripts is now either derived in-script or
carry-forward with explicit source anchor.

Checkpoint stage `069` falls in III.3's range; returned
**clean (0 findings)** under v2 with no hardcoded constants — every
constant in 069's scripts is either derived (`Cres2` from definition,
`PresGap` via `Solve` in Mathematica) or carried forward with source
anchor (`Cres2 ≈ 0.99441883...` and `Pres ≈ 1.005612487...` carried
from stage 067's sech-Gaussian benchmark with 60-dps mpmath provenance,
`Pe_req` carried from stages 048/049/052). The III.3 v2 sweep surfaced
0 new hardcoded_result findings across the batch. The 2 substantive
paper_misalignment items at stage 062 added: (a) closed-form definition
`C_sp_sq := Osp²/(Nss·Npp)` (no new constant) and (b) Cauchy
parameterization `Osp = cos(θ)·√(Nss·Npp)` (introduces only the
declarative symbol `θ`); the σφ coupling sign flip changed only the sign
of an existing expression. The two banner-relabel items (067, 072) are
pure string changes. III.2 v2 update text remains accurate for that range;
I.1, I.2, II.1, III.1 v2 close text remains accurate for those ranges;
III.4 update remains accurate for batch III.4.)

## Audit Rule

Every constant used in a checkpoint audit should fall into one of three buckets:

- `derived in audit`
- `carried forward with source anchor`
- `probe-only numeric value labeled`

Anything outside those buckets is treated as an audit defect.

## Completed Entries

### Stage 001

- Canonical stage: `paper/stages/stage_001.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `1/(2 sqrt(pi))`
  derived as the normalized real monopole harmonic `Y_00`
- `2 sqrt(pi)`
  derived from the mouth-average extraction rule
  `delta a = q_00 / (2 sqrt(pi))`
- `6`
  derived as the specialization `ell(ell+1)` at `ell = 2`

Audit note:

- The earlier surface-measure ambiguity is now explicit rather than hidden.
  The Stage 001 audits check both the weighted wall-action form and the
  densitized convention actually used by the ledger.
- The current canonical Stage 001 card now states that the ledger adopts the
  densitized one-dimensional convention, rather than leaving that choice
  implicit.
- The added Mathematica mirror rechecks the harmonic normalization, chain-rule
  sign, and weighted-vs-densitized wall-action split in a second CAS.

### Stage 002

- Canonical stage: `paper/stages/stage_002.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage002_breathing_reduction_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `1/(2 sqrt(pi))`
  carried from the normalized `Y_00` convention fixed in Stage 001
- `2 sqrt(pi)`
  carried from the Stage 001 mouth-average bridge
- `4 pi`
  derived in the Stage 002 audit as `(2 sqrt(pi))^2`
- `6`
  derived as the specialization `ell(ell+1)` at `ell = 2`

Audit note:

- The audits check the conservative `(a, L)` matrix reduction and the grouped
  real `P_2` degeneracy without introducing any free numerical coefficients.
- The current Stage 002 note and canonical stage card now state explicitly that
  the Stage 001 surface weight has already been absorbed into the effective
  axial coefficients and profiles before the reduced overlaps are written.
- The added Mathematica mirror rederives the `Y_00` bridge, the `4 pi`
  reduction factor, the conservative two-mode matrix system, and the `ell=2`
  restoring shift independently of the SymPy implementation.

### Stage 003

- Canonical stage: `paper/stages/stage_003.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage003_bdg_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage003_bdg_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage003_021_foundational_stress.py`
  and
  `mathematica/numerical/stage003_021_foundational_stress.wl`
- Verdict: `clean for current symbolic + stress scope`

Constants reviewed:

- none in the theorem path
  the Schur-complement kernels, low-frequency moments, pole formulas, and
  grouped trace / anomaly identities are all checked symbolically in both CAS
  layers
- JSON sample values in
  `scripts/numerical/stage003_021_foundational_samples.json`
  probe-only numeric values labeled

Audit note:

- The current Stage 003 checkpoint is not resting on hidden literals.
- The symbolic theorem path is exact in both CAS layers.
- The shared numerical-stress harness now resolves its repo-local sample JSON
  after the `research/pde_ledger/` move, so the stress layer is runnable rather
  than merely listed in the stage card.
- The sample JSON is used only for perturbative-validity and scaling probes; it
  does not supply any constant used to derive the Stage 003 formulas.
- Red-team batch I.1 (2026-05-21) identified and patched a Wolfram Language
  multi-line continuation defect in the original `.wl`: the `lRed = ...`
  assignment spanning several lines only captured the kinetic terms, missing
  the potential and coupling additions. Downstream results were unaffected
  because the dispersion derivations flowed through `mMat`, `kMat`, `cMat`,
  and `oMat` rather than from `lRed`. The corrected `.wl` now adds the missing
  terms via parenthesised `lRed = lRed + (...)` and verifies all four EL
  residuals against the SymPy convention.

### Stage 022

- Canonical stage: `paper/stages/stage_022.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `27`, `9`, `81`
  derived in the compact outgoing `l=2` fingerprint
  `a^5/(27 c_s^5)`, `a^2/(9 c_s^2)`, and `4 a^4/(81 c_s^4)`
- `54`, `6`, `8`, `15`, `5`
  derived by solving the invariant normalization product against
  `2 G / (5 c^5)` and then collapsing the `K_2`, `K_4` target formulas
- there are no probe-only decimals in the theorem path

Audit note:

- The current SymPy audit explicitly performs the Stage-021 dictionary
  back-substitution round-trip for `N0`, `N2`, and `N4`; the earlier review
  note claiming that gap was stale.
- The Mathematica mirror independently replays the same round-trip and the
  normalization-product solve.
- The checkpoint therefore no longer depends on unverified dictionary
  substitutions or unexplained literals.

### Stage 023

- Canonical stage: `paper/stages/stage_023.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `5`, `20`, `4`
  derived in the grouped `(1,2,2)` metric as the exact `Ggrp` norms of
  `ebar`, `ea`, and `eb`
- `1/5`, `2/5`, `1/10`, `1/2`
  derived from the corresponding rank-one projectors and grouped inverse-map
  formulas under the same weighted metric
- `2`, `3`
  derived in the prefactor-cancellation solve for the exact `N2` and `N4`
  constant-prefactor targets
- `54`, `27`, `5`
  carried forward with source anchor from the Stage-022 normalization bridge;
  Stage `023` now rebuilds `Gamma5_port` through the Stage-021 exact outgoing
  `l=2` branch before using the invariant product
- there are no probe-only decimals in the theorem path

Audit note:

- The earlier “assembled reduced coefficients only” caveat is now stale.
  Both CAS layers explicitly reconstruct representative one-port `Z_n` and
  `N_n` formulas from the underlying `(\Delta, S, Q, H, P)` lane data before
  assembling the grouped packet.
- The Stage-022 outgoing odd coefficient is no longer a bare literal at the
  bundle level. Both CAS layers now rebuild `Gamma5_port` through the same
  Stage-021 exact outgoing route used in Stage `022`.
- The remainder of the checkpoint then stays symbolic: grouped decomposition of
  `D_A0`, `D_A2`, `D_A4`, `N_A0`, isotropic reduction, prefactor constraints,
  anisotropy transport, and monotonicity derivatives.
- The current checkpoint therefore does not depend on hidden literals or on an
  unverified black-box coefficient handoff.

### Stage 024

- Canonical stage: `paper/stages/stage_024.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `sqrt(15/(8 pi))`
  derived as the normalization factor for the real STF `\ell = 2` harmonics
- `4 pi / 15`, `4 pi / 122`
  derived unit-sphere fourth and sixth moments used in the exact overlap
  contractions
- `sqrt(5)/(7 sqrt(pi))`
  derived from the exact `Y_20` triple-overlap matrix
- `1`, `1/2`, `-1`
  derived grouped axisymmetric splitting signature from that overlap matrix
- `1/4`, `3/4`
  derived grouped defect weights implied by the `20/21/22` signature
- there are no probe-only decimals in the theorem path

Audit note:

- The earlier notation caveat is stale. The current note explicitly says that
  `H_r = G_{U,r}^2 + G_{W,r}^2` is just the Stage-6 combined gauge/mixed
  strength written with a new letter to avoid collision with Newton's `G`.
- The earlier “tautological Section II” caveat is also stale. Both CAS layers
  now include unequal-lane witness checks showing that grouped defects become
  nonzero off the isotropic locus.
- The stage still treats the radial/axial overlap layer as carried reduced data
  from Stage 6, but that scope limit is explicit in the note and canonical
  stage card. The checkpoint claim itself is the angular closure and splitting
  law, and that claim is now clean.

### Stage 036

- Canonical stage: `paper/stages/stage_036.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `9`, `11`, `8`
  carried forward with source anchor from the Stage-18 selected-branch
  D/N-normal-form coefficients and the `8 / (pi^2 A)` dimensionless scaling
- `18`, `81`
  derived algebraically inside the derivative and endpoint reductions from the
  same closed form `G(xi,delta)`
- `2/9`
  derived in the near-onset expansion coefficient `-2 xi^2 / (9 delta)`
- there are no probe-only decimals in the theorem path

Audit note:

- The earlier boundary-assumption caveat is stale. The current SymPy audit uses
  `xi >= 0` via `nonnegative=True`, and the Mathematica mirror assumes
  `0 <= xi < 1` explicitly.
- Both CAS layers now check the same live theorem packet: exact `G`, exact
  loading split, manifestly positive derivative form, onset value, upper
  endpoint, and near-onset series.
- The checkpoint therefore no longer depends on a hidden `xi > 0` restriction
  that would narrow the onset-boundary claim.

### Stage 089

- Canonical stage: `paper/stages/stage_089.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `4/3`, `1/3`
  carried exact minimal-isotropic loading and support demands from the Stage-71
  precursor packet
- `3.46622291347846`, `3.46752913273870`, `3.46752922945601`
  carried Family-1 loading-ratio window markers from the upstream support-window
  packet
- `2.46752922945601`
  carried Family-1 hard support ceiling from the Stage 63/64 ceiling packet
- `1.00005192880220`
  carried zero-bias Family-1 baseline from the Stage-62 transport map
- there are no probe-only decimals in the theorem path

Audit note:

- This checkpoint is a closed arithmetic theorem conditional on the upstream
  minimal-isotropic module. It does not need to rederive the earlier support
  windows in order to verify the explicit zero-bias Family-1 verdict.
- Both CAS layers keep the carried thresholds explicit and check the exact
  ordering and margin claims without hidden literals.
- The resulting theorem claim is exactly the one stated in the stage card:
  the explicit Family-1 branch succeeds already at `Pe_req = 0` for the
  minimal isotropic passive/outgoing demand.

### Stage 090

- Canonical stage: `paper/stages/stage_090.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage090_updated_reduced_status_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage090_updated_reduced_status_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- `3/4`, `1/4`
  carried from the minimal isotropic conservative module fixed upstream in the
  grouped-`P_2` packet
- `4/3`
  derived in the audit as `1 / (3/4)`
- `1/3`
  derived in the audit as `(1/4) / (3/4)`
- `3.46622291347846`
  carried forward as `rho_suff^(chi)` from the Stage 69 support-window audit
- `2.46752922945601`
  carried forward as `zeta_max^(F1)` from the Stage 63/64 support ceiling
- `1.00005192880220`
  carried forward as `A_F1` from the Stage 62 transport map

Audit note:

- This script is intentionally a status-consistency audit. It does not smuggle
  those decimals in as unexplained literals; it declares them as carried
  threshold data and uses them only for the branch-ordering checks that Stage 73
  is supposed to summarize.
- Because the checkpoint claim is itself only that reduced status boundary, and
  the carried inputs are explicit and source-anchored, this narrow audit surface
  is sufficient for citation support.

### Stage 096

- Canonical stage: `paper/stages/stage_096.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `6`
  derived in the audit as the Laplace-sphere eigenvalue `ell(ell+1)` at
  `ell = 2`
- `1/4`
  derived from the obstruction formula once `eps_2 = eps_4 = 0`
- `3/4`
  derived as `1 - c_pole`
- `4/3`
  derived as `1 / c_geom`
- `1/3`
  derived as `c_pole / c_geom`

Audit note:

- The script rechecks the isotropic `l=0 <-> l=2` decoupling directly instead
  of treating `eps_2 = eps_4 = 0` as a free status assertion.

### Stage 163

- Canonical stage: `paper/stages/stage_163.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl`
- Verdict: `clean for current symbolic scope`

Constants reviewed:

- no free constants in the symbolic theorem path
  the parent compensation defect, `delta_perp` transport, microscopic
  parent-variable formula, outlet-defect transport, and tangent/normal split
  are all checked symbolically in both CAS layers
- `1.77799353547498`, `0.758035078944663`, `4.651033550168876`,
  `0.6703621156734617`
  carried Family-1 readbacks for `r_*`, `g_*`, `Sigma0_can`, and `S_can`,
  used only in the final numerical coefficient banner
- the derived readback coefficients
  `4 sqrt(1+r_*^2)`, `-1/sqrt(1+r_*^2)`, `1/(4 sqrt(1+r_*^2))`,
  `Sigma0_can S_can / sqrt(1+r_*^2)`, and `16 / sqrt(1+r_*^2)`
  are explanatory numeric readbacks, not theorem inputs

Audit note:

- The current symbolic theorem path does not depend on the Family-1 readback
  packet. Those values are only used to print the final numerical transport
  coefficients for the canonical point.
- Because the readbacks are explanatory rather than proof-critical, they are a
  provenance item to track, not a remaining trust defect.
- The checkpoint therefore no longer needs a separate hardening gate before it
  can count as citation support.

### Stage 185

- Canonical stage: `paper/stages/stage_185.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage185_187_orbit_stress.py`
  and
  `mathematica/numerical/stage185_187_orbit_stress.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1 + deltaU_*`, `1 + chi0_*`
  carried symbolic exponents from the Stage 183 tracking compiler
- `E_*`, `F_*`
  carried symbolic coefficients from the Stage 183 nontracking compiler; not
  numerically specialized in either CAS layer
- `2`, `11`, `9`, `4`
  source-anchored coefficients appearing in the carried Stage 183 formulas for
  `E_*` and `F_*`, not unexplained literals
- there are no free decimal literals in the symbolic theorem path

Audit note:

- The Stage 185 audits now reconstruct the primitive microscopic ratios
  `(gamma, c_{etaU}, T_U, K_U, K_eta^{(eff)}, K_W^{(eff)}, lambda_W, mu_W)`
  before rebuilding `chi_0`, `delta_U`, `epsilon_W`,
  `Z_W/Omega_W^2`, and `epsilon_eta`.
- The tracking, nontracking, and dressing monomial laws are then checked both
  from those primitive-ratio compilers and from the direct monomial ratios,
  removing the earlier tautology concern.
- The existing `185--187` numerical stress layer remains secondary; the main
  theorem path is now symbolic in both CAS layers.

### Stage 200

- Canonical stage: `paper/stages/stage_200.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1`
  structural zero-defect target carried by the Packet-A finish line
  `chi_Q = 1`
- `2`
  structural packet-length and pairing coefficient carried by the exact
  two-point orbit transport / cocycle formulas
- `3`, `5`, `9`
  carried from the Stage 200 definition
  `chi_Q = 3 (S beta^5 + 9 Sigma_5) / (3 S - Sigma_0)` and therefore source
  anchored rather than inserted ad hoc

Audit note:

- The Stage 200 audits are symbolic and contain no free decimal literals.
- The added Mathematica mirror checks the same reference-free packet identities,
  mismatch conversions, and linearized compiler in a second CAS.

### Stage 203

- Canonical stage: `paper/stages/stage_203.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage203_free_quintuple_scalar_closure_slice_and_crossing_theorem_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1`
  structural target value carried by the scalar closure condition
  `widehat chi_Q = 1`
- `2`
  structural coefficient carried by the Stage 192 quotient matrix and the
  graph-tangent / repair formulas
- `(1 + deltaU_*) / (1 + chi0_*)`
  carried symbolically from the Stage 192 same-free-quintuple graph formulas
- `E_*`, `F_*`
  carried symbolic coefficients from the monomial compiler; not numerically
  specialized in either audit

Audit note:

- The Stage 203 audits are symbolic and contain no free decimal literals.
- The added Mathematica mirror checks the graph-kernel theorem, graph-error
  packet compiler, inverse compiler, and repair vector in a second CAS.

### Stage 218

- Canonical stage: `paper/stages/stage_218.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `5`
  derived in the audit as the number of primitive free axes
  `{\lambda, c, \gamma, U, W}`
- `30`
  derived in the audit as the total count of nonempty proper support strata
  `Sum[Binomial[5,k], {k,1,4}] = 2^5 - 2`
- `1140`
  carried forward as the Stage 215 support-`<=4` global ledger budget
- `179`
  carried forward as the Stage 217 lifted support-5 per-envelope Bezout bound
- `750`
  carried forward as the Stage 217 projected-chart fallback per-envelope bound
- `324`, `1500`, `1464`, `2640`
  derived in the audit from the carried Stage 215 and Stage 217 budget packets

Audit note:

- The Stage 218 audits are symbolic / combinatorial and contain no free decimal
  literals.
- The added Mathematica mirror checks the boundary-identification counts,
  splice theorems, and carried budget arithmetic in a second CAS.

### Stage 221

- Canonical stage: `paper/stages/stage_221.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1/2`
  derived in the audit as the exact maximum of the conservative shape factor
  `r / (1 + r^2)` at `r = 1`
- `eta / (1 + eta^2)`
  derived in the audit as the exact low-loss envelope factor on the domain
  `r >= 1 / eta`
- `1 / (2 Q_* eta)`
  derived in the audit by substituting `gamma_* = omega_* / (2 Q_*)` into the
  low-loss detuning threshold
- `7`, `2`, `3`, `1/4`, `5`, `11`, `40`, `80`, `3/5`
  labeled probe-only numerical values used only for the constructive sign/scale
  sanity slice, not for the symbolic proof path

Audit note:

- The theorem path is symbolic in both CAS layers; there are no free decimal
  literals.
- The added Mathematica mirror checks the simple-pole normal form, the carried
  Stage 220 derivative identity, the wall-like specialization, the exact
  line-shape tradeoff laws, and the linear survival-window formulas.

### Stage 239

- Canonical stage: `paper/stages/stage_239.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `I_2`
  derived in the audit as the exact diagonal rigid-mouth packet compiler in the
  physical `(U,V)` chart
- `(0, -V, U - V)`
  derived in the audit from the Stage 236 dependent-plane carrier after the
  substitutions `q_nt = U`, `q_eta = V`
- `(0, 1, 1)`
  carried forward as the Stage 236 equal-drift dressing ray and used explicitly
  as a sourced packet direction rather than as a hidden literal
- `U = 0`, `V = 0`
  derived in the audit as the exact Cartesian orbit-lock conditions in the
  physical logarithmic chart

Audit note:

- The Stage 239 audits are symbolic and contain no probe numerics or free
  decimal literals.
- The added Mathematica mirror checks the diagonal physical chart, the exact
  physical-to-microscopic correction compiler, the support-blindness statements,
  and the Cartesian orbit-lock theorem in a second CAS.

### Stage 242

- Canonical stage: `paper/stages/stage_242.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage242_actual_twin_support_placement_and_coherent_orbit_lock_compiler_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `2 / 11`
  carried forward from the coherent local D/N placement map used in the Stage
  242 physical support variable
- `4 / 3`
  carried forward from the Stage 240 selected support demand
  `Pi_tr = (4 / 3) C_mix`
- `8 / Pi^2`
  carried forward from the Stage 240 mixed-support carrier
  `C_mix = 8 Lambda (1 - epsilon) / Pi^2`
- `2 / 3`
  derived in the audit by substituting the selected support demand into
  `varrho_phys = Pi^2 Pi_tr / (16 Lambda)`
- `1 / (2 + beta^2)` and `beta / (1 + beta + beta^2)`
  derived in the audit by rewriting the Stage 241 thresholds in the realized
  support variable `epsilon`
- `2 deltaU / ((1 + deltaU) (11 + 9 deltaU))`
  derived in the audit by differentiating the carried support variable
  `epsilon = epsilon_W (1 - (2 / 11) deltaU / (1 + deltaU))`
- `3 / 2`, `2 / 3`, `13 / 17`, `1 / 3`, `1 / 5`, `7 / 11`, `1`
  labeled probe-only numerical values used only in the rational sanity sample,
  not in the symbolic proof path

Audit note:

- The Stage 242 theorem path is symbolic in both CAS layers; the explicit
  rational sample point is probe-only.
- The added Mathematica mirror checks the realized support coordinate, the
  threshold rewrites, the support-blind orbit packet, the infinitesimal
  observable compilers, and the exact support/orbit split in a second CAS.

### Stage 243

- Canonical stage: `paper/stages/stage_243.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `-Sqrt[2] / 4` and `Sqrt[2] Sqrt[Pi] / 8`
  derived in the audit by evaluating the exact Gaussian leakage and work
  integrals
- `1 / 2`
  used only as the standard quadratic normalization in the non-rigid free
  energy and carried transparently in both CAS scripts
- `chi_lam / k_V`
  derived in the audit as the exact non-rigid ratio `V / U`
- `-a / (4 b)` and `1 - b - a^2 / (8 b)`
  derived in the audit as the interior stationary point and vertex value of the
  compensated source quadratic
- `r^-6`, `Exp[-2 kappa r] / r^4`, `Exp[-4 kappa r] / r^2`
  carried forward from the one-port short-range same-charge kernel verdict and
  checked explicitly for vanishing `r * V(r)` tails

Audit note:

- The Stage 243 audits are symbolic in both CAS layers and contain no probe
  numerics.
- The added Mathematica mirror checks the exact Gaussian leakage/work channel,
  the linear `(U,V)` solve, the recovery slice, and the short-range limit
  firewall in a second CAS.

### Stage 248

- Canonical stage: `paper/stages/stage_248.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage248_event_chain_stress.py`
- Numerical stress (Mathematica):
  `mathematica/numerical/stage248_event_chain_stress.wl`
- Verdict: `clean for current checkpoint scope`

Constants reviewed:

- `1 / 2`
  carried transparently as the standard kinetic-energy normalization in the
  reduced one-dimensional event chain
- `1 / E`
  derived in the audit as the exact pure-Coulomb outer turning point
- `Pi`
  derived in the audit through the exact near-top parabolic action integral
- `2.5413906350657705`, `3.272783388968954`, `2.1447620194324593`, `0.4`,
  `0.673752615`, `0.546377065`, `1.23312756`, `23.3128`
  derived benchmark outputs checked against the declared Session-II benchmark
  specialization
- `5.0`, `0.18`, `2.5`, `0.19999794`, `3.42933112`, `0.23944389`,
  `0.39096144`, `0.19039548`, `0.19744614`, `0.30222297`, `0.34437471`,
  `0.42826825`, `2.59221845`, `0.28091705`
  labeled benchmark-only numeric inputs copied from the declared Session-II
  readback, not used as part of the symbolic theorem path
- In `stage248_event_chain_samples.json`, the probe families are parameterized
  by barrier gap, subbarrier fraction, cross-window fraction, `I_new /
  I_Coul`, `\Xi_{\rm turn}`, and `|V'(r_{\rm turn})|`.
  Those are probe-only stress inputs.
- The stress harness then derives `V_{\rm peak}`, `E_{\rm sub}`,
  `v_{\rm crit,new}`, `v_{0,\rm sub}`, `v_{\rm cross}`,
  `r_{\rm turn,Coul}`, `\lambda_{\rm th}`, and the transmission ratio from the
  exact Stage 248 formulas instead of hard-coding those outputs.
- The near-top stress block uses probe-only `(\Delta E, K_{\rm peak}, m_s,
  \hbar_{\rm eff})` tuples and compares direct quadrature against the exact
  parabolic top-action formula.

Audit note:

- The Stage 248 theorem path is symbolic in both CAS layers up through energy
  conservation, threshold-speed compilation, Coulomb reference formulas, and
  near-top action.
- The Session-II numbers are intentionally confined to a benchmark-only
  specialization layer; they are tracked as declared readback values rather
  than hidden derivation inputs.
- The new shared numerical-stress layer now checks three admissible event-chain
  probe families plus three near-top probe families in both CAS layers.
- The Coulomb WKB action is numerically integrated and compared against the
  exact closed form in both CAS layers; the near-top action is also compared
  against direct quadrature in both CAS layers.

### Stage 253

- Canonical stage: `paper/stages/stage_253.tex`
- SymPy audit:
  `scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py`
- Mathematica audit:
  `mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- Numerical stress:
  `scripts/numerical/stage253_material_threshold_stress.py`
- Numerical stress (Mathematica):
  `mathematica/numerical/stage253_material_threshold_stress.wl`
- The lattice-turnover compiler, legacy recovery slice, harmonic
  trigger/stiffness compiler, Korringa ceiling, and screening ratios are all
  derived symbolically in both CAS layers.
- `\Upsilon_{\rm lat}` remains an explicit calibration parameter; it is not
  solved for or hidden inside a hard-coded literal.
- `4.79562976` is a carried forward legacy Session-V reduced lattice-rate input
  used only in the declared legacy-slice recovery.
- `0.5489386551062235`, `6.94311167`, `0.75`, and `1.0` are benchmark-only
  Stage 252 slice inputs for `s_c`, `s_0`, `f_{\rm lat}`, and `\mu_\eta`.
- `0.39096144`, `0.42826825`, and `2.73855812` are benchmark-only
  turning-point / force-match inputs carried from the declared Session-V
  benchmark slice.
- `65.45193925961132`, `13.64824695299483`, `136.23361317476524`,
  `8.736185210116078`, `0.9128891530016525`, `2.1908464937104797`,
  `10.95423248`, `5.837462857946154`, and `0.5489386551062235` are
  benchmark-derived outputs computed from those declared inputs.
- The Stage 253 theorem path is symbolic in both CAS layers. All explicit
  decimals are benchmark-only readbacks, not hidden theorem inputs.
- In `stage253_material_threshold_samples.json`, the declared
  `Pi_ep/Pi_chi/Pi_k/Pi_t` targets are probe-only screening margins. The
  stress harnesses derive `\lambda_{\rm ep}\omega_D`, `r_{\rm turn}`,
  `k_{\rm eff}`, and `T` from the exact Stage 253 thresholds rather than
  hard-coding those candidate values directly.
- `\Upsilon_{\rm lat}` in the stress layer is either the raw microscopic slice
  `1`, the exact legacy replay `\gamma_{\rm lat,safe}^{\rm eq} /
  \gamma_{\rm lattice}^{\rm red}`, or an explicit custom calibration probe.

## Open Follow-up

- Symbolic parity is now closed across the checkpoint support set.
- The stage-level trust baseline now lives in `CHECKPOINT_TRUST_AUDIT.md`.
- Dedicated numerical hardening is now in place for checkpoint stages `248`
  and `253`.
- The next checkpoint pass should decide whether to widen numerical stress
  beyond the current four stress-backed checkpoints (`003`, `185`, `248`,
  `253`) and how to expose that verification surface compactly in the paper.
