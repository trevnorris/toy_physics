# PathA 22b minimal combination xi

## Gate 0

### 0a mhat reconciliation

- Outcome: `MHAT_DIMENSIONFUL_CONFIRMED`.
- Reasoning: The normalization law itself forces a dimensionful source-map scale. The later natural-map sentence is dimensionless, so it cannot be the whole object in that law. The consistent reading is mhat = mhat0*g_mhat, with mhat0 carrying units and g_mhat the dimensionless profile/source-map factor.
- Effect on pathA_22a: No flip to the pathA_22a decomposition: mhat0^2 remains the scale carrier and g_mhat remains a dimensionless branch/source-map residual until a natural-source-map derivation pins it.
- Provenance: research/pde/paper/pde.tex:2053-2062; research/pde/paper/pde.tex:2075-2083; research/pde/paper/pde.tex:2095-2099; software/stage1_solver/src/stage1_solver/dimensional_check.py:4280-4288.

| ingredient | dimension | SymPy monomial | derivation | provenance |
| --- | --- | --- | --- | --- |
| `Gamma5` | `T^5` | `T_exp**5` | Odd outgoing coefficient dimension from the explicit a^5/c_s^5 carrier. | (c_s/a)^2-normalized (dimensionless) P0 per pathA_22a; research/pde/paper/pde.tex:2053-2062 |
| `G/c^5` | `L^-2 T^3 M^-1` | `T_exp**3/(L_exp**2*M_exp)` | Right side of the invariant odd-coefficient normalization. | research/pde/paper/pde.tex:2075-2079 |
| `mhat0` | `L^-1 T^-1 M^-1/2` | `1/(L_exp*sqrt(M_exp)*T_exp)` | Scale carrier forced by [mhat0]^2*[Gamma5]=[G/c^5]. | research/pde/paper/pde.tex:2075-2083; software/stage1_solver/src/stage1_solver/dimensional_check.py:4280-4288 |
| `natural source-map factor` | `1` | `1` | Dimensionless profile/source-map factor; the correction a^2/r^2 is dimensionless. | research/pde/paper/pde.tex:2095-2099 |

0a checks:
- `Gamma5 outgoing coefficient`: **CONSISTENT** (expected `T^5`, actual `T^5`). Gamma5=chi_Q*P0*a^5/(const*c_s^5) with (c_s/a)^2-normalized (dimensionless) P0 per pathA_22a and dimensionless chi_Q.
- `BT odd-coefficient right side G/c^5`: **CONSISTENT** (expected `L^-2 T^3 M^-1`, actual `L^-2 T^3 M^-1`). The invariant odd-coefficient target carries the dimension that mhat^2*Gamma5 must match.
- `required mhat from mhat^2*Gamma5=G/c^5`: **CONSISTENT** (expected `L^-1 T^-1 M^-1/2`, actual `L^-1 T^-1 M^-1/2`). Solving the unit equation fixes the source-map scale carrier.
- `factorized normalization right side`: **CONSISTENT** (expected `L^-2 T^-2 M^-1`, actual `L^-2 T^-2 M^-1`). After substituting Gamma5, the right side has the dimension of mhat0^2.
- `dimensionful law`: **CONSISTENT** (expected `L^-2 T^3 M^-1`, actual `L^-2 T^3 M^-1`). This is the non-tautological unit reconciliation for the odd-coefficient law. terms=mhat0^2*Gamma5:L^-2 T^3 M^-1, G/c^5:L^-2 T^3 M^-1
- `direct dimensionless mhat reading in odd-coefficient law`: **INCONSISTENT** (expected negative) (expected `L^-2 T^3 M^-1`, actual `T^5`; factor needed `L^-2 T^-2 M^-1`). This is expected to fail if the dimensionless natural-map sentence is read as the whole mhat.
- `natural source-map correction a^2/r^2`: **CONSISTENT** (expected `1`, actual `1`). The natural-map sentence is dimensionless and therefore belongs to the profile/source-map factor.

### 0b Z-kernel cancellation source assessment

- Outcome: `DOES_NOT_CANCEL (NOT_ESTABLISHED — sources do not establish either cancellation route; a later Gate-4 action-level derivation could still find one)`.
- Scope: Gate 0 is a source-availability assessment; the action-level derivation of the `g_G`/`g_mhat²` kernels is deferred to Gate 4.
- Exact structural condition: Z-independence of g_mhat^2/g_G can be established by either route (a) shared factored scalar, where both g_G and g_mhat^2 are I_Z times separate w-independent field-content kernels with the same I_Z=int sqrt(g_w)*Z(w) dw, or route (b) pointwise-proportional kernels, where K_stress(w)=const*K_source(w) for all w. Route (b) is equivalently the identity K_stress(w)*K_source(v)=K_stress(v)*K_source(w) for all w,v, making the weighted-average ratio independent of the Z/W_eff profile.
- Do the sources establish either route? `False`. Gate 0 is a source-availability assessment; the action-level derivation of the g_G/g_mhat^2 kernels is deferred to Gate 4. The parent action places Z(w) in the localized Maxwell kinetic term and Maxwell equation, while the GNLS matter sector, exact current, source-map statements, and brane projection W(w) are distinct. The cited sources do not establish route (a), a shared factored I_Z functional for g_G and g_mhat^2, nor route (b), proportional stress/source kernels.
- Gate 4 implication: Before declaring W_eff/full transverse profile irreducible, Gate 4 must test both cancellation routes: route (a) shared factored scalar and route (b) pointwise-proportional kernels. Unless Gate 4 proves one of those routes from the action-level kernels, W_eff/full transverse profile remains on the critical path for the ratio.
- Provenance: research/pde/paper/pde.tex:277; research/pde/paper/pde.tex:289-295; research/pde/paper/pde.tex:357-416; research/pde/paper/pde.tex:496-565; software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md:80-104; software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md.

0b checks and controls:
- `common scalar functional algebra`: **CONSISTENT** (expected `K_stress/K_source`, actual `K_stress/K_source`). If both factors are proven to contain the same factored I_Z, the scalar cancels algebraically.
- `non-factorizing weighted negative control`: **CONSISTENT** (expected `DOES_NOT_CANCEL`, actual `DOES_NOT_CANCEL`). Distinct kernels in int Z*K_stress / int Z*K_source leave dependence on Z.
- `weighted proportional control`: **CONSISTENT** (expected `CANCELS`, actual `CANCELS`). Pointwise proportional weighted kernels make int Z*K_stress / int Z*K_source independent of Z; this is route (b), distinct from shared scalar route (a).
- Negative control: `w + 1` vs `w**2 + 1` gives `DOES_NOT_CANCEL`; condition residual `(v - w)*(v*w + v + w - 1)`.
- Proportional-kernel control: `2*w**2 + 2` vs `w**2 + 1` gives `CANCELS` by route (b); condition residual `0`.
- Positive algebra control: common scalar functional gives `CANCELS` with ratio `K_stress/K_source`.

### Dimensional checks

- `mhat0^2*Gamma5 vs G/c^5`: **CONSISTENT** (expected `L^-2 T^3 M^-1`, actual `L^-2 T^3 M^-1`). terms=mhat0^2*Gamma5:L^-2 T^3 M^-1, G/c^5:L^-2 T^3 M^-1
- `I_Z=int sqrt(g_w) Z(w) dw scalar carrier`: **CONSISTENT** (expected `L`, actual `L`). The exact length dimension depends on the w-coordinate convention; it cancels only after a proven common-factor derivation.

### Target-blindness

- `TARGET_BLIND_PASS` over the new Gate-0 module and Mathematica cross-check.

### Residual ledger

- The natural dimensionless source-map factor g_mhat is not dynamically pinned by Gate 0.
- The parent sources do not establish either route (a) a shared factored I_Z for g_G and g_mhat^2 or route (b) pointwise-proportional kernels.
- Gate 4 cannot use a reduced field-content ratio unless a later action-level derivation proves one of the two cancellation routes; otherwise it needs W_eff/full transverse kernels.

## Gate 1

### Measure determination

- PDE definition used: `flat_dw_no_sqrt_g_w`. `pde.tex:289-295` defines `Z_int=int Z(w) dw`; it does not include a `sqrt(g_w)` factor.
- Discrepancy flag: Gate 0's structural carrier included sqrt_g_w, but pde.tex defines Z_int with flat dw. The frozen export provides no independent sqrt_g_w profile; in the flat/densitized convention sqrt_g_w=1 and the numeric quadrature is identical.
- Numeric `sqrt(g_w)` variant: `NOT_INDEPENDENTLY_COMPUTABLE_FROM_EXPORTED_Z_W`. If the export is read in the flat/densitized convention, the variant equals `2.03111437236 L`.
- Exported grid measure label: `4*pi*r^2 dr dw`; w quadrature used the exported `w_faces` widths.

### Quadrature result

- Headline: `Z_int ~= 2.03111437236 L` (finite-box, floor-dominated, domain-dependent; ideal +/-infinity integral divergent).
- Rule: cell-centered midpoint sum over exported w_faces widths. Primary background: `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json`.
- Domain/grid: `w in [0, 1.85]`, `nw=8`, `dw=0.23125`.
- Resolution diagnostic: nearest-grid change `0.00111874036259 L` (grid-resolution delta only -- NOT total uncertainty); full 4/6/8-point ladder spread `0.00429934874052 L`.
- Tail/truncation diagnostic: edge values `Z_left=0.94042125771`, `Z_right=0.94042125771`, max edge/peak ratio `0.75875254043`.
- Edge-cell contribution diagnostic: the two boundary cells contribute `0.434944831691 L`, fraction `0.214140984678` of the finite-domain integral.
- Ideal infinite-domain status: `BLOCKED_NEEDS_DECAYING_OR_ZERO_FLOOR_EXTENSION`. The exported nonzero floor means the finite-domain result is the faithful exported-grid integral, but the omitted ideal tail is not bounded as a small correction by this export.
- Floor decomposition: `localization_floor=0.8` over domain width `1.85 L` contributes `1.48 L`, which is `72.9%` of finite-box `Z_int`; the localized Gaussian part contributes `0.551114372359 L` (`27.1%`), so the Gaussian is the minority.
- Gaussian-only omitted-tail diagnostic, ignoring the nonzero floor extension: `0.04852910805 L`.

### Dimensions and scope

- `exported Z_w localization weight`: **CONSISTENT** (expected `1`, actual `1`). The code exports Z_w from dimensionless floor/amplitude constants and a dimensionless Gaussian argument.
- `Z_int=int Z(w) dw`: **CONSISTENT** (expected `L`, actual `L`). pde.tex defines Z_int with flat dw, so the integral carries the w-coordinate length dimension.
- Scope: Gate 1 reports Z_int only as a coupling-normalization artifact: mu0_eff=mu0/Z_int and q_eff=q_star/sqrt(Z_int). Z_int is not a factor in P0*chi_Q*g_mhat^2*lambda_gamma^5/g_G, does not gate the xi verdict, and does not promote Z_int to lambda_gamma or alter photon-cone speed. If needed downstream, carry it symbolically as mu0_eff=mu0/Z_int and q_eff=q_star/sqrt(Z_int), never as the numeric finite-box value.

### Provenance

- Measure: research/pde/paper/pde.tex:289-295.
- Localized Maxwell source of `Z(w)`: research/pde/paper/pde.tex:357-416.
- Coupling reduction: research/pde/paper/pde.tex:541-563.
- Exported `Z_w`: software/stage1_solver/src/stage1_solver/m1c_background_export.py:166-170.
- Frozen data: `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json`.

### Target-blindness

- `TARGET_BLIND_PASS` over the new Gate-1 module and Mathematica cross-check.

### Residual ledger

- The exported primary grid is only nw=8; the nearest-grid finite-domain change is recorded as a grid-resolution delta only, not as the total uncertainty.
- The exported Z_w is not small at the domain edges, so the finite interval is not a controlled small-tail approximation to an ideal infinite integral.
- Because the frozen source has a nonzero localization_floor, an infinite continuation of that same formula would not give a finite Z_int; an external decaying/zero-floor continuation would be needed to bound the omitted tail.
- Floor provenance: localization_floor=0.8 is an undocumented solver config constant (coupled_branch.py:188-192; m1c_physical_run.py:121-123), differs across presets (patha_closed_newton.py:61-63), and has no source support; pde.tex:277,289-295 uses a localized Z(w) over (-infinity,+infinity). Next step to unblock: export/derive a genuinely decaying Z(w), or obtain documented physical provenance for the floor plus a physical w-extent.
- No photon-cone speed or lambda_gamma normalization is derived or modified in Gate 1.

- Gate 1 outcome: `BLOCKED_NEEDS_DECAYING_Z_PROFILE_OR_FLOOR_PROVENANCE`.

## Gate 2

### Existing pathA_20b result

- Group/module consumed: `software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:24-50`. It established the transverse gauge principal block `P_T=C_E*omega^2-C_B*k^2`, hence `c_bulk^2=C_B/C_E`, and carried `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`.
- Gate-2 scope: reduce the localized zero-mode Maxwell action onto the brane, form the photon cone, and treat the missing speed normalization explicitly. No Gate 3-5 work is used.

### C_E/C_B reduction

- Metric finding: `FLAT_UNWARPED_BULK_METRIC` with `eta_MN=diag(-1,+1,+1,+1,+1)`. Provenance: `research/pde/paper/pde.tex:242-244`.
- Source action: `L_EM=-(Z(w)/(4*mu0))*F_MN F^MN` plus gauge fixing and external/background source terms. Provenance: `pde.tex:357-416`.
- Reduction assumptions: `A_w approximately 0, partial_w A_mu approximately 0, J^w approximately 0, F_mu_w approximately 0`. Provenance: `pde.tex:543-552`.
- Measure: `flat_dw_no_sqrt_g_w` because `pde.tex:289-295` defines `Z_int=int Z(w) dw`; no independent `sqrt(g_w)` factor is introduced.
- Flat-metric expansion: `F_MN F^MN = -2*F_tx**2 - 2*F_ty**2 - 2*F_tz**2 + 2*F_xy**2 + 2*F_xz**2 + 2*F_yz**2` after the zero-mode assumptions, i.e. `-2*F_tx**2 - 2*F_ty**2 - 2*F_tz**2 + 2*F_xy**2 + 2*F_xz**2 + 2*F_yz**2`.
- Reduced Maxwell density: `F_tx**2*Z_w/(2*mu0) + F_ty**2*Z_w/(2*mu0) + F_tz**2*Z_w/(2*mu0) - F_xy**2*Z_w/(2*mu0) - F_xz**2*Z_w/(2*mu0) - F_yz**2*Z_w/(2*mu0)`; after the flat `dw` integration: `F_tx**2*Z_int/(2*mu0) + F_ty**2*Z_int/(2*mu0) + F_tz**2*Z_int/(2*mu0) - F_xy**2*Z_int/(2*mu0) - F_xz**2*Z_int/(2*mu0) - F_yz**2*Z_int/(2*mu0)`.
- Reduced coefficients: `C_E=Z_int/mu0`, `C_B=Z_int/mu0`.
- Cone ratio: `C_B/C_E=1`. This is a computed flat-metric result, not a second free knob; both coefficients are `Z_int/mu0`, so the common localization factor cancels from the characteristic speed.
- Transverse operator: `-Z_int*(k2 - omega**2)/mu0`.

### Speed Normalization

- Carried speed map: `beta_bulk_to_brane` with status `UNPINNED_BY_SOURCES`.
- Genuine open knobs after the F^2 computation: `1` (`bulk_metric_to_acoustic_speed_identification`).
- Physical photon cone carried by Gate 2: `c_gamma^2 = beta_bulk_to_brane^2` because `C_B/C_E=1` in bulk-metric units.
- GNLS sound speed: `c_s^2 = 5*K*rho0**4/m_GNLS`.
- Result: `lambda_gamma = beta_bulk_to_brane*sqrt(m_GNLS/(5*K*rho0**4))`.
- Outcome: `STILL_TUNABLE_LAMBDAGAMMA`. This is still tunable through one named speed-normalization residual, not a two-condition geometry result.
- Gate-5 protection: carry `lambda_gamma` as an explicit open knob until the speed map is pinned; the downstream xi verdict is `FAIL_ABLE_PENDING_LAMBDAGAMMA_SPEED_MAP`, not `REAL_MATCH`.

### Negative Control

- Verdict: `GUARD_AGAINST_FORCING_LAMBDAGAMMA_TO_ONE`.
- Forced equality would require `beta_bulk_to_brane^2=5*K*rho0^4/m_GNLS`.
- Symbolic residual carried: `-5*K*rho0**4/m_GNLS + beta_bulk_to_brane**2`.
- Residual is an identity zero? `False`. This is a guard against forcing `lambda_gamma=1`, not evidence from a discriminating source equation.

### Dimensional Checks

- `C_E electric principal coefficient`: **CONSISTENT** (expected `L^-4 T^2 M^-1`, actual `L^-4 T^2 M^-1`). Defined by the coefficient multiplying E_i E_i in the reduced Maxwell action.
- `C_B magnetic principal coefficient`: **CONSISTENT** (expected `L^-4 T^2 M^-1`, actual `L^-4 T^2 M^-1`). The flat bulk metric gives the same scalar-weighted coefficient as the electric term.
- `C_B/C_E flat bulk cone ratio`: **CONSISTENT** (expected `1`, actual `1`). In bulk-metric units the localized Maxwell zero mode has c_bulk^2=1.
- `speed-normalization beta_bulk_to_brane`: **CONSISTENT** (expected `L T^-1`, actual `L T^-1`). This is the remaining source-unpinned map from unit bulk-metric speed to physical brane speed.
- `c_gamma=beta*sqrt(C_B/C_E)`: **CONSISTENT** (expected `L T^-1`, actual `L T^-1`). Physical brane photon speed after the flat reduction fixes C_B/C_E=1.
- `c_s=sqrt(5*K*rho0^4/m_GNLS)`: **CONSISTENT** (expected `L T^-1`, actual `L T^-1`). GNLS sound speed from the already-derived pathA_19/20 EOS law.
- `lambda_gamma=c_gamma/c_s`: **CONSISTENT** (expected `1`, actual `1`). The single remaining speed-normalization dial is dimensionless only as beta_bulk_to_brane/c_s.
- `bulk-metric transverse Maxwell principal operator`: **CONSISTENT** (expected `L^-5 T`, actual `L^-5 T`). The flat parent metric uses x^0 and spatial coordinates in the same bulk-metric units. terms=C_E*partial_0^2 A_T:L^-5 T, C_B*laplacian A_T:L^-5 T
- `photon and sound dispersions`: **CONSISTENT** (expected `T^-2`, actual `T^-2`). This confirms comparable units only; it does not assert equal coefficients. terms=omega^2:T^-2, beta^2*(C_B/C_E)*k^2:T^-2, c_s^2*k^2:T^-2

### Provenance

- research/pde/paper/pde.tex:242-244.
- research/pde/paper/pde.tex:289-295.
- research/pde/paper/pde.tex:357-416.
- research/pde/paper/pde.tex:543-552.
- research/pde/paper/pde.tex:553-558.
- software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:24-50.

### Target-Blindness

- `TARGET_BLIND_PASS` over the new Gate-2 module, tests, and Mathematica cross-check. No final-comparison constants are used in the derivation.

### Residual Ledger

- The flat localized Maxwell reduction computes C_B/C_E=1, so C_B/C_E is no longer an open Gate-2 knob.
- The bulk-to-brane speed normalization beta_bulk_to_brane is not pinned by pde.tex:357-416 or by the zero-mode reduction at pde.tex:543-558.
- Z_int cancels from the cone and remains a coupling-normalization artifact, not lambda_gamma.
- A definite numerical lambda_gamma requires a source-derived identification of the unit bulk Maxwell cone with the acoustic brane speed scale.
- Gate 5 must carry lambda_gamma as an explicit open knob; the overall xi verdict remains FAIL_ABLE_PENDING_LAMBDAGAMMA_SPEED_MAP and must not be folded into a REAL_MATCH.

- Gate 2 outcome: `STILL_TUNABLE_LAMBDAGAMMA`.

## Gate 2 addendum -- beta status

- Classification: `BETA_GENUINE_GAP`.
- Decisive reason: the declared parent action fixes the localized Maxwell zero-mode cone only in the independent flat bulk metric, while the GNLS acoustic cone is fixed by the medium equation of state. The sources do not provide a postulate or field equation identifying the bulk metric's unit null speed with the GNLS acoustic speed scale.
- Maxwell side: the bulk metric is declared as `eta_MN=diag(-1,+1,+1,+1,+1)`, so the Maxwell kinetic term `-(Z(w)/(4*mu0))*F_MN F^MN` has equal electric and magnetic principal coefficients after the controlled zero-mode reduction. The reduction gives `C_B/C_E=1` in bulk-metric units, with `Z_int` cancelling from the cone and remaining only in `mu0_eff`/`q_eff`. Provenance: research/pde/paper/pde.tex:242-248, research/pde/paper/pde.tex:357-416, research/pde/paper/pde.tex:541-565, and Gate-2 reduction above at this report:118-126.
- GNLS side: the matter sector is a gauged GNLS medium with `P(rho)=K*rho^5` and `c_s^2=(1/m)*dP/drho=5*K*rho^4/m`. This derives the acoustic speed from the medium's EoS and mass parameter, not from the Maxwell metric. Provenance: research/pde/paper/pde.tex:316-352.
- Coupling status: the parent action states a combined matter--gauge system, but the common coupling is minimal charge/current coupling plus the localized Maxwell equation; it does not add a constitutive relation equating the gauge kinetic metric to the acoustic metric or tying `beta_bulk_to_brane^2` to `5*K*rho0^4/m_GNLS`. Provenance: research/pde/paper/pde.tex:318-340, research/pde/paper/pde.tex:370-416, research/pde/paper/pde.tex:903-931. The "fixed/not fixed" firewall also leaves later branch realization data open rather than promoting a metric-acoustic identification to parent-action content: research/pde/paper/pde.tex:567-577.
- Legacy-scope provenance: the emergent-acoustic EM identification in research/em_fields/paper/em_fields.tex:160,172,475-482, which would force `c_gamma=c_s`, is the superseded legacy construction; the canonical Path-A parent action re-postulates EM as an independent gauge field on the flat bulk metric, research/pde/paper/pde.tex:257-262, per software/stage1_solver/reports/pathA_20_velocity_constants.md:35 and software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:68, so `BETA_GENUINE_GAP` is scoped to the declared parent action and is not overturned by the legacy emergent route.
- Refinement of pathA_20b: pathA_20b's `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED` and `PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING` are not resolved by the Gate-2 `F^2` reduction. Gate 2 removes the extra `C_B/C_E` ambiguity by proving `C_B/C_E=1` in bulk-metric units, but it does not derive the remaining speed map to acoustic units. Provenance: software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:47-69 and this report:128-136,170-176.
- Not a mere convention: choosing units can set one speed's numerical value, but the dimensionless ratio `c_gamma/c_s` is invariant under common brane unit rescalings. Declaring the EM unit cone to equal the acoustic cone would therefore be an added physical identification, not a prediction already contained in the action.
- Missing input needed to derive it later: an explicit `bulk_metric_to_acoustic_metric_identification` postulate or branch-level constitutive law, for example a source-derived equation fixing `beta_bulk_to_brane^2 = 5*K*rho0^4/m_GNLS` for equality, or more generally fixing the EM bulk metric's time/length scale as a function of the GNLS medium parameters.
- Implication: `lambda_gamma = beta_bulk_to_brane*sqrt(m_GNLS/(5*K*rho0^4)) = beta_bulk_to_brane/c_s` is a calibration input until that missing physical input is supplied; it is not a parent-action prediction.

## Gate 3 — χ_Q

### Definition

- Scope: linear-response solve around the frozen M1c finite-core branch. No nonlinear profile re-solve and no deep/empty-throat solve were run.
- Frozen branch: `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json` with grid `{'nr': 10, 'nw': 8, 'r_range': [0.0, 2.0], 'w_range': [0.0, 1.85]}` and stationary residual `4.5892669664482355e-09`.
- Boundary labels: freeze sheet `mouth Dirichlet / open-impedance(Robin) exit / no hard cap`; packet `open_impedance` with exit model `impedance_mismatch_open_exit`.
- Hankel check: `Lambda_2^out=I*z**5/9 + z**4/9 + z**2/3 - 3` and `Y_out=I*z**5/27 + 4*z**4/81 + z**2/9 + 1`; the `z^5` coefficient is `1/27`.
- Anti-tautology: cached `Gamma_port=0.14814814814814817` is present only as legacy provenance and is not used by the solve.
- Source adapter: frozen derived BdG wall modes from the physical packet; the outgoing DtN/self-energy sweep itself is recomputed from the frozen background.
- Extraction: for each domain, grid, and frequency window, fit `E = C5` from `-Im(Sigma_out/Sigma_cons)/omega^5 = C5 + C7*omega^2 + C9*omega^4`, then set `chi_Q = E_defect / E_reference`.
- Physical reference: the branch-geometry reference keeps the exported `Z_w` and `R0_w`, and undefects only the condensate/gauge sector (`rho` uniform, `A -> 0`). The flat reference (`Z_w=1`, `R0_w=a`) is retained only to quantify the old reference choice.
- Exterior caveat: `radial_scale>1` uses the existing radial exterior taper, not a nonlinear far-field re-solve. The sweep tests/bounds boundary-placement sensitivity on this frozen continuation; lack of plateau still blocks non-characterized inputs.

### r_max Domain Sweep

- `dr` control: nr is recomputed from the scale-1 base dr so r_max changes do not deliberately coarsen the radial cell size.
- Branch-geometry reference tail mean `7.420107259365566e-01`, tail spread `1.394e-03`, relative tail spread `1.879e-03`, plateau `True`.
- Plateau onset: `r_max=6.0` by the `1.5e-03` tail-delta rule; strict `1e-3` final-lock onset `r_max=7.0`.
- Flat `Z_w=1` reference is a plateau-tail diagnostic over scales `[3.0]`: tail mean `7.038601320016022e-01`, tail spread `0.000e+00`, relative tail spread `0.000e+00`.

Branch-geometry reference:

| radial_scale | r_max | grid | dr | E_defect | E_reference | chi_Q | delta chi_Q | max residual |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 2 | 10x8 | 1.618182e-01 | 4.399618116589e-01 | 2.416783352202e+00 | 1.820443736746e-01 |  | 9.382e-16 |
| 1.5 | 3 | 16x8 | 1.635294e-01 | 1.936856040594e+00 | 5.588811643100e+00 | 3.465595486628e-01 | 1.645151749882e-01 | 9.319e-16 |
| 2 | 4 | 22x8 | 1.643478e-01 | 7.845910949754e+00 | 1.171598578464e+01 | 6.696756972891e-01 | 3.231161486263e-01 | 1.002e-15 |
| 2.5 | 5 | 29x8 | 1.593333e-01 | 2.384749254423e+01 | 3.253017802069e+01 | 7.330882889440e-01 | 6.341259165495e-02 | 1.065e-15 |
| 3 | 6 | 35x8 | 1.605556e-01 | 5.931640551817e+01 | 8.003232505849e+01 | 7.411555952525e-01 | 8.067306308464e-03 | 1.087e-15 |
| 3.5 | 7 | 41x8 | 1.614286e-01 | 1.282190805257e+02 | 1.727259937405e+02 | 7.423264891925e-01 | 1.170893940022e-03 | 1.101e-15 |
| 4 | 8 | 47x8 | 1.620833e-01 | 2.500120710280e+02 | 3.366938786516e+02 | 7.425500933647e-01 | 2.236041721781e-04 | 1.001e-15 |

Flat `Z_w=1` plateau-tail diagnostic:

| radial_scale | r_max | grid | dr | E_defect | E_reference | chi_Q | delta chi_Q | max residual |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 6 | 35x8 | 1.605556e-01 | 5.931640551817e+01 | 8.427300087233e+01 | 7.038601320016e-01 |  | 1.116e-15 |

### Z_w Reference Comparison

- Physically correct zero for the defect's radiative coupling: `branch_geometry` (`Z_w/R0_w` retained, condensate/gauge undefected).
- Matched converged-grid reference shift, branch minus flat: `6.081952123083167e-02` (relative `8.54%`). This is a one-sided definitional systematic because branch-geometry is the physical zero and flat `Z_w=1` is the less-physical alternative; it is not folded into the numerical error bar.

Matched converged even-nw rows:

| nr | nw | r_max | branch chi_Q | flat Z_w=1 chi_Q | branch-flat shift | relative shift |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 48 | 64 | 6 | 7.127073077159e-01 | 6.517756787792e-01 | 6.093162893678e-02 | 8.549% |
| 48 | 80 | 6 | 7.125120276544e-01 | 6.518009360986e-01 | 6.071109155583e-02 | 8.521% |
| 48 | 96 | 6 | 7.122748130763e-01 | 6.514589698764e-01 | 6.081584319989e-02 | 8.538% |

### Grid Convergence at Post-Plateau r_max

- Fixed post-plateau `r_max=6` (`radial_scale=3`) branch grid sequence values: `[0.7411555952524825, 0.7220737682611891, 0.7159043960619648]`.
- Observed delta ratios: `[0.3233114000058398]`.
- Richardson/geometric-tail limit: `7.129567649611814e-01`; finest value `7.159043960619648e-01`; old two-grid mean `7.189890821615770e-01`.
- Grid extrapolation uncertainty `6.032e-03`; raw finest-limit correction `2.948e-03`; model uncertainty `6.032e-03`; old two-grid mean bias `6.032e-03`.
- Grid convergence: `True` with observed delta-ratio `3.233e-01` and relative tail correction `4.134e-03`.

| grid | r_max | exterior dr | dw | E_defect | E_reference | chi_Q | fit RMS max | max residual |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 35x8 | 6 | 1.605556e-01 | 2.055556e-01 | 5.931640551817e+01 | 8.003232505849e+01 | 7.411555952525e-01 | 4.098e-13 | 1.087e-15 |
| 48x10 | 6 | 1.179592e-01 | 1.681818e-01 | 5.661117583814e+01 | 7.840082042374e+01 | 7.220737682612e-01 | 4.382e-13 | 1.766e-15 |
| 61x12 | 6 | 9.322581e-02 | 1.423077e-01 | 5.592800392480e+01 | 7.812216859185e+01 | 7.159043960620e-01 | 4.547e-13 | 2.442e-15 |

### Even-nw Characterization

- Canonical characterization JSON: `_scratch/pathA_22b_gate3_nw_characterization.json`.
- Outcome class: `NW_CONVERGES_chi_inf=0.711956789348_pm_0.00418_order=1.525`; reason: finite power-law limit with stable coarsest-point jackknife and second-nr confirmation.
- Tail-supported central uses the `nw>=56` fit: raw `7.120688981355322e-01`, reported as `0.712`.
- Numerical bar `0.0008` comes from `sqrt(max(tail RMS 4.838e-04, tail jackknife 7.404e-04)^2 + nr offset 3.177e-04^2)`.
- Conservative full-sweep fit uncertainty is kept separately: `0.0042` from the `nw=8..160` fit.
- Evidence: jackknife stable `True`, nr-independent `True`, stored flat tail through `nw=160`; review verification extends the flat tail to `nw=320`.

### Tiny-Defect Linearity

- Definition: branch-geometry reference plus a small fraction of the frozen condensate/gauge defect.
- Slope mean `-7.654165346717084e-02`; slope relative spread `8.774e-03`; linear toward zero `True`.

| defect fraction | chi_Q | chi_Q - 1 | (chi_Q - 1)/fraction | max residual |
| ---: | ---: | ---: | ---: | ---: |
| 1.000000e-03 | 9.999237316117e-01 | -7.626838828689e-05 | -7.626838828689e-02 | 4.443e-16 |
| 3.000000e-03 | 9.997707501453e-01 | -2.292498546721e-04 | -7.641661822403e-02 | 4.594e-16 |
| 1.000000e-02 | 9.992306004611e-01 | -7.693995389059e-04 | -7.693995389059e-02 | 4.549e-16 |

### Trivial Uniform Consistency

- The uniform self-ratio is demoted to a trivial consistency check: `uniform no-defect control coefficient divided by separately built uniform vacuum reference coefficient`.
- Uniform no-defect max `|chi_Q-1|`: `0.000e+00`; pass `True`. This is not evidence for the physical magnitude.

### Dimensional Check

- `defect omega^5 coefficient`: **CONSISTENT** (expected `T^5`, actual `T^5`). The defect coefficient multiplying omega^5 carries T^5.
- `vacuum omega^5 coefficient`: **CONSISTENT** (expected `T^5`, actual `T^5`). The same extraction on the uniform background carries T^5.
- `closure placement factor (R_exit/c_s)^5`: **CONSISTENT** (expected `T^5`, actual `T^5`). The outgoing-Hankel boundary-placement factor is common to defect and vacuum extractions.
- `radius-consistent chi_Q ratio`: **CONSISTENT** (expected `1`, actual `1`). chi_Q is the ratio of two same-radius omega^5 coefficients.

### Provenance

- research/pde/paper/pde.tex:1964.
- research/pde/paper/pde.tex:1985-1988.
- research/pde/paper/pde.tex:2034-2069.
- software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls:686-753.
- software/stage1_solver/mathematica/mt15_06_m1c_prep_crossengine.wls:896-985.

### Target-Blindness

- `TARGET_BLIND_PASS` over the new Gate-3 module, tests, and Mathematica cross-check. The final comparison constants are not used in the derivation.

### Residual Ledger

- Frozen finite-core background loaded; no nonlinear profile re-solve was performed.
- Linear source adapter uses the frozen derived BdG wall modes; the outgoing DtN/self-energy sweep is recomputed.
- Retarded operator uses exact spherical-Hankel Y_out at the truncation radius; cached Gamma_port is not used in the solve.
- Branch-reference r_max tail relative spread 1.879e-03; onset r_max=6.0; plateau True.
- Post-plateau fixed-r_max Richardson limit 7.129567649611814e-01; observed delta-ratio 3.233e-01.
- Matched even-nw Z_w reference shift 6.082e-02 (8.54% of the branch value); carried separately as a one-sided definitional systematic.
- Tiny-defect slope relative spread 8.774e-03; linear True.

### OUTCOME

- OUTCOME: `DELTA_Q_NE_0_0.712_pm_0.0008_NW_CONVERGES`.
- Converged `chi_Q = 0.712 +/- 0.0008` (numerical, tail-supported), branch-geometry reference.
- Conservative full-sweep fit uncertainty: `+/- 0.0042`; this is retained as a coarse-transient diagnostic, while the flat nw-tail plus nr-offset gives the honest numerical bar.
- Separate `Z_w`-reference definitional systematic: one-sided `0.061` absolute (`8.5%`) toward the flat `Z_w=1` alternative; not folded symmetrically into the numerical error bar.
- Convergence evidence: jackknife-stable `True`, nr-independent `True`, flat even-nw tail, and matched branch-geometry reference.
- Reason: The physically correct defect zero keeps the branch Z_w/R0_w geometry and removes only the condensate/gauge defect. The production magnitude is earned from the unfiltered even-nw branch-geometry characterization when that stored sweep has a stable finite limit, nr confirmation, and a matched flat-Zw reference comparison. The flat-Zw shift is a separate one-sided definitional systematic, not part of the numerical error bar. If the even-nw evidence is unavailable or non-convergent, the magnitude reverts to a calibration knob.

