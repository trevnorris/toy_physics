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

