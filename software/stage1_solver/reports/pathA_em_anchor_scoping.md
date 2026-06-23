# Path-A EM Anchor Scoping: What Can Pin `lambda_gamma = c_gamma/c_s`?

Date: 2026-06-22

Scope: scoping derivation/analysis only. No canonical decision or paper files are edited here.

## Executive Verdict

The canonical parent action admits many EM observables, but almost none isolate `lambda_gamma`.

The EM sector cleanly constructs:

- a zero-mode photon cone in flat bulk-metric units, with `C_B/C_E = 1`;
- a brane coupling reduction, `mu0_eff = mu0/Z_int`;
- a canonical charge reduction, `q_eff = q_*/sqrt(Z_int)`;
- a Coulomb plus Gaussian-KK Yukawa pattern;
- time-domain KK tails and massive-mode dispersion.

Those observables pin EM coupling strength, localization/profile data, or photon speed. They do not, by themselves, pin the ratio to the GNLS acoustic speed unless the observable is explicitly a speed comparison against `c_s`.

The honest minimal anchor for the xi verdict is therefore a speed anchor:

```text
lambda_gamma = c_gamma/c_s.
```

A fine-structure analog exists after restored-unit brane normalization, but it pins the product

```text
alpha_toy = mu0_eff q_eff^2 c_gamma/(4 pi hbar)
```

and therefore does not isolate `lambda_gamma` unless the EM coupling combination `mu0_eff q_eff^2` has already been independently fixed. The direct `lambda_gamma = 1` route closes the count, but it is not derived by the current action. It is either a calibration to the observed GW/light speed equality or an added metric-acoustic postulate.

## Source Facts Used

The parent action separates the GNLS medium from the localized Maxwell field. The matter sector gives

```text
c_s^2(rho) = (1/m) dP/drho = 5 K rho^4/m
```

from the stiff EOS (research/pde/paper/pde.tex:342-352). The Maxwell sector is an independent localized gauge action

```text
L_EM = - Z(w) F_MN F^MN/(4 mu0) - ... - A_M J_ext^M
```

with exact Maxwell equation `partial_M(Z F^{MN}) + ... = mu0 J_tot^N` (research/pde/paper/pde.tex:355-416). The charge ontology and zero-mode reductions are

```text
q_eff = q_*/sqrt(Z_int),    Z_int = int Z(w) dw
mu0_eff = mu0/Z_int
```

(research/pde/paper/pde.tex:289-295, research/pde/paper/pde.tex:553-563). The zero-mode reduction is explicitly a controlled reduction, not microscopic ontology (research/pde/paper/pde.tex:541-565).

Gate 2 computed the flat-metric zero-mode Maxwell cone:

```text
C_E = Z_int/mu0,    C_B = Z_int/mu0,    C_B/C_E = 1
```

so `Z_int` cancels from the cone, but the speed map remains unpinned (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:118-135). The beta addendum classifies this as `BETA_GENUINE_GAP`: the parent action does not identify the Maxwell bulk metric with the GNLS acoustic metric (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:180-191; software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:200-208).

The counting status is settled: after calibrating `g_G` on Newtonian `G`, the quadrupole leaves `g_mhat^2` and `lambda_gamma^5` underdetermined by one; an EM-sector anchor is load-bearing (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:107-115, software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:237-240).

## D1. Dimensionless EM Observables From The Canonical Action

Here "isolates `lambda_gamma`" means the observable alone fixes `c_gamma/c_s`, with `c_s` already derived from the GNLS medium. Observables that pin products involving `q_*`, `mu0`, or `Z_int` do not isolate it.

| Observable class | Dimensionless form | Source / derivation | Pins | Isolates `lambda_gamma`? |
| --- | --- | --- | --- | --- |
| Bulk zero-mode cone ratio | `C_B/C_E = 1` | Flat metric and zero-mode reduction give equal electric/magnetic principal coefficients (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:118-126). | Removes one cone-shape knob; pins only the isotropic Maxwell cone in bulk-metric units. | No. The bulk-to-brane/acoustic speed map remains `beta_bulk_to_brane` (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:128-136). |
| Photon/acoustic speed ratio | `lambda_gamma = c_gamma/c_s` | `c_s` is derived from EOS; `c_gamma` is the Maxwell cone speed after the unpinned speed map (research/pde/paper/pde.tex:342-352; software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:130-135). | Directly pins the missing xi factor. | Yes, if the observable is actually a speed comparison or an independently measured photon speed plus derived `c_s`. |
| Fine-structure analog | `alpha_toy = mu0_eff q_eff^2 c_gamma/(4 pi hbar)` | Restored-unit brane Maxwell coupling; see D2. | Pins `mu0_eff q_eff^2 c_gamma/hbar`, equivalently an EM coupling-speed product. | No, unless `mu0_eff q_eff^2` is already independently calibrated. |
| Coulomb strength made dimensionless with the photon quantum scale | `K_C/(hbar c_gamma) = alpha_toy`, where `K_C` is the Coulomb energy-length coefficient in the same brane normalization | Static zero mode gives Coulomb plus Yukawa; `A_0(r) = mu0_eff q_*/(4 pi r)[1+...]` in the unrecanonically normalized field (research/4d_em_fields/paper/4d_em_fields.tex:991-1008). Canonical normalization gives `q_eff` (research/4d_em_fields/paper/4d_em_fields.tex:822-856). | Same as `alpha_toy`: coupling-speed product. | No, for the same reason. |
| Coulomb strength made dimensionless with the acoustic quantum scale | `K_C/(hbar c_s) = alpha_toy * lambda_gamma` | Same Coulomb sector, but the denominator uses `c_s` instead of `c_gamma`. | Product `alpha_toy lambda_gamma`. | No. This is worse as a lambda anchor because `alpha_toy` and `lambda_gamma` are entangled. |
| Gaussian localization range | `a/lambda = sqrt(pi) a/Z_int` or `r/lambda = sqrt(pi) r/Z_int` | Gaussian `Z(w)=exp(-w^2/lambda^2)` gives `Z_int=lambda sqrt(pi)` (research/4d_em_fields/paper/4d_em_fields.tex:290-303). | Localization length/profile scale relative to an external length such as `a` or an observed separation. | No. It pins `Z_int/a` or `Z_int/r`, not `c_gamma/c_s`. |
| Yukawa amplitude ratios | `c_{2m}/c_0 = 4^{-m} binom(2m,m)`; leading correction coefficient `1/2` | Mode coupling formula and parity rule (research/4d_em_fields/paper/4d_em_fields.tex:938-972); leading potential correction (research/4d_em_fields/paper/4d_em_fields.tex:1017-1022). | Shape of the Gaussian KK tower. | No. Independent of `c_gamma` and `c_s`. |
| KK mass ratios | `m_n^2 = 2n/lambda^2`; `m_{2m}/m_2 = sqrt(m)` | Gaussian Sturm-Liouville spectrum (research/4d_em_fields/paper/4d_em_fields.tex:883-919). | Gaussian spectrum shape and localization length if one mass scale is measured. | No. Ratios are speed-free; absolute thresholds need `c_gamma/lambda`. |
| Massive-mode threshold products | `omega_n Z_int/c_gamma = sqrt(2 pi n)` for the Gaussian tower | Restoring speed in `omega_n = c_gamma m_n` and using `Z_int=lambda sqrt(pi)` from the spectrum above. | Product/ratio `c_gamma/Z_int`. | Not alone. With an independent `Z_int` range anchor and known `c_s`, it can infer `lambda_gamma`; otherwise it pins a two-constant combination. |
| Massive-mode group velocity | `v_g/c_gamma = k/sqrt(k^2+m_n^2)` | Above threshold, modes have `omega^2=k^2+m_n^2` and `v_g=k/omega<1` in units with light speed one (research/4d_em_fields/paper/4d_em_fields.tex:1415-1426). Restored units replace `1` by `c_gamma`. | Dispersion shape plus `c_gamma/lambda` if absolute units are used. | Only if the massless wavefront speed `c_gamma` is measured against `c_s`, or if `Z_int` is separately known. |
| Retarded tail shape | Bessel argument `m_n sqrt(c_gamma^2 Delta t^2-r^2)` and wavefront speed `c_gamma` | Massive retarded Green functions have support inside the light cone (research/4d_em_fields/paper/4d_em_fields.tex:1120-1141). | KK content and wavefront speed. | The wavefront speed can isolate `lambda_gamma`; the tail shape alone pins `c_gamma/Z_int`. |
| Charge sign | `eta_Q = sign(q_*) = +/-1` | Charge sign ontology (research/pde/paper/pde.tex:279-312). | Branch sign. | No. It carries no speed information. |
| Mixed-sector invariants | Dimensionless ratios can be formed from measured `E_w`, `C_a`, and brane fields after choosing scales | Mixed invariants are exact observables (research/pde/paper/pde.tex:473-490), and they survive in the linearized PDE (research/pde/paper/pde.tex:894-950). | Mixed-core response/profile data. | Not from the current action alone. A future branch law could make a speed observable, but no current closed form isolates `lambda_gamma`. |

Bottom line for D1: the only canonical observable that isolates `lambda_gamma` is a speed-cone observable, either directly `c_gamma/c_s` or indirectly a photon-speed measurement plus the already derived `c_s`. Charge, Coulomb, fine-structure, and KK-profile observables pin multi-constant EM combinations unless supplemented by separate anchors.

## D2. Fine-Structure Analog

A dimensionless EM coupling can be constructed, but the current source stack does not already use it as a declared observable.

### Restored-Unit Form

In a brane SI-like normalization, restore independent charge dimension `Q`:

```text
[q_eff]    = Q
[mu0_eff]  = M L Q^-2
[c_gamma] = L T^-1
[hbar]    = M L^2 T^-1
```

Then

```text
alpha_toy = mu0_eff q_eff^2 c_gamma/(4 pi hbar)
```

is dimensionless:

```text
[mu0_eff q_eff^2 c_gamma/hbar]
= (M L Q^-2)(Q^2)(L T^-1)/(M L^2 T^-1)
= 1.
```

I ran one dedicated dimensional check for this report. It returned:

```text
mu_eff*q_eff^2*c/hbar: 1
after q_eff=q_star/sqrt(Z), mu_eff=mu0/Z: 1
mu_eff*q_eff^2/(hbar*c): L^-2 T^2
```

So the tentative multiplicative-`c_gamma` form is the homogeneous restored-unit form; the inverse-`c_gamma` alternative is not homogeneous in the brane SI-like convention.

Substituting the reductions used by the sources,

```text
q_eff = q_*/sqrt(Z_int),    mu0_eff = mu0/Z_int
```

gives

```text
alpha_toy
= mu0 q_*^2 c_gamma/(4 pi hbar Z_int^2).
```

This expression assumes the same restored-unit convention in which the higher-dimensional `q_*` carries the corresponding `sqrt(Z_int)` charge normalization. The key physical point is invariant: `alpha_toy` is one dimensionless product involving EM charge normalization, EM kinetic normalization/localization, and photon speed.

### Normalization Caveat

The source formulas use two equivalent packages for the zero-mode coupling:

- unrecanonical Maxwell field: `mu0_eff = mu0/Z_int`, source charge `q_*`;
- canonically normalized brane zero mode: `q_eff = q_*/sqrt(Z_int)`, kinetic coefficient `mu0`.

The EM paper explicitly gives both packages (research/4d_em_fields/paper/4d_em_fields.tex:768-779, research/4d_em_fields/paper/4d_em_fields.tex:822-856). The measured Coulomb strength is the invariant product, not an independent measurement of both field normalizations. One must not spend both `mu0_eff` and `q_eff` as if they were independent reductions of the same field normalization.

### What `alpha_toy = alpha_physical` Pins

Equating to physical `alpha` pins

```text
mu0_eff q_eff^2 c_gamma/hbar = 4 pi alpha_physical,
```

or, after substitution,

```text
mu0 q_*^2 c_gamma/(Z_int^2 hbar) = 4 pi alpha_physical.
```

This does not isolate `lambda_gamma` alone. It is one equation for a product containing `q_*`, `mu0`, `Z_int`, and `c_gamma`. It can recover

```text
lambda_gamma
= (4 pi alpha_physical hbar)/(mu0_eff q_eff^2 c_s)
```

only if `mu0_eff q_eff^2` has already been independently fixed by an EM coupling/charge anchor. Without that extra EM calibration, `alpha_toy` is not the missing one-unknown closure for the xi verdict.

## D3. Crux: What Speed Does The Toy Quadrupole Radiation Use?

### D3(a). Propagation Speed In The Toy Quadrupole DtN

The toy outgoing quadrupole radiation used for `chi_Q` propagates on the GNLS medium/outgoing-port speed `c_s`, not on the Maxwell photon speed `c_gamma`.

The direct source is the compact outgoing `l=2` fingerprint:

```text
Y_out(omega)
= 1
  + a^2 omega^2/(9 c_s^2)
  + 4 a^4 omega^4/(81 c_s^4)
  + i a^5 omega^5/(27 c_s^5)
  + O(omega^6).
```

This is written explicitly with `c_s` (research/pde/paper/pde.tex:1954-1967). The one-pole retarded package uses

```text
Omega_Q = 3 c_s/(2a),
sigma_Q^can = 4 a^5/(27 c_s^5)
```

(research/pde/paper/pde.tex:1969-1988). The leading odd bridge is

```text
Gamma_5 = chi_Q P_0 a^5/(27 c_s^5)
```

(research/pde/paper/pde.tex:2040-2063). Gate 3 also extracted the `omega^5` coefficient from a spherical-Hankel outgoing DtN and recorded the closure placement factor as `(R_exit/c_s)^5` (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:193-204, software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:284-289).

The GR comparison target, however, is the Burke-Thorne coefficient

```text
mhat^2 Gamma_5 = 2 G/(5 c^5),
```

with `c` the relativistic light speed in the target theory (research/pde/paper/pde.tex:2073-2085). Path-A then carries `c = lambda_gamma c_s` in the dimensionless audit (software/stage1_solver/reports/pathA_22a_dimensional_skeleton.md:56-61), yielding the verdict

```text
P0 chi_Q g_mhat^2 lambda_gamma^5/g_G = 54/5
```

(software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:37-39).

So the speed split is:

- toy outgoing quadrupole DtN/radiation carrier: `c_s`;
- target GR light-speed denominator: `c = c_gamma`;
- verdict conversion: `lambda_gamma^5 = (c_gamma/c_s)^5`.

### D3(b). Is `lambda_gamma = 1` Forced?

Not by the canonical action.

The declared parent action treats EM as an independent localized gauge field on the flat bulk metric, minimally coupled to the GNLS medium. Gate 2 proves only `C_B/C_E = 1` in bulk-metric units; it does not identify that unit cone with the acoustic cone (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:118-136). The beta addendum explicitly states that the action supplies no relation tying `beta_bulk_to_brane^2` to `5 K rho0^4/m_GNLS` (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:182-190). Decision 14 therefore classifies `lambda_gamma` as a genuine input/gap, not a convention and not a derivation (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:54-59, software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:200-208).

If the toy outgoing quadrupole wave is identified with real gravitational radiation, then physical GR faithfulness wants its propagation speed to equal the light speed. That is the statement

```text
c_s = c_gamma,    lambda_gamma = 1.
```

But this is an extra physical identification relative to the current parent action. The verdict's `lambda_gamma^5` can absorb a mismatch algebraically; it does not prove the cones are equal.

### D3(c). Classification Of `lambda_gamma = 1`

Using the beta taxonomy:

- Not `CONVENTION`: the ratio `c_gamma/c_s` is dimensionless and invariant under common unit rescalings (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:189).
- Not `DERIVABLE`: the parent action lacks the metric-acoustic identification (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:182-191).
- Current status: `BETA_GENUINE_GAP`.
- If set from the measured physical fact "GW speed equals light speed": calibration of a genuine-gap input to an external measured fact.
- If inserted as an internal law "the Maxwell cone and acoustic cone are the same": added postulate, specifically a `bulk_metric_to_acoustic_metric_identification` postulate.

Those two choices have the same numerical consequence, `lambda_gamma=1`, but they change the model's character differently.

## D4. Consequences Of Adopting `lambda_gamma = 1`

### Count Closure

Yes. With `g_G` calibrated on Newtonian `G`, `P0` and `chi_Q` derived, and `lambda_gamma=1`, the verdict becomes

```text
P0 chi_Q g_mhat^2/g_G = 54/5.
```

Therefore

```text
g_mhat = sqrt((54/5) g_G/(P0 chi_Q)).
```

So the quadrupole pins `g_mhat` directly. The count closes, but the GR quadrupole remains an absorbed calibration anchor rather than a prediction, as Decision 14 already states for the post-Gate-4 status (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:237-240).

### Held-Out Surplus Under `lambda_gamma = 1`

The `lambda_gamma`-carrying tails become numerically trivial as speed-ratio discriminants, but the held-out observables do not all become trivial.

Decision 14 maps the held-out dependencies as follows:

- `g-2`: depends on `c_s`, `chi_Q`, and `Xi_1=P_1/P_0`; no explicit `lambda_gamma` entry is settled (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:127).
- `5PN`: carries a `lambda_gamma^3` tail and depends on `P0,P2,P4,D0`, `chi_Q`, and `g_mhat^2` (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:128).
- `ringdown/QNM`: depends on `c_s` and branch moments, with `lambda_gamma` still marked uncertain (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:129).
- `moving-throat`: carries a `lambda_gamma^3` tail plus the full bundle and explicit EM sector `mu0,q_*,Z_int` (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:130).

Therefore:

- The `5PN` and moving-throat `lambda_gamma^3` factors become `1`. The speed-ratio part of those tests is no longer discriminating.
- The `5PN` observable can still be sharp if the branch-derived `P2`, `P4`, `D0`, and shared `chi_Q` are computed target-blind.
- The moving-throat observable still tests the full branch bundle and EM-sector constants/profile data.
- `g-2` retains predictive content through `Xi_1` and `chi_Q`.
- Ringdown retains whatever branch/QNM content survives its still-aspirational status.

So `lambda_gamma=1` is not globally trivializing. It sharpens the model by removing one free factor, but it kills the specific discriminating power of any held-out whose only novel dependence was the `lambda_gamma^3` tail.

## D5. Honest EM Anchor Menu

| Option | Physical input | Combination pinned | Isolates `lambda_gamma`? | Cost | Consequence for held-outs |
| --- | --- | --- | --- | --- | --- |
| Fine-structure analog | Match `alpha_toy` to measured `alpha` | `mu0_eff q_eff^2 c_gamma/(4 pi hbar)` | No, unless `mu0_eff q_eff^2` is already fixed. | Uses a real EM measurement, but spends it on a multi-constant product. Also requires adopting the fine-structure construction, which is not currently a declared canonical observable. | Preserves nontrivial `lambda_gamma` tails only if additional anchors recover `c_gamma`; by itself leaves count open. |
| Direct photon speed / photon-acoustic comparison | Measure or calibrate the photon cone `c_gamma` and compare to derived `c_s` | `lambda_gamma = c_gamma/c_s` | Yes. | Calibration to a measured speed fact; no new internal postulate if treated as external calibration. | Preserves maximal held-out content if measured `lambda_gamma` is not forced to 1; tail factors remain nontrivial predictions using the fixed value. |
| `lambda_gamma=1` from GW speed equals light speed | Identify toy quadrupole wave speed with real GW speed and use empirical GW/light equality | `lambda_gamma=1` | Yes. | If external: calibration to measured fact. If internal: added metric-acoustic postulate. Not derived by the action. | Count closes. `lambda_gamma^3` tails become unity, so speed-ratio discrimination is lost; other branch/bundle predictions remain. |
| Photon dispersion / massive KK thresholds | Measure massless wavefront plus massive thresholds/group velocities | Massless front pins `c_gamma`; thresholds pin `c_gamma/Z_int`; static range can pin `Z_int` | Yes only if the massless front is compared to `c_s`, or if `Z_int` is separately anchored and `c_s` known. | Constructible from current localized Maxwell action. Practically, massive KK modes may be only bounded, not detected. | If it pins a nontrivial `lambda_gamma`, it preserves speed-tail discriminants. If only bounds exist, verdict count is not exactly closed. |
| Coulomb/charge strength | Coulomb law normalization | Invariant Coulomb strength, e.g. `mu0_eff q_*^2` or canonical `mu0 q_eff^2` depending on field normalization | No. | Real EM measurement, but dimensionful and normalization-convention sensitive. Needs `hbar` and a speed to become `alpha`. | Anchors EM coupling/profile sector, not the xi speed gap. |
| Static Yukawa deviations | Precision Coulomb-law departures | `Z_int/a` or profile length/range, plus fixed Gaussian amplitude ratios | No. | Constructible if `Z(w)` is genuinely decaying. Current exported `Z_int` is blocked by nondecaying floor/finiteness issues (software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md:69-85). | Useful for EM profile predictions and moving-throat EM sector; does not close xi count alone. |
| Bulk `C_B/C_E=1` | None; derived from action | Equal principal coefficients in flat metric units | No. | Free derived result, already done. | Removes a fake cone knob but leaves `beta_bulk_to_brane` open. |
| Charge sign / ontology | Branch sign assignment | `eta_Q = +/-1` | No. | Already part of canonical ontology (research/pde/paper/pde.tex:279-312). | No lambda consequence. |
| Legacy emergent-acoustic photon identification | Declare photon emergent on acoustic cone | `c_gamma=c_s`, hence `lambda_gamma=1` | Yes. | Added postulate/change of parent-action character. Decision 14 identifies this as the superseded route, not canonical Path-A evidence (software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md:205-208). | Same numerical consequence as `lambda_gamma=1`, but with stronger model-character cost. |
| Reverse-engineer `lambda_gamma` from `54/5` | GR quadrupole target | `g_mhat^2 lambda_gamma^5` product | No independent isolation. | Forbidden target-facing calibration; it consumes the very verdict being tested. | Destroys the independent-anchor discipline. |

## D6. Recommendation

Use a speed anchor, not a charge anchor.

The cleanest next directive should treat `lambda_gamma` as the directly calibrated photon/acoustic cone ratio:

```text
lambda_gamma := c_gamma/c_s.
```

This closes the verdict count with the fewest new theoretical assumptions because `c_s` is already derived from the GNLS EOS and the photon cone is the actual missing EM speed input. It is also the most honest toy-analog framing: the parent action has independent cones, so the model must either measure/calibrate their ratio or explicitly add a cone-identification postulate.

I do not recommend using `alpha_toy` as the minimal lambda anchor. It is a legitimate dimensionless EM coupling, but it pins `mu0_eff q_eff^2 c_gamma/hbar`, not `c_gamma/c_s` alone. It is useful for the broader EM-sector calibration, especially charge/coupling normalization, but it is not the one-equation closure for the current xi underdetermination unless paired with another EM coupling anchor.

The conceptual decision is a clean either/or:

1. **Analog-calibration route:** Calibrate `lambda_gamma` from the photon cone relative to the derived acoustic speed. This is a calibration to a measured speed ratio. It keeps the parent action's independent-cone character and preserves nontrivial `lambda_gamma^3` tail factors if the calibrated ratio is not exactly one.

2. **GR-faithful equality route:** Set `lambda_gamma=1` because the toy quadrupole wave is required to share the physical light cone. This closes the count fastest, but it is not derived. It must be labeled either as a measured-fact calibration or as an added metric-acoustic postulate. It makes the `lambda_gamma^3` held-out tails numerically trivial while leaving branch/bundle held-outs such as `g-2`, `P2/P4` 5PN structure, ringdown, and moving-throat EM/profile content as the remaining predictive surplus.

For the next gate directive, the most defensible wording is:

```text
Anchor lambda_gamma by an independent speed-cone input.
Default: calibrate c_gamma/c_s directly.
If the user elects GR-faithfulness, set lambda_gamma=1 and classify it explicitly as CALIBRATION_TO_MEASURED_FACT or ADDED_POSTULATE, not DERIVED.
```
