# S11c-c2 N6: native material face-slot construction (astra)

Authored 2026-09-06 against `be542cb2`. This is a construction specification for review and subsequent folding into §5c and the build directive. It does not supply an expected N6 residual or report a completed N6 computation.

All script references below are under `research/pde_ledger_v3/scripts/`. Abbreviations: **a** = `S11c_a_interface_geometry_sympy_audit.py`; **b** = `S11c_b_brane_operator_sympy_audit.py`; **c2** = `S11c_c2_selfenergy_fold_sympy_audit.py`. The checked source SHA-256 values are:

```text
a   f4da0f2224257381300e5e78cc74a1341485b82446ac98251ca25d4b3c1a6be7
b   2a1fef275636dd210c8864ab1762859f63a4dfb054008fa21620e54b14e47fd2
c2  24d3ecb40c999604d6dc529cd004d62a9db9cb4186d50896663d35d46499d2e3
```

**1. The object and the binding decision**

For each fixed pair

\[
\alpha\in\{\mathrm{LAB\_HELD},\mathrm{MATERIAL\_ADVECTED}\},\qquad
\rho\in\{\mathrm{RHO4\_CONSTANT},\mathrm{RHOBR\_CONSTANT}\},
\]

construct

\[
 I_E^{\alpha,\rho}=\mathsf{extract}\bigl(\mathsf{close}_E(S_{P,E})-S_{P,E}\bigr),
 \qquad
 I_{M\to E}^{\alpha,\rho}=\mathsf{extract}\bigl(\mathsf{close}_M(S_{P,M})-S_{P,M}\bigr),
 \qquad R_{N6}^{\alpha,\rho}=I_E^{\alpha,\rho}-I_{M\to E}^{\alpha,\rho}.
\]

Here `SLAB_M` means the row-valued face carrier \(S_{P,M}\), not a material bulk action or a fully rebuilt slab. `close_E` and `close_M` are the **same imported same-(α,ρ) c1 response**, with their respective source arguments. There is **no separate transformation T** on either carrier, closed increment, or extracted object. Both inputs to the single Eulerian `extract` already use Eulerian face covectors, row components, fields, test functions, and pressure-slot identities.

**Decision: route 2 binds material μ_θ**, obtained by varying the material-route bulk density. Route 1 binds the imported Eulerian `mu_theta_operator`. Bind the chosen amplitude consistently in the face substrate **and in the c1 response source**. Binding material μ only in the open face work is insufficient: its pressure-independent contribution cancels from the increment.

The reasoning has two parts. The parent `a.task_rep_invariance` (1481–1496) constructs native material geometry and differences it directly against `RAW_PRIMARY`; `build_material_face_source` performs the covector conversion internally. It does not authorize a transformation of an increment. However, the parent leaves the same anchoring-specific `mu_theta_branch[branch]` as a constitutive placeholder in both face routes; **that geometry test alone does not determine which bulk operator c2 should bind**. For the requested native material *constitutive-and-face* construction, instantiate that placeholder with the material-route constitutive derivative, using the existing `b.build_operator` density-routing convention (2762–2783). Otherwise one tests only geometry representation with an Eulerian constitutive source in both routes. Those are different controls. The choice here includes the N4 constitutive channel without adding a transformation to an already mapped face object. Equality of the resulting increments remains a finding to be adjudicated, not a premise used to simplify or repair them.

In particular, `mu_theta_M` versus `mu_theta_L` names **anchoring**, not the coordinate route. At `LAB_HELD`, material-route μ is bound to the placeholder `mu_theta_L`; at `MATERIAL_ADVECTED`, it is bound to `mu_theta_M`. Neither anchoring nor density changes between operands.

**2. Construct the material constitutive source before taking variations**

Let \(\mathcal E_\alpha\) denote `b.construct_energy(alpha, ...).density`, with the same energy coefficients and live background-jet convention used by the imported Eulerian slab. Define

\[
 g_i=\frac{D_i\rho_4}{\rho_4},\qquad a_\rho=u_i g_i,\qquad
 h_\alpha=\begin{cases}u_iD_iW_{bg}/W_{bg},&\alpha=\mathrm{LAB\_HELD},\\0,&\alpha=\mathrm{MATERIAL\_ADVECTED},\end{cases}
\]

\[
 \mathcal E_{M,\alpha,\rho}
 =\bigl[(1+\operatorname{tr}\nabla u)\,
 \mathcal E_\alpha(\theta+a_\rho,e_W+h_\alpha,
 D_i(\theta+a_\rho),D_i(e_W+h_\alpha),u,\nabla u)\bigr]_{\text{wave degree }2},
\]

\[
 \bar\mu_M^{\alpha,\rho}
 =\frac{\partial\mathcal E_M}{\partial\theta}
   -\sum_iD_i\frac{\partial\mathcal E_M}{\partial(\partial_i\theta)}.
\]

These are precisely the typed uses of `b.material_pullback` (1916–1981) and the second return value of `b.operator_from_density` (2290–2384). The latter is an **amplitude**; its exported physical row has an additional `epsilon_shape`. For route 1,

\[
 \bar\mu_E^{\alpha,\rho}
 =\texttt{inputs.mu[(alpha,rho)][1]}/\epsilon.
\]

Use the live `total_derivative` chain and the same coefficient provenance; do not use `committed_material_pullback`, `frozen_derivative`, a Hessian freeze, or a new background-density law. The slab convention has `STRONG_ROW_JET_DEPTH=3`; the additional weak spatial derivative must retain the corresponding depth-4 information. Build termwise derivative/grade circuits so this definition does not require expansion of the whole bulk operator. Independently PIT-check the unpulled derivative of the reconstructed energy against the imported μ amplitude. This check detects coefficient or grade drift and is not a prescription to replace the independently constructed material result with the imported result.

The density Jacobian belongs **only in this quadratic scalar construction**. For a homogeneous quadratic input, its linear correction produces cubic wave terms and is removed by the quadratic projection. Do not force a Jacobian contribution to survive. Do not apply `material_pullback` to linear rows or to pressure-slot coefficients.

The exact representative laws from `b.density_pair` are

\[
 \rho_4=\rho_{br}/W_0\quad\Rightarrow\quad g_i=0
 \quad (\mathrm{RHO4\_CONSTANT}),
\]

\[
 \rho_4=\rho_{br}/W_{bg}\quad\Rightarrow\quad
 g_i=-D_iW_{bg}/W_{bg}
 \quad (\mathrm{RHOBR\_CONSTANT}).
\]

Do not confuse this gradient of **four-density** with the separate `BACKGROUND_ADVECTION` origin in the projected mass row, which involves the integrated density and cancels under pressure-slot projection.

**3. SLAB_M: the native face fold, with exact pressure identities**

Define the four amplitude slots

\[
 P=(\texttt{delta_p_plus},\texttt{d_w_delta_p_plus},
     \texttt{delta_p_minus},\texttt{d_w_delta_p_minus}),
 \qquad \mathcal P_P F=F-F|_{P=0}.
\]

They are the **same SymPy symbols**, including assumptions, used by the imported route-1 slab and `c2.build_face` (509–515). `a` inherits the two pressures from S11b and declares the two normal jets (163–164); `b` binds those four inherited identities (480–492). A small direct SymPy check during authoring compared all four `a` symbols with `S11c_b_exports.LEDGER[name]['value']`: all four equalities were `True`. No pressure-symbol mismatch was found. The implementation must print a further identity/assumption census against `inputs.a(name)` at the c2 boundary; do not silently create new symbols by spelling.

Create a **new bundle adapter** with the schema accepted by `b.selected_substrate_axes` and `b.face_generalized_force_rows`. Its values come from the following material builders, filtered to the fixed α,ρ and `DELTA_W`. The dictionary builders specify the authoritative quantities; for performance, use their exact single-case source functions instead of looping over every case.

| Bundle key / dependency | Native material builder and selected axes |
|---|---|
| `traction` | `a.build_geometry_quantity('TRACTION', route='MATERIAL')`, `(α,s,'DELTA_W',ρ)` |
| `closure_shape_deriv` | `a.build_geometry_quantity('CLOSURE_SHAPE_DERIV', route='MATERIAL')`, `(α,s,'DELTA_W',ρ)` |
| `virtual_work_shape_deriv` | `a.build_geometry_quantity('VIRTUAL_WORK_SHAPE_DERIV', route='MATERIAL')`, `(α,'DELTA_W','DELTA_W',ρ)` and `(α,'DELTA_W','ZETA_C',ρ)` |
| `virtual_constraint` | `a.build_virtual_constraint_raw(route='MATERIAL')`, `(α,'DELTA_W',ρ)` |
| `evolution_mass_balance` | `a.build_evolution_raw('EVOLUTION_MASS_BALANCE', route='MATERIAL')`, `(α,'DELTA_W',ρ)` |
| `evolution_term_origins` | `a.build_evolution_raw('EVOLUTION_TERM_ORIGINS', route='MATERIAL')`, same axes |
| face velocity \(V_{M,s}\) | `a.build_geometry_quantity('FACE_VELOCITY', route='MATERIAL')`, `(α,s,'DELTA_W')` |

The single-case equivalents are `build_face_source(..., route='MATERIAL')` → `traction_raw`, `closure_raw`, `face_velocity_raw`, `true_area_flux_raw`; `virtual_work_cases(..., route='MATERIAL')`; `virtual_constraint_route`; and `evolution_route`. Apply the same profile finalization and retained grades as their dictionary builders. `FACE_VELOCITY` has no density axis in `build_geometry_quantity`; its default representative is `RHO4_CONSTANT`. This is a density-independent geometric velocity, not a direction to change the case's constitutive density.

Inside `build_material_face_source` (703–813), flattening is differentiated at fixed material X, its covector is mapped with `material_inverse_transpose` (742–755), the Eulerian graph conormal is normalized, and virtual displacement and face velocity are constructed from the material kinematics (772–794). `face_velocity_raw` returns \(\epsilon\,\partial_\lambda(\hat n\cdot v_s)|_0\). Reuse this result; a second inverse transpose, area Jacobian, θ substitution, or thickness remap on these outputs is forbidden.

Use the actual fold in `b.build_operator`, **only its face and pressure-slot portion**:

1. Bind \(\bar\mu_M\) into the material closure, traction, and virtual work with `b.bind_mu_theta_operand`. In `face_generalized_force_rows` (2150–2179), binding precedes `diff(..., delta_v_u[a])` and `diff(..., delta_v_e_W)`. Preserve that order and the separate center-work operand.
2. Let \(Q_{M,s}\) be the bound `closure_shape_deriv` for face s, \(B_M\) the material evolution mass balance, and \(O_M\) its origins. Form \(\widehat B_M=B_M-\sum_sQ_{M,s}\), and subtract the same sum from `O_M['TRUE_AREA_FACE_FLUX']`, exactly as `b.build_operator` (2814–2851). This replaces the supplied relative-flux contribution by its closure-law contribution **once**. Do not add another mass equation or another area factor.
3. Call `face_generalized_force_rows(material_bundle, α, ρ, corrected_origins, mu_M_amplitude)` or its termwise equivalent. Form the row dictionary
   \[
   F_M=\{U:F_{M,u},\ E_W:F_{M,e},\ \mathrm{THETA}:\widehat B_M\},
   \qquad S_{P,M}=\mathcal P_P F_M.
   \]
   Keep the traction/closure/work source operands as provenance; traction is already consumed through virtual work and must not be added again.
4. Carry the virtual constraint as native material provenance. Its constraint-reaction rows and all bulk/kinetic rows are pressure independent: do not build them merely to subtract them. Print their pressure-dependency census and the slot guard below. The center row is provenance, not a fourth c2 balance.
5. Route 1 is \(S_{P,E}=\mathcal P_P\texttt{c2.expanded_rows(inputs.slab[(α,ρ)])}\). Preserve U's three components plus the THETA and E_W rows and both faces. No imported `coupling_kernel` is an ingredient of either increment.

The adapter must distinguish raw dictionaries from a's exported `(axes, metadata)` rows and b's `(axes, value)` bundle rows; do not feed raw dictionaries to `selected_substrate_axes`. There is also an existing non-pressure jet naming distinction: `a.grad_theta` uses `theta_d1` etc., while b's constitutive operator uses `grad_theta_1` etc. These must be identified as derivatives of the same field at the jet/physical-field boundary, with an explicit table. This is not a mismatch of P, and it must not be handled by indiscriminate symbol renaming.

**4. What survives pressure projection: a checked limitation**

The full virtual work contains both μ and pressure, but additively. Its schematic traction amplitude is

\[
 p+\Lambda_X(\bar\mu/\rho_{br}^{bg}-p/\rho_m).
\]

Consequently the **open pressure coefficients themselves carry neither θ nor μ**. A direct small SymPy calculation of `finalize(virtual_work_cases(...))[1]`, followed by pressure projection and virtual differentiation, confirmed this for both anchorings and both coordinate routes at `RHOBR_CONSTANT`. It also found zero mechanical U pressure carrier for `LAB_HELD`, and an order-σ_W U pressure carrier for `MATERIAL_ADVECTED`. These are local face calculations, not an evaluation of R_N6.

For orientation, let \(\Pi_s=p_s+s\eta W_0w_1\,p_{s,w}/2\), \(\beta=1-\Lambda_X/\rho_m\), and use the retained profile convention \(D_iW_{bg}=\sigma_W w_{1,i}\). The checked mechanical carriers reduce to

\[
 (S_P)_{E_W}=-\epsilon\frac{W_0\beta}{2}(\Pi_++\Pi_-),\qquad
 (S_P)_{U_i}=\begin{cases}
 0,&\alpha=\mathrm{LAB\_HELD},\\
 \epsilon\frac{\beta D_iW_{bg}}2(\Pi_++\Pi_-),&\alpha=\mathrm{MATERIAL\_ADVECTED}.
 \end{cases}
\]

Use these only as a compact description of the check, **not as a shared hand-coded replacement for the two native constructions**. The THETA pressure row must still be obtained from the corrected material evolution/closure fold. N4 enters the self-energy through the **closed pressure source μ**, even though its direct μ terms in the open face rows cancel. A θ shift of \(S_P\) cannot reproduce that channel.

**5. close: identical opaque c1 operator, different native source arguments**

For each route r and face s, use `inputs.response[(α,s,ρ)]`, i.e. the imported `s11c_c1_face_response`. Resolve the live integrated density from `background_density_map` before composing `DELTA_P`; this is `c2.build_face` stages 0–2 (476–500). The source amplitude is, equivalently to the imported expression,

\[
 b_{r,s}=(1+\Lambda_V/\rho_m)\bar V_{r,s}
       +\frac{\Lambda_A}{\rho_m\rho_{br}^{bg}}\bar\mu_r,
 \qquad \bar V_{r,s}=V_{r,s}/\epsilon,
\]

\[
 p_{r,s}=\mathcal R_{\alpha,s,\rho}Z_{\alpha,s}[b_{r,s}],
 \qquad p_{r,s,w}=\text{the same outgoing normal continuation applied to }p_{r,s}.
\]

These formulas identify the slots of the **imported** expression; they do not authorize rederiving Z or its inverse. Preserve its noncommutative order and all retained scattering terms. Use the same `kernel_bridge` definitions and the same resolvent in both routes. `build_face` differentiates the outgoing factor \(\exp(i s q_{out}(w-sW_0/2))\) before setting the reference face, so the normal-jet multiplier uses **s and the output leg**, not the input momentum or an arbitrary sign.

Bind \(\bar\mu_M\) and material \(\bar V_{M,s}\) directly to the c1 source slots. `build_face`'s `mu_override` and `velocity_override` callbacks can replace those amplitudes, but the stock function also expands `composed_source` before applying the kernel. A finite-field/retained-grade low-level evaluator is needed to obey the early-evaluation mandate; calling this stock path with giant symbolic replacements is not a tractability solution. The response is opaque with respect to coordinate mapping: do not substitute θ or multiply a Jacobian into `DELTA_P`, Z, or the resolvent.

Write the carrier as \(S_{P,r}=\sum_{p\in P}C_{r,p}p\) after its slot-linearity check. Then compute

\[
 J_r=\Bigl[\sum_{p\in P}C_{r,p}\bigl(\chi_{r,p}-p\bigr)\Bigr]_{\rm retained},
 \qquad I_r=\mathsf{extract}(\mathsf{physical\_fields}(J_r)).
\]

Here \(\chi_{r,p}\) is the appropriate pressure or normal-jet closure amplitude. This retains the bare subtraction even if later sector restriction kills some bare terms. Amplitude normalization supplies ε exactly once, through the original carrier. Apply the same weak restriction and grade truncation to both operands; never evaluate a Fourier/DtN integral and never set \(k_{in}=k_{out}\).

**6. Controls and the reverse-u issue**

**Advection control.** Introduce a dimensionless source tag t in the material **bulk density substitution** only:

\[
 \theta\mapsto\theta+t a_\rho,\qquad
 \partial_i\theta\mapsto D_i(\theta+t a_\rho).
\]

Recompute \(\mathcal E_M(t)\) and \(\bar\mu_M(t)\) before binding. Baseline is t=1; the corrupted material construction uses t=0. Keep the LAB_HELD e_W shift, density law, coefficients, material geometry, material V, c1 response, Z and resolvent unchanged. Bind the recomputed μ at **both** face and c1-source sites. Do not omit `BACKGROUND_ADVECTION` from the pressure-independent mass base and call that this control; do not omit a term after constructing the increment.

For `RHOBR_CONSTANT`, \(a_\rho=-u\cdot\nabla W_{bg}/W_{bg}\) is live. Emit the tagged μ contribution, pressure-source contribution, carrier contraction, and all six extracted-block sensitivities, with grades. A simple scalar check \(\mathcal E=B(\theta+g u)^2/2+C(\theta+g u)(e+h u)\) gives \(\bar\mu(1)-\bar\mu(0)=Bgu\): it demonstrates a pressure-source path at first gradient order; it is not a computation of the full N6 residual or a guarantee that every extracted block responds.

For `RHO4_CONSTANT`, evaluate \(D_i(\rho_{br}/W_0)\), \(a_\rho\), and their tagged dependency in the material μ circuit, and print **computed structural absence**. Do not generate an identical corrupted copy and emit A−A as an independence test. Material and Eulerian μ may still differ through the LAB_HELD thickness shift; absence of this θ-advection probe says nothing about the complete N6 residual. Had Eulerian μ been chosen for both routes, this N4 source channel would instead be absent in **both** density representatives; that is not this specification.

**Tilt control.** Supply an injectable Eulerian face-carrier factory, using the same fold in §3 with `route='EULERIAN'`. First PIT-compare its unmodified output with the imported \(S_{P,E}\), row by row, face by face, pressure and normal-jet coefficient separately. Print both operands and the reconstruction residual. Then mutate the **upper face normal's selected first background-slope contribution in direction 1**, inside the Eulerian face source before traction, virtual work, relative flux and closure are formed. One precise implementation introduces a sign tag on the upper background graph-conormal contribution \(-D_1W_{bg}/2\); take tag +1 for baseline and −1 for corruption, normalize the modified conormal, and propagate that normal to the face fold. Hold the virtual-displacement kinematics fixed for this deliberate normal error. Other face and direction contributions remain available and generic.

Use fresh source objects and mutation-aware caches. Do not edit a shared symbol/profile globally. For a carrier-only tilt probe, keep the c1-source μ and V at their baseline Eulerian values as well as keeping Z/resolvent fixed; the mutation then tests the normal-derived carrier coefficients directly. This breaks, among other things, the LAB_HELD cancellation between normal tilt and virtual vertical displacement. It is a deliberately inconsistent normal used as a diagnostic, not a new physical anchoring.

A small check using `a.build_face_source('LAB_HELD',1,'DELTA_W','RHOBR_CONSTANT')` extracted the coefficient of `delta_v_u_1 * W_bg_d1` in the background-normal work contribution `-probe_pressure * n.dot(virtual_displacement)`: it was 0 before the normal-only reversal and `-probe_pressure` afterward. This is a source-level witness for the injection, not a PIT result for the complete closed increment.

The existing `a.reverse_upper_x1` hook changes `source_jet_scales` throughout the Eulerian source, including virtual kinematics, so it is broader than this normal-only injection. The new factory must expose the actual normal operand being changed. The c2 `FLIP_FACE_SLOPE` override (384–387) modifies a c1 DtN Fourier jet and **does not implement this carrier control**.

**Reverse-u qualification: do not manufacture an off-diagonal.** `c2.extract` (337–341) already takes the curl of the U rows restricted to THETA, E_W, and longitudinal trials. Its existence or nonzero value is not evidence for N4. Genuine route provenance is the material inverse-transpose geometry construction and the material-density-derived μ passed into the pressure response. Genuine N4 sensitivity is a computed dependence on the tag t, traced from that source to a particular block, not the presence of a curl instruction.

Variation of the material **bulk scalar** does contain the reverse chain-rule contribution \(g_i\bar\mu\) in its U Euler–Lagrange row (plus the appropriate derivative terms). In the small quadratic check above, the mixed derivative is \(Bg+Ch\), versus \(Ch\) when θ advection is omitted. But that bulk U row has no pressure slots and cancels in this increment. Adding \(g_i\) times a face mass row to \(S_P\), or mapping a virtual θ slot again, would invent a different fold.

There is an additional retained-grade limitation: the N4-only change in μ from a quadratic density is linear in u and its jets, of at least first background-gradient order. The uncorrupted LAB_HELD U carrier is zero; the MATERIAL_ADVECTED U carrier already has one σ_W. Thus their N4-only U-carrier products are absent or at least σ_W² and are outside the retained rectangle. Scalar-trial reverse blocks cannot be credited with an N4 signal merely because they were already present without it. Print tag sensitivities before and after truncation and permit computed absence in reverse blocks. **An unconditional requirement that a retained N4-induced reverse U contribution survive in this pressure-slot increment is incompatible with this existing face fold.** Replace that old requirement with provenance and blockwise sensitivity reporting; do not satisfy it vacuously or add uncancelled bulk rows. The live RHOBR_CONSTANT advection test is at the μ/pressure source and its surviving forward couplings; its numerical responses remain to be computed.

**Independence evidence.** The main Eulerian operand remains imported; the main material carrier is rebuilt from native material sources and material constitutive variation. Neither is constructed by subtraction from, substitution into, or copying the other. For the tilt probe print the unmodified/corrupted Eulerian factory increments and their difference, with the material operand evaluated unchanged on the same samples. For the live advection probe print unmodified/corrupted material increments and their difference, with Eulerian unchanged. Record source fingerprints and evaluate the unaffected operand too. A nonzero change from a one-sided source corruption, alongside unchanged opposite-route data, supplies independence evidence. If no change is found, report that result and the probe's support; never repair the output to make a control bite. A failed Eulerian factory/import check must be reported as reconstruction drift, so the tilt difference is not misattributed.

**7. Tractability and exact finite-field PIT contract**

No `build_case` end-to-end for either route. No full `b.build_operator` rebuild merely to obtain pressure rows. No full-symbolic expansion or zero test over all retained grades of the closed operands. The F/G `carrier()` algebra (99–112) is reusable; its `pinned()` path depending on a completed `build_case` is not.

Build row/face/slot coefficient circuits first. Keep ε, η, σ_W as formal grade indices; evaluate generic scalar coefficients and truncate the grade algebra **before closure expansion**. Retain η⁰/η¹ × σ_W⁰/σ_W¹, including ησ_W and all corresponding allowed scattering terms. Differentiate through the live jet circuits before specializing their inputs, or use derivative-aware finite-field jets; substituting numbers and then taking a spatial derivative would incorrectly freeze the background. Preserve mixed-partial identifications, density relations, kernel relations, and the input/output distinction. Do not replace η or σ_W by random numbers before extracting grades.

Print these slot guards using the same PIT:

\[
 G_{lin,r}=S_r-S_r|_{P=0}
       -\sum_{p\in P}p\,\left.\partial_pS_r\right|_{P=0},
 \qquad G_{cross,r,pq}=\partial_p\partial_q S_r,
\]

with r=E using the imported expanded-row circuit and r=M using the face-fold circuit. Also evaluate the full circuit twice, at P and χ(P), numerically in the retained algebra, and compare their difference with the coefficient-carrier evaluation. This tests cancellation of the pressure-independent base without materializing a second symbolic closed slab. Keep both face families generic, test their separate contributions and their sum, and include denominators' pressure-dependency checks. Guards are emitted measurements, not assertions.

For residuals, guards, reconstruction checks and controls:

1. Split by fixed α,ρ; row/face as applicable; retained (ε,η,σ_W) grade; six weak blocks; alpha-normalized formal-integral signature; and non-integral remainder. Normalize only the allowed compact-support interior in-plane IBP. Keep Fourier phases, integration measures and support/distribution signatures formal. Test the rational coefficient circuits of this normal form; do not integrate or treat exponentials as floating-point functions.
2. Use several exact prime fields, with multiple independent generic sample draws for each prime and cell. For example, validate and record primes 1000000007, 1000000009 and 998244353, and use at least eight valid draws per prime per cell, increasing the count if the reported degree bound requires it. Represent i in an appropriate extension field when necessary. Sample both momentum legs independently; keep each corresponding coordinate, field jet, coefficient, normal momentum and response definition shared across routes and controls. Do not independently randomize an inverse: evaluate it from the same Z and denominator definition in both routes.
3. Clear rational denominators within each block/cell, preserving their nonzero conditions. Reject singular samples **jointly** across the compared operands. Retained inverse series require an invertible grade-zero denominator. Select real `Piecewise`/outgoing-branch cells before modular reduction; finite fields have no order. Record the cell and use a certified rational parametrization or algebraic-extension evaluation respecting the dispersion/normal-momentum relations. Do not sample q and k as unrelated symbols in an on-shell coefficient. Cover each nonsingular branch cell and report excluded transition/singular loci separately. An unsupported branch is an explicit coverage gap, not an all-zero result.
4. Derive a conservative numerator degree bound from each arithmetic circuit after denominator clearing (and after any cell parametrization). Record the sample sets, rejected counts, seeds, primes and bounds. For independent coordinate draws from sets of size N, with rejection probability bounded by E/N (for example using the sum E of the excluded-denominator degree bounds), a usable per-valid-draw upper bound is min(1,D/(N−E)), provided E<N. Multiply across independent valid draws for a fixed nonzero polynomial over a certified good prime, and union-bound the declared family of tests. If a different sampler is used, give its actual bound. Handle bad primes explicitly: justify nonvanishing reduction by content/height bounds or state the FN bound as conditional on at least one chosen prime being good. In the latter case, when the good prime is unidentified, use the worst per-prime bound, not a product assuming every prime is good. Multiple primes alone do not eliminate that condition.
5. A certified nonzero modular numerator is evidence of a nonzero coefficient in the declared formal weak-kernel normal form. State that scope; do not claim to have evaluated an integral or proved more equivalences than the allowed normalizer supplies. All-zero samples mean only “no nonzero found; conditional false-negative probability ≤ δ,” with conditions and coverage printed. Missing degree/cell certification means no certified δ for that part. Neither result is an exit predicate on R_N6.

The branch sampler, formal-kernel normalization and degree certificates are implementation work, not features already supplied by `retained_shape` or the current F/G diagnostic. Keep the transcript compact: counts, numeric residues, fingerprints, degree/cell metadata and small support tables; no giant expressions.

**8. Emit contract and concrete gaps to carry into the build directive**

PRINT the two operands **before** their residual, with full α,ρ labels:

```text
S11CC2_REP_INVARIANCE_EULERIAN_OPERAND
S11CC2_REP_INVARIANCE_MATERIAL_OPERAND
S11CC2_REP_INVARIANCE_RESIDUAL
```

Use equivalent operand-first triples for reconstruction, slot guards and live controls, and an explicit computed-absence record for the inactive advection representative. Every mathematical object, including intermediate μ, V, slot coefficients, pressures, normal jets, row components, weak blocks, guards and zero objects, carries `[L,T,M]` and its (ε,η,σ_W) support. Metadata-only tags are dimensionless. For zero objects retain the declared target units separately from empty computed support.

Use the inherited dimensional schema and derivative/product/integral rules, but compute consistency independently of the desired label. In particular, μ has `[-1,-2,1]`, V `[1,-1,0]`, pressure `[-2,-2,1]`, its normal jet `[-3,-2,1]`; U balance has `[-2,-2,1]`, THETA mass balance `[-3,-1,1]`, E_W balance `[-1,-2,1]`. Slot coefficients have row dimension minus slot dimension; the weak blocks include their test-field and derivative dimensions and must be checked separately. No stamping all blocks with an energy dimension. Require an able-to-fail dimensional diagnostic: for example, separately feed a deliberately extra W_0 into one pressure summand and print the resulting incompatible addition dimensions. Numeric specialization must not erase the dimensional circuit. Report unknown units and inconsistencies instead of coercing them to the target. No assert, no residual-zero builder exit, no expected R_N6 value.

Findings about present code support:

- The parent has all required material face, work, velocity, virtual-constraint and evolution builders. There is **no ready-made material c2 carrier bundle/factory**. `b.build_operator(route='MATERIAL')` changes the bulk density but still calls `substrate_bundle`, which reads the exported Eulerian face substrate; it is not SLAB_M. A thin raw-case/bundle adapter and face-only fold are required.
- No mismatch was found for the four pressure symbols in the direct a-to-b export comparison. The c2 boundary census and explicit non-pressure jet bridge remain required implementation checks.
- The μ route and source binding must be added explicitly. The published b μ is Eulerian (`task_slab_operator` calls `build_operator` with the default route). The existing low-level c2 callbacks can bind replacement amplitudes, but do not provide native material construction or early finite-field evaluation.
- `build_face_shift_raw` has no route parameter. It is not needed for this carrier: the material source already calls `physical_trace_fields`, including both pressure normal jets. If a future diagnostic requires a separately exported material `face_shift` table, that table needs an adapter; do not substitute an Eulerian table under a material label.
- The normal-only Eulerian carrier injection and its imported-carrier PIT check are missing. The existing broader a slope hook and c2 DtN-jet override do not supply that exact control.
- The old mandatory retained N4 reverse-u survival clause is unsupported for the pressure-slot increment, for the cancellation/grade reasons in §6. Preserve the authentic face fold and report channel support; do not add a new reverse row.
- The companion needs a finite-field circuit evaluator with branch handling, denominator/degree certificates, slot guards and dimensional failure evidence. Existing whole-symbolic close/extract helpers are not that evaluator.

The obsolete cross-anchoring c2 `REP_INVARIANCE_*` loop and the current T-on-increment language must be replaced when the reviewed spec is folded into the build directive. This authoring change edits only this specification. **WITHHOLD: R_N6 has no supplied expected value. Its computed difference and controls are adjudicated externally.**
