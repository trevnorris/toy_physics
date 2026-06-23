# pathA_23 Stage 0 action and contracts

## VERDICT

`ACTION_SPECIFIED_CLASSIFIED` - a candidate brane-elastic extension is specified and classified as `NEW_PARENT_ACTION`. The in-plane constitutive law is not derived from microstructure in Stage 0; it is carried as a postulated menu, so any later positive Maxwell result would be conditional. The brane-to-bulk coupling is `SYMMETRY_ALLOWED_POSTULATED`.

Stage 0 stops here. No no-leak result, constitutive winner, spectrum result, Maxwell dictionary, charge normalization, or cone-lock statement is claimed.

Failure-token audit:

| token | Stage-0 status |
| --- | --- |
| `FAIL_ILLDEFINED_BRANE_EMBEDDING` | not triggered; embedding, induced metric, normal, curvature assumptions, thickness profile, and BC contract are specified below |
| `FAIL_NO_CONSISTENT_COUPLING` | not triggered; a normal-work term and conserved-current source term are specified with units and provenance |
| `FAIL_UNSPECIFIED_SUBSTRUCTURE` | not triggered; the microstructure path is explicitly `POSTULATED` |
| `FAIL_MICROSTRUCTURE_REVERSE_ENGINEERED` | not triggered; no gyrostat/Cosserat substructure is claimed as derived |
| `FAIL_NO_ADMISSIBLE_CURRENT` | not triggered; the conserved-current plus \(\Phi_{\rm br}\) completion passes the cheap pre-check |
| `FAIL_PREFERRED_W_ANISOTROPY_AT_PRINCIPAL_LEVEL` | not triggered at Stage 0; declared brane principal blocks contain no \(k_w\), while later preferred-frame tests remain open |
| `FAIL_PROJECTION_REDUCTION_FIREWALL` | not triggered; \(Z(w)\), \(W(w)\), and \(B_\ell(w)\) are kept distinct |

## Explicit action

The Stage-0 candidate is a family
\[
S_{\rm cand}^{(X)} =
S_{\psi}^{\rm keep}
+ S_{\rm brane}^{(X)}
+ S_{\rm cpl},
\qquad
X\in\{\mathrm{Cauchy},\mathrm{rot},\mathrm{Cosserat}\}.
\]

Reconciliation with the frozen parent: `pde.tex` still contains the localized \(A_M\) Maxwell sector. This report does not edit or reinterpret that paper. The brane-elastic fork is a new parent action, not a variational rewrite of the \(A_M\) action. If it were later integrated side-by-side with the frozen \(A_M\) sector, that would introduce an explicit double-counting problem unless a later gate removes, identifies, or demotes one EM carrier.

### 1. Bulk scalar sector kept from `pde.tex`

\[
S_{\psi}^{\rm keep}
=\int dt\,d^4X\,\mathcal L_\psi,
\]
\[
\mathcal L_\psi
=
\frac{i\hbar}{2}
\left(\psi^*D_t\psi-\psi D_t\psi^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
-V_{\rm conf}({\bf X};\Sigma)\rho
-U(\rho),
\]
\[
D_t\psi=\partial_t\psi+\frac{i q_\star}{\hbar}A_0\psi,
\qquad
D_i\psi=\partial_i\psi-\frac{i q_\star}{\hbar}A_i\psi,
\]
\[
\rho=|\psi|^2,\qquad
U(\rho)=\frac{K}{4}\rho^5,\qquad
h(\rho)=\frac{5K}{4}\rho^4.
\]

Units: \([\mathcal L_\psi]=E/L^4=M L^{-2}T^{-2}\), so
\([dt\,d^4X\,\mathcal L_\psi]=E T\). This is the frozen scalar sector from `research/pde/paper/pde.tex:316-348`; the flat metric and index conventions are from `pde.tex:230-255`.

### 2. Brane-elastic sector

Use material brane coordinates \(\xi^a\), \(a=1,2,3\), and time \(t\). The displacement fields are
\[
u^a(t,\xi)\quad\text{in-plane physical displacement},
\qquad
u_w(t,\xi)\quad\text{off-brane displacement}.
\]

For any constitutive candidate \(X\),
\[
S_{\rm brane}^{(X)}
=
\int dt\,d^3\xi\,\sqrt{\gamma}\,
\left[
\frac{1}{2}\rho_\parallel\,\gamma_{ab}\dot u^a\dot u^b
+\frac{1}{2}\rho_w\,\dot u_w^2
-\mathcal U_\parallel^{(X)}
-\mathcal U_w
\right].
\]

The off-brane sector carried into later gates is
\[
\mathcal U_w
=
\frac{1}{2}\tau_w\,D_a u_w D^a u_w
+\frac{1}{2}\kappa_w\,(D^aD_a u_w)^2
+\frac{1}{2}k_w u_w^2.
\]

Units:

| quantity | dimension |
| --- | --- |
| \(u^a,u_w,\xi^a,w\) | \(L\) |
| \(\rho_\parallel,\rho_w\) | \(M/L^3\) |
| \(\mathcal U_\parallel,\mathcal U_w\) | \(E/L^3=M L^{-1}T^{-2}\) |
| \(\tau_w,\mu_{\rm br},\mu_R\) | \(E/L^3\) |
| \(\kappa_w\) | \(E/L\) |
| \(k_w\) | \(E/L^5\) |

The coefficient \(k_w\) is not assumed positive by Stage 0. If \(u_w\) is massless and unsuppressed, Stage 4 has the named failure path `FAIL_BENDING_MASSLESS_FIFTH_FORCE`.

### 3. Candidate in-plane constitutive menu

Stage 0 does not choose the winner.

**Cauchy/Navier candidate**
\[
e_{ab}=\frac{1}{2}(D_a u_b+D_b u_a),
\]
\[
\mathcal U_\parallel^{\rm Cauchy}
=
\frac{1}{2}\lambda_{\rm br}(D_a u^a)^2
+\mu_{\rm br}e_{ab}e^{ab}.
\]

**Rotational/MacCullagh candidate**
\[
r_a=\epsilon_{abc}D^b u^c,
\qquad
\mathcal U_\parallel^{\rm rot}
=
\frac{1}{2}\mu_R r_a r^a.
\]

**Cosserat/micropolar candidate**

This candidate adds an internal micro-rotation \(\varpi_a\). If Stage 2 refuses extra internal variables, this option is not available.
\[
\mathcal U_\parallel^{\rm Cosserat}
=
\frac{1}{2}\lambda_c(D_a u^a)^2
+\mu_c e_{\langle ab\rangle}e^{\langle ab\rangle}
+\frac{1}{2}\kappa_c(2\varpi_a-r_a)(2\varpi^a-r^a)
+\frac{1}{2}A_c D_a\varpi_bD^a\varpi^b+\cdots .
\]
Here \(\varpi_a\) is dimensionless, \(\kappa_c\) has units \(E/L^3\), and \(A_c\) has units \(E/L\).

### 4. Brane-to-bulk coupling

The coupling is
\[
S_{\rm cpl}
=
\int dt\,d^3\xi\,\sqrt{\gamma}\,
\left[
u_w\,\Pi_n[\psi,\Sigma]
+J_{\rm br}^a A_a^{\rm br}
-J_{\rm br}^0 \Phi_{\rm br}
\right]
+S_{\rm mouth}.
\]

Definitions:
\[
A_a^{\rm br}=\alpha_u u_a.
\]
\(\Phi_{\rm br}\) is an auxiliary scalar-potential analog included only to make a generic conserved current structurally admissible. It is not a Stage-0 proof of Maxwell gauge structure; the C5 obstruction remains a Stage-5 test.

\(\Pi_n\) is the normal generalized force density supplied by the bulk scalar/throat sector, with
\[
[\Pi_n]=E/L^4,
\qquad
[u_w\Pi_n]=E/L^3.
\]
Operationally it may be represented by the normal projection of the bulk stress or by the variation of the confinement geometry, but Stage 0 does not derive its coefficient.

The source-current units are
\[
[J_{\rm br}^0]=Q/L^3,\qquad
[J_{\rm br}^a]=Q/(L^2T),
\]
\[
[A_a^{\rm br}]=M L/(TQ),
\qquad
[\Phi_{\rm br}]=M L^2/(T^2Q),
\qquad
[\alpha_u]=M/(TQ).
\]
Thus \(J_{\rm br}^aA_a^{\rm br}\) and \(J_{\rm br}^0\Phi_{\rm br}\) both have units \(E/L^3\).

Provenance label: `SYMMETRY_ALLOWED_POSTULATED`. It is scalar under brane spatial rotations/diffeomorphisms in the chosen static gauge and uses the defect/throat current allowed by the charge ontology, but it is not derived from `pde.tex`.

### 5. Finite-thickness representation

The brane is treated as a regulated finite-thickness layer with normalized material profile
\[
\int_{-\infty}^{+\infty}B_\ell(w)\,dw=1,
\qquad
[B_\ell]=L^{-1}.
\]
The brane integral can be written as a bulk integral with \(B_\ell(w)\):
\[
\int dt\,d^3\xi\,\sqrt{\gamma}\,\mathcal L_{\rm brane}
\quad\leftrightarrow\quad
\int dt\,d^4X\,B_\ell(w)\,\mathcal L_{\rm brane}
\]
in the flat, small-slope chart. Coefficients distribute as
\[
\rho_{\parallel}^{(4)}(w)=\rho_\parallel B_\ell(w),
\qquad
\mu_{\rm br}^{(4)}(w)=\mu_{\rm br}B_\ell(w),
\qquad
\mu_R^{(4)}(w)=\mu_R B_\ell(w),
\]
and similarly for \(\rho_w,\tau_w,\kappa_w,k_w\). The thin-brane limit is \(B_\ell(w)\to\delta(w)\), but Stage 0 keeps \(B_\ell\) regulated to avoid hidden \(w\)-Jacobian and delta-squared errors.

## Embedding and symmetry contract

**Embedding.** The reference brane is \(w=0\) in the flat bulk chart
\[
X^i=(x,y,z,w),\qquad \eta_{MN}=\mathrm{diag}(-1,+1,+1,+1,+1)
\]
using the frozen parent convention. A perturbed brane material point is
\[
Y^i(t,\xi)=\bigl(\xi^a+u^a(t,\xi),\,u_w(t,\xi)\bigr).
\]

**Induced metric and normal.**
\[
\gamma_{ab}=\delta_{ij}\partial_aY^i\partial_bY^j.
\]
At linear order,
\[
\gamma_{ab}=\delta_{ab}+\partial_a u_b+\partial_bu_a+O(u^2,\partial u_w\partial u_w).
\]
The unit normal is
\[
n_i=\frac{(-D_a u_w,\,1)}{\sqrt{1+D_cu_wD^cu_w}}+O(u_\parallel D u_w),
\]
and the small-slope extrinsic curvature is
\[
K_{ab}=D_aD_bu_w+O((Du_w)^2,u_\parallel).
\]

**Extrinsic-curvature assumptions.** The in-plane principal operator is evaluated on the flat reference brane. Extrinsic curvature enters only through the \(u_w\) bending sector and through throat/mouth boundary data. Background curvature is zero away from defects.

**Physical status of \(u_\parallel\).** \(u^a\) is a physical in-plane displacement of brane material, not merely a material relabeling. A later rotational candidate may make \(u^a\to u^a+\alpha_u^{-1}D^a\chi\) a gauge-type redundancy only if the full kinetic, source, and constraint system supports it. Stage 0 does not assume that.

**Bulk flow intersection.** The bulk velocity from the scalar sector intersects the brane through
\[
v_n=v^w-v^aD_a u_w
\]
at the mouth/layer. The normal coupling uses normal pressure/force data. No bulk tangential shear modulus is introduced.

**Defect/throat boundary conditions.** The brane domain is the reference brane with throat mouths excised. At each mouth boundary \(\partial M_A\), the variational boundary work is
\[
\delta W_A=t_A^a\delta u_a+f_A\delta u_w+m_A^{ab}\delta K_{ab}
\]
or a Dirichlet anchoring limit of the same data. Current conservation uses either closed defect worldlines on the brane or a matched boundary flux into the throat. At spatial infinity, variations are compactly supported or radiative so the source-current boundary term vanishes.

**Residual symmetries.** Away from mouth data the candidate keeps time translations, brane Euclidean symmetry, and internal material homogeneity/isotropy in the brane directions. The embedding singles out the normal \(w\) direction, but the declared brane principal blocks contain only \(k_x,k_y,k_z\), not \(k_w\). Preferred-frame effects from the \(w\) embedding and the bulk flow remain later-stage tests.

**Projection firewall.** The old `Z(w)`/`W(w)` separation is retained and extended:

| profile | role |
| --- | --- |
| \(Z(w)\) | frozen localized-Maxwell kinetic profile in `pde.tex`; belongs to the old \(A_M\) dynamics |
| \(W(w)\) | brane-observer projection kernel in `pde.tex` |
| \(B_\ell(w)\) | new brane material/thickness profile for elastic density and modulus |

No equality among \(Z\), \(W\), and \(B_\ell\) is allowed without a later controlled choice. The mixed fields \(A_w,J^w,F_{\mu w},E_w,C_a\) from the frozen parent are not erased by this Stage-0 brane action.

## Microstructure contract

Path chosen: `POSTULATED`, not independently derived.

Reason: the available physical record motivates "deeper substructure/cohesion" beneath the GP/NLS mean field, but it does not specify independent internal variables, a symmetry class, a coarse-graining rule, or an invariant list that is independently fixed before seeing the desired MacCullagh form. Choosing spinning gyrostats solely to obtain curl-only elasticity would be reverse-engineered for this directive. Stage 0 therefore does not claim a microstructure derivation.

Consequence: the conditional-verdict rule is active. If later stages recover Maxwell structure using the rotational candidate, the honest statement would be "Maxwell structure follows from a postulated rotational elasticity," not "the medium derives electromagnetism."

Allowed invariant list carried to Stage 2:

| candidate | variables | invariants |
| --- | --- | --- |
| Cauchy/Navier | \(u^a\) | \((D_au^a)^2\), \(e_{ab}e^{ab}\) |
| rotational/MacCullagh | \(u^a\) | \(r_ar^a=(\nabla_3\times u)^2\) |
| Cosserat/micropolar | \(u^a,\varpi_a\) | \(e_{\langle ab\rangle}e^{\langle ab\rangle}\), \((2\varpi_a-r_a)^2\), \(D_a\varpi_bD^a\varpi^b\), parity-allowed higher terms if declared |
| bending | \(u_w\) | \(D_au_wD^au_w\), \((D^aD_au_w)^2\), \(u_w^2\) |

## Coupling provenance

| term | provenance | status |
| --- | --- | --- |
| \(u_w\Pi_n[\psi,\Sigma]\) | symmetry-allowed normal work term | `SYMMETRY_ALLOWED_POSTULATED`; coefficient and exact bulk functional not derived |
| \(J_{\rm br}^aA_a^{\rm br}-J_{\rm br}^0\Phi_{\rm br}\) | gauge-admissible source structure for conserved defect current | `SYMMETRY_ALLOWED_POSTULATED`; not a charge derivation |
| \(J_{\rm br}^\mu\) conservation | defect worldline continuity or matched throat flux | structurally derived from no-endpoint current assumption, but charge normalization is deferred |
| \(A_a^{\rm br}=\alpha_u u_a\) | admissibility parameterization | postulated; not the Stage-5 EM dictionary |

## Current-admissibility pre-check

Result: `ADMISSIBLE_WITH_CONSERVED_DEFECT_CURRENT_AND_PHI_COMPLETION`. Therefore Stage 0 does not trigger `FAIL_NO_ADMISSIBLE_CURRENT`.

The bar cleared here is only that a conserved current is structurally writable in the coupling. \(\Phi_{\rm br}\) is a
non-dynamical auxiliary in this pre-check; its mechanical origin and the C5 kinetic-term gauge obstruction are deferred to
Stage 5. Whether the defect/throat sector actually supplies such a conserved current is the able-to-fail Stage-6a test, so
this pre-check cannot fire for that failure mode by design.

For the source term
\[
S_J=\int dt\,d^3\xi\,\sqrt{\gamma}\,
\left(J_{\rm br}^aA_a^{\rm br}-J_{\rm br}^0\Phi_{\rm br}\right),
\]
with
\[
\delta A_a^{\rm br}=D_a\chi,\qquad
\delta\Phi_{\rm br}=-\partial_t\chi,
\]
the variation is
\[
\delta S_J
=
\int dt\,d^3\xi\,\sqrt{\gamma}\,
\left(J_{\rm br}^aD_a\chi+J_{\rm br}^0\partial_t\chi\right)
=
-\int dt\,d^3\xi\,\sqrt{\gamma}\,
\chi(\partial_tJ_{\rm br}^0+D_aJ_{\rm br}^a)
+\text{boundary}.
\]
With \(\partial_tJ_{\rm br}^0+D_aJ_{\rm br}^a=0\) and the boundary condition above, the source coupling is gauge-invariant.

A defect current can have this structure:
\[
J_{\rm br}^\mu(x,t)
=
\sum_A q_A\int d\tau\,
\dot X_A^\mu(\tau)\,
\frac{\delta^{(4)}(x-X_A(\tau))}{\sqrt{\gamma}},
\qquad
q_A=\eta_{Q,A}e_\star
\]
up to the later normalization map. Closed worldlines or matched throat flux give the conservation law. The charge sign is carried by \(\eta_Q\), not circulation.

Negative controls from both scripts:

- \(J^a u_a\) without \(\Phi_{\rm br}\) is not gauge-admissible for a generic time-dependent charge current.
- The \(\Phi_{\rm br}\)-completed term is not gauge-invariant for an unconserved current.

## Constitutive-form menu and DOF audit

Constitutive menu:

| candidate | raw in-plane content | known Stage-0 risk |
| --- | --- | --- |
| Cauchy/Navier | 3 in-plane displacement components | generic longitudinal elastic wave |
| rotational/MacCullagh | 3 in-plane displacement components with curl-only potential | C5 kinetic obstruction: curl-only potential is not full-action gauge |
| Cosserat/micropolar | \(u^a\) plus micro-rotation \(\varpi_a\) | extra micro-rotation modes and gap-tuning risk |

Raw brane field count:

| sector | raw variables | Stage-0 count |
| --- | --- | --- |
| in-plane \(u^a\) | \(u_x,u_y,u_z\) | 3 |
| off-brane \(u_w\) | scalar bending | 1 |
| source completion \(\Phi_{\rm br}\) | auxiliary temporal/scalar-potential analog | 1 auxiliary |

The brane-elastic sector therefore has 4 raw elastic configuration fields, or 5 raw variables if the current-admissibility scalar \(\Phi_{\rm br}\) is included. This is not yet a propagating-mode count. A later rotational/gauge result could reduce the in-plane sector to two transverse modes, but Stage 0 does not claim that. Cauchy generically carries a longitudinal mode. Cosserat may add modes.

Comparison to frozen \(A_M\):

| parent field | components | physical count before brane reduction |
| --- | --- | --- |
| \(A_M\) in 4+1 | 5 | massless gauge field has \(D-2=3\) physical polarizations |
| strict brane zero-mode \(A_\mu\) | 4 | Maxwell gauge/constraint count gives 2 physical polarizations |
| brane \(u^a\) candidate | 3 | no gauge reduction proved in Stage 0 |

Classification: `NEW_PARENT_ACTION`. The brane-elastic EM sector is not derived from the frozen \(A_M\) action and is not equivalent to it by field redefinition at Stage 0.

## Assumptions and provenance

| item | Stage-0 assumption/status | provenance |
| --- | --- | --- |
| Flat 4+1 bulk chart | kept | `research/pde/paper/pde.tex:242-255` |
| Bulk scalar sector | kept as written | `pde.tex:316-348` |
| Frozen localized Maxwell sector | reference/comparison, not edited | `pde.tex:257-261`, `pde.tex:357-365`, `pde.tex:410-427` |
| Charge sign | \(\eta_Q\), not circulation | `pde.tex:279-312` |
| \(v\leftrightarrow A\), \(\Omega\leftrightarrow F\) link | frozen-parent context only | `pde.tex:435-490` |
| Projection firewall | \(Z\), \(W\), and \(B_\ell\) kept distinct | `pde.tex:277-278`, `pde.tex:492-565` |
| No bulk shear constraint | brane elasticity only; no bulk tangential modulus introduced | `software/stage1_solver/decisions/15_em_medium_native_physical_picture.md:43-48`, `:63-74` |
| MacCullagh template and C5 obstruction | acknowledged, not resolved | `decisions/15_em_medium_native_physical_picture.md:128-179` |
| Conceptual costs C1/C9/C10 | carried | `decisions/15_em_medium_native_physical_picture.md:216-239` |
| Stage-0 acceptance items | action, contracts, provenance, current pre-check, DOF/classification | `software/stage1_solver/directives/pathA_23_em_medium_native.md:130-149` |
| Discipline | able-to-fail, dimensional units, no paper/decision edits | `pathA_23_em_medium_native.md:230-254` |

## Open questions handed to later stages

1. Stage 1: Does the coupling structure leak into bulk transverse shear, Magnus, or gravity-flow channels, including through \(u_w\)?
2. Stage 2: Can any independently motivated substructure derive Cauchy, rotational, or Cosserat elasticity, or is the constitutive law only postulated?
3. Stage 2/4: Are the Cauchy longitudinal mode, Cosserat micro-rotation modes, and \(u_w\) scalar gapped, constrained, confined, or fatal?
4. Stage 3: Does the post-constitutive traction/couple-stress closure source bulk transverse stress?
5. Stage 5: Does \(\Phi_{\rm br}\) have a mechanical origin that fixes the C5 kinetic-term obstruction, or does the full action fail gauge structure?
6. Stage 6a: Can \(J_{\rm br}^\mu\) be derived from the defect/throat coupling with the \(\eta_Q e_\star\) charge ontology and no source double-counting?
7. Stage 6b: If the frozen \(A_M\) sector remains in the canonical parent, how is double-counting avoided?
8. Stage 6c: Are \(c_\gamma^2=\mu_{\rm br}/\rho_\parallel\) and the bulk \(c_s^2\) tied by common microparameters, or not?

## Script pointers

Primary Mathematica:

- `software/stage1_solver/tools/pathA_23_stage0_action_contracts.wl`
- Output: `software/stage1_solver/_scratch/pathA_23_stage0_action_contracts_mathematica.json`
- Run: `timeout 600 math -script software/stage1_solver/tools/pathA_23_stage0_action_contracts.wl`
- Result: 25/25 checks, exit 0.

Independent SymPy cross-check:

- `software/stage1_solver/tools/pathA_23_stage0_action_contracts_sympy.py`
- Output: `software/stage1_solver/_scratch/pathA_23_stage0_action_contracts_sympy.json`
- Run: `timeout 600 python3 software/stage1_solver/tools/pathA_23_stage0_action_contracts_sympy.py`
- Result: 25/25 checks, exit 0.

What the scripts verify: restored-unit homogeneity for bulk, brane, thickness, bending, normal coupling, and current coupling; source-current gauge variation with negative controls; no explicit \(k_w\) in the declared brane principal blocks; raw DOF arithmetic against the frozen \(A_M\) count.

Honesty note: of the 25 dual-engine checks, roughly half are bookkeeping. The hardcoded DOF integer arithmetic, the
no-\(k_w\) transcription checks built from operators that contain no \(k_w\) by construction, and several definition-echo
dimensional checks are not load-bearing. The genuinely physics-discriminating content is the 3 gauge/current-admissibility
checks with falsifiable negative controls plus the roughly 8 substantive dimensional compositions. This is the appropriate
posture for a specification stage, but "25/25" should not be read as deep physics verification.

What the scripts do not verify: no-leak, constitutive provenance, stress/couple-stress closure, Hamiltonian positivity, mode spectrum, Maxwell equations, Lorentz/E/B mixing, charge normalization, or cone lock.
