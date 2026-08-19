# Opposite-Orientation Throats, Projected Coincidence, and Core Contact

## Status and purpose

This document records a geometric extension and the calculations it requires. It is a derivation plan, not a claim that the required solutions already exist. It supplements, rather than replaces, the canonical [ontology and closure ledger](toy_model_ontology_summary.md). Definitions of flux, gravity-source symbols, projected and full coordinates, and epistemic status are imported from that ontology; the equations repeated here are only those needed to specify the focused two-interface calculation.

This plan is conceptually frozen with the ontology: revisions should now be driven by the slab, spectrum, Green-matrix, source-profile, and nonlinear-throat calculations rather than by adding untested mechanisms.

Definitions repeated from the ontology are guarded by a [mechanical synchronization check](../software/check_ontology_shared_definitions.py); that check detects documentation drift and does not validate the physics.

The central observation is simple:

> Two defects can occupy the same projected three-dimensional position while remaining distinct in the normal direction of a finite-thickness brane.

This matters for opposite electric orientations, particle–antiparticle identification, short-range electric behavior, annihilation, neutral composites, and the adequacy of a single normal-displacement field.

The framework below preserves the following distinctions:

- charge orientation is not particle species;
- projected overlap is not four-dimensional core overlap;
- cancellation of a long-range electric disturbance is not disappearance of the nonlinear material cores;
- a same-species particle–antiparticle pair is not the same problem as two oppositely charged objects with different internal or composite structure;
- linear theory can determine mode content and finite-thickness kernels, but contact, reconnection, and annihilation require a nonlinear time-dependent calculation.

## 1. Particle branch, orientation, and antiparticle

Let a throat species be labeled by \(\alpha\). The label includes every orientation-independent part of the object:

- trapped support-mode family;
- core topology and boundary data;
- equilibrium radius and shape;
- constituent-transport and order-conversion structure;
- any additional conserved or metastable internal structure.

Let

\[
\mathcal T_{\alpha,s},
\qquad s\in\{+1,-1\},
\]

denote species \(\alpha\) with normal orientation \(s\). Reversing \(s\) on the same solved internal branch defines the leading antiparticle candidate only if the complete operation transforms every required orientation-odd internal datum and preserves every orientation-even scalar datum. The target map is

\[
\mathcal C_{\rm th}
=\mathcal I_{\rm internal}\circ\mathcal R_w:
\mathcal T_{\alpha,+}
\longleftrightarrow
\mathcal T_{\alpha,-}.
\]

Here \(\mathcal R_w\) is geometric reflection of the parent fields, while \(\mathcal I_{\rm internal}\) transforms any additional odd phase, circulation, chirality, support-mode convention, or topological datum. Whether \(\mathcal I_{\rm internal}\) is trivial is an output for each solved species branch; geometric reflection alone is not assumed to be the full antiparticle operation.

For a reversible equilibrium branch, orientation-even magnitudes evaluated with reflected environmental data must include

\[
E_\alpha^{(+)}=E_\alpha^{(-)},
\qquad
Q_{n,\alpha}^{(+)}=Q_{n,\alpha}^{(-)},
\qquad
Q_{\chi,\alpha}^{(+)}=Q_{\chi,\alpha}^{(-)},
\]

while all required orientation-odd data reverse. Chirality, circulation, knot orientation, support-mode phase convention, topological index, or conversion-current data may therefore have to transform as part of the complete map. Electric-orientation reversal alone has not yet established that map.

For a conservative periodic branch, the complete map must additionally transform support-mode action, frequency, and phase data. Compare a derived cycle-averaged invariant or the cycle-averaged stress and momentum flux in the same stated ensemble; two throats with arbitrary relative phase are not automatically reflection partners.

For a driven relaxational or mixed branch, equality of an equilibrium energy may be unavailable or insufficient. The reflected stationary tests must instead include, as applicable, equality of positive throughput, dissipation rate, entropy-production rate, lifetime, stationary stress and momentum fluxes, and the complete retarded response kernels under reflected environmental data.

The environment must transform with the throat. For any source-side gravity comparison the target covariance imported from the ontology is

\[
\mathbf g_{\rm eff}
[\mathcal C_{\rm th}\mathcal T;
 \mathcal R_w\mathcal B_\infty]
=
\mathbf g_{\rm eff}
[\mathcal T;\mathcal B_\infty].
\]

Degeneracy in one unchanged environment additionally requires \(\mathcal R_w\mathcal B_\infty=\mathcal B_\infty\). Unequal bulk-facing reservoirs may produce a calculable environmental splitting and do not by themselves disprove the intrinsic candidate throat operation.

Therefore, reversing an electron-like throat produces a positron-like candidate, not a proven complete antiparticle solution. It does not produce a proton-like object. A proton-like object must instead occupy a distinct internal branch, possess composite throat structure, or carry additional topology not removed by electric-orientation cancellation.

This distinction prevents the model from predicting that every positive and negative object must annihilate merely because their electric orientations are opposite.

Flux notation follows the ontology, including \(\partial_t n+\nabla_4\cdot\mathbf J_n=0\). Let \(\mathcal A_\alpha\) be its selected \(w\)-normal three-dimensional reference cut through an axis-aligned throat. With \(\hat{\mathbf n}_\alpha=s\hat{\mathbf w}\),

\[
\widetilde Q_{n,\alpha}^{(w)}
=\int_{\mathcal A_\alpha}
\mathbf J_n\cdot\hat{\mathbf w}\,d^3\Sigma,
\qquad
Q_{n,\alpha}^{\rm net}
=\int_{\mathcal A_\alpha}
\mathbf J_n\cdot\hat{\mathbf n}_\alpha\,d^3\Sigma.
\]

Thus \(\widetilde Q_{n,\alpha}^{(w)}\) is signed in the global \(w\) convention and \(Q_{n,\alpha}^{\rm net}\) is signed local net outward flux. On this reference cut the normal definitions give the **identity**

\[
\boxed{
\widetilde Q_{n,\alpha}^{(w)}
=sQ_{n,\alpha}^{\rm net}
}.
\]

For a general curved cut with normal field \(\hat{\mathbf n}(\mathbf X)\), use the invariant surface flux; the simple global-\(w\) identity then requires a derived projected-area relation and is not automatic.

For each independently solved stationary outward-drain branch, define \(Q_{n,\alpha}\equiv Q_{n,\alpha}^{\rm net}>0\) only after its sign is verified. Under reflection-symmetric boundary and reservoir data, the actual physical targets for \(\mathcal C_{\rm th}\)-related partners are

\[
\boxed{
Q_{n,\alpha}^{(+)}=Q_{n,\alpha}^{(-)}
},
\qquad
\boxed{
\widetilde Q_{n,\alpha}^{(w)}(+)
=-\widetilde Q_{n,\alpha}^{(w)}(-)
}.
\]

The quantity \(Q_{\chi,\alpha}\ge0\) remains the separate gross local ordered-to-bulk drain magnitude. In the fraction/relaxational branch it is \(\int_{\Omega_{{\rm throat},\alpha}}n\Gamma_{\rm drain}\,d^4X\), with the same number-per-time dimensions as \(Q_{n,\alpha}\); in an inertial order-field branch the corresponding gross diagnostic must be derived rather than imported.

## 2. Finite slab and two interface coordinates

Let the unperturbed ordered slab have thickness \(H_0\equiv H_{\rm br}\), where \(H_{\rm br}\) is the homogeneous background value used by the ontology. Its interfaces are

\[
w_+^{(0)}=+\frac{H_0}{2},
\qquad
w_-^{(0)}=-\frac{H_0}{2}.
\]

Introduce independent interface displacements

\[
w_+(\mathbf x,t)=+\frac{H_0}{2}+h_+(\mathbf x,t),
\qquad
w_-(\mathbf x,t)=-\frac{H_0}{2}+h_-(\mathbf x,t).
\]

The fields \(h_+\) and \(h_-\) are collective coordinates for the upper and lower interfaces only while those interfaces can be represented as graphs over the brane coordinates. They are appropriate for background modes, far fields, finite-thickness response, and weakly to moderately nonlinear separated-interface configurations. They do not by themselves represent a complete throat core once the geometry develops overhangs, necks, tubes, folded sheets, reconnection, pinch-off, or topology change. Explicitly,

\[
\boxed{
h_\pm
=\text{far- and intermediate-field interface coordinates}
}
\]

whereas

\[
\boxed{
\begin{gathered}
\chi_B(\mathbf x,w),\ n(\mathbf x,w),\ \theta(\mathbf x,w),\\
\text{support fields, and core level sets}
\end{gathered}
=\text{full nonlinear throat geometry}
}.
\]

Thus two interface coordinates extend the earlier one-height description; they do not replace the full four-dimensional core solution.

The natural parity combinations are

\[
h_m=\frac{h_++h_-}{2},
\qquad
h_t=\frac{h_+-h_-}{2}.
\]

They have direct geometric meanings:

\[
w_{\rm mid}=h_m,
\qquad
H(\mathbf x,t)=H_0+2h_t.
\]

Under reflection \(w\to-w\),

\[
h_+\to-h_-,
\qquad
h_-\to-h_+,
\]

and consequently

\[
h_m\to-h_m,
\qquad
h_t\to+h_t.
\]

Thus \(h_m\) is the reflection-odd mid-surface displacement and \(h_t\) is the reflection-even thickness disturbance.

Translation invariance in the parent \(w\) direction requires a uniform \(h_m\) shift to cost no energy on an isolated homogeneous background. A dynamically selected slab thickness may give \(h_t\) a positive restoring curvature. This already shows why one height field is generally insufficient: it suppresses either the second interface or one of these parity channels by construction.

## 3. Oriented interface protrusions

For illustration, let \(a_c>0\) be a core displacement amplitude. A one-sided outward protrusion can have the schematic boundary data

\[
s=+1:
\qquad
(h_+,h_-)=(a_c,0),
\]

and

\[
s=-1:
\qquad
(h_+,h_-)=(0,-a_c).
\]

In the parity basis these become

\[
s=+1:
\qquad
(h_m,h_t)=\left(+\frac{a_c}{2},+\frac{a_c}{2}\right),
\]

\[
s=-1:
\qquad
(h_m,h_t)=\left(-\frac{a_c}{2},+\frac{a_c}{2}\right).
\]

The orientation-odd component reverses, while the thickness component does not.

For an equal-amplitude opposite-orientation pair at the same projected position,

\[
h_m^{(+)}+h_m^{(-)}=0,
\qquad
h_t^{(+)}+h_t^{(-)}=a_c.
\]

Thus the illustrative total pair has

\[
h_m^{\rm tot}=0,
\qquad
h_t^{\rm tot}=a_c.
\]

This is the core of the proposal. The signed far disturbance may cancel, while an even thickness deformation, stored stress, conversion zone, and two support-mode configurations remain. Electrical neutralization is therefore not identical to material disappearance.

These boundary values are illustrative rather than frozen and are not a solved throat profile. The nonlinear parent-field solution must determine whether a physical core is one-sided, spans the slab, displaces both interfaces, or requires a geometry that cannot be represented as two single-valued graphs.

## 4. Projected position versus core support

A brane observer assigns a throat the projected position

\[
\mathbf x_a=(x_a,y_a,z_a).
\]

A general full four-dimensional spatial point is \(\mathbf X=(\mathbf x,w)\).

The canonical arena metric is the ambient Euclidean \(\delta_{AB}\), so \(\|\mathbf Y_+-\mathbf Y_-\|_4\) below is Euclidean distance. It is independent of the interface-graph parameterization, not of arbitrary changes to the ambient metric; a nonflat parent branch would require geodesic distance instead.

Because the parent fields may have noncompact tails, “core” cannot mean the literal support of every disturbance. Let \(\mathcal I_a(\mathbf X)\) be a solved localized diagnostic such as order depletion, defect density, conversion intensity, or another invariant selected by the parent solution, and define

\[
\boxed{
\mathcal C_a(\epsilon)
=\left\{
\mathbf X:\mathcal I_a(\mathbf X)\ge\epsilon
\right\}
\subset\mathbb R^3_{\mathbf x}\times\mathbb R_w.
\]

The diagnostic and threshold range must give compact core sets, or at least closed localized sets for which the distance infimum is attained. The separate \(+/-\) labels must also be derived rather than imposed. Use one of two procedures:

1. **Orientation-resolved diagnostics:** derive genuinely distinct \(\mathcal I_+\) and \(\mathcal I_-\), for example from an oriented topological or conversion-density field. Their threshold sets may be intersected directly.
2. **Connected-component continuation:** use one total diagnostic \(\mathcal I(\mathbf X)\), identify two thresholded connected components at large \(R\), and continue them toward smaller \(R\). Contact or merger occurs when their closures first touch or their component count or topology changes.

Below, \(\mathcal C_a\) abbreviates a physical set obtained by one of these procedures unless \(\epsilon\) is shown explicitly. Once a merger leaves only one component and no orientation-resolved diagnostic, separate \(\mathcal C_+\) and \(\mathcal C_-\) cease to exist; they must not be reconstructed by an analyst-chosen decomposition of the total field.

Two cores can therefore satisfy

\[
\mathbf x_+=\mathbf x_-
\]

without satisfying

\[
\mathcal C_+\cap\mathcal C_-\ne\varnothing.
\]

In the interface-graph regime, let the total local interface positions in the two-core configuration be

\[
w_+^{\rm tot}(\mathbf x)
=+\frac{H_0}{2}+h_+^{\rm tot}(\mathbf x),
\qquad
w_-^{\rm tot}(\mathbf x)
=-\frac{H_0}{2}+h_-^{\rm tot}(\mathbf x).
\]

Let \(\ell_+^{\rm in}(\mathbf x)\) be the inward reach of a core anchored on the upper interface and \(\ell_-^{\rm in}(\mathbf x)\) the inward reach of a core anchored on the lower interface. The graph-regime normal gap is

\[
\begin{aligned}
d_w^{\rm core}(\mathbf x)
&=\left[w_+^{\rm tot}(\mathbf x)-\ell_+^{\rm in}(\mathbf x)\right]
-\left[w_-^{\rm tot}(\mathbf x)+\ell_-^{\rm in}(\mathbf x)\right]\\[1mm]
&=H_0+h_+^{\rm tot}(\mathbf x)-h_-^{\rm tot}(\mathbf x)
-\ell_+^{\rm in}(\mathbf x)-\ell_-^{\rm in}(\mathbf x)\\[1mm]
&=H_0+2h_t^{\rm tot}(\mathbf x)
-\ell_+^{\rm in}(\mathbf x)-\ell_-^{\rm in}(\mathbf x).
\end{aligned}
\]

At zero projected separation:

- \(d_w^{\rm core}>0\): the projected mouths overlap, but the four-dimensional cores remain disjoint;
- \(d_w^{\rm core}=0\): the cores first touch;
- \(d_w^{\rm core}<0\): the separate-core graph ansatz has overlapped and must be replaced by a merged nonlinear parent-field solution.

The quantities \(\ell_\pm^{\rm in}\) must be outputs of the core solution. They cannot be chosen merely to enforce or avoid annihilation.

For the illustrative outward pair, \(h_+^{\rm tot}=a_c\), \(h_-^{\rm tot}=-a_c\), and \(h_t^{\rm tot}=a_c\). With negligible inward reach,

\[
d_w^{\rm core}\simeq H_0+2a_c,
\]

not merely \(H_0\). In that restricted geometry the two cores cannot touch. Contact requires the solved cores or deformation zones to extend sufficiently toward the mid-surface, one core to span the slab, or the interfaces to undergo a topology change. This is a conditional geometric result, not yet a solved property of the model's throats.

Once graph coordinates fail, use the interface-graph-independent ambient core-set distance while two cores remain identifiable:

\[
d_{\rm core}(\epsilon)
=\inf_{\substack{
\mathbf Y_+\in\mathcal C_+(\epsilon)\\
\mathbf Y_-\in\mathcal C_-(\epsilon)}}
\left\|\mathbf Y_+-\mathbf Y_-\right\|_4.
\]

This definition remains meaningful when the interfaces develop overhangs, necks, or folds. For orientation-resolved sets, actual contact is

\[
\boxed{
\mathcal C_+(\epsilon)
\cap
\mathcal C_-(\epsilon)
\ne\varnothing
}.
\]

For compact orientation-resolved sets this is equivalent to \(d_{\rm core}(\epsilon)=0\); for nonclosed or nonlocalized sets, zero infimum alone would not prove intersection. Under connected-component continuation, closure contact and the component-count or topology change define merger, after which separate set distance is no longer meaningful. The inferred contact point and qualitative outcome must be robust under reasonable variation of \(\epsilon\). Neither contact nor set distance decides whether the subsequent evolution is reconnection, annihilation, pinch-off, or another conserving topology change.

## 5. Quadratic two-interface action

The following is a minimal local, time-reversal-even, second-order conservative quadratic interface template:

\[
S_{\rm int}^{(2)}
=\frac12\int dt\,d^3x\,
\left[
\dot{\mathbf h}^{\,T}\mathbf M\dot{\mathbf h}
-(\partial_i\mathbf h)^{T}\mathbf K^{ij}(\partial_j\mathbf h)
-\mathbf h^T\mathbf U\mathbf h
\right],
\]

where

\[
\mathbf h=
\begin{pmatrix}
h_+\\
h_-
\end{pmatrix}.
\]

On an isotropic reflection-symmetric background, transformation to \((h_m,h_t)\) should diagonalize reflection parity at quadratic order. A minimal form is

\[
S_{mt}^{(2)}
=\frac12\int dt\,d^3x\,
\left[
\rho_m\dot h_m^2-\kappa_m|\nabla h_m|^2
+\rho_t\dot h_t^2-\kappa_t|\nabla h_t|^2
-\mu_t^2h_t^2
\right].
\]

Within this minimal template, parent translation invariance forbids a local \(\mu_m^2h_m^2\) term. A healthy realization of the displayed thickness-restoring branch would require

\[
\rho_m>0,
\quad
\rho_t>0,
\quad
\kappa_m\ge0,
\quad
\kappa_t\ge0,
\quad
\mu_t^2>0,
\]

after all coupled fields and constraints are included.

This action is only a diagnostic template. A more general effective response may contain mixed time-space derivatives, first-order gyroscopic or time-reversal-odd terms, higher-gradient bending terms, nonlocal kernels generated by integrating out bulk modes, frequency- and wave-number-dependent coefficients, and constrained density and phase variables. The actual spectrum must also include coupling to the order profile \(\chi_B(w)\), longitudinal motion \(u_L\), bulk density and phase, shear motion, and any other parent variables retained by the reduction. Reflection parity may block some mixing on the exact symmetric background, while throats and curvature can reactivate it.

If the selected order-conversion branch is relaxational or explicitly mixed inertial-dissipative, the physical linear calculation is a retarded response problem with frequency-dependent damping and every required internal reservoir variable. The parity decomposition and stability tests still apply, but the dynamic object is the complete constrained retarded operator \(\mathcal O_R(\omega,\mathbf k)\), not a free-energy Hessian. A conservative, relaxational, or mixed branch is permitted only when its material, momentum, energy, inertia/drag, and required entropy/reservoir ledgers are explicit; an unspecified blend is not.

## 6. Linear source and Green-matrix problem

Let \(\boldsymbol\Phi\) collect the relevant linear normal/order variables, beginning with

\[
\boldsymbol\Phi=(h_m,h_t,u_L,\delta\chi_B,\ldots)^T.
\]

### 6.1 Conservative or equilibrium static problem

Only a conservative branch or a genuine equilibrium reduction supplies a free-energy Hessian,

\[
\mathcal H(\mathbf k)
=\frac{\delta^2F}
{\delta\boldsymbol\Phi^\dagger\delta\boldsymbol\Phi},
\qquad
\mathcal G_{\rm stat}(\mathbf k)
=\mathcal H^{-1}(\mathbf k),
\]

on the physical constrained subspace. Its quadratic static functional may be written

\[
F_{\rm stat}^{(2)}
=\frac12\int\frac{d^3k}{(2\pi)^3}
\left[
\boldsymbol\Phi^\dagger
\mathcal H(\mathbf k)
\boldsymbol\Phi
+2\operatorname{Re}\!\left(
\boldsymbol J^\dagger\boldsymbol\Phi
\right)
\right],
\]

provided the chosen boundary ensemble admits that thermodynamic potential. The physical pair potential must be obtained from the correct fixed-value, fixed-source, fixed-flux, mixed, or topological functional, with one-body self terms subtracted.

### 6.2 Dynamic retarded response

Poles, damping, phase lag, radiation, and time-dependent response are instead governed by the retarded operator

\[
\mathcal O_R(\omega,\mathbf k),
\qquad
\mathcal G_R(\omega,\mathbf k)
=\mathcal O_R^{-1}(\omega,\mathbf k).
\]

For a conservative theory, \(\mathcal O_R\) is related to the action's linearized equations and its zero-frequency limit must be consistent with \(\mathcal H\). For a relaxational or mixed theory, \(\mathcal O_R\) may be non-Hermitian and must include every internal reservoir variable needed to close energy, momentum, and entropy accounting. Its inverse is a response function, not automatically the Hessian of an interaction energy.

### 6.3 Conservative intrinsically periodic pair

A self-consistent supported throat may be conservative yet intrinsically time periodic. Such a solution is not automatically governed by the static free energy above. A phase-locked periodic pair may instead admit a cycle-averaged Hamiltonian, Routhian, quasienergy, or averaged action, but only in a derived ensemble that states what is held fixed—for example support-mode action, frequency, relative phase, or the appropriate conjugate data. If no scalar functional is established, derive the mean force from the cycle-averaged Noether stress and momentum flux. Do not substitute a static Hessian merely because the cycle-averaged geometry is stationary.

### 6.4 Driven dissipative stationary pair

For a maintained dissipative nonequilibrium pair, do not assume that a scalar \(E_{ab}(R)\) exists. Derive the quasistatic force from the complete stress tensor, material and order momentum flux, conversion terms, and reservoir response. A potential may be introduced only after the force is shown to be conservative and path independent over the relevant configuration space.

### 6.5 Source parity and pair kernels

For a static oriented throat—or for a declared cycle-averaged source of a phase-locked periodic throat—the source vector should separate orientation-odd and orientation-even data:

\[
\boldsymbol J_{\alpha,s}
=s\,\boldsymbol J_\alpha^{\rm odd}
+\boldsymbol J_\alpha^{\rm even},
\]

where

\[
\boldsymbol J_\alpha^{\rm odd}
=J_{m,\alpha}\mathbf e_m
+\boldsymbol J_{\alpha,\rm other}^{\rm odd},
\qquad
\boldsymbol J_\alpha^{\rm even}
=J_{t,\alpha}\mathbf e_t
+\boldsymbol J_{\alpha,\rm other}^{\rm even}.
\]

The \(h_m\)-like and \(h_t\)-like entries are leading interface projections. Neither parity sector is assumed one-dimensional, and circulation, chirality, phase, topology, order, density, longitudinal, bulk, or reservoir data must enter the appropriate complete block. The orientation-even entries may depend on throat species and support state; \(s\) multiplies the complete odd source.

When orientation-labeled parity contributions are useful, write

\[
\boldsymbol J_{\alpha,s}^{\rm odd}
\equiv s\boldsymbol J_\alpha^{\rm odd},
\qquad
\boldsymbol J_{\alpha,s}^{\rm even}
\equiv\boldsymbol J_\alpha^{\rm even}.
\]

For two sources, the raw quadratic cross-kernel is schematically

\[
\mathcal K_{ab}(R)
=\int\frac{d^3k}{(2\pi)^3}
e^{i\mathbf k\cdot\mathbf R}
\boldsymbol J_a^\dagger(\mathbf k)
\mathcal G_{\rm stat}(\mathbf k)
\boldsymbol J_b(\mathbf k).
\]

This is an equilibrium/conservative static kernel only. The physical interaction energy is not automatically equal to the bare expression: it must use the correct conjugate functional for the boundary ensemble and subtract one-body self terms. Fixed-value, fixed-source, fixed-flux, and mixed boundary data can change the interaction sign. A conservative periodic pair requires the cycle-averaged treatment above, including its action/frequency/phase ensemble; the driven dissipative analog uses \(\mathcal G_R\) to obtain response and the stress/momentum ledger to obtain force. Neither is assigned a static energy by analogy.

On a reflection-symmetric background, the propagator separates into odd and even parity blocks unless a parity-breaking datum is present. The leading kernel then has the schematic block form

\[
\mathcal K_{ab}(R)
\sim
s_as_b
(\boldsymbol J_a^{\rm odd})^\dagger
\mathcal G_{{\rm stat},{\rm odd}}
\boldsymbol J_b^{\rm odd}
+(\boldsymbol J_a^{\rm even})^\dagger
\mathcal G_{{\rm stat},{\rm even}}
\boldsymbol J_b^{\rm even}
+\cdots.
\]

This formula displays a target parity organization only. It does not assign the final force sign or establish the range of either sector. The physical force sign must be obtained from the correctly Legendre-transformed or constrained energy appropriate to the throat boundary data, after one-body self terms are subtracted. Fixed-value, fixed-source, fixed-flux, mixed, or topological boundary conditions need not produce the same on-shell interaction sign.

After all constraints and coupled fields are included, \(\mathcal H_{\rm odd}\) and \(\mathcal H_{\rm even}\) are equilibrium static parity blocks, while \(\mathcal O_{R,{\rm odd}}\) and \(\mathcal O_{R,{\rm even}}\) are their dynamic retarded counterparts. Each may contain multiple eigenbranches rather than one scalar mode. For a conservative/equilibrium Coulomb potential, the desired low-wave-number outputs include at least one healthy odd Hessian eigenvalue

\[
\lambda_E(k)
\in\operatorname{spec}\mathcal H_{\rm odd}(k),
\qquad
\lambda_E(k)
=\kappa_Ek^2+O(k^4),
\qquad
\kappa_E>0,
\]

with nonzero throat coupling for the Coulomb-carrying branch, and at least one primarily thickness-like even eigenvalue

\[
\lambda_t(k)
\in\operatorname{spec}\mathcal H_{\rm even}(k),
\qquad
\lambda_t(k)
=\mu_t^2+\kappa_tk^2+O(k^4),
\qquad
\mu_t^2>0,
\]

for a screened thickness-like branch if the intended construction requires it. These are target static spectral outputs for selected eigenbranches, not consequences of parity alone. The retarded blocks must separately supply poles, damping, and radiation. Stable thickness selection makes a primarily thickness-like mode eligible to be gapped, but even variables may mix with density, order, longitudinal, reservoir, or other gapless variables; multiple odd and even branches may exist. The calculation must inventory every additional long-range branch and its throat coupling. The actual eigenvectors may mix \(h_m,h_t,\delta\chi_B\), density, phase, and other retained fields.

A Coulomb-like static potential in three brane dimensions specifically requires leading \(k^2\) static stiffness:

\[
G_E(0,k)\sim\frac{1}{\kappa_Ek^2},
\qquad
G_E(R)\sim\frac{1}{4\pi\kappa_E R}.
\]

Gaplessness alone is insufficient. A flexural branch with

\[
E_{\rm flex}
\sim\frac{\kappa_4}{2}
\int d^3x\,(\nabla^2\phi)^2
\]

has a \(k^4\) static kernel. With ordinary local inertia it gives \(\omega^2\propto k^4\), but it does not give the Coulomb Green function. The calculation must therefore distinguish a Coulomb-carrying branch with \(k^2\) static stiffness, a flexural branch with \(k^4\) static stiffness, and any nonstandard kinetic structure only if that structure is independently derived.

If the desired even eigenbranch is in fact gapped and approximately local, it may have

\[
G_t(R)\sim\frac{e^{-R/\lambda_t}}{4\pi\kappa_tR},
\qquad
\lambda_t=\frac{\sqrt{\kappa_t}}{\mu_t}.
\]

Thus odd-long-range and even-short-range behavior is a target to be derived from the full coupled spectrum, not an assumption attached to the parity labels.

## 7. Finite-thickness form factors

The two interfaces cannot be replaced by delta functions in \(w\) when \(R\) approaches the slab thickness. In the general multicomponent calculation, let \(\boldsymbol f_a(\mathbf k,w)\) be the solved throat source profile across all retained fields and let \(\boldsymbol{\mathscr G}_{\rm stat}(\mathbf k;w,w')\) be the full matrix parent static Green operator. The conservative/equilibrium kernel is

\[
\mathcal K_{ab}^{\rm 4D}(R)
=\int\frac{d^3k}{(2\pi)^3}e^{i\mathbf k\cdot\mathbf R}
\int dw\,dw'\,
\boldsymbol f_a^\dagger(\mathbf k,w)
\boldsymbol{\mathscr G}_{\rm stat}(\mathbf k;w,w')
\boldsymbol f_b(\mathbf k,w').
\]

The corresponding dynamic calculation uses the full retarded matrix \(\boldsymbol{\mathscr G}_R(\omega,\mathbf k;w,w')\). Neither expression may discard source components merely because the leading interface projection is simple.

For that leading scalar interface projection only, let \(f_+(w)\) and \(f_-(w)\) be normalized interface profiles. Then

\[
G_{+-}(R)
=\int\frac{d^3k}{(2\pi)^3}e^{i\mathbf k\cdot\mathbf R}
\int dw\,dw'\,
f_+(w)\mathscr G_{\rm stat}^{\rm interface}(\mathbf k;w,w')f_-(w').
\]

The same-interface kernels \(G_{++}\) and \(G_{--}\) are defined analogously.

The calculation must distinguish three regimes:

1. **Far field, \(R\gg H_0\).** Only the lightest projected eigenmodes survive. This is where a \(1/R\) electric potential and possible leading factorization into one-throat magnitudes can emerge.
2. **Finite-thickness crossover, \(R\sim H_0\).** Interface profiles and massive modes modify the kernel. These terms belong to finite-size and environmental corrections rather than the leading universal charge.
3. **Projected coincidence, \(R\to0\).** Linear Green functions can diagnose whether the quadratic response remains finite or singular, but they cannot replace the graph-regime gap, the nonlinear set distance, or the parent-field analysis of merger, reconnection, conversion, or annihilation.

No particular soft-core formula should be assumed in advance. The parent operator and solved profiles must determine the crossover.

## 8. Same-species particle–antiparticle pair

Assuming the brane is globally two-sided and co-orientable, or that an equivalent oriented intersection number exists, the candidate count for a same-species pair \(\mathcal T_{\alpha,+}\) and \(\mathcal T_{\alpha,-}\) is

\[
N_s=+1-1=0.
\]

Within the linear source description, exact projected coincidence of a complete same-species reflection pair gives cancellation of the entire odd source vector,

\[
\boldsymbol J_{\alpha,+}^{\rm odd}
+\boldsymbol J_{\alpha,-}^{\rm odd}=0.
\]

For intrinsically periodic throats, this statement is exact only when \(\boldsymbol J^{\rm odd}\) denotes the cycle-averaged stationary source or when the complete \(\mathcal C_{\rm th}\) operation supplies the required instantaneous phase map. An arbitrary relative support-mode phase may leave transient odd multipoles, phase-dependent stresses, or sideband radiation even when the cycle-averaged Coulomb monopole cancels.

The complete even data add rather than cancel:

\[
\boldsymbol J_{\alpha,+}^{\rm even}
+\boldsymbol J_{\alpha,-}^{\rm even}
=2\boldsymbol J_\alpha^{\rm even}.
\]

If this oriented count is derived and conserved, its zero value permits the pair to relax into neutral configurations. It does not determine which of the following occurs:

- elastic scattering;
- a metastable neutral bound state analogous to a particle–antiparticle atom;
- core contact followed by reconnection;
- annihilation into propagating transverse, longitudinal, order, density, or bulk modes;
- partial annihilation with a neutral remnant;
- exclusion caused by an even core deformation or another conserved quantity.

For a conservative or genuine equilibrium branch, a useful target organization—not a derived potential—is

\[
E_{\alpha\bar\alpha}(R)
=E_{\rm odd}(R)
+E_{\rm even}(R)
+E_{\rm mix}(R)
+E_{\rm core}^{\rm nl}(R).
\]

On an exactly reflection-symmetric background, the mixed term should vanish at quadratic order. After the physical electric sign has been derived from the correct boundary ensemble, the leading long-range odd contribution should be attractive for opposite orientation. The even contribution may be attractive, repulsive, or sign-changing and must be calculated. The nonlinear core term decides contact, reconnection, annihilation, saturation, or exclusion; no desired short-range barrier may be inserted by hand. For a conservative periodic pair, replace this static energy by the derived cycle-averaged functional at fixed support-mode action, frequency, relative phase, or other stated conjugate data; absent such a functional, use cycle-averaged stress and momentum flux. For a driven dissipative pair, use the same labels only to organize response channels and derive force from the stress, momentum-flux, conversion, and reservoir ledgers unless conservativity is independently established.

A viable electron-like branch must possess an allowed electron–positron annihilation channel, but annihilation need not be instantaneous or compulsory in every encounter. The time-dependent calculation must conserve total constituent content, momentum, and energy; in a relaxational or mixed realization it must also account for entropy and internal reservoir excitation.

## 9. Opposite charge with different species

Consider \(\mathcal T_{\alpha,-}\) and \(\mathcal T_{\beta,+}\) with \(\alpha\ne\beta\). Their leading electric monopoles may cancel, but neither every odd internal component nor the even source data need cancel:

\[
\boldsymbol J_{\alpha,-}
+\boldsymbol J_{\beta,+}
=-\boldsymbol J_\alpha^{\rm odd}
+\boldsymbol J_\beta^{\rm odd}
+\boldsymbol J_\alpha^{\rm even}
+\boldsymbol J_\beta^{\rm even}.
\]

Species-dependent odd multipoles or internal structure can remain even when the leading electric monopole vanishes. The even data generally add, and the pair may retain a nonzero neutral core, distinct support-mode content, or additional topology. Electric neutrality alone must not force annihilation.

In particular, a proton-like object cannot be modeled solely as the \(s=+1\) orientation of the electron-like branch. That construction is the positron-like branch. A proton-like sector requires separate internal or composite structure, and the theory must explain why electron–proton overlap can form a bound neutral system rather than following the same annihilation channel as electron–positron overlap.

## 10. Linear derivation program

The background and linear calculations should proceed in the following order.

### 10.1 Solve the finite slab

Determine a stationary \(n_0(w),\chi_0(w)\) profile with two interfaces and a selected thickness \(H_0\). Verify

\[
\chi_0(w)=\chi_0(-w),
\qquad
\chi_0'(0)=0,
\qquad
\chi_0(0)=\max_w\chi_0(w),
\]

and establish whether the solution is stable or metastable.

Record the asymptotic data \(\mathcal B_\infty\), test whether the ordinary vacuum satisfies \(\mathcal R_w\mathcal B_\infty=\mathcal B_\infty\), and retain an explicitly asymmetric environment as a separate calculable comparison rather than folding it into intrinsic throat parity.

### 10.2 Derive the interface collective coordinates

Project perturbations onto upper- and lower-interface translation modes to obtain \(h_+\) and \(h_-\) in the graph regime. Transform to \(h_m,h_t\), derive their kinetic and gradient matrices, and verify reflection parity mechanically rather than assigning it. State explicitly where these collective coordinates cease to resolve the parent-field core geometry.

### 10.3 Diagonalize the complete coupled spectrum

Include \(\mathbf u_T\), \(u_L\), \(\delta\chi_B\), density/phase perturbations, normal bulk motion, and every internal variable required by the selected conservative, relaxational, or mixed branch. For an equilibrium branch, diagonalize the complete odd and even Hessian blocks \(\mathcal H_{\rm odd/even}\). Separately analyze the retarded blocks \(\mathcal O_{R,{\rm odd/even}}\) for dynamics; do not assume one scalar branch in either parity sector. Determine:

- which combination is the gapless normal/electric candidate;
- whether a healthy odd eigenmode has the required \(k^2\) static stiffness;
- whether a thickness-like even eigenmode is actually gapped after all mixing is included;
- whether any additional odd or even long-range eigenbranches have nonzero throat coupling;
- whether any \(k^4\) flexural branch is distinct from the Coulomb carrier;
- whether ghosts or unstable eigenvalues occur;
- whether the transverse light modes remain isolated;
- whether defects or curvature produce unacceptable odd–even or longitudinal–transverse mixing.

### 10.4 Derive throat source profiles

For each species and orientation, derive the interface/core form factors and the complete decomposition \(\boldsymbol J_{\alpha,s}=s\boldsymbol J_\alpha^{\rm odd}+\boldsymbol J_\alpha^{\rm even}\). Verify the full candidate operation \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\), not only its leading \(h_m\)- and \(h_t\)-like projections. If asymptotic electric factorization is obtained, the electric sign appears only through \(s\), while one-throat coupling magnitudes remain nonnegative after branch selection.

### 10.5 Compute the Green matrix

Compute the full multicomponent \(\boldsymbol{\mathscr G}_{\rm stat}(\mathbf k;w,w')\) for equilibrium interactions and \(\boldsymbol{\mathscr G}_R(\omega,\mathbf k;w,w')\) for dynamic response. Treat scalar \(G_{++},G_{--},G_{+-}\), or \(h_m,h_t\) propagators, only as leading interface projections. Extract:

- the \(R\gg H_0\) asymptotics;
- the \(R\sim H_0\) crossover;
- the behavior of the quadratic kernel as \(R\to0\);
- the size of odd–even mixing away from the symmetric background.

### 10.6 Test electric factorization

For multiple independently solved throat species and support states, compute the asymptotic coefficient matrix of the reflection-odd Coulomb monopole channel only. After removing \(s_as_b\), determine whether its leading part is positive semidefinite and approximately rank one:

\[
C_{E,ab}(R)
=g_ag_b+\delta C_{ab}^{\rm odd}(R),
\qquad
g_a\ge0.
\]

For pairs with \(g_ag_b>0\), require both \(\delta C_{ab}^{\rm odd}/(g_ag_b)\to0\) and \(R\partial_R\delta C_{ab}^{\rm odd}/(g_ag_b)\to0\) in the far-field regime before inferring the leading force sign. Compute reflection-even, parity-mixed, and nonlinear core interactions separately; none belongs inside \(C_{E,ab}\). Common slab thickness alone does not establish factorization.

### 10.7 Accelerating-source and radiation test

If the Coulomb-carrying odd eigenmode is dynamical rather than constrained, an accelerating oriented throat may radiate into it. The complete moving and accelerating throat calculation must measure

\[
P_{u_T},
\quad
P_{\rm odd},
\quad
P_{\rm even},
\quad
P_{u_L},
\quad
P_{\rm bulk},
\quad
P_\chi,
\]

corresponding to the two transverse light modes, the odd normal/electric sector, the even thickness sector, the in-plane longitudinal mode, density or bulk-sound modes, and order-conversion or reservoir modes. The odd electric mode must not be assumed nondynamical before the complete operator is derived.

### 10.8 Support and gravity-source tests

On each provisional one-core geometry, solve the complete variable-coefficient transverse eigenproblem as a spectral seed. Require finite weighted norm, finite energy, and acceptable localization for a true bound mode. For a resonance, impose outgoing conditions and derive a complex frequency or equivalent leakage rate and lifetime rather than also imposing bound-state normalizability. Identify any rim, collar, sheath, or cavity wall from the solved energy and stress localization rather than imposing a sharp \(\mu_R>0\) domain.

Then solve the nonlinear free-boundary or periodic-orbit problem in which the support field's time-averaged stress, throat shape, throughput, conversion, and reservoir loading are mutually consistent. Harmonic balance, shooting, continuation, or direct time evolution are candidate numerical routes. Only afterward compute the Floquet or linear-response spectrum about that self-consistent supported throat.

Independently derive the source-side environmental gravity field before introducing a later probe. Track the raw \(\widetilde Q_g^{(s)}(r)\), test covariance under \((\mathcal T,\mathcal B_\infty)\mapsto(\mathcal C_{\rm th}\mathcal T,\mathcal R_w\mathcal B_\infty)\), and determine whether one orientation-even but radially signed \(\overline Q_g(r)\) results. Test \(\overline Q_g(r_{\rm local})>0\) separately from reflection parity. The fixed operational calibration and passive response remain canonical ontology tasks and must not be collapsed into this source-side comparison.

For each ordinary supported branch carried into the pair problem, import and verify the ontology's separate requirements of positive-definite quasistatic passive and inertial response and retarded passivity for an unexcited internal reservoir. Common wrong-sign response is not rescued by species universality.

## 11. Nonlinear contact and annihilation program

Linear graph-regime theory must hand off to the full parent fields and core level sets before overhangs, necks, reconnection, pinch-off, or core contact invalidate \(h_\pm\).

### 11.1 One-core solutions

Obtain both orientations of each species from the same parent equations under reflected environmental data, and separately test degeneracy in a reflection-symmetric unchanged environment. Verify the complete candidate map \(\mathcal C_{\rm th}\), including equal orientation-even scalar magnitudes and reversal of every required orientation-odd internal datum. Use the linear transverse mode only as a seed, continue to a self-consistent nonlinear periodic supported throat, and then determine its Floquet or response stability. For a bound seed record its norm and finite energy; for a resonant seed record outgoing conditions, complex frequency or decay rate, and lifetime. Identify diagnostic support regions and full core level sets rather than identifying the interfaces alone with the throat.

### 11.2 Two-core stationary or phase-locked periodic families

Continue same- and opposite-orientation solutions from large projected separation through \(R\sim H_0\) toward \(R=0\). For periodic throats, construct families at explicitly fixed support-mode action, frequency, relative phase, or other derived conjugate data rather than silently treating the averaged geometry as static. Track:

- the graph-regime gap \(d_w^{\rm core}=H_0+2h_t^{\rm tot}-\ell_+^{\rm in}-\ell_-^{\rm in}\), while valid;
- the thresholded full set distance \(d_{\rm core}(\epsilon)\) while two cores remain identifiable, orientation-resolved intersection or component-closure/topology transition, and robustness under reasonable \(\epsilon\) variation after the graph description fails;
- whether orientation-resolved diagnostics or connected-component continuation supplies physical core labels, and the parameter value at which separate labels cease to exist after merger;
- interface topology;
- throat radii and support-mode amplitudes;
- support-mode actions or frequencies, relative phase, phase locking, and any phase-dependent sidebands;
- even and odd field energy;
- signed local net outward flux \(Q_n^{\rm net}\), selected positive throughput \(Q_n\), the reference-cut identity \(\widetilde Q_n^{(w)}=sQ_n^{\rm net}\) or derived curved-cut replacement, and the branch-specific gross order-conversion rate or diagnostic \(Q_\chi\);
- the provisional raw projected gravity source, its reflected-environment covariance, and the orientation-even but radially signed \(\overline Q_g(r)\);
- reflection-symmetric and asymmetric \(\mathcal B_\infty\), including return-halo overlap and any environmental splitting;
- the support-mode norm, energy/stress localization, leakage rate, and lifetime;
- projected positions \(\mathbf x_a\) and the full set distance defined with four-dimensional dummy points \(\mathbf Y_\pm\);
- negative and zero modes.

### 11.3 Time-dependent collision families

Evolve same-species particle–antiparticle collisions over impact parameter, relative velocity, support-mode phase, and internal excitation using the full parent fields. Measure every outgoing channel, including \(P_{u_T},P_{\rm odd},P_{\rm even},P_{u_L},P_{\rm bulk},P_\chi\), and close the material, momentum, energy, and possible entropy ledgers.

### 11.4 Different-species overlap

Repeat the contact calculation for oppositely charged but internally different or composite throats. Determine which additional invariants prevent or modify annihilation and whether a stable neutral composite exists.

## 12. Observable consequences and failure conditions

The geometry can be useful only if its consequences are acceptable.

The two-interface construction fails if any of the following occurs:

1. no stationary drain branch has \(Q_n^{\rm net}>0\), or independently solved reflection partners fail to have equal positive throughput under reflection-symmetric environmental data;
2. throat and environmental data fail to transform covariantly, or the ordinary vacuum is insufficiently reflection symmetric to avoid unacceptable orientation-dependent gravitational splitting;
3. the reflected solutions fail to produce one orientation-even signed \(\overline Q_g(r)\), or the ordinary local branch fails the separate attractive-sign test;
4. no complete odd/even decomposition of the solved throat source exists; an equilibrium branch has no healthy odd Hessian eigenbranch with nonzero throat coupling and the required \(k^2\) static stiffness; or additional long-range odd/even modes produce unacceptable forces or radiation;
5. every primarily thickness-like even branch remains gapless or produces unacceptable long-range species-dependent forces when the intended construction requires screening;
6. no normalizable bound seed or acceptably long-lived outgoing resonance exists, or no seed continues to a self-consistent nonlinear periodic supported throat with positive throughput and acceptable Floquet/response stability;
7. the selected conservative, relaxational, or mixed order branch cannot close material, momentum, energy, inertia/drag, and required entropy/internal-reservoir ledgers;
8. no complete \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\) maps a solved species branch into an equilibrium equal-energy partner, a conservative periodic partner with consistently mapped action/frequency/phase data and cycle-averaged invariants, or a driven partner with equal reflected stationary throughput, dissipation, entropy production, lifetime, stress/momentum fluxes, and response kernels as applicable;
9. the interface graph description loses validity without a consistent handoff from projected \(\mathbf x_a\) to full \(\mathbf X=(\mathbf x,w)\), parent fields, physically labeled localized threshold sets, and the ambient four-dimensional set distance using \(\mathbf Y_\pm\);
10. no well-defined graph-regime gap, compact or suitably closed core diagnostic, orientation-resolved or connected-component identification procedure, threshold-robust contact criterion, or rule for retiring separate labels after merger exists;
11. projected coincidence necessarily produces an unphysical singularity, arbitrary branch selection, loss of conservation, or unacceptable short-distance scattering or bound-state behavior;
12. same-species opposite-orientation pairs have no conservation-respecting annihilation or neutralization channel;
13. different-species opposite charges are forced to annihilate solely because their leading electric monopoles cancel despite unmatched odd or even internal data;
14. the required odd electric branch produces unacceptable extra radiation or drag, including unsuppressed radiation incompatible with the intended light/electromagnetic phenomenology;
15. a conservative periodic pair is assigned a static potential without a derived cycle-averaged ensemble or stress force; a driven dissipative force is assigned an unjustified potential; the odd coefficient \(C_{E,ab}\) hides even/core forces or non-subleading derivative corrections; or a desired flux, gravity, support, spectrum, or short-range result requires a separately tuned coefficient or absolute-value convention;
16. the ordinary supported branch has a non-positive passive or inertial response direction after conventions are fixed, or its retarded response violates passivity for an unexcited internal reservoir.

Finite-thickness corrections, any neutral remnant or hard core, leakage into transverse light, and every additional radiative branch are therefore calculable predictions rather than adjustable exceptions.

## 13. Compact interpretation

The finite brane has two interfaces in the canonical flat Euclidean parent space. Charge orientation can select which side carries the dominant protrusion, while a brane observer records only \(\mathbf x_a\); full core points are \(\mathbf X=(\mathbf x,w)\). Opposite orientations may therefore overlap in projection while remaining separated across the slab. For a complete same-species candidate pair, the cycle-averaged or correctly phase-mapped full odd source vector can cancel even though even source data, thickness deformation, and internal support structures remain; arbitrary periodic relative phase need not cancel instantaneously. The corrected graph-regime separation uses the locally deformed thickness \(H_0+2h_t^{\rm tot}\), not the undeformed \(H_0\) alone. Nonlinear cores are compact or suitably closed localized diagnostic level sets identified through orientation-resolved fields or connected-component continuation. Their intersection or component-topology transition defines contact while two labels remain physical; \(d_{\rm core}(\epsilon)\) is a threshold-robust ambient Euclidean precontact measure using full points \(\mathbf Y_\pm\), not a post-merger decomposition.

For the same solved internal throat branch, \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\) defines the antiparticle candidate, so electron-like plus consistently transformed electron-like means electron–positron, not electron–proton. Whether \(\mathcal I_{\rm internal}\) is trivial, whether the partner tests close, and what happens in collision remain to be derived. A same-species opposite pair may annihilate or form a temporary neutral state. An oppositely charged object of a different or composite species can remain materially distinct after electric cancellation.

The logical separation is

\[
\text{odd electric cancellation}
\not\Rightarrow
\text{even material cancellation}
\not\Rightarrow
\text{core contact}
\not\Rightarrow
\text{annihilation}.
\]

Linear two-interface theory must determine the complete odd/even source sectors, equilibrium Hessian and retarded-response blocks, the selected eigenbranches' static stiffnesses, the multicomponent finite-thickness Green matrices, the support seed and nonlinear periodic throat, reflected-environment gravity behavior, and every radiation channel. Equilibrium kernels may define an ensemble-qualified potential; a conservative periodic pair requires a derived cycle-averaged action/frequency/relative-phase ensemble or averaged stress force; driven dissipative forces require stress, momentum-flux, conversion, and reservoir accounting unless conservativity is derived. The coefficient \(C_{E,ab}\) covers only the odd Coulomb monopole channel. Parity does not by itself determine range, and gaplessness does not by itself produce a Coulomb kernel. Electric neutrality for different species need not cancel unmatched odd multipoles or even material data. The interface variables apply only while the surfaces remain graphs; full parent-field, nonlinear, time-dependent theory must determine contact, reconnection, annihilation, exclusion, and neutral composite formation.
