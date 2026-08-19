# Guided Light and the Photon-Soliton Research Program

## Status and purpose

This document records the current light ontology, the guided-mode hypothesis, the proposed route to a self-bound classical photon analog, and the calculations required to test it. It supplements, rather than replaces, the canonical [ontology and closure ledger](toy_model_ontology_summary.md) and the focused [opposite-orientation throat plan](opposite_orientation_throat_coincidence.md).

It is a **research and derivation plan**, not a claim that the required guided modes, nonlinear couplings, or photon-like solutions have already been obtained.

The model is a classical toy analog. The medium and its constitutive laws may be postulated so that they possess the properties required by the force and light sectors. The test is not whether those properties emerge accidentally from the few equations already written. The test is whether one frozen, mathematically coherent medium can support all required sectors at once without independent retuning.

The following epistemic labels apply throughout:

- **Postulated:** part of the candidate medium's definition.
- **Operational definition:** a quantity defined by how a brane observer would measure it.
- **Target relation:** a result the selected medium must produce; it may not be inserted as an answer.
- **Derived under stated assumptions:** follows from a specified action, approximation, branch, or boundary condition.
- **Open hypothesis:** a plausible physical route that remains to be tested.
- **Failure criterion:** a result that rejects the selected light branch or, where stated, the broader candidate architecture.

The conceptual freeze rule remains in force:

> A constitutive term or coefficient may be postulated as part of the common medium, but once selected it must be propagated through the stable-brane spectrum, far-field forces, photon sector, throat support, radiation, and every dependent calculation. It may not be adjusted only for the photon calculation and then ignored elsewhere.

## Central proposal

The leading light hypothesis is:

> **The photon analog is a gapless, two-polarization, transverse guided carrier of the ordered finite-thickness brane. Its profile is confined across the brane's normal thickness, and its finite three-dimensional envelope is maintained by finite-thickness dispersion together with a reflection-even, nonradiative longitudinal/thickness/order dressing that is slaved to the transverse intensity.**

This hypothesis uses only structures already present in the ontology:

- the ordered finite-thickness brane;
- the globally defined normal axis associated with the extra spatial direction;
- transverse brane shear;
- the in-plane longitudinal mode;
- the two interface variables and their even thickness combination;
- density and material-order response;
- the de-structured bulk as the surrounding phase of the same medium.

It does **not** introduce a separate photon substance or a fundamental electromagnetic field.

The intended causal picture is

```text
ordered finite-thickness brane
        │
        ├── supports a lowest gapless guided transverse branch
        │       │
        │       ├── two independent transverse polarizations
        │       ├── bound profile across the normal direction w
        │       └── net energy transport along the brane
        │
        ├── finite thickness supplies dispersion and a material length scale
        │
        └── finite transverse intensity excites an even material response
                │
                ├── longitudinal compression or rarefaction
                ├── local thickness change
                ├── density change
                └── order-strength change
                         │
                         └── self-created moving guide for the transverse packet
```

The final object would therefore be predominantly transverse but materially dressed:

```text
observable transverse light carrier
        +
bound, determined longitudinal/thickness/order response
        =
classical photon-packet candidate
```

## 1. Scope of the photon claim

### 1.1 Linear light is not yet a photon

The current quadratic brane action conditionally supports two transverse shear modes with speed

\[
c_\gamma^2=\frac{\mu_R}{\rho_{\rm br}},
\]

where \(\mu_R\) is the effective brane shear modulus and \(\rho_{\rm br}\) is its effective inertia density.

That establishes a candidate **classical light-wave sector** under the supplied assumptions. It does not yet establish a localized photon object.

A linear wave equation can support plane waves and ordinary wave packets. It does not normally select:

- one localized unit of energy;
- a preferred packet width;
- a preferred amplitude;
- resistance to three-dimensional spreading;
- stable collisions;
- a discrete energy-frequency law.

The photon program is therefore a nonlinear extension of the background light sector, not a restatement of the existing two-mode count.

### 1.2 Derived homogeneous transverse–longitudinal separation

Under the selected homogeneous, isotropic, parity-even quadratic brane response in three brane dimensions,

\[
\mathcal M(\mathbf k)
=
\mu_Rk^2P_T
+
B_{\rm comp}k^2P_L,
\]

where

\[
P_T
=
I-\frac{\mathbf k\mathbf k^T}{k^2},
\qquad
P_L
=
\frac{\mathbf k\mathbf k^T}{k^2}.
\]

The audited spectrum is

\[
\omega_T^2
=
\frac{\mu_R}{\rho_{\rm br}}k^2,
\qquad
\text{multiplicity }2,
\]

and

\[
\omega_L^2
=
\frac{B_{\rm comp}}{\rho_{\rm br}}k^2,
\qquad
\text{multiplicity }1.
\]

The characteristic polynomial is cubic in \(\omega^2\), with nonzero leading coefficient proportional to \(-\rho_{\rm br}^3\). For \(\rho_{\rm br}\ne0\), the selected homogeneous system therefore has exactly three finite roots counted with multiplicity. Degeneracy can merge the known transverse and longitudinal eigenvalues but cannot create a hidden additional finite mode. At \(B_{\rm comp}=\mu_R\), all three roots coincide; the eigenvalue alone then does not distinguish transverse from longitudinal character, but no new mode appears.

At generic coefficients,

\[
\boxed{
P_T\mathcal M P_L=0.
}
\]

The transverse light root is independent of \(B_{\rm comp}\), while the longitudinal root is independent of \(\mu_R\). Equivalently:

> The homogeneous \(D=3\) quadratic transverse and longitudinal eigenbranches have zero linear cross-block under the selected action.

This is a **derived result under the stated homogeneous quadratic assumptions**. It does not establish decoupling at nonlinear order, on a nonuniform slab, near interfaces or defects, in a parity- or time-reversal-breaking background, or after additional reservoir variables are introduced.

### 1.3 Wave packet, solitary wave, soliton, and photon analog

These terms should remain distinct.

A **wave packet** is a localized superposition of linear waves. It may spread.

A **solitary wave** travels while retaining its shape for a long time because nonlinearity, dispersion, geometry, or another response suppresses spreading.

A **soliton** is a stronger object: a dynamically stable solitary wave that survives interactions in a recognizable form, often with only phase or position shifts.

A **classical photon analog** in this model should satisfy more than temporary localization. At minimum it should possess:

- finite localized energy and momentum;
- exactly two freely selectable asymptotic polarization degrees of freedom;
- a dominant transverse carrier;
- propagation at or extremely near the light speed \(c_\gamma\);
- no freely adjustable rest configuration;
- negligible unacceptable leakage into bulk, longitudinal, electric-odd, thickness-even, or conversion modes;
- a stable or sufficiently long-lived three-dimensional envelope;
- a meaningful frequency, wavelength, energy, and momentum relation;
- sufficiently weak or soliton-like interaction with other light packets.

A classical soliton does not automatically supply quantum mechanics. Even a successful packet may form a continuous family in amplitude, width, or action. The later emergence of a universal action scale, an \(E\)-versus-\(\omega\) quantum law, discrete emission and absorption, photon number, and statistics would remain separate closure problems.

### 1.4 How serious is failure?

The linear light sector could remain mathematically valid even if no free self-bound packet exists. A throat could also trap an otherwise non-self-bound wave as a cavity or defect mode.

However, if the intended ontology is required to explain a freely propagating photon as a localized classical object, then failure to find **any** viable self-bound packet branch is a major photon-ontology failure. It would leave the model with an analog of classical radiation but no analog of an individual free photon.

For that reason, this program should be treated as an early, high-priority falsification track rather than postponed until after the full throat solve.

## 2. The normal direction and the geometry of light

### 2.1 Signed normal versus unsigned normal axis

The finite brane is assumed to be globally two-sided and co-orientable. Let \(N^A\) denote a unit normal to the brane in the parent four-dimensional space.

The **sign** of the normal distinguishes \(+w\) from \(-w\) and is used by the throat sector to represent electric orientation.

Light should ordinarily depend only on the **normal axis**, not on which end is called positive. Its leading constitutive law should therefore depend on reflection-even combinations such as

\[
N^A N^B,
\]

which are invariant under

\[
N^A\rightarrow -N^A.
\]

This gives two related but distinct uses of the same geometry:

\[
\boxed{
\begin{aligned}
N^A
&:\ \text{signed orientation available to the charge sector},\\
N^A N^B
&:\ \text{unsigned normal axis available to the light-guiding sector}.
\end{aligned}
}
\]

The leading free-light equations should therefore be symmetric under \(w\to-w\). A photon should not acquire electric charge merely because its transverse mode is guided by the normal structure of the slab.

### 2.2 Tangential projector

The normal defines a projector onto directions tangent to the brane:

\[
P^{AB}=\delta^{AB}-N^A N^B.
\]

A guided light displacement should be tangent:

\[
N_A u_T^A=0.
\]

For a wave with brane-tangent propagation vector \(k^A\),

\[
N_Ak^A=0,
\]

transversality also requires

\[
k_Au_T^A=0.
\]

The parent medium has four spatial directions. Removing the one normal direction and the one propagation direction leaves two independent polarization directions:

\[
4-1-1=2.
\]

This gives a geometric interpretation of the two-polarization target:

> A light displacement is forbidden from pointing through the normal axis and is transverse to its direction of travel inside the brane, leaving two independent tangential polarization directions.

This count remains conditional on the selected action and constraints. It must be verified by the complete guided-mode spectrum rather than accepted from counting alone.

### 2.3 Transverse isotropy

The most natural material design is **transversely isotropic**:

- no preferred direction among \(x,y,z\) in a homogeneous brane;
- materially different response along the normal axis \(w\).

This permits strong normal structure without selecting a fixed compass direction inside observed space.

A permanently preferred direction within the brane could produce unacceptable anisotropy. By contrast, a unique normal axis is already part of the brane ontology.

## 3. What it means for the light disturbance to “want to spread in w”

The phrase should not mean that the photon's net energy propagation is primarily into the bulk. A photon observed in the brane must carry energy and momentum along the brane.

The more precise proposal is:

> The microscopic shear deformation associated with light naturally develops a strong profile, relaxation tendency, or leakage tendency across the normal axis, while the packet's group velocity and net energy transport remain tangent to the brane. The finite ordered slab and the different bulk phase constrain that normal structure into a guided mode.

A guided wave has two different directions built into it:

- **normal mode structure:** how the field varies across the guide;
- **propagation direction:** where energy travels along the guide.

The light field may therefore have a strongly structured or nearly standing profile in \(w\) while translating along a brane direction such as \(z\).

A schematic mode is

\[
\mathbf u_T(\mathbf x,w,t)
=
\operatorname{Re}\!\left[
\left(
\boldsymbol\epsilon_1 A_1(\mathbf x,t)
+
\boldsymbol\epsilon_2 A_2(\mathbf x,t)
\right)
f_0(w)
e^{i(\mathbf k\cdot\mathbf x-\omega t)}
\right],
\]

where

- \(\boldsymbol\epsilon_1,\boldsymbol\epsilon_2\) are the two tangential polarizations;
- \(f_0(w)\) is the lowest guided profile across the slab;
- \(A_1,A_2\) are slowly varying packet envelopes;
- \(\mathbf k\) is tangent to the brane.

The profile \(f_0(w)\) is part of the photon solution even though a brane observer may only resolve the projected three-dimensional envelope.

## 4. Finite slab as a waveguide

### 4.1 Why the ordered slab can guide shear

The ordered brane supports shear. The de-structured bulk is not assumed to support the same resolved transverse elastic branch.

A minimal variable-coefficient transverse action may have the form

\[
S_T^{(2)}
=
\frac12\int dt\,d^3x\,dw
\left[
\rho_T(w)|\dot{\mathbf u}_T|^2
-
\mu_\parallel(w)|\nabla_\parallel\mathbf u_T|^2
-
\mu_w(w)|\partial_w\mathbf u_T|^2
-
\mathbf u_T\cdot\mathbf V_T(w)\mathbf u_T
\right],
\]

where the coefficient profiles descend from the solved background \(n_0(w),\chi_0(w)\) and the selected shear-supporting constitutive structure.

This expression is a template, not a frozen final action. The actual reduction may include:

- mixed derivatives;
- nonlocal kernels;
- higher-gradient terms;
- coupling to density, order, thickness, or longitudinal fields;
- branch-specific damping and internal reservoirs.

The normal eigenproblem has schematic form

\[
\mathcal L_w f_n(w)
=
\lambda_n\,\mathcal W_w f_n(w).
\]

Projection onto a normal branch then produces

\[
\omega_n^2(k)
=
\lambda_n+c_n^2k^2+a_{4,n}k^4+a_{6,n}k^6+\cdots.
\]

### 4.2 Gapless lowest branch

The photon branch must be gapless in the intended long-wavelength regime:

\[
\boxed{
\lambda_0=0,
\qquad
\omega_0(k)\rightarrow0
\quad\text{as}\quad k\rightarrow0.
}
\]

Higher normal branches may be gapped:

\[
\lambda_{n>0}>0.
\]

A naive thickness standing wave with a fixed normal frequency would instead give

\[
\omega^2=\Omega_w^2+c^2k^2,
\]

which behaves like a massive mode at low \(k\). That is not the desired free-light branch.

The lowest branch must therefore be an acoustic-like guided zero mode, a symmetry-protected mode, a coupled-interface mode, or another branch whose frequency vanishes appropriately even though its profile is confined across \(w\).

### 4.3 Confinement is not guaranteed by “zero bulk shear” alone

It is tempting to say that the mode is automatically trapped because the bulk does not support shear. That is only a design intuition.

The full normal operator must show one of the following:

- exponential or otherwise normalizable decay into the bulk;
- exact confinement by boundary conditions and material order;
- an acceptably long-lived outgoing resonance with a calculated leakage rate;
- another mathematically closed nonpropagating bulk response.

A vanishing bulk shear modulus by itself may create a floppy or strongly coupled region rather than a clean bound state. The spectral problem must decide.

The localization test must not be reduced to normalizing \(f_0\) to unit norm, because that normalization can hide a profile that spreads arbitrarily far into the bulk. At each fixed in-brane wave number \(k\), and within the same complete constrained parity and polarization block as the candidate mode, let \(\omega_{\rm ess}(k)\) denote the lower edge of the relevant essential spectrum and define

\[
\boxed{
\Delta_{\rm spec}(k)
=
\omega_{\rm ess}(k)-\omega_0(k).
}
\]

If the generalized eigenproblem is naturally posed in \(\omega^2\), also track

\[
\Delta_{\rm spec}^{(2)}(k)
=
\omega_{\rm ess}^2(k)-\omega_0^2(k),
\]

so the conclusion does not depend on a square-root branch convention. A true isolated bound mode has positive spectral separation and zero resonance width at fixed \(k\). For an outgoing resonance, write its pole conventionally as

\[
\omega_{\rm pole}(k)
=
\Omega_0(k)-\frac{i}{2}\Gamma_0(k),
\qquad
\Gamma_0(k)\ge0,
\]

and calculate the lifetime and leakage from \(\Gamma_0(k)\).

The calculation must also track a normal localization length \(\ell_w(k)\), obtained from the asymptotic decay or a stated positive energy- or norm-weighted width, and a brane-participation fraction. Let \(b_{\rm br}(w)\) be a smooth window derived from the solved ordered profile, with \(0\le b_{\rm br}\le1\), rather than an arbitrarily selected hard slab boundary. For a healthy bound branch, define

\[
\boxed{
\mathcal P_{\rm br}(k)
=
\frac{
\displaystyle\int_{\mathbb R}dw\,
b_{\rm br}(w)
f_0^*(w;k)\mathcal W_w(k)f_0(w;k)
}{
\displaystyle\int_{\mathbb R}dw\,
f_0^*(w;k)\mathcal W_w(k)f_0(w;k)
}.
}
\]

For a leaky resonance, whose outgoing profile is not ordinarily square-integrable, any analogous participation measure must use a stated finite-window, flux, pole-residue, complex-scaling, or absorbing-boundary prescription and must be stable as that regulator is varied.

To expose delocalization that a unit-norm convention would conceal, fix a nonzero on-brane amplitude—for an even lowest mode with nonzero mid-surface value, for example,

\[
f_0(0;k)=1,
\]

and track the resulting effective three-dimensional inertia

\[
\boxed{
\mathcal I_{\rm 3D}(k)
=
\int_{\mathbb R}dw\,
f_0^*(w;k)\mathcal W_w(k)f_0(w;k).
}
\]

If the mode has a mid-surface node, the normalization must instead fix another nonzero on-brane field or energy diagnostic. Once an explicit brane source or observable has been defined, the associated pole residue should also be tracked.

No universal requirement is imposed that \(\Delta_{\rm spec}\) remain bounded away from zero or that \(\ell_w\) approach a finite constant as \(k\to0\). A gapless guided or surface branch may have a closing separation or a penetration depth that scales with wavelength. The required result is an acceptable derived asymptotic regime: the leakage ratio or lifetime remains acceptable, \(\mathcal P_{\rm br}\) does not vanish unacceptably, \(\mathcal I_{\rm 3D}\) remains finite and nonzero under the fixed on-brane normalization, \(\ell_w\) remains compatible with an effective three-dimensional description, and the brane-observable pole residue does not disappear.

### 4.4 Finite thickness supplies a length scale and dispersion

The slab thickness \(H_{\rm br}\) introduces a natural material length scale. After normal-mode projection, the effective dispersion may contain

\[
\omega_T^2(k)
=
c_\gamma^2k^2
+a_4H_{\rm br}^2k^4
+a_6H_{\rm br}^4k^6
+\cdots.
\]

These corrections matter because a self-bound packet generally requires a balance between spreading and nonlinear response.

The signs and magnitudes of \(a_4,a_6,\ldots\) are outputs of the slab reduction. They may be useful, negligible, destabilizing, or incompatible with other sectors.

## 5. The canal-soliton analogy: useful and limited

The shallow-channel solitary-wave analogy captures three relevant ideas:

\[
\boxed{
\text{restricted geometry}
+
\text{dispersion}
+
\text{nonlinearity}
}
\]

A finite-depth channel introduces dispersion. Nonlinear wave speed changes can oppose that spreading. Side walls suppress lateral expansion and make the problem nearly one-dimensional.

The brane analogy is:

\[
\begin{aligned}
\text{channel depth}
&\leftrightarrow H_{\rm br},\\
\text{finite-depth dispersion}
&\leftrightarrow\text{guided-mode dispersion from slab thickness},\\
\text{water nonlinearity}
&\leftrightarrow\text{nonlinear shear/order/material response}.
\end{aligned}
\]

The limitation is equally important:

\[
\boxed{
\text{confinement across }w
\neq
\text{localization in all three brane directions}.
}
\]

The slab supplies boundaries across its thickness, but no fixed side walls surround a photon in the two brane directions transverse to its motion.

For a packet moving along \(z\), ordinary slab guidance can prevent loss into \(w\) while the envelope still spreads in \(x\), \(y\), or along \(z\).

The missing analog of the canal walls may therefore need to be a **self-induced channel carried by the photon itself**.

## 6. Linear directional spreading and the role of substructure

### 6.1 Normal-dominant spreading is a permissible constitutive target

A guided envelope before normal projection might contain

\[
i(\partial_t+v_g\partial_z)A
+
D_z\partial_z^2A
+
D_\perp(\partial_x^2+\partial_y^2)A
+
D_w\partial_w^2A
+\cdots
=0.
\]

The proposed substructure could favor

\[
|D_w|\gg |D_\perp|,
\]

so most of the mode's natural transverse variation lies across \(w\), where the finite slab constrains it.

This is physically allowable because the brane already distinguishes normal and tangential response. The medium need not spread a disturbance isotropically through all four spatial directions.

### 6.2 Why exact linear self-collimation is not automatic

Inside a homogeneous and isotropic brane, an ordinary dispersion relation depends only on \(|\mathbf k|\):

\[
\omega=F(|\mathbf k|).
\]

A packet centered on \(k_0\hat{\mathbf z}\) generally has nonzero curvature of \(\omega\) in the transverse momentum directions. That curvature produces diffraction.

Therefore, exact suppression of in-brane spreading for every propagation direction usually requires additional structure beyond a generic local isotropic linear wave:

- special anisotropic or self-collimating dispersion;
- nonlocal response;
- coupled modes whose diffraction cancels;
- or a nonlinear self-created guide that follows the packet.

A permanent microscopic lattice with preferred directions inside the brane might create self-collimation but would also risk observable spatial anisotropy. A packet-created guide is more naturally compatible with rotational symmetry because it aligns itself with the packet's own direction of travel.

### 6.3 Linear self-collimation may still help

The selected constitutive law may be designed so that the lowest guided branch has unusually small transverse diffraction over a broad wave-number range.

This would reduce the burden placed on nonlinear binding. It need not provide perfect localization by itself.

A useful target is therefore not necessarily

\[
D_\perp=0,
\]

but rather

\[
|D_\perp|\ll \text{ordinary elastic scale}
\]

through the intended photon frequency range, followed by nonlinear cancellation of the residual spreading.

## 7. The proposed photon structure

### 7.1 Transverse carrier with localized envelope

The most photon-like candidate is an oscillatory guided carrier inside a localized envelope:

\[
\boxed{
\mathbf u_T
=
\operatorname{Re}\!\left[
\left(
\boldsymbol\epsilon_1A_1
+
\boldsymbol\epsilon_2A_2
\right)
f_0(w)
e^{i(k_0z-\omega_0t)}
\right].
}
\]

This separates several roles:

- \(e^{i(k_0z-\omega_0t)}\): carrier frequency and wavelength;
- \(f_0(w)\): guided normal profile;
- \(A_1,A_2\): localized three-dimensional envelope;
- \(\boldsymbol\epsilon_1,\boldsymbol\epsilon_2\): two polarization directions.

This route is preferred over a single nonoscillating hump because it naturally retains:

- frequency;
- wavelength;
- phase;
- polarization;
- a moving localized energy packet.

### 7.2 What should be localized

The displacement field itself need not have compact support. A constant displacement can carry no local strain energy.

The physically important localized quantities are:

- shear strain;
- material velocity;
- stress;
- energy density;
- momentum density;
- the helper-field disturbance.

A finite-energy target has schematic form

\[
E_\gamma
=
\int d^3x\,dw
\left[
\frac12\rho_T|\dot{\mathbf u}_T|^2
+
\frac12\mu|\nabla\mathbf u_T|^2
+
E_{\rm nl}
+
E_{\rm disp}
+
E_{\rm helper}
\right]
<\infty.
\]

The energy and stress should remain concentrated even if the displacement possesses a harmless offset or weak tail.

### 7.3 No independent rest state

A photon-like branch should not generically admit a stationary localized lump with the same conserved data. Its natural solutions should travel at a speed locked near the transverse light branch.

The detailed requirement may be expressed through a dispersion or energy-momentum relation rather than imposed as a primitive rule. A target low-energy relation is approximately

\[
E_\gamma\simeq c_\gamma |\mathbf P_\gamma|
\]

for the freely propagating branch, with calculable corrections.

### 7.4 Route A and Route B

Two nonlinear reductions should be retained, but they have different roles.

**Route A — KdV-like solitary shear pulse.** Finite-thickness dispersion and amplitude-dependent shear propagation may produce an effectively one-dimensional equation of the form

\[
\partial_t A
+c_\gamma\partial_z A
+\alpha A\partial_z A
+\beta H_{\rm br}^2\partial_z^3 A
=0.
\]

A localized traveling hump would be a useful demonstration that the brane's nonlinearity and dispersion can balance. It is the closest direct analog of a shallow-channel solitary wave. It is not the preferred final photon model because it is naturally one-dimensional, generally has amplitude-dependent speed, does not by itself solve lateral localization, and does not naturally encode an oscillatory carrier with two polarization states.

**Route B — guided transverse carrier with a localized vector envelope.** The preferred photon candidate is the oscillatory two-polarization carrier described above, with its envelope bound by finite-thickness dispersion and a slaved material response. Route B retains frequency, phase, wavelength, and polarization and can in principle produce a fully three-dimensional packet.

Route A should therefore be used as an inexpensive precursor and diagnostic. Success would establish a nonlinear-dispersive balance but would not close the photon ontology. Failure would constrain the constitutive signs and make the stronger Route B problem less promising. Route B remains the principal target.

## 8. The stray longitudinal mode as a binding field

### 8.1 The longitudinal mode has not yet been assigned a final role

The canonical ontology already contains an in-plane longitudinal mode \(u_L\). It is not charge and it is not the \(\pm w\) throat orientation.

Its physical dispersion, confinement, leakage, and coupling to transverse waves remain open.

This makes it a legitimate candidate participant in photon localization.

### 8.2 Linear decoupling, nonlinear coupling

In the selected homogeneous \(D=3\) weak-field theory, the transverse and longitudinal branches are derived to have zero quadratic cross-block. This preserves two freely selectable transverse light polarizations at leading order. The result is not assumed to survive an arbitrary nonuniform slab or an enlarged quadratic action; Stage L2 must test the complete operator.

At finite amplitude, however, transverse shear can produce compression, density change, order change, or thickness response. Nonlinear coupling is therefore natural even if linear mixing vanishes.

A natural weakly nonlinear ordering is

\[
\boxed{
u_T=O(\epsilon),
\qquad
\eta_{\rm helper}=O(\epsilon^2),
\qquad
\text{feedback on }u_T=O(\epsilon^3).
}
\]

This hierarchy is a reduction target, not an independent postulate. It follows only if the first allowed even helper source is quadratic in the transverse amplitude. A derived constitutive law with a different consistent leading order must state and propagate that ordering.

The desired feedback is

\[
\boxed{
\text{transverse intensity}
\rightarrow
\text{even material compression/order response}
\rightarrow
\text{local change in shear propagation}
\rightarrow
\text{transverse self-guiding}.
}
\]

### 8.3 The helper field is likely a mixed eigenmode

The actual helper should not be assumed to be pure \(u_L\). In a finite slab, the relevant reflection-even material eigenmode may mix

- longitudinal compression \(\nabla\cdot\mathbf u_L\);
- thickness disturbance \(h_t\);
- order variation \(\delta\chi_B\);
- density variation \(\delta n\);
- bulk or reservoir response.

Introduce a schematic helper coordinate

\[
\boxed{
\eta\equiv\eta_{\rm helper}
=
a_L\nabla\cdot\mathbf u_L
+a_t h_t
+a_\chi\delta\chi_B
+a_n\delta n
+\cdots.
}
\]

The complete slab spectrum must derive the coefficients and identify whether one combination produces the required response.

The research question is therefore:

> Which non-transverse, reflection-even eigenmode is sourced by transverse-wave intensity, and does its backreaction self-focus the guided carrier without radiating or destabilizing the brane?

### 8.4 Why the helper does not create a third photon polarization

Let

\[
I=|A_1|^2+|A_2|^2
\]

be the total transverse intensity.

The desired helper is slaved to \(I\):

\[
\eta=\eta[I].
\]

Then \(A_1\) and \(A_2\) remain the two independent polarization amplitudes. The helper amplitude is determined by the packet rather than freely selectable.

Thus the photon candidate contains

\[
\boxed{
\text{two free transverse polarization components}
+
\text{one determined material dressing}.
}
\]

It does not contain three freely selectable propagating polarizations.

The background brane may still support an independent free longitudinal wave. The claim is only that the longitudinal component inside the photon packet is a bound response rather than an independent photon label.

### 8.5 Reflection parity protects electrical neutrality

The transverse intensity \(I\) is reflection even under \(w\to-w\). It can therefore source even responses such as

- \(u_L\) compression;
- \(h_t\) thickness change;
- \(\delta n\);
- \(\delta\chi_B\).

It should not directly source the reflection-odd electric mid-surface channel at leading order on a symmetric brane.

This gives a useful division:

\[
\boxed{
\begin{aligned}
\text{photon intensity}
&\rightarrow\text{even helper dressing},\\
\text{oriented throat}
&\rightarrow\text{odd electric source}.
\end{aligned}
}
\]

A successful photon packet can therefore use thickness and longitudinal structure without acquiring a leading electric monopole.

## 9. Nonradiative versus radiative longitudinal response

### 9.1 Three separate speed comparisons

The homogeneous audit introduces three speeds that must not be collapsed into one inequality:

\[
c_\gamma
:\ \text{low-amplitude transverse guided speed},
\]

\[
c_L^2
=
\frac{B_{\rm comp}}{\rho_{\rm br}}
:\ \text{homogeneous in-plane longitudinal speed},
\]

and

\[
c_{s0}
:\ \text{asymptotic bulk-sound speed in the simple phase-matching model}.
\]

Let \(v_{\rm packet}=|\mathbf v|\) denote the translational speed of the nonlinear photon packet.

The first comparison is

\[
v_{\rm packet}
\quad\text{versus}\quad
c_L.
\]

For a simple acoustic longitudinal branch, \(v_{\rm packet}<c_L\) helps prevent a DC comoving source from resonantly emitting a free in-plane longitudinal wave. It is neither necessary nor sufficient for the complete mixed helper to be nonradiative.

The second comparison is

\[
c_L
\quad\text{versus}\quad
c_{s0}.
\]

It places the **free homogeneous longitudinal eigenbranch** below, at, or inside the simple bulk-sound continuum.

The third comparison is

\[
v_{\rm packet}
\quad\text{versus}\quad
c_{s0}.
\]

It places the DC comoving helper source below, at, or above the corresponding simple bulk-sound grazing condition. Higher harmonics and sidebands require their own frequency-dependent tests.

If a later sector enforces

\[
v_{\rm packet}\simeq c_\gamma=c_{s0},
\]

the DC comoving helper may lie near a bulk grazing threshold for collinear Fourier components. This is not automatically fatal: interface overlap, polarization, finite-thickness dispersion, a gap, or a mixed helper eigenmode may suppress radiation. It is nevertheless a cross-sector compatibility condition that must be calculated rather than omitted.

### 9.2 Free longitudinal branch versus the bulk-sound continuum

For the simple homogeneous phase-matching comparison, let

\[
\omega=c_Lk,
\qquad
c_L^2=\frac{B_{\rm comp}}{\rho_{\rm br}},
\]

and let the isotropic bulk-sound branch satisfy

\[
\omega^2
=
c_{s0}^2(k^2+k_w^2),
\]

where \(k\) is the conserved tangential wave number and \(k_w\) is the bulk normal wave number. Matching the free brane-longitudinal branch to bulk sound gives

\[
\boxed{
k_w^2
=
k^2
\left(
\frac{c_L^2}{c_{s0}^2}-1
\right).
}
\]

The three kinematic regimes are:

| Regime | Normal wave number | Kinematic meaning |
|---|---:|---|
| \(c_L<c_{s0}\) | \(k_w^2<0\) | No propagating bulk-sound channel for the matched free longitudinal mode; an evanescent profile is possible |
| \(c_L=c_{s0}\) | \(k_w=0\) | Grazing or continuum threshold |
| \(c_L>c_{s0}\) | \(k_w^2>0\) | A propagating bulk-sound channel is kinematically open |

The audited `KW_ZERO_LOCUS` is

\[
\boxed{
B_{\rm comp}
=
\rho_{\rm br}c_{s0}^{\,2},
\qquad
c_L=c_{s0}.
}
\]

This phase-matching test determines whether a propagating bulk channel exists in the simple homogeneous comparison. It does not determine whether the brane mode couples to that channel, whether the overlap vanishes by symmetry, whether a bound state in the continuum exists, or whether an evanescent solution is a normalizable bound eigenstate. Those questions require the complete brane–bulk interface operator, boundary conditions, kinetic norm, and spectral classification.

Thus the three regimes must be called **bulk channel kinematically closed**, **grazing threshold**, and **bulk channel kinematically open**—not simply bound, threshold, and radiating.

### 9.3 Forced-helper radiation and harmonic audit

The photon helper is not necessarily an on-shell free longitudinal wave. For a packet traveling with velocity \(\mathbf v\) and possessing comoving internal frequency \(\Omega\), write its nonlinear source spectrum as

\[
\boxed{
\omega_{m,\mathbf q}
=
\mathbf q\cdot\mathbf v+m\Omega,
\qquad
m\in\mathbb Z.
}
\]

The free in-plane longitudinal resonance test is

\[
\omega_L(\mathbf q)
=
\mathbf q\cdot\mathbf v+m\Omega.
\]

For the simple isotropic bulk-sound model, a forced source component has

\[
\boxed{
k_w^2(m,\mathbf q)
=
\frac{
\left(
\mathbf q\cdot\mathbf v+m\Omega
\right)^2
}{
c_{s0}^2
}
-|\mathbf q|^2.
}
\]

Kinematic intersection with a receiving dispersion branch is not enough. Radiation into branch \(j\) additionally requires nonzero projected source overlap:

\[
\boxed{
\Pi_j(\mathbf q)\,
\mathbf S_m(\mathbf q)
\ne0.
}
\]

The complete audit must include every longitudinal, bulk, thickness, order, density, electric, and conversion branch. The KW threshold places the **free longitudinal eigenbranch** relative to the simple bulk-sound continuum. The Fourier-support and projector audit decides whether the **nonlinearly forced helper** emits.

A slaved off-shell helper can remain nonradiative even when the corresponding free eigenmode lies in an open continuum, provided its Fourier support does not hit a receiving pole or its interface overlap is sufficiently suppressed. Conversely, a DC helper may be localized while a second harmonic, sideband, or higher harmonic radiates.

### 9.4 Comoving helper equation

A diagnostic helper equation may be written

\[
\partial_t^2\eta
-c_\eta^2\nabla^2\eta
+\Omega_\eta^2\eta
=\alpha I+\beta\nabla^2 I+\cdots.
\]

In a packet frame \(\xi=z-vt\), the operator changes character depending on the relative dispersion and speed.

A bound helper requires an elliptic, screened, or otherwise nonradiative comoving response over the packet's spectrum. A hyperbolic outgoing solution indicates wake production.

This is a direct calculation and a major early gate.

### 9.5 A gapped helper may be useful

The helper need not itself be gapless. A moderately gapped even mode can generate a localized, finite-range response around the photon envelope.

After being integrated out, it may produce a nonlocal effective interaction

\[
\eta(\mathbf x)
=-\alpha\int d^3x'\,G_\eta(\mathbf x-\mathbf x')I(\mathbf x'),
\]

which feeds back into the transverse equation.

A finite-range or nonlocal response can stabilize a three-dimensional packet more effectively than a simple local cubic self-focusing term.

The same gapped helper must remain compatible with the electric-sector requirement that unwanted even thickness forces be short-ranged.

### 9.6 Current radiation status

| Channel or question | Status |
|---|---|
| Homogeneous quadratic \(T\leftrightarrow L\) conversion | Derived to vanish in \(D=3\) under the selected action |
| Number of finite homogeneous roots | Exactly three under the stated assumptions; no hidden finite mode |
| Free longitudinal branch versus bulk continuum | Simple kinematic threshold derived at \(c_L=c_{s0}\) |
| Actual free-longitudinal leakage | Open; requires interface coupling and spectral boundary conditions |
| Nonlinear DC helper response | Desired if localized and nonradiative; not yet derived |
| Nonlinear \(2\omega_0\), sideband, and higher-harmonic leakage | Open |
| Defect-, curvature-, or slab-gradient-induced conversion | Open |
| Longitudinal-to-bulk leakage in the complete slab | S11b / open |

## 10. Candidate reduced envelope system

The following system is a diagnostic target, not a postulated final equation.

Let \(A_a\), \(a=1,2\), denote the two polarization envelopes and

\[
I=\sum_{a=1}^2|A_a|^2.
\]

A projected guided-mode equation may take the schematic form

\[
\begin{aligned}
i(\partial_t+v_g\partial_z)A_a
&+D_\parallel\partial_z^2A_a
+D_\perp\nabla_\perp^2A_a\\
&+\gamma_{\rm dir}I A_a
+\gamma_\eta\eta A_a
+\mathcal N_a^{\rm pol}(A_1,A_2)
+\cdots
=0,
\end{aligned}
\]

where

- \(D_\parallel\) controls longitudinal envelope dispersion;
- \(D_\perp\) controls in-brane diffraction;
- \(\gamma_{\rm dir}\) is direct transverse nonlinearity;
- \(\gamma_\eta\eta\) is helper-mediated self-guiding;
- \(\mathcal N_a^{\rm pol}\) contains polarization-dependent nonlinearities allowed by symmetry.

The helper equation may take the form

\[
\mathcal L_\eta(\partial_t,\nabla)\eta
=-\alpha I+\cdots.
\]

In a quasistatic nonradiative regime,

\[
\eta=-\alpha\mathcal L_\eta^{-1}I,
\]

and the transverse field experiences a nonlocal effective nonlinearity.

The signs are decisive:

- focusing response can oppose diffraction;
- defocusing response accelerates spreading;
- an unbounded focusing energy can cause collapse;
- a saturating, competing, or nonlocal response may permit stable three-dimensional packets.

No sign should be chosen after the fact. Once a common constitutive law is selected, these coefficients are calculated projections of that law.

## 11. Minimal constitutive requirements

The photon branch may postulate a medium with the following shared properties.

### 11.1 Ordered-state shear support

The ordered phase must possess positive, bounded transverse shear energy and the required two transverse polarizations.

### 11.2 Normal-tangential anisotropy

The homogeneous brane may be isotropic within \(x,y,z\) while having distinct coefficients across \(w\).

### 11.3 Guided gapless branch

The full variable-coefficient normal problem must contain a lowest gapless branch with a normalizable or acceptably long-lived profile.

### 11.4 Finite-thickness dispersion

The guided branch must possess enough higher-order dispersion or nonlocality to participate in envelope stabilization without producing unacceptable observed vacuum dispersion.

### 11.5 Nonlinear transverse-material coupling

Finite transverse intensity must couple to one or more reflection-even material variables.

Candidate constitutive invariants include schematic terms such as

\[
(\nabla\cdot\mathbf u_L)\,\mathcal S_T^2,
\qquad
h_t\mathcal S_T^2,
\qquad
\delta\chi_B\mathcal S_T^2,
\qquad
\delta n\mathcal S_T^2,
\]

where \(\mathcal S_T\) denotes an appropriate transverse shear strain invariant.

### 11.6 Bounded energy

The nonlinear energy must be bounded below on the physical constrained space. A self-focusing term that permits catastrophic collapse is not sufficient.

Possible stabilizers may include, if selected as common constitutive structure:

- saturating order response;
- competing higher-order nonlinearities;
- positive higher-gradient terms;
- nonlocal helper response;
- finite-thickness mode coupling.

### 11.7 Shared use across sectors

Every selected coefficient must also enter, where relevant:

- the stable slab;
- the longitudinal spectrum;
- electric even-sector response;
- throat support-mode spectrum;
- moving-throat radiation;
- drag thresholds;
- force-sector Green matrices.

No photon-only hidden material law is permitted.

## 12. Research path independent of the throat

This program can proceed in parallel with the far-field force audit because its primary inputs are properties of the throat-free slab.

The shared upstream problem is

```text
parent constitutive law
        ↓
stable finite slab profile
        ↓
complete background spectrum
```

The downstream tracks then separate:

```text
                    stable slab and spectrum
                             │
              ┌──────────────┴──────────────┐
              │                             │
      low-k far-field audit        finite-k nonlinear light audit
              │                             │
  force carriers and ranges       guided photon-packet existence
```

The force track emphasizes

- \(k\to0\);
- static Hessians and Green matrices;
- long-range tails;
- screening;
- source parity;
- characteristic speeds.

The photon track emphasizes

- finite carrier wave number;
- normal guided profiles;
- dispersion curvature;
- nonlinear projection coefficients;
- self-localization;
- stability and collisions.

Both must use the same frozen slab and constitutive coefficients.

## 13. Detailed derivation program

### 13.1 Stage L0 — freeze the light capability contract

Before solving, define the exact output required from the medium.

Freeze the following notation and distinctions:

- \(c_\gamma\): low-amplitude transverse guided speed;
- \(c_L\): homogeneous in-plane longitudinal speed;
- \(c_{s0}\): asymptotic bulk-sound speed entering the simple phase-matching comparison;
- \(v_{\rm packet}\): translational velocity of the nonlinear photon packet;
- \(k_w\): bulk normal wave number;
- `KW_ZERO_LOCUS`: the simple homogeneous grazing locus \(c_L=c_{s0}\);
- \(\omega_0\): carrier frequency in the selected frame;
- \(\Omega\): comoving internal periodic frequency used to organize source harmonics and sidebands;
- free-mode continuum placement versus forced-source radiation as distinct questions.

At minimum:

- one lowest gapless guided transverse branch;
- two degenerate tangential polarizations on the homogeneous brane;
- a profile confined across \(w\);
- a sufficiently nondispersive observational regime;
- calculable higher-order dispersion outside that regime;
- a candidate nonlinear self-guiding mechanism;
- no leading electric-odd source from a free photon;
- acceptable coupling to longitudinal, thickness, density, order, and bulk modes;
- an explicit audit of DC, harmonic, and sideband radiation channels.

This document is the first version of that contract.

### 13.2 Stage L1 — solve the stable slab

Obtain a stationary background

\[
n_0(w),\qquad \chi_0(w),
\]

with selected finite thickness \(H_{\rm br}\).

Verify:

- reflection symmetry for the ordinary vacuum branch;
- positive or sufficiently metastable interface energy;
- bounded energy;
- the origin and sign of shear support;
- the branch choice for order dynamics;
- consistent bulk boundary conditions.

No guided-light conclusion is trustworthy without a stable background.

### 13.3 Stage L2 — derive the complete linear guided spectrum

Linearize every relevant field about the slab:

\[
\delta\boldsymbol\Phi
=
(\mathbf u_T,u_L,h_m,h_t,\delta n,\delta\theta,\delta\chi_B,\ldots)^T.
\]

Solve the normal eigenproblem for each brane wave vector \(\mathbf k\).

Determine:

- whether the complete \(w\)-dependent operator reproduces the exact homogeneous \(D=3\) two-transverse/one-longitudinal roots in the appropriate limit;
- whether the zero homogeneous quadratic transverse–longitudinal cross-block survives or is replaced by slab-gradient, interface, anisotropic, parity-breaking, or reservoir-induced mixing;
- the number and parity of guided branches;
- the normal profiles \(f_{a,n}(w;k)\);
- the dispersions \(\omega_{a,n}(k)\);
- energy signs and norm signs;
- group velocities;
- the relevant essential-spectrum edges and \(\Delta_{\rm spec}(k)\), or \(\Delta_{\rm spec}^{(2)}(k)\) where natural;
- bulk leakage and resonance widths \(\Gamma_0(k)\);
- normal localization lengths \(\ell_w(k)\);
- brane-participation fractions \(\mathcal P_{\rm br}(k)\);
- fixed-on-brane-normalization effective inertia \(\mathcal I_{\rm 3D}(k)\), and pole residues once a brane observable is specified;
- the bulk essential spectrum and the placement of the free longitudinal branch relative to it;
- whether the complete slab reproduces or replaces the simple `KW_ZERO_LOCUS`;
- whether the longitudinal branch is a bound eigenstate, threshold state, bound state in the continuum, or outgoing resonance;
- interface overlaps and leakage widths for every kinematically open channel;
- whether the \(k_w=0\) locus produces long normal tails, singular normalization, anomalous density of states, or enhanced mode conversion;
- transverse-longitudinal mixing;
- interface and order participation;
- whether the lowest transverse branch is truly gapless;
- whether the two polarizations remain degenerate.

The homogeneous KW formula is a baseline and consistency limit, not a substitute for the complete variable-coefficient operator.

### 13.4 Stage L3 — construct the projected transverse effective action

Project the parent action or response operator onto the lowest transverse branch.

Record:

- \(c_\gamma\);
- \(a_4,a_6,\ldots\);
- frequency dependence;
- nonlocal kernels;
- damping if the selected branch is mixed or relaxational;
- the range over which the mode is approximately linear and nondispersive.

The projected action must reproduce the conditional quadratic light sector already present in the ontology.

### 13.5 Stage L4 — derive nonlinear couplings

Expand the same parent constitutive law to the first nontrivial nonlinear order.

Project all allowed terms onto:

- the two lowest transverse polarization amplitudes;
- each candidate even helper mode;
- any potentially dangerous odd, bulk, or conversion mode.

Derive coefficients for:

\[
I^2,
\qquad
\eta I,
\qquad
h_tI,
\qquad
\delta\chi_BI,
\qquad
\delta nI,
\]

and polarization-dependent invariants.

Promote the cubic-invariant audit explicitly. Among the candidate terms to derive or exclude are

\[
(\nabla\cdot\mathbf u_L)\mathcal S_T^2,
\qquad
h_t\mathcal S_T^2,
\qquad
\delta n\,\mathcal S_T^2,
\qquad
\delta\chi_B\,\mathcal S_T^2.
\]

For every candidate coupling, determine:

- whether it is independent or a total derivative;
- its reflection and time-reversal parity;
- whether it preserves the homogeneous quadratic block separation;
- its DC, \(2\omega_0\), sideband, and higher-harmonic source content;
- its projector and interface overlap with free longitudinal, bulk, thickness, order, density, electric, and conversion continua;
- its focusing or defocusing sign;
- whether the complete constrained energy remains bounded below.

Determine:

- focusing or defocusing sign;
- saturation;
- energy boundedness;
- whether circular and linear polarizations remain equivalent;
- whether the intensity sources the electric-odd block;
- whether the helper is local, nonlocal, inertial, or relaxational.

### 13.6 Stage L5 — derive the weakly nonlinear envelope system

Use a multiple-scale expansion around a transverse carrier \((k_0,\omega_0)\).

Unless the derived constitutive law produces a different consistent ordering, use

\[
A=O(\epsilon),
\qquad
\eta=O(\epsilon^2),
\]

so the helper first responds to transverse intensity and feeds back on the carrier at \(O(\epsilon^3)\).

The output should be a **derived** vector envelope/helper system, not an assumed nonlinear Schrödinger equation.

The reduction must state:

- scale separation;
- carrier frequency range;
- amplitude ordering;
- retained helper modes;
- neglected branches;
- conservation or dissipation law;
- validity time and length scales.

### 13.7 Stage L6 — perform inexpensive no-go and existence tests

Before large simulations, test:

- energy boundedness;
- virial or scaling identities;
- collapse criteria;
- transverse instability of one-dimensional pulses;
- free longitudinal eigenmode confinement and leakage;
- DC forced-helper localization;
- harmonic and sideband radiation into every receiving branch;
- threshold behavior at `KW_ZERO_LOCUS`;
- full-slab corrections to the simple homogeneous phase-matching formula;
- whether a finite-energy three-dimensional branch is allowed;
- whether helper-field elimination produces focusing of the required sign.

This stage should reject unsuitable constitutive branches early.

### 13.8 Stage L7 — solve reduced solitary-wave families

Proceed in increasing difficulty.

1. **One-dimensional traveling pulse.** Determine whether nonlinearity and dispersion can balance at all.
2. **Transverse stability test.** Perturb the one-dimensional pulse in the two lateral brane directions.
3. **Two-dimensional localized envelope.** Test one lateral direction plus propagation localization.
4. **Full three-dimensional brane envelope.** Seek localization in \(x,y,z\) with the normal profile already guided in \(w\).
5. **Vector polarization families.** Continue linear, circular, and general elliptical polarization states.

A one-dimensional soliton is a useful precursor but is not yet a photon.

### 13.9 Stage L8 — continue to the full parent fields

A reduced envelope solution is only a seed.

Reconstruct or directly solve the parent fields in four spatial dimensions and verify:

- the predicted normal profile;
- finite total energy;
- absence of hidden secular tails;
- consistency of the helper dressing;
- actual leakage;
- agreement with reduced coefficients;
- nonlinear stability.

### 13.10 Stage L9 — long-time stability and Floquet/response analysis

If the photon candidate is internally oscillatory, analyze perturbations around the traveling periodic solution.

Determine:

- translational zero modes;
- phase mode;
- polarization modes;
- width and shape modes;
- longitudinal/helper modes;
- growing instabilities;
- radiation continua;
- lifetime for a metastable resonance.

### 13.11 Stage L10 — collision tests

Evolve two photon candidates over:

- relative direction;
- impact parameter;
- polarization;
- phase;
- carrier frequency;
- amplitude.

Measure:

- shape recovery;
- phase and position shifts;
- frequency conversion;
- helper-field radiation;
- bulk leakage;
- fusion or breakup;
- energy and momentum conservation.

Weak interaction or soliton-like passage is required for a vacuum-like light sector.

### 13.12 Stage L11 — throat coupling, later

Only after a free photon branch is understood should the throat calculation ask:

- whether the same transverse mode can form the trapped support wave;
- how a throat emits and absorbs the packet;
- how acceleration sources both polarizations;
- whether the same coupling normalizes electricity, magnetism, and radiation;
- whether the packet can enter, remain in, or escape a throat cavity.

The free-light program therefore reduces uncertainty in the later supported-throat problem.

## 14. Acceptance criteria for a classical photon analog

A successful candidate should satisfy all of the following within a declared regime.

### 14.1 Background and spectral criteria

1. The stable slab exists and has bounded energy.
2. The lowest transverse guided branch is gapless.
3. The appropriate homogeneous limit reproduces the closed two-transverse/one-longitudinal finite-mode census with zero selected quadratic cross-block, while the complete slab retains exactly two healthy freely selectable transverse photon polarizations.
4. The normal profile is spectrally isolated and localized, or is an acceptably long-lived calculated resonance, as judged by the derived asymptotics of \(\Delta_{\rm spec}\), \(\Gamma_0\), \(\ell_w\), \(\mathcal P_{\rm br}\), \(\mathcal I_{\rm 3D}\), and the brane-observable pole residue where defined.
5. The relevant frequency family lies on one connected physical branch.

Criterion 4 is not satisfied merely by rescaling the eigenfunction to unit norm. Nor does it require a \(k\)-independent positive spectral gap or a \(k\)-independent localization length. The asymptotic scaling must preserve the declared effective three-dimensional regime and an appreciable coupling to brane observables.

The complete slab calculation must also classify the free longitudinal branch relative to the bulk essential spectrum and determine whether `KW_ZERO_LOCUS` survives, shifts, splits, or is replaced. Kinematic opening of a bulk channel is not itself evidence of leakage without nonzero interface overlap and outgoing flux.

### 14.2 Localization criteria

6. Total energy and momentum are finite.
7. Energy and stress are localized across \(w\) and within all three brane directions.
8. The packet does not spread secularly over the tested lifetime.
9. The packet does not collapse to the cutoff scale.
10. The helper response remains localized and attached.

### 14.3 Polarization and symmetry criteria

11. There are only two freely selectable photon polarizations.
12. The helper field is slaved rather than independently selectable.
13. The packet has no leading electric-odd monopole on a reflection-symmetric background.
14. Linear and circular polarization branches do not acquire unacceptable speed or energy splitting.

### 14.4 Propagation criteria

15. The packet propagates at or extremely near \(c_\gamma\).
16. Speed dependence on amplitude, width, or frequency is acceptably small in the intended regime.
17. The energy-momentum relation is approximately light-like.
18. There is no persistent longitudinal, bulk, electric, thickness, or conversion wake from the DC helper, carrier harmonics, sidebands, or higher nonlinear source components.
19. Any metastable leakage lifetime is sufficiently long.

### 14.5 Interaction criteria

20. Two packets pass with weak or soliton-like interaction in the ordinary regime.
21. Collisions do not generically fuse, collapse, or radiate a large fraction of energy.
22. Polarization-dependent interactions remain acceptable.

### 14.6 Cross-sector criteria

23. The same coefficients leave the far-field electric, magnetic, gravity, and unwanted-mode audits viable.
24. The helper mode does not create an unacceptable new long-range force.
25. The photon branch does not force a zero critical velocity or vacuum drag for ordinary throats.
26. No coefficient is retuned only for the photon branch.

## 15. Failure conditions

The selected photon branch fails if any required result proves impossible, including:

1. no stable finite-thickness shear-supporting slab exists;
2. no gapless guided transverse branch exists;
3. the lowest guided branch is unacceptably leaky into the bulk, merges into the relevant essential spectrum without an acceptable resonance, delocalizes so that \(\mathcal P_{\rm br}\) or the brane-observable pole residue vanishes unacceptably, or has divergent or vanishing fixed-on-brane-normalization effective inertia incompatible with the intended three-dimensional mode;
4. the background supports the wrong number of transverse polarizations;
5. polarization degeneracy is unacceptably broken;
6. the finite slab supplies no usable dispersion or nonlocality and every packet spreads;
7. every allowed nonlinear response is defocusing when focusing is required;
8. every focusing branch collapses or has energy unbounded below;
9. no reflection-even helper mode can bind to the packet;
10. the DC helper or any unavoidable carrier harmonic or sideband necessarily radiates an unacceptable persistent longitudinal, bulk, thickness, density, order, electric, or conversion wake;
11. the helper becomes a third freely selectable photon polarization;
12. the packet necessarily carries a leading electric-odd monopole;
13. localization exists only at a single isolated frequency with no acceptable light spectrum;
14. the packet speed varies unacceptably with intensity or frequency;
15. photon-photon interactions are generically strong, inelastic, or destructive;
16. a viable packet requires constitutive coefficients incompatible with the stable slab or force sectors;
17. a viable result exists only after introducing a photon-specific fitted term not propagated through the common medium;
18. the parent-field solution contradicts the reduced envelope solution;
19. no sufficiently long-lived three-dimensional packet exists despite exhausting the allowed shared constitutive branches.

Failure of the self-bound branch would not erase the conditional linear two-transverse-mode result. It would mean that this medium supplies classical light waves but does not close the intended localized free-photon ontology.

## 16. Interaction with the far-field force program

The photon and force programs should share a **medium capability matrix**.

For every background eigenmode, record:

- field composition;
- reflection parity;
- transverse or longitudinal character;
- normal profile across \(w\);
- static small-\(k\) stiffness;
- dynamic dispersion;
- gap or screening length;
- essential-spectrum placement and any grazing threshold;
- group velocity;
- energy sign;
- damping or leakage;
- allowed source couplings;
- interface and mode-projector overlaps;
- DC, harmonic, and sideband source support;
- role in light, electricity, magnetism, gravity, inertia, or unwanted radiation.

The light branch needs the finite-\(k\), nonlinear part of this matrix. The force branch needs the low-\(k\), static and retarded part.

A common parameter search should determine whether one region simultaneously provides:

- stable slab thickness;
- healthy transverse guided light;
- acceptable longitudinal/helper behavior;
- a healthy odd \(k^2\) Coulomb carrier;
- screened unwanted even modes;
- acceptable gravity response;
- no ghosts;
- no low-threshold drag;
- photon localization.

A failure to find a common region is more informative than separate success in individually tuned sectors.

## 17. Frozen light ledger

The following entries should be added to the broader provenance ledger as this program advances.

| Claim or quantity | Current status | Required evidence | Downstream dependencies | Failure condition |
|---|---|---|---|---|
| Homogeneous \(D=3\) mode census and quadratic block separation | Derived under stated assumptions | [S11](../research/pde_ledger_v3/steps/S11_stray_longitudinal.md): \(\omega_T^2=(\mu_R/\rho_{\rm br})k^2\) with multiplicity 2, \(\omega_L^2=(B_{\rm comp}/\rho_{\rm br})k^2\) with multiplicity 1, cubic characteristic polynomial with nonzero leading coefficient for \(\rho_{\rm br}\ne0\), and zero selected \(T\)-\(L\) cross-block | Free light, helper ordering, support mode, magnetism, radiation | Wrong mode count, instability, or use outside the homogeneous isotropic parity-even quadratic scope without reopening the result |
| Stable finite slab | Open | Solved \(n_0(w),\chi_0(w)\) with selected thickness and spectrum | Every light result | No stable or metastable slab |
| Lowest guided transverse profile \(f_0(w;k)\) | Target | Variable-coefficient normal eigenproblem with \(\Delta_{\rm spec}\), \(\Gamma_0\), \(\ell_w\), \(\mathcal P_{\rm br}\), fixed-on-brane \(\mathcal I_{\rm 3D}\), regulator tests for resonances, and a pole residue once an observable is specified | Confinement, effective inertia, dispersion, photon packet | No bound or long-lived resonance branch, unacceptable delocalization or leakage, vanishing brane participation/residue, or singular effective inertia |
| Gapless guided photon branch | Target | \(\lambda_0=0\), \(\omega_0(k)\to0\) | Light-like propagation | Lowest branch remains gapped |
| Two guided polarizations | Target reduction | Complete constrained guided spectrum | Photon identity | Extra or missing polarization |
| Finite-thickness dispersion | Open | Coefficients \(a_4,a_6,\ldots\) from projection | Soliton balance | No compatible dispersion |
| KW grazing threshold | Derived kinematic threshold under the simple homogeneous phase-matching model | \(k_w^2=k^2(c_L^2/c_{s0}^2-1)\), with **KW_ZERO_LOCUS** \(B_{\rm comp}=\rho_{\rm br}c_{s0}^2\) or \(c_L=c_{s0}\) | Free-longitudinal continuum placement and full-slab threshold audit | Threshold is treated as proof of coupling, leakage, boundness, or absence of a bound state in the continuum |
| Full interface leakage closure | Open / S11b | Complete brane–bulk operator, kinetic norm, essential spectrum, interface matching, bound/threshold/BIC/resonance classification, source overlap, outgoing flux, and leakage width | Light confinement, longitudinal observability, drag, packet lifetime | Unacceptable leakage or no mathematically closed spectral classification |
| Longitudinal/even helper mode | Open hypothesis | Derived mixed eigenmode, coupling to \(I\), and placement relative to all receiving continua | Self-guiding | No localized nonradiative helper |
| Nonlinear helper radiation closure | Open | Complete DC, harmonic, and sideband source spectrum; mode projectors; interface overlaps; retarded response; emitted powers; and threshold behavior | Photon localization and lifetime | Any unavoidable source component radiates unacceptably |
| Nonlinear envelope coefficients | Target | Projection of shared constitutive law | Existence and stability | Wrong sign, collapse, or independent tuning |
| Three-dimensional localized packet | Target relation | Reduced and parent-field traveling solution | Classical photon ontology | No finite-energy stable branch |
| Packet stability | Target | Linear, Floquet, or retarded spectrum | Lifetime and interactions | Growing mode or unacceptable leakage |
| Photon collision behavior | Target | Time-dependent two-packet simulations | Vacuum linearity analog | Strong generic inelasticity |
| Energy-frequency and energy-momentum laws | Open | Conserved charges of packet family | Photon phenomenology | Unacceptable non-light-like behavior |
| Quantized action or \(E=\hbar\omega\) analog | Outside present classical closure / later target | Additional mechanism not yet specified | Quantum photon behavior | Not used to judge the classical packet stage |

## 18. Recommended computational artifacts

Each stage should produce reproducible symbolic and numerical outputs.

### Symbolic outputs

- parent quadratic operator;
- normal-mode reduction;
- polarization projectors;
- homogeneous characteristic polynomial and finite-root census;
- simple and full-slab \(k_w=0\) threshold conditions;
- small- and finite-\(k\) dispersion expansions;
- nonlinear invariant inventory;
- DC, harmonic, and sideband source spectra;
- receiving-mode projector and interface-overlap formulas;
- projected coupling coefficients;
- envelope/helper equations;
- energy, momentum, and virial identities;
- resonance and emission conditions.

### Numerical outputs

- slab profiles;
- guided mode profiles across \(w\);
- essential-spectrum edges and spectral-separation curves;
- free-longitudinal and forced-source continuum-placement maps;
- parameter sweeps through and around **KW_ZERO_LOCUS**;
- resonance-width and lifetime curves;
- localization-length, brane-participation, fixed-on-brane-inertia, and applicable pole-residue curves;
- regulator-convergence tests for outgoing resonances;
- dispersion surfaces;
- helper Green functions;
- one-dimensional solitary branches;
- transverse-instability spectra;
- full three-dimensional packets;
- long-time evolution;
- collision scans;
- parameter compatibility maps shared with the far-field audit.

### Documentation outputs

For each derived claim, record:

- exact equations and branch assumptions;
- frozen parameters;
- boundary conditions;
- solver and convergence tests;
- conserved or dissipated quantities;
- status: precursor, conditional derivation, target, or failure;
- downstream sectors reopened by any change.

## 19. Immediate next questions

The first practical questions are deliberately narrower than “solve the photon.”

1. What transversely isotropic shear action is being postulated for the ordered state?
2. How do its coefficients depend on \(n\) and \(\chi_B\)?
3. Does a symmetric finite slab possess a gapless bound transverse branch or an acceptably long-lived outgoing resonance?
4. What are the branch's exact normal profile \(f_0(w;k)\), spectral separation, width, localization length, brane participation, and fixed-on-brane effective inertia?
5. What higher-order dispersion follows from finite thickness?
6. How do \(v_{\rm packet}\), \(c_\gamma\), \(c_L\), and \(c_{s0}\) compare, and how does the complete slab behave at and near the simple **KW_ZERO_LOCUS**?
7. Which even eigenmode is most strongly sourced by transverse intensity?
8. Does eliminating that mode produce focusing, defocusing, or nonlocal response, and do any of its DC, harmonic, or sideband source components radiate?
9. Is the nonlinear energy bounded below?
10. Can a one-dimensional traveling envelope exist?
11. Is that envelope stable to the two lateral brane directions?
12. Can one common coefficient set satisfy both this photon audit and the low-\(k\) force audit?

Answering these questions would determine whether a full three-dimensional photon solve is promising before committing to the hardest numerical continuation.

## 20. Compact interpretation

One conserved four-dimensional medium forms an ordered finite-thickness brane. The brane has a globally defined normal axis but remains isotropic within its three observed directions. The normal **sign** is available to oriented throats and electric charge, while the unsigned normal axis organizes light guidance.

Light is postulated to be the lowest gapless transverse shear branch of that slab. Its displacement is tangent to the brane and transverse to its direction of propagation, leaving two independent polarizations. Its field has a bound profile across \(w\), while its energy travels along the brane.

The finite slab can prevent normal leakage and supply a material length scale and dispersion, but normal guidance alone does not localize a packet within all three brane directions. The photon candidate therefore carries an oscillatory transverse carrier inside a finite envelope. Its intensity excites a reflection-even material response—most likely a mixture of longitudinal compression, thickness, density, and order change. That response is slaved to the two transverse polarization amplitudes and creates a moving guide that opposes in-brane diffraction.

The central candidate structure is

\[
\boxed{
\text{gapless guided transverse carrier}
+
\text{finite-thickness dispersion}
+
\text{slaved even material dressing}
\longrightarrow
\text{self-guided classical photon packet}.
}
\]

The longitudinal mode is therefore not redefined as charge or promoted to a third photon polarization. It is investigated as a possible bound support field for the packet. If its dispersion permits a comoving nonradiative response, it may be essential. If it necessarily forms a wake, causes drag, or destabilizes the brane, it becomes a failure channel.

The homogeneous \(D=3\) quadratic theory now supplies a closed three-mode census and exact linear transverse–longitudinal separation under the selected action. It also supplies a genuine free-longitudinal/bulk-sound grazing threshold at \(c_L=c_{s0}\). These results strengthen the weak-field light sector but do not close photon confinement or longitudinal leakage. The photon-soliton program therefore proceeds with a purely transverse linear carrier and a nonlinearly induced even helper, while the complete slab and retarded-response calculations determine whether the helper and all of its harmonics remain bound or radiate.

This program is largely independent of the full throat solve and should run in parallel with the far-field force audit. Both depend on the same stable slab and complete background spectrum, and both must use the same frozen constitutive coefficients. The photon track tests finite-wave-number nonlinear localization; the force track tests long-wavelength carriers, ranges, and couplings.

The decisive question is not whether the current equations accidentally contain a photon soliton. The medium may be designed to contain the required guided and nonlinear structures. The decisive question is whether one such frozen medium remains stable and supports the photon, force, throat, radiation, and reservoir sectors simultaneously without incompatible requirements or independent retuning.
