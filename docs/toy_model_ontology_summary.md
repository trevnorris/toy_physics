# Ontology and Closure Ledger of the Toy Analog Physics Model

## Purpose and scope

This document is the canonical statement of what the model postulates, defines operationally, derives under stated assumptions, leaves open, and treats as a failure condition. It also serves as the frozen derivation-provenance and closure ledger. It describes a **classical toy analog**, not a claim that the real universe is literally made this way.

The focused [opposite-orientation derivation plan](opposite_orientation_throat_coincidence.md) imports definitions and epistemic status from this ontology and develops the detailed two-interface contact calculation; it is not a second canonical source of definitions.

**Conceptual freeze rule.** After the present technical cleanup, new mechanisms, coefficients, or ontology branches should be added only when an explicit derivation exposes a missing assumption or a failed closure test. The next default step is calculation within this ledger, not further speculative expansion.

The project asks whether one internally consistent medium can support gravity-like, light-like, electric, magnetic, inertial, and cosmological behavior at the same time. The medium and its constitutive laws may be postulated. The test is whether the resulting system is mathematically coherent, whether its sectors can coexist without independent retuning, and which departures from GR or Maxwell it predicts.

The following epistemic labels apply throughout:

- **Postulated:** part of the candidate model's definition.
- **Operational definition:** a calibrated source or response quantity used by a brane observer; it is not by itself a microscopic derivation.
- **Target relation:** a result the completed model must produce, not an input that may be imposed to obtain the desired phenomenology.
- **Derived under stated assumptions:** follows from postulated equations only within the named ansatz, branch, boundary ensemble, or idealized limit.
- **Open hypothesis:** physically motivated, but not yet established by the model's derivations.
- **Failure criterion:** a result whose impossibility or unacceptable consequence rejects the candidate architecture or a specified branch.

No mathematical variable should be mistaken for a separate substance merely because it is convenient to describe one aspect of the medium.

## One-sentence description

The model is a classical \(4+1\)-dimensional, one-medium analog in which our three-dimensional space is an ordered, shear-supporting brane; particles are finite, oriented throats through it; a trapped transverse brane-shear standing mode helps hold each throat open; throat transport, order conversion, stress, return, and reservoir response are the proposed material mechanism underlying a calibrated gravity observable; and freely propagating light, electric coupling, magnetism, radiation, passive response, inertia, and possible cosmic expansion are target behaviors of the same complete material configuration.

## 1. The arena

The model assumes four spatial coordinates. Write a projected brane point as \(\mathbf x=(x,y,z)\) and a full spatial point as

\[
\mathbf X=(\mathbf x,w),
\]

and one time coordinate \(t\). The familiar directions \(x,y,z\) lie along the brane. The coordinate \(w\) is normal to it and points into the bulk.

The postulated parent arena is

\[
\boxed{
\mathbb R_t\times(\mathbb R^4,\delta_{AB})
},
\]

with preferred Newtonian time \(t\) and flat Euclidean spatial metric \(\delta_{AB}\). The operators \(\nabla_4\), spatial dot products, reflection \(\mathcal R_w\), and norms such as \(\|\mathbf Y_+-\mathbf Y_-\|_4\) use this ambient metric. Constitutive coefficients or a solved background may be anisotropic without changing this metric postulate. If a future branch replaces the flat ambient geometry, Euclidean set distance must be replaced by the corresponding geodesic distance and every dependent Green function and wave operator reopened.

This background arena is part of the model's setup. The model does not currently attempt to derive space or time themselves. Instead, it asks whether effective \(3+1\)-dimensional physics can emerge for observers confined to a brane inside this \(4+1\)-dimensional arena.

## 2. The one fundamental substance

Everything in the model is made from one compressible condensate. At mean-field level it is represented by a complex field

\[
\psi=\sqrt n\,e^{i\theta},
\]

where:

- \(n=|\psi|^2\) is the constituent number density;
- \(\theta\) is the condensate phase;
- gradients of \(\theta\) generate the four-dimensional material velocity \(\mathbf v_{\rm med}\) in the conservative Madelung branch;
- the density and velocity live throughout all four spatial dimensions.

The total constituent current \(\mathbf J_n\) is defined by exact parent-medium conservation,

\[
\boxed{
\partial_t n+\nabla_4\cdot\mathbf J_n=0
}.
\]

In the conservative Madelung branch,

\[
\mathbf J_n=n\mathbf v_{\rm med}.
\]

A coarse-grained relaxational or mixed branch may instead have

\[
\mathbf J_n
=n\mathbf v_{\rm med}+\mathbf J_n^{\rm rel},
\]

but any relative or diffusive current \(\mathbf J_n^{\rm rel}\) must be derived with its momentum, energy, and entropy partners. It is not an additional material source. Thus \(n\) has constituent-number-per-four-volume units, while an aperture integral of \(\mathbf J_n\) has constituent-number-per-time units.

The working substrate is a generalized nonlinear-Schrödinger or Gross–Pitaevskii medium with a stiff polytropic equation of state, presently written

\[
P=K n^5.
\]

The present conservative substrate action, equation of state, and microscopic constants are postulates of the analog. The full parent dynamical framework is not closed and may require dissipative or reservoir variables rather than a conventional action. Any such reservoir variables must represent unresolved or coarse-grained degrees of freedom of the same globally conserved substrate; they are not a second substance or an external source of material or energy. The model does not yet derive either branch from a deeper constituent theory.

The medium may possess density, momentum, pressure, effective inertia densities, and kinetic coefficients. These properties must not be interpreted as an inventory of little conventional massive particles whose rest masses add up to an observed particle's mass. A substrate coefficient such as a condensate mass parameter is not any of the throat's active, passive, or inertial mass coefficients. Those are distinct emergent source and response properties of the complete dressed throat.

The medium is intended to be globally conserved. Apparent sources and sinks seen by a brane observer must ultimately be projections of four-dimensional transport or changes in the material's order, not creation or destruction of a second substance.

The bulk is therefore an internal reservoir, not an external supply. In a closed-loop realization, localized throat drainage transfers material, momentum, and energy from ordered brane degrees of freedom into de-structured bulk degrees of freedom, while distributed return transfers material back into the ordered state. The bulk may store the corresponding response as compression, flow, internal energy, entropy, or other unresolved excitations of the same medium. Any sustained cycle must include the evolution of that reservoir rather than treating it as an inexhaustible battery.

## 3. Two material states of the same medium

An order variable \(\chi_B\) distinguishes two states of the medium:

\[
\chi_B\simeq1
\quad\Longrightarrow\quad
\text{ordered, brane-like, shear-supporting material},
\]

\[
\chi_B\simeq0
\quad\Longrightarrow\quad
\text{de-structured, bulk-like material}.
\]

The brane and bulk are therefore not two substances. They are two states of one substance. Material can de-structure from the brane state into the bulk state and can reorder from the bulk state into the brane state.

The symbol \(\chi_B\) always denotes the material order variable. It is interpreted literally as an ordered fraction only on the fraction/relaxational branch. If an inertial or mixed branch requires both an order field and a separate coarse-grained ordered material fraction, reserve \(f_B\in[0,1]\) for that fraction and derive its coupling to \(\chi_B\).

The existence of this two-state structure and the potential that supports it are major postulates. The final theory must also make a major branch choice about the dynamics of \(\chi_B\).

### 3.1 Dissipative branch

If \(\chi_B\) is a dissipative material fraction:

- drainage and return produce entropy;
- a permanent throat may require a maintained nonequilibrium state;
- the energy source and return pathway must be identified;
- particle longevity must be reconciled with continual conversion.

Any heat, bath, or reservoir variables used to close this branch must be internal degrees of freedom of the one medium. Its complete ledger must still close material, momentum, energy, and entropy.

In this branch, de-structuring may be understood as loss of **coarse-grained structural memory**: material that ceases to support shear can no longer retain the resolved elastic stress or geometric imprint carried by the ordered brane. That resolved information is not fundamentally destroyed; it is transferred into heat, disorder, and other unresolved internal degrees of freedom and appears as entropy production. This is a candidate physical interpretation of the dissipative branch, not a derivation of the branch choice or a claim that reversible order conversion is impossible.

### 3.2 Inertial or reversible branch

If \(\chi_B\) is an inertial field:

- order conversion carries canonical energy and momentum;
- additional conversion waves or oscillations should exist;
- an apparent drain may be part of coherent transport rather than irreversible consumption.

### 3.3 Explicitly mixed inertial-dissipative branch

A physical effective theory may contain both reversible order dynamics and damping into identified internal degrees of freedom. Schematically, but not as a selected final equation,

\[
\ddot\chi_B
+\Gamma_\chi\dot\chi_B
+\frac{\delta E}{\delta\chi_B}
=\text{explicit internal-reservoir coupling}.
\]

The final parent theory must explicitly select and close a conservative, relaxational, or mixed inertial-dissipative order branch. What is forbidden is leaving the branch unspecified or using damping without its internal reservoir, energy, momentum, and entropy partners. Every selected branch must close constituent conservation, momentum and energy accounting, a consistent inertia/drag formalism, and entropy production where damping is present. If the effective description requires fluctuation/noise partners or records their omission, those consequences must also be stated. The branch choice affects gravity, throat stability, leakage, cosmology, and the global reservoir loop.

Branch choice also controls which mathematical objects are legitimate. A conservative or genuine equilibrium static problem may use a free-energy Hessian \(\mathcal H\), its inverse \(\mathcal G_{\rm stat}\), and an ensemble-appropriate interaction potential. Dynamic response uses a retarded operator \(\mathcal O_R(\omega,\mathbf k)\) and \(\mathcal G_R=\mathcal O_R^{-1}\) to determine poles, damping, phase lag, and radiation. A conservative but intrinsically periodic supported throat or pair is not automatically an equilibrium static problem: its interaction may require a cycle-averaged Hamiltonian, Routhian, quasienergy, or averaged action in an explicitly derived ensemble with fixed support-mode action, frequency, phase relation, or other conjugate data. If no such scalar functional is established, its mean force must instead be derived from the cycle-averaged Noether stress and momentum flux. For a maintained dissipative stationary pair, force must be obtained from stress, momentum flux, conversion, and reservoir response; a scalar pair potential may be introduced only after that force is shown to be conservative and path independent.

## 4. Field and variable dictionary

The model uses several descriptions of the same medium at different levels:

| Level | Variables | Meaning |
|---|---|---|
| Parent arena and medium | \(t,\delta_{AB},\mathbf x,\mathbf X=(\mathbf x,w),\psi,n,\theta,\chi_B,f_B,\mathbf v_{\rm med},\mathbf J_n\) | Preferred time, ambient Euclidean spatial metric, projected and full spatial points, condensate, density, phase, order field, optional distinct coarse-grained ordered fraction, material velocity, and total constituent current |
| Background brane | \(n_0(w),\chi_0(w),H_{\rm br}\) | Brane profile and local selected slab thickness |
| Brane collective modes | \(\mathbf u_T,u_L,h_+,h_-,h_m,h_t,\phi_E\) | Transverse and longitudinal motion; upper- and lower-interface displacement; mid-surface and thickness combinations; and schematic notation for the Coulomb-carrying eigencombination if one is derived |
| Linear operators | \(\mathcal H,\mathcal G_{\rm stat},\mathcal O_R,\mathcal G_R\) | Equilibrium free-energy Hessian and static Green matrix versus dynamic retarded operator and retarded response matrix |
| Throat coordinates | \(s_a,\mathbf x_a,R_a\) | Orientation, projected three-dimensional position, and characteristic size of throat \(a\) |
| Reflection and antiparticle candidate | \(\mathcal R_w,\mathcal I_{\rm internal},\mathcal C_{\rm th}\) | Geometric \(w\)-reflection, any additional transformation of odd internal data, and the candidate complete throat map \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\) |
| Interface versus core geometry | \(h_\pm,h_m,h_t\) versus \(\chi_B,n,\theta\), support fields, and core level sets | Graph-regime interface coordinates versus the full nonlinear four-dimensional throat geometry |
| Transport and conversion | \(\widetilde Q_n^{(w)},Q_n^{\rm net},Q_n,Q_\chi,Q_{\rm ap}\) | Signed global-\(w\) coordinate flux, signed local net outward flux, positive net outward throughput on a selected drain branch, branch-specific nonnegative gross local order-conversion drain, and optional gross-aperture shorthand |
| Gravity source and calibration | \(\widetilde Q_g^{(s)}(r),\overline Q_g(r),\boldsymbol{\mathfrak g}_{\rm src},\mathcal C_{\rm ref},\mathbf g_{\rm eff},C_g,G_{\rm eff}\) | Raw orientation-labeled source coefficient; reflection-even but radially signed enclosed coefficient; source-side environmental field; fixed operational calibration; calibrated observable; inverse-square source normalization; and effective gravitational constant in the observer's chosen units |
| Charge quantities | \(s_a,N_s,g_a,C_{E,ab}\) | Throat orientation, net oriented count, possible derived nonnegative leading one-throat coupling magnitude, and coefficient of the reflection-odd Coulomb monopole channel only |
| Other observables | \(\mathbf B_{\rm eff},M_{\rm active}^{\rm local},M_{\rm active}^{\rm eff}(r),M_{\rm passive},M_{\rm inertial}\) | Magnetic response, local and scale-dependent active mass, and the passive and inertial mass roles |

This hierarchy preserves the one-medium ontology. The collective coordinates are effective descriptions of material motion or order; they are not additional fundamental substances. The symbol \(\phi_E\) is schematic notation for whichever eigencombination of the coupled \(h_+\)-\(h_-\)-\(\chi_B\) and other retained fields has the required Coulomb-carrying \(k^2\) static stiffness, if such a branch exists. It is not an additional fundamental field. Likewise,

\[
C_{E,ab}(R)=g_ag_b+\delta C_{ab}^{\rm odd}(R)
\]

is a target decomposition of the reflection-odd Coulomb monopole coefficient, not a derived identity or a container for even thickness and nonlinear core forces. In any such factorized convention, the \(g_a\) are nonnegative magnitudes; electric sign remains in \(s_a\).

Let \(\mathcal A_a\) be the selected \(w\)-normal three-dimensional reference cut through an axis-aligned throat, with one global orientation \(\hat{\mathbf w}\). The constituent-flux convention on this cut is

\[
\widetilde Q_{n,a}^{(w)}
\equiv
\int_{\mathcal A_a}
\mathbf J_n\cdot\hat{\mathbf w}\,d^3\Sigma,
\]

where \(\widetilde Q_{n,a}^{(w)}\) is signed in the global \(w\) coordinate. With local outward throat normal

\[
\hat{\mathbf n}_a=s_a\hat{\mathbf w},
\]

the signed local net outward flux is

\[
Q_{n,a}^{\rm net}
\equiv
\int_{\mathcal A_a}
\mathbf J_n\cdot\hat{\mathbf n}_a\,d^3\Sigma.
\]

The coordinate relation

\[
\boxed{
\widetilde Q_{n,a}^{(w)}=s_aQ_{n,a}^{\rm net}
}
\]

is an **identity** following from the two normal conventions on the selected reference cut. For a general curved cut with normal field \(\hat{\mathbf n}(\mathbf X)\), the invariant flux is \(\int\mathbf J_n\cdot\hat{\mathbf n}(\mathbf X)\,d^3\Sigma\); the simple global-\(w\) identity is then conditional on the geometry of that cut and may require a projected-area relation. For a solved stationary outward-drain branch, define the selected positive throughput only after verifying the branch sign:

\[
Q_{n,a}\equiv Q_{n,a}^{\rm net}>0.
\]

Positivity is a property of the selected solution, not a consequence of the surface-integral notation. Under reflection-symmetric boundary and reservoir data, the actual physical targets for independently solved \(\mathcal C_{\rm th}\)-related partners are

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

In the fraction/relaxational branch, when \(\Gamma_{\rm drain}\) is a local nonnegative conversion rate per unit time, define

\[
\boxed{
Q_{\chi,a}
=\int_{\Omega_{{\rm throat},a}}
n\Gamma_{\rm drain}\,d^4X
\ge0
}.
\]

It has the same constituent-number-per-time dimensions as \(Q_{n,a}\) but need not equal it. It excludes distributed return and is not a signed net rate. In an inertial order-field branch, an equivalent gross conversion diagnostic and its units must be derived from that dynamics rather than importing this constitutive formula. Signed net conversion remains in \(\Gamma_{\rm drain}\), \(\Gamma_{\rm return}\), and \(\Gamma_B\) only on branches where those rates are defined.

## 5. The brane: our effective three-dimensional space

Our observed space is identified with a finite ordered slab or wall region near \(w=0\), bounded by bulk material on both sides. Having bulk on both sides is important because a throat may puncture toward either \(+w\) or \(-w\).

**Postulated global orientation.** The present candidate assumes that the brane is globally two-sided and co-orientable, with a continuous normal convention inherited from the global \(w\) coordinate. This allows \(+w/-w\) throat orientations to be compared consistently at separated points. If later solutions admit non-co-orientable topology or no global normal bundle, the charge definition must be reformulated rather than used unchanged.

The symbol \(\mathcal R_w\) denotes geometric reflection of the parent fields and environment. The complete antiparticle candidate may additionally require \(\mathcal I_{\rm internal}\), so the candidate throat operation is

\[
\boxed{
\mathcal C_{\rm th}
=\mathcal I_{\rm internal}\circ\mathcal R_w
}.
\]

The solved species branch must determine whether \(\mathcal I_{\rm internal}\) is trivial. Reflection covariance compares \(\mathcal C_{\rm th}\mathcal T\) with the geometrically reflected environment \(\mathcal R_w\mathcal B_\infty\). Partners are degenerate in one unchanged vacuum only when that environment is itself symmetric, \(\mathcal R_w\mathcal B_\infty=\mathcal B_\infty\). Unequal bulk-facing reservoirs may instead produce a calculable orientation-dependent environmental splitting; that is distinct from failure of the intrinsic candidate throat map.

The coordinate \(w=0\) is not an externally preferred absolute plane. It is defined relationally as the local reflection-symmetric mid-surface of the ordered slab between its two bulk-facing interfaces. In the working formation picture, this is the meeting and symmetry surface of the counterposed bulk regions at which the ordered slab forms. The \(+w\) and \(-w\) bulk regions are two sides of the same globally conserved medium, not separate substances. For a symmetric stable slab solution,

\[
\chi_0(w)=\chi_0(-w),
\qquad
\chi_0(0)=\max_w\chi_0(w),
\]

so the ordered state is strongest at the mid-surface and weakens toward both interfaces. Translating the whole brane merely redefines the coordinate origin; only displacement relative to its local mid-surface is physical. This symmetry identifies the center once a slab exists, but does not by itself derive brane formation, finite thickness, or stability.

A finite slab generally requires two nonlinear interface coordinates rather than one single-valued height. Write

\[
w_+(\mathbf x,t)=+\frac{H_{\rm br}}2+h_+(\mathbf x,t),
\qquad
w_-(\mathbf x,t)=-\frac{H_{\rm br}}2+h_-(\mathbf x,t).
\]

Their mid-surface and half-thickness combinations are

\[
h_m=\frac{h_++h_-}{2},
\qquad
h_t=\frac{h_+-h_-}{2},
\qquad
\delta H=2h_t.
\]

Under reflection through the slab,

\[
\boxed{
h_+\to-h_-,
\qquad
h_-\to-h_+
},
\]

and therefore

\[
h_m\to-h_m,
\qquad
h_t\to+h_t.
\]

Both interface displacements acquire the required minus sign while their upper/lower labels are exchanged. Uniform \(h_m\) translation must remain unpinned when the parent medium is translation invariant, whereas stable thickness selection may give a primarily thickness-like eigenmode a restoring scale. Parity alone does not determine range because these coordinates may mix with order, density, phase, longitudinal, bulk, or reservoir variables.

The fields \(h_+\) and \(h_-\) describe the interfaces only while each remains a single-valued graph over the brane coordinates. They are useful background, far-field, and separated-interface collective coordinates, not a complete nonlinear throat topology. Overhangs, necks, tubes, folded sheets, reconnection, pinch-off, and topology change must be represented by the parent fields \(\chi_B(\mathbf x,w),n(\mathbf x,w),\theta(\mathbf x,w)\), the support fields, and core level sets. The earlier single coordinate \(h\) should be read only as schematic notation for a derived normal eigencombination.

A viable background brane must possess:

- finite, dynamically selected thickness;
- an energy functional bounded below, with interfaces that are dynamically stable or sufficiently metastable;
- in-plane shear support;
- a well-posed perturbation spectrum;
- the required transverse and longitudinal collective modes.

These are prerequisites for a stable throat, not secondary details. The current order potential admits interfaces but does not by itself select the width of a finite ordered slab. A complex condensate plus a scalar order fraction also does not automatically generate solid-like transverse shear. The shear-supporting constitutive action and its modulus must therefore either be explicitly postulated or derived from a richer order structure.

**Status:** the ordered two-state ontology and present wall/shear constitutive structure are postulated; interface profiles can be derived from those postulates; finite slab-thickness selection, microscopic shear provenance, and the complete stable background spectrum remain open.

Any future thickness-selection or shear-generation mechanism must be propagated through every dependent sector. It may change:

- the light speed \(c_\gamma\);
- the longitudinal spectrum;
- the electric-sector stiffness;
- the throat radius;
- leakage and conversion;
- electric and magnetic couplings.

Observers and excitations bound to the brane see effective three-dimensional fields obtained by projecting the full four-dimensional medium through a finite window in \(w\).

This projection has physical consequences. When four-dimensional continuity is integrated across a finite window, its edges leave a boundary-flux remainder. This projected term is called \(S_{\rm leak}\). Leakage is therefore a native possibility of the dimensional reduction rather than an arbitrary cosmological source added afterward.

## 6. Particles are finite supported throats

There are no point particles in the ontology. A particle is a finite, nonlinear **throat**, **puncture**, or **defect** in which the brane extends into the bulk. A throat is an oriented, open site of phase conversion or transport between ordered brane material and de-structured bulk material.

Write a solved throat branch as

\[
\mathcal T_{\alpha,s},
\qquad
s\in\{+1,-1\},
\]

where \(\alpha\) contains orientation-independent species data and \(s\) is the electric orientation. A change of \(s\) does not by itself change \(\alpha\) into another particle species.

The working decomposition of a particle is:

- **Particle identity:** the complete nonlinear throat, its internal wave, its near fields, and its surrounding flow.
- **Charge sign:** the direction in which it punctures, \(+w\) or \(-w\).
- **Antiparticle relation:** reversing \(s\) on the same solved internal throat branch gives the leading antiparticle candidate. Geometric reflection \(\mathcal R_w\) reverses the embedding orientation, while a possibly nontrivial \(\mathcal I_{\rm internal}\) must transform every additional odd phase, chirality, circulation, knot orientation, support-mode phase convention, topological index, or conversion-current datum and preserve the required even scalars. The complete candidate is \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\). An electron-like throat related by a complete consistent \(\mathcal C_{\rm th}\) is positron-like, not proton-like; neither \(\mathcal I_{\rm internal}\) nor the complete map is yet derived.
- **Species identity:** distinct particles may require different internal support modes, core topology, or composite throat structure. A proton-like positive object may not be represented merely by reversing the orientation of an electron-like elementary throat.
- **Structural support:** energy \(E_{\rm support}\) in a trapped transverse brane-shear standing mode—the model's light or photon mode—helps hold the aperture open.
- **Finite size:** a dynamical balance among standing-wave support, throat flow, brane tension, and bulk backpressure.
- **Constituent transport:** \(Q_n^{\rm net}\) is signed local net outward flux and obeys the coordinate identity \(\widetilde Q_n^{(w)}=sQ_n^{\rm net}\). On a solved stationary drain branch, \(Q_n\equiv Q_n^{\rm net}>0\) is the selected net outward throughput used by critical-flow estimates.
- **Order conversion:** \(Q_\chi\ge0\), the positive branch-specific gross ordered-to-bulk drain rate or derived diagnostic in the throat region, excluding distributed return.
- **Gross aperture process:** \(Q_{\rm ap}\), an optional umbrella shorthand for the complete aperture process when the distinction between \(Q_n\) and \(Q_\chi\) is not at issue. It does not assert \(Q_n=Q_\chi\).
- **Gravity-source strength:** \(\widetilde Q_g^{(s)}(r)\) is the provisional raw orientation-labeled projected coefficient. After reflection closure, \(\overline Q_g(r)\in\mathbb R\) is orientation-even but retains its radial sign after enclosed return. Attraction on the ordinary local branch is a separate sign test.
- **Active gravitational mass:** a local source property and a possibly scale-dependent effective source inferred from the gravity field.
- **Passive gravitational response:** the throat's response coefficient in an external gravity flow.
- **Inertial response:** the resistance encountered when the complete dressed throat configuration is accelerated.

For a reversible equilibrium branch, candidate antiparticle degeneracy includes equality of equilibrium energy. For a conservative periodic branch, the complete \(\mathcal C_{\rm th}\) test must also map support-mode action, frequency, and phase data and compare the appropriate cycle-averaged invariant or stress/momentum flux; an arbitrary relative phase is not a reflection partner. For a driven relaxational or mixed branch, the reflected stationary comparison must instead or additionally test equality of positive throughput, dissipation and entropy-production rates, lifetime, stationary stress and momentum fluxes, and retarded response kernels under reflected environmental data. Equality of a scalar energy is not assumed where no equilibrium or periodic variational functional exists.

### 6.1 Spectrally localized ordered support mode

The support mode is defined spectrally on the complete variable-coefficient parent geometry. It is not assumed to propagate as a brane-shear wave through fully de-structured material with \(\chi_B\simeq0\). The **target bound-state problem** is

\[
\mathcal L_T[\chi_B,n,\text{geometry}]\phi_T
=\omega^2
\mathcal W_T[\chi_B,n]\phi_T,
\]

For a true bound-state branch, require finite weighted norm, for example

\[
\boxed{
\int d^4X\,
\phi_T^*\mathcal W_T\phi_T<\infty
},
\]

acceptable asymptotic decay, and finite energy. A leaky resonance is a different spectral problem: impose outgoing conditions and derive a complex frequency or an equivalent decay rate and lifetime. It is not required to satisfy the bound state's normalizability condition. In an open, relaxational, or mixed realization, a resonance is acceptable only if its leakage and lifetime are derived and sufficiently small.

The ordered rim, collar, sheath, cavity wall, or other region in which the solved mode's energy and supporting stress are concentrated may be denoted diagnostically by \(\Omega_{\rm support}\). It is identified after solving the eigenproblem—for example as a region containing most of the modal energy—not imposed as a sharp PDE domain by \(\mu_R>0\). The mode's energy localization, supporting stress, coexistence with throughput, leakage into longitudinal, bulk, electric/order, or conversion modes, and behavior for opposite-interface cores at projected coincidence are **not yet derived**.

The linear eigenproblem is only the first stage of the supported-particle calculation:

1. A linear mode on a provisional throat geometry supplies a spectral seed and identifies candidate bound or resonant branches.
2. The actual particle must be a self-consistent nonlinear free-boundary or time-periodic solution in which the support field's time-averaged stress, throat geometry, flow, conversion, and reservoir loading agree. Harmonic balance, shooting, continuation, or direct time evolution are possible numerical methods, not additional physical postulates.
3. Stability is the Floquet or linear-response spectrum about that self-consistent periodic solution, not merely the seed eigenvalue problem.

The stationary solution is a coupled feedback system,

\[
\boxed{
E_{\rm support},\omega_{\rm standing}
\longleftrightarrow
R_{\rm throat},\text{shape}
\longleftrightarrow
\widetilde Q_n^{(w)},Q_n^{\rm net},Q_n,Q_\chi
\longleftrightarrow
p_{\rm bulk},\text{brane tension}
}.
\]

The resulting flow and conversion pattern enters the provisional orientation-dependent gravity-source map,

\[
\widetilde Q_g^{(s)}(r)
=\widetilde{\mathcal F}_r\!\left[
\widetilde Q_n^{(w)},Q_n^{\rm net},Q_n,Q_\chi,
\Gamma_{\rm return},\mathbf J_\chi,
\text{bulk and brane response}
\right].
\]

No form for \(\widetilde{\mathcal F}_r\) is assumed. The reflected solutions must derive covariance with reflected environmental data,

\[
\boxed{
\mathbf g_{\rm eff}
[\mathcal C_{\rm th}\mathcal T;
 \mathcal R_w\mathcal B_\infty]
=
\mathbf g_{\rm eff}
[\mathcal T;\mathcal B_\infty]
},
\]

and, for same-environment degeneracy, \(\mathcal R_w\mathcal B_\infty=\mathcal B_\infty\). Independently solved drain branches must have equal positive throughput even though

\[
\widetilde Q_n^{(w)}(+)
=-\widetilde Q_n^{(w)}(-).
\]

After reflection closure, the common coefficient \(\overline Q_g(r)\) remains radially signed. The ordinary local particle branch must separately establish \(\overline Q_g(r_{\rm local})>0\); screening, return overshoot, or a possible sign change at larger radius must remain visible. The arrows identify conceptual roles, not a one-way monotonic sequence. The support mode, geometry, signed transport, throughput, conversion, brane tension, bulk backpressure, and global return must be solved together.

**Candidate offset-driven conversion and entrainment.** Conditional on the symmetric slab profile above, a throat carries ordered material away from the local mid-surface and toward regions where the ordered state is less favorable. The resulting de-structuring, chemical-potential change, and pressure response may entrain surrounding brane material into the aperture. This provides a candidate mechanical driver for \(Q_\chi\) and, through the induced flow, \(Q_n\). It is not suction into empty space, need not affect \(Q_n\) and \(Q_\chi\) identically, and must be derived from the coupled throat solution. Reflection between \(+w\) and \(-w\) reverses the geometric orientation but not the scalar drain magnitude.

**Candidate critical-flow saturation.** A first idealized bound can be calculated for a fixed mouth radius and fixed upstream state. If the flow is steady, inviscid, isentropic, and nozzle-like, write \(P=Kn^\gamma\) with \(\gamma=5\), constituent mass parameter \(m_\star\), and

\[
c_s^2=\frac{\partial P}{\partial(m_\star n)}
=\frac{5K}{m_\star}n^4.
\]

At an ideal sonic throat, the standard polytropic critical point gives

\[
n_*=n_0\,3^{-1/4},
\qquad
c_{s*}=c_{s0}\,3^{-1/2},
\qquad
j_{n,\max}^{\rm ideal}=n_0c_{s0}\,3^{-3/4}.
\]

If the mouth is a three-ball of radius \(R_{\rm throat}\) in the brane, its trans-\(w\) cross-sectional three-volume is \(\mathcal A_3=4\pi R_{\rm throat}^3/3\), so

\[
Q_{n,\max}^{\rm ideal}
=\frac{4\pi R_{\rm throat}^3}{3}
n_0c_{s0}\,3^{-3/4}.
\]

This is a calculable upper envelope for \(Q_n\) under the stated assumptions, not an established property of the full condensate throat. Quantum pressure, order conversion, nonuniform geometry, critical-current instabilities, and the solved pressure boundary conditions may lower or replace it. The relation can be inverted to estimate the minimum radius needed to carry a target \(Q_n\), but it does not by itself bound \(Q_\chi\), stabilize the aperture, or impose a universal maximum particle mass.

The observed position \(\mathbf x_a=(x_a,y_a,z_a)\) is a projection onto the brane. Oppositely oriented throats may therefore have the same projected position while their nonlinear cores remain supported on opposite interfaces and separated in \(w\). In the interface-graph regime their local normal gap is

\[
d_w^{\rm core}
=H_{\rm br}+2h_t^{\rm tot}
-\ell_+^{\rm in}-\ell_-^{\rm in}.
\]

For the illustrative equal-amplitude outward pair, \(h_t^{\rm tot}=a_c\), so before inward core reach is included,

\[
d_w^{\rm core}\simeq H_{\rm br}+2a_c.
\]

Projected mouth footprints can therefore overlap without full four-dimensional core contact. The fields \(h_\pm\) cease to be sufficient if overhangs, necks, reconnection, pinch-off, or topology change occurs; the parent fields and core level sets must then supply the full geometry and an interface-graph-independent ambient core-set distance. Whether exact projected coincidence yields separated cores, contact, reconnection, annihilation, exclusion, or a neutral composite is an open short-range question developed in the [opposite-orientation derivation plan](opposite_orientation_throat_coincidence.md).

Because parent fields may have noncompact tails, a nonlinear core is not the literal support of every disturbed field. Choose a solved, invariant or operationally specified localized core diagnostic \(\mathcal I_a(\mathbf X)\)—for example order depletion, defect density, or conversion intensity—and define a thresholded core set

\[
\boxed{
\mathcal C_a(\epsilon)
=\left\{
\mathbf X:\mathcal I_a(\mathbf X)\ge\epsilon
\right\}
}.
\]

The admissible diagnostic and threshold range must make \(\mathcal C_a(\epsilon)\) compact, or at least closed and localized so that the distance infimum is attained. Individual core labels must themselves be outputs of the nonlinear solution. Two acceptable procedures are:

1. **Orientation-resolved diagnostics:** derive distinct physical \(\mathcal I_+\) and \(\mathcal I_-\), for example from an oriented topological or conversion-density field.
2. **Connected-component continuation:** use one total diagnostic \(\mathcal I(\mathbf X)\), identify its two thresholded connected components at large separation, and continue those components as \(R\) decreases. Contact or merger is the parameter value at which their closures first touch or the number or topology of components changes.

The symbols \(\mathcal C_+\) and \(\mathcal C_-\) may be used only while one of these procedures still identifies two physical cores. After merger, a decomposition of the total field into two analyst-assigned contributions is not permitted. While two sets remain identifiable, define the precontact separation

\[
d_{\rm core}(\epsilon)
=\inf_{\substack{
\mathbf Y_+\in\mathcal C_+(\epsilon)\\
\mathbf Y_-\in\mathcal C_-(\epsilon)}}
\left\|\mathbf Y_+-\mathbf Y_-\right\|_4,
\]

using the ambient Euclidean norm fixed in Section 1. For orientation-resolved sets, define actual contact by

\[
\boxed{
\mathcal C_+(\epsilon)
\cap
\mathcal C_-(\epsilon)
\ne\varnothing
}.
\]

For compact orientation-resolved sets, contact is equivalent to \(d_{\rm core}(\epsilon)=0\); without the required closure/localization, zero infimum alone is insufficient. Under connected-component continuation, use closure contact and the component-count or topology change as the merger criterion rather than forcing a post-merger \(+/-\) decomposition. Numerical conclusions must be robust over a reasonable range of \(\epsilon\). Neither contact nor distance selects the subsequent merger, reconnection, annihilation, or exclusion outcome.

The trapped support mode is not itself gravitational mass, and the model does not assume \(E_{\rm support}=M_{\rm active}^{\rm local}c^2\). The mode is structural support, not the direct gravitational source. The supported open aperture permits transport and conversion; the source-side flow, pressure, stress, conversion, return, and reservoir state must first produce an environmental gravity field independently of any later probe. A fixed calibration then defines the operational gravity observable, after which arbitrary test-throat response is derived separately.

Throat orientation reverses \(\widetilde Q_n^{(w)}\) when reflection partners have equal selected throughput. The quantity \(Q_n^{\rm net}\) is signed before branch selection; \(Q_n>0\) is the selected drain-branch throughput; and \(Q_\chi\) and \(Q_{\rm ap}\) are local nonnegative gross magnitudes. Same-sign gravity for reflected orientations and attractive radial sign are distinct requirements: reflection closure produces the signed \(\overline Q_g(r)\), and local attraction requires \(\overline Q_g(r_{\rm local})>0\).

The word **stationary** applies to the throat's geometry, time-averaged stress, and mean constituent-transport and order-conversion rates. Its supporting field remains dynamically oscillatory, for example

\[
\Phi(\mathbf x,t)
=\operatorname{Re}\!\left[\phi(\mathbf x)e^{-i\omega t}\right].
\]

The complete object may therefore be a periodic nonlinear bound state rather than a strictly static field configuration. Its stability may require Floquet analysis, not merely the Hessian of a static energy functional.

**Status — open hypothesis:** a finite oriented throat and its support mechanism are working ontological hypotheses. A self-consistent stationary geometry, spectrally localized trapped transverse standing mode, separate \(\widetilde Q_n^{(w)}\), \(Q_n^{\rm net}\), selected positive \(Q_n\), and \(Q_\chi\), reflection-covariant signed gravity map, and perturbative or Floquet stability have not yet been jointly derived.

## 7. Active, passive, and inertial mass

The word “mass” combines three observational roles that the ontology must initially keep separate.

### 7.1 Aperture processes, net gravity strength, and active mass

Before a drain branch is selected, the throat region has signed local net outward constituent flux

\[
Q_n^{\rm net}
=\int_{\mathcal A_a}
\mathbf J_n\cdot\hat{\mathbf n}_a\,d^3\Sigma.
\]

On a solved stationary outward-drain branch,

\[
Q_n\equiv Q_n^{\rm net}>0,
\]

is its positive net outward throughput. The same throat may also carry

\[
Q_\chi\equiv\text{gross ordered-to-bulk conversion rate at the throat},
\qquad Q_\chi\ge0.
\]

The corresponding coordinate-normal flux obeys the identity \(\widetilde Q_n^{(w)}=sQ_n^{\rm net}\) on the selected reference cut. Independently solved outward-drain partners have opposite coordinate-flux signs by orientation; their coordinate fluxes are equal and opposite only when the separately tested positive throughputs are equal. The quantity \(Q_\chi\) is a positive gross local drain magnitude or derived inertial-branch diagnostic. It does not include distributed return and is not itself a signed net conversion rate. In the fraction/relaxational branch, net conversion signs remain encoded by \(\Gamma_{\rm drain}\), \(\Gamma_{\rm return}\), and \(\Gamma_B\); an inertial branch must derive its own signed bookkeeping from its equations.

Their relationship depends on the selected \(\chi_B\) dynamics. When their distinction is not important, \(Q_{\rm ap}\) denotes the complete gross aperture process, but it is only shorthand and does not assert that \(Q_n\) and \(Q_\chi\) are numerically identical.

A sonic or other critical-flow calculation may bound the magnitude \(Q_n\) for a fixed throat geometry and upstream state. Turning that into a bound on \(M_{\rm active}^{\rm local}\) additionally requires the derived gravity map and the equilibrium relation that fixes \(R_{\rm throat}\). Conversely, once a provisional branch supplies those relations, an observed target mass can be used to infer the minimum mouth radius required by that branch. Such an inference is a calibration or necessary-condition test until the radius and source map are independently derived.

Before reflection parity is known, distributed return and the surrounding brane/bulk response define only the provisional raw coefficient

\[
\widetilde Q_g^{(s)}(r)
=\widetilde{\mathcal F}_r\!\left[
\widetilde Q_n^{(w)},Q_n^{\rm net},Q_n,Q_\chi,
\Gamma_{\rm return},\mathbf J_\chi,
\text{bulk and brane response}
\right].
\]

The final projected reduction must derive \(\widetilde{\mathcal F}_r\) and test reflection covariance with the environmental data transformed as well. If that test closes, introduce the orientation-even but still radially signed coefficient

\[
\boxed{
\overline Q_g^{(+)}(r)
=\overline Q_g^{(-)}(r)
\equiv\overline Q_g(r)\in\mathbb R
}
\]

under the properly reflected environment. Reflection evenness does not establish radial attraction. Let \(r_{\rm local}\) enclose the throat and its immediate conversion region but little distributed return, and require the ordinary local branch to pass the separate sign test

\[
\overline Q_g(r_{\rm local})>0.
\]

Only on that verified attractive local branch is the **operational definition** of local active mass

\[
M_{\rm active}^{\rm local}
\equiv
\frac{C_g\overline Q_g(r_{\rm local})}{G_{\rm eff}},
\qquad
\overline Q_g(r_{\rm local})>0.
\]

The three gravity normalizations have distinct operational roles:

\[
\begin{aligned}
\mathcal C_{\rm ref}
&:\ \text{parent-field units}\longrightarrow\text{acceleration units},\\
C_g
&:\ \text{radial/source geometric normalization in the inverse-square form},\\
G_{\rm eff}
&:\ \text{acceleration-source units}\longrightarrow\text{mass units}.
\end{aligned}
\]

These are not three freely adjustable physical couplings. Only their conventionally independent combinations may be calibrated once, after which they are frozen. None may be retuned by species, radius, or environment to enforce a desired sign or universality ratio.

A well-defined local active mass requires either a matching region

\[
R_{\rm throat}\ll r\ll L_{\rm return},
\]

in which \(\overline Q_g(r)\) is approximately constant and positive, or another scale-independent prescription for extracting the throat's local source coefficient. Here \(L_{\rm return}\) is the characteristic scale over which distributed return significantly changes the enclosed gravity source. If no such plateau exists, \(M_{\rm active}^{\rm local}\) must be obtained through a near-to-far matching procedure rather than chosen at an arbitrary observation radius.

At a general observation scale, retain a signed effective source:

\[
M_{\rm active}^{\rm eff}(r)
\equiv\frac{C_g\overline Q_g(r)}{G_{\rm eff}}
=f_{\rm ret}(r;\mathcal E)M_{\rm active}^{\rm local},
\]

where \(\mathcal E\) denotes any permitted environmental dependence. No form or sign restriction for \(f_{\rm ret}\) is assumed: it must expose screening, return overshoot, or any sign reversal rather than hide it.

After calibration, the corresponding spherically symmetric effective observable is parameterized as

\[
\mathbf g_{\rm eff}(r)
=-\frac{C_g\overline Q_g(r)}{r^2}\,\hat{\mathbf r}.
\]

This does not say that \(\overline Q_g\) is conventional mass. It says that effective active mass is the observer's signed parameterization of the enclosed gravity-source strength.

The return-profile map from \(M_{\rm active}^{\rm local}\) to \(M_{\rm active}^{\rm eff}(r)\) must be universal, environmentally calculable, or otherwise observationally acceptable; it cannot vary arbitrarily by throat species. If distributed return completely compensates localized drainage inside a sufficiently large region, then \(\overline Q_g(r)\) may approach zero and gravity may be screened beyond a return scale. Return overshoot or a sign reversal is also an open dynamical possibility that must remain explicit and observationally acceptable rather than being removed by an absolute value.

### 7.2 Passive gravitational mass

Passive gravitational mass is the leading isotropic, quasistatic coefficient determining how an arbitrary test throat responds to the externally generated, fixed-calibration gravity observable defined in Section 8:

\[
F_i^{\rm test}(\omega\to0)
=\mathcal M_{ij}^{\rm passive}(0)
g_{{\rm eff},j}+\cdots.
\]

More generally, the response may be tensorial, frequency dependent, and history dependent:

\[
F_i(\omega)
=\mathcal M^{\rm passive}_{ij}(\omega)g_{{\rm eff},j}(\omega)+\cdots.
\]

The scalar \(M_{\rm passive}\) must emerge in the low-frequency isotropic limit. It is a response property of the complete dressed throat, is not automatically equal to the active source strength, and must be derived from the full moving-defect response.

After the gravity sign and calibration conventions are frozen, the ordinary stable branch must have positive passive response. Writing the symmetric quasistatic part as \(\mathcal M_{(ij)}^{\rm passive}(0)\), require

\[
\boxed{
\xi_i\mathcal M_{(ij)}^{\rm passive}(0)\xi_j>0
\quad
\text{for every nonzero }\boldsymbol\xi
},
\]

or \(M_{\rm passive}>0\) in the isotropic limit. Universality of a wrong-sign or singular passive response is not acceptable.

The source-side environmental field is derived independently of the later test species. A single reference convention may set the units of \(\mathbf g_{\rm eff}\), but that calibration is then frozen; it may not absorb species dependence and make universal free fall tautological.

### 7.3 Composite inertial response

Inertial mass governs the resistance of the complete dressed throat to acceleration:

\[
\mathbf F_{\rm ext}=M_{\rm inertial}\mathbf a+\cdots.
\]

It may contain several coupled contributions.

**Standing-wave reconfiguration.** The trapped transverse brane-shear standing mode—the model's light or photon mode—need not accelerate as a rigid object. As the throat changes velocity:

- its counterpropagating components may acquire different timing;
- phases and nodes must readjust;
- the standing-pattern frequency or round-trip cycle may change;
- the resonant throat geometry may change;
- momentum must be transferred into the internal wave.

The underlying transverse wave speed \(c_\gamma\) need not slow. What may reorganize is the standing pattern, its internal clock, or its resonant cycle.

**Near-field reorganization.** A moving throat carries an inflow pattern, pressure and density disturbances, brane deformation, an order-conversion region, and possibly an upstream compression and downstream wake. Acceleration requires this dressed configuration to be rebuilt. That rebuilding can produce a reversible added-inertia-like response even if no microscopic constituent carries the observed particle mass.

The theory must distinguish reversible near-field reorganization, which contributes to inertia, from outgoing wave production, which produces radiation or damping.

The correct inertia formalism depends on the \(\chi_B\) branch. For a closed reversible branch,

\[
M_{ij}^{\rm eff}
=\left.\frac{\partial P_i}{\partial v_j}\right|_{\mathbf v=0},
\]

where \(P_i\) includes internal standing-wave momentum, near-field medium momentum, conversion-field momentum, and throat-shape and brane-deformation contributions.

On the ordinary stable branch, the symmetric zero-frequency inertial response must be positive definite:

\[
\boxed{
\xi_iM_{(ij)}^{\rm eff}(0)\xi_j>0
\quad
\text{for every nonzero }\boldsymbol\xi
}.
\]

Zero or negative inertial directions signal a marginal or unstable branch rather than an acceptable universal mass.

For a dissipative branch, a conserved canonical momentum cannot be assumed. The response to a small projected throat displacement \(\delta\mathbf x_a\) is instead frequency dependent:

\[
F_i(\omega)
=\left[-\omega^2M_{ij}(\omega)-i\omega\Gamma_{ij}(\omega)+\cdots\right]\delta x_{a,j}(\omega),
\]

where \(M_{ij}(\omega)\) is the dynamical inertial response and \(\Gamma_{ij}(\omega)\) is damping. Its ledger must include entropy, heat, or explicit reservoir variables, all representing degrees of freedom internal to the same substrate.

The retarded response of an ordinary unexcited internal reservoir must also be passive: for every small periodic perturbation, the cycle-averaged external power entering the complete throat-plus-reservoir system must satisfy

\[
\boxed{
\left\langle
\mathbf F_{\rm ext}\cdot\dot{\mathbf x}_a
\right\rangle_{\rm cycle}
\ge0
}.
\]

Equality is allowed for a lossless channel. A response that extracts net work from an unexcited reservoir is unacceptable unless an explicitly prepared energetic state and its depletion are included in the ledger.

For an explicitly mixed inertial-dissipative branch, both reversible momentum storage and the retarded damping response must be derived from the same closed system with its internal reservoir variables. The conservative and dissipative pieces may coexist, but neither may be inserted without the associated material, momentum, energy, and entropy accounting.

**Status:** the possible components of inertia are identified, but neither the reversible momentum, positive-definite ordinary response, dissipative kernel, nor retarded passivity has been derived from a supported-throat solution.

### 7.4 Radiation thresholds and uniform motion

A stationary comoving distortion around a uniformly moving throat may be acceptable, but the drag test depends on the \(\chi_B\) branch.

**Reversible branch.** On a homogeneous isotropic background, the Landau-like kinematic emission threshold is determined by the complete conservative coupled spectrum,

\[
v_c=\min_a\inf_{k>0}\frac{\omega_a(k)}{k},
\]

where \(a\) runs over every branch coupled to the throat: transverse, longitudinal, bulk-sound, flexural, electric/order, and conversion modes.

For an anisotropic background or defect environment, the threshold is directional. With motion along \(\hat{\mathbf v}\), the corresponding kinematic form is

\[
v_c(\hat{\mathbf v})
=\inf_{\substack{a,\mathbf k\\
\mathbf k\cdot\hat{\mathbf v}>0}}
\frac{\omega_a(\mathbf k)}
{\mathbf k\cdot\hat{\mathbf v}},
\]

equivalently testing the resonance condition \(\omega_a(\mathbf k)=\mathbf k\cdot\mathbf v\). The isotropic expression must not be used after a throat or environment breaks rotational symmetry without checking this directional generalization.

If a branch has \(\omega\propto k^2\), then \(\omega/k\rightarrow0\) at long wavelength. Avoiding drag would then require symmetry suppression or vanishing coupling rather than merely a low velocity.

**Dissipative or mixed branch.** Absence of propagating radiation is insufficient to guarantee drag-free uniform motion. The zero- and low-frequency friction, order-relaxation, and conversion kernels must vanish or remain below observational constraints. Schematically, the additional test is

\[
\Gamma_{ij}(\omega\rightarrow0)\rightarrow0
\]

or sufficient suppression for uniform translation. Exact zero damping is not assumed.

### 7.5 Velocity-dependent aperture processes

Lateral motion may alter constituent transport, order conversion, and the effective net gravity source. For an isotropic throat, reversing the velocity should not ordinarily change either scalar magnitude, so the provisional low-speed forms are

\[
Q_A(v)
=Q_{A0}\left[1+\alpha_A\frac{v^2}{c_*^2}+O(v^4)\right],
\qquad A\in\{n,\chi\}.
\]

The notation \(Q_{\rm ap}(v)\) may summarize the combined aperture process without identifying the two coefficients. The selected drain branch must continue to satisfy \(Q_n(v)=Q_n^{\rm net}(v)>0\). Distributed return then determines the signed \(\overline Q_g(r;v)\). Possible causes include upstream pressure, throat widening, support-mode reconfiguration, modified conversion, or coupling between lateral and normal flow. Any increase, decrease, screening, or reversal of the signed effective source must remain visible.

**Status:** velocity-dependent constituent transport, order conversion, and net gravity strength are hypotheses and outputs sought from the moving-throat family, not established results. Any dependence on velocity relative to the substrate is a preferred-medium effect that must be calculated and constrained.

### 7.6 The equivalence requirements

Active, passive, and inertial mass need not arise from identical microscopic mechanisms. The relevant family of stable supported-throat solutions must instead produce species-universal ratios in a defined low-energy operational regime:

\[
\left.
\frac{M_{\rm passive}}{M_{\rm inertial}}
\right|_{\substack{
v\to0,\ \omega\to0\\
\text{weak external field}\\
\text{homogeneous environment}}}
=\eta_{\rm PI},
\]

and

\[
\left.
\frac{M_{\rm active}^{\rm local}}{M_{\rm inertial}}
\right|_{\substack{
v\to0,\ \omega\to0\\
R_{\rm throat}\ll r\ll L_{\rm return}\\
\text{or equivalent matching}}}
=\eta_{\rm AI}.
\]

The constants may be normalized only after units and source/response conventions are fixed. Universality of passive/inertial response is the analog of universal free fall. Universality of local active/inertial response is a source-strength reciprocity or source-universality requirement. Neither result alone establishes the full geometric equivalence principle of GR. At finite velocity, frequency, field strength, or environmental inhomogeneity, the tensorial and dispersive response kernels remain the relevant quantities.

These ratios are tested only after the ordinary branch has passed the separate positive active-source, passive-response, inertial-response, and retarded-passivity conditions. Equality of equally wrong-signed coefficients is not universality.

Such low-energy universality would be an emergent property of the same complete object: local active mass arises from the matched source-region gravity strength, passive mass from response to \(\mathbf g_{\rm eff}\), and inertia from the coupled dynamical response of the support mode, geometry, near fields, and conversion sector. The separate return profile \(f_{\rm ret}(r;\mathcal E)\) must then map local active mass into scale-dependent effective mass in a universal or calculable way.

**Status:** the three-role distinction is an ontological correction; the local active source map, passive response, inertial kernel, both universality ratios, and the return-profile map to effective active mass remain open outputs.

## 8. Gravity

**Postulated mechanism.** Gravity is attributed to the quasi-steady material flow, pressure, stress, conversion, and return configuration associated with an open throat. Signed coordinate flux \(\widetilde Q_n^{(w)}\), signed local net outward flux \(Q_n^{\rm net}\), selected positive throughput \(Q_n\), and nonnegative gross local conversion drain \(Q_\chi\) are distinct quantities; \(Q_{\rm ap}\) is only an umbrella shorthand for gross aperture processes.

**Open source-side reduction.** Reserve \(\mathbf v_{\rm med}\) for four-dimensional material velocity. Before selecting any later test probe, define a source-side environmental field

\[
\boldsymbol{\mathfrak g}_{\rm src}
=\mathcal G_{\rm src}\!\left[
\mathbf J_n,
\mathbf v_{\rm med},
D_t\mathbf v_{\rm med},
n,P,\mu,\chi_B,
\boldsymbol\tau_{\rm order},
\Gamma_B,\mathbf J_\chi,
\text{return and reservoir state};
\mathcal B_\infty
\right].
\]

The list contains candidate parent variables; it does not require every one to survive the final reduction, and no formula for \(\mathcal G_{\rm src}\) is assumed. Crucially, no test-throat response is an input to this source-side functional.

**Fixed operational calibration.** One fixed reference convention sets units through

\[
\boxed{
\mathbf g_{\rm eff}
=\mathcal C_{\rm ref}\,
\boldsymbol{\mathfrak g}_{\rm src}
}.
\]

The calibration \(\mathcal C_{\rm ref}\) may be established using one reference throat or another operational standard, but once fixed it cannot depend on the species of a later test body. It maps the derived parent-field quantity into acceleration units; \(C_g\) separately records the chosen radial/source geometric normalization, and \(G_{\rm eff}\) converts the calibrated source coefficient into mass units. Only conventionally independent combinations may be fixed once. The source-side gravity field is derived from the environmental medium configuration independently of the later probe species. Passive mass is then a separately derived response coefficient of an arbitrary throat to that calibrated field. The calibration may set units, but it may not make universal free fall tautological by absorbing species-dependent response. “Gravity is inflow” identifies the proposed underlying material mechanism; it does not make raw \(\mathbf v_{\rm med}\) numerically identical to \(\mathbf g_{\rm eff}\).

**Target reflection covariance and radial sign.** Electric orientation parity and radial gravitational sign are separate tests. Let \(\mathcal R_w\) act geometrically on the parent fields and environmental data, and let \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\) denote the complete candidate throat operation. The general target is

\[
\boxed{
\mathbf g_{\rm eff}
[\mathcal C_{\rm th}\mathcal T;
 \mathcal R_w\mathcal B_\infty]
=
\mathbf g_{\rm eff}
[\mathcal T;\mathcal B_\infty]
}.
\]

The corresponding source target is an orientation-even but radially signed coefficient,

\[
\boxed{
\overline Q_g^{(+)}(r)
=\overline Q_g^{(-)}(r)
},
\qquad
\overline Q_g(r)\in\mathbb R.
\]

For the two orientations to be degenerate in the same vacuum environment, one additionally requires

\[
\boxed{
\mathcal R_w\mathcal B_\infty
=\mathcal B_\infty
}.
\]

If the two bulk-facing reservoirs differ, an orientation-dependent environmental splitting is a calculable environmental effect and may be an unacceptable charge-correlated gravitational signal; it is not automatically a failure of the intrinsic candidate throat operation. In the chosen radial convention,

\[
\mathbf g_{\rm eff}(r)
=-\frac{C_g\overline Q_g(r)}{r^2}\hat{\mathbf r}.
\]

The ordinary local particle branch must independently satisfy \(\overline Q_g(r_{\rm local})>0\). Distributed return may drive \(\overline Q_g(r)\) toward zero, produce screening or overshoot, or possibly reverse its sign. Such behavior must be derived and constrained, not hidden by an absolute-value convention.

Gravity is not identified with light. The leading local compressional response propagates on the bulk sound branch with characteristic speed \(c_s\). The complete time-dependent gravity response may also contain order-conversion, interface, return, reservoir, or dissipative relaxation modes. Its poles must be derived from the coupled response operator; schematically,

\[
g_i(\omega,\mathbf k)
=\mathcal R_{ij}^{g}(\omega,\mathbf k)
\mathcal S_j(\omega,\mathbf k).
\]

The intended ordinary local and observational far-field branch is attractive and inverse-square on the brane. Existing calculations support the scaling only under effective sink/return assumptions. The model must still derive \(\widetilde{\mathcal F}_r\), \(\mathcal G_{\rm src}\), \(\mathcal C_{\rm ref}\), reflection covariance, and the signed \(\overline Q_g(r)\). It must also show that a supported throat produces stable transport and conversion and that the return distribution gives an acceptable \(M_{\rm active}^{\rm eff}(r)\), including any screening, overshoot, or sign-change scale.

Gravitomagnetism is the vorticity or frame-dragging component of this gravity flow. It is distinct from electromagnetic magnetism, even though both may be described informally as swirl.

**Status:** gravity as throat-driven material response is postulated; inverse-square scaling and leading attraction are derived only within prior effective sink/return assumptions. The raw source map \(\widetilde{\mathcal F}_r\), source-side reduction \(\mathcal G_{\rm src}\), fixed calibration \(\mathcal C_{\rm ref}\), reflected-environment covariance, orientation-even signed \(\overline Q_g(r)\), local sign, response poles, return-profile universality, and large-scale behavior are **not yet derived**.

## 9. Light

Light is an in-plane shear wave of the brane. If \(\mathbf u\) is the brane displacement, its transverse part \(\mathbf u_T\) satisfies

\[
\nabla\cdot\mathbf u_T=0.
\]

The supplied quadratic brane action supports two transverse modes in three brane dimensions. Their characteristic speed is

\[
c_\gamma^2=\frac{\mu_R}{\rho_{\rm br}},
\]

where \(\mu_R\) is the effective brane shear modulus and \(\rho_{\rm br}\) is its effective inertia density.

The two-polarization count is conditional on the assumed action, isotropic inertia, and the relevant structural premises. The model has not yet derived the shear modulus from the underlying medium.

Freely propagating light and the throat-support mode must be distinguished. Free light is a background-brane \(\mathbf u_T\) excitation. The support mode is a spectrally normalizable bound state or acceptably long-lived resonance of the complete variable-coefficient transverse operator; it is not a second photon substance and is not assumed to propagate as a shear wave through fully de-structured bulk material. Any \(\Omega_{\rm support}\) is a diagnostic energy/stress-localization region extracted from that solved mode, not a hard PDE domain imposed by \(\mu_R>0\).

Whether transverse light is sufficiently confined to the brane is an interface question. Uniform transverse coupling to the bulk vanishes under the supplied assumptions, but the physically relevant nonuniform, gradient-driven coupling remains to be completed. The model must show that bulk leakage, longitudinal conversion, or defect-induced mixing does not spoil ordinary transverse propagation. It must also derive how an accelerating oriented throat excites the two far-zone transverse polarizations through the same coupling structure as static electricity and magnetism.

**Status:** the conditional two-transverse-mode count follows from the supplied quadratic action; the origin of the shear modulus, nonuniform interface isolation, and propagation through defect backgrounds remain open.

## 10. The stray longitudinal mode

The same brane can also carry an in-plane longitudinal displacement \(u_L\). In an ordinary symmetric elastic response, the medium consistently supports both:

- two transverse modes; and
- one longitudinal mode.

The longitudinal mode is not a mathematical inconsistency and is not something the model intends to tune away merely to imitate exact Maxwell theory. It is a characterized physical departure from a transverse-only electromagnetic field.

It is essential to keep three objects separate:

- \(u_L\): in-plane longitudinal displacement;
- \(h_+,h_-\), or equivalently \(h_m,h_t\): the two interface displacements and their mid-surface and thickness combinations in the normal \(w\) direction;
- \(s=\pm1\): the orientation of a throat and the sign of charge.

The longitudinal mode is not charge and does not become the throat's \(\pm w\) orientation. Its eventual physical role depends on the brane–bulk interface law: it may remain bound, radiate into bulk sound, acquire its own long-wavelength dispersion, or mix with transverse waves near defects. Moving and accelerating throat calculations must compute any emitted longitudinal power rather than assume that only the transverse light sector radiates.

**Status:** simultaneous transverse and longitudinal mode support is established for the supplied homogeneous quadratic response. The physical dispersion, brane–bulk leakage, and defect-induced mixing of the longitudinal branch remain open.

## 11. Electric charge as an oriented defect

Electric charge is identified with the orientation of a puncture:

\[
s=+1 \quad\text{for a \(+w\) throat},
\qquad
s=-1 \quad\text{for a \(-w\) throat}.
\]

The throat is not modeled as a freely adjustable ordinary scalar source. It imposes an oriented geometric and material boundary condition on the finite-thickness brane. Schematically, a signed material disturbance \(\sigma\) obeys

\[
\sigma\big|_{\partial\Omega_a}=s_a\sigma_0
\]

at the core of throat \(a\), and the external brane configuration minimizes its material energy subject to that fixed oriented defect. The quantity \(\sigma_0\) is a core boundary amplitude or branch condition, not an intrinsic scalar charge magnitude.

### 11.1 The two-interface electric sector

For a finite slab, the normal embedding variables are \(h_+(\mathbf x)\) and \(h_-(\mathbf x)\), or equivalently \(h_m\) and \(h_t\). Any single coordinate \(h\) is shorthand for a derived eigencombination of these variables. The physical electric disturbance should not be identified with the absolute value of any bare height alone. If the parent medium is translation invariant in \(w\), uniformly shifting the whole brane cannot by itself cost energy.

The electric energy is therefore expected to depend on a derived signed stability variable

\[
\sigma=\sigma[h_+,h_-,\chi_B],
\]

which may contain:

- gradients or curvature of the interface or parity coordinates;
- asymmetry between the two brane interfaces;
- local thinning or thickening;
- normal strain through the slab;
- departure of \(\chi_B\) from its stable profile;
- a signed moment of the order distribution through \(w\).

Under reflection,

\[
w\rightarrow-w,
\qquad
s\rightarrow-s,
\qquad
\sigma\rightarrow-\sigma,
\]

while the stored disturbance energy remains symmetric:

\[
W_{\rm stress}(-\sigma)=W_{\rm stress}(\sigma).
\]

A schematic local stored-stress energy is

\[
W_{\rm stress}(\sigma)
=\frac{K_\sigma}{2}\sigma^2
+\frac{\lambda_\sigma}{4}\sigma^4+\cdots,
\qquad K_\sigma>0.
\]

In an overlap region, \(\sigma=\sigma_1+\sigma_2\). Same orientations increase \(|\sigma|\) and the stored stress energy, while opposite orientations partially cancel it. Positive even powers alone do not destabilize the ordered material state.

Actual approach to de-structuring requires a coupled material potential,

\[
V_{\rm mat}(\chi_B,\sigma)
=V_0(\chi_B)
+\frac12K(\chi_B)\sigma^2
+\frac14\Lambda(\chi_B)\sigma^4+\cdots,
\]

whose \(\sigma\)-dependence changes the stability of the ordered \(\chi_B\) branch. Possible signals include vanishing ordered-branch curvature,

\[
\left.
\frac{\partial^2V_{\rm mat}}{\partial\chi_B^2}
\right|_{\chi_{\rm ordered}}
\rightarrow0,
\]

a falling barrier between ordered and de-structured states, or equality of their effective free energies. No particular \(K(\chi_B)\), \(\Lambda(\chi_B)\), or threshold is assumed. Same-oriented overlap may both store more stress energy and reduce the ordered material's stability margin; these are related but distinct effects.

This local stability energy does not by itself guarantee a long-range electric field. If \(\sigma\) were an independent field with a local quadratic restoring term, its response would ordinarily be gapped and screened. The actual eigenvectors may mix \(h_m,h_t,\delta\chi_B\), density, phase, longitudinal, bulk, and reservoir variables. After all constraints and coupled fields are included, \(\mathcal H_{\rm odd}\) and \(\mathcal H_{\rm even}\) denote the complete equilibrium static Hessian blocks, while \(\mathcal O_{R,{\rm odd}}\) and \(\mathcal O_{R,{\rm even}}\) denote the dynamic retarded-response blocks. Each may contain multiple eigenbranches. The bare \(h_m\) and \(h_t\) coordinates are only leading odd and even interface projections, not guaranteed complete eigenmodes. For a conservative or genuine equilibrium Coulomb potential, the desired static spectral outputs include at least one odd Hessian eigenvalue

\[
\lambda_E(k)
\in\operatorname{spec}\mathcal H_{\rm odd}(k),
\qquad
\lambda_E(k)
\stackrel{?}{=}
\kappa_Ek^2+O(k^4),
\qquad
\kappa_E>0,
\]

for a Coulomb-carrying odd eigenbranch with nonzero throat coupling, and at least one primarily thickness-like even eigenvalue

\[
\lambda_t(k)
\in\operatorname{spec}\mathcal H_{\rm even}(k),
\qquad
\lambda_t(k)
\stackrel{?}{=}
\mu_t^2+\kappa_tk^2+O(k^4),
\qquad
\mu_t^2>0,
\]

for a screened thickness-like eigenbranch. These are **target relations** for selected equilibrium eigenvalues, not identities for the full blocks. The retarded blocks must separately determine poles, damping, phase lag, and radiation. Parity alone does not establish range: an even variable may mix with another gapless field, multiple odd or even branches may exist, and stable thickness selection does not guarantee that every even eigenmode is screened. Every additional long-range mode and its throat coupling must be calculated and constrained.

The working requirement is that the nonlinear two-interface/order core selects the physical orientation branch and interaction sign, while an eigenmode with the required static stiffness carries the long-range tail. Schematically,

\[
E_{\rm far}
=\frac{\kappa_E}{2}\int d^3x\,|\nabla\phi_E|^2.
\]

Equivalently,

\[
E_{\rm far}
=\frac12
\int\frac{d^3k}{(2\pi)^3}
\kappa_Ek^2
|\phi_E(\mathbf k)|^2,
\]

so the target static Green function is

\[
G_E(0,k)\sim\frac{1}{\kappa_Ek^2},
\qquad
G_E(R)\sim\frac{1}{4\pi\kappa_ER}.
\]

A Coulomb tail therefore requires leading \(k^2\) static stiffness; gaplessness alone is insufficient. With ordinary local inertia that branch normally has \(\omega_E\simeq c_Ek\). A flexural branch with \(k^4\) static stiffness can have \(\omega\propto k^2\), but it has a different static Green function. The two branches must not be identified merely because both are gapless. A nonstandard kinetic operator may alter the dynamic scaling only if it is explicitly derived.

The symbol \(\phi_E\) denotes whichever eigencombination diagonalization identifies as the healthy Coulomb carrier with the required \(k^2\) static stiffness, if one exists; it is not an additional fundamental field.

### 11.2 Target derived electric strength and required interaction

The intended force mechanism is not an ordinary healthy scalar field linearly sourced by charge. After such a field is solved on shell, a conventional linear source coupling generally favors like-source attraction. That is the wrong mechanism here.

Instead:

\[
s=\pm1
\longrightarrow
\text{oriented throat boundary condition}
\longrightarrow
\text{signed brane strain/order disturbance}
\longrightarrow
\text{stability energy}
\longrightarrow
\text{force}.
\]

Same-oriented throat disturbances must compound the departure from the preferred brane state and raise the configurational energy. Opposite-oriented disturbances must partially cancel and restore the brane toward that state. The throats are not passive marbles pulled together by a deformable sheet; they are defects imposing material states whose compatibility depends on orientation.

Charge sign belongs to an individual oriented throat. Interaction strength is to be derived from the shared finite-thickness brane boundary-value problem rather than assigned as an independent scalar property by hand. Whether the leading long-distance coefficient factorizes into one-throat couplings is an output. A local healthy \(k^2\)-stiffness Coulomb carrier and a conserved oriented current favor asymptotic factorization, while nonlinear core, polarization, environmental, finite-size, and many-body corrections may remain relational.

Define the target coefficient of the **reflection-odd Coulomb monopole channel only** by

\[
C_{E,ab}(R)
=g_ag_b+\delta C_{ab}^{\rm odd}(R),
\]

where \(g_a\) is a possible derived nonnegative leading coupling magnitude of throat \(a\) to the Coulomb-carrying electric eigenmode. The correction is restricted to finite-size, polarization, environmental, or nonlinear corrections within that same odd monopole channel. It must not absorb reflection-even thickness forces or nonlinear core interactions.

For a conservative or genuine equilibrium branch, the target full pair-energy organization is

\[
E_{ab}(R)
=s_as_b\frac{g_ag_b+\delta C_{ab}^{\rm odd}(R)}{R}
+E_{ab}^{\rm even}(R)
+E_{ab}^{\rm mix}(R)
+E_{\rm core}^{\rm nl}(R).
\]

On an exactly reflection-symmetric background, \(E_{ab}^{\rm mix}\) should vanish at quadratic order, but even and nonlinear terms may remain, especially near \(R\sim H_{\rm br}\). For a conservative intrinsically periodic pair, the corresponding scalar object—if one exists—must be a derived cycle-averaged Hamiltonian, Routhian, quasienergy, or averaged action in a stated fixed-action, fixed-frequency, fixed-relative-phase, or other conjugate ensemble. Otherwise its mean force must come from cycle-averaged Noether stress and momentum flux. In a driven dissipative stationary branch, the same decomposition may organize response channels, but the force must be obtained from stress, momentum flux, conversion, and reservoir response unless a conservative potential is independently established.

For the odd equilibrium channel, let \(A_{ab}(R)=g_ag_b+\delta C_{ab}^{\rm odd}(R)\). Its radial force contribution is

\[
F_{R,ab}^{\rm odd}
=s_as_b
\left[
\frac{A_{ab}(R)}{R^2}
-\frac{A_{ab}'(R)}{R}
\right].
\]

Here positive \(F_R\) means force toward increasing separation. Thus a positive coefficient alone does not fix the force sign at finite \(R\). For \(g_ag_b>0\), the leading far-field like-repels/opposite-attracts target additionally requires

\[
\frac{\delta C_{ab}^{\rm odd}(R)}{g_ag_b}\to0,
\qquad
\frac{R\,\partial_R\delta C_{ab}^{\rm odd}(R)}{g_ag_b}\to0
\]

in the asymptotic regime. More generally,

\[
g_a=g\!\left(H_a,\mathcal B_a,\text{shared constitutive data}\right),
\]

where \(H_a\) is the local slab thickness and \(\mathcal B_a\) denotes relevant core boundary data.

For ordinary particles on a nearly homogeneous slab,

\[
H_a\simeq H_b\simeq H_{\rm br}
\]

is helpful but insufficient by itself to make charge universal. The supported-throat solution must also show that \(\mathcal B_a\) lies on a universal, quantized, or topology-controlled orientation branch, or that the leading Coulomb-carrier coupling is insensitive to the internal support-mode state. Near-universality requires

\[
g_a\simeq g_0,
\qquad g_0\ge0
\]

for ordinary unit-charge throats. Different standing-wave configurations must not cause unacceptable species-dependent electric coupling unless that nonuniversality is retained as a falsifiable prediction. Thickness, core, environmental, and many-body corrections must all be calculated and constrained.

The **leading far-field odd channel** gives repulsion for \(s_as_b=+1\) and attraction for \(s_as_b=-1\) only if the coefficient and derivative conditions above are derived. In an equilibrium branch, this requires a bounded-below, ghost-free two-interface/order boundary-value problem whose physical ensemble has the required interaction-energy ordering after one-body self terms are subtracted. Fixed-value, fixed-source, fixed-flux, mixed, or topological boundary conditions need not give the same sign. In a conservative periodic branch, the sign must follow from the correctly constrained cycle-averaged functional or averaged stress force. In a driven branch, it must come from the complete stress and momentum-flux calculation. No sign is guaranteed by calling a height coordinate a scalar or by choosing a source sign.

The nonlinear throat must fix the core boundary branch and normalization and supply the solved source profile seen by the linear far field. The exterior may still linearize to a \(1/R\) potential and a \(1/R^2\) force; strong nonlinearity need not persist to arbitrarily large distance.

Finite slab thickness introduces a further distinction between far-field cancellation and core contact. For a same-species pair with opposite orientation, the reflection-odd disturbance can cancel at projected coincidence while reflection-even thickness disturbance and internal support configurations remain nonzero. For periodic throats this means the cycle-averaged odd source, or the instantaneous source only when \(\mathcal C_{\rm th}\) fixes the required relative phase; arbitrary relative phase may leave transient odd response and sidebands. In the graph regime,

\[
d_w^{\rm core}
=H_{\rm br}+2h_t^{\rm tot}
-\ell_+^{\rm in}-\ell_-^{\rm in}.
\]

The illustrative equal-amplitude outward pair gives \(d_w^{\rm core}\simeq H_{\rm br}+2a_c\) before inward reach. Thus vanishing leading electric monopole does not imply disappearance or contact of the material cores. The finite-thickness Green matrix must determine the \(R\sim H_{\rm br}\) crossover. If the interfaces cease to be graphs, full parent fields, core level sets, and the nonlinear set distance must replace \(h_\pm\). The exact-coincidence outcome remains open.

If the leading asymptotic interaction factorizes, the convention is

\[
g_a\ge0,
\qquad
q_a=s_ag_a.
\]

The orientation \(s_a=\pm1\) carries the entire electric sign. The derived quantity \(g_a\) is conventionally defined as a nonnegative coupling magnitude. Any apparent sign in the derived coupling should therefore be absorbed into the orientation branch or treated as evidence that the proposed branch assignment is inconsistent.

Nonfactorizing corrections then represent departures from ideal Coulomb universality rather than preventing the definition of a leading effective charge. Leading factorization is favored by, but does not follow merely from naming, a local \(k^2\)-stiffness Coulomb carrier and conserved oriented current.

**Status:** oriented throat sign and the boundary-value ontology are postulated; a signed far-field form is conditionally motivated. The healthy \(k^2\)-stiffness eigenmode, leading coupling \(g_a\), odd-channel correction \(\delta C_{ab}^{\rm odd}\) and its derivative, equilibrium or stress-derived force sign, even/mixed/core interactions, core-data universality, and degree of same-brane universality remain open.

### 11.3 Additive charge and conservation

The binary label \(s=\pm1\) does not automatically provide additive charge conservation or a Gauss law. The global sum also requires the brane's postulated co-orientability or an equivalent oriented intersection construction. A candidate route is

\[
N_s=\sum_a s_a.
\]

More generally, the **target topological definition** is

\[
N_s
=\sum_a
I(\text{throat}_a,\text{oriented brane}),
\]

where the intersection number \(I\) remains to be defined. If the brane becomes non-co-orientable or loses a global normal bundle, \(N_s\) must be reformulated rather than used unchanged.

This net oriented throat number is logically distinct from the nonnegative magnitude with which the brane couples a throat to the long-range electric eigenmode. If the leading asymptotic coupling factorizes and is sufficiently universal, one may define the conventional long-distance additive charge,

\[
Q_{\rm eff}=g_0N_s.
\]

A same-species throat–antithroat pair has zero net oriented count, so charge conservation permits annihilation or relaxation into neutral propagating and bound modes; it does not require instantaneous annihilation on contact. Oppositely oriented objects belonging to different internal or composite species may instead form a neutral bound state or retain additional conserved structure. The model must derive these distinctions rather than identifying every positive orientation with the antiparticle of every negative one.

To establish this, the model must define:

- the topological or geometric quantity being counted;
- why continuous evolution conserves it;
- how oppositely oriented throats can be created or annihilated;
- how the count appears as a projected Gauss-like flux law.

**Status:** conservation of oriented throat number, leading asymptotic factorization, and universality of \(g_a\) are distinct open requirements. Neither an additive conventional charge nor a Gauss-like projected law has yet been derived.

## 12. Magnetism and the moving oriented throat

Electromagnetic magnetism is the sign-dependent transverse response of the same oriented throat when it moves. A moving \(\pm w\) puncture is expected to source a response of schematic form

\[
\mathbf J_T\propto s\,\mathbf V.
\]

Electricity and magnetism are only structurally related until their static and moving responses follow from one conserved oriented throat current and one shared brane coupling structure. A target continuity law is

\[
\partial_t\rho_s+\nabla\cdot\mathbf J_s=0.
\]

The model may not independently choose the leading electric coupling magnitude \(g_a\), magnetic moving-source normalization, relative sign, velocity coupling, and light-sector response. Deriving them from one current and one coupling structure is the actual electromagnetic unification test. If the leading asymptotic electric interaction factorizes, the same derived nonnegative \(g_a\) magnitudes must enter the moving current, with electric sign carried by \(s_a\); nonlinear or environmental corrections remain separately testable.

A moving throat can also produce two other effects that must not be called magnetism:

1. **Electromagnetic magnetism:** the sign-dependent transverse response proportional to \(s\mathbf V\).
2. **Inertial medium response:** sign-independent rebuilding or compression of the dressed substrate.
3. **Velocity-dependent gravity:** any change in \(\widetilde Q_n^{(w)}(v)\), \(Q_n^{\rm net}(v)\), selected \(Q_n(v)\), \(Q_\chi(v)\), the raw gravity source, or the resulting signed \(\overline Q_g(r;v)\) after reflection closure.

All three may coexist, but they are different observables and require separate bookkeeping.

The direct moving-throat work has supported the expected magnetic tensor structure, distance scaling, and velocity order at a structural level. Its normalization and sign remain tied to the unresolved throat interior and electric boundary condition.

This is not yet exact Maxwell theory. The present quadratic system does not generate a first-class \(U(1)\) Gauss constraint, and the magnetic variable has a characterized time-reversal departure from Maxwell's magnetic field. The honest target is a Maxwell-like analog with explicit additional structure, not a hidden claim of exact electromagnetism.

**Status:** the sign-dependent current form and magnetic tensor/falloff structure are conditionally supported; the shared electric–magnetic normalization, physical sign, conserved oriented current, and complete moving-throat response remain open.

### 12.1 Dynamic electromagnetic closure target

For a slowly moving and accelerating oriented throat, one conserved oriented current and one derived coupling structure must control

\[
\begin{array}{ll}
O(v^0):&\text{the quasistatic electric response},\\[1mm]
O(v):&\text{the magnetostatic transverse response},\\[1mm]
O(a):&\text{outgoing transverse light radiation}.
\end{array}
\]

The continuity target

\[
\partial_t\rho_s+
\nabla\cdot\mathbf J_s=0
\]

is necessary but not sufficient. The parent coupling must show how acceleration transfers energy and momentum into both \(\mathbf u_T\) polarizations. If leading electric factorization is derived, the same nonnegative one-throat coupling \(g_a\) must enter the static electric asymptotic coefficient, the moving magnetic-current normalization, and the transverse radiation amplitude, with all electric sign carried by \(s_a\).

The moving and accelerating throat calculation must also determine the powers

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

emitted into the two transverse light modes, the odd normal/electric mode, the even thickness mode, the in-plane longitudinal mode, bulk density or sound modes, and order-conversion or reservoir modes. These extra channels are calculable predictions; none may be removed by an independently tuned coefficient.

**Status — target relation:** static, magnetic, and radiative normalization has **not yet been unified or derived**.

## 13. Brane–bulk conversion and leakage

Total constituent density obeys the exact four-dimensional conservation law defined in Section 2. The order-sector evolution must remain branch specific.

In the fraction/relaxational branch, where \(\chi_B\) itself is a coarse-grained ordered material fraction, the constitutive evolution equation for \(\chi_Bn\) is

\[
\partial_t(\chi_B n)
+\nabla_4\!\cdot\!\left(\chi_B n\mathbf v_{\rm med}+\mathbf J_\chi\right)
=n\Gamma_B,
\qquad
\Gamma_B=\Gamma_{\rm return}-\Gamma_{\rm drain}.
\]

Here:

- \(\Gamma_{\rm drain}\) converts ordered brane material into de-structured bulk material at throats;
- \(\Gamma_{\rm return}\) converts bulk material back into the ordered brane state;
- \(\mathbf J_\chi\) transports order independently of total material transport.

In an inertial branch, \(\chi_B\) instead obeys the selected Euler–Lagrange equation. The balance above must not be imposed as a second independent equation for the same variable. Any \(\mathbf J_\chi\), \(\Gamma_B\), drain, or return balance must be derived from the inertial dynamics or assigned to a distinct coarse-grained ordered-fraction variable, which may be denoted \(f_B\). A mixed branch must state explicitly whether its order field and material fraction are the same effective variable or separate coupled variables and must provide exactly one nonredundant closed evolution system.

On the fraction/relaxational branch, the quantity

\[
Q_\chi
=\int_{\Omega_{\rm throat}}
n\Gamma_{\rm drain}\,d^4X
\ge0
\]

denotes the positive gross rate of ordered-to-bulk conversion associated with a throat. It has constituent-number-per-time units, does not include distributed return, and is not itself a signed net conversion rate. Net conversion signs remain encoded by \(\Gamma_{\rm drain}\), \(\Gamma_{\rm return}\), and \(\Gamma_B\). On an inertial branch, the corresponding nonnegative gross diagnostic is **not yet derived** and this rate formula is not assumed.

This is phase conversion, not suction into empty space and not transfer between two substances. It is not automatically identical to signed local net constituent flux \(Q_n^{\rm net}\), selected positive throughput \(Q_n\), or signed coordinate flux \(\widetilde Q_n^{(w)}\). A complete theory must supply the material, momentum, and energy partners of the selected branch-specific order evolution and must make them consistent with the explicitly chosen conservative, relaxational, or mixed inertial-dissipative \(\chi_B\) branch.

In the closed-loop interpretation, throat drainage populates de-structured bulk degrees of freedom and distributed return repopulates the ordered brane state. Bulk backpressure or a chemical-potential imbalance is a candidate driver of return, not an assumed external pump. The global solution must track where drained material, momentum, energy, and, in the dissipative branch, entropy reside throughout this cycle.

The same projected law contributes to, but does not alone define, \(\mathbf g_{\rm eff}\). It must first derive the provisional source coefficient

\[
\widetilde Q_g^{(s)}(r)
=\widetilde{\mathcal F}_r\!\left[
\widetilde Q_n^{(w)},Q_n^{\rm net},Q_n,Q_\chi,
\Gamma_{\rm return},\mathbf J_\chi,
\text{bulk and brane response}
\right].
\]

No particular \(\widetilde{\mathcal F}_r\) or source-side map \(\mathcal G_{\rm src}\) is assumed. Reflection closure must reduce the raw coefficients to the orientation-even but radially signed \(\overline Q_g(r)\), after which the local attractive sign is tested separately and \(\mathcal C_{\rm ref}\) is fixed without reference to later probe species. In a dissipative or mixed effective description, any entropy, heat, or reservoir variables required to close these maps remain unresolved degrees of freedom of the same globally conserved medium, not an external source.

**Status:** projected leakage exists structurally, but the relation among \(\widetilde Q_n^{(w)}\), \(Q_n^{\rm net}\), selected \(Q_n\), and \(Q_\chi\), the source/controller functions, the raw and reflection-even signed gravity maps, the source-side reduction and fixed operational calibration, and the complete material, momentum, energy, and possible entropy budgets remain open.

## 14. Cosmological expansion hypothesis

In the fraction/relaxational branch—or in an inertial/mixed branch only after an equivalent coarse-grained balance is derived—the cosmological hypothesis is that the conversion balance contains a slow DC component for which distributed return exceeds localized throat drainage. The causal proposal is

\[
\text{distributed return excess}
\longrightarrow
\text{added ordered medium}
\longrightarrow
\text{change in intrinsic brane three-volume}
\longrightarrow
\text{possible expansion}.
\]

For a three-dimensional brane, the relevant intrinsic measure is three-volume, or hyperarea when discussing its measure as an embedded hypersurface—not ordinary two-dimensional area.

Because return acts across the brane while drainage tracks matter defects, the balance could in principle yield an early drain-dominated era followed by a late return-dominated era. This remains a cosmological expansion hypothesis, not a derived dark-energy sector.

This proposal is redistribution within a closed medium, not injection from an external energy source. If return changes the brane's intrinsic three-volume over cosmological time, the bulk reservoir's density, pressure, internal energy, and entropy must evolve consistently. A sustained expansion cannot appeal to an unchanged inexhaustible bulk battery; its duration and asymptotic behavior must follow from the global reservoir ledger.

The model must still establish:

- the sign and magnitude of the net conversion rate;
- global material, momentum, and energy conservation;
- whether returned material increases intrinsic three-volume rather than density or thickness;
- the expansion law \(\dot a\);
- whether \(\ddot a>0\);
- the corresponding effective stress-energy or cosmological equation.

Expansion does not automatically imply acceleration, and acceleration does not automatically reproduce observed dark-energy phenomenology.

**Status:** distributed return and expansion form a coherent hypothesis; neither the net sign, the three-volume response, the scale-factor dynamics, nor dark-energy-like acceleration has been derived.

## 15. Multiple characteristic speeds

The ontology contains several physically distinct propagation branches:

- \(c_s\): the bulk density or sound speed governing the leading local compressional response, not automatically every pole of the full gravity response;
- \(c_\gamma\): transverse brane-shear speed, identified with the speed of light;
- \(\omega_{E,a}(k)\): the as-yet undiagonalized electric/order branches of the coupled \(h_+\)-\(h_-\)-\(\chi_B\) sector.

Only after a healthy electric eigenbranch has been shown to possess leading \(k^2\) static stiffness and, with ordinary local inertia, a dispersion

\[
\omega_E(k)\simeq c_Ek
\]

does \(c_E\) become a well-defined electric-sector speed. A flexural branch with \(k^4\) static stiffness can instead have \(\omega\propto k^2\), but it does not generate the same Coulomb Green function. A nonstandard kinetic structure can alter this connection only if it is derived. Massive combinations may be screened, while parity alone does not guarantee which sector is massive.

The complete time-dependent gravity response may additionally contain order-conversion, interface, return, reservoir, or dissipative-relaxation poles. Its propagation cannot be assigned the single speed \(c_s\) before the coupled response operator is solved.

Equality of these speeds is not implied by dimensional analysis or by their coexistence in one medium. In particular, \(c_\gamma=c_s\) is presently a calibration target rather than a derived consequence. A future reduction of the brane moduli and throat dynamics may either derive a cone lock or expose a genuine incompatibility.

Lorentz symmetry is therefore not fundamental in the ontology. It can only be an emergent property of the effective long-wavelength sectors if their characteristic cones and couplings align closely enough. Velocity-dependent constituent transport or order conversion, low-threshold mode emission, or vacuum drag would be additional preferred-medium effects that must be calculated and constrained.

## 16. Frozen-parameter rule

Because the project is built backward from several required phenomena, consistency depends on a frozen ledger:

> Once a parameter, constitutive relation, sign, or normalization has been fixed by an earlier sector, it is frozen before proceeding. It may not be independently retuned to repair a later sector. A later modification is acceptable only if it follows from omitted parent dynamics and its consequences are propagated consistently through every completed sector.

The ledger should record:

- the quantity or constitutive choice;
- where it was fixed;
- whether it was postulated or derived;
- which later sectors depend on it;
- whether a later revision has been propagated through all dependencies.

Without this rule, separate effective knobs could make each sector appear successful without establishing a single medium.

### 16.1 Derivation-provenance table

This table records claim maturity rather than merely whether a script ran. “Derived under stated assumptions” never means that the upstream postulates or nonlinear closure have been established. A path identifies a repository artifact carrying the stated evidence; a precursor is labeled as such and does not close the current formulation. `not yet derived` means that no repository derivation closes the stated claim. The **last provenance audit** column records only a path-and-status inspection; it is not mathematical validation or a derivation rerun.

| Claim or quantity | Status | Assumptions | Derivation or script | Frozen inputs | Downstream dependencies | Last provenance audit | Failure condition |
|---|---|---|---|---|---|---|---|
| Parent arena and metric | **Postulated** | Preferred Newtonian time and flat Euclidean four-space | No derivation is needed to make the arena a postulate | \(t,\delta_{AB}\) | Every operator, Green function, wave speed, reflection, and core distance | 2026-08-18; definition checked | A different ambient metric is used without reopening and propagating all dependencies |
| \(P=Kn^5\) | **Postulated** | Stiff polytropic substrate | Historical provenance only: `research/pde_ledger_v2/notes/stage005_pathA20_source_map.md` explicitly says not to consume it as current closure; no derivation is needed to make the EOS a postulate | EOS exponent and \(K\) | \(c_s\), critical-flow tests, throat and gravity reductions | 2026-08-18; path exists, historical warning checked | EOS gives instability or cannot participate in a consistent parent action |
| \(\chi_B\) order dynamics | **Postulated branch choice / open closure** | One-medium order sector; fraction balance used only in its branch | Conservative/relaxational precursor: `research/pde_ledger_v2/paper/stages/stage_006.tex`; a fully closed conservative, relaxational, or mixed branch is **not yet derived** | Order potential, transport, reservoir variables, and whether a separate fraction \(f_B\) exists | Conversion, support, gravity, drag, cosmology | 2026-08-18; precursor scope checked | Chosen branch overdetermines \(\chi_B\) or cannot close material, momentum, energy, and entropy/reservoir ledgers |
| Finite-thickness slab existence and selection | **Open hypothesis** | Two-state order sector and suitable boundary data | Interface precursor: `research/pde_ledger_v2/paper/stages/stage_006.tex`; failed alternative wall branch: `software/stage1_solver/reports/pathA_24_T1_wall.md`; width selection **not yet derived** | Order potential and parent constitutive data | Every brane, throat, electric, and light sector | 2026-08-18; both paths and caveats checked | No stable or metastable selected slab exists |
| Shear modulus and two transverse modes | **Derived under stated assumptions** for mode count; modulus provenance open | Supplied homogeneous quadratic brane action and isotropic inertia | Current requirement/action record: `research/pde_ledger_v3/steps/S9_light_requires_shear.md`; current conditional count: `research/pde_ledger_v3/steps/S10_two_transverse_photons.md`; v2 precursor: `research/pde_ledger_v2/paper/stages/stage_003.tex` | \(\mu_R,\rho_{\rm br}\) and supplied action | Light speed, support mode, magnetism, radiation | 2026-08-18; paths and conditional scope checked | No healthy two-polarization shear sector or no acceptable modulus origin |
| Longitudinal mode | **Derived under stated assumptions** | Supplied homogeneous compressional response | Current record: `research/pde_ledger_v3/steps/S11_stray_longitudinal.md`; v2 precursor: `research/pde_ledger_v2/paper/stages/stage_003.tex` | Homogeneous elastic coefficients | Leakage, drag, extra radiation, Maxwell departure | 2026-08-18; paths and conditional scope checked | Instability or unacceptable mixing/radiation |
| Two-interface parity variables | **Operational definition / target reduction** | Reflection-symmetric graph-regime slab | **not yet derived** | \(H_{\rm br}\), solved interface profiles | Electric spectrum, finite-thickness kernels, contact geometry | 2026-08-18; no closing path identified | No consistent graph-regime reduction or parity assignment |
| Constituent conservation, flux convention, and reflection throughput | **Postulated conservation law; identity plus target relation** | Selected \(w\)-normal reference cut, or derived geometric replacement for a curved cut, and reflected environmental data | Historical flux-law audit: `software/stage1_solver/reports/pathA_20_velocity_constants.md`; accepted stationary throat flux and reflected-throughput equality **not yet derived** | Current law, EOS, mouth geometry, boundary/reservoir data | Critical flow, conversion, gravity reduction | 2026-08-18; path exists and underdetermined verdict checked | Conservation fails, cut geometry invalidates the identity without replacement, no positive branch exists, or reflected throughput differs |
| Ideal critical-flow estimate | **Derived only in an idealized limit** | Steady inviscid isentropic nozzle flow, fixed geometry/upstream state, no quantum-pressure correction | Precursor branch audit: `software/stage1_solver/reports/pathA_20_velocity_constants.md`; accepted throat flux law **not yet derived** | EOS and chosen upstream data | Necessary radius/throughput bounds | 2026-08-18; path and caveat checked | Full throat has no corresponding critical branch or violates the assumed limit |
| Support seed and nonlinear periodic throat | **Target relation** | Provisional variable-coefficient geometry followed by self-consistent free-boundary/periodic continuation | **not yet derived** | Solved throat, shear law, flow/conversion, and open-system boundary data | Stability, particle lifetime, radiation | 2026-08-18; no closing path identified | No bound/outgoing-resonance seed, nonlinear continuation, or acceptable Floquet/response stability with nonzero throughput |
| Inverse-square gravity under sink/return assumptions | **Derived under stated assumptions** | Effective localized sink/return completion | Conditional precursor: `software/stage1_solver/reports/pathA_29_brane_bulk_return.md`; full gravity reduction **not yet derived** | Effective return ansatz and background data | Gravity range and matching program | 2026-08-18; path and conditional scope checked | Full reduction changes range or sign unacceptably |
| Gravity source, reflection reduction, and calibration | **Target relation / operational calibration** | Full projected continuity, stress, conversion, return, reservoir, and reflected environment | \(\widetilde Q_g^{(s)}\to\overline Q_g\), \(\mathcal G_{\rm src}\), and \(\mathcal C_{\rm ref}\): **not yet derived** | Parent dynamics and solved local/global state | Active mass, passive response, cosmology | 2026-08-18; no closing path identified | No probe-independent source field, signed reflection closure, or fixed calibration exists |
| Active, passive, and inertial mass | **Operational definitions**; extraction open | Defined low-energy weak-field homogeneous regime, valid local matching, positive ordinary response, and passive reservoir | **not yet derived** | Fixed gravity normalization and solved moving throat | Source universality and free-fall analog | 2026-08-18; no closing path identified | Passive/inertial response is non-positive, retarded response is active without a prepared reservoir, or low-energy ratios are not species universal |
| Complete odd/even source, Hessian, and retarded blocks | **Target relation** | Equilibrium \(\mathcal H_{\rm odd/even}\) and dynamic \(\mathcal O_{R,{\rm odd/even}}\) on the complete constrained system | One-height precursor only: `research/pde_ledger_v2/paper/stages/stage_030.tex`; current parity blocks, source profiles, and support coupling **not yet derived** | Stable slab and all coupled fields | Coulomb range, thickness forces, damping, radiation | 2026-08-18; path and precursor scope checked | No healthy equilibrium odd \(k^2\) branch, inconsistent retarded spectrum, or unacceptable additional long range |
| \(1/R\) electric far field | **Conditional precursor**; current claim open | Postulated scalar/interface branch and selected boundary data | One-height precursor: `research/pde_ledger_v2/paper/stages/stage_031.tex`; current two-interface carrier **not yet derived** | Candidate operator and source branch | Electric force and factorization | 2026-08-18; path and precursor scope checked | Full carrier lacks \(k^2\) stiffness or gives wrong range |
| Electric force sign | **Target relation** | Correct equilibrium ensemble and self-energy subtraction, derived conservative-periodic averaged ensemble or stress force, or complete driven stress/momentum/reservoir force ledger | Boundary-ensemble precursor/audit: `research/pde_ledger_v2/paper/stages/stage_032.tex`; sign **not yet derived** | Solved core boundary, support action/frequency/phase, and branch data | Attraction/repulsion and bound states | 2026-08-18; path and open status checked | Desired sign requires an unsupported static/periodic ensemble, unjustified potential, or fitted mechanism |
| Leading odd-channel factorization \(g_ag_b\) | **Target relation** | Local Coulomb carrier plus universal solved throat couplings | **not yet derived** | Multiple solved throat species and frozen normalization | Charge universality, magnetic/radiative normalization | 2026-08-18; no closing path identified | Odd monopole coefficient is not acceptably rank one/positive or its value/derivative corrections are not asymptotically subleading |
| Conserved oriented current | **Target relation** | Co-orientable brane and conserving throat topology | **not yet derived** | Complete \(\mathcal C_{\rm th}\)/topological map | Additive charge, magnetism, radiation | 2026-08-18; no closing path identified | No conserved count/current or global orientation exists |
| Magnetic tensor structure | **Derived under stated assumptions** at structural level | Rigid moving profile and postulated transverse response | `software/em_charge_attribute/magnetism_moving_throat_result.md`; normalization and electric tie remain open | Supplied transverse operator and provisional source | Dynamic electromagnetic closure | 2026-08-18; path and conditional/R1 scope checked | Full moving throat changes tensor/falloff or needs independent normalization |
| Transverse radiation from acceleration | **Target relation** | Solved accelerating throat and far-zone outgoing conditions | **not yet derived** | Conserved current and shared coupling | Electromagnetic closure and energy ledger | 2026-08-18; no closing path identified | No tied transverse radiation or unacceptable extra channels |
| Projected coincidence and core-contact criterion | **Operational geometric definition** | Graph regime for \(d_w^{\rm core}\); orientation-resolved diagnostics or connected-component continuation with compact or suitably closed localized level sets otherwise | **not yet derived** | Solved diagnostic, threshold range, component labels, level sets, and interface geometry | Collision, annihilation, neutral binding | 2026-08-18; no closing path identified | No graph-to-core handoff, physical label-continuation rule, threshold-robust contact test, or ambient four-dimensional set distance |
| Same-species annihilation and different-species neutral binding | **Open hypothesis / target dynamics** | Complete \(\mathcal C_{\rm th}\) partners and nonlinear time evolution | **not yet derived** | Full internal map, core topology, and conservation laws | Particle phenomenology | 2026-08-18; no closing path identified | False universal annihilation or no conserving antiparticle channel |
| Local/global drain-return fixed point and cosmological expansion | **Open hypothesis** | Closed one-medium reservoir and cosmological reduction | Conditional return precursor: `software/stage1_solver/reports/pathA_29_brane_bulk_return.md`; fixed point and expansion **not yet derived** | Global material, momentum, energy, and entropy data | Return scale, active mass, expansion history | 2026-08-18; path and conditional scope checked | No local–global fixed point or no acceptable expansion law |

Any change to a frozen input requires updating its row, rerunning every real downstream derivation listed here, and marking dependent claims as reopened until that propagation is complete.

The [mechanical synchronization check](../software/check_ontology_shared_definitions.py) guards the repeated definitions of \(h_m,h_t\), the two-minus interface reflection map, the flux identity, \(\mathcal C_{\rm th}\), and the graph-regime core gap across this ontology and its focused companion. A passing result establishes documentation consistency only; it is not mathematical or physical validation.

## 17. Explicit failure conditions

The present candidate architecture fails if any required result proves impossible, including:

1. no finite-thickness, shear-supporting brane exists whose energy is bounded below and whose interfaces and perturbation spectrum are stable or sufficiently metastable;
2. no normalizable bound seed or acceptably long-lived outgoing resonance exists, or no seed continues to a self-consistent nonlinear periodic supported throat with nonzero throughput and acceptable Floquet/response stability;
3. no stationary drain branch exists with \(Q_n^{\rm net}>0\), or independently solved reflection partners fail to have equal positive throughput under reflection-symmetric environmental data;
4. exact constituent conservation cannot be maintained; the aperture cut does not support the claimed flux identity and no correct geometric replacement is supplied; or no consistent relation among \(\widetilde Q_n^{(w)}\), \(Q_n^{\rm net}\), selected \(Q_n\), gross local conversion drain \(Q_\chi\), distributed return, raw gravity source, and signed \(\overline Q_g(r)\) can be derived;
5. no source-side gravity field can be derived independently of the later probe species, or the fixed operational calibration absorbs species-dependent response and makes passive universality tautological;
6. reflected solutions fail to transform covariantly with reflected reservoir data, or the ordinary vacuum branch is not sufficiently reflection symmetric to avoid unacceptable charge-correlated gravitational splitting;
7. the reflected solutions fail to produce one orientation-even gravity-source coefficient, or the ordinary local particle branch fails to satisfy \(\overline Q_g(r_{\rm local})>0\);
8. screening, return overshoot, or a sign reversal at larger radius is hidden by an absolute-value convention, is unacceptably species dependent, or violates observational constraints;
9. no matching plateau or equivalent scale-independent near-to-far extraction exists for local active mass, or the map to the signed scale-dependent effective source is neither universal nor environmentally calculable;
10. a conservative/equilibrium coupled system has no healthy odd Hessian eigenbranch with nonzero throat coupling and the \(k^2\) static stiffness required for a \(1/R\) potential, a conservative periodic branch cannot derive the corresponding acceptable averaged force, or a driven branch cannot derive it from stress and momentum-flux response;
11. no complete odd/even decomposition of the solved throat source exists, additional odd/even long-range modes produce unacceptable forces or radiation, or the full parity-block spectrum contains ghosts or instabilities;
12. the coupled \(V_{\rm mat}(\chi_B,\sigma)\) boundary problem cannot provide bounded energy, a stable brane, a physical oriented core branch, and freedom from unacceptable modes; an equilibrium calculation cannot obtain the required interaction sign from the correct energy ensemble; a conservative periodic calculation assigns a static potential without deriving the correct cycle-averaged ensemble or stress force; or a driven dissipative calculation assigns a potential without first deriving a conservative force from stress, material/order momentum flux, conversion, and reservoir response;
13. the odd-monopole coefficient \(C_{E,ab}\) is used to hide even or nonlinear forces, its value or derivative corrections are not asymptotically subleading, or common slab conditions plus derived core data cannot produce sufficiently universal leading electric coupling among ordinary unit-charge throats;
14. the brane geometry does not support a global co-orientation or equivalent topological construction needed for additive signed throat number, or no conserved oriented current/Gauss-like projected law exists;
15. the ordinary stable branch has a non-positive passive or inertial response direction after conventions are fixed, or a retarded response extracts net work from an unexcited internal reservoir;
16. either low-energy \(M_{\rm active}^{\rm local}/M_{\rm inertial}\) or low-energy \(M_{\rm passive}/M_{\rm inertial}\) cannot be sufficiently species universal in its defined operational regime;
17. in a conservative branch, unavoidable sub-threshold emission or zero-critical-velocity coupling produces unacceptable drag;
18. in a relaxational or mixed branch, low-frequency friction, relaxation, or conversion lag produces unacceptable drag during uniform motion;
19. the chosen conservative, relaxational, or mixed order dynamics overdetermines \(\chi_B\) or cannot close its nonredundant material, momentum, energy, inertia/drag, and required entropy/internal-reservoir ledgers;
20. transverse light leaks or converts excessively into longitudinal, bulk, odd-electric, even-thickness, or conversion modes;
21. an accelerating oriented throat cannot source the two transverse light modes with normalization tied to the same leading electric coupling as static electricity and magnetism, or it necessarily produces unacceptable unsuppressed additional radiation;
22. no self-consistent fixed point exists between the local supported throat and the global drain-return/reservoir boundary conditions;
23. no consistent handoff exists from projected positions \(\mathbf x_a\) and two-interface graphs to full points \(\mathbf X=(\mathbf x,w)\), compact or suitably closed localized diagnostic core sets measured with \(\mathbf Y_\pm\), and nonlinear parent-field topology; neither orientation-resolved diagnostics nor connected-component continuation supplies physical core labels through contact; separate labels are imposed after merger; the inferred contact is not robust under reasonable threshold variation; or projected coincidence causes an unphysical singularity, arbitrary branch selection, unacceptable short-range force, or loss of particle identity;
24. no complete \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\) produces an equilibrium equal-energy partner, a conservative periodic partner with consistently mapped action/frequency/phase data and cycle-averaged invariants, or a driven partner with equality of the applicable stationary throughput, dissipation, entropy production, lifetime, stress/momentum fluxes, and retarded response kernels;
25. a same-species throat–antithroat pair has no material-, momentum-, and energy-conserving annihilation or neutralization channel, or the theory incorrectly forces distinct-species opposite charges to annihilate solely because their odd electric disturbance cancels;
26. any desired gravity, electric, contact, radiation, support, or cosmological result can be obtained only by introducing a separately fitted coefficient or retuning a frozen upstream input without propagating the change.

These are not admissions that the model has already failed. They define the computational tests that prevent the ontology from being protected by post-hoc adjustment.

## 18. What the model does not contain

The current ontology contains:

- no fundamental point particles;
- no identification of observed particle mass with trapped-wave energy;
- no identification of observed particle mass with a sum of constituent rest masses;
- no requirement that active, passive, and inertial mass have the same microscopic mechanism;
- no electric interaction magnitude assigned independently by hand—charge sign is intrinsic to orientation, while a leading one-throat coupling may emerge only from the shared brane boundary-value problem;
- no second material substance;
- no empty bulk vacuum that actively sucks material away;
- no ordinary linearly sourced scalar explanation of electric repulsion;
- no assignment of a static interaction energy to a conservative periodic pair without a derived cycle-averaged ensemble or averaged stress-force construction;
- no assignment of an equilibrium interaction potential to a maintained dissipative pair unless its quasistatic force is first shown to be conservative;
- no use of the odd Coulomb coefficient \(C_{E,ab}\) to absorb reflection-even thickness or nonlinear core forces;
- no identification of a proton-like object with an electron-like throat whose orientation has merely been reversed—that operation defines the positron-like branch;
- no claim that electric-orientation reversal alone completes the antiparticle map;
- no identification of raw \(\mathbf v_{\rm med}\) with the calibrated \(\mathbf g_{\rm eff}\);
- no treatment of graph-regime interface coordinates as complete nonlinear core topology;
- no assumption that a gapless or reflection-odd mode is automatically Coulomb-like;
- no assumption that a reflection-even mode is automatically screened;
- no fundamental electromagnetic gauge field inserted solely to reproduce Maxwell;
- no identification of the longitudinal mode with electric charge;
- no automatic equality of gravity, light, and electric propagation speeds;
- no claim that the model is general relativity or exact electromagnetism;
- no claim that the stable brane or nonlinear particle interior has already been solved.

The model instead seeks a single material system whose effective sectors approximate selected GR-like and Maxwell-like behaviors while producing explicit, falsifiable departures where exact agreement is impossible.

## 19. Closure sequence

The ontology must be closed in dependency order:

\[
\boxed{
\begin{aligned}
\text{parent dynamical framework}
&\longrightarrow \text{stable finite brane}\\
&\longrightarrow \text{background-brane spectrum}\\
&\longrightarrow \text{stationary supported throat}\\
&\longrightarrow \text{throat stability and internal spectrum}\\
&\longrightarrow \text{moving-throat family}\\
&\longrightarrow \text{two-throat interactions}\\
&\longrightarrow \text{global conversion balance}.
\end{aligned}
}
\]

This sequence is the practical dependency order, not a one-way physical closure. The local and global problems form the feedback loop

\[
\boxed{
\text{local supported throat under asymptotic reservoir data}
\;\longleftrightarrow\;
\text{global drain-return and reservoir solution}
}.
\]

Introduce provisional environmental data

\[
\mathcal B_\infty
=\left\{
n_\infty,
P_\infty,
\mu_\infty,
\chi_\infty,
T_\infty,
\text{return profile},
\text{reservoir state}
\right\},
\]

retaining only the variables appropriate to the selected conservative, relaxational, or mixed branch. Closure requires the following iteration:

1. solve the local supported throat for provisional \(\mathcal B_\infty\);
2. extract \(Q_n^{\rm net}\), selected positive \(Q_n\), \(\widetilde Q_n^{(w)}\), conversion, stress, momentum, energy, entropy where applicable, and emitted disturbances;
3. solve the global bulk, return, and reservoir problem;
4. update backpressure and \(\mathcal B_\infty\), reflect environmental data when comparing orientations, and update the raw source, signed \(\overline Q_g(r)\), and source-side field;
5. iterate until a self-consistent fixed point is obtained.

No local result is fully closed if it relies on externally frozen reservoir data that the resulting global solution would change.

### 19.1 Background brane

The throat-free slab problem must determine:

- the selected local thickness \(H_{\rm br}\);
- whether a symmetric solution has a relational mid-surface with \(\chi_0(w)\) maximized there and how the two bulk-facing regions participate in thickness selection;
- the two-interface variables \(h_+,h_-\), their leading parity projections \(h_m,h_t\), the equilibrium Hessian blocks \(\mathcal H_{\rm odd/even}\), and the dynamic retarded blocks \(\mathcal O_{R,{\rm odd/even}}\), rather than one assumed scalar branch per parity;
- the graph-regime validity domain and handoff from interface collective coordinates to full parent-field topology;
- the origin and value of the shear modulus;
- energy bounds and interface stability;
- two transverse shear modes;
- the in-plane longitudinal mode;
- normal displacement or flexural modes;
- whether a healthy odd eigenbranch has \(k^2\) static stiffness and whether any flexural \(k^4\) branch is distinct;
- whether a thickness-like even eigenbranch is actually screened after all mixing is included;
- whether additional odd or even long-range branches have nonzero throat coupling and acceptable forces and radiation;
- order-variable and bulk-sound modes;
- brane–bulk leakage and mixing;
- all relevant dispersions \(\omega_a(k)\);
- the background coefficients and reservoir data needed by the source-side gravity reduction \(\mathcal G_{\rm src}\);
- the asymptotic localization properties of the transverse support operator needed for a bound-state or resonance test;
- the constitutive data needed to evaluate \(V_{\rm mat}(\chi_B,\sigma)\).

### 19.2 Stationary supported throat

The stationary or time-periodic supported-throat solution must determine:

- a linear bound-state or outgoing-resonance seed on a provisional complete variable-coefficient geometry;
- for a bound seed, its weighted norm, finite energy, and asymptotic decay; for a resonant seed, its outgoing conditions, complex frequency or equivalent leakage rate, and lifetime;
- continuation of that seed to a self-consistent nonlinear free-boundary or periodic supported throat whose time-averaged support stress, geometry, flow, conversion, and reservoir loading agree;
- whether the diagnostic localization region \(\Omega_{\rm support}\) is a rim, collar, sheath, cavity wall, or another connected solved geometry;
- how its averaged pressure or stress holds the aperture open;
- the equilibrium throat radius and shape;
- the support of each oriented core across the slab and whether it is anchored to one interface, spans the slab, or requires a multi-sheet geometry;
- signed local net outward flux \(Q_n^{\rm net}\), selected positive drain-branch throughput \(Q_n\), and signed coordinate-normal flux \(\widetilde Q_n^{(w)}=sQ_n^{\rm net}\) on the selected \(w\)-normal reference cut, or the derived geometric replacement for a curved cut;
- nonnegative gross local order-conversion drain \(Q_\chi\), including its branch-specific definition and number-per-time units;
- gross shorthand \(Q_{\rm ap}\), if useful, without conflating those processes;
- whether displacement relative to the brane mid-surface produces the proposed conversion-and-entrainment feedback;
- whether a sonic or other critical-flow condition bounds \(Q_n\) for the solved mouth geometry and upstream state, and whether that bound participates in the equilibrium-radius relation;
- the provisional raw gravity coefficient \(\widetilde Q_g^{(s)}(r)\), its covariance under \((\mathcal T,\mathcal B_\infty)\mapsto(\mathcal C_{\rm th}\mathcal T,\mathcal R_w\mathcal B_\infty)\), and the orientation-even signed \(\overline Q_g(r)\);
- the ordinary local sign test \(\overline Q_g(r_{\rm local})>0\), with local active mass extracted from a derived plateau \(R_{\rm throat}\ll r\ll L_{\rm return}\) or an equivalent scale-independent near-to-far prescription;
- the probe-independent source-side field \(\boldsymbol{\mathfrak g}_{\rm src}\) and one fixed calibration \(\mathcal C_{\rm ref}\) to \(\mathbf g_{\rm eff}\);
- the charge-orientation boundary branch;
- the coupled two-interface/order core structure;
- the complete odd/even source vector rather than only bare \(h_m\)- and \(h_t\)-like entries;
- the complete candidate map \(\mathcal C_{\rm th}=\mathcal I_{\rm internal}\circ\mathcal R_w\) of every orientation-even and orientation-odd internal datum, including whether \(\mathcal I_{\rm internal}\) is trivial;
- consistency with provisional asymptotic reservoir data \(\mathcal B_\infty\);
- whether the core boundary amplitude and leading electric coupling are universal, quantized, topology controlled, or species dependent.

Its geometry, averaged stress, and mean constituent-transport and order-conversion rates may be stationary while its internal standing mode remains periodic.

### 19.3 Throat stability and internal spectrum

Perturbations around the supported throat must determine:

- translation zero modes;
- orientation-changing channels;
- same-species throat–antithroat contact, neutralization, reconnection, and annihilation channels;
- radius and shape oscillations;
- trapped-standing-wave perturbations;
- unstable or negative-energy modes;
- radiation channels into every background branch;
- normalization or outgoing-resonance data and energy localization of the spectral seed while throughput remains nonzero;
- the Floquet or linear-response spectrum about the self-consistent nonlinear periodic throat, rather than only stability of the seed eigenproblem.

A viable throat must have no unphysical negative-energy modes, no growing linear or Floquet modes, no spontaneous collapse, widening, or orientation flip, and sufficiently long-lived metastability if exact stability is unavailable.

### 19.4 Moving-throat family

The moving family must determine:

- standing-mode reconfiguration;
- near-field rebuilding;
- for the reversible branch, total canonical momentum \(P_i(\mathbf v)\) and \(\partial P_i/\partial v_j\);
- for a relaxational or mixed branch, the frequency-dependent inertia, conversion, and low-frequency damping kernels with complete internal-reservoir accounting;
- local active, passive, and inertial mass relations;
- passive response kept separate from the fixed source-field calibration, so later test-species dependence remains measurable;
- positive-definite symmetric passive and inertial response on the ordinary stable branch;
- nonnegative cycle-averaged absorbed power for perturbations of an unexcited dissipative reservoir;
- \(\widetilde Q_n^{(w)}(v)\), \(Q_n^{\rm net}(v)\), selected \(Q_n(v)\), the nonnegative gross drain magnitude \(Q_\chi(v)\), the raw gravity source, and signed \(\overline Q_g(r;v)\) after reflection closure;
- electromagnetic magnetism;
- one-current normalization across \(O(v^0)\) electricity, \(O(v)\) magnetism, and \(O(a)\) transverse radiation;
- far-zone powers \(P_{u_T},P_{\rm odd},P_{\rm even},P_{u_L},P_{\rm bulk},P_\chi\);
- reversible-branch mode-emission thresholds;
- dissipative-branch uniform-motion friction and conversion lag;
- whether acceleration or uniform motion radiates into any branch.

### 19.5 Two-throat interactions

The two-body problem must determine:

- gravitational force normalization;
- whether the pair is equilibrium static, conservative phase-locked periodic, or driven dissipative, and the corresponding static ensemble, fixed support-mode action/frequency/relative-phase data, or stress/momentum/reservoir force prescription;
- reflected or asymmetric environmental data and any resulting orientation-dependent splitting;
- the effect of return-halo overlap on the signed \(\overline Q_g(r)\) rather than on an absolute source magnitude;
- projected throat positions \(\mathbf x_a\) and the full four-dimensional core-set distance using dummy points \(\mathbf Y_\pm\);
- the general odd/even source blocks and all unmatched species-dependent components;
- the odd-monopole decomposition \(C_{E,ab}(R)=g_ag_b+\delta C_{ab}^{\rm odd}(R)\), with \(g_a,g_b\ge0\), all leading electric sign carried by \(s_a,s_b\), and both \(\delta C_{ab}^{\rm odd}\) and \(R\partial_R\delta C_{ab}^{\rm odd}\) asymptotically subleading;
- separate reflection-even, parity-mixed, and nonlinear-core forces rather than hiding them inside \(C_{E,ab}\);
- whether leading factorization occurs;
- whether core data are universal, quantized, or topology controlled;
- whether thickness, environmental, finite-size, and polarization corrections are acceptably small;
- whether many-body effects preserve approximate additivity;
- the electric force sign from the correct equilibrium ensemble, the derived conservative-periodic averaged functional or stress force, or, for a driven pair, the complete stress, material/order momentum-flux, conversion, and reservoir response;
- the finite-thickness cross-interface Green matrix and its \(R\sim H_{\rm br}\) corrections;
- whether oppositely oriented throats can coincide in projected three-space while remaining separated in \(w\);
- the graph-regime gap \(d_w^{\rm core}=H_{\rm br}+2h_t^{\rm tot}-\ell_+^{\rm in}-\ell_-^{\rm in}\), a solved localized core diagnostic, thresholded sets, orientation-resolved intersection or component-closure/topology transition, and the threshold-robust distance \(d_{\rm core}(\epsilon)\) while two cores remain identifiable after graph breakdown;
- orientation-resolved diagnostics or connected-component continuation that keeps the two core labels physical up to contact and stops assigning separate labels after merger;
- whether the exact-coincidence limit gives exclusion, saturation, annihilation, reconnection, or a stable neutral composite, with same-species antiparticles distinguished from different or composite species;
- overlap of return halos and the modification of each throat's \(Q_n^{\rm net}\), selected throughput, conversion, raw gravity source, signed effective source, and active mass by the other's environment;
- conditions under which isolated one-throat profiles approximately superpose;
- momentum carried by the medium and return flow;
- environmental corrections to electric coupling and active mass;
- whether a regime \(R_{\rm throat}\ll R\ll L_{\rm return}\) actually exists rather than being assumed;
- magnetic velocity-dependent interaction;
- nonlinear short-range behavior.

### 19.6 Global conversion

The global problem must determine:

- localized drainage and distributed return;
- a nonredundant branch-specific order evolution: the \(\chi_Bn\) balance for a fraction/relaxational branch, an Euler–Lagrange equation for an inertial order field, or explicitly coupled \(\chi_B\) and \(f_B\) variables where both are needed;
- the \(\widetilde{\mathcal F}_r\) relation among \(\widetilde Q_n^{(w)}\), \(Q_n^{\rm net}\), selected \(Q_n\), \(Q_\chi\), return, and the raw gravity coefficient;
- the reflected-environment reduction to signed \(\overline Q_g(r)\), the probe-independent source-side field \(\boldsymbol{\mathfrak g}_{\rm src}\), and fixed calibration \(\mathcal C_{\rm ref}\) to \(\mathbf g_{\rm eff}\);
- the return or gravity-screening scale, if one exists;
- whether \(\overline Q_g(r)\) approaches zero, a constant, changes sign, or has another asymptotic form;
- whether the return profile is universal or environmentally controlled;
- the evolving density, pressure, internal energy, and entropy of the bulk reservoir throughout any sustained drain-return cycle;
- material conservation;
- momentum and energy conservation;
- entropy production and internal reservoir accounting if the selected branch is relaxational or mixed;
- long-term brane evolution.

The global solution must update \(\mathcal B_\infty\) and be iterated with the local throat until the same reservoir state, backpressure, return profile, transport, conversion, and emitted fluxes are recovered on both sides of the feedback loop.

## 20. Compact mental picture

One globally conserved four-dimensional medium lives in preferred Newtonian time and flat Euclidean parent space \((\mathbb R^4,\delta_{AB})\) and forms an ordered, finite-thickness, shear-supporting three-dimensional brane. The present candidate assumes the brane is globally two-sided and co-orientable, with its local mid-surface defined relationally. A finite oriented throat at projected position \(\mathbf x_a\) extends the ordered structure into a bulk-facing defect. Its complete geometry lives at full points \(\mathbf X=(\mathbf x,w)\). A linear normalizable bound mode or outgoing resonance on a provisional geometry is only a spectral seed. The particle must be a self-consistent nonlinear periodic supported throat whose averaged support stress, geometry, flow, conversion, and reservoir loading close together, followed by Floquet or response stability. A rim, collar, sheath, cavity wall, or other \(\Omega_{\rm support}\) is diagnosed from the solved energy and stress localization. Support energy is structural and is not itself observed mass.

A throat has signed local net outward flux \(Q_n^{\rm net}\), signed global-coordinate flux \(\widetilde Q_n^{(w)}=sQ_n^{\rm net}\) on its selected \(w\)-normal reference cut, and, on a selected stationary drain branch, positive throughput \(Q_n=Q_n^{\rm net}>0\). Exact total conservation defines \(\mathbf J_n\); a curved aperture needs the invariant surface flux and a derived projected-area relation. Its transport, branch-specific nonnegative gross conversion drain \(Q_\chi\), stress, distributed return, and reservoir response must first determine the probe-independent source field \(\boldsymbol{\mathfrak g}_{\rm src}\). One fixed \(\mathcal C_{\rm ref}\) calibrates that field into \(\mathbf g_{\rm eff}\), which is distinct from raw \(\mathbf v_{\rm med}\); passive response is then derived separately. Reflection covariance compares the complete candidate throat operation \(\mathcal C_{\rm th}\) with geometrically reflected environmental data, and same-environment degeneracy requires \(\mathcal R_w\mathcal B_\infty=\mathcal B_\infty\). Reflection closure yields orientation-even but radially signed \(\overline Q_g(r)\); local active mass is defined only where the ordinary branch has \(\overline Q_g(r_{\rm local})>0\), while larger-scale screening, overshoot, or sign reversal remains visible. Ordinary passive and inertial response must be positive, dissipative response must be passive for an unexcited reservoir, and mass-ratio universality is required only in the defined weak-field, low-velocity, low-frequency regime.

The two brane interfaces have leading reflection-odd mid-surface and reflection-even thickness projections, but the complete source and response decompose into parity blocks that may contain multiple eigenbranches. Electric orientation imposes a signed material boundary condition, while even thickness and nonlinear core dressing may remain after the cycle-averaged or correctly phase-mapped odd source cancels for a candidate pair. A Coulomb-like equilibrium tail requires a healthy odd Hessian eigenbranch with nonzero throat coupling and leading \(k^2\) static stiffness; parity or gaplessness alone is insufficient, and a flexural \(k^4\) branch is a different static problem. The coefficient \(C_{E,ab}\) belongs only to that odd Coulomb monopole channel; even, mixed, and core interactions remain separate, and both the odd correction and its radial derivative must be asymptotically subleading for the leading force sign. Equilibrium interactions use \(\mathcal H\) and an ensemble-appropriate potential; a conservative periodic pair requires a derived cycle-averaged action/frequency/relative-phase ensemble or averaged stress force; damping and radiation use \(\mathcal O_R\); and a driven dissipative pair gets its force from stress and momentum/reservoir ledgers unless conservativity is derived. Opposite orientations can have \(\mathbf x_+=\mathbf x_-\) while their cores remain separated. In the graph regime their gap is \(H_{\rm br}+2h_t^{\rm tot}-\ell_+^{\rm in}-\ell_-^{\rm in}\); interface coordinates hand off to full parent fields and thresholded localized core sets. Orientation-resolved diagnostics or connected-component continuation must supply physical core labels through contact and retire them after merger; \(d_{\rm core}(\epsilon)\) is an ambient Euclidean precontact measure while two cores remain identifiable. The complete \(\mathcal C_{\rm th}\), not geometric reflection alone, defines only the antiparticle candidate until the branch-appropriate partner tests are derived.

Electromagnetic closure requires one conserved oriented current and one derived coupling structure to connect the \(O(v^0)\) static electric response, the \(O(v)\) magnetic response, and \(O(a)\) emission into the two transverse light polarizations. Radiation into odd-normal, even-thickness, longitudinal, bulk, and conversion/reservoir modes is a calculable prediction. The local supported throat and the global drain-return/reservoir problem must be iterated to a common fixed point with frozen parameters and complete material, momentum, energy, and, where applicable, entropy accounting. Cosmological expansion through distributed return remains an open closed-reservoir hypothesis, not a derived dark-energy sector.

The central test is whether one parent framework and one self-consistent family of brane, throat, moving, interacting, and global reservoir solutions can produce all of these behaviors without independent retuning. Every unclosed relation remains a target or failure test in the provenance ledger rather than an assumed success.
