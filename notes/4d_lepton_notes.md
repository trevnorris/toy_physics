## 1. Session frame and target

### 1.1 Why we started with a falsification-first lepton audit

The immediate purpose of this phase is **not** to jump straight into the full moving-throat PDE. The handoff explicitly asks us to use the same research style that produced the strongest earlier results in the toy-model series: choose a destination with a hard known answer, work backward to the physical requirements that must hold, and then test whether the current ontology can host those requirements before attempting a full first-principles derivation. This is the project’s established “known-destination, falsification-first” method, and the handoff argues that particle-identity questions should be approached the same way. 

That pivot makes sense in the current project state. The 1PN program is already frozen through its declared conservative closure hierarchy, but it also explicitly leaves spin couplings, dissipative leakage/radiation reaction, strong-field completion, and the fully solved moving-throat PDE outside its solved scope. At the same time, the 2PN notes show that pushing the model harder keeps revealing new particle/throat structure rather than merely refining old coefficients. So the fastest way to learn whether the model can really support particle physics is still to ask what a particle would need **before** we bury those requirements inside a much larger PDE problem.  

### 1.2 What counts as success for this phase

Success in this phase does **not** mean a full Standard Model derivation. The handoff narrows the target to a much smaller and sharper question: can the toy-model defect coherently host the **minimal charged-lepton package**? That means identifying whether the model can support a genuine intrinsic low-energy two-state sector, a rest-frame spin-like observable that is not ordinary orbital motion, an intrinsic magnetic moment with the Pauli/Dirac scaling (\mu \sim (q/m)S), and ideally some sign of a (2\pi) versus (4\pi) identity distinction. A successful phase result could therefore be either a constructive route to that package or an early no-go showing exactly which required structure is missing. 

Failure is just as informative here as success. The handoff is explicit that the audit should reject candidates early if they reduce to ordinary circulation, rigid-body rotation, or some unrelated ad hoc patch that does not reuse the throat ontology. In other words, this phase is supposed to answer a hard structural question: does the current defect picture even have room for a lepton-like internal identity, or does it stall before that? 

### 1.3 The minimal charged-lepton package as the destination

The external anchor is the ordinary charged-lepton low-energy package. At the field-equation level, that means keeping the Dirac and Pauli structures in view:
[
(i\gamma^\mu D_\mu - m)\psi = 0,
\qquad
\mu = g,\frac{q}{2m},S,
]
with the leading Dirac target (g=2). But the handoff is clear that the real issue is not formal symbolism; it is whether the toy model can produce the **physical content** that those equations encode: an intrinsic two-state internal degree of freedom, an intrinsic magnetic dipole law, and a reason that the resulting identity is not just classical vorticity or orbital motion in disguise. 

That makes the “minimal charged-lepton package” the right first particle target for the program. It is close enough to the existing ontology to be meaningful—because the model already has charge-like circulation, internal throat geometry, and a nontrivial mouth/wake response structure—but still sharp enough to falsify the model quickly if the necessary internal structure simply cannot exist.  

### 1.4 Endpoint equations and observables to aim at

For this audit, the destination is best formulated as a small set of endpoint equations and observables. First, the low-energy object should behave like a charged spin-(\tfrac12) defect rather than a scalar monopole, so the target magnetic coupling is
[
H_{\rm Pauli}\sim \cdots + q\phi - \mu!\cdot!B,
\qquad
\mu = g,\frac{q}{2m},S.
]
Second, the model should contain a **rest-frame** internal observable (S) that survives after ordinary orbital motion is removed. Third, the handoff asks us to look, where possible, for a (2\pi/4\pi) identity signal rather than stopping at classical circulation alone. These are the equations and observables against which the next-stage toy-model structures have to be judged. 

## 2. Starting map from the frozen papers

### 2.1 What the earlier program already fixes

Before the 4D program, the earlier papers had already frozen several nontrivial viability conditions. The Newtonian sector is carried by an instantaneous Poisson structure. The orbital 1PN program fixed the perihelion-precession ledger through an effective inertia sector with (\beta=3), together with the decomposition that later became (\kappa_\rho=1), (\kappa_{\rm add}=1/2), and (\kappa_{\rm PV}=3/2) in the declared adiabatic closure. The weak-field optical sector fixed the barotropic exponent to
[
P(\rho)=K_{\rm EOS}\rho^n,\qquad n=5,
]
and the vector/wake sector fixed the parity-even wake data
[
\alpha^2=\frac34,\qquad a_H=0,\qquad K_{\rm vec}=\frac{2}{\pi^2}.
]
The EM/cavity branch also selected the familiar throat aspect ratio
[
\frac{L}{a}\approx \frac{\sqrt2,\pi}{x_{01}}\approx 1.85.
]
These are not loose preferences anymore; they are carry-forward constraints the later papers deliberately preserve.    

The full conservative 1PN assembly then repackaged these earlier results inside the action-based 4D framework and closed the two-body Einstein–Infeld–Hoffmann ledger through order (c^{-2}). In that paper the values
[
q=1,\qquad n=5,\qquad \kappa_{\rm add}=\frac12,\qquad \kappa_{\rm PV}=\frac32,\qquad \beta_{\rm 1PN}=3
]
are not ad hoc inserts anymore; they are treated as the frozen one-body and two-body constants of the conservative near-zone program. 

### 2.2 Charge and mass dictionary already present

The project also already has a working geometric/topological dictionary for “particle” properties. In the EM/throat paper, the effective charge-like quantity is tied to circulation and throat area:
[
q = \kappa_q,\rho_0,\pi a^2 \Gamma,
]
while the mass-like quantity is tied to throat geometry:
[
m_G = \kappa_m,\rho_0,\pi a^2 L.
]
This is one of the strongest reasons a lepton audit makes sense now: the model is not starting from a featureless point particle, but from a defect whose charge and inertia are already expressed in terms of circulation, size, and internal support. 

### 2.3 The 4D EM / projection / leakage structure already exists

The action-based 4D program already contains the key field-theory ingredients needed for particle-identity questions. On the matter side, it defines bulk density (\rho=|\psi|^2), a gauge-invariant velocity
[
v_i=\frac{\hbar}{m}\partial_i\theta-\frac{q}{m}A_i,
]
and the exact vorticity–gauge identity
[
\Omega_{ij}=\partial_i v_j-\partial_j v_i=-\frac{q}{m}F_{ij}.
]
On the observer side, it defines brane quantities by a projection map (\mathcal P_W) and proves that projected continuity is generically open:
[
\partial_t\rho_{\rm brane}+\nabla_3!\cdot!\mathbf J_{\rm brane}=S_{\rm leak}.
]
So the brane is already a **filter on a bulk system**, not a literal hard wall. 

The full 4D EM sector is equally important. The localized 4+1 Maxwell action yields a genuine bulk gauge theory
[
\partial_M!\left(Z(w)F^{MN}\right)+\frac{1}{\xi}\partial^N(\partial!\cdot!A)=\mu_0 J^N,
]
and under controlled zero-mode assumptions it reduces to an effective 3+1 Maxwell sector with
[
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
]
Outside that strict zero-mode regime, the theory explicitly keeps the mixed channels
[
A_w,\qquad E_w=F_{w0},\qquad C_a=F_{aw},\qquad J^w
]
as real dynamical objects. Those mixed channels later became central in the plasma paper’s beyond-MHD/open-system story, which is one reason they became natural candidates in the lepton audit.   

### 2.4 Why the 1PN program is strong but still incomplete

The 1PN assembly is strong because it closes the conservative two-body EIH structure cleanly and shows that the 4D toy ontology can reproduce the known weak-field gravitational ledger. But the paper is equally explicit about what it does **not** yet solve. Spin couplings are not part of the frozen conservative result, and neither are dissipative leakage/radiation-reaction effects, strong-field completion, or a theorem for the fully solved moving-throat PDE. So the 1PN result should be read as a very strong closure **within its declared hierarchy**, not as a completed particle theory. 

### 2.5 What 2PN added

The 2PN notes are where the particle picture visibly changes. The conservative cross sector is no longer just “more coefficients.” It becomes a small throat-response theory: a carried-forward odd dipole wake sector, a genuinely new even support bundle built from
[
P_0\oplus P_2,
]
and a separate geometry-closure channel with deficit
[
\Delta_{\rm geom}=\frac{281}{80}.
]
In the zero-frequency language this operator carries split odd residues
[
R_{1\perp}=\frac72,\qquad R_{10}=4,
]
an even support metric that is effectively the identity on the six real (P_0\oplus P_2) ports, and a scalar source vector
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right).
]
This is the crucial frozen starting point for the lepton audit: the defect is already more than a scalar monopole, but it is not yet obviously a spinor-like particle either. 

## 3. First lepton falsifiers

### 3.1 The first gate: can the frozen 2PN operator host a genuine intrinsic two-state sector?

Given the handoff, the first question was not about family structure or weak chirality. It was simpler and more brutal: does the current frozen ontology even contain a candidate for a **rest-frame intrinsic two-state sector**? The handoff makes that the first meaningful particle-physics gate because anything more elaborate—magnetic moments, (g=2), (4\pi) behavior—depends on having such a sector first. 

### 3.2 Why ordinary circulation or orbital motion does not count as spin

The handoff is explicit about what should fail. If a candidate “spin” variable reduces to ordinary circulation, rigid-body rotation, or some relabeling of orbital vorticity, it does **not** count as success. That matters because the earlier program already contains spin-labeled classical structures, such as the vortical far field used to reproduce weak-field Kerr/Lense–Thirring behavior,
[
v_\phi(r,\theta)=\frac{D}{r^2}\sin\theta,\qquad D=\frac{4G}{c^2}J,
]
but those are still frame-dragging/orbital-response analogs. The handoff deliberately refuses to accept that kind of object as intrinsic particle spin.  

### 3.3 The transverse dipole pair as the first serious internal candidate

Within the frozen 2PN operator, the only serious initial candidate for an internal low-energy doublet is the split odd (\ell=1) sector. The notes distinguish a degenerate transverse pair from a separate longitudinal channel:
[
R_{1\perp}=\frac72,\qquad R_{10}=4.
]
That transverse (|m|=1) pair is the first place the model naturally contains a rest-frame internal **two-plane** rather than a single scalar amplitude. This is why it became the first target of the spin audit. 

### 3.4 Why the frozen conservative operator is oscillator-like, not Pauli-like

The first serious no-go came from the operator form itself. In the minimal low-frequency completion now frozen in the 2PN notes, each channel is described by a one-pole admittance such as
[
Y_{1\perp}(\omega)=\frac{7/2}{1-\omega^2/\Omega_{1\perp}^2},
\qquad
Y_{10}(\omega)=\frac{4}{1-\omega^2/\Omega_{10}^2},
]
with analogous forms for the (P_0\oplus P_2) support channels and the geometry mode. That means the conservative 2PN throat is presently encoded as a set of **even-in-(\omega)** oscillator-like channels. This is an inference from the frozen operator form, but it is a strong one: the structure looks like internal second-order mode mechanics, not a first-order Pauli spinor sector with an intrinsic Berry or precession term already built in. 

### 3.5 Initial no-go: richer than a scalar monopole, but not yet spin-(\tfrac12)

The first falsifier therefore landed in a very specific place. The current 2PN throat is plainly richer than a scalar monopole: it already has an internal odd dipole pair, even (P_0\oplus P_2) support channels, and a separate geometry-closure channel. But at the frozen conservative level, that structure still looks like a small internal response theory built from integer-(\ell) oscillator channels, not yet like a genuine intrinsic spin-(\tfrac12) sector. So the initial conclusion of the audit was: **the model has enough structure to host an internal two-plane, but not enough yet to claim a true same-charge spin-(\tfrac12) doublet.**  

## 4. Berry-curvature and intrinsicness tests

### 4.1 Exact Berry-curvature test for collective coordinates

The handoff’s first hard filter is that any acceptable spin candidate must be a **rest-frame intrinsic structure**, not just ordinary circulation or orbital motion relabeled. That makes the right first test a symplectic one: given the parent 4D matter action, does a candidate internal coordinate actually carry a nontrivial first-order Berry/Magnus structure once center-of-mass motion is removed? The parent GNLS sector is the correct starting point because it already contains the first-order matter term
[
\mathcal L_\psi=
\frac{i\hbar}{2}\bigl(\psi^*D_t\psi-\psi(D_t\psi)^*\bigr)-\cdots,
]
and the current lepton target is explicitly the minimal charged-lepton package rather than a full PDE-first derivation.  

For any collective-coordinate family (\psi(\mathbf X;q)), the reduced first-order term is
[
L_{\rm red}^{(1)}=\mathcal A_a(q),\dot q^a,
\qquad
\mathcal A_a(q)=\frac{i\hbar}{2}\int d^4X,
\bigl(\psi^*\partial_a\psi-\psi,\partial_a\psi^*\bigr),
]
and the corresponding Berry curvature is
[
F_{ab}
======

# \partial_a\mathcal A_b-\partial_b\mathcal A_a

2\hbar,\Im!\int d^4X,\partial_a\psi^*,\partial_b\psi.
]
This gives a direct, frame-clean test of whether an internal coordinate pair really carries its own intrinsic symplectic area.

### 4.2 Result for the fixed-center internal dipole pair: zero curvature

The natural first candidate inside the frozen 2PN throat was the split odd (\ell=1) sector, because the conservative 2PN notes already show that the defect is no longer a scalar monopole but a small response theory with a carried-forward dipole wake plus a new (P_0\oplus P_2) support/closure layer. In particular, the odd sector is already split into a transverse pair and a distinct longitudinal channel, which is the first place the model naturally offers an internal two-plane.  

Take a fixed-center transverse dipole pair ((D_x,D_y)) in the simplest real (|m|=1) basis, with a common background phase,
[
\psi(\mathbf X;D_x,D_y)
=======================

e^{i\theta_0(\mathbf X)}
\Big[\Psi_0(\mathbf X)

* D_x,f(\mathbf X)\cos\phi
* D_y,f(\mathbf X)\sin\phi
* O(D^2)\Big].
  ]
  Then at the origin of the dipole plane,
  [
  \partial_{D_x}\psi=e^{i\theta_0}f\cos\phi,
  \qquad
  \partial_{D_y}\psi=e^{i\theta_0}f\sin\phi,
  ]
  so the symplectic curvature vanishes:
  [
  F_{D_xD_y}(0)
  =
  2\hbar,\Im!\int d^4X,\partial_{D_x}\psi^*,\partial_{D_y}\psi
  =0.
  ]
  The common background phase drops out, and the real transverse pair carries no intrinsic symplectic area of its own. That already tells us that the frozen conservative dipole plane is an internal oscillator plane, not yet a Pauli-like two-state sector.

### 4.3 Result for literal core translation: nonzero curvature proportional to circulation

The same test changes completely if the collective coordinates are not internal polarizations but the **literal transverse position** of a circulating core,
[
\psi_{\mathbf R}(\mathbf X)=\psi_0(\mathbf X-\mathbf R),
\qquad
\mathbf R=(X,Y).
]
Then
[
F_{XY}
======

2\hbar,\Im!\int d^4X,\partial_x\psi_0^*\partial_y\psi_0,
]
and for a circulating branch this reduces to a boundary/topological term proportional to the circulation label (\Gamma). Since the charge-like quantity in the current ontology is already tied to circulation by
[
q=\kappa_q,\rho_0,\pi a^2\Gamma,
]
the nonzero curvature in this case is exactly the Magnus/transport structure one would expect from moving a charged/circulating defect through the medium. 

### 4.4 Interpretation: the model has Magnus transport, but not yet intrinsic spin

This comparison gives a clean conceptual split. The parent 4D theory absolutely does contain gyroscopic structure, but in the present freeze it appears when one moves a **circulating object through space**, not when one isolates a fixed-center internal doublet. That is perfectly consistent with the status of the conservative program: the 1PN assembly explicitly leaves spin couplings unresolved, while the 2PN notes say only that the particle has acquired a richer dipole-plus-(P_0\oplus P_2)-plus-closure response structure.  

The conservative 2PN operator also fits that reading. Its minimal low-frequency completion is even in (\omega), so the frozen channel data look like oscillator-type response functions rather than first-order spin-precession kernels. In other words, the current files already support internal orientations and anisotropic response, but not yet a genuine first-order intrinsic spin sector.

### 4.5 Consequence: the first candidate fails the handoff’s intrinsic-spin bar

The result of the first Berry-curvature gate is therefore negative in exactly the way the handoff wanted us to be willing to find. The fixed-center transverse dipole pair is a real internal structure, but by itself it does **not** carry the intrinsic symplectic two-form a same-charge spin-like degree of freedom would need. The only nonzero gyroscopic invariant at this stage is tied to actual transport of a circulating core. By the handoff’s standard, that means the first candidate comes up short: it gives a particle-shaped hydrodynamic internal oscillator, not yet an intrinsic charged-lepton-like spin variable.  

## 5. Mixed (a!-!w) twist as the first same-charge internal candidate

### 5.1 Why the frozen parity-even wake alone cannot rescue spin

Once the fixed-center dipole pair failed the intrinsicness test, the next question was where the missing chirality could still live inside the frozen ontology. The parity-even 1PN wake sector was an obvious place to check, but it turns out not to be enough. The 1PN wake basis includes transverse, longitudinal, and helical pieces, yet the full conservative EIH match collapses the parity-even solution to the unique minimal real point
[
\alpha^2=\frac34,\qquad a_H=0,\qquad K_{\rm vec}=\frac{2}{\pi^2},
]
and the paper is explicit that the wake construction is restricted to the parity-even conservative sector. So the already-frozen wake algebra does not contain a retained helical sign that could serve as a same-charge internal spin bit.  

### 5.2 Why the mixed (A_w/F_{\mu w}/J^w) sector is the natural next place to look

The full 4+1 ontology, however, already contains a genuinely mixed sector that survives beyond the strict 3+1 Maxwell limit. The electromagnetic/plasma program keeps
[
A_w,\qquad F_{\mu w},\qquad J^w,\qquad E_w=F_{w0},\qquad C_a=F_{aw}
]
as real dynamical channels, and the reduced force law on the brane explicitly acquires the extra mixed term
[
q_s\big(E_a+(\mathbf v_s\times\mathbf B)_a+v_s^w C_a\big).
]
The same papers also emphasize that these channels are precisely what are suppressed when one takes the controlled zero-mode Maxwell/MHD limit. So the mixed (a!-!w) sector is not an invention added for the lepton audit; it is already the project’s designated place where hidden bulk structure leaks back into brane-facing dynamics. 

### 5.3 Same-charge internal handedness proposal

That led to the next live candidate: a same-charge internal **mixed-twist** degree of freedom built from the transverse odd pair together with the first odd (w)-mode. A minimal collective-coordinate ansatz is
[
\psi(b_x,b_y,\chi)
==================

e^{i\theta_0}\Big[
\Psi_0
+b_x,u_x(\mathbf x)\phi_0(w)
+b_y,u_y(\mathbf x)\phi_0(w)
+i\chi,\zeta\big(-b_yu_x(\mathbf x)+b_xu_y(\mathbf x)\big),\partial_w\phi_1(w)
\Big],
]
where (\chi=\pm1) labels the sign of the mixed (a!-!w) quadrature, (\phi_0) is the even localized mode, and (\phi_1) is the first odd one. In this reduced picture the Berry curvature becomes
[
F_{b_xb_y}
==========

2\hbar,\Im!\int d^4X,
\partial_{b_x}\psi^*\partial_{b_y}\psi
======================================

2\hbar,\chi,\zeta,\mathcal I,
]
with
[
\mathcal I \sim
\int d^3x,(u_x^2+u_y^2),
\int dw,\phi_0(w),\partial_w\phi_1(w).
]
So the mixed-twist sign (\chi) is precisely the sign of the Berry coefficient.

### 5.4 A same-charge, same-static-mass doublet is kinematically possible

The project’s Gaussian/Hermite localization machinery makes this concrete. For
[
Z(w)=e^{-w^2/\lambda^2},
\qquad
\phi_n(w)=\frac{H_n(w/\lambda)}{\sqrt{\lambda\sqrt\pi,2^n n!}},
]
odd modes decouple from a strictly brane-centered source because
[
H_{2m+1}(0)=0.
]
But the mixed derivative overlap does **not** vanish:
[
\int_{-\infty}^{\infty}
e^{-w^2/\lambda^2},\phi_0(w),\partial_w\phi_1(w),dw
===================================================

\frac{\sqrt2}{\lambda}\neq 0.
]
That is the first controlled unit-test showing that a mixed-sector Berry coefficient can survive even when the odd mode is hidden from the leading brane source. 

At the same time, the charge branch remains fixed by circulation and throat geometry,
[
q=\kappa_q,\rho_0,\pi a^2\Gamma,
]
so flipping (\chi) need not flip (q). And to this order the static density correction is (\chi)-even, since
[
|u+i\chi v|^2=u^2+v^2.
]
So ((\Gamma,\chi=+1)) and ((\Gamma,\chi=-1)) are naturally same-charge, same-static-mass candidates at the kinematic level. This was the first point in the session where a genuine **same-charge** internal sign became plausible. 

### 5.5 New no-go inside this route: scalar/even channels cannot distinguish the chirality

The mixed-twist route still contains a sharp no-go, though. The 2PN even support/closure sector is built from (P_0\oplus P_2) support channels and a separate geometry-closure direction, and in the current conservative freeze those channels couple through density/enthalpy/support variables that are **even** in the internal sign (\chi). That means the scalar bundle can renormalize support, stiffness, and closure, but it cannot by itself tell (\chi=+1) from (\chi=-1). In practical terms, the same-charge chiral bridge cannot come from the already-frozen parity-even support layer; it must come from the genuinely mixed (a!-!w) sector itself.  

So Section 5’s conclusion is mixed but important: the current ontology does contain a plausible same-charge internal chirality slot, but only in the mixed (a!-!w) sector, not in the frozen parity-even wake or scalar support bundle.

## 6. Magnetic moment and (g)-factor audit

### 6.1 Zeeman coupling test for the mixed-twist doublet

Once the same-charge mixed-twist doublet became kinematically plausible, the next gate was magnetic: does the same structure couple to an external brane magnetic field like an intrinsic dipole moment? Using the parent minimal coupling, take a weak external field in the usual axial gauge,
[
A_a^{\rm ext}=\frac12(\mathbf B\times\mathbf x)*a,
\qquad
A_w^{\rm ext}=0.
]
The linear magnetic interaction comes from the same minimal-coupling Hamiltonian already present in the parent action,
[
\delta H^{(1)}=-\int d^4X,J^A A_A.
]
Evaluated on the mixed-twist doublet, the Zeeman shift takes the form
[
\delta H_B=-\mu_z^{(\chi)}B_z,
]
with
[
S_z^{(\chi)}=\chi,\hbar,\zeta,\mathcal I,s,
\qquad
\mu_z^{(\chi)}=\chi,\frac{q\hbar,\zeta,\mathcal I}{2m*\psi},s,
]
where (s=b_x^2+b_y^2). So the same overlap that generated the Berry curvature also generates an intrinsic magnetic moment tied to the same internal sign.

### 6.2 What worked: a genuine rest-frame same-charge magnetic moment is reachable

This is a real success of the audit. The mixed-twist route does produce a **rest-frame** Zeeman coupling of the Pauli form
[
\mu_z^{(\chi)}=\frac{q}{2m_\psi},S_z^{(\chi)}.
]
So the magnetic-moment part of the handoff’s target is not empty: the toy model can reach an intrinsic dipole law scaling like (qS/m), and it can do so with the same same-charge internal label that survived the Section 5 audit. The parent minimal-coupling structure therefore contains enough room for an intrinsic magnetic response, even though it does not yet contain a full spinor sector.  

### 6.3 What failed: minimal coupling gives (g=1)-type behavior, not automatic Dirac (g=2)

What the route does **not** automatically give is the Dirac factor. At this stage the minimal relation is
[
\mu=\frac{q}{2m_\psi}S,
]
which is (g=1) relative to the parent GNLS mass parameter (m_\psi), not (g=2). This is the central negative result of the magnetic audit. The files do not contain a built-in Pauli term in the parent action, and the simple mixed-twist construction behaves like an intrinsic current loop rather than an already-Dirac object. The handoff explicitly identified (g=2) as part of the target, so on that criterion the current route is still short.  

### 6.4 Why the (v^wC_a) term is not a second independent Zeeman contribution

A tempting loophole was to try to count the mixed-sector force term (v^wC_a) as a second magnetic contribution and hope it doubles the coefficient. But the parent minimal-coupling structure does not allow that interpretation. The plasma reduction shows that the brane force law contains
[
q\big(E_a+(\mathbf v\times\mathbf B)*a+v^w C_a\big),
]
yet this mixed term is part of the force extracted from the **same** minimal interaction, not a second independent Zeeman Hamiltonian. In the point-particle analogue with (A_w=0),
[
L*{\rm int}=q,v^aA_a,
]
and the Euler–Lagrange force derived from it contains both the ordinary ((\mathbf v\times\mathbf B)_a) piece and the mixed (v^wC_a) piece. So adding a “brane current magnetic moment” and then adding a second “mixed (v^wC) magnetic moment” would double-count the same linear interaction. 

### 6.5 The remaining loophole: a mass-ratio renormalization is possible, but not derived

There is still one honest opening left in this section. The parent action itself states that the GNLS mass parameter (m_\psi) is the bulk matter mass scale and is **unrelated** to the reduced point-particle masses that appear later in the brane-facing particle description. So even though the minimal mixed-twist route gives
[
\mu=\frac{q}{2m_\psi}S,
]
the observable reduced-particle form would be
[
\mu = g_{\rm eff},\frac{q}{2M_{\rm part}},S,
\qquad
g_{\rm eff}=\frac{M_{\rm part}}{m_\psi}.
]
That means a Dirac-like (g=2) could still appear if the reduction from the parent matter scale to the observable particle scale enforces the right mass ratio. The project’s next-steps notes also leave room for geometry-coupled effective coefficients and more elaborate charged/neutral matter splits, so this loophole is structurally real. But the current files do **not** derive the needed ratio. So the magnetic-moment audit ends with a partial positive result and a sharp missing piece: intrinsic same-charge magnetic moment is reachable, automatic Dirac (g=2) is not yet.  

## 7. Induced Pauli-term and scalar-bundle bridge tests

### 7.1 Could the frozen (P_0\oplus P_2) support bundle generate the missing Pauli term?

After the mixed (a!-!w) twist route produced a same-charge internal sign and a (g=1)-type magnetic moment law, the next natural loophole was the new even 2PN response layer itself. The handoff’s magnetic target is not just “some dipole response,” but a Pauli-type intrinsic Zeeman structure with a leading Dirac target (g=2), so the right question became: can the already-frozen (P_0\oplus P_2) support/closure bundle mediate an extra (B,S) term for the same-charge internal sector?   

The 2PN notes make this a well-posed question. The conservative cross sector is already factorized into a carried-forward odd dipole wake plus a new even support bundle and a separate geometry-closure channel. In the dynamic completion, the scalar source vector still points only along the axisymmetric slots,
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
]
and the pure geometry closure enters through
[
\Delta_{\rm geom}=\frac{281}{80}.
]
So if the even bundle is going to repair the magnetic sector, it must do so through a very small scalar subspace, not through an arbitrary new operator family.  

### 7.2 Generic response formula and the scalar subspace

The bridge paper already supplies the generic response logic needed for this test. If a small internal coordinate set (q^a) is integrated out adiabatically, the effective energy shift is controlled by the quadratic response matrix (K_{ab}) and the mixed source vector (f_a):
[
H_{\rm eff}(\epsilon)
=====================

H_0+H_\epsilon\epsilon
+\frac12\left(H_{\epsilon\epsilon}-f_a(K^{-1})^{ab}f_b\right)\epsilon^2+\cdots.
]
Applying the same structure to an external magnetic source (B_n) and a candidate internal spin source (S_n), the bridge contribution takes the schematic form
[
\delta H_{BS}=-,\Lambda_B^{,T}K^{-1}\Lambda_S,B_nS_n.
]
So the entire 2PN-induced Pauli-term question reduces to whether the fixed-charge internal sign has a nonzero linear source vector (\Lambda_S) into the scalar even bundle. 

Because the current static source vector has nonzero entries only in the (P_0) and (P_{20}) directions, and the remaining negative piece is a pure geometry-closure channel, the relevant scalar subspace is just
[
Q_s=(P_0,;P_{20},;g).
]
That is already a dramatic narrowing. Whatever magnetic repair the current closure can provide must live in this three-channel scalar block, not in the full (P_2) multiplet. 

### 7.3 Strict frozen answer: the scalar bundle does not carry a same-charge chiral bridge

The strict frozen closure gives a clean negative result. The 1PN paper explicitly leaves spin couplings out of scope, and the 2PN even bundle is a support/closure sector built from scalar and quadrupole support channels, not from an explicitly chiral or parity-odd source structure. In other words, the current even bundle already knows about support, closure, and static source loading, but it does not contain an odd-even bridge term carrying a same-charge internal sign into (P_0), (P_{20}), or (g).  

The geometry side reinforces that interpretation. In the unified 4D program, the matter contribution to geometry forces enters through density-weighted confinement loading,
[
-\frac{\partial H_\psi}{\partial a}
===================================

-\int d^4X,\rho,\partial_a V_{\rm conf},
\qquad
-\frac{\partial H_\psi}{\partial L}
===================================

-\int d^4X,\rho,\partial_L V_{\rm conf},
]
so the support/closure bundle is reading the state primarily through density, confinement, and scalar response channels. Those are naturally even under a same-charge internal sign flip. The audit therefore concluded that in the **strict frozen conservative closure**, the scalar bundle can renormalize support and breathing, but it cannot distinguish the two members of a same-charge chiral pair. In that sense,
[
\Lambda_S^{(P_0,P_{20},g)}=0
]
for the closure already frozen in the papers.  

### 7.4 The two-field neutrality split remains a real loophole

The next-steps document keeps one important loophole open: charge neutrality may ultimately require either a neutralizing background or a two-field matter model with a neutral condensate plus a charged excitation field. That matters here because separating support-carrying mass from charge-carrying current is exactly the kind of move that could, in principle, change the magnetic package without rewriting the whole ontology.  

So the audit did not treat the two-field split as dead. It treated it as a legitimate future model extension. But it also recognized that the current files keep that split at the level of neutrality/source bookkeeping; they do **not** yet derive the specific derivative or nonminimal couplings that would be needed to turn it into a same-charge chiral bridge.

### 7.5 Why the minimal two-field derivative extensions still fail

The reason the first two-field fixes still fail is structural. The project’s open-system and mixed-sector files already classify (j^w), (J^w), (A_w), (F_{\mu w}), and higher-mode couplings as the channels where genuinely new odd or leakage-like structure lives. By contrast, the current scalar support and geometry sectors are density-driven and chirality-even. So the simplest scalar two-field couplings can shift mass, support, and closure, but they still do not create the missing same-charge Pauli bridge on their own.   

The audit therefore concluded that a two-field split may still be valuable later for neutrality and mass bookkeeping, but the **minimal** two-field extension does not rescue the magnetic target by itself. The missing structure is still an explicitly chiral or parity-odd internal coupling, not just a neutral-support/charged-current bookkeeping split.

### 7.6 Bottom line: current closure does not rescue the Dirac magnetic package

This part of the session ended with a real no-go. The current conservative 2PN closure, even when paired with the most obvious scalar-bundle and two-field loopholes, does **not** generate the missing Dirac magnetic package. It can support a same-charge internal sign in the mixed (a!-!w) sector and it can support a (g=1)-type intrinsic magnetic moment, but it does not produce the extra Pauli bridge needed to close the (g=2) target. That is exactly the kind of early failure the handoff asked us to accept quickly if the structure is absent.  

## 8. Pivot: step back and review the physical particle model

### 8.1 Why we stopped chasing spin directly

The reason for the pivot was methodological, not evasive. The handoff insists that the next particle-physics step should still follow the old known-destination falsification workflow, and the files themselves show why. The conservative 1PN program is frozen within a declared closure hierarchy but explicitly leaves spin couplings unresolved, while the 2PN notes say that pushing the defect harder keeps revealing genuinely new internal structure rather than simply refining earlier coefficients. That is a sign that the particle ontology itself is still incomplete.   

So instead of continuing to force a spin interpretation onto an incomplete particle model, the audit stepped back and asked a simpler question: what physical object is the model already describing, and what obvious particle properties are still missing from that description?

### 8.2 Particle picture: a 3D-facing mouth attached to a 4D corridor

The files already describe a defect as something more specific than a point particle. The carry-forward ontology says the brane sees an approximately spherical mouth/monopole face, while the interior supports cylindrical or cavity-like structure, and charge and mass are already tied to throat geometry and circulation rather than assumed fundamental. In plain language, the current particle is best thought of as a **3D-facing mouth attached to a 4D corridor**.  

The next-steps document makes that picture even more explicit by upgrading the master variables to a fully dynamic 4D bulk with a throat gate, 4+1 gauge fields, and geometry DOFs ((a,L)). So the audit was not inventing a new ontology here; it was simply consolidating what the project files already say the particle has become. 

### 8.3 Outside circulation versus inside 4D recirculation

The project already has equations for both the “outside” and the “inside” parts of the flow story, but they live in different sectors. On the brane side, there is already a wake/circulation picture through the parity-even 1PN wake and the ordinary projected current/velocity fields. In the full 4+1 theory, the missing inside piece lives in the mixed channels:
[
E_w=F_{w0},\qquad C_a=F_{aw},
]
together with the exact kinematic identity
[
\partial_A v_B-\partial_B v_A=-\frac{q}{m}F_{AB}.
]
So the “circulation after it goes in” naturally belongs to the mixed (a!-!w) part of the vorticity/gauge sector, not just to ordinary 3D vorticity on the brane.  

What is still missing is a solved moving-throat recirculation/plumbing law that closes this loop explicitly. The files already contain the tensor slots where that physics would live, but they do not yet give a frozen particle-level rule for how exterior circulation, throat intake, and 4D outflow close into one full recirculation story.

### 8.4 Hidden bulk interaction between particles

The answer to “can particles interact through hidden 4D channels we do not directly observe on the brane?” is yes in principle, but only partly in the current reduced mechanics. The localized 4+1 Maxwell sector, mixed fields (A_w) and (F_{\mu w}), transverse currents (J^w), leakage terms, and the Hermite/KK tower all already exist in the full ontology. Those channels generate Yukawa-like static corrections, mixed force terms, leakage, and extra work/energy reservoirs that a 3D observer would only see indirectly after projection.   

But the conservative 1PN point-particle reduction intentionally closes before those channels are allowed to dominate. So the hidden bulk interaction story is already present in the exact theory, but not yet reduced to a standard two-particle brane law.

### 8.5 The brane is a filter, not a brick wall

One of the clearest conceptual gains of the audit was the recognition that the brane is best viewed as a **projection filter** on a bulk system, not as a literal wall. The unified 4D papers insist on the distinction between projection and reduction:
[
\mathcal P_W[Q](t,\mathbf x)
============================

\int W(w),Q(t,\mathbf x,w),dw,
]
and the exact projected continuity law is open,
[
\partial_t \overline{\rho_s}+\partial_a\overline{j_s^a}=S_{\rm leak}^{(s)}.
]
So even if bulk evolution is conservative, the brane subsystem is generically open. That is why 4D behavior can alter 3D motion without the brane observer ever seeing a local 3D source for it.   

### 8.6 The (+w) versus (-w) puncture question

The user’s question about a throat opening toward (+w) or (-w) turned out to isolate a real missing variable. In the centered Gaussian localization picture, the odd Hermite modes satisfy
[
H_{2m+1}(0)=0,
]
so they decouple from a strictly brane-centered source. In the standard zero-mode brane reduction, one also suppresses (A_w), (F_{\mu w}), and (J^w). Put together, that means the present frozen reduction largely projects away “which side of the brane the throat opens toward.” The current files have the right channels for such a distinction to matter, but they do not yet freeze an explicit side/orientation variable.   

### 8.7 What the current ontology still does not explicitly represent

The particle audit therefore identified three missing physical layers. First, the model lacks an explicit **side/orientation variable** distinguishing (+w) from (-w) punctures. Second, it lacks a fully solved **recirculation plumbing** law that tracks how exterior intake, interior bulk transport, and back-pressure close into one particle-scale loop. Third, it lacks a reduced **bulk pair kernel** telling us how two actual defects interact through the mixed (A_w/F_{\mu w}/J^w) sector. The next-steps document also makes clear that the full dynamic model will need those pieces anyway if it is to handle throat dynamics, back-pressure, and bulk/brane exchange honestly.   

## 9. Oriented particle package

### 9.1 New particle labels introduced in the audit

The session’s answer to the missing-particle-properties problem was to define a minimal **oriented particle package** before going back to spin. In the audit, a single particle state was enlarged schematically to
[
\mathcal X_i=
\big(\mathbf X_i,\mathbf V_i,\ a_i,L_i,\Gamma_i,\eta_i,\tau_i,\ Q_i^\alpha\big),
]
where ((a,L,\Gamma,Q^\alpha)) are already present in the project’s geometry/circulation/response ontology, while
[
\eta_i=\pm1,\qquad \tau_i=\pm1
]
are **new proposed labels**, not yet frozen in the project files. The files already contain the structures they would have to live in—dynamic 4D bulk variables, 4+1 gauge fields, projection kernels, and geometry DOFs—but not the labels themselves.  

The idea is to separate two things the old particle model was conflating: a side/orientation distinction of the throat itself, and a same-charge internal twist/chirality distinction living in the mixed (a!-!w) sector.

### 9.2 How the side label (\eta=\pm1) would be wired in

The cleanest proposed implementation of the side bit was an **odd-(w)** deformation of objects that are already part of the frozen ontology, such as the confinement potential or the EM localization profile:
[
V_{\rm conf}^{(\eta)}
=====================

V_{\rm even}
+\eta,V_{\rm odd},
\qquad
Z^{(\eta)}(w)
=============

Z_{\rm even}(w)
+\eta,Z_{\rm odd}(w),
]
with the odd pieces changing sign under (w\mapsto -w). This gives the throat an orientation relative to the brane without changing its leading even brane-facing monopole character.   

That proposal is also consistent with the current blindness of the centered reduction. In the zero-mode/centered regime, odd modes are largely projected away from direct brane observables, so an (\eta)-odd structure can exist in the bulk without automatically ruining the carry-forward 1PN/2PN far-field behavior. The audit therefore treated (\eta) as a plausible extension of the present ontology rather than a contradiction of it.

### 9.3 How the same-charge twist label (\tau=\pm1) would be wired in

The proposed same-charge internal twist label (\tau) lives in a different slot: the mixed (a!-!w) vorticity/gauge structure already present in the 4+1 theory. The key exact identity is
[
\partial_A v_B-\partial_B v_A=-\frac{q}{m}F_{AB},
]
so in particular the mixed component
[
\Omega_{aw}\equiv \partial_a v_w-\partial_w v_a
]
is directly tied to the mixed field strength (F_{aw}). The audit therefore proposed that a same-charge internal sign could be built from a transport-subtracted mixed-twist moment, schematically
[
\Xi_i
\sim
n_i^a\int_{\mathcal T_i}\rho_i,\widetilde{\Omega}_{aw}^{(i)},d^4X,
\qquad
\tau_i=\operatorname{sign}(\Xi_i).
]
That would make (\tau) a rest-frame internal label attached to mixed brane–bulk twist rather than to overall charge sign. The exact formula is still a proposal from the audit, but the tensor slot it uses is already present in the project files.  

### 9.4 Bulk pair kernel idea

Once (\eta) and (\tau) exist, the natural next reduced object is a bulk pair kernel describing how two oriented/thickened defects exchange through the hidden sector. At the parent level the obvious schematic form is a 5D current-exchange action plus the open-system correction terms already emphasized in the 4D program:
[
\Delta S_{AB}^{\rm bulk}
========================

\frac12\int d^5x,d^5x',
J_A^M(x),G^{(5)}*{MN}(x,x'),J_B^N(x')
+\Delta S*{AB}^{\rm leak/cov}.
]
After reduction, the audit proposed a symmetry-allowed kernel of the form
[
\Delta L_{AB}^{\rm bulk}
========================

K_0
+\eta_A\eta_B,K_{\eta\eta}
+\tau_A\tau_B,K_{\tau\tau}
+(\eta_A\tau_B+\eta_B\tau_A),K_{\eta\tau}
+\cdots.
]
This is not yet a file-derived equation, but it is the natural reduced pair law suggested by the existing localized 4+1 Maxwell, KK tower, and mixed-sector force structure.   

### 9.5 Why the oriented particle package comes before the lepton target

The purpose of the oriented particle package is not to replace the minimal charged-lepton target; it is to make it testable. The handoff’s destination still requires a genuine same-charge internal two-state sector, an intrinsic magnetic moment, and possibly a (2\pi/4\pi) identity signal. But after the first audit, it became clear that the current particle ontology is missing the side/orientation and mixed-twist ingredients needed even to state those questions cleanly. So the audit treated the oriented particle package as a prerequisite layer:
[
\text{oriented particle package}
;\longrightarrow;
\text{minimal charged-lepton package}.
]
Only after the particle model carries explicit side/orientation and mixed-twist structure does it make sense to return to the lepton audit and ask whether any of those new labels survive the handoff’s intrinsicness test.  

## 10. Five checks for the oriented particle package

### 10.1 Far-field neutrality / invisibility — pass, with a condition

The first check was whether adding an oriented particle layer would automatically spoil the already-frozen far-field brane physics. The current files say it need not. The next-steps program already treats charge neutrality as a blocking requirement and keeps brane observables as projections of bulk fields rather than stitched boundary values. In the same framework, the Gaussian/Hermite localization tower shows that odd transverse modes decouple from a strictly brane-centered source because odd Hermite modes vanish at (w=0). That means a side/orientation structure can, in principle, live in an odd-(w) sector while remaining invisible to the leading even brane projection. In that narrow sense, the model already has a mathematically clean hiding place for a side bit without wrecking the frozen 1PN/2PN far-field picture.   

### 10.2 Same-charge two-state separation — open, but not contradicted

The second check was whether the ontology has room for a same-charge internal pair rather than merely a particle/antiparticle sign flip. The handoff’s target is explicit: a viable lepton-like sector must contain a genuine internal two-state degree of freedom at fixed charge branch, not just reversed circulation relabeled as spin. The current project state does not yet derive such a pair, but it also does not forbid it. Charge is still carried by the circulation/geometry dictionary
[
q=\kappa_q,\rho_0,\pi a^2\Gamma,
\qquad
m_G=\kappa_m,\rho_0,\pi a^2L,
]
while the full 4+1 ontology already contains separate mixed channels (A_w), (F_{\mu w}), and (J^w) that are not exhausted by (\Gamma) alone. So a same-charge internal sign is still kinematically plausible, but it is not yet frozen as an actual derived branch of the model.   

### 10.3 Recovery of the current frozen model when the new pieces are turned off — pass

This check came out strongly positive. The current 4D/1PN reduction already identifies a controlled brane Maxwell limit in which
[
A_w=0,\qquad \partial_wA_\mu=0,
]
and the remaining brane field obeys the standard reduced Maxwell equation with (\mu_0^{\rm eff}=\mu_0/Z_{\rm int}). The EM summary gives the same message from the localized-Maxwell side: once (w)-dependence and transverse current are suppressed, one recovers the standard brane sector and the odd tower disappears from direct brane sourcing. So an oriented-particle extension that vanishes when its odd-(w) and mixed-sector pieces are set to zero is a true extension of the frozen model, not a rewrite of it.  

### 10.4 Hidden bulk exchange already has explicit ledgers — pass strongly

The fourth check was whether hidden 4D interaction would be “mysterious” if it existed. It is not. The 4D/plasma program already gives exact projected continuity with leakage, explicit mixed force terms, and named diagnostics for flow support and back-pressure. In particular, the next-steps file already promotes through-flow and back-pressure to concrete observables,
[
\Phi_w(t)=\int_{\Sigma_w}J_{{\rm mass},w},d^3S,
\qquad
\Delta h(t)=\langle h(\rho)\rangle_{\rm near\ mouth}-\langle h(\rho)\rangle_{\rm far\ bulk},
]
while the plasma ledger keeps mixed-sector work channels such as (J^wE_w) and leakage terms (S_{\rm leak}^{(s)}). So if two oriented defects really exchange through the hidden sector, the current ontology already has exact places where that difference must appear.   

### 10.5 Non-orbital rest-frame survival — not yet passed

This is the only check that still blocks a return to the lepton target. The handoff is explicit that the next candidate must survive the “intrinsic versus orbital” filter: a valid spin-like variable must remain meaningful in the defect rest frame after ordinary circulation is removed. The 1PN conservative assembly explicitly leaves spin couplings unresolved, and the 2PN notes, while much richer, still only guarantee a carried-forward dipole wake plus a new (P_0\oplus P_2) support layer and a separate geometry-closure channel. That is enough to justify a deeper particle model, but not enough yet to claim a rest-frame same-charge two-state invariant.   

### 10.6 Audit result

So the oriented particle package passed the particle-ontology checks more strongly than it passed the lepton checks. Far-field invisibility is viable, the frozen model is recovered when the odd/mixed sector is switched off, and hidden bulk exchange already has exact ledgers. Same-charge internal separation remains open rather than blocked. But non-orbital rest-frame survival is still unproved, which is why the session next turned to overlap-integral tests for the new odd/mixed structures rather than directly claiming progress on spin.  

## 11. Overlap-integral tests for the new particle ingredients

### 11.1 Centered-throat parity structure

The first overlap test used the standard centered reduction that is already frozen in the localized-Maxwell sector. For Gaussian localization,
[
Z(w)=e^{-w^2/\lambda^2},
\qquad
\phi_n(w)=\frac{H_n(w/\lambda)}{\sqrt{\lambda\sqrt\pi,2^n n!}},
]
and the brane coupling of a strict (\delta(w)) source is controlled by (f_n(0)). Because odd Hermite modes satisfy
[
H_{2m+1}(0)=0,
]
they decouple from a strictly centered brane source. That is the exact parity fact behind the session’s conclusion that the strict centered baseline has no visible side/orientation label built into its leading brane observables.  

### 11.2 Why the strict baseline does not already contain a side bit

Once that parity structure is taken seriously, the current centered baseline cannot by itself distinguish “opens toward (+w)” from “opens toward (-w)” in any leading projected scalar observable. That is not because the physics is impossible; it is because the standard reduction was intentionally built to preserve the centered brane-compatible sector and suppress odd-(w) information. This is why the session treated the side label (\eta=\pm1) as a **new proposed variable** rather than as an overlooked quantity already frozen in the files. The smallest extension consistent with the current logic is therefore an odd-(w) deformation of the already-existing confinement or localization structures, for example
[
V_{\rm conf}^{(\eta)}=V_{\rm even}+\eta,V_{\rm odd},
\qquad
Z^{(\eta)}=Z_{\rm even}+\eta,Z_{\rm odd},
]
with (V_{\rm odd}) and (Z_{\rm odd}) odd under (w\mapsto -w). This is a proposal from the session, not a file-derived frozen equation. It is motivated by the parity facts above.  

### 11.3 Mixed derivative overlaps of the first odd mode

The key discovery of the overlap analysis was that “odd mode invisible to direct brane coupling” does **not** mean “odd mode irrelevant.” In the same Gaussian/Hermite basis, the direct even–odd scalar overlap vanishes by parity, but the mixed derivative overlap does not. In the session’s normalized unit test,
[
\int_{-\infty}^{\infty}e^{-w^2/\lambda^2},\phi_0(w),\partial_w\phi_1(w),dw
==========================================================================

\frac{\sqrt2}{\lambda}\neq 0.
]
That calculation used only the frozen Hermite tower and its weighted inner product. So the first odd mode is hidden from a direct brane source, but not from the derivative structures that appear naturally in the mixed (a!-!w) sector. 

### 11.4 What this means physically

This was the decisive reason the session stayed interested in the mixed sector. In the strict zero-mode Maxwell limit,
[
A_w\approx 0,\qquad \partial_wA_\mu\approx 0,\qquad J^w\approx 0,
]
the mixed field (F_{aw}) vanishes and the candidate fixed-center twist invariant collapses. But once the first odd (w)-mode and the mixed (F_{aw})-type channel are kept alive, the centered throat can support a nonzero same-charge odd coordinate without immediately becoming visible as a direct brane scalar. In other words: centered geometry kills the direct scalar side label, but it does **not** kill a mixed-sector odd internal coordinate.  

### 11.5 The next diagnostic promoted by this overlap result

The overlap analysis therefore promoted a very specific linear test for the new odd variable. In the session’s notation, the right quantities to check were:

[
\Xi'(0),\qquad
\partial_b q\big|*{0},\qquad
\partial_b m\big|*{0},
]
where (b) is the centered odd amplitude, (\Xi) is the transport-subtracted mixed-twist observable, and (q,m) are the particle branch charge and mass. This was the exact falsifier that then fed into Section 12: if (\Xi'(0)) vanishes, the odd mode is irrelevant; if (\partial_b q|_0\neq 0), it is not same-charge; if (\partial_b m|_0\neq 0), it is not a clean internal branch. Those conditions were not lifted from a prior paper; they were the session’s reduced-model diagnostic distilled from the file-based parity and projection structure. 

## 12. Same-charge odd mode: success and limitation

### 12.1 Linear same-charge test — passes

The first outcome of that diagnostic was positive. In the session’s reduced centered-throat analysis, the mixed derivative overlap made it possible for the transport-subtracted internal twist observable (\Xi) to acquire a nonzero linear slope once the first odd mode was retained:
[
\Xi'(0)\neq 0.
]
At the same time, the charge and mass dictionaries stay fixed at linear order so long as the circulation branch (\Gamma) and the background geometry ((a,L)) are held fixed:
[
q=\kappa_q,\rho_0,\pi a^2\Gamma,
\qquad
m_G=\kappa_m,\rho_0,\pi a^2L.
]
So the centered odd mode passed the first same-charge test in the narrow linear sense:
[
\partial_b q\big|_0=0,
\qquad
\partial_b m\big|_0=0.
]
That is the clearest positive result of this part of the session: the model can carry a **same-charge nonorbital odd coordinate** once the mixed sector is retained.  

### 12.2 Why that still falls short of a genuine two-state sector

But the success was only linear and kinematic. The conservative 2PN throat still encodes the odd sector as a stable response theory, not a bifurcated internal branch. The 2PN notes freeze positive odd residues
[
R_{1\perp}=\frac72,\qquad R_{10}=4,
]
and package the corresponding channels into pole-completed kernels such as
[
Z_{1\perp}(\omega)=\frac{1-\omega^2/\Omega_{1\perp}^2}{7/2},
\qquad
Z_{10}(\omega)=\frac{1-\omega^2/\Omega_{10}^2}{4}.
]
That is the algebra of a stable oscillator sector. So the same odd coordinate that passed the linear same-charge test is still a **continuous amplitude** (b), not yet a discrete (\pm) pair or a true low-energy doublet.  

### 12.3 What the result means

This gave the session a very useful narrowing. The model does contain a candidate same-charge internal coordinate that is not immediately reducible to overall circulation. But it does **not** yet contain the mechanism that would turn that coordinate into a genuine two-state internal sector of the kind the handoff requires for the minimal charged-lepton package. In the language of the audit, the model has crossed from “no internal same-charge sign at all” to “same-charge internal odd coordinate exists,” but it is still short of “intrinsic two-state particle identity.”  

### 12.4 The missing ingredient exposed by this limitation

By the end of this step, the missing ingredient had become much more specific. The next candidate structure had to be something that **locks or quantizes** the odd coordinate without collapsing it back into ordinary orbital flow. That is why the subsequent parts of the session moved to two possible mechanisms: nonlinear wall/geometry bifurcation of the odd branch, and a mixed-sector Berry rotor that might quantize the surviving transverse internal plane. The point of Sections 10–12 is that neither of those later moves would have made sense without first isolating the result recorded here: the model already supports a same-charge nonorbital odd coordinate, but not yet a protected two-state sector.  

## 13. Family-1 wall/geometry search for sign locking

### 13.1 Why the Family-1 wall sector became the next target

By the time the audit reached this stage, the model had already passed one important test and failed another. It had produced a plausible same-charge, nonorbital **odd internal coordinate** in the mixed (a!-!w) sector, but that coordinate was still continuous rather than a genuine low-energy (\pm) doublet. The natural next place to look was therefore the throat’s own support and geometry layer: could the wall/support structure turn that continuous odd amplitude into a sign-locked branch? This was not a random shift in focus. The 2PN notes had already reorganized the particle into a carried-forward dipole wake plus a new (P_0\oplus P_2) support bundle and a separate geometry-closure channel, while the next-steps/inner-throat notes explicitly identified the missing layer as a true inner-throat/wall-response reduction rather than another abstract coefficient fit.   

In other words, once the mixed-sector audit revealed “same-charge internal sign exists, but is not yet discrete,” the Family-1 wall/geometry sector became the first honest place where a (\pm) bifurcation could emerge **without** abandoning the throat ontology the earlier papers had already frozen.

### 13.2 The support-minus-closure picture at 2PN

The 2PN notes now support a very concrete structural reading of the throat. On the even side, the conservative static cross sector is no longer just a list of coefficients; it is a positive support bundle carried by (P_0\oplus P_2) ports plus a separate geometry-closure deficit. In the canonical real basis, the scalar source vector is
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
]
while the pure geometry closure appears through
[
\Delta_{\rm geom}=\frac{281}{80}.
]
The same notes also show that the static support operator can be written in a compact Family-1 wall form with local profiles
[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
z_{\rm curv}(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4,
]
and a source profile
[
S(\mu)=10-\frac{63}{2}z_{\rm base}(\mu)=\frac{7}{16}+\frac{45}{16}\mu^2.
]
So the static 2PN throat is already a genuine low-order Family-1 flare/support theory in the basis ({1,\mu^2,\mu^4}), equivalently (P_0\oplus P_2) at the level relevant here. 

That support-minus-closure structure is the reason the wall sector became so attractive. It is already a small, positive, physically interpretable response bundle with one distinguished negative closure direction. If any conservative sign-locking mechanism is going to appear before we invoke a genuinely topological/Berry term, it should show up as a nonlinear interaction between the odd mode and this support/closure layer.

### 13.3 The search for a double-well / quartic softening mechanism

The concrete question was whether a centered odd amplitude (b) could soften the Family-1 wall sector strongly enough to generate a double-well. Because the centered background is even under (b\to -b), the first allowed couplings into the even bundle are quadratic in (b). That leads to the minimal reduced Landau system
[
E(b,Q,g)
========

\frac12\kappa_b b^2+\frac14 u_b b^4+\frac16 v_b b^6
+\frac12 Q^T K_Q Q+\frac{g^2}{2Y_g(0)}
-b^2(\lambda_Q^TQ+\lambda_g g)+\cdots,
]
where (Q) is the even support bundle, here minimally the scalar pair ((q_0,q_{20})), and (g) is the monopole geometry/breathing channel. Eliminating the even variables gives
[
E_{\rm eff}(b)=\frac12\kappa_b b^2+\frac14 u_{\rm eff} b^4+\frac16 v_{\rm eff} b^6+\cdots,
\qquad
u_{\rm eff}=u_b-2,\lambda^T K^{-1}\lambda.
]

This was the first real nonlinear wall/geometry falsifier of the session. If (u_{\rm eff}) never goes negative on a physical Family-1 branch, then the wall route dies cleanly. If it does, the wall sector at least opens the door to a (\pm) branch. The appeal of this route was that it reuses the exact 2PN support/closure hierarchy the notes had already frozen rather than adding a new particle structure by hand.

### 13.4 Reduction to a small set of coefficients

The structure simplifies dramatically because the Family-1 support and source are already locked together. If the odd amplitude shifts the base wall profile as
[
\delta z_{\rm base}(\mu)=b^2\bigl[\alpha_0+\alpha_2P_2(\mu)\bigr]+O(b^4),
]
then the frozen source relation
[
S(\mu)=10-\frac{63}{2}z_{\rm base}(\mu)
]
forces
[
\delta S(\mu)=-\frac{63}{2}b^2\bigl[\alpha_0+\alpha_2P_2(\mu)\bigr].
]
Projecting onto the canonical scalar ports gives the induced couplings
[
\lambda_0=-\frac{21}{\sqrt5}(3\alpha_0+\alpha_2),
\qquad
\lambda_{20}=-21\alpha_2.
]

So the whole scalar support bridge collapses to just two coefficients ((\alpha_0,\alpha_2)). The geometry side contributes one more coefficient (\lambda_g) through the breathing channel, and the bare odd branch contributes the self-interaction coefficients (u_b) and (v_b). If one also tracks the first (b^2)-shift of the static monopole stiffness, the remaining geometry input is
[
\kappa_{00}^{(2)}\equiv \partial_{b^2}K_{00}\big|*0.
]
That is how the original “large unknown nonlinear wall sector” reduced to the small coefficient set
[
\alpha_0,\ \alpha_2,\ \lambda_g,\ \kappa*{00}^{(2)},\ u_b,\ v_b.
]
The 2PN notes already fix the static monopole geometry susceptibility to
[
Y_g(0)=\frac{109}{280},
]
so once these coefficients are known the quartic and sextic shifts are fully determined. 

### 13.5 What the known 2PN data already fix, and what remains open

The 2PN files already fix a large part of the wall/geometry background against which the odd branch must be tested. On the odd side the conservative residues are positive,
[
R_{1\perp}=\frac72,\qquad R_{10}=4,
]
and the explicit axisymmetric wall/DtN scaffold gives positive odd poles,
[
z_{1\perp}\approx 2.561,\qquad z_{10}\approx 2.531.
]
On the even side, the support/source/closure data are already frozen through the Family-1 profiles quoted above, the source vector (J), and the geometry-closure deficit (\Delta_{\rm geom}=281/80). The monopole geometry channel is also positive/passive on the worked Family-1 branch, with a low-frequency closure dominated by one breathing mode. 

What the current files do **not** yet give are the actual odd-deformed Family-1 coefficients themselves. There is no frozen derivation yet of
[
\alpha_0,\ \alpha_2,\ \lambda_g,\ \kappa_{00}^{(2)},\ u_b,\ v_b
]
because that would require solving the odd-deformed throat branch rather than only the even/static support problem. So the wall route was narrowed to a handful of coefficients, but not numerically closed.

### 13.6 Key threshold logic: support softening versus geometry stabilization

Even before those coefficients are fully known, the 2PN structure already tells us how the competition must work. The quartic correction is always negative-semidefinite:
[
u_{\rm eff}
===========

u_b-\frac{2646}{5}\bigl(3\alpha_0^2+2\alpha_0\alpha_2+2\alpha_2^2\bigr)-\frac{109}{140}\lambda_g^2.
]
So the Family-1 support/geometry sector can only **soften** the odd mode; it cannot stiffen it. For the transverse branch later identified as the more relevant same-charge candidate, one finds (\alpha_2=-\alpha_0), so
[
u_{\rm eff}^{(|m|=1)}
=====================

u_b-\frac{7938}{5}\alpha_0^2-\frac{109}{140}\lambda_g^2.
]
This immediately shows that support-side softening is much more efficient than pure geometry breathing at driving the quartic negative.

But the known odd quadratic backbone stays positive. Because the current files keep the odd poles positive, the centered odd branch has
[
\kappa_b>0.
]
So quartic softening alone cannot produce the clean continuous pitchfork that one would ideally want for an intrinsic (\pm) sector. At best, if (u_{\rm eff}<0) and (v_{\rm eff}>0), the wall route can generate **outer minima** through a first-order/metastable mechanism. Geometry breathing then naturally looks more useful as a **sextic stabilizer** than as the driver of the bifurcation itself. This is why the wall route remained interesting but never became a clean “spin solved” result in the session. 

## 14. Why the axisymmetric wall route comes up short

### 14.1 The transverse (|m|=1) branch is the right candidate for a same-charge internal plane

Once the wall-softening question was posed, it mattered which odd branch was being softened. The 2PN notes split the odd sector into a transverse pair and a distinct longitudinal channel. That split is already present in the solved residues
[
R_{1\perp}=\frac72,\qquad R_{10}=4.
]
The transverse (|m|=1) pair is the only odd sector that already provides an internal **two-plane**. The axisymmetric (m=0) odd branch is only one real amplitude. So if the goal is a same-charge internal degree of freedom that could later be tested against the Pauli/Dirac package, the transverse pair is the right live candidate. The (m=0) branch remains useful as an instability probe, but by itself it can never be more than an Ising-like sign choice.  

### 14.2 The current axisymmetric wall package depends only on (s=b_x^2+b_y^2)

Write the transverse odd branch in real amplitudes
[
\mathbf b=(b_x,b_y),\qquad s=b_x^2+b_y^2.
]
Because the static Family-1 wall/source package is axisymmetric, all of its explicit input depends only on (\mu=\cos\theta), not on the azimuth (\phi). The same is true of the solved source vector:
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
]
which carries no (P_{21}) or (P_{22}) drive. So after integrating out the even channels, the static reduced energy of the transverse pair must be an (O(2))-invariant function of (s), not of the angle in the ((b_x,b_y)) plane. In the reduced language of the audit,
[
E_{\rm stat}(b_x,b_y)=E_{\rm stat}(s).
]
That is already enough to predict trouble for any attempt to get a discrete same-charge doublet purely from the current wall law. 

### 14.3 Including (P_{22}) support does not help, because (|Q_2|^2=s^2)

A possible loophole was the existence of the real (m=\pm2) support channels in the 2PN operator. The transverse odd pair necessarily backreacts into a real quadrupole
[
Q_2=(b_x^2-b_y^2,\ 2b_xb_y).
]
So one might hope that integrating out the (P_{22}) channels would create a genuine anisotropy. But the algebra closes against that hope:
[
|Q_2|^2
=======

# (b_x^2-b_y^2)^2+(2b_xb_y)^2

# (b_x^2+b_y^2)^2

s^2.
]
So the (P_{22}) channels do appear, but in the current axisymmetric support package they contribute only another (s^2) term. They soften or stiffen the **radius** of the transverse branch; they do not select an angle in the plane. This is why the session eventually concluded that (P_{22}) is a passive support channel in the current freeze, not a genuine anisotropy source. The source vector has no (P_{22}) entries, so there is no static (\cos2\phi)-type term to split the rotor. 

### 14.4 Outcome: if the branch condenses, it condenses into a ring, not a discrete pair

This symmetry result forces a very specific outcome. If higher-order wall softening ever does create a nonzero transverse radius (s_*), the minima are not a pair
[
b=\pm b_*,
]
but a whole circle
[
b_x^2+b_y^2=s_*.
]
So the current axisymmetric Family-1 wall/geometry package can produce, at best, a **transverse rotor**. That is already richer than a scalar monopole, but it is still not the handoff’s required “genuine internal two-state sector.” A ring of minima is a continuous internal orientation manifold, not a protected same-charge doublet.  

### 14.5 The (m=0) branch can make an Ising-like (\pm) pair, but not a Pauli-like doublet

The axisymmetric (m=0) odd branch behaves differently. Because it is one real amplitude, if its quartic softening is strong enough it can generate a true pair
[
b=\pm b_*.
]
And in the wall audit it actually softens more strongly than the transverse branch. But even in the best case this only yields an **Ising-like sign label**. It does not produce an internal two-plane, and it carries no natural route to a Pauli spinor structure. So the axisymmetric branch may later become a useful extra internal quantum number, but it cannot by itself satisfy the handoff’s minimal charged-lepton target.  

## 15. Last major fork: static side bit versus Berry route

### 15.1 A static (+w/-w) side label by itself cannot split the transverse rotor

Once the axisymmetric wall route stalled at a rotor, the next question was whether a throat “side” label—opening preferentially toward (+w) or (-w)—could statically break the degeneracy. The answer inside the current freeze is no. A side label (\eta=\pm1) certainly breaks (w\to -w), but it is still a **scalar under in-plane (O(2))**. So in the current axisymmetric throat package, it can only enter static terms like
[
\eta,s,\qquad \eta^2 s,\qquad \eta,s^2,\ldots
]
which change the radial potential but do not pick an angle in the ((b_x,b_y)) plane. Since the static wall/source package depends only on (\mu) and the current centered reduction suppresses odd modes from direct brane sourcing, a side bit alone cannot convert the transverse ring into a discrete (\pm) pair.  

### 15.2 Why a genuine (m=2) anisotropy source would be needed

To split the transverse rotor **statically**, one needs a real in-plane quadrupolar background tensor,
[
\Sigma_2=(\Sigma_{2c},\Sigma_{2s}),
]
entering the reduced energy as
[
\delta E_{\rm aniso}
====================

# -\lambda_{22},\Sigma_2!\cdot!Q_2

-V_2\cos2(\phi-\phi_0),
]
where
[
Q_2=s_*(\cos2\phi,\sin2\phi).
]
That is the minimal way to break the in-plane (O(2)) down to a discrete twofold choice. Without such an (m=2) source, the static energy of the transverse branch can only depend on (s). This is why the session’s fork became so sharp: either find a real (P_{22}) drive, or abandon the static side-bit rescue and follow the mixed-sector Berry route instead.

### 15.3 No file-derived (P_{22}) drive exists in the current freeze

The session then checked whether such a source was already hiding in the frozen files. It is not. The even 2PN operator does contain (P_{22c}) and (P_{22s}) as support channels, but the actual source data are still purely axisymmetric:
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right).
]
The Family-1 support/source law is equally axisymmetric:
[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
z_{\rm curv}(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4,
\qquad
S(\mu)=\frac{7}{16}+\frac{45}{16}\mu^2,
]
and even the minimal strict PDE-side repair identified in the notes remains axisymmetric,
[
B_{\rm tan}(\mu)=B_0+B_2\mu^2.
]
So the (P_{22}) slots are present, but only as passive support capacity, not as an active anisotropy drive. There is no frozen (m=2) background source that could statically split the transverse rotor. 

### 15.4 So the static side-bit route fails as the main rescue

This gave the session its last major bifurcation. The static side/orientation route failed cleanly inside the current conservative freeze: a pure (+w/-w) side label does not split the transverse rotor, and the current files do not already contain the (P_{22}) anisotropy source that would be required. That left one live corridor: the **mixed-sector Berry route**. The reason it stayed alive is structural. The parent matter action is already first order in time, and the full 4+1 Maxwell/plasma program keeps the mixed channels
[
A_w,\qquad F_{\mu w},\qquad J^w,\qquad E_w,\qquad C_a=F_{aw}
]
as genuine physical sectors beyond the strict zero-mode Maxwell limit. Those are exactly the kinds of channels that can supply a first-order internal phase law even when the static axisymmetric wall package cannot. So the session’s conclusion at this fork was precise: the static side-bit rescue fails, and any remaining route to a same-charge internal doublet must come from a mixed-sector Berry or topological mechanism rather than from the present conservative wall law alone.   

## 16. Mixed-sector Berry rotor

### 16.1 Why the parent theory allows a first-order internal symplectic term

Once the static side-bit and purely conservative wall routes had been pushed as far as they could go, the last genuinely live corridor was the mixed (a!-!w) sector. The reason is structural. The parent matter theory is already first order in time,
[
\mathcal L_\psi=
\frac{i\hbar}{2}\bigl(\psi^*D_t\psi-\psi(D_t\psi)^*\bigr)-\cdots,
]
so a collective-coordinate Berry/Magnus term is kinematically allowed from the beginning. At the same time, the 4+1 Maxwell/plasma sector explicitly retains the mixed channels
[
A_w,\qquad F_{\mu w},\qquad J^w,\qquad E_w=F_{w0},\qquad C_a=F_{aw}
]
outside the strict zero-mode reduction, and the bulk vorticity identity ties those mixed fields directly to matter kinematics through
[
\Omega_{ij}=-\frac{q}{m}F_{ij}.
]
So the ontology already contains the exact ingredients needed for a first-order internal rotor; they simply disappear in the special 3+1 Maxwell limit where the mixed sector is turned off.   

### 16.2 Why the mixed (a!-!w) channel can generate a nonzero Berry coefficient

The localized 4+1 Maxwell analysis also explains why this route was invisible for so long. In the Gaussian localization tower,
[
Z(w)=e^{-w^2/\lambda^2},
\qquad
f_n(w)=H_n(w/\lambda),
\qquad
m_n^2=\frac{2n}{\lambda^2},
]
odd Hermite modes obey
[
H_{2m+1}(0)=0,
]
so they do not couple to a strictly brane-centered source at leading order. That is why the zero-mode brane Maxwell limit suppresses the very channels we need. But this is a parity statement about direct sourcing, not a theorem that the odd tower cannot appear in mixed derivative overlaps. Once the mixed (F_{aw})-type channel is retained, the first odd mode becomes a legitimate internal quadrature partner for the transverse (|m|=1) mouth pair rather than disappearing from the problem.  

### 16.3 Minimal collective-coordinate model of the transverse mixed rotor

The session’s minimal live collective-coordinate ansatz was therefore built from the transverse odd pair (b_x,b_y) together with the first odd (w)-mode and a same-charge chirality sign (\chi=\pm1). In reduced form, the most general axisymmetric low-energy Lagrangian is
[
L_{\rm tr}
==========

\frac{M_b}{2}\dot b_i\dot b_i
+
\frac{\chi,\nu_0}{2},\epsilon_{ij}b_i\dot b_j
---------------------------------------------

V_{\rm stat}(b_x^2+b_y^2),
]
where (M_b) is the inertial coefficient of the transverse odd branch and (\nu_0) is the mixed-sector Berry coefficient. The physical point is that the sign (\chi) need not be tied to the circulation/charge branch (\Gamma). Charge in the current ontology still lives on the older geometric/circulation dictionary,
[
q=\kappa_q,\rho_0,\pi a^2\Gamma,
]
so flipping (\chi) can, in principle, leave (q) unchanged. This is the first point in the whole audit where a same-charge internal sign and a genuine first-order internal symplectic term coexist in one reduced model.   

### 16.4 If the transverse radius pins, the mixed sector becomes an internal rotor

If the static wall/support sector ever succeeds in pinning the radius
[
s_* = b_x^2+b_y^2,
]
then the transverse odd pair reduces to a rotor:
[
b_x=\sqrt{s_*}\cos\phi,\qquad
b_y=\sqrt{s_*}\sin\phi.
]
Substituting this into the reduced Lagrangian gives
[
L_{\rm rot}
===========

\frac{I}{2}\dot\phi^2+\chi,\kappa_0,\dot\phi,
\qquad
I=M_b s_*,
\qquad
\kappa_0=\frac{\nu_0 s_*}{2}.
]
Quantization then yields
[
E_n^{(\chi)}=\frac{(n\hbar-\chi\kappa_0)^2}{2I},
\qquad n\in\mathbb Z.
]
This was the first route in the session that could simultaneously address three handoff requirements at once: a same-charge internal sector, a nonorbital first-order internal law, and a possible (2\pi/4\pi) identity distinction if the Berry flux happens to land on a half-integer. The handoff explicitly identified that last point as one of the things worth testing rather than assuming away.   

### 16.5 Why this became the last live corridor

At this stage of the session, the mixed-sector Berry rotor became the only route still plausibly compatible with the handoff’s minimum charged-lepton package. The scalar-bundle rescue of (g=2) had already failed, the purely static side-bit route had already failed, and the axisymmetric wall route could at best give a transverse rotor or an Ising-like (m=0) sign, but not a protected same-charge doublet. By contrast, the mixed-sector Berry construction reused only structures already present in the 4D program: the first-order matter action, the mixed gauge components, and the projection/open-system logic. That is why the session treated it as the last serious corridor rather than as another speculative branch.   

## 17. Why the Berry route still does not close the lepton target

### 17.1 The transverse branch is not yet pinned

The first missing ingredient is static pinning. The conservative 2PN throat already fixes positive odd dipole residues and positive odd poles, so the centered odd branch remains a stable normal mode in the current freeze. The even Family-1 wall/source law is also axisymmetric and support-dominated. So while the wall sector can soften the transverse radius (s=b_x^2+b_y^2), the current files do **not** yet derive a physical branch with
[
s_*>0.
]
If the transverse branch ever does pin, the axisymmetric wall package says it pins first to a **ring** of minima,
[
b_x^2+b_y^2=s_*,
]
not to a discrete two-state pair. That is already enough to show that the Berry route is not yet closed by the frozen 2PN data alone. 

### 17.2 There is no file-derived static (P_{22}) splitter

The second missing ingredient is angular splitting inside the transverse plane. The 2PN even operator already contains the full real (P_0\oplus P_2) support bundle, including the (P_{22c},P_{22s}) channels, but the solved static source vector is still purely axisymmetric:
[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right).
]
The Family-1 support/source law is likewise built only from (\mu=\cos\theta) through low-order flare polynomials, and even the strict PDE-side repair of the wall law remains axisymmetric. So the current freeze contains (P_{22}) as a **support capacity**, but not as a driven anisotropy source. There is therefore no file-derived static term of the form
[
-V_2\cos2(\phi-\phi_0)
]
that could split the transverse rotor into a genuine doublet. This is why the static side-bit rescue failed in the session. 

### 17.3 The current internal direction (w) is noncompact, not a compact spinor bundle

The third missing ingredient is topological/compactness structure. The frozen EM localization program uses a noncompact Gaussian direction
[
Z(w)=e^{-w^2/\lambda^2},
]
with Hermite/KK tower
[
f_n(w)=H_n(w/\lambda),
\qquad
m_n^2=\frac{2n}{\lambda^2}.
]
That is a line with localized modes, not a compact internal circle or double-cover bundle. The only explicit topology discriminator already frozen in the earlier bridge/full 1PN program is the added-mass split between a (w)-uniform throat and a compact (4)-ball,
[
\kappa_{\rm add}^{\rm throat}=\frac12,
\qquad
\kappa_{\rm add}^{(B^4)}=\frac13.
]
That is a real geometry/topology result, but it is not a spinor-like (\mathbb Z_2) or (4\pi)-return structure. So the present files do **not** already contain a compactness law that would force the Berry flux to a half-integer.   

### 17.4 The half-integer Berry-flux condition remains unforced

In the reduced rotor model, a protected low-energy doublet appears only when
[
\alpha\equiv \frac{\kappa_0}{\hbar}=\frac{\nu_0 s_*}{2\hbar}\in \mathbb Z+\frac12.
]
At that point the ground-state pair is degenerate and a (2\pi) rotation contributes a phase (-1) while (4\pi) returns (+1). But in the current frozen files, (\nu_0) is a continuous mixed-sector overlap coefficient and (s_*) is, at best, a yet-to-be-derived pinned radius. There is no built-in quantization law that forces (\alpha) to a half-integer rather than leaving it as a tunable continuous throat parameter. So the session’s most important negative result at this stage was precise: the ontology supports, at best, a same-charge Berry rotor, but **not yet** a protected charged-lepton doublet.   

### 17.5 Bottom line of the Berry audit

The Berry route therefore remains alive but incomplete. It solves the right kind of problem—it gives a same-charge, first-order, nonorbital internal structure—but it still lacks two closing ingredients: a way to pin the transverse branch at (s_*>0), and a way to discretize or quantize the resulting internal phase so the handoff’s (2\pi/4\pi) target is met for structural rather than accidental reasons. That is why the next exact falsifier became so narrow: either find a real (m=2) anisotropy source or a compact/topological phase law, or accept that this branch stalls at “same-charge Berry rotor without protected doublet.”  

## 18. Session endpoint: what is alive, what is dead, what is missing

### 18.1 What is still alive

By the end of the session, three things were clearly still alive inside the current ontology. First, the model can support a **same-charge nonorbital internal odd coordinate** once the mixed (a!-!w) sector is kept. Second, the mixed-sector route can support a **nonzero Berry/gyroscopic coefficient**, because the parent matter action is first order in time and the first odd (w)-mode survives in derivative overlaps even when it is hidden from a direct brane source. Third, if a transverse radius (s_*) is ever pinned, the resulting low-energy dynamics are those of a genuine **internal rotor** rather than a mere relabeling of orbital circulation. This is already much more than the scalar-monopole picture the project started from, and it is exactly the kind of “little extra” structure the handoff encouraged us to follow if it appeared.    

### 18.2 What is dead or stalled inside the current freeze

Several routes are now cleanly ruled out or stalled. Ordinary circulation/orbital relabelings fail the handoff’s intrinsic-spin bar outright. The scalar-bundle rescue of the magnetic sector fails to produce the missing Dirac piece in the current closure. The static side/orientation route fails because a pure (+w/-w) side bit is a scalar under the same in-plane (O(2)) symmetry and the current files do not provide a driven (P_{22}) anisotropy source. The axisymmetric Family-1 wall route can, at best, produce a continuous rotor for the transverse branch or an Ising-like (\pm) sign for the axisymmetric odd branch, but not a Pauli-like internal doublet. And the present Gaussian noncompact localization scheme gives no built-in half-integer/topological quantization of Berry flux. These are all genuine no-gos or stalls inside the current frozen ontology, not merely missing calculations.    

### 18.3 What the session isolated as the missing ingredients

The session narrowed the missing physics to a very small set. To keep the minimal charged-lepton route alive, the model now needs **either** a genuine non-axisymmetric (m=2) source that can split the transverse rotor, **or** a compact/topological phase law that quantizes the Berry flux, **or** some new first-order internal structure beyond the present conservative closure that effectively does one of those jobs. All three options are much sharper than the session’s starting point, because the earlier branches have already been pushed to their failure conditions. That is why the session’s endpoint should be viewed as progress rather than impasse: the model is no longer “missing everything.” It is missing something very specific.   

## 19. Best restart point for a fresh session

### 19.1 Start from the oriented-particle package, not from raw spin claims

A fresh session should **not** reopen the whole spin audit from scratch. The clean restart point is the oriented-particle package that the session built after the first no-gos: a 3D-facing mouth attached to a 4D corridor, with an explicit side/orientation question, a mixed-twist internal sector, projection/open-system diagnostics, and a clear distinction between what is already frozen and what is only proposed. That particle package is the right base layer because it is the first place where the hidden bulk channels, back-pressure, and internal mixed-sector structure are all present in one picture. Only after that layer is fixed does it make sense to ask again whether any surviving internal coordinate can meet the handoff’s charged-lepton target.   

### 19.2 Treat the mixed-sector Berry rotor as the only live corridor

The next session should treat the mixed-sector Berry rotor as the **only** live corridor unless new file-based evidence appears. The correct question is no longer “can the model somehow make spin?” It is: can the current oriented particle, with the mixed (A_w/F_{\mu w}/J^w) sector kept alive, produce either a real transverse pinning (s_*>0) or a genuine quantization/discretization law for the internal phase? That framing prevents a new session from wandering back into already-killed routes such as static scalar-bundle rescues or purely axisymmetric side-bit splitting.   

### 19.3 Frame the next work as identifying the missing discretizer/quantizer

The cleanest “fresh start” prompt is therefore not “derive spin,” but rather: identify the missing **discretizer** or **quantizer** of the mixed-sector rotor. Concretely, the next session should ask one of two sharply framed questions. First: is there any honest non-axisymmetric source, beyond the current axisymmetric Family-1 freeze, that can activate the passive (P_{22}) slots and generate a real (\cos2\phi) splitter? Second: is there any compact or topological internal phase law—perhaps associated with a modified localization structure, an internal bundle, or a throat-specific cover-space variable—that can force
[
\frac{\nu_0 s_*}{2\hbar}\in \mathbb Z+\frac12?
]
If the answer to both is no, the minimal charged-lepton route should be declared stalled. If either answer is yes, then the next phase finally has a concrete nonorbital doublet target. That is the most productive restart point because it is both honest about the current no-gos and maximally sensitive to the one piece of new physics the model still seems to need.   
