# From 4D Action to Atomic Structure and the Lepton g-Factor: Reduced-Sector Derivations in the 4D Superfluid Model

# Part I — Foundations

*Draft write-up / handoff artifact for the 4D superfluid program.*

This part establishes the shared dictionary and claim boundaries used by the later atomic and lepton sections. Its job is not to prove the hydrogen or lepton results by itself. Its job is to make clear what the underlying ontology is, what is exact in the current stack, what is only available through controlled reduction, and what remains outside present closure.

The most important discipline carried through this document is that the project does **not** treat every equation on equal footing. Some statements follow directly from the declared 4D action. Some arise only after explicit brane reduction. Some are fixed only within a stated internal-response protocol. The later atomic and lepton derivations only make sense if that hierarchy stays visible from the beginning.

---

## 4. Ontology of the 4D superfluid defect model

### 4.1 Brane, bulk, throats, and localized defects

The current program is built on a fully action-based model in which the fundamental arena is a 
\((4+1)\)-dimensional spacetime with coordinates
\[
X^M=(t,x,y,z,w),
\qquad
\mathbf X=(x,y,z,w)\in\mathbb R^4,
\qquad
\mathbf x=(x,y,z)\in\mathbb R^3.
\]
The extra coordinate \(w\) is not decorative. It is the direction along which the theory distinguishes the full bulk description from what a three-dimensional observer on the brane actually measures. The familiar three-dimensional world is therefore not the fundamental stage of the model but a reduced or projected sector of a higher-dimensional fluid-like substrate.

Within that ontology, localized matter-like objects are modeled as defect structures with a throat interpretation. The brane-facing side of such an object appears as a localized mouth or puncture near \(w=0\), while the interior extends into the bulk as a tube- or cavity-like structure. The older geometric language of the project therefore survives the move to the full 4D action: a defect is still a **brane–bulk throat**, but it is now embedded in a precise bulk field theory rather than being discussed only through effective-side bookkeeping.

The geometry sector is intentionally kept explicit. Two effective geometric degrees of freedom,
\[
a(t),\qquad L(t),
\]
represent the throat radius and axial extent. They are not mere labels; they belong to the dynamical closure of the model and can exchange energy with the matter and gauge sectors. In the parent paper this appears through a geometry action and damped evolution laws for \(a\) and \(L\), with quasi-static behavior recovered only as a controlled limit rather than as the default interpretation.

A crucial consequence of this ontology is that the brane observer does not automatically see the full defect. What the brane observer sees is always a reduced image of bulk structure. That is why the same object can look approximately monopolar on the brane while still supporting nontrivial interior support, wake, and response structure once the theory is pushed beyond the lowest-order scalar sector.

### 4.2 Matter sector as gauged 4D GNLS

The matter/order-parameter sector is a complex bulk field
\[
\psi(\mathbf X,t),\qquad \rho(\mathbf X,t)=|\psi|^2,
\]
with a gauged four-spatial-dimensional nonlinear Schrödinger structure. The exact matter action is written in terms of gauge-covariant derivatives,
\[
D_t\psi=\partial_t\psi+\frac{i q_*}{\hbar}A_0\psi,
\qquad
D_i\psi=\partial_i\psi-\frac{i q_*}{\hbar}A_i\psi,
\]
where the signed microscopic branch coupling is
\[
q_*\equiv \eta_Q e_*,\qquad \eta_Q=\pm 1,\qquad e_*>0.
\]
The matter Lagrangian density is
\[
\mathcal L_\psi
=
\frac{i\hbar}{2}(\psi^*D_t\psi-\psi D_t\psi^*)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
-V_{\rm conf}(\mathbf X;a,L)|\psi|^2
-U(\rho).
\]

The self-interaction is not arbitrary. It is frozen to the stiff-polytropic choice already selected earlier in the program:
\[
P(\rho)=K\rho^5,
\qquad
U(\rho)=\frac K4\rho^5,
\qquad
h(\rho)=\frac{dU}{d\rho}=\frac{5K}{4}\rho^4.
\]
This keeps the equation-of-state exponent \(n=5\) visible as a carry-forward constraint rather than a tunable parameter. It also means the matter sector has a definite hydrodynamic interpretation: density, current, sound speed, and pressure are inherited from one fixed GNLS medium rather than inserted ad hoc for each new application.

Varying the action yields the exact bulk matter equation
\[
i\hbar D_t\psi
=
\left[-\frac{\hbar^2}{2m}D_iD_i+V_{\rm conf}(\mathbf X;a,L)+h(\rho)\right]\psi,
\]
together with the exact conserved number current
\[
j^i=\frac{\hbar}{m}\,\mathrm{Im}(\psi^*D_i\psi),
\qquad
\partial_t\rho+\partial_i j^i=0.
\]
When \(\rho>0\), one may rewrite the current as \(j^i=\rho v^i\), which is the starting point of the Madelung/hydrodynamic picture used repeatedly throughout the program.

For the purposes of the write-up, the key point is not merely that a nonlinear Schrödinger field exists. The key point is that the project’s later atomic and lepton derivations start from a matter sector that is already exact, already gauged, already endowed with a fixed EOS, and already coupled to a dynamical throat geometry. They do not begin from a point-particle Bohr model.

### 4.3 Electromagnetic sector and localized Maxwell reduction

The electromagnetic sector is supplied by a localized Maxwell theory in \(4+1\) dimensions. The gauge field is
\[
A_M=(A_0,A_i),\qquad F_{MN}=\partial_MA_N-\partial_NA_M,
\]
and its action includes a transverse localization profile \(Z(w)\), gauge fixing, and external sources. At the level of the unified theory, the electromagnetic Lagrangian density is
\[
\mathcal L_{\rm EM}
=
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
-A_MJ_{\rm ext}^M.
\]
Because the matter sector is already minimally coupled, variation of the matter action automatically generates the dynamical matter current. The explicit source term in the electromagnetic action is therefore reserved for external/background sources rather than duplicating the matter current by hand.

This sector matters for two reasons. First, it gives an exact bulk gauge theory rather than a purely kinematic “EM-from-flow” dictionary. Second, it produces ordinary \(3+1\) Maxwell theory on the brane as a **controlled reduction** of the localized \(4+1\) action. Under the explicit zero-mode assumptions used in the Maxwell limit, integration over the transverse direction yields an effective brane coupling
\[
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}},
\qquad
Z_{\rm int}=\int_{-\infty}^{\infty}Z(w)\,dw.
\]
The same reduction implies the thickness-controlled observable charge strength
\[
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}}.
\]

That reduced Maxwell law is the entry point for the later atomic program. It is what allows a Coulomb sector to appear on the brane while still keeping the microscopic charge bookkeeping inside the higher-dimensional ontology.

The electromagnetic paper is equally explicit about what the reduction does **not** do. It suppresses the mixed core channels
\[
(A_w,\,J^w,\,F_{\mu w})
\]
only in a controlled far-field brane limit. It does not erase them from the microscopic theory. That distinction becomes essential later, because the atomic finite-size program and the lepton spin/anomaly program both eventually lean on precisely the channels that the zero-mode Maxwell limit leaves dormant.

---

## 5. Reduction principles needed for the atomic and lepton programs

### 5.1 Reduction versus projection

The parent 4D papers make a distinction that the later derivations must preserve: **projection** and **reduction** are not the same operation.

Projection defines what a brane observer measures. Given a fixed normalized weight \(W(w)\), a projected brane observable is a weighted average,
\[
\langle f\rangle_W(\mathbf x,t)=\int_{-\infty}^{\infty}W(w)f(\mathbf x,w,t)\,dw.
\]
This is an operational definition. It is exact once \(W\) is declared.

Reduction, by contrast, is a controlled elimination of the extra coordinate under an explicit ansatz or scaling regime. It is what produces effective brane equations such as the Poisson hook or the zero-mode Maxwell law. Reduction changes the dynamical description. Projection does not.

This distinction is central for the later atomic and lepton sections. The hydrogenic Coulomb problem is a **reduced** brane-effective sector, not merely a projected snapshot of the full bulk system. Conversely, questions about mixed core structure, throat response, or lepton-like internal modes cannot be settled by projection alone, because projection hides the very channels that may matter microscopically.

A second reason this distinction matters is methodological. The project’s strongest results come from preserving the difference between what is exact in the parent theory and what only appears after controlled brane-effective reduction. If that line is blurred, later claims begin to look stronger than they really are.

### 5.2 Zero-mode Maxwell limit and its controlled domain

The effective brane Maxwell law is obtained only after a specific reduction regime is declared. In practice, the far-field brane limit suppresses the transverse/mixed sector through assumptions of the form
\[
A_w\approx 0,
\qquad
\partial_w A_\mu\approx 0,
\qquad
J^w\approx 0,
\qquad
F_{\mu w}\approx 0,
\]
with brane-localized sources and zero-mode dominance. In that regime the localized \(4+1\) Maxwell system collapses to ordinary \(3+1\) Maxwell theory with renormalized coupling.

This is exactly the regime needed for the first-pass atomic reduction. It legitimizes using an effective Coulomb interaction on the brane and makes the hydrogenic sector mathematically well-posed within the existing stack.

But the domain of validity is just as important as the result itself. The zero-mode reduction is a long-distance, low-frequency, far-field brane law. It is not the final word on what happens near the defect core, at finite throat size, or when transverse excitations become relevant. The EM paper is explicit that beyond-zero-mode structure appears as higher transverse modes and Yukawa-suppressed corrections, with the mixed sector remaining part of the microscopic ontology. That is why the later write-up must treat the zero-mode Coulomb law as a clean leading sector, not as the whole story.

### 5.3 Charge sign as topological branch label

One of the most important corrections in the newer stack is the charge-ontology firewall. Electric charge sign is no longer identified with circulation, throat radius, or breathing variables. It is carried by a fixed topological puncture orientation,
\[
\eta_Q=\pm 1,
\qquad
q_*=\eta_Q e_*.
\]
This matters because the older effective-side language could tempt later sections into treating circulation as if it directly defined electric sign. The corrected ontology forbids that move.

The same firewall also cleans up the notation inherited from the gravity-side papers. The historical scalar coefficient often written as bare \(q=1\) in the Newtonian/1PN gravity ledger is **not** electric charge. In the corrected notation it is written \(\kappa_\rho=1\). This is more than a cosmetic fix. Without it, later lepton or atomic arguments would risk mixing gravitational bookkeeping coefficients with electromagnetic branch data.

The positive consequence of the correction is that the model now has a cleaner division of labor. Electric sign is topological. Electric strength is determined by localization. Circulation survives, but it belongs to the magnetic/vortical sector rather than to the electric-charge dictionary. That separation is one of the enabling conditions for the later reduced derivations of magnetic moment and gyromagnetic ratio.

### 5.4 Observable charge magnitude as thickness-controlled coupling

Once charge sign is separated from circulation, the next key principle is that the **observable** brane charge is thickness-controlled. The microscopic coupling \(e_*\) belongs to the branch itself, but the charge measured on the brane is weakened by the localization thickness:
\[
e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}},
\qquad
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}}.
\]
For a Gaussian localization profile,
\[
Z(w)=e^{-w^2/\lambda^2},
\qquad
Z_{\rm int}=\lambda\sqrt\pi.
\]
So thicker localization means a larger \(Z_{\rm int}\) and therefore a weaker effective brane coupling.

This principle is one of the most important bridges into the later reduced-sector work. In the hydrogenic sector it means the orbital scale can be expressed directly in terms of the localization thickness once the reduced Coulomb coefficient is inserted. In the lepton anomaly program it provides the first clean structural reason that support redistribution and boundary-layer effects can renormalize the observable electromagnetic response without changing the topological branch identity itself.

For the write-up, the right way to state the principle is simple: topology fixes the sign; localization thickness controls the measured strength.

---

## 6. Scope of the present derivations

The current document is not built on one undifferentiated claim type. It stands on a hierarchy of exact statements, controlled reductions, protocol closures, and reduced-sector consequences. The safest way to continue the write-up is to keep that hierarchy explicit.

| Status label | Meaning in this document | Typical examples |
|---|---|---|
| **Exact** | Follows directly from the declared parent action, projection map, or exact algebra | Bulk GNLS equation, exact projected continuity with leakage, exact longitudinal identity |
| **Controlled reduction** | Requires an explicit regime assumption or reduced ansatz | Poisson hook, zero-mode Maxwell reduction, worldtube/small-body limit |
| **Protocol closure** | Fixed only within a declared internal-response or support hierarchy | Adiabatic breathing closure at 1PN, low-frequency throat-response packaging at 2PN |
| **Reduced-sector consequence** | A later derivation built from the previous layers without claiming a full solved bulk-defect theorem | Hydrogenic Bohr-scale reduction, finite-size \(P_{22}\) response, reduced \(g=2\) bridge |

### 6.1 What is derived from the current action stack

The parent 4D action stack already supplies several pieces that are solid enough to support later derivations.

First, the project now has a fully action-based bulk parent theory with three explicit sectors: the gauged 4D GNLS matter sector, the localized \(4+1\) Maxwell sector, and the geometry sector with effective throat variables \((a,L)\). That alone is a major shift from the earlier purely effective-side language.

Second, several exact structures carry forward independently of later closure choices: the exact Euler–Lagrange equations, exact conserved current and bulk continuity, exact projected continuity with a leakage source \(S_{\rm leak}\), and the exact longitudinal identity that underlies the gravitational Poisson hook.

Third, the program already has two high-value controlled reductions that later parts repeatedly reuse. One is the quasi-static Poisson hook that produces a brane-facing inverse-square scalar sector. The other is the zero-mode Maxwell reduction that produces a brane Maxwell/Coulomb sector with thickness-controlled coupling.

Fourth, the conservative post-Newtonian gravity side is already much stronger than a collection of heuristics. The full conservative 1PN assembly reproduces the standard EIH two-body ledger exactly **within a declared closure hierarchy**, while the 2PN assembly reproduces the standard generic-frame ADM Hamiltonian through \(2\)PN, again within a declared closure hierarchy. The 2PN papers also show that once the theory is pushed harder, the defect can no longer be treated as a pure scalar monopole: a carried-forward dipole wake is joined by a new even \(P_0\oplus P_2\) mouth/support sector and a separate geometry-closure channel.

Those are the backbone results that later sections are allowed to reuse without apology. In other words, the later atomic and lepton programs are not being built on a vacuum. They are being built on an exact parent action plus a frozen set of controlled brane-effective and post-Newtonian ledgers.

### 6.2 What remains outside present closure

At the same time, the current stack still has hard boundaries, and the write-up should say so plainly.

The project does **not** yet claim a fully solved moving-throat PDE theorem in the 4D bulk. The particle limit remains a controlled reduction rather than a numerically completed defect solution. The open-system objects such as leakage, boundary forces, and related exchange channels exist exactly in the parent description, but the conservative 1PN and 2PN derivations do not amount to a finished dynamical solution of the full open-system bulk problem.

Several sectors also remain explicitly unresolved in the present papers: spin couplings, dissipative leakage and radiation reaction, strong-field completion, higher-PN sectors beyond the solved conservative hierarchy, and the completed dynamic derivation of all appendix-side pole or compliance data. Even where the 2PN notes reveal throat-response structure, not every response datum is claimed as an assumption-free bulk theorem.

That limitation directly constrains how the later sections should be read. The hydrogenic, finite-size, \(P_{22}\), and lepton \(g\)-factor results developed in the companion notes are best described as **reduced-sector derivations inside the current hierarchy**, not as final first-principles theorems of the fully solved 4D defect. They are legitimate because they respect the existing action, reduction rules, and ontology firewall. But they must still be labeled honestly as closure-level results where appropriate.

For the full write-up, the safest summary is this: the foundations are strong enough to support serious reduced-sector atomic and lepton derivations, but the complete moving-throat and full spin-coupled theory is still open.

---

# Part II — Hydrogen from the Existing 4D Action

*Draft write-up / handoff artifact for the 4D superfluid program.*

This part isolates the first atomic target that can already be attacked with the current 4D stack: the hydrogenic bound-state scale. The aim is deliberately narrow. It is **not** to claim a finished first-principles theorem of atomic structure from the full moving-throat PDE. It is to show that the present action, together with its already-declared reduction rules, supports a controlled hydrogenic sector in which the familiar Bohr scale appears as an **energy minimum** rather than as a hand-inserted Bohr circulation rule.

That distinction matters. In the old language, one could always recover the Bohr radius by rephrasing the Bohr/de Broglie condition in fluid vocabulary. The present program is stronger only if the scale emerges from the reduced 4D action itself. The derivation below shows that such a path already exists, first in a one-electron / fixed-source sector and then in a genuine two-defect sector where the reduced mass appears automatically.

The claims in this part are therefore all of the same kind: they are **reduced-sector derivations built from the exact parent action and controlled zero-mode reduction**, not yet full moving-throat theorems.

---

## 7. The one-electron hydrogenic reduction

The hydrogenic program starts from the exact structures already established in Part I:

- the gauged 4D GNLS matter sector,
- the localized Maxwell sector with thickness-controlled brane coupling,
- and the distinction between exact bulk dynamics and controlled brane reduction.

The narrow question of this section is simple:

> If one places a light negative branch in the Coulomb field of a heavy positive branch and stays inside the controlled zero-mode Maxwell regime, what effective three-dimensional bound-state problem does the 4D action produce?

### 7.1 Separation of the lowest transverse mode

In the first-pass hydrogenic sector, take the positive branch as a static heavy source located at the brane origin. Work in the same zero-mode/far-field Maxwell regime already used in the EM reduction:
\[
A_w\approx 0,
\qquad
\partial_w A_\mu\approx 0,
\qquad
J^w\approx 0,
\qquad
F_{\mu w}\approx 0.
\]
The matter field is then factorized into a lowest transverse profile and a brane wavefunction,
\[
\psi(\mathbf x,w,t)=\chi_0(w)\,\phi(\mathbf x,t),
\qquad
\int_{-\infty}^{\infty}|\chi_0(w)|^2\,dw=1.
\]

This is the natural first reduction because the atomic target is not trying to resolve transverse excited matter modes. The lowest mode is enough to ask whether a hydrogenic scale exists at all.

After substituting this ansatz into the exact matter action and integrating over the transverse coordinate, two reduced quantities appear immediately:
\[
E_\perp
\equiv
\int dw\,\chi_0^*(w)
\left[
-\frac{\hbar^2}{2m}\partial_w^2+V_{\rm conf}^{(w)}(w;a,L)
\right]\chi_0(w),
\]
which is the transverse confinement offset, and
\[
\Gamma_{10}\equiv \int dw\,|\chi_0(w)|^{10},
\]
which controls the inherited stiff-EOS nonlinear correction after reduction.

The structural point is important. The matter reduction does **not** delete the nonlinear GNLS sector. It carries it forward into the brane theory through a definite overlap factor. So even the first hydrogenic reduction already knows that the parent medium is not ordinary linear Schrödinger mechanics.

### 7.2 Effective 3D Coulomb sector on the brane

The gauge side supplies the other half of the hydrogen problem. In the static zero-mode Maxwell sector, the positive branch generates a brane scalar potential with the usual Coulomb zero mode plus a KK/Yukawa tower. Writing the reduced attractive coefficient as \(g_C>0\), the light opposite-charge branch feels the effective hydrogenic potential
\[
V_{\rm H}(r)
= -\frac{g_C}{r}
-\frac{g_C}{2r}e^{-2r/\lambda}
+\mathcal O\!\left(\frac{e^{-2\sqrt2\,r/\lambda}}{r}\right).
\]
The leading term is Coulombic. The next term is the first even transverse correction.

At this stage it is best to keep \(g_C\) symbolic. Under standard brane-Maxwell matching,
\[
g_C = \frac{e_{\rm eff}^2}{4\pi\epsilon_0}
= \frac{e_*^2}{4\pi\epsilon_0 Z_{\rm int}},
\qquad
Z_{\rm int}=\int Z(w)\,dw.
\]
So the observable Coulomb strength is already thickness-controlled. The microscopic branch parameter \(e_*\) is not the same thing as the charge measured on the brane.

For Gaussian localization,
\[
Z(w)=e^{-w^2/\lambda^2},
\qquad
Z_{\rm int}=\lambda\sqrt\pi,
\qquad
e_{\rm eff}=\frac{e_*}{\sqrt{\lambda\sqrt\pi}}.
\]
This relation is what later turns the hydrogenic scale into a direct diagnostic of localization thickness.

### 7.3 KK/Yukawa and GNLS correction channels

Once the reduction is written honestly, two correction channels are already visible before any precision fitting begins.

The first is the KK/Yukawa correction from the localized Maxwell tower,
\[
\delta V_Y(r)=-\frac{g_C}{2r}e^{-2r/\lambda}+\cdots,
\]
which strengthens attraction relative to pure Coulomb at short range.

The second is the inherited GNLS self-pressure term. After transverse reduction, the one-electron brane action contains
\[
-\frac{K\Gamma_{10}}{4}|\phi|^{10},
\]
so the reduced equation of motion becomes
\[
i\hbar\partial_t\phi
=
\left[
-\frac{\hbar^2}{2m}\nabla^2
+E_\perp
+V_{\rm H}(r)
+\frac{5K\Gamma_{10}}{4}|\phi|^8
\right]\phi.
\]
This is the first genuine hydrogenic field equation extracted from the current stack. It is not yet full hydrogen as a finished theorem, but it is already a nontrivial reduced consequence of the parent 4D action.

Two transverse structures should be kept distinct from the outset:

- \(Z(w)\), the EM localization profile that controls \(e_{\rm eff}\), and
- \(\chi_0(w)\), the matter transverse profile that controls \(E_\perp\) and \(\Gamma_{10}\).

A future geometry closure may tie both to the same throat data, but that relation is not yet derived in the present stack. Part II therefore keeps them separate on purpose.

---

## 8. The hydrogenic energy functional

With the reduced brane equation in hand, the next question is whether it really supports a Bohr-type minimum. The cleanest way to test that is variationally.

### 8.1 Variational \(1s\) ansatz

Normalize the reduced brane field by
\[
\int d^3x\,|\phi|^2=1.
\]
For a first-pass hydrogenic probe, use the standard normalized \(1s\)-type trial state
\[
\phi_a(r)=\frac{1}{\sqrt{\pi a^3}}e^{-r/a},
\qquad a>0.
\]
This ansatz does not assume Bohr quantization. It simply asks whether the reduced action prefers a particular bound-state size.

The needed expectation values are elementary:
\[
\langle T\rangle = \frac{\hbar^2}{2ma^2},
\qquad
\left\langle \frac1r \right\rangle = \frac1a,
\qquad
\left\langle \frac{e^{-2r/\lambda}}{r} \right\rangle
= \frac{\lambda^2}{a(a+\lambda)^2},
\]
and for the inherited nonlinearity,
\[
\int d^3x\,|\phi_a|^{10} = \frac{1}{125\pi^4 a^{12}}.
\]

### 8.2 Reduced energy as a function of orbital scale

The resulting one-parameter energy is
\[
\boxed{
E(a)
=
E_\perp
+\frac{\hbar^2}{2ma^2}
-\frac{g_C}{a}
-\frac{g_C}{2a(1+a/\lambda)^2}
+\frac{K\Gamma_{10}}{500\pi^4 a^{12}}
+\cdots
}
\]
where the omitted terms are higher KK/Yukawa contributions beyond the leading even mode.

This formula is the first real hydrogenic checkpoint in the program. It says the current 4D action already produces a brane effective energy with exactly the structure needed for a bound-state minimum:

- a kinetic term that prefers delocalization,
- an attractive Coulomb sector that prefers localization,
- a short-range KK correction that increases attraction,
- and a stiff nonlinear term that resists overcompression.

That is already enough to ask whether the Bohr scale appears dynamically.

### 8.3 Clean decoupling limit

The clean hydrogenic limit is the regime in which both correction channels are small at the bound-state scale:
\[
\lambda\ll a,
\qquad
\frac{K\Gamma_{10}}{a^{12}}\ll \frac{g_C}{a}.
\]
In that regime the reduced energy collapses to
\[
E(a)\approx E_\perp+\frac{\hbar^2}{2ma^2}-\frac{g_C}{a}.
\]
This is the precise sense in which the present program reaches a hydrogenic Coulomb problem. It is not claimed as an exact microscopic theorem at all scales. It is claimed as a controlled decoupling sector of the parent action.

The full reduced functional also shows immediately where the model can fail. If the KK correction or the GNLS stiffness is not negligible at hydrogenic scale, the minimum will move. Those shifts are not nuisances; they are falsifiers.

---

## 9. Bohr radius without imposed Bohr quantization

The central result of Part II is that the Bohr scale appears as an energy minimum of the reduced action. No separate Bohr/de Broglie circulation rule is needed.

### 9.1 Action-based minimization

Minimizing the clean-limit energy,
\[
E(a)\approx E_\perp+\frac{\hbar^2}{2ma^2}-\frac{g_C}{a},
\]
with respect to \(a\) gives
\[
\frac{dE}{da}=-\frac{\hbar^2}{ma^3}+\frac{g_C}{a^2}=0.
\]
So the unique positive minimizer is
\[
\boxed{a_* = \frac{\hbar^2}{m g_C}}.
\]
This is the exact point at which the present program stops being merely a rephrased Bohr model. The quantized length scale appears from a variational competition inside the reduced field theory.

### 9.2 Recovery of the Bohr-scale radius

Under standard brane-Maxwell matching,
\[
g_C=\frac{e_{\rm eff}^2}{4\pi\epsilon_0},
\]
so the minimum becomes
\[
\boxed{
a_* = \frac{4\pi\epsilon_0\hbar^2}{m e_{\rm eff}^2}.
}
\]
Once the light reduced mass parameter is identified with the electron mass parameter in the one-body sector, this is exactly the standard Bohr-radius expression in its fixed-source form.

The same minimum also gives the clean hydrogenic binding scale,
\[
E_* - E_\perp = -\frac{m g_C^2}{2\hbar^2}
= -\frac{m e_{\rm eff}^4}{2(4\pi\epsilon_0)^2\hbar^2}.
\]
So the reduced action is already pointing to both the familiar length scale and the familiar binding scale.

### 9.3 Thickness law for effective charge and orbital size

The action-based hydrogen minimum becomes especially interesting after inserting the thickness-controlled brane charge law,
\[
e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}}.
\]
Then
\[
\boxed{
a_* = \frac{4\pi\epsilon_0\hbar^2 Z_{\rm int}}{m e_*^2}.
}
\]
For Gaussian localization,
\[
\boxed{
a_* = \frac{4\pi^{3/2}\epsilon_0\hbar^2}{m e_*^2}\,\lambda.
}
\]
This is the clearest “charge-from-thickness” consequence available in the present hydrogen sector. Thicker localization weakens the observable brane coupling and therefore inflates the orbital scale linearly.

That is a real structural output of the current theory, not a verbal gloss. It says hydrogen is not just a test of Coulomb attraction. It is already a test of how localization in the transverse direction renormalizes observable charge on the brane.

The first correction shifts are also straightforward to read off around the clean minimum \(a_0=\hbar^2/(mg_C)\):
\[
\delta E_Y(a)= -\frac{g_C}{2a(1+a/\lambda)^2},
\qquad
\delta E_{\rm NL}(a)=\frac{K\Gamma_{10}}{500\pi^4 a^{12}}.
\]
The KK sector pulls the minimum inward; the stiff GNLS term pushes it outward. A future precision hydrogen program would have to show that both are small in the physical regime.

---

## 10. Two-body upgrade

The fixed-source hydrogenic result is already meaningful, but it is not yet the cleanest formulation. The next obvious upgrade is to treat the atom as a true two-defect problem rather than a light branch in an externally imposed field.

### 10.1 Genuine two-defect reduction

Let the negative and positive defect centers on the brane be \(\mathbf x_e(t)\) and \(\mathbf x_p(t)\). In the coherent-defect / small-worldtube regime, the reduced conservative action is taken to be
\[
S_{ep}^{\rm red}
=
\int dt\,
\Bigg[
\frac12 m_e \dot{\mathbf x}_e^2
+\frac12 m_p \dot{\mathbf x}_p^2
-E_e^{\rm int}
-E_p^{\rm int}
-V_{ep}(r)
+\delta L_{\rm fs}
\Bigg],
\]
with
\[
\mathbf r\equiv \mathbf x_e-\mathbf x_p,
\qquad
r=|\mathbf r|.
\]
The inter-defect potential takes the same reduced Maxwell form as before,
\[
V_{ep}(r)
=
-\frac{g_C}{r}
-\frac{g_C}{2r}e^{-2r/\lambda}
+\mathcal O\!\left(\frac{e^{-2\sqrt2\,r/\lambda}}{r}\right)
+\delta V_{\rm fs}(r),
\]
where \(\delta V_{\rm fs}\) now explicitly collects finite-size and throat-response effects beyond the point-defect sector.

This is already conceptually cleaner than the one-body note because the internal GNLS and confinement structure is reclassified as defect self-energy rather than being left as a leading orbital mean-field term.

### 10.2 Emergent reduced mass

Introduce the center-of-mass and relative coordinates,
\[
M\equiv m_e+m_p,
\qquad
\mu\equiv \frac{m_e m_p}{m_e+m_p},
\]
\[
\mathbf R\equiv \frac{m_e\mathbf x_e+m_p\mathbf x_p}{M},
\qquad
\mathbf r\equiv \mathbf x_e-\mathbf x_p.
\]
Then the kinetic term splits exactly:
\[
\frac12 m_e \dot{\mathbf x}_e^2+\frac12 m_p \dot{\mathbf x}_p^2
=
\frac12 M\dot{\mathbf R}^2+\frac12 \mu\dot{\mathbf r}^2.
\]
So the reduced mass is not inserted by hand. It emerges automatically from the two-body reduction.

The relative Hamiltonian is therefore
\[
H_{\rm rel}
=
-\frac{\hbar^2}{2\mu}\nabla_r^2
+E_{\rm int}^{\rm tot}
-\frac{g_C}{r}
-\frac{g_C}{2r}e^{-2r/\lambda}
+\delta V_{\rm fs}(r),
\]
with
\[
E_{\rm int}^{\rm tot}=E_e^{\rm int}+E_p^{\rm int}.
\]

That is the correct hydrogenic central-force sector in the present reduced ontology.

### 10.3 Recovery of the fixed-source limit

Using the same normalized \(1s\)-type trial state for the relative coordinate,
\[
\phi_a(r)=\frac{1}{\sqrt{\pi a^3}}e^{-r/a},
\]
one finds the relative energy
\[
\boxed{
E_{\rm rel}(a)
=
E_{\rm int}^{\rm tot}
+\frac{\hbar^2}{2\mu a^2}
-\frac{g_C}{a}
-\frac{g_C}{2a(1+a/\lambda)^2}
+\langle \delta V_{\rm fs}\rangle_a.
}
\]
In the clean two-body hydrogenic limit,
\[
\lambda\ll a,
\qquad
|\langle \delta V_{\rm fs}\rangle_a|\ll \frac{g_C}{a},
\]
the minimum is
\[
\boxed{
a_* = \frac{\hbar^2}{\mu g_C}
= \frac{4\pi\epsilon_0\hbar^2}{\mu e_{\rm eff}^2}
= \frac{4\pi\epsilon_0\hbar^2 Z_{\rm int}}{\mu e_*^2}.
}
\]
For Gaussian localization,
\[
\boxed{
a_* = \frac{4\pi^{3/2}\epsilon_0\hbar^2}{\mu e_*^2}\,\lambda.
}
\]
So the thickness law survives intact, but the kinematic mass is now the correct reduced mass.

The fixed-source result is recovered smoothly in the heavy-source limit,
\[
m_p\gg m_e
\qquad\Longrightarrow\qquad
\mu\to m_e.
\]
This is an important consistency check. The one-body hydrogen note was not an artifact. It was the heavy-source limit of a cleaner pair reduction.

---

## What Part II establishes

Part II does establish several things already.

1. The current 4D action stack supports a controlled hydrogenic reduction.
2. In that reduction, the Bohr scale appears as a variational minimum of the reduced action.
3. The observable atomic scale is thickness-controlled through the brane charge renormalization \(e_{\rm eff}=e_*/\sqrt{Z_{\rm int}}\).
4. The two-body upgrade naturally produces the reduced mass \(\mu\).
5. The same reduced functional already isolates the first correction channels that can shift the hydrogenic radius: KK/Yukawa attraction and finite-size/throat-response structure.

## What Part II does not yet establish

Part II does **not** yet provide a full first-principles theorem of hydrogen from the complete moving-throat PDE.

It does not yet derive:

- the microscopic constants \(m\), \(m_e\), \(m_p\), or \(e_*\) from throat geometry alone,
- the full finite-size interaction \(\delta V_{\rm fs}(r)\),
- the complete relation between the EM profile \(Z(w)\) and the matter profile \(\chi_0(w)\),
- or the angular/orbital mode structure of excited hydrogenic states.

Those missing pieces matter, but they do not erase the central result of this part: there is already a mathematically clean route from the existing 4D action to the Bohr-scale hydrogenic minimum.

---

# Part III — Finite-Size Atomic Response and Core Regulation

*Draft write-up / handoff artifact for the 4D superfluid program.*

Part II ended with one unresolved object:
\[
\delta V_{\rm fs}(r).
\]
Once hydrogen is written as a genuine two-defect problem, the remaining atomic question is no longer vague. It becomes precise: what internal perturbation does one finite throat feel in the field of another, and how does that response modify the short-distance sector of the bound pair?

This part develops the first answer already supported by the current stack. It does four things.

First, it shows that the internally active atomic load is **not** the raw Coulomb well depth itself. In the frozen zero-mode Maxwell sector, the first shape-sensitive load begins at the **Hessian** of the partner potential across the finite defect. That gives a natural
\[
P_0\oplus P_2
\]
load multiplet consisting of a scalar trace and a quadrupolar tidal tensor.

Second, it uses the 1PN bridge response logic and the already-frozen 2PN support hierarchy to identify the first usable finite-size correction channels. At an interim bridge level, a scalar mouth-load embedding already yields a practical inverse-square correction. But the direct action-level perturbation map then sharpens the picture: pure Coulomb drives the quadrupole channel and leads to an inverse-sixth far-field response, while the scalar breathing channel enters only through the finite-thickness Yukawa tower.

Third, it shows that the same quadrupole tide that modifies hydrogen is also the first physical source for the dormant \(P_{22}\) mouth anisotropy already present in the 2PN throat-response bundle.

Fourth, it resolves the apparent short-distance divergence by replacing the point-Hessian approximation with a finite-throat tidal kernel. The regulator is not imposed by hand. It is supplied by the physical waist of the throat itself.

As in the earlier parts, every result below should be read at the same claim level: these are **reduced-sector consequences of the exact parent action plus controlled reduction and frozen support/response hierarchy**, not yet full theorems of the completely solved moving-throat PDE.

---

## 11. Atomic perturbation map from the reduced action

The first task is to identify what the partner defect actually loads.

Part II already produced a reduced two-body hydrogenic sector with an inter-defect potential of the form
\[
V_{ep}(r)
=
-\frac{g_C}{r}
-\frac{g_C}{2r}e^{-2r/\lambda}
+\cdots.
\]
But that formula by itself only gives the center-to-center interaction energy. It does not yet tell us what the field of one defect does to the **internal coordinates** of the other.

### 11.1 Why the internal load is not the raw Coulomb scalar

Take one defect, call it body \(A\), with brane center \(R_A\) and internal coordinates \(Q^a\). Let body \(B\) generate an external static potential energy
\[
V_B(\mathbf x)
\]
felt by the normalized internal density of \(A\). Write the internal coordinates of \(A\) relative to its center as
\[
\boldsymbol\xi=(\xi^1,\xi^2,\xi^3,\xi^w),
\qquad
\int d^3\xi\,dw\;\rho_A(\boldsymbol\xi;Q)=1.
\]
Then the external interaction energy of \(A\) is
\[
H_A^{\rm ext}[Q;V_B]
=
\int d^3\xi\,dw\;\rho_A(\boldsymbol\xi;Q)
V_B(R_A+\boldsymbol\xi).
\]

Now expand the partner potential around the center of \(A\):
\[
V_B(R_A+\boldsymbol\xi)
=
V_B(R_A)
+\xi^i\partial_iV_B(R_A)
+\frac12\xi^i\xi^j\partial_i\partial_jV_B(R_A)
+\cdots,
\]
where in the frozen zero-mode sector the relevant spatial derivatives are brane derivatives \(i,j\in\{x,y,z\}\).

Define the internal brane moments
\[
d_A^i(Q)
\equiv
\int d^3\xi\,dw\;\rho_A\,\xi^i,
\qquad
M_A^{ij}(Q)
\equiv
\int d^3\xi\,dw\;\rho_A\,\xi^i\xi^j.
\]
Then
\[
H_A^{\rm ext}
=
V_B(R_A)
+\partial_iV_B(R_A)\,d_A^i(Q)
+\frac12\partial_i\partial_jV_B(R_A)\,M_A^{ij}(Q)
+\cdots.
\]

This already shows why the raw Coulomb well depth is not the first internal load.

Because the density is normalized, the zeroth-order term
\[
V_B(R_A)
\]
shifts only the center-of-mass energy. It is independent of the internal coordinates once normalization is fixed. And because internal coordinates are taken in a centered frame, the first dipole moment is set to zero,
\[
d_A^i(Q_0)=0,
\]
so the first-order term does not drive the internal deformation sector either.

Therefore the first internally active load begins at second spatial derivative order:
\[
\boxed{
H_{A,\rm int}^{\rm drive}
=
\frac12\,T_{ij}(R_A)\,M_A^{ij}(Q)+\cdots,
\qquad
T_{ij}\equiv \partial_i\partial_jV_B.
}
\]
The atomic perturbation variable is thus not the potential itself but the **local Hessian of the partner field across the finite throat**.

### 11.2 Hessian/tidal loading and the first active multipoles

To expose the response channels, decompose the Hessian into trace and tracefree parts:
\[
T_{ij}=T_0\,\delta_{ij}+T_{2,ij},
\qquad
T_0\equiv \frac13\delta^{ij}T_{ij}=\frac13\nabla_3^2V_B,
\qquad
\delta^{ij}T_{2,ij}=0.
\]
Likewise decompose the internal second-moment derivatives around equilibrium:
\[
M_{A,a}^{ij}
\equiv
\left.\partial_{Q^a}M_A^{ij}\right|_{Q_0}
=
\frac13\,\mathcal M^{(0)}_{A,a}\,\delta^{ij}
+\mathcal M^{(2)\,ij}_{A,a},
\qquad
\delta_{ij}\mathcal M_{A,a}^{(2)\,ij}=0.
\]
Then the linear internal drive becomes
\[
\delta H_{A,\rm drive}
=
G_{A,a}^{(0)}T_0\,\delta Q^a
+G_{A,a\alpha}^{(2)}T_{2,\alpha}\,\delta Q^a,
\]
where \(T_{2,\alpha}\) denotes any real quadrupole basis decomposition of \(T_{2,ij}\), and the response couplings are the moment overlaps
\[
G_{A,a}^{(0)}=\frac12\mathcal M_{A,a}^{(0)},
\qquad
G_{A,a\alpha}^{(2)}
=\frac12 E_{ij}^{(\alpha)}\mathcal M_{A,a}^{(2)\,ij}.
\]

So the directly derived atomic load is not a single scalar number. It is the multiplet
\[
\boxed{
\epsilon_{\rm atom}\sim (T_0,\,T_{2,\alpha}).
}
\]
This is the clean atomic realization of the already-frozen 2PN support structure:

- \(T_0\) excites the **monopole / breathing** channel,
- \(T_{2,\alpha}\) excites the **quadrupole** support channel.

### 11.3 Scalar versus quadrupolar response channels

Let the equilibrium stiffness matrix of the internal coordinates be
\[
K_{ab}
=
\left.\frac{\partial^2H_{A,0}}{\partial Q^a\partial Q^b}\right|_{Q_0}.
\]
Adiabatic elimination gives
\[
\delta Q^a_{\rm ad}
=
-(K^{-1})^{ab}
\left(
G_{A,b}^{(0)}T_0+G_{A,b\alpha}^{(2)}T_{2,\alpha}
\right),
\]
and therefore the bodywise response energy
\[
\delta H_{A,\rm resp}
=
-\frac12\Lambda_A^{(0)}T_0^2
-\frac12\Lambda_{A,\alpha\beta}^{(2)}T_{2,\alpha}T_{2,\beta}
-\Lambda_{A,\alpha}^{(02)}T_0T_{2,\alpha},
\]
with
\[
\Lambda_A^{(0)}
=
G_{A,a}^{(0)}(K^{-1})^{ab}G_{A,b}^{(0)},
\]
\[
\Lambda_{A,\alpha\beta}^{(2)}
=
G_{A,a\alpha}^{(2)}(K^{-1})^{ab}G_{A,b\beta}^{(2)},
\]
\[
\Lambda_{A,\alpha}^{(02)}
=
G_{A,a}^{(0)}(K^{-1})^{ab}G_{A,b\alpha}^{(2)}.
\]

For an isotropic equilibrium defect, the scalar–quadrupole mixing term vanishes and the quadrupole block reduces to
\[
\Lambda_{A,\alpha\beta}^{(2)}=\Lambda_A^{(2)}\delta_{\alpha\beta}.
\]
The full finite-size pair correction is then the sum of the two bodywise response energies.

This is the real output of the perturbation-map step. The current Maxwell + GNLS reduction does not merely say that “finite size matters.” It says exactly what the first internally active loads are and how the response coefficients are built: as moment-overlap / stiffness contractions of the finite throat profile.

---

## 12. Finite-size response of the bound pair

Once the atomic load multiplet is identified, the next question is what radial laws it produces in hydrogen.

### 12.1 First usable finite-size correction to the hydrogenic potential

Before the direct atomic perturbation map was fully sorted, the 1PN bridge already guaranteed one important structural fact: once a scalar perturbation variable \(\epsilon\) is declared, the leading internal-response correction is quadratic in that load,
\[
\delta H_{\rm resp}(\epsilon)
= -\frac12\chi\,\epsilon^2.
\]
That observation, combined with the already-frozen scalar 2PN susceptibility package, gives a useful interim hydrogenic embedding:
\[
\delta V_{\rm fs}(r,\omega)
=
-\Lambda_{\rm fs}\chi_{\rm fs}(\omega)\,\epsilon_A(r,\omega)\epsilon_B(r,\omega).
\]
Under the minimal scalar mouth-load choice
\[
\epsilon_A(r)=\epsilon_B(r)=\frac{g_C}{E_{\rm th}r},
\]
one obtains the first usable pair correction
\[
\delta V_{\rm fs}^{\rm(min)}(r)
=
-\frac{C_{\rm fs}g_C^2}{r^2},
\qquad
C_{\rm fs}
\equiv
\frac{\Lambda_{\rm fs}\chi_{\rm fs}(0)}{E_{\rm th}^2}.
\]
With a normalized \(1s\)-type trial state, this yields
\[
\delta E_{\rm fs}^{\rm(min)}(a)
=-\frac{2C_{\rm fs}g_C^2}{a^2},
\]
and the clean Coulombic minimum shifts to
\[
a_* = \frac{\hbar^2}{\mu g_C}-4C_{\rm fs}g_C.
\]

This provisional inverse-square correction was useful for two reasons. It showed that the current 1PN/2PN support hierarchy already supports a definite hydrogenic response channel, and it fixed the sign in the right direction for a softening of the effective bound well.

But it is only an interim embedding. Once the direct atomic perturbation map is applied, the load structure becomes sharper.

### 12.2 Why the naive short-distance response diverges

Use the direct action-level load map on the hydrogenic potential
\[
V_B(r)
=
-\frac{g_C}{r}
-\frac{g_C}{2r}e^{-mr}+\cdots,
\qquad
m\equiv \frac{2}{\lambda}.
\]

For the pure Coulomb term,
\[
\nabla_3^2\left(\frac1r\right)=0
\qquad (r>0),
\]
so away from the source the Coulomb sector produces **no scalar trace load**. The scalar breathing channel begins only once the finite-thickness Yukawa tower is included:
\[
T_0(r)
=
\frac13\nabla_3^2V_B(r)
=
-\frac{g_Cm^2}{6}\frac{e^{-mr}}{r}+\cdots
=
-\frac{2g_C}{3\lambda^2}\frac{e^{-2r/\lambda}}{r}+\cdots.
\]

The quadrupole channel is different. For the pure Coulomb piece,
\[
\boxed{
T_{2,ij}^{\rm C}
=
-g_C\frac{3n_in_j-\delta_{ij}}{r^3}.
}
\]
So the pure Coulomb field already drives the quadrupolar response. Under isotropy and neglecting scalar–quadrupole mixing, the direct finite-size pair correction is
\[
\delta V_{\rm fs}(r)
=
-\frac12\Lambda^{(0)}T_0(r)^2
-\frac12\Lambda^{(2)}T_{2,ij}(r)T_2^{ij}(r)+\cdots.
\]
For pure Coulomb,
\[
T_{2,ij}^{\rm C}T_{\rm C}^{2,ij}=6\frac{g_C^2}{r^6},
\]
so the genuine direct Coulomb response is
\[
\boxed{
\delta V_{\rm fs}^{\rm C}(r)
=
-3\Lambda^{(2)}\frac{g_C^2}{r^6}.
}
\]
The Yukawa scalar trace contributes the first directly derived inverse-square-type law,
\[
\boxed{
\delta V_{00}^{\rm Y}(r)
=
-\frac{\Lambda^{(0)}g_C^2m^4}{72}
\frac{e^{-2mr}}{r^2}+\cdots
=
-\frac{2\Lambda^{(0)}g_C^2}{9\lambda^4}
\frac{e^{-4r/\lambda}}{r^2}+\cdots.
}
\]
So the direct atomic map says something stronger than the interim scalar ansatz:

- the raw Coulomb field does **not** load the scalar breathing mode,
- pure Coulomb drives the quadrupole channel and yields an inverse-sixth far-field response,
- the first direct inverse-square correction is Yukawa-dressed and therefore thickness-suppressed.

Now the short-distance issue is obvious. If one continues the point-source Hessian law
\[
\delta V_{\rm fs}^{\rm C}(r)\propto -\frac1{r^6}
\]
all the way to the origin and inserts it into a hydrogenic \(1s\) expectation value, one gets
\[
\langle r^{-6}\rangle_{1s}
\sim
\int_0^\infty dr\,r^2|\phi_{1s}|^2\,\frac1{r^6}
\sim
\int_0 dr\,r^{-4},
\]
which diverges.

So the direct perturbation map does not merely produce a new radial law. It exposes the exact place where the point-defect approximation has been pushed beyond its regime of validity.

### 12.3 Identification of the core-regulation problem

The divergence has a very specific origin. It comes from combining two approximations outside their common domain:

1. a **point-source Coulomb kernel**, and
2. a **local Hessian expansion** valid only when the external field varies slowly across the finite defect.

Once the inter-defect separation becomes comparable to the physical throat waist, both approximations fail together. The defect can no longer be treated as a point source, and the partner field can no longer be treated as slowly varying across an infinitesimal core.

That identifies the real next task. The theory does not need an ad hoc short-distance cutoff. It needs a physical core regulator rooted in the finite geometry of the throat itself.

---

## 13. Constant-area throat deformation and the \(P_{22}\) mouth mode

The same tidal tensor that generates the atomic finite-size correction also tells us what deformation the throat mouth wants to adopt.

### 13.1 Mass conservation as a constant-area constraint

Write the circular mouth boundary in polar angle \(\theta\) as
\[
r(\theta)=a[1+u(\theta)],
\qquad |u|\ll 1.
\]
Its cross-sectional area is
\[
A
=
\frac12\int_0^{2\pi}r(\theta)^2\,d\theta
=
\pi a^2+a^2\int_0^{2\pi}u(\theta)\,d\theta+O(u^2).
\]
So fixed area implies the linear constraint
\[
\int_0^{2\pi}u(\theta)\,d\theta=0.
\]
The monopole \(m=0\) Fourier mode is therefore forbidden.

Within the present ontology this is not an arbitrary geometric preference. It is the effective statement of mass/inflow conservation at the mouth. A free isotropic collapse of the entire cross section would choke the supported throat rather than merely deform it.

### 13.2 Area-preserving quadrupole deformation of the throat mouth

Centering removes the dipole sector as well:
\[
\int_0^{2\pi}u(\theta)e^{\pm i\theta}\,d\theta=0.
\]
Now choose the principal-axis angle \(\phi_0\) of the planar tidal tensor so that the traceless planar load is represented by
\[
T_{xx}-T_{yy}=T_2^\perp\cos2\phi_0,
\qquad
2T_{xy}=T_2^\perp\sin2\phi_0.
\]
The linear tidal coupling to the boundary deformation is then
\[
\delta U_{\rm tid}^{(1)}
\propto
-a^2T_2^\perp
\int_0^{2\pi}
 u(\theta)\cos2(\theta-\phi_0)
\,d\theta.
\]
So the first linearly driven non-axisymmetric harmonic is exactly
\[
\boxed{
u(\theta)=\epsilon_2\cos2(\theta-\phi_0).}
\]
Area preservation removes \(m=0\), centering removes \(m=1\), and the traceless tidal Hessian couples first to \(m=2\). Therefore the first allowed tidal squish is the \(P_{22}\) mode.

The natural exact completion of this linear mode is the area-preserving ellipse with semiaxes
\[
a_1=ae^{\sigma},
\qquad
a_2=ae^{-\sigma},
\]
so that the area remains fixed,
\[
A=\pi a_1a_2=\pi a^2.
\]
Its exact boundary can be written as
\[
r_b(\theta)
=
\frac{a}{\sqrt{e^{-2\sigma}\cos^2(\theta-\phi)+e^{2\sigma}\sin^2(\theta-\phi)}}.
\]
Expanding for small \(\sigma\) gives
\[
\frac{r_b(\theta)}{a}
=
1+\sigma\cos2(\theta-\phi)
+\frac{\sigma^2}{4}\bigl(3\cos4(\theta-\phi)-1\bigr)+O(\sigma^3),
\]
so the exact nonlinear mouth remains an ellipse whose leading non-axisymmetric harmonic is indeed the required \(P_{22}\) geometry.

### 13.3 Why the first driven deformation is an ellipse

For a uniform mouth profile over the ellipse, the real quadrupole moments are
\[
Q_c\equiv \langle x^2-y^2\rangle
=
\frac{a^2}{2}\sinh2\sigma\cos2\phi,
\]
\[
Q_s\equiv 2\langle xy\rangle
=
\frac{a^2}{2}\sinh2\sigma\sin2\phi.
\]
Thus the mouth carries a genuine real \(P_{22c},P_{22s}\) source pair with amplitude
\[
Q_{22}=\frac{a^2}{2}\sinh2\sigma.
\]
Coupling the planar tidal tensor to this quadrupole gives the exact static mouth energy
\[
\boxed{
U_{\rm tid}^{\rm mouth}(\sigma,\phi;r)
=
-\frac{a^2}{4}T_2^{\rm eff}(r;a)
\sinh2\sigma\cos2(\phi-\phi_0).
}
\]
For weak field,
\[
U_{\rm tid}^{\rm mouth}
=
-\frac{a^2}{2}T_2^{\rm eff}(r;a)
\sigma\cos2(\phi-\phi_0)+O(\sigma^3).
\]

So the proton-side tide does exactly what the geometry suggests: it compresses the mouth along one axis, fixed area forces compensating stretch along the orthogonal axis, and the first resulting non-axisymmetric deformation is the ellipse.

With a support-sector stiffness
\[
U_{\rm supp}(\sigma)=\frac12k_{22}\sigma^2+O(\sigma^4),
\qquad k_{22}>0,
\]
minimization gives the weak-field ellipticity
\[
\sigma_*(r,\phi)
=
\chi_{22}T_2^{\rm eff}(r;a)\cos2(\phi-\phi_0),
\qquad
\chi_{22}\equiv \frac{a^2}{2k_{22}}.
\]
Substituting back yields the mouth-response energy
\[
\delta V_{\rm mouth}(r,\phi)
=
-\frac12\chi_{22}a^2\,[T_2^{\rm eff}(r;a)]^2\cos^22(\phi-\phi_0).
\]
In the far field,
\[
T_2^{\rm eff}(r;a)\to -\frac{3g_C}{r^3},
\]
so the \(r^{-6}\) law is recovered as a genuine weak-field asymptotic. That is the important refinement: the inverse-sixth law is not wrong, but it only holds while the source can still be treated as effectively pointlike compared with the throat waist.

---

## 14. Finite-throat core regulator

Once the same tidal tensor is treated as acting between **finite** throats rather than mathematical points, the short-distance divergence disappears.

### 14.1 Core-resolved effective tidal load

The correct object near the core is not the point-Hessian but the finite-body convolution
\[
T_{ij}^{\rm eff}(R;a)
=
\int d^3\xi\;\rho_a(\xi)
\,\partial_i\partial_j
\left(-\frac{g_C}{|R+\xi|}\right),
\]
where \(\rho_a\) is a normalized source profile of characteristic width or waist \(a\).

For any smooth finite profile, \(T_{ij}^{\rm eff}\) is analytic at \(R=0\) and its natural scale is at most
\[
T_{ij}^{\rm eff}=O\!\left(\frac{g_C}{a^3}\right)
\]
rather than singular. The point-Hessian divergence is therefore not a real prediction of the defect ontology. It is an artifact of extrapolating a far-field approximation into the overlap regime.

A controlled explicit example is obtained from a Gaussian-smoothed source. The regularized Coulomb potential is
\[
\Phi_a(r)
=
-g_C\frac{\operatorname{erf}(r/\sqrt2a)}{r}.
\]
Its quadrupole tidal amplitude can be written as
\[
T_2^{\rm eff}(r;a)
=
\Phi_a''(r)-\frac{\Phi_a'(r)}{r}
=
\frac{g_C}{a^3}\,\mathcal F\!\left(\frac{r}{a}\right),
\]
with
\[
\boxed{
\mathcal F(x)
=
-\frac{3\,\operatorname{erf}(x/\sqrt2)}{x^3}
+\sqrt{\frac{2}{\pi}}e^{-x^2/2}\left(1+\frac{3}{x^2}\right).
}
\]
This function has the exact asymptotics
\[
\mathcal F(x)\sim -\frac{3}{x^3}
\qquad (x\gg 1),
\]
so the familiar point-source tide is recovered in the far field, while near the core
\[
\boxed{
\mathcal F(x)
=
-\sqrt{\frac{2}{\pi}}\frac{x^2}{5}+O(x^4)
\qquad (x\to 0).
}
\]
Thus the tidal load does not diverge at overlap. It actually softens to zero, because once the partner source is resolved as a finite smooth throat, the field becomes locally isotropic at the center.

### 14.2 Removal of the short-distance divergence

Because the mouth-response energy is quadratic in the effective tide,
\[
\delta V_{\rm mouth}(r)\propto -[T_2^{\rm eff}(r;a)]^2,
\]
we now obtain the correct asymptotic behavior in both regimes:
\[
\delta V_{\rm mouth}(r)
\sim -\frac{g_C^2}{r^6}
\qquad (r\gg a),
\]
but
\[
\delta V_{\rm mouth}(r)
\sim -\frac{g_C^2}{a^{10}}r^4
\qquad (r\ll a).
\]
The old inverse-sixth law therefore survives only as a **far-field asymptotic**. Near the core, the response softens and becomes integrable in the hydrogenic variational problem.

The regulator is not an inserted cutoff. It is the physical statement that a finite throat of waist \(a\) cannot be probed by a point-Hessian beyond that scale.

### 14.3 Atomic meaning of the regulated response

Two support mechanisms already visible in the ontology explain why the smoothing scale \(a\) is physically protected rather than arbitrary.

First, as the minor axis of the area-preserving ellipse shrinks,
\[
a_- = ae^{-\sigma},
\]
a trapped EM zero-mode or cavity contribution rises roughly like an inverse squared confinement scale,
\[
U_\gamma\sim \frac{K_\gamma}{a_-^2}=\frac{K_\gamma}{a^2}e^{2\sigma}.
\]
So narrow squeezing rapidly costs support energy.

Second, the parent matter sector carries the stiff bulk law
\[
P(\rho)=K\rho^5,
\]
which heavily penalizes concentrated compression. At the effective support level this is naturally represented by a rapidly rising even function of the ellipticity, for example
\[
U_{\rm bulk}(\sigma)\sim K_B(\cosh4\sigma-1),
\]
though the exact coefficient remains open.

The net consequence is that the throat has a supported finite waist and a supported finite-mouth deformation sector. That same finite support both regulates the atomic core and prepares the geometric background that later becomes important for the lepton-side \(P_{22}\) story.

---

## What Part III establishes

Part III does establish several things.

1. The first internally active atomic load is the partner-field Hessian across the finite defect, not the raw Coulomb scalar.
2. The direct reduced atomic load multiplet is naturally
   \[
   (T_0,\,T_{2,\alpha}),
   \]
   matching the already-frozen \(P_0\oplus P_2\) support structure.
3. A minimal bridge-level scalar embedding already yields a practical first inverse-square finite-size correction.
4. The direct action-level map sharpens that statement: pure Coulomb drives a quadrupolar inverse-sixth far-field response, while the scalar inverse-square piece enters only through finite-thickness Yukawa structure.
5. The short-distance divergence is not fundamental. It is the sign that the point-source Hessian approximation has been extended past its valid regime.
6. Under constant area and centering, the first tidally driven non-axisymmetric mouth deformation is uniquely the \(P_{22}\) ellipse.
7. Treating the partner as a finite throat of waist \(a\) replaces the divergent point-Hessian by a smooth finite tidal kernel and removes the apparent core singularity.

## What Part III does not yet establish

Part III does **not** yet provide a full solved atomic core theorem from the moving-throat PDE.

It does not yet derive:

- the exact microscopic core profile of the throat,
- the exact numerical response coefficients \(\Lambda^{(0)}\), \(\Lambda^{(2)}\), \(k_{22}\), or the full mixed-sector reduction map,
- the complete relation between the atomic \(P_{22}\) mouth background and the lepton-side internal rotor,
- or the full many-electron extension beyond hydrogen.

Those remain open. But the important point is that the atomic program now has a real response map, a real geometric deformation law, and a real physical core regulator.

---

# Part IV — From Atomic \(P_{22}\) Forcing to the Lepton Doublet

*Draft write-up / handoff artifact for the 4D superfluid program.*

Part III ended with a result that matters far beyond hydrogen. The finite-size atomic response did not merely add a correction to the bound-state potential. It uncovered a concrete non-axisymmetric mouth degree of freedom: under centering and constant-area flow, the first tidally driven deformation of the throat mouth is the real \(P_{22}\) ellipse.

That fact immediately touches the parallel lepton program. The earlier lepton audit had already narrowed the live same-charge possibilities to one very specific corridor: a mixed-sector internal rotor carried by the surviving \((a\text{-}w)\) branch. But that rotor remained continuous in the frozen model. What it lacked was a real static splitter capable of turning a continuous internal angle into a protected two-state sector.

Part III supplied the first physical candidate for exactly that missing ingredient.

This part explains the bridge in four steps.

First, it recalls why the mixed-sector rotor was the only same-charge corridor left alive and why that corridor still stalled without a real \(m=2\) anisotropy source.

Second, it shows that the atomic tide derived in Part III is precisely the missing splitter at the effective level: once the mouth is forced into a static \(P_{22}\) background, the mixed-sector rotor sees the required \(\pi\)-periodic potential.

Third, it addresses the isolated-particle limit. The atomic tide vanishes as the partner is removed, so an atomically induced ellipse alone cannot support a free lepton. To keep the lepton corridor alive, one must show that the \(P_{22}\) sector persists intrinsically. The controlled topological handoff from the parallel lepton program does exactly that by locking the internal Berry flux to the half-integer lattice and, in the minimal branch, to
\[
\alpha=\frac12.
\]

Fourth, it replaces the previously symbolic intrinsic bracing coefficient \(h_\alpha\) with a concrete mixed-sector core stress calculation. The result is that the same mixed branch that already carried the first-order rotor also produces a real static quadrupole stress, and therefore a real intrinsic \(P_{22}\) mouth bracer.

As in the earlier parts, the claim level is the same throughout: these are **reduced-sector consequences of the existing action stack plus controlled closure inputs**, not yet the final theorem of the fully solved moving-throat PDE.

---

## 15. The lepton-side need for a static \(P_{22}\) splitter

By the end of the lepton audit, most naive same-charge routes had already failed.
Ordinary circulation was not intrinsic spin. Static side labels collapsed back into the opposite-charge branch bookkeeping. Purely axisymmetric wall modes never produced a protected Pauli-like doublet. Once those routes were removed, only one corridor remained structurally alive: the **mixed \((a\text{-}w)\) sector**.

### 15.1 Continuous mixed-sector rotor in the frozen model

The surviving reduced ansatz uses the lowest transverse odd pair
\[
(b_x,b_y)
\]
and a same-charge mixed-twist sign
\[
\tau=\pm1.
\]
At low energy the mixed sector reduces to the transverse collective-coordinate model
\[
L_{\rm tr}
=
\frac{M_b}{2}\dot b_i\dot b_i
+\frac{\tau\nu_0}{2}\epsilon_{ij}b_i\dot b_j
-V_{\rm stat}(b_x^2+b_y^2).
\]
If the transverse radius pins at
\[
s_*=b_x^2+b_y^2,
\]
write
\[
b_x=\sqrt{s_*}\cos\varphi,
\qquad
b_y=\sqrt{s_*}\sin\varphi.
\]
Then the dynamics collapse to an internal rotor,
\[
L_{\rm rot}
=
\frac{I}{2}\dot\varphi^2+\tau\kappa_0\dot\varphi,
\qquad
I=M_b s_*,
\qquad
\kappa_0=\frac{\nu_0 s_*}{2}.
\]
Its spectrum is
\[
E_n^{(\tau)}=\frac{(n\hbar-\tau\kappa_0)^2}{2I}.
\]

This was already an important milestone. It was the first place in the corrected ontology where a same-charge internal sign and a genuine first-order internal phase law coexisted in one reduced model without being smuggled in as ordinary orbital motion.

But by itself it was still only a **continuous rotor**.

### 15.2 Why a protected doublet requires a real splitter

A continuous rotor is not yet a lepton-like doublet.
Even if the Berry coefficient is nonzero, the angle \(\varphi\) still lives on a continuous circle unless something breaks the accidental \(O(2)\) degeneracy down to the headless \(m=2\) structure appropriate to a Pauli-like two-state package.

What the frozen model still lacked was a real static anisotropy of the form
\[
-V_2\cos 2(\varphi-\varphi_0).
\]
The lepton audit had already isolated this as the missing discretizer:

1. the mixed sector provides the first-order internal rotor,
2. but without a genuine \(P_{22}\) source, the same-charge sector remains continuously tunable,
3. so no protected doublet follows.

That is why the passive \(P_{22}\) capacity already visible in the 2PN support layer mattered so much. The geometry bundle could already *host* a real quadrupole anisotropy. What it could not yet do was *generate* one.

Part III changed that situation.

---

## 16. Atomic tide as the missing splitter mechanism

The finite-size atomic response derived in Part III supplied the first concrete source for a static \(P_{22}\) mouth background.
The key point was not just that the partner field perturbs the defect. The key point was that, after centering and constant-area flow are imposed, the first tidally driven non-axisymmetric mouth mode is uniquely the \(P_{22}\) ellipse.

### 16.1 Static \(P_{22}\) forcing from the bound-state environment

Let the mouth ellipse be parameterized by an ellipticity \(\sigma\) and an orientation angle \(\phi_m\), with semiaxes
\[
a_+=a e^{\sigma},
\qquad
 a_-=a e^{-\sigma},
\]
so the area \(\pi a_+a_- = \pi a^2\) stays fixed.
The corresponding mouth quadrupole tensor can be written as
\[
Q^{\rm mouth}_{ij}(\sigma,\phi_m)
=
Q_m(\sigma)\,\mathcal N_{ij}(\phi_m),
\qquad
Q_m(\sigma)=\frac{a^2}{2}\sinh 2\sigma,
\]
with the headless quadrupole basis
\[
\mathcal N_{ij}(\phi)
\equiv
n_i n_j-\frac12\delta_{ij},
\qquad
n_i=(\cos\phi,\sin\phi).
\]

In Part III the atomic tide was shown to couple linearly to this mouth quadrupole through the core-resolved tidal load \(T_2^{\rm eff}(r;a)\). The reduced mouth energy in the bound-state environment is therefore
\[
U_{\rm mouth}^{\rm atom}(\sigma,\phi_m;r)
=
\frac12 k_{22}\sigma^2
+\frac14 u_{22}\sigma^4
-\frac{a^2}{4}T_2^{\rm eff}(r;a)\sinh 2\sigma\cos2(\phi_m-\phi_0),
\]
where \(\phi_0\) is the principal axis of the atomic tidal field.

In the weak-forcing limit this gives
\[
\sigma_{\rm atom}(r)
\approx
\frac{a^2}{2k_{22}}T_2^{\rm eff}(r;a),
\]
so the proton-side bound state really does generate a static, finite, non-axisymmetric mouth background.

### 16.2 Coupling between mouth ellipticity and the mixed-sector rotor

The mouth ellipse and the mixed-sector rotor are distinct variables.
The mouth ellipse is a support/background deformation. The rotor angle \(\varphi\) is a low-energy internal collective coordinate carried by the mixed sector.

The lowest symmetry-allowed coupling between them is
\[
U_{\rm mix-mouth}
=
-g_{\rm mix}Q_m(\sigma)\cos2(\varphi-\phi_m).
\]
Once the atomic tide has frozen the mouth into a static \(P_{22}\) background, that coupling becomes a real anisotropy potential for the rotor.
If the mouth responds adiabatically compared with the internal rotor motion, then the rotor sees an effective splitting potential inherited from the mouth background itself.

### 16.3 Effective rotor-plus-splitter dynamics

After adiabatic elimination of the slow mouth variable, the mixed-sector rotor acquires the reduced atomic form
\[
L_{\rm rot}^{\rm atom}
=
\frac{I}{2}\dot\varphi^2
+\hbar\alpha\dot\varphi
-V_2(r)\cos2(\varphi-\phi_0),
\]
with
\[
V_2(r)\propto g_{\rm mix}Q_m\bigl(\sigma_{\rm atom}(r)\bigr).
\]
This is exactly the structural term the earlier lepton audit had been missing.

So Part III did more than regulate the atomic core. It supplied the first physical mechanism that can turn the passive \(P_{22}\) support slots of the throat into a real same-charge splitter.

At this stage, however, the argument is still only a **bound-state** argument.
As the partner is removed,
\[
T_2^{\rm eff}(r;a)\to 0,
\qquad
V_2(r)\to 0.
\]
So the atomic tide alone cannot yet explain a free lepton. Something intrinsic must keep the \(P_{22}\) sector alive as \(r\to\infty\).

---

## 17. Isolated-particle persistence of the \(P_{22}\) sector

The next step is the controlled topological handoff from the parallel lepton program.
The purpose of that handoff is very specific: to prevent the isolated throat from relaxing back into the axisymmetric singlet sector after the atomic partner is removed.

### 17.1 Topological half-flux lock \(\alpha=1/2\)

The lepton-side closure treats the isolated particle as an autonomous self-reproducing GNLS soliton.
Under that closure, the common scalar phase is forced into the trivial class while the mixed sector lives on a non-orientable internal bundle with central sign
\[
U_{\rm loop}=-I_2.
\]
The resulting Berry flux is then constrained to the half-integer lattice,
\[
\alpha\in \mathbb Z+\frac12.
\]
The minimal isolated branch is therefore
\[
\boxed{\alpha=\frac12.}
\]

For the present write-up, this should be read at the same honesty level as the earlier closure inputs: it is a controlled topological result imported from the parallel lepton program, not yet a direct theorem of the fully solved moving-throat PDE in the uploaded core stack.

What it supplies is crucial. The mixed-sector rotor is no longer an arbitrary continuous Berry-shifted circle. It is now the minimal half-flux rotor.

### 17.2 Intrinsic core stress as a permanent quadrupole bracer

Once the mixed sector carries a headless half-flux structure, the first static observable it can generate is not a vector but a traceless rank-2 tensor. At the reduced level this is naturally written as an intrinsic quadrupole bracer,
\[
T^{(\alpha)}_{ij}(\rho)
=
\mathcal T_\alpha f_\alpha(\rho/a)
\left(n_i n_j-\frac12\delta_{ij}\right),
\]
where \(f_\alpha\) is localized near the throat core and \(\mathbf n\sim -\mathbf n\) is a headless orientation.

Coupling this intrinsic stress to the mouth quadrupole gives the isolated mouth energy
\[
U_{\rm iso}(\sigma,\phi_m)
=
\frac12 k_{22}\sigma^2
+\frac14 u_{22}\sigma^4
-h_\alpha\sinh2\sigma\cos2(\phi_m-\phi_\alpha),
\]
where \(\phi_\alpha\) is the intrinsic bracer direction.

The same notation makes it easy to unify the atomic and intrinsic sources. Define the real \(P_{22}\) vectors
\[
\mathbf Q_m(\sigma,\phi_m)
=
Q_m(\sigma)(\cos2\phi_m,\sin2\phi_m),
\]
\[
\mathbf H_\alpha
=
h_\alpha(\cos2\phi_\alpha,\sin2\phi_\alpha),
\]
\[
\mathbf H_{\rm atom}(r)
=
\frac{a^2}{4}T_2^{\rm eff}(r;a)(\cos2\phi_0,\sin2\phi_0).
\]
Then the full mouth energy is
\[
U_{\rm tot}(\sigma,\phi_m;r)
=
\frac12 k_{22}\sigma^2
+\frac14 u_{22}\sigma^4
-
\bigl(\mathbf H_\alpha+\mathbf H_{\rm atom}(r)\bigr)\cdot
\mathbf Q_m(\sigma,\phi_m).
\]
So the atomic tide no longer has to create \(P_{22}\) from nothing. It only modulates a pre-existing intrinsic \(P_{22}\) sector.

### 17.3 Why the isolated throat does not relax back to an axisymmetric singlet

The isolated limit is obtained by removing the partner:
\[
\mathbf H_{\rm atom}(r)\to 0.
\]
The mouth energy then reduces to
\[
U_{\rm iso}(\sigma,\phi_m)
=
\frac12 k_{22}\sigma^2
+\frac14 u_{22}\sigma^4
-h_\alpha\sinh2\sigma\cos2(\phi_m-\phi_\alpha).
\]
If \(h_\alpha\neq0\), the circular configuration \(\sigma=0\) is no longer the minimum.
Minimization in the orientation gives
\[
\phi_m=\phi_\alpha\pmod{\pi},
\]
and in the weak-bracing regime the scalar minimum is
\[
\sigma_\infty
\approx
\frac{2h_\alpha}{k_{22}}.
\]
So the isolated particle remains elliptic even when the external tide vanishes.

That establishes persistence of the \(P_{22}\) mouth sector at the geometric level.
But it also has a sharper dynamical consequence. Once the isolated mouth anisotropy is present, the internal rotor sees the reduced isolated Lagrangian
\[
L_{\rm rot}^{\rm iso}
=
\frac{I}{2}\dot\varphi^2
+\hbar\alpha\dot\varphi
-V_\infty\cos2(\varphi-\phi_\alpha),
\qquad
\alpha=\frac12.
\]
Gauge-transforming
\[
\psi(\varphi)=e^{i\varphi/2}\tilde\psi(\varphi)
\]
removes the Berry shift from the Hamiltonian,
\[
\tilde H
=
-\frac{\hbar^2}{2I}\partial_\varphi^2
-V_\infty\cos2(\varphi-\phi_\alpha),
\]
but the boundary condition becomes antiperiodic:
\[
\tilde\psi(\varphi+2\pi)=-\tilde\psi(\varphi).
\]
Because the intrinsic splitter is \(\pi\)-periodic, the half-turn operator
\[
(S\tilde\psi)(\varphi)=\tilde\psi(\varphi+\pi)
\]
commutes with \(\tilde H\), while antiperiodicity implies
\[
S^2=-1.
\]
The eigenvalues of \(S\) are therefore \(\pm i\), and complex conjugation maps the \(+i\) and \(-i\) sectors into one another at equal energy. Every energy level of the reduced isolated rotor is therefore exactly paired.

So inside the reduced half-flux + \(P_{22}\) model, the isolated same-charge doublet is no longer accidental. It is symmetry-protected.

---

## 18. Microscopic derivation of the intrinsic \(P_{22}\) bracing coefficient

The previous section still left one possible caveat: the intrinsic bracer coefficient \(h_\alpha\) had not yet been shown to be nonzero from a concrete core profile. This section closes that gap.

### 18.1 Mixed-sector core ansatz

Work in a local static gauge near the isolated core,
\[
A_a\approx 0,
\qquad
\partial_t A_w\approx 0,
\]
so the physical mixed field is carried by
\[
F_{aw}=\partial_aA_w.
\]
The lowest nontrivial localized transverse odd representative is
\[
\boxed{
A_w(x,y,w)=b_i x_i\,\chi(r,w),
\qquad
r=\sqrt{x^2+y^2},
}
\]
where \(b_i\) is a constant planar amplitude and \(\chi(r,w)\) is a smooth localized profile.

Write
\[
b_i=\sqrt{s}\,n_i(\phi_\alpha),
\qquad
s=b_kb_k,
\qquad
n_i=(\cos\phi_\alpha,\sin\phi_\alpha).
\]
Under
\[
\phi_\alpha\to \phi_\alpha+\pi,
\]
the representative changes sign, but all quadratic observables remain unchanged. So the first physical static observable must lie in the headless \(m=2\) sector.

### 18.2 Integrated planar stress theorem

The static mixed-sector energy density is
\[
\mathcal E_{\rm mix}
=
\frac{Z(w)}{2\mu_0}
\left[F_{aw}F_{aw}+(\partial_wA_w)^2\right].
\]
The anisotropic planar stress is carried by the traceless part of the \(F_{aw}F_{bw}\) term. Define the integrated traceless core stress
\[
\Pi^{\rm mix}_{ij}
\equiv
\int d^2x_\perp\,dw\;
\frac{Z(w)}{\mu_0}
\left(
\partial_iA_w\partial_jA_w
-\frac12\delta_{ij}\partial_kA_w\partial_kA_w
\right).
\]

For the ansatz above,
\[
\partial_iA_w
=
b_i\chi+(b_kx_k)\frac{x_i}{r}\partial_r\chi.
\]
Angular averaging in the mouth plane yields the exact theorem
\[
\boxed{
\Pi^{\rm mix}_{ij}
=
C_{\rm mix}
\left(b_i b_j-\frac12 s\,\delta_{ij}\right),
}
\]
with
\[
\boxed{
C_{\rm mix}
=
\frac{\pi}{2\mu_0}
\int_{-\infty}^{\infty}dw\,Z(w)
\int_0^{\infty}r\,dr\,(2\chi+r\chi_r)^2.
}
\]
Because the integrand is a square and \(Z(w)>0\) on the localized Maxwell support, one has
\[
C_{\rm mix}>0
\]
for every nontrivial localized core profile.

So the intrinsic mixed-sector static stress is automatically a real traceless rank-2 tensor — exactly the \(P_{22}\) object needed for the mouth bracer.

### 18.3 Derivation of \(h_\alpha>0\)

The reduced mixed-sector dynamics already imply a Berry coefficient
\[
\kappa_0=\frac{\nu_0 s}{2},
\qquad
\alpha=\frac{\kappa_0}{\hbar}=\frac{\nu_0 s}{2\hbar}.
\]
Imposing the topological half-flux lock
\[
\alpha=\frac12
\]
fixes the baseline amplitude to
\[
\boxed{
s_\alpha=\frac{\hbar}{\nu_0}.
}
\]
Substituting this into the stress theorem gives the isolated intrinsic bracer
\[
\Pi^{(\alpha)}_{ij}
=
C_{\rm mix}\frac{\hbar}{\nu_0}\,\mathcal N_{ij}(\phi_\alpha).
\]

Now couple this to the area-preserving mouth quadrupole,
\[
Q^{\rm mouth}_{ij}(\sigma,\phi_m)
=
Q_m(\sigma)\,\mathcal N_{ij}(\phi_m),
\qquad
Q_m(\sigma)=\frac{a^2}{2}\sinh2\sigma,
\]
through the lowest symmetry-allowed scalar transfer
\[
U_{\rm coup}=-\Lambda_Q\Pi^{(\alpha)}_{ij}Q^{\rm mouth}_{ij}.
\]
Using
\[
\mathcal N_{ij}(\phi_1)\mathcal N_{ij}(\phi_2)=\frac12\cos2(\phi_1-\phi_2),
\]
one obtains
\[
U_{\rm coup}
=
-h_\alpha\sinh2\sigma\cos2(\phi_m-\phi_\alpha),
\]
with the explicit microscopic coefficient
\[
\boxed{
h_\alpha
=
\frac{\Lambda_Q a^2 C_{\rm mix}\hbar}{4\nu_0}.
}
\]
Thus, provided the core profile is nontrivial and the mouth-core transfer coefficient \(\Lambda_Q\) does not vanish accidentally,
\[
\boxed{h_\alpha>0.}
\]

So the intrinsic \(P_{22}\) source is no longer a placeholder. It is a concrete static stress carried by the same mixed sector that already carried the first-order rotor.

### 18.4 Permanent isolated ellipticity

With the microscopic coefficient in hand, the isolated mouth energy is
\[
U_{\rm iso}(\sigma,\phi_m)
=
\frac12 k_{22}\sigma^2
+\frac14 u_{22}\sigma^4
-h_\alpha\sinh2\sigma\cos2(\phi_m-\phi_\alpha).
\]
Minimization gives
\[
\phi_m=\phi_\alpha\pmod{\pi},
\]
and in the weak-bracing regime,
\[
\boxed{
\sigma_\infty
=
\frac{2h_\alpha}{k_{22}}+O(h_\alpha^2)
=
\frac{\Lambda_Q a^2 C_{\rm mix}\hbar}{2k_{22}\nu_0}+O(h_\alpha^2).
}
\]
So the isolated throat remains elliptic even after the atomic partner is removed.

This sharpens the conceptual picture from the end of Part III. The atomic sector is not required to create \(P_{22}\) from nothing. It only excites and modulates an intrinsic \(P_{22}\) branch that is already present once the topological half-flux closure and the mixed-sector core stress are taken into account.

---

## What Part IV establishes

Part IV does establish several things.

1. The mixed \((a\text{-}w)\) sector remains the only live same-charge corridor in the corrected frozen model.
2. The earlier lepton obstruction was precise: the mixed rotor existed, but it lacked a real static \(P_{22}\) splitter.
3. The finite-size atomic tide derived in Part III supplies the first concrete effective splitter by forcing the mouth into a static \(P_{22}\) ellipse.
4. A controlled topological half-flux closure keeps the isolated rotor on the minimal branch
   \[
   \alpha=\frac12.
   \]
5. The isolated mixed core produces a real intrinsic traceless quadrupole stress, and therefore a real intrinsic bracing coefficient \(h_\alpha>0\).
6. As a result, the isolated throat does not relax back into the axisymmetric singlet sector; it remains in the \(P_{22}\) branch.
7. In the reduced half-flux + \(\pi\)-periodic splitter model, the isolated same-charge sector supports exact spectral doublets.

## What Part IV does not yet establish

Part IV does **not** yet prove a complete lepton theorem from the fully solved parent PDE.

In particular, it does not yet derive:

- the full moving-throat proof of the topological half-flux lock,
- the exact numerical values of \(\nu_0\), \(\Lambda_Q\), \(k_{22}\), \(u_{22}\), \(I\), or \(V_\infty\),
- the full dissipative/open-system stability of the isolated doublet,
- or the magnetic-moment package itself.

Those are still open. But the crucial structural gap has now been closed: the same-charge corridor is no longer missing its discretizer.

---

# Part V — Dirac \(g=2\) from Geometric Gear Reduction

*Draft write-up / handoff artifact for the 4D superfluid program.*

Part IV ended with the first viable isolated lepton package in the corrected ontology.
The same-charge sector now contains three ingredients that had not previously coexisted in a single reduced model:

1. an ordinary observable charged brane circulation mode;
2. an intrinsic mixed-sector rotor locked to the minimal half-flux branch
   \[
   \alpha=\frac12;
   \]
3. and a real static \(P_{22}\) mouth bracer that protects the reduced doublet without destroying the low-energy loop geometry.

That makes the next question unavoidable.
Once the isolated defect carries both a charged brane loop and an intrinsic half-flux spin package, does the mismatch between those two structures force the tree-level Pauli / Dirac target
\[
\mu = g\,\frac{q}{2m}S,
\qquad
 g=2?
\]

This part shows that, within the current reduced closure, the answer is yes.
The factor of two does not need to be inserted by hand.
It emerges from a geometric gear reduction:

- the **external magnetic moment** is carried by an ordinary integer \(2\pi\) brane circulation mode;
- the **internal spin** is carried by the half-flux mixed-sector tether;
- and the last remaining factor is fixed by a common-carrier theorem equating the magnetic circulation inertia with the particle’s reduced inertial ledger.

So the arc of this part is simple.

First, we compute the magnetic moment of the charged \(P_{22}\) ellipse and show that the leading moment law is independent of the ellipticity itself.

Second, we identify the intrinsic spin in the isolated lepton package with the already-locked half-flux branch.

Third, we compare the external loop law to the Pauli form and isolate the exact reduced gyromagnetic ratio.

Fourth, we derive the mass-ledger matching condition
\[
M_\mu=M_{\rm part}
\]
from a common-carrier closure on the moving throat support.
That closes the tree-level Dirac bridge and simultaneously exposes the single mismatch parameter
\[
\zeta
\]
that will become the anomaly channel in Part VI.

As in the earlier parts, the claim level remains the same throughout: these are **reduced-sector consequences of the existing action stack plus controlled closure assumptions**, not yet final theorems of the fully solved moving-throat PDE.

---

## 19. Magnetic moment of the charged elliptical loop

The first task is to separate the role of the static \(P_{22}\) bracer from the role of the observable charged circulation.
The bracer stabilizes the isolated mouth shape and protects the reduced doublet, but it is not itself the magnetic moment.
The magnetic moment comes from the ordinary brane charge moving around the mouth.

### 19.1 Area integral of the \(P_{22}\) ellipse

Take a uniform applied magnetic field perpendicular to the mouth plane,
\[
\mathbf B = B_z\,\hat{\mathbf z}.
\]
In symmetric gauge,
\[
\mathbf A = \frac12\,\mathbf B\times \mathbf r,
\]
the minimal coupling of the observable brane charge \(q_{\rm eff}\) to a closed planar trajectory \(\mathbf r(t)=(x(t),y(t))\) is
\[
L_B
=
q_{\rm eff}\,\mathbf A\cdot \dot{\mathbf r}
=
\frac{q_{\rm eff}B_z}{2}(x\dot y-y\dot x).
\]
Averaging over one period \(T\), the Zeeman energy is
\[
U_B=-\mu_z B_z,
\]
with magnetic moment
\[
\mu_z
=
\frac{q_{\rm eff}}{2T}\int_0^T (x\dot y-y\dot x)\,dt
=
\frac{q_{\rm eff}}{2T}\oint (x\,dy-y\,dx).
\]
So the moment is controlled by the oriented area swept by the charged loop.

Now evaluate this for the area-preserving \(P_{22}\) ellipse selected in Part IV:
\[
x(\Theta)=a e^{\sigma}\cos\Theta,
\qquad
y(\Theta)=a e^{-\sigma}\sin\Theta,
\qquad 0\le \Theta<2\pi.
\]
Then
\[
x\,dy-y\,dx = a^2\,d\Theta,
\]
so the loop area integral is exactly
\[
\boxed{
\oint (x\,dy-y\,dx)=2\pi a^2.
}
\]
The result is exact and does not depend on the ellipticity \(\sigma\).

### 19.2 Ellipticity independence of the leading magnetic coefficient

The exact area law immediately gives
\[
\boxed{
\mu_z=\frac{q_{\rm eff}\pi a^2}{T}.
}
\]
This is the crucial separation of roles:

- the static \(P_{22}\) background changes the shape of the mouth,
- but as long as the area-preserving constraint holds, it does **not** renormalize the leading current-area coefficient.

So the ellipse is not a correction to the loop law. It is the background geometry on which the charged loop moves.

To make the loop structure explicit, define the time-averaged external circulation angular momentum using the inertial ledger \(M_\mu\) carried by the observable charged loop:
\[
\bar L_{{\rm ext},z}
=
\frac{M_\mu}{T}\int_0^T (x\dot y-y\dot x)\,dt
=
\frac{M_\mu}{T}\oint (x\,dy-y\,dx).
\]
Comparing with the magnetic moment gives the exact charged-loop identity
\[
\boxed{
\mu_z
=
\frac{q_{\rm eff}}{2M_\mu}\,\bar L_{{\rm ext},z}.
}
\]
This is just the usual class-1 charged-loop law, but now written in the defect notation.

Because the observable brane circulation closes on an ordinary \(2\pi\) loop, the circulation sector is an integer mode.
For the \(m\)th loop branch,
\[
\bar L_{{\rm ext},z}^{(m)} = m\hbar,
\qquad m\in\mathbb Z,
\]
so
\[
\boxed{
\mu_z^{(m)}
=
\frac{q_{\rm eff}}{2M_\mu}\,m\hbar.
}
\]
The lowest nontrivial charged branch is the minimal integer mode
\[
m=1.
\]
That ordinary integer brane loop will supply the “orbital” side of the later gear reduction.

---

## 20. Internal spin from the topological half-flux lock

The second ingredient is the isolated mixed-sector rotor derived in Part IV.
There the same-charge sector reduced to a headless internal director with a Berry-offset term and a real \(\pi\)-periodic splitter. The minimal topological branch is locked to
\[
\alpha=\frac12.
\]

### 20.1 External circulation quantum \(m=1\)

Before comparing external and internal sectors, it helps to state the minimal observable branch explicitly.
The external charged circulation is an ordinary brane loop, so its lowest nontrivial mode is the integer quantum
\[
m=1.
\]
This is not yet intrinsic spin. It is the ordinary loop quantum of the observable charged circulation around the mouth.
Its role is to determine the loop magnetic moment through
\[
\mu_z^{(1)}
=
\frac{q_{\rm eff}}{2M_\mu}\,\hbar.
\]

### 20.2 Internal rotor spin \(S=\hbar/2\)

The isolated internal rotor, by contrast, is not an ordinary \(2\pi\) orbital loop.
In reduced form it is described by
\[
L_{\rm rot}^{\rm iso}
=
\frac I2\dot\varphi^2 + \hbar\alpha\dot\varphi
-
V_\infty\cos 2(\varphi-\phi_\alpha),
\]
with the locked half-flux value
\[
\alpha=\frac12.
\]
Because the internal carrier is a headless mixed-sector director rather than a full vector rotor, the Berry offset directly sets the intrinsic spin magnitude:
\[
S_0 = |\alpha|\hbar.
\]
So the minimal isolated branch carries
\[
\boxed{
S_0=\frac{\hbar}{2}.
}
\]
The two reduced spin projections are therefore
\[
S_z=s\,\alpha\hbar,
\qquad s=\pm1,
\]
so on the minimal branch
\[
S_z=\pm \frac{\hbar}{2}.
\]

This is the key structural mismatch that the previous program had been circling around:

- the observable charged loop is an ordinary integer \(2\pi\) circulation mode,
- while the intrinsic spin is a half-quantum mixed-sector tether.

That is where the factor of two will come from.

---

## 21. Reduced gyromagnetic ratio

Once the external loop law and the internal half-flux spin are both written in the same reduced language, the magnetic package becomes almost automatic.

### 21.1 Pauli-form comparison

Start from the exact charged-loop result on the \(m\)th branch:
\[
\mu_z^{(m)}
=
\frac{q_{\rm eff}}{2M_\mu}\,m\hbar.
\]
Now compare this to the target Pauli form for an intrinsic spin package,
\[
\mu_z
=
g_{\rm red}\,\frac{q_{\rm eff}}{2M_{\rm part}}\,S_z.
\]
Using the intrinsic spin magnitude
\[
|S_z|=|\alpha|\hbar,
\]
we obtain the reduced gyromagnetic ratio
\[
\boxed{
g_{\rm red}^{(m)}
=
\frac{|m|}{|\alpha|}\,\frac{M_{\rm part}}{M_\mu}.
}
\]
This is the cleanest formula in the whole bridge.
It isolates the entire problem into two pieces:

1. a **topological/geometric factor**
   \[
   \frac{|m|}{|\alpha|},
   \]
2. a **mass-ledger factor**
   \[
   \frac{M_{\rm part}}{M_\mu}.
   \]

### 21.2 Factor-of-two from \(m=1\) versus \(\alpha=1/2\)

On the minimal charged branch,
\[
m=1,
\qquad
\alpha=\frac12,
\]
so the geometric factor is immediately
\[
\frac{|m|}{|\alpha|}=2.
\]
Therefore the reduced magnetic law becomes
\[
\boxed{
g_{\rm red}=2\,\frac{M_{\rm part}}{M_\mu}.
}
\]
This is the sharp version of the “gear reduction” picture.
The factor of two does not come from a lucky convention. It comes from comparing:

- an **external** integer charged loop with period \(2\pi\),
- to an **internal** half-flux tether with spin magnitude \(\hbar/2\).

So the geometry already wants a Dirac factor.
What remains is the mass-ledger match.

### 21.3 Conditions for the tree-level Dirac target

The exact tree-level Dirac hit follows if the circulation inertia of the observable loop is the same reduced particle mass ledger that defines the one-particle dynamics:
\[
\boxed{
M_\mu=M_{\rm part}.
}
\]
Under that closure,
\[
\boxed{
g_{\rm red}=2.
}
\]
Equivalently, the reduced magnetic law becomes
\[
\mu_z
=
\frac{q_{\rm eff}}{2M_{\rm part}}(2S_z)
=
2\,\frac{q_{\rm eff}}{2M_{\rm part}}S_z,
\]
which is exactly the tree-level Pauli / Dirac structure.

So Part V has now reduced the tree-level question to one precise microscopic closure:

> Does the same throat support branch carry both the observable charge and the reduced inertial mass of the charged loop?

That is the point of the next section.

---

## 22. Inertia matching and the common-carrier theorem

The remaining factor
\[
\frac{M_{\rm part}}{M_\mu}
\]
should not be left as a floating constant.
The current reduced ontology already has enough structure to express it in terms of the support measures carried by the moving mouth/collar branch.

### 22.1 Charge and inertia on a shared throat support measure

Keep the static isolated \(P_{22}\) background with fixed waist \(a\), ellipticity \(\sigma\), and orientation \(\phi\).
Parameterize the mouth support by the similar-ellipse family
\[
\mathbf u(\lambda,\vartheta;\sigma,\phi)
=
R(\phi)
\begin{pmatrix}
\lambda a e^{\sigma}\cos\vartheta\\[4pt]
\lambda a e^{-\sigma}\sin\vartheta
\end{pmatrix},
\qquad 0\le \lambda\le 1,
\]
where \(R(\phi)\) is the planar rotation matrix.
The circulation mode is a phase shift
\[
\vartheta\to \vartheta+\Theta(t),
\]
so the moving support point is
\[
\mathbf x(t)=\mathbf X(t)+\mathbf u(\lambda,\vartheta+\Theta(t);\sigma,\phi).
\]

The exact geometry gives two useful identities.
Because the ellipse is area-preserving,
\[
\left|\partial_\lambda \mathbf u\times \partial_\vartheta \mathbf u\right|_z
=
\lambda a^2,
\]
so the transverse area element is
\[
dA_\perp = \lambda a^2\,d\lambda\,d\vartheta,
\]
independent of \(\sigma\).
And because \(\partial_\Theta=\partial_\vartheta\),
\[
u_x\,\partial_\Theta u_y-u_y\,\partial_\Theta u_x
=
\lambda^2 a^2.
\]
So the oriented area swept per loop phase is again independent of the ellipticity.

Now let the reduced one-particle support carry normalized charge and inertial measures
\[
dq = q_{\rm eff}\,c_q(\lambda,\vartheta,w)\,d\Sigma,
\qquad
d\mu_{\rm eff}=M_{\rm part}\,c_M(\lambda,\vartheta,w)\,d\Sigma,
\]
with
\[
\int c_q\,d\Sigma = 1,
\qquad
\int c_M\,d\Sigma = 1.
\]
Here:

- \(dq\) is the observable brane charge carried by the circulation branch,
- \(d\mu_{\rm eff}\) is the dressed inertial support measure of the same reduced particle,
- \(d\Sigma\) is the carrier volume element,
- and \(M_{\rm part}\) is the translational particle mass ledger.

The reduced translational kinetic term is then
\[
T_{\rm trans}=
\frac12\int |\dot{\mathbf X}|^2\,d\mu_{\rm eff}
=
\frac12 M_{\rm part}\dot{\mathbf X}^2.
\]
If the support is centered,
\[
\int \partial_\Theta \mathbf u\,d\mu_{\rm eff}=0,
\]
so the same inertial measure that defines the particle mass ledger also defines the circulation observables.

### 22.2 Derivation of \(M_\mu=M_{\rm part}\)

For the circulation mode \(\Theta(t)\), each carrier element has planar velocity
\[
\dot{\mathbf u}=\partial_\Theta \mathbf u\,\dot\Theta.
\]
The differential magnetic moment is
\[
d\mu_z
=
\frac12(\mathbf u\times \dot{\mathbf u})_z\,dq
=
\frac12\lambda^2 a^2\dot\Theta\,dq,
\]
so after integration,
\[
\mu_z
=
\frac{q_{\rm eff}a^2\dot\Theta}{2}\,C_{2,q},
\qquad
C_{2,q}
\equiv
\int \lambda^2 c_q\,d\Sigma.
\]
Likewise the differential angular momentum is
\[
dL_z
=
(\mathbf u\times \dot{\mathbf u})_z\,d\mu_{\rm eff}
=
\lambda^2 a^2\dot\Theta\,d\mu_{\rm eff},
\]
so
\[
L_z
=
M_{\rm part}a^2\dot\Theta\,C_{2,M},
\qquad
C_{2,M}
\equiv
\int \lambda^2 c_M\,d\Sigma.
\]
Therefore the exact reduced relation is
\[
\boxed{
\mu_z
=
\frac{q_{\rm eff}}{2M_{\rm part}}\frac{C_{2,q}}{C_{2,M}}L_z.
}
\]
So the effective circulation inertia is
\[
\boxed{
M_\mu = M_{\rm part}\,\frac{C_{2,M}}{C_{2,q}}.
}
\]
This is the honest general answer.
The only possible mismatch comes from a difference between the charge-support profile and the inertial-support profile.

Now impose the **common-carrier closure**:
\[
c_q(\lambda,\vartheta,w)=c_M(\lambda,\vartheta,w)\equiv c(\lambda,\vartheta,w),
\]
or equivalently,
\[
\boxed{
dq=\frac{q_{\rm eff}}{M_{\rm part}}\,d\mu_{\rm eff}
}
\]
pointwise on the moving carrier.
Then
\[
C_{2,q}=C_{2,M}\equiv C_2,
\]
and the second-moment factor cancels identically:
\[
\boxed{
\mu_z = \frac{q_{\rm eff}}{2M_{\rm part}}L_z.
}
\]
Therefore
\[
\boxed{
M_\mu=M_{\rm part}.
}
\]
The nice feature is that this result is profile-independent. It does not care whether the support is disk-like, collar-like, or spread in \(w\), and it survives the \(P_{22}\) ellipticity because the geometric \(\lambda^2 a^2\) factor appears on both sides and cancels.

### 22.3 Tree-level closure \(g=2\)

Once the common-carrier theorem is combined with the half-flux spin package, the reduced tree-level hit is immediate.
Take
\[
m=1,
\qquad
\alpha=\frac12,
\qquad
M_\mu=M_{\rm part}.
\]
Then
\[
\mu_z
=
\frac{q_{\rm eff}}{2M_{\rm part}}\,\hbar
=
\frac{q_{\rm eff}}{2M_{\rm part}}(2S_z)
=
2\,\frac{q_{\rm eff}}{2M_{\rm part}}S_z,
\]
so
\[
\boxed{
g_{\rm red}=2.
}
\]
Within the current reduced closure, that is the tree-level theorem.
The factor of two still comes from the geometric mismatch between the ordinary external loop and the half-flux internal tether, but the last mass-ledger step is no longer left floating. It is fixed by a concrete support theorem.

### 22.4 Definition of the mismatch parameter \(\zeta\)

The same derivation also exposes the exact failure channel.
Define
\[
\boxed{
\zeta \equiv \frac{C_{2,q}}{C_{2,M}}.
}
\]
Then the general reduced magnetic law is
\[
\mu_z
=
\frac{q_{\rm eff}}{2M_{\rm part}}\,\zeta\,L_z,
\qquad
M_\mu=M_{\rm part}\,\zeta^{-1}.
\]
On the minimal half-flux branch,
\[
m=1,
\qquad
\alpha=\frac12,
\]
so
\[
\boxed{
g_{\rm red}=2\zeta.
}
\]
Thus:

- if \(\zeta=1\), the tree-level Dirac target is hit exactly;
- if \(\zeta>1\), the charge carrier sits effectively farther out than the inertial support;
- if \(\zeta<1\), the inertial support sits effectively farther out than the charge carrier.

This is the most useful structural output of Part V.
The tree-level problem is solved, and the anomaly problem is now localized in a single dimensionless support-mismatch factor.

---

## What Part V establishes

At the level of the current reduced closure hierarchy, Part V establishes five things.

1. The leading magnetic moment of the charged mouth loop is controlled by the enclosed area and is independent of the \(P_{22}\) ellipticity itself.
2. The isolated mixed-sector rotor carries intrinsic spin
   \[
   S=\hbar/2
   \]
   on the minimal half-flux branch.
3. Comparing the external loop law to the internal spin law gives the reduced gyromagnetic ratio
   \[
   g_{\rm red}^{(m)} = \frac{|m|}{|\alpha|}\frac{M_{\rm part}}{M_\mu}.
   \]
4. The minimal branch
   \[
   m=1,
   \qquad
   \alpha=\frac12
   \]
   already supplies the topological factor of two.
5. If charge and reduced inertia are carried by the same throat-support branch, then
   \[
   M_\mu=M_{\rm part}
   \]
   and the tree-level Dirac target follows exactly:
   \[
   g_{\rm red}=2.
   \]

This is the first point in the program where the lepton-side structural ingredients finally close into a real tree-level magnetic theorem.

## What Part V does not yet establish

Part V still does **not** solve the full lepton problem.
It does not yet derive:

- the full moving-throat PDE for the isolated defect,
- the exact microscopic transport law beyond the common-carrier closure,
- the anomaly \(g-2\),
- or the full relativistic Dirac equation.

What it does do is sharpen the entire anomaly problem to one precise channel:
\[
\zeta\neq 1.
\]
That is exactly the opening needed for the next part.

---

# Part VI — The Electron Anomaly as a Continuum Self-Dressing Problem

*Draft write-up / handoff artifact for the 4D superfluid program.*

Part V closed the tree-level bridge as tightly as the present reduced hierarchy allows.
Once the charge and inertia ledgers are forced onto a common throat-support branch, the model gives
\[
\zeta=1,
\qquad
M_\mu=M_{\rm part},
\qquad
g_{\rm tree}=2.
\]
So the lepton problem changes character at this point.
The question is no longer whether the corrected 4D ontology can reach the Dirac target.
Within reduced closure, it can.

The next question is subtler and more interesting:

> If tree level is already exact, what microscopic redistribution of support produces the very small measured deviation \(g-2\)?

This part develops the answer as a sequence of increasingly local closures.
The anomaly is interpreted not as a new particle sector and not as a failure of the tree-level bridge, but as a **small self-dressing of the carrier support**:

1. a first Schwinger-scale transfer of visible charge support toward the mouth boundary;
2. a quadratic inertial backreaction from the same \(P_{22}\) self-loop that protects the isolated doublet;
3. a finite-thickness softening of that inertial loading;
4. a local support-transport PDE that generates the blur scale instead of fitting it;
5. and a final charge-side local transport correction that nearly saturates the measured electron value.

The result is the strongest reduced anomaly chain reached in the present program:

- tree level gives a rigid
  \[
  g=2;
  \]
- the first support transfer lands naturally on the Schwinger scale;
- the quadratic \(P_{22}\) self-loop fixes a rational inertial coefficient
  \[
  \eta_1=\frac{11}{36};
  \]
- local finite-thickness and local charge transport then drive the reduced prediction to within
  \[
  2.27\times 10^{-12}
  \]
  of the measured electron magnitude, using no new fit parameter.

As before, the claim level remains disciplined: this part is a **reduced-sector derivation within the declared closure hierarchy**, not yet a final theorem of the fully solved moving-throat PDE.

---

## 23. From exact tree level to the measured anomaly

The anomaly story begins by isolating exactly what must move once the tree-level bridge has been closed.

### 23.1 Why \(\zeta=1\) at tree level

Part V established the common-carrier theorem:
if the observable charge and the reduced inertial ledger are carried by the same throat-support measure, then
\[
C_{2,q}=C_{2,M}
\qquad\Longrightarrow\qquad
\zeta\equiv \frac{C_{2,q}}{C_{2,M}}=1.
\]
The reduced Pauli comparison then gives
\[
g_{\rm red}=2\,\zeta,
\]
so tree level closes to
\[
\boxed{
\zeta_{\rm tree}=1,
\qquad
 g_{\rm tree}=2.
}
\]
That result is already structurally important.
It means the anomaly is not a correction to the factor-of-two mechanism itself.
The factor of two is already present in the mismatch between the ordinary integer brane loop \(m=1\) and the intrinsic half-flux rotor \(\alpha=1/2\).

So the measured electron cannot be telling us that the geometric gear reduction is wrong.
It can only be telling us that the exact common-carrier identity is slightly softened once self-dressing of the support is included.

### 23.2 Experimental target for the electron

Using the same numerical inputs adopted in the derivation notes,
\[
\alpha = 7.2973525643\times 10^{-3},
\qquad
|g_e| = 2.00231930436092,
\]
we define the usual anomaly magnitude
\[
a_e\equiv \frac{|g_e|-2}{2}=1.15965218046\times10^{-3}.
\]
In the present language that corresponds to
\[
\boxed{
\zeta_e = \frac{|g_e|}{2}=1+a_e=1.00115965218046.
}
\]
So the observed electron requires only
\[
\boxed{
\delta\zeta_e \equiv \zeta_e-1 = 1.15965218046\times10^{-3}.
}
\]
The numerical target is therefore very small.
The whole anomaly is a part-per-thousand effect riding on top of an exact tree-level \(g=2\) core.

### 23.3 Size and meaning of the required correction

This smallness matters physically.
In support-language terms, the measured electron does **not** demand a drastic rearrangement of the defect.
It demands only a tiny redistribution of where the visible charge weight sits relative to the effective inertial weight.

That is what makes a self-dressing interpretation plausible.
The target is not “replace the tree-level mechanism.”
It is “slightly separate the charge second moment from the inertia second moment.”

So Part VI proceeds by asking the smallest honest sequence of questions:

1. what is the first carrier-support correction compatible with the frozen ontology?
2. what part of that correction should live on the charge side, and what part on the inertia side?
3. can the first few local transport layers account for the observed part-per-thousand deviation?

The answer turns out to be yes, at least within the present reduced closure.

---

## 24. Schwinger-scale geometric support transfer

The first non-common-carrier correction should be as small and as geometric as possible.
The frozen stack already points away from a bulk-volume renormalization and toward a boundary-support correction: the EM zero mode is brane-visible, the mixed sector survives microscopically, and the first nontrivial response structure is a mouth/support layer rather than a pure scalar monopole.

### 24.1 Collar transfer model for the charge carrier

Take the normalized reference mouth disk with second moment
\[
C_2^{(0)}.
\]
The smallest conservative support deformation is to transfer a tiny fraction \(f\) of the visible charge support from the interior toward a boundary collar while leaving the total carried charge fixed.

In the sharpened transferred-collar geometry used later in this part, the exact charge-side ratio is
\[
\boxed{
Q_{\rm sharp}(f)
\equiv
\frac{C_{2,q}}{C_2^{(0)}}
=1+f-f^2.
}
\]
At leading order,
\[
Q_{\rm sharp}(f)=1+f+O(f^2),
\]
so the first correction is linear in the transferred support fraction.

The natural reduced identification is
\[
\boxed{
f\equiv \frac{\alpha}{2\pi}.
}
\]
This is the smallest loop-sized support transfer available once the topological half-flux sector is already present.

### 24.2 First anomaly correction at order \(\alpha/2\pi\)

If the inertia side is still held fixed at this first step, then
\[
\zeta = 1+f+O(f^2),
\]
so the magnetic moment becomes
\[
\boxed{
g=2\Bigl(1+f+O(f^2)\Bigr)
=2\left(1+\frac{\alpha}{2\pi}+O(\alpha^2)\right).
}
\]
So the first support-transfer correction lands automatically on the Schwinger scale.

The exact transferred-collar geometry sharpens this to
\[
\zeta = 1+f-f^2+O(f^3)
\]
before inertial backreaction is even included.
That refinement becomes important once the quadratic self-loop is added.

### 24.3 Why the one-loop scale appears naturally

The structural point is that the one-loop scale does not appear here as an imported QED diagram.
It appears because the first visible departure from exact common-carrier support is a tiny boundary transfer, and the natural size of that transfer is the already-available loop fraction
\[
f=\frac{\alpha}{2\pi}.
\]

So the first part of the anomaly is not mysterious inside the reduced model.
What remains is to understand why the measured electron is **slightly below** the sharp one-loop charge transfer. That is where inertial backreaction enters.

---

## 25. Quadratic inertial backreaction from the \(P_{22}\) self-loop

The same \(P_{22}\) support layer that protected the isolated doublet in Parts III–IV also gives the first honest inertial self-dressing.
The anomaly therefore stops being a purely charge-side problem. The defect must surf its own induced support wake.

### 25.1 Continuum self-loop interpretation

Once the charge support is pushed slightly outward, the mixed-sector throat response does not remain inert.
The induced \(22\)-support mode loads the mouth and feeds back into the inertial ledger.
That is the continuum analogue of a self-loop: the defect propagates a support disturbance and then responds to the disturbance it created.

The first nontrivial inertial correction is quadratic rather than linear, because the support loading is built from the square of the same \(P_{22}\) mode that already appears in the protected-doublet package.
So the inertial second moment takes the form
\[
\boxed{
C_{2,M}=C_2^{(0)}\left[1+\eta_1 f^2+O(f^3)\right].
}
\]

### 25.2 Canonical \(22\)-support normalization

In the minimal reduced self-loop closure, the \(22\)-support mode is canonically normalized and the inertial share carried by that mode is fixed by the same response bookkeeping already used throughout the throat-response program.
The result is that the quadratic inertial coefficient is not free.
It reduces to a definite rational:
\[
\boxed{
\eta_1=\frac{11}{36}.
}
\]

This is one of the most structurally satisfying outputs of the anomaly program.
The first genuinely nontrivial coefficient is not fitted from the electron.
It drops out of the reduced \(P_{22}\) self-loop once the minimal support identifications are made.

### 25.3 Derivation of \(\eta_1=11/36\)

The derivation note resolves the coefficient into two factors:

1. the canonical \(22\)-self-loop normalization contributes the overall factor
   \[
   \frac12;
   \]
2. the reduced inertial-loading share contributes
   \[
   \frac{11}{18}.
   \]

Multiplying them gives
\[
\boxed{
\eta_1=\frac12\cdot \frac{11}{18}=\frac{11}{36}.
}
\]

With the sharpened charge collar
\[
Q_{\rm sharp}(f)=1+f-f^2,
\]
the reduced anomaly law becomes
\[
\boxed{
\zeta
=1+f-\frac{47}{36}f^2+O(f^3),
}
\]
and therefore
\[
\boxed{
g_{\rm red}
=2\left[1+f-\frac{47}{36}f^2+O(f^3)\right].
}
\]

### 25.4 Resulting reduced anomaly law

Numerically, this already lands very close to the electron target.
Using the same \(\alpha\) and \(|g_e|\) values as above, the sharp self-loop prediction is
\[
|g_e|_{\rm sharp}\approx 2.00231929740804.
\]
The remaining miss is only
\[
\boxed{
\Delta g_{\rm sharp}\approx 6.95\times 10^{-9}.
}
\]
So the anomaly problem is no longer one of order \(10^{-3}\).
By this stage it has already been compressed to a sub-ppb residual on top of the correct tree-level and Schwinger-scale structure.

That remaining miss is the opening for finite-thickness softening.

---

## 26. Boundary leakage and finite-thickness softening

The sharp \(11/36\) self-loop closure still assumes a razor-thin mouth boundary.
But the ontology never required the throat mouth to be perfectly sharp.
Once leakage and healing are admitted, the inertial visibility of the self-loop must be smeared across a narrow collar.

### 26.1 Healing collar and blurred boundary visibility

Let
\[
\delta_{\rm leak}\equiv \frac{\ell_{\rm heal}}{a}
\]
be the fractional thickness of the mouth blur layer.
Instead of seeing the \(22\)-support mode with full sharp visibility, the inertial ledger sees it through a visibility profile
\[
\chi_\delta(\rho).
\]
Then the quadratic inertial coefficient is softened only through the overlap factor
\[
\boxed{
\eta_1(\delta)=\frac{11}{36}\,\mathcal B_\delta,
\qquad
\mathcal B_\delta\equiv \int \chi_\delta\,\varphi_{22}^2\,d\mu_0.
}
\]

For a centered blur collar, the first-order law is universal:
\[
\boxed{
\mathcal B_\delta = 1-3\delta + O(\delta^2).
}
\]
So any small positive boundary thickness softens the inertial loading.

### 26.2 Softening of the inertial-loading share

For the exact linear-ramp collar used in the detailed note,
\[
\boxed{
\mathcal B_{\rm lin}(\delta)
=
1-3\delta+5\delta^2-5\delta^3+3\delta^4-\delta^5+\frac{\delta^6}{7}.
}
\]
The electron target fixes the required softened inertial coefficient to
\[
\eta_1^{(e)}\approx 0.3029782629,
\]
so the required overlap factor is
\[
\boxed{
\mathcal B_e
=\frac{\eta_1^{(e)}}{11/36}
\approx 0.9915652240755.
}
\]
This is only a very small departure from the sharp-boundary value \(1\).

### 26.3 Small boundary thickness needed to reach the electron target

In the linear-collar model, the required blur is
\[
\boxed{
\delta_{\rm leak}^{(e)} \approx 2.82485413698\times10^{-3}.
}
\]
A centered cosine collar gives essentially the same value.
So the remaining discrepancy after the sharp \(11/36\) closure can be absorbed by a healing collar only about
\[
0.282\%
\]
of the throat radius.

This is a physically modest correction, not a structural overhaul.
It tells us that the razor-sharp closure was already almost right.
The remaining work is to stop treating \(\delta_{\rm leak}\) as an effective parameter and derive it from a local transport law.

---

## 27. Local self-transport closure for the blur scale

The next step is to derive the blur dynamically instead of fitting it from the measured electron.
The successful route here is local rather than cosmological: the isolated defect surfs its own support leakage along the mouth collar.

### 27.1 Lowest D/N support mode as the transport clock

The frozen throat stack already contains the lowest mixed D/N half-shifted support scale
\[
\boxed{
k_{1/2}=\frac{\pi}{2L}.
}
\]
This scale is the natural support clock for a local collar-transport problem.
If \(c_s\) is the effective support speed along the collar, then the relaxation rate of the lowest support mode is
\[
\omega_{1/2}=c_s k_{1/2}.
\]
At the same time, the first local mixed transport speed is taken to be the one-loop fraction of that same support speed,
\[
\boxed{
v_{\rm mix}=f c_s,
\qquad
f=\frac{\alpha}{2\pi}.
}
\]

### 27.2 Local PDE for collar transport

Let
\[
x=a(1-\rho)\ge 0
\]
measure inward distance from the mouth edge, and let \(\chi(x,t)\) be the inertial visibility of the induced \(22\)-support mode.
The minimal stationary collar-transport equation is
\[
\boxed{
\partial_t\chi + v_{\rm mix}\,\partial_x\chi
=
\omega_{1/2}(1-\chi),
}
\]
with boundary conditions
\[
\chi(0,t)=0,
\qquad
\chi(\infty,t)=1.
\]
In the steady state, the solution is exponential:
\[
\boxed{
\chi(x)=1-e^{-x/\ell_{\rm blur}},
\qquad
\ell_{\rm blur}=\frac{v_{\rm mix}}{\omega_{1/2}}.
}
\]
Using the support clock above gives
\[
\boxed{
\ell_{\rm blur}=\frac{f}{k_{1/2}}=\frac{2L}{\pi}f.
}
\]
So the blur thickness becomes
\[
\boxed{
\delta_{\rm loc}
\equiv
\frac{\ell_{\rm blur}}{a}
=
\frac{2L}{\pi a}f.
}
\]

### 27.3 Derivation of the blur thickness from local support dynamics

With the frozen aspect ratio
\[
\frac{L}{a}\approx 1.85,
\]
this gives
\[
\boxed{
\delta_{\rm loc}
\approx 1.367846338647\times 10^{-3}.
}
\]
Feeding that into the exact exponential overlap law yields
\[
\mathcal B_{\exp}(\delta_{\rm loc})\approx 0.991848746224,
\qquad
\eta_1^{\rm loc}
=\frac{11}{36}\mathcal B_{\exp}(\delta_{\rm loc})
\approx 0.303064894679.
\]
The resulting local-inertia prediction is
\[
\boxed{
g_{\rm loc}^{\rm (inertia)}
\approx 2.00231930412721,
}
\]
with residual
\[
\boxed{
|g_e|-g_{\rm loc}^{\rm (inertia)}
\approx 2.34\times10^{-10}.
}
\]
So once the blur is derived from a local support PDE, the remaining miss collapses by about another factor of thirty.

But one piece is still being carried forward rather than locally derived: the sharp charge-side collar law.
That is the last local closure needed in this part.

---

## 28. Charge-side transport closure

The inertia side now has a genuine local PDE, so the remaining discrepancy must live on the charge side if the reduced model is to close cleanly.
The smallest honest possibility is a mean-zero local shape mode inside the already-transferred collar.

### 28.1 Mean-zero shape mode inside the transferred collar

Define
\[
\kappa\equiv \frac{2L}{\pi a},
\qquad
\tau(f)=1-\sqrt{1-f},
\]
and keep the already-established global transferred-collar geometry
\[
Q_{\rm sharp}(f)=1+f-f^2.
\]
Inside that collar, resolve the first local charge-shape mode as a D/N half-wave with zero collar mean:
\[
\phi_q(\rho)
=
\cos\!\left(\frac{\pi(1-\rho)}{2\tau}\right)-\bar c(f),
\qquad
\langle \phi_q\rangle_{\rm col}=0.
\]
Because the mode is mean-zero, it cannot disturb the already-fixed tree-level or Schwinger-scale carrier totals.
It only changes how the transferred charge is distributed *within* the collar.

The minimal amplitude is tied to the same support scale already used on the inertia side:
\[
\boxed{
A_q(f)=\frac{\tau}{\kappa}.
}
\]

### 28.2 Local correction to the charge second moment

With this local collar mode, the global charge second moment becomes
\[
\boxed{
Q_{\rm loc}(f)
=
1+f-f^2+2fA_q(f)\,\Xi(f),
}
\]
where \(\Xi(f)\) is the collar-weighted radius-squared overlap of the shape mode.
Its small-\(f\) expansion begins at cubic order,
\[
\boxed{
Q_{\rm loc}(f)
=
1+f-f^2
+
\frac{4-\pi}{\pi^2\kappa}f^3
+O(f^4).
}
\]
So the first genuinely local charge correction arrives at exactly the order that was still missing.
It introduces no new fit parameter.

Combining this with the local inertia factor
\[
\eta_1(f)=\frac{11}{36}\,\mathcal B_{\exp}(\kappa f)
\]
gives the improved reduced anomaly law
\[
\boxed{
g_{\rm loc}(f)
=
2\Bigl[Q_{\rm loc}(f)-\eta_1(f)f^2\Bigr].
}
\]
Its perturbative form begins
\[
\boxed{
g_{\rm loc}(f)
=
2\left[
1+f-\frac{47}{36}f^2
+
\left(
\frac{11}{6}\kappa+\frac{4-\pi}{\pi^2\kappa}
\right)f^3
+O(f^4)
\right].
}
\]

### 28.3 Near-closure of the measured electron \(g\)-factor

With the frozen aspect ratio
\[
\frac{L}{a}=1.85,
\qquad
\kappa\approx 1.17774657888,
\]
the local charge-plus-inertia closure gives
\[
\boxed{
g_{\rm loc}
\approx 2.00231930435865.
}
\]
Compared against the same electron target used throughout this part, the remaining residual is
\[
\boxed{
|g_e|-g_{\rm loc}
\approx 2.27\times10^{-12}.
}
\]
This is the strongest anomaly result reached so far.
Within the current reduced closure, the measured electron is essentially saturated.
The remaining gap now naturally lives at the unresolved
\[
O(f^4)
\]
common charge-inertia transport layer rather than requiring a new mechanism.

---

## What Part VI establishes

Part VI does establish a coherent reduced anomaly chain.

1. The measured electron deviation can be localized to the support-ratio channel
   \[
   \zeta-1.
   \]
2. The first carrier-support transfer naturally appears at the Schwinger scale
   \[
   f=\frac{\alpha}{2\pi}.
   \]
3. The first inertial backreaction is a quadratic \(P_{22}\) self-loop with derived coefficient
   \[
   \eta_1=\frac{11}{36}.
   \]
4. A small finite-thickness boundary layer softens that inertial loading in exactly the direction required by the measured electron.
5. A local support-transport PDE generates the blur scale without introducing a new fit parameter.
6. A final mean-zero local charge-shape mode inside the transferred collar pushes the reduced prediction to within
   \[
   2.27\times10^{-12}
   \]
   of the measured electron magnitude.

So within the present closure hierarchy, the anomaly is no longer a loose phenomenological add-on. It is a concrete continuum self-dressing problem with an explicit reduced solution path.

## What Part VI does not yet establish

Part VI still does **not** prove the full measured electron anomaly as an exact theorem of the fully solved parent PDE.

In particular, it does not yet derive:

- the full moving-throat support transport without reduced closure assumptions,
- the exact common charge-inertia transport layer at
  \[
  O(f^4),
  \]
- the complete spin-coupled and dissipative/open-system completion of the isolated lepton,
- or the full relativistic/QED equivalence class beyond the reduced magnetic moment.

So the correct claim is not “the full electron has been solved exactly.”
The correct claim is stronger and more careful:
within the current reduced hierarchy, the anomaly has been driven into a tiny and sharply identified higher-order remainder.

---

# Part VII — Consolidated Results and Interpretation

*Draft write-up / handoff artifact for the 4D superfluid program.*

Parts II through VI developed the main reduced derivation chain of the current program.
The route began with a hydrogenic reduction of the existing 4D action, passed through finite-size tidal response and throat-mouth regulation, then moved into the isolated lepton package, the tree-level Dirac bridge, and finally the electron anomaly as a small support-redistribution problem.

By this point the project has crossed an important threshold.
The question is no longer whether the corrected ontology can reproduce any familiar physical targets at all.
Within reduced closure, it clearly can.
The more delicate question is how to sort the results into three categories:

1. what now looks structurally locked within the current reduced hierarchy;
2. what remains dependent on closure choices, controlled ansatzes, or not-yet-derived transport coefficients;
3. and what physical interpretation ties the whole sequence together.

That is the purpose of Part VII.
It does not introduce a new mechanism.
Instead it consolidates the discoveries into a single map of the present state of the program.

The central message is this:

> Within the current reduced-sector hierarchy, the 4D superfluid model now supports a coherent chain from the existing 4D action to hydrogenic length scales, finite-throat atomic regulation, a protected lepton doublet, a tree-level Dirac factor, and a near-closure of the measured electron anomaly through local support transport.

At the same time, the present work still stops short of a final first-principles theorem of the full moving-throat theory.
The distinction between those two statements must remain explicit throughout.

---

## 29. Main derived results

This section gathers the main outputs of the preceding parts into their cleanest final forms.
Each result is stated at the same claim level used in the derivations themselves: a reduced-sector consequence of the present action stack plus declared closure assumptions.

### 29.1 Bohr-radius path from the existing 4D action

The first major success of the program was to show that a Bohr-scale orbital radius can arise from the **existing 4D action** without inserting Bohr’s circulation rule by hand.

The starting point is the gauged 4D GNLS matter sector combined with the localized Maxwell reduction.
In the controlled zero-mode regime, the matter field factorizes into a lowest transverse mode and a brane-visible orbital factor,
\[
\psi(\mathbf x,w,t)=\chi_0(w)\,\phi(\mathbf x,t),
\]
which produces an effective three-dimensional hydrogenic variational problem.
The reduced one-parameter energy has the form
\[
E(a)
=
E_\perp
+\frac{\hbar^2}{2\mu a^2}
-\frac{g_C}{a}
+\delta E_{\rm KK}(a)
+\delta E_{\rm NL}(a),
\]
where \(\mu\) is the reduced mass, \(g_C\) is the effective Coulomb coefficient inherited from the reduced Maxwell sector, and the remaining terms encode KK / finite-thickness and inherited GNLS-pressure corrections.

In the clean decoupling limit,
\[
\delta E_{\rm KK},\,\delta E_{\rm NL}\to 0,
\]
minimization gives
\[
\boxed{
 a_* = \frac{\hbar^2}{\mu g_C}.
}
\]
Once the Coulomb coefficient is written in standard brane form,
\[
 g_C = \frac{e_{\rm eff}^2}{4\pi\epsilon_0},
\]
one finds
\[
\boxed{
 a_* = \frac{4\pi\epsilon_0\hbar^2}{\mu e_{\rm eff}^2}.
}
\]

This is the key hydrogenic result of Part II.
The Bohr-radius scale is recovered as the minimum of a reduced action-derived energy, not as a separately imposed phase-winding law.
In that sense, the program established the existence of a real action-based path to the Bohr radius.

The derivation also exposed the first thickness law for charge:
\[
 e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}},
\qquad
 a_* = \frac{4\pi\epsilon_0\hbar^2 Z_{\rm int}}{\mu e_*^2}.
\]
So within the present reduction, the observable orbital scale grows with the effective localization thickness.
This became the first concrete realization of the earlier intuition that charge magnitude on the brane is controlled by defect localization rather than by an independent inserted constant.

### 29.2 Physical origin of the \(P_{22}\) mouth mode

The second major success was not an atomic spectrum result but a geometric one.
Once the finite-throat structure of the bound pair is taken seriously, the first nontrivial atomic loading is not a raw scalar Coulomb well but a **tidal Hessian load** across the finite mouth of the defect.

That shift changes the whole structure of the short-distance problem.
Instead of asking how a pointlike electron moves in a central potential, the reduced theory asks how a finite mouth responds to the spatial variation of its partner’s field.
The relevant load decomposes into scalar and quadrupolar pieces, and the leading shape-sensitive channel is quadrupolar.

The constant-area throat constraint then becomes decisive.
If the mouth is compressed along one direction while preserving the net inflow area of the throat, the first admissible deformation is not a uniform contraction but an area-preserving ellipse.
That is exactly the \(m=2\) mouth mode:
\[
\boxed{
\text{constant area + quadrupole tide }\Longrightarrow P_{22}\text{ mouth deformation.}
}
\]

This did two things at once.
First, it supplied the first physical origin for the \(P_{22}\) geometry that the lepton-side program had been missing.
Second, it turned the apparent short-distance divergence of the naive quadrupole response into a finite-throat problem rather than a pathology.

The resulting core-resolved response showed that the far-field singular law survives only asymptotically.
Near the finite throat waist, the effective tidal load softens and the short-distance plunge is capped.
So the finite-throat regulator was not added externally; it emerged from taking the defect geometry and support thickness seriously.

### 29.3 Protected lepton doublet from intrinsic half-flux bracing

The third major success was the synthesis between the atomic \(P_{22}\) forcing and the isolated lepton program.
Before that synthesis, the mixed-sector rotor still had a continuous internal orientation and therefore lacked a structurally protected two-state package.
The missing ingredient was a real static \(P_{22}\) splitter.

The atomic tide supplied the first environmental splitter.
But the deeper question was whether the isolated defect could keep a \(P_{22}\) sector even after the external partner was removed.
That question was answered by combining the atomic mouth mechanics with the topological half-flux lock
\[
\alpha=\frac12.
\]

Under that lock, the mixed-sector core carries a permanent traceless quadrupolar stress.
The mouth energy takes the reduced form
\[
U_{\rm iso}(\sigma,\phi_m)
=
\frac12 k_{22}\sigma^2 + \frac14 u_{22}\sigma^4
-h_\alpha\sinh 2\sigma\cos 2(\phi_m-\phi_\alpha),
\]
with
\[
 h_\alpha>0
\]
in the adopted microscopic core closure.
So the circular throat is no longer the ground state.
Instead the isolated defect remains in a permanently braced \(P_{22}\) sector,
\[
\boxed{
\sigma_\infty>0.
}
\]

This is the clean reduced origin of the protected lepton doublet.
The doublet is not produced by adding discrete quantum labels by hand.
It is the low-energy remnant of a permanently braced, headless quadrupolar throat with a half-flux mixed-sector tether.

### 29.4 Tree-level Dirac \(g=2\)

Once the isolated defect carries both an ordinary brane-visible charged loop and an intrinsic half-flux rotor, the tree-level gyromagnetic question becomes unavoidable.
The reduced answer turned out to be strikingly simple.

The external magnetic moment is carried by the ordinary observable loop.
For the area-preserving \(P_{22}\) ellipse,
\[
\oint (x\,dy-y\,dx)=2\pi a^2,
\]
so the leading magnetic moment is independent of the ellipticity itself.
The external loop therefore behaves as an ordinary charged circulation mode with integer winding,
\[
\bar L_{{\rm ext},z}^{(m)} = m\hbar,
\qquad m=1\text{ on the minimal branch}.
\]

The internal spin is instead set by the half-flux rotor,
\[
S = |\alpha|\hbar = \frac{\hbar}{2}.
\]
So before any inertia matching, the reduced gyromagnetic ratio reads
\[
 g_{\rm red}^{(m)} = \frac{|m|}{|\alpha|}\frac{M_{\rm part}}{M_\mu}.
\]
For the minimal charged branch \(m=1\) and the locked half-flux \(\alpha=1/2\), the entire problem collapses to the inertia ratio.

The common-carrier theorem then closes that ratio.
If the same throat-support measure carries both the magnetic circulation inertia and the particle’s reduced inertial ledger, then
\[
 M_\mu = M_{\rm part}.
\]
At that point the tree-level result is exact within the reduced closure:
\[
\boxed{
 g_{\rm tree}=2.
}
\]

This is one of the main conceptual outputs of the whole write-up.
The factor of two is not put in by hand and is not treated as a mysterious external quantum rule.
It emerges from a geometric gear reduction between an ordinary integer brane loop and a half-twisted internal topological tether, closed by a common-carrier inertia theorem.

### 29.5 Continuum-mechanical route to the electron anomaly

Once tree level is fixed at \(g=2\), the measured electron anomaly becomes a much more refined question.
It no longer asks for a new mechanism that creates the factor of two.
Instead it asks what slightly separates the charge-support ledger from the inertial-support ledger.

This was expressed through the mismatch parameter
\[
\zeta = \frac{C_{2,q}}{C_{2,M}},
\qquad
 g=2\zeta.
\]
Tree level gives
\[
\zeta_{\rm tree}=1.
\]
So the anomaly is exactly the statement that the visible charge and the effective inertia are not carried by perfectly identical support distributions once self-dressing is included.

The reduced derivation then unfolded in four steps.

First, a Schwinger-scale collar transfer of charge support was introduced with size
\[
 f=\frac{\alpha}{2\pi}.
\]
That already gives the natural one-loop anomaly scale.

Second, the same \(P_{22}\) self-loop that protects the isolated doublet generates a quadratic inertial backreaction,
\[
\eta_1=\frac{11}{36},
\]
so that the leading anomaly law becomes
\[
 g
=
2\left[1+f-\left(1+\frac{11}{36}\right)f^2+\cdots\right].
\]

Third, the sharp-boundary closure was softened by a finite healing collar.
The resulting visibility factor decreases the inertial loading slightly and moves the reduced value much closer to the experimental target.

Fourth, both inertia-side and charge-side local transport were resolved through collar transport equations and local support modes.
That yielded a nearly closed reduced value,
\[
 g_{\rm loc}\approx 2.00231930435865,
\]
leaving only a tiny residual at the next common transport layer.

So the anomaly, in the final reduced interpretation of Part VI, is not a mysterious ultraviolet effect.
It is a **tiny redistribution of carrier support** around a throat whose tree-level magnetic / spin structure is already exact.

---

## 30. What is already a theorem within reduced closure

The phrase “theorem within reduced closure” is important here.
It means the following: once a given reduction and closure set has been explicitly declared, the result follows exactly from that reduced system rather than being fitted ad hoc.

Several of the results reached in Parts II–VI now fall into that category.

### 30.1 Hydrogenic action minimum

Within the zero-mode Maxwell reduction and lowest-transverse-mode factorization, the hydrogenic energy minimum follows directly from the reduced action.
The emergence of the Bohr-scale length
\[
 a_* = \frac{\hbar^2}{\mu g_C}
\]
is therefore already a theorem of the reduced hydrogenic sector.
What is *not* yet a theorem is the full geometry-only determination of the couplings entering that formula.

### 30.2 Reduced-mass emergence in the two-body upgrade

Once the hydrogenic problem is upgraded from a fixed source to a genuine two-defect reduction, the split into center-of-mass and relative coordinates is exact within the reduced kinematics.
So the appearance of the reduced mass \(\mu\) is not fitted; it follows from the two-body structure itself.

### 30.3 Constant-area \(P_{22}\) forcing under quadrupolar load

Within the finite-throat mouth mechanics, a quadrupolar tide together with area preservation forces the first driven mouth mode into the \(P_{22}\) sector.
That link between quadrupole load and area-preserving ellipse is already locked in the reduced mouth geometry.

### 30.4 Finite-throat core regularization

Once the mouth is treated as a finite object rather than a point, the naive short-distance singularity is replaced by a core-resolved effective load.
The removal of the divergence is therefore a theorem of the finite-throat reduction, not an artificial external cutoff.

### 30.5 Reduced isolated doublet under half-flux bracing

Accepting the autonomous-solition half-flux lock
\[
\alpha=\frac12
\]
as the intrinsic topological input, and accepting the adopted microscopic mixed-core stress closure, the existence of a nonzero intrinsic bracing coefficient
\[
 h_\alpha>0
\]
and therefore of a permanent isolated \(P_{22}\) sector is already exact within that reduced model.
So the protected two-state package is no longer only a qualitative hope.
It is a consequence of the reduced lepton mouth-plus-rotor system.

### 30.6 Tree-level Dirac factor under common-carrier closure

Given the minimal charged branch \(m=1\), the half-flux rotor \(\alpha=1/2\), and the common-carrier theorem
\[
M_\mu=M_{\rm part},
\]
the tree-level result
\[
 g=2
\]
is exact within reduced closure.
This is one of the strongest statements reached in the whole program.

### 30.7 Rational quadratic backreaction coefficient

Within the minimal canonical \(22\)-self-loop closure, the inertial backreaction coefficient
\[
\eta_1 = \frac{11}{36}
\]
is a derived rational, not a fit.
The model may eventually modify it once more complete transport physics is included, but within the declared reduced hierarchy it is already fixed.

### 30.8 Local blur scale and local cubic transport law

Once the collar transport PDE is adopted,
\[
\partial_t\chi + v_{\rm mix}\partial_x\chi = \omega_{1/2}(1-\chi),
\]
with the minimal support identifications used in Part VI, the local blur length and the first explicit cubic anomaly coefficient are exact outputs of that local PDE closure.
The same is true of the final charge-side local mode correction used to produce the near-closure of the measured electron magnitude.

---

## 31. What remains approximate or model-dependent

The strengths of the reduced program only become more useful when the open pieces are stated just as sharply.
This section collects the main remaining caveats.

### 31.1 The full moving-throat PDE is not yet solved

Every major result in this manuscript has been derived inside a declared reduction and closure hierarchy.
That is not a weakness if it is kept explicit, but it does mean the full moving-throat PDE has not yet been solved and shown to reproduce all the reduced sectors automatically.

So whenever the text says “derived,” the correct reading is:

> derived within the present reduced hierarchy from the existing action stack plus the explicitly declared closures.

The difference between that and a final theorem of the full dynamical 4D system remains substantial.

### 31.2 Suppressed mixed and off-zero-mode sectors

The hydrogenic reduction suppresses the mixed channels
\[
A_w,\quad J^w,\quad F_{\mu w},
\]
and focuses on the controlled zero-mode Maxwell sector.
That is a clean first step, but it is not the whole microscopic theory.
In particular, several of the later transport and bracing effects are best understood as low-energy remnants of channels that are suppressed in the simplest brane-visible reduction.

### 31.3 Localization profiles and thickness bookkeeping

The reduced couplings depend on localization data such as
\[
Z(w),\quad Z_{\rm int},\quad \chi_0(w),
\]
and several derivation notes used Gaussian or similarly simple explicit profiles to produce closed formulas.
Those profiles are physically reasonable and mathematically useful, but they are not yet derived uniquely from the full throat geometry.
So the thickness-to-charge law is structurally established, while the full geometry-to-profile closure remains open.

### 31.4 Microscopic transport coefficients are not fully closed

By the end of Part VI the anomaly chain depends on several transport-scale identifications that are highly structured but not yet fully microscopic theorems of the parent system.
Examples include:

- the exact mixed-support drift speed used in the local collar PDE;
- the precise local charge-side support mode amplitude;
- the final common charge-inertia transport layer beyond the orders explicitly computed.

The fact that these ingredients can be chosen in a parameter-free and nearly target-saturating way is encouraging.
But they should still be classified as reduced closure relations rather than final microscopic derivations.

### 31.5 The remaining \(O(f^4)\) layer

The final local charge-transport closure drove the reduced electron value to within roughly
\[
2.27\times10^{-12}
\]
of the adopted target.
At that point the remaining discrepancy sits naturally at the next transport layer.
So it would be premature to describe the measured electron as already fully derived exactly.
The honest statement is stronger and narrower:

> the present reduced transport chain nearly saturates the measured electron anomaly without new fit parameters, and the remaining gap is small enough to reside naturally in the next common charge-inertia layer.

### 31.6 Constants not yet derived from pure geometry alone

The present program has demonstrated paths to the Bohr radius and to the electron anomaly, but not yet a complete geometry-only derivation of all microscopic constants such as the particle mass scale or the full observable charge scale.
Those remain long-range targets of the larger program rather than completed outputs of the present manuscript.

### 31.7 Multi-electron and chemistry sectors remain early-stage

Although the atomic finite-size program clarified the geometric origin of the \(P_{22}\) mouth response and made hydrogenic progress concrete, the helium and many-electron sectors are not yet closed at the same standard.
So this manuscript should not yet be read as a derivation of full chemical orbital structure.
It has established the first serious one-electron and finite-size foundations on which that later program could be built.

---

## 32. Physical interpretation of the anomaly as support redistribution

The most useful conceptual summary of Parts V and VI is that the anomaly is **not** the factor of two.
The factor of two belongs to tree level.
The anomaly is what happens when the visible charge support and the effective inertia support stop being carried by exactly the same second-moment measure.

### 32.1 Tree level as exact common-carrier geometry

At tree level, the geometry is simple:

- the external brane-visible charge runs on an ordinary integer loop;
- the internal mixed-sector tether carries a half-flux intrinsic spin;
- the same throat support measure carries both the observable magnetic inertia and the reduced particle inertia.

That gives
\[
\zeta=1,
\qquad
 g=2.
\]
In this interpretation, Dirac’s tree-level result is not a quantum postulate but a statement of exact common-carrier geometry.

### 32.2 The anomaly as \(\zeta>1\)

Once self-dressing begins, the support distributions separate slightly:
\[
\zeta = \frac{C_{2,q}}{C_{2,M}} > 1.
\]
That means the visible charge support sits, on average, just slightly farther out in second-moment space than the effective inertial support.
The anomaly therefore measures a tiny charge-over-inertia bias in the dressed throat.

This interpretation is physically economical.
The measured electron does not require a large distortion of the defect.
It requires only a part-per-thousand support mismatch and, by the end of the present program, only a very tiny next-order correction beyond the already-computed local transport layers.

### 32.3 Why the Schwinger scale appears

In standard field-theory language, the first anomaly correction is often represented as a loop contribution.
In the present language, the same scale appears because the smallest visible departure from perfect common-carrier support is a tiny collar transfer of order
\[
 f=\frac{\alpha}{2\pi}.
\]
So the Schwinger scale is reinterpreted here as the natural magnitude of the first support transfer rather than as the bookkeeping output of a virtual-particle diagram.

### 32.4 Why the quadratic backreaction is tied to the \(P_{22}\) sector

The same throat-response layer that protects the isolated doublet also provides the first nontrivial inertial self-dressing.
That is conceptually satisfying: the anomaly is not grafted onto the spin package from outside.
It is produced by the same finite-throat support mechanics that made the isolated doublet possible in the first place.

In that sense, the anomaly is not an independent add-on to the spin story.
It is a higher-order support consequence of the same \(P_{22}\)-braced defect.

### 32.5 Boundary blur and local transport as physical softening

The remaining refinement of the anomaly came from finite-thickness boundary mechanics and local collar transport.
This matters conceptually because it replaces what could have become abstract fitting freedom with a recognizable physical picture:

- support is transferred toward the boundary;
- the boundary is not infinitely sharp, but carries a small healing collar;
- the mixed-support disturbance transports locally on that collar;
- and the charge-side support also resolves its first local mode instead of remaining perfectly sharp.

So the near-closure of the measured electron is interpreted as the output of a very small, very local support-dressing cascade on top of an already exact tree-level geometric core.

### 32.6 Final interpretive summary

The cleanest one-line summary of the whole lepton chain is therefore:
\[
\boxed{
\text{tree-level }g=2\text{ is exact common-carrier gear reduction,}
\quad
 g-2\text{ is tiny local support redistribution.}
}
\]

That is the conceptual endpoint reached by Parts V and VI.
It is also the most useful bridge into the final part of the write-up, because the remaining work is now easy to state plainly:

- compute the next common charge-inertia transport layer;
- strengthen the microscopic derivation of the local transport coefficients;
- and eventually lift the whole chain out of reduced closure into a fuller moving-throat solution.

---

# Part VIII — Open Problems

*Draft write-up / handoff artifact for the 4D superfluid program.*

The previous seven parts established the strongest reduced-sector chain reached so far in the current program.
Starting from the existing 4D action stack, the derivations produced:

- a controlled hydrogenic reduction with a Bohr-scale variational minimum;
- a finite-throat atomic response theory in which the first active shape mode is the constant-area \(P_{22}\) ellipse;
- an isolated lepton package stabilized by intrinsic half-flux bracing;
- a tree-level geometric derivation of the Dirac target \(g=2\);
- and a local support-transport program that drives the reduced electron result to
  \[
  g_{\rm loc}\approx 2.00231930435865,
  \]
  leaving only
  \[
  |g_e|-g_{\rm loc}\approx 2.27\times 10^{-12}
  \]
  relative to the current NIST/CODATA target.

That is enough progress that the remaining tasks are no longer vague.
They are concrete, localized, and sharply testable.
The purpose of Part VIII is to collect those tasks in the order that now seems most rational.

The important theme is that the open problems are **not all of the same type**.
Some are the next natural layer of the same reduced transport hierarchy.
Some are deeper closure problems that require the full moving-throat PDE rather than another clever reduced ansatz.
And some concern the extension of the present derivation chain beyond hydrogen and the isolated lepton package into genuinely many-body structure.

So the final part is organized as a research program rather than a loose list of caveats.
Each section asks three questions:

1. what remains unfinished;
2. why it matters for the present claims;
3. and what would count as a genuine closure rather than another provisional fit.

---

## 33. Remaining \(O(f^4)\) common charge-inertia transport layer

The most immediate technical open problem is also the most constrained one.
By the end of Part VI, the reduced electron anomaly law had been refined to the point that the remaining discrepancy is naturally of the same order as
\[
 f^4,
 \qquad
 f=\frac{\alpha}{2\pi}.
\]
Numerically,
\[
 f\approx 1.1614\times 10^{-3},
 \qquad
 f^4\approx 1.82\times 10^{-12},
\]
which is exactly the scale of the remaining miss.
That means the next open term is no longer speculative bookkeeping.
It is the first term that is actually required if the reduced anomaly chain is to be promoted from a near-closure to a real measured-value derivation.

The reason this term matters conceptually is that the previous layers were still **partially separated**.
The charge support was first transferred outward through a collar mechanism.
The inertia side was then dressed by the quadratic \(P_{22}\) self-loop and softened by the blurred boundary collar.
Finally, the charge side acquired its own local shape mode inside the collar.
That was enough to reduce the residual to the \(10^{-12}\) scale.
But the next term is expected to be the first one that is genuinely **common**: a transport layer in which the charge and inertia visibility fields can no longer be treated as sequentially corrected ledgers.
They must be propagated together.

### 33.1 Why the residual now points specifically to \(O(f^4)\)

The present reduced law can be written schematically as
\[
 g(f)=2\Bigl[1+f-c_2 f^2+c_3 f^3+O(f^4)\Bigr],
\]
where the coefficients \(c_2\) and \(c_3\) are no longer free.
In the current closure,
\[
 c_2=\frac{47}{36},
\]
and the cubic coefficient is fixed by the local collar transport and local charge-shape mode.
That is why the last visible residual is small enough to sit naturally in the quartic layer.

This is a healthy sign.
A derivation program is in good shape when the uncomputed remainder is roughly the size expected from the next neglected order, not when it is orders of magnitude larger.
Here the remaining discrepancy is not screaming for a new mechanism.
It is pointing toward the next order in the same mechanism.

### 33.2 What the missing quartic layer probably represents

The likely content of the quartic layer is a **common charge-inertia transport equation** rather than one more independent correction on either side.
In the previous steps, the charge and inertia sectors were still updated in a staggered way:

- the charge side was transported outward into the collar;
- the inertia side was then backreacted through the \(22\)-self-loop and blurred by the healing collar;
- and the charge side was finally given its own local mode.

At \(O(f^4)\), this separation is expected to break down.
The next correction likely contains one or more of the following ingredients:

1. a shared collar profile in which the charge and inertial support densities obey coupled rather than independent local transport equations;
2. a common visibility kernel in which the same support perturbation contributes simultaneously to both \(C_{2,q}\) and \(C_{2,M}\);
3. mixed cubic-quadratic interference between the charge-side local mode and the blurred inertial \(P_{22}\) self-loop;
4. and possibly the first finite departure from the still-carried sharp support ratio used for the undressed interior.

In reduced notation, this means replacing the current sequential structure by something like
\[
 \partial_t \mathbf u + \mathcal V(\mathbf u)\partial_x \mathbf u = \mathcal R(\mathbf u),
 \qquad
 \mathbf u=(u_q,u_M),
\]
with \(u_q\) and \(u_M\) the local charge and inertia support fields.
The quartic anomaly coefficient would then be a genuine output of the first coupled local transport system.

### 33.3 What would count as closure here

A genuine closure of this section would not merely be “get the last few digits right.”
The standard should be stricter.
The program should derive a quartic coefficient from a coupled local transport law whose coefficients are already fixed by the existing reduced throat data.

The clean target is:
\[
 g_e = 2\Bigl[1+f-c_2 f^2+c_3 f^3-c_4 f^4+O(f^5)\Bigr]
\]
with
\[
 c_4
\]
computed rather than fitted.
If that happens and the result lands on the measured electron value within the expected next neglected order, then the present anomaly program can honestly be described as a derived reduced-sector result rather than a near miss.

If instead the quartic coefficient must be inserted by hand, or if the coupled charge-inertia PDE forces a coefficient of the wrong sign or wrong order of magnitude, then the present near-closure would have to be reclassified as a numerically suggestive but incomplete closure.

### 33.4 Why this is the best immediate next step

Among all remaining tasks, this one is the most local and the most disciplined.
It does not require solving the full moving-throat PDE.
It does not require introducing a new physical sector.
And it tests the current anomaly chain exactly where it is now weakest: at the first order where truly common support transport should appear.

For that reason, this is the highest-priority next problem if the immediate goal is to convert the present electron result into a full reduced derivation.

---

## 34. Full moving-throat PDE closure

The quartic transport layer is the next tactical step, but the deepest strategic open problem is still the same one that has been present from the beginning:

> the full moving-throat PDE has not yet been solved.

Everything in the present write-up lives inside a declared closure hierarchy.
That hierarchy has been remarkably productive.
It has produced hydrogenic reductions, finite-size throat response, a protected lepton doublet, and a near-closed electron anomaly.
But it is still a hierarchy of reduced sectors, support layers, and carefully separated response channels.
It is not yet a single fully dynamical bulk-brane-defect solution.

### 34.1 Why the reduced hierarchy has reached its natural limit

The reduced derivations were possible because the current program repeatedly separated scales and sectors:

- zero-mode Maxwell on the brane;
- lowest transverse matter mode;
- frozen or adiabatically slaved throat shapes;
- isolated response channels such as \(P_0\), \(P_2\), and the mixed-sector rotor;
- and local collar transport approximations that use fixed mode clocks and effective support lengths.

This is not a flaw.
It is what made the derivation chain possible at all.
But it also means that many of the quantities used in the reduced notes still enter as inherited profile data or reduced transport clocks rather than as direct outputs of a single PDE solution.

The program therefore now sits at a familiar threshold:
reduced closure has probably extracted most of the easy structure from the current action stack.
The next qualitative leap will likely require solving, or at least much more tightly constraining, the genuinely dynamical throat problem.

### 34.2 What “full moving-throat closure” would have to include

A real closure would have to evolve, in a common system, at least the following:

1. the bulk-brane matter field \(\psi(\mathbf x,w,t)\) without freezing it to a single transverse mode from the outset;
2. the throat geometry variables, including the waist scale, length scale, and mouth deformation channels;
3. the localized electromagnetic and mixed fields, including the channels suppressed in the far-field reduction;
4. the support-layer response that currently appears only through effective coefficients such as \(k_{22}\), \(u_{22}\), \(\Lambda_Q\), transport clocks, and visibility profiles;
5. and the continuity / leakage bookkeeping that currently enters through reduced source terms.

In a symbolic shorthand, the desired object is not another effective law for one visible coefficient.
It is something closer to
\[
 \mathcal E[\psi, A_\mu, A_w, g_{\rm throat}, \Sigma_{\rm mouth}] = 0
\]
with all reduced coefficients derived by projection or asymptotic analysis from that common dynamical solution.

### 34.3 What the full PDE would settle immediately

A successful moving-throat closure would do more than simply make the write-up look more complete.
It would settle several of the present reduced ambiguities in one shot.
For example, it would provide direct answers to questions such as:

- Why is the lowest support clock the relevant one in the local collar transport problem?
- Are the adopted visibility profiles exact tails of the true throat boundary PDE or only first asymptotic surrogates?
- Does the same microscopic carrier really support both the charge and inertia ledgers exactly at tree level?
- Is the half-flux locked \(P_{22}\) bracing dynamically selected by the full defect solution, or only admitted by it?
- And what higher support channels first enter beyond the \(P_0\oplus P_2\) bundle?

Those are not cosmetic issues.
They are the questions that separate a very successful reduced program from a fully closed field-theoretic one.

### 34.4 What would count as meaningful progress before a full solution

A full exact solution may be too much to demand immediately.
But there are intermediate achievements that would already count as major progress:

1. deriving the first nontrivial moving-throat mode equations around the isolated or hydrogenic background;
2. deriving the support clocks and visibility tails from that dynamical linearization instead of choosing them by reduced closure;
3. and obtaining a controlled asymptotic expansion that demonstrates which reduced notes are asymptotic theorems of the full PDE rather than clever closures.

So while full moving-throat closure is the deepest open problem, it also contains a ladder of intermediate milestones.
The program does not need to leap directly from reduced notes to an exact global solution.
But it does need to begin collapsing the reduced coefficients back into a common dynamical origin.

---

## 35. Microscopic derivation of all carrier and stiffness coefficients

Even inside reduced closure, several key coefficients still enter as partially inherited, partially modeled data rather than as outputs of a single microscopic throat profile.
This is the third major open problem.

The issue is not that every coefficient is equally uncertain.
Some are already sharply anchored inside the current hierarchy.
Others are only fixed once one chooses a particular collar profile, core ansatz, or local transport clock.
The goal of this section is to separate those cases clearly.

### 35.1 Coefficients that are already relatively well anchored

Several quantities are already much better constrained than they were at the beginning of the project.
Examples include:

- the tree-level support matching condition behind \(M_\mu=M_{\rm part}\);
- the half-flux lock \(\alpha=1/2\) in the reduced isolated lepton closure;
- the canonical quadratic self-loop coefficient
  \[
  \eta_1=\frac{11}{36}
  \]
  in the minimal \(P_{22}\) inertial backreaction closure;
- the lowest D/N half-shifted support length that enters the local collar transport picture;
- and the exact geometric statement that the leading magnetic moment is independent of \(P_{22}\) ellipticity.

Those are not fully microscopic in the sense of coming from a solved throat profile, but they are not loose phenomenological parameters either.
They are already tightly connected to the present reduced hierarchy.

### 35.2 Coefficients that remain only partially derived

Other quantities are still only partially closed.
These include, in one form or another:

- the full microscopic values of the mouth elastic coefficients such as \(k_{22}\) and \(u_{22}\);
- the bracing and coupling amplitudes that enter the static \(P_{22}\) sector;
- the visibility-tail parameters associated with the blurred boundary and local collar transport;
- the precise support-transfer amplitudes that map mixed-sector activity into charge or inertia redistribution;
- the exact throat aspect ratio and its stability across the relevant reduced sectors;
- and any higher support-channel coefficients beyond the first local transport modes already used.

In the anomaly program, the remaining quartic layer is the most visible reminder of this fact.
But the same issue also appears in the hydrogenic and lepton sections whenever a support stiffness, response amplitude, or transport clock is inherited from a prior reduced note rather than solved directly.

### 35.3 Why unifying the coefficients matters

The point of deriving these coefficients microscopically is not aesthetic tidiness.
It is structural economy.
Right now the program succeeds because the same handful of themes keep reappearing:

- finite throat support rather than point defects;
- the \(P_{22}\) mouth sector as the first nontrivial shape mode;
- localized charge determined by thickness or support concentration;
- and local transport on short support collars.

A truly satisfying closure would show that all of the recurring coefficients in those stories descend from a **single microscopic throat profile** and a single equation of state.
That would sharply reduce the sense that the program is moving from one tailored reduced note to another.
It would reveal that the reduced notes are just different asymptotic faces of the same common defect solution.

### 35.4 The clearest next target in this category

The cleanest next target is not “derive every constant at once.”
It is to pick the few coefficients that appear repeatedly in the strongest derivation chain and derive them from a common microscopic profile.
The best candidates are:

1. the support length / transport clock used in the local collar PDE;
2. the blurred-boundary visibility scale entering the inertial softening;
3. the charge-side local collar mode amplitude;
4. and the stiffness ratio that sets the strength of the \(P_{22}\) bracing and response.

If those four pieces can be shown to descend from one common microscopic defect profile, then the whole reduced anomaly program will become much more unified.

---

## 36. Extension from hydrogen to helium and multi-electron structure

The present document achieved something real for hydrogen.
It established that the existing 4D action stack supports a plausible and explicit path to the Bohr-radius scale without inserting Bohr quantization by hand.
It also identified the first physical finite-size correction channels and the first nontrivial mouth deformation under a bound-state tidal load.

But that is still only the beginning of an atomic theory.
The next truly atomic open problem is not another refinement of hydrogen.
It is **helium and beyond**.

### 36.1 What hydrogen has actually solved

The hydrogen program solved a very specific question:
can a single light defect bound to a heavy opposite branch in the reduced Maxwell + GNLS stack produce a Bohr-scale variational minimum with finite-throat correction channels?

The answer appears to be yes.
That is important, but it does not yet amount to a theory of atomic orbitals in the full many-electron sense.
Hydrogen is still the one-body or effective two-body sector of the larger problem.

### 36.2 Why helium is the real next atomic test

Helium is the first place where the reduction can no longer hide behind central symmetry and a single Coulomb partner.
A real helium derivation would have to address, at minimum:

- two light same-charge defects in the same nuclear well;
- finite-throat repulsion between those defects;
- support-mediated shape response when both throats are active simultaneously;
- and a physically honest replacement for any fitted effective screening parameter.

In other words, helium is the first place where the theory must stop being merely hydrogenic and become genuinely many-body.
That is why it is the right next atomic milestone.

### 36.3 What is still missing before helium can be claimed

Several ingredients are not yet ready for a trustworthy helium derivation.
The most important are:

1. a stronger same-charge defect-defect interaction kernel;
2. a clearer spin / exchange-like completion of the lepton package;
3. and a multi-throat extension of the finite-mouth response theory.

Without those pieces, helium would risk becoming a disguised fit.
The program would simply invent a new geometric coefficient to replace the effective screening constant it was trying to transcend.
That would be a step backward.

### 36.4 What a successful helium program would look like

A real helium success would not merely be hitting a single number.
It would show that the same finite-throat machinery already used in hydrogen and the lepton sector can produce:

- a non-ad hoc two-electron ground-state geometry;
- a total energy that lands in the right range without fitting an effective screening constant by hand;
- and a clear explanation of why multi-electron structure requires a three-dimensional throat-mouth organization rather than a naive planar picture.

That is also the first place where the present \(P_{22}\) machinery might begin to connect back to the earlier intuition that three-dimensional orbital structure could emerge from finite-throat repulsion and support packing, rather than from probabilistic clouds taken as primitive.

### 36.5 Why this section should wait until after the local transport closure

Even though helium is a major target, it is probably **not** the next best task.
The anomaly program is now extremely close to a reduced closure, while helium still lacks several sectors outright.
So the correct order of work is likely:

1. finish the local charge-inertia transport closure for the electron anomaly;
2. strengthen same-charge and spin-coupling completion;
3. then return to helium with those tools in hand.

That ordering avoids forcing the atomic many-body program to run ahead of the lepton-side structures it almost certainly needs.

---

## 37. Stronger same-charge and spin-coupling completion

The final major open problem is the one that sits underneath both the helium program and the deeper lepton claims: a fuller same-charge and spin-coupling completion.

The present write-up reached impressive results for the isolated lepton package.
It found a protected reduced doublet, a tree-level Dirac factor, and a highly accurate anomaly chain.
But those achievements still live in a controlled reduced setting.
They do not yet amount to a full theory of all same-charge interactions, all spin couplings, or all many-body spin structure.

### 37.1 Same-charge interactions remain only partially closed

The lepton package solved the isolated defect and the self-dressed defect better than the defect-defect same-charge problem.
That is a very specific asymmetry in the present program.
It means the theory currently knows more about how one defect dresses itself than about how two like defects fully interact.

A stronger completion would need a same-charge kernel that includes, in a common language:

- the mixed-sector rotor and its static \(P_{22}\) background;
- finite-mouth deformation under mutual loading;
- local transport and support overlap;
- and any repulsive or avoided-overlap mechanisms that play the role of a true same-charge exclusion structure.

That kernel is likely to matter not only for helium, but also for whether the isolated lepton package really scales into a many-lepton theory.

### 37.2 Spin couplings to external fields are not yet fully derived

The current chain produced the reduced gyromagnetic ratio and a local self-dressing of the magnetic moment.
But that is not the same thing as a full spin-coupling theory.
Several obvious targets remain:

- Zeeman-type response in an applied external field;
- spin-orbit coupling in a hydrogenic or multi-electron bound state;
- precession dynamics of the mixed-sector rotor under a genuinely time-dependent perturbation;
- and the response of the protected doublet once the static reduced closures are relaxed.

These are not optional embellishments.
If the program aims to claim a real replacement-level account of lepton spin, then it must eventually produce the standard field-response structures, not only the static isolated value of \(g\).

### 37.3 Statistics and many-body spin structure remain beyond the present closure

A protected two-state package is a major milestone, but it is not yet a derivation of full fermionic many-body behavior.
There is still a significant gap between:

- an isolated defect with a topologically braced reduced doublet;

and

- a many-defect theory with the correct exchange structure, same-charge repulsion, and spin-dependent state organization.

The current write-up should therefore be careful not to overstate what has been achieved.
The program has established a plausible continuum origin for a protected two-state lepton package.
It has **not yet** established the whole many-body statistics story.

### 37.4 What would count as decisive progress here

The strongest next tests in this category would be:

1. derive a same-charge two-defect kernel that includes the mixed-sector and finite-mouth response without new ad hoc coefficients;
2. derive the response of the reduced doublet to an external magnetic field and recover the expected low-field spin splitting;
3. show that the same package remains stable when embedded in a genuinely interacting two-defect or few-defect problem;
4. and identify whether the reduced doublet extends naturally into a full exchange / exclusion structure or whether an additional topological ingredient is required.

If those targets are met, then the lepton program will move from “isolated reduced package” to a real same-charge spin theory.
If they fail, then the present isolated results may still be valuable, but they would have to be reinterpreted as a partial sector rather than a full lepton completion.

---

## Final research priority ladder

The open problems above are not equally urgent.
The clearest order of operations now appears to be:

1. **Finish the \(O(f^4)\) common charge-inertia transport layer.**  
   This is the nearest, cleanest, most falsifiable step toward turning the near-closed electron result into a full reduced derivation.

2. **Collapse more of the transport and stiffness coefficients into a common microscopic throat profile.**  
   This strengthens the anomaly chain without yet demanding a full global PDE solution.

3. **Begin lifting the strongest reduced sectors into the moving-throat PDE.**  
   The goal is not an instant full exact solution, but a controlled dynamical derivation of the support clocks, visibility tails, and response coefficients that the reduced notes currently inherit.

4. **Strengthen same-charge and spin-coupling completion.**  
   This is likely a prerequisite for a trustworthy multi-electron atomic program.

5. **Return to helium and multi-electron structure with those tools in hand.**  
   Only then can the program test whether finite-throat response, same-charge support exclusion, and three-dimensional mouth organization really generate multi-electron structure without ad hoc screening fits.

This priority ladder also clarifies the status of the present write-up.
The strongest immediate opportunity is still the electron anomaly, because the local reduced derivation is already very close to closure.
Hydrogen has already crossed the “path exists” threshold.
Helium and many-body structure remain longer-range goals that depend on solving the same-charge and spin-coupling problems more fully.

---

## Closing summary

The most honest summary of the current state is:

> The program has moved from speculative ontology to a disciplined reduced derivation chain, but the final conversion of that chain into a full field-theoretic closure still depends on one last local anomaly layer, a stronger microscopic coefficient unification, and eventually the moving-throat PDE itself.

That is not a weak endpoint.
It is a good endpoint for a handoff.
The central derivation chain is already in place.
The remaining work is now organized into concrete, falsifiable problems rather than diffuse aspirations.

---

# Appendices — Technical Reference and Handoff Map

*Draft write-up / handoff artifact for the 4D superfluid program.*

These appendices collect the technical material that sits underneath Parts I–VIII of the write-up.
The goal is not to introduce a new mechanism.
It is to give the next session a compact reference for notation, reduced equations, variational identities, geometric lemmas, transport closures, benchmark coefficients, and file-to-section provenance.

Throughout, the same claim level used in the main text is maintained:

> all formulas below are stated at the **reduced-sector / declared-closure** level reached in the current program, unless explicitly marked otherwise.

A recurring notational warning is important enough to state up front.

## A note on the two different alphas

Earlier notes used the symbol \(\alpha\) in two different ways:

1. **Topological half-flux / Berry lock**
   \[
   \alpha_{\rm top} = \frac{\kappa_0}{\hbar},
   \]
   with the isolated lepton closure fixing
   \[
   \alpha_{\rm top} = \frac12.
   \]

2. **Fine-structure constant**
   \[
   \alpha_{\rm fs} \approx 0.0072973525643.
   \]

To avoid ambiguity, these appendices use \(\alpha_{\rm top}\) and \(\alpha_{\rm fs}\) explicitly.
The reduced anomaly expansion parameter is then
\[
f \equiv \frac{\alpha_{\rm fs}}{2\pi}.
\]

---

## Appendix A — Notation and symbol table

### A.1 Coordinates, geometry, and support variables

| Symbol | Meaning |
|---|---|
| \(t\) | Time coordinate. |
| \(\mathbf x=(x,y,z)\) | Three visible brane coordinates. |
| \(w\) | Extra / transverse coordinate used in the 4D bulk-brane localization picture. |
| \(r\) | Ordinary radial distance in the brane-visible sector. |
| \(R\) | Two-body center-of-mass coordinate. |
| \(\mathbf r = \mathbf x_e-\mathbf x_p\) | Two-body relative coordinate. |
| \(a\) | Throat-mouth radius / reference mouth scale. In Part II this also labels the hydrogenic variational scale. Context determines usage. |
| \(L\) | Throat-support length used in the D/N half-shifted support ladder. |
| \(\sigma\) | Area-preserving mouth ellipticity variable. |
| \(\phi_m\) | Orientation angle of the mouth ellipse. |
| \(\phi_\alpha\) | Orientation angle of the intrinsic mixed-sector quadrupole bracer. |
| \(\rho\) | Density in the GNLS matter sector. |
| \(\rho_\*\) | Boundary / support reference density in healing-length estimates. |
| \(\lambda\) | Localization thickness scale of the Maxwell zero mode. |
| \(\delta\) | Boundary blur thickness in units of the throat radius, \(\delta = \ell_{\rm blur}/a\). |
| \(\kappa\) | Reduced blur coefficient \(\kappa = 2L/(\pi a)\). |
| \(\tau(f)\) | Collar width fraction \(\tau(f)=1-\sqrt{1-f}\). |

### A.2 Fields, couplings, and reduction data

| Symbol | Meaning |
|---|---|
| \(\psi(\mathbf x,w,t)\) | Full matter field in the gauged 4D GNLS sector. |
| \(\phi(\mathbf x,t)\) | Brane-visible reduced matter field after transverse-mode separation. |
| \(\chi_0(w)\) | Lowest transverse localization mode used in the hydrogenic reduction. |
| \(A_\mu\) | Brane-visible gauge field components. |
| \(A_w\) | Mixed / transverse gauge component. |
| \(F_{\mu\nu}\), \(F_{\mu w}\) | Gauge field strengths. |
| \(q_*\) | Microscopic branch charge carried by the localized defect. |
| \(q_{\rm eff}\) | Observable brane charge after localization reduction. |
| \(e_*\) | Microscopic unit charge scale in the branch sector. |
| \(e_{\rm eff}\) | Observable electric charge magnitude after reduction. |
| \(Z(w)\) | Localization profile used in the reduced Maxwell sector. |
| \(Z_{\rm int}\) | Integrated localization weight \(\int Z(w)\,dw\). |
| \(\Gamma_{10}\) | Reduced tenth-power transverse overlap \(\Gamma_{10}=\int |\chi_0|^{10}dw\). |
| \(g_C\) | Effective Coulomb coefficient on the brane. In standard brane units \(g_C=e_{\rm eff}^2/(4\pi\epsilon_0)\). |
| \(S_{\rm leak}\) | Leakage / source term in the projected continuity law. |
| \(k_{1/2}\) | Lowest D/N half-shifted support wavenumber, \(k_{1/2}=\pi/(2L)\). |
| \(K_{22}\) | Wall/support coefficient in the \(22\)-channel; the reduced working value used in the healing note is \(K_{22}=8/3\). |

### A.3 Atomic and finite-size response quantities

| Symbol | Meaning |
|---|---|
| \(m_e, m_p\) | Electron and proton masses in the two-body upgrade. |
| \(M=m_e+m_p\) | Total two-body mass. |
| \(\mu = m_e m_p/(m_e+m_p)\) | Reduced mass. |
| \(E_\perp\) | Transverse-mode energy offset after separation. |
| \(\phi_a(r)\) | Hydrogenic \(1s\)-type variational ansatz \(\phi_a=(\pi a^3)^{-1/2}e^{-r/a}\). |
| \(T_2^{\rm eff}(r;a)\) | Finite-throat effective quadrupolar tidal load. |
| \(\mathcal F(r/a)\) | Core-resolved shape function appearing in \(T_2^{\rm eff}=(g_C/a^3)\mathcal F(r/a)\). |
| \(\delta V_{\rm fs}(r)\) | Finite-size correction to the reduced pair interaction. |

### A.4 Lepton, rotor, and anomaly quantities

| Symbol | Meaning |
|---|---|
| \(\alpha_{\rm top}\) | Topological Berry / half-flux parameter. Current isolated closure fixes \(\alpha_{\rm top}=1/2\). |
| \(\kappa_0\) | Berry coefficient of the internal rotor. |
| \(\Pi^{\rm mix}_{ij}\) | Integrated traceless planar stress generated by the mixed-sector core. |
| \(C_{\rm mix}\) | Positive mixed-sector stress coefficient in the \(P_{22}\) theorem. |
| \(h_\alpha\) | Intrinsic \(P_{22}\) bracing coefficient sourced by the mixed core. |
| \(k_{22},u_{22}\) | Quadratic and quartic mouth stiffness coefficients in the reduced \(P_{22}\) mouth energy. |
| \(I\) | Effective rotor inertia in the isolated lepton reduction. |
| \(M_\mu\) | Magnetic circulation inertia. |
| \(M_{\rm part}\) | Reduced particle inertia in the lepton package. |
| \(C_{2,q}\) | Charge-carrier second moment. |
| \(C_{2,M}\) | Inertial-carrier second moment. |
| \(\zeta\) | Charge/inertia mismatch factor, \(\zeta=C_{2,q}/C_{2,M}\). |
| \(\eta_1\) | Quadratic inertial backreaction coefficient in the anomaly law. |
| \(\mathcal B_\delta\) | Overlap softening factor generated by a blurred boundary visibility profile. |
| \(\chi_\delta(\rho)\) | Radial boundary visibility profile across the mouth. |
| \(f\) | Reduced anomaly expansion parameter \(f=\alpha_{\rm fs}/(2\pi)\). |

---

## Appendix B — Exact reduced equations used in the derivations

### B.1 Parent 4D matter sector

The current action stack uses the gauged 4D GNLS matter Lagrangian
\[
\mathcal L_\psi
=
\frac{i\hbar}2\left(\psi^*D_t\psi-\psi D_t\psi^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
-V_{\rm conf}(\mathbf X;a,L)|\psi|^2
-U(\rho),
\]
with stiff EOS potential
\[
U(\rho)=\frac K4\rho^5,
\qquad
h(\rho)=\frac{5K}4\rho^4.
\]

The corresponding reduced GNLS equation is
\[
i\hbar D_t\psi
=
\left[
-\frac{\hbar^2}{2m}D_iD_i
+V_{\rm conf}
+\frac{5K}4\rho^4
\right]\psi.
\]

### B.2 Projected continuity and leakage

The projected plasma-sector continuity identity carried throughout the later notes is
\[
\partial_t\rho_{\rm brane} + \partial_a j^a_{\rm brane} = S_{\rm leak}.
\]
This is the reduced source law behind the later healing-length and collar-transport notes.
In the successful local anomaly chain, the dark-energy route was dropped, but the **local** support / leakage structure remained essential.

### B.3 Localized Maxwell reduction and observable charge

In the static zero-mode regime the reduced Maxwell sector gives an effective Coulomb potential of the form
\[
A_0(r)
=
\frac{\mu_0^{\rm eff} q_*}{4\pi r}
\left[
1+\frac12 e^{-2r/\lambda}+\cdots
\right].
\]

The observable charge is localization-controlled:
\[
q_{\rm eff} = \frac{q_*}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff} = \frac{e_*}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int} = \int Z(w)\,dw.
\]
For a Gaussian localization profile one has
\[
Z_{\rm int} = \sqrt{\pi}\,\lambda.
\]

### B.4 Lowest-mode hydrogenic reduction

The hydrogenic reduction uses
\[
\psi(\mathbf x,w,t)=\chi_0(w)\,\phi(\mathbf x,t),
\qquad
\int |\chi_0|^2dw=1.
\]
This produces the reduced one-body brane energy
\[
E[\phi]
=
\int d^3x\,
\left[
\frac{\hbar^2}{2m}|\nabla\phi|^2
+\bigl(E_\perp + V_{\rm H}(r)\bigr)|\phi|^2
+\frac{K\Gamma_{10}}4|\phi|^{10}
\right],
\]
with
\[
\Gamma_{10} = \int |\chi_0(w)|^{10}dw
\]
and
\[
V_{\rm H}(r)
=
-\frac{g_C}r
-\frac{g_C}{2r}e^{-2r/\lambda}
+\cdots.
\]

### B.5 Two-body upgrade

The genuine pair reduction replaces the fixed-source approximation by
\[
\frac12 m_e\dot x_e^2 + \frac12 m_p\dot x_p^2
=
\frac12 M\dot R^2 + \frac12 \mu\dot r^2,
\qquad
M=m_e+m_p,
\qquad
\mu=\frac{m_e m_p}{m_e+m_p}.
\]

In the clean decoupling limit this gives the Bohr-scale minimum
\[
a_* = \frac{\hbar^2}{\mu g_C}
= \frac{4\pi\epsilon_0\hbar^2}{\mu e_{\rm eff}^2}.
\]

### B.6 Finite-size response and core regulation

The first successful atomic perturbation map showed that the active internal load is not the raw \(1/r\) Coulomb scalar but the finite-throat tidal / Hessian load.
At the reduced level this is summarized by
\[
T_2^{\rm eff}(r;a)=\frac{g_C}{a^3}\,\mathcal F(r/a),
\]
with asymptotic behavior
\[
\mathcal F(x)\sim -\frac3{x^3}
\quad (x\gg 1),
\qquad
\mathcal F(x)\sim -\sqrt{\frac2\pi}\frac{x^2}5
\quad (x\ll 1).
\]
The far-field divergence is therefore replaced by a finite core-resolved response near the throat waist.

### B.7 Mouth energy and rotor equations

The reduced isolated mouth energy is
\[
U_{\rm iso}(\sigma,\phi_m)
=
\frac12 k_{22}\sigma^2
+\frac14 u_{22}\sigma^4
-h_\alpha\sinh 2\sigma\cos 2(\phi_m-\phi_\alpha).
\]

The associated isolated rotor-plus-splitter Lagrangian is
\[
L_{\rm rot}^{\rm iso}
=
\frac I2\dot\varphi^2
+\hbar\alpha_{\rm top}\dot\varphi
-
V_\infty\cos 2(\varphi-\phi_\alpha),
\qquad
\alpha_{\rm top}=\frac12.
\]

### B.8 Common-carrier theorem and anomaly skeleton

The circulation-level magnetic relation is
\[
\mu_z
=
\frac{q_{\rm eff}}{2M_{\rm part}}
\frac{C_{2,q}}{C_{2,M}}
L_z.
\]
So
\[
\zeta \equiv \frac{C_{2,q}}{C_{2,M}},
\qquad
g_{\rm red} = 2\zeta.
\]

The common-carrier theorem is the pointwise closure
\[
dq = \frac{q_{\rm eff}}{M_{\rm part}}\,d\mu_{\rm eff},
\]
which implies
\[
C_{2,q}=C_{2,M},
\qquad
M_\mu=M_{\rm part},
\qquad
g_{\rm tree}=2.
\]

The anomaly program then organizes the mismatch factor as
\[
\zeta = 1 + f - (1+\eta_1)f^2 + O(f^3),
\qquad
f=\frac{\alpha_{\rm fs}}{2\pi}.
\]

---

## Appendix C — Variational integrals for the hydrogenic sector

### C.1 Trial state

The normalized \(1s\)-type trial state is
\[
\phi_a(r)=\frac1{\sqrt{\pi a^3}}e^{-r/a},
\qquad a>0,
\]
with
\[
\int d^3x\,|\phi_a|^2 = 1.
\]

### C.2 Radial integral identities

For later use,
\[
\int_0^\infty r^n e^{-\beta r}dr = \frac{n!}{\beta^{n+1}}
\qquad
(\Re\beta>0).
\]

Using spherical symmetry,
\[
\int d^3x\,F(r)=4\pi\int_0^\infty r^2F(r)\,dr.
\]

### C.3 Kinetic term

For the \(1s\) ansatz one finds
\[
|\nabla \phi_a|^2 = \frac1{a^2}|\phi_a|^2,
\]
hence
\[
\int d^3x\,|\nabla\phi_a|^2 = \frac1{a^2},
\qquad
\langle T\rangle = \frac{\hbar^2}{2ma^2}.
\]

### C.4 Coulomb expectation value

The Coulomb integral is
\[
\left\langle \frac1r \right\rangle
=
\int d^3x\,|\phi_a|^2\frac1r
=
\frac1a.
\]

### C.5 Yukawa-dressed Coulomb expectation value

For the first KK/Yukawa correction used in the hydrogen note,
\[
\left\langle \frac{e^{-2r/\lambda}}r \right\rangle
=
\int d^3x\,|\phi_a|^2\frac{e^{-2r/\lambda}}r
=
\frac{\lambda^2}{a(a+\lambda)^2}
=
\frac1{a(1+a/\lambda)^2}.
\]

### C.6 Inherited GNLS tenth-power moment

The stiff-EOS term requires
\[
\int d^3x\,|\phi_a|^{10}
=
\frac1{125\pi^4 a^{12}}.
\]

### C.7 Reduced one-parameter energy

Putting the preceding pieces together gives
\[
E(a)
=
E_\perp
+\frac{\hbar^2}{2ma^2}
-\frac{g_C}a
-\frac{g_C}{2a(1+a/\lambda)^2}
+\frac{K\Gamma_{10}}{500\pi^4 a^{12}}
+\cdots.
\]
In the clean hydrogenic decoupling limit this reduces to
\[
E(a)\approx E_\perp + \frac{\hbar^2}{2ma^2} - \frac{g_C}a,
\]
whose unique positive minimizer is
\[
a_* = \frac{\hbar^2}{m g_C}
\]
or, in the two-body upgrade,
\[
a_* = \frac{\hbar^2}{\mu g_C}.
\]

---

## Appendix D — \(P_{22}\) geometry and constant-area ellipse identities

### D.1 Area-preserving ellipse parametrization

The area-preserving \(P_{22}\) mouth ellipse is written as
\[
x(\Theta)=a e^\sigma \cos\Theta,
\qquad
y(\Theta)=a e^{-\sigma}\sin\Theta.
\]
Its semiaxes are
\[
a_x = a e^\sigma,
\qquad
a_y = a e^{-\sigma},
\]
so
\[
a_x a_y = a^2.
\]

### D.2 Constant area

The enclosed area is therefore exactly
\[
A = \pi a_x a_y = \pi a^2.
\]
This is the clean algebraic statement behind the constant-area throat constraint used in the atomic core regulator and in the lepton mouth mechanics.

### D.3 Magnetic area integral

The charged-loop area integral used in the Dirac bridge is
\[
\oint (x\,dy-y\,dx)=2\pi a^2.
\]
So the leading magnetic moment of the charged loop is **independent of the ellipticity** \(\sigma\).
This is why the \(P_{22}\) deformation protects the doublet without renormalizing the leading magnetic coefficient.

### D.4 Mouth quadrupole tensor

Define the headless director tensor
\[
\mathcal N_{ij}(\phi)=n_i n_j - \frac12\delta_{ij},
\qquad
n=(\cos\phi,\sin\phi).
\]
Then the mouth quadrupole tensor takes the reduced form
\[
Q^{\rm mouth}_{ij}(\sigma,\phi_m)
=
Q_m(\sigma)\,\mathcal N_{ij}(\phi_m),
\qquad
Q_m(\sigma)=\frac{a^2}2\sinh 2\sigma.
\]

### D.5 Contraction identity

The useful contraction identity is
\[
\mathcal N_{ij}(\phi_1)\mathcal N_{ij}(\phi_2)
=
\frac12\cos 2(\phi_1-\phi_2).
\]
This is the algebra that turns the scalar mouth-core coupling into the reduced
\[
-h_\alpha\sinh 2\sigma\cos 2(\phi_m-\phi_\alpha)
\]
form.

### D.6 Small-ellipticity expansion

For \(|\sigma|\ll 1\),
\[
\sinh 2\sigma = 2\sigma + \frac{(2\sigma)^3}6 + O(\sigma^5),
\]
so
\[
Q_m(\sigma)=a^2\sigma + O(\sigma^3).
\]
This is the regime used in the weak-bracing estimate
\[
\sigma_\infty \approx \frac{2h_\alpha}{k_{22}}.
\]

---

## Appendix E — Mixed-sector stress integrals and intrinsic \(P_{22}\) bracing

### E.1 Minimal mixed-core ansatz

The microscopic lepton bracing note used the representative mixed-sector core ansatz
\[
A_w(x,y,w)=b_i x_i\,\chi(r,w),
\qquad
r=\sqrt{x^2+y^2}.
\]
This is the minimal localized core profile that can carry a nontrivial planar traceless stress.

### E.2 Gauge-invariant traceless planar stress

The integrated mixed-core stress is
\[
\Pi^{\rm mix}_{ij}
\equiv
\int d^2x_\perp\,dw\;
\frac{Z(w)}{\mu_0}
\left(
\partial_iA_w\partial_jA_w
-\frac12\delta_{ij}\partial_kA_w\partial_kA_w
\right).
\]

For the ansatz above, exact angular averaging gives
\[
\Pi^{\rm mix}_{ij}
=
C_{\rm mix}
\left(
b_i b_j - \frac12 s\,\delta_{ij}
\right),
\qquad
s=b_k b_k,
\]
with
\[
C_{\rm mix}
=
\frac{\pi}{2\mu_0}
\int_{-\infty}^\infty dw\,Z(w)
\int_0^\infty r\,dr\,
\bigl(2\chi+r\partial_r\chi\bigr)^2.
\]
Because the integrand is a square and \(Z(w)>0\) in the localized sector,
\[
C_{\rm mix}>0
\]
for every nontrivial localized core profile.

### E.3 Director form

Using the headless tensor \(\mathcal N_{ij}\),
\[
\Pi^{\rm mix}_{ij}
=
C_{\rm mix}\,s\,\mathcal N_{ij}(\phi_\alpha).
\]
So the mixed-sector core automatically produces a true \(P_{22}\) planar stress rather than a scalar or vector bracer.

### E.4 Topological half-flux lock

The reduced mixed-sector dynamics use
\[
\alpha_{\rm top} = \frac{\nu_0 s}{2\hbar}.
\]
With the isolated closure
\[
\alpha_{\rm top}=\frac12,
\]
the baseline amplitude becomes
\[
s_\alpha = \frac{\hbar}{\nu_0}.
\]
Substituting this into the stress theorem gives
\[
\Pi^{(\alpha)}_{ij}
=
C_{\rm mix}\frac{\hbar}{\nu_0}\,\mathcal N_{ij}(\phi_\alpha).
\]

### E.5 Mouth coupling and intrinsic bracing coefficient

The lowest symmetry-allowed scalar mouth-core coupling is
\[
U_{\rm coup}
=
-\Lambda_Q\,\Pi^{(\alpha)}_{ij}Q^{\rm mouth}_{ij}.
\]
Using the contraction identity from Appendix D, this becomes
\[
U_{\rm coup}
=
-h_\alpha\sinh 2\sigma\cos 2(\phi_m-\phi_\alpha),
\]
with
\[
h_\alpha
=
\frac{\Lambda_Q a^2 C_{\rm mix}\hbar}{4\nu_0}.
\]

### E.6 Weak-bracing isolated ellipticity

The isolated mouth energy then yields
\[
\sigma_\infty
=
\frac{2h_\alpha}{k_{22}} + O(h_\alpha^2)
=
\frac{\Lambda_Q a^2 C_{\rm mix}\hbar}{2k_{22}\nu_0}+O(h_\alpha^2).
\]
So the isolated throat remains in the \(P_{22}\) sector whenever the mixed profile is nontrivial and \(\Lambda_Q\neq 0\).

### E.7 Gaussian worked example

For the separable Gaussian representative
\[
\chi(r,w)=\frac1a e^{-r^2/(2a^2)}W_1(w),
\]
one finds
\[
C_{\rm mix}^{\rm G}
=
\frac{\pi}{2\mu_0}
\int_{-\infty}^\infty Z(w)W_1(w)^2\,dw.
\]
With the simple odd-mode choice
\[
Z(w)=e^{-w^2/\lambda^2},
\qquad
W_1(w)=\frac w\lambda e^{-w^2/(2\lambda^2)},
\]
this evaluates to
\[
C_{\rm mix}^{\rm G}
=
\frac{\pi\sqrt{2\pi}}{16\mu_0}\lambda,
\]
and therefore
\[
h_\alpha^{\rm G}
=
\frac{\Lambda_Q\pi\sqrt{2\pi}}{64\mu_0}
\frac{a^2\lambda\hbar}{\nu_0}.
\]

---

## Appendix F — Local transport PDEs and visibility profiles

### F.1 Inertial collar transport PDE

The successful local blur closure used the stationary collar transport equation
\[
\partial_t\chi + v_{\rm mix}\partial_x\chi = \omega_{1/2}(1-\chi),
\qquad
x=a(1-\rho).
\]
The minimal local closure takes
\[
v_{\rm mix} = f c_s,
\qquad
\omega_{1/2} = c_s k_{1/2},
\qquad
k_{1/2}=\frac{\pi}{2L}.
\]

Its stationary solution is
\[
\chi(x)=1-e^{-x/\ell_{\rm blur}},
\qquad
\ell_{\rm blur}=\frac{v_{\rm mix}}{\omega_{1/2}}
=\frac f{k_{1/2}}
=\frac{2L}\pi f.
\]
Hence the local reduced blur thickness is
\[
\delta_{\rm loc}
=
\frac{\ell_{\rm blur}}a
=
\frac{2L}{\pi a}f
=
\kappa f.
\]

### F.2 Exponential visibility profile

The corresponding radial visibility profile is
\[
\chi_{\exp}(\rho)=1-e^{-(1-\rho)/\delta}.
\]
The overlap factor used to soften the inertial \(22\)-self-loop is
\[
\mathcal B_\delta
=
6\int_0^1 \rho^5\chi_\delta(\rho)\,d\rho.
\]
For the exponential profile,
\[
\mathcal B_{\exp}(\delta)
=
1-6\delta+30\delta^2-120\delta^3+360\delta^4-720\delta^5
+720\delta^6\bigl(1-e^{-1/\delta}\bigr).
\]

The blurred inertial coefficient is then
\[
\eta_1(f)
=
\frac{11}{36}\,\mathcal B_{\exp}(\kappa f).
\]

### F.3 Linear collar profile used as an earlier geometric check

Before the local PDE closure, the linear collar profile
\[
\chi_{\rm lin}(\rho)
=
\begin{cases}
1, & 0\le \rho \le 1-\delta,\\
\dfrac{1-\rho}\delta, & 1-\delta \le \rho \le 1
\end{cases}
\]
gave
\[
\mathcal B_{\rm lin}(\delta)
=
1-3\delta+5\delta^2-5\delta^3+3\delta^4-\delta^5+\frac{\delta^6}7.
\]
This profile is retained here because it gives a useful check that the needed electron-scale softening is geometrically tiny and profile-stable.

### F.4 Charge-side collar mode

The final successful anomaly closure added a mean-zero local charge-shape mode inside the already-transferred collar:
\[
\phi_q(\rho)
=
\cos\!\left(\frac{\pi(1-\rho)}{2\tau}\right)-\bar c(f),
\qquad
\tau(f)=1-\sqrt{1-f},
\]
with collar average
\[
\langle Y\rangle_{\rm col}
=
\frac1f \int_{\sqrt{1-f}}^1 2\rho Y(\rho)\,d\rho.
\]
The mean-zero condition is
\[
\bar c(f)=\left\langle
\cos\!\left(\frac{\pi(1-\rho)}{2\tau}\right)
\right\rangle_{\rm col},
\qquad
\langle \phi_q\rangle_{\rm col}=0.
\]

The minimal amplitude uses the same D/N support scale as the inertia-side blur:
\[
A_q(f)=\frac{a\tau}{2L/\pi}=\frac\tau\kappa.
\]

### F.5 Local charge second moment

Define
\[
\Xi(f)=\langle \rho^2\phi_q\rangle_{\rm col}.
\]
Then the global charge second-moment ratio becomes
\[
Q_{\rm loc}(f)
=
1+f-f^2+2fA_q(f)\Xi(f).
\]
Its leading small-\(f\) behavior is
\[
\Xi(f)=\left(\frac4{\pi^2}-\frac1\pi\right)f+O(f^2),
\qquad
A_q(f)=\frac f{2\kappa}+O(f^2),
\]
so
\[
Q_{\rm loc}(f)
=
1+f-f^2
+
\frac{4-\pi}{\pi^2\kappa}f^3
+O(f^4).
\]

---

## Appendix G — Series coefficients for the reduced anomaly law

This appendix collects the anomaly coefficients in the order they appeared.

### G.1 Tree level

The common-carrier theorem gives
\[
\zeta_{\rm tree}=1,
\qquad
g_{\rm tree}=2.
\]

### G.2 Schwinger-scale support transfer

The first successful geometric support-transfer layer gives the natural one-loop scale
\[
f=\frac{\alpha_{\rm fs}}{2\pi}.
\]
At this level the reduced anomaly law begins as
\[
g = 2\bigl[1+f+O(f^2)\bigr].
\]

### G.3 Sharp self-loop closure

The canonically normalized \(22\)-self-loop yields
\[
\eta_1^{\rm min} = \frac{11}{36},
\]
so the sharp-boundary reduced law is
\[
g_{\rm sharp}
=
2\left[
1+f-\frac{47}{36}f^2+O(f^3)
\right].
\]

### G.4 Local inertia blur

With the local PDE blur \(\delta_{\rm loc}=\kappa f\), one has
\[
\eta_1(f)=\frac{11}{36}\mathcal B_{\exp}(\kappa f),
\]
and the exact reduced form
\[
g_{\rm blur}
=
2\Bigl[
1+f-(1+\eta_1(f))f^2
\Bigr].
\]

Expanding the exponential overlap gives
\[
g_{\rm blur}
=
2\left[
1+f-\frac{47}{36}f^2
+\frac{11}6\kappa f^3
-\frac{55}6\kappa^2 f^4
+O(f^5)
\right].
\]

### G.5 Local charge-side transport closure

After adding the mean-zero local charge-shape mode,
\[
g_{\rm loc}(f)
=
2\Bigl[Q_{\rm loc}(f)-\eta_1(f)f^2\Bigr].
\]
Its cubic expansion is
\[
g_{\rm loc}(f)
=
2\left[
1+f-\frac{47}{36}f^2
+
\left(
\frac{11}6\kappa + \frac{4-\pi}{\pi^2\kappa}
\right)f^3
+O(f^4)
\right].
\]

The cubic coefficients are therefore
\[
c_{3,\rm inertia} = \frac{11}6\kappa,
\qquad
c_{3,q} = \frac{4-\pi}{\pi^2\kappa},
\qquad
c_{3,\rm total}=c_{3,\rm inertia}+c_{3,q}.
\]

### G.6 Residual structure

At the current best reduced closure, the remaining miss sits naturally at
\[
O(f^4),
\]
which is exactly where a coupled common charge-inertia transport PDE should first become unavoidable.
So the present anomaly hierarchy is internally consistent in the sense that the remaining discrepancy is of the same order as the next omitted transport layer.

---

## Appendix H — Numerical values and benchmark outputs

The tables below collect the numerical inputs reused across the successful derivation chain.
They are not a substitute for rerunning the SymPy helpers, but they are the frozen working values used in the present write-up.

### H.1 External numerical inputs

| Quantity | Value |
|---|---:|
| Fine-structure constant \(\alpha_{\rm fs}\) | 0.0072973525643 |
| Reduced anomaly parameter \(f=\alpha_{\rm fs}/(2\pi)\) | 0.001161409732093 |
| Electron \(|g_e|\) target used in the write-up | 2.00231930436092 |
| Electron anomaly \(a_e=|g_e|/2-1\) | 0.00115965218046 |
| Frozen throat aspect ratio \(L/a\) | 1.85 |
| Reduced support coefficient \(\kappa=2L/(\pi a)\) | 1.177746578880 |
| Wall/support coefficient \(K_{22}\) used in the healing note | \(8/3\) |

### H.2 Reduced coefficients

| Quantity | Value |
|---|---:|
| Minimal self-loop coefficient \(\eta_1^{\rm min}=11/36\) | 0.305555555555556 |
| Electron-target effective \(\eta_1^{(e)}\) | 0.302978262909392 |
| Local blur \(\delta_{\rm loc}=\kappa f\) | 0.001367846338650 |
| Exponential overlap \(\mathcal B_{\exp}(\delta_{\rm loc})\) | 0.991848746223625 |
| Local inertial coefficient \(\eta_1^{\rm loc}=(11/36)\mathcal B_{\exp}\) | 0.303064894679441 |
| Charge-side cubic coefficient \(c_{3,q}\) | 0.073848525604098 |
| Total cubic coefficient \(c_{3,\rm total}\) | 2.233050586884145 |

### H.3 Benchmark anomaly predictions

| Closure level | Reduced prediction |
|---|---:|
| Sharp \(11/36\) self-loop | 2.00231929740804 |
| Local inertia blur only | 2.00231930412721 |
| Final local charge-transport closure | 2.00231930435865 |
| Remaining miss \(|g_e|-g_{\rm loc}\) | 2.270e-12 |

### H.4 Interpretation of the benchmark ladder

The benchmark sequence shows the internal logic of the successful anomaly program.

1. **Tree level** closes at \(g=2\).
2. **Support transfer** supplies the Schwinger-scale one-loop magnitude.
3. **Quadratic inertial backreaction** fixes the leading negative correction through \(\eta_1=11/36\).
4. **Local inertia blur** converts the sharp mouth into a collar transport problem.
5. **Local charge-shape transport** removes essentially all of the remaining visible residual without introducing a new fit parameter.

The final remaining discrepancy is therefore no longer of a size that suggests a new mechanism.
It is naturally of the order expected from the next unresolved transport layer.
