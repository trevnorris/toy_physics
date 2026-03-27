# Lepton / Electron Program — Sections 1–3 Draft

## 1. Frozen background and source anchors

The work so far sits inside a deliberately frozen reduced hierarchy rather than a completed first-principles moving-throat solution. The carry-forward assumptions are:

- EOS exponent fixed at \(n=5\).
- Added-inertia coefficient fixed at \(\kappa_{\rm add}=1/2\).
- Pressure-volume response coefficient fixed at \(\kappa_{\rm PV}=3/2\) within the declared adiabatic one-degree-of-freedom closure.
- Conservative 1PN bookkeeping fixed with \(\beta_{\rm 1PN}=3\) and exact conservative EIH equality in the solved scope.
- Preferred cavity/throat aspect ratio selected at \(L/a\approx 1.85\).
- The 2PN conservative particle/throat structure sharpened into a carried-forward dipole wake sector plus a genuine \(P_0\oplus P_2\) mouth/support response layer and a separate geometry-closure channel.

Those statements matter because everything in the mass program below is conditional on them. We are not claiming a universal theorem for every possible protocol of the parent model. We are claiming theorems inside the present reduction hierarchy.

Just as important are the explicit non-claims. The current stack does **not** yet solve the full moving-throat bulk PDE. It also leaves unresolved the dissipative/open-system sectors, radiation reaction, strong-field completion, spin couplings, and higher-PN extensions beyond the solved conservative regime. In the 1PN full summary, even some background objects remain symbolic, including the EOS normalization \(K_{\rm EOS}\) and the detailed localization profiles such as \(W(w)\) and \(Z(w)\). So the derivations below should be read as reduced-closure theorems, not final completed defect simulations.

At the ontology level, the model still treats a particle as a brane–bulk throat defect. The brane sees a roughly monopolar mouth field, the interior carries throat/cavity structure, and charge and mass are emergent book-keeping layers after reduction rather than primitive ingredients. Earlier carry-forward dictionaries already encoded this perspective through geometric quantities like
\[
q\sim \kappa_q\rho_0\pi a^2\Gamma,
\qquad
m_G\sim \kappa_m\rho_0\pi a^2L.
\]
Those are not the final particle-identity dictionary, but they are the right source anchors for the present program.

The bulk parent theory already contains the energy ledgers needed to motivate a mass-from-defect-energy approach. In the 4D GNLS/Maxwell stack, the gauge-invariant matter energy density is of the form
\[
u_s=
\frac{\hbar^2}{2m_s}(D_A\psi_s)^*(D_A\psi_s)+V_{{\rm conf},s}\rho_s+U_s(\rho_s),
\]
while the localized electromagnetic energy density is
\[
u_{\rm EM}=
\frac{Z(w)}{2\mu_0}\left(E_AE_A+\tfrac12F_{AB}F_{AB}\right),
\]
and the total bulk energy is the integral of matter plus field contributions over the three brane coordinates and the transverse \(w\) direction. That is the parent ledger from which the reduced defect-mass program was abstracted.

For carry-forward work, the source anchors worth keeping visible are:

- `lepton_target_handoff.md` for the frozen hierarchy and methodology;
- `4d_1pn_full_summary.md` for what is fixed versus what remains symbolic;
- `4d_plasma.tex` / `4d_plasma_summary.md` for the parent bulk energy ledger;
- `4d_2pn_summary.md` and the 2PN notes for the \(P_0\oplus P_2\) mouth/support structure;
- `atom_work.md` and `lepton_work.md` for the later lepton corridor and magnetic-moment targets.

## 2. Core isolated-electron mass theorem

The first durable mass result is that, within the present closure, the isolated defect can be represented by a reduced one-parameter static energy functional of the form
\[
F(a,\rho)=\frac{A(\rho)}{a}+\frac{B(\rho)}{a^2}+C(\rho)a^3.
\]
Here the three terms are interpreted as:

- a trapped support-wave ledger \(E_w=A/a\),
- a self-flow or feed ledger \(E_f=B/a^2\),
- a pressure-volume or compressional ledger \(E_{\rm PV}=Ca^3\).

This is the cleanest reduced rest-energy ansatz that matches the scalings already fixed by the bridge hierarchy. It is not the full 4D PDE. It is the reduced defect energy ledger consistent with the current carry-forward closure.

At fixed ambient density \(\rho\), stationarity with respect to the throat size \(a\) gives
\[
\partial_aF=0
\quad\Longrightarrow\quad
-\frac{A}{a^2}-\frac{2B}{a^3}+3Ca^2=0.
\]
Multiplying by \(a\) and rewriting in terms of the three energy ledgers gives the virial identity
\[
-E_w-2E_f+3E_{\rm PV}=0,
\]
or equivalently
\[
\boxed{E_w+2E_f=3E_{\rm PV}.}
\]
This is the first theorem worth carrying forward: the reduced defect mass ledger is not arbitrary. At equilibrium, its three sectors are tightly locked.

The second step is to combine that virial law with the already frozen 1PN density response target. Define
\[
x\equiv \frac{E_f}{E_w}.
\]
Then the virial identity implies
\[
\frac{E_{\rm PV}}{E_w}=\frac{1+2x}{3},
\qquad
F=E_w\frac{4+5x}{3}.
\]
Inside the adiabatic one-degree-of-freedom closure, the bridge program fixes \(n=5\), \(\kappa_\rho=1\), and \(\kappa_{\rm PV}=3/2\). With those inputs, the density-slope ledger forces
\[
\frac{E_f}{E_w}=x=\frac{2}{11}.
\]
Once that ratio is known, the whole internal partition is fixed:
\[
\frac{E_f}{E_w}=\frac{2}{11},
\qquad
\frac{E_{\rm PV}}{E_w}=\frac{5}{11},
\]
and therefore
\[
\boxed{E_w:E_f:E_{\rm PV}=11:2:5.}
\]
This is one of the strongest reduced results of the entire program. It means the present hierarchy already fixes the **shape** of the isolated rest-mass ledger.

A second direct consequence is the total-to-wave relation. If we write
\[
E_w=11u,\qquad E_f=2u,\qquad E_{\rm PV}=5u,
\]
then the total stationary rest energy is
\[
F_*=18u=\frac{18}{11}E_w.
\]
So the zeroth-order isolated particle mass is not an extra input: it is exactly \(18/11\) times the trapped support-wave ledger on the same branch.

The same result can be written directly in terms of the reduced coefficients. Imposing the partition on
\[
F(a)=\frac{A}{a}+\frac{B}{a^2}+Ca^3
\]
gives
\[
a_*=\frac{11B}{2A},
\qquad
C=\frac{80A^5}{11^5B^4},
\]
and therefore
\[
E_w^*=\frac{2A^2}{11B},
\qquad
E_f^*=\frac{4A^2}{121B},
\qquad
E_{\rm PV}^*=\frac{10A^2}{121B}.
\]
Adding them gives the exact stationary reduced mass coefficient
\[
\boxed{F_*=\frac{36A^2}{121B}.}
\]
This is the core isolated-electron mass theorem in reduced form.

One more identity from the same closure is worth keeping. Differentiating the stationary condition with respect to ambient density yields the predicted breathing slope
\[
\boxed{\frac{d\ln a}{d\ln \rho}=-\frac{57}{64}.}
\]
So the current hierarchy fixes not only the internal partition, but also how the equilibrium throat radius responds to changes in the background density.

What Section 2 does **not** do is close the absolute mass scale. The theorem fixes the reduced structure,
\[
F_*\propto \frac{A^2}{B},
\]
but it does not by itself determine the microscopic value of \(A^2/B\). That is exactly why the next section matters.

## 3. Microscopic identification of the reduced coefficients

The main gain of the next stage was to stop treating \(A\) and \(B\) as abstract reduction coefficients. Within a minimal isolated-throat ansatz, both can be identified microscopically.

### 3.1 The support-wave coefficient

Assume the isolated throat carries a trapped support mode with frequency
\[
\omega_{\rm supp}=c_s k_{\rm supp},
\qquad
k_{\rm supp}=\frac{\chi_w}{a},
\]
where \(\chi_w\) is a pure dimensionless mode number fixed by the chosen support branch. Let the adiabatic support action be
\[
\mathcal I_w\equiv \frac{E_w}{\omega_{\rm supp}}.
\]
Then the support-wave energy is
\[
E_w=\mathcal I_w\omega_{\rm supp}
=\mathcal I_w c_s\frac{\chi_w}{a}.
\]
Comparing with the reduced ledger \(E_w=A/a\) immediately gives
\[
\boxed{A(\rho)=\mathcal I_w\chi_w c_s(\rho).}
\]
So the wave coefficient is not mysterious: it is the support action times the relevant support-mode wavenumber coefficient times the local sound speed.

On the lowest mixed D/N branch, the mode number is
\[
\chi_w=\frac{\pi}{2\Lambda},
\qquad \Lambda\equiv \frac{L}{a},
\]
so the coefficient specializes to
\[
\boxed{A(\rho)=\frac{\pi\mathcal I_w}{2\Lambda}c_s(\rho).}
\]
That is the first microscopic closure of the reduced mass theorem.

### 3.2 The self-flow coefficient

Now take the feed/self-flow ledger to come from a compact isotropic radial flux in four spatial dimensions. If \(\Phi\) is the conserved radial throughput, then in \(d\) spatial dimensions the static radial speed scales as
\[
u(r)=\frac{\Phi}{\rho S_{d-1}r^{d-1}},
\]
where \(S_{d-1}\) is the area of the unit \((d-1)\)-sphere. The static kinetic energy is then
\[
E_f^{(d)}=\frac{\rho}{2}\int_a^\infty \nu(r)^2 S_{d-1}r^{d-1}\,dr
=\frac{\Phi^2}{2S_{d-1}(d-2)\rho a^{d-2}}.
\]
The sharp selection result is that:
\[
E_f^{(3)}=\frac{\Phi^2}{8\pi\rho a},
\qquad
E_f^{(4)}=\frac{\Phi^2}{8\pi^2\rho a^2}.
\]
Only the compact **4D** self-flow reproduces the bridge scaling \(1/(\rho a^2)\). So the current reduced hierarchy is selecting a 4D compact self-flow branch, not a 3D slice-flow branch.

Comparing the 4D result with the reduced ledger \(E_f=B/a^2\) gives
\[
\boxed{B(\rho)=\frac{\Phi^2}{8\pi^2\rho}.}
\]
That is the second microscopic closure.

### 3.3 The support-quantum form of the isolated mass

Once the microscopic \(A\) and \(B\) are substituted back into the stage-2 reduced mass theorem, the reference-branch rest energy becomes
\[
F_0=\frac{288\pi^2\mathcal I_w^2\chi_w^2\rho_0c_{s0}^2}{121\Phi^2}.
\]
But the local reference branch also obeys the fixed partition \(E_f/E_w=2/11\), which gives
\[
a_0=\frac{11\Phi^2}{16\pi^2\mathcal I_w c_{s0}\chi_w\rho_0}.
\]
Eliminating \(\Phi\) between those two expressions collapses the whole result to
\[
\boxed{F_0=\frac{18}{11}\frac{\mathcal I_w\chi_w c_{s0}}{a_0}.}
\]
Since
\[
E_{w,0}=\frac{\mathcal I_w\chi_w c_{s0}}{a_0},
\]
this is simply
\[
\boxed{F_0=\frac{18}{11}E_{w,0}.}
\]
So the support-wave interpretation of the reduced theorem is exact: the total isolated rest energy is \(18/11\) times the trapped support-wave energy on that branch.

Specializing again to the D/N half-mode,
\[
\chi_w=\frac{\pi}{2\Lambda},
\qquad
\omega_{1/2}=\frac{\pi c_s}{2L},
\]
we get
\[
\boxed{F_0=\frac{18}{11}\mathcal I_w\omega_{1/2}.}
\]
If one later imposes the minimal support-action quantization
\[
\mathcal I_w=\hbar,
\]
then the reference-branch theorem becomes
\[
\boxed{F_0=\frac{18}{11}\hbar\omega_{1/2}.}
\]

That is as far as the present reduced hierarchy goes cleanly. It does **not** yet derive the absolute electron mass, because the absolute branch scale \(a_0\) or equivalently \(L_0\) is still unresolved. But it does isolate the last obstruction very precisely: the program no longer lacks a microscopic energy functional, and it no longer lacks a clean reference-branch mass theorem. What it still lacks is the final scale-fixing closure for the isolated defect.

# Lepton / Electron Program — Sections 4–6 Draft

## 4. Exact D/N spectral closure

The next durable step was to remove the last ambiguity in the support channel. Earlier stages had already shown that the reduced rest-energy theorem can be written as
\[
F_0=\frac{18}{11}\,\mathcal I_w\,\omega_{\rm supp},
\]
but the “half-mode” support clock was still being used as a plausible microscopic branch rather than an exact throat eigenmode. The D/N closure closes that gap.

### 4.1 Exact finite-throat D/N support ladder

On the exact finite-throat branch, the longitudinal support field lives on the interval
\[
z\in[0,L]
\]
and obeys
\[
\psi''+k^2\psi=0.
\]
The mouth is Dirichlet,
\[
\psi(0)=0,
\]
and the bottom is Neumann,
\[
\psi'(L)=0.
\]

Writing
\[
\psi(z)=A\sin(kz)+B\cos(kz),
\]
the mouth condition sets \(B=0\), so
\[
\psi(z)=A\sin(kz).
\]
Then the bottom condition requires
\[
Ak\cos(kL)=0.
\]
For a nontrivial mode \(A\neq0\), one must have
\[
\cos(kL)=0.
\]
Therefore
\[
\boxed{
k_j=\frac{\pi}{L}\left(j+\frac12\right),
\qquad
j=0,1,2,\dots
}
\]
and
\[
\boxed{
\omega_j=c_s k_j
=\frac{\pi c_s}{L}\left(j+\frac12\right).
}
\]

So the half-shifted support branch is not an ansatz. It is the exact D/N free-eigenmode spectrum of the finite throat.

This same branch is visible directly in the exact mouth operator. Solving the boundary-value problem with a prescribed mouth datum \(\psi_m\) gives
\[
\psi(z,\omega)=\psi_m\frac{\cos(k(L-z))}{\cos(kL)},
\qquad
k=\frac{\omega}{c_s},
\]
so the mouth derivative defines
\[
\boxed{
Z_{00}(\omega)=-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L}{c_s}\right).
}
\]
Its poles occur at
\[
\boxed{
\omega_j^{\rm pole}=\frac{\pi c_s}{L}\left(j+\frac12\right),
}
\]
exactly the same ladder.

### 4.2 Exact round-trip phase closure

The later same-charge analysis showed that the lepton corridor depends not only on the existence of the support spectrum, but on whether the common scalar loop phase is actually forced to cancel.

Write the scalar round-trip factor as
\[
R_{\rm rt}=r_0\,r_L\,e^{2ikL},
\qquad
\phi_0=\delta_0+\delta_L+2kL.
\]
For the D/N branch, the reflection factors are fixed by the boundary conditions themselves:
\[
r_D=-1,
\qquad
r_N=+1.
\]
So the full round-trip multiplier is
\[
R_{\rm rt}=r_Dr_N e^{2ikL}=(-1)(+1)e^{2ikL}.
\]
On the exact D/N ladder,
\[
2k_jL=(2j+1)\pi,
\]
hence
\[
R_{\rm rt}(k_j)=(-1)e^{i(2j+1)\pi}=1.
\]
Therefore every trapped D/N mode obeys
\[
\boxed{
R_{\rm rt}=1,
\qquad
\phi_0\equiv0\pmod{2\pi}.
}
\]

This result is stronger than a mere “steady state.” A pump-balanced state can keep the amplitude constant without forcing the intrinsic round-trip coefficient to equal one. The cycle-map form
\[
A_{n+1}=\Lambda A_n + D
\]
makes that explicit. Only the autonomous eigenmode closure
\[
D=0,
\qquad
A_{n+1}=A_n\neq0
\]
forces
\[
\Lambda=1.
\]
The exact D/N trapped branch realizes exactly that self-reproducing scalar support law.

A useful cross-check is that the D/N eigenfunctions
\[
\psi_j(z,t)=A_j\sin(k_jz)e^{-i\omega_j t}
\]
carry zero longitudinal standing-wave current:
\[
J_z\propto \Im(\psi_j^*\partial_z\psi_j)=0.
\]
So this branch is not a one-way transport state. It is a genuine trapped support channel.

### 4.3 Fixed-geometry tower versus self-consistent particle tower

Once the D/N ladder is exact, two different excitation questions must be separated.

If the geometry is externally held fixed, then the stage-2 rest-energy theorem gives
\[
F_j^{\rm(fixed)}
=
\frac{18}{11}\mathcal I_w\omega_j
=
\frac{9\pi}{11L}(2j+1)\,\mathcal I_wc_s.
\]
With \(L=\Lambda a\),
\[
\boxed{
F_j^{\rm(fixed)}
=
\frac{9\pi}{11\Lambda}(2j+1)\frac{\mathcal I_wc_s}{a}.
}
\]
This is the correct ladder for support excitations on one and the same throat.

But isolated particles do not live at fixed geometry. The stage-1 stationary theorem says
\[
F_*=\frac{36}{121}\frac{A^2}{B},
\]
and on the D/N branch
\[
A_j=\mathcal I_wc_s\chi_j,
\qquad
\chi_j=\frac{\pi}{\Lambda}\left(j+\frac12\right).
\]
So the stationary radius itself responds to the chosen support mode:
\[
a_j=\frac{11B}{2A_j}
=
\frac{11\Lambda\Phi^2}{8\pi^3\mathcal I_wc_s\rho\,(2j+1)}.
\]
Thus
\[
a_j\propto \frac{1}{2j+1}.
\]
Once the throat re-equilibrates, the self-consistent support frequency scales as
\[
\omega_j^{\rm(eq)}\propto (2j+1)^2,
\]
and so does the rest energy:
\[
\boxed{
F_j^{\rm(eq)}
=
\frac{72\pi^4\mathcal I_w^2c_s^2\rho}{121\Lambda^2\Phi^2}(2j+1)^2.
}
\]
Therefore the self-consistent isolated-particle tower is
\[
\boxed{
\frac{F_j^{\rm(eq)}}{F_0^{\rm(eq)}}=(2j+1)^2.
}
\]

This is one of the sharpest falsifiers in the whole chain. If the only family label were the D/N longitudinal support index \(j\), then the first two excited branches would predict
\[
1:9:25:\dots
\]
for the charged-lepton mass ladder. But the observed muon/electron and tau/electron ratios are about
\[
206.77
\qquad\text{and}\qquad
3477.37,
\]
so the naïve support-only family picture is decisively ruled out.

### 4.4 What Section 4 really establishes

Section 4 closes three things worth remembering.

First, the half-shifted support branch is exact, not heuristic.

Second, the scalar loop phase on the trapped D/N branch really does cancel:
\[
\phi_0\equiv0\pmod{2\pi}.
\]

Third, the simplest family picture already fails. The isolated support tower can organize an electron-like ground branch, but it cannot by itself generate the charged-lepton mass hierarchy.

That is why the next two sections had to split apart the remaining mass bottleneck from the already-closed phase structure.

## 5. Charge / circulation / throughput audit

After the D/N spectral closure, the isolated-electron rest-energy theorem was already very compact:
\[
F_0^{\rm(eq)}
=
\frac{72\pi^4}{121\Lambda^2}
\frac{\mathcal I_w^2c_s^2\rho}{\Phi^2},
\qquad
\Lambda=\frac{L}{a}\approx1.85.
\]
So the remaining obstacle was precise. Everything left was packed into the single amplitude combination
\[
\frac{\mathcal I_w^2c_s^2\rho}{\Phi^2}.
\]

The natural next question was whether \(\Phi\) could be fixed topologically, by identifying it with some charge/circulation invariant. The audit showed that this is **not** what the corrected ontology gives.

### 5.1 What the exact topology law actually fixes

In the corrected 4D GNLS ontology, the gauge-invariant bulk velocity is
\[
v_i=\frac{\hbar}{m_\psi}\left(\partial_i\theta-\frac{q_*}{\hbar}A_i\right),
\]
where \(q_*=\eta_Qe_*\) is the signed microscopic branch label and \(m_\psi\) is the GNLS matter parameter.

For any loop \(\mathcal C\) surrounding the mouth, phase single-valuedness gives the exact fluxoid law
\[
\oint_{\mathcal C}\left(\partial_i\theta-\frac{q_*}{\hbar}A_i\right)d\ell^i=2\pi n,
\qquad
n\in\mathbb Z.
\]
Equivalently,
\[
\Gamma_{\mathcal C}
\equiv
\oint_{\mathcal C}v_i\,d\ell^i
=
n\frac{h}{m_\psi}
-\frac{q_*}{m_\psi}\Phi_{\mathcal C},
\]
with \(\Phi_{\mathcal C}\) the enclosed gauge flux.

This is a rigid topological law, but it is a law for the **circulation/fluxoid sector**. It quantizes the tangential vortical class. It does not yet determine the radial feed amplitude \(\Phi\).

That distinction matters because the corrected ontology keeps several objects separate:

- electric sign is the puncture orientation \(\eta_Q=\pm1\),
- microscopic charge branch is \(q_*=\eta_Qe_*\),
- observable brane charge is thickness-dressed,
  \[
  q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
  \]
- circulation belongs to the magnetic/vortical sector,
- and radial throughput is a separate 4D flux amplitude.

So “charge,” “circulation,” and “throughput” are not interchangeable labels for one quantity.

### 5.2 What exact continuity fixes

The parent bulk continuity law is exact:
\[
\partial_t\rho+\partial_i j^i=0,
\qquad
j^i=\rho v^i.
\]
For a stationary isolated defect,
\[
\partial_t\rho=0,
\qquad
\nabla_4\cdot\mathbf j=0.
\]

Split the current into a compact 4D radial part and a tangential circulation part:
\[
\mathbf j=j_r(r)\hat{\mathbf r}+\mathbf j_{\rm circ}.
\]
The tangential circulation is divergence-free by itself, so the radial equation becomes
\[
\nabla_4\cdot\bigl(j_r(r)\hat{\mathbf r}\bigr)=0.
\]
In four spatial dimensions,
\[
\nabla_4\cdot\bigl(f(r)\hat{\mathbf r}\bigr)
=
\frac1{r^3}\frac{d}{dr}\bigl(r^3f(r)\bigr).
\]
Hence
\[
\frac{d}{dr}\bigl(r^3j_r(r)\bigr)=0
\]
and therefore
\[
j_r(r)=\frac{\Phi}{2\pi^2r^3},
\qquad
u_r(r)=\frac{\Phi}{2\pi^2\rho\,r^3}.
\]

This is exact and useful, but its meaning is very specific:
\[
\Phi=\int_{S_r^3}\mathbf j\cdot d\mathbf S
\]
is an integration constant of stationary continuity. Continuity proves that it is conserved. It does **not** quantize it.

### 5.3 The theorem of the audit

The sharp result of the audit is therefore
\[
\boxed{
\text{Topology fixes the circulation/holonomy class, but not the throughput amplitude } \Phi.
}
\]

This can be restated in geometric language. The current exact equations give:

- a 1-cycle line integral for circulation,
- a 2-surface fluxoid law for the gauge sector,
- and a 3-sphere flux for the radial throughput.

The present stack contains no exact constitutive theorem that turns the 1-cycle winding data into the 3-cycle throughput amplitude.

That is why the mass bottleneck survives unchanged:
\[
\boxed{
F_0^{\rm(eq)}
=
\frac{72\pi^4}{121\Lambda^2}
\frac{\mathcal I_w^2c_s^2\rho}{\Phi^2}.
}
\]

### 5.4 What Section 5 really establishes

Section 5 is easy to summarize and easy to misread.

It does **not** say topology is irrelevant. Topology already fixes the external circulation class, the internal mixed holonomy class, and the same-charge half-flux corridor.

What it does say is narrower and more important:

\[
\boxed{
\text{the remaining missing theorem is not another winding law;}
}
\]
\[
\boxed{
\text{it is an amplitude law for }\Phi.
}
\]

That conclusion was the reason the project turned next to dynamic mouth-output, open-system balance, and eventually the nonlinear turbine tests.

## 6. Dynamic turbine scaling

The next handoff changed the family question in a useful way. Instead of keeping the radial throughput fixed across branches, it treated the support wave as a “turbine” that could dynamically pump more 4D fluid on higher branches. That means the geometry must be allowed to re-equilibrate with the chosen throughput, rather than keeping the aspect ratio \(\Lambda=L/a\) frozen by hand.

This does not rescue the naïve support tower, but it produces a much sharper family theorem.

### 6.1 Dynamic re-equilibration laws

Keep the already-closed reduced mass functional
\[
F_j(a_j)=\frac{A_j}{a_j}+\frac{B_j}{a_j^2}+C_j a_j^3
\]
together with the frozen partition
\[
E_w:E_f:E_{\rm PV}=11:2:5.
\]
For the \(j\)-th D/N support branch define
\[
\nu_j\equiv 2j+1.
\]
Then a convenient reduced parameterization is
\[
A_j=U\,\frac{\nu_j}{\Lambda_j},
\qquad
B_j=V\,\Phi_j^2,
\qquad
C_j=W\,\Lambda_j,
\qquad
\Lambda_j=\frac{L_j}{a_j}.
\]
Here:

- \(A_j/a_j\) is the support-wave ledger,
- \(B_j/a_j^2\) is the compact 4D radial self-flow ledger,
- \(C_ja_j^3\) is the displaced-vacuum / pressure-volume ledger.

Imposing the already-closed ratio
\[
\frac{E_f}{E_w}=\frac2{11}
\]
gives the same stationary formulas as before,
\[
a_j=\frac{11B_j}{2A_j},
\qquad
F_j=\frac{36}{121}\frac{A_j^2}{B_j}.
\]
Substituting the dynamic forms of \(A_j,B_j,C_j\) yields the exact branch laws
\[
\boxed{
R_j\equiv\frac{F_j}{F_0}
=
\nu_j^{1/3}\phi_j^{2/3},
\qquad
\phi_j\equiv\frac{\Phi_j}{\Phi_0},
}
\]
and
\[
\boxed{
\frac{a_j}{a_0}=\nu_j^{-1/6}\phi_j^{2/3},
\qquad
\frac{L_j}{L_0}=\nu_j^{2/3}\phi_j^{-2/3},
\qquad
\frac{\Lambda_j}{\Lambda_0}=\nu_j^{5/6}\phi_j^{-4/3}.
}
\]

So once dynamic throughput is allowed, the family problem is no longer a fixed-geometry support ladder. The geometry and the pumping strength are now locked together.

### 6.2 Why “heavier means deeper” fails inside the reduced closure

These equations immediately imply a no-go.

If a higher branch is genuinely deeper than the electron branch, then
\[
L_j>L_0
\quad\Longrightarrow\quad
\phi_j<\nu_j.
\]
Substituting that into the mass law gives
\[
R_j<\nu_j.
\]
So any deeper branch obeys
\[
\boxed{
R_1<3,
\qquad
R_2<5.
}
\]

That is a very strong statement. Inside the present reduced closure, a branch that is truly **deeper** than the electron can never approach the charged-lepton mass ratios. So the intuitive story “heavier families are deeper throats” is ruled out at this level.

### 6.3 Exact geometry required by any target mass ratio

The dynamic law can be inverted without specifying the turbine model. Solving for \(\phi_j\) gives
\[
\phi_j=\frac{R_j^{3/2}}{\sqrt{\nu_j}}.
\]
Substituting back,
\[
\boxed{
\frac{\Phi_j}{\Phi_0}=\frac{R_j^{3/2}}{\sqrt{\nu_j}},
\qquad
\frac{a_j}{a_0}=\frac{R_j}{\sqrt{\nu_j}},
\qquad
\frac{L_j}{L_0}=\frac{\nu_j}{R_j},
\qquad
\frac{\Lambda_j}{\Lambda_0}=\frac{\nu_j^{3/2}}{R_j^2}.
}
\]

This means that once a target mass ratio is specified, the required dynamic geometry is already fixed by the reduced closure.

For the observed charged-lepton targets,
\[
R_\mu\approx206.77,
\qquad
R_\tau\approx3477.37,
\]
the resulting geometries are not deep. They are extremely **wide and shallow**:
\[
\frac{a_\mu}{a_0}\approx119.38,
\qquad
\frac{L_\mu}{L_0}\approx1.45\times10^{-2},
\qquad
\frac{\Lambda_\mu}{\Lambda_0}\approx1.22\times10^{-4},
\]
and
\[
\frac{a_\tau}{a_0}\approx1555.13,
\qquad
\frac{L_\tau}{L_0}\approx1.44\times10^{-3},
\qquad
\frac{\Lambda_\tau}{\Lambda_0}\approx9.25\times10^{-7}.
\]

So if one insists on reproducing the muon and tau masses inside this dynamic closure, the heavier branches become flatter, not deeper.

### 6.4 Minimal scale-free turbine law

To encode the handoff’s “support wave acts as a turbine” idea without adding a new dimensionful scale, the simplest response law is
\[
\Phi_j
=
\Phi_0\left(\frac{\mathcal P_j}{\mathcal P_0}\right)^s,
\]
where \(s\) is a dimensionless turbine exponent and
\[
\mathcal P_j\equiv E_{w,j}\omega_j
\]
is the natural pumping power.

At equilibrium,
\[
E_{w,j}=\frac{11}{18}F_j,
\]
so
\[
\frac{\mathcal P_j}{\mathcal P_0}
=
\left(\frac{F_j}{F_0}\right)^2
=
R_j^2.
\]
Hence
\[
\phi_j=R_j^{2s}.
\]
Substituting into the dynamic mass law gives the exact self-consistent ladder
\[
\boxed{
R_j=\nu_j^{1/(3-4s)},
\qquad
s<\frac34.
}
\]
This exposes a critical response point
\[
\boxed{
s_c=\frac34.
}
\]

There is a second no-go hidden here. For branches to get deeper with \(j\), one needs
\[
s<\frac12.
\]
But in that regime the mass exponent never exceeds one, so
\[
R_1\le3,
\qquad
R_2\le5.
\]
Therefore any turbine law capable of reaching the real charged-lepton scale must lie in the opposite regime
\[
s>\frac12,
\]
where the higher branches become wider and shallower.

### 6.5 What the turbine fit actually says

Fitting the minimal monomial law to the physical \(\mu/e\) and \(\tau/e\) ratios gives
\[
s_\mu\approx0.69849,
\qquad
s_\tau\approx0.70066,
\]
with a common best fit near
\[
\boxed{
s\approx0.70.
}
\]
So the nonlinear turbine idea can bridge the gap in **order of magnitude**, but only in a regime very close to the critical point \(s=3/4\). It does not naturally produce the real lepton masses with a simple universal power law, and the branches it does produce are geometrically shallow.

### 6.6 What Section 6 really establishes

Section 6 changed the family discussion in an important way.

It showed that the old support-only \(1:9:25\) falsifier was too simple, because once throughput and geometry are allowed to respond dynamically, the family ladder changes.

But it also produced an even stronger negative result:

\[
\boxed{
\text{inside the present reduced closure, heavier cannot mean deeper.}
}
\]

A dynamic turbine law is still logically allowed, but only in a near-critical regime, and in that regime the higher branches are flatter and wider rather than deeper and more cable-like.

So Section 6 does not rescue the naive family picture. It sharpens exactly where it fails and what kind of additional ingredient would be required if one wants to keep a family interpretation alive.

# Lepton / Electron Program — Sections 7–9 Draft

## 7. Exact phase closure versus missing amplitude closure

After the D/N spectral theorem and the family-ladder falsifier, the program split cleanly into two different closure problems.

One of them is a **phase** problem: does the internal support wave really reproduce itself after one full loop, so that the common scalar phase is forced into the trivial class?

The other is an **amplitude** problem: what fixes the actual throughput / mouth-output strength that normalizes the isolated rest energy?

Those two problems are often conflated, but they are not the same. Section 7 records the point where the project stopped treating them as one question.

### 7.1 The exact cycle-map distinction

The support channel can always be written in one-cycle form as
\[
A_{n+1}=\Lambda A_n + D,
\qquad
\Lambda=\rho e^{i\phi_0},
\]
where:

- \(\Lambda\) is the intrinsic round-trip coefficient of the support loop,
- \(\phi_0\) is the common scalar loop phase,
- and \(D\) is an external replenishment / pump term.

This simple equation already shows that three notions of “stability” are different.

A driven or leakage-balanced steady state only requires the amplitude to settle to a fixed point,
\[
A_* = \frac{D}{1-\Lambda}.
\]
That does **not** force
\[
\Lambda=1.
\]
It only means pumping and loss have balanced in whatever open-system way is available.

If one weakens the condition further to mere non-decay of the modulus,
\[
|A_{n+1}|=|A_n|,
\]
one gets only
\[
|\Lambda|=1,
\]
so the phase can still drift.

Only the autonomous self-reproducing eigenmode condition
\[
D=0,
\qquad
A_{n+1}=A_n\neq0
\]
forces
\[
\Lambda=1,
\qquad
\phi_0=2\pi N.
\]
So the exact issue was never “is the amplitude steady?” The exact issue was “is the isolated particle a true autonomous eigenmode?”

### 7.2 What the exact D/N branch contributes

The exact finite-throat D/N branch already supplies a large part of that closure.

Its support spectrum is
\[
\omega_j=\frac{\pi c_s}{L}\left(j+\frac12\right),
\]
and the exact reflection structure is
\[
r_D=-1,
\qquad
r_N=+1.
\]
Therefore the scalar round-trip factor is
\[
R_{\rm rt}=r_Dr_N e^{2ikL}.
\]
On the exact D/N ladder,
\[
2k_jL=(2j+1)\pi,
\]
so
\[
R_{\rm rt}(k_j)=(-1)e^{i(2j+1)\pi}=1.
\]
Equivalently,
\[
\boxed{
\phi_0\equiv0\pmod{2\pi}
}
\]
on the exact trapped D/N eigenmodes.

There is also a useful current check. For
\[
\psi_j(z,t)=A_j\sin(k_j z)e^{-i\omega_j t},
\]
the longitudinal standing-wave current is
\[
J_z\propto \Im\!\left(\psi_j^*\partial_z\psi_j\right)=0.
\]
So the exact D/N branch is a true trapped standing-wave support channel rather than a one-way transport state.

### 7.3 What this does and does not prove

The phase result is strong, but it comes with one important honesty condition.

The exact D/N boundary-value problem proves that **if** the isolated particle lives on that trapped eigenmode branch, then the common scalar phase is killed exactly. But the uploaded project files still do not fully derive, from the complete moving-throat dynamics alone, that an isolated free particle must realize the strong autonomous-eigenmode condition. That is still one of the explicit unresolved sectors of the current stack.

So the carry-forward status is:
\[
\boxed{
\text{phase closure is essentially solved on the exact D/N branch,}
}
\]
but
\[
\boxed{
\text{the model still lacks a first-principles amplitude closure.}
}
\]

### 7.4 Why the mass bottleneck survives

The isolated electron rest-energy theorem had already been reduced to
\[
F_0^{\rm(eq)}
=
\frac{72\pi^4}{121\Lambda^2}
\frac{\mathcal I_w^2 c_s^2\rho}{\Phi^2},
\qquad
\Lambda=\frac{L}{a}.
\]
The exact D/N branch, even with its full phase closure, does not determine the remaining amplitude combination
\[
\frac{\mathcal I_w^2 c_s^2\rho}{\Phi^2}.
\]
So the support-phase theorem does not yet close the electron mass.

This was the point where the program cleanly separated:

- the **phase side**, which is essentially solved on the exact trapped branch,
- from the **amplitude side**, which still needs an open-system mouth-output or throughput law.

### 7.5 What Section 7 really establishes

Section 7 records the first major split in the program.

The exact D/N support problem already does most of the work needed to trivialize the common scalar phase. What remains missing is not “more spectral closure.” It is the actual amplitude law that fixes the size of the mouth-output / throughput channel.

So the correct carry-forward slogan is:
\[
\boxed{
\text{phase closure: nearly in hand; amplitude closure: still open.}
}
\]

## 8. Open-system mouth-output channel

Once the phase problem was separated from the amplitude problem, the next step was to ask what physical channel could actually fix the remaining mass normalization.

The answer was not another internal winding law. The answer had to come from an **open-system mouth-output channel**.

### 8.1 Three flux objects that must not be identified

The corrected ontology makes an immediate distinction between three different flux quantities:

1. the **internal 4D radial throughput** \(\Phi\), which enters the self-flow energy ledger through
   \[
   E_f=\frac{\Phi^2}{8\pi^2\rho a^2},
   \]
2. the **exact projected brane–bulk exchange source** built from the transverse bulk current \(J_w\),
3. and the **brane-side mouth-output flux** measured through a mouth surface and governed only operationally by an effective response law.

Those are not the same object. The old temptation to set
\[
\Phi=S_{\rm leak}
\]
by fiat was therefore conceptually wrong.

### 8.2 Exact D/N ground mode: no trans-brane current, but nonzero mouth gradient

On the exact ground D/N branch,
\[
\chi_0(z)=\sqrt{\frac{2}{L}}\sin\!\left(\frac{\pi z}{2L}\right),
\qquad
\omega_0=\frac{\pi c_s}{2L}.
\]
This mode obeys
\[
\chi_0(0)=0,
\qquad
\chi_0'(0)=\frac{\pi}{\sqrt2\,L^{3/2}},
\qquad
|\chi_0'(0)|^2=\frac{\pi^2}{2L^3}.
\]
So the trapped support wave vanishes at the mouth, but it has a nonzero mouth gradient.

For the exact single trapped mode
\[
\psi_0(z,t)=A\chi_0(z)e^{-i\omega t},
\]
the transverse bulk current is
\[
J_w=
\frac{\hbar}{2im_\psi}
\left(
\psi^*\partial_w\psi-
\psi\,\partial_w\psi^*
\right).
\]
The exact single-mode result is
\[
\boxed{J_w=0.}
\]
Therefore, if one identifies the projected source strictly with the exact current-based leakage channel, then the isolated trapped D/N branch gives
\[
\boxed{S_{\rho,\rm brane}=0.}
\]

This was a major clarification. The exact D/N support wave is not a direct trans-brane current injector. It is a **boundary hammer**.

### 8.3 Why the linear mouth channel is the wrong variable

The exact mouth operator on the D/N branch is
\[
Z_{00}(\omega)
=
-\frac{\omega}{c_s}
\tan\!\left(\frac{\omega L}{c_s}\right).
\]
At the exact ground support frequency
\[
\omega_0=\frac{\pi c_s}{2L},
\]
this operator sits on a pole. So the linear mouth-value channel is not the right constitutive variable to use on the trapped support branch.

That conclusion also follows geometrically: the exact trapped mode has
\[
\psi(0)=0,
\]
so a linear mouth-displacement or mouth-value source term vanishes anyway.

The first nontrivial mouth datum is therefore not the mouth value. It is the **quadratic normal stress** generated by the exact boundary gradient.

### 8.4 Exact quantum-geometric hammer stress

Take the reduced 1D support Hamiltonian density
\[
\mathcal H_{\rm supp}
=
\frac{\rho_\parallel}{2}(\partial_t X)^2
+
\frac{\rho_\parallel c_s^2}{2}(\partial_z X)^2,
\]
with exact standing-wave field
\[
X_0(z,t)=A_0\chi_0(z)\cos(\omega_0 t).
\]
At the mouth \(z=0\), the Dirichlet condition gives
\[
X_0(0,t)=0,
\qquad
\partial_t X_0(0,t)=0,
\]
so the boundary traction comes entirely from the gradient term:
\[
T_{nn}(t)
=
\frac{\rho_\parallel c_s^2}{2}
\big(\partial_zX_0(0,t)\big)^2.
\]
Substituting the exact D/N mode gives
\[
T_{nn}(t)
=
\frac{\rho_\parallel c_s^2A_0^2\pi^2}{4L^3}\cos^2(\omega_0 t)
=
\bar T_{nn}\bigl(1+\cos2\omega_0 t\bigr),
\]
with
\[
\bar T_{nn}
=
\frac{\rho_\parallel c_s^2A_0^2\pi^2}{8L^3}.
\]

Now impose the support-action normalization
\[
\mathcal I_w=\frac{\bar E_w}{\omega_0}=\hbar,
\qquad
\bar E_w=\frac{\rho_\parallel A_0^2\omega_0^2}{2}.
\]
This fixes
\[
A_0^2=\frac{2\hbar}{\rho_\parallel\omega_0}.
\]
Substituting into \(\bar T_{nn}\) cancels the internal support normalization completely and gives the exact hammer theorem
\[
\boxed{
\bar T_{nn}^{(0)}=\frac{\pi\hbar c_s}{2L^2}.
}
\]
So the isolated D/N ground mode, once anchored by \(\mathcal I_w=\hbar\), exerts a parameter-free cycle-averaged mouth hammer stress.

### 8.5 The correct mouth constitutive law

The bridge-response framework only requires a generalized port law of the form
\[
j_i(\omega)=\sum_j Z^{\rm eff}_{ij}(\omega)u_j(\omega),
\]
with the effort variable declared by protocol.

For the present jackhammer problem, the natural effort is the **normal mouth stress**, not the vanishing mouth value. So the correct constitutive variable is
\[
u_\sigma(t)=T_{nn}(t).
\]
The corresponding stress-port law is
\[
j_{\rm mouth}(\omega)=Z_\sigma^{\rm eff}(\omega)u_\sigma(\omega).
\]
Since the hammer stress contains both a DC and a \(2\omega_0\) part,
\[
u_\sigma(t)=\bar T_{nn}^{(0)}+\bar T_{nn}^{(0)}\cos(2\omega_0 t),
\]
the first coarse-grained mouth injector is
\[
\boxed{
j_{\rm mouth}^{\rm(dc)}
=
Z_\sigma^{\rm eff}(0)
\frac{\pi\hbar c_s}{2L^2}.
}
\]
This was the first explicit analytical mouth-output law available within the present stack.

### 8.6 What still remains unresolved

The hammer stress is exact, but the constitutive coefficient is not. The project files still do not derive the actual low-frequency mouth response coefficient
\[
Z_\sigma^{\rm eff}(0)
\]
from first principles.

So Section 8 closes the **load**, but not the **response**. It narrows the amplitude bottleneck without removing it.

### 8.7 What Section 8 really establishes

Section 8 records the transition from a vague “maybe there is a leakage term” picture to a much sharper statement:

- the exact trapped D/N mode is not a trans-brane current source,
- it is a resonant mouth hammer,
- the correct mouth datum is the quadratic normal stress,
- and that stress has the exact cycle-averaged amplitude
  \[
  \bar T_{nn}^{(0)}=\frac{\pi\hbar c_s}{2L^2}.
  \]

So the remaining amplitude problem is no longer “find some source somehow.” It is:
\[
\boxed{
\text{derive the mouth response to the exact hammer stress.}
}
\]

## 9. Static brane response: compliance, not DC flow

Section 8 suggested a natural next step: perhaps the exact DC hammer stress could be converted into a nonzero static mouth source by the static GNLS compressibility of the brane.

Section 9 showed that this is **not** what the conservative theory does.

The conservative static theory produces a scalar compliance. It does **not** produce a DC flux mobility.

### 9.1 Static GNLS compliance law

Work around the uniform GNLS brane background with the already-frozen stiff-polytropic EOS
\[
P(\rho)=K\rho^5,
\qquad
U(\rho)=\frac{K}{4}\rho^5.
\]
Then
\[
U''(\rho)=5K\rho^3,
\qquad
c_s^2(\rho)=\frac{1}{m}\frac{dP}{d\rho}=\frac{5K}{m}\rho^4,
\]
so at the background density \(\rho_0\),
\[
\boxed{U''(\rho_0)=\frac{mc_s^2}{\rho_0}.}
\]

Write
\[
\rho=\rho_0(1+\eta),
\qquad
|\eta|\ll1,
\]
and keep only the static scalar density channel. The quadratic conservative energy is then
\[
\delta^2E
=
\int d^3x\left[
\frac{\rho_0mc_s^2}{2}\eta^2
+
\frac{\rho_0\hbar^2}{8m}|\nabla\eta|^2
-
 s(\mathbf x)\eta
\right],
\]
where \(s(\mathbf x)\) is the effective mouth load.

Varying gives the static Euler–Lagrange equation
\[
\rho_0mc_s^2\eta
-
\frac{\rho_0\hbar^2}{4m}\nabla^2\eta
=
 s(\mathbf x),
\]
or
\[
\boxed{
(1-\ell^2\nabla^2)\eta
=
\frac{s(\mathbf x)}{\rho_0mc_s^2},
\qquad
\ell^2=\frac{\hbar^2}{4m^2c_s^2}.
}
\]
With the standard acoustic healing-length convention
\[
\xi_h=\frac{\hbar}{\sqrt2mc_s},
\]
one has
\[
\ell=\frac{\xi_h}{\sqrt2}.
\]
So the exact conservative static response is a Yukawa / Helmholtz **compliance** problem.

### 9.2 Why the static flux mobility vanishes

The GNLS current is
\[
\mathbf j=\rho\mathbf v=\frac{\hbar}{m}\rho\nabla\theta.
\]
In a static conservative minimization, the phase appears only through the positive term
\[
\frac{\hbar^2\rho}{2m}|\nabla\theta|^2.
\]
For a bounded localized static response with no imposed phase winding in the brane-response channel, the minimum is therefore
\[
\nabla\theta=0,
\qquad
\mathbf j=0.
\]
Equivalently, the steady continuity equation
\[
\nabla\cdot(\rho\nabla\theta)=0
\]
with finite energy and no topological constraint admits only the constant-phase branch.

So the quantity closed by the conservative EOS is **not** a flux admittance. It is a compliance. In particular,
\[
\boxed{Z_{\sigma,\rm flux}^{\rm eff}(0)=0}
\]
in the conservative static sector.

This is one of the clearest theorems worth carrying forward. It means the DC mouth source cannot come from static EOS mobility alone.

### 9.3 What the \(P_{22}\) branch can still do

The pre-existing \(P_{22}\) mouth branch does matter, but only in a limited way at zero frequency.

Model the scalar / quadrupole response near equilibrium by a reduced Hessian
\[
H=
\begin{pmatrix}
A & \varepsilon B\\
\varepsilon B & C
\end{pmatrix},
\qquad
f=
\begin{pmatrix}
f_0\\0
\end{pmatrix},
\]
where \(f_0\) is the scalar load and \(\varepsilon\sim\sigma_\infty\) measures the intrinsic quadrupolar pre-stress.

Then the scalar compliance is
\[
\chi_{00}^{\rm eff}
\equiv
\frac{q_0}{f_0}
=
\frac{C}{AC-\varepsilon^2B^2}.
\]
Expanding for small \(\varepsilon\),
\[
\boxed{
\chi_{00}^{\rm eff}
=
\frac1A+
\frac{\varepsilon^2B^2}{A^2C}+O(\varepsilon^4).
}
\]
So the \(P_{22}\) branch can renormalize the scalar compliance, but only at **even order** in the pre-stress. There is no linear scalar-to-quadrupole generation of DC flow, and still no static mouth source follows.

### 9.4 Consequence for the mass theorem

At the end of Section 8, the mass theorem looked temptingly close to closure because the hammer load was exact:
\[
\bar T_{nn}^{(0)}=\frac{\pi\hbar c_s}{2L^2}.
\]
Section 9 shows that if one tries to interpret the brane response coefficient as a static **flux** mobility, then the conservative theory gives
\[
Z_{\sigma,\rm flux}^{\rm eff}(0)=0.
\]
So the static conservative route cannot close the electron mass.

The correct conservative zero-frequency law is instead a compliance law of the form
\[
\delta V = \chi_{00}^{\rm eff}(0)\,\bar T_{nn}^{(0)},
\]
which can renormalize geometry or scalar stiffness but does not itself fix the absolute throat scale.

### 9.5 What Section 9 really establishes

Section 9 is the theorem that protects the GNLS vacuum from “melting” under static hammering.

The conservative brane behaves like an elastic mattress, not like a fluid with a DC bleed channel. The exact static object is a real scalar compliance, and the corresponding static flux mobility vanishes:
\[
\boxed{Z_{\sigma,\rm flux}^{\rm eff}(0)=0.}
\]

That result narrows the remaining bottleneck dramatically. The missing constitutive law is no longer “some static mouth mobility from EOS.” It has to be a genuinely **dynamic** rectification law — radiation, leakage, streaming, memory, or another open-system AC\(\to\)DC conversion mechanism.

So the correct carry-forward statement is:
\[
\boxed{
\text{static mouth mobility is the wrong target; dynamic rectification is the real one.}
}
\]

# Lepton / Electron Program — Sections 10–12 Draft

## 10. Dynamic AC→DC rectification law

Section 9 closed the conservative static response problem and showed that the zero-frequency mouth response is a **compliance**, not a DC source law:
\[
Z^{\rm eff}_{\sigma,\mathrm{flux}}(0)=0.
\]
That immediately changes the logic of the open-system sector. If a nonzero coarse-grained mouth injector exists,
\[
j_{\rm mouth}^{\rm(dc)}\neq0,
\]
it cannot come from linear static GNLS elasticity. It must come from a genuinely **dynamic** rectification mechanism.

The exact D/N ground branch already supplies the necessary periodic hammer. From the earlier mouth-stress theorem,
\[
T_{nn}(t)=\bar T\bigl(1+\cos 2\omega_0 t\bigr),
\qquad
\bar T=\frac{\pi\hbar c_s}{2L^2},
\qquad
\omega_0=\frac{\pi c_s}{2L}.
\]
So the oscillatory drive seen by the brane is
\[
T_{nn}^{\rm(ac)}(t)=T_a\cos(\Omega t),
\qquad
T_a=\frac{\pi\hbar c_s}{2L^2},
\qquad
\Omega=2\omega_0=\frac{\pi c_s}{L}.
\]
The question of Section 10 is therefore precise: can second-order GNLS acoustics convert this AC hammer into a nonzero cycle-averaged DC mass current?

### 10.1 Second-order GNLS acoustics

Take a one-dimensional normal coordinate \(x\) away from the mouth and write the standard compressible GNLS hydrodynamics in conservative form,
\[
\partial_t\rho+\partial_x(\rho u)=0.
\]
Expand in a small acoustic amplitude parameter \(\varepsilon\):
\[
\rho=\rho_0+\varepsilon\rho_1+\varepsilon^2\rho_2+\cdots,
\qquad
u=\varepsilon u_1+\varepsilon^2u_2+\cdots.
\]
At first order,
\[
\partial_t\rho_1+\rho_0\partial_xu_1=0.
\]
At second order,
\[
\partial_t\rho_2+\rho_0\partial_xu_2+\partial_x(\rho_1u_1)=0.
\]
Cycle-averaging over one period kills the explicit time derivative, so the exact DC conservation law becomes
\[
\partial_xJ_{\rm dc}=0,
\]
with
\[
\boxed{
J_{\rm dc}=\rho_0\overline{u_2}+\langle\rho_1u_1\rangle.
}
\]
This is the fundamental rectification object. It says that second-order acoustics can carry a genuine DC mass current even when the linear acoustics have zero time average.

That is the first theorem of the section:
\[
\boxed{
\text{the first place a DC mouth current can appear is second order in the acoustic response.}
}
\]

### 10.2 Dispersive linear acoustics and the normal-stress impedance

The linearized GNLS momentum / Bernoulli equation gives the standard Bogoliubov dispersion law,
\[
\Omega^2=c_s^2k^2+\frac{\hbar^2k^4}{4m^2}.
\]
Using the healing-length convention
\[
\xi_h=\frac{\hbar}{\sqrt2\,m c_s},
\]
this becomes
\[
\boxed{
\Omega^2=c_s^2k^2\left(1+\frac12\xi_h^2k^2\right).
}
\]
The corresponding phase speed is
\[
c_{\rm ph}=\frac{\Omega}{k}.
\]
The linear normal stress carried by the acoustic field is then
\[
\Pi_1=\rho_0\frac{\Omega}{k}u_1=\rho_0c_{\rm ph}u_1.
\]
At the exact hammer frequency
\[
\Omega=\frac{\pi c_s}{L},
\]
solving the dispersion relation for \(c_{\rm ph}\) gives
\[
\boxed{
c_{\rm ph}^2
=
\frac{c_s^2}{2}
\left(
1+\sqrt{1+\frac{2\pi^2\xi_h^2}{L^2}}
\right).
}
\]
So the exact hammer is not driving a perfectly nondispersive acoustic channel; the healing length explicitly modifies the effective stress-to-velocity conversion.

### 10.3 Closed standing branch versus outgoing radiative branch

The crucial split comes here.

For a **closed conservative standing wave**, take
\[
u_1^{\rm st}=U\sin(kx)\sin(\Omega t),
\qquad
\rho_1^{\rm st}=\frac{\rho_0U}{c_{\rm ph}}\cos(kx)\cos(\Omega t).
\]
Then
\[
\langle \rho_1^{\rm st}u_1^{\rm st}\rangle=0.
\]
So the exact trapped conservative branch carries no rectified DC mouth flux.

By contrast, for a minimal **outgoing / radiative branch**, take
\[
u_1^{\rm out}=U\cos(kx-\Omega t),
\qquad
\rho_1^{\rm out}=\frac{\rho_0U}{c_{\rm ph}}\cos(kx-\Omega t).
\]
Then the quadratic transport term survives:
\[
\langle \rho_1^{\rm out}u_1^{\rm out}\rangle
=
\frac{\rho_0U^2}{2c_{\rm ph}}.
\]
Matching the normal stress to the hammer amplitude,
\[
T_a=\rho_0c_{\rm ph}U,
\qquad
U=\frac{T_a}{\rho_0c_{\rm ph}},
\]
gives the local cycle-averaged mouth-flux density
\[
\boxed{
j_{\rm mouth}^{\rm(dc,rad)}
=
\frac{T_a^2}{2\rho_0c_{\rm ph}^3}.
}
\]
Substituting the exact hammer amplitude,
\[
T_a=\frac{\pi\hbar c_s}{2L^2},
\]
one gets
\[
\boxed{
j_{\rm mouth}^{\rm(dc,rad)}
=
\frac{\pi^2\hbar^2c_s^2}{8\rho_0c_{\rm ph}^3L^4}.
}
\]
Or, eliminating \(c_{\rm ph}\) explicitly,
\[
\boxed{
j_{\rm mouth}^{\rm(dc,rad)}(L)
=
\frac{\pi^2\hbar^2}{8\rho_0c_sL^4}
\left[
\frac{2}{1+\sqrt{1+2\pi^2\xi_h^2/L^2}}
\right]^{3/2}.
}
\]

That is the second theorem of the section:
\[
\boxed{
\text{a nonzero DC mouth source exists only on an explicitly open / radiative branch.}
}
\]

### 10.4 Integrated mouth source versus internal sink

The conservative throat geometry still gives the compact internal sink law
\[
\Phi^2=
\frac{8\pi^3\rho_0\hbar c_s}{11\Lambda^2}
L,
\qquad
\Phi\propto L^{1/2},
\]
with
\[
\Lambda=\frac{L}{a}.
\]
For the mouth output, the relevant area is the scalar mouth area
\[
A_{\rm mouth}=\pi a^2=\pi\frac{L^2}{\Lambda^2}.
\]
So the integrated radiative mouth source is
\[
\dot M_{\rm mouth}^{\rm(dc,rad)}
=
A_{\rm mouth}\,j_{\rm mouth}^{\rm(dc,rad)}
=
\frac{\pi^3\hbar^2c_s^2}{8\Lambda^2\rho_0c_{\rm ph}^3L^2}.
\]
Using the explicit phase speed again,
\[
\boxed{
\dot M_{\rm mouth}^{\rm(dc,rad)}(L)
=
\frac{\pi^3\hbar^2}{8\Lambda^2\rho_0c_sL^2}
\left[
\frac{2}{1+\sqrt{1+2\pi^2\xi_h^2/L^2}}
\right]^{3/2}.
}
\]
The scaling comparison is immediate:
\[
\dot M_{\rm mouth}^{\rm(dc,rad)}\propto L^{-2},
\qquad
\Phi\propto L^{1/2}.
\]
So the two channels necessarily cross at some positive depth. That observation is the bridge to Section 11.

### 10.5 What Section 10 really establishes

Section 10 did not finish the electron mass theorem by itself, but it closed a crucial missing component.

It proved that:

1. static GNLS elasticity cannot generate a DC mouth current,
2. second-order acoustics can,
3. the closed standing-wave branch still carries zero DC mouth current,
4. the first nonzero mouth injector is the **open/radiative** branch,
5. and its exact scaling law is known.

The carry-forward formula worth remembering is
\[
\boxed{
j_{\rm mouth}^{\rm(dc,rad)}(L)
=
\frac{\pi^2\hbar^2}{8\rho_0c_sL^4}
\left[
\frac{2}{1+\sqrt{1+2\pi^2\xi_h^2/L^2}}
\right]^{3/2}.
}
\]
This formula is what finally made it meaningful to search for a unique selected electron scale.

## 11. Unique equilibrium crossover and selected electron scale

Once the dynamic mouth-source law is known, the electron scale-selection problem becomes sharp. The internal sink increases with depth, while the radiative mouth source decreases with depth. A stable isolated particle should therefore sit at the unique crossover point where the two exactly balance.

### 11.1 The equilibrium crossover condition

Take the two exact stage-10 channels:
\[
\dot M_{\rm mouth}^{\rm(dc,rad)}(L)
=
\frac{\pi^3\hbar^2}{8\Lambda^2\rho_0c_sL^2}
\left[
\frac{2}{1+\sqrt{1+2\pi^2\xi_h^2/L^2}}
\right]^{3/2},
\]
\[
\Phi(L)=
\sqrt{\frac{8\pi^3\rho_0\hbar c_s}{11\Lambda^2}}\,L^{1/2}.
\]
The equilibrium selection principle is simply
\[
\boxed{
\dot M_{\rm mouth}^{\rm(dc,rad)}(L_*)=\Phi(L_*).
}
\]
This is the first real scale-selection principle in the chain that does not insert \(L\) by hand.

### 11.2 Uniqueness of the selected scale

Define
\[
s=\sqrt{1+\frac{2\pi^2\xi_h^2}{L^2}}>1.
\]
Then the mouth-source logarithmic derivative is
\[
\frac{d\ln \dot M_{\rm mouth}}{dL}
=
-\frac{s+3}{2Ls}<0,
\]
while the sink derivative is
\[
\frac{d\Phi}{dL}
=
\frac12\sqrt{\frac{8\pi^3\rho_0\hbar c_s}{11\Lambda^2}}\,L^{-1/2}>0.
\]
So:

- \(\dot M_{\rm mouth}\) is strictly decreasing,
- \(\Phi\) is strictly increasing.

At small depth,
\[
L\to0^+:
\qquad
\dot M_{\rm mouth}\to\infty,
\qquad
\Phi\to0.
\]
At large depth,
\[
L\to\infty:
\qquad
\dot M_{\rm mouth}\to0,
\qquad
\Phi\to\infty.
\]
Therefore there is exactly one positive crossover scale:
\[
\boxed{
\text{the equilibrium-selection equation has a unique positive solution }L_*.
}
\]

### 11.3 Exact crossover equation and dimensionless reduction

After algebraic simplification, the exact crossover equation is
\[
\boxed{
64\Lambda^2\rho_0^3c_s^3L_*^2
\left(L_*+\sqrt{L_*^2+2\pi^2\xi_h^2}\right)^3
=
11\pi^3\hbar^3.
}
\]
This can be put into a cleaner dimensionless form by defining
\[
z=L+\sqrt{L^2+2\pi^2\xi_h^2},
\qquad
z=\sqrt2\,\pi\xi_h y,
\qquad y>1.
\]
Then the exact selection law becomes
\[
\boxed{
y(y^2-1)^2=\mathcal C,
\qquad
\mathcal C=
\frac{11\hbar^3}{64\sqrt2\,\pi^2\Lambda^2\rho_0^3c_s^3\xi_h^5}.
}
\]
Since
\[
\frac{d}{dy}\Big[y(y^2-1)^2\Big]=(y^2-1)(5y^2-1)>0
\qquad (y>1),
\]
this dimensionless equation also has a unique root \(y_*>1\).

The selected depth is therefore
\[
\boxed{
L_*=
\frac{\pi\xi_h}{\sqrt2}
\frac{y_*^2-1}{y_*}.
}
\]

### 11.4 Selected electron mass law

The stage-8 isolated rest-energy theorem had already reduced the electron mass on any chosen branch to
\[
F_0(L)=\frac{9\pi\hbar c_s}{11L}.
\]
Substituting the selected depth gives
\[
\boxed{
F_0(L_*)
=
\frac{9\sqrt2\,\hbar c_s}{11\xi_h}
\frac{y_*}{y_*^2-1},
}
\]
where \(y_*\) is fixed implicitly by the crossover law above.

There is also a useful long-wave approximation. If
\[
L_*\gg\xi_h,
\]
then the exact equation reduces to
\[
L_*
\approx
\left(
\frac{11\pi^3\hbar^3}{512\Lambda^2\rho_0^3c_s^3}
\right)^{1/5}.
\]
Substituting that back gives
\[
\boxed{
F_0(L_*)
\approx
\frac{18\,2^{4/5}}{11^{6/5}}
\pi^{2/5}\Lambda^{2/5}\hbar^{2/5}\rho_0^{3/5}c_s^{8/5}.
}
\]

### 11.5 What Section 11 really establishes

Section 11 is easy to state carefully.

It **does** prove that the dynamic mouth-source law plus the conservative sink law select a unique electron depth. So the scale \(L\) is no longer arbitrary.

But it does **not** produce a parameter-free numerical electron mass in the strongest possible sense, because the selected scale still inherits the GNLS background data
\[
\rho_0,
\qquad c_s,
\qquad \xi_h.
\]
That is not a flaw. A dimensionful rest energy must inherit some dimensionful vacuum rulers. What Section 11 actually closes is the following theorem:
\[
\boxed{
\text{the electron branch has a unique dynamically selected scale }L_*,
\text{ but that scale is still vacuum-background dependent.}
}
\]
That conclusion set up the next pivot: if absolute mass depends on the background, perhaps the **family ratios** are the cleaner place to search for pure geometric structure.

## 12. Family-ratio theorems and no-gos

Section 12 generalized the crossover principle from the ground branch to the whole D/N tower and asked whether the heavier charged leptons could be obtained by applying the same equilibrium-selection rule to higher support harmonics.

The answer is one of the strongest negative results in the whole program.

### 12.1 The generalized crossover family equation

Keep the same support-action anchor on each isolated branch,
\[
\mathcal I_w=\hbar,
\]
and define
\[
\nu_j\equiv 2j+1,
\qquad
R_j\equiv\frac{F_j}{F_0}.
\]
The exact D/N support frequency is
\[
\omega_j=\frac{\pi c_s\nu_j}{2L_j},
\]
so the isolated rest energy on branch \(j\) is
\[
F_j=\frac{9\pi\hbar c_s\nu_j}{11L_j}.
\]
Using the already-derived re-equilibration identities,
\[
\frac{L_j}{L_0}=\frac{\nu_j}{R_j},
\qquad
\frac{\Lambda_j}{\Lambda_0}=\frac{\nu_j^{3/2}}{R_j^2},
\]
and reimposing the same dynamic crossover selection on every branch yields the exact family equation
\[
\boxed{
R_j^3
=
\nu_j^{5/3}
\frac{1+\sqrt{1+\alpha R_j^2}}
     {1+\sqrt{1+\alpha}},
\qquad
\alpha=2\pi^2\left(\frac{\xi_h}{L_{*,0}}\right)^2.
}
\]
All dimensionful background scales cancel. The entire family problem collapses to one universal dimensionless vacuum parameter \(\alpha\).

### 12.2 Exact universal bounds

The ratio factor obeys the sharp inequality
\[
1
\le
\frac{1+\sqrt{1+\alpha R_j^2}}
     {1+\sqrt{1+\alpha}}
\le R_j
\qquad
(\alpha\ge0,\ R_j\ge1).
\]
Substituting this into the exact family equation yields the universal bounds
\[
\boxed{
\nu_j^{5/9}\le R_j\le \nu_j^{5/6}.
}
\]
For the first two excited D/N branches,
\[
\nu_1=3,
\qquad
\nu_2=5,
\]
one gets
\[
\boxed{
1.841057547\le R_1\le2.498049533,
}
\]
\[
\boxed{
2.445212848\le R_2\le3.823622457.
}
\]
These are exact consequences of the current crossover-family mechanism.

### 12.3 Comparison with the physical charged-lepton ratios

The charged-lepton targets are
\[
R_\mu=\frac{m_\mu}{m_e}\approx206.768,
\qquad
R_\tau=\frac{m_\tau}{m_e}\approx3477.365.
\]
They are nowhere near the allowed D/N crossover windows
\[
R_1=O(1),
\qquad
R_2=O(1).
\]
So the central theorem of the section is
\[
\boxed{
\text{the simple }e,\mu,\tau\leftrightarrow j=0,1,2\text{ D/N support-tower hypothesis is decisively falsified.}
}
\]
This is not a soft “not yet derived” statement. It is a real no-go theorem for the simplest family mechanism.

### 12.4 The parameter-free endpoint ladders

The exact family equation has two useful pure-geometric endpoint limits.

If
\[
\alpha\to0,
\]
then
\[
\boxed{
R_j=\nu_j^{5/9}.
}
\]
If
\[
\alpha\to\infty,
\]
then
\[
\boxed{
R_j=\nu_j^{5/6}.
}
\]
So the parameter-free endpoint values are:
\[
R_1^{\rm IR}=3^{5/9}\approx1.8411,
\qquad
R_2^{\rm IR}=5^{5/9}\approx2.4452,
\]
\[
R_1^{\rm UV}=3^{5/6}\approx2.4980,
\qquad
R_2^{\rm UV}=5^{5/6}\approx3.8236.
\]
Any exact solution of the present D/N crossover-family equation must lie between these two ladders.

There is also a robustness check: if one freezes \(\Lambda\) instead of letting the geometry re-equilibrate, the family equation changes, but the failure becomes even stronger. The corresponding bound tightens to
\[
\nu_j^{2/5}\le R_j\le\nu_j^{2/3},
\]
which still leaves the physical muon and tau completely out of range.

### 12.5 What Section 12 really establishes

Section 12 did two things at once.

First, it made the family problem more elegant by showing that the crossover ratios depend only on a universal dimensionless vacuum parameter \(\alpha\), not on all the background units separately.

Second, it proved that the most obvious family interpretation already fails:
\[
\boxed{
\text{charged-lepton families are not the plain D/N harmonic tower plus the same crossover selector.}
}
\]
So Section 12 is one of the key places where the project gains information by falsification. It does not leave the family problem shapeless. It tells us exactly what sort of family mechanism cannot be right, which is why the later stages moved toward reverse-engineering, integer-action sectors, and multi-constraint lattice intersections instead.

# Lepton / Electron Program — Sections 13–15 Draft

## 13. Reverse-engineered lepton families under the universal 11:2:5 ledger

Sections 10–12 showed that the isolated-electron branch can be selected by a dynamic crossover condition, but they also showed that the plain D/N family ladder fails badly if one tries to identify
\[
e,\mu,\tau \leftrightarrow j=0,1,2
\]
directly. The next move was therefore to invert the problem phenomenologically while **keeping the stationary reduced ledger fixed**.

The carry-forward ledger is
\[
E_w:E_f:E_{\rm PV}=11:2:5,
\qquad
F=E_w+E_f+E_{\rm PV}=\frac{18}{11}E_w.
\]
This is the crucial lock. If a heavier branch still lives on the same reduced stationary ontology, then the observed mass ratio must equal the total support-wave energy ratio:
\[
R\equiv \frac{F}{F_0},
\qquad
Q\equiv \frac{E_w}{E_{w,0}}.
\]
Since
\[
F=\frac{18}{11}E_w,
\qquad
F_0=\frac{18}{11}E_{w,0},
\]
one gets the exact identification
\[
\boxed{R=Q.}
\]
This means the empirical lepton mass ratio is not just another observable. Inside the frozen reduced hierarchy, it **is** the required total support-energy multiplier.

### 13.1 Action, frequency, and the combined support multiplier

To split that required support-energy multiplier into microscopic pieces, introduce two factors.

First, let
\[
W\equiv \frac{\mathcal I_w}{\mathcal I_{w,0}}
\]
be the support-action or topological multiplier.

Second, let
\[
\nu
\]
be the D/N spatial-frequency multiplier. For the original D/N tower,
\[
\nu=2j+1.
\]

The useful combined quantity is then
\[
\boxed{\Sigma\equiv W\nu.}
\]

This is the quantity that actually enters the stationary throat geometry. The support coefficient scales as
\[
A\propto \frac{\mathcal I_w\,\nu}{\Lambda},
\qquad
\Lambda\equiv \frac{L}{a},
\]
so the reduced support energy scales as
\[
Q=\frac{E_w}{E_{w,0}}=\frac{\Sigma}{\lambda},
\qquad
\lambda\equiv \frac{L}{L_0}.
\]
Because the ledger already fixes
\[
R=Q,
\]
the stationary reduced geometry becomes
\[
\boxed{
\frac{L}{L_0}=\frac{\Sigma}{R},
\qquad
\frac{a}{a_0}=\frac{R}{\sqrt{\Sigma}},
\qquad
\frac{\Lambda}{\Lambda_0}=\frac{\Sigma^{3/2}}{R^2}.
}
\]

This is one of the most important family formulas in the whole program. It says that the stationary geometry is controlled primarily by the **combined support multiplier** \(\Sigma=W\nu\), not by topology and frequency separately.

### 13.2 Exact action–frequency split

The D/N support law still gives
\[
\omega=\frac{\pi c_s\nu}{2L}.
\]
Therefore
\[
\frac{\omega}{\omega_0}
=
\frac{\nu}{L/L_0}
=
\frac{\nu}{\Sigma/R}
=
\frac{R}{W}.
\]
So the support-energy multiplier factors exactly as
\[
R=Q
=\left(\frac{\mathcal I_w}{\mathcal I_{w,0}}\right)
 \left(\frac{\omega}{\omega_0}\right)
=W\,\frac{R}{W}.
\]
This means:

- the action/topological multiplier is \(W\),
- the frequency multiplier is \(\omega/\omega_0=R/W\),
- the stationary geometry only sees \(\Sigma=W\nu\).

So once the mass ratio is fixed by observation, increasing \(W\) lowers the required support frequency, and decreasing \(W\) forces the support frequency upward.

### 13.3 Exact crossover inversion

The dynamic crossover law from Sections 10–11 can be rewritten in terms of the support multipliers. Let
\[
\alpha\equiv 2\pi^2\left(\frac{\xi_h}{L_0}\right)^2.
\]
Then the generalized branch equation is
\[
\boxed{
R^3
=
(W\nu)^{5/3}
\frac{1+\sqrt{1+\alpha (R/W)^2}}
     {1+\sqrt{1+\alpha}}.
}
\]
Using \(\Sigma=W\nu\), the same equation can also be written as
\[
\boxed{
R^3
=
\Sigma^{5/3}
\frac{1+\sqrt{1+\alpha (R\nu/\Sigma)^2}}
     {1+\sqrt{1+\alpha}}.
}
\]
This is the master inversion formula for the heavy-lepton branches inside the frozen reduced hierarchy.

### 13.4 Frequency-dominant versus action-dominant branches

There are two qualitatively different regimes.

#### Frequency-dominant branch: \(W<R\)

If \(W<R\), then
\[
\frac{R}{W}>1.
\]
In that regime the crossover equation implies the lower bound
\[
W\ge \frac{R^3}{\nu^{5/2}}.
\]
Consistency with \(W<R\) therefore requires
\[
\boxed{\nu>R^{4/5}.}
\]
This is a very strong restriction. It says that a frequency-dominant interpretation is possible only if the D/N support harmonic is already extremely high.

#### Action-dominant / low-harmonic branch: \(W\ge R\)

If \(W\ge R\), then \(R/W\le1\), and the crossover factor is bounded above by one. Therefore
\[
R^3\le (W\nu)^{5/3},
\]
which gives the exact lower bound
\[
\boxed{W\nu\ge R^{9/5}.}
\]
This is the regime that matters for low odd harmonics such as \(\nu=1,3,5\). In that regime, the action/topology multiplier must do almost all of the heavy lifting.

### 13.5 Universal low-harmonic benchmark

On the minimal low-harmonic branch, the lower bound saturates:
\[
W\nu=R^{9/5}.
\]
Substituting into the exact geometry formulas gives a universal benchmark:
\[
\boxed{
\frac{L}{L_0}=R^{4/5},
\qquad
\frac{a}{a_0}=R^{1/10},
\qquad
\frac{\Lambda}{\Lambda_0}=R^{7/10}.
}
\]
The support-frequency ratio on this benchmark branch is
\[
\frac{\omega}{\omega_0}=\nu R^{-4/5}.
\]
So on a low-harmonic heavy branch the frequency is actually **lower** than on the electron branch, not higher.

### 13.6 Muon and tau benchmarks

Using the charged-lepton mass anchors
\[
R_\mu\approx 206.7682829877,
\qquad
R_\tau\approx 3477.3652666018,
\]
the frequency-dominant consistency thresholds are
\[
\nu_\mu>R_\mu^{4/5}\approx 71.1848,
\qquad
\nu_\tau>R_\tau^{4/5}\approx 680.7697.
\]
So the first few D/N harmonics cannot possibly carry the muon or tau by frequency uplift alone.

On the minimal low-harmonic branch, the resulting benchmark geometries are:

For the muon:
\[
\frac{L_\mu}{L_0}\approx 71.1848,
\qquad
\frac{a_\mu}{a_0}\approx 1.7043,
\qquad
\frac{\Lambda_\mu}{\Lambda_0}\approx 41.7675.
\]

For the tau:
\[
\frac{L_\tau}{L_0}\approx 680.7697,
\qquad
\frac{a_\tau}{a_0}\approx 2.2601,
\qquad
\frac{\Lambda_\tau}{\Lambda_0}\approx 301.2140.
\]

The corresponding minimal topological-action multipliers are enormous:
\[
W_{\mu,\min}=\frac{R_\mu^{9/5}}{\nu},
\qquad
W_{\tau,\min}=\frac{R_\tau^{9/5}}{\nu}.
\]
For low harmonics \(\nu=1,3,5\), this puts the muon in the \(W\sim10^3\)–\(10^4\) range and the tau in the \(W\sim10^5\)–\(10^6\) range.

### 13.7 What Section 13 really establishes

Section 13 is the point where the family story changes qualitatively.

The exact results worth carrying forward are:
\[
\boxed{R=Q,}
\qquad
\boxed{\Sigma=W\nu,}
\qquad
\boxed{
\frac{L}{L_0}=\frac{\Sigma}{R},
\;
\frac{a}{a_0}=\frac{R}{\sqrt{\Sigma}},
\;
\frac{\omega}{\omega_0}=\frac{R}{W}.
}
\]

The physical meaning is equally important:

1. the universal stationary ledger locks mass ratio and support-energy ratio together exactly;
2. low-order D/N harmonics cannot explain the muon or tau by frequency uplift alone;
3. any low-harmonic heavy-lepton branch must be **action/topology dominated**;
4. and the resulting heavy branches are predicted to be **very deep**, **only mildly wider**, and of **much larger aspect ratio** than the electron branch.

This is the first place where the heavy leptons start to look like what later handoffs called “deep needles” or “wrapped cables,” rather than just higher support tones.

## 14. Multi-constraint intersection picture

Section 13 produced an exact inversion, but by itself it still allowed a continuum of candidate branches. The next conceptual step was to recognize that a true particle family cannot be generated by scaling one knob continuously. A stable branch must sit at a **simultaneous intersection** of several independent constraints.

Inside the present reduced hierarchy, the relevant constraints are:

1. the universal stationary ledger \(E_w:E_f:E_{\rm PV}=11:2:5\),
2. the dynamic crossover condition from Sections 10–11,
3. an integer action/topology multiplier \(W\in\mathbb Z_{>0}\),
4. an odd D/N support multiplier \(\nu\in2\mathbb N+1\),
5. and ultimately a universal background parameter \(\alpha\) shared by all branches.

This turns the family problem into a **lattice problem** rather than a continuum fit.

### 14.1 Exact branch equation in lattice form

With
\[
\Sigma=W\nu,
\qquad
\alpha=2\pi^2\left(\frac{\xi_h}{L_0}\right)^2,
\]
the exact multi-constraint branch law is
\[
R^3
=
\Sigma^{5/3}
\frac{1+\sqrt{1+\alpha (R\nu/\Sigma)^2}}
     {1+\sqrt{1+\alpha}}.
\]
The geometry remains
\[
\frac{L}{L_0}=\frac{\Sigma}{R},
\qquad
\frac{a}{a_0}=\frac{R}{\sqrt{\Sigma}},
\qquad
\frac{\Lambda}{\Lambda_0}=\frac{\Sigma^{3/2}}{R^2},
\qquad
\frac{\omega}{\omega_0}=\frac{R}{W}.
\]
So \(\Sigma\) still controls geometry, while \(W\) separately controls the action–frequency split.

### 14.2 Exact admissible strip on the low-harmonic branch

The physically interesting regime is the action-dominant / low-harmonic strip,
\[
\frac{R\nu}{\Sigma}\le1.
\]
In that regime the crossover factor is bounded by
\[
\frac{R\nu}{\Sigma}
\le
\frac{1+\sqrt{1+\alpha (R\nu/\Sigma)^2}}
     {1+\sqrt{1+\alpha}}
\le 1.
\]
This yields the exact strip
\[
\boxed{R^{9/5}\le \Sigma \le \frac{R^3}{\nu^{3/2}}.}
\]
Since \(\Sigma=W\nu\), integer topology turns this into the explicit lattice condition
\[
\boxed{W\in\mathbb Z_{>0},\quad \nu\in2\mathbb N+1,\quad R^{9/5}\le W\nu\le \frac{R^3}{\nu^{3/2}}.}
\]
This is the cleanest mathematical expression of the “simultaneous intersection” ontology.

### 14.3 First admissible low-harmonic integer candidates

Using the observed charged-lepton mass ratios,
\[
R_\mu\approx 206.7682829877,
\qquad
R_\tau\approx 3477.3652666018,
\]
one can ask for the first admissible integer values of \(W\) at low odd harmonics \(\nu=1,3,5\).

For the muon, the first admissible low-harmonic integer candidates are:

| \(\nu\) | first integer \(W\) | \(\Sigma=W\nu\) | \(L/L_0\) | \(a/a_0\) | \(\Lambda/\Lambda_0\) | \(\omega/\omega_0\) | required \(\alpha\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 14719 | 14719 | 71.18597 | 1.70430 | 41.76857 | 0.0140477 | \(1.1127\times10^{-4}\) |
| 3 | 4907 | 14721 | 71.19564 | 1.70418 | 41.77708 | 0.0421374 | \(1.0192\times10^{-3}\) |
| 5 | 2944 | 14720 | 71.19080 | 1.70424 | 41.77283 | 0.0702338 | \(5.6707\times10^{-4}\) |

For the tau, the first admissible low-harmonic integer candidates are:

| \(\nu\) | first integer \(W\) | \(\Sigma=W\nu\) | \(L/L_0\) | \(a/a_0\) | \(\Lambda/\Lambda_0\) | \(\omega/\omega_0\) | required \(\alpha\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2367286 | 2367286 | 680.77001 | 2.26009 | 301.21422 | 0.00146893 | \(2.5837\times10^{-6}\) |
| 3 | 789096 | 2367288 | 680.77059 | 2.26008 | 301.21460 | 0.00440677 | \(8.2162\times10^{-6}\) |
| 5 | 473458 | 2367290 | 680.77116 | 2.26008 | 301.21498 | 0.00734461 | \(1.3849\times10^{-5}\) |

### 14.4 What these first integer ceilings mean

The most striking feature of these tables is that, for each species separately, the first admissible low-harmonic integers have almost the **same** \(\Sigma\), and therefore almost the same geometry, across \(\nu=1,3,5\).

So once the observed mass ratio is fixed, the geometry is nearly fixed as well. Changing \(\nu\) mostly just trades support budget between action and frequency:

- larger \(\nu\) lowers the required integer \(W\),
- but the combined geometric multiplier \(\Sigma=W\nu\) barely moves.

This is why the family problem starts to look like a sparse lattice problem. The observed mass ratio selects a narrow geometric site, and the integer tumblers \((W,\nu)\) must land on that site closely enough to satisfy the dynamic crossover equation.

### 14.5 Why there is no clean small-winding interpretation

Section 13 already showed that a frequency-dominant branch requires
\[
\nu>R^{4/5}.
\]
For the observed charged-lepton ratios this means
\[
\nu_\mu>71.1848,
\qquad
\nu_\tau>680.7697.
\]
So there is no clean explanation in terms of small \(W\) and a modest harmonic uplift.

Likewise, on the low-harmonic branch the first admissible integer ceilings already require
\[
W_\mu\sim 10^3\text{--}10^4,
\qquad
W_\tau\sim 10^5\text{--}10^6.
\]
So no simple law such as \(W=2\), \(W=3\), or \(W=j+1\) emerges. The heavy leptons are not small-winding variants of the electron branch.

### 14.6 First common-background test

A true family theorem would require the same background invariant \(\alpha\) for all branches. The first low-harmonic integer ceilings fail that test badly. For example:

- the muon \(\nu=1\) ceiling needs \(\alpha\approx1.1127\times10^{-4}\),
- the tau \(\nu=1\) ceiling needs \(\alpha\approx2.5837\times10^{-6}\).

The same mismatch persists for \(\nu=3\) and \(\nu=5\).

So the first admissible integer ceilings are **not** yet a family theorem. They are only the first indication that the family problem is sparse and highly constrained.

### 14.7 What Section 14 really establishes

Section 14 sharpens the ontology of family replication.

The exact carry-forward statement is:
\[
\boxed{
\text{a stable branch is not selected by one scale factor, but by a simultaneous intersection of}
}
\]
\[
\boxed{
\text{stationary ledger} + \text{dynamic crossover} + \text{integer topology} + \text{common background.}
}
\]

The mathematical form of that statement is the admissible strip
\[
\boxed{R^{9/5}\le W\nu\le \frac{R^3}{\nu^{3/2}}.}
\]

This is the point where the family story stops being a one-parameter scaling exercise and becomes a genuine sparse-lattice problem.

## 15. Universal-\(\alpha\) lattice search

Section 14 made the next task inevitable. If the electron, muon, and tau live in the same GNLS background, then the dimensionless vacuum parameter
\[
\alpha=2\pi^2\left(\frac{\xi_h}{L_0}\right)^2
\]
must be the same for all three. The problem is therefore to search the integer lattice \((W,\nu)\) for simultaneous charged-lepton branches with one common \(\alpha\).

### 15.1 Exact inversion for \(\alpha\)

The branch equation can be solved exactly for \(\alpha\). Define
\[
k\equiv \frac{R^3}{\Sigma^{5/3}},
\qquad
\beta\equiv \frac{R\nu}{\Sigma}.
\]
Then the exact inversion is
\[
\boxed{
\alpha(R,\Sigma,\nu)
=
\frac{4k\,(k-\beta^2)(1-k)}{(k^2-\beta^2)^2}.
}
\]
This is important because it turns the family problem into a direct arithmetic comparison problem. One no longer has to solve a transcendental equation numerically at every lattice point.

### 15.2 Universal-vacuum variable and near-integer structure

A second useful reformulation uses
\[
x\equiv \sqrt{1+\alpha},
\qquad
t\equiv \left(\frac{\Sigma}{R^{9/5}}\right)^{1/3},
\qquad
\gamma\equiv \frac{\nu}{R^{4/5}}.
\]
Then the branch equation reduces to
\[
\boxed{2t^5+\gamma^2(x-1)t^4-(x+1)=0.}
\]
So a universal vacuum means:

1. choose one common \(x=\sqrt{1+\alpha}\),
2. solve the species-specific quintic for \(t\),
3. and require that
   \[
   W=\frac{R^{9/5}t^3}{\nu}
   \in \mathbb Z_{>0}
   \]
   for both species.

This shows that the family problem is really a **simultaneous near-integer problem**.

### 15.3 Search domain and admissible counts

The explicit search was carried out over the action-dominant strip with:

- muon odd harmonics
  \[
  1\le \nu_\mu\le 71,
  \]
- tau odd harmonics
  \[
  1\le \nu_\tau\le 679,
  \]
- and all integer \(W\) such that the implied \(\alpha\le1\).

Inside that domain, the admissible counts are:

- **4404** muon tuples,
- **1,028,044** tau tuples.

Every admissible muon tuple was then compared against the nearest tau tuples in \(\alpha\).

### 15.4 No exact equality found in the searched box

Within this full search box, no exact equality
\[
\alpha_\mu=\alpha_\tau
\]
was found to numerical precision.

That is an important negative result. It means the present reduced equations do **not** yet pick out a unique exact family theorem on their own inside the searched domain.

But the search also found very tight near-intersections.

### 15.5 Best precision near-family

The tightest pair in the search box is
\[
(\nu_\mu,W_\mu)=(3,5025),
\qquad
(\nu_\tau,W_\tau)=(1,2424686).
\]
These give
\[
\alpha_\mu=0.16959239403167644\ldots,
\qquad
\alpha_\tau=0.16959239404618672\ldots,
\]
so
\[
\boxed{|\Delta\alpha|\approx1.4510\times10^{-11}.}
\]
This is the best precision family found inside the search domain.

Its common-vacuum estimate is
\[
\alpha_*\approx0.16959239403893158,
\qquad
\frac{L_0}{\xi_h}\approx10.7885.
\]
The corresponding branch geometries are:

Muon:
\[
\Sigma_\mu=15075,
\quad
\frac{L_\mu}{L_0}\approx72.9077,
\quad
\frac{a_\mu}{a_0}\approx1.68405,
\quad
\frac{\Lambda_\mu}{\Lambda_0}\approx43.2930,
\quad
\frac{\omega_\mu}{\omega_0}\approx0.04115.
\]

Tau:
\[
\Sigma_\tau=2424686,
\quad
\frac{L_\tau}{L_0}\approx697.277,
\quad
\frac{a_\tau}{a_0}\approx2.23317,
\quad
\frac{\Lambda_\tau}{\Lambda_0}\approx312.236,
\quad
\frac{\omega_\tau}{\omega_0}\approx0.001434.
\]

This is the clearest realization of the “Leviathan” picture: very large action, very low support frequency, and very deep bulk profiles.

### 15.6 Best balanced near-family

The search also found a much more balanced near-family:
\[
(\nu_\mu,W_\mu)=(7,2126),
\qquad
(\nu_\tau,W_\tau)=(115,20809).
\]
These give
\[
\alpha_\mu=0.07635277420936142\ldots,
\qquad
\alpha_\tau=0.07635277408891019\ldots,
\]
so
\[
\boxed{|\Delta\alpha|\approx1.2045\times10^{-10}.}
\]
This is looser by about one order of magnitude, but the required integers are dramatically smaller.

Its common-vacuum estimate is
\[
\alpha_*\approx0.07635277414913581,
\qquad
\frac{L_0}{\xi_h}\approx16.0788.
\]
The corresponding geometries are:

Muon:
\[
\frac{L_\mu}{L_0}\approx71.9743,
\qquad
\frac{a_\mu}{a_0}\approx1.69494,
\qquad
\frac{\Lambda_\mu}{\Lambda_0}\approx42.4643,
\qquad
\frac{\omega_\mu}{\omega_0}\approx0.09726.
\]

Tau:
\[
\frac{L_\tau}{L_0}\approx688.175,
\qquad
\frac{a_\tau}{a_0}\approx2.24789,
\qquad
\frac{\Lambda_\tau}{\Lambda_0}\approx306.142,
\qquad
\frac{\omega_\tau}{\omega_0}\approx0.16711.
\]

This shows that the universal-vacuum search does **not** uniquely force the huge-\(W\), ultra-slow-frequency corner. The family lattice is richer than that.

### 15.7 Least-complexity near-family

A third representative point is the least-complexity candidate found with a relatively loose but still small mismatch:
\[
(\nu_\mu,W_\mu)=(23,664),
\qquad
(\nu_\tau,W_\tau)=(653,3637),
\]
with
\[
|\Delta\alpha|\approx3.6002\times10^{-6}.
\]
This point is less precise, but it shows that once larger odd harmonics are allowed, the integer action need not be astronomically large.

### 15.8 What Section 15 really establishes

Section 15 changes the family picture in two important ways.

First, it shows that the simultaneous-intersection ontology is not empty. The integer lattice contains extremely tight common-background near-families, including one with
\[
|\Delta\alpha|\sim10^{-11}.
\]
So the family problem is genuinely sparse and structured rather than hopelessly diffuse.

Second, it shows that the universal-vacuum search does **not** yet select a unique family theorem. The current reduced equations allow several quite different near-family corners:

- a “Leviathan” corner with huge \(W\) and tiny \(\omega/\omega_0\),
- a more balanced corner with moderate \(W\) and higher odd harmonics,
- and still other less precise corners at larger \(\nu\).

What survives across all of them is the broad geometric trend:

- heavy leptons are **much deeper** than the electron branch,
- they are only **mildly wider**,
- and they have a **much larger aspect ratio**.

So the honest carry-forward conclusion is
\[
\boxed{
\text{the present theory has reached a sparse simultaneous near-integer family lattice,}
}
\]
not yet
\[
\boxed{
\text{a unique exact Diophantine family theorem.}
}
\]

That means one more selector is still missing. The obvious candidates are a secondary torsion/parity rule, a decay-width or stability selector, or some additional topological congruence that collapses the near-family lattice to one exact family.

# Lepton / Electron Program — Sections 16–18 Draft

## 16. Resonance and phase-lock no-gos

Sections 13–15 produced a sparse family lattice of deep heavy-branch candidates, but they did not yet explain **why** only rare branches would be stable. The next natural idea was a coupled-oscillator selector: perhaps a charged lepton survives only when the 4D longitudinal support wave phase-locks to a natural 3D mouth rebound mode.

That program produced two leading-order no-gos worth carrying forward.

The first is a symmetry veto: the exact scalar mouth hammer does **not** linearly drive the area-preserving \(P_{22}\) ellipse.

The second is a direct resonance veto: the representative deep-needle family candidates do **not** directly phase-lock even to the scalar \(P_0\) breathing branch.

So the known mouth oscillators do not presently supply the missing family selector.

### 16.1 The exact mouth hammer oscillates at \(2\omega_j\)

For the exact D/N support branch,
\[
X_j(z,t)\propto \sin(k_j z)\cos(\omega_j t),
\qquad
\omega_j=\frac{\pi c_s\nu_j}{2L_j},
\qquad
\nu_j=2j+1.
\]
The mouth normal traction is quadratic in the support gradient:
\[
T_{nn}(t)\propto (\partial_z X_j)^2.
\]
Therefore the oscillatory part of the hammer is
\[
\boxed{
T_{nn}^{\rm(ac)}(t)=T_a\cos(2\omega_j t).
}
\]
This matters because any direct forced resonance must be tested against the drive frequency \(2\omega_j\), not \(\omega_j\).

### 16.2 The scalar \(P_0\) hammer does zero linear work on the \(P_{22}\) mouth mode

Parameterize the headless quadrupole mouth by
\[
a_+=ae^{\sigma},
\qquad
 a_-=ae^{-\sigma}.
\]
Its area is exactly
\[
A_{\rm mouth}=\pi a_+a_- = \pi a^2,
\]
so it is independent of the ellipticity \(\sigma\).

A spatially isotropic pressure load contributes
\[
U_P(t)=-P(t)A_{\rm mouth}=-P(t)\pi a^2,
\]
and therefore
\[
\boxed{
\frac{\partial U_P}{\partial \sigma}=0.
}
\]
So the exact scalar mouth hammer cannot linearly excite the area-preserving quadrupole:
\[
\boxed{P_0\text{ hammer} \not\to P_{22}\text{ oscillator at linear order.}}
\]
This is one of the cleanest structural results in the entire program. It means the simplest 3D/4D impedance-matching picture is not just numerically unsupported; it is symmetry-forbidden at leading order.

### 16.3 Natural \(P_{22}\) rebound frequency from dynamic GNLS compliance

Although the direct forcing vanishes, the natural quadrupole rebound frequency can still be computed from the dynamic extension of the GNLS compliance law.

The scalar linearized GNLS branch obeys the Bogoliubov dispersion
\[
\Omega^2=c_s^2k^2\left(1+\frac12\xi_h^2k^2\right).
\]
For the headless quadrupole boundary deformation
\[
\delta r(\theta,t)=a\,\delta\sigma(t)\cos2(\theta-\phi),
\]
the azimuthal wavenumber is
\[
k_2=\frac{2}{a}.
\]
Hence the natural mouth rebound frequency is
\[
\boxed{
\Omega_{22}(a)=\frac{2c_s}{a}\sqrt{1+2\frac{\xi_h^2}{a^2}}.
}
\]
The intrinsic quadrupolar bracer from the isolated-mouth potential,
\[
U_{\rm iso}(\sigma)=\frac12k_{22}\sigma^2+\frac14u_{22}\sigma^4-h_\alpha\sinh2\sigma,
\]
shifts the equilibrium but does not change the oscillation frequency at first order. So the direct lock criterion is
\[
2\omega_j=\Omega_{22,j},
\]
or equivalently
\[
\boxed{
\frac{2\omega_j}{\Omega_{22,j}}
=
\frac{\pi\nu_j}{2\Lambda_j\sqrt{1+2(\xi_h/a_j)^2}},
\qquad
\Lambda_j\equiv\frac{L_j}{a_j}.
}
\]
Since direct lock requires this ratio to be near one, deep high-aspect-ratio states are already in trouble.

### 16.4 Representative deep-needle families fail direct \(P_{22}\) lock badly

The representative near-family candidates from the universal-\(\alpha\) search all sit far below direct quadrupole lock.

For the best-precision near-family,
\[
\frac{2\omega_\mu}{\Omega_{22,\mu}}\approx 0.0582,
\qquad
\frac{2\omega_\tau}{\Omega_{22,\tau}}\approx 0.00270.
\]
For the best balanced near-family,
\[
\frac{2\omega_\mu}{\Omega_{22,\mu}}\approx 0.139,
\qquad
\frac{2\omega_\tau}{\Omega_{22,\tau}}\approx 0.318.
\]
So none of the Stage-15 deep needles naturally align with direct \(P_{22}\) resonance.

This can also be expressed as a threshold. Direct lock needs roughly
\[
\Lambda_j^{\rm lock}\lesssim 1.57\,\nu_j,
\]
while the deep-needle candidates have
\[
\Lambda_j\sim 10^2\text{ to }10^3.
\]
They are simply too deep.

### 16.5 The symmetry-correct retest: the scalar \(P_0\) breathing mode

Because the exact D/N hammer is scalar, the correct direct-lock target is the scalar breathing branch rather than the quadrupole.

For an axisymmetric scalar mouth mode on a circular mouth of radius \(a\), the regular linear GNLS mode is
\[
\eta(r)=J_0(kr),
\]
with the same Bogoliubov dispersion
\[
\Omega^2=c_s^2k^2\left(1+\frac12\xi_h^2k^2\right).
\]
Two natural rim conditions were tested:

- free-edge / Neumann breathing, \(J_1(ka)=0\), so \(k=j_{1,1}/a\);
- Dirichlet robustness check, \(J_0(ka)=0\), so \(k=j_{0,1}/a\).

Writing \(\chi\in\{j_{1,1},j_{0,1}\}\), the scalar breathing frequency is
\[
\boxed{
\Omega_{00}^{(\chi)}(a)
=
\frac{\chi c_s}{a}
\sqrt{1+\frac12\chi^2\frac{\xi_h^2}{a^2}}.
}
\]
The direct scalar lock condition is then
\[
2\omega_j = k\,\Omega_{00,j},
\qquad k\in\mathbb Z_{>0},
\]
or
\[
\boxed{
\frac{2\omega_j}{\Omega_{00,j}^{(\chi)}}
=
\frac{\pi\nu_j}
{\chi\Lambda_j\sqrt{1+\frac12\chi^2(\xi_h/a_j)^2}}.
}
\]
For any direct harmonic lock, this ratio must be at least one.

### 16.6 The deep-needle families also fail direct \(P_0\leftrightarrow P_0\) lock

The same representative candidates remain far below direct scalar lock.

For the best-precision near-family:
\[
\frac{2\omega_\mu}{\Omega_{00,\mu}^{(D)}}\approx 0.0482,
\qquad
\frac{2\omega_\mu}{\Omega_{00,\mu}^{(N)}}\approx 0.0296,
\]
\[
\frac{2\omega_\tau}{\Omega_{00,\tau}^{(D)}}\approx 0.00224,
\qquad
\frac{2\omega_\tau}{\Omega_{00,\tau}^{(N)}}\approx 0.00139.
\]
For the balanced near-family:
\[
\frac{2\omega_\mu}{\Omega_{00,\mu}^{(D)}}\approx 0.1156,
\qquad
\frac{2\omega_\mu}{\Omega_{00,\mu}^{(N)}}\approx 0.0719,
\]
\[
\frac{2\omega_\tau}{\Omega_{00,\tau}^{(D)}}\approx 0.2643,
\qquad
\frac{2\omega_\tau}{\Omega_{00,\tau}^{(N)}}\approx 0.1649.
\]
Even the strongest case only reaches about \(0.264\), still far below the \(k=1\) fundamental lock.

A generic threshold estimate on the minimal action-dominant branch gives
\[
\nu_j \gtrsim \frac{\chi\Lambda_0}{\pi}R_j^{7/10}.
\]
For the observed charged-lepton ratios this means approximately
\[
\nu_\mu \gtrsim 59.15\text{ (Dirichlet)},\quad 94.24\text{ (Neumann)},
\]
\[
\nu_\tau \gtrsim 426.56\text{ (Dirichlet)},\quad 679.66\text{ (Neumann)}.
\]
So the low-harmonic deep-needle branches found earlier do not meet the scalar breathing lock threshold either.

### 16.7 The tearing story also fails in both channels

The tearing logic turned out to go in the opposite direction from the original handoff.

For the quadrupole channel, the exact scalar hammer does not linearly drive \(\sigma\) at all, so there is no current-order resonant growth of the ellipticity. Even if one inserts a hypothetical anisotropic transfer coefficient \(\varepsilon_{02}\), the driven response behaves like
\[
|\delta\sigma|\sim
\frac{\varepsilon_{02}T_a}
{I_{22}\sqrt{(\Omega_{22}^2-4\omega_j^2)^2+\cdots}},
\qquad
T_a\propto L^{-2},
\]
so off resonance
\[
|\delta\sigma|\propto L^{-2}.
\]
Thus deeper branches are **less** likely to tear through the \(P_{22}\) channel.

For the scalar breathing mode, the hammer does couple linearly, but away from resonance the fractional strain scales as
\[
|b_j|\sim \frac{T_a}{\rho_0c_s^2}
\frac{1}{|1-(2\omega_j/\Omega_{00,j})^2|}.
\]
Since the deep-needle candidates all satisfy \(2\omega_j\ll \Omega_{00,j}\), the enhancement factor stays near one and the strain scales essentially as
\[
\boxed{|b_j|/|b_0|\sim (L_0/L_j)^2.}
\]
So the deeper states are again **safer**, not more fragile.

Therefore neither direct mouth-oscillator channel currently produces a natural upper-depth tearing cutoff that would stop the family ladder after three generations.

### 16.8 What Section 16 really establishes

Section 16 leaves three sharp conclusions.

First,
\[
\boxed{\text{the exact scalar D/N hammer does not directly }P_0\to P_2\text{ lock at leading order.}}
\]

Second,
\[
\boxed{\text{the known deep-needle family candidates do not directly }P_0\leftrightarrow P_0\text{ lock either.}}
\]

Third,
\[
\boxed{\text{simple mouth-mode resonance is not the missing family selector.}}
\]

So if a resonance selector exists, it must be more subtle than a direct harmonic match. The remaining possibilities are things like anisotropic transfer, parametric / Floquet conversion, mixed two-step coupling, or an additional topological selector layered on top of the sparse heavy-lepton lattice.

## 17. Highest-confidence conclusions to carry forward

By the end of the program, some results had become much more reliable than others. This section separates the carry-forward theorems from the weaker speculative branches.

### 17.1 Exact or near-exact reduced theorems that survived every pivot

The following results are the strongest parts of the chain and should be treated as the main carry-forward backbone.

#### (a) The reduced isolated-defect mass ledger is rigid

Within the frozen reduced hierarchy, the stationary isolated defect obeys
\[
F(a,\rho)=\frac{A(\rho)}{a}+\frac{B(\rho)}{a^2}+C(\rho)a^3,
\]
with virial relation
\[
E_w+2E_f=3E_{\rm PV}.
\]
The fixed closure then forces
\[
\boxed{E_w:E_f:E_{\rm PV}=11:2:5,}
\]
and therefore
\[
\boxed{F=\frac{18}{11}E_w.}
\]
This is the single most important thermodynamic theorem in the entire project.

#### (b) The reduced coefficients have clean microscopic meanings

The support-wave coefficient and the 4D self-flow coefficient are
\[
\boxed{A=\mathcal I_w\chi_w c_s,}
\qquad
\boxed{B=\frac{\Phi^2}{8\pi^2\rho}.}
\]
So the reduced mass functional is no longer just algebraic bookkeeping; it has a direct support-plus-self-flow interpretation.

#### (c) The exact D/N support branch is real, not heuristic

The finite-throat D/N spectrum is
\[
\boxed{k_j=\frac{\pi}{L}\left(j+\frac12\right),
\qquad
\omega_j=\frac{\pi c_s}{L}\left(j+\frac12\right).}
\]
And its round-trip coefficient is exactly
\[
\boxed{R_{\rm rt}=1,\qquad \phi_0\equiv0\pmod{2\pi}.}
\]
So the scalar support phase is essentially closed.

#### (d) The trapped D/N mode is a mouth hammer, not a trans-brane current

For the exact ground branch, the linear transverse current vanishes,
\[
J_w=0,
\]
but the mouth gradient is nonzero, producing a cycle-averaged hammer stress
\[
\boxed{\bar T_{nn}^{(0)}=\frac{\pi\hbar c_s}{2L^2}.}
\]
That is the cleanest closed mouth theorem we obtained.

#### (e) Static GNLS brane response is compliance, not DC flow

The conservative zero-frequency brane problem is a Yukawa/Helmholtz compliance problem, not a transport law. Therefore
\[
\boxed{Z_{\sigma,\mathrm{flux}}^{\mathrm{eff}}(0)=0.}
\]
So any nonzero coarse-grained mouth injector must be dynamic and second order or higher.

#### (f) The first nonzero mouth injector is the open/radiative second-order branch

The explicit AC\(\to\)DC mouth law is
\[
\boxed{
j_{\rm mouth}^{\rm(dc,rad)}(L)
=
\frac{\pi^2\hbar^2}{8\rho_0c_sL^4}
\left[
\frac{2}{1+\sqrt{1+2\pi^2\xi_h^2/L^2}}
\right]^{3/2}.
}
\]
That is the first genuine dynamical source law in the whole mouth-output chain.

#### (g) The electron crossover scale is unique once the radiative branch is admitted

Balancing the radiative mouth output against the inward sink selects a unique crossover depth \(L_*\). The selected electron branch is therefore not arbitrary at the reduced level.

### 17.2 Strong no-gos that are now probably permanent within the present hierarchy

Several family-generation ideas were not just left unfinished. They were actively falsified inside the current reduced stack.

#### (a) Plain support-harmonic family replication fails

If one identifies
\[
e,\mu,\tau \leftrightarrow j=0,1,2
\]
on the exact D/N support tower, the self-consistent isolated-particle ladder is
\[
\frac{F_j}{F_0}=(2j+1)^2.
\]
So the plain harmonic-family prediction is
\[
1:9:25:\dots
\]
which is nowhere near the real charged-lepton ratios.

#### (b) Topology fixes circulation class, not throughput amplitude

The exact fluxoid law quantizes the vortical / holonomy sector, but the stationary radial throughput \(\Phi\) remains an integration constant of continuity, not a topologically fixed quantity.

So the missing mass theorem is an amplitude law, not another winding law.

#### (c) Simple dynamic-turbine scaling does not rescue the family ladder

Allowing \(\Phi_j\) to scale dynamically changes the geometry, but within the reduced ladder it also yields the sharp result that a genuinely deeper branch must satisfy
\[
R_j<\nu_j.
\]
So “heavier means deeper” fails in the simplest turbine interpretation.

#### (d) Direct phase lock to the known mouth oscillators fails

The exact scalar hammer does not linearly drive the area-preserving quadrupole, and the deep-needle candidates also fail direct scalar breathing lock.

So the presently known mouth resonances do not select the lepton family lattice.

### 17.3 The best current heavy-lepton picture

If one insists on reverse-engineering the observed charged-lepton ratios while preserving the universal stationary ledger, the most robust conclusion is that heavy branches are controlled by the combined support multiplier
\[
\boxed{\Sigma=W\nu,}
\]
with
\[
\boxed{R=Q=\frac{E_w}{E_{w,0}}.}
\]
Inside that inversion, low-harmonic heavy branches are action/topology dominated rather than frequency dominated. Their characteristic geometry is:

- much deeper than the electron branch,
- only mildly wider,
- and much larger in aspect ratio.

The universal-\(\alpha\) search then turns the family problem into a sparse near-integer lattice problem. That is genuine progress, but it is not yet a unique exact family theorem.

### 17.4 What still remains open

The clean open problems are now narrow.

1. **Amplitude closure for the support/throat feed.**
   The exact support phase is basically closed, but the amplitude law that fixes \(\Phi\) from first principles is still missing.

2. **Open-system radiation / transport microphysics.**
   The radiative mouth branch exists, but the actual outgoing fraction and its microscopic selection law remain open.

3. **A second selector for the family lattice.**
   The heavy-lepton lattice still needs one more exact condition — topological, parity-like, mixed-sector, or decay-based — to collapse the sparse near-family lattice to a unique theorem.

4. **The full moving-throat PDE and strong-field completion.**
   Many background and transport objects are still symbolic because the deeper PDE sector remains unsolved.

### 17.5 What Section 17 really establishes

Section 17 is the carry-forward filter.

The exact results worth remembering are not “everything we tried.” They are:
\[
\boxed{E_w:E_f:E_{\rm PV}=11:2:5,}
\qquad
\boxed{F=\frac{18}{11}E_w,}
\qquad
\boxed{R_{\rm rt}=1,}
\qquad
\boxed{\bar T_{nn}^{(0)}=\frac{\pi\hbar c_s}{2L^2},}
\qquad
\boxed{Z_{\sigma,\mathrm{flux}}^{\mathrm{eff}}(0)=0.}
\]
And the strongest falsifiers worth remembering are:
\[
\boxed{1:9:25\text{ harmonic families fail,}}
\qquad
\boxed{\text{throughput is not fixed by topology alone,}}
\qquad
\boxed{\text{direct mouth-mode phase lock fails for the deep needles.}}
\]

## 18. Best chapter split for the next session

The material is now large enough that the next session should not treat it as one continuous note. The cleanest carry-forward structure is a three-chapter split plus a short appendix.

### 18.1 Chapter I — The isolated-electron theorem chain

This chapter should contain the tightest, most theorem-like part of the project: the reduced isolated-particle mass program.

Recommended contents:

1. frozen reduced hierarchy and source anchors,
2. reduced mass functional and virial relation,
3. exact \(11{:}2{:}5\) partition,
4. microscopic identification of \(A\) and \(B\),
5. exact D/N support spectrum and phase closure,
6. mouth-hammer theorem,
7. static-compliance theorem,
8. second-order radiative mouth law,
9. unique crossover selection of the electron branch.

This chapter is the closest thing the current stack has to a finished theorem chain.

### 18.2 Chapter II — Family no-gos and protected exclusions

This chapter should collect all the mechanisms that were genuinely tested and found not to work inside the present hierarchy.

Recommended contents:

1. failure of the plain support-only \((2j+1)^2\) tower,
2. failure of charge/circulation arguments to fix \(\Phi\),
3. failure of simple dynamic-turbine/deeper-branch scaling,
4. symmetry veto on direct \(P_0\to P_2\) lock,
5. failure of direct \(P_0\leftrightarrow P_0\) breathing lock,
6. failure of simple tearing-cutoff explanations.

This chapter matters because it preserves the project’s falsification-first discipline. It makes clear what future sessions do **not** need to re-test unless the underlying ontology changes.

### 18.3 Chapter III — Heavy-lepton reverse engineering and sparse lattice structure

This chapter should contain the most speculative but still useful carry-forward part: the reverse-engineered heavy-lepton geometry program.

Recommended contents:

1. exact lock \(R=Q\) under the universal stationary ledger,
2. combined multiplier \(\Sigma=W\nu\),
3. action-dominant versus frequency-dominant branches,
4. deep-needle benchmark geometry,
5. multi-constraint intersection picture,
6. universal-\(\alpha\) lattice inversion,
7. best near-family candidates and their geometry,
8. why a second selector is still needed.

This chapter is the right place to continue if the next session wants to keep searching for a family theorem.

### 18.4 Recommended writing order for the next session

The most efficient order is not the chronological one.

The best sequence is:

1. **Write Chapter I first.** It contains the strongest closed results and gives the cleanest standalone theorem chain.
2. **Write Chapter II second.** It records the strongest no-gos and prevents the project from re-running already-falsified routes.
3. **Write Chapter III last.** It is the exploratory branch and should be framed explicitly as a reverse-engineering / lattice program rather than a finished derivation.

### 18.5 Best appendix split

A short appendix should carry the narrower technical branches that are important but not central to the main theorem chain.

Suggested appendix items:

- exact D/N mouth operator and pole structure,
- the circulation / fluxoid audit,
- the detailed \(\alpha\)-inversion formulas,
- the near-family candidate tables,
- the phase-lock ratio tables,
- and the unresolved-symbol ledger \((\rho_0,c_s,\xi_h,\Phi, Z^{\rm eff}, \ldots)\).

### 18.6 What Section 18 really establishes

The chapter split is not just organizational convenience. It reflects the actual confidence structure of the work.

- **Chapter I** is the reduced electron theorem chain.
- **Chapter II** is the protected no-go ledger.
- **Chapter III** is the heavy-lepton search frontier.

That is probably the cleanest way to carry this program into another session without losing which results are hard-won theorems, which are falsifiers, and which are still exploratory.
