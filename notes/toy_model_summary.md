1pn_orbital_dynamics.tex

## 1) What this paper is trying to accomplish

The paper asks a very specific question:

**Can a (mostly) scalar, superfluid/hydrodynamic toy model reproduce the *GR 1PN perihelion precession* while keeping the *0PN (Newtonian) gravity exactly instantaneous* in the near zone?**

The architectural stance is:

* **0PN / Newtonian gravity** is enforced by an *instantaneous Poisson sector* (solves the aberration problem by construction).
* **Finite propagation speed** enters via an additional “lag/retarded” scalar sector, but the paper argues it *does not* generate the observed 1PN precession for a central field.
* Therefore, to match GR’s 1PN perihelion advance, the missing piece is put into an **effective inertia / kinetic prefactor** of the orbiting defect.

---

## 2) Core objects and notation (things you’ll keep using)

### Matter / sources

Mass density is represented as a sum over moving “defects” (smoothed point masses):
[
\rho(\mathbf{x},t)=\sum_{i=1}^N m_i,W!\left(\mathbf{x}-\mathbf{x}_i(t)\right).
]

### Potentials: split into instantaneous + lag

Total scalar potential is split as
[
\Phi(\mathbf{x},t)=\Phi_{\mathrm{P}}(\mathbf{x},t)+\Phi_{\mathrm{L}}(\mathbf{x},t),
]
where (\Phi_{\mathrm{P}}) is the instantaneous Poisson piece, and (\Phi_{\mathrm{L}}) carries the retarded / radiative corrections.

### Key constants

* (G): effective gravitational constant in the toy model.
* (\cS): propagation speed of the scalar sector (interpretable as a “sound speed” in the medium). When matching GR, the identification used is (\cS=c).
* (\mu \equiv GM): the usual Kepler parameter for a central mass.

---

## 3) Governing field equations (the backbone)

### Instantaneous “0PN” sector

[
\nabla^2 \Phi_{\mathrm{P}}(\mathbf{x},t)=4\pi G,\rho(\mathbf{x},t).
]

### Lag / “retarded” sector (driven by time-variation of (\Phi_{\mathrm{P}}))

[
\frac{\partial^2 \Phi_{\mathrm{L}}}{\partial t^2}
=\cS^2\nabla^2\Phi_{\mathrm{L}}-\frac{\partial^2\Phi_{\mathrm{P}}}{\partial t^2}.
]

A numerically convenient equivalent is to evolve the **total** potential by the sourced wave equation:
[
\frac{\partial^2 \Phi}{\partial t^2}
=\cS^2\left[\nabla^2\Phi-4\pi G,\rho\right],
]
and then compute (\Phi_{\mathrm{P}}) by FFT Poisson solve each step, defining (\Phi_{\mathrm{L}}=\Phi-\Phi_{\mathrm{P}}) in post-processing.

### Force law on defects

[
\mathbf{a}(\mathbf{x},t)=-\nabla\Phi(\mathbf{x},t)
= -\nabla\Phi_{\mathrm{P}}-\nabla\Phi_{\mathrm{L}}.
]

---

## 4) The Static Limit Theorem (why “instantaneous gravity” is consistent here)

A central claim is:

> For time-independent sources, the lag field decays away and the total potential reduces to the instantaneous Poisson solution.

Concretely, in the static limit:
[
\nabla^2\Phi(\mathbf{x})=4\pi G,\rho(\mathbf{x}),
\qquad
\Phi(\mathbf{x})=\Phi_{\mathrm{P}}(\mathbf{x}),\qquad
\Phi_{\mathrm{L}}(\mathbf{x})\to 0.
]

This is the paper’s built-in solution to the classic “aberration” worry: **the near-zone, leading-order gravitational field is exactly instantaneous** in this model.

---

## 5) 0PN orbit dynamics: Newtonian potential from a point sink

For a point-like central source:
[
\rho(\mathbf{x})=M,\delta^{(3)}(\mathbf{x}),
]
Poisson gives the usual:
[
\Phi(r)=-\frac{\mu}{r},
\qquad \mu\equiv GM,
]
and therefore
[
a_r(r)=-\frac{d\Phi}{dr}=-\frac{\mu}{r^2}.
]

So: **the toy model reproduces exactly Newtonian central-force orbits at 0PN.**

---

## 6) Retarded scalar sector: why it gives *no* 1PN perihelion precession (key result)

The paper derives a scalar Liénard–Wiechert–like retarded potential for a moving point source using the retarded Green function. But the crucial reduction is:

* For **central-field dynamics** (static source at the origin), the retarded solution collapses back to the Newtonian (-\mu/r) form, i.e.
  [
  \Phi_{\mathrm{ret}}(t,r(t))=-\frac{\mu}{r(t)},\qquad \Phi_{\mathrm{L}}(r)=0.
  ]

**Therefore the scalar lag sector contributes no 1PN perihelion advance:**
[
\Delta\varphi_{\mathrm{scalar}}(a,e)=0.
]

This is the paper’s “null result” that forces the inertia-based mechanism.

---

## 7) The mechanism that *does* produce GR-like 1PN precession: position-dependent inertia

The move is: keep the potential Newtonian, but modify the **kinetic prefactor** in the effective test-body Lagrangian.

### Effective inertia ansatz

[
m_{\mathrm{eff}}(r)=m\left[1+\sigma(r)\right],
\qquad
\sigma(r)=\beta,\frac{\mu}{\cS^2 r}.
]

### Effective Lagrangian used for orbit reduction

[
L=\frac12\left[1+\sigma(r)\right]\left(\dot r^2+r^2\dot\varphi^2\right)-\Phi_{\mathrm{scalar}}(r),
\qquad
\Phi_{\mathrm{scalar}}(r)=-\frac{\mu}{r}.
]

### Resulting perihelion advance per orbit

The final compact result is:
[
\Delta\varphi_{\mathrm{tot}}
\simeq 2\pi\beta,\frac{\mu}{\cS^2,a(1-e^2)}.
]

To match Schwarzschild 1PN:
[
\Delta\varphi_{\mathrm{GR}} = 6\pi,\frac{\mu}{c^2,a(1-e^2)},
]
the paper sets (\cS=c) and fixes:
[
\beta=3.
]

**Interpretation inside the paper:** all the GR-like 1PN precession is being “carried” by the position dependence of the kinetic term (an effective spatial-metric / inertia effect), not by a retarded correction to (\Phi).

---

## 8) Hydrodynamic interpretation of (\beta): the (\kappa)-decomposition (very important “physics meaning”)

The paper wants (\beta=3) to be *derivable* from superfluid/throat microphysics later, so it decomposes:
[
\beta=\kappa_\rho+\kappa_{\mathrm{add}}+\kappa_{\mathrm{PV}}.
]

And then argues for the size of each contribution:

### (i) Density-driven cavitation mass: (\kappa_\rho=1)

Heuristic: defect mass scales like displaced/evacuated volume times ambient density,
[
m\sim \rho_0 V_{\mathrm{cav}},
]
and since (\rho(r)) is perturbed by the gravitational potential (barotropic response), you get a leading-order inertial renormalization (\propto \mu/(\cS^2 r)) with coefficient ~1.

### (ii) Added mass: (\kappa_{\mathrm{add}}=\tfrac12)

Classical hydrodynamics result: a moving spherical cavity in potential flow entrains fluid, adding an inertial term equal to **half** the displaced mass (in the simplest incompressible idealization).

### (iii) Pressure–volume (PV) inertia must supply the rest

Given (\beta=3), the paper concludes:
[
\kappa_{\mathrm{PV}}
\simeq 3 - 1 - \frac12
= \frac32.
]
This is framed as a strong constraint: **the throat’s internal compressible dynamics / PV work must contribute a large effective inertia** if the toy model is to reproduce GR.

---

## 9) What the numerics actually do (implementation-level memory)

The “Numerical experiments” section is mainly sanity/consistency checks of the analytic story via two styles:

### A) Full PDE evolution on a 3D grid

* Evolve (\Phi) via the sourced wave equation (finite differences, Courant-limited timestep).
* Recompute (\Phi_{\mathrm{P}}) each step by FFT Poisson solve.
* Define (\Phi_{\mathrm{L}}=\Phi-\Phi_{\mathrm{P}}) for diagnostics.
* Defects move under (-\nabla\Phi) with a symplectic integrator.
* Static-limit test: hold (\rho) fixed, evolve, and verify (\Phi_{\mathrm{L}}\to 0).

### B) Reduced orbit integration (central potential)

* Scalar-only: integrate orbits under (\Phi=-\mu/r) and confirm null precession.
* Effective 1PN: integrate using the modified kinetic prefactor with (\beta=3) and verify the perihelion advance matches the GR formula.

---

## 10) The “mental checklist” for future work with this paper

If I’m going to work with/extend this model later, here are the *non-negotiable* frozen takeaways:

1. **Field split is structural:** (\Phi=\Phi_{\mathrm{P}}+\Phi_{\mathrm{L}}) with (\Phi_{\mathrm{P}}) enforcing instantaneous 0PN gravity.
2. **Lag sector is retarded but does not generate 1PN precession** for a static central source in the test-mass limit.
3. **GR-like 1PN perihelion precession is put into the inertia sector** via
   [
   \sigma(r)=\beta,\frac{\mu}{\cS^2 r},\quad \beta=3\ (\cS=c).
   ]
4. **(\beta) is not treated as magic:** it’s targeted as a hydrodynamic sum
   [
   \beta=\kappa_\rho+\kappa_{\mathrm{add}}+\kappa_{\mathrm{PV}},
   ]
   with (\kappa_\rho=1), (\kappa_{\mathrm{add}}=1/2), (\kappa_{\mathrm{PV}}=3/2) required.
5. **Numerical solver philosophy:** evolve (\Phi) as a wave, solve Poisson by FFT every step, and treat (\Phi_{\mathrm{L}}) diagnostically as “what’s left.”

===

1pn_optics.tex

## 1) What `1pn_optics.tex` is doing (one sentence)

It builds the **1PN optical / clock sector** of the toy model by treating the vacuum as a **stiff polytropic superfluid** whose density/sound-speed profile around a flux-tube “mass defect” produces a **refractive index (N(r))**; that single (N(r)) is then shown to reproduce **GR light bending + Shapiro delay**, while **gravitational redshift** is obtained via **local-density scaling of rest mass / oscillator frequency**—and the “optical metric alone” mismatch with orbital precession motivates a **bi-metric / dressed-soliton** viewpoint.

---

## 2) Inputs it assumes from Paper 1 (orbital sector)

The optics paper treats these as already established:

* **Effective orbital Lagrangian form** (same structure as Paper 1):
  [
  L
  = \frac{1}{2}\bigl[1 + \sigma(r)\bigr]
  \bigl(\dot r^2 + r^2 \dot\varphi^2\bigr)

  * \Phi_{\mathrm{eff}}(r).
    \tag{`eq:SG-L-eff-again`}
    ]
* (\sigma(r)) is the **inertia/kinetic prefactor** (hydrodynamic “dressing”), and Paper 1 fixed the orbital coefficient to **(\beta=3)** (reviewed in this optics paper’s intro).

---

## 3) Vacuum model: stiff polytrope + hydrostatic balance

### Equation of state (the vacuum “material law”)

[
P = K \rho^n.
\tag{`eq:eos-polytrope`}
]

### Background sound speed

[
c_0^2 \equiv \left.\frac{\partial P}{\partial \rho}\right|_{\rho=\rho_0}
= n K \rho_0^{,n-1}.
\tag{`eq:c0-def`}
]

### Hydrostatic balance around a central mass defect

[
\frac{1}{\rho(r)},\frac{\mathrm{d}P}{\mathrm{d}r}
= \frac{\mathrm{d}\Phi}{\mathrm{d}r}.
\tag{`eq:hydrostatic-balance`}
]

Linearized (weak field, (\rho\simeq\rho_0), (\Phi=-\mu/r), (\mu=GM)):
[
\frac{1}{\rho_0},\frac{\mathrm{d}P}{\mathrm{d}r}
\simeq \frac{\mu}{r^2}.
\tag{`eq:hydrostatic-linear`}
]

### Integrating the pressure perturbation

Define (\Delta P(r)=P(r)-P_0), integrate from (r) to (\infty):
[
\Delta P(r)
= - \frac{\mu \rho_0}{r}
= - \frac{GM \rho_0}{r}.
\tag{`eq:deltaP`}
]

### Density perturbation (via (c_0^2 = \partial P/\partial\rho))

[
\frac{\Delta\rho(r)}{\rho_0}
= \frac{\Delta P(r)}{\rho_0 c_0^2}
= -,\frac{GM}{c_0^2 r}.
\tag{`eq:deltarho`}
]

---

## 4) From density → sound speed → refractive index (N(r))

Sound speed from the EOS:
[
c_s^2(\rho)=\frac{\mathrm{d}P}{\mathrm{d}\rho}=nK\rho^{,n-1}.
]

Linearized fractional change (about (\rho_0)):
[
\frac{\Delta c_s}{c_0}
\simeq \frac{n-1}{2},\frac{\Delta\rho}{\rho_0}.
\tag{`eq:deltacs-rho`}
]

Substitute (\Delta\rho/\rho_0) from above:
[
\frac{\Delta c_s}{c_0}
\simeq -,\frac{n-1}{2},\frac{GM}{c_0^2 r}.
\tag{`eq:deltacs`}
]

Define refractive index as
[
N(r)\equiv \frac{c_0}{c_s(r)}.
\tag{`eq:app-ppn-N` (definition)}
]

Then to leading order:
[
N(r)
\simeq 1 - \frac{\Delta c_s}{c_0}
\simeq 1 + \frac{n-1}{2},\frac{GM}{c_0^2 r}.
\tag{`eq:N-general-n`}
]

### Special choice: (n=5) (the paper’s key fixed value)

[
N_{n=5}(r)
\simeq 1 + 2,\frac{GM}{c_0^2 r}.
\tag{`eq:N-n5`}
In later “GR matching” statements the paper effectively identifies (c_0 \leftrightarrow c) (asymptotic signal speed).

---

## 5) Gravitational lensing from refraction (light bending)

For (n=5), the paper just uses
[
N(r)\simeq 1 + 2,\frac{GM}{c_0^2 r},
\tag{`eq:lensing-N`}
so
[
\ln N(r)\simeq 2,\frac{GM}{c_0^2 r}.
\tag{`eq:lensing-lnN`}
]

Geometric-optics deflection approximation:
[
\Delta\theta
\simeq \int_{-\infty}^{+\infty} \nabla_\perp \ln N\bigl(r(z)\bigr),\mathrm{d}z.
\tag{`eq:lensing-deflection-general`}
]

Carrying out the standard integral gives:
[
\Delta\theta
= \frac{4GM}{b c^2}.
\tag{`eq:lensing-dtheta-result`}
(again using (c_0\to c) in the final GR comparison).

### General (n) result (useful for “what if we change stiffness?”)

Define (\alpha_n\equiv (n-1)/2) (introduced in the appendix), then:
[
\Delta\theta_n
= \frac{2\alpha_n GM}{b c^2}
= \frac{(n-1)GM}{b c^2}.
\tag{`eq:app-dtheta-n`}
So the GR coefficient (4GM/(bc^2)) selects (n=5) in this “pure refraction” construction.

---

## 6) Shapiro time delay from the same (N(r))

Coordinate travel time along a (nearly straight) path is modeled as:
[
t
\simeq \frac{1}{c_0}
\int_{-Z_\mathrm{E}}^{+Z_\mathrm{R}}
\left[1 + 2,\frac{GM}{c_0^2 r(z)}\right]\mathrm{d}z.
\tag{`eq:shapiro-t-integral`}
]

Define (\Delta t=t-t_0), then:
[
\Delta t
\simeq \frac{2GM}{c_0^3}
\int_{-Z_\mathrm{E}}^{+Z_\mathrm{R}}
\frac{\mathrm{d}z}{r(z)}.
\tag{`eq:shapiro-dt-integral`}
Evaluating the log integral yields:
[
\Delta t
\simeq \frac{2GM}{c_0^3},
\ln!\left(\frac{4 r_\mathrm{E} r_\mathrm{R}}{b^2}\right).
\tag{`eq:shapiro-dt-result`}
General (n):
[
\Delta t_n
= \alpha_n \frac{GM}{c^3}
\ln!\left(\frac{4 r_\mathrm{E} r_\mathrm{R}}{b^2}\right).
\tag{`eq:app-dt-n`}
(So lensing + Shapiro both key off the same (\alpha_n=(n-1)/2) in this setup.)

---

## 7) Gravitational redshift: clocks via local density / mass scaling

This paper’s redshift mechanism is **not** “(g_{tt}) time dilation of the optical metric” (since the optical metric takes (g_{tt}^{(\text{opt})}=-1) below). Instead:

1. Use the same density profile:
   [
   \frac{\delta\rho(r)}{\rho_0}
   = -\frac{GM}{r c^2}

* \mathcal{O}!\left(\frac{G^2M^2}{c^4 r^2}\right).
  \tag{`eq:redshift-deltarho`}

2. Assume local rest mass scales with local vacuum density:
   [
   m(r) \propto \rho_{\text{local}}(r)
   \ \Rightarrow
   \frac{\delta m(r)}{m_0}
   = \frac{\delta\rho(r)}{\rho_0}.
   \tag{`eq:redshift-dm`}
   So:
   [
   \frac{\delta m(r)}{m_0}
   = -\frac{GM}{r c^2}+\cdots
   \tag{`eq:redshift-dm-final`}
3. Oscillator frequency scales with mass:
   [
   \omega \propto m
   \ \Rightarrow
   \frac{\delta\omega}{\omega_0}
   = \frac{\delta m}{m_0}.
   \tag{`eq:redshift-omega-m`}
   Final chain:
   [
   \frac{\delta\omega(r)}{\omega_0}
   = \frac{\delta m(r)}{m_0}
   = \frac{\delta\rho(r)}{\rho_0}
   = -\frac{GM}{r c^2}

* \mathcal{O}!\left(\frac{G^2M^2}{c^4 r^2}\right).
  \tag{`eq:redshift-full-chain`}
  So the leading redshift matches the usual GR weak-field scaling.

---

## 8) Metric viewpoint and the “optical metric ≠ matter metric” point

### Optical metric (what light “sees” in geometric optics)

[
\mathrm{d}s^2
= -c^2 \mathrm{d}t^2

* N^2(r),\bigl(\mathrm{d}r^2 + r^2 \mathrm{d}\Omega^2\bigr).
  \tag{`eq:SG-optical-metric`}
  Null geodesics of this metric reproduce the lensing + Shapiro integrals used earlier.

### Post-Newtonian metric form used to connect to orbital dynamics

[
\mathrm{d}s^2
= -\bigl[1 + 2\Phi_{\text{eff}}(r)/c^2 + \dots\bigr] c^2 \mathrm{d}t^2

* \bigl[1 + 2\Psi(r)/c^2 + \dots\bigr]
  \bigl(\mathrm{d}r^2 + r^2 \mathrm{d}\varphi^2\bigr).
  \tag{`eq:SG-metric-PN`}
  The corresponding geodesic Lagrangian (slow-motion expansion):
  [
  L_{\text{geo}}
  = \frac{1}{2}\bigl[1 + 2\Psi(r)/c^2\bigr]
  \bigl(\dot r^2 + r^2 \dot\varphi^2\bigr)

- \Phi_{\text{eff}}(r)

* \mathcal{O}!\left(\frac{v^4}{c^2}\right).
  \tag{`eq:SG-L-geo`}
  Comparing this with the orbital-sector effective Lagrangian implies the key identification:
  [
  1+\sigma(r)\ \longleftrightarrow\ 1+\frac{2\Psi(r)}{c^2}.
  ]
  So (\sigma(r)) is effectively the matter-sector “spatial potential” piece.

### The crucial warning: using the optical metric as the *full* spacetime metric breaks perihelion precession

In the appendix, the paper shows that if you naively interpret the optical metric as the full spacetime metric for massive bodies, you get (in PPN language) (\gamma_{\text{opt}}=2), (\beta_{\text{opt}}=1), leading to:
[
\Delta\varphi_{\text{opt,PN}}
= 10,\frac{\pi GM}{a(1-e^2)c^2},
]
i.e. **coefficient 10 instead of GR’s 6** (overestimate by (5/3)). (See the “Precession from the optical metric” appendix discussion around `subsec:app-optical-precession`.)

**Resolution proposed by the paper:** massive defects are **hydrodynamically dressed solitons**, so matter does *not* follow the bare optical metric; light probes the refractive/acoustic projection, while matter probes an emergent/dressed metric (bi-metric structure).

---

## 9) Fixed key values / “model knobs that got frozen”

* **Polytropic index:** (n=5) is singled out by the optical tests in the “pure refraction” setup.
* **Index coefficient:** for general (n), (\alpha_n=(n-1)/2); for (n=5), (\alpha=2).
* **Refractive index profile (the workhorse):**
  [
  N(r)\simeq 1 + 2\frac{GM}{c^2 r}.
  ]
* **Optical-sector inferred PPN:** matching light deflection implies (\gamma=1) (same as GR), but the paper emphasizes this is being carried by the refractive/spatial part with (g_{tt}^{(\text{opt})}=-1), hence the need for the dressed-matter sector for orbits.

---

## 10) Numerical implementation essentials (ray tracing / time-of-flight)

The appendix gives a concrete algorithm that’s worth remembering:

### Hamiltonian for rays in (N(\mathbf{x}))

Dispersion:
[
\omega(\mathbf{x},\mathbf{k})
= c_s(\mathbf{x}),|\mathbf{k}|
= \frac{c}{N(\mathbf{x})},|\mathbf{k}|.
]
Hamiltonian:
[
H(\mathbf{x},\mathbf{k})
= c_s(\mathbf{x}),|\mathbf{k}|
= \frac{c}{N(\mathbf{x})},|\mathbf{k}|.
]

Ray equations (affine parameter (\lambda)):
[
\dot{\mathbf{x}}
= \frac{\partial H}{\partial \mathbf{k}}
= c_s(\mathbf{x})\frac{\mathbf{k}}{|\mathbf{k}|},
\qquad
\dot{\mathbf{k}}
= -\frac{\partial H}{\partial \mathbf{x}}
= -|\mathbf{k}|,\nabla c_s(\mathbf{x}).
]

Practical note from the paper: initialize with (|\mathbf{k}|=1) and **renormalize occasionally** to control drift.

### Extracting the numerical deflection angle

Compute from incoming/outgoing asymptotic wavevectors:
[
\Delta\theta_{\text{num}}(b)
= \arccos!\left(
\frac{\mathbf{k}*\mathrm{in}\cdot\mathbf{k}*\mathrm{out}}
{|\mathbf{k}*\mathrm{in}||\mathbf{k}*\mathrm{out}|}
\right).
]

Time-of-flight is essentially (t\sim (1/c)\int N,d\ell) (consistent with the Shapiro integral setup in the main text).

===

1pn_spin_and_nbody.tex

## 1) What this paper adds to Papers I (orbits) + II (optics)

It extends the toy model’s **1PN coverage** to:

1. **Spin / gravitomagnetism**: reproduce the weak-field Kerr **(g_{0i})** and the **Lense–Thirring (LT) precession** using a *vortical (vector) flow* carried by a rotating “dyon” defect.

2. **(N)-body 1PN dynamics**: reproduce the **Einstein–Infeld–Hoffmann (EIH)** Lagrangian structure for multiple point masses by showing how the toy model generates:

* the **(v^4)** kinetic corrections,
* the **static nonlinear (“gravity gravitates”) 3-body term** (\sim G^2 m_A m_B m_C/(r_{AB} r_{AC})),
* the **velocity-dependent pairwise terms** (the “EIH cross-tensor” structure).

The punchline: **matching EIH fixes specific hydrodynamic parameters** in the vector/wake sector (notably (\alpha) and (K)).

---

## 2) Inputs it assumes from the earlier papers (reused definitions)

It imports the earlier “dictionary” objects and treats them as established:

* **Effective scalar potential** (central-field limit)
  [
  \Phi_{\text{eff}}(r)\approx -\frac{\mu}{r} + \text{(small 1PN pieces)},\qquad \mu\equiv GM.
  ]

* **Position-dependent inertia / kinetic prefactor**
  [
  \sigma(r)=\beta,\frac{\mu}{\mathcal{S}^2 r},
  ]
  with Paper I’s GR perihelion match giving **(\beta=3)** (and typically (\mathcal{S}\to c) when comparing to GR).

* **Optics side** (Paper II) is referenced mainly as “already matched,” not rederived here.

---

## 3) Spin sector: how LT precession is produced

### 3.1 GR target (what must be matched)

Weak-field Kerr gives a gravitomagnetic precession vector for a gyroscope:
[
\boldsymbol{\Omega}*{\text{LT}}(\mathbf{r})
= \frac{G}{c^2 r^3}\Big(3(\mathbf{J}!\cdot!\hat{\mathbf{r}})\hat{\mathbf{r}}-\mathbf{J}\Big),
]
with the associated weak-field Kerr off-diagonal metric structure summarized by
[
g*{0i}^{\rm (GR)} \propto \epsilon_{ijk}\frac{J^j x^k}{r^3}.
]

### 3.2 Toy-model metric dictionary (what “spin” couples through)

The paper uses an “acoustic / effective” weak-field form where the **off-diagonal piece is sourced by a vector potential tied to bulk flow**:
[
ds^2 \supset -\frac{4}{c^3},\mathbf{A}*{\rm eff}\cdot d\mathbf{x},dt,
\quad \text{with}\quad \mathbf{A}*{\rm eff}\propto \rho_0,\mathbf{v}.
]
(Exact proportionality is a calibration constant.)

So: **find a far-field vortical flow (\mathbf{v}(\mathbf{r}))** whose induced (\mathbf{A}*{\rm eff}) matches Kerr’s (g*{0i}).

### 3.3 Required far-field flow (“dyon” vortex)

A pure “rotating sphere” backflow dies too fast:
[
\mathbf{v}_{\text{backflow}}(\mathbf{r})
\sim \frac{R^3}{r^3},\boldsymbol{\Omega}\times\mathbf{r}
\quad (r\gg R),
]
so the paper instead uses a **dyon**: a defect that carries both a **scalar sink (mass)** and a **vortex structure (spin)**.

The key far-field vortical profile is taken as (spherical coordinates)
[
v_\phi(r,\theta)=\frac{D}{r^2}\sin\theta\qquad (r\gg a),
]
which has the right “dipole-like” scaling to reproduce the Kerr-like (g_{0i}) and hence (\Omega_{\rm LT}\propto 1/r^3).

### 3.4 Spin calibration (fixed coefficient)

Matching the toy-model vector potential to GR fixes:
[
D=\frac{4G}{c^2},J.
]
This is the “spin → flow strength” calibration you carry forward.

**What to remember:**
If you want spin effects (frame dragging, LT precession) in the toy model, you attach a **vortical dyon component** with far-field (v_\phi\propto J/r^2), normalized by (D=4GJ/c^2).

---

## 4) (N)-body 1PN sector: EIH structure and how the toy model generates it

### 4.1 EIH decomposition (what must be reproduced)

The paper treats the standard 1PN structure as:
[
L_{\text{EIH}} = L_{\text{N}} + \frac{1}{c^2}L_{1\rm PN} + \mathcal{O}(c^{-4}),
]
and conceptually splits (L_{1\rm PN}) into three pieces:

1. **Kinetic** (v^4)-type corrections ((L_{\rm kin}))
2. **Static nonlinear** 3-body term ((L_{\rm stat}))
3. **Velocity-dependent pairwise** term ((L_{\rm vec}))

The “cross-tensor” part of the EIH velocity-dependent interaction between bodies (A,B) is written (this is the calibration target):
[
L_{\text{vec}}^{(AB)}
=\frac{Gm_A m_B}{r_{AB}}
\left[
\frac{3}{2}(v_A^2+v_B^2)
-\frac{7}{2}\mathbf{v}_A\cdot\mathbf{v}*B
-\frac{1}{2}(\mathbf{v}*A\cdot\mathbf{n}*{AB})(\mathbf{v}*B\cdot\mathbf{n}*{AB})
\right],
]
with (\mathbf{n}*{AB}=(\mathbf{x}_A-\mathbf{x}*B)/r*{AB}).

### 4.2 Where each piece comes from in the toy model

#### (i) (L_{\rm kin}): relativistic expansion of defect kinetic energy

This comes from treating defects as moving in the effective metric already fixed by Papers I/II (the same “matter metric” behind (\sigma(r))).

#### (ii) (L_{\rm stat}): “gravity gravitates” via mass depending on local density/pressure

Key ansatz:
[
m_A(\mathbf{x}*A)=m*{A,0}\left[1+\kappa_\rho,\frac{\Phi_{\rm loc}(\mathbf{x}_A)}{c^2}

* \mathcal{O}!\left(\frac{\Phi^2}{c^4}\right)\right],
  ]
  with local potential from other bodies approximately
  [
  \Phi_{\rm loc}(\mathbf{x}*A)\simeq -\sum*{C\neq A}\frac{G m_{C,0}}{r_{AC}}.
  ]

Plugging this into the Newtonian pair potential energy
[
V_{\rm N}= -\frac{1}{2}\sum_{A\neq B}\frac{G,m_A(\mathbf{x}*A)m_B(\mathbf{x}*B)}{r*{AB}}
]
generates an emergent **3-body** static term of the schematic form
[
V*{\rm stat}^{(3)} \sim \kappa_\rho,\frac{G^2}{c^2}\sum_{A,B,C}^{\prime}\frac{m_A m_B m_C}{r_{AB}r_{AC}},
]
up to symmetrization conventions (the paper rewrites it into the standard EIH-looking symmetric prime-sum form and denotes the overall coefficient as (C_{\rm stat}\propto \kappa_\rho)).

**What to remember:** (\kappa_\rho) is the knob that controls the EIH **static nonlinear** coefficient.

#### (iii) (L_{\rm vec}): overlap energy of moving dyon flows

This is the key “new mechanism” for the EIH velocity-dependent tensor structure:

Define a vector interaction energy between two moving dyons as an overlap integral of their induced flow fields:
[
V_{\text{vec}}^{(AB)}=\rho_0\int d^3x,\mathbf{u}_A(\mathbf{x})\cdot\mathbf{u}_B(\mathbf{x}).
]

To compute this cleanly, the paper uses a **Fourier-space wake ansatz** for the translational flow of a moving dyon. It decomposes the response into transverse/longitudinal pieces with projectors:
[
\mathcal{P}_L^{ij}(\mathbf{k})=\frac{k^i k^j}{k^2},\qquad
\mathcal{P}_T^{ij}(\mathbf{k})=\delta^{ij}-\frac{k^i k^j}{k^2}.
]

and takes the moving wake field to be (core formula):
[
\mathbf{u}_{\text{trans}}(\mathbf{k};\mathbf{v})
================================================

i\frac{K}{k},\mathcal{P}_T(\mathbf{k}),\mathbf{v}
;+;
i\frac{K a_H}{k^2}(\mathbf{k}\times\mathbf{v})
;+;
i\frac{K\alpha}{k^3},\mathbf{k}(\mathbf{k}\cdot\mathbf{v}).
]

Interpretation of parameters:

* **(K)**: overall coupling (strength of wake/flow response)
* **(\alpha)**: longitudinal (compressible) wake admixture
* **(a_H)**: helical transverse component (set to 0 for the minimal EIH match)

Evaluating the overlap integral yields the EIH-like pair interaction
[
V_{\text{vec}}^{(AB)}
=\frac{Gm_A m_B}{c^2 r_{AB}}
\left[
C_\parallel,\mathbf{v}_A\cdot\mathbf{v}_B
+
C_L,(\mathbf{v}*A\cdot\mathbf{n}*{AB})(\mathbf{v}*B\cdot\mathbf{n}*{AB})
+\cdots
\right].
]

The derived coefficients are (crucial closed-form result):
[
C_\parallel(\alpha,a_H)=K\pi^2(-1+a_H^2-\alpha^2),
\qquad
C_L(\alpha,a_H)=K\pi^2(-1+a_H^2+\alpha^2).
]

### 4.3 Exact EIH cross-tensor match fixes (\alpha) and (K)

To match the EIH cross coefficients
[
C_\parallel \stackrel{!}{=} -\frac{7}{2},\qquad
C_L \stackrel{!}{=} -\frac{1}{2},
]
the paper selects the **minimal** choice (a_H=0) and solves, obtaining:
[
\alpha^2=\frac{3}{4},
\qquad
K=\frac{2}{\pi^2}.
]

This is one of the most important “frozen” outputs of the paper:
**EIH matching forces a specific transverse/longitudinal mix in the wake response**.

---

## 5) The “carry-forward constants” from this paper

If you only memorize a few things from `1pn_spin_and_nbody.tex`, make it these:

### Spin / gravitomagnetism

* Required far-field spin flow:
  [
  v_\phi(r,\theta)=\frac{D}{r^2}\sin\theta,
  ]
* Calibration:
  [
  D=\frac{4G}{c^2},J,
  ]
* This reproduces Kerr-like (g_{0i}) and LT precession.

### (N)-body EIH match (vector sector)

* Wake ansatz parameters:

  * **set (a_H=0)** (minimal clean match),
  * **(\alpha^2=3/4)**,
  * **(K=2/\pi^2)**,
* Then:
  [
  C_\parallel=-\frac{7}{2},\qquad C_L=-\frac{1}{2},
  ]
  which are exactly the hard-to-get EIH cross-tensor coefficients.

### Static nonlinear sector

* Defect mass responds to local potential:
  [
  m_A(\mathbf{x}*A)=m*{A,0}\Bigl[1+\kappa_\rho,\Phi_{\rm loc}(\mathbf{x}_A)/c^2+\cdots\Bigr],
  ]
* This generates the EIH-style **3-body static** term with an overall coefficient proportional to **(\kappa_\rho)**.

---

## 6) Practical “mental checklist” for future work

When you later build sims / derivations that include spin and multiple bodies, this paper’s logic is:

1. Use the **scalar sector** for Newtonian + whatever you already use for orbital 1PN.
2. For spin:

   * treat each spinning defect as a **dyon with a calibrated vortex dipole** ((D=4GJ/c^2)).
3. For (N)-body 1PN:

   * generate the static nonlinearity from (m(\Phi_{\rm loc})) (tied to (\kappa_\rho)),
   * generate EIH velocity-dependent terms from the **flow-overlap** energy,
   * enforce EIH cross-tensor structure by fixing (\alpha^2=3/4), (K=2/\pi^2), (a_H=0).

===

em_fields.tex

This paper is doing for **electromagnetism** what the earlier 1PN papers did for **gravity**: it defines a **hydrodynamic dictionary** (superfluid variables → EM potentials/fields), shows that **Maxwell’s equations** emerge from an **acoustic wave operator**, and shows the **Lorentz force** emerges from **Magnus + pressure/enthalpy forces** on a “dyon” throat defect. It also ties EM strength and the **EM/Gravity hierarchy** to **defect geometry** ((a,L,\Gamma)) and a **cavity-mode energy minimization** that selects the preferred **aspect ratio** (L/a).

---

## 1) The core “hydrodynamic → EM” dictionary (the thing to cache)

### Fluid variables

* (\rho(\mathbf{x},t)): density, (\rho_0) background density
* (\mathbf{v}(\mathbf{x},t)): bulk velocity
* (h(\mathbf{x},t)): specific enthalpy (Bernoulli-like scalar)
* (\boldsymbol{\omega} \equiv \nabla\times\mathbf{v}): vorticity
* (c) is identified with the **acoustic speed** (c_s) in the weak-field limit.

### EM potentials (key definitions)

Introduce a normalization constant (\lambda) and define
[
\phi_{\mathrm{EM}} \equiv \lambda \left(h + \frac{1}{2} v^2\right),
\qquad
\mathbf{A} \equiv \lambda,\mathbf{v}.
\tag{EM dictionary}
]

### EM fields

[
\mathbf{B} \equiv \nabla \times \mathbf{A}
\quad\Rightarrow\quad
\mathbf{B} = \lambda,\boldsymbol{\omega},
]
[
\mathbf{E} \equiv -\nabla \phi_{\mathrm{EM}} - \partial_t \mathbf{A}
\quad\Rightarrow\quad
\mathbf{E} = -\lambda\left[\nabla!\left(h+\frac12 v^2\right)+\partial_t\mathbf{v}\right].
]

**Immediate payoff:** the homogeneous Maxwell equations are identities:
[
\nabla\cdot\mathbf{B}=0,\qquad \nabla\times\mathbf{E}=-\partial_t\mathbf{B}.
]

---

## 2) Gauge invariance (and what it means in the fluid)

The standard gauge transform
[
\phi_{\mathrm{EM}}'=\phi_{\mathrm{EM}}-\partial_t\chi,\qquad
\mathbf{A}'=\mathbf{A}+\nabla\chi
]
leaves (\mathbf{E},\mathbf{B}) invariant.

Fluid interpretation (the key idea): gauge freedom corresponds to **reparameterizing the velocity potential** in the irrotational sector ((\mathbf{v}=\nabla\Phi)) and compensating in the Bernoulli combination (h+\tfrac12 v^2).

---

## 3) Inhomogeneous Maxwell equations from an acoustic wave equation

### Flat-space wave operator

In the weak-field/flat limit the acoustic d’Alembertian becomes
[
\Box ;\rightarrow; -\frac{1}{c^2}\partial_t^2 + \nabla^2.
]

### Postulate: sourced wave equation for the 4-potential

Assemble the four-potential
[
A^\mu=\left(\frac{\phi_{\mathrm{EM}}}{c},,\mathbf{A}\right),
]
and promote the vacuum acoustic equation to
[
\Box A^\mu = -\mu_0 J^\mu,
\qquad
\partial_\mu A^\mu = 0 \quad \text{(Lorenz gauge)}.
]

In components:
[
\Box \phi_{\mathrm{EM}} = -\mu_0 c^2 \rho_e,
\qquad
\Box \mathbf{A} = -\mu_0 \mathbf{J},
]
with the flat Lorenz gauge
[
\frac{1}{c^2}\partial_t \phi_{\mathrm{EM}} + \nabla\cdot\mathbf{A} = 0.
]

From these, the paper derives the inhomogeneous Maxwell equations:
[
\nabla \cdot \mathbf{E} = \frac{\rho_e}{\epsilon_0},
\qquad
\nabla\times\mathbf{B}

* \frac{1}{c^2}\partial_t \mathbf{E}
  = \mu_0 \mathbf{J},
  ]
  together with current conservation (\partial_\mu J^\mu=0), and the standard constraint
  [
  \epsilon_0 \mu_0 = \frac{1}{c^2}.
  ]

---

## 4) What counts as “charge” and “current” (dyon defects)

A moving defect worldline (\mathbf{x}_d(t)) is modeled as a point source:
[
\rho_e(\mathbf{x},t)=q,\delta^3(\mathbf{x}-\mathbf{x}_d(t)),
\qquad
\mathbf{J}(\mathbf{x},t)=\rho_e,\mathbf{u}_d(t).
]

### Effective charge from geometry + circulation

Charge is not fundamental here; it is **geometric/topological**:
[
q ;=; \kappa_q,\rho_0 \Gamma A
= \kappa_q,\rho_0 \pi a^2 \Gamma.
]

* (a): throat radius
* (\Gamma): circulation (treated as a **topologically conserved integer label** / “charge state”)
* (\kappa_q): dimensionless calibration constant.

### Coulomb limit (breathing-mode argument)

In the static limit:
[
\nabla^2 \phi_{\mathrm{EM}} = -\mu_0 c^2 q,\delta^3(\mathbf{x}),
]
so outside the core
[
\phi_{\mathrm{EM}}(r)=\frac{\mu_0 c^2 q}{4\pi r}
=\frac{q}{4\pi\epsilon_0 r}.
]

---

## 5) Lorentz force = Magnus + “electric-like” pressure/enthalpy force

This is the other headline result: the paper shows the standard Lorentz force emerges as the *macroscopic* force on a throat defect.

### Magnus force (magnetic term)

For a vortex line (per unit length):
[
\mathbf{f}*M
= \rho_0 \Gamma,\hat{\mathbf{t}} \times (\mathbf{u} - \mathbf{v}*{\text{fluid}}).
]
For a straight throat along (\hat{\mathbf{z}}), the (u)-dependent part is
[
\mathbf{f}_{M,u}=\rho_0\Gamma,\hat{\mathbf{z}}\times\mathbf{u}.
]

The paper defines the magnetic Lorentz force per unit length as
[
\mathbf{f}*{L,\mathrm{mag}}
= \frac{q}{L},\mathbf{u}\times\mathbf{B}.
]
Matching (\mathbf{f}*{L,\mathrm{mag}}=\mathbf{f}_{M,u}) fixes how the *coarse-grained* throat field (\mathbf{B}) encodes circulation/geometry (in that simplified matching it yields an effective (B_0) proportional to (-L/(\kappa_q\pi a^2))).

### Electric-like force (pressure gradients + unsteady acceleration)

A small defect in a slowly varying background feels a net pressure-gradient force approximated as
[
\mathbf{F}*{\mathrm{press}}\simeq -V*{\mathrm{eff}},\nabla p(\rho)
\simeq -\rho_0 V_{\mathrm{eff}},\nabla h,
]
and breathing/unsteady effects contribute a term (\propto \rho_0 V_{\mathrm{eff}}\partial_t\mathbf{v}). Combining:
[
\mathbf{F}*E
\simeq \rho_0 V*{\mathrm{eff}}
\left[-\nabla\left(h+\frac12 v^2\right)-\partial_t\mathbf{v}\right]
= \frac{\rho_0 V_{\mathrm{eff}}}{\lambda},\mathbf{E}.
]
This motivates the identification
[
q \equiv \frac{\rho_0 V_{\mathrm{eff}}}{\lambda}
\quad\Rightarrow\quad
\mathbf{F}_E = q,\mathbf{E}.
]

### Net force statement (final packaged result)

The paper summarizes the dyon dynamics as
[
\mathbf{F}_{\mathrm{net}}
= m_G \mathbf{g}

* q\left(\mathbf{E} + \mathbf{u}_d \times \mathbf{B}\right)
* \mathbf{F}*{\mathrm{self}} + \mathbf{F}*{\mathrm{PN}},
  ]
  and thus
  [
  m_G \frac{d\mathbf{u}_d}{dt}
  = m_G \mathbf{g}
* q\left(\mathbf{E} + \mathbf{u}_d \times \mathbf{B}\right)
* \cdots
  ]

---

## 6) Defect “mass” and the EM/Gravity hierarchy

### Gravitational mass from throat volume

[
m_G ;=; \kappa_m,\rho_0 V(a,L)
= \kappa_m,\rho_0 \pi a^2 L.
]

### Force ratio (clean diagnostic)

[
\frac{F_{\mathrm{elec}}}{F_{\mathrm{grav}}}
= \frac{1}{4\pi \epsilon_0 G} \frac{q^2}{m_G^2}.
]

Using the cavity-selected geometry (below), the paper writes the key scaling as
[
\frac{F_{\mathrm{elec}}}{F_{\mathrm{grav}}}
= \left[\frac{1}{4\pi \epsilon_0 G}
\frac{\kappa_q^2}{\kappa_m^2}
\frac{x_{01}^2}{2\pi^2}\right]
\frac{\Gamma^2}{a^2}.
]

**Interpretation you’ll reuse later:** once the aspect ratio is fixed, the EM/Gravity hierarchy is controlled mainly by the **circulation-to-radius** ratio (\Gamma/a) (up to the calibration constants (\kappa_q,\kappa_m)).

---

## 7) Resonant cavity enthalpy minimization (why (L/a) gets selected)

The paper uses a “bubble/cavity mode” energy for the fundamental Bessel/standing-wave mode (with (x_{01}) = first zero of (J_0)):

[
E_{\mathrm{mode}}(a,L)
= \alpha_{\mathrm{cav}},\rho_0 c_s^2,\pi a^2 L
\left(\frac{x_{01}^2}{a^2} + \frac{\pi^2}{L^2}\right)
= \alpha_{\mathrm{cav}},\rho_0 c_s^2 \pi
\left(L x_{01}^2 + a^2 \frac{\pi^2}{L}\right).
]

Minimizing this at fixed “charge content” yields the preferred aspect ratio:
[
\frac{L}{a} = \frac{\sqrt{2},\pi}{x_{01}}
\simeq 1.85.
]

That number shows up repeatedly in the later brane/bulk and hybrid papers.

---

## 8) What’s *actually* free vs. fixed (calibration picture)

This paper is very explicit that the normalization (\lambda), (\epsilon_0), (\mu_0) are **not** arbitrary once you decide how to calibrate:

* (c) is the **acoustic signal speed**.
* (\epsilon_0\mu_0=1/c^2).
* (\epsilon_0) is fixed operationally by matching the **breathing-mode** (Coulomb law).
* (\lambda) is fixed by matching the **electric-force** identification (q=\rho_0 V_{\mathrm{eff}}/\lambda) together with the chosen defect geometry/calibration of (q).
* (\kappa_q,\kappa_m) are the remaining dimensionless “microphysics” knobs capturing geometric/non-ideal corrections.

---

## 9) Minimal cache for future work (what I’d keep in RAM)

If we’re building plasma sims / extending EM later, you mostly need:

1. **Dictionary**
   [
   \phi_{\mathrm{EM}}=\lambda(h+\tfrac12 v^2),\quad \mathbf{A}=\lambda\mathbf{v},
   \quad \mathbf{B}=\lambda(\nabla\times\mathbf{v}),\quad
   \mathbf{E}=-\nabla\phi_{\mathrm{EM}}-\partial_t\mathbf{A}.
   ]

2. **Field equations**
   [
   \Box A^\mu = -\mu_0 J^\mu,\quad \partial_\mu A^\mu=0,
   \quad \epsilon_0\mu_0=1/c^2,
   ]
   giving
   [
   \nabla\cdot\mathbf{E}=\rho_e/\epsilon_0,\quad
   \nabla\times\mathbf{B}-\frac{1}{c^2}\partial_t\mathbf{E}=\mu_0\mathbf{J}.
   ]

3. **Source model**
   [
   \rho_e=q\delta^3(\mathbf{x}-\mathbf{x}_d),\quad \mathbf{J}=\rho_e\mathbf{u}_d,
   \quad q=\kappa_q\rho_0\pi a^2\Gamma,\quad
   m_G=\kappa_m\rho_0\pi a^2 L.
   ]

4. **Lorentz force emergence**
   [
   \mathbf{F}=m_G\mathbf{g}+q(\mathbf{E}+\mathbf{u}\times\mathbf{B})+\cdots
   ]
   with magnetic part from **Magnus** and electric part from **pressure/enthalpy + unsteady acceleration**.

5. **Geometry selection**
   [
   L/a=\sqrt2\pi/x_{01}\simeq 1.85
   ]
   and hierarchy scaling (F_E/F_G\propto \Gamma^2/a^2).

===

brane_bulk_ontology.tex

## 1) What this paper is for

It resolves the **sphere–cylinder tension** between earlier sectors of the toy model:

* **Gravity / 1PN (Papers I–III):** defects must look like **spherical sinks** to reproduce the correct far-field monopole behavior used in orbital dynamics.
* **EM / cavity modes (Paper IV):** stable charged defects require a **cylindrical resonant cavity** (radius (a), length (L)) with Bessel-mode structure, selecting a specific **aspect ratio**.

**Resolution (this paper’s core move):**
Defects are **brane–bulk throats**: a localized “mouth” on the **3D brane** that opens into a **4D superfluid bulk**.

* Far from the mouth on the brane, you see an **almost perfect monopole** (spherical sink).
* Inside the throat (in the bulk direction), you naturally get **cylindrical cavity modes**.

---

## 2) Ontology and geometry you should keep in mind

### Brane–bulk setup

* Bulk is **4D** with coordinates ((x,y,z,w)).
* Observable world is a **3D brane** at (w=0).
* A “particle/defect” is not a 3D object; it’s a **throat** extending into the bulk.

### Throat domain (idealized)

[
\mathcal{T}
\simeq \left{
(\mathbf{x},w),\middle|,
0 \le w \le L,;
\sqrt{x^2+y^2+z^2}\lesssim a
\right}.
\tag{eq:throat_domain}
]

### A concrete “rounded funnel” toy geometry (used for finite-size estimates)

[
R(w) = a\left[1+\frac{1}{2}\exp!\left(-\frac{5w}{a}\right)\right].
\tag{eq:funnel_radius}
]

Interpretation: near (w=0) the mouth is a bit flared; deeper in (w) it approaches radius (\sim a).

---

## 3) Bulk superfluid model (the governing PDE backbone)

### Polytropic equation of state

[
P = K \rho^n.
\tag{eq:poly_eos}
]

### Bulk continuity + Euler (4D)

[
\partial_t \rho + \nabla_4 \cdot (\rho,\mathbf{v}_4)=0,
\tag{eq:bulk_continuity}
]
[
\rho\Big(\partial_t \mathbf{v}_4 + (\mathbf{v}_4\cdot\nabla_4)\mathbf{v}*4\Big)
= -\nabla_4 P + \mathbf{f}*{\mathrm{ext}}.
\tag{eq:bulk_euler}
]

### Enthalpy field and linearization

[
h(\rho)\equiv \int^\rho \frac{dP}{\rho'}
= \frac{nK}{n-1}\rho^{n-1},
\tag{eq:enthalpy_def}
]
[
h \simeq \frac{c_s^2}{\rho_0},\delta\rho.
\tag{eq:enthalpy_perturbation}
]

### 4D acoustic wave equation (linear bulk modes)

[
\partial_t^2 h - c_s^2\nabla_4^2 h = 0,
\qquad \nabla_4^2=\partial_x^2+\partial_y^2+\partial_z^2+\partial_w^2.
\tag{eq:4d_acoustic}
]

This is the “mode engine” used later to show why the *interior* of the throat wants cylindrical Bessel profiles.

---

## 4) Dimensional reduction: how the brane sees a 3D monopole

To connect to the earlier 3D effective descriptions, the paper defines **brane-projected fields** by integrating over (w) with a kernel (K(w)) peaked at (w=0):

[
\rho_{\mathrm{3D}}(\mathbf{x},t)
\equiv \int_{-\infty}^{+\infty}dw;K(w),\rho(\mathbf{x},w,t),
\tag{eq:rho3d_def}
]
[
\mathbf{v}*{\mathrm{3D}}(\mathbf{x},t)
\equiv \frac{1}{N}\int*{-\infty}^{+\infty}dw;K(w),\Pi_{\parallel}\mathbf{v}_4(\mathbf{x},w,t),
\qquad N=\int dw,K(w).
\tag{eq:v3d_def}
]

Then the induced 3D flow is decomposed as
[
\mathbf{v}*{\mathrm{3D}}=\nabla\Phi+\mathbf{v}*T,
\qquad \nabla\cdot \mathbf{v}*T=0,
\tag{eq:v3d_decomposition}
]
leading to an effective 3D continuity + Poisson structure on the brane:
[
\partial_t\rho*{\mathrm{3D}}+\nabla\cdot(\rho*{\mathrm{3D}}\mathbf{v}*{\mathrm{3D}})=0,
\tag{eq:3d_continuity}
]
[
\nabla^2\Phi = 4\pi G_{\mathrm{eff}}\rho_{\mathrm{3D}}.
\tag{eq:3d_poisson}
]

**Key modeling point:** (G_{\mathrm{eff}}) is treated as a calibrated effective constant here (not derived from first principles in this paper).

---

## 5) Far-field structure: “why 1PN gravity sees a sphere”

### Multipole expansion of the brane potential

[
\Phi(r,\theta)\approx -\frac{GM}{r} - G,\frac{Q}{r^3}P_2(\cos\theta)+\cdots.
\tag{eq:phi_multipole}
]

And more generally:
[
\Phi(r,\theta)
= -\frac{GM}{r}\Bigg[1+\sum_{\ell\ge2}\alpha_\ell\left(\frac{a}{r}\right)^\ell P_\ell(\cos\theta)\Bigg].
\tag{eq:phi_multipole_general}
]

### Core scaling claim

Even if the throat geometry has anisotropy (\varepsilon), the angular correction is suppressed as:
[
\Phi(r,\theta)\approx -\frac{GM}{r}\left[1+\mathcal{O}!\left(\varepsilon\frac{a^2}{r^2}\right)P_2(\cos\theta)+\cdots\right].
\tag{eq:phi_farfield_scaling}
]

So: **the brane observer gets a clean monopole at (r\gg a)**, with finite-size “shape” corrections down by ((a/r)^2). This is the bridge between “throat interior is cylindrical” and “far field is spherical.”

### Concrete Gaussian toy model (used to sanity-check scaling)

[
\rho_4(r,\theta,w)=\rho_0 e^{-r^2/a^2}e^{-w^2/L^2}\big[1+\varepsilon P_2(\cos\theta)\big],
\tag{eq:rho4_gaussian}
]
[
\rho_{\mathrm{3D}}(r,\theta)=\int dw,\rho_4
=\sqrt{\pi}L\rho_0 e^{-r^2/a^2}\big[1+\varepsilon P_2(\cos\theta)\big].
\tag{eq:rho3d_gaussian}
]
Total mass:
[
M=\int d^3x,\rho_{\mathrm{3D}}=\pi^2 L a^3\rho_0.
\tag{eq:gaussian_mass}
]
Quadrupole moment scaling:
[
Q_{20}=\frac{3}{10}\pi^2 L a^5\varepsilon\rho_0,
\tag{eq:q20_result}
]
so
[
\frac{Q_{20}}{M}=\frac{3}{10}\varepsilon a^2.
\tag{eq:q_over_m_scaling}
]
which is exactly what drives the ((a/r)^2) suppression in the potential.

### Numerical “rounded funnel” estimate: quadrupole coefficient is small

The abstract reports a representative hard-bounded funnel integration giving a leading quadrupole coefficient:
[
\alpha_2 \approx 1.4\times 10^{-2}
]
for (L/a=2) and a (10%) quadrupolar anisotropy (this is used as evidence that the spherical-sink approximation for 1PN gravity is robust).

**Interpretation for later work:** (\alpha_2) (and higher (\alpha_\ell)) are where **finite-size / 2PN-ish** falsifiable corrections live.

---

## 6) Inside the throat: why the EM sector wants a cylinder (mode separation)

This paper explicitly connects the throat interior to the cylindrical cavity modes from the EM construction.

### Separation ansatz for acoustic enthalpy modes

[
h(t,r,w)=\Re{H(r)W(w)e^{-i\omega t}}.
\tag{eq:separated_ansatz}
]

Resulting separated equations:
[
\nabla_\perp^2 H+k_r^2 H=0,
\qquad
\partial_w^2 W+k_w^2 W=0,
\tag{eq:radial_eq_general, eq:axial_eq_general}
]
with dispersion
[
\omega^2=c_s^2(k_r^2+k_w^2).
\tag{eq:dispersion_relation}
]

### Fundamental cylindrical mode (Bessel × sine standing wave)

The fundamental is written in Bessel form (with (x_{01}) the first zero of (J_0)):
[
h_1(t,r,w)\propto J_0!\left(\frac{x_{01}r}{a}\right)\sin!\left(\frac{\pi w}{L}\right)\cos(\omega t).
\tag{eq:fundamental_mode}
]
with
[
\omega^2=c_s^2\left(\frac{x_{01}^2}{a^2}+\frac{\pi^2}{L^2}\right).
\tag{eq:fundamental_dispersion}
]

### Why a specific aspect ratio is selected: energy minimization at fixed “charge”

The paper uses:
[
\mathcal{E}[h]\sim \int_{\mathcal{T}} d^3x,dw;\rho_0\left[\frac{1}{2c_s^2}(\partial_t h)^2+\frac{1}{2}(\nabla_4 h)^2\right],
\tag{eq:enthalpy_functional}
]
and shows mode scalings of the form
[
\mathcal{E}\propto A^2(k_r^2+k_w^2)\int_{\mathcal{T}}H^2W^2,
\tag{eq:energy_kw}
]
[
\mathcal{Q}\propto A^2\int_{\mathcal{T}}H^2W^2.
\tag{eq:charge_functional}
]
So at fixed (\mathcal{Q}), minimizing (\mathcal{E}) reduces to minimizing ((k_r^2+k_w^2)), yielding the preferred geometry:
[
\frac{L}{a}=\frac{\sqrt{2}\pi}{x_{01}}\approx 1.85.
\tag{eq:aspect_ratio_result}
]

**This is the paper’s key “geometry unification” statement:**

* **Gravity** probes the *mouth* → spherical sink.
* **EM** probes the *throat interior* → cylindrical cavity with a fixed (L/a).

---

## 7) EM ontology update: magnetostatics from brane transverse wakes

A big “cleanup” claim here is that you can get long-range magnetostatics **without requiring bulk vorticity**, by identifying the EM vector potential with the **brane transverse wake** component (\mathbf{v}_T).

### Definitions and equations

[
\mathbf{A}\equiv \kappa_A,\mathbf{v}_T,
\tag{eq:A_from_vT}
]
with Coulomb gauge and a vector Poisson equation
[
\nabla\cdot\mathbf{A}=0,
\qquad
\nabla^2\mathbf{A}=-\kappa_A,\mathcal{Q},\mathbf{u},\delta^{(3)}(\mathbf{r}).
\tag{eq:vector_poisson_A}
]
Solution for a moving point “charge”:
[
\mathbf{A}(\mathbf{r})=\frac{\kappa_A\mathcal{Q}}{4\pi}\frac{\mathbf{u}}{r},
\tag{eq:A_u_over_r}
]
and hence
[
\mathbf{B}=\nabla\times\mathbf{A}
=\frac{\kappa_A\mathcal{Q}}{4\pi}\frac{\mathbf{u}\times\mathbf{r}}{r^3},
\tag{eq:biot_savart_point}
]
i.e. a Biot–Savart-like far field.

**Interpretation:**

* **Charge** is reinterpreted as a **vorticity flux / circulation** property of the throat.
* **Magnetostatics** comes from the correct completion of the wake basis on the brane (the transverse part).

---

## 8) Parameter constraint that shows up again: (\alpha^2=3/4)

This paper also records:
[
\alpha^2=\frac{3}{4}.
\tag{eq:alpha_squared_value}
]

This is the same “wake-mixing” constraint that appears in the EIH-matching story (velocity-dependent (N)-body structure) in the spin/(N)-body paper. Here it’s part of the ontology unification: completing the wake basis consistently pins the same parameter.

---

## 9) The “carry-forward checklist” (what to remember)

If we’re working with this paper later, the minimal cache is:

1. **Defects are brane–bulk throats**: spherical mouth on the brane, cylindrical cavity in the bulk.
2. **Dimensional reduction** defines (\rho_{\mathrm{3D}}) and (\mathbf{v}*{\mathrm{3D}}) by a kernel (K(w)) near (w=0); then (\mathbf{v}*{\mathrm{3D}}=\nabla\Phi+\mathbf{v}*T) and (\nabla^2\Phi=4\pi G*{\mathrm{eff}}\rho_{\mathrm{3D}}).
3. **Far-field geometry corrections are suppressed by ((a/r)^2)**; multipole form encoded by (\alpha_\ell) in
   [
   \Phi=-\frac{GM}{r}\Big[1+\sum_{\ell\ge2}\alpha_\ell(a/r)^\ell P_\ell(\cos\theta)\Big].
   ]
   Representative estimate: (\alpha_2\sim 1.4\times10^{-2}) for a funnel geometry case.
4. **Interior throat modes** are cylindrical Bessel × axial standing waves; minimizing energy at fixed charge selects
   [
   \frac{L}{a}=\frac{\sqrt2\pi}{x_{01}}\approx 1.85.
   ]
5. **Magnetostatics from brane transverse wake:** (\mathbf{A}=\kappa_A\mathbf{v}_T) and ( \mathbf{B}=\nabla\times\mathbf{A}) gives Biot–Savart form for a moving charge.

===

1pn_hybrid.tex

Below is the same “essentials-only / future working memory” extraction for **`1pn_hybrid.tex`**. This one is a *capstone* paper: it **stitches scalar + vector sectors together**, shows a **uniqueness result** (forces **(n=5)** and **(q=1)**), and then pushes into **strong-field structure** (acoustic horizon, photon-sphere-like orbit, shadow/lensing predictions) with the **nonlinear inflow** and **throat impedance** as the control knobs for deviations.

---

## 0) One-paragraph purpose

**Goal:** build a single “hybrid” 1PN effective theory where:

* **Scalar sector (density starvation)** makes defect mass depend on the local potential/density.
* **Vector/optical sector (polytropic + flow)** makes propagation speed vary and introduces flow couplings that generate EIH-like velocity-dependent terms.
* When you expand the combined test-body Lagrangian to 1PN and demand agreement with **GR/EIH**, the paper claims the model itself selects:

  * **polytropic index (n=5)**,
  * **equivalence exponent (q=1)**,
  * and carries forward the already-fixed wake-mixing **(\alpha^2=3/4)**.

Then it treats the **nonlinear inflow** as a physical completion that produces an **acoustic horizon** and **strong-field optical observables** (photon sphere / shadow scale), with deviations governed by **throat size (a)**, aspect ratio (\Lambda=L/a), and **throat impedance**.

---

## 1) Frozen parameter “ledger” (this paper assumes earlier fixes)

Right in the intro, the paper consolidates the series’ parameters into one table:

* (\beta = 3) (Paper I) — total scalar inertia renormalization (perihelion)
* (\kappa_\rho=1), (\kappa_{\rm PV}=3/2), (\kappa_{\rm add}=1/2) (Paper I) — decomposition of (\beta)
* (n=5) (Paper II) — stiff polytrope (optical matching)
* (\alpha^2=3/4) (Paper III) — wake mixing (EIH vector sector)
* (\Lambda_\star=L/a\simeq 1.85) (Paper V) — preferred throat geometry
* (\Gamma\in\mathbb{Z}) (Paper IV) — circulation/charge control, (Q\propto\Gamma)

**Notation note:** it uses (\beta_{\rm idx}) for the *optical index* (1/r) coefficient, to avoid collision with the orbital (\beta).

---

## 2) Thermodynamic vacuum model (core assumptions)

### Polytropic EOS + sound speed

[
P(\rho)=K\rho^n
]
[
\mathcal{S}^2(\rho)=\frac{dP}{d\rho}=Kn\rho^{n-1}
]
and the background is fixed so that
[
\mathcal{S}(\rho_0)=c.
]

Small density contrast:
[
\delta \equiv \frac{\rho-\rho_0}{\rho_0}, \qquad |\delta|\ll 1
]
so
[
\mathcal{S}^2(\rho)=c^2\left[1+(n-1)\delta+\mathcal{O}(\delta^2)\right].
]

### “Equivalence assumption”: defect mass tracks vacuum density

[
M(\rho)\propto \rho^q.
]
Interpretation in-text: (q=1) means mass scales with the feeding fluid density; (q=0) would be a rigid mass unaffected by the medium.

Weak-field density response to potential:
[
\rho(\Phi)=\rho_0[1+\delta(\Phi)],\qquad \delta(\Phi)=\gamma\frac{\Phi}{c^2}+\mathcal{O}!\left(\frac{\Phi^2}{c^4}\right),
]
so
[
M(\Phi)=M_0\left[1+q\gamma\frac{\Phi}{c^2}+\mathcal{O}!\left(\frac{\Phi^2}{c^4}\right)\right].
]

The hybrid “master schematic” Lagrangian is introduced as
[
L = -M(\Phi)c^2\sqrt{1-\frac{v^2}{\mathcal{S}^2(\Phi)}} ;+; L_{\text{int}}[\Phi,\mathbf{v}].
]

---

## 3) Scalar sector (density starvation → variable mass → + sign in (\Phi v^2/c^2))

Central potential:
[
\Phi(r)=-\frac{GM_B}{r}.
]

The scalar-sector depletion statement is made explicitly as
[
\frac{\delta\rho}{\rho_0}\approx \frac{\Phi}{c^2},
]
which then drives
[
M(\Phi)=M_0\left(1+q\frac{\Phi}{c^2}+\mathcal{O}!\left(\frac{\Phi^2}{c^4}\right)\right).
]

Relativistic scalar Lagrangian:
[
L_{\text{sc}}=-M(\Phi)c^2\sqrt{1-\frac{v^2}{c^2}}.
]

Using
[
\sqrt{1-\frac{v^2}{c^2}} = 1-\frac12\frac{v^2}{c^2}-\frac18\frac{v^4}{c^4}+\cdots,
]
the key *hybrid-relevant* output is: **scalar sector contributes a (+) coefficient to (\Phi v^2/c^2)**, summarized later as
[
L_{\text{sc}} \supset +\frac12,q,M_0,\frac{\Phi v^2}{c^2}.
]

Also: the Newtonian coupling from this form is
[
L^{(0)} \supset -M_0 q,\Phi,
]
which is why the Newtonian limit will force (q=1) (next section).

---

## 4) Vector sector (polytrope → (c_s(\Phi)) → refractive index + kinetic sign flip + flow terms)

From the polytrope expansion, in weak field:
[
\frac{c_s^2(\Phi)}{c^2} = 1 + (n-1)\frac{\Phi}{c^2} + \mathcal{O}!\left(\frac{\Phi^2}{c^4}\right),
]
so
[
\frac{c_s(\Phi)}{c} = 1 + \frac12(n-1)\frac{\Phi}{c^2} + \mathcal{O}!\left(\frac{\Phi^2}{c^4}\right).
]

Define the index:
[
N(\Phi)=\frac{c}{c_s(\Phi)}
= 1 - \frac12(n-1)\frac{\Phi}{c^2} + \mathcal{O}!\left(\frac{\Phi^2}{c^4}\right).
]

Optical metric form used:
[
ds_{\text{opt}}^2 = -\frac{c^2}{N^2(\Phi)}dt^2 + N^2(\Phi),d\mathbf{x}^2.
]

Vector kinetic piece:
[
L_{\text{kin,vec}}=-M_0 c^2\sqrt{1-\frac{v^2}{c_s^2(\Phi)}},
]
and its key 1PN contribution is the **opposite sign**:
[
L_{\text{kin,vec}} \supset -\frac12(n-1),M_0,\frac{\Phi v^2}{c^2}.
]

Then there is an explicit **flow coupling** term (used for the EIH-like velocity-dependent structure):
[
L_{\text{flow}}
= M_0,\mathbf{v}\cdot\mathbf{u}(\mathbf{x})
= M_0,\mathbf{v}\cdot\mathbf{u}_{\text{trans}}

* M_0,\mathbf{v}\cdot\mathbf{u}_{\text{vort}}.
  ]

And an (N)-body vector template is written in the EIH style:
[
L_{\text{vec}}^{(AB)}
= \frac{G M_A M_B}{r_{AB}}
\left[\cdots + \text{(terms in }\mathbf{v}_A\cdot\mathbf{v}_B,,
(\mathbf{v}*A\cdot\mathbf{n}*{AB})(\mathbf{v}*B\cdot\mathbf{n}*{AB})\text{)}\right],
]
with the detailed coefficient bookkeeping deferred to the appendices/tables.

---

## 5) The hybrid Lagrangian and the uniqueness result (q=1,; n=5)

The hybrid is defined as:
[
L_{\text{hyb}}=L_{\text{sc}}+L_{\text{kin,vec}}.
]

Expanded to the needed order, the “famous” term is the (\Phi v^2/c^2) coefficient:
[
L_{\text{hyb}}\supset
C_{\text{hyb}}(n,q),M_0,\frac{\Phi v^2}{c^2},
\qquad
C_{\text{hyb}}(n,q)=\frac12 q-\frac12(n-1).
]

### Step 1: Newtonian limit fixes (q)

The paper isolates the Newtonian piece as
[
L_{\text{hyb}}^{(0)}=\frac12 M_0 v^2 - M_0 q,\Phi,
]
so matching the Newtonian gravitational coupling requires
[
q=1.
]

Then
[
C_{\text{hyb}}(n)\equiv C_{\text{hyb}}(n,q=1)=\frac{2-n}{2}.
]

### Step 2: 1PN/EIH matching fixes (n) (and uses the flow sector + (\alpha))

The GR/EIH “target template” is written schematically as
[
L_{\text{EIH}}
= \frac12 M_0 v^2 - M_0 \Phi

* B_{\text{EIH}},M_0\frac{\Phi v^2}{c^2}
* D_{\text{EIH}},M_0\frac{\Phi^2}{c^2}
* \cdots
  ]
  and the matching condition is presented as
  [
  \frac{2-n}{2} + C_{\text{flow}}(n)=B_{\text{EIH}},
  ]
  where (C_{\text{flow}}(n)) is the extra (\Phi v^2/c^2) contribution induced by the flow coupling.

In the “hybrid uniqueness” section, the paper states the solution as:
[
n=5,\qquad q=1,
]
and reiterates the already-fixed wake-mixing constraint
[
\alpha^2=\frac34.
]

**Conceptual memory hook:**

* Scalar density-starvation gives (+\tfrac12 q).
* Vector kinetic gives (-\tfrac12(n-1)).
* Flow terms add (C_{\text{flow}}(n)) and also supply the EIH cross-velocity structures.
* Demanding the GR/EIH pattern forces **(q=1)** and **(n=5)** (and uses **(\alpha^2=3/4)** for the wake basis).

---

## 6) Strong-field completion: nonlinear inflow → acoustic horizon

The strong-field “engine” is a steady, spherically symmetric inflow (u(r)) feeding the throat.

Mass flux conservation:
[
\dot{M}_{\text{flux}} = 4\pi r^2 \rho(r) u(r).
]

Radial Euler equation (effective potential included):
[
u\frac{du}{dr} = -\frac{1}{\rho}\frac{dP}{dr} - \frac{d\Phi}{dr}.
]

Combined into the standard transonic form (plus a throat/finite-size forcing term in this paper):
[
\left(u^2-c_s^2\right)\frac{du}{dr}
= u\left[\frac{2c_s^2}{r}-\frac{d\Phi_{\text{eff}}}{dr}\right]
\quad\text{(and in the horizon section: }+\mathcal{F}_{\text{throat}}(r;a,n)\text{)}.
]

**Acoustic horizon definition:**
[
r=r_H \quad\text{such that}\quad |u(r_H)|=c_s(r_H).
]

Transonic “regularity” conditions are written as:
[
u^2(r_H)=c_s^2(r_H),\qquad
\frac{2c_s^2(r_H)}{r_H}
-\left.\frac{d\Phi}{dr}\right|*{r_H}
+\mathcal{F}*{\text{throat}}(r_H;a,n)=0.
]

Scaling statement (used later for phenomenology):
[
r_H \sim \kappa,\frac{GM}{c^2},\qquad r_H \gtrsim a.
]

---

## 7) Photon sphere, lensing, and “shadow” scale from the optical metric

The spherical optical metric is written as:
[
ds_{\text{opt}}^2
= -\frac{c^2}{N^2(r)}dt^2 + N^2(r)\left(dr^2+r^2 d\Omega^2\right).
]

Weak-field index profile:
[
N(r)=1-\frac12(n-1)\frac{\Phi(r)}{c^2}+\cdots,
\qquad \Phi(r)\simeq -\frac{GM}{r},
]
and for the hybrid-selected (n=5):
[
N(r)\simeq 1+2\frac{GM}{rc^2}.
]

Ray equation used repeatedly:
[
\frac{d\phi}{dr}
= \frac{b}{r^2}\frac{N^2(r)}{\sqrt{N^2(r)-b^2/r^2}}.
]

Photon-sphere-like orbit conditions are expressed in terms of the effective potential:
[
V_{\text{eff}}(r_{\text{ph}};b_{\text{ph}})=0,\qquad
\partial_r V_{\text{eff}}(r_{\text{ph}};b_{\text{ph}})=0,
]
equivalently:
[
\frac{b_{\text{ph}}^2}{r_{\text{ph}}^2}=N^2(r_{\text{ph}}),\qquad
\frac{d}{dr}\left[\frac{N^2(r)}{r^2}\right]*{r=r*{\text{ph}}}=0,
]
which implies the compact slope condition:
[
\frac{N'(r_{\text{ph}})}{N(r_{\text{ph}})}=\frac{1}{r_{\text{ph}}}.
]

For a pure (1/r) index model (N(r)=1+\beta_{\rm idx}GM/(rc^2)), it gives the scaling estimate:
[
r_{\text{ph}}\simeq \frac{3\beta_{\rm idx}GM}{c^2},\qquad
b_{\text{ph}}\simeq \sqrt{3},r_{\text{ph}}.
]

But the paper emphasizes the **full nonlinear inflow index (N(r)=c/c_s(r))** deviates from (1/r) inside/near (r_H), so it summarizes the numerical outcome as:
[
r_{\text{ph}}\simeq \kappa_{\text{ph}}\frac{GM}{c^2},
\qquad
b_{\text{ph}}\simeq \eta_{\text{ph}}\frac{GM}{c^2},
]
with (\kappa_{\text{ph}}) “slightly larger than 3” (Schwarzschild) for the (n=5) inflow solutions they compute.

Bending angle expression used in the numerical appendix:
[
\Delta\phi(b)
=2\int_{r_{\min}}^\infty
\frac{b}{r^2}\frac{N^2(r)}{\sqrt{N^2(r)-b^2/r^2}},dr-\pi,
\qquad
N^2(r_{\min})=\frac{b^2}{r_{\min}^2}.
]

---

## 8) Throat impedance + mass–radius scaling (where deviations come from)

Throat impedance is defined as:
[
Z_{\text{th}}\equiv \frac{\Delta P}{\dot{M}_{\text{flux}}}.
]

A geometric mass scaling is written as:
[
M \sim \rho_0 a^3,\mathcal{F}(\Lambda),\qquad \Lambda\equiv \frac{L}{a}.
]

Define the logarithmic mass–radius exponent:
[
k \equiv \frac{d\log M}{d\log a}.
]

The paper summarizes the nonlinear solutions by an approximate power-law:
[
M \simeq \mu_*,a^{k_*},\qquad k_*\lesssim 1,
]
and therefore the horizon radius scales as:
[
r_H \simeq \kappa_*\frac{GM}{c^2}
\simeq \kappa_*\frac{G\mu_*}{c^2},a^{k_*}.
]

**Interpretation that matters for later use:**
Once you’re past 1PN matching, **strong-field differences** (shadow size, photon-sphere shift, etc.) are controlled by:

* (a) and (\Lambda=L/a),
* the impedance (Z_{\text{th}}),
* and the effective exponent (k_*) (how mass grows with throat size).

---

## 9) “Predictions / tests” section: the paper’s preferred deviation parameters

This section is intentionally schematic (it’s about how to *look for failures*).

Key dimensionless knobs it defines:

* Aspect ratio:
  [
  \Lambda=\frac{L}{a}.
  ]
* Finite-size deviation scaling:
  [
  \epsilon_{\text{fs}}\sim \left(\frac{a}{r_H}\right)^p.
  ]
* Compactness:
  [
  \mathcal{C}\equiv \frac{GM}{r_H c^2}\sim \frac{1}{\kappa}.
  ]

The idea is: in regimes where (a/r_H) isn’t tiny (or impedance is large), you should expect systematic departures from GR-like strong-field relations even if 1PN is nailed.

---

## 10) What I’d want “cached” to work with this paper later

If we’re going to build on `1pn_hybrid.tex`, the minimum memory set is:

1. **Hybrid matching logic:** scalar gives (+\tfrac12 q) in (\Phi v^2/c^2), vector kinetic gives (-\tfrac12(n-1)), flow adds (C_{\text{flow}}(n)).
   Enforcing Newtonian + EIH selects:
   [
   q=1,\quad n=5,\quad \alpha^2=\frac34.
   ]

2. **Optics dictionary:** (N(r)=c/c_s(r)); for (n=5) weak field (N\simeq 1+2GM/(rc^2)), but near/inside (r_H) you must use the **full nonlinear inflow** (c_s(r)).

3. **Strong-field structure comes from inflow:** acoustic horizon (r_H) where (|u|=c_s); photon sphere is determined by (N'/N=1/r) using the same (N(r)).

4. **Deviations are not “free”:** they’re controlled by **throat impedance** and **mass–radius scaling** (M\simeq \mu_* a^{k_*}) with (k_*\lesssim 1).