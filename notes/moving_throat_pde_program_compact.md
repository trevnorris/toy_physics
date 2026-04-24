# Moving-Throat PDE Program Compact — V2 Full

## 0. Purpose

This document is the **V2-full compact working master** for the moving-throat
PDE program. It is meant for AI consumption: a new session should be able to
work across the parent field theory, moving-throat PDE scaffold, reduced
support/mouth/core machinery, grouped real `P2` finish line, same-charge and
relaxed/open-system extensions, and the new V2 branch-realization audit without
re-reading the full derivation logs.

This version deliberately keeps the broad mathematical coverage of the original
compact master while integrating the V2 audit corrections. The earlier V2
rewrite was useful as a short field guide, but it was too aggressive as the sole
memory document. This V2-full version restores the older math cards and branch
compiler material, then adds the audit-corrected parent-action, gauge, boundary,
target-surface, branch-freeze, and solver-pipeline layers.

The document is therefore organized as a **single compact ledger** for future
work on:

- the exact parent `4+1` GNLS/Maxwell theory and projection/reduction hooks,
- the corrected status of the moving wall as either an effective closure or a
  promoted throat action,
- the Maxwell gauge-localization safety patch,
- the open finite-radius throat boundary protocol,
- finite-throat support, D/N AC impedance, and mouth DtN data,
- wall/BdG/Maxwell/mixed Schur complements and stability gates,
- Family-1 mouth/core and co-evolving mouth/core branches,
- grouped real `P2` conservative and outgoing quadrupole finish lines,
- coherent local-kernel monomials, similarity orbit, free-quintuple search, and
  selected-support/orbit-lock packets,
- relaxed/open-system leakage, work, non-rigid mouth, export, heat-partition,
  and material-screening companions,
- the reduced 5PN / actual-branch finish-line packet,
- the V2-19 through V2-23 executable branch-realization protocol and first
  target-blind open-throat residual extraction.

This is **not** a paper draft. It is not a claim that the full nonlinear
moving-throat PDE has been solved. It is a compact program ledger whose job is
to preserve the exact/reduced/effective/open status of every ingredient while
keeping enough formulas available for future derivation work.

The extension is integrated thematically rather than as a new stage-by-stage
appendix. Stage numbers remain useful provenance markers, but the document is
organized around the actual theorem structure that AI or a human derivation
session needs to carry forward.

---

## 1. Reading Rules and Status Firewall

### 1.1 Claim-status tags

Every major statement in this document should be read under one of these tags.

- **Exact**: follows directly from the declared action, exact definitions, or exact algebra.
- **Exact Within Closure**: exact inside a stated reduced closure family, branch, or hierarchy.
- **Reduced / Controlled Reduction**: follows only after a stated ansatz, low-frequency reduction, projection, or branch restriction.
- **Effective Closure**: a physically motivated closure choice that is not yet a unique theorem of the completed PDE.
- **Numerically Located**: defined exactly but realized at values currently obtained by numerical solve rather than closed form.
- **Open**: still depends on the actual completed moving-throat PDE branch.

### 1.2 Non-negotiable notation firewall

The following separations are structural.

1. Electric charge is carried by
   \[
   \eta_Q,\qquad q_\star,\qquad q_{\rm eff},
   \]
   not by circulation.
2. The historical gravity-side bare `q=1` is the mass-dressing coefficient
   \[
   \kappa_\rho=1,
   \]
   not electric charge.
3. Grouped labels `20/21/22` denote grouped real `P2` lanes, not spacetime indices.
4. The mixed channels
   \[
   A_w,\qquad J^w,\qquad F_{\mu w},\qquad E_w,\qquad C_a
   \]
   are suppressed only in the strict far-field brane reduction. They remain part
   of the microscopic ontology and are required for the honest outgoing bridge.

### 1.3 Present theorem-status summary

Current best reading of the program:

- the parent `4+1` field-theory block is fixed at the exact declared-action level,
- the moving-throat geometry lift and reduced wall/support/outgoing program are explicit and stable enough for theorem work inside the carried closures,
- the grouped real `P2` Packet-A finish line is no longer an eight-slot residual ledger at theorem order; on the natural isotropic point-particle branch it has collapsed to the single outgoing-normalization scalar
  \[
  \chi_Q-1,
  \]
- the one-port same-charge audit no longer leaves the mixed sector as a generic placeholder: the static bundle creates no new long-range attractive law, and the linear dynamic bundle creates no new kernel class beyond resonant dispersive enhancement of the already-short-range families,
- the actual weak-axisymmetric same-charge bottleneck is now the transported prefactor packet
  \[
  (\Delta_{\rm norm},a_{P_0},b_{P_0}),
  \]
  equivalently the pair
  \[
  (\Delta_{\rm norm},\Xi_1),
  \qquad
  \Xi_1=\frac{P_1}{P_0},
  \]
  on the weak-axisymmetric line
  \[
  b_{P_0}=3a_{P_0},
  \]
- on the rigid-mouth actual branch the post-static same-charge problem is already diagonal in the physical logarithmic chart
  \[
  (U,V)=\left(
  \ln\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2},
  \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}
  \right),
  \]
  and orbit lock is the Cartesian codimension-two condition `U=V=0`,
- the passive/outgoing support-selection ambiguity has collapsed to the exact lowest-twin slice
  \[
  \rho_\alpha=\frac43,
  \qquad
  \zeta_{\rm req}=\frac13,
  \]
  with actual placement then reduced to one coordinate on the selected twin-support curve,
- the post-225 barrier companion is now a codimension-three relaxed branch around the Stage-225 slice, with exact recovery map
  \[
  \ell_w=f_U=a=b=0,
  \]
  and with the same same-charge long-range verdict still frozen,
- on that relaxed branch the stationary lowered front end is carried by
  \[
  V_{\rm eff}^{(230)}(r)
  =
  V_{\rm short}^{(1p)}(r)
  -
  \lambda_L S_{\rm leak}(r)
  -
  \lambda_W\mathcal W_w^{\rm sess}(r)
  -
  \Delta E_{UV}(r)
  -
  \mathcal M_{\sigma,\rm red}(r),
  \]
  while the dynamic continuation already has explicit peak, threshold, turning-point,
  WKB, and Goldilocks compilers,
- the microscopic active-leg export law is now super-Ohmic,
  \[
  K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5,
  \]
  with an exact event-safe half-plane and a channel-resolved vacuum/lattice
  heat-partition law,
- the materials companion is now organized by one calibration dictionary
  \[
  (t_*,\lambda_{\rm phys},E_*,\Upsilon_{\rm lat})
  \]
  and four explicit screening ratios
  \[
  \Pi_{\rm ep},\ \Pi_\chi,\ \Pi_k,\ \Pi_T,
  \]
- the full reduced back end is still the exact four-scalar packet
  \[
  \Delta_{\rm full}=(\Delta_Q,q_{\rm tr},q_{\rm nt},q_\eta),
  \]
  together with its explicit graph-error realization form,
- the local free-quintuple mixed-ray search is still finite through support-cardinality `<=5`,
- and the main remaining theorem gap is now actual branch realization of either
  the strict selected support/orbit/outgoing packet or, if the relaxed barrier
  companion is being used, the realized leakage/non-rigid/source/export/calibration
  packet rather than more algebraic compression.

So the active bottleneck is no longer “derive another compiler.” It is whether
the actual moving-throat branch returns the isotropic one-pole grouped-`P2`
carrier, the canonical outgoing normalization, an admissible weak-axisymmetric
prefactor slope `\Xi_1`, the rigid-mouth Cartesian orbit-lock packet, the
selected twin-support placement state, and—if one needs the back-end repair
audit rather than the front-end branch test—the local Hessian packets required
by the now-completed search sieve.


### 1.4 V2 audit overlay: status corrections and branch-realization discipline

The V2 audit is an **overlay** on the older compact master. It does not delete
the legacy compilers. It corrects their status and makes the branch-realization
workflow executable.

#### Parent-action status correction

The current strict parent action is GNLS plus localized Maxwell, with geometry
entering through the confinement argument:
\[
S_{\rm current}=S_\psi[\psi,A;\Sigma]+S_{\rm EM}[A].
\]
In that action, \(\Sigma/R\) is parent-level as a **confinement-coupling
argument**, but it is not yet parent-level as an autonomous dynamical throat
field. Varying \(V_{\rm conf}(\mathbf X;\Sigma)\) supplies a wall force/source;
it does not by itself produce a wall PDE with \(\eta_{tt}\),
\(\partial_w(T_w\eta_w)\), or \(\Delta_{S^2}\eta\).

A parent-complete moving-throat statement must instead be
\[
\boxed{
S_{\rm total}=S_\psi[\psi,A,\Sigma]+S_{\rm EM}[A]+S_\Sigma[\Sigma;\mathcal C_\Sigma].
}
\]
At quadratic order around a stationary branch,
\[
S_\Sigma\to S_\eta^{(2)}+O(\eta^3).
\]
Until this promotion is made, the distributed wall PDE should be labeled an
**effective linear wall closure**, not an already-derived strict parent field.

#### Maxwell gauge-localization correction

The unweighted Lorenz gauge-fixing term is safe as a formal bulk gauge device,
but unsafe for naive noncompact zero-mode action reduction. Two safe readings
are allowed:

1. keep the unweighted term in the bulk, impose \(\partial\!\cdot A=0\) before
   zero-mode reduction, and choose any convenient 3+1 gauge fixing afterward;
2. localize the gauge fixing with finite \(H_{\rm loc}(w)\), preferably
   \(H_{\rm loc}=Z\), so that
   \[
   \mathcal L_{\rm gf}^{\rm loc}
   =-\frac{H_{\rm loc}(w)}{2\xi\mu_0}(\partial\!\cdot A)^2.
   \]

The general gauge-fixed equation is
\[
\boxed{
\partial_M(ZF^{MN})+\frac{1}{\xi}\partial^N(H\,\partial\!\cdot A)=\mu_0J^N.
}
\]
For \(H=Z\), the zero-mode gauge-fixed action is finite and \(\xi_4=\xi\). The
mixed observables \(E_w\) and \(C_a\) remain exact gauge invariants and are not
artifacts of the gauge choice.

#### Open-throat endpoint correction

The actual branch-realization geometry is an **open finite-radius conduit**:
\[
\boxed{R(0)=a,
\qquad R(L)>0.}
\]
The old hard-cap statement \(R(L)=0\) should be treated as an obsolete toy
idealization unless explicitly declared as such. The D/N half-ladder survives
as an AC support-coordinate impedance/reflection result:
\[
q(0)=0,
\qquad T_wq_w(L,\omega)+Y_L(\omega)q(L,\omega)=0,
\]
with the Neumann end recovered when \(Y_L\to0\). DC/background throughput exits
through the open finite-radius junction and is tracked by leakage/work terms.

#### Branch-freeze rule

Before target residuals are evaluated, freeze:

- parent action and gauge convention,
- wall/interface action status,
- open-exit boundary protocol,
- projection/source map,
- support family,
- mode/port list,
- stability/passivity gates,
- extraction formulas,
- solver/output schema.

No post-residual refit is allowed. A stable open branch that fails the target
packet is falsification data, not a pipeline error.

#### V2 target-realization status

The symbolic grouped-`P2`/outgoing-quadrupole compiler is no longer the main
gap. The active V2 PDE question is whether a frozen, target-blind, open-throat
branch can output stable non-dark full-bundle data satisfying the isotropic
residual packet
\[
R_{\rm pole}=R_{\rm norm}=R_{P2}=R_{P4}=0,
\]
and, when included, the tail gate
\[
R_{\rm tail}=0.
\]
The first V2-23 reduced open-throat branch is a useful negative control: it
passes open/stability/outgoing-transfer gates but fails the target packet.

---

## 2. Parent Theory and Exact Bulk Equations

### 2.1 Arena, coordinates, and indices

**Status:** `Exact`

The fundamental arena is a `(4+1)`-dimensional spacetime with coordinates
\[
x^M=(t,x,y,z,w),
\qquad
M,N\in\{0,1,2,3,4\}.
\]

Bulk spatial coordinates are
\[
\mathbf X=(x,y,z,w)\in\mathbb R^4,
\]
while brane spatial coordinates are
\[
\mathbf x=(x,y,z)\in\mathbb R^3.
\]

Brane indices are
\[
\mu,\nu\in\{0,1,2,3\},
\]
bulk spatial indices are
\[
i,j\in\{x,y,z,w\},
\]
and brane spatial indices are
\[
a,b\in\{x,y,z\}.
\]

The bulk metric is
\[
\eta_{MN}=\mathrm{diag}(-1,+1,+1,+1,+1),
\]
with d'Alembertian
\[
\Box_5=-\partial_t^2+\nabla_3^2+\partial_w^2.
\]

### 2.2 Core fields and corrected dynamical status

**Status:** `Exact` for \(\psi,A_M\); `Exact confinement argument / Effective or promoted dynamical field` for \(\Sigma/R\)

The strict parent dynamical fields already present in the declared action are
\[
\psi(\mathbf X,t),\qquad A_M(x).
\]
The throat geometry enters the exact parent matter sector through the
confinement argument
\[
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
\]
where
\[
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r\in S^2.
\]

The shape field \(R(\Omega,w,t)\) is therefore already exact as the argument of
\(V_{\rm conf}(\mathbf X;\Sigma)\). It becomes an autonomous parent dynamical
field only if a throat action \(S_\Sigma\), or its quadratic approximation
\(S_\eta^{(2)}\), is included in \(S_{\rm total}\). Without that promotion, the
moving-wall PDE is an effective closure sourced by the exact parent fields.

The bulk density is
\[
\rho=|\psi|^2.
\]

The localized gauge field is
\[
A_M=(A_0,A_i),
\qquad
F_{MN}=\partial_MA_N-\partial_NA_M.
\]

The old finite-dimensional throat variables
\[
a(t),\qquad L(t)
\]
survive only as collective moments of the moving-throat field, not as the
fundamental geometry variables.

### 2.3 Charge ontology

**Status:** `Exact`

The corrected electric-charge bookkeeping is
\[
\eta_Q=\pm 1,
\qquad
q_\star=\eta_Q e_\star,
\qquad
e_\star>0.
\]

After canonical zero-mode brane normalization,
\[
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_\star}{\sqrt{Z_{\rm int}}},
\qquad
Z_{\rm int}=\int_{-\infty}^{+\infty} Z(w)\,dw.
\]

The firewall is:

- electric-charge sign is carried by \(\eta_Q\),
- circulation belongs to the magnetic/vortical sector,
- the historical gravity-side `q=1` is really \(\kappa_\rho=1\).

### 2.4 Exact parent action and promoted-throat option

**Status:** `Exact` for GNLS + localized Maxwell; `Effective Closure` for wall dynamics unless \(S_\Sigma\) is included

The strict parent theory currently fixed by the program is
\[
S_{\rm current}=\int dt\,d^4X\,(\mathcal L_\psi+\mathcal L_{\rm EM}),
\]
with geometry encoded through the confinement coupling
\[
V_{\rm conf}(\mathbf X;\Sigma).
\]

This action gives exact matter and Maxwell equations plus a wall force from
\(V_{\rm conf}\). It does **not** by itself give autonomous wall inertia or
stiffness. The parent-complete moving-throat action, when promoted, is
\[
\boxed{
S_{\rm total}=S_\psi[\psi,A,\Sigma]+S_{\rm EM}[A]+S_\Sigma[\Sigma;\mathcal C_\Sigma].
}
\]
The quadratic wall action used later should be read as
\[
S_\Sigma\to S_\eta^{(2)}+O(\eta^3)
\]
unless and until a nonlinear \(S_\Sigma\) is declared.

#### 2.4.1 Matter sector: gauged GNLS

\[
\boxed{
\mathcal L_\psi
=
\frac{i\hbar}{2}\left(\psi^*D_t\psi-\psi D_t\psi^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
-V_{\rm conf}(\mathbf X;\Sigma)\,\rho
-U(\rho).
}
\]

The covariant derivatives are
\[
D_t\psi=\partial_t\psi+\frac{i q_\star}{\hbar}A_0\psi,
\qquad
D_i\psi=\partial_i\psi-\frac{i q_\star}{\hbar}A_i\psi.
\]

The frozen stiff-polytropic equation of state is
\[
P(\rho)=K\rho^5,
\qquad
U(\rho)=\frac{K}{4}\rho^5,
\qquad
h(\rho)=\frac{dU}{d\rho}=\frac{5K}{4}\rho^4.
\]

The corresponding bulk sound speed is
\[
c_s^2(\rho)=\frac{1}{m}\frac{dP}{d\rho}=\frac{5K}{m}\rho^4.
\]

#### 2.4.2 Localized Maxwell sector and gauge-fixing status

The localized Maxwell kinetic term is
\[
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}.
\]
The historically displayed gauge-fixing term is unweighted:
\[
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2.
\]
So the current compact action may still be written as
\[
\boxed{
\mathcal L_{\rm EM}
=
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
-A_MJ_{\rm ext}^M.
}
\]
with the following safety rule: the unweighted Lorenz term is a bulk gauge
fixing device and should not be naively integrated over a noncompact zero mode.
For parent-clean gauge-fixed PDE work, use the weighted form
\[
\boxed{
\mathcal L_{\rm gf}
=-\frac{H(w)}{2\xi\mu_0}(\partial\!\cdot A)^2,
\qquad H_{\rm int}=\int H(w)dw<\infty,
}
\]
with the preferred structural choice \(H=Z\). Then \(\xi_4=\xi Z_{\rm int}/H_{\rm int}\), and for \(H=Z\), \(\xi_4=\xi\).

The full source entering Maxwell's equation is
\[
J_{\rm tot}^M=J_\psi^M+J_{\rm ext}^M.
\]

Important bookkeeping rule:
varying the covariant matter action already generates the dynamical matter
current \(J_\psi^M\), so it must not be double-counted in the explicit Maxwell
source term.

### 2.5 Exact field equations already fixed by the parent theory

#### 2.5.1 Gauged GNLS equation

**Status:** `Exact`

Variation with respect to \(\psi^\ast\) gives
\[
\boxed{
i\hbar D_t\psi
=
\left[
-\frac{\hbar^2}{2m}D_iD_i
+V_{\rm conf}(\mathbf X;\Sigma)
+h(\rho)
\right]\psi.
}
\]

#### 2.5.2 Exact current and continuity

**Status:** `Exact`

The bulk number current is
\[
\boxed{
j^i=\frac{\hbar}{m}\,\Im(\psi^\ast D_i\psi).
}
\]

Exact continuity is
\[
\boxed{
\partial_t\rho+\partial_i j^i=0.
}
\]

Where \(\rho>0\), define the bulk velocity by
\[
j^i=\rho v^i.
\]

#### 2.5.3 Exact localized Maxwell equation

**Status:** `Exact`, with gauge-fixing convention explicit

With a general gauge-fixing weight \(H(w)\), variation with respect to \(A_N\)
gives
\[
\boxed{
\partial_M\!\left(Z(w)F^{MN}\right)
+\frac{1}{\xi}\partial^N\!\left(H(w)\,\partial\!\cdot\!A\right)
=
\mu_0 J_{\rm tot}^N.
}
\]
The legacy displayed equation is the special case \(H=1\). The localized
parent-clean option is \(H=Z\).

The Bianchi identities are
\[
\partial_{[L}F_{MN]}=0.
\]

Because \(ZF^{MN}\) is antisymmetric,
\[
\partial_N\partial_M(ZF^{MN})=0,
\]
so the divergence consistency equation is
\[
\boxed{
\frac{1}{\xi}\Box_5\big(H\,\partial\!\cdot A\big)=\mu_0\partial_NJ_{\rm tot}^N.
}
\]
Gauge invariance requires current conservation:
\[
\partial_MJ_{\rm tot}^M=0.
\]

#### 2.5.4 Madelung rewrite and Euler-like form

**Status:** `Exact`

Write
\[
\psi=\sqrt{\rho}\,e^{i\theta}.
\]

The gauge-invariant bulk velocity is
\[
\boxed{
v_i=\frac{\hbar}{m}\partial_i\theta-\frac{q_\star}{m}A_i.
}
\]

The quantum potential is
\[
\boxed{
Q(\rho)=-\frac{\hbar^2}{2m}\frac{\nabla_4^2\sqrt{\rho}}{\sqrt{\rho}}.
}
\]

The exact Euler-like equation is
\[
\boxed{
m(\partial_t+v_j\partial_j)v_i
=
q_\star(E_i+v_jB_{ij})
-\partial_i\!\left(V_{\rm conf}+h(\rho)+Q(\rho)\right).
}
\]

Here
\[
E_i=-\partial_tA_i-\partial_iA_0,
\qquad
B_{ij}=F_{ij}.
\]

#### 2.5.5 Exact vorticity–gauge identity

**Status:** `Exact`

Define the bulk vorticity 2-form
\[
\Omega_{ij}=\partial_i v_j-\partial_j v_i.
\]

Away from phase singularities,
\[
\boxed{
\Omega_{ij}=-\frac{q_\star}{m}F_{ij}.
}
\]

So circulation belongs to the magnetic/vortical sector rather than the electric-charge dictionary.

#### 2.5.6 Exact mixed-sector gauge invariants

**Status:** `Exact`

The mixed fields
\[
E_w=F_{w0}=-\partial_tA_w-\partial_wA_0,
\]
\[
C_a=F_{aw}=\partial_aA_w-\partial_wA_a,
\]
are exact gauge invariants under
\[
A_0\to A_0-\partial_t\chi,
\qquad
A_a\to A_a+\partial_a\chi,
\qquad
A_w\to A_w+\partial_w\chi.
\]

These mixed channels are suppressed only in the strict far-field zero-mode brane
reduction. They remain part of the microscopic ontology and are essential to the
honest outgoing bridge.

#### 2.5.7 Cold-start projection and zero-mode reduction hooks

**Status:** `Reduced / Controlled Reduction`

For a normalized brane weight \(W(w)\),
\[
\int W(w)\,dw=1,
\]
the exact projected brane observables are
\[
\rho_{\rm brane}(\mathbf x,t)=\int W(w)\rho(\mathbf x,w,t)\,dw,
\]
\[
\mathbf j_{\rm brane}(\mathbf x,t)=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\qquad
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane}.
\]

Projected continuity is exact:
\[
\boxed{
\partial_t\rho_{\rm brane}+\nabla_3\cdot \mathbf j_{\rm brane}=S_{\rm leak},
}
\]
with leakage term
\[
\boxed{
S_{\rm leak}
=
-\left[Wj^w\right]_{-\infty}^{+\infty}
+\int W'(w)j^w\,dw.
}
\]

Under the Helmholtz split
\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,
\qquad
\nabla_3\cdot\mathbf v_T=0,
\]
the exact longitudinal identity is
\[
\boxed{
\rho_{\rm brane}\,\nabla_3^2\varphi
=
S_{\rm leak}
-\partial_t\rho_{\rm brane}
-(\nabla_3\rho_{\rm brane})\cdot(\nabla_3\varphi+\mathbf v_T).
}
\]

In the quasi-static longitudinal regime this becomes the brane Poisson hook for
\(\varphi\).

Under the controlled far-field zero-mode assumptions
\[
A_w\approx 0,
\qquad
\partial_w A_\mu\approx 0,
\qquad
J^w\approx 0,
\qquad
F_{\mu w}\approx 0,
\]
integration over \(w\) gives the effective brane Maxwell sector
\[
\boxed{
\partial_\mu F^{\mu\nu}=\mu_0^{\rm eff}J_{\rm eff}^\nu,
\qquad
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
}
\]

This is a controlled reduction, not a denial of the mixed-core structure. It is
the short cold-start bridge between the exact parent theory and the reduced brane
language used later in the document.

Gauge-fixing caution: if the five-dimensional gauge-fixing term is unweighted
\((H=1)\), impose Lorenz gauge before zero-mode reduction and choose any desired
finite 3+1 gauge fixing afterward. If \(H=Z\), the zero-mode gauge-fixed action
reduces directly with the same \(Z_{\rm int}\) normalization as the kinetic
term. In both readings, the gauge-invariant mixed-core fields are suppressed
only by the far-field ansatz, not removed from the microscopic theory.


### 2.6 What is fixed and what is not at this level

**Fixed at this level**

- the strict parent field content \(\psi,A_M\),
- \(\Sigma/R\) as a confinement-coupling argument,
- the exact gauged GNLS plus localized Maxwell action,
- the exact bulk equations and exact mixed-sector observables,
- the corrected charge ontology,
- the Maxwell gauge-localization safety rule.

**Not fixed at this level**

- autonomous throat dynamics unless \(S_\Sigma\) is promoted,
- the detailed moving-throat branch geometry,
- the reduced wall/support/outgoing hierarchy,
- the actual outgoing quadrupole normalization branch,
- the full PDE realization of the reduced coherent closures.

---

## 3. Moving-Throat Geometry and Throat-Mode Decomposition

### 3.1 Moving throat as a level-set / shape-field lift

**Status:** `Effective Closure`

The smallest moving-throat lift currently used by the program is the hybrid
level-set / shape-field representation
\[
\boxed{
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
}
\]
with
\[
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r\in S^2.
\]

The finite throat surface is
\[
\Sigma(\mathbf X,t)=0.
\]

Sign convention:

- exterior region: \(\Sigma>0\),
- interior/support region: \(\Sigma<0\).

This lift is not yet claimed as a unique theorem of the parent PDE.
It is the smallest effective geometry choice that:

- keeps the `4+1` ontology intact,
- makes the mouth \(S^2\) multipoles explicit,
- and lets the old collective variables emerge as moments rather than as
  fundamental geometry degrees of freedom.

### 3.2 Reference stationary throat and recovery of \((a,L)\)

**Status:** `Exact Within Closure`, with V2 open-endpoint correction

The reference stationary throat is
\[
\Sigma_0(\mathbf X)=r-R_0(w),
\]
with:

- mouth at \(w=0\),
- finite depth \(0\le w\le L_0\),
- mouth radius
  \[
  a_0=R_0(0),
  \]
- open finite-radius exit
  \[
  \boxed{R_0(L_0)>0.}
  \]

The old hard-cap phrase \(R_0(L_0)=0\) is no longer the physical branch-realization
geometry. It may be used only as an explicitly declared toy cap. The V2 branch
geometry is an open conduit: DC/background transport can exit the finite-radius
junction, while AC support coordinates may still reflect through an impedance
mismatch.

The old collective variables are recovered as geometry moments:
\[
a(t)=\frac{1}{4\pi}\int_{S^2}R(\Omega,0,t)\,d\Omega.
\]

If \(L(\Omega,t)\) denotes the chosen open-exit section of the conduit, then
\[
L(t)=\frac{1}{4\pi}\int_{S^2}L(\Omega,t)\,d\Omega.
\]

So the old \((a,L)\) closure survives as the lowest geometry moments of the
distributed open-throat field.

### 3.3 Promoted confinement coupling

**Status:** `Exact Within Closure`

The old confinement potential is promoted to a moving-surface coupling by
\[
\boxed{
V_{\rm conf}(\mathbf X;\Sigma)=V_{\rm wall}\!\left(\frac{\Sigma(\mathbf X,t)}{\ell_c}\right),
}
\]
with \(\ell_c\) a wall-thickness scale.

Linearize around the reference throat by writing
\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]
Then
\[
\delta\Sigma=-\eta,
\]
so the direct wall-bulk coupling is
\[
\boxed{
\delta V_{\rm conf}
=
-\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
}
\]

This is the basic linear source through which the moving wall drives the matter
and gauge sectors in the reduced hierarchy.

### 3.4 Harmonic decomposition on the mouth sphere

**Status:** `Exact Within Closure`

Expand the wall displacement in real spherical harmonics:
\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+
\sum_{m\in P_2({\rm real})} q_{2m}(w,t)\,Y_{2m}^{\rm real}(\Omega)
+
\eta_{\ge 3}(\Omega,w,t).
}
\]

The grouped real `P2` set is
\[
\{20,\ 21c,\ 21s,\ 22c,\ 22s\}.
\]

With
\[
Y_{00}=\frac{1}{2\sqrt\pi},
\]
the physical mouth-average shift \(\delta a\) and the normalized monopole
coefficient are related by
\[
\boxed{
q_{00}(0,t)=2\sqrt\pi\,\delta a(t).
}
\]

A useful axisymmetric split is
\[
\eta_0(w,t)=\alpha_a(w)\,\delta a(t)+\alpha_L(w)\,\delta L(t)+g(w,t),
\]
where \(g(w,t)\) is the residual axisymmetric geometry lane orthogonal to the
collective \((a,L)\) directions.

So the geometry sector separates into:

- collective \(l=0\) throat motion \((\delta a,\delta L)\),
- residual axisymmetric geometry,
- grouped real `l=2` lanes,
- higher \(l\ge 3\) lanes.

### 3.5 Why the grouped real `P2` bundle matters

**Status:** `Reduced / Controlled Reduction`

The grouped real `P2` sector is the first nontrivial harmonic family beyond the
monopole. In the moving-throat program it is the literal geometry/support
realization of the conservative grouped quadrupole payload already exposed by
the earlier reduced PN hierarchy.

It is therefore not decorative harmonic bookkeeping.
It is the first throat-localized multipole bundle that the full PDE must supply
honestly.

### 3.6 Geometry-sector status at this point

**Fixed within the present geometry lift**

- level-set / shape-field throat definition,
- recovery of \((a,L)\) as collective moments,
- promoted confinement coupling,
- mouth-sphere harmonic decomposition.

**Still open**

- the fully realized moving-throat branch geometry of the completed PDE,
- whether the true branch preserves the same reduced harmonic separation cleanly
  beyond the present closure.

---

## 4. Reduced Wall and Finite-Throat Support Engine

### 4.1 Finite-throat D/N support branch

**Status:** `Exact Within Closure`

On the finite throat interval
\[
z\in[0,L],
\]
the minimal internal support equation is
\[
\psi''+k^2\psi=0.
\]

The selected finite-throat support branch imposes:

- Dirichlet at the mouth,
  \[
  \boxed{\psi(0)=0,}
  \]
- Neumann at the bottom,
  \[
  \boxed{\psi'(L)=0.}
  \]

The eigenmodes are
\[
\psi_j(z)=A_j\sin(k_j z),
\qquad
k_j=\frac{\pi}{L}\left(j+\frac12\right),
\qquad
j=0,1,2,\dots
\]
with frequencies
\[
\boxed{
\omega_j=c_s k_j
=
\frac{\pi c_s}{L}\left(j+\frac12\right).
}
\]

This is the exact trapped-support ladder used downstream inside the chosen
finite-throat D/N closure. V2 status: this ladder is an AC support-coordinate
impedance/reflection result, not evidence for a hard geometric cap. The open
finite-radius conduit remains \(R(L)>0\); DC/background throughput exits and is
tracked by leakage/work terms.

### 4.2 Mouth DtN operator and trapped round-trip closure

**Status:** `Exact Within Closure`

For prescribed mouth datum \(\psi_m\), the finite-throat D/N branch gives the
exact mouth derivative
\[
\boxed{
Z_{00}(\omega)
=
-\frac{\omega}{c_s}\tan\!\left(\frac{\omega L}{c_s}\right).
}
\]

Its poles are exactly the D/N ladder above.

The scalar round-trip factor is
\[
R_{\rm rt}=r_0r_Le^{2ikL}.
\]

For the D/N branch,
\[
r_D=-1,
\qquad
r_N=+1,
\]
so on the exact D/N ladder
\[
R_{\rm rt}=1,
\qquad
\phi_0\equiv 0 \pmod{2\pi}.
\]

This is the trapped-support closure currently carried into the reduced program.
With the open-exit patch, the Neumann end is equivalently the low-impedance AC
limit of
\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0,
\]
while the physical throat exit remains finite-radius and open.

### 4.3 Minimal distributed wall action

**Status:** `Effective Closure`

The minimal passive distributed wall action presently used is the quadratic
effective closure
\(S_\eta^{(2)}\). It is parent-complete only if promoted as the leading
expansion of a throat action \(S_\Sigma\). In its current status, it supplies a
consistent linear wall PDE and branch data, but its coefficients are not to be
post-hoc refit. The action is
\[
\boxed{
S_\eta^{(2)}
=
\frac12\int dt\,dw\,d\Omega\,\sqrt{\gamma_0}
\left[
\mu_\eta(w)(\partial_t\eta)^2
-T_w(w)(\partial_w\eta)^2
-T_\Omega(w)\,\eta(-\Delta_{S^2})\eta
-K_\eta(w)\eta^2
\right].
}
\]

Here:

- \(\mu_\eta(w)\) is the effective wall inertia density,
- \(T_w(w)\) is axial wall stiffness,
- \(T_\Omega(w)\) is angular stiffness on the mouth sphere,
- \(K_\eta(w)\) is a local restoring potential.

These are fixed effective constitutive functions of the chosen reference throat.
They are not to be refit stage by stage to rescue downstream normalization.

From this point onward the wall amplitudes are written in a densitized
one-dimensional convention: after integrating over the reference sphere, the
remaining surface weight is absorbed into the effective axial coefficients and
modal amplitudes.

The resulting modal operator is
\[
\mu_\eta q_{lm,tt}
-\partial_w(T_w\partial_w q_{lm})
+\bigl[K_\eta+l(l+1)T_\Omega\bigr]q_{lm}
=
S_{lm}^{(\psi,A)}+f_{lm}^{\rm ext}.
\]

So the scalar lane \(l=0\) and grouped real `P2` lane \(l=2\) are already split
before any additional matter/gauge closure is imposed.

### 4.4 Axisymmetric reduction back to \((a,L)\)

**Status:** `Reduced / Controlled Reduction`

Using the two-mode axisymmetric truncation
\[
\eta_0(w,t)=2\sqrt\pi\,[\alpha_a(w)\delta a(t)+\alpha_L(w)\delta L(t)],
\]
the distributed wall theory reduces to
\[
L_{\rm red}^{(0)}
=
\frac12 M_{AB}\dot Q^A\dot Q^B
-\frac12 K_{AB}Q^AQ^B,
\qquad
Q^A=(\delta a,\delta L),
}
\]
with
\[
M_{AB}=4\pi\int dw\,\mu_\eta\,\alpha_A\alpha_B,
\]
\[
K_{AB}=4\pi\int dw\,[T_w\alpha_A'\alpha_B'+K_0\alpha_A\alpha_B].
\]

This is the conservative reduction back to the old \((a,L)\) matrix system.
So the distributed wall is a lift of the old closure, not a replacement that
forgets it.

### 4.5 Grouped real `P2` reduction and isotropic degeneracy

**Status:** `Reduced / Controlled Reduction`

For one grouped real quadrupole component,
\[
\eta_{2m}(\Omega,w,t)=\beta_2(w)q_{2m}(t)Y_{2m}^{\rm real}(\Omega),
\]
the reduced one-mode Lagrangian is
\[
L_{2m}=\frac12 M_2\dot q_{2m}^2-\frac12 K_2 q_{2m}^2,
\]
with
\[
M_2=\int dw\,\mu_\eta\beta_2^2,
\]
\[
K_2=\int dw\,[T_w(\beta_2')^2+(K_\eta+6T_\Omega)\beta_2^2].
\]

Before symmetry breaking or matter/gauge coupling is turned on, the grouped real
`P2` channels are degenerate on the isotropic reference throat.

That is the microscopic reason the grouped quadrupole block appears as a
degenerate bundle before additional anisotropy is introduced.

### 4.6 Wall/support-engine status at this point

**Fixed within the present reduced engine**

- finite-throat D/N support ladder,
- exact D/N mouth DtN operator,
- trapped round-trip closure on that branch,
- minimal distributed wall action,
- conservative reduction back to \((a,L)\),
- grouped real `P2` reduction and isotropic degeneracy.

**Still open**

- the full coupled matter/gauge renormalization of these reduced lanes,
- the outgoing odd response and quadrupole-normalization branch,
- the actual branch realized by the completed moving-throat PDE.

### 4.7 Linearized coupled bulk/interface skeleton

**Status:** `Reduced / Controlled Reduction`

Take a stationary reference solution
\[
\psi(\mathbf X,t)=e^{-i\mu_0 t/\hbar}\psi_0(\mathbf X),
\qquad
A_M(\mathbf X,t)=A_{M0}(\mathbf X),
\qquad
R(\Omega,w,t)=R_0(w).
\]

At the reduced linearized level:

- the matter sector is BdG-like,
- the geometry sector enters through \(\delta V_{\rm conf}\),
- and the mixed channels
  \[
  \delta A_w,\qquad \delta J^w,\qquad \delta F_{\mu w}
  \]
  remain active.

The linearized matter sector is schematically
\[
i\hbar \partial_t
\begin{bmatrix}
\delta\psi\\
\delta\psi^\ast
\end{bmatrix}
=
L_{\rm BdG}
\begin{bmatrix}
\delta\psi\\
\delta\psi^\ast
\end{bmatrix}
+
C_A[\delta A_M]
+
C_\eta[\eta],
\]
with wall source entering through
\[
\delta V_{\rm conf}
=
-\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c}\,\eta.
\]

The linearized Maxwell sector is
\[
\boxed{
\partial_M\!\left(Z(w)\delta F^{MN}\right)
+\frac{1}{\xi}\partial^N(\partial\!\cdot\!\delta A)
=
\mu_0\,\delta J^N.
}
\]

The linearized geometry sector is
\[
\boxed{
\mu_\eta\partial_t^2\eta
-\partial_w(T_w\partial_w\eta)
-T_\Omega\Delta_{S^2}\eta
+K_\eta\eta
=
S_\eta^{(\psi)}+S_\eta^{(A)}+f_{\rm ext}.
}
\]

This is still a reduced-sector scaffold, not the fully solved coupled PDE.
But it is the point at which the wall, support, and mixed gauge sectors enter one
common linearized hierarchy.

### 4.8 Reduced conservative grouped-lane bundle

**Status:** `Exact Within Closure`

The moving-throat PDE is expected to reduce, lane by lane, to the grouped real
`P2` bundle.

For each grouped lane \(A\in\{20,21,22\}\), define:

- wall/worldtube amplitude \(q_A\),
- stable BdG support modes \(X_{A\alpha}\) with frequencies \(\varpi_{A\alpha}\),
- localized brane-like gauge coordinates \(U_{A,r}\) with frequencies \(\Omega_{U,A,r}\),
- mixed coordinates \(W_{A,r}\) with frequencies \(\Omega_{W,A,r}\),
- internal mixed-sector couplings \(R_{A,r}\).

On the stable separated-pole branch, eliminating the stable BdG modes gives the
exact conservative moments
\[
B_{A,0}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\qquad
B_{A,2}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\qquad
B_{A,4}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

For each Maxwell/mixed port \(r\), define
\[
\Delta_{A,r}=\Omega_{U,A,r}^2\Omega_{W,A,r}^2-R_{A,r}^2,
\qquad
S_{A,r}=\Omega_{U,A,r}^2+\Omega_{W,A,r}^2,
\]
\[
Q_{A,r}=g_{U,A,r}^2\Omega_{W,A,r}^2+2g_{U,A,r}g_{W,A,r}R_{A,r}+g_{W,A,r}^2\Omega_{U,A,r}^2,
\qquad
G_{A,r}=g_{U,A,r}^2+g_{W,A,r}^2.
\]

Then
\[
Z_{A,0}^{(r)}=\frac{Q_{A,r}}{\Delta_{A,r}},
\qquad
Z_{A,2}^{(r)}=\frac{Q_{A,r}S_{A,r}-G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^2},
\]
\[
Z_{A,4}^{(r)}=
\frac{Q_{A,r}(S_{A,r}^2-\Delta_{A,r})-S_{A,r}G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^3}.
\]

Summing over ports gives \(Z_{A,n}\).

The full conservative grouped-lane operator is
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6),
\]
with
\[
\boxed{
D_{A,0}=K_A-B_{A,0}-Z_{A,0},
}
\qquad
\boxed{
D_{A,2}=-\big(M_A+B_{A,2}+Z_{A,2}\big),
}
\]
\[
\boxed{
D_{A,4}=-\big(B_{A,4}+Z_{A,4}\big).
}
\]

So at this level the completed PDE is expected to supply a grouped conservative
bundle of static, quadratic, and quartic low-frequency coefficients lane by lane.

### 4.9 Conservative-lane status at this point

**Fixed within the present reduced hierarchy**

- the linearized coupled bulk/interface scaffold,
- the exact conservative BdG and Maxwell/mixed moment formulas on the selected
  stable reduced branch,
- the grouped conservative operator coefficients \((D_{A,0},D_{A,2},D_{A,4})\).

**Still open**

- the actual values of the branch data on the completed moving-throat solution,
- whether the physical branch lands on the isotropic/passive/outgoing route,
- the outgoing odd normalization data.

---

## 5. Family-1 Core/Mouth Closure Stack

### 5.1 Core/mouth variables and normalized Family-1 ratios

**Status:** `Exact Within Closure`

The core/mouth branch is parameterized by the microscopic throat-core quantities
\[
K_s,\qquad K_q,\qquad \lambda,\qquad g_s,\qquad g_q.
\]

The normalized Family-1 ratios are
\[
\boxed{
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g_c=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
}
\]

The parent compensation condition is
\[
1+\mathfrak r^2=4(\mathfrak g_c-\mathfrak r)^2,
\]
so on the compensated family
\[
\mathfrak g_c
=
\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}.
\]

The lower compensated branch is
\[
\boxed{
\mathfrak g_-(\mathfrak r)=\mathfrak r-\frac12\sqrt{1+\mathfrak r^2}.
}
\]

### 5.2 Positive localized mouth-source theorem

**Status:** `Exact Within Closure`

Inside the first one-lane positive localized-source closure on the first D/N
interval, take a nonnegative normalized axial source profile
\[
\sigma(z)\ge 0,
\qquad
\int_0^L \sigma(z)\,dz=1.
\]

The normalized mouth-bias factor is
\[
\boxed{
\mathfrak g[\sigma]
=
\int_0^L \sigma(z)\cos\!\left(\frac{\pi z}{2L}\right)\,dz.
}
\]

Because
\[
0\le \cos\!\left(\frac{\pi z}{2L}\right)\le 1,
\]
every positive normalized source law satisfies
\[
\boxed{
0\le \mathfrak g[\sigma]\le 1.
}
\]

So within this one-lane positive localized-source setup:

- the upper compensated Family-1 branch is impossible,
- the lower compensated branch is the only admissible compensated candidate.

This theorem does **not** cover arbitrary sign-changing, multimode, or
nonlocalized mouth data.

### 5.3 Explicit mouth boundary-layer law

**Status:** `Effective Closure`

The explicit mouth boundary layer is modeled by the positive source-density free energy
\[
\boxed{
F_{\rm mouth}[\sigma]
=
\int_0^L dz\,
\Big[
\Theta_\sigma\,\sigma\!\left(\ln\frac{\sigma}{\sigma_*}-1\right)
+
V_m(z)\sigma
\Big].
}
\]

Near the mouth, the effective potential is linearized as
\[
\boxed{
V_m(z)\approx V_1 z,
\qquad
V_1=\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-q_\star\left.\partial_zA_0\right|_{\rm m}.
}
\]

That linear potential is a reduced mouth-layer closure on the active interval
\([0,L]\), not a claim that the full throat potential is globally linear.

Using the positive Onsager current
\[
J_\sigma=-M_\sigma\,\sigma\,\partial_z\mu_\sigma^{\rm chem},
\]
the stationary zero-flux branch is exactly
\[
\boxed{
\sigma_\Pi(z)
=
\frac{\Pi e^{-\Pi z/L}}{L(1-e^{-\Pi})},
\qquad
\Pi=\frac{V_1L}{\Theta_\sigma}>0.
}
\]

So the previously ad hoc truncated exponential source family becomes the exact
zero-flux equilibrium branch of the reduced mouth-layer model.

### 5.4 Explicit core-to-mouth gain map

**Status:** `Exact Within Closure`

The explicit throat-core ansatz gives the actual coupled mouth gains
\[
\boxed{
M_s=\frac{Lg_s^2}{K_s\Theta_\sigma},
\qquad
M_q=-\frac{L(K_sg_q-\lambda g_s)^2}{K_s(K_sK_q+\lambda^2)\Theta_\sigma}.
}
\]

In normalized variables, define
\[
\boxed{
\Sigma_0:=\frac{Lg_s^2}{K_s\Theta_\sigma}.
}
\]
Then
\[
\boxed{
M_s=\Sigma_0,
\qquad
M_q=-\Sigma_0\,R_q,
\qquad
R_q=\frac{(\mathfrak g_c-\mathfrak r)^2}{1+\mathfrak r^2}.
}
\]

So the coupled Family-1 mouth law is
\[
\boxed{
\Pi=\Sigma_0\Big[1-R_q\,\mathcal S_q(\Pi)\Big].
}
\]

On the exact compensated branch,
\[
\boxed{
R_q=\frac14.
}
\]

So the outlet-consistent one-parameter mouth closure is derived rather than assumed.

### 5.5 Self-matched mouth susceptibility closure

**Status:** `Effective Closure`

The first self-matched mouth susceptibility closure identifies the source channel
with the same active shell layer as the shell/compliance mode:
\[
\boxed{
\Theta_\sigma=H_wJ_s.
}
\]

This is a same-layer closure.
It removes an otherwise free susceptibility scale rather than introducing a new
fit parameter. If the actual mouth source lives on a materially different layer,
the resulting prefactor need not stay unchanged.

Under this closure,
\[
\boxed{
\Sigma_0=\frac{20}{9}\,\widehat T_m^2,
\qquad
\widehat T_m:=\frac{\rho_w\ell\sqrt L\,\mathcal T_m}{\hbar c_{s,w}}.
}
\]

So the shell gain is explicit in terms of a normalized mouth-traction amplitude.

### 5.6 Actual Family-1 gains on the canonical point

**Status:** `Numerically Located`

For the concrete Family-1 branch,
\[
\mathfrak r_{F1}\approx 1.77799353547498.
\]

At the canonical compensation point \(\Pi_*\),
\[
\Pi_*\approx 1.50882951349316,
\qquad
\mathcal S_q(\Pi_*)\approx 0.658075937605429.
\]

The natural equal-normalized branch \(\mathfrak g_c=1\) gives
\[
M_s^{\rm nat,*}\approx 1.66854252965624,
\qquad
M_q^{\rm nat,*}\approx -0.242696939724365.
\]

The exact compensated branch gives
\[
M_s^{\rm comp,*}\approx 1.80594111095636,
\qquad
M_q^{\rm comp,*}\approx -0.451485277739090.
\]

So:

- the shell gain changes by about `8.23%`,
- the mixed-gain magnitude increases by about a factor `1.86`,
- while under the self-matched susceptibility closure the normalized mouth
  traction differs by only about `4.04%`.

This means the natural branch is on the correct sign side and not far from the
canonical compensated branch, but it does not reproduce it exactly.

### 5.7 Self-consistent explicit mouth branch

**Status:** `Exact Within Closure`

On the explicit positive mouth family, identify
\[
\mathfrak g_c=\mathfrak g_\Pi,
\]
where
\[
\boxed{
\mathfrak g_\Pi
=
\frac{2\Pi(2\Pi e^\Pi+\pi)}{(4\Pi^2+\pi^2)(e^\Pi-1)}.
}
\]

Then the mixed-to-shell ratio becomes
\[
\boxed{
R_q(\Pi)=\frac{(\mathfrak g_\Pi-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
}
\]
and the Family-1 mouth law closes to
\[
\boxed{
\Pi=\Sigma_0\Big[1-R_q(\Pi)\mathcal S_q(\Pi)\Big].
}
\]

Equivalently,
\[
\boxed{
\Sigma_0(\Pi)=\frac{\Pi}{1-R_q(\Pi)\mathcal S_q(\Pi)}.
}
\]

Under the self-matched closure,
\[
\boxed{
\widehat T_m(\Pi)=
\sqrt{\frac{9\Pi}{20\left[1-R_q(\Pi)\mathcal S_q(\Pi)\right]}}.
}
\]

So the explicit mouth branch is governed by a single scalar bias parameter
\(\Pi\), not a free gain pair.

### 5.8 Singular equal-normalized limit and unique regular canonical branch

**Status:** `Exact Within Closure`

Inside the explicit positive exponential mouth family:

1. the equal-normalized branch \(\mathfrak g_c=1\) is **not** a finite-bias branch;
2. it is reached only in the singular point-source limit \(\Pi\to\infty\);
3. and in the same limit the normalized mouth traction diverges.

So the equal-normalized branch is not a regular finite branch of the explicit
mouth-layer dynamics.

The lower compensated Family-1 value is
\[
\mathfrak g_-^{F1}\approx 0.758035078944663.
\]

Since
\[
\frac{2}{\pi}<\mathfrak g_-^{F1}<1,
\]
there exists a unique positive finite \(\Pi_*\) such that
\[
\mathfrak g_{\Pi_*}=\mathfrak g_-^{F1}.
\]

Numerically,
\[
\boxed{
\Pi_*\approx 1.50882951349316,
\qquad
\Sigma_0(\Pi_*)\approx 1.80594111095636,
\qquad
\widehat T_m(\Pi_*)\approx 0.901484054174205.
}
\]

So inside the explicit Family-1 positive exponential mouth-layer closure:

- the upper compensated branch is impossible,
- the equal-normalized branch is singular,
- the lower compensated branch is the unique regular finite-bias / finite-traction branch.

### 5.9 Family-1 mouth-stack status at this point

**Fixed within the explicit mouth closure stack**

- one-lane positive-source branch-selection theorem,
- explicit exponential mouth boundary layer,
- explicit core-to-mouth gain map,
- normalized gain ratio \(R_q\),
- self-matched susceptibility closure,
- self-consistent explicit mouth branch,
- unique regular canonical branch inside the explicit positive family.

**Still open**

- whether the actual moving-throat mouth layer is realized closely enough by this
  one-lane positive exponential closure,
- whether the same-layer susceptibility closure is the physically realized one,
- and how the full co-evolving mouth/core branch shifts the canonical point once
  backreaction is turned back on.

---

## 6. Co-Evolving Mouth/Core Branch

### 6.1 Exact canonical full-profile residual

**Status:** `Exact Within Closure`

On the explicit canonical Family-1 mouth branch, the exponential source
\[
\Sigma_*(x)=\frac{\Pi_*e^{-\Pi_*x}}{1-e^{-\Pi_*}},
\qquad x\in[0,1],
\]
generates exact shell and mixed profiles
\[
T_s(x;\Pi_*)
=
\frac{1-e^{-\Pi_* x}}{\Pi_*(1-e^{-\Pi_*})}
-
\frac{x e^{-\Pi_*}}{1-e^{-\Pi_*}},
\]
\[
T_q(x;\Pi_*)
=
A_q\sinh\!\left(\frac{\pi x}{2}\right)
-
C_q\cosh\!\left(\frac{\pi x}{2}\right)
+
C_q e^{-\Pi_* x},
\]
with
\[
C_q=
\frac{\Pi_*}{(1-e^{-\Pi_*})(\pi^2/4-\Pi_*^2)},
\]
\[
A_q=
\frac{C_q\left(\frac{\pi}{2}\sinh(\pi/2)+\Pi_*e^{-\Pi_*}\right)}
{\frac{\pi}{2}\cosh(\pi/2)}.
\]

Using the canonical compensated outlet data
\[
M_s^*=4\Sigma_m^*,
\qquad
M_q^*=-\Sigma_m^*,
\qquad
\Sigma_m^*\approx 0.451485277739090,
\]
the full mouth potential is
\[
\Phi_*(x)=4\Sigma_m^*\,T_s(x;\Pi_*)-\Sigma_m^*\,T_q(x;\Pi_*).
\]

Relative to the tangent exponential potential \(\Pi_*x\), define
\[
\boxed{
R_*(x):=\Phi_*(x)-\Pi_*x.
}
\]

Then the canonical branch satisfies
\[
\boxed{
R_*(0)=0,
\qquad
R_*'(0)=0,
\qquad
R_*''(0)
=
-\,3\Sigma_m^*\frac{\Pi_*}{1-e^{-\Pi_*}}<0.
}
\]

So the explicit exponential source is tangent-matched at the mouth, but the full
coupled potential is locally sublinear there. Within the mouth-only closure this
means the true self-consistent source broadens relative to \(\Sigma_*\).

### 6.2 First-order mouth-only correction and one-step nonlinear check

**Status:** `Reduced / Controlled Reduction`

The exact full positive source on the mouth-only branch is
\[
\Sigma_{\rm full}(x)\propto e^{-\Phi_*(x)}=e^{-\Pi_*x-R_*(x)}.
\]

Expanding to first order around \(\Sigma_*\) gives
\[
\boxed{
\delta\Sigma_{\rm act}(x)
=
-\Sigma_*(x)\,\widetilde R_*(x),
\qquad
\widetilde R_*(x)=R_*(x)-\langle R_*\rangle_*.
}
\]

Only two projected moment shifts matter:
\[
\delta\mathfrak g_{\rm act}
=
-\operatorname{Cov}_*(c,R_*),
\qquad
\delta\mathcal S_{\rm act}
=
-\operatorname{Cov}_*(K_q,R_*),
\]
where
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

On the canonical explicit branch,
\[
\boxed{
\delta\mathfrak g_{\rm act}\approx -0.0648069687666328,
\qquad
\delta\mathcal S_{\rm act}\approx -0.0388718368650403.
}
\]

Using the previously frozen canonical rigidity data,
\[
\mathfrak g_*' \approx 0.0714453558083195,
\qquad
A_T\approx -4.27263956256927,
\qquad
B_T\approx 0.134875005736706,
\]
the first-order retuning is
\[
\boxed{
\delta\Pi_{\rm act}\approx 0.907084414842908,
\qquad
\delta\widehat T_{m,\rm act}\approx 0.271653979462338.
}
\]

So the mouth-only corrected point is
\[
\boxed{
\Pi_{\rm corr}\approx 2.41591392833607,
\qquad
\widehat T_{m,\rm corr}\approx 1.17313803363654.
}
\]

This is a controlled first-order mouth-only correction, not an exact finite
nonlinear law. A one-step nonlinear Picard probe gives
\[
\Pi_1\approx 2.53914847609768,
\qquad
\widehat T_{m,1}\approx 1.21036942084359,
\]
which moves in the same direction and on the same scale, but does not by itself
prove full convergence of the mouth-only iteration.

### 6.3 Exact co-evolving Family-1 map

**Status:** `Exact Within Closure`

The next reduction lets the positive mouth source and the core loading ratio
co-evolve. For any normalized positive profile
\[
\Sigma(x)\ge 0,
\qquad
\int_0^1 \Sigma(x)\,dx=1,
\]
define the moments
\[
\mathfrak g[\Sigma]
=
\int_0^1 \Sigma(x)\cos\!\left(\frac{\pi x}{2}\right)\,dx,
\]
\[
\mathcal S[\Sigma]
=
\int_0^1 \Sigma(x)
\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}\,dx.
\]

On the explicit Family-1 core branch, the shell/mixed ratio is no longer free:
\[
\boxed{
\mathcal R[\Sigma]
=
\frac{(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak r_{F1}\approx 1.77799353547498.
}
\]

The shell and mixed kernels are
\[
\mathcal T_s[\Sigma](x)=\int_0^1 \min(x,y)\,\Sigma(y)\,dy,
\]
\[
\mathcal T_q[\Sigma](x)
=
\int_0^1
\frac{\sinh\!\left(\frac{\pi}{2}\min(x,y)\right)
\cosh\!\left(\frac{\pi}{2}(1-\max(x,y))\right)}{(\pi/2)\cosh(\pi/2)}
\Sigma(y)\,dy.
\]

So the exact co-evolving mouth potential is
\[
\boxed{
\Phi_{\Sigma_0}[\Sigma](x)
=
\Sigma_0\Big[\mathcal T_s[\Sigma](x)-\mathcal R[\Sigma]\mathcal T_q[\Sigma](x)\Big],
}
\]
and the reduced fixed-point equation is
\[
\boxed{
\Sigma(x)=
\frac{e^{-\Phi_{\Sigma_0}[\Sigma](x)}}
{\int_0^1 e^{-\Phi_{\Sigma_0}[\Sigma](y)}\,dy}.
}
\]

The lower compensated branch is still characterized by
\[
\mathfrak g[\Sigma]=\mathfrak g_*,
\qquad
\mathfrak g_*\approx 0.758035078944663,
\]
which is exactly equivalent to
\[
\boxed{
\mathcal R[\Sigma]=\frac14.
}
\]

Near the canonical branch,
\[
\boxed{
\delta\mathcal R
=
-\frac{\delta\mathfrak g}{\sqrt{1+\mathfrak r_{F1}^2}}
+O(\delta\mathfrak g^2),
}
\]
so source broadening \((\delta\mathfrak g<0)\) drives the mixed loading ratio
above \(1/4\).

For any self-consistent profile on this reduced branch,
\[
\boxed{
\Pi[\Sigma]
=
\Phi'_{\Sigma_0}[\Sigma](0)
=
\Sigma_0\Bigl[1-\mathcal R[\Sigma]\mathcal S[\Sigma]\Bigr].
}
\]

Under the self-matched susceptibility closure,
\[
\boxed{
\Sigma_0=\frac{20}{9}\widehat T_m^2.
}
\]

### 6.4 Frozen-traction fixed point

**Status:** `Numerically Located`

Holding the old canonical traction fixed at
\[
\Sigma_0^*\approx 1.80594111095636,
\qquad
\widehat T_{m,*}\approx 0.901484054174204,
\]
the exact co-evolving map on the analyzed positive branch window converges to a
positive fixed point with
\[
\boxed{
\mathfrak g_{\rm fp}\approx 0.693352419668063,
\qquad
\mathcal S_{\rm fp}\approx 0.6216013167514007,
\qquad
\mathcal R_{\rm fp}\approx 0.2827139049082381,
\qquad
\Pi_{\rm fp}\approx 1.4885734438300713.
}
\]

So the branch survives and the mouth bias stays close to \(\Pi_*\), but the fixed
point no longer lands exactly on the compensated value \(\mathcal R=1/4\).

This is strong anti-tuning evidence inside the reduced closure: once backreaction
is allowed, the old canonical point does not protect itself.

### 6.5 Renormalized canonical branch under full co-evolution

**Status:** `Numerically Located`

The exact restoration condition is
\[
\mathfrak g_{\rm fp}(\Sigma_0)=\mathfrak g_*.
\]
On the analyzed positive branch window, this condition has a unique numerically
located root
\[
\boxed{
\Sigma_0^{\rm can}\approx 4.651033550168876,
\qquad
\widehat T_{m,\rm can}\approx 1.446708366456762.
}
\]

At that renormalized traction, the co-evolving reduced fixed point satisfies
\[
\boxed{
\mathfrak g_{\rm can}=\mathfrak g_*,
\qquad
\mathcal R_{\rm can}=\frac14,
\qquad
\mathcal S_{\rm can}\approx 0.6703621156734617,
\qquad
\Pi_{\rm can}\approx 3.8715643774790087.
}
\]

Relative to the original canonical point, exact reduced compensation costs roughly
\[
\boxed{
\frac{\Sigma_0^{\rm can}}{\Sigma_0^*}-1\approx 1.5754070949223031,
}
\]
\[
\boxed{
\frac{\Pi_{\rm can}}{\Pi_*}-1\approx 1.5659389234213572,
}
\]
\[
\boxed{
\frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1\approx 0.6048074946616844.
}
\]

So the compensated Family-1 branch survives the analyzed co-evolving closure, but
only as a renormalized finite-bias, finite-traction branch rather than at the old
canonical point.

### 6.6 Co-evolving mouth/core status at this point

**Fixed within the analyzed reduced closure**

- the exact co-evolving Family-1 fixed-point map,
- the first-order sign of the compensation defect transport,
- failure of exact compensation at the old canonical traction,
- the numerically located renormalized compensated branch on the analyzed
  positive window.

**Still open**

- whether the full moving-throat PDE dynamically realizes this reduced
  co-evolving branch closely enough for referee-level purposes,
- whether other microscopic branches exist outside the analyzed positive reduced
  window,
- and how the co-evolving reduced branch feeds the final outgoing grouped-`P2`
  normalization on the completed PDE solution.

---

## 7. Grouped Real `P2` Conservative and Retarded Finish Line

### 7.1 Normalized grouped conservative response and weighted trace/anomaly map

**Status:** `Exact Within Closure`

Define the normalized grouped conservative response
\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6),
\qquad A\in\{20,21,22\}.
\]
Then
\[
\boxed{
 u_2^{(A)}=-\frac{D_{A,2}}{D_{A,0}},
}
\qquad
\boxed{
 u_4^{(A)}=\frac{D_{A,2}^2-D_{A,0}D_{A,4}}{D_{A,0}^2}.
}
\]

For any grouped triple `\(x_A\)`, define the weighted trace/anomaly variables
\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]
\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2}.
\]
Applied to `\(u_2^{(A)}\)` and `\(u_4^{(A)}\)`, the exact grouped conservative
isotropy gate is
\[
\boxed{a_2=b_2=a_4=b_4=0.}
\]

### 7.2 Exact isotropic one-pole surface and one-parameter conservative carrier

**Status:** `Exact Within Closure`

On the isotropic grouped branch,
\[
D_{20,n}=D_{21,n}=D_{22,n}=:D_n,
\qquad
N_{20,n}=N_{21,n}=N_{22,n}=:N_n.
\]
Then the grouped conservative moments collapse to one common carrier
\[
D_Q^{\rm(cons)}(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]
and the common response moments are
\[
\bar u_2=-\frac{D_2}{D_0},
\qquad
\bar u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

Define the exact conservative one-pole defect
\[
\boxed{\Delta_{\rm pole}:=\bar u_4-4\bar u_2^{\,2}.}
\]
A short algebraic reduction gives
\[
\boxed{
\Delta_{\rm pole}=-\frac{3D_2^2+D_0D_4}{D_0^2}.
}
\]
So the isotropic one-pole surface is exactly
\[
\boxed{D_0D_4+3D_2^2=0.}
\]
If one rewrites the conservative moments in the wall/BdG/Maxwell bundle,
\[
D_0=K-B_0-Z_0,
\qquad
D_2=-(M+B_2+Z_2),
\qquad
D_4=-(B_4+Z_4),
\]
this becomes
\[
\boxed{D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.}
\]

Define
\[
\boxed{\Omega_Q^2:=-\frac{D_0}{4D_2}.}
\]
Then on the isotropic one-pole surface the conservative grouped-real `P2`
carrier is forced to the one-parameter form
\[
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
}
\]
So the exact conservative front-end theorem target is
\[
\boxed{a_2=b_2=a_4=b_4=0,\qquad \Delta_{\rm pole}=0.}
\]

### 7.3 Exact scalar/geometry firewall

**Status:** `Exact Within Closure`

At the isotropic quadratic wall/support level, the operator is block-diagonal in
angular momentum:
\[
\mathcal D^{(0)}(\omega)=\mathcal D_{l=0}(\omega)\oplus\mathcal D_{l=2}(\omega)\oplus\cdots.
\]
So the `l=0` scalar/geometry lane does not linearly enter the grouped `l=2`
carrier on the exact isotropic background.

If a small anisotropy parameter `\(\varepsilon_{\rm mix}\)` turns on the first
`l=0\leftrightarrow l=2` mixing,
\[
\mathcal D(\omega,\varepsilon_{\rm mix})
=
\begin{pmatrix}
\mathcal D_0(\omega) & \varepsilon_{\rm mix} C(\omega)^T\\[4pt]
\varepsilon_{\rm mix} C(\omega) & \mathcal D_2(\omega)I_3
\end{pmatrix},
\]
then the exact Schur complement gives
\[
\boxed{
\mathcal D_{2,\rm eff}(\omega,\varepsilon_{\rm mix})
=
\mathcal D_2(\omega)I_3
-
\varepsilon_{\rm mix}^2 C(\omega)\mathcal D_0(\omega)^{-1}C(\omega)^T.
}
\]
Therefore
\[
\boxed{
\partial_{\varepsilon_{\rm mix}}\mathcal D_{2,\rm eff}(\omega,\varepsilon_{\rm mix})\big|_{\varepsilon_{\rm mix}=0}=0,
}
\]
and scalar/geometry contamination of the grouped `l=2` carrier begins only at
`O(\varepsilon_{\rm mix}^2)`.

### 7.4 Exact outgoing `l=2` fingerprint and isotropic deformation algebra

**Status:** `Exact Within Closure`

Define the dimensionless outgoing argument
\[
\boxed{z:=\frac{a\omega}{c_s}.}
\]
The exact outgoing spherical `l=2` DtN operator is
\[
\Lambda_2^{\rm out}(z)=z\,\frac{h_2^{(1)\prime}(z)}{h_2^{(1)}(z)}
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6).
\]
Normalize by the static slot:
\[
\boxed{\widehat Y_2^{\rm out}(z):=-\frac{3}{\Lambda_2^{\rm out}(z)}.}
\]
Then
\[
\boxed{
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]

Match this to the isotropic retarded grouped-`P2` module
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+
\frac14\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5}
+O(\omega^6),
\]
with
\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\sigma_Q^{\rm can}=\frac{9}{8\Omega_Q^5}=\frac{4a^5}{27c_s^5}.
\]
The exact compact outgoing branch therefore fixes
\[
\boxed{\chi_Q=1.}
\]

A convenient isotropic DtN deformation family preserving the same even branch is
\[
\boxed{
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6).
}
\]
Canonical-even matching fixes
\[
\boxed{
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
}
\]
The exact outgoing-normalization scalar is then
\[
\boxed{
\chi_Q=
\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}.
}
\]
Linearizing about the canonical branch,
\[
\beta=1+\varepsilon_\beta,
\qquad
\Sigma_0=\delta\Sigma_0,
\qquad
\Sigma_5=\delta\Sigma_5,
\]
one gets
\[
\boxed{
\chi_Q-1
=
5\varepsilon_\beta
+\frac{\delta\Sigma_0}{3S}
+\frac{9\,\delta\Sigma_5}{S}
+O(2).
}
\]
So the entire isotropic retarded realization problem is already localized in the
three deformation coordinates `(\beta,\Sigma_0,\Sigma_5)`.

### 7.5 Source-map reduction, higher-odd irrelevance, and exact Packet-A finish line

**Status:** `Exact Within Closure`

The exact observable odd closure factorizes as
\[
\boxed{m_{\hat 0}^{\,2}\chi_Q N_Q=1.}
\]
On the natural point-particle source-map branch,
\[
\boxed{m_{\hat 0}\to 1,\qquad N_Q=\chi_Q^{-1}.}
\]
Therefore the Packet-A normalization slot is
\[
\boxed{
\Delta_{\rm norm}=P_0^{\rm target}(\chi_Q^{-1}-1),
\qquad
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
}
\]

Any extra isotropic odd retarded correction whose first nonzero term is
`O(\omega^7)` leaves `\chi_Q`, `N_Q`, and `\Delta_{\rm norm}` unchanged at theorem
order. So there is no remaining higher-odd loophole in the reduced point-particle
`2.5`PN finish line.

Once the isotropic one-pole conservative surface and the isotropic outgoing lane
are imposed,
\[
a_2=b_2=a_4=b_4=a_{P_0}=b_{P_0}=\Delta_{\rm pole}=0,
\]
and the full Packet-A residual collapses to one scalar only:
\[
\boxed{\Delta_{\rm branch}=0\iff\chi_Q=1.}
\]
Equivalently, on the natural point-particle source-map branch,
\[
\boxed{\Delta_{\rm branch}=0\iff N_Q=1\iff \Delta_Q:=\chi_Q-1=0.}
\]
The exact isotropic DtN realization gate may be written as
\[
\boxed{\chi_Q=1\iff 3S(\beta^5-1)+\Sigma_0+27\Sigma_5=0.}
\]
So Packet A is now completely localized: the whole branch-side retarded finish
line is one scalar equation.


## 7A. V2 Full-Bundle Target Surface, 5PN Finish Line, and Executable Branch Pipeline

This section is the V2 audit layer that should be used whenever a future PDE
branch or solver export is being evaluated. It does not replace the older
Packet-A/Packet-B and coherent-kernel machinery; it gives the frozen residual
surface and extraction protocol that the actual branch must satisfy.

### 7A.1 Isotropic full-bundle target card

On the isotropic grouped real `P2` branch, define
\[
D(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]
\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6),
\]
with
\[
\boxed{D_0=K-B_0-Z_0,}
\qquad
\boxed{D_2=-(M+B_2+Z_2),}
\qquad
\boxed{D_4=-(B_4+Z_4).}
\]
The normalized conservative response and outgoing prefactor are
\[
\boxed{u_2=-D_2/D_0,}
\qquad
\boxed{u_4=(D_2^2-D_0D_4)/D_0^2,}
\]
\[
\boxed{P_0=N_0/D_0,}
\]
\[
\boxed{P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},}
\]
\[
\boxed{
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
}
\]

The one-pole conservative gate is
\[
\boxed{
R_{\rm pole}:=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2=0.
}
\]
The constant-prefactor outgoing gate is
\[
\boxed{P_2=0,
\qquad P_4=0,}
\]
equivalently
\[
\boxed{N_2=2D_2N_0/D_0,}
\]
\[
\boxed{N_4=N_0(D_2^2+2D_0D_4)/D_0^2.}
\]
The universal quadrupole normalization is
\[
\boxed{
R_{\rm norm}:=\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-\frac{54Gc_s^5}{5a^5c^5}=0.
}
\]
Equivalently,
\[
\boxed{
\gamma_{\rm quad}^{\rm eff}
=\widehat m_0^{\,2}\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}
=\frac{2G}{5c^5}.
}
\]
If the 4PN tail scalar has not already been derived independently, include
\[
\boxed{R_{\rm tail}:=\Theta_{\rm tail}(c/c_s)^3-1=0.}
\]

### 7A.2 Local target-sheet formulas

On the target sheet,
\[
\boxed{D_0=\frac{3(B_2+M+Z_2)^2}{B_4+Z_4},}
\]
\[
\boxed{K=B_0+Z_0+\frac{3(B_2+M+Z_2)^2}{B_4+Z_4},}
\]
\[
\boxed{P_0^{\rm target}=\frac{54Gc_s^5}{5\mathcal S_{\rm port}a^5c^5\widehat m_0^{\,2}},}
\]
\[
\boxed{
N_0=\frac{162Gc_s^5(B_2+M+Z_2)^2}
{5\mathcal S_{\rm port}a^5c^5\widehat m_0^{\,2}(B_4+Z_4)},
}
\]
\[
\boxed{
N_2=-\frac{108Gc_s^5(B_2+M+Z_2)}
{5\mathcal S_{\rm port}a^5c^5\widehat m_0^{\,2}},
}
\]
\[
\boxed{
N_4=-\frac{18Gc_s^5(B_4+Z_4)}
{\mathcal S_{\rm port}a^5c^5\widehat m_0^{\,2}}.
}
\]
The target packet is locally codimension five in the output slots
\((K,N_0,N_2,N_4,\Theta_{\rm tail})\), so the branch-freeze rule is essential.

### 7A.3 Grouped projectors and angular source-map memory card

The grouped metric is
\[
G_{\rm grp}=\mathrm{diag}(1,2,2).
\]
For a grouped triple \(x=(x_{20},x_{21},x_{22})\),
\[
\boxed{x=\bar x(1,1,1)+a_x(4,-1,-1)+b_x(0,1,-1),}
\]
where
\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\qquad
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}2.
\]
The grouped anisotropy norm is
\[
\boxed{A_x^2=4a_x^2+\frac45b_x^2.}
\]
A weak axisymmetric quadrupole perturbation has fixed signature
\[
(20,21,22)\sim\left(1,\frac12,-1\right),
\qquad
\boxed{b=3a.}
\]
The canonical real-STF angular source map is identity on the isotropic branch:
\[
\boxed{\widehat m_{\rm ang}=1.}
\]
Thus the remaining normalization gap is radial/axial, source-amplitude,
denominator \(D_0\), port factor \(\mathcal S_{\rm port}\), or mixed-transfer
\(N_0\), not an angular STF-basis mismatch.

### 7A.4 Stability, passivity, and non-dark-port gates

The branch must pass at least these reduced gates before target residuals are
scientifically meaningful:

- open finite-radius exit: \(R(L)>0\), no hard cap;
- wall positivity: \(\mu_\eta>0\), \(T_w>0\), \(K_\eta+6T_\Omega\ge0\) in the grouped lane;
- BdG support positivity: retained BdG modes have positive stable frequencies;
- Schur static stability: \(D_0=K-B_0-Z_0>0\);
- Maxwell/mixed block positivity:
  \[
  \boxed{\Delta=\Omega_U^2\Omega_W^2-R^2>0;}
  \]
- outgoing transfer non-darkness: \(N_0\ne0\), preferably \(N_0>0\);
- no ghost/Krein contamination in the retained positive-energy reduction;
- target-blind pre-freeze: no post-residual retuning of \(K,M,B_n,Z_n,N_n\).

### 7A.5 Reduced 5PN / actual-branch finish-line card

The reduced 5PN package is not a standard full conservative 5PN two-body
assembly. It is the carry-forward finish-line compression of the 5PN / 2.5PN /
4PN moving-throat branch problem.

On the actual coherent branch, the reduced finish line is
\[
\boxed{
 d\ln R_{\rm tr}=0,
 \qquad
 d\ln R_{\rm target}=0,
 \qquad
 d\ln\epsilon_\eta=0,
 \qquad
 N_Q=1.
}
\]
On the natural source-map branch,
\[
\boxed{N_Q=\chi_Q^{-1},
\qquad N_Q=1\iff \chi_Q=1.}
\]
The weak-axisymmetric grouped normalization defect is
\[
\boxed{\Xi_1=\frac{P_1}{P_0}.}
\]
The exact coherent normal form uses branch-adapted monomial drifts
\[
\Sigma_{\rm tr}=\delta\ln\mathfrak C_{{\rm tr},*},
\qquad
\Sigma_{\rm nt}=\delta\ln\mathfrak C_{{\rm nt},*},
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta,
\]
with triangular observable structure
\[
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
\]
The exact similarity-orbit theorem says that zero weak-axisymmetric defect is
equivalent to tangent motion along the monomial-preserving similarity orbit:
\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\mathbf x\in T_{\rm id}\mathcal G_*.
\]
So the reduced 5PN continuation point is not another coefficient search. It is
whether the actual moving-throat branch lands on the isotropic state, the
weak-axisymmetric orbit-lock tangent, and the outgoing normalization condition
on one and the same branch.

### 7A.6 Solver/extraction pipeline card

The executable V2 chain is:
\[
\boxed{
\text{solver export}
\to
\text{V2-22B validation}
\to
\text{V2-22A profile adapter}
\to
\text{V2-21 grouped observable extraction}
\to
\text{residual/tolerance report}.
}
\]

Required frozen solver metadata:

```text
pre_target_freeze = true
target_blind = true
no_post_residual_refit = true
boundary_class = open_impedance
R_exit > 0
no_hard_cap = true
at_least_one_mixed_Aw_port = true
```

Required extracted data:

```text
K, M
B0, B2, B4
Z0, Z2, Z4
N0, N2, N4
mhat0
S_port
Theta_tail    # if tail gate is active
```

Failure classes:

1. handoff/validation failure;
2. open-throat/stability failure;
3. target-realization failure.

A stable branch may fail the target packet. That is useful negative data, not a
pipeline error.

### 7A.7 V2-23 first target-blind open-throat branch result

The first reduced open-throat branch is a frozen one-dimensional branch:

```text
branch_id: v2_23_minimal_open_throat_frozen_demo
pre_target_freeze: true
target_blind: true
no_post_residual_refit: true
boundary_class: open_impedance
boundary_protocol: open_impedance_AC_reflecting_DC_leaking
```

It solves an open finite-radius throat with
\[
R(0)=a,
\qquad
R(L)>0,
\]
and reports
\[
R_{\rm mouth}=1,
\qquad
R_{\rm exit}=0.452938042901,
\qquad
R_{\rm min}=0.452938042901.
\]
The branch passes open/stability/outgoing-transfer mechanics, but fails the
full target packet:

```text
target_packet_pass: false
P0 / P0_target = 0.00183195525067
one_pole_ratio = 0.0147641804366
```

Interpretation: the branch is a useful reduced negative control and solver
fixture. It is not a PDE realization of the universal quadrupole target.


## 8. Coherent Local-Kernel Reduction, Realization Compiler, and Finite Local Search

### 8.1 Coherent-kernel effective ratios and tracking branch

**Status:** `Exact Within Closure`

Inside the first coherent local D/N support closure, the mixed and support lanes
couple through the same local wall/\(U\) density, so the coherent branch is
controlled by the effective stiffnesses
\[
K_{U1}=K_U(1+\delta_U),
\qquad
K_\eta^{(\mathrm{eff})}=K_\eta+6T_\Omega,
\]
\[
K_W^{(\mathrm{eff})}=K_W+\frac{\pi^2T_W}{4L^2},
\qquad
K_\phi^{(\mathrm{eff})}=K_\phi+\frac{\pi^2T_\phi}{4L^2},
\]
and the dimensionless ratios
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\]
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\]
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\zeta=\frac{\lambda_\phi^2K_W^{(\mathrm{eff})}}{\lambda_W^2K_\phi^{(\mathrm{eff})}},
\]
\[
\Lambda=\frac{27\pi^2Gc_s^5K_W^{(\mathrm{eff})}}{20a^5c^5\mu_W}.
\]

The coherent local-kernel hypothesis forces
\[
\rho_0=\sigma_0=\chi_0,
\qquad
\epsilon_\phi=\zeta\epsilon_W,
\qquad
Z_\phi=\zeta Z_W,
\]
so the reduced Stage-27 branch lands exactly on the tracking surface
\[
\boxed{
R_{\rm tr}=R_U=R_\phi
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
}
\]

On the constructive coherent branch \(\chi_0>0\), \(\delta_U>0\),
\[
\boxed{
\frac{1}{1+\delta_U}<R_{\rm tr}<1.
}
\]

### 8.2 Coherent placement map and support-compensation theorem

**Status:** `Exact Within Closure`

Define the split mixed blocking ratio
\[
\epsilon
=
\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
\]
Then the coherent branch obeys the exact reduced placement map
\[
\boxed{
M_{\rm mix}
=
\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
}
\]
\[
\boxed{
M_{\rm supp}
=
\frac{8\zeta Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\zeta\epsilon)},
}
\]
\[
\boxed{
M_{\rm tr}=M_{\rm mix}+M_{\rm supp}=M_{\rm mix}S(\zeta;\epsilon),
}
\]
with support-enhancement factor
\[
\boxed{
S(\zeta;\epsilon)
=
1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}.
}
\]

The outgoing demand ratio is
\[
\boxed{
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
}
\]
so \(R_{\rm target}\) is independent of \(\zeta\). The support lane changes only
the available baseline, not the target.

The coherent tracking branch then obeys
\[
M_{\rm tr}=G_{\rm tr}(\xi,\delta;R_{\rm tr}),
\qquad
R_{\rm target}=F_{\rm tr}(\xi,\delta;R_{\rm tr}),
\]
where
\[
G_{\rm tr}(\xi,\delta;R)
=
\frac{9\xi(\xi+\delta)}{9\delta+(9+2R^2)\xi},
\]
\[
F_{\rm tr}(\xi,\delta;R)
=
\frac{[9\delta+(9+2R^2)\xi]^2[9\delta+(9+2R)\xi]^2}
{81(1-\xi)\,[9\delta^2+18\delta\xi+(9+2R^2)\xi^2]^2}.
\]

At fixed \((\xi,\delta)\), the coherent tracking branch is ordered by \(R\):
\[
\frac{dG_{\rm tr}}{dR}<0,
\qquad
\frac{dF_{\rm tr}}{dR}>0.
\]
So because \(R_{\rm tr}<1\) on the constructive split-\(U\) branch, the coherent
local kernel requires more total loading and gives less normalized response than
the old flat branch.

That deficit is still not a reduced-level no-go. On the stable-side domain
\[
0<\epsilon<1,
\qquad
0\le \zeta<1/\epsilon,
\qquad
0<\xi<1,
\]
the support-enhancement factor is strictly increasing and invertible:
\[
\boxed{
\frac{dS}{d\zeta}=\frac{1-\epsilon}{(1-\zeta\epsilon)^2}>0.
}
\]

If the mixed-only baseline lies below the tracking critical load
\[
M_{\rm crit}(\delta,R)
=
G_{\rm tr}(1,\delta;R)
=
\frac{9(1+\delta)}{9\delta+9+2R^2},
\]
and the mixed-only branch has not already reached the target, then the unique
required support ratio is
\[
\boxed{
\zeta_{\rm req}
=
\frac{S_{\rm req}-1}{1+\epsilon(S_{\rm req}-2)},
\qquad
S_{\rm req}:=\frac{M_{\rm req}}{M_{\rm mix}},
}
\]
with \(\zeta_{\rm req}<\zeta_{\rm crit}<1/\epsilon\). So there is no
reduced-level support no-go on the coherent tracking branch.

This theorem remains a reduced existence statement. The open PDE question is
whether the realized physical branch actually supplies \(\zeta\ge\zeta_{\rm req}\).

### 8.3 Microscopic slippage ledger and triangular defect normal form

**Status:** `Exact Within Closure`

The positive microscopic coherent state is
\[
x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U),
\]
with grouped weak-axisymmetric log-drift vector
\[
\delta\mathbf x
=
\begin{pmatrix}
\lambda_1\\
c_1\\
\gamma_1\\
\kappa_U\\
\kappa_\eta\\
\kappa_W\\
\mu_1\\
\tau_1
\end{pmatrix}
=
\begin{pmatrix}
\delta\ln\lambda_W\\
\delta\ln c_{\eta U}\\
\delta\ln\gamma\\
\delta\ln K_U\\
\delta\ln K_\eta^{(\mathrm{eff})}\\
\delta\ln K_W^{(\mathrm{eff})}\\
\delta\ln\mu_W\\
\delta\ln T_U
\end{pmatrix}_{\!\rm grp}.
\]

The five direct microscopic slippages are
\[
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
\]
\[
\Sigma_\chi=\gamma_1+c_1-\kappa_U=\delta\ln\chi_0,
\qquad
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta=\delta\ln\epsilon_\eta,
\]
\[
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W=\delta\ln\epsilon_W,
\qquad
\Sigma_\delta=\tau_1-\kappa_U=\delta\ln\delta_U.
\]

The branch-adapted coordinates compress these to
\[
\boxed{
\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
}
\]
\[
\boxed{
\Sigma_{\rm nt}
=
\Sigma_Z
+
\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
}
\]

Then the three reduced observable drifts take the exact triangular form
\[
\boxed{
\Theta_1
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}\,\Sigma_{\rm tr},
}
\]
\[
\boxed{
\Xi_1
=
\frac{2\chi_0}{(1+\chi_0)(1+\delta_U)}\,\Sigma_{\rm tr}
+\Sigma_{\rm nt},
}
\]
\[
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\Sigma_\eta.
}
\]

So the coherent grouped defect decomposes into:

- a tracking sector \(\Sigma_{\rm tr}\),
- a genuine nontracking transfer-shape sector \(\Sigma_{\rm nt}\),
- and a dressing sector \(\Sigma_\eta\).

### 8.4 Final direct monomial invariants

**Status:** `Reduced / Controlled Reduction`

At first grouped weak-axisymmetric/reference-branch order, the three defect
coordinates are the first logarithmic drifts of three direct microscopic
monomials.

First define the reference-branch exponents
\[
\boxed{
E_*=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
}
\]
\[
\boxed{
F_*=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}
{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
}
\]

The final direct monomials are
\[
\boxed{
\mathfrak C_{{\rm tr},*}
:=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
:=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*},
}
\]
\[
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
}
\]

Their first log drifts are exactly
\[
\boxed{
\delta\ln \mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
\]

So at this reference-branch order the coherent zero-defect ledger is equivalent to
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln \mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln \mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]

This is the sharpest honest linearized invariant statement: zero defect is
equivalent, at first grouped weak-axisymmetric/reference-branch order, to
preservation of the three direct microscopic monomials above.

### 8.5 Similarity orbit and finite quotient closure

**Status:** `Exact Within Closure`

Inside the positive coherent microscopic sector
\[
\mathcal M_+
=
\bigl\{
(\lambda_W,c_{\eta U},\gamma,K_U,K_\eta^{(\mathrm{eff})},K_W^{(\mathrm{eff})},\mu_W,T_U)>0
\bigr\}.
\]
define the invariant map
\[
\boxed{
\mathcal I(x)=\bigl(\mathfrak C_{{\rm tr},*}(x),\mathfrak C_{{\rm nt},*}(x),\epsilon_\eta(x)\bigr).
}
\]

The monomial-drift map is
\[
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\[2pt]
\delta\ln\mathfrak C_{{\rm nt},*}\\[2pt]
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*\,\delta\mathbf x,
\]
with
\[
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\[4pt]
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\[4pt]
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
\]

On the constructive coherent branch,
\[
\det M_*^{(\tau,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0,
\]
so \(M_*\) has rank \(3\) and kernel dimension \(5\).

Choose free finite co-scalings
\[
(\Lambda,\ C,\ \Gamma,\ U,\ W)
\]
for
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}),
\]
and determine the remaining three by exact monomial preservation:
\[
K_\eta^{(\mathrm{eff})}\mapsto e^{\,2C-U}K_\eta^{(\mathrm{eff})},
\]
\[
T_U\mapsto
e^{\,U-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U)}T_U,
\]
\[
\mu_W\mapsto e^{\,M(\Lambda,C,\Gamma,U,W)}\mu_W,
\]
where
\[
M(\Lambda,C,\Gamma,U,W)
=
2C-U+2W-2\Lambda
-
E_*(2\Gamma+2\Lambda-U-W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Gamma+C-U).
\]

This defines the exact five-parameter similarity orbit \(\mathcal G_*\). The
finite quotient theorem is:
\[
\boxed{
\mathcal I(\widetilde x)=\mathcal I(x)
\iff
\widetilde x\in \mathcal G_*\!\cdot x.
}
\]

Equivalently,
\[
\boxed{
\mathcal M_+/\mathcal G_*\cong(\mathbb R_{>0})^3
\quad\text{with quotient coordinates}\quad
(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_\eta).
}
\]

So within the positive coherent microscopic sector, finite weak-axisymmetric
defect motion is exactly motion in this three-dimensional quotient. This is an
invariant-structure theorem for the reduced coherent hierarchy, not yet a theorem
that the full PDE dynamically remains on one similarity orbit.

### 8.6 Coherent/invariant status at this point

**Fixed within the coherent reduced hierarchy**

- the coherent local-kernel tracking reduction,
- the exact coherent placement map and support-enhancement factor,
- the stable-side reduced support-compensation theorem,
- the microscopic slippage ledger and triangular defect normal form,
- the final direct monomial invariant coordinates,
- and the exact similarity-orbit / quotient closure in the positive coherent sector.

**Still open**

- whether the actual moving-throat PDE realizes the coherent local-kernel
  hypothesis closely enough on the physical branch,
- whether the realized branch reaches the required support ratio
  \(\zeta_{\rm req}\) before softening,
- and whether the full PDE branch preserves the three reduced coherent invariants
  closely enough that the quotient description is the right physical branch
  language rather than only the right reduced one.

---

### 8.7 Branch-observable packet and transfer-shape / outgoing-prefactor compiler

**Status:** `Exact Within Closure`

The first exact PDE-facing observable packet is
\[
\boxed{
\Delta_{\rm obs}^{(1)}
:=
\begin{pmatrix}
\delta\ln R_{\rm tr}\\[2pt]
\delta\ln \mathfrak N_*\\[2pt]
\delta\ln \epsilon_\eta
\end{pmatrix},
\qquad
\mathfrak N_*:=\mathcal T^2R_{\rm tr}^{B_*},
\qquad
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
}
\]
The exact coefficient identity
\[
\boxed{A_{{\rm tr},*}=B_*C_{{\rm tr},*}}
\]
makes the tracking feed-through subtraction exact. The observable packet compiles
to the defect scalars as
\[
\boxed{\Theta_1=\delta\ln R_{\rm tr},}
\]
\[
\boxed{\Xi_1=\delta\ln\mathfrak N_*-B_*\,\delta\ln R_{\rm tr},}
\]
\[
\boxed{\mathcal R_1=-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta-\Xi_1.}
\]
So the first-order coherent zero-defect theorem is already
\[
\boxed{\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=\delta\ln\mathfrak N_*=\delta\ln\epsilon_\eta=0.}
\]

The exact direct coherent transfer shape is
\[
\boxed{
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}
=
\Lambda_0\,\frac{1-\epsilon_\eta}{R_{\rm target}},
\qquad
\Lambda_0:=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
}
\]
Its exact first-order identities are
\[
\boxed{\delta\ln\mathcal T^2=\Xi_1,}
\qquad
\boxed{\delta\ln(1-\epsilon_\eta)=\mathcal R_1+\Xi_1,}
\qquad
\boxed{\delta\ln R_{\rm target}=\mathcal R_1.}
\]

The grouped outgoing compiler is
\[
\boxed{K_0=P_0,}
\qquad
\boxed{K_2=P_2+\frac{a^2}{9c_s^2}P_0,}
\qquad
\boxed{K_4=P_4+\frac{a^2}{9c_s^2}P_2+\frac{4a^4}{81c_s^4}P_0,}
\]
\[
\boxed{\Gamma_5=\frac{a^5}{27c_s^5}P_0.}
\]
So the first grouped weak-axisymmetric scalar is already the outgoing-prefactor
slope:
\[
\boxed{\Xi_1=\frac{\delta\ln\mathcal T_A^2}{\epsilon\lambda_A}=\frac{P_1}{P_0}.}
\]

### 8.8 Direct-defect / dressing split and weak-axisymmetric filters

**Status:** `Exact Within Closure`

The direct transfer-shape defect depends only on the four microscopic slippages
\[
(\Sigma_Z,\Sigma_\chi,\Sigma_\epsilon,\Sigma_\delta),
\]
while the selected-branch dressing residual introduces only the additional
slippage
\[
\Sigma_\eta=\delta\ln\epsilon_\eta.
\]
The exact block-triangular split is
\[
\boxed{
\Theta_1=-C_{\rm tr}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1=-\frac{\epsilon_\eta}{1-\epsilon_\eta}\Sigma_\eta.
}
\]
So the full first-order rigidity theorem is
\[
\boxed{\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\Sigma_{\rm tr}=\Sigma_{\rm nt}=\Sigma_\eta=0.}
\]

Two exact filters are worth carrying forward.

First, the direct transfer shape is support-blind:
\[
\boxed{
\partial_\zeta\mathcal T^2
=
\partial_\zeta\mathfrak N_*
=
\partial_\zeta R_{\rm target}
=
\partial_\zeta\Xi_1=0.
}
\]
So coherent support enhancement can move the baseline loading, but it cannot
cancel the direct weak-axisymmetric transfer-shape defect.

Second, a pure grouped real `P2` anisotropy cannot linearly feed the scalar
off-bundle slippages. Its scalar feed-down begins only at quadratic order
through grouped invariants. So the remaining **linear** theorem gate must live in
the direct outlet / transfer-shape data themselves.

### 8.9 Reference-free final verdict packet

**Status:** `Exact Within Closure`

Let the target similarity orbit be
\[
\boxed{\mathcal O_*:=\mathcal G_*\!\cdot\mathbf x_*.}
\]
The full reduced verdict is now four scalars only. In additive form,
\[
\boxed{
\Delta_{\rm full}^{(x\mid \mathcal O_*)}
:=
\bigl(
\Delta_Q(x),
q_{\rm tr}^{(x\leftarrow \mathcal O_*)},
q_{\rm nt}^{(x\leftarrow \mathcal O_*)},
q_\eta^{(x\leftarrow \mathcal O_*)}
\bigr),
\qquad
\Delta_Q:=\chi_Q-1.
}
\]
The equivalent multiplicative chart is
\[
\boxed{
\mathcal V_{\rm full}^{(x\mid \mathcal O_*)}
:=
\bigl(
\chi_Q(x),
\mathfrak R_{\rm tr}^{(x\leftarrow \mathcal O_*)},
\mathfrak R_{\rm nt}^{(x\leftarrow \mathcal O_*)},
\mathfrak R_\eta^{(x\leftarrow \mathcal O_*)}
\bigr),
}
\]
and the mismatch chart is
\[
\boxed{
\mathcal M_{\rm full}^{(x\mid \mathcal O_*)}
:=
\bigl(
\chi_Q(x),
 m_T^{(x\leftarrow \mathcal O_*)},
 m_K^{(x\leftarrow \mathcal O_*)},
 m_\mu^{(x\leftarrow \mathcal O_*)}
\bigr).
}
\]
The exact conversion laws are
\[
\boxed{\mathfrak R_{\rm tr}=e^{q_{\rm tr}}=(m_T)^{1+\chi_{0,*}},}
\qquad
\boxed{\mathfrak R_\eta=e^{q_\eta}=\frac{1}{m_K},}
\]
\[
\boxed{\mathfrak R_{\rm nt}=e^{q_{\rm nt}}=\frac{m_\mu}{m_Km_T^{F_*}}.}
\]
So the exact reference-free full home-stretch theorem is
\[
\boxed{
\Delta_{\rm full}^{(x\mid \mathcal O_*)}=0
\iff
\chi_Q(x)=1\ \text{and}\ x\in\mathcal O_*.
}
\]
There is no privileged orbit representative left anywhere in the reduced endgame.

### 8.10 Canonical orbit projection and free-quintuple target graph

**Status:** `Exact Within Closure`

Let the positive microscopic state be
\[
\boxed{
\mathbf x:=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
}
\]
Carry forward the exact direct microscopic monomials
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
}
\]
\[
\boxed{\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.}
\]
Fix target values
\[
\mathfrak C_{{\rm tr},*}^{\rm target},\qquad
\mathfrak C_{{\rm nt},*}^{\rm target},\qquad
\epsilon_{\eta,*}^{\rm target},
\]
so the target orbit is
\[
\boxed{
\mathcal O_*=
\{x>0:\mathfrak C_{{\rm tr},*}=\mathfrak C_{{\rm tr},*}^{\rm target},\ \mathfrak C_{{\rm nt},*}=\mathfrak C_{{\rm nt},*}^{\rm target},\ \epsilon_\eta=\epsilon_{\eta,*}^{\rm target}\}.
}
\]

Define the free quintuple
\[
\boxed{\mathbf y:=(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).}
\]
Then the exact same-free-quintuple target graph is
\[
\boxed{
\delta_{U,*}^{\rm graph}(\mathbf y)
:=
\left[
\frac{\mathfrak C_{{\rm tr},*}^{\rm target}}{
\left(\dfrac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}}
\right]^{1/(1+\chi_{0,*})},
}
\]
\[
\boxed{T_U^{\rm graph}(\mathbf y)=\frac{L^2K_U}{\pi^2}\,\delta_{U,*}^{\rm graph}(\mathbf y),}
\]
\[
\boxed{K_{\eta,*}^{\rm graph}(\mathbf y)=\frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}^{\rm target}},}
\]
\[
\boxed{
\mu_W^{\rm graph}(\mathbf y)
=
\frac{\mathfrak C_{{\rm nt},*}^{\rm target}\,c_{\eta U}^2\,(K_W^{(\mathrm{eff})})^2}
{\epsilon_{\eta,*}^{\rm target}\,K_U\,\lambda_W^2}
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}}\right)^{-E_*}
\left(\delta_{U,*}^{\rm graph}(\mathbf y)\right)^{F_*}.
}
\]
So the target orbit is the exact five-dimensional graph
\[
\boxed{
\mathcal O_*=
\{\mathbf x_*^{\rm graph}(\mathbf y):\mathbf y\in(\mathbb R_{>0})^5\}.
}
\]

### 8.11 Graph-error packet and same-free-quintuple closure surface

**Status:** `Exact Within Closure`

The exact graph projection is simply “freeze the free quintuple and replace the
dependent triple by the graph target values”:
\[
\boxed{
\Pi_{\mathcal O_*}^{\rm graph}(\mathbf x)
:=
\mathbf x_*^{\rm graph}(\pi_{\rm free}(\mathbf x))
=
\Pi_{\mathcal O_*}^{\rm can}(\mathbf x).
}
\]
Define the dependent-triple graph errors
\[
\boxed{E_T(\mathbf x):=\ln\!\left(\frac{T_U}{T_U^{\rm graph}(\pi_{\rm free}(\mathbf x))}\right),}
\]
\[
\boxed{E_K(\mathbf x):=\ln\!\left(\frac{K_\eta^{(\mathrm{eff})}}{K_{\eta,*}^{\rm graph}(\pi_{\rm free}(\mathbf x))}\right),}
\]
\[
\boxed{E_\mu(\mathbf x):=\ln\!\left(\frac{\mu_W}{\mu_W^{\rm graph}(\pi_{\rm free}(\mathbf x))}\right).}
\]
These are exactly the Packet-B mismatch logs:
\[
\boxed{E_T=\frac{q_{\rm tr}}{1+\chi_{0,*}}=\ln m_T,}
\qquad
\boxed{E_K=-q_\eta=\ln m_K,}
\]
\[
\boxed{E_\mu=q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}=\ln m_\mu.}
\]
So the exact same-free-quintuple realization packet is
\[
\boxed{
\Delta_{\rm real}^{\rm graph}(\mathbf x\mid\mathcal Z_*)
=
\bigl(\chi_Q(\mathbf x)-1,\ E_T(\mathbf x),\ E_K(\mathbf x),\ E_\mu(\mathbf x)\bigr).
}
\]
The unique same-free-quintuple repair vector is simply the negative graph-error
triple on the dependent coordinates:
\[
\boxed{\Delta\mathbf x_{\rm rep}=(0,0,0,0,-E_K,0,-E_\mu,-E_T)^T.}
\]

### 8.12 Graph-slice theorem and one-scalar graph search

**Status:** `Exact Within Closure`

Define the graph-closure scalar
\[
\boxed{\widehat\chi_Q(\mathbf y):=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y)),}
qquad
\boxed{\widehat\Delta_Q(\mathbf y):=\widehat\chi_Q(\mathbf y)-1.}
\]
Since every graph point already lies on `\(\mathcal O_*\)`, the fully reduced
closure set is the codimension-one graph slice
\[
\boxed{
\mathcal Z_*
=
\{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\chi_Q(\mathbf y)=1\}
=
\{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\Delta_Q(\mathbf y)=0\}.
}
\]

For any graph-lifted free-quintuple family `\(\mathbf y(\tau)\)`, the Stage-175
monomial map annihilates the graph tangent:
\[
\boxed{M_*\,\dot{\Delta\mathbf x}_{\rm graph}=0.}
\]
So Packet B vanishes identically on graph-aligned families, and the home-stretch
search drops from four scalars to one scalar.

If a continuous graph-aligned one-parameter family has
\[
\widehat\Delta_Q(\tau_-)\,\widehat\Delta_Q(\tau_+)<0,
\]
then there exists `\(\tau_*\)` with
\[
\widehat\Delta_Q(\tau_*)=0,
\qquad
\mathbf x_{\rm graph}(\tau_*)\in\mathcal Z_*.
\]
So the first exact existence test on the graph is already a one-scalar crossing
theorem.

### 8.13 Explicit log-ray compiler and local scalar predictors

**Status:** `Exact Within Closure`

The smallest explicit one-parameter family on the free quintuple is the log-ray
\[
\boxed{\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s},}
\]
with base point `\(\mathbf y_\circ\)` and logarithmic direction
`\(\mathbf s=(s_\lambda,s_c,s_\gamma,s_U,s_W)\)`. Its graph lift stays on the
target orbit for all `\(\tau\)`.

The scalarized closure function on the ray is
\[
\boxed{\Phi_{\mathbf s}(\tau):=\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau)).}
\]
The first and second local directional packets are
\[
\Phi_0:=\Phi_{\mathbf s}(0),
\qquad
L_0:=\left.\frac{d}{d\tau}\ln\Phi_{\mathbf s}(\tau)\right|_{\tau=0},
\qquad
L_1:=\left.\frac{d^2}{d\tau^2}\ln\Phi_{\mathbf s}(\tau)\right|_{\tau=0}.
\]
Stage `187` gives the log-linear predictor
\[
\boxed{\tau_{\log}:=-\frac{\ln\Phi_0}{L_0},}
\]
while Stage `188` gives the quadratic log predictor through the discriminant
\[
\boxed{\Delta_{\log}:=L_0^2-2L_1\ln\Phi_0.}
\]
So the local graph-aligned search depends only on the free-ray direction and the
directional Hessian data of `\(\ln\widehat\chi_Q\)`.

### 8.14 Certified local ray screens and support-cardinality search

**Status:** `Exact Within Closure`

Stages `189-198` turn the free-quintuple local search into a certified sieve.
The primitive rays are controlled by the exact local packet
\[
(H_0,K_0,[\underline K_1,\overline K_1]),
\]
which determines admissibility, root brackets, and primitive ranking. Pairwise
mixed rays reduce to finite quartic stationary problems. Triple and quadruple
simplices add genuine interior content, but each closes to finite certified
interior packets and exact boundary splices.

By Stage `198`, the whole support-`<=4` search is already a finite certified
ledger with global interval splices and a preferred evaluation budget of
\[
\boxed{1140}
\]
exact candidate evaluations.

### 8.15 Final support-`<=5` local mixed-ray sieve

**Status:** `Exact Within Closure`

There are only five primitive free axes,
\[
\mathfrak I_5=\{\lambda,c,\gamma,U,W\},
\]
so the exact support-cardinality ceiling is `5`. The unique positive
five-coordinate simplex contributes the only genuinely new support-`5` interior
packet, and its entire boundary is exactly the already-finished support-`<=4`
search set.

Therefore the full local mixed-ray search reduces to one exact splice between the
support-`<=4` global ledger and the unique support-`5` interior packet:
\[
\boxed{
\tau_{\le 5,*}^{\rm best}
=
\min\bigl(\tau_{\le 4,*}^{\rm best},\ \tau_{5,*}^{\rm best,int}\bigr).
}
\]
The preferred exact total evaluation budget is
\[
\boxed{1464,}
\]
with fallback chart budget
\[
\boxed{2640.}
\]
So there is no remaining support-cardinality theorem to build. The next honest
move is to insert actual PDE-derived Hessian envelopes and branch data into this
completed ledger.


### 8.16 One-port same-charge mixed-kernel verdict and static/dynamic narrowing

**Status:** `Exact Within Closure`

The post-`201` same-charge audit inserts actual one-port wall/BdG/Maxwell/mixed
bundle data into the already-completed local mixed-ray ledger.

After the stable BdG support mode is integrated out,
\[
\boxed{
K_*:=K-\frac{C^2}{\varpi^2},
\qquad
\Delta:=\Omega_U^2\Omega_W^2-R^2,
}
\]
\[
\boxed{
Q:=G_U^2\Omega_W^2+2G_UG_WR+G_W^2\Omega_U^2,
\qquad
P:=\Omega_U^2G_W+RG_U,
}
\]
\[
\boxed{
D_0:=K_*-\frac{Q}{\Delta}.
}
\]
The exact reduced static `3\times 3` bundle has determinant
\[
\boxed{\det\mathcal K_{\rm red}=\Delta D_0.}
\]
So the natural admissible static branch is
\[
\boxed{
\Omega_U^2>0,
\qquad
\Delta>0,
\qquad
D_0>0.
}
\]

On that branch the exact static susceptibility shift is
\[
\boxed{
\delta V_{\rm mix}(r)
=
-\frac12 J(r)^T\mathcal K_{\rm red}^{-1}J(r)\le 0,
}
\]
so the induced same-charge static correction is always attractive or neutral at
quadratic order.

The wall–mixed bridge is already tied to the outgoing-prefactor chain:
\[
\boxed{
\Lambda:=\frac{P}{\Delta},
\qquad
N_0=\Lambda^2,
\qquad
P_0=\frac{N_0}{D_0},
}
\]
\[
\boxed{
\chi_{qW}=\frac{\Lambda}{D_0},
\qquad
\chi_{qW}^2=\frac{P_0}{D_0}.
}
\]
So static same-charge softening and outgoing quadrupole normalization are not
independent knobs inside the one-port bundle.

For the first primitive same-charge source families
\[
\mathcal S_Q(x)=\frac{1}{x^3},
\qquad
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
\]
the exact static mixed correction can generate only
\[
\boxed{
\frac{1}{x^6},
\qquad
\frac{e^{-2\kappa x}}{x^4},
\qquad
\frac{e^{-4\kappa x}}{x^2}.
}
\]
So the minimal static one-port mixed bundle creates **no new long-range
attractive law**. It survives only as short-range coefficient renormalization /
hybridization / near-softening.

The linear dynamic lift keeps the same narrowing. With
\[
K_B(\omega)
=
K-M\omega^2-\frac{C^2}{\varpi^2-\omega^2},
\qquad
D_\Pi(\omega)=K_B(\omega)-\frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)},
\]
the one-port mixed sector makes those same short-range families frequency
dependent, but it does **not** create a new conservative kernel class. The first
outgoing correction is phase-lag / pumping rather than linear barrier bypass.
After the residue/linewidth audit and the exact `5`PN compatibility surface,
the only linear dynamic corridor left is resonant dispersive enhancement of the
already-existing short-range families inside a finite branch-compatible window.

So the post-`201` same-charge bundle verdict is now:

- no new static long-range law,
- no new linear dynamic kernel class,
- one surviving short-range resonant corridor tied to the same bundle
  denominators that already appear in the `5`PN / `2.5`PN normalization chain.

### 8.17 Actual-branch same-charge ceiling, `\Xi_1`, and the static-first theorem

**Status:** `Exact Within Closure`

The later same-charge stages then convert the Stage-206 compatibility surface
into an actual-branch ceiling test.

For any prefactor ceiling `P_{\rm crit}`, the actual weak-axisymmetric branch is
carried by the grouped signature
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
and by the exact outgoing-prefactor slope
\[
\boxed{\Xi_1=\frac{P_1}{P_0}.}
\]
The branch prefactors are
\[
\boxed{
P_A=\bar P_0\bigl(1+\epsilon\lambda_A\Xi_1\bigr),
}
\]
so the grouped anomaly packet becomes
\[
\boxed{
a_{P_0}=\frac{\epsilon\bar P_0\Xi_1}{4},
\qquad
b_{P_0}=\frac{3\epsilon\bar P_0\Xi_1}{4},
}
\]
hence
\[
\boxed{b_{P_0}=3a_{P_0}.}
\]

Therefore the exact all-lane ceiling collapses to one scalar inequality:
\[
\boxed{
\bar P_0\bigl(1+|\epsilon\Xi_1|\bigr)\le P_{\rm crit}.
}
\]
Using the normalization compiler,
\[
\boxed{
\bar P_0=\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}},
}
\]
this becomes
\[
\boxed{
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}
\bigl(1+|\epsilon\Xi_1|\bigr)
\le P_{\rm crit}.
}
\]

So the actual weak-axisymmetric same-charge bottleneck is no longer generic
mixed-sector anisotropy. It is the explicit packet
\[
(\Delta_{\rm norm},a_{P_0},b_{P_0}),
\]
or, on the weak-axisymmetric line, just
\[
\boxed{(\Delta_{\rm norm},\Xi_1).}
\]

The strict `5`PN even-gate package, the pure-transfer sieve, the selected-branch
rigid-split classifier, and its continuum pullback then sharpen the verdict
further:

- the surviving first-order outlet is exactly the transported prefactor slope
  `\(\Xi_1\)`,
- the wall-like dynamic ceilings can be classified exactly along the selected
  branch,
- but on the carried branch family they are always weaker than the universal
  transported static ceilings.

So the first global kill condition is still
\[
\boxed{\text{the transported static }\Xi_1\text{ budget, not the wall-like dynamic window}.}
\]

### 8.18 Rigid-mouth physical normal form and Cartesian orbit lock

**Status:** `Exact Within Closure`

On the actual coherent branch, the physical direct observables are
\[
R_{\rm tr},
\qquad
\mathcal T^2,
\qquad
\epsilon_\eta,
\]
with exact identity
\[
\boxed{
R_{\rm target}\,\mathcal T^2
=
\Lambda_0(1-\epsilon_\eta),
\qquad
\Lambda_0:=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
}
\]

On the rigid-mouth slice `\(q_{\rm tr}=0\)`, the surviving finite packet is
already diagonal in the physical logarithmic chart
\[
\boxed{
U:=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V:=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
}
\]
with
\[
\boxed{
q_{\rm nt}=U,
\qquad
q_\eta=V.
}
\]

So the rigid-mouth branch is the exact Cartesian product of

1. a pure transfer-shape leg
   \[
   \mathcal T^2=\mathcal T_{\rm ref}^2 e^U,
   \qquad
   \epsilon_\eta=\epsilon_{\eta,\rm ref},
   \]
2. a pure dressing leg
   \[
   \mathcal T^2=\mathcal T_{\rm ref}^2,
   \qquad
   \epsilon_\eta=\epsilon_{\eta,\rm ref}e^V.
   \]

The exact physical-to-microscopic dependent-plane compiler is
\[
\boxed{
\begin{pmatrix}
\Delta_T\\
\Delta_{K_\eta}\\
\Delta_\mu
\end{pmatrix}
=
\begin{pmatrix}
0\\
-V\\
U-V
\end{pmatrix}.
}
\]
So:

- clearing the static transfer-shape defect changes only `\(\mu_W\)`,
- the post-static dressing correction is the equal
  `\(K_\eta^{(\mathrm{eff})}\)`–`\(\mu_W\)` shift,
- and the full orbit-restoring correction is the sum of those two pieces.

The rigid-mouth orbit-lock theorem is now completely sharp:
\[
\boxed{
U=0,\ V=0
\iff
\mathcal T^2=\mathcal T_{\rm ref}^2,\ \epsilon_\eta=\epsilon_{\eta,\rm ref}
\iff
R_{\rm target}=R_{\rm target,ref},\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]

The whole rigid-mouth packet and its microscopic correction are support-blind:
coherent support enhancement can move the support baseline, but it cannot alter
either the `(U,V)` packet or the orbit-restoring correction it requires.

### 8.19 Selected twin-support slice and primitive ranking

**Status:** `Exact Within Closure`

The passive/outgoing support-selection ambiguity is now fixed by the minimal
isotropic conservative quadrupole precursor
\[
Y_Q^{\rm cons}(\omega)
=
c_0+\frac{c_1}{1-\omega^2/\Omega_Q^2},
\qquad
c_0+c_1=1.
\]
On the natural isotropic branch
\[
\boxed{
c_0=\frac34,
\qquad
c_1=\frac14,
}
\]
so the exact support loading ratio is
\[
\boxed{
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\rho_\alpha-1=\frac13,
\qquad
\Pi_{\rm tr}=\frac43\,C_{\rm mix}.
}
\]

Equivalently, the selected branch is the exact one-parameter twin-support curve
\[
\boxed{
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma = \frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23,
}
\]
and it lies strictly inside the symmetric lowest-twin support window for the
whole allowed range. So the mixed-only and non-twin branches are gone from the
live same-charge closure.

On that selected curve only two primitive thresholds survive:
\[
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\qquad
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
\]
The exact primitive ranking is then
\[
\boxed{
\begin{aligned}
&0<\varrho<\varrho_{W\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_\Lambda > q_W > |q_U|,\\[1mm]
&\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_W > q_\Lambda > |q_U|,\\[1mm]
&\varrho_{U\Lambda}<\varrho<\frac23
&&\Longrightarrow&&
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
\end{aligned}
}
\]

So the quartic same-charge repair is **not** generically a large split-`U`
effect. Along the selected branch,

- `\(q_\chi\)` is always dominant,
- `\(q_Z\)` is always second,
- `\(q_W>|q_U|\)` everywhere,
- and the live quartic interpretation is generically
  interference / overlap / wall-blocking / outgoing-scale.

### 8.20 Actual twin-support placement, support/orbit split, and branch packet

**Status:** `Exact Within Closure`

On the actual coherent local D/N branch, the support coordinate that places the
branch on the selected twin-support curve is
\[
\boxed{
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
}
\]
So the realized selected point is
\[
\boxed{
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}=\frac{2\epsilon}{1-\epsilon}
=\frac{4}{3\varrho_{\rm phys}}-2.
}
\]
Once `\(\epsilon\)` is known, the Stage-224 ranking region is fixed with no
further support solve.

The coherent orbit packet can then be read directly from the physical placement
variables:
\[
\boxed{
q_{\rm tr}
=
(1+\delta_{U,*})\ln\frac{\chi_0}{\chi_{0,\rm ref}}
+
(1+\chi_{0,*})\ln\frac{\delta_U}{\delta_{U,\rm ref}},
}
\]
\[
\boxed{
q_{\rm nt}
=
\ln\frac{Z_W}{Z_{W,\rm ref}}
-
\ln\frac{\Lambda}{\Lambda_{\rm ref}}
+
E_*\ln\frac{\epsilon_W}{\epsilon_{W,\rm ref}}
-
F_*\ln\frac{\delta_U}{\delta_{U,\rm ref}},
}
\]
\[
\boxed{
q_\eta=\ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
}
\]

At infinitesimal order, the direct-observable compiler is
\[
\boxed{\Theta_1=d\ln R_{\rm tr},}
\]
\[
\boxed{
\Xi_1
=
-d\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,d\ln\epsilon_\eta,
}
\]
and the coherent orbit-lock theorem is
\[
\boxed{
q_{\rm tr}=q_{\rm nt}=q_\eta=0
\iff
d\ln R_{\rm tr}=d\ln R_{\rm target}=d\ln\epsilon_\eta=0.
}
\]
The outgoing finish line remains separate:
\[
\boxed{
N_Q=1,
\qquad
\text{equivalently }\chi_Q=1
\text{ on the natural source-map branch.}
}
\]

The exact support/orbit separation is now also explicit. The coherent support
ratio `\(\zeta\)` affects only the support packet
\[
\boxed{
\bigl(\zeta,\ M_{\rm mix},\ S(\zeta;\epsilon),\ M_{\rm tr}\bigr),
}
\]
while the orbit packet depends only on
\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda).
\]
So support compensation cannot rescue orbit-lock failure, and orbit lock does
not determine support enhancement.

The smallest exact front-end packet that the completed moving-throat branch now
needs to return is therefore
\[
\boxed{
\Bigl(
\epsilon,\,
\varrho_{\rm phys},\,
\sigma_{\rm phys},\,
\text{ranking region},\,
R_{\rm tr},\,
R_{\rm target},\,
\epsilon_\eta,\,
d\ln R_{\rm tr},\,
d\ln R_{\rm target},\,
d\ln\epsilon_\eta,\,
N_Q-1
\Bigr).
}
\]
This is the most concrete actual-branch evaluation checklist the compact program
currently has.


## 9. Relaxed-Constraint / Open-System Barrier Extension (`226-236`)

### 9.1 Relaxed branch declaration and exact recovery slice

**Status:** `Exact Within Closure`

The strict Stage-225 front end is now carried as the standard-recovery slice
\[
\boxed{
\mathfrak B_{\rm std}
:=
\Bigl\{
J^w=0,\ U=0,\ V=0,\ \varsigma(z)\equiv1
\Bigr\}
\subset
(\mathcal P_{225},\mathcal S_{225}).
}
\]

The post-225 relaxed branch is the codimension-three lift
\[
\boxed{
\mathfrak B_{226}^{\rm relax}
=
\Bigl\{
(\mathcal P_{225},\mathcal S_{225});
L_w,\ L_{UV},\ L_\varsigma
\Bigr\},
}
\]
with
\[
L_w=(S_{\rm leak},\mathcal W_w),
\qquad
L_{UV}=(U,V,\mathcal D_{UV}),
\qquad
L_\varsigma=(\varsigma,\varsigma_{\min},\mathrm{signchg}).
\]

The exact recovery slice is
\[
\boxed{
\ell_w=0,
\qquad
f_U=0,
\qquad
a=b=0
\iff
S_{\rm leak}=\mathcal W_w=U=V=\mathcal D_{UV}=0,
\quad
\varsigma\equiv1.
}
\]

The short-range/open-system firewall is unchanged on the relaxed branch.
Both the static and linear dynamic conservative mixed-sector corrections remain in
the exact kernel span
\[
\boxed{
\mathcal K_{\rm SR}
=
\mathrm{span}\!\left\{
\frac{1}{x^6},
\frac{e^{-2\kappa x}}{x^4},
\frac{e^{-4\kappa x}}{x^2}
\right\},
}
\]
so lifting `\(J^w\)`, lifting rigid-mouth factorization, or allowing a
compensated sign-changing source does **not** adjoin a new Coulomb-like or
Yukawa-`\(1/x\)` same-charge attraction.

The minimal relaxed packet is therefore
\[
\boxed{
\mathcal P_{226}^{\rm relax}
=
\Bigl(
\mathcal P_{225};
S_{\rm leak},
\mathcal W_w,
U,
V,
\mathcal D_{UV},
\varsigma(z),
\varsigma_{\min},
\mathrm{signchg}
\Bigr).
}
\]

### 9.2 Selected-branch leakage/work compiler

**Status:** `Exact Within Closure`

The selected-support demand carried from Stage 225 is
\[
\Pi_{\rm tr}
=
\frac43\,C_{\rm mix}
=
\frac{32\Lambda(1-\epsilon)}{3\pi^2}
=
\frac{16}{\pi^2}\Lambda\,\varrho,
\qquad
\varrho=\frac23(1-\epsilon).
\]

With the minimal selected-branch pullback
\[
\mathcal E_0(r)=\eta_{\rm leak}\,\Pi_{\rm tr}(r),
\]
the first exact open-system observables are
\[
\boxed{
S_{\rm leak}(r)
=
\frac{16\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{3\pi^{5/2}\lambda^3}
\,\Lambda(r)(1-\epsilon(r))
=
\frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
\,\Lambda(r)\varrho(r),
}
\]
\[
\boxed{
\mathcal W_w^{\rm bulk}(r)
=
\frac{512\sqrt2\,\eta_{\rm leak}^2\mu_w q\rho_0}{9\pi^{9/2}\lambda^3}
\,\Lambda(r)^2(1-\epsilon(r))^2,
}
\]
\[
\boxed{
\mathcal W_w^{\rm sess}(r)
=
\frac{2048\,\eta_{\rm leak}^2\mu_w q\rho_0}{9\pi^4\lambda^2}
\,\Lambda(r)^2(1-\epsilon(r))^2.
}
\]

This whole lane is support-side, not orbit-side:
\[
\partial_{R_{\rm tr}}S_{\rm leak}
=
\partial_{R_{\rm target}}S_{\rm leak}
=
\partial_{\epsilon_\eta}S_{\rm leak}
=
0,
\]
and likewise for both work scalars.

Its exact parity is
\[
S_{\rm leak}(-\eta_{\rm leak})=-S_{\rm leak}(\eta_{\rm leak}),
\qquad
\mathcal W_w^{\rm bulk/sess}(-\eta_{\rm leak})
=
\mathcal W_w^{\rm bulk/sess}(\eta_{\rm leak}).
\]

A practical packet is
\[
\boxed{
\mathcal P_{227}^{\rm leak}
=
\Bigl(
\mathcal P_{225};
\Lambda,\epsilon,\eta_{\rm leak},
S_{\rm leak},
\mathcal W_w^{\rm bulk},
\mathcal W_w^{\rm sess}
\Bigr).
}
\]

### 9.3 Non-rigid mouth/dressing packet

**Status:** `Exact Within Closure`

The orbit-side non-rigid packet is generated by the reduced free energy
\[
\boxed{
\mathcal F_{\rm nr}(U,V;r)
=
\frac12 a_U U^2
+
\frac12 a_V V^2
-
\chi_{UV}(r)\,UV
-
f_U(r)\,U.
}
\]

Define
\[
\boxed{
\Delta_{UV}:=a_U a_V-\chi_{UV}^2.
}
\]
Admissibility is exactly
\[
\boxed{\Delta_{UV}>0.}
\]

The stationary packet is
\[
\boxed{
U(r)=\frac{a_V f_U(r)}{\Delta_{UV}(r)},
\qquad
V(r)=\frac{\chi_{UV}(r)f_U(r)}{\Delta_{UV}(r)},
\qquad
\frac{V}{U}=\frac{\chi_{UV}}{a_V}.
}
\]

The exact finite physical compiler is
\[
\boxed{
\mathcal T^2(r)
=
\mathcal T_{\rm ref}^2 e^{U(r)},
\qquad
\epsilon_\eta(r)
=
\epsilon_{\eta,\rm ref}e^{V(r)},
}
\]
\[
\boxed{
\frac{R_{\rm target}(r)}{R_{\rm target,ref}}
=
\frac{
1-\epsilon_{\eta,\rm ref}e^{V(r)}
}{
1-\epsilon_{\eta,\rm ref}
}
e^{-U(r)}.
}
\]

The dependent rigid-mouth plane therefore induces the exact microscopic
correction
\[
\boxed{
\mathbf y_{\rm nr}^{\rm dep}
=
\begin{pmatrix}
0\\
-V\\
U-V
\end{pmatrix},
}
\]
and the induced transfer-to-dressing drain is
\[
\boxed{
\mathcal D_{UV}
=
\chi_{UV}UV
=
\frac{\chi_{UV}^2 a_V f_U^2}{\Delta_{UV}^2}
\ge 0.
}
\]

At weak-axisymmetric order,
\[
U=\varepsilon u_1+O(\varepsilon^2),
\qquad
V=\varepsilon v_1+O(\varepsilon^2),
\]
the exact first-order packet is
\[
\boxed{\Xi_1^{\rm nr}=u_1,}
\]
\[
\boxed{
\mathcal R_1^{\rm nr}+\Xi_1^{\rm nr}
=
-\frac{\epsilon_{\eta,\rm ref}}{1-\epsilon_{\eta,\rm ref}}\,v_1,
}
\qquad
\boxed{
\mathcal R_1^{\rm nr}
=
-u_1
-\frac{\epsilon_{\eta,\rm ref}}{1-\epsilon_{\eta,\rm ref}}\,v_1.
}
\]
So the rigid-mouth front-end pair `\((\Delta_{\rm norm},\Xi_1)\)` is the
projection of a larger two-leg packet once the dressing leg is activated.

### 9.4 Compensated multimode source compiler

**Status:** `Exact Within Closure`

The first explicit compensated source branch is
\[
\boxed{
\varsigma(x;r)
=
1+a(r)\cos(\pi x)+b(r)\cos(2\pi x),
}
\qquad
\int_0^1 \varsigma(x;r)\,dx = 1.
\]

Its two carried source moments remain exact:
\[
\boxed{
\mathfrak g(a,b)
=
\frac{2}{\pi}\left(1+\frac{a}{3}-\frac{b}{15}\right),
}
\]
\[
\boxed{
\mathcal S(a,b)
=
\frac{2\tanh(\pi/2)}{\pi}
\left(1+\frac{a}{5}+\frac{b}{17}\right).
}
\]

With normalized moments
\[
\widetilde g:=\frac{\pi}{2}\mathfrak g-1,
\qquad
\widetilde S:=\frac{\pi}{2\tanh(\pi/2)}\mathcal S-1,
\]
the two-mode source family is exactly invertible:
\[
\boxed{
a=\frac{85}{42}\widetilde S+\frac{25}{14}\widetilde g,
\qquad
b=\frac{425}{42}\widetilde S-\frac{85}{14}\widetilde g.
}
\]

The carried mixed-to-shell loading law extends as
\[
\boxed{
\mathcal R(a,b)
=
\frac{(\mathfrak g(a,b)-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2}.
}
\]

For the transported Session-I branch,
\[
a(r)=a_0 s(r),
\qquad
b(r)=b_0 s(r),
\qquad
s(r)=\frac{r_\sigma^2}{r^2+r_\sigma^2},
\]
the sign-change threshold on the sampled orientation `\(a_0>0,\ b_0<0\)` is
\[
\boxed{
\mathrm{signchg}(r)
\iff
s(r)>\frac{1}{a_0-b_0}
\iff
r<r_\sigma\sqrt{a_0-b_0-1}.
}
\]

The compact packet to carry forward is
\[
\boxed{
\mathcal M_\sigma^{\rm pack}(r)
=
\bigl(
\varsigma(x;r),\,
\varsigma_{\min}(r),\,
\mathrm{signchg}(r),\,
\mathfrak g(r),\,
\mathcal S(r),\,
\mathcal R(r)
\bigr).
}
\]

### 9.5 Carried Stage-230 stationary lowered front end

**Status:** `Reduced / Controlled Reduction`

The extension does not restate a standalone Stage-230 theorem block, but Stage 231
carries forward its exact reduced stationary front end as
\[
\boxed{
V_{\rm eff}^{(230)}(r)
=
V_{\rm short}^{(1p)}(r)
-
\lambda_L S_{\rm leak}(r)
-
\lambda_W\mathcal W_w^{\rm sess}(r)
-
\Delta E_{UV}(r)
-
\mathcal M_{\sigma,\rm red}(r).
}
\]

Here the compact master writes

- `\(\Delta E_{UV}(r)\)` for the stationary non-rigid lowering,
  \[
  \boxed{
  \Delta E_{UV}(r)
  =
  -\mathcal F_{\rm nr}\bigl(U(r),V(r);r\bigr)
  =
  \frac{a_V f_U(r)^2}{2\Delta_{UV}(r)},
  }
  \]
- and `\(\mathcal M_{\sigma,\rm red}(r)\)` for the scalar barrier-lowering
  contribution built from the Stage-229 source packet
  `\(\mathcal M_\sigma^{\rm pack}(r)\)`.

This is only a short-range/open-system reshaping of the already-audited one-port
same-charge front end. It does not change the kernel-span theorem in `§9.1`.

### 9.6 Dynamic event-chain, thresholds, and WKB

**Status:** `Exact Within Closure`

On the carried lowered front end `\(V(r):=V_{\rm eff}^{(230)}(r)\)`, the reduced
one-dimensional dynamics is
\[
\boxed{
m_s\ddot r=-V'(r),
\qquad
E=\frac{m_s}{2}\dot r^{\,2}+V(r).
}
\]

If the lowered branch has one peak
\[
V'(r_{\rm peak})=0,
\qquad
V''(r_{\rm peak})<0,
\qquad
V_{\rm peak}=V(r_{\rm peak}),
\]
then the finite-radius classical threshold speed is
\[
\boxed{
v_{\rm crit,new}
=
\sqrt{\frac{2}{m_s}\bigl(V_{\rm peak}-V(r_0)\bigr)}.
}
\]

For the pure Coulomb comparison,
\[
V_{\rm Coul}(r)=\frac{1}{r},
\]
the contact threshold at the same launch radius is
\[
\boxed{
v_{\rm contact,Coul}
=
\sqrt{\frac{2}{m_s}\left(\frac1{r_{\rm contact}}-\frac1{r_0}\right)}.
}
\]

Hence the exact classical window theorem is
\[
\boxed{
v_{\rm crit,new}<v_0<v_{\rm contact,Coul}
\Longrightarrow
\text{lowered branch reaches contact while Coulomb still turns back.}
}
\]

For a subbarrier energy `\(V(r_0)<E<V_{\rm peak}\)`, the turning points are
\[
\boxed{
V\bigl(r_-(E)\bigr)=E,
\qquad
V\bigl(r_+(E)\bigr)=E,
}
\]
with exact transport law
\[
\boxed{
\frac{dr_\pm}{dE}=\frac{1}{V'(r_\pm(E))}.
}
\]

The WKB action and transmission factor are
\[
\boxed{
I_{\rm new}(E)
=
\frac{1}{\hbar_{\rm eff}}
\int_{r_-(E)}^{r_+(E)}
\sqrt{2m_s\bigl(V(r)-E\bigr)}\,dr,
\qquad
T_{\rm new}(E)=e^{-2I_{\rm new}(E)}.
}
\]

Against the pure-Coulomb reference,
\[
\boxed{
\frac{T_{\rm new}(E)}{T_{\rm Coul}(E)}
=
\exp\!\bigl[-2\bigl(I_{\rm new}(E)-I_{\rm Coul}(E)\bigr)\bigr].
}
\]

The turning-point event chain also carries two reduced diagnostics:
\[
\boxed{
\Xi_{\rm turn}(E)=\Xi_1\bigl(r_+(E)\bigr),
}
\qquad
\boxed{
\lambda_{\rm th}(E)
=
\left|\frac{E}{V'(r_+(E))}\right|.
}
\]

Near the peak,
\[
V(r)
=
V_{\rm peak}
-\frac{K_{\rm peak}}{2}(r-r_{\rm peak})^2+\cdots,
\qquad
K_{\rm peak}:=-V''(r_{\rm peak})>0,
\]
so
\[
\boxed{
I_{\rm top}(E)
=
\frac{\pi\Delta E}{\hbar_{\rm eff}}
\sqrt{\frac{m_s}{K_{\rm peak}}}
+O(\Delta E^{3/2}),
\qquad
\Delta E:=V_{\rm peak}-E.
}
\]

### 9.7 Conditional helicity-export diagnostic

**Status:** `Exact Within Closure`

The hidden-channel diagnostic carried by the parent projection algebra is
\[
\boxed{
\frac{dH_{\rm sub}}{dt}
=
-\Phi_h-2\mathcal C_h,
}
\]
where `\(H_{\rm sub}\)` is the unresolved projected helicity on the Stage-231
event chain.

Introduce the orientation closure label
\[
\sigma\in\{+1,-1\},
\]
and write
\[
\boxed{
\dot H_{{\rm sub},\sigma}(t;E)
=
\Gamma_0(t;E)+\sigma\,\Gamma_1(t;E)
=
\Gamma_0(t;E)\bigl[1+\sigma\,\alpha_h(t;E)\bigr],
}
\]
with
\[
\alpha_h:=\frac{\Gamma_1}{\Gamma_0}.
\]

Then the exact instantaneous and integrated comparison ratios are
\[
\boxed{
R_{\rm inst}(t;E)
=
\frac{1+\alpha_h(t;E)}{1-\alpha_h(t;E)},
\qquad
\alpha_h=\frac{R_{\rm inst}-1}{R_{\rm inst}+1},
}
\]
\[
\boxed{
R_{\rm int}(E)
=
\frac{1+\bar\alpha_h(E)}{1-\bar\alpha_h(E)},
\qquad
\bar\alpha_h=\frac{R_{\rm int}-1}{R_{\rm int}+1}.
}
\]

So within the declared mixed/vortical closure the aligned-vs-anti-aligned
preference is a one-scalar asymmetry problem. This remains a diagnostic
unloading theorem, not yet a microscopic spin theorem.

### 9.8 Crossing-vs-collapse / Goldilocks compiler

**Status:** `Exact Within Closure`

From the carried event chain,
\[
\dot r
=
-\sqrt{\frac{2}{m_s}\bigl(E-V(r)\bigr)},
\]
the exact path-time is
\[
\boxed{
T_{\rm traj}(E;r_a\to r_b)
=
\sqrt{\frac{m_s}{2}}
\int_{r_b}^{r_a}\frac{dr}{\sqrt{E-V(r)}}.
}
\]

The Stage-233 characteristic-width closure compresses this to
\[
\boxed{
t_{\rm cross}(E)
=
\lambda_{\rm eff}
\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
}
\]

On the active dressing leg, the local unstable normal form
\[
\mu_\eta\ddot V=g_{UV}\chi_{\rm peak}V
\]
gives the undamped collapse time
\[
\boxed{
t_{\rm collapse}
=
\sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
}
\]

Hence the survival ratio is
\[
\boxed{
\mathcal S(E)
=
\frac{t_{\rm cross}(E)}{t_{\rm collapse}}
=
\lambda_{\rm eff}
\sqrt{\frac{m_s g_{UV}\chi_{\rm peak}}{2\mu_\eta(E-V_{\rm peak})}}.
}
\]

The exact lower-edge theorem is
\[
\boxed{
\mathcal S(E)<1
\iff
E>E_{\rm edge}
:=
V_{\rm peak}
+
\frac{m_s}{\mu_\eta}
\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
}
\]

At fixed launch radius,
\[
E=\frac{m_s}{2}v_0^2+V(r_0),
\]
so the speed-space compiler becomes
\[
\boxed{
v_{\rm safe,min}^2
=
v_{\rm crit,new}^2
+
\frac{\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}}{\mu_\eta}.
}
\]

On heavy-throat scaling
\[
\mu_\eta=\alpha m_s,
\]
the explicit mass dependence cancels:
\[
\boxed{
E_{\rm edge}
=
V_{\rm peak}
+
\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
}
\]

So the wall-timescale story is controlled by
\[
\lambda_{\rm eff},\qquad
\chi_{\rm peak},\qquad
g_{UV},\qquad
\mu_\eta/m_s,
\]
not by mass scaling alone.

### 9.9 Microscopic cubic/quintic export kernel and safe surface

**Status:** `Exact Within Closure`

The first active-leg microscopic odd export channels are

- a derivative-coupled scalar outlet, giving
  \[
  \boxed{
  \Gamma_3
  =
  \Pi_{V0}^2\,
  \gamma_1\eta_0^2\frac{\Omega_{U,0}^4}{\Delta_0^2},
  }
  \]
- and the selected mixed quadrupole outlet, giving
  \[
  \boxed{
  \Gamma_5
  =
  \Pi_{V-}^2\frac{a^5}{27c_s^5}P_{0,-}.
  }
  \]

So the first microscopic export kernel is
\[
\boxed{
\Sigma_{\rm exp}(\omega)
=
-i\Gamma_3\omega^3-i\Gamma_5\omega^5+O(\omega^7),
}
\]
equivalently
\[
\boxed{
K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5.
}
\]

The homogeneous active-leg equation is
\[
\boxed{
\mu_\eta\ddot V-\kappa_V V+\Gamma_3 V^{(3)}-\Gamma_5 V^{(5)}=0,
\qquad
\kappa_V:=g_{UV}\chi_{\rm peak}.
}
\]

Its odd power input is exact:
\[
\boxed{
\dot V F_{\rm odd}
=
\frac{d}{dt}\mathcal S_{\rm odd}
-\Gamma_3\ddot V^{\,2}
-\Gamma_5\dddot V^{\,2},
}
\]
so the exported power after Schott subtraction is
\[
\boxed{
\mathcal P_{\rm exp}
=
\Gamma_3\ddot V^{\,2}+\Gamma_5\dddot V^{\,2}\ge0.
}
\]

For exponential growth `\(V\sim e^{st}\)`, the characteristic polynomial is
\[
\boxed{
F(s)=\Gamma_5 s^5+\Gamma_3 s^3+\mu_\eta s^2-\kappa_V.
}
\]
Because `\(F'(s)>0\)` for `\(s>0\)`, there is exactly one positive real growth
root for any finite passive `\(\Gamma_3,\Gamma_5\ge0\)`. So the minimal
microscopic kernel has no finite unconditional analogue of the Session-IV
envelope threshold.

What it **does** have is the exact event-safe surface
\[
\boxed{
\Gamma_3 s_c^3+\Gamma_5 s_c^5
\ge
\mu_\eta(s_0^2-s_c^2),
\qquad
s_0:=\sqrt{\kappa_V/\mu_\eta},
\qquad
s_c:=\frac{1}{t_{\rm cross}}.
}
\]
In normalized form,
\[
\boxed{
\widehat\Gamma_3+s_c^2\widehat\Gamma_5
\ge
\frac{s_0^2-s_c^2}{s_c^3},
\qquad
\widehat\Gamma_n:=\Gamma_n/\mu_\eta.
}
\]

### 9.10 Vacuum/lattice heat partition and cold-survival energy

**Status:** `Exact Within Closure`

Split the microscopic coefficients as
\[
\Gamma_3=\Gamma_3^{\rm vac}+\Gamma_3^{\rm lat},
\qquad
\Gamma_5=\Gamma_5^{\rm vac}+\Gamma_5^{\rm lat}.
\]

For an arbitrary event window,
\[
\boxed{
E_{\rm vac}(T)=\Gamma_3^{\rm vac}\mathcal I_2(T)+\Gamma_5^{\rm vac}\mathcal I_3(T),
}
\]
\[
\boxed{
E_{\rm lat}(T)=\Gamma_3^{\rm lat}\mathcal I_2(T)+\Gamma_5^{\rm lat}\mathcal I_3(T),
}
\]
with
\[
\mathcal I_2(T)=\int_0^T \ddot V^{\,2}dt,
\qquad
\mathcal I_3(T)=\int_0^T \dddot V^{\,2}dt.
\]

Define the exact shape quotient
\[
\boxed{
\mathfrak r_V:=\frac{\mathcal I_3}{\mathcal I_2}.
}
\]
Then the heat partition fractions are
\[
\boxed{
f_{\rm vac}(\mathfrak r_V)
=
\frac{\Gamma_3^{\rm vac}+\Gamma_5^{\rm vac}\mathfrak r_V}
{\Gamma_3+\Gamma_5\mathfrak r_V},
\qquad
f_{\rm lat}(\mathfrak r_V)
=
\frac{\Gamma_3^{\rm lat}+\Gamma_5^{\rm lat}\mathfrak r_V}
{\Gamma_3+\Gamma_5\mathfrak r_V}.
}
\]

On the single-growth cold branch `\(V(t)=V_{\rm in}e^{st}\)`,
\[
\boxed{
\mathfrak r_V=s^2,
}
\]
and the event-equivalent rates are
\[
\boxed{
\gamma_{\rm vac}^{\rm eq}(s)=\Gamma_3^{\rm vac}s^2+\Gamma_5^{\rm vac}s^4,
\qquad
\gamma_{\rm lat}^{\rm eq}(s)=\Gamma_3^{\rm lat}s^2+\Gamma_5^{\rm lat}s^4.
}
\]

On the exact cold-survival edge,
\[
\Gamma_3 s_c^3+\Gamma_5 s_c^5=\mu_\eta(s_0^2-s_c^2),
\]
the minimum total exported energy is
\[
\boxed{
E_{\rm exp,min}^{\rm safe}
=
\frac{V_{\rm in}^2}{2}(e^2-1)\,\mu_\eta(s_0^2-s_c^2),
}
\]
with channel-resolved shares
\[
\boxed{
E_{\rm vac/lat,min}^{\rm safe}
=
f_{\rm vac/lat}(s_c)\,
\frac{V_{\rm in}^2}{2}(e^2-1)\,\mu_\eta(s_0^2-s_c^2).
}
\]

The safe-edge event-equivalent lattice rate is therefore
\[
\boxed{
\gamma_{\rm lat,safe}^{\rm eq}
=
f_{\rm lat}(s_c)\,\mu_\eta\frac{s_0^2-s_c^2}{s_c}.
}
\]

The old Session-IV `\(3:1\)` split is now an exact coefficient surface:
\[
\boxed{
f_{\rm lat}(s_c)=\frac34
\iff
\Gamma_3^{\rm lat}+\Gamma_5^{\rm lat}s_c^2
=
3\bigl(\Gamma_3^{\rm vac}+\Gamma_5^{\rm vac}s_c^2\bigr).
}
\]
A speed-independent `\(3:1\)` family exists iff the same fraction holds in both
channels,
\[
\Gamma_3^{\rm lat}=\frac34\Gamma_3,
\qquad
\Gamma_5^{\rm lat}=\frac34\Gamma_5.
\]

### 9.11 Physical calibration and material-screening companion

**Status:** `Reduced / Controlled Reduction`

Use the physical calibration dictionary
\[
\boxed{
t^{\rm phys}=t_*\,t,
\qquad
r^{\rm phys}=\frac{\lambda_{\rm phys}}{\lambda_{\rm ref}}\,r,
\qquad
E^{\rm phys}=E_*\,E.
}
\]

The exact lattice-turnover threshold is then
\[
\boxed{
(\lambda_{\rm ep}\omega_D)_{\min}^{(236)}
=
\frac{\gamma_{\rm lat,safe}^{\rm eq}}{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,t_*},
}
\]
or equivalently
\[
\boxed{
(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys})_{\min}^{(236)}
=
\frac{\gamma_{\rm lat,safe}^{\rm eq}}{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,s_c}.
}
\]

At the turning point,
\[
\boxed{
\chi_{\lambda,\rm lattice}(r_{\rm turn})
=
\frac{2\lambda_{\rm ref}}{r_{\rm turn}}.
}
\]
This is a pure geometry ratio; it does **not** determine a stiffness by itself.

The force-matched stiffness compiler is
\[
\boxed{
k_{\rm eff,req}
=
\mathcal K_{\rm turn}\,\frac{E_*}{\lambda_{\rm phys}^2},
\qquad
\mathcal K_{\rm turn}
:=
\frac{\lambda_{\rm ref}^2}{r_{\rm turn}}|V_{\rm red}'(r_{\rm turn})|.
}
\]

The Korringa-limited spin-survival ceiling is
\[
\boxed{
T_{\max}
=
\frac{\mathcal K_{\rm corr}}{t_{\rm cross}^{\rm phys}}
=
\frac{s_c\,\mathcal K_{\rm corr}}{t_*}.
}
\]

Stage 236 then packages the condensed-matter screen into four exact ratios:
\[
\boxed{
\Pi_{\rm ep}
:=
\frac{\Upsilon_{\rm lat}\,\zeta_{\rm ep}\,\lambda_{\rm ep}\,\omega_D\,t_*}
{\gamma_{\rm lat,safe}^{\rm eq}}
=
\frac{\lambda_{\rm ep}\omega_D}{(\lambda_{\rm ep}\omega_D)_{\min}^{(236)}},
}
\]
\[
\boxed{
\Pi_\chi
:=
\chi_{\lambda,\rm lattice}(r_{\rm turn}),
}
\]
\[
\boxed{
\Pi_k
:=
\frac{k_{\rm eff}\,\lambda_{\rm phys}^2}{\mathcal K_{\rm turn}E_*}
=
\frac{k_{\rm eff}}{k_{\rm eff,req}},
}
\]
\[
\boxed{
\Pi_T(T)
:=
\frac{\mathcal K_{\rm corr}}{T\,t_{\rm cross}^{\rm phys}}
=
\frac{T_{\max}}{T}.
}
\]

A candidate host survives the Stage-236 companion screen only if
\[
\boxed{
\Pi_{\rm ep}\ge1,
\qquad
\Pi_\chi\ge1,
\qquad
\Pi_k\ge1,
\qquad
\Pi_T(T)\ge1.
}
\]

### 9.12 Global reading of the `226-236` extension

**Status:** `Open`

The post-225 extension does **not** overturn the earlier same-charge verdict.
It leaves intact:

- the no-new-long-range-kernel theorem,
- the rigid-mouth / selected-support front-end packet,
- the Packet-A finish line and actual-prefactor ceiling,
- and the back-end graph/orbit compiler.

What it does add is a controlled companion branch around that front end:

1. a codimension-three relaxed lift
   \[
   (J^w,\ U/V,\ \varsigma),
   \]
2. an explicit lowered stationary front end `\(V_{\rm eff}^{(230)}\)`,
3. a complete reduced event chain with threshold, turning-point, WKB, and
   Goldilocks compilers,
4. a hidden-channel unloading diagnostic,
5. a microscopic cubic/quintic export kernel with exact event-safe surface,
6. a channel-resolved heat-partition law,
7. and a four-inequality material-screening companion.

So the extension should be read as a **short-range/open-system barrier and
materials companion**, not as a replacement for the strict selected-branch
Packet-A / Packet-B realization program.

On that companion branch the live open objects are now the realized coefficients
and calibration data
\[
\eta_{\rm leak},
\qquad
f_U,\ \chi_{UV},
\qquad
(a,b)\ \text{or}\ (\mathfrak g,\mathcal S,\mathcal R),
\qquad
\Gamma_3,\Gamma_5,
\]
\[
\Gamma_3^{\rm vac/lat},\ \Gamma_5^{\rm vac/lat},
\qquad
t_*,\ \lambda_{\rm phys},\ E_*,\ \Upsilon_{\rm lat},\ \mathcal K_{\rm corr},
\]
rather than more missing symbolic compilers.


## 10. Final Theorem Ledger

### 10.1 Exact parent-theory statements

**Status:** `Exact`

- The bulk arena, strict parent fields `psi,A_M`, charge ontology, GNLS action,
  localized Maxwell action, GNLS equation, continuity law, Maxwell equation,
  Madelung rewrite, and mixed gauge-invariant identities in Section `2` are the
  exact starting theory of this program.
- `Sigma/R` is exact as a confinement-coupling argument. It becomes an
  autonomous parent dynamical field only after `S_Sigma` is included; otherwise
  the wall PDE in Sections `3-4` is an effective closure.
- Nothing in the later reduced hierarchy changes the strict parent equations;
  all later reductions are branch-specific or closure-specific consequences.

### 10.2 Exact results within explicit reduced closures

**Status:** `Exact Within Closure`

- The finite-throat D/N support problem, the minimal distributed wall action, and
  the grouped conservative lane coefficients are exact inside the reduced
  moving-throat wall/support hierarchy.
- The one-lane positive localized-source theorem, the explicit exponential
  mouth-layer law, the exact Family-1 gain map, the self-consistent explicit
  mouth branch, and the singular exclusion of the equal-normalized branch are
  exact inside the stated positive-mouth closure.
- The co-evolving Family-1 fixed-point map, the equivalence
  `\(\mathfrak g=\mathfrak g_*\iff \mathcal R=1/4\)`, and the mouth-slope identity
  are exact inside the analyzed co-evolving reduced closure.
- The coherent local-kernel tracking reduction, the coherent placement map, the
  support-enhancement factor, the microscopic slippage normal form, and the
  finite similarity-orbit / quotient theorem are exact inside the coherent
  reduced hierarchy and the positive coherent microscopic sector where stated.
- The post-175 extension is now also exact inside the carried hierarchy:
  the isotropic grouped-real `P2` conservative one-pole surface, the scalar/
  geometry firewall, the exact outgoing `l=2` DtN fingerprint, the source/outgoing
  factorization, higher-odd irrelevance, the Packet-A finish-line theorem
  `\(\chi_Q=1\)`, the reference-free four-scalar verdict packet, the canonical orbit
  projection, the free-quintuple target graph, the graph-error packet, the
  graph-slice theorem, the explicit log-ray compiler, and the finite support-`<=5`
  local mixed-ray sieve.
- The post-201 same-charge / rigid-mouth / selected-branch extension is now also
  exact inside the carried one-port, weak-axisymmetric, rigid-mouth, and
  selected-branch closures:
  the one-port static susceptibility kernel and square-law suppression theorem,
  the dynamic mixed-port / phase-lag no-go and resonant-survival window,
  the actual-branch prefactor-ceiling compiler in `\((\Delta_{\rm norm},\Xi_1)\)`,
  the selected-branch static-first theorem, the rigid-mouth physical normal form
  `(U,V)`, the physical-to-microscopic dependent compiler, the exact support
  slice `\(\rho_\alpha=4/3\)`, the selected twin-support ranking theorem, and the
  actual placement / support-versus-orbit packet compiler.
- The post-225 relaxed/open-system barrier companion is now also exact within the
  carried reduced hierarchy where stated:
  the codimension-three relaxed-branch declaration and recovery slice,
  the selected-branch leakage/work compiler, the finite non-rigid `U/V`
  compiler, the compensated two-mode source compiler, the dynamic
  event-chain/turning-point/WKB compiler, the conditional unresolved-helicity
  export diagnostic, the Goldilocks lower-edge theorem, the microscopic
  cubic/quintic export kernel with exact event-safe surface, the channel-resolved
  vacuum/lattice heat-partition law, and the Stage-236 material-screening
  companion invariants.

### 10.3 Numerically located branch data

**Status:** `Numerically Located`

- The explicit Family-1 canonical mouth point is
  \[
  \Pi_*\approx 1.50882951349316,
  \qquad
  \widehat T_{m,*}\approx 0.901484054174205.
  \]
- The mouth-only first-order/full-profile correction gives the shifted comparison
  point
  \[
  \Pi_{\rm corr}\approx 2.41591392833607,
  \qquad
  \widehat T_{m,\rm corr}\approx 1.17313803363654,
  \]
  with a one-step nonlinear probe at
  \[
  \Pi_1\approx 2.53914847609768,
  \qquad
  \widehat T_{m,1}\approx 1.21036942084359.
  \]
- At the old canonical traction, the co-evolving reduced fixed point survives but
  shifts to
  \[
  \mathfrak g_{\rm fp}\approx 0.693352419668063,
  \qquad
  \mathcal R_{\rm fp}\approx 0.2827139049082381,
  \qquad
  \Pi_{\rm fp}\approx 1.4885734438300713.
  \]
- Exact reduced compensation is restored only on the renormalized co-evolving
  branch
  \[
  \Sigma_0^{\rm can}\approx 4.651033550168876,
  \qquad
  \widehat T_{m,\rm can}\approx 1.446708366456762,
  \qquad
  \Pi_{\rm can}\approx 3.8715643774790087.
  \]

### 10.4 Reduced endpoint statements now fixed by the updated compiler

**Status:** `Reduced / Controlled Reduction`

- The grouped real `P2` conservative front end still has the exact isotropic
  one-pole target surface
  \[
  a_2=b_2=a_4=b_4=0,
  \qquad
  \Delta_{\rm pole}=0,
  \qquad
  \widehat Y_Q^{\rm cons}(\omega)=\frac34+\frac14(1-\omega^2/\Omega_Q^2)^{-1}.
  \]
- The entire Packet-A retarded finish line is still one scalar equation only:
  \[
  \boxed{\chi_Q=1.}
  \]
  On the natural source-map branch this is equivalent to
  \[
  N_Q=1,
  \qquad
  \Delta_Q:=\chi_Q-1=0,
  \qquad
  m_{\hat 0}^{\,2}\chi_QN_Q=1.
  \]
- The one-port same-charge bundle verdict is now sharp:
  the static mixed sector generates no new long-range attractive family, and the
  linear dynamic mixed sector generates no new kernel class beyond resonant
  dispersive enhancement of the already-short-range families.
- The actual weak-axisymmetric same-charge ceiling is now the packet
  \[
  (\Delta_{\rm norm},a_{P_0},b_{P_0}),
  \]
  equivalently
  \[
  (\Delta_{\rm norm},\Xi_1),
  \qquad
  \Xi_1=\frac{P_1}{P_0},
  \qquad
  b_{P_0}=3a_{P_0},
  \]
  with exact robust ceiling
  \[
  \boxed{
  \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}
  \bigl(1+|\epsilon\Xi_1|\bigr)
  \le P_{\rm crit}.
  }
  \]
- The first global kill condition on the carried same-charge branch is therefore
  still the transported static `\(\Xi_1\)` budget, not the wall-like dynamic
  window.
- On the rigid-mouth actual branch the physical packet is already diagonal in
  \[
  U=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
  \qquad
  V=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
  \]
  with exact dependent-plane compiler
  \[
  (\Delta_T,\Delta_{K_\eta},\Delta_\mu)=(0,-V,U-V),
  \]
  so orbit lock is the Cartesian condition
  \[
  U=V=0.
  \]
- The passive/outgoing support-selection problem has collapsed to the exact
  lowest-twin slice
  \[
  \rho_\alpha=\frac43,
  \qquad
  \zeta_{\rm req}=\frac13,
  \qquad
  \Pi_{\rm tr}=\frac43\,C_{\rm mix},
  \]
  equivalently the one-parameter selected twin-support curve
  \[
  \epsilon_* = 1-\frac{3\varrho}{2},
  \qquad
  \sigma = \frac{4}{3\varrho}-2,
  \qquad
  0<\varrho<\frac23.
  \]
- The selected-branch primitive hierarchy is exact in the three regions bounded
  by
  \[
  \varrho_{W\Lambda}
  =
  \frac{2(1+\beta^2)}{3(2+\beta^2)},
  \qquad
  \varrho_{U\Lambda}
  =
  \frac{2(1+\beta^2)}{3(1+\beta+\beta^2)},
  \]
  and it is generically an interference / overlap / wall-blocking /
  outgoing-scale correction rather than a large split-`U` effect.
- The smallest actual front-end branch packet is now
  \[
  \Bigl(
  \epsilon,\,
  \varrho_{\rm phys},\,
  \sigma_{\rm phys},\,
  \text{ranking region},\,
  R_{\rm tr},\,
  R_{\rm target},\,
  \epsilon_\eta,\,
  d\ln R_{\rm tr},\,
  d\ln R_{\rm target},\,
  d\ln\epsilon_\eta,\,
  N_Q-1
  \Bigr),
  \]
- The relaxed/open-system continuation now has an explicit codimension-three lift
  of the selected-branch front end,
  \[
  \mathfrak B_{226}^{\rm relax}
  =
  \{(\mathcal P_{225},\mathcal S_{225});L_w,L_{UV},L_\varsigma\},
  \]
  with exact recovery slice
  \[
  \ell_w=f_U=a=b=0.
  \]
- Its carried stationary lowered front end is
  \[
  V_{\rm eff}^{(230)}(r)
  =
  V_{\rm short}^{(1p)}(r)
  -
  \lambda_L S_{\rm leak}(r)
  -
  \lambda_W\mathcal W_w^{\rm sess}(r)
  -
  \Delta E_{UV}(r)
  -
  \mathcal M_{\sigma,\rm red}(r),
  \]
  and its dynamic continuation already supplies
  \[
  (r_{\rm peak},V_{\rm peak}),\quad
  v_{\rm crit,new},\quad
  r_\pm(E),\quad
  I_{\rm new}(E),\quad
  \Xi_{\rm turn}(E),\quad
  \lambda_{\rm th}(E).
  \]
- On the relaxed branch the hidden-channel unloading diagnostic is the one-scalar
  asymmetry problem
  \[
  \dot H_{{\rm sub},\sigma}=\Gamma_0(1+\sigma\alpha_h),
  \]
  while the active `V`-leg export law is the super-Ohmic kernel
  \[
  K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5,
  \]
  with exact event-safe surface
  \[
  \Gamma_3 s_c^3+\Gamma_5 s_c^5
  \ge
  \mu_\eta(s_0^2-s_c^2).
  \]
- The Stage-236 materials companion is now organized by the exact inequality stack
  \[
  \Pi_{\rm ep}\ge1,
  \qquad
  \Pi_\chi\ge1,
  \qquad
  \Pi_k\ge1,
  \qquad
  \Pi_T(T)\ge1,
  \]
  together with the calibration dictionary
  \[
  (t_*,\lambda_{\rm phys},E_*,\Upsilon_{\rm lat}).
  \]
  with separate support packet
  \[
  (\zeta,\ M_{\rm mix},\ S(\zeta;\epsilon),\ M_{\rm tr}).
  \]
- The full reduced back end is still the exact four-scalar packet
  \[
  \Delta_{\rm full}=(\Delta_Q,q_{\rm tr},q_{\rm nt},q_\eta),
  \]
  together with the equivalent multiplicative, mismatch, graph-error, and local
  mixed-ray charts. The local search sieve remains exact, but it is no longer
  the front-line bottleneck.

### 10.5 Practical global verdict

**Status:** `Open`

- None of the post-225 relaxed/open-system stages changes the earlier one-port
  same-charge verdict: there is still no new static long-range attractive family
  and no new linear dynamic conservative kernel class. The relaxed branch remains
  a short-range/open-system companion around the strict selected-support front end.
- There is currently no known algebraic blocker, no known stale-notation
  obstruction, no missing Packet-A compiler step, no remaining support-selection
  ambiguity on the passive/outgoing side, no hidden same-charge long-range kernel
  class left un-audited at the one-port level, and no unfinished
  support-cardinality theorem in the local search sieve.
- The strongest remaining theorem gap is now branch realization of one combined
  front-end packet, not another symbolic closure:
  whether the actual moving-throat branch returns
  \[
  N_Q=1,
  \qquad
  \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}
  \bigl(1+|\epsilon\Xi_1|\bigr)\le P_{\rm crit},
  \]
  together with the rigid-mouth orbit condition
  \[
  d\ln R_{\rm tr}=d\ln R_{\rm target}=d\ln\epsilon_\eta=0,
  \]
  on one and the same selected twin-support placement state.
- The free-quintuple graph/orbit packet and the local mixed-ray sieve are still
  the exact back-end repair and audit language, but the active front-end
  bottleneck is now: selected support placement, actual `\Xi_1` placement, and
  rigid-mouth orbit lock on the realized branch.

- If one follows the Stage-226-236 companion branch rather than staying on the
  strict Stage-225 slice, the active open objects are no longer symbolic closure
  formulas but realized coefficients and calibration data:
  \[
  \eta_{\rm leak},
  \quad
  f_U,\chi_{UV},
  \quad
  (a,b)\ \text{or}\ (\mathfrak g,\mathcal S,\mathcal R),
  \quad
  \Gamma_3,\Gamma_5,
  \]
  \[
  \Gamma_3^{\rm vac/lat},\ \Gamma_5^{\rm vac/lat},
  \qquad
  t_*,\ \lambda_{\rm phys},\ E_*,\ \Upsilon_{\rm lat},\ \mathcal K_{\rm corr}.
  \]
So the active bottleneck is no longer compression. It is realization of the
same combined support / orbit / outgoing packet on the true moving-throat branch.


## 11. Open Realization Gap

The remaining gap is now concentrated and explicit.

### 11.1 What is no longer the main risk

- not dropped factors, sign mistakes, or stale symbolic algebra,
- not disagreement about the grouped trace/anomaly projectors,
- not ambiguity in the Packet-A source/outgoing factorization,
- not a hidden higher-odd loophole in the point-particle `2.5`PN theorem,
- not an unfinished support-selection problem on the passive/outgoing side,
- not a missing static or linear dynamic same-charge kernel class at the one-port level,
- not uncertainty about whether the first global kill is the wall-like dynamic window,
- not ambiguity in the reference-free orbit packet or the same-free-quintuple
  repair map,
- not unfinished support-cardinality combinatorics in the local mixed-ray search.

### 11.2 What is still genuinely open

- If the relaxed/open-system companion is kept, the realized Stage-226-230
  branch data are also genuinely open:
  \[
  \eta_{\rm leak},
  \qquad
  f_U,\ \chi_{UV},
  \qquad
  (a,b)\ \text{or}\ (\mathfrak g,\mathcal S,\mathcal R),
  \]
  together with the resulting lowered front end
  \[
  V_{\rm eff}^{(230)}(r)
  \quad\text{and its dynamic tags}\quad
  (r_{\rm peak},V_{\rm peak},r_\pm,I_{\rm new},\Xi_{\rm turn},\lambda_{\rm th}).
  \]
- The realized microscopic export coefficients and their physical split remain open:
  \[
  \Gamma_3,\ \Gamma_5,
  \qquad
  \Gamma_3^{\rm vac/lat},\ \Gamma_5^{\rm vac/lat},
  \]
  as do the channel-resolved heat fractions and the exact cold-event shape quotient
  on the realized branch.
- The physical calibration dictionary of the material-screening companion is still open:
  \[
  t_*,\ \lambda_{\rm phys},\ E_*,\ \Upsilon_{\rm lat},\ \mathcal K_{\rm corr},
  \]
  and therefore so are the final host-screening ratios
  \[
  \Pi_{\rm ep},\ \Pi_\chi,\ \Pi_k,\ \Pi_T.
  \]
- The actual one-port reduced branch data, or any exactly equivalent physical
  packet, on the realized moving-throat branch:
  \[
  (\Delta,D_0,P_0,\Xi_1),
  \qquad
  \text{or equivalently }
  (\Delta_{\rm norm},\Xi_1).
  \]
- Whether the realized weak-axisymmetric branch satisfies the transported static
  ceiling
  \[
  \frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}
  \bigl(1+|\epsilon\Xi_1|\bigr)\le P_{\rm crit}.
  \]
- The actual rigid-mouth physical packet
  \[
  (R_{\rm tr},\mathcal T^2,\epsilon_\eta),
  \qquad
  \text{equivalently }(U,V),
  \]
  and whether it lands at
  \[
  U=V=0
  \quad\text{or}\quad
  d\ln R_{\rm tr}=d\ln R_{\rm target}=d\ln\epsilon_\eta=0
  \]
  on the same realized branch.
- The actual selected twin-support placement state
  \[
  (\epsilon,\varrho_{\rm phys},\sigma_{\rm phys},\text{ranking region}),
  \]
  together with the separate support packet
  \[
  (\zeta,M_{\rm mix},S(\zeta;\epsilon),M_{\rm tr}).
  \]
- Whether one and the same realized branch satisfies
  \[
  N_Q=1,
  \qquad
  \Xi_1\text{-ceiling admissibility},
  \qquad
  \text{and rigid-mouth orbit lock},
  \]
  rather than satisfying them only on nearby or repaired branches.
- The actual directional-Hessian / admissibility packets needed to populate the
  now-finite local mixed-ray search ledger if the back-end graph repair audit is
  still needed after the front-end branch test.

### 11.3 How to read the current program honestly

- The document now supports the claim that there is a mathematically coherent
  reduced derivation chain from the parent PDE setup to an exact Packet-A finish
  line, an exact Packet-B graph/orbit compiler, and a finite local search sieve.
- It does not yet justify collapsing every reduced closure into an unconditional
  full-PDE realization theorem.
- So the right referee-facing language is: exact parent theory, exact statements
  inside explicit closures, numerically located reduced branches, exact verdict /
  graph / projector compilers, a finite local search ledger, and an explicitly
  isolated realization gap.


## 11A. V2 Branch-Extraction Working Protocol

This section is the practical AI protocol for using the V2-full compact when a
candidate branch, numerical output, or analytic reduced profile is supplied.

### 11A.1 Branch-freeze checklist

Before computing any target residuals, require a manifest containing:

```text
parent_action_status        # strict current action or promoted S_Sigma
gauge_fixing_convention     # H=1 pre-Lorenz or localized H=Z
boundary_class              # must be open_impedance for V2 branch tests
R_exit                      # must be > 0
projection_source_map       # mhat0 and S_port convention
support_family              # D/N AC support, Robin/impedance, etc.
mode_port_list              # wall, BdG, U/W mixed ports
stability_gates             # wall positivity, BdG positivity, Delta>0, D0>0
extraction_formulas         # definitions of K,M,B,Z,N,P,u
pre_target_freeze=true
target_blind=true
no_post_residual_refit=true
```

If any of these are missing, the branch is not ready for target-realization
judgment.

### 11A.2 Observable extraction steps

1. Validate open geometry and profile arrays.
2. Normalize or explicitly record the axial/radial profile measures.
3. Extract wall coefficients \(K,M\).
4. Extract BdG moments \(B_0,B_2,B_4\) from positive stable modes.
5. Extract Maxwell/mixed conservative moments \(Z_0,Z_2,Z_4\).
6. Extract outgoing-transfer moments \(N_0,N_2,N_4\).
7. Compute
   \[
   D_0=K-B_0-Z_0,
   \quad
   D_2=-(M+B_2+Z_2),
   \quad
   D_4=-(B_4+Z_4).
   \]
8. Compute \(u_2,u_4,P_0,P_2,P_4\) and grouped trace/anomaly data.
9. Evaluate residuals \(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4}\), plus
   \(R_{\rm tail}\) when active.
10. Report pass/fail without changing branch data.

### 11A.3 Common wrong moves to avoid

- Do not treat \(\Sigma/R\) as a strict parent dynamical field unless
  \(S_\Sigma\) has been included.
- Do not integrate the unweighted five-dimensional gauge-fixing term over the
  noncompact zero mode.
- Do not use a hard cap \(R(L)=0\) for branch-realization tests.
- Do not erase \(A_w,J^w,F_{\mu w},E_w,C_a\) because the far-field brane limit
  suppresses them.
- Do not tune \(K,N_0,N_2,N_4,\Theta_{\rm tail}\) after seeing target residuals.
- Do not call a stable target-failing branch a pipeline failure.
- Do not treat grouped-lane anisotropy as arbitrary if the perturbation is weak
  axisymmetric: it must obey \(b=3a\).

### 11A.4 Minimal branch-output schema

```yaml
branch_id: string
branch_freeze_hash: string
pre_target_freeze: true
target_blind: true
no_post_residual_refit: true
boundary_class: open_impedance
R_mouth: number
R_exit: number    # > 0
wall:
  K: number
  M: number
  positivity_pass: bool
bdg:
  modes:
    - coupling: number
      varpi: number   # > 0
maxwell_mixed:
  ports:
    - Omega_U: number
      Omega_W: number
      R: number
      g_U: number
      g_W: number
      Delta: number   # > 0
source_map:
  mhat0: number
  S_port: number
constants:
  G: number
  c: number
  c_s: number
  a: number
optional_tail:
  Theta_tail: number
```

The compact formulas above are enough to turn this schema into the V2 residual
packet.


## 12. Quick Translation Dictionary

### 12.1 Core fields and brane reduction hooks

The exact parent fields are
\[
\psi(\mathbf X,t),
\qquad
\rho=|\psi|^2,
\qquad
A_M=(A_0,A_i),
\qquad
F_{MN}=\partial_MA_N-\partial_NA_M.
\]

The exact brane projection hooks are
\[
\rho_{\rm brane}=\int W(w)\rho(\mathbf x,w,t)\,dw,
\qquad
\mathbf j_{\rm brane}=\int W(w)\mathbf j_{xyz}(\mathbf x,w,t)\,dw,
\]
\[
\mathbf v_{\rm brane}=\mathbf j_{\rm brane}/\rho_{\rm brane},
\qquad
\partial_t\rho_{\rm brane}+\nabla_3\cdot\mathbf j_{\rm brane}=S_{\rm leak}.
\]

In the quasi-static longitudinal regime,
\[
\mathbf v_{\rm brane}=\nabla_3\varphi+\mathbf v_T,
\qquad
\nabla_3\cdot\mathbf v_T=0,
\]
and \(\varphi\) is the brane velocity-potential / Poisson-hook variable.

### 12.2 Moving-throat geometry variables

The moving throat is represented by the level-set / shape field
\[
\Sigma(\mathbf X,t)=r-R(\Omega,w,t),
\qquad
r=\sqrt{x^2+y^2+z^2},
\qquad
\Omega=\mathbf x/r.
\]

The throat surface is \(\Sigma=0\), and the wall displacement is
\[
R(\Omega,w,t)=R_0(w)+\eta(\Omega,w,t).
\]

The old collective variables
\[
a(t),\qquad L(t)
\]
are collective moments of the moving-throat field, not the fundamental geometry
variables.

### 12.3 Harmonic and grouped real `P2` variables

The mouth-sphere harmonic decomposition is
\[
\eta(\Omega,w,t)
=
\eta_0(w,t)Y_{00}(\Omega)
+\sum_{m\in P_2({\rm real})}q_{2m}(w,t)Y_{2m}^{\rm real}(\Omega)+\eta_{\ge 3}.
\]

The grouped real `P2` lanes are the weighted triplet
\[
A\in\{20,21,22\},
\]
with the full five-mode real bundle understood underneath.

The grouped conservative response variables are
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\omega^2+D_{A,4}\omega^4+O(\omega^6),
\]
\[
Y_A(\omega)=\frac{D_{A,0}}{D_A^{\rm(cons)}(\omega)}
=
1+u_2^{(A)}\omega^2+u_4^{(A)}\omega^4+O(\omega^6).
\]

The outgoing prefactor data are
\[
P_0=\frac{N_0}{D_0},
\qquad
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\qquad
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

### 12.4 Family-1 mouth/core variables

The throat-core closure variables are
\[
K_s,\qquad K_q,\qquad \lambda,\qquad g_s,\qquad g_q.
\]

Their normalized Family-1 ratios are
\[
\mathfrak r=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}}.
\]

The mouth/core control variables used later are
\[
\Sigma_0=\frac{20}{9}\widehat T_m^2,
\qquad
\mathcal R=\frac{(\mathfrak g-\mathfrak r)^2}{1+\mathfrak r^2},
\]
\[
\mathcal S[\Sigma]=\int_0^1 K_q(x)\Sigma(x)\,dx,
\qquad
\Pi=\Sigma_0\,[1-\mathcal R\,\mathcal S].
\]

So:

- \(\widehat T_m\) is the normalized mouth traction,
- \(\Sigma_0\) is its exact reduced quadratic proxy,
- \(\mathcal R\) is the mixed-to-shell loading ratio,
- \(\mathcal S\) is the mouth-profile quadrupole overlap,
- \(\Pi\) is the mouth slope / bias variable.

### 12.5 Coherent local-kernel variables

The coherent reduced hierarchy is organized by
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\]
\[
Z_W=\frac{\lambda_W^2}{K_\eta^{(\mathrm{eff})}K_W^{(\mathrm{eff})}},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\]
\[
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\zeta=\frac{\lambda_\phi^2K_W^{(\mathrm{eff})}}{\lambda_W^2K_\phi^{(\mathrm{eff})}}.
\]

The key coherent tracking factor is
\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
and the coherent support-enhancement factor is
\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\qquad
\epsilon=\epsilon_W\!\left[1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right].
\]

### 12.6 Microscopic grouped defect variables

The positive coherent microscopic state is
\[
x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]

Its grouped weak-axisymmetric log-drift vector is
\[
\delta\mathbf x
=
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1)^T.
\]

The direct microscopic slippages are
\[
\Sigma_Z=\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right),
\qquad
\Sigma_\chi=\delta\ln\chi_0,
\qquad
\Sigma_\eta=\delta\ln\epsilon_\eta,
\]
\[
\Sigma_\epsilon=\delta\ln\epsilon_W,
\qquad
\Sigma_\delta=\delta\ln\delta_U.
\]

The branch-adapted coordinates are
\[
\Sigma_{\rm tr}=(1+\chi_0)\Sigma_\delta+(1+\delta_U)\Sigma_\chi,
\]
\[
\Sigma_{\rm nt}
=
\Sigma_Z
+
\frac{2\epsilon_W}{1-\epsilon}\frac{11+9\delta_U}{11(1+\delta_U)}\,\Sigma_\epsilon
-\left[
\frac{2\chi_0}{1+\delta_U}
+
\frac{4\epsilon_W\delta_U}{11(1-\epsilon)(1+\delta_U)^2}
\right]\Sigma_\delta.
\]

The observable defect variables are

- \(\Theta_1\): tracking-factor drift
- \(\Xi_1\): grouped transfer-shape drift
- \(\mathcal R_1\): selected-branch demand-ratio drift

and they satisfy the exact triangular normal form
\[
\Theta_1 \leftrightarrow \Sigma_{\rm tr},
\qquad
\Xi_1 \leftrightarrow \Sigma_{\rm nt},
\qquad
\mathcal R_1+\Xi_1 \leftrightarrow \Sigma_\eta.
\]

### 12.7 Final invariant variables

The final reduced coherent invariant coordinates are
\[
\mathfrak C_{{\rm tr},*}
=
\chi_0^{\,1+\delta_{U,*}}\,
\delta_U^{\,1+\chi_{0,*}},
\]
\[
\mathfrak C_{{\rm nt},*}
=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*},
\]
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
\]

Inside the reduced coherent hierarchy, these are the three quotient coordinates
carried by the final invariant theorem.


### 12.8 Branch-observable, transfer-shape, and rigid-mouth physical packet

The Stage-171 PDE-facing observable packet is
\[
\Delta_{\rm obs}^{(1)}
=
(\delta\ln R_{\rm tr},\ \delta\ln\mathfrak N_*,\ \delta\ln\epsilon_\eta)^T,
\qquad
\mathfrak N_*:=\mathcal T^2R_{\rm tr}^{B_*}.
\]

Here
\[
B_*=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}},
\qquad
A_{{\rm tr},*}=B_*C_{{\rm tr},*}.
\]

The exact transfer-shape / selected-branch variables are
\[
\mathcal T^2=\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\qquad
R_{\rm target},
\qquad
1-\epsilon_\eta,
\]
with
\[
R_{\rm target}\mathcal T^2
=
\Lambda_0(1-\epsilon_\eta),
\qquad
\Lambda_0=\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]

At first grouped weak-axisymmetric order,
\[
\delta\ln\mathcal T^2=\Xi_1,
\qquad
\delta\ln(1-\epsilon_\eta)=\mathcal R_1+\Xi_1,
\qquad
\delta\ln R_{\rm target}=\mathcal R_1.
\]

On the rigid-mouth slice `\(q_{\rm tr}=0\)`, the exact physical packet is already
diagonal in
\[
\boxed{
U:=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V:=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
}
\]
with
\[
\boxed{
q_{\rm nt}=U,
\qquad
q_\eta=V.
}
\]
So the actual rigid-mouth branch is charted directly by the physical pair
`(U,V)` rather than by a triangular quotient representation.

### 12.9 Packet-A finish-line and actual prefactor-ceiling variables

Packet A is no longer best thought of as an eight-slot branch residual at theorem
order. Its exact conservative front end is the isotropic one-pole grouped-real
`P2` surface
\[
a_2=b_2=a_4=b_4=0,
\qquad
\Delta_{\rm pole}=\bar u_4-4\bar u_2^{\,2}=0,
\]
with common conservative carrier
\[
\widehat Y_Q^{\rm cons}(\omega)=\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
\]
The exact outgoing scalar is
\[
\chi_Q=\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0},
\]
and on the natural point-particle source-map branch one has
\[
N_Q=\chi_Q^{-1},
\qquad
m_{\hat 0}^{\,2}\chi_QN_Q=1.
\]
So the Packet-A finish line is the single scalar
\[
\Delta_Q:=\chi_Q-1,
\]
equivalently
\[
\Delta_{\rm norm}=P_0^{\rm target}(\chi_Q^{-1}-1),
\qquad
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
\]

On the actual weak-axisymmetric branch, the first same-charge ceiling is carried
by the prefactor packet
\[
(\Delta_{\rm norm},a_{P_0},b_{P_0}),
\]
or equivalently by
\[
\boxed{
(\Delta_{\rm norm},\Xi_1),
\qquad
\Xi_1=\frac{P_1}{P_0},
}
\]
with
\[
\boxed{
a_{P_0}=\frac{\epsilon\bar P_0\Xi_1}{4},
\qquad
b_{P_0}=\frac{3\epsilon\bar P_0\Xi_1}{4}.
}
\]
The robust all-lane ceiling is
\[
\boxed{
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}
\bigl(1+|\epsilon\Xi_1|\bigr)\le P_{\rm crit}.
}
\]
So Packet A now has two front-end faces:

- the isotropic outgoing finish line `\(\chi_Q=1\)`,
- the transported actual-branch prefactor ceiling in `\((\Delta_{\rm norm},\Xi_1)\)`.

### 12.10 Packet-B orbit / graph variables and actual placement packet

Important notation firewall:
\[
R_{\rm tr}
\]
is the coherent **tracking factor**, while
\[
\mathfrak R_{\rm tr}
\]
is the finite **orbit/invariant ratio** used in Packet B.

The additive Packet-B orbit packet is
\[
(q_{\rm tr},q_{\rm nt},q_\eta).
\]
The multiplicative orbit packet is
\[
(\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta),
\qquad
\mathfrak R_{\rm tr}=e^{q_{\rm tr}},
\quad
\mathfrak R_{\rm nt}=e^{q_{\rm nt}},
\quad
\mathfrak R_\eta=e^{q_\eta}.
\]
The mismatch packet is
\[
(m_T,m_K,m_\mu),
\]
with exact conversions
\[
m_T=e^{q_{\rm tr}/(1+\chi_{0,*})},
\qquad
m_K=e^{-q_\eta},
\qquad
m_\mu=e^{q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}}.
\]
On the same-free-quintuple graph, the dependent-triple graph errors are
\[
E_T=\frac{q_{\rm tr}}{1+\chi_{0,*}},
\qquad
E_K=-q_\eta,
\qquad
E_\mu=q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}.
\]
So the graph-error packet
\[
(E_T,E_K,E_\mu)
\]
is just the Packet-B mismatch packet written directly on the dependent triple.

On the actual coherent placement branch, the finite packet is also available
directly in the physical placement variables:
\[
\boxed{
q_{\rm tr}
=
(1+\delta_{U,*})\ln\frac{\chi_0}{\chi_{0,\rm ref}}
+
(1+\chi_{0,*})\ln\frac{\delta_U}{\delta_{U,\rm ref}},
}
\]
\[
\boxed{
q_{\rm nt}
=
\ln\frac{Z_W}{Z_{W,\rm ref}}
-
\ln\frac{\Lambda}{\Lambda_{\rm ref}}
+
E_*\ln\frac{\epsilon_W}{\epsilon_{W,\rm ref}}
-
F_*\ln\frac{\delta_U}{\delta_{U,\rm ref}},
}
\]
\[
\boxed{
q_\eta=\ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
}
\]
The coherent orbit-lock theorem is
\[
\boxed{
q_{\rm tr}=q_{\rm nt}=q_\eta=0
\iff
d\ln R_{\rm tr}=d\ln R_{\rm target}=d\ln\epsilon_\eta=0.
}
\]

The actual placement / support front-end packet is
\[
\boxed{
\Bigl(
\epsilon,\,
\varrho_{\rm phys},\,
\sigma_{\rm phys},\,
\text{ranking region},\,
R_{\rm tr},\,
R_{\rm target},\,
\epsilon_\eta,\,
d\ln R_{\rm tr},\,
d\ln R_{\rm target},\,
d\ln\epsilon_\eta,\,
N_Q-1
\Bigr),
}
\]
with separate support packet
\[
\boxed{
(\zeta,\ M_{\rm mix},\ S(\zeta;\epsilon),\ M_{\rm tr}).
}
\]
So Packet B still has the exact graph/orbit backend, but it now also has a
concrete actual-branch front-end placement chart.

### 12.11 Free-quintuple graph and local search variables

The free quintuple is
\[
\mathbf y=(\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\mathrm{eff})}).
\]
The target orbit is the exact graph
\[
\mathcal O_*=
\{\mathbf x_*^{\rm graph}(\mathbf y):\mathbf y\in(\mathbb R_{>0})^5\},
\]
and the graph-closure scalar is
\[
\widehat\chi_Q(\mathbf y)=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y)).
\]
So the full reduced closure set is the codimension-one graph slice
\[
\mathcal Z_*=
\{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\chi_Q(\mathbf y)=1\}.
\]
The explicit one-parameter free-family compiler is the log-ray
\[
\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s},
\]
with scalarized ray function
\[
\Phi_{\mathbf s}(\tau)=\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau)).
\]
Its local directional packet is
\[
\Phi_0,
\qquad
L_0=\frac{d}{d\tau}\ln\Phi_{\mathbf s}\Big|_0,
\qquad
L_1=\frac{d^2}{d\tau^2}\ln\Phi_{\mathbf s}\Big|_0.
\]
The finite local mixed-ray search is organized by support-cardinality on the
five primitive axes `\(\{\lambda,c,\gamma,U,W\}\)`, with exact ceiling `5` and final
splice
\[
\tau_{\le 5,*}^{\rm best}=
\min\bigl(\tau_{\le 4,*}^{\rm best},\tau_{5,*}^{\rm best,int}\bigr).
\]
The preferred full evaluation budget is `1464`, with fallback `2640`.

### 12.12 Relaxed/open-system barrier companion variables

The post-225 relaxed branch is organized by three lifted lanes:
\[
L_w=(S_{\rm leak},\mathcal W_w),
\qquad
L_{UV}=(U,V,\mathcal D_{UV}),
\qquad
L_\varsigma=(\varsigma,\varsigma_{\min},\mathrm{signchg}).
\]

The selected-support leakage/work compiler uses
\[
\Pi_{\rm tr},
\qquad
\eta_{\rm leak},
\qquad
S_{\rm leak},
\qquad
\mathcal W_w^{\rm bulk},
\qquad
\mathcal W_w^{\rm sess}.
\]

The orbit-side non-rigid packet uses
\[
a_U,\qquad
a_V,\qquad
\chi_{UV},\qquad
\Delta_{UV}=a_Ua_V-\chi_{UV}^2,
\]
\[
U=\ln\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2},
\qquad
V=\ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}},
\qquad
\mathcal D_{UV}=\chi_{UV}UV.
\]

The compensated source packet uses
\[
\varsigma(x;r)=1+a(r)\cos(\pi x)+b(r)\cos(2\pi x),
\]
\[
\mathfrak g[\varsigma],\qquad
\mathcal S[\varsigma],\qquad
\mathcal R[\varsigma],
\qquad
\mathcal M_\sigma^{\rm pack}(r).
\]

The carried lowered stationary front end is
\[
V_{\rm eff}^{(230)}(r)
=
V_{\rm short}^{(1p)}(r)
-\lambda_L S_{\rm leak}(r)
-\lambda_W\mathcal W_w^{\rm sess}(r)
-\Delta E_{UV}(r)
-\mathcal M_{\sigma,\rm red}(r).
\]

Its dynamic event-chain variables are
\[
r_{\rm peak},\qquad
V_{\rm peak},\qquad
v_{\rm crit,new},\qquad
r_\pm(E),\qquad
I_{\rm new}(E),\qquad
\Xi_{\rm turn}(E),\qquad
\lambda_{\rm th}(E).
\]

The wall-timescale / unloading / export variables are
\[
t_{\rm cross},\qquad
t_{\rm collapse},\qquad
\mathcal S(E),\qquad
E_{\rm edge},\qquad
v_{\rm safe,min},
\]
\[
\alpha_h,\qquad
R_{\rm inst},\qquad
R_{\rm int},
\qquad
\Gamma_3,\qquad
\Gamma_5,
\qquad
K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5.
\]

The heat-partition / materials-companion variables are
\[
\Gamma_3^{\rm vac/lat},\qquad
\Gamma_5^{\rm vac/lat},\qquad
\mathfrak r_V,\qquad
f_{\rm vac},\qquad
f_{\rm lat},
\]
\[
\gamma_{\rm vac}^{\rm eq}(s),\qquad
\gamma_{\rm lat}^{\rm eq}(s),\qquad
\gamma_{\rm lat,safe}^{\rm eq},
\]
\[
t_*,\qquad
\lambda_{\rm phys},\qquad
E_*,\qquad
\Upsilon_{\rm lat},\qquad
\mathcal K_{\rm turn},\qquad
\mathcal K_{\rm corr},
\]
\[
\Pi_{\rm ep},\qquad
\Pi_\chi,\qquad
\Pi_k,\qquad
\Pi_T(T).
\]


---

## 13. Ultra-Compact Recovery Cards

These cards are intentionally redundant. They are the fastest way for a new AI
session to recover the working state without searching the document.

### 13.1 Parent / gauge / boundary card

```text
Strict parent: S_current = S_psi[psi,A;Sigma] + S_EM[A]
Wall PDE strict parent status: requires + S_Sigma[Sigma]
Quadratic wall closure: S_eta^(2), effective unless promoted
Gauge H=1: impose Lorenz before zero-mode reduction
Gauge H=Z: finite localized gauge fixing, xi_4 = xi
Physical branch endpoint: R(L) > 0
D/N ladder: AC impedance/reflection, not a hard cap
Mixed channels: A_w, J^w, F_{mu w}, E_w, C_a stay alive
```

### 13.2 Full-bundle coefficient card

\[
D_0=K-B_0-Z_0,
\quad
D_2=-(M+B_2+Z_2),
\quad
D_4=-(B_4+Z_4).
\]
\[
u_2=-D_2/D_0,
\qquad
u_4=(D_2^2-D_0D_4)/D_0^2.
\]
\[
P_0=N_0/D_0,
\quad
P_2=(D_0N_2-2D_2N_0)/D_0^2.
\]
\[
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

### 13.3 V2 target card

\[
R_{\rm pole}=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2=0.
\]
\[
R_{\rm norm}=\widehat m_0^{\,2}\mathcal S_{\rm port}N_0/D_0
-54Gc_s^5/(5a^5c^5)=0.
\]
\[
P_2=P_4=0.
\]
\[
R_{\rm tail}=\Theta_{\rm tail}(c/c_s)^3-1=0\quad\text{if active}.
\]

### 13.4 5PN actual-branch card

\[
d\ln R_{\rm tr}=d\ln R_{\rm target}=d\ln\epsilon_\eta=0,
\qquad
N_Q=1.
\]
\[
\Xi_1=P_1/P_0,
\qquad
b=3a\quad\text{for weak axisymmetric grouped anisotropy}.
\]

### 13.5 Solver card

```text
freeze -> validate -> profile adapter -> observable extraction -> residual packet
mechanical pass does not imply target realization
first V2-23 reduced open branch: stable/open but target-failing
```
