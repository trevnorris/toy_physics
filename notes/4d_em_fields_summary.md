4d_em_fields.tex

## 1) What this paper is trying to accomplish

This paper ("Paper VIII" in the toy-model sequence) isolates the **electromagnetic / gauge-field sector** of the unified brane-bulk toy model and shows how ordinary **3+1 Maxwell theory on the brane** arises from a **localized 4+1 Maxwell action**.

It answers the referee-style question left open by the earlier EM-from-flow papers:

> Where do the **inhomogeneous** Maxwell equations, propagators, and retarded structure come from dynamically, without being postulated by hand?

### 1.1 Charge ontology update

This paper now uses the same charge ontology as the corrected `4d.tex`.

1. Electric-charge sign is **not** defined by circulation, throat radius, or breathing variables.
2. The fixed topological puncture orientation is
   \[
   \eta_Q=\pm1.
   \]
3. A positive microscopic charge scale `e_*` and signed defect-branch coupling are introduced by
   \[
   e_*>0,
   \qquad
   q_*\equiv \eta_Q e_*.
   \]
4. The observable brane charge strength is controlled by localization thickness. After canonical zero-mode normalization,
   \[
   e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}},
   \qquad
   q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}}=\eta_Q e_{\rm eff}.
   \]
5. Circulation belongs to the quantized magnetic/vortical sector, not to the electric-charge dictionary.
6. The zero-mode Maxwell reduction suppresses the mixed core sector \((A_w,J^w,F_{\mu w})\) only as a controlled far-field brane limit. It does not remove those channels from the microscopic ontology.

### Core deliverables (as stated/derived)

1. **Bulk dynamics (exact):** starting from a localized 4+1 action, derive the bulk Euler-Lagrange equations
   \[
   \partial_M\big(Z(w)F^{MN}\big)+\frac{1}{\xi}\,\partial^N(\partial\!\cdot\!A)=\mu_0 J^N,
   \]
   together with the Bianchi identities and the divergence consistency identity tying gauge fixing to current conservation.

2. **Controlled brane reduction (Maxwell limit):** under explicit reduction assumptions (axial gauge, zero-mode dominance, brane-localized sources), integrate over the transverse direction and obtain standard 3+1 Maxwell on the brane with computable effective coupling
   \[
   \mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}},
   \qquad
   Z_{\rm int}=\int_{-\infty}^{\infty} Z(w)\,dw.
   \]
   For a Gaussian profile \(Z(w)=e^{-w^2/\lambda^2}\), this becomes
   \[
   \mu_0^{\rm eff}=\frac{\mu_0}{\lambda\sqrt\pi}.
   \]

3. **Thickness-controlled brane charge strength:** the same reduced action also implies the equivalent canonical-normalization result
   \[
   q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
   \qquad
   e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}}.
   \]
   Thus the microscopic branch label \(q_*\) is fixed, while the observable brane coupling is weakened by thicker localization.

4. **KK corrections (quantitative, falsifiable):** for Gaussian localization, solve the transverse Sturm-Liouville problem in closed form. The result is a discrete KK tower with a rigid mass/coupling pattern that yields Yukawa-suppressed departures from Coulomb's law. For a fixed defect branch,
   \[
   A_0(r)=\frac{\mu_0^{\rm eff}q_*}{4\pi r}
   \Big[1+\tfrac12 e^{-2r/\lambda}+\cdots\Big].
   \]
   Equivalently, after canonical normalization the same brane field may be written in terms of \(q_{\rm eff}\).

5. **Moving sources (causal + brane-Lorentz covariant):** each KK mode propagates with a covariant retarded Green function. Massive modes produce causal tails inside the light cone; summing the tower preserves causality and depends on brane momenta only through the Lorentz scalar \(k^2\).

6. **Source consistency (derived, not imposed):** a minimal brane matter model with local \(U(1)\) phase symmetry produces a conserved Noether current, which supplies the conserved source required by the Maxwell divergence identity. In the updated bookkeeping,
   \[
   J_\psi^\mu=q_* j^\mu,
   \qquad
   J_{\rm tot}^M=J_\psi^M+J_{\rm ext}^M.
   \]


## 2) Frozen conventions, objects, and notation

### 2.1 Coordinates, indices, and metric

Bulk coordinates:
\[
x^M=(t,x,y,z,w),
\qquad
M,N\in\{0,1,2,3,4\}.
\]

Brane indices:
\[
\mu,\nu\in\{0,1,2,3\},
\qquad
x^\mu=(t,x,y,z).
\]

Brane spatial indices:
\[
i,j\in\{1,2,3\}.
\]

Flat bulk metric (mostly plus):
\[
\eta_{MN}=\mathrm{diag}(-1,1,1,1,1),
\qquad
\eta^{MN}=(\eta_{MN})^{-1}.
\]

Bulk d'Alembertian:
\[
\Box\equiv \eta^{MN}\partial_M\partial_N
= -\partial_t^2+\nabla^2+\partial_w^2.
\]

### 2.2 Localization profile

A nonnegative localization profile \(Z(w)\) multiplies the Maxwell kinetic term. This paper chooses a Gaussian:
\[
Z(w)=\exp\!\left(-\frac{w^2}{\lambda^2}\right),
\]
with finite integral
\[
Z_{\rm int}=\int_{-\infty}^{\infty} Z(w)\,dw = \lambda\sqrt\pi.
\]

All integrations by parts in \(w\) assume the relevant boundary terms vanish as \(|w|\to\infty\).

### 2.3 Charge labels

The paper uses the fixed topological sign label
\[
\eta_Q\in\{+1,-1\},
\qquad
q_*=\eta_Q e_*,
\qquad
e_*>0.
\]

The canonically normalized brane coupling is
\[
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
\qquad
e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}}.
\]

Interpretation:

* \(q_*\) is the fixed microscopic defect-branch charge label.
* \(q_{\rm eff}\) is the observable brane coupling after zero-mode normalization.
* Charge sign is carried by \(\eta_Q\), not by circulation or throat geometry.

### 2.4 Gauge field and field strength

Bulk gauge potential \(A_M(x)\) is treated as a **covector**. Field strength:
\[
F_{MN}=\partial_M A_N-\partial_N A_M,
\qquad
F^{MN}=\eta^{MP}\eta^{NQ}F_{PQ}.
\]

Divergence:
\[
\partial\!\cdot\!A\equiv \partial_M A^M,
\qquad
A^M\equiv \eta^{MN}A_N.
\]

### 2.5 Gauge fixing parameter

A Lorenz-type gauge-fixing term is included with parameter \(\xi\) (called \(\xi_{\rm gf}\) in the paper). Gauge-invariant observables do not depend on \(\xi\) when the total source is conserved.

### 2.6 Source/current bookkeeping

The compact standalone Maxwell-sector notation uses a contravariant bulk current \(J^M(x)\) with source coupling
\[
\mathcal L_{\rm int}=-A_M J^M.
\]

In the fully coupled embedding, however, the split is
\[
J_{\rm tot}^M=J_\psi^M+J_{\rm ext}^M,
\]
where:

* \(J_\psi^M\) is the dynamical matter current generated by minimal coupling;
* \(J_{\rm ext}^M\) contains external/background sources, including neutralizing jellium if used;
* the explicit Maxwell action term is most cleanly read as \(-A_M J_{\rm ext}^M\) once matter is coupled separately.

So the bare symbol \(J^M\) in the standalone Maxwell derivation should be read as a shorthand total current, and one must not double-count \(J_\psi^M\).


## 3) Full localized 4+1 Maxwell sector (the action)

The localized Maxwell action is
\[
S_{\rm EM}[A;Z]
=\int dt\,d^3x\,dw\;
\left[
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
-A_M J^M
\right].
\]

Notes:

* Only the Maxwell kinetic term is multiplied by \(Z(w)\).
* In the coupled reading, the explicit source term should be interpreted with the bookkeeping of §2.6.
* The reduction result \(\mu_0^{\rm eff}=\mu_0/Z_{\rm int}\) remains primary, but the equivalent charge-normalization statement is also recorded explicitly in this updated version.


## 4) Exact Euler-Lagrange equations and identities (bulk statements)

Everything in this section is **exact**.

### 4.1 Bulk Maxwell equations (with localization + gauge fixing)

Varying \(A_N\) gives
\[
\partial_M\big(Z(w)F^{MN}\big)
+\frac{1}{\xi}\,\partial^N(\partial\!\cdot\!A)
=\mu_0 J^N.
\]
Define the Maxwell operator
\[
\mathcal M^N[A]\equiv
\partial_M\big(ZF^{MN}\big)
+\frac{1}{\xi}\,\partial^N(\partial\!\cdot\!A),
\qquad
\mathcal M^N[A]=\mu_0 J^N.
\]

### 4.2 Bianchi identities

From the definition of \(F\),
\[
\partial_{[L}F_{MN]}=0
\quad\Leftrightarrow\quad
\partial_LF_{MN}+\partial_MF_{NL}+\partial_NF_{LM}=0.
\]

### 4.3 Divergence consistency identity

Taking \(\partial_N\) of the EOM,
\[
\partial_N\partial_M(ZF^{MN})=0,
\]
so one obtains
\[
\partial_N(\mathcal M^N[A]-\mu_0 J^N)
=
\frac{1}{\xi}\,\Box(\partial\!\cdot\!A)-\mu_0\,\partial_N J^N.
\]

Consequences:

* In Lorenz gauge \(\partial\!\cdot\!A=0\), the Maxwell equations enforce bulk current conservation.
* Conversely, conserved sources are exactly what make the gauge-fixed description consistent.


## 5) Gauge structure and brane Lorentz covariance

### 5.1 Gauge transformations

Infinitesimal gauge transformation:
\[
\delta_\chi A_M=\partial_M\chi.
\]

Then \(\delta_\chi F_{MN}=0\), while the interaction term varies as
\[
\delta_\chi \mathcal L_{\rm int}
=-(\partial_M\chi)J^M
= -\partial_M(\chi J^M)+\chi\,\partial_M J^M.
\]
So gauge invariance of the coupling is equivalent to \(\partial_M J^M=0\).

### 5.2 Brane Lorentz symmetry + transformation types

Because \(Z(w)\) depends only on \(w\), the theory preserves exact **3+1 Poincare symmetry** along the brane.

Transformation types:

* \(A_\mu\) transforms as a **covector**.
* \(J^\mu\) transforms as a **vector**.

Under a brane Lorentz transform \(x'^\mu=\Lambda^\mu{}_{\nu}x^\nu\):
\[
A'_\mu(x')=(\Lambda^{-1})^\nu{}_{\mu}A_\nu(x),
\qquad
J'^\mu(x')=\Lambda^\mu{}_{\nu}J^\nu(x).
\]

### 5.3 Propagator depends only on \(k^2\)

Fourier transforming along the brane, \(\partial_\mu\to i k_\mu\), the response kernel becomes
\[
D_{\rm eff}(k^2)\sim\sum_n \frac{c_n}{k^2+m_n^2+i\epsilon},
\qquad
k^2\equiv \eta^{\mu\nu}k_\mu k_\nu=-\omega^2+\mathbf k^2.
\]
This is the momentum-space form matched by the reduced mode operator
\((\square_4+m_n^2)\).

This dependence on brane momentum only through \(k^2\) makes the reduced response manifestly brane-Lorentz covariant.


## 6) Controlled brane reduction and the effective coupling \(\mu_0^{\rm eff}\)

Everything here is **controlled** and requires explicit reduction assumptions.

### 6.1 Axial gauge reachability and residual 3+1 gauge freedom

Under \(\delta A_M=\partial_M\chi\), the transverse component shifts as
\[
A_w\mapsto A_w+\partial_w\chi.
\]
One may reach axial gauge by choosing
\[
\chi_{\rm Ax}(x^\mu,w)\equiv -\int_0^w A_w(x^\mu,w')\,dw',
\]
so that \(\partial_w\chi_{\rm Ax}=-A_w\).

Residual gauge freedom in axial gauge satisfies \(\partial_w\chi=0\), hence
\[
\chi(x^\mu,w)=\chi_b(x^\mu),
\]
which is ordinary 3+1 gauge freedom on the brane.

### 6.2 Brane-localized sources

A conserved brane source is modeled in the controlled far-field Maxwell reduction as
\[
J^\mu(x,w)=J_b^\mu(x)\,\delta(w),
\qquad
J^w(x,w)=0.
\]
Then \(\partial_M J^M=0\) reduces distributionally to \(\partial_\mu J_b^\mu=0\). Here
\(J^w=0\) is the reduction ansatz used to isolate the brane Maxwell sector, not a
microscopic claim that the charged core lacks mixed \((w)\)-structure.

In the updated ontology, the brane source may itself be split as
\[
J_{b,\rm tot}^\mu=J_{\psi,b}^\mu+J_{{\rm ext},b}^\mu,
\qquad
J_{\psi,b}^\mu=q_* j_b^\mu.
\]
So the sign of the dynamical electric branch is carried by \(q_* = \eta_Q e_*\), while neutralizing backgrounds belong in \(J_{{\rm ext},b}^0\).

Define the integrated current
\[
J_{\rm eff}^\mu(x)\equiv \int_{-\infty}^{\infty} J^\mu(x,w)\,dw,
\]
so for \(\delta(w)\) sources, \(J_{\rm eff}^\mu=J_b^\mu\).

### 6.3 Exact integrated brane equation

Write the brane components of the EOM as
\[
\partial_M\big(ZF^{M\nu}\big)
+\frac{1}{\xi}\,\partial^\nu(\partial\!\cdot\!A)
=\mu_0 J^\nu.
\]
Integrating over \(w\) gives
\[
\partial_\mu\Big(\int Z F^{\mu\nu}dw\Big)
+\Big[ZF^{w\nu}\Big]_{-\infty}^{+\infty}
+\frac{1}{\xi}\,\partial^\nu\Big(\int (\partial\!\cdot\!A)dw\Big)
=\mu_0 J_{\rm eff}^\nu.
\]

For localized fields, the boundary term vanishes. Define the weighted brane field strength
\[
\bar F^{\mu\nu}(x)\equiv \frac{1}{Z_{\rm int}}\int Z(w)F^{\mu\nu}(x,w)\,dw.
\]
Then
\[
\partial_\mu\bar F^{\mu\nu}
+\frac{1}{\xi Z_{\rm int}}\,\partial^\nu\Big(\int (\partial\!\cdot\!A)dw\Big)
=\frac{\mu_0}{Z_{\rm int}}\,J_{\rm eff}^\nu.
\]
This is exact, but not yet the standard brane Maxwell system.

### 6.4 Zero-mode (Maxwell) limit and \(\mu_0^{\rm eff}\)

Impose the controlled zero-mode ansatz in axial gauge:
\[
A_\mu(x,w)\simeq a_\mu(x),
\qquad
\partial_w A_\mu\simeq 0,
\qquad
A_w=0.
\]
Then
\[
F^{\mu\nu}(x,w)\simeq f^{\mu\nu}(x),
\qquad
f_{\mu\nu}=\partial_\mu a_\nu-\partial_\nu a_\mu,
\]
and \(\bar F^{\mu\nu}\simeq f^{\mu\nu}\). Because the five-dimensional
gauge-fixing term is unweighted in \(w\), the clean zero-mode brane equation is
obtained by imposing Lorenz gauge before reduction, so that the explicit
integrated gauge-fixing term vanishes. One then gets
\[
\partial_\mu f^{\mu\nu}
=\mu_0^{\rm eff} J_b^\nu,
\qquad
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
\]
For Gaussian localization,
\[
\mu_0^{\rm eff}=\frac{\mu_0}{\lambda\sqrt\pi}.
\]
A gauge-fixed \(3+1\) form may be chosen afterward inside the reduced theory; it
is not obtained directly by replacing the unweighted five-dimensional
gauge-fixing integral with \(Z_{\rm int}\,\partial_\mu a^\mu\).

### 6.5 Canonical normalization and thickness-controlled charge strength

The reduced zero-mode action has schematic form
\[
S_{\rm EM}^{\rm eff}
\sim -\frac{Z_{\rm int}}{4\mu_0}\int f_{\mu\nu}f^{\mu\nu}
+\int a_\mu J_b^\mu.
\]
Rescale the brane zero mode by
\[
a_\mu=\frac{a_\mu^{\rm can}}{\sqrt{Z_{\rm int}}}.
\]
Then the canonically normalized source coupling reads
\[
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
\qquad
 e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}}.
\]
For the Gaussian profile,
\[
Z_{\rm int}=\sqrt\pi\,\lambda,
\qquad
 e_{\rm eff}\propto \lambda^{-1/2}.
\]
So thicker localization weakens the brane-observable electric coupling. The microscopic label \(q_*\) itself does **not** vary.

### 6.6 Consistency of the \(N=w\) equation

The transverse component of the bulk EOM is
\[
\partial_M\big(ZF^{Mw}\big)
+\frac{1}{\xi}\,\partial^w(\partial\!\cdot\!A)
=\mu_0 J^w.
\]
In axial gauge with the zero-mode ansatz, \(F^{\mu w}\simeq 0\) and \(\partial^w(\partial\!\cdot\!A)\simeq 0\), so consistency requires \(J^w\simeq 0\).

This is a **far-field Maxwell-limit condition**. It suppresses the mixed sector \((A_w,J^w,F_{\mu w})\) near the brane far field, but it should not be read as a microscopic claim that a charged defect core has no mixed or odd-\(w\) structure.

### 6.7 What causes deviations from 3+1 Maxwell?

The exact integrated equation shows deviations appear when any of the following fail:

1. \(A_\mu\) has significant \(w\)-dependence (KK excitations).
2. \(ZF^{w\nu}\) does not decay fast enough for the boundary term to vanish.
3. \(J^w\neq0\) (transverse current / leakage).


## 7) Gaussian KK tower and corrections to Coulomb (static sector)

### 7.1 KK mode expansion and transverse Sturm-Liouville problem

In axial gauge and Lorenz gauge for brane components, expand
\[
A_\mu(x,w)=\sum_{n\ge0} a_\mu^{(n)}(x)f_n(w).
\]
The transverse profiles solve
\[
-\frac{d}{dw}\Big(Z(w)\frac{df_n}{dw}\Big)=m_n^2 Z(w)f_n(w).
\]

### 7.2 Closed-form solution for Gaussian localization

For \(Z(w)=e^{-w^2/\lambda^2}\), one gets
\[
f_n(w)=H_n\!\left(\frac{w}{\lambda}\right),
\qquad
m_n^2=\frac{2n}{\lambda^2},
\qquad
n=0,1,2,\ldots
\]
with weighted norm
\[
\|f_n\|^2=\int Z f_n^2 dw = \lambda\sqrt\pi\,2^n n!.
\]

### 7.3 Brane couplings and parity selection rule

Define the brane coupling weights
\[
c_n\equiv \frac{f_n(0)^2}{\|f_n\|^2}.
\]
Hermite parity implies
\[
c_{2m+1}=0.
\]
For even modes,
\[
c_{2m}=\frac{1}{\lambda\sqrt\pi}\,\frac{1}{4^m}\binom{2m}{m},
\qquad
c_0=\frac{1}{\lambda\sqrt\pi}=\frac{1}{Z_{\rm int}}.
\]

Interpretation: odd transverse structure decouples from centered brane sources in the leading Maxwell reduction, but it is not thereby removed from the microscopic ontology.

### 7.4 Effective brane propagator (momentum space)

The brane-to-brane propagator is
\[
D_{\rm eff}(k^2)=\sum_{n\ge0}\frac{c_n}{k^2+m_n^2+i\epsilon},
\qquad
m_n^2=\frac{2n}{\lambda^2}.
\]

### 7.5 Coulomb + Yukawa potential for a static defect branch

For a fixed brane defect branch with
\[
J_b^0=q_*\delta^{(3)}(\mathbf x),
\]
the brane scalar potential is
\[
A_0(r)=\frac{\mu_0 q_*}{4\pi r}\sum_{n\ge0} c_n e^{-m_n r}.
\]
Separating the zero mode,
\[
A_0(r)=\frac{\mu_0^{\rm eff} q_*}{4\pi r}
\left[
1+\sum_{m\ge1}\frac{c_{2m}}{c_0}e^{-m_{2m}r}
\right],
\qquad
m_{2m}=\frac{2\sqrt m}{\lambda}.
\]
Leading correction:
\[
A_0(r)=\frac{\mu_0^{\rm eff} q_*}{4\pi r}
\left[
1+\frac12 e^{-2r/\lambda}+\mathcal O\!\left(e^{-2\sqrt2\,r/\lambda}\right)
\right].
\]

Equivalently, after canonical normalization one may rewrite the same brane field in terms of \(q_{\rm eff}\). In either normalization, the Yukawa tower modifies the field of a fixed topological charge branch; it does **not** imply that electric charge varies with circulation, throat radius, or breathing.


## 8) Moving sources and retarded structure (time-dependent sector)

### 8.1 Mode-by-mode brane equations

For conserved brane sources, each KK mode satisfies
\[
(\square_4+m_n^2)a_\mu^{(n)}(x)=\mu_0 c_n J_{b\,\mu}(x),
\qquad
\square_4\equiv \partial_t^2-\nabla^2.
\]

### 8.2 Retarded Green functions and effective retarded kernel

Let \(G^{(m)}_{\rm ret}\) solve
\[
(\square_4+m^2)G^{(m)}_{\rm ret}(x-x')=\delta^{(4)}(x-x'),
\qquad
G^{(m)}_{\rm ret}=0\ \text{for}\ t<t'.
\]
Then
\[
a_\mu^{(n)}(x)=\mu_0 c_n\int d^4x'\,G^{(m_n)}_{\rm ret}(x-x')\,J_{b\,\mu}(x').
\]

Define the KK-summed retarded kernel
\[
G^{\rm eff}_{\rm ret}(x)=\sum_{n\ge0} c_n G^{(m_n)}_{\rm ret}(x),
\]
so
\[
a_\mu(x)=\mu_0\int d^4x'\,G^{\rm eff}_{\rm ret}(x-x')\,J_{b\,\mu}(x').
\]

### 8.3 Explicit causal support: light cone + tail

Massless retarded Green function:
\[
G^{(0)}_{\rm ret}(t,r)=\theta(t)\frac{\delta(t-r)}{4\pi r}.
\]

Massive retarded Green function:
\[
G^{(m)}_{\rm ret}(t,r)
=\theta(t)\left[
\frac{\delta(t-r)}{4\pi r}
-\theta(t-r)\frac{m}{4\pi}
\frac{J_1\!\left(m\sqrt{t^2-r^2}\right)}{\sqrt{t^2-r^2}}
\right].
\]
Thus massive modes add a causal tail inside the forward light cone.

### 8.4 Moving defect-branch worldline

A moving defect on a fixed charge branch \(q_*\) with brane worldline \(z^\mu(\tau)\) has current
\[
J_b^\mu(x)=q_*\int d\tau\,u^\mu(\tau)\,\delta^{(4)}(x-z(\tau)),
\qquad
u^\mu\equiv \frac{dz^\mu}{d\tau}.
\]
Then
\[
a_\mu(x)=\mu_0 q_*\sum_{n\ge0} c_n\int d\tau\,
u_\mu(\tau)
G^{(m_n)}_{\rm ret}\big(x-z(\tau)\big).
\]
The zero mode reproduces the usual Lienard-Wiechert structure, while the massive modes add Yukawa/tail corrections.


## 9) Matter current conservation and Maxwell consistency

### 9.1 Minimal brane matter model

A representative minimally coupled brane matter Lagrangian is
\[
\mathcal L_\psi
=
\frac{i\hbar}{2}\left(\psi^* D_t\psi-\psi(D_t\psi)^*\right)
-\frac{\hbar^2}{2m}(D_i\psi)^*(D_i\psi)
-U(|\psi|^2),
\]
with Paper VII sign conventions
\[
D_t=\partial_t+\frac{i q_*}{\hbar}A_0,
\qquad
D_i=\partial_i-\frac{i q_*}{\hbar}A_i.
\]

Gauge transformation:
\[
A_0\mapsto A_0-\partial_t\chi,
\qquad
A_i\mapsto A_i+\partial_i\chi,
\qquad
\psi\mapsto e^{i q_*\chi/\hbar}\psi.
\]

This matter model is used here as a **source-consistency closure**, not as a claim that the entire ambient vacuum condensate is uniformly charged. In the defect reading, either:

* \(\psi\) is a localized, unit-normalized defect profile carrying fixed total charge \(q_*\); or
* background neutrality is implemented separately in \(J_{\rm ext}^0\).

### 9.2 Noether current and charge current

The brane Noether current is
\[
j^0=|\psi|^2,
\]
\[
j^i=\frac{\hbar}{2mi}\big(\psi^* D^i\psi-(D^i\psi)^*\psi\big)
=\frac{\hbar}{m}\,\mathrm{Im}(\psi^* D^i\psi).
\]
Define the charge current
\[
J_\psi^\mu\equiv q_* j^\mu.
\]

### 9.3 Off-shell identity and on-shell conservation

The local phase symmetry implies the off-shell Noether identity
\[
\partial_\mu j^\mu
+\frac{i}{\hbar}\big(\psi^*\mathrm{EL}_{\psi^*}-\psi\,\mathrm{EL}_\psi\big)=0.
\]
On shell,
\[
\partial_\mu j^\mu=0
\qquad\Leftrightarrow\qquad
\partial_\mu J_\psi^\mu=0.
\]

### 9.4 Embedding as a bulk source and closing the loop

Use the brane-localized bulk current model
\[
J^\mu(x,w)=J_\psi^\mu(x)\,\delta(w),
\qquad
J^w(x,w)=0.
\]
Then
\[
\partial_M J^M
=\partial_\mu(J_\psi^\mu\delta(w))
=\delta(w)\,\partial_\mu J_\psi^\mu.
\]
So brane conservation implies bulk conservation, exactly the condition required by the Maxwell divergence identity.


## 10) Validity regime, leading deviations, and falsifiable tests

### 10.1 When the 3+1 Maxwell limit holds

A sufficient set of conditions is:

1. **Zero-mode dominance:** \(A_\mu(x,w)\approx a_\mu(x)\) over the support of \(Z\), with \(\partial_w A_\mu\approx0\).
2. **No transverse flux:** \(J^w=0\).
3. **Localized boundary behavior:** \(Z(w)F^{w\nu}\to0\) fast enough to drop the boundary term.
4. **Conserved brane current:** \(\partial_\mu J_b^\mu=0\).

These are controlled **far-field Maxwell assumptions**. They suppress \((A_w,J^w,F_{\mu w})\), but do not erase those channels from the microscopic defect ontology.

### 10.2 Leading deviation structure (Gaussian profile)

Static potential:
\[
A_0(r)=\frac{\mu_0^{\rm eff}q_*}{4\pi r}
\left[1+\tfrac12 e^{-2r/\lambda}+\cdots\right].
\]
Rigid Gaussian predictions:

* the first brane-coupled mass is \(m_2=2/\lambda\);
* the leading Yukawa amplitude is fixed to \(1/2\).

### 10.3 Time-dependent thresholds and dispersion

Each coupled KK mode has dispersion
\[
\omega^2=k^2+m_n^2,
\]
so massive modes have subluminal group velocity. The first brane-coupled threshold is \(\omega=m_2=2/\lambda\).

### 10.4 Tail signatures

Massive modes add causal after-response inside the light cone. Time-domain pulse experiments constrain or potentially detect this structure.

### 10.5 Charge leakage / transverse current

The Maxwell regime assumes \(J^w=0\). The framework does not forbid \(J^w\neq0\), but then one has left the controlled brane Maxwell limit and must keep the mixed equations explicitly.

### 10.6 Profile-shape dependence

All KK predictions depend on \(Z(w)\) only through the Sturm-Liouville problem and the resulting \(\{m_n,c_n\}\). Different localization profiles would give different discrete patterns.


## 11) Referee verification suite and common pitfalls

### 11.1 Wolfram Language (WL) referee suite

The paper includes WL scripts that verify the nontrivial identities by checking that residuals simplify to zero under the frozen conventions.

Scripts summarized in the paper:

* `maxwell_from_4d.wl` -- vary the localized action, derive the EOM for all components, check the Bianchi identities and divergence identity, verify the brane reduction, and print the updated bookkeeping \(q_* = \eta_Q e_*\), \(J_{\rm tot}=J_\psi+J_{\rm ext}\), and \(q_{\rm eff}=q_*/\sqrt{Z_{\rm int}}\).
* `maxwell_symmetry_checks.wl` -- check gauge variation of the interaction term, the shift of \(\partial\!\cdot\!A\), the brane Lorentz transformation types, and the Paper VII matter-gauge convention.
* `maxwell_kk_coulomb_propagator.wl` -- solve the Gaussian SL problem, verify the Hermite spectrum, norms, couplings \(c_n\), parity rule, Coulomb+Yukawa structure for a fixed defect branch, and the \(k^2\)-only propagator dependence.
* `maxwell_axial_gauge_brane_sources.wl` -- verify axial gauge reachability, residual 3+1 gauge freedom, brane \(\delta(w)\) sources, integrated equations, and the controlled-far-field reading of \(J^w=0\).
* `maxwell_moving_charges_propagator.wl` -- verify that the time-dependent response remains Lorentz covariant and depends only on \(k^2\), with source language interpreted as a fixed defect branch rather than an arbitrary charge.
* `matter_current_conservation.wl` -- derive the brane Noether current with the Paper VII sign convention, verify the off-shell identity, and record \(J_\psi^\mu=q_* j^\mu\).

### 11.2 Common pitfalls explicitly flagged

1. Treating \(A_\mu\) as a vector instead of a covector breaks Lorentz checks.
2. Mixing covariant and contravariant components without explicit metric raising/lowering produces sign errors.
3. WL `EulerEquations` output is safest when converted to residual form before simplification.
4. The Noether identity has a fixed relative sign; flipping it generates fake failures.
5. Hermite norms are most robustly checked by integer spot-checks plus validated closed forms.


## 12) Minimal cache for future work

If you want the minimal working set, keep the following.

1. **Action**
   \[
   S=\int d^5x\left[-\frac{Z}{4\mu_0}F_{MN}F^{MN}-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2-A_MJ^M\right].
   \]
   In the coupled reading, treat \(J^M\) as shorthand total source and split \(J_{\rm tot}=J_\psi+J_{\rm ext}\).

2. **Charge ontology**
   \[
   q_*=\eta_Q e_*,
   \qquad
   q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}},
   \qquad
   e_{\rm eff}=\frac{e_*}{\sqrt{Z_{\rm int}}}.
   \]
   Sign is fixed by \(\eta_Q\); thickness controls the reduced brane coupling.

3. **Bulk EOM + divergence identity**
   \[
   \partial_M(ZF^{MN})+\xi^{-1}\partial^N(\partial\!\cdot\!A)=\mu_0J^N,
   \qquad
   \xi^{-1}\Box(\partial\!\cdot\!A)=\mu_0\partial_NJ^N.
   \]

4. **Controlled brane Maxwell limit**
   * axial gauge: \(A_w=0\), residual \(\chi=\chi_b(x)\),
   * brane source: \(J^\mu=J_b^\mu\delta(w)\), \(J^w=0\),
   * zero mode: \(\partial_w A_\mu\approx0\),
   giving
   \[
   \partial_\mu f^{\mu\nu}=\mu_0^{\rm eff}J_b^\nu,
   \qquad
   \mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}=\frac{\mu_0}{\lambda\sqrt\pi}.
   \]
   The same reduction gives \(q_{\rm eff}=q_*/\sqrt{Z_{\rm int}}\).

5. **Gaussian KK data**
   \[
   m_n^2=\frac{2n}{\lambda^2},
   \qquad
   c_{2m+1}=0,
   \qquad
   c_{2m}=\frac{1}{\lambda\sqrt\pi}\frac{1}{4^m}\binom{2m}{m}.
   \]
   Static potential:
   \[
   A_0(r)=\frac{\mu_0 q_*}{4\pi r}\sum_n c_n e^{-m_n r}
   =\frac{\mu_0^{\rm eff}q_*}{4\pi r}\left[1+\frac12 e^{-2r/\lambda}+\cdots\right].
   \]

6. **Matter current source**
   \[
   j^0=|\psi|^2,
   \qquad
   j^i=\frac{\hbar}{m}\,\mathrm{Im}(\psi^*D^i\psi),
   \qquad
   J_\psi^\mu=q_* j^\mu.
   \]
   Embed as \(J^\mu=J_\psi^\mu\delta(w)\), \(J^w=0\), with any neutralizing background placed in \(J_{\rm ext}^0\).

That is enough to recover: the localized Maxwell sector, the brane reduction, the thickness-controlled brane charge strength, the Gaussian Yukawa correction pattern, the causal retarded structure, and the matter-current consistency closure.
