# Parent Throat Action — Runtime Monitor and Hard Falsifiers

## Purpose

Once the Branch-B patch is treated as a hardcoded CFD boundary condition, the
right question is no longer how to derive it.

The right question is:

> what exactly should be monitored in the grid, and what signatures kill the
> gravity analogy immediately?

This step compiles that monitor suite from the exact 4D brane-projection and
bridge identities.

---

## 1. What to monitor directly

The primary raw monitors are the exact projected brane quantities:

\[
\rho_{\rm brane}=\int W(w)\rho\,dw,
\qquad
\mathbf J_{\rm brane}=\int W(w)\mathbf J_{xyz}\,dw,
\]

and especially the projected leakage source

\[
S_\rho
=
-[WJ^w]_{-\infty}^{+\infty}
+\int W'(w)J^w\,dw.
\]

If the boundary term vanishes,

\[
S_\rho=\int W'(w)J^w\,dw.
\]

This is the raw source that drives the gravity hook. It is the one to monitor.

You should also record

\[
\mathbf v_{\rm brane}=\frac{\mathbf J_{\rm brane}}{\rho_{\rm brane}},
\]

along with the projection non-commutativity stress

\[
R_{ij}=\Pi_{ij}-\rho_{\rm brane}v_i v_j,
\]

and the momentum leakage correction

\[
S_{J_i}
=
-\Big[W\frac{J_iJ_w}{\rho}\Big]
+\int W'(w)\frac{J_iJ_w}{\rho}\,dw.
\]

These are regime guards: if they are large, the clean 3D gravity-like reduction
is not operating.

---

## 2. The exact and approximate gravity monitors

The exact open-system continuity law is

\[
\partial_t\rho_{\rm brane}+\nabla_3\cdot \mathbf J_{\rm brane}=S_\rho.
\]

So the exact continuity residual is

\[
R_{\rm cont}
:=
\partial_t\rho_{\rm brane}+\nabla_3\cdot \mathbf J_{\rm brane}-S_\rho.
\]

Using

\[
\mathbf v_{\rm brane}=\nabla_3\phi_3+\mathbf v_T,
\qquad \nabla_3\cdot\mathbf v_T=0,
\]

the exact longitudinal identity is

\[
\rho_{\rm brane}\nabla_3^2\phi_3
=
S_\rho-\partial_t\rho_{\rm brane}
-(\nabla_3\rho_{\rm brane})\cdot(\nabla_3\phi_3+\mathbf v_T).
\]

So the exact Poisson residual is

\[
R_{\rm Pois}^{\rm exact}
:=
\rho_{\rm brane}\nabla_3^2\phi_3
-S_\rho
+\partial_t\rho_{\rm brane}
+(\nabla_3\rho_{\rm brane})\cdot(\nabla_3\phi_3+\mathbf v_T).
\]

In the linearized quasi-static regime,

\[
\rho_0\nabla_3^2\phi_3\approx S_\rho,
\]

so the reduced runtime check is

\[
R_{\rm Pois}^{\rm lin}:=\rho_0\nabla_3^2\phi_3-S_\rho.
\]

These two residuals should stay small before any gravity-like interpretation is
taken seriously.

---

## 3. Massive-scalar hard fail

The clean Newtonian exterior candidate is

\[
\phi_N(r)=-\frac{A}{r},
\]

for which

\[
Q_r(r):=4\pi r^2\partial_r\phi_N=4\pi A
\]

is exactly constant.

By contrast, a Yukawa-screened exterior

\[
\phi_Y(r)=-\frac{A e^{-\mu r}}{r}
\]

gives

\[
Q_r(r)=4\pi A e^{-\mu r}(1+\mu r),
\]

which is not constant.

The corresponding exterior mass diagnostic is

\[
\mu_{\rm eff}^2(r)
:=
\frac{\nabla_3^2\phi_3-S_\rho/\rho_0}{\phi_3}.
\]

For the Newtonian exterior,

\[
\mu_{\rm eff}^2=0.
\]

For the Yukawa exterior,

\[
\mu_{\rm eff}^2=\mu^2.
\]

The script also checks non-Newtonian impostors explicitly:

\[
\phi(r)=-\frac{A}{r^2}
\quad\Rightarrow\quad
Q_r=\frac{8\pi A}{r},
\qquad
\mu_{\rm eff}^2=\frac{2}{r^2},
\]

and

\[
\phi(r)=A\log r
\quad\Rightarrow\quad
Q_r=4\pi A r,
\qquad
\mu_{\rm eff}^2=\frac{1}{r^2\log r}.
\]

Both fail the flux-plateau/massless-exterior checks.

So the first hard fail is:

> if \(Q_r(r)\) does not plateau, or \(\mu_{\rm eff}^2(r)\) approaches a
> nonzero exterior constant, the simulation is not mimicking Newtonian gravity.

---

## 4. Optical/redshift hard fail

The 4D bridge summary fixed the weak-field optics coefficient to

\[
\alpha_n=\frac{n-1}{2}=2
\qquad\text{since }n=5.
\]

So the weak-field probe sector must obey the coefficient-2 law. A clean runtime
fit is

\[
\alpha_{\rm fit}(r)
:=
-\frac{c_{\rm probe}^2\,[N_{\rm probe}(r)-1]}{\Phi_{\rm eff}(r)},
\]

where \(\Phi_{\rm eff}\) is the effective exterior potential inferred from the
same source sector used to define the longitudinal field.

Then the hard redshift/optics condition is

\[
\alpha_{\rm fit}\approx 2.
\]

Equivalently, the analog bending and Shapiro coefficients must match the same
weak-field number:

\[
\Delta\theta=\frac{4GM_{\rm eff}}{bc_{\rm probe}^2},
\qquad
\Delta t \propto 2\frac{GM_{\rm eff}}{c_{\rm probe}^3}.
\]

So the second hard fail is:

> if the fitted optical/redshift coefficient does not plateau near \(2\), the
> simulation is not reproducing the \(n=5\) gravity-like optics sector.

---

## 5. Tail check

If the 4PN tail scalar is available independently, the compact note adds

\[
R_{\rm tail}:=\Theta_{\rm tail}(c/c_s)^3-1.
\]

This is optional, but if measured it is another direct failure channel.

---

## 6. Concrete synthetic exhibit

The symbolic checks are backed by concrete grid snapshots from the runtime
postprocessor self-test:

- periodic exact-consistency snapshot:
  - `max |R_cont| = 1.4210854715202004e-14`
  - `max |R_Pois_exact| = 0.0`
  - `alpha_fit tail = 2.0000000000000004`
- Newton exterior:
  - `Q_r` tail coefficient of variation
    `= 0.00048311949183278685`
  - `mu_eff^2` tail median `= -0.003597967054167384`
- Yukawa exterior:
  - `Q_r` tail coefficient of variation
    `= 0.25427545643672567`
  - `mu_eff^2` tail median `= 1.954292666547206`
- bad-optics exterior:
  - `alpha_fit tail = 1.4`
- `1/r^2` impostor exterior:
  - `Q_r` tail coefficient of variation
    `= 0.10349929749560405`
  - `mu_eff^2` tail median `= 0.3602269975484812`
  - fail-fast classifier verdict: `FAIL`
- `log r` impostor exterior:
  - `Q_r` tail coefficient of variation
    `= 0.10301948878406167`
  - `mu_eff^2` tail median `= 0.21968065668811232`
  - fail-fast classifier verdict: `FAIL`

This makes the runtime monitor test more than a set of definition
substitutions: it exercises the same numerical path that will be used on CFD
snapshots.

---

## 7. Interpretation

The monitor hierarchy is now clear.

1. **Raw source:** \(S_\rho\) from projected leakage.
2. **Derived gravity-like field:** the longitudinal brane potential \(\phi_3\).
3. **Regime guards:** \(R_{ij}\), \(S_{J_i}\), \(R_{\rm cont}\), \(R_{\rm Pois}\).
4. **Hard falsifiers:** Yukawa screening and the wrong optical coefficient.

So the CFD runtime should be built around the projected leakage source and the
Helmholtz-extracted longitudinal potential, not around the Stage-45 support
field \(\phi\).
