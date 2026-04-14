# Moving-Throat PDE — Stage 186: Free-Quintuple Scalar Closure Slice, Exact Graph-Family Tangency, and the One-Parameter Crossing Theorem

## Status

**Exact within the carried Stage-175 orbit/quotient projector calculus and the Stage-185 exact free-quintuple target graph.**

This stage does **not** introduce a new constitutive law.
It takes the Stage-185 target-orbit graph and proves that, once a candidate family is aligned with that graph, the full reduced home-stretch test collapses to a **single scalar closure function** on the five free microscopic coordinates.

---

## Purpose

Stage 185 solved the target orbit `\(\mathcal O_*\)` explicitly as a graph over the free quintuple
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
\]
That already removed the abstractness from the Stage-184 canonical projection. But one practical question was still left open:

> once the target orbit is written as an explicit graph, what is the actual reduced closure problem for a one-parameter branch family in that free-quintuple space?

This stage answers that exactly.

The main outputs are:

1. the exact **graph-closure scalar**
   \[
   \widehat\chi_Q(\mathbf y):=\chi_Q\bigl(\mathbf x_*^{\rm graph}(\mathbf y)\bigr),
   \]
2. the exact theorem that the full reduced closure set is the codimension-one graph slice
   \[
   \mathcal Z_*=
   \{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\chi_Q(\mathbf y)=1\},
   \]
3. the exact formulas for the dependent-triple log tangents of any graph-lifted free family,
4. the exact theorem that those graph tangents lie in the Stage-175 similarity-orbit kernel,
5. the exact same-free-quintuple decomposition of an arbitrary candidate family into
   - the graph-lifted orbit piece, and
   - the three dependent graph errors \((E_T,E_K,E_\mu)\),
6. and the exact **one-parameter crossing theorem** reducing the graph-aligned closure search to a single scalar sign-change problem.

So Stage 186 is the first point where a reduced moving-throat branch search becomes a genuine **one-scalar spectral-placement problem** on the free-quintuple target graph.

---

## 1. Carry-forward exact free-quintuple target graph

Let the positive microscopic state be
\[
\boxed{
\mathbf x=
\bigl(
\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U
\bigr).
}
\]

Carry forward the Stage-185 free-quintuple projection
\[
\boxed{
\mathbf y:=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
}
\]

The dependent triple is
\[
(T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W).
\]

The exact Stage-185 target-graph map is
\[
\boxed{
\mathbf x_*^{\rm graph}(\mathbf y)
:=
\bigl(
\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,
K_{\eta,*}^{\rm graph}(\mathbf y),
K_W^{(\mathrm{eff})},
\mu_W^{\rm graph}(\mathbf y),
T_U^{\rm graph}(\mathbf y)
\bigr),
}
\]
with
\[
\boxed{
\delta_{U,*}^{\rm graph}(\mathbf y)
:=
\left[
\frac{\mathfrak C_{{\rm tr},*}^{\rm target}}
{\left(\dfrac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}}
\right]^{1/(1+\chi_{0,*})},
}
\]
\[
\boxed{
T_U^{\rm graph}(\mathbf y)
=
\frac{L^2K_U}{\pi^2}\,\delta_{U,*}^{\rm graph}(\mathbf y),
}
\]
\[
\boxed{
K_{\eta,*}^{\rm graph}(\mathbf y)
=
\frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}^{\rm target}},
}
\]
\[
\boxed{
\mu_W^{\rm graph}(\mathbf y)
=
\frac{\mathfrak C_{{\rm nt},*}^{\rm target}\,c_{\eta U}^2\,(K_W^{(\mathrm{eff})})^2}
{\epsilon_{\eta,*}^{\rm target}\,K_U\,\lambda_W^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{-E_*}
\left(\delta_{U,*}^{\rm graph}(\mathbf y)\right)^{F_*}.
}
\]

By Stage 185,
\[
\boxed{
\mathcal O_*=
\{\mathbf x_*^{\rm graph}(\mathbf y):\mathbf y\in(\mathbb R_{>0})^5\}.
}
\]

---

## 2. Exact graph-closure scalar and the codimension-one closure slice

Define the exact graph-closure scalar
\[
\boxed{
\widehat\chi_Q(\mathbf y)
:=
\chi_Q\bigl(\mathbf x_*^{\rm graph}(\mathbf y)\bigr).
}
\]
Equivalently, define the graph residual
\[
\boxed{
\widehat\Delta_Q(\mathbf y):=\widehat\chi_Q(\mathbf y)-1.
}
\]

Since every graph state already lies on `\(\mathcal O_*\)`, the full reduced closure set is obtained by imposing only the Packet-A scalar target on the graph.

\[
\boxed{\textbf{Theorem (Stage 186 scalar graph-slice theorem).}}
\]

The fully reduced closure set is exactly
\[
\boxed{
\mathcal Z_*
=
\Bigl\{
\mathbf x_*^{\rm graph}(\mathbf y):
\widehat\chi_Q(\mathbf y)=1
\Bigr\}
=
\Bigl\{
\mathbf x_*^{\rm graph}(\mathbf y):
\widehat\Delta_Q(\mathbf y)=0
\Bigr\}.
}
\]

### Proof

By Stage 185, `\(\mathbf x_*^{\rm graph}(\mathbf y)\in\mathcal O_*\)` for every positive free quintuple `\(\mathbf y\)`.
Therefore the only remaining reduced closure condition is the Packet-A scalar condition
\[
\chi_Q=1.
\]
Evaluating that on the graph gives exactly the displayed slice condition.

So once the target orbit has been solved as a graph, the reduced closure set is a codimension-one scalar slice of that graph.

---

## 3. Exact graph-family tangent formulas

Let
\[
\boxed{\mathbf y(\tau)}
\]
be a smooth one-parameter free-quintuple family, and define the corresponding graph-lifted microscopic family
\[
\boxed{
\mathbf x_{\rm graph}(\tau):=\mathbf x_*^{\rm graph}(\mathbf y(\tau)).
}
\]

Introduce the free log-tangent components
\[
\boxed{
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_W)
:=
\left(
\frac{d\ln\lambda_W}{d\tau},
\frac{d\ln c_{\eta U}}{d\tau},
\frac{d\ln\gamma}{d\tau},
\frac{d\ln K_U}{d\tau},
\frac{d\ln K_W^{(\mathrm{eff})}}{d\tau}
\right).
}
\]

Then the graph formulas imply the exact dependent-triple log tangents
\[
\boxed{
\frac{d\ln \delta_{U,*}^{\rm graph}}{d\tau}
=
-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
\bigl(\gamma_1+c_1-\kappa_U\bigr),
}
\]
\[
\boxed{
\tau_1^{\rm graph}
:=
\frac{d\ln T_U^{\rm graph}}{d\tau}
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
\bigl(\gamma_1+c_1-\kappa_U\bigr),
}
\]
\[
\boxed{
\kappa_{\eta}^{\rm graph}
:=
\frac{d\ln K_{\eta,*}^{\rm graph}}{d\tau}
=
2c_1-\kappa_U,
}
\]
\[
\boxed{
\mu_1^{\rm graph}
:=
\frac{d\ln\mu_W^{\rm graph}}{d\tau}
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-E_*\bigl(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W\bigr)
-F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}\bigl(\gamma_1+c_1-\kappa_U\bigr).
}
\]

So the full graph-family log tangent is
\[
\boxed{
\dot{\Delta\mathbf x}_{\rm graph}
=
\begin{pmatrix}
\lambda_1\\
c_1\\
\gamma_1\\
\kappa_U\\
\kappa_{\eta}^{\rm graph}\\
\kappa_W\\
\mu_1^{\rm graph}\\
\tau_1^{\rm graph}
\end{pmatrix}.
}
\]

---

## 4. Exact graph-family tangency to the Stage-175 orbit kernel

Carry forward the exact Stage-175 monomial-drift map
\[
\boxed{
\mathbf q=M_*\,\Delta\mathbf x,
}
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

Then the graph-family tangent satisfies the exact kernel identity
\[
\boxed{
M_*\,\dot{\Delta\mathbf x}_{\rm graph}=0.
}
\]

\[
\boxed{\textbf{Corollary (graph-family orbit tangency).}}
\]

Every graph-lifted free-quintuple family lies identically in the Stage-175 similarity-orbit tangent:
\[
\boxed{
\dot{\Delta\mathbf x}_{\rm graph}(\tau)\in\ker M_*.
}
\]

So the entire Packet-B part of the reduced home-stretch theorem vanishes identically on graph-aligned families. That is the exact reason the closure search drops from four scalars to one scalar on the graph.

---

## 5. Same-free-quintuple candidate decomposition and exact graph errors

Let
\[
\boxed{
\mathbf x_{\rm cand}(\tau)
}
\]
be an arbitrary candidate family sharing the same free-quintuple path `\(\mathbf y(\tau)\)` as the graph family.
Write its dependent triple against the graph lift as
\[
\boxed{
T_U(\tau)=T_U^{\rm graph}(\tau)e^{E_T(\tau)},
\qquad
K_\eta^{(\mathrm{eff})}(\tau)=K_{\eta,*}^{\rm graph}(\tau)e^{E_K(\tau)},
\qquad
\mu_W(\tau)=\mu_W^{\rm graph}(\tau)e^{E_\mu(\tau)}.
}
\]

Then the exact same-free-quintuple error vector is
\[
\boxed{
\Delta\mathbf x_{\rm err}
=
\begin{pmatrix}
0\\
0\\
0\\
0\\
E_K\\
0\\
E_\mu\\
E_T
\end{pmatrix}.
}
\]
Applying the Stage-175 monomial map gives the exact quotient packet carried by those graph errors:
\[
\boxed{
q_{\rm tr}=(1+\chi_{0,*})E_T,
}
\]
\[
\boxed{
q_{\rm nt}=E_\mu-E_K-F_*E_T,
}
\]
\[
\boxed{
q_\eta=-E_K.
}
\]

The inverse formulas are therefore
\[
\boxed{
E_T=\frac{q_{\rm tr}}{1+\chi_{0,*}},
\qquad
E_K=-q_\eta,
\qquad
E_\mu=q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}.
}
\]

This is exactly the Stage-185 graph-error packet, now written as the unique same-free-quintuple decomposition of any candidate family into

- the graph-lifted orbit piece, and
- the genuine quotient-failure piece.

---

## 6. Exact family packet on and off the graph

For a general same-free-quintuple candidate family, the exact reduced four-scalar packet is
\[
\boxed{
\Delta_{\rm fam}^{\rm graph}(\tau)
=
\bigl(
\chi_Q(\mathbf x_{\rm cand}(\tau))-1,
E_T(\tau),
E_K(\tau),
E_\mu(\tau)
\bigr).
}
\]

So the full reduced closure criterion remains
\[
\boxed{
\mathbf x_{\rm cand}(\tau)\in\mathcal Z_*
\iff
\chi_Q(\mathbf x_{\rm cand}(\tau))=1,
\ E_T(\tau)=0,
\ E_K(\tau)=0,
\ E_\mu(\tau)=0.
}
\]

But on the graph-lifted family itself one has
\[
E_T\equiv E_K\equiv E_\mu\equiv 0,
\]
so the packet collapses to
\[
\boxed{
\Delta_{\rm fam}^{\rm graph,\,lift}(\tau)
=
\bigl(
\widehat\chi_Q(\mathbf y(\tau))-1,
0,
0,
0
\bigr).
}
\]

This is the exact one-scalar reduction of the home-stretch family test.

---

## 7. One-parameter graph crossing theorem

Let
\[
\boxed{
\mathbf y:[\tau_-,\tau_+]\to(\mathbb R_{>0})^5
}
\]
be a continuous free-quintuple path, and let
\[
\mathbf x_{\rm graph}(\tau)=\mathbf x_*^{\rm graph}(\mathbf y(\tau)).
\]
Define the graph scalar residual
\[
\boxed{
\widehat\Delta_Q(\tau):=
\widehat\chi_Q(\mathbf y(\tau))-1.
}
\]

\[
\boxed{\textbf{Theorem (Stage 186 one-parameter graph crossing theorem).}}
\]

If
\[
\widehat\Delta_Q(\tau_-)\,\widehat\Delta_Q(\tau_+)<0,
\]
then there exists
\[
\boxed{
\tau_*\in(\tau_-,\tau_+)
}
\]
such that
\[
\boxed{
\widehat\Delta_Q(\tau_*)=0,
\qquad
\mathbf x_{\rm graph}(\tau_*)\in\mathcal Z_*.
}
\]

### Proof

The graph residual `\(\widehat\Delta_Q(\tau)\)` is continuous because both the free path and the graph lift are continuous and `\(\chi_Q\)` is continuous on the carried reduced branch. So the conclusion follows immediately from the intermediate value theorem.

This is the first exact existence theorem for a closed reduced state along a one-parameter graph-aligned family.

---

## 8. Simple crossing and tangency on the graph

Let `\(\tau_*\)` be a graph-family closure point with
\[
\widehat\Delta_Q(\tau_*)=0.
\]
Then:

### 8.1 Transverse graph crossing
If
\[
\boxed{
\frac{d\widehat\Delta_Q}{d\tau}(\tau_*)\neq 0,
}
\]
then the closure point is a simple transverse crossing of the scalar graph slice. In particular, it is locally isolated on the one-parameter family.

### 8.2 Graph tangency
If
\[
\boxed{
\frac{d\widehat\Delta_Q}{d\tau}(\tau_*)=0,
}
\]
then the graph family is tangent to the codimension-one closure slice at first order.

Because the Packet-B part vanishes identically on graph families, these are the only two first-order possibilities.

So the first-order family audit on the graph is reduced completely to the derivative of one scalar function.

---

## 9. Operational meaning for the PDE program

Stage 185 reduced the target orbit from an abstract similarity set to an explicit graph over five free coordinates.
Stage 186 goes one decisive step further.

### 9.1 Graph-aligned PDE searches are one-scalar searches

If a reduced moving-throat solve is organized as a graph-aligned family in the free quintuple, then the whole home-stretch packet becomes
\[
(\widehat\chi_Q-1,0,0,0).
\]
So the closure search is no longer a four-scalar search. It is a one-scalar search.

### 9.2 Same-free-quintuple candidate families are still fully explicit

If the candidate family is not exactly graph-aligned, there is still no abstract quotient problem left.
The failure is exactly the dependent graph-error triple
\[
(E_T,E_K,E_\mu),
\]
and those errors generate the exact quotient packet
\[
(q_{\rm tr},q_{\rm nt},q_\eta)

automatically.
\]

### 9.3 The remaining reduced theorem gap is now sharply executable

The completed moving-throat PDE no longer needs to be asked for a vague “final closure.”
It needs to do one of two concrete things:

1. either return a graph-aligned free-quintuple family and let one solve
   \[
   \widehat\chi_Q(\mathbf y)=1,
   \]
2. or return a same-free-quintuple candidate family and let one solve the explicit four-scalar packet
   \[
   (\chi_Q-1,E_T,E_K,E_\mu)=0.
   \]

That is the smallest clean family-level theorem gate reached so far.
