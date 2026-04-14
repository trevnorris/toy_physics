# Moving-Throat PDE — Stage 187: Explicit Free-Quintuple Log-Ray Sweep, Finite Graph Invariance, and the Scalar Root Predictor

## Status

**Exact within the carried Stage-185 free-quintuple target graph and the Stage-186 scalar graph-slice theorem.**

This stage does **not** introduce a new constitutive law.
It turns the abstract one-parameter family of Stage 186 into an explicit, audit-ready family generator on the five free microscopic coordinates and reduces the local branch search to a scalar root problem on a logarithmic ray.

---

## Purpose

Stage 186 proved that once the target orbit is written as the exact free-quintuple graph
\[
\mathbf x_*^{\rm graph}(\mathbf y),
\]
all reduced closure points are exactly the scalar slice
\[
\widehat\chi_Q(\mathbf y)=1.
\]
That is already a complete existence theorem, but it still leaves one practical question open:

> what is the smallest explicit one-parameter family in free-quintuple space that can be followed, audited, and differentiated without rebuilding the whole orbit geometry each time?

This stage answers that exactly.

The main outputs are:

1. the exact **free-quintuple logarithmic ray**
   \[
   \mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s},
   \]
2. the exact finite graph lift of that ray, with the dependent triple carried by four explicit exponents,
3. the exact theorem that the full graph-lifted ray stays on the target orbit `\(\mathcal O_*\)` for **all** `\(\tau\)`, not only at first order,
4. the explicit **primitive direction table** for the five free microscopic directions,
5. the exact scalarized closure function
   \[
   \Phi_{\mathbf s}(\tau):=\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau)),
   \]
6. the monotone-ray uniqueness theorem,
7. and the exact affine and log-linear root predictors that turn the Stage-186 crossing test into a concrete local branch-audit tool.

So Stage 187 is the first fully explicit **family compiler** for the reduced home stretch.

---

## 1. Carry-forward free quintuple and target graph

The free quintuple is
\[
\boxed{
\mathbf y:=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
}
\]
The exact Stage-185 target graph gives the dependent triple
\[
(T_U^{\rm graph}(\mathbf y),\ K_{\eta,*}^{\rm graph}(\mathbf y),\ \mu_W^{\rm graph}(\mathbf y))
\]
and therefore the full target-orbit graph point
\[
\boxed{
\mathbf x_*^{\rm graph}(\mathbf y)
=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,
K_{\eta,*}^{\rm graph}(\mathbf y),
K_W^{(\mathrm{eff})},
\mu_W^{\rm graph}(\mathbf y),
T_U^{\rm graph}(\mathbf y)).
}
\]
By Stage 186,
\[
\boxed{
\mathcal Z_*
=
\{\mathbf x_*^{\rm graph}(\mathbf y):\widehat\chi_Q(\mathbf y)=1\},
\qquad
\widehat\chi_Q(\mathbf y):=\chi_Q(\mathbf x_*^{\rm graph}(\mathbf y)).
}
\]
So the only missing practical ingredient is an explicit one-parameter family in `\(\mathbf y\)`.

---

## 2. Exact free-quintuple logarithmic ray

Fix a positive base point
\[
\boxed{
\mathbf y_\circ
=(\lambda_{W,\circ},\ c_{\eta U,\circ},\ \gamma_\circ,\ K_{U,\circ},\ K_{W,\circ}^{(\mathrm{eff})})
\in(\mathbb R_{>0})^5,
}
\]
and a logarithmic direction
\[
\boxed{
\mathbf s=(s_\lambda,s_c,s_\gamma,s_U,s_W)\in\mathbb R^5.
}
\]
The exact log-ray through `\(\mathbf y_\circ\)` is
\[
\boxed{
\mathbf y_{\mathbf s}(\tau)
:=
\mathbf y_\circ\odot e^{\tau\mathbf s},
}
\]
meaning
\[
\lambda_W(\tau)=\lambda_{W,\circ}e^{s_\lambda\tau},
\qquad
c_{\eta U}(\tau)=c_{\eta U,\circ}e^{s_c\tau},
\qquad
\gamma(\tau)=\gamma_\circ e^{s_\gamma\tau},
\]
\[
K_U(\tau)=K_{U,\circ}e^{s_U\tau},
\qquad
K_W^{(\mathrm{eff})}(\tau)=K_{W,\circ}^{(\mathrm{eff})}e^{s_W\tau}.
\]
So the free-quintuple logarithmic slopes are exactly constant:
\[
\frac{d\ln\lambda_W}{d\tau}=s_\lambda,
\qquad
\frac{d\ln c_{\eta U}}{d\tau}=s_c,
\qquad
\frac{d\ln\gamma}{d\tau}=s_\gamma,
\qquad
\frac{d\ln K_U}{d\tau}=s_U,
\qquad
\frac{d\ln K_W^{(\mathrm{eff})}}{d\tau}=s_W.
\]

### 2.1 Local completeness of the log-ray library

Let `\(\mathbf y(\tau)\)` be any smooth positive free-quintuple path, and fix `\(\tau_0\)`.
Define
\[
\boxed{
\mathbf s_{\tau_0}:=
\left.
\frac{d\ln \mathbf y}{d\tau}
\right|_{\tau_0}.
}
\]
Then
\[
\boxed{
\mathbf y(\tau_0+\Delta\tau)
=
\mathbf y(\tau_0)\odot e^{\Delta\tau\mathbf s_{\tau_0}}+O(\Delta\tau^2).
}
\]
So every smooth free-quintuple path is locally represented, to first order, by a unique log-ray.
This is why the log-ray family is the natural Stage-186 continuation for local branch auditing.

---

## 3. Exact finite graph lift of the log-ray

Introduce the carried ratio
\[
\boxed{
\mathfrak a_*:=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.
}
\]
Define the four dependent exponents
\[
\boxed{
\sigma_\delta(
\mathbf s)
:=
-\mathfrak a_*(s_\gamma+s_c-s_U),
}
\]
\[
\boxed{
\sigma_T(\mathbf s):=s_U+\sigma_\delta(\mathbf s),
\qquad
\sigma_{K_\eta}(\mathbf s):=2s_c-s_U,
}
\]
\[
\boxed{
\sigma_\mu(\mathbf s)
:=
2s_c-s_U+2s_W-2s_\lambda
-E_*(2s_\gamma+2s_\lambda-s_U-s_W)
+F_*\sigma_\delta(\mathbf s).
}
\]
Then the exact Stage-185 target graph along the ray is
\[
\boxed{
\delta_{U,*}^{\rm graph}(\tau)
=
\delta_{U,*}^{\rm graph}(0)
\,e^{\sigma_\delta\tau},
}
\]
\[
\boxed{
T_U^{\rm graph}(\tau)
=
T_{U,\circ}^{\rm graph}
\,e^{\sigma_T\tau},
}
\]
\[
\boxed{
K_{\eta,*}^{\rm graph}(\tau)
=
K_{\eta,\circ}^{\rm graph}
\,e^{\sigma_{K_\eta}\tau},
}
\]
\[
\boxed{
\mu_W^{\rm graph}(\tau)
=
\mu_{W,\circ}^{\rm graph}
\,e^{\sigma_\mu\tau}.
}
\]
So the full graph-lifted microscopic ray is
\[
\boxed{
\mathbf x_{\mathbf s}^{\rm graph}(\tau)
=
\bigl(
\lambda_{W,\circ}e^{s_\lambda\tau},
 c_{\eta U,\circ}e^{s_c\tau},
 \gamma_\circ e^{s_\gamma\tau},
 K_{U,\circ}e^{s_U\tau},
 K_{\eta,\circ}^{\rm graph}e^{\sigma_{K_\eta}\tau},
 K_{W,\circ}^{(\mathrm{eff})}e^{s_W\tau},
 \mu_{W,\circ}^{\rm graph}e^{\sigma_\mu\tau},
 T_{U,\circ}^{\rm graph}e^{\sigma_T\tau}
\bigr).
}
\]

The corresponding constant microscopic log-tangent is
\[
\boxed{
\dot{\Delta\mathbf x}_{\mathbf s}^{\rm graph}
=
(s_\lambda,s_c,s_\gamma,s_U,\sigma_{K_\eta},s_W,\sigma_\mu,\sigma_T)^T.
}
\]

---

## 4. Exact finite monomial invariance theorem on the graph ray

The three Stage-185 target monomials are
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta.
\]
By direct substitution of the log-ray graph lift,
\[
\boxed{
\mathfrak C_{{\rm tr},*}(\mathbf x_{\mathbf s}^{\rm graph}(\tau))
=
\mathfrak C_{{\rm tr},*}^{\rm target},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}(\mathbf x_{\mathbf s}^{\rm graph}(\tau))
=
\mathfrak C_{{\rm nt},*}^{\rm target},
}
\]
\[
\boxed{
\epsilon_\eta(\mathbf x_{\mathbf s}^{\rm graph}(\tau))
=
\epsilon_{\eta,*}^{\rm target},
}
\]
for **all** `\(\tau\in\mathbb R\)`.
Therefore
\[
\boxed{
\mathbf x_{\mathbf s}^{\rm graph}(\tau)\in\mathcal O_*
\qquad\text{for all }\tau.
}
\]

This is the finite version of the Stage-186 kernel theorem.
Equivalently, the exact Stage-175 quotient map annihilates the constant ray tangent:
\[
\boxed{
M_*\dot{\Delta\mathbf x}_{\mathbf s}^{\rm graph}=0.
}
\]
So the Packet-B part of the home-stretch verdict vanishes identically on every free-quintuple log-ray graph lift.

---

## 5. Primitive free-direction table

The five primitive free rays are the coordinate directions
\[
\mathbf e_\lambda,
\mathbf e_c,
\mathbf e_\gamma,
\mathbf e_U,
\mathbf e_W.
\]
For each one, the dependent exponents are exact.
Let
\[
\sigma_\delta\sim \delta_{U,*}^{\rm graph},
\qquad
\sigma_T\sim T_U^{\rm graph},
\qquad
\sigma_{K_\eta}\sim K_{\eta,*}^{\rm graph},
\qquad
\sigma_\mu\sim \mu_W^{\rm graph}.
\]
Then:

| Primitive ray | `\(\sigma_\delta\)` | `\(\sigma_T\)` | `\(\sigma_{K_\eta}\)` | `\(\sigma_\mu\)` |
|---|---:|---:|---:|---:|
| `\(\mathbf e_\lambda\)` | `\(0\)` | `\(0\)` | `\(0\)` | `\(-2-2E_*\)` |
| `\(\mathbf e_c\)` | `\(-\mathfrak a_*\)` | `\(-\mathfrak a_*\)` | `\(2\)` | `\(2-F_*\mathfrak a_*\)` |
| `\(\mathbf e_\gamma\)` | `\(-\mathfrak a_*\)` | `\(-\mathfrak a_*\)` | `\(0\)` | `\(-2E_*-F_*\mathfrak a_*\)` |
| `\(\mathbf e_U\)` | `\(+\mathfrak a_*\)` | `\(1+\mathfrak a_*\)` | `\(-1\)` | `\(-1+E_*+F_*\mathfrak a_*\)` |
| `\(\mathbf e_W\)` | `\(0\)` | `\(0\)` | `\(0\)` | `\(2+E_*\)` |

Two immediate consequences are worth recording.

1. The `\(\lambda_W\)` and `\(K_W^{(\mathrm{eff})}\)` rays do **not** move
   \(\delta_{U,*}^{\rm graph}\), \(T_U^{\rm graph}\), or
   \(K_{\eta,*}^{\rm graph}\) at all; they act only through `\(\mu_W^{\rm graph}\)`.
2. The direct tracking quantity `\(\delta_{U,*}^{\rm graph}\)` is controlled only by the combination
   \(s_\gamma+s_c-s_U\).

So this table is the first exact primitive audit library for the free-quintuple branch search.

---

## 6. Scalarized closure function on a log-ray

Define the scalarized graph-ray closure function
\[
\boxed{
\Phi_{\mathbf s}(\tau)
:=
\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau)),
}
\]
and the scalar residual
\[
\boxed{
\Delta_{\mathbf s}(\tau):=\Phi_{\mathbf s}(\tau)-1.
}
\]
Then Stage 186 immediately becomes
\[
\boxed{
\mathbf x_{\mathbf s}^{\rm graph}(\tau)\in\mathcal Z_*
\iff
\Delta_{\mathbf s}(\tau)=0.
}
\]
So once a base point and a free direction are chosen, the whole reduced closure search is a one-variable root problem in `\(\tau\)`.

Introduce the exact directional derivative operator on free-quintuple log space:
\[
\boxed{
\mathcal D_{\mathbf s}
:=
s_\lambda\partial_{\ln\lambda_W}
+s_c\partial_{\ln c_{\eta U}}
+s_\gamma\partial_{\ln\gamma}
+s_U\partial_{\ln K_U}
+s_W\partial_{\ln K_W^{(\mathrm{eff})}}.
}
\]
Then
\[
\boxed{
\Phi_{\mathbf s}'(0)
=
(\mathcal D_{\mathbf s}\widehat\chi_Q)(\mathbf y_\circ),
}
\]
and, whenever `\(\Phi_{\mathbf s}(0)>0\)`,
\[
\boxed{
L_{\mathbf s}(0)
:=
\left.\frac{d}{d\tau}\ln\Phi_{\mathbf s}(\tau)\right|_{\tau=0}
=
(\mathcal D_{\mathbf s}\ln\widehat\chi_Q)(\mathbf y_\circ)
=
\frac{\Phi_{\mathbf s}'(0)}{\Phi_{\mathbf s}(0)}.
}
\]

---

## 7. Monotone-ray uniqueness theorem

Assume `\(\Phi_{\mathbf s}(\tau)>0\)` on an interval `[\tau_-,\tau_+]`.
If the logarithmic slope
\[
L_{\mathbf s}(\tau)=\frac{d}{d\tau}\ln\Phi_{\mathbf s}(\tau)
\]
has a fixed sign on that interval, then `\(\Phi_{\mathbf s}(\tau)\)` is monotone there.
So:

\[
\boxed{\textbf{Theorem (Stage 187 monotone-ray uniqueness theorem).}}
\]

If
\[
\Phi_{\mathbf s}(\tau_-)<1<\Phi_{\mathbf s}(\tau_+)
\qquad\text{or}\qquad
\Phi_{\mathbf s}(\tau_+)<1<\Phi_{\mathbf s}(\tau_-),
\]
and `\(L_{\mathbf s}(\tau)\)` has fixed sign on `[\tau_-,\tau_+]`, then there exists a **unique**
\[
\boxed{\tau_*\in(\tau_-,\tau_+)}
\]
with
\[
\boxed{\Phi_{\mathbf s}(\tau_*)=1.}
\]

So the scalar graph-slice theorem of Stage 186 becomes a practical uniqueness test as soon as one can certify monotonicity along the chosen free-quintuple ray.

---

## 8. Exact local root predictors

Let
\[
\Phi_0:=\Phi_{\mathbf s}(0),
\qquad
\Phi_1:=\Phi_{\mathbf s}'(0),
\qquad
L_0:=\frac{\Phi_1}{\Phi_0}.
\]
Assume `\(\Phi_1\neq0\)`.
Then the exact affine root predictor is
\[
\boxed{
\tau_{\rm aff}:=\frac{1-\Phi_0}{\Phi_1}.
}
\]
If `\(\Phi_0>0\)` and `\(L_0\neq0\)`, the exact log-linear predictor is
\[
\boxed{
\tau_{\log}:=-\frac{\ln\Phi_0}{L_0}.
}
\]

These are the two natural first local predictors for the scalarized branch closure point.

### 8.1 First-order agreement of the two predictors

Write
\[
\Phi_0=1+\varepsilon,
\qquad |\varepsilon|\ll1.
\]
Then
\[
\boxed{
\tau_{\log}-\tau_{\rm aff}
=
-\frac{\varepsilon^2}{2L_0}+O(\varepsilon^3).
}
\]
So both predictors agree to first order in the closure defect and differ only at quadratic order.

This is useful operationally:

- `\(\tau_{\rm aff}\)` is the direct linear predictor for `\(\widehat\chi_Q\)`,
- `\(\tau_{\log}\)` is the direct linear predictor for `\(\ln\widehat\chi_Q\)`.

Near the scalar slice, they are equivalent to first order.

---

## 9. Best current reading after Stage 187

The Stage-186 scalar graph-slice theorem is now fully operational.

The reduced closure search no longer needs a vague “family” in the free quintuple.
It now has an explicit audit-ready one-parameter generator:
\[
\mathbf y_{\mathbf s}(\tau)=\mathbf y_\circ\odot e^{\tau\mathbf s}.
\]
Its graph lift is exact, its monomial invariants are exact, its Packet-B quotient defect vanishes identically, and its remaining reduced closure problem is the single scalar root equation
\[
\widehat\chi_Q(\mathbf y_{\mathbf s}(\tau))=1.
\]

So after Stage 187 the next honest step is no longer to invent a family parameterization.
It is to compute, on candidate free-quintuple rays,

1. the scalar `\(\Phi_{\mathbf s}(0)=\widehat\chi_Q(\mathbf y_\circ)\)`,
2. the directional logarithmic slope `\(L_0=(\mathcal D_{\mathbf s}\ln\widehat\chi_Q)(\mathbf y_\circ)\)`,
3. the local predictor `\(\tau_{\log}\)` (or `\(\tau_{\rm aff}\)`),
4. and then verify whether the actual scalarized branch crosses `1` uniquely before leaving the controlled branch.

That is the sharpest explicit one-parameter continuation point reached so far.
