# Moving-Throat PDE — Stage 207: Primitive-Ray Hessian Envelopes, Exact Certified Ray Table, and the Local Elimination Theorem

## Status

**Exact within the carried Stage 204 primitive free-ray family and the Stage 206 certified local bracket / ranking framework, once a local diagonal Hessian envelope of `\(\ln\widehat\chi_Q\)` is supplied on the chosen primitive intervals.**

This stage does **not** introduce a new constitutive law.
It specializes the generic Stage 206 search sieve to the five primitive free-quintuple rays and isolates the smallest actual differential data the completed PDE must return before the primitive local search is fully decided.

---

## Purpose

Stage 204 already gave the exact primitive free-direction table
\[
\mathbf e_\lambda,\ \mathbf e_c,\ \mathbf e_\gamma,\ \mathbf e_U,\ \mathbf e_W,
\]
and Stage 206 upgraded the scalarized branch search to a certified local bracket theorem on an arbitrary oriented ray.

That left one sharp continuation point:

> what does the Stage 206 sieve become when we apply it to the **actual primitive candidate rays**, and how much Hessian data does the PDE really need at that level?

This stage answers that exactly.

The main outputs are:

1. the exact **canonical orientation rule** for each primitive ray,
2. the exact theorem that a primitive ray uses only one diagonal Hessian restriction
   \[
   \partial_{\ell_i\ell_i}\ln\widehat\chi_Q,
   \]
   with no off-diagonal entry,
3. the exact **primitive certified ray table** carrying both the Stage 204 microscopic drift data and the Stage 206 bracket data,
4. the exact **primitive elimination theorem**,
5. the exact **primitive winner theorem** when one certified bracket sits strictly to the left of all others,
6. and the exact statement that off-diagonal Hessian entries first appear only for genuine two-coordinate mixed rays.

So Stage 207 is the first point where the final scalar branch search becomes a concrete five-row primitive audit table rather than a generic ray library.

---

## 1. Carry-forward free log coordinates and the scalar closure function

Keep the free-quintuple logarithmic coordinates
\[
\boxed{
\boldsymbol\ell:=(\ell_\lambda,\ell_c,\ell_\gamma,\ell_U,\ell_W)
=(\ln\lambda_W,\ln c_{\eta U},\ln\gamma,\ln K_U,\ln K_W^{(\mathrm{eff})}).
}
\]
Define the carried scalar closure function
\[
\boxed{
h(\boldsymbol\ell):=\ln\widehat\chi_Q(\boldsymbol\ell).
}
\]
Fix a positive base point `\(\boldsymbol\ell_\circ\)` and write
\[
\boxed{
h_0:=h(\boldsymbol\ell_\circ),
\qquad
\sigma_0:=\operatorname{sgn}(h_0),
\qquad
H_0:=|h_0|>0.
}
\]
The oriented logarithmic residual of Stage 206 is therefore
\[
H_{\mathbf s}(\tau)=\sigma_0\,h(\boldsymbol\ell_\circ+\tau\mathbf s).
\]
Its initial directional slope and curvature are
\[
K_0(\mathbf s)=\sigma_0\,(\mathbf s\cdot\nabla_\ell h)(\boldsymbol\ell_\circ),
\qquad
H_1(\mathbf s;\tau)=\sigma_0\,\mathbf s^T H_h(\boldsymbol\ell_\circ+\tau\mathbf s)\mathbf s,
\]
where `\(H_h\)` is the Hessian of `\(h\)` in the log coordinates.

So the Stage 206 data for any candidate ray are still:

1. the positive initial defect `\(H_0\)`,
2. the oriented initial slope `\(K_0\)`,
3. a local curvature envelope for `\(H_1\)`.

---

## 2. Primitive candidate rays and their canonical orientation

The five primitive free rays are
\[
\boxed{
\mathbf e_\lambda,
\quad
\mathbf e_c,
\quad
\mathbf e_\gamma,
\quad
\mathbf e_U,
\quad
\mathbf e_W.
}
\]
Define the five oriented gradient components at the base point
\[
\boxed{
\Gamma_i:=\sigma_0\,\partial_{\ell_i}h(\boldsymbol\ell_\circ),
\qquad
i\in\{\lambda,c,\gamma,U,W\}.
}
\]
These are the Stage 206 initial slopes on the **positive** primitive axes.

### 2.1 Monotone primitive rays

If `\(\Gamma_i\neq0\)`, define the canonical sign
\[
\boxed{
\varepsilon_i:=-\operatorname{sgn}(\Gamma_i),
\qquad
\widehat{\mathbf s}_i:=\varepsilon_i\mathbf e_i.
}
\]
Then the oriented initial slope is exactly
\[
\boxed{
K_i:=K_0(\widehat{\mathbf s}_i)=\varepsilon_i\Gamma_i=-|\Gamma_i|<0.
}
\]
So every nonflat primitive axis admits a unique canonical forward orientation that points toward local closure.

### 2.2 Turning primitive rays

If `\(\Gamma_i=0\)`, the primitive axis is a turning candidate. Its sign is irrelevant at first order, and Stage 206 says the forward root question is decided entirely by the sign and envelope of the second derivative.

So the primitive search problem splits cleanly:

- `\(\Gamma_i\neq0\)`  → monotone primitive candidate,
- `\(\Gamma_i=0\)`     → turning primitive candidate.

---

## 3. Exact primitive Hessian-envelope theorem

For a primitive candidate ray,
\[
\widehat{\mathbf s}_i=\varepsilon_i\mathbf e_i,
\qquad
\varepsilon_i^2=1.
\]
Therefore the exact second-order directional operator is
\[
\boxed{
\mathcal H_{\widehat{\mathbf s}_i}
=
(\varepsilon_i\mathbf e_i)^T H_h (\varepsilon_i\mathbf e_i)
=
\partial_{\ell_i\ell_i}.
}
\]
So the primitive ray curvature depends only on the corresponding diagonal Hessian restriction:
\[
\boxed{
H_1^{(i)}(\tau)
=
\sigma_0\,\partial_{\ell_i\ell_i}h(\boldsymbol\ell_\circ+\tau\widehat{\mathbf s}_i).
}
\]
This is the first main theorem of the stage.

### Primitive Hessian-envelope theorem

If, on an interval `\([0,T_i]\)`, the oriented diagonal Hessian restriction obeys
\[
\boxed{
\underline\kappa_i
\le
\sigma_0\,\partial_{\ell_i\ell_i}h(\boldsymbol\ell_\circ+\tau\widehat{\mathbf s}_i)
\le
\overline\kappa_i,
\qquad
0\le\tau\le T_i,
}
\]
then that interval is the **full** Stage 206 curvature envelope for the primitive ray.

So no off-diagonal Hessian entry is needed to certify any primitive ray.
The completed PDE must provide only the five diagonal restrictions
\[
\partial_{\ell_\lambda\ell_\lambda}h,
\quad
\partial_{\ell_c\ell_c}h,
\quad
\partial_{\ell_\gamma\ell_\gamma}h,
\quad
\partial_{\ell_U\ell_U}h,
\quad
\partial_{\ell_W\ell_W}h
\]
along the corresponding oriented primitive intervals.

---

## 4. Sign-adapted primitive microscopic drift table

Stage 204 already fixed the dependent graph exponents in terms of the free log-ray direction.
Write
\[
\boxed{
\mathfrak a_*:=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.
}
\]
For a primitive ray `\(\widehat{\mathbf s}_i=\varepsilon_i\mathbf e_i\)`, the dependent graph exponents are just the Stage 204 primitive values multiplied by the orientation sign `\(\varepsilon_i\)`.

| Primitive ray | `\(\sigma_\delta\)` | `\(\sigma_T\)` | `\(\sigma_{K_\eta}\)` | `\(\sigma_\mu\)` |
|---|---:|---:|---:|---:|
| `\(\widehat{\mathbf s}_\lambda=\varepsilon_\lambda\mathbf e_\lambda\)` | `\(0\)` | `\(0\)` | `\(0\)` | `\(\varepsilon_\lambda(-2-2E_*)\)` |
| `\(\widehat{\mathbf s}_c=\varepsilon_c\mathbf e_c\)` | `\(-\varepsilon_c\mathfrak a_*\)` | `\(-\varepsilon_c\mathfrak a_*\)` | `\(2\varepsilon_c\)` | `\(\varepsilon_c(2-F_*\mathfrak a_*)\)` |
| `\(\widehat{\mathbf s}_\gamma=\varepsilon_\gamma\mathbf e_\gamma\)` | `\(-\varepsilon_\gamma\mathfrak a_*\)` | `\(-\varepsilon_\gamma\mathfrak a_*\)` | `\(0\)` | `\(\varepsilon_\gamma(-2E_*-F_*\mathfrak a_*)\)` |
| `\(\widehat{\mathbf s}_U=\varepsilon_U\mathbf e_U\)` | `\(+\varepsilon_U\mathfrak a_*\)` | `\(\varepsilon_U(1+\mathfrak a_*)\)` | `\(-\varepsilon_U\)` | `\(\varepsilon_U(-1+E_*+F_*\mathfrak a_*)\)` |
| `\(\widehat{\mathbf s}_W=\varepsilon_W\mathbf e_W\)` | `\(0\)` | `\(0\)` | `\(0\)` | `\(\varepsilon_W(2+E_*)\)` |

So Stage 207 does not replace the Stage 204 primitive family compiler.
It adds the exact certified search data to it.

---

## 5. Exact certified primitive ray table

For each primitive index `\(i\in\{\lambda,c,\gamma,U,W\}\)`, define the primitive data packet
\[
\boxed{
\mathcal R_i^{\rm prim}
:=
(H_0,\Gamma_i,\underline\kappa_i,\overline\kappa_i,T_i).
}
\]
The Stage 206 root map is
\[
\mathcal T(H_0,k;c)= -\frac{2H_0}{-k-\sqrt{k^2-2cH_0}} = \frac{2H_0}{k+\sqrt{k^2-2cH_0}},
\qquad k>0.
\]

### 5.1 Monotone primitive rows

If `\(\Gamma_i\neq0\)`, set
\[
\boxed{k_i:=|\Gamma_i|>0.}
\]
The two discriminants are
\[
\boxed{
\Delta_i^{\rm lo}:=k_i^2-2\underline\kappa_iH_0,
\qquad
\Delta_i^{\rm hi}:=k_i^2-2\overline\kappa_iH_0.
}
\]
If both are nonnegative, define the certified bracket
\[
\boxed{
\tau_{i,\rm lo}:=\mathcal T(H_0,k_i;\underline\kappa_i),
\qquad
\tau_{i,\rm hi}:=\mathcal T(H_0,k_i;\overline\kappa_i),
}
\]
with width
\[
\boxed{W_i:=\tau_{i,\rm hi}-\tau_{i,\rm lo}.}
\]
If also `\(\tau_{i,\rm hi}\le T_i\)`, then the primitive ray is locally admissible and Stage 206 guarantees one unique true primitive closure point
\[
\boxed{\tau_*^{(i)}\in[\tau_{i,\rm lo},\tau_{i,\rm hi}].}
\]

### 5.2 Turning primitive rows

If `\(\Gamma_i=0\)`, the primitive ray is turning-admissible only when
\[
\boxed{\overline\kappa_i<0.}
\]
Then the exact turning comparison roots are
\[
\boxed{
\tau_{i,\rm lo}^{\rm(tp)}:=\sqrt{-\frac{2H_0}{\underline\kappa_i}},
\qquad
\tau_{i,\rm hi}^{\rm(tp)}:=\sqrt{-\frac{2H_0}{\overline\kappa_i}},
}
\]
with turning width
\[
\boxed{W_i^{\rm(tp)}:=\tau_{i,\rm hi}^{\rm(tp)}-\tau_{i,\rm lo}^{\rm(tp)}.}
\]
If `\(\tau_{i,\rm hi}^{\rm(tp)}\le T_i\)`, then Stage 206 again gives one unique primitive closure point on that axis.

### 5.3 Compact certified primitive-ray table

So each primitive row is decided by the exact data

| Row type | Input data | Certified output |
|---|---|---|
| monotone primitive | `\(H_0,\ k_i,\ \underline\kappa_i,\ \overline\kappa_i,\ T_i\)` | `\([\tau_{i,\rm lo},\tau_{i,\rm hi}],\ W_i\)` |
| turning primitive | `\(H_0,\ \underline\kappa_i,\ \overline\kappa_i,\ T_i\)` with `\(\overline\kappa_i<0\)` | `\([\tau_{i,\rm lo}^{\rm(tp)},\tau_{i,\rm hi}^{\rm(tp)}],\ W_i^{\rm(tp)}\)` |

This is the first fully explicit certified ray table for the moving-throat free-quintuple search.

---

## 6. Exact primitive elimination theorem

A primitive candidate ray is **eliminated** from the certified local search sieve if any of the following exact conditions hold.

### 6.1 Monotone primitive elimination

For `\(\Gamma_i\neq0\)`, the primitive row is eliminated if

1. one of the discriminants is negative,
   \[
   \Delta_i^{\rm lo}<0\quad\text{or}\quad\Delta_i^{\rm hi}<0,
   \]
   so no certified local quadratic comparison root exists; or
2. the certified upper bracket exceeds the local validity interval,
   \[
   \tau_{i,\rm hi}>T_i.
   \]

### 6.2 Turning primitive elimination

For `\(\Gamma_i=0\)`, the primitive row is eliminated if

1. `\(\overline\kappa_i\ge0\)`, so the turning branch is not locally downward-bending; or
2. `\(\tau_{i,\rm hi}^{\rm(tp)}>T_i\)`.

### 6.3 Primitive elimination theorem

Therefore:
\[
\boxed{\textbf{A primitive ray survives the certified local sieve iff it satisfies the corresponding Stage 206 admissibility test on its own primitive row.}}
\]

So after Stage 207 the primitive local search is no longer heuristic at all.
It is a finite exact elimination problem on five rows.

---

## 7. Exact primitive winner theorem and the local priority pair

Suppose primitive rows `\(i\)` and `\(j\)` are both admissible.
Then Stage 206 gives the exact pairwise ordering theorem:
\[
\boxed{
\tau_{i,\rm hi}<\tau_{j,\rm lo}
\quad\Longrightarrow\quad
\tau_*^{(i)}<\tau_*^{(j)}.
}
\]
So if there exists one primitive row `\(i\)` with
\[
\boxed{
\tau_{i,\rm hi}<\min_{j\neq i}\tau_{j,\rm lo},
}
\]
then `\(i\)` is the unique earliest primitive closure direction.

When certified brackets overlap, the exact theorem-level winner is no longer fixed, but the Stage 206 local priority pair remains
\[
\boxed{
\mathcal Q_i:=(\tau_{i,\rm hi},W_i)
}
\]
(or the turning analogue when `\(\Gamma_i=0\)`).
So the primitive local ranking rule is:

1. eliminate every nonadmissible primitive row,
2. certify strict winner status whenever one upper bracket lies below all other lower brackets,
3. among overlapping surviving rows, sort lexicographically by `\(\mathcal Q_i\)`.

This is the exact Stage 207 primitive winner theorem.

---

## 8. Exact diagonal-Hessian sufficiency theorem and the mixed-ray preview

The primitive Hessian-envelope theorem shows that primitive certification uses only the five diagonal Hessian restrictions.
The reason is exact:
\[
\boxed{
\mathcal H_{\widehat{\mathbf s}_i}=\partial_{\ell_i\ell_i}.
}
\]
But if we move to a genuine two-coordinate ray
\[
\mathbf s=a\mathbf e_i+b\mathbf e_j,
\qquad i\neq j,
\]
then the exact second directional derivative becomes
\[
\boxed{
\mathcal H_{\mathbf s}
=
a^2\partial_{\ell_i\ell_i}
+2ab\partial_{\ell_i\ell_j}
+b^2\partial_{\ell_j\ell_j}.
}
\]
So the off-diagonal Hessian entries `\(\partial_{\ell_i\ell_j}h\)` first appear **only** for mixed rays.

This gives the exact sufficiency theorem:
\[
\boxed{
\textbf{Stage 207 exhausts all certified primitive-ray information using only the diagonal Hessian envelopes.}
}
\]
No cross-Hessian datum can change any primitive certified bracket or primitive winner ordering.
It can only matter when the search is enlarged beyond the primitive family.

---

## 9. Best current reading after Stage 207

Stage 204 gave the exact primitive free-direction family.
Stage 206 turned the scalarized branch search into a certified local bracket / ranking theorem.
Stage 207 now compresses the primitive local search to the smallest exact data packet the completed PDE still has to provide:

1. the base defect magnitude `\(H_0=|\ln\widehat\chi_Q(\boldsymbol\ell_\circ)|\)`,
2. the five oriented primitive gradient components `\(\Gamma_i\)`,
3. the five oriented primitive diagonal Hessian envelopes `\([\underline\kappa_i,\overline\kappa_i]\)`,
4. the five local validity lengths `\(T_i\)`.

Those data determine the entire primitive certified ray table, the primitive eliminations, and—when brackets separate—the exact primitive winner.

So the next honest continuation is no longer to discuss primitive rays abstractly.
It is to move one step beyond them:

> evaluate the first genuine two-coordinate mixed rays, where the off-diagonal Hessian entries of `\(\ln\widehat\chi_Q\)` can finally enter and possibly beat the primitive certified brackets.

That is the natural Stage 208 gate.
