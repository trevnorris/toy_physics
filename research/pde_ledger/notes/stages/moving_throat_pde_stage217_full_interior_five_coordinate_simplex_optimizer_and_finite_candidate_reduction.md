# Moving-Throat PDE — Stage 217: Unique Five-Coordinate Mixed Simplex, Exact Face Reduction, Full Interior Optimizer, and the Finite Algebraic Candidate Set

## Status

**Exact within the carried Stage 249 support-`<=4` certified ledger, once a compact interior ratio window and the corresponding validity map are supplied on the chosen five-coordinate patch.**

Stage 250 now supplies the support-cardinality-`5` gate and the exact five-face boundary splice. Because there is only **one** primitive five-coordinate simplex,
\[
\{\lambda,c,\gamma,U,W\},
\]
this note can focus entirely on the only genuinely new remaining content: the **four-parameter interior** problem.

---

## Purpose

Stage 249 finished the full support-`<=4` search, and Stage 250 then built the unique positive spherical five-simplex, reduced its full codimension-one boundary back to the Stage 249 quadruple packets, and isolated the two canonical five-way interior screens.

That leaves one exact continuation point:

> can the **unique interior five-coordinate** mixed branch beat the already-finished support-`<=4` ledger, and can that interior optimization again be reduced to a finite algebraic candidate set rather than a free four-parameter continuum scan?

This stage answers that.

The main outputs are:

1. the exact positive spherical five-simplex and its five-face reduction back to the Stage 249 quadruple packets,
2. the exact four-ratio interior chart and certified root functional,
3. the exact four-component stationary numerator system,
4. the exact lifted polynomial stationary system with degree pattern `(3,3,3,3,2)`,
5. the exact finite algebraic candidate-set theorem with lifted Bézout bound `162`,
6. the direct square-root-free elimination into three quintic cross-consistency polynomials plus one sextic square condition on a chosen interior chart, giving the projected bound `750`,
7. two exact special reductions explaining why the gradient-optimal and equal-mix interior screens are the right canonical screens,
8. and the exact local improvement / non-improvement theorems against the imported Stage 249 five-face boundary ledger.

So Stage 251 turns the last genuinely new local mixed-simplex problem into a finite candidate problem.

---

## 1. The unique positive spherical five-simplex and exact face reduction

Keep the oriented primitive free axes
\[
\widehat{\mathbf e}_\lambda,\quad
\widehat{\mathbf e}_c,\quad
\widehat{\mathbf e}_\gamma,\quad
\widehat{\mathbf e}_U,\quad
\widehat{\mathbf e}_W,
\]
with positive oriented slope magnitudes
\[
k_\lambda,\quad k_c,\quad k_\gamma,\quad k_U,\quad k_W>0.
\]

The unique support-cardinality-`5` search set is the positive spherical five-simplex
\[
\boxed{
\Delta_5^+
:=
\Bigl\{
\mathbf a=(a_\lambda,a_c,a_\gamma,a_U,a_W)\in\mathbb R_{\ge 0}^5:
a_\lambda^2+a_c^2+a_\gamma^2+a_U^2+a_W^2=1
\Bigr\}.
}
\]
Each simplex point generates the oriented mixed ray
\[
\boxed{
\widehat{\mathbf s}_5(\mathbf a)
=
a_\lambda\widehat{\mathbf e}_\lambda
+a_c\widehat{\mathbf e}_c
+a_\gamma\widehat{\mathbf e}_\gamma
+a_U\widehat{\mathbf e}_U
+a_W\widehat{\mathbf e}_W.
}
\]

Its five codimension-one faces are exactly the five Stage 249 primitive quadruple simplices:
\[
Q_{\widehat\lambda}=\{c,\gamma,U,W\},
\qquad
Q_{\widehat c}=\{\lambda,\gamma,U,W\},
\]
\[
Q_{\widehat\gamma}=\{\lambda,c,U,W\},
\qquad
Q_{\widehat U}=\{\lambda,c,\gamma,W\},
\qquad
Q_{\widehat W}=\{\lambda,c,\gamma,U\}.
\]

So the entire codimension-one boundary is already closed by the Stage 249 ledger.

Let the imported Stage 249 closed-simplex intervals be
\[
\mathcal I_{\widehat\lambda}^{\square},
\quad
\mathcal I_{\widehat c}^{\square},
\quad
\mathcal I_{\widehat\gamma}^{\square},
\quad
\mathcal I_{\widehat U}^{\square},
\quad
\mathcal I_{\widehat W}^{\square},
\]
with lower and upper face minima
\[
\boxed{
\beta_5^{\rm lo}
:=
\min\bigl(
\tau_{\widehat\lambda,\min}^{\rm lo,\square},
\tau_{\widehat c,\min}^{\rm lo,\square},
\tau_{\widehat\gamma,\min}^{\rm lo,\square},
\tau_{\widehat U,\min}^{\rm lo,\square},
\tau_{\widehat W,\min}^{\rm lo,\square}
\bigr),
}
\]
\[
\boxed{
\beta_5^{\rm hi}
:=
\min\bigl(
\tau_{\widehat\lambda,\min}^{\rm hi,\square},
\tau_{\widehat c,\min}^{\rm hi,\square},
\tau_{\widehat\gamma,\min}^{\rm hi,\square},
\tau_{\widehat U,\min}^{\rm hi,\square},
\tau_{\widehat W,\min}^{\rm hi,\square}
\bigr).
}
\]

So from this point onward, the only genuinely new content is the **interior** of `\(\Delta_5^+\)`.

---

## 2. Interior ratio chart and the exact certified objective

On the interior chart `a_\lambda>0`, define the positive ratio coordinates
\[
r:=\frac{a_c}{a_\lambda}>0,
\qquad
s:=\frac{a_\gamma}{a_\lambda}>0,
\qquad
t:=\frac{a_U}{a_\lambda}>0,
\qquad
u:=\frac{a_W}{a_\lambda}>0.
\]
Then
\[
\boxed{
\mathbf a(r,s,t,u)
=
\frac{(1,r,s,t,u)}{\sqrt{1+r^2+s^2+t^2+u^2}}.
}
\]
The oriented initial slope magnitude is
\[
\boxed{
k_5(r,s,t,u)
=
\frac{k_\lambda+r k_c+s k_\gamma+t k_U+u k_W}
{\sqrt{1+r^2+s^2+t^2+u^2}}.
}
\]

For either envelope label
\[
\star\in\{{\rm lo},{\rm hi}\},
\]
write the `5\times 5` symmetric Hessian-envelope block entries as
\[
\nu_{\lambda\lambda,\star},\ \nu_{\lambda c,\star},\ \nu_{\lambda\gamma,\star},\ \nu_{\lambda U,\star},\ \nu_{\lambda W,\star},
\]
\[
\nu_{cc,\star},\ \nu_{c\gamma,\star},\ \nu_{cU,\star},\ \nu_{cW,\star},
\]
\[
\nu_{\gamma\gamma,\star},\ \nu_{\gamma U,\star},\ \nu_{\gamma W,\star},
\]
\[
\nu_{UU,\star},\ \nu_{UW,\star},\ \nu_{WW,\star}.
\]

Define the exact discriminant coefficients
\[
A_\star:=k_\lambda^2-2H_0\nu_{\lambda\lambda,\star},
\qquad
B_\star:=2k_\lambda k_c-4H_0\nu_{\lambda c,\star},
\]
\[
C_\star:=2k_\lambda k_\gamma-4H_0\nu_{\lambda\gamma,\star},
\qquad
D_\star:=2k_\lambda k_U-4H_0\nu_{\lambda U,\star},
\qquad
E_\star:=2k_\lambda k_W-4H_0\nu_{\lambda W,\star},
\]
\[
F_\star:=k_c^2-2H_0\nu_{cc,\star},
\qquad
G_\star:=2k_c k_\gamma-4H_0\nu_{c\gamma,\star},
\qquad
H_\star:=2k_c k_U-4H_0\nu_{cU,\star},
\qquad
I_\star:=2k_c k_W-4H_0\nu_{cW,\star},
\]
\[
J_\star:=k_\gamma^2-2H_0\nu_{\gamma\gamma,\star},
\qquad
K_\star:=2k_\gamma k_U-4H_0\nu_{\gamma U,\star},
\qquad
L_\star:=2k_\gamma k_W-4H_0\nu_{\gamma W,\star},
\]
\[
M_\star:=k_U^2-2H_0\nu_{UU,\star},
\qquad
N_\star:=2k_U k_W-4H_0\nu_{UW,\star},
\qquad
O_\star:=k_W^2-2H_0\nu_{WW,\star}.
\]

Then the exact interior discriminant numerator is
\[
\boxed{
\Delta^\sharp_\star(r,s,t,u)
=
A_\star+B_\star r+C_\star s+D_\star t+E_\star u
+F_\star r^2+G_\star rs+H_\star rt+I_\star ru
+J_\star s^2+K_\star st+L_\star su
+M_\star t^2+N_\star tu+O_\star u^2.
}
\]

The exact certified root function is
\[
\boxed{
\tau_\star(r,s,t,u)
=
\frac{2H_0\sqrt{1+r^2+s^2+t^2+u^2}}
{k_\lambda+r k_c+s k_\gamma+t k_U+u k_W+\sqrt{\Delta^\sharp_\star(r,s,t,u)}}.
}
\]
Equivalently define
\[
\boxed{
\Phi_\star(r,s,t,u)
:=
\frac{k_\lambda+r k_c+s k_\gamma+t k_U+u k_W+\sqrt{\Delta^\sharp_\star(r,s,t,u)}}
{\sqrt{1+r^2+s^2+t^2+u^2}},
\qquad
\tau_\star=\frac{2H_0}{\Phi_\star}.
}
\]

Let the compact interior ratio window be
\[
\boxed{
\mathcal W_5:=[0,R_5]\times[0,S_5]\times[0,T_5]\times[0,U_5],
\qquad
0<R_5,S_5,T_5,U_5<\infty,
}
\]
with validity map `\(\mathcal T_5(r,s,t,u)\)`. The admissible interior set is
\[
\boxed{
\mathcal A_{5,\star}^{\rm int}
:=
\Bigl\{
(r,s,t,u)\in(0,\infty)^4\cap \mathcal W_5:
\Delta^\sharp_\star\ge 0,\ 
\tau_\star(r,s,t,u)\le \mathcal T_5(r,s,t,u)
\Bigr\}.
}
\]

---

## 3. Exact four-component stationary numerator theorem

Write
\[
S:=1+r^2+s^2+t^2+u^2,
\qquad
K_{\rm lin}:=k_\lambda+r k_c+s k_\gamma+t k_U+u k_W.
\]

Define the exact slope numerators
\[
\boxed{
M_r:=S k_c-r K_{\rm lin}
=
k_c(1+s^2+t^2+u^2)-r(k_\lambda+s k_\gamma+t k_U+u k_W),
}
\]
\[
\boxed{
M_s:=S k_\gamma-s K_{\rm lin}
=
k_\gamma(1+r^2+t^2+u^2)-s(k_\lambda+r k_c+t k_U+u k_W),
}
\]
\[
\boxed{
M_t:=S k_U-t K_{\rm lin}
=
k_U(1+r^2+s^2+u^2)-t(k_\lambda+r k_c+s k_\gamma+u k_W),
}
\]
\[
\boxed{
M_u:=S k_W-u K_{\rm lin}
=
k_W(1+r^2+s^2+t^2)-u(k_\lambda+r k_c+s k_\gamma+t k_U).
}
\]

Also define the discriminant-transport numerators
\[
\boxed{
L_{r,\star}:=S\,\partial_r\Delta^\sharp_\star-2r\Delta^\sharp_\star,
}
\qquad
\boxed{
L_{s,\star}:=S\,\partial_s\Delta^\sharp_\star-2s\Delta^\sharp_\star,
}
\]
\[
\boxed{
L_{t,\star}:=S\,\partial_t\Delta^\sharp_\star-2t\Delta^\sharp_\star,
}
\qquad
\boxed{
L_{u,\star}:=S\,\partial_u\Delta^\sharp_\star-2u\Delta^\sharp_\star.
}
\]

Then the exact stationary numerators are
\[
\boxed{
\mathcal N_{r,\star}=2M_r\sqrt{\Delta^\sharp_\star}+L_{r,\star},
\qquad
\mathcal N_{s,\star}=2M_s\sqrt{\Delta^\sharp_\star}+L_{s,\star},
}
\]
\[
\boxed{
\mathcal N_{t,\star}=2M_t\sqrt{\Delta^\sharp_\star}+L_{t,\star},
\qquad
\mathcal N_{u,\star}=2M_u\sqrt{\Delta^\sharp_\star}+L_{u,\star}.
}
\]

### Exact stationary numerator theorem

For every interior admissible point with `\(\Delta^\sharp_\star>0\)`,
\[
\boxed{
\partial_r\Phi_\star=0,\quad
\partial_s\Phi_\star=0,\quad
\partial_t\Phi_\star=0,\quad
\partial_u\Phi_\star=0
\iff
\mathcal N_{r,\star}=\mathcal N_{s,\star}=\mathcal N_{t,\star}=\mathcal N_{u,\star}=0.
}
\]
So every interior optimizer is an admissible common zero of four exact stationary numerators.

---

## 4. Exact lifted polynomial stationary system and finite candidate set

Introduce the auxiliary square-root variable
\[
y:=\sqrt{\Delta^\sharp_\star(r,s,t,u)}\ge 0.
\]
Then the lifted stationary system is
\[
\boxed{
\mathcal F_{r,\star}:=2M_r\,y+L_{r,\star}=0,
\qquad
\mathcal F_{s,\star}:=2M_s\,y+L_{s,\star}=0,
}
\]
\[
\boxed{
\mathcal F_{t,\star}:=2M_t\,y+L_{t,\star}=0,
\qquad
\mathcal F_{u,\star}:=2M_u\,y+L_{u,\star}=0,
}
\]
together with
\[
\boxed{
\mathcal F_{\Delta,\star}:=y^2-\Delta^\sharp_\star(r,s,t,u)=0.
}
\]

### Exact degree ledger

Each stationary equation `\(\mathcal F_{r,\star},\mathcal F_{s,\star},\mathcal F_{t,\star},\mathcal F_{u,\star}\)` has total degree `3` in `(r,s,t,u,y)`, while `\(\mathcal F_{\Delta,\star}\)` has total degree `2`. So the exact lifted degree pattern is
\[
\boxed{(3,3,3,3,2).}
\]

Hence the lifted Bézout candidate bound is
\[
\boxed{3\cdot 3\cdot 3\cdot 3\cdot 2 = 162.}
\]

Define the lifted admissible stationary set
\[
\boxed{
\mathcal C_{5,\star}^{\rm int,lift}
:=
\Bigl\{
(r,s,t,u,y)\in \mathcal W_5\times\mathbb R_{\ge 0}:
\mathcal F_{r,\star}=\mathcal F_{s,\star}=\mathcal F_{t,\star}=\mathcal F_{u,\star}=\mathcal F_{\Delta,\star}=0,
\ \tau_\star(r,s,t,u)\le \mathcal T_5(r,s,t,u)
\Bigr\}.
}
\]

### Exact lifted finite candidate-set theorem

Assume the optimizer does not sit on an artificial outer boundary of the chosen ratio window and that the lifted stationary set is zero-dimensional on that window. Then every interior optimizer of `\(\tau_\star\)` belongs to the finite admissible set `\(\mathcal C_{5,\star}^{\rm int,lift}\)`.

So the interior optimized brackets are finite evaluations:
\[
\boxed{
\tau_{5,\min}^{\star,\rm int}
=
\min_{(r,s,t,u,y)\in\mathcal C_{5,\star}^{\rm int,lift}}
\tau_\star(r,s,t,u),
\qquad
\star\in\{{\rm lo},{\rm hi}\}.
}
\]

This is the preferred exact candidate compiler for the unique support-cardinality-`5` interior problem.

---

## 5. Direct square-root-free elimination and the projected candidate bound

The lifted system is the cleanest compiler, but the direct elimination remains useful.

### 5.1 Quintic cross-consistency polynomials

From the four stationary equations,
\[
2M_r y + L_{r,\star}=0,\quad
2M_s y + L_{s,\star}=0,\quad
2M_t y + L_{t,\star}=0,\quad
2M_u y + L_{u,\star}=0,
\]
eliminate `\(y\)` pairwise and define
\[
\boxed{
\mathcal C_{rs,\star}:=M_sL_{r,\star}-M_rL_{s,\star},
}
\qquad
\boxed{
\mathcal C_{rt,\star}:=M_tL_{r,\star}-M_rL_{t,\star},
}
\qquad
\boxed{
\mathcal C_{ru,\star}:=M_uL_{r,\star}-M_rL_{u,\star}.
}
\]
Each is **quintic** in `(r,s,t,u)`.

### 5.2 Sextic square condition

Squaring the `\(r\)`-stationary numerator gives
\[
\boxed{
\mathcal S_{r,\star}:=L_{r,\star}^2-4M_r^2\Delta^\sharp_\star=0,
}
\]
which is **sextic** in `(r,s,t,u)`.

So on a nondegenerate `r`-chart every interior stationary point lies in the projected algebraic pre-candidate set
\[
\boxed{
\widetilde{\mathcal C}_{5,\star;r}^{\rm int}
:=
\Bigl\{
(r,s,t,u)\in\mathcal W_5:
\mathcal C_{rs,\star}=0,\ 
\mathcal C_{rt,\star}=0,\ 
\mathcal C_{ru,\star}=0,\ 
\mathcal S_{r,\star}=0
\Bigr\}.
}
\]

Because the three cross-consistency polynomials are quintic and the square condition is sextic, Bézout gives the projected one-chart bound
\[
\boxed{
5\cdot 5\cdot 5\cdot 6 = 750.
}
\]

This projected bound is useful, but the lifted `162`-point system above remains the preferred exact compiler because it avoids chart-dependent square-root degeneracies.

---

## 6. Two exact special reductions

These justify the canonical interior screens inherited from the skipped gate logic.

### 6.1 Diagonal-isotropic curvature reduction

Suppose the interior Hessian envelope is diagonal-isotropic:
\[
\nu_{\lambda c,\star}=\nu_{\lambda\gamma,\star}=\nu_{\lambda U,\star}=\nu_{\lambda W,\star}=0,
\]
\[
\nu_{c\gamma,\star}=\nu_{cU,\star}=\nu_{cW,\star}=\nu_{\gamma U,\star}=\nu_{\gamma W,\star}=\nu_{UW,\star}=0,
\]
and
\[
\nu_{\lambda\lambda,\star}=\nu_{cc,\star}=\nu_{\gamma\gamma,\star}=\nu_{UU,\star}=\nu_{WW,\star}=: \nu_\star.
\]
Then the unique interior gradient-optimal ray
\[
\boxed{
\mathbf a_{\rm grad}
=
\frac{(k_\lambda,k_c,k_\gamma,k_U,k_W)}
{\sqrt{k_\lambda^2+k_c^2+k_\gamma^2+k_U^2+k_W^2}}
}
\]
is an exact stationary ray.

In the `a_\lambda>0` ratio chart this is
\[
\boxed{
r_{\rm grad}=\frac{k_c}{k_\lambda},
\qquad
s_{\rm grad}=\frac{k_\gamma}{k_\lambda},
\qquad
t_{\rm grad}=\frac{k_U}{k_\lambda},
\qquad
u_{\rm grad}=\frac{k_W}{k_\lambda}.
}
\]

### 6.2 Full fivefold symmetry reduction

Suppose the free slopes are all equal,
\[
k_\lambda=k_c=k_\gamma=k_U=k_W=:k,
\]
and the Hessian envelope is fully permutation-symmetric over the five primitive axes. Then the equal-mix barycenter
\[
\boxed{
\mathbf a_{\rm eq}=\frac{(1,1,1,1,1)}{\sqrt5}
}
\]
is an exact stationary ray. In the `a_\lambda>0` chart this is simply
\[
\boxed{
r=s=t=u=1.
}
\]

So the two canonical interior screens are not arbitrary. They are the exact special reductions of the full five-coordinate interior problem.

---

## 7. Exact local five-coordinate improvement / non-improvement theorems

Let
\[
\tau_{5,\min}^{\rm lo,int}
\le
\tau_{5,*}^{\rm best,int}
\le
\tau_{5,\min}^{\rm hi,int}
\]
be the exact interior certified interval supplied by the lifted finite candidate compiler.

Because the whole codimension-one boundary is already closed by the imported Stage 249 five-face ledger, the local five-coordinate comparison is now exact.

### Genuine five-coordinate interior improvement theorem

If
\[
\boxed{
\tau_{5,\min}^{\rm hi,int}<\beta_5^{\rm lo},
}
\]
then every admissible interior five-coordinate winner beats every admissible support-`<=4` boundary winner. So the unique support-cardinality-`5` simplex carries a genuine interior improvement.

### No genuine five-coordinate interior improvement theorem

If
\[
\boxed{
\tau_{5,\min}^{\rm lo,int}>\beta_5^{\rm hi},
}
\]
then no admissible interior five-coordinate optimizer beats the already-finished support-`<=4` boundary ledger on the chosen window.

So the support-cardinality-`5` question is now reduced to a finite local interval comparison.

---

## 8. Exact incremental evaluation budget

The lifted interior candidate compiler contributes at most
\[
\boxed{162}
\]
stationary candidates per envelope on the chosen window, hence at most
\[
\boxed{324}
\]
candidate evaluations across the `{lo,hi}` envelopes.

This is the only genuinely new finite search work beyond the already-finished support-`<=4` ledger.

---

## 9. Best current summary after Stage 251

The support-cardinality-`5` problem is no longer a continuum search.

- Its entire codimension-one boundary is already closed by Stage 249.
- Its genuinely new interior optimizer is governed by an exact lifted polynomial system with degree pattern `(3,3,3,3,2)`.
- Every interior optimizer lies in a finite admissible candidate set with lifted Bézout bound `162`.
- And the local improvement/non-improvement verdict against the full support-`<=4` ledger is now an exact interval comparison.

That means the natural next move is Stage 252 = the full support-`<=5` completion theorem, i.e. splice the unique five-coordinate interior packet to the already-finished support-`<=4` ledger and close the entire local mixed-ray search.
