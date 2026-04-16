# Moving-Throat PDE — Stage 226: Relaxed-Constraint Branch Declaration and Short-Range Open-System Compiler

## Status

**Exact within the carried Stage-225 placement/orbit packet, the exact projected leakage identity, the one-port short-range same-charge kernel verdict, and the declared reduced non-rigid / compensated-source closure used to translate the barrier session into the PDE stack.**

This stage does **not** yet compute the lowered barrier itself.
It declares the first honest post-Stage-225 branch family in which three standard-recovery suppressions are lifted in a controlled way:

1. transverse-current suppression `J^w \approx 0`,
2. rigid-mouth factorization / orbit slice `U = V = 0`,
3. positive-source Family-1 mouth/core closure.

The main outcome is that the Session-I barrier corridor now enters the derivation stack as a **codimension-three lift** of the Stage-225 front end, while the same-charge long-range verdict remains unchanged:

- no new static long-range attractive law,
- no new linear dynamic kernel class,
- only short-range/open-system reshaping of the already-carried one-port families.

---

## Purpose

Stage 225 ended with a concrete front-end branch test:

- actual selected twin-support placement,
- actual coherent orbit packet,
- separate outgoing finish line.

But the barrier session deliberately stepped **outside** the strict standard-recovery slice by relaxing three assumptions that the compact program had been carrying as part of the far-field / rigid-mouth / positive-source corridor.

So the next derivation move is not yet another barrier formula.
It is a branch declaration.

Stage 226 therefore does four things:

1. fixes the precise **base slice** inherited from Stage 225,
2. declares the three lifted closure lanes,
3. proves the exact **recovery map** back to the Stage-225 slice,
4. and hard-wires the crucial short-range theorem so the relaxed branch cannot later be misread as a hidden long-range same-charge force.

---

## Provenance

This stage is the moving-throat translation of the first barrier-session stress test into the PDE program.

Conceptually it sits directly after:

- **Stage 225**, which turned selected support placement and coherent orbit lock into an actual branch packet,
- the compact same-charge verdict, which already narrowed the one-port mixed sector to short-range static families and resonant short-range linear dynamics,
- and the barrier-session Section I, which relaxed `J^w \approx 0`, rigid-mouth factorization, and positive-source closure simultaneously.

So Stage 226 is the declaration step that makes those relaxations part of the formal derivation tree rather than a side calculation.

---

## 0. Why this stage is needed

Without this step, the barrier-session branch sits outside the audited PDE stack in an awkward way.
The compact program already says that the actual same-charge bottleneck is now the front-end packet

\[
(\Delta_{\rm norm},\Xi_1),
\qquad
U = V = 0,
\qquad
\text{selected twin-support placement},
\]

and that the mixed channels

\[
A_w,\qquad J^w,\qquad F_{\mu w},\qquad E_w,\qquad C_a
\]

are suppressed only in the strict far-field brane reduction, not deleted from the ontology.

So if we want to compare the standard branch to the relaxed session branch honestly, we need one intermediate theorem saying:

> what is being lifted, what stays fixed, and what does **not** change about the same-charge kernel class?

That is exactly the job of Stage 226.

---

## 1. Carried base branch from Stage 225

Stage 225 returns the concrete front-end packet

\[
\mathcal P_{225}
=
\Bigl(
\epsilon,
\varrho_{\rm phys},
\sigma_{\rm phys},
\text{ranking region},
R_{\rm tr},
R_{\rm target},
\epsilon_\eta,
 d\ln R_{\rm tr},
 d\ln R_{\rm target},
 d\ln \epsilon_\eta,
N_Q-1
\Bigr),
\]

with separate support packet

\[
\mathcal S_{225}
=
\bigl(
\zeta,
M_{\rm mix},
S(\zeta;\epsilon),
M_{\rm tr}
\bigr).
\]

The compact rigid-mouth normal form then identifies the physical post-static variables

\[
U := \ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V := \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\]

and on the strict rigid-mouth orbit slice one asks for

\[
U = 0,
\qquad
V = 0.
\]

For the present stage, the standard-recovery slice is therefore taken to be

\[
\boxed{
\mathfrak B_{\rm std}
:=
\Bigl\{
J^w = 0,
\ U = 0,
\ V = 0,
\ \varsigma(z) \equiv 1
\Bigr\}
\subset
\bigl(\mathcal P_{225},\mathcal S_{225}\bigr).
}
\]

Here `\varsigma(z)` is the Session-I mouth/core source profile.
It is written with a new symbol on purpose so it is **not** confused with the Stage-225 selected-support coordinate `\sigma_{\rm phys}`.

So Stage 226 starts from the precise standard slice

- no transverse leakage,
- rigid-mouth orbit lock,
- uncompensated positive-source mouth/core input.

---

## 2. The three lifted closure lanes

We now declare three independent reduced lifts.

### 2.1 Open-system leakage/work lane

The exact projected continuity identity already gives

\[
S_{\rm leak}
=
-\bigl[W(w)j^w\bigr]_{-\infty}^{+\infty}
+
\int_{-\infty}^{+\infty} W'(w)j^w(w)\,dw.
\]

To make the lift concrete, choose the normalized Gaussian projector

\[
W(w)=\frac{e^{-w^2}}{\sqrt\pi}
\]

and the first parity-odd reduced transverse current profile

\[
j^w(w)=\ell_w j_0\,w e^{-w^2}.
\]

Then the boundary term vanishes and the exact leakage source is

\[
\boxed{
S_{\rm leak}
=
\int_{-\infty}^{+\infty}W'(w)j^w(w)\,dw
=
-\frac{\sqrt2}{4}\,\ell_w j_0.
}
\]

So the leakage channel is nonzero as soon as the transverse-current lift `\ell_w` is turned on.

To accompany it, take the first parity-matched transverse electric field profile

\[
E_w(w)=E_0\,w e^{-w^2}.
\]

Then the reduced scalar-photon work channel becomes

\[
\boxed{
\mathcal W_w
:=
\int_{-\infty}^{+\infty} j^w(w)E_w(w)\,dw
=
\frac{\sqrt{2\pi}}{8}\,\ell_w j_0 E_0.
}
\]

So the open-system lift is carried by the exact pair

\[
\boxed{
L_w := \bigl(S_{\rm leak},\mathcal W_w\bigr),
}
\]

with both components vanishing on `\ell_w = 0`.

The sign of `S_{\rm leak}` depends on the orientation of the chosen odd profile.
Its magnitude is the physically relevant point here: the branch is no longer closed once `J^w` is allowed.

---

### 2.2 Non-rigid mouth / dressing lane

The compact program already diagonalized the strict rigid-mouth packet in the physical logarithmic variables `(U,V)`.
To declare the first non-rigid lift, use the minimal reduced free energy

\[
\mathcal F_{UV}(U,V)
=
\frac12 k_U U^2
+
\frac12 k_V V^2
-
\chi_\lambda U V
-
f_U U.
\]

Here:

- `k_U > 0` is the transfer-shape stiffness,
- `k_V > 0` is the dressing-leg stiffness,
- `\chi_\lambda` is the first non-rigid mouth coupling,
- `f_U` is the reduced same-charge forcing on the transfer leg.

The stationarity equations are

\[
k_U U - \chi_\lambda V = f_U,
\qquad
k_V V - \chi_\lambda U = 0.
\]

Solving exactly gives

\[
\boxed{
U = \frac{k_V f_U}{k_U k_V - \chi_\lambda^2},
\qquad
V = \frac{\chi_\lambda f_U}{k_U k_V - \chi_\lambda^2}.
}
\]

So the rigid/factorized and non-rigid cases separate sharply:

- if `f_U = 0`, then `U = V = 0`,
- if `f_U \neq 0` and `\chi_\lambda = 0`, then `V = 0`,
- if `f_U \neq 0` and `\chi_\lambda \neq 0`, then `V \neq 0` exactly.

The exact response ratio is

\[
\boxed{
\frac{V}{U} = \frac{\chi_\lambda}{k_V}.
}
\]

So the dressing leg is activated linearly once rigid-mouth factorization is lifted.

The stability Hessian is

\[
H_{UV}=
\begin{pmatrix}
k_U & -\chi_\lambda \\
-\chi_\lambda & k_V
\end{pmatrix},
\qquad
\det H_{UV}=k_U k_V - \chi_\lambda^2.
\]

Hence the reduced non-rigid branch is admissible precisely when

\[
\boxed{k_U k_V - \chi_\lambda^2 > 0.}
\]

On that admissible branch, the induced drain term is

\[
\boxed{
\mathcal D_{UV}:=\chi_\lambda U V
=
\frac{\chi_\lambda^2 k_V f_U^2}{(k_U k_V - \chi_\lambda^2)^2}
\ge 0.
}
\]

So the minimal reduced closure already reproduces the qualitative Session-I statement:

> once the rigid-mouth factorization is broken, the dressing leg is no longer inert and the transfer forcing drains positive energy into it.

At infinitesimal order, this lane is the relaxed continuation of the compact transfer-shape scalar:

\[
\Xi_1 = \delta U.
\]

So the later explicit `\Xi_1` bookkeeping belongs to the same lifted lane.

---

### 2.3 Compensated sign-changing source lane

To lift the positive-source Family-1 corridor without changing the mean source strength, take the minimal compensated profile

\[
\boxed{
\varsigma(z)=1+a\cos(\pi z)+b\cos(2\pi z),
\qquad z\in[0,1].
}
\]

Its mean is exactly preserved:

\[
\boxed{
\int_0^1 \varsigma(z)\,dz = 1.
}
\]

So this branch changes source *shape* without changing total normalized source weight.

Now write

\[
y := \cos(\pi z),
\qquad y\in[-1,1],
\]

so that

\[
\cos(2\pi z)=2y^2-1.
\]

Then the compensated source becomes the exact quadratic

\[
\boxed{
\varsigma(y)=1-b+a y+2b y^2,
\qquad y\in[-1,1].
}
\]

The interior stationary point is

\[
\boxed{
y_*=-\frac{a}{4b}}
\]

when `b \neq 0`, and the corresponding vertex value is

\[
\boxed{
\varsigma(y_*) = 1-b-\frac{a^2}{8b}.
}
\]

So a sign-changing compensated branch exists whenever the quadratic takes a negative value somewhere on `[-1,1]`.
In particular:

- the boundary candidates are
  \[
  \varsigma(1)=1+a+b,
  \qquad
  \varsigma(-1)=1-a+b,
  \]
- if `b>0` and `|a|\le 4b`, the interior minimum is exact and sign change occurs whenever
  \[
  1-b-\frac{a^2}{8b}<0.
  \]

So the compensated-source lane is not a hand-waving statement.
It is an explicit zero-mean shape deformation with an exact sign-change test.

We package this lift as

\[
\boxed{
L_\varsigma := \bigl(\varsigma(z),\ \varsigma_{\min},\ \mathrm{signchg}\bigr),
}
\]

where `signchg` is the boolean condition `\varsigma_{\min}<0`.

---

## 3. Exact relaxed-branch declaration

We can now declare the relaxed family directly.

\[
\boxed{
\mathfrak B_{226}^{\rm relax}
:=
\Bigl\{
(\mathcal P_{225},\mathcal S_{225});
L_w,\ L_{UV},\ L_\varsigma
\Bigr\}
}
\]

with the explicit lifted lanes

\[
L_w = (S_{\rm leak},\mathcal W_w),
\qquad
L_{UV}=(U,V,\mathcal D_{UV}),
\qquad
L_\varsigma=(\varsigma,\varsigma_{\min},\mathrm{signchg}).
\]

So the relaxed branch is **not** a replacement for the Stage-225 packet.
It is an **augmentation** of that packet by three new channels:

1. open-system export,
2. non-rigid mouth/dressing response,
3. compensated source-shape response.

The selected twin-support placement packet from Stage 225 is carried through unchanged.
The relaxed branch only adds new reduced lanes around it.

---

## 4. Exact recovery map back to the Stage-225 slice

The Stage-225/compact front end is recovered on the codimension-three slice

\[
\boxed{
\ell_w = 0,
\qquad
f_U = 0,
\qquad
a = b = 0.
}
\]

Indeed:

\[
\ell_w=0
\implies
S_{\rm leak}=0,
\quad
\mathcal W_w=0,
\]

\[
f_U=0
\implies
U=0,
\quad
V=0,
\quad
\mathcal D_{UV}=0,
\]

\[
a=b=0
\implies
\varsigma(z)\equiv 1,
\quad
\varsigma_{\min}=1,
\quad
\mathrm{signchg}=\mathrm{False}.
\]

So the relaxed branch declaration is honest in both directions:

- it extends the carried front end,
- and it reduces exactly back to the Stage-225 standard-recovery slice.

---

## 5. Short-range kernel-span theorem for the lifted branch

This is the most important structural statement of the stage.

The lifted branch is **not** allowed to change the already-audited same-charge long-range verdict.
So we now record that theorem directly in the derivation stack.

Take the same primitive reduced source profiles already isolated by the one-port same-charge audit:

\[
\mathcal S_Q(x)=\frac{1}{x^3},
\qquad
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
\qquad \kappa>0.
\]

For the one-port static mixed bundle, the exact conservative correction is of the form

\[
\boxed{
\delta V_{\rm stat}(x)
=
-\frac12\left[
\frac{\mathcal C_6}{x^6}
+
2\mathcal C_4\frac{e^{-2\kappa x}}{x^4}
+
\mathcal C_2\frac{e^{-4\kappa x}}{x^2}
\right],
}
\]

and the linear monochromatic dynamic bundle keeps the **same** spatial span with frequency-dependent coefficients,

\[
\boxed{
\Re\,\mathfrak V_{\rm dyn}(x,\omega)
=
-\frac12\left[
\frac{\mathcal C_6(\omega)}{x^6}
+
2\mathcal C_4(\omega)\frac{e^{-2\kappa x}}{x^4}
+
\mathcal C_2(\omega)\frac{e^{-4\kappa x}}{x^2}
\right].
}
\]

Therefore both the static and the linear dynamic conservative corrections lie in the exact kernel span

\[
\boxed{
\mathcal K_{\rm SR}
=
\mathrm{span}\!\left\{
\frac{1}{x^6},
\frac{e^{-2\kappa x}}{x^4},
\frac{e^{-4\kappa x}}{x^2}
\right\}.
}
\]

In particular,

\[
\boxed{
\lim_{x\to\infty} x\,\delta V_{\rm stat}(x)=0,
\qquad
\lim_{x\to\infty} x\,\Re\mathfrak V_{\rm dyn}(x,\omega)=0.
}
\]

So the relaxed branch declaration hard-wires the following invariant:

> turning on `J^w`, turning on non-rigid `U/V` response, or allowing compensated sign-changing source data does **not** adjoin a new `1/x`-type or Yukawa-`1/x` same-charge attraction.

What the lift can do is only:

1. change the coefficients of the already-carried short-range families,
2. open energy-export channels through `S_{\rm leak}` and `\mathcal W_w`,
3. move energy between the transfer and dressing legs through `\mathcal D_{UV}`.

That is the exact short-range/open-system meaning of the Session-I corridor.

---

## 6. Practical Stage-226 output packet

The smallest reduced packet that Stage 226 asks the later derivation tree to carry is

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

If the linear weak-axisymmetric packet is being used, then `\Xi_1` is simply the infinitesimal transfer-shape component of the same `U` lane.

So the correct reading of the relaxed branch is:

- Stage 225 still fixes the **placement/orbit/outgoing** front end,
- Stage 226 adds the first **open-system / non-rigid / compensated-source** lanes,
- and the short-range kernel class remains frozen.

---

## 7. What this stage achieves

Stage 226 closes four concrete bookkeeping gaps.

### 7.1 It imports the Session-I relaxations without breaking the compact firewall

The three relaxed assumptions are now represented directly inside the derivation stack:

- `J^w \neq 0`,
- non-rigid `(U,V)` dynamics,
- compensated sign-changing `\varsigma(z)`.

### 7.2 It proves exact recovery of the carried Stage-225 slice

The relaxed family is not a different theory.
It is a codimension-three lift of the already-carried selected-support / orbit packet.

### 7.3 It records the reduced energetic meaning of the new lanes

- `S_{\rm leak}` and `\mathcal W_w` are the open-system export channels,
- `\mathcal D_{UV}` is the transfer-to-dressing drain,
- `\varsigma(z)` is a zero-mean source-shape deformation that can change sign.

### 7.4 It hard-wires the short-range theorem

This is the key conceptual safeguard.
The relaxed same-charge branch is now formally declared to be a **short-range/open-system bypass** branch.
It is **not** allowed to be reinterpreted later as a hidden long-range attractive law.

---

## 8. Immediate next derivation step

The next honest continuation is now sharply defined.

Stage 227 should take the newly declared open-system lane and compile the reduced observables that actually turn on once `J^w` is allowed:

\[
S_{\rm leak},
\qquad
J^w E_w,
\qquad
\text{and their pullback into the same selected placement/orbit packet.}
\]

That is the right first quantitative continuation of the relaxed branch.
Only after that does it make sense to fold in the full non-rigid mouth packet and the compensated source compiler.

---

## 9. SymPy-backed status

The accompanying audit script verifies all of the algebra used here:

1. the exact Gaussian/odd-profile leakage reduction
   \[
   S_{\rm leak}=-\frac{\sqrt2}{4}\,\ell_w j_0,
   \]
2. the exact Gaussian/odd-profile work channel
   \[
   \mathcal W_w=\frac{\sqrt{2\pi}}{8}\,\ell_w j_0 E_0,
   \]
3. the exact non-rigid mouth response
   \[
   U = \frac{k_V f_U}{k_Uk_V-\chi_\lambda^2},
   \qquad
   V = \frac{\chi_\lambda f_U}{k_Uk_V-\chi_\lambda^2},
   \]
4. the exact drain identity
   \[
   \mathcal D_{UV}=\chi_\lambda U V
   =
   \frac{\chi_\lambda^2 k_V f_U^2}{(k_Uk_V-\chi_\lambda^2)^2},
   \]
5. the exact compensated-source normalization
   \[
   \int_0^1 \varsigma(z)\,dz = 1,
   \]
6. the exact quadratic rewrite
   \[
   \varsigma(y)=1-b+a y+2b y^2,
   \qquad y=\cos(\pi z),
   \]
7. the interior stationary point and candidate minimum
   \[
   y_*=-\frac{a}{4b},
   \qquad
   \varsigma(y_*)=1-b-\frac{a^2}{8b},
   \]
8. the exact standard-recovery slice
   \[
   \ell_w=0,
   \quad f_U=0,
   \quad a=b=0
   \implies
   S_{\rm leak}=\mathcal W_w=U=V=\mathcal D_{UV}=0,
   \quad \varsigma\equiv1,
   \]
9. and the strict short-range limits
   \[
   \lim_{x\to\infty}x\,\delta V_{\rm stat}(x)=0,
   \qquad
   \lim_{x\to\infty}x\,\Re\mathfrak V_{\rm dyn}(x,\omega)=0.
   \]

Supporting file:
- `moving_throat_pde_stage226_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
