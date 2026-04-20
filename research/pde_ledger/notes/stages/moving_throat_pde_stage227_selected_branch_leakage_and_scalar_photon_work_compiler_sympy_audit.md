# Moving-Throat PDE — Stage 227: Selected-Branch Leakage and Scalar-Photon Work Compiler

## Status

**Exact within the carried Stage-226 relaxed open-system branch declaration, the exact projected continuity / EM energy identities, the exact Stage-225 selected-support pullback, and the declared one-mode scalar-photon amplitude closure used to translate Session I into the PDE stack.**

This stage does **not** yet lower the barrier by itself.
It compiles the first two observables that actually switch on once the `J^w \neq 0` lane is allowed:

1. the projected leakage source `S_{\rm leak}`,
2. the scalar-photon work channel `J^w E_w` in the reduced one-mode sense used by the barrier session.

The main outcome is that these observables are now attached directly to the same **selected-support packet** that Stage 225 already fixed:

- they pull back through the selected support demand `\Pi_{\rm tr}`,
- they depend only on the support-side variables `(\Lambda,\epsilon)` or `(\Lambda,\varrho)`,
- and they remain separate from the orbit-lock packet.

---

## Purpose

Stage 226 declared the open-system branch family but stopped at the declaration step.
It told us that once `J^w` is no longer suppressed, the PDE stack must carry explicit leakage and work observables, but it did not yet compile those observables onto the actual Stage-225 selected branch.

Stage 227 closes that gap.

It does five things:

1. restates the exact parent-theory identities that make `S_{\rm leak}` and `J^wE_w` legitimate observables,
2. fixes a minimal odd Gaussian scalar-photon profile that gives an exact leakage compiler,
3. defines the exact one-mode bulk work scalar and the reduced Session-I work scalar,
4. pulls the amplitude back through the actual selected-support demand `\Pi_{\rm tr}` from Stage 225,
5. proves that the whole leakage/work lane is **support-side**, not orbit-side.

So this stage is the first quantitative continuation of the relaxed branch declaration.

---

## Provenance

This stage sits directly after:

- **Stage 225**, which fixed the actual selected-support packet and separated support placement from orbit lock,
- **Stage 226**, which declared the open-system lift as a codimension-three extension of the standard branch,
- the exact projection and mixed-sector energy identities from the parent 4D / plasma stack,
- and the Session-I barrier run, which tracked explicit `S_{\rm leak}` and `J^wE_w` channels once `J^w` was allowed.

So Stage 227 is the first point where the relaxed open-system lane is no longer just “allowed.”
It is explicitly compiled into the moving-throat derivation tree.

---

## 0. Why this stage is needed

The compact program already distinguishes two front-end packets:

1. the **selected-support** packet,
2. the **orbit-lock** packet.

Stage 225 made that split concrete.
It showed that the selected support point is fixed by the coherent support variables, while the coherent orbit-lock test is carried by the separate support-blind orbit packet.

Stage 226 then declared the relaxed branch but did not yet say where the new leakage/work observables attach.

That point matters.
If these observables are really part of the same same-charge front end, then they should not float freely.
They should either:

- attach to the selected-support packet,
- attach to the orbit packet,
- or mix the two.

Stage 227 proves that, in the first honest reduced compiler, they attach to the **selected-support packet alone**.

---

## 1. Exact parent-theory identities that turn on when `J^w` is allowed

The exact projected continuity identity already gives

\[
\boxed{
S_{\rm leak}
=
-\bigl[W(w)j^w\bigr]_{-\infty}^{+\infty}
+
\int_{-\infty}^{+\infty}W'(w)j^w(w)\,dw.
}
\]

So once `j^w` is nonzero, the projected brane subsystem is explicitly open.

The exact localized-Maxwell energy ledger gives

\[
\boxed{
\partial_t u_{\rm EM}+\partial_A S_{\rm EM}^A
=
-J^A E_A
=
-(J^aE_a+J^wE_w).
}
\]

So `J^wE_w` is an explicit scalar-photon work channel as soon as the mixed field `E_w` is active.

In projected form one has

\[
\partial_t\overline{u_{\rm EM}}+\partial_a\overline{S_{\rm EM}^a}
=
-\overline{J^A E_A}
+\mathcal L_{\rm EM}^{(w)},
\]

with

\[
\mathcal L_{\rm EM}^{(w)}
=
-\bigl[W S_{\rm EM}^w\bigr]_{-\infty}^{+\infty}
+
\int W'(w)S_{\rm EM}^w\,dw.
\]

Stage 227 keeps the leakage source `S_{\rm leak}` and the scalar-photon work channel `J^wE_w` as the first reduced observables in the newly opened lane.

---

## 2. Minimal odd Gaussian leakage profile and exact projected leakage compiler

Choose the normalized Gaussian projector

\[
\boxed{
W_\lambda(w)
=
\frac{e^{-w^2/\lambda^2}}{\lambda\sqrt\pi}.
}
\]

Introduce the first odd scalar-photon profile

\[
\boxed{
\phi_\lambda(w)
=
\frac{2w}{\sqrt\pi\,\lambda^3}
\,e^{-w^2/\lambda^2}.
}
\]

Now define the reduced mixed field and current by

\[
\boxed{
E_w(w;r)=-\mathcal E_0(r)\,\phi_\lambda(w),
\qquad
j^w(w;r)=\mu_w\rho_0\,E_w(w;r),
\qquad
J^w=q\,j^w.
}
\]

Here:

- `\mathcal E_0(r)` is the selected-branch scalar-photon amplitude,
- `\mu_w` is the reduced transverse mobility,
- `\rho_0` is the carried density scale,
- `q` is the reduced charge label used in the session barrier closure.

Because both `W_\lambda` and `\phi_\lambda` decay rapidly, the boundary term vanishes exactly:

\[
\bigl[W_\lambda j^w\bigr]_{-\infty}^{+\infty}=0.
\]

The leakage source is therefore

\[
S_{\rm leak}(r)
=
\int_{-\infty}^{+\infty}W_\lambda'(w)j^w(w;r)\,dw.
\]

A direct evaluation gives

\[
\boxed{
S_{\rm leak}(r)
=
\frac{\sqrt2\,\mu_w\rho_0}{2\sqrt\pi\,\lambda^3}
\,\mathcal E_0(r).
}
\]

So the first honest selected-branch leakage scalar is linear in the reduced scalar-photon amplitude `\mathcal E_0(r)`.

This is the exact projected-continuity compiler for the chosen one-mode Gaussian closure.

---

## 3. Exact one-mode bulk work scalar and the reduced Session-I work scalar

The same odd Gaussian lane carries an exact one-mode bulk work scalar

\[
\boxed{
\mathcal W_w^{\rm bulk}(r)
:=
q\int_{-\infty}^{+\infty}j^w(w;r)E_w(w;r)\,dw.
}
\]

For the chosen profile this evaluates exactly to

\[
\boxed{
\mathcal W_w^{\rm bulk}(r)
=
\frac{\sqrt2\,\mu_w q\rho_0}{2\sqrt\pi\,\lambda^3}
\,\mathcal E_0(r)^2.
}
\]

So the one-mode bulk work scalar obeys the exact relation

\[
\boxed{
\mathcal W_w^{\rm bulk}(r)=q\,\mathcal E_0(r)\,S_{\rm leak}(r).
}
\]

Thus, once the selected-branch amplitude is fixed, the one-mode bulk work lane is not independent of the leakage lane.

### 3.1 Reduced Session-I scalar-photon work compiler

The barrier session tracked a reduced scalar-photon work scalar of the form

\[
\boxed{
\mathcal W_w^{\rm sess}(r)
:=
\frac{2\mu_w q\rho_0}{\lambda^2}\,\mathcal E_0(r)^2.
}
\]

Within the present one-mode Gaussian closure this is exactly the thickness-rescaled bulk work scalar

\[
\boxed{
\mathcal W_w^{\rm sess}(r)
=
2\sqrt{2\pi}\,\lambda\,\mathcal W_w^{\rm bulk}(r).
}
\]

So the Session-I work scalar is not a separate physical lane.
It is the wall-thickness rescaling of the same one-mode bulk scalar-photon work channel.

Combining this with the leakage compiler gives the exact quadratic law

\[
\boxed{
\mathcal W_w^{\rm sess}(r)
=
\frac{4\pi q\lambda^4}{\mu_w\rho_0}
\,S_{\rm leak}(r)^2.
}
\]

So the reduced Session-I work scalar is quadratic in the projected leakage amplitude on the chosen one-mode branch.

---

## 4. Pullback into the Stage-225 selected-support packet

Stage 225 already fixed the exact selected-support demand

\[
\Pi_{\rm tr}=
\frac43\,C_{\rm mix},
\qquad
C_{\rm mix}=
\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

Hence

\[
\boxed{
\Pi_{\rm tr}(r)
=
\frac{32\Lambda(r)(1-\epsilon(r))}{3\pi^2}.
}
\]

Using the actual selected-branch coordinate

\[
\varrho(r)=\frac23(1-\epsilon(r)),
\]

this can be rewritten exactly as

\[
\boxed{
\Pi_{\rm tr}(r)=\frac{16}{\pi^2}\,\Lambda(r)\,\varrho(r).
}
\]

Stage 227 now chooses the simplest selected-branch amplitude pullback:

\[
\boxed{
\mathcal E_0(r):=\eta_{\rm leak}\,\Pi_{\rm tr}(r),
}
\]

where `\eta_{\rm leak}` is the reduced leakage-strength parameter already used in the session barrier closure.

### 4.1 Leakage source on the selected branch

Substituting into the exact leakage compiler gives

\[
\boxed{
S_{\rm leak}(r)
=
\frac{16\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{3\pi^{5/2}\lambda^3}
\,\Lambda(r)(1-\epsilon(r)).
}
\]

Equivalently,

\[
\boxed{
S_{\rm leak}(r)
=
\frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
\,\Lambda(r)\varrho(r).
}
\]

### 4.2 One-mode bulk work scalar on the selected branch

The exact one-mode bulk work scalar becomes

\[
\boxed{
\mathcal W_w^{\rm bulk}(r)
=
\frac{512\sqrt2\,\eta_{\rm leak}^2\mu_w q\rho_0}{9\pi^{9/2}\lambda^3}
\,\Lambda(r)^2(1-\epsilon(r))^2.
}
\]

Equivalently,

\[
\boxed{
\mathcal W_w^{\rm bulk}(r)
=
\frac{128\sqrt2\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^{9/2}\lambda^3}
\,\Lambda(r)^2\varrho(r)^2.
}
\]

### 4.3 Reduced Session-I work scalar on the selected branch

The reduced Session-I scalar becomes

\[
\boxed{
\mathcal W_w^{\rm sess}(r)
=
\frac{2048\,\eta_{\rm leak}^2\mu_w q\rho_0}{9\pi^4\lambda^2}
\,\Lambda(r)^2(1-\epsilon(r))^2.
}
\]

Equivalently,

\[
\boxed{
\mathcal W_w^{\rm sess}(r)
=
\frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
\,\Lambda(r)^2\varrho(r)^2.
}
\]

So the selected-branch leakage/work packet is now fully explicit in the realized support variables.

---

## 5. Exact support-versus-orbit split and recovery slice

The key structural result of the stage is now immediate.

All three selected-branch formulas

\[
S_{\rm leak}(r),
\qquad
\mathcal W_w^{\rm bulk}(r),
\qquad
\mathcal W_w^{\rm sess}(r)
\]

depend only on

\[
\Lambda(r),\qquad \epsilon(r)
\quad\text{or equivalently}\quad
\Lambda(r),\qquad \varrho(r).
\]

They do **not** depend on the separate coherent orbit variables

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta.
\]

So the exact support-side factorization is

\[
\boxed{
\partial_{R_{\rm tr}} S_{\rm leak}
=
\partial_{R_{\rm target}} S_{\rm leak}
=
\partial_{\epsilon_\eta} S_{\rm leak}=0,
}
\]

and likewise for both work scalars.

This proves that Stage 227 is a **selected-support compiler**, not an orbit-lock compiler.

### 5.1 Orientation parity

The pullback parameter `\eta_{\rm leak}` enters the compiled observables with exact parity

\[
\boxed{
S_{\rm leak}(-\eta_{\rm leak})=-S_{\rm leak}(\eta_{\rm leak}),
}
\]

while both work scalars are even:

\[
\boxed{
\mathcal W_w^{\rm bulk}(-\eta_{\rm leak})
=
\mathcal W_w^{\rm bulk}(\eta_{\rm leak}),
\qquad
\mathcal W_w^{\rm sess}(-\eta_{\rm leak})
=
\mathcal W_w^{\rm sess}(\eta_{\rm leak}).
}
\]

So changing the transport orientation flips the sign of net leakage but leaves the exported-work magnitudes unchanged.

### 5.2 Recovery of the standard slice

The Stage-226 standard-recovery slice is recovered at

\[
\boxed{
\eta_{\rm leak}=0.
}
\]

Then

\[
\mathcal E_0(r)=0,
\qquad
S_{\rm leak}(r)=0,
\qquad
\mathcal W_w^{\rm bulk}(r)=0,
\qquad
\mathcal W_w^{\rm sess}(r)=0.
\]

So the Stage-227 compiler reduces exactly back to the carried Stage-225 / Stage-226 front end on the closed-branch slice.

---

## 6. Practical Stage-227 output packet

The smallest packet that this stage asks the later derivation tree to carry is

\[
\boxed{
\mathcal P_{227}^{\rm leak}
=
\Bigl(
\mathcal P_{225};
\Pi_{\rm tr},
\mathcal E_0,
S_{\rm leak},
\mathcal W_w^{\rm bulk},
\mathcal W_w^{\rm sess}
\Bigr).
}
\]

Equivalently, on the realized selected branch one may use the packet

\[
\boxed{
\mathcal P_{227}^{\rm leak}
=
\Bigl(
\mathcal P_{225};
\Lambda,
\epsilon,
\eta_{\rm leak},
S_{\rm leak},
\mathcal W_w^{\rm bulk},
\mathcal W_w^{\rm sess}
\Bigr).
}
\]

So the first quantitative open-system lane is now fully attached to the same selected-support packet that Stage 225 had already fixed.

---

## 7. What this stage achieves

Stage 227 closes five concrete bookkeeping gaps.

### 7.1 It turns the Stage-226 open-system declaration into actual observables

The leakage/work lane is no longer only “allowed.”
It now carries explicit scalars:

- `S_{\rm leak}`,
- `\mathcal W_w^{\rm bulk}`,
- `\mathcal W_w^{\rm sess}`.

### 7.2 It ties the open-system lane to the selected-support packet

The amplitude pullback

\[
\mathcal E_0=\eta_{\rm leak}\Pi_{\rm tr}
\]

forces the entire leakage/work compiler through the same selected-support demand that Stage 225 already fixed.

### 7.3 It proves the support-versus-orbit split

The leakage/work lane depends only on the support-side variables `(\Lambda,\epsilon)` or `(\Lambda,\varrho)`.
It does not depend on the separate orbit packet.

### 7.4 It clarifies the relation between the exact one-mode bulk work scalar and the Session-I work scalar

The Session-I scalar-photon work quantity is the thickness-rescaled form of the same one-mode bulk work channel.
So the session quantity is not a different lane; it is a reduced compiler of the same mixed-sector work observable.

### 7.5 It records the exact parity structure of the leakage/work lane

- leakage is odd in the transport-orientation parameter,
- work magnitude is even.

That is the first honest sign theorem for the selected-branch open-system packet.

---

## 8. Immediate next derivation step

The next honest continuation is now sharply defined.

Stage 228 should take the already-declared non-rigid mouth lane `(U,V)` from Stage 226 and compile it onto the same selected-support branch, so that one can study:

\[
U(r),\qquad V(r),\qquad \mathcal D_{UV}(r)
\]

side by side with

\[
S_{\rm leak}(r),\qquad \mathcal W_w^{\rm sess}(r).
\]

That is the right point to join the open-system export lane to the dressing-leg drain lane in one packet.

---

## 9. SymPy-backed status

The accompanying audit script verifies all of the algebra used here:

1. the exact Gaussian leakage compiler
   \[
   S_{\rm leak}(r)
   =
   \frac{\sqrt2\,\mu_w\rho_0}{2\sqrt\pi\,\lambda^3}
   \,\mathcal E_0(r),
   \]
2. the exact one-mode bulk work scalar
   \[
   \mathcal W_w^{\rm bulk}(r)
   =
   \frac{\sqrt2\,\mu_w q\rho_0}{2\sqrt\pi\,\lambda^3}
   \,\mathcal E_0(r)^2,
   \]
3. the exact relation
   \[
   \mathcal W_w^{\rm bulk}=q\,\mathcal E_0\,S_{\rm leak},
   \]
4. the exact reduced Session-I compiler
   \[
   \mathcal W_w^{\rm sess}
   =
   2\sqrt{2\pi}\,\lambda\,\mathcal W_w^{\rm bulk}
   =
   \frac{2\mu_w q\rho_0}{\lambda^2}\,\mathcal E_0^2,
   \]
5. the exact quadratic leakage-to-work law
   \[
   \mathcal W_w^{\rm sess}
   =
   \frac{4\pi q\lambda^4}{\mu_w\rho_0}
   \,S_{\rm leak}^2,
   \]
6. the exact selected-support pullback
   \[
   \Pi_{\rm tr}
   =
   \frac{32\Lambda(1-\epsilon)}{3\pi^2}
   =
   \frac{16\Lambda\varrho}{\pi^2},
   \qquad
   \mathcal E_0=\eta_{\rm leak}\Pi_{\rm tr},
   \]
7. the resulting compiled selected-branch formulas for
   `S_{\rm leak}`, `\mathcal W_w^{\rm bulk}`, and `\mathcal W_w^{\rm sess}`,
8. the exact support-versus-orbit split
   \[
   \partial_{R_{\rm tr}}=
   \partial_{R_{\rm target}}=
   \partial_{\epsilon_\eta}=0
   \quad\text{on the compiled leakage/work lane},
   \]
9. the exact orientation parity
   \[
   S_{\rm leak}(-\eta_{\rm leak})=-S_{\rm leak}(\eta_{\rm leak}),
   \quad
   \mathcal W(-\eta_{\rm leak})=\mathcal W(\eta_{\rm leak}),
   \]
10. and exact recovery of the closed standard slice at `\eta_{\rm leak}=0`.

Supporting file:
- `moving_throat_pde_stage227_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py`
