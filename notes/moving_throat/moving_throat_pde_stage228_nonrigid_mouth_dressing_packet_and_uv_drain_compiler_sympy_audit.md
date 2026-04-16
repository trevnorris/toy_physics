# Moving-Throat PDE — Stage 228: Non-Rigid Mouth/Dressing Packet and `U/V` Drain Compiler

## Status

**Exact within the carried rigid-mouth actual-branch chart, the Stage-226 non-rigid reduced free-energy lane, and the declared reduced orbit-side forcing closure that lifts the rigid-mouth slice into a two-coordinate physical packet `(U,V)`.**

This stage does **not** yet lower the barrier by itself.
It compiles the second lifted lane from Stage 226 into the actual physical observables that the compact program already uses:

1. transfer-shape motion `\mathcal T^2`,
2. dressing motion `\epsilon_\eta`,
3. target-ratio motion `R_{\rm target}`,
4. the dependent microscopic correction,
5. and the positive `U/V` drain term.

The main outcome is that once rigid-mouth factorization is broken, the direct transfer packet is still carried by `U`, but the finite selected-branch ratio acquires a second exact leg through `V`; the dressing variable becomes active, and the resulting compiler is **orbit-side / support-blind** rather than a new support-placement lane.

---

## Purpose

Stage 226 declared the non-rigid mouth/dressing lane abstractly.
It introduced a minimal reduced `U/V` free energy and showed that a transfer-side forcing can activate the dressing leg once the mouth is allowed to flex.

But that declaration step stopped short of the actual physical compiler.
The compact program and the rigid-mouth barrier audit had already fixed the exact physical chart

\[
U=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\]

together with the exact selected-branch identity relating `\mathcal T^2`, `\epsilon_\eta`, and `R_{\rm target}`.

So the missing step is now sharp:

> take the abstract Stage-226 non-rigid `U/V` lane and compile it all the way into the actual finite physical packet, the weak-axisymmetric first-order packet, and the dependent microscopic correction.

That is what Stage 228 does.

---

## Provenance

This stage sits directly after:

- **Stage 226**, which declared the non-rigid mouth/dressing lift through a reduced `U/V` free energy,
- **Stage 227**, which proved the leakage/work compiler is selected-support side rather than orbit side,
- the compact rigid-mouth normal form, which already diagonalized the actual branch in `(U,V)`,
- and the barrier-session Section I, which reported a nonzero `V`, active `\epsilon_\eta`, and a positive `U/V` drain once the rigid-mouth assumption was relaxed.

So Stage 228 is the first place where the barrier session’s `U/V` lane is no longer just a reduced add-on.
It becomes part of the actual derivation stack.

---

## 0. Why this stage is needed

The compact program already says that, on the rigid-mouth actual branch, orbit lock is the Cartesian codimension-two condition

\[
U=0,\qquad V=0.
\]

Stage 226 then showed that a non-rigid mouth coupling and a transfer-side forcing produce a nonzero solution `(U,V)` on an admissible reduced branch.
But as long as that result remains in abstract free-energy form, the derivation stack is still missing three things:

1. the exact finite update of the physical observables `(\mathcal T^2,\epsilon_\eta,R_{\rm target})`,
2. the exact dependent microscopic correction induced by that finite packet,
3. the precise relation between the new finite packet and the earlier weak-axisymmetric scalar `\Xi_1`.

Without those compilers, one cannot tell whether the Session-I `U/V` branch is:

- a support-placement effect,
- an orbit-side deformation,
- or an uncontrolled mixture of the two.

Stage 228 closes that ambiguity.

---

## 1. Carried exact rigid-mouth chart

The rigid-mouth actual branch is already charted by the diagonal logarithmic pair

\[
\boxed{
U:=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V:=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]

So the direct observables are

\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2 e^U,
\qquad
\epsilon_\eta=\epsilon_{\eta,\rm ref}e^V.
}
\]

The selected-branch identity
\[
R_{\rm target}\,\mathcal T^2=\Lambda_0(1-\epsilon_\eta)
\]
then gives the exact target-ratio compiler
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}\,e^{-U}.
}
\]

So the rigid-mouth branch already tells us something important before any new closure is added:

- `U` is the direct transfer-shape leg,
- `V` is the direct dressing leg,
- and the target ratio is the commuting product of those two finite legs.

Stage 228 keeps this exact chart and only changes the law that determines the finite point `(U,V)`.

---

## 2. Reduced non-rigid mouth/dressing free energy

Stage 226 already declared the minimal reduced non-rigid lane.
For the present stage it is convenient to rewrite the same closure in the session notation

\[
\chi_{UV}(r):=g_{UV}\,\chi_\lambda(r),
\]

so the reduced free energy becomes

\[
\boxed{
\mathcal F_{\rm nr}(U,V;r)
=
\frac12 a_U U^2
+
\frac12 a_V V^2
-
\chi_{UV}(r)\,U V
-
f_U(r)\,U.
}
\]

Here:

- `a_U>0` is the transfer-shape stiffness,
- `a_V>0` is the dressing-leg stiffness,
- `\chi_{UV}(r)` is the signed non-rigid mouth coupling,
- `f_U(r)` is the orbit-side transfer forcing.

Two important limiting cases are immediate.

### 2.1 Rigid-mouth orbit slice

If the forcing is absent,

\[
f_U(r)=0,
\]

then the unique stationary point is

\[
U=0,\qquad V=0.
\]

So the non-rigid compiler still recovers the exact rigid-mouth orbit slice when the forcing is shut off.

### 2.2 Factorized-but-forced slice

If the mouth coupling is shut off,

\[
\chi_{UV}(r)=0,
\]

then the transfer leg survives but the dressing leg stays inert:

\[
U=\frac{f_U}{a_U},
\qquad
V=0.
\]

So the dressing activation is a genuinely non-rigid effect rather than a rephrasing of the old transfer leg.

---

## 3. Exact stationary `U/V` packet

The stationarity equations are linear:

\[
a_U U-\chi_{UV}V=f_U,
\qquad
a_V V-\chi_{UV}U=0.
\]

Define the non-rigid determinant

\[
\boxed{
\Delta_{UV}:=a_U a_V-\chi_{UV}^2.
}
\]

Then the exact stationary packet is

\[
\boxed{
U(r)=\frac{a_V f_U(r)}{\Delta_{UV}(r)},
\qquad
V(r)=\frac{\chi_{UV}(r)f_U(r)}{\Delta_{UV}(r)}.
}
\]

The exact response ratio is therefore

\[
\boxed{
\frac{V}{U}=\frac{\chi_{UV}}{a_V}.
}
\]

So the sign of `V/U` is the sign of the signed mouth coupling.
The Session-I branch with `U>0` and `V<0` is therefore the branch on which `\chi_{UV}<0`.

The admissibility condition is simply positive-definiteness of the Hessian:

\[
H_{\rm nr}=
\begin{pmatrix}
a_U & -\chi_{UV}\\
-\chi_{UV} & a_V
\end{pmatrix},
\qquad
\det H_{\rm nr}=\Delta_{UV}.
\]

Hence the reduced non-rigid packet is admissible precisely when

\[
\boxed{
\Delta_{UV}>0.
}
\]

This is the same condition already anticipated in Stage 226, now carried in the session notation.

---

## 4. Exact finite physical compiler

Substituting the exact stationary packet into the carried physical chart gives

\[
\boxed{
\mathcal T^2(r)
=
\mathcal T_{\rm ref}^2
\exp\!\left(\frac{a_V f_U(r)}{\Delta_{UV}(r)}\right),
}
\]

\[
\boxed{
\epsilon_\eta(r)
=
\epsilon_{\eta,\rm ref}
\exp\!\left(\frac{\chi_{UV}(r)f_U(r)}{\Delta_{UV}(r)}\right),
}
\]

\[
\boxed{
\frac{R_{\rm target}(r)}{R_{\rm target,ref}}
=
\frac{
1-\epsilon_{\eta,\rm ref}\exp\!\left(\chi_{UV}f_U/\Delta_{UV}\right)
}{
1-\epsilon_{\eta,\rm ref}
}
\exp\!\left(-a_V f_U/\Delta_{UV}\right).
}
\]

This is the exact finite compiler that Stage 226 did not yet provide.

It shows that once rigid-mouth factorization is broken, the reduced branch is no longer characterized by a single transfer coordinate alone:

- `U` moves `\mathcal T^2` directly,
- `V` moves `\epsilon_\eta` directly,
- and the selected-branch target ratio sees both legs simultaneously.

So the non-rigid packet is already larger than the original one-scalar rigid-mouth front end.

---

## 5. Exact dependent microscopic correction and positive drain

The rigid-mouth actual branch already carries the exact dependent-plane compiler

\[
\boxed{
\mathbf y_{\rm rm}^{\rm dep}(U,V)
=
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

So the first non-rigid packet induces the exact microscopic correction

\[
\boxed{
\mathbf y_{\rm nr}^{\rm dep}(r)
=
\begin{pmatrix}
0\\
-\chi_{UV}f_U/\Delta_{UV}\\
(a_V-\chi_{UV})f_U/\Delta_{UV}
\end{pmatrix}.
}
\]

This gives the physical interpretation immediately.

- The transfer leg changes `\mu_W` directly.
- The dressing activation changes `K_\eta^{(\rm eff)}` directly and shifts `\mu_W` by the same amount with opposite sign.
- The `T_U` slot remains untouched on the rigid-mouth dependent plane.

The exact `U/V` drain term is

\[
\boxed{
\mathcal D_{UV}:=\chi_{UV}UV
=
\frac{\chi_{UV}^2\,a_V\,f_U^2}{\Delta_{UV}^2}
\ge 0.
}
\]

So even when `\chi_{UV}` is negative and `V` is opposite in sign to `U`, the reduced drain is still positive.
This is the exact version of the Session-I statement that a relaxed mouth can drain positive energy into the dressing leg.

---

## 6. Weak-axisymmetric first-order packet and prefactor shift

The compact program already gives the exact first-order identities

\[
\delta\ln\mathcal T^2=\Xi_1,
\qquad
\delta\ln(1-\epsilon_\eta)=\mathcal R_1+\Xi_1,
\qquad
\delta\ln R_{\rm target}=\mathcal R_1.
\]

To connect the finite Stage-228 packet to that weak-axisymmetric front end, set

\[
U=\varepsilon u_1+O(\varepsilon^2),
\qquad
V=\varepsilon v_1+O(\varepsilon^2).
\]

Then the exact finite compiler gives

\[
\delta\ln\mathcal T^2=u_1,
\]

\[
\delta\ln(1-\epsilon_\eta)
=
-\frac{\epsilon_{\eta,\rm ref}}{1-\epsilon_{\eta,\rm ref}}\,v_1,
\]

\[
\delta\ln R_{\rm target}
=
-u_1
-\frac{\epsilon_{\eta,\rm ref}}{1-\epsilon_{\eta,\rm ref}}\,v_1.
\]

So the non-rigid first-order packet is

\[
\boxed{
\Xi_1^{\rm nr}=u_1,
}
\]

\[
\boxed{
\mathcal R_1^{\rm nr}+\Xi_1^{\rm nr}
=
-\frac{\epsilon_{\eta,\rm ref}}{1-\epsilon_{\eta,\rm ref}}\,v_1,
}
\]

\[
\boxed{
\mathcal R_1^{\rm nr}
=
-u_1
-\frac{\epsilon_{\eta,\rm ref}}{1-\epsilon_{\eta,\rm ref}}\,v_1.
}
\]

This is the key first-order compiler of the stage.

It says:

1. the direct prefactor-shape scalar is still carried by `U`,
2. the dressing leg enters only through the selected-branch residual at first order,
3. and once `V` is activated, the target-ratio packet is no longer determined by `\Xi_1` alone.

So the compact front-end packet `(\Delta_{\rm norm},\Xi_1)` is the rigid-mouth projection of a larger two-leg packet on the non-rigid branch.

---

## 7. Orbit-side / support-blind split

Stage 227 attached the leakage/work compiler to the **selected-support** packet.
Stage 228 is the complementary statement for the non-rigid packet.

Within the present closure, the quantities

\[
U,\qquad
V,\qquad
\epsilon_\eta,\qquad
\frac{R_{\rm target}}{R_{\rm target,ref}},
\qquad
\mathbf y_{\rm nr}^{\rm dep},
\qquad
\mathcal D_{UV},
\]

depend only on

\[
a_U,\qquad
a_V,\qquad
\chi_{UV},\qquad
f_U,\qquad
\epsilon_{\eta,\rm ref},
\]

and **not** on the selected-support coordinates from Stage 225.

So if `(\Lambda,\varrho)` denote the support-side placement variables, then

\[
\partial_\Lambda U
=
\partial_\varrho U
=
\partial_\Lambda V
=
\partial_\varrho V
=
0,
\]

and similarly for the finite observable compiler above, so long as the orbit-side forcing `f_U` and signed mouth coupling `\chi_{UV}` are carried independently.

That is the right structural split:

- Stage 227 is the support-side open-system compiler,
- Stage 228 is the orbit-side non-rigid packet compiler.

---

## 8. Session-I readback

The session write-up records the representative relaxed-rigid-mouth point

\[
U(r_{\rm eval})\approx 0.14313458,
\qquad
V(r_{\rm eval})\approx -0.03619791,
\qquad
\epsilon_{\eta,\rm ref}=0.3.
\]

Feeding those numbers into the exact Stage-228 compiler gives

\[
\epsilon_\eta(r_{\rm eval})
=
0.3\,e^{-0.03619791}
\approx 0.28933482,
\]

which matches the Session-I value.

The same finite compiler gives

\[
\frac{R_{\rm target}(r_{\rm eval})}{R_{\rm target,ref}}
\approx 0.87984149.
\]

And the dependent microscopic correction becomes

\[
\mathbf y_{\rm nr}^{\rm dep}(r_{\rm eval})
=
\begin{pmatrix}
0\\
0.03619791\\
0.17933249
\end{pmatrix}.
\]

So the sampled branch is a clear transfer-up / dressing-down branch:

- transfer-shape grows,
- the dressing variable falls,
- and the target ratio shifts accordingly through the exact selected-branch identity.

At first weak-axisymmetric order, the same point gives

\[
\Xi_1^{\rm nr}\approx 0.14313458,
\qquad
\mathcal R_1^{\rm nr}\approx -0.12762119
\quad
(\epsilon_{\eta,\rm ref}=0.3).
\]

So the session’s `\Xi_1` readback remains the direct transfer scalar, while the newly active dressing leg contributes the additional selected-branch residual.

---

## 9. What this stage accomplishes

Stage 228 closes the main gap left after Stage 226.

1. It turns the abstract Stage-226 `U/V` free-energy lane into an **exact finite compiler** for the physical observables `(\mathcal T^2,\epsilon_\eta,R_{\rm target})`.
2. It shows that the dressing variable `\epsilon_\eta` is no longer a frozen bystander once rigid-mouth factorization is broken.
3. It proves that the induced `U/V` drain is always nonnegative on the admissible branch.
4. It identifies the exact microscopic correction carried by the finite non-rigid packet.
5. It shows that, at first weak-axisymmetric order, the direct prefactor scalar is still `U`, but the active dressing leg adds a second selected-branch residual.
6. It cleanly separates the non-rigid packet from the support-side lane compiled in Stage 227.

So after Stage 228 the barrier-session `U/V` story is no longer floating outside the derivation stack.
It is an explicit orbit-side packet compiler.

---

## 10. Immediate next step

The next stage should now use the same non-rigid packet but relax the **source profile** itself beyond the positive Family-1 corridor.

That means:

1. promote the mouth/core source to the first compensated multimode branch,
2. compile its mouth-bias and shell-loading functionals,
3. and then feed those functionals together with the Stage-227 leakage/work lane and the Stage-228 non-rigid packet into the first honest reduced stationary barrier compiler.

That is the natural continuation to Stage 229.
