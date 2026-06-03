# Moving-Throat PDE — Stage 236: Rigid-Mouth Microscopic Dependent-Plane Projectors, the Equal-Drift Dressing Ray, and the Static-Only Restoration Gap

## Status

**Exact within the carried Stage-252 rigid-mouth packet split and the later microscopic quotient section on the dependent triple.**

This stage does **not** yet compute the actual magnitude of the rigid-mouth dressing coordinate on the completed PDE branch.
It does something narrower and sharper:

> it identifies the exact **microscopic carrier** of the static-blind residue `q_eta` once the first static same-charge strip has already been cleared.

The main new result is stricter than the Stage-252 direct-space wording:

> on the rigid-mouth slice the surviving dressing residue is not a generic three-coordinate failure.
> It lives on the exact diagonal ray
> \[
> (\Delta_T,\Delta_{K_\eta},\Delta_\mu)
> =
> -q_\eta\,(0,1,1),
> \]
> inside the dependent triple.
> So after the static gate is cleared, the unresolved same-charge orbit defect is exactly an **equal-drift `K_\eta`–`\mu` dressing shift at fixed `T_U`**.

In other words, Stage 252 identified the static-blind direct-space line.
Stage 253 identifies its exact microscopic image.

---

## Purpose

Stage 252 showed that on the rigid-mouth direct-observable plane
\[
(R_1,E_1):=
(\delta\ln R_{\rm target},\,\delta\ln\epsilon_\eta),
\]
the surviving quotient packet is
\[
(q_{\rm nt},q_\eta)=(\Xi_1,E_1),
\]
and that the first static same-charge gate only constrains
\[
q_{\rm nt}=\Xi_1,
\]
while the dressing coordinate
\[
q_\eta=\delta\ln\epsilon_\eta
\]
remains completely invisible to that gate.

But the later quotient-program compilers also say that every quotient failure is represented microscopically only on the dependent triple
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
via the exact dependent section.
So the right next step is now immediate:

> what is the exact microscopic dependent-coordinate image of the rigid-mouth packet \((q_{\rm nt},q_\eta)\)?

This stage answers that exactly.

---

## 1. Carry-forward rigid-mouth packet and dependent section

From Stage 252 the rigid-mouth packet map is
\[
\boxed{
\mathbf q_{\rm rm}:=
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_{\rm rm}
\begin{pmatrix}
R_1\\
E_1
\end{pmatrix},
\qquad
M_{\rm rm}=
\begin{pmatrix}
-1 & -c_\eta\\
0 & 1
\end{pmatrix},
}
\]
with
\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]
So explicitly,
\[
q_{\rm nt}=-R_1-c_\eta E_1=\Xi_1,
\qquad
q_\eta=E_1.
\]

From the exact microscopic quotient section on the dependent triple,
\[
\Delta_T^{(q)}=\frac{q_{\rm tr}}{1+\chi_{0,*}},
\qquad
\Delta_{K_\eta}^{(q)}=-q_\eta,
\qquad
\Delta_\mu^{(q)}=\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}+q_{\rm nt}-q_\eta.
\]
On the rigid-mouth slice,
\[
q_{\rm tr}=0,
\]
so the dependent correction reduces exactly to
\[
\boxed{
\mathbf y_{\rm rm}:=
\begin{pmatrix}
\Delta_T^{(q)}\\
\Delta_{K_\eta}^{(q)}\\
\Delta_\mu^{(q)}
\end{pmatrix}
=
S_{\rm rm}^{\rm dep}
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix},
\qquad
S_{\rm rm}^{\rm dep}:=
\begin{pmatrix}
0 & 0\\
0 & -1\\
1 & -1
\end{pmatrix}.
}
\]
So explicitly,
\[
\boxed{
\Delta_T^{(q)}=0,
\qquad
\Delta_{K_\eta}^{(q)}=-q_\eta,
\qquad
\Delta_\mu^{(q)}=q_{\rm nt}-q_\eta.
}
\]

The rigid-mouth quotient failure therefore fills the full plane
\[
\Delta_T=0
\]
inside the dependent triple.

---

## 2. Exact direct-observable-to-dependent compiler

Composing the direct rigid-mouth packet map with the dependent microscopic section gives
\[
\boxed{
\mathbf y_{\rm rm}
=
C_{\rm rm}^{\rm dep}
\begin{pmatrix}
R_1\\
E_1
\end{pmatrix},
\qquad
C_{\rm rm}^{\rm dep}:=
S_{\rm rm}^{\rm dep}M_{\rm rm}
=
\begin{pmatrix}
0 & 0\\
0 & -1\\
-1 & -\dfrac{1}{1-\epsilon_{\eta,*}}
\end{pmatrix}.
}
\]
So the exact microscopic dependent correction in direct rigid-mouth variables is
\[
\boxed{
\Delta_T^{(q)}=0,
}
\]
\[
\boxed{
\Delta_{K_\eta}^{(q)}=-E_1=-\delta\ln\epsilon_\eta,
}
\]
\[
\boxed{
\Delta_\mu^{(q)}
=
-R_1-\frac{1}{1-\epsilon_{\eta,*}}E_1
=
-\delta\ln R_{\rm target}
-\frac{1}{1-\epsilon_{\eta,*}}\delta\ln\epsilon_\eta.
}
\]

So even before any projector algebra is applied, one exact fact is already clear:

> on the rigid-mouth slice, the microscopic quotient correction never uses `T_U` at first order.
> The entire failure is carried only by the `(K_\eta,\mu)` plane.

---

## 3. Exact dependent-plane packet projectors

The rigid-mouth dependent compiler has an exact left inverse on the plane `\Delta_T=0`:
\[
\boxed{
L_{\rm rm}^{\rm dep}:=
\begin{pmatrix}
0 & -1 & 1\\
0 & -1 & 0
\end{pmatrix},
\qquad
L_{\rm rm}^{\rm dep}S_{\rm rm}^{\rm dep}=I_2.
}
\]
So the packet coordinates can be recovered directly from the dependent correction by
\[
\boxed{
q_{\rm nt}=\Delta_\mu^{(q)}-\Delta_{K_\eta}^{(q)},
\qquad
q_\eta=-\Delta_{K_\eta}^{(q)}.
}
\]

Push the packet projectors back to the dependent plane:
\[
P_{\rm nt}^{\rm dep}:=
S_{\rm rm}^{\rm dep}
\begin{pmatrix}1&0\\0&0\end{pmatrix}
L_{\rm rm}^{\rm dep},
\qquad
P_\eta^{\rm dep}:=
S_{\rm rm}^{\rm dep}
\begin{pmatrix}0&0\\0&1\end{pmatrix}
L_{\rm rm}^{\rm dep}.
\]
They are explicit:
\[
\boxed{
P_{\rm nt}^{\rm dep}=
\begin{pmatrix}
0&0&0\\
0&0&0\\
0&-1&1
\end{pmatrix},
\qquad
P_\eta^{\rm dep}=
\begin{pmatrix}
0&0&0\\
0&1&0\\
0&1&0
\end{pmatrix}.
}
\]
These are exact complementary projectors on the rigid-mouth dependent plane:
\[
(P_{\rm nt}^{\rm dep})^2=P_{\rm nt}^{\rm dep},
\qquad
(P_\eta^{\rm dep})^2=P_\eta^{\rm dep},
\qquad
P_{\rm nt}^{\rm dep}P_\eta^{\rm dep}=P_\eta^{\rm dep}P_{\rm nt}^{\rm dep}=0,
\]
and
\[
P_{\rm nt}^{\rm dep}+P_\eta^{\rm dep}
=
\begin{pmatrix}
0&0&0\\
0&1&0\\
0&0&1
\end{pmatrix},
\]
which is the identity on the plane `\Delta_T=0`.

So every rigid-mouth dependent correction decomposes uniquely as
\[
\boxed{
\mathbf y_{\rm rm}=\mathbf y_{\rm nt}+\mathbf y_\eta,
\qquad
\mathbf y_{\rm nt}:=P_{\rm nt}^{\rm dep}\mathbf y_{\rm rm},
\qquad
\mathbf y_\eta:=P_\eta^{\rm dep}\mathbf y_{\rm rm}.
}
\]
Explicitly,
\[
\boxed{
\mathbf y_{\rm nt}=
\begin{pmatrix}
0\\
0\\
q_{\rm nt}
\end{pmatrix},
\qquad
\mathbf y_\eta=
-q_\eta
\begin{pmatrix}
0\\
1\\
1
\end{pmatrix}.
}
\]

So the rigid-mouth dependent plane contains two exact, complementary, physically meaningful pieces:

- `y_nt`: the pure `\mu` mismatch seen by the first static strip,
- `y_eta`: the equal-drift dressing ray that remains after that strip is cleared.

---

## 4. The equal-drift dressing ray and the microscopic meaning of the static strip

Stage 252 showed that the static strip is
\[
q_{\rm nt}=\Xi_1=0.
\]
In the dependent microscopic plane, this becomes
\[
\boxed{
q_{\rm nt}=0
\iff
\Delta_\mu^{(q)}=\Delta_{K_\eta}^{(q)}.
}
\]
So the entire static-blind line is exactly the diagonal ray
\[
\boxed{
\mathbf y_\eta=-q_\eta(0,1,1)^T.
}
\]
Equivalently, using the direct observable compiler,
\[
q_{\rm nt}=0
\iff
R_1=-c_\eta E_1,
\]
and then
\[
\boxed{
\mathbf y_{\rm rm}
=
\begin{pmatrix}
0\\
-E_1\\
-E_1
\end{pmatrix}
=
\frac{R_1}{c_\eta}
\begin{pmatrix}
0\\
1\\
1
\end{pmatrix}.
}
\]
So the Stage-252 static-blind line in the direct observable plane maps exactly to an **equal-drift `K_\eta`–`\mu` ray** in the dependent microscopic plane.

Its microscopic norm is exact:
\[
\boxed{
\|\mathbf y_\eta\|^2 = 2q_\eta^2 = 2E_1^2 = \frac{2R_1^2}{c_\eta^2}.
}
\]
Therefore clearing the first static same-charge ceiling does **not** force the microscopic quotient correction to be small.
It only forces that correction to lie on a one-dimensional diagonal ray.

The full rigid-mouth orbit-lock point remains the endpoint of that ray:
\[
\boxed{
q_{\rm nt}=0,\ q_\eta=0
\iff
\Delta_T^{(q)}=\Delta_{K_\eta}^{(q)}=\Delta_\mu^{(q)}=0.
}
\]

---

## 5. Exact microscopic correction compilers

The dependent-plane projectors immediately give the two natural microscopic corrections.

### 5.1 Static-only microscopic correction

Remove only the static component:
\[
\boxed{
\Delta\mathbf y_{\rm static}:=-\mathbf y_{\rm nt}=
\begin{pmatrix}
0\\
0\\
-q_{\rm nt}
\end{pmatrix}.
}
\]
After this correction,
\[
\mathbf y_{\rm rm}+\Delta\mathbf y_{\rm static}=\mathbf y_\eta,
\qquad
q_{\rm nt}\to 0,
\qquad
q_\eta\to q_\eta.
\]
So the first static ceiling is cleared by changing only `\mu_W` inside the dependent triple.

### 5.2 Full orbit-lock correction

Remove the entire rigid-mouth dependent correction:
\[
\boxed{
\Delta\mathbf y_{\rm orbit}:=-\mathbf y_{\rm rm}.
}
\]
Equivalently,
\[
\boxed{
\Delta\mathbf y_{\rm orbit}
=
\Delta\mathbf y_{\rm static}
+
q_\eta
\begin{pmatrix}
0\\
1\\
1
\end{pmatrix}.
}
\]
So the extra step beyond the static gate is again completely sharp:

> once the `\mu` mismatch `q_{\rm nt}` has been removed, the only remaining orbit-restoring correction is the equal `K_\eta`–`\mu` dressing shift generated by `q_\eta`.

This is the exact rigid-mouth meaning of the static-only restoration gap.

---

## 6. What changes physically after Stage 253

Stage 252 already said that the first static same-charge strip is not the whole rigid-mouth orbit-lock problem.
Stage 253 now says something more microscopic and more useful:

1. on the rigid-mouth slice, the quotient-failure image is the full plane `\Delta_T=0`,
2. the first static gate tests only the `\mu-K_\eta` difference,
3. the static-blind residue is exactly the diagonal equal-drift ray `\Delta_\mu=\Delta_{K_\eta}`,
4. and the surviving same-charge obstruction after the static strip is cleared is therefore not generic throat motion but one scalar `K_\eta`–`\mu` dressing amplitude.

So the next honest theorem gate is now even sharper than at Stage 252:

> compute the actual dressing coordinate `q_\eta=\delta\ln\epsilon_\eta`, because after the first static gate is cleared that single scalar is exactly the amplitude of the remaining equal-drift microscopic obstruction.

---

## 7. Best current verdict after Stage 253

The same-charge corridor is still alive, but the rigid-mouth bottleneck has narrowed again.

It is no longer enough to know that the branch lies inside the Stage-252 static strip.
On the rigid-mouth slice that strip only kills the `q_{\rm nt}` packet component.
The surviving dressing coordinate
\[
q_\eta=\delta\ln\epsilon_\eta
\]
maps microscopically to the exact equal-drift ray
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
=
-q_\eta(0,1,1),
\]
which remains completely unconstrained by the first static ceiling except for its one-dimensional direction.

So the sharp rigid-mouth statement is now:

> the first static same-charge ceiling is a necessary codimension-one condition, but the true rigid-mouth orbit-lock problem is still codimension two, and the entire post-static microscopic obstruction is the single dressing coordinate `q_\eta` carried by an equal-drift `K_\eta`–`\mu` ray at fixed `T_U`.

That is the cleanest continuation point into Stage 237.

---

## 8. SymPy-backed status

The accompanying SymPy audit verifies:

- the rigid-mouth dependent section obtained from the later microscopic quotient formulas,
- the exact direct-observable-to-dependent compiler `C_{\rm rm}^{\rm dep}`,
- the exact left inverse `L_{\rm rm}^{\rm dep}`,
- the complementary dependent-plane projectors `P_{\rm nt}^{\rm dep}` and `P_\eta^{\rm dep}`,
- the decomposition `y_{\rm rm}=y_{\rm nt}+y_\eta`,
- the static-strip equivalence `q_{\rm nt}=0 \iff \Delta_\mu^{(q)}=\Delta_{K_\eta}^{(q)}`,
- the equal-drift dressing ray and its exact norm,
- and the static-only and full orbit-lock correction compilers.

Supporting file:
- `moving_throat_pde_stage236_rigid_mouth_microscopic_dependent_plane_projectors_equal_drift_dressing_ray_and_static_only_restoration_gap_sympy_audit.py`
