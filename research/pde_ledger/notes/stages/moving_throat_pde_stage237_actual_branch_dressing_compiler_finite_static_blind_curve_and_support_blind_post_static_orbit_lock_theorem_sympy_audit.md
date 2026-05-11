# Moving-Throat PDE — Stage 237: Actual-Branch Dressing Compiler, the Finite Static-Blind Curve, and the Support-Blind Post-Static Orbit-Lock Theorem

## Status

**Exact within the carried Stage-253 rigid-mouth dependent-plane split and the later direct-branch / microscopic packet compilers on the actual coherent branch.**

This stage does **not** solve the full moving-throat PDE branch.
It does something narrower and more useful:

> it computes the surviving rigid-mouth dressing coordinate `q_eta` exactly in the actual branch observables and in the actual microscopic variables, and shows that this post-static obstruction is completely blind to the coherent support-enhancement sector.

So after Stage 237, the next same-charge gate is no longer vague.
Once the first static ceiling `q_nt = 0` has been cleared, the remaining obstruction is exactly the single actual-branch scalar
\[
q_\eta = \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\mathrm{ref}}}\right),
\]
which may also be read directly from the target observable `R_{\rm target}` on the static-blind curve.

---

## Purpose

Stage 253 proved that on the rigid-mouth slice the remaining quotient failure is carried microscopically by the equal-drift dependent-plane ray
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
=
-q_\eta(0,1,1).
\]
So after the first static gate is cleared, the entire unresolved same-charge obstruction is one scalar: the dressing coordinate `q_eta`.

But Stage 253 still left one obvious unfinished task:

> how do we compute `q_eta` on the actual coherent branch, rather than only locating its microscopic carrier?

The later branch-packet compilers already contain the answer. On the actual coherent branch the finite quotient packet is charted exactly by the three physical observables
\[
(R_{\rm tr},\,R_{\rm target},\,\epsilon_\eta),
\]
and the dressing coordinate is simply the logarithmic ratio of `\epsilon_\eta` itself. This stage isolates that fact, pushes it through the rigid-mouth static gate, and shows that the surviving dressing obstruction is rigorously support-blind.

---

## 1. Exact rigid-mouth finite packet on the actual branch

Relative to a coherent reference branch
\[
(R_{\rm tr,ref},\,R_{\rm target,ref},\,\epsilon_{\eta,\rm ref}),
\]
the exact finite quotient coordinates are
\[
q_{\rm tr}
=
-C_* \ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right),
\]
\[
q_{\rm nt}
=
B_* \ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right)
+\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\]
\[
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]

Here `B_*` and `C_*` are the carried direct-branch packet coefficients from the earlier observable compiler.
On the rigid-mouth slice,
\[
q_{\rm tr}=0,
\qquad
R_{\rm tr}=R_{\rm tr,ref},
\]
so the surviving packet is exactly
\[
\boxed{
q_{\rm nt}
=
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\qquad
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]

So on the actual rigid-mouth branch the two finite surviving coordinates are already completely explicit:

- `q_nt` is the target–dressing mismatch of the finite branch observables,
- `q_eta` is the pure dressing logarithmic ratio.

---

## 2. Exact finite static-blind curve and direct computation of `q_eta`

The first static same-charge ceiling is
\[
q_{\rm nt}=0.
\]
On the rigid-mouth branch this becomes the exact finite relation
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}.
}
\]

So the static-blind set is not just a tangent line.
It is an exact one-parameter finite curve in the direct branch observables.

Parameterizing that curve by `q_eta`, we use
\[
\epsilon_\eta=\epsilon_{\eta,\rm ref}e^{q_\eta},
\]
which gives
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^{q_\eta}}{1-\epsilon_{\eta,\rm ref}}.
}
\]

Conversely, once the static gate has been cleared, `q_eta` can be computed directly from the actual branch target observable alone:
\[
\boxed{
q_\eta
=
\ln\!\left(
\frac{1-(1-\epsilon_{\eta,\rm ref})\,R_{\rm target}/R_{\rm target,ref}}
{\epsilon_{\eta,\rm ref}}
\right).
}
\]

So after the static gate is passed, the direct same-charge branch is exactly one-dimensional and can be charted either by

- the dressing observable `\epsilon_\eta`, or
- the target observable `R_{\rm target}`.

The full rigid-mouth orbit-lock point remains the endpoint of that curve:
\[
q_{\rm nt}=0,\ q_\eta=0
\iff
R_{\rm target}=R_{\rm target,ref},
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
\]

---

## 3. Exact first-order compiler and the tangent of the finite static-blind curve

Linearizing around the coherent reference branch gives
\[
q_\eta=\delta\ln\epsilon_\eta,
\]
and on the rigid-mouth slice
\[
q_{\rm nt}
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta.
\]

Define the carried dressing coefficient
\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]
Then the exact first-order direct packet is
\[
\boxed{
q_{\rm nt}=-R_1-c_\eta E_1,
\qquad
q_\eta=E_1,
}
\]
with
\[
R_1:=\delta\ln R_{\rm target},
\qquad
E_1:=\delta\ln\epsilon_\eta.
\]

So once the static gate is cleared,
\[
q_{\rm nt}=0
\iff
R_1=-c_\eta q_\eta,
\qquad
q_\eta=-\frac{R_1}{c_\eta}.
\]

This is exactly the tangent relation of the finite curve from Section 2.
Indeed,
\[
\frac{d}{dq_\eta}
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right)\Bigg|_{q_\eta=0}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
=
-c_\eta.
\]

So the earlier rigid-mouth static-blind line is precisely the tangent of the finite actual-branch static-blind curve.

---

## 4. Exact microscopic compiler for `q_eta`

### 4.1 Direct microscopic coherent variables

On the actual coherent branch,
\[
\epsilon_\eta
=
\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
so the finite dressing coordinate is
\[
\boxed{
q_\eta
=
\ln\!\left(
\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}
\frac{K_{U,\rm ref}K_{\eta,\rm ref}^{(\mathrm{eff})}}{c_{\eta U,\rm ref}^2}
\right).
}
\]

Equivalently,
\[
\boxed{
q_\eta
=
2\ln\!\left(\frac{c_{\eta U}}{c_{\eta U,\rm ref}}\right)
-
\ln\!\left(\frac{K_U}{K_{U,\rm ref}}\right)
-
\ln\!\left(\frac{K_\eta^{(\mathrm{eff})}}{K_{\eta,\rm ref}^{(\mathrm{eff})}}\right).
}
\]

### 4.2 First-order microscopic drift packet

For the weak-axisymmetric microscopic drifts
\[
(c_1,\kappa_U,\kappa_\eta),
\]
the exact first-order extractor is
\[
\boxed{
q_\eta
=
\delta\ln\epsilon_\eta
=
2c_1-\kappa_U-\kappa_\eta.
}
\]

So the post-static same-charge obstruction is directly the combined logarithmic drift of

- the wall–U coupling `c_{ηU}`,
- the `U` stiffness `K_U`,
- and the effective wall stiffness `K_η^{(\mathrm{eff})}`.

No additional transport bookkeeping is needed.

---

## 5. Exact support-blindness theorem for the dressing coordinate

Introduce the coherent support ratio in the usual form
\[
\zeta
=
\frac{\lambda_\phi^2 K_W^{(\mathrm{eff})}}{\lambda_W^2 K_\phi^{(\mathrm{eff})}},
\]
and the total coherent support-enhanced baseline
\[
M_{\rm tr}
=
M_{\rm mix}
\left[
1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}
\right].
\]

Then the exact dressing coordinate satisfies
\[
\boxed{
\partial_\zeta q_\eta = 0,
\qquad
\partial_{M_{\rm tr}} q_\eta = 0.
}
\]

At the microscopic parameter level this is even sharper:
\[
\boxed{
\partial_{\lambda_\phi} q_\eta = 0,
\qquad
\partial_{K_\phi^{(\mathrm{eff})}} q_\eta = 0.
}
\]

So the coherent support lane is exactly blind to the dressing obstruction.

This has a strong physical consequence:

> support enhancement can rescue the steady normalization side of the branch, but it cannot change the post-static same-charge dressing obstruction.

So if the static gate has been cleared and `q_eta` is still nonzero, no adjustment of the coherent support ratio can remove that failure.

---

## 6. Post-static orbit-lock theorem on the actual rigid-mouth branch

Stage 253 showed that after the static gate is cleared, the remaining microscopic correction is the equal-drift dependent-plane ray
\[
\mathbf y_\eta = -q_\eta(0,1,1)^T.
\]

Stage 237 now identifies the amplitude of that ray exactly:
\[
\boxed{
\mathbf y_\eta
=
-\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right)
\begin{pmatrix}
0\\
1\\
1
\end{pmatrix}.
}
\]
Equivalently, on the finite static-blind curve,
\[
\boxed{
\mathbf y_\eta
=
-\ln\!\left(
\frac{1-(1-\epsilon_{\eta,\rm ref})\,R_{\rm target}/R_{\rm target,ref}}
{\epsilon_{\eta,\rm ref}}
\right)
\begin{pmatrix}
0\\
1\\
1
\end{pmatrix}.
}
\]

Therefore the full rigid-mouth post-static orbit-lock theorem is
\[
\boxed{
q_{\rm nt}=0\ \text{and}\ q_\eta=0
\iff
R_{\rm target}=R_{\rm target,ref}
\ \text{and}\
\epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]

At first order this becomes
\[
\boxed{
q_{\rm nt}=0\ \text{and}\ q_\eta=0
\iff
R_1=0
\ \text{and}\
E_1=0.
}
\]

So after the first static ceiling is cleared, the same-charge barrier reduces to one exact test:

> is `\epsilon_\eta` invariant on the actual rigid-mouth branch?

If yes, the post-static dressing ray collapses and orbit lock is restored.
If no, the same-charge obstruction remains, and its exact amplitude is `q_eta`.

---

## 7. Best current summary after Stage 237

The continuation from Stage 253 is now complete.

1. The static-blind microscopic carrier is still the equal-drift dependent-plane ray
   \[
   -q_\eta(0,1,1).
   \]
2. Its exact actual-branch amplitude is now explicit:
   \[
   q_\eta = \ln(\epsilon_\eta/\epsilon_{\eta,\rm ref}).
   \]
3. Once the first static gate is cleared, that same amplitude can be read directly from `R_{\rm target}`:
   \[
   q_\eta
   =
   \ln\!\left(
   \frac{1-(1-\epsilon_{\eta,\rm ref})R_{\rm target}/R_{\rm target,ref}}
   {\epsilon_{\eta,\rm ref}}
   \right).
   \]
4. The dressing coordinate is exactly support-blind.
5. Therefore the next actual same-charge theorem gate is no longer the full quotient packet.
   It is simply:

> compute `\epsilon_\eta` on the actual rigid-mouth branch and test whether it is invariant.

That is the sharpest post-static same-charge criterion reached so far.

---

## 8. SymPy-backed status

The accompanying SymPy audit verifies:

- the carried finite packet formulas and their rigid-mouth reduction,
- the exact finite static-blind curve and its inverse `q_eta(R_{\rm target})`,
- the first-order packet compiler and the tangent slope `-c_\eta`,
- the exact microscopic coherent compiler for `q_eta`,
- the first-order microscopic drift extractor `q_\eta = 2c_1-\kappa_U-\kappa_\eta`,
- the exact support-blindness identities,
- the post-static equal-drift dependent-plane ray with its actual-branch amplitude,
- and the codimension-two rigid-mouth orbit-lock theorem in both finite and first-order form.

Supporting file:
- `moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_sympy_audit.py`
