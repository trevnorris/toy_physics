
# Same-Charge Barrier Audit — Stage 017: Direct Branch-Observable Static Gate and the Two-Observable Kill Test

## Status

**Exact within the carried Stage-016 same-charge reduction and the later `5PN` direct-branch observable compiler.**

This stage does not introduce a new physical mechanism.
It rewrites the unresolved rigid-mouth orbit-lock problem in the sharpest actual-branch variables currently available:

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta.
\]

The main new result is stronger than the Stage-016 wording:

> at first weak-axisymmetric order, the surviving static same-charge gate is already independent of the tracking observable \(R_{\rm tr}\). Once the coherent branch is written in direct observables, the remaining static bottleneck lives entirely on the two-observable plane \((R_{\rm target},\epsilon_\eta)\).

So the “rigid mouth” picture remains useful as an interpretation, but the exact reduced kill test is now a **target–dressing mismatch theorem**, not a direct mouth-motion theorem.

---

## 0. Why this stage is needed

Stage 016 used the Stage-171 observable compiler to show that on a track-locked branch
\[
\delta\ln R_{\rm tr}=0
\quad\Longrightarrow\quad
\Xi_1=\delta\ln \mathfrak N_*,
\qquad
\mathfrak N_*=\mathcal T^2 R_{\rm tr}^{B_*}.
\]

That was already a useful rigid-mouth translation. But the later `5PN` notes go further.
They show that the coherent weak-axisymmetric orbit packet can be charted directly by the exact branch observables
\[
(R_{\rm tr},R_{\rm target},\epsilon_\eta),
\]
and that in this chart the first-order defect packet becomes completely triangular.

So the next honest step is to eliminate the remaining quotient bookkeeping and write the same-charge static gate directly in those observables.

---

## 1. Exact finite quotient chart from direct branch observables

Relative to a coherent reference branch
\[
(R_{\rm tr,ref},R_{\rm target,ref},\epsilon_{\eta,\rm ref}),
\]
the exact finite quotient coordinates are

\[
q_{\rm tr}
=
-\,C_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right),
\]

\[
q_{\rm nt}
=
B_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right)
+
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\]

\[
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\]

with the carried branch constants

\[
C_*=
\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
{\chi_{0,*}\delta_{U,*}},
\qquad
B_*=
\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
\]

The inverse map is exact:

\[
R_{\rm tr}
=
R_{\rm tr,ref}\,e^{-q_{\rm tr}/C_*},
\]

\[
\epsilon_\eta
=
\epsilon_{\eta,\rm ref}\,e^{q_\eta},
\]

\[
R_{\rm target}
=
R_{\rm target,ref}\,
e^{-q_{\rm nt}-(B_*/C_*)q_{\rm tr}}
\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}.
\]

So once the actual weak-axisymmetric branch is known, it can be tested either by the finite quotient packet \((q_{\rm tr},q_{\rm nt},q_\eta)\) or directly by the three physical observables \((R_{\rm tr},R_{\rm target},\epsilon_\eta)\).

---

## 2. Exact first-order defect compiler in direct branch language

Linearizing the finite chart gives

\[
q_{\rm tr}
=
-\,C_*\,\delta\ln R_{\rm tr},
\]

\[
q_{\rm nt}
=
B_*\,\delta\ln R_{\rm tr}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln\epsilon_\eta
-
\delta\ln R_{\rm target},
\]

\[
q_\eta
=
\delta\ln\epsilon_\eta.
\]

Composing this with the exact quotient-to-defect compiler yields the triangular first-order direct-branch map

\[
\boxed{\Theta_1=\delta\ln R_{\rm tr},}
\]

\[
\boxed{
\Xi_1
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln\epsilon_\eta,
}
\]

\[
\boxed{\mathcal R_1=\delta\ln R_{\rm target}.}
\]

So the physical first-order defect packet is exactly equivalent to the three direct observable drifts
\[
(\delta\ln R_{\rm tr},\delta\ln R_{\rm target},\delta\ln\epsilon_\eta).
\]

The inverse map is equally simple:

\[
\delta\ln R_{\rm tr} = \Theta_1,
\]

\[
\delta\ln R_{\rm target} = \mathcal R_1,
\]

\[
\delta\ln \epsilon_\eta
=
-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\,
(\mathcal R_1+\Xi_1).
\]

---

## 3. Two exact cancellations

Define the selected-branch dressing coefficient

\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]

Then
\[
\Xi_1 = -\,\delta\ln R_{\rm target} - c_\eta\,\delta\ln\epsilon_\eta.
\]

This yields two exact structural cancellations.

### 3.1 `R_tr` cancels out of `Xi_1`

Although the finite nontracking quotient contains the tracking observable,
\[
q_{\rm nt}
=
B_*\,\delta\ln R_{\rm tr}
-
c_\eta\,\delta\ln\epsilon_\eta
-
\delta\ln R_{\rm target},
\]
the physical defect is
\[
\Xi_1
=
q_{\rm nt}
+
\frac{B_*}{C_*}\,q_{\rm tr}.
\]

So the tracking contribution cancels exactly, and
\[
\boxed{\partial_{\delta\ln R_{\rm tr}}\Xi_1=0.}
\]

This is stronger than the Stage-016 wording.
It means the surviving first-order static same-charge gate is **already mouth-blind** in the direct observable chart.

### 3.2 Support/source rescaling cancels out of the entire first-order defect packet

The later `5PN` notes also show that the direct branch observables \((R_{\rm tr},R_{\rm target},\epsilon_\eta)\), and therefore the defect packet \((\Theta_1,\Xi_1,\mathcal R_1)\), are exactly blind to the total support-compensation baseline and to the coherent support-enhancement ratio:
\[
\partial_{M_{\rm tr}}(\Theta_1,\Xi_1,\mathcal R_1)=0,
\qquad
\partial_\zeta(\Theta_1,\Xi_1,\mathcal R_1)=0.
\]

So the Stage-015 support/source safety verdict and the Stage-016 static gate really are distinct theorem pieces.

---

## 4. Exact rigid-mouth specialization revisited

The rigid-mouth interpretation still has one exact role:
\[
\delta\ln R_{\rm tr}=0
\quad\Longleftrightarrow\quad
\Theta_1=0.
\]

On that branch,
\[
q_{\rm tr}=0,
\qquad
q_{\rm nt}
=
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right).
\]

And at first order,
\[
\boxed{q_{\rm nt}=\Xi_1.}
\]

So the Stage-016 statement
\[
\Xi_1=\delta\ln\mathfrak N_*
\]
on a track-locked branch sharpens here to:

> under rigid-mouth lock, the surviving first-order defect is exactly the tangent of the finite **two-observable** nontracking quotient built from \(R_{\rm target}\) and \(\epsilon_\eta\).

That is the cleanest current form of the “inside can still repackage while the mouth stays locked” picture.

---

## 5. The direct two-observable static gate

Let
\[
R_1:=\delta\ln R_{\rm target},
\qquad
E_1:=\delta\ln\epsilon_\eta.
\]

Then the entire first-order same-charge placement scalar is

\[
\boxed{\Xi_1 = -R_1 - c_\eta E_1.}
\]

So the Stage-007 / Stage-016 static ceilings become exact two-observable band conditions.

### 5.1 Robust static gate

\[
\boxed{
|\epsilon(-R_1-c_\eta E_1)|
\le
0.367930328492646.
}
\]

Equivalently,
\[
R_1
\in
\left[
-c_\eta E_1-\frac{0.367930328492646}{|\epsilon|},
\,
-c_\eta E_1+\frac{0.367930328492646}{|\epsilon|}
\right].
\]

### 5.2 Nonempty-corridor gate

\[
\boxed{
|\epsilon(-R_1-c_\eta E_1)|
\le
0.737619063660757.
}
\]

Equivalently,
\[
R_1
\in
\left[
-c_\eta E_1-\frac{0.737619063660757}{|\epsilon|},
\,
-c_\eta E_1+\frac{0.737619063660757}{|\epsilon|}
\right].
\]

So the first unresolved same-charge kill test is now literally a strip in the \((R_1,E_1)\) plane.

---

## 6. Canonical direct-branch families

For a prescribed first-order defect value
\[
\Xi_1=\xi,
\]
the direct branch families are simple.

### 6.1 Pure target-drift family

\[
E_1=0,
\qquad
R_1=-\xi.
\]

So if dressing is frozen, the whole defect is a pure selected-target drift.

### 6.2 Pure dressing-drift family

\[
R_1=0,
\qquad
E_1=-\frac{\xi}{c_\eta}.
\]

So if the target is frozen, the whole defect is a pure dressing drift.

### 6.3 Balanced minimal-norm family

Minimizing
\[
R_1^2+E_1^2
\]
subject to
\[
-R_1-c_\eta E_1=\xi
\]
gives the exact minimum-norm branch

\[
\boxed{
R_1=-\frac{\xi}{1+c_\eta^2},
\qquad
E_1=-\frac{c_\eta\,\xi}{1+c_\eta^2}.
}
\]

So even before the full PDE branch is known, the direct observable chart already gives a canonical least-deformation family that realizes any chosen \(\Xi_1\).

---

## 7. What changes in the physical interpretation

The user’s “hyper-trumpet / internal choke” picture is still a useful heuristic, but Stage 17 sharpens it in an important way.

Stage 16 suggested:
- rigid mouth \(\Rightarrow\) internal transfer observable \(\mathfrak N_*\) carries the first unresolved scalar.

Stage 17 now shows:
- at first order, the same unresolved scalar is already
  \[
  \Xi_1=-\delta\ln R_{\rm target}-c_\eta\,\delta\ln\epsilon_\eta,
  \]
  so it is **independent of the mouth tracking variable \(R_{\rm tr}\)**.

That means the first static failure mode is not a mouth-motion theorem even in disguised form.
It is a **selected-target / dressing mismatch theorem**.

A “hyper-trumpet choke” reading is still consistent with the math if the physical reason \(R_{\rm target}\) or \(\epsilon_\eta\) drifts is internal throat repackaging.
But the reduced observable that actually decides the branch is not the mouth variable. It is the two-observable combination above.

---

## 8. Best current verdict after Stage 017

The same-charge corridor is still alive, but the unresolved static bottleneck has narrowed again.

It is no longer:

- support/source strength,
- dynamic wall-window survival,
- or even mouth tracking at first order.

It is now exactly:

\[
\boxed{
\Xi_1
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln\epsilon_\eta,
}
\]
together with the Stage-016 robust / nonempty ceilings.

So the next honest theorem gate is:

> compute the actual coherent weak-axisymmetric branch drifts of \(R_{\rm target}\) and \(\epsilon_\eta\), then test whether their exact linear combination clears the direct static gate above.

That is the sharpest branch-level same-charge kill test reached so far.
