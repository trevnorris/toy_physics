# Moving-Throat PDE — Stage 188: Branch-Observable Completion and the Exact First-Order Observable Compiler

## Status

**Exact within the coherent local D/N weak-axisymmetric reference-branch closure.**

This stage does **not** replace the exact finite quotient closure of Stage 187.
It identifies the exact **first-order PDE-facing observable packet** whose vanishing is equivalent to vanishing first grouped weak-axisymmetric defect on the selected coherent branch.

---

## Purpose

Stage 187 finished the coherent weak-axisymmetric problem at the exact finite quotient level:

a) the similarity orbit `\(\mathcal G_*\)` is exact,

b) the finite quotient coordinates are the three direct microscopic invariants
\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta,
\]

and

c) the first grouped weak-axisymmetric defect is the infinitesimal motion of the actual branch in that exact three-dimensional quotient.

That is the right reduced finish for the orbit side, but it is not yet the most useful **PDE-facing** statement.
The actual moving-throat solve will not first return abstract quotient coordinates. It will return direct branch observables.

The missing bridge is therefore:

> convert the Stage 187 tangent quotient packet into a direct observable packet built from the exact branch quantities
> \[
> R_{\rm tr},
> \qquad
> \mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
> \qquad
> \epsilon_\eta,
> \]
> and then compile that observable packet directly into the defect triple
> \[
> (\Theta_1,\Xi_1,\mathcal R_1).
> \]

That is what this stage does.

The main outputs are:

1. an exact coefficient identity
   \[
   A_{{\rm tr},*}=B_* C_{{\rm tr},*},
   \]
2. an exact first-order compiler from the observable packet
   \[
   \Delta_{\rm obs}^{(1)}:=
   \begin{pmatrix}
   \delta\ln R_{\rm tr}\\
   \delta\ln \mathfrak N_*\\
   \delta\ln \epsilon_\eta
   \end{pmatrix}
   \]
   to the Stage 187 tangent quotient packet,
3. an exact first-order compiler from `\(\Delta_{\rm obs}^{(1)}\)` to the defect packet
   \[
   \Delta_{\rm def}^{(1)}:=(\Theta_1,\Xi_1,\mathcal R_1)^T,
   \]
4. and the sharp branch-observable zero-defect theorem
   \[
   \Theta_1=\Xi_1=\mathcal R_1=0
   \iff
   \delta\ln R_{\rm tr}=0,
   \quad
   \delta\ln \mathfrak N_*=0,
   \quad
   \delta\ln\epsilon_\eta=0.
   \]

So this stage is the clean PDE-facing completion of the Stage 187 quotient theorem.

---

## 1. Carry-forward exact branch observables

### 1.1 Tracking factor

The coherent tracking factor is the exact branch observable
\[
\boxed{
R_{\rm tr}
=
\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)}.
}
\]
By definition of the grouped weak-axisymmetric tracking defect,
\[
\boxed{
\delta\ln R_{\rm tr}=\Theta_1.
}
\]

### 1.2 Corrected nontracking transfer-shape observable

The direct grouped weak-axisymmetric defect `\(\Xi_1\)` is the logarithmic drift of the effective transfer shape,
\[
\boxed{
\delta\ln \mathcal T^2=\Xi_1,
\qquad
\mathcal T^2=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
}
\]
But the genuine nontracking coordinate is the tracking-feed-through-subtracted composite
\[
\boxed{
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
}
\]
Its first grouped weak-axisymmetric logarithmic drift is
\[
\boxed{
\delta\ln \mathfrak N_*=\Sigma_{\rm nt}.
}
\]

### 1.3 Dressing observable and complementary selected-branch observable

The third direct branch observable is already microscopic:
\[
\boxed{
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
\]
The exact selected-branch identity
\[
\boxed{
\frac{R_{\rm target}\mathcal T^2}{\Lambda_0}=1-\epsilon_\eta,
\qquad
\Lambda_0:=\frac{27\pi^2Gc_s^5}{20a^5c^5},
}
\]
shows that the complementary selected-branch observable
\[
\boxed{
\mathfrak E:=1-\epsilon_\eta
}
\]
has logarithmic drift
\[
\boxed{
\delta\ln\mathfrak E
=
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta.
}
\]
So the selected-branch residual is not an extra independent datum. It is the complementary drift of the same dressing variable.

---

## 2. Exact coefficient identity

Carry forward the exact coherent reference-branch coefficients
\[
\boxed{
C_{{\rm tr},*}
=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
}
\]
\[
\boxed{
A_{{\rm tr},*}
=
\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}.
}
\]
Then with
\[
\boxed{
B_*:=\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}},
}
\]
we have the exact identity
\[
\boxed{
A_{{\rm tr},*}=B_* C_{{\rm tr},*}.
}
\]

This identity is the algebraic hinge that lets the corrected nontracking composite `\(\mathfrak N_*\)` subtract the universal tracking feed-through exactly.

---

## 3. Exact first-order observable packet and its quotient-tangent image

Define the first-order branch-observable packet
\[
\boxed{
\Delta_{\rm obs}^{(1)}:=
\begin{pmatrix}
\delta\ln R_{\rm tr}\\[2pt]
\delta\ln \mathfrak N_*\\[2pt]
\delta\ln \epsilon_\eta
\end{pmatrix}.
}
\]

Stage 187 already fixed the tangent quotient packet as
\[
\Delta_{\rm quot}^{(1)}:=
\begin{pmatrix}
\delta\ln \mathfrak C_{{\rm tr},*}\\[2pt]
\delta\ln \mathfrak C_{{\rm nt},*}\\[2pt]
\delta\ln \epsilon_\eta
\end{pmatrix}
=
\begin{pmatrix}
\Sigma_{\rm tr}\\[2pt]
\Sigma_{\rm nt}\\[2pt]
\Sigma_\eta
\end{pmatrix}.
\]

Using the Stage 184/185 identities,
\[
\delta\ln \mathfrak C_{{\rm tr},*}=-\frac{1}{C_{{\rm tr},*}}\,\delta\ln R_{\rm tr},
\qquad
\delta\ln \mathfrak C_{{\rm nt},*}=\delta\ln \mathfrak N_*,
\qquad
\delta\ln\epsilon_\eta=\delta\ln\epsilon_\eta,
\]
so the exact first-order observable-to-quotient compiler is
\[
\boxed{
\Delta_{\rm quot}^{(1)}
=
\underbrace{
\begin{pmatrix}
-\dfrac{1}{C_{{\rm tr},*}} & 0 & 0\\[8pt]
0 & 1 & 0\\[2pt]
0 & 0 & 1
\end{pmatrix}}_{\mathcal C_{\rm obs\to quot}}
\Delta_{\rm obs}^{(1)}.
}
\]
Its determinant is
\[
\det \mathcal C_{\rm obs\to quot}=-\frac{1}{C_{{\rm tr},*}}\neq 0.
\]
So the observable packet and the tangent quotient packet have exactly the same zero set.

### Interpretation

Stage 187 still remains the exact **finite** quotient closure.
Stage 188 does not replace it by a new finite quotient statement.
It shows something narrower and more useful for the actual PDE solve:

> at first grouped weak-axisymmetric/reference-branch order, the direct branch-observable packet is an exact invertible compiler for the tangent packet of the Stage 187 quotient.

That is the precise sense in which `\(R_{\rm tr}\)`, `\(\mathfrak N_*\)`, and `\(\epsilon_\eta\)` become first-class PDE outputs.

---

## 4. Exact first-order compiler from branch observables to the defect triple

The observable defect packet is
\[
\boxed{
\Delta_{\rm def}^{(1)}:=
\begin{pmatrix}
\Theta_1\\[2pt]
\Xi_1\\[2pt]
\mathcal R_1
\end{pmatrix}.
}
\]
Using
\[
\Theta_1=\delta\ln R_{\rm tr},
\]
\[
\delta\ln\mathfrak N_*=
\delta\ln\mathcal T^2+B_*\delta\ln R_{\rm tr}
=\Xi_1+B_*\Theta_1,
\]
and
\[
\mathcal R_1+\Xi_1
=-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta,
\]
we obtain the exact first-order observable compiler
\[
\boxed{
\Delta_{\rm def}^{(1)}
=
\underbrace{
\begin{pmatrix}
1 & 0 & 0\\[2pt]
-B_* & 1 & 0\\[2pt]
B_* & -1 & -\dfrac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\end{pmatrix}}_{\mathcal C_{\rm obs\to def}}
\Delta_{\rm obs}^{(1)}.
}
\]
Equivalently, componentwise,
\[
\boxed{
\Theta_1=\delta\ln R_{\rm tr},
}
\]
\[
\boxed{
\Xi_1=\delta\ln\mathfrak N_*-B_*\,\delta\ln R_{\rm tr},
}
\]
\[
\boxed{
\mathcal R_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta
-\Xi_1.
}
\]

The inverse compiler is also exact:
\[
\boxed{
\Delta_{\rm obs}^{(1)}
=
\underbrace{
\begin{pmatrix}
1 & 0 & 0\\[2pt]
B_* & 1 & 0\\[2pt]
0 & -\dfrac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}} & -\dfrac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}
\end{pmatrix}}_{\mathcal C_{\rm def\to obs}}
\Delta_{\rm def}^{(1)}.
}
\]
So
\[
\det \mathcal C_{\rm obs\to def}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\neq 0
\qquad (0<\epsilon_{\eta,*}<1).
\]
Therefore the defect packet and the observable packet also have exactly the same zero set.

---

## 5. Factorization through the Stage 187 tangent packet

The Stage 183/185 triangular normal form already gives
\[
\Theta_1=-C_{{\rm tr},*}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{{\rm tr},*}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1+
\Xi_1=-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\Sigma_\eta.
\]
With
\[
\Delta_{\rm quot}^{(1)}=(\Sigma_{\rm tr},\Sigma_{\rm nt},\Sigma_\eta)^T,
\]
this is the exact tangent-packet compiler
\[
\Delta_{\rm def}^{(1)}
=
\underbrace{
\begin{pmatrix}
-C_{{\rm tr},*} & 0 & 0\\[2pt]
A_{{\rm tr},*} & 1 & 0\\[2pt]
-A_{{\rm tr},*} & -1 & -\dfrac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\end{pmatrix}}_{\mathcal C_{\rm quot\to def}}
\Delta_{\rm quot}^{(1)}.
\]

Stage 188 then shows that this factorizes exactly as
\[
\boxed{
\mathcal C_{\rm obs\to def}=
\mathcal C_{\rm quot\to def}\,\mathcal C_{\rm obs\to quot}.
}
\]
The first-column entries match because of the identity
\[
A_{{\rm tr},*}=B_* C_{{\rm tr},*}.
\]

So the Stage 188 observable compiler is not a new independent rule. It is exactly the Stage 187 tangent quotient compiler rewritten in the natural direct branch observables.

---

## 6. Exact first-order branch-observable zero-defect theorem

Because both compilers above are invertible, we obtain the sharp observable theorem:
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=0,
\quad
\delta\ln \mathfrak N_*=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]
Equivalently,
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=0,
\quad
\delta\ln(\mathcal T^2R_{\rm tr}^{B_*})=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]
And because the tangent quotient packet is an exact compiler image of the same observable packet,
\[
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
\]
So the Stage 184 branch-composite theorem and the Stage 187 tangent quotient theorem are now joined cleanly.

### Interpretation

The actual moving-throat branch no longer has to be tested against a broad microscopic slippage ledger to answer the first-order defect question.
It only has to answer three direct branch-observable questions:

1. is the coherent tracking factor `\(R_{\rm tr}\)` invariant?
2. is the corrected nontracking composite `\(\mathfrak N_*\)` invariant?
3. is the microscopic dressing ratio `\(\epsilon_\eta\)` invariant?

If yes, then the full coherent first grouped weak-axisymmetric defect vanishes.

---

## 7. Best current theorem statement after Stage 188

Stage 187 already reduced the coherent weak-axisymmetric problem to an exact finite quotient.
Stage 188 now adds the direct PDE-facing front end.

At the first grouped weak-axisymmetric/reference-branch order:

1. the exact finite quotient tangent packet
   \[
   (\delta\ln\mathfrak C_{{\rm tr},*},\;\delta\ln\mathfrak C_{{\rm nt},*},\;\delta\ln\epsilon_\eta)
   \]
   is an exact invertible compiler image of the branch-observable packet
   \[
   (\delta\ln R_{\rm tr},\;\delta\ln\mathfrak N_*,\;\delta\ln\epsilon_\eta),
   \]
2. the defect packet
   \[
   (\Theta_1,\Xi_1,\mathcal R_1)
   \]
   is an exact invertible compiler image of that same branch-observable packet,
3. and therefore the first grouped weak-axisymmetric zero-defect theorem can be written directly in the exact branch observables
   \[
   R_{\rm tr},\qquad \mathfrak N_*=\mathcal T^2R_{\rm tr}^{B_*},\qquad \epsilon_\eta.
   \]

So the next honest theorem gate for the completed PDE is now completely sharp:

> compute the first grouped weak-axisymmetric drifts of
> \[
> R_{\rm tr},
> \qquad
> \mathfrak N_*=
> \mathcal T^2R_{\rm tr}^{B_*},
> \qquad
> \epsilon_\eta,
> \]
> on the actual moving-throat branch.

That is the smallest branch-specific front end of the coherent defect theorem reached so far.

---

## 8. Immediate next stage

The right next stage after this one is the first one that uses these exact branch observables as the PDE-facing inputs of the later grouped-bundle/outgoing compiler.

In particular, the next clean move is to keep the present Stage 188 observable packet fixed and then ask how the completed grouped real `P2` bundle packages those observables into the finite branch packet that controls:

- grouped isotropy,
- one-pole conservative closure,
- and outgoing quadrupole normalization.

So Stage 188 is the natural hinge between the exact orbit–quotient closure and the later endgame packet compiler.
