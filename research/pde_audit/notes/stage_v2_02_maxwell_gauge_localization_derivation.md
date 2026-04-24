# Stage V2-02 — Maxwell Gauge-Localization Audit

## Purpose

This stage audits the localized Maxwell sector used by the 4D / moving-throat program, focusing on the fact that the Maxwell kinetic term is localized by a transverse profile \(Z(w)\), while the displayed Lorenz-type gauge-fixing term is not localized.

The current action form is

\[
S_{\rm EM}
=
\int dt\,d^3x\,dw\,
\left[
-\frac{Z(w)}{4\mu_0}F_{MN}F^{MN}
-\frac{1}{2\xi\mu_0}(\partial\!\cdot\!A)^2
-A_MJ^M
\right].
\]

The audit asks whether this is safe for the bulk equations, the zero-mode brane reduction, and the moving-throat mixed-sector response stack.

---

## 1. General weighted gauge-fixing functional

To separate the issue cleanly, introduce a general gauge-fixing weight \(H(w)\):

\[
\mathcal L_{\rm gf}
=
-\frac{H(w)}{2\xi\mu_0}\,B^2,
\qquad
B\equiv\partial_MA^M.
\]

The current papers correspond to

\[
H(w)=1.
\]

The localized patch corresponds to something like

\[
H(w)=Z(w)
\]

or, more generally, any finite localized profile \(H_{\rm loc}(w)\) with finite integral.

Varying the gauge-fixing term gives

\[
\delta S_{\rm gf}
=
-\frac{1}{\xi\mu_0}
\int H B\,\partial_M\delta A^M
=
\frac{1}{\xi\mu_0}
\int \partial_M(HB)\,\delta A^M,
\]

so the Maxwell equation becomes

\[
\boxed{
\partial_M\bigl(ZF^{MN}\bigr)
+
\frac{1}{\xi}\partial^N\bigl(HB\bigr)
=
\mu_0J^N.
}
\]

For \(H=1\), this reduces to the current displayed equation:

\[
\partial_M\bigl(ZF^{MN}\bigr)
+
\frac{1}{\xi}\partial^N B
=
\mu_0J^N.
\]

For \(H=Z\), it becomes

\[
\partial_M\bigl(ZF^{MN}\bigr)
+
\frac{1}{\xi}\partial^N\bigl(ZB\bigr)
=
\mu_0J^N.
\]

---

## 2. Divergence consistency

Because \(ZF^{MN}\) is antisymmetric in \(M,N\), while \(\partial_M\partial_N\) is symmetric,

\[
\partial_N\partial_M\bigl(ZF^{MN}\bigr)=0.
\]

The script verifies this explicitly in a symbolic antisymmetric subblock.

Taking the divergence of the general equation gives

\[
\boxed{
\frac{1}{\xi}\Box_5\bigl(HB\bigr)
=
\mu_0\partial_NJ^N.
}
\]

So for conserved sources,

\[
\Box_5(HB)=0.
\]

For \(H=1\), this is

\[
\Box_5B=0.
\]

For \(H=Z\), it is

\[
\Box_5(ZB)=0.
\]

The two are not identical as gauge-propagation equations. Since \(H=Z\) depends on \(w\),

\[
\Box_5(HB)
=
H\Box_5B+2H'\partial_wB+H''B.
\]

The script verifies the identity

\[
\Box(HB)-\left(H\Box B+2H'B_w+H''B\right)=0.
\]

This is not a fatal issue. Gauge fixing is not a physical observable. But it means that switching from \(H=1\) to \(H=Z\) is not a purely cosmetic rewrite; it changes the gauge-condition propagation while preserving gauge-invariant observables when handled consistently.

---

## 3. Zero-mode brane reduction

Take the usual zero-mode ansatz

\[
A_\mu(x,w)\simeq a_\mu(x),
\qquad
A_w=0,
\qquad
\partial_wA_\mu\simeq0.
\]

Then

\[
F_{\mu\nu}(x,w)=f_{\mu\nu}(x),
\qquad
F_{\mu w}=0,
\qquad
B=\partial_\mu a^\mu.
\]

The Maxwell kinetic term reduces to

\[
S_{\rm kin}^{(0)}
=
-\frac{Z_{\rm int}}{4\mu_0}
\int d^4x\,f_{\mu\nu}f^{\mu\nu},
\qquad
Z_{\rm int}=\int_{-\infty}^{\infty}Z(w)\,dw.
\]

For the Gaussian profile

\[
Z(w)=e^{-w^2/\lambda^2},
\]

SymPy verifies

\[
\boxed{Z_{\rm int}=\sqrt\pi\,\lambda.}
\]

The effective coupling is therefore

\[
\mu_0^{\rm eff}=\frac{\mu_0}{Z_{\rm int}}.
\]

### 3.1 Unweighted gauge fixing

For the current unweighted term \(H=1\), the zero-mode gauge-fixing contribution on a finite regulator interval \([-R,R]\) is

\[
S_{\rm gf,unweighted}^{(0)}(R)
=
-\frac{2R}{2\xi\mu_0}
\int d^4x\,\left(\partial_\mu a^\mu\right)^2.
\]

Thus

\[
\frac{H_{\rm int}}{Z_{\rm int}}
=
\frac{2R}{\sqrt\pi\lambda}
\longrightarrow \infty
\qquad(R\to\infty).
\]

So the unweighted gauge-fixing action diverges for a noncompact zero mode unless

\[
\partial_\mu a^\mu=0
\]

is imposed before reduction.

Equivalently, if one tries to interpret the finite-box unweighted term as a 4D gauge-fixing term with

\[
-\frac{Z_{\rm int}}{2\xi_4\mu_0}(\partial\!\cdot a)^2,
\]

then matching gives

\[
\xi_4=\xi\frac{Z_{\rm int}}{2R}
\longrightarrow 0.
\]

So the noncompact unweighted zero-mode limit is a singular Landau-type gauge limit, not an ordinary finite 4D gauge-fixed action.

### 3.2 Weighted localized gauge fixing

For \(H=Z\), the zero-mode gauge-fixing term is instead

\[
S_{\rm gf,weighted}^{(0)}
=
-\frac{Z_{\rm int}}{2\xi\mu_0}
\int d^4x\,\left(\partial_\mu a^\mu\right)^2.
\]

This is finite and has the same \(Z_{\rm int}\) normalization as the Maxwell kinetic term.

For any finite gauge weight \(H_{\rm int}=\int H(w)dw\), the 4D gauge parameter obeys

\[
\boxed{
\xi_4=\xi\frac{Z_{\rm int}}{H_{\rm int}}.
}
\]

The script verifies the matching identity

\[
\frac{H_{\rm int}}{\xi}=\frac{Z_{\rm int}}{\xi_4}.
\]

For \(H=Z\), one gets \(\xi_4=\xi\).

---

## 4. Mixed-sector gauge invariants

The moving-throat program must keep the mixed sector alive:

\[
A_w,
\qquad
F_{\mu w},
\qquad
J^w,
\qquad
E_w,
\qquad
C_a.
\]

Using the convention

\[
A_0\to A_0-\partial_t\chi,
\qquad
A_a\to A_a+\partial_a\chi,
\qquad
A_w\to A_w+\partial_w\chi,
\]

the mixed fields are

\[
E_w=-\partial_tA_w-\partial_wA_0,
\]

\[
C_a=\partial_aA_w-\partial_wA_a.
\]

The script verifies

\[
\delta E_w=0,
\qquad
\delta C_a=0.
\]

So the mixed fields are exact gauge-invariant observables. They should not be treated as artifacts of the gauge-fixing choice.

---

## 5. Stage verdict

### 5.1 What passes

The bulk equation with the current unweighted \(H=1\) Lorenz term is mathematically consistent as a formal gauge-fixed bulk equation:

\[
\partial_M(ZF^{MN})+rac1\xi\partial^N(\partial\!\cdot A)=\mu_0J^N.
\]

Its divergence identity is also consistent when the total current is conserved:

\[
\Box_5(\partial\!\cdot A)=0.
\]

The mixed-sector observables \(E_w\) and \(C_a\) are gauge invariant and survive independently of this choice.

### 5.2 What fails

The unweighted gauge-fixing term fails as a noncompact zero-mode gauge-fixed action. For \(A_\mu(x,w)=a_\mu(x)\), the Maxwell kinetic coefficient is finite:

\[
Z_{\rm int}=\sqrt\pi\lambda,
\]

but the unweighted gauge-fixing coefficient is

\[
\int_{-\infty}^{\infty}dw=\infty.
\]

Therefore the current unweighted term cannot be naively reduced to a standard finite 3+1 gauge-fixed Maxwell action by replacing the integral with \(Z_{\rm int}\).

### 5.3 Safe interpretations

There are two safe ways to state the program.

#### Option A — keep the current unweighted gauge fixing as a bulk gauge device

Use \(H=1\) only at the bulk gauge-fixed equation level, impose Lorenz gauge

\[
\partial\!\cdot A=0
\]

before the zero-mode brane reduction, and then choose any convenient 3+1 gauge fixing after the reduced Maxwell theory has been obtained.

In this reading, the clean zero-mode Maxwell equation is obtained from the gauge-invariant part:

\[
\partial_\mu f^{\mu\nu}=\frac{\mu_0}{Z_{\rm int}}J_b^\nu.
\]

The unweighted gauge-fixing term is not part of the reduced finite action.

#### Option B — localize the gauge-fixing functional

Replace the gauge-fixing term by

\[
\mathcal L_{\rm gf}^{\rm loc}
=
-\frac{H_{\rm loc}(w)}{2\xi\mu_0}(\partial\!\cdot A)^2,
\]

with finite

\[
H_{\rm int}=\int H_{\rm loc}(w)dw.
\]

The simplest choice is

\[
H_{\rm loc}=Z.
\]

Then the zero-mode gauge-fixed action is finite and has the same localization normalization as the kinetic term. The cost is that the gauge-condition propagation becomes

\[
\Box_5(Z\partial\!\cdot A)=0
\]

for conserved sources.

---

## 6. Recommended Volume 2 patch

The cleanest patch for the derivation ledger is to state the Maxwell sector in one of these two forms:

### Conservative patch

Keep the already-published equation, but add the warning:

> The unweighted Lorenz term is a bulk gauge-fixing device. In the zero-mode brane reduction, the Lorenz condition is imposed before integrating over the noncompact transverse coordinate. A finite 3+1 gauge-fixing term may then be chosen in the reduced brane Maxwell theory. The unweighted five-dimensional gauge-fixing term should not be integrated over the noncompact zero mode unless \(\partial_\mu a^\mu=0\).

### Structural patch

Use the localized gauge-fixing action

\[
\boxed{
\mathcal L_{\rm gf}^{\rm loc}
=
-\frac{Z(w)}{2\xi\mu_0}(\partial\!\cdot A)^2.
}
\]

Then the brane zero mode reduces cleanly to a finite 3+1 gauge-fixed Maxwell action with \(\xi_4=\xi\).

My recommendation is to use the conservative patch immediately for consistency with the existing papers, then mark the structural patch as the preferred parent-action cleanup for any future fully gauge-fixed PDE derivation.

---

## 7. Final status line

\[
\boxed{
\text{The unweighted Maxwell gauge-fixing term is not fatal, but it is unsafe for naive zero-mode action reduction.}
}
\]

A correct Volume 2 statement is:

\[
\boxed{
\text{Either impose Lorenz gauge before reducing, or localize the gauge-fixing functional.}
}
\]
