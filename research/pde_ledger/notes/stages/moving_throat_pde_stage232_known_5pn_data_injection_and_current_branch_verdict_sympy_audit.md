
# Moving-Throat PDE — Stage 249: Known 5PN Data Injection and Current Branch Verdict

## Status

**Exact within the carried same-charge barrier chain plus the numerically located `5`PN Family-1 support/source branch on the refreshed `\Lambda_{\rm EM}` geometry.**

This stage does **not** solve the full moving-throat PDE.
It takes the support/source data that are already numerically present in the `5`PN stack, injects them into the same-charge audit language, and asks whether the corridor dies there or whether the unresolved bottleneck stays where the earlier audit said it was.

---

## Purpose

The earlier same-charge barrier audit had already narrowed the live corridor to a small ordered chain:

1. the **dynamic** wall-like window is *not* the first kill condition;
2. the **support/source** side must still be large enough;
3. the first unresolved kill test is the **static** placement / orbit-lock side.

The `5`PN continuation then did something new: it actually **numerically located** the explicit Family-1 support/source branch on the refreshed exact `\Lambda_{\rm EM}` geometry. That means the next honest stage is no longer more support algebra. It is simply:

> inject the known `5`PN support/source data into the same-charge audit chain and decide whether the live bottleneck moves, or stays on the unresolved static orbit-lock / coherent placement packet.

The main result is clear:

> after the known `5`PN data are injected, the same-charge corridor is still alive.
> The support/source side is numerically safe by a large margin, while the first unresolved gate remains the actual PDE-selected static orbit-lock / placement point.

---

## 1. What is already numerically known from the `5`PN stack

### 1.1 Refreshed exact `\Lambda_{\rm EM}` geometry branch

With the carried thin-wall Family-1 lock
\[
\frac{\ell}{a}=\frac1{20},
\]
the exact EM-branch geometry refresh gives
\[
\Lambda_\ell:=\frac{L}{\ell}
=
20\,\Lambda_{\rm EM}
=
\frac{20\sqrt2\,\pi}{x_{01}}
\approx 36.94973154240256,
\]
where \(x_{01}\) is the first positive zero of \(J_0\).

The same branch fixes
\[
\eta=\Lambda_\ell\approx 36.94973154240256,
\]
\[
\chi_s=\frac{L}{2\ell}=\frac{\Lambda_\ell}{2}\approx 18.47486577120128,
\]
\[
\kappa
=
4\chi_s^2+\frac45\Lambda_\ell^2
=
\frac95\Lambda_\ell^2
\approx 2457.5087899001137.
\]

On the Robin support side, the lowest support eigenvalue solves
\[
y\tan y=\eta,
\qquad
0<y<\frac{\pi}{2},
\]
so the exact support-softening factor is
\[
A_K
=
\frac{\kappa+\pi^2/4}{\kappa+y^2}
\approx 1.0000521380385143,
\]
and the Family-1 support ceiling is
\[
\boxed{
\zeta_{\max}
=
A_K\frac{\pi^2}{4}
\approx 2.4675297457259358.
}
\]

### 1.2 Exact support-drop kernel and fixed-point equation

The same support/source branch is selected by the exact fixed-point law
\[
Pe=\Xi\,\Delta(Pe;\kappa,\eta),
\]
with stationary zero-flux source profile
\[
\Sigma_{Pe}(x)=\frac{Pe\,e^{Pe x}}{e^{Pe}-1},
\qquad x\in[0,1],
\]
and support-drop kernel
\[
K_{\kappa,\eta}(x)
=
\frac{\cosh(\alpha x)+(\eta/\alpha)\sinh(\alpha x)-\cosh(\alpha(1-x))}
{\alpha\sinh\alpha+\eta\cosh\alpha},
\qquad
\alpha=\sqrt\kappa.
\]
So
\[
\Delta(Pe;\kappa,\eta)
=
\int_0^1 K_{\kappa,\eta}(x)\,\Sigma_{Pe}(x)\,dx.
\]

Its exact endpoint values are
\[
\Delta_0(\kappa,\eta)
=
\frac{\eta(\cosh\alpha-1)}
{\alpha^2(\alpha\sinh\alpha+\eta\cosh\alpha)}
\approx 1.7377393923469950\times 10^{-4},
\]
\[
\Delta_\infty(\kappa,\eta)
=
\frac{\cosh\alpha+(\eta/\alpha)\sinh\alpha-1}
{\alpha\sinh\alpha+\eta\cosh\alpha}
\approx 2.0172162594593645\times 10^{-2}.
\]

So every constructive root obeys the exact branch bracket
\[
\Xi\,\Delta_0
\le Pe_*
\le
\Xi\,\Delta_\infty.
\]

### 1.3 Two explicit wall-depth extractions

On the carried Family-1 branch, the two explicit wall-depth extractions are
\[
\Theta_w^{(\chi)}\approx 4.06863235008162,
\qquad
\Theta_w^{(J)}\approx 0.927552032539308
\]
for the benchmark \(\lambda_\mu=1\).

The exact wall/source figures of merit are then
\[
\Xi_\chi = 168\,\Theta_w^{(\chi)}\Lambda_\ell^2
\approx 5.5548332017764099\times 10^5,
\]
\[
\Xi_J = 168\,\Theta_w^{(J)}\Lambda_\ell^2
\approx 1.2663707072528143\times 10^5.
\]

The numerically located fixed-point roots are:

#### \(\chi\)-weighted extraction
\[
Pe_*^{(\chi)}\approx 11155.7265863205869,
\]
\[
\zeta_{\rm phys}^{(\chi)}
=
A_K\,\Omega_{Pe_*^{(\chi)}}^2
\approx 2.4675296478814376,
\]
\[
\rho_{\alpha,\max}^{(\chi)}
=
1+\zeta_{\rm phys}^{(\chi)}
\approx 3.4675296478814376.
\]

#### \(J\)-weighted extraction
\[
Pe_*^{(J)}\approx 2504.9703142859238,
\]
\[
\zeta_{\rm phys}^{(J)}
=
A_K\,\Omega_{Pe_*^{(J)}}^2
\approx 2.4675278051675084,
\]
\[
\rho_{\alpha,\max}^{(J)}
=
1+\zeta_{\rm phys}^{(J)}
\approx 3.4675278051675084.
\]

Here the exact overlap-boost law is
\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe\,e^{Pe}+\pi\bigr)}
{(4Pe^2+\pi^2)(e^{Pe}-1)},
\]
or, in numerically stable form,
\[
\Omega_{Pe}
=
\frac{\pi Pe\bigl(2Pe+\pi e^{-Pe}\bigr)}
{(4Pe^2+\pi^2)(1-e^{-Pe})}.
\]

So the support/source side is no longer just “probably okay.”
It is numerically located on the actual explicit Family-1 operator branch.

---

## 2. Exact margins after plugging in the known `5`PN data

The natural isotropic passive/outgoing grouped-\(P_2\) branch still requires only
\[
\zeta_{\rm req}=\frac13,
\qquad
\rho_\alpha^{\rm req}=\frac43.
\]

So the injected support/source safety margins are:

### \(\chi\)-weighted branch
\[
\zeta_{\rm phys}^{(\chi)}-\zeta_{\rm req}
\approx 2.1341963145481043,
\]
\[
\rho_{\alpha,\max}^{(\chi)}-\frac43
\approx 2.1341963145481043.
\]

### \(J\)-weighted branch
\[
\zeta_{\rm phys}^{(J)}-\zeta_{\rm req}
\approx 2.1341944718341751,
\]
\[
\rho_{\alpha,\max}^{(J)}-\frac43
\approx 2.1341944718341751.
\]

Useful ratio form:

### \(\chi\)-weighted branch
\[
\frac{\zeta_{\rm phys}^{(\chi)}}{\zeta_{\rm req}}
\approx 7.402588943644313,
\qquad
\frac{\rho_{\alpha,\max}^{(\chi)}}{4/3}
\approx 2.600647235911078.
\]

### \(J\)-weighted branch
\[
\frac{\zeta_{\rm phys}^{(J)}}{\zeta_{\rm req}}
\approx 7.402583415502525,
\qquad
\frac{\rho_{\alpha,\max}^{(J)}}{4/3}
\approx 2.600645853875631.
\]

So the explicit support/source branch overshoots the canonical isotropic demand by a factor of about \(7.4\) in the \(\zeta\) variable.

At the same time, both operator-selected points sit extremely close to the Family-1 ceiling:

### \(\chi\)-weighted branch
\[
\zeta_{\max}-\zeta_{\rm phys}^{(\chi)}
\approx 9.784449817674381\times 10^{-8}.
\]

### \(J\)-weighted branch
\[
\zeta_{\max}-\zeta_{\rm phys}^{(J)}
\approx 1.940558427350484\times 10^{-6}.
\]

That is not a contradiction.
It simply means that the numerically selected Family-1 support/source point nearly saturates the **Family-1 support ceiling** while still sitting far above the much smaller isotropic demand \(\zeta_{\rm req}=1/3\).

---

## 3. Transported consequence for the same-charge audit chain

The earlier same-charge audit had already established two load-bearing orderings before any `5`PN injection:

1. the **dynamic** selected-branch wall-window is not the first kill condition;
2. the first real bottleneck is the **static** transported placement / \(\Xi_1\) side.

Injecting the numerically located `5`PN support/source data does **not** change that ordering.
It sharpens it.

The carried consequences are now:

- the support/source side is numerically safe;
- the canonical passive/outgoing normalization side stays exact on the natural isotropic branch;
- so the first unresolved kill condition remains exactly where the reduced `5`PN finish-line notes already said it was:
  the actual PDE-selected orbit-lock / coherent placement point.

In the exact `5`PN finish-line language, the still-missing numerical object is the actual branch point satisfying
\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln \epsilon_\eta=0,
\]
together with the canonical outgoing-normalization condition
\[
N_Q=1.
\]

So after the known `5`PN data are injected, the support/source side is no longer the place where the same-charge idea should die inside the current reduced theorem stack.

---

## 4. Current best verdict after plugging in the known numbers

The same-charge corridor is still alive.

More precisely:

1. **Static same-sign Maxwell shaping** is still not the answer.
2. **Dynamic resonance** is still not the first bottleneck.
3. **Support/source enhancement** is now numerically located and strongly non-bottlenecked.
4. The first unresolved gate is still the **actual PDE-selected static orbit-lock / placement point**.

So after Stage 249 the question is no longer:

> can the reduced support/source side possibly be large enough?

It already is.

The question is now:

> when the completed moving-throat branch is actually selected, does its realized orbit packet / static placement verdict land inside the surviving same-charge window, or does the route finally die there?

That is the cleanest current stopping point.

---

## 5. SymPy-backed status

The accompanying audit script verifies all of the following:

1. the exact refreshed `\Lambda_{\rm EM}` geometry formulas
   \[
   \Lambda_\ell=\frac{20\sqrt2\,\pi}{x_{01}},
   \qquad
   \kappa=\frac95\Lambda_\ell^2;
   \]
2. the exact Robin support equation
   \[
   y\tan y=\eta
   \]
   and the derived Family-1 support ceiling
   \[
   \zeta_{\max}=A_K\pi^2/4;
   \]
3. the exact support-drop endpoint formulas
   \[
   \Delta_0(\kappa,\eta),
   \qquad
   \Delta_\infty(\kappa,\eta);
   \]
4. the exact fixed-point equation
   \[
   Pe=\Xi\,\Delta(Pe;\kappa,\eta)
   \]
   for both wall-depth extractions;
5. the numerically located roots
   \[
   Pe_*^{(\chi)},\qquad Pe_*^{(J)};
   \]
6. the exact transported support/source values
   \[
   \zeta_{\rm phys}^{(\chi)},\quad
   \zeta_{\rm phys}^{(J)},\quad
   \rho_{\alpha,\max}^{(\chi)},\quad
   \rho_{\alpha,\max}^{(J)};
   \]
7. the exact margins above
   \[
   \zeta_{\rm req}=\frac13,
   \qquad
   \rho_\alpha^{\rm req}=\frac43;
   \]
8. and the proximity of both explicit Family-1 points to the ceiling \(\zeta_{\max}\).

So the note is not just a verbal carry-forward statement.
It is backed by an executable reconstruction of the known `5`PN data injection step.

---

## 6. What the next honest stage should do

The next stage should **not** invent more support algebra.

The honest next theorem gate is now:

1. take the exact unresolved coherent placement packet
   \[
   (R_{\rm tr},R_{\rm target},\epsilon_\eta),
   \]
   or equivalently
   \[
   (d\ln R_{\rm tr},d\ln R_{\rm target},d\ln\epsilon_\eta),
   \]
2. express its weak-axisymmetric verdict as the actual static \(\Xi_1\) / transported placement scalar used in the same-charge chain,
3. and test whether the realized branch clears the already-carried static ceiling.

That is where the present stack says the real answer now lives.
