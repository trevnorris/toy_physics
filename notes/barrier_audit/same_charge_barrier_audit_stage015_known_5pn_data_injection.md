# Same-Charge Barrier Audit — Stage 015: Known 5PN Data Injection and Current Branch Verdict

## Purpose

The earlier same-charge barrier audit reduced the live corridor to a very narrow condition set:

1. the **dynamic** wall-like window is *not* the first kill condition;
2. the **support/source** side must supply enough coherent enhancement;
3. the first unresolved kill test is the **static** placement / orbit-lock side.

This stage plugs the **actually numerically located** 5PN support/source data into that audit chain.

The point is not to solve the full moving-throat branch. The current 5PN notes are explicit that the support/source side has already been numerically located, while the actual PDE-selected orbit-lock / coherent placement point is still not numerically present in the files. So the clean question is:

> once the numerically located 5PN support/source data are inserted, does the same-charge corridor die, or does the unresolved bottleneck stay where the earlier audit said it was?

---

## 1. What is already numerically known from the 5PN stack

The numerically located Family-1 support/source branch gives the following exact values on the refreshed geometry branch:

\[
\Lambda_\ell \approx 36.94973154240256,
\qquad
\kappa \approx 2457.5087899001137,
\qquad
\zeta_{\max} \approx 2.4675297457259358.
\]

On the two explicit extraction branches:

### \(\chi\)-weighted extraction
\[
Pe_*^{(\chi)} \approx 11155.7265863205869,
\]
\[
\zeta_{\rm phys}^{(\chi)} \approx 2.4675296478814376,
\]
\[
\rho_{\alpha,\max}^{(\chi)} \approx 3.4675296478814376.
\]

### \(J\)-weighted extraction
\[
Pe_*^{(J)} \approx 2504.9703142859238,
\]
\[
\zeta_{\rm phys}^{(J)} \approx 2.4675278051675084,
\]
\[
\rho_{\alpha,\max}^{(J)} \approx 3.4675278051675084.
\]

The natural isotropic passive/outgoing grouped-\(P_2\) branch still requires only

\[
\zeta_{\rm req}=\frac13,
\qquad
\rho_\alpha^{\rm req}=\frac43.
\]

So the support/source branch is no longer just “probably okay.” It is numerically located and strongly above the canonical isotropic demand.

---

## 2. Exact margins after plugging in the known 5PN data

The support/source safety margins are:

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
\approx 2.600647235911438.
\]

### \(J\)-weighted branch
\[
\frac{\zeta_{\rm phys}^{(J)}}{\zeta_{\rm req}}
\approx 7.402583415502525,
\qquad
\frac{\rho_{\alpha,\max}^{(J)}}{4/3}
\approx 2.600645853875631.
\]

So the explicit support/source branch overshoots the canonical isotropic requirement by a factor of about \(7.4\) in the \(\zeta\) variable.

At the same time, the numerically selected Family-1 points sit very close to the Family-1 ceiling:

\[
\zeta_{\max}-\zeta_{\rm phys}^{(\chi)}
\approx 9.784449820573488\times 10^{-8},
\]
\[
\zeta_{\max}-\zeta_{\rm phys}^{(J)}
\approx 1.9405584274474645\times 10^{-6}.
\]

That is not a contradiction. It simply means the numerically selected Family-1 branch nearly saturates the *Family-1 support* ceiling while still sitting far above the much smaller isotropic demand \(\zeta_{\rm req}=1/3\).

---

## 3. Transported consequence for the same-charge audit chain

The earlier same-charge audit already proved two things before any 5PN injection:

1. the **dynamic** selected-branch wall-window is not the first kill condition;
2. the first real bottleneck is the **static** transported placement / \(\Xi_1\) side.

Injecting the numerically located 5PN support/source data does not change that ordering.
Instead it sharpens it:

- the support/source side is numerically safe;
- the canonical passive/outgoing normalization side remains exact on the natural branch;
- so the first unresolved kill condition stays where the 5PN notes say it is:
  the actual PDE-selected orbit-lock / coherent placement point.

In the 5PN finish-line language, the missing numerical object is still the actual branch point satisfying

\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln \epsilon_\eta=0,
\]
together with the canonical outgoing normalization condition
\[
N_Q=1.
\]

So the support/source side is no longer the place where the same-charge idea should die inside the current reduced theorem stack.

---

## 4. Current best verdict after plugging in the known numbers

The same-charge corridor is still alive.

More precisely:

1. **Static same-sign Maxwell shaping** is still not the answer.
2. **Dynamic resonance** is still not the first bottleneck.
3. **Support/source enhancement** is now numerically located and strongly non-bottlenecked.
4. The first unresolved gate is still the **actual PDE-selected static orbit-lock / placement point**.

So after Stage 015 the question is no longer:

> can the reduced support/source side possibly be large enough?

It already is.

The question is now:

> when the completed moving-throat branch is actually selected, does its realized orbit packet / static placement verdict land inside the surviving same-charge window, or does the route finally die there?

That is the cleanest current stopping point.

---

## 5. What the next honest stage should do

The next stage should not invent more support algebra. It should compile the actual unresolved packet into the same-charge audit language.

Concretely, the next theorem gate is:

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
