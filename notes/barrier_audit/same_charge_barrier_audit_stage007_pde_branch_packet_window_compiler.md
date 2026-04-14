# Same-Charge Barrier Audit — Stage 007: PDE Branch-Packet Compiler, Weak-Axisymmetric Ceiling Transport, and the First Actual-Branch Kill Test

## 0. Purpose

Stage 006 forced the primitive finite-throat one-port branch onto the exact isotropic `5`PN target surface and found that the dynamic same-charge corridor survives only inside a **finite normalization window**.

But that still left one gap:

> how do we test the **actual** moving-throat / `5`PN branch, rather than the primitive family, against that survival window?

The next honest step is therefore not another primitive scan. It is to compile the real PDE-selected branch data into the same variables the Stage-006 window cares about.

That is what this stage does.

The main outputs are:

1. the exact compiler from the final `5`PN branch packet to the lane prefactors
   \[
   P_0^{(20)},\qquad P_0^{(21)},\qquad P_0^{(22)},
   \]
2. the exact isotropic window test in terms of the normalization defect `\Delta_{\rm norm}`,
3. the exact weak-axisymmetric transported ceiling test in terms of the grouped prefactor defects `a_{P_0}, b_{P_0}`,
4. the axisymmetric specialization
   \[
   b_{P_0}=3a_{P_0}
   \]
   and the equivalent one-scalar amplitude law
   \[
   \Xi_1=\frac{P_1}{P_0},
   \]
5. and the explicit anisotropy headroom left at the Stage-006 compatibility point.

So after Stage 007, the problem is no longer

> “extract the normalization target from the PDE somehow.”

It is now

> “does the actual branch packet land inside one explicit finite corridor in `\Delta_{\rm norm}` and weak-axisymmetric prefactor slope?”

---

## 1. Frozen input carried forward

### 1.1 Exact final branch packet from the `5`PN endgame

The exact reduced branch verdict packet is
\[
\Delta_{\rm branch}
=
(a_2,b_2,a_4,b_4,a_{P_0},b_{P_0},\Delta_{\rm pole},\Delta_{\rm norm}).
\]
Its normalization slot is
\[
\Delta_{\rm norm}
=
\hat m_0^{\,2}\,\bar P_0-
\frac{54Gc_s^5}{5a^5c^5}.
\]
So the actual isotropic mean prefactor seen by the same-charge window test is already fixed exactly by
\[
\boxed{
\bar P_0=
\frac{\Delta_{\rm norm}+\dfrac{54Gc_s^5}{5a^5c^5}}{\hat m_0^{\,2}}.
}
\]

### 1.2 Exact grouped inverse map for the lane prefactors

The grouped prefactor anomalies compile back to the three lane prefactors by
\[
\boxed{P_0^{(20)}=\bar P_0+4a_{P_0},}
\]
\[
\boxed{P_0^{(21)}=\bar P_0-a_{P_0}+b_{P_0},}
\]
\[
\boxed{P_0^{(22)}=\bar P_0-a_{P_0}-b_{P_0}.}
\]
So the actual moving-throat branch packet already tells us exactly which lane-wise static prefactors have to be compared against the transported same-charge window.

### 1.3 Finite survival ceilings carried from Stage 006

Stage 006 found four useful primitive-family ceilings.

At the stricter `10%`-loss benchmark:

- **both wall-like poles survive** up to
  \[
  P_{\rm both}^{(10)}
  =0.0028313316855593175,
  \]
- a **nonempty wall-like corridor** survives up to
  \[
  P_{\rm one}^{(10)}
  =0.0035965105896846573.
  \]

At the looser `30%`-loss benchmark:

- **both wall-like poles survive** up to
  \[
  P_{\rm both}^{(30)}
  =0.00817339430971383,
  \]
- a **nonempty wall-like corridor** survives up to
  \[
  P_{\rm one}^{(30)}
  =0.0116633929790174.
  \]

These are still transported primitive-family ceilings, not final full-PDE theorems. But they are the first exact dynamic windows we have, so they are the right actual-branch kill test to carry forward.

---

## 2. Exact compiler from the actual branch packet to the transported window test

The branch packet now enters the same-charge audit in the most direct possible way.

For any chosen ceiling `P_\mathrm{crit}`, define the three actual lane prefactors
\[
P_{20}=\bar P_0+4a_{P_0},
\qquad
P_{21}=\bar P_0-a_{P_0}+b_{P_0},
\qquad
P_{22}=\bar P_0-a_{P_0}-b_{P_0}.
\]
Then the strongest transported sufficient condition that **all grouped lanes** stay inside that primitive-family window is
\[
\boxed{
\max\{P_{20},P_{21},P_{22}\}
\le P_{\rm crit}.
}
\]
Equivalently, in packet variables,
\[
\boxed{
\max\Bigl\{
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}+4a_{P_0},
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}-a_{P_0}+b_{P_0},
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}-a_{P_0}-b_{P_0}
\Bigr\}
\le P_{\rm crit},
}
\]
where
\[
T_{\rm quad}:=\frac{54Gc_s^5}{5a^5c^5}.
\]

So the Stage-006 window has now been converted into a direct inequality on the actual branch packet.

---

## 3. Exact isotropic kill test in terms of `\Delta_{\rm norm}`

If the actual branch is exactly isotropic at the prefactor level,
\[
a_{P_0}=b_{P_0}=0,
\]
then every lane sees the same static prefactor
\[
P_{20}=P_{21}=P_{22}=\bar P_0.
\]
So the transported isotropic ceiling test is simply
\[
\boxed{\bar P_0\le P_{\rm crit}.}
\]
Using the exact normalization compiler,
\[
\boxed{
\Delta_{\rm norm}
\le
\hat m_0^{\,2}P_{\rm crit}-T_{\rm quad}.
}
\]
This is the first actual-branch same-charge kill test written directly in the endgame residual language.

### 3.1 Exact calibrated-branch lower bound on `\hat m_0`

If the real branch already hits the universal quadrupole normalization exactly,
\[
\Delta_{\rm norm}=0,
\]
then isotropic survival at ceiling `P_{\rm crit}` is equivalent to
\[
\boxed{
\hat m_0^{\,2}
\ge
\frac{T_{\rm quad}}{P_{\rm crit}}.
}
\]
So the same-charge corridor now places a direct lower bound on the source-map factor of the actual calibrated branch.

That is a much sharper statement than anything before Stage 006.

---

## 4. Weak-axisymmetric specialization from the exact grouped signature

The later moving-throat grouped notes collapse the weak-axisymmetric outgoing slippage bundle to one scalar amplitude with grouped signature
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
and identify that amplitude with the physical outgoing-prefactor slope
\[
\boxed{\Xi_1=\frac{P_1}{P_0}.}
\]
So the weak-axisymmetric prefactor lanes take the exact first-order form
\[
\boxed{P_A=\bar P_0\bigl(1+\epsilon\lambda_A\Xi_1\bigr).}
\]
Explicitly,
\[
P_{20}=\bar P_0(1+\epsilon\Xi_1),
\qquad
P_{21}=\bar P_0\Bigl(1+\frac12\epsilon\Xi_1\Bigr),
\qquad
P_{22}=\bar P_0(1-\epsilon\Xi_1).
\]

The grouped trace/anomaly compiler then gives
\[
\boxed{a_{P_0}=\frac{\epsilon\bar P_0\Xi_1}{4},}
\qquad
\boxed{b_{P_0}=\frac{3\epsilon\bar P_0\Xi_1}{4}.}
\]
So the exact weak-axisymmetric branch law is
\[
\boxed{b_{P_0}=3a_{P_0}.}
\]

This is important because it means the actual same-charge window is not generic in the full `(a_{P_0},b_{P_0})` plane. On the weak-axisymmetric branch it collapses to a one-dimensional axisymmetric line.

---

## 5. Exact transported ceiling test in terms of `\Xi_1`

On the axisymmetric weak-anisotropy line,
\[
P_{20}=\bar P_0+4a_{P_0},
\qquad
P_{21}=\bar P_0+2a_{P_0},
\qquad
P_{22}=\bar P_0-4a_{P_0}.
\]
So the worst surviving lane is always the one with the larger sign of `a_{P_0}`. Therefore the exact robust all-lane ceiling test collapses to
\[
\boxed{
\bar P_0+4|a_{P_0}|
\le
P_{\rm crit}.
}
\]
Equivalently, in the one-scalar `\Xi_1` language,
\[
\boxed{
\bar P_0\bigl(1+|\epsilon\Xi_1|\bigr)
\le
P_{\rm crit}.
}
\]
Using the exact normalization compiler,
\[
\boxed{
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}igl(1+|\epsilon\Xi_1|\bigr)
\le
P_{\rm crit}.
}
\]

This is the cleanest same-charge continuation reached so far.

The actual moving-throat branch now lives or dies, at this transported level, by one explicit inequality in

- the normalization defect `\Delta_{\rm norm}`,
- the source-map factor `\hat m_0`,
- and the weak-axisymmetric outgoing-prefactor slope `\Xi_1=P_1/P_0`.

### 5.1 Exact calibrated-branch bound with weak anisotropy included

If the actual branch is exactly calibrated,
\[
\Delta_{\rm norm}=0,
\]
then the robust transported ceiling becomes
\[
\boxed{
\hat m_0^{\,2}
\ge
\frac{T_{\rm quad}(1+|\epsilon\Xi_1|)}{P_{\rm crit}}.
}
\]
So weak-axisymmetric prefactor loading raises the lower bound on the source-map factor linearly in the absolute outgoing slope.

That is a very concrete actual-branch falsification test.

---

## 6. Explicit headroom at the Stage-006 compatibility point

Stage 006 found the concrete compatibility point
\[
\bar P_0 = P_{0,\rm target,compat}
\approx 0.002069792318062885.
\]
Substituting this into the robust weak-axisymmetric ceiling law
\[
|\epsilon\Xi_1|\le \frac{P_{\rm crit}}{\bar P_0}-1,
\qquad
|a_{P_0}|\le \frac{P_{\rm crit}-\bar P_0}{4},
\]
gives the following explicit budgets.

### 6.1 Stricter `10%`-loss benchmark

For **both wall-like poles** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 0.367930328492646,
\qquad
|a_{P_0}|
\lesssim 1.90384841874108\times 10^{-4}.
\]
For a **nonempty wall-like corridor** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 0.737619063660757,
\qquad
|a_{P_0}|
\lesssim 3.81679567905443\times 10^{-4}.
\]

### 6.2 Looser `30%`-loss benchmark

For **both wall-like poles** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 2.94889585703134,
\qquad
|a_{P_0}|
\lesssim 1.52590049791274\times 10^{-3}.
\]
For a **nonempty wall-like corridor** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 4.63505472371892,
\qquad
|a_{P_0}|
\lesssim 2.39840016523863\times 10^{-3}.
\]

So the Stage-006 compatibility point still has finite weak-axisymmetric headroom.
But the stricter `10%` robust budget is not huge. That is exactly the sort of narrow corridor we would expect if the idea is real but difficult.

---

## 7. What Stage 007 changes

Stage 007 does not yet prove success or failure.
But it does convert the remaining ambiguity into an actual-branch compiler.

Before this stage, the next step was still phrased loosely as

> “extract the actual branch-compatible normalization target from the PDE.”

After this stage, we know the exact thing to compute:

1. the real branch packet gives `\Delta_{\rm norm}`, `a_{P_0}`, and `b_{P_0}`,
2. those compile to the actual lane prefactors `P_{20},P_{21},P_{22}`,
3. on the weak-axisymmetric branch they collapse to one scalar `\Xi_1=P_1/P_0`,
4. and the same-charge corridor survives only if those data satisfy one explicit finite inequality.

So the best current summary is:

> the idea is still alive, but the actual moving-throat branch now has to land inside a sharply delimited corridor in `(\Delta_{\rm norm},\Xi_1)` or equivalently `(\Delta_{\rm norm},a_{P_0},b_{P_0})`.

That is the first genuine PDE-selected same-charge kill test.
