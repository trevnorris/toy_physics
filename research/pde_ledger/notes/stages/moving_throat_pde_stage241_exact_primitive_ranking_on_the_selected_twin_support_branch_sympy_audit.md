# Moving-Throat PDE — Stage 241: Exact Primitive Ranking on the Selected Twin-Support Branch

## Status

**Exact within the carried selected-branch support identities, the support-carried primitive coherent closure, and the constructive coherent bound**
\[
0<\beta<\frac{2}{11}.
\]

This stage does **not** yet place the actual moving-throat branch at a unique point on the selected support curve.
It closes a narrower reduced question:

> once Stage 240 has fixed the passive/outgoing support slice to the symmetric lowest-twin branch, what primitive quartic-carrier hierarchy survives along that exact slice?

The answer is completely sharp.

The selected same-charge branch is the one-parameter curve

\[
\boxed{
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma = \frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23,
}
\]

and on that curve only **two** primitive crossover thresholds remain,

\[
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\qquad
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
\]

So the complete selected-branch primitive ranking is

\[
\boxed{
\begin{aligned}
&0<\varrho<\varrho_{W\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_\Lambda > q_W > |q_U|,\\[1mm]
&\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_W > q_\Lambda > |q_U|,\\[1mm]
&\varrho_{U\Lambda}<\varrho<\frac23
&&\Longrightarrow&&
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
\end{aligned}
}
\]

So the selected same-charge anomaly branch now has a completely explicit primitive ranking phase diagram.

---

## Purpose

Stage 240 compressed the support-selection side of the passive/outgoing same-charge branch to one exact support ratio:

\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43\,C_{\rm mix}.
\]

Equivalently, the support selector

\[
\varrho:=\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
\]

was no longer free. It was tied to the same selected branch by

\[
\epsilon_* = 1-\frac{3\varrho}{2}.
\]

That killed the old mixed-only and non-twin support branches, but it still left one real microscopic ambiguity:

> along the exact selected twin-support slice, which primitive coherent carriers actually dominate the quartic same-charge correction?

This stage answers that by taking the primitive weight formulas from the coherent branch, pulling them back to the selected twin-support curve, and solving the surviving carrier crossovers exactly.

---

## Provenance

This note is the moving-throat Stage-241 version of the selected-branch primitive-ranking computation carried by the anomaly-side support track.

Conceptually it sits directly after:

- **Stage 240**, which fixed the selected support ratio to the exact lowest-twin slice
  \[
  \rho_\alpha=\frac43,
  \qquad
  \Pi_{\rm tr}=\frac43\,C_{\rm mix},
  \]
- and the carried coherent primitive-weight formulas in the same anomaly ledger.

So the real job of Stage 241 is not to reopen support selection.
It is to determine what primitive hierarchy survives **after** support selection has already collapsed.

---

## 0. Why this stage is needed

Before this step, the selected same-charge branch was already known to lie strictly inside the symmetric lowest-twin support window.
But that still did **not** tell us whether the quartic same-charge repair should be read physically as dominated by

- outgoing-scale motion `q_\Lambda`,
- wall blocking `q_W`,
- overlap/interference `q_Z, q_\chi`,
- or the split-`U` companion `q_U`.

The coherent primitive side had already narrowed the ambiguity substantially:

- `q_\chi` is always above `q_Z`,
- `q_Z` is always above `q_W`,
- `q_W` is always above `|q_U|`,
- `q_Z` is always above `q_\Lambda`,

so the only unresolved comparisons were

\[
q_W \ \text{vs.}\ q_\Lambda,
\qquad
|q_U| \ \text{vs.}\ q_\Lambda.
\]

Stage 241 shows that, on the exact selected twin-support curve, those two comparisons are controlled by only two exact threshold coordinates in `\varrho`.

---

## 1. The selected branch is an exact one-parameter twin-support curve

Stage 240 fixed the selected support ratio to

\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43,
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
\]

So

\[
\varrho
=
\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\frac{2(1-\epsilon_*)}{3}.
\]

Hence the selected branch satisfies the exact support law

\[
\boxed{
\epsilon_* = 1-\frac{3\varrho}{2}.
}
\]

Since `0<\epsilon_*<1`, this gives the exact selected-branch range

\[
\boxed{
0<\varrho<\frac23.
}
\]

Now convert to the blocking variable

\[
\sigma:=\frac{2\epsilon_*}{1-\epsilon_*}.
\]

Substituting the selected-branch law gives

\[
\sigma
=
\frac{2(1-3\varrho/2)}{3\varrho/2}
=
\boxed{
\frac{4}{3\varrho}-2.
}
\]

So the selected branch is not a two-parameter region in `(\epsilon_*,\varrho)`.
It is one exact curve, and all surviving primitive carrier comparisons have to be read along that curve.

---

## 2. The selected curve lies strictly inside the symmetric lowest-twin window

The support windows already carried forward are

\[
0<\sigma\le \frac{1}{\varrho}-2
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
\]

\[
\frac{1}{\varrho}-2 < \sigma \le \frac{2}{\varrho}-2
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
\sigma > \frac{2}{\varrho}-2
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
\]

On the selected branch,

\[
\sigma_{\rm sel}=\frac{4}{3\varrho}-2.
\]

Then

\[
\sigma_{\rm sel} - \left(\frac{1}{\varrho}-2\right)
=
\frac{1}{3\varrho}>0,
\]

and

\[
\left(\frac{2}{\varrho}-2\right)-\sigma_{\rm sel}
=
\frac{2}{3\varrho}>0.
\]

So for every allowed point on the selected branch,

\[
\boxed{
\frac{1}{\varrho}-2
<
\sigma_{\rm sel}
<
\frac{2}{\varrho}-2.
}
\]

Therefore the selected same-charge branch lies **strictly inside** the symmetric lowest-twin regime for all `0<\varrho<2/3`.

So mixed-only and non-twin branches are gone from the live closure before we even start ranking the primitive quartic carriers.

---

## 3. Primitive coherent weights on the selected curve

On the support-carried minimum-norm coherent closure, introduce the positive normalization factor

\[
N(\epsilon_*,\beta)
=
5(1-\epsilon_*)^2 + 6\epsilon_*^2(1+\beta^2)>0.
\]

The primitive carrier weights are

\[
\boxed{
w_\Lambda
=
\frac{\epsilon_*^2(1+\beta^2)}{N},
}
\]

\[
\boxed{
w_Z
=
\frac{1-2\epsilon_*+(2+\beta^2)\epsilon_*^2}{N},
}
\]

\[
\boxed{
w_\chi
=
\frac{2\bigl[1-2\epsilon_*+(2+\beta^2)\epsilon_*^2\bigr]}{N},
}
\]

\[
\boxed{
w_W
=
\frac{\epsilon_*(1-\epsilon_*)}{N},
}
\]

\[
\boxed{
|w_U|
=
\frac{\beta\,\epsilon_*(1-\epsilon_*)}{N}.
}
\]

With the common coherent scalar `\Lambda_1` factored out, the primitive amplitudes are

\[
q_\Lambda=\Lambda_1 w_\Lambda,\qquad
q_Z=\Lambda_1 w_Z,\qquad
q_\chi=\Lambda_1 w_\chi,\qquad
q_W=\Lambda_1 w_W,\qquad
q_U=-\Lambda_1|w_U|.
\]

So the primitive hierarchy reduces to the ordering of the weights.

### 3.1 Branch-independent identities

Before imposing the selected-branch curve, the coherent branch already gives the exact identities

\[
\boxed{
w_\chi = 2w_Z,
}
\]

\[
w_Z-w_\Lambda
=
\frac{(1-\epsilon_*)^2}{N}>0,
\]

\[
w_Z-w_W
=
\frac{\beta^2\epsilon_*^2+3(\epsilon_*-1/2)^2+1/4}{N}>0,
\]

\[
w_W-|w_U|
=
\frac{\epsilon_*(1-\epsilon_*)(1-\beta)}{N}>0,
\]

because the constructive coherent branch obeys

\[
\boxed{
0<\beta<\frac{2}{11}<1.
}
\]

So the branch-independent ranking is

\[
\boxed{
q_\chi > q_Z > q_W > |q_U|,
\qquad
q_Z > q_\Lambda.
}
\]

Therefore the only surviving comparisons are

\[
q_W \ \text{vs.}\ q_\Lambda,
\qquad
|q_U| \ \text{vs.}\ q_\Lambda.
\]

Those are exactly the two thresholds that remain on the selected curve.

---

## 4. Surviving threshold 1: `q_W` versus `q_\Lambda`

The exact crossover condition from the coherent primitive side is

\[
q_W=q_\Lambda
\iff
\epsilon_*=\frac{1}{2+\beta^2}.
\]

Insert the selected-branch law

\[
\epsilon_* = 1-\frac{3\varrho}{2}.
\]

Then

\[
1-\frac{3\varrho}{2}=\frac{1}{2+\beta^2},
\]

so the selected-branch threshold is

\[
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)}.
}
\]

Equivalently,

\[
\epsilon_*-\frac{1}{2+\beta^2}
=
\frac32\bigl(\varrho_{W\Lambda}-\varrho\bigr).
\]

So:

- if
  \[
  0<\varrho<\varrho_{W\Lambda},
  \]
  then `q_\Lambda > q_W`;
- if
  \[
  \varrho>\varrho_{W\Lambda},
  \]
  then `q_W > q_\Lambda`.

A useful exact factorization on the selected curve is

\[
\boxed{
w_\Lambda-w_W
=
\frac{(2-3\varrho)(2+\beta^2)(\varrho_{W\Lambda}-\varrho)}
{18\beta^2\varrho^2-24\beta^2\varrho+8\beta^2+33\varrho^2-24\varrho+8}.
}
\]

So the sign transfer is completely explicit:
the only zero occurs at `\varrho=\varrho_{W\Lambda}`.

---

## 5. Surviving threshold 2: `|q_U|` versus `q_\Lambda`

The second coherent primitive crossover is

\[
|q_U|=q_\Lambda
\iff
\epsilon_*=\frac{\beta}{1+\beta+\beta^2}.
\]

Again insert the selected-branch law:

\[
1-\frac{3\varrho}{2}
=
\frac{\beta}{1+\beta+\beta^2}.
\]

Then the second selected-branch threshold is

\[
\boxed{
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
\]

Equivalently,

\[
\epsilon_*-\frac{\beta}{1+\beta+\beta^2}
=
\frac32\bigl(\varrho_{U\Lambda}-\varrho\bigr).
\]

So:

- if
  \[
  \varrho<\varrho_{U\Lambda},
  \]
  then `q_\Lambda > |q_U|`;
- if
  \[
  \varrho>\varrho_{U\Lambda},
  \]
  then `|q_U| > q_\Lambda`.

The exact sign factorization on the selected curve is

\[
\boxed{
w_\Lambda-|w_U|
=
\frac{(2-3\varrho)(1+\beta+\beta^2)(\varrho_{U\Lambda}-\varrho)}
{18\beta^2\varrho^2-24\beta^2\varrho+8\beta^2+33\varrho^2-24\varrho+8}.
}
\]

Again the sign transfer is explicit and controlled by one threshold only.

---

## 6. Ordering and numerical size of the two thresholds

The thresholds are not independent.
Their difference is

\[
\varrho_{U\Lambda}-\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)(1-\beta)}{3(1+\beta+\beta^2)(2+\beta^2)} > 0,
\]

because `0<\beta<2/11<1`.

And

\[
\frac23-\varrho_{U\Lambda}
=
\frac{2\beta}{3(1+\beta+\beta^2)} > 0.
\]

So the exact threshold ordering on the selected branch is

\[
\boxed{
0<\varrho_{W\Lambda}<\varrho_{U\Lambda}<\frac23.
}
\]

This means the selected twin-support curve always splits into **three** ranking regions and never fewer.

### 6.1 Numerical windows from the constructive coherent bound

Because

\[
0<\beta<\frac{2}{11},
\]

the threshold windows are very narrow.

First,

\[
\frac{d\varrho_{W\Lambda}}{d\beta}
=
\frac{4\beta}{3(\beta^2+2)^2}>0,
\]

so `\varrho_{W\Lambda}` increases across the constructive branch.
Therefore

\[
\boxed{
\frac13 < \varrho_{W\Lambda} < \frac{125}{369}\approx 0.338753.
}
\]

Second,

\[
\frac{d\varrho_{U\Lambda}}{d\beta}
=
-\frac{2(1-\beta^2)}{3(1+\beta+\beta^2)^2}<0
\qquad (0<\beta<1),
\]

so `\varrho_{U\Lambda}` decreases across the constructive branch.
Therefore

\[
\boxed{
\frac{250}{441}\approx 0.566893 < \varrho_{U\Lambda} < \frac23.
}
\]

So the selected curve has a very clean geometry:

- only the **low-`\varrho`** end allows `q_\Lambda` to beat `q_W`,
- across the middle of the selected curve, `q_W` beats `q_\Lambda` but `q_\Lambda` still beats `|q_U|`,
- only near the **large-`\varrho` / very weak-blocking** end does `|q_U|` overtake `q_\Lambda`.

---

## 7. Exact primitive ranking theorem on the selected twin-support branch

Combining the branch-independent ordering

\[
q_\chi > q_Z > q_W > |q_U|,
\qquad
q_Z > q_\Lambda,
\]

with the two selected-branch thresholds above gives the complete selected-branch ranking.

### Region I — low `\varrho`, strong blocking

If

\[
0<\varrho<\varrho_{W\Lambda},
\]

then

\[
\boxed{
q_\chi > q_Z > q_\Lambda > q_W > |q_U|.
}
\]

### Region II — intermediate `\varrho`

If

\[
\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda},
\]

then

\[
\boxed{
q_\chi > q_Z > q_W > q_\Lambda > |q_U|.
}
\]

### Region III — large `\varrho`, very weak blocking

If

\[
\varrho_{U\Lambda}<\varrho<\frac23,
\]

then

\[
\boxed{
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
}
\]

So the selected anomaly branch now has a completely explicit primitive ranking phase diagram.

---

## 8. Physical reading

This stage is the cleanest answer the selected support branch can currently support.

The quartic same-charge repair is **not** generically driven by large split-`U` motion.
Instead:

1. `q_\chi` is always the dominant primitive carrier;
2. `q_Z` is always the second-largest carrier;
3. `q_W` dominates over `|q_U|` everywhere;
4. `q_\Lambda` overtakes `q_W` only in the low-`\varrho` / high-blocking corner;
5. and `|q_U|` overtakes `q_\Lambda` only in the large-`\varrho` / very weak-blocking corner.

So the selected same-charge quartic layer is generically an

\[
\text{interference / overlap / wall-blocking / outgoing-scale correction,}
\]

not a large axial split correction.

That is already a strong sharpening of the old broader support-phase picture.

---

## 9. What this stage fixes and what it leaves open

### What is now fixed

- the selected branch is an exact one-parameter twin-support curve;
- the selected branch lies strictly inside the symmetric lowest-twin support window;
- only two primitive crossover thresholds survive on that curve;
- the complete primitive ranking is exact in three regions.

### What is still open

The actual moving-throat branch is still not pinned to one unique value of `\varrho` (equivalently `\epsilon_*`) inside the selected curve.

So the next honest task is now extremely narrow:

> determine the actual physical position of the moving-throat branch on the selected twin-support curve.

Once that single coordinate is known, the quartic same-charge hierarchy stops being a phase diagram and becomes one definite branch verdict.

That is exactly what the next stage has to do.

---

## 10. SymPy-backed status

The accompanying SymPy audit verifies:

1. the exact selected-branch reduction
   \[
   \epsilon_* = 1-\frac{3\varrho}{2},
   \qquad
   \sigma = \frac{4}{3\varrho}-2,
   \]
2. the strict inclusion of the selected curve in the symmetric lowest-twin support window,
3. the exact primitive weight identities
   \[
   w_\chi=2w_Z,
   \qquad
   w_Z>w_\Lambda,\ w_Z>w_W,\ w_W>|w_U|,
   \]
4. the exact selected-branch thresholds `\varrho_{W\Lambda}` and `\varrho_{U\Lambda}`,
5. the factorized sign-transfer laws for `w_\Lambda-w_W` and `w_\Lambda-|w_U|`,
6. the exact ordering
   \[
   0<\varrho_{W\Lambda}<\varrho_{U\Lambda}<\frac23,
   \]
7. the numerical threshold windows implied by `0<\beta<2/11`,
8. and representative region checks confirming all three ranking regimes.

Supporting file:
- `moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py`
