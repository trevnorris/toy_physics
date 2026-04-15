# Moving-Throat PDE — Stage 212: Selected-Branch Numerator/Denominator Signature and the Softening-Depth Crossover Theorem

## Status

**Exact within the explicit finite-throat D/N selected-branch closure built on the Stage-024/025 selected-mode normalization product and compared against the Stage-211 pure-transfer numerator/denominator split.**

This stage does **not** solve the full moving-throat PDE.
It takes the two exact rigid same-charge survivors from Stage 211,

- **numerator-rigid** \((\pi_1=0)\),
- **denominator-rigid** \((\delta_1=0)\),

and asks the next natural microscopic question:

> once the actual selected quadrupole branch is written in its exact finite-throat form, which of those two rigid signatures does the real PDE-selected mixed branch most closely resemble?

---

## Purpose

Stage 211 already reduced the surviving same-charge pure-transfer corridor to the exact static split
\[
\Xi_1 = 2(\pi_1-\delta_1),
\qquad
\pi_1:=\frac{P_{01}}{P},
\qquad
\delta_1:=\frac{\Delta_{01}}{\Delta},
\]
and then isolated the two rigid 1D subcorridors

- numerator-rigid \((\pi_1=0)\),
- denominator-rigid \((\delta_1=0)\).

But the later moving-throat selected-branch chain already replaced the free one-port static factor by the exact selected-mode normalization product
\[
N_-(x)=\frac{\beta_0\,s_-(x)^2}{\kappa_0^2\,(A-x)},
\qquad 0\le x < A,
\]
with exact selected overlap
\[
s_-(x)=
\frac{\bigl[\kappa_0^2(x+\Delta K_{\rm ax})+\kappa_1^2x\bigr]^2}
{\kappa_0^2(x+\Delta K_{\rm ax})^2+\kappa_1^2x^2}.
\]

So the next honest question is now sharper:

> when the actual selected branch is written in its exact softening-depth language, is it numerator-like, denominator-like, or something genuinely different from both rigid Stage-211 subcorridors?

The main outputs of this stage are:

1. the exact factorization of the selected-branch normalization product into a **numerator-like** and a **denominator-like** piece;
2. the exact log-slope classifier that decides which rigid Stage-211 signature the selected branch most closely resembles at any point on the stable branch;
3. the universal crossover theorem in the dimensionless variables
   \[
   \xi=\frac{x}{A},
   \qquad
   \delta=\frac{\Delta K_{\rm ax}}{A};
   \]
4. the exact threshold
   \[
   \delta=\frac89
   \]
   separating always-denominator branches from branches that begin numerator-dominant;
5. and the sharp conclusion that the actual selected branch is **never** literally one of the rigid Stage-211 subcorridors — it is an exact co-loading product — but it becomes unambiguously **denominator-like** near softening, and for all \(\delta\ge 8/9\) it is denominator-like on the entire stable branch.

So after Stage 212 the continuation point is no longer “numerator-rigid or denominator-rigid?” in the abstract.
It is:

> what are the actual selected-branch ratios \((\xi,\delta)\) on the physical moving-throat branch, and where do they land on this universal classifier map?

---

## 1. Frozen input carried forward

### 1.1 Stage-211 static split

Stage 211 isolated the exact static pure-transfer identity
\[
\Xi_1 = 2(\pi_1-\delta_1),
\qquad
\pi_1:=\frac{P_{01}}{P},
\qquad
\delta_1:=\frac{\Delta_{01}}{\Delta},
\]
and then split the two-dimensional pure-transfer corridor into

- the 1D numerator-rigid branch \(\pi_1=0\),
- the 1D denominator-rigid branch \(\delta_1=0\).

That gave the right static classifiers, but it was still a classifier internal to the primitive mixed-slope space.

### 1.2 Actual selected-branch product

The later moving-throat selected-branch chain already replaced the free static load-factor language by the exact selected-mode normalization product
\[
N_-(x)=\frac{\beta_0\,s_-(x)^2}{\kappa_0^2\,(A-x)},
\]
with exact D/N constants
\[
\kappa_0^2=\frac{8}{\pi^2},
\qquad
\kappa_1^2=\frac{16}{9\pi^2}.
\]

The point of Stage 212 is to compare the Stage-211 numerator/denominator split against **this** actual selected-branch object rather than against a free one-port static factor.

---

## 2. Exact dimensionless selected-branch factorization

Introduce the dimensionless stable-branch variables
\[
\xi:=\frac{x}{A},
\qquad
\delta:=\frac{\Delta K_{\rm ax}}{A},
\qquad
0\le \xi < 1,
\qquad
\delta>0.
\]

With
\[
x=A\xi,
\qquad
\Delta K_{\rm ax}=A\delta,
\]
the selected-branch normalization product factors exactly as
\[
N_-(x)=\frac{8\beta_0}{\pi^2 A}\,F(\xi,\delta),
\]
where
\[
\boxed{
F(\xi,\delta)
=
\frac{(9\delta+11\xi)^4}
{81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2}.
}
\]

So the same universal selected-branch function \(F\) still controls the normalization product; the extra factor \(8/\pi^2=\kappa_0^2\) is a fixed D/N overlap constant and does not affect the classifier below.

Now split this into
\[
\boxed{
F(\xi,\delta)=F_{\rm num}(\xi,\delta)\,F_{\rm den}(\xi),
}
\]
with
\[
\boxed{
F_{\rm num}(\xi,\delta)
=
\frac{(9\delta+11\xi)^4}
{81(9\delta^2+18\delta\xi+11\xi^2)^2},
}
\qquad
\boxed{
F_{\rm den}(\xi)=\frac{1}{1-\xi}.
}
\]

This is the key exact factorization.

It says:

- the selected-branch **numerator-like** gain is the overlap / source-map / internal-transfer factor \(F_{\rm num}\),
- the selected-branch **denominator-like** gain is the explicit softening factor \((1-\xi)^{-1}\).

So the actual PDE-selected branch already comes with a built-in numerator/denominator split.
But unlike Stage 211, the split is not “either/or.” It is an exact product.

---

## 3. Exact log-slope classifier

To decide which rigid Stage-211 branch the selected branch most closely resembles, the right invariant is the log-slope split of \(F\) along the physical softening coordinate \(\xi\).

Define
\[
L_{\rm num}:=\partial_\xi\ln F_{\rm num},
\qquad
L_{\rm den}:=\partial_\xi\ln F_{\rm den},
\qquad
L_{\rm tot}:=\partial_\xi\ln F=L_{\rm num}+L_{\rm den}.
\]

The exact derivatives are
\[
\boxed{
L_{\rm num}(\xi,\delta)
=
\frac{72\delta^2}
{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)},
}
\]
\[
\boxed{
L_{\rm den}(\xi)=\frac{1}{1-\xi}.
}
\]

So the exact selected-branch numerator/denominator classifier is
\[
\boxed{
\mathcal R_{ND}(\xi,\delta)
:=
\frac{L_{\rm num}}{L_{\rm den}}
=
\frac{72\delta^2(1-\xi)}{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}.
}
\]

Interpretation:

- \(\mathcal R_{ND}>1\): the selected branch is **numerator-like** at that point;
- \(\mathcal R_{ND}<1\): it is **denominator-like**;
- \(\mathcal R_{ND}=1\): exact crossover.

This is the Stage-212 replacement for the rigid Stage-211 subcorridors.
It is no longer a statement about a free primitive mixed slope.
It is a statement about the actual selected normalization product.

---

## 4. Exact onset and near-softening limits

### 4.1 Onset

At zero softening,
\[
\xi=0,
\]
the classifier is
\[
\boxed{
\mathcal R_{ND}(0,\delta)=\frac{8}{9\delta}.
}
\]

So the selected branch begins

- numerator-like if \(0<\delta<8/9\),
- exactly balanced if \(\delta=8/9\),
- denominator-like if \(\delta>8/9\).

### 4.2 Near softening

As \(\xi\to1^-\),
\[
L_{\rm den}(\xi)=\frac{1}{1-\xi}\to+\infty,
\]
while
\[
L_{\rm num}(\xi,\delta)\to
\frac{72\delta^2}{(9\delta+11)(9\delta^2+18\delta+11)},
\]
which is finite.

Therefore
\[
\boxed{
\lim_{\xi\to1^-}\mathcal R_{ND}(\xi,\delta)=0.
}
\]

So the actual selected branch is always denominator-like sufficiently close to the softening edge.
This is already a strong answer to the Stage-211 continuation question:

> whatever the selected branch does near onset, it becomes denominator-like before softening.

---

## 5. Exact crossover theorem

The numerator/denominator crossover condition is
\[
\mathcal R_{ND}(\xi,\delta)=1.
\]
Equivalently,
\[
L_{\rm num}=L_{\rm den}.
\]

Clearing denominators gives the exact cubic
\[
\boxed{
\mathcal P(\xi,\delta)
=
121\xi^3+297\delta\xi^2+333\delta^2\xi+81\delta^3-72\delta^2
=0.
}
\]

The derivative is
\[
\boxed{
\partial_\xi\mathcal P
=
363\xi^2+594\delta\xi+333\delta^2 > 0
}
\]
for every \(\xi\ge0\), \(\delta>0\).

So \(\mathcal P\) is strictly increasing in \(\xi\).
That yields the exact theorem.

### 5.1 Always-denominator regime

If
\[
\delta\ge \frac89,
\]
then
\[
\mathcal P(0,\delta)=9\delta^2(9\delta-8)\ge0.
\]
Since \(\mathcal P\) is strictly increasing, it follows that
\[
\mathcal P(\xi,\delta)>0
\qquad
\text{for all }0<\xi<1,
\]
so
\[
\boxed{
\delta\ge\frac89
\quad\Longrightarrow\quad
\mathcal R_{ND}(\xi,\delta)<1
\ \text{for the entire stable branch.}
}
\]

This means:

> if the physical axial gap ratio satisfies \(\delta\ge 8/9\), the selected PDE branch is denominator-like from the start.

### 5.2 Mixed regime

If
\[
0<\delta<\frac89,
\]
then
\[
\mathcal P(0,\delta)<0,
\qquad
\lim_{\xi\to1^-}\mathcal P(\xi,\delta)>0.
\]
Because \(\mathcal P\) is strictly increasing, there exists a **unique**
\[
\xi_*(\delta)\in(0,1)
\]
such that
\[
\mathcal P(\xi_*,\delta)=0.
\]

So
\[
\boxed{
0<\delta<\frac89
\quad\Longrightarrow\quad
\begin{cases}
\mathcal R_{ND}>1,& 0\le \xi<\xi_*(\delta),\\[4pt]
\mathcal R_{ND}=1,& \xi=\xi_*(\delta),\\[4pt]
\mathcal R_{ND}<1,& \xi_*(\delta)<\xi<1.
\end{cases}
}
\]

This is the exact universal crossover theorem.
It says the actual selected branch interpolates from numerator-like to denominator-like whenever the axial gap ratio is sufficiently small.

---

## 6. Sample crossover depths

The exact crossover root \(\xi_*(\delta)\) is algebraic but not especially transparent in radicals, so the most useful quick reading is numerical.

For a few representative gap ratios:

- \(\delta=\frac14\):
  \[
  \xi_*\approx 0.107223051105697;
  \]
- \(\delta=\frac12\):
  \[
  \xi_*\approx 0.081847937860074;
  \]
- \(\delta=\frac34\):
  \[
  \xi_*\approx 0.032505121082825.
  \]

So even when the selected branch begins numerator-like, that window is usually quite short.
The denominator-like regime takes over early.

That is important for the same-charge audit, because the actual normalization hit is not expected at infinitesimal loading.
It happens deeper on the stable branch, where denominator dominance is the natural expectation.

---

## 7. What this says about the Stage-211 rigid subcorridors

Stage 211 asked which rigid static signature the real PDE-selected mixed branch most closely resembles.
Stage 212 gives the precise answer.

### 7.1 It is not literally either rigid subcorridor

The selected branch carries
\[
F(\xi,\delta)=F_{\rm num}(\xi,\delta)\,F_{\rm den}(\xi),
\]
so both factors move simultaneously.
Therefore the actual selected branch is **not** literally numerator-rigid or denominator-rigid.

It is an exact co-loading branch.

### 7.2 It becomes denominator-like near the physical target window

Because the denominator factor diverges and the numerator factor stays finite as \(\xi\to1^-\), the selected branch always becomes denominator-like sufficiently near softening.

So if the physical branch hits the universal target at appreciable softening depth, the right Stage-211 proxy is the denominator-rigid one, not the numerator-rigid one.

### 7.3 Only very early on can it look numerator-like

Numerator-like behavior is confined to the small-softening regime
\[
0\le \xi < \xi_*(\delta),
\qquad
0<\delta<\frac89.
\]

So the numerator-rigid Stage-211 branch is best read as an **onset-side local proxy**, not as the global selected-branch signature.

---

## 8. Best current verdict after Stage 212

The continuation question from Stage 211 now has a clean answer.

1. The real selected PDE branch is not one of the rigid Stage-211 subcorridors.
   It is an exact numerator/denominator **co-loading** product.
2. The exact classifier is
   \[
   \mathcal R_{ND}(\xi,\delta)
   =
   \frac{72\delta^2(1-\xi)}
   {(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}.
   \]
3. If \(\delta\ge 8/9\), the selected branch is denominator-like on the whole stable branch.
4. If \(0<\delta<8/9\), the selected branch begins numerator-like but crosses uniquely to denominator-like at \(\xi=\xi_*(\delta)\).
5. Near softening — and therefore near any large selected-branch normalization gain — the selected branch is always denominator-like.

So the next honest stage is now very specific:

> feed the actual moving-throat selected-branch data into \((\xi,\delta)\), place the physical branch on this universal classifier map, and then compare the resulting denominator-vs-numerator signature against the concrete Stage-211 dynamic ceilings.

That is the clean continuation point.

---

## 9. SymPy-backed status

The accompanying audit verifies:

- the exact reduction of the selected-mode product \(N_-(x)\) to the universal dimensionless function \(F(\xi,\delta)\);
- the exact factorization \(F=F_{\rm num}F_{\rm den}\);
- the exact log-slope formulas for \(L_{\rm num}\), \(L_{\rm den}\), and \(\mathcal R_{ND}\);
- the exact onset law \(\mathcal R_{ND}(0,\delta)=8/(9\delta)\);
- the near-softening limit \(\lim_{\xi\to1^-}\mathcal R_{ND}=0\);
- the exact crossover cubic \(\mathcal P(\xi,\delta)=0\) and the positivity of \(\partial_\xi\mathcal P\);
- the always-denominator threshold \(\delta=8/9\);
- and the sample crossover depths for \(\delta=1/4,1/2,3/4\).

Supporting file:
- `moving_throat_pde_stage212_selected_branch_numerator_denominator_signature_and_softening_depth_crossover_theorem_sympy_audit.py`

---

## 10. Immediate next step

The next clean move is now well defined.

1. Keep the exact Stage-212 classifier map.
2. Insert the actual selected-branch moving-throat ratios \((\xi,\delta)\).
3. Decide whether the physical branch sits in the numerator-like onset window or the denominator-like deeper-softening window.
4. Then compile that placement into the carried Stage-211 wall-like dynamic ceilings.

That is the smallest next theorem gate that directly connects the selected-branch normalization geometry to the already-audited same-charge dynamic windows.
