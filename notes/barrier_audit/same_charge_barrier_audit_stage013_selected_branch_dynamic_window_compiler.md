# Same-Charge Barrier Audit — Stage 013: Selected-Branch Classifier-to-Dynamic Window Compiler and the Static-First Theorem

## 0. Purpose

Stage 012 gave the exact selected-branch numerator/denominator classifier
\[
\mathcal R_{ND}(\xi,\delta)
=
\frac{72\delta^2(1-\xi)}{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)},
\]
and Stage 011 gave the first actual wall-like dynamic responses of the two rigid pure-transfer survivors.

So the next honest question is now completely sharp:

> if the **actual selected branch** is a numerator/denominator co-loading product rather than one rigid split, does its wall-like dynamic window ever become the first kill condition, or is the first real ceiling still the transported static `\Xi_1` budget?

This stage answers that inside the exact rigid-split compiler built from Stages 011 and 012.

The main outputs are:

1. the exact selected-branch share weights
   \[
   w_{\rm num}=\frac{\mathcal R_{ND}}{1+\mathcal R_{ND}},
   \qquad
   w_{\rm den}=\frac{1}{1+\mathcal R_{ND}},
   \]
   which split the selected branch into its numerator-carried and denominator-carried parts;
2. the exact selected-branch wall-like dynamic slopes per unit `\Xi_1` as affine mixtures of the carried Stage-011 rigid-branch slopes;
3. the exact sign theorem that the upper wall-like pole always worsens, while the lower wall-like pole flips sign only at one finite classifier threshold
   \[
   \mathcal R_*\approx 1.229255438463336;
   \]
4. the associated onset threshold
   \[
   \delta_*^{(\rm dyn)}=\frac{8}{9\mathcal R_*}\approx 0.723111617875019;
   \]
5. and the central verdict that the **dynamic ceilings are everywhere weaker than the universal transported static ceilings**:
   the first kill condition on the selected branch is still the static `\Xi_1` budget, not the wall-like dynamic window.

So after Stage 013, the continuation point is no longer to wonder whether the selected-branch dynamic window kills the same-charge corridor. It does not, at least inside this exact rigid-split compiler on the concrete compatibility branch. The remaining question is still the static placement of `\Xi_1` on the actual moving-throat branch.

---

## 1. Frozen input carried forward

### 1.1 Stage-011 rigid dynamic data

Stage 011 isolated two exact 1D pure-transfer survivors on the concrete compatibility branch.

The first is the **numerator-rigid** branch
\[
\pi_1=0,
\]
so the same-charge signal is carried entirely by the denominator:
\[
\Xi_1=-2\delta_1.
\]
Its carried Stage-011 unit direction has
\[
\Xi_1\approx 1.73611235,
\]
and wall-like dynamic slopes
\[
\delta\ln \mathcal R_{Q,+}\approx -0.52346582,
\qquad
\delta\ln \mathcal R_{Q,-}\approx +0.71358484.
\]
So per unit `\Xi_1`, the denominator-carried dynamic slopes are
\[
\boxed{
 s_{+}^{(\rm den)}\approx -0.301516097158113,
 \qquad
 s_{-}^{(\rm den)}\approx +0.411024574532864.
}
\]

The second is the **denominator-rigid** branch
\[
\delta_1=0,
\]
so the same-charge signal is carried entirely by the numerator:
\[
\Xi_1=2\pi_1.
\]
Its carried Stage-011 unit direction has
\[
\Xi_1\approx 0.69293215,
\]
and wall-like dynamic slopes
\[
\delta\ln \mathcal R_{Q,+}\approx -0.35245541,
\qquad
\delta\ln \mathcal R_{Q,-}\approx -0.23169484.
\]
So per unit `\Xi_1`, the numerator-carried dynamic slopes are
\[
\boxed{
 s_{+}^{(\rm num)}\approx -0.508643465308977,
 \qquad
 s_{-}^{(\rm num)}\approx -0.334368725711457.
}
\]

These four numbers are the only dynamic inputs needed in Stage 013.

### 1.2 Stage-012 selected-branch classifier

Stage 012 showed that the actual selected branch is not literally numerator-rigid or denominator-rigid. It is an exact co-loading product with classifier
\[
\mathcal R_{ND}(\xi,\delta)
=
\frac{72\delta^2(1-\xi)}{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}.
\]

The derivative is
\[
\partial_\xi \mathcal R_{ND}
=
-
\frac{72\delta^2\Bigl(81\delta^3+261\delta^2+297\delta\xi(2-\xi)+\xi^2(363-242\xi)\Bigr)}{(9\delta+11\xi)^2(9\delta^2+18\delta\xi+11\xi^2)^2},
\]
which is strictly negative on the stable interval
\[
0\le \xi<1,
\qquad
\delta>0.
\]
So the classifier decreases monotonically along the selected branch.

That monotonicity is the key input that lets us turn Stage-012 signature data into exact Stage-013 ceiling statements.

---

## 2. Exact rigid-split share compiler

The selected branch has numerator-like and denominator-like log-slope contributions
\[
L_{\rm num},\qquad L_{\rm den},
\qquad
\mathcal R_{ND}=\frac{L_{\rm num}}{L_{\rm den}}.
\]
So its exact contribution shares are
\[
\boxed{
 w_{\rm num}=\frac{L_{\rm num}}{L_{\rm num}+L_{\rm den}}=
 \frac{\mathcal R_{ND}}{1+\mathcal R_{ND}},
 \qquad
 w_{\rm den}=\frac{L_{\rm den}}{L_{\rm num}+L_{\rm den}}=
 \frac{1}{1+\mathcal R_{ND}}.
}
\]

Inside the exact rigid-split compiler, these are the weights with which the selected branch samples the carried Stage-011 dynamic responses.

So the selected-branch wall-like dynamic slopes per unit `\Xi_1` are
\[
\boxed{
 S_+(\mathcal R_{ND})
 =
 \frac{\mathcal R_{ND}s_+^{(\rm num)}+s_+^{(\rm den)}}{1+\mathcal R_{ND}},
}
\]
\[
\boxed{
 S_-(\mathcal R_{ND})
 =
 \frac{\mathcal R_{ND}s_-^{(\rm num)}+s_-^{(\rm den)}}{1+\mathcal R_{ND}}.
}
\]

These are exact affine mixtures of the rigid-branch per-unit-`\Xi_1` slopes.

This is the first Stage-013 compression:

> once the Stage-012 classifier is known, the selected-branch dynamic response is completely determined inside the rigid-split compiler.

---

## 3. Exact sign theorem for the selected wall-like poles

### 3.1 Upper wall-like pole

Both carried upper-pole slopes are negative:
\[
 s_+^{(\rm num)}<0,
 \qquad
 s_+^{(\rm den)}<0.
\]
So for every classifier value
\[
\boxed{S_+(\mathcal R_{ND})<0.}
\]

The upper wall-like pole always worsens.

### 3.2 Lower wall-like pole

The carried lower-pole slopes have opposite sign:
\[
 s_-^{(\rm num)}<0,
 \qquad
 s_-^{(\rm den)}>0.
\]
So the lower wall-like pole flips sign exactly once, at
\[
S_-(\mathcal R_*)=0.
\]
This gives the exact threshold
\[
\boxed{
\mathcal R_*=rac{s_-^{(\rm den)}}{-s_-^{(\rm num)}}
\approx 1.229255438463336.
}
\]

Therefore
\[
\boxed{
\mathcal R_{ND}<\mathcal R_*
\Longrightarrow
S_-(\mathcal R_{ND})>0,
}
\]
\[
\boxed{
\mathcal R_{ND}=\mathcal R_*
\Longrightarrow
S_-(\mathcal R_{ND})=0,
}
\]
\[
\boxed{
\mathcal R_{ND}>\mathcal R_*
\Longrightarrow
S_-(\mathcal R_{ND})<0.
}
\]

So the selected branch has a universal dynamic-sign split:

- if the classifier is not too numerator-dominant, the lower wall-like pole actually **improves**;
- only once the selected branch is strongly numerator-dominant do both wall-like poles worsen.

This is already enough to show that the denominator-like part of the classifier map is dynamically safer than the numerator-like part.

---

## 4. Immediate consequences for the Stage-012 classifier map

### 4.1 Every denominator-like point has infinite nonempty dynamic ceiling

Stage 012 already proved that denominator-like means
\[
\mathcal R_{ND}\le 1.
\]
Since
\[
1<\mathcal R_*,
\]
we get immediately
\[
\boxed{
\mathcal R_{ND}\le 1
\Longrightarrow
S_-(\mathcal R_{ND})>0.
}
\]

So every denominator-like selected-branch point has the same split-sign dynamic response as the Stage-011 denominator-carried rigid branch:

- upper wall-like pole worsens,
- lower wall-like pole improves.

That means its **nonempty** dynamic ceiling is infinite.

In particular, the whole always-denominator regime from Stage 012,
\[
\delta\ge \frac89,
\]
inherits an infinite nonempty dynamic ceiling on the whole stable branch.

### 4.2 A stronger onset threshold

At onset,
\[
\mathcal R_{ND}(0,\delta)=\frac{8}{9\delta}.
\]
Requiring onset to stay below the sign-flip threshold gives
\[
\frac{8}{9\delta}\le \mathcal R_*,
\]
so
\[
\boxed{
\delta\ge \delta_*^{(\rm dyn)}:=\frac{8}{9\mathcal R_*}
\approx 0.723111617875019.
}
\]

Because `\mathcal R_{ND}` decreases monotonically with `\xi`, every branch with
\[
\delta\ge \delta_*^{(\rm dyn)}
\]
stays below `\mathcal R_*` on the whole stable interval. Therefore
\[
\boxed{
\delta\ge 0.723111617875019
\Longrightarrow
\text{nonempty dynamic ceiling is infinite on the entire selected branch.}
}
\]

This is stronger than the Stage-012 denominator-like theorem. It says even a substantial subset of the onset-side numerator-like branches still never lose their nonempty dynamic window.

---

## 5. Exact dynamic ceilings in `|\epsilon\Xi_1|`

Use the carried Stage-011 wall-like dynamic figures
\[
\mathcal R_{Q,-}\approx 30.199907560250075,
\qquad
\mathcal R_{Q,+}\approx 36.171186483269487,
\]
and the same stricter `10%`-loss requirement
\[
\mathcal R_{Q,\rm req}\approx 21.8545662963584.
\]
Define the dynamic margins
\[
\ell_-:=\ln\frac{\mathcal R_{Q,-}}{\mathcal R_{Q,\rm req}}
\approx 0.323428979934714,
\]
\[
\ell_+:=\ln\frac{\mathcal R_{Q,+}}{\mathcal R_{Q,\rm req}}
\approx 0.503852964869151.
\]

Then the selected-branch **robust** dynamic ceiling on `|\epsilon\Xi_1|` is
\[
\boxed{
B_{\rm dyn}^{(\rm both)}(\mathcal R_{ND})
=
\min\!\left(
\frac{\ell_+}{-S_+(\mathcal R_{ND})},
\frac{\ell_-}{-S_-(\mathcal R_{ND})}
\right),
}
\]
with the second term understood as `+\infty` whenever `S_-\ge 0`.

The **nonempty** dynamic ceiling is
\[
\boxed{
B_{\rm dyn}^{(\rm nonempty)}(\mathcal R_{ND})
=
\begin{cases}
+\infty, & S_-(\mathcal R_{ND})\ge 0,\\[4pt]
\max\!\left(
\dfrac{\ell_+}{-S_+(\mathcal R_{ND})},
\dfrac{\ell_-}{-S_-(\mathcal R_{ND})}
\right), & S_-(\mathcal R_{ND})<0.
\end{cases}
}
\]

The endpoint values are already enough to understand the whole story:
\[
\boxed{
B_{\rm dyn}^{(\rm both)}(0)
\approx 1.671064893775584,
}
\]
\[
\boxed{
\lim_{\mathcal R_{ND}\to\infty} B_{\rm dyn}^{(\rm both)}
\approx 0.967282389363822,
}
\]
\[
\boxed{
\lim_{\mathcal R_{ND}\to\infty} B_{\rm dyn}^{(\rm nonempty)}
\approx 0.990581810705233.
}
\]

So even the worst selected-branch robust dynamic ceiling is still close to one full unit of `|\epsilon\Xi_1|`.

That is already much looser than the transported static budgets below.

---

## 6. Universal transported static ceilings in `|\epsilon\Xi_1|`

Stage 011 gave the transported static ceilings in the rigid branch parameter `t`. Converting them to `|\epsilon\Xi_1|` by multiplying with the carried rigid-branch `\Xi_1` values yields the same universal numbers from both rigid splits:
\[
\boxed{
B_{\rm stat}^{(\rm both)}\approx 0.367930328492646,
}
\]
\[
\boxed{
B_{\rm stat}^{(\rm nonempty)}\approx 0.737619063660757.
}
\]

This universality is expected: on the pure-transfer corridor,
\[
\delta\ln P_0 = \Xi_1,
\]
so the transported Stage-007/011 static budgets naturally live in the branch-invariant variable `|\epsilon\Xi_1|`.

Now compare the worst selected-branch dynamic ceilings with these universal static budgets:
\[
\boxed{
\inf_{\mathcal R_{ND}\ge 0} B_{\rm dyn}^{(\rm both)}
\approx 0.967282389363822
>
0.367930328492646
= B_{\rm stat}^{(\rm both)}.
}
\]
\[
\boxed{
\inf_{\mathcal R_{ND}\ge 0,\;B_{\rm dyn}^{(\rm nonempty)}<\infty}
B_{\rm dyn}^{(\rm nonempty)}
\approx 0.990581810705233
>
0.737619063660757
= B_{\rm stat}^{(\rm nonempty)}.
}
\]

So the selected-branch dynamic window is **everywhere weaker** than the universal transported static ceiling.

This is the central Stage-013 theorem.

---

## 7. Sample classifier points

For a few representative classifier values:

### 7.1 Pure denominator-carried limit `\mathcal R_{ND}=0`
\[
S_+\approx -0.301516097158113,
\qquad
S_-\approx +0.411024574532864,
\]
so
\[
B_{\rm dyn}^{(\rm both)}\approx 1.671064893775584,
\qquad
B_{\rm dyn}^{(\rm nonempty)}=+\infty.
\]

### 7.2 Exact numerator/denominator balance `\mathcal R_{ND}=1`
\[
S_+\approx -0.405079781233545,
\qquad
S_-\approx +0.038327924410703,
\]
so
\[
B_{\rm dyn}^{(\rm both)}\approx 1.243836370541187,
\qquad
B_{\rm dyn}^{(\rm nonempty)}=+\infty.
\]

### 7.3 Sign-flip threshold `\mathcal R_{ND}=\mathcal R_*`
\[
S_- = 0,
\]
so the lower wall-like pole is exactly neutral and the nonempty dynamic ceiling is still infinite.

### 7.4 Strong numerator-like point `\mathcal R_{ND}=10`
\[
S_+\approx -0.489813704567990,
\qquad
S_-\approx -0.266605698416519,
\]
so
\[
B_{\rm dyn}^{(\rm both)}\approx 1.028662448947899,
\qquad
B_{\rm dyn}^{(\rm nonempty)}\approx 1.213136035184892.
\]

Even here, the robust dynamic ceiling is still far above the universal static robust budget.

---

## 8. Best current verdict after Stage 013

Stage 013 does not kill the same-charge corridor.
It sharpens the verdict instead.

Inside the exact rigid-split compiler built from Stages 011 and 012 on the concrete compatibility branch:

1. the selected-branch wall-like dynamic response is completely controlled by the single classifier `\mathcal R_{ND}`;
2. the upper wall-like pole always worsens, but the lower one improves whenever
   \[
   \mathcal R_{ND}\le \mathcal R_*\approx 1.229255438463336;
   \]
3. every denominator-like point therefore has infinite nonempty dynamic ceiling;
4. if
   \[
   \delta\ge \frac89,
   \]
   the whole selected branch is denominator-like and hence has infinite nonempty dynamic ceiling on the whole stable interval;
5. even more strongly, if
   \[
   \delta\ge 0.723111617875019,
   \]
   the entire selected branch stays below the sign-flip threshold `\mathcal R_*`, so its nonempty dynamic ceiling is still infinite everywhere;
6. and for **all** selected-branch signatures,
   \[
   B_{\rm dyn}^{(\rm both)} > B_{\rm stat}^{(\rm both)},
   \qquad
   B_{\rm dyn}^{(\rm nonempty)} > B_{\rm stat}^{(\rm nonempty)}.
   \]

So the first kill condition on the selected same-charge branch is still the transported static `\Xi_1` budget, not the wall-like dynamic window.

That is the right continuation point after Stage 013.
