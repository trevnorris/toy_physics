# Moving-Throat PDE — Stage 231: Continuum Placement Pullback of the Selected-Branch Dynamic-Class Map

## Status

**Exact within the carried universal selected-branch D/N geometry, the exact continuum placement map, and the Stage-230 selected-branch classifier-to-dynamic-window compiler.**

This stage does **not** solve the full moving-throat PDE.
It takes the exact Stage-230 dynamic classifier
\[
\mathcal R_{ND}(\xi,\delta),
\]
and pulls it all the way back through the actual continuum selected-branch placement map, so the same-charge dynamic verdict is expressed directly in the physical kernel ratios of the moving-throat branch.

---

## Purpose

Stage 230 solved the exact dynamic-sign problem on the full selected-branch classifier half-line:

- the wall-like dynamic window never became the first kill condition,
- the only live global kill still came from the transported static \(\Xi_1\) budget.

But Stage 230 still spoke in the abstract classifier coordinate
\[
\mathcal R_{ND}.
\]

The next honest step is therefore very specific:

> pull that classifier and its dynamic sign thresholds back through the **actual continuum selected-branch placement map**, so the same-charge verdict is stated directly in the physical kernel ratios of the moving-throat branch.

The main outputs of this stage are:

1. the exact physical classifier
   \[
   \mathcal R_{\rm phys}(\delta,R_{\rm target})
   :=
   \mathcal R_{ND}\!\bigl(\xi_{\rm req}(\delta,R_{\rm target}),\delta\bigr),
   \]
   where \(\xi_{\rm req}\) is the unique stable selected-branch point solving
   \[
   F(\xi,\delta)=R_{\rm target};
   \]
2. the exact monotonicity theorem
   \[
   \partial_{R_{\rm target}}\mathcal R_{\rm phys}<0,
   \]
   so larger normalization demand ratio always pushes the physical selected branch in the denominator-like / dynamically safer direction;
3. the exact pulled-back target thresholds
   \[
   R_{\rm flip}(\delta),
   \qquad
   R_{\rm den}(\delta),
   \]
   corresponding respectively to
   \[
   \mathcal R_{\rm phys}=\mathcal R_*,
   \qquad
   \mathcal R_{\rm phys}=1;
   \]
4. the equivalent exact inequalities on the continuum kernel ratios
   \[
   (\epsilon_\eta,\epsilon_W,\rho,Z_W,\delta_0,\Lambda)
   \]
   and on the mixed-baseline coordinate \(M_{\rm mix}\);
5. and the refined verdict that even on the **actual continuum-selected branch** the dynamic window is still not the first kill condition.

So after Stage 231, the same-charge corridor is still alive.
What changes is that the classifier map is no longer a sample-branch statement.
It is now a physical continuum-kernel statement.

---

## 1. Frozen input carried forward

### 1.1 Universal selected-branch geometry

From the carried universal D/N selected-branch geometry, the stable selected branch is controlled by the exact functions
\[
F(\xi,\delta)
=
\frac{(9\delta+11\xi)^4}{81(1-\xi)(9\delta^2+18\delta\xi+11\xi^2)^2},
\]
\[
G(\xi,\delta)
=
\frac{9\xi(\xi+\delta)}{9\delta+11\xi},
\]
with
\[
0\le \xi<1,\qquad \delta>0.
\]

Their exact monotonicities are
\[
\partial_\xi F
=
\frac{(9\delta+11\xi)^3\bigl(81\delta^3+189\delta^2\xi+72\delta^2+297\delta\xi^2+121\xi^3\bigr)}
{81(1-\xi)^2(9\delta^2+18\delta\xi+11\xi^2)^3}
>0,
\]
\[
\partial_\xi G
=
\frac{9(9\delta^2+18\delta\xi+11\xi^2)}{(9\delta+11\xi)^2}
>0.
\]

The endpoint data are
\[
F(0,\delta)=1,
\qquad
F\to+\infty\quad(\xi\to1^-),
\]
\[
G(0,\delta)=0,
\qquad
G_{\max}(\delta)=\frac{9(1+\delta)}{9\delta+11}.
\]

So for fixed \(\delta\), the physical selected branch is placed by a unique normalization locus
\[
F(\xi,\delta)=R_{\rm target},
\]
together with the support-feasibility frontier
\[
M_{\rm mix}\le G(\xi,\delta).
\]

### 1.2 Exact selected-branch classifier and Stage-230 dynamic thresholds

From Stages 229–230, the exact selected-branch numerator/denominator classifier is
\[
\mathcal R_{ND}(\xi,\delta)
=
\frac{72\delta^2(1-\xi)}
{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)},
\]
with strict monotonicity
\[
\partial_\xi \mathcal R_{ND}
=
-
\frac{72\delta^2\bigl(81\delta^3+261\delta^2+297\delta\xi(2-\xi)+\xi^2(363-242\xi)\bigr)}
{(9\delta+11\xi)^2(9\delta^2+18\delta\xi+11\xi^2)^2}
<0.
\]

The carried Stage-230 sign threshold is
\[
\mathcal R_*\approx 1.229255438463336,
\]
and the denominator-like threshold is
\[
\mathcal R_{ND}=1.
\]

The associated onset thresholds are
\[
\delta_*^{(\rm dyn)}=\frac{8}{9\mathcal R_*}\approx 0.723111617875019,
\qquad
\delta_{\rm den}=\frac89.
\]

So Stage 230 already gave two exact global statements:

- if \(\delta\ge \delta_*^{(\rm dyn)}\), the nonempty dynamic ceiling is infinite on the whole selected branch;
- if \(\delta\ge 8/9\), the whole selected branch is denominator-like.

### 1.3 Continuum placement map

From the carried continuum kernel extraction, the actual moving-throat branch is placed by the exact dimensionless ratios
\[
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{\rm eff}},
\qquad
\epsilon_W=\frac{c_{UW}^2\sigma}{K_UK_W^{\rm eff}},
\]
\[
\rho=\frac{c_{UW}c_{\eta U}}{K_Uc_{\eta W}},
\qquad
Z_W=\frac{c_{\eta W}^2}{K_\eta^{\rm eff}K_W^{\rm eff}},
\qquad
\delta_0=\frac{\pi^2T_w}{L^2K_\eta^{\rm eff}},
\]
together with the radiative demand scale
\[
\Lambda=\frac{27\pi^2Gc_s^5K_W^{\rm eff}}{20a^5c^5\mu_W}.
\]

The exact placement formulas are
\[
\delta=\frac{\delta_0}{1-\epsilon_\eta},
\]
\[
M_{\rm mix}
=
\frac{8Z_W(1+\rho)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon_W)},
\]
\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon_W)^2}{Z_W(1+\rho)^2},
\]
with exact product law
\[
R_{\rm target}M_{\rm mix}
=
\frac{8\Lambda(1-\epsilon_W)}{\pi^2}.
\]

So the actual physical selected branch is obtained by

1. computing the continuum ratios,
2. mapping them to \((\delta,M_{\rm mix},R_{\rm target})\),
3. solving the unique equation \(F(\xi,\delta)=R_{\rm target}\),
4. and checking \(M_{\rm mix}\le G(\xi,\delta)\).

---

## 2. Exact pullback of the classifier to the physical placement locus

Define the physical selected-branch point \(\xi_{\rm req}(\delta,R_{\rm target})\) by the unique stable solution of
\[
F(\xi_{\rm req},\delta)=R_{\rm target},
\qquad
R_{\rm target}\ge 1.
\]

Then the actual physical classifier is
\[
\boxed{
\mathcal R_{\rm phys}(\delta,R_{\rm target})
:=
\mathcal R_{ND}\!\bigl(\xi_{\rm req}(\delta,R_{\rm target}),\delta\bigr).
}
\]

This is the exact continuum pullback of the Stage-230 classifier.

Because
\[
\partial_\xi F>0,
\qquad
\partial_\xi \mathcal R_{ND}<0,
\]
implicit differentiation gives
\[
\frac{\partial \xi_{\rm req}}{\partial R_{\rm target}}
=
\frac{1}{\partial_\xi F}>0,
\]
and therefore
\[
\boxed{
\frac{\partial \mathcal R_{\rm phys}}{\partial R_{\rm target}}
=
\frac{\partial_\xi \mathcal R_{ND}}{\partial_\xi F}
<
0.
}
\]

So for fixed \(\delta\):

> larger normalization demand ratio \(R_{\rm target}\) always pushes the physical selected branch in the denominator-like / dynamically safer direction.

Using the exact product law
\[
R_{\rm target}M_{\rm mix}
=
\frac{8\Lambda(1-\epsilon_W)}{\pi^2},
\]
the same statement can be rewritten as
\[
\boxed{
\frac{\partial \mathcal R_{\rm phys}}{\partial M_{\rm mix}}>0
\quad
\text{(at fixed }\delta,\Lambda,\epsilon_W\text{)}.
}
\]

So larger mixed baseline drives the physical selected branch in the numerator-like direction, while smaller mixed baseline drives it in the denominator-like direction.

This is the first real physical reading of the Stage-230 classifier map.

---

## 3. Exact threshold compiler for any classifier cap

Fix any classifier cap \(c>0\).
Define the threshold polynomial
\[
P_c(\xi,\delta)
:=
c(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)-72\delta^2(1-\xi).
\]

Its exact derivative is
\[
\partial_\xi P_c
=
3\bigl(87c\delta^2+198c\delta\xi+121c\xi^2+24\delta^2\bigr)>0
\]
for
\[
\xi\ge 0,\qquad \delta>0,\qquad c>0.
\]

So if a threshold root exists, it is unique.

At onset,
\[
P_c(0,\delta)=9\delta^2(9c\delta-8),
\]
which is equivalent to the carried onset classifier
\[
\mathcal R_{ND}(0,\delta)=\frac{8}{9\delta}.
\]

Therefore the exact onset threshold for the classifier cap \(c\) is
\[
\boxed{
\delta_c=\frac{8}{9c}.
}
\]

This gives the exact pullback theorem:

- if \(\delta\ge \delta_c\), then the physical selected branch already satisfies \(\mathcal R_{\rm phys}\le c\) at onset, so the pulled-back target threshold is simply
  \[
  R_{\rm target}\ge 1;
  \]
- if \(0<\delta<\delta_c\), there is a unique \(\xi_c(\delta)\in(0,1)\) with
  \[
  \mathcal R_{ND}(\xi_c,\delta)=c,
  \]
  and the pulled-back target threshold is
  \[
  \boxed{
  R_c(\delta)=F(\xi_c(\delta),\delta)>1.
  }
  \]

So every dynamic-class condition on the abstract selected branch has a unique continuum-placement target curve.

---

## 4. Two pulled-back dynamic-class surfaces

### 4.1 Lower-pole sign-flip surface

Take
\[
c=\mathcal R_*\approx 1.229255438463336.
\]
Then
\[
\delta_c=\delta_*^{(\rm dyn)}=\frac{8}{9\mathcal R_*}\approx 0.723111617875019.
\]

Define the exact pulled-back sign-flip target curve
\[
\boxed{
R_{\rm flip}(\delta):=R_{\mathcal R_*}(\delta).
}
\]

Then the physical selected branch has
\[
\mathcal R_{\rm phys}\le \mathcal R_*
\iff
R_{\rm target}\ge R_{\rm flip}(\delta).
\]

Equivalently:

> once the physical normalization demand ratio exceeds \(R_{\rm flip}(\delta)\), the lower wall-like pole improves and the nonempty dynamic ceiling becomes infinite.

For
\[
\delta\ge \delta_*^{(\rm dyn)},
\]
the threshold collapses to onset:
\[
\boxed{
R_{\rm flip}(\delta)=1.
}
\]

### 4.2 Denominator-like surface

Take
\[
c=1.
\]
Then
\[
\delta_c=\frac89.
\]

Define the exact pulled-back denominator target curve
\[
\boxed{
R_{\rm den}(\delta):=R_1(\delta).
}
\]

Then
\[
\mathcal R_{\rm phys}\le 1
\iff
R_{\rm target}\ge R_{\rm den}(\delta).
\]

So once the physical target ratio exceeds \(R_{\rm den}(\delta)\), the actual selected branch is fully denominator-like.

For
\[
\delta\ge \frac89,
\]
this again collapses to onset:
\[
\boxed{
R_{\rm den}(\delta)=1.
}
\]

Because
\[
\mathcal R_*>1,
\]
the sign-flip surface is weaker than the denominator surface:
\[
R_{\rm flip}(\delta)\le R_{\rm den}(\delta).
\]

So as \(R_{\rm target}\) increases, the physical branch crosses the “lower pole improves” threshold before it reaches the fully denominator-like regime.

---

## 5. Sample pulled-back thresholds

The exact pulled-back thresholds on a few representative \(\delta\)-slices are:

| \(\delta\) | \(\xi_{\rm flip}\) | \(R_{\rm flip}\) | \(\xi_{\rm den}\) | \(R_{\rm den}\) |
|---:|---:|---:|---:|---:|
| \(0.25\) | \(0.087442106\) | \(1.330868539\) | \(0.107223051\) | \(1.393832566\) |
| \(0.50\) | \(0.051428579\) | \(1.139956630\) | \(0.081847938\) | \(1.221087062\) |
| \(0.75\) | \(0\) | \(1\) | \(0.032505121\) | \(1.071471867\) |

So the physical selected branch becomes dynamically sign-safe quite early.
For \(\delta=0.75\), the nonempty dynamic ceiling is already infinite from onset, even though the branch is not yet denominator-like from onset because \(0.75<8/9\).

This is a strong continuation of the Stage-230 picture.

---

## 6. Exact continuum-kernel inequalities

Because
\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon_W)^2}{Z_W(1+\rho)^2},
\]
the pulled-back sign-flip condition
\[
R_{\rm target}\ge R_{\rm flip}(\delta)
\]
is exactly equivalent to
\[
\boxed{
Z_W(1+\rho)^2
\le
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon_W)^2}{R_{\rm flip}(\delta)}.
}
\]

Likewise, the denominator-like condition
\[
R_{\rm target}\ge R_{\rm den}(\delta)
\]
is equivalent to
\[
\boxed{
Z_W(1+\rho)^2
\le
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon_W)^2}{R_{\rm den}(\delta)}.
}
\]

Using the exact product law
\[
R_{\rm target}M_{\rm mix}
=
\frac{8\Lambda(1-\epsilon_W)}{\pi^2},
\]
the same thresholds become upper bounds on the mixed baseline:
\[
\boxed{
M_{\rm mix}
\le
\frac{8\Lambda(1-\epsilon_W)}{\pi^2R_{\rm flip}(\delta)}
\quad\Longleftrightarrow\quad
\text{nonempty dynamic ceiling infinite},
}
\]
\[
\boxed{
M_{\rm mix}
\le
\frac{8\Lambda(1-\epsilon_W)}{\pi^2R_{\rm den}(\delta)}
\quad\Longleftrightarrow\quad
\text{physical branch denominator-like}.
}
\]

So at fixed product scale \(8\Lambda(1-\epsilon_W)/\pi^2\):

- lowering \(M_{\rm mix}\) first drives the actual branch across the dynamic sign-flip threshold,
- and lowering it further drives the branch fully into the denominator-like regime.

This is the first exact continuum-kernel dynamic-class map.

---

## 7. The static-first theorem survives the pullback

Stage 230 already proved the global inequalities
\[
\inf B_{\rm dyn}^{(\rm both)}
\approx 0.967282389363822
>
0.367930328492646
=
B_{\rm stat}^{(\rm both)},
\]
\[
\inf B_{\rm dyn}^{(\rm nonempty)}
\approx 0.990581810705233
>
0.737619063660757
=
B_{\rm stat}^{(\rm nonempty)}.
\]

The continuum placement map only restricts the full classifier half-line to a physical subset
\[
\mathcal R_{\rm phys}(\delta,R_{\rm target})\subseteq [0,\infty).
\]
So the same strict inequalities survive under pullback:
\[
\boxed{
B_{\rm dyn}^{(\rm both)}(\text{physical})
>
B_{\rm stat}^{(\rm both)},
}
\]
\[
\boxed{
B_{\rm dyn}^{(\rm nonempty)}(\text{physical})
>
B_{\rm stat}^{(\rm nonempty)}.
}
\]

Therefore:

> even on the actual continuum-selected branch, the first kill condition remains the transported static \(\Xi_1\) budget, not the wall-like dynamic window.

This is the central Stage-231 theorem.

---

## 8. Best current verdict after Stage 231

Stage 231 does not kill the same-charge corridor.
It sharpens the physical classification instead.

The exact selected-branch classifier is now pulled all the way back to the continuum kernel ratios:

1. for fixed \(\delta\), larger \(R_{\rm target}\) always makes the actual selected branch more denominator-like;
2. equivalently, at fixed product scale, larger \(M_{\rm mix}\) makes it more numerator-like;
3. there is an exact pulled-back sign-flip threshold \(R_{\rm flip}(\delta)\) beyond which the nonempty dynamic ceiling is infinite;
4. there is a second exact pulled-back denominator threshold \(R_{\rm den}(\delta)\) beyond which the physical branch is denominator-like;
5. these thresholds translate directly into exact inequalities on
   \[
   (\epsilon_\eta,\epsilon_W,\rho,Z_W,\delta_0,\Lambda)
   \]
   or on \(M_{\rm mix}\) through the exact product law;
6. but nowhere on the continuum placement map does the dynamic window become the first kill condition.

So after Stage 231, the same-charge idea is still alive.
The remaining first kill condition is still the **static placement of \(\Xi_1\) on the actual moving-throat branch**, not the wall-like dynamic window.

---

## 9. SymPy-backed status

The accompanying audit script verifies:

- the exact monotonicity formulas for \(F(\xi,\delta)\), \(G(\xi,\delta)\), and \(\mathcal R_{ND}(\xi,\delta)\);
- the implicit-derivative compiler
  \[
  \partial_{R_{\rm target}}\mathcal R_{\rm phys}
  =
  \frac{\partial_\xi\mathcal R_{ND}}{\partial_\xi F}
  <0;
  \]
- the exact threshold polynomial
  \[
  P_c(\xi,\delta)
  \]
  and its unique-root / onset law
  \[
  \delta_c=\frac{8}{9c};
  \]
- the pulled-back sign-flip and denominator thresholds \(R_{\rm flip}(\delta)\), \(R_{\rm den}(\delta)\) on representative slices;
- the exact continuum product law
  \[
  R_{\rm target}M_{\rm mix}
  =
  \frac{8\Lambda(1-\epsilon_W)}{\pi^2};
  \]
- and the Stage-230 static-first inequalities after pullback.

The note and script are therefore aligned: the stage is a clean compiler from the abstract selected-branch classifier to the actual continuum placement coordinates of the moving-throat branch.
