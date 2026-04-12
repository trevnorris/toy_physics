# 5PN continuation notes — Stages 262–264

This batch takes the exact blocked Family-1 corridor from Stages 258–261 and applies it directly to the **actual selected moving-throat grouped-`P2` branch** rather than to a generic one-pole carrier.

The key input from the earlier 5PN notes is the actual grouped-`P2` carrier formula

\[
c_{\rm pole}=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\]

with the actual isotropic branch at

\[
\epsilon_2=\epsilon_4=0
\quad\Longrightarrow\quad
c_{\rm pole}=\frac14.
\]

So the continuation point after Stage 261 is no longer “what happens for a generic one-pole branch?” It is:

> what does the **actual** moving-throat grouped-`P2` branch do inside the exact blocked Family-1 corridor, and how robust is that branch under the first weak-anisotropy contamination?

---

## Stage 262 — exact actual grouped-`P2` branch map

Using

\[
c_{\rm pole}=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\qquad
\zeta_{\rm req}
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]

we get the exact blocked demand map in the **actual** grouped-`P2` variables:

\[
\rho_\alpha
=
\frac{4(1+\epsilon_2)^2}{4\epsilon_2^2+8\epsilon_2-\epsilon_4+3},
\]

\[
\zeta_{\rm req}
=
\frac{1+\epsilon_4}
{4(1+\epsilon_2)^2(1-\epsilon_{\rm blk})-(1-2\epsilon_{\rm blk})(1+\epsilon_4)}.
\]

At the actual isotropic point,

\[
\epsilon_2=\epsilon_4=0,
\qquad
c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac{1}{3-2\epsilon_{\rm blk}}.
\]

So the exact actual isotropic grouped-`P2` branch carries precisely the same demand point already isolated earlier in generic one-pole language.

### Exact twin-safety numerator

The universal lowest-twin boundary becomes the exact branch numerator

\[
M_{\rm twin}
:=
2(1+\epsilon_2)^2-(1+\epsilon_4).
\]

It satisfies

\[
c_{\rm pole}-\frac12
=
-\frac{M_{\rm twin}}{4(1+\epsilon_2)^2},
\]

and

\[
\zeta_{\rm req}-1
=
-\frac{2(1-\epsilon_{\rm blk})M_{\rm twin}}
{4(1+\epsilon_2)^2(1-\epsilon_{\rm blk})-(1-2\epsilon_{\rm blk})(1+\epsilon_4)}.
\]

So the actual grouped-`P2` branch is

- twin-safe iff `M_twin >= 0`,
- exactly on the twin boundary iff `M_twin = 0`,
- non-twin iff `M_twin < 0`.

### Exact contamination monotonicity

The geometry-lane contamination acts asymmetrically:

\[
\frac{\partial c_{\rm pole}}{\partial \epsilon_2}
=
-\frac{1+\epsilon_4}{2(1+\epsilon_2)^3}<0,
\qquad
\frac{\partial c_{\rm pole}}{\partial \epsilon_4}
=
\frac{1}{4(1+\epsilon_2)^2}>0,
\]

and on every admissible blocked branch,

\[
\frac{\partial \zeta_{\rm req}}{\partial \epsilon_2}<0,
\qquad
\frac{\partial \zeta_{\rm req}}{\partial \epsilon_4}>0.
\]

So positive `epsilon_2` contamination **softens** the support demand, while positive `epsilon_4` contamination **hardens** it.

That is already more informative than the generic one-pole picture.

---

## Stage 263 — exact blocked Family-1 corridor on the actual branch

The exact blocked Family-1 ceiling from Stage 261 is

\[
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})
=
\frac{\zeta_{\max}^{(F1)}(1-\epsilon_{\rm blk})}
{1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}}.
\]

Rewriting directly in the actual grouped-`P2` variables gives the second exact numerator

\[
M_{F1}
:=
4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})(1+\epsilon_2)^2-(1+\epsilon_4).
\]

Then

- `M_F1 > 0` is exact Family-1 admissibility,
- `M_twin >= 0` is exact lowest-twin sufficiency.

So the whole blocked corridor is now expressed directly in the actual branch variables.

### Actual isotropic margins

At the physical isotropic point,

\[
M_{\rm twin}(0,0)=1,
\qquad
M_{F1}(0,0)=4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})-1.
\]

So the actual isotropic branch sits exactly **one full unit** below the universal twin/non-twin boundary.

And because `zeta_max^(F1) > 1`,

\[
\frac{d}{d\epsilon_{\rm blk}}M_{F1}(0,0)
=
\frac{4\zeta_{\max}^{(F1)}(\zeta_{\max}^{(F1)}-1)}{\bigl(1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}\bigr)^2}>0,
\]

the exact Family-1 admissibility margin at the isotropic point **grows** with blocking.

On the exact `Lambda_EM`-refreshed branch we get numerically:

- twin margin at the isotropic point: `1`,
- Family-1 admissibility margin at `epsilon_blk = 0`:
  `4 c_pole,max^(F1)(0) - 1 ≈ 1.84644104209044`,
- Family-1 admissibility margin at the critical blocking ceiling:
  `3`.

So the actual isotropic grouped-`P2` branch is not merely admissible. It sits comfortably inside the exact blocked corridor.

---

## Stage 264 — exact second-order weak-anisotropy tolerance theorem

Stage 78 in the notes says the first nonzero geometry contamination enters only at

\[
\epsilon_2,\epsilon_4 = O(\chi^2).
\]

Write this as

\[
\epsilon_2 = a_2 \chi^2,
\qquad
\epsilon_4 = a_4 \chi^2,
\qquad
y := \chi^2.
\]

Then the exact twin-safe and Family-1 conditions become explicit quadratics in `y`.

### Exact twin-safe quadratic

\[
M_{\rm twin}(y)
=
2a_2^2 y^2 + (4a_2-a_4)y + 1.
\]

So the actual branch remains twin-safe iff `M_twin(y) >= 0`.

### Exact Family-1 admissibility quadratic

\[
M_{F1}(y)
=
4c_{\rm pole,max}^{(F1)}a_2^2 y^2
+
(8c_{\rm pole,max}^{(F1)}a_2-a_4)y
+
(4c_{\rm pole,max}^{(F1)}-1).
\]

So the exact blocked Family-1 corridor persists iff `M_F1(y) > 0`.

### Three immediate theorem-level consequences

1. **Finite safe neighborhood around the isotropic point**

   Because
   \[
   M_{\rm twin}(0)=1,
   \qquad
   M_{F1}(0)=4c_{\rm pole,max}^{(F1)}-1>0,
   \]
   every weak-anisotropy branch has a finite safe neighborhood around `chi = 0`.

2. **Initial drift controlled only by `a4 - 2 a2`**

   The exact initial slopes are
   \[
   \left.\frac{dc_{\rm pole}}{d(\chi^2)}\right|_{\chi=0}
   =
   \frac{a_4-2a_2}{4},
   \]
   \[
   \left.\frac{d\zeta_{\rm req}}{d(\chi^2)}\right|_{\chi=0}
   =
   \frac{4(1-\epsilon_{\rm blk})(a_4-2a_2)}{(3-2\epsilon_{\rm blk})^2}.
   \]

   So:
   - `a4 - 2 a2 > 0` pushes the actual branch toward larger pole weight and larger support demand,
   - `a4 - 2 a2 < 0` pushes it toward smaller pole weight and smaller support demand.

3. **Exit can occur only through exact quadratic roots**

   Any eventual loss of twin-safety or Family-1 admissibility must occur through the corresponding positive roots of `M_twin(y)` or `M_F1(y)`. It cannot happen arbitrarily close to the actual isotropic point.

So the Stage-78 `O(chi^2)` statement is now turned into an explicit algebraic tolerance theorem.

---

## Net result after Stages 262–264

The “next move” after the generic one-pole corridor is now complete.

1. The exact blocked support corridor has been rewritten directly in the **actual** selected moving-throat grouped-`P2` branch variables `epsilon_2, epsilon_4`.
2. The actual isotropic branch sits at
   \[
   c_{\rm pole}=\frac14,
   \qquad
   M_{\rm twin}=1,
   \qquad
   M_{F1}=4c_{\rm pole,max}^{(F1)}-1>0,
   \]
   so it is safely inside both the universal twin-safe strip and the exact blocked Family-1 corridor.
3. The first weak-anisotropy contamination has been promoted from the qualitative statement `O(chi^2)` to exact quadratics in `chi^2`.

So the next honest theorem gate is now very sharp:

> extract the actual weak-anisotropy coefficients `(a_2,a_4)` — equivalently the first concrete `l=0 <-> l=2` induced geometry-mixing law — from the moving-throat branch, then test those coefficients against the exact quadratic twin-safe and Family-1 admissibility conditions above.
