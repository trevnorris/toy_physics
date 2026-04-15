# Moving-Throat PDE — Stage 209: Strict `5`PN Even-Gate Package, the Surviving Mixed Corridor, and the Pure-Transfer Subcorridor

## Status

**Exact within the explicit finite-throat one-port weak-axisymmetric logarithmic-slope closure built on the Stage-206 compatibility branch, the Stage-208 primitive compiler, and the imported strict `5`PN even-gate package.**

This stage does **not** solve the full moving-throat PDE.
It takes the explicit mixed-sector survivor from Stage 208 and carries it into the stricter imported `5`PN first-order package
\[
\Xi_{\rm load},\qquad K_1,\qquad H_{\rm even},
\]
then determines exactly whether the corridor dies or survives on the concrete noncanonical compatibility branch.

---

## Purpose

Stage 208 showed that the explicit finite-throat one-port branch still supports a nontrivial mixed-sector corridor after imposing the **first-order conservative compensation surface** that keeps the actual grouped response fixed on that branch.

But the later `5`PN notes impose a stricter first-order package than Stage 208 used. In that later language, the surviving weak-axisymmetric grouped problem is organized by three scalars:

\[
\Xi_{\rm load}=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0},
\]
\[
K_1=D_{21}+\frac{D_{01}}{9},
\]
\[
H_{\rm even}=D_{41}-\frac{2}{3}D_{21}-\frac{D_{01}}{27}.
\]

So the next honest step is not another loose mechanism scan. It is:

> carry the Stage-208 mixed-sector survivor into the stricter imported `5`PN even-gate package and see whether the corridor collapses or survives.

The main outputs are:

1. the exact bridge
   \[
   \Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0},
   \]
   so the same-charge scalar from Stage 208 is already the imported `5`PN load defect;
2. the exact comparison between the Stage-208 compensation surface and the stricter `5`PN even-gate package;
3. the explicit mixed-sector-only strict-even-gate solve on the concrete Stage-206 compatibility point;
4. the sharper **pure-transfer** subcorridor that survives once one also enforces the Stage-208 conservative-shape preservation on this noncanonical sample branch;
5. the transported same-charge ceiling budgets on both the strict mixed corridor and the pure-transfer subcorridor.

So after Stage 209, the question is no longer

> does some mixed anisotropy survive?

It is now

> does the actual moving-throat branch realize the surviving same-charge effect primarily as mixed-sector outgoing-transfer enhancement, with the conservative one-pole bundle frozen at first order?

---

## 1. Frozen input carried forward

### 1.1 Explicit one-port compatibility branch from Stages 206–208

Keep the same finite-throat branch:

- wall and brane-like internal coordinate on the lowest N/N zero mode,
- trapped support and mixed channel on the lowest D/N half-wave,
- exact overlap constant
  \[
  \kappa=\frac{2\sqrt2}{\pi}.
  \]

With primitive parameters
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
\left(\frac12,\frac3{10},\frac25,\frac14,1,\frac75,2,1\right),
\]
Stage 208 fixed the isotropic one-pole base data
\[
D_0\approx 24.2373099886223,
\qquad
D_2\approx -1.18562046858190,
\qquad
D_4\approx -0.173991572849491,
\]
\[
u_2\approx 0.0489171640391802,
\qquad
u_4\approx 0.00957155575054425,
\qquad
\frac{D_4}{D_0}\approx -0.00717866681290820,
\]
with exact one-pole identity
\[
u_4-4u_2^2=0.
\]
The same branch also carries
\[
P_{0,\rm compat}=\frac{N_0}{D_0}\approx 0.002069792318062885.
\]

### 1.2 Stage-208 first-order primitive compiler

Stage 208 already compiled the primitive weak-axisymmetric slopes
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W})
\]
into
\[
D_{01},\qquad D_{21},\qquad D_{41},\qquad N_{01},\qquad \Xi_1.
\]
On the mixed-only sector, the Stage-208 conservative compensation surface was
\[
D_{21}=-u_2D_{01},
\qquad
D_{41}=\frac{D_4}{D_0}D_{01},
\]
and the transported same-charge scalar was already
\[
\Xi_1=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]

### 1.3 Stage-207 transported same-charge ceilings

At this same compatibility point the robust carried budgets were
\[
|\epsilon\Xi_1|\lesssim 0.367930328492646
\]
for the stricter `10%`-loss “both wall-like poles survive” criterion, and
\[
|\epsilon\Xi_1|\lesssim 0.737619063660757
\]
for the stricter `10%`-loss “nonempty wall-like corridor” criterion.
The looser `30%` windows remain
\[
|\epsilon\Xi_1|\lesssim 2.94889585703134,
\qquad
|\epsilon\Xi_1|\lesssim 4.63505472371892.
\]

---

## 2. Exact bridge to the imported `5`PN load defect

The imported weak-axisymmetric `5`PN load defect is
\[
\Xi_{\rm load}:=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]
But Stage 208 already gave
\[
\Xi_1=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]
So coefficientwise,
\[
\boxed{\Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0}.}
\]

That matters because we are not introducing a new scalar here. We are simply re-reading the Stage-208 same-charge scalar in the stricter imported `5`PN language.

So the same-charge scalar already is the `5`PN loading defect.

---

## 3. Exact comparison: Stage-208 compensation is weaker than the strict `5`PN package

The Stage-208 first-order conservative compensation surface on an arbitrary one-pole base branch is
\[
D_{21}=-u_2D_{01},
\qquad
D_{41}=\frac{D_4}{D_0}D_{01}.
\]
Insert this into the stricter imported `5`PN even gates:
\[
K_1=D_{21}+\frac{D_{01}}{9},
\qquad
H_{\rm even}=D_{41}-\frac{2}{3}D_{21}-\frac{D_{01}}{27}.
\]
Then exactly,
\[
\boxed{K_1=\left(\frac19-u_2\right)D_{01},}
\]
\[
\boxed{H_{\rm even}=\left(\frac{D_4}{D_0}+\frac{2u_2}{3}-\frac1{27}\right)D_{01}.}
\]
Using
\[
\frac{D_4}{D_0}=u_2^2-u_4,
\]
and the one-pole identity
\[
u_4=4u_2^2,
\]
one also gets
\[
\boxed{H_{\rm even}=\left(-3u_2^2+\frac{2u_2}{3}-\frac1{27}\right)D_{01}.}
\]

So the Stage-208 conservative-shape preservation and the stricter imported `5`PN even gates are **not** the same condition. They agree only on the canonical branch for which the coefficients above vanish.

On the explicit Stage-206 compatibility branch the coefficients are numerically
\[
K_1\approx 0.0621939470719309\,D_{01},
\qquad
H_{\rm even}\approx -0.0116042611571584\,D_{01}.
\]
Both are nonzero. Therefore, on this noncanonical sample branch,
\[
\boxed{\text{Stage-208 compensation + strict `5`PN even gates} \iff D_{01}=0.}
\]

This is the first real sharpening of the corridor.

It means that once the stricter imported `5`PN gates are imposed on top of the Stage-208 compensation surface, the surviving same-charge corridor can no longer use first-order conservative static loading. It must pass through a branch with
\[
D_{01}=0.
\]

---

## 4. Mixed-sector-only strict even-gate corridor

Now restrict attention to the mixed/U primitive family
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}).
\]
The strict even-gate matrix on the concrete branch is
\[
\begin{pmatrix}
-0.255028994532 & -0.132167046465 & -0.067875763349 & 0.568483003085 & 0.300375205864 \\
-0.086801409924 & -0.010480267714 & -0.045759251298 & 0.510038362482 & 0.131455867026
\end{pmatrix}.
\]
It has rank `2`, hence nullity `3`.

So the mixed sector **survives** the strict imported `5`PN even-gate package as a three-dimensional corridor.

A convenient raw null basis is
\[
w_1\approx(-0.606454972136,
            0.656652628212,
            1,
            0,
            0),
\]
\[
w_2\approx(6.983614208603,
           -9.174307357027,
            0,
            1,
            0),
\]
\[
w_3\approx(1.616693986742,
           -0.846872492318,
            0,
            0,
            1).
\]
The corresponding same-charge scalars are
\[
\Xi_1(w_1)\approx 1.33691841376792,
\]
\[
\Xi_1(w_2)\approx -13.9944400566810,
\]
\[
\Xi_1(w_3)\approx -5.02163500066813.
\]
So the strict even-gate package does **not** kill the mixed-sector corridor. It only deforms it relative to the weaker Stage-208 compensation surface.

Let `\Pi_{\rm even}` denote the ambient-Euclidean orthogonal projector onto this nullspace. Then the induced operator norm of the same-charge functional on the strict corridor is
\[
\sigma_{\rm even}=\|\Pi_{\rm even}\,\Xi_1\|_2.
\]
Numerically,
\[
\boxed{\sigma_{\rm even}\approx 2.67386816837173.}
\]
That gives a canonical same-charge gain scale for unit microscopic mixed-sector drift amplitude.

---

## 5. The pure-transfer subcorridor

Now impose the full intersection of

1. Stage-208 conservative-shape preservation, and
2. the stricter imported `5`PN even-gate package.

On this noncanonical sample branch, Section 3 showed that this is equivalent to solving
\[
D_{01}=0,
\qquad
D_{21}=0,
\qquad
D_{41}=0.
\]
In the mixed-only sector, that is the `3 x 5` linear system built from
\[
\text{eq1}=D_{21}+u_2D_{01},
\qquad
\text{eq2}=D_{41}-\frac{D_4}{D_0}D_{01},
\qquad
D_{01}=0.
\]
The intersection matrix has rank `3`, hence nullity `2`.

So a **two-dimensional** mixed-sector corridor still survives even after imposing both the Stage-208 compensation surface and the imported strict `5`PN even gates.

A convenient raw basis is
\[
t_1\approx(-4.359222794718,
           3.107402039105,
           18.703510605854,
           1,
           0),
\]
\[
t_2\approx(1.909256655687,
          -1.163651238154,
          -0.482414494705,
           0,
           1).
\]
On these directions,
\[
D_{01}(t_i)=D_{21}(t_i)=D_{41}(t_i)=0,
\]
while the same-charge scalar remains nonzero:
\[
\Xi_1(t_1)\approx 11.0106276743889,
\qquad
\Xi_1(t_2)\approx -5.66658382170817.
\]
At the same time,
\[
N_{01}(t_1)\approx 0.552361328292489,
\qquad
N_{01}(t_2)\approx -0.284270966124842.
\]
So on this subcorridor the whole effect is carried purely by outgoing-transfer loading:
\[
\boxed{D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\Xi_{\rm load}=\frac{N_{01}}{N_0}.}
\]

This is the cleanest surviving mechanism so far.

It means the same-charge corridor can survive all gates reached up to this stage with the conservative one-pole bundle frozen at first order. The remaining effect is purely a mixed-sector transfer enhancement.

Let `\Pi_{\rm transfer}` denote the ambient-Euclidean orthogonal projector onto this two-dimensional subspace. Then the induced same-charge norm is
\[
\sigma_{\rm transfer}=\|\Pi_{\rm transfer}\,\Xi_1\|_2.
\]
Numerically,
\[
\boxed{\sigma_{\rm transfer}\approx 2.31561904386057.}
\]

---

## 6. Transported same-charge ceiling budgets on the strict corridors

Interpret the ambient microscopic mixed-sector drift amplitude as
\[
\|x_{\rm mixed}\|_2=t.
\]
If the operator norm of `\Xi_1` on a corridor is `\sigma`, then the transported Stage-207 ceiling law becomes
\[
|\epsilon|t \le \frac{\text{budget}}{\sigma}.
\]

### 6.1 Strict three-dimensional even-gate corridor

Using
\[
\sigma_{\rm even}\approx 2.67386816837173,
\]
the robust carried budgets become:

for the stricter `10%`-loss “both wall-like poles survive” test,
\[
|\epsilon|t \lesssim 0.137602269567650;
\]

for the stricter `10%`-loss “nonempty wall-like corridor” test,
\[
|\epsilon|t \lesssim 0.275862165676603.
\]

The looser `30%` budgets are
\[
|\epsilon|t \lesssim 1.10285760977778,
\qquad
|\epsilon|t \lesssim 1.73346419189450.
\]

### 6.2 Pure-transfer two-dimensional subcorridor

Using
\[
\sigma_{\rm transfer}\approx 2.31561904386057,
\]
the robust carried budgets become:

for the stricter `10%`-loss “both wall-like poles survive” test,
\[
|\epsilon|t \lesssim 0.158890698998242;
\]

for the stricter `10%`-loss “nonempty wall-like corridor” test,
\[
|\epsilon|t \lesssim 0.318540765855427.
\]

The looser `30%` budgets are
\[
|\epsilon|t \lesssim 1.27348056877049,
\qquad
|\epsilon|t \lesssim 2.00164821411704.
\]

So the corridor narrows again, but it does not die.

And in fact the pure-transfer subcorridor leaves slightly **more** same-charge headroom than the larger strict even-gate corridor, because the surviving `\Xi_1` functional is a little more concentrated there.

---

## 7. What Stage 209 changes

Before this stage, the strongest positive statement was only

> a mixed-sector corridor survives the first-order conservative compensation surface from Stage 208.

After this stage, the picture is much sharper.

1. The imported strict `5`PN loading defect is **exactly** the same scalar already isolated in Stage 208:
   \[
   \Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0}.
   \]
2. The Stage-208 compensation surface is weaker than the stricter imported `5`PN even-gate package.
3. On the concrete noncanonical sample branch, imposing both structures together forces
   \[
   D_{01}=0.
   \]
4. Even after that sharpening, a nontrivial mixed-sector corridor still survives.
5. The sharpest surviving subcorridor is the **pure-transfer** family with
   \[
   D_{01}=D_{21}=D_{41}=0,
   \qquad
   \Xi_1=\frac{N_{01}}{N_0}.
   \]

So the best current summary is:

> the idea survives Stage 209, but no longer as generic mixed anisotropy. The sharpest surviving mechanism is mixed-sector outgoing-transfer enhancement with the conservative one-pole bundle frozen at first order.

That is a real narrowing, and it is the kind of narrowing we want.

---

## 8. Script-backed status

The accompanying SymPy audit verifies:

- the exact bridge
  \[
  \Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0};
  \]
- the exact formulas obtained by inserting the Stage-208 compensation surface into the imported strict `5`PN even gates;
- the concrete noncanonical compatibility-point coefficients in front of `D_{01}`;
- the strict mixed-only even-gate matrix, its rank-`2` nullity-`3` solve, the displayed raw null basis, the corresponding `\Xi_1` values, and the induced corridor norm `\sigma_{\rm even}`;
- the pure-transfer `3 x 5` intersection system, its rank-`3` nullity-`2` solve, the displayed raw basis, the corresponding `\Xi_1` and `N_{01}` values, and the induced norm `\sigma_{\rm transfer}`;
- and the transported same-charge ceiling budgets on both surviving strict corridors.

Supporting file:
- `moving_throat_pde_stage209_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.py`
