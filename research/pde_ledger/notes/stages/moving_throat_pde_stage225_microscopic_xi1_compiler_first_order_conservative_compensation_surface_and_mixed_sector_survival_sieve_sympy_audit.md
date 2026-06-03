# Moving-Throat PDE — Stage 225: Microscopic `\Xi_1` Compiler, First-Order Conservative Compensation Surface, and the Mixed-Sector Survival Sieve

## Status

**Exact within the explicit finite-throat one-port weak-axisymmetric logarithmic-slope closure built on the Stage-223 compatibility branch and the transported Stage-224 same-charge ceiling test.**

This stage does **not** solve the full moving-throat PDE.
It takes the first actual finite-throat base branch that already satisfies the isotropic `5`PN compatibility surface, perturbs it by weak-axisymmetric primitive microscopic drifts, and determines exactly which first-order mechanism families can preserve the conservative grouped response while still carrying a nonzero same-charge scalar
\[
\Xi_1=\frac{P_1}{P_0}.
\]

---

## Purpose

Stage 224 converted the actual-branch same-charge test into one explicit transported inequality in
\[
\Delta_{\rm norm},
\qquad
\Xi_1=\frac{P_1}{P_0}.
\]

So the next honest step is no longer another abstract branch-packet manipulation.
It is to compute `\Xi_1` **microscopically** on the explicit finite-throat one-port branch, while keeping track of the conservative grouped-`P2` conditions that the `5`PN branch still has to respect.

The main outputs are:

1. the exact arbitrary-base first-order formulas for
   \[
   u_2^{(1)},\qquad u_4^{(1)},\qquad \Xi_1,
   \]
   together with the exact compensation surface that preserves the conservative grouped response on a one-pole base branch;
2. the specialization of those formulas to the explicit Stage-223 compatibility point;
3. the exact primitive compiler from microscopic logarithmic slopes
   \[
   (x_K,x_M,x_{\lambda_B},x_{\varpi},x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W})
   \]
   to
   \[
   D_{01},\qquad D_{21},\qquad D_{41},\qquad N_{01},\qquad \Xi_1;
   \]
4. the first mechanism sieve:
   - wall-only,
   - pure BdG-only,
   - mixed-sector-only;
5. one explicit mixed-sector compensated family that survives and carries nonzero `\Xi_1`;
6. the direct translation of that surviving family into the transported Stage-224 same-charge headroom budgets.

So after Stage 225, the question is no longer

> can some microscopic anisotropy produce a useful `\Xi_1`?

It is now

> which microscopic families survive the conservative first-order compensation surface, and what same-charge headroom do they actually leave?

---

## 1. Frozen input carried forward

### 1.1 Explicit isotropic one-port base branch from Stage 223

Keep the same explicit finite-throat branch used in Stages 222–224:

- lowest N/N zero mode for the wall and the brane-like internal coordinate,
- lowest D/N half-wave for the trapped support and mixed coordinate,
- exact overlap constant
  \[
  \kappa = \frac{2\sqrt2}{\pi}.
  \]

With primitive couplings
\[
C=\kappa\lambda_B,
\qquad
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R,
\]
the static bundle data are
\[
\Delta = \Omega_U^2\Omega_W^2-R^2,
\qquad
S_2=\Omega_U^2+\Omega_W^2,
\qquad
H=G_U^2+G_W^2,
\]
\[
Q = G_U^2\Omega_W^2 + 2G_UG_WR + G_W^2\Omega_U^2,
\qquad
P = \Omega_U^2G_W + RG_U.
\]
The primitive moments are
\[
B_0=\frac{C^2}{\varpi^2},
\qquad
B_2=\frac{C^2}{\varpi^4},
\qquad
B_4=\frac{C^2}{\varpi^6},
\]
\[
Z_0=\frac{Q}{\Delta},
\qquad
Z_2=\frac{QS_2-H\Delta}{\Delta^2},
\qquad
Z_4=\frac{Q(S_2^2-\Delta)-S_2H\Delta}{\Delta^3},
\]
\[
N_0=\frac{P^2}{\Delta^2}.
\]

For the concrete Stage-223 compatibility sample
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
\left(\frac12,\frac3{10},\frac25,\frac14,1,\frac75,2,1\right),
\]
the exact compatibility wall stiffness is
\[
K_{\mathrm{compat}}\approx 24.4737548792910.
\]
The corresponding isotropic one-pole base values are
\[
D_0 \approx 24.2373099886223,
\qquad
D_2 \approx -1.18562046858190,
\qquad
D_4 \approx -0.173991572849491,
\]
\[
u_2 \approx 0.0489171640391802,
\qquad
u_4 \approx 0.00957155575054425,
\qquad
\frac{D_4}{D_0}\approx -0.00717866681290820,
\]
with exact one-pole identity
\[
u_4-4u_2^2=0.
\]
The same branch reproduces the carried compatibility prefactor
\[
P_{0,\mathrm{target,compat}}
=
\frac{N_0}{D_0}
\approx 0.00206979231806289.
\]

### 1.2 Stage-224 transported same-charge ceilings at that point

At the same compatibility point, the stricter weak-axisymmetric budgets carried from Stage 224 are
\[
|\epsilon\Xi_1| \lesssim 0.367930328492646
\]
for the `10%`-loss “both wall-like poles survive” criterion, and
\[
|\epsilon\Xi_1| \lesssim 0.737619063660757
\]
for the `10%`-loss “nonempty wall-like corridor” criterion.

The corresponding looser `30%` transported budgets are
\[
|\epsilon\Xi_1| \lesssim 2.94889585703134,
\qquad
|\epsilon\Xi_1| \lesssim 4.63505472371892.
\]

---

## 2. Exact arbitrary-base first-order formulas

Let a one-pole isotropic base branch be described by
\[
D_0,\qquad D_2,\qquad D_4,\qquad N_0.
\]
Then
\[
u_2=-\frac{D_2}{D_0},
\qquad
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\qquad
P_0=\frac{N_0}{D_0}.
\]
Introduce weak-axisymmetric first-order slopes
\[
D_{A0}=D_0+\epsilon\lambda_A D_{01},
\qquad
D_{A2}=D_2+\epsilon\lambda_A D_{21},
\qquad
D_{A4}=D_4+\epsilon\lambda_A D_{41},
\]
\[
N_{A0}=N_0+\epsilon\lambda_A N_{01},
\qquad
\lambda_{20}=1,
\quad
\lambda_{21}=\frac12,
\quad
\lambda_{22}=-1.
\]
Then the exact first-order physical slopes are
\[
\boxed{
u_2^{(1)}=
\frac{-D_0D_{21}+D_2D_{01}}{D_0^2}
= -\frac{D_{21}+u_2 D_{01}}{D_0},
}
\]
\[
\boxed{
u_4^{(1)}=
\frac{-D_0(D_0D_{41}+D_{01}D_4-2D_2D_{21})+2D_{01}(D_0D_4-D_2^2)}{D_0^3},
}
\]
\[
\boxed{
\Xi_1=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
\]

So on any one-pole base branch the same-charge scalar is already a **static loading mismatch** between

- outgoing-transfer strengthening `N_{01}/N_0`, and
- conservative static loading `D_{01}/D_0`.

---

## 3. Exact conservative first-order compensation surface

If the conservative grouped response is to remain fixed to first order,
\[
u_2^{(1)}=0,
\qquad
u_4^{(1)}=0,
\]
then the exact compensation surface is
\[
\boxed{D_{21}=-u_2 D_{01},}
\]
\[
\boxed{D_{41}=\frac{D_4}{D_0}D_{01}.}
\]
Using
\[
\frac{D_4}{D_0}=u_2^2-u_4,
\]
this can also be written as
\[
D_{41}=(u_2^2-u_4)D_{01}.
\]
On a one-pole branch,
\[
u_4=4u_2^2,
\]
so the second equation reduces to
\[
\boxed{D_{41}=-3u_2^2 D_{01}.}
\]

This is the exact arbitrary-base continuation of the canonical `5`PN even-preserving surface.
Once it is imposed, the only remaining first-order outlet is `\Xi_1`.

For the concrete Stage-223 compatibility point, the compensation surface becomes
\[
D_{21}\approx -0.0489171640391802\,D_{01},
\qquad
D_{41}\approx -0.00717866681290820\,D_{01}.
\]

---

## 4. Primitive logarithmic-slope compiler

Parameterize primitive weak-axisymmetric microscopic drifts by logarithmic slopes
\[
(x_K,x_M,x_{\lambda_B},x_{\varpi},x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}),
\]
so that each positive primitive parameter is dressed as
\[
p_A = p\,e^{\epsilon\lambda_A x_p}.
\]

Then the exact first-order primitive moment drifts are:

### 4.1 BdG sector
\[
B_{0,1}=B_0(2x_{\lambda_B}-2x_{\varpi}),
\]
\[
B_{2,1}=B_2(2x_{\lambda_B}-4x_{\varpi}),
\]
\[
B_{4,1}=B_4(2x_{\lambda_B}-6x_{\varpi}).
\]

### 4.2 Conservative Maxwell/mixed sector
\[
\Delta_1 = 2\Omega_U^2\Omega_W^2(x_{\Omega_U}+x_{\Omega_W})-2R^2x_{\lambda_R},
\]
\[
S_{2,1}=2\Omega_U^2x_{\Omega_U}+2\Omega_W^2x_{\Omega_W},
\]
\[
H_1=2G_U^2x_{\lambda_U}+2G_W^2x_{\lambda_W},
\]
\[
Q_1
=
2G_U^2\Omega_W^2(x_{\lambda_U}+x_{\Omega_W})
+
2G_UG_WR(x_{\lambda_U}+x_{\lambda_W}+x_{\lambda_R})
+
2G_W^2\Omega_U^2(x_{\lambda_W}+x_{\Omega_U}),
\]
\[
P_1^{\rm raw}
=
\Omega_U^2G_W(2x_{\Omega_U}+x_{\lambda_W})
+
RG_U(x_{\lambda_R}+x_{\lambda_U}).
\]
So
\[
Z_{0,1}=\frac{Q_1\Delta-Q\Delta_1}{\Delta^2},
\]
\[
Z_{2,1}=
\frac{\Delta(-\Delta H_1-H\Delta_1+QS_{2,1}+S_2Q_1)+2\Delta_1(\Delta H-QS_2)}{\Delta^3},
\]
\[
Z_{4,1}=
-\frac{\Delta^2HS_{2,1}+\Delta^2S_2H_1+\Delta^2Q_1-2\Delta HS_2\Delta_1-2\Delta QS_2S_{2,1}-2\Delta Q\Delta_1-\Delta S_2^2Q_1+3QS_2^2\Delta_1}{\Delta^4},
\]
\[
N_{0,1}=\frac{2PP_1^{\rm raw}}{\Delta^2}-\frac{2P^2\Delta_1}{\Delta^3}.
\]

### 4.3 First-order bundle compiler
\[
\boxed{D_{01}=Kx_K-B_{0,1}-Z_{0,1},}
\]
\[
\boxed{D_{21}=-(Mx_M+B_{2,1}+Z_{2,1}),}
\]
\[
\boxed{D_{41}=-(B_{4,1}+Z_{4,1}),}
\]
\[
\boxed{N_{01}=N_{0,1}.}
\]

### 4.4 Concrete microscopic `\Xi_1` compiler at the compatibility point

On the concrete Stage-223 compatibility point, the same-charge scalar compiles numerically to
\[
\boxed{
\begin{aligned}
\Xi_1 \approx{}&
-1.00975540977030\,x_K
+0.00418038073077834\,x_{\lambda_B}
-0.00418038073077834\,x_{\varpi} \\
&+0.324464020216766\,x_{\lambda_U}
+1.69086641859305\,x_{\lambda_W}
+0.423379354382463\,x_{\lambda_R} \\
&-0.747843374599229\,x_{\Omega_U}
-4.11424577297551\,x_{\Omega_W}.
\end{aligned}
}
\]
The `x_M` coefficient vanishes identically on this base branch.
So the strongest positive same-charge leverage sits in the mixed-channel loading `x_{\lambda_W}`, while the strongest negative leverage sits in the mixed frequency drift `x_{\Omega_W}`.

---

## 5. Mechanism sieve

### 5.1 Wall-only family — exact generic no-go

If only
\[
(x_K,x_M)
\]
are active, then
\[
D_{01}=Kx_K,
\qquad
D_{21}=-Mx_M,
\qquad
D_{41}=0.
\]
The conservative first-order compensation equations become
\[
Ku_2 x_K - Mx_M = 0,
\qquad
-\frac{D_4}{D_0}Kx_K = 0.
\]
So whenever
\[
\frac{D_4}{D_0}\neq 0,
\]
the second equation forces
\[
x_K=0,
\]
and then the first forces
\[
x_M=0.
\]
Therefore
\[
\boxed{\text{wall-only compensated deformations are generically trivial.}}
\]

This already kills the naive pure-wall route on the Stage-223 compatibility point, since there
\[
\frac{D_4}{D_0}\approx -0.00717866681290820\neq 0.
\]

### 5.2 Pure BdG family — exact sample-point no-go

If only
\[
(x_{\lambda_B},x_{\varpi})
\]
are active, then the compensation equations are
\[
\begin{pmatrix}
-(B_2+u_2B_0) & 2B_2+u_2B_0 \\
-(B_4-\tfrac{D_4}{D_0}B_0) & 3B_4-\tfrac{D_4}{D_0}B_0
\end{pmatrix}
\binom{x_{\lambda_B}}{x_{\varpi}}=0.
\]
Its exact determinant is
\[
\boxed{
\Delta_{\rm BdG}
=
-B_0B_2\frac{D_4}{D_0}-2B_0B_4u_2-B_2B_4.
}
\]
On the Stage-223 compatibility point,
\[
\Delta_{\rm BdG}\approx -5.11886996120011\times 10^{-5}\neq 0.
\]
So
\[
\boxed{\text{the pure BdG compensated family is also trivial on the concrete branch.}}
\]

This is the same structural conclusion reached later in the `5`PN continuation notes in a stricter language: neither pure wall nor pure support/BdG anisotropy carries the live weak-axisymmetric corridor by itself.

### 5.3 Mixed-sector-only family — explicit surviving corridor

Now activate only the mixed/U family
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}).
\]
On the Stage-223 compatibility point, the compensation matrix is
\[
\begin{pmatrix}
-0.241952861865934 & -0.122133861432532 & -0.0656784156312263 & 0.553209522700447 & 0.288144673113677 \\
-0.250543086743604 & -0.0937748521387244 & -0.0899548469020231 & 0.881694465041011 & 0.325834311088034
\end{pmatrix}.
\]
It has rank `2`, hence nullity `3`.

A convenient null basis is
\[
v_1\approx(-0.610255553634424,\ 0.671187016268095,\ 1,\ 0,\ 0),
\]
\[
v_2\approx(7.05469842496522,\ -9.44615143817664,\ 0,\ 1,\ 0),
\]
\[
v_3\approx(1.61486053113911,\ -0.839860892848583,\ 0,\ 0,\ 1).
\]
The corresponding same-charge slopes are
\[
\Xi_1(v_1)\approx 1.36026097049402,
\qquad
\Xi_1(v_2)\approx -14.4310278139755,
\qquad
\Xi_1(v_3)\approx -5.01037421295998.
\]
So the mixed/U family retains a **nontrivial compensated nullspace** and can still carry nonzero `\Xi_1`.

Therefore
\[
\boxed{\text{the same-charge idea survives this stage only in a constrained mixed-sector corridor.}}
\]

---

## 6. Direct same-charge headroom on the first surviving mixed family

Choose the first surviving mixed basis vector `v_1` and write its microscopic amplitude as `t`.
Then
\[
\Xi_1 = \sigma_1 t,
\qquad
\sigma_1\approx 1.36026097049402.
\]
The transported Stage-224 ceiling law is
\[
|\epsilon\Xi_1| \le \text{budget},
\]
so on this family it becomes
\[
|\epsilon t| \le \frac{\text{budget}}{\sigma_1}.
\]

The explicit headroom windows are:

### `10%`-loss, both wall-like poles survive
\[
|\epsilon t| \lesssim 0.270485102839510.
\]

### `10%`-loss, nonempty wall-like corridor
\[
|\epsilon t| \lesssim 0.542262903708006.
\]

### `30%`-loss, both wall-like poles survive
\[
|\epsilon t| \lesssim 2.16788978070904.
\]

### `30%`-loss, nonempty wall-like corridor
\[
|\epsilon t| \lesssim 3.40747461278373.
\]

So the strict transported window is not huge, but it is definitely nonzero.
That is exactly what we would expect if the corridor is real but narrow.

---

## 7. What Stage 225 changes

Before this stage, the same-charge continuation was still phrased as

> compute `\Xi_1` somehow from the branch and compare it against the Stage-224 ceiling.

After this stage, the problem is materially sharper.

1. The microscopic same-charge scalar is now compiled directly from primitive one-port weak-axisymmetric slopes.
2. The conservative first-order compensation surface is exact on any one-pole base branch:
   \[
   D_{21}=-u_2 D_{01},
   \qquad
   D_{41}=\frac{D_4}{D_0}D_{01}.
   \]
3. Pure wall and pure BdG deformations are killed once that conservative surface is imposed.
4. The mixed/U family survives and carries an explicit nonzero `\Xi_1`.
5. The transported same-charge window can now be read directly as a bound on one microscopic mixed-sector amplitude.

So the best current summary is:

> the idea survives this stage, but only as a constrained mixed-sector corridor.
> Pure wall-only or pure support/BdG anisotropy does not survive the first conservative same-charge compensation test.

That is exactly the narrowing we wanted.

---

## 8. SymPy-backed status

The paired audit script verifies all of the following:

1. the exact arbitrary-base formulas for
   \[
   u_2^{(1)},\qquad u_4^{(1)},\qquad \Xi_1;
   \]
2. the exact first-order conservative compensation surface;
3. the full primitive logarithmic-slope compiler;
4. the Stage-223 compatibility-point values
   \[
   D_0,\qquad D_2,\qquad D_4,\qquad u_2,\qquad u_4,\qquad P_{0,\mathrm{target,compat}};
   \]
5. the concrete compatibility-point `\Xi_1` linear form;
6. the wall-only and pure-BdG no-go results;
7. the mixed/U compensation matrix, its rank-`2` / nullity-`3` structure, the convenient null basis above, and the corresponding `\Xi_1` values;
8. the transported Stage-224 amplitude windows on the first surviving mixed family.

Supporting file:
- `moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py`
