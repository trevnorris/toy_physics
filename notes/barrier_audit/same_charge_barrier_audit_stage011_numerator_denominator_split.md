# Same-Charge Barrier Audit — Stage 011: Numerator/Denominator Split of the Pure-Transfer Corridor and the First Actual Dynamic-Window Test

## 0. Purpose

Stage 010 showed that the strict same-charge survivor is no longer a generic mixed-sector anisotropy. On the concrete compatibility branch it has already collapsed to the **pure-transfer** subcorridor
\[
D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\frac{N_{01}}{N_0}=2\,\delta\ln\Lambda,
\qquad
\Lambda=\frac{P}{\Delta}.
\]

So the next honest question is even sharper:

> once the conservative one-pole bundle is frozen at first order, is the surviving same-charge effect being carried by the **load numerator** or by the **load denominator**, and does either option die once the actual wall-like dynamic window is imposed?

That is the job of this stage.

The main outputs are:

1. the exact static split
   \[
   \Xi_1 = 2(\pi_1-\delta_1),
   \qquad
   \pi_1:=\frac{P_{01}}{P},
   \qquad
   \delta_1:=\frac{\Delta_{01}}{\Delta};
   \]
2. the exact subcorridor theorem:
   - pure-transfer + \(\pi_1=0\) leaves a 1D **numerator-rigid** survivor,
   - pure-transfer + \(\delta_1=0\) leaves a 1D **denominator-rigid** survivor,
   - imposing both \(\pi_1=0\) and \(\delta_1=0\) kills the corridor exactly;
3. the explicit positive-\(\Xi_1\) unit directions on the concrete sample branch;
4. the first actual wall-like dynamic-window audit on those two 1D survivors;
5. and the sharp verdict that **neither** rigid split is killed by the dynamic window on the concrete compatibility point — the first ceiling is still the transported static one.

So after Stage 011 the live question is not whether the pure-transfer corridor can survive either rigidity split at all. It is which split the real PDE-selected mixed branch actually resembles.

---

## 1. Frozen input carried forward

### 1.1 Concrete compatibility branch

Keep the same explicit finite-throat one-port branch used in Stages 006–010:

- lowest N/N zero mode for the wall and brane-like internal coordinate,
- lowest D/N half-wave for the trapped support and mixed coordinate,
- exact overlap constant
  \[
  \kappa=\frac{2\sqrt2}{\pi},
  \]
- primitive parameters
  \[
  (\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
  =\left(\frac12,\frac3{10},\frac25,\frac14,1,\frac75,2,1\right),
  \]
- and the exact isotropic compatibility wall stiffness
  \[
  K_{\rm compat}\approx 24.473754879290965.
  \]

The static one-port primitives remain
\[
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R,
\]
\[
\Delta=\Omega_U^2\Omega_W^2-R^2,
\qquad
P=\Omega_U^2G_W+RG_U,
\qquad
N_0=\frac{P^2}{\Delta^2}.
\]

### 1.2 Stage-010 pure-transfer corridor

Stage 010 already reduced the strict same-charge survivor to
\[
D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\frac{N_{01}}{N_0}.
\]

So all that remains here is to split the static outgoing load factor
\[
\Lambda=\frac{P}{\Delta}
\]
into its numerator and denominator pieces.

---

## 2. Exact numerator/denominator split theorem

Define the exact static slopes
\[
\pi_1:=\frac{P_{01}}{P},
\qquad
\delta_1:=\frac{\Delta_{01}}{\Delta}.
\]
Then on the pure-transfer corridor,
\[
\Xi_1=\frac{N_{01}}{N_0}=2\frac{P_{01}}{P}-2\frac{\Delta_{01}}{\Delta},
\]
so
\[
\boxed{
\Xi_1 = 2(\pi_1-\delta_1).
}
\]

On the concrete sample branch, the exact row formulas in the mixed primitive slope space
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W})
\]
are
\[
\pi_1
=
\frac{3}{19}x_{\lambda_U}
+\frac{16}{19}x_{\lambda_W}
+\frac{3}{19}x_{\lambda_R}
+\frac{32}{19}x_{\Omega_U},
\]
\[
\delta_1
=
\frac{50}{25-98\pi^2}x_{\lambda_R}
+\frac{196\pi^2}{98\pi^2-25}x_{\Omega_U}
+\frac{196\pi^2}{98\pi^2-25}x_{\Omega_W}.
\]

So the numerator and denominator are not sampling the same microscopic slots:

- the **numerator** sees \(\lambda_U,\lambda_W,\lambda_R,\Omega_U\),
- the **denominator** sees \(\lambda_R,\Omega_U,\Omega_W\),
- and only \(\lambda_R\) and \((\Omega_U,\Omega_W)\) are shared.

That already tells us the two rigidities are physically distinct, not just algebraically complementary.

### 2.1 Exact subcorridor counts

Inside the five-dimensional mixed primitive slope space, the audit gives:

- pure-transfer:
  \[
  \operatorname{rank}=3,
  \qquad
  \operatorname{nullity}=2;
  \]
- pure-transfer + numerator rigidity \(\pi_1=0\):
  \[
  \operatorname{rank}=4,
  \qquad
  \operatorname{nullity}=1;
  \]
- pure-transfer + denominator rigidity \(\delta_1=0\):
  \[
  \operatorname{rank}=4,
  \qquad
  \operatorname{nullity}=1;
  \]
- pure-transfer + both rigidities:
  \[
  \operatorname{rank}=5,
  \qquad
  \operatorname{nullity}=0.
  \]

Equivalently, on a basis of the pure-transfer nullspace, the exact reduced determinant is
\[
\det[(\pi_1,\delta_1)|_{\rm pure\ transfer}]
=
\frac{196(200+147\pi^2)(80000+343225\pi^2+43218\pi^4)}{475(8670000+14894275\pi^2+2117682\pi^4)}
\neq 0.
\]

So the split theorem is exact:
\[
\boxed{
\text{pure-transfer} + \pi_1=0 \Rightarrow \text{1D survivor},
}
\]
\[
\boxed{
\text{pure-transfer} + \delta_1=0 \Rightarrow \text{1D survivor},
}
\]
\[
\boxed{
\text{pure-transfer} + \pi_1=0 + \delta_1=0 \Rightarrow \text{only the trivial drift}.
}
\]

This is the Stage-011 analogue of the Stage-010 co-loading no-go.

---

## 3. Positive-\(\Xi_1\) unit survivors on the concrete branch

Orient the surviving directions so that \(\Xi_1>0\), since that is the sign relevant to the same-charge barrier-softening corridor.

### 3.1 Numerator-rigid branch \(\pi_1=0\)

A Euclidean unit generator is
\[
 v_{\rm num}
 \approx
 (-0.55551149,
 \,0.31814576,
 \, -0.65766801,
 \, -0.04533730,
 \, -0.39447126),
\]
with
\[
\pi_1(v_{\rm num})=0,
\qquad
\delta_1(v_{\rm num})\approx -0.86805617,
\qquad
\Xi_1(v_{\rm num})\approx 1.73611235.
\]
So this branch carries the same-charge signal entirely through the denominator:
\[
\Xi_1=-2\delta_1.
\]

### 3.2 Denominator-rigid branch \(\delta_1=0\)

A Euclidean unit generator is
\[
 v_{\rm den}
 \approx
 (-0.26583993,
 \,0.18448137,
 \,0.94454459,
 \,0.04984499,
 \, -0.02543112),
\]
with
\[
\delta_1(v_{\rm den})=0,
\qquad
\pi_1(v_{\rm den})\approx 0.34646608,
\qquad
\Xi_1(v_{\rm den})\approx 0.69293215.
\]
So this branch carries the same-charge signal entirely through the numerator:
\[
\Xi_1=2\pi_1.
\]

### 3.3 Immediate static reading

The numerator-rigid survivor produces a larger same-charge scalar per unit mixed drift:
\[
\sigma_{\rm num}\approx 1.73611235,
\qquad
\sigma_{\rm den}\approx 0.69293215.
\]
So at fixed ambient microscopic amplitude it is the **stronger static lever**. The denominator-rigid branch is the **gentler** one.

At this point alone one might be tempted to prefer the numerator-rigid branch. Stage 011 shows why the dynamic test is needed before making that judgment.

---

## 4. First actual dynamic-window split on the wall-like poles

Now carry those two unit directions into the actual pole census of the concrete compatibility branch.

At the undeformed compatibility point, the wall-like poles are
\[
\omega_-\approx 1.997535678933614,
\qquad
\mathcal R_{Q,-}\approx 30.199907560250075,
\]
\[
\omega_+\approx 4.949054323643126,
\qquad
\mathcal R_{Q,+}\approx 36.171186483269487,
\]
with
\[
P_0\approx 0.002069792318062885.
\]

Using symmetric finite difference on the full pole census gives the first-order log-slopes.

### 4.1 Numerator-rigid positive-\(\Xi_1\) motion

The audit finds
\[
\delta\ln P_0 \approx +1.73611235,
\]
\[
\delta\ln \mathcal R_{Q,+} \approx -0.52346582,
\qquad
\delta\ln \mathcal R_{Q,-} \approx +0.71358484,
\]
with only negligible wall-pole frequency drift.

So the numerator-rigid branch has a very specific dynamic signature:

- it **hurts** the upper wall-like dynamic figure,
- but it **improves** the lower wall-like dynamic figure.

This is a split-sign dynamic response.

### 4.2 Denominator-rigid positive-\(\Xi_1\) motion

The audit finds
\[
\delta\ln P_0 \approx +0.69293215,
\]
\[
\delta\ln \mathcal R_{Q,+} \approx -0.35245541,
\qquad
\delta\ln \mathcal R_{Q,-} \approx -0.23169484,
\]
again with negligible wall-pole frequency drift.

So the denominator-rigid branch has the opposite qualitative pattern:

- it **hurts both** wall-like dynamic figures,
- but it does so more mildly than the numerator-rigid branch hurts the upper wall pole.

This is a same-sign dynamic penalty.

That is the first genuinely physical difference between the two rigid subcorridors.

---

## 5. Comparison with the actual dynamic window

Use the same local wall-like survival threshold carried from the earlier stages.
At the stricter `10%`-loss benchmark,
\[
\mathcal R_{Q,\rm req}\approx 21.8545662963584.
\]

### 5.1 Dynamic ceilings

The first-order wall-like dynamic ceilings are:

#### Numerator-rigid
\[
|\epsilon|t \lesssim 0.96253269
\qquad
(\text{both wall-like poles survive}),
\]
\[
|\epsilon|t \lesssim \infty
\qquad
(\text{nonempty wall-like corridor}).
\]
The nonempty dynamic ceiling is infinite at first order because one wall-like pole improves while the other worsens.

#### Denominator-rigid
\[
|\epsilon|t \lesssim 1.39592653
\qquad
(\text{both wall-like poles survive}),
\]
\[
|\epsilon|t \lesssim 1.42955095
\qquad
(\text{nonempty wall-like corridor}).
\]
So the denominator-rigid branch is the **only** one with a genuinely finite nonempty dynamic ceiling on the concrete branch.

### 5.2 Static ceilings from the carried Stage-007 transport

But the transported static ceilings are still much tighter:

#### Numerator-rigid
\[
|\epsilon|t \lesssim 0.21192772
\qquad
(\text{both wall-like poles}),
\]
\[
|\epsilon|t \lesssim 0.42486828
\qquad
(\text{nonempty wall-like corridor}).
\]

#### Denominator-rigid
\[
|\epsilon|t \lesssim 0.53097598
\qquad
(\text{both wall-like poles}),
\]
\[
|\epsilon|t \lesssim 1.06448959
\qquad
(\text{nonempty wall-like corridor}).
\]

So on the actual sample compatibility point,
\[
\boxed{
\text{dynamic ceiling} > \text{transported static ceiling}
}
\]
for both rigid splits.

That is the decisive Stage-011 result.

---

## 6. What Stage 011 changes

Before this stage, the best statement was only:

> the Stage-010 pure-transfer corridor survives unless both the interference and hybridization pieces are rigid simultaneously.

After this stage, the picture is much sharper.

1. The pure-transfer corridor splits cleanly into two exact 1D branches:
   - numerator-rigid,
   - denominator-rigid.
2. Imposing both rigidities kills the corridor exactly.
3. The two surviving branches are dynamically different:
   - numerator-rigid is a stronger static lever, but it dynamically **splits** the two wall poles;
   - denominator-rigid is a weaker static lever, but it dynamically **penalizes both** wall poles.
4. On the concrete compatibility point, however, **neither** split is killed by the actual dynamic window.
5. The first true ceiling is still the transported **static** one.

So the next clean question is no longer:

> does numerator rigidity or denominator rigidity kill the mechanism?

The answer is no.

The sharper next question is:

> which of those two structural load splits does the real PDE-selected mixed branch most closely realize?

That is the right continuation point after Stage 011.
