# Moving-Throat PDE — Stage 247: Relaxed Stationary Barrier Compiler from One-Port Short-Range, Leakage, `U/V` Drain, and Compensated Source Packets

## Status

**Exact within**

1. the carried one-port static mixed-bundle theorem that restricts the admissible same-charge short-range families to
   
   
\[
   r^{-6},
   \qquad
   e^{-2\kappa r}/r^4,
   \qquad
   e^{-4\kappa r}/r^2,
   \]
2. the exact Stage-244 leakage / work compiler,
3. the exact Stage-245 non-rigid `U/V` drain compiler,
4. the exact Stage-246 compensated source packet,
5. and the declared **stationary barrier-unit embedding** in which the exported, drained, and source-shifted scalars enter the reduced same-charge barrier with nonnegative weights.

The **Session-I benchmark specialization** used below additionally sets the work-channel embedding to unit weight and reads the weighted `U/V` lowering directly from the recorded stationary run.

This stage does **not** derive a new same-charge kernel class.
It is the first honest downstream compiler that takes the already-proved one-port short-range audit and the three relaxed packets from Stages 244–246 and assembles them into one explicit stationary lowered-barrier law.

---

## Purpose

Stage 244 compiled the open-system support-side packet
\[
(S_{\rm leak},\,\mathcal W_w),
\]
Stage 245 compiled the orbit-side non-rigid packet
\[
(U,\,V,\,\mathcal D_{UV}),
\]
and Stage 246 compiled the source-side compensated packet
\[
(\sigma_{\min},\,\mathfrak g[\sigma],\,\mathcal S[\sigma],\,\mathcal R[\sigma]).
\]

What was still missing was the **stationary compiler** that tells the derivation stack how these packets lower the reduced same-charge barrier **after** the one-port audit has already fixed the admissible short-range families.

That is exactly what Stage 247 does.

It has four jobs:

1. restate the carried one-port short-range baseline in a form ready for barrier work,
2. import the Stage-244 leakage / work packet into barrier units,
3. import the Stage-245 weighted `U/V` drain and the Stage-246 compensated source response,
4. assemble the minimal reduced stationary compiler and benchmark it against the recorded Session-I softening point.

---

## Provenance

This stage sits directly after:

- **Stage 244**, which attached leakage and scalar-photon work to the selected-support packet,
- **Stage 245**, which attached a positive `U/V` drain to the non-rigid orbit-side packet,
- **Stage 246**, which converted the sign-changing source branch into explicit mouth/source observables,
- the one-port same-charge barrier audit, which already proved that the static mixed bundle does **not** generate a new long-range attractive family,
- and the Session-I stationary relaxed-constraint run summarized in Appendix B.2 of the barrier write-up.

So Stage 247 is the place where the three relaxed lanes stop floating as separate diagnostics and become one stationary reduced barrier law.

---

## 0. Why this stage is needed

Before this step, the derivation stack had four disconnected pieces:

1. a one-port short-range static baseline,
2. an open-system support packet,
3. an orbit-side non-rigid drain packet,
4. a compensated sign-changing source packet.

That was no longer enough.
The barrier session had already shown a substantial stationary lowering,
\[
\frac{V_{\rm eff}(r_{\rm soft})}{V_{\rm Coul}(r_{\rm soft})}\approx 0.31446203,
\]
but there was not yet a single derivation-stack formula saying **how** the already-compiled packets combine to produce that lowered branch.

Stage 247 closes that gap.

---

## 1. Carried one-port short-range baseline

Let the reduced one-port static bundle be controlled by
\[
\Delta=\Omega_U^2\Omega_W^2-R_{\rm mix}^2,
\]
\[
Q=G_U^2\Omega_W^2+2G_UG_WR_{\rm mix}+G_W^2\Omega_U^2,
\]
\[
P=\Omega_U^2G_W+R_{\rm mix}G_U,
\]
\[
D_0=K_* - \frac{Q}{\Delta}.
\]

The exact static susceptibilities are
\[
\chi_{qq}=\frac1{D_0},
\qquad
\chi_{qU}=\frac{P_U}{\Delta D_0},
\qquad
\chi_{qW}=\frac{P}{\Delta D_0},
\]
\[
\chi_{UU}=\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
\qquad
\chi_{UW}=\frac{K_*R_{\rm mix}+G_UG_W}{\Delta D_0},
\qquad
\chi_{WW}=\frac{K_*\Omega_U^2-G_U^2}{\Delta D_0},
\]
with
\[
P_U=G_U\Omega_W^2+R_{\rm mix}G_W.
\]

Using the primitive short-range source amplitudes
\[
\beta_Q,\qquad \beta_U,\qquad \beta_W,
\]
the exact one-port product-family coefficients are
\[
\mathcal C_6=\chi_{qq}\beta_Q^2,
\]
\[
\mathcal C_4=\chi_{qU}\beta_Q\beta_U+\chi_{qW}\beta_Q\beta_W,
\]
\[
\mathcal C_2=\chi_{UU}\beta_U^2+2\chi_{UW}\beta_U\beta_W+\chi_{WW}\beta_W^2.
\]

If the carried pre-existing core terms are written as
\[
-3\alpha_6 r^{-6},
\qquad
-\alpha_2 e^{-4\kappa r}/r^2,
\]
then the exact one-port short-range baseline may be collected as
\[
A_6=3\alpha_6+\frac12\mathcal C_6,
\qquad
A_4=\mathcal C_4,
\qquad
A_2=\alpha_2+\frac12\mathcal C_2,
\]
so that
\[
\boxed{
V_{\rm short}^{(1p)}(r)
=
\frac1r\left(1+\frac12 e^{-2\kappa r}\right)
-
\frac{A_6}{r^6}
-
A_4\frac{e^{-2\kappa r}}{r^4}
-
A_2\frac{e^{-4\kappa r}}{r^2}.
}
\]

This is the carried short-range baseline that Stage 247 lowers further.

### Structural reading

This formula is already under the barrier-audit firewall:

- it contains the admissible one-port short-range families,
- it does **not** contain a new attractive `1/r` law,
- and it is already tied to the same denominator pair `(\Delta,D_0)` that also feeds the outgoing-normalization stack.

So Stage 247 is downstream of the one-port verdict, not an attempt to reopen it.

---

## 2. Imported relaxed packets

### 2.1 Stage-244 leakage and reduced Session-I work scalar

For compactness, write
\[
\mathfrak L(r):=\Lambda(r)\varrho(r).
\]
Then Stage 244 gives the exact selected-branch support packet
\[
\boxed{
S_{\rm leak}(r)
=
\frac{8\sqrt2\,\eta_{\rm leak}\mu_w\rho_0}{\pi^{5/2}\lambda^3}
\,\mathfrak L(r),
}
\]
\[
\boxed{
\mathcal W_w^{\rm sess}(r)
=
\frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
\,\mathfrak L(r)^2.
}
\]

This stage uses the **reduced Session-I work scalar** `\mathcal W_w^{\rm sess}` because that is the object actually tabulated in the Session-I stationary run under the label `J^wE_w`.

### 2.2 Stage-245 weighted `U/V` lowering scalar

Stage 245 gives the exact drain
\[
\mathcal D_{UV}(r)=\frac{\chi_{UV}(r)^2 a_V f_U(r)^2}{\Delta_{UV}(r)^2},
\qquad
\Delta_{UV}=a_Ua_V-\chi_{UV}^2.
\]

The natural stationary barrier-lowering scalar is its weighted version
\[
\boxed{
\Delta E_{UV}(r):=\eta_{UV}\mathcal D_{UV}(r),
}
\]
where `\eta_{UV}` is the same Session-I barrier-lowering weight carried in the stationary script.

### 2.3 Stage-246 compensated source response

For the transported compensated two-mode source
\[
\sigma(x;r)=1+a(r)\cos(\pi x)+b(r)\cos(2\pi x),
\qquad
a(r)=a_0 s(r),
\qquad
b(r)=b_0 s(r),
\]
with
\[
s(r)=\frac{r_\sigma^2}{r^2+r_\sigma^2},
\]
Stage 246 gives the exact mouth-bias and loading ratio. For the stationary compiler it is convenient to use the monotone source-response scalar
\[
\boxed{
\mathcal M_\sigma(r):=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r)\bigr],
}
\]
with
\[
\mathcal R(r)=\frac{(\mathfrak g(r)-\mathfrak r_{F1})^2}{1+\mathfrak r_{F1}^2},
\qquad
\mathfrak g(r)=\frac{2}{\pi}\left(1+\frac{a(r)}{3}-\frac{b(r)}{15}\right),
\]
\[
\mathcal R_\infty=
\frac{\left(\frac{2}{\pi}-\mathfrak r_{F1}\right)^2}{1+\mathfrak r_{F1}^2}.
\]

On the Session-I orientation
\[
a_0>0,
\qquad
b_0<0,
\]
one has
\[
\mathfrak g(r)\ge \frac{2}{\pi},
\]
and, throughout the sampled branch, still
\[
\mathfrak g(r)<\mathfrak r_{F1}.
\]
Hence `\mathcal R(r)` decreases below its far-field value, so
\[
\boxed{
\mathcal M_\sigma(r)\ge 0
\quad\text{on the Session-I compensated branch.}
}
\]

This is the exact source-side lowering scalar used below.

---

## 3. Exact relaxed stationary barrier compiler

The minimal stationary compiler is now
\[
\boxed{
V_{\rm eff}^{(247)}(r)
=
V_{\rm short}^{(1p)}(r)
-
\lambda_L S_{\rm leak}(r)
-
\lambda_W \mathcal W_w^{\rm sess}(r)
-
\Delta E_{UV}(r)
-
\mathcal M_\sigma(r).
}
\]

Here:

- `\lambda_L\ge 0` converts the open-subsystem leakage source into stationary barrier units,
- `\lambda_W\ge 0` converts the reduced Session-I work scalar into stationary barrier units,
- `\Delta E_{UV}` is already the weighted orbit-side drain scalar,
- `\mathcal M_\sigma` is already the weighted source-response scalar.

### 3.1 Exact lowering identity

By direct subtraction,
\[
\boxed{
V_{\rm short}^{(1p)}(r)-V_{\rm eff}^{(247)}(r)
=
\lambda_L S_{\rm leak}(r)
+
\lambda_W \mathcal W_w^{\rm sess}(r)
+
\Delta E_{UV}(r)
+
\mathcal M_\sigma(r).
}
\]

So if the imported packet scalars are nonnegative on the branch, then
\[
\boxed{
V_{\rm eff}^{(247)}(r)\le V_{\rm short}^{(1p)}(r).
}
\]

That is the exact stationary lowering theorem of the stage.

### 3.2 Support / orbit / source split

The compiler stays block-structured:

- `V_{\rm short}^{(1p)}` is the carried one-port short-range baseline,
- `(S_{\rm leak},\mathcal W_w^{\rm sess})` are support-side / open-system objects,
- `\Delta E_{UV}` is the orbit-side non-rigid drain,
- `\mathcal M_\sigma` is the source-side compensated mouth response.

So Stage 247 does **not** mix the three relaxed lanes algebraically. It composes them downstream of the one-port theorem.

### 3.3 No-new-kernel reading

This is important enough to say plainly.

The one-port audit already proved that the admissible static kernel families are exhausted by the short-range baseline `V_{\rm short}^{(1p)}`. The new imported terms in Stage 247 are **branch-local packet scalars**, not new convolution kernels. So the compiler lowers the barrier by:

- exporting energy into the open-system leakage/work lane,
- draining energy into the dressing leg,
- and shifting the compensated mouth/source branch,

but **not** by introducing any new attractive `1/r`-type law.

So the live corridor remains a **short-range/open-system bypass**.

---

## 4. Session-I one-point benchmark at the strongest softening point

Appendix B.2 of the barrier-session write-up records the strongest stationary softening point at
\[
r_{\rm soft}=0.18,
\qquad
V_{\rm eff}(r_{\rm soft})=1.74701126,
\qquad
V_{\rm Coul}(r_{\rm soft})=5.55555556,
\]
so that
\[
\frac{V_{\rm eff}(r_{\rm soft})}{V_{\rm Coul}(r_{\rm soft})}=0.31446203.
\]

Using the recorded one-port parameters
\[
K_*=4.0,
\quad
\Omega_U^2=9.0,
\quad
\Omega_W^2=16.0,
\quad
G_U=1.0,
\quad
G_W=1.25,
\quad
R_{\rm mix}=1.35,
\]
\[
\beta_Q=0.03,
\quad
\beta_{U0}=0.15,
\quad
\beta_{W0}=0.20,
\quad
\kappa=1,
\]
and taking the session specialization `\alpha_6=\alpha_2=0` for this benchmark slice, the exact Stage-247 baseline gives
\[
\Delta=142.17750000,
\qquad
D_0=3.76481862,
\]
\[
\boxed{
V_{\rm short}^{(1p)}(r_{\rm soft})=3.74163698.
}
\]

So the one-port short-range baseline alone lowers the Coulomb barrier to about
\[
\frac{V_{\rm short}^{(1p)}(r_{\rm soft})}{V_{\rm Coul}(r_{\rm soft})}
\approx 0.67349466,
\]
which is **not yet** the full Session-I softening.

That remaining drop is exactly what the imported Stage-244/245/246 packets must supply.

### 4.1 Imported packet values on the benchmark slice

The Session-I table also records
\[
\mathcal W_w^{\rm sess}(r_{\rm soft})=1.51632107,
\qquad
\Delta E_{UV}(r_{\rm soft})=0.21064278.
\]

From the exact Stage-244 work law,
\[
\mathcal W_w^{\rm sess}(r)
=
\frac{512\,\eta_{\rm leak}^2\mu_w q\rho_0}{\pi^4\lambda^2}
\,\mathfrak L(r)^2,
\]
with the Session-I parameters
\[
\lambda=1,
\qquad
\eta_{\rm leak}=0.03,
\qquad
\mu_w=0.8,
\qquad
q=1,
\qquad
\rho_0=1,
\]
one infers
\[
\boxed{
\mathfrak L(r_{\rm soft})=\Lambda(r_{\rm soft})\varrho(r_{\rm soft})
=20.01677473,
}
\]
and therefore
\[
\boxed{
S_{\rm leak}(r_{\rm soft})=0.31069599.
}
\]

For the compensated source packet, using
\[
a_0=2.2,
\qquad
b_0=-0.6,
\qquad
r_\sigma=0.8,
\qquad
\xi_R=0.9,
\qquad
\mathfrak r_{F1}=1.77799353547498,
\]
one finds
\[
\boxed{
\mathcal M_\sigma(r_{\rm soft})=\xi_R\bigl[\mathcal R_\infty-\mathcal R(r_{\rm soft})\bigr]
=0.18386120.
}
\]

### 4.2 Exact benchmark decomposition

Choose the natural Session-I work specialization
\[
\lambda_W=1.
\]
Then the Stage-247 compiler requires only one additional leakage embedding coefficient,
\[
\lambda_L,
\]
which is fixed exactly by the recorded softening point:
\[
\boxed{
\lambda_L
=
\frac{V_{\rm short}^{(1p)}(r_{\rm soft})
-
\mathcal W_w^{\rm sess}(r_{\rm soft})
-
\Delta E_{UV}(r_{\rm soft})
-
\mathcal M_\sigma(r_{\rm soft})
-
V_{\rm eff}^{\rm sess}(r_{\rm soft})}
{S_{\rm leak}(r_{\rm soft})}.
}
\]

Numerically,
\[
\boxed{
\lambda_L=0.26971918.
}
\]

So the benchmark decomposition is
\[
3.74163698
-
0.26971918\times 0.31069599
-
1.51632107
-
0.21064278
-
0.18386120
=
1.74701126.
\]

This is the most useful concrete result of the stage:

> once the one-port short-range baseline is fixed, the recorded Session-I stationary lowering can be decomposed exactly into
> 
> 1. a carried one-port short-range contribution,
> 2. a reduced Session-I work contribution,
> 3. a weighted `U/V` drain contribution,
> 4. a compensated source-response contribution,
> 5. and one remaining leakage-to-barrier embedding coefficient.

### 4.3 Residual structure of the benchmark

The benchmark also makes the packet split numerically transparent:
\[
V_{\rm short}^{(1p)}(r_{\rm soft})-
\mathcal W_w^{\rm sess}(r_{\rm soft})-
\Delta E_{UV}(r_{\rm soft})
=2.01467313,
\]
so work plus orbit-side drain alone do **not** yet reach the recorded softening.
After adding the compensated source packet,
\[
V_{\rm short}^{(1p)}(r_{\rm soft})-
\mathcal W_w^{\rm sess}(r_{\rm soft})-
\Delta E_{UV}(r_{\rm soft})-
\mathcal M_\sigma(r_{\rm soft})
=1.83081193,
\]
leaving a final gap of
\[
0.08380067,
\]
which is exactly what the leakage term closes when
\[
\lambda_L=0.26971918.
\]

So the stationary readback is now algebraically sharp.

---

## 5. What this stage achieves physically

Stage 247 closes the first true stationary endgame gap on the relaxed same-charge branch.

1. It keeps the one-port long-range firewall intact.
   The short-range baseline remains the only carrier of admissible same-charge kernel families.
2. It converts the Stage-244, Stage-245, and Stage-246 packets into one barrier-level compiler.
3. It shows that the Session-I lowered branch can be decomposed into short-range, open-system, orbit-drain, and compensated-source pieces without inventing a new asymptotic attraction.
4. It reduces the remaining stationary embedding freedom, on the benchmark slice, to a single leakage-to-barrier coefficient once the work channel is taken in the Session-I normalization.

So the same-charge stationary corridor is no longer just “the barrier got lower in a script.”
It is now a derivation-stack object with an explicit packet decomposition.

---

## 6. Result

The first honest relaxed stationary barrier compiler is now in place.

The exact formula is
\[
\boxed{
V_{\rm eff}^{(247)}(r)
=
V_{\rm short}^{(1p)}(r)
-
\lambda_L S_{\rm leak}(r)
-
\lambda_W \mathcal W_w^{\rm sess}(r)
-
\Delta E_{UV}(r)
-
\mathcal M_\sigma(r).
}
\]

Its lowering identity is exact,
\[
V_{\rm short}^{(1p)}-V_{\rm eff}^{(247)}
=
\lambda_L S_{\rm leak}+\lambda_W\mathcal W_w^{\rm sess}+\Delta E_{UV}+\mathcal M_\sigma,
\]
and the Session-I strongest-softening point is reproduced by the concrete benchmark decomposition above.

So Stage 247 gives the correct stationary continuation of the relaxed branch:

- **one-port short-range baseline**,
- **support-side open-system export**,
- **orbit-side non-rigid drain**,
- **source-side compensated mouth response**,
- and **no new long-range same-charge law**.

---

## 7. Immediate next step

The next stage should promote this stationary compiler into the dynamic event chain.

That means:

1. take `V_{\rm eff}^{(247)}(r)` as the lowered stationary front end,
2. compute the peak, turning point, threshold speed, and WKB action on that lowered branch,
3. then carry the same short-range/open-system discipline into the dynamic scattering problem.

That is exactly the Stage-248 job.
