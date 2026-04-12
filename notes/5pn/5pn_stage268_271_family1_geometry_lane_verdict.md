# 5PN continuation — Stages 268–271

## Goal

Continue from the Stage-267 corridor theorem without oversimplifying the scalar/geometry lane.

The actual next question was:

> what does the first *known* scalar/geometry lane from the moving-throat / Family-1 program look like **at the level of the Stage-78 contamination coefficients** `(a_2,a_4)`, and does it really look dangerous enough to threaten the grouped-`P2` twin-safe / Family-1-safe corridor?

The answer from this batch is much sharper than a generic “maybe.” The explicit Family-1 scalar lane is a **coupled two-pole monopole breathing channel**. For the Stage-78 contamination coefficients, that two-pole channel has an **exact contamination-equivalent one-pole diagnostic** with

\[
\lambda_{\rm coeff}=\frac{G_2}{G_4}\approx 28.1365336267,
\qquad
c_{\rm eff}=\frac{G_2^2}{G_0G_4}\approx 0.6072383642.
\]

So the known scalar lane is neither pure pole nor pure contact at the contamination level. But on the actual isotropic moving-throat branch the geometry lane is dynamically inert through \(O(\omega^4)\), so this scalar lane only re-enters after an explicit \(l=0\leftrightarrow l=2\) mixing source is turned on. The actual isotropic branch itself remains exactly safe and the active reduced bottleneck remains the outgoing quadrupole normalization defect \(N_Q\).

---

## Stage 268 — exact effective one-pole diagnostic for arbitrary scalar lanes

Start from the exact Stage-78 contamination pair
\[
a_2=\frac{\Omega_Q^2M_0^2G_2}{K_{\rm pole}G_0^2},
\qquad
a_4=\frac{\Omega_Q^4M_0^2(G_0G_4-G_2^2)}{K_{\rm pole}G_0^3}.
\]

These depend only on the first three scalar-lane moments \((G_0,G_2,G_4)\). That means there is an exact “moment-equivalent” one-pole diagnostic for *any* scalar lane:

\[
\Omega_{\rm eff}^2=\frac{G_2}{G_4},
\qquad
c_{\rm eff}=\frac{G_2^2}{G_0G_4},
\qquad
G_{{\rm pole},\rm eff}=\frac{G_2^2}{G_4},
\qquad
G_{{\rm contact},\rm eff}=G_0-\frac{G_2^2}{G_4}.
\]

Then
\[
a_2
=
\left(\frac{M_0^2}{K_{\rm pole}G_0}\right)
\left(\frac{\Omega_Q^2}{\Omega_{\rm eff}^2}\right)c_{\rm eff},
\]
\[
a_4
=
\left(\frac{M_0^2}{K_{\rm pole}G_0}\right)
\left(\frac{\Omega_Q^2}{\Omega_{\rm eff}^2}\right)^2
c_{\rm eff}(1-c_{\rm eff}),
\]
**exactly**.

So the Stage-266 corridor theorem applies to any scalar/geometry lane through the exact effective one-pole diagnostic \((\Omega_{\rm eff}^2,c_{\rm eff})\). This is not a fit; it is an identity for the contamination pair \((a_2,a_4)\).

### Exact two-pole decomposition

For a genuine contact-plus-two-pole scalar lane
\[
D_g = G_c+\frac{R_1}{1-s/L_1}+\frac{R_2}{1-s/L_2},
\]
the obstruction term splits exactly as
\[
G_0G_4-G_2^2
=
G_c\!\left(\frac{R_1}{L_1^2}+\frac{R_2}{L_2^2}\right)
+
R_1R_2\!\left(\frac1{L_1}-\frac1{L_2}\right)^2.
\]

So a two-pole lane can pick up “effective contact” both from the literal contact slot and from pole separation.

---

## Stage 269 — the actual Family-1 monopole breathing channel

The 2PN constructive appendix already fixed the scalar monopole breathing channel as
\[
K_{00}(s)
=
-\frac{757}{2520}
+\frac{R_-}{1-s/\lambda_-}
+\frac{R_+}{1-s/\lambda_+},
\]
with
\[
\lambda_- \approx 6.405572392138922,
\qquad
\lambda_+ \approx 254.444968136936126,
\]
\[
R_- \approx 0.002552474771738,
\qquad
R_+ \approx 0.386733239513976,
\]
and exact static value
\[
K_{00}(0)=\frac{4}{45}.
\]

So the actual constructive scalar lane is explicitly **two-pole**.

From this we get the first three low-frequency moments
\[
G_0=\frac{4}{45},
\qquad
G_2=\frac{R_-}{\lambda_-}+\frac{R_+}{\lambda_+}
\approx 1.91838640121275\times10^{-3},
\]
\[
G_4=\frac{R_-}{\lambda_-^2}+\frac{R_+}{\lambda_+^2}
\approx 6.81813341566730\times10^{-5}.
\]

The exact contamination-equivalent one-pole diagnostic is therefore
\[
\lambda_{\rm coeff}=\frac{G_2}{G_4}\approx 28.1365336267,
\qquad
c_{\rm eff}=\frac{G_2^2}{G_0G_4}\approx 0.6072383642.
\]

So for the Stage-78 contamination coefficients, the actual Family-1 two-pole breathing lane behaves exactly like a **mixed contact-plus-pole** scalar lane with pole fraction about \(0.61\).

### Important non-equivalence

The 2PN note also quotes an interval-accuracy Padé reduction
\[
\lambda_{\rm pade}\approx 202.9235163675.
\]

This is **not** the same thing as \(\lambda_{\rm coeff}\). The Padé pole is useful for fitting \(K_{00}(s)\) over a finite \(s\)-interval. The coefficient-equivalent pole is the exact one-pole reduction that preserves \((a_2,a_4)\). Conflating them would be an unnecessary oversimplification.

---

## Stage 270 — exact corridor comparison for the actual Family-1 breathing lane

Using the exact coefficient-equivalent pole fraction \(c_{\rm eff}\), the scalar-lane danger variable becomes

\[
u_{\rm eff}=r_{\rm eff}(1-c_{\rm eff}),
\qquad
r_{\rm eff}:=\frac{\Omega_Q^2}{\lambda_{\rm coeff}},
\]
with
\[
1-c_{\rm eff}\approx 0.3927616358.
\]

So the Stage-266 thresholds become explicit pole-ratio thresholds:

### Initial-drift thresholds
\[
r_{\rm eff}>\frac{2}{1-c_{\rm eff}}
\approx 5.09214703737
\quad\Longrightarrow\quad
\text{support demand grows initially},
\]

\[
r_{\rm eff}>\frac{4}{1-c_{\rm eff}}
\approx 10.1842940747
\quad\Longrightarrow\quad
\text{the twin margin shrinks initially}.
\]

### Actual failure thresholds
\[
r_{\rm eff}\ge
\frac{4+2\sqrt2}{1-c_{\rm eff}}
\approx 17.3856774766
\quad\Longrightarrow\quad
\text{actual twin failure can occur},
\]

and, using the exact Lambda\(_{\rm EM}\)-refreshed hard Family-1 ceiling \(c_\*=0.7116102605\ldots\),
\[
r_{\rm eff}\ge
\frac{8c_\*+4\sqrt{c_\*(4c_\*-1)}}{1-c_{\rm eff}}
\approx 26.1684980784
\quad\Longrightarrow\quad
\text{actual Family-1 failure can occur}.
\]

So the first known constructive scalar lane is not dangerous by composition alone. It becomes dangerous only if the grouped quadrupole pole is **much faster** than the contamination-equivalent breathing pole in the same reduced spectral variable.

---

## Stage 271 — actual branch verdict

Now combine this with the later moving-throat theorem ledger.

On the actual isotropic branch the geometry lane is dynamically inert through \(O(\omega^4)\) with respect to the grouped real `P2` carrier:
\[
\epsilon_2=\epsilon_4=0.
\]

So the conservative grouped-`P2` branch is exactly
\[
c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13
\]
on the actual isotropic branch, and the explicit Family-1 support/source side is automatic there.

That means:

1. the actual isotropic branch is **exactly safe** against scalar/geometry contamination;
2. the first known scalar-lane correction enters only after an explicit \(l=0\leftrightarrow l=2\) mixing source is turned on;
3. that first correction is controlled by the **Family-1 two-pole monopole breathing lane**;
4. and even then it only threatens the corridor if the physical pole ratio \(r_{\rm eff}\) is large enough to exceed the thresholds above.

So the scalar/geometry lane is not the natural source of 5PN failure on the actual isotropic branch.

---

## Best current reading after Stages 268–271

The next theorem gate is now much smaller.

### What is now settled

- Any scalar/geometry lane has an exact contamination-equivalent one-pole diagnostic \((\Omega_{\rm eff}^2,c_{\rm eff})\) for the pair \((a_2,a_4)\).
- The actual constructive scalar lane from the 2PN Family-1 module is a **two-pole monopole breathing channel**, not a literal one-pole lane.
- For contamination purposes, that actual scalar lane is exactly equivalent to
  \[
  \lambda_{\rm coeff}\approx 28.1365,
  \qquad
  c_{\rm eff}\approx 0.60724.
  \]
- On the actual isotropic moving-throat branch,
  \[
  \epsilon_2=\epsilon_4=0
  \]
  exactly, so the scalar lane does not contaminate the grouped `P2` carrier there at all.
- The explicit Family-1 support/source side is already automatic on that same actual isotropic branch.

### What remains genuinely open

The scalar/geometry lane only re-enters after an explicit \(l=0\leftrightarrow l=2\) mixing source is turned on. At that point the only remaining scalar-lane uncertainty is the physical pole ratio
\[
r_{\rm eff}=\frac{\Omega_Q^2}{\lambda_{\rm coeff}}.
\]

So the next honest continuation is now:

> determine, on the actual weakly anisotropic moving-throat reduced operator, whether the physical grouped-quadrupole pole is anywhere near the explicit danger window
> \[
> r_{\rm eff}\gtrsim 17.39 \quad\text{(twin failure)},
> \qquad
> r_{\rm eff}\gtrsim 26.17 \quad\text{(Family-1 failure)}.
> \]

Unless it is, the first known scalar/geometry lane is a controlled correction, not the active 5PN bottleneck.
