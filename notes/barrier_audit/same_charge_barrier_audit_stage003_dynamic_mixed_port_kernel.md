# Same-Charge Barrier Audit — Stage 003: Dynamic Mixed-Port Kernel, Phase-Lag No-Go, and the Resonant-Survival Gate

## 0. Purpose

Stage 002 closed the first honest **static** mixed-bundle audit.
It showed that the minimal one-port wall / Maxwell / mixed bundle does not create a new long-range attractive law. It only generates short-range product families built from the primitive source profiles.

That left one clear escape hatch:

> perhaps a genuinely time-dependent / non-adiabatic mixed-sector drive can do something qualitatively stronger than the static bundle.

This stage is the first clean test of that possibility.

The job here is not to solve the full driven two-throat PDE.
The job is narrower:

1. promote the Stage-002 one-port bundle to a **frequency-dependent** mixed-port kernel,
2. derive its exact complex susceptibility,
3. ask whether linear monochromatic driving creates a genuinely new spatial attractive family,
4. separate the in-phase barrier reshaping from the out-of-phase pumping / leakage channel,
5. and identify the only live dynamic corridor that survives the linear audit.

The result is again restrictive:

> at linear order, the dynamic mixed bundle does **not** create a new spatial kernel class.
> It only makes the same short-range families frequency dependent, and the first outgoing-port correction is purely phase-lag / pumping rather than conservative barrier lowering.

So if the idea survives dynamically, it survives through a **resonant dispersive** corridor or a genuinely nonlinear / parametric corridor, not through a linear monochromatic barrier bypass.

---

## 1. Frozen input carried forward

### 1.1 Static one-port bundle from Stage 002

After the stable BdG support mode is integrated into the static wall stiffness,

\[
K_* = K - \frac{C^2}{\varpi^2},
\]

the static one-port bundle is controlled by

\[
\Delta = \Omega_U^2\Omega_W^2 - R^2,
\]
\[
Q = G_U^2\Omega_W^2 + 2 G_U G_W R + G_W^2\Omega_U^2,
\]
\[
P = \Omega_U^2 G_W + R G_U,
\]
\[
D_0 = K_* - \frac{Q}{\Delta},
\]

with exact outgoing-prefactor bridge

\[
N_0 = \frac{P^2}{\Delta^2},
\qquad
P_0 = \frac{N_0}{D_0}.
\]

Stage 002 proved that the static mixed correction is always attractive or neutral on the admissible branch and that, for the primitive reduced source families,

\[
\mathcal S_Q(x)=\frac1{x^3},
\qquad
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
\]

the bundle can produce only the three product families

\[
\frac1{x^6},
\qquad
\frac{e^{-2\kappa x}}{x^4},
\qquad
\frac{e^{-4\kappa x}}{x^2}.
\]

So the linear dynamic audit must be judged against that static baseline.

### 1.2 Why the dynamic mixed sector is the right next lane

The 4D plasma extension keeps the mixed fields

\[
A_w,
\qquad
F_{\mu w},
\qquad
J^w,
\]

alive once one relaxes the strict far-field brane reduction, and it treats brane-facing non-ideality as conservative transport into the hidden `w` sector and higher localization modes rather than as purely local dissipation. So if a field-driven corridor survives at all, it has to pass through the dynamic mixed sector rather than through the static zero-mode Maxwell law.

### 1.3 Outgoing-port motivation

The moving-throat program already showed where a genuine passive/outgoing channel first appears in the reduced wall language: the mixed-sector port. On the outgoing `l=2` branch the first odd fingerprint starts at `i omega^5` rather than at an even static correction. So the same reduced sector that can carry dynamic pumping is also the place where the radiative quadrupole bridge lives.

That is exactly why the dynamic barrier audit has to distinguish **in-phase real** response from **out-of-phase imaginary** response.

---

## 2. Minimal dynamic one-port mixed bundle

Take the same reduced coordinates as in Stage 002:

- wall/worldtube amplitude `q`,
- brane-like internal gauge coordinate `U`,
- mixed `A_w / F_(mu w) / J^w` coordinate `W`.

Now keep the frequency dependence explicitly.

### 2.1 Dynamic reduced coefficients

Define

\[
K_B(\omega)
:=
K - M\omega^2 - \frac{C^2}{\varpi^2-\omega^2},
\]
\[
A(\omega):=\Omega_U^2-\omega^2,
\qquad
W(\omega):=\Omega_W^2-\omega^2-\Pi(\omega),
\]

where `Pi(omega)` is the effective outgoing / port self-energy carried by the mixed coordinate.

Then define

\[
\Delta_\Pi(\omega):=A(\omega)W(\omega)-R^2,
\]
\[
Q_\Pi(\omega):=G_U^2W(\omega)+2G_UG_WR+G_W^2A(\omega),
\]
\[
D_\Pi(\omega):=K_B(\omega)-\frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)}.
\]

### 2.2 Dynamic reduced stiffness matrix

The frequency-domain reduced bundle is

\[
\mathcal K_{\rm dyn}(\omega)
=
\begin{pmatrix}
K_B(\omega) & -G_U & -G_W \\
-G_U & A(\omega) & -R \\
-G_W & -R & W(\omega)
\end{pmatrix}.
\]

Its exact determinant identity is

\[
\boxed{
\det \mathcal K_{\rm dyn}(\omega)
=
\Delta_\Pi(\omega)\,D_\Pi(\omega).
}
\]

So the same two denominator factors govern the dynamic bundle:

- the internal `(U,W)` block denominator `Delta_Pi`,
- the dressed wall denominator `D_Pi`.

### 2.3 Static reduction check

At zero frequency with no outgoing dressing,

\[
\omega=0,
\qquad
\Pi(0)=0,
\]

one has

\[
K_B(0)=K-\frac{C^2}{\varpi^2}=K_*,
\qquad
A(0)=\Omega_U^2,
\qquad
W(0)=\Omega_W^2,
\]

so the dynamic bundle reduces exactly to the Stage-002 static bundle.

So Stage 003 is a genuine continuation, not a different reduced model.

---

## 3. Exact dynamic susceptibility law

Define the complex reduced quadratic response functional

\[
\mathfrak V_{\rm mix}(x,\omega)
:=
-\frac12 J(x,\omega)^T\mathcal K_{\rm dyn}(\omega)^{-1}J(x,\omega).
\]

This is the direct dynamic analogue of the Stage-002 static on-shell energy shift.

- `Re \mathfrak V_mix` is the in-phase conservative barrier reshaping.
- `Im \mathfrak V_mix` is the quadrature / pumping / leakage channel.

### 3.1 Exact inverse entries

The inverse has the same form as in Stage 002, but with the dynamic coefficients substituted:

\[
\chi_{qq}(\omega)=\frac1{D_\Pi(\omega)},
\]
\[
\chi_{qU}(\omega)=\frac{P_U(\omega)}{\Delta_\Pi(\omega)D_\Pi(\omega)},
\qquad
P_U(\omega):=G_UW(\omega)+RG_W,
\]
\[
\chi_{qW}(\omega)=\frac{P(\omega)}{\Delta_\Pi(\omega)D_\Pi(\omega)},
\qquad
P(\omega):=A(\omega)G_W+RG_U,
\]
\[
\chi_{UU}(\omega)=\frac{K_B(\omega)W(\omega)-G_W^2}{\Delta_\Pi(\omega)D_\Pi(\omega)},
\]
\[
\chi_{UW}(\omega)=\frac{K_B(\omega)R+G_UG_W}{\Delta_\Pi(\omega)D_\Pi(\omega)},
\]
\[
\chi_{WW}(\omega)=\frac{K_B(\omega)A(\omega)-G_U^2}{\Delta_\Pi(\omega)D_\Pi(\omega)}.
\]

So every dynamic same-charge correction is still controlled by the same denominator pair `Delta_Pi D_Pi`.

### 3.2 Collinear-source theorem

If the reduced source is collinear in source space,

\[
J(x,\omega)=\mathcal S(x,\omega)
\begin{pmatrix}
s_q\\ s_U\\ s_W
\end{pmatrix},
\]

then the response still factorizes exactly:

\[
\boxed{
\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\,\chi_s(\omega)\,\mathcal S(x,\omega)^2,
\qquad
\chi_s(\omega)=\frac{\mathcal N_s(\omega)}{\Delta_\Pi(\omega)D_\Pi(\omega)}.
}
\]

Here

\[
\mathcal N_s(\omega)
=
\Delta_\Pi s_q^2
+2P_U s_q s_U
+2P s_q s_W
+\bigl(K_BW-G_W^2\bigr)s_U^2
+2\bigl(K_BR+G_UG_W\bigr)s_Us_W
+\bigl(K_BA-G_U^2\bigr)s_W^2.
\]

So the dynamic bundle still acts by a scalar susceptibility on each collinear reduced source family.

---

## 4. Dynamic product-family theorem for the primitive same-charge loads

Now return to the same primitive reduced source profiles used in Stage 002:

\[
\mathcal S_Q(x)=\frac1{x^3},
\qquad
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x}.
\]

Take the reduced source vector

\[
J(x,\omega)
=
\begin{pmatrix}
\beta_Q\mathcal S_Q(x)\\
\beta_U\mathcal S_Y(x)\\
\beta_W\mathcal S_Y(x)
\end{pmatrix}.
\]

Then the exact dynamic correction is

\[
\boxed{
\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\left[
\frac{\mathcal C_6(\omega)}{x^6}
+
2\mathcal C_4(\omega)\frac{e^{-2\kappa x}}{x^4}
+
\mathcal C_2(\omega)\frac{e^{-4\kappa x}}{x^2}
\right],
}
\]

with

\[
\mathcal C_6(\omega)=\chi_{qq}(\omega)\,\beta_Q^2,
\]
\[
\mathcal C_4(\omega)=\chi_{qU}(\omega)\beta_Q\beta_U+
\chi_{qW}(\omega)\beta_Q\beta_W,
\]
\[
\mathcal C_2(\omega)=\chi_{UU}(\omega)\beta_U^2
+2\chi_{UW}(\omega)\beta_U\beta_W
+\chi_{WW}(\omega)\beta_W^2.
\]

This is the first strong dynamic theorem of the audit:

> the linear monochromatic mixed bundle does **not** create a new spatial kernel family.
> It keeps exactly the same three short-range spatial families as the static one-port bundle, but promotes their coefficients to complex frequency-dependent susceptibilities.

So linear time dependence does **not** yet buy a qualitatively new radial law.

It only buys:

1. frequency-dependent renormalization of the same short-range families,
2. complex phase structure,
3. and the possibility of resonant amplification through the denominator pair `Delta_Pi D_Pi`.

---

## 5. Exact outgoing-port derivative identity

The real dynamic novelty is not the source-product structure. It is the outgoing port.

Let

\[
e_W=(0,0,1)^T,
\]

so that the mixed coordinate is the `W` lane. Because `Pi(omega)` enters only in the `WW` slot,

\[
\partial_\Pi \mathcal K_{\rm dyn}(
\omega)
=
- e_W e_W^T.
\]

Matrix calculus then gives the exact identity

\[
\boxed{
\partial_\Pi \mathcal K_{\rm dyn}(\omega)^{-1}
=
\mathcal K_{\rm dyn}(\omega)^{-1}
\,e_W e_W^T\,
\mathcal K_{\rm dyn}(\omega)^{-1}.
}
\]

Therefore the dynamic response functional obeys

\[
\boxed{
\partial_\Pi \mathfrak V_{\rm mix}(x,\omega)
=
-\frac12
\bigl[e_W^T\mathcal K_{\rm dyn}(\omega)^{-1}J(x,\omega)\bigr]^2.
}
\]

So for a small outgoing port around the conservative branch `Pi = 0`,

\[
\boxed{
\delta\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\Pi(\omega)
\,\mathcal T_J(\omega)^2
+O\bigl(\Pi(\omega)^2\bigr),
}
\]

with exact transfer amplitude

\[
\mathcal T_J(\omega)
:=
 e_W^T\mathcal K_{\rm cons}(\omega)^{-1}J(x,\omega),
\qquad
\mathcal K_{\rm cons}:=\mathcal K_{\rm dyn}\big|_{\Pi=0}.
\]

That is the exact dynamic generalization of the Stage-004 outgoing-transfer factor from the moving-throat quadrupole bridge.

---

## 6. Linear phase-lag no-go theorem

Now impose the minimal passive/outgoing form

\[
\Pi(\omega)=i\,\Gamma_{\rm out}(\omega),
\qquad
\Gamma_{\rm out}(\omega)\ge 0,
\]

with the understanding that on the outgoing quadrupole branch one has

\[
\Gamma_{\rm out}(\omega)
\sim
\frac{a^5}{27c_s^5}\,\omega^5
\qquad (l=2\ \text{fingerprint}).
\]

On the conservative branch the reduced matrix is real for real `omega` away from its poles, so `\mathcal T_J(\omega)^2` is real. Therefore

\[
\boxed{
\delta\mathfrak V_{\rm mix}(x,\omega)
=
-\frac{i}{2}\Gamma_{\rm out}(\omega)\,\mathcal T_J(\omega)^2
+O\bigl(\Gamma_{\rm out}^2\bigr).
}
\]

Hence the first outgoing correction satisfies

\[
\boxed{
\Re\,\delta\mathfrak V_{\rm mix}(x,\omega)=0
\qquad\text{to first order in the passive/outgoing port.}
}
\]

But the quadrature load is nonzero:

\[
\Im\,\delta\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\Gamma_{\rm out}(\omega)\,\mathcal T_J(\omega)^2.
\]

For a monochromatic drive, the average absorbed power is therefore

\[
\boxed{
\overline{P}_{\rm abs}^{(1)}(x,\omega)
=
-\omega\,\Im\delta\mathfrak V_{\rm mix}(x,\omega)
=
\frac{\omega}{2}\Gamma_{\rm out}(\omega)\,\mathcal T_J(\omega)^2
\ge 0.
}
\]

This is the core dynamic no-go statement of the stage:

> the leading passive/outgoing correction feeds **pumping / leakage / phase lag**, not conservative barrier lowering.

So a linear outgoing phase term by itself cannot be counted as a same-charge barrier-softening mechanism.
It is an energy-flow channel, not a static or in-phase screening channel.

---

## 7. What still survives dynamically

The phase-lag no-go does **not** kill the whole dynamic corridor.
It only kills the easiest version.

### 7.1 What is dead at linear order

The following picture is now strongly disfavored:

> apply a monochromatic mixed-sector drive, let the leading outgoing phase term turn on, and obtain a first-order conservative barrier reduction from that phase term alone.

The first outgoing term is purely imaginary at linear order.
So that route does not lower the barrier. It only loads the system dynamically.

### 7.2 What is still alive

There are still two linear dynamic possibilities left.

#### (a) Dispersive even-frequency reshaping

The real part of the dynamic susceptibilities

\[
\mathcal C_6(\omega),
\qquad
\mathcal C_4(\omega),
\qquad
\mathcal C_2(\omega)
\]

changes through the even frequency dependence of

\[
K_B(\omega),
\qquad
A(\omega),
\qquad
W(\omega)
\]

and therefore through the pole structure of

\[
\Delta_0(\omega)=\bigl(\Omega_U^2-\omega^2\bigr)\bigl(\Omega_W^2-\omega^2\bigr)-R^2,
\]
\[
D_0(\omega)=K_B(\omega)-\frac{Q_0(\omega)}{\Delta_0(\omega)}.
\]

Away from poles this is just an even analytic deformation of the static coefficients.
Near poles it can be parametrically amplified.

#### (b) Resonant dispersive enhancement

Only in a neighborhood where `Delta_0(omega)` or `D_0(omega)` becomes small can the real in-phase part become much larger than its static value.
But those are exactly the regions where the same transfer factors also amplify the imaginary pumping channel once the passive/outgoing port is restored.

So the dynamic survival test is no longer a kernel-class question.
It is a **resonance / linewidth / quality-factor** question:

> can the real dispersive enhancement of the existing short-range attractive families become large enough **before** the absorptive / leakage channel simply turns the mechanism into disguised heating or branch loss?

### 7.3 What remains beyond this stage

A genuinely new spatial kernel class would require something beyond linear monochromatic response, for example

- nonlinear parametric mixing,
- sideband generation,
- multi-frequency beating,
- or explicit time-dependent modulation of the reduced coefficients themselves.

None of those are part of the current linear one-port audit.

---

## 8. Updated reduced barrier audit after Stage 003

At linear monochromatic order the dynamic same-charge audit potential is best written as the **real** part of the driven response correction:

\[
\widetilde V_{\rm audit}^{(3)}(x,\omega)
=
\frac{1}{x}\Bigl(1+\frac12 e^{-2\kappa x}\Bigr)
-
3\alpha_6\,\mathcal G_6(x)
-
\alpha_2\frac{e^{-4\kappa x}}{x^2}
+
\Re\,\mathfrak V_{\rm mix}(x,\omega),
\]

with

\[
\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\left[
\frac{\mathcal C_6(\omega)}{x^6}
+
2\mathcal C_4(\omega)\frac{e^{-2\kappa x}}{x^4}
+
\mathcal C_2(\omega)\frac{e^{-4\kappa x}}{x^2}
\right].
\]

The dynamic pumping diagnostic is

\[
\overline{P}_{\rm abs}(x,\omega)
=
-\omega\,\Im\mathfrak V_{\rm mix}(x,\omega).
\]

So any linear dynamic claim now has to clear **both** tests:

1. `Re \mathfrak V_mix` must materially reduce the barrier,
2. `P_abs` must stay small enough that the effect is not just covert heating / leakage.

That is the first honest dynamic kill test.

---

## 9. Current verdict after Stage 003

Stage 001 killed the naive Maxwell-only story and showed that pure inverse-sixth attraction cannot erase the barrier by itself.

Stage 002 showed that the first honest static mixed bundle creates no new long-range attractive family — only short-range product kernels tied to the same one-port bundle that already feeds the 5PN / 2.5PN normalization chain.

Stage 003 now sharpens the dynamic side.

1. The linear dynamic one-port bundle still produces only the same three spatial families
   
   \[
   x^{-6},
   \qquad
   e^{-2\kappa x}/x^4,
   \qquad
   e^{-4\kappa x}/x^2.
   \]

2. Time dependence only makes their coefficients complex and frequency dependent.
3. The first passive/outgoing correction is purely phase-lag / pumping at linear order and therefore does **not** lower the barrier conservatively.
4. The only surviving linear dynamic corridor is a **resonant dispersive enhancement** of the already-known short-range families, and that corridor must beat its own absorptive / leakage load to count as real barrier engineering rather than disguised heating.

So the idea is still alive, but now in a much narrower form:

> not “dynamic mixed driving creates a new attractive law,”
> but
> “a narrow resonant mixed-sector window may amplify the existing short-range attractive channels enough to matter before absorption and branch loss dominate.”

---

## 10. Immediate next step

The next theorem gate is now very sharp.

1. Keep the exact Stage-003 dynamic one-port kernel.
2. Choose the first explicit admissible branch data for
   
   \[
   K,
   \quad
   M,
   \quad
   C,
   \quad
   \varpi,
   \quad
   \Omega_U,
   \quad
   \Omega_W,
   \quad
   G_U,
   \quad
   G_W,
   \quad
   R.
   \]

3. Compute the real dispersive enhancement and the absorptive load near the first internal poles.
4. Ask whether any resonance window gives a materially larger
   
   \[
   -\Re\mathfrak V_{\rm mix}
   \]

   before
   
   \[
   \overline{P}_{\rm abs}
   \]

   becomes the dominant effect.

That is the smallest honest continuation point after the linear dynamic audit.
