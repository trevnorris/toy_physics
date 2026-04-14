same_charge_barrier_audit_stage002_one_port_mixed_bundle_kernel.md

# Same-Charge Barrier Audit — Stage 002: One-Port Mixed-Bundle Static Kernel and the Square-Law Suppression Test

## 0. Purpose

Stage 001 left the same-charge survival corridor in the form

\[
-\gamma_{\rm mix}\,\mathcal G_{\rm mix}(x),
\]

with the explicit instruction that the next honest step was to replace the placeholder kernel by the **first kernel actually implied by the moving-throat mixed bundle**.

That is what this stage does.

The goal is not yet to solve the full two-throat PDE.
The goal is narrower and cleaner:

1. take the already isolated isotropic one-port wall/BdG/Maxwell/mixed bundle,
2. derive its exact **static susceptibility kernel** under an external same-charge load,
3. identify how that kernel is controlled by the same quantities that already appear in the 5PN / 2.5PN normalization stack,
4. and ask whether the static mixed sector really creates a new barrier-softening family, or only renormalizes short-range kernels that were already present.

The main result is strong:

> in the minimal static one-port setting, the mixed bundle does **not** create a new long-range attractive family.
> It produces only quadratic product kernels built from the primitive source profiles.

So the static mixed corridor survives, but it is already much narrower than the generic Stage-001 placeholder made it look.

---

## 1. Frozen input from the moving-throat one-port bundle

Work with the same isotropic one-port reduction already used in the moving-throat normalization chain.
After the stable BdG support mode is integrated out, the effective wall stiffness is

\[
K_*:=K-\frac{C^2}{\varpi^2}.
\]

The one-port Maxwell/mixed block is then controlled by

\[
\Delta:=\Omega_U^2\Omega_W^2-R^2,
\]

\[
Q:=G_U^2\Omega_W^2+2G_UG_WR+G_W^2\Omega_U^2,
\]

\[
P:=\Omega_U^2G_W+RG_U,
\]

and the static conservative wall operator is

\[
D_0
:=
K_* - \frac{Q}{\Delta}.
\]

These are exactly the same bundle quantities that already feed the isotropic 5PN / 2.5PN ledger:

\[
P_0=\frac{N_0}{D_0},
\qquad
N_0=\frac{P^2}{\Delta^2},
\qquad
m_{\hat 0}^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5}.
\]

So the barrier-softening audit is now directly tied to the same reduced bundle data that the normalization program already cares about.

---

## 2. Exact static reduced bundle and admissibility conditions

Take three reduced static coordinates:

- wall/worldtube amplitude `q`,
- brane-like internal gauge coordinate `U`,
- mixed `A_w/F_{\mu w}/J^w` coordinate `W`.

Their static quadratic energy is

\[
V_{\rm stat}(q,U,W;r)
=
\frac12 K_* q^2
+\frac12 \Omega_U^2 U^2
+\frac12 \Omega_W^2 W^2
-RUW
-G_U qU
-G_W qW
-J_q(r)q-J_U(r)U-J_W(r)W.
\]

In matrix form,

\[
V_{\rm stat}
=
\frac12 X^T\mathcal K_{\rm red}X-J^TX,
\qquad
X=(q,U,W)^T,
\qquad
J=(J_q,J_U,J_W)^T,
\]

with

\[
\mathcal K_{\rm red}
=
\begin{pmatrix}
K_* & -G_U & -G_W \\
-G_U & \Omega_U^2 & -R \\
-G_W & -R & \Omega_W^2
\end{pmatrix}.
\]

A direct computation gives the exact determinant identity

\[
\boxed{
\det \mathcal K_{\rm red}=\Delta D_0.
}
\]

So the natural static admissibility conditions are the same ones that already appeared in the normalization program:

\[
\Omega_U^2>0,
\qquad
\Delta>0,
\qquad
D_0>0.
\]

Under these hypotheses the internal `(U,W)` block is positive definite and the full reduced static matrix is positive definite by the Schur-complement criterion.

That matters because it fixes the sign of the induced static same-charge correction.

---

## 3. Exact static susceptibility kernel

Minimizing the quadratic energy gives

\[
X_*(r)=\mathcal K_{\rm red}^{-1}J(r),
\]

and the exact on-shell energy shift is

\[
\boxed{
\delta V_{\rm mix}(r)
=
-\frac12 J(r)^T\mathcal K_{\rm red}^{-1}J(r).
}
\]

Because `\mathcal K_{\rm red}` is positive definite on the admissible branch, the quadratic form is nonnegative and therefore

\[
\boxed{
\delta V_{\rm mix}(r)\le 0.
}
\]

So the first one-port mixed-bundle result is simple but important:

> the static mixed bundle is always attractive or neutral at second order in the external load.

### 3.1 Exact inverse entries

The inverse matrix entries are

\[
\chi_{qq}:=(\mathcal K_{\rm red}^{-1})_{qq}=\frac{1}{D_0},
\]

\[
\chi_{qU}:=(\mathcal K_{\rm red}^{-1})_{qU}=\frac{P_U}{\Delta D_0},
\qquad
P_U:=G_U\Omega_W^2+RG_W,
\]

\[
\chi_{qW}:=(\mathcal K_{\rm red}^{-1})_{qW}=\frac{P}{\Delta D_0},
\]

\[
\chi_{UU}:=(\mathcal K_{\rm red}^{-1})_{UU}
=
\frac{K_*\Omega_W^2-G_W^2}{\Delta D_0},
\]

\[
\chi_{UW}:=(\mathcal K_{\rm red}^{-1})_{UW}
=
\frac{K_*R+G_UG_W}{\Delta D_0},
\]

\[
\chi_{WW}:=(\mathcal K_{\rm red}^{-1})_{WW}
=
\frac{K_*\Omega_U^2-G_U^2}{\Delta D_0}.
\]

So every static same-charge correction is controlled by the same denominator `\Delta D_0`.

### 3.2 Collinear-source theorem

If the external load is collinear in the reduced source space,

\[
J(r)=\mathcal S(r)
\begin{pmatrix}
s_q\\ s_U\\ s_W
\end{pmatrix},
\]

then the induced correction factorizes exactly:

\[
\boxed{
\delta V_{\rm mix}(r)
=
-\frac12\chi_s\,\mathcal S(r)^2,
\qquad
\chi_s=\frac{\mathcal N_s}{\Delta D_0},
}
\]

with

\[
\mathcal N_s
=
\Delta s_q^2
+2P_U s_q s_U
+2P s_q s_W
+(K_*\Omega_W^2-G_W^2)s_U^2
+2(K_*R+G_UG_W)s_U s_W
+(K_*\Omega_U^2-G_U^2)s_W^2.
\]

So the first PDE-constrained mixed kernel is not a free function.
It is the square of the source profile times an exact bundle susceptibility.

That already places a strong restriction on how flexible the barrier-softening corridor can be.

---

## 4. Exact bridge to the 5PN / 2.5PN outgoing-prefactor data

Define the outgoing-load factor

\[
\Lambda:=\frac{P}{\Delta}.
\]

Then the static outgoing-transfer and prefactor quantities are

\[
N_0=\Lambda^2,
\qquad
P_0=\frac{\Lambda^2}{D_0}.
\]

The wall-mixed cross susceptibility becomes

\[
\boxed{
\chi_{qW}=\frac{\Lambda}{D_0}.
}
\]

So the same outgoing-load factor that feeds the 5PN / 2.5PN normalization chain also controls the static wall–mixed barrier response.

There is an especially useful identity:

\[
\chi_{qW}^2=\frac{P_0}{D_0}.
\]

So static same-charge softening and outgoing quadrupole normalization are not independent knobs inside the one-port bundle.

If one tries to enhance the static barrier-softening route by:

- increasing `\Lambda`, or
- softening `D_0`,

one is simultaneously pushing on the same variables that govern the normalization side of the 5PN / 2.5PN story.

That is the first exact bridge from the barrier audit back into the existing bundle ledger.

---

## 5. Product-kernel theorem for the first primitive source families

The next question is whether the static mixed bundle creates a genuinely new radial family.

Use the first two primitive reduced source profiles that are already natural from the frozen files:

- quadrupolar/Coulomb-Hessian drive
  \[
  \mathcal S_Q(x)=\frac{1}{x^3},
  \]
- Yukawa/mixed-sector drive
  \[
  \mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
  \qquad
  \kappa:=\frac{a}{\lambda}.
  \]

Take the reduced source vector to be

\[
J(x)=
\begin{pmatrix}
\beta_Q\,\mathcal S_Q(x)\\
\beta_U\,\mathcal S_Y(x)\\
\beta_W\,\mathcal S_Y(x)
\end{pmatrix}.
\]

Then the exact static correction is

\[
\boxed{
\widetilde{\delta V}_{\rm mix}(x)
=
-\frac12
\left[
\frac{\mathcal C_6}{x^6}
+
2\mathcal C_4\,\frac{e^{-2\kappa x}}{x^4}
+
\mathcal C_2\,\frac{e^{-4\kappa x}}{x^2}
\right],
}
\]

where

\[
\mathcal C_6=\chi_{qq}\,\beta_Q^2,
\]

\[
\mathcal C_4=
\chi_{qU}\,\beta_Q\beta_U
+
\chi_{qW}\,\beta_Q\beta_W,
\]

\[
\mathcal C_2=
\chi_{UU}\,\beta_U^2
+
2\chi_{UW}\,\beta_U\beta_W
+
\chi_{WW}\,\beta_W^2.
\]

This is the key structural result of the stage.

### 5.1 What the theorem says

The minimal static one-port mixed bundle can generate only the three product families

\[
\frac{1}{x^6},
\qquad
\frac{e^{-2\kappa x}}{x^4},
\qquad
\frac{e^{-4\kappa x}}{x^2}.
\]

So it does **not** create a new slower-than-source attractive family.
In particular, it does not create a new long-range attractive term of the form

\[
-\frac{1}{x},
\qquad
-\frac{e^{-2\kappa x}}{x},
\qquad
\text{or anything longer-ranged than the primitive source profiles themselves.}
\]

What it does do is:

1. renormalize the already-known inverse-sixth channel,
2. renormalize the already-known Yukawa-square channel,
3. and add one hybrid cross family
   \[
   \frac{e^{-2\kappa x}}{x^4},
   \]
   which is still shorter-ranged than a Yukawa tail and still much shorter-ranged than Coulomb.

So the static mixed bundle is already much less exotic than the generic Stage-001 placeholder suggested.

---

## 6. Updated reduced barrier audit potential

The natural Stage-002 replacement for the generic mixed term is therefore

\[
\widetilde V_{\rm audit}^{(2)}(x)
=
\frac{1}{x}\Bigl(1+\frac12 e^{-2\kappa x}\Bigr)
-
3\alpha_6\,\mathcal G_6(x)
-
\alpha_2\frac{e^{-4\kappa x}}{x^2}
+
\widetilde{\delta V}_{\rm mix}(x),
\]

with `\widetilde{\delta V}_{\rm mix}` fixed by the one-port susceptibility theorem above.

Equivalently, after collecting like families,

\[
\widetilde V_{\rm audit}^{(2)}(x)
=
\frac{1}{x}\Bigl(1+\frac12 e^{-2\kappa x}\Bigr)
-
\frac{\alpha_{6,\rm eff}}{x^6}
-
\alpha_{4,\rm mix}\frac{e^{-2\kappa x}}{x^4}
-
\alpha_{2,\rm eff}\frac{e^{-4\kappa x}}{x^2}
+
\text{(core-resolved replacement of }x^{-6}\text{ near overlap)}.
\]

So the minimal static mixed bundle does **not** open an entirely new barrier corridor.
It only changes the coefficients of short-range attractive families, with one extra hybrid family in between.

That is a large narrowing of the survival space.

---

## 7. What this kills, and what still survives

### 7.1 What is now strongly disfavored

The following picture is now much harder to sustain inside the present static audit:

> a modest mixed-sector static deformation produces a qualitatively new long-range attractive channel that materially changes same-charge approach before the usual short-range response terms matter.

The one-port bundle does not do that.
Its static correction is a quadratic response built from source products, and those products are all short-range.

### 7.2 What still survives

A static same-charge mixed corridor is **not** dead yet, but it has become much more specific.
It can survive only through one of the following:

1. unusually large coefficient renormalization of the already-known short-range families,
2. a sizable hybrid `e^{-2\kappa x}/x^4` cross family,
3. or a near-instability softening of `D_0` that magnifies the static susceptibility.

But the third option is not a free win, because the same `D_0` appears in the outgoing-prefactor and normalization chain.
So trying to exploit near-instability for barrier softening pushes directly on the same bundle objects that the 5PN / 2.5PN program already constrains.

### 7.3 Practical reading

At the static one-port level, the surviving question is no longer

> can the mixed sector create a magical new attractive law?

It is now

> can the mixed bundle generate **large enough coefficients** on the short-range families it already knows how to produce, while staying on an admissible branch?

That is a much sharper kill test.

---

## 8. Current verdict after Stage 002

Stage 001 already killed the naive Maxwell-only story and showed that pure inverse-sixth attraction cannot erase the barrier by itself.

Stage 002 narrows the static mixed corridor further.

1. The first PDE-constrained static mixed correction is an exact quadratic susceptibility of the one-port bundle.
2. Its denominator is `\Delta D_0`, the same bundle object that already controls the normalization chain.
3. For the first primitive reduced source families, the bundle can only produce
   \[
   x^{-6},\qquad e^{-2\kappa x}/x^4,\qquad e^{-4\kappa x}/x^2.
   \]
4. So the static corridor is no longer about discovering a new long-range attractive law.
   It is about whether the actual moving-throat branch can supply sufficiently large **short-range** mixed-bundle coefficients without violating the same admissibility conditions the 5PN / 2.5PN program already needs.

That means the idea is still alive, but the static version is already looking much more like

- coefficient engineering of short-range response,

than like

- genuine barrier bypass.

---

## 9. Immediate next step

The next honest continuation is now even narrower.

1. Keep the Stage-002 static susceptibility kernel.
2. Test whether the actual branch data can make its coefficients large while keeping
   \[
   \Delta>0,
   \qquad
   D_0>0,
   \qquad
   \Delta_{\rm branch}=0,
   \qquad
   \Delta_{\rm orbit}=0.
   \]
3. If the static corridor stays too weak, move to the next nontrivial possibility:
   a genuinely **time-dependent / non-adiabatic** mixed-sector audit, where dynamic pumping might change the effective kernel class rather than merely its coefficients.

That is the cleanest continuation point now.

same_charge_barrier_audit_stage003_dynamic_mixed_port_kernel.md

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

same_charge_barrier_audit_stage004_resonance_linewidth_tradeoff.md


# Same-Charge Barrier Audit — Stage 004: Resonance/Linewidth Tradeoff, the Dispersive No-Free-Lunch Theorem, and the Linear Survival Window

## 0. Purpose

Stage 003 killed the naive dynamic story but left one narrow escape hatch alive:

> perhaps a **resonant dispersive** mixed-sector window can amplify the already-known short-range attractive families enough to matter before the same pole simply turns into absorptive pumping / leakage.

This stage turns that possibility into an exact line-shape audit.

The job is not to solve the full driven two-throat PDE.
It is much sharper:

1. isolate the local simple-pole normal form of the dynamic one-port mixed bundle,
2. derive the exact relation between the conservative dispersive gain and the absorptive load,
3. specialize that relation to the wall-like pole of the Stage-003 bundle,
4. and state the first honest linear survival criterion.

The result is again restrictive:

> at linear order, the **largest possible conservative barrier reshaping** occurs exactly where the absorptive load is of the same order.
> If one insists on a cleaner low-loss window, the maximum conservative enhancement is suppressed by a simple exact factor.

So after Stage 004, the dynamic corridor is no longer “maybe resonance helps.”
It is:

> only if the actual PDE branch supplies a pole with a sufficiently large residue-to-linewidth ratio can the linear mixed-sector route survive.

---

## 1. Frozen input carried forward

### 1.1 Dynamic one-port bundle from Stage 003

The Stage-003 dynamic reduced correction has the exact collinear-source form

\[
\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\,\chi_s(\omega)\,\mathcal S(x,\omega)^2,
\qquad
\chi_s(\omega)=\frac{\mathcal N_s(\omega)}{\Delta_\Pi(\omega)\,D_\Pi(\omega)}.
\]

For the primitive same-charge source families
\[
\mathcal S_Q(x)=\frac1{x^3},
\qquad
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
\]
the spatial kernel class is already frozen:

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

So the dynamic audit is **not** about inventing a new spatial family.
It is about what the pole structure can do to the coefficients of the existing short-range families.

### 1.2 Why the pole analysis is the right next move

The 4D/plasma ontology keeps the mixed-sector variables

\[
A_w,\qquad F_{\mu w},\qquad J^w
\]

alive beyond the strict far-field brane reduction, and the moving-throat / 5PN stack has already narrowed the remaining weak-axisymmetric grouped payload to the outgoing prefactor side. So once the static and first dynamic no-gos are accepted, the only linear corridor left is a mixed-sector **resonant** one.

---

## 2. Exact simple-pole normal form

Let

\[
F_\Pi(\omega):=\Delta_\Pi(\omega)\,D_\Pi(\omega)
\]

be the full dynamic denominator entering the reduced susceptibility.
Suppose the conservative branch has a simple real pole at \(\omega=\omega_*\):

\[
F_0(\omega_* )=0,
\qquad
F_0'(\omega_*)\neq 0,
\qquad
F_0:=F_\Pi\big|_{\Pi=0}.
\]

Define the outgoing-port sensitivity
\[
Z_* := -\partial_\Pi F_\Pi(\omega)\big|_{(\omega_*,\,\Pi=0)}.
\]

Then to first order in detuning
\[
\delta:=\omega-\omega_*,
\]
and to first order in the passive/outgoing port,
\[
F_\Pi(\omega)
=
F_0'(\omega_*)\,\delta
-
\Pi(\omega_*)\,Z_*
+
O(\delta^2,\Pi\delta,\Pi^2).
\]

If the port is passive,
\[
\Pi(\omega_*) = i\,\Gamma_*,
\qquad
\Gamma_*\ge 0,
\]
then every reduced coefficient near that pole has the universal Breit–Wigner form

\[
\chi_s(\omega)
\approx
\frac{A_*}{\delta-i\gamma_*},
\qquad
A_*:=\frac{\mathcal N_s(\omega_*)}{F_0'(\omega_*)},
\qquad
\gamma_*:=\frac{\Gamma_* Z_*}{|F_0'(\omega_*)|}.
\]

So the whole linear resonance problem collapses to two real scalars:

- the signed residue scale \(A_*\),
- the positive linewidth \(\gamma_*\).

---

## 3. Specialization to the wall-like pole

The Stage-003 one-port bundle is controlled by the reduced wall denominator
\[
D_\Pi(\omega)
=
K_B(\omega) - \frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)}.
\]

The exact Stage-003 derivative identity is
\[
\partial_\Pi D_\Pi(\omega) = -N(\omega),
\]
where
\[
N(\omega)=\frac{P(\omega)^2}{\Delta_\Pi(\omega)^2}.
\]

So at a **wall-like** conservative pole
\[
D_0(\omega_*)=0,
\qquad
\Delta_0(\omega_*)\neq 0,
\]
one has the local normal form

\[
D_\Pi(\omega)
=
D_0'(\omega_*)\,\delta
-
\Pi(\omega_*)\,N_*
+
O(\delta^2,\Pi\delta,\Pi^2),
\qquad
N_*:=N(\omega_*).
\]

Therefore the wall susceptibility itself is

\[
\chi_{qq}(\omega)
=
\frac1{D_\Pi(\omega)}
\approx
\frac{1}{D_0'(\omega_*)}\,
\frac1{\delta-i\gamma_{\rm wall}},
\qquad
\gamma_{\rm wall}:=\frac{\Gamma_*\,N_*}{|D_0'(\omega_*)|}.
\]

This is the first exact dynamic self-limiting statement of the audit:

> the same transfer factor \(N_*\) that helps the outgoing quadrupole bridge also widens the wall pole once the passive/outgoing port is restored.

So a larger transfer strength is not a free win. It simultaneously increases the linewidth that weakens the dispersive amplification.

---

## 4. Exact dispersive/absorptive tradeoff theorem

Now write the universal local line shape as

\[
\chi_*(\omega)=\frac{A_*}{\delta-i\gamma_*},
\qquad
\gamma_*>0.
\]

Rationalizing the denominator gives

\[
\chi_*(\omega)
=
A_*\frac{\delta+i\gamma_*}{\delta^2+\gamma_*^2}.
\]

So the conservative and absorptive pieces are exactly

\[
\Re \chi_*(\omega)=A_*\frac{\delta}{\delta^2+\gamma_*^2},
\qquad
\Im \chi_*(\omega)=A_*\frac{\gamma_*}{\delta^2+\gamma_*^2}.
\]

Introduce the dimensionless detuning ratio

\[
r:=\frac{|\delta|}{\gamma_*}.
\]

Then the line-shape magnitudes collapse to

\[
|\Re\chi_*|
=
\frac{|A_*|}{\gamma_*}\,
\frac{r}{1+r^2},
\qquad
|\Im\chi_*|
=
\frac{|A_*|}{\gamma_*}\,
\frac{1}{1+r^2},
\]
and therefore

\[
\boxed{
\frac{|\Re\chi_*|}{|\Im\chi_*|}=r.
}
\]

So the ratio of useful conservative reshaping to absorptive loading is **nothing but the detuning in linewidth units**.

### 4.1 Exact maximum of the conservative line shape

The dispersive factor
\[
f(r):=\frac{r}{1+r^2}
\]
has derivative
\[
f'(r)=\frac{1-r^2}{(1+r^2)^2}.
\]

So the exact maximum occurs at
\[
r=1,
\]
that is,
\[
|\delta|=\gamma_*.
\]

At that point
\[
\boxed{
\max_r |\Re\chi_*|
=
\frac{|A_*|}{2\gamma_*},
}
\]
and simultaneously
\[
\boxed{
|\Re\chi_*|=|\Im\chi_*|.
}
\]

This is the cleanest no-free-lunch theorem in the linear audit:

> the largest possible conservative dispersive enhancement appears precisely where the absorptive load is of the same size.

### 4.2 Exact low-loss bound

Suppose one demands that the absorptive part be at most a fraction \(\eta\) of the conservative part:

\[
|\Im\chi_*| \le \eta\,|\Re\chi_*|,
\qquad
0<\eta\le 1.
\]

Using \(|\Re|/|\Im|=r\), this is equivalent to
\[
r\ge \frac1\eta.
\]

Since \(f(r)\) decreases for \(r\ge 1\), the largest conservative magnitude allowed by that low-loss condition occurs at the boundary \(r=1/\eta\). Therefore

\[
\boxed{
\sup_{\,|\Im\chi_*|\le \eta|\Re\chi_*|}\,|\Re\chi_*|
=
\frac{|A_*|}{\gamma_*}\,
\frac{\eta}{1+\eta^2}.
}
\]

For small \(\eta\),
\[
\frac{\eta}{1+\eta^2}
=
\eta-\eta^3+O(\eta^5).
\]

So in a genuinely low-loss window, the best linear conservative enhancement scales only **linearly** with the allowed loss fraction.

This is the central Stage-004 theorem.

---

## 5. Barrier language and absorbed-power language

Take one of the already-frozen spatial families \(S_j(x)\) from Stage 003.
Near a simple pole its coefficient has the local form
\[
\chi_j(\omega)\approx \frac{A_j}{\delta-i\gamma_*}.
\]

Then the in-phase barrier reshaping contribution is

\[
U_j^{\rm disp}(x,\omega)
=
-\frac12\,\Re\chi_j(\omega)\,S_j(x)^2,
\]

while the out-of-phase absorbed-power diagnostic is

\[
\overline P_j(x,\omega)
=
-\omega\,\Im\mathfrak V_j(x,\omega)
=
\frac{\omega}{2}\,|\Im\chi_j(\omega)|\,S_j(x)^2.
\]

So the exact tradeoff becomes

\[
\boxed{
\frac{\overline P_j(x,\omega)}{\omega\,|U_j^{\rm disp}(x,\omega)|}
=
\frac{|\Im\chi_j|}{|\Re\chi_j|}
=
\frac1r
=
\frac{\gamma_*}{|\delta|}.
}
\]

This makes the physical reading completely transparent:

- **on resonance** \((\delta=0)\), the conservative reshaping vanishes and the response is purely absorptive;
- **at the dispersive optimum** \((|\delta|=\gamma_*)\), the power load and the conservative barrier term are locked one-to-one;
- **in a low-loss window** \((|\delta|\gg\gamma_*)\), the absorptive fraction is smaller, but the conservative gain is proportionally suppressed.

---

## 6. Quality-factor form of the same theorem

Introduce the local quality factor
\[
Q_*:=\frac{\omega_*}{2\gamma_*}.
\]

Then
\[
\frac{|\delta|}{\omega_*}=\frac{r}{2Q_*}.
\]

So if one imposes the same low-loss condition
\[
|\Im\chi_*| \le \eta |\Re\chi_*|,
\]
the exact detuning requirement becomes

\[
\boxed{
\frac{|\omega-\omega_*|}{\omega_*}
\ge
\frac{1}{2Q_*\eta}.
}
\]

This is useful because it tells us how demanding the “clean dispersive” corridor really is.

- If \(Q_*\) is modest, the required low-loss detuning is not a tiny perturbation of the resonance.
- If \(Q_*\) is very large, the corridor is narrow in absolute frequency, but the residue-to-linewidth scale \(|A_*|/\gamma_*\) can also become large enough to matter.

So the surviving linear question is not “does resonance exist?”
It is:

> does the actual moving-throat branch support a pole with sufficiently large \(Q_*\) and sufficiently large residue scale \(|A_*|\) that the bound above can beat the same-charge barrier threshold?

---

## 7. What survives after Stage 004

Stage 003 already proved that linear monochromatic driving never creates a new spatial family.
Stage 004 now adds the exact resonance theorem.

### 7.1 What is dead

The following picture is now strongly disfavored:

> tune a monochromatic mixed-sector drive to resonance and obtain a large conservative barrier reduction without paying a comparable pumping / leakage price.

That does not happen at linear order.

- On resonance the effect is purely absorptive.
- At the point of maximal conservative gain the absorptive and conservative pieces are equal in magnitude.
- In lower-loss windows the best achievable conservative gain is bounded by
  \[
  \frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2}.
  \]

So there is no linear resonant “free lunch.”

### 7.2 What is still alive

A linear dynamic corridor still survives, but only in the narrow form:

1. one of the already-known short-range families must already be spatially relevant in the tunneling region,
2. the actual PDE branch must supply a pole with sufficiently large \(|A_*|/\gamma_*|\),
3. and that branch must remain admissible under the same 5PN / outgoing constraints already frozen earlier.

So after Stage 004, the dynamic same-charge route lives or dies on a **residue-to-linewidth** question, not on a generic resonance slogan.

---

## 8. Best current verdict after Stage 004

The idea is still alive, but it has narrowed again.

- Stage 001 killed the naive Coulomb/KK story.
- Stage 002 killed the hope for a brand-new static mixed-sector spatial law.
- Stage 003 killed the linear outgoing-phase shortcut.
- Stage 004 now kills the naive “resonance fixes it” shortcut.

What survives is only this:

> a high-quality mixed-sector pole could still amplify one of the already-existing short-range attractive families enough to matter, but only if its residue-to-linewidth ratio is large enough to overcome the exact low-loss tradeoff bound.

That is a real theorem gate, not a slogan.

---

## 9. Immediate next step

The continuation point is now very sharp.

1. Keep the exact Stage-003 dynamic one-port bundle.
2. Choose an explicit primitive branch family for
   \[
   K,\ M,\ C,\ \varpi,\ \Omega_U,\ \Omega_W,\ G_U,\ G_W,\ R.
   \]
3. Compute the first actual pole locations \(\omega_*\), residues \(A_*\), and linewidths \(\gamma_*\).
4. Compare
   \[
   \frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2}
   \]
   against the barrier-softening size demanded by the reduced same-charge problem.

That is the smallest honest continuation point after the resonance/linewidth audit.

same_charge_barrier_audit_stage005_concrete_branch_residue_linewidth_test.md

# Same-Charge Barrier Audit — Stage 005: Concrete Finite-Throat Primitive Branch, Pole Census, and the Residue/Linewidth Survival Test

## 0. Purpose

Stage 004 reduced the linear dynamic corridor to one exact question:

> can an admissible pole of the mixed one-port bundle produce a large enough residue-to-linewidth ratio to matter in the tunneling region before the same pole is dominated by absorptive loading?

The next honest move is therefore not another generic resonance slogan.
It is to pick one **explicit primitive branch family**, compute the first actual poles, and write the survival test in direct branch data.

That is what this stage does.

The main outputs are:

1. one concrete finite-throat primitive branch built from the lowest N/N wall/U profile and the lowest D/N half-wave support/W profile,
2. the exact quartic pole polynomial of the conservative one-port bundle,
3. the exact simple-pole cancellation
   \[
   \mathcal R_{Q,*}:=\frac{|A_{Q,*}|}{\gamma_*}
   =
   \frac{27c_s^5}{a^5\omega_*^5 N_*(\omega_*)},
   \]
   for the pure quadrupolar same-charge family,
4. the exact low-loss survival inequality against a required local barrier reduction,
5. one explicit numerical pole census,
6. and the first clean evidence of a **static/dynamic tension**: strengthening the static outgoing prefactor makes the linear dynamic residue/linewidth corridor worse on the same branch family.

So after Stage 005, the idea is still alive, but it has become much more specific.

---

## 1. Frozen input carried forward

The parent 4D/plasma ontology keeps the mixed channels
\[
A_w,\qquad F_{\mu w},\qquad J^w
\]
alive outside the strict far-field brane reduction, so the linear dynamic same-charge corridor still has to live in the mixed sector rather than in pure brane Maxwell shaping. The 5PN / moving-throat side independently says that the explicit isotropic finite-throat overlap model is already the correct place to package the grouped conservative/prefactor data, and that the surviving weak-axisymmetric grouped scalar is the outgoing-prefactor slope
\[
\Xi_1=\frac{P_1}{P_0}.
\]
So the correct next move is a concrete pole test on that same mixed/outgoing bundle. fileciteturn21file2turn24file4turn23file18

---

## 2. Concrete primitive finite-throat branch

Fix the finite throat interval
\[
s\in[0,L].
\]
Use the lowest N/N zero mode for the wall and the brane-like internal coordinate,
\[
u_0(s)=\frac1{\sqrt L},
\]
and the lowest D/N half-wave for the trapped support and mixed coordinate,
\[
f_0(s)=\sqrt{\frac2L}\sin\frac{\pi s}{2L}.
\]
Then the exact overlap constant is
\[
\kappa:=\int_0^L u_0(s)f_0(s)\,ds=\frac{2\sqrt2}{\pi}.
\]
On this branch the overlap-renormalized one-port couplings are
\[
C=\kappa\lambda_B,
\qquad
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R.
\]

The dynamic one-port wall/BdG/Maxwell/mixed bundle is therefore
\[
K_B(\omega)=K-M\omega^2-\frac{C^2}{\varpi^2-\omega^2},
\]
\[
A(\omega)=\Omega_U^2-\omega^2,
\qquad
W(\omega)=\Omega_W^2-\omega^2,
\]
\[
\Delta(\omega)=A(\omega)W(\omega)-R^2,
\]
\[
Q(\omega)=G_U^2W(\omega)+2G_UG_WR+G_W^2A(\omega),
\]
\[
D(\omega)=K_B(\omega)-\frac{Q(\omega)}{\Delta(\omega)}.
\]

The static admissibility conditions remain the same as earlier:
\[
\Delta_0>0,
\qquad
D_0>0,
\]
with
\[
\Delta_0=\Omega_U^2\Omega_W^2-R^2,
\qquad
D_0=K-\frac{C^2}{\varpi^2}-\frac{Q_0}{\Delta_0}.
\]

---

## 3. Exact quartic pole polynomial

Writing
\[
y=\omega^2,
\]
the conservative pole condition is exactly the quartic equation
\[
F(y)=0,
\]
with
\[
F(y)=
\Bigl((K-My)(\varpi^2-y)-C^2\Bigr)
\Bigl((\Omega_U^2-y)(\Omega_W^2-y)-R^2\Bigr)
-
(\varpi^2-y)
\Bigl(G_U^2(\Omega_W^2-y)+2G_UG_WR+G_W^2(\Omega_U^2-y)\Bigr).
\]

So the primitive finite-throat one-port branch has exactly four conservative poles in the generic admissible case.

This is already a useful narrowing:
we are no longer talking about a vague resonance landscape, but about one explicit quartic pole census.

---

## 4. Exact residue/linewidth cancellation

For a simple conservative pole \(\omega_*\) with \(\Delta(\omega_*)\neq 0\), attach the passive outgoing quadrupole port
\[
\Pi_{\rm out}(\omega)=i\Gamma_5\omega^5,
\qquad
\Gamma_5=\frac{a^5}{27c_s^5}.
\]

The exact transfer factor is
\[
N(\omega)=\frac{P(\omega)^2}{\Delta(\omega)^2},
\qquad
P(\omega)=A(\omega)G_W+RG_U.
\]

For the **wall susceptibility**
\[
\chi_{qq}(\omega)=\frac1{D(\omega)},
\]
the simple-pole residue and linewidth satisfy
\[
|A_{qq,*}|=\frac1{|D'(\omega_*)|},
\qquad
\gamma_*=rac{\Gamma_5\omega_*^5N(\omega_*)}{|D'(\omega_*)|},
\]
so the derivative cancels out exactly:
\[
\boxed{
\mathcal R_{qq,*}:=\frac{|A_{qq,*}|}{\gamma_*}
=
\frac{1}{\Gamma_5\omega_*^5N(\omega_*)}
=
\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)}.
}
\]

For the **pure quadrupolar same-charge family**
\[
S_Q(x)=\frac1{x^3},
\]
Stage 002/003 give
\[
\mathfrak V_Q(x,\omega)=-\frac12\chi_{qq}(\omega)\frac1{x^6},
\]
so the same exact ratio controls the useful conservative coefficient:
\[
\boxed{
\mathcal R_{Q,*}=\mathcal R_{qq,*}
=
\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)}.
}
\]

This is the first truly branch-level simplification of the dynamic same-charge problem.
The residue-to-linewidth figure is controlled only by

1. the pole frequency \(\omega_*\), and
2. the exact outgoing transfer factor \(N(\omega_*)\).

Everything else cancels out.

---

## 5. Exact low-loss survival inequality

Let \(\Delta V_{\rm req}(x)\) be the local barrier reduction required at radius \(x\) and impose the same low-loss condition as in Stage 004,
\[
|\Im\chi|\le \eta |\Re\chi|,
\qquad 0<\eta\le 1.
\]
Then the maximum conservative line shape allowed in that window is
\[
|\Re\chi|_{\max,\eta}
=
\frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2}.
\]
For the pure quadrupolar family this gives the exact local survival criterion
\[
\frac12\,\mathcal R_{Q,*}\,\frac{\eta}{1+\eta^2}\,\frac1{x^6}
\ge
\Delta V_{\rm req}(x).
\]
Equivalently,
\[
\boxed{
\mathcal R_{Q,*}
\ge
2\,\Delta V_{\rm req}(x)
\frac{1+\eta^2}{\eta}
\,x^6.
}
\]
Substituting the exact pole ratio gives the branch-level theorem gate
\[
\boxed{
\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)}
\ge
2\,\Delta V_{\rm req}(x)
\frac{1+\eta^2}{\eta}
\,x^6.
}
\]
So the linear dynamic corridor lives or dies on one explicit inequality in the primitive branch data.

---

## 6. Explicit numerical primitive slice

Take the admissible sample branch
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,K,M)
=(0.5,\,0.3,\,0.4,\,0.25,\,1.0,\,1.4,\,2.0,\,3.0,\,1.0)
\]
with \(a=c_s=1\).

The overlap-renormalized couplings are
\[
C\approx 0.450158158078553,
\quad
G_U=0.3,
\quad
G_W\approx 0.360126526462843,
\quad
R\approx 0.225079079039277.
\]
The static branch data are
\[
\Delta_0\approx 1.90933940817883,
\qquad
D_0\approx 2.76355510933127,
\qquad
N_0\approx 0.0501661980249591,
\qquad
P_0\approx 0.0181527764203329.
\]
So the sample slice is statically admissible.

### 6.1 Uncoupled roots

The uncoupled wall/BdG roots are
\[
\omega_{\rm wall/BdG}^{(0)}\approx 1.68143182591478,
\qquad
2.04274007519334,
\]
and the uncoupled internal U/W roots are
\[
\omega_{UW}^{(0)}\approx 0.974601723746314,
\qquad
1.41779810977117.
\]

### 6.2 Full conservative pole census

The full quartic pole census is
\[
\omega_*\approx 0.938272741746753 \quad (\text{internal-like}),
\]
\[
\omega_*\approx 1.39141087653805 \quad (\text{internal-like}),
\]
\[
\omega_*\approx 1.72045371048003 \quad (\text{wall-like}),
\]
\[
\omega_*\approx 2.04539948783659 \quad (\text{wall-like}).
\]

The corresponding pure-Q residue/linewidth figures are
\[
\mathcal R_{Q,*}\approx 18.7069287828307,
\quad
0.380740659074003,
\quad
16.0250330226177,
\quad
32.0025481088465.
\]

So on this concrete slice the strongest candidate is the **upper wall-like pole**.

---

## 7. Comparison to the Stage-001 barrier benchmark

Use the Stage-001 illustrative reduced barrier benchmark at \(x=1\):
\[
V_{\rm known}(1)\approx 1.181909222592,
\qquad
\epsilon=0.1,
\qquad
\Delta V_{\rm req}(1)=V_{\rm known}(1)-\epsilon\approx 1.081909222592.
\]
Since \(S_Q(1)^2=1\), the exact low-loss thresholds are
\[
\mathcal R_{Q,*}^{\rm req}(\eta=0.1)
\approx 21.8545662963584,
\]
\[
\mathcal R_{Q,*}^{\rm req}(\eta=0.3)
\approx 7.86187368416853.
\]

### 7.1 Sample-slice verdict

For the lower wall-like pole,
\[
\mathcal R_{Q,*}\approx 16.0250,
\]
so it

- **fails** the \(10\%\)-loss benchmark,
- but **passes** the \(30\%\)-loss benchmark.

For the upper wall-like pole,
\[
\mathcal R_{Q,*}\approx 32.0025,
\]
so it passes both.

So the concrete primitive slice does **not** kill the dynamic corridor.
But it also does **not** make it generic.
Only some poles on the same branch are actually strong enough.

---

## 8. First static/dynamic tension

Now scan the outgoing-leg coupling \(\lambda_W\) while keeping all other sample parameters fixed.
The resulting static prefactor \(P_0\) and upper-wall residue/linewidth figure \(\mathcal R_{Q,*}\) are:

| \(\lambda_W\) | \(P_0\) | \(D_0\) | upper wall \(\omega_*\) | upper wall \(\mathcal R_{Q,*}\) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.005947405318 | 2.827234421584 | 2.044022723028 | 145.483858657863 |
| 0.4 | 0.018152776420 | 2.763555109331 | 2.045399487837 | 32.002548108846 |
| 0.6 | 0.038000163140 | 2.665913497210 | 2.047932775068 | 13.688535635681 |
| 0.8 | 0.067170782681 | 2.534309585220 | 2.051906688892 | 7.580971267466 |
| 1.0 | 0.108473308110 | 2.368743373361 | 2.057783390355 | 4.827389255647 |

So on this explicit family:

- stronger outgoing coupling raises the static prefactor \(P_0\),
- but the same move lowers the dynamic residue/linewidth figure \(\mathcal R_{Q,*}\).

That is the first clean evidence that the **static outgoing-normalization corridor** and the **linear dynamic low-loss corridor** are in real tension on the same branch family.

This is exactly the kind of tradeoff Stage 004 suggested but could not yet show concretely.

---

## 9. What Stage 005 changes

Stage 005 does not prove that the same-charge idea works.
But it does change the status materially.

Before this stage, the dynamic corridor lived or died on a generic slogan:
“find a pole with a large enough residue-to-linewidth ratio.”

After this stage, we have a much sharper statement:

1. the primitive finite-throat branch gives an exact quartic pole census,
2. the pure-Q same-charge family has the exact branch-level figure
   \[
   \mathcal R_{Q,*}=\frac{27c_s^5}{a^5\omega_*^5N(\omega_*)},
   \]
3. the survival test is one explicit inequality in \(\omega_*\) and \(N(\omega_*)\),
4. and on a concrete admissible slice the corridor is **non-empty**, but only for part of the pole spectrum and only with a real static/dynamic tradeoff.

So the idea survives this stage, but narrowly.

---

## 10. Immediate next step

The next clean move is now even sharper than before:

1. keep the exact primitive finite-throat branch,
2. replace the sample numerical couplings by a PDE-grounded overlap extractor,
3. evaluate \(N(\omega_*)\) and the pole census on the actual moving-throat branch,
4. and compare the resulting \(\mathcal R_{Q,*}\) against the same low-loss survival threshold.

That is the first point where the surviving same-charge dynamic corridor stops being a generic mixed-sector hope and becomes a concrete branch test.

same_charge_barrier_audit_stage006_5pn_target_surface_dynamic_window.md

# Same-Charge Barrier Audit — Stage 006: 5PN Isotropic Target Surface, Primitive-Branch Compatibility, and the Dynamic Survival Window

## 0. Purpose

Stage 005 showed that the linear same-charge corridor does **not** immediately die on the explicit primitive finite-throat branch. But that result was still too loose, because the primitive slice had not yet been forced onto the exact isotropic target surface already isolated by the 5PN / 2.5PN / 4PN moving-throat endgame.

So the next honest question is:

> if the same primitive one-port branch is required to satisfy the exact isotropic one-pole and outgoing-normalization conditions, does the dynamic same-charge corridor survive or die?

That is what this stage answers.

The main outputs are:

1. the exact symbolic compatibility equation between the isotropic one-pole condition and the isotropic outgoing-normalization condition,
2. its specialization to the explicit primitive finite-throat one-port family,
3. one concrete compatibility branch on the Stage-005 sample slice,
4. the resulting compatibility-branch pole census,
5. and a finite **dynamic survival window** in the branch-compatible target parameter.

So after Stage 006, the problem is no longer merely “find a good pole.” It is:

> can the same branch support a good pole **while also lying on the 5PN-compatible isotropic target surface**?

---

## 1. Frozen input carried forward

### 1.1 Primitive finite-throat one-port branch from Stage 005

Keep the same explicit finite-throat branch:

- lowest N/N zero mode for the wall and brane-like internal coordinate,
- lowest D/N half-wave for the trapped support and mixed coordinate,
- overlap constant
  \[
  \kappa = \frac{2\sqrt2}{\pi}.
  \]

With reduced couplings
\[
C=\kappa\lambda_B,
\qquad
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R,
\]
we still have
\[
\Delta = \Omega_U^2\Omega_W^2-R^2,
\qquad
Q = G_U^2\Omega_W^2 + 2G_UG_WR + G_W^2\Omega_U^2,
\qquad
P = \Omega_U^2G_W + RG_U.
\]

The primitive bundle moments are
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
where
\[
S_2=\Omega_U^2+\Omega_W^2,
\qquad
H=G_U^2+G_W^2,
\]
and
\[
N_0=\frac{P^2}{\Delta^2}.
\]

### 1.2 Exact isotropic 5PN / 2.5PN target surface

On the isotropic one-port bundle,
\[
D_0 = K-B_0-Z_0,
\qquad
D_2 = -(M+B_2+Z_2),
\qquad
D_4 = -(B_4+Z_4),
\]
with normalized conservative response
\[
u_2 = -\frac{D_2}{D_0},
\qquad
u_4 = \frac{D_2^2-D_0D_4}{D_0^2},
\]
and outgoing prefactor
\[
P_0 = \frac{N_0}{D_0}.
\]

The exact isotropic one-pole condition is
\[
\nu_4 = 4\nu_2^2
\iff
D_0(B_4+Z_4) = 3(M+B_2+Z_2)^2.
\]
The isotropic outgoing-normalization condition is
\[
P_0 = P_{0,\mathrm{target}},
\]
where for the fully calibrated moving-throat branch
\[
P_{0,\mathrm{target}} = \frac{54Gc_s^5}{5a^5c^5\,\hat m_0^{\,2}}.
\]

Stage 006 keeps this target symbolic as \(P_{0,\mathrm{target}}\), because on the primitive reduced branch the important question is compatibility first.

---

## 2. Exact compatibility equation

The one-pole condition solves for the wall stiffness as
\[
K_{\mathrm{pole}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4} + B_0 + Z_0.
\]
The outgoing-normalization condition solves for the same wall stiffness as
\[
K_{\mathrm{norm}}
=
\frac{N_0}{P_{0,\mathrm{target}}} + B_0 + Z_0.
\]
So simultaneous isotropic one-pole success and isotropic normalization success require
\[
K_{\mathrm{pole}} = K_{\mathrm{norm}},
\]
which is equivalent to the exact compatibility equation
\[
\boxed{
\frac{N_0}{P_{0,\mathrm{target}}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
}
\]
Equivalently, the primitive branch itself induces the unique branch-compatible target
\[
\boxed{
P_{0,\mathrm{target,compat}} = \frac{N_0(B_4+Z_4)}{3(M+B_2+Z_2)^2}.
}
\]
This is the first exact place where the same-charge branch is forced to talk directly to the 5PN isotropic surface.

Two points are worth emphasizing.

First, this is **not** a generic fit condition. It is an exact consequence of trying to satisfy the same isotropic target surface from two sides.

Second, the compatibility equation does **not** determine every coupling. It tells us whether the primitive branch wants a normalization target that is even compatible with its own conservative one-pole structure.

---

## 3. Primitive specialization of the compatibility equation

Substituting the explicit primitive one-port data gives
\[
P_{0,\mathrm{target,compat}}
=
\frac{\dfrac{P^2}{\Delta^2}\left(\dfrac{C^2}{\varpi^6}+Z_4\right)}{3\left(M+\dfrac{C^2}{\varpi^4}+Z_2\right)^2},
\]
or equivalently
\[
\boxed{
\frac{P^2/\Delta^2}{P_{0,\mathrm{target}}}
=
\frac{3\left(M+\dfrac{C^2}{\varpi^4}+Z_2\right)^2}{\dfrac{C^2}{\varpi^6}+Z_4}.
}
\]
So on the primitive family the isotropic 5PN-compatible surface is already a single explicit algebraic relation in the radial/axial couplings and frequencies.

That is a substantial tightening compared with Stage 005.

---

## 4. Concrete sample compatibility branch

Now specialize to the same Stage-005 sample values
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
(0.5,0.3,0.4,0.25,1.0,1.4,2.0,1.0),
\]
with \(a=c_s=1\).

The overlap-renormalized primitive data are
\[
C \approx 0.450158158078553,
\qquad
G_U = 0.3,
\qquad
G_W \approx 0.360126526462843,
\qquad
R \approx 0.225079079039277.
\]
The static bundle quantities are
\[
\Delta \approx 1.90933940817883,
\qquad
Q \approx 0.354725283210515,
\qquad
P \approx 0.427650250174625,
\]
\[
B_0 \approx 0.0506605918211689,
\quad
B_2 \approx 0.0126651479552922,
\quad
B_4 \approx 0.00316628698882306,
\]
\[
Z_0 \approx 0.185784298847558,
\quad
Z_2 \approx 0.172955320626603,
\quad
Z_4 \approx 0.170825285860668,
\]
\[
N_0 \approx 0.0501661980249591.
\]

The exact compatibility target on this primitive slice is
\[
\boxed{
P_{0,\mathrm{target,compat}}
\approx 0.00206979231806289.
}
\]
The corresponding compatibility wall stiffness is
\[
\boxed{
K_{\mathrm{compat}}
\approx 24.4737548792910.
}
\]
So the compatibility branch is much stiffer than the loose Stage-005 sample branch. Its compatible static denominator is
\[
D_{0,\mathrm{compat}} = K_{\mathrm{compat}}-B_0-Z_0
\approx 24.2373099886222.
\]

This is already informative.

The Stage-005 sample branch had
\[
P_0 \approx 0.0181527764203329,
\]
while the same primitive family, when forced onto the exact isotropic one-pole/normalization compatibility surface, wants
\[
P_{0,\mathrm{target,compat}} \approx 0.00207.
\]
So the 5PN-compatible branch lives at a much lower static prefactor and a much higher wall stiffness on this particular primitive slice.

---

## 5. Pole census on the compatibility branch

Using \(K=K_{\mathrm{compat}}\), the conservative pole census is
\[
\omega_* \approx 0.971575315129468 \quad (\text{internal-like}),
\]
\[
\omega_* \approx 1.41651290122561 \quad (\text{internal-like}),
\]
\[
\omega_* \approx 1.99753567893361 \quad (\text{wall-like}),
\]
\[
\omega_* \approx 4.94905432364313 \quad (\text{wall-like}).
\]

The exact pure-quadrupolar residue/linewidth figure remains
\[
\mathcal R_{Q,*} = \frac{27c_s^5}{a^5\omega_*^5 N(\omega_*)}.
\]
On the compatibility branch the four values are
\[
\mathcal R_{Q,*} \approx 0.159888393135835 \quad (\text{internal-like}),
\]
\[
\mathcal R_{Q,*} \approx 0.000806281535937178 \quad (\text{internal-like}),
\]
\[
\mathcal R_{Q,*} \approx 30.1999075602499 \quad (\text{wall-like}),
\]
\[
\mathcal R_{Q,*} \approx 36.1711864832695 \quad (\text{wall-like}).
\]
So on this concrete 5PN-compatible branch the dynamic same-charge corridor is **not** carried by the internal poles at all. It is carried entirely by the wall-like poles.

This is already a nontrivial structural simplification.

---

## 6. Dynamic survival window on the compatibility surface

Carry forward the same illustrative local barrier benchmark from Stage 005 at \(x=1\):
\[
V_{\mathrm{known}}(1) \approx 1.181909222592,
\qquad
\epsilon = 0.1,
\qquad
\Delta V_{\mathrm{req}}(1) \approx 1.081909222592.
\]
Then the required residue/linewidth thresholds are
\[
\mathcal R_{Q,*}^{\mathrm{req}}(\eta=0.1)
\approx 21.8545662963584,
\]
\[
\mathcal R_{Q,*}^{\mathrm{req}}(\eta=0.3)
\approx 7.86187368416853.
\]

Now scan \(\lambda_W\) **along the exact compatibility surface**, i.e. always resetting \(K\) to \(K_{\mathrm{compat}}(\lambda_W)\). The resulting branch-compatible target and wall-like residue/linewidth figures are:

| \(\lambda_W\) | \(P_{0,\mathrm{target,compat}}\) | \(K_{\mathrm{compat}}\) | lower wall \(\mathcal R_Q\) | upper wall \(\mathcal R_Q\) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.000576970879843 | 29.3158464872314 | 138.814136942081 | 137.502546600713 |
| 0.4 | 0.002069792318063 | 24.4737548792910 | 30.1999075602499 | 36.1711864832695 |
| 0.6 | 0.004865681200486 | 21.1544287401845 | 12.8348600273988 | 16.7575510327116 |
| 0.8 | 0.009169913681573 | 19.0298300900561 | 7.06074242207991 | 9.69035785242054 |
| 1.0 | 0.014981190324091 | 17.7824591822917 | 4.45922850098679 | 6.30111094469551 |

This scan shows a clean monotonic tradeoff on the explicit compatibility family:

- increasing \(\lambda_W\) raises the branch-compatible static target \(P_{0,\mathrm{target,compat}}\),
- the same move lowers the required wall stiffness \(K_{\mathrm{compat}}\),
- and both wall-like dynamic figures \(\mathcal R_Q\) fall monotonically.

So the static/dynamic tension from Stage 005 survives **even after** the branch is forced onto the exact 5PN isotropic target surface.

### 6.1 Finite survival windows

At the stricter \(10\%\)-loss benchmark,

- the **lower wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.00283133168555932,
  \]
- the **upper wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.00359651058968466.
  \]

At the looser \(30\%\)-loss benchmark,

- the **lower wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.00817339430971383,
  \]
- the **upper wall** pole survives only up to
  \[
  P_{0,\mathrm{target,compat}} \lesssim 0.0116633929790174.
  \]

So the dynamic corridor is not generically open on the explicit primitive family. It survives only inside a finite interval of the same branch-compatible target that the isotropic 5PN surface itself wants.

That is the sharpest same-charge compatibility statement reached so far.

---

## 7. What Stage 006 changes

Stage 006 does **not** prove the same-charge idea works.
But it does materially tighten the status.

Before this stage, the dynamic corridor lived on a primitive branch that had not yet been forced to satisfy the exact isotropic 5PN surface.

After this stage, we know four much stronger things.

### 7.1 The primitive branch can be put on the exact isotropic 5PN surface

The isotropic one-pole and normalization conditions do not conflict abstractly. They reduce to one exact compatibility equation.

### 7.2 The dynamic corridor is not killed automatically by that calibration

On the concrete sample slice, once the branch is moved to the compatibility wall stiffness, both wall-like poles still clear the stricter \(10\%\)-loss benchmark.

### 7.3 The same branch develops a genuine target window

The dynamic corridor survives only for a finite range of the branch-compatible static target. So the eventual PDE-selected normalization cannot be arbitrarily large on this primitive family if the dynamic same-charge corridor is to remain alive.

### 7.4 The wall-like poles are now the only relevant survivors

The internal poles are dynamically irrelevant on the compatibility branch. So the surviving same-charge corridor is a wall-like corridor, not a generic mixed-pole corridor.

---

## 8. Best current verdict after Stage 006

The idea is still alive.
But it is now alive in a much narrower form:

> an explicit primitive finite-throat branch can satisfy the exact isotropic 5PN compatibility surface **and** still retain a wall-like dynamic same-charge corridor, but only inside a finite branch-compatible normalization window.

So the next honest move is no longer another generic resonance scan.
It is:

1. extract the actual branch-compatible normalization target from the moving-throat PDE,
2. compare it against the finite survival window above,
3. and see whether the real branch lands inside or outside that window.

That is the next clean kill test.

same_charge_barrier_audit_stage007_pde_branch_packet_window_compiler.md

# Same-Charge Barrier Audit — Stage 007: PDE Branch-Packet Compiler, Weak-Axisymmetric Ceiling Transport, and the First Actual-Branch Kill Test

## 0. Purpose

Stage 006 forced the primitive finite-throat one-port branch onto the exact isotropic `5`PN target surface and found that the dynamic same-charge corridor survives only inside a **finite normalization window**.

But that still left one gap:

> how do we test the **actual** moving-throat / `5`PN branch, rather than the primitive family, against that survival window?

The next honest step is therefore not another primitive scan. It is to compile the real PDE-selected branch data into the same variables the Stage-006 window cares about.

That is what this stage does.

The main outputs are:

1. the exact compiler from the final `5`PN branch packet to the lane prefactors
   \[
   P_0^{(20)},\qquad P_0^{(21)},\qquad P_0^{(22)},
   \]
2. the exact isotropic window test in terms of the normalization defect `\Delta_{\rm norm}`,
3. the exact weak-axisymmetric transported ceiling test in terms of the grouped prefactor defects `a_{P_0}, b_{P_0}`,
4. the axisymmetric specialization
   \[
   b_{P_0}=3a_{P_0}
   \]
   and the equivalent one-scalar amplitude law
   \[
   \Xi_1=\frac{P_1}{P_0},
   \]
5. and the explicit anisotropy headroom left at the Stage-006 compatibility point.

So after Stage 007, the problem is no longer

> “extract the normalization target from the PDE somehow.”

It is now

> “does the actual branch packet land inside one explicit finite corridor in `\Delta_{\rm norm}` and weak-axisymmetric prefactor slope?”

---

## 1. Frozen input carried forward

### 1.1 Exact final branch packet from the `5`PN endgame

The exact reduced branch verdict packet is
\[
\Delta_{\rm branch}
=
(a_2,b_2,a_4,b_4,a_{P_0},b_{P_0},\Delta_{\rm pole},\Delta_{\rm norm}).
\]
Its normalization slot is
\[
\Delta_{\rm norm}
=
\hat m_0^{\,2}\,\bar P_0-
\frac{54Gc_s^5}{5a^5c^5}.
\]
So the actual isotropic mean prefactor seen by the same-charge window test is already fixed exactly by
\[
\boxed{
\bar P_0=
\frac{\Delta_{\rm norm}+\dfrac{54Gc_s^5}{5a^5c^5}}{\hat m_0^{\,2}}.
}
\]

### 1.2 Exact grouped inverse map for the lane prefactors

The grouped prefactor anomalies compile back to the three lane prefactors by
\[
\boxed{P_0^{(20)}=\bar P_0+4a_{P_0},}
\]
\[
\boxed{P_0^{(21)}=\bar P_0-a_{P_0}+b_{P_0},}
\]
\[
\boxed{P_0^{(22)}=\bar P_0-a_{P_0}-b_{P_0}.}
\]
So the actual moving-throat branch packet already tells us exactly which lane-wise static prefactors have to be compared against the transported same-charge window.

### 1.3 Finite survival ceilings carried from Stage 006

Stage 006 found four useful primitive-family ceilings.

At the stricter `10%`-loss benchmark:

- **both wall-like poles survive** up to
  \[
  P_{\rm both}^{(10)}
  =0.0028313316855593175,
  \]
- a **nonempty wall-like corridor** survives up to
  \[
  P_{\rm one}^{(10)}
  =0.0035965105896846573.
  \]

At the looser `30%`-loss benchmark:

- **both wall-like poles survive** up to
  \[
  P_{\rm both}^{(30)}
  =0.00817339430971383,
  \]
- a **nonempty wall-like corridor** survives up to
  \[
  P_{\rm one}^{(30)}
  =0.0116633929790174.
  \]

These are still transported primitive-family ceilings, not final full-PDE theorems. But they are the first exact dynamic windows we have, so they are the right actual-branch kill test to carry forward.

---

## 2. Exact compiler from the actual branch packet to the transported window test

The branch packet now enters the same-charge audit in the most direct possible way.

For any chosen ceiling `P_\mathrm{crit}`, define the three actual lane prefactors
\[
P_{20}=\bar P_0+4a_{P_0},
\qquad
P_{21}=\bar P_0-a_{P_0}+b_{P_0},
\qquad
P_{22}=\bar P_0-a_{P_0}-b_{P_0}.
\]
Then the strongest transported sufficient condition that **all grouped lanes** stay inside that primitive-family window is
\[
\boxed{
\max\{P_{20},P_{21},P_{22}\}
\le P_{\rm crit}.
}
\]
Equivalently, in packet variables,
\[
\boxed{
\max\Bigl\{
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}+4a_{P_0},
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}-a_{P_0}+b_{P_0},
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}-a_{P_0}-b_{P_0}
\Bigr\}
\le P_{\rm crit},
}
\]
where
\[
T_{\rm quad}:=\frac{54Gc_s^5}{5a^5c^5}.
\]

So the Stage-006 window has now been converted into a direct inequality on the actual branch packet.

---

## 3. Exact isotropic kill test in terms of `\Delta_{\rm norm}`

If the actual branch is exactly isotropic at the prefactor level,
\[
a_{P_0}=b_{P_0}=0,
\]
then every lane sees the same static prefactor
\[
P_{20}=P_{21}=P_{22}=\bar P_0.
\]
So the transported isotropic ceiling test is simply
\[
\boxed{\bar P_0\le P_{\rm crit}.}
\]
Using the exact normalization compiler,
\[
\boxed{
\Delta_{\rm norm}
\le
\hat m_0^{\,2}P_{\rm crit}-T_{\rm quad}.
}
\]
This is the first actual-branch same-charge kill test written directly in the endgame residual language.

### 3.1 Exact calibrated-branch lower bound on `\hat m_0`

If the real branch already hits the universal quadrupole normalization exactly,
\[
\Delta_{\rm norm}=0,
\]
then isotropic survival at ceiling `P_{\rm crit}` is equivalent to
\[
\boxed{
\hat m_0^{\,2}
\ge
\frac{T_{\rm quad}}{P_{\rm crit}}.
}
\]
So the same-charge corridor now places a direct lower bound on the source-map factor of the actual calibrated branch.

That is a much sharper statement than anything before Stage 006.

---

## 4. Weak-axisymmetric specialization from the exact grouped signature

The later moving-throat grouped notes collapse the weak-axisymmetric outgoing slippage bundle to one scalar amplitude with grouped signature
\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1,
\]
and identify that amplitude with the physical outgoing-prefactor slope
\[
\boxed{\Xi_1=\frac{P_1}{P_0}.}
\]
So the weak-axisymmetric prefactor lanes take the exact first-order form
\[
\boxed{P_A=\bar P_0\bigl(1+\epsilon\lambda_A\Xi_1\bigr).}
\]
Explicitly,
\[
P_{20}=\bar P_0(1+\epsilon\Xi_1),
\qquad
P_{21}=\bar P_0\Bigl(1+\frac12\epsilon\Xi_1\Bigr),
\qquad
P_{22}=\bar P_0(1-\epsilon\Xi_1).
\]

The grouped trace/anomaly compiler then gives
\[
\boxed{a_{P_0}=\frac{\epsilon\bar P_0\Xi_1}{4},}
\qquad
\boxed{b_{P_0}=\frac{3\epsilon\bar P_0\Xi_1}{4}.}
\]
So the exact weak-axisymmetric branch law is
\[
\boxed{b_{P_0}=3a_{P_0}.}
\]

This is important because it means the actual same-charge window is not generic in the full `(a_{P_0},b_{P_0})` plane. On the weak-axisymmetric branch it collapses to a one-dimensional axisymmetric line.

---

## 5. Exact transported ceiling test in terms of `\Xi_1`

On the axisymmetric weak-anisotropy line,
\[
P_{20}=\bar P_0+4a_{P_0},
\qquad
P_{21}=\bar P_0+2a_{P_0},
\qquad
P_{22}=\bar P_0-4a_{P_0}.
\]
So the worst surviving lane is always the one with the larger sign of `a_{P_0}`. Therefore the exact robust all-lane ceiling test collapses to
\[
\boxed{
\bar P_0+4|a_{P_0}|
\le
P_{\rm crit}.
}
\]
Equivalently, in the one-scalar `\Xi_1` language,
\[
\boxed{
\bar P_0\bigl(1+|\epsilon\Xi_1|\bigr)
\le
P_{\rm crit}.
}
\]
Using the exact normalization compiler,
\[
\boxed{
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}igl(1+|\epsilon\Xi_1|\bigr)
\le
P_{\rm crit}.
}
\]

This is the cleanest same-charge continuation reached so far.

The actual moving-throat branch now lives or dies, at this transported level, by one explicit inequality in

- the normalization defect `\Delta_{\rm norm}`,
- the source-map factor `\hat m_0`,
- and the weak-axisymmetric outgoing-prefactor slope `\Xi_1=P_1/P_0`.

### 5.1 Exact calibrated-branch bound with weak anisotropy included

If the actual branch is exactly calibrated,
\[
\Delta_{\rm norm}=0,
\]
then the robust transported ceiling becomes
\[
\boxed{
\hat m_0^{\,2}
\ge
\frac{T_{\rm quad}(1+|\epsilon\Xi_1|)}{P_{\rm crit}}.
}
\]
So weak-axisymmetric prefactor loading raises the lower bound on the source-map factor linearly in the absolute outgoing slope.

That is a very concrete actual-branch falsification test.

---

## 6. Explicit headroom at the Stage-006 compatibility point

Stage 006 found the concrete compatibility point
\[
\bar P_0 = P_{0,\rm target,compat}
\approx 0.002069792318062885.
\]
Substituting this into the robust weak-axisymmetric ceiling law
\[
|\epsilon\Xi_1|\le \frac{P_{\rm crit}}{\bar P_0}-1,
\qquad
|a_{P_0}|\le \frac{P_{\rm crit}-\bar P_0}{4},
\]
gives the following explicit budgets.

### 6.1 Stricter `10%`-loss benchmark

For **both wall-like poles** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 0.367930328492646,
\qquad
|a_{P_0}|
\lesssim 1.90384841874108\times 10^{-4}.
\]
For a **nonempty wall-like corridor** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 0.737619063660757,
\qquad
|a_{P_0}|
\lesssim 3.81679567905443\times 10^{-4}.
\]

### 6.2 Looser `30%`-loss benchmark

For **both wall-like poles** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 2.94889585703134,
\qquad
|a_{P_0}|
\lesssim 1.52590049791274\times 10^{-3}.
\]
For a **nonempty wall-like corridor** to remain alive,
\[
|\epsilon\Xi_1|
\lesssim 4.63505472371892,
\qquad
|a_{P_0}|
\lesssim 2.39840016523863\times 10^{-3}.
\]

So the Stage-006 compatibility point still has finite weak-axisymmetric headroom.
But the stricter `10%` robust budget is not huge. That is exactly the sort of narrow corridor we would expect if the idea is real but difficult.

---

## 7. What Stage 007 changes

Stage 007 does not yet prove success or failure.
But it does convert the remaining ambiguity into an actual-branch compiler.

Before this stage, the next step was still phrased loosely as

> “extract the actual branch-compatible normalization target from the PDE.”

After this stage, we know the exact thing to compute:

1. the real branch packet gives `\Delta_{\rm norm}`, `a_{P_0}`, and `b_{P_0}`,
2. those compile to the actual lane prefactors `P_{20},P_{21},P_{22}`,
3. on the weak-axisymmetric branch they collapse to one scalar `\Xi_1=P_1/P_0`,
4. and the same-charge corridor survives only if those data satisfy one explicit finite inequality.

So the best current summary is:

> the idea is still alive, but the actual moving-throat branch now has to land inside a sharply delimited corridor in `(\Delta_{\rm norm},\Xi_1)` or equivalently `(\Delta_{\rm norm},a_{P_0},b_{P_0})`.

That is the first genuine PDE-selected same-charge kill test.

same_charge_barrier_audit_stage008_microscopic_xi1_primitive_compiler.md

# Same-Charge Barrier Audit — Stage 008: Microscopic `\Xi_1` Compiler, First-Order Conservative Compensation Surface, and the Mixed-Sector Survival Sieve

## 0. Purpose

Stage 007 converted the actual-branch same-charge test into one explicit transported inequality in
\[
\Delta_{\rm norm},\qquad \Xi_1=\frac{P_1}{P_0}.
\]

So the next honest step is no longer another abstract branch-packet manipulation.
It is to compute `\Xi_1` **microscopically** on the explicit finite-throat one-port branch, while keeping track of the conservative grouped-`P2` conditions that the `5`PN branch still has to respect.

That is what this stage does.

The main outputs are:

1. the exact arbitrary-base first-order formulas for
   \[
   u_2^{(1)},\qquad u_4^{(1)},\qquad \Xi_1,
   \]
   together with the exact compensation surface that preserves the conservative grouped response on a one-pole base branch;
2. the specialization of those formulas to the explicit Stage-006 compatibility point;
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
6. the direct translation of that family into the transported Stage-007 same-charge ceiling budgets.

So after Stage 008, the question is no longer

> can some microscopic anisotropy produce a useful `\Xi_1`?

It is now

> which microscopic families survive the conservative first-order compensation surface, and what same-charge headroom do they actually leave?

---

## 1. Frozen input carried forward

### 1.1 Explicit isotropic one-port base branch from Stage 006

Keep the same explicit finite-throat branch used in Stages 005–007:

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

For the concrete Stage-006 compatibility sample
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
\left(\frac12,\frac3{10},\frac25,\frac14,1,\frac75,2,1\right),
\]
the explicit base values are
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
and the exact one-pole condition holds:
\[
u_4-4u_2^2=0.
\]
The same branch also reproduces the carried Stage-006 compatibility-point prefactor
\[
P_{0,\rm target,compat}
=
\frac{N_0}{D_0}
\approx 0.002069792318062885.
\]

### 1.2 Stage-007 transported same-charge ceilings at that point

At the same compatibility point, the robust weak-axisymmetric budget from Stage 007 is
\[
|\epsilon\Xi_1| \lesssim 0.367930328492646
\]
for the stricter `10%`-loss “both wall-like poles survive” criterion, and
\[
|\epsilon\Xi_1| \lesssim 0.737619063660757
\]
for the stricter `10%`-loss “nonempty wall-like corridor” criterion.
The looser `30%` budgets are much larger and are carried forward unchanged.

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
\Xi_1=rac{P_1}{P_0}=rac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
}
\]

So on any base branch the same-charge scalar is already a **static loading mismatch** between

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

This is the exact arbitrary-base continuation of the canonical `5`PN even-preserving surface. Once it is imposed, the only remaining first-order outlet is `\Xi_1`.

For the concrete Stage-006 compatibility point, the compensation surface becomes
\[
D_{21}\approx -0.0489171640391802\,D_{01},
\qquad
D_{41}\approx -0.00717866681290820\,D_{01}.
\]

---

## 4. Primitive logarithmic slope compiler

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

On the concrete Stage-006 compatibility point, the microscopic same-charge scalar compiles numerically to
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
So on this branch the strongest positive same-charge leverage sits in the mixed-channel loading `x_{\lambda_W}`, while the strongest negative leverage sits in the mixed frequency drift `x_{\Omega_W}`.

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
Therefore:
\[
\boxed{\text{wall-only compensated deformations are generically trivial.}}
\]

This is already enough to kill the naive “pure wall anisotropy” route on the actual Stage-006 base point, since there
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
On the Stage-006 compatibility point,
\[
\Delta_{\rm BdG}\approx -5.11886996120011\times 10^{-5}\neq 0.
\]
So:
\[
\boxed{\text{the pure BdG compensated family is also trivial on the concrete branch.}}
\]

This is the same structural conclusion the later `5`PN continuation notes reached in a different language: neither pure wall nor pure BdG anisotropy carries the live weak-axisymmetric corridor by itself.

### 5.3 Mixed-sector-only family — explicit surviving corridor

Now activate only the mixed/U family
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}).
\]
On the Stage-006 compatibility point, the compensation matrix is
\[
\begin{pmatrix}
-0.241952861865934 & -0.122133861432532 & -0.0656784156312263 & 0.553209522700447 & 0.288144673113677 \\
-0.250543086743604 & -0.0937748521387244 & -0.0899548469020231 & 0.881694465041011 & 0.325834311088034
\end{pmatrix}.
\]
It has rank `2`, hence nullity `3`.
A convenient null basis is:
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

Therefore:
\[
\boxed{\text{the idea survives this stage only in a constrained mixed-sector corridor.}}
\]

---

## 6. Direct same-charge headroom on the first surviving mixed family

Choose the first surviving mixed basis vector `v_1` and write its microscopic amplitude as `t`. Then
\[
\Xi_1 = \sigma_1 t,
\qquad
\sigma_1\approx 1.36026097049402.
\]
The transported Stage-007 ceiling law is
\[
|\epsilon\Xi_1| \le \text{budget},
\]
so on this family it becomes
\[
|\epsilon t| \le \frac{\text{budget}}{\sigma_1}.
\]
The explicit budgets are:

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

So the strict transported window is not huge, but it is definitely nonzero. That is exactly the sort of result we would expect if the corridor is real but narrow.

---

## 7. What Stage 008 changes

Before this stage, the same-charge continuation was still phrased as

> compute `\Xi_1` somehow from the branch and compare it against the Stage-007 ceiling.

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

> the idea survives this stage, but only as a constrained mixed-sector corridor. Pure wall-only or pure support/BdG anisotropy does not survive the first conservative same-charge compensation test.

That is a real narrowing, and it is exactly the kind of narrowing we wanted.

same_charge_barrier_audit_stage009_strict_5pn_even_gate_mixed_corridor.md

# Same-Charge Barrier Audit — Stage 009: Strict `5`PN Even-Gate Package, the Surviving Mixed Corridor, and the Pure-Transfer Subcorridor

## 0. Purpose

Stage 008 showed that the explicit finite-throat one-port branch still supports a nontrivial mixed-sector corridor after imposing the **first-order conservative compensation surface** that keeps the actual grouped response fixed on that branch.

But the later `5`PN notes impose a stricter first-order package than Stage 008 used. In that later language, the surviving weak-axisymmetric grouped problem is organized by three scalars:

a) the load defect
\[
\Xi_{\rm load}=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0},
\]

b) the first even gate
\[
K_1=D_{21}+\frac{D_{01}}{9},
\]

c) the hidden-even gate
\[
H_{\rm even}=D_{41}-\frac{2}{3}D_{21}-\frac{D_{01}}{27}.
\]

So the next honest step is not another loose mechanism scan. It is:

> carry the Stage-008 mixed-sector survivor into the stricter imported `5`PN even-gate package and see whether the corridor collapses or survives.

That is what this stage does.

The main outputs are:

1. the exact bridge
   \[
   \Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0},
   \]
   so the same-charge scalar from Stage 008 is already the imported `5`PN load defect;
2. the exact comparison between the Stage-008 compensation surface and the stricter `5`PN even-gate package;
3. the explicit mixed-sector-only strict-even-gate solve on the concrete Stage-006 compatibility point;
4. the sharper **pure-transfer** subcorridor that survives once one also enforces the Stage-008 conservative-shape preservation on this noncanonical sample branch;
5. the transported same-charge ceiling budgets on both the strict mixed corridor and the pure-transfer subcorridor.

So after Stage 009, the question is no longer

> does some mixed anisotropy survive?

It is now

> does the actual moving-throat branch realize the surviving same-charge effect primarily as mixed-sector outgoing-transfer enhancement, with the conservative one-pole bundle frozen at first order?

---

## 1. Frozen input carried forward

### 1.1 Explicit one-port compatibility branch from Stages 006–008

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
Stage 008 fixed the isotropic one-pole base data
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

### 1.2 Stage-008 first-order primitive compiler

Stage 008 already compiled the primitive weak-axisymmetric slopes
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W})
\]
into
\[
D_{01},\qquad D_{21},\qquad D_{41},\qquad N_{01},\qquad \Xi_1.
\]
On the mixed-only sector, the Stage-008 conservative compensation surface was
\[
D_{21}=-u_2 D_{01},
\qquad
D_{41}=\frac{D_4}{D_0}D_{01},
\]
and the transported same-charge scalar was already
\[
\Xi_1=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]

### 1.3 Stage-007 transported same-charge windows

At this same compatibility point the robust carried budgets were
\[
|\epsilon\Xi_1|\lesssim 0.367930328492646
\]
for the stricter `10%`-loss “both wall-like poles survive” criterion, and
\[
|\epsilon\Xi_1|\lesssim 0.737619063660757
\]
for the stricter `10%`-loss “nonempty wall-like corridor” criterion. The looser `30%` windows are kept as secondary reference budgets.

---

## 2. Exact bridge to the imported `5`PN load defect

The imported weak-axisymmetric `5`PN load defect is
\[
\Xi_{\rm load}:=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]
But Stage 008 already gave
\[
\Xi_1=\frac{P_1}{P_0}=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]
So coefficientwise,
\[
\boxed{\Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0}.}
\]

That matters because it means we are not introducing a new scalar here. We are simply re-reading the Stage-008 same-charge scalar in the stricter imported `5`PN language.

So the same-charge scalar already is the `5`PN loading defect.

---

## 3. Exact comparison: Stage-008 compensation is weaker than the strict `5`PN package

The Stage-008 first-order conservative compensation surface on an arbitrary one-pole base branch is
\[
D_{21}=-u_2 D_{01},
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

So the Stage-008 conservative-shape preservation and the stricter imported `5`PN even gates are **not** the same condition. They agree only on the canonical branch for which the coefficients above vanish.

On the explicit Stage-006 compatibility branch the coefficients are numerically
\[
K_1\approx 0.0621939470719309\,D_{01},
\qquad
H_{\rm even}\approx -0.0116042611571584\,D_{01}.
\]
Both are nonzero. Therefore, on this noncanonical sample branch,
\[
\boxed{\text{Stage-008 compensation + strict `5`PN even gates} \iff D_{01}=0.}
\]

This is the first real sharpening of the corridor.

It means that once the stricter imported `5`PN gates are imposed on top of the Stage-008 compensation surface, the surviving same-charge corridor can no longer use first-order conservative static loading. It must pass through a branch with
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
So the strict even-gate package does **not** kill the mixed-sector corridor. It only deforms it relative to the weaker Stage-008 compensation surface.

The ambient-Euclidean operator norm of `\Xi_1` on this strict three-dimensional corridor is
\[
\boxed{\sigma_{\rm even}\approx 2.67386816837173.}
\]
That gives a canonical same-charge gain scale for unit microscopic mixed-sector drift amplitude.

---

## 5. The pure-transfer subcorridor

Now impose the full intersection of

1. Stage-008 conservative-shape preservation, and
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

So a **two-dimensional** mixed-sector corridor still survives even after imposing both the Stage-008 compensation surface and the imported strict `5`PN even gates.

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

The ambient-Euclidean operator norm of `\Xi_1` on this pure-transfer two-dimensional corridor is
\[
\boxed{\sigma_{\rm transfer}\approx 2.31561904386057.}
\]

---

## 6. Transported same-charge ceiling budgets on the strict corridors

Interpret the ambient microscopic mixed-sector drift amplitude as
\[
\|x_{\rm mixed}\|_2=t.
\]
If the operator norm of `\Xi_1` on a corridor is `\sigma`, then the transported Stage-007 ceiling law becomes
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

## 7. What Stage 009 changes

Before this stage, the strongest positive statement was only

> a mixed-sector corridor survives the first-order conservative compensation surface from Stage 008.

After this stage, the picture is much sharper.

1. The imported strict `5`PN loading defect is **exactly** the same scalar already isolated in Stage 008:
   \[
   \Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0}.
   \]
2. The Stage-008 compensation surface is weaker than the stricter imported `5`PN even-gate package.
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

> the idea survives Stage 009, but no longer as generic mixed anisotropy. The sharpest surviving mechanism is mixed-sector outgoing-transfer enhancement with the conservative one-pole bundle frozen at first order.

That is a real narrowing, and it is the kind of narrowing we want.

same_charge_barrier_audit_stage010_pure_transfer_load_factor_rigidity.md

# Same-Charge Barrier Audit — Stage 010: Pure-Transfer Load Factor, Outgoing-Rigidity Sieve, and the First Co-Loading No-Go

## 0. Purpose

Stage 009 narrowed the surviving same-charge mechanism to a very specific subcorridor:
\[
D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\Xi_{\rm load}=\frac{N_{01}}{N_0}.
\]

So the next honest step is no longer another generic mixed-anisotropy scan. The surviving effect is already telling us what it is:

> it is a **pure outgoing-load** problem.

That means the right next question is:

> once the conservative one-pole bundle is frozen at first order, what microscopic one-port load factor is actually carrying \(\Xi_1\), and does that mechanism survive the first natural rigidity filters?

This stage answers that.

The main outputs are:

1. the exact pure-transfer theorem
   \[
   \Xi_1=\frac{N_{01}}{N_0}=2\frac{P_{01}}{P}-2\frac{\Delta_{01}}{\Delta},
   \]
   on the Stage-009 pure-transfer corridor;
2. the exact one-port load-factor decomposition
   \[
   \Lambda:=\frac{P}{\Delta}
   =
   \frac{G_W}{\Omega_W^2}\,\frac{1+I}{1-H},
   \]
   with
   \[
   I=\frac{RG_U}{\Omega_U^2G_W},
   \qquad
   H=\frac{R^2}{\Omega_U^2\Omega_W^2};
   \]
3. the exact pure-transfer identity
   \[
   \Xi_1=2\,\delta\ln\Lambda;
   \]
4. the outgoing-rigidity sieve:
   - \(i=0\),
   - \(h=0\),
   - \(m=0\),
   - and the combined pair \(i=h=0\);
5. the first exact same-charge **co-loading no-go**:
   \[
   i=h=0
   \quad\Longrightarrow\quad
   \text{only the trivial pure-transfer drift survives;}
   \]
6. and the transported dynamic ceilings on the remaining one-dimensional rigid survivors.

So after Stage 010, the question is no longer

> does some pure-transfer effect survive?

It is now

> which outgoing-load channels are allowed to move together, and which rigidity assumptions kill the same-charge corridor outright?

---

## 1. Frozen input carried forward

### 1.1 Concrete finite-throat one-port sample branch

Keep the same explicit branch from Stages 006–009:

- wall and brane-like internal coordinate on the lowest N/N zero mode,
- trapped support and mixed coordinate on the lowest D/N half-wave,
- exact overlap constant
  \[
  \kappa=\frac{2\sqrt2}{\pi}.
  \]

The primitive parameters remain
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
\left(\frac12,\frac3{10},\frac25,\frac14,1,\frac75,2,1\right).
\]

The one-port mixed primitives are
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

### 1.2 Stage-009 pure-transfer subcorridor

Stage 009 showed that on this noncanonical sample branch the full intersection of

1. the Stage-008 conservative-shape preservation, and
2. the stricter imported `5`PN even-gate package

is equivalent to
\[
D_{01}=D_{21}=D_{41}=0.
\]

So on that two-dimensional mixed-sector subcorridor the same-charge scalar becomes
\[
\boxed{
\Xi_1=\frac{N_{01}}{N_0}.
}
\]

That is the starting point of this stage.

---

## 2. Exact pure-transfer load theorem

Differentiate
\[
N_0=\frac{P^2}{\Delta^2}.
\]
Then on any weak first-order deformation,
\[
\frac{N_{01}}{N_0}=2\frac{P_{01}}{P}-2\frac{\Delta_{01}}{\Delta}.
\]
Because the pure-transfer subcorridor has
\[
D_{01}=0,
\]
the surviving same-charge scalar is exactly
\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
=
2\frac{P_{01}}{P}
-
2\frac{\Delta_{01}}{\Delta}.
}
\]

So the surviving same-charge mechanism is literally the logarithmic slope of the one-port load factor
\[
\Lambda:=\frac{P}{\Delta}.
\]

Equivalently,
\[
\boxed{
\Xi_1 = 2\,\delta\ln\Lambda.
}
\]

This is already a major sharpening of the Stage-009 verdict.

The mechanism is no longer “mixed-sector enhancement” in general. It is one exact load-factor slope.

---

## 3. Exact one-port factorization

Write
\[
P=\Omega_U^2G_W+RG_U,
\qquad
\Delta=\Omega_U^2\Omega_W^2-R^2.
\]
Then
\[
P=\Omega_U^2G_W\left(1+\frac{RG_U}{\Omega_U^2G_W}\right),
\qquad
\Delta=\Omega_U^2\Omega_W^2\left(1-\frac{R^2}{\Omega_U^2\Omega_W^2}\right).
\]
So the load factor becomes
\[
\boxed{
\Lambda
=
\frac{P}{\Delta}
=
\frac{G_W}{\Omega_W^2}\,
\frac{1+I}{1-H},
}
\]
with the exact one-port invariants
\[
\boxed{
I=\frac{RG_U}{\Omega_U^2G_W},
\qquad
H=\frac{R^2}{\Omega_U^2\Omega_W^2}.
}
\]

Define the microscopic logarithmic drifts
\[
m:=\delta\ln\!\left(\frac{G_W}{\Omega_W^2}\right),
\qquad
i:=\delta\ln I,
\qquad
h:=\delta\ln H.
\]
Then on the Stage-009 pure-transfer subcorridor,
\[
\boxed{
\Xi_1
=
2\left[
m+\frac{I}{1+I}\,i+\frac{H}{1-H}\,h
\right].
}
\]

So the outgoing-load scalar has three transparent pieces:

1. a direct mixed-leg / port-weight piece \(m\),
2. an interference-ratio piece \(i\),
3. a hybridization-ratio piece \(h\).

### 3.1 Exact sample-branch coefficients

On the concrete sample branch,
\[
I=\frac{3}{16},
\qquad
H=\frac{25}{98\pi^2}.
\]
So the exact Stage-010 load law is
\[
\boxed{
\Xi_1
=
2m+\frac{6}{19}i+\frac{50}{98\pi^2-25}\,h
}
\]
on the pure-transfer subcorridor.

This is the first exact microscopic decomposition of the surviving same-charge scalar.

---

## 4. The outgoing-rigidity sieve

Now impose the first natural rigidity filters.

### 4.1 Combined interference and hybridization rigidity

If both
\[
i=0,
\qquad
h=0,
\]
are imposed on the pure-transfer corridor, the exact reduced `2 x 2` rigidity matrix on the Stage-009 basis has nonzero determinant. The audit factors it exactly and finds
\[
\det(i,h)|_{\rm pure\ transfer}\neq 0.
\]

So the combined rigidity system has only the trivial solution:
\[
\boxed{
D_{01}=D_{21}=D_{41}=0,\quad i=0,\quad h=0
\ \Longrightarrow\
x_{\rm mixed}=0.
}
\]

This is the first exact same-charge co-loading no-go of the audit.

It says:

> if both the interference ratio and the hybridization ratio are frozen, the pure-transfer same-charge mechanism dies on this concrete branch.

### 4.2 Single-rigidity survivors

The situation is different if only one rigidity is imposed.

The audit finds:

- \(i=0\) leaves a one-dimensional survivor,
- \(h=0\) leaves a one-dimensional survivor,
- \(m=0\) also leaves a one-dimensional survivor.

So the corridor is narrower, but it does not die under a **single** outgoing-rigidity filter.

This is useful physically. It means the same-charge branch still has room to live if one of the outgoing subchannels is rigid, but not if both \(i\) and \(h\) are rigid simultaneously.

---

## 5. Concrete unit directions and same-charge gain scales

Using ambient Euclidean normalization in the mixed primitive space
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}),
\]
the audit finds the following unit directions.

### 5.1 `i = 0` survivor
\[
v_i\approx
(0.45280825,\,-0.29424612,\,-0.82815170,\,-0.04054866,\,0.14458380),
\]
with
\[
|\Xi_1(v_i)|\approx 1.26576248.
\]

### 5.2 `h = 0` survivor
\[
v_h\approx
(0.66561963,\,-0.38941932,\,0.46712837,\,0.03609301,\,0.43103536),
\]
with
\[
|\Xi_1(v_h)|\approx 2.04509123.
\]

### 5.3 `m = 0` survivor
\[
v_m\approx
(0.13386239,\,-0.10586713,\,-0.98242900,\,-0.05389175,\,-0.05293356),
\]
with
\[
|\Xi_1(v_m)|\approx 0.29342952.
\]

The `m = 0` result is especially interesting. It means the direct mixed-leg factor can be frozen while the same-charge signal is carried entirely by a correlated interference/hybridization deformation.

So the pure-transfer corridor is not synonymous with a raw mixed-leg effect.

---

## 6. Transported dynamic ceilings on the rigid one-dimensional survivors

Carry forward the same Stage-007 transported budgets used in Stage 009. Interpreting the ambient microscopic drift amplitude as
\[
\|x_{\rm mixed}\|_2=t,
\]
the dynamic ceiling is still
\[
|\epsilon|t \le \frac{\text{budget}}{\sigma},
\]
where \(\sigma\) is the \(\Xi_1\) operator norm on the relevant corridor.

### 6.1 Reference pure-transfer corridor
\[
\sigma_{\rm transfer}\approx 2.31561904.
\]
So the stricter `10%`-loss ceilings are
\[
|\epsilon|t \lesssim 0.15889070
\qquad
(\text{both wall-like poles survive}),
\]
\[
|\epsilon|t \lesssim 0.31854077
\qquad
(\text{nonempty wall-like corridor}).
\]

### 6.2 `i = 0` survivor
\[
\sigma_i\approx 1.26576248.
\]
So the stricter `10%`-loss ceilings become
\[
|\epsilon|t \lesssim 0.29067881,
\qquad
|\epsilon|t \lesssim 0.58274682.
\]

### 6.3 `h = 0` survivor
\[
\sigma_h\approx 2.04509123.
\]
So the stricter `10%`-loss ceilings become
\[
|\epsilon|t \lesssim 0.17990900,
\qquad
|\epsilon|t \lesssim 0.36067783.
\]

### 6.4 `m = 0` survivor
\[
\sigma_m\approx 0.29342952.
\]
So the stricter `10%`-loss ceilings become
\[
|\epsilon|t \lesssim 1.25389678,
\qquad
|\epsilon|t \lesssim 2.51378617.
\]

The `m = 0` branch leaves the most headroom simply because the same-charge scalar is much smaller per unit ambient microscopic drift there.

That does **not** make it automatically the physically best mechanism, but it does make it the least constrained by the transported dynamic ceiling.

---

## 7. What Stage 010 changes

Before this stage, the strongest statement was only

> the strict even-gate survivor is a pure-transfer mixed corridor.

After this stage, the picture is much sharper.

1. The surviving same-charge scalar is exactly one **outgoing load-factor slope**:
   \[
   \Xi_1=2\,\delta\ln\Lambda.
   \]
2. On the concrete branch, that load factor splits into three pieces:
   \[
   m,\qquad i,\qquad h.
   \]
3. Freezing both \(i\) and \(h\) kills the corridor outright.
4. Freezing only one of \(i,h,m\) still leaves a one-dimensional same-charge survivor.
5. Even with \(m=0\), the same-charge effect can survive through correlated interference/hybridization motion.

So the best current summary is:

> the idea survives Stage 010, but only as a very structured outgoing co-loading effect. Pure transfer is real, but it is not generically “just the mixed leg.” And the first exact no-go now says that simultaneous interference and hybridization rigidity kills the mechanism on this concrete branch.

That is a real narrowing, and it is the right kind of narrowing.

same_charge_barrier_audit_stage011_numerator_denominator_split.md

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

same_charge_barrier_audit_stage012_selected_branch_signature.md

# Same-Charge Barrier Audit — Stage 012: Selected-Branch Numerator/Denominator Signature and the Softening-Depth Crossover Theorem

## 0. Purpose

Stage 011 split the surviving same-charge pure-transfer corridor into two exact 1D static subcorridors:

- **numerator-rigid** (\(\pi_1=0\)),
- **denominator-rigid** (\(\delta_1=0\)),

and showed that both survive the first actual wall-like dynamic window on the concrete compatibility point.

So the next honest question is sharper:

> which of those two rigid static signatures does the **actual PDE-selected mixed branch** most closely resemble?

The moving-throat selected-branch notes already reduced the real quadrupole normalization branch to the exact scalar softening-depth problem
\[
N_-(x)=\frac{\beta_0\,s_-(x)^2}{\kappa_0^2\,(A-x)},
\qquad
0\le x < A,
\]
with the exact selected overlap
\[
s_-(x)=\frac{\bigl[\kappa_0^2(x+\Delta K_{\rm ax})+\kappa_1^2x\bigr]^2}
{\kappa_0^2(x+\Delta K_{\rm ax})^2+\kappa_1^2x^2}.
\]

This stage carries the Stage-011 split into that actual selected-branch language.

The main outputs are:

1. the exact factorization of the selected-branch normalization product into a **numerator-like** and a **denominator-like** piece;
2. the exact log-slope classifier that decides which split the selected branch resembles at any point on the stable branch;
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
   separating always-denominator-dominant branches from branches that begin numerator-dominant;
5. and the conclusion that the real selected branch is **never** literally one of the rigid Stage-011 subcorridors — it is an exact co-loading product — but it becomes unambiguously **denominator-like** near softening, and for all \(\delta\ge 8/9\) it is denominator-like on the entire stable branch.

So after Stage 012, the next question is no longer “numerator-rigid or denominator-rigid?” in the abstract.
It is:

> what are the actual selected-branch ratios \((\xi,\delta)\) on the physical moving-throat branch, and where do they land on this universal classifier map?

---

## 1. Frozen input carried forward

### 1.1 Stage-011 static split

Stage 011 isolated the exact static pure-transfer identity
\[
\Xi_1 = 2(\pi_1-\delta_1),
\qquad
\pi_1:=\frac{P_{01}}{P},
\qquad
\delta_1:=\frac{\Delta_{01}}{\Delta},
\]
and then split the two-dimensional pure-transfer corridor into:

- the 1D numerator-rigid branch \(\pi_1=0\),
- the 1D denominator-rigid branch \(\delta_1=0\).

That gave the right static classifiers, but it was still a classifier internal to the primitive mixed-slope space.

### 1.2 Actual selected-branch product

The later moving-throat selected-branch chain already replaced the free static load-factor language by the exact selected-mode normalization product
\[
N_-(x)=\frac{\beta_0\,s_-(x)^2}{\kappa_0^2\,(A-x)},
\]
with the exact D/N constants
\[
\kappa_0^2=\frac{8}{\pi^2},
\qquad
\kappa_1^2=\frac{16}{9\pi^2}.
\]

The point of Stage 012 is to compare the Stage-011 numerator/denominator split against **this** actual selected-branch object rather than against a free one-port static factor.

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

So the same universal Stage-25 branch function \(F\) still controls the selected normalization product; the extra factor \(8/\pi^2=\kappa_0^2\) is a fixed D/N overlap constant and does not affect the numerator/denominator classifier below.

Now split this into

\[
\boxed{
F(\xi,\delta)
=
F_{\rm num}(\xi,\delta)\,F_{\rm den}(\xi),
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

But unlike Stage 011, the split is not “either/or.” It is an exact product.

---

## 3. Exact log-slope classifier

To decide which Stage-011 rigid branch the selected branch most closely resembles, the right invariant is the log-slope split of \(F\) along the physical softening coordinate \(\xi\).

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
\frac{72\delta^2(1-\xi)}
{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}.
}
\]

Interpretation:

- \(\mathcal R_{ND}>1\): the selected branch is **numerator-like** at that point;
- \(\mathcal R_{ND}<1\): it is **denominator-like**;
- \(\mathcal R_{ND}=1\): exact crossover.

This is the Stage-012 replacement for the Stage-011 rigid subcorridors.

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

As \(\xi\to 1^-\),
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

This is already a strong answer to the Stage-011 continuation question:

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

## 7. What this says about the Stage-011 rigid subcorridors

Stage 011 asked which rigid static signature the real PDE-selected mixed branch most closely resembles.

Stage 012 gives the precise answer.

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

So if the physical branch hits the universal target at appreciable softening depth, the right Stage-011 proxy is the denominator-rigid one, not the numerator-rigid one.

### 7.3 Only very early on can it look numerator-like

Numerator-like behavior is confined to the small-softening regime
\[
0\le \xi < \xi_*(\delta),
\qquad
0<\delta<\frac89.
\]

So the numerator-rigid Stage-011 branch is best read as an **onset-side local proxy**, not as the global selected-branch signature.

---

## 8. Best current verdict after Stage 012

The continuation question from Stage 011 now has a clean answer.

1. The real selected PDE branch is not one of the rigid Stage-011 subcorridors.
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

> feed the actual moving-throat selected-branch data into \((\xi,\delta)\), place the physical branch on this universal classifier map, and then compare the resulting denominator-vs-numerator signature against the concrete Stage-011 dynamic ceilings.

That is the clean continuation point.

same_charge_barrier_audit_stage013_selected_branch_dynamic_window_compiler.md

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

same_charge_barrier_audit_stage014_continuum_placement_dynamic_class_map.md


# Same-Charge Barrier Audit — Stage 014: Continuum Placement Pullback of the Selected-Branch Dynamic-Class Map

## 0. Purpose

Stage 013 solved the exact dynamic-sign problem on the full selected-branch classifier half-line:
the wall-like dynamic window never became the first kill condition, and the only live global kill still came from the transported static `\Xi_1` budget.

But Stage 013 still spoke in the abstract classifier variable
\[
\mathcal R_{ND}.
\]
The next honest step is therefore very specific:

> pull that classifier and its dynamic sign thresholds back through the **actual continuum selected-branch placement map**, so the same-charge verdict is expressed directly in the physical kernel ratios of the moving-throat branch.

This stage does exactly that.

The main outputs are:

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

So after Stage 014, the same-charge corridor is still alive.
What changes is that the classifier map is no longer a sample-branch statement.
It is now a physical continuum-kernel statement.

---

## 1. Frozen input carried forward

### 1.1 Universal selected-branch geometry

From the carried Stage-025 universal D/N branch geometry, the stable selected branch is controlled by the exact functions
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
\partial_\xi F>0,
\qquad
\partial_\xi G>0,
\]
and the endpoint data are
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

### 1.2 Exact selected-branch classifier

From Stages 012–013, the exact selected-branch numerator/denominator classifier is
\[
\mathcal R_{ND}(\xi,\delta)
=
\frac{72\delta^2(1-\xi)}
{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)},
\]
with strict monotonicity
\[
\partial_\xi \mathcal R_{ND}<0.
\]

The carried Stage-013 sign threshold is
\[
\mathcal R_*\approx 1.229255438463336,
\]
and the denominator-like threshold is
\[
\mathcal R_{ND}=1.
\]

The associated onset thresholds in `\delta` are

\[
\delta_*^{(\rm dyn)}=\frac{8}{9\mathcal R_*}\approx 0.723111617875019,
\]
\[
\delta_{\rm den}=\frac89.
\]

So Stage 013 already gave two exact global statements:

- if \(\delta\ge \delta_*^{(\rm dyn)}\), the nonempty dynamic ceiling is infinite on the whole selected branch;
- if \(\delta\ge 8/9\), the whole selected branch is denominator-like.

### 1.3 Continuum placement map

From the carried Stage-026 continuum kernel extraction, the actual moving-throat branch is placed by the exact dimensionless ratios
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

This is the exact continuum pullback of the Stage-013 classifier.

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

This is the first real physical reading of the Stage-013 classifier map.

---

## 3. Exact threshold compiler for any classifier cap

Fix any classifier cap \(c>0\).
Define the threshold polynomial
\[
P_c(\xi,\delta)
:=
c(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)-72\delta^2(1-\xi).
\]

The exact derivative is
\[
\partial_\xi P_c
=
3\Bigl(87c\delta^2+198c\delta\xi+121c\xi^2+24\delta^2\Bigr)>0
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

## 5. Sample threshold values

The exact pulled-back thresholds on a few representative `\delta` slices are:

| \(\delta\) | \(\xi_{\rm flip}\) | \(R_{\rm flip}\) | \(\xi_{\rm den}\) | \(R_{\rm den}\) |
|---:|---:|---:|---:|---:|
| \(0.25\) | \(0.087442106\) | \(1.330868539\) | \(0.107223051\) | \(1.393832566\) |
| \(0.50\) | \(0.051428579\) | \(1.139956630\) | \(0.081847938\) | \(1.221087062\) |
| \(0.75\) | \(0\) | \(1\) | \(0.032505121\) | \(1.071471867\) |

So the physical selected branch becomes dynamically sign-safe quite early.
For \(\delta=0.75\), the nonempty dynamic ceiling is already infinite from onset, even though the branch is not yet denominator-like from onset because \(0.75<8/9\).

This is a strong continuation of the Stage-013 picture.

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
\frac{8\Lambda(1-\epsilon_W)}{\pi^2 R_{\rm flip}(\delta)}
\quad\Longleftrightarrow\quad
\text{nonempty dynamic ceiling infinite},
}
\]
\[
\boxed{
M_{\rm mix}
\le
\frac{8\Lambda(1-\epsilon_W)}{\pi^2 R_{\rm den}(\delta)}
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

Stage 013 already proved the global inequalities
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

> even on the actual continuum-selected branch, the first kill condition remains the transported static `\Xi_1` budget, not the wall-like dynamic window.

This is the central Stage-014 theorem.

---

## 8. Best current verdict after Stage 014

Stage 014 does not kill the same-charge corridor.
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

So after Stage 014, the same-charge idea is still alive.
The remaining first kill condition is still the **static placement of `\Xi_1` on the actual moving-throat branch**, not the wall-like dynamic window.

same_charge_barrier_audit_stage015_known_5pn_data_injection.md

# Same-Charge Barrier Audit — Stage 015: Known 5PN Data Injection and Current Branch Verdict

## Purpose

The earlier same-charge barrier audit reduced the live corridor to a very narrow condition set:

1. the **dynamic** wall-like window is *not* the first kill condition;
2. the **support/source** side must supply enough coherent enhancement;
3. the first unresolved kill test is the **static** placement / orbit-lock side.

This stage plugs the **actually numerically located** 5PN support/source data into that audit chain.

The point is not to solve the full moving-throat branch. The current 5PN notes are explicit that the support/source side has already been numerically located, while the actual PDE-selected orbit-lock / coherent placement point is still not numerically present in the files. So the clean question is:

> once the numerically located 5PN support/source data are inserted, does the same-charge corridor die, or does the unresolved bottleneck stay where the earlier audit said it was?

---

## 1. What is already numerically known from the 5PN stack

The numerically located Family-1 support/source branch gives the following exact values on the refreshed geometry branch:

\[
\Lambda_\ell \approx 36.94973154240256,
\qquad
\kappa \approx 2457.5087899001137,
\qquad
\zeta_{\max} \approx 2.4675297457259358.
\]

On the two explicit extraction branches:

### \(\chi\)-weighted extraction
\[
Pe_*^{(\chi)} \approx 11155.7265863205869,
\]
\[
\zeta_{\rm phys}^{(\chi)} \approx 2.4675296478814376,
\]
\[
\rho_{\alpha,\max}^{(\chi)} \approx 3.4675296478814376.
\]

### \(J\)-weighted extraction
\[
Pe_*^{(J)} \approx 2504.9703142859238,
\]
\[
\zeta_{\rm phys}^{(J)} \approx 2.4675278051675084,
\]
\[
\rho_{\alpha,\max}^{(J)} \approx 3.4675278051675084.
\]

The natural isotropic passive/outgoing grouped-\(P_2\) branch still requires only

\[
\zeta_{\rm req}=\frac13,
\qquad
\rho_\alpha^{\rm req}=\frac43.
\]

So the support/source branch is no longer just “probably okay.” It is numerically located and strongly above the canonical isotropic demand.

---

## 2. Exact margins after plugging in the known 5PN data

The support/source safety margins are:

### \(\chi\)-weighted branch
\[
\zeta_{\rm phys}^{(\chi)}-\zeta_{\rm req}
\approx 2.1341963145481043,
\]
\[
\rho_{\alpha,\max}^{(\chi)}-\frac43
\approx 2.1341963145481043.
\]

### \(J\)-weighted branch
\[
\zeta_{\rm phys}^{(J)}-\zeta_{\rm req}
\approx 2.1341944718341751,
\]
\[
\rho_{\alpha,\max}^{(J)}-\frac43
\approx 2.1341944718341751.
\]

Useful ratio form:

### \(\chi\)-weighted branch
\[
\frac{\zeta_{\rm phys}^{(\chi)}}{\zeta_{\rm req}}
\approx 7.402588943644313,
\qquad
\frac{\rho_{\alpha,\max}^{(\chi)}}{4/3}
\approx 2.600647235911438.
\]

### \(J\)-weighted branch
\[
\frac{\zeta_{\rm phys}^{(J)}}{\zeta_{\rm req}}
\approx 7.402583415502525,
\qquad
\frac{\rho_{\alpha,\max}^{(J)}}{4/3}
\approx 2.600645853875631.
\]

So the explicit support/source branch overshoots the canonical isotropic requirement by a factor of about \(7.4\) in the \(\zeta\) variable.

At the same time, the numerically selected Family-1 points sit very close to the Family-1 ceiling:

\[
\zeta_{\max}-\zeta_{\rm phys}^{(\chi)}
\approx 9.784449820573488\times 10^{-8},
\]
\[
\zeta_{\max}-\zeta_{\rm phys}^{(J)}
\approx 1.9405584274474645\times 10^{-6}.
\]

That is not a contradiction. It simply means the numerically selected Family-1 branch nearly saturates the *Family-1 support* ceiling while still sitting far above the much smaller isotropic demand \(\zeta_{\rm req}=1/3\).

---

## 3. Transported consequence for the same-charge audit chain

The earlier same-charge audit already proved two things before any 5PN injection:

1. the **dynamic** selected-branch wall-window is not the first kill condition;
2. the first real bottleneck is the **static** transported placement / \(\Xi_1\) side.

Injecting the numerically located 5PN support/source data does not change that ordering.
Instead it sharpens it:

- the support/source side is numerically safe;
- the canonical passive/outgoing normalization side remains exact on the natural branch;
- so the first unresolved kill condition stays where the 5PN notes say it is:
  the actual PDE-selected orbit-lock / coherent placement point.

In the 5PN finish-line language, the missing numerical object is still the actual branch point satisfying

\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln \epsilon_\eta=0,
\]
together with the canonical outgoing normalization condition
\[
N_Q=1.
\]

So the support/source side is no longer the place where the same-charge idea should die inside the current reduced theorem stack.

---

## 4. Current best verdict after plugging in the known numbers

The same-charge corridor is still alive.

More precisely:

1. **Static same-sign Maxwell shaping** is still not the answer.
2. **Dynamic resonance** is still not the first bottleneck.
3. **Support/source enhancement** is now numerically located and strongly non-bottlenecked.
4. The first unresolved gate is still the **actual PDE-selected static orbit-lock / placement point**.

So after Stage 015 the question is no longer:

> can the reduced support/source side possibly be large enough?

It already is.

The question is now:

> when the completed moving-throat branch is actually selected, does its realized orbit packet / static placement verdict land inside the surviving same-charge window, or does the route finally die there?

That is the cleanest current stopping point.

---

## 5. What the next honest stage should do

The next stage should not invent more support algebra. It should compile the actual unresolved packet into the same-charge audit language.

Concretely, the next theorem gate is:

1. take the exact unresolved coherent placement packet
   \[
   (R_{\rm tr},R_{\rm target},\epsilon_\eta),
   \]
   or equivalently
   \[
   (d\ln R_{\rm tr},d\ln R_{\rm target},d\ln\epsilon_\eta),
   \]
2. express its weak-axisymmetric verdict as the actual static \(\Xi_1\) / transported placement scalar used in the same-charge chain,
3. and test whether the realized branch clears the already-carried static ceiling.

That is where the present stack says the real answer now lives.

same_charge_barrier_audit_stage016_rigid_mouth_orbit_lock_gate.md

# Same-Charge Barrier Audit — Stage 016: Rigid-Mouth Orbit-Lock Compiler and the Static Turbulence Gate

## Status

**Exact within the carried Stage-007 / Stage-015 same-charge reduction plus the Stage-171 coherent branch-observable compiler.**

This stage does **not** claim that hydrodynamic turbulence or cavitation has already been derived from the full PDE.
It shows something narrower and more useful:

> if one adopts the physical reading that the brane mouth is rigid while the internal branch can still repackage loading, then the unresolved same-charge kill test is exactly a bound on the internal orbit/transfer packet, and the Stage-007 static ceiling can be read as a first reduced “choke / turbulence” gate.

The theorem content is exact inside the reduced branch language.
The words “turbulence”, “choked flow”, and “collapse” are **interpretive labels** for leaving that surviving static branch window.

---

## 0. Why this stage is needed

Stage 015 already showed that the explicit `5PN` support/source branch is numerically safe and that the first unresolved same-charge bottleneck is still the **PDE-selected orbit-lock / coherent placement point**. The same notes are explicit that the numerically missing object is the actual point satisfying the coherent placement conditions
\[
d\ln R_{\rm tr}=0,\qquad
d\ln R_{\rm target}=0,\qquad
d\ln \epsilon_\eta=0.
\]

So the natural next question is no longer whether support/source is large enough.
It is:

> how should the unresolved orbit packet be interpreted if the defect mouth is geometrically rigid and the remaining load can only reorganize internally?

This stage compiles that idea into exact reduced formulas.

---

## 1. Interpretive firewall

Two distinctions matter.

### 1.1 Rigid mouth is **not** identical to `D_{01}=0`

The geometric statement “the brane entrance is rigid” is a statement about the mouth geometry.
The algebraic condition
\[
D_{01}=0
\]
is a stronger statement: it says the **effective static grouped lane operator** does not drift at first weak-axisymmetric order.

So throughout this stage:

- **mouth rigidity / tracking lock** is represented first by
  \[
  \delta\ln R_{\rm tr}=0,
  \]
- while
  \[
  D_{01}=0
  \]
  is introduced only later as an additional **operator-rigidity hypothesis**.

### 1.2 The carried `L/a` lock is still a branch freeze, not a full-PDE theorem

The preferred aspect ratio
\[
L/a \approx 1.85
\]
is a carried branch value / reference freeze in the present stack, not yet a theorem that every loaded moving-throat branch must obey exactly.

So the physical reading “the throat does not deepen with more loading” should currently be read as:

> on the branch family we are actually carrying, the large-scale geometry is treated as effectively locked, and the first remaining freedom is internal transfer / placement rather than wholesale mouth growth.

That is strong enough for the present reduced audit, but it is not yet a finished nonlinear PDE theorem.

---

## 2. Exact branch-observable compiler from Stage 171

Stage 171 already gives the first-order observable compiler
\[
\Theta_1=\delta\ln R_{\rm tr},
\]
\[
\Xi_1=\delta\ln \mathfrak N_* - B_*\,\delta\ln R_{\rm tr},
\]
\[
\mathcal R_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln \epsilon_\eta
-\Xi_1,
\]
where
\[
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*}
\]
is the corrected nontracking branch observable.

So the full coherent weak-axisymmetric defect packet is already an exact linear compiler image of the direct observable packet
\[
\bigl(\delta\ln R_{\rm tr},\ \delta\ln\mathfrak N_*,\ \delta\ln\epsilon_\eta\bigr).
\]

This gives the first sharp reduced translation of the user’s picture:

- `R_tr` measures tracking / mouth-side placement,
- `\mathfrak N_*` measures corrected internal transfer / nontracking load,
- `\epsilon_\eta` measures selected-branch dressing.

---

## 3. Rigid-mouth / track-locked specialization

If the mouth-side branch is rigid at the first observable level, impose
\[
\boxed{\delta\ln R_{\rm tr}=0.}
\]

Then the exact Stage-171 compiler collapses to
\[
\boxed{\Theta_1=0,}
\]
\[
\boxed{\Xi_1=\delta\ln\mathfrak N_*,}
\]
\[
\boxed{
\mathcal R_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta
-\delta\ln\mathfrak N_*.
}
\]

So on a track-locked branch the direct same-charge placement scalar is **literally**
the logarithmic drift of the corrected internal transfer observable `\mathfrak N_*`.

That is the precise reduced sense in which “the entrance is rigid but the inside can still load” is already visible in the math.

---

## 4. Stronger operator-rigidity specialization

Stage 007 and Stage 008 already showed that the weak-axisymmetric static load scalar is
\[
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}
=
\frac{P_1}{P_0}.
\]

If one imposes the stronger operator-rigidity hypothesis
\[
\boxed{D_{01}=0,}
\]
then
\[
\boxed{
\Xi_{\rm load}
=
\frac{N_{01}}{N_0}.
}
\]

So on a branch where

1. the mouth/track side is rigid,
2. the static entrance operator is also rigid,

the entire weak-axisymmetric load defect is reinterpreted as **pure internal outgoing-transfer strengthening**.

This is the closest exact reduced analogue of the “internal choke opens while the entrance does not” picture.

---

## 5. Transported static ceiling as an internal choke / turbulence gate

Stage 007 already transported the primitive-family same-charge window onto the actual branch packet:
\[
\bar P_0(1+|\epsilon\Xi_1|)\le P_{\rm crit},
\]
equivalently
\[
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}(1+|\epsilon\Xi_1|)\le P_{\rm crit},
\qquad
T_{\rm quad}:=\frac{54Gc_s^5}{5a^5c^5}.
\]

So for any actual branch packet, the first robust static same-charge survival gate is
\[
\boxed{
|\epsilon\Xi_1|
\le
\frac{P_{\rm crit}\hat m_0^{\,2}}{\Delta_{\rm norm}+T_{\rm quad}}-1.
}
\]

On a calibrated branch with
\[
\Delta_{\rm norm}=0,
\]
this becomes
\[
\boxed{
|\epsilon\Xi_1|
\le
\frac{P_{\rm crit}\hat m_0^{\,2}}{T_{\rm quad}}-1.
}
\]

Now add the rigid-mouth specialization from Section 3:
\[
\Xi_1=\delta\ln\mathfrak N_*.
\]
Then the same survival gate becomes
\[
\boxed{
|\epsilon\,\delta\ln\mathfrak N_*|
\le
\frac{P_{\rm crit}\hat m_0^{\,2}}{\Delta_{\rm norm}+T_{\rm quad}}-1.
}
\]

And with the stronger operator-rigidity closure `D_{01}=0`,
\[
\Xi_{\rm load}=\Xi_1=\frac{N_{01}}{N_0},
\]
so the same gate is
\[
\boxed{
\left|\epsilon\,\frac{N_{01}}{N_0}\right|
\le
\frac{P_{\rm crit}\hat m_0^{\,2}}{\Delta_{\rm norm}+T_{\rm quad}}-1.
}
\]

This is the exact reduced formula behind the proposed “choked-flow / turbulence threshold” reading.

---

## 6. The Stage-007 numerical budgets in this language

At the Stage-006 compatibility point
\[
\bar P_0 \approx 0.002069792318062885,
\]
Stage 007 gave the strict `10%` robust budget
\[
|\epsilon\Xi_1|\lesssim 0.367930328492646
\]
for **both wall-like poles** to remain alive, and the looser nonempty-corridor budget
\[
|\epsilon\Xi_1|\lesssim 0.737619063660757.
\]

So under the rigid-mouth specialization these are directly reinterpreted as

### Robust static gate
\[
\boxed{
|\epsilon\,\delta\ln\mathfrak N_*|
\lesssim 0.367930328492646.
}
\]

### Nonempty static gate
\[
\boxed{
|\epsilon\,\delta\ln\mathfrak N_*|
\lesssim 0.737619063660757.
}
\]

Under the additional operator-rigidity closure `D_{01}=0`, the same become

### Robust internal-transfer gate
\[
\boxed{
\left|\epsilon\,\frac{N_{01}}{N_0}\right|
\lesssim 0.367930328492646.
}
\]

### Nonempty internal-transfer gate
\[
\boxed{
\left|\epsilon\,\frac{N_{01}}{N_0}\right|
\lesssim 0.737619063660757.
}
\]

So the Stage-007 `36.8%` figure really is the first exact reduced scalar at which the surviving static same-charge branch stops being robust.

---

## 7. Why this is compatible with the support/source verdict

Stage 015 already showed that the explicit `5PN` support/source branch overshoots the canonical isotropic demand by a large margin and is not the active bottleneck. What is still missing is the actual PDE-selected orbit-lock / coherent placement point.

So if one adopts the rigid-mouth reading, the situation is:

- support/source can feed the branch,
- the dynamic wall-window is not the first failure,
- the first unresolved gate is whether the **internal transfer / placement observable** `\mathfrak N_*` remains inside the Stage-007 static ceiling.

That is exactly the mathematical form of the “internal choke versus rigid entrance” picture.

---

## 8. Best current verdict after Stage 016

The user’s physical interpretation is **mostly compatible** with the reduced math, but with two important corrections:

1. **Good match:** it is reasonable to read the surviving same-charge problem as a rigid-mouth / internal-load competition, because on the exact observable compiler
   \[
   \delta\ln R_{\rm tr}=0
   \quad\Longrightarrow\quad
   \Xi_1=\delta\ln\mathfrak N_*,
   \]
   so the unresolved scalar is literally an internal transfer/placement drift.

2. **Important correction:** the statement
   \[
   D_{01}=0
   \]
   is **not** the same thing as “the entrance radius cannot change.”
   It is a stronger effective-static-operator rigidity condition.
   It is useful for the reduced thought experiment, but it should not be confused with the geometric mouth lock itself.

3. **Interpretive caution:** the Stage-007 `36.8%` ceiling is an exact reduced static branch bound. Calling it a “turbulence threshold” is a plausible physical interpretation, but not yet a derivation of literal hydrodynamic turbulence from the full PDE.

So the next real falsification step is:

> compute or model the actual rigid-mouth orbit packet
> \[
> \bigl(\delta\ln R_{\rm tr},\delta\ln\mathfrak N_*,\delta\ln\epsilon_\eta\bigr)
> \]
> strongly enough to decide whether the realized branch clears or exceeds the static gate above.

That is where the present stack now says the answer lives.

same_charge_barrier_audit_stage017_direct_branch_observable_static_gate.md


# Same-Charge Barrier Audit — Stage 017: Direct Branch-Observable Static Gate and the Two-Observable Kill Test

## Status

**Exact within the carried Stage-016 same-charge reduction and the later `5PN` direct-branch observable compiler.**

This stage does not introduce a new physical mechanism.
It rewrites the unresolved rigid-mouth orbit-lock problem in the sharpest actual-branch variables currently available:

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta.
\]

The main new result is stronger than the Stage-016 wording:

> at first weak-axisymmetric order, the surviving static same-charge gate is already independent of the tracking observable \(R_{\rm tr}\). Once the coherent branch is written in direct observables, the remaining static bottleneck lives entirely on the two-observable plane \((R_{\rm target},\epsilon_\eta)\).

So the “rigid mouth” picture remains useful as an interpretation, but the exact reduced kill test is now a **target–dressing mismatch theorem**, not a direct mouth-motion theorem.

---

## 0. Why this stage is needed

Stage 016 used the Stage-171 observable compiler to show that on a track-locked branch
\[
\delta\ln R_{\rm tr}=0
\quad\Longrightarrow\quad
\Xi_1=\delta\ln \mathfrak N_*,
\qquad
\mathfrak N_*=\mathcal T^2 R_{\rm tr}^{B_*}.
\]

That was already a useful rigid-mouth translation. But the later `5PN` notes go further.
They show that the coherent weak-axisymmetric orbit packet can be charted directly by the exact branch observables
\[
(R_{\rm tr},R_{\rm target},\epsilon_\eta),
\]
and that in this chart the first-order defect packet becomes completely triangular.

So the next honest step is to eliminate the remaining quotient bookkeeping and write the same-charge static gate directly in those observables.

---

## 1. Exact finite quotient chart from direct branch observables

Relative to a coherent reference branch
\[
(R_{\rm tr,ref},R_{\rm target,ref},\epsilon_{\eta,\rm ref}),
\]
the exact finite quotient coordinates are

\[
q_{\rm tr}
=
-\,C_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right),
\]

\[
q_{\rm nt}
=
B_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right)
+
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\]

\[
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\]

with the carried branch constants

\[
C_*=
\frac{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
{\chi_{0,*}\delta_{U,*}},
\qquad
B_*=
\frac{2(1+\chi_{0,*}+\delta_{U,*})}{\delta_{U,*}}.
\]

The inverse map is exact:

\[
R_{\rm tr}
=
R_{\rm tr,ref}\,e^{-q_{\rm tr}/C_*},
\]

\[
\epsilon_\eta
=
\epsilon_{\eta,\rm ref}\,e^{q_\eta},
\]

\[
R_{\rm target}
=
R_{\rm target,ref}\,
e^{-q_{\rm nt}-(B_*/C_*)q_{\rm tr}}
\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}.
\]

So once the actual weak-axisymmetric branch is known, it can be tested either by the finite quotient packet \((q_{\rm tr},q_{\rm nt},q_\eta)\) or directly by the three physical observables \((R_{\rm tr},R_{\rm target},\epsilon_\eta)\).

---

## 2. Exact first-order defect compiler in direct branch language

Linearizing the finite chart gives

\[
q_{\rm tr}
=
-\,C_*\,\delta\ln R_{\rm tr},
\]

\[
q_{\rm nt}
=
B_*\,\delta\ln R_{\rm tr}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln\epsilon_\eta
-
\delta\ln R_{\rm target},
\]

\[
q_\eta
=
\delta\ln\epsilon_\eta.
\]

Composing this with the exact quotient-to-defect compiler yields the triangular first-order direct-branch map

\[
\boxed{\Theta_1=\delta\ln R_{\rm tr},}
\]

\[
\boxed{
\Xi_1
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln\epsilon_\eta,
}
\]

\[
\boxed{\mathcal R_1=\delta\ln R_{\rm target}.}
\]

So the physical first-order defect packet is exactly equivalent to the three direct observable drifts
\[
(\delta\ln R_{\rm tr},\delta\ln R_{\rm target},\delta\ln\epsilon_\eta).
\]

The inverse map is equally simple:

\[
\delta\ln R_{\rm tr} = \Theta_1,
\]

\[
\delta\ln R_{\rm target} = \mathcal R_1,
\]

\[
\delta\ln \epsilon_\eta
=
-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}\,
(\mathcal R_1+\Xi_1).
\]

---

## 3. Two exact cancellations

Define the selected-branch dressing coefficient

\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]

Then
\[
\Xi_1 = -\,\delta\ln R_{\rm target} - c_\eta\,\delta\ln\epsilon_\eta.
\]

This yields two exact structural cancellations.

### 3.1 `R_tr` cancels out of `Xi_1`

Although the finite nontracking quotient contains the tracking observable,
\[
q_{\rm nt}
=
B_*\,\delta\ln R_{\rm tr}
-
c_\eta\,\delta\ln\epsilon_\eta
-
\delta\ln R_{\rm target},
\]
the physical defect is
\[
\Xi_1
=
q_{\rm nt}
+
\frac{B_*}{C_*}\,q_{\rm tr}.
\]

So the tracking contribution cancels exactly, and
\[
\boxed{\partial_{\delta\ln R_{\rm tr}}\Xi_1=0.}
\]

This is stronger than the Stage-016 wording.
It means the surviving first-order static same-charge gate is **already mouth-blind** in the direct observable chart.

### 3.2 Support/source rescaling cancels out of the entire first-order defect packet

The later `5PN` notes also show that the direct branch observables \((R_{\rm tr},R_{\rm target},\epsilon_\eta)\), and therefore the defect packet \((\Theta_1,\Xi_1,\mathcal R_1)\), are exactly blind to the total support-compensation baseline and to the coherent support-enhancement ratio:
\[
\partial_{M_{\rm tr}}(\Theta_1,\Xi_1,\mathcal R_1)=0,
\qquad
\partial_\zeta(\Theta_1,\Xi_1,\mathcal R_1)=0.
\]

So the Stage-015 support/source safety verdict and the Stage-016 static gate really are distinct theorem pieces.

---

## 4. Exact rigid-mouth specialization revisited

The rigid-mouth interpretation still has one exact role:
\[
\delta\ln R_{\rm tr}=0
\quad\Longleftrightarrow\quad
\Theta_1=0.
\]

On that branch,
\[
q_{\rm tr}=0,
\qquad
q_{\rm nt}
=
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right).
\]

And at first order,
\[
\boxed{q_{\rm nt}=\Xi_1.}
\]

So the Stage-016 statement
\[
\Xi_1=\delta\ln\mathfrak N_*
\]
on a track-locked branch sharpens here to:

> under rigid-mouth lock, the surviving first-order defect is exactly the tangent of the finite **two-observable** nontracking quotient built from \(R_{\rm target}\) and \(\epsilon_\eta\).

That is the cleanest current form of the “inside can still repackage while the mouth stays locked” picture.

---

## 5. The direct two-observable static gate

Let
\[
R_1:=\delta\ln R_{\rm target},
\qquad
E_1:=\delta\ln\epsilon_\eta.
\]

Then the entire first-order same-charge placement scalar is

\[
\boxed{\Xi_1 = -R_1 - c_\eta E_1.}
\]

So the Stage-007 / Stage-016 static ceilings become exact two-observable band conditions.

### 5.1 Robust static gate

\[
\boxed{
|\epsilon(-R_1-c_\eta E_1)|
\le
0.367930328492646.
}
\]

Equivalently,
\[
R_1
\in
\left[
-c_\eta E_1-\frac{0.367930328492646}{|\epsilon|},
\,
-c_\eta E_1+\frac{0.367930328492646}{|\epsilon|}
\right].
\]

### 5.2 Nonempty-corridor gate

\[
\boxed{
|\epsilon(-R_1-c_\eta E_1)|
\le
0.737619063660757.
}
\]

Equivalently,
\[
R_1
\in
\left[
-c_\eta E_1-\frac{0.737619063660757}{|\epsilon|},
\,
-c_\eta E_1+\frac{0.737619063660757}{|\epsilon|}
\right].
\]

So the first unresolved same-charge kill test is now literally a strip in the \((R_1,E_1)\) plane.

---

## 6. Canonical direct-branch families

For a prescribed first-order defect value
\[
\Xi_1=\xi,
\]
the direct branch families are simple.

### 6.1 Pure target-drift family

\[
E_1=0,
\qquad
R_1=-\xi.
\]

So if dressing is frozen, the whole defect is a pure selected-target drift.

### 6.2 Pure dressing-drift family

\[
R_1=0,
\qquad
E_1=-\frac{\xi}{c_\eta}.
\]

So if the target is frozen, the whole defect is a pure dressing drift.

### 6.3 Balanced minimal-norm family

Minimizing
\[
R_1^2+E_1^2
\]
subject to
\[
-R_1-c_\eta E_1=\xi
\]
gives the exact minimum-norm branch

\[
\boxed{
R_1=-\frac{\xi}{1+c_\eta^2},
\qquad
E_1=-\frac{c_\eta\,\xi}{1+c_\eta^2}.
}
\]

So even before the full PDE branch is known, the direct observable chart already gives a canonical least-deformation family that realizes any chosen \(\Xi_1\).

---

## 7. What changes in the physical interpretation

The user’s “hyper-trumpet / internal choke” picture is still a useful heuristic, but Stage 17 sharpens it in an important way.

Stage 16 suggested:
- rigid mouth \(\Rightarrow\) internal transfer observable \(\mathfrak N_*\) carries the first unresolved scalar.

Stage 17 now shows:
- at first order, the same unresolved scalar is already
  \[
  \Xi_1=-\delta\ln R_{\rm target}-c_\eta\,\delta\ln\epsilon_\eta,
  \]
  so it is **independent of the mouth tracking variable \(R_{\rm tr}\)**.

That means the first static failure mode is not a mouth-motion theorem even in disguised form.
It is a **selected-target / dressing mismatch theorem**.

A “hyper-trumpet choke” reading is still consistent with the math if the physical reason \(R_{\rm target}\) or \(\epsilon_\eta\) drifts is internal throat repackaging.
But the reduced observable that actually decides the branch is not the mouth variable. It is the two-observable combination above.

---

## 8. Best current verdict after Stage 017

The same-charge corridor is still alive, but the unresolved static bottleneck has narrowed again.

It is no longer:

- support/source strength,
- dynamic wall-window survival,
- or even mouth tracking at first order.

It is now exactly:

\[
\boxed{
\Xi_1
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,
\delta\ln\epsilon_\eta,
}
\]
together with the Stage-016 robust / nonempty ceilings.

So the next honest theorem gate is:

> compute the actual coherent weak-axisymmetric branch drifts of \(R_{\rm target}\) and \(\epsilon_\eta\), then test whether their exact linear combination clears the direct static gate above.

That is the sharpest branch-level same-charge kill test reached so far.

same_charge_barrier_audit_stage018_rigid_mouth_packet_projectors.md

# Same-Charge Barrier Audit — Stage 018: Rigid-Mouth Packet Projectors, the Static-Blind Dressing Line, and the Codimension-Two Orbit-Lock Point

## Status

**Exact within the carried Stage-017 direct branch-observable compiler and the later `5PN` orbit/quotient projector calculus.**

This stage does not add a new physical mechanism.
It upgrades the Stage-017 rigid-mouth strip picture into an exact **two-coordinate projector calculus** on the direct observable plane
\[
(R_1,E_1)
:=
(\delta\ln R_{\rm target},\,\delta\ln\epsilon_\eta).
\]

The main new result is sharper than the Stage-017 wording:

> on the rigid-mouth branch, the first static same-charge scalar sees only one quotient coordinate, `q_nt = Xi_1`, while a second exact dressing coordinate `q_eta` survives completely invisible to that static gate. So the static strip is a codimension-one test inside a codimension-two orbit-lock problem.

In other words, clearing the first static ceiling is necessary but not sufficient even after the mouth side is frozen.

---

## 0. Why this stage is needed

Stage 017 already proved that on the direct coherent branch
\[
\Theta_1=\delta\ln R_{\rm tr},
\qquad
\Xi_1=-\delta\ln R_{\rm target}-c_\eta\,\delta\ln\epsilon_\eta,
\qquad
\mathcal R_1=\delta\ln R_{\rm target},
\]
with
\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]
So under rigid-mouth lock
\[
\delta\ln R_{\rm tr}=0,
\]
the surviving static scalar is already
\[
\Xi_1=-R_1-c_\eta E_1.
\]

But the later orbit/quotient program also shows that full orbit lock is an exact quotient-packet statement: the quotient packet must vanish, not just one scalar projection of it. In the full microscopic language, orbit lock is exactly the zero-set condition
\[
Q_{\rm quot}\,\Delta x=0
\iff
M_*\Delta x=0,
\]
where the finite Packet-B coordinates are `q = (q_tr,q_nt,q_eta)`. So the right next step is to reduce that projector picture onto the rigid-mouth slice and see what the static strip is actually testing. This is exactly what Stage 175 isolates at the microscopic level. 

---

## 1. Exact rigid-mouth direct packet

On the rigid-mouth branch define the direct observable vector
\[
\boxed{
\mathbf x_{\rm rm}:=
\begin{pmatrix}
R_1\\
E_1
\end{pmatrix}
=
\begin{pmatrix}
\delta\ln R_{\rm target}\\
\delta\ln\epsilon_\eta
\end{pmatrix}.
}
\]
The surviving quotient packet is
\[
\boxed{
\mathbf q_{\rm rm}:=
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
\begin{pmatrix}
\Xi_1\\
E_1
\end{pmatrix}.
}
\]
Using the Stage-017 direct compiler,
\[
q_{\rm nt}=-R_1-c_\eta E_1,
\qquad
q_\eta=E_1,
\]
so
\[
\boxed{
\mathbf q_{\rm rm}=M_{\rm rm}\,\mathbf x_{\rm rm},
\qquad
M_{\rm rm}:=
\begin{pmatrix}
-1 & -c_\eta\\
0 & 1
\end{pmatrix}.
}
\]
This matrix is an involution:
\[
\boxed{M_{\rm rm}^2=I_2.}
\]
So the direct observable plane and the rigid-mouth quotient packet are related by an exact one-step involutive compiler.

The inverse map is therefore identical:
\[
\boxed{
\mathbf x_{\rm rm}=M_{\rm rm}\,\mathbf q_{\rm rm},
}
\]
that is,
\[
\boxed{R_1=-q_{\rm nt}-c_\eta q_\eta,}
qquad
\boxed{E_1=q_\eta.}
\]

---

## 2. Exact canonical packet projectors on the rigid-mouth plane

In packet space, the two obvious complementary projectors are
\[
Q_{\rm nt}:=
\begin{pmatrix}
1&0\\
0&0
\end{pmatrix},
\qquad
Q_\eta:=
\begin{pmatrix}
0&0\\
0&1
\end{pmatrix},
\qquad
Q_{\rm nt}+Q_\eta=I_2.
\]
Push them back to direct-observable space with the exact section `S_rm = M_rm`:
\[
\boxed{P_{\rm nt}:=S_{\rm rm}Q_{\rm nt}M_{\rm rm}},
\qquad
\boxed{P_\eta:=S_{\rm rm}Q_\eta M_{\rm rm}}.
\]
Because `M_rm` is its own inverse, these are explicit:
\[
\boxed{
P_{\rm nt}=
\begin{pmatrix}
1 & c_\eta\\
0 & 0
\end{pmatrix},
\qquad
P_\eta=
\begin{pmatrix}
0 & -c_\eta\\
0 & 1
\end{pmatrix}.
}
\]
They are exact complementary projectors:
\[
P_{\rm nt}^2=P_{\rm nt},
\qquad
P_\eta^2=P_\eta,
\qquad
P_{\rm nt}P_\eta=P_\eta P_{\rm nt}=0,
\qquad
P_{\rm nt}+P_\eta=I_2.
\]

So every rigid-mouth direct branch point splits uniquely as
\[
\boxed{
\mathbf x_{\rm rm}=\mathbf x_{\rm nt}+\mathbf x_\eta,
\qquad
\mathbf x_{\rm nt}:=P_{\rm nt}\mathbf x_{\rm rm},
\qquad
\mathbf x_\eta:=P_\eta\mathbf x_{\rm rm}.
}
\]
Explicitly,
\[
\boxed{
\mathbf x_{\rm nt}=
\begin{pmatrix}
R_1+c_\eta E_1\\
0
\end{pmatrix}
=
\begin{pmatrix}
-\Xi_1\\
0
\end{pmatrix},
}
\]
\[
\boxed{
\mathbf x_\eta=
\begin{pmatrix}
-c_\eta E_1\\
E_1
\end{pmatrix}
=
\begin{pmatrix}
-c_\eta q_\eta\\
q_\eta
\end{pmatrix}.
}
\]
So the rigid-mouth direct plane already contains two exact, complementary, physically meaningful pieces:

- `x_nt`: the piece seen by the first static scalar,
- `x_eta`: the pure dressing/selected-target drift that the static scalar cannot see.

---

## 3. Exact codimension-two orbit-lock theorem on the rigid-mouth branch

Because the rigid-mouth quotient packet is
\[
\mathbf q_{\rm rm}=(q_{\rm nt},q_\eta)^T=(\Xi_1,E_1)^T,
\]
full rigid-mouth orbit lock is exactly
\[
\boxed{
\mathbf q_{\rm rm}=0
\iff
q_{\rm nt}=0\ \text{and}\ q_\eta=0.
}
\]
Equivalently,
\[
\boxed{
R_1=0,
\qquad
E_1=0.
}
\]
Using the direct compiler again,
\[
\boxed{
\mathbf q_{\rm rm}=0
\iff
\Xi_1=0\ \text{and}\ R_1=0.
}
\]
Since `q_eta = E_1`, the condition `Xi_1 = 0` alone is not enough.

So on the rigid-mouth branch:

- the static strip `Xi_1 = 0` is codimension one,
- the true orbit-lock point is codimension two.

This is the sharpest version yet of “the static gate is not the whole orbit-lock problem.”

---

## 4. The static-blind dressing line

The compensated strip from Stage 017 is now simply the direct-space image of the packet projector `Q_eta`:
\[
\Xi_1=0
\iff
q_{\rm nt}=0
\iff
\mathbf x_{\rm rm}=\mathbf x_\eta=
\begin{pmatrix}
-c_\eta q_\eta\\
q_\eta
\end{pmatrix}.
\]
So the entire static-blind line is parameterized by the single dressing coordinate `q_eta`.

Its direct-space norm is exact:
\[
\boxed{
\|\mathbf x_\eta\|^2=(1+c_\eta^2)\,q_\eta^2.
}
\]
Therefore the static strip contains points arbitrarily far from the orbit-lock point.
For any prescribed size `L > 0`, choose
\[
q_\eta = \frac{L}{\sqrt{1+c_\eta^2}},
\]
then
\[
\Xi_1=0,
\qquad
\|\mathbf x_{\rm rm}\|=L.
\]

So the first static same-charge gate does **not** bound the true rigid-mouth orbit-failure amplitude.
It only bounds the `q_nt` component.

This is the exact static-blindness theorem on the rigid-mouth slice.

---

## 5. Canonical correction compilers

The projector formulas immediately give the two natural rigid-mouth corrections.

### 5.1 Static-only correction

Project to the static strip by removing only the `q_nt` component:
\[
\boxed{
\Delta\mathbf x_{\rm static}:=-\mathbf x_{\rm nt}
=
\begin{pmatrix}
\Xi_1\\
0
\end{pmatrix}.
}
\]
After this correction,
\[
\mathbf x_{\rm rm}+\Delta\mathbf x_{\rm static}=\mathbf x_\eta,
\qquad
\Xi_1\to 0,
\qquad
q_\eta\to q_\eta.
\]
So the branch clears the first static ceiling but generally still fails orbit lock.

### 5.2 Full orbit-lock correction

Project all the way to the orbit-lock point by removing both packet components:
\[
\boxed{
\Delta\mathbf x_{\rm orbit}:=-\mathbf x_{\rm rm}
=
\begin{pmatrix}
q_{\rm nt}+c_\eta q_\eta\\
-q_\eta
\end{pmatrix}.
}
\]
This is exactly the sum of the static correction and the static-blind dressing correction:
\[
\boxed{
\Delta\mathbf x_{\rm orbit}
=
\Delta\mathbf x_{\rm static}+
\begin{pmatrix}
 c_\eta q_\eta\\
-q_\eta
\end{pmatrix}.
}
\]
So the extra step beyond the static gate is not mysterious. It is the exact removal of the packet component `q_eta`.

---

## 6. What changes physically after Stage 018

Stage 017 already showed that first-order same-charge survival on the rigid-mouth branch is governed by the strip
\[
|\epsilon\Xi_1|\le B,
\qquad
\Xi_1=-R_1-c_\eta E_1.
\]
Stage 018 now sharpens that into an exact orbit-packet statement:

1. the static gate constrains only `q_nt = Xi_1`,
2. the dressing coordinate `q_eta = E_1` survives completely outside that gate,
3. and the full rigid-mouth orbit-lock problem is exactly the vanishing of both packet coordinates.

So the next honest theorem gate is no longer “does the branch clear the static strip?”
It is:

> compute the actual rigid-mouth dressing coordinate `q_eta = \delta\ln\epsilon_\eta` (equivalently `R_1` once `Xi_1` is known), because that is the exact static-blind residue that still blocks orbit lock after the first static same-charge ceiling is cleared.

same_charge_barrier_audit_stage019_microscopic_dependent_plane_projectors.md

# Same-Charge Barrier Audit — Stage 019: Rigid-Mouth Microscopic Dependent-Plane Projectors, the Equal-Drift Dressing Ray, and the Static-Only Restoration Gap

## Status

**Exact within the carried Stage-018 rigid-mouth packet split and the Stage-175/176 microscopic quotient section on the dependent triple.**

This stage does not yet compute the actual magnitude of the rigid-mouth dressing coordinate on the completed PDE branch.
It does something narrower and sharper:

> it identifies the exact **microscopic carrier** of the static-blind residue `q_eta` once the first static same-charge strip has already been cleared.

The main new result is that on the rigid-mouth slice the surviving dressing residue is not a generic three-coordinate failure.
It lives on the exact diagonal ray
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
=
-q_\eta\,(0,1,1),
\]
inside the dependent triple.
So after the static gate is cleared, the unresolved same-charge orbit defect is exactly an **equal-drift `K_\eta`–`\mu` dressing shift at fixed `T_U`**.

---

## 0. Why this stage is needed

Stage 018 showed that on the rigid-mouth direct-observable plane
\[
(R_1,E_1):=(\delta\ln R_{\rm target},\,\delta\ln\epsilon_\eta),
\]
the surviving quotient packet is
\[
(q_{\rm nt},q_\eta)=(\Xi_1,E_1),
\]
and that the first static same-charge gate only constrains
\[
q_{\rm nt}=\Xi_1,
\]
while the dressing coordinate
\[
q_\eta=\delta\ln\epsilon_\eta
\]
remains completely invisible to that gate.

But the later quotient-projector program already says that every quotient failure is represented microscopically only on the dependent triple
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu),
\]
via the exact section `S_{(T,K_\eta,\mu)}`. So the right next step is now obvious:

> what is the exact microscopic dependent-coordinate image of the rigid-mouth packet `(q_{\rm nt},q_\eta)`?

This stage answers that exactly.

---

## 1. Carry-forward rigid-mouth packet and dependent section

From Stage 018 the rigid-mouth packet map is
\[
\boxed{
\mathbf q_{\rm rm}:=
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_{\rm rm}
\begin{pmatrix}
R_1\\
E_1
\end{pmatrix},
\qquad
M_{\rm rm}=
\begin{pmatrix}
-1 & -c_\eta\\
0 & 1
\end{pmatrix},
}
\]
with
\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]
So explicitly,
\[
q_{\rm nt}=-R_1-c_\eta E_1=\Xi_1,
\qquad
q_\eta=E_1.
\]

From the exact microscopic quotient section on the dependent triple,
\[
\Delta_T^{(q)}=\frac{q_{\rm tr}}{1+\chi_{0,*}},
\qquad
\Delta_{K_\eta}^{(q)}=-q_\eta,
\qquad
\Delta_\mu^{(q)}=\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}+q_{\rm nt}-q_\eta.
\]
On the rigid-mouth slice,
\[
q_{\rm tr}=0,
\]
so the dependent correction reduces exactly to
\[
\boxed{
\mathbf y_{\rm rm}:=
\begin{pmatrix}
\Delta_T^{(q)}\\
\Delta_{K_\eta}^{(q)}\\
\Delta_\mu^{(q)}
\end{pmatrix}
=
S_{\rm rm}^{\rm dep}
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix},
\qquad
S_{\rm rm}^{\rm dep}:=
\begin{pmatrix}
0 & 0\\
0 & -1\\
1 & -1
\end{pmatrix}.
}
\]
So explicitly,
\[
\boxed{
\Delta_T^{(q)}=0,
\qquad
\Delta_{K_\eta}^{(q)}=-q_\eta,
\qquad
\Delta_\mu^{(q)}=q_{\rm nt}-q_\eta.
}
\]

The rigid-mouth quotient failure therefore fills the full plane `\Delta_T=0` inside the dependent triple.

---

## 2. Exact direct-observable-to-dependent compiler

Composing the direct rigid-mouth packet map with the dependent microscopic section gives
\[
\boxed{
\mathbf y_{\rm rm}
=
C_{\rm rm}^{\rm dep}
\begin{pmatrix}
R_1\\
E_1
\end{pmatrix},
\qquad
C_{\rm rm}^{\rm dep}:=S_{\rm rm}^{\rm dep}M_{\rm rm}
=
\begin{pmatrix}
0 & 0\\
0 & -1\\
-1 & -\dfrac{1}{1-\epsilon_{\eta,*}}
\end{pmatrix}.
}
\]
So the exact microscopic dependent correction in direct rigid-mouth variables is
\[
\boxed{
\Delta_T^{(q)}=0,
}
\]
\[
\boxed{
\Delta_{K_\eta}^{(q)}=-E_1=-\delta\ln\epsilon_\eta,
}
\]
\[
\boxed{
\Delta_\mu^{(q)}=-R_1-\frac{1}{1-\epsilon_{\eta,*}}E_1
=-\delta\ln R_{\rm target}-\frac{1}{1-\epsilon_{\eta,*}}\delta\ln\epsilon_\eta.
}
\]

So even before any projector algebra is applied, one exact fact is already clear:

> on the rigid-mouth slice, the microscopic quotient correction never uses `T_U` at first order.
> The entire failure is carried only by the `(K_\eta,\mu)` plane.

---

## 3. Exact dependent-plane packet projectors

The rigid-mouth dependent compiler has an exact left inverse on the plane `\Delta_T=0`:
\[
\boxed{
L_{\rm rm}^{\rm dep}:=
\begin{pmatrix}
0 & -1 & 1\\
0 & -1 & 0
\end{pmatrix},
\qquad
L_{\rm rm}^{\rm dep}S_{\rm rm}^{\rm dep}=I_2.
}
\]
So the packet coordinates can be recovered directly from the dependent correction by
\[
\boxed{
q_{\rm nt}=\Delta_\mu^{(q)}-\Delta_{K_\eta}^{(q)},
\qquad
q_\eta=-\Delta_{K_\eta}^{(q)}.
}
\]

Push the packet projectors back to the dependent plane:
\[
P_{\rm nt}^{\rm dep}:=S_{\rm rm}^{\rm dep}
\begin{pmatrix}1&0\\0&0\end{pmatrix}
L_{\rm rm}^{\rm dep},
\qquad
P_\eta^{\rm dep}:=S_{\rm rm}^{\rm dep}
\begin{pmatrix}0&0\\0&1\end{pmatrix}
L_{\rm rm}^{\rm dep}.
\]
They are explicit:
\[
\boxed{
P_{\rm nt}^{\rm dep}=
\begin{pmatrix}
0&0&0\\
0&0&0\\
0&-1&1
\end{pmatrix},
\qquad
P_\eta^{\rm dep}=
\begin{pmatrix}
0&0&0\\
0&1&0\\
0&1&0
\end{pmatrix}.
}
\]
These are exact complementary projectors on the rigid-mouth dependent plane:
\[
(P_{\rm nt}^{\rm dep})^2=P_{\rm nt}^{\rm dep},
\qquad
(P_\eta^{\rm dep})^2=P_\eta^{\rm dep},
\qquad
P_{\rm nt}^{\rm dep}P_\eta^{\rm dep}=P_\eta^{\rm dep}P_{\rm nt}^{\rm dep}=0,
\]
and
\[
P_{\rm nt}^{\rm dep}+P_\eta^{\rm dep}=
\begin{pmatrix}
0&0&0\\
0&1&0\\
0&0&1
\end{pmatrix},
\]
which is just the identity on the plane `\Delta_T=0`.

So every rigid-mouth dependent correction decomposes uniquely as
\[
\boxed{
\mathbf y_{\rm rm}=\mathbf y_{\rm nt}+\mathbf y_\eta,
\qquad
\mathbf y_{\rm nt}:=P_{\rm nt}^{\rm dep}\mathbf y_{\rm rm},
\qquad
\mathbf y_\eta:=P_\eta^{\rm dep}\mathbf y_{\rm rm}.
}
\]
Explicitly,
\[
\boxed{
\mathbf y_{\rm nt}=
\begin{pmatrix}
0\\0\\q_{\rm nt}
\end{pmatrix},
\qquad
\mathbf y_\eta=
-q_\eta
\begin{pmatrix}
0\\1\\1
\end{pmatrix}.
}
\]

This is the microscopic dependent-plane analogue of the Stage-018 direct-space projector split.

---

## 4. The equal-drift dressing ray and the microscopic meaning of the static strip

Stage 018 showed that the static strip is
\[
q_{\rm nt}=\Xi_1=0.
\]
In the dependent microscopic plane, this becomes
\[
\boxed{
q_{\rm nt}=0
\iff
\Delta_\mu^{(q)}=\Delta_{K_\eta}^{(q)}.
}
\]
So the entire static-blind line is exactly the diagonal ray
\[
\boxed{
\mathbf y_\eta=-q_\eta(0,1,1)^T.
}
\]
Equivalently, using the direct observable compiler,
\[
q_{\rm nt}=0
\iff
R_1=-c_\eta E_1,
\]
and then
\[
\boxed{
\mathbf y_{\rm rm}=
\begin{pmatrix}
0\\-E_1\\-E_1
\end{pmatrix}
=
\frac{R_1}{c_\eta}
\begin{pmatrix}
0\\1\\1
\end{pmatrix}.
}
\]
So the Stage-018 static-blind line in the direct observable plane maps exactly to an **equal-drift `K_\eta`–`\mu` ray** in the dependent microscopic plane.

Its microscopic norm is exact:
\[
\boxed{
\|\mathbf y_\eta\|^2 = 2 q_\eta^2 = 2 E_1^2 = \frac{2R_1^2}{c_\eta^2}.
}
\]
Therefore clearing the first static same-charge ceiling does not force the microscopic quotient correction to be small.
It only forces that correction to lie on a one-dimensional diagonal ray.

The full rigid-mouth orbit-lock point remains the endpoint of that ray:
\[
\boxed{
q_{\rm nt}=0,
\ q_\eta=0
\iff
\Delta_T^{(q)}=\Delta_{K_\eta}^{(q)}=\Delta_\mu^{(q)}=0.
}
\]

---

## 5. Exact microscopic correction compilers

The dependent-plane projectors immediately give the two natural microscopic corrections.

### 5.1 Static-only microscopic correction

Remove only the static component:
\[
\boxed{
\Delta\mathbf y_{\rm static}:=-\mathbf y_{\rm nt}=
\begin{pmatrix}
0\\0\\-q_{\rm nt}
\end{pmatrix}.
}
\]
After this correction,
\[
\mathbf y_{\rm rm}+\Delta\mathbf y_{\rm static}=\mathbf y_\eta,
\qquad
q_{\rm nt}\to 0,
\qquad
q_\eta\to q_\eta.
\]
So the first static ceiling is cleared by changing only `\mu_W` inside the dependent triple.

### 5.2 Full orbit-lock correction

Remove the entire rigid-mouth dependent correction:
\[
\boxed{
\Delta\mathbf y_{\rm orbit}:=-\mathbf y_{\rm rm}.
}
\]
Equivalently,
\[
\boxed{
\Delta\mathbf y_{\rm orbit}
=
\Delta\mathbf y_{\rm static}
+
q_\eta
\begin{pmatrix}
0\\1\\1
\end{pmatrix}.
}
\]
So the extra step beyond the static gate is again completely sharp:

> once the `\mu` mismatch `q_{\rm nt}` has been removed, the only remaining orbit-restoring correction is the equal `K_\eta`–`\mu` dressing shift generated by `q_\eta`.

---

## 6. What changes physically after Stage 019

Stage 018 already said that the first static same-charge strip is not the whole rigid-mouth orbit-lock problem.
Stage 019 now says something more microscopic and more useful:

1. on the rigid-mouth slice, the quotient-failure image is the full plane `\Delta_T=0`,
2. the first static gate tests only the `\mu-K_\eta` difference,
3. the static-blind residue is exactly the diagonal equal-drift ray `\Delta_\mu=\Delta_{K_\eta}`,
4. and the surviving same-charge obstruction after the static strip is cleared is therefore not generic throat motion but one scalar `K_\eta`–`\mu` dressing amplitude.

So the next honest theorem gate is now even sharper than at Stage 018:

> compute the actual dressing coordinate `q_\eta = \delta\ln\epsilon_\eta`, because after the first static gate is cleared that single scalar is exactly the amplitude of the remaining equal-drift microscopic obstruction.

same_charge_barrier_audit_stage020_actual_branch_dressing_compiler.md

# Same-Charge Barrier Audit — Stage 020: Actual-Branch Dressing Compiler, the Finite Static-Blind Curve, and the Support-Blind Post-Static Orbit-Lock Theorem

## Status

**Exact within the carried Stage-019 rigid-mouth dependent-plane split and the later direct-branch / microscopic packet compilers on the actual coherent branch.**

This stage does not solve the full moving-throat PDE branch.
It does something narrower and more useful:

> it computes the surviving rigid-mouth dressing coordinate `q_eta` exactly in the actual branch observables and in the actual microscopic variables, and shows that this post-static obstruction is completely blind to the coherent support-enhancement sector.

So after Stage 020, the next same-charge gate is no longer vague.
Once the first static ceiling `q_nt = 0` has been cleared, the remaining obstruction is exactly the single actual-branch scalar
\[
q_\eta = \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\mathrm{ref}}}\right),
\]
which may also be read directly from the target observable `R_target` on the static-blind curve.

---

## 0. Why this stage is needed

Stage 019 proved that on the rigid-mouth slice the remaining quotient failure is carried microscopically by the equal-drift dependent-plane ray
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)
=
-q_\eta(0,1,1).
\]
So after the first static gate is cleared, the entire unresolved same-charge obstruction is one scalar: the dressing coordinate `q_eta`.

But Stage 019 still left one obvious unfinished task:

> how do we compute `q_eta` on the actual coherent branch, rather than only locating its microscopic carrier?

The later branch-packet notes already contain the answer. On the actual coherent branch the finite quotient packet is charted exactly by the three physical observables
\[
(R_{\rm tr},\,R_{\rm target},\,\epsilon_\eta),
\]
and the dressing coordinate is simply the logarithmic ratio of `\epsilon_\eta` itself. This stage isolates that fact, pushes it through the rigid-mouth static gate, and shows that the surviving dressing obstruction is rigorously support-blind.

---

## 1. Exact rigid-mouth finite packet on the actual branch

Relative to a coherent reference branch
\[
(R_{\rm tr,ref},\,R_{\rm target,ref},\,\epsilon_{\eta,\rm ref}),
\]
the exact finite quotient coordinates are
\[
q_{\rm tr}=-C_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right),
\]
\[
q_{\rm nt}=B_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right)
+\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\]
\[
q_\eta=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]

On the rigid-mouth slice,
\[
q_{\rm tr}=0,
\qquad
R_{\rm tr}=R_{\rm tr,ref},
\]
so the surviving packet is exactly
\[
\boxed{
q_{\rm nt}
=
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\qquad
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]

So on the actual rigid-mouth branch the two finite surviving coordinates are already completely explicit:

- `q_nt` is the target–dressing mismatch of the finite branch observables,
- `q_eta` is the pure dressing logarithmic ratio.

---

## 2. Exact finite static-blind curve and direct computation of `q_eta`

The first static same-charge ceiling is
\[
q_{\rm nt}=0.
\]
On the rigid-mouth branch this becomes the exact finite relation
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}.
}
\]
So the static-blind set is not just a tangent line. It is an exact one-parameter finite curve in the direct branch observables.

Parameterizing that curve by `q_eta`, we use
\[
\epsilon_\eta = \epsilon_{\eta,\rm ref}e^{q_\eta},
\]
which gives
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^{q_\eta}}{1-\epsilon_{\eta,\rm ref}}.
}
\]
Conversely, once the static gate has been cleared, `q_eta` can be computed directly from the actual branch target observable alone:
\[
\boxed{
q_\eta
=
\ln\!\left(
\frac{1-(1-\epsilon_{\eta,\rm ref})\,R_{\rm target}/R_{\rm target,ref}}
{\epsilon_{\eta,\rm ref}}
\right).
}
\]
So after the static gate is passed, the direct same-charge branch is exactly one-dimensional and can be charted either by

- the dressing observable `\epsilon_\eta`, or
- the target observable `R_{\rm target}`.

The full rigid-mouth orbit-lock point remains the endpoint of that curve:
\[
q_{\rm nt}=0,
\ q_\eta=0
\iff
R_{\rm target}=R_{\rm target,ref},
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
\]

---

## 3. Exact first-order compiler and the tangent of the finite static-blind curve

Linearizing around the coherent reference branch gives
\[
q_\eta = \delta\ln\epsilon_\eta,
\]
and on the rigid-mouth slice
\[
q_{\rm nt}
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta.
\]
Define the carried dressing coefficient
\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]
Then the exact first-order direct packet is
\[
\boxed{
q_{\rm nt}=-R_1-c_\eta E_1,
\qquad
q_\eta=E_1,
}
\]
with
\[
R_1:=\delta\ln R_{\rm target},
\qquad
E_1:=\delta\ln\epsilon_\eta.
\]

So once the static gate is cleared,
\[
q_{\rm nt}=0
\iff
R_1=-c_\eta q_\eta,
\qquad
q_\eta=-\frac{R_1}{c_\eta}.
\]
This is exactly the tangent relation of the finite curve from Section 2.
Indeed,
\[
\frac{d}{dq_\eta}
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right)\Bigg|_{q_\eta=0}
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
=-c_\eta.
\]
So the Stage-018 static-blind line is precisely the tangent of the finite actual-branch static-blind curve.

---

## 4. Exact microscopic compiler for `q_eta`

### 4.1 Direct microscopic coherent variables

On the actual coherent branch,
\[
\epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}},
\]
so the finite dressing coordinate is
\[
\boxed{
q_\eta
=
\ln\!\left(\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}
\frac{K_{U,\rm ref}K_{\eta,\rm ref}^{(\mathrm{eff})}}{c_{\eta U,\rm ref}^2}
\right).
}
\]
Equivalently,
\[
\boxed{
q_\eta
=
2\ln\!\left(\frac{c_{\eta U}}{c_{\eta U,\rm ref}}\right)
-
\ln\!\left(\frac{K_U}{K_{U,\rm ref}}\right)
-
\ln\!\left(\frac{K_\eta^{(\mathrm{eff})}}{K_{\eta,\rm ref}^{(\mathrm{eff})}}\right).
}
\]

### 4.2 First-order microscopic drift packet

For the weak-axisymmetric microscopic drifts
\[
(c_1,\,\kappa_U,\,\kappa_\eta),
\]
the exact first-order extractor is
\[
\boxed{
q_\eta
=
\delta\ln\epsilon_\eta
=
2c_1-\kappa_U-\kappa_\eta.
}
\]
So the post-static same-charge obstruction is directly the combined logarithmic drift of

- the wall–U coupling `c_{ηU}`,
- the `U` stiffness `K_U`,
- and the effective wall stiffness `K_η^{(eff)}`.

No additional transport bookkeeping is needed.

---

## 5. Exact support-blindness theorem for the dressing coordinate

Introduce the coherent support ratio in the usual form
\[
\zeta
=
\frac{\lambda_\phi^2 K_W^{(\mathrm{eff})}}{\lambda_W^2 K_\phi^{(\mathrm{eff})}},
\]
and the total coherent support-enhanced baseline
\[
M_{\rm tr}=M_{\rm mix}
\left[
1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}
\right].
\]
Then the exact dressing coordinate satisfies
\[
\boxed{
\partial_\zeta q_\eta = 0,
\qquad
\partial_{M_{\rm tr}} q_\eta = 0.
}
\]
At the microscopic parameter level this is even sharper:
\[
\boxed{
\partial_{\lambda_\phi} q_\eta = 0,
\qquad
\partial_{K_\phi^{(\mathrm{eff})}} q_\eta = 0.
}
\]
So the coherent support lane is exactly blind to the dressing obstruction.

This has a strong physical consequence:

> support enhancement can rescue the steady normalization side of the branch, but it cannot change the post-static same-charge dressing obstruction.

So if the static gate has been cleared and `q_eta` is still nonzero, no adjustment of the coherent support ratio can remove that failure.

---

## 6. Post-static orbit-lock theorem on the actual rigid-mouth branch

Combine the results of Stages 019 and 020.

- Stage 019 showed that after the static gate is cleared, the remaining microscopic correction is the equal-drift dependent-plane ray
  \[
  \mathbf y_\eta = -q_\eta(0,1,1)^T.
  \]
- Stage 020 shows that the amplitude of that ray is exactly
  \[
  q_\eta = \ln(\epsilon_\eta/\epsilon_{\eta,\rm ref}),
  \]
  equivalently
  \[
  q_\eta = -R_1/c_\eta
  \quad\text{on the post-static branch.}
  \]

Therefore the full rigid-mouth post-static orbit-lock theorem is
\[
\boxed{
q_{\rm nt}=0\ \text{and}\ q_\eta=0
\iff
R_{\rm target}=R_{\rm target,ref}
\ \text{and}\
\epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]
At first order this becomes
\[
\boxed{
q_{\rm nt}=0\ \text{and}\ q_\eta=0
\iff
R_1=0
\ \text{and}\
E_1=0.
}
\]
So after the first static ceiling is cleared, the same-charge barrier reduces to one exact test:

> is `\epsilon_\eta` invariant on the actual rigid-mouth branch?

If yes, the post-static dressing ray collapses and orbit lock is restored.
If no, the same-charge obstruction remains, and its exact amplitude is `q_eta`.

---

## 7. Best current summary after Stage 020

The continuation from Stage 019 is now complete.

1. The static-blind microscopic carrier is still the equal-drift dependent-plane ray
   \[
   -q_\eta(0,1,1).
   \]
2. Its exact actual-branch amplitude is now explicit:
   \[
   q_\eta = \ln(\epsilon_\eta/\epsilon_{\eta,\rm ref}).
   \]
3. Once the first static gate is cleared, that same amplitude can be read directly from `R_target`:
   \[
   q_\eta
   =
   \ln\!\left(
   \frac{1-(1-\epsilon_{\eta,\rm ref})R_{\rm target}/R_{\rm target,ref}}
   {\epsilon_{\eta,\rm ref}}
   \right).
   \]
4. The dressing coordinate is exactly support-blind.
5. Therefore the next actual same-charge theorem gate is no longer the full quotient packet.
   It is simply:

> compute `\epsilon_\eta` on the actual rigid-mouth branch and test whether it is invariant.

That is the sharpest post-static same-charge criterion reached so far.

same_charge_barrier_audit_stage021_physical_branch_transfer_shape_compiler.md

# Same-Charge Barrier Audit — Stage 021: Physical-Branch Transfer-Shape Compiler, Packet Factorization, and the Post-Static Dressing-Invariance Theorem

## Status

**Exact within the carried Stage-020 actual-branch dressing compiler and the later coherent-branch observable compiler already frozen in the moving-throat notes.**

This stage does not solve the full moving-throat PDE branch.
It does the next sharp reduction after Stage 020:

> it rewrites the actual same-charge packet directly in the physical coherent branch variables 
> \((R_{\rm tr},\,\mathcal T^2,\,\epsilon_\eta)\),
> shows that the static same-charge ceiling is exactly transfer-shape rigidity,
> and isolates the post-static barrier as pure dressing invariance.

So after this stage the same-charge chain factors into three clean gates:

1. tracking rigidity,
2. transfer-shape rigidity,
3. dressing rigidity.

And once the first two have been cleared, the barrier reduces to one scalar test only:
\[
\delta\ln\epsilon_\eta=0.
\]

---

## 0. Why this stage is needed

Stage 020 proved that on the rigid-mouth slice the surviving post-static obstruction is
\[
q_\eta = \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\]
and that this scalar is exactly support-blind.

But Stage 020 still left one useful bridge implicit:

> what does the static same-charge gate `q_nt = 0` look like in the actual coherent branch variables before the dressing scalar is read off?

The later coherent-branch observable compiler already supplies the missing object. On that branch the target observable is not primitive; it is tied to one wall-normalized transfer shape
\[
\mathcal T^2.
\]
Once that identity is inserted, the same-charge packet becomes triangular in the physical variables themselves.

---

## 1. Exact coherent-branch observables and the transfer-shape identity

On the actual coherent branch the direct observables are
\[
R_{\rm tr} = \frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=\frac{1+\chi_0+\delta_U}{(1+\chi_0)(1+\delta_U)},
\]
\[
\mathcal T^2 = \frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2},
\]
\[
R_{\rm target}
=
\Lambda_0\,\frac{\Omega_W^2(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]
with the exact selected-branch identity
\[
\boxed{
R_{\rm target}\,\mathcal T^2 = \Lambda_0(1-\epsilon_\eta).
}
\]
So the actual coherent branch itself supplies the finite packet
\[
(R_{\rm tr},\,\mathcal T^2,\,\epsilon_\eta),
\]
from which the quotient coordinates can be reconstructed exactly.

---

## 2. Exact finite same-charge packet in physical branch variables

Relative to a coherent reference branch
\[
(R_{\rm tr,ref},\,\mathcal T_{\rm ref}^2,\,\epsilon_{\eta,\rm ref}),
\]
the carried finite packet is
\[
q_{\rm tr}=-C_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right),
\]
\[
q_{\rm nt}
=
B_*\ln\!\left(\frac{R_{\rm tr}}{R_{\rm tr,ref}}\right)
+
\ln\!\left(\frac{1-\epsilon_\eta}{1-\epsilon_{\eta,\rm ref}}\right)
-
\ln\!\left(\frac{R_{\rm target}}{R_{\rm target,ref}}\right),
\]
\[
q_\eta = \ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]
Using
\[
R_{\rm target}\,\mathcal T^2 = \Lambda_0(1-\epsilon_\eta),
\qquad
R_{\rm target,ref}\,\mathcal T_{\rm ref}^2 = \Lambda_0(1-\epsilon_{\eta,\rm ref}),
\]
one gets the exact finite factorization
\[
\boxed{
q_{\rm nt} + \frac{B_*}{C_*}q_{\rm tr}
=
\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right).
}
\]
So the finite nontracking packet is already the transfer-shape ratio up to the universal tracking feed-through.

On the rigid-mouth slice,
\[
q_{\rm tr}=0,
\qquad
R_{\rm tr}=R_{\rm tr,ref},
\]
so the finite surviving packet becomes
\[
\boxed{
q_{\rm nt}=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
q_\eta=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]
Therefore the first static same-charge ceiling is exactly
\[
\boxed{
q_{\rm nt}=0
\iff
\mathcal T^2 = \mathcal T_{\rm ref}^2.
}
\]
So the finite static-blind set is transfer-shape rigidity, not an additional hidden quotient condition.

---

## 3. Exact first-order physical drift compiler

Linearizing the coherent observables gives
\[
\delta\ln R_{\rm tr}
=
-
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\bigl[(1+\delta_U)\,d\ln\chi_0 + (1+\chi_0)\,d\ln\delta_U\bigr],
\]
\[
\delta\ln \mathcal T^2
=
 d\ln Z_W - d\ln \Omega_W^2
 + \frac{2\chi_0}{1+\chi_0}d\ln\chi_0
 + \frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
\]
\[
\delta\ln R_{\rm target}
=
 d\ln\Omega_W^2 - d\ln Z_W
 - \frac{2\chi_0}{1+\chi_0}d\ln\chi_0
 - \frac{2\epsilon}{1-\epsilon}d\ln\epsilon
 - \frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta,
\]
\[
\delta\ln\epsilon_\eta = d\ln\epsilon_\eta.
\]
So the first-order same-charge packet is
\[
\boxed{
q_{\rm tr} = -C_*\,\delta\ln R_{\rm tr},
\qquad
q_{\rm nt} + \frac{B_*}{C_*}q_{\rm tr} = \delta\ln\mathcal T^2,
\qquad
q_\eta = d\ln\epsilon_\eta.
}
\]
On the rigid-mouth branch,
\[
q_{\rm tr}=0,
\]
so the surviving first-order packet is simply
\[
\boxed{
q_{\rm nt} = \delta\ln\mathcal T^2,
\qquad
q_\eta = d\ln\epsilon_\eta.
}
\]
Thus the first-order static same-charge ceiling is
\[
\boxed{
q_{\rm nt}=0
\iff
\delta\ln\mathcal T^2=0.
}
\]
And once that ceiling has been cleared, the remaining post-static obstruction is exactly
\[
\boxed{
q_\eta = d\ln\epsilon_\eta.
}
\]

---

## 4. Exact support-blindness factorization of the physical packet

The coherent support-enhancement sector enters only through the baseline factor
\[
M_{\rm tr}=M_{\rm mix}
\left[1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon}\right].
\]
But the direct same-charge observables satisfy
\[
\partial_\zeta R_{\rm tr}=0,
\qquad
\partial_\zeta \mathcal T^2=0,
\qquad
\partial_\zeta \epsilon_\eta=0,
\]
\[
\partial_{M_{\rm mix}} R_{\rm tr}=0,
\qquad
\partial_{M_{\rm mix}} \mathcal T^2=0,
\qquad
\partial_{M_{\rm mix}} \epsilon_\eta=0.
\]
Therefore the full first-order packet is support-blind:
\[
\boxed{
\partial_\zeta q_{\rm tr}=\partial_\zeta q_{\rm nt}=\partial_\zeta q_\eta=0,
\qquad
\partial_{M_{\rm mix}} q_{\rm tr}=\partial_{M_{\rm mix}} q_{\rm nt}=\partial_{M_{\rm mix}} q_\eta=0.
}
\]
So support enhancement may rescue the steady normalization side of the branch, but it cannot change the actual same-charge packet at first weak-axisymmetric order.

---

## 5. Post-static dressing-invariance theorem on the actual branch

The same-charge chain on the actual coherent rigid-mouth branch now factors exactly into three gates.

### 5.1 Tracking gate
\[
q_{\rm tr}=0
\iff
R_{\rm tr}=R_{\rm tr,ref}
\iff
(1+\delta_U)d\ln\chi_0 + (1+\chi_0)d\ln\delta_U = 0
\quad\text{(at first order).}
\]

### 5.2 Static-blind transfer-shape gate
\[
q_{\rm nt}=0
\iff
\mathcal T^2 = \mathcal T_{\rm ref}^2
\iff
\delta\ln\mathcal T^2 = 0
\quad\text{(at first order).}
\]

### 5.3 Post-static dressing gate
After the first two have been cleared, the remaining obstruction is exactly
\[
q_\eta
=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right),
\qquad
q_\eta^{(1)} = d\ln\epsilon_\eta.
\]
So the final post-static criterion is
\[
\boxed{
q_\eta = 0
\iff
\epsilon_\eta = \epsilon_{\eta,\rm ref}
\iff
 d\ln\epsilon_\eta = 0
\quad\text{(at first order).}
}
\]
This is the exact post-static dressing-invariance theorem.

---

## 6. Best current summary after Stage 021

The continuation from Stage 020 is now complete.

1. The actual same-charge packet is charted directly by
   \[
   (R_{\rm tr},\,\mathcal T^2,\,\epsilon_\eta).
   \]
2. The finite static same-charge ceiling is exactly
   \[
   q_{\rm nt}=0 \iff \mathcal T^2 = \mathcal T_{\rm ref}^2.
   \]
3. The first-order static ceiling is exactly
   \[
   q_{\rm nt} = \delta\ln\mathcal T^2.
   \]
4. The whole direct packet is support-blind.
5. Therefore, once tracking rigidity and transfer-shape rigidity have both been imposed, the same-charge barrier reduces to one exact test only:

> compute `\epsilon_\eta` on the actual rigid-mouth branch and check whether it is invariant.

That is the sharpest post-static physical-variable criterion reached so far.

same_charge_barrier_audit_stage022_rigid_mouth_physical_normal_form.md

# Same-Charge Barrier Audit — Stage 022: Rigid-Mouth Physical Normal Form, Exact Physical-to-Microscopic Correction Compiler, and the Cartesian Orbit-Lock Theorem

## Status

**Exact within the carried Stage-019 rigid-mouth dependent-plane projector calculus and the Stage-021 physical-branch transfer-shape compiler.**

This stage does not solve the full moving-throat PDE branch.
It does the next clean reduction after Stage 021:

> it diagonalizes the rigid-mouth same-charge packet in the actual physical logarithmic variables 
> \((\mathcal T^2,\epsilon_\eta)\),
> converts that diagonal packet into an exact microscopic dependent-plane compiler,
> and shows that rigid-mouth orbit lock is already a Cartesian product of transfer-shape rigidity and dressing rigidity.

So after this stage the surviving rigid-mouth same-charge geometry is no longer triangular.
It is exactly a two-axis normal form:

1. a **pure transfer-shape axis**,
2. a **pure dressing axis**.

And the corresponding microscopic correction splits just as sharply:

- clearing the static transfer-shape defect changes only `\mu_W`,
- clearing the post-static dressing defect adds the equal `K_\eta^{(\mathrm{eff})}`–`\mu_W` shift.

---

## 0. Why this stage is needed

Stage 019 already proved that on the rigid-mouth slice the quotient-failure image in the dependent microscopic plane is
\[
\mathbf y_{\rm rm}
=
\mathbf y_{\rm nt}+\mathbf y_\eta,
\qquad
\mathbf y_{\rm nt}=\begin{pmatrix}0\\0\\q_{\rm nt}\end{pmatrix},
\qquad
\mathbf y_\eta=-q_\eta\begin{pmatrix}0\\1\\1\end{pmatrix}.
\]
So after the first static gate is cleared, the remaining obstruction is already the equal-drift `K_\eta`–`\mu` ray.

Stage 021 then rewrote the rigid-mouth finite packet directly in the physical branch observables and proved
\[
q_{\rm nt}=
\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
q_\eta=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]
So the remaining obvious question is:

> what happens if we use these physical logarithmic variables themselves as the rigid-mouth coordinates?

This stage shows that doing so diagonalizes the packet completely and turns the microscopic correction problem into an exact two-axis compiler.

---

## 1. Exact rigid-mouth physical logarithmic chart

On the rigid-mouth branch define the physical logarithmic coordinates
\[
\boxed{
U:=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V:=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]
By Stage 021 these are exactly the surviving finite packet coordinates:
\[
\boxed{
q_{\rm nt}=U,
\qquad
q_\eta=V.
}
\]
So the rigid-mouth physical packet compiler is already diagonal:
\[
\boxed{
\mathbf q_{\rm rm}^{\rm phys}
:=
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_{\rm phys}
\begin{pmatrix}
U\\
V
\end{pmatrix},
\qquad
M_{\rm phys}=I_2.
}
\]

The third direct observable on the rigid-mouth branch, `R_{\rm target}`, is then recovered exactly from the selected-branch identity
\[
R_{\rm target}\,\mathcal T^2=
\Lambda_0(1-\epsilon_\eta),
\qquad
R_{\rm target,ref}\,\mathcal T_{\rm ref}^2=
\Lambda_0(1-\epsilon_{\eta,\rm ref}),
\]
which gives
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}\,e^{-U}.
}
\]
So the full rigid-mouth actual branch is already charted exactly by the diagonal logarithmic pair `(U,V)`.

---

## 2. Exact physical projectors and the two commuting finite legs

Because the rigid-mouth packet is diagonal in `(U,V)`, the canonical physical packet projectors are simply
\[
\boxed{
P_{\mathcal T}:=
\begin{pmatrix}1&0\\0&0\end{pmatrix},
\qquad
P_\eta:=
\begin{pmatrix}0&0\\0&1\end{pmatrix}.
}
\]
They are exact complementary projectors:
\[
P_{\mathcal T}^2=P_{\mathcal T},
\qquad
P_\eta^2=P_\eta,
\qquad
P_{\mathcal T}P_\eta=P_\eta P_{\mathcal T}=0,
\qquad
P_{\mathcal T}+P_\eta=I_2.
\]
So every rigid-mouth physical point decomposes uniquely as
\[
\boxed{
\begin{pmatrix}U\\V\end{pmatrix}
=
\begin{pmatrix}U\\0\end{pmatrix}
+
\begin{pmatrix}0\\V\end{pmatrix}.
}
\]

### 2.1 Pure transfer-shape leg

The image of `P_{\mathcal T}` is the exact finite transfer-shape leg
\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2 e^U,
\qquad
\epsilon_\eta=\epsilon_{\eta,\rm ref},
\qquad
\frac{R_{\rm target}}{R_{\rm target,ref}}=e^{-U}.
}
\]
So pure transfer-shape motion changes only `\mathcal T^2` directly and compensates `R_{\rm target}` inversely.

### 2.2 Pure dressing leg

The image of `P_\eta` is the exact finite dressing leg
\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2,
\qquad
\epsilon_\eta=\epsilon_{\eta,\rm ref}e^V,
\qquad
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}.
}
\]
So pure dressing motion is exactly the finite static-blind curve from Stage 020, now recognized as one coordinate axis of the physical chart.

### 2.3 Exact commutativity

Because the target ratio factorizes as
\[
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\underbrace{e^{-U}}_{\text{transfer leg}}
\underbrace{\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}}_{\text{dressing leg}},
\]
the two finite legs commute exactly.
So the rigid-mouth branch is an exact Cartesian product of

- transfer-shape motion, and
- dressing motion.

This is the physical normal-form version of the earlier packet projector calculus.

---

## 3. Exact physical-to-microscopic dependent-plane compiler

Stage 019 gives the dependent-plane packet carriers
\[
\mathbf y_{\rm nt}=\begin{pmatrix}0\\0\\q_{\rm nt}\end{pmatrix},
\qquad
\mathbf y_\eta=-q_\eta\begin{pmatrix}0\\1\\1\end{pmatrix}.
\]
Substituting the physical coordinates `q_{\rm nt}=U`, `q_\eta=V` gives the exact rigid-mouth dependent correction
\[
\boxed{
\mathbf y_{\rm rm}^{\rm dep}(U,V)
=
\begin{pmatrix}
\Delta_T\\
\Delta_{K_\eta}\\
\Delta_\mu
\end{pmatrix}
=
\begin{pmatrix}
0\\
-V\\
U-V
\end{pmatrix}.
}
\]
So the physical-to-microscopic compiler matrix is
\[
\boxed{
C_{\rm phys}^{\rm dep}
:=
\begin{pmatrix}
0&0\\
0&-1\\
1&-1
\end{pmatrix},
\qquad
\mathbf y_{\rm rm}^{\rm dep}=C_{\rm phys}^{\rm dep}
\begin{pmatrix}U\\V\end{pmatrix}.
}
\]
A left inverse is immediate:
\[
\boxed{
L_{\rm phys}^{\rm dep}
:=
\begin{pmatrix}
0&-1&1\\
0&-1&0
\end{pmatrix},
\qquad
L_{\rm phys}^{\rm dep}C_{\rm phys}^{\rm dep}=I_2.
}
\]
So the physical packet coordinates can be recovered directly from the dependent microscopic correction by
\[
\boxed{
U=\Delta_\mu-\Delta_{K_\eta},
\qquad
V=-\Delta_{K_\eta}.
}
\]

This is the cleanest rigid-mouth compiler obtained so far:

- `U` is the microscopic `\mu-K_\eta` difference,
- `V` is minus the `K_\eta` drift itself.

---

## 4. Exact microscopic images of the two physical axes

Push the physical projectors through `C_{\rm phys}^{\rm dep}`.

### 4.1 Pure transfer-shape image

Applying `P_{\mathcal T}` gives
\[
\boxed{
\mathbf y_{\mathcal T}^{\rm dep}
=
C_{\rm phys}^{\rm dep}
\begin{pmatrix}U\\0\end{pmatrix}
=
\begin{pmatrix}
0\\
0\\
U
\end{pmatrix}.
}
\]
So a pure transfer-shape defect is carried microscopically by a `\mu_W` shift only.

### 4.2 Pure dressing image

Applying `P_\eta` gives
\[
\boxed{
\mathbf y_\eta^{\rm dep}
=
C_{\rm phys}^{\rm dep}
\begin{pmatrix}0\\V\end{pmatrix}
=
-V\begin{pmatrix}0\\1\\1\end{pmatrix}.
}
\]
So a pure dressing defect is carried microscopically by the exact equal-drift `K_\eta^{(\mathrm{eff})}`–`\mu_W` ray.

This is the physical version of the Stage-019 dependent-plane theorem.
The difference is that the two microscopic pieces are now the direct images of the actual physical axes, not only of abstract quotient coordinates.

---

## 5. Exact correction compilers

Because the rigid-mouth packet is diagonal in `(U,V)`, the orbit-restoring corrections are immediate.

### 5.1 Static-only correction

To clear only the transfer-shape defect, subtract the pure transfer image:
\[
\boxed{
\Delta\mathbf y_{\rm static}
:=-\mathbf y_{\mathcal T}^{\rm dep}
=
\begin{pmatrix}0\\0\\-U\end{pmatrix}.
}
\]
After this correction,
\[
\mathbf y_{\rm rm}^{\rm dep}+
\Delta\mathbf y_{\rm static}
=
-V\begin{pmatrix}0\\1\\1\end{pmatrix},
\]
so the branch is moved exactly onto the pure dressing ray.

Thus the first static ceiling is cleared by changing only `\mu_W`.

### 5.2 Post-static dressing correction

Once the static ceiling has been cleared, the remaining orbit-restoring correction is just the opposite of the dressing image:
\[
\boxed{
\Delta\mathbf y_{\eta,\rm rest}
:=+V\begin{pmatrix}0\\1\\1\end{pmatrix}.
}
\]
So the extra step beyond the static gate is the equal positive shift in

- `K_\eta^{(\mathrm{eff})}`,
- `\mu_W`.

### 5.3 Full orbit-lock correction

Removing the whole rigid-mouth dependent correction gives
\[
\boxed{
\Delta\mathbf y_{\rm orbit}
:=-\mathbf y_{\rm rm}^{\rm dep}
=
\begin{pmatrix}
0\\
V\\
V-U
\end{pmatrix}
=
\Delta\mathbf y_{\rm static}+
\Delta\mathbf y_{\eta,\rm rest}.
}
\]
So the full orbit-restoring correction is literally the sum of

1. the static-only `\mu_W` correction,
2. the post-static equal `K_\eta`–`\mu_W` dressing correction.

This is the sharpest correction split reached so far.

---

## 6. Exact support-blindness of the physical normal form

Stage 021 already showed that the direct physical observables are support-blind:
\[
\partial_\zeta \mathcal T^2=0,
\qquad
\partial_\zeta \epsilon_\eta=0,
\qquad
\partial_{M_{\rm mix}}\mathcal T^2=0,
\qquad
\partial_{M_{\rm mix}}\epsilon_\eta=0.
\]
Therefore the physical logarithmic coordinates themselves satisfy
\[
\boxed{
\partial_\zeta U=
\partial_{M_{\rm mix}}U=
\partial_\zeta V=
\partial_{M_{\rm mix}}V=0.
}
\]
And because the microscopic compiler is linear in `(U,V)`, the dependent correction and all three correction compilers are support-blind as well:
\[
\boxed{
\partial_\zeta \mathbf y_{\rm rm}^{\rm dep}=
\partial_{M_{\rm mix}}\mathbf y_{\rm rm}^{\rm dep}=0,
}
\]
\[
\boxed{
\partial_\zeta \Delta\mathbf y_{\rm static}=
\partial_\zeta \Delta\mathbf y_{\rm orbit}=
0,
\qquad
\partial_{M_{\rm mix}} \Delta\mathbf y_{\rm static}=
\partial_{M_{\rm mix}} \Delta\mathbf y_{\rm orbit}=0.
}
\]
So coherent support enhancement cannot alter either

- the rigid-mouth physical packet, or
- the microscopic orbit-restoring correction required by that packet.

---

## 7. Cartesian orbit-lock theorem and first-order form

On the rigid-mouth actual branch,
\[
\boxed{
q_{\rm nt}=0,
\ q_\eta=0
\iff
U=0,
\ V=0
\iff
\mathcal T^2=\mathcal T_{\rm ref}^2,
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]
Because
\[
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}e^{-U},
\]
this is also equivalent to
\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2,
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}
\iff
R_{\rm target}=R_{\rm target,ref},
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]
So the rigid-mouth orbit-lock point is already a Cartesian codimension-two point in the physical logarithmic chart.

At first order the same statement becomes simply
\[
U=\delta\ln\mathcal T^2,
\qquad
V=\delta\ln\epsilon_\eta,
\]
so
\[
\boxed{
\mathbf y_{\rm rm}^{\rm dep,(1)}
=
\begin{pmatrix}
0\\
-\delta\ln\epsilon_\eta\\
\delta\ln\mathcal T^2-\delta\ln\epsilon_\eta
\end{pmatrix}.
}
\]
The static-blind line from Stages 018–020 is just the axis
\[
U=0,
\]
and its microscopic image is the pure equal-drift dressing ray.

So the earlier direct-space triangular packet was not the final normal form.
The later coherent transfer-shape compiler diagonalizes it completely.

---

## 8. Best current summary after Stage 022

The continuation from Stage 021 is now complete.

1. On the rigid-mouth actual branch, the physical logarithmic variables
   \[
   \left(\ln\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2},\ \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right)
   \]
   are already the exact finite packet coordinates.
2. In that chart, the packet projector calculus is diagonal.
3. The microscopic dependent correction is exactly
   \[
   (\Delta_T,\Delta_{K_\eta},\Delta_\mu)=(0,-V,U-V).
   \]
4. The static-only correction changes only `\mu_W`.
5. The post-static dressing correction is the equal `K_\eta`–`\mu_W` shift.
6. The whole rigid-mouth physical normal form and its microscopic correction compiler are support-blind.

So the same-charge barrier is now fully factorized on the rigid-mouth actual branch:

> first clear the pure transfer-shape axis `U`,
> then test whether the remaining pure dressing axis `V` collapses.

That is the sharpest Cartesian rigid-mouth statement reached so far.

same_charge_barrier_audit_stage023_selected_branch_loading_ratio_from_minimal_isotropic_quadrupole_precursor.md

# Same-Charge Barrier Audit — Stage 023: Selected-Branch Loading Ratio from the Minimal Isotropic Quadrupole Precursor

## Provenance

This standalone working note is extracted from the corresponding step in `g2_full_output.md`
so the stage-numbered same-charge barrier audit artifacts match the actual calculations that were done.



## Goal

Step 22 reduced the last support-selection ambiguity to the scalar demand parameter

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}.
```

But that still left one question open:

> **what does the selected passive/outgoing normalization side actually choose for
> `\Pi_{\rm tr}`?**

The later moving-throat notes answer that in two stages.

1. First, the selected-branch quadrupole-demand product **cancels** all separate
   outgoing-normalization amplitudes and depends only on the loading ratio
   ```math
   \rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}}.
   ```
2. Second, the natural minimal isotropic conservative quadrupole precursor fixes
   that loading ratio exactly through its contact/pole split.

So Step 23 is the point where the support selector `\varrho` is finally tied back
into the selected-branch normalization side.

The main result is:

```math
\boxed{
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43\,C_{\rm mix}.
}
```

Equivalently,

```math
\boxed{
\varrho = \frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43.
}
```

So the natural minimal isotropic passive/outgoing branch is **not** mixed-only and
**not** non-twin. It sits exactly on the symmetric-lowest-twin support slice.

---

## Inputs from the selected-branch normalization side

The later moving-throat notes give the exact product identities

```math
\Pi_{\rm tr}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm req},
\qquad
C_{\rm mix}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm mix}.
```

So immediately

```math
\boxed{
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}
=:\rho_\alpha.
}
```

This is the key simplification: once the outgoing quadrupole branch is normalized,
all the separate selected-mode amplitudes drop out of the support test.

In the spectral notation of the selected branch,

```math
N_Q^{(\rm target)} = \hat m_-^2\,\beta_0\,\frac{s_-}{\lambda_-},
```

so the same identities can also be written as

```math
\Pi_{\rm tr} = \hat m_-^2\frac{s_-}{\lambda_-}\alpha_{\rm req},
\qquad
C_{\rm mix} = \hat m_-^2\frac{s_-}{\lambda_-}\alpha_{\rm mix}.
```

Again the ratio is just `\rho_\alpha`.

---

## Step 23A — Exact contact-plus-pole inverse formulas

The natural minimal conservative quadrupole precursor is written as

```math
Y_Q^{\rm cons}(\omega)
=
 c_0 + \frac{c_1}{1-\omega^2/\Omega_Q^2},
```

with normalized static limit

```math
c_0 + c_1 = 1.
```

On the explicit support/source branch, the natural reading is:

- the mixed baseline carries the static contact fraction,
- the extra support lane carries the finite conservative pole.

So the same precursor can be written as

```math
Y_Q^{\rm cons}(\omega)
=
\frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
+
\frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Introducing

```math
\rho_\alpha := \frac{\alpha_{\rm req}}{\alpha_{\rm mix}},
```

this becomes

```math
Y_Q^{\rm cons}(\omega)
=
\frac{1}{\rho_\alpha}
+
\frac{\rho_\alpha-1}{\rho_\alpha}
\frac{1}{1-\omega^2/\Omega_Q^2}.
```

Therefore the exact contact/pole data are

```math
\boxed{c_0 = \frac{1}{\rho_\alpha}},
\qquad
\boxed{c_1 = \frac{\rho_\alpha-1}{\rho_\alpha}},
```

and the inverse formulas are

```math
\boxed{\rho_\alpha = \frac{1}{c_0} = \frac{1}{1-c_1}},
```

```math
\boxed{\zeta_{\rm req} := \rho_\alpha-1 = \frac{c_1}{c_0}.}
```

So the support/source loading ratio is directly encoded in the static contact /
pole split of the conservative quadrupole precursor.

---

## Step 23B — Matching to the minimal isotropic quadrupole module

The 2.5PN quadrupole audit already fixed the smallest viable isotropic
conservative precursor to

```math
c_0 = \frac34,
\qquad
c_1 = \frac14,
\qquad
\Omega_Q = \frac{3c_s}{2a}.
```

Inserting those values into the exact inverse formulas gives immediately

```math
\boxed{\rho_\alpha = \frac{1}{3/4} = \frac43,}
```

```math
\boxed{\zeta_{\rm req} = \frac{1/4}{3/4} = \frac13.}
```

Then the selected demand product is

```math
\boxed{
\Pi_{\rm tr}
=
\rho_\alpha C_{\rm mix}
=
\frac43 C_{\rm mix}.
}
```

So the selected branch is no longer carrying an arbitrary support demand. The
natural minimal isotropic passive/outgoing branch fixes it exactly.

---

## Step 23C — Exact support-selector form of the selected branch

Step 22 defined

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\qquad
C_{\rm mix} = \frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

Substituting

```math
\Pi_{\rm tr} = \frac43 C_{\rm mix}
```

gives

```math
\varrho
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\boxed{\frac{2(1-\epsilon_*)}{3}.}
```

And the required support enhancement is simply

```math
S_{\rm req}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\boxed{\frac43.}
```

So the selected branch is no longer scanning all support-demand sectors. It is
locked to one exact support ratio.

---

## Step 23D — Regime meaning

Stage 22 already split the support regimes by

```math
\Pi_{\rm tr} \le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
```

```math
C_{\rm mix} < \Pi_{\rm tr} \le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
```

```math
\Pi_{\rm tr} > 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
```

Because the selected branch gives

```math
\Pi_{\rm tr} = \frac43 C_{\rm mix},
```

it follows exactly that

```math
\boxed{
C_{\rm mix} < \Pi_{\rm tr} < 2C_{\rm mix}.
}
```

So:

- mixed-only is **not** enough,
- the symmetric lowest twin **is** enough,
- and non-twin asymmetry is **not** required.

This is already a real simplification of the anomaly bridge.

---

## Main result of the step

The selected-branch normalization side has now fixed the support ratio carried by
the natural minimal isotropic passive/outgoing quadrupole branch:

```math
\boxed{
\rho_\alpha = \frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.
}
```

Equivalently,

```math
\boxed{
\varrho = \frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac43.
}
```

So the last support ambiguity has collapsed from three sectors

- mixed-only,
- symmetric lowest twin,
- non-twin asymmetry,

to exactly **one** selected support slice:

```math
\text{symmetric lowest twin, with demand ratio } \Pi_{\rm tr}/C_{\rm mix}=4/3.
```

---

## What the next step should be

The next honest move is now very sharp:

> restrict the Step-21 primitive quartic ranking problem to this selected
> twin-support branch, and ask how much of the remaining `q_W` versus
> `q_\Lambda` ambiguity survives there.

That is the smallest next derivation that still genuinely pushes the anomaly
closure forward.
same_charge_barrier_audit_stage024_exact_primitive_ranking_on_selected_twin_support_branch.md

# Same-Charge Barrier Audit — Stage 024: Exact Primitive Ranking on the Selected Twin-Support Branch

## Provenance

This standalone working note is extracted from the corresponding step in `g2_full_output.md`
so the stage-numbered same-charge barrier audit artifacts match the actual calculations that were done.



## Goal

Step 23 made the main support-side simplification:

```math
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \frac43.
```

So the natural minimal isotropic passive/outgoing quadrupole branch is **not**
allowed to roam over all three support sectors anymore. It lives on exactly one
selected support slice:

```math
\text{symmetric lowest twin, with } \Pi_{\rm tr}/C_{\rm mix}=4/3.
```

That means the old Step-21/22 phase diagram is now overkill. The real next
question is narrower:

> **once we restrict to that selected twin-support curve, how much of the
> primitive quartic ranking ambiguity survives?**

This step answers that exactly.

The main result is that the whole selected branch is one exact curve

```math
\boxed{\epsilon_* = 1 - \frac{3\varrho}{2}},
\qquad
\boxed{\sigma = \frac{4}{3\varrho}-2},
\qquad
0<\varrho<\frac23,
```

and along that curve only **two** ranking thresholds remain:

```math
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
}
```

```math
\boxed{
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
```

So the full selected-branch primitive ranking is:

```math
\boxed{
\begin{aligned}
&0<\varrho<\varrho_{W\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_\Lambda > q_W > |q_U|,\\[1mm]
&\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_W > q_\Lambda > |q_U|,\\[1mm]
&\varrho_{U\Lambda}<\varrho<\frac23
&&\Longrightarrow&&
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
\end{aligned}
}
```

That is the cleanest anomaly ranking statement reached so far.

---

## Step 24A — The selected branch is an exact one-parameter twin-support curve

Step 22 defined

```math
\varrho := \frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\qquad
\sigma = \frac{2\epsilon_*}{1-\epsilon_*}.
```

Step 23 then fixed the selected support ratio to

```math
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \frac43,
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
```

So

```math
\varrho
=
\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\frac{2(1-\epsilon_*)}{3}.
```

Hence

```math
\boxed{\epsilon_* = 1 - \frac{3\varrho}{2}.}
```

Since `0<\epsilon_*<1`, this gives the exact selected-branch range

```math
\boxed{0<\varrho<\frac23.}
```

Now convert to `\sigma`:

```math
\sigma
=
\frac{2\epsilon_*}{1-\epsilon_*}
=
\frac{2\bigl(1-3\varrho/2\bigr)}{3\varrho/2}
=
\boxed{\frac{4}{3\varrho}-2.}
```

So the selected branch is not a 2D region in `(epsilon_*,\varrho)` at all. It
is a single exact curve.

---

## Step 24B — The selected curve sits strictly inside the twin window

Step 22 gave the exact support windows in `\sigma`:

```math
0<\sigma\le \frac1\varrho-2
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
```

```math
\frac1\varrho-2 < \sigma \le \frac2\varrho-2
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
```

```math
\sigma > \frac2\varrho-2
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
```

On the selected branch,

```math
\sigma_{\rm sel} = \frac{4}{3\varrho}-2.
```

Then

```math
\sigma_{\rm sel} - \left(\frac1\varrho-2\right)
=
\frac{1}{3\varrho} > 0,
```

and

```math
\left(\frac2\varrho-2\right) - \sigma_{\rm sel}
=
\frac{2}{3\varrho} > 0.
```

So for every allowed point on the selected branch,

```math
\boxed{
\frac1\varrho-2 < \sigma_{\rm sel} < \frac2\varrho-2.
}
```

That is the exact proof that the selected branch lies **strictly inside** the
symmetric-lowest-twin regime.

So mixed-only and non-twin branches are gone from the live anomaly closure.

---

## Step 24C — Surviving threshold 1: `q_W` versus `q_\Lambda`

Step 21 gave the exact crossover

```math
q_W = q_\Lambda
\iff
\epsilon_* = \frac{1}{2+\beta^2}.
```

Insert the selected-branch law

```math
\epsilon_* = 1 - \frac{3\varrho}{2}.
```

Then the threshold on the selected branch is

```math
1 - \frac{3\varrho}{2} = \frac{1}{2+\beta^2},
```

hence

```math
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)}.
}
```

Therefore:

- if
  ```math
  0<\varrho<\varrho_{W\Lambda},
  ```
  then `q_\Lambda > q_W`;
- if
  ```math
  \varrho>\varrho_{W\Lambda},
  ```
  then `q_W > q_\Lambda`.

So the outgoing-scale lane overtakes the wall-blocking lane only in the **low-`\varrho` / high-blocking** corner of the selected curve.

---

## Step 24D — Surviving threshold 2: `|q_U|` versus `q_\Lambda`

Step 21 also gave the exact crossover

```math
|q_U| = q_\Lambda
\iff
\epsilon_* = \frac{\beta}{1+\beta+\beta^2}.
```

Again insert

```math
\epsilon_* = 1 - \frac{3\varrho}{2},
```

and solve for `\varrho`:

```math
\boxed{
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
```

So:

- if
  ```math
  \varrho<\varrho_{U\Lambda},
  ```
  then `q_\Lambda > |q_U|`;
- if
  ```math
  \varrho>\varrho_{U\Lambda},
  ```
  then `|q_U| > q_\Lambda`.

This is the selected-branch version of Step 21’s “very weak blocking” corner.

---

## Step 24E — Ordering of the two thresholds

The two thresholds are not independent. Their difference is

```math
\varrho_{U\Lambda} - \varrho_{W\Lambda}
=
\frac{2(1+\beta^2)(1-\beta)}{3(1+\beta+\beta^2)(2+\beta^2)} > 0
```

because `0<\beta<2/11<1`.

And

```math
\frac23 - \varrho_{U\Lambda}
=
\frac{2\beta}{3(1+\beta+\beta^2)} > 0.
```

So the exact threshold ordering on the selected branch is

```math
\boxed{0 < \varrho_{W\Lambda} < \varrho_{U\Lambda} < \frac23.}
```

That means the selected twin-support curve always splits into **three** ranking
regions and never fewer.

---

## Step 24F — Full primitive ranking on the selected twin-support branch

Step 21 already proved the branch-independent ordering facts

```math
q_\chi > q_Z,
\qquad
q_Z > q_W,
\qquad
q_W > |q_U|.
```

So only `q_\Lambda` moves relative to `q_W` and `|q_U|`.
Using the two selected-branch thresholds above, the complete ranking is now exact.

### Region I — low `\varrho`, strong blocking

If

```math
0<\varrho<\varrho_{W\Lambda},
```

then

```math
\boxed{q_\chi > q_Z > q_\Lambda > q_W > |q_U|.}
```

### Region II — intermediate `\varrho`

If

```math
\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda},
```

then

```math
\boxed{q_\chi > q_Z > q_W > q_\Lambda > |q_U|.}
```

### Region III — large `\varrho`, very weak blocking

If

```math
\varrho_{U\Lambda}<\varrho<\frac23,
```

then

```math
\boxed{q_\chi > q_Z > q_W > |q_U| > q_\Lambda.}
```

So the selected anomaly branch now has a completely explicit primitive ranking
phase diagram.

---

## Step 24G — Numerical size of the surviving thresholds

Using the constructive coherent bound

```math
0<\beta<\frac{2}{11},
```

one finds

```math
\boxed{
\frac13 < \varrho_{W\Lambda} < \frac{125}{369} \approx 0.338753,
}
```

and

```math
\boxed{
\frac{250}{441} \approx 0.566893 < \varrho_{U\Lambda} < \frac23.
}
```

So the selected twin-support curve has a very clean geometry.

- Only the **low-`\varrho`** end allows `q_\Lambda` to beat `q_W`.
- Across the middle of the selected curve, `q_W` beats `q_\Lambda` but
  `q_\Lambda` still beats `|q_U|`.
- Only near the **large-`\varrho` / very weak-blocking** end does `|q_U|`
  overtake `q_\Lambda`.

That is already much sharper than the full constructive-branch picture from Step 21.

---

## Main result of the step

The natural minimal isotropic passive/outgoing branch has collapsed the old
support-selection problem to one exact twin-support curve:

```math
\boxed{\epsilon_* = 1 - \frac{3\varrho}{2}},
\qquad
\boxed{\sigma = \frac{4}{3\varrho}-2},
\qquad
0<\varrho<\frac23.
```

On that curve, the primitive quartic hierarchy is controlled by only two exact
thresholds:

```math
\boxed{
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\qquad
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
}
```

So the complete selected-branch ranking is

```math
\boxed{
\begin{aligned}
&0<\varrho<\varrho_{W\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_\Lambda > q_W > |q_U|,\\[1mm]
&\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda}
&&\Longrightarrow&&
q_\chi > q_Z > q_W > q_\Lambda > |q_U|,\\[1mm]
&\varrho_{U\Lambda}<\varrho<\frac23
&&\Longrightarrow&&
q_\chi > q_Z > q_W > |q_U| > q_\Lambda.
\end{aligned}
}
```

This is the strongest quartic anomaly ranking statement reached so far from the
selected moving-throat branch side.

---

## What the next step should be

The next honest move is now very narrow:

> derive the actual physical position of the moving-throat branch on this
> selected twin-support curve — equivalently, pin `\epsilon_*` or `\varrho`
> rather than leaving it parametric.

Once that single coordinate is known, the quartic carrier hierarchy stops being a
phase diagram and becomes one definite anomaly prediction.
# Same-Charge Barrier Audit — Stage 025
## Actual Twin-Support Placement and Coherent Orbit-Lock Compiler

## Purpose

Stage 024 closed the reduced support-selection algebra by restricting the selected
same-charge branch to the exact one-parameter twin-support curve

\[
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma = \frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23.
\]

The remaining task is no longer another ranking theorem. It is to place the
**actual coherent moving-throat branch** on that curve and then evaluate the
separate orbit-lock packet.

This note compiles those two jobs into one exact front-end stage.

---

## 1. Actual selected-twin placement

On the coherent local D/N branch, the stationary placement packet is

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta).
\]

The exact coherent reduction uses

\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
\]

The selected same-charge branch already carries the Stage-023 support demand

\[
\Pi_{\rm tr}=\frac43 C_{\rm mix},
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

Therefore the actual selected-branch coordinate on the Stage-024 twin-support
curve is

\[
\varrho_{\rm phys}
:=
\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=
\frac23(1-\epsilon).
\]

Hence

\[
\epsilon_{*,\rm phys}=\epsilon,
\qquad
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}=\frac{4}{3\varrho_{\rm phys}}-2=
\frac{2\epsilon}{1-\epsilon}.
\]

So the parametric Stage-024 curve is converted into a realized support point as
soon as the completed PDE returns the coherent packet
\((\delta_U,\epsilon_W)\).

---

## 2. Exact threshold rewrite in the realized variable \(\epsilon\)

Stage 024 found the two primitive thresholds

\[
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\qquad
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
\]

Using
\(\epsilon_*=1-\tfrac32\varrho\),
these become

\[
\epsilon_{W\Lambda}=\frac{1}{2+\beta^2},
\qquad
\epsilon_{U\Lambda}=\frac{\beta}{1+\beta+\beta^2}.
\]

So the realized primitive ranking is decided directly by the actual coherent
support variable \(\epsilon\):

\[
\epsilon>\epsilon_{W\Lambda}
\iff
0<\varrho<\varrho_{W\Lambda},
\]

\[
\epsilon_{U\Lambda}<\epsilon<\epsilon_{W\Lambda}
\iff
\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda},
\]

\[
0<\epsilon<\epsilon_{U\Lambda}
\iff
\varrho_{U\Lambda}<\varrho<\frac23.
\]

---

## 3. Support-lane classifier on the realized branch

The coherent support classifier is still

\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\text{mixed-only enough},
\]

\[
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\text{lowest symmetric twin enough},
\]

\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\text{non-twin asymmetry required}.
\]

On the selected same-charge branch,

\[
\Pi_{\rm tr}=\frac43 C_{\rm mix},
\]

so the realized support slice is **automatically** in the lowest symmetric twin
window:

\[
C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}.
\]

Important notation point:

- the Stage-023 demand-side ratio \(\zeta_{\rm req}=1/3\) belongs to the reduced
  contact/pole loading problem,
- the physical coherent local D/N support ratio \(\zeta\) is a separate support
  variable,
- and on the lowest symmetric twin branch the physical coherent support value is
  \(\zeta=1\).

These two uses of \(\zeta\) must not be conflated.

---

## 4. Exact orbit-side observables

The coherent orbit packet is carried by

\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]

\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]

together with \(\epsilon_\eta\).

The support ratio \(\zeta\) does **not** enter
\(\epsilon\),
\(R_{\rm tr}\),
or
\(R_{\rm target}\).
So the orbit packet is exactly support-blind.

---

## 5. Infinitesimal compiler

Let the actual coherent placement drifts be

\[
(d\ln\chi_0,
 d\ln\delta_U,
 d\ln Z_W,
 d\ln\epsilon_W,
 d\ln\epsilon_\eta,
 d\ln\Lambda).
\]

Then

\[
d\ln\epsilon
=
d\ln\epsilon_W
-
\frac{2\delta_U}{(1+\delta_U)(11+9\delta_U)}d\ln\delta_U.
\]

Also

\[
d\ln R_{\rm tr}
=
-
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\Big[(1+\delta_U)d\ln\chi_0+(1+\chi_0)d\ln\delta_U\Big].
\]

And

\[
d\ln R_{\rm target}
=
d\ln\Lambda-d\ln Z_W
-
\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-
\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-
\frac{2\epsilon}{1-\epsilon}d\ln\epsilon.
\]

The direct observable weak-axisymmetric defect packet is therefore

\[
\Theta_1=d\ln R_{\rm tr},
\]

\[
\Xi_1
=
-d\ln R_{\rm target}-\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
=
-d\ln\Lambda+d\ln Z_W+
\frac{2\chi_0}{1+\chi_0}d\ln\chi_0+
\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
\]

\[
\mathcal R_1=d\ln R_{\rm target}.
\]

So \(\Xi_1\) is explicitly blind to both the support lane \(\zeta\) and the
selected-branch demand variable \(\epsilon_\eta\).

---

## 6. Exact coherent orbit-lock theorem gate

The actual coherent local D/N branch satisfies orbit lock iff

\[
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
\]

equivalently iff

\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln\epsilon_\eta=0.
\]

The outgoing finish line remains separate:

\[
N_Q=1,
\qquad
\text{or equivalently }\chi_Q=1
\text{ on the natural source-map branch.}
\]

So the actual remaining 5PN / moving-throat endgame is exactly:

1. extract the stationary coherent placement state,
2. compute the weak-axisymmetric tangent,
3. place the branch on the selected twin-support curve via \(\epsilon\),
4. evaluate
   \(d\ln R_{\rm tr}\),
   \(d\ln R_{\rm target}\),
   \(d\ln\epsilon_\eta\),
   and \(N_Q-1\).

---

## 7. Practical Stage-025 output packet

The smallest exact compiler packet to return from the actual moving-throat branch is

\[
\Bigl(
\epsilon,
\varrho_{\rm phys},
\sigma_{\rm phys},
\text{ranking region},
R_{\rm tr},
R_{\rm target},
\epsilon_\eta,
 d\ln R_{\rm tr},
 d\ln R_{\rm target},
 d\ln\epsilon_\eta,
 N_Q-1
\Bigr).
\]

That is the sharp front edge after Stage 024.
