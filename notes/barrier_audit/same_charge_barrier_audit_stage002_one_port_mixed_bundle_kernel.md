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
