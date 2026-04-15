# Moving-Throat PDE — Stage 203: Dynamic Mixed-Port Kernel, Phase-Lag No-Go, and the Resonant-Survival Gate

## Status

**Exact within the carried isotropic one-port wall/BdG/Maxwell/mixed closure, away from the internal poles**
\[
\varpi^2-\omega^2\neq 0,
\qquad
\Delta_\Pi(\omega)\neq 0,
\qquad
D_\Pi(\omega)\neq 0,
\]
and exact to first order in the outgoing-port dressing when the phase-lag result is invoked.

This is the first post-Stage-202 dynamic insertion of actual mixed-port bundle data into the now-closed local mixed-ray ledger.
It does **not** solve the full driven two-throat moving PDE.
It computes the first honest frequency-dependent same-charge kernel already implied by the one-port mixed bundle and identifies the only linear dynamic corridor left alive after the static audit.

---

## Purpose

Stage 202 closed the first honest static same-charge mixed-bundle audit.
It showed that the minimal one-port wall/BdG/Maxwell/mixed bundle does **not** create a new long-range attractive law.
It only generates short-range product families built from the primitive source profiles.

That leaves one clean escape hatch:

> perhaps a genuinely time-dependent / non-adiabatic mixed-port drive can do something qualitatively stronger than the static bundle.

The job of Stage 203 is narrower than solving the full driven PDE.
It is to:

1. promote the Stage-202 one-port bundle to a frequency-dependent mixed-port kernel,
2. derive its exact complex susceptibility,
3. test whether linear monochromatic driving creates a genuinely new spatial attractive family,
4. separate the in-phase barrier reshaping from the out-of-phase pumping / leakage channel,
5. and identify the only surviving linear dynamic corridor.

The main outputs are:

1. the exact dynamic reduced `3 x 3` bundle and determinant identity,
2. the exact dynamic susceptibility law and inverse-entry formulas,
3. the exact collinear-source and primitive-source product-family theorems,
4. the exact outgoing-port derivative identity,
5. the exact linear phase-lag theorem,
6. and the resulting survival gate: **the only linear dynamic corridor left is resonant dispersive enhancement of the already-known short-range families, and that corridor must beat its own absorptive load to count as real barrier engineering.**

So Stage 203 keeps the dynamic same-charge idea alive, but only in a much narrower form than the raw mixed-port intuition suggested.

---

## 1. Frozen input carried forward

### 1.1 Static one-port bundle from Stage 202

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

Stage 202 proved that the static mixed correction is always attractive or neutral on the admissible branch and that, for the primitive reduced source families,
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

So the linear dynamic audit has to be judged against that already-frozen static baseline.

### 1.2 Why the dynamic mixed sector is the right next lane

The parent `4+1` / plasma ontology keeps the mixed variables
\[
A_w,
\qquad
F_{\mu w},
\qquad
J^w,
\]
alive once one relaxes the strict far-field brane reduction, and it treats brane-facing non-ideality as conservative transport into the hidden `w` sector and higher localization modes rather than as purely local dissipation.
So if a field-driven corridor survives at all, it has to pass through the dynamic mixed sector rather than through the static zero-mode Maxwell law.

### 1.3 Outgoing-port motivation

The moving-throat program already isolated the mixed-sector port as the first place where a genuine passive/outgoing channel can enter the reduced wall language.
On the outgoing `l=2` branch the first odd fingerprint begins at `i\omega^5`, not at an even static correction.
So the same reduced sector that can carry dynamic pumping is also the sector that carries the radiative quadrupole bridge.

That is exactly why the dynamic barrier audit has to distinguish **real in-phase reshaping** from **imaginary out-of-phase pumping**.

---

## 2. Minimal dynamic one-port mixed bundle

Keep the same reduced coordinates as in Stage 202:

- wall/worldtube amplitude `q`,
- brane-like internal gauge coordinate `U`,
- mixed `A_w/F_{\mu w}/J^w` coordinate `W`.

Now retain the frequency dependence explicitly.

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
where `\Pi(\omega)` is the effective outgoing/self-energy carried by the mixed coordinate.
Then define
\[
\Delta_\Pi(\omega):=A(\omega)W(\omega)-R^2,
\]
\[
Q_\Pi(\omega):=G_U^2W(\omega)+2G_UG_WR+G_W^2A(\omega),
\]
\[
\boxed{
D_\Pi(\omega):=K_B(\omega)-\frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)}.
}
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

- the internal `(U,W)` block denominator `\Delta_\Pi`,
- the dressed wall denominator `D_\Pi`.

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
so the dynamic bundle reduces exactly to the Stage-202 static one-port bundle.

So Stage 203 is a genuine continuation, not a different reduced model.

---

## 3. Exact dynamic susceptibility law

Define the complex reduced quadratic response functional
\[
\mathfrak V_{\rm mix}(x,\omega)
:=
-\frac12 J(x,\omega)^T\mathcal K_{\rm dyn}(\omega)^{-1}J(x,\omega).
\]

This is the direct dynamic analogue of the Stage-202 static on-shell shift.

- `\Re\,\mathfrak V_{\rm mix}` is the in-phase conservative barrier reshaping.
- `\Im\,\mathfrak V_{\rm mix}` is the quadrature / pumping / leakage channel.

### 3.1 Exact inverse entries

The inverse entries are the direct dynamic lifts of the Stage-202 formulas:
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

So every dynamic same-charge correction is still controlled by the same denominator pair `\Delta_\Pi D_\Pi`.

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

Return to the same primitive reduced source profiles used in Stage 202:
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

This is the first strong dynamic theorem of the stage:

> the linear monochromatic mixed bundle does **not** create a new spatial kernel family.
> It keeps exactly the same three short-range spatial families as the static one-port bundle, but promotes their coefficients to complex frequency-dependent susceptibilities.

So linear time dependence does **not** yet buy a qualitatively new radial law.
It buys only:

1. frequency-dependent renormalization of the same short-range families,
2. complex phase structure,
3. and the possibility of resonant amplification through the denominator pair `\Delta_\Pi D_\Pi`.

---

## 5. Exact outgoing-port derivative identity

Let
\[
e_W=(0,0,1)^T,
\]
so that the mixed coordinate is the `W` lane.
Because `\Pi(\omega)` enters only in the `WW` slot,
\[
\partial_\Pi \mathcal K_{\rm dyn}(\omega)
=
- e_W e_W^T.
\]

Matrix calculus therefore gives the exact identity
\[
\boxed{
\partial_\Pi \mathcal K_{\rm dyn}(\omega)^{-1}
=
\mathcal K_{\rm dyn}(\omega)^{-1}
\,e_W e_W^T\,
\mathcal K_{\rm dyn}(\omega)^{-1}.
}
\]
Hence the dynamic response functional obeys
\[
\boxed{
\partial_\Pi \mathfrak V_{\rm mix}(x,\omega)
=
-\frac12
\bigl[e_W^T\mathcal K_{\rm dyn}(\omega)^{-1}J(x,\omega)\bigr]^2.
}
\]

So for a small outgoing port around the conservative branch `\Pi=0`,
\[
\boxed{
\delta\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\Pi(\omega)\,\mathcal T_J(\omega)^2
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

This is the exact dynamic lift of the Stage-004 outgoing-transfer factor, now written directly for the driven same-charge source vector.

---

## 6. Phase-lag no-go at linear order

Take a passive/outgoing mixed port at linear order:
\[
\Pi(\omega)=i\Gamma(\omega)+O\bigl(\Gamma(\omega)^2\bigr),
\qquad
\Gamma(\omega)\ge 0.
\]
On the conservative off-pole branch, `\mathcal K_{\rm cons}(\omega)` and the reduced source data are real, so `\mathcal T_J(\omega)` is real.
Therefore the first outgoing correction becomes
\[
\boxed{
\delta\mathfrak V_{\rm mix}^{(1)}(x,\omega)
=
-\frac{i}{2}\Gamma(\omega)\,\mathcal T_J(\omega)^2.
}
\]
Its real and imaginary parts are then
\[
\boxed{
\Re\,\delta\mathfrak V_{\rm mix}^{(1)}=0,
\qquad
\Im\,\delta\mathfrak V_{\rm mix}^{(1)}=-\frac12\Gamma(\omega)\,\mathcal T_J(\omega)^2.
}
\]
So the first outgoing correction is **purely phase-lag / pumping** at linear order.
It does **not** lower the barrier conservatively.

Define the absorbed-power diagnostic by
\[
\boxed{
\overline P_{\rm abs}(x,\omega)
:=
-\omega\,\Im\mathfrak V_{\rm mix}(x,\omega).
}
\]
Then at first outgoing order,
\[
\boxed{
\overline P_{\rm abs}^{(1)}(x,\omega)
=
\frac{\omega\Gamma(\omega)}{2}\,\mathcal T_J(\omega)^2
\ge 0.
}
\]

So the first passive/outgoing mixed-port correction is not barrier softening.
It is dissipative loading in exactly the sense the 4D/plasma ontology suggested.

This is the Stage-203 phase-lag no-go theorem.

---

## 7. What still survives: the resonant dispersive corridor

The previous section kills the naive dynamic story:

> the first passive/outgoing dynamic correction is not a new attractive contribution.
> It is pure quadrature at linear order.

So the only linear dynamic corridor left is the **even dispersive** dependence already present in the conservative coefficients through
\[
K_B(\omega),
\qquad
A(\omega),
\qquad
W(\omega),
\]
and therefore through the conservative denominator pair
\[
\Delta_0(\omega)=\bigl(\Omega_U^2-\omega^2\bigr)\bigl(\Omega_W^2-\omega^2\bigr)-R^2,
\]
\[
D_0(\omega)=K_B(\omega)-\frac{Q_0(\omega)}{\Delta_0(\omega)},
\qquad
Q_0(\omega):=Q_\Pi(\omega)\big|_{\Pi=0}.
\]

Away from poles this is only an even analytic deformation of the static coefficients.
Only near zeros of `\Delta_0` or `D_0` can the real in-phase response be parametrically amplified.
But those are exactly the regions where the same transfer factors also amplify `\mathcal T_J` and therefore the absorptive channel once the passive/outgoing port is restored.

So the linear dynamic survival gate is no longer a kernel-class question.
It is the following resonance test:

> can the real dispersive enhancement of the already-known short-range attractive families become large enough **before** the absorptive / leakage channel simply turns the mechanism into disguised heating or branch loss?

That is the correct Stage-203 meaning of “resonant survival.”
The full residue/linewidth tradeoff is deferred to the next stage.

---

## 8. Updated reduced barrier audit after Stage 203

At linear monochromatic order the same-charge audit potential is best written as the **real** part of the driven response correction:
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
The pumping / leakage diagnostic is
\[
\boxed{
\overline P_{\rm abs}(x,\omega)
=
-\omega\,\Im\mathfrak V_{\rm mix}(x,\omega).
}
\]

So any linear dynamic claim now has to clear **both** tests:

1. `\Re\,\mathfrak V_{\rm mix}` must materially reduce the barrier,
2. `\overline P_{\rm abs}` must stay small enough that the effect is not just covert heating / leakage.

That is the first honest dynamic kill test.

---

## 9. Best current summary after Stage 203

Stage 203 closes the first honest **linear dynamic** same-charge audit of the one-port moving-throat mixed bundle.

- The dynamic reduced bundle is an exact `3 x 3` frequency-domain system with determinant `\Delta_\Pi D_\Pi`.
- At zero frequency with `\Pi(0)=0`, it reduces exactly to the Stage-202 static bundle.
- For the first primitive same-charge source families, the linear dynamic bundle still produces only
  \[
  x^{-6},\qquad e^{-2\kappa x}/x^4,\qquad e^{-4\kappa x}/x^2.
  \]
- Time dependence only makes their coefficients complex and frequency dependent; it does **not** introduce a new spatial kernel family.
- The first passive/outgoing correction is exactly pure phase lag:
  \[
  \Re\,\delta\mathfrak V_{\rm mix}^{(1)}=0,
  \qquad
  \overline P_{\rm abs}^{(1)}=\frac{\omega\Gamma}{2}\mathcal T_J^2\ge0.
  \]
- Therefore the only surviving linear dynamic corridor is resonant dispersive enhancement of the already-known short-range families, and that corridor must defeat its own absorptive load to count as real barrier engineering.

So the dynamic mixed-port idea survives, but only as a narrow resonance / quality-factor problem, not as a linear monochromatic barrier-bypass law.

---

## 10. Script-backed status

The accompanying SymPy audit verifies:

1. the determinant identity
   \[
   \det \mathcal K_{\rm dyn}=\Delta_\Pi D_\Pi,
   \]
2. the exact static reduction back to the Stage-202 one-port bundle,
3. the exact inverse-entry formulas,
4. the exact dynamic susceptibility formula
   \[
   \mathfrak V_{\rm mix}=-\tfrac12 J^T\mathcal K_{\rm dyn}^{-1}J,
   \]
5. the collinear-source factorization theorem,
6. the primitive-source product-family theorem preserving the `x^{-6}`, `e^{-2\kappa x}/x^4`, and `e^{-4\kappa x}/x^2` families,
7. the exact outgoing-port derivative identity
   \[
   \partial_\Pi \mathfrak V_{\rm mix}=-\tfrac12\mathcal T_J^2,
   \]
8. the linear outgoing correction
   \[
   \delta\mathfrak V_{\rm mix}^{(1)}=-\tfrac{i}{2}\Gamma\mathcal T_J^2,
   \]
9. and the phase-lag consequence
   \[
   \Re\,\delta\mathfrak V_{\rm mix}^{(1)}=0,
   \qquad
   \overline P_{\rm abs}^{(1)}>0
   \]
   on a constructive off-pole slice.

Supporting file:
- `moving_throat_pde_stage203_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py`

---

## 11. Immediate next step

The next theorem gate is now very sharp.

1. Keep the exact Stage-203 dynamic one-port kernel.
2. Move to the local simple-pole / linewidth normal form near the first admissible internal resonances.
3. Quantify the residue-to-linewidth tradeoff between
   \[
   -\Re\,\mathfrak V_{\rm mix}
   \]
   and
   \[
   \overline P_{\rm abs}.
   \]
4. Ask whether any resonance window gives a materially larger conservative reshaping before absorptive loading takes over.

That is the cleanest continuation point after the phase-lag no-go theorem and the first resonant-survival gate.
