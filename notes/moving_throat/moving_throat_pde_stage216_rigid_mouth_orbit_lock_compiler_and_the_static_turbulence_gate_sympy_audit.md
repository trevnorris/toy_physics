# Moving-Throat PDE — Stage 216: Rigid-Mouth Orbit-Lock Compiler and the Static Turbulence Gate

## Status

**Exact within the carried Stage-207 / Stage-215 same-charge reduction plus the Stage-171 coherent branch-observable compiler.**

This stage does **not** claim that hydrodynamic turbulence or cavitation has already been derived from the full PDE.
It shows something narrower and more useful:

> if one adopts the physical reading that the brane mouth is rigid while the internal branch can still repackage loading, then the unresolved same-charge kill test is exactly a bound on the internal orbit/transfer packet, and the Stage-207 static ceiling can be read as a first reduced “choke / turbulence” gate.

The theorem content is exact inside the reduced branch language.
The words “turbulence”, “choked flow”, and “collapse” are **interpretive labels** for leaving that surviving static branch window.

---

## Purpose

Stage 215 already showed that the explicit `5`PN support/source branch is numerically safe and that the first unresolved same-charge bottleneck is still the **PDE-selected orbit-lock / coherent placement point**. The same notes are explicit that the numerically missing object is the actual point satisfying the coherent placement conditions
\[
d\ln R_{\rm tr}=0,\qquad
d\ln R_{\rm target}=0,\qquad
d\ln \epsilon_\eta=0.
\]

So the natural next question is no longer whether support/source is large enough.
It is:

> how should the unresolved orbit packet be interpreted if the defect mouth is geometrically rigid and the remaining load can only reorganize internally?

This stage compiles that idea into exact reduced formulas.

The main result is sharp:

> on a track-locked branch, the surviving static same-charge scalar is literally the drift of the corrected internal transfer observable \(\mathfrak N_*\), and with the stronger operator-rigidity hypothesis \(D_{01}=0\) the same static gate becomes a pure internal outgoing-transfer bound.

---

## 1. Interpretive firewall

Two distinctions matter.

### 1.1 Rigid mouth is **not** identical to `D_{01}=0`

The geometric statement “the brane entrance is rigid” is a statement about the mouth geometry.
The algebraic condition
\[
D_{01}=0
\]
is a stronger statement: it says the **effective static grouped-lane operator** does not drift at first weak-axisymmetric order.

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

This gives the first sharp reduced translation of the rigid-mouth picture:

- \(R_{\rm tr}\) measures tracking / mouth-side placement,
- \(\mathfrak N_*\) measures corrected internal transfer / nontracking load,
- \(\epsilon_\eta\) measures selected-branch dressing.

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
the logarithmic drift of the corrected internal transfer observable \(\mathfrak N_*\).

That is the precise reduced sense in which “the entrance is rigid but the inside can still load” is already visible in the math.

---

## 4. Stronger operator-rigidity specialization

Stages 207 and 208 already showed that the weak-axisymmetric static load scalar is
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

Stage 207 already transported the primitive-family same-charge window onto the actual branch packet:
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

And with the stronger operator-rigidity closure \(D_{01}=0\),
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

## 6. The Stage-207 numerical budgets in this language

At the Stage-206 compatibility point
\[
\bar P_0 \approx 0.002069792318062885,
\]
Stage 207 gave the strict `10%` robust budget
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

Under the additional operator-rigidity closure \(D_{01}=0\), the same become

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

So the Stage-207 `36.8%` figure really is the first exact reduced scalar at which the surviving static same-charge branch stops being robust.

---

## 7. Why this is compatible with the support/source verdict

Stage 215 already showed that the explicit `5`PN support/source branch overshoots the canonical isotropic demand by a large margin and is not the active bottleneck. What is still missing is the actual PDE-selected orbit-lock / coherent placement point.

So if one adopts the rigid-mouth reading, the situation is:

- support/source can feed the branch,
- the dynamic wall-window is not the first failure,
- the first unresolved gate is whether the **internal transfer / placement observable** \(\mathfrak N_*\) remains inside the Stage-207 static ceiling.

That is exactly the mathematical form of the “internal choke versus rigid entrance” picture.

---

## 8. Best current verdict after Stage 216

The rigid-mouth physical interpretation is **mostly compatible** with the reduced math, but with two important corrections:

1. **Good match:** it is reasonable to read the surviving same-charge problem as a rigid-mouth / internal-load competition, because on the exact observable compiler
   \[
   \delta\ln R_{\rm tr}=0
   \quad\Longrightarrow\quad
   \Xi_1=\delta\ln\mathfrak N_*,
   \]
   so the unresolved scalar is literally an internal transfer / placement drift.

2. **Important correction:** the statement
   \[
   D_{01}=0
   \]
   is **not** the same thing as “the entrance radius cannot change.”
   It is a stronger effective-static-operator rigidity condition.
   It is useful for the reduced thought experiment, but it should not be confused with the geometric mouth lock itself.

3. **Interpretive caution:** the Stage-207 `36.8%` ceiling is an exact reduced static branch bound. Calling it a “turbulence threshold” is a plausible physical interpretation, but not yet a derivation of literal hydrodynamic turbulence from the full PDE.

So the next real falsification step is:

> compute or model the actual rigid-mouth orbit packet
> \[
> \bigl(\delta\ln R_{\rm tr},\delta\ln\mathfrak N_*,\delta\ln\epsilon_\eta\bigr)
> \]
> strongly enough to decide whether the realized branch clears or exceeds the static gate above.

That is where the present stack now says the answer lives.

---

## 9. SymPy-backed status

The accompanying audit script verifies all of the following:

1. the exact Stage-171 observable compiler
   \[
   \Theta_1=\delta\ln R_{\rm tr},
   \qquad
   \Xi_1=\delta\ln\mathfrak N_* - B_*\,\delta\ln R_{\rm tr},
   \qquad
   \mathcal R_1=
   -\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta-\Xi_1;
   \]
2. the track-locked specialization
   \[
   \delta\ln R_{\rm tr}=0
   \Longrightarrow
   \Theta_1=0,\quad \Xi_1=\delta\ln\mathfrak N_*;
   \]
3. the exact prefactor identity
   \[
   \Xi_{\rm load}
   =
   \frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}
   =
   \frac{P_1}{P_0};
   \]
4. the stronger operator-rigidity specialization
   \[
   D_{01}=0
   \Longrightarrow
   \Xi_{\rm load}=\frac{N_{01}}{N_0};
   \]
5. the transported static ceiling
   \[
   |\epsilon\Xi_1|
   \le
   \frac{P_{\rm crit}\hat m_0^{\,2}}{\Delta_{\rm norm}+T_{\rm quad}}-1;
   \]
6. its calibrated-branch simplification
   \[
   \Delta_{\rm norm}=0
   \Longrightarrow
   |\epsilon\Xi_1|
   \le
   \frac{P_{\rm crit}\hat m_0^{\,2}}{T_{\rm quad}}-1;
   \]
7. the equivalent \(\bar P_0\)-form
   \[
   \bar P_0=\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}
   \Longrightarrow
   |\epsilon\Xi_1|
   \le
   \frac{P_{\rm crit}}{\bar P_0}-1;
   \]
8. and the numerical recovery of the two carried Stage-207 budgets
   \[
   0.367930328492646,
   \qquad
   0.737619063660757
   \]
   from the Stage-206 compatibility-point value
   \[
   \bar P_0 \approx 0.002069792318062885.
   \]

So the note is not just a verbal reinterpretation.
It is backed by an executable reconstruction of the rigid-mouth packet compiler and its static gate.

---

## 10. What the next honest stage should do

The next stage should **not** invent more support algebra.

The honest next theorem gate is now:

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
