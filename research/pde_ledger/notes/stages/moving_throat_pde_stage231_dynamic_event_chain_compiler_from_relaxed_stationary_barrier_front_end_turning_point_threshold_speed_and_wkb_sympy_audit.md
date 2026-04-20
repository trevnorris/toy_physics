# Moving-Throat PDE — Stage 231: Dynamic Event-Chain Compiler from the Relaxed Stationary Barrier Front End, Turning Points, Threshold Speed, and WKB

## Status

**Exact within**

1. the carried Stage-230 stationary lowered-barrier front end
   \[
   V_{\rm eff}^{(230)}(r),
   \]
2. the reduced one-dimensional separation dynamics
   \[
   m_s\ddot r=-\partial_r V_{\rm eff}^{(230)}(r),
   \]
3. the standard turning-point / WKB reduction for a smooth single-peak barrier,
4. and the declared Session-II benchmark specialization.

This stage does **not** introduce a new same-charge kernel class.
It is the dynamic continuation of Stage 230, not a replacement for the static-first barrier audit. The short-range/open-system firewall remains unchanged: the dynamic event chain is built on the already-lowered short-range branch rather than on a new long-range attractive law or a new linear dynamic kernel class.

---

## Purpose

Stage 230 assembled the lowered stationary branch
\[
V_{\rm eff}^{(230)}(r)
=
V_{\rm short}^{(1p)}(r)
-
\lambda_L S_{\rm leak}(r)
-
\lambda_W\mathcal W_w^{\rm sess}(r)
-
\Delta E_{UV}(r)
-
\mathcal M_\sigma(r).
\]

What was still missing was the **dynamic event chain** that turns that stationary front end into the objects actually used by the Session-II scattering test:

1. the barrier peak \((r_{\rm peak},V_{\rm peak})\),
2. the finite-radius classical threshold speed,
3. the subbarrier turning points,
4. the WKB action and transmission factor,
5. the dynamic turning-point diagnostics carried by the same lowered branch,
6. and the exact window condition that says when the lowered branch reaches contact while the pure Coulomb reference still turns back.

That is exactly what Stage 231 does.

Script-backed status:
- `scripts/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
  checks the symbolic event-chain conservation laws, threshold-speed and WKB
  compilers, the Coulomb reference formulas, the near-top action normal form,
  and the declared Session-II benchmark specialization.
- `mathematica/moving_throat_pde_stage231_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
  mirrors the same symbolic compiler path in the second CAS and keeps the
  Session-II benchmark numbers explicitly confined to the benchmark-only
  readback layer.

---

## Provenance

This stage sits directly after Stage 230 and uses the same claim-status discipline.

- The stationary lowering already came from the one-port short-range baseline plus the Stage-227 leakage/work packet, the Stage-228 non-rigid `U/V` drain packet, and the Stage-229 compensated-source packet.
- The same-charge barrier audit had already shown that the linear dynamic mixed bundle does **not** generate a new conservative kernel class at large distance. Linear time dependence only makes the already-admissible short-range families frequency dependent, while the first outgoing correction is phase lag / pumping rather than direct conservative barrier lowering.
- Session II then tested what the lowered stationary front end does dynamically inside a reduced scattering problem.

So Stage 231 is the correct next step in the audit trail: it keeps the static-first theorem logic intact and promotes the already-lowered front end into explicit event-chain formulas.

---

## 0. Why this stage is needed

Before this step, the derivation stack had a fully explicit **stationary** lowered branch but no corresponding dynamic compiler.

That meant the stack still lacked direct formulas for statements of the form

- “the lowered branch becomes classically traversable at finite launch speed,”
- “the outer turning point moves inward,”
- “the WKB exponent is reduced relative to Coulomb,”
- and “the turning-point branch carries a definite dynamic gradient trigger.”

Stage 231 closes that gap.

---

## 1. Dynamic lift of the Stage-230 stationary front end

Take the Stage-230 reduced barrier as the one-dimensional same-charge front end
\[
V(r):=V_{\rm eff}^{(230)}(r).
\]

The reduced radial dynamics is
\[
\boxed{
 m_s\ddot r=-V'(r).
}
\]

### 1.1 Exact energy integral

Define the reduced mechanical energy
\[
\boxed{
E=\frac{m_s}{2}\dot r^{\,2}+V(r).
}
\]
Then
\[
\frac{dE}{dt}=m_s\dot r\ddot r+V'(r)\dot r
=\dot r\bigl(m_s\ddot r+V'(r)\bigr)=0.
\]
So the event chain is completely controlled by the stationary front end once the launch energy is fixed.

### 1.2 Launch energy from finite initial radius

If the reduced trajectory is launched at radius \(r_0\) with inward speed \(v_0\), its energy is
\[
\boxed{
E_0(v_0)=\frac{m_s}{2}v_0^2+V(r_0).
}
\]
This is the quantity used throughout the rest of the compiler.

---

## 2. Barrier peak and classical-threshold compiler

Assume the lowered branch has one barrier peak between the launch point and the chosen contact scale:
\[
V'(r_{\rm peak})=0,
\qquad
V''(r_{\rm peak})<0,
\qquad
V_{\rm peak}:=V(r_{\rm peak}).
\]

### 2.1 Finite-radius classical threshold speed on the lowered branch

The minimum launch speed required to reach the top of the reduced barrier is determined by
\[
E_0(v_{\rm crit,new})=V_{\rm peak}.
\]
So
\[
\boxed{
 v_{\rm crit,new}
 =
 \sqrt{\frac{2}{m_s}\bigl(V_{\rm peak}-V(r_0)\bigr)}.
}
\]

This is the first exact dynamic compiler output of the lowered branch.

### 2.2 Pure-Coulomb contact threshold at the same launch radius

For the reference potential
\[
V_{\rm Coul}(r)=\frac1r,
\]
the energy needed to reach a chosen contact radius \(r_{\rm contact}\) from the same \(r_0\) is
\[
E_{\rm Coul}(v_0)=\frac{m_s}{2}v_0^2+\frac1{r_0}.
\]
The corresponding threshold speed is therefore
\[
\boxed{
 v_{\rm contact,Coul}
 =
 \sqrt{\frac{2}{m_s}\left(\frac1{r_{\rm contact}}-\frac1{r_0}\right)}.
}
\]

### 2.3 Exact classical-window theorem

Let the lowered branch be single-peaked on \([r_{\rm contact},r_0]\), so that crossing the peak is sufficient to reach contact. Then if
\[
\boxed{
 v_{\rm crit,new}<v_0<v_{\rm contact,Coul},
}
\]
one has simultaneously
\[
E_0(v_0)>V_{\rm peak},
\qquad
E_{\rm Coul}(v_0)<\frac1{r_{\rm contact}}.
\]
So the lowered branch reaches the contact scale while the pure Coulomb branch still turns back outside it.

This is the cleanest exact statement of the Session-II classical corridor.

---

## 3. Turning-point compiler and subbarrier WKB event chain

Now fix a subbarrier energy \(E\) satisfying
\[
V(r_0)<E<V_{\rm peak}.
\]

### 3.1 Turning points

The classically forbidden interval is bounded by the turning points
\[
\boxed{
V\bigl(r_-(E)\bigr)=E,
\qquad
V\bigl(r_+(E)\bigr)=E,
\qquad
r_-(E)<r_+(E).
}
\]
Here \(r_+(E)\) is the outer turning point and \(r_-(E)\) is the inner turning point.

Differentiating the turning-point equation gives the exact transport law
\[
\boxed{
\frac{dr_\pm}{dE}=\frac{1}{V'\bigl(r_\pm(E)\bigr)}.
}
\]
So on the usual barrier geometry,

- \(V'(r_+(E))<0\), hence \(dr_+/dE<0\): the outer turning point moves inward as energy rises;
- \(V'(r_-(E))>0\), hence \(dr_-/dE>0\): the inner turning point moves outward as energy rises.

Therefore the forbidden width shrinks as \(E\) approaches the peak from below.

### 3.2 Launch speed for a chosen subbarrier energy

If one wants the reduced trajectory to have energy \(E_{\rm sub}\) at launch radius \(r_0\), then
\[
\boxed{
 v_{0,\rm sub}(E_{\rm sub})
 =
 \sqrt{\frac{2}{m_s}\bigl(E_{\rm sub}-V(r_0)\bigr)}.
}
\]

### 3.3 Exact WKB action and transmission factor

The reduced WKB action is
\[
\boxed{
 I_{\rm new}(E)
 =
 \frac{1}{\hbar_{\rm eff}}
 \int_{r_-(E)}^{r_+(E)}
 \sqrt{2m_s\bigl(V(r)-E\bigr)}\,dr.
}
\]
The transmission factor is
\[
\boxed{
 T_{\rm new}(E)=e^{-2I_{\rm new}(E)}.
}
\]

Differentiating under the integral sign and using the vanishing boundary terms at the turning points gives
\[
\boxed{
\frac{dI_{\rm new}}{dE}
=
-\frac{\sqrt{m_s/2}}{\hbar_{\rm eff}}
\int_{r_-(E)}^{r_+(E)}
\frac{dr}{\sqrt{V(r)-E}}
<0.
}
\]
So any mechanism that lowers the barrier or raises the incident energy decreases the action and increases the transmission exponentially.

### 3.4 Pure-Coulomb reference action

For the pure Coulomb reference
\[
V_{\rm Coul}(r)=\frac1r,
\]
the outer turning point is exact:
\[
\boxed{
 r_{\rm turn,Coul}(E)=\frac1E.
}
\]
If the inner comparison point is the chosen contact radius \(r_{\rm contact}\), then
\[
\boxed{
 I_{\rm Coul}(E)
 =
 \frac{\sqrt{2m_s}}{\hbar_{\rm eff}}
 \int_{r_{\rm contact}}^{1/E}
 \sqrt{\frac1r-E}\,dr.
}
\]
The integral closes in elementary form for \(Er_{\rm contact}<1\):
\[
\boxed{
I_{\rm Coul}(E)
=
\frac{\sqrt{2m_s}}{\hbar_{\rm eff}}
\left[
\frac{\pi}{2\sqrt{E}}
-
\sqrt{r_{\rm contact}\bigl(1-Er_{\rm contact}\bigr)}
-
\frac{\arcsin\sqrt{Er_{\rm contact}}}{\sqrt{E}}
\right].
}
\]

### 3.5 Exact transmission-ratio compiler

The reduced-vs-Coulomb tunneling ratio is therefore
\[
\boxed{
\frac{T_{\rm new}(E)}{T_{\rm Coul}(E)}
=
\exp\!\bigl[-2\bigl(I_{\rm new}(E)-I_{\rm Coul}(E)\bigr)\bigr].
}
\]
So once the lowered branch gives \(I_{\rm new}<I_{\rm Coul}\), the enhancement is completely fixed.

---

## 4. Near-top normal form

For energies just below the peak it is useful to expand
\[
V(r)=V_{\rm peak}-\frac{K_{\rm peak}}{2}(r-r_{\rm peak})^2+O\bigl((r-r_{\rm peak})^3\bigr),
\qquad
K_{\rm peak}:=-V''(r_{\rm peak})>0.
\]
Define
\[
\Delta E:=V_{\rm peak}-E>0.
\]
Then the turning points are
\[
\boxed{
 r_\pm(E)
 =
 r_{\rm peak}
 \pm
 \sqrt{\frac{2\Delta E}{K_{\rm peak}}}
 +O\bigl(\Delta E\bigr).
}
\]
The leading WKB action is
\[
\boxed{
 I_{\rm top}(E)
 =
 \frac{\pi\Delta E}{\hbar_{\rm eff}}\sqrt{\frac{m_s}{K_{\rm peak}}}
 +O\bigl(\Delta E^{3/2}\bigr).
}
\]

So near the top the event chain collapses to one curvature scale \(K_{\rm peak}\).

This normal form is not needed for the Session-II benchmark itself, but it is the right local compiler for later near-threshold analysis.

---

## 5. Dynamic turning-point diagnostics carried forward

Even before the magnetic/helicity branch is added, the turning-point event chain naturally carries two reduced diagnostics that are useful later.

### 5.1 Dynamic barrier scalar at the outer turning point

The same weak-axisymmetric barrier scalar from the stationary front end may be sampled dynamically on the event chain:
\[
\boxed{
\Xi_{\rm turn}(E):=\Xi_1\bigl(r_+(E)\bigr).
}
\]
This quantity is not a new theorem by itself, but it tells the later weak-axisymmetric audit where the subbarrier path is sitting on the transported front-end packet.

### 5.2 Trigger-width diagnostic

The Session-II gradient trigger was defined by
\[
\chi_\lambda\equiv \lambda\,\bigl|\partial_r\ln V(r)\bigr|.
\]
Solving the threshold condition \(\chi_\lambda=1\) at the outer turning point gives
\[
\boxed{
\lambda_{\rm th}(E)
=
\left|\frac{V\bigl(r_+(E)\bigr)}{V'\bigl(r_+(E)\bigr)}\right|
=
\left|\frac{E}{V'\bigl(r_+(E)\bigr)}\right|.
}
\]
So the same turning-point event chain also determines the first dynamic confinement-width trigger.

### 5.3 Scope boundary

Stage 231 stops here.
The aligned-vs-anti-aligned helicity-export diagnostic belongs to the next stage because it is not part of the scalar scattering compiler itself; it is an additional mixed/vortical diagnostic layered on top of this event chain.

---

## 6. Session-II benchmark specialization

The Session-II run used
\[
m_s=1,
\qquad
\hbar_{\rm eff}=1,
\qquad
r_0=5,
\qquad
r_{\rm contact}=0.18,
\qquad
E_{\rm sub}=2.5.
\]
The reported lowered-branch dynamic observables were
\[
r_{\rm peak}=0.23944389,
\qquad
V_{\rm peak}=3.42933112,
\qquad
V(r_0)=0.19999794,
\]
\[
r_{\rm turn,new}=0.39096144,
\qquad
r_{\rm inner}=0.19039548,
\qquad
I_{\rm new}=0.19744614,
\]
with the additional turning-point diagnostics
\[
\Xi_{\rm turn}=0.34437471,
\qquad
\lambda_{\rm th}=0.42826825.
\]

### 6.1 Peak-to-threshold compiler

Using the exact finite-radius threshold law,
\[
v_{\rm crit,new}
=
\sqrt{2\bigl(V_{\rm peak}-V(r_0)\bigr)}
=
\sqrt{2\bigl(3.42933112-0.19999794\bigr)}
\approx 2.54139063.
\]
This reproduces the Session-II threshold value.

For the pure Coulomb contact comparison,
\[
v_{\rm contact,Coul}
=
\sqrt{2\left(\frac1{0.18}-\frac15\right)}
\approx 3.27278339.
\]
So the reduced classical corridor is explicit:
\[
\boxed{
2.54139063<v_0<3.27278339.
}
\]

### 6.2 Subbarrier launch speed and Coulomb turning point

At the fixed subbarrier energy \(E_{\rm sub}=2.5\), the launch speed is
\[
v_{0,\rm sub}
=
\sqrt{2\bigl(2.5-0.19999794\bigr)}
\approx 2.14476202.
\]
The pure-Coulomb outer turning point is
\[
r_{\rm turn,Coul}=\frac1{2.5}=0.4,
\]
which agrees with the reported Coulomb reference value \(0.40000141\) up to the expected numerical rounding of the session scan.

### 6.3 WKB enhancement

The reported Coulomb comparison action was
\[
I_{\rm Coul}\approx 0.30222297.
\]
Using the exact transmission law,
\[
T_{\rm new}=e^{-2I_{\rm new}}=e^{-2(0.19744614)}\approx 0.67375262,
\]
\[
T_{\rm Coul}=e^{-2I_{\rm Coul}}=e^{-2(0.30222297)}\approx 0.54637707.
\]
So
\[
\frac{T_{\rm new}}{T_{\rm Coul}}
=
\exp\!\bigl[-2(I_{\rm new}-I_{\rm Coul})\bigr]
\approx 1.23312756.
\]
That is a transmission increase of
\[
\boxed{
\left(\frac{T_{\rm new}}{T_{\rm Coul}}-1\right)\times 100\%
\approx 23.3128\%.
}
\]

### 6.4 Above-threshold contact demonstration

The Session-II report also gave an explicit above-threshold demonstration speed
\[
v_{0,\rm cross}=2.59221845.
\]
This lies inside the exact window:
\[
2.54139063<2.59221845<3.27278339.
\]
So the lowered branch is above its own barrier while the pure Coulomb branch is still below the energy needed for the same contact radius.

Using the Coulomb energy at that speed,
\[
E_{\rm Coul}(v_{0,\rm cross})
=
\frac12 v_{0,\rm cross}^2+\frac15,
\]
the pure-Coulomb turning point is
\[
r_{\rm turn,Coul}(v_{0,\rm cross})
=
\frac{1}{E_{\rm Coul}(v_{0,\rm cross})}
\approx 0.2809,
\]
consistent with the reported value \(0.28091705\).

So the benchmark does exactly what Stage 231 needs it to do: it exhibits a concrete speed range where the lowered branch reaches the chosen contact scale and the Coulomb comparison does not.

### 6.5 Coulomb closed-form check

The exact Coulomb WKB formula from §3.4 gives, with \(m_s=\hbar_{\rm eff}=1\),
\[
I_{\rm Coul}(2.5;0.18)
=
\sqrt2\left[
\frac{\pi}{2\sqrt{2.5}}
-
\sqrt{0.18\,(1-2.5\times 0.18)}
-
\frac{\arcsin\sqrt{2.5\times 0.18}}{\sqrt{2.5}}
\right]
\approx 0.30230580.
\]
That is within the numerical tolerance implied by the reported Session-II discrete reference value \(0.30222297\), so the closed-form Coulomb side of the compiler is consistent with the dynamic benchmark.

---

## 7. What Stage 231 proves

Stage 231 upgrades the relaxed barrier branch from a stationary softening plot into a real dynamic event chain.

It proves, within the declared reduced closure:

1. the exact reduced motion is still energy-conserving once the Stage-230 front end is fixed,
2. the barrier peak determines a finite launch threshold
   \[
   v_{\rm crit,new}=\sqrt{2\bigl(V_{\rm peak}-V(r_0)\bigr)/m_s},
   \]
3. there is an exact comparison threshold for the same contact scale on the pure Coulomb branch,
4. the subbarrier branch is controlled by explicit turning points and an explicit WKB action,
5. the transmission improvement is exactly
   \[
   \exp\!\bigl[-2(I_{\rm new}-I_{\rm Coul})\bigr],
   \]
6. and the turning-point path naturally carries the dynamic diagnostics
   \[
   \Xi_{\rm turn},\qquad \lambda_{\rm th},
   \]
   needed by the next stages.

So the lowered same-charge branch is no longer just “a smaller plotted barrier.”
It is a complete reduced scattering object with explicit classical and subbarrier compilers.

At the same time, this stage **does not** reopen the barrier-audit verdict.
The event-chain improvement still sits on a short-range/open-system front end, and the linear dynamic mixed bundle still contributes phase-lag / pumping rather than a new conservative kernel class.

---

## 8. Immediate next step

The next stage should attach the magnetic/helicity diagnostic to this event chain.
That means:

1. keep the dynamic turning-point and WKB compiler derived here,
2. add the aligned-vs-anti-aligned mixed/vortical export observable,
3. and test whether the preferred branch is the one that most effectively unloads unresolved repulsive structure while traversing the same lowered event chain.

That is exactly the Stage-232 job.
