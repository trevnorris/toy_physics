# Closed-book 0PN referee review

Scoped examination of the repaired 0PN requirements, wall identities, verification record, and two-mouth release procedure. Forces on the open mouth, a wall patch, and the carrier are kept separate. No gravitational law or interface closure is assumed.

---

## 1. Independent derivation (wall patch only)

**Constitutive setup (given).** Normal matching \(v_n=V_n\) on the wall. Tangential balance
\[
K v+\gamma(v-V_\parallel)=T,
\]
with \(K>0\), \(\gamma>0\), and Rayleigh density \(\gamma(v-V_\parallel)^2/2\).

**Solved material trace**
\[
v=\frac{T+\gamma V_\parallel}{K+\gamma}.
\]
Matches the printed solved trace.

**Wall traction from the slip term**
\[
f_\parallel=\gamma(v-V_\parallel)=\frac{\gamma}{K+\gamma}(T-K V_\parallel)=-\frac{\gamma(KV_\parallel-T)}{K+\gamma}.
\]
Matches the printed solved wall force. On the archived \(T=0\) branch,
\[
f_\parallel=-\frac{K\gamma}{K+\gamma}V_\parallel
\]
(opposite the imposed wall motion).

**Combined slip power on wall plus material.** With action–reaction \(f_{\mathrm{mat}}=-f_\parallel\),
\[
f_\parallel V_\parallel+f_{\mathrm{mat}}v=-\gamma(v-V_\parallel)^2.
\]
Matches the printed power identity. Dissipation is non-positive as required by the Rayleigh form.

**Local slopes (wall component traction)**
\[
\partial_U f_\parallel=\frac{\gamma\,T_U}{K+\gamma},\qquad
\partial_V f_\parallel=-\frac{\gamma(K-T_V)}{K+\gamma},
\]
\[
\partial_U f_\parallel+\partial_V f_\parallel=-\frac{\gamma(K-T_U-T_V)}{K+\gamma}.
\]
Sum vanishes iff \(T_U+T_V=K\) at fixed coefficients. The printed relative-flow control
\[
\bigl[K\gamma/(K+\gamma),\;-K\gamma/(K+\gamma),\;0\bigr]
\]
is exactly the case \(T_U=K\), \(T_V=0\). The changed-incident case with nonzero sum is the correct contrast. Wall residuals at zero are consistent with the algebra.

**Falsification attempts on these identities.** No algebraic counterexample under the stated constitutive assumptions. The identities hold for a **wall patch**. They do not by themselves fix a force on an open aperture or on the geometric center.

---

## 2. Separation of forces (mouth / wall / carrier)

| Object | What the sourced law fixes | What remains open |
| --- | --- | --- |
| Wall patch (E5) | Normal matching; tangential \(f_\parallel\) and slip power as above | Mapping from distant companion flow to local \(T\); transmission of \(f_\parallel\) to the rim |
| Open aperture | Required \(D(a)>0\) (nonzero relative normal current) | Sideways momentum transfer; stress jump; rim reaction; storage/exchange outside the translating body |
| Geometric center / carrier | Native solver can prescribe center motion | Net drive \(F_x(U,V)\); inertia response \(\chi_{ij}\); absolute supports \(H_i\) |

**Hard scope fact.** Applying E5 normal matching on the **open** aperture forces zero relative normal current there, which contradicts \(D(a)>0\). E5 can describe a wall beside a draining opening. It does not specify the sideways force retained by the draining aperture. Permeable option E3 leaves crossing velocity free relative to the pattern and likewise does not close the aperture force. Prior boundary adjudication selected no endpoint. **No sourced term in the record fixes the net physical force on the open mouth or its center.**

**Carried momentum.** In the mouth frame the ambient-frame carried influx shifts by \(-D(a)V\) at linear order. For equal crossing mass rates that shift cancels between the two carried contributions. On the stationary equal-current single-carrier branch with one continuous differentiable tangential phase, the carried tangential difference vanishes. Cancellation removes the carried piece only. Stress difference, rim forces, and exterior exchange remain. A purely stress-based center force is therefore admissible only after that cancellation is stated; it still does not determine the value of the stress/rim contribution. The record’s conditioning on carried cancellation is correct and incomplete as a closure.

---

## 3. 0PN response coefficients \(A\), \(B\)

Definitions as given: ambient-frame companion velocity \(U\) at the target location with target absent; center velocity \(V\); drive \(F_x\) on the center (at steady prescribed translation, minus the external holding force). Preparations: held center for \(A=\partial_U F_x(0,0)\); isolated translation through resting ambient for \(B=\partial_V F_x(0,0)\); odd difference over twice the speed; half-speed repeat; inward drain throughout.

**Conditional calculus (sound).** If \(F_x\) depends only on \(U-V\) at leading order, then \(B=-A\). If also \(A>0\), an isolated mover feels \(-A\mathbf{V}\) and loses mechanical energy at rate \(A|\mathbf{V}|^2\) to leading order. That conditional drag would enter orbital energy/angular momentum bookkeeping. The Newtonian mechanism still needs positive radial response and only subleading net tangential acceleration; subtracting an isolated radial control does not remove a physical pair torque.

**What does not force \(B=-A\).**
- A common boost of mouth, incident material, and distant reservoir preserves velocities relative to that reservoir and supplies **no** identity between \(A\) and \(B\).
- Held periodic feed sets a phase gradient; isolated drift translates the intake with the outer bath fixed. Those are different preparations.
- Writing carried momentum in the mouth frame cancels as above and still supplies no value for total drive.
- Wall slope sum \(\partial_U f+\partial_V f=0\) is a **component wall** statement under \(T_U+T_V=K\). Inferring \(B=-A\) for the geometric center additionally requires rim transmission and the full opening momentum balance. None of that is sourced.
- Archived isolated exterior (brane-side) drift loads changed sign across successive fixed-thickness refinements (minus, plus, minus). They do not fix either physical slope.
- Present boundary data determine neither slope nor \(A+B=0\).

**Falsification of the kinematic “relative-flow” claim for the mouth.** The claim that the mouth must sense only \(U-V\) is not entailed by observer equivalence, carried cancellation, or the wall slip law. It remains an open physical question. The wall relative-flow discriminant is a local check on \(T\), not a mouth-center identity.

---

## 4. Two-mouth release from rest

Protocol structure (positions, velocities, independent field data, history, intake/reservoir laws preserved; supports removed without impulse; accelerations, \(\partial_t\) fields, and constraint stresses re-solved; unsteady pressure allowed to jump) is formally consistent with instantaneous release. Axis from mouth 1 to 2. Release formula using **absolute** total supports \(H_i\) and instantaneous coupled response \(\chi_{ij}\) is the correct linear onset relation
\[
\alpha_i(0^+)=-\sum_j\chi_{ij}H_j
\]
under smooth force-to-acceleration response and zero acceleration while the same \(H_j\) are applied. Freezing constraints at release would be wrong; the record’s unfrozen independent initial data are appropriate. At the prepared rest instant the \(BV\) drift piece vanishes, so release cannot settle the later drag question.

**Link status (all still open for physical mouths).**

1. **Radius → inward rate \(D(a)\).** Archived boundary data admit different prescribed rates at fixed radius. Radius does not select the flow. Open.
2. **Source → incident flow.** Conditional exterior monopole \(U_{2\to1}=+Q_2/(4\pi d^2)\), \(U_{1\to2}=-Q_1/(4\pi d^2)\) under steady spherical flux and uniform density. Sign convention along the 1→2 axis is consistent with inflow toward a drain. Intervening material budget linking physical intake to exterior strength remains open.
3. **Incident flow → retained force.** Need total support forces \(H_i\) or integrated net transfer (carried + stress + rim) on each translating opening. Archived unequal pinned mouths supply brane-side incoming loads on both mouths and **neither** total support reaction. Open.
4. **Force → acceleration.** Need \(\chi_{ij}\) from independent weak probes with both centers free. Native solver prescribes center motion. A separate finite-companion benchmark shows the readout is possible when outlet and inertia are specified. Physical \(\chi_{ij}\) unmeasured. Open.

Under the further **ifs** (vanishing isolated gravitational drive by symmetry, \(\chi_{ii}\to 1/M_i\), subleading cross-acceleration, local incident velocity supplying the leading interaction), one would obtain
\[
\alpha_1=+A(a_1)\frac{Q_2}{4\pi M_1 d^2},\qquad
\alpha_2=-A(a_2)\frac{Q_1}{4\pi M_2 d^2}
\]
at leading order. With present inputs, **neither** initial acceleration has a physical numerical prediction. Comparison of \(d^2\alpha_1/M_2\) and \(-d^2\alpha_2/M_1\) to one positive \(G\) is not available.

---

## 5. Verification record (as printed)

Reported checks: source surface and normal matching; reconstruction of \(K v+\mathrm{d}R/\mathrm{d}v\); archived 3-by-3 homogeneous tangent coefficient; forced wall reaction vs \(T-Kv\); implicit differentiation of stationarity in \(U\) and \(V\) with slope and sum checks; each actual/reference pair consumed by the guard.

Printed residuals are identically zero on every listed wall identity. Forcing discriminant values match the independent algebra. **Coverage is wall outputs only.** Independent local form mutations (including deletions) are out of scope for this review, as stated.

**Defects that affect claims.**
- The record may not claim mouth-center \(A\), \(B\), or \(A+B=0\) from wall residuals.
- Conditioning a stress-based center expression on carried cancellation is necessary and not sufficient.
- Release language using absolute supports and coupled response is procedurally sound; it does not supply the missing \(H_i\) or \(\chi_{ij}\).
- No printed evidence addresses opening momentum balance, rim transmission, \(D(a)\) selection, or pair supports.

---

## 6. Findings by severity

### Critical
1. **No sourced net force on the open mouth / geometric center.** Wall E5 and permeable E3 do not close aperture sideways force. Exact correction: treat wall \(f_\parallel\) as wall-patch only; require an explicit mouth-level momentum-transfer law or independent total-force readings in both \(A\) and \(B\) preparations before any center-force claim.
2. **E5 normal matching must not be imposed on the open aperture.** Exact correction: restrict normal matching to the wall beside the opening; keep \(D(a)>0\) on the aperture.
3. **Present data do not determine \(A\), \(B\), or \(A+B=0\).** Archived drift signs are refinement-unstable; boundary data leave both slopes open. Exact correction: delete or gate any claim that the mouth force is pure relative-flow drag until (i) physical \(T(U,V)\) or total-force measurements exist, (ii) rim transmission is specified, and (iii) opening stress/carried balance is closed.
4. **Two-mouth numerical \(\alpha_i\) are unavailable.** \(D(a)\), \(H_i\) (or transfer law), and \(\chi_{ij}\) are all open. Exact correction: strip any implication of an evaluated Newtonian \(G\) comparison; keep release formulas symbolic with open coefficients.

### Major
5. **Wall slope sum \(\leftrightarrow\) mouth \(B=-A\) is a non sequitur without transmission and opening balance.** Mutation: any executable or prose bridge from \(T_U+T_V=K\) to mouth \(A+B=0\) must be removed or gated on those extra identities.
6. **Brane-side loads ≠ total supports \(H_i\).** Mutation: unequal-mouth archive citations must not be used as \(H_i\) inputs to the release formula.
7. **Verification scope creep.** Mutation: keep guard language and pass/fail claims strictly as “wall outputs”; do not mark mouth readiness from wall residuals.
8. **Relative-flow discriminant is wall-local.** It validates algebra under assumed \(T_U,T_V\); it does not calibrate companion-driven \(T\) or mouth \(A\).

### Minor / procedural
9. Odd-difference and half-speed protocols for future \(A,B\) readings are sound.
10. Finite-companion weak-probe method is the right template for \(\chi_{ij}\) once free centers are allowed.
11. Release-from-rest correctly leaves the later drag question open.
12. Monopole \(Q/(4\pi d^2)\) must remain labeled conditional (steady spherical flux, uniform density, material budget).

---

## 7. Exact executable / mutation instructions

1. Retain the printed wall guard suite unchanged in scope; do not add mouth-center pass criteria to it.
2. Add (or require before readiness) either: total-force odd-difference jobs for \(A(a)\) and \(B(a)\) at matched far boundary, radius, and intake law; or a declared opening momentum-transfer law used by the center ODE.
3. For release: record absolute held-pair supports \(H_i\); measure \(\chi_{ij}\) with separate weak probes on both free centers in the same prepared state; then apply \(\alpha_i(0^+)=-\sum_j\chi_{ij}H_j\). Do not substitute brane-side incoming loads for \(H_i\).
4. Enforce a single physical \(D(a)\) so each radius selects its drain rate; reject archives that prescribe rate independently of radius as physical \(D(a)\).
5. Delete or hard-gate any code path or assertion that sets mouth \(B=-A\) from wall \(T_U+T_V=K\) alone.
6. Keep carried-cancellation as an explicit precondition before dropping carried terms from a stress-only center force; do not treat cancellation as determination of that force.
7. Leave orbital integration deferred; do not emit a numerical \(G\) test until \(A(a)\), masses/inertia, and \(Q_i(a)\) are fixed inputs.

---

## 8. What cannot be verified from this record

- Numerical values of \(A\), \(B\), \(K\), \(\gamma\), \(T_U\), \(T_V\), \(H_i\), \(\chi_{ij}\), or \(D(a)\) (only symbolic wall residuals are shown).
- That production code matches the identities beyond the printed actual/reference pairs.
- Local form mutations or deletions (assigned to another reviewer).
- Physical rim-transmission pathway and exterior material budget.
- Whether a Newtonian-sector comparison would succeed after inputs are fixed.
- Any executable evidence outside wall outputs (opening balance, pair supports, free-center probes, orbit).

---

## 9. Separate scoped verdicts

**Wall derivation.** Sound for the wall patch under the stated slip law, Rayleigh density, and normal matching. Force, power, solved trace, and slope-sum algebra check independently. The derivation correctly does **not** determine the net force on the open aperture or geometric center, and correctly does not force mouth \(A+B=0\) from kinematics or carried cancellation alone.

**Executable evidence as reported.** Internally consistent with the wall identities; printed residuals are zero; forcing discriminant matches the independent algebra. Evidence is **wall-scoped only**. It does not underwrite mouth \(A\), \(B\), rim transmission, total supports, or two-mouth release forces.

**Conditional two-mouth readiness.** **Not ready.** Missing physical \(D(a)\) selection, mouth-level force law or calibrated \((A,B)\) / absolute \(H_i\), and free-center inertia response \(\chi_{ij}\). The release protocol’s structure (unfrozen constraints, absolute supports, coupled instantaneous response) is acceptable as a future procedure once those inputs exist. No physical numerical initial accelerations and no \(G\) comparison are licensed by the present record.
