## Closed-book physics peer review: 0PN requirements and two-mouth release record

### Direct answers

**Force conjugate to wall velocity.**  
From the Rayleigh density \(R=\gamma(v-V_\parallel)^2/2\), the traction conjugate to wall velocity is
\[
f_\parallel=\gamma(v-V_\parallel)=\frac{\gamma}{K+\gamma}\bigl(T-K V_\parallel\bigr).
\]
(The material feels the opposite traction \(-\,f_\parallel\).) On the archived branch \(T=0\) this is \(f_\parallel=-[\gamma K/(K+\gamma)]V_\parallel\), opposite the imposed wall motion.

**Mechanical power of the two opposite forces.**  
Wall power \(f_\parallel V_\parallel\) plus material power \((-f_\parallel)v\) sums to
\[
-\gamma(v-V_\parallel)^2=-2R.
\]
This is an identity of the slip term alone.

**Relation to force on the geometric center.**  
None is fixed by the supplied equations. Wall-patch tractions are component stresses on the throat wall. Net drive on the geometric center still requires transmission of that traction to the rim plus the full momentum balance on the open aperture (carried momentum, two-sided stress difference, independent rim reaction). No wall-to-aperture matching relation is supplied, so \(f_\parallel\) does not determine \(F_x\), \(A(a)\), or \(B(a)\).

**Signed two-body release rule.**  
Lab axis from mouth 1 to mouth 2; supports removed without impulse at a settled held state. Linear response:
\[
\alpha_i(0^+)=-\sum_j\chi_{ij}H_j.
\]
Under the additional far-zone hypotheses (incident velocity supplies the whole leading drive, \(\chi_{ii}\to 1/M_i\), cross-acceleration subleading, \(V=0\) at release so the \(BV\) term drops),
\[
\alpha_1=+\frac{A(a_1)Q_2}{4\pi M_1 d^2},\qquad
\alpha_2=-\frac{A(a_2)Q_1}{4\pi M_2 d^2},
\]
with the exterior assignments \(U_{2\to1}=+Q_2/(4\pi d^2)\) and \(U_{1\to2}=-Q_1/(4\pi d^2)\) conditional on steady spherical monopole flow at uniform density.

**Independent inputs that remain unknown.**  
1. Physical intake law \(D(a)\) (radius does not select flow in the archived boundary data).  
2. Mouth-intake to exterior monopole strength \(Q_i\) (intervening material budget).  
3. Local exterior traction map \(T(U_\parallel,V_\parallel)\) (or \(\partial_{U_\parallel}T\), \(\partial_{V_\parallel}T\)).  
4. Whole-mouth linear responses \(A(a)=\partial_U F_x(0,0)\) and \(B(a)=\partial_V F_x(0,0)\), or an equivalent mouth-level momentum-transfer law.  
5. Rim transmission / opening momentum closure from wall traction to center force.  
6. Held-pair total support forces \(H_i\) (archives report brane-side loads only).  
7. Inertia and cross-response \(\chi_{ij}\) (and draining-state masses \(M_i\)) for free centers.  
8. Surface law and tangential exchange with receiving material on the open aperture.

---

### Independent derivation (wall patch)

Declared scalar wall equation:
\[
(K+\gamma)v-\gamma V_\parallel=T.
\]
Solution:
\[
v=\frac{T+\gamma V_\parallel}{K+\gamma},\qquad
v-V_\parallel=\frac{T-K V_\parallel}{K+\gamma}.
\]
Rayleigh derivatives:
\[
\frac{\partial R}{\partial v}=\gamma(v-V_\parallel),\qquad
\frac{\partial R}{\partial V_\parallel}=-\gamma(v-V_\parallel).
\]
Dissipative forces (standard Rayleigh sign: force \(=-\partial R/\partial\dot q\)) recover fluid traction \(-\gamma(v-V_\parallel)\) and wall traction \(+\gamma(v-V_\parallel)\), i.e. equal-and-opposite pair with power sum \(-2R\).

Substitute \(T=T_U U+T_V V\) (component notation suppressed). Wall-force slopes:
\[
\partial_U f_\parallel=\frac{\gamma}{K+\gamma}\,T_U,\qquad
\partial_V f_\parallel=\frac{\gamma}{K+\gamma}\,(T_V-K).
\]
Their sum vanishes iff
\[
T_U+T_V=K
\]
at fixed \(K,\gamma\). That is a condition on the omitted forcing, not a consequence of the slip operator. Boost of mouth + material + reservoir together leaves both reservoir-relative velocities unchanged and supplies no identity \(A=-B\) for the center force. Carried-momentum rewrite in the mouth frame shifts the ambient reading by \(-D(a)V\) at linear order; for equal crossing mass rates that shift cancels between incoming and transmitted contributions and does not fix the net drive. Stress difference and rim force remain open.

Falsification attempts that fail to break the slip identities: (i) change of observer; (ii) carried-momentum frame shift at equal rates; (iii) enforcing normal material–wall velocity agreement on the open aperture (that step contradicts \(D(a)>0\) and so cannot be used to close the mouth force). Falsification that succeeds against any claim that E5 fixes the mouth: absence of rim/aperture closure.

---

### Findings by severity

**Critical**

1. **Wall traction ≠ center force.**  
   Exact defect: treating \(f_\parallel\) or its slopes as \(F_x\), \(A\), or \(B\).  
   Correction: keep wall-patch identities scoped to the wall; require an independent mouth momentum law or total-force readings in both preparations before quoting \(A(a)\), \(B(a)\).  
   Executable: in any code path that prints “force slopes,” assert that the reported object is wall traction; refuse to export it as center \(F_x\) unless a separate rim-transmission module is present and tested.

2. **No total support force \(H_i\).**  
   Archives give brane-side loads on unequal pinned mouths, not the mechanical support reaction needed for release.  
   Correction: measure or define \(H_i\) as the total external force maintaining each held center; then \(\alpha_i(0^+)=-\sum_j\chi_{ij}H_j\).  
   Executable: add paired readouts (total support force on each center) in the held-pair preparation; reject release predictions that ingest only brane-side loads.

3. **\(D(a)\) and radius→intake unset.**  
   Archived boundary data admit different prescribed rates at fixed radius; radius does not select the flow.  
   Correction: mark every chain from \((a_1,a_2)\) to \(\alpha_i\) as open until a unique physical \(D(a)\) is fixed and matched to exterior \(Q_i\).

4. **\(A+B=0\) is not implied.**  
   Present boundary data determine neither slope nor \(A+B=0\). The slip-level condition \(T_U+T_V=K\) is necessary for opposite *wall* slopes and is still insufficient for the center.  
   Correction: strike any inference that E5 or boost invariance yields \(B=-A\) for the mouth.  
   Executable discriminant: hold \(K,\gamma\) fixed; mutate only \((T_U,T_V)\) so that \(T_U+T_V\neq K\); recompute wall slopes from the formula above; re-run the guarded checks (trace equation, opposite wall/fluid forces, power \(=-2R\)). Those three checks still pass. Separately mutate a reported center-slope pair to enforce \(A=-B\) while leaving \(T_U+T_V\) untouched: the slip checks still pass. Therefore passage of the guarded checks does not validate \(A=-B\) or any center-force claim.

**High**

5. **Guarded checks are slip-internal only.**  
   Enforcing the scalar trace, opposite wall/fluid forces, and power \(=-2R\) validates the Rayleigh–slip algebra. It does not validate mouth \(A\), \(B\), rim transmission, or release accelerations.  
   Correction: split test suites into “wall-patch identities” vs “center-force / release” and require the latter to fail closed until \(H_i\), \(\chi_{ij}\), and a mouth law are supplied.

6. **Normal matching on the open aperture is inconsistent with drainage.**  
   Option E5 over the aperture forces zero relative normal current, contrary to \(D(a)>0\). E5 can only describe wall beside an opening; E3 leaves crossing velocity free and still does not specify sideways aperture force.  
   Correction: do not apply E5 as an aperture closure; keep sideways mouth response as an independent requirement.

7. **Archived isolated drift sign flips (minus, plus, minus) under lateral refinement.**  
   Brane-side drift readouts are not a stable determination of \(B\).  
   Correction: treat archived drift slopes as inconclusive; replace with the prescribed preparation for \(B(a)\) (isolated translation through resting ambient, odd difference, half-speed repeat) once total center force is defined.

**Medium**

8. **Exterior monopole link is conditional.**  
   \(U=Q/(4\pi d^2)\) assumes steady dominant spherical flux and uniform density, and still needs intake↔exterior strength matching. Sign convention for \(Q_i\) must stay consistent with inward drainage and the stated \(U_{2\to1}=+Q_2/(4\pi d^2)\).  
   Correction: label the \(\alpha_i\) formulas as conditional on that exterior closure.

9. **Inertia matrix \(\chi_{ij}\) unmeasured for physical mouths.**  
   The finite-companion benchmark shows the readout is possible for specified boundaries; the native solver that *prescribes* center motion does not replace a free-center probe.  
   Correction: require independent weak-probe response with both centers free before converting \(H_j\) into \(\alpha_i\).

10. **Carried-momentum cancellation assumes equal crossing mass rates.**  
    If inflow and transmitted rates differ, the frame shift need not cancel; the record’s “supplies no value for the total drive” remains correct a fortiori, but the cancellation step should be stated as conditional on equal rates.

**Low / clarifications**

11. Half-speed odd-difference protocol for \(A\) and \(B\) is appropriate for linear response.  
12. Release from rest cannot constrain later drag (\(BV\) vanishes at the prepared instant); that statement is correct.  
13. Subtracting an isolated control in a radial comparison does not remove a physical pair torque; correctly flagged.  
14. Newtonian comparison needs positive radial response and only subleading *net* tangential acceleration; deferred orbit run cannot be used to close missing 0PN inputs.

---

### What cannot be verified from this record

- Numerical values or code output of the original guarded calculation (only that it passed the three slip identities).  
- Whether any printed slopes in that calculation satisfied \(T_U+T_V=K\).  
- Actual solver implementation of boundary options E3/E5, feed initialization, or support removal.  
- Sign-stable continuum limit of the archived drift fluxes beyond the reported minus/plus/minus sequence.  
- Any wall-to-aperture matching (explicitly absent).  
- Free two-mouth trajectories, work, and torque (orbit deferred; inputs incomplete).  
- Equality of the Newtonian \(G\) extracted from \(d^2\alpha_1/M_2\) and \(-d^2\alpha_2/M_1\) across radii and separation.

---

### Scoped verdicts

**Verdict on the derivation (wall-patch slip algebra).**  
Sound. Independent reconstruction confirms: solution for \(v\); wall traction \(f_\parallel=\gamma(v-V_\parallel)=\gamma(T-KV_\parallel)/(K+\gamma)\); equal-and-opposite fluid/wall tractions; power sum \(-2R\); slope formulas; and the precise condition \(T_U+T_V=K\) for opposite wall slopes. No algebraic defect found inside that scope. The derivation does not, and on the given equations cannot, deliver the force on the geometric center.

**Verdict on the record (0PN requirements and release chain).**  
Physically cautious and correctly open. It properly separates wall-patch identities from mouth-center response, rejects boost and carried-momentum arguments as closures for \(A+B=0\), flags E5/E3 as insufficient for aperture side force, and lays out the release chain with each link’s status. The signed conditional release rule is consistent with the stated monopole and linear-response assumptions. The defect of the deliverable is incompleteness of physical inputs (\(D(a)\), \(Q_i\), mouth force law or \(A,B\), \(H_i\), \(\chi_{ij}\), rim closure)—which the record itself declares—rather than a false claim that those quantities are known. No justified numerical two-mouth acceleration or Newtonian \(G\) comparison is available from the supplied material.

**Overall.**  
Pass the wall-patch mathematics as identities of option E5 on the wall. Fail any inference that those identities, the guarded checks, or the archived brane-side loads determine mouth \(A\), \(B\), net center force, or release accelerations. Proceed to a free two-mouth run only after supplying the independent inputs listed above and re-validating with preparations that read total center force in both the companion-incident and isolated-translation tests.
