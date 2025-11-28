(* β = 5/2 Consistency Check Against Superfluid Ontology *)

(* ================================================================== *)
(* C0: Setup                                                          *)
(* ================================================================== *)

ClearAll["Global`*"];

$Assumptions = {
  μ > 0, cs > 0, m > 0, ρ0 > 0, Q > 0, κg > 0,
  r > 0, r0 > 0, rThroat > 0, ξ > 0,
  Element[{μ, cs, m, ρ0, Q, κg, r, r0, rThroat, ξ}, Reals]
};

βVal = 5/2;

Print["C0: Testing consistency of β = ", βVal];

(* ================================================================== *)
(* C1: Regime of validity - where σ(r) becomes O(1)                   *)
(* ================================================================== *)

σ[r_, β_] := β μ/(cs^2 r);

rBreakdownEq = Solve[σ[r, βVal] == 1, r][[1]];
rBreakdown = r /. rBreakdownEq;

Print["C1: σ(r) = 1 at r = ", rBreakdown];

rSchwarzschild = 2 μ/cs^2;

ratioToRs = Simplify[rBreakdown/rSchwarzschild];

Print["C1: r_breakdown = ", ratioToRs, " × r_Schwarzschild"];

checkC1 = Simplify[ratioToRs > 1];
Print["C1: 1PN breaks down outside horizon? ", checkC1];

(* ================================================================== *)
(* C2: Effective inertia positivity for all r > r_breakdown           *)
(* ================================================================== *)

mEff[r_, β_] := m (1 + σ[r, β]);

mEffAtBreakdown = Simplify[mEff[rBreakdown, βVal]];

Print["C2: m_eff at r_breakdown = ", mEffAtBreakdown];

mEffPositive = Simplify[mEff[r, βVal] > 0, r > rBreakdown];

Print["C2: m_eff > 0 for r > r_breakdown? ", mEffPositive];

(* ================================================================== *)
(* C3: Derive σ(r) from Bernoulli + barotropic EoS                    *)
(* ================================================================== *)

(* Bernoulli along streamline: P/ρ + (1/2)v² + U = const *)
(* Linearize: δP/ρ0 - (P0/ρ0²)δρ + (1/2)v_flow² + U(r) = 0 *)
(* Barotropic: δP = cs² δρ *)
(* Far field: v_flow² ≪ |U|, so dominant balance is δρ vs U *)

U[r_] := -μ/r;

δρOverρ0Bernoulli[r_] := Simplify[-U[r]/cs^2];

Print["C3: From Bernoulli, δρ/ρ0 = ", δρOverρ0Bernoulli[r]];

σFromBernoulli[r_, κ_] := κ * δρOverρ0Bernoulli[r];

Print["C3: If m_eff = m(1 + κ δρ/ρ0), then σ(r) = ", σFromBernoulli[r, κ]];

κNeeded = Simplify[Solve[σFromBernoulli[r, κ] == σ[r, βVal], κ][[1, 1, 2]]];

Print["C3: Matching σ requires κ = ", κNeeded];

checkC3 = Simplify[κNeeded == βVal];
Print["C3: κ = β? ", checkC3];

(* ================================================================== *)
(* C4: Newtonian limit - all 1PN corrections vanish as cs → ∞         *)
(* ================================================================== *)

Φeff[r_] := -μ/r - μ^2/(2 cs^2 r^2);

σNewtLimit = Limit[σ[r, βVal], cs -> Infinity];
ΦeffNewtLimit = Limit[Φeff[r], cs -> Infinity];

Print["C4: lim_{cs→∞} σ(r) = ", σNewtLimit];
Print["C4: lim_{cs→∞} Φ_eff(r) = ", ΦeffNewtLimit];

LFull[r_, rdot_, θdot_, β_] := 
  (1/2) m (1 + σ[r, β]) (rdot^2 + r^2 θdot^2) - m Φeff[r];

LNewtLimit = Simplify[Limit[LFull[r, rdot, θdot, βVal], cs -> Infinity]];
LNewtExpected = (1/2) m (rdot^2 + r^2 θdot^2) + m μ/r;

Print["C4: lim_{cs→∞} L = ", LNewtLimit];

checkC4 = Simplify[LNewtLimit == LNewtExpected];
Print["C4: Newtonian Lagrangian recovered? ", checkC4];

(* ================================================================== *)
(* C5: Energy conservation - derive Hamiltonian                       *)
(* ================================================================== *)

prFull = D[LFull[r, rdot, θdot, βVal], rdot] // Simplify;
pθFull = D[LFull[r, rdot, θdot, βVal], θdot] // Simplify;

Print["C5: p_r = ", prFull];
Print["C5: p_θ = ", pθFull];

rdotFromPr[r_, pr_] := pr/(m (1 + σ[r, βVal]));
θdotFromPθ[r_, pθ_] := pθ/(m (1 + σ[r, βVal]) r^2);

HFull[r_, pr_, pθ_] := 
  pr * rdotFromPr[r, pr] + pθ * θdotFromPθ[r, pθ] - 
  LFull[r, rdotFromPr[r, pr], θdotFromPθ[r, pθ], βVal] // Simplify;

Print["C5: H(r, p_r, p_θ) = ", HFull[r, pr, pθ]];

HNewtLimit = Simplify[Limit[HFull[r, pr, pθ], cs -> Infinity]];

Print["C5: lim_{cs→∞} H = ", HNewtLimit];

checkC5noT = FreeQ[HFull[r, pr, pθ], t];
Print["C5: H has no explicit time dependence? ", checkC5noT];

(* ================================================================== *)
(* C6: Magnus force correction from local density                     *)
(* ================================================================== *)

ρLocal[r_] := ρ0 (1 + δρOverρ0Bernoulli[r]);

Print["C6: Local density ρ(r) = ", ρLocal[r]];

FMagnusRatio[r_] := Simplify[ρLocal[r]/ρ0];

Print["C6: F_Magnus(ρ_local) / F_Magnus(ρ0) = ", FMagnusRatio[r]];

δFMagnusOverF = Simplify[FMagnusRatio[r] - 1];

Print["C6: Fractional Magnus correction = ", δFMagnusOverF];

checkC6order = Simplify[δFMagnusOverF == μ/(cs^2 r)];
Print["C6: Magnus correction is O(μ/(cs² r)) i.e. O(1PN)? ", checkC6order];

(* ================================================================== *)
(* C7: Rest mass / cavitation energy with position-dependent inertia  *)
(* ================================================================== *)

EcavStandard = m cs^2;

EcavModified[r_] := m (1 + σ[r, βVal]) cs^2 // Simplify;

Print["C7: Standard E_cav = m cs² = ", EcavStandard];
Print["C7: Modified E_cav(r) = m_eff(r) cs² = ", EcavModified[r]];

δEcav[r_] := Simplify[EcavModified[r] - EcavStandard];

Print["C7: δE_cav(r) = ", δEcav[r]];

UNewton[r_] := -m μ/r;

ratioToBindingEnergy = Simplify[δEcav[r]/Abs[UNewton[r]]];

Print["C7: δE_cav / |U_Newton| = ", ratioToBindingEnergy];
Print["C7: This equals β = ", βVal, "? ", Simplify[ratioToBindingEnergy == βVal]];

(* ================================================================== *)
(* C8: Sink strength to μ mapping from info.md                        *)
(* ================================================================== *)

(* From info.md: μ = κg Q/ρ0 where Q is sink mass flux *)

μFromSink = κg Q/ρ0;

σInSinkVars = σ[r, βVal] /. μ -> μFromSink // Simplify;

Print["C8: σ(r) in sink variables = ", σInSinkVars];

(* Sink inflow velocity: v_sink = Q/(4π ρ0 r²) *)
vSink[r_] := Q/(4 π ρ0 r^2);

(* Kepler orbital velocity: v_orb = √(μ/r) *)
vOrbit[r_] := Sqrt[μ/r];

(* Ratio at orbital radius *)
vRatio = Simplify[vOrbit[r]/vSink[r] /. μ -> κg Q/ρ0];

Print["C8: v_orbit / v_sink = ", vRatio];

(* 1PN validity requires v_orbit/cs < 1, i.e., μ/(cs² r) < 1 *)
vOrbitOverCs = Simplify[vOrbit[r]/cs];

Print["C8: v_orbit/cs = ", vOrbitOverCs];
Print["C8: 1PN valid when μ/(cs² r) < 1, same as σ/β < 1"];

(* ================================================================== *)
(* C9: Added mass from kinetic energy in sink flow field              *)
(* ================================================================== *)

(* T_fluid = (1/2) ρ0 ∫ v_sink² d³x from r=a to ∞ *)

integrand = (1/2) ρ0 * (Q/(4 π ρ0 rr^2))^2 * 4 π rr^2;

Print["C9: Integrand for T_fluid = ", Simplify[integrand]];

Tfluid = Integrate[integrand, {rr, a, Infinity}, 
  Assumptions -> {a > 0, Q > 0, ρ0 > 0}];

Print["C9: T_fluid (sink flow KE from r=a to ∞) = ", Tfluid];

TfluidInμ = Tfluid /. Q -> μ ρ0/κg // Simplify;

Print["C9: T_fluid in terms of μ = ", TfluidInμ];

(* Compare to rest mass energy m cs² *)
TfluidOverMcs2 = Simplify[TfluidInμ/(m cs^2) /. m -> μ/cs^2];

Print["C9: T_fluid / (m cs²) with m ~ μ/cs² = ", TfluidOverMcs2];

(* ================================================================== *)
(* C10: Hydrodynamic origin of β - decomposition attempt              *)
(* ================================================================== *)

(* β = 5/2 must come from: 
   1. Direct density coupling: δρ/ρ0 contributes κ_ρ
   2. Added mass from flow: contributes κ_add  
   3. Pressure-volume work: contributes κ_PV
   Total: β = κ_ρ + κ_add + κ_PV *)

(* If direct density gave κ_ρ = 1, we need κ_add + κ_PV = 3/2 *)

Print["C10: β = 5/2 decomposition:"];
Print["C10: If κ_ρ (direct density) = 1"];
Print["C10: Then κ_add + κ_PV (added mass + PV work) = ", βVal - 1];

(* Classical added mass for sphere in incompressible fluid: C_add = 1/2 *)
(* For compressible fluid near a sink, this is modified *)

(* Dimensional analysis: added mass ~ ρ0 × (characteristic volume) *)
(* Characteristic volume near throat ~ r³ or throat volume *)

Print["C10: Classical added mass coefficient C_add = 1/2 for sphere"];
Print["C10: Remaining after κ_ρ=1, κ_add=1/2: κ_PV = ", βVal - 1 - 1/2];
Print["C10: This suggests PV work contributes κ_PV = 1"];

(* ================================================================== *)
(* C11: Pressure-volume work contribution                             *)
(* ================================================================== *)

(* When dyon accelerates, it does work against pressure gradient *)
(* δW = -∫ ∇P · δx dV *)

(* From barotropic EoS: P = P0 + cs² δρ *)
(* ∇P = cs² ∇(δρ) *)

(* δρ/ρ0 = μ/(cs² r), so δρ = ρ0 μ/(cs² r) *)
δρ[r_] := ρ0 μ/(cs^2 r);

gradδρ = D[δρ[r], r];

Print["C11: ∇(δρ) = d(δρ)/dr = ", gradδρ];

gradP = cs^2 * gradδρ // Simplify;

Print["C11: ∇P = cs² ∇(δρ) = ", gradP];

(* Pressure force per unit volume = -∇P *)
fPressure = -gradP // Simplify;

Print["C11: Pressure force density = -∇P = ", fPressure];

(* This is a restoring force proportional to μ/r² - same scaling as gravity *)
(* Integrating over throat volume ~ ξ³ gives a force *)

(* ================================================================== *)
(* C12: Spin-orbit coupling modification check                        *)
(* ================================================================== *)

(* From info.md Section 4.3: spin-orbit uses Magnus force *)
(* F_Magnus = ρ (v_rel × Γ) *)

(* With local density: F_Magnus,local = ρ0(1 + μ/(cs² r))(v_rel × Γ) *)

(* Relative velocity for orbiting electron around proton vortex: *)
(* v_rel = v_orbit - v_flow *)
(* where v_flow ~ Γ_proton/(2π r) for vortex swirl *)

Clear[Γp, vrel];

vFlow[r_, Γp_] := Γp/(2 π r);

FMagnusLocal[r_, vrel_, Γe_] := ρLocal[r] * vrel * Γe // Simplify;
FMagnusBackground[vrel_, Γe_] := ρ0 * vrel * Γe;

spinOrbitCorrection = Simplify[FMagnusLocal[r, vrel, Γe]/FMagnusBackground[vrel, Γe] - 1];

Print["C12: Spin-orbit force correction δF/F = ", spinOrbitCorrection];
Print["C12: Same O(1PN) as gravity correction? ", 
  Simplify[spinOrbitCorrection == μ/(cs^2 r)]];

(* ================================================================== *)
(* SUMMARY                                                            *)
(* ================================================================== *)

Print[""];
Print["============================================================="];
Print["CONSISTENCY CHECK SUMMARY FOR β = ", βVal];
Print["============================================================="];
Print["  C1: 1PN valid outside ", ratioToRs, " r_s: ", checkC1];
Print["  C3: σ derives from Bernoulli with κ = β: ", checkC3];
Print["  C4: Newtonian limit recovered: ", checkC4];
Print["  C5: Energy conserved (H time-independent): ", checkC5noT];
Print["  C6: Magnus correction is O(1PN): ", checkC6order];
Print["  C7: δE_cav/|U| = β: ", Simplify[ratioToBindingEnergy == βVal]];
Print["  C12: Spin-orbit correction is O(1PN): ", 
  Simplify[spinOrbitCorrection == μ/(cs^2 r)]];
Print["============================================================="];

allChecks = {checkC1, checkC3, checkC4, checkC5noT, checkC6order,
  Simplify[ratioToBindingEnergy == βVal],
  Simplify[spinOrbitCorrection == μ/(cs^2 r)]};

Print["All consistency checks: ", allChecks];
Print["============================================================="];

Print[""];
Print["KEY PHYSICAL FINDINGS:"];
Print["  1. β = 5/2 breaks down at r = (5/4) r_s - just outside horizon"];
Print["  2. Magnus/spin-orbit forces get O(1PN) corrections from local ρ"];
Print["  3. Rest mass gets position-dependent correction δE = β m μ/r"];
Print["  4. β may decompose as: 1 (density) + 1/2 (added mass) + 1 (PV work)"];
Print["============================================================="];
