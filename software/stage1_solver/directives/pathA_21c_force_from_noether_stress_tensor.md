# Directive pathA_21c — Inter-defect force + SIGN from the GNLS momentum-BALANCE tensor (calibrate-predict framed)

**Status:** Design-reviewed (Codex `gpt-5.5` → SOUND-WITH-FIXES; all 5 must-fix + 4 should-fix applied 2026-06-19, plus
the calibrate-predict reframe per the user); pending Codex confirm-pass → execution (user already chose pathA_21c). Follows
the pathA_21b 4-agent review, which confirmed pathA_21b's drain-velocity Gauss solve + the closed G1 BVP but caught that
the inter-defect force **COEFFICIENT** was still a heuristic product `F = m_GNLS·N_∞,3·Q2·v1` dressed as a Π_cross surface
integral (the tensor written as prose, never integrated; the parent supplies only an Euler force-PER-VOLUME identity).

**Date:** 2026-06-19
**Owner:** Codex (DERIVES + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** pathA_21b review (decision-13 §11/§12); user chose "pathA_21c: derive force + sign" (2026-06-19). Chain:
`…pathA_21`→`pathA_21b`→**this**→ option C (throat-profile SOLVE) → `pathA_22`.

## Why this step
A force between two sources is correctly obtained by integrating the conserved/balanced **momentum-flux (stress) tensor**
over a control surface enclosing one source (the fluid analog of the Maxwell stress tensor for charges) — PLUS any
body-force volume term where the medium is not translation-invariant. pathA_21b skipped the tensor and asserted the
coefficient. Doing it properly (a) closes the overclaim honestly, (b) confirms/corrects the leading far-field coefficient
and `r`-power as RESULTS, and (c) **attempts to derive the SIGN** (attractive vs repulsive) from the integrated flux —
a genuinely target-blind result (we never tuned "attractive"). The momentum-flux tensor is a STANDARD Noether object; this
is a derivation, not new physics.

## Calibrate-predict framing (the success criterion — read first)
This is a TOY ANALOG validated by **calibrate-on-a-subset-then-predict-the-surplus** ([[feedback-calibrate-predict-methodology]],
[[project-analog-framework-goal]]). Sort every output of this step into exactly one bucket and LABEL it:
- **Target-blind PREDICTIONS (derive these; do not tune):** the force STRUCTURE (bilinear in the two sources,
  `∝ Q1 Q2 / r²` in the compact reduced-3D lane), the `r`-power, and the **SIGN** (attractive/repulsive). These were never
  calibrated to anything, so deriving them is real surplus.
- **CALIBRATION KNOB (allowed; count it, do NOT call it "derived"):** the force's overall NORMALIZATION (the dimensionless
  `I_F,12^full` / `Θ_Q` / branch-data piece), which a downstream step fixes by matching a trusted anchor (Newtonian `G` /
  the GR quadrupole `54/5`). Leaving this un-derived ab initio is NOT a flaw — it is the bootstrap anchor.
The DERIVED-FORM GATE forbids ONE thing: **calling a calibrated/posited quantity "derived."** It does NOT forbid
calibration. The win here is the derived structure + sign (or honest `SIGN_RESIDUAL`) with the normalization labeled a
knob — feeding downstream surplus (`G`, then g−2 / 5PN / multi-defect). Forbidden: (a) claim the SAME observable as both
calibrated-to and predicted; (b) knobs ≥ independent predictions. Count degrees of freedom honestly.

## Scope & stance
ANALYTIC far-field derivation (forms + a definite sign where the far field allows it; otherwise a precisely-named
residual). It does NOT need the solved interior profile: it uses the pathA_21b Gauss-solved far-field drain fields
`v1, v2` and the parent action. Carry the interior profile correction as a named residual; do NOT manufacture interior
data. Mass-symbol discipline: `m_GNLS` (constituent) vs `m_defect` (throat). Extend `dimensional_check.py` side-by-side
(NEW group `--patha21c-*`; leave `dimensional_dictionary()`, global `D`, D1–D3, and the pathA_18/19/20/20b/21/21b groups
intact — additive only). Same infra: read-only `research/pde_audit/simulation/`; never touch `physical_export_permitted`;
`research/pde_audit/{scripts,notes}` READ-ONLY; no `$RT exec`; standalone `python3` + `math -script`; each script under
`timeout 600` (exit 124 = failure → reformulate); ≤2 concurrent MMA seats; YAML/markdown human output, JSON only for
machine artifacts.

**The medium is NOT globally translation-invariant.** The parent action carries explicit backgrounds — `V_conf(X;Σ0)`,
the gauge localization `Z(w)`, and the external source `J_ext` — which break exact spatial translation invariance. So the
Noether construction yields a momentum-**BALANCE** law, NOT a pure conservation law:
`∂_t g_i + ∂_j Π_ij = f_i^body`, where `Π_ij` collects the translation-invariant field stresses and `f_i^body` collects
the explicit-background body forces. Every explicit-background term must be either carried in `f_i^body` (and accounted
for in the force) or carried as a NAMED residual — none may be silently dropped.

**Carried-forward UNCHANGED (do NOT re-derive or weaken):** the pathA_21b G1 stationary BVP + its BC table; the
drain-velocity Gauss solve (`v_r=−Θ_Q J/(4π N_∞,3 r²)` reduced-3D, `v_R∝R⁻³` bulk with `Ω_3=2π²`); the honest residuals
G2/G4/G5/G6 (do NOT close `W_eff`/G5); the pathA_21b `I_F,12^full` definition (this step may add `ΔI_F,12^profile` or
state the supersession explicitly, but must not silently weaken the prior residual); P2 `MASS_BRIDGE_FORM_NOT_DERIVED` /
`EP_NOT_DERIVED`, P3 `{L,T,M}`, P4 `NEWTON_G_FORM_NOT_DERIVED`; all pathA_20/20b residuals.

**DERIVED-FORM GATE (anti-restatement teeth).** The force is accepted as derived ONLY when it is
`F_12 = −∮_{∂U₂} Π_ij n̂_j dS + ∫_{U₂} f_i^body dV` with `Π` and `f^body` **derived from the parent action via Noether**
and the integrals actually carried out, so the STRUCTURE, `r`-power, and SIGN come OUT. Explicit FAILs: (1) asserting
`F = m_GNLS·N·Q2·v1` (or any heuristic momentum-uptake product) and labeling it a surface integral; (2) positing `Π` or
`f^body` instead of deriving them; (3) a sign from a positivity convention (`positive=True`, choosing `Q_i>0`); (4) a
harness "check" that is `x==x`, `expr===expr`, or `Solve[...]` compared to its own output standing in for a derivation
(the pathA_21b Gauss-velocity checks were this — do NOT repeat); (5) silently dropping an explicit-background term instead
of carrying it in `f^body` or as a residual; (6) choosing a stress-tensor improvement that changes the cross integral or
the sign without proving it. A piece the far field genuinely cannot fix → a PRECISELY NAMED residual (valid), never the
heuristic shortcut and never a calibrated value mislabeled "derived." Dual-engine `.wl`/SymPy validates DIMENSIONS/ALGEBRA
only; the derivation is the Noether construction + the explicit integrals.

## Work items

### P1c-A — derive the momentum-BALANCE tensor `Π_ij` + the body forces `f_i^body` (Noether)
- From the parent gauged-GNLS action (`pde.tex:318-391`), derive via Noether (spatial translation of the
  translation-invariant field terms) the canonical momentum-flux tensor `Π_ij`, written term-by-term: convective
  `m_GNLS ρ v_i v_j`; EOS pressure `δ_ij P(ρ)`, `P=Kρ⁵`; the Madelung **quantum stress** `σ_Q,ij` from the `−ħ²/2m ∇²`
  term; and the gauge/Maxwell stress where active.
- **Stress-representative discipline:** derive the CANONICAL tensor first. If a symmetrization/Belinfante improvement is
  used, state it explicitly and PROVE its closed-surface cross integral is zero or strictly subleading; otherwise carry
  `STRESS_IMPROVEMENT_AMBIGUITY_RESIDUAL`.
- **Explicit-background body forces:** `V_conf(X;Σ0)`, `Z(w)`, and `J_ext` break translation invariance → collect their
  contributions in `f_i^body` (e.g. `f_i^Vconf = −ρ ∂_i V_conf`), NOT inside `Π`. List each; any not put in closed form is
  a named residual (`VCONF_BODY_FORCE_RESIDUAL`, etc.).
- **VERIFY** the balance law `∂_t g_i + ∂_j Π_ij = f_i^body` reproduces the parent Euler force-per-volume identity
  (`pde.tex:445-451`; `part01:280-286`) — dimension + algebra dual-engine checked. This is the anchor that the tensor is
  the RIGHT one.

### P1c-B — the inter-defect force as an explicit integral (structure + `r`-power as RESULTS; normalization = knob)
- **Sign convention (pin it):** `n̂` outward from the control region `U₂`; `F_12` = force ON defect 2 BY defect 1; in the
  stationary case `F_12 = −∮_{∂U₂} Π_ij n̂_j dS + ∫_{U₂} f_i^body dV`, with `∂U₂` a small closed reduced-3D surface
  enclosing defect 2 and excluding defect 1.
- Evaluate `Π_cross` (the part of `Π[v1+v2]` bilinear in the two sources) on the carried Gauss far-field fields. The
  force STRUCTURE (`∝ Q1 Q2 / r²`) and `r`-power must EMERGE — do NOT posit `F = m_GNLS·N·Q2·v1`. The **lane**: the `4π`
  reduced-3D coefficient uses the carried reduced-3D Gauss solve; any bulk-4D result must use `Ω_3=2π²` and stay separate.
- **Availability discipline (derive or residualize each):** the density response `δρ` is NOT given by the Gauss `v` — derive
  it from far-field Bernoulli/EOS (state assumptions) or carry `DENSITY_RESPONSE_RESIDUAL`; the quantum stress
  contribution must be derived and shown subleading/integrable in the far field or carried as
  `QUANTUM_STRESS_FAR_FIELD_RESIDUAL`; the `V_conf` body-force volume term must be computed or carried as
  `VCONF_BODY_FORCE_RESIDUAL`; the Maxwell stress enters the coefficient ONLY after proving the `Z(w)`/gauge-fixing/`J_ext`
  terms vanish or cancel in the selected lane, else residualize.
- **Normalization = CALIBRATION KNOB:** the overall dimensionless prefactor (`I_F,12^full` / `Θ_Q`) is an allowed knob,
  to be fixed downstream against an anchor — label it as such, do NOT call it derived, and do NOT close `W_eff`/G5. Confirm
  or correct the pathA_21b structural form `C_F ∝ m_GNLS N_∞,3 Q1 Q2/(4π)`.

### P1c-C — the SIGN verdict + the calibrate-predict ledger
- From the computed integral (surface flux + body-force volume term, with the pinned convention), determine whether `F_12`
  points toward or away from defect 1. The sign must come from the integrated convective + pressure (+ Maxwell, if
  legitimately included) flux + the body-force term — NOT a convention. Verdict: `FORCE_ATTRACTIVE_DERIVED` /
  `FORCE_REPULSIVE_DERIVED`, OR — if the far field genuinely cannot fix it without a specific interior piece — a
  precisely-named `SIGN_RESIDUAL_<piece>` (e.g. `SIGN_RESIDUAL_PROFILE_TRACTION_ORIENTATION`,
  `SIGN_RESIDUAL_QUANTUM_VCONF_PROFILE`). A derived sign is the headline target-blind result; an honest `SIGN_RESIDUAL` is
  an acceptable landing — do NOT manufacture a sign to "win."
- **Calibrate-predict ledger:** state explicitly which outputs are target-blind PREDICTIONS (structure, `r`-power, sign)
  vs CALIBRATION KNOBS (normalization, branch data), count them, and name the downstream surplus this feeds (`G`, then
  g−2 / 5PN / multi-defect). Replace the pathA_21b P1b verdict label with the honest pathA_21c label reflecting all of
  the above.

## Deliverables
- A new report `reports/pathA_21c_force_from_noether_stress_tensor.md` (the Noether derivation trace, the balance law, the
  force integral, the sign verdict, the calibrate-predict ledger, the carried items).
- A new additive harness group `--patha21c-*` in `dimensional_check.py` (dims + algebra; NO tautological self-checks) + an
  independent `.wl` cross-check.
- Consumed inputs (read first): the pathA_21b report + directive; decision-13 §11/§12; the parent anchors
  (action `pde.tex:318-391`, Euler `:445-451`, EOS `:342-352`, continuity `:396-406`; `part01:275-286`).

## Acceptance criteria (PASS/FAIL with NAMED RESIDUALS; exit-0 NECESSARY not SUFFICIENT)
1. **P1c-A:** `Π_ij` + `f_i^body` DERIVED from the action via Noether (term-by-term), the balance law
   `∂_t g_i + ∂_j Π_ij = f_i^body` shown to reproduce the parent Euler force-per-volume identity (dim + algebra checked);
   stress-improvement ambiguity proven-zero/subleading or residualized; every explicit-background term placed in `f^body`
   or named.
2. **P1c-B:** `F_12 = −∮ Π_ij n̂_j dS + ∫ f_i^body dV` carried out with the pinned sign convention; the STRUCTURE and
   `r`-power are RESULTS; `δρ`/quantum-stress/`V_conf`/Maxwell each derived or precisely residualized; the overall
   normalization labeled a CALIBRATION KNOB (not "derived"); reduced-3D vs bulk lane explicit; `W_eff`/G5 not closed. NO
   heuristic product.
3. **P1c-C:** an explicit SIGN verdict — `FORCE_ATTRACTIVE_DERIVED` / `FORCE_REPULSIVE_DERIVED` from the integral, OR a
   precisely-named `SIGN_RESIDUAL_<piece>`; no positivity-convention sign; the calibrate-predict ledger (predictions vs
   knobs, counted) stated; the pathA_21b P1b label corrected.
4. **Carried items unchanged:** the pathA_21b G1 BVP + Gauss solve + residuals (incl. `I_F,12^full`); P2/EP/P3/P4;
   pathA_20/20b residuals — verbatim, neither re-opened nor weakened.
5. **Harness:** new `--patha21c-*` group passes; dual-engine `.wl` agreement for ALGEBRAIC/DIMENSIONAL claims only, with
   NO `x==x` / `expr===expr` / `Solve`-compared-to-its-own-output posing as a derivation; scripts exit 0 within
   `timeout 600`; prior groups + `dimensional_dictionary()` untouched (additive only).

**Fail conditions (explicit):** the heuristic product `F=m·N·Q2·v1` labeled a surface integral; positing `Π`/`f^body`
instead of deriving them; a sign from a positivity convention; a tautological harness check standing in for a derivation;
silently dropping an explicit-background term; an undeclared stress-improvement that changes the cross integral/sign;
calling the calibrated normalization "derived"; claiming the same observable as both calibrated-to and predicted;
re-opening/weakening a carried item; mutating a prior group or the dictionary. VALID expected outcomes: a derived
far-field structure + `r`-power with the normalization a labeled knob and interior pieces named residuals; a DERIVED
attractive/repulsive SIGN (the win), OR a precisely-named `SIGN_RESIDUAL_<piece>`.

## Out of scope
The throat-profile SOLVE itself (option C — pathA_21c is analytic far-field); the scale-map → `m̂0²·S_port` → B2c rerun
(`pathA_22`); fixing the calibration normalization to an anchor (that is the downstream calibrate step); `α_J` /
`ħ`-provenance / the `2π` classification (carried new-physics walls); any freeze amendment.

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script + the `.wl`); an adversarial distrust-restated-target pass with
the SAME teeth that caught pathA_21 AND pathA_21b: is `Π_ij` + `f^body` actually Noether-DERIVED from the action (show the
variation) or posited? does the balance law reproduce the parent Euler identity? is `F_12` a REAL integral whose
structure/`r`-power/sign come OUT, with the pinned sign convention, or the heuristic product relabeled again? is the SIGN
from the integrated flux or a convention (or honestly residualized)? is the overall normalization correctly labeled a
calibration KNOB (not "derived"), and is the calibrate-predict ledger's knob/prediction count honest? are the harness
checks substantive (no Solve-vs-own-output)? Verify carried items verbatim. Claude reads only residuals. Then gate to
option C.
