# Directive pathA_21b — Close the inter-defect force HONESTLY + write the CODEABLE stationary-throat BVP (the option-C spec)

**Status:** Design-reviewed (Codex `gpt-5.5` → SOUND-WITH-FIXES; all 5 must-fix + 3 should-fix applied 2026-06-19);
pending Codex confirm-pass → user execution gate. This is a
corrective + spec-hardening step that follows the pathA_21 execution + 5-agent review (decision-13 §11). pathA_21's
NEGATIVE verdicts (P2 `MASS_BRIDGE_FORM_NOT_DERIVED`, `EP_NOT_DERIVED`, P3 `{L,T,M}` retained, P4
`NEWTON_G_FORM_NOT_DERIVED`) were ALL adversarially re-confirmed HONEST and are CARRIED FORWARD UNCHANGED — do NOT re-open
them. pathA_21b fixes the two things the review found short:
1. **P1 was a RESTATEMENT, not a derivation.** The `1/r²` drain field was hand-inserted as a literal
   (`velocity_from_1 = -q1/(4π r²)`), never solved from continuity; the "derivation" check reduced to `x==x`; and the
   attractive SIGN came from declaring symbols `positive=True` + choosing `Q_i>0` (a convention the directive forbids),
   not from the pressure/compressibility response. The verdict label `..._DERIVED_CONDITIONAL` OVERCLAIMS.
2. **The P5 profile-solve spec is honest but NOT computable.** It is a residual ledger, not a closed boundary-value
   problem: no stationary field-equation set (no time-independent GNLS, no stationary Maxwell, no equation that SELECTS
   `R0(w)`), the `I_F,12` cross-stress integrand `Π_cross` is never written out, and `J`/`W_eff`/the reduction kernel
   have no closing conditions. An option-C solver engineer would have to invent equations.

**Date:** 2026-06-19
**Owner:** Codex (DERIVES + codes; iterates until scripts exit 0). Claude reviews afterward.
**Trigger:** pathA_21 review (decision-13 §11). Chain: `pathA_19`→`pathA_20`→`pathA_20b`→`pathA_21`→**this
(`pathA_21b`: honest force + codeable BVP)**→ **option C** (the throat-profile SOLVE, now actually specified) →
`pathA_22` (scale-map → `m̂0²·S_port` → re-run B2c).

## Why this step
The convergence finding (pathA_20 §9, pathA_21 §11) is real and now sharpened: EVERY remaining number bottlenecks on the
SOLVED stationary throat profile. pathA_21 was supposed to (a) derive the force FORM and (b) produce a spec concrete
enough to drive that solve. It did neither rigorously: (a) it assumed the force, and (b) the spec names outputs but gives
no closed problem. pathA_21b makes both real: it closes the transcribable sector (the stationary field equations) so
option C can be coded for it, and it replaces every remaining gap with a NAMED branch-realization residual that option C
must choose or solve — without faking any missing physics (the genuine new-physics walls stay residuals).

## Scope & stance
SYMBOLIC DERIVATION + SPECIFICATION (forms + a closed BVP statement; NOT the numerical solve — that is option C). Carry
`λ_γ`, `α_J`, the `J`-value, and the genuinely-unsolved profile functionals SYMBOLIC; do NOT manufacture missing
profile/operator data. No model-formula changes, no freeze touch, no `m̂0²·S_port` un-pin (`pathA_22`). Mass-symbol
discipline: `m_GNLS` (constituent: EOS/healing/Madelung/action) vs `m_defect` (throat). Extend `dimensional_check.py`
side-by-side (NEW group `--patha21b-*`; leave `dimensional_dictionary()`, global `D`, D1–D3, and the
pathA_18/19/20/20b/21 groups intact — additive only). Same infra constraints: read-only `research/pde_audit/simulation/`;
never touch `physical_export_permitted`; `research/pde_audit/{scripts,notes}` READ-ONLY; no `$RT exec`; standalone
`python3` + `math -script`; each script under `timeout 600` (exit 124 = failure → reformulate); ≤2 concurrent MMA seats;
YAML/markdown for human/LLM output, JSON only for machine artifacts.

**Carried-forward UNCHANGED (do NOT re-derive or weaken):** P2 `MASS_BRIDGE_FORM_NOT_DERIVED` + `EP_NOT_DERIVED`; P3
`{L,T,M}` retained (`HBAR_FREE_SUBSTRATE_RELATION_MISSING`); P4 `NEWTON_G_FORM_NOT_DERIVED`; all pathA_20/20b residuals
(`STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`, `H_2PI_RATE_CLASSIFICATION_UNDETERMINED`,
`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`, `HBAR_PROVENANCE_UNDETERMINED`). These were re-confirmed honest by three
independent adversarial agents; pathA_21b neither re-opens nor "improves" them.

**DERIVED-FORM GATE (anti-restatement — still the core teeth).** A positive force/field form is accepted ONLY with a
SOURCE-EQUATION CHAIN from the parent equations to the result, with every intermediate step shown (no hand-inserted
target field, no `x==x` self-checks standing in for a derivation). A sign asserted via `positive=True`/symbol-positivity
declarations is NOT a derived sign. Dual-engine `.wl`/SymPy agreement validates DIMENSIONS/ALGEBRA ONLY; it cannot upgrade
a restatement into a derivation. A missing chain → a NAMED RESIDUAL (a VALID outcome), never a PASS or a "conditional
derivation."

## Work items

### P1b — DERIVE the inter-defect force from the stationary drain field (no hand-inserted `1/r²`, no convention sign)
- **Solve, do not assume, the far-field drain velocity.** From stationary continuity `∇·(ρ0 𝐯)=0` with a localized
  drain/source of strength fixed by the flux `J`, SOLVE for the far-field radial velocity in BOTH the 4D-bulk geometry and
  the reduced-3D geometry. The `r`-power must EMERGE from the dimensionality of the divergence operator (Gauss on the
  enclosing sphere), not be inserted. Show the explicit step `∮ ρ0 𝐯·d𝐒 = (source) ⇒ v_r(r)`. If the background `ρ0(r)`
  is not asymptotically constant, carry the resulting correction symbolically.
- **Compute the force via the full momentum-flux stress tensor.** Write the cross-term integrand `Π_cross` EXPLICITLY —
  which terms enter (convective `m_GNLS ρ 𝐯⊗𝐯`, pressure `h(ρ)`, quantum stress from `Q(ρ)`, `V_conf`) and on what
  control surface `∂U₂` around defect 2 — then derive the force `F_12` by integrating `Π_cross·𝐧` over `∂U₂`. The force
  coefficient `C_F` and the `r`-power are RESULTS of this integral, not posited.
- **Derive the SIGN (attractiveness) from the pressure/compressibility response.** Show attractiveness follows from the
  Bernoulli pressure drop in the entrained flow (`sign(dh/dρ)` / the density-response sign), NOT from declaring symbols
  positive. Any EOS/stability sign used (e.g. `dh/dρ>0`) must be CITED from the parent EOS (`P=Kρ⁵ ⇒ h=(5K/4)ρ⁴`) or a
  stated stability condition — NOT introduced as a CAS `positive=True` assumption; the defect coupling ORIENTATION (which
  way each drain pulls) remains a profile-residual if not derived. If the sign genuinely cannot be fixed without the
  solved interior profile, that is acceptable ONLY as the explicit residual `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL` (name
  the exact missing profile piece) — never via a positivity convention.
- **Honest outcome menu:** either (i) a genuinely-derived far-field force law `F_12(C_F, r-power)` with a derived sign and
  the remaining interior factors flagged to P5b, OR (ii) a NAMED RESIDUAL with the specific missing equation. Either way,
  REPLACE the overclaiming pathA_21 verdict label with the honest one; mark
  `P1_INVERSE_SQUARE_FIELD_ASSUMED_NOT_SOLVED` RESOLVED only if the continuity/Gauss solve actually succeeds, otherwise
  `SUPERSEDED_BY_<named residual>`.

### P5b — the stationary-throat BVP: closed where transcribable, NAMED RESIDUALS where branch-realized
Produce the boundary-value problem the option-C solve will implement. For EACH equation/condition: the PDE/condition
itself, its domain, the fields it acts on, every boundary condition, the measure, the dimension, the frame
(4D-bulk/brane/reduced-3D), and a real source anchor (file:line or eq label). **Per-item bar:** for each of G1–G5, EITHER
(i) transcribe a SOURCE-ANCHORED equation/BC/closure from the parent action, OR (ii) emit a NAMED branch-realization
residual (and, where useful, list candidate branch choices as explicit ASSUMPTIONS, not closure). A residual for branch
data the parent does NOT select is a VALID outcome, NOT a failure — faking closure to avoid one is the failure. The win is
that the transcribable sector (notably G1) is genuinely closed and every remaining unknown is a named, enumerated
branch-realization residual rather than a silent gap. The six review gaps:
- **G1 — stationary field-equation set (EXPECTED CLOSEABLE).** Write the time-independent GNLS for `ψ0(X)` (the stationary
  limit of the parent gauged GNLS: chemical potential / eigenvalue, the `−ħ²/2m_GNLS ∇₄²`, `V_conf`, `h(ρ0)`, gauge
  coupling), with its asymptotic + throat BCs; and the stationary Maxwell equation for `A_{M0}` with its gauge condition +
  BCs. These ARE transcribed from the parent action's stationary limit — derive them and anchor them; this is the closed
  core of the BVP.
- **G2 — the condition on `R0(w)`.** ATTEMPT to derive a free-boundary / Euler–Lagrange condition by varying the action
  w.r.t. the surface `Σ0 = r − R0(w)`. If the parent provides such a selecting equation, state it (no "topology BC"
  placeholder); if it does NOT, record `R0_FREE_BOUNDARY_CONDITION_UNDERIVED` and list the candidate bottom-topology BCs
  as explicit ASSUMPTIONS (enumerated branches), not closure.
- **G3 — the `Π_cross` momentum-flux integrand + control surface** (shared with P1b). DERIVE `Π_cross` from the action /
  Noether momentum stress (or the Euler stress balance), written out term-by-term on a specified control surface `∂U₂`,
  giving the computable functional `I_F,12`. If it is not source-anchorable, record `PI_CROSS_STRESS_TENSOR_UNDERIVED`
  and mark `I_F,12` a residual (status `new-physics`/branch), NOT `profile-solve`.
- **G4 — the `J`-value.** State the conservation / no-leakage BCs, then identify any SEPARATE source-anchored
  regularity / choking / energy condition that fixes the VALUE of `J` (no-leakage alone conserves flux but need not fix
  `J`). If such a selector exists, state it; if none does, record `J_VALUE_BRANCH_PARAMETER` / `J_SELECTOR_UNDERIVED` and
  say what would fix it downstream. No silent free parameter.
- **G5 — the reduction kernel.** Give the explicit projection FORMULAS for the declared transverse weights `W(w)`,
  `χ_N(w)` and the integrals that yield `W_eff` and `N_∞,3` (the parent DECLARES/controls these weights). If the kernel
  SHAPE is not source-SELECTED, keep `W_KERNEL_UNDERDECLARED` / `W_EFF_REDUCTION_UNDERIVED`, or label a candidate kernel
  as an explicit option-C branch assumption — not a placeholder definition.
- **G6 — brane `c_γ/c_s` (`λ_γ`).** Needs the brane zero-mode reduction, which pathA_20b proved underived. KEEP the
  explicit residual `BRANE_ZERO_MODE_REDUCTION_UNDERIVED` — do NOT fake closure. State precisely what the option-C solve
  would have to compute to close it.
- **Genuine new-physics stays residual.** `α_J` as an independent profile functional, the `ħ`-provenance, and the
  cycle-vs-angular `2π` classification are NOT in scope to close here (Agent D re-confirmed these walls). Keep them as
  named residuals; do NOT introduce a restatement (e.g. `α_J := m_defect c²/(ħJ)`, or `α_H,ω := H_throat/(ħJ_ω)` accepted
  as the bridge) to make the BVP look complete.
- **Anchor-fix:** correct the `α_H,ω` source anchor — `pde.tex:318-391` is the parent ACTION, not a canonical
  Hamiltonian; re-label as action-level and flag "canonical Hamiltonian must be constructed" as the residual.

### P5b output form
Emit the BVP as (a) a human-readable spec doc (`reports/` or `notes/`) with the equations + BC table, AND (b) an updated
machine table that supersedes the pathA_21 P5 table, with the SAME schema columns
(`symbol`/`definition`/`dimension`/`frame`/`source anchor`/`closes-which-output`/`status`/`residual-if-absent`/
`downstream consumer`) but where every `profile-solve` row whose closing equation now exists points to that equation
(status `profile-solve` with a CONCRETE BVP reference), and only the true new-physics rows remain
`new-physics`/placeholder-free residuals. Include the `𝔅` items (`pde.tex:2515-2566`); keep the
pathA_21b-needed vs pathA_22-deferred split.

## Acceptance criteria (PASS/FAIL with NAMED RESIDUALS; exit-0 NECESSARY not SUFFICIENT)
1. **P1b force:** the far-field drain velocity is SOLVED from continuity (the `r`-power emerges from Gauss, not inserted);
   `F_12` follows from an EXPLICIT `Π_cross` surface integral; the sign is DERIVED from the pressure/compressibility
   response OR carried as `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL`. NO hand-inserted `1/r²`, NO `x==x` check standing in
   for a derivation, NO positivity-convention sign. The pathA_21 P1 verdict label is corrected.
2. **P5b BVP:** for each of G1–G5, EITHER a source-anchored equation/condition (PDE + domain + BCs + measure + anchor) OR
   a NAMED branch-realization residual (candidate assumptions enumerated where useful); G1 (the stationary field
   equations) is expected CLOSED; G6 + the new-physics items are explicit residuals. NO P5b row has an UNNAMED placeholder
   definition (a named residual is NOT a placeholder). The `α_H,ω` anchor is corrected.
3. **Carried negatives unchanged:** P2/EP/P3/P4 verdicts and the pathA_20/20b residuals appear verbatim, neither re-opened
   nor weakened.
4. **Harness:** new `--patha21b-*` group passes; dual-engine `.wl` agreement for ALGEBRAIC/DIMENSIONAL claims only (NOT a
   proof of any derivation, and NO `expr===expr` self-comparison passed off as a check); scripts exit 0 within
   `timeout 600`; pathA_18/19/20/20b/21 groups + `dimensional_dictionary()` untouched (additive only).

**Fail conditions (explicit):** a hand-inserted drain field or `r`-power; a sign from a positivity convention; an `x==x` /
`expr===expr` check presented as verifying a non-trivial relation; any G1–G5 left as an UNNAMED placeholder (a named
branch-realization residual is a VALID outcome, NOT a failure); faking a G1–G5 closure the parent does not actually
provide; faking G6 / `α_J` / the `2π` classification with a restatement; re-opening or weakening a carried negative;
mutating a prior group or the dictionary. VALID expected outcomes: a derived far-field force with interior factors
flagged; an honest `ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL`; G2/G4/G5/G6/`α_J`/`ħ` remaining named residuals; a BVP that is
closed for the transcribable sector (G1) and explicitly, enumerably residual for the branch-realization + new-physics
sectors.

## Out of scope
The numerical profile SOLVE itself (option C — pathA_21b only SPECIFIES it as a closed BVP); the scale-map →
`m̂0²·S_port` → B2c rerun (`pathA_22`); any freeze amendment; closing `α_J`/`ħ`/brane-`c_γ` (genuine new-physics walls).

## Review (orchestrator, after Codex)
Transliteration-fidelity (one clean agent per new script + the `.wl`); an adversarial pass with the SAME
distrust-restated-target teeth that caught pathA_21 — specifically: is the drain velocity actually SOLVED from continuity
(show the Gauss step) or re-inserted? does the `r`-power EMERGE? is `Π_cross` written out and the force a real surface
integral? is the sign from the pressure response or a convention? are G1–G5 genuinely codeable (would a solver engineer
guess anywhere)?; plus a completeness-critic pass on the BVP (what's still missing to run option C?). Verify the carried
negatives are verbatim. Claude reads only residuals. Then gate to option C.
