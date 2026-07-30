# I-2 (ledger_stage005) source map — pathA_20 sound speed + light/sound ratio

> Prep note for the NEXT stage (I-2), captured 2026-07-07 from the Part-I Explore fan-out so a fresh session can author the
> reshape directive without re-running discovery. Verify against the source reports before finalizing the directive.
> Companion: `part1_medium_atomic_split.md` (the finalized Part I split). Build-order id **005**, Part I.
> **Headline verdict token:** `C_GAMMA_RATIO_UNDERDETERMINED`. **Nature:** EARNED (c_s) + **CALIBRATED (λγ free — the honest landing)**.
> Not a FAIL-headline; the honest landing is "the light/sound ratio is unpinned by the parent action → a free calibration knob."

## File inventory
- **Directives:** `software/stage1_solver/directives/pathA_20_emergent_constants_cs_c.md`,
  `software/stage1_solver/directives/pathA_20b_cgamma_cs_linearization.md`
- **Reports:** `software/stage1_solver/reports/pathA_20_velocity_constants.md`,
  `software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md`
- **Tools (.wl):** `tools/pathA_20_velocity_constants_crosscheck.wl`, `tools/pathA_20b_cgamma_cs_crosscheck.wl`
- **SymPy source-of-truth (extract, do NOT import):** shared harness `software/stage1_solver/src/stage1_solver/dimensional_check.py`
  — `run_patha20_velocity_constants` (CLI `--patha20-velocity`; 21/21 dim + 5/5 alg), `run_patha20b_cgamma_cs`
  (CLI `--patha20b-cgamma-cs`; 11/11 dim + 7/7 alg). **NO standalone per-gate `.py`** → the reshape extracts one (exactly as I-1 did;
  pathA_20/20b are NOT truly `.wl`-only — the SymPy already exists+passes in the harness).
- **Verdict tokens:** pathA_20 = `C_GAMMA_RATIO_UNDERDETERMINED` (headline) + `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA` (flux)
  + `HBAR_PROVENANCE_UNDETERMINED` (ħ). pathA_20b sharpens the ratio → `C_GAMMA_BULK_UNDERDETERMINED` (bulk) +
  `C_GAMMA_RATIO_STILL_UNDERDETERMINED` (brane).

## Atomic claim-set (distilled; the directive's §2)
1. **Sound speed from EOS (EARNED rel. to the imposed EOS).** `c_s²(ρ) = (1/m)·dP/dρ = 5Kρ⁴/m` from `P=Kρ⁵`; `[c_s]=LT⁻¹`;
   state dependence `c_s ∝ ρ²`, `dln c_s/dln ρ = 2`. EOS is POSTULATED (`EOS_CLOSURE_IMPOSED`).
2. **Three velocity scales (dims EARNED, ratios not derived).** `v_b=(ħ/m)∇θ` (condensate flow), `c_s` (phonon), `c_γ` (photon/gauge).
3. **`c=c_γ` terminal ceiling — DERIVED from the wave sector, NOT `E=mc²`.** `ω²=c_γ²k²`; trapped transverse mode
   `ω²=c_γ²(k∥²+k⊥²)`; bound-mode clock `ω₀=c_γ k⊥`, boosted → internal clock `ω₀/γ` (time dilation). Non-circular. EARNED.
4. **λγ = c_γ/c_s is a FREE calibration input (the headline / CALIBRATED landing).** Carried symbolic; tail `(c/c_s)³=λγ³`.
   `C_GAMMA_RATIO_UNDERDETERMINED` — the parent action does NOT pin it; `c_γ=c_s` from shared dims or weak-field prose is REJECTED.
   (Downstream pathA_22a labels λγ a "TRUE FREE calibration knob"; it enters `54/5` to the 5th power.)
5. **Mass-bridge candidate (POSTULATED → pathA_21).** `m_defect=α_J·ħJ/c_γ²`; recorded only, does NOT collapse `M`, `α_J` not derived.
6. **Flux law (UNDETERMINED).** Steady no-leakage → `J=∫ρv_b·dΣ` surface-independent but `v_b` accelerates (nozzle); no accepted `J_crit`
   → `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`. Conditional ideal Euler-nozzle numbers (`c_{s,*}/c_s0=3^{−1/2}`, etc.) are
   **recorded, NOT accepted**.
7. **ħ provenance (UNDETERMINED).** `HBAR_PROVENANCE_UNDETERMINED`; anti-tautology: `ħ=m c_s0 a` is a pin rearrangement unless `a` is
   independently fixed by an ħ-free relation; `h/2π` meaningful only in winding/cycle bookkeeping.

**pathA_20b (coupled linearization — sharpens claim 4; reference-only per blueprint §4, but its negative control is the able-to-fail teeth):**
8. Coupled principal symbol on the neutralized homogeneous background: phonon block `det P_ph = ħ(ω² − (5Kρ0⁴/m)k²)`; gauge transverse
   block `P_T = C_E ω² − C_B k² = C_E(ω² − c_bulk²k²)` with `c_bulk²=C_B/C_E`; coupled `det = P_ph·P_T²`; **off-diagonal principal terms
   VANISH** on the neutralized bg. Two dispersions: phonon `ω²=c_s0²k²` (+`k⁴` Bogoliubov), gauge `ω²=c_bulk²k²`. `[c_s]=[c_γ]=LT⁻¹`
   explicitly **non-evidentiary for equality**.
9. Two-layer verdict: bulk `C_GAMMA_BULK_UNDERDETERMINED` (`BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`); brane
   `C_GAMMA_RATIO_STILL_UNDERDETERMINED` (`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`, `BRANE_PHOTON_CONE_REQUIRES_PROFILE`). If `c_bulk` is
   ρ0-independent, `λγ ∝ ρ0⁻²`.
10. **Negative control (the able-to-fail teeth):** the forbidden forced equality `C_B/C_E = 5Kρ0⁴/m` →
    `FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION` (independent symbols `c_bulk`, `c_s0`; no source equation found in the parent action).

## Scope rollup
- **EARNED:** `c_s²=5Kρ⁴/m` + state-dependence (1); velocity dims (2); `c=c_γ` ceiling + time dilation (3); coupled principal symbol +
  both dispersions (8).
- **POSTULATED (labeled):** EOS `P=Kρ⁵` closure (1); mass-bridge candidate (5); the neutralizing external source (background of 8).
- **CALIBRATED / free:** `λγ=c_γ/c_s` + `λγ³` tail (4, 9) — the headline, unpinned by the action, ρ0-dependent.
- **UNDETERMINED / open:** flux law `J_crit` (6); ħ provenance + `h/2π` split (7); bulk + brane cone normalization (9).

## Carried residuals (first-class, not repaired — the directive's §3)
pathA_20: `EOS_CLOSURE_IMPOSED`, `C_GAMMA_RATIO_UNDERDETERMINED` (carry λγ + λγ³), `STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA`,
`NO_NET_ACCRETION_BC_UNDERIVED`, `HBAR_FREE_SUBSTRATE_RELATION_MISSING` (BLOCKS_HBAR_EMERGENT), `H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED`.
pathA_20b: `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`, `PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING` (BLOCKS_BULK_EQUALS_C_S — no
EQUALS verdict absent a source equation), `BRANE_ZERO_MODE_REDUCTION_UNDERIVED`, `BRANE_PHOTON_CONE_REQUIRES_PROFILE`.

## Reshape notes (how I-2 differs from I-1)
- **Consuming stage** (like stage002 consumed stage001): consumes I-1 (`ledger_stage004`) dictionary + the `EOS_FROM_GNLS_FACTOR`
  handoff (`h0=(m c_s0²)/4`, `ξ_h=√2 ħ/(m c_s0)`) as CITED symbolic inputs — do NOT re-derive them here.
- **Extract** a standalone `_sympy_audit.py` from the harness + **re-author an independent `.wl`** (different route — e.g. native
  `Solve`/`Series` for the dispersion determinants).
- **Able-to-fail teeth:** the pathA_20b `FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION` negative control (assert that forcing
  `c_bulk=c_s` without a source equation is rejected — proves λγ is genuinely free, not smuggled to 1) + dimensional firewall + a
  corrupted-EOS probe that breaks `c_s²=5Kρ⁴/m`.
- **Represent λγ as a labeled free calibration input** (like I-1's labeled non-derivations) — the honest CALIBRATED landing, first-class.
- **Decision reaffirmed (blueprint §4):** pathA_20b is folded as **reference-only support in the ONE stage I-2** (its coupled-linearization
  negative control pulled in as teeth), NOT a separate stage.
- Target stem: `ledger_stage005_sound_speed_light_ratio` (or similar; confirm slug when authoring).

## Process (unchanged, calibrated — same as I-1)
Author reshape directive → Codex xhigh design-review → fold to `DIRECTIVE_CLEAN` (no GLM on Parts I–VI) → **pre-exec user gate** →
Codex build (`--sandbox danger-full-access`, background, `< /dev/null`, xhigh) → dual-engine both exit 0 → orchestrator arbiter re-run
(via runners) → full tri-review on fresh agents (arbiter + fidelity + adversarial-scoped-to-reshape-integrity, with ablation) →
remediate → bump counts 4→5 → note/card/LaTeX(Part I appendix already exists → just add `\input{stages/stage_005}`)/registration →
rebuild PDF → commit + docs/memory sync. ⛔ The I-1 exemplars that served as directive/prompt templates are **not retained** —
they lived in gitignored `_scratch/ledger_stage004_*` and no copy survives.
