# Ledger exclusions — the "approaches that didn't work" paper backlog

**What this is.** The tracked backlog of material **deliberately EXCLUDED from the v2 ledger** under the standing rule
(user, 2026-07-22): *the ledger shows the SURVIVING solution — clean, direct, simple, groundable; post-mortems of
tried-and-discarded/retired approaches go to a separate future "approaches that didn't work" paper.* Nothing here is
lost — the source material is preserved (permanent `software/…` reports + banked seeds); this note is the index so the
future paper can be assembled and so the ledger stays clean.

**The line** (see memory `feedback_ledger_surviving_solution_only`): KEEP in the ledger = a *feature/departure of the
surviving model*; MOVE here = a *post-mortem of a discarded/retired mechanism*. A discarded approach's one-line
*surviving residue* stays in the ledger (e.g. "μ_R is postulated"); its detailed apparatus comes here.

---

## Exclusion 1 — the couple-stress no-go (`pathA_35` gateL): why the retired polar field `P` cannot source light's `μ_R`

- **Excluded 2026-07-22** (first application of the standing rule). Was slated as ledger stages 030/031 (III-2/III-3);
  the per-Part gate had approved a 3-stage Part III, but on review the couple-stress no-go is a **post-mortem of the
  retired `P`** (Decision 16), not a feature of the surviving model → excluded.
- **Surviving residue that STAYS in the ledger** (stage003 + parameter register): *"`μ_R` is a postulated MacCullagh
  modulus; a derivation from a deeper polar substructure was attempted and failed, so `μ_R` remains an input knob."*
  One clean honest-scope line — NOT the apparatus.
- **What moves to the failures paper** (the apparatus): the retired-`P` couple-stress analysis — `ARROWS_SUPPLY_TRACTION`
  (rank-2 cut traction), the surface projections `J_P=mρa²ℓ_g` / `Γ_P=mρc_s²a²ℓ_g` / `M_gap_P`, the 7×7 principal
  symbol + signed-root branch counts, the three `P` configs (live → 3 hidden spin waves + unbounded minor
  `k²(Γ_P μ_R k² − λ_Pu²)`; gapped → residual `2 μ_R Ω_u`; slaved → `μ_eff=μ_R−2λ_Pu` + `POSTULATED_SURFACE_ELASTICITY`),
  the C5/`φ` obstruction, and the aggregated `FAIL_COUPLE_STRESS_NOGO`. This is Decision-16's "light" payoff-failure
  evidence — perfect failures-paper content, wrong for the surviving-solution ledger.
- **Preserved source material (nothing lost):**
  - Permanent: `software/stage1_solver/reports/pathA_35_gateL_light.md` + `tools/pathA_35_gateL_{sympy.py,.wl}`.
  - Decision record: `software/stage1_solver/decisions/16_retire_brane_polar_field.md` (evidence items 2 "Light" + 4).
  - **Durable seed (tracked):** the full source distillation is folded into the tracked
    `notes/stage030_pathA35_gateL_source_map.md` (the derivation content the paper needs). *(A session-transient
    physics-verified reshape directive was also produced in gitignored `_scratch/` — cleared Codex→Grok→Codex
    `DIRECTIVE_SOUND` + `GROK_COMPUTE_CLEAN` — but `_scratch/` is NOT tracked, so it is a bonus, not the preservation:
    the tracked source map + the permanent `pathA_35` reports are the durable record.)*

---

## Exclusion 2 — superseded EM (charge + magnetism) routes (stub; full write-up when Parts IV/V are re-scoped)

- **Excluded 2026-07-22** (model_map §3.4 + blueprint Part IV/V rows point here). These are **discarded EM mechanisms**
  superseded by the surviving puncture-deflection (charge) + moving-throat (magnetism) builds; carried as one-line
  pointers in the ledger, apparatus → failures paper.
- **Surviving residue in the ledger:** the surviving charge mechanism is the static ±w throat (puncture-deflection,
  `R1_REQUIRED(bc_selection)`); magnetism is the moving throat. These predecessors are named as superseded only.
- **What moves to the failures paper (preserved source, permanent):**
  - **Charge:** leftover-scalar `u_L`-clamp (`NO_NATIVE_CLAMP`) → `software/em_charge_attribute/leftover_scalar_electric_sign_result.md`;
    defect/antidefect wall-healing → `docs/em_phaseC_force_decomposition.md`; the old `pathA_38` throat-body Coulomb
    anchor (`THROAT_ELECTRIC_LOCALIZED_COULOMB`) → `software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md`.
  - **Magnetism:** the old `pathA_39` route resting on the barred `j∝sV` → `software/stage1_solver/reports/pathA_39_*` (superseded by the natively-derived `J_T=q_T sηV` moving-throat build).
- Full write-up deferred to the Part IV/V re-scopes; this stub makes the model_map §3.4 + blueprint pointers resolve.

## Exclusion 3 — the retired-`P` action machinery (stages 006/007; stub)

- **Excluded 2026-07-22.** The action-level machinery of the retired polar field `P` (Decision 16), beyond the
  couple-stress no-go (Exclusion 1): stage006's `α_aniso` easy-plane anisotropy (`α_aniso χ_B (P·ŵ)²`) and stage007's
  `λ_Pu`/`L_pol`/`L_Pu` terms + the historical `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`/`DOF=8` dual-tier drift apparatus.
- **Surviving residue in the ledger:** stage006 operative `DRIFT(5)`; stage007 operative `DOF=4` / `POST_D16_DRIFT(7)`
  + operative action `{S_GNLS,gL_Mac,gL_uw}`. The `α_aniso`/`λ_Pu` pointers in stage006/007/register resolve here.
- **Preserved source (nothing lost):** the authoritative decision `software/stage1_solver/decisions/16_retire_brane_polar_field.md`;
  the stage006/007 **audit scripts still compute the historical tier** as verification provenance (unchanged); git
  history retains the full pre-trim notes/cards. Full write-up deferred; this stub makes the retirement pointers resolve.

## Graveyard audit — DONE (2026-07-22)
The focused audit of the already-built stages/docs ran (covering stages 001–029 + the maps/trackers) and the
user-approved **"docs-only, keep hash"** cleanup landed: reader-facing prose trimmed to operative-only + one-line
provenance across 11 files (heaviest = the `stage007` note+card retired-`P` dual-tier apparatus → operative freeze;
`model_map` §3.3/§4/§3.4/§6; `parameter_register`; `stage006` α_aniso; blueprint stale EM rows; two stale gravity
cross-refs), PDF rebuilt. **No script edited, no hash changed** — the audit scripts keep computing the historical tier
as verification, flagged by a one-line "the script also computes the historical tier" note in each trimmed doc.
**Confirmed KEEP (surviving, not graveyard):** the θ-as-Maxwell-φ `FATAL_FLAW`, `FAIL_CAUCHY_STRAY_LONGITUDINAL`,
gravity's honest-negatives, the density-port vector-freedom proof, stages 001/002. **Gravity (Part II) came back
essentially clean** — its honest-negatives + retirement-proofs are surviving methodology/features.
