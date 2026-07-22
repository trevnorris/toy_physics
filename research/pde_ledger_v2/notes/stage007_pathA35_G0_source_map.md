# I-4 (ledger_stage007) source map — pathA_35 G0 shear-surface freeze

> Running-start prep for the NEXT stage (I-4), captured 2026-07-07 from a source-mapping fan-out so a fresh session can
> author the reshape directive without re-discovery. Verify against the cited sources before finalizing the directive.
> Companion: `part1_medium_atomic_split.md` (row 007). Build-order id **007**, Part I — the LAST Part-I stage
> (after it: the Part II gravity user gate).
> **Reshape of an EXISTING dual-engine pair** (NOT fresh-authored): `software/stage1_solver/tools/pathA_35_G0_sympy.py`
> (~645 lines) + `pathA_35_G0.wl` (~501 lines), `ENGINE_AGREE`.
> **Headline verdicts (verbatim):** `T0_SHEAR_FROZEN(d9520d3819c3)` + `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`.
> Full frozen-action SHA-256: `d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3`.
> **⚠ AMENDED 2026-07-21 (Decision 16):** the ledger stage has since retired the polar field `P` — `L_pol`, `L_Pu`,
> `λ_Pu`, and structural postulates 3/4/5. The freeze + verdicts here are the immutable HISTORICAL tier (still
> hash-verified, 11/DOF=8 STAND); the operative post-D16 stage carries DOF=4 and `POST_D16_DRIFT(7)`. See the stage
> note's read-first Decision-16 layer.

## File inventory
- **Freeze report (the authoritative source):** `software/stage1_solver/reports/pathA_35_G0_freeze.md` (17 KB) —
  scope line 9: "reports-only plus dual-engine dimensional and flat-brane DOF count checks; **no Gate L verdict
  computed**". Machine ledger: `reports/pathA_35_G0_results.yaml`.
- **Scripts to reshape:** `tools/pathA_35_G0_sympy.py` + `tools/pathA_35_G0.wl` (see §2 for their structure).
- **Directive (methodology to carry):** `directives/pathA_35_shear_surface_brane_gates.md` (43 KB; §1 binding
  methodology :38–62, §7 G0 protocol :249–282, §10 dimensional firewall :344–360).
- **The 2026-07-04 erratum:** appended to the freeze report itself (:323–336, append-only) — see §4.
- **Gate-L artifacts (EXCLUDED from I-4 — name them in the directive's What-NOT-to-do):**
  `tools/pathA_35_gateL_sympy.py` (54 KB) + `pathA_35_gateL.wl` (43 KB); `reports/pathA_35_gateL_light.md`
  (`FAIL_COUPLE_STRESS_NOGO`) + `_results.yaml`. The G0 scripts compute NO Gate-L content (only prose strings naming
  the exposures that "remain able to fire"; scope field `.py:610`, `classification_guard` `.py:579`).

## §1 The frozen content (what the stage represents)
- **Kept (0 new):** GNLS + the T0 polar-OP action `L_pol = ½mρa²(D_t^v Pⁱ)² − ½mρc_s²a²(∂_j Pⁱ)² − ¼mρc_s²(PⁱPⁱ−1)²`
  (T0 block hash `8fa41ac51e88…`, byte-embedded in the G0 block — kept byte-for-byte from `pathA_24_T0_freeze.md`).
- **Frozen brane/light action** (canonical block, report :200–261):
  - `L_Mac = ½ρ_br(∂_t uᵃ)² − ½μ_R Ω_uᵃΩ_uᵃ` (MacCullagh rotational shear);
  - `L_Pu = −λ_Pu ϖ_a Ω_uᵃ = −λ_Pu ŵ·(P_∥×(∇_∥×u))`, `ϖ_a := (ŵ×P_∥)_a` — parity-EVEN (the direct `P_∥·Ω_u` is
    parity-ODD, excluded); ⚠ preserve verbatim: this operator "re-admits the ε-contracted/chiral class excluded by
    T0" and REQUIRES the conceded axis `ŵ` (report :58, :227) — a structural-postulate cost, not a free choice;
  - `L_uw = ½ρ_br(∂_t u_w)² − ½ρ_br Ω_w² u_w²` (the u_w gap);
  - `S_brane = ∫dt d⁴X g_ℓ(w)[L_Mac + L_Pu + L_uw]`, `g_ℓ(w) = exp(−(w/ℓ_g)²)/(√π ℓ_g)`, `∫g dw = 1`.
- **What the hash is OF:** the exact bytes inside the single `freeze-action` fenced block (excluding fences), 1-based
  byte range 8111–13310, length 5200 (report :280; extraction command :268–269). The scripts re-extract + SHA-256 +
  assert equality with the expected hash (freeze-fidelity check).

## §2 The script pair — what the reshape strips/keeps
- **STRIP (blueprint §5):** `argparse --compare` (`.py:11,641–644`); the `agreement_payload` JSON writes
  (`SCRATCH/pathA_35_G0_sympy.json`; `.wl:478–490,500` `WriteString`) — the engines currently agree via a `--compare`
  JSON diff = a **payload mirror**. Re-architect to the v2 convention: each engine standalone, print-only, asserts
  independently in-process, both exit 0. **This is the single biggest structural change.**
- **KEEP (the load-bearing checks, per block):**
  - `check_freeze_fidelity` (`.py:108–125`): re-extract the fenced block → SHA-256 == `d9520d3819c3…` + the T0 hash
    `8fa41ac51e88…` byte-embedded.
  - `build_dimension_payload` (`.py:185–433`): inline dim audit of EVERY frozen term — kept GNLS, `L_pol`,
    profile+measure, brane inertia, MacCullagh, reused couple-stress coeffs, parity-repaired P–u, u_w gap, the full
    projected traction `T_na` incl. `T_wa = mρ v_w v_a`, action measures, and all G0.5 linearization quantities
    (`O_u`, `c_γ²`, `ω_P²`, `ω_radial²`, `ω_uw²`).
  - `compute_flat_brane_dof` (`.py:439–501`) + `flat_brane_modes` (`.py:517–580`): DOF=8 **genuinely rank-computed**
    (curl projector, tangent/radial P, u_w, φ) and wired reported→computed (assert total==8 at `.py:532`).
    Breakdown (report G0.5 :110–129): in-plane uᵃ curl-transverse 2 + kinetic−curl 1; T0 Pⁱ tangent 3 + radial 1;
    u_w 1; C5 φ 0.
  - Able-to-fail ablations (all FIRED, report :301–308): dim — `drop_m_from_T_wa`, `MacCullagh_without_curl`;
    DOF — `drop_u_w_gap_term`, `drop_P_soft_spin_radial_term`, `zero_u_longitudinal_component` (each 8→7).
- The `.wl` is currently a MIRROR of the `.py` (same payload construction) — the reshape must re-author it as a
  genuinely independent route (mirror policy), e.g. its own dim machinery + its own rank/projector construction.

## §3 The "11" drift enumeration (report G0.4 :68–108) — represent verbatim
- **4 constants:** `ρ_br` [M L⁻³] surface inertia · `μ_R` [M L⁻¹T⁻²] MacCullagh modulus · `λ_Pu` [M L⁻¹T⁻²]
  parity-repaired P–u coupling · `Ω_w` [T⁻¹] bare u_w gap scale. (T0 couple-stress coeffs = KEPT, 0 new.)
- **1 function:** `g_ℓ(w)` — fixed Gaussian shape + one width `ℓ_g`; no free-form profile/kernel.
- **6 structural postulates:** (1) imposed `ŵ` axis + `w=0` surface (conceded-wall); (2) `uᵃ` same-medium surface
  collective, tangentially free-slip (`u̇ᵃ ≠ vᵃ`); (3) T0 `Pⁱ` reused as the Cosserat micro-rotation reservoir
  (0 new DOF); (4) baseline `Pⁱ` spin-wave status = `massless` (alternates `gapped`/`slaved-rigid` named-inactive);
  (5) the `ŵ`-dependent parity-even P–u operator (re-admits the ε-contracted/chiral class); (6) no C5 `φ` analog /
  no longitudinal constraint.
- Field/DOF sub-count new at G0 = 4 (uᵃ=3, u_w=1) — SEPARATE from the 11-input drift count.

## §4 The 2026-07-04 erratum (carry as first-class; report :323–336)
Operative correction (quote): "**The `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` count … is NOT inflated by a `ρ_br`
overcount — this freeze's count stands.**" The earlier GLM claim that `ρ_br` duplicated pathA_25's derived
`varrho_br[ρ]` is superseded: pathA_25's object belongs to the CLOSED density-smectic candidate (`FAIL_NOT_CODIM1`),
`OUT_OF_ACTIVE_NG5`. This `ρ_br`/`μ_R` = genuine postulated shear-surface inertia/modulus with a registered-pending
pathA_40 **Route-A** reduction (`ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` = register R10). "**The honest
cross-sector drift** (per pathA_41) is instead `{ρ_B0, χ_c, C_hu}`" — a Part-VI item, NOT absorbed into the "11".
Corroboration: `reports/pathA_41_ng5_second_medium_drift.md:1–7,22–26,48–52` (`NO_OVERCOUNT_ROUTE_A_PENDING`).

## §5 Downstream consumers (mark as load-bearing exports)
- **`ledger_stage003`** consumes the physics: `μ_R`, `ρ_br` symbols (its audit `.py:118–119,148–149`),
  `c_γ² = μ_R/ρ_br` (`:588`), the transverse dispersion (`:587`) — the frozen `L_Mac` IS stage003's starting
  Lagrangian; the token is cited in its note (`ledger_stage003_...md:34,53`) and results. I-4 is the token's formal
  home (split note :35–36).
- **`ledger_stage006`** cites `c_γ² = μ_R/ρ_br` + the χ_B relationship (see trip-up 5).
- **Parameter register:** rows already present — `μ_R` (:verified stage 003), `ρ_br`, edge R4, Route-A R10, the
  pathA_41 trio. **Stage007 ADDS:** `λ_Pu` [M L⁻¹T⁻²], `Ω_w` [T⁻¹], `g_ℓ` [L⁻¹] (function) + `ℓ_g` [L] rows, the
  6 structural postulates as a noted structural block, and re-homes `μ_R`/`ρ_br` provenance to I-4 (keeping the
  stage-003 dim verification attribution).

## §6 Dimensional targets (as the scripts assert; sympy `build_dimension_payload` :204–219, :273, :372)
`Pⁱ` dimensionless · `uᵃ, u_w` = L · `ρ_br` = M L⁻³ · `μ_R, λ_Pu` = M L⁻¹T⁻² · `Ω_w` = T⁻¹ · `g_ℓ(w)` = L⁻¹ ·
`ℓ_g` = L · `c_γ² = μ_R/ρ_br` = L²T⁻².

## §7 Directive methodology to carry (binding, from `pathA_35_shear_surface_brane_gates.md`)
- §1 :38–62: postulate the FULL structure up front; "postulating an *ingredient* is allowed; postulating an *outcome*
  is not"; late ingredient = `AD_HOC_RESCUE` → fresh G0; every knob counted; ≥2 new inputs ⇒ `SECOND_MEDIUM_DRIFT`
  reported plainly; a clean all-pass is suspicious.
- G0.2 target-blind rule: `g(w)` admitted on locality/minimality grounds ONLY, never because it helps a gate.
- Anti-impose: G0 freezes *terms*, NOT gate answers (no bounded-below / traction / longitudinal-is-gauge assertions).
- §10: dual-engine dim checks on REAL expressions, able-to-fail, tautology-free.

## §8 Reshape trip-ups (pin in the directive)
1. **The "11" is a HARDCODED literal** — `ledger()` (`.py:583–592`) returns subcounts `4/4/1/6/11/8` + the verdict
   string as literals; the 11 is NOT recomputed from the enumeration (the DOF=8 IS genuinely rank-computed). The
   reshape must either DERIVE the subcounts from an enumerated ledger table (preferred — matches stage006's computed
   DRIFT tally) or explicitly label them bookkeeping literals. Do not let the verdict be a stamp.
2. **Payload-mirror → independent engines** (§2 above) — the biggest structural change; the `.wl` needs a genuinely
   independent route, not the shared payload.
3. **Hash fidelity is byte-range-sensitive:** the `freeze-action` block byte range (8111–13310, len 5200) is
   report-revision-specific; the standalone stage script must re-extract by fence-parsing (not fixed offsets) or
   embed the frozen bytes. Any edit ABOVE the block in the report breaks fixed offsets (the erratum was appended
   below, safely).
4. **μ_R notational collision:** I-4's `μ_R` is the 3D BRANE modulus (M L⁻¹T⁻²) — distinct from stage006's 4D gate
   density `μ_R⁽⁴⁾` (M L⁻²T⁻²), related by register edge R17 `μ_R = ∫χ_B μ_R⁽⁴⁾ dw` (PENDING). Do not conflate.
5. **χ_B-wall vs frozen g_ℓ(w) supersession (stage006 relationship):** the stage006 χ_B order-field wall superseded
   the fixed-shape `g_ℓ(w)` profile as the *material-state* closure, but the G0 freeze REMAINS the light-sector
   **constitutive** freeze (MacCullagh/P–u/u_w action). The directive must state both facts and not claim `g_ℓ(w)`
   as the material wall (nor let the χ_B wall retro-invalidate the freeze — stage003 consumes it as-is).
6. **Do NOT absorb `{ρ_B0, χ_c, C_hu}` into the "11"** — Part-VI (pathA_41) item per the erratum.
7. **Gate-L exclusion:** I-4 stops at the freeze; no gate verdict recomputed/imported (the exposure-name strings in
   the G0 scripts are prose, not computations — keep them as printed provenance only).

## Verdict tokens + honest scope
Headline `T0_SHEAR_FROZEN(d9520d3819c3)` + `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`. EARNED: freeze-fidelity (hash),
flat-brane linear DOF=8 (rank-computed, able-to-fail), dimensional firewall. POSTULATED/CALIBRATED: the 11 (labeled,
enumerated, count computed in the reshape). Carried: Gate-L exposures (prose provenance), the 2026-07-04 erratum,
Route-A pending. Does NOT earn light (stage003's job), does NOT compute any gate verdict.

## Process (reshape variant — same calibrated per-stage flow)
Author directive (reshape contract + the §8 pins + faithful-cover of §2's check blocks) → Codex xhigh design-review →
fold to `DIRECTIVE_CLEAN` (no GLM on Parts I–VI) → pre-exec USER GATE → Codex builds the two scripts
(`--sandbox danger-full-access`, background, xhigh) → dual-engine both exit 0 → arbiter re-run via runners → full
tri-review on fresh agents (arbiter + fidelity + adversarial-scoped-to-reshape-integrity, with ablation) → remediate →
bump counts 6→7 → update + Codex-verify the parameter register → note/card/`\input{stages/stage_007}`/registration →
PDF → commit + docs/memory sync. Target stem: `ledger_stage007_shear_surface_g0_freeze` (confirm slug on authoring).
**After I-4: Part I COMPLETE → the Part II (gravity) atomic-split user gate.**
