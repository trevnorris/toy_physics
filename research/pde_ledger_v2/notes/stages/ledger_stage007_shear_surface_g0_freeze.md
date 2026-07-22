# ledger_stage007_shear_surface_g0_freeze

## ⚠ Decision-16 amendment (2026-07-21) — the P-retirement layer (READ FIRST)

The operative frozen action is `{S_GNLS, gL_Mac, gL_uw}` under `g_ℓ`, with operative **DOF = 4**,
**`POST_D16_DRIFT(7)`** (= 3 constants `{ρ_br, μ_R, Ω_w}` + 1 function `g_ℓ` + 3 survivor postulates `{1, 2, 6}`),
and `c_γ² = μ_R/ρ_br`. A polar-field `P` substructure (`L_pol`/`L_Pu`/`λ_Pu`, postulates 3/4/5) was tried as a
light-stiffness source and RETIRED (Decision 16, `software/stage1_solver/decisions/16_retire_brane_polar_field.md`);
its analysis lives in the failures-paper backlog (`notes/ledger_exclusions_failures_paper_backlog.md`). The freeze
hash `d9520d3819c3` byte-includes those now-retired `P` terms: the block is **unchanged on disk** and both engines
still fence-parse and SHA-256-verify it (the retirement is a symbolic action-summand set partition over the immutable
hash anchor, no byte surgery), so the historical freeze-as-run record STANDS immutable. `c_γ² = μ_R/ρ_br`, `L_Mac`,
`L_uw`, `g_ℓ`, `{ρ_br, μ_R, Ω_w, ℓ_g}`, and the R10 Route-A reduction are untouched; Part III/light consumes `L_Mac`
as-is (no `P` sector). The `χ_B = |P_∥|²` route (c) stays a named, high-risk, Part-VII-adjacent future gate needing a
new T0 freeze (obsolete-as-carried, not foreclosed).

Ledger earned-label for the amended stage: `G0_FREEZE_FIDELITY_PLUS_POST_D16_LAYER_VERIFIED`. Dual-engine after the
amendment: **SymPy 142 PASS / Mathematica 140 PASS**, both exit 0, CWD-independent.

## Status

**Part I — The medium. I-4 (build-order 007; the LAST Part-I stage).** Reshape of the existing dual-engine gate
`pathA_35` **G0** (the shear-surface brane/light-sector **constitutive freeze**). Surviving headline:

- **`T0_SHEAR_FROZEN(d9520d3819c3)`** — the frozen brane/light action block, SHA-256
  `d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3`, with the T0 polar-OP block (short hash
  `8fa41ac51e88`) kept byte-for-byte and byte-embedded inside it.
- **`POST_D16_DRIFT(7)`** — the honest operative drift ledger: 7 independent new inputs, **computed in both engines
  from an enumerated table** (3 constants + 1 function + 3 survivor postulates), never a typed literal.

**EARNED** (what this stage genuinely computes): freeze-fidelity (fence-parsed byte-exact hash re-verification of
both freeze blocks + the T0 embedding), the flat-brane linear operative **DOF = 4** (rank-computed, able-to-fail), and
the dimensional firewall over every operative frozen term. **POSTULATED/CALIBRATED** (the honest landing,
first-class): the 7 operative freeze inputs. The freeze freezes **terms, not gate answers** (anti-impose); this stage
computes **no gate verdict** (Gate-L excluded, below) and does **not** earn light (stage003's job — stage003 derives
the 2 transverse photons at `c_γ² = μ_R/ρ_br` from this stage's frozen `L_Mac`).

*(The audit script additionally computes the historical freeze-as-run tier — `DOF=8`,
`SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` — as verification provenance; not reproduced here.)*

Ledger-local earned-label (NOT a source verdict token): `G0_FREEZE_FIDELITY_DOF_DIM_VERIFIED`.

## Purpose

G0 is where the program **postulated the full shear-surface brane structure up front and paid for it in public**:
every new constant, the one profile function, and every structural concession counted into a single drift number,
frozen *before* any gate was run (the pathA_35 directive's binding rule — "postulating an *ingredient* is allowed;
postulating an *outcome* is not"; a late ingredient = `AD_HOC_RESCUE` → fresh G0; a clean all-pass is suspicious).
The rebuilt ledger needs this stage as the **formal home** of the freeze token that Part III consumes, and as the
Part-I closing statement of what the medium's brane sector *costs*.

## The frozen action (canonical block content; hash-guarded)

The historical freeze declares five action summands; the **operative** action is `{S_GNLS, gL_Mac, gL_uw}`. The
`P`-dependent complement `{L_pol, gL_Pu}` (with the constant `λ_Pu`) is RETIRED by Decision 16 (→ failures-paper
backlog); it remains byte-embedded in the frozen block as a hash-verified historical record only. The byte-level
SHA-256 anchor is unchanged — the partition is symbolic.

**⟨operative, SURVIVES⟩** the GNLS parent action `S_GNLS` (I-1/I-2), 0 new.

**Frozen brane/light blocks** (`S_brane = ∫ dt d⁴X g_ℓ(w) [L_Mac + L_uw]`, operative):

```
L_Mac = ½ ρ_br (∂_t uᵃ)² − ½ μ_R Ω_uᵃ Ω_uᵃ          (MacCullagh rotational shear; Ω_u = ∇_∥×u)   ⟨operative⟩
L_uw  = ½ ρ_br (∂_t u_w)² − ½ ρ_br Ω_w² u_w²          (the u_w gap)                                  ⟨operative⟩
g_ℓ(w) = exp(−(w/ℓ_g)²)/(√π ℓ_g),   ∫ g_ℓ dw = 1      (fixed Gaussian profile, one width ℓ_g)        ⟨operative⟩
```

The operative `ŵ` is the **intrinsic wall normal** (postulate 1, annotation softened; no longer a concession *for*
the retired P–u operator).

## Freeze fidelity (EARNED)

Both engines re-extract the `freeze-action` fenced blocks from the two canonical reports
(`software/stage1_solver/reports/pathA_35_G0_freeze.md`, `reports/pathA_24_T0_freeze.md`) **by byte-exact
fence-parsing** — never by fixed byte offsets — and assert SHA-256 equality with the frozen constants plus the
byte-level embedding of the T0 block inside the G0 block. Each engine uses its own extractor (the SymPy route is an
exact-match line scanner; the Mathematica route is an independent byte-pattern scanner over the raw byte stream).
Teeth: a one-byte in-memory corruption of the extracted block trips the SHA mismatch; a nonexistent fence tag trips
the extractor's missing-fence path; adversarial ablation additionally confirmed a one-byte edit *outside* the block
(in the erratum region) leaves the stage green — fence-parsing, not whole-file hashing. The short hash in the
headline token is asserted to be a prefix of the verified full constant.

*Byte-range convention remark:* the scripts print the freshly-computed content byte range informatively; it reads
length **5201** where the committed report metadata says **5200** (1-based 8111–13310). The report's historical
range convention excluded the block's final newline while the hash always included it; the reshaped scripts report
the range of the actually-hashed bytes. Never asserted either way — the SHA-256 match against the frozen constant
is the fidelity anchor.

## Dimensional firewall (EARNED; full source surface carried)

Exponent-triple `{L,T,M}` audit of every operative frozen term, faithful to the source gate's surface (fidelity-leg
coverage-diffed, nothing dropped): kept GNLS (`c_s²=5Kρ⁴/m`, `U=(K/4)ρ⁵`, quantum pressure, `mρvv`), profile +
measure (with `∫g_ℓ dw = 1` genuinely integrated in both engines), the two operative brane blocks `L_Mac` and `L_uw`
(each term `M L⁻¹ T⁻²`, the 3D surface density; `g_ℓ·[brane]` restores the 4D density), the action measures, the full
projected traction `T_na` including `T_wa = mρ v_w v_a`, `O_u`, `c_γ² = μ_R/ρ_br`, and the surviving linearization
mode `ω_uw,bare² = Ω_w²`.

The retired `P`-sector terms (`L_pol`, `L_Pu`, `λ_Pu`, and the retired `P`-mode frequencies) are not live-checked as
survivors; their dim analysis lives in the failures-paper backlog and the audit script's historical tier.

Targets: `uᵃ, u_w` = `L` · `ρ_br` = `M L⁻³` · `μ_R` = `M L⁻¹ T⁻²` · `Ω_w` = `T⁻¹` · `g_ℓ` = `L⁻¹` · `ℓ_g` = `L`.
Teeth (both fired in the source gate and re-fire here): `drop_m_from_T_wa`, `MacCullagh_without_curl`.

**Notational firewall (register edge R22):** this stage's `μ_R` (3D brane modulus, `M L⁻¹ T⁻²`) and stage006's
`μ_R⁽⁴⁾` (4D shear-stiffness density, `M L⁻² T⁻²`) are asserted **dimensionally distinct** (a computed
exponent-triple inequality, with its own tooth), and the R17 projection `μ_R = ∫ χ_B μ_R⁽⁴⁾ dw` is asserted
dim-consistent (`[μ_R⁽⁴⁾]·L = [μ_R]`) — status **PENDING**. Two symbols, one pending reduction; never conflated.

## Flat-brane linear DOF = 4 (operative, EARNED, rank-computed)

**Operative post-Decision-16 DOF = 4**, computed from rank/nullity bookkeeping in both engines (never typed) on the
retired field set (`Pⁱ` removed): `uᵃ` curl-transverse 2 + kinetic−curl 1 + `u_w` 1 + `φ` 0 = 4. The two engines use
genuinely different linear algebra (SymPy: the source gate's k-aligned curl-rank machinery; Mathematica: a generic
`k=(1,2,3)` outer-product transverse projector, ranked). The absent `φ` (no longitudinal-constraint analog —
structural postulate 6) is printed plainly: it is the fact the later C5 crux attacks (stage006's θ-as-φ no-go; the
material-state pivot). Operative teeth: reinject one `Pⁱ` mode → 5, the full `Pⁱ` block → 8, and the survivor
ablations `drop_u_w_gap_term` / `zero_u_longitudinal_component` → 3 (each ≠ 4). The new-field content at G0 (`uᵃ` 3 +
`u_w` 1 = 4 DOF) is computed separately from an enumerated field list and kept **out** of the drift count.

*(The audit script additionally computes the historical freeze-as-run tier — `DOF=8`,
`SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` — as verification provenance; not reproduced here.)*

## The drift ledger: `POST_D16_DRIFT(7)`, computed (POSTULATED/CALIBRATED — the honest landing)

The operative post-Decision-16 drift **`POST_D16_DRIFT(7)`** is **computed** (not typed) in both engines from an
explicit enumeration table and asserted equal to the frozen token — the source gate's hardcoded literals are retired:

- **3 constants:** `ρ_br` [M L⁻³] surface inertia · `μ_R` [M L⁻¹T⁻²] MacCullagh modulus · `Ω_w` [T⁻¹] bare u_w gap
  scale. (Table dims are asserted equal to the independently-verified firewall dims — the table cannot silently
  misinform.)
- **1 function:** `g_ℓ(w)` — fixed Gaussian shape, ONE width knob `ℓ_g`; admitted on locality/minimality grounds
  only (target-blind G0.2); no free-form profile.
- **3 survivor structural postulates** `{1, 2, 6}`: (1) imposed `ŵ` axis + `w=0` surface (annotation softened to
  "intrinsic wall normal"); (2) `uᵃ` same-medium surface collective, tangentially free-slip (`u̇ᵃ ≠ vᵃ`); (6) no C5
  `φ` analog / no longitudinal constraint.

Postulates 3/4/5 and the constant `λ_Pu` are RETIRED with `P` (Decision 16 → failures-paper backlog);
`POST_D16_DRIFT(7)` is derived as the exact set partition `historical ∖ {λ_Pu, postulates 3/4/5}`. T0 couple-stress
coefficients: asserted 0 new (kept).

Teeth: leave `λ_Pu` live → n=8; leave any retired postulate → n=8; a **same-cardinality** `Ω_w`↔`λ_Pu` swap keeps
n=7 but **fails the set-partition assert**; drop a survivor → n=6; miscategorize `Ω_w`; corrupt `ρ_br` table-dim;
inject `ρ_B0` → the anti-absorption guard fires; corrupt n before verdict assembly → the verdict-string equality
fires. Adversarial ablation confirmed the checks are anchored to in-engine derivations and the frozen token
(dual-corruption of both engines fails both) — no cross-engine comparison anywhere; the old `--compare` payload
mirror is dead.

## Anti-absorption (2026-07-04 erratum residue)

`ρ_br`/`μ_R` are genuine postulated shear-surface inertia/modulus; Route-A reduction pending (R10,
`ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT`; corroboration token `NO_OVERCOUNT_ROUTE_A_PENDING`); the
cross-sector drift `{ρ_B0, χ_c, C_hu}` is a Part-VI (pathA_41) item, not absorbed here (guarded in-engine by the
anti-absorption assert). (The earlier `ρ_br`-overcount claim was superseded — `varrho_br` belongs to the CLOSED
density-smectic candidate, `OUT_OF_ACTIVE_NG5`, a different structure.)

## Supersession (both facts, exactly)

1. The stage006 χ_B order-field wall **superseded the fixed-shape `g_ℓ(w)` profile as the *material-state*
   closure** (the wall is a postulated field, not a fixed profile).
2. The G0 freeze **REMAINS the light-sector *constitutive* freeze** — the operative MacCullagh/`u_w` action; stage003
   consumes the frozen `L_Mac` as-is.

Neither artifact retro-invalidates the other (register edge R21 — a scope split, not a reduction; `ℓ_g` stays a
counted knob of the constitutive freeze).

## Gate-L exclusion

Scope carried verbatim: **"G0 freeze only; no Gate L verdict computed"**, plus the classification guard (counts
only; no gate verdict, no boundedness, no gauge, no leak claim). The exposure-name strings (the gates that "remain
able to fire") are printed as prose provenance only — nothing is computed from them; nothing is imported from the
Gate-L artifacts (`pathA_35_gateL_sympy.py`/`.wl`, `reports/pathA_35_gateL_light.md` = `FAIL_COUPLE_STRESS_NOGO`).
Gate-L is the retired-`P` couple-stress no-go — NOT a ledger stage; it lives in the failures-paper backlog
(Exclusion 1), not Part III.

## Downstream consumers (load-bearing exports)

- **`ledger_stage003` (Part III, light):** consumes `μ_R`, `ρ_br`, `c_γ² = μ_R/ρ_br`, and the frozen `L_Mac` — all
  **operative survivors**, so light is unaffected by the P-retirement (`pathA_36` derives both photons with no `P`
  sector). It cites the historical `T0_SHEAR_FROZEN(d9520d3819c3)` **and** the operative subset
  `POST_D16_ACTION{S_GNLS,gL_Mac,gL_uw}`. This stage is the token's formal home.
- **`ledger_stage006` (I-3):** cites `c_γ² = μ_R/ρ_br` and the supersession relationship above; its `μ_R⁽⁴⁾` is
  R17-related, R22-distinct.
- **Parameter register:** rows `Ω_w`, `g_ℓ`/`ℓ_g` live; `λ_Pu` **retired** (Decision 16); the structural block
  reduced to survivor postulates `{1,2,6}`; `μ_R`/`ρ_br` provenance re-homed to I-4 (stage-003 dim-verification
  attribution kept); edges R21 (scope split) and R22 (dimensional distinctness) retained.

## Verification

- **⭐ Decision-16 amendment (2026-07-21):** the P-retirement layer applied per
  `_scratch/decision16_amendment_directive.md` (directive cleared the Codex→Grok→Codex bookend — Codex
  `DIRECTIVE_SOUND`, Grok `GROK_COMPUTE_CLEAN`). Dual-engine after the amendment: **SymPy 142 PASS / Mathematica
  140 PASS**, both exit 0, CWD-independent; the four tracked transcripts regenerated under `scripts/output/` +
  `mathematica/output/`. ⏳ Fresh-agent tri-review of the amended scripts + docs is the next gate. The original-build
  record below is retained.
- **Reshape (blueprint §5):** argparse/`--compare`/JSON payload-mirror stripped; each engine standalone,
  print-only, in-process asserts, no file writes (file READS of the two canonical freeze reports only, paths
  resolved from the script location — verified from a foreign CWD); float-free symbolic payload; the `.wl`
  re-authored as a genuinely independent route (own fence scanner, own dim machinery, own projector/rank DOF
  construction, own enumeration tallies).
- **Dual-engine (original build):** SymPy audit **96 PASS / 0 FAIL**, exit 0 · Mathematica audit **94 PASS / 0
  FAIL**, exit 0 (both from repo root and from a foreign CWD). Superseded by the amendment tallies above.
- **Tri-review (fresh agents):** orchestrator arbiter re-run via the runners (7/7 PASS both engines);
  **`FIDELITY_CLEAN`** (block-by-block coverage diff vs the source `build_dimension_payload` — no dropped check, no
  changed target; all dims independently reconfirmed); **`ADVERSARIAL_CLEAN`** (26-run mutation matrix:
  computed-11 genuine, hash genuinely fence-parsed with the outside-the-block control, dual-corruption class fails
  both engines, teeth recompute, `.wl` independent, containment clean).
- **Remediation (tri-review nits folded):** the two inert X≡X DOF checks deleted; drift-table dims made
  load-bearing (asserted against firewall dims + new tooth); short-hash-prefix assert added; Route-A status token
  printed. (An interrupted remediation session transiently dropped the `.wl` drift section; caught by tally
  arithmetic + check-name diff and repaired; the section re-verified live.)

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_35_G0_{sympy.py,.wl}` (reshaped; sources unchanged).
- Frozen reports (hash-anchored at runtime): `software/stage1_solver/reports/pathA_35_G0_freeze.md` (incl. the
  2026-07-04 erratum) + `reports/pathA_35_G0_results.yaml`; `software/stage1_solver/reports/pathA_24_T0_freeze.md`.
- Methodology: `software/stage1_solver/directives/pathA_35_shear_surface_brane_gates.md` (§1, G0.2, §7, §10).
- Reshape directive + tri-review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage007_*` +
  `_scratch/adv_stage007/`.
- **Decision-16 amendment:** `software/stage1_solver/decisions/16_retire_brane_polar_field.md` (the decision) +
  `research/pde_ledger_v2/_scratch/decision16_amendment_directive.md` (the amendment directive + Codex/Grok review
  trail).
- Source map: `research/pde_ledger_v2/notes/stage007_pathA35_G0_source_map.md`.
