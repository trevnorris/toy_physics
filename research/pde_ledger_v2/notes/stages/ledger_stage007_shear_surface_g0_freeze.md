# ledger_stage007_shear_surface_g0_freeze

## ⚠ Decision-16 amendment (2026-07-21) — the P-retirement layer (READ FIRST)

Decision 16 (user, 2026-07-13; `software/stage1_solver/decisions/16_retire_brane_polar_field.md`) retires the brane
polar field `P` and its couplings — `L_pol` (the T0 polar-OP block), `L_Pu` (the parity-even P–u coupling), the
constant `λ_Pu`, and structural postulates 3/4/5 — from the frozen medium action. Every `P` payoff failed
independently (charge `NATIVE_P_NO_EMERGENT_GAUSS`; light `pathA_35 FAIL_COUPLE_STRESS_NOGO`, and `pathA_36` derives
both photons with **no** `P` sector; the polar-vector wall self-localization falsified, `pathA_24` T1), and the frozen
massless-`P`+`λ_Pu` baseline is structurally unstable — a long-wavelength helical Lifshitz mode
`det = k²(2 A_P μ_R k − λ_Pu²) < 0` for any `λ_Pu ≠ 0` (`INSTABILITY_CONFIRMED_STRUCTURAL`, the pre-registered
`pathA_35` G0.3 `FAIL_NOT_BOUNDED_BELOW` exposure firing; retained finding, commit `a5c079eb`).

**Route A — the freeze is retained; retirement is a computed layer, not a re-freeze.** This stage now carries TWO
tiers (both computed in both engines):

- **HISTORICAL (freeze-as-run, immutable).** The G0 freeze is **unchanged on disk**; both engines still fence-parse
  and SHA-256-verify the frozen `d9520d3819c3` block (T0 `8fa41ac51e88` byte-embedded) against the untouched reports,
  and still compute the freeze's honest cost *as run* — `SECOND_MEDIUM_DRIFT_AT_FREEZE(11)` and **DOF = 8**. History
  is **not** falsified downward: a dedicated tooth fires if anyone tries to compute the historical drift as 7 or the
  historical DOF as 4. The freeze is the anchor the instability finding is *about* (Decision 16 §Scope); re-introducing
  orientational order would need a **new T0 freeze** with its own gauntlet — not un-retiring this.
- **OPERATIVE (post-Decision-16, live).** The medium's live Part-I action is the retired subset, computed as an exact
  **action-summand set partition** of the historical summands:
  `operative {S_GNLS, gL_Mac, gL_uw}` ⊎ `retired {L_pol, gL_Pu}` == `historical {S_GNLS, L_pol, gL_Mac, gL_Pu, gL_uw}`
  (token `POST_D16_ACTION{S_GNLS,gL_Mac,gL_uw}_OF(d9520d3819c3)`; a symbolic partition over the immutable hash anchor,
  **no byte-substring surgery**). Operative **DOF = 4** (the removed `Pⁱ` block = tangent 3 + radial 1 = 4 takes
  8 → 4), and operative **`POST_D16_DRIFT(7)`** derived as `historical ∖ {λ_Pu, postulates 3/4/5}` = 3 constants
  `{ρ_br, μ_R, Ω_w}` + 1 function `g_ℓ` + 3 postulates `{1, 2, 6}`. Postulate (1)'s annotation softens to "intrinsic
  wall normal" (no longer an axis conceded *for* `L_Pu`); its membership and count are unchanged.

Retirement is **able-to-fail**: action-summand partition teeth (mutate a survivor / admit a retired term / drop a
survivor), operative-DOF teeth (one reinjected `Pⁱ` mode → 5, full block → 8), operative-drift teeth (leave `λ_Pu` → 8,
leave a retired postulate → 8, a same-cardinality `Ω_w`↔`λ_Pu` swap → set-partition failure at unchanged count), the
rebased survivor drift-enumeration teeth, and the historical-integrity teeth. **Untouched:** `c_γ² = μ_R/ρ_br`,
`L_Mac`, `L_uw`, `g_ℓ`, `{ρ_br, μ_R, Ω_w, ℓ_g}`, and the R10 Route-A reduction; Part III/light consumes `L_Mac` as-is
(no `P` sector). The **`χ_B = |P_∥|²` route (c) is NOT dead** — it stays a named, high-risk, Part-VII-adjacent future
gate needing a new T0 freeze (obsolete-as-carried, not foreclosed). Provenance stamp:
`DECISION16_PROVENANCE retired={L_pol,L_Pu,λ_Pu,postulates_3/4/5}; reason=P_RETIRED_ALL_PAYOFFS_FAILED_PLUS_LIFSHITZ_INSTABILITY`.

Ledger earned-label for the amended stage: `G0_FREEZE_FIDELITY_PLUS_POST_D16_LAYER_VERIFIED`. Dual-engine after the
amendment: **SymPy 142 PASS / Mathematica 140 PASS**, both exit 0, CWD-independent.

## Status

> ⚠ **Amended by Decision 16 (2026-07-21) — see the P-retirement layer above.** The text below describes the
> HISTORICAL freeze-as-run (the immutable tier); the OPERATIVE live action drops `L_pol`/`L_Pu`/`λ_Pu` + postulates
> 3/4/5 (DOF 8 → 4, drift 11 → 7).

**Part I — The medium. I-4 (build-order 007; the LAST Part-I stage).** Reshape of the existing dual-engine gate
`pathA_35` **G0** (the shear-surface brane/light-sector **constitutive freeze**). Source headline verdicts, carried
verbatim:

- **`T0_SHEAR_FROZEN(d9520d3819c3)`** — the frozen brane/light action block, SHA-256
  `d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3`, with the T0 polar-OP block (short hash
  `8fa41ac51e88`) kept byte-for-byte and byte-embedded inside it.
- **`SECOND_MEDIUM_DRIFT_AT_FREEZE(11)`** — the honest drift ledger: 11 independent new inputs, **computed in both
  engines from an enumerated table** (4 constants + 1 function + 6 structural postulates), never a typed literal.

**EARNED** (what this stage genuinely computes): freeze-fidelity (fence-parsed byte-exact hash re-verification of
both freeze blocks + the T0 embedding), the flat-brane linear **DOF = 8** (rank-computed, able-to-fail), and the
dimensional firewall over every frozen term. **POSTULATED/CALIBRATED** (the honest landing, first-class): the 11
freeze inputs. The freeze freezes **terms, not gate answers** (anti-impose); this stage computes **no gate verdict**
(Gate-L excluded, below) and does **not** earn light (stage003's job — stage003 derives the 2 transverse photons at
`c_γ² = μ_R/ρ_br` from this stage's frozen `L_Mac`).

Ledger-local earned-label (NOT a source verdict token): `G0_FREEZE_FIDELITY_DOF_DIM_VERIFIED`.

## Purpose

G0 is where the program **postulated the full shear-surface brane structure up front and paid for it in public**:
every new constant, the one profile function, and every structural concession counted into a single drift number,
frozen *before* any gate was run (the pathA_35 directive's binding rule — "postulating an *ingredient* is allowed;
postulating an *outcome* is not"; a late ingredient = `AD_HOC_RESCUE` → fresh G0; a clean all-pass is suspicious).
The rebuilt ledger needs this stage as the **formal home** of the freeze token that Part III consumes, and as the
Part-I closing statement of what the medium's brane sector *costs*.

## The frozen action (canonical block content; hash-guarded)

The historical freeze declares **five action summands** `{S_GNLS, L_pol, gL_Mac, gL_Pu, gL_uw}`. Decision 16 retires
`{L_pol, gL_Pu}` (the `P`-dependent complement); the **operative** action is `{S_GNLS, gL_Mac, gL_uw}` (each summand
labeled ⟨operative⟩ / ⟨retired⟩ below). The byte-level SHA-256 anchor is unchanged — the partition is symbolic.

**⟨historical, kept 0 new⟩** the GNLS parent action `S_GNLS` (I-1/I-2) — ⟨operative, SURVIVES⟩.

**⟨RETIRED by Decision 16⟩** the T0 polar-OP action `L_pol`, byte-for-byte from `pathA_24_T0_freeze.md` (T0 block hash
`8fa41ac51e88…`, byte-embedded in the G0 block; still hash-verified as a *historical* record, dim-audited as a
retired-historical term — it *was* dimensionally homogeneous, gone by decision not defect):

```
L_pol = ½ m ρ a² (D_t^v Pⁱ)² − ½ m ρ c_s² a² (∂_j Pⁱ)² − ¼ m ρ c_s² (PⁱPⁱ−1)²          ⟨RETIRED⟩
```

**Frozen brane/light blocks** (`S_brane = ∫ dt d⁴X g_ℓ(w) [L_Mac + L_Pu + L_uw]`):

```
L_Mac = ½ ρ_br (∂_t uᵃ)² − ½ μ_R Ω_uᵃ Ω_uᵃ          (MacCullagh rotational shear; Ω_u = ∇_∥×u)   ⟨operative⟩
L_Pu  = −λ_Pu ϖ_a Ω_uᵃ,   ϖ_a := (ŵ×P_∥)_a           (parity-EVEN P–u coupling)                    ⟨RETIRED⟩
L_uw  = ½ ρ_br (∂_t u_w)² − ½ ρ_br Ω_w² u_w²          (the u_w gap)                                  ⟨operative⟩
g_ℓ(w) = exp(−(w/ℓ_g)²)/(√π ℓ_g),   ∫ g_ℓ dw = 1      (fixed Gaussian profile, one width ℓ_g)        ⟨operative⟩
```

⚠ Preserved verbatim from the freeze report (now a **retired-historical parity record**): the `L_Pu` operator
"re-admits the ε-contracted/chiral class excluded by T0" and REQUIRED the conceded axis `ŵ` — a structural-postulate
cost, not a free choice (the direct `P_∥·Ω_u` is parity-ODD, excluded). With `L_Pu` retired, the operative `ŵ` is just
the **intrinsic wall normal** (postulate 1, annotation softened; no longer a concession *for* the P–u operator).

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

Exponent-triple `{L,T,M}` audit of every frozen term, faithful to the source gate's full surface (fidelity-leg
coverage-diffed, nothing dropped): kept GNLS (`c_s²=5Kρ⁴/m`, `U=(K/4)ρ⁵`, quantum pressure, `mρvv`), kept T0
`L_pol` (each term `M L⁻² T⁻²`, the 4D action density; `Pⁱ` dimensionless), the kept T0 couple-stress coefficients
(labeled kept, 0 new), profile + measure (with `∫g_ℓ dw = 1` genuinely integrated in both engines), the three brane
blocks (each term `M L⁻¹ T⁻²`, the 3D surface density; `g_ℓ·[brane]` restores the 4D density), the action measures,
the full projected traction `T_na` including `T_wa = mρ v_w v_a`, and the G0.5 linearization mode frequencies.

**Post-Decision-16 firewall split:** the **operative live** surface audits only the survivors — kept GNLS, `L_Mac`,
`L_uw`, `g_ℓ`, `T_na`, `O_u`, `c_γ² = μ_R/ρ_br`, and the surviving mode `ω_uw,bare² = Ω_w²`. The `P`-sector terms
(`L_pol`, `L_Pu`, `λ_Pu`, and the mode frequencies **`ω_P²` and `ω_radial²` — both retired P-modes**, freeze report
:251 / :252–253) stay dim-audited but only as **retired-historical** freeze-as-run records, never live-checked as
survivors.

Targets: `Pⁱ` dimensionless (retired) · `uᵃ, u_w` = `L` · `ρ_br` = `M L⁻³` · `μ_R` = `M L⁻¹ T⁻²` · `λ_Pu` =
`M L⁻¹ T⁻²` (retired) · `Ω_w` = `T⁻¹` · `g_ℓ` = `L⁻¹` · `ℓ_g` = `L`. Teeth (both fired in the source gate and re-fire
here): `drop_m_from_T_wa`, `MacCullagh_without_curl`.

**Notational firewall (register edge R22):** this stage's `μ_R` (3D brane modulus, `M L⁻¹ T⁻²`) and stage006's
`μ_R⁽⁴⁾` (4D shear-stiffness density, `M L⁻² T⁻²`) are asserted **dimensionally distinct** (a computed
exponent-triple inequality, with its own tooth), and the R17 projection `μ_R = ∫ χ_B μ_R⁽⁴⁾ dw` is asserted
dim-consistent (`[μ_R⁽⁴⁾]·L = [μ_R]`) — status **PENDING**. Two symbols, one pending reduction; never conflated.

## Flat-brane linear DOF = 8 (EARNED, rank-computed)

The total is **computed from rank/nullity bookkeeping** in both engines and asserted equal to the frozen 8 — never
typed. Breakdown (report G0.5): in-plane `uᵃ` curl-transverse **2** + kinetic−curl **1**; T0 `Pⁱ` tangent **3** +
radial **1**; `u_w` **1**; C5 `φ` **0**. The two engines use genuinely different linear algebra: the SymPy keeps
the source gate's k-aligned curl-rank machinery; the Mathematica builds a generic `k=(1,2,3)` outer-product
transverse projector and ranks that. The absent `φ` (no longitudinal-constraint analog — structural postulate 6) is
printed plainly: it is the fact the later C5 crux attacks (stage006's θ-as-φ no-go; the material-state pivot).
Teeth: `drop_u_w_gap_term`, `drop_P_soft_spin_radial_term`, `zero_u_longitudinal_component` — each recomputes
total 7 ≠ 8. The new-field content at G0 (`uᵃ` 3 + `u_w` 1 = 4 DOF) is computed separately from an enumerated field
list and kept **out** of the 11-input drift count.

**Operative post-Decision-16 DOF = 4** (rank-computed on the retired field set, `Pⁱ` removed): `uᵃ` curl-transverse 2
+ kinetic−curl 1 + `u_w` 1 + `φ` 0 = 4; the removed `Pⁱ` block (tangent 3 + radial 1 = 4) is exactly what takes
historical 8 → operative 4. Operative teeth: reinject one `Pⁱ` mode → 5, the full `Pⁱ` block → 8, and the survivor
ablations `drop_u_w_gap_term` / `zero_u_longitudinal_component` → 3 (each ≠ 4). The `drop_P_soft_spin_radial_term`
ablation stays HISTORICAL-only (no operative analogue — `P` is gone); a historical-integrity tooth blocks falsifying
DOF=8 down to 4.

## The drift ledger: the 11, computed (POSTULATED/CALIBRATED — the honest landing)

The source gate's `ledger()` returned the subcounts as hardcoded literals; the reshape replaces this with the
stage006 computed-tally pattern. Both engines build an explicit enumeration table and **compute** the subcounts and
`n = 4+1+6 = 11` by counting categories; the verdict string is **built from the computed n** and asserted equal to
the frozen token.

- **4 constants:** `ρ_br` [M L⁻³] surface inertia · `μ_R` [M L⁻¹T⁻²] MacCullagh modulus · `λ_Pu` [M L⁻¹T⁻²]
  parity-repaired P–u coupling · `Ω_w` [T⁻¹] bare u_w gap scale. (Table dims are asserted equal to the
  independently-verified firewall dims — the table cannot silently misinform.)
- **1 function:** `g_ℓ(w)` — fixed Gaussian shape, ONE width knob `ℓ_g`; admitted on locality/minimality grounds
  only (target-blind G0.2); no free-form profile.
- **6 structural postulates** (verbatim): (1) imposed `ŵ` axis + `w=0` surface (conceded-wall); (2) `uᵃ`
  same-medium surface collective, tangentially free-slip (`u̇ᵃ ≠ vᵃ`); (3) T0 `Pⁱ` reused as the Cosserat
  micro-rotation reservoir (0 new DOF); (4) baseline `Pⁱ` spin-wave status = `massless` (alternates
  `gapped`/`slaved-rigid` named-inactive); (5) the `ŵ`-dependent parity-EVEN P–u operator (re-admits the
  ε-contracted/chiral class excluded by T0; requires the conceded `ŵ`); (6) no C5 `φ` analog / no longitudinal
  constraint.

T0 couple-stress coefficients: asserted 0 new (kept). Teeth: drop-one-entry → n=10 fires; miscategorize → subcount
asserts fire; inject `ρ_B0` → the anti-absorption guard fires; corrupt n before verdict assembly → the
verdict-string equality fires; corrupt a table dim → the dim-consistency assert fires. Adversarial ablation
confirmed the checks are anchored to in-engine derivations and the frozen token (dual-corruption of both engines
fails both) — no cross-engine comparison anywhere; the old `--compare` payload mirror is dead.

**Operative post-Decision-16 drift `POST_D16_DRIFT(7)`** — derived (not merely recounted) as the exact set partition
`historical ∖ {λ_Pu, postulates 3/4/5}` = 3 constants `{ρ_br, μ_R, Ω_w}` + 1 function `g_ℓ` + 3 survivor postulates
`{1, 2, 6}` (postulate 1's annotation softened to "intrinsic wall normal"). Operative teeth: leave `λ_Pu` live → n=8;
leave any retired postulate → n=8; a **same-cardinality** `Ω_w`↔`λ_Pu` swap keeps n=7 but **fails the set-partition
assert**; the rebased survivor teeth (drop a survivor → n=6, miscategorize `Ω_w`, corrupt `ρ_br` table-dim); inject
`ρ_B0` → operative anti-absorption; corrupt n → token equality. The historical 11 stays computed; a
historical-integrity tooth blocks falsifying it down to 7.

## The 2026-07-04 erratum (carried first-class) + anti-absorption

**The 11 STANDS — no `ρ_br` overcount.** An earlier GLM-based claim that `ρ_br` duplicated pathA_25's derived
`varrho_br[ρ]` was superseded by the pathA_41 framing adjudication: `varrho_br` belongs to the **CLOSED**
density-smectic candidate (`FAIL_NOT_CODIM1`), `OUT_OF_ACTIVE_NG5` — a different structure from this active
shear-surface brane. This `ρ_br`/`μ_R` is genuine postulated shear-surface inertia/modulus with a
**registered-pending pathA_40 Route-A reduction** (`ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT`; register
edge R10; corroboration token `NO_OVERCOUNT_ROUTE_A_PENDING`). **The honest cross-sector drift (per pathA_41) is
`{ρ_B0, χ_c, C_hu}` — a Part-VI item, never absorbed into the historical 11 _or_ the operative 7** (guarded in-engine
by the anti-absorption assert on both tables). The historical 11 STANDS as the freeze-as-run cost; the operative 7 is
the Decision-16 retirement, **not** an overcount correction.

## Supersession (both facts, exactly)

1. The stage006 χ_B order-field wall **superseded the fixed-shape `g_ℓ(w)` profile as the *material-state*
   closure** (the wall is a postulated field, not a fixed profile).
2. The G0 freeze **REMAINS the light-sector *constitutive* freeze** — the MacCullagh/P–u/u_w action; stage003
   consumes the frozen `L_Mac` as-is.

Neither artifact retro-invalidates the other (register edge R21 — a scope split, not a reduction; `ℓ_g` stays a
counted knob of the constitutive freeze).

## Gate-L exclusion

Scope carried verbatim: **"G0 freeze only; no Gate L verdict computed"**, plus the classification guard (counts
only; no gate verdict, no boundedness, no gauge, no leak claim). The exposure-name strings (the gates that "remain
able to fire") are printed as prose provenance only — nothing is computed from them; nothing is imported from the
Gate-L artifacts (`pathA_35_gateL_sympy.py`/`.wl`, `reports/pathA_35_gateL_light.md` = `FAIL_COUPLE_STRESS_NOGO`,
a Part-III matter).

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
