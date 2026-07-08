# II-G2a (ledger_stage013) source map — pathA_31 harmonic profiles + M/K projection

> Running-start prep captured 2026-07-08 (post stage012, in the /clear prep) so a fresh session can author the reshape
> directive without re-discovery. Verify line refs against the sources before finalizing (recon agent line counts:
> report 112, `_sympy.py` 1359, `.wl` 428, directive 289).
> Companion: `part2_gravity_atomic_split.md` (rows 013/014/015 + the Cluster-A reshape-cost map + the pathA_31 trip-ups).
> Build-order id **013**, Part II. **pathA_31 SPLITS 3-way: 013 (II-G2a harmonic profiles + M/K projection) + 014
> (II-G2b truncation consistency) + 015 (II-G2c legacy-Hessian + Hellmann–Feynman force).**
> **Source top-line: `BREATHING_CALIBRATED`** — the joint 3-stage verdict; **013 carries the M/K-projection component**
> (the (a,L) collective closure's LHS), 014 the truncation-consistency component, 015 the structure + HF component.

## ⭐ The two headline differences from 011/012 (READ FIRST)
1. **A 3-WAY split, not 2-way.** The clean cuts are TIGHT — 013 = the harmonic-lift profiles + the `M_AB`/`K_AB`
   operator-projection + the dynamical-EOM LHS assembly. It must **NOT** recompute the generalized-eigenproblem /
   β-sweep / N-convergence overlaps (**014**) or the legacy-Hessian pattern-match / Hellmann–Feynman force (**015**).
   See §1 for the exact cut lines.
2. **⭐ 013 is the FIRST Part-II CALIBRATED stage — it ADDS calibration knobs.** Unlike stages 008–012 (all "zero new
   counted knobs"), pathA_31 downgrades to `BREATHING_CALIBRATED` precisely because the wall constants `{μ_η, T_w, K_η}`
   are **calibration inputs** (`stiffness_provenance` = "calibration"; `.py` L988–992). This is a genuine register
   departure — see §6. The EARNED content is the *structure* (`M_AB`/`K_AB` from real ∫dw operator projection, profiles
   derived as harmonic lifts); the *values* `μ_η,T_w,K_η` are calibrated.

## The `.wl` is a HYBRID (unlike 012's clean transfer-matrix route)
The pathA_31 `.wl` is a genuine independent symbolic route (native `DSolveValue` BVP + `Integrate` + `NIntegrate`
generalized `Eigensystem`) **that nonetheless imports the SymPy expr export and cross-checks against it** (`Get[sympyExprFile]`
L19; `checks` assoc differences `FullSimplify[... − sympy... == 0]` L171–207; numeric overlaps vs `sympy*` L291–311). So
the reshape is HEAVIER than 012's: severing the bridge means the `.wl` must derive its own 013 `M_AB`/`K_AB` and assert
against its OWN values (not against imported `sympy*` symbols). This is the Cluster-A reshape-cost note's "it must derive
its own inputs instead of Importing the `.py`'s expressions" applied to a genuinely-hybrid engine.

## File inventory (all under `software/stage1_solver/`)
- **Report:** `reports/pathA_31_scalar_breathing.md` — 013 content: `## Operator, BCs, Inner Product` :3–10, `## Profiles
  And Projection` :11–30, `## Dynamical EOM` :31–36 (LHS only; RHS `F_A^(HF)` = 015), `## Dimensional Check` :87–98
  (013-anchored — walks `M_AB`/`K_AB`), `## Reduction Certificate` :99–102 (shared provenance). **014:** `## Truncation
  Consistency` :37–59. **015:** `## Structure` :60–67, `## Hellmann-Feynman Force` :68–74, `## Static-Dynamic Limit`
  :75–78. Shared: `## Counterfactual Guard` :79–82 (degenerate `M_det→0` = 013; `wrong_o_k` = 014; `hf_mismatch` = 015),
  `## Engine Agreement` :83–86 (dropped bridge), `## Files` :103–112.
- **`.py`:** `tools/pathA_31_scalar_breathing_sympy.py`. **`.wl`:** `tools/pathA_31_scalar_breathing.wl`.
  Directive `directives/pathA_31_scalar_breathing.md`; `pathA_31_results.yaml`.

## §1 The 013 slice (`.py` line ranges) — the CLEAN CUTS
- **013 CORE = inside `symbolic_engine()` (def L610), lines 630–679:**
  - **harmonic-lift profiles** L630–644: `general = C1·sinh(β w) + C2·cosh(β w)` (L632); `coeff_a` from
    `α_a(0)=1, α_a(L0)=0` (L633–637) → `α_a = sinh(L0·β − β w)/sinh(L0·β)` (L643); `coeff_L` from
    `α_L(0)=0, α_L(L0)=r_AL` (L638–642) → `α_L = r_AL·sinh(β w)/sinh(L0·β)` (L644).
  - **harmonic residual assertion** L646–649: `L0_alpha_a = (−T_w·α'' + K_η·α)/μ_η`; **raises if ≠0** (a 013 tooth —
    proves α genuinely solves `L₀α=0`).
  - **M integrands** L651–655 (`μ_η·α_A·α_B`); **K integrands** L656–660 (`T_w·α_A'·α_B' + K_η·α_A·α_B`).
  - **the ∫dw projection** L661–666: `M_aa..M_LL`, `K_aa..K_LL = 4π·integrate(<integrand>, (w,0,L0))`.
  - `M`,`K` matrices + `M_det`,`K_det` L667–670; `build_dimensional_check(...)` call L671–679.
  - Payload 013 sub-dicts (in the return L878–1006): `modal_operator` L886–898, `profiles` L899–912,
    `projection_integrals` L913–920, `dynamical_EOM` L922–931 (LHS; note `K_AB_provenance="operator_projection_not_static_Hessian"`
    L930), `stiffness_provenance`/`calibration_inputs` L988–993, `dimensional_check` key L994.
  - **Dim helpers (013-anchored, shared machinery):** `dim_of`/`build_dimensional_check` L106–215, L314–444; `compact()`
    L58–59.
- **⭐ CLEAN CUT — 013 symbolic content ENDS at L679. Do NOT recompute (all inside/after `symbolic_engine`):**
  - L681–687 `E_geom` + `H_legacy` + `build_structure_gate(...)` → **015**;
  - L689–703 `source_density` + distributed HF `F_dist_a/L` + legacy-parametric HF `F_legacy_a/L` + `hf_*_ok` → **015**;
  - L760–809 `build_truncation(...)` calls + the numeric-WL append → **014**.
- **Standalone functions that are 014/015 territory (013 must NOT touch):** `build_structure_gate` L217–311 (**015**);
  `_baseline_functions` L447–473, `galerkin_overlap` L476–557, `build_truncation` L560–607 (**014**).

## §2 The 013 claim-set (derive + assert; report quotes)
- **The ℓ=0 wall operator + inner product (CONSUMED/POSTULATED; report :3–10):** `L₀ = μ_η⁻¹(−∂_w(T_w ∂_w ·) + K_η ·)`
  on `w∈[0,L0]`; axisymmetric breathing (the `T_Ω` angular term drops at `ℓ(ℓ+1)=0`); μ_η-weighted inner product
  `⟨f,g⟩_μ = 4π∫₀^{L0} μ_η f g dw` (`.py` L896). Print the ℓ=0 restriction + `Y00=1/(2√π)` (cert L843–847).
- **Harmonic-lift profiles (EARNED; report :11–30):** `α_a`,`α_L` solved as harmonic lifts of `L₀` under the collective
  BCs (`profile_provenance="derived"`, cert L861–862). The residual assertion `L₀α_a = L₀α_L = 0` (L646–649) is the
  able-to-fail check that they genuinely annihilate the operator.
- **`M_AB`/`K_AB` by real ∫dw operator projection (EARNED — the 013 crux):** `M_AB = 4π∫μ_η α_A α_B dw`,
  `K_AB = 4π∫[T_w α_A' α_B' + K_η α_A α_B] dw` — COMPUTED by `sp.integrate`, NOT typed from the legacy Hessian
  (`K_AB_provenance="operator_projection_not_static_Hessian"`; `forbidden_fit_flags` all `False`: `K_from_static_hessian`,
  `M_or_K_typed_from_legacy_values`, `full_matrix_fit`, L869–875). This is the (a,L) collective-closure LHS.
- **The dynamical-EOM LHS (EARNED; report :31–36):** assemble `M_AB Q̈ + K_AB Q = F_A^(HF)` — 013 emits the LHS
  (`M`,`K`) + the EOM structure; the RHS `F_A^(HF)` is filled by **015** (do NOT compute the HF force here).
- **The 013 dimensional leg + corrupt-`[T_w]` probe (EARNED; report :87–98):** walk `[M_AB]`,`[K_AB]` from the sourced
  dims (`μ_η:(-1,1,0)`, `T_w:(1,1,-2)`, `β:(-1,0,0)`, `L0:(1,0,0)`, `r_AL:0`, L336–342); the corrupt-`[T_w]` probe
  (adds one power of L → `BREATHING_FAIL_DIMENSIONAL`, L358–374) + the two-verdict self-ablation (L1089–1104,
  `fail_suppressed`). This is the 013 analog of 011's corrupt-`[K]` probe.
- **The 013-scoped verdict + joint composition:** compute the 013 component (the M/K-projection certificate) and print
  the joint `BREATHING_CALIBRATED = (013: M/K operator-projection + (a,L) closure) ∧ (014: truncation consistency) ∧
  (015: legacy-structure + HF force)`. Do NOT type `BREATHING_CALIBRATED` as 013-earned alone (the pattern from
  009/010 + 011/012).

## §3 Reshape cost (the bridge to sever) — HEAVIER than 012 (hybrid `.wl`)
Scratch-yaml payload-mirror (Cluster-A variant) — SAME family as 011/012 but the `.wl` is a hybrid. Strip:
- **`.py` → MMA expr export `_scratch/pathA_31_sympy_exprs.wl`:** `SYM_EXPR_WL` def L46; `expressions_for_mma` +
  `digest_mapping`→`expr_digest`/`sympyExpressionDigest` L719–752; `SYM_EXPR_WL.write_text(...)` L754–758; the numeric
  Galerkin append L808–809.
- **`.py` → own scratch YAML:** `SYM_YAML` L45; `yaml_write(SYM_YAML, payload)` L1346.
- **`.wl` READS + cross-checks (SEVER — make the `.wl` assert against its OWN 013 values):** `sympyExprFile` L16;
  `Get[sympyExprFile]` L19–20; the `checks` assoc `FullSimplify[... − sympy... == 0]` L171–207; the numeric-overlap
  comparison vs `sympy*` L291–311; the digest write-back L316.
- **`.wl` → own scratch YAML:** `yamlOut` L17; YAML build L313–423; `Export[yamlOut,...]` L425.
- **`.py` re-reads the MMA scratch:** `MMA_YAML` L47; `engine_status()` L1009–1040 (`digest_matches`/`engine_agreement`).
- **Results-YAML + report + feed writers:** `result_payload()` L1107–1155 + `yaml_write(RESULTS_YAML,...)` L1351;
  `write_report()` L1158–1273; `write_feed_note()` L1276–1304; `print_summary()` L1307–1341; `main()` L1344–1355.
- **⭐ The `.wl` 013 KEEP (make it fully independent):** `α_A`/`α_L` via `DSolveValue` L29–34; `mass*`/`stiff*` via
  `Integrate` L36–42; dets L43–44; the dimensional block + corrupt-`[T_w]` probe L101–169. Re-target its asserts from
  `sympy*`-differencing to its OWN native `M`/`K`. **Zero file I/O.** (015 `.wl` = L46–99, 014 `.wl` = L209–306 — EXCLUDE.)
- **Arity discipline (standing, stage007–012 lesson):** def/call arity scan + unevaluated-leakage transcript scan.

## §4 Consumed / exported
- **Consumes (cite, dual-site integrity, don't re-derive — §2f pattern):** the **frozen throat packet** from Gate-1
  (pathA_30 / **ledger_stage011+012**): the domain `[0,L0]`, the frozen const-coeff wall packet `L0=37/20, T_w=1,
  K_η=1, β=1` (`BETA_L0_FROM_R0=37/20` L55, `beta_from_R0` L588–595), and the ℓ=0 restriction. **`L0=37/20` is the
  branch-determinable `L/a=37/20`** (see `[[project-ansatz-value-catalog]]` route (b)) — cite it, don't refit. ⚠ `c_S`
  does NOT appear as a live 013 symbol — the matter-sector `c_s`/BdG `k⁴` is explicitly DEFERRED (`phonon_limit_caveat`
  L868), so this stage does NOT consume the pathA_30 speed. (Dual-site integrity applies to the consumed geometry packet;
  design the two independent sites at authoring — e.g. `β` from `√(K_η/T_w)` vs from the `R0` taper.)
- **Exports (split-note flow):** the `(a,L)` collective closure + `M_AB`/`K_AB` → **014** (truncation consistency
  consumes the M/K matrices for the generalized eig) + **022/023** (the ℓ=0 map). 013's share = the harmonic-lift
  profiles, the `M_AB`/`K_AB` operator projection, and the dynamical-EOM LHS.

## §5 Teeth candidates (013-specific)
Keep/assign to 013: **harmonic-lift residual assertion** (`L₀α_a = L₀α_L = 0`, L648–649 — corrupt a profile → residual
≠0 → raise); **M/K-produced-by-∫dw** (`forbidden_fit_flags` all False + `operator_projection_not_static_Hessian` — a
mutant that types `K_AB` from the legacy Hessian must be caught); **degenerate counterfactual `M_det→0`** (`wrong_zero_M`
L707–712, 817); **corrupt-`[T_w]` dimensional probe** (L358–374 → `BREATHING_FAIL_DIMENSIONAL` + two-verdict
self-ablation L1089–1104); **consumed-geometry corruption** (dual-site, fires both engines). ⚠ **NOT 013 (do not pull in):**
the HF `x−x` two-route enforcement (`unsimplified_routes_identical` L977 → `BREATHING_FAIL_HF_FORCE`) = **015**; the
overlap-floor `o_k≥0.9` + β-sweep anti-gaming = **014**. PLUS the `.wl` def/call arity scan + unevaluated-leakage scan.

## §6 Register expectation ⭐ (the departure from 008–012)
**013 ADDS calibration knobs** — the first Part-II stage to do so. The wall constants `{μ_η, T_w, K_η}` are
`stiffness_provenance="calibration"` (L988–992; cert `input_provenance` L863–866). At registration, determine the
**counted** set:
- `μ_η` (wall inertia, `[μ_η]=(-1,1,0)` = M L⁻¹), `T_w` (wall tension, `[T_w]=(1,1,-2)` = M L T⁻²) — likely **COUNTED
  calibration knobs** (they are tuned, not derived).
- `K_η` (wall stiffness) = `T_w·β²` (L628) — **tied to `T_w` and `β`**, so likely a MANIFESTATION of `T_w` (not
  independently counted); `β = √(K_η/T_w)` is fixed by the geometry (`β L0 = 37/20` via `beta_from_R0`), so `β` is
  ACTION-geometry (tracked, not counted, like stage011's `L0`).
- So the likely irreducible add = **{μ_η, T_w}** (2 counted), with `K_η` a manifestation and `β` geometry — but CONFIRM
  against the scripts at authoring (this is the register's first Part-II calibration entry; get the provenance classes
  right per `[[feedback-parameter-register-every-stage]]` — an `IMPOSED`/calibration knob dressed as `DERIVED` would
  falsely shrink the count). No new reduction edge expected (013 is a calibrated closure, not a reduction); a structural
  edge may record the operator-projection provenance if useful.
- `r_AL` (the BC ratio) is dimensionless (`r_AL:0`); `ℓ_c`, legacy `κ,χ,σ_a,σ_L` appear but are 015/legacy-comparison
  symbols — check whether they enter 013 (they should NOT).

## Verdict tokens + honest scope
013 carries the M/K-projection component of **`BREATHING_CALIBRATED`**: the harmonic-lift profiles + the `M_AB`/`K_AB`
operator ∫dw projection + the (a,L) dynamical-EOM LHS are EARNED-structure; the wall constants `{μ_η,T_w,K_η}` are
CALIBRATION inputs (→ `CALIBRATED`, not `PASS`). Caveats to carry (mostly 014/015's, but note them): the modal overlap
does NOT guard profile-correctness (014); clean 2-mode truncation only for order-unity wall stiffness `K_η/T_w≲2.6`
(014); the phonon-limit / BdG `k⁴` is DEFERRED (`kξ≪1`). CITED: the frozen geometry packet `L0=37/20` (branch-determinable),
the ℓ=0 restriction. Nothing here re-derives `c_S` (matter sector deferred).

## Process (unchanged, calibrated)
Author reshape directive (§1 the 3-way clean cut + KEEP-013-only + §2 faithful cover + §3 bridge-strip incl. the
hybrid-`.wl` re-independence + §5 013 teeth + §6 the FIRST-calibration register note) → Codex xhigh design-review →
fold to `DIRECTIVE_CLEAN` (no GLM on Parts I–VI; 013 is repackaging-of-calibrated — but ⚠ per the clarified GLM policy,
if the calibration-knob accounting reads as "might be new physics" reconsider; it is not — it repackages pathA_31's
already-earned calibrated closure) → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`,
background, stdin, absolute paths, xhigh) → dual-engine exit 0 (repo root + foreign CWD) → arbiter re-run → tri-review
(fidelity + adversarial-with-ablation, incl. the M/K-operator-projection genuineness + the 3-way cut integrity + arity
scan + tally spot-check) → remediate → fresh-agent re-verify → registration 12→13 + parameter register (the FIRST
calibration knobs — Codex-verify the counted set) → note/card/`\input{stages/stage_013}` → PDF → commit + docs/memory
sync. Target stem: `ledger_stage013_breathing_harmonic_mk_projection` (confirm slug at authoring).
