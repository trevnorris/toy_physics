# II-G2c (ledger_stage015) source map — pathA_31 legacy-Hessian structure recovery + Hellmann–Feynman force (COMPLETES BREATHING_CALIBRATED)

> Running-start prep captured 2026-07-08 (post stage014, in the /clear prep) so a fresh session can author the reshape
> directive without re-discovery. **All line refs below VERIFIED against the current sources 2026-07-08** (`_sympy.py` 1359,
> `.wl` 428 — but note the `.wl` was cleaned in stage014's build; the pathA_31 SOURCE `.wl` under `software/stage1_solver/`
> is unchanged, 428 lines).
> Companion: `part2_gravity_atomic_split.md` (rows 013/014/015 + the pathA_31 trip-ups); the SIBLING maps
> `stage013_pathA31_breathing_source_map.md` (013 = M/K projection, built) + `stage014_pathA31_truncation_source_map.md`
> (014 = truncation, built). Build-order id **015**, Part II. **pathA_31 SPLITS 3-way: 013 (M/K projection, DONE count 13)
> + 014 (truncation consistency, DONE count 14) + 015 (THIS: legacy-Hessian recovery + Hellmann–Feynman force).**
> **⭐ 015 is the COMPLETING leg** — it COMPLETES the joint `BREATHING_CALIBRATED = (013) ∧ (014) ∧ (015)` (the 010-completes-009,
> 012-completes-011 pattern). Source top-line: **`BREATHING_CALIBRATED`**.

## ⭐ The THREE headline differences from 014 (READ FIRST)
1. **⚠ 015 is SYMBOLIC / float-free again — 014 was numeric.** 015 returns to the exact-symbolic discipline of 013 (the
   structure gate + the HF force are `sp.integrate`/`FullSimplify` closed-form algebra, no `eigh`/`quad`). So the contract
   is `float-free` (exact rationals; `L0=37/20` never `1.85`), `expect_zero`-style asserts — the INVERSE of 014's
   float-bearing exception. (014 was the one numeric stage; 015 is back to symbolic.)
2. **⭐⭐ THE pathA_31 v1 REJECTION scar's OTHER locus lives HERE — the HF `x−x` two-route tautology.** pathA_31 v1 was
   tri-review-REJECTED for (a) a gamed truncation threshold (014's scar, absent there) AND **(b) an HF `x−x` tautology** —
   "both HF routes" was ONE expression typed twice. **015 must emit the TWO HF routes as GENUINELY DIFFERENT
   constructions** (a distributed source-projection vs a Hellmann–Feynman parametric derivative) that agree only after
   simplification; the `unsimplified_routes_identical` guard (verdict L1059–1060) FAILS if the two UNSIMPLIFIED expression
   trees are literally identical. This is the reshape's central able-to-fail obligation (contrast 013's M/K-projection
   genuineness, 014's anti-gaming threshold).
3. **⭐ 015 is the COMPLETING stage — print the FULL joint verdict as now-earned.** Unlike 013/014 (which printed
   `BREATHING_CALIBRATED` as a partial component with siblings PENDING), 015 lands all three components → print the joint
   `BREATHING_CALIBRATED = (013) ∧ (014) ∧ (015)` as COMPLETE (the 010/012 completing pattern). It also emits the EOM RHS
   `F_A^(HF)` that 013 deferred (013 emitted the LHS only).

## §1 The 015 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)
- **`build_structure_gate(M, K, H_legacy, L0, beta, mu_eta, T_w, r_AL, kappa, chi, sigma_a, sigma_L)` L217–311 — the
  legacy-Hessian structure recovery.** Pattern-matches the operator-projected `M_AB`/`K_AB` (013's, CONSUMED) against the
  legacy static Hessian `H_legacy`: `m_posdef` (leading minor + det vs the positive closed forms, L244–251); `k_structure`
  (K symmetric + `K_aL<0` + `K_rank` matches `legacy_rank` + `K_zero_pattern` matches `legacy_zero_pattern` →
  `K_structure_ok`, L253–276); the legacy checks `legacy_symmetric`/`legacy_offdiag_negative`/`legacy_det_positive`
  L238–242; the structure_counterfactual probes `M_aa→−M_aa` (non-posdef) + `K_aL→−K_aL` (sign-flip) L280–283/L298–310.
- **`E_geom` + `H_legacy` + the `build_structure_gate(...)` call L681–687:** `E_geom = ½κ(δ_L−χδ_a)² + ½σ_a δ_a² + ½σ_L δ_L²`
  L681–685; `H_legacy = sp.hessian(E_geom, (δ_a, δ_L))` L686; the gate call L687. ⚠ **The legacy structure constants
  `{κ, χ, σ_a, σ_L}` enter HERE (015's territory)** — they are the legacy energy-Hessian parameterization the operator
  projection recovers; they must be ABSENT from 013's `M_AB`/`K_AB` (013's free-symbol-name ancestry firewall).
- **The Hellmann–Feynman force (BOTH routes) L689–703 — the RHS `F_A^(HF)` 013 deferred:**
  - `source_density = ρ*·Vp0/ell_c` L689 (the linearized wall confinement force density).
  - **Route 1 — distributed source projection:** `F_dist_A = 4π∫₀^{L0} source_density·α_A dw` (`F_dist_a_uns`/`F_dist_L_uns`
    L690–691 → `.doit()` L692–693).
  - **Route 2 — Hellmann–Feynman parametric derivative:** `V_conf_linear = (Vp0/ell_c)(r − R_param)`,
    `R_param = R0(w) + qa·α_a + qL·α_L` L695–697; `F_legacy_A = −4π∫ ρ*·∂V_conf/∂q_A |_{q=0} dw` (`F_legacy_a_uns`/`_L_uns`
    L698–699 → `.doit()` L700–701).
  - **The two-route agreement:** `hf_a_ok = (F_dist_a − F_legacy_a) == 0`, `hf_L_ok` L702–703. ⚠ **The routes are genuinely
    DIFFERENT unsimplified trees** — route 2 carries the `−∂/∂q` double-negative of the Hellmann–Feynman derivative
    (report :72 shows `−4π∫(−Vp0 ρ* …)`) while route 1 is the direct projection (report :71 `4π∫(Vp0 ρ* …)`) — they agree
    only after simplification (this is the anti-`x−x` genuineness).
- **The HF slices of the counterfactual guard (015-owned) L824–826 + L834–838:** in the `degenerate` block,
  `F_a_dist`/`F_a_legacy_frozen`/`hf_mismatch` L824–827; in the `nontrivial` block, `F_a_dist`/`F_a_legacy_frozen`/
  `hf_mismatch`/`fails` L834–838 (⚠ `nontrivial["fails"]` is driven by `hf_mismatch` — the constant profile is REJECTED by
  the HF, NOT the overlap; this is the 014↔015 boundary: 014 recorded constant PASSES the overlap, 015 REJECTS it via HF).
  The `wrong_zero_F_a`/`wrong_const_F_a` counterfactual HF values L716–717.
- **The 015 verdict rungs L1053–1060:** structure rung L1053–1058 (`not (M_posdef and K_symmetric and K_structure_ok)` →
  `BREATHING_FAIL_STRUCTURE`); **HF rung L1059–1060** (`not hf_force_reduces OR unsimplified_routes_identical` →
  `BREATHING_FAIL_HF_FORCE` — ⭐ the `unsimplified_routes_identical` anti-`x−x` guard).
- **The static-dynamic limit (report :75–77):** `Q̈→0` in `M_AB Q̈ + K_AB Q = F_A^(HF)` → `K_AB Q = F_A` → after the legacy
  pattern comparison, the static `∂E_geom/∂Q = 0` limit. (Prose/payload assembly — locate the exact payload key at build.)
- **⭐ CLEAN CUT — 015 owns the structure gate + the HF force + the static-dynamic limit + the HF guard slices ONLY. Do NOT
  recompute (013/014 territory):**
  - 013: the harmonic profiles + `M_AB`/`K_AB` `∫dw` projection L627–679 + the dim block L326–444 + the native degenerate
    `M_det` slice. **013 is DONE — CITE its `M_AB`/`K_AB` + profiles for the structure pattern-match + the HF projection.**
  - 014: `galerkin_overlap`/`build_truncation`/`_baseline_functions` L447–607 + the module constants L50–55 + the calls
    L760–762 + the `wrong_o_k`/`overlap_passes` guard slices. **014 is DONE — CITE the truncation certificate.**

## §1b The `.wl` 015 slice (VERIFIED) — L27 + L46–50 + L59–60 + L62–99
- `src` L27 (HF source density); `forceDistA`/`forceDistL` L46–47 (Route 1); `vConf` L48 + `forceLegacyA`/`forceLegacyL`
  L49–50 (Route 2); `wrongZeroFA` L59 / `wrongConstFA` L60 (HF counterfactual slices); `egeom` L62 + `legacyH` L63
  (legacy Hessian); the structure gate L65–99 (`structureMPosdef`, `structureKStructureOk`, the M/K-vs-legacy rank/zero-pattern
  match L80–86, the probes L87–99). KEEP this as the native independent route (native `Integrate` HF + `D[egeom,{{da,dL},2}]`
  Hessian + `Det`/`MatrixRank` structure).
- **⚠ The bridge to sever (same as 013/014):** `Get[sympyExprFile]` L20, the `sympy*`-comparison `checks` assoc (L171–207,
  which cross-checks `wrongZeroFA`/`wrongConstFA`/`forceDist`/`forceLegacy` vs `sympy*` at L186–189), the digest write-back
  L316, `Export` L425. SEVER: the `.wl` asserts against its OWN native HF/structure values. (013 `.wl` = M/K + dim block
  L26/L29–44/L101–169; 014 `.wl` = numeric Galerkin L209–306 — EXCLUDE both.)

## §2 The 015 claim-set (derive + assert; report quotes)
- **Legacy-Hessian structure recovery (EARNED; report `## Structure` :60–66).** `E_geom = ½κ(δ_L−χδ_a)² + ½σ_a δ_a² +
  ½σ_L δ_L²` (:62); `H_legacy = ∂²E_geom/∂Q²` = `{aa: χ²κ+σ_a, aL: −χκ, LL: κ+σ_L}` (:63). Pattern check (:64): the
  operator-projected `M`/`K` (013's) are a "symmetric full-rank ratio-plus-support Hessian with negative off-diagonal for
  χ>0" — `M_posdef=True`, `K_structure_ok=True`, `K_offdiag_negative=True`; `structure_from_computed_matrices=True` (:65).
  The able-to-fail probes `M_aa→−M_aa` + `K_aL→−K_aL` (:66) flip the gate. ⭐ **This certifies the operator-projected
  closure has the SAME structural form as the legacy energy Hessian** — the (a,L) closure is recovered, not re-postulated.
- **The Hellmann–Feynman force, BOTH routes independently emitted (EARNED — the CRUX; report `## Hellmann-Feynman Force`
  :68–73).** Measure `4π∫dw` (not μ_η-weighted, :70). Route 1 (distributed, :71) vs Route 2 (legacy-parametric HF, :72);
  Simplified differences `{F_a: 0, F_L: 0}` (:73) — the two GENUINELY-DIFFERENT routes agree. ⚠ Emit `unsimplified_routes_identical
  = (F_dist_a_uns == F_legacy_a_uns as raw trees)` = **False** (the anti-`x−x` guard) + `hf_force_reduces = (hf_a_ok and
  hf_L_ok)` = True. **This fills the EOM RHS `F_A^(HF)` that 013 deferred.**
- **The static-dynamic limit (EARNED; report :75–77).** `Q̈→0` → `K_AB Q = F_A` = the static `∂E_geom/∂Q=0` limit (ties the
  dynamical closure back to the legacy static energy).
- **The HF guard slices (EARNED able-to-fail; report :81).** The `degenerate`/`nontrivial` `hf_mismatch` — ⚠ the
  `constant_one` profile is REJECTED here via `hf_mismatch=True` (the HF is the REAL profile guard 014's overlap could not
  supply; the 014↔015 boundary — carry it: 014 showed constant PASSES the overlap, 015 shows it FAILS the HF).
- **The 015-scoped landing + the COMPLETED joint verdict.** Land at the 015 component AND print the joint as COMPLETE:
  `BREATHING_CALIBRATED = (013: harmonic profiles + M/K operator-projection + (a,L) EOM LHS) ∧ (014: truncation consistency)
  ∧ (015: legacy-Hessian structure recovery + Hellmann–Feynman force, computed here) — COMPLETE`. `CALIBRATED ⇐ {μ_η,T_w,β}
  (013) [+ any 015 HF-force knob, §6]`.

## §3 Reshape cost (the bridge to sever) — hybrid `.wl`, SYMBOLIC block
Same scratch-YAML payload-mirror family as 013 (NOT the numeric 014). Strip the `.py`'s scratch-YAML/`_sympy_exprs.wl`
export + the MMA re-read + the report/feed/summary writers (`main` L1344–1355); and the `.wl`'s `Get`/`sympy*`-checks
(L171–207, the HF/structure cross-checks at L186–207)/digest/`Export`. **KEEP the `.wl`'s native HF `Integrate` + Hessian
`D` + structure `Det`/`MatrixRank` route, re-targeted to its OWN values.** ⚠ **Route-2 genuineness in the `.wl` too:**
`forceLegacyA` (L49, the `−4π∫ρ*(D[vConf,qa]/.q→0)`) must be a genuine parametric-derivative construction, NOT the same
expression as `forceDistA` (L46) — the `.wl`'s own anti-`x−x`. **Zero file I/O.** Arity discipline (standing).

## §4 Consumed / exported
- **Consumes (cite, dual-site integrity, don't re-derive):**
  1. **013's `M_AB`/`K_AB`** — for the structure-gate pattern-match (the gate compares 013's M/K against `H_legacy`). ⚠
     Same normalization note as 014: 013's closed-form M/K are the symbolic objects; cite them, do NOT re-run the `∫dw`
     projection. 2. **013's profiles `α_a`/`α_L`** — for the HF source projection `4π∫ source_density·α_A dw`. 3. **013's
     EOM LHS** — 015 fills the RHS `F_A^(HF)`. 4. The frozen packet + `ρ*` (frozen density, consumed from stage005/011).
  ⚠ **`c_S` NOT consumed** (matter-sector deferred, `kξ≪1`) — as 013/014.
- **Exports:** the COMPLETED `(a,L)` breathing closure (LHS 013 + RHS 015) + the legacy-structure recovery → the ℓ=0
  cross-ℓ map (022/023). Completes pathA_31.

## §5 Teeth candidates (015-specific)
1. **⭐ The HF two-route anti-`x−x` genuineness (the CRUX):** `unsimplified_routes_identical` must be a COMPUTED comparison
   of the two RAW route trees (False) — ablation: force the two routes to be the literal same expression → the guard fires
   → `BREATHING_FAIL_HF_FORCE`. AND `hf_force_reduces` (the two routes agree after simplification) — ablation: corrupt one
   route (e.g. drop a sign / a profile factor) → `hf_a_ok` fails.
2. **Structure-gate probes** (`M_aa→−M_aa` non-posdef; `K_aL→−K_aL` sign-flip) → the gate flips → `BREATHING_FAIL_STRUCTURE`;
   the legacy rank/zero-pattern match must be a genuine comparison (a mutated legacy Hessian breaks it).
3. **The HF `hf_mismatch` guard slices** (degenerate + constant profiles) — the constant profile's `hf_mismatch=True`
   genuinely REJECTS it (the 014↔015 boundary: constant passed the overlap, fails the HF). Ablation: a profile that
   matched the HF would not be rejected.
4. **Consumed-013 dual-site** (M/K + profiles cited; a corruption fires both engines) + `.wl` arity scan.
⚠ **NOT 015 (do not pull in):** the harmonic-residual / M/K `∫dw` projection / forbidden_fit_flags (013); the
generalized-eig / overlaps / floor / β_L0 sweep / N-convergence (014).

## §6 Register expectation — ⭐ THE KEY 015 QUESTION (may add a calibration knob; CONFIRM)
Unlike 014's clean zero, **015 introduces the HF driving force**, which brings new parameters — determine at registration
+ Codex-verify:
- **`Vp0` (wall confinement potential radial slope) + `ell_c` (confinement length)** — the HF force density is `ρ*·Vp0/ell_c`.
  `ell_c` was INERT + tracked-not-counted at stage011 (`δV_conf=0`); HERE `δV_conf≠0` (the HF force IS the confinement
  drive), so `ell_c` becomes LIVE, and `Vp0` is new. ⭐ **Whether `Vp0` (or the force-density group `ρ*Vp0/ell_c`) is a NEW
  counted CALIB knob** (the breathing driving-force scale — 015 would be the SECOND calibration-adding Part-II stage after
  013) **or is absorbed/construction/consumed** MUST be determined against the scripts + Codex-verified. Do NOT pre-decide;
  the honest read is likely a new CALIB force-scale, but confirm whether it cancels in the observable or is set by the
  frozen packet. `ρ*` is consumed (stage005/011).
- **The legacy structure constants `{κ, χ, σ_a, σ_L}`** — these are the legacy energy-Hessian PARAMETERIZATION the operator
  projection RECOVERS (a consistency mapping to `{μ_η, T_w, β}`), NOT (likely) independent new knobs — they are the
  comparison basis of the structure gate. CONFIRM they map onto 013's counted set (a re-parameterization), not counted
  afresh (that would double-count the closure).
- Likely a structural edge recording the legacy-Hessian-recovery provenance (the operator projection reproduces the static
  energy Hessian) + possibly a reduction-debt edge if `Vp0` is a new CALIB (the same R30 nonlinear-throat reduction that
  would derive the wall response would also derive the confinement drive). CONFIRM at registration.

## Verdict tokens + honest scope
015 carries the **structure-recovery + HF-force component (3/3) of `BREATHING_CALIBRATED`** and COMPLETES the joint verdict
(013 ∧ 014 ∧ 015). EARNED-structure = the operator-projected closure recovers the legacy energy-Hessian form + the HF force
is derived by TWO independent routes that agree (the RHS `F_A^(HF)`). CALIBRATED = the wall constants `{μ_η,T_w,β}` (013)
[+ likely the HF drive `Vp0`/force-scale, §6]. Caveats: the HF is the profile guard the overlap could not supply (the
014↔515 boundary — constant profile fails the HF, passes the overlap); `c_S` deferred (`kξ≪1`). ⭐ The v1 HF `x−x` rig is
015's locus — both routes GENUINELY DIFFERENT, `unsimplified_routes_identical` computed.

## Process (unchanged, calibrated — the UPDATED per-stage pipeline)
Author the II-G2c reshape directive (§1 the clean 015 slice + KEEP-015-only + §2 faithful cover + §3 bridge-strip incl. the
hybrid-`.wl` HF/structure re-independence + §5 the HF anti-`x−x` teeth + §6 the KEY `Vp0` register question) → **Codex xhigh
design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless compute-verify pass** (assess + independently
verify each catch — Grok is noisy) → fold → **Codex confirm-pass on the folds** → **pre-exec USER GATE** → Codex builds the
two scripts (`--sandbox danger-full-access`, background, stdin, xhigh) → dual-engine both exit 0 (repo root + foreign CWD) →
arbiter re-run → full tri-review (fidelity + adversarial-with-**per-tooth ablation**, ⭐ hunt the HF two-route genuineness +
any vacuous able-to-fail) → remediate → fresh-agent re-verify → bump counts 14→15 → parameter register (⭐ resolve the `Vp0`
question) + Codex-verify → note/card/`\input{stages/stage_015}` + registration → rebuild PDF → commit + docs/memory sync.
Orchestrator authors notes/cards/LaTeX/registration; Codex codes. Target stem: `ledger_stage015_breathing_structure_hf_force`
(confirm slug at authoring).
