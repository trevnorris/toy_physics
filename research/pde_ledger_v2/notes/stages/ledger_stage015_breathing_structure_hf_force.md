# ledger_stage015_breathing_structure_hf_force

## Status

**Part II — Gravity. II-G2c (build-order 015).** Reshape of the **legacy-Hessian structure-recovery + Hellmann–Feynman-force
leg** of gate **`pathA_31`** (the scalar "breathing" of the frozen finite throat). Source top-line, verbatim:
**`BREATHING_CALIBRATED`** — a JOINT 3-stage verdict, and **this stage carries its legacy-structure + HF-force component (3/3)
and COMPLETES the joint.** ⭐ **pathA_31 splits 3-way** per the finalized Part-II split: **013** = the harmonic-lift profiles +
the `M_AB`/`K_AB` operator projection + the dynamical-EOM **LHS** (DONE, count 13); **014** = truncation consistency (the
combined-basis generalized eigenproblem, the `β_L0` sweep, the N-convergence) (DONE, count 14); **015 (this stage)** =
legacy-Hessian structure recovery + the Hellmann–Feynman force + the static-dynamic limit. 015 returns to the
**exact-symbolic / float-free** discipline of 013 (the INVERSE of 014, the one numeric stage) — the structure gate and the HF
force are `sp.integrate`/`FullSimplify` closed-form algebra with `expect_zero`-style asserts and exact rationals
(`L0 = 37/20`, never a decimal).

Stage 013 delivered the `(a,L)` collective closure `M_AB Q̈ + K_AB Q = F_A` — the **left-hand side** of the breathing
dynamics — by projecting the frozen wall operator onto two derived harmonic-lift profiles, and **deferred the RHS.** This
stage answers the two remaining questions: **(Q1)** does that operator-projected closure have the SAME structural form as the
legacy static energy Hessian `H_legacy` of `E_geom`? and **(Q2)** what is the EOM driving force `F_A^(HF)` (the deferred RHS),
and is it recovered by two genuinely-different constructions that agree? It then takes the **static-dynamic limit** and lands
the **completed** joint verdict.

- **CONSUMED (cited, dual-site integrity, NOT re-derived):** stage 013's **collective closure** — the harmonic-lift profiles
  `α_a, α_L` (for the HF source projection) and the operator-projected matrices `M_AB, K_AB` (for the structure-gate
  pattern-match) — and the **frozen wall packet** `{L0 = 37/20, β = 1, β·L0 = 37/20}`; and stage 014's **truncation
  certificate** (the validity window `β_L0 ∈ [0.1,3.0]` / `K_η/T_w ≲ 2.6`), cited in one line. ⚠ **013's `∫dw` operator
  projection / M/K closed-form derivation / profile BVP-solve and 014's generalized eig are NOT recomputed** — 015 declares
  the closed forms as cited literals, guards their citation-integrity, and runs its OWN 015 content (`H_legacy`, the structural
  comparison, the HF `4π∫ source_density·α_A dw` projection which 013 never did). ⚠ **`c_S` is NOT consumed** — 015 is a
  static structure + static force test; the matter-sector `c_s`/BdG `k⁴` is DEFERRED under `kξ≪1` (`phonon_limit_caveat`).
- **NEW SYMBOLS (015-owned):** the **legacy structure constants `{κ, χ, σ_a, σ_L}`** (the `E_geom` energy-Hessian
  parameterization the operator projection recovers) — ⚠ ABSENT from 013's `M_AB`/`K_AB` (013's free-symbol-name firewall
  proves the exclusion; they become LIVE HERE as the comparison basis `H_legacy`); and the **HF drive** `Vp0` (wall
  confinement potential radial slope) + `ell_c` (confinement length, now LIVE — INERT / tracked-not-counted at stage 011 where
  `δV_conf = 0`; here `δV_conf ≠ 0`, the HF force IS the confinement drive) + `ρ*` (frozen density, CONSUMED from
  stage005/011). In the force they appear ONLY as the ratio `Vp0/ell_c`.
- **EARNED (structure):** `H_legacy` own-built by `sp.hessian(E_geom)`; the genuine STRUCTURAL SIGNATURE MATCH of the cited
  013 `M_AB`/`K_AB` against `H_legacy` (M pos-definite by exact-identity Sylvester certificates; K symmetric; `K_aL < 0`;
  rank and zero-pattern match) → the `(a,L)` operator-projected closure is RECOVERED, not re-postulated; the
  **Hellmann–Feynman force `F_A^(HF)`** derived by TWO genuinely-different routes that agree only after simplification (the
  RHS 013 deferred); the static-dynamic limit.
- **CALIBRATED (values):** 015 adds **one new counted CALIB knob** — the breathing driving-force scale `Vp0/ell_c` (with
  `Vp0` + the now-live `ell_c` as manifestations of the one ratio-scale; `ρ*` consumed) — making 015 the **second
  calibration-adding Part-II stage** after 013. The legacy `{κ, χ, σ_a, σ_L}` are NOT counted afresh (they are the
  legacy-Hessian pattern basis, a re-parameterization of the already-calibrated `{μ_η, T_w, β}` closure). The verdict is the
  joint `BREATHING_CALIBRATED`, not `..._PASS`.

Ledger-local earned-label (NOT a source verdict token): `BREATHING_STRUCTURE_HF_FORCE_EARNED`. The joint verdict composes,
now **COMPLETE**, as `BREATHING_CALIBRATED = (013: harmonic-lift profiles + M_AB/K_AB by ∫dw operator projection + (a,L) EOM
LHS) ∧ (014: truncation consistency) ∧ (015: legacy-Hessian structure recovery + Hellmann–Feynman force + static-dynamic
limit)`, and the EOM RHS `F_A^(HF)` that 013 deferred is filled here.

## Purpose

Stage 013 reduced the ℓ=0 breathing to the collective closure `M_AB Q̈ + K_AB Q = F_A` and stopped at the LHS; stage 014
certified the 2-mode truncation is faithful. Two things remain to complete the calibrated breathing sector. First, the
operator-projected `M_AB`/`K_AB` should be recognizable as the same object a physicist would write from a static energy
functional — otherwise the "closure" is an unmoored algebraic artifact. Stage 015 builds the legacy static energy
`E_geom = ½κ(δ_L − χ δ_a)² + ½σ_a δ_a² + ½σ_L δ_L²`, forms its Hessian `H_legacy`, and shows the cited `M_AB`/`K_AB` carry the
SAME structural signature (a symmetric, full-rank, negative-off-diagonal, positive-definite ratio-plus-support form) — the
`(a,L)` closure is *recovered*, not re-postulated. Second, the dynamics need a right-hand side: the wall confinement drive.
Stage 015 supplies `F_A^(HF)` — and, because the pathA_31 v1 rejection's other locus was an **HF `x−x` tautology** ("both HF
routes" was one expression typed twice), it emits the two routes as GENUINELY DIFFERENT constructions — a distributed
source-projection `F_dist` and a Hellmann–Feynman parametric derivative `F_legacy` — that agree only after simplification, with
the raw unsimplified trees computed-distinct. The static-dynamic limit `Q̈→0 ⇒ K_AB Q = F_A` then ties the dynamical closure
back to the legacy static `∂E_geom/∂Q = 0`, completing the joint `BREATHING_CALIBRATED`.

## The derivation (both engines, own routes)

- **The legacy energy and its own-built Hessian (015-owned).** `E_geom = ½κ(δ_L − χ δ_a)² + ½σ_a δ_a² + ½σ_L δ_L²`;
  `H_legacy = ∂²E_geom/∂Q² = {aa: χ²κ+σ_a, aL: −χκ, LL: κ+σ_L}` (built via `sp.hessian` / `D[egeom,{{da,dL},2}]`), with
  `legacy_symmetric`, `legacy_offdiag_negative` (`−χκ < 0` for χ,κ>0), and `legacy_det_positive`
  (`det = κσ_a + κχ²σ_L + σ_aσ_L`). The legacy constants `{κ, χ, σ_a, σ_L}` enter here; a firewall echo confirms the cited
  `M_AB`/`K_AB` free symbols exclude them (they are `{L0, T_w, β, μ_η, r_AL}`).

- **The legacy-Hessian structure recovery (EARNED).** The structure gate reads the CITED `M_AB`/`K_AB` and checks a genuine
  STRUCTURAL SIGNATURE against `H_legacy`: **M positive-definite** by Sylvester — `M_aa > 0` and `det M > 0` for
  `B = β·L0 > 0`, each discharged by an **exact-identity certificate** (not `sp.ask`/`is_positive`, which return `None`/raise
  on the transcendentals, and NOT the source's `M_aa − m_aa_positive_form == 0` form-equality, which is `X≡X` once M is
  cited): `f1 = sinh B cosh B − B` has `f1(0)=0` and `f1' = cosh(2B) − 1 = 2 sinh²B` (a manifest square) ⇒ `f1 > 0`;
  `f2 = sinh²B − B² = (sinh B − B)(sinh B + B)` with `sinh B − B` monotone-positive (`g(0)=0`, `g' = cosh B − 1 =
  2 sinh²(B/2)`) ⇒ `f2 > 0`; **K symmetric**; **`K_aL < 0`** for B>0 (`K_aL ∝ −1/sinh B`); **rank(K) = rank(H_legacy) = 2**;
  **zero-pattern(K) = zero-pattern(H_legacy)** → `K_structure_ok`; `structure_from_computed_matrices = True` (firewall-tied);
  `full_matrix_fit = False` (the recovery is structural — `M_aa ≠ H_aa` entrywise — not a full numeric fit). The able-to-fail
  probes `M_aa → −M_aa` (non-posdef) and `K_aL → −K_aL` (sign-flip) flip the gate. **The operator-projected closure has the
  same structural form as the legacy energy Hessian.**

- **The Hellmann–Feynman force, both routes independently emitted (EARNED — the crux).** With the linearized wall confinement
  force density `source_density = ρ*·Vp0/ell_c` and the action measure `4π∫dw` (NOT μ_η-weighted):
  - **Route 1 (distributed source projection):** `F_dist_A = 4π∫₀^{L0} source_density·α_A dw`.
  - **Route 2 (Hellmann–Feynman parametric derivative):** `F_legacy_A = −4π∫₀^{L0} ρ*·∂V_conf/∂q_A|_{q=0} dw` with
    `V_conf = (Vp0/ell_c)(r − R0(w) − q_a α_a − q_L α_L)` — built by a genuine `sp.diff` / `D[vConf,q_A]`, carrying the
    `−∂/∂q` double-negative.
  The two agree after simplification: `hf_force_reduces = (F_dist_A − F_legacy_A == 0) = True`, with the closed form
  `F_dist_a = 4π Vp0 ρ*(cosh(L0β) − 1)/(β ell_c sinh(L0β))`. ⭐ The **anti-`x−x` guard** is computed:
  `unsimplified_routes_identical = (raw tree F_dist_A_uns == raw tree F_legacy_A_uns) = False` — the two UNSIMPLIFIED trees
  are literally different (`+4π∫(+ρ*Vp0α/ell_c)` vs `−4π∫(−ρ*Vp0α/ell_c)`), agreeing only after `.doit()`. This fills the EOM
  RHS `F_A^(HF)` that 013 deferred. The verdict rung `if not hf_force_reduces OR unsimplified_routes_identical:
  BREATHING_FAIL_HF_FORCE` places the anti-`x−x` in the OR — typing Route 2 as a copy of Route 1 trips it.

- **The static-dynamic limit (EARNED).** `Q̈ → 0` in the same `M_AB Q̈ + K_AB Q = F_A` system gives `K_AB Q = F_A`; after the
  legacy-pattern comparison this is the static `∂E_geom/∂Q = 0` limit (verified as the exact `Q̈→0` residual identity of the
  EOM, `(M Q̈ + K Q − F)|_{Q̈=0} = K Q − F`, not a fabricated separate static solve).

- **The HF profile guard (EARNED able-to-fail; the 014↔015 boundary).** 015 owns only the `F_a`/`hf_mismatch` slices of the
  source's counterfactual guard. The `degenerate` profile (`α_a ≡ 0`) gives `wrong_zero_F_a = 0` → `hf_mismatch =
  (wrong_zero_F_a − F_legacy_a ≠ 0) = True`; the `nontrivial` profile (`α_a ≡ 1`) gives `wrong_const_F_a = 4π L0 Vp0 ρ*/ell_c`
  → `hf_mismatch = True`. **Both wrong profiles are REJECTED by the HF** (their distributed force ≠ the frozen legacy force).
  ⚠ **The 014↔015 boundary:** 014 recorded that the constant profile PASSES the modal overlap (`overlap_passes = True`); 015
  shows it FAILS the HF. **The Hellmann–Feynman mismatch is the profile guard the modal overlap could not supply** — together
  with 013's `𝓛₀[α]=0` residual it certifies profile-correctness, which 014's overlap (truncation-consistency only) does not.

- **The completed landing.** Computed from the 015 rungs (the STRUCTURE rung `not (M_posdef ∧ K_symmetric ∧
  K_structure_ok)` → `BREATHING_FAIL_STRUCTURE`; the HF rung with the anti-`x−x` in the OR; the two `hf_mismatch` able-to-fail
  assertions), landing at the 015 component of `BREATHING_CALIBRATED` and printing the joint 3-stage composition as
  **COMPLETE** (013 done ∧ 014 done ∧ 015 computed here), with the EOM RHS `F_A^(HF)` filled and the three carried caveats.

## Consumed inputs

**Cited — no file reads; genuine DUAL-SITE citation-integrity with an anti-tautology guard.** stage 013's collective profiles
`α_a, α_L`, its operator-projected `M_AB, K_AB`, and the frozen packet `{L0 = 37/20, β = 1, β·L0 = 37/20}` are guarded by
independently-corruptible relations: **site A (branch anchor)** `β·L0 − 37/20 ≡ 0`, reading an independent cited datum;
**site B (profile-consumption)** the consumed profiles satisfy BOTH the harmonic residual `−α'' + β²α = 0` AND the collective
BC values `{α_a(0)=1, α_a(L0)=0, α_L(0)=0, α_L(L0)=r_AL}` (⚠ residual alone misses a kernel-preserving corruption — a rescale
leaves the residual 0 but breaks the BC — so residual ∧ BC are both able-to-fail; this VERIFIES the citation, it does NOT
re-derive the profiles); and **the consumed-M/K citation-integrity** covering every entry — the two **det-identities**
(`det M = 4π²μ_η²r_AL²(sinh²B − B²)/(β²sinh²B)`, `det K = 16π²T_w²β²r_AL²`) catch diagonal-value corruptions, and because the
determinant carries `M_aL²` it is BLIND to an off-diagonal SIGN flip, so explicit sign checks (`M_aL > 0` via `h = B − tanh B`,
`h' = tanh²B`; `K_aL < 0`) close that hole. Every one-value corruption (the `37/20` anchor, a kernel-preserving or non-kernel
profile, a cited M/K diagonal value or off-diagonal sign) fails BOTH engines. `c_S` is NOT consumed (matter sector deferred).
Stage 014's truncation certificate is a one-line narrative cite (no recomputation).

## Exports

- The **completed `(a,L)` breathing closure** (013's LHS + 015's RHS `F_A^(HF)`) + the legacy-structure recovery →
  **stages 022/023** (the ℓ=0 cross-ℓ map). This COMPLETES pathA_31.
- Register: **one new counted CALIB knob** — the breathing driving-force scale `Vp0/ell_c` (`Vp0` + the now-live `ell_c` its
  manifestations; `ρ*` consumed from stage005/011) — 015 is the **second calibration-adding Part-II stage** after 013. The
  legacy `{κ, χ, σ_a, σ_L}` are the legacy-Hessian PATTERN basis (a re-parameterization mapping to the calibrated `{μ_η, T_w,
  β}` closure; NOT counted afresh — that would double-count the closure). Structural edge **R32** (the legacy-Hessian
  structure recovery — the operator-projected M/K reproduce the static energy Hessian's structural signature; discharges
  nothing) + reduction-debt edge **R33** (the confinement-drive `Vp0/ell_c` reduction — a sibling of R30's wall-response debt;
  the same deferred nonlinear-throat interior that would derive `{μ_η, T_w, β}` would derive the confinement drive; PENDING).

## Verification

- **Reshape (blueprint §5) — bridge-severing + hybrid-`.wl` symbolic re-independence:** stripped the `.py`'s
  `expressions_for_mma` `sympy*` export dict, the scratch-WL numeric append, the scratch-YAML/results writers, the MMA-scratch
  re-read, the report/feed/summary writers, and `sympy_expression_digest`; and the `.wl`'s `Get[sympyExprFile]`, the whole
  `checks` `sympy*`-comparison assoc, the digest write-back, and the `Export`. The `.wl` KEEPS its native HF `Integrate` +
  Hessian `D` + structure `Det`/`MatrixRank` route but its asserts were **re-targeted to its OWN computed values**; it
  declares the cited profiles + M/K closed forms as native literals (NOT via `DSolveValue`, NOT by re-`Integrate`), builds its
  OWN `H_legacy` and HF force (Route 2 a genuine `D[vConf,qa]`, not a copy of Route 1), and emits its OWN
  `unsimplifiedRoutesIdentical` raw-tree comparison. Both engines standalone, print-only, **zero file I/O**, SYMBOLIC /
  float-free (exact rationals — `L0 = 37/20`; `expect_zero`/`expect_bool` asserts; no scipy/numpy/floats/tolerances), ledger
  idioms. **Clean 3-way cut:** no 013 content (harmonic profiles / `M_AB`/`K_AB` `∫dw` projection / EOM LHS — cited) and no
  014 content (generalized eig / overlaps / β_L0 sweep / N-convergence — cited); `c_S` is not a live symbol.
- **Dual-engine agreement (transcript-level):** the two engines independently compute the structure verdict (M pos-def / K
  structure OK / probes flip), the two HF routes' simplified forms + `hf_force_reduces = True` +
  `unsimplified_routes_identical = False`, and the two `hf_mismatch = True` guard slices — SymPy vs Wolfram
  `Integrate`/`D`/`Det` (the `.wl` `Tanh[βL0/2]` form equals the `.py` `(cosh L0β − 1)/sinh L0β` form via the half-angle
  identity, each asserted against its own expected form); the scripts do NOT cross-read (the arbiter compares transcripts).
  The `.wl` carries an arity self-check + unevaluated-leakage guard (the stage-007 silent-skip lesson).
- **Dual-engine:** SymPy **95 PASS / 0 FAIL** · Mathematica **102 PASS / 0 FAIL**, both exit 0, CWD-independent (repo root +
  foreign CWD); runner transcripts under `scripts/output/` + `mathematica/output/`.
- **Directive review — Codex → Grok → Codex bookend (blueprint §6):** Codex `DIRECTIVE_NEEDS_REVISION` → 1 BLOCKING folded
  (the consumed-M/K det-identity was blind to the off-diagonal sign flip `M_aL → −M_aL` — the determinant carries `M_aL²`;
  added explicit `M_aL > 0` / `K_aL < 0` sign checks) → Codex confirm `DIRECTIVE_CLEAN`; then a **Grok-4.5 compute-verify
  pass** caught 1 BLOCKING (the Sylvester positivity `M_aa > 0` is NOT float-free dischargeable — `sp.ask(Q.positive)` returns
  `None`, `bool(expr>0)` raises — so, left unspecified, it would silently regress to the vacuous `M_aa − m_aa_positive_form ==
  0` form-equality; required explicit exact-identity certificates + a form-equality prohibition) + 1 nit, both folded → Codex
  final confirm `DIRECTIVE_CLEAN`. Grok independently confirmed the Hessian/det forms, the HF raw-tree difference, the
  Sylvester certificates, both `hf_mismatch` values, and the det-blind-to-sign fact.
- **Tri-review (fresh agents):** arbiter re-run via the runners (both engines, both CWDs, reproduced); **`FIDELITY_CLEAN`**
  (an independent read-and-reason audit hand-verified every value — all seven positivity-certificate identities and the four
  entry-to-certificate residuals reduce to 0, the two HF integrals, the Hessian/det, both `hf_mismatch` values; same claims,
  no drop/soften/upgrade; clean consume-from-013/014 boundary; the HF two-route genuineness and the certificate positivity
  real); **`ADVERSARIAL_ISSUES`** — the three central burdens proven GENUINE by live copy-mutations (the HF two-route
  anti-`x−x` — typing Route 2 = Route 1 trips `BREATHING_FAIL_HF_FORCE`; the structure-gate exact-identity certificate — no
  disguised form-equality; the consumed-013 dual-site incl. the det-blind off-diagonal; the `.wl` a genuine independent
  engine), but **1 BLOCKING falsification-strength defect: tooth 4 was a vacuous `x−x` stamp in both engines**
  (`expect_fail(compact(F_legacy_a − F_legacy_a) != 0)` — a self-compare, empirically blind to a mutation of its object;
  the real profile-guard teeth backstop it), plus 3 lower-severity unable-to-fail literal/scaffolding flags.
- **Remediation (make-genuine ×4, de-count ×2, both engines):** tooth 4 rewired to read the REAL wrong-profile distributed
  forces (`expect_nonzero(wrong_{zero,const}_F_a − F_legacy_a)`, a mutation making a wrong force equal `F_legacy_a` fires it);
  the static-dynamic flag made a genuine `Q̈→0` EOM-residual identity; `full_matrix_fit` made `expect_nonzero(M_aa − H_aa)`;
  `structure_from_computed_matrices` tied to the free-symbol firewall; the trivially-True typed-twice `expect_bool` and the
  two vacuous "bypassing probe" meta-asserts de-counted (the genuine verdict + entry-corruption teeth retained). A
  **fresh-agent re-verify** confirmed **`REVERIFY_CLEAN`** via the coupling meta-test — every one of the 6 remediated
  constructs fires under a mutation of its named object in both engines; the three central burdens + the physics/verdict
  untouched; no new vacuous tooth. New tallies 95/102 (net −4 each from the honest de-counts).
- **Teeth (all fire, mutations on copies):** (1a) type Route 2 = Route 1 → `unsimplified_routes_identical` True → the HF rung
  fires `BREATHING_FAIL_HF_FORCE`; (2) corrupt one HF route → `hf_force_reduces` False → the rung fires; (3) structure
  genuineness — each positivity traces to a certificate identity (NOT a form-equality); corrupt cited `M_aa` sign →
  `M_posdef` fails, corrupt `K_aL` sign → `K_structure_ok` fails, mutate `H_legacy` (drop χ) → the zero-pattern match fails;
  the consumed-M/K citation ablations — a cited diagonal entry → the det-identity fires, an `M_aL` sign flip → the `M_aL > 0`
  check fires (the det is blind), a `K_aL` sign flip → `K_aL < 0` fires; (4) replace a wrong profile's distributed force with
  `F_legacy_a` → `hf_mismatch` False → the guard no longer rejects; (5) the structure probes flip the gate; (6) corrupt
  `E_geom` (drop the `σ_a` support or the `χ` cross term) → `H_legacy` entry/`legacy_det_positive` fires; (7) any one-value
  corruption of the consumed 013 packet (the `37/20` anchor; a kernel-preserving profile → BC leg; a non-kernel profile →
  residual leg; a cited M/K entry) fails both engines; (8) the `.wl` arity self-check catches a def/call mismatch.

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_31_scalar_breathing_{sympy.py,.wl}` (015 slice = `build_structure_gate`
  L217–311; `E_geom`/`H_legacy` L681–687; the Hellmann–Feynman force L689–703; the `hf_force`/`static_dynamic` payloads
  L960–984, `unsimplified_routes_identical = hstr(F_dist_a_uns) == hstr(F_legacy_a_uns)` L977; the `F_a`/`hf_mismatch` guard
  slices L716–717/L824–838; the 015 verdict rungs L1053–1060; sources unchanged);
  `software/stage1_solver/reports/pathA_31_scalar_breathing.md` (`## Structure` :60–66, `## Hellmann-Feynman Force` :68–73,
  `## Static-Dynamic Limit` :75–77; `## Counterfactual Guard` :79–82 shared). SIBLING (cited, not recomputed):
  `## Operator, BCs, Inner Product` / `## Profiles And Projection` / `## Dynamical EOM` / `## Dimensional Check` :3–36, :87–98
  = stage 013; `## Truncation Consistency` :37–58 = stage 014.
- Reshape directive + review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage015_*` (directive; Codex→Grok→Codex
  design-review + confirm logs; execute/remediation/re-verify logs; fidelity + adversarial + re-verify prompts). The
  directive design-review used the Codex → Grok → Codex bookend (blueprint §6): the two BLOCKING folds were the det-blind
  off-diagonal sign hole (Codex) and the float-free positivity-certificate requirement (Grok, compute-verified).
- Running-start source map: `research/pde_ledger_v2/notes/stage015_pathA31_structure_hf_source_map.md`. Split row:
  `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 015). ⭐ **COMPLETES the joint `BREATHING_CALIBRATED`** =
  (013: M/K projection + (a,L) LHS) ∧ (014: truncation consistency) ∧ (015: legacy-Hessian structure recovery +
  Hellmann–Feynman force RHS + static-dynamic limit).
