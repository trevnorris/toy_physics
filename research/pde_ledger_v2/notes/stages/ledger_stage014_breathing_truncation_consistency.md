# ledger_stage014_breathing_truncation_consistency

## Status

**Part II — Gravity. II-G2b (build-order 014).** Reshape of the **truncation-consistency leg** of gate **`pathA_31`** (the
scalar "breathing" of the frozen finite throat). Source top-line, verbatim: **`BREATHING_CALIBRATED`** — a JOINT 3-stage
verdict, and **this stage carries its truncation-consistency component (2/3).** ⭐ **pathA_31 splits 3-way** per the
finalized Part-II split: **013** = the harmonic-lift profiles + the `M_AB`/`K_AB` operator projection + the dynamical-EOM
LHS (DONE, count 13); **014 (this stage)** = truncation consistency (the combined-basis generalized eigenproblem, the
`β_L0` sweep, the N-convergence); **015** = legacy-Hessian structure recovery + the Hellmann–Feynman force. This is the
**FIRST NUMERIC / float-bearing stage in the rebuild** — stages 004–013 were exact/symbolic/float-free; 014 runs
`scipy.linalg.eigh` + `scipy.integrate.quad` (`.py`) and `Eigensystem` + `NIntegrate` (`.wl`).

Stage 013 reduced the ℓ=0 breathing to a **2-mode collective-coordinate closure** `Q = (δ_a, δ_L)` and projected the
frozen wall operator onto the two collective harmonic-lift profiles `α_a, α_L` to build `M_AB, K_AB`. **This stage tests
that the 2-mode truncation is faithful, not an artifact:** it enlarges the basis to `B = {α_a, α_L, g_1..g_N}` (the two
collective modes PLUS `N` sine "lane" modes `g_n = sin((n−½)πw/L0)`), solves the combined-basis generalized eigenproblem
`K v = ω² M v`, and asks whether the two lowest generalized modes are captured by the collective `span{α_a, α_L}` to a
predeclared overlap floor, across a `β_L0` sweep and as `N→16`.

- **CONSUMED (cited, dual-site integrity, NOT re-derived):** stage 013's **collective closure** — the harmonic-lift
  profiles `α_a, α_L` used as the collective basis of `B` — and the **frozen wall packet** `{L0 = 37/20, β = 1, β·L0 =
  37/20}` (the branch-determinable `L/a`). ⚠ **013's `∫dw` operator projection, `M_AB`/`K_AB` closed forms, and EOM LHS
  are NOT recomputed** — 014 builds its OWN numeric Galerkin Gram over the full basis; the citation-integrity is on the
  consumed profiles + packet. ⚠ **`c_S` is NOT consumed** — 014 is a static-spectrum truncation test; the matter-sector
  `c_s`/BdG `k⁴` is DEFERRED under `kξ≪1` (`phonon_limit_caveat`).
- **POSTULATED (labeled):** the `g`-lane basis `g_n = sin((n−½)πw/L0)` (BCs `g(0)=0, T_w g'(L0)=0` — the half-shifted
  Neumann-at-cap ladder, the same one stage 012's DtN produced) is the IMPOSED enlarging basis. The truncation controls —
  `FLOOR = 0.9` (`EPS_TRUNC = 0.1`), `N_FINAL = 16`, `N_CONVERGENCE = [4,8,12,16]`, and the `BETA_L0_SWEEP` grid — are
  predeclared **method/tolerance** parameters (tracked, not counted — convergence/solver knobs, not physics DOFs).
- **EARNED (structure):** the combined-basis generalized eigenproblem genuinely solved (`eigh`/`Eigensystem`); the modal
  overlaps `o_1, o_2` (mass-Gram projection of the two lowest generalized modes onto `span{α_a, α_L}`); the floor-gated
  truncation certificate; the COMPUTED `β_L0` validity window (with a real upper edge); the N-convergence.
- **CALIBRATED (values):** 014 **adds no new counted knobs** — the underlying wall constants `{μ_η, T_w, K_η, β}` are
  stage-013's calibration inputs. The verdict inherits the joint `BREATHING_CALIBRATED`, not `..._PASS`.

Ledger-local earned-label (NOT a source verdict token): `BREATHING_TRUNCATION_CONSISTENT_EARNED`. The joint verdict
composes as `BREATHING_CALIBRATED = (013: harmonic-lift profiles + M_AB/K_AB by ∫dw operator projection + (a,L) EOM LHS)
∧ (014: truncation consistency — combined-basis generalized eig / β_L0 window / N-convergence, computed here) ∧ (015:
legacy-structure + Hellmann–Feynman force)`.

## Purpose

Stage 013 delivered the `(a,L)` collective closure `M_AB Q̈ + K_AB Q` — the left-hand side of the breathing dynamics —
by projecting the frozen wall operator onto two derived harmonic-lift profiles. But a 2-mode truncation of an
infinite-dimensional wall field is only meaningful if those two modes genuinely capture the low end of the spectrum. This
stage supplies that certificate. It embeds the two collective modes in a much larger Galerkin basis (adding `N` sine lane
modes), solves the full generalized eigenproblem `K v = ω² M v`, and measures how much of each of the two lowest
generalized eigenvectors lives in `span{α_a, α_L}` — the modal overlaps `o_1, o_2`. When both clear a predeclared floor,
the 2-mode truncation is clean. The earned deliverable is a **honest validity window**: the truncation is clean only for
order-unity wall stiffness (`K_η/T_w ≲ 2.6`); sharp walls genuinely fail and would need more modes. ⚠ **The honest scope
limit — carried, not softened:** the modal overlap certifies *truncation-consistency*, NOT *profile-correctness* (a
constant, wrong profile still clears the overlap floor); the profile guard is 013's `𝓛₀[α]=0` residual and 015's
Hellmann–Feynman mismatch. The **pathA_31 v1 rejection scar was a gamed truncation threshold**, so the reshape's central
burden is that the floor and window are genuinely COMPUTED (not typed to pass) and able-to-fail. The source pair computed
this through a scratch-YAML bridge with a hybrid `.wl` (a genuine `NIntegrate`+`Eigensystem` Galerkin route that
nonetheless imported and float-diffed the `.py`'s numeric export); the reshape's burden is to sever that bridge and
re-target the `.wl` to apply its OWN floor to its OWN computed overlaps.

## The derivation (both engines, own routes)

- **The combined basis + the generalized eigenproblem (POSTULATED basis + CITED collective modes).** The collective modes
  `α_a = sinh(β(1−x))/sinh β`, `α_L = sinh(βx)/sinh β` (normalized `x = w/L0 ∈ [0,1]`, `r_AL = 1` representative) are cited
  from stage 013; the `g`-lane basis is `g_n = sin((n−½)πx)`. Over `B = {α_a, α_L, g_1..g_N}` the mass/stiff Gram is built
  by numeric quadrature — `M_ij = ∫₀¹ φ_i φ_j`, `K_ij = ∫₀¹ [φ_i' φ_j' + β² φ_i φ_j]` — rank-deficient rows are dropped,
  and the generalized eigenproblem `K v = ω² M v` is solved (`eigh(stiff, mass)` / `Eigensystem[{K, M}]`) and sorted
  ascending.

- **The modal overlaps + the truncation floor (EARNED — the anti-gaming crux).** Each of the two lowest generalized
  eigenvectors is projected onto `span{α_a, α_L}` via the mass-Gram (`o_k = ‖P_α v_k‖_M / ‖v_k‖_M`, `k=1,2`). The
  truncation certificate is `pass = o_1 ≥ FLOOR ∧ o_2 ≥ FLOOR ∧ min(ω_1², ω_2²) > 0`, `FLOOR = 0.9` predeclared — a genuine
  predicate on the COMPUTED overlaps. At the physical `β_L0 = 37/20, N=16`: `o_1 = 0.99311, o_2 = 0.98776, min(ω²) =
  3.42252, gap = 2.22787`, `pass = True`. All three conjuncts are separately toothed (the `o_2` floor, the `min(ω²)>0`
  positivity leg, and — via the sweep — the window edge).

- **The β_L0 sweep + the COMPUTED validity window (EARNED — the anti-gaming crux).** Sweeping `BETA_L0_SWEEP = {0.1, 0.2,
  0.5, 1.0, 1.85, 2.0, 3.0, 5.0, 10.0, 18.0, 30.0, 50.0}` (which spans well BEYOND the pass region, to `o_1 ≈ 0.27` at
  `β_L0=50`), the certificate passes for `β_L0 ∈ {0.1 … 3.0}` and **FAILS for `β_L0 ∈ {5, 10, 18, 30, 50}`** (`o_1 = 0.860
  < FLOOR` already at `β_L0=5`). The window `beta_window = {min: 0.1, max: 3.0}` is COMPUTED from the passing set — it has
  a **real upper edge**, not an everywhere-pass. The honest caveat: clean 2-mode truncation holds only for order-unity
  wall stiffness `K_η/T_w ≲ 2.6` (the emergent edge between `β_L0=3` pass and `β_L0=5` fail; `β² = K_η/T_w`, and
  `(3/1.85)² ≈ 2.63`) — an EMERGENT number from the computed sweep, not a typed criterion. The physical `β_L0 = 37/20` is
  the CITED Gate-1 anchor, not the β that maximizes overlap.

- **N-convergence (EARNED).** At the physical `β_L0 = 37/20`, running `N ∈ {4, 8, 12, 16}` gives a stable `o_1`
  (0.99301 → 0.99310 → 0.99311 → 0.99311, max |Δ| ≈ 1e-4), `pass = True` at every N — so the truncation certificate is not
  an artifact of a lucky small N. Each row is asserted at its declared N (a label check), and the stability check is
  genuinely able-to-fail (a deliberately drifting series fires it). The `mass_condition` growth with N (3.1e3 → 1.9e5) is a
  benign conditioning artifact of the enlarging basis, noted as such.

- **The two counterfactual overlaps (EARNED able-to-fail + the honest caveat).** Only the OVERLAP slices of the source's
  counterfactual guard are 014's (`M_det`/`M_posdef` are 013's; `hf_mismatch`/`F_a` are 015's). The `degenerate_zero`
  profile (`α_a ≡ 0`) collapses the collective span → `rank_deficient_basis = True`, `o_2 = 0.2227 < FLOOR` → **FAILS** the
  floor (a genuine 014 tooth). The `constant_one` profile (`α_a ≡ 1`, a WRONG profile) → `o_1 = 1.0, o_2 = 0.9738 ≥ FLOOR`
  → **PASSES** the overlap (`overlap_passes = True`). ⚠ This is the honest scope limit: the overlap certifies
  truncation-consistency, NOT profile-correctness — 014 records that the constant profile PASSES and does **not** claim it
  is rejected (its rejection is 015's Hellmann–Feynman mismatch + 013's residual).

- **The 013→014 M/K seam (no naïve equality).** 014's generalized-eig 2×2 collective block IS 013's `M_AB`/`K_AB` (same
  operator, same profiles), but the numeric Galerkin integrates over the normalized `x = w/L0 ∈ [0,1]` with unit weight,
  so its collective sub-block equals 013's closed forms only up to the `4π`, `L0`-Jacobian, and `μ_η`-weight prefactors —
  which cancel in the generalized eig and the mass-normalized overlaps (both normalization-invariant). The seam is
  therefore guarded via the consumed profiles + packet (below), **not** a naïve `numeric M_aa == symbolic M_aa` (which
  would false-fail). The numeric `r_AL = 1` is the representative point (the overlaps are `r_AL`-invariant).

- **The 014-scoped landing.** Computed from the 014 rungs — `BREATHING_FAIL_TRUNCATION_INCONSISTENT` is the truncation
  verdict rung (`not truncation_status or min(ω²) ≤ 0`); the degenerate-fails, window-has-edge, and N-converged checks are
  able-to-fail ASSERTIONS — landing at the 014 component of `BREATHING_CALIBRATED`. The joint 3-stage composition is
  printed (013 done ∧ 014 here ∧ 015 sibling), with the three carried caveats (overlap ≠ profile-correctness; `K_η/T_w ≲
  2.6`; BdG `k⁴` deferred); the verdict is NOT typed as 014-alone-earned.

## Consumed inputs

**Cited — no file reads; genuine DUAL-SITE citation-integrity with an anti-tautology guard:** stage 013's collective
profiles `α_a, α_L` + the frozen packet `{L0 = 37/20, β = 1, β·L0 = 37/20}` are guarded by two independently-corruptible
relations — **site A (branch anchor)** `β·L0 − 37/20 ≡ 0`, reading an independent CITED `betaL0_cited` datum (NOT the
`β, L0` the sweep uses, so a lone anchor corruption genuinely breaks it — verified non-vacuous); and **site B
(profile-consumption integrity)** the consumed profiles satisfy BOTH the harmonic residual `−α'' + β² α = 0` AND the
collective BC values `{α_a(0)=1, α_a(1)=0, α_L(0)=0, α_L(1)=1}`. ⚠ **The residual alone is insufficient** — a
kernel-preserving corruption (a rescale / wrong endpoint normalization, e.g. `sinh(β(1−x))/cosh β`) stays in `ker(𝓛₀)`
(residual `= 0`) yet has the wrong BC; only the BC leg catches it, while a non-kernel corruption (`α+1`) is caught by the
residual leg — so site B = residual ∧ BC, both able-to-fail. This VERIFIES the citation (the profiles received from 013
are the harmonic lifts); it does NOT re-derive them. Every one-value corruption (the `37/20` anchor, a kernel-preserving
profile, or a non-kernel profile) fails BOTH engines. `c_S` is NOT consumed (matter sector deferred).

## Exports

- The **validity window + the truncation certificate** → **stage 015** (its Hellmann–Feynman force and structure recovery
  run on the same 2-mode closure this stage certifies) + **stages 022/023** (the ℓ=0 cross-ℓ map). Distinct from the ℓ=2
  grouped-`P2` port kernel of stages 016/017.
- Register: **zero new counted knobs** — the numeric truncation controls `{FLOOR = 0.9 (EPS_TRUNC = 0.1), N_FINAL = 16,
  N_CONVERGENCE, BETA_L0_SWEEP}` are method/tolerance parameters (tracked, not counted); the wall packet `{μ_η, T_w, β,
  K_η}` is CONSUMED from stage 013 (no double-count). One structural edge: **R31** (the truncation validity window —
  `K_η/T_w ≲ 2.6`; a validity certificate, dual-engine + able-to-fail, discharging NOTHING).

## Verification

- **Reshape (blueprint §5) — bridge-severing + hybrid-`.wl` numeric re-independence:** stripped the `.py`'s numeric-WL
  scratch export (`SYM_EXPR_WL` append), the scratch-YAML/results writers, the MMA-scratch re-read, and the report/feed/
  summary writers; and the `.wl`'s `Get[sympyExprFile]`, the sweep-grid read `sympyBetaSweep[[All,1]]`, the four `*Delta`
  float-diffs, the `numericChecks` assoc, and the `Export`. The `.wl` KEEPS its native `NIntegrate`+`Eigensystem`
  `galerkinRow` route but its asserts were **re-targeted from `sympy*`-float-differencing to its OWN computed overlaps +
  floor + window**; it declares its OWN sweep grid / floor / N-grid and the cited profile closed forms natively (NOT via
  `DSolveValue`). Both engines standalone, print-only, **zero file I/O**, float-BEARING (numpy/scipy · NIntegrate/
  Eigensystem; predeclared anchors exact — `β_L0 = 37/20`; assertions numeric-tolerance, `5e-6` selected / `5e-4`
  counterfactual / `1e-3` N-stability, all with wide pass margins), ledger idioms. **Clean 3-way cut:** no 013 content
  (harmonic profiles / `M_AB`/`K_AB` `∫dw` projection / EOM LHS — cited) and no 015 content (legacy `E_geom` Hessian /
  `build_structure_gate` / the Hellmann–Feynman force / static-dynamic limit / `hf_mismatch`); `c_S` is not a live symbol.
- **Dual-engine agreement (transcript-level):** the two engines independently compute the sweep pass-pattern, the window
  `[0.1, 3.0]`, the selected `o_k`, the N-convergence, and the two counterfactual overlaps — `scipy.eigh`+`quad` vs
  `Eigensystem`+`NIntegrate`, agreeing far tighter (`max |Δ| ≈ 1e-13`) than the `0.9` floor margin; the scripts do NOT
  cross-read (the arbiter compares transcripts). The `.wl` correctly M-normalizes its own (un-normalized) `Eigensystem`
  vectors; it carries an arity self-check (the stage-007 silent-skip lesson).
- **Dual-engine:** SymPy/SciPy **93 PASS / 0 FAIL** · Mathematica **100 PASS / 0 FAIL** (the +7 = the `.wl`'s native arity
  self-check + unevaluated-leakage block), both exit 0, CWD-independent (repo root + foreign CWD); runner transcripts under
  `scripts/output/` + `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (both engines, both CWDs, reproduced); **`FIDELITY_CLEAN`**
  (an independent from-scratch numpy Galerkin re-derived every value — `o_1 = 0.99310910`, `o_2 = 0.98776370`, `min(ω²) =
  3.42251945`, window `[0.1, 3.0]`, edge-fail `o_1 = 0.85984718` at `β_L0=5`, degenerate `o_2 = 0.22268966`, constant
  passes `o_1=1.0/o_2=0.97384719`, `K_η/T_w` edge `2.630` — all matched; coverage symmetric across engines; 3-way cut,
  anti-gaming threshold, dual-site consumption, and M/K seam faithful); **`ADVERSARIAL_ISSUES`** — the anti-gaming
  floor/window (the v1 scar) proven CLEAN (23/23 genuine ablations fired at their own assert; the `min(ω²)>0` leg LIVE; the
  identity-sub-Gram projection tooth material, `[0.557, 0.382]` vs `[0.993, 0.988]`; dual-site residual∧BC both able-to-fail),
  but **2 BLOCKING falsification-strength defects: teeth 4 and 5 were vacuous `x := <const>; expect_fail(x)` stamps in both
  engines** (dead able-to-fail legs, decoupled from the physics — though the underlying facts stayed genuinely guarded by
  the neighboring baseline `expect_bool` reads).
- **Remediation (teeth 4/5 made genuine, both engines):** rewired each tooth's `expect_fail` to a **shared predicate**
  (`constant_overlap_caveat` / `degenerate_guard_catches`) that the baseline `expect_bool` also uses, applied to a
  **mutated copy** of the real computed row — so the leg fires because the status is corrupted, not because a literal was
  flipped. A **fresh-agent re-verify** confirmed `REVERIFY_CLEAN` via the deciding coupling meta-test: corrupting the
  shared predicate to always-`True` made each tooth's `expect_fail` stop firing → the audit correctly FAILED, proving
  genuine coupling. No regression (93/0, 100/0; all other teeth intact; the vacuous strings gone). A follow-up deleted
  three now-dead `Module` locals from the `.wl` (re-run 100/0).
- **Teeth (9, all fire):** (1a) drop the `o_2 ≥ FLOOR` conjunct → the degenerate no longer rejected; (1b) force `min(ω²) ≤
  0` → the positivity conjunct fires → `BREATHING_FAIL_TRUNCATION_INCONSISTENT`; (2) hardcode a high-β row to pass / trim
  the sweep → the computed window loses its real edge; (3) identity-sub-Gram (a material projection change, NOT the inert
  "skip mass-normalization") → the physical-point overlap assertion fails; (4) a constant that fell below the floor breaks
  the overlap-passes caveat (via the shared predicate on a mutated copy); (5) a non-collapsing degenerate is not caught by
  the guard (shared predicate on a mutated copy); (6) a deliberately non-converging N-series fires the stability assert +
  the N-label check catches an all-N=4 mutant; (7) any one-value corruption of the consumed 013 packet (the `37/20`
  anchor, a kernel-preserving profile → BC leg, a non-kernel profile → residual leg) fails both engines; (8) the `.wl`
  arity self-check catches a def/call mismatch.

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_31_scalar_breathing_{sympy.py,.wl}` (014 slice = `galerkin_overlap`
  L476–557, `build_truncation` L560–607, `_baseline_functions` L447–473, the module constants L50–55, the calls L760–762,
  the `wrong_o_k`/`overlap_passes` guard slices; sources unchanged); `software/stage1_solver/reports/pathA_31_scalar_breathing.md`
  (`## Truncation Consistency` :37–58; `## Counterfactual Guard` :79–82 shared). SIBLING (cited, not recomputed):
  `## Operator, BCs, Inner Product` / `## Profiles And Projection` / `## Dynamical EOM` / `## Dimensional Check` :3–36,
  :87–98 = stage 013; `## Structure` / `## Hellmann-Feynman Force` / `## Static-Dynamic Limit` :60–77 = stage 015.
- Reshape directive + review artifacts — ⛔ **not retained** (they lived in gitignored `_scratch/`; no copy survives, so the names that follow record what existed rather than citing it): `research/pde_ledger_v2/_scratch/ledger_stage014_*` (directive; Codex→Grok→Codex
  design-review logs; execute/remediation/re-verify logs). The directive design-review used the **Codex → Grok → Codex**
  bookend (blueprint §6): Codex `DIRECTIVE_CLEAN` after 3 BLOCKING + 2 nits folded (site-B residual∧BC, the `min(ω²)>0`
  tooth, the N-convergence non-converging ablation, the identity-sub-Gram projection tooth, the native cited-closed-form
  `.wl`), a Grok-4.5 compute-verification pass (`DIRECTIVE_CLEAN`, all anchors independently re-derived), then a Codex
  confirm-pass on the folds.
- Running-start source map: `research/pde_ledger_v2/notes/stage014_pathA31_truncation_source_map.md`. Split row:
  `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 014). Carries the truncation-consistency component of
  the joint `BREATHING_CALIBRATED`; the M/K-projection + (a,L)-closure (013) and the legacy-structure + HF-force (015)
  components complete the fold.
