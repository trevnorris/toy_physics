# ledger_stage013_breathing_harmonic_mk_projection

## Status

**Part II — Gravity. II-G2a (build-order 013).** Reshape of the **harmonic-profile + M/K-projection leg** of gate
**`pathA_31`** (the scalar "breathing" of the frozen finite throat). Source top-line, verbatim: **`BREATHING_CALIBRATED`**
— a JOINT 3-stage verdict, and **this stage carries its M/K-projection + (a,L)-closure component.** ⭐ **pathA_31 splits
3-way** per the finalized Part-II split: **013 (this stage)** = the harmonic-lift profiles + the `M_AB`/`K_AB` operator
projection + the dynamical-EOM LHS; **014** = truncation consistency (generalized eig / β_L0 sweep / N-convergence);
**015** = legacy-Hessian structure recovery + the Hellmann–Feynman force. This is the **FIRST Part-II CALIBRATED stage** —
the gravity-sector algebra was knob-free through stages 008–012; the calibration enters here with the breathing-mode wall
response.

The ℓ=0 (axisymmetric) linearized breathing of the frozen throat reduces to a **2-mode collective-coordinate closure**
`Q = (δ_a, δ_L)`. Given Gate-1's frozen const-coefficient wall operator on the throat interval `[0, L0]`, this stage
solves the two collective harmonic-lift profiles `α_a, α_L` (from `𝓛₀ α = 0` under the collective boundary conditions),
projects the wall action's inertia and stiffness onto them by **real `∫dw` operator projection** to build the
mass/stiffness matrices `M_AB, K_AB`, and assembles the dynamical-EOM **LHS** `M_AB Q̈ + K_AB Q`.

- **CONSUMED (cited, dual-site integrity, NOT re-derived):** the **frozen throat packet** from Gate-1
  (`ledger_stage011`+`ledger_stage012`): the domain `[0, L0]` (cap `R0(L0)=0`), the frozen const-coefficient wall packet
  `L0 = 37/20, T_w = 1, K_η = 1, β = 1` (`K_η = T_w β²`, `β·L0 = 37/20` = the branch-determinable `L/a`), and the ℓ=0
  restriction (`Y00 = 1/(2√π)`; the `T_Ω` angular term drops at `ℓ(ℓ+1)=0`). ⚠ **`c_S` is NOT consumed** — 013's symbolic
  content is speed-free; the matter-sector `c_s`/BdG `k⁴` is explicitly DEFERRED under `kξ≪1` (`phonon_limit_caveat`).
- **POSTULATED (labeled):** the collective boundary conditions `{α_a(0)=1, α_a(L0)=0}`, `{α_L(0)=0, α_L(L0)=r_AL}` (the
  IMPOSED collective-coordinate normalizations of the straight-throat lift). `r_AL` (the `α_L` cap ratio, `[r_AL]=1`) is a
  dimensionless control-construction ratio that parameterizes the collective length mode, not the physics.
- **EARNED (structure):** the harmonic-lift profiles `α_a, α_L` (DERIVED, proven to annihilate `𝓛₀`); the `M_AB`/`K_AB`
  matrices by real `∫dw` operator projection (NOT typed from the legacy static Hessian — `forbidden_fit_flags` computed
  `False` via free-symbol-**name** ancestry, `K_AB_provenance = operator_projection_not_static_Hessian`); the
  dynamical-EOM LHS.
- **CALIBRATED (values):** the wall constants `μ_η` (inertia), `T_w` (tension), and `β` (inverse-length scale) are
  calibration inputs, with `K_η = T_w β²` a derived manifestation. **The STRUCTURE is EARNED; the VALUES are CALIBRATED**
  → the top-line is `BREATHING_CALIBRATED`, not `..._PASS`. Geometry alone does not derive the wall constants (the source
  is explicit: `beta_from_R0: "geometry alone does not derive it"`).

Ledger-local earned-label (NOT a source verdict token): `BREATHING_HARMONIC_MK_PROJECTION_EARNED`. The joint verdict
composes as `BREATHING_CALIBRATED = (013: harmonic-lift profiles + M_AB/K_AB by ∫dw operator projection + (a,L) EOM LHS,
computed here) ∧ (014: truncation consistency) ∧ (015: legacy-structure + Hellmann–Feynman force)`.

## Purpose

Stages 011/012 pinned the frozen-throat interior operator and its D/N resonator response. This stage asks the next
question: how does the throat *breathe* — how does the ℓ=0 collective deformation of the frozen wall respond
dynamically? The physical content is a reduction of the distributed wall field `η(w,t)` to two collective coordinates
`(δ_a, δ_L)` (a "mouth-radius" mode and a "length" mode), and the earned deliverable is the mass/stiffness matrices of
that reduced dynamics — **built by projecting the wall operator onto the derived collective profiles, not read off a
static energy Hessian.** The reduced `M_AB Q̈ + K_AB Q` is the left-hand side of the collective equation of motion; its
right-hand side (the Hellmann–Feynman driving force) and the consistency of the 2-mode truncation are the sibling stages
015 and 014. The result — a genuine operator-projected `(a,L)` closure — is banked forward to stage 014 (which consumes
`M_AB`/`K_AB` for the generalized eigenproblem) and stages 022/023 (the ℓ=0 cross-ℓ map). The source pair computed this
through a scratch-YAML bridge with a hybrid `.wl` (a genuine `DSolveValue`+`Integrate` route that nonetheless imported
and cross-checked the `.py`'s expressions); the reshape's burden is to sever that bridge, re-target the `.wl` to assert
against its OWN native M/K, and keep the operator-projection genuineness un-stampable.

## The derivation (both engines, own routes)

- **The operator, BCs, inner product, ℓ=0 restriction (CITED / POSTULATED).** The frozen wall operator is
  `𝓛₀ = μ_η⁻¹(−∂_w(T_w ∂_w ·) + K_η ·)` on `w ∈ [0, L0]`, with `K_η = T_w β²` and `β = √(K_η/T_w)`. The reduction is
  axisymmetric (ℓ=0): `η(w,t) = η₀₀(w,t)·Y00` with `Y00 = 1/(2√π)` and `∫_{S²} Y00² dΩ = 1`, and the angular `T_Ω` term
  drops because `ℓ(ℓ+1)=0`. The projection carries the μ-weighted inner product `⟨f,g⟩_μ = 4π ∫₀^{L0} μ_η f g dw`.
  (To avoid the source's `L0` name-overload, the interior operator's application is written `Lop[α]` and `L0` denotes the
  domain length.)

- **The harmonic-lift profiles (EARNED).** The general solution of `𝓛₀ α = 0` (⟺ `α'' − β² α = 0`) is
  `α = C₁ sinh(βw) + C₂ cosh(βw)`. Solving for `(C₁,C₂)` under each collective BC pair gives
  `α_a = sinh(L0 β − β w)/sinh(L0 β)` (`α_a(0)=1, α_a(L0)=0`) and `α_L = r_AL sinh(β w)/sinh(L0 β)`
  (`α_L(0)=0, α_L(L0)=r_AL`); `profile_provenance = derived`. The **harmonic-residual assertion**
  `Lop[α] = (−T_w α'' + K_η α)/μ_η ≡ 0` (a RAISE on failure) proves the profiles genuinely annihilate the operator — and
  ⚠ its able-to-fail tooth uses a **non-kernel** corruption (a wrong wavenumber `sinh(L0β−2βw)` or `α+1`, NOT a rescale or
  argument sign-flip, which stay in `ker(𝓛₀)` and leave the residual `0`; the tri-review's decisive AB2b ablation
  confirmed a kernel-preserving corruption is correctly caught only by the *separate* BC assert, proving the residual
  check is genuinely about the residual). The BC values are asserted separately.

- **`M_AB` / `K_AB` by real `∫dw` operator projection (EARNED — the crux).** With M integrands `μ_η α_A α_B` and K
  integrands `T_w α_A' α_B' + K_η α_A α_B`, the projections are
  `M_AB = 4π ∫₀^{L0} μ_η α_A α_B dw` and `K_AB = 4π ∫₀^{L0} [T_w α_A' α_B' + K_η α_A α_B] dw`, computed by `sp.integrate`
  / native `Integrate`. Closed forms (both engines, agreeing):
  - `K_aa = 4π T_w β / tanh(L0 β)`, `K_aL = −4π T_w β r_AL / sinh(L0 β)`, `K_LL = 4π T_w β r_AL² / tanh(L0 β)`,
    `det(K) = 16π² T_w² β² r_AL²`;
  - `M_aa = −2π μ_η (L0 β tanh(L0 β) − sinh²(L0 β)) / (β sinh²(L0 β) tanh(L0 β))` (and the `aL`, `LL`, `det(M)` forms),
  - `[M_AB] = M`, `[K_AB] = M T⁻²`, `[K/M] = T⁻²`.
  **The genuineness spine** (the anti-"typed-M/K" burden, since M/K would otherwise be hardcodable — one could type the
  legacy static Hessian `{χ²κ+σ_A, −χκ, κ+σ_L}`):
  1. **Independent re-integration (NOT a self-identity).** The re-integrated full integrands are compared against the
     **hardcoded report closed-forms** (not against the same-function projection output — that would be a vacuous `X≡X`).
     A corrupted integrand (dropping the `T_w α'α'` gradient term) changes the re-integral and the comparison FAILS
     (tooth-3, confirmed firing in both engines). *(This retarget-to-closed-forms was a tri-review remediation of the
     one non-blocking nit: the baseline comparison had originally compared two identical projection calls.)*
  2. **Free-symbol-NAME ancestry.** The `forbidden_fit_flags` (`K_from_static_hessian`, `M_or_K_typed_from_legacy_values`,
     `full_matrix_fit`) are COMPUTED from the symbol **names** of `M_AB`/`K_AB` (`{s.name …}`, robust to SymPy's
     assumption-fragility where `Symbol('kappa') ≠ Symbol('kappa', positive=True)` and an object-set intersection would
     silently miss): the names must lie in `{muEta, Tw, beta, L0, rAL, pi}` and contain NONE of the legacy structure names
     `{kappa, chi, sigmaA, sigmaL}`. Typing `K_AB` from the legacy Hessian introduces those names → the flag flips `True`
     → CAUGHT in both engines (tri-review AB1/WL-AB1).

- **The dynamical-EOM LHS (EARNED).** The collective EOM LHS is `M_AB Q̈^B + K_AB Q^B = F_A^(HF)` with `Q = (δ_a, δ_L)`,
  emitted as the two expanded rows `M_aa δ̈_a + M_aL δ̈_L + K_aa δ_a + K_aL δ_L = F_a` (+ the `L`-row twin). ⚠ This stage
  emits the LHS + the EOM structure ONLY; the RHS `F_A^(HF)` (the Hellmann–Feynman force) is a symbolic placeholder filled
  by stage 015 — it is not computed here (and the HF `x−x` two-route enforcement, the pathA_31 v1 rig locus, is 015's).

- **The 013 dimensional legs + corrupt-`[T_w]` probe (EARNED).** From the sourced dims (order `{L,M,T}`) `L0:(1,0,0)`,
  `β:(-1,0,0)`, `μ_η:(-1,1,0)`, `T_w:(1,1,-2)`, `r_AL:(0,0,0)`, and derived `K_η = T_w β²:(-1,1,-2)`, the walk gives
  `[M_AB]=(0,1,0)=M`, `[K_AB]=(0,1,-2)=M T⁻²`, `[K/M]=(0,0,-2)`. The corrupt-`[T_w]` probe adds one power of `L` to
  `[T_w]`, shifting `[K_η]` and hence `[K_AB]` to `(1,1,-2)=L M T⁻²` (off `M T⁻²`, while `[M_AB]` stays `M`) →
  `BREATHING_FAIL_DIMENSIONAL` (`mutation_fires=True`), with a two-verdict self-ablation (with mutation → the fail;
  without → `BREATHING_CALIBRATED`, `fail_suppressed=True`).

- **The degenerate `M_det→0` counterfactual (EARNED, native).** Setting the degenerate profile `α_a → 0` and recomputing
  `M` via the SAME `∫dw` projection gives `M_det ≡ 0` (the closure degenerates) — asserted CAUGHT via a **native**
  non-degeneracy test (`M_det ≡ 0`), NOT the 015 structure-gate `M_posdef` (importing that would breach the 3-way cut).

- **The 013-scoped landing.** Computed from the 013 rungs — `BREATHING_FAIL_DIMENSIONAL` is the dimensional verdict rung;
  the harmonic-residual, the M/K-projection-genuineness (free-symbol-name ancestry + the independent re-integration), and
  the degenerate-`M_det` checks are able-to-fail ASSERTIONS that RAISE — landing at the 013 component of
  `BREATHING_CALIBRATED`. The joint 3-stage composition is printed (013 ∧ 014 ∧ 015), with `CALIBRATED ⇐ {μ_η, T_w, K_η}
  calibration`; the verdict is NOT typed as 013-alone-earned.

## Consumed inputs

**Cited — no file reads; genuine DUAL-SITE citation-integrity with an explicit anti-tautology guard:** the frozen wall
packet `{L0 = 37/20, T_w = 1, K_η = 1, β = 1}` (from Gate-1, stages 011/012) is guarded by two independently-corruptible
relations — **site A (constitutive)** `β − √(K_η/T_w) ≡ 0` and **site B (geometric/branch)** `β·L0 − 37/20 ≡ 0` (the
branch-determinable `L/a`) — plus the frozen-export anchor `(L0, T_w, K_η, β) = (37/20, 1, 1, 1)`. ⚠ **Anti-tautology
(the crux of the guard):** the packet values are held as independently-corruptible CITED datums (`K_eta_cited`,
`Tw_cited`, `beta_cited`, `L0_cited`) — the guard does **not** re-expand `K_η` through the operator alias `K_η = T_w β²`,
because `β − √(T_w β²/T_w) = β − β ≡ 0` would be a vacuous `X≡X` (the disguised-tautology rig class that bit stages
006/007/009). A lone `K_eta_cited` corruption genuinely breaks site A; corrupting the `37/20`/`L0` breaks site B; every
one-value corruption fails BOTH engines (tri-review AB3/WL confirmed, and AB3b confirmed that substituting the alias makes
the non-vacuity tooth itself fail). This checks the packet's internal consistency; it does **not** derive the wall
constants (they are calibration). The domain `[0, L0]` is cited from stage 011 (`L0` already registered as
ACTION-geometry — not re-counted). `c_S` is NOT consumed (matter sector deferred).

## Exports

- The `(a,L)` collective closure — the profiles `α_a, α_L`, the `M_AB`/`K_AB` matrices, and the dynamical-EOM LHS — →
  **stage 014** (which consumes `M_AB`/`K_AB` for the combined-basis generalized eigenproblem) + **stages 022/023** (the
  ℓ=0 cross-ℓ map). Distinct from the ℓ=2 grouped-`P2` port kernel of stages 016/017 — do not conflate the ℓ=0 breathing
  closure with the ℓ=2 sector.
- Register: **⭐ the FIRST Part-II counted calibration knobs — 3 counted CALIB `{μ_η, T_w, β}`** (the frozen-throat-wall
  constitutive packet: inertia `[μ_η]=M L⁻¹`, tension `[T_w]=M L T⁻²`, inverse-length scale `[β]=L⁻¹`), with **`K_η = T_w
  β²` a DERIVED manifestation (edge R29)** — so 3 counted, not 4. ⚠ `β` is counted, **not** dressed as geometry (the
  source is explicit that geometry alone does not derive the wall constants; `β·L0=37/20` is `β(calibrated)·L0(=the
  stage-011 L/a geometry)`, not an independent branch-pin of `β` — the conservative, non-shrinking read). `r_AL`
  (dimensionless collective BC ratio) is a control-construction `CANDIDATE` (tracked, not counted, like `k_warp`/`α`).
  One **reduction-debt** edge: **R30** (a nonlinear-throat solve would derive the wall response `{μ_η, T_w, β}` from the
  medium — `PENDING`, the same deferred throat as R10/R23).

## Verification

- **Reshape (blueprint §5) — bridge-severing + hybrid-`.wl` re-independence:** stripped the `.py`'s scratch-YAML /
  `_sympy_exprs.wl` export, the numeric-Galerkin append, the MMA-YAML re-read (`engine_status`), the report/feed/summary
  writers, and the `RESULTS_YAML`/`SYM_YAML` writes inside `main`; and the `.wl`'s `Get[sympyExprFile]`, the
  `sympy*`-comparison `checks` assoc, the numeric-overlap-vs-`sympy*`, the digest write-back, and the `Export`. The `.wl`
  KEEPS its native `DSolveValue`+`Integrate` route but its asserts were **re-targeted from `sympy*`-differencing to its
  OWN native M/K/profiles**. Both engines standalone, print-only, **zero file I/O**, float-free (`L0 = 37/20` as an exact
  rational), ledger idioms. **Clean 3-way cut:** no 014 content (generalized eig / `galerkin_overlap` / `build_truncation`
  / β_L0 sweep / N-convergence / `beta_from_R0`) and no 015 content (legacy `E_geom` Hessian / `build_structure_gate` /
  the Hellmann–Feynman force / static-dynamic limit); the EOM RHS `F_A^(HF)` is a symbolic placeholder; `c_S` is not a
  live symbol.
- **The `.wl` independent route:** the KEPT `DSolveValue`+`Integrate` engine solves `−y'' + β² y = 0` under the collective
  BCs and integrates the operator integrands natively (Mathematica-native `Csch`/`Coth` closed forms, e.g.
  `K_aa = 4β π T_w Coth[β L0]`, algebraically equal to the SymPy `tanh/sinh` forms) — a genuinely different route from the
  `.py`'s `sp.solve`-on-general-solution + `sp.integrate`; both agree. Its own dual-site packet guard, free-symbol-name
  ancestry, native degenerate-`M_det` slice, dim legs + corrupt-`[T_w]`, and an arity self-check (which genuinely catches
  a def/call mismatch — the stage-007 silent-skip lesson).
- **Dual-engine:** SymPy **78 PASS / 0 FAIL** · Mathematica **84 PASS / 0 FAIL** (the +6 = the `.wl`'s native arity
  self-check block), both exit 0, CWD-independent (repo root + foreign CWD); runner transcripts under `scripts/output/` +
  `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (both engines, repo root + foreign CWD, reproduced);
  **`FIDELITY_CLEAN`** (every profile/M/K/det independently re-derived to ~1e-124 via a *different* normalization path;
  no dropped check; `expand_hyperbolic_args` confirmed a math-preserving identity — `sinh(β(L0−w)) ≡ sinh(βL0−βw)`, same
  function — that only lets `sp.integrate` evaluate; genuine dual-engine; 3-way cut clean; FIRST-calibration honesty
  intact); **`ADVERSARIAL_CLEAN`** (10-ablation mutation matrix, all fired at the intended assert: the M/K free-symbol-name
  ancestry catches a legacy-Hessian mutant in both engines; the non-kernel residual tooth holds (AB2b decisive); the
  `K_eta_cited` anti-tautology is real in both engines (AB3b: substituting the alias makes the non-vacuity tooth fail);
  the native degenerate-`M_det` guard fires; the corrupt-`[T_w]` probe propagates; the `.wl` is a genuine independent
  engine, not a mirror; arity self-check catches a planted mismatch; tallies reconcile 78 = 72 + 6-arity-net, 84).
- **Remediation (one non-blocking nit, both legs CLEAN):** both legs flagged that the baseline "independent
  re-integration" comparison was value-wise `X−X` (a second `project_from_profiles` call compared to the first), passing
  by construction — genuineness carried by the closed-form assertions + tooth-3. Retargeted (both engines) to compare the
  re-integration against the hardcoded report closed-forms, making each of the 6 asserts genuinely able-to-fail; a
  **fresh-agent re-verify** confirmed `REVERIFY_CLEAN` (closed-forms are hardcoded literals — not a relabeled `X−X`; no
  self-compare remains; tooth-3 still fires in both engines; baselines 78/84; nothing else changed).
- **Teeth (6, all fire):** (1) harmonic-residual — a NON-KERNEL profile corruption (`sinh(L0β−2βw)` or `α+1`) → residual
  ≠ 0 → raise; (2) M/K typed from the legacy Hessian → `{kappa,chi,sigmaA,sigmaL}` in the name set → flag flips → caught;
  (3) drop the `T_w α'α'` gradient integrand → the independent re-integration mismatches the closed-forms → fail; (4)
  corrupt-`[T_w]` → `[K_AB]` off `M T⁻²` → `BREATHING_FAIL_DIMENSIONAL` + two-verdict self-ablation; (5) dual-site
  consumed-packet — any one-value corruption (`K_eta_cited`, `37/20`, `L0`, `T_w`) fails both engines (site A verified
  non-vacuous); (6) degenerate — bypass the native `M_det≡0` test so `α_a=0` is not caught → the guard assertion fails.

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_31_scalar_breathing_{sympy.py,.wl}` (013 slice = `symbolic_engine()`
  L627–679; sources unchanged); `software/stage1_solver/reports/pathA_31_scalar_breathing.md` (`## Operator, BCs, Inner
  Product` :3–10, `## Profiles And Projection` :11–30, `## Dynamical EOM` :31–36, `## Dimensional Check` :87–98,
  `## Reduction Certificate` :99–102). SIBLING (cited, not recomputed): `## Truncation Consistency` :37–58 = stage 014;
  `## Structure` / `## Hellmann-Feynman Force` / `## Static-Dynamic Limit` :60–77 = stage 015.
- Reshape directive + review artifacts — ⛔ **not retained** (they lived in gitignored `_scratch/`; no copy survives, so the names that follow record what existed rather than citing it): `research/pde_ledger_v2/_scratch/ledger_stage013_*` (directive; Codex→Grok→Codex
  design-review logs; execute/remediation logs). The directive design-review used the new **Codex → Grok → Codex**
  bookend (blueprint §6): Codex `DIRECTIVE_CLEAN`, a Grok-4.5 compute-verification pass (which caught the kernel-preserving
  residual-tooth defect folded into §4.1 + five hardening nits), then a Codex confirm-pass on the folds.
- Running-start source map: `research/pde_ledger_v2/notes/stage013_pathA31_breathing_source_map.md`. Split row:
  `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 013). Carries the M/K-projection + (a,L)-closure
  component of the joint `BREATHING_CALIBRATED`; the truncation-consistency (014) and legacy-structure + HF-force (015)
  components complete the fold.
