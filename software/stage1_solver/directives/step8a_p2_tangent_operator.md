# Build Directive — §7 Step 8a: ℓ=2 P2 weak-axisymmetric tangent operator + static (ω=0) tangent solve

**Status:** build directive (Claude drafts requirements + acceptance criteria; Codex designs the discretization and writes the code; Claude reviews with clean agents + a transliteration-fidelity audit + an independent arbiter re-run).

**Plan anchor:** `docs/branch_realization_execution_plan.md` §7 step 8 ("WP3 — P2 tangent on the converged stationary branch"). Step 8 is decomposed (user-approved, mirrors the step-3 3a/3b split) into:

- **8a (THIS directive)** — the ℓ=2 P2 weak-axisymmetric *tangent operator* + a *static* (ω=0) driven tangent solve + MMS certification of the new operator pieces. No waves, no absorber, no physical response packet.
- **8b (later)** — the genuine outgoing-wave absorber (CAP / matter-wave PML) upgrading the Robin exit + a measured wave-reflection coefficient + the *driven* (ω≠0) response. (Deferred from step 5.)
- **8c (later)** — the genuinely current-carrying conservation test (deferred from step 6) + the low-frequency-expansion surrogate response coefficients.

This step stays **engineering-smoke / target-blind**, exactly like steps 3–7: placeholder parameters, NO §G/§H observable extraction, NO frozen `free_choice` values (GATE A not done), export guard untouched (GATE B not done).

---

## 1. What the ℓ=2 P2 tangent is (transliterate from the ledger — do NOT invent)

The WP3 sector linearizes the converged WP1 stationary isotropic branch in the **ℓ=2 (grouped real P2) weak-axisymmetric angular sector**. Canonical source — `notes/moving_throat_pde_program_compact.md`:

- **Linearized coupled bulk/interface skeleton — §4.7, L1383–1455.** Around a stationary reference `ψ = e^{−iμ₀t/ℏ}ψ₀(X)`, `A_M = A_{M0}(X)`, `R = R₀(w)`:
  - **Matter sector (BdG-like):** `iℏ ∂_t [δψ, δψ*]ᵀ = L_BdG [δψ, δψ*]ᵀ + C_A[δA_M] + C_η[η]` (L1406–1423).
  - **Maxwell sector (boxed, L1431–1439):** `∂_M(Z(w) δF^{MN}) + (1/ξ) ∂^N(∂·δA) = μ₀ δJ^N`.
  - **Geometry/wall sector (boxed, L1441–1451):** `μ_η ∂_t²η − ∂_w(T_w ∂_w η) − T_Ω Δ_{S²}η + K_η η = S_η^(ψ) + S_η^(A) + f_ext`.
- **Angular split — §4.3, L1298–1308.** The modal wall operator is `μ_η q_{lm,tt} − ∂_w(T_w ∂_w q_{lm}) + [K_η + l(l+1) T_Ω] q_{lm} = S_{lm}^(ψ,A) + f_{lm}^ext`. With `Δ_{S²} Y_{lm} = −l(l+1) Y_{lm}`, the scalar lane `l=0` and the grouped real P2 lane `l=2` are split *before* any matter/gauge closure. **For l=2, l(l+1)=6.**
- **Bulk angular reduction.** The 3D radial Laplacian is `Δ_3D = ∂_r² + (2/r)∂_r + (1/r²)Δ_{S²}`; on a `Y_{2m}` mode `Δ_{S²} → −6`, so the matter (and scalar-potential) bulk Laplacian acquires the **centrifugal term `−l(l+1)/r² = −6/r²` × field**. This is the principal NEW operator piece versus the isotropic (`tensor_laplacian`, no centrifugal) operators of steps 1–7.
- **Grouped real P2 / isotropic degeneracy — §4.5, L1340–1364.** On the isotropic reference throat the five grouped real P2 channels `{20, 21c, 21s, 22c, 22s}` (L1107–1109) are **degenerate**: the radial tangent operator depends on `ℓ` only, not on `m`. So 8a builds the ℓ=2 *radial* tangent operator (m-independent for an isotropic background); the m-multiplicity is a downstream extraction/bookkeeping detail, not solved five times here.
- **Pinned wall→matter coupling — §3.3, L1075–1089 (boxed, "Exact Within Closure"):** `δV_conf = −(V_wall'(Σ₀/ℓ_c)/ℓ_c) η` — "the basic linear source through which the moving wall drives the matter and gauge sectors." This direction of coupling IS canonically pinned and should be transliterated.

**The tangent operator must BE the linearization of the WP1 residual** in the ℓ=2 sector — i.e. the Jacobian of the same coupled residual used in step 3, restricted to ℓ=2 angular modes (which differs from the ℓ=0 Jacobian *only* by the angular/centrifugal terms above and the now-live `l(l+1)T_Ω` wall term). It must NOT be an independently re-posited operator. This is enforced by acceptance gate **A** (tangent-vs-finite-difference-of-residual).

---

## 2. Pinned vs. OPEN coupling — the scope boundary (critical)

`§4.6 (L1377–1381)` lists **"the full coupled matter/gauge renormalization of these reduced lanes"** and **"the outgoing odd response and quadrupole-normalization branch"** as **STILL OPEN**. Concretely:

- **PINNED (transliterate, use in 8a):** the matter BdG-like linear operator `L_BdG`; the gauge coupling `C_A[δA]` (the same covariant `D_i` coupling already in `coupled_pde_residual`); the boxed Maxwell operator; the boxed wall geometry operator with the `l(l+1)T_Ω` term (l=2); the wall→matter coupling `δV_conf` (L1080–1086); the centrifugal `−6/r²` bulk term.
- **OPEN / schematic (do NOT invent a kernel):** the matter/gauge → wall *source* `S_η^(ψ)`, `S_η^(A)` / `S_{lm}^(ψ,A)` (named at L1304, L1449 but given no explicit kernel; §4.6 marks the coupled renormalization open). This is the physical "induced P₂₂ source." It is the engine of the *physical* response packet, which belongs to 8c and may itself require a Claude+Codex methodology call (or a user-scope decision, like the effective-closure decision) if it is genuinely not derivable from the frozen ledger.

**8a resolution:** drive the tangent with a **target-blind surrogate external forcing `f_ext` / `f_{lm}^ext`** (a generic, non-physical source that excites all coupled sectors so the convergence study is a real test of the coupled ℓ=2 operator), NOT the physical matter/gauge→wall backreaction source. If, on reading the source, Codex finds the matter→wall source `S_η^(ψ,A)` is in fact explicitly pinned somewhere it can cite verbatim, it may use it — but absent a verbatim canonical kernel it must **STOP-flag** rather than posit one (mirrors the step-2 wall-constitutive STOP and the step-3 wall↔branch-balance STOP). Inventing the matter→wall coupling would be a fidelity violation.

The wall constitutive packet (`μ_η, T_w, T_Ω, K_η`) remains `free_choice` with **MMS-only / engineering-smoke PLACEHOLDER** values frozen per-run and labelled in config/code/report (identical convention to step 2). The static balance `R₀(w)` is the prescribed step-3 background; no invented force law.

---

## 3. In scope — deliverables (8a)

1. **ℓ=2 tangent operator** (new module, e.g. `p2_tangent.py`): a linear operator that is the ℓ=2 linearization of the step-3 coupled residual around a converged WP1 background `(ψ₀, A₀, R₀)`, including the centrifugal `−6/r²` matter/scalar-potential term, the ℓ=2 angular structure of the Maxwell sector, the live `[K_η + 6 T_Ω]` wall term, and the pinned `δV_conf` wall→matter coupling. The WP1 background is the **step-3 engineering-smoke** converged branch (`run_branch_continuation`, placeholder params), NOT a frozen physical branch.
2. **Static (ω=0) driven tangent solve**: solve `T_{ℓ=2} δu = f` (elliptic at ω=0) with a target-blind surrogate forcing `f` that excites all coupled sectors, under the existing Robin/sponge open-exit BCs (no absorber yet), reusing the existing Newton/GMRES + colored-sparse-Jacobian preconditioner path (the tangent is linear, so this is a single linear solve or a one-Newton-step solve; either is acceptable — Codex chooses).
3. **MMS certification of the NEW operator pieces** (added to `mms_benchmarks.py` / `run_all_mms_benchmarks`): non-circular manufactured-solution benchmarks for the centrifugal matter term and the ℓ=2 Maxwell angular structure, converging at order ~2. (The wall `l(l+1)T_Ω` term is already MMS-covered at ℓ=2.) Forcing derived from the continuum operator via sympy, NEVER from the discrete stencil under test.
4. **Target-blind surrogate tangent observables**: raw-field functionals of `δu` (e.g. response-norm `‖δu‖`, per-sector response norms, a mouth ℓ=2 deformation surrogate) — NO §G/§H observables (no `D2/D4/N2/N4/P2/P4`, no `R_pole/R_norm`, no `χ_Q`, no extraction map).
5. **Grid-convergence of the surrogate observables** (reuse the step-4 ladder pattern): each surrogate observable converges at ~order 2 under refinement; honest floor/null labels in the step-4 vocabulary.
6. **Harness** (`p2_tangent_harness.py`): thin CLI — obtain the WP1 background → build the tangent → static solve → write report + machine-readable table + deterministic digest.
7. **Report** `reports/step8a_p2_tangent.md` + tests in `tests/test_stage1_solver.py`.

---

## 4. Out of scope (→ 8b / 8c / gated)

- The driven (ω≠0) response and its low-frequency expansion → `R_pole`/`R_norm` (compact response sections ~L2542, L2677, L2728). → **8c**.
- The genuine outgoing-wave absorber (CAP / PML) + measured reflection coefficient (deferred from step 5). → **8b**.
- The genuinely current-carrying conservation test (deferred from step 6). → **8c**.
- The physical matter/gauge→wall source `S_η^(ψ,A)` / the physical induced-P₂₂ response. → **8c** (likely a STOP / methodology call; see §2).
- ANY §G/§H observable extraction, the field→coefficient extraction map, the frozen physical run. → gated behind USER GATE A (freeze free_choice values target-blind — none frozen) + GATE B (export-guard flip) + the extraction map (still TODO).

---

## 5. Honesty & target-blindness constraints (review axis #1)

- **Target-blindness is the #1 review axis.** Grep-clean: no `R_norm`/`R_pole`/`χ_Q`/`N_Q` targets, no GR quadrupole constant `54 G c_s^5/(5 a^5 c^5)` or `54/5`, no `χ_Q=1`, no §H target values, no extraction map, no benchmark/reference import in the solve path. `ℓ=2` is a STRUCTURAL choice fixed by prereg §B (grouped real P2), NOT a tuned parameter, NOT chosen to hit a target.
- **Engineering smoke.** Placeholder/natural O(1) parameters; no frozen free_choice values anywhere; export guard untouched; emit NO physical packet.
- **Surrogate, not physical.** State loudly in the report that the observables are target-blind numerical surrogates demonstrating operator + solve machinery, NOT the physical WP3 observables (`d ln R_tr`, `N_Q`, the coefficient bundle) — those need the extraction map + GATES A/B + (per §2) the open matter→wall source.
- **No can't-fail pass-checks** (the recurring step-6/step-7 sin). Every counted gate must be a genuinely computed boolean that can FAIL on a broken build. Any asserted-by-construction claim (target-blind, export-guard-false, prior-step-passed) goes in a separate `asserted_checks` section with the `_not_a_physics_gate` suffix and is EXCLUDED from `passed`. Do NOT "make the export-guard check real" by importing the firewalled `research/pde_audit/simulation/physical_nonlinear_model.py`.
- **MMS-clean ≠ fidelity** (the load-bearing lesson). A wrong ℓ=2 operator still converges at order 2 against its own wrong sympy forcing. The transliteration-fidelity audit (Claude, post-build) is the load-bearing correctness check, not the convergence table.

---

## 6. Reuse / firewall

- Reuse: `coupled_branch.py` (residual + pack/unpack + background solve), `operators.py`, `newton.py`, `preconditioners.py`, `mms.py`/`mms_benchmarks.py`, `boundaries.py`, `config.py` (config_hash payload + exclusions), `manifest.py`, the report/digest patterns.
- **Firewall:** never import from / write under `research/pde_audit/simulation/`. Run output stays under the gitignored `runs/`/`data/`/`figures/`/`_scratch/`.
- **Sponge×MMS gap:** 8a adds NO differential absorber, so the sponge stays a pointwise term and the gap stays LOW. Flag in the report that 8b's absorber, if it enters a differential stencil, MUST be added to the MMS forcing first.

---

## 7. Acceptance criteria (the gate sequence will check these)

- **A. Tangent IS the Jacobian (load-bearing).** `T_{ℓ=2}·δu` matches the central finite difference of the full ℓ=2 residual, `[R_{ℓ=2}(u₀+εδu) − R_{ℓ=2}(u₀−εδu)]/(2ε)`, to ≲1e-9 on a small grid (mirrors the existing JVP-vs-FD probe). Proves the tangent is the linearization, not an independent posit.
- **B. MMS order ~2, non-circular**, for the centrifugal matter term and the ℓ=2 Maxwell angular structure; pytest + harness.
- **C. Static tangent solve converges** to the solver floor via the existing Newton/GMRES + preconditioner path.
- **D. Surrogate observables converge** at ~order 2 across a ≥3-level refinement ladder; honest floor/null labels.
- **E. ℓ=2 tangent well-posed at ω=0** — demonstrate the static ℓ=2 operator is non-singular (the ℓ=0 phase/normalization zero-mode is absent in the ℓ=2 sector); e.g. a conditioning/smallest-singular-value check on a small grid, or clean Newton/GMRES convergence with no null-space stagnation.
- **F. Target-blind** — grep-clean per §5; a test asserts no target constants in the solve path.
- **G. Deterministic** — diagnostics digest (sha256[:16]) reproduced across processes.
- **H. Honest pass-checks** — counted gates are genuinely computed; asserted items separated per §5.

A directive item Codex cannot satisfy without inventing open physics (per §2) is a **STOP-flag**, not a guess.

---

## 8. Tests (add to `tests/test_stage1_solver.py`)

At minimum: tangent-vs-FD-of-residual consistency (gate A); MMS-order for each new operator piece (gate B); ℓ=2 well-posedness (gate E); a target-blindness assertion (gate F); digest determinism (gate G); a non-tautological test that the centrifugal term actually changes the residual versus ℓ=0 (so the ℓ=2 structure is provably live, mirroring step-3's "cross-sector terms LIVE" test).

---

## 9. Report requirements (`reports/step8a_p2_tangent.md`)

Scope framing (engineering-smoke, target-blind, surrogate-not-physical); the canonical-source citations (compact §4.3/§4.5/§4.7 line ranges) for each transliterated operator piece; the pinned-vs-open coupling boundary (§2) stated explicitly, including the matter→wall source deferral; the MMS table; the tangent-vs-FD check; the static-solve convergence; the surrogate-observable convergence study with honest floor/null labels; counted pass-checks vs asserted-by-construction checks; provenance (background commit/config_hash); reproduction commands; the sponge×MMS forward note for 8b.
