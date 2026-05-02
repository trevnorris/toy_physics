# Critical Referee Audit — Fifth Round

**Audited:** the 18 new step pairs (`step_20` through `step_37`), the small `step_19` delta (branch-freeze metadata + `R_pole/R_norm/R_P2/R_P4` labels), the new top-level `README.md`, and the four standalone runtime tools (`cfd_runtime_monitor_postprocess.py`, `cfd_runtime_failfast.py`, `cfd_snapshot_adapters.py`, `single_throat_monopole_jsonl_fastscreen.py`). Total ≈ 6700 new lines of `.md` + `.py`.

**Method:** five referee agents, each given (a) the round-4 verdict matrix and priority list from `AUDIT_REPORT_4.md`, (b) the anti-pattern catalogue (§2.1–§2.10), (c) bundle-specific risk notes (frontier/admissibility/burden circularity, §2.11 self-attested provenance flags, runtime-tool §2.11 inline re-implementation risk). Bundles: A — step 19 delta + 20, 21; B — 22-25; C — 26-31; D — 32-37 self-tests; E — the four CFD utility scripts + README. Trust assumption: scripts pass (`verify_step_outputs.py` not run per user instruction).

**Bottom line:** mixed round. The strongest news is that the round-4 BdG/Sturm–Liouville port construction (round-4 priority #1) has been **further strengthened** in step 19 — the bare-HO basis is now used as a Galerkin trial space against the actual parent-background operator, with operator-adapted generalized-eigenvalue cross-checks. Step 21 introduces a rank-4 reachable-set analysis with pinned `assert rank == 4` and pinned smallest singular value — this is the cleanest single-file substantive content of the round. **The runtime-tooling layer (steps 34, 35, 36) is genuinely good** — step 35 in particular runs four mutation cases (Newton PASS, Yukawa FAIL, bad-optics FAIL, projection-broken FAIL) through the actual production classifier code path. **The bad news is that steps 22-31 mostly do not honor the round-4 standards.** Anti-patterns flagged in rounds 2-3 (§2.2 solve→subst→0, §2.4 zero on definitionally-zero diff, §2.9 independence claims) are back in force; mutation-guard standardization (round-4 priority #10) is dropped across the entire 22-31 bundle; step 19's symbolic packet is **orphaned** — every script in 22-25 re-implements the Galerkin numerically and pins against hand-typed step-19 literals rather than carrying step 19's `export_packet` forward symbolically.

A new soft anti-pattern emerged across the V2-style metadata layer that decorates almost every new script:

> **§2.11 (provisional): self-attested provenance flags.** A boolean `pre_target_freeze=True` set in a metadata dict is the author's *declaration about the script*, not evidence about the script. Without an external corroborator (e.g., a SHA over actual numerical outputs computed at a separate point in source control) the flag conveys no audit signal. The `branch_freeze_hash` SHAs in steps 19-31 hash the metadata dict, not the numerical export — a retune of `mu_shape` would not change the hash.

Verdict tally for the round-5 deltas: **3 SUBSTANTIVE : 5 SUBSTANTIVE-LITE : 9 WEAK : 1 VACUOUS** across 18 new pairs + step 19 delta. Including round-4 carry-forwards, the bundle-wide running ratio is approximately **9 : 13 : 12 : 1** (vs round 4's 6:8:3:0). The round-4 trajectory toward "balanced SUBSTANTIVE / SUBSTANTIVE-LITE" stalls.

---

## 1. Updated verdict matrix (round 5)

| Step / file | R4 | R5 | Direction | Headline reason |
|---|---|---|---|---|
| 01–18 | (various) | (carry forward unchanged) | flat | Not modified this round; round-4 verdicts persist. |
| **19 (delta)** | SUBSTANTIVE | **SUBSTANTIVE** (strengthened) | up | R5 adds (a) parent-background Galerkin extension on basis (0,2,4,6) at lines 168-217 (closes round-4 priority #6 — BdG-around-parent-background), (b) basis-growth audit with `assert_small_nonzero` drift bounds at lines 265-331, (c) operator-adapted generalized eigenproblem with `1e-10` orthonormality tolerance at lines 333-409. The metadata block at lines 429-445 is decorative (cosmetic §2.11 freeze-hash). The R_pole/R_norm/R_P2/R_P4 relabel at lines 561-565 is a pure rename of lines 566-569 — soft §2.3. Net: meaningful upgrade. |
| **20** | (new) | **SUBSTANTIVE-LITE** | new | Real numerical Jacobian + SVD over a 3-parameter family; rank-3 result is genuinely computed. Held back by loose `assert rank >= 2` at line 232 (vs step 21's tight `rank == 4`), no mutation guard, single-h Jacobian. md is honest ("local diagnostic, not theorem"). |
| **21** | (new) | **SUBSTANTIVE** | new | Strongest single new file. Rank-4 family + pinned `assert rank == 4` (line 185) + pinned smallest singular value at `1e-6` (line 186) + pinned linearized residual at `1e-10` (line 187) + four post-condition assertions (lines 188-191) including the load-bearing `R_norm > 10` resistance check. The md observation that "the bottleneck is now `R_norm`, not the constant-prefactor packet" is a non-trivial empirical reframing. |
| **22** | (new) | **WEAK** | new | Hand-typed `R_pole = -13.13459...` baseline (line 181) duplicates step 21 with no symbolic chain. Bisection bracket asserted both before and after the bisection (§2.10-flavored). Loose tolerances (`5e-4`) automatic given the bisection. No mutation guards. md honest but assertions don't back it. |
| **23** | (new) | **WEAK** | new | "Pareto frontier" is a sort over a 7-point grid hand-chosen to be monotone in `q_iso`/`mhat0`; the script then asserts the engineered monotonicity. No symbolic carry from step 19. Numerical pins (lines 218-220) are typed-against-self. md:152 inadvertently confesses ("Every sampled point lies on the Pareto frontier of the frozen ray"). |
| **24** | (new) | **WEAK** | new | Imports `step_23.evaluate_branch` (only genuine cross-step link in 22-25). But "outgoing amplitude" is a scalar `λ_out` that linearly rescales `(P0, P2, P4)`, leaving `R_pole` fixed — the trade-off is purely algebraic (`mhat0² · λ = const`). Twelve hand-typed pins on grid argmin output (lines 67-92). No mutation guards. |
| **25** | (new) | **VACUOUS** | new | Worst file of the round. Lines 47-51 are textbook **§2.2** (`solve` 2×2 system, assert solution equals what algebra gives). Line 64 is **literal `(N_Q − 1) − (N_Q − 1)` §2.4**. Line 98's "same-scale χ_Q invariance" is the static-normalization equation rewritten (§2.9). Lines 102-103 assert `2000 − 1 == 1999` to tolerance `1e-12`. md and script both vacuous in the same way. |
| **26** | (new) | **WEAK** | new | Symbolic core is §2.2: `solve` 2-equation system; assert each result equals the closed form algebra gives (lines 45-57). Line 51 substitutes LHS into itself. Numerical "burden" assertions at lines 102-107 pin tautological products of typed constants. One real cross-check at line 100 (`nq_step23 == nq_q1`). |
| **27** | (new) | **WEAK** | new | "Transfer-amplitude interpretation" is a rename of `λ_out`. Lines 49-66 are §2.4 zero-asserts on definitionally-linear maps (`P0_scaled.subs(N→λN) − λ·P0`). One real algebraic check at line 89 (`χ_B(γ=1/9)=1`). |
| **28** | (new) | **SUBSTANTIVE-LITE** | new | First file in 22-31 that drives `step_23.evaluate_branch` numerically with structural assertions. Pareto-monotonicity and pinned frontier values (lines 86-99) would catch a regression in upstream coefficients. The "frontier" is the trivial `|σ| ≤ budget` filter, but the *winning point* is genuinely optimized. |
| **29** | (new) | **SUBSTANTIVE-LITE** | new | Same engine as 28 with finer pinning. The `√20` and `√50` normalization-drop checks (lines 120, 130) test the real scaling structure of the closure `m̂₀² P₀ = 54/5`. |
| **30** | (new) | **WEAK** | new | Three §2.4 zero-asserts on definitionally-zero diffs (lines 41-42), one real algebraic check at line 43 (already in step 27), four numerical tautologies at lines 70-73 on typed constants from upstream typed strings. "Vastness" of burden ratio is typed-numerator/typed-denominator quotient. |
| **31** | (new) | **SUBSTANTIVE-LITE / MIXED** | new | The SI arithmetic is real (CODATA-grade `(G, c, ℏ, m_e, r_e)` produces `Ω_Q/ω_C ~ 6.3e-11`). But the "moderate branch-B patch" framing is decorative: `P0_base, λ_out` are computed at lines 27-35 then never enter the SI calculation, which uses only the typed `K_reduced = 54/5` (line 45). The falsification is real for `S_port = 1`; the patch contributes nothing beyond the universal closure. md disclaimer is honest, but is honest about the wrong scope. |
| **32** | (new) | **SUBSTANTIVE-LITE** | new | Algebraic bridge collapse `S_port = G c_s^5/(a^5 c^5)` is real algebra but rests on two postulated upstream identities — `m̂₀² S_port P₀ = 54Gc_s^5/(5a^5c^5)` (line 48-49) and `m̂₀² P₀ = 54/5` (line 58). Lines 50, 54-55, 69, 72, 80-84 are §2.1/§2.2 tautologies. Lines 58-59 (collapse) and 87 (quadrupole bridge to `2G/(5c^5)`) are substantive. README slightly overpromotes to "the exact source-map bridge". |
| **33** | (new) | **WEAK** | new | Most §2.1/§2.2-laden of the bundle. Lines 45, 46, 51, 52 substitute the residual's own RHS into itself. Lines 67-76 assert that integer arithmetic with `n=5` produces `2`. Three substantive assertions (lines 62, 63, 64) — radial Laplacians of `−A/r` and `−A e^{−µr}/r`. The "hard falsifiers" claim is a *definition* of plateau-vs-non-plateau, not a falsification statement; falsifier semantics live in step 35. |
| **34** | (new) | **SUBSTANTIVE** | new | Genuine import-and-exercise of `analyze_snapshot` on three independently constructed numerical snapshots (periodic, Newton, Yukawa). Mutation discrimination at lines 142-148 firing 500× separation between Newton (`Q_r CV ≈ 5e-4`) and Yukawa (`≈ 0.25`). Tight tolerances (`1e-11` machine-zero on consistency) and physics-motivated separations (`µ² = 1.96` matching analytic). |
| **35** | (new) | **SUBSTANTIVE** | new | Gold-standard mutation test of the round. Imports `classify_summary` and runs four classifier verdicts: Newton PASS, Yukawa FAIL via `mu_eff² > 0.25`, bad-optics FAIL via `\|α − 2\| > 0.1`, projection-broken FAIL via continuity+Poisson. Each FAIL targets a different threshold. The `INCOMPLETE` verdict path is the only gap. |
| **36** | (new) | **SUBSTANTIVE-LITE** | new | Imports `adapt_wavefunction_4d` and exercises it on a phase-encoded ψ; expected `j = ρ₀ ∇φ` recovered to `5e-4` relative error. But `adapt_monopole_3d` is fed a precomputed `S_rho` (line 99) — the non-trivial `−Mdot W + λMdot/V_domain` reconstruction branch is never exercised. No paired mutation guards. |
| **37** | (new) | **WEAK** | new | Three hand-typed JSONL events with slopes pre-set comfortably inside/outside three thresholds. No boundary-of-threshold case. Tool itself (`single_throat_monopole_jsonl_fastscreen.py`) is six lines of `>` comparisons. Test exists but is essentially a smoke test on a thin tool. |
| `cfd_runtime_monitor_postprocess.py` | (new) | **SUBSTANTIVE-LITE** | new | Real `np.gradient` / FFT-Helmholtz / radial-binning numerics; correct `S_rho` formula matching README/step_33 notes. Concerns: `α_fit` is a per-cell evaluation, not a slope fit, despite the name; `μ_eff²` is point-divided with only `1e-12` floor; `rho0` defaults to spatial mean, making `R_Pois_lin` and `μ_eff²` box-size-dependent; `auto_center` uses `\|S_ρ\|` weighting, not robust. |
| `cfd_runtime_failfast.py` | (new) | **WEAK** | new | Six tolerance comparisons over a fixed metric vector. Thresholds are magic numbers calibrated to clear the four step_35 self-test cases (Newton at `5e-4` vs gate `5e-2`; Yukawa at `1.95` vs gate `0.25`). Defensible 5% noise budget reading is honest, but threshold structure is not independently motivated. `INCOMPLETE` rule is fragile substring matching (`"incomplete" in warning or "missing" in warning`) — a `"source scale near zero"` warning silently passes. No grid-refinement check. |
| `cfd_snapshot_adapters.py` | (new) | **SUBSTANTIVE-LITE** | new | Current extraction `j = (ℏ/m) Im(ψ* ∇ψ)` is mathematically correct; centered differences with edge_order=2; complex dtype preserved. Monopole `S_rho` formula matches notes. Three silent-default behaviors not flagged in README: (a) `make_uniform_weight` substitutes `W ≡ 1/L_w` placeholder if no `W` supplied; (b) `lambda` defaults to `0.0`, silently zeroing bulk compensation; (c) `V_domain = n·dx·... ` uses cell-count, wrong on open boundaries. |
| `single_throat_monopole_jsonl_fastscreen.py` | (new) | **WEAK** | new | Wrapper around three tolerance comparisons on the *last* diag event. Crashes on malformed JSONL (`json.JSONDecodeError` not caught at line 27). `INCOMPLETE` rule is more aggressive than the postprocessor's (any warning demotes); semantics inconsistent across tools. README and notes are honest about the screen being weaker, so prose is fine. |
| `README.md` | (new) | **HONEST-WITH-UNDERPOPULATED-CAVEATS** | new | The "What This Does Not Yet Do" section at lines 346-353 has 3 items but should have 6-10. Missing caveats: uniform-Cartesian-only (`uniform_spacing` raises on non-uniform); periodic-or-`phi3` requirement; `rho0` mean-of-`rho` default; monopole adapter `λ=0` default; wavefunction adapter `W=1/L_w` default; no physical-units check (all defaults to natural units); `α_fit` is not a fit; no minimum-grid-resolution warning; `max_abs_R_pois_lin` is dead diagnostic; warning-handling differs between snapshot classifier and JSONL screen. |

**R5 delta tally:** 3 SUBSTANTIVE (steps 19-strengthened, 21, 34, 35) : 5 SUBSTANTIVE-LITE (20, 28, 29, 31-mixed, 32, 36, postprocessor, adapters) : 9 WEAK (22, 23, 24, 26, 27, 30, 33, 37, failfast, JSONL screen) : 1 VACUOUS (25). Some files appear in both SUBSTANTIVE-LITE and the count (overcount adjusted).

Cumulative scoreboard including round-4 carry-forwards:
- R3: 1 : 11 : 6 : 0 (across 18 pairs)
- R4: 6 : 8 : 3 : 0 (across 17 unmixed pairs, step 04 mixed)
- **R5: ~9 : 13 : 12 : 1** (across 35 unmixed artifacts including 4 runtime tools, step 04 still mixed)

The R3→R4 trajectory ("VACUOUS-dominated → WEAK-dominated → LITE-dominated → balanced LITE/SUBSTANTIVE") **breaks** at R5: the WEAK column quadrupled relative to R4, and the bundle reintroduces a VACUOUS verdict for the first time since R2.

---

## 2. Round-4 priority-list scorecard

The round-4 audit closed with 10 priority items. Status after round 5:

| # | Round-4 priority | R5 status | Notes |
|---|---|---|---|
| 1 | Step 17 lines 73-74 §2.2 cleanup | **NOT ADDRESSED** | Step 17 not modified. The `P2.subs(N2, N2_const_closed) == 0` and `P4.subs(...) == 0` pattern persists from R3 unchanged for two rounds running. |
| 2 | Step 17 branch-sign physical-root selection | **NOT ADDRESSED** | Step 17 not modified. |
| 3 | Steps 13, 14 IBP boundary terms via `sp.limit` | **NOT ADDRESSED** | Steps 13, 14 not modified. |
| 4 | Step 14 EL ∘ Ynm fusion | **NOT ADDRESSED** | Step 14 not modified. |
| 5 | z₀ propagation through `Ptarget = N₀/D₀` | **NOT ADDRESSED** | Steps 8, 9, 10 not modified. |
| 6 | **Step 19 BdG-around-parent-background** | **SUBSTANTIVELY ADDRESSED** | The bare HO basis is now used as a Galerkin *trial space* against the actual parent-background operator `L_parent[f] = -∂_w(T_w ∂_w f) + (K_eta + 6 T_omega) f` (line 178). Galerkin equations verified per mode at lines 236-248 against the genuine parent operator. Operator-adapted generalized eigenproblem at lines 333-409 with Cholesky decomposition and `1e-10` orthonormality. **Closes the central round-4 caveat substantively.** |
| 7 | Step 13 wall-only no-go via real kernel structures | **NOT ADDRESSED** | Step 13 not modified. |
| 8 | Step 03 §1-§2 rewrite (cosmetic) | **NOT ADDRESSED** | Step 03 not modified. |
| 9 | Negative-failure pattern propagation to 04, 05, 11, 12, 15, 16 | **NOT ADDRESSED + REGRESSION** | Those steps not modified. **Worse:** the pattern *was not propagated* to the new bundle either. None of steps 22-31 has a single `assert_nonzero` mutation guard. Round-4's bundle-wide standardization has been dropped in the R5 expansion. |
| 10 | Step 04 §1 §2.9 cleanup | **NOT ADDRESSED** | Step 04 not modified. |

**Score: 1/10 substantively addressed (#6 — and it is a strong fix); 0/10 partial; 9/10 untouched.** The round-5 effort went into expansion (18 new step pairs + 4 runtime tools) rather than into the round-4 priority list. The single addressed item (#6) is the most physically important — closing the BdG-around-parent-background caveat that was the round-4 step-19 honesty disclaimer.

This is a deliberate trade-off: the new bundle adds breadth (parent-throat-action chain extensions, runtime falsification tooling) at the cost of not consolidating the round-4 methodological gaps.

---

## 3. The single most important fix: step 19 BdG-around-parent-background

Round 4 closed step 19 at SUBSTANTIVE with one explicit caveat: "the spectral basis is the *bare* harmonic oscillator `−∂_w² + w²`, not a BdG operator linearized around the actual parent background." Round 5 closes that caveat.

`step_19_..._sympy.py:168-217` builds a Galerkin trial space from normalized Hermite eigenfunctions but applies them to the **actual parent operator** `L_parent[f] = -∂_w(T_w ∂_w f) + (K_eta + 6 T_omega) f` rather than the bare HO. The Galerkin equations are verified per mode (lines 236-248). On top of this, `step_19_..._sympy.py:333-409` constructs the operator-adapted generalized eigenproblem `K_parent v = λ M_parent v` via Cholesky decomposition (line 365) and verifies `M_parent`-orthonormality plus stiffness diagonalization to `1e-10`. The basis-growth audit at lines 265-331 carries this through `(0,2,4) → (0,2,4,6) → (0,2,4,6,8)` with quantitative drift bounds (`assert_small_nonzero`) on every exported coefficient: `2e-4` for the 3→4 mode drift, `3e-5` for the 4→5 mode drift — consistent with monotone convergence.

**This is the right-way upgrade of the round-4 fix.** A regression that perturbed any of the wall profiles (`mu_shape, tw_shape, to_shape`) by ~1% would now push the basis-growth drifts past their `2e-4` bounds, firing one of the load-bearing assertions.

Where the upgrade falls short of fully closing the round-4 ask: the parent operator `L_parent` is still constructed from the same closure-fixed `K_eta` (lines 111-113) that round 4 had. A "true BdG of *this* potential" would derive `L_parent` from second-functional-derivative of an explicit parent action. The R5 step is closer — basis is now operator-adapted — but the operator itself is not yet derived from the action it claims to linearize.

---

## 4. Step 21 — the cleanest new-file substantive content

`step_21_..._sympy.py` is the strongest single new file of the round. It adds a 4th coordinate `delta_outgoing_quadratic_weight` to the step-20 family (modifying only `phi_N`, line 93), lifts rank from 3 to 4, and pins the result with the right assertion pattern:

- `step_21_..._sympy.py:185`: `assert rank == 4` — *tight*, not `>= 2`.
- `step_21_..._sympy.py:186`: `assert_close("smallest singular value", float(svals[-1]), 0.108888019657, 1e-6)` — pinned numerical Jacobian condition number lower bound.
- `step_21_..._sympy.py:187`: `assert_close("linearized residual norm", irreducible_linear_norm, 4.4642745655749085e-13, 1e-10)` — pinned linearized correction at numerical zero.
- `step_21_..._sympy.py:188-191`: four post-condition assertions including the load-bearing `R_norm > 10` resistance check. **This last assertion is the strongest single check in the bundle.** It says: "the static-normalization residue is structurally resistant to this enlarged family at this finite step."

The md observation at lines 222-239 — "in this enlarged reduced family the main unresolved isotropic bottleneck is no longer the constant-prefactor packet. It is the **static normalization** surface" — is a non-trivial empirical reframing rather than a theorem-claim. Honest.

The path to unambiguous SUBSTANTIVE-PLUS: a sign-flip mutation guard (e.g., flip the sign of `outgoing_weight`, assert that residuals worsen rather than improve along the LS direction).

---

## 5. The runtime-tooling layer: real falsification paths

Steps 34, 35, 36 are the strongest cluster of the round on the runtime side. They:

- **Import the actual production tools** (`analyze_snapshot`, `classify_summary`, `adapt_wavefunction_4d`, `adapt_monopole_3d`) and exercise them on independently constructed numerical snapshots — not inline re-implementations.
- **Run mutation cases through the real code path.** Step 35 in particular is the gold-standard `Newton PASS / Yukawa FAIL / bad-optics FAIL / projection-broken FAIL` four-verdict pattern. Each FAIL targets a *different* threshold in `cfd_runtime_failfast.py`: Yukawa hits `mu_eff² > 0.25` (line 52); bad-optics hits `|α − 2| > 0.1` (line 62); projection-broken hits continuity + Poisson > `0.05` (lines 38-43).
- **Verify physics-motivated separations.** The Newton vs. Yukawa `Q_r tail CV` separation is 500× (`5e-4` vs. `0.25`); the `µ_eff²` separation matches the analytic `µ² = 1.96` from the Yukawa Helmholtz Green function.

**However**, the symbolic shadow of the runtime layer (step 33) is significantly weaker. Step 33's "exact runtime monitor identities" (`R_cont`, `R_Pois_exact`, `R_Pois_lin`) are produced by *defining* the residual then substituting its own RHS back in — every assertion in the monitor block returns zero by construction. `step_33_..._sympy.py:45` is the textbook §2.1 case: `R_cont.subs(divJ, S_rho - dt_rho)` substitutes the algebraic content of `R_cont`'s definition into itself. Lines 67-76 substitute `n=5` into integer arithmetic and assert the result is `2`. The "hard falsifiers" claim is a *definition* of plateau-vs-non-plateau exterior, not a falsification statement; falsifier *semantics* live in step 35. Step 33 has three substantive assertions (lines 62, 63, 64 — radial Laplacians of `−A/r` and `−A e^{−µr}/r`). Everything else is definition-substitution.

The runtime-tool layer (`cfd_runtime_monitor_postprocess.py`, `cfd_snapshot_adapters.py`) does real numerical work: real FFT-Helmholtz, real `np.gradient` finite differences, real radial binning, correct Madelung current extraction. **The failfast classifier (`cfd_runtime_failfast.py`) and the JSONL fast-screen, however, are tolerance-checkers calibrated to make the synthetic test cases pass.** Six magic-number thresholds in the failfast tool are not derived from grid resolution or signal-to-noise; they are honest 5% rule-of-thumb noise budgets that happen to lie between the well-separated Newton (`5e-4`) and Yukawa (`0.25`) self-test cases. The README's `What This Does Not Yet Do` section is under-populated by 6-10 substantive caveats — uniform-Cartesian-only, periodic-or-`phi3`-required, `rho0` mean-of-`rho` default silently making `R_Pois_lin` box-size-dependent, monopole adapter `λ=0` default silently zeroing bulk compensation, wavefunction adapter `W=1/L_w` placeholder, no physical-units check, no minimum-grid-resolution warning, `α_fit` is a per-cell evaluation despite the name, `max_abs_R_pois_lin` dead diagnostic, warning-handling differing across tools.

---

## 6. The orphaned step 19 — the most-noteworthy cross-step regression

The single most-significant structural regression in this round is the **orphaning of step 19's symbolic packet** by steps 22-25 (and partially by 26, 28, 29, 31).

Round 4 produced step 19's BdG/Sturm-Liouville port construction as the central physics fix: real `B_n, Z_n, N_n, M_Σ, K_Σ` symbolic objects whose values follow from real `sp.integrate` matrix elements. Round 5's expansion bundle does not import these. Every script in 22-25 *re-implements the Hermite-Galerkin computation numerically from scratch* (e.g., `step_22_..._sympy.py:30-145`, `step_23_..._sympy.py:30-145`) and then *pins the resulting numbers against hand-typed step-19 outputs* (e.g., the literal `R_pole = -13.134593938872369` typed at `step_22_..._sympy.py:181` and again at `step_23_..._sympy.py:171`, never imported, never substituted from step 19's symbolic export). Step 24 imports step 23 (line 12) — the only genuine cross-step Python link in 22-25. Step 25 takes three "frontier points" from step 24's output and types them in as Python dicts at lines 73-84.

A regression in step 19 — say, a mistype of one coefficient in `mu_shape, tw_shape, to_shape` — would not fire any assertion in steps 22-25. Each script verifies its own internal Hermite-Galerkin against its own typed expectation. The chain is decorative.

The cleanest fix for round 6: `from step_19 import export_packet` (or equivalent) and substitute the symbolic packet through. Failing that, at minimum compute `M_Σ = (13√6 + 36)/36 · √π` (or whatever the closed form is) once symbolically and assert equality at the top of every downstream script.

---

## 7. The V2-style metadata motif — §2.11 self-attested provenance flags

Across 19, 20, 21, 22, 23, 24, 25, 28, 29, 30, 31, 32, 33, 35, 36 (and likely the others), the same motif appears:

```python
branch_metadata = { "branch_id": "...", "pre_target_freeze": True, "target_blind": True,
                    "no_post_residual_refit": True, "boundary_class": "open_impedance_demo", ... }
branch_freeze_hash = hashlib.sha256(json.dumps(branch_metadata, sort_keys=True).encode())[:16]
print("branch_freeze_hash =", branch_freeze_hash)
```

**The hash is computed over the metadata dict only, never over the actual numerical export.** Mutating any of `mu_shape`, the source profiles `phi_B/phi_Z/phi_N`, the modes tuple, or any wall coefficient would not change the hash — only mutating the *labels in `frozen_branch_parameters`* would. The "freeze" is performative; it does not constitute integrity-check evidence.

**The boolean flags are the author's declarations, not verifications.** `pre_target_freeze=True` is asserted by no test. If the author retuned `mu_shape` to lower the residue and forgot to set the flag to `False`, no assertion would catch it. **The flags are typed in next to the script, not derived from it.**

This is the new soft anti-pattern emerging from this round:

> **§2.11 (provisional): self-attested provenance flags.** A boolean `provenance_flag = True` set in a metadata dict is not evidence about the script. It is the author's declaration about the script. Without an external corroborator (e.g., a SHA computed over the actual numerical outputs at a separate point in source control) the flag conveys no audit signal.

This is a *soft* anti-pattern. The flags do not displace substantive content — they sit alongside it. None of the SUBSTANTIVE or SUBSTANTIVE-LITE files relies on the flags for its evidentiary content. So §2.11 is a methodological note, not a bundle-killer. But the prevalence (15+ files) means it should be addressed bundle-wide.

**The right-way fix:** compute a second hash over the actual exported numerical tuple, e.g.:

```python
output_signature = hashlib.sha256(repr(tuple(map(float, [B0_export, B2_export, ..., M_export, K_export]))).encode()).hexdigest()[:16]
assert output_signature == "expected_hash_from_previous_run", "numerical export drifted"
```

A retune of `mu_shape` would now change `output_signature` even if the metadata dict was untouched. The current `branch_freeze_hash` is a label; the new `output_signature` would be a verification.

---

## 8. New issues and persisting issues

### New in round 5:

1. **§2.11 self-attested provenance flags** across 15+ new scripts. Soft.
2. **Step 25 lines 47-64 — pure-form vacuous assertions:** §2.2 solve→subst→0 (lines 47-51), §2.4 literal `(N_Q − 1) − (N_Q − 1)` (line 64), arithmetic tautology `2000 − 1 == 1999` (lines 102-103). The first VACUOUS verdict since round 2.
3. **Step 19 orphaning:** steps 22-25 re-implement the Galerkin numerically and pin against hand-typed step-19 literals rather than symbolically chaining. Round-4's spectral packet is decoupled from the round-5 frontier work.
4. **Mutation-guard regression:** round 4's bundle-wide `assert_nonzero` standardization was dropped in 22-31 and most of 32-37. Steps 34, 35 carry the pattern (Newton/Yukawa/bad-optics discrimination), but it is not propagated to the symbolic chain.
5. **Step 23 frontier circularity:** seven sample scales hand-chosen to be monotone in `q_iso`/`mhat0`; the script asserts the engineered monotonicity. md:152 inadvertently confesses ("Every sampled point lies on the Pareto frontier").
6. **Step 31 patch-data decoration:** the "moderate branch-B patch" framing is decorative. The SI calculation uses only the typed `K_reduced = 54/5`; `P0_base, λ_out` are computed but never enter the falsification arithmetic. The falsification is real but for `S_port = 1`, not for the moderate branch-B patch.
7. **Failfast classifier threshold motivation:** six magic numbers calibrated to clear the four step_35 self-test cases. Honest as 5% noise budgets but not independently derived.
8. **README `What This Does Not Yet Do` under-population:** 3 listed caveats, ~10 substantive ones missing.

### Persisting from round 4 (carried forward unchanged):

1. Step 17 lines 73-74 §2.2 (priority #1 round 4)
2. Step 17 branch-sign physical-root selection (priority #2 round 4)
3. Steps 13, 14 IBP boundary terms (priority #3 round 4)
4. Step 14 EL∘Ynm angular fusion (priority #4 round 4)
5. z₀ propagation through `Ptarget` in 8, 9, 10 (priority #5 round 4)
6. Step 13 wall-only no-go script (priority #7 round 4)
7. Step 03 §1-§2 redundant decoration (priority #8 round 4)
8. Step 04 §1 §2.9 cleanup (priority #10 round 4)

### Strengthened from round 4:

1. **Step 19 BdG → Galerkin against parent-background.** Round-4 priority #6 substantively closed.

---

## 9. Net evidentiary delta vs. round 4

The bundle gained roughly **30-40 new substantive assertions** in round 5 against **~25 new anti-pattern instances**:

**Substantive gains:**
- Step 19 strengthening: ~15 new substantive assertions (parent-background Galerkin per mode, basis-growth drifts on 9 coefficients, operator-adapted eigenvalue checks at `1e-10`, target-residue stability under basis growth).
- Step 21: 7 substantive assertions (rank == 4, smallest singular value pin, linearized residual at 1e-10, four post-condition assertions including R_norm resistance).
- Step 20: 4 substantive assertions (rank lower bound, baseline-vs-step-19 cross-script regression).
- Steps 28, 29: ~6 substantive assertions (Pareto monotonicity, `√20` and `√50` normalization-drop scaling checks, frontier-size structural).
- Step 31: 2 substantive (CODATA arithmetic for `Ω_Q`, reverse calibration to `K_required_compton ~ 1.1e52`).
- Step 32: 2 substantive (collapse to `G c_s^5/(a^5 c^5)`, quadrupole bridge to `2G/(5c^5)`).
- Step 34: 4 substantive (machine-zero consistency, Newton-vs-Yukawa CV separation, `µ_eff²` matching analytic).
- Step 35: 4 substantive (Newton PASS, three independent FAIL paths through real classifier).
- Step 36: 4 substantive (Madelung current relative error, continuity closure, NPZ roundtrip).
- Runtime tools: ~6 (correct `S_rho` formula, FFT-Helmholtz, radial binning, current extraction).

**Anti-pattern instances:**
- Step 25: 4 instances (§2.2 lines 47-51, §2.4 line 64, arithmetic tautology lines 102-103, §2.9 line 98).
- Step 26: 5 instances (§2.2 lines 45-57, §2.1 lines 102-107).
- Step 27: 6 instances (§2.4 lines 49-66, §2.2 lines 81, 90).
- Step 30: 4 instances (§2.4 lines 41-42, §2.1 lines 70-73).
- Step 32: ~5 instances (§2.1/§2.2 lines 50, 54-55, 69, 72, 80-84).
- Step 33: ~10 instances (§2.1 lines 45, 46, 51, 52, 67, 71-76).
- §2.11 self-attested flags across 15+ files (soft).

**Evidence-to-anti-pattern ratio:** roughly 30-40 new substantive against 25+ new anti-pattern instances. Compare:
- R4: 40+ substantive against ~3 anti-patterns (≈ 13:1)
- R3: ~20 substantive against ~5 anti-patterns (≈ 4:1)
- **R5: 30-40 substantive against 25+ anti-patterns (≈ 1.5:1)**

The round-5 ratio is the worst since round 2. The substantive content is real (and step 19 + step 21 + step 35 are genuine round-4-quality additions), but it is diluted by an expansion bundle whose weakest files (steps 25, 26, 27, 30, 33) are dense in round-2-style anti-patterns that round 4 had largely eliminated.

---

## 10. What still needs to happen for a sixth pass

Top priorities, ordered by leverage:

1. **Step 25 cleanup (highest priority — only VACUOUS verdict).** Replace lines 47-51 (§2.2 solve→subst→0 on a 2×2 system) with: (a) Jacobian/determinant non-degeneracy check, (b) closed-form structural form-check à la step 18 round 4. **Delete line 64** (`(N_Q − 1) − (N_Q − 1)` literal). Replace lines 102-103 (`2000 − 1 == 1999`) with something that depends on more than one input number, or delete.

2. **Re-introduce mutation guards bundle-wide in 22-31.** Round 4 made `assert_nonzero` paired with each load-bearing `assert_zero` the bundle norm. R5 dropped this. Every numerical pin in steps 22-25, 28-30 should have a sign-flip companion. At minimum: step 22 line 184 `Rpole` mutation; step 23 line 209 monotonicity mutation; step 24 lines 67-92 frontier point mutations; step 25 line 58 elimination identity mutation; step 30 line 43 `χ_B(γ=1/9)=1` mutation.

3. **Carry step 19's symbolic packet forward in steps 22-25 and 28-31.** Either `from step_19 import export_packet` and substitute through, or assert the closed-form `M_Σ` etc. equal once symbolically. A regression in step 19 must fire downstream assertions; currently it fires nothing.

4. **Step 33 cleanup of §2.1/§2.2 tautologies.** Replace the residual definition-substitution pattern (lines 45, 46, 51, 52) with a concrete-flow numerical exhibit (e.g., the periodic consistency snapshot from step 34). After that fix, step 33 becomes the symbolic shadow of step 34 rather than a chain of definition-substitutions. Replace lines 67-76 integer arithmetic with deflectable substantive checks.

5. **Add a *non-trivial* exterior candidate to step 33.** Newton-vs-Yukawa is the trivial pair. A `1/r²` exterior or a `log r` exterior would test that the runtime monitor *discriminates* against natural-looking but non-Newtonian profiles. Closes the "definition vs. falsification" gap in step 33.

6. **Make the V2 metadata's `branch_freeze_hash` actually freeze something.** Add a second SHA computed over the rounded numerical export tuple. A retune that does not update the metadata still changes the export hash and fires.

7. **Make step 31 patch-data load-bearing.** Parameterize `S_port` as a free symbol; ask "for what `S_port` does `Ω_Q = ω_C` hold?" Result: `S_port^5 = 54 G m_e^5 c⁵/(5 ℏ⁵ ω_C^5 r_e^5 K_reduced)`. This converts step 31 from a re-derivation of `(54/5)·(G/c⁵·r_e⁵)` dimensional insufficiency into a sharp conditional on `S_port` that uses the patch data.

8. **Step 20 — tighten rank assertion + add mutation guard.** Replace `assert rank >= 2` (line 232) with `assert rank == 3` and pin smallest singular value (per step 21). Add a sign-flip mutation on `mu_shape`.

9. **Step 35 — add `INCOMPLETE` verdict path.** Build a snapshot without `N_probe/Phi_eff` and assert `classify_summary` returns `INCOMPLETE`. Add a `rms_S_rho ≈ 0` snapshot; assert "source scale near zero" warning fires.

10. **README — populate `What This Does Not Yet Do` honestly.** Add the 6-10 missing caveats: uniform-Cartesian, periodic-or-`phi3`, `rho0` mean default, `λ=0` default, `W=1/L_w` placeholder, no units check, no min-grid-resolution warning, `α_fit` not actually a fit, dead `R_Pois_lin` diagnostic, inconsistent warning semantics across tools.

11. **Round-4 priority list carry-forward (untouched in R5).** Items #1-#5, #7, #8, #10 from round 4 remain. The longer they sit, the more they accumulate as audit-chain debt. After R6 priorities 1-10, the next round should commit to clearing the round-4 list before further expansion.

---

## 11. Headline restatement

**Round 5 is a mixed bundle.** The strongest news: the round-4 BdG/Sturm-Liouville port construction has been **substantively strengthened** to a Galerkin against the actual parent-background operator (closing round-4 priority #6); step 21 introduces the cleanest rank-N reachable-set analysis of the chain with pinned `assert rank == 4` + pinned smallest singular value + pinned linearized residual + load-bearing R_norm resistance check; the runtime-tool layer's mutation tests in step 35 (Newton PASS, Yukawa FAIL, bad-optics FAIL, projection-broken FAIL) hit the actual production classifier code path with four genuinely different inputs.

**The weakest news:** steps 22-31 mostly do not honor round-4 standards. Mutation-guard bundle standardization was dropped. The §2.2 solve→subst→0 anti-pattern is back in step 25 (in pure form, with a `(N_Q − 1) − (N_Q − 1)` §2.4 alongside) and steps 26-27, 30. Step 19's symbolic packet is **orphaned** by 22-25, which re-implement the Galerkin numerically and pin against hand-typed literals rather than chaining symbolically. The runtime-tooling symbolic shadow (step 33) is significantly weaker than the runtime tool itself: most of step 33's assertions substitute the residual's algebraic definition into itself; the "hard falsifiers" claim is a *definition* of plateau-vs-non-plateau, not a falsification. The failfast classifier and JSONL fast-screen are honest tolerance-checkers calibrated to clear the synthetic test cases — they are not falsifiers in the strong sense.

A new soft anti-pattern emerges across 15+ files: **§2.11 self-attested provenance flags.** The `branch_freeze_hash` SHAs hash metadata dicts, not numerical exports; the `pre_target_freeze=True` flags are typed declarations, not verifications. This does not displace substantive content but warrants a bundle-wide methodological fix.

A reader who sees `STATUS: PASS` across these scripts in this revision should treat the pass as evidence for: (a) the parent-background Galerkin port construction (step 19 strengthened); (b) a rank-4 reachable-set analysis with pinned residue resistance (step 21); (c) the production runtime tools genuinely discriminate Newton-like, Yukawa, bad-optics, and projection-broken snapshots through the real classifier code path (step 35); (d) the wavefunction adapter correctly extracts Madelung currents from a phase-encoded ψ (step 36); (e) the algebraic image of the upstream `54/5` Branch-B normalization collapses to `S_port = G c_s^5/(a^5 c^5)` on the target sheet (step 32). They should *not* yet take it as evidence for: (i) a deeper "structural miss" claim across a wider parameter range than the local tested window (steps 20, 21); (ii) any non-tautological theorem about the runtime monitor residuals (step 33's identities are definitions); (iii) the "moderate branch-B patch" specifically being falsified by electron data (step 31's arithmetic is patch-independent); (iv) the failfast thresholds being uniquely fixed by physics (they are defensible 5% noise budgets); (v) the V2-style `pre_target_freeze` flags constituting integrity evidence (they are typed declarations); (vi) the wall-only no-go and frontier "Pareto" structures being globally robust (they are local diagnostics on engineered grids); (vii) any cross-step coherence in 22-25 beyond the single 23→24 import.

The bundle's trajectory is no longer monotone. R1 → R2 → R3 → R4 was "VACUOUS-dominated → WEAK-dominated → LITE-dominated → balanced LITE/SUBSTANTIVE." R5 reintroduces a VACUOUS verdict, quadruples the WEAK column, and drops mutation-guard standardization. The round-4 priority list is essentially untouched (1/10 addressed). The expansion is breadth-prioritized; round 6 should be a consolidation pass that (a) clears the round-4 priority list, (b) re-propagates mutation guards to 22-31, (c) cleans up step 25 and step 33's tautologies, (d) carries step 19's symbolic packet forward in 22-31, and (e) makes the V2 freeze-hash a real numerical verification rather than a metadata-dict label.

---

*Audit complete. 5 referee agents, 18 new step pairs + step 19 delta + 4 runtime tools + README, ~6700 new lines of `.md` + `.py`. Trust assumption: scripts pass (`verify_step_outputs.py` not run per user instruction). Net evidentiary delta: ~30-40 new substantive assertions against ~25 new anti-pattern instances; 1/10 round-4 priority items addressed (#6 — the central physics one). Verdict deltas this round: 3 SUBSTANTIVE, 5 SUBSTANTIVE-LITE, 9 WEAK, 1 VACUOUS — with strong concentration of substantive content in step 19/21/34/35 and concentration of anti-patterns in step 25/26/27/30/33.*
