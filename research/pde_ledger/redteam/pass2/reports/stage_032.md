---
unit_id: 032
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage032_source_map_from_mode_integrals.md]
  paper_appendix: present
---

# Audit unit 032 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_032.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage032_source_map_from_mode_integrals.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 54)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.txt`

## What the paper claims

Stage 032 removes the abstract selected-branch source-map factor `m̂_-` by showing that on the natural D/N source branch it is fixed entirely by the same wall–D/N overlap that already controls the loading. The card's boxed `\stagefield{Output}` deliverables are: (1) the source-map factor `\widehat m_- = (v·e_-)/\kappa_0`, hence `\widehat m_-^2 = s_-/\kappa_0^2` (eq:app-stage032-mhat-minus, lines 36–41); (2) the exact bound `1 ≤ \widehat m_-^2 < σ/\kappa_0^2 = 11/9` (eq:app-stage032-mhat-bound, lines 44–47); and (3) the invariant selected product `N_-(α_0) = \widehat m_-^2 P_{0,-} = β_0 s_-^2/(\kappa_0^2 λ_-)` set equal to `N_Q^{target}` (eq:app-stage032-Nminus / -Nminus-target, lines 54–63). The `\stagefield{Checks}` enumerate the structural supports: same overlap vector `v=(\kappa_0,\kappa_1)^T`; denominator `\kappa_0` at zero loading (`\widehat m_-=1` at onset); and `σ/\kappa_0^2 = 1+2/9 = 11/9`. The notes add the explicit closed values that the terse card omits: `\kappa_0 = 2√2/π`, `\kappa_1 = -4/(3π)`, `σ = 88/(9π²)`, `\kappa_0^2 = 8/π²`, the local-kernel reductions `λ_U=g_U` etc., the Schur structure `Σ_wall = Ξ I_2 + α v v^T`, and the eventual 2.5PN target literal `54 G c_s^5/(5 a^5 c^5)` (notes §6, line 252) — none of which the scripts compute.

## What the script claims to verify

Both engines verify five sections that mirror the card/notes. §1: the exact axial mode integrals — orthonormality of the N/N basis and D/N half-wave, and the overlaps `\kappa_0 = 2√2/π`, `\kappa_1 = -4/(3π)`, `σ = 88/(9π²)`, `σ/\kappa_0^2 = 11/9`. §2: the local-kernel mode reductions, confirming each coupling collapses to the same overlap vector `v` (`L_src = g_Q Q (v·q)` etc.). §3: the Schur-complement decomposition `Σ = Ξ I + α vv^T` with `Ξ = g_U^2/A_U` and the stated `α` coefficient (Mathematica also re-derives the Schur image via an independent `LinearSolve` route and cross-checks against the direct inversion). §4: the natural D/N source map — the closed forms `λ_-`, `s_-`, and `\widehat m_-^2`, independently cross-checked against an actual eigensolve of the loaded operator `M(α_0) = diag(A, A+DK) - α_0 vv^T`, plus the endpoint checks `\widehat m_-^2(α_0=0)=1` and `lim_{α_0→∞}\widehat m_-^2 = 11/9`. §5: the eliminated product `N_-` reproduces `β_0 \kappa_0^2/A` at `α_0=0` and vanishes as `α_0→∞`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `\widehat m_-^2 = s_-/\kappa_0^2` (lines 36–41) | sympy 163 `mhat_sq = s_minus_nat/kappa0**2` + independent eigensolve check 187–189; wl 148, 167–170 | match |
| Bound `1 ≤ \widehat m_-^2 < 11/9` (lines 44–47) | sympy 195–196 (`=1` at 0, `→11/9` at ∞); wl 175–176 | match |
| `N_- = β_0 s_-^2/(\kappa_0^2 λ_-)` (lines 54–57) | sympy 203–215 (limits 0 and ∞); wl 179–196 (+ independent-eigenvector path 181–188) | match |
| `v=(\kappa_0,\kappa_1)^T` same overlap (check 1) | sympy 91–110 reductions; wl 81–85 | match |
| `\widehat m_-=1` at onset (check 2) | sympy 195; wl 175 | match |
| `σ/\kappa_0^2 = 11/9` (check 3) | sympy 65; wl 58 | match |
| `\kappa_0=2√2/π`, `\kappa_1=-4/(3π)`, `σ=88/(9π²)` (notes 67–75) | sympy 62–64; wl 55–57 | match |
| `N_-(α_0) = N_Q^{target}`, `=54 G c_s^5/(5 a^5 c^5)` (tex 62; notes 252) | not computed by either script | n/a (out of script scope; see note) |

The 2.5PN target equality is a downstream identification the card writes only as `N_Q^{target}`; the explicit `54 G c_s^5/(5 a^5 c^5)` literal lives only in the notes (§6) and is deferred to Stage 033 ("Stage~033 uses `N_-` to write the microscopic target"). Neither script emits any of `54, G, c_s, a, c`, so it is not a script-emitted deliverable for this unit — no reconciliation obligation. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58–61 | `expect_zero` basis orthonormality | structural (mode integrals) | yes |
| A2 | sympy | 62–65 | `expect_zero` κ_0, κ_1, σ, σ/κ_0²=11/9 | notes overlaps + check 3 | yes |
| A3 | sympy | 91–110 | `expect_zero` 5 local-kernel reductions | check 1 (same `v`) | yes |
| A4 | sympy | 143 | `expect_zero` `Σ - [Ξ I + α vv^T]` | notes §3 Schur structure | yes |
| A5 | sympy | 187–189 | `expect_zero` `s_check - s_minus_nat` (eigensolve vs closed form) | `\widehat m_-^2 = s_-/\kappa_0^2` (deliverable 1) | yes |
| A6 | sympy | 191–194 | `expect_zero` `λ_- ` eigensolve vs closed form | supports deliverable 1/3 | yes |
| A7 | sympy | 195 | `expect_zero` `\widehat m_-^2(0) - 1` | bound lower edge | yes |
| A8 | sympy | 196 | `expect_zero` `lim_∞ \widehat m_-^2 - 11/9` | bound upper edge (deliverable 2) | yes |
| A9 | sympy | 208–211 | `expect_zero` `Nprod(0) - β_0 κ_0²/A` | deliverable 3 endpoint | yes |
| A10 | sympy | 215 | `expect_zero` `lim_∞ Nprod` | deliverable 3 endpoint | yes |
| B1 | wl | 51–58 | `expectZero` basis + overlaps | = A1/A2 | yes |
| B2 | wl | 81–85 | `expectZero` reductions | = A3 | yes |
| B3 | wl | 121 | `expectZero` Schur via Inverse vs LinearSolve | independent Schur route | yes |
| B4 | wl | 134–135 | `expectZero` `Σ - [Ξ I + α vv^T]` (both routes) | = A4 | yes |
| B5 | wl | 167–170 | `expectZero` `s_check - s_minus_nat` | = A5 | yes |
| B6 | wl | 171–174 | `expectZero` `λ_-` eigensolve vs closed | = A6 | yes |
| B7 | wl | 175–176 | `expectZero` endpoints 1 and 11/9 | = A7/A8 | yes |
| B8 | wl | 185–188 | `expectZero` `nProdNat - nProdIndep` (eigenvector path) | deliverable 3, extra independent route | yes |
| B9 | wl | 191–196 | `expectZero` Nprod endpoints | = A9/A10 | yes |

No assertion is tautological: in every case one side is a hand-written closed form and the other is computed (an integral, a matrix inverse/LinearSolve, or an `eigenvects`/`Eigensystem` solve), so a wrong closed form would make the residual nonzero and trip `Exit[1]`/`AssertionError`.

## Findings

### F1 — stale_output

**Severity:** low (informational, non-blocking)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage032_source_map_from_mode_integrals_sympy_audit.txt` (mtime 1779777841 = 2026-05-26 00:44:01)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage032_source_map_from_mode_integrals_mathematica_audit.txt` (mtime 1779777860 = 2026-05-26 00:44:20)
- vs sympy script mtime 1780523951 = 2026-06-03 15:59:11
- vs mathematica script mtime 1780523951 = 2026-06-03 15:59:11

**What's wrong:**
Both saved `.txt` outputs predate their scripts by ~8 days. I checked git to see what changed in the scripts after the outputs were captured: the only commit touching either script after the output's last commit (c6ee7d7, 2026-05-26) is `e2a4780` ("numbering reconciliation Phase 1 (deterministic): 273 doc-only stage-label fixes"). `git show e2a4780` on both scripts is a 10-line cosmetic banner relabel only — every changed `+` line is a `banner("STAGE 15.x …")` → `banner("STAGE 32.x …")` rename, with NO equation, value, variable, or logic change (confirmed: filtering the added lines for anything other than the `STAGE n` token returns empty for both files). Consequently the captured outputs still display the stale `STAGE 15.x` section banners (e.g. sympy output line 3 "STAGE 15.1 …") and the trailing "All Stage 15 checks passed." line, while the live scripts now print "STAGE 32.x". The **mathematics** in both transcripts is current and correct: I algebraically verified the two engines' `mhat_-^2` printouts (sympy output line 105; wl output line 65) are identical (both reduce to `11/18 + (968·α_0 + 63·dK·π²)/(18·√(7744·α_0² + 1008·α_0·dK·π² + 81·dK²·π⁴))`), and every `expectZero`/`expect_zero` residual is 0 / PASS.

Secondary (same numbering-drift class, source-side, not in the output): the sympy docstring line 3 ("Moving-throat PDE Stage 15 SymPy audit.") and the final `print("All Stage 15 checks passed.")` (sympy line 217 / wl line 198) are also stale `Stage 15` tokens that Phase 1's matcher did not catch (it targeted `banner(...)` strings, not this docstring wording or the trailing print). These are pure label strings with zero effect on any assertion or emitted value.

**Why this matters:**
Nothing mathematical breaks; this is purely a freshness/label artifact. It matters only because the verifier's standard re-run will regenerate the transcripts and the banner/print labels will flip to `STAGE 32.x` / "All Stage 32 checks passed.", so the committed outputs should be refreshed to match. The residual `Stage 15` docstring/print tokens belong to the orchestrator's deterministic numbering-reconciliation pipeline (label-only), not to the math red-team's assertion scope.

**Required change:**
Informational. The orchestrator's standard exec re-run will refresh both `.txt` outputs (this is the documented `stale_output` workflow). No assertion or value change is required for the math to be correct. If the orchestrator's numbering pipeline picks these up, the source-side label fixes are: sympy line 3 `Stage 15` → `Stage 032` (or `32`); sympy line 217 / wl line 198 `All Stage 15 checks passed.` → `All Stage 032 checks passed.` — label-only, no logic. (Out of scope for a Codex math directive; listed for orchestrator awareness only.)

**Verification:**
After re-run, the refreshed `.txt` outputs should show `STAGE 32.1 …` banners and (if the label tokens are fixed) "All Stage 032 checks passed."; every residual line should remain `0`/`PASS`, and the two engines' `mhat_-^2` printouts should remain algebraically identical.

## Independent-derivation check (Mathematica)

The `.wl` is parallel in structure (same section banners, same check names) but performs genuine independent re-derivation at the load-bearing steps, so it is not a transliteration:

- **Schur (§3):** sympy uses `Bmat * Kint.inv() * Bmat.T` (sympy line 132). Mathematica deliberately does NOT just invert: it solves `kInt . y == Transpose[bMat] . z` for arbitrary external `z` via `LinearSolve`, reads off the Schur matrix from the linear coefficients of `z1,z2` (wl lines 104–117), AND separately does the direct `Inverse` route (wl line 120), then cross-checks the two against each other (wl line 121, `expectZero["Schur via Inverse vs LinearSolve", …]`). This is an extra independent route the sympy script lacks.
- **Source map (§4):** both engines independently construct the loaded operator and call an eigensolver (`M_loaded.eigenvects()` sympy 172; `Eigensystem[mLoaded]` wl 157), select the lower branch by numeric probe, normalize, and verify `(v·e_-)^2` equals the hand-written closed form `s_minus`. Independent solver, same target.
- **Elimination (§5):** Mathematica adds `nProdNat - nProdIndep` (wl 181–188), routing the product through the eigenvector-derived `sCheck`/`lamMinusIndep` rather than the closed forms — a second independent path with no sympy counterpart.

Conclusion: independent derivation confirmed; no `mathematica_transliteration` finding.

## Engine cross-check

Both engines emit the same closed form for the central result. SymPy output line 105:
`(63*pi**2*DeltaK + 968*alpha0 + 11*sqrt(4608*alpha0**2 + (9*pi**2*DeltaK + 56*alpha0)**2))/(18*sqrt(4608*alpha0**2 + (9*pi**2*DeltaK + 56*alpha0)**2))`.
Mathematica output line 65:
`(11 + (968*alpha0 + 63*dK*Pi^2)/Sqrt[7744*alpha0^2 + 1008*alpha0*dK*Pi^2 + 81*dK^2*Pi^4])/18`.
Expanding `(9π²DK + 56α_0)² = 81π⁴DK² + 1008π²DK·α_0 + 3136α_0²`, the sympy radicand `4608α_0² + that = 7744α_0² + 1008π²DK·α_0 + 81π⁴DK²`, identical to the Mathematica radicand; the two surface forms are algebraically equal. All other residuals are 0 in sympy and PASS in Mathematica. The Mathematica `Limit::alimv` warnings on the two `α_0→∞` limits are benign (Mathematica ignores limit-variable assumptions; the limits still evaluate to 11/9 and 0 and PASS — these are finite limits, not pole evaluations, so the MEMORY pole-non-determinism caveat does not apply). Engines agree.

## Verdict justification

`verdict: findings` with a single low-severity, non-blocking `stale_output` finding. Attacks attempted and failed to break the unit: (a) tautology — every `expectZero` pits a hand-written closed form against an independently computed integral / matrix solve / eigensolve, so none can pass vacuously; (b) endpoint correctness — I hand-verified `s_-(α_0=0)=κ_0^2` (⇒ `\widehat m_-^2=1`) and `s_-(α_0→∞)=σ` (⇒ `11/9`), matching the card's bound exactly; (c) engine disagreement — the two `mhat_-^2` surface forms are algebraically identical; (d) symbol domains — the §4/§5 positivity assumptions (`A,DK,σ,Kprod,β_0 > 0`, `α_0 ≥ 0`) are consistent with the physical setup (loading parameter `α_0 ≥ 0`, positive stiffnesses) and with the card; (e) Mathematica independence — confirmed by the LinearSolve and eigenvector-path cross-checks. I read the card, the notes, and the appendix row: the script's verified claims (`\widehat m_-^2 = s_-/\kappa_0^2`, bound `1 ≤ \widehat m_-^2 < 11/9`, eliminated product `N_- = β_0 s_-^2/(\kappa_0^2 λ_-)` with correct endpoints) match the paper deliverables exactly. The only defect is output freshness (and same-class stale `Stage 15` label tokens in the docstring/print), which is cosmetic and resolved by the verifier's standard re-run; no directive is warranted because there is no script-side math change for Codex to apply.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value the scripts emit, located in the docs. The `.tex` card is terse (symbolic only); the explicit closed values legitimately live in the `.md` notes — per the augmentation guards a value living correctly in the notes is a MATCH.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `κ_0 = 2√2/π` | py 62 / wl 55; sympy out 9, wl out 9 | notes:67 `kappa_0 = 2 sqrt(2) / pi` | MATCH |
| `κ_1 = -4/(3π)` | py 63 / wl 56; sympy out 10, wl out 10 | notes:69 `kappa_1 = - 4 / (3 pi)` | MATCH |
| `σ = κ_0²+κ_1² = 88/(9π²)` | py 64 / wl 57; sympy out 11, wl out 11 | notes:75 `sigma = ... = 88 / (9 pi^2)`; tex:42 (`σ=κ_0²+κ_1²` symbolic) | MATCH |
| `κ_0² = 8/π²` | implicit denom in py 65 / wl 58 (`σ/κ_0²`) | notes:81 `kappa_0^2 = 8 / pi^2` | MATCH |
| `σ/κ_0² = 11/9` | py 65 / wl 58; sympy out 12, wl out 12 | tex:46 `=\frac{11}{9}`, tex:70 `=1+2/9`; notes:83 `sigma / kappa_0^2 = 11 / 9` | MATCH |
| `\widehat m_-^2 = s_-/κ_0²` (closed form) | py 163/165 / wl 148/150; sympy out 105, wl out 65 | tex:40 `\widehat m_-^{\,2}=\frac{s_-}{\kappa_0^2}`; notes:205 | MATCH |
| `\widehat m_-^2(α_0=0) = 1` | py 195 / wl 175; sympy out 108, wl out 70 | tex:46 lower bound, tex:69 onset; notes:220/224 | MATCH |
| `\widehat m_-^2(α_0→∞) = 11/9` | py 196 / wl 176; sympy out 109, wl out 74 | tex:46 upper bound; notes:220/224 | MATCH |
| `Ξ = g_U²/A_U` (Schur) | py 135 / wl 124; sympy out 95, wl out 55 | notes:161 `Xi_0 = g_U^2 / Omega_U^2` (kernel form); notes:151 `Xi(omega)=g_U^2/A_U(omega)` | MATCH |
| `α` Schur coefficient (closed form) | py 136 / wl 125; sympy out 96, wl out 56 | notes:152–156 `alpha(omega)=...`, `Delta_UW=A_U A_W - g_R^2 σ` | MATCH |
| `Σ = Ξ I + α vv^T` (decomposition) | py 143 / wl 134; sympy out 97–100, wl out 57 | notes:147 `Sigma_wall = Xi I_2 + alpha v v^T`; tex (Schur implicit) | MATCH |
| local-kernel reductions `L_src=g_Q Q (v·q)` etc. | py 91–110 / wl 81–85; sympy out 25–34, wl out 33–47 | notes:110–116, 182; tex:21–28 | MATCH |
| `N_- = β_0 s_-²/(κ_0² λ_-)` | py 203–204 / wl 179–180 (endpoints out 114–115, wl out 82–87) | tex:56 `=\frac{\beta_0 s_-^2}{\kappa_0^2 λ_-}`; notes:248 | MATCH |
| `N_-(α_0=0) = β_0 κ_0²/A` | py 207–211 / wl 191–192; sympy out 114, wl out 82 | derived endpoint; consistent with tex:56 + onset | MATCH (internal endpoint of deliverable) |

INTERNAL scaffolding (accounted for, no finding): basis orthonormality residuals `u_i·u_j=δ_ij`, `f_0·f_0=1` (verification of the basis, not deliverables); `λ_-` closed form (intermediate, cross-checked vs eigensolve); `R` (radical), `s_minus`/`lamMinus` symbolic intermediates; probe-point sort rule (`α_0=A=DK=1`); `Schur via Inverse vs LinearSolve` and `nProdNat - nProdIndep` engine-internal cross-route residuals; the Mathematica unused symbol declarations `gConst, cs, radius, cSpeed, mhat` (wl lines 138–141) — declared but never used in any computation or assertion, benign dead code, emits no value. The 2.5PN target literal `54 G c_s^5/(5 a^5 c^5)` (notes:252, tex:62 as `N_Q^{target}`) is NOT script-emitted (no `54/G/c_s/a/c` in either script) → out of reconciliation scope for this unit.

reconciliation: complete; 14 deliverable values checked, 0 misaligned.

## Self-test notes

I checked the traps relevant here. (1) Variable independence: no spurious `diff`/`D` of a constant — the only derivatives are SymPy `sp.limit`/Mathematica `Limit` of `α_0`-dependent expressions, and `mhat_sq` genuinely depends on `α_0`, so the limits are non-trivial (verified by hand: 1 at 0, 11/9 at ∞). (2) Symmetry/parity: the §1 integrals are over `[0,L]` (not a symmetric domain) and the printed nonzero overlaps `κ_0,κ_1` confirm they don't vanish; no false "vanishes" claim. (3) Trivial-case pre-check: I substituted `α_0=0` into `s_minus` by hand and got `κ_0²` (⇒ `\widehat m_-^2=1`), and `α_0→∞` gives `σ` (⇒ `11/9`), matching both the assertions and the card. No directive is written (the sole finding is informational `stale_output`, resolved by the verifier's re-run), so the path-specification and paper-round-trip self-tests are not applicable.
