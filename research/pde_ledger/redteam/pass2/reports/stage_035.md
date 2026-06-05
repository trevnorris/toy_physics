---
unit_id: 035
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
  notes_stage_files: [moving_throat_pde_stage035_dimensionless_normalization_locus.md]
  paper_appendix: present
---

# Audit unit 035 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_035.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage035_dimensionless_normalization_locus.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 60 references stage 035)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.txt`

## What the paper claims

Stage 035 non-dimensionalizes the Stage-034 softening-depth target. Inserting the exact D/N finite-throat constants `kappa_0^2 = 8/pi^2`, `kappa_1^2 = 16/(9 pi^2)` (ratio `eta = 2/9`) and dimensionless variables `xi = x/A`, `delta = DeltaK_ax/A` (with `0 <= xi < 1`, `delta > 0`), the normalization product factors as `N_-(x) = N_-(0) F(xi,delta)` with `N_-(0) = beta_0 kappa_0^2/A`. The `\stagefield{Output}` is verbatim: *"Stage~035 outputs the shape function \eqref{eq:app-stage035-F}, the target ratio \eqref{eq:app-stage035-Rtarget}, the monotonicity derivative \eqref{eq:app-stage035-F-derivative}, the existence/uniqueness theorem \eqref{eq:app-stage035-existence}, and the required total loading \eqref{eq:app-stage035-alpha-req}."* Distinct deliverables: (1) closed-form `F(xi,delta) = (9 delta + 11 xi)^4 / [81 (1 - xi)(9 delta^2 + 18 delta xi + 11 xi^2)^2]`; (2) `R_target = N_Q^target A/(beta_0 kappa_0^2)` and locus `F = R_target`; (3) the manifestly-positive `dF/dxi`; (4) endpoints `F(0,delta) = 1`, `lim_{xi->1^-} F = +infinity` with leading constant `C(delta) = (9 delta + 11)^4/[81(9 delta^2 + 18 delta + 11)^2]` (this drives the existence/uniqueness trichotomy); (5) `alpha_req(xi,delta) = 9 pi^2 A xi (xi + delta)/[8(9 delta + 11 xi)]` and its `xi->1^-` limit `alpha_crit = 9 pi^2 A (1 + delta)/[8(11 + 9 delta)]`. The notes additionally give the near-onset series of `F` through `O(xi^2)`, the support-coupling split `g_B,req^2/varpi^2 = alpha_req - alpha_mix` with `alpha_mix = Chi^2/(Omega_U^2 Delta_0)`, and the two asymptotic locus formulas. The appendix row (line 60) summarizes: "Shape function F(xi,delta), monotonicity, endpoint theorem, and total loading."

## What the script claims to verify

Both scripts construct the dimensional ingredients from the physical premises — `x = A xi`, `DeltaK = A delta`, the two kappa constants, the directional-loading ratio `alpha_x`, and the normalization product `N_x` (a degree-4-over-degree-5 rational in the kappa-weighted blocks) — and then reduce. Section 35.1 forms `F = N_x/(beta_0 kappa_0^2/A)` from the physical construction and asserts `F - F_target == 0` against the independently typed closed form; it also pins `R_target` (Mathematica adds an independent `R_target == Pi^2 A NQ/(8 beta_0)` check). Section 35.2 asserts `dF/dxi` equals the manifestly-positive factored form, `F(0,delta) = 1`, the softening constant `C(delta)` (via a `(1-xi) F` limit at `xi->1^-` AND an independent `eps_soft` series), all against typed closed forms. Section 35.3 asserts `alpha_req` and its `alpha_crit` limit against typed forms and prints the support-coupling combination `alpha_req - alpha_mix`. Section 35.4 asserts the near-onset series of `F` (through `O(xi^2)`) and of `alpha_req` against typed forms. Each `expect_zero`/`expectZero` subtracts a from-construction expression from an independently typed target, so the checks are non-tautological.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| (1) `F` closed form, eq app-stage035-F | sympy L53–58 `F - F_target == 0`; wl L50–61 `f - fTarget == 0` | match |
| (2) `R_target`, eq app-stage035-Rtarget; locus `F = R_target` | sympy L60–62 prints `R_target`; wl L55–62 asserts `R_target == Pi^2 A NQ/(8 beta_0)` | match |
| (3) `dF/dxi` positive form, eq app-stage035-F-derivative | sympy L65–73; wl L64–71 `dF - dFTarget == 0` | match |
| (4a) `F(0,delta) = 1` | sympy L74; wl L72 | match |
| (4b) `lim F = +inf` via `C(delta)`, eq app-stage035-Cdelta | sympy L75–82 (limit + eps_soft series); wl L74–90 | match |
| (4c) existence/uniqueness trichotomy, eq app-stage035-existence | implied by monotonicity + endpoints + `C(delta)` (no separate numeric assertion) | match (theorem follows from verified ingredients) |
| (5a) `alpha_req`, eq app-stage035-alpha-req | sympy L85–89; wl L92–105 `alpha_req - target == 0` | match |
| (5b) `alpha_crit`, eq app-stage035-alpha-crit | sympy L91–94; wl L97–106 | match |
| notes: support coupling `g_B,req^2/varpi^2 = alpha_req - alpha_mix` | sympy L96–99; wl L108–110 (printed, not asserted — definitional) | match |
| notes §6.1 near-onset `F` series O(xi^2) | sympy L102–105; wl L112–122 | match |
| notes §6.2 `1-xi_req ~ C(delta)/R_target` | covered by verified `C(delta)`; algebraic rearrangement | match |

All paper-side deliverables map to a faithful script-side check; no `mismatch`, `missing`, or `extra` rows. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `simplify(F - F_target) == 0` | (1) F closed form | yes |
| A2 | sympy | 73 | `simplify(dF - dF_target) == 0` | (3) monotonicity derivative | yes |
| A3 | sympy | 74 | `F_target(xi=0) - 1 == 0` | (4a) onset endpoint | yes |
| A4 | sympy | 78 | `soft_const - C(delta) == 0` (limit) | (4b) softening constant | yes |
| A5 | sympy | 82 | `eps_soft series coeff - C(delta) == 0` | (4b) softening constant (indep route) | yes |
| A6 | sympy | 89 | `alpha_req - target == 0` | (5a) total loading | yes |
| A7 | sympy | 94 | `alpha_crit - target == 0` (limit) | (5b) crit threshold | yes |
| A8 | sympy | 105 | near-onset F series == target | notes §6.1 | yes |
| A9 | sympy | 110 | near-onset alpha series == target | notes §6.1 (leading term) | yes |
| B1 | mathematica | 61 | `f - fTarget == 0` | (1) F closed form | yes |
| B2 | mathematica | 62 | `rTarget - Pi^2 A NQ/(8 beta0) == 0` | (2) R_target | yes |
| B3 | mathematica | 71 | `dF - dFTarget == 0` | (3) monotonicity | yes |
| B4 | mathematica | 72 | `fTarget(xi=0) - 1 == 0` | (4a) onset | yes |
| B5 | mathematica | 83 | `softConst - target == 0` (limit) | (4b) softening | yes |
| B6 | mathematica | 90 | `softScaledSeries - target == 0` | (4b) softening (indep route) | yes |
| B7 | mathematica | 105 | `alphaReq - target == 0` | (5a) total loading | yes |
| B8 | mathematica | 106 | `alphaCrit - target == 0` (limit) | (5b) crit | yes |
| B9 | mathematica | 122 | near-onset F series == target | notes §6.1 | yes |
| B10 | mathematica | 123 | near-onset alpha series == target | notes §6.1 | yes |

Every assertion subtracts a from-construction expression (`F`/`f`, `dF`, `alpha_req`, limits, series) from an independently typed target; none is `x == x` by construction. The shape function `F` is reduced from `N_x` (built from the kappa constants) and compared against `F_target` (typed by hand) — a genuine algebraic identity that fails if the construction were wrong. The two limit checks (A4/A5 and B5/B6) cross-validate `C(delta)` by two independent routes (one-sided limit vs. epsilon series), strengthening section 4.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.txt` (mtime 2026-05-21 17:35:55)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.txt` (mtime 2026-05-21 17:36:00)
- scripts: `..._sympy_audit.py` and `..._mathematica_audit.wl` (both mtime 2026-06-03 15:59:11)

**What's wrong:**
Both saved `.txt` outputs are older than their scripts. The single commit that touched these scripts after the outputs were captured is `e2a4780` ("numbering reconciliation Phase 1 (deterministic): 273 doc-only stage-label fixes"), which is **label-only**: it changed banner strings (`STAGE 018` → `STAGE 035`, `STAGE 18.1` → `STAGE 35.1`, etc.) and nothing in the equations, values, variables, or assertions. The visible consequence is that the saved outputs still print the OLD banners: the sympy `.txt` shows `STAGE 18.1 — EXACT D/N DIMENSIONLESS SHAPE FUNCTION` (line 3) and `All Stage 18 checks passed.` (line 66); the Mathematica `.txt` shows `STAGE 018 — DIMENSIONLESS NORMALIZATION LOCUS` (line 3). The computed results in both transcripts are still correct (every `... = 0` residual and every `PASS:` line) because the algebra was never touched.

**Why this matters:**
The content of the captured outputs is trustworthy; only the banner labels are stale. This is informational — the orchestrator's re-run will refresh the transcripts (and update the banners to "STAGE 035"). Not blocking. (Note: the sympy script's module docstring at `.py:3` still says "Stage 18 SymPy audit" and its final print at `.py:112` says "All Stage 18 checks passed." — cosmetic label residuals outside the value-correctness scope; flag here for completeness, no math impact.)

**Required change:**
None to the script math. The verifier should re-run both scripts to regenerate fresh transcripts (which will carry the corrected "STAGE 035" banners). Optionally, Codex may relabel the residual "Stage 18" string at `.py:3` (docstring) and `.py:112` (final print) to "Stage 035" for consistency, but this is cosmetic and not load-bearing.

**Verification:**
After re-run, the sympy `.txt` header line should read `STAGE 35.1 — ...` and the Mathematica `.txt` line 3 should read `STAGE 035 — DIMENSIONLESS NORMALIZATION LOCUS`; all residual lines remain `0` and all `PASS:` lines remain.

## Independent-derivation check (Mathematica)

The `.wl` is structurally parallel to the `.py`: both define `x = A xi`, `deltaKSub = A delta`, the same `alphaX`, the same `nX`, then `f = nX/(beta0 kappa0Sq/A)`, the same typed targets, the same limits and series. This is the same variable choreography. However, the parallelism is at the level of the **shared physical premises** (the two kappa constants and the `N_x`/`alpha_x` construction), which is the legitimate common starting point — not one engine echoing the other's intermediate algebraic manipulations. Each engine then reduces `nX` to closed form independently using its own simplifier (`sp.simplify`/`sp.series`/`sp.limit` vs. `FullSimplify`/`Series`/`Limit[..., Direction->"FromBelow"]`). The Mathematica script is also not a pure mirror: it adds an assertion the sympy script lacks — `expectZero["R_target - Pi^2 A NQ/(8 beta0)", ...]` (wl:62) verifies the `R_target` closed form, whereas sympy only prints `R_target` (py:60–61). Verdict: acceptable second-engine independence; not a `mathematica_transliteration` finding. Quoted corresponding sections: sympy `N_x` (py:47–50) vs. Mathematica `nX` (wl:44–48) — same physical construction; sympy `F = sp.simplify(N_x/(beta0*kappa0_sq/A))` (py:53) vs. `f = FullSimplify[nX/(beta0*kappa0Sq/A), ...]` (wl:50) — same reduction performed by independent engines.

## Engine cross-check

Both engines agree on every deliverable:
- `eta = 2/9` (sympy out L5; wl out L5).
- `F`: sympy out L6–12 `-(9δ+11ξ)^4 / [81(ξ-1)(2ξ²+9(δ+ξ)²)]` ≡ wl out L6 `-1/81*(9*delta+11*xi)^4/((-1+xi)*(9*delta^2+18*delta*xi+11*xi^2)^2)` (note `2ξ²+9(δ+ξ)² = 9δ²+18δξ+11ξ²`; both carry the `(ξ-1)` sign that cancels the leading minus to give the positive paper form). `F - closed D/N form = 0` in both.
- `R_target = Pi^2 A NQ/(8 beta0)` (sympy out L14; wl out L7).
- `dF/dxi - manifestly positive form = 0`, `F(0,delta) - 1 = 0` in both.
- `C(delta) = (9δ+11)^4/[81(9δ²+18δ+11)^2]` (sympy out L36; wl out L17 `(11+9*delta)^4/(81*(11+9*delta*(2+delta))^2)` — same expression) — `= 0` residual in both, by both the limit and the eps_soft routes.
- `alpha_req = 9π²Aξ(δ+ξ)/(8(9δ+11ξ))` (sympy out L44–48; wl out L22 `(9*A*Pi^2*xi*(delta+xi))/(72*delta+88*xi)` — same, `8·9=72`, `8·11=88`).
- `alpha_crit = 9π²A(1+δ)/(8(11+9δ))` (sympy out L50; wl out L23 `(9*A*(1+delta)*Pi^2)/(88+72*delta)` — same).
- `g_B,req^2/varpi^2 = alpha_req - Chi²/(Omega_U² Delta_0)` (sympy out L52–57; wl out L28).
- near-onset F series and alpha series residuals `= 0` in both (sympy out L62–65; wl out L29–34).

Engines agree. `engines_agree: true`.

## Verdict justification

Verdict is `findings` solely because of the low-severity `stale_output` (F1): the saved `.txt` transcripts predate a label-only relabel commit and still print "STAGE 18" banners, while the computed content remains correct. There is no script-side math defect and no paper misalignment. Attacks tried and failed: (a) traced the full `N_x → F` reduction by hand and confirmed it cancels all `pi` and `A` factors to exactly `(9δ+11ξ)^4/[81(1-ξ)(9δ²+18δξ+11ξ²)²]` — the assertion is a genuine identity, not tautological; (b) verified `alpha_x → 9π²Aξ(ξ+δ)/(8(9δ+11ξ))` by hand; (c) checked symbol domains (`xi` nonnegative, `delta` positive, `0<=xi<1` in Mathematica) match the paper's stated branch `0<=xi<1, delta>0`; (d) confirmed the softening constant is cross-validated by two independent routes; (e) confirmed the second engine is not a pure transliteration and adds an independent `R_target` assertion. I read the paper card, the notes, and the appendix row; the script's verified claims match the paper's `\stagefield{Output}` deliverables exactly.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 9 deliverable values checked, 0 misaligned.

Every result/deliverable value the scripts emit reconciles with the paper card and/or notes. Output `.txt` files are stale-by-label only (see F1); the reconciliation is based on the script source plus the saved outputs, whose computed content is unchanged from the current scripts.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `eta = 2/9` | py:35 / wl:36; sympy out L5, wl out L5 | tex:24 `\eta=...=\frac29`; md:33 `eta := ... = 2/9` | MATCH |
| `F(xi,delta) = (9δ+11ξ)^4/[81(1-ξ)(9δ²+18δξ+11ξ²)²]` | py:54 / wl:51–54; sympy out L6–12, wl out L6 | tex:46–48 (eq app-stage035-F); md:61–63 | MATCH |
| `R_target = Pi^2 A NQ/(8 beta0)` ( = NQ A/(beta0 kappa0²)) | py:60 / wl:55–56; sympy out L14, wl out L7 | tex:53–55 (eq app-stage035-Rtarget); md:71–72 | MATCH |
| `dF/dxi = (9δ+11ξ)^3(81δ³+72δ²+189δ²ξ+297δξ²+121ξ³)/[81(1-ξ)²(9δ²+18δξ+11ξ²)³]` | py:66–69 / wl:65–67; sympy out L20–33, wl out L12 | tex:69–72 (eq app-stage035-F-derivative); md:84–88 | MATCH |
| `F(0,delta) = 1` | py:74 / wl:72; sympy out L35, wl out L15 | tex:77; md:94 | MATCH |
| `C(delta) = (9δ+11)^4/[81(9δ²+18δ+11)²]` | py:76 / wl:78–79; sympy out L36/L38, wl out L17 | tex:86 (eq app-stage035-Cdelta); md:104, md:191 | MATCH |
| `alpha_req = 9π²Aξ(ξ+δ)/[8(9δ+11ξ)]` | py:86 / wl:93–94; sympy out L44–48, wl out L22 | tex:104–105 (eq app-stage035-alpha-req); md:120–121 | MATCH |
| `alpha_crit = 9π²A(1+δ)/[8(11+9δ)]` | py:92 / wl:101; sympy out L50, wl out L23 | tex:110 (eq app-stage035-alpha-crit); md:125 | MATCH |
| near-onset `F` series `1 + (1+8/(9δ))ξ + (1+8/(9δ)-28/(27δ²))ξ²` | py:103 / wl:113–114; sympy out L62, wl out L29 | tex:117–118 (eq app-stage035-onset-expansion); md:167–169 | MATCH |

INTERNAL items (genuine scaffolding / intermediate, no finding):
- `g_B,req^2/varpi^2 = alpha_req - Chi²/(Omega_U² Delta_0)` — printed, not asserted; definitional combination matching notes §5 (md:149–151); reported as a derived gate, value reconciles (informational MATCH-equivalent, not a numeric deliverable).
- `alpha_mix = Chi²/(Omega_U² Delta_0)` — definitional input (notes md:141), not a stage deliverable value.
- `alpha_req` near-onset series O(ξ²) coefficient `-π²A/(36δ)` — script-internal asymptotic; notes §6.1 (md:181) state only the leading term `π²A ξ/8`, which matches; the O(ξ²) term is intermediate detail, not a stated deliverable.
- `x = A xi`, `DeltaK_sub = A delta`, `N_x`, `alpha_x`, `kappa0_sq = 8/pi²`, `kappa1_sq = 16/(9pi²)` — construction inputs (kappa constants are stated in tex:21–24/md:27–29, MATCH; the rest are intermediate).

## Self-test notes

Checked the traps relevant to this unit. (1) Variable independence: the only derivatives/series/limits are in `xi` of expressions (`F_target`, `alpha_req_target`) that genuinely depend on `xi`, so `dF/dxi` and the series are not identically zero — confirmed by the nonzero printed forms in the outputs. (2) Symmetry/parity: no unbounded-domain integrals in this unit, n/a. (3) Trivial-case pre-check: hand-substituting `xi=0` gives `F=1` (matches A3) and the leading near-onset coefficient `(1+8/(9δ))` is nonzero; the `N_x → F` reduction cancels to a `pi`-free, `A`-free rational exactly equal to `F_target`, so `F - F_target` is a real (non-trivial) zero. (5) Paper round-trip: F1 prescribes no math change (re-run only), so it introduces no new misalignment; all nine deliverables reconcile to the card/notes verbatim.
