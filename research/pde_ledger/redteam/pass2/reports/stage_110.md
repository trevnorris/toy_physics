---
unit_id: 110
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage110_robin_outlet_model.md]
  paper_appendix: present
---

# Audit unit 110 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_110.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage110_robin_outlet_model.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: §"Outlet DtN and Robin outlet tests" L345, "Robin and mixed-pole tests" L399–412, eq:app-part04-robin-outlet L405, eq:app-part04-chi-robin L410; canonical outgoing fingerprint L305–322)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.txt`

## What the paper claims

Stage 110 takes the canonical outgoing `l=2` DtN expansion `Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9 + O(z^6)` (carried forward; identical to the appendix eq. L308–313 and to stage 108's `Lambda_out`) and adds an isotropic Robin throat-core shift `Lambda_2^R(z) = Lambda_2^out(z) + rho_R`. Normalizing by the actual static slot `(-3 + rho_R)`, the notes box the exact low-frequency response `Y_2^R = 1 + z^2/(9-3 rho_R) + (4-rho_R)/(9(3-rho_R)^2) z^4 + i z^5/(27-9 rho_R) + O(z^6)`. Relative to the canonical outgoing fingerprint `Y_2^out = 1 + z^2/9 + 4 z^4/81 + i z^5/27` (appendix L317–321), the deliverable normalization factor is `\chi_Q^R = 3/(3-rho_R)` (card quote L16: "Pure Robin shift gives \(\chi_Q^{\rm R}=3/(3-\rho_R)\) and generically spoils the branch"; appendix eq:app-part04-chi-robin L410), with linearization `chi_Q^R = 1 + rho_R/3 + rho_R^2/9 + O(rho_R^3)` (notes L60). The card's `\stagefield{Purpose}`/`\stagefield{Verification}` make the verification transcript the audit target; the `Checks` checklist asks for pure-scale/argument deformations, Robin/mixed-pole limits, and even+odd preservation. Distinct deliverables in this stage's scripts: the four series coefficients (c2, c4, c5, and the implicit constant 1), `chi_Q^R`, and its linearization. (The Stage-92 branch-triple `(b,a_0,a_5)=(0,rho_R,0)` and `delta chi_Q^R = rho_R/3` are notes-side mappings, not script outputs.)

## What the script claims to verify

Both engines construct `Lambda_R = Lambda_out + rho`, form `Y_R = (-3+rho)/Lambda_R`, series-expand to z^5, extract the z^2/z^4/z^5 coefficients (the z^5 coefficient divided by `I`), define `chi_R = c5/(1/27)` (ratio of the Robin odd coefficient to the canonical `Y_2^out` odd coefficient `1/27`), and linearize `chi_R` in `rho` to order rho^2. Five assertions then confirm `c2 = 1/(9-3 rho)`, `c4 = (4-rho)/(9(3-rho)^2)`, `c5 = 1/(27-9 rho)`, `chi_R = 3/(3-rho)`, and `chi_R_lin = 1 + rho/3 + rho^2/9`. The RHS literals are the paper's boxed results; the LHS are re-derived from the Robin DtN by symbolic series, so the checks are non-tautological. The SymPy script declares `z, rho` real; the `.wl` adds `rho != 3` (the pole) to its assumptions.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Y_2^R` constant term = 1 (canonical normalization) | implicit: series starts at 1 (printed `Y_R(z) = 1 + ...`) | match (printed, not asserted, but follows from c2..c5 + structure) |
| `c2 = 1/(9-3 rho_R)` (even fingerprint shift) | py:26 / wl:49 `c2 - 1/(9-3 rho) == 0` | match |
| `c4 = (4-rho_R)/(9(3-rho_R)^2)` | py:27 / wl:50 `c4 - (4-rho)/(9(3-rho)^2) == 0` | match |
| `c5 = 1/(27-9 rho_R)` (odd normalization shift) | py:28 / wl:51 `c5 - 1/(27-9 rho) == 0` | match |
| `chi_Q^R = 3/(3-rho_R)` (Output / boxed) | py:29 / wl:52 `chi_R - 3/(3-rho) == 0` | match |
| `chi_Q^R linearized = 1 + rho_R/3 + rho_R^2/9` | py:30 / wl:53 | match |
| Stage-92 triple `(0, rho_R, 0)`; `delta chi_Q^R = rho_R/3` | not emitted (notes-only mapping) | not a script deliverable — no finding |

All script-side checks trace to a specific paper deliverable; no orphan assertions. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 26 | `simplify(c2 - 1/(9-3 rho)) == 0` | even-branch coeff c2 | yes |
| A2 | sympy | 27 | `simplify(c4 - (4-rho)/(9(3-rho)^2)) == 0` | even-branch coeff c4 | yes |
| A3 | sympy | 28 | `simplify(c5 - 1/(27-9 rho)) == 0` | odd-normalization coeff c5 | yes |
| A4 | sympy | 29 | `simplify(chi_R - 3/(3-rho)) == 0` | `chi_Q^R` (Output) | yes |
| A5 | sympy | 30 | `simplify(chi_R_lin - (1+rho/3+rho^2/9)) == 0` | `chi_Q^R` linearization | yes |
| A6 | mathematica | 49 | `expectZero[c2 - 1/(9-3 rho)]` | c2 | yes |
| A7 | mathematica | 50 | `expectZero[c4 - (4-rho)/(9(3-rho)^2)]` | c4 | yes |
| A8 | mathematica | 51 | `expectZero[c5 - 1/(27-9 rho)]` | c5 | yes |
| A9 | mathematica | 52 | `expectZero[chiR - 3/(3-rho)]` | `chi_Q^R` | yes |
| A10 | mathematica | 53 | `expectZero[chiRLinear - (1+rho/3+rho^2/9)]` | linearization | yes |

All ten rows are anchored: the LHS coefficients are produced by symbolic series of the Robin DtN, and the RHS are the paper's claimed values — the assertion fails if the series produces anything else. No tautology, no hardcoded LHS.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl:31-53`
- (parallel to) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.py:8-30`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:

1. Identical premise construction:
   - py:8–9 `Lambda_out = -3 + z**2/3 + z**4/9 + I*z**5/9` ; `Lambda_R = Lambda_out + rho`
   - wl:31–32 `lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9` ; `lambdaR = lambdaOut + rho`

2. Identical normalization + series + coefficient choreography:
   - py:10–16 `Y_R = simplify((-3+rho)/Lambda_R)`; `series(Y_R,z,0,6).removeO()`; `c2=coeff(z,2)`; `c4=coeff(z,4)`; `c5=coeff(z,5)/I`; `chi_R = c5/Rational(1,27)`
   - wl:33–39 `yR = FullSimplify[(-3+rho)/lambdaR]`; `Series[yR,{z,0,5}]//Normal//Expand`; `Coefficient[...,z,2]`; `...z,4`; `Coefficient[...,z,5]/I`; `chiR = c5/(1/27)`

3. Identical linearization and identical assertion order/content:
   - py:17 + py:26–30 `series(chi_R,rho,0,3)`; then the five `assert simplify(LHS - RHS) == 0`
   - wl:40 + wl:49–53 `Series[chiR,{rho,0,2}]`; then the five `expectZero[LHS - RHS]` in the SAME order with the SAME RHS literals.

Both engines `Series`-expand the *same* expression `(-3+rho)/(Lambda_out+rho)`, extract coefficients via the same `/I` and `/(1/27)` operations, and assert the same five identities in the same sequence. The only differences are syntactic (`series` vs `Series`, order 6 vs 5, `simplify` vs `FullSimplify`, `Together[Expand[]]` inside `expectZero`). There is no second, structurally distinct route to `chi_Q^R` (e.g., a direct closed-form `chi_Q^R = -3 L_5_canonical / L_5_Robin` style derivation, a residue/limit computation, or an algebraic factorization of `(-3+rho)/Lambda_R` that does not pass through the same `Series`). This is the 105–175 first-pass orchestrator-direct band flagged for transliteration watch, and the `.wl` re-types the `.py`.

**Why this matters:**
The second-engine policy requires Mathematica to independently corroborate the result from the physical premises, so that a single algebra error or SymPy-specific series convention is caught. A transliteration cannot catch such an error — both engines share the identical derivation path, so they would share the identical mistake. The result is correct (verified by hand below), so this is an independence defect, not a correctness defect; severity medium.

**Required change:**
Re-author the `.wl` so it reaches `chi_Q^R = 3/(3-rho_R)` (and the c2/c4/c5 coefficients) by a route that does NOT mirror the SymPy `Series`-then-`Coefficient` pipeline. One viable independent route: compute `chi_Q^R` directly as the ratio of the leading odd outgoing coefficients without series-expanding `Y_R` — e.g. form the canonical `Y_2^out = -3/Lambda_out` and the Robin `Y_2^R = (-3+rho)/Lambda_R`, take `chiR := Limit[Coefficient-free ratio]` via `Series[yR/yOut, {z,0,5}]` leading odd term, or derive `chi_Q^R = (-3+rho)/(-3) * (Lambda_out z^5-coeff)/(Lambda_R z^5-coeff)` algebraically from the closed forms; and obtain the c2/c4 via `SeriesCoefficient[yR,{z,0,k}]` (a different extraction primitive than the `.py`'s `coeff`) or via direct `Residue`/`Normal[Series[...]]` of `Y_2^R - Y_2^out`. The acceptance criterion is that the `.wl` does not perform the identical `Series[(-3+rho)/lambdaR,{z,0,5}]` → `Coefficient[...,z,k]` → `/(1/27)` → `Series[...,{rho,0,2}]` sequence in the same order with the same RHS literals. (Designing the concrete independent route is Codex's job per the contract; the directive states the requirement and the no-mirror acceptance test, not a pre-designed script.)

**Verification:**
The refreshed `.wl` must (a) still emit `chi_Q^R = -3/(rho-3)` and `chi_Q^R linearized = 1 + rho/3 + rho^2/9` and the same c2/c4/c5 in the output `.txt`, (b) exit 0, and (c) no longer contain the line-for-line mirror of py:8–30 — specifically the verifier should confirm the coefficient-extraction primitive and/or the route to `chiR` differs structurally from `Series[...]`→`Coefficient[...]/(1/27)`.

## Independent-derivation check (Mathematica)

Not independent. See F1. The `.wl` mirrors the `.py` premise-by-premise, coefficient-by-coefficient, and assertion-by-assertion. Note that sibling stage 108's engines use the same `(Coefficient[yArg,z,5]/I)/(1/27)` band idiom, so this is the established band pattern — but the transliteration-watch directive for 105–175 explicitly asks that a line-by-line port be flagged even when correct, which this is.

## Engine cross-check

Both engines agree exactly. SymPy output: `c2 = -1/(3*rho - 9)`, `c4 = (4 - rho)/(9*(rho**2 - 6*rho + 9))`, `c5 = -1/(9*rho - 27)`, `chi_Q^R = -3/(rho - 3)`, `chi_Q^R linearized = rho**2/9 + rho/3 + 1`, `stage110: PASS`. Mathematica output: `c2 = (9 - 3*rho)^(-1)`, `c4 = -1/9*(-4 + rho)/(-3 + rho)^2`, `c5 = (27 - 9*rho)^(-1)`, `chi_Q^R = -3/(-3 + rho)`, `chi_Q^R linearized = 1 + rho/3 + rho^2/9`, all five `PASS`. These are identical expressions in different normal forms (`-3/(rho-3) == 3/(3-rho)`, `(4-rho)/(9(rho^2-6rho+9)) == (4-rho)/(9(3-rho)^2)`). No `engine_disagreement`.

## Verdict justification

`verdict: findings` with one finding: `mathematica_transliteration` (medium). The mathematics is correct — I hand-verified the geometric-series expansion of `Y_R = 1/(1+u)` with `u = [z^2/3 + z^4/9 + i z^5/9]/(rho-3)` and reproduced c2, c4, c5, `chi_Q^R = 3/(3-rho)`, and the `1 + rho/3 + rho^2/9` linearization; the carried-forward `Lambda_2^out` matches the appendix (L308–313) and stage 108 (`Lambda_out`) exactly; the assertions are non-tautological (RHS are the paper's boxed values, LHS re-derived by series), the `rho != 3` pole domain is correctly handled, and all six emitted deliverables reconcile with the card/notes/appendix (0 misaligned). Paper alignment is exact. The only defect is that the `.wl` is a line-by-line port of the `.py` rather than an independent second-engine derivation — flagged per the 105–175 transliteration watch. This is not `paper_misalignment` (no user gate), not stop-cold (correct result, no downstream propagation), and `engines_agree: true`. The independence direction (re-author vs. accept) is a user-level call per the dual-engine policy.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation (every RESULT value the scripts emit):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Y_2^R(z) = 1 + z^2/(9-3 rho) + (4-rho)/(9(3-rho)^2) z^4 + i z^5/(27-9 rho)` | py out L1; wl out L5 | notes L31–42 (boxed) | MATCH |
| `c2 = 1/(9-3 rho_R)` | py out L2; wl out L6 | notes L36 (`z^2/(9-3 rho_R)`) | MATCH |
| `c4 = (4-rho_R)/(9(3-rho_R)^2)` | py out L3; wl out L7 | notes L37–38 | MATCH |
| `c5 = 1/(27-9 rho_R)` | py out L4; wl out L8 | notes L39 (`i z^5/(27-9 rho_R)`) | MATCH |
| `chi_Q^R = 3/(3-rho_R)` | py out L5; wl out L9 | card L16 (`\chi_Q^{\rm R}=3/(3-\rho_R)`); notes L55; appendix L410 (eq:app-part04-chi-robin) | MATCH |
| `chi_Q^R linearized = 1 + rho_R/3 + rho_R^2/9` | py out L6; wl out L10 | notes L60 | MATCH |

Internal scaffolding (accounted for, no finding): pass/fail flags (`stage110: PASS`, `PASS:`/`FAIL:`), the `expectZero` residuals (all `= 0`), the banner string, `fmt`/`pass`/`fail`/`banner` helper definitions, the `$Assumptions`/`real=True` declarations.

Note on the Stage-92 mapping `(b,a_0,a_5)=(0,rho_R,0)` and `delta chi_Q^R = rho_R/3` (notes L68, L72): these are notes-side derivations to the stage-92 linearized notation; the scripts do not emit them, so they are not script deliverables and require no reconciliation.

Output freshness: both `.txt` outputs are newer than their scripts (py output 15:18 > py script 15:08; wl output 15:24 > wl script 15:08) — fresh, no `stale_output`. The reconciliation is based on script source + committed outputs, no execution performed.

reconciliation: complete; 6 values checked, 0 misaligned.

## Self-test notes

Checked: (1) Variable independence — no `diff`/`D` in this stage; the only differentiation-like operations are `Series` in `z` and `rho`, and both `Y_R` and `chi_R` genuinely depend on those variables, so no identically-zero-derivative trap. (2) Pole/domain — `rho = 3` is the genuine pole of `1/(3-rho)`; SymPy declares `rho` real (unrestricted, safe), `.wl` adds `rho != 3` (correct). (3) Trivial-case pre-check — at `rho=0`, all coefficients reduce to the canonical `1/9, 4/81, 1/27` and `chi_Q^R → 1` (confirms the Robin shift correctly degenerates to the canonical branch); at small `rho`, `chi_Q^R ≈ 1 + rho/3` matches the linearization. (4) The transliteration finding's proposed remedy does not introduce a new constant or change any RHS literal, so no new `paper_misalignment` is created (paper round-trip clean). I tried to break the assertions as tautological (LHS is series-derived, not hardcoded — fails to break) and as paper-misaligned (all six deliverables reconcile — fails to break); the only surviving attack is the independence/transliteration one.
