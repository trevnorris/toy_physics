---
unit_id: 227
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T19:13:50Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 227 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_227.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 66, 747–804 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card (`stage_227.tex:15`) gives `\stagefield{Output}`: "Pure-transfer theorem: $\Xi_1=2\delta\ln(P/\Delta)$. If the outgoing transfer shape is rigid while interference and hybridization ratios are also rigid, then nontrivial co-loading is impossible." The card's `Inputs` (line 9) imports the one-port numerator `P=Ω_U²G_W+RG_U`, denominator `Δ=Ω_U²Ω_W²-R²`, and load factor `Λ=P/Δ`. The notes (the authoritative detail source) enumerate six deliverables: (1) the pure-transfer theorem `Ξ_1 = N_01/N_0 = 2 P_01/P - 2 Δ_01/Δ`; (2) the one-port factorization `Λ = (G_W/Ω_W²)(1+I)/(1-H)` with `I = R G_U/(Ω_U² G_W)`, `H = R²/(Ω_U²Ω_W²)`; (3) `Ξ_1 = 2 δln Λ`; (4) the microscopic slope law `Ξ_1 = 2[m + I/(1+I) i + H/(1-H) h]` with sample specialization `2m + (6/19)i + 50/(98π²-25) h`, plus `I=3/16`, `H=25/(98π²)`; (5) the combined `i=h=0` rigidity determinant being nonzero (the co-loading no-go); and (6) the one-dimensional `i=0/h=0/m=0` unit survivors, their gain scales `σ_i,σ_h,σ_m`, the reference `σ_transfer`, and the transported 10%-loss ceilings. The appendix (lines 750–804) restates the same factorization (eq. app-part07-lambda-factorization), the `Ξ_1=2δln Λ` identity, the co-loading no-go hypothesis, and the rigidity-sieve theorem.

## What the script claims to verify

The SymPy script builds the one-port compiler on the explicit finite-throat sample branch, computes the Stage-226 pure-transfer corridor basis `T` (2D nullspace of the `D01/D21/D41` mixed matrix, line 149), and then asserts: the quotient-rule pure-transfer theorem (line 171), the `Λ` factorization identity (line 176), the microscopic `m,i,h` slope law and its sample specialization (lines 187, 190) plus `I=3/16`, `H=25/(98π²)` (lines 192–193), the combined `i=h` reduced determinant equals a hardcoded closed form and is nonzero (lines 216–217), the three one-dimensional unit survivors against pre-stored numeric vectors and their `|Ξ_1|` gain scales (lines 237–252), and the corridor norm `σ_transfer` plus the four transported 10%-loss ceiling pairs (lines 269, 288–291). All assertions are substantive and the derived intermediate quantities (`N01_mixed`, `P01_mixed`, `Delta01_mixed`, the nullspaces, the reduced rows) are independently differentiated/computed rather than pre-constructed, so the checks are non-tautological.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| (1) Ξ_1 = N_01/N_0 = 2 P_01/P − 2 Δ_01/Δ | line 169–171 | match |
| (2) Λ = (G_W/Ω_W²)(1+I)/(1−H), I, H defs | lines 173–176 | match |
| (3) Ξ_1 = 2 δln Λ | covered transitively via (1)+(2)+(4); Λ slope = m+… checked at 187 | match |
| (4) Ξ_1 = 2[m + I/(1+I)i + H/(1−H)h] + sample law + I=3/16, H=25/(98π²) | lines 179–193 | match |
| (5a) combined i=h determinant nonzero (no-go) | lines 214–217 | match (nonzero established) |
| (5b) exact closed-form value of that determinant | line 215 literal `(200+147π²)` vs notes `(251+215π²)` | mismatch (notes only) |
| (6) i=0/h=0/m=0 survivors, σ_i/σ_h/σ_m, σ_transfer, ceilings | lines 225–291 | match |
| Mathematica second-engine derivation | (none) | missing |

`paper_alignment: partial` — the stage card and appendix are fully aligned with the script; the only paper-side conflict is the determinant's *exact* form, which appears only in the notes (the .tex states "nonzero," which the script does verify), and the second engine is absent.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 147 | `M_transfer.rank() == 3` | corridor setup (input) | yes |
| A2 | sympy | 150 | `len(transfer_null) == 2` | corridor 2D (deliverable 0/input) | yes |
| A3 | sympy | 157–160 | `assert_close` basis vs `expected_t` | corridor basis | yes (independent nullspace) |
| A4 | sympy | 171 | `simplify(Xi_transfer − 2(P01/P − Δ01/Δ)) == 0` | claim 1 | yes (N01 diff'd independently) |
| A5 | sympy | 176 | `simplify(Λ − G_W/Ω_W²·(1+I)/(1−H)) == 0` | claim 2 | yes (algebraic identity) |
| A6 | sympy | 187 | `simplify(Xi_transfer − Xi_mih) == 0` | claim 4 | yes |
| A7 | sympy | 190 | `Xi_sample_specialized == 0` | claim 4 (sample law) | yes |
| A8 | sympy | 192–193 | `I_s == 3/16`, `H_s == 25/(98)/π²` | claim 4 (coeffs) | yes |
| A9 | sympy | 216 | `simplify(det_ih − expected_det) == 0` | claim 5b (exact value) | yes but value disputed by notes |
| A10 | sympy | 217 | `det_ih != 0` | claim 5a (no-go) | yes |
| A11 | sympy | 237–242 | `assert_close` v_i/v_h/v_m vs stored | claim 6 (survivors) | yes (independent nullspaces) |
| A12 | sympy | 250–252 | `assert_close(abs(Xi_v…), …)` | claim 6 (gain scales) | yes |
| A13 | sympy | 269 | `assert_close(σ_transfer, 2.31561904…)` | claim 6 (norm) | yes (projector-derived) |
| A14 | sympy | 288–291 | `assert_close` four ceiling pairs | claim 6 (ceilings) | yes |
| — | mathematica | — | (no script) | all | missing |

All "yes" rows are non-tautological: every numeric expectation is checked against a value the engine derives independently in-script (differentiation, `nullspace()`, `det()`, projector norm), not against the same literal used to build it.

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.md:294`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.py:215`

**What's wrong:**
The notes state the combined `i=h` rigidity determinant as
> (notes:294) `-\frac{19(-25+98\pi^2)(251+215\pi^2)(441\pi^2+4400)}{6(8670000+14894275\pi^2+2117682\pi^4)}`

while the script's `expected_det` (which the run confirms via `assert simplify(det_ih − expected_det) == 0` and the saved output line 24 reproduces) is
> (script:215) `-19*(-25 + 98*pi**2)*(200 + 147*pi**2)*(441*pi**2 + 4400)/(6*(8670000 + 14894275*pi**2 + 2117682*pi**4))`

Every factor matches except the middle one: notes `(251+215π²)` vs script `(200+147π²)`. The difference is exactly `51 + 68π²`. This is the same notes-side drift pattern flagged for sibling stages 222/223 (notes = script + 68 in the π² coefficient): here the π² coefficient differs by exactly 68 (215 vs 147) and the constant by 51. The script's value is what SymPy actually computes (confirmed in the committed output), and the determinant is independently derived in-script from `red_i`, `red_h`, and the nullspace basis `T`, so the script side is internally self-consistent. The stage card and the part-07 appendix do **not** state the determinant's exact form (they only assert "nonzero," which the script does verify at line 217), so the conflict is confined to the notes prose.

**Why this matters:**
The no-go conclusion (`det ≠ 0 ⇒ only trivial pure-transfer drift`) is unaffected — both candidate polynomials are manifestly positive-definite in π², so the determinant is nonzero either way. But the notes are the authoritative detail source for this stage, and a reader comparing the notes' closed form against the script (or against a future Mathematica re-derivation) will hit a contradiction. Per protocol this is routed to the user; the auditor does not silently pick a side.

**Required change:**
None by Codex. See `## Resolve before fix_loop` in the directive. Codex must not edit notes/ or the script literal until the user chooses a direction.

**Verification:**
After user resolution: if the script is correct (expected, given the committed run confirms it), the notes line 294 middle factor is corrected to `(200+147π²)` by Claude (notes are prose, Claude-owned per the file-ownership policy) and no script change is made; the script continues to exit 0. If the notes are correct, the script's `expected_det` literal at line 215 is changed and re-run, but this would only pass if SymPy actually produces `(251+215π²)` — which the current committed output (line 24) shows it does not.

### F2 — missing_verification_script

**Severity:** medium
**Subtype:** missing_mathematica
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `*stage227*` file exists)

**What's wrong:**
There is no Mathematica audit script for unit 227. The stage card itself records "Mathematica audit: none yet" (`stage_227.tex:11`). The dual-engine rule requires a `.wl` wherever Mathematica *can* independently verify the stage — the test is possibility, not necessity. Stage 227 is entirely closed-form symbolic linear algebra over exact rationals and `π`: differentiation of rational functions w.r.t. a scalar deformation parameter, a 3×5 mixed matrix and its 2D null space, a 2×2 reduced determinant, one-dimensional null spaces of the reduced interference/hybridization/mixed-leg rows, and a Gram-projector operator norm. Every one of these maps onto native Mathematica primitives (`D`, `Solve`/`NullSpace`, `Det`, `FullSimplify`, `Eigensystem`/`Norm`, exact `Pi`). Independent verification is clearly POSSIBLE, so the absence is a finding.

**Why this matters:**
Stage 227 carries the first exact same-charge co-loading no-go (an MTDC-T11.3 anchor cited downstream by the same-charge sieve, rigid-mouth normal form, and twin-support branch per `stage_227.tex:26`). A single-engine result for a load-bearing no-go has no cross-check; the F1 notes/script determinant disagreement is precisely the kind of thing a second engine would catch.

**Required change:**
Codex authors the `.wl` per the claim manifest in the directive, deriving each result with native Mathematica primitives and a *different decomposition* than the SymPy script (anti-transliteration guard spelled out in the directive). Do not port the `.py` line-by-line.

**Verification:**
`redteam exec-mathematica 227` produces the new `.wl`'s output; all in-file checks (M1–M7 below) pass and the script exits 0. The verifier confirms the determinant's nonzero conclusion and the survivor gain scales match the SymPy values, and confirms the file is an independent derivation (not a transliteration).

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. The directive's claim manifest mandates an independent decomposition (anti-transliteration guard) so the future `.wl` is not a port.

## Engine cross-check

Only SymPy is present; no cross-check possible (`engines_agree: n/a`). The SymPy output is fresh (output mtime 2026-05-11 12:51 > script mtime 2026-05-11 11:58) and reports `EXIT_CODE: 0` with all checks passing.

## Verdict justification

The SymPy script is a strong, non-tautological audit: every numeric expectation is checked against an independently derived quantity, the hand-written `m/i/h` linear forms (lines 179–181) correctly match the `δln` of the definitions, the quotient-rule theorem genuinely re-derives `N01` rather than assuming it, and the factorization/sample-law/survivor/ceiling checks all trace to a paper deliverable. I tried to break the determinant assertion (A9) as a hardcoded_result — it is not, because `det_ih` is computed from the reduced rows and nullspace and the literal is what SymPy actually emits (output line 24). I tried to break the pure-transfer theorem (A4) as tautological — it is not, because `N01_mixed` is an independent derivative of `N0d`. Two findings stand: (F1) a notes-vs-script disagreement on the determinant's middle factor `(251+215π²)` vs `(200+147π²)`, matching the sibling 222/223 +68-in-π² notes-drift pattern, routed to the user; and (F2) the absent second engine for a stage where Mathematica verification is plainly possible. Verdict: `findings`, no stop-cold (the no-go conclusion is robust to F1 since both polynomials are positive-definite, so no math is broken and nothing downstream is invalidated pending resolution).

## Self-test notes

Variable-independence: the only derivatives are `sp.diff(..., eps)` of dressed quantities that genuinely depend on `eps` through every `exp(eps·x)` factor, so no derivative is identically zero. Parity/integral traps: no integrals in this stage — all checks are algebraic/linear-algebraic. Trivial-case: the `det_ih != 0` conclusion holds for both candidate polynomials since `200+147π²` and `251+215π²` are each strictly positive; I confirmed the `m/i/h` linear forms reproduce the correct `δln` of `G_W/Ω_W²`, `R G_U/(Ω_U²G_W)`, and `R²/(Ω_U²Ω_W²)`. Path spec: the F2 target `.wl` path lives under `mathematica/` and is named verbatim in the directive. Paper round-trip: F2 prescribes no constant that conflicts with the paper; F1 is held for user resolution and prescribes no Codex edit.
