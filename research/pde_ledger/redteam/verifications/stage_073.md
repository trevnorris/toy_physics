---
unit_id: 073
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 073

## Per-finding outcomes

### F1 — tautological_check (Mathematica precedence bug)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:57` now reads
```
expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];
```
The parentheses bind `- 37` outside the Rule, so `len/ell -> 37` substitutes first and the residual is `37 - 37`.

**Assessment:**
Edit matches the directive exactly (diff line `+expectZero["eta(reference) - 37", (eta /. (len/ell) -> 37) - 37];`). The saved Mathematica output now prints `eta(reference) - 37 = 0` with `PASS: eta(reference) - 37`. Combined with F2's edit, `eta` is now `len/ell` (the simplified output shows `eta under K_m = T_X/ell -> len/ell`), so the substitution and subtraction operate on the genuine reference number. No collateral edits.

### F2 — tautological_check (eta - L/ell cancellation)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:60-62` introduces `eta_sym = K_m * L / T_X` and then `eta = sp.simplify(eta_sym.subs(K_m, T_X / ell))` before the assertion.
- `mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:53-54` introduces `etaSym = km*len/tx;` and `eta = FullSimplify[etaSym /. km -> tx/ell, ...]` before the assertion.

**Assessment:**
Both engines now build `eta` symbolically in `K_m`/`km` first, then perform the Robin-closure substitution. The assertion `eta - L/ell == 0` now exercises the substitution rather than a trivial cancellation of `T_X`. Saved sympy output: `eta under K_m = T_X/ell -> L/ell`, `eta - L/ell = 0`. Saved Mathematica output: `eta under K_m = T_X/ell -> len/ell`, `eta - L/ell = 0`, `PASS: eta - L/ell`. Edit matches the directive verbatim. No collateral edits.

### F3 — hardcoded_result (geometry-map literals lack symbolic identity)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.py:36-42` adds a symbolic block defining `L_sym, a_sym, ell_sym`, `Lambda_star_sym = L_sym/a_sym`, `ell_over_a_sym = ell_sym/a_sym`, `Lambda_ell_sym = sp.simplify(Lambda_star_sym/ell_over_a_sym)`, and `expect_zero("Lambda_ell - L/ell (symbolic)", Lambda_ell_sym - L_sym/ell_sym)`.
- `mathematica/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.wl:28-35` adds a `Module[...]` block defining `lambdaStarSym = lSym/aSym`, `ellOverASym = ellSym/aSym`, `lambdaEllSym = FullSimplify[lambdaStarSym/ellOverASym, ...]`, and `expectZero["Lambda_ell - L/ell (symbolic)", lambdaEllSym - lSym/ellSym]`.

**Assessment:**
Both new blocks precede the existing numerical specialization, leaving the original `epsilon_r = 1/20` / `Lambda_star = 37/20` block intact. Saved sympy output shows `Lambda_ell - L/ell (symbolic) = 0` (line 5) before `Lambda_ell - 37 = 0` (line 10). Saved Mathematica output shows `Lambda_ell - L/ell (symbolic) = 0` and `PASS: Lambda_ell - L/ell (symbolic)` (lines 5-6) before `Lambda_ell - 37 = 0` and `PASS: Lambda_ell - 37` (lines 11-12). The new assertion is non-tautological: it tests the symbolic identity `(L/a)/(ell/a) = L/ell` independent of the chosen reference rationals. A perturbation to `2 L_sym/a_sym` would yield residual `L_sym/ell_sym ≠ 0`. Edit matches directive. No collateral edits.

### F4 — mathematica_transliteration (INFORMATIONAL)

**Classification:** resolved

**What changed:**
No files changed (the directive explicitly instructed Codex to skip F4). The directive's `## Applied: F4` block records `files_changed: []`, `summary: informational only — no changes per directive`, `deviation: none`.

**Assessment:**
Compliant with the directive. The informational-only nature is acknowledged. The Mathematica script remains a transliteration; this will be flagged in `MATHEMATICA_MIRROR_POLICY.md` per the post-batch tracker policy.

## Exec log assessment

**SymPy:** exit=n/a. The orchestrator-captured `redteam/exec_logs/stage_073_sympy.log` does not exist on disk. The canonical saved output `scripts/output/moving_throat_pde_stage073_family1_geometry_map_sympy_audit_refresh.txt` is present, was regenerated 2026-05-22 23:09:49 (after script mtime 23:08:30), and shows all four assertion lines passing (printed residuals all `0`; the script uses `raise AssertionError` on nonzero residual, so the clean output proves exit 0). Notable lines:
- `Lambda_ell - L/ell (symbolic) = 0`
- `Lambda_ell - 37 = 0`
- `eta under K_m = T_X/ell -> L/ell`
- `eta - L/ell = 0`
- `eta(reference) - 37 = 0`

**Mathematica:** exit=n/a. The orchestrator-captured `redteam/exec_logs/stage_073_mathematica.log` does not exist on disk. The canonical saved output `mathematica/output/moving_throat_pde_stage073_family1_geometry_map_mathematica_audit.txt` is present, was regenerated 2026-05-22 23:09:58 (after script mtime 23:08:40), and ends with `Stage 073 Mathematica audit passed.` (which the script prints only just before `Exit[0]`). Notable lines:
- `Lambda_ell - L/ell (symbolic) = 0` / `PASS: Lambda_ell - L/ell (symbolic)`
- `Lambda_ell - 37 = 0` / `PASS: Lambda_ell - 37`
- `eta under K_m = T_X/ell -> len/ell`
- `eta - L/ell = 0` / `PASS: eta - L/ell`
- `eta(reference) - 37 = 0` / `PASS: eta(reference) - 37`

**Output freshness:** Confirmed. Both `.txt` outputs (mtimes 23:09:49 and 23:09:58) are newer than their respective scripts (mtimes 23:08:30 and 23:08:40). The fix_loop re-ran both engines after the edits.

## Material-change assessment

`material_change`: false.

The unit's printed downstream ledger values (`Lambda_ell = 37`, `eta = 37`) are unchanged. The edits hardened the assertion residuals (made them non-tautological) and added a new symbolic identity check but did not alter any derived symbolic content or numeric value that downstream units consume. The intermediate symbol `eta` is now printed as `L/ell` (sympy) / `len/ell` (Mathematica) instead of being formed by direct `T_X/ell * L/T_X` cancellation, but the simplified value is the same.

## Side observations (non-blocking)

- The Mathematica output line `eta under K_m = T_X/ell -> len/ell` uses the local Mathematica symbol `len` rather than `L` — this is purely cosmetic (the .wl uses `len` because `L` collides with built-ins) and is consistent with the rest of the script.
- The orchestrator's exec_logs for this unit are absent on disk; I relied on the canonical saved `.txt` outputs (per the prompt's "if the captured exec_logs are missing... read the .txt outputs" guidance) and confirmed both engines emit their terminal "passed" / final-ledger markers, which the scripts only reach when all `expectZero`/`expect_zero` checks succeed.

## Verdict justification

All three mechanically applicable findings (F1, F2, F3) were applied verbatim against the directive's required-change blocks; the resulting `.txt` outputs are fresh, contain the new and updated assertion lines, and reach their respective end-of-script markers (which the scripts only emit on full pass). F4 was correctly skipped per directive instruction. No collateral edits or regressions in the diff, and the printed downstream values (`Lambda_ell = 37`, `eta = 37`) are unchanged.

stage 073: verified
