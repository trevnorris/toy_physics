---
unit_id: 099
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 099

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Both engines replaced the four `R_n - (N_Q - 1) == 0` checks with three non-tautological structural identities driven by a FREE positive symbol `K0_sym` / `k0Sym`:
- sympy `scripts/.../_sympy_audit.py:63-80` introduces `K0_sym` and asserts (i) `K0 K4 - 4 K2^2 = 0`, (ii) `9*K2^{5/2}/K0^{3/2} - Gamma5_struct = 0`, (iii) `Gamma5_struct.subs(K0_sym, N_Q*K0_target) / (2G/(5c^5)) - N_Q = 0`.
- Mathematica `.../_mathematica_audit.wl:45-56` mirrors the same three identities with `k0Sym` as a free positive symbol.

**Assessment:**
Correct and matches the directive's verbatim code blocks. The new identities are non-tautological in the directive's sense: `K0_sym` is a free positive symbol (no upstream `K0_sym = N_Q * K0_target` substitution before the branch identity / sqrt-form checks). The branch identity `K0 K4 - 4 K2^2` algebraically tests the canonical `1/(4 Omega_Q^{2n})` structural ratios; the sqrt-form check tests the canonical odd-coefficient form; the third check exercises the appendix's `chi_Q=1` factorization. No collateral edits beyond the cited line ranges.

### F2 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
Per Cluster A direction (a) in `redteam/resolutions/batch_IV1_paper_alignment.md`, the sympy docstring (`_sympy_audit.py:3-19`) annotates the three card `\stagefield{Checks}` as upstream carry-ins (orthogonality at 094; static-limit at 091/092/094/096; minimal-module at 088/089/090) and notes that the local Yhat_Q^cons static-slot and pole-residue assertions act as substantive anchors here (the F4 implementation).

**Assessment:**
The orchestrator-direct Applied: F2 block legitimately closes this paper_misalignment via the documented user-resolution route. The docstring carry-in annotation is in place and the local-side substantive anchors land via F4. Matches the resolution policy.

### F3 — paper_misalignment (notes_contradicts_script)

**Classification:** resolved

**What changed:**
Per Cluster C direction (a), the script banners and final-line prints were relabeled to STAGE 099:
- sympy `_sympy_audit.py:3` ("Stage 099 SymPy audit..."), line 40 (`banner("STAGE 099 — REDUCED FINISH LINE")`), line 85 (`"STAGE 099 AUDIT PASSED"`).
- Mathematica `_mathematica_audit.wl:26` (`banner["STAGE 099 — REDUCED FINISH LINE"]`), line 62 (`"STAGE 099 AUDIT PASSED"`).
- Output transcripts echo the new "STAGE 099" labels (sympy txt:3,16; Mathematica txt:3,22).

**Assessment:**
Banner-number convention aligned to audit-unit ID 099 in both engines. The paper card title `Stage~116` is explicitly deferred to PAPER_CLEANUP_TRACKER (within-scope deviation noted in directive Applied block). No stale "STAGE 82/082" remains in either script or transcript.

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
Both engines now exercise the conservative module form `Yhat_Q^cons` non-trivially:
- sympy `_sympy_audit.py:49-57` asserts (i) `Yhat_Q_cons.subs(omega, 0) - 1 = 0` (static slot) and (ii) `sp.residue(Yhat_Q_cons, omega, Omega_Q) - (-Omega_Q/8) = 0` (pole residue).
- Mathematica `_mathematica_audit.wl:37-39` mirrors: `(yhatCons /. omega -> 0) - 1` and `Residue[yhatCons, {omega, omegaQ}] - (-omegaQ/8)`.

**Assessment:**
Both new assertions land before the `K0_target` checks and exercise the conservative module's partial-fraction structure directly — these would fail under perturbation (e.g., `1/2 + 1/(2(1 - ...))` gives residue `-Omega_Q/4`, not `-Omega_Q/8`). Non-tautological, exact match to the directive's verbatim code.

## Exec log assessment

**SymPy:** exit=0 (inferred from "STAGE 099 AUDIT PASSED" at line 16, no traceback). 6 expected `<name> = 0` lines counted at txt lines 6, 7, 8, 9, 10, 11 — matches 6 `expect_zero` calls in the script (static slot, pole residue, K0_target geom, branch identity, Gamma_5 sqrt form, Gamma_5 normalization). Notable lines:
- `Yhat_Q^cons pole residue at omega=Omega_Q is -Omega_Q/8 = 0` (txt:7)
- `branch identity K0 K4 = 4 K2^2 = 0` (txt:9)
- `Gamma_5 normalization equals N_Q on chi_Q=1 branch = 0` (txt:11)

**Mathematica:** exit=0 (inferred from final "STAGE 099 AUDIT PASSED"; no `FAIL:` lines and no `Exit[1]`). 6 expected PASS lines counted at txt lines 7, 9, 11, 13, 15, 17 — matches 6 `expectZero` calls. Notable lines:
- `PASS: Yhat_Q^cons pole residue at omega=omegaQ is -omegaQ/8` (txt:9)
- `PASS: branch identity K0 K4 = 4 K2^2` (txt:13)
- `PASS: Gamma_5 normalization equals N_Q on chi_Q=1 branch` (txt:17)

**Output freshness:** confirmed. sympy script mtime 1779902197 < sympy output 1779913722; Mathematica script mtime 1779902218 < Mathematica output 1779913844. Outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The replacement substituted tautological-by-construction checks with non-tautological structural identities that exercise the same canonical `1/(4 Omega_Q^{2n})` and `9 K_0/(32 Omega_Q^5)` relations the paper relies on. No derived numerical or symbolic result that downstream stages consume has changed — the algebraic content (Yhat_Q^cons structure, K0_target two-form, `Gammabar_5/(2G/5c^5) = N_Q` on the chi_Q=1 branch) is identical to before, just verified more rigorously. Banner relabeling and docstring annotation are cosmetic / documentation. No downstream-stale flag warranted from a math-content standpoint.

## Side observations (non-blocking)

- The card's title "Stage~116" remains misaligned with the audit-unit ID 099; this is explicitly tracked in PAPER_CLEANUP_TRACKER per Cluster C (a) and is intentional within this batch. Not a verification failure.
- The sympy transcript does not emit explicit "PASS:" lines (the sympy `expect_zero` only prints `<name> = <expr>` and raises on non-zero), unlike the Mathematica `expectZero` which prints both `<name> = ...` and `PASS:`. Cosmetic only; the absence of an AssertionError plus the printed `= 0` is equivalent to a pass.

## Verdict justification

All four findings (F1 tautological_check, F2 paper_misalignment via Cluster A (a), F3 paper_misalignment via Cluster C (a), F4 insufficient_verification) are resolved: the directive's verbatim code edits are in place at the cited line ranges in both engines; both transcripts show the expected 6 zero residuals / 6 PASS lines with no failures; outputs are fresh (mtimes post-edit); no regressions; no introduction of new tautologies (`K0_sym` is genuinely free, not pre-substituted). Material content is unchanged downstream.
