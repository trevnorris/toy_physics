---
unit_id: 240
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: false
diff_empty: true
relgate_sympy_exit: 0
relgate_mma_exit: 0
engines_agree: true
findings_count: 1
---

# Verification — unit 240 (pass-2, ZERO-SCRIPT-CORRECTION)

## Inputs reviewed
- report: `redteam/pass2/reports/stage_240.md` (verdict: findings, 1 low, paper-side)
- directive: `redteam/pass2/directives/stage_240.md` (applied: false; needs_user_resolution: true; Codex applies nothing)
- captured diff: `redteam/pass2/exec_logs/stage_240_diff.patch` — EMPTY (0 bytes), confirms no Codex script edit
- reliability-gate re-run: `redteam/pass2/exec_logs/relgate/sympy_240.txt`, `mma_240.txt` — both exit 0

## Disposition

### F1 — paper_misalignment (paper_missing_script_claim / card-text-lag)
- **Class:** CARD-TEXT-LAG. Card `paper/stages/stage_240.tex:11` says `Mathematica audit: none yet` while a committed, passing `.wl` (M1-M6) exists.
- **Disposition:** DEFERRED to user / paper batch P4-51. Documentation-side only; no result value affected; Codex correctly made no edit. Not a script fix.

### Variable-independence self-test-trap fix — CONFIRMED (re-check holds)
- The extracted weights are pulled from `Y_support`, which carries `Omega_Q` through its pole (`.wl` out L9-10 shows the explicit `Apart` partial-fraction with `omega±omegaQ` simple poles; `.py` multiplies by the pole denominator pre-limit). The fix exercises independence against a genuinely VARIABLE-bearing object, not a constant.
- Guards (py L115-116/119-120/124-125; wl L77/L84) force the pre-limit object to contain `Omega_Q`, so the subsequent `diff(c0_static,Omega_Q)==0` / `diff(c1_static,Omega_Q)==0` (sympy out L11-12; mma out L25-28) tests the real notes §2.1 claim that `Omega_Q` drops out of the static extraction — non-tautological.
- Values reconcile exactly: ρ_α = 4/3, ζ_req = 1/3, Π_tr = (4/3)C_mix, S_req = 4/3, ϱ = 2(1-ε_*)/3, regime C_mix < Π_tr < 2C_mix. All MATCH across card / notes / appendix and both engine outputs.

## Engine sanity
- sympy_240.txt: 23 `[ok]` lines + "All Stage 240 symbolic checks passed", exit 0.
- mma_240.txt: PASS M1-M6 (all sub-checks) + "Stage 240 Mathematica audit passed", exit 0.
- Both fresh relgate outputs end all-checks-passed; engines agree (extracted c0 = alphaMix/alphaReq, c1 = 1 - alphaMix/alphaReq mirrors the SymPy extraction). `.wl` is NOT a pure transliteration (Apart route vs direct pole-multiply) — no transliteration finding.

## Verdict
- overall_verdict: **verified**
- material_change: **false** (zero Codex edits, empty diff, relgate byte-identical & exit 0)
- Only finding is paper-side card-text-lag → deferred to P4-51; F1 self-test-trap fix confirmed holding.
