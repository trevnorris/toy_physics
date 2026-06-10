---
unit_id: 249
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_empty: true
reliability_gate: passed
engines_agree: true
findings_count: 1
findings_resolved_in_script: 0
findings_deferred: 1
---

# Verification — unit 249 (pass-2)

## Summary

Zero-script-correction unit. Codex made no edits; the captured diff
`redteam/pass2/exec_logs/stage_249_diff.patch` is empty (0 bytes). The reliability
gate passed: both fresh relgate outputs exit 0 and end all-checks-passed
(`exec_logs/relgate/sympy_249.txt` → "All symbolic and numerical checks passed";
`exec_logs/relgate/mma_249.txt` → "All Mathematica checks passed"). The sole finding
is a card-text-lag-class `paper_misalignment`, deferred to P4-51 (not a script fix).

## Reliability gate (fresh relgate outputs)

- SymPy `exec_logs/relgate/sympy_249.txt`: exit 0; closes with "All symbolic and
  numerical checks passed." Symbolic identities (S1 subtracted-RHS = -2*S_cov,
  S2 closure factorization, S3/S4 Möbius forms) reduce as expected; benchmark
  values reconstructed (alpha_pk=0.663669919..., R_int=4.109209231..., abar=0.608549990...).
- Mathematica `exec_logs/relgate/mma_249.txt`: exit 0; closes with "All Mathematica
  checks passed." M1-M5 all PASS; M1-M4 symbolic residuals = 0; eta_h cancels (True);
  M5 benchmark diffs all ~1e-9 < 1e-7 tol; asymmetry ordering + positive packet True.
- Engines agree to displayed precision (alpha_pk, R_int, abar identical; all
  symbolic residuals 0). No engine_disagreement.

## Independence / round-trip re-confirmation

Confirmed the pass-1 round-trip fix holds: the Möbius inverses are genuinely
Solve-derived, not tautological. SymPy obtains `alpha_from_Rpk` and `abar`
via `sp.solve` and compares the solved result to `(R-1)/(R+1)` — a real inversion,
not `x==x`; the `eta_h` cancellation is a falsifiable `free_symbols` test. The
Mathematica `.wl` mirrors with an independent `Solve[...]` + `FreeQ`, reached via a
structurally distinct choreography (M2 defines `Hdot[s_]:=Gamma0+s Gamma1` directly,
M2/M3 fused with a fresh `RpkExpr`), so it is a genuine second engine, not a
transliteration. No surviving tautology or self-referential round-trip.

## Findings disposition

| # | class | severity | disposition |
|---|---|---|---|
| F1 | paper_misalignment (card-text-lag) | low | DEFERRED → P4-51 (paper-side; stage verifies) |

F1: the stage card's `\stagefield{Verification}` line at `paper/stages/stage_249.tex:4`
still says "Mathematica audit: none yet," but the passing independent `.wl` exists
(added pass-1). Pure card-text lag — no math impact, no script change. Routed to the
P4-51 paper-side card-text-lag sweep; the directive carries only a `## Resolve before
fix_loop` block and Codex made no edits. All 11 deliverable values reconcile against
card + notes at matching precision.

## Verdict

overall_verdict: verified. material_change: false. The math holds under attack: both
engines present, independent, and agreeing; all symbolic identities → 0; benchmark
reconstruction within tol; round-trip concern resolved. The lone finding is the
card-text-lag `paper_misalignment` deferred to P4-51 — the stage itself verifies.
