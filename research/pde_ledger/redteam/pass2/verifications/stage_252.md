---
unit_id: 252
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_patch: /var/projects/toy_physics/research/pde_ledger/redteam/pass2/exec_logs/stage_252_diff.patch (empty, 0 bytes)
relgate:
  sympy: /var/projects/toy_physics/research/pde_ledger/redteam/pass2/exec_logs/relgate/sympy_252.txt (exit 0, all checks passed)
  mathematica: /var/projects/toy_physics/research/pde_ledger/redteam/pass2/exec_logs/relgate/mma_252.txt (exit 0, M1-M9 all PASS)
  engines_agree: true
findings_count: 1
needs_user_resolution: true
---

# Verification — unit 252 (pass-2)

## Summary

Zero-script-correction stage. Codex applied no edits; the captured diff
(`exec_logs/stage_252_diff.patch`) is empty (0 bytes). The reliability-gate
re-runs both exit 0:

- SymPy (`relgate/sympy_252.txt`) ends "All symbolic and numerical checks passed."
- Mathematica (`relgate/mma_252.txt`) ends "All Mathematica checks passed for
  M1-M9." All M1-M9 sections PASS; the two nonzero residuals (M9 safe-prefactor
  2.84e-14 within tol 1e-9; M9 lattice-energy 8.67e-19 within tol 1e-12) are
  inside tolerance.

The engines agree exactly (f_vac/f_lat forms, drift law, I1/I2/I3, safe-energy
reduced form, 3:1 surface, and the Session-IV benchmark numerics all match).

## Finding disposition

### F1 — paper_misalignment (subtype: paper_missing_script_claim) — DEFERRED

- **Class:** card-text-lag. The stage-252 card verification line
  (`paper/stages/stage_252.tex:4`) still reads "Mathematica audit: none yet,"
  but a complete passing Mathematica audit `.wl` (M1-M9) now exists.
- **Disposition:** DEFERRED to the paper-side doc-sync sweep (P4-51 class).
  Not a script defect — the math is sound on both engines. Codex applies
  nothing on this unit; resolution is a user-gated paper-side edit (replace the
  verification line to cite the `.wl`). `needs_user_resolution: true`.

## Pass-1 round-trip re-check (the load-bearing concern)

Confirmed the pass-1 `gamma_safe_eq` X−X tautology is resolved. The
load-bearing rate identity is now `safe_combo/sc == Gamma3 sc^2 + Gamma5 sc^4`
(sympy L140 / wl M7-L182, "M7 safe rate raw quotient = 0"), a genuine algebraic
identity that fails if the kernel powers are wrong (sc^3/sc → sc^2,
sc^5/sc → sc^4). The downstream bridge (sympy L148 / wl M7-L185) is connective —
it documents the safe-edge reduction to mu_eta(s0^2-sc^2)/sc and fires
genuinely (the `G3 sc^3 + G5 sc^5` Add survives as a structural node under
division by sc). Other attacks (variable-independence on `D[flat,rV]`, subs-fired
check on the safe-rule substitution, benchmark calibration pins) all cleared in
the audit report and are consistent with the relgate output.

## Verdict

`verified` — `material_change: false`. No script change; both engines pass at
exit 0 with byte-identical-to-report results. The sole finding is a
low-severity card-text-lag, deferred to the P4-51 doc-sync sweep pending user
direction.
