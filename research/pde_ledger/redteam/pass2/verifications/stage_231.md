---
unit_id: 231
batch: VII.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
overall_verdict: verified
material_change: false
scripts_changed_this_batch: false
diff_empty: true
relgate_sympy_exit: 0
relgate_mma_exit: 0
engines_agree: true
findings_disposition:
  - finding: F1
    type: paper_misalignment
    subtype: paper_missing_script_claim
    disposition: deferred
    target: PAPER_CLEANUP P4-51
    is_script_fix: false
---

# Stage 231 verification — pass 2 (zero-script-correction batch VII.2)

## Verdict

`verified` — `material_change: false`. No script was changed this batch. The
captured diff `redteam/pass2/exec_logs/stage_231_diff.patch` is empty (0 bytes);
Codex applied nothing. The reliability-gate re-run passed on both engines
(`relgate/sympy_231.txt` and `relgate/mma_231.txt`, both exit 0, byte-identical
to committed outputs).

## Finding disposition (one line each)

- **F1 — paper_misalignment / paper_missing_script_claim (low):** card
  `paper/stages/stage_231.tex:11` Verification line reads "Mathematica audit:
  none yet" while a passing, committed `.wl` exists and confirms every
  deliverable. This is the standing CARD-TEXT-LAG class — a stale STATUS
  annotation, NOT a value/identity mismatch — and is **deferred to
  PAPER_CLEANUP P4-51** per the standing user decision. It is not a script fix
  and is correctly left for the later paper pass (same situation as sibling
  Stage 230 F1). The stage still goes `verified`.

## Script sanity confirmation

- `relgate/sympy_231.txt` ends "All Stage 231 symbolic and numerical audits
  passed." (exit 0).
- `relgate/mma_231.txt` ends "All Stage 231 Mathematica audits passed." (exit 0),
  every M1–M7 line `PASS`.
- The pass-1 dF/dξ coefficient fix (240/189→189, with 121·ξ³) is reflected in
  both fresh outputs:
  - SymPy out L6: `... + 189*delta**2*xi + 72*delta**2 + 297*delta*xi**2 + 121*xi**3 ...`
  - Mathematica out L3: `M1 ... D[F,x] numerator polynomial factor = 72*d^2 + 81*d^3 + 189*d^2*x + 297*d*x^2 + 121*x^3`
  - notes L84/98/162 read 189/121; reconciled — fix held.
- Engines agree exactly: `R_* = 1.229255438463336`, product-law residual = 0 in
  both, sample R_flip/R_den rows match to 9-digit precision (sympy out L30-32 vs
  mma out L52-72). No engine disagreement.

## Note

This is a zero-script-correction unit: no `.py` or `.wl` edit, no Codex change,
empty diff, 0 Mathematica seats consumed beyond the orchestrator reliability
re-run. The lone finding is paper-side prose only and prescribes no script edit,
so no Codex change could introduce a new paper_misalignment.
