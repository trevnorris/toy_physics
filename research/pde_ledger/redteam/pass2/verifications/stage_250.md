---
unit_id: 250
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-10T00:00:00Z
overall_verdict: verified
material_change: false
codex_edits: none
diff_empty: true
relgate:
  sympy: exit_0_all_checks_passed
  mathematica: exit_0_all_checks_passed
  byte_identical_rerun: true
findings_count: 1
findings_disposition:
  - id: F1
    class: card_text_lag
    subtype: paper_missing_script_claim
    severity: low
    status: deferred
    routed_to: P4-51
    note: >-
      Card stage_250.tex:4 \stagefield{Verification} reads "Mathematica audit: none yet"
      but a full independent passing .wl (M1–M7 all PASS) exists. Paper-side text update,
      user-gated; NOT a script fix. Stage itself verifies. Deferred to the dedicated
      numbering/coverage pass P4-51.
load_bearing_claim_recheck:
  claim: global goldilocks-window / monotonicity (one-sided safe band)
  established_globally: true
  mechanism: Resolve[ForAll] / Reduce over the full positive domain (M1 global dt_cross/dE<0,
    M3 unique positive-domain edge + Reduce, M4 S<1 iff E>E_edge globally)
  point_evaluation: false
  pass1_fix_intact: true
transliteration_check: independent (.wl uses quantifier elimination ForAll/Reduce with no
  .py counterpart for the M4 biconditional window certificate; not a transliteration)
---

# Verification — unit 250

## Summary

ZERO-SCRIPT-CORRECTION unit. Codex applied no edits; captured diff
`redteam/pass2/exec_logs/stage_250_diff.patch` is empty (0 bytes, confirmed).
Reliability-gate re-run passed clean on both engines, byte-identical to the
committed outputs.

## Relgate confirmation

- SymPy (`redteam/pass2/exec_logs/relgate/sympy_250.txt`): ends
  "All symbolic and numerical checks passed." — exit 0.
- Mathematica (`redteam/pass2/exec_logs/relgate/mma_250.txt`): M1–M7 each PASS,
  ends "All Stage 250 Mathematica checks passed." — exit 0.

## Load-bearing claim re-check (holds)

The global goldilocks-window / monotonicity claim (one-sided safe band:
`t_collapse` energy-independent AND `t_cross(E)` monotone decreasing ⟹ a single
lower edge `E_edge`) is established GLOBALLY, not by a single sample point:

- M1: `Resolve[ForAll[..., D[t_cross, En] < 0]]` over the full positive domain
  → "M1 global dt_cross/dE < 0 = True".
- M3: `Reduce` for the edge + `Resolve[ForAll[..., S^2==1 ⟹ En==E_edge]]`
  → "M3 unique positive-domain edge = True".
- M4: `Resolve[ForAll[..., Equivalent[S<1, En>E_edge]]]`
  → "M4 S(E)<1 iff E>E_edge globally = True".

The pass-1 single-sample-point → global fix is present and intact. This is a
genuine quantifier-elimination certificate over the whole positive domain, not a
point evaluation.

## Transliteration check

The `.wl` is an independent re-derivation. M4's global window-equivalence
certificate has no counterpart in the `.py` (SymPy only asserts the edge formula
and a pointwise derivative sign). No `mathematica_transliteration` finding.

## Finding disposition

F1 (paper_misalignment / paper_missing_script_claim, low): the card
`\stagefield{Verification}` line says "Mathematica audit: none yet" while a full
passing `.wl` exists. This is the CARD-TEXT-LAG class — a paper-side, user-gated
text update, NOT a script fix. The stage itself verifies. Deferred to P4-51.

## Verdict

overall_verdict: verified — math holds, both engines agree, global certificates
intact, the lone finding is a deferred card-text-lag with no script impact.
