---
unit_id: 081
batch: III.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T22:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 081

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
A two-line `expect_zero` call was inserted into
`scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:40-41`,
immediately after the `print("Q(zeta;eps_blk) =", Q)` line (L39) and before the
existing `Q(0)-1` anchor check (now L42):

```python
expect_zero("Q-closedform",
    Q - (1 + zeta - 2*eps_blk*zeta)/(1 - eps_blk*zeta))
```

**Assessment:**
Correct and complete. The inserted RHS is the verbatim notes/paper closed form
`Q = [1 + (1-2 eps_blk) zeta]/[1 - eps_blk zeta]`; the constants `1`, `2`, and the
`(1-2 eps_blk)` coefficient are unchanged, and it exactly mirrors the full-identity
assertion the Mathematica engine already carries at `wl:54`
(`qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)`). The assertion is
**non-tautological**: the LHS `Q` is produced by `sp.solve` of the support-demand
law and then `sp.simplify`'d (out L11: `(2*eps_blk*zeta - zeta - 1)/(eps_blk*zeta - 1)`),
while the RHS is an independently-written literal closed form; the residual reduces to
`0` only because the two forms are genuinely identical (numerator and denominator both
sign-negated). This closes the coverage gap: SymPy now load-bearingly verifies the
stage's primary deliverable (the full closed-form inversion) across `zeta` and
`eps_blk`, not merely at the two interpolation points `zeta=0,1`. No collateral edit:
the captured diff is exactly the 2-line insertion, nothing else touched.

### F2 — numbering cross-reference (DEFERRED, not in scope)

Not a categorized finding for this loop. The directive explicitly deferred the two
comment-only cross-stage labels (`:32` `Stage-35`, `:46` `Stage 63`) to the dedicated
numbering pass because they are CROSS-references (to stages 052/080), not stage-081
self-labels. Confirmed unchanged in the diff and present in the current file (L32, L46)
as expected. Their presence is intentional and is not treated as a defect.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Q-closedform = 0` (new assertion passes — `expect_zero` would have raised otherwise)
- `Q(0)-1 = 0`, `Q(1)-2 = 0` (pre-existing anchors still pass)
- `Q(zeta;eps_blk) = (2*eps_blk*zeta - zeta - 1)/(eps_blk*zeta - 1)` (solved form, the LHS of the new check)
- `# exit_code: 0`

**Mathematica:** exit=n/a — no `.wl` change was required or made (the Mathematica engine
already asserted the full identity at `wl:54`; F1 only brought SymPy up to parity). Not
re-run for this finding.

**Output freshness:** confirmed. The committed
`scripts/output/...stage081...sympy_audit.txt` shows `Q-closedform = 0` (line 7) and
its mtime (2026-06-05 15:46:36) is newer than the edited script mtime
(2026-06-05 15:37:51), so it was regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit is a purely additive verification assertion — it
strengthens coverage of an already-correct, already-asserted (in Mathematica) identity
and changes no derived value, symbol, or output number. The `Q` closed form, the five
`Pi/C` thresholds, and the blocking ceiling are all unchanged in the output. No
downstream unit depends on anything that moved.

## Side observations (non-blocking)

None beyond the deferred F2 cross-references already accounted for above.

## Verdict justification

The sole finding F1 is resolved: the directed `expect_zero("Q-closedform", ...)`
assertion was inserted verbatim with the exact notes/paper closed form, the SymPy run
exits 0 with the new `Q-closedform = 0` line passing, the committed output was freshly
regenerated, and the captured diff shows no collateral edits. The new check is
non-tautological (solve-derived LHS vs. literal RHS) and brings SymPy to parity with the
Mathematica engine. The deferred F2 numbering cross-references are intentionally untouched
per orchestrator policy. Verdict: verified.
