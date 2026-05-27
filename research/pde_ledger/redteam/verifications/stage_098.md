---
unit_id: 098
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 098

## Per-finding outcomes

### F1 — paper_misalignment

**Classification:** resolved

**What changed:**
Per orchestrator's `## Applied: F1` block (Cluster C direction (a) from `redteam/resolutions/batch_IV1_paper_alignment.md`), script-side banners and docstrings were updated to `STAGE 098` to match the audit-unit number:
- `scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:3` now reads `Stage 098 SymPy audit: ...`
- `scripts/moving_throat_pde_stage098_family1_support_is_automatic_sympy_audit.py:42` now reads `print('\nSTAGE 098 AUDIT PASSED')`
- `mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:38` now reads `banner["STAGE 098 — FAMILY-1 SUPPORT IS AUTOMATIC"];`
- `mathematica/moving_throat_pde_stage098_family1_support_is_automatic_mathematica_audit.wl:74` already read `Print["Stage 098 Mathematica audit passed."];` (was already 098 pre-fix).

**Assessment:**
All four script-side stage-number labels now uniformly read `Stage 098 / STAGE 098`, eliminating the original `Stage 81 / 081` divergence on the script side. The paper card title "Stage~115" and notes "Stages 80–81" reference are deliberately deferred to a future paper-cleanup pass per the resolution document; the verifier scope is scripts-only, so the script-side relabeling is complete for the verifier's lane. The exec logs reflect the new banners (`STAGE 098 — FAMILY-1 SUPPORT IS AUTOMATIC`, `STAGE 098 AUDIT PASSED`, `Stage 098 Mathematica audit passed.`). No collateral edits beyond the four cited sites.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/.../stage098_..._sympy_audit.py:24-25` now carries a two-line provenance comment naming `notes/stages/moving_throat_pde_stage098_family1_support_is_automatic.md` as the source of `zeta_max^(F1) = 2.46752922945601`, directly above the unchanged `zmax_F1 = sp.N('2.46752922945601')` line (now line 26).
- `mathematica/.../stage098_..._mathematica_audit.wl:61-63` carries the parallel `(* ... *)` Mathematica comment, directly above the unchanged `zMaxF1 = SetPrecision[2.46752922945601, 20];` line (now line 64).

**Assessment:**
Matches the directive's "Required change" verbatim. The numeric literal `2.46752922945601` is preserved unchanged, as are the `expectApprox` targets `0.456730991107963169017835980412` and `2.01079823834804688464927835412`. The exec-log outputs remain bit-identical to pre-fix in their numeric content. The provenance comment is informational (not a derivation), which is consistent with the directive's "informational unless the carry-forward source cannot be located" language.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/.../stage098_..._sympy_audit.py:30-34` adds four new lines after the existing `assert gap_F1 > 0`:

```
# Numeric pin (matches Mathematica's expectApprox targets):
zeta_edge_F1_target = sp.Float('0.456730991107963169017835980412', 30)
gap_F1_target = sp.Float('2.01079823834804688464927835412', 30)
assert abs(zeta_edge_F1 - zeta_edge_F1_target) < sp.Float('1e-15', 30)
assert abs(gap_F1 - gap_F1_target) < sp.Float('1e-15', 30)
```

The existing `assert gap_F1 > 0` line is kept.

**Assessment:**
Matches the directive's "After" block verbatim. The two new assertions are non-tautological: they compare `zeta_edge_F1` (derived from `zmax_F1` via the symbolic `zeta_edge.subs(zmax, zmax_F1)` formula) against an independently-quoted literal from notes §3 / Mathematica's `expectApprox` targets. A typo in `2.46752922945601` would now produce a mismatch larger than `1e-15` and the script would fail — confirming the original audit's intent. The SymPy exec log shows the printed values `Family-1 zeta_edge = 0.456730991107963169017835980412` and `Family-1 margin = 2.01079823834804688464927835412` matching the pinned targets exactly to 30 digits.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `zeta_edge = zmax/(3*zmax - 2)`
- `zmax - zeta_edge = 3*zmax*(zmax - 1)/(3*zmax - 2)`
- `Family-1 zeta_edge = 0.456730991107963169017835980412`
- `Family-1 margin = 2.01079823834804688464927835412`
- `STAGE 098 AUDIT PASSED`

All 6 sympy asserts (lines 12, 15, 21, 29, 33, 34) must have passed for the final banner to print; the final banner is present. The numeric pin targets are reproduced exactly in the printed output.

**Mathematica:** exit=0. Notable lines (6 PASS lines, matching the 6 expect* calls at WL lines 56-59 and 70-71):
- `PASS: zeta_req - 1/(3-2 eps)`
- `PASS: d zeta_req / d eps exact formula`
- `PASS: gap factorization`
- `PASS: automatic-support gap is positive for zmax > 1`
- `PASS: Family-1 zeta_edge numeric check` (diff ≈ 4.11e-18)
- `PASS: Family-1 margin numeric check` (diff ≈ 4.11e-18)
- `Stage 098 Mathematica audit passed.`

PASS-line count (6) matches expect-call count (6). Both diffs are well below the 1e-15 tolerance.

**Output freshness:** confirmed. Sympy script mtime 1779902165 < sympy output mtime 1779913721; Mathematica script mtime 1779902169 < Mathematica output mtime 1779913837. Both outputs were re-generated after the fix.

## Material-change assessment

`material_change`: false.

The edits are: (a) docstring/banner string relabeling, (b) inline provenance comments, (c) added SymPy assertions that pin numeric values already produced. No derived numeric or symbolic quantity changes. No downstream unit depends on a numerical or symbolic result from this stage that wasn't already present pre-fix.

## Side observations (non-blocking)

- The Mathematica `expectApprox` target `2.01079823834804688464927835412` differs from the actually-computed `gapF1` value at the 17th digit (the printed `2.01079823834804688876172914845137997165` shows the rounding of the F1 propagation in 20-digit precision). The diff is `4.11e-18`, well below the `1e-15` tolerance, so this is not a finding — just noting that the targets pinned in both engines round at ~16 decimal places and the 30-digit string is partly cosmetic. The SymPy pin uses the same target literal and passes the same way.
- The paper card title "Stage~115" is explicitly out-of-scope for this verifier (paper/notes excluded by scope and by Cluster C direction (a)); the orchestrator's `## Applied: F1` block correctly flags it for a future paper-side pass via `PAPER_CLEANUP_TRACKER`.

## Verdict justification

All three findings are `resolved`. F1's script-side relabeling is complete and consistent across all four cited script sites; the paper-side deferral is documented in the `## Applied: F1` block per the user-resolved Cluster C direction. F2's provenance comments are in place at both engine sites with no change to the numeric literal. F3's SymPy numeric pin mirrors the Mathematica `expectApprox` targets and is non-tautological — it would catch a typo in `zmax_F1`. Both exec logs exit 0, both are fresher than their scripts, and the Mathematica PASS-line count (6) matches the expect-call count (6). No regressions visible in the current file state. Verdict: `verified`.
