---
unit_id: 167
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 167

## Per-finding outcomes

### F1 — hardcoded_result (mislabeled banner)

**Classification:** resolved

**What changed:**
The banner literal was corrected in both engines:
- `scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py:29` — `banner("STAGE 150 — …")` → `banner("STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION")`.
- `mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl:26` — `banner["STAGE 150 — …"];` → `banner["STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION"];`.

The git diff shows exactly one line changed per file — no collateral edits. The inline stage-number comments (sympy L33/L52/L81), explicitly placed out of scope by the directive, were correctly left untouched.

**Assessment:**
The edit matches the directive's required change verbatim. The corrected banner propagates into both regenerated transcripts: sympy output line 3 and mathematica output line 3 now read `STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION` (the report/directive said "line 11" but the banner sits at output line 3 in both transcripts; the substantive requirement — that the captured evidence read STAGE 167, not 150 — is satisfied unambiguously). This is a pure string/provenance fix with zero math impact, so there is no tautology hazard and no possibility of introducing a paper_misalignment. Finding fully closed.

## Exec log assessment

**SymPy:** exit=0 (inferred). The `expect_zero` harness `raise`s `AssertionError` on any non-zero; the transcript contains no traceback/AssertionError text and prints all six zero-checks cleanly: `delta ln r_c = 0`, `delta ln frak r = 0`, `delta ln frak g = 0` (lines 26-28), and the two channels + `delta_perp = 0` (lines 33-35). 7 `expect_zero(` occurrences in the .py minus the 1 definition = 6 live asserts, matching the report's inventory.

**Mathematica:** exit=0. 13 `PASS:` lines, 0 `FAIL`, and the closing line `Stage 167 Mathematica audit passed.` (output line 70). 14 `expectZero[` occurrences in the .wl minus the 1 function definition = 13 live asserts, all PASS. (The directive's verification note said "15 mathematica PASS lines"; the correct count is 13 — assertion-inventory rows A7-A15 expand to 13 individual `expectZero` calls, since A14 and A15 are 3-assert grouped rows. The "15" was an over-count in the note, not a missing assertion; coverage is complete.)

**Output freshness:** confirmed. Script mtimes are 15:54:48 (.py) / 15:54:49 (.wl); transcript mtimes are 16:10:14 (sympy .txt) / 16:11:29 (mathematica .txt). Both outputs are newer than their scripts, so the transcripts were regenerated after the banner fix.

## Material-change assessment

`material_change`: false. The only edit is a non-load-bearing display-banner string. No derived expression, assertion, or numeric result changed — all derived log-drifts, the three parent invariants, both off-family channels, and `delta_perp` are byte-identical to the pre-fix transcripts. No downstream unit depends on a print banner.

## Side observations (non-blocking)

- The directive's "line 11" / "15 PASS lines" figures in the verification note are both slightly off (banner is at output line 3; mathematica has 13 asserts). These are documentation inaccuracies in the directive, not defects in the applied fix; they do not affect closure of F1.
- The auditor noted but did not flag the inline stage-number comments referencing Stages 166/165/163 vs. the notes' 251/250/248. Correctly left out of scope; not a verification concern.

## Verdict justification

The single low-severity finding (mislabeled STAGE 150 banner) is fully resolved in both engines with a one-line edit per file and no collateral changes; the git diff confirms surgical scope. Both transcripts were regenerated post-fix (mtimes newer than scripts), now display STAGE 167, and both runs are clean — 6 SymPy zero-checks with no AssertionError and 13 Mathematica PASS with 0 FAIL plus the "audit passed" sentinel. The fix is provenance-only with no math impact. Verdict: verified.
