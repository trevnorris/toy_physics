---
unit_id: 050
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T18:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 050

## Per-finding outcomes

The original audit raised two informational stale_output findings (F1 = stale
banner/transcript; F2 = stale numbering self-labels in the SymPy script). The
orchestrator folded these into a single directive finding (F1) covering the two
unambiguous SELF-labels (docstring + close line), with the import-comment
cross-ref explicitly DEFERRED. I verify against the directive's F1.

### F1 — stale_output (unambiguous self-labels; number-only)

**Classification:** resolved

**What changed:**
Per `stage_050_diff.patch`, exactly two single-line edits in
`scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`:
- `py:3` docstring: `Stage 33 SymPy audit.` → `Stage 50 SymPy audit.`
- `py:118`: `print("\nAll Stage-33 symbolic checks passed.")` →
  `print("\nAll Stage-50 symbolic checks passed.")`
Current file confirms `py:3` now reads `Stage 50 SymPy audit.`. The committed
SymPy transcript (`...sympy_audit.txt:42`) now closes `All Stage-50 symbolic
checks passed.` and (`:3`) banners `STAGE 50 — PHYSICAL ZETA VS ZETA_REQ`.

**Assessment:**
Correct and complete. The change is label-only — strip the stage number from
both sides and the lines are identical to HEAD. No assertion, symbol, numeric
expression, import, or control flow was touched (diff is exactly the two label
lines; no other hunks).

(a) Directive's required change applied — yes, both targeted lines.
(b) Label-only / strip-the-number identical to HEAD — yes; the diff contains
    only `33`→`50` substitutions inside otherwise-identical text.
(c) Deferred cross-ref at `py:61` (`# Imported from Stage 32's explicit D/N
    overlap extraction ...`) was correctly LEFT UNTOUCHED — confirmed present
    verbatim. It is a cross-reference to the imported stage (049), not a
    self-label, and is deferred to the dedicated SCRIPT/OUTPUT-band numbering
    plan; its presence is correct, not a defect. The import statement at
    `py:17` (`...stage049_dn_overlap_zeta_sympy_audit import twin_support_ratio`)
    is unchanged and already canonical.
(d) No math/assertion/value changed — confirmed; both engines reproduce
    identical symbolic results (zeta_req, S, zeta_n, x_max, S_n^(max),
    S_1^(max)=2(eps-5)/(eps-9), S_2^(max)=2(eps-13)/(eps-25)) with all
    residuals 0 / all PASS.
(e) material_change false — confirmed.

`grep -nE "Stage 33|Stage-33"` on the py now returns nothing.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 50 — PHYSICAL ZETA VS ZETA_REQ` (banner)
- `zeta_0^(twin) - 1 (anchors doubling to stage 049 import) = 0`
- `x_max - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)] = 0`
- `All Stage-50 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- `STAGE 050 — PHYSICAL ZETA VS ZETA_REQ` (banner)
- `PASS: zeta_n^(twin)(x_max) - zeta_req`
- `PASS: Numerator of (S_n^(max) - S_n^(twin)) - (1-eps)(2n+1)^2 n(n+1) x`
- `Stage 050 Mathematica audit passed.`

**Output freshness:** confirmed. Both committed `.txt` outputs have mtime
1780683731, newer than the py (1780681733) and the wl (1780523951). Logs
are dated 2026-06-05, post-fix.

## Material-change assessment

`material_change`: false. The edit is a pure label correction in the SymPy
script's docstring and closing print; no derived result, assertion, or numeric
value changed. The Mathematica engine still independently re-derives the same
quantities and agrees. No downstream unit is affected.

## Side observations (non-blocking)

None. The deferred import-comment cross-ref ("Stage 32") at py:61-62 is the
known SCRIPT/OUTPUT-band item already tracked in the dedicated numbering plan;
flagging only for traceability, not as a defect.

## Verdict justification

The single directive finding F1 is resolved: the two stale "Stage 33"
self-labels are now "Stage 50", the change is label-only and identical to HEAD
modulo the stage number, the deferred cross-ref at py:61 and the canonical
import at py:17 are untouched, no math/assertion/value changed, both engines
exit 0 with refreshed transcripts whose mtimes postdate the scripts, and the
grep backstop confirms no residual stale self-labels. material_change is false.
Verdict: verified.
