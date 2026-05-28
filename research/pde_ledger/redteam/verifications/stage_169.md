---
unit_id: 169
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 169

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy `moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:133-144`: after the existing `print` of `Xi_perp_num` (line 131), Codex added per-coefficient extraction (`Xi_perp.coeff(XiT/Xiv/XiL)`), numeric substitution of `g_num`/`r_num`, and a loop that compares each computed coefficient to the paper literals `0.758035078944663` (XiT), `1.00314310113848` (Xiv), `1.88373219118005` (XiL), raising `AssertionError` if `|got - want| > 1e-12`.
- Mathematica `..._mathematica_audit.wl:105-119`: analogous block after the numeric print (line 103) using `Coefficient[xiPerp /. {g->gNum, r->rNum}, ...]` and the `fail/pass` helpers with the same three paper literals and `10^-12` tolerance.
- The pre-existing tautological `eps_perp - Xi_perp Igrp == 0` line was left in place per the directive ("harmless").

**Assessment:**
Correct and addresses the finding. The new checks are genuine, non-tautological tolerance assertions: `got` is computed by extracting the symbolic coefficient from the independently-defined `Xi_perp` (`g`, `g+1/(2s)`, `2g+3/(4s)` with `s=sqrt(1+r^2)`) and substituting numeric `g`/`r`, while `want` is a fixed paper literal. An incorrect Stage 253 weighting coefficient would now make `|got - want|` exceed tolerance and fail. The exec logs confirm this is a real comparison: the SymPy transcript reports nonzero residuals `|diff| = 3.97e-15` (Xiv) and `4.46e-15` (XiL) and exactly `0` for XiT (paper quoted `g` to full precision) — exactly the signature of comparing an independently-computed value against a truncated paper literal, not an identity. Mathematica shows three PASS lines for the same checks. Tolerance `1e-12` correctly accommodates the paper's truncation of Xiv/XiL while still catching any real coefficient error. The edit matches the directive verbatim with no collateral changes.

### F2 — mathematica_transliteration

**Classification:** resolved (policy-accepted)

**What changed:**
No `.wl` edit. The directive's `## Applied: F2` records an orchestrator policy decision to accept the Mathematica file as a documented mirror (per MATHEMATICA_MIRROR_POLICY), consistent with the disposition applied to prior-batch `mathematica_transliteration` findings. The algebra in the mirror is correct; no independent re-derivation was forced.

**Assessment:**
This is a deliberate, legitimate policy disposition, not a skipped finding. The original report itself flagged F2 as "a policy call" with "no mechanical code edit mandated by this finding alone," and the directive instructed acceptance under the established mirror policy. The F1 numeric coefficient checks (now real assertions in both engines) restore substantive second-engine coverage of the load-bearing Stage 253 weights, mitigating the transliteration concern. Classified as resolved per the orchestrator's stated policy disposition.

### F3 — stale_output (banner mislabel)

**Classification:** resolved

**What changed:**
- SymPy line 30: `banner("STAGE 152 — ...")` → `banner("STAGE 169 — ...")`.
- Mathematica line 31: `banner["STAGE 152 — ..."]` → `banner["STAGE 169 — ..."]`.

**Assessment:**
Correct. The git diff shows exactly the single-line banner change in each file. Both regenerated transcripts now read `STAGE 169 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE` at line 3 (the banner emits a leading blank line, so the title lands on transcript line 3, not line 11 as the directive loosely estimated; the substance — correct stage number in the header — is satisfied). No collateral edits.

## Exec log assessment

**SymPy:** exit=0 (no traceback; all assertions printed `= 0` / matched, script ran to the carry-forward block). Notable lines:
- `I[x,y] - [4 a_x a_y + 4/5 b_x b_y] = 0`
- `b_x - 3 a_x = 0` ; `I_axis - (7/10) eps^2 x1 y1 = 0`
- `<Y20> = 0` ; `<Y20^2> = 1/(4*pi)` ; `linear coefficient in averaged log-observable = 0`
- `Xi_perp coeff on XiT = 0.758... (paper 0.758..., |diff| = 0)`
- `Xi_perp coeff on Xiv = 1.0031... (paper 1.0031..., |diff| = 3.97E-15)`
- `Xi_perp coeff on XiL = 1.8837... (paper 1.8837..., |diff| = 4.46E-15)`

**Mathematica:** exit=0 (final line `Stage 169 Mathematica audit passed.` and `Exit[0]`; all `expectZero` reported PASS, no FAIL/Exit[1]). Notable lines:
- `PASS: I[x,y] - [4 a_x a_y + 4/5 b_x b_y]`
- `PASS: b_x - 3 a_x` ; `PASS: I_axis - (7/10) eps^2 x1 y1`
- `PASS: average Y20` ; `PASS: average Y20^2 - 1/(4 pi)` ; `PASS: linear coefficient in averaged log-observable`
- `PASS: Xi_perp coeff on xiT` / `xiv` / `xiL`

**Output freshness:** confirmed re-generated post-fix. SymPy script mtime 2026-05-28 15:55:10, its output 16:10:16 (newer). Mathematica script mtime 15:55:21, its output 16:11:42 (newer). Both transcripts carry the corrected `STAGE 169` banner, confirming the runs reflect the post-edit scripts.

## Material-change assessment

`material_change`: false.

The edits are: (F1) added assertions that test pre-existing computed quantities against fixed paper literals — no derived value changed, the script merely now checks what it previously only printed; (F3) a cosmetic banner string. The numeric Family-1 combination forwarded downstream (`0.758.../1.003.../1.884...`) is unchanged — it was already printed correctly pre-fix. No downstream-consumed result was altered, so units > 169 are not substantively affected by these edits.

## Side observations (non-blocking)

- The directive estimated the corrected banner would land on transcript "line 11"; because `banner` prints a leading blank line, the title is actually on transcript line 3. Cosmetic only; the F3 substance (correct stage number) is met. Not a blocker.
- The SymPy `s_num` (line 129) is computed but unused (the coefficient checks substitute symbolic `r` and let SymPy form `s` internally). Harmless dead assignment; pre-existing, not introduced by this fix.

## Verdict justification

All three findings are closed: F1's vacuous block-4 coverage is replaced by genuine per-coefficient tolerance assertions in both engines (confirmed non-tautological by the nonzero ~4e-15 residuals against truncated paper literals); F2 is a documented, legitimate policy-accepted mirror with F1 restoring substantive second-engine coverage; F3's banner mislabel is corrected and propagated into fresh transcripts. Both exec logs exit 0 with all assertions passing, outputs are newer than their scripts, and the git diff shows only the directed edits with no regressions. Verdict: verified.
