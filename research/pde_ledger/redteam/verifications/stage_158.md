---
unit_id: 158
batch: IV.6
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-28T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 158

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
Paper card `paper/stages/stage_158.tex:23-24` now reads:
- item 2: `Even-preservation of the canonical gate is verified downstream: see \ref{stage:159}.`
- item 3: `Tangent motion on the parent compensation family giving \(\delta_\perp=0\) is verified downstream: see \ref{stage:162} and \ref{stage:163}.`

Scripts untouched for F1 (per user-gate resolution direction Cluster C: forward-carry citations).

**Assessment:**
Per orchestrator note, the user resolution was Cluster C (forward-carry citations to downstream stages). The card now clearly tags items 2 and 3 as downstream-verified rather than claiming Stage 158 owns them. Resolution direction respected; no script-side additions were required for F1.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy (`scripts/.../stage158_..._sympy_audit.py`): the `Ms = Sigma0 + dSigma0` line and the `expect_zero("delta Ms law", ...)` line are gone. Lines 48-52 now define only `Sigma0, dSigma0, Rstar, dR` symbols, `Mq`, `Mq_lin`, `Mq0`, and the single `delta Mq law` check.
- Mathematica (`mathematica/.../stage158_..._mathematica_audit.wl`): mirror deletion at lines 38-44 — `mS` definition and `expectZero["delta Ms law", ...]` removed; only the `delta Mq law` check remains.

**Assessment:**
Matches directive option (a) verbatim. The transcripts no longer contain `delta Ms law` lines (sympy:5-10, mathematica:5-16), confirming the tautology was removed. No collateral damage to other assertions.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy lines 64-79: new section "3b. Composed boxed identities (notes §3-§4)" with fresh symbols `dg_sym`, `r_sym`, the `dR_from_dg = -dg_sym/sqrt(1+r_sym^2)` carrier, and two `expect_zero` assertions: `composed delta Mq law` and `composed delta Pi law`.
- Mathematica lines 55-65: mirror block with `dgSym`, `rSym`, `dRFromDg`, and two `expectZero` calls for `composed delta Mq law` and `composed delta Pi law`.

**Assessment:**
The new assertions match the directive's symbolic forms verbatim (R_* = 1/4 substituted as `Rational(1,4)` / literal `1/4`; S_* kept symbolic as `Sstar`/`sStar`; Sigma_0 kept symbolic). They are non-tautological: the boxed RHS is constructed independently of `dR_from_dg`, so a sign flip on the dg coefficient would produce a nonzero residual. Both transcripts now show `composed delta Mq law = 0` and `composed delta Pi law = 0` (sympy lines 8-9; mathematica lines 11-14 with `PASS:` confirmations).

## Exec log assessment

**SymPy:** exit=0. All six assertions zero — `linear delta R law`, `delta Mq law`, `delta Pi law`, `composed delta Mq law`, `composed delta Pi law`, `linear Delta_Q law`. No `delta Ms law` line (F2 verified). Numerical coefficients block unchanged.

**Mathematica:** exit=0. All six assertions show `= 0` followed by `PASS:` (lines 5-16). No `delta Ms law` line. Final "Stage 158 Mathematica audit passed." banner present.

**Output freshness:** Confirmed. Sympy script mtime May 28 10:04; output txt mtime May 28 11:30. Mathematica wl mtime May 28 10:04; output txt mtime May 28 11:31. Both outputs newer than their scripts.

## Material-change assessment

`material_change`: false.

No derived numerical result changed. The numerical-coefficients block (sympy lines 96-127; mathematica lines 79-110) is byte-identical in structure and values to the pre-fix run; the printed coefficients match the auditor's report cross-check table exactly. F2 removed a redundant assertion (no value); F3 added new assertions that confirm a symbolic identity already implicit in the printed composed forms. F1 was a paper-side citation rewrite that does not touch any numerics. Downstream consumers see unchanged carry-forward values.

## Side observations (non-blocking)

- Banner still reads "STAGE 158 — LINEAR DEFECT TRANSPORT..." (correct; the auditor noted an old "STAGE 141" residue but the current files already have 158, so either it was fixed earlier or the auditor's note was stale).
- The Mathematica `dPi/dThat` line at txt:28 reports precision `19.69897000433602` rather than `20.` — minor precision-tracking artifact of the multiplication chain, not a correctness issue.

## Verdict justification

All three findings are resolved. F1's paper card was edited per the user-gate direction (Cluster C forward-carry citations) and the orchestrator correctly left scripts untouched for it. F2's tautological `delta Ms law` is gone from both engines. F3's two composed boxed identities are now asserted symbolically in both engines, are non-tautological under the directive's self-test, and both transcripts pass with exit 0. Output mtimes confirm post-fix re-runs. No material change to downstream-consumed numerics. Verdict: `verified`.
