---
unit_id: 215
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 215

The single finding in the original report is non-script (F1 = paper-side card-text lag: the `\stagefield{Verification}` line still reads "Mathematica audit: none yet" though a passing dual-engine `.wl` was added in the pass-1 retrofit). No Codex source edit was directed (`applied: false`, no `## Applied:` block) and none was needed — the directive holds before invoking Codex and routes the card edit to user resolution (P4-51). Verification therefore confirms (A) the committed outputs are clean and fresh, (B) the audit disposition still holds on those artifacts (independent `.wl`, non-tautological checks), (C) the card-text-lag deferral to P4-51 is correctly classed, and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (card-text lag: "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to P4-51)

**What changed:**
Nothing in scripts or the card. paper.tex is off-limits to the red-team, and per file-ownership the card edit is applied by Codex only after the user authorizes a follow-up directive (Claude reviews). The directive (`directives/stage_215.md`) carries only the `## Resolve before fix_loop` (P4-51 deferral) block — `applied: false`, no Codex edit, no diff patch produced (`exec_logs/stage_215_diff.patch` absent, as expected for a no-edit unit). The finding routes `stage_215.tex:11` to the paper-cleanup gate to cite the present passing `.wl` and drop "none yet".

**Assessment:**
Legitimately deferred and correctly classed. This is a pure prose/status discrepancy, not a math error: the `.wl` (M1-M7, 256 lines) exists and passes — committed output ends "Stage 215 Mathematica audit passed." with every `PASS:` present and no `FAIL`. The card text understates verification coverage (claims single-engine where the unit is now dual-engine) but contradicts no result. Non-blocking; outside the scripts-only scope. Deferral to P4-51 is the prescribed remedy.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed. The load-bearing interval/order theorems (claims 2-6) are proved by two methodologically distinct routes — the `.py` enumerates ordered integer tuples and checks the inequality pointwise (e.g. the full-simplex splice over 3136 samples, the unique-winner over 542760 min-over-others samples), while the `.wl` decides the universally-quantified ordering over the continuum via `Resolve[ForAll[{reals}, Implies[premise, conclusion]], Reals]` real quantifier elimination (M3-M6). Enumerate-and-check vs. decide-over-reals is a genuinely different extraction on every load-bearing object — INDEPENDENT, not a port. The shared combinatorial counting (M1 vs py I) and the one-line `Min`-flattening associativity (M2 vs py II) are shared *definitional* premises, which is permitted.
- **Non-tautological checks:** confirmed. A4-A7 / M3-M6 test order relations of independently-defined `Min` expressions that could fail for a wrong definition; A1-A2 / M1 count an independently constructed combinatorial object (10 triples / 5 quadruples / 4 faces / 2 incidence / 4 axis); A8 / M7 factorize the literal budgets (`3*3*3*2=54`, `10*12+10*48=600`, `5*2*54=540`, `600+540=1140`) rather than asserting a number against itself.
- **Engines agree:** confirmed. Both committed outputs land on identical numeric deliverables (per-envelope 54, support-≤3 600, interior 540, full 1140; combinatorics 10/5/4/2/4) and finish with their pass banners.

## Exec log assessment

**SymPy:** exit=0. Notable lines: "nested Min flattening (lo) = True" / "(hi) = True"; "verified local full-simplex interval theorem on 3136 ordered integer samples"; "verified unique certified winner theorem on 542760 min-over-others samples"; "full support<=4 budget - 1140 = 0"; closes "All Stage 215 identities and interval theorems verified."

**Mathematica:** exit=0. Notable lines: "M5 unique-winner quantified results = {True, True, True, True, True}"; "PASS: M5 unique certified winner over five quadruple intervals"; "PASS: M6 support<=4 splice theorem"; "M7 per-envelope candidate bound = 54" → "PASS: M7 full support<=4 budget - 1140"; closes "Stage 215 Mathematica audit passed." Every check prints PASS; no FAIL.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime 2026-06-09 16:51, newer than the `.py` and `.wl` (both 2026-06-02 11:0x). Banners are canonical "STAGE 215" (no stale renumber label), `git status --porcelain` is empty for all four paths (refreshed outputs were already byte-identical / committed), and no `FAIL` appears in either committed output.

## Material-change assessment

`material_change`: false. No source code changed — the only orchestrator action was the independent re-run that confirmed the committed outputs were already fresh (byte-identical). No derived result changed, so no downstream unit is affected. The pending card-text edit (P4-51) is paper prose only and carries no math impact.

## Side observations (non-blocking)

None beyond the single reported finding. The combinatorial block (M1 vs py I) and the `Min`-flattening identity (M2 vs py II) share a definitional premise, but the load-bearing extraction (the certified-interval orderings) is genuinely independent across engines; I concur with the auditor's INDEPENDENT verdict and add nothing.

## Verdict justification

`verified`. The lone finding is non-script and resolved: F1 is a paper-side card-text lag ("Mathematica audit: none yet" vs. a passing independent `.wl`), correctly classed as a `paper_misalignment` outside the scripts-only scope and USER-DEFERRED to P4-51 — no Codex edit was directed or needed (`applied: false`, no diff). The audit disposition holds on the committed artifacts: both engines exit 0 with every check passing and no FAIL, the `.wl` is genuinely independent (`Resolve[ForAll, …, Reals]` real-QE vs. `.py` exhaustive integer enumeration on every load-bearing interval theorem), and the assertions are non-tautological. Outputs are clean and fresh (mtime 2026-06-09 16:51 > scripts 2026-06-02; canonical STAGE 215 banners; clean git tree). No regressions; `material_change: false`.
