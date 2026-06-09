---
unit_id: 202
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

# Verification — unit 202

The single finding is a `paper_misalignment` (F1 = card/notes say "Mathematica audit: none yet" though a complete passing `.wl` exists). It is the card-text-lag class, USER-DEFERRED to PAPER_CLEANUP P4-51 — no Codex script edit was directed and none was needed (directive `applied: false`, `needs_user_resolution: true`; diff patch absent, as expected). A separate stale committed SymPy `.txt` banner (was "STAGE 185") was refreshed by the orchestrator's independent re-run. Verification therefore confirms (A) the output refresh landed clean, (B) the audit disposition still holds on the refreshed artifacts, (C) the deferral to P4-51 is correctly classed (stale STATUS annotation, not a value/identity mismatch), and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (card-text lag: "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to PAPER_CLEANUP P4-51)

**What changed:**
Nothing in scripts. `paper/` and `notes/` are off-limits to the red-team and the directive holds for user resolution (the `## Resolve before fix_loop` block, no `## Applied:`). There is no `stage_202_diff.patch` and no source-file change — consistent with a finding that prescribes no Codex edit. The card line `stage_202.tex:11` ("Mathematica audit: none yet") and the notes Supporting-files omission of the `.wl` are routed to PAPER_CLEANUP P4-51 to be reconciled against the present passing `.wl`.

**Assessment:**
Correctly classed and legitimately deferred. This is a stale STATUS annotation: the prose understates coverage (claims single-engine) while a full dual-engine audit exists and passes. It is NOT a value or identity mismatch — every symbolic deliverable reconciles 11/11 MATCH per the report's value table, and both engines reach all-zero residuals on the shared identities. The `.wl` is genuinely independent of the `.py` (the `.py` POSITS the boxed closed forms `deltaU_graph/T_graph/Keta_graph/mu_graph` at py:92–104 then verifies by back-substitution; the `.wl` DERIVES them via `LinearSolve` of the 3x3 log-linear dependent-triple system at wl:112–135, an output of the solve — derive-vs-posit on the load-bearing object). The first-pass iter-2 rework (transcendental `Solve` → log-linearized `LinearSolve`) lives ONLY in the `.wl`; the `.py` never solves, so there is no shared solve to be a port of. Non-blocking and outside scripts-only scope.

## Disposition re-confirmation (post-refresh)

- **Output refresh clean (A):** the committed SymPy `.txt` now carries the canonical banner `STAGE 202 — EXACT FREE-QUINTUPLE TARGET GRAPH AND THE FIRST REDUCED-FAMILY TEST` (line 3) and `STAGE 202 AUDIT COMPLETE` (line 620) — the stale "STAGE 185" banner is gone. Its mtime is Jun 9 16:51, newer than the `.py` (Jun 3 15:59) and `.wl` (Jun 2 09:52). The Mathematica `.txt` banner reads `STAGE 202 - FREE-QUINTUPLE TARGET GRAPH MATHEMATICA AUDIT` (line 3) / `STAGE 202 MATHEMATICA AUDIT PASSED` (line 63), mtime Jun 9 16:51.
- **Audit disposition holds (B):** SymPy committed output shows every result line collapsing to zero — `log(target Ctr/Cnt/eps reconstructed by graph)=0` (out L84–86), `L in mu_graph.free_symbols? no` (L87), the 8-component `graph projection equals canonical orbit projection = [0…0]` (L250–265), the three graph-error identities + `E-log(m)` rewrites all `=0` (L317–322), `repair vector rewrite = [0…0]` (L476–490), `family packet on target graph = [0,0,0,0]` (L552–558), and `multiplicative family chart on target graph = [0,0,0,0]` (L616–622). Mathematica committed output prints PASS on all eleven checks (graph log-system residual M0, M1–M3 monomials, M4 `mu_graph` free of L and Pi + `d log(mu_graph)/dL=0`, M5 quotient packet, M6 8-comp projection + threaded equations, M7 repair + reduced-family packet) with no FAIL.
- **Deferral correctly classed (C):** F1 is a stale `\stagefield{Verification}` annotation reconciled to P4-51, not a value/identity mismatch — there is no numeric or symbolic disagreement anywhere; the report's reconciliation table is 11/11 MATCH, 0 misaligned.

## Exec log assessment

**SymPy:** exit=0 (exec log header `# argv: timeout 600 python3 …stage202…sympy_audit.py`, footer `# exit_code: 0`). Notable lines: "log(target Ctr reconstructed by graph / Ctr_target) = 0" (L84); "L in mu_graph.free_symbols? no" (L87); "E_T - q_tr/(1+chi0_*) = 0", "E_K + q_eta = 0", "E_mu - (q_nt - q_eta + F_* q_tr/(1+chi0_*)) = 0" (L317–319); graph-projection / repair / family vectors all `[0]` (L250–265, L476–490, L552–558). No FAIL/traceback.

**Mathematica:** exit=0 (footer `# exit_code: 0`; the re-run took 14s per the orchestrator note — well under the 600s timeout, confirming the iter-2 `LinearSolve` rework is deterministic). Notable lines: "graph log-system residual = {{0,0,0}}" + "PASS: graph log-system residual" (L18–20); "PASS: M1/M2/M3 … log residual" (L25–29); "M4 mu_graph is free of L and Pi = True" + "PASS: M4 d log(mu_graph)/dL" (L31–34); "M5 graph-error quotient packet = {{0,0,0}}" PASS (L39–41); "M6 projection component log residuals = {{0×8}}" + threaded equations True PASS (L46–50); "M7 repair-vector rewrite residuals = {{0×8}}" + "M7 reduced-family packet … = {{0,0,0,0}}" PASS (L55–60). Every check PASS; zero FAIL.

**Output freshness:** confirmed. Both committed `.txt` carry mtime Jun 9 16:51, newer than the `.py` (Jun 3 15:59) and `.wl` (Jun 2 09:52). SymPy banner is now canonical "STAGE 202" (the stale-185 banner is cleared). No `stage_202_diff.patch` exists — correct, since no source file was edited.

## Material-change assessment

`material_change`: false. No source code changed (no diff patch, no `## Applied:` block; F1 is a user-deferred paper-prose lag). The only edits were the regenerated committed `.txt` outputs (banner relabel + transcript refresh). No derived result changed, so no downstream unit (> 202) is affected.

## Side observations (non-blocking)

None beyond the reported finding. I concur with the auditor that the substitution/identity checks are non-tautological (a posited/derived graph is fed into independently-written monomials, so a wrong graph formula would yield a nonzero residual) and that the `(1+δ_{U,*})` exponent is a carried constant `deltaUs`, not the solved `δ_U` — no circularity. The `.wl` carries two genuinely-extra coverage items the `.py` lacks (the M0 `LinearSolve` consistency residual and the explicit `D[logMuGraph,L]=0` derivative on top of `FreeQ`), strengthening rather than echoing.

## Verdict justification

`verified`. The single finding is a non-script paper_misalignment correctly classed as the card-text-lag class and USER-DEFERRED to PAPER_CLEANUP P4-51 (stale "Mathematica audit: none yet" STATUS annotation, not a value/identity mismatch); no Codex edit was directed and none exists (diff patch absent). The orchestrator's independent re-run cleared the stale SymPy "STAGE 185" banner — the refreshed `.txt` now carries the canonical "STAGE 202" banner with every SymPy result `= 0`, and the Mathematica output prints PASS on all eleven checks (M0/M1–M7) with no FAIL, both exit 0. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (LinearSolve-derives vs `.py`-posits the load-bearing graph object; the iter-2 linearization lives only in the `.wl`). No regressions; `material_change: false`.
