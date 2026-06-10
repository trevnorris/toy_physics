---
unit_id: 224
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 224

## Per-finding outcomes

### F1 — paper_misalignment (card verification-status staleness)

**Classification:** resolved

**What changed:**
Nothing in scripts — correctly. The directive (`redteam/pass2/directives/stage_224.md`) contains ONLY a `paper_misalignment` finding (card `paper/stages/stage_224.tex:11` reads "Mathematica audit: none yet" while a complete passing `.wl` exists). The directive explicitly instructs Codex to do nothing (`applied: false`, `needs_user_resolution: true`) and routes the card-text fix to the user. The diff patch `redteam/pass2/exec_logs/stage_224_diff.patch` is empty; both `.py` and `.wl` are unchanged from HEAD.

**Assessment:**
Correct handling. This is a stale STATUS annotation OUTSIDE scripts-only scope — the card understates the dual-engine coverage that a genuinely-independent, passing `.wl` already provides. It is a paper-side text lag, not a script defect: there is no script change required or possible. Per the standing user decision it is deferred to PAPER_CLEANUP P4-51. Scripts-only verdict is unaffected; this finding alone does NOT roll the unit to needs_rework/blocked. F1 = deferred card-text-lag, resolved at the scripts level (no action correct).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- "Verified exact grouped inverse map: ... Recovered Pbar = (Delta_norm + T_quad)/mhat0**2; Recovered aP = aP; Recovered bP = bP" — inverse map recovers the free symbols (anchored, non-tautological).
- "Exact relation: bP = 3 aP" and the robust collapse "Common robust form = Pbar * (1 + |eps Xi1|) = Pbar + 4 |aP|".
- "All Stage 224 symbolic and numerical checks passed." `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines:
- M1: "solved inverse rule = {pbarUnknown -> ... p20Lane/5 + (2*p21Lane)/5 + (2*p22Lane)/5 ...}" via `Solve` then three "PASS: M1 solved {Pbar,aP,bP} recovers input" (residual = 0) — independently derived inverse, not hardcoded.
- M3: "LinearSolve recovered {Pbar, aP, bP} = {(deltaNorm + tQuad)/mhat0^2, (eps*(deltaNorm + tQuad)*xi1)/(4*mhat0^2), (3*eps*(deltaNorm + tQuad)*xi1)/(4*mhat0^2)}" — distinct inversion mechanism agreeing with SymPy.
- M5: "PASS: M5 positive sign P20 dominates" / "PASS: M5 negative sign P22 dominates" — a domination proof the `.py` does not perform; `Abs` resolved via `Refine`.
- M6: all eight "M6 ... defining relation residual = 0 / PASS" (budgets checked back against the ceiling `pcrit`).
- "All Stage 224 Mathematica checks passed." `# exit_code: 0`.

No FAIL lines in either log.

**Output freshness:** The orchestrator just re-ran both engines directly; both exec logs are dated 2026-06-09T19:21:18-06:00 and are deterministic (the committed `.txt` outputs are byte-identical to the fresh runs per the batch close-out). Not failing on committed `.txt` mtime, per the freshness directive.

## Material-change assessment

`material_change`: false. No script or derived result changed — the diff patch is empty and both engines reproduce the same symbolic forms and numeric budgets as the committed outputs. No downstream unit is affected.

## Side observations (non-blocking)

- Spot-checked the load-bearing asserts for tautology: the M6 headroom checks are genuinely load-bearing — budgets are COMPUTED (`eps_xi_budget = Pcrit_val/barP0_compat - 1`, `a_budget = (Pcrit_val - barP0_compat)/4`, py L145-146; `.wl` L197-198) and then checked back against the **ceiling** via defining relations `(budget+1)*barP0 == Pcrit` and `barP0 + 4*a_budget == Pcrit` (py L155-156; `.wl` L204-210). A wrong rearrangement (e.g. `/2` instead of `/4`) would fail. The first-pass strengthening (budgets → ceiling defining-relations) HOLDS; not value-vs-itself.
- Independence confirmed at the source level: `.py` hardcodes the inverse-projection coefficients (L30-32, L73-74) while `.wl` derives them via `Solve` (L95-100) and `LinearSolve` (L123) on the forward lane matrix `{{1,4,0},{1,-1,1},{1,-1,-1}}`. Agreement between the two routes is cross-engine confirmation, not a shared derivation. The only shared items are physical premises (lane matrix and signature `{1,1/2,-1}`). Independence call: independent.

## Verdict justification

The single directive finding (F1) is a paper-card verification-status lag ("Mathematica audit: none yet" while a passing independent `.wl` exists) that is outside scripts-only scope and correctly deferred to PAPER_CLEANUP P4-51 by standing user decision; no script change is required or possible, so Codex's no-op (empty diff) is correct. On the scripts-only mandate the unit is clean: both engines re-ran fresh to exit 0 with all-PASS bodies and no FAIL lines, the load-bearing M6 budget asserts are non-tautological (computed budgets checked against the ceiling), and the `.wl` is a genuine independent recomputation (Solve/LinearSolve-derived inverse map and a `.wl`-only domination proof) rather than a transliteration. Verdict: `verified`.
