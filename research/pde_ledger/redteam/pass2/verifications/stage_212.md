---
unit_id: 212
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

# Verification — unit 212

The sole finding (F1) is a paper-side stale verification-status line ("Mathematica audit: none yet" on the card while a passing `.wl` exists). It is non-script and USER-DEFERRED to P4-51; the directive explicitly holds Codex (`applied: false`, `needs_user_resolution: true`) and prescribes no script edit. Verification therefore confirms (A) outputs clean/fresh, (B) the audit disposition holds on the refreshed artifacts, (C) the card-text-lag deferral is correctly classed and the pass-1 notes-typo corrections still hold, and (D) `material_change: false`.

## Per-finding outcomes

### F1 — paper_misalignment (stale verification-status line: card "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED to P4-51)

**What changed:**
Nothing in scripts. There is no `stage_212_diff.patch` (confirmed absent), consistent with the directive holding Codex: `applied: false`, `needs_user_resolution: true`, "Codex applies NOTHING here." The card-text edit (`paper/stages/stage_212.tex:11` cite the present `.wl`, drop "none yet") is the user's call, deferred to paper-cleanup P4-51. paper.tex is off-limits to the red-team.

**Assessment:**
Legitimately deferred and correctly classed as a paper-prose status lag with no math impact. A full Mathematica audit (M1-M4) is present and passes — the refreshed transcript ends "Stage 212 Mathematica audit passed." with every M1-M4 check PASS and no FAIL — so the card text understates coverage but contradicts no result. Non-blocking; outside the scripts-only scope. The directive's expected resolution (direction (a): card stale) matches the deferral.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed. The load-bearing objects are the interval/order theorems. SymPy proves block-III by finite integer-grid enumeration (3136/924/462/896 ordered samples; out L54-60), which alone samples a finite subset. Mathematica proves the same theorems by symbolic quantifier elimination over the reals — `Resolve[ForAll[..., Implies[...], ...], Reals]` (M3a-M3d, out L39-52 all `True`) — a real-closed-field decision procedure that decides the universal statement for ALL reals. Categorically different proof strategies (sample-the-grid vs. decide-over-Reals); not a port. The shared M1 combinatorics (`combinations`/`Subsets`), M2 Min auto-flattening, and M4 arithmetic are non-load-bearing scaffolding/definitions where share-of-premise is allowed.
- **Non-tautological assertions:** confirmed. Counts are computed one route (`combinations`/`Subsets`) and checked against an independent route (`binomial`/`Binomial`) — `#pairs - binomial(5,2) = 0`, `#triples - binomial(5,3) = 0` (sympy out L19-20); the interval theorems are falsifiable quantified inequalities; the budget literals are checked against the product decomposition `10*12=120`, `10*48=480`, `120+480=600` (wl out L60/L57/L64), not against themselves.
- **Notes-typo corrections HOLD:** confirmed. Grep of both committed outputs for the stale tokens `188` / `243` / `245` / `246` / `247` / stale `STAGE 24x` / `STAGE 18x` returns NONE. The renumber to 209/211/212/213 is present (the script self-labels "STAGE 212"), and budget `120/480/600` is present and correct in both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines: "#pairs - binomial(5,2) = 0" (L19), "#triples - binomial(5,3) = 0" (L20); "nested Min flattening (lo) = True / (hi) = True" (L48-49); "verified local full-simplex interval theorem on 3136 ordered integer samples" (L54) through "global three-coordinate no-improvement theorem on 462 samples" (L60); "pairwise/triple interior/full budget - 120/-480/-600 = 0" (L68-70); "All Stage 212 identities and interval theorems verified." (L72).

**Mathematica:** exit=0. Notable lines: "PASS: M1 #pairs" (L20), "PASS: M1 #triples" (L22); "PASS: M1 every primitive pair lies in three triples" (L24); "PASS: M2 lo/hi nested Min identity" (L32/L34); M3a-M3d all `True` with "PASS" (L39-52); "PASS: M4 10*12 / 10*48 / 120+480" (L61/L63/L65); "Stage 212 Mathematica audit passed." (L67). Every check PASS; no FAIL.

**Output freshness:** confirmed. Both committed `.txt` outputs carry mtime 2026-06-09T16:51:54, newer than the `.py`/`.wl` sources (both 2026-06-02T11:38:51). Committed `.txt` tails match the exec-log transcripts (budget 120/480/600 and the two "passed" banners). The orchestrator's independent re-run refreshed them post-fix (212's outputs were already fresh — byte-identical refresh).

## Material-change assessment

`material_change`: false. No source code changed (no diff patch; Codex held for user resolution); the only artifact touch was the orchestrator's transcript refresh (byte-identical). No derived result changed, so no downstream unit is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The auditor's INDEPENDENT verdict and the non-tautology checks re-confirm cleanly on the refreshed artifacts.

## Verdict justification

`verified`. The lone finding F1 is a paper-side stale card status line ("Mathematica audit: none yet" while a passing M1-M4 `.wl` exists), correctly classed as a documentation-status lag with no math impact and USER-DEFERRED to paper-cleanup P4-51; the directive held Codex (no diff patch, `applied: false`), which is the prescribed handling. Outputs are clean and fresh (both `.txt` mtimes 16:51:54 newer than the 11:38:51 sources; budget 120/480/600 present in both engines; no stale 188/243/245/246/247 tokens). The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (`Resolve[ForAll,...,Reals]` symbolic QE vs `.py` finite integer-grid enumeration), assertions are non-tautological (counts vs `binomial`, budget vs product decomposition), and the pass-1 notes-typo corrections (188→120, renumber to 209/211/212/213) still hold. Both engines exit 0 with every check PASS and no FAIL. `material_change: false`.
