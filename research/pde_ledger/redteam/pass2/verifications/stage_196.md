---
unit_id: 196
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 196

Both findings in the original report are non-script (F1 = stale committed SymPy `.txt`, refreshed by the orchestrator's independent re-run; F2 = paper-side card-text lag, USER-DEFERRED to paper-cleanup). No Codex source edit was directed and none was needed. Verification therefore confirms (A) the output refresh landed clean, (B) the audit disposition still holds on the refreshed artifacts, and (C) the deferral is correctly classed.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 179" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run regenerated the committed outputs. `scripts/output/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.txt` now reads `STAGE 196 — EXACT HIGHER-ODD IRRELEVANCE THEOREM` at line 3 and `STAGE 196 LEDGER` at line 165 (the stale "STAGE 179" banner is gone). Its mtime is Jun 9 14:07, newer than the `.py` (Jun 3 15:59). The Mathematica output was already current and was likewise refreshed (mtime Jun 9 14:07), banner `STAGE 196 - HIGHER-ODD IRRELEVANCE THEOREM` (line 3) and `STAGE 196 MATHEMATICA LEDGER` (line 72).

**Assessment:**
Correct and complete. The refresh is the prescribed remedy. Every SymPy result line is an equality `= 0`: §I lines 30, 48, 49, 50; §II lines 65, 66, 79; §III lines 131, 132, 133; §IV lines 151, 152, 153; §V line 162. The Mathematica output prints an explicit `PASS:` for every check (lines 13, 15, 17, 19, 26, 28, 30, 32, 43, 45, 47, 49, 51, 53, 55, 63, 65, 67, 69) with no `FAIL` anywhere. Output refresh is clean.

### F2 — paper_misalignment (card-text lag: "Mathematica audit: none yet")

**Classification:** resolved (correctly routed — non-script, USER-DEFERRED)

**What changed:**
Nothing in scripts (paper.tex is off-limits to the red-team). The finding routes `stage_196.tex:11` to the paper-cleanup tracker to cite the present passing `.wl`. This is correctly classified as a paper-prose lag with no math impact.

**Assessment:**
Legitimately deferred. The `.wl` exists and passes (Mathematica output lines 13-69 all PASS), so the card text understates coverage but contradicts no result. Non-blocking; outside the scripts-only scope.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed on the refreshed artifacts. §III regenerates the canonical-even slots from the native outgoing operator — `lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2,x],x]/SphericalHankelH1[2,x]]` (wl:134) and `Solve` for the even slots (wl:146) — a route the `.py` does not contain (grep for `SphericalHankel`/`FunctionExpand` in the `.py` returns nothing). The native-Hankel window prints in the refreshed Mathematica output: "outgoing l=2 DtN window through z^5 = -3 + x^2/3 + x^4/9 + (I/9)*x^5" (line 37), and the Solve-derived matching laws reconcile to the SymPy targets (lines 42-51). Both engines land on the same `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)` (py out L121-130 / wl out L41).
- **`l7` derivative non-vacuous:** confirmed. `l7` genuinely appears at the z⁷ slot — SymPy output line 96 shows `729·i·L7·z^7/(2187S−729Σ₀)` in the matched normalized DtN response, and Mathematica output line 24 shows `- (I*l7*z^7)/l0` in the generic DtN response. So `d chi_Q/dL7 = 0` (py out L133 / wl out L52) and the z^0..z^5 independence checks (wl out L54-55 `{0,0,0,0,0,0}`) are real can-fail statements, not trivially-true.
- **0 reconciliation misalignments:** confirmed. The report's deliverable table reconciles 13/13 values MATCH against the notes; every reconciled value (`sigma_Q^can`, `Sigma_2`, `Sigma_4`, `chi_Q`, `N_Q`, `P0_target`, `Delta_norm`, L7 first-entry coeff `-i/L0`, `Delta_Q`) is present and identical in the refreshed outputs (e.g. `N_Q` py out L138-144 / wl out L60; `Delta_norm` py out L145-153 / wl out L61; `-I/L0` py out L159 / wl out L31).

## Exec log assessment

**SymPy:** exit=0 (inferred from refreshed committed output). Notable lines: "exact response-side difference identity = 0" (L30); "first higher-odd correction really starts at omega^7 = 0" (L50); "chi_Q extractor - deformation algebra formula = 0" (L132); "d chi_Q^(series extractor) / dL7 = 0" (L133); "Delta_norm - P0_target*(1/chi_Q - 1) = 0" (L153). All equalities hold to zero.

**Mathematica:** exit=0 (inferred from refreshed committed output). Notable lines: "PASS: exact response-side difference identity - SymPy target" (L13); "PASS: matched coefficients z^0..z^5 independent of L7" (L55); "PASS: chi_Q extractor - SymPy deformation formula" (L51); "PASS: Delta_norm - P0_target*(1/chi_Q - 1)" (L69). Every check prints PASS; no FAIL.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime Jun 9 14:07, newer than `.py` (Jun 3 15:59) and `.wl` (Jun 1 12:00). SymPy banner is now canonical "STAGE 196" (the F1 stale-179 banner is cleared).

## Material-change assessment

`material_change`: false. No source code changed; the only edits were regenerated committed `.txt` outputs (banner relabel + transcript refresh). No derived result changed, so no downstream unit is affected.

## Side observations (non-blocking)

None beyond the two reported findings. The §I/§II "forced-form" identities were re-examined by the auditor and correctly not flagged as transliteration (a difference of two rational functions is algebraically forced to its unique closed form); the load-bearing independence lives in §III and is real. I concur and add nothing.

## Verdict justification

`verified`. Both findings are non-script and resolved: F1's stale SymPy banner is cleared by the orchestrator re-run — the refreshed `.txt` carries the canonical "STAGE 196" banner with every SymPy check `= 0`, and the Mathematica output prints PASS on every check with no FAIL; F2 is a paper-side card-text lag correctly USER-DEFERRED to paper-cleanup. The audit disposition holds on the refreshed artifacts: the `.wl` is genuinely independent (native `SphericalHankelH1`/`Solve` §III route absent from the `.py`), the `l7` derivative is non-vacuous (`l7` present at the z⁷ slot in both outputs), and reconciliation is 13/13 MATCH with 0 misalignments. No regressions; `material_change: false`.
