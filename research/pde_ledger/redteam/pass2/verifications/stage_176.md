---
unit_id: 176
batch: V.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 176

This stage was audited **clean** in pass-2 (0 findings, no Codex changes). Scripts were
not modified; the orchestrator's reliability re-run refreshed the committed outputs. This is
a LIGHT confirmation per the verify prompt: (a) audit verdict + independent-route, (b)
deliverable values still verify, (c) banner/exit/no-FAIL, (d) empty script diff.

## Per-finding outcomes

No findings in the original report (`findings_count: 0`). Nothing to resolve. The
confirmation checks below stand in for the per-finding pass.

### (a) Audit verdict holds — `.wl` is a genuinely independent route

**Classification:** resolved (verdict confirmed clean).

The Mathematica audit is an independent derivation, not a port of the SymPy script. The
load-bearing object `Sigma_exact` — the first-order log-drift of `Lambda^2/K` that A2/A3/A4
and B2/B3/B4 all test against the claimed forms — is extracted by **different machinery in
each engine**:

- SymPy: truncated Taylor series, reads the `eps^1` coefficient —
  `sp.series(sp.log((Lambdap**2/Kp)/(Lambda**2/K)), eps, 0, 2).removeO()/eps` (py:62-72).
- Mathematica: symbolic derivative evaluated at the origin —
  `D[Log[((lambdaP^2/kP)/(lambda^2/k))], eps] /. eps -> 0` (wl:62-65),
  explicitly documented as divergent in the `.wl` comment at wl:60-61.

The primitives (`lambda`, `mCal`, `iCal`, `hCal`) are each rebuilt from the same six port
symbols in Mathematica (wl:32-35), not imported from a pre-simplified SymPy form;
simplifier is `FullSimplify[Together[Expand[...]]]` (wl:20-24) vs SymPy `simplify(expand(...))`.
The candidate forms `Sigma_factored`/`Sigma_expanded`/`2·dlnM` are necessarily identical in
both engines because they *are* the paper's claimed results, but the independence lives in
`Sigma_exact`. Independent route confirmed; no transliteration concern.

### (b) Deliverable values still verify in the refreshed outputs

**Classification:** resolved.

All four residuals are `0` in both refreshed transcripts:
- factorization of `Lambda^2/K` = 0 (sympy out L10; mathematica out L10-11 PASS)
- factored first-order defect = 0 (sympy L18; mathematica L19-20 PASS)
- expanded primitive-variable transport = 0 (sympy L19; mathematica L21-22 PASS)
- rigidity reduction to `2 d ln M` = 0 (sympy L30; mathematica L33-34 PASS)

### (c) Banner / exit / no FAIL

**Classification:** resolved.

Both exec logs carry banner `STAGE 176 — OUTGOING LOAD-FACTOR FACTORIZATION`. SymPy
`exit_code: 0`; Mathematica `exit_code: 0` with `Stage 176 Mathematica audit passed.`
(out L42). No `FAIL:` line in either transcript; all four Mathematica checks print `PASS:`.

### (d) Empty script diff

**Classification:** resolved.

`redteam/pass2/exec_logs/stage_176_diff.patch` is 0 lines. `git diff --stat HEAD` against the
`.py` and `.wl` returns nothing. The `.py` (2026-05-30 01:10:19) and `.wl` (2026-05-30
01:10:19) are untouched since pass-1; only outputs were regenerated. Scripts unchanged.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `exact factorization of Lambda^2/K = 0` (L10)
- `factored first-order defect formula = 0` (L18)
- `expanded primitive-variable transport = 0` (L19)
- `rigidity reduction of Sigma_exact to 2 d ln M = 0` (L30)

**Mathematica:** exit=0. Notable lines:
- `exact factorization of Lambda^2/K = 0` / `PASS:` (L10-11)
- `factored first-order defect formula = 0` / `PASS:` (L19-20)
- `expanded primitive-variable transport = 0` / `PASS:` (L21-22)
- `rigidity reduction of Sigma_exact to 2 d ln M = 0` / `PASS:` (L33-34); `Stage 176 Mathematica audit passed.` (L42)

**Output freshness:** confirmed fresh. `.py` 2026-05-30 01:10:19 → `.txt` 2026-06-09
00:19:04; `.wl` 2026-05-30 01:10:19 → `.txt` 2026-06-09 00:19:04. Both transcripts post-date
their scripts by the orchestrator's reliability re-run. No staleness.

## Material-change assessment

`material_change`: false. No script edit; no derived result changed. Nothing downstream is
affected by this verification.

## Side observations (non-blocking)

The directive file `redteam/pass2/directives/stage_176.md` does not exist — expected for a
0-finding clean stage with no Codex pass. No issue.

## Verdict justification

`verified`. The stage was clean in pass-2 with zero findings; this light confirmation holds
on all four axes. The `.wl` is a genuinely independent route — `Sigma_exact` is obtained by
symbolic `D[Log,eps]` at `eps->0` in Mathematica vs a Taylor-series `eps^1` read in SymPy,
with primitives rebuilt from premises in both. All four deliverable residuals are `0` and
every Mathematica check prints `PASS:`; both engines exit 0 under the `STAGE 176` banner with
no `FAIL`. The committed outputs are fresh (mtimes 2026-06-09, newer than the 2026-05-30
scripts), and the script diff is empty.
