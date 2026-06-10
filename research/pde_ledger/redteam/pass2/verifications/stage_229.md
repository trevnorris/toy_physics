---
unit_id: 229
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 229

## Per-finding outcomes

The auditor report for unit 229 carries `verdict: clean` with `findings_count: 0`,
and accordingly NO directive file was written (`redteam/pass2/directives/stage_229.md`
does not exist — confirmed by `ls`, expected). There are therefore no findings to
classify. The verifier's job here reduces to confirming that the `clean` verdict still
holds: re-reading both scripts, confirming both exec logs pass, and spot-checking that
the load-bearing assertions are non-tautological and that the `.wl` is a genuine
independent recomputation rather than a transliteration.

### (no findings) — confirmation that `clean` holds

**Classification:** resolved (vacuously — nothing to fix; `clean` re-confirmed)

**What changed:** Nothing. The git diff patch
(`exec_logs/stage_229_diff.patch`) is 0 bytes — legitimately empty, not a
truncation. `git status --porcelain` on both script files returns empty: the
`.py` and `.wl` are byte-identical to HEAD. This is the expected signature of the
first pass-2 batch (VII.1) that required zero script corrections; Codex applied no
edits because the auditor found nothing to fix.

**Assessment:** Confirmed clean. Detail below.

1. **Scripts re-read.** Both the SymPy
   (`scripts/...stage229...sympy_audit.py`) and Mathematica
   (`mathematica/...stage229...mathematica_audit.wl`) scripts were read in full.
   They match the auditor's description.

2. **Both exec logs pass.** SymPy log ends `# exit_code: 0` (line 31); body emits
   "All Stage 229 symbolic and numerical audits passed." Mathematica log ends
   `# exit_code: 0` (line 77); body shows 24 `PASS` lines spanning M1–M10
   (M1 reduction, M2 factorization, M3 log-slopes, M4 classifier, M5 onset,
   M6 softening limits + pole reciprocal, M7 four cubic coefficients,
   M8 three derivative coefficients + positivity, M9 P(0,δ) + threshold,
   M10 three sample roots + side classifications + δ=1 slice) and "Stage 229
   Mathematica audit passed." No FAIL lines, no warnings.

3. **Load-bearing asserts are non-tautological.** Each `.py` assertion compares a
   closed form against an *independently derived* expression, not against itself:
   - py:49 — `N_-` substitution residual vs `(8β₀/π²A)·F`, where `F` is built
     from the dimensionless ansatz, then `simplify(...) == 0`.
   - py:53 — `F − F_num·F_den == 0` (factorization).
   - py:71-73 — `L_num`, `L_den`, `R_ND` obtained by *symbolic differentiation*
     of `log F_num`/`log F_den` (py:63-65), compared to the boxed closed forms.
   - py:97 — `numerator = expand(-together(R_ND − 1).as_numer_denom()[0])`
     (the crossover cubic is *derived from* `R_ND − 1`, not pinned), asserted
     equal to `expected_P` whose leading term is `121*xi**3`.
   - py:104 — `P(0,δ) − 9δ²(9δ−8) == 0`.
   - py:124-134 — `select_stable_real_root` requires a *unique* real root in
     (0,1) (matching the strict-monotonicity theorem), and the side-probe checks
     require `R_ND > 1` left of the root and `< 1` right of it (genuine
     classification, would fail if the signature were wrong).
   None of these is self-referential.

4. **`.wl` is a genuine independent recomputation, not a transliteration.** The
   Mathematica path uses a different simplification stack (`Cancel∘Together` +
   `FullSimplify`, wl:18-23, 72-74) versus SymPy's `simplify`/`together`, and
   re-derives `F` from the constants (wl:74, M1 "dimensionless F from constants").
   Decisively, the `.wl` adds a universal-quantifier *proof* the `.py` lacks:
   `Resolve[ForAll[{xiPos,deltaPos}, Implies[xiPos>=0 && deltaPos>0,
   (dP/.…)>0]], Reals]` (wl:170-180), output "M8 derivative positivity = True".
   The `.py` only checks the *coefficient form* of `∂_ξP` (py:99-101) and never
   proves positivity. A CAD-based `Resolve[ForAll]` decision has no SymPy
   counterpart — the second engine does strictly *more* than the first, the
   opposite of a transliteration. Root-finding also differs (`sp.nroots` vs
   `NSolve[…, Reals]`).

5. **Crossover-cubic leading coefficient `121ξ³` emitted; no surviving `189ξ³`.**
   - SymPy: py:96 `expected_P = 121*xi**3 + …`, derived-and-asserted at py:95-97;
     output line 20 emits `… + 121*xi**3`.
   - Mathematica: `targetCoeffs = {…, 121}` (wl:145); `P` is derived from
     `Numerator[Together[RND − 1]]` (wl:137), `CoefficientList` taken, and the
     `ξ³` slot matched — log lines 35 (`… 121*xi^3`) and 42-43
     (M7 coefficient 3 residual 0, PASS).
   - No `189` appears in either script or either exec log. The −68 correction is
     fully landed on both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `F(xi,delta) = (9*delta + 11*xi)**4/((81 - 81*xi)*(9*delta**2 + 18*delta*xi + 11*xi**2)**2)` (line 10)
- `P(xi,delta)   = 81*delta**3 + 333*delta**2*xi - 72*delta**2 + 297*delta*xi**2 + 121*xi**3` (line 20 — leading `121*xi**3`)
- `R_ND(0,delta) = 8/(9*delta)`, `lim_{xi->1^-} R_ND = 0` (lines 16-17)
- `All Stage 229 symbolic and numerical audits passed.` (line 30)

**Mathematica:** exit=0. Notable lines:
- `M7 coefficient list = {-72*delta^2 + 81*delta^3, 333*delta^2, 297*delta, 121}` (line 35) + `PASS M7 coefficient 3` (line 42-43)
- `M8 derivative positivity = True` / `PASS M8 derivative positivity` (lines 51-52) — the `Resolve[ForAll]` universal proof
- `M9 P(0,delta) = 9*delta^2*(-8 + 9*delta)` / `PASS M9 onset polynomial factor` + `PASS M9 threshold root` (lines 53-57)
- Sample roots match SymPy to ~1e-16 (lines 58-71); `M10 delta=1 always denominator slice = True` (line 74)
- `Stage 229 Mathematica audit passed.` (line 76)

Both engines agree on every load-bearing emitted value (F, R_ND, onset, softening
limit, crossover cubic `121` coefficient, P(0,δ), all three sample roots, δ=1 slice).

**Output freshness:** Per the batch context, the orchestrator re-ran both engines
directly; the exec logs are dated 2026-06-09T19:21:18 and exit 0. The run is
deterministic — the committed `.txt` outputs are byte-identical to the fresh run.
Empty diff is consistent with no script change. Not failing on committed-`.txt` mtime,
per the freshness note.

## Material-change assessment

`material_change`: false. No script was edited (empty diff, scripts byte-identical to
HEAD), so no derived result changed. Downstream units (>229) are unaffected by this
verification; nothing to flag.

## Side observations (non-blocking)

The auditor noted (non-blocking) that the stage card's `\stagefield{Verification}`
line still reads "Mathematica audit: none yet" while a substantive `.wl` now exists —
a stale prose line in the card (doc-pointer-sync class), not a math defect and outside
the scripts-only scope. Recorded here only to carry it forward to the doc sync; it does
NOT affect this verdict.

## Verdict justification

`verified`. The auditor's `clean`/0-findings verdict holds. There were no findings and
no directive, the git diff is legitimately empty, and both scripts are byte-identical to
HEAD. Both exec logs exit 0 with full PASS coverage (M1–M10 on the `.wl`; the single
"all passed" line plus emitted values on the `.py`). The load-bearing assertions are
non-tautological (each closed form is checked against an independently derived
expression — substitution+simplify, symbolic differentiation, or the numerator of
`R_ND − 1`), and the `.wl` is a genuine independent recomputation that does strictly
more than the `.py` via the `Resolve[ForAll]` universal positivity proof. The crossover
cubic leading coefficient is `121ξ³` on both engines with no surviving `189ξ³`. No
regressions, no material change.
