---
unit_id: 227
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 227

## Per-finding outcomes

### F1 — paper_misalignment (notes_contradicts_script; user-resolved, direction (a))

**Classification:** resolved

**What changed:**
`notes/stages/moving_throat_pde_stage227_..._sympy_audit.md:294` — the combined `i=h`
rigidity determinant middle factor was corrected from `(251+215\pi^2)` to
`(200+147\pi^2)`. Context line 294 now reads:
`-\frac{19(-25+98\pi^2)(200+147\pi^2)(441\pi^2+4400)}{6(8670000+14894275\pi^2+2117682\pi^4)} \neq 0.`
No other notes edit. The SymPy script `.py` was NOT touched (git shows it unmodified),
matching the user-authorized notes-only direction (a).

**Assessment:**
Correct. This is the authorized `+51+68π²` slip correction (same sibling 222/223
notes-drift signature). I confirmed the stale `251` and `215` literals no longer appear
anywhere in the notes file, and the new `200`/`147` factor is present at the determinant
line. Cross-engine corroboration: the NEW `.wl` independently computes the determinant
(M4, output line 66) and emits
`(-19*(-25 + 98*Pi^2)*(200 + 147*Pi^2)*(4400 + 441*Pi^2))/(6*(8670000 + 14894275*Pi^2 + 2117682*Pi^4))`
— byte-for-byte the SymPy value (`200+147π²` middle factor), so a second engine now
corroborates that the corrected notes value is the true determinant, not the stale one.
The SymPy output `.txt` shows as modified in git, but the diff is purely a header-block
format change from the orchestrator's independent re-run; the substantive body
(determinant `200+147π²`, basis, survivors, ceilings) is identical.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage227_..._mathematica_audit.wl` (321 lines)
plus its output capture. It verifies M1–M7 with 34 asserting checks, all passing,
exit 0.

**Assessment — independence (genuinely NOT a transliteration):**
The `.wl` uses a fundamentally different decomposition than the `.py` on every axis the
anti-transliteration guard named:
- **δ mechanism:** `.py` dresses each primitive with `exp(eps*x)` then
  `sp.diff(..., eps).subs(eps, 0)` (an exponential ansatz differentiated at 0).
  `.wl` instead applies the analytic log-derivative rule directly to the closed forms:
  `deltaMixed[expr] = Total[ x·param·D[expr, param] ]` over `logPairs` (lines 110-114).
  Different mechanism, same first-order response.
- **Corridor basis:** `.py` compares the raw `M_transfer.nullspace()` vectors. `.wl`
  takes `NullSpace[transferMatrix]` then RE-PARAMETERIZES via
  `LinearSolve[tailSolveMatrix, {1,0}]`/`{0,1}` to pin tail coords (lines 228-236) — a
  distinct normalization route, not the `.py`'s `nullspace()[0]` index choreography.
- **Norm:** `.wl` builds the Gram projector natively
  (`Transpose . Inverse[gram] . Transpose`, lines 288-292) and evaluates the quadratic
  form numerically, comparing only the scalar `σ_transfer`, not symbolic intermediates.
- **K compatibility:** re-derived from the Stage-242 relation
  `B0+Z0+3(1+B2+Z2)^2/(B4+Z4)` (lines 149-153), not transliterated.
M1–M7 are each genuine and non-tautological: every numeric/symbolic expectation is
checked against an in-`.wl`-derived quantity (independent log-variation, `NullSpace`,
`Det`, `MatrixRank`, projector norm), not against the literal used to construct it.

**Assessment — M4 no-hardcode (det≠0 structural and robust):**
Confirmed. The M4 ASSERTING checks are:
- line 252: `expectTrue["...structurally nonzero", FullSimplify[detIH] =!= 0]`
- line 253: `expectTrue["...reduced matrix has rank 2", MatrixRank[reducedIH] === 2]`
The disputed factor `200+147π²` appears ONLY in a `Print` (line 251) for the verifier's
eyes; it is never an assertion operand. The no-go is established structurally (symbolic
nonzero `=!=` plus rank-2), so the conclusion is robust to F1: both `200+147π²` and the
old `251+215π²` are strictly positive in π², hence det≠0 either way. The literal is not
load-bearing.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `det[(i,h)|_pure transfer] = -19*(-25 + 98*pi**2)*(200 + 147*pi**2)*(441*pi**2 + 4400)/(6*(...))`
- `sigma_transfer = 2.3156190438605546`
- `Stage 227 audit completed successfully.`

**Mathematica:** exit=0. Notable lines:
- `OK  M4 combined i=h determinant is structurally nonzero` (condition = True)
- `OK  M4 combined i=h reduced matrix has rank 2`
- `det[(i,h)|pure transfer] = (-19*(-25 + 98*Pi^2)*(200 + 147*Pi^2)*(4400 + 441*Pi^2))/(6*(...))`
- 34 `OK` lines, 0 `FAIL`, `All Stage 227 Mathematica checks passed.`

Token convention is non-vacuous: every helper (`expectZero`, `expectTrue`, `expectClose`,
`expectVectorClose`, `expectDirectionClose`) terminates with `Exit[1]` on the failure
branch (lines 30/39/51/64/77); the trailing `Exit[0]` (line 320) fires only after all 34
checks emit `OK`. A failed check would abort with a nonzero exit, so the 34/0 OK/FAIL
tally is genuine, not print-only.

**Output freshness:** `.wl` mtime 2026-06-02 17:17:11 < its output mtime 2026-06-02
17:25:19 — output is post-fix and post-script. SymPy `.txt` was refreshed by the
orchestrator's independent re-run (header-only diff; substantive body unchanged).

## Material-change assessment

`material_change`: false.

The only edits are (a) a notes-prose correction of a determinant literal to match the
already-verified SymPy value (no derived quantity changed; the script and its output are
unchanged in substance), and (b) a NEW second-engine `.wl` that corroborates the existing
SymPy results without altering any of them. No downstream-consumed result changed value;
the MTDC-T11.3 co-loading no-go conclusion is unchanged (and now cross-checked). Nothing
downstream goes stale on account of unit 227.

## Side observations (non-blocking)

- `redteam/exec_logs/stage_227_diff.patch` is 0 bytes. The substantive changes are
  visible via `git status`/`git diff` (notes M, sympy output M header-only, new untracked
  `.wl`+output), so this did not impede verification; flagging only so the orchestrator
  knows the captured patch was empty for this unit.
- The `.wl` flips the sign convention on `v_m` implicitly (output M5 prints the negated
  direction vs the `.py`), but `expectDirectionClose` compares up to sign, so this is
  correct and intentional, not a discrepancy.

## Verdict justification

Both findings are resolved. F1's notes literal now reads `200+147π²` (the verified SymPy
value), the stale `251`/`215` are gone, the `.py` is untouched, and the new `.wl`
independently corroborates the corrected determinant cross-engine. F2's `.wl` is a genuine
second engine — a different δ-mechanism (analytic log-derivative vs exp-dressing), a
re-parameterized corridor basis, a native Gram projector, and a re-derived K — not a
port; all 34 M1–M7 checks assert against independently-derived quantities and exit 0. M4
verifies the co-loading no-go structurally (symbolic nonzero + rank-2) without hardcoding
the F1-disputed factor, so the no-go is robust regardless of F1's resolution. Helpers
genuinely `Exit[1]` on failure, so the passing run is non-vacuous. Verdict: verified.
