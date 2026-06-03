---
unit_id: 224
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 224

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
New `mathematica/moving_throat_pde_stage224_..._mathematica_audit.wl` (218 lines) plus its
committed `mathematica/output/...txt` transcript. Covers M1–M6:
- M1 (`.wl:93-110`): grouped inverse derived by `Solve` of the forward 3-vector system
  `{p20,p21,p22} == laneMatrix.{Pbar,aP,bP}`, then confirms the solved `Pbar/aP/bP` recover
  the input seeds.
- M2 (`:112-117`): isotropic ceiling rearrangement residual == 0.
- M3 (`:119-131`): WA compiler — but recovered via `LinearSolve[laneMatrix, waLanes]`, then
  `aPwa = eps Pbar Xi1/4`, `bPwa = 3 eps Pbar Xi1/4`, `bPwa == 3 aPwa`.
- M4 (`:133-138`): `laneMatrix.{Pbar,aPwa,bPwa}` re-expands to the WA lanes.
- M5 (`:140-181`): sign-case `Refine[Abs[eps xi1], eps xi1 > 0/<0]` (NO Abs dictionary), plus
  `expectTrueUnder` lane-dominance and the `Pbar(1+|eps xi1|) == Pbar + 4|aP|` robust form for
  both signs.
- M6 (`:183-214`): carried compat point + four ceilings (exact `Rationalize` via `exactDecimal`),
  computes budgets, asserts the two defining relations per ceiling.

**Assessment:**
Genuinely independent, not a transliteration. The `.py` asserts the inverse projector
coefficients `(P20+2P21+2P22)/5` etc. as pre-chosen forms (`.py:30`); the `.wl` instead derives
the inverse by `Solve`/`LinearSolve` on the forward system and reads the solution off
(`.wl:95-100, 123`) — a different decomposition. The `.py` resolves `Abs` by a
`subs({Abs(eps*Xi1): zabs})` dictionary (the fragile trick the auditor flagged); the `.wl`
resolves it by `Refine` under a sign assumption (`.wl:144-145`), exactly the anti-transliteration
guard the directive demanded, and the log shows `Refine[Abs[eps xi1], eps xi1>0] -> eps*xi1` /
`-(eps*xi1)`, confirming the case-split actually fired. All six manifest items map to PASS lines;
none is a `True==True` placeholder — M1/M3/M4/M5 residuals are simplified symbolic expressions that
reduce to 0, M5 dominance is a genuine inequality `FullSimplify`d to `True` under the sign
assumption, and `bPwa==3aPwa` / the robust-form identity are non-trivial. 24 PASS / 0 FAIL.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
`.py:152-159` (old): eight `assert_close(float(budgets[k][i]), <hardcoded budget literal>)`.
`.py:153-156` (new): a loop that, for each ceiling, asserts
`(eps_xi_budget + 1)*barP0_compat == Pcrit_val` and `barP0_compat + 4*a_budget == Pcrit_val`,
with `Pcrit_val = ceilings[key]`. The eight pre-baked budget literals are gone (confirmed by the
diff: 8 deletions, replaced by the 4-line loop). A `# carried from Stage 223` comment was added
at `:133`. The four ceiling literals, the compat point, and the printed headroom numbers
(`:149-151`) are unchanged — confirmed against the SymPy exec log (lines 36-39 print the same
`0.36793.../0.73761.../2.94889.../4.63505...` budgets).

**Assessment — F2 fix is GENUINE, not still-hardcoded.**
The old check was `budget == <second copy of the answer>`: a budget-formula typo could pass as
long as the duplicated literal carried the same typo, and the literal could silently drift from
upstream. The new check ties each budget back to the ceiling via its *defining relation*:
`eps_xi_budget` and `a_budget` are computed from `Pcrit_val`/`barP0_compat`, and the assertion
re-derives `Pcrit_val` from them. A wrong budget *formula* (e.g. forgetting the `-1`, or `/3`
instead of `/4`) now fails the assertion; only the correct formula closes the loop. It is
`f(budget, inputs) == ceiling`, not `budget == budget`. This is exactly the form the directive
specified and the auditor pre-validated in its self-test (report line 127). The residual concern
— that the relation is an algebraic identity for these definitions, so it cannot catch a wrong
*ceiling/compat value* — is acknowledged and acceptable: F2 is a low-severity hardening to remove
the duplicate-answer gap, not to re-verify upstream numbers (those are stage 223's job). No new
constant introduced; no value altered.

## Notes-renumber (authorized, notes-only)

Out of scripts-only scope to read the notes file directly, but the directive's
`## Applied: notes-renumber` block records: only
`notes/stages/...stage224..._sympy_audit.md` touched; stale Stage 240/241 labels → canonical
Stage 223/224; math/values unchanged; scripts/paper/appendix untouched. The diff patch contains
no notes-file hunk and no script/paper hunk beyond the F2 `.py` edit, consistent with a
notes-only, content-neutral renumber. No script-side impact.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `bP = 3*Xi1*eps*(Delta_norm + T_quad)/(4*mhat0**2)` / `Exact relation: bP = 3 aP` (log 26-27).
- Printed budgets unchanged (`both_10: |eps Xi1| <= 0.36793...`, log 36-39).
- `All Stage 224 symbolic and numerical checks passed.` / `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines:
- `M1 solved Pbar recovers input residual = 0 / PASS` (Solve route fired; log 14-20).
- `LinearSolve recovered {Pbar, aP, bP} = {(deltaNorm+tQuad)/mhat0^2, (eps...xi1)/(4 mhat0^2), (3 eps...xi1)/(4 mhat0^2)}` then `M3 bP = 3 aP ... PASS` (log 32-38).
- `Refined Abs[eps xi1] for eps xi1 > 0: eps*xi1` (Refine, not dict; log 53-54).
- 24× `PASS`, 0× `FAIL`; `All Stage 224 Mathematica checks passed.` / `# exit_code: 0`.

**Output freshness:** confirmed. `.py` mtime 16:30:28, `.wl` mtime 16:32:26; both committed `.txt`
outputs mtime 16:38:06 — newer than their scripts. Both transcripts regenerated post-fix.

## Fabricated-literal / tail-noise assessment

Confirmed genuine floating-point/bisection tail noise, not a masked discrepancy. Stage 224's
hardcoded literals vs upstream stage 223's saved output:

| key | stage 224 literal | stage 223 source | first divergence |
|---|---|---|---|
| compat | 0.002069792318062885 | 0.002069792318062883 | 16th sig fig |
| both_10 | 0.0028313316855593175 | 0.0028313316855593336 | ~14th–15th sig fig |
| one_10 | 0.0035965105896846573 | 0.003596510589684656 | ~16th sig fig |
| both_30 | 0.00817339430971383 | 0.008173394309713826 | ~17th sig fig |
| one_30 | 0.0116633929790174 | 0.011663392979017383 | ~17th sig fig |

Largest absolute divergence (`both_10`) ≈ 1.6e-17, i.e. ~6e-15 relative — far inside the
`assert_close` default `tol=1e-12` (`.py:7`) and stage 223's own `tol=5e-13`. These are
double-precision/bisection-tail artifacts of re-typing upstream numbers, exactly as the auditor
characterized; they do NOT constitute a value conflict. Importantly, the F2-rewritten checks do
not compare 224's literals against 223's — they only verify each literal is internally consistent
with its own derived budget, so the tail noise is not even exercised as a cross-stage tolerance;
no real discrepancy is masked. (Side note: tighter traceability would re-type the exact upstream
strings, but the directive explicitly instructed Codex NOT to alter the carried values, so this is
correct behavior, not a defect.)

## PASS-count check

Committed Mathematica output: `grep -c "^PASS:"` = 24, `grep -c FAIL` = 0. Matches the orchestrator's
reported "24 PASS, 0 FAIL." SymPy block ends "All Stage 224 ... checks passed" with exit 0.

## Material-change assessment

`material_change`: **false**.

No derived result changed. F1 adds a second engine that re-derives identical symbolic content
(all residuals == 0; recovered forms match the `.py` to the symbol). F2 swaps one set of passing
assertions for another, equally-passing set; the four carried ceiling/compat literals and all
printed headroom numbers are byte-for-byte unchanged. The notes renumber is content-neutral.
Nothing a downstream unit could depend on (no constant, no carried numeric, no symbolic identity)
was altered. Downstream units do not need re-audit on account of unit 224.

## Side observations (non-blocking)

- The F2 inverse-relations are algebraic identities of the budget definitions, so they cannot
  detect a corrupted *ceiling* or *compat* literal (only a corrupted budget formula). This is the
  intended, directive-sanctioned scope of a low-severity hardening; flagged only for completeness,
  not blocking.
- Minor traceability gap persists: 224's literals are re-typed copies of 223's, divergent in the
  trailing 1–2 ULPs. Within tolerance and explicitly out of scope for this fix (directive forbade
  altering carried values). Not a finding.

## Verdict justification

Both findings are `resolved`. F1: a genuinely independent Mathematica `.wl` derives M1–M6 via
`Solve`/`LinearSolve` and sign-case `Refine` (a different decomposition than the SymPy round-trip,
no Abs-substitution dictionary), 24 PASS / 0 FAIL, exit 0. F2: the self-referential
`budget == hardcoded-literal` checks are replaced by defining-relation assertions tying each budget
to the ceiling via `Pcrit_val`/`barP0_compat`, so a wrong budget formula now fails — genuine, not
still-hardcoded. The literal-vs-upstream differences are confirmed 14th–17th-sig-fig tail noise
well inside tolerance and not even cross-compared by the new checks. Outputs regenerated post-fix.
No regressions in the diff, no material change. Verdict: **verified**.
