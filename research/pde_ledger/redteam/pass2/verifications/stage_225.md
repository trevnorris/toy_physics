---
unit_id: 225
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 225

This unit's auditor report carried verdict `clean` with **0 findings**, so there is no
directive file (expected — its absence is not a defect) and the captured diff patch is
legitimately **empty**: Codex applied no edits, both `.py` and `.wl` are byte-identical to
HEAD. Per the verifier contract this is not a re-audit; my job is to confirm that the
`clean` verdict holds — i.e. (a) the scripts still pass under their own assertions, (b) both
exec logs are exit-0 with genuine PASS bodies, and (c) the load-bearing checks are
non-tautological and the `.wl` is an independent recomputation, not a transliteration.

## Per-finding outcomes

No findings were raised by the auditor (`findings_count: 0`), so there is nothing to
classify per-finding. The remainder of this report records the confirmation checks that the
`clean` verdict survives scrutiny.

### Confirmation C1 — empty diff is legitimate, not a missed/failed apply

**Classification:** resolved (n/a — no finding)

**What changed:** Nothing. `stage_225_diff.patch` is a 1-line/empty file. There was no
directive (`directives/stage_225.md` does not exist — expected for a 0-finding unit). Both
script files were re-read in their current on-disk state.

**Assessment:** Correct and expected for a zero-correction unit. An empty diff here is the
*designed* outcome of the auditor's `clean` verdict, not a swallowed edit or a runner
failure. Nothing to address.

### Confirmation C2 — the M2 one-pole `expectNonZero` negative control is genuine

**Classification:** resolved (n/a — no finding; spot-check per batch directive)

**What changed:** Pre-existing (first pass strengthened it). `.wl` lines 165-168 add
`expectNonZero["M2 wrong -2 coefficient is not an identity", onePoleD41 - (-2 u2Base^2) D01]`,
complementing the substantive positive check at wl:164
(`expectZero[..., onePoleD41 - (-3 u2Base^2) D01]`). SymPy mirrors the positive identity at
py:72 (`assert sp.simplify(one_pole_D41 - (-3 * u2**2) * D01) == 0`).

**Assessment:** The control is **genuine and discriminating**, not a no-op. The helper
`expectNonZero` (wl:49-53) *fails* (`fail` → `Exit[1]`) when its argument cleans to zero, and
*passes* only when the residual is nonzero. Fed the wrong `-2` coefficient, the residual is
`-((D01*D2^2)/D0^2)` (mathematica log line 32) = `-u2^2 D01 ≠ 0`, so the control passes
*because it correctly detects the discrepancy*. Critically, if the correct `-3` coefficient
were substituted, the residual would be exactly `0` and `expectNonZero` would FAIL — which is
precisely what pins the `-3` uniquely and proves the companion positive check at wl:164 is not
`0==0`. The `one_pole_D41`/`onePoleD41` quantity is itself derived (substituting
`D4 -> -3 D0 u2^2` into the independently `Solve`/`solve`-derived `d41Comp`/`D41_comp`), so the
positive identity is load-bearing on both engines. Genuine negative control confirmed.

### Confirmation C3 — `.wl` is an independent recomputation, not a transliteration

**Classification:** resolved (n/a — no finding; spot-check)

**Assessment:** Independent on three distinct subsystems, consistent with the auditor's
independence call:
- First-order slopes: `.wl` extracts via series-coefficient `epsSlope[expr] =
  Coefficient[Normal[Series[expr,{eps,0,1}]], eps]` (wl:92, 134-136) vs `.py` analytic
  `sp.diff(...).subs(eps,0)/lam` (py:36-38) — different algorithms, same forms.
- Surviving-family nullspace: `.wl` uses built-in `NullSpace[mixedMatrix]` then re-gauges via
  `LinearSolve[Transpose[freeBlock], UnitVector[3,j]]` and additionally asserts
  `mixedMatrix . v == 0` (wl:398-416) vs `.py` hand block-elimination `A11.inv()` (py:354-363).
- Wall no-go: `.wl` proves structurally with `MatrixRank==2 && NullSpace==={}` (wl:360-361) vs
  `.py` `sp.solve(...)` asserting the lone `{xK:0,xM:0}` solution (py:297-300).
- Distinct symbol names (`Kwall`/`mass`/`omU` vs `K`/`M`/`OmU`) rule out a token-level port.
Not a transliteration.

### Confirmation C4 — load-bearing asserts are externally anchored (non-tautological)

**Classification:** resolved (n/a — no finding)

**Assessment:** The numeric checkpoints are anchored to literal expected values rather than
self-referential: the seven Stage-223 numbers (py:233-239 / wl:323-329), the nine Xi1
coefficients (py:273-281 / wl:340-348), the BdG sample determinant
`-5.118...e-5` (py:320 / wl:382), the null basis (py:366-369 / wl:421-431), the three
`Xi1(v_i)` (py:381-383 / wl:440-442), and the four transported windows (py:406-409 /
wl:458) all compare against hardcoded expected literals via `assert_close`/`expectClose`.
The symbolic-identity asserts (`expectZero`) reduce a derived-minus-expected residual to 0,
where the "expected" side is an independently written closed form (e.g. py:43-50 / wl:143-150
for `u4^(1)`). None collapse to `0==0`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- "Verified arbitrary-base first-order formulas" / "Verified exact first-order conservative
  compensation surface" / "On a one-pole base branch: D41 = -3*D01*D2**2/D0**2" (log 8-16).
- "All Stage 225 symbolic and numerical checks passed." (log 64), `# exit_code: 0` (log 65).
- Every `assert`/`assert_close` is reached (the script raises on failure), and the printed
  numeric deliverables match the .wl side (D0=24.2373..., Delta_BdG=-5.1189e-5, etc.).

**Mathematica:** exit=0. Notable lines:
- All M1-M8 blocks emit `PASS:` on every line; no `FAIL:`. 
- Negative control discriminates: "M2 wrong -2 coefficient is not an identity residual =
  -((D01*D2^2)/D0^2)" → "PASS" (log 32-33) — nonzero residual correctly accepted.
- `expectZero` residuals are exact `0`; `expectClose` residuals are ~1e-14..1e-18 (under tol).
- "All Stage 225 Mathematica checks passed." (log 165), `# exit_code: 0` (log 166).

**Output freshness:** Both exec logs are dated 2026-06-09T19:21:18-06:00 (the orchestrator's
direct re-run), exit 0, deterministic. Per batch context the committed `.txt` outputs are
byte-identical to the fresh run; I do not fail on committed-`.txt` mtime. Output is fresh.

## Material-change assessment

`material_change`: **false**. No script edits were applied (empty diff), so no derived result
changed. No downstream unit can be affected by a non-change. Nothing for the orchestrator to
mark stale on account of unit 225.

## Side observations (non-blocking)

- The auditor noted a pre-existing documentation lag: the card's `\stagefield{Verification}`
  line reads "Mathematica audit: none yet" and notes §8 lists only the SymPy file, both
  predating the now-present `.wl`. This is scripts-out-of-scope for me and is correctly routed
  to the orchestrator's ordinary post-batch doc-pointer sync; it understates coverage rather
  than misstating any value. Not a verification blocker.
- The `.wl` carries an extra genuine control beyond the negative one above: M6
  `expectNonZero["M6 pure-BdG determinant nonzero at sample", bdgSampleDet]` (wl:383), which
  the `.py` does not have — additional strengthening, not a defect.

## Verdict justification

The unit was audited `clean` with zero findings; there was therefore no directive and no edit
to apply, and the empty diff patch is the designed outcome rather than a failure. I re-read
both current scripts, confirmed both exec logs end `# exit_code: 0` with genuine PASS bodies
and no FAIL, confirmed the M2 one-pole `expectNonZero` negative control is a real
discriminating check (it would FAIL on the correct `-3` coefficient and passes only because
the wrong `-2` yields the nonzero residual `-u2^2 D01`), confirmed the load-bearing asserts
are externally anchored and non-tautological, and confirmed the `.wl` is an independent
recomputation across three subsystems rather than a transliteration. No script-side defect.
Verdict: **verified**, exits 0/0, `material_change: false`.
