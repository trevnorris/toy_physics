---
unit_id: 001
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:35:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 001

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**

`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`:

- Line 3: `Needs["VariationalMethods\`"];` inserted after `$HistoryLength = 0;`, exactly as the directive specified.
- The hand-rolled `eulerLagrange2D` helper (former lines 59-67) is fully deleted from the file (verified by absence of any `eulerLagrange2D` token).
- Section III.1 (lines 169-173): the modal-wall EL is computed via
  `elDensEq = EulerEquations[ldens, q[t, w], {t, w}];`
  `elDens = FullSimplify[elDensEq[[1]] - elDensEq[[2]], ...]`
  No manual sign flip; the `targetDens` form is unchanged from the original claim. (Iter-2 removed the erroneous `First[]` wrapper from iter-1 — the current form subtracts the two sides of the returned equation.)
- Section III.2 (lines 175-179): identical pattern applied to `lweighted`.
- Section III.4 (lines 188-192): identical pattern applied to `ldensForced`.
- Section IV (lines 205-211): `elAx` and `elAw` are now computed by `-VariationalD[lmax, axField, {x, w}]` and `-VariationalD[lmax, awField, {x, w}]`. A two-line inline comment (lines 206-207) names the convention: "VariationalD uses the standard variational derivative dL/dA - D_i(dL/d(D_i A)); the existing Maxwell residual target is the opposite-side equation residual." The unary minus is therefore a documented convention adjustment, not a silent fingerprint of porting.
- New Section I.3b (lines 128-139): an independent eigenvalue check using `SphericalHarmonicY[l, 0, theta, phi]` from Mathematica's library, with `expectZero` calls for `l=0` and `l=2`.
- Section II (lines 143-145): the documented-parallelism comment is present verbatim as required by the directive.

**Assessment:**

The directive's intent landed. Specifically:

- `VariationalMethods\`` is loaded (yes).
- `EulerEquations[...]` calls replace the deleted `eulerLagrange2D` at every site the directive named (III.1, III.2, III.4); the helper itself is gone from the file.
- `VariationalD[...]` is used for the Maxwell-component variations (IV).
- The new `I.3b` section using `SphericalHarmonicY` exists and exercises the eigenvalue claim non-tautologically: the LHS is built from Mathematica's library `SphericalHarmonicY`, not from the bespoke real-basis polynomials `y00, y20, ...`; a sign/factor error in the bespoke polynomials in I.3 would not propagate into I.3b, so I.3b is a genuine second witness for the eigenvalue claim.
- The Section II parallelism comment is present.

Non-tautology audit of the new assertions:

- `EulerEquations` returns `lhs == rhs` equation objects; the script extracts both sides and subtracts. Subtraction against the independently-constructed `targetDens / targetWeighted / targetForced` is a genuine check — the LHS comes from a different symbolic engine path (Mathematica's variational package) than the RHS (hand-written canonical form).
- `VariationalD` evaluates the standard variational derivative; the unary minus is documented as a convention adjustment relative to the `target*` form used previously. The `target*` expressions are unchanged from the original report, so the residual = 0 result is checking that Mathematica's `VariationalD` agrees with the canonical Maxwell residual form `d_mu(Z F^{mu nu}) - (1/xi) d^nu(d.A) - mu0 J^nu`, which is non-tautological.
- `lapEig[l]` applies the explicit angular Laplacian operator to `SphericalHarmonicY[l, 0, theta, phi]` and checks the eigenvalue. The check would fail if Mathematica's built-in `SphericalHarmonicY` had a different convention than the eigenvalue formula `-l(l+1)`; passing it is non-trivial.

No collateral edits beyond the directive's instructions: the sympy script is untouched (mtime April 21), only the `.wl` file was modified, and no `target*` expressions were altered.

## Exec log assessment

**SymPy:** exit=0. The sympy script was not touched; the prior log still reflects the unchanged sympy behavior. Notable lines:
- `densitized Euler-Lagrange equation = 0`
- `sourced densitized Euler-Lagrange equation = 0`
- `localized-Maxwell x-component = 0` / `localized-Maxwell w-component = 0`

**Mathematica:** exit=0. Notable lines (post-fix, captured 2026-05-21T00:26:48):
- `SphericalHarmonicY[0,0]: lap eigenvalue = 0 = 0` / `PASS`
- `SphericalHarmonicY[2,0]: lap eigenvalue = -6 = 0` / `PASS`
- `densitized Euler-Lagrange equation = 0` / `PASS` (now via `EulerEquations`)
- `weighted Euler-Lagrange equation = 0` / `PASS`
- `sourced densitized Euler-Lagrange equation = 0` / `PASS`
- `localized-Maxwell x-component = 0` / `PASS` (now via `VariationalD`)
- `localized-Maxwell w-component = 0` / `PASS`
- Final: `All Stage-001 symbolic checks passed.` `# exit_code: 0`

**Output freshness:** the Mathematica saved `.txt` output at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage001_geometry_lift_mathematica_audit.txt` has mtime 2026-05-11 12:40, OLDER than the post-fix script mtime 2026-05-21 00:22. The orchestrator-captured exec log (`exec_logs/stage_001_mathematica.log`, mtime 2026-05-21 00:28) reflects the post-fix run, so the verification is sound, but the saved `.txt` output the auditor would re-read is stale relative to the current script. Flagging this as a freshness gap for the orchestrator to refresh if downstream tooling reads the `.txt` rather than the log. (Not a verification blocker — the exec log is the authoritative post-fix record.)

## Material-change assessment

`material_change`: false.

Rationale: every `target*` expression is unchanged. Only the LHS computation route was swapped (hand-rolled EL → `EulerEquations`/`VariationalD`, plus the new I.3b second-witness check). All residuals remain identically zero. No new claim is derived; no existing claim is modified. Downstream units that depend on the Stage 001 results (modal-wall EL form, Maxwell component equations, harmonic eigenvalues) continue to see the same content. No specific downstream unit needs targeted re-audit on account of this change.

## Side observations (non-blocking)

- The saved `.txt` output file is stale relative to the patched `.wl` script (May 11 vs May 21). Not a verification blocker for this unit, but worth a one-off regeneration so any downstream consumer that reads the `.txt` sees the iter-2 content.
- The `VariationalD` unary minus is documented inline; it is now a convention adjustment with a named justification rather than a silent sign flip, which is the correct fix for the "manual sign flip is a fingerprint of porting" observation in the original report.
- Section II remains a parallel check by design, with the required documenting comment in place. The original report explicitly carved this out as acceptable.
- The orchestrator-captured `stage_001_diff.patch` reflects the iter-1 state (it contains `First[EulerEquations[...]]` calls). The actual file on disk reflects the iter-2 state (the `First[]` wrappers are gone). The exec log timestamp post-dates iter-2, so the passing residuals are the iter-2 residuals. The diff patch is therefore not current — flagging for orchestrator awareness; not a blocker.

## Verdict justification

Codex addressed all five required changes in F1 (Needs, EulerEquations substitution, VariationalD substitution, new I.3b SphericalHarmonicY check, Section II parallelism comment) and the iter-2 delta corrected the `First[]` extraction bug that prevented the first attempt from running. The post-fix Mathematica script exits 0 with every `expectZero` PASSing, every previously asserted residual is still zero, and the new I.3b check exercises the eigenvalue claim through Mathematica's built-in spherical harmonics — a genuine second witness to the bespoke real-basis derivation, not a relabeling. No regressions in the diff, no `target*` tampering, no collateral edits. F1 is resolved; overall verdict `verified`.
