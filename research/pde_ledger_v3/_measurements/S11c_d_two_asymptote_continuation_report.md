# Two-asymptote continuation: end lifts saved, native binding failed

The authorized continuation finished on 2026-09-28 at 07:23:39 UTC with a
code error in `bind-native-action`. It saved both local regular inverse
matrices and all eight first-grade end fields before stopping. The intended
`END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED` stage was **not completed**: native
forcing, carrier reconstruction, mutations and final domain assembly were
not reached. No retry, new scientific job or external review was launched.

The [completion checkpoint](S11c_d_two_asymptote_continuation_checkpoint.json)
pins the actual results, logs and complete indices. Original files remain at
`_scratch/s11c/s11c-d-two-asymptote-20260927/continuation/complete`.
This inspection used JSON, source text and file hashes only; it restored no
scientific pickle and did not independently recompute mathematical results.

## Saved progress

- **536 completed journal operations:** 61 original and 12 diagnostic returns
  restored byte-for-byte; 413 new operations; 50 exact-operand reuse records.
  Every completed input, return, started receipt and completed receipt was
  checked against the operation index. The 26 copied supporting artifacts
  match their saved originals. No completed parent function was replayed.
- **50 denominator certificates**, covering every entry in both pole blocks.
  All original `is_zero` properties were `None`. The actual certificates use
  exact real-component reconstruction and two same-sign terms: 41 positive,
  nine negative. Their five distinct displayed expressions have the stated
  signs; all saved reconstruction residuals are zero. The first entry reuses
  the diagnostic series and matches its actual saved coefficient.
- **Two regular inverse matrices and eight end-field bundles** are saved.
  Sixteen end-equation power checks and eight inverse/residue checks report
  exact zero 5-by-5 matrices, with no unresolved entries. The 17 inherited
  source/grade/end joins also remain zero. Their complete operands and
  returns are indexed and hashed. These are worker-reported checks on partial
  outputs, not an independent result-validation verdict.
- **2,366 files / 9,597,308 bytes** are preserved, including the incomplete
  operation input. All 447 source/input records and 15 source snapshots
  match their pinned hashes. The original 279 production artifacts and
  59 diagnostic artifacts remain unchanged.

The checkpoint lists each end-field and regular-matrix artifact with its
exact route and SHA256. It also pins the complete operation and artifact
indices without duplicating their large payloads into the report.

## Actual failure and small repair

The sole incomplete operation is `bind-native-action`. Its saved input is
`complete/operations/0536-bind-native-action/input.pickle` (79,669 bytes),
SHA256 `8e1cf051c7d8b47ae8ae775d63f93845ee995dfa939ff7062dfbc91cff92562c`.
No return or completed receipt exists for that operation.

The exception is an unbound profile limit:

```text
Limit(s11cdWProfile(_liftBound_s11cdProfileCoordinate),
      _liftBound_s11cdProfileCoordinate, -oo, dir='+')
```

Source inspection explains the mismatch. `prepare_bound` calls
`bind(alpha_safe(v))`. `alpha_safe` renames an integration variable to a
Dummy, including its occurrence inside this endpoint Limit. `bind` looks up
endpoint values by exact expression identity. The renamed Limit misses that
lookup, then correctly reaches the unknown-Limit refusal. The
[existing source summary](S11c_d_uniform_source_focused.json) already records
the corresponding unrenamed W-profile endpoint as zero. No new endpoint
evaluation, regulator removal or physics assumption is needed.

The proposed local correction is to substitute the already accepted exact
endpoint map before alpha-renaming, inside `alpha_safe`:

```python
return visit(expression.xreplace(endpoint_values), {})
```

This replaces its current `return visit(expression,{})`. It does not evaluate
an Integral or Limit, alter an endpoint value, or relax the unknown-Limit
guard. This is a source-level repair proposal; the failed worker and all its
pinned files remain unchanged. The correction has not run on scientific
objects.

A separately authorized continuation should restore all 536 completed
returns, both regular matrices, the eight complete end-field bundles and
their supporting outputs, then resume **only** the saved native-binding
input and subsequent unfinished forcing/control/domain construction. It
must not reconstruct the completed Laurent coefficients or end fields.
Keep the ordinary 900s outer/840s native guard and a new result directory
with the existing session's completion hook. This is the next concrete
dependency, not a new general certification or review campaign.

## Actual resources and completion checks

Worker time was **640.024100 seconds**; whole-guard time was 641.232967 seconds.
Peak whole-job memory was **139,350,016 bytes** (about 133 MiB). All 322 resource
samples show zero swap and zero memory events; at most three tasks were
sampled. Minimum host available memory was 17,934,106,624 bytes. The actual
2 GiB memory / zero-swap / one-CPU / nice15 / 32-task / one-native-thread
limits and both deadlines were verified. This was not a resource or timeout
failure; the guard's stop reason is null.

Supervisor, child, guard and coordinator all record exit 1. Scientific stdout
is empty and no checks JSON exists, so stdout/checks identity is unavailable.
Scientific stderr contains 4,424 bytes: the saved failure traceback plus the
outer `main()` frame. Surrounding guard, coordinator and hook stderr are
empty. The completion hook was armed before science and finished normally.
These facts preserve the failure rather than converting process termination
into scientific acceptance.

## Scope retained

The native census remains 25 cells, 71 local terms and 160 nonlocal terms.
The copied plane-jet premise remains explicitly unverified. No forcing or
formal-carrier mutation passed in this run, and no nonzero integrated
response is claimed. Full two-asymptote response, Green/FORM, A11/A12,
forcing pairing, coincident-point extension and retarded equivalence remain
open. The unchanged physical slice remains evanescent.

Literal reviews remain Claude **NEEDS REVISION** and Grok **CLEAR FOR THIS
BOUNDED INGREDIENT**, with the earlier local dispositions; there is no fresh
independent CLEAR. Earlier accepted reference/outgoing results, all failures,
incident history, Lean work and shared guard were preserved. The protected
builder suffix still hashes to
`f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2`.
The concise new checkpoint is below 1 MiB; large binary/index scratch
artifacts remain in place.
