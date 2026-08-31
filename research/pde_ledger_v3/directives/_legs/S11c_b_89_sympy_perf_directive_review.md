# Decision-list review — S11c-b #89 PY performance-fix directive (before any builder)

You are one of two independent legs reviewing an **orchestrator-written** directive before it goes to a
builder. The directive asks a builder to make the #89-repaired SymPy engine's **full run complete** in
tractable time **without changing any computed object** (a pure performance fix). Find defects in the
DIRECTIVE — ways it could let the builder change the physics, miss the real bottleneck, or leak the withheld
answer. Report a numbered list (file:line — problem — correction), then a one-line verdict.

## Artifact under review
- Directive: `research/pde_ledger_v3/directives/S11c_b_89_sympy_perf_directive.md`
- Premises: `research/pde_ledger_v3/directives/_measurements/S11c_b_89_sympy_perf_directive.md`

## What you are handed (all readable — you are a reviewer)
- The repaired engine (working tree): `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
  — note the uncommitted #89 repair is applied; review against the current bytes.
- The projection code (`retained_grade` L881, `first_shape_series` L713, `map_object`, `PROFILE_GRADE_SUBS`),
  the control tasks (L3033–3480), the builder caches (L70–74), and `BRANCHES`/`DENSITY_REPS` (L50–51).

## Required method
Verify every premise against the actual bytes; for any contested claim show the command (`sed`/`grep`/short
`python3`) and its literal output. A prose assertion with no command is discarded. Copy anything you execute
to `/tmp`; never modify the working tree.

## Checks (report any failure)

1. **Is the named bottleneck the real one?** Confirm `retained_grade`/`first_shape_series` is uncached and is
   invoked at volume by the controls (count the call sites and the loop cardinality). If a *different* path
   dominates the cost (e.g. a specific control's rank/`sp.series`/`simplify` unrelated to `retained_grade`,
   or a builder cache that is silently missing for a variant), the directive is aimed wrong — that is the
   finding.

2. **Could a permitted technique change the physics?** The directive allows memoization (pure-function, safe)
   and, conditionally, replacing `sp.series` with a bounded denominator expansion. Scrutinize the second:
   is the "only `eta_bg`-denominators come from the profile substitutions" premise actually true (check
   `PROFILE_GRADE_SUBS` and every `eta_bg` denominator the projection can see)? If a bounded expansion could
   silently differ from `sp.series(·,eta_bg,0,2)` on some retained-grade term, the directive must forbid it
   unless the §4 equivalence check covers that case — verify the §4 check is strong enough (does one operator
   row + one kernel scalar exercise the risky denominators?).

3. **Is the equivalence self-check load-bearing and non-tautological?** §4 requires keeping the current
   `first_shape_series` as a reference and emitting `reference(x) − new(x) = 0`. Confirm this is a genuine
   two-route check (the reference must be the *unmodified* original, not the new function under another name)
   and that a memoization-only change would still exercise it. Flag if the check could pass while the
   projection changed (e.g. sample inputs that miss the denominators, or a reference that is not actually the
   old code).

4. **Answer leak (rule 5).** The corrected live basis count must not appear as a target; only the public
   frozen `26` may. Grep the directive and premises. Confirm the saved partial primaries (which contain the
   corrected value) are NOT handed to the builder as an acceptance anchor.

5. **Does it forbid weakening coverage for speed?** The controls must keep all cases/anchorings/reps. Confirm
   the directive forbids reducing control cases, and that none of its suggested techniques implicitly drops
   coverage (e.g. "reuse cached objects" must not skip a genuinely-distinct variant).

6. **Completeness / termination.** Will the suggested techniques plausibly make the full run finish, or is
   there a second bottleneck the directive ignores that would just move the stall (e.g. if memoization fixes
   `retained_grade` but a control still rebuilds an expensive kernel per case with a distinct cache key)? If
   the fix as scoped would not actually let the run complete, say so.

## Physics filter
Report a finding only if it catches a way the fix would change the physics, mislead the builder, leak the
answer, or fail to make the run complete. Not style.

## Output
Numbered discrepancies, then `VERDICT: CLEAR` or `VERDICT: N defects`. Nothing else.
