# S11 WL engine fix round 1 — decision list (folded once after two legs)

Target: `research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl` at HEAD
of this commit. Companion measurements (generated): `../_measurements/S11_wl_engine_fix_round1_brief.md`.

## The measured wall (every claim here is twin-backed)

- XKIN_ANISO D4 was killed at 79,127 s wall, the last 16.2 h with zero emissions, inside the
  `Reduce[And @@ equations, variables, Complexes]` at `emitLocus` (wl:211) — between the
  `_IDENTICALLY_SATISFIED` emission (wl:205) and the `_INCONSISTENT` emission (wl:216).
- The engine has three quantifier-elimination call sites: the wl:211 `Reduce`, the per-branch
  `Resolve[Exists[…], Reals]` real-existence test (wl:137), and the rank-decision
  `Resolve[…, Reals]` (wl:679). The completed D3 record implies 700+ executions of this class per
  XKIN cell (164 `_INCONSISTENT` + 164 `_REAL_STATUS` records + 379 `_REAL_ADMISSIBLE` branch
  entries; 140 of the paired records are `NOT_APPLICABLE` and ran nothing). XKIN_ANISO D2's cell
  was killed by the memory guard at 891 MB available, 1,275 s from a cold kernel, with its last
  emission inside its strata blocks — the rank-decision site's territory.
- The diagnosis (two legs, scripts + literal stdout archived) measured: the stuck system's
  operands grow modestly D3→D4 (rationalized minors 9–12 → 52–74 terms, degree 10 → 14) while the
  K-direction Gröbner basis of the rationalized system completes in 1–10 s where `Reduce` did not
  return in 15 h. ⚠ The COEFF-direction Gröbner probe itself timed out at 600 s in SymPy: the wall
  is NOT uniformly cheap, so the repair must be shaped as bounded attempts plus certified
  fallbacks, not a one-for-one primitive swap.

## Out of scope — two spec-compliance defects found by this round's review legs

Registered in `research/pde_ledger_v3/DEFECT_REGISTER.md` (entries dated 2026-08-16): the
pointwise `_IDENTICALLY_SATISFIED` route (wl:203) and the joint-system `ROOT_COINCIDENCE` locus
(wl:887–890) versus the spec's per-pair obligation. ⛔ This round must NOT touch either — they are
committed-record defects in BOTH the regression baseline and the wall repair's scope would grow;
they get their own round. Obligation 5's byte-identity is measured against the committed records
as they are.

## What must be true after the fix

1. **Objects unchanged, routes free, nothing weakened.** Every locus-protocol record remains the
   object `directives/S11_SHARED_PHYSICS.md` §5 names, at least as strong as its committed form,
   with the pinned payload field sequences and token vocabularies. `_IDENTICALLY_SATISFIED` and
   `_INCONSISTENT` remain computed from `_EQUATIONS` (spec line 251: never read off any solver's
   return, including this cell's own `_SOLUTION`).
2. **One uniform route at every QE site.** The repaired logic is the same at all three call sites
   (wl:211, wl:137, wl:679) and selects within itself only on measured budget expiry —
   ⛔ never on package, dimension, root index, stratum index, or tag identity. The diff is the
   evidence; a route reachable only by some cells' identity is a build failure.
3. **Honest incompleteness, but only after a real attempt.** An undecided token may stand only
   after the route's bounded deciding attempts have all run and failed within their budgets, and
   the record carries the live operands of the failed attempts. No token is ever chosen to make a
   cell complete. ⭐ Acceptance probe (able to fail): every undecided-class locus record in a newly
   completing cell is re-attacked post-build by a bounded, independent, exact-rational decision
   probe on its own emitted `_EQUATIONS` and premises; ⛔ any such record the cheap probe decides
   is a build failure.
4. **Completion is defined by progress, not just a clock.** XKIN_ANISO D3 and D4 run to cell
   completion under the per-cell contract (fresh kernel), with no inter-emission silent gap
   exceeding 1,800 s and an outer per-cell wall of 14,400 s as this round's stop: a cell still
   emitting at the outer wall is a measurement to report, never a failure to paper over. **D2 is
   attempted under the same contract with the memory guard armed.** A memory-guard death is a
   measurement, not a build failure — ⛔ but its death record must be analyzed: if the last-emitted
   tag places the death inside a QE call the fix claims to have bounded, the fix did not apply and
   the build fails.
5. **The 19 completed cells are the regression bar, run for real.** Re-run every completed cell
   with the fixed engine and compare ALL emitted `WL_S11_*` tags — including the `LOCAL_` premise
   inventories, which the spec names as objects — against the committed `.out` baseline:
   byte-identical. Reading the committed file runs no engine; the comparison operand is a fresh
   run's output.
6. **Newly completing cells face a manifest census, not a suffix count.** The expected tag surface
   for a new XKIN cell is derived from the completed XKIN D3 record: every tag FAMILY present at
   D3 (names reduced modulo root index, stratum index, branch index, and dimension-dependent
   counts) must be present, plus every §8 named object for the cell. A family present at D3 and
   absent in the new cell is a build failure. Emitted objects are additionally spot-verified by
   exact rational specialization against their own `_EQUATIONS` — an oracle independent of the
   route that produced them.
7. **The measurements twin is generated by script**, never transcribed, and every measured claim
   in the build report carries its command and literal output.

## Builder operational constraints

- One kernel at a time (two-seat licence; nothing else running). Absolute paths everywhere.
  Transcript outside the repository. Never run the full sweep in one kernel — per-cell only.
- Long runs must show observable progress; long + silent is the failure mode under repair.
- The working tree is committed before the build touches it.
