# Build review (leg 2, re-run) — S11c-b #89 PY §3a repair engine (+ lever-C projection)

You are ONE of two independent legs reviewing a built SymPy engine. Derive everything from first principles;
do NOT trust the engine's own emitted conclusions un-recomputed. You have no prior context and must not assume
any expected outcome.

## Read and execute the FULL review spec
The complete artifact description and the checks (A1–A4 physics, B5–B7 lever-C, C8 soundness) are in:
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_b_89_sympy_3a_repair_build_review.md`
Read it in full and execute EVERY check in it. Working dir: `/var/projects/toy_physics`.

## ⛔ COMPLETION & EVIDENCE CONTRACT — this is why leg-1 was rejected
- Do NOT output any VERDICT until EVERY check (A1, A2, A3, A4, B5, B6, B7, C8) has run to completion and you
  have SAVED its script AND its literal stdout to `/home/trevnorris/.s11_build/S11c_b_89_buildleg_grok2/`
  (mkdir it first). An early VERDICT with unfinished checks is a FAILED review, not a CLEAR.
- In your final report, list each check with (a) the absolute paths of its saved script+stdout and (b) the
  literal key number(s) it produced. A claim with no saved script+stdout is discarded — a prose "I derived X"
  counts for nothing here.

## ⚠ EFFICIENCY — do not time yourself out (this is what killed leg-1)
- Part A form ablation (MANDATORY): copy the engine `.py` to /tmp, zero the Hessian / second-background-jet
  map, re-run, and report the literal basis count BOTH ways. The count MUST move under the ablation; if it
  does not, the repair did nothing — that is a finding.
- Part B (lever-C physics identity): do NOT rebuild the full energy density — it is expensive and will time you
  out. Instead extract a BOUNDED sample of already-constructed operator and kernel scalars and compare
  `first_shape_series_reference(x)` vs `first_shape_series_fast(x)` NUMERICALLY (substitute random rationals;
  symbolic `simplify` crashes on these). Any genuinely nonzero difference is a physics bug — report it with the
  scalar. Also confirm Integral/Derivative-bearing scalars route to the exact reference.

## ⛔ Do not modify the working tree
The `.out` is a read-only git-annex pointer (content is present locally — read/grep it directly). Ablate a
/tmp COPY of the engine `.py`, never the working tree.

## Output
Numbered discrepancies (file:line — problem — how found, each with your saved script+stdout absolute paths and
the literal number), then a final line `VERDICT: CLEAR` or `VERDICT: N defects`. Nothing else.
