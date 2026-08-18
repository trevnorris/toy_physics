# Measurements — S11 census-instrument repair round 4 directive

Round-3 scoped legs (2026-08-18): fresh-agent report (scripts `/tmp/s11_r3_leg/`, archived) and
`~/.s11_build/census_repair3_grok_leg.txt` (scripts `/tmp/s11_r3_leg/` grok set). Both legs
certified all six round-3 repairs live (per-class ablation diffs), the round-3 calibration
able-to-fail at `89ed80c9`, the reducer recount exact (917/348/815), all 104 surviving witness
failures membership-driven, and converged on ONE finding.

## The false OMITTED (both legs; orchestrator-confirmed)

`WL_S11_XKIN_ANISO_D2_ROOT_COINCIDENCE_COEFF_SOLUTION`, round-3 `containment_wl.stdout`
verdict `OMITTED_BRANCH`, both candidates `NOT_COVERED_SAMPLED`. The candidates are inverse
charts (`bComp` in terms of `muR`) of the same quadratic whose `muR`-roots the record emits.

```
$ python3  # emitted muR-roots r1,r2 from the committed record; census candidate cand0;
           # union product (muR-r1)(muR-r2) with bComp -> cand0, generic rational point
union product at candidate0, generic point: 0.e-137
```

## The mechanism (agent leg; orchestrator-reproduced)

At the sampler's point `{muR:1, sRho:2, k1:-1, k2:-1/2}` the covering branch's constraint is an
exact algebraic zero that the instrument classifies NONZERO:

```
$ python3  # all four candidate-root pairings through s11_census_math.simplify_residual at fd9a5835
cand0 root0: status=NONZERO numeric=0.54869684 + 0.8923707*I     # genuinely nonzero
cand0 root1: status=NONZERO numeric=-7.112828e-161 - 1.778207e-161*I   # EXACT ZERO, misclassified
cand1 root0: status=NONZERO numeric=0.54869684 - 0.8923707*I
cand1 root1: status=NONZERO numeric=-7.112828e-161 + 1.778207e-161*I   # EXACT ZERO, misclassified
```

Structural test `simplified != 0` at `s11_census_math.py:550`, `:570`; sampled fallback
`_sample_union_coverage` `:689-717` reached because the symbolic union product exceeds the
2000-char simplify cap (`:544`); promotion at `:829`.

## D3/D4 candidate quality — CORRECTED by both round-4 directive legs

The round-3 agent leg's "cand[1] is covered (union product zero 6/6)" claim was an artifact of
multiplying the constraints WITHIN the single emitted branch `{bComp -> 0, muR -> 0}` —
`bComp * 0 = 0` — which turns the branch's AND into an OR. Both round-4 directive legs refuted
it and the orchestrator confirmed:

```
$ python3  # D3 equations from the record; candidate {muR: 0, sRho: -k1^2/(k2^2+k3^2)}, bComp=1/2
candidate on the variety? residuals at generic point (bComp=1/2):  0 / 0 / 0
emitted branch {bComp:0, muR:0} covers this point? bComp there = 1/2 (needs 0) -> NOT covered
```

Both D3/D4 missing candidates are genuinely uncovered: the records each carry TWO omitted
families (the registered `{muR == bComp, sRho == 1}` class and the `{muR -> 0,
sRho -> -k1^2/K_perp^2}` family). The record-level OMITTED verdicts are genuine and
STRONGER than round 3 reported.

## Bounding (both legs)

Sampled spurious/probe refutations re-decided in the round-3 reviews (6 spurious, 6+8 probes)
were all backed by genuinely nonzero values and upheld — the misclassification was caught only
where nested-radical constants arise (coverage of radical inverse charts).


---

# Round-4 directive-leg measurements (Codex + Grok, 2026-08-18)

Leg reports: `~/.s11_build/repair4_directive_codex_leg.txt` (scripts
`/tmp/s11_r4_directive_review/`), `~/.s11_build/repair4_directive_grok_leg.txt` (scripts
`/tmp/s11_r4_review/`). Convergent findings, folded into the directive at this commit:

1. Required 1 certificate-strength nondeterminacy (Codex): two compliant implementations
   (fixed-precision interval vs minimal-polynomial) differ on
   `Sqrt[10^200 + 1] - 10^100` (exact nonzero lower bound `1/(2*10^100+1)`; a 50-digit
   rigorous interval still contains zero) — POLICY_A SPURIOUS_BRANCH vs POLICY_B
   BRANCH_MEMBERSHIP_UNDECIDED on the same byte-shaped payload. Fold: exact route mandatory;
   new calibration plant with that shape.
2. Reading-B union object (Grok finding 1 = Codex finding 2): per-branch AND semantics named
   explicitly; within-branch constraint products forbidden. Grok's D3 evidence:
   `missing1_prod_all=0 status=ZERO` (false cover) vs one-from-each combination
   `product=bComp status=NONZERO_SAMPLED` (uncovered).

Codex verification that stands: D2 minpoly certificates (`_x` for both zero pairings,
`729*_x**2 - 800*_x + 800` with certified lower bound 1/2 for the nonzero pair); D2 symbolic
union products simplify to zero exactly; UNDECIDED propagation truth-table sound; round-3
STRATUM5/6 witness conjuncts unaffected by Required 1 (five TRUE, three contingent, no FALSE);
the genuine-refutation plant feasible (`NOT_COVERED_SAMPLED` at a certified-nonzero point).
