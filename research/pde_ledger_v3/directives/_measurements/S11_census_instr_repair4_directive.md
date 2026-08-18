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

## Secondary contamination (agent leg)

The D3/D4 `ROOT_COINCIDENCE_COEFF` OMITTED verdicts stand on their `cand[0]` (verified
genuinely uncovered: union product nonzero 6/6 generic points, equation residual zero 6/6,
defined 6/6) but their `cand[1]` is covered (union product zero 6/6) — one false missing
candidate on each line; the record-level verdicts remain genuine.

## Bounding (both legs)

Sampled spurious/probe refutations re-decided in the round-3 reviews (6 spurious, 6+8 probes)
were all backed by genuinely nonzero values and upheld — the misclassification was caught only
where nested-radical constants arise (coverage of radical inverse charts).
