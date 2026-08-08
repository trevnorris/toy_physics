# W0 — integrator statement: R4 binding and coverage

⛔ Restricted to the orchestrator and review legs. Do not give this document to an engine builder.

⚠ Do not read `research/pde_ledger_v3/steps/S9_PILOT_ADJUDICATION.md`.

## Launch status — blocked

⛔ W0 cannot launch in the current worktree. S10-Python reaches its registry comparison and exits 1 because
it attempts to iterate a declared dimension that is now a `BoundDimensionLaw`:

```text
TypeError: 'BoundDimensionLaw' object is not iterable
```

The first builder regression therefore cannot be established for that engine. Do not issue any of the four
builder assignments until the work described by `W3_fix_round2.md` has landed and all four unchanged
engines have regenerated their committed deterministic baseline projections.

## R4 and its current wired residuals

The registry defines R4 as

```yaml
residual:
  - Sub
  - [Q, Q.brane.c_gamma]
  - [Sqrt, [Div, [Q, Q.brane.mu_R], [Q, Q.brane.rho_br]]]
```

`Q.brane.c_gamma` has LTM dimension `[1,-1,0]`. The two current bindings do not supply that object:

- S9-WL binds all three QIDs through `WL_S9_CANDIDATE_SPEED_SQUARED1`; the committed value bound as
  `Q.brane.c_gamma` is `muR/rhoBr`. The instantiated residual is
  `muR/rhoBr - Sqrt[muR/rhoBr]`.
- S10-WL binds `Q.brane.rho_br` to `WL_S10_MAIN_D3_Q6_INERTIAL_COEFFICIENTS` with
  `sequence_first`, `Q.brane.mu_R` to `WL_S10_MAIN_D3_Q6_STIFFNESS_COEFFICIENTS` with `list`, and
  `Q.brane.c_gamma` to `WL_S10_MAIN_D3_Q3_DISTINCT_ROOTS` with `sequence_second`. The selected committed
  value is `((k1^2+k2^2+k3^2)*muR)/rhoBr`, so the instantiated residual is
  `(k1^2+k2^2+k3^2)*(muR/rhoBr) - Sqrt[muR/rhoBr]`.

S9 is wired to a squared speed. S10 is wired to a squared frequency carrying
`k1^2+k2^2+k3^2`. These are different binding defects.

## What each engine already emits

| Step/engine | Existing squared-speed object | Missing |
|---|---|---|
| S9-WL | `WL_S9_CANDIDATE_SPEED_SQUARED1 = muR/rhoBr` | linear-speed object and two defining links |
| S9-Python | `PY_S9_MAIN_SPEED_SQUARED_CANDIDATES = [mu_R/rho_br]` | linear-speed object, two defining links, and R4 row |
| S10-WL | `WL_S10_MAIN_D{2,3,4,5}_ROOT2_Q6_RATE_COEFFICIENT = muR/rhoBr` | selected linear-speed object and two defining links |
| S10-Python | `PY_S10_MAIN_D{2,3,4,5}_ROOT2_Q6_ROOT_OVER_WAVENUMBER_NORM = mu_R/rho_br` | selected linear-speed object, two defining links, and R4 rows |

For both S10 engines and every listed MAIN cell, the corresponding root-1 quotient is `0`. Selection is
therefore load-bearing.

## Builder-artifact adjudication

After the launch hold is cleared, review all four build artifacts before changing a harness configuration.

1. Treat every name-and-value manifest as disclosure and review input, not as an acceptance oracle. Check
   every manifested value against the required object, placement grain, computed operands, and premises.
2. In each of the ten engine/cells—one S9 MAIN cell per engine and S10 MAIN `D=2,3,4,5` per engine—verify
   that membership equals the emitted predicate truth set. Then compute the cardinality verdict here. The
   integrator requires exactly one member in every cell; that number is not builder feedback.
3. For the sole member of each accepted cell, verify the root association, both operands and guard of the
   quotient link, the scalar `[1,-1,0]` dimension witness, the sign premises, and both operands and guard of
   the square link. Reject a pinned failure or any missing association.
4. Verify provenance from the upstream-quotient sentinel artifact and verify selection from the
   fixed-spectrum selection-only artifact. Four engines agreeing on a value is not a substitute.

The equality `L_T=0` carries no information about whether the right object was selected, derived, or bound:
any implementation that defines `c_T` as a square root of its own `v_T^2` makes it zero. The quotient link
is likewise a defining consistency object, and the in-scope quotient is algebraically indistinguishable
from a recombination; the upstream mutation supplies the provenance evidence. The instantiated **R4
residual** is the residual here that can be nonzero on a wrong harness binding and adjudicates the registry
relation. It remains in this integrator-only statement.

## Current coverage holes

There are exactly two `registry_residual` rows repository-wide: one in `checks_S9.yaml` and one in
`checks_S10.yaml`. Both name R4 and both declare `engine: wl`.

- Python does not run R4 at either step.
- S9 reaches only its MAIN three-component cell.
- S10 reaches only MAIN at `D=3`; MAIN `D=2,4,5` do not enter R4.

The four emission builds do not make R4 tested. They only create eligible linear-speed objects. Repointing
and coverage expansion are a separate reviewed change to two closed-step configurations.

## Required repoint

Only after all ten builder placements pass the adjudication above:

1. In S9, repoint the existing WL `Q.brane.c_gamma` binding from the squared-speed emission to the scalar
   `c_T` emission associated with that cell's sole predicate-true root. Add the corresponding Python R4
   row, binding the two material QIDs from computed engine emissions and `Q.brane.c_gamma` from Python's
   corresponding scalar `c_T` emission.
2. In S10, create R4 rows for both engines at MAIN `D=2,3,4,5`. Bind the material QIDs from each cell's
   emitted coefficient objects and bind `Q.brane.c_gamma` directly to the scalar `c_T` emission associated
   with that cell's sole predicate-true root. Do not bind a certificate collection, root list, selector
   into a root list, squared frequency, or squared speed as `Q.brane.c_gamma`.
3. Resolve actual tag spelling only from the reviewed name-and-value manifests. The placement is fixed by
   cell and selected-root association; spelling may differ by engine.
4. The resulting reviewed coverage is two S9 rows (one per engine) and eight S10 rows (two engines across
   four dimensions). For every row, emit the instantiated R4 operands and residual before its guard, and
   permit a nonzero residual to fail visibly.
5. Keep the configuration repoint separate from the engine-emission changes. Do not describe a successful
   emission build as closing or testing R4.

## Explicit deferral

R5 has the analogous linear-speed shape through `Q.brane.B_comp` and currently has no harness row. The S11
engines emit longitudinal squared-frequency roots but no longitudinal linear-speed object. R5/S11
emission, binding, and coverage are deferred from W0 and require their own reviewed work.
