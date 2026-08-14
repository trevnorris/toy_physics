# Measurements — S11_wl_fix_round1_brief.md

Every factual claim in the brief, with the command that produced it and its literal output. Run from
`/var/projects/toy_physics` at `46ba77c2` (the committed as-reviewed engine; the working-tree `.wl` is
byte-identical to that blob). Leg evidence preserved under `~/.s11_build/wl/evidence_r1/` (copied from
the fresh-agent leg's scratch `/tmp/s11wls_leg_tIDD/`); leg reports at `~/.s11_build/wl/`.

This brief was folded once after two review legs on it (Codex `fix1_leg_codex.log`, Grok
`fix1_leg_grok.log`); the pre-fold defects they found — an acceptance green path for a non-deciding
fix, an expected-value regression anchor, an undeclared point-failure payload demand, a reachability
overreach, dependency-only walker probes, and the promotion that became item 4 — are recorded in those
logs, not restated here.

## Item 1 — `RationalQ` undefined, not a builtin, and the emitted payloads are inert

```
$ rg -c "RationalQ" research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
1

$ sed -n '905,910p' research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
        DimensionalProduct[dimensions]],
    Head[expression] === Power &&
        (IntegerQ[expression[[2]]] || RationalQ[expression[[2]]]),
      baseDimension = dimensionOfScalar[expression[[1]], atomDimensions, assumptions];
      If[ListQ[baseDimension], expression[[2]] baseDimension,
        DimensionalPower[baseDimension, expression[[2]]]],

$ timeout 60 /usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel -noprompt \
    -run 'Print[{IntegerQ[1/2], Head[1/2] === Rational, RationalQ[1/2]}]; Quit[]'
{False, True, RationalQ[1/2]}
```

The unevaluated echo in third position is the kernel's own statement that `RationalQ` is not a builtin.
(A `wolframscript -code` form of this probe works in the orchestrator's environment but exits 255 in
the Codex sandbox; the direct-kernel command above reproduces everywhere and is the twin's probe.)

Captured live consequence (leg capture, preserved; truncation is part of the pipe):

```
$ grep -n "RationalQ" ~/.s11_build/wl/evidence_r1/xkin2.out | head -2 | cut -c1-200
208:WL_S11_XKIN_ANISO_D2_ROOT1_DIM_OVER_KSQ: DimensionalProduct[{{0, 0, 0}, {2, 0, 0}, {2*dimensionSymbol, 0, -2}, {0, 0, 0}, DimensionalAlternatives[{-2*dimensionSymbol, -2, 2}, DimensionalProduct[{{
209:WL_S11_XKIN_ANISO_D2_ROOT2_DIM_OVER_KSQ: DimensionalProduct[{{0, 0, 0}, {2, 0, 0}, {2*dimensionSymbol, 0, -2}, {0, 0, 0}, DimensionalAlternatives[{-2*dimensionSymbol, -2, 2}, Which[RationalQ[1/2],
```

## Item 2 — status by structure, literal certificate

```
$ sed -n '641,648p' research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
emitComponentCount[package_, dimension_, baseQuantity_String, value_,
    freeParameters_List] := Module[{status, certificate},
  status = If[freeParameters === {}, tokenConstant, tokenUndecided];
  certificate = If[SameQ[status, tokenConstant], {
      "FREE_PARAMETERS" -> freeParameters,
      "TEST_OBJECT" -> True,
      "OPERANDS" -> {"VALUE" -> value, "FREE_PARAMETERS" -> freeParameters}
    }, tokenNotApplicable];

$ grep -c "CHANGE_LOCUS" ~/.s11_build/wl/evidence_r1/main2_single.out
72
$ grep "CHANGE_LOCUS" ~/.s11_build/wl/evidence_r1/main2_single.out | grep -vc "NOT_APPLICABLE"
0
```

"No token excluded by the shape of the code": the assignment has two arms for three declared tokens.

## Item 3 — value-conditional stratum drop; no declared failure payload

```
$ sed -n '1151,1155p' research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
        component["Variables"]]]],
    mergeCandidates[candidates]];
  strata = Select[strata, ListQ[#["Point"]] &];
  ordering = Lookup[strata, "Branch", {}];
  emitCell[package, dimension, "STRATUM_ORDERING", ordering];

$ rg -n "POINT_OUTCOME|POINT.*FAIL" research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md; echo "exit=$?"
exit=1
```

Latent: no captured stream shows the drop firing (every admissible component produced a point). The
spec declares an exact point and nothing else — hence the item's two honest outcomes.

## Item 4 — a solver non-answer becomes PROVED_FALSE

```
$ sed -n '91,95p;198,207p' research/pde_ledger_v3/mathematica/S11_stray_longitudinal_mathematica_audit.wl
truthStatus[testObject_] := Which[
  TrueQ[testObject], tokenProvedTrue,
  SameQ[testObject, False], tokenProvedFalse,
  True, tokenUndecided
];
  unrestrictedReduction = Quiet[Reduce[And @@ equations, variables, Complexes]];
  inconsistentOperands = {"EQUATIONS" -> equations,
    "SOLVE_VARIABLES" -> variables,
    "UNRESTRICTED_REDUCTION" -> unrestrictedReduction};
  inconsistentTest = SameQ[unrestrictedReduction, False];
  emitCell[package, dimension, baseQuantity <> "_INCONSISTENT", {
    "STATUS_TOKEN" -> truthStatus[inconsistentTest],
    "TEST_OBJECT" -> inconsistentTest,
    "OPERANDS" -> inconsistentOperands
  }];

$ timeout 60 /usr/local/Wolfram/Wolfram/15.0/Executables/WolframKernel -noprompt \
    -run 'u=Inactive[Reduce][x==0,x,Complexes]; Print[{u, SameQ[u,False]}]; Quit[]'
{Inactive[Reduce][x == 0, x, Complexes], False}
```

`SameQ[<unevaluated reduction>, False]` is the decided boolean `False`, which `truthStatus` maps to
the proved token — the control demonstrates the path without any engine run.

The branch is live, not hypothetical — the engine's own captured stream reaches this site with radical
systems:

```
$ awk 'BEGIN{n=0;s=0} /_INCONSISTENT:/{n++; if(index($0,"Sqrt[")>0)s++} END{print "inconsistent_tags=" n; print "with_sqrt=" s}' ~/.s11_build/wl/evidence_r1/xkin2.out
inconsistent_tags=79
with_sqrt=14
```

Spec's pinned payload form:

```
$ sed -n '253,256p' research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
⭐⭐ **THE BOOLEAN-TEST PAYLOAD FORM IS PINNED.** Each `_IDENTICALLY_SATISFIED` and `_INCONSISTENT`
payload has this ordered field sequence:

```
```

## Acceptance item 2's feasibility claim

The `XKIN_ANISO` D2 cell emits the `_DIM_` family early enough to observe within budget — the leg's
capture reached those tags at stream lines 208–209 (grep above), and its runner record is preserved:

```
$ cat ~/.s11_build/wl/evidence_r1/xkin2.out.runner
exit=137 elapsed=267s killed=RSS_OVER_6GB rss_kb=6299572
```

## Not measurements — decisions, marked as such

- Which findings became items, and the promotion of item 4 from a below-filter footnote after a leg
  showed its unreachability premise false: orchestrator adjudication (rule 13; every number above
  re-measured by the orchestrator, not transcribed).
- Acceptance shapes (identity residuals; parse-as-declared-object; per-family certificate ablation;
  live corrupted-extractor diff; mapping unit probe): orchestrator decisions incorporating both legs'
  counterexample constructions.
- Regression-scope review deliberately moved from builder-facing acceptance to the script legs: the
  committed stream must not be a builder target.
