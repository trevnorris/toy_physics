# Measurements behind S11b_wl_engine_fix_round2_directive.md (rule 2)
# The two dead-control / chi5-tuning defects (found by both WL-repair legs), at these lines:
```
1237:unrelatedCausalityAggregation = aggregateNamedTestRecords[
1255:  "UNRELATED_OBJECT_CONTROL" -> unrelatedCausalityAggregation,
1261:      unrelatedCausalityAggregation["AGGREGATE_TEST_OBJECT"] &&
1608:    stratumRule = First[Solve[dispersion == 0, chi5]];
1614:      "CLASSIFICATION" -> classifyLeadingBlock[stratumBlock]|>
1648:thresholdUnrelatedClassification = thresholdLeadingClassification;
```
F-WL-3c: CLASSIFICATION -> classifyLeadingBlock[stratumBlock], stratumBlock = block /. First[Solve[
  dispersion==0, chi5]] -> chi5 tuned to force Det=0 -> tautological rank; thresholdUnrelatedClassification
  = thresholdLeadingClassification (alias) -> A-A unrelated control (residuals {0,0} by construction).
F-WL-3b: unrelatedCausalityAggregation re-aggregates the UNMODIFIED kernelOrientationIdentities -> A-A.
