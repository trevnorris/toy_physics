# Measurements behind S11b_comparator_fix_round1_directive.md (rule 2)
# The three baseline holes, each a genuine difference reported as non-DISAGREE:
```
$ python3 research/pde_ledger_v3/directives/_measurements/s11b_comparator_holes_probe.py
D1 function-head collision: AGREE (reason=None)
D2 status buries sibling: UNDECIDED (reason=ENGINE_STATUS PY=STATUS_TOKEN=UNDECIDED WL=STATUS_TOKEN=UNDECIDED)
D3 symbol-pair tuple promoted: AGREE (reason=None)
```
D1: function heads f_a/fA collapse under transliteration -> AGREE (should be UNCOMPARED collision).
D2: STATUS_TOKEN present -> UNDECIDED before residualing RESULT a vs b (should be DISAGREE).
D3: Symbol-keyed pair-tuple promoted to Association -> AGREE (should be STRUCTURE DISAGREE).
