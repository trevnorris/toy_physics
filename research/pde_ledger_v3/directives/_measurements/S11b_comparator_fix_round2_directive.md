# Measurements behind S11b_comparator_fix_round2_directive.md (rule 2)

## Defect 1 — same-name function+symbol falsely collides (hides a real difference as UNCOMPARED)
```
$ python3 -c "<load comparator>; compare fA(x)+fA vs g(x)+h"
fA(x)+fA vs g(x)+h  -> UNCOMPARED (should be DISAGREE)
  reason: $: RESIDUAL_FAILURE TRANSLITERATION_COLLISION ENGINE=PY TARGET=fA SOURCES=('FUNCTION:fA', 'SYMBOL:fA')
```

## Defect 2 — classification emitted only AFTER full-operand sstr (render hang; code)
L1268: print('\\n'.join(render_comparison(comparison)), flush=True) — the flush follows the full
render_comparison, which sstr-renders PY_PARSED_OBJECT/WL_PARSED_OBJECT. Grok measured
S11B_COMPRESSIONAL_RESPONSE: compare_records 8.3s -> DISAGREE, render_comparison >240s on the
120571/97753-char operands; the per-leaf residual budget does not enclose render.
