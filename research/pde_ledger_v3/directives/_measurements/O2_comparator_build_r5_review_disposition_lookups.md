# Measurements — O2 comparator build r5 review dispositions (generated 2026-10-07 21:41)

Generator: `_scratch/s9b_build/gen/o2_comparator_build_r5_lookups.sh` (sed/grep/sha256sum only; lines cut at 300 characters). Leg evidence quoted verbatim from the legs' filed stdout, copied to `_scratch/s9b_build/o2_comparator_build_review_r5_claude_evidence/` and `_scratch/s9b_build/o2_comparator_build_review_r5_grok_evidence/`.

```
$ sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r5.sha256
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/ablate_O2_cross_engine_comparator.py: OK
_scratch/s11c/o2-comparator-build-r5/comparison.jsonl: OK
_scratch/s11c/o2-comparator-build-r5/accounting.jsonl: OK
_scratch/s11c/o2-comparator-build-r5/catalog.json: OK
_scratch/s11c/o2-comparator-build-r5/validation.json: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r5_report.md: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r5_measurements.txt: OK
```

## Round-5 findings W1–W3 in the r5 artifact
W1: literal-match counts of the four WL provenance / held-operator names in the r5 production output.
```
$ for t in '"wl::6"' '"wl::3,8"' '"wl::D"' '"wl::Map"'; do printf '%s lines=' "$t"; grep -cF "$t" _scratch/s11c/o2-comparator-build-r5/comparison.jsonl; done
"wl::6" lines=0
"wl::3,8" lines=0
"wl::D" lines=14
"wl::Map" lines=14
```

W1 detail: where the remaining wl::D / wl::Map occur (compared inventory, difference pairs, raw Symbol leaves), and the census field that holds them.
```
$ cat _scratch/s9b_build/gen/o2_r5_w1_detail.sh
#!/bin/bash
# Literal-match counts: where wl::D and wl::Map remain in the r5 production comparison output.
f=/var/projects/toy_physics/_scratch/s11c/o2-comparator-build-r5/comparison.jsonl
for t in D Map; do
  printf 'named_OPEN_operands with wl::%s: ' $t; grep -c "\"named_OPEN_operands\":{[^}]*\"wl::$t\"" $f
  printf 'difference pairs ["wl::%s",n]: ' $t; grep -c "\[\"wl::$t\",-\?[0-9]*\]" $f
  printf 'Symbol leaves "Symbol","wl::%s": ' $t; grep -c "\"Symbol\",\"wl::$t\"" $f
  printf 'census "named_operands_and_heads" objects holding wl::%s: ' $t; grep -o "\"named_operands_and_heads\":{[^{}]*\"wl::$t\"" $f | wc -l
done
```

```
$ bash _scratch/s9b_build/gen/o2_r5_w1_detail.sh
named_OPEN_operands with wl::D: 0
difference pairs ["wl::D",n]: 0
Symbol leaves "Symbol","wl::D": 14
census "named_operands_and_heads" objects holding wl::D: 8
named_OPEN_operands with wl::Map: 0
difference pairs ["wl::Map",n]: 0
Symbol leaves "Symbol","wl::Map": 14
census "named_operands_and_heads" objects holding wl::Map: 8
```

```
$ sed -n '838,839p;859p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
def structure(n):
    names, heads, options, orientations = Counter(), Counter(), Counter(), Counter()
    return {'named_operands_and_heads':dict(sorted(names.items())),
```

```
$ grep -n 'raw_syntax_policy' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
996:    'raw_syntax_policy':'raw trees, censuses and serialization hashes are diagnostics, not additional OPEN comparisons',
1410:                         'serialization_scope':OPEN_SCOPE['raw_syntax_policy'],
```

W2: harness lines that name term_orientation.
```
$ grep -c 'term_orientation' research/pde_ledger_v3/scripts/ablate_O2_cross_engine_comparator.py
4
```

W3: the printed not_compared list.
```
$ grep -n -A6 "'not_compared'" research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py | head -10
992:    'not_compared':['where an object sits among OPEN arguments','how many times an object occurs',
993-                    'content outside these inventories in OPEN occurrences and balance entries: profile-free velocity or metric algebra; held derivative or variation variables; held aggregate sets and binders; scalar coefficients beyond their orientation sign; provenance, constructor options and
994-    'occurrence_indices':'local engine enumeration only; no cross-engine argument or occurrence slot pairing',
995-    'leaf_policy':'parsed terminal nodes and constructor options; consumed leaves construct a role/head wherever the action sits, a named operand/label, or a complete live-object key. Options, binder declarations outside a live key, provenance and other syntax are outside. Inactive operator spel
996-    'raw_syntax_policy':'raw trees, censuses and serialization hashes are diagnostics, not additional OPEN comparisons',
997-}
998-
```

## N1 — sign reading of a negative non-integer rational coefficient
```
$ grep -n "startswith('-')" research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
793:                if child.head == 'Number' and child.value.startswith('-'):
1151:    if n.head=='Number' and n.value.startswith('-'):
```

```
$ sed -n '1,14p' _scratch/s9b_build/o2_comparator_build_review_r5_claude_evidence/rational_probe/probe.stdout
CASE integer -1 both
   entry orientation py [-1]
   entry orientation wl [-1]
   balance orientation difference {"OPEN_free": [], "[\"OPEN_MomentumDensity_0\"]": []}
   action_comparison orientation [{"differences_orientation": [], "py": {"named_OPEN_operands": {"I_br_live": 1}, "live_arguments": {}, "head": {"OPEN_MomentumDensity_0": 1}, "orientation": {"-1": 1}, "role": {"OPEN_MomentumDensity_0": 1}}, "wl": {"named_OPEN_operands": {"I_br_live": 1}, "live_argum
CASE Rational(-1,2) vs -(1/2) [same sign]
   entry orientation py [1]
   entry orientation wl [-1]
   balance orientation difference {"OPEN_free": [], "[\"OPEN_MomentumDensity_0\"]": [["-1", -1], ["1", 1]]}
   action_comparison orientation [{"differences_orientation": [["-1", -1], ["1", 1]], "py": {"named_OPEN_operands": {"I_br_live": 1}, "live_arguments": {}, "head": {"OPEN_MomentumDensity_0": 1}, "orientation": {"1": 1}, "role": {"OPEN_MomentumDensity_0": 1}}, "wl": {"named_OPEN_operands": {"I_br_liv
CASE Rational(-1,2) vs +(1/2) [opposite sign]
   entry orientation py [1]
   entry orientation wl [1]
   balance orientation difference {"OPEN_free": [], "[\"OPEN_MomentumDensity_0\"]": []}
```

```
$ grep -h 'terms_with_Rational_on_sign_path\|differing' _scratch/s9b_build/o2_comparator_build_review_r5_claude_evidence/sign_path_probe.stdout _scratch/s9b_build/o2_comparator_build_review_r5_claude_evidence/rational_probe/real_stream_exposure.stdout | head -12
hold_inplane py components 3 additive_terms 24 terms_with_Rational_on_sign_path 0
hold_inplane wl components 3 additive_terms 24 terms_with_Rational_on_sign_path 0
hold_bulk py components 1 additive_terms 8 terms_with_Rational_on_sign_path 0
hold_bulk wl components 1 additive_terms 8 terms_with_Rational_on_sign_path 0
hold_normal py components 1 additive_terms 32 terms_with_Rational_on_sign_path 0
hold_normal wl components 1 additive_terms 32 terms_with_Rational_on_sign_path 0
energy_balance py components 1 additive_terms 5 terms_with_Rational_on_sign_path 0
energy_balance wl components 1 additive_terms 5 terms_with_Rational_on_sign_path 0
join rows with differing action_comparison: 0
balance rows with differing printed output: 0
```

```
$ grep -h 'sign_mismatches\|negative_rational_terms' _scratch/s9b_build/o2_comparator_build_review_r5_grok_evidence/recompute_stdout.txt | head -4
BALANCE sign_mismatches=0 negative_rational_terms=0
```

## Measurement for the record — WL-only xi_w second derivative in the in-plane current entries
```
$ head -12 _scratch/s9b_build/o2_comparator_build_review_r5_claude_evidence/xi2/wl_inplane_terms.stdout
component head Plus length 8
TERM 1 head=Times roles={InternalForce[1]} XiW orders {{1, 12}} n>=2 positions 0
TERM 2 head=OpenAction roles={OutwardAdditionalMomentumPartner[1]} XiW orders {{1, 12}} n>=2 positions 0
TERM 3 head=OpenFirstVariation roles={MomentumCurrent[1, 1]} XiW orders {{1, 43}, {2, 12}} n>=2 positions 12
    head chain to first: {List, Identity, Identity, Identity, Identity, Identity}
    enclosing expr: (x1^3*VR[Sqrt[x1^2 + x2^2 + x3^2]]*Derivative[2][XiW][Sqrt[x1^2 + x2^2 + x3^2]])/(x1^2 + x2^2 + x3^2)^(3/2)
TERM 4 head=OpenFirstVariation roles={MomentumCurrent[1, 2]} XiW orders {{1, 43}, {2, 12}} n>=2 positions 12
    head chain to first: {List, Identity, Identity, Identity, Identity, Identity}
    enclosing expr: (x1^2*x2*VR[Sqrt[x1^2 + x2^2 + x3^2]]*Derivative[2][XiW][Sqrt[x1^2 + x2^2 + x3^2]])/(x1^2 + x2^2 + x3^2)^(3/2)
TERM 5 head=OpenFirstVariation roles={MomentumCurrent[1, 3]} XiW orders {{1, 43}, {2, 12}} n>=2 positions 12
    head chain to first: {List, Identity, Identity, Identity, Identity, Identity}
    enclosing expr: (x1^2*x3*VR[Sqrt[x1^2 + x2^2 + x3^2]]*Derivative[2][XiW][Sqrt[x1^2 + x2^2 + x3^2]])/(x1^2 + x2^2 + x3^2)^(3/2)
```

```
$ head -12 _scratch/s9b_build/o2_comparator_build_review_r5_claude_evidence/xi2/py_flux_orders.stdout
PY_O2_MOMENTUM_FLUX occurrences= 12
    OPEN_MomentumFlux_0_0 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_0_1 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_0_2 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_1_0 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_1_1 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_1_2 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_2_0 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_2_1 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_2_2 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_3_0 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
    OPEN_MomentumFlux_3_1 xi_w derivative orders {1: 105} xi_w undifferentiated (non-dummy arg) 4
```

