# Measurements — O2 comparator build r3 review dispositions (generated 2026-10-07 18:25)

Generator: `_scratch/s9b_build/gen/o2_comparator_build_r3_lookups.sh` (sed/grep/sha256sum only; lines cut at 300 characters).

```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r3.sha256)
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/ablate_O2_cross_engine_comparator.py: OK
_scratch/s11c/o2-comparator-build-r3/comparison.jsonl: OK
_scratch/s11c/o2-comparator-build-r3/accounting.jsonl: OK
_scratch/s11c/o2-comparator-build-r3/catalog.json: OK
_scratch/s11c/o2-comparator-build-r3/validation.json: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r3_report.md: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r3_measurements.txt: OK
```

## V1 — the coupled_embedding join, its two engine constructions, and the count-token reason
```
$ grep -n "coupled_embedding" scripts/O2_cross_engine_comparator.py | head -4
407:    J('coupled_embedding','COUPLED_INPUTS/embedding','B_HOLD_LIVE/O4Identity',395,239,'§3.2'),
```

```
$ grep -n "O4EquationIdentityCount" scripts/O2_cross_engine_comparator.py
1265:        'COUPLED_INPUTS_MODEL_POINT/O4EquationIdentityCount':'Additional WL unsettled-count token; PY carries it within the joined coupled-embedding object, with no separate count token (§3.2).',
```

```
$ sed -n '395,396p' scripts/O2_live_balance_sympy_audit.py
            embedding=sp.Tuple(op['E_h_live'], sp.Eq(xi, ell*h, evaluate=False),
                               Str('normal_equation_identity_and_count_unsettled')),
```

```
$ sed -n '22p;56p;235p;239p;340p' mathematica/O2_live_balance_mathematica_audit.wl
operand[name_, section_, contract_] := OPEN[name, origin[section, contract]];
embeddingRelation = operand[EhLive, "3.2", "5"];
holdBalance = account[momentumEntries];
  "Status" -> ConditionalNamedBalance, "O4Identity" -> UnresolvedIdentification[embeddingRelation,holdBalance],
  "O4EquationIdentityCount" -> Unsettled,
```

## V2 — the OPEN-occurrence feature inventory: presence sets, Name heads, Str→Text
```
$ sed -n '984,1010p' scripts/O2_cross_engine_comparator.py
def features(n,bindings,engine):
    """Inventories of objects, not counts of engine-specific syntax nodes."""
    mapped = bindings.apply(n,engine)
    names,live,binders = Counter(),Counter(),Counter()
    def visit(v,head=False):
        if v.head=='Name' and not head and bindings.kinds.get(v.value) not in ('coordinate','parameter','profile','binder'):
            names[v.value] = 1
        if v.head in ('LiveProfile','ProfileDerivative'):
            live[object_key(v)] = 1
        if v.head in ('ProfileDerivative','Derivative','Subs','Lambda','FirstVariation','HeldDerivative'):
            binders[object_key(v)] = 1
        for i,x in enumerate(v.args):
            visit(x,head=v.head=='Apply' and i==0)
    visit(feature_tree(n,bindings,engine))
    return {'head':Counter({json.dumps([mapped.head,mapped.value,
                data(mapped.args[0]) if mapped.head=='Apply' else None],separators=(',',':')):1}),
            'named_OPEN_operands':names,'live_arguments':live,'binders':binders,
            'argument_trees':Counter({fingerprint(mapped):1})}


def merge_fields(current,fields):
    for field,values in fields.items():
        inventory=current.setdefault(field,Counter())
        if field in OBJECT_FIELDS:
            inventory.update({key:1 for key in values if key not in inventory})
        else:
            inventory.update(values)
```

```
$ grep -n "'Text'" scripts/O2_cross_engine_comparator.py | head -6
65:                        'Text' if isinstance(n.value, str) else 'Number', n.value)
95:            if len(a) != 1 or a[0].head != 'Text':
97:            return Node({'Symbol': 'Symbol', 'Dummy': 'Dummy', 'Str': 'Text',
105:                         and v.args[0].head == 'Text' for v in a):
171:            if any(e.head != 'Rule' or e.args[0].head != 'Text' for e in entries):
175:            left = atom('Text', ast.literal_eval(tok))
```

## V3 — the Sqrt translation and the harness knife that targets it
```
$ sed -n '645,650p' scripts/O2_cross_engine_comparator.py
        if head.head in ('FunctionName','Symbol'):
            key = bindings.name(head.value, engine)
            if engine == 'wl' and head.value == 'Sqrt' and len(a) == 2:
                return sp.sqrt(algebra(a[1], bindings, engine))
            if bindings.kinds.get(key) == 'profile':
                if len(a) < 2:
```

```
$ grep -n 'canonical_sqrt_removed' scripts/ablate_O2_cross_engine_comparator.py scripts/O2_cross_engine_comparator.py | head -6
scripts/ablate_O2_cross_engine_comparator.py:58:    ('canonical_sqrt_removed',
```

