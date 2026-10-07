# Measurements — O2 comparator build r2 review dispositions (generated 2026-10-07 17:16)

Generator: `_scratch/s9b_build/gen/o2_comparator_build_r2_lookups.sh` (sed/grep/sha256sum only; lines cut at 300 characters).

```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r2.sha256)
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
_scratch/s11c/o2-comparator-build-r2/comparison.jsonl: OK
_scratch/s11c/o2-comparator-build-r2/accounting.jsonl: OK
_scratch/s11c/o2-comparator-build-r2/catalog.json: OK
_scratch/s11c/o2-comparator-build-r2/summary-complete.json: OK
_scratch/s11c/o2-comparator-build-r2/action-coverage.json: OK
_scratch/s11c/o2-comparator-build-r2/balance-evidence.json: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r2_report.md: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r2_measurements.txt: OK
```

## T1 — the OPEN-action feature keys, the Derivative clause, and the tests that read them
```
$ sed -n '863,886p' scripts/O2_cross_engine_comparator.py
def features(n,bindings,engine):
    """Computed multisets; repeated occurrences are counted, never zipped."""
    mapped = bindings.apply(n,engine)
    names,live,binders = Counter(),Counter(),Counter()
    def visit(v):
        if v.head in ('Symbol','Dummy') and bindings.kinds.get(v.value) not in ('coordinate','parameter','profile','binder'):
            names[v.value] += 1
        if v.head=='Apply' and (bindings.kinds.get(v.args[0].value)=='profile' or
                (v.args[0].head=='FunctionName' and not v.args[0].value.startswith('py::OPEN_'))):
            live[json.dumps(data(v),separators=(',',':'))] += 1
        if v.head in ('Derivative','Subs','Lambda') or (v.head=='Apply' and
                (v.args[0].value in ('wl::Function','wl::OpenFirstVariation') or
                 'wl::Derivative' in {x.value for x in v.args[0].args})):
            binders[fingerprint(v)] += 1
        arguments=v.args
        if v.head=='Apply':
            arguments=v.args[2:] if v.args[0].value=='wl::OpenAction' else v.args[1:]
        for x in arguments:
            visit(x)
    visit(mapped)
    return {'head':Counter({json.dumps([mapped.head,mapped.value,
                data(mapped.args[0]) if mapped.head=='Apply' else None],separators=(',',':')):1}),
            'named_OPEN_operands':names,'live_arguments':live,'binders':binders,
            'argument_trees':Counter({fingerprint(mapped):1})}
```

```
$ grep -n "wl::Derivative" scripts/O2_cross_engine_comparator.py
875:                 'wl::Derivative' in {x.value for x in v.args[0].args})):
```

```
$ sed -n '430,460p' scripts/test_O2_cross_engine_comparator.py
    def test_action_field_differences(self):
        py=fn('OPEN_MomentumDensity_0',sy('I_br_live')+', '+fn('V_r',sy('x1')))
        wl='OpenAction[MomentumDensity[1],{IBrLive},VR[x1]]'
        base=self.fixture(py,wl)
        for field,mutation in (
            ('role',wl.replace('MomentumDensity[1]','InternalForce[1]')),
            ('head',wl.replace('MomentumDensity[1]','InternalForce[1]')),
            ('named_OPEN_operands',wl.replace('IBrLive','NBrLive')),
            ('live_arguments',wl.replace('VR[x1]','VR[x2]')),
            ('orientation','-('+wl+')'),
            ('argument_trees',wl.replace('VR[x1]','Wrapper[VR[x1]]'))):
            with self.subTest(field=field):
                self.assertNotEqual(self.action_deltas(base,field),
                                    self.action_deltas(self.fixture(py,mutation),field))

    def test_balance_entry_field_differences(self):
        py='Add('+fn('OPEN_MomentumDensity_0',sy('I_br_live')+', '+fn('V_r',sy('x1')))+',Integer(3))'
        wl='OpenAction[MomentumDensity[1],{IBrLive},VR[x1]]+3'
        def delta(text,field):
            product=self.balance_products(self.fixture(py,text,balance=True))[0]
            return {role:fields.get(field,[]) for role,fields in product['entry_differences'].items()}
        for field,changed in (
            ('role',wl.replace('MomentumDensity[1]','InternalForce[1]')),
            ('head',wl.replace('MomentumDensity[1]','InternalForce[1]')),
            ('named_OPEN_operands',wl.replace('IBrLive','NBrLive')),
            ('live_arguments',wl.replace('VR[x1]','VR[x2]')),
            ('orientation','-('+wl+')'),
            ('argument_trees',wl.replace('VR[x1]','Wrapper[VR[x1]]'))):
            with self.subTest(field=field):
                self.assertNotEqual(delta(wl,field),delta(changed,field))

```

## T2 — the relation branch and the closed-part extraction
```
$ sed -n '1050,1062p' scripts/O2_cross_engine_comparator.py
        return {'outcome':'not_formed','reason':'text, native boolean or name is not a subtractable operand'},(0,0)
    # Relational operands retain the relation head, with each side compared.
    relations = {'Equality','Unequality','StrictGreaterThan','StrictLessThan','GreaterThan','LessThan'}
    if a.head in relations or b.head in relations:
        if a.head != b.head or len(a.args) != len(b.args):
            return {'outcome':'not_formed','reason':'relational structure differs'},(0,0)
        result, counts = [], [0,0]
        for x,y in zip(a.args,b.args):
            r,c = compare(x,y,bindings)
            result.append(r)
            counts = [counts[i]+c[i] for i in range(2)]
        return {'outcome':'relation','head':a.head,'operands':result},tuple(counts)
    if a.head == b.head == 'Symbol':
```

```
$ sed -n '985,996p' scripts/O2_cross_engine_comparator.py
            entries.append({'role':role,'orientation':term_orientation(term),
                            'open_free':not is_open,'operand':data(mapped)})
            current=group.setdefault(role,{})
            for field,values in fields.items():
                current.setdefault(field,Counter()).update(values)
        sides[engine]=entries
        groups[engine]=group
        closed[engine]=Node('Add',tuple(closed_terms)) if closed_terms else atom('Number',0)
    result={'outcome':'balance','entries':sides,'entry_differences':{
        role:field_deltas(groups['py'].get(role,{}),groups['wl'].get(role,{}))
        for role in sorted(groups['py'].keys() | groups['wl'].keys())},
        'closed_operands':{e:data(n) for e,n in closed.items()}}
```

```
$ grep -n 'rel1\|right operand\|relation' scripts/test_O2_cross_engine_comparator.py | head -12
```

