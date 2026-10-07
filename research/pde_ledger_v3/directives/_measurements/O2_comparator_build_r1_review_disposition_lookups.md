# Measurements — O2 comparator build r1 review dispositions (generated 2026-10-07 16:05)

Generator: `_scratch/s9b_build/gen/o2_comparator_build_r1_lookups.sh` (sed/grep/sha256sum only; lines cut at 300 characters).

```
$ (cd /var/projects/toy_physics && sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r1.sha256)
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
_scratch/s11c/o2-comparator-build-r1/comparison.jsonl: OK
_scratch/s11c/o2-comparator-build-r1/accounting.jsonl: OK
_scratch/s11c/o2-comparator-build-r1/catalog.json: OK
_scratch/s11c/o2-comparator-build-r1/summary-complete.json: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r1_report.md: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r1_measurements.txt: OK
```

## S1 — orientation census: which heads get an orientation; the OPEN-free label; the carried entry in both engines
```
$ grep -n 'OPEN-free' scripts/O2_cross_engine_comparator.py
```

```
$ grep -n 'def signed_actions' -A30 scripts/O2_cross_engine_comparator.py
679:    def signed_actions(v,sign=1,dependency=False):
680-        # Linear syntax carries the enclosing expression's orientation. This is
681-        # the engines' held calculus/aggregation syntax, not a constitutive law:
682-        # PY audit 134-140, 242-244; WL audit 142-143, 181-193.
683-        if v.head == 'Mul':
684-            for child in v.args:
685-                if child.head == 'Number' and child.value.startswith('-'):
686-                    sign *= -1
687-            for child in v.args:
688-                if child.head != 'Number':
689-                    signed_actions(child,sign,dependency)
690-        elif v.head == 'Apply':
691-            head = v.args[0]
692-            linear_arguments = ()
693-            if (head.head == 'Apply' and len(head.args) == 2
694-                    and head.args[0].value == 'wl::Inactive'):
695-                operator = head.args[1].value
696-                if operator in ('wl::Total','wl::Map','wl::D'):
697-                    linear_arguments = (1,)
698-            elif head.value in ('wl::Function','wl::OpenFirstVariation'):
699-                linear_arguments = (2,) if head.value == 'wl::Function' else (1,)
700-            elif head.value == 'py::OPEN_SumOverAllNativeFaces':
701-                linear_arguments = (2,)
702-            role = None
703-            if head.value.startswith('py::OPEN_'):
704-                role = head.value
705-            elif head.value == 'wl::OpenAction' and len(v.args) > 1:
706-                role = json.dumps(data(v.args[1]),separators=(',',':'))
707-            if role is not None:
708-                if dependency:
709-                    argument_actions[role] += 1
```

```
$ sed -n '259,260p' scripts/O2_live_balance_sympy_audit.py
    entries = ((1, storage), (1, transport), (1, carry), (1, partners),
               (-1, internal), (-1, reduced_load)) + body_entries
```

```
$ grep -n 'CarriedExchange' mathematica/O2_live_balance_mathematica_audit.wl
232:  entry[CarriedExchange,1,carriedMomentum,origin["5","1,7,10"]],
```

## S2 — derivative evaluation point: the comparator branches and the stripped-argument control
```
$ grep -n 'Subs\|Derivative' scripts/O2_cross_engine_comparator.py | head -30
597:        # InputForm Derivative[n][F][argument] keeps the applied argument live.
599:                and head.args[0].args[0] == atom('Symbol','Derivative')):
605:                derivative = sp.Derivative(sp.Function(key)(var),(var,int(orders[0].value)),evaluate=False)
606:                return sp.Subs(derivative,var,algebra(a[1],bindings,engine)).doit()
608:    if h == 'Derivative':
616:        return sp.Derivative(v,*specs,evaluate=False)
617:    if h == 'Subs' and len(a) == 3:
637:        return sp.Subs(algebra(dummy(a[0]),local,engine),
720:        elif v.head in ('Derivative','Lambda'):
721:            body = 0 if v.head == 'Derivative' else 1
882:    if engine == 'wl' and path[:2] == ('MATERIAL_MOMENTUM','MaterialProfileDerivatives'):
```

```
$ grep -n 'def test_applied_argument_stripped' -A6 scripts/test_O2_cross_engine_comparator.py
97:    def test_applied_argument_stripped(self):
98-        self.changes(fn('V_r','Symbol("x1")'),'VR[x1]','VR')
99-
100-    def test_live_profile_frozen(self):
101-        self.changes(fn('o2_rho_br_live','Symbol("x1")'),'RhoBr[x1]','11')
102-
103-    def test_open_head(self):
```

```
$ grep -c 'Subs' scripts/test_O2_cross_engine_comparator.py
2
```

## S3 — orientation through a held SymPy Derivative: controls
```
$ grep -n 'Derivative' scripts/test_O2_cross_engine_comparator.py
124:        py = 'Subs(Derivative(Function("V_r")(Dummy("z")), Tuple(Dummy("z"), Integer(1))), Tuple(Dummy("z")), Tuple(Symbol("x1")))'
125:        self.changes(py,'Derivative[1][VR][x1]','Derivative[2][VR][x1]')
209:        py = 'Subs(Derivative(Function("V_r")(Dummy("z")), Tuple(Dummy("z"), Integer(1))), Tuple(Dummy("z")), Tuple(Symbol("x1")))'
210:        result = self.fixture(py,'Derivative[1][VR][x1]')[0][0]
```

## S4 — item-4 differences and the constant argument field
```
$ grep -n 'live_arguments_and_binders' scripts/O2_cross_engine_comparator.py
739:            'live_arguments_and_binders':'complete mapped operand tree; no argument erasure'}
```

```
$ grep -n 'def differences' -A25 scripts/O2_cross_engine_comparator.py
651:def differences(a, b, path=()):
652-    """Tree differences retain mismatched subtrees by reference to printed operands.
653-
654-    Descend through paired children even when their parent heads differ. Unpaired
655-    subtrees are emitted in full so mutations below a mismatch remain observable.
656-    """
657-    if (a.head,a.value,a.options,len(a.args)) != (b.head,b.value,b.options,len(b.args)):
658-        yield {'path':path,'py_node':[a.head,a.value,a.options,len(a.args)],
659-               'wl_node':[b.head,b.value,b.options,len(b.args)],
660-               'subtrees':'see mapped operands at this path'}
661-    for i,(x,y) in enumerate(zip(a.args,b.args)):
662-        yield from differences(x,y,(*path,i))
663-    for i in range(min(len(a.args),len(b.args)),max(len(a.args),len(b.args))):
664-        yield {'path':(*path,i),'py_unpaired':data(a.args[i]) if i < len(a.args) else None,
665-               'wl_unpaired':data(b.args[i]) if i < len(b.args) else None}
666-
667-
668-def structure(n):
669-    names, heads, options, orientations = Counter(), Counter(), Counter(), Counter()
670-    argument_actions = Counter()
671-    def visit(v):
672-        heads[v.head] += 1
673-        if v.head in ('Symbol','FunctionName','Dummy'):
674-            names[v.value] += 1
675-        options.update(v.options)
676-        for child in v.args:
```

