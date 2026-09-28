#!/usr/bin/env python3
"""First-grade end lifts and native whole-line forcing, from accepted returns.

No producer imports, root/mode solver, inverse reconstruction, numerical
response, integral evaluation or automatic continuation. A fresh method gate
and the shared whole-job guard are required before importing SymPy.
"""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import pickle
import re
import resource
import signal
import time
import traceback

ROOT = Path('/var/projects/toy_physics')
STORE = ROOT / '_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
GRADES = ((0, 0), (1, 0), (0, 1))


def require(test, message):
    if not test:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def route(path):
    path = Path(path)
    return {'path': str(path), 'canonicalPath': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'sha256': digest(path)}


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def containment():
    group = next(s[3:] for s in Path('/proc/self/cgroup').read_text().splitlines()
                 if s.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    result = {k: (root/k).read_text().strip()
              for k in ('memory.max', 'memory.swap.max', 'pids.max')}
    result.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={k: os.environ.get(k) for k in THREADS})
    require(result['memory.max'] == str(2*1024**3) and result['memory.swap.max'] == '0'
            and result['pids.max'] == '32' and result['nice'] >= 15
            and len(result['affinity']) == 1
            and all(v == '1' for v in result['threads'].values()),
            'required whole-job containment absent; no science imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise TimeoutError('840-second native limit; preserve incomplete input')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(840)
    return {**result, 'nativeWallSeconds': 840}


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.'))
                or module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError(('unapproved saved class', module, name))
        return super().find_class(module, name)


class Journal:
    def __init__(self, out):
        self.out, self.records = out, []
        self.active = None

    def value(self, name, value):
        path = self.out/name; path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            pickle.dump(value, stream, protocol=4)
            stream.flush(); os.fsync(stream.fileno())
        return route(path)

    def op(self, name, function, *args):
        folder = 'operations/%04d-%s' % (len(self.records), name)
        record = {'name': name, 'startedUtc': datetime.now(timezone.utc).isoformat(),
                  'input': self.value(folder+'/input.pickle', args)}
        save(self.out/folder/'started.json', record)
        self.active = record
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle', result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json', record)
        self.records.append(record)
        self.active = None
        return result


def construct(spec, out, journal):
    import sympy as sp
    from sympy.core.function import AppliedUndef
    j = journal
    actual = {name: route(r['path']) for name, r in spec['inputs'].items()}
    save(out/'prehashes.json', actual)
    require(actual == spec['inputs'], 'pinned source/input bytes changed')
    require(actual['reduction']['sha256'] == actual['ownReferenceReduction']['sha256'],
            'native and accepted own-reference reduction bytes differ')
    cache = {}
    def saved(name):
        require(name in actual, ('unmanifested input', name))
        if name not in cache:
            def read(path):
                with Path(path).open('rb') as stream:
                    return SavedCodec(stream).load()
            cache[name] = j.op('restore-'+name, read, actual[name]['path'])
        return cache[name]
    def text_json(name):
        return json.loads(Path(actual[name]['path']).read_text())
    def simplify(value):
        # Local scalar algebra only. Native integrals are never sent here.
        require(not value.has(sp.Integral, sp.Limit), 'algebraic reducer received an action integral')
        return sp.cancel(sp.simplify(value))
    def matrix_reduce(value):
        return sp.ImmutableMatrix(value).applyfunc(simplify)
    def zero_check(name, value):
        result = j.op(name, matrix_reduce, value)
        save(out/(name+'.json'), {'shape': list(result.shape),
             'zero': result == sp.zeros(*result.shape),
             'unresolvedEntries': [[i, k] for i in range(result.rows)
                                  for k in range(result.cols) if result[i, k] != 0]})
        require(result == sp.zeros(*result.shape), ('exact residual unresolved or nonzero', name))
        return result

    cp = text_json('candidateCheckpoint')
    require(cp['status'] == 'VALIDATED_FIXED_INPUT_SEPARATED_POINT_OUTGOING_PRESCRIPTION_CANDIDATE',
            'accepted candidate checkpoint required')
    production = text_json('productionCheckpoint')
    require(cp['productionCheckpoint']['sha256'] == actual['productionCheckpoint']['sha256']
            and cp['validationCheckpoint']['sha256'] == actual['validationCheckpoint']['sha256'],
            'accepted checkpoint chain identity')
    require(cp['primaryArtifacts']['operation-index.json']['sha256'] == actual['candidateOperationIndex']['sha256']
            and cp['primaryArtifacts']['prescription-domain.pickle']['sha256'] == actual['domain']['sha256'],
            'accepted operation inventory and domain identities')
    for name, relative in spec['candidateArtifacts'].items():
        require(production['artifactSha256'][relative] == actual[name]['sha256'],
                ('accepted candidate route changed', name))
    candidate = saved('candidate')
    context, units = saved('context'), saved('units')
    uniform, reduction = saved('uniform'), saved('reduction')
    action, assembly_packet, grades = saved('actions'), saved('assembly'), saved('grades')
    assembly = assembly_packet['result']; state = reduction['reductionState']
    physical = text_json('physicalInput')
    require(candidate['fixedPhysicalInput'] == physical, 'same fixed physical input')
    require(candidate['effectiveReferenceGrades'] == {'eta_bg': 0, 'sigma_W': 0},
            'actual reference grades must remain zero')
    require(candidate['pointwiseDomain'] == sp.Ne(state['z'], state['zp'], evaluate=False),
            'accepted separated-point domain')
    require(action['reducedActionSha256'] == actual['reduction']['sha256'], 'native action/reduction source join')
    require(assembly_packet['provenance']['ACTIONS_SHA256'] == actual['actions']['sha256'],
            'saved assembly/action source join')
    require(grades['inputPackets'][actual['actions']['path']] == actual['actions']['sha256'],
            'saved grades/action source join')
    require(grades['inputPackets'][actual['assembly']['path']] == actual['assembly']['sha256'],
            'saved grades/assembly source join')
    z, zp, xi, alpha = (state[n] for n in ('z', 'zp', 'xi', 'regulator'))
    eps, eta, sigma = (state['symbols'][n] for n in ('epsilon_shape', 'eta_bg', 'sigma_W'))
    kn = context['normalMomentum']; generators = (eps, eta, sigma)
    require(tuple(grades['generators']) == generators, 'native grade ordering')
    require(len(action['fields']) == len(action['probes']) == len(action['columns']) == 5,
            'complete native five-field action')
    parameters = {name: sp.Rational(value) for name, value in physical['parameters'].items()}
    profile_formulas = {name: sp.sympify(value, locals={'xi': xi})
                        for name, value in physical['profiles'].items()}
    require(all(v.free_symbols <= {xi} for v in profile_formulas.values()), 'fixed source profile binding')
    endpoint_values = uniform['profileEndpoints']
    p0, inverse = candidate['fixedSymbol'], candidate['fixedInverse']
    require(p0 == saved('fixedSymbol') and inverse == saved('fixedInverse'), 'candidate/source operand join')
    known = dict(reduction['dimensionState']['known'])
    for packet in (action, assembly_packet, grades, uniform):
        known.update(packet['dimensionState']['known'])
    fields = list(map(tuple, units['fieldUnits'])); rows = list(map(tuple, units['rowUnits']))
    ku = tuple(units['spectralMeasureUnit']); length_unit = tuple(-v for v in ku)
    require(all(tuple(action['columnUnits'][(k, i)]) == rows[i]
                for i in range(5) for k in range(5)), 'native equation-unit join')
    require([tuple(known[f]) for f in action['fields']] == fields, 'native field-unit join')
    column_units = [tuple(a-b for a, b in zip(row, ku)) for row in rows]
    frame_units = [[tuple(a-b for a, b in zip(field, c)) for c in column_units] for field in fields]
    force_units = [[tuple(a-b for a, b in zip(row, c)) for c in column_units] for row in rows]
    require(frame_units == [[tuple(v) for v in row] for row in units['kernelEntryUnits']],
            'residue-frame units from inverse spectral measure')
    j.value('units-and-binding.pickle', {'fieldUnits': fields, 'rowUnits': rows,
        'residueColumnUnits': column_units, 'endFieldUnits': frame_units,
        'forcingUnits': force_units, 'physicalInput': physical,
        'profileFormulas': profile_formulas, 'profileEndpoints': endpoint_values,
        'profileRegulator': alpha, 'nativeLimitsUnchanged': True})

    coordinate_symbols = {z, zp, xi, alpha, kn, state['transfer'],
                          *state['normal_map'].values()}
    coordinate_symbols = {v for v in coordinate_symbols if isinstance(v, sp.Symbol)}
    @lru_cache(maxsize=None)
    def bind(expression):
        if isinstance(expression, sp.MatrixBase):
            return sp.ImmutableMatrix(expression).applyfunc(bind)
        if expression in endpoint_values:
            return endpoint_values[expression]
        if isinstance(expression, sp.Symbol):
            if expression in generators or expression in coordinate_symbols or isinstance(expression, sp.Dummy):
                return expression
            require(expression.name in parameters, ('unbound material input', expression))
            return parameters[expression.name]
        if isinstance(expression, AppliedUndef) and expression.func in state['profiles'].values():
            name = next(n for n, f in state['profiles'].items() if f == expression.func)
            return profile_formulas[name].subs(xi, bind(expression.args[0]))
        if not expression.args:
            return expression
        if isinstance(expression, sp.Integral):
            value = bind(expression.function)
            return sp.Integral(value, *(sp.Tuple(v, bind(a), bind(b)) for v,a,b in expression.limits)) if value != 0 else sp.S.Zero
        if isinstance(expression, sp.Derivative):
            value = bind(expression.expr)
            require(not value.has(sp.Integral, sp.Limit), 'derivative of nonlocal integral needs explicit disposition')
            return sp.diff(value, *expression.variable_count)
        if isinstance(expression, sp.Subs):
            return bind(expression.expr).subs(list(zip(expression.variables, map(bind, expression.point))), simultaneous=True)
        if isinstance(expression, sp.Limit):
            raise ValueError(('unbound profile limit', expression))
        return expression.func(*(bind(v) for v in expression.args))

    class UnitAudit:
        """Strict source-unit checker: no guessed units or new constraint solve."""
        def same(self, a, b):
            require(a is None or b is None or a == b, ('unit mismatch', a, b))

        @lru_cache(maxsize=None)
        def measure(self, v):
            if v == 0:
                return None
            if v in known:
                return tuple(known[v])
            if isinstance(v, sp.Symbol) and re.fullmatch(r'[wm]1_profile(?:_d[123](?:d[123])*)?', v.name):
                return (0,0,0)  # exact native DimensionAnalysis rule
            if isinstance(v, AppliedUndef):
                require(v.func in known, ('unknown source function units', v.func))
                return tuple(known[v.func])
            if v.is_number or v in (sp.true, sp.false):
                return (0, 0, 0)
            if isinstance(v, sp.Add):
                values = [self.measure(a) for a in v.args if self.measure(a) is not None]
                for other in values[1:]:
                    self.same(values[0], other)
                return values[0] if values else None
            if isinstance(v, sp.Mul):
                values = [self.measure(a) for a in v.args if self.measure(a) is not None]
                return tuple(sum(a[i] for a in values) for i in range(3))
            if isinstance(v, sp.Pow):
                self.same(self.measure(v.exp), (0,0,0)); base = self.measure(v.base)
                if v.exp.is_number:
                    return tuple(v.exp*a for a in base) if base else None
                self.same(base, (0,0,0)); return (0,0,0)
            if isinstance(v, sp.Derivative):
                unit = self.measure(v.expr)
                return tuple(unit[i]-sum(n*self.measure(x)[i] for x,n in v.variable_count) for i in range(3)) if unit else None
            if isinstance(v, sp.Subs):
                if isinstance(v.expr, sp.Derivative):
                    bound = dict(zip(v.variables, v.point))
                    measured = self.measure(v.expr.expr)
                    return tuple(measured[i]-sum(n*self.measure(bound.get(x,x))[i]
                        for x,n in v.expr.variable_count) for i in range(3)) if measured else None
                for a,b in zip(v.variables, v.point):
                    self.same(self.measure(a), self.measure(b))
                return self.measure(v.expr)
            if isinstance(v, sp.Integral):
                unit = self.measure(v.function)
                return tuple(unit[i]+sum(self.measure(l[0])[i] for l in v.limits) for i in range(3)) if unit else None
            if isinstance(v, sp.Limit):
                return self.measure(v.args[0])
            if v.func == sp.DiracDelta:
                n = v.args[1] if len(v.args)>1 else 0
                return tuple(-(n+1)*a for a in self.measure(v.args[0]))
            if v.func in (sp.exp, sp.sin, sp.cos, sp.tanh, sp.log, sp.Heaviside, sp.erf, sp.erfc):
                self.same(self.measure(v.args[0]), (0,0,0)); return (0,0,0)
            if v.func == sp.sign:
                return (0,0,0)
            if v.func in (sp.conjugate, sp.re, sp.im, sp.Abs):
                return self.measure(v.args[0])
            if isinstance(v, sp.Piecewise):
                values = [self.measure(a.expr) for a in v.args]
                for value in values[1:]: self.same(values[0], value)
                return values[0]
            raise ValueError(('unsupported unit node', v.func, v))
    unit = UnitAudit()
    for generator in generators:
        unit.same(unit.measure(generator), (0,0,0))
    # Numerical binding removes dimensional symbols. Keep an explicit typed
    # matrix-product audit, tied to the accepted row/field/inverse units.
    def add_units(*values):
        return tuple(sum(v[n] for v in values) for n in range(3))
    def neg_units(value):
        return tuple(-v for v in value)
    typed_products = []
    for i in range(5):
        for col in range(5):
            for a in range(5):
                for b in range(5):
                    punit = tuple(units['symbolEntryUnits'][a][b])
                    unit.same(punit, add_units(rows[a], neg_units(fields[b])))
                    au = frame_units[i][a]; bu = tuple(units['inverseEntryUnits'][i][a])
                    ar = frame_units[b][col]; br = tuple(units['inverseEntryUnits'][b][col])
                    units_c2 = add_units(au, punit, ar)
                    terms_c1 = (add_units(bu, punit, ar), add_units(au, punit, br),
                                add_units(au, punit, neg_units(ku), ar))
                    for measured in (*terms_c1, add_units(length_unit, units_c2)):
                        unit.same(measured, frame_units[i][col])
                    typed_products.append((i,col,a,b,units_c2,terms_c1))
    j.value('typed-principal-products.pickle', {'products':typed_products,
        'sourceUnits':units, 'numericBindingAfterUnitAudit':True,
        'scope':'Matrix contraction units; not a numerical channel normalization'})

    # Reuse the exact grade coefficients but retain native infinite integrals.
    def prepare_native(assembly, grades, action):
        records = {tuple(v['address']): v['record'] for v in grades['records'].values()}
        inventory = text_json('gradeInventory')
        selected = {tuple(v['address']):v for v in inventory if v['address'][0] in ('factor','source')}
        require(len(selected) == 115 and all(v['grades'] == [[0,0,0]] for v in selected.values()),
                'saved 80 factor/35 source records have only zero grade')
        require(set(selected) == {a for a in records if a[0] in ('factor','source')},
                'complete factor/source inventory-to-packet address join')
        for address, item in selected.items():
            record = records[address]
            require(set(record['COEFFICIENTS']) == {(0,0,0)}
                    and not record['ORIGINAL'].has(*generators),
                    ('actual saved factor/source record is not grade-free', address))
        j.value('factor-source-grade-inventory-join.pickle',
                {'inventory':selected,'records':{a:records[a] for a in selected},
                 'actualPacketGradesAllZero':True,'unboundedInterchangeClaimed':False})
        cells = {(v['ROW'],v['COLUMN']):v for v in assembly['ROWS']}
        require(set(cells) == {(i,k) for i in range(5) for k in range(5)}, 'all native cells')
        native = {}; census = []; unit_records = []
        for i in range(5):
            for k in range(5):
                cell = cells[i,k]; probe = action['probes'][k]
                require(cell['RECONSTRUCTION_RESIDUAL'] == 0 and cell['NONLINEAR_OR_AFFINE_REMAINDER'] == 0,
                        'accepted assembly residual status')
                local = {g:[] for g in GRADES}; nonlocal_terms = {g:[] for g in GRADES}
                for order, coefficient in cell['LOCAL'].items():
                    record = records['local',order,i,k]
                    require(record['ORIGINAL'] == coefficient, 'saved local coefficient identity')
                    expected = tuple(a-b+order*c for a,b,c in zip(rows[i],fields[k],length_unit))
                    unit.same(unit.measure(coefficient), expected)
                    require(tuple(record['UNIT']) == expected, 'saved local grade units')
                    for g in GRADES:
                        c = record['COEFFICIENTS'].get((0,*g), sp.S.Zero)
                        unit.same(unit.measure(c), expected)
                        local[g].append((order, c))
                    unit_records.append({'kind':'local','row':i,'column':k,'order':order,'unit':expected})
                for term, (integral, coefficient) in enumerate(cell['NONLOCAL']):
                    record = records['cell',i,k,term]
                    require(record['ORIGINAL'] == coefficient and not integral.has(*generators),
                            'saved grade-free integral/coefficient identity')
                    require(all(len(limit) == 3 for node in integral.atoms(sp.Integral)
                                for limit in node.limits), 'original definite ordered integral limits')
                    expected = tuple(a-b for a,b in zip(rows[i],unit.measure(integral)))
                    unit.same(unit.measure(coefficient), expected)
                    require(tuple(record['UNIT']) == expected, 'saved nonlocal grade units')
                    for g in GRADES:
                        c = record['COEFFICIENTS'].get((0,*g),sp.S.Zero)
                        unit.same(unit.measure(c), expected)
                        nonlocal_terms[g].append((term, c, integral))
                    unit_records.append({'kind':'nonlocal','row':i,'column':k,'term':term,'unit':expected})
                for g in GRADES:
                    native[g,i,k] = {'local':local[g], 'nonlocal':nonlocal_terms[g], 'probe':probe}
                census.append({'row':i,'column':k,'localTerms':len(cell['LOCAL']),
                               'nonlocalTerms':len(cell['NONLOCAL'])})
        require(sum(v['nonlocalTerms'] for v in census) == len(grades['termJoins']), 'all saved nonlocal terms retained')
        require(not any(g[0] for record in records.values() for g in record['COEFFICIENTS']),
                'linear field epsilon already removed in source actions')
        return native, census, unit_records
    native, census, unit_records = j.op('join-native-grades-and-units', prepare_native, assembly, grades, action)
    # This join transfers the already accepted source-to-symbol calculation.
    origin = {eta:sp.S.Zero,sigma:sp.S.Zero}
    def canonical_origin(node):
        # Structural grade selection only; no Fourier/operator construction.
        if isinstance(node, sp.MatrixBase):
            return sp.ImmutableMatrix(node).applyfunc(canonical_origin)
        if node in origin: return origin[node]
        if not node.args: return node
        args = tuple(canonical_origin(v) for v in node.args)
        if args == node.args: return node
        if isinstance(node, sp.Integral):
            return sp.Integral(*args) if args[0] != 0 else sp.S.Zero
        if isinstance(node, sp.Derivative):
            require(not args[0].has(sp.Integral, sp.Limit), 'origin derivative needs disposition')
            return sp.diff(args[0], *args[1:])
        if isinstance(node, sp.Subs):
            require(not args[0].has(sp.Integral,sp.Limit), 'origin Subs needs nonlocal disposition')
            return sp.Subs(*args).doit(deep=False)
        return node.func(*args)
    strong_zero = j.op('native-origin-selection', canonical_origin, action['strong'])
    reference_strong = j.op('reference-origin-normal-form', canonical_origin,
        uniform['records']['REFERENCE']['background']['strong'])
    j.value('native-reference-action-join.pickle', {
        'nativeStrongAtZero':strong_zero,'acceptedReferenceBackground':reference_strong,
        'acceptedReferenceSymbol':uniform['records']['REFERENCE']['strong'],
        'nativeGradeActions':native, 'sourceActionColumns':action['columns'],
        'sourceToSymbolEvidence':'accepted uniform-source strong_matrix route; no replay',
        'newPlaneJetRule':'differentiate the accepted harmonic-action identity in normal momentum',
        'integralEvaluationPerformed':False})
    require(strong_zero == reference_strong, 'native grade-zero/reference background structural join')
    plane_jet_premise = {
        'status':'UNVERIFIED_NONLOCAL_PLANE_JET_EXTENSION',
        'zeroGradeContainsAbelRegulator':bool(strong_zero.has(alpha)),
        'zeroGradeProfileOccurrences':tuple(sorted((v for v in strong_zero.atoms(AppliedUndef)
            if v.func in state['profiles'].values()),key=sp.default_sort_key)),
        'translationInvarianceEstablished':False,
        'momentumDerivativePassesOrderedIntegralsAndAbelLimitEstablished':False,
        'required':'Justify the native regulated harmonic identity on finite momentum derivatives near each inherited pole.',
        'scope':'The absence checks are necessary source diagnostics, not a translation/interchange theorem.'}
    j.value('plane-jet-premise.pickle',plane_jet_premise)
    save(out/'plane-jet-premise.json', {k:v for k,v in plane_jet_premise.items()
        if k != 'zeroGradeProfileOccurrences'} | {
        'zeroGradeProfileOccurrenceCount':len(plane_jet_premise['zeroGradeProfileOccurrences'])})
    reference_bound = j.op('bind-reference-symbol',bind,uniform['records']['REFERENCE']['strong'])
    zero_check('reference-symbol-join',reference_bound-p0)
    j.value('native-grades-and-units.pickle', {'actions':native,'census':census,'units':unit_records})
    save(out/'native-term-census.json',census)

    # Distinguish local meromorphic products from distribution products.
    branch = saved('branchEvidence'); radical_square = branch['minusSquare']
    blocks = {b['index']:b for b in candidate['blocks']}
    require(set(blocks) == {16,17}, 'inherited real-block set changed')
    end_tables = {}; end_symbols = {}
    for end in ('LEFT','RIGHT'):
        table = saved(end.lower()+'Pencil')
        # Names are source-defined coordinate identities, not new roots.
        def map_end(expression):
            mapping = {}
            for s in expression.free_symbols:
                if s.name == 's11cdCurrentRightMomentum': mapping[s] = kn
                elif s.name == 's11cdAcousticRightNormalMomentum': mapping[s] = sp.I*sp.sqrt(radical_square)
                elif s.name in ('s11cdPairingRightFrequency','s11cdPairingLeftFrequency'): mapping[s] = parameters['omega']
                elif s in (eta,sigma): mapping[s] = s
                elif s.name in parameters: mapping[s] = parameters[s.name]
                else: raise ValueError(('unmapped end-pencil coordinate',end,s))
            return expression.xreplace(mapping)
        symbolic = j.op('bind-'+end+'-symbol',bind,uniform['records'][end]['strong'])
        table_pencil = j.op('map-'+end+'-pencil',map_end,table['pencil'])
        zero_check(end+'-uniform-pencil-join',symbolic-table_pencil)
        zero_check(end+'-reference-join',table_pencil.xreplace(origin)-p0)
        end_symbols[end] = symbolic
        end_tables[end] = {}
        for g in GRADES:
            value = table['pencilCoefficients'][(*g,0)]
            mapped = j.op(end+'-coefficient-'+str(g[0])+str(g[1]),map_end,value)
            end_tables[end][g,0] = mapped
            # Select from the saved symbol for this new consumer join only.
            selected = j.op(end+'-symbol-grade-'+str(g[0])+str(g[1]),
                lambda a,g:a.diff(eta,g[0]).diff(sigma,g[1]).xreplace(origin),symbolic,g)
            zero_check(end+'-grade-join-'+str(g[0])+str(g[1]),mapped-selected)
            derivative = j.op(end+'-derivative-map-'+str(g[0])+str(g[1]),map_end,
                              table['pencilCoefficients'][(*g,1)])
            end_tables[end][g,1] = derivative
            zero_check(end+'-branch-derivative-join-'+str(g[0])+str(g[1]),
                       derivative-mapped.diff(kn))

    def convolution(a,b,degree):
        return [simplify(sum(a[t]*b[n-t] for t in range(n+1))) for n in range(degree+1)]
    def polynomial_series(expression,variable,radical,k0,rseries,degree):
        poly = sp.Poly(expression,variable,radical,domain='EX')
        result = [sp.S.Zero]*(degree+1)
        for (a,b),coefficient in poly.terms():
            require(not coefficient.has(variable,radical), 'unresolved polynomial coefficient')
            kseries = [sp.binomial(a,n)*k0**(a-n) if n<=a else sp.S.Zero for n in range(degree+1)]
            power = [sp.S.One]+[sp.S.Zero]*degree
            for _ in range(b): power = convolution(power,rseries,degree)
            term = convolution(kseries,power,degree)
            result = [simplify(x+coefficient*y) for x,y in zip(result,term)]
        return result
    def regular_part(chart,args,residue,k0,q0):
        value,variable,square,radical = args
        require(chart['backJoin'] == 0 and chart['rationalInCoordinates'], 'accepted inverse chart identity')
        require(variable == kn and square == radical_square, 'same local branch chart')
        r0 = simplify(q0/sp.I)
        require(r0.is_positive is True and simplify(r0*r0-square.subs(kn,k0)) == 0,
                'inherited nonzero positive radical at pole')
        degree = 3  # denominator order <=2; constant inverse term needs d+1.
        radpoly = sp.Poly(square,kn)
        rhs = [simplify(sp.diff(square,kn,n).subs(kn,k0)/sp.factorial(n)) for n in range(degree+1)]
        require(radpoly.degree() == 2, 'bounded quadratic branch')
        rseries = [r0]
        for n in range(1,degree+1):
            rseries.append(simplify((rhs[n]-sum(rseries[t]*rseries[n-t] for t in range(1,n)))/(2*r0)))
        numerator,denominator = sp.fraction(chart['rational'])
        nseries = polynomial_series(numerator,kn,radical,k0,rseries,degree)
        dseries = polynomial_series(denominator,kn,radical,k0,rseries,degree)
        order = next((i for i,v in enumerate(dseries) if v != 0),None)
        require(order is not None and order <= 2, 'unsupported local denominator order')
        require(dseries[order].is_zero is False, 'denominator leading coefficient nonzero undecided')
        negative = nseries[:max(0,order-1)]
        residue_join = simplify((nseries[order-1] if order else 0)-dseries[order]*residue)
        constant = simplify((nseries[order]-dseries[order+1]*residue)/dseries[order])
        return {'source':value,'chart':chart,'k0':k0,'q0':q0,'branchSeries':rseries,
                'numeratorSeries':nseries,'denominatorSeries':dseries,'denominatorOrder':order,
                'forbiddenHigherPoleCoefficients':negative,'residueJoin':residue_join,
                'savedResidue':residue,'regularCoefficient':constant}

    lifts = {}; exact_checks = []; local_data = {}
    for index in (17,16):
        block = blocks[index]; k0 = block['k']; residue = block['residue']
        require(residue == saved('residue'+str(index)), 'saved residue/frame join')
        q0 = saved('q'+str(index)); source_at_pole = saved('p'+str(index))
        regular = sp.zeros(5)
        for i in range(5):
            for k in range(5):
                chart = saved('chart%d%d'%(i,k)); args = saved('chartInput%d%d'%(i,k))
                require(args[0] == inverse[i,k], 'chart/actual inverse-entry join')
                coefficient = j.op('block%d-regular-%d-%d'%(index,i,k),regular_part,
                                   chart,args,residue[i,k],k0,q0)
                require(all(v == 0 for v in coefficient['forbiddenHigherPoleCoefficients'])
                        and coefficient['residueJoin'] == 0,'local quotient/residue source check')
                regular[i,k] = coefficient['regularCoefficient']
        regular = sp.ImmutableMatrix(regular)
        j.value('block%d-regular-matrix.pickle'%index,regular)
        pprime = j.op('block%d-source-derivative'%index,
                     lambda a,b:a.diff(kn).subs(kn,b),p0,k0)
        zero_check('block%d-residue-left'%index,source_at_pole*residue)
        zero_check('block%d-residue-right'%index,residue*source_at_pole)
        zero_check('block%d-inverse-constant-left'%index,source_at_pole*regular+pprime*residue-sp.eye(5))
        zero_check('block%d-inverse-constant-right'%index,regular*source_at_pole+residue*pprime-sp.eye(5))
        for end in ('LEFT','RIGHT'):
            for g in GRADES[1:]:
                tag='%s-block%d-grade%d%d'%(end,index,*g)
                ph = j.op(tag+'-end-coefficient',lambda a,b:a.subs(kn,b),end_tables[end][g,0],k0)
                dh = j.op(tag+'-end-derivative',lambda a,b:a.subs(kn,b),end_tables[end][g,1],k0)
                def principal(a,b,h,dh):
                    c2 = -a*h*a
                    c1 = -(b*h*a+a*h*b+a*dh*a)
                    return {'doublePole':matrix_reduce(c2),'simplePole':matrix_reduce(c1)}
                c = j.op(tag+'-principal-part',principal,residue,regular,ph,dh)
                polynomial = c['simplePole']+sp.I*z*c['doublePole']
                # Separate application of the actual source symbol to the new
                # degree-one field, not a copy of the construction formula.
                residual = j.op(tag+'-end-action',lambda p,dp,h,a,f:
                    sp.ImmutableMatrix(p*f-sp.I*dp*f.diff(z)+h*a),
                    source_at_pole,pprime,ph,residue,polynomial)
                for power in (0,1):
                    selected = residual.applyfunc(lambda v:sp.diff(v,z,power).subs(z,0))
                    zero_check(tag+'-end-action-power'+str(power),selected)
                field = sp.exp(sp.I*k0*z)*polynomial
                lifts[end,index,g] = field
                j.value(tag+'-field.pickle',{'polynomial':polynomial,'field':field,
                         'principalPart':c,'endCoefficient':ph,'branchDerivative':dh,
                         'fieldEntryUnits':frame_units,'columnUnits':column_units})
                exact_checks.append({'end':end,'block':index,'grade':list(g),'exactResidualZero':True})
        local_data[index] = {'k0':k0,'residue':residue,'sourceAtPole':source_at_pole,'sourceDerivative':pprime}

    # Rename each native integral once, before inserting any lifted field.
    # Bound identifiers never capture an observation coordinate from the lift.
    @lru_cache(maxsize=None)
    def alpha_safe(expression):
        def visit(node,env):
            if node in env: return env[node]
            if isinstance(node,sp.Integral):
                renames = {lim[0]:sp.Dummy('liftBound_'+lim[0].name,**lim[0].assumptions0) for lim in node.limits}
                for old,new in renames.items(): known[new] = unit.measure(old)
                limits = []
                for n,(v,a,b) in enumerate(node.limits):
                    outer = {lim[0]:renames[lim[0]] for lim in node.limits[n+1:]}
                    scope = {**env,**outer}
                    limits.append((renames[v],visit(a,scope),visit(b,scope)))
                return sp.Integral(visit(node.function,{**env,**renames}),*limits)
            if isinstance(node,sp.Subs):
                renames = {v:sp.Dummy('liftSubs_'+v.name,**v.assumptions0) for v in node.variables}
                for old,point in zip(node.variables,node.point):
                    known[renames[old]] = unit.measure(point)
                return sp.Subs(visit(node.expr,{**env,**renames}),tuple(renames.values()),
                               tuple(visit(v,env) for v in node.point))
            if not node.args: return node
            return node.func(*(visit(a,env) for a in node.args))
        return visit(expression,{})

    def substitute_field(expression,probe,field):
        @lru_cache(maxsize=None)
        def visit(node):
            if isinstance(node,AppliedUndef) and node.func == probe:
                require(len(node.args)==1, 'single normal-coordinate field slot')
                return field.subs(z,node.args[0],simultaneous=True)
            if isinstance(node,sp.Integral):
                value = visit(node.function)
                return sp.Integral(value,*node.limits) if value != 0 else sp.S.Zero
            if isinstance(node,sp.Derivative):
                value = visit(node.expr)
                require(not value.has(sp.Integral,sp.Limit), 'unsupported derivative of integral in native action')
                return sp.diff(value,*node.variable_count)
            if isinstance(node,sp.Subs):
                return visit(node.expr).subs(list(zip(node.variables,node.point)),simultaneous=True)
            if not node.args: return node
            return node.func(*(visit(a) for a in node.args))
        return visit(expression)

    # Formal linearity at the retained regulator only. Extract coefficients
    # independent of ALL bound variables, retaining the original limit order.
    # No integral/limit is evaluated or moved through another limiting action.
    def linear_form(expression):
        @lru_cache(maxsize=None)
        def visit(node):
            if isinstance(node,sp.Integral):
                inner = sp.expand(visit(node.function),mul=True,power_base=False,power_exp=False,log=False)
                terms=[]
                bound=tuple(limit[0] for limit in node.limits)
                for term in sp.Add.make_args(inner):
                    if term == 0: continue
                    coefficient, dependent = term.as_independent(*bound, as_Add=False)
                    terms.append(coefficient*sp.Integral(dependent,*node.limits))
                return sp.Add(*terms)
            if not node.args: return node
            return node.func(*(visit(v) for v in node.args))
        return sp.expand(visit(expression),mul=True,power_base=False,power_exp=False,log=False)

    def formal_coefficients(matrix):
        normalized=sp.ImmutableMatrix(matrix).applyfunc(linear_form)
        carriers=set()
        def collect(node):
            if isinstance(node,sp.Integral):
                carriers.add(node);return  # nested carriers remain inside this operand
            for arg in node.args: collect(arg)
        for entry in normalized: collect(entry)
        ordered=tuple(sorted(carriers,key=sp.default_sort_key))
        symbols=tuple(sp.Dummy('nativeIntegralCarrier'+str(i)) for i in range(len(ordered)))
        binding=dict(zip(ordered,symbols));positions=[];coefficients=[]
        for i in range(normalized.rows):
            for col in range(normalized.cols):
                encoded=normalized[i,col].xreplace(binding)
                require(not encoded.has(sp.Integral,sp.Limit), 'unencoded integral or limit in formal carrier algebra')
                terms=sp.Poly(encoded,*symbols,domain='EX').terms() if symbols else [((),encoded)]
                for powers,coefficient in terms:
                    require(not coefficient.has(*symbols,sp.Integral,sp.Limit), 'carrier-dependent scalar coefficient')
                    positions.append({'row':i,'column':col,'carrierPowers':powers})
                    coefficients.append(coefficient)
        return {'normalized':normalized,'carriers':ordered,'symbols':symbols,
                'positions':positions,'coefficients':sp.ImmutableMatrix(len(coefficients),1,coefficients),
                'scope':'Free formal algebra on unchanged ordered integral carriers; no integrated-value claim.'}

    ell = parameters['L_W']; chi = {'RIGHT':(1+sp.tanh(z/ell))/2,'LEFT':(1-sp.tanh(z/ell))/2}
    require(simplify(chi['RIGHT']+chi['LEFT']-1)==0,'smooth partition')
    j.value('auxiliary-partition.pickle', {'chi':chi,'physicalInputChanged':False,
        'heldFixedUnderBackgroundGrades':True,'profileAbelRegulator':alpha,
        'sourceCanonicalHeavisideSubtractionUnchanged':True})
    def prepare_bound(native):
        prepared = {}
        for (g,i,k),cell in native.items():
            prepared[g,i,k] = {'probe':cell['probe'],
                'local':[(n,bind(c)) for n,c in cell['local']],
                'nonlocal':[(t,bind(c),bind(alpha_safe(v))) for t,c,v in cell['nonlocal']]}
        return prepared
    prepared = j.op('bind-native-action', prepare_bound, native)
    def apply_action(g,field):
        local = sp.zeros(5); nonlocal_matrix = sp.zeros(5); terms=[]
        for i in range(5):
            for k in range(5):
                cell = prepared[g,i,k]
                for col in range(5):
                    f = field[k,col]
                    local[i,col] += sum(c*sp.diff(f,z,n) for n,c in cell['local'])
                    for term,c,integral in cell['nonlocal']:
                        value = c*substitute_field(integral,cell['probe'],f) if c != 0 and f != 0 else sp.S.Zero
                        nonlocal_matrix[i,col] += value
                        terms.append({'row':i,'sourceColumn':k,'frameColumn':col,'term':term,
                                      'grade':g,'coefficient':c,'nativeIntegral':integral,'value':value})
        return {'local':sp.ImmutableMatrix(local),'nonlocal':sp.ImmutableMatrix(nonlocal_matrix),
                'total':sp.ImmutableMatrix(local+nonlocal_matrix),'terms':terms}

    def decompose(forcing,end_plane,reference_plane,side_actions,omission=None):
        commutators={};omitted=[]
        for end, action_record in side_actions.items():
            total=sp.MutableDenseMatrix(action_record['local'])
            for term in action_record['terms']:
                address=(end,term['row'],term['sourceColumn'],term['frameColumn'],term['term'])
                if address == omission:
                    omitted.append(term)
                else:
                    total[term['row'],term['frameColumn']] += term['value']
            commutators[end]=sp.ImmutableMatrix(total)-chi[end]*reference_plane[end]
        require(len(omitted) == (1 if omission is not None else 0), 'exactly the specified native addend omitted')
        value=-(forcing['total']-sum((chi[e]*end_plane[e] for e in ('LEFT','RIGHT')),sp.zeros(5))) \
              -sum(commutators.values(),sp.zeros(5))
        return {'decomposed':value,'conditionalCommutators':commutators,'omission':omission,'omittedTerms':omitted,
                'planeJetExtensionStatus':'UNVERIFIED_NONLOCAL_PLANE_JET_EXTENSION'}

    force_summaries=[]; domain_records=[]; responsive=False; formal_nonlocal_responsive=False
    for index in (17,16):
        data=local_data[index];k0=data['k0'];residue=data['residue']
        seed=sp.exp(sp.I*k0*z)*residue
        for g in GRADES[1:]:
            tag='block%d-grade%d%d'%(index,*g)
            lift=sum((chi[e]*lifts[e,index,g] for e in ('LEFT','RIGHT')),sp.zeros(5))
            forcing=j.op(tag+'-native-forcing',apply_action,g,seed)
            lifted=j.op(tag+'-native-lift-action',apply_action,(0,0),lift)
            side_actions={e:j.op(tag+'-'+e+'-cutoff-action',apply_action,(0,0),chi[e]*lifts[e,index,g])
                          for e in ('LEFT','RIGHT')}
            end_plane={e:sp.exp(sp.I*k0*z)*end_tables[e][g,0].subs(kn,k0)*residue for e in ('LEFT','RIGHT')}
            reference_plane={e:sp.exp(sp.I*k0*z)*(data['sourceAtPole']*(lifts[e,index,g]/sp.exp(sp.I*k0*z))
                              -sp.I*data['sourceDerivative']*(lifts[e,index,g]/sp.exp(sp.I*k0*z)).diff(z))
                             for e in ('LEFT','RIGHT')}
            assembly_result=j.op(tag+'-assemble-decomposition',decompose,forcing,end_plane,reference_plane,side_actions)
            decomposed=assembly_result['decomposed']
            direct=-forcing['total']-lifted['total']
            j.value(tag+'-forcing-operands.pickle',{'seed':seed,'lift':lift,'forcing':forcing,
                'liftedAction':lifted,'sideActions':side_actions,'endPlaneActions':end_plane,
                'referencePlaneJetActions':reference_plane,'conditionalCommutators':assembly_result['conditionalCommutators'],
                'direct':direct,'decomposed':decomposed,'forcingUnits':force_units,
                'planeJetReductionPremise':plane_jet_premise,'wholeLinePairingEstablished':False})
            # Both sides of the unverified plane-jet identification must survive
            # any later reconstruction failure. Do not move this below a guard.
            raw_plane={e:j.op(tag+'-'+e+'-native-plane-action',apply_action,(0,0),lifts[e,index,g])
                       for e in ('LEFT','RIGHT')}
            j.value(tag+'-plane-jet-reduction-operands.pickle',{
                'nativePlaneActions':raw_plane,'sourceSymbolPlaneActions':reference_plane,
                'referenceSourceJoin':'native-reference-action-join.pickle',
                'numericalOrIntegralEvaluation':False,'premise':plane_jet_premise})
            reconstruction=j.op(tag+'-reconstruction-carriers',formal_coefficients,direct-decomposed)
            zero_check(tag+'-reconstruction-coefficients-zero',reconstruction['coefficients'])
            if not formal_nonlocal_responsive:
                candidates=[]
                for end in ('LEFT','RIGHT'):
                    for term in side_actions[end]['terms']:
                        if term['value'] != 0:
                            candidates.append((end,term['row'],term['sourceColumn'],term['frameColumn'],term['term']))
                j.value(tag+'-nonlocal-mutation-candidates.pickle',candidates)
                for number,address in enumerate(candidates):
                    # Reassemble from the actual source-term list while omitting
                    # one negative addend of the correct decomposition. Hold the
                    # independently computed direct forcing fixed throughout.
                    changed=j.op(tag+'-nonlocal-omission-'+str(number),decompose,
                                 forcing,end_plane,reference_plane,side_actions,address)
                    raw_residual=direct-changed['decomposed']
                    formal=j.op(tag+'-mutated-carriers-'+str(number),formal_coefficients,raw_residual)
                    reduced=j.op(tag+'-mutated-coefficients-'+str(number),matrix_reduce,formal['coefficients'])
                    surviving=[i for i,(position,value) in enumerate(zip(formal['positions'],reduced))
                               if any(position['carrierPowers']) and value.is_zero is False]
                    j.value(tag+'-nonlocal-mutation-'+str(number)+'.pickle',{
                        'address':address,'removedActualTerms':changed['omittedTerms'],
                        'baselineDirect':direct,'baselineDecomposition':decomposed,
                        'mutatedDecomposition':changed['decomposed'],'actualResidual':raw_residual,
                        'carrierReduction':formal,'reducedCoefficients':reduced,
                        'exactNonzeroCarrierCoefficientPositions':surviving,
                        'formalControlResponsive':bool(surviving),
                        'integratedResponseNonzeroEstablished':False})
                    if surviving:
                        formal_nonlocal_responsive=True;break
            if not responsive:
                # Independent exact end-lift control, still not an assertion of
                # a nonzero evaluated nonlocal integral or outgoing response.
                for end in ('LEFT','RIGHT'):
                    action_change=j.op(tag+'-'+end+'-end-mutation',matrix_reduce,
                        end_tables[end][g,0].subs(kn,k0)*residue)
                    flags=[v.is_zero for v in action_change]
                    if any(flag is False for flag in flags):
                        j.value(tag+'-end-mutation.pickle',{'removedEndLift':lifts[end,index,g],
                            'heldFixedEndForcing':end_plane[end],'residualWithoutLift':action_change,
                            'responsiveExactCoefficient':True,'nonlocalControlSubstitute':False})
                        responsive=True;break
            local_part=-forcing['local']-lifted['local']
            # These source-specific local limits are diagnostic. Local and
            # nonlocal tails can cancel; no conclusion about the full defect
            # follows from checking the local summand in isolation.
            local_limits=[]
            j.value(tag+'-local-tail-input.pickle', {'localPart':local_part,'phase':sp.exp(sp.I*k0*z),
                'weight':1+(z/ell)**2,'coordinate':z,'ends':(-sp.oo,sp.oo)})
            for i in range(5):
                for col in range(5):
                    stripped=local_part[i,col]/sp.exp(sp.I*k0*z)
                    for end in (-sp.oo,sp.oo):
                        operand=sp.Limit((1+(z/ell)**2)*stripped,z,end)
                        local_limits.append({'row':i,'column':col,'end':str(end),
                            'operand':operand,'status':'UNEVALUATED_LOCAL_SUMMAND_NOT_FULL_DEFECT'})
            obligations=[]
            for end, action_record in [('PROFILE_FORCING', forcing), *side_actions.items()]:
                for term in action_record['terms']:
                    if term['value']==0:continue
                    obligations.append({'side':end,'address':{k:term[k] for k in ('row','sourceColumn','frameColumn','term')},
                        'operand':term['value'],'orderedIntegralLimits':[v.limits for v in sorted(term['value'].atoms(sp.Integral),key=sp.default_sort_key)],
                        'hasProfileAbelRegulator':term['value'].has(alpha),
                        'required':'combine with the other side and profile forcing before the source Abel weak limit; establish weighted tail/pole-pairing regularity',
                        'status':'UNRESOLVED_UNBOUNDED_NONLOCAL_PAIRING'})
            record={'block':index,'grade':g,'directForcing':direct,'profileForcingTerms':forcing['terms'],
                'localWeightedTailLimits':local_limits,'nonlocalObligations':obligations,
                'fullDefectEndLimitOperands':{str(v):direct.applyfunc(lambda a:sp.Limit(a,z,v)) for v in (-sp.oo,sp.oo)},
                'profileAbelLimitOrder':'after combination/convolution; unchanged source convention',
                'realPoleExclusionIndependentOfProfileAbel':True,
                'schwartzMembershipEstablished':False,'wholeLineOutgoingActionEstablished':False,
                'status':'FORCING_DOMAIN_UNRESOLVED'}
            j.value(tag+'-domain.pickle',record)
            summary={'block':index,'grade':list(g),'exactEndEquations':True,
                     'formalNativeActionReconstructionZero':True,'nonlocalObligations':len(obligations),
                     'localWeightedLimitsEvaluated':0,
                     'localWeightedLimits':len(local_limits),'forcingDomainEstablished':False}
            save(out/(tag+'-summary.json'),summary);force_summaries.append(summary)
            domain_records.append({'block':index,'grade':list(g),'artifact':tag+'-domain.pickle',
                                   'status':'FORCING_DOMAIN_UNRESOLVED'})
    require(formal_nonlocal_responsive, 'no exact nonzero formal nonlocal carrier mutation; preserve all control operands')
    require(responsive,'no exact responsive end-lift mutation found; preserve actual outputs')
    save(out/'domain-status.json',{'status':'FORCING_DOMAIN_UNRESOLVED','records':domain_records,
        'preciseDependency':'source-specific full unbounded nonlocal cancellation, weighted tails and prescription pairing; local limits alone do not settle it',
        'planeJetExtensionStatus':'UNVERIFIED_NONLOCAL_PLANE_JET_EXTENSION',
        'formalNonlocalMutationResponsive':formal_nonlocal_responsive,'nonlocalMutationIntegratedResponseEstablished':False,
        'globalProfileTheoremRequired':False,'newResponseOrGreenOperatorClaimed':False})
    return {'status':'END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED','case':'LAB_HELD__RHO4_CONSTANT',
            'exactEndChecks':exact_checks,'forcingSummaries':force_summaries,
            'nativeCells':len(census),'nativeNonlocalTerms':sum(v['nonlocalTerms'] for v in census),
            'exactEndMutationResponsive':responsive,'formalNonlocalMutationResponsive':formal_nonlocal_responsive,
            'nonlocalIntegratedMutationEstablished':False,'planeJetExtensionEstablished':False,
            'wholeLineForcingDomainEstablished':False,'responseConstructed':False,'formCompleted':False,
            'a11Cleared':False,'a12Cleared':False,'newRoots':0,'newModeSolves':0,'newLUSolves':0,
            'producerCalls':0,'integralsEvaluated':False,'freshIndependentClearClaimed':False,
            'completedOperations':len(j.records),'automaticRetry':False}


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--input-manifest',type=Path,required=True)
    parser.add_argument('--gate-receipt',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args()
    spec=json.loads(args.input_manifest.read_text());gate=json.loads(args.gate_receipt.read_text())
    require(gate['status']=='READY_FOR_ONE_GUARDED_TWO_ASYMPTOTE_LIFT_STAGE','method gate incomplete')
    require(gate['workerSha256']==digest(Path(__file__)) and gate['inputManifestSha256']==digest(args.input_manifest),
            'method gate worker/input identity')
    require(gate['scopeExplicitlyApproved'] is True and gate['substantiveReviewFindingsClosed'] is True,
            'scope or independent method-review disposition missing')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False,
            'bounded stage duration/retry contract')
    require(gate['sharedGuardSha256']==digest(ROOT/'scripts/s11c_guarded_run.py'), 'shared guard identity')
    observed=containment()
    out=args.run_directory.resolve();out.relative_to(STORE);out.mkdir(parents=True,exist_ok=False)
    save(out/'native-containment.json',observed);save(out/'input-manifest.json',spec);save(out/'gate-receipt.json',gate)
    journal=Journal(out);started=time.monotonic()
    try:
        result=construct(spec,out,journal)
        posthashes={name:route(v['path']) for name,v in spec['inputs'].items()}
        save(out/'posthashes.json',posthashes)
        require(posthashes==spec['inputs'],'input posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started,allSourceHashesUnchanged=True)
        save(out/'operation-index.json',journal.records)
        save(out/'checks.json',result)
        save(out/'artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        print(json.dumps(result,indent=2,allow_nan=False))
    except BaseException:
        failure_traceback = traceback.format_exc()
        signal.alarm(0)  # bookkeeping only; no further science after a failure
        post = {}
        for name, record in spec['inputs'].items():
            try: post[name] = route(record['path'])
            except OSError as error: post[name] = {'error':str(error),'path':record['path']}
        save(out/'failure-posthashes.json',post)
        save(out/'partial-operation-index.json',journal.records)
        save(out/'failure.json',{'traceback':failure_traceback,'completedOperations':len(journal.records),
                                'incompleteOperation':journal.active,
                                'wallSeconds':time.monotonic()-started,'automaticRetry':False})
        save(out/'failure-artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        raise


if __name__=='__main__':
    main()
