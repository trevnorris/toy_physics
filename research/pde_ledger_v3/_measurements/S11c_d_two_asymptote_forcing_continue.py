#!/usr/bin/env python3
"""Resume only unfinished two-asymptote end-lift/forcing work.

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
import shutil
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
    """Restore prior returns without calling their functions; journal new work."""
    def __init__(self, out):
        self.out, self.records, self.stack = out, [], []
        self.next_number = 0
        self.restored_parent = self.restored_diagnostic = 0
        self.restored_artifacts = []

    @property
    def active(self):
        return self.stack[-1] if self.stack else None

    def value(self, name, value):
        path = self.out/name; path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            pickle.dump(value, stream, protocol=4)
            stream.flush(); os.fsync(stream.fileno())
        return route(path)

    def copy(self, name, record):
        path = Path(record['path'])
        require(route(path) == record, ('saved bytes changed', path))
        target = self.out/name; target.parent.mkdir(parents=True, exist_ok=True)
        with path.open('rb') as src, target.open('xb') as dest:
            shutil.copyfileobj(src, dest); dest.flush(); os.fsync(dest.fileno())
        result = route(target)
        require(result['sha256'] == record['sha256'], 'copied return identity')
        return result

    def begin(self, name):
        folder = 'operations/%04d-%s' % (self.next_number, name)
        self.next_number += 1
        return folder, {'name': name, 'startedUtc': datetime.now(timezone.utc).isoformat()}

    def restore(self, previous, diagnostic=False):
        name = ('diagnostic-' if diagnostic else '')+previous['name']
        folder, record = self.begin(name)
        record.update(input=self.copy(folder+'/input.pickle', previous['input']),
                      execution=('RESTORED_DIAGNOSTIC_COMPLETE_RETURN' if diagnostic
                                 else 'RESTORED_PRIOR_COMPLETE_RETURN'), priorOperation=previous)
        save(self.out/folder/'started.json', record); self.stack.append(record)
        record['value'] = self.copy(folder+'/value.pickle', previous['value'])
        with Path(record['value']['path']).open('rb') as stream:
            result = SavedCodec(stream).load()
        record['finishedUtc'] = datetime.now(timezone.utc).isoformat()
        save(self.out/folder/'completed.json', record)
        self.records.append(record); self.stack.pop()
        if diagnostic: self.restored_diagnostic += 1
        else: self.restored_parent += 1
        return result

    def artifact(self, name, record, restore=False):
        result = self.copy(name, record)
        self.restored_artifacts.append({'source':record,'copy':result,
                                       'execution':'COPIED_PRIOR_COMPLETE_ARTIFACT'})
        if restore:
            with Path(result['path']).open('rb') as stream:
                return SavedCodec(stream).load()
        return result

    def op(self, name, function, *args):
        folder, record = self.begin(name)
        record.update(input=self.value(folder+'/input.pickle', args),
                      execution='NEW_UNFINISHED_OPERATION')
        save(self.out/folder/'started.json', record); self.stack.append(record)
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle', result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json', record)
        self.records.append(record); self.stack.pop()
        return result

    def reuse(self, name, previous, result, args):
        folder, record = self.begin(name)
        record.update(input=self.value(folder+'/input.pickle',args),
                      execution='REUSED_COMPLETE_RETURN_IDENTICAL_OPERANDS',priorOperation=previous)
        save(self.out/folder/'started.json',record);self.stack.append(record)
        record.update(value=self.copy(folder+'/value.pickle',previous['value']),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',record)
        self.records.append(record);self.stack.pop()
        return result

    def resume_input(self, name, function, previous_input):
        # The failed call never returned. Copy its exact saved input before
        # executing only that unfinished function under the existing guard.
        folder, record = self.begin(name)
        record.update(input=self.copy(folder+'/input.pickle',previous_input),
                      execution='NEW_UNFINISHED_OPERATION_FROM_SAVED_INPUT',
                      priorIncompleteInput=previous_input)
        save(self.out/folder/'started.json',record);self.stack.append(record)
        with Path(record['input']['path']).open('rb') as stream:
            args = SavedCodec(stream).load()
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle',result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',record)
        self.records.append(record);self.stack.pop()
        return result


def construct(spec, out, journal):
    import sympy as sp
    from sympy.core.function import AppliedUndef
    j = journal
    actual = {name:route(r['path']) for name,r in spec['inputs'].items()}
    save(out/'prehashes.json',actual)
    require(actual == spec['inputs'],'pinned source/input bytes changed')
    prior = {r['name']:j.restore(r) for r in spec['completedOperations']}
    require(j.restored_parent == 536 and len(prior) == 536,
            'all completed continuation returns must be restored')
    cache = {name[len('restore-'):]:value for name,value in prior.items()
             if name.startswith('restore-')}
    load_names = {'units-and-binding.pickle','plane-jet-premise.pickle',
                  'auxiliary-partition.pickle'} | set(spec['endFieldArtifacts'])
    bundles = {name:j.artifact(name,actual[key],restore=name in load_names)
               for name,key in spec['carryArtifacts'].items()}
    save(out/'restored-artifact-index.json',j.restored_artifacts)
    save(out/'restored-completed-stage.json',{
        'completedReturnsRestored':536,'priorFunctionsExecuted':0,
        'copiedSupportingArtifacts':len(bundles),
        'regularMatricesRestoredWithoutConstruction':2,
        'endFieldsRestoredWithoutConstruction':8,
        'firstUnfinishedOperation':'bind-native-action',
        'priorIncompleteInput':actual['incompleteNativeBindingInput'],
        'parentFailurePreserved':True})
    candidate, context, units = (cache[n] for n in ('candidate','context','units'))
    uniform, reduction, action, assembly_packet, grades = (
        cache[n] for n in ('uniform','reduction','actions','assembly','grades'))
    state = reduction['reductionState']; binding = bundles['units-and-binding.pickle']
    physical = binding['physicalInput']; parameters = {
        name:sp.Rational(value) for name,value in physical['parameters'].items()}
    profile_formulas = binding['profileFormulas']; endpoint_values = binding['profileEndpoints']
    frame_units, column_units, force_units = (binding[n] for n in
        ('endFieldUnits','residueColumnUnits','forcingUnits'))
    z, zp, xi, alpha = (state[n] for n in ('z','zp','xi','regulator'))
    eps, eta, sigma = (state['symbols'][n] for n in ('epsilon_shape','eta_bg','sigma_W'))
    generators = (eps,eta,sigma); kn = context['normalMomentum']
    p0, inverse = cache['fixedSymbol'], cache['fixedInverse']
    known = dict(reduction['dimensionState']['known'])
    for packet in (action,assembly_packet,grades,uniform):
        known.update(packet['dimensionState']['known'])
    coordinate_symbols = {z,zp,xi,alpha,kn,state['transfer'],*state['normal_map'].values()}
    coordinate_symbols = {v for v in coordinate_symbols if isinstance(v,sp.Symbol)}
    native,census,unit_records = prior['join-native-grades-and-units']
    plane_jet_premise = bundles['plane-jet-premise.pickle']
    radical_square = cache['branchEvidence']['minusSquare']
    blocks = {b['index']:b for b in candidate['blocks']}
    end_tables = {end:{(g,order):prior[end+('-coefficient-' if order==0 else '-derivative-map-')
                    +str(g[0])+str(g[1])] for g in GRADES for order in (0,1)}
                  for end in ('LEFT','RIGHT')}


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
    lifts = {}; exact_checks = []; local_data = {}
    for index in (17,16):
        block = blocks[index]
        local_data[index] = {'k0':block['k'],'residue':block['residue'],
            'sourceAtPole':cache['p'+str(index)],
            'sourceDerivative':prior['block%d-source-derivative'%index]}
        for end in ('LEFT','RIGHT'):
            for g in GRADES[1:]:
                tag = '%s-block%d-grade%d%d'%(end,index,*g)
                field = bundles[tag+'-field.pickle']
                lifts[end,index,g] = field['field']
                for power in (0,1):
                    check = json.loads((out/(tag+'-end-action-power'+str(power)+'.json')).read_text())
                    require(check == {'shape':[5,5],'zero':True,'unresolvedEntries':[]},
                            ('saved end-equation check missing',tag,power))
                exact_checks.append({'end':end,'block':index,'grade':list(g),
                                     'exactResidualZero':True})

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
        return visit(expression.xreplace(endpoint_values),{})

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

    ell = parameters['L_W']; partition = bundles['auxiliary-partition.pickle']
    chi = partition['chi']
    require(partition['physicalInputChanged'] is False
            and partition['heldFixedUnderBackgroundGrades'] is True
            and partition['sourceCanonicalHeavisideSubtractionUnchanged'] is True,
            'saved partition scope changed')

    def endpoint_evidence(native_action, endpoints):
        limits = set()
        for cell in native_action.values():
            for _,_,integral in cell['nonlocal']:
                limits.update(integral.atoms(sp.Limit))
        ordered = tuple(sorted(limits,key=sp.default_sort_key))
        return {'sourceEndpointMap':endpoints,'actualNativeLimits':ordered,
                'matched':[(v,endpoints[v]) for v in ordered if v in endpoints],
                'unknown':[v for v in ordered if v not in endpoints],
                'integralOrLimitEvaluation':False}
    evidence = j.op('endpoint-binding-source-evidence',endpoint_evidence,native,endpoint_values)
    save(out/'endpoint-binding-source-evidence.json',{
        'actualNativeLimitCount':len(evidence['actualNativeLimits']),
        'matched':[{'limit':sp.sstr(k),'savedValue':sp.sstr(v)} for k,v in evidence['matched']],
        'unknownLimits':[sp.sstr(v) for v in evidence['unknown']],
        'substitutionBeforeAlphaRenaming':True,'unknownLimitGuardRetained':True,
        'integralsOrLimitsEvaluated':False})

    def prepare_bound(native):
        prepared = {}
        for (g,i,k),cell in native.items():
            prepared[g,i,k] = {'probe':cell['probe'],
                'local':[(n,bind(c)) for n,c in cell['local']],
                'nonlocal':[(t,bind(c),bind(alpha_safe(v))) for t,c,v in cell['nonlocal']]}
        return prepared
    def resume_bound(saved_native):
        require(saved_native == native,'saved unfinished native input differs from restored source')
        return prepare_bound(saved_native)
    prepared = j.resume_input('bind-native-action',resume_bound,actual['incompleteNativeBindingInput'])
    save(out/'endpoint-binding-complete.json',{
        'status':'KNOWN_ENDPOINT_SUBSTITUTION_AND_NATIVE_BINDING_COMPLETED',
        'priorIncompleteInput':actual['incompleteNativeBindingInput'],
        'knownEndpointValuesChanged':False,'unknownLimitGuardRetained':True,
        'integralsOrLimitsEvaluated':False})
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
    require(gate['status']=='READY_FOR_ONE_GUARDED_TWO_ASYMPTOTE_FORCING_CONTINUATION','method gate incomplete')
    require(gate['workerSha256']==digest(Path(__file__)) and gate['inputManifestSha256']==digest(args.input_manifest),
            'method gate worker/input identity')
    require(gate['scopeExplicitlyApproved'] is True and gate['substantiveReviewFindingsClosed'] is True,
            'scope or independent method-review disposition missing')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False,
            'bounded stage duration/retry contract')
    require(gate['sharedGuardSha256']==digest(ROOT/'scripts/s11c_guarded_run.py'), 'shared guard identity')
    require(gate['supervisorSha256']==digest(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'), 'supervisor identity')
    approval_path=Path(gate['scopeApprovalPath'])
    require(gate['scopeApprovalSha256']==digest(approval_path),'continuation approval identity')
    approval=json.loads(approval_path.read_text())
    require(approval['status']=='AUTHORIZED_ONE_GUARDED_TWO_ASYMPTOTE_FORCING_CONTINUATION'
            and approval['workerSha256']==digest(Path(__file__))
            and approval['inputManifestSha256']==digest(args.input_manifest), 'continuation approval scope')
    require(gate['correctionDispositionSha256']==digest(Path(gate['correctionDispositionPath'])),
            'bounded correction disposition identity')
    require(str(args.run_directory.resolve())==gate['resultDirectory'], 'continuation result directory')
    observed=containment()
    out=args.run_directory.resolve();out.relative_to(STORE);out.mkdir(parents=True,exist_ok=False)
    save(out/'native-containment.json',observed);save(out/'input-manifest.json',spec);save(out/'gate-receipt.json',gate)
    journal=Journal(out);started=time.monotonic()
    try:
        result=construct(spec,out,journal)
        posthashes={name:route(v['path']) for name,v in spec['inputs'].items()}
        save(out/'posthashes.json',posthashes)
        require(posthashes==spec['inputs'],'input posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started,allSourceHashesUnchanged=True,
                      restoredParentOperations=journal.restored_parent,
                      restoredDiagnosticOperations=journal.restored_diagnostic,
                      priorFunctionsExecuted=0)
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
                                'incompleteOperation':journal.active,'incompleteOperationStack':journal.stack,
                                'restoredParentOperations':journal.restored_parent,
                                'restoredDiagnosticOperations':journal.restored_diagnostic,
                                'restoredArtifacts':journal.restored_artifacts,
                                'wallSeconds':time.monotonic()-started,'automaticRetry':False})
        save(out/'failure-artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        raise


if __name__=='__main__':
    main()
