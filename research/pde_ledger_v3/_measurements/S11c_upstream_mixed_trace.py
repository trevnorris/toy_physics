#!/usr/bin/env python3
"""Candidate saved-term closed-response trace; no producer or production edits.

Scientific imports occur only after pooled containment and immutable input pins.
This prints operands, identities and controls. Independent review and a separate
record must decide their meaning; successful execution is not physics acceptance.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import resource
import re
import sys
import textwrap
import time
import traceback
from types import SimpleNamespace

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(value, message):
    if value is not True:
        raise ValueError(message)


def save(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def extract_function(source, name):
    tree = ast.parse(source)
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name]
    require(len(nodes) == 1, 'unique native function ' + name)
    node = nodes[0]
    return ast.get_source_segment(source, node)


def literal_record(source, key):
    found = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Dict):
            continue
        for k, v in zip(node.keys, node.values):
            if isinstance(k, ast.Constant) and k.value == key:
                fields = {ast.literal_eval(a): b for a, b in zip(v.keys, v.values)}
                call = fields['value']
                require(isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                        and call.func.id == '_restore' and len(call.args) == 1,
                        'literal restore source ' + key)
                found.append(ast.literal_eval(call.args[0]))
    require(len(found) == 1, 'unique literal record ' + key)
    return found[0]


def containment():
    memory = 4 * 1024**3
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines()
                 if x.startswith('0::'))
    base = Path('/sys/fs/cgroup') / group.lstrip('/')
    actual = {k: (base/k).read_text().strip() for k in ('memory.max', 'memory.swap.max', 'pids.max')}
    actual.update(cgroup=str(base), affinity=sorted(os.sched_getaffinity(0)),
                  threads={k: os.environ.get(k) for k in THREADS})
    require(actual['memory.max'] == str(memory) and actual['memory.swap.max'] == '0'
            and actual['pids.max'] == '32' and len(actual['affinity']) == 1
            and all(v == '1' for v in actual['threads'].values()), 'pooled resource limits')
    require('S11C_POOLED_GUARD_MANIFEST' in os.environ, 'pooled guard required')
    require(resource.getrlimit(resource.RLIMIT_CPU) == (resource.RLIM_INFINITY,) * 2,
            'no inherited CPU deadline')
    resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    return actual | dict(nativeAddressSpace=memory, durationLimits=None)


class Journal:
    def __init__(self, out):
        self.out = out
        self.completed = []
        self.active = None

    def encode(self, x):
        if isinstance(x, (sp.Basic, sp.MatrixBase)):
            return dict(text=str(x), srepr=sp.srepr(x))
        if isinstance(x, dict):
            require(all(isinstance(k, str) for k in x), 'literal string journal keys')
            return {k: self.encode(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [self.encode(v) for v in x]
        return x

    def emit(self, name, value):
        path = self.out/(name + '.json')
        save(path, self.encode(value))
        return dict(path=path.name, sha256=sha(path), bytes=path.stat().st_size)

    def stage(self, name, inputs, function):
        self.active = name
        operand = self.emit(name + '-input', inputs)
        result = function()
        receipt = self.emit(name + '-return', result)
        self.completed.append(dict(name=name, input=operand, result=receipt))
        self.active = None
        return result

    def zero(self, name, left, right):
        residual = sp.factor(sp.cancel(sp.simplify(left-right)))
        self.emit(name, dict(left=left, right=right, residual=residual))
        require(residual == 0, name)



def assignment(source, function, target, occurrence=0):
    """Verbatim native assignment; call sites declare any diagnostic replacement."""
    body = ast.parse(extract_function(source, function)).body[0].body
    nodes = [n for n in body if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == target for t in n.targets)]
    require(len(nodes) > occurrence, 'native assignment ' + target)
    return ast.get_source_segment(extract_function(source, function), nodes[occurrence])


def selected_constructor(source, key, case, outer=None):
    raw = literal_record(source, key)
    node = ast.parse(raw, mode='eval').body
    def label(x):
        require(isinstance(x, ast.Call) and len(x.args) == 1, 'literal case label')
        return ast.literal_eval(x.args[0])
    if outer:
        node = next(x.args[1] for x in node.args if label(x.args[0]) == outer)
    matches = [x.args[1] for x in node.args if [label(k) for k in x.args[0].args] == list(case)]
    require(len(matches) == 1, 'unique selected native case')
    value = next(x.args[1] for x in matches[0].args if label(x.args[0]) == 'VALUE')
    return ast.get_source_segment(raw, value)


def science(manifest, J):
    scripts = ROOT/'research/pde_ledger_v3/scripts'
    c2 = (scripts/'S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
    native_json = json.loads(Path(manifest['nativeRecords']).read_text())
    native = {}
    for name, record in native_json.items():
        original = selected_constructor(Path(record['source']).read_text(), record['recordKey'],
                                        record['case'], record['outer'])
        J.emit('native-' + name + '-source-join', dict(original=original,
               selected=record['constructor'], byteIdentical=original == record['constructor']))
        require(original == record['constructor'], 'native constructor identity ' + name)
        native[name] = sp.sympify(original, locals={'Str': Str})
    named = lambda v, name: next(x for k, x in v if str(k) == name)
    response, kernel = native['response'], native['dtn']
    all_symbols = set().union(*(v.atoms(sp.Symbol) for v in native.values()))
    def atom(name):
        matches = [x for x in all_symbols if x.name == name]
        require(len(matches) == 1, 'unique native symbol ' + name)
        return matches[0]
    physical = json.loads(Path(manifest['physicalInput']).read_text())
    values = {k: sp.sympify(v) for k, v in physical['parameters'].items()}
    values['omega'] = sp.Integer(3)
    # Contrast variables remain independent; never substitute their saved numbers.
    bindings = {a: values[a.name] for a in all_symbols
                if a.name in values and a.name not in ('eta_bg', 'sigma_W', 'epsilon_shape')}
    eta, sigma, epsilon = (atom(x) for x in ('eta_bg', 'sigma_W', 'epsilon_shape'))
    prior = Path(manifest['priorDirectory'])
    saved_point = json.loads((prior/'physical-point.json').read_text())
    require(saved_point['originalParameters'] == physical['parameters'], 'saved physical input identity')
    def restore_json(name):
        raw = json.loads((prior/name).read_text())
        J.emit('restore-' + name.removesuffix('.json'), dict(original=raw,
               path=str(prior/name), sha256=sha(prior/name), functionCalled=False))
        return raw
    saved = restore_json('profile-action-return.json')
    saved_dimensions = restore_json('dimensions.json')
    point = restore_json('physical-point.json')
    # Only named scalar returns are decoded. No prior producer/function is called.
    restore = lambda v: sp.sympify(v['srepr'], locals={'Str': Str})
    bare_integrand = restore(saved['integrand'])
    qi_value, qo_value = restore(point['inputDepth']), restore(point['outputDepth'])
    qin, qout = atom('s11cc1_q_out_input'), atom('s11cc1_q_out_output')
    kin = tuple(atom('s11cc1_k_input_' + str(i)) for i in (1,2,3))
    kout = tuple(atom('s11cc1_k_output_' + str(i)) for i in (1,2,3))
    middle = tuple(sp.Symbol('trace_middle_k' + str(i), real=True) for i in (1,2,3))
    qm = sp.Symbol('trace_middle_q', nonzero=True)
    Q = restore(point['outputProfileMomentum'])
    J.zero('saved-frequency-binding',restore(point['omega']),values['omega'])
    J.zero('saved-input-momentum',restore(point['inputProfileMomentum']),0)
    J.zero('saved-output-momentum',Q,sp.Rational(1,10))
    edge = tuple(restore(x) for x in point['edge'])
    endpoint_binding = dict(zip(kin, (0,*edge))) | dict(zip(kout, (Q,*edge)))
    endpoint_binding |= {qin: qi_value, qout: qo_value, middle[1]:edge[0], middle[2]:edge[1]}
    z = atom('s11cc1_dtn_operator_lab_held_plus')
    identity = atom('s11cc1_identity_operator')
    resolvent = named(response, 'RESOLVENT')
    definition = named(response, 'RESOLVENT_DEFINITION')[1]
    coefficient = sp.expand(definition).coeff(z).subs(bindings)
    J.zero('resolvent-definition', definition.subs(bindings), identity + coefficient*z)
    ns = dict(sp=sp, NEW_DIMENSIONS={})
    exec(compile(extract_function(c2, 'fourier_profiles'), '<native-fourier_profiles>', 'exec'), ns)
    inputs = SimpleNamespace(a=atom)
    first = named(kernel, 'FIRST_SHAPE').subs(bindings)
    transfer = ns['fourier_profiles'](inputs, first, kout, kin)
    flat = named(kernel, 'FLAT_DIAGONAL')
    flat = flat.xreplace({d:sp.S.One for d in flat.atoms(sp.DiracDelta)}).subs(bindings)
    left = dict(zip(kin,middle)) | {qin:qm}
    right = dict(zip(kout,middle)) | {qout:qm}
    # Execute the native three-leg constructor verbatim: the original direct slot is zero.
    native_matrix = assignment(c2, 'kernel_bridge', 'z_three')
    tree = ast.parse(native_matrix).body[0]
    slot = tree.value.args[0].elts[0].elts[2]
    require(isinstance(slot,ast.Constant) and slot.value == 0, 'unchanged native bare direct slot')
    z_matrix = sp.Matrix([[flat, first],[0,flat.xreplace(dict(zip(kout,kin))|{qout:qin})]])
    ns.update(z_matrix=z_matrix, transfer=transfer, leftmap=left, rightmap=right,
              z_middle=flat.xreplace(dict(zip(kout,middle))|{qout:qm}))
    exec(native_matrix, ns)
    Z = ns['z_three'].subs(endpoint_binding)
    D = sp.Symbol('saved_bare_mixed_operator_slot')
    # Full operator slot includes the two conserved-edge deltas. They are
    # factored only in the final reduced integrand, never evaluated at zero.
    reduced_dimension = restore(saved_dimensions['reducedOneDimensionalKernel'])
    full_dimensions = restore(saved_dimensions['nativeKernel'])
    J.emit('kernel-measure',dict(reducedDimension=reduced_dimension,
        fullDimensions=full_dimensions,factoredEdgeDeltas=saved_dimensions['factoredEdgeDeltas'],
        restoredFullDimension=reduced_dimension+sp.Matrix([2,0,0]),
        note='Native first-shape products are formal composition integrands, not evaluated convolutions.'))
    require(saved_dimensions['factoredEdgeDeltas'] == 2,'inherited edge-delta count')
    require(all(x == reduced_dimension+sp.Matrix([2,0,0]) for x in full_dimensions),
            'full/reduced kernel dimensional join')
    tag = sp.Symbol('direct_slot_tag', real=True)
    added = sp.zeros(3); added[0,2] = tag*eta*sigma*D
    J.emit('closure-input',dict(nativeConstructor=native_matrix, bareBaseline=Z, insertion=added,
        coefficient=coefficient, definition=definition, bindings={str(k):v for k,v in bindings.items()},
        savedIntegrand=bare_integrand, selectedMomentum={str(k):v for k,v in endpoint_binding.items()},
        note='Only the missing direct slot changes. Existing first-shape iteration is present once.'))
    def close(matrix):
        context = dict(sp=sp, coefficient=coefficient, z_three=matrix)
        statement = assignment(c2,'kernel_bridge','three_inverse')
        exec(statement, context)
        return context['three_inverse']*matrix
    def grade(value):
        return sp.cancel(sp.diff(value,eta,sigma).subs({eta:0,sigma:0},simultaneous=True))
    def matzero(name, left_operand, right_operand):
        residual = (left_operand-right_operand).applyfunc(sp.cancel)
        J.emit(name,dict(left=left_operand,right=right_operand,residual=residual))
        require(residual == sp.zeros(*residual.shape),name)
    def closure():
        old, new = close(Z), close(Z+added)
        J.emit('closure-raw',dict(old=old,new=new,oldMixed=grade(old[0,2]),
               newMixed=grade(new[0,2]), difference=new-old))
        matzero('old-boundary-closure',(sp.eye(3)+coefficient*Z)*old,Z)
        matzero('new-boundary-closure',(sp.eye(3)+coefficient*(Z+added))*new,Z+added)
        delta = (new-old).applyfunc(sp.cancel)
        mixed = grade(delta[0,2])
        selected = sp.cancel(mixed.subs(tag,1))
        factor = sp.cancel(sp.diff(selected,D))
        J.zero('direct-linearity',selected,factor*D)
        # Independent resolvent-difference identity; both external response factors are retained.
        check = (sp.eye(3)+coefficient*Z).upper_triangular_solve(added)
        right_inverse = (sp.eye(3)+coefficient*(Z+added)).upper_triangular_solve(sp.eye(3))
        matzero('resolvent-difference-identity',delta,check*right_inverse)
        return dict(old=old,new=new,delta=delta,mixed=mixed,factor=factor)
    C = J.stage('closed-response',dict(baseline=Z,insertion=added,coefficient=coefficient),closure)

    def reference():
        source = native['trace'][0]
        p, jet = atom('delta_p_plus'), atom('d_w_delta_p_plus')
        trace = sp.cancel(source/epsilon).subs(bindings)
        value_coefficient = sp.diff(trace,p)
        height = sp.diff(trace,jet)
        constant = trace.subs({p:0,jet:0}, simultaneous=True)
        profile = atom('w1_profile')
        height_constant = height.subs(profile,0)
        height_coefficient = sp.diff(height,profile)
        height_hat = height_coefficient*atom('s11cc1_w1_profile_hat_transfer')
        height_kernel = ns['fourier_profiles'](inputs,height_hat,kout,kin)
        normal_coordinate = sp.Symbol('trace_normal_coordinate',real=True)
        ref = values['W_0']/2
        extension = sp.exp(sp.I*qout*(normal_coordinate-ref))
        normal_output = sp.diff(extension,normal_coordinate).subs(normal_coordinate,ref)
        normal_middle,normal_input = normal_output.xreplace({qout:qm}),normal_output.xreplace({qout:qin})
        trace_statement = assignment(c2,'reference_pressure_kernels','trace_three')
        context = dict(sp=sp,value_coefficient=value_coefficient,height_constant=height_constant,
            normal_output=normal_output,normal_middle=normal_middle,normal_input=normal_input,
            height_kernel=height_kernel,left=left,right=right)
        exec(trace_statement,context)
        T = context['trace_three'].subs(endpoint_binding)
        old = T.upper_triangular_solve(C['old'])
        new = T.upper_triangular_solve(C['new'])
        delta_ref = (new-old).applyfunc(sp.cancel)
        normal_matrix = sp.diag(qo_value,qm,qi_value)*sp.I
        delta_jet = normal_matrix*delta_ref
        J.emit('reference-before-guards',dict(nativeTrace=source,amplitudeTrace=trace,
            valueCoefficient=value_coefficient,height=height,constant=constant,
            heightConstant=height_constant,heightCoefficient=height_coefficient,
            nativeMatrixSource=trace_statement,traceOperator=T,oldReference=old,newReference=new,
            deltaReference=delta_ref,normalExtension=extension,normalMatrix=normal_matrix,
            deltaNormalJet=delta_jet,referenceMixed=grade(delta_ref[0,2]),
            jetMixed=grade(delta_jet[0,2])))
        J.zero('native-trace-reconstruction',trace,value_coefficient*p+height*jet+constant)
        J.zero('height-reconstruction',height,height_constant+height_coefficient*profile)
        J.zero('trace-constant',constant,0)
        matzero('trace-difference-identity',T*delta_ref,C['delta'])
        matzero('normal-jet-identity',delta_jet,normal_matrix*delta_ref)
        return dict(trace=T,referenceMixed=grade(delta_ref[0,2]),jetMixed=grade(delta_jet[0,2]))
    T = J.stage('reference-pressure',dict(deltaPhysical=C['delta'],nativeTrace=native['trace'][0]),reference)
    # Direct source factors, not a new slab source or incoming transverse mode.
    source = named(response,'DELTA_P').subs({resolvent:1,z:1},simultaneous=True)/epsilon
    source = source.subs(bindings)
    velocity_atom = atom('s11cc1_V_lab_held_plus')
    velocity_coefficient = sp.diff(source,velocity_atom)
    response_factor = sp.cancel(C['factor']*velocity_coefficient)
    ref_factor = sp.cancel(sp.diff(T['referenceMixed'].subs(tag,1),D)*velocity_coefficient)
    jet_factor = sp.cancel(sp.diff(T['jetMixed'].subs(tag,1),D)*velocity_coefficient)
    J.emit('external-factor-dependence',dict(physical=sorted(response_factor.free_symbols,key=str),
        reference=sorted(ref_factor.free_symbols,key=str),normalJet=sorted(jet_factor.free_symbols,key=str)))
    require(not (response_factor.free_symbols | ref_factor.free_symbols | jet_factor.free_symbols),
            'fixed external factors independent of middle momentum and shape grades')
    # Persist physical factors and exact complex domain evidence before acceptance checks.
    denominators = [1+coefficient*Z[0,0],1+coefficient*Z[2,2]]
    certificates=[]
    for d in denominators:
        real,imag = (sp.simplify(x) for x in sp.expand_complex(d).as_real_imag())
        norm=sp.simplify(real**2+imag**2)
        certificates.append(dict(value=d,real=real,imaginary=imag,normSquared=norm,
                                 positive=norm.is_positive,finite=norm.is_finite))
    factor_norm = sp.simplify(sp.expand_complex(response_factor*sp.conjugate(response_factor)))
    J.emit('physical-factors',dict(source=source,velocityCoefficient=velocity_coefficient,
        physicalPressureFactor=response_factor,referencePressureFactor=ref_factor,
        normalJetFactor=jet_factor,externalDenominators=certificates,factorNormSquared=factor_norm,
        physicalIntegrand=response_factor*bare_integrand,referenceIntegrand=ref_factor*bare_integrand,
        jetIntegrand=jet_factor*bare_integrand,
        integralEvaluated=False,sourceIsTransverseMode=False,fullSlabContraction=False))
    for cert in certificates:
        require(cert['positive'] is True and cert['finite'] is True,'external closure domain')
    depth = sp.Symbol('positive_middle_depth', positive=True)
    middle_certificates = []
    for label, physical_q in [('propagating', depth), ('evanescent', sp.I*depth)]:
        middle_denominator = (1+coefficient*Z[1,1]).subs(qm,physical_q)
        real_part = sp.simplify(sp.re(middle_denominator))
        middle_certificates.append(dict(branch=label,denominator=middle_denominator,
                                        realPart=real_part,positive=real_part.is_positive))
    J.emit('middle-domain',dict(certificates=middle_certificates,
        endpoints='q=0 is excluded from these inverses; saved bare-action endpoint limits remain inherited.'))
    require(all(x['positive'] is True for x in middle_certificates),'native outgoing middle closure domain')
    require(factor_norm.is_finite is True,'finite response factor')
    J.emit('response-status',dict(zero=response_factor.is_zero,
        nonzeroNorm=factor_norm.is_positive,
        note='An exact zero is a possible result, not a failed expected-answer check.'))
    # Addressed direct-slot omission, sign reversal, and one-sided closure corruption.
    omit = C['mixed'].subs(tag,0)
    flipped = C['mixed'].subs(tag,-1)
    one_side = grade(((sp.eye(3)+coefficient*Z).upper_triangular_solve(added))[0,2]).subs(tag,1)
    right = C['mixed'].subs(tag,1)
    movement = sp.cancel(sp.diff(one_side-right,D))
    J.emit('controls',dict(omitted=omit,baseline=right,flipped=flipped,
        oneSided=one_side,missingInputResponseMovement=movement,
        missingInputResponseNorm=sp.simplify(sp.expand_complex(movement*sp.conjugate(movement)))))
    J.zero('omission-control',omit,0)
    J.zero('sign-control',flipped,-right)
    require(sp.simplify(sp.expand_complex(movement*sp.conjugate(movement))).is_positive is True,
            'input-resolvent omission responds')
    # Native grade convention remains rectangular; no pure-second-order completion claimed.
    J.emit('scope',dict(eta=eta,sigma=sigma,epsilon=epsilon,
        outputs='difference in selected physical/reference pressure and normal-jet kernels',
        deferred=['slab-row contraction','two-face cancellations','on-shell transverse excitation',
                  'full second-order shape','production regeneration','defect solve','power or loss'],
        inheritedBareEvidence='saved selected tanh action; no new boundary/profile/integral calculation',
        independentBuildClearance=False,scientificAcceptance=False))
    return dict(operationCount=len(J.completed),candidateEvidenceOnly=True,
                productionChanges=False,integralOrLossValue=False)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    manifest=json.loads(args.inputs.read_text())
    for path,pin in manifest['sourcePins'].items():
        require(sha(path)==pin,'source pin '+path)
    require(manifest['sourcePins'][str(Path(__file__).resolve())]==sha(__file__),'worker pin')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c')
    args.out.mkdir(exist_ok=False)
    started=time.monotonic(); J=None; result={}; code=1
    try:
        save(args.out/'containment.json',containment())
        global sp,Str
        import sympy as sp
        from sympy.core.symbol import Str
        J=Journal(args.out)
        result=science(manifest,J)
        result['executionStatus']='COMPLETED_CANDIDATE_OPERANDS'
        code=0
    except BaseException:
        result=dict(executionStatus='FAILED_PRESERVED',traceback=traceback.format_exc(),
                    incompleteOperation=None if J is None else J.active,automaticRetry=False)
        save(args.out/'failure.json',result)
    finally:
        post={p:dict(expected=h,actual=sha(p)) for p,h in manifest['sourcePins'].items()}
        save(args.out/'posthashes.json',post)
        if any(v['expected']!=v['actual'] for v in post.values()):
            result['integrityFailure']=True; code=1
        save(args.out/'operation-index.json',[] if J is None else J.completed)
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result)
        sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':
    sys.exit(main())
