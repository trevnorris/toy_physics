#!/usr/bin/env python3
"""Candidate upper-face mixed-grade instrument; no producer or production edits.

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


def science(manifest, J):
    scripts = ROOT/'research/pde_ledger_v3/scripts'
    text = (scripts/'S11c_c1_bulk_closure_sympy_audit.py').read_text()
    exports = (scripts/'S11c_c1_exports.py').read_text()
    physical = json.loads(Path(manifest['physicalInput']).read_text())
    parameters = {k: sp.sympify(v) for k, v in physical['parameters'].items()}
    parameters['omega'] = sp.Integer(3)  # Explicit diagnostic frequency, as in approved pilot.
    w0, length, rho, omega, cs = [parameters[k] for k in ('W_0','L_W','rho_m','omega','c_s0')]
    edge = [parameters[k] for k in ('s11cdTangentialMomentum1','s11cdTangentialMomentum2')]
    cutoff_squared = omega**2/cs**2 - sum(k**2 for k in edge)
    Q = sp.Rational(manifest['outputProfileMomentum'])
    q0, qQ = sp.sqrt(cutoff_squared), sp.sqrt(cutoff_squared-Q**2)
    J.emit('physical-point', dict(originalParameters=physical['parameters'], omega=omega,
        edge=edge, inputProfileMomentum=sp.S.Zero, outputProfileMomentum=Q,
        cutoffSquared=cutoff_squared, inputDepth=q0, outputDepth=qQ,
        anchoring='LAB_HELD', density='RHO4_CONSTANT', bulkDrain=0,
        scope='Upper-face bare impedance action; no on-shell slab mode or loss calculation.'))
    require(q0.is_positive is True and qQ.is_positive is True, 'nongrazing incident/output bulk legs')

    def native():
        raw = literal_record(exports, 'dtn_kernel')
        full = sp.sympify(raw, locals={'Str': Str})
        row = next(v for k,v in full if str(k[0]) == 'LAB_HELD' and k[1] == 1)
        labels = {str(k):v for k,v in row}
        parts = {str(k):v for k,v in labels['VALUE']}
        first = parts['FIRST_SHAPE']
        symbols = {s.name:s for s in first.free_symbols}
        kin = tuple(symbols['s11cc1_k_input_'+str(i)] for i in (1,2,3))
        jets = tuple(symbols['s11cc1_w1_profile_jet_hat_'+str(i)] for i in (1,2,3))
        ns = {'sp':sp, 'Mapping':dict, 'Inputs':object, 'DIRECTIONS':(1,2,3),
              'k_in':kin, 'w1_hat':symbols['s11cc1_w1_profile_hat_transfer'], 'w1_jet_hat':jets}
        for name in ('shape_source','dtn_flat_symbol','dtn_first_kernel'):
            exec(compile(extract_function(text,name),str(scripts/'S11c_c1_bulk_closure_sympy_audit.py'),'exec'), ns)
        ctx = SimpleNamespace(**{key:symbols[value] for key,value in
            {'eta':'eta_bg','sigma_W':'sigma_W','W0':'W_0','rho_m':'rho_m','omega':'omega','c_s0':'c_s0'}.items()})
        shape = ns['shape_source'](ctx)
        built = ns['dtn_first_kernel'](ctx,shape,symbols['s11cc1_q_out_output'],symbols['s11cc1_q_out_input'])
        J.zero('native-linear-source-join',built,first)
        z0 = ns['dtn_flat_symbol'](ctx,symbols['s11cc1_q_out_output'])
        diagonal = parts['FLAT_DIAGONAL'].xreplace({d:sp.S.One for d in parts['FLAT_DIAGONAL'].atoms(sp.DiracDelta)})
        J.zero('native-flat-source-join',z0,diagonal)
        return dict(literal=raw, first=first, dimension=labels['DIMENSION_L_T_M'],
                    shape=shape, symbols=symbols, kin=kin, ctx=ctx, nativeFunction=ns['dtn_first_kernel'])

    # The callable/context are ephemeral; the full operands are emitted separately.
    N = native()
    J.emit('native-operands', {k:v for k,v in N.items() if k not in ('ctx','nativeFunction')})
    h,s = sp.symbols('h s', real=True)
    k,H,S = sp.symbols('k H S', real=True)
    qi,qh,qs,qo = sp.symbols('q_i q_h q_s q_o', nonzero=True)
    om,rm = N['ctx'].omega,N['ctx'].rho_m
    amps = sp.symbols('A00 A10 A01 A11')
    modes = [(0,0,k,qi,amps[0]),(1,0,k+H,qh,amps[1]),
             (0,1,k+S,qs,amps[2]),(1,1,k+H+S,qo,amps[3])]

    def boundary():
        geometry=(scripts/'S11c_a_interface_geometry_sympy_audit.py').read_text()
        assignments={}
        for node in ast.walk(ast.parse(geometry)):
            if isinstance(node,ast.Assign):
                for target in node.targets:
                    if isinstance(target,ast.Name) and target.id in ('denominator','normal_exact'):
                        fragment=ast.get_source_segment(geometry,node)
                        if 'grad_h' in fragment:
                            assignments.setdefault(target.id,[]).append(fragment)
        require(all(len(assignments.get(k,[]))==1 for k in ('denominator','normal_exact')),
                'unique native graph normal assignments')
        normal_ns={'sp':sp,'grad_h':(s,sp.S.Zero,sp.S.Zero),'face':1,
                   'dot':lambda a,b:sum(x*y for x,y in zip(a,b))}
        exec('\n'.join(assignments[k][0] for k in ('denominator','normal_exact')),normal_ns)
        normal=sp.Matrix(normal_ns['normal_exact'])
        J.emit('native-graph-normal',dict(source=assignments,normal=normal))
        # Unit-normal normalization has no linear-s term and no h*s term.
        J.zero('normal-through-slope',sp.diff(normal[0],s).subs(s,0),-1)
        J.zero('normal-vertical-through-slope',sp.diff(normal[3],s).subs(s,0),0)
        velocity, pressure = 0,0
        slope_factor=sp.Symbol('native_tangential_normal_multiplier',real=True)
        for a,b,momentum,depth,amplitude in modes:
            for n in range(2-a):
                shifted = h**(a+n)*s**b*amplitude*(sp.I*depth)**n/sp.factorial(n)
                velocity += sp.I*depth*shifted
                pressure += sp.I*om*rm*shifted
                if b == 0:
                    velocity += slope_factor*normal[0].diff(s).subs(s,0)*s*sp.I*momentum*shifted
        equations = [sp.expand(velocity-1).coeff(h,a).coeff(s,b) for a,b,_,_,_ in modes]
        solution = {}
        for equation,amplitude in zip(equations,amps):
            solved = sp.solve(equation.subs(solution),amplitude)
            require(len(solved)==1,'unique triangular boundary amplitude')
            solution[amplitude] = sp.factor(solved[0])
        coefficients = [sp.factor(sp.expand(pressure.subs(solution)).coeff(h,a).coeff(s,b)) for a,b,_,_,_ in modes]
        J.emit('boundary-before-guards', dict(normal=normal, modes=modes, velocity=velocity,
            pressure=pressure, equations=equations, amplitudes={str(a):v for a,v in solution.items()},
            coefficients=coefficients, omittedGrades=['h^2','s^2'],
            note='h and s represent the independent eta and sigma slots; normalizations restored below.'))
        for i,equation in enumerate(equations):
            J.zero('boundary-residual-'+str(i),equation.subs(solution),0)
        epsilon=sp.Symbol('epsilon',real=True)
        J.zero('wave-amplitude-linearity',pressure.subs({a:epsilon*v for a,v in solution.items()}),
               epsilon*pressure.subs(solution))
        # Native linear correspondence retains arbitrary profile-direction k.
        sx=N['symbols']; mapping=dict(zip(N['kin'],(k,*edge)))
        mapping.update({sx['s11cc1_q_out_input']:qi,sx['s11cc1_q_out_output']:qo})
        native_h=N['nativeFunction'](N['ctx'],{'height_hat':1,'tilt_hat':(0,0,0)},qo,qi).subs(mapping)
        native_h=sp.expand(native_h).subs(k**2,om**2/N['ctx'].c_s0**2-sum(e**2 for e in edge)-qi**2)
        native_s=N['nativeFunction'](N['ctx'],{'height_hat':0,'tilt_hat':(1,0,0)},qo,qi).subs(mapping)
        J.zero('native-height-coefficient',coefficients[1].subs(qh,qo),native_h)
        J.zero('native-slope-coefficient',coefficients[2].subs({qs:qo,slope_factor:1}),native_s)
        mixed=coefficients[3].subs(slope_factor,1)
        selected=sp.factor(mixed.subs(k,0))
        J.zero('selected-mixed-reduction',selected,-sp.I*om*rm*H*qi/(qh*qo))
        J.zero('rigid-height-translation',selected.subs(H,0),0)
        J.zero('constant-end-zero-jet',(h*s*mixed).subs(s,0),0)
        return dict(mixed=mixed,selected=selected,height=coefficients[1],slope=coefficients[2],
                    controlCoefficient=coefficients[3],controlMultiplier=slope_factor,
                    amplitudes={str(a):v for a,v in solution.items()})

    B=J.stage('boundary-coefficient',dict(modes=modes,h=h,s=s,omega=om,rho=rm),boundary)

    def profile():
        y=sp.Symbol('y',real=True); L=sp.Symbol('L',positive=True); t=sp.Symbol('t',real=True)
        xi=sp.Symbol('xi',real=True)
        w=sp.sympify(physical['profiles']['w'],locals={'xi':xi}).subs(xi,y/L)
        derivative=sp.diff(w,y); jet=L*derivative
        c2=(scripts/'S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
        coordinates=(y,*sp.symbols('edge_y1 edge_y2',real=True))
        functions={'sp':sp,'re':re,'Y':coordinates,'NEW_DIMENSIONS':{},
                   'DIMENSION_SCHEMA':{'w1_profile':(0,0,0)}}
        node=next(n for n in ast.walk(ast.parse(c2)) if isinstance(n,ast.FunctionDef) and n.name=='at_source')
        exec(compile(textwrap.dedent(ast.get_source_segment(c2,node)),'native-at_source','exec'),functions)
        for name in ('fourier_profiles','profile_bindings'):
            exec(compile(extract_function(c2,name),'native-'+name,'exec'),functions)
        atoms=dict(N['symbols'])
        def atom(name):
            if name not in atoms:atoms[name]=sp.Symbol(name,real=True)
            return atoms[name]
        source_ctx=SimpleNamespace(a=atom,wave=set(),values={'L_W':L})
        source_ctx.at_source=lambda x:functions['at_source'](source_ctx,x)
        original_bindings=functions['profile_bindings'](source_ctx)
        w_function=source_ctx.at_source(atom('w1_profile'))
        field_sub={w_function:w}
        # Do NOT evaluate the Fourier integrals. Bind only their local operands.
        local_profile=source_ctx.at_source(atom('w1_profile')).subs(field_sub)
        local_jet=source_ctx.at_source(atom('w1_profile_d1')).subs(field_sub).doit()
        J.zero('native-profile-operand',local_profile,w)
        J.zero('native-profile-jet-operand',local_jet,jet)
        J.emit('native-fourier-operands',dict(original=original_bindings,profile=local_profile,
            jet=local_jet,coordinates=coordinates,
            edgeReduction='Each unchanged edge integral exp(-i delta_k y)/(2 pi) gives its Dirac delta; the profile transform retains 1/(2 pi).'))
        u=sp.Symbol('u',positive=True)
        logistic=(1+ (u-1)/(u+1))/2
        J.zero('profile-logistic-derivative',sp.diff(logistic,u),1/(1+u)**2)
        J.zero('profile-native-jet',jet,sp.sech(y/L)**2/2)
        J.zero('profile-reflection',w+w.subs(y,-y),1)
        # x=exp(2y/L): w'(y)dy=dx/(1+x)^2. Euler's beta integral
        # with exponent 1-i L t/2 and its conjugate has real part 1.
        z=L*t/2
        beta_product=sp.gamma(1-sp.I*z)*sp.gamma(1+sp.I*z)/(2*sp.pi)
        A=L*t/(4*sp.sinh(sp.pi*L*t/2))
        gamma_reduced=sp.simplify(sp.expand_func(beta_product))
        J.zero('fourier-beta-reflection',gamma_reduced,A)
        J.zero('fourier-derivative',(sp.I*t)*(L/(4*sp.I*sp.sinh(sp.pi*L*t/2))),A)
        atzero=sp.limit(A,t,0)
        J.zero('fourier-derivative-mass',atzero,1/(2*sp.pi))
        return dict(profile=w,derivative=derivative,jet=jet,substitution=sp.exp(2*y/L),
            betaProduct=beta_product,A=A,removableValue=atzero,
            heightTransform='delta(t)/2 + PV[L/(4 i sinh(pi L t/2))]',
            jetTransform=L*A,
            distributionIdentity='t*delta(t)=0; t*PV(1/t)=1; fixes t*w_hat(t)=A(t)/i.',
            convention='Forward exp(-i k y)/(2 pi), inverse exp(+i k y); two conserved edge deltas factored.',
            L=L,t=t)

    P=J.stage('profile-transform',dict(profile=physical['profiles']['w'],nativeConvention='c2 profile_bindings'),profile)

    def action():
        t=P['t']; A=P['A'].subs(P['L'],length)
        At=A; AQ=A.subs(t,Q-t)
        q=sp.Symbol('q_t',nonzero=True)
        hweighted=w0*At/(2*sp.I)  # t*h_hat/eta, including removal of constant delta.
        slope=length*AQ/2        # s_hat/sigma; derivative of delta is not inserted.
        selected=B['selected'].subs({H:t,qi:q0,qh:q,qo:qQ,om:omega,rm:rho})
        # selected/t is polynomially cancelled BEFORE the transfer t=0 is used.
        integrand=sp.factor(sp.cancel(selected/t)*hweighted*slope)
        prefactor=-w0*omega*rho*q0*length/(4*qQ)
        J.zero('native-profile-weighting',integrand,prefactor*At*AQ/q)
        symmetrized=(integrand+integrand.xreplace({t:Q-t,q:sp.Symbol('q_Q_minus_t',nonzero=True)}))/2
        point=Q/2; depth=sp.sqrt(cutoff_squared-point**2)
        sample=sp.simplify(integrand.subs({t:point,q:depth}))
        wrong_sheet=integrand.subs({t:point,q:-depth})
        control_sub={k:sp.S.Zero,H:point,qi:q0,qh:depth,qo:qQ,om:omega,rm:rho}
        control_weight=(hweighted*slope/t).subs(t,point)
        omitted_slope=sp.simplify(B['controlCoefficient'].subs(control_sub).subs(B['controlMultiplier'],0)*control_weight)
        reversed_slope=sp.simplify(B['controlCoefficient'].subs(control_sub).subs(B['controlMultiplier'],-1)*control_weight)
        J.emit('responsive-controls',dict(t=point,physicalDepth=depth,baseline=sample,
            wrongSheet=wrong_sheet,wrongSheetMovement=sp.simplify(wrong_sheet-sample),
            omittedSlope=omitted_slope,omissionMovement=sp.simplify(omitted_slope-sample),
            slopeReversal=reversed_slope,reversalMovement=sp.simplify(reversed_slope-sample)))
        require(sample.is_negative is True,'actual propagating sample sign')
        require(sp.simplify(wrong_sheet-sample).is_zero is False,'sheet control must move')
        require(sp.simplify(omitted_slope-sample).is_zero is False,'addressed native slope omission must move')
        # Whole-interval sign is a factor argument, not a sampled assertion.
        pos=sp.Symbol('positive_transfer',positive=True)
        positive_A=P['A'].subs({P['L']:length,t:pos})
        J.zero('A-even',At.subs(t,-t),At)
        require(positive_A.is_positive is True,'A positive for positive real transfer')
        require(P['removableValue'].is_positive is True,'A positive at removable zero')
        require(prefactor.is_negative is True,'real prefactor sign')
        a=q0
        numerator=sp.simplify(prefactor*At*AQ)
        endpoints=[]
        for end in (-a,a):
            value=sp.simplify(numerator.subs(t,end))
            endpoints.append(dict(t=end,numerator=value,finite=value.is_finite,negative=value.is_negative))
            require(value.is_finite is True and value.is_negative is True,'finite nonzero branch numerator')
        u=sp.Symbol('u',positive=True)
        J.zero('interior-endpoint-square',cutoff_squared-(a-u**2)**2,u**2*(2*a-u**2))
        J.zero('exterior-endpoint-square',(a+u**2)**2-cutoff_squared,u**2*(2*a+u**2))
        tail=sp.limit(At*sp.exp(sp.pi*length*t/2)/t,t,sp.oo)
        J.zero('profile-tail',tail,length/2)
        g=sp.Symbol('positive_depth',positive=True)
        interior=sp.simplify(integrand.subs(q,g))
        exterior=sp.simplify(integrand.subs(q,sp.I*g))
        J.zero('exterior-real-part',sp.re(exterior),0)
        # Full t convolution already includes both height/slope assignments.
        # The explicitly symmetrized form integrates identically by t -> Q-t;
        # it includes a factor 1/2 and has branch sets for BOTH middle legs.
        return dict(integrand=integrand,prefactor=prefactor,A=At,otherA=AQ,
            symmetrized=symmetrized,interior=interior,exterior=exterior,
            branchPoints=[-a,a],inputDispersionResidual=sp.simplify(q0**2-cutoff_squared),
            outputDispersionResidual=sp.simplify(qQ**2+Q**2-cutoff_squared),
            endpoints=endpoints,tailCoefficient=tail,
            signIngredients=dict(A_even=True,A_positive_on_positive_argument=True,
                A_zero=P['removableValue'],prefactor=prefactor),
            integrationNotPerformed=True,
            reasoningForReview='Interior has two positive A factors / positive outgoing depth and negative prefactor; exterior denominator is positive imaginary. Square-root endpoints locally integrable; exponential tails integrable. Delta/PV singularities removed by exact transfer multiplication.',
            scope='Coefficient of independent eta*sigma only, not the complete physical eta^2 shape correction.')

    R=J.stage('profile-action',dict(coefficient=B['selected'],profile=P,physical=parameters),action)
    # Dimensions are bound to native Z0 and the reduced Fourier convention.
    # [rho_m]=M/L^4 in this four-spatial-dimensional model.
    native_dim=N['dimension']
    dimZ=sp.ImmutableMatrix([-3,-1,1]); dimL=sp.ImmutableMatrix([1,0,0])
    # rho*omega, W*L, two A (each dimensionless), 1/q and dt cancel:
    # rho*omega + (W L) = (-2,-1,1). A is the transform of w' and is dimensionless.
    reduced=sp.ImmutableMatrix([-4,-1,1])+2*dimL
    J.emit('dimensions',dict(nativeKernel=native_dim,nativeImpedance=dimZ,
        reducedOneDimensionalKernel=reduced,expectedReducedKernel=dimZ+dimL,
        residual=reduced-(dimZ+dimL),factoredEdgeDeltas=2))
    require(reduced==dimZ+dimL,'reduced kernel dimensions')
    require(all(tuple(d)==(0,-1,1) for d in native_dim),'native three-dimensional kernel dimensions')
    # Exact source fact, not execution of the c2 producer.
    bridge=extract_function((scripts/'S11c_c2_selfenergy_fold_sympy_audit.py').read_text(),'kernel_bridge')
    tree=ast.parse(bridge)
    assignment=next(n for n in ast.walk(tree) if isinstance(n,ast.Assign)
                    and any(isinstance(t,ast.Name) and t.id=='z_three' for t in n.targets))
    slot=assignment.value.args[0].elts[0].elts[2]
    require(isinstance(slot,ast.Constant) and slot.value==0,'native direct slot source')
    J.emit('native-direct-slot',dict(source=ast.get_source_segment(bridge,assignment),literalSlot=0,
        productionChanged=False,closedResolventComputed=False))
    return dict(operationCount=len(J.completed),candidateEvidenceOnly=True,
        independentClearance=False,productionChanges=False,integralOrLossValue=False)


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
