#!/usr/bin/env python3
"""Candidate selected direct-term source and slab-consumer diagnostic.

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
        self.active = name
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



def selected_increment(row, eta, sigma, D, slots, ref_factor, jet_factor):
    """The same selected-slot substitution for the physical and ablated rows."""
    pplus,pminus,jplus,jminus=slots
    increment=row.subs({pplus:eta*sigma*D*ref_factor,jplus:eta*sigma*D*jet_factor,
                        pminus:0,jminus:0},simultaneous=True)
    mixed=sp.cancel(sp.diff(increment,eta,sigma).subs({eta:0,sigma:0},simultaneous=True))
    return increment,mixed


def exact_nonzero(value):
    """Same exact zero predicate; retain unknown and never use a numeric tolerance."""
    direct = value.is_zero
    evidence = dict(value=value, directZero=direct, finite=value.is_finite,
                    decision=None, route='direct-symbolic-flag')
    if value.is_finite is False:
        return evidence
    if direct is not None:
        evidence['decision'] = direct is False
        return evidence
    if value.free_symbols or value.is_number is not True:
        return evidence
    expanded = sp.expand_complex(value)
    real, imaginary = (sp.cancel(x) for x in expanded.as_real_imag())
    evidence.update(route='exact-constant-real-imaginary', expanded=expanded,
                    real=real, imaginary=imaginary,
                    componentZero=[real.is_zero, imaginary.is_zero],
                    componentFinite=[real.is_finite, imaginary.is_finite])
    if real.is_finite is True and imaginary.is_finite is True:
        if real.is_zero is False or imaginary.is_zero is False:
            evidence['decision'] = True
        elif real.is_zero is True and imaginary.is_zero is True:
            evidence['decision'] = False
    return evidence


def forbidden_profile_keys(profiles):
    return sorted(k.name for k in profiles if k.name in ('eta_bg','sigma_W')
                  or k.name.startswith(('w1_profile','m1_profile','gamma_')))


def zero_grade_leftovers(value, profiles):
    names = {k.name for k in profiles} | {'eta_bg','sigma_W'}
    return sorted(a.name for a in value.free_symbols if a.name in names
                  or a.name.startswith(('w1_profile','m1_profile')))


def science(manifest, J):
    """New source/consumer algebra only; prior closure/integral returns are inputs."""
    extractor_path = Path(manifest['sourceExtractor'])
    context = {'__file__': str(extractor_path), '__name__': 'source_metadata'}
    exec(compile(extractor_path.read_text(), str(extractor_path), 'exec'), context)
    records = json.loads(Path(manifest['nativeRecords']).read_text())
    actual = context['build_records']()
    J.emit('native-source-census-join', dict(saved=records, actual=actual, identical=records == actual))
    require(records == actual, 'complete native source/census identity')
    decode = lambda s: sp.sympify(s, locals={'Str': Str})
    mu_pair = decode(records['chemicalSource']['valueConstructorText'])
    density_case = records['geometry']['background_density_map']['cases'][0]
    density = decode(density_case['valueConstructorText'])
    density_context = decode(density_case['caseConstructorText'])
    velocities = {int(r['case'][1]): decode(r['valueConstructorText'])
                  for r in records['geometry']['face_velocity']['cases']}
    responses = {int(r['case'][1]): decode(r['valueConstructorText'])
                 for r in records['faceResponseSources']['cases']}
    row_source = {}
    for r in records['slabConsumers']['rows']:
        address = r['address']
        name = 'U' + str(address[-1]) if address[-2] == 'EXPANDED' else address[-2]
        row_source[name] = sp.Add(*(decode(c['constructorText']) for c in r['selectedPressureJetChildren']))
    objects = [mu_pair, density_context, *velocities.values(), *responses.values(), *row_source.values()]
    symbols = set().union(*(x.atoms(sp.Symbol) for x in objects))
    def atom(name):
        found = [a for a in symbols if a.name == name]
        require(len(found) == 1, 'unique native atom ' + name)
        return found[0]
    eta, sigma, epsilon = (atom(n) for n in ('eta_bg', 'sigma_W', 'epsilon_shape'))
    physical = json.loads(Path(manifest['physicalInput']).read_text())
    values = {k: sp.sympify(v) for k,v in physical['parameters'].items()}
    values['omega'] = sp.Integer(3)
    numeric = {a: values[a.name] for a in symbols if a.name in values
               and a not in (eta,sigma,epsilon)}
    # Use the same actual equality family as Inputs.profiles, keeping sigma independent.
    profiles = {eq.lhs: eq.rhs for eq in density_context.atoms(sp.Equality)
                if isinstance(eq.lhs, sp.Symbol) and eq.lhs != sigma}
    J.emit('profile-key-contract',dict(keys=sorted(a.name for a in profiles),
        forbiddenKeys=forbidden_profile_keys(profiles), independentGrades=[eta,sigma]))
    require(not forbidden_profile_keys(profiles),'independent grade/profile keys are not replaced')
    density_map = {atom('rho_br_bg_rho4_constant'): density[1]}
    def bind(expr):
        expr = expr.subs(density_map, simultaneous=True)
        for _ in range(len(profiles)+2):
            changed = expr.xreplace(profiles)
            if changed == expr:
                break
            expr = changed
        else:
            raise ValueError('unresolved background equality cycle')
        return expr.subs(numeric, simultaneous=True)
    J.emit('binding-context',dict(physicalInput=physical, frequency=values['omega'],
        density=density, densityMap={str(k):v for k,v in density_map.items()},
        profileEqualities={str(k):v for k,v in profiles.items()},
        numeric={str(k):v for k,v in numeric.items()},
        independentGrades=[eta,sigma], epsilon=epsilon))
    def regular(name, expr):
        reduced = sp.cancel(expr)
        numerator, denominator = sp.fraction(reduced)
        den0 = sp.simplify(denominator.subs({eta:0,sigma:0}, simultaneous=True))
        J.emit(name+'-domain',dict(original=expr,reduced=reduced,numerator=numerator,
            denominator=denominator,denominatorAtZero=den0,finite=den0.is_finite,zero=den0.is_zero))
        require(den0.is_finite is True and den0.is_zero is False, name+' regular grade domain')
        return reduced
    def flat(name, expr):
        reduced = regular(name,expr)
        zeroth = reduced.subs({eta:0,sigma:0}, simultaneous=True)
        direct = sp.diff(eta*sigma*reduced,eta,sigma).subs({eta:0,sigma:0}, simultaneous=True)
        J.zero(name+'-grade-product',direct,zeroth)
        J.emit(name+'-grade-split',dict(full=reduced,zeroGrade=zeroth,higherGrades=reduced-zeroth,
            meaning='Grade selection with a regularity certificate; product identity, not an independent native-truncation run.',
            hiddenBackgroundSymbols=sorted([a for a in reduced.free_symbols if a.name in
                ('W_bg','mu_R_bg','rho_br_bg_rho4_constant')], key=str),
            zeroGradeUnresolvedSymbols=zero_grade_leftovers(zeroth,profiles)))
        require(not any(a.name in ('W_bg','mu_R_bg','rho_br_bg_rho4_constant')
                        for a in reduced.free_symbols),'complete background substitution')
        require(not zero_grade_leftovers(zeroth,profiles),'zero grade has no unresolved background/profile/grade symbols')
        return zeroth
    # Load only selected readable returns; never call an earlier function.
    prior = Path(manifest['priorDirectory'])
    restored = {}
    for name in manifest['restoreFiles']:
        raw = json.loads((prior/name).read_text())
        J.emit('restored-'+name.removesuffix('.json'),dict(path=str(prior/name),sha256=sha(prior/name),
            saved=raw,functionCalled=False))
        restored[name] = raw
    F = restored['physical-factors.json']
    restore = lambda x: decode(x['srepr'])
    v_coefficient = restore(F['velocityCoefficient'])
    J.emit('inherited-source-normalization',dict(velocityCoefficient=v_coefficient,
        status='inherited return; not independent revalidation'))
    require(v_coefficient.is_zero is False and v_coefficient.is_finite is True,'source normalization finite/nonzero')
    ref_factor = sp.cancel(restore(F['referencePressurePerUnitV'])/v_coefficient)
    jet_factor = sp.cancel(restore(F['normalJetPerUnitV'])/v_coefficient)
    D = sp.Symbol('inherited_whole_bare_mixed_kernel')
    named = lambda value,key: next(v for k,v in value if str(k)==key)
    sources, source_parts, flats = {}, {}, {}
    mu = bind(mu_pair[1]/epsilon)
    J.emit('native-chemical-amplitude',dict(raw=mu_pair,amplitude=mu))
    mu0 = flat('chemical-amplitude',mu)
    for face in (1,-1):
        label = 'plus' if face == 1 else 'minus'
        response = responses[face]
        rz = named(response,'RESOLVENT')
        dtn = atom('s11cc1_dtn_operator_lab_held_'+label)
        raw = named(response,'DELTA_P').subs({rz:1,dtn:1}, simultaneous=True)/epsilon
        V = atom('s11cc1_V_lab_held_'+label)
        M = atom('s11cc1_mu_theta_lab_held_'+label)
        # Saved and original upper source must agree before a new channel is used.
        if face == 1:
            J.zero('saved-upper-c1-source',bind(raw),bind(restore(F['source'])))
        velocity = bind(velocities[face]/epsilon)
        vc = bind(sp.diff(raw,V)); mc = bind(sp.diff(raw,M))
        if face == 1:
            J.zero('native-saved-velocity-normalization',vc,v_coefficient)
        reconstruction = bind(raw.subs({V:velocity,M:mu},simultaneous=True))
        parts = dict(velocity=vc*velocity,chemical=mc*mu)
        J.emit(label+'-source-input',dict(raw=raw,nativeVelocity=velocities[face],
            velocityAmplitude=velocity,chemicalAmplitude=mu,velocityCoefficient=vc,
            chemicalCoefficient=mc,combined=reconstruction,parts=parts))
        J.zero(label+'-source-linearity',reconstruction,parts['velocity']+parts['chemical'])
        sources[face] = reconstruction
        source_parts[face] = {k:flat(label+'-'+k,v) for k,v in parts.items()}
        flats[face] = flat(label+'-source',reconstruction)
    # Execute only the unchanged native trial-field helpers, never the producer.
    c2 = (ROOT/'research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py').read_text()
    X = tuple(sp.Symbol('s11cc2X'+str(i),real=True) for i in (1,2,3))
    time_symbol = sp.Symbol('s11cc2Time',real=True)
    ns = dict(sp=sp,re=re,X=X,TIME=time_symbol,NEW_DIMENSIONS={})
    for name in ('field','wave_jet'):
        exec(compile(extract_function(c2,name),'<native-'+name+'>','exec'),ns)
    pattern = re.compile(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)')
    def wave_atoms(expr):
        return sorted([a for a in expr.free_symbols if pattern.fullmatch(a.name)
                       or a.name.startswith('grad_theta_')],key=str)
    def restrict(expr, sector):
        return sp.expand(expr.xreplace({a:ns['wave_jet'](a,sector=sector) for a in wave_atoms(expr)}))
    def operator_form(expr):
        from sympy.core.function import AppliedUndef
        jets=sorted(expr.atoms(sp.Derivative,AppliedUndef),key=sp.default_sort_key)
        dummies={j:sp.Dummy('independentJet'+str(i)) for i,j in enumerate(jets)}
        polynomial=expr.xreplace(dummies)
        coefficients={str(j):sp.cancel(sp.diff(polynomial,d)) for j,d in dummies.items()}
        remainder=sp.expand(polynomial-sum(coefficients[str(j)]*d for j,d in dummies.items()))
        require(remainder==0,'linear arbitrary-trial jet reconstruction')
        vals=list(coefficients.values())
        decisions={name:exact_nonzero(v) for name,v in coefficients.items()}
        state='zero' if all(v==0 for v in vals) else (
            'nonzero' if any(v['decision'] is True for v in decisions.values()) else 'unresolved')
        return dict(expression=expr,independentJetCoefficients=coefficients,
                    coefficientNonzeroEvidence=decisions,remainder=remainder,state=state)
    restrictions = {}
    for face,source0 in flats.items():
        label = 'plus' if face == 1 else 'minus'
        restrictions[face]={sector:restrict(source0,sector) for sector in
                            ('TRANSVERSE','LONGITUDINAL','THETA','E_W')}
        J.emit(label+'-source-restrictions',dict(source00=source0,
            gradeScope='Only direct kernel(1,1) times source(0,0) times consumer(0,0); higher source grades are preserved but not restricted.',
            velocity00=source_parts[face]['velocity'],chemical00=source_parts[face]['chemical'],
            restrictions=restrictions[face],transverseParts={k:restrict(v,'TRANSVERSE')
            for k,v in source_parts[face].items()},
            operatorForms={k:operator_form(v) for k,v in restrictions[face].items()},
            ansatz=extract_function(c2,'wave_jet'),onShellModeConstructed=False))
    # All row coefficients precede restriction. No new pressure equation or traction sign.
    pplus, pminus = atom('delta_p_plus'),atom('delta_p_minus')
    jplus, jminus = atom('d_w_delta_p_plus'),atom('d_w_delta_p_minus')
    slots = (pplus,pminus,jplus,jminus)
    coefficients, row_factors = {}, {}
    for name,raw in row_source.items():
        row = bind(raw)
        coeff = {str(s):sp.diff(row,s) for s in slots}
        J.zero(name+'-slot-linearity',row,sum(coeff[str(s)]*s for s in slots))
        coefficients[name] = {str(s):flat(name+'-'+str(s),coeff[str(s)]) for s in slots}
        increment,mixed=selected_increment(row,eta,sigma,D,slots,ref_factor,jet_factor)
        expected = D*(coefficients[name][str(pplus)]*ref_factor+
                      coefficients[name][str(jplus)]*jet_factor)
        J.emit(name+'-consumer-input',dict(raw=raw,bound=row,slotCoefficients=coeff,
            incrementPerCombinedSource=increment,mixedPerCombinedSource=mixed))
        J.zero(name+'-consumer-grade-join',mixed,expected)
        row_factors[name]=sp.cancel(mixed/D)
    J.emit('zero-grade-jet-consumers',dict(
        coefficients={name:coeff[str(jplus)] for name,coeff in coefficients.items()},
        statuses={name:('zero' if coeff[str(jplus)]==0 else
            'nonzero' if coeff[str(jplus)].is_zero is False else 'unresolved')
            for name,coeff in coefficients.items()},
        note='A zero coefficient means this inherited jet factor is unused at the retained mixed grade.'))
    # Selected Fourier amplitudes are source at kin and consumer at kout, not products
    # of spatially varying coefficients disguised as a Fourier convolution.
    point = json.loads(Path(manifest['physicalPoint']).read_text())
    edge = tuple(restore(x) for x in point['edge'])
    kin = (restore(point['inputProfileMomentum']),*edge)
    kout = (restore(point['outputProfileMomentum']),*edge)
    J.zero('inherited-frequency',restore(point['omega']),values['omega'])
    require(point['originalParameters']==physical['parameters'],'physical input inherited identity')
    amplitudes={name:sp.Symbol('amplitude_'+name) for name in ('u_1','u_2','u_3','theta','e_W')}
    def plane(expr,momentum):
        mapping={}
        for a in wave_atoms(expr):
            name=a.name
            if name.startswith('grad_theta_'):name='theta_d'+name.rsplit('_',1)[1]
            match=pattern.fullmatch(name)
            require(match is not None,'native plane-wave jet '+name)
            base,suffix=match.groups(); value=amplitudes[base]
            for axis in re.findall(r'd([123])',suffix):value*=sp.I*momentum[int(axis)-1]
            value*=(-sp.I*values['omega'])**(2 if '_tt' in suffix else 1 if '_t' in suffix else 0)
            mapping[a]=value
        return sp.expand(expr.xreplace(mapping))
    source_plane=plane(flats[1],kin)
    require(epsilon not in flats[1].free_symbols and epsilon not in flats[-1].free_symbols,
            'native source amplitudes normalized exactly once')
    outputs={name:sp.cancel(factor*source_plane) for name,factor in row_factors.items()}
    source_channels={name:sp.diff(source_plane,a) for name,a in amplitudes.items()}
    J.zero('source-amplitude-linearity',source_plane,sum(source_channels[n]*a for n,a in amplitudes.items()))
    J.emit('selected-fourier-contraction',dict(kin=kin,kout=kout,omega=values['omega'],
        sourceAmplitude=source_plane,sourceChannels=source_channels,
        rowFactorsPerCombinedSource=row_factors,rowAmplitudesPerWholeKernel=outputs,
        fullMixedIncrement='eta*sigma_W times whole reduced kernel D times row amplitude; no new integral',
        epsilonConvention='native row coefficients keep their epsilon; normalized source carries no epsilon',
        incidentTransverseMode=False,lowerFaceCorrectionConstructed=False))
    # Reverse route: report native slot absence separately from the form-only curl.
    H=sp.Function('selectedKernelAppliedSource')(*X,time_symbol)
    vector=[row_factors['U'+str(i)]*H for i in range(3)]
    curl=[sp.diff(vector[(i+2)%3],X[(i+1)%3])-sp.diff(vector[(i+1)%3],X[(i+2)%3]) for i in range(3)]
    J.emit('weak-directions',dict(vectorPerSourceImage=vector,reverseCurl=curl,
        reverseEvidence='Native pressure/jet consumer census; H/curl is a form-only routing test, not a solved reverse mode.',
        nativeUCensus=[dict(address=r['address'],counts=r['pressureSymbolOccurrenceCounts'],
            selectedChildren=r['selectedPressureJetChildren'])
            for r in records['slabConsumers']['rows'] if r['address'][-2]=='EXPANDED'],
        nativeGeneralizedForceU=records['physicalGeneralizedForceULiteral'],
        forwardTransverseSource=restrictions[1]['TRANSVERSE'],
        forwardAction=({name:sp.S.Zero for name in row_factors}
            if operator_form(restrictions[1]['TRANSVERSE'])['state']=='zero'
            else 'UNRESOLVED: nonzero/unknown source needs its actual kernel action; saved off-shell witness is insufficient'),
        note='Source restriction precedes any kernel action; zero does not imply whole-operator decoupling.'))
    # Responsive controls: native pieces removed, not expected answers inserted.
    vel0=source_parts[1]['velocity']; chem0=source_parts[1]['chemical']
    e_t=atom('e_W_t')
    velocity_piece=sp.diff(vel0,e_t)*e_t
    velocity_movement=sp.cancel(plane(velocity_piece,kin).subs(amplitudes['e_W'],1))
    theta_movement=sp.cancel(sp.diff(plane(chem0,kin),amplitudes['theta']))
    choices=[]
    for jet in wave_atoms(chem0):
        if jet.name.startswith('u_'):
            piece=sp.diff(chem0,jet)*jet
            movement=restrict(piece,'TRANSVERSE')
            if operator_form(movement)['state']=='nonzero':
                choices.append((jet,piece,movement))
    J.emit('source-controls',dict(velocityAddress=str(e_t),removedVelocityPiece=velocity_piece,
        velocityMovement=velocity_movement,removedChemicalChannel=chem0,thetaMovement=theta_movement,
        divergenceOmissions=[dict(jet=j,piece=p,curlSourceMovement=m,
            operatorForm=operator_form(m)) for j,p,m in choices]))
    # Carry the actual omissions through the actual scalar consumers. The full
    # row epsilon is saved; divide it only in explicitly labeled amplitude checks.
    scalar_rows=('THETA_BALANCE','E_W_BALANCE')
    J.emit('end-to-end-controls-input',dict(rowFactors=row_factors,upperSource00=flats[1],
        sourcePlane=source_plane,velocityPiece=velocity_piece,chemicalChannel=chem0,
        divergenceChoices=[dict(jet=j,piece=p,restrictedPiece=r) for j,p,r in choices],
        scalarConsumers={name:row_source[name] for name in scalar_rows},
        upperPressureSlot=pplus,referenceFactor=ref_factor,jetFactor=jet_factor,
        grades=[eta,sigma],epsilon=epsilon))
    downstream=[]
    for name in scalar_rows:
        factor=row_factors[name]
        for channel,source_change in [('velocity',velocity_movement),('chemical',theta_movement)]:
            selected_amplitude='e_W' if channel=='velocity' else 'theta'
            probe={a:sp.Integer(1 if n==selected_amplitude else 0) for n,a in amplitudes.items()}
            original_source=sp.cancel(source_plane.subs(probe,simultaneous=True))
            damaged_source=sp.cancel(original_source-source_change)
            before=sp.cancel(factor*original_source)
            after=sp.cancel(factor*damaged_source)
            row_change=sp.cancel(after-before)
            amplitude=sp.cancel(row_change/epsilon)
            downstream.append(dict(row=name,sourceChannel=channel,consumerFactor=factor,
                removedSourceContribution=source_change,originalSource=original_source,
                ablatedSource=damaged_source,baselineRowPerD=before,ablatedRowPerD=after,
                rowMovementPerD=row_change,rowAmplitudeMovementPerD=amplitude,
                amplitudeFreeOfEpsilon=epsilon not in amplitude.free_symbols,
                nonzeroEvidence=exact_nonzero(amplitude),
                nonzero=exact_nonzero(amplitude)['decision'] is True))
        for jet,piece,movement in choices:
            damaged_source=flats[1]-piece
            damaged_restriction=restrict(damaged_source,'TRANSVERSE')
            before=sp.expand(factor*restrictions[1]['TRANSVERSE'])
            after=sp.expand(factor*damaged_restriction)
            row_change=sp.expand(after-before)
            amplitude=sp.expand(row_change/epsilon)
            downstream.append(dict(row=name,sourceChannel='divergence-omission',jet=jet,
                omittedPiece=piece,originalSource=flats[1],ablatedSource=damaged_source,
                baselineRestriction=restrictions[1]['TRANSVERSE'],ablatedRestriction=damaged_restriction,
                consumerFactor=factor,removedCurlSourceContribution=movement,
                baselineRowPerD=before,ablatedRowPerD=after,
                rowMovementPerD=row_change,rowAmplitudeMovementPerD=amplitude,
                amplitudeFreeOfEpsilon=epsilon not in amplitude.free_symbols,
                operatorForm=operator_form(amplitude)))
    J.emit('source-omissions-through-consumers',dict(
        normalization='Full row movements keep epsilon; amplitude movements divide that one native row epsilon. D is the inherited whole mixed kernel.',
        evidenceRole='Sensitivity/wiring controls, not independent correctness proofs; divergence omission tests the same curl restriction path.',
        records=downstream))
    # Remove an addressed native upper pressure slot, then use the same selected
    # substitution and mixed-coefficient routine as the unmodified rows.
    probe={a:sp.Integer(1 if name=='e_W' else 0) for name,a in amplitudes.items()}
    probe_source=sp.cancel(source_plane.subs(probe,simultaneous=True))
    consumer_controls=[]
    for name in scalar_rows:
        raw=row_source[name]
        damaged_raw=raw.subs(pplus,0)
        damaged_bound=bind(damaged_raw)
        damaged_increment,damaged_mixed=selected_increment(
            damaged_bound,eta,sigma,D,slots,ref_factor,jet_factor)
        damaged_factor=sp.cancel(damaged_mixed/D)
        before=sp.cancel(row_factors[name]*probe_source)
        after=sp.cancel(damaged_factor*probe_source)
        change=sp.cancel((after-before)/epsilon)
        consumer_controls.append(dict(row=name,removedSlot=pplus,native=raw,
            ablatedNative=damaged_raw,bound=damaged_bound,probeAmplitudes={str(k):v for k,v in probe.items()},
            sourceAmplitude=probe_source,selectedIncrement=damaged_increment,
            mixedPerSource=damaged_mixed,baselineRowPerD=before,ablatedRowPerD=after,
            rowMovementPerD=after-before,
            amplitudeMovementPerD=change,amplitudeFreeOfEpsilon=epsilon not in change.free_symbols,
            nonzeroEvidence=exact_nonzero(change),
            nonzero=exact_nonzero(change)['decision'] is True))
    J.emit('native-pressure-consumer-omission',consumer_controls)
    # A misplaced scalar consumer tests the actual reverse curl implementation.
    routed=row_factors['E_W_BALANCE']*H
    wrong_vector=[routed,0,0]
    wrong_curl=[sp.diff(wrong_vector[(i+2)%3],X[(i+1)%3])-
                sp.diff(wrong_vector[(i+1)%3],X[(i+2)%3]) for i in range(3)]
    J.emit('routing-control',dict(actualVector=vector,wrongVector=wrong_vector,
        actualCurl=curl,wrongCurl=wrong_curl,movement=[sp.expand(a-b) for a,b in zip(wrong_curl,curl)],
        role='FORM misrouting control, not another physical operator'))
    # Infer consumer dimensions before binding physical numbers. Source dimensions
    # are the supplied native amplitude declarations, not a fresh energy derivation.
    schema_node=next(n for n in ast.parse(c2).body if isinstance(n,ast.Assign) and
        any(isinstance(t,ast.Name) and t.id=='DIMENSION_SCHEMA' for t in n.targets))
    schema=ast.literal_eval(schema_node.value)
    def dimension(expr):
        if expr.is_number:return sp.zeros(3,1)
        if isinstance(expr,sp.Symbol):
            require(expr.name in schema,'native dimension '+expr.name)
            return sp.Matrix(schema[expr.name])
        if isinstance(expr,sp.Add):
            vals=[dimension(a) for a in expr.args]
            require(all(v==vals[0] for v in vals),'homogeneous native sum')
            return vals[0]
        if isinstance(expr,sp.Mul):
            return sum((dimension(a) for a in expr.args),sp.zeros(3,1))
        if isinstance(expr,sp.Pow) and expr.exp.is_number:return dimension(expr.base)*expr.exp
        raise ValueError('unsupported dimension constructor '+sp.srepr(expr))
    unit_rows=[]
    targets={'THETA_BALANCE':sp.Matrix([-3,-1,1]),'E_W_BALANCE':sp.Matrix([-1,-2,1])}
    for name,raw in row_source.items():
        for slot in slots:
            coeff=sp.diff(raw,slot)
            if coeff!=0:
                total=dimension(coeff)+dimension(slot)
                unit_rows.append(dict(row=name,slot=slot,coefficient=coeff,
                    coefficientDimension=dimension(coeff),total=total,expected=targets.get(name)))
    J.emit('consumer-unit-joins',unit_rows)
    J.emit('units-and-scope',dict(nativeChemical=records['chemicalSource']['caseConstructorText'],
        nativeConsumerDimensions=records['slabConsumers']['caseDimensionsConstructorText'],
        inheritedReducedKernelDimension=[-2,-1,1],normalizedSourceDimension=[1,-1,0],
        reducedPressureFourierDimension=[-1,-2,1],
        declarationStatus='Source/kernel dimensions and Fourier measure convention are inherited declarations; consumer dimensions above are computed against supplied expected row dimensions.',
        note='Two conserved-edge deltas are factored in inherited D. The remaining input-plane measure convention is inherited, not a new integral verification; no delta(0) is evaluated.',
        scopes=['selected upper-face direct retained addition','arbitrary-curl flat source restriction',
                'selected off-shell scalar/longitudinal Fourier witness'],
        deferred=['full finite inverse','lower-face correction and cancellations','on-shell scattering',
                  'production repair','loss','primitive calibration','drain flow','defect sweep']))
    # All available control and unit evidence precedes responsiveness decisions.
    require(velocity_movement.is_zero is False and theta_movement.is_zero is False,'native scalar source controls respond')
    require(bool(choices),'native divergence-piece omission is applicable and responds')
    require(all(r['amplitudeFreeOfEpsilon'] and
                (r['nonzero'] if 'nonzero' in r else r['operatorForm']['state']=='nonzero')
                for r in downstream),'addressed source omissions reach actual scalar consumers')
    require(all(r['amplitudeFreeOfEpsilon'] and r['nonzero'] for r in consumer_controls),
            'native pressure consumer omission responds')
    require(any(sp.expand(a-b)!=0 for a,b in zip(wrong_curl,curl)),'wrong scalar-to-vector route responds')
    require(all(x['expected'] is not None and x['total']==x['expected'] for x in unit_rows),
            'native consumer dimensions')
    return dict(candidateEvidenceOnly=True,sourceAndConsumerCheck=True,productionChanges=False,
        priorClosureOrIntegralReplayed=False,integralOrLossValue=False)


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
        result=J.stage('source-consumer',dict(sourcePins=manifest['sourcePins'],
            priorDirectory=manifest['priorDirectory'],restoreFiles=manifest['restoreFiles']),
            lambda:science(manifest,J))
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
        save(args.out/'artifact-index.json',[
            dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size)
            for p in sorted(args.out.iterdir()) if p.is_file()])
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result)
        sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':
    sys.exit(main())
