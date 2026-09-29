#!/usr/bin/env python3
"""Unexecuted transverse-end/REFERENCE-face build candidate.

No producer import, thickness spectrum, permeability continuation, forced
response, Born coefficient, or power-loss calculation. A new explicit science
approval after the costed-plan stop is required even after build clearance.
"""
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time
import tempfile
import traceback

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
ENDS = ('REFERENCE', 'LEFT', 'RIGHT')
DRIVES = ('AMPLITUDE', 'PRESSURE', 'OUTWARD_VELOCITY', 'RELATIVE_MASS_FLUX',
          'AFFINITY', 'BULK_VELOCITY')


def require(value, reason):
    if not value:
        raise ValueError(reason)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def route(path):
    path = Path(path)
    return dict(path=str(path), canonicalPath=str(path.resolve(strict=True)),
                bytes=path.stat().st_size, sha256=sha(path))


class IntegrityError(RuntimeError):
    """A source or stored-object pin changed; never a local math fallback."""


def publish(path, payload):
    """Publish a complete fsynced byte string with no overwrite window."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.'+path.name+'.', suffix='.partial', dir=path.parent)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(payload); stream.flush(); os.fsync(stream.fileno())
    os.link(temporary, path)  # Atomic exclusive publication; existing path fails.
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)
    os.unlink(temporary)


def save(path, value):
    publish(path, (json.dumps(value, indent=2, allow_nan=False)+'\n').encode())


class OperationBudget(BaseException):
    """Deliberately escapes ordinary algebra/domain exception handlers."""


class NativeDeadline(BaseException):
    """Fatal whole-job deadline; disarm before failure bookkeeping."""


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.')) or
                module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError((module, name))
        return super().find_class(module, name)


class Journal:
    """Complete operands, atomic content-addressed values, local math failures."""
    def __init__(self, out, deadline):
        self.out, self.deadline = out, deadline
        self.records, self.active, self.objects = [], None, {}
        def timeout(*_):
            signal.setitimer(signal.ITIMER_REAL, 0)
            if time.monotonic() >= self.deadline:
                raise NativeDeadline('native whole-job deadline')
            raise OperationBudget('bounded mathematical operation')
        signal.signal(signal.SIGALRM, timeout)
        self.arm_native()

    def arm_native(self):
        signal.setitimer(signal.ITIMER_REAL, 0)
        remaining = self.deadline-time.monotonic()
        if remaining <= 0:
            raise NativeDeadline('native whole-job deadline')
        signal.setitimer(signal.ITIMER_REAL, remaining)

    def blob(self, value):
        existing = self.objects.get(id(value))
        if existing is not None and existing[0] is value:
            record = existing[1]
            if route(record['path']) != record:
                raise IntegrityError('previous object changed')
            return record
        payload = pickle.dumps(value, protocol=4)
        digest = hashlib.sha256(payload).hexdigest()
        path = self.out/'objects'/(digest+'.pickle')
        if path.exists():
            if path.stat().st_size != len(payload) or sha(path) != digest:
                raise IntegrityError('existing content-addressed object is incomplete or changed')
        else:
            publish(path, payload)
        record = route(path)
        if record['sha256'] != digest:
            raise IntegrityError('published object hash mismatch')
        self.objects[id(value)] = (value, record)
        return record

    def op(self, name, function, *args, seconds=20, required=False):
        path = self.out/'operations'/('%04d-%s' % (len(self.records), name))
        path.mkdir(parents=True, exist_ok=False)
        self.active = dict(name=name, operands=[], startedUtc=datetime.now(timezone.utc).isoformat())
        save(path/'started.json', self.active)
        for arg in args:
            self.active['operands'].append(self.blob(arg))
        save(path/'input-references.json', self.active)
        remaining = self.deadline-time.monotonic()
        if remaining <= 0:
            signal.setitimer(signal.ITIMER_REAL, 0)
            raise NativeDeadline('native whole-job deadline')
        signal.setitimer(signal.ITIMER_REAL, min(seconds, remaining))
        failed = None
        try:
            value = function(*args)
            signal.setitimer(signal.ITIMER_REAL, 0)
        except NativeDeadline:
            signal.setitimer(signal.ITIMER_REAL, 0)
            raise
        except OperationBudget:
            signal.setitimer(signal.ITIMER_REAL, 0)
            if required:
                signal.setitimer(signal.ITIMER_REAL, 0)
                raise
            failed = dict(status='UNRESOLVED', reason='OPERATION_BUDGET', traceback=traceback.format_exc())
        except Exception as error:
            signal.setitimer(signal.ITIMER_REAL, 0)
            if required or isinstance(error, (IntegrityError, OSError, pickle.UnpicklingError)):
                signal.setitimer(signal.ITIMER_REAL, 0)
                raise
            failed = dict(status='UNRESOLVED', reason='LOCAL_OPERATION_EXCEPTION',
                          exceptionType=type(error).__name__, traceback=traceback.format_exc())
        # A short operation timer never covers serialization or error receipts.
        # Only the remaining whole-job deadline is armed here, with no 1ms retry.
        if failed is not None:
            self.active['localFailureBeforePersistence'] = failed
        self.arm_native()
        if failed is not None:
            failed['operationReceipt'] = str(path/'unresolved.json')
            record = dict(self.active, status='UNRESOLVED_NO_RETRY', failure=failed)
            save(path/'unresolved.json', record)
            self.records.append(record); self.active = None
            return failed
        record = dict(self.active, status='COMPLETE', value=self.blob(value))
        save(path/'complete.json', record)
        self.records.append(record); self.active = None
        return value


def containment():
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines()
                 if x.startswith('0::'))
    cgroup = Path('/sys/fs/cgroup')/group.lstrip('/')
    limits = {key: (cgroup/key).read_text().strip()
              for key in ('memory.max', 'memory.swap.max', 'pids.max')}
    limits.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={key: os.environ.get(key) for key in THREADS})
    require(limits['memory.max'] == '2147483648' and limits['memory.swap.max'] == '0'
            and limits['pids.max'] == '32' and limits['nice'] >= 15
            and len(limits['affinity']) == 1
            and all(x == '1' for x in limits['threads'].values()),
            'ordinary containment missing; science not imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise OperationBudget('bounded operation/native deadline')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 840)
    return limits


def construct(spec, out, journal):
    import sympy as sp
    import numpy as np

    w = sp.Symbol('premiseFrequency', positive=True)
    v = sp.Symbol('premiseRayCoordinate', nonnegative=True)
    k = sp.Symbol('premiseNormalMomentum', real=True)
    q, qb = sp.symbols('premiseDepthMomentum premiseLeftDepthMomentum', complex=True)
    params = {name: sp.Rational(value) for name, value in spec['physicalInput']['parameters'].items()}

    def show(value):
        if isinstance(value, dict): return {str(a): show(b) for a, b in value.items()}
        if isinstance(value, (list, tuple)): return [show(x) for x in value]
        if isinstance(value, sp.MatrixBase): return show(value.tolist())
        if isinstance(value, np.ndarray): return show(value.tolist())
        if isinstance(value, complex): return dict(real=value.real, imag=value.imag)
        if isinstance(value, (np.integer, np.floating, np.bool_)): return value.item()
        if isinstance(value, sp.Basic): return str(value)
        return value

    def op(name, function, *args, seconds=20, required=False, render=True):
        value = journal.op(name, function, *args, seconds=seconds, required=required)
        if render:
            save(out/(name+'.json'), show(value))
        else:
            save(out/(name+'.json'),dict(operation=name,completeObjectReceipt=journal.records[-1],
                 wholeRestoredPayloadRendered=False))
        return value

    def restored(name, key):
        return op(name, restore, key, required=True, render=False)

    def restore(name):
        record = spec['inputs'][name]
        if route(record['path']) != record: raise IntegrityError('changed source '+name)
        with Path(record['path']).open('rb') as stream:
            return SavedCodec(stream).load()

    def flat(value):
        if isinstance(value, sp.MatrixBase): return list(value)
        if isinstance(value, dict): return [z for x in value.values() for z in flat(x)]
        if isinstance(value, (tuple, list)): return [z for x in value for z in flat(x)]
        return [value]

    def zero(value): return all(sp.cancel(x) == 0 for x in flat(value))
    def clean(matrix): return sp.ImmutableMatrix(matrix).applyfunc(sp.cancel)

    def raw_denominators(values):
        return tuple(dict.fromkeys(
            [sp.fraction(sp.together(value))[1] for value in values]+
            [power.base for value in values for power in value.atoms(sp.Pow)
             if power.exp.is_negative is True]))

    def bind(value, label, extra=None):
        mapping = dict(extra or {})
        # Endpoints come from the accepted saved packet, never re-evaluated.
        value = value.xreplace(uniform['profileEndpoints'])
        for symbol in value.free_symbols-set(mapping):
            name = symbol.name
            if name in ('omega', 's11cdFrequency'): mapping[symbol] = w
            elif name == 's11cdTangentialMomentum1': mapping[symbol] = 2*v
            elif name == 's11cdTangentialMomentum2': mapping[symbol] = v
            elif name == 's11cdSpectralNormalMomentum': mapping[symbol] = k
            elif name == 'eta_bg': mapping[symbol] = 0 if label == 'REFERENCE' else params[name]
            elif name == 'sigma_W':
                mapping[symbol] = 0 if label == 'REFERENCE' else params['eta_bg']*params['W_0']/params['L_W']
            elif name in params: mapping[symbol] = params[name]
            elif symbol not in (w, v, k, q, qb): raise ValueError('unsupported source symbol '+name)
        bound = value.xreplace(mapping)
        require(not bound.atoms(sp.Integral, sp.Limit, sp.Derivative), 'unresolved source atom')
        return bound

    def grade(matrix):
        epsilon = [s for s in matrix.free_symbols if s.name == 'epsilon_shape']
        if not epsilon and zero(matrix): return matrix, matrix
        require(len(epsilon) == 1, 'quadratic current amplitude carrier')
        e = epsilon[0]
        coefficient = matrix.applyfunc(lambda x: sp.Poly(x, e).nth(2))
        return coefficient, matrix-e**2*coefficient

    uniform = restored('restore-uniform','uniform')
    reduction = restored('restore-branch-context','reduction')
    units = restored('restore-unit-context','unitContext')
    common = restored('restore-uniform-common','uniformCommon')
    require(all(isinstance(x, dict) and x.get('reason') != 'OPERATION_BUDGET'
                for x in (uniform, reduction, units)), 'required source restoration incomplete')
    save(out/'unit-context.json', show(dict(unitFrame=spec['physicalInput']['unit_frame'],
         fields=units['fieldUnits'], current=units['currentUnit'], strong=uniform['units']['strong'])))
    face_index = json.loads(Path(spec['inputs']['faceIndex']['path']).read_text())
    source_tree = ast.parse(Path(spec['inputs']['sourceEngine']['path']).read_text())
    acoustic_class = next(node for node in source_tree.body
                          if isinstance(node,ast.ClassDef) and node.name=='ClosedAcousticEnergy')
    acoustic_method = next(node for node in acoustic_class.body
                           if isinstance(node,ast.FunctionDef) and node.name=='construct')
    sign_loops = [node for node in ast.walk(acoustic_method)
                  if isinstance(node,ast.For) and isinstance(node.target,ast.Name)
                  and node.target.id=='sign' and isinstance(node.iter,ast.Tuple)]
    require(len(sign_loops)==1, 'source face-order schema')
    face_order = ast.literal_eval(sign_loops[0].iter)
    law_nodes = {name:next(node for node in acoustic_method.body if isinstance(node,ast.Assign)
                 and any(isinstance(target,ast.Name) and target.id==name for target in node.targets))
                 for name in ('pressure','velocity','local_flux','harmonic')}
    acoustic_law_provenance = {name:dict(line=node.lineno,assignment=ast.unparse(node)) for name,node in law_nodes.items()}
    save(out/'acoustic-law-source.json',dict(source=spec['inputs']['sourceEngine'],laws=acoustic_law_provenance))
    save(out/'native-face-order-source.json',dict(order=face_order,
        source=spec['inputs']['sourceEngine'],sourceLine=sign_loops[0].lineno,
        expression=ast.unparse(sign_loops[0].iter)))

    def branch_inventory(state):
        equations, mapping = state['branch_equations'], state['branch_map']
        residuals = [equation.rhs-mapping[equation.lhs] for equation in equations]
        return dict(equations=equations, mapping=mapping, residuals=residuals,
                    allEquationKeysPresent=all(eq.lhs in mapping for eq in equations),
                    joins=zero(residuals),scope='SUPPLIED_BRANCH_INVENTORY_INTEGRITY_NOT_INDEPENDENT_DERIVATION', tangents=state['tangents'],
                    groups=state['momentum_groups'], normalMap=state['normal_map'])

    branches = op('source-branch-map-joins', branch_inventory, reduction['reductionState'])
    require(branches.get('joins') is True, 'saved branch equation/map inconsistency')

    def carrier_law_check(label, result, mapping, depth):
        time_coordinate, depth_coordinate = sp.symbols('carrierTime carrierOutwardDepth',real=True)
        harmonic = sp.Symbol('carrierHarmonic',nonzero=True)
        epsilon = sp.Symbol('carrierRealFieldScale',real=True)
        minus_amplitude,plus_amplitude = sp.symbols('carrierConjugateAmplitude carrierAmplitude')
        plus = plus_amplitude*harmonic*sp.exp(sp.I*(q*depth_coordinate-w*time_coordinate))
        minus = minus_amplitude/harmonic*sp.exp(-sp.I*(qb*depth_coordinate-w*time_coordinate))
        real_ansatz = epsilon*(plus+minus)/2
        # Source-native pressure and outward velocity laws, applied to one
        # real harmonic carrier pair, not independent counterpropagating waves.
        pressure = -params['rho_m']*sp.diff(real_ansatz,time_coordinate)
        velocity = sp.diff(real_ansatz,depth_coordinate)
        product = sp.expand(pressure*velocity)
        average = product.coeff(harmonic,0)
        normalization = epsilon**2*minus_amplitude*plus_amplitude
        normalized_flux = sp.cancel(average/normalization).subs(depth_coordinate,0)
        pressure_scale = sp.cancel(-params['rho_m']*sp.diff(plus,time_coordinate)/plus)
        velocity_scale = sp.cancel(sp.diff(plus,depth_coordinate)/plus)
        saved_pressure = bind(result['OPEN_BULK_PRESSURE_SCALES'][0],label,mapping)
        saved_velocity = bind(result['OPEN_BULK_VELOCITY_SCALES'][0],label,mapping)
        minus_pressure_scale = sp.cancel(-params['rho_m']*sp.diff(minus,time_coordinate)/minus)
        minus_velocity_scale = sp.cancel(sp.diff(minus,depth_coordinate)/minus)
        saved_minus_pressure = bind(result['OPEN_BULK_PRESSURE_SCALES'][1],label,mapping)
        saved_minus_velocity = bind(result['OPEN_BULK_VELOCITY_SCALES'][1],label,mapping)
        residuals = dict(pressure=sp.cancel(saved_pressure-pressure_scale),
                         velocity=sp.cancel(saved_velocity-velocity_scale),
                         minusPressure=sp.cancel(saved_minus_pressure-minus_pressure_scale),
                         minusVelocity=sp.cancel(saved_minus_velocity-minus_velocity_scale),
                         depthFlux=sp.cancel(depth-normalized_flux))
        return dict(scope='SOURCE_LAW_AND_REAL_CARRIER_NORMALIZATION_CHECK_NOT_INDEPENDENT_CLOSURE',
            provenance=acoustic_law_provenance,plus=plus,conjugateLeg=minus,
            realAnsatz=real_ansatz,pressure=pressure,velocity=velocity,product=product,
            zeroHarmonicAverage=average,amplitudeNormalization=normalization,
            normalizedFlux=normalized_flux,pressureScale=pressure_scale,velocityScale=velocity_scale,
            savedPressure=saved_pressure,savedVelocity=saved_velocity,savedFlux=depth,
            minusPressureScale=minus_pressure_scale,minusVelocityScale=minus_velocity_scale,
            savedMinusPressure=saved_minus_pressure,savedMinusVelocity=saved_minus_velocity,
            residuals=residuals,joins=zero(residuals))

    def bound_end(label, packet, pair_tuple):
        pair, dimensions = pair_tuple  # Authoritative end_pairing_check.py:165-166.
        result = pair['result']
        wl, wr = result['FREQUENCY_LEGS']; kl, kr = result['NORMAL_LEGS']
        ql, qr = result['BULK_LEGS']
        mapping = {wl:w, wr:w, kl:k, kr:k, ql:qb, qr:q}
        raw_matrix = bind(result['CLOSED_PENCIL_LEGS'][0], label, mapping)
        matrix = clean(raw_matrix)
        wave = sp.expand(bind(result['ACOUSTIC_WAVE_ROWS'][0], label, mapping))
        wp = sp.Poly(wave, q)
        require(wp.degree() == 2 and wp.nth(1) == 0, 'unsupported depth equation')
        q2 = sp.cancel(-wp.nth(0)/wp.nth(2))
        oldq = packet['radical']
        original = bind(packet['originalAlgebraic'], label, {oldq:oldq})
        relation = bind(packet['originalRelation'], label, {oldq:oldq})
        scale = sp.sqrt(sp.cancel(sp.Poly(relation, oldq).nth(2)/wp.nth(2)))
        source_join = clean(matrix.subs(q, scale*oldq)-original)
        wave_join = sp.cancel(wave.subs(q, scale*oldq)-relation)
        curl = bind(uniform['curl'], label)
        common_curl_join = clean(bind(common['curl'],label)-curl)
        gauge = tuple(bind(value,label) for value in common['gauge'])
        longitudinal_probes = tuple(sp.ImmutableMatrix.vstack(value,sp.zeros(2,1)) for value in gauge)
        lift = sp.ImmutableMatrix.vstack(curl[:, 1:3], sp.zeros(2, 2))
        raw_gram = lift.H*lift
        gram = clean(raw_gram)
        # A 2x2 coordinate inverse only; never invert the strong physical pencil.
        raw_weak = lift.H*raw_matrix*lift
        raw_restricted = gram.inv()*raw_weak
        weak, restricted = clean(raw_weak), clean(raw_restricted)
        raw_t_denominators = tuple(dict.fromkeys(
            [sp.fraction(sp.together(x))[1] for x in (*raw_weak,*raw_restricted)] +
            [power.base for x in (*raw_weak,*raw_restricted) for power in x.atoms(sp.Pow)
             if power.exp.is_negative is True]))
        invariant = clean(matrix*lift-lift*restricted)
        currents, grade_residuals = {}, {}
        for name in ('SLAB_CURRENT_MATRIX', 'BULK_NORMAL_CURRENT_DENSITY_MATRIX',
                     'BULK_DEPTH_CURRENT_MATRIX', 'INTERFACE_POWER_MATRIX'):
            value, residual = grade(result[name])
            currents[name] = bind(value, label, mapping)
            grade_residuals[name] = residual
        amplitude_legs = tuple(tuple(sorted((symbol for symbol in dimensions
            if getattr(symbol,'name','').startswith('s11cdCurrent'+side+'Amplitude')),
            key=lambda symbol:symbol.name)) for side in ('Plus','Minus'))
        require(all(len(amplitudes)==5 for amplitudes in amplitude_legs),
                'five physical amplitude columns in each harmonic leg')
        rows, reconstruction, raw_face_expressions = [], [], []
        flat_exterior = []
        for face_index_value,face in enumerate(result['FACE_LEG_OBJECTS']):
            require(len(face)==2,'both saved harmonic face legs required')
            face_rows,face_residuals,face_expressions = [],[],[]
            for leg,amplitudes in enumerate(amplitude_legs):
                row,residual,expressions = {},{},{}
                for name in DRIVES:
                    expression = face[leg][name]
                    coefficients = sp.ImmutableMatrix(1,5,
                        lambda i,j:sp.diff(expression,amplitudes[j]))
                    row[name] = bind(coefficients,label,mapping)
                    residual[name] = expression-(coefficients*sp.ImmutableMatrix(amplitudes))[0]
                    expressions[name] = expression
                face_rows.append(row);face_residuals.append(residual);face_expressions.append(expressions)
                pressure_scale = bind(result['OPEN_BULK_PRESSURE_SCALES'][leg],label,mapping)
                velocity_scale = bind(result['OPEN_BULK_VELOCITY_SCALES'][leg],label,mapping)
                flat_exterior.append(dict(face=face_index_value,harmonicLeg=leg,
                    pressure=row['PRESSURE'],amplitude=row['AMPLITUDE'],bulkVelocity=row['BULK_VELOCITY'],
                    pressureScale=pressure_scale,velocityScale=velocity_scale,
                    pressureResidual=clean(row['PRESSURE']-pressure_scale*row['AMPLITUDE']),
                    velocityResidual=clean(row['BULK_VELOCITY']-velocity_scale*row['AMPLITUDE'])))
            rows.append(tuple(face_rows));reconstruction.append(tuple(face_residuals))
            raw_face_expressions.append(tuple(face_expressions))
        depth = result['OPEN_BULK_CURRENT_COEFFICIENTS'][3]
        e = next(s for s in depth.free_symbols if s.name == 'epsilon_shape')
        free_amps = {s:sp.S.One for s in depth.free_symbols
                     if s.name in ('s11cdAcousticLeftAmplitude', 's11cdAcousticRightAmplitude')}
        depth = bind(sp.Poly(depth,e).nth(2), label, {**mapping, **free_amps})
        carrier = carrier_law_check(label,result,mapping,depth)
        denominator_operands = [*raw_matrix,*raw_restricted,*raw_weak,*raw_gram,
            *[value for current in currents.values() for value in current],
            *[value for face in rows for row in face for coefficients in row.values() for value in coefficients]]
        raw_source_values = [*result['CLOSED_PENCIL_LEGS'][0],
            *[value for name in currents for value in result[name]],
            *[expression for face in raw_face_expressions for leg in face for expression in leg.values()]]
        carried_symbols = {symbol:symbol for columns in amplitude_legs for symbol in columns}
        carried_symbols.update({symbol:symbol for value in raw_source_values for symbol in value.free_symbols
                                if symbol.name=='epsilon_shape'})
        saved_source_factors = raw_denominators(raw_source_values)
        bound_saved_source_factors = tuple(bind(factor,label,{**mapping,**carried_symbols})
                                          for factor in saved_source_factors)
        denominators = tuple(dict.fromkeys((*raw_denominators(denominator_operands),*bound_saved_source_factors)))
        coupling = [uniform['records'][label]['coupling'][name] for name in ('TH','HT')]
        checks = dict(sourceJoin=zero(source_join), waveJoin=zero(wave_join),
            suppliedCurlJoin=zero(common_curl_join), carrierLawAndNormalization=carrier['joins'],
            originalBranch=zero(packet['branchResiduals']), pairingBranch=zero(result['SOURCE_BRANCH_JOINS']),
            savedCoupling=zero(coupling), invariantT=zero(invariant),
            weakHermitian=zero(weak-weak.H), grade=zero(grade_residuals),
            driveLinear=zero(reconstruction),
            flatExteriorJoin=zero([x[key] for x in flat_exterior for key in ('pressureResidual','velocityResidual')]),
            sourceFrequencyLive=any(s.name == 'omega' for s in packet['originalAlgebraic'].free_symbols),
            sourceTangentsLive=all(any(s.name == name for s in packet['originalAlgebraic'].free_symbols)
                for name in ('s11cdTangentialMomentum1','s11cdTangentialMomentum2')),
            actualMatrixFrequencyLive=matrix.has(w), actualMatrixTangentialLive=matrix.has(v),
            restrictedFrequencyLive=restricted.has(w), restrictedTangentialLive=restricted.has(v),
            declaredSymbols=not matrix.free_symbols-{w,v,k,q})
        return dict(end=label, checks=checks, originalSource=packet['originalAlgebraic'],
            originalRelation=packet['originalRelation'], sourceJoin=source_join, waveJoin=wave_join,
            savedSourceBranch=packet['branchResiduals'], pairingSourceBranch=result['SOURCE_BRANCH_JOINS'],
            sourceCoupling=coupling, matrix=matrix, rawMatrix=raw_matrix, wave=wave, q2=q2, scale=scale,
            lift=lift, gram=gram, suppliedCurlJoin=common_curl_join, savedGauge=gauge,
            longitudinalDirectionProbes=longitudinal_probes, rawWeak=raw_weak, rawRestricted=raw_restricted,
            rawTDenominators=raw_t_denominators, chart=sp.factor(gram.det()), restricted=restricted, weak=weak,
            invariantResidual=invariant, weakHermitianResidual=weak-weak.H,
            currents=currents, gradeResiduals=grade_residuals, drives=rows,
            harmonicAmplitudeColumns=amplitude_legs,rawSavedFaceExpressions=raw_face_expressions,
            sourceDenominatorOperands=denominator_operands,
            rawSavedSourceDenominatorFactors=saved_source_factors,
            boundSavedSourceDenominatorFactors=bound_saved_source_factors,
            driveReconstruction=reconstruction, flatExteriorJoins=flat_exterior, sourceDenominators=denominators,
            outgoingDepthCoefficient=depth, carrierLawCheck=carrier,
            flatExteriorJoinScope='SOURCE_INTERNAL_CONSISTENCY_NOT_INDEPENDENT_CLOSURE',dimensions=dimensions)

    def faces(bound):
        plus_lift = bound['lift']; records = []
        literals = []
        if bound['end'] == 'REFERENCE':
            for item in face_index['records']:
                if item['tag'] != 'PY_S11CC2_FOLD_SYMBOL_MAP_LAB_HELD_RHO4_CONSTANT': continue
                for identity in item['velocityIdentifications']:
                    text = identity['savedValueSrepr']
                    if hashlib.sha256(text.encode()).hexdigest() != identity['literalSha256']:
                        raise IntegrityError('c2 index literal changed')
                    expression = sp.sympify(text)
                    field = next(symbol for symbol in expression.free_symbols if symbol.name=='e_W_t')
                    literals.append(dict(indexFaceLabel=item['face'],source=expression,
                        coefficientsByHarmonicLeg=tuple(bind(expression,'REFERENCE',{field:rate})
                            for rate in (-sp.I*w,sp.I*w)),literalSha256=identity['literalSha256']))
        for index,face_rows in enumerate(bound['drives']):
            for leg,row in enumerate(face_rows):
                lift = plus_lift if leg==0 else sp.conjugate(plus_lift)
                contractions = {name:clean(row[name]*lift) for name in DRIVES}
                loaded_controls = []
                for name in DRIVES:
                    urow = row[name][:,:3]
                    if zero(urow):
                        loaded_controls.append(dict(drive=name,status='NOT_APPLICABLE_NO_U_DEPENDENCE',
                            originalRow=row[name],uRow=urow,baseline=contractions[name]))
                        continue
                    terms = [(slot,term) for slot in range(3)
                             for term in sp.Add.make_args(sp.expand(row[name][0,slot])) if term!=0]
                    columns = []
                    for column in range(lift.cols):
                        contributions = [(slot,term,sp.cancel(term*lift[slot,column])) for slot,term in terms]
                        rational = all(value.is_rational_function(w,v,k,q,qb) is True
                                       for _,_,value in contributions)
                        # Exact rational identity, not an assumptions-based
                        # nonzero-at-every-point claim on the k/v domain.
                        loaded = [item for item in contributions if item[2]!=0] if rational else []
                        if not loaded:
                            structural = rational and all(value==0 for _,_,value in contributions)
                            columns.append(dict(column=column,contributions=contributions,
                                status='NOT_APPLICABLE_STRUCTURALLY_UNLOADED_COLUMN' if structural else 'UNRESOLVED'))
                            continue
                        slot,term,contribution = loaded[0]
                        mutated = sp.MutableDenseMatrix(row[name]); mutated[0,slot]-=term
                        before = row[name]*lift[:,column]
                        after = mutated*lift[:,column]; raw_movement = after-before
                        movement = clean(raw_movement)
                        columns.append(dict(column=column,slot=slot,omittedNativeUTerm=term,
                            sourceContributions=contributions,originalRow=row[name],mutatedRow=mutated,
                            fixedLoadedLift=lift[:,column],before=before,after=after,
                            rawMovement=raw_movement,reducedMovement=movement,
                            status='GENERIC_RATIONAL_NONIDENTITY' if not zero(movement) else 'UNRESOLVED',
                            pointwiseNonzeroClaim=False))
                    supported = all(item['status']!='UNRESOLVED' for item in columns)
                    loaded_controls.append(dict(drive=name,status='SUPPORTED_GENERIC_DEPENDENCE' if supported else 'UNRESOLVED',
                        originalRow=row[name],baseline=contractions[name],columns=columns))
                slot_diagnostics = []
                for name in DRIVES:
                    for slot,field_name in ((3,'Theta'),(4,'E')):
                        original = row[name]; probe = sp.eye(5)[:,slot]
                        changed = sp.MutableDenseMatrix(original);changed[0,slot]=0
                        before,after = original*probe,changed*probe
                        raw_difference = after-before; reduced = clean(raw_difference)
                        coefficient = sp.cancel(original[0,slot])
                        rational = coefficient.is_rational_function(w,v,k,q,qb) is True
                        status = ('NOT_APPLICABLE_STRUCTURAL_ZERO' if coefficient==0 else
                                  'GENERIC_COEFFICIENT_PRESENT' if rational else 'UNRESOLVED')
                        slot_diagnostics.append(dict(drive=name,field=field_name,slot=slot,status=status,
                            scope='NATIVE_SLOT_PRESENCE_ONLY_NOT_T_CANCELLATION_OR_POWER',
                            originalRow=original,changedRow=changed,probe=probe,before=before,after=after,
                            rawDifference=raw_difference,reducedDifference=reduced,coefficient=coefficient))
                ew = next(item for item in slot_diagnostics
                          if item['drive']=='OUTWARD_VELOCITY' and item['field']=='E')
                scaled_velocity_coefficient = sp.cancel(ew['coefficient']/w)
                ew_calibrated = (not scaled_velocity_coefficient.free_symbols
                                 and scaled_velocity_coefficient!=0 and w.is_positive is True)
                ew = dict(ew,frequencyNormalizedCoefficient=scaled_velocity_coefficient,
                          frequencyDomain=sp.Gt(w,0,evaluate=False),presenceCalibrationSupported=ew_calibrated)
                probe_vectors = [probe if leg==0 else sp.conjugate(probe)
                                 for probe in bound['longitudinalDirectionProbes']]
                probes = [dict(kind='SAVED_CURL_GAUGE_DIRECTION_EMBEDDED_AS_U_PROBE',vector=probe,
                    normalization='Arbitrary coordinate probe; not a normalized physical longitudinal mode',
                    values={name:row[name]*probe for name in DRIVES}) for probe in probe_vectors]
                # Index values are orientation-blind. The two harmonic signs
                # come from their actual time characters, not face orientation.
                original_velocity = row['OUTWARD_VELOCITY']
                joins = [sp.cancel(original_velocity[0,4]-item['coefficientsByHarmonicLeg'][leg]) for item in literals]
                records.append(dict(faceOrdinal=index,nativeSourceOrientation=face_order[index],harmonicLeg=leg,
                    harmonicMeaning='POSITIVE_CHARACTER' if leg==0 else 'NEGATIVE_CONJUGATE_CHARACTER',
                    amplitudeColumns=bound['harmonicAmplitudeColumns'][leg],nativeRows=row,lift=lift,
                    contractions=contractions,driveZero=zero(contractions),loadedUTermControls=loaded_controls,
                    nontransverseSourceProbes=probes,nativeThetaESlotDiagnostics=slot_diagnostics,
                    eWPresenceCalibration=ew,c2IndexLiteralJoins=joins,
                    c2IndexLiteralSupported=(len(literals)==2) if bound['end']=='REFERENCE' else None,
                    c2OrientationIndependentlyChecked=False))
        forms = {name:clean(plus_lift.H*matrix*plus_lift) for name,matrix in bound['currents'].items()
                 if name!='SLAB_CURRENT_MATRIX'}
        checks = dict(allHarmonicDrivesZero=all(item['driveZero'] for item in records),
            nativeEVelocityPresence=all(item['eWPresenceCalibration']['presenceCalibrationSupported'] for item in records),
            loadedUControls=all(control['status']!='UNRESOLVED' for item in records for control in item['loadedUTermControls']),
            c2IndexLiteral=(bound['end']!='REFERENCE' or all(item['c2IndexLiteralSupported'] and
                zero(item['c2IndexLiteralJoins']) for item in records)),lossSideFormsZero=zero(forms))
        return dict(end=bound['end'],faces=records,physicalForms=forms,checks=checks,
            harmonicLegScope='BOTH_SAVED_REAL_FIELD_LEGS_NOT_TWO_ACOUSTIC_INCIDENCE_WAVES',
            c2LiteralEvidence=literals,c2EvidenceScope='ORIENTATION_BLIND_INDEX_LITERAL_ONLY',
            flatExteriorSourceInternalJoins=bound['flatExteriorJoins'],
            c2PressureTraceFirstShapeJoin='NOT_SUPPORTED_BY_THIS_INPUT_INDEX')

    def determinant_data(matrix):
        # Preserve row-clearing operands before cancellation; no H/full pencil
        # determinant, radical resultant, or discriminant is constructed.
        raw_denominators = tuple(sp.fraction(sp.together(x))[1] for x in matrix)
        row_denominators = tuple(sp.prod(raw_denominators[i*matrix.cols:(i+1)*matrix.cols])
                                 for i in range(matrix.rows))
        cleared = sp.ImmutableMatrix(matrix.rows,matrix.cols,
            lambda i,j:sp.cancel(row_denominators[i]*matrix[i,j]))
        cleared_determinant = sp.expand(cleared.det())
        determinant = sp.cancel(cleared_determinant/sp.prod(row_denominators))
        numerator, denominator = sp.fraction(determinant)
        return dict(matrix=matrix, determinant=determinant, numerator=numerator, denominator=denominator,
                    rawEntryDenominators=raw_denominators,rowDenominators=row_denominators,
                    rowClearedMatrix=cleared,rowClearedDeterminant=cleared_determinant)

    def real_exclusion_certificate(expression, w0, v0):
        specialized = sp.together(expression.subs({w:w0,v:v0}))
        numerator,denominator = sp.fraction(specialized)
        parts = []
        for part in (numerator,denominator):
            real,imag = map(sp.expand,part.as_real_imag())
            rp,ip = sp.Poly(real,k,domain=sp.QQ),sp.Poly(imag,k,domain=sp.QQ)
            if rp.is_zero and ip.is_zero:
                parts.append(dict(real=real,imag=imag,identicallyZero=True,realZeroCount=None))
                continue
            gcd = sp.gcd(rp,ip)
            parts.append(dict(real=real,imag=imag,gcd=gcd.as_expr(),identicallyZero=False,
                realZeroCount=gcd.count_roots(-sp.oo,sp.oo) if gcd.degree()>0 else 0))
        return dict(expression=expression,specialized=specialized,parts=parts,
                    noRealExcludedPoint=all(x['realZeroCount']==0 for x in parts))

    def census(bound, w0, v0):
        specialized = bound['rawRestricted'].subs({w:w0,v:v0})
        data = determinant_data(specialized)
        if specialized.has(q, qb): return dict(status='UNRESOLVED', reason='T_RESTRICTION_RETAINS_DEPTH_ROOT', **data)
        # These are exclusively transverse row/coordinate exclusions. Their
        # global real-k check licenses pointwise T absence/completeness only.
        exclusions = [real_exclusion_certificate(x,w0,v0) for x in
            (*bound['rawTDenominators'],bound['chart'],*data['rawEntryDenominators'],
             *data['rowDenominators'],data['denominator'])]
        real, imag = map(sp.expand, data['rowClearedDeterminant'].as_real_imag())
        rp, ip = sp.Poly(real,k,domain=sp.QQ), sp.Poly(imag,k,domain=sp.QQ)
        if rp.is_zero and ip.is_zero:
            return dict(status='UNRESOLVED', reason='DEGENERATE_T_DETERMINANT', exclusions=exclusions, **data)
        gcd = sp.gcd(rp, ip)
        # Keep multiplicities. A twofold transverse root can contain two valid
        # polarizations; actual lifted kernel/current ranks decide support.
        intervals = gcd.intervals(eps=sp.Rational(1,10**28)) if gcd.degree()>0 else []
        distinct = gcd.sqf_part().count_roots(-sp.oo,sp.oo) if gcd.degree()>0 else 0
        multiplicity_count = sum(mult for _,mult in intervals)
        factorization = sp.factor_list(gcd.as_expr(),k)
        return dict(status='REAL_T_CANDIDATES', omega=w0, ray=v0, kappa=sp.sqrt(5)*v0,
            realPolynomial=rp.as_expr(), imaginaryPolynomial=ip.as_expr(), gcd=gcd.as_expr(),
            realRemainder=sp.rem(rp,gcd).as_expr(), imaginaryRemainder=sp.rem(ip,gcd).as_expr(),
            intervals=intervals, count=distinct, totalRootMultiplicity=multiplicity_count,
            factorization=factorization, exclusions=exclusions,
            noRealExcludedPoint=all(x['noRealExcludedPoint'] for x in exclusions),
            coverage=(distinct==len(intervals)), **data)

    def no_pole_interval(expression, bound, point, interval, depth_sign):
        lo,hi = interval
        mapping = {w:point['omega'],v:point['ray'],qb:depth_sign*q}
        specialized = expression.subs(mapping)
        q2_raw = bound['q2'].subs({w:point['omega'],v:point['ray']})
        q2 = sp.cancel(q2_raw)
        stages,zero_tests,coefficient_factors = [],[],[]
        pending = [('source',specialized)]
        queued = {specialized}
        seen_factors = set()

        def retain_factor(value, stage):
            if value.has(q):
                if value not in queued:
                    queued.add(value);pending.append((stage,value))
            elif value not in seen_factors:
                seen_factors.add(value);coefficient_factors.append((stage,value))

        def split(value, stage):
            # Preserve the source expression and its uncancelled negative-power
            # factors before together/cancel can merge removable factors.
            literal = raw_denominators([value])
            together = sp.together(value)
            numerator,denominator = sp.fraction(together)
            for factor in (*literal,*raw_denominators([together])):
                retain_factor(factor,stage+'-domain-factor')
            stages.append(dict(stage=stage,raw=value,together=together,
                numerator=numerator,denominator=denominator,uncancelledFactors=literal))
            return numerator,denominator

        def reduce_polynomial(value, stage):
            polynomial = sp.Poly(value,q,domain='EX')
            modulus = sp.Poly(q*q-q2,q,domain='EX')
            before_factors = raw_denominators([*polynomial.all_coeffs(),*modulus.all_coeffs()])
            for factor in before_factors:
                retain_factor(factor,stage+'-input-coefficient-denominator')
            reduced_poly = sp.rem(polynomial,modulus)
            after_factors = raw_denominators(reduced_poly.all_coeffs())
            for factor in after_factors:
                retain_factor(factor,stage+'-introduced-coefficient-denominator')
            reduced = reduced_poly.as_expr()
            stages.append(dict(stage=stage,polynomial=polynomial.as_expr(),modulus=modulus.as_expr(),
                reduced=reduced,inputCoefficientDenominators=before_factors,
                outputCoefficientDenominators=after_factors))
            return reduced

        try:
            for factor in raw_denominators([q2_raw]):
                retain_factor(factor,'physical-quadratic-coefficient-denominator')
            cursor = 0
            while cursor<len(pending):
                name,value = pending[cursor];cursor+=1
                numerator,denominator = split(value,name+'-fraction')
                for side,part in (('numerator',numerator),('denominator',denominator)):
                    tag = name+'-'+side
                    reduced = reduce_polynomial(part,tag+'-quadratic-reduction')
                    # This product is a rational expression. Both its numerator
                    # AND denominator are reduced and retained independently.
                    product = sp.Mul(reduced,reduced.xreplace({q:-q}),evaluate=False)
                    product_n,product_d = split(product,tag+'-opposite-sheet-product')
                    norm_n = reduce_polynomial(product_n,tag+'-norm-numerator')
                    norm_d = reduce_polynomial(product_d,tag+'-norm-denominator')
                    stages.append(dict(stage=tag+'-norm-rational-pair',
                        numerator=norm_n,denominator=norm_d,uncancelledProduct=product))
                    if norm_n.has(q) or norm_d.has(q):
                        raise ValueError('unsupported residual depth dependence in denominator norm')
                    retain_factor(norm_n,tag+'-norm-numerator-nonzero')
                    retain_factor(norm_d,tag+'-norm-denominator-nonzero')
            # Every retained scalar factor is tested on this root interval only.
            # Complex rational coefficients use a real/imaginary gcd; neither
            # an identically zero numerator nor an unsupported domain passes.
            for stage,value in coefficient_factors:
                numerator,denominator = sp.fraction(sp.together(value))
                tests = []
                for side,part in (('numerator',numerator),('denominator',denominator)):
                    real,imag = map(sp.expand,part.as_real_imag())
                    rp,ip = sp.Poly(real,k,domain=sp.QQ),sp.Poly(imag,k,domain=sp.QQ)
                    if rp.is_zero and ip.is_zero:
                        tests.append(dict(side=side,real=real,imaginary=imag,identicallyZero=True,intervalRootCount=None))
                    else:
                        gcd = sp.gcd(rp,ip)
                        tests.append(dict(side=side,real=real,imaginary=imag,gcd=gcd.as_expr(),
                            identicallyZero=False,intervalRootCount=gcd.count_roots(lo,hi) if gcd.degree()>0 else 0))
                zero_tests.append(dict(stage=stage,factor=value,numerator=numerator,denominator=denominator,tests=tests))
        except (sp.PolynomialError,sp.polys.polyerrors.CoercionFailed,ValueError) as error:
            return dict(status='UNRESOLVED',expression=expression,specialized=specialized,
                quadratic=q2_raw,stages=stages,retainedCoefficientFactors=coefficient_factors,
                tests=zero_tests,unsupportedType=type(error).__name__,unsupportedReason=str(error),excluded=True)
        supported = all(test['intervalRootCount']==0 for entry in zero_tests for test in entry['tests'])
        return dict(status='CERTIFIED_ON_INTERVAL' if supported else 'UNRESOLVED',
            expression=expression,specialized=specialized,quadratic=q2_raw,stages=stages,
            retainedCoefficientFactors=coefficient_factors,tests=zero_tests,excluded=not supported)

    def numeric(matrix, mapping):
        return np.asarray(matrix.subs(mapping).evalf(40).tolist(), dtype=complex)

    def numerical_face_legs(bound,mapping,basis):
        records = []
        for face_index_value,face_rows in enumerate(bound['drives']):
            for leg,row in enumerate(face_rows):
                loaded_basis = basis if leg==0 else basis.conj()
                native_rows = {name:numeric(row[name],mapping) for name in DRIVES}
                contractions = {name:value@loaded_basis for name,value in native_rows.items()}
                e_coefficient = native_rows['OUTWARD_VELOCITY'][0,4]
                records.append(dict(faceOrdinal=face_index_value,harmonicLeg=leg,
                    nativeRows=native_rows,loadedBasis=loaded_basis,contractions=contractions,
                    nativeEVelocityCoefficient=e_coefficient,
                    nativeEVelocityPresent=abs(e_coefficient)>1e-12))
        return records

    def loaded_control(bound, symbol, mapping, basis, scale):
        matrix = bound['matrix']; records = []
        baseline_matrix = numeric(matrix,mapping)
        original = (baseline_matrix/scale[:,None])@basis
        for column in range(basis.shape[1]):
            selected,best = None,0.
            for i in range(matrix.rows):
                for j in range(matrix.cols):
                    if abs(basis[j,column])<=1e-12 or not matrix[i,j].has(symbol): continue
                    for term in sp.Add.make_args(sp.expand(matrix[i,j])):
                        if not term.has(symbol): continue
                        movement = abs(complex(term.subs(mapping).evalf(40))*basis[j,column])/scale[i]
                        if movement>best:
                            best,selected = movement,(i,j,term)
            if selected is None:
                records.append(dict(column=column,status='UNRESOLVED',reason='NO_LOADED_DEPENDENT_SOURCE_TERM'))
                continue
            i,j,term = selected; changed = sp.MutableDenseMatrix(matrix); changed[i,j]-=term
            altered_matrix = numeric(changed,mapping)
            altered = (altered_matrix/scale[:,None])@basis
            difference = altered-original
            column_norm = float(np.linalg.norm(difference[:,column]))
            records.append(dict(column=column,status='RESPONSIVE' if column_norm>1e-8 else 'UNRESOLVED',
                selectedEntry=(i,j),omittedNativeTerm=term,mutatedMatrix=altered_matrix,
                mutatedResidual=altered,difference=difference,loadedColumnMovement=column_norm))
        return dict(status='RESPONSIVE' if records and all(x['status']=='RESPONSIVE' for x in records) else 'UNRESOLVED',
            scope='PER_LOADED_POLARIZATION_SOURCE_DEPENDENCE_NOT_INDEPENDENT_PHYSICS',symbol=symbol,
            baselineMatrix=baseline_matrix,fixedLoadedBasis=basis,baselineResidual=original,columns=records)

    def depth_selection(bound, q2, w0, v0, k0):
        candidates = (sp.sqrt(q2),-sp.sqrt(q2))
        records = []
        for candidate in candidates:
            mapping = {w:w0,v:v0,k:k0,q:candidate,qb:sp.conjugate(candidate)}
            native_flux = sp.cancel(bound['outgoingDepthCoefficient'].subs(mapping))
            carrier_flux = sp.cancel(bound['carrierLawCheck']['normalizedFlux'].subs(mapping))
            residual = sp.cancel(native_flux-carrier_flux)
            numeric_flux = complex(native_flux.evalf(40))
            if q2.is_positive is True:
                selected = residual==0 and abs(numeric_flux.imag)<1e-10 and numeric_flux.real>0
                criterion = 'REAL_DEPTH_POSITIVE_SOURCE_AND_CARRIER_OUTWARD_FLUX'
            else:
                selected = residual==0 and sp.im(candidate).is_positive is True
                criterion = 'IMAGINARY_DEPTH_DECAY_NOT_PROPAGATING_FLUX'
            records.append(dict(candidate=candidate,nativeFlux=native_flux,carrierFlux=carrier_flux,
                normalizationResidual=residual,selected=bool(selected),criterion=criterion))
        selected = [item['candidate'] for item in records if item['selected']]
        return dict(status='SUPPORTED' if len(selected)==1 else 'UNRESOLVED',
                    candidates=records,selected=selected[0] if len(selected)==1 else None)

    def inspect_root(bound, point, interval):
        (lo,hi), multiplicity = interval; km = (lo+hi)/2
        base = dict(status='UNRESOLVED', interval=interval, multiplicity=multiplicity)
        if lo<=0<=hi: return dict(base, reason='ZERO_NORMAL_CURRENT_THRESHOLD')
        q2point = sp.cancel(bound['q2'].subs({w:point['omega'],v:point['ray']}))
        qpolys = [sp.Poly(x,k,domain=sp.QQ) for x in sp.fraction(q2point)]
        if any(p.is_zero or (p.degree()>0 and p.count_roots(lo,hi)>0) for p in qpolys):
            return dict(base, reason='BULK_BRANCH_OR_POLE_INTERVAL', q2=q2point)
        q2 = sp.cancel(q2point.subs(k,km))
        if q2.is_positive is not True and q2.is_negative is not True:
            return dict(base, reason='DEPTH_SIGN_UNRESOLVED', q2=q2)
        sign = 1 if q2.is_positive else -1
        try:
            exclusions = [no_pole_interval(x,bound,point,(lo,hi),sign)
                          for x in (*bound['sourceDenominators'],bound['chart'])]
        except (sp.PolynomialError,sp.polys.polyerrors.CoercionFailed):
            return dict(base, reason='DENOMINATOR_CERTIFICATE_DOMAIN_UNSUPPORTED')
        if any(x['excluded'] for x in exclusions):
            return dict(base, reason='SOURCE_POLE_OR_CHART_NOT_EXCLUDED', exclusions=exclusions)
        selection = depth_selection(bound,q2,point['omega'],point['ray'],km)
        if selection['status']!='SUPPORTED':
            return dict(base,reason='DEPTH_SELECTION_UNRESOLVED',depthSelection=selection,exclusions=exclusions)
        q0 = selection['selected']
        mapping = {w:point['omega'],v:point['ray'],k:km,q:q0,qb:sp.conjugate(q0)}
        matrix = numeric(bound['matrix'],mapping); lift = numeric(bound['lift'],mapping)
        restricted = numeric(bound['restricted'],mapping)
        row_scale = np.maximum(np.linalg.norm(matrix,axis=1),1.)
        _, singular, vh = np.linalg.svd(restricted)
        tolerance = 1e-10*max(1.,float(np.linalg.norm(restricted)))
        null = singular<tolerance
        base.update(physicalDepth=q0, depthSquared=q2, exclusions=exclusions,depthSelection=selection,
                    matrix=matrix, lift=lift, restricted=restricted, singularValues=singular)
        if not null.any() or np.any((singular>=tolerance)&(singular<100*tolerance)):
            return dict(base, reason='T_RANK_UNRESOLVED')
        physical = lift@vh.conj().T[:,null]
        basis, _ = np.linalg.qr(physical)
        residual = (matrix/row_scale[:,None])@basis
        if np.linalg.norm(residual)>1e-8:
            return dict(base, reason='PHYSICAL_FULL_PENCIL_RESIDUAL', basis=basis, residual=residual)
        forms = {name:basis.conj().T@numeric(value,mapping)@basis
                 for name,value in bound['currents'].items()}
        face = numerical_face_legs(bound,mapping,basis)
        current = forms['SLAB_CURRENT_MATRIX']; hermitian = current-current.conj().T
        eig, rotation = np.linalg.eigh((current+current.conj().T)/2)
        depth_coefficient = complex(bound['outgoingDepthCoefficient'].subs(mapping).evalf(40))
        controls = [loaded_control(bound,symbol,mapping,basis,row_scale) for symbol in (w,v)]
        base.update(basis=basis,geometricNullity=basis.shape[1],
                    rootMultiplicityAccounted=(basis.shape[1]==multiplicity),
                    rankInterpretation='NUMERICAL_SUPPORTED_CURRENT_RANK_LOWER_BOUND',
                    fullResidual=residual, rowScale=row_scale, forms=forms,
                    face=face, current=current, hermitianResidual=hermitian, currentEigenvalues=eig,
                    rotation=rotation, controls=controls, outgoingDepthCoefficient=depth_coefficient)
        if q2.is_positive and (abs(depth_coefficient.imag)>1e-10 or depth_coefficient.real<=0):
            return dict(base, reason='OUTGOING_DEPTH_SIGN_UNRESOLVED')
        if any(not item['nativeEVelocityPresent'] for item in face) or any(
                np.linalg.norm(value)>1e-8 for item in face for value in item['contractions'].values()) or any(
                np.linalg.norm(value)>1e-8 for name,value in forms.items() if name!='SLAB_CURRENT_MATRIX'):
            return dict(base, reason='LOSS_SIDE_T_DRIVE_OR_FORM_NOT_ZERO')
        if np.linalg.norm(hermitian)>1e-8*max(1.,np.linalg.norm(current)) or np.min(np.abs(eig))<1e-9:
            return dict(base, reason='CURRENT_REALITY_OR_RANK_UNRESOLVED')
        if any(control['status']!='RESPONSIVE' for control in controls):
            return dict(base, reason='LOADED_SOURCE_CONTROL_UNRESOLVED')
        return dict(base, status='AVAILABLE', positiveCurrentRank=int(sum(eig>0)),
                    negativeCurrentRank=int(sum(eig<0)), currentDomain='EXACTLY_UNDRIVEN_T_EXTERIOR')

    def saved_seed(bound, old):
        info = old['info']
        k0, q0 = complex(info['K']), complex(info['PHYSICAL_Q'])
        w0 = sp.Rational(str(info['OMEGA'].real)) if isinstance(info['OMEGA'],complex) else sp.Rational(str(info['OMEGA']))
        v0 = params['s11cdTangentialMomentum2']
        mapping = {w:w0,v:v0,k:k0,q:q0,qb:q0.conjugate()}
        matrix = numeric(bound['matrix'],mapping)
        lift = numeric(bound['lift'],mapping)
        # Projection through the 2-column geometric lift, not a new modal solve.
        gram = lift.conj().T@lift
        dual = np.linalg.inv(gram)@lift.conj().T
        right = np.asarray(old['right'],complex)
        t_residual = right-lift@dual@right
        scale = np.maximum(np.linalg.norm(matrix,axis=1),1.)
        pencil_join = (matrix-np.asarray(old['pencil'],complex))/scale[:,None]
        kernel = (matrix/scale[:,None])@right
        face = numerical_face_legs(bound,mapping,right)
        slab = numeric(bound['currents']['SLAB_CURRENT_MATRIX'],mapping)
        bulk = numeric(bound['currents']['BULK_NORMAL_CURRENT_DENSITY_MATRIX'],mapping)
        current = right.conj().T@slab@right
        projected_bulk = right.conj().T@bulk@right
        operand_joins = {'slab':slab-np.asarray(old['currentOperands']['CURRENT_SLAB'],complex),
                         'bulk':bulk-np.asarray(old['currentOperands']['CURRENT_BULK'],complex)}
        controls = [loaded_control(bound,symbol,mapping,right,scale) for symbol in (w,v)]
        wave = complex(bound['wave'].subs(mapping).evalf(40))
        normal_reality_residual = k0-k0.conjugate()
        normal_is_real = k0.imag == 0.0
        face_zero = all(np.linalg.norm(value)<1e-8 for item in face
                        for value in item['contractions'].values())
        bulk_zero = np.linalg.norm(projected_bulk)<1e-8
        hermitian_residual = current-current.conj().T
        hermitian_part = (current+current.conj().T)/2
        current_eigenvalues = np.linalg.eigvalsh(hermitian_part)
        hermitian_supported = np.linalg.norm(hermitian_residual)<1e-8*max(1.,np.linalg.norm(current))
        current_rank_supported = bool(current_eigenvalues.size and np.min(np.abs(current_eigenvalues))>1e-9)
        transport = dict(scope='RECOMPUTED_SLAB_FORM_ONLY_ON_UNDRIVEN_T_ZERO_PROJECTED_BULK_DOMAIN',
            domainSupported=bool(face_zero and bulk_zero),right=right,slabMatrix=slab,bulkMatrix=bulk,
            currentForm=current,projectedBulk=projected_bulk,hermitianResidual=hermitian_residual,
            hermitianPart=hermitian_part,eigenvaluesOfHermitianPart=current_eigenvalues,
            hermitianTolerance=1e-8*max(1.,float(np.linalg.norm(current))),rankTolerance=1e-9,
            positiveCurrentRank=int(sum(current_eigenvalues>1e-9)),
            negativeCurrentRank=int(sum(current_eigenvalues < -1e-9)),
            interpretation='Tolerance-based current check on the saved basis; no new normalized basis')
        history = dict(currentDefined=old.get('currentDefined'),
            sourceSheetMembership=info.get('SHEET_MEMBERSHIP'),sourceExactRealNormal=info.get('EXACT_REAL_NORMAL'),
            sourceBulkDecayDiskCertified=info.get('BULK_DECAY_DISK_CERTIFIED'),
            sourcePhysicalCurrentNormalization=info.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'),
            scope='HISTORICAL_METADATA_NOT_RECOMPUTED_AVAILABILITY_GATE')
        missing = [name for name in ('depthIntegral','currentGram') if name not in old]
        if missing:
            historical_comparison = dict(status='UNAVAILABLE',missingFields=missing,
                availableSavedDepthIntegral=old.get('depthIntegral'),availableSavedGram=old.get('currentGram'),
                availabilityGated=False)
        else:
            depth = complex(old['depthIntegral'])
            weighted_matrix = slab+depth*bulk
            weighted_gram = right.conj().T@weighted_matrix@right
            saved_gram = np.asarray(old['currentGram'],complex)
            weighted_residual = weighted_gram-saved_gram
            comparison_supported = np.linalg.norm(weighted_residual)<1e-8
            historical_comparison = dict(status='COMPARED',right=right,slabMatrix=slab,bulkMatrix=bulk,
                savedDepthIntegral=depth,weightedMatrix=weighted_matrix,recomputedWeightedGram=weighted_gram,
                savedGram=saved_gram,residual=weighted_residual,residualNorm=float(np.linalg.norm(weighted_residual)),
                tolerance=1e-8,supported=bool(comparison_supported),availabilityGated=True,
                source='S11c_d_uniform_response.py:74-88; no depth integral is re-evaluated')
        depth_evidence = []
        for candidate in (q0,-q0):
            dm = dict(mapping); dm.update({q:candidate,qb:candidate.conjugate()})
            native_flux = complex(bound['outgoingDepthCoefficient'].subs(dm).evalf(40))
            carrier_flux = complex(bound['carrierLawCheck']['normalizedFlux'].subs(dm).evalf(40))
            normalization_residual = native_flux-carrier_flux
            decay = bool(abs(candidate.real)<1e-12 and candidate.imag>0)
            real_outward = bool(abs(candidate.imag)<1e-12 and abs(native_flux.imag)<1e-10 and native_flux.real>0)
            depth_evidence.append(dict(candidate=candidate,nativeFlux=native_flux,carrierFlux=carrier_flux,
                normalizationResidual=normalization_residual,decayingImaginaryDepth=decay,
                realOutwardDepth=real_outward,coordinateRealityTolerance=1e-12,
                historicalDecayDisk=info.get('BULK_DECAY_DISK_CERTIFIED'),
                supported=(abs(normalization_residual)<1e-10 and (decay or real_outward))))
        checks = dict(savedFrequency=(w0==params['omega']),
            transverseMembership=np.linalg.norm(t_residual)<1e-8,
            physicalPencilJoin=np.linalg.norm(pencil_join)<1e-8,
            physicalKernel=np.linalg.norm(kernel)<1e-8, physicalWave=abs(wave)<1e-8,
            recomputedRealNormal=normal_is_real,
            actualPhysicalDepth=(depth_evidence[0]['supported'] and not depth_evidence[1]['supported']),
            slabCurrentOperands=np.linalg.norm(operand_joins['slab'])<1e-8,
            bulkCurrentOperands=np.linalg.norm(operand_joins['bulk'])<1e-8,
            nativeEVelocityPresence=all(item['nativeEVelocityPresent'] for item in face),
            undriven=face_zero,projectedBulk=bulk_zero,
            transverseCurrentHermitian=hermitian_supported,transverseCurrentNonzeroRank=current_rank_supported,
            controls=all(item['status']=='RESPONSIVE' for item in controls))
        if historical_comparison['status']=='COMPARED':
            checks['historicalWeightedGramComparison'] = historical_comparison['supported']
        return dict(status='AVAILABLE' if all(checks.values()) else 'UNRESOLVED',
            scope='SELECTED_SAVED_SEED_SUBSPACE_JOIN_NOT_COMPLETE_CENSUS',info=info,
            historicalEligibility=history,historicalWeightedGramComparison=historical_comparison,
            checks=checks,sourceMatrix=matrix,savedPencil=old['pencil'],sourceLift=lift,
            savedRight=right,transverseResidual=t_residual,pencilResidual=pencil_join,
            kernelResidual=kernel,waveResidual=wave,face=face,recomputedTransverseTransport=transport,
            actualNormalMomentum=k0,normalRealityResidual=normal_reality_residual,
            normalRealityCriterion='Zero imaginary part of the saved numeric coefficient, checked directly',
            currentOperandResiduals=operand_joins,projectedBulk=projected_bulk,controls=controls,
            savedDepthSignEvidence=depth_evidence)

    def classify(label, bound, w0, v0, prefix):
        if v0 == 0:
            record = dict(end=label,status='UNRESOLVED',coverage='UNRESOLVED',reason='CURL_CHART_ZERO_RAY', omega=str(w0),ray=str(v0),kappa='0')
            save(out/(prefix+'-summary.json'),record); return record
        point = op(prefix+'-census',census,bound,w0,v0,seconds=12)
        results = []
        if point.get('status')=='REAL_T_CANDIDATES' and point['coverage'] and zero(
                (point['realRemainder'],point['imaginaryRemainder'])):
            for index, interval in enumerate(point['intervals']):
                if time.monotonic()>journal.deadline-60: break
                results.append(op(prefix+'-root-%02d'%index,inspect_root,bound,point,interval,seconds=12))
            complete = (point['noRealExcludedPoint'] and len(results)==len(point['intervals'])
                        and all(x.get('status')=='AVAILABLE' and x.get('rootMultiplicityAccounted') is True for x in results))
            found = any(x.get('status')=='AVAILABLE' for x in results)
            state = 'AVAILABLE' if found else ('ABSENT' if not point['intervals'] and point['noRealExcludedPoint'] else 'UNRESOLVED')
            # A partially classified point retains existential availability but
            # never receives a complete inventory or inferred absence.
            coverage = 'COMPLETE' if complete else 'PARTIAL_UNRESOLVED'
        else:
            state, coverage = 'UNRESOLVED', 'UNRESOLVED'
        record = dict(end=label,status=state,coverage=coverage,omega=str(w0),ray=str(v0),
            kappa=str(sp.sqrt(5)*v0),candidateCount=point.get('count'),candidateMultiplicity=point.get('totalRootMultiplicity'),
            transverseDomainCertified=point.get('noRealExcludedPoint',False),
            candidateStatuses=[x.get('status','UNRESOLVED') for x in results],
            positiveCurrentRank=sum(x.get('positiveCurrentRank',0) for x in results),
            negativeCurrentRank=sum(x.get('negativeCurrentRank',0) for x in results),
            physicalLightCalibration='NOT_ESTABLISHED',lossFractionComputed=False)
        save(out/(prefix+'-summary.json'),record)
        return record

    def loci(bound):
        data = determinant_data(bound['restricted'])
        if data['numerator'].has(q,qb): return dict(status='UNRESOLVED', reason='DEPTH_DEPENDENT_T_LOCUS', **data)
        # No discriminant campaign: only zero-normal candidates and source
        # denominator/chart exclusions are proposed as exact equations.
        zero_normal = sp.factor(data['numerator'].subs(k,0))
        qp = sp.Poly(bound['q2'],k)
        bulk = qp.nth(0) if qp.degree()==2 and qp.nth(1)==0 and qp.nth(2).is_negative is True else None
        return dict(status='CANDIDATE_EQUATIONS', zeroNormal=zero_normal,
            sourceDenominator=data['denominator'], chart=bound['chart'],
            bulkMaximum=bulk, bulkBoundary=sp.factor(bulk) if bulk is not None else None,
            bulkAvailabilityDomain=sp.Gt(bulk,0) if bulk is not None else None,
            coordinateRelation=sp.Eq(sp.Symbol('kappa',nonnegative=True)**2,5*v*v),
            completePartition=False, **data)

    def slice_candidates(equation, v0, lower, upper):
        expression = sp.cancel(equation.subs(v,v0))
        numerator, denominator = sp.fraction(expression)
        real,imag = map(sp.expand,numerator.as_real_imag())
        rp,ip = sp.Poly(real,w,domain=sp.QQ),sp.Poly(imag,w,domain=sp.QQ)
        if rp.is_zero and ip.is_zero: return dict(status='UNRESOLVED', reason='IDENTICALLY_ZERO_SLICE')
        gcd = sp.gcd(rp,ip)
        intervals = gcd.intervals(eps=sp.Rational(1,10**18)) if gcd.degree()>0 else []
        intervals = [(interval,m) for interval,m in intervals if interval[1]>=lower and interval[0]<=upper]
        samples = []
        for index,((lo,hi),multiplicity) in enumerate(intervals):
            left_edge = lower if index==0 else intervals[index-1][0][1]
            right_edge = upper if index+1==len(intervals) else intervals[index+1][0][0]
            if not left_edge<lo<=hi<right_edge: continue
            samples.append(dict(interval=(lo,hi),multiplicity=multiplicity,
                                below=(left_edge+lo)/2,above=(hi+right_edge)/2))
        return dict(status='SLICE_CANDIDATES', equation=equation,ray=v0,
                    specialized=expression,denominator=denominator,intervals=intervals,samples=samples)

    frequencies = [sp.Rational(x) for x in spec['grid']['frequencies']]
    prefix_frequencies = [sp.Rational(x) for x in spec['grid']['prefixFrequencies']]
    require(all(x in frequencies for x in prefix_frequencies),'prefix uses existing frequency rows')
    transition_budget = spec['grid']['transitionPairBudget']
    require(isinstance(transition_budget,int) and 0<=transition_budget<=64,'bounded transition pairs')
    rays = [sp.Rational(x) for x in spec['grid']['rayCoordinates']]
    require(len(frequencies)*len(rays)<=256 and all(x>0 for x in frequencies)
            and all(x>=0 for x in rays), 'bounded rational-ray window')
    fixed = params['s11cdTangentialMomentum2']
    require(params['s11cdTangentialMomentum1']==2*fixed, 'saved azimuth join')
    ends, end_states, summaries, face_results, curves = {}, {}, [], {}, {}
    point_serial = 0
    visited_pairs = set()
    seed_results = {}
    save(out/'planned-grid.json',dict(ends=ENDS, pairs=[dict(omega=str(x),ray=str(y))
         for x in frequencies for y in rays], defaultUnvisitedState='UNRESOLVED',
         maximumDistinctPairs=spec['grid']['maxPairs']))
    def sample(label,bound,w0,v0):
        nonlocal point_serial
        pair_key = (str(w0),str(v0))
        if pair_key not in visited_pairs and len(visited_pairs)>=spec['grid']['maxPairs']:
            return dict(status='UNRESOLVED',reason='DISTINCT_PAIR_BUDGET')
        visited_pairs.add(pair_key)
        point_serial += 1
        if w0==params['omega'] and v0==fixed:
            selected = seed_results.get(label,[])
            found = any(item.get('status')=='AVAILABLE' for item in selected)
            result = dict(end=label,status='AVAILABLE' if found else 'UNRESOLVED',
                coverage='SELECTED_SAVED_SUBSPACES_ONLY',omega=str(w0),ray=str(v0),
                kappa=str(sp.sqrt(5)*v0),producerReplayed=False,newRootCensus=False,
                seedReceipts=[label+'-seed-join-'+str(index)+'.json' for index in (16,17)])
            save(out/('point-%04d-%s-summary.json'%(point_serial,label)),result)
        else:
            result = classify(label,bound,w0,v0,'point-%04d-%s'%(point_serial,label))
        summaries.append(result)
        return result

    # Useful pointwise output precedes all optional generic locus algebra.
    # No mandatory H polynomial or symbolic plane elimination exists here.
    for label in ENDS:
        if time.monotonic()>journal.deadline-150:
            end_states[label] = 'UNRESOLVED_BUDGET_NOT_ATTEMPTED'; continue
        packet = restored(label+'-restore-original-symbol',label+'Frequency')
        pairing = restored(label+'-restore-current',label+'Pairing')
        if any(isinstance(value,dict) and value.get('reason')=='OPERATION_BUDGET'
               for value in (packet,pairing)):
            end_states[label] = 'UNRESOLVED_SOURCE_RESTORE_BUDGET'; continue
        bound = op(label+'-source-bind',bound_end,label,packet,pairing,seconds=35)
        if 'checks' not in bound or not all(bound['checks'].values()):
            end_states[label] = 'UNRESOLVED_SOURCE_JOIN'; continue
        face = op(label+'-face-premise',faces,bound,seconds=20)
        face_results[label] = face
        if 'checks' not in face or not all(face['checks'].values()):
            end_states[label] = 'UNRESOLVED_LOSSLESS_T_PREMISE'; continue
        ends[label] = bound; end_states[label] = 'SUPPORTED_SOURCE_PREMISE'
        seed_results[label] = []
        for index in (16,17):
            key = label+'Seed'+str(index)
            if key not in spec['inputs']:
                seed_results[label].append(dict(status='UNRESOLVED',reason='NO_PINNED_SAVED_SEED_OPERAND',key=key))
                continue
            saved = restored(label+'-restore-seed-'+str(index),key)
            if isinstance(saved,dict) and saved.get('reason')=='OPERATION_BUDGET':
                seed_results[label].append(saved); continue
            seed_results[label].append(op(label+'-seed-join-'+str(index),saved_seed,bound,saved,seconds=12))
        for w0 in prefix_frequencies:
            sample(label,bound,w0,fixed)
    save(out/'source-end-states.json',end_states)
    save(out/'reference-face-drive.json',show(face_results.get('REFERENCE',{'status':'UNRESOLVED'})))
    seed_hashes = {}
    for name,record in spec['inputs'].items():
        if any(name==label+'Seed'+str(index) for label in ENDS for index in (16,17)):
            seed_hashes.setdefault(record['sha256'],[]).append(name)
    save(out/'saved-seed-comparison.json',show(dict(selectedRecords=seed_results,
        savedSeedProvenance=spec.get('savedSeedProvenance'),sharedAcceptedHashes=seed_hashes,
        sharedBytesAreIndependentEvidence=False,completeCensusEstablished=False,producerReplayed=False)))
    transition_count = 0
    for label,bound in ends.items():
        if time.monotonic()>journal.deadline-150:
            curves[label] = {'status':'UNRESOLVED_BUDGET_NOT_ATTEMPTED'}; continue
        curve = op(label+'-optional-loci',loci,bound,seconds=20); curves[label] = curve
        if curve.get('status')!='CANDIDATE_EQUATIONS': continue
        for key in ('zeroNormal','bulkBoundary'):
            if curve.get(key) is None: continue
            crossings = op(label+'-'+key+'-slice',slice_candidates,curve[key],fixed,
                           min(frequencies),max(frequencies),seconds=10)
            if crossings.get('status')!='SLICE_CANDIDATES': continue
            for crossing in crossings['samples']:
                if transition_count+2>transition_budget or time.monotonic()>journal.deadline-110: break
                below = sample(label,bound,crossing['below'],fixed)
                above = sample(label,bound,crossing['above'],fixed)
                values = {side:sp.cancel(curve[key].subs({w:crossing[side],v:fixed})) for side in ('below','above')}
                # Actual sides are persisted even when no physical crossing is
                # found; a candidate equation is never promoted by its name.
                op('transition-%03d'%transition_count,lambda x:x,
                   dict(end=label,locus=key,candidate=crossing,below=below,above=above,
                        equationValues=values,physicalCutoffEstablished=False))
                transition_count += 2
    save(out/'threshold-loci.json',show(curves))
    stop = None
    ordered = [(x,fixed) for x in frequencies if x not in prefix_frequencies]+[(x,y) for x in frequencies for y in rays if y!=fixed]
    for w0,v0 in ordered:
        if time.monotonic()>journal.deadline-100:
            stop = 'BUDGET_SAVED_PREFIX'; break
        for label,bound in ends.items():
            if time.monotonic()>journal.deadline-60:
                stop = 'BUDGET_SAVED_PREFIX'; break
            sample(label,bound,w0,v0)
    seen = {(item['end'],item['omega'],item['ray']) for item in summaries if 'end' in item}
    unvisited = [dict(end=label,omega=str(x),ray=str(y),status='UNRESOLVED',reason='NOT_VISITED')
                 for label in ENDS for x in frequencies for y in rays
                 if (label,str(x),str(y)) not in seen]
    save(out/'unvisited-grid.json',unvisited)
    save(out/'grid-summary.json',summaries)
    return dict(status='TRANSVERSE_FACE_SAVED_CANDIDATE_REQUIRES_INSPECTION',case=spec['case'],
        endStates=end_states,points=len(summaries),gridStop=stop,completePlanePartition=False,
        referenceOnlyFaceScope=True,matchedAllGradeEndPremiseEstablished=False,
        firstOrderBoundaryMapsEstablished=False,totalLossComputed=False,
        thicknessClassification='DEFERRED',physicalLightCalibration='NOT_ESTABLISHED',
        azimuth='2:1; kappa=sqrt(5)*v',sourceEndMismatchMayRemain=True)


def main():
    parser = argparse.ArgumentParser(__doc__)
    for name in ('input-manifest','gate-receipt','run-directory'):
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    spec = json.loads(args.input_manifest.read_text()); gate = json.loads(args.gate_receipt.read_text())
    require(gate.get('status')=='READY_FOR_ONE_GUARDED_TRANSVERSE_FACE_PREMISE_JOB', 'fresh execution gate required')
    require(gate.get('scienceAuthorizationAfterCostedPlanGo') is True
            and gate.get('userScienceApprovalRecord'), 'current preparation approval authorizes no science')
    approval = gate['userScienceApprovalRecord']
    require(route(approval['path'])==approval, 'post-plan explicit science approval pin')
    require(gate['workerSha256']==sha(__file__) and gate['inputManifestSha256']==sha(args.input_manifest), 'worker/input pins')
    require(gate['independentBuildClearance'] is True and gate['scienceJobOrdinal']==1, 'review/one-job gate')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False, 'normal duration/no retry')
    require(gate['guardSha256']==sha(ROOT/'scripts/s11c_guarded_run.py') and gate['supervisorSha256']==sha(
        ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'), 'shared guard/supervisor pins')
    for name,record in spec['inputs'].items(): require(route(record['path'])==record,'changed input '+name)
    require(json.loads(Path(spec['inputs']['physicalInputFile']['path']).read_text())==spec['physicalInput'], 'physical inputs unchanged')
    for join in spec['checkpointJoins']:
        checkpoint = json.loads(Path(join['checkpoint']).read_text())
        require(checkpoint['artifacts'][join['artifact']]['sha256']==spec['inputs'][join['inputKey']]['sha256'], 'accepted source checkpoint')
    limits = containment(); started = time.monotonic()
    out = args.run_directory.resolve(); out.relative_to(ROOT/'_scratch/s11c'); out.mkdir(parents=True,exist_ok=False)
    journal = Journal(out,started+840)
    save(out/'native-limits.json',limits); save(out/'manifest.json',spec); save(out/'gate.json',gate)
    try:
        result = construct(spec,out,journal)
        post = {name:route(record['path']) for name,record in spec['inputs'].items()}
        save(out/'posthashes.json',post); require(post==spec['inputs'],'posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started,operations=len(journal.records),allSourcesUnchanged=True)
        save(out/'operation-index.json',journal.records); save(out/'checks.json',result)
        print(json.dumps(result,indent=2,allow_nan=False))
    except BaseException:
        signal.setitimer(signal.ITIMER_REAL,0)
        save(out/'failure.json',dict(traceback=traceback.format_exc(),incompleteOperation=journal.active,
             wallSeconds=time.monotonic()-started,noRetry=True))
        save(out/'operation-index.json',journal.records)
        save(out/'failure-posthashes.json',{name:route(record['path']) for name,record in spec['inputs'].items()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        save(out/'artifact-index.json',{str(path.relative_to(out)):route(path)
             for path in sorted(out.rglob('*')) if path.is_file()})


if __name__ == '__main__':
    main()
