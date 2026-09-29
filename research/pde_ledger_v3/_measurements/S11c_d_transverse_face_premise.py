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


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


class OperationBudget(BaseException):
    """Deliberately escapes ordinary algebra/domain exception handlers."""


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.')) or
                module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError((module, name))
        return super().find_class(module, name)


class Journal:
    """Each operation has complete operands; repeated objects use saved refs.

    The identity table holds strong references, so Python id reuse cannot alias
    an operand. It references only complete immutable-by-convention returns.
    Every attempted operation has its own receipt even when values coincide.
    """
    def __init__(self, out, deadline):
        self.out, self.deadline = out, deadline
        self.records, self.active, self.objects = [], None, {}

    def blob(self, value):
        existing = self.objects.get(id(value))
        if existing is not None and existing[0] is value:
            return existing[1]
        payload = pickle.dumps(value, protocol=4)
        digest = hashlib.sha256(payload).hexdigest()
        path = self.out/'objects'/(digest+'.pickle')
        path.parent.mkdir(exist_ok=True)
        if not path.exists():
            with path.open('xb') as stream:
                stream.write(payload); stream.flush(); os.fsync(stream.fileno())
        record = route(path)
        self.objects[id(value)] = (value, record)
        return record

    def op(self, name, function, *args, seconds=20):
        path = self.out/'operations'/('%04d-%s' % (len(self.records), name))
        path.mkdir(parents=True, exist_ok=False)
        self.active = dict(name=name, operands=[self.blob(arg) for arg in args],
                           startedUtc=datetime.now(timezone.utc).isoformat())
        save(path/'input-references.json', self.active)
        remaining = self.deadline-time.monotonic()
        if remaining <= 0:
            raise OperationBudget('native whole-job deadline')
        signal.setitimer(signal.ITIMER_REAL, min(seconds, remaining))
        try:
            value = function(*args)
            record = dict(self.active, status='COMPLETE', value=self.blob(value))
            save(path/'complete.json', record)
            self.records.append(record); self.active = None
            return value
        except OperationBudget:
            record = dict(self.active, status='BUDGET_STOP_UNRESOLVED_NO_RETRY')
            save(path/'budget-stop.json', record)
            self.records.append(record); self.active = None
            if time.monotonic() >= self.deadline:
                raise
            return {'status': 'UNRESOLVED', 'reason': 'OPERATION_BUDGET',
                    'operationReceipt': str(path/'budget-stop.json')}
        finally:
            signal.setitimer(signal.ITIMER_REAL, max(.001, self.deadline-time.monotonic()))


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

    def op(name, function, *args, seconds=20):
        value = journal.op(name, function, *args, seconds=seconds)
        save(out/(name+'.json'), show(value))
        return value

    def restore(name):
        record = spec['inputs'][name]
        require(route(record['path']) == record, 'changed source '+name)
        with Path(record['path']).open('rb') as stream:
            return SavedCodec(stream).load()

    def flat(value):
        if isinstance(value, sp.MatrixBase): return list(value)
        if isinstance(value, dict): return [z for x in value.values() for z in flat(x)]
        if isinstance(value, (tuple, list)): return [z for x in value for z in flat(x)]
        return [value]

    def zero(value): return all(sp.cancel(x) == 0 for x in flat(value))
    def clean(matrix): return sp.ImmutableMatrix(matrix).applyfunc(sp.cancel)

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

    uniform = op('restore-uniform', restore, 'uniform')
    reduction = op('restore-branch-context', restore, 'reduction')
    units = op('restore-unit-context', restore, 'unitContext')
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
    save(out/'native-face-order-source.json',dict(order=face_order,
        source=spec['inputs']['sourceEngine'],sourceLine=sign_loops[0].lineno,
        expression=ast.unparse(sign_loops[0].iter)))

    def branch_inventory(state):
        equations, mapping = state['branch_equations'], state['branch_map']
        residuals = [equation.rhs-mapping[equation.lhs] for equation in equations]
        return dict(equations=equations, mapping=mapping, residuals=residuals,
                    allEquationKeysPresent=all(eq.lhs in mapping for eq in equations),
                    joins=zero(residuals), tangents=state['tangents'],
                    groups=state['momentum_groups'], normalMap=state['normal_map'])

    branches = op('source-branch-map-joins', branch_inventory, reduction['reductionState'])
    require(branches.get('joins') is True, 'saved branch equation/map inconsistency')

    def bound_end(label, packet, pair_tuple):
        pair, dimensions = pair_tuple  # Authoritative end_pairing_check.py:165-166.
        result = pair['result']
        wl, wr = result['FREQUENCY_LEGS']; kl, kr = result['NORMAL_LEGS']
        ql, qr = result['BULK_LEGS']
        mapping = {wl:w, wr:w, kl:k, kr:k, ql:qb, qr:q}
        matrix = clean(bind(result['CLOSED_PENCIL_LEGS'][0], label, mapping))
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
        lift = sp.ImmutableMatrix.vstack(curl[:, 1:3], sp.zeros(2, 2))
        gram = clean(lift.H*lift)
        # A 2x2 coordinate inverse only; never invert the strong physical pencil.
        weak = clean(lift.H*matrix*lift)
        restricted = clean(gram.inv()*weak)
        invariant = clean(matrix*lift-lift*restricted)
        currents, grade_residuals = {}, {}
        for name in ('SLAB_CURRENT_MATRIX', 'BULK_NORMAL_CURRENT_DENSITY_MATRIX',
                     'BULK_DEPTH_CURRENT_MATRIX', 'INTERFACE_POWER_MATRIX'):
            value, residual = grade(result[name])
            currents[name] = bind(value, label, mapping)
            grade_residuals[name] = residual
        amplitudes = sorted((s for s in dimensions if getattr(s, 'name', '').startswith(
            's11cdCurrentPlusAmplitude')), key=lambda s:s.name)
        require(len(amplitudes) == 5, 'five physical amplitude columns required')
        rows, reconstruction = [], []
        for face in result['FACE_LEG_OBJECTS']:
            row, residual = {}, {}
            for name in DRIVES:
                expression = face[0][name]
                coefficients = sp.ImmutableMatrix(1, 5, lambda i,j:sp.diff(expression, amplitudes[j]))
                row[name] = bind(coefficients, label, mapping)
                residual[name] = expression-(coefficients*sp.ImmutableMatrix(amplitudes))[0]
            rows.append(row); reconstruction.append(residual)
        pressure_scale = bind(result['OPEN_BULK_PRESSURE_SCALES'][0],label,mapping)
        velocity_scale = bind(result['OPEN_BULK_VELOCITY_SCALES'][0],label,mapping)
        flat_exterior = [dict(pressure=row['PRESSURE'], amplitude=row['AMPLITUDE'],
            bulkVelocity=row['BULK_VELOCITY'], pressureScale=pressure_scale,
            velocityScale=velocity_scale,
            pressureResidual=clean(row['PRESSURE']-pressure_scale*row['AMPLITUDE']),
            velocityResidual=clean(row['BULK_VELOCITY']-velocity_scale*row['AMPLITUDE']))
            for row in rows]
        depth = result['OPEN_BULK_CURRENT_COEFFICIENTS'][3]
        e = next(s for s in depth.free_symbols if s.name == 'epsilon_shape')
        free_amps = {s:sp.S.One for s in depth.free_symbols
                     if s.name in ('s11cdAcousticLeftAmplitude', 's11cdAcousticRightAmplitude')}
        depth = bind(sp.Poly(depth,e).nth(2), label, {**mapping, **free_amps})
        denominators = tuple(dict.fromkeys(sp.denom(sp.cancel(x)) for x in
            [*matrix, *restricted, *gram, *[z for m in currents.values() for z in m],
             *[z for row in rows for m in row.values() for z in m]]))
        coupling = [uniform['records'][label]['coupling'][name] for name in ('TH','HT')]
        checks = dict(sourceJoin=zero(source_join), waveJoin=zero(wave_join),
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
            sourceCoupling=coupling, matrix=matrix, wave=wave, q2=q2, scale=scale,
            lift=lift, gram=gram, chart=sp.factor(gram.det()), restricted=restricted, weak=weak,
            invariantResidual=invariant, weakHermitianResidual=weak-weak.H,
            currents=currents, gradeResiduals=grade_residuals, drives=rows,
            driveReconstruction=reconstruction, flatExteriorJoins=flat_exterior, sourceDenominators=denominators,
            outgoingDepthCoefficient=depth, dimensions=dimensions)

    def faces(bound):
        lift = bound['lift']; records = []
        literals = []
        if bound['end'] == 'REFERENCE':
            for item in face_index['records']:
                if item['tag'] != 'PY_S11CC2_FOLD_SYMBOL_MAP_LAB_HELD_RHO4_CONSTANT': continue
                for identity in item['velocityIdentifications']:
                    text = identity['savedValueSrepr']
                    require(hashlib.sha256(text.encode()).hexdigest() == identity['literalSha256'],
                            'c2 literal hash differs')
                    expression = sp.sympify(text)
                    field = next(s for s in expression.free_symbols if s.name == 'e_W_t')
                    literals.append(dict(face=item['face'], source=expression,
                        coefficient=bind(expression, 'REFERENCE', {field:-sp.I*w}),
                        literalSha256=identity['literalSha256']))
        for index, row in enumerate(bound['drives']):
            contractions = {name:clean(row[name]*lift) for name in DRIVES}
            original = row['OUTWARD_VELOCITY']
            changed = sp.MutableDenseMatrix(original); changed[0,4] = 0
            probe = sp.eye(5)[:,4]
            difference = clean((changed-original)*probe)
            # Tilting a supplied transverse vector into each native dependent
            # physical slot tests the zero-drive restriction, not bulk power.
            dependencies = []
            for name in DRIVES:
                for slot in range(5):
                    coefficient = sp.cancel(row[name][0,slot])
                    if coefficient.is_zero is False:
                        tilted = lift[:,0]+sp.eye(5)[:,slot]
                        dependencies.append(dict(drive=name, slot=slot, coefficient=coefficient,
                            originalVector=lift[:,0], tiltedVector=tilted,
                            baseline=clean(row[name]*lift[:,0]), changed=clean(row[name]*tilted),
                            movement=clean(row[name]*(tilted-lift[:,0]))))
            native_controls = []
            for name in DRIVES:
                terms = [(slot, term) for slot in range(5)
                         for term in sp.Add.make_args(sp.expand(row[name][0,slot]))
                         if term != 0]
                if terms:
                    slot, term = terms[0]
                    mutated = sp.MutableDenseMatrix(row[name]); mutated[0,slot] -= term
                    loaded = sp.eye(5)[:,slot]
                    native_controls.append(dict(drive=name, slot=slot, omittedTerm=term,
                        originalRow=row[name], changedRow=mutated, loadedProbe=loaded,
                        original=row[name]*loaded, changed=mutated*loaded,
                        difference=clean((mutated-row[name])*loaded)))
            addressed_literals = [item for item in literals if item['face']==face_order[index]]
            joins = [sp.cancel(original[0,4]-item['coefficient']) for item in addressed_literals]
            records.append(dict(face=index,sourceFaceOrientation=face_order[index], nativeRows=row, lift=lift, contractions=contractions,
                driveZero=zero(contractions), eWControl=dict(original=original, changed=changed,
                    probe=probe, difference=difference, responsive=not zero(difference)),
                nativeTermControls=native_controls, dependentSlotProbes=dependencies,
                dependencyResponsive=bool(dependencies) and all(not zero(x['movement']) for x in dependencies),
                c2VelocityJoins=joins, c2VelocityJoinSupported=(len(addressed_literals)==1) if bound['end']=='REFERENCE' else None))
        forms = {name:clean(lift.H*matrix*lift) for name,matrix in bound['currents'].items()
                 if name != 'SLAB_CURRENT_MATRIX'}
        checks = dict(allDrivesZero=all(x['driveZero'] for x in records),
            eWControls=all(x['eWControl']['responsive'] for x in records),
            dependentSlotControls=all(x['dependencyResponsive'] for x in records),
            nativeControls=all(x['nativeTermControls'] and all(not zero(z['difference'])
                for z in x['nativeTermControls']) for x in records),
            c2Velocity=(bound['end']!='REFERENCE' or
                all(x['c2VelocityJoinSupported'] and zero(x['c2VelocityJoins']) for x in records)),
            lossSideFormsZero=zero(forms))
        return dict(end=bound['end'], faces=records, physicalForms=forms, checks=checks,
                    c2LiteralEvidence=literals, flatReferenceExteriorJoins=bound['flatExteriorJoins'],
                    c2PressureTraceFirstShapeJoin='NOT_SUPPORTED_BY_THIS_INPUT_INDEX')

    def determinant_data(matrix):
        # Specialize the tiny T matrix before arithmetic. No full/H determinant,
        # radical elimination, resultants, or discriminants are computed.
        determinant = sp.cancel(matrix.det())
        numerator, denominator = sp.fraction(determinant)
        return dict(matrix=matrix, determinant=determinant, numerator=numerator, denominator=denominator)

    def census(bound, w0, v0):
        specialized = clean(bound['restricted'].subs({w:w0,v:v0}))
        data = determinant_data(specialized)
        if specialized.has(q, qb): return dict(status='UNRESOLVED', reason='T_RESTRICTION_RETAINS_DEPTH_ROOT', **data)
        real, imag = map(sp.expand, data['numerator'].as_real_imag())
        try:
            rp, ip = sp.Poly(real,k,domain=sp.QQ), sp.Poly(imag,k,domain=sp.QQ)
        except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed):
            return dict(status='UNRESOLVED', reason='UNSUPPORTED_POINT_POLYNOMIAL_DOMAIN', **data)
        if rp.is_zero and ip.is_zero: return dict(status='UNRESOLVED', reason='DEGENERATE_T_DETERMINANT', **data)
        gcd = sp.gcd(rp, ip).sqf_part()
        intervals = gcd.intervals(eps=sp.Rational(1,10**28)) if gcd.degree()>0 else []
        count = gcd.count_roots(-sp.oo,sp.oo) if gcd.degree()>0 else 0
        return dict(status='REAL_T_CANDIDATES', omega=w0, ray=v0, kappa=sp.sqrt(5)*v0,
            realPolynomial=rp.as_expr(), imaginaryPolynomial=ip.as_expr(), gcd=gcd.as_expr(),
            realRemainder=sp.rem(rp,gcd).as_expr(), imaginaryRemainder=sp.rem(ip,gcd).as_expr(),
            intervals=intervals, count=count,
            coverage=(count == sum(m for _,m in intervals)), **data)

    def no_pole_interval(expression, bound, point, interval, depth_sign):
        lo, hi = interval
        specialized = sp.cancel(expression.subs({w:point['omega'],v:point['ray'],qb:depth_sign*q}))
        numerator, denominator = sp.fraction(specialized)
        certificates = []
        q2 = sp.cancel(bound['q2'].subs({w:point['omega'],v:point['ray']}))
        for value in (numerator,denominator):
            reduced = sp.rem(sp.Poly(value,q),sp.Poly(q*q-q2,q)).as_expr()
            norm = sp.cancel(reduced*reduced.xreplace({q:-q}))
            norm = sp.rem(sp.Poly(sp.fraction(norm)[0],q),sp.Poly(q*q-q2,q)).as_expr()
            real_norm = sp.cancel(norm*sp.conjugate(norm))
            polys = [sp.Poly(x,k,domain=sp.QQ) for x in sp.fraction(real_norm)]
            counts = [None if p.is_zero else (p.count_roots(lo,hi) if p.degree()>0 else 0) for p in polys]
            certificates.append(dict(reduced=reduced, norm=norm, realNorm=real_norm,
                                     polynomialNumerator=polys[0].as_expr(), polynomialDenominator=polys[1].as_expr(),
                                     intervalRootCounts=counts))
        return dict(expression=expression, specialized=specialized, certificates=certificates,
                    excluded=any(c != 0 for item in certificates for c in item['intervalRootCounts']))

    def numeric(matrix, mapping):
        return np.asarray(matrix.subs(mapping).evalf(40).tolist(), dtype=complex)

    def loaded_control(bound, symbol, mapping, basis, scale):
        matrix = bound['matrix']; selected = None; best = 0.
        for i in range(matrix.rows):
            for j in range(matrix.cols):
                if np.linalg.norm(basis[j,:])<=1e-12 or not matrix[i,j].has(symbol): continue
                for term in sp.Add.make_args(sp.expand(matrix[i,j])):
                    if not term.has(symbol): continue
                    amount = complex(term.subs(mapping).evalf(40))
                    movement = abs(amount)*float(np.linalg.norm(basis[j,:]))/scale[i]
                    if movement > best:
                        best = movement; selected = (i,j,term)
        if selected is None: return dict(status='UNRESOLVED', reason='NO_LOADED_DEPENDENT_SOURCE_TERM')
        i,j,term = selected; changed = sp.MutableDenseMatrix(matrix); changed[i,j] -= term
        baseline_matrix, changed_matrix = numeric(matrix,mapping), numeric(changed,mapping)
        original, altered = (baseline_matrix/scale[:,None])@basis, (changed_matrix/scale[:,None])@basis
        return dict(status='RESPONSIVE' if best>1e-8 else 'UNRESOLVED', symbol=symbol,
            selectedEntry=(i,j), omittedNativeTerm=term, baselineMatrix=baseline_matrix,
            mutatedMatrix=changed_matrix, baselineResidual=original, mutatedResidual=altered,
            difference=altered-original, norm=float(np.linalg.norm(altered-original)))

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
        q0 = sp.sqrt(q2)
        mapping = {w:point['omega'],v:point['ray'],k:km,q:q0,qb:sp.conjugate(q0)}
        matrix = numeric(bound['matrix'],mapping); lift = numeric(bound['lift'],mapping)
        restricted = numeric(bound['restricted'],mapping)
        row_scale = np.maximum(np.linalg.norm(matrix,axis=1),1.)
        _, singular, vh = np.linalg.svd(restricted)
        tolerance = 1e-10*max(1.,float(np.linalg.norm(restricted)))
        null = singular<tolerance
        base.update(physicalDepth=q0, depthSquared=q2, exclusions=exclusions,
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
        face = [{name:numeric(row[name],mapping)@basis for name in DRIVES} for row in bound['drives']]
        current = forms['SLAB_CURRENT_MATRIX']; hermitian = current-current.conj().T
        eig, rotation = np.linalg.eigh((current+current.conj().T)/2)
        depth_coefficient = complex(bound['outgoingDepthCoefficient'].subs(mapping).evalf(40))
        controls = [loaded_control(bound,symbol,mapping,basis,row_scale) for symbol in (w,v)]
        base.update(basis=basis, fullResidual=residual, rowScale=row_scale, forms=forms,
                    face=face, current=current, hermitianResidual=hermitian, currentEigenvalues=eig,
                    rotation=rotation, controls=controls, outgoingDepthCoefficient=depth_coefficient)
        if q2.is_positive and (abs(depth_coefficient.imag)>1e-10 or depth_coefficient.real<=0):
            return dict(base, reason='OUTGOING_DEPTH_SIGN_UNRESOLVED')
        if any(np.linalg.norm(x)>1e-8 for row in face for x in row.values()) or any(
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
        face = [{name:numeric(row[name],mapping)@right for name in DRIVES}
                for row in bound['drives']]
        slab = numeric(bound['currents']['SLAB_CURRENT_MATRIX'],mapping)
        bulk = numeric(bound['currents']['BULK_NORMAL_CURRENT_DENSITY_MATRIX'],mapping)
        current = right.conj().T@slab@right
        current_join = current-np.asarray(old['currentGram'],complex) if 'currentGram' in old else None
        operand_joins = {'slab':slab-np.asarray(old['currentOperands']['CURRENT_SLAB'],complex),
                         'bulk':bulk-np.asarray(old['currentOperands']['CURRENT_BULK'],complex)}
        controls = [loaded_control(bound,symbol,mapping,right,scale) for symbol in (w,v)]
        wave = complex(bound['wave'].subs(mapping).evalf(40))
        projected_bulk = right.conj().T@bulk@right
        eligible = bool(info['SHEET_MEMBERSHIP'] and info['EXACT_REAL_NORMAL'] and old['currentDefined'])
        checks = dict(savedEligibility=eligible, savedFrequency=(w0==params['omega']),
            transverseMembership=np.linalg.norm(t_residual)<1e-8,
            physicalPencilJoin=np.linalg.norm(pencil_join)<1e-8,
            physicalKernel=np.linalg.norm(kernel)<1e-8, physicalWave=abs(wave)<1e-8,
            slabCurrentOperands=np.linalg.norm(operand_joins['slab'])<1e-8,
            bulkCurrentOperands=np.linalg.norm(operand_joins['bulk'])<1e-8,
            undriven=all(np.linalg.norm(value)<1e-8 for row in face for value in row.values()),
            projectedBulk=np.linalg.norm(projected_bulk)<1e-8,
            currentGram=('currentGram' in old and np.linalg.norm(current_join)<1e-8),
            controls=all(item['status']=='RESPONSIVE' for item in controls))
        return dict(status='AVAILABLE' if all(checks.values()) else 'UNRESOLVED',
            scope='SELECTED_SAVED_SEED_SUBSPACE_JOIN_NOT_COMPLETE_CENSUS', info=info,
            checks=checks, sourceMatrix=matrix, savedPencil=old['pencil'], sourceLift=lift,
            savedRight=right, transverseResidual=t_residual, pencilResidual=pencil_join,
            kernelResidual=kernel, waveResidual=wave, face=face, current=current,
            savedCurrent=old.get('currentGram'), currentResidual=current_join,
            currentOperandResiduals=operand_joins, projectedBulk=projected_bulk, controls=controls)

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
            complete = len(results)==len(point['intervals']) and all(x.get('status')=='AVAILABLE' for x in results)
            found = any(x.get('status')=='AVAILABLE' for x in results)
            state = 'AVAILABLE' if found else ('ABSENT' if not point['intervals'] else 'UNRESOLVED')
            # A partially classified point retains existential availability but
            # never receives a complete inventory or inferred absence.
            coverage = 'COMPLETE' if complete else 'PARTIAL_UNRESOLVED'
        else:
            state, coverage = 'UNRESOLVED', 'UNRESOLVED'
        record = dict(end=label,status=state,coverage=coverage,omega=str(w0),ray=str(v0),
            kappa=str(sp.sqrt(5)*v0),candidateCount=point.get('count'),
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
        gcd = sp.gcd(rp,ip).sqf_part()
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
        packet = op(label+'-restore-original-symbol',restore,label+'Frequency')
        pairing = op(label+'-restore-current',restore,label+'Pairing')
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
            saved = op(label+'-restore-seed-'+str(index),restore,key)
            if isinstance(saved,dict) and saved.get('reason')=='OPERATION_BUDGET':
                seed_results[label].append(saved); continue
            seed_results[label].append(op(label+'-seed-join-'+str(index),saved_seed,bound,saved,seconds=12))
        for w0 in frequencies[:3]:
            sample(label,bound,w0,fixed)
    save(out/'source-end-states.json',end_states)
    save(out/'reference-face-drive.json',show(face_results.get('REFERENCE',{'status':'UNRESOLVED'})))
    save(out/'saved-seed-comparison.json',show(dict(selectedRecords=seed_results,
        suppliedFiniteResponseSummary=spec.get('savedSeedSummary'),
        completeCensusEstablished=False, producerReplayed=False)))
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
                if transition_count>=24 or time.monotonic()>journal.deadline-110: break
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
    ordered = [(x,fixed) for x in frequencies[3:]]+[(x,y) for x in frequencies for y in rays if y!=fixed]
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
