#!/usr/bin/env python3
"""Bounded reference inverse and regular spectral density from saved operands.

Imports no scientific producer. Production requires the whole-job resource
guard and an explicit implementation-gate receipt. No profile response,
integration, mode construction, pole search, or automatic retry is performed.
An inverse symbol is never labelled a complete outgoing Green operator.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import pickle
import resource
import signal
import sys
import time
import traceback

M = Path(__file__).resolve().parent
PROJECT = M.parents[2]
STORE = PROJECT / '_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            value.update(block)
    return value.hexdigest()


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def route(path):
    path = Path(path)
    return {'path': str(path), 'canonicalPath': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'sha256': digest(path)}


def containment():
    group = next(s[3:] for s in Path('/proc/self/cgroup').read_text().splitlines()
                 if s.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    observed = {n: (root / n).read_text().strip()
                for n in ('memory.max', 'memory.swap.max', 'pids.max')}
    observed.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                    affinity=sorted(os.sched_getaffinity(0)),
                    threads={n: os.environ.get(n) for n in THREADS})
    if not (observed['memory.max'] == str(2 * 1024**3)
            and observed['memory.swap.max'] == '0'
            and observed['pids.max'] == '32' and observed['nice'] >= 15
            and len(observed['affinity']) == 1
            and all(v == '1' for v in observed['threads'].values())):
        raise RuntimeError('Required whole-job containment missing; no science imported')
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise TimeoutError('900-second native limit; keep completed operations')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(900)
    return observed


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.'))
                or module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError(('unapproved codec module', module, name))
        return super().find_class(module, name)


class Journal:
    def __init__(self, out):
        self.out = out
        self.records = []

    def value(self, relative, value):
        path = self.out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            pickle.dump(value, stream, protocol=4)
            stream.flush()
            os.fsync(stream.fileno())
        return route(path)

    def op(self, name, function, *args):
        folder = 'operations/%03d-%s' % (len(self.records), name)
        item = {'name': name, 'input': self.value(folder + '/input.pickle', args),
                'startedUtc': datetime.now(timezone.utc).isoformat()}
        save_json(self.out / folder / 'started.json', item)
        # Persist inputs before invoking the actual operation and its result
        # before any later guard; a failed call leaves its input and start record.
        value = function(*args)
        item['value'] = self.value(folder + '/value.pickle', value)
        item['finishedUtc'] = datetime.now(timezone.utc).isoformat()
        save_json(self.out / folder / 'completed.json', item)
        self.records.append(item)
        return value


def sign(permutation):
    return -1 if sum(a > b for i, a in enumerate(permutation)
                     for b in permutation[i+1:]) % 2 else 1


def determinant(rows):
    """Explicit finite scalar expression, with no matrix-inverse placeholder."""
    n = len(rows)
    return sum(sign(p) * math.prod(rows[i][p[i]] for i in range(n))
               for p in itertools.permutations(range(n)))


def adjugate_entry(rows, i, j):
    # Transposed cofactor: remove source row j and source column i.
    minor = [[v for b, v in enumerate(row) if b != i]
             for a, row in enumerate(rows) if a != j]
    return (-1)**(i+j) * determinant(minor)


def require(test, label):
    if not test:
        raise ValueError(label)


def construct(spec, out):
    import sympy as sp
    import mpmath as mp
    journal = Journal(out)
    actual = []
    for expected in spec['files']:
        record = route(expected['path'])
        save_json(out / 'consumed' / ('%03d.json' % len(actual)), record)
        actual.append(record)
        require(record == expected, ('consumed file changed', record['path']))
    by_path = {r['path']: r for r in actual}

    def read_saved(key):
        path = spec[key]
        require(path in by_path, ('unmanifested packet', key))
        with Path(path).open('rb') as stream:
            return SavedCodec(stream).load()

    packet = journal.op('restore-uniform', read_saved, 'uniformPacket')
    reduction = journal.op('restore-reduction', read_saved, 'reductionPacket')
    physical = json.loads(Path(spec['physicalInput']).read_text())
    checkpoint = json.loads(Path(spec['endCheckpoint']).read_text())
    require(checkpoint['status'] == 'ACCEPTED_CASE_END_AND_CURRENT_INPUT_SOURCES',
            'accepted case-end source status')
    key = 'cases/LAB_HELD__RHO4_CONSTANT/uniform-source.pickle'
    require(str(Path(checkpoint['runDirectory']) / key) == spec['uniformPacket'],
            'checkpoint-selected own baseline packet')
    require(checkpoint['artifacts'][key]['sha256'] == by_path[spec['uniformPacket']]['sha256'],
            'accepted uniform hash')
    require(checkpoint['inputPackets'][spec['reductionPacket']] ==
            by_path[spec['reductionPacket']]['sha256'], 'accepted reduction hash')
    journal.value('selected-inputs.pickle', {
        'reference': packet['records']['REFERENCE'], 'units': packet['units'],
        'constantFourierMass': packet['constantFourierMass'],
        'reductionState': reduction['reductionState'], 'physicalInput': physical})
    source = packet['records']['REFERENCE']['strong']
    require(isinstance(source, sp.MatrixBase) and source.shape == (5, 5),
            'actual physical strong symbol, not six-component weak gauge matrix')
    require(not any(packet['records']['REFERENCE']['unresolved'].values()),
            'saved reference symbol is resolved')
    state = reduction['reductionState']
    kns = [s for s in source.free_symbols if s.name == 's11cdSpectralNormalMomentum']
    require(len(kns) == 1, 'unique actual spectral coordinate')
    kn = kns[0]
    require(state['z'] != state['zp'], 'distinct source and observation coordinates')
    known = packet['dimensionState']['known']
    field_units = []
    for name in ('u1', 'u2', 'u3', 'theta', 'eW'):
        matches = [u for f, u in known.items()
                   if getattr(f, '__name__', None) == 's11cdReducedField' + name]
        require(len(matches) == 1, ('native field-unit binding', name))
        field_units.append(tuple(matches[0]))
    entry_units = [[tuple(packet['units']['strong'][(5*i+j,)]) for j in range(5)]
                   for i in range(5)]
    row_units = [tuple(a+b for a, b in zip(entry_units[i][0], field_units[0]))
                 for i in range(5)]
    unit_residuals = [[tuple(a-b-c for a, b, c in
                            zip(row_units[i], field_units[j], entry_units[i][j]))
                       for j in range(5)] for i in range(5)]
    journal.value('unit-join-operands.pickle',
                  (field_units, row_units, entry_units, unit_residuals))
    require(all(v == 0 for row in unit_residuals for u in row for v in u),
            'consistent physical row-field unit differences')
    inverse_units = [[tuple(a-b for a, b in zip(field_units[i], row_units[j]))
                      for j in range(5)] for i in range(5)]

    # Entry symbols are temporary algebraic carriers, not output placeholders.
    generic = [[sp.Symbol('s11cdKernelEntry_%d_%d' % (i, j))
                for j in range(5)] for i in range(5)]
    binding = {generic[i][j]: source[i, j] for i in range(5) for j in range(5)}
    journal.value('entry-carriers.pickle', {'generic': generic, 'sourceBinding': binding})
    det_generic = journal.op('generic-determinant', determinant, generic)
    det_source = journal.op('source-determinant', lambda a, b: a.xreplace(b),
                            det_generic, binding)
    require(det_source != 0, 'reference determinant is not an identically zero expression')
    adj_generic = [[None]*5 for _ in range(5)]
    adj_source = [[None]*5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            adj_generic[i][j] = journal.op('cofactor-%d-%d' % (i, j),
                                           adjugate_entry, generic, i, j)
            adj_source[i][j] = journal.op('source-cofactor-%d-%d' % (i, j),
                lambda a, b: a.xreplace(b), adj_generic[i][j], binding)
    inverse = journal.op('actual-inverse', lambda a, d:
        sp.ImmutableMatrix([[sp.Mul(v, sp.Pow(d, -1, evaluate=False), evaluate=False)
                             for v in row] for row in a]), adj_source, det_source)
    require(not inverse.has(*binding), 'no temporary carriers in physical inverse')
    # Generic polynomial certificates lift through the explicitly saved exact
    # scalar-entry binding. This is algebraic recipe validation, not physics review.
    proofs = []
    for side in ('left', 'right'):
        for i in range(5):
            for j in range(5):
                raw = sum((generic[i][k]*adj_generic[k][j] if side == 'left'
                           else adj_generic[i][k]*generic[k][j]) for k in range(5))
                raw -= det_generic if i == j else 0
                residual = journal.op('%s-polynomial-%d-%d' % (side, i, j), sp.expand, raw)
                proofs.append(residual)
    journal.value('algebraic-certificate.pickle',
        {'genericIdentityResiduals': proofs, 'sourceBinding': binding,
         'identityDomain': sp.Ne(det_source, 0, evaluate=False),
         'inverseEntryUnits': inverse_units})
    require(all(v == 0 for v in proofs), 'left/right explicit cofactor identities')
    for side in ('left', 'right'):
        residual = sp.ImmutableMatrix(5, 5, lambda i, j: sp.Add(
            *(sp.Mul(source[i, k], inverse[k, j], evaluate=False) if side == 'left'
              else sp.Mul(inverse[i, k], source[k, j], evaluate=False) for k in range(5)),
            -int(i == j), evaluate=False))
        journal.value(side + '-physical-residual.pickle', residual)

    # Preserve supplied real-frequency branch functions and Fourier data.
    # Never infer a contour prescription from a pointwise inverse identity.
    fourier = state['fourier_mass']
    spectral_unit = tuple(known[kn])
    kernel_units = [[tuple(a+b for a, b in zip(unit, spectral_unit))
                     for unit in row] for row in inverse_units]
    journal.value('physical-units.pickle', {
        'fieldUnits': field_units, 'rowUnits': row_units,
        'symbolEntryUnits': entry_units, 'inverseEntryUnits': inverse_units,
        'spectralMeasureUnit': spectral_unit, 'kernelEntryUnits': kernel_units,
        'numericFrame': physical['unit_frame']})
    phase = sp.exp(sp.I*kn*(state['z']-state['zp']))
    density = journal.op('regular-spectral-density',
        lambda a, b, c: a.applyfunc(lambda v: sp.Mul(b, v, sp.Pow(c, -1), evaluate=False)),
        inverse, phase, fourier)
    branch_nodes = sorted(source.atoms(sp.Piecewise), key=sp.default_sort_key)
    radicals = sorted([v for v in source.atoms(sp.Pow) if v.exp.is_Rational
                       and v.exp.q != 1], key=sp.default_sort_key)
    journal.value('outgoing-source-context.pickle', {
        'timeCharacter': sp.exp(-sp.I*state['omega']*state['t']),
        'normalCharacter': phase, 'normalMomentum': kn,
        'normalizationOperands': state['normalization_operands'],
        'fourierMass': fourier, 'constantFourierMass': packet['constantFourierMass'],
        'sourceBranchEquations': state['branch_equations'],
        'sourceBranchMap': state['branch_map'], 'actualPiecewiseNodes': branch_nodes,
        'actualRadicalNodes': radicals, 'regularSpectralDensity': density,
        'profileRegulator': state['regulator'], 'referenceCoupling': packet['records']['REFERENCE']['coupling']})

    mp.mp.dps = 80
    parameters = physical['parameters']
    live = sorted(source.free_symbols - {kn}, key=sp.default_sort_key)
    require(all(s.name in parameters for s in live),
            ('all physical coefficient bindings supplied', [s.name for s in live if s.name not in parameters]))
    bindings = {s: sp.sympify(parameters[s.name]) for s in live}
    length = sp.sympify(parameters['L_W'])
    probe_records = []
    # Three fixed algebraic probes, not modes or new physical incident inputs.
    for number, k in enumerate((sp.S.Zero, 1/length, 2/length)):
        point = {**bindings, kn: k}
        numeric_source = journal.op('probe-source-%d' % number,
            lambda a, b: a.xreplace(b).evalf(80), source, point)
        numeric_inverse = journal.op('probe-inverse-%d' % number,
            lambda a, b: a.xreplace(b).evalf(80), inverse, point)
        def convert(matrix):
            return mp.matrix([[mp.mpc(str(sp.re(v)), str(sp.im(v))) for v in row]
                              for row in matrix.tolist()])
        def native(matrix):
            # Persist decimal values in the same standard SymPy codec, never
            # context-generated mpmath matrix classes that may not unpickle.
            return sp.ImmutableMatrix(matrix.rows, matrix.cols,
                lambda i, j: sp.Float(str(mp.re(matrix[i, j])), 80)
                             + sp.I*sp.Float(str(mp.im(matrix[i, j])), 80))
        a, r = convert(numeric_source), convert(numeric_inverse)
        def direct_solve(matrix):
            converted = convert(matrix)
            columns = [mp.lu_solve(converted, mp.eye(5)[:, j]) for j in range(5)]
            return native(mp.matrix([[columns[j][i, 0] for j in range(5)] for i in range(5)]))
        independent_value = journal.op('probe-independent-%d' % number,
                                       direct_solve, numeric_source)
        independent = convert(independent_value)
        condition = mp.mnorm(a, mp.inf)*mp.mnorm(independent, mp.inf)
        left, right = a*r-mp.eye(5), r*a-mp.eye(5)
        difference = mp.mnorm(r-independent, mp.inf)/(1+mp.mnorm(independent, mp.inf))
        entry = next(((i, j) for i in range(5) for j in range(5) if a[i, j] != 0), None)
        require(entry is not None, 'nonzero source entry for one-sided mutation')
        mutated = a.copy()
        mutated[entry[0], entry[1]] *= -1
        mutation = mutated*r-mp.eye(5)
        record = {'k': str(k), 'conditionEstimate': str(condition),
                  'leftResidualNorm': str(mp.mnorm(left, mp.inf)),
                  'rightResidualNorm': str(mp.mnorm(right, mp.inf)),
                  'independentRelativeDifference': str(difference),
                  'mutationEntry': entry, 'mutationResponseNorm': str(mp.mnorm(mutation-left, mp.inf))}
        journal.value('probes/%d/complete-operands.pickle' % number,
            {'bindings': point, 'source': numeric_source, 'inverse': numeric_inverse,
             'independent': independent_value, 'workingDecimalPrecision': 80,
             'leftResidual': native(left), 'rightResidual': native(right),
             'mutatedSource': native(mutated), 'mutationResidual': native(mutation),
             'record': record})
        save_json(out/'probes'/str(number)/'summary.json', record)
        probe_records.append(record)
        require(mp.isfinite(condition) and condition < mp.mpf('1e12'), 'regular bounded-condition algebraic probe')
        require(max(mp.mnorm(left, mp.inf), mp.mnorm(right, mp.inf), difference) < mp.mpf('1e-40'),
                'high-precision binding of actual inverse and independent solve')
        require(mp.mnorm(mutation-left, mp.inf) > mp.mpf('1e-20'), 'one-sided symbol mutation moves residual')
    # This construction deliberately does not prescribe how to cross real
    # determinant singularities or continue the complete coupled symbol.
    prescription = {
        'status': 'REGULAR_SPECTRAL_DENSITY_BUILT_OUTGOING_PRESCRIPTION_PENDING',
        'constructed': ['full saved reference symbol inverse', 'regular spectral density',
                        'saved bulk branch and Fourier dependency context'],
        'unconstructed': ['whole-line outgoing singular-integral prescription',
                          'profile/two-asymptote response', 'radiating-domain witness'],
        'reason': 'A pointwise inverse on det(P0) != 0 is not a prescription at real-axis singularities. Source bulk branch data are retained for the separate coupled outgoing-contour adjudication.',
        'formExport': False, 'newPhysicalBinding': False, 'radiatingCoverage': False,
        'centreEliminated': False, 'independentPhysicsClearance': False}
    save_json(out/'prescription-status.json', prescription)
    post = [route(r['path']) for r in actual]
    save_json(out/'posthashes.json', post)
    require(actual == post, 'all consumed routes/bytes unchanged')
    save_json(out/'operation-index.json', journal.records)
    return {'status': prescription['status'], 'case': 'LAB_HELD__RHO4_CONSTANT',
            'sourceSha256': digest(spec['uniformPacket']), 'probes': probe_records,
            'completedOperations': len(journal.records), 'posthashesUnchanged': True,
            'recipeCertificate': 'explicit cofactor identity plus exact source substitution',
            'scientificScope': 'New inverse-symbol ingredient only; no completed science replay',
            'outgoingGreenOperatorCompleted': False, 'formCompleted': False}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--input-manifest', type=Path, required=True)
    parser.add_argument('--gate-receipt', type=Path, required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    observed = containment()
    spec = json.loads(args.input_manifest.read_text())
    gate = json.loads(args.gate_receipt.read_text())
    require(gate['status'] == 'READY_FOR_ONE_GUARDED_REFERENCE_KERNEL_STAGE', 'implementation gate not closed')
    require(gate['workerSha256'] == digest(Path(__file__)), 'reviewed worker identity')
    require(gate['inputManifestSha256'] == digest(args.input_manifest), 'reviewed input identity')
    require(gate['scopeAuthorization'] == 'Approved', 'bounded stage authorization')
    out = args.run_directory.resolve()
    out.relative_to(STORE)
    out.mkdir(parents=True, exist_ok=False)
    save_json(out/'native-containment.json', observed)
    save_json(out/'input-manifest.json', spec)
    save_json(out/'gate-receipt.json', gate)
    started = time.monotonic()
    try:
        result = construct(spec, out)
    except BaseException as exc:
        save_json(out/'failure.json', {'type': type(exc).__name__, 'message': str(exc),
                  'traceback': traceback.format_exc(), 'wallSeconds': time.monotonic()-started,
                  'automaticRetry': False, 'completedArtifactsPreserved': True})
        raise
    result['wallSeconds'] = time.monotonic()-started
    save_json(out/'checks.json', result)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
