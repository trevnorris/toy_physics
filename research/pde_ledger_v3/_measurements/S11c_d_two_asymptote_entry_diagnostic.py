#!/usr/bin/env python3
"""Diagnose only the incomplete block17-regular-0-0 entry; no continuation.

Restores the failed operation's saved input. Recomputes only its unreturned
quadratic-branch/Taylor intermediates using unchanged reviewed helpers, saves
them before observing zero flags, and stops before quotient/end-lift work.
All execution requires separate explicit approval and the shared guard.
"""
import argparse
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
STORE = ROOT/'_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


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


def diagnose(spec, out, j):
    # main verified all pinned bytes and actual containment before this import.
    import sympy as sp


    def simplify(value):
        # Local scalar algebra only. Native integrals are never sent here.
        require(not value.has(sp.Integral, sp.Limit), 'algebraic reducer received an action integral')
        return sp.cancel(sp.simplify(value))


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


    def restore(path):
        with Path(path).open('rb') as stream:
            return SavedCodec(stream).load()

    started = json.loads(Path(spec['inputs']['incompleteStarted']['path']).read_text())
    require(started['name'] == 'block17-regular-0-0'
            and started['input'] == spec['inputs']['incompleteInput'],
            'incomplete input/operation receipt mismatch')
    parent = Path(spec['inputs']['incompleteInput']['path']).parent
    require(not (parent/'value.pickle').exists() and not (parent/'completed.json').exists(),
            'the parent entry now has a return; no replay permitted')
    operands = j.op('restore-incomplete-entry-input', restore,
                    spec['inputs']['incompleteInput']['path'])
    require(isinstance(operands, tuple) and len(operands) == 5, 'saved entry schema')
    chart, args, residue, k0, q0 = operands
    require(isinstance(args, tuple) and len(args) == 4, 'saved chart argument schema')
    value, kn, square, radical = args
    require(chart['backJoin'] == 0 and chart['rationalInCoordinates'],
            'saved accepted inverse chart identity')
    require(isinstance(kn, sp.Symbol) and isinstance(radical, sp.Symbol),
            'saved scalar coordinate schema')
    j.value('unchanged-entry-operands.pickle', operands)
    # The following is the failed entry's existing series calculation split into
    # durable operations. None of the parent's 61 complete operations is called.
    r0 = j.op('inherited-positive-radical', lambda q: simplify(q/sp.I), q0)
    branch_join = j.op('inherited-radical-square-join',
                      lambda r, s, k, point: simplify(r*r-s.subs(k, point)),
                      r0, square, kn, k0)
    require(r0.is_positive is True and branch_join == 0,
            'inherited nonzero positive radical at pole')
    degree = 3
    radpoly = j.op('quadratic-branch-polynomial', lambda s, k: sp.Poly(s, k), square, kn)
    rhs = j.op('quadratic-branch-taylor',
        lambda s, k, point, d: [simplify(sp.diff(s, k, n).subs(k, point)/sp.factorial(n))
                               for n in range(d+1)], square, kn, k0, degree)
    require(radpoly.degree() == 2, 'bounded quadratic branch')
    rseries = [r0]
    for n in range(1, degree+1):
        coefficient = j.op('radical-series-coefficient-'+str(n),
            lambda rhs, prior, r, order:
                simplify((rhs[order]-sum(prior[t]*prior[order-t]
                         for t in range(1, order)))/(2*r)), rhs, rseries, r0, n)
        rseries.append(coefficient)
    j.value('branch-series.pickle', rseries)
    numerator, denominator = j.op('saved-rational-fraction', sp.fraction, chart['rational'])
    nseries = j.op('numerator-series', polynomial_series,
                   numerator, kn, radical, k0, rseries, degree)
    dseries = j.op('denominator-series', polynomial_series,
                   denominator, kn, radical, k0, rseries, degree)
    # Full series returns are on disk before asking SymPy for zero status.
    def observe(series):
        records = []
        for index, coefficient in enumerate(series):
            structural = coefficient != 0
            flag = coefficient.is_zero
            records.append({'index': index, 'coefficient': coefficient,
                'expressionText': sp.sstr(coefficient), 'expressionSrepr': sp.srepr(coefficient),
                'structurallyUnequalZero': structural, 'isZero': flag,
                'isZeroType': type(flag).__module__+'.'+type(flag).__qualname__,
                'isZeroIsPythonFalse': flag is False,
                'isZeroIsPythonTrue': flag is True, 'isZeroIsNone': flag is None})
        order = next((i for i, v in enumerate(series) if v != 0), None)
        order_supported = order is not None and order <= 2
        guard_passes = order_supported and series[order].is_zero is False
        return {'coefficients': records, 'selectedOrder': order,
                'originalOrderSupported': order_supported,
                'originalNonzeroGuardWouldPass': guard_passes}
    flags = j.op('observe-original-denominator-guard', observe, dseries)
    summary = {k: v for k, v in flags.items() if k != 'coefficients'}
    summary['coefficients'] = [{k: v for k, v in r.items()
                              if k not in ('coefficient', 'isZero', 'structurallyUnequalZero')}
                              | {'structurallyUnequalZero': bool(r['structurallyUnequalZero']),
                                 'isZeroRepr': repr(r['isZero'])}
                              for r in flags['coefficients']]
    summary.update(status='INCOMPLETE_ENTRY_SERIES_DIAGNOSTIC_ONLY',
                   sourceInput=spec['inputs']['incompleteInput'],
                   originalFailureReproduced=not flags['originalNonzeroGuardWouldPass'],
                   numeratorSeriesLength=len(nseries), denominatorSeriesLength=len(dseries),
                   quotientComputed=False, regularCoefficientAccepted=False,
                   diagnosticIsNewZeroProof=False, parentCompletedOperationsReplayed=0)
    save(out/'denominator-observation.json', summary)
    return {'status': 'INCOMPLETE_ENTRY_SERIES_DIAGNOSTIC_ONLY',
            'incompleteOperation': 'block17-regular-0-0',
            'sourceInputSha256': spec['inputs']['incompleteInput']['sha256'],
            'parentCompletedOperationsReplayed': 0, 'completedOperations': len(j.records),
            'originalGuardWouldPass': flags['originalNonzeroGuardWouldPass'],
            'selectedOrder': flags['selectedOrder'], 'quotientComputed': False,
            'endLiftConstructed': False, 'forcingConstructed': False,
            'newRoots': 0, 'newModeSolves': 0, 'newLUSolves': 0, 'producerCalls': 0,
            'integralsEvaluated': False, 'numericalProbes': 0,
            'scientificStageAccepted': False, 'freshIndependentClearClaimed': False,
            'automaticRetry': False}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--input-manifest', type=Path, required=True)
    parser.add_argument('--gate-receipt', type=Path, required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    spec = json.loads(args.input_manifest.read_text())
    gate = json.loads(args.gate_receipt.read_text())
    require(gate['status'] == 'READY_FOR_ONE_GUARDED_TWO_ASYMPTOTE_ENTRY_DIAGNOSTIC',
            'diagnostic approval gate incomplete')
    require(gate['workerSha256'] == digest(Path(__file__))
            and gate['inputManifestSha256'] == digest(args.input_manifest), 'pinned diagnostic identity')
    require(gate['scopeExplicitlyApproved'] is True and gate['automaticRetry'] is False,
            'explicit diagnostic approval absent')
    require(gate['seconds'] == 900 and gate['nativeSeconds'] == 840,
            'diagnostic duration contract')
    require(gate['sharedGuardSha256'] == digest(ROOT/'scripts/s11c_guarded_run.py')
            and gate['supervisorSha256'] == digest(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),
            'guard/supervisor identity')
    approval_path = Path(gate['scopeApprovalPath'])
    require(digest(approval_path) == gate['scopeApprovalSha256'], 'approval identity')
    approval = json.loads(approval_path.read_text())
    require(approval['status'] == 'EXPLICITLY_APPROVED_ONE_GUARDED_TWO_ASYMPTOTE_ENTRY_DIAGNOSTIC'
            and approval['workerSha256'] == digest(Path(__file__))
            and approval['inputManifestSha256'] == digest(args.input_manifest)
            and approval['pendingGateSha256'] == gate['pendingGateSha256'],
            'approval scope mismatch')
    actual = {name: route(r['path']) for name, r in spec['inputs'].items()}
    require(actual == spec['inputs'], 'pinned parent/input bytes changed; no science imported')
    observed = containment()
    out = args.run_directory.resolve(); out.relative_to(STORE)
    require(str(out) == gate['resultDirectory'], 'wrong diagnostic result directory')
    out.mkdir(parents=True, exist_ok=False)
    save(out/'native-containment.json', observed)
    save(out/'input-manifest.json', spec); save(out/'gate-receipt.json', gate)
    save(out/'prehashes.json', actual)
    journal = Journal(out); started = time.monotonic()
    try:
        result = diagnose(spec, out, journal)
        signal.alarm(0)  # diagnostic completed; only hash/index bookkeeping remains
        posthashes = {name: route(v['path']) for name, v in spec['inputs'].items()}
        save(out/'posthashes.json', posthashes)
        require(posthashes == spec['inputs'], 'input posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started, allSourceHashesUnchanged=True)
        save(out/'operation-index.json', journal.records); save(out/'checks.json', result)
        save(out/'artifact-index.json', {str(p.relative_to(out)): route(p)
             for p in sorted(out.rglob('*')) if p.is_file()})
        print(json.dumps(result, indent=2, allow_nan=False))
    except BaseException:
        failure_traceback = traceback.format_exc()
        signal.alarm(0)
        post = {}
        for name, record in spec['inputs'].items():
            try: post[name] = route(record['path'])
            except OSError as error: post[name] = {'error': str(error), 'path': record['path']}
        save(out/'failure-posthashes.json', post)
        save(out/'partial-operation-index.json', journal.records)
        save(out/'failure.json', {'traceback': failure_traceback,
             'completedOperations': len(journal.records), 'incompleteOperation': journal.active,
             'wallSeconds': time.monotonic()-started, 'automaticRetry': False})
        save(out/'failure-artifact-index.json', {str(p.relative_to(out)): route(p)
             for p in sorted(out.rglob('*')) if p.is_file()})
        raise


if __name__ == '__main__':
    main()
