#!/usr/bin/env python3
"""Validate only saved reference-kernel artifacts; never import a producer."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import traceback

THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def require(value, message):
    if not value:
        raise ValueError(message)


def containment():
    group = next(v[3:] for v in Path('/proc/self/cgroup').read_text().splitlines()
                 if v.startswith('0::'))
    root = Path('/sys/fs/cgroup')/group.lstrip('/')
    result = {n: (root/n).read_text().strip()
              for n in ('memory.max', 'memory.swap.max', 'pids.max')}
    result.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
        affinity=sorted(os.sched_getaffinity(0)), threads={n: os.environ.get(n) for n in THREADS})
    require(result['memory.max'] == str(2*1024**3) and result['memory.swap.max'] == '0'
        and result['pids.max'] == '32' and result['nice'] >= 15
        and len(result['affinity']) == 1 and all(v == '1' for v in result['threads'].values()),
        'required validator containment missing')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise TimeoutError('native validator limit')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(900)
    return result


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.'))
                or module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError((module, name))
        return super().find_class(module, name)


def validate(source, manifest, out):
    import sympy as sp
    import mpmath as mp
    mp.mp.dps = 80
    pins = json.loads(manifest.read_text())
    require(pins['sourceRoot'] == str(source), 'exact saved artifact root')
    require(pins['validatorSha256'] == sha(Path(__file__)), 'validator source pin')
    actual = {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*')) if p.is_file()}
    save(out/'prehashes.json', actual)
    require(actual == pins['files'], 'all saved output files match inventory')
    save(out/'input-manifest.json', pins)
    def load(name):
        require(name in actual, 'unmanifested saved operand: '+name)
        with (source/name).open('rb') as f:
            return SavedCodec(f).load()
    operations = json.loads((source/'operation-index.json').read_text())
    named = {r['name']: r for r in operations}
    def value(name):
        p = Path(named[name]['value']['path'])
        return load(str(p.relative_to(source)))
    selected = load('selected-inputs.pickle')
    symbol = selected['reference']['strong']
    carriers = load('entry-carriers.pickle')
    certificate = load('algebraic-certificate.pickle')
    inverse = value('actual-inverse')
    det = value('source-determinant')
    context = load('outgoing-source-context.pickle')
    units = load('physical-units.pickle')
    fields, rows, entries, residuals = load('unit-join-operands.pickle')
    checks = {
        'symbolShape': symbol.shape == (5, 5),
        'inverseShape': inverse.shape == (5, 5),
        'exactSourceBinding': all(carriers['sourceBinding'][carriers['generic'][i][j]] == symbol[i,j]
            for i in range(5) for j in range(5)),
        'certificateSourceBinding': certificate['sourceBinding'] == carriers['sourceBinding'],
        'fiftySavedZeroCertificates': len(certificate['genericIdentityResiduals']) == 50
            and all(v == 0 for v in certificate['genericIdentityResiduals']),
        'determinantDomain': certificate['identityDomain'] == sp.Ne(det, 0, evaluate=False),
        'inverseScalarFactors': all(inverse[i,j] == sp.Mul(value('source-cofactor-%d-%d' % (i,j)),
            sp.Pow(det,-1)) for i in range(5) for j in range(5)),
        'noUnboundCarriers': not inverse.has(*carriers['sourceBinding']),
        'savedUnitResidualsZero': all(v == 0 for row in residuals for entry in row for v in entry),
        'actualUnitJoins': all(tuple(a-b for a,b in zip(rows[i],fields[j])) == entries[i][j]
            for i in range(5) for j in range(5)),
        'inverseUnitOrientation': all(tuple(a-b for a,b in zip(fields[i],rows[j])) == units['inverseEntryUnits'][i][j]
            for i in range(5) for j in range(5)),
        'savedDensityContext': context['regularSpectralDensity'] == value('regular-spectral-density'),
        'densityScalarFactors': all(context['regularSpectralDensity'][i,j] == sp.Mul(
            context['normalCharacter'], inverse[i,j], sp.Pow(context['fourierMass'],-1))
            for i in range(5) for j in range(5)),
        'savedFourierMass': context['fourierMass'] == selected['reductionState']['fourier_mass'],
        'distinctProfileOperandRetained': context['constantFourierMass'] == selected['constantFourierMass']}
    for side in ('left', 'right'):
        residual = load(side+'-physical-residual.pickle')
        expected = sp.ImmutableMatrix(5,5,lambda i,j: sp.Add(
            *(sp.Mul(symbol[i,k],inverse[k,j]) if side == 'left'
              else sp.Mul(inverse[i,k],symbol[k,j]) for k in range(5)),
            -int(i==j)))
        checks[side+'SavedResidualExpression'] = residual == expected
    save(out/'structural-checks.json', checks)
    require(all(checks.values()), 'saved structural checks failed')
    numeric = []
    def convert(m):
        return mp.matrix([[mp.mpc(str(sp.re(v)),str(sp.im(v))) for v in row] for row in m.tolist()])
    def native(m):
        return sp.ImmutableMatrix(m.rows,m.cols,lambda i,j:
            sp.Float(str(mp.re(m[i,j])),80)+sp.I*sp.Float(str(mp.im(m[i,j])),80))
    for i in range(3):
        saved = load('probes/%d/complete-operands.pickle' % i)
        require(saved['source'] == value('probe-source-%d' % i)
            and saved['inverse'] == value('probe-inverse-%d' % i)
            and saved['independent'] == value('probe-independent-%d' % i), 'numeric operand joins')
        a,r,d,mut = (convert(saved[k]) for k in ('source','inverse','independent','mutatedSource'))
        left,right,mutation = a*r-mp.eye(5),r*a-mp.eye(5),mut*r-mp.eye(5)
        condition = mp.mnorm(a,mp.inf)*mp.mnorm(d,mp.inf)
        difference = mp.mnorm(r-d,mp.inf)/(1+mp.mnorm(d,mp.inf))
        movement = mp.mnorm(mutation-left,mp.inf)
        joins = max(mp.mnorm(left-convert(saved['leftResidual']),mp.inf),
            mp.mnorm(right-convert(saved['rightResidual']),mp.inf),
            mp.mnorm(mutation-convert(saved['mutationResidual']),mp.inf))
        entry = tuple(saved['record']['mutationEntry'])
        mutation_matches = all(mut[x,y] == (-a[x,y] if (x,y)==entry else a[x,y])
            for x in range(5) for y in range(5))
        result = {'probe': i, 'conditionEstimate': str(condition),
            'leftResidualNorm': str(mp.mnorm(left,mp.inf)), 'rightResidualNorm': str(mp.mnorm(right,mp.inf)),
            'independentRelativeDifference': str(difference), 'mutationResponseNorm': str(movement),
            'savedResidualRoundTripDifference': str(joins), 'oneSidedMutationMatches': mutation_matches}
        with (out/('probe-%d-recomputed.pickle' % i)).open('xb') as f:
            pickle.dump({'sourceArtifact': 'probes/%d/complete-operands.pickle' % i,
                'leftResidual':native(left),'rightResidual':native(right),'mutationResidual':native(mutation),
                'result':result},f,protocol=4)
            f.flush();os.fsync(f.fileno())
        save(out/('probe-%d.json' % i),result)
        require(condition < mp.mpf('1e12') and mp.isfinite(condition), 'condition bound')
        require(max(mp.mnorm(left,mp.inf),mp.mnorm(right,mp.inf),difference) < mp.mpf('1e-40'), 'inverse residuals')
        require(joins < mp.mpf('1e-70'), 'stored residual roundtrip')
        require(movement > mp.mpf('1e-20') and mutation_matches, 'mutation sensitivity')
        numeric.append(result)
    post = {str(p.relative_to(source)): sha(p) for p in sorted(source.rglob('*')) if p.is_file()}
    save(out/'posthashes.json', post)
    require(post == actual, 'saved results unchanged')
    result = {'status':'SAVED_REFERENCE_INGREDIENT_VALIDATED','structuralChecks':checks,
        'probes':numeric,'sourceRoot':str(source),'savedFiles':len(actual),
        'allSourceHashesUnchanged':True,'newConstructors':0,'newLUSolves':0,
        'outgoingGreenOperatorCompleted':False,'formCompleted':False,
        'finishedUtc':datetime.now(timezone.utc).isoformat()}
    save(out/'checks.json',result)
    print(json.dumps(result,indent=2))


def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--run-directory',type=Path,required=True)
    args=p.parse_args()
    observed=containment()
    root=Path('/var/projects/toy_physics/_scratch/s11c')
    source=args.source.resolve(); source.relative_to(root)
    out=args.run_directory.resolve(); out.relative_to(root)
    require(out != source and source not in out.parents, 'separate fresh validation output')
    out.mkdir(parents=True,exist_ok=False)
    save(out/'native-containment.json',observed)
    try:
        validate(source,args.manifest,out)
    except BaseException:
        save(out/'failure.json',{'traceback':traceback.format_exc(),'automaticRetry':False})
        raise


if __name__ == '__main__':
    main()
