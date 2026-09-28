#!/usr/bin/env python3
"""One saved-pencil feasibility diagnostic, not a scattering/flux constructor."""
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
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def route(path):
    path = Path(path)
    return dict(path=str(path), canonicalPath=str(path.resolve(strict=True)),
                bytes=path.stat().st_size, sha256=sha(path))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.')) or
                module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError((module, name))
        return super().find_class(module, name)


class Journal:
    def __init__(self, out):
        self.out, self.records, self.active = out, [], None

    def op(self, name, function, *args):
        folder = self.out/'operations'/('%03d-%s' % (len(self.records), name))
        folder.mkdir(parents=True, exist_ok=False)
        for filename, value in [('input.pickle', args)]:
            with (folder/filename).open('xb') as stream:
                pickle.dump(value, stream, protocol=4)
                stream.flush(); os.fsync(stream.fileno())
        self.active = dict(name=name, input=route(folder/'input.pickle'),
                           startedUtc=datetime.now(timezone.utc).isoformat())
        save(folder/'started.json', self.active)
        result = function(*args)
        with (folder/'value.pickle').open('xb') as stream:
            pickle.dump(result, stream, protocol=4)
            stream.flush(); os.fsync(stream.fileno())
        self.active.update(value=route(folder/'value.pickle'),
                           finishedUtc=datetime.now(timezone.utc).isoformat())
        save(folder/'completed.json', self.active)
        self.records.append(self.active); self.active = None
        return result


def containment():
    group = next(v[3:] for v in Path('/proc/self/cgroup').read_text().splitlines()
                 if v.startswith('0::'))
    cgroup = Path('/sys/fs/cgroup')/group.lstrip('/')
    limits = {k:(cgroup/k).read_text().strip()
              for k in ('memory.max', 'memory.swap.max', 'pids.max')}
    limits.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={k:os.environ.get(k) for k in THREADS})
    require(limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0'
            and limits['pids.max']=='32' and limits['nice']>=15
            and len(limits['affinity'])==1 and all(v=='1' for v in limits['threads'].values()),
            'required whole-job containment absent; no science imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise TimeoutError('840-second native stop; no retry')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(840)
    return {**limits, 'nativeSeconds':840}


def construct(spec, out, j):
    import sympy as sp
    import mpmath as mp
    mp.mp.dps = 70

    def display(value):
        if isinstance(value, dict):return {str(k):display(v) for k,v in value.items()}
        if isinstance(value, (tuple, list)):return [display(v) for v in value]
        if isinstance(value, sp.MatrixBase):return [[str(value[i,h]) for h in range(value.cols)] for i in range(value.rows)]
        if isinstance(value, sp.Basic):return str(value)
        return value

    def op(name, function, *args):
        value = j.op(name, function, *args)
        save(out/(name+'.json'), display(value))
        return value

    def restore(record):
        require(route(record['path'])==record, 'saved input identity')
        with Path(record['path']).open('rb') as stream:return SavedCodec(stream).load()

    def bind(packet, target):
        k=packet['momentum'];q=packet['radical'];w=packet['frequency']
        matrix=packet['livePencil'].subs(w,target)
        wave=sp.expand(packet['wave'].subs(w,target))
        polynomial=sp.Poly(wave,q)
        supported=(polynomial.degree()==2 and polynomial.nth(1)==0 and polynomial.nth(2)!=0)
        square=sp.cancel(-polynomial.nth(0)/polynomial.nth(2)) if supported else None
        return {'matrix':matrix,'wave':wave,'momentum':k,'radical':q,
                'targetFrequency':target,'radicalSquare':square,
                'originalFrequency':w,'originalOrigin':packet['origin'],
                'savedUnitConvention':packet['unitConvention'],
                'checks':{'fiveFields':matrix.shape==(5,5),
                          'onlyBoundMomentumAndRadical':not (matrix.free_symbols|wave.free_symbols)-{k,q},
                          'realMomentumCoordinate':k.is_real is True,
                          'quadraticRadicalRelation':supported,
                          'savedBranchJoinsZero':all(v==0 for v in packet['branchResiduals']),
                          'savedSeedJoinsZero':all(v==0 for v in packet['referencePencilResidual'])
                                               and packet['referenceWaveResidual']==0}}

    def determinant(bound):
        matrix=bound['matrix'];denominators=[];rows=[]
        for i in range(5):
            entries=[sp.cancel(matrix[i,h]) for h in range(5)]
            denominator=sp.lcm([sp.denom(v) for v in entries])
            denominators.append(denominator)
            rows.append([sp.cancel(v*denominator) for v in entries])
        cleared=sp.ImmutableMatrix(rows)
        numerator, denominator=sp.fraction(sp.cancel(cleared.det(method='domain-ge')/sp.prod(denominators)))
        return {'originalMatrix':matrix,'clearedMatrix':cleared,'rowDenominators':denominators,
                'numerator':numerator,'denominator':denominator,
                'checks':{'nonzeroDeterminantNumerator':numerator!=0}}

    def eliminate(det, bound):
        k=bound['momentum'];q=bound['radical']
        result=sp.expand(sp.resultant(det['numerator'],bound['wave'],q))
        real,imag=result.as_real_imag()
        re_poly=sp.Poly(real,k,domain=sp.QQ);im_poly=sp.Poly(imag,k,domain=sp.QQ)
        common=sp.gcd(re_poly,im_poly)
        return {'resultant':result,'realPart':re_poly.as_expr(),'imaginaryPart':im_poly.as_expr(),
                'commonPolynomial':common.as_expr(),'degree':int(common.degree()) if not common.is_zero else None,
                'realDivisionRemainder':sp.rem(re_poly,common).as_expr() if not common.is_zero else None,
                'imagDivisionRemainder':sp.rem(im_poly,common).as_expr() if not common.is_zero else None,
                'checks':{'nonzeroResultant':result!=0,'boundedDegree':not common.is_zero and common.degree()<=64}}

    def isolate(elimination, bound):
        k=bound['momentum'];poly=sp.Poly(elimination['commonPolynomial'],k,domain=sp.QQ).sqf_part()
        intervals=poly.intervals(eps=sp.Rational(1,10**45)) if poly.degree()>0 else []
        return {'squareFreePolynomial':poly.as_expr(),'intervals':intervals,
                'scope':'All real zeros of the necessary elimination polynomial at this ONE frequency; not physical modes.',
                'checks':{'realRootCountMatches':sum(m for _,m in intervals)==poly.count_roots(-sp.oo,sp.oo),
                          'boundedCandidateCount':len(intervals)<=32}}

    def number(value):
        return mp.mpc(str(sp.re(value).evalf(75)),str(sp.im(value).evalf(75)))

    def numeric_matrix(matrix, substitutions):
        return mp.matrix([[number(matrix[i,h].subs(substitutions)) for h in range(matrix.cols)]
                          for i in range(matrix.rows)])

    def encoded(matrix):
        return [[{'real':mp.nstr(mp.re(matrix[i,h]),75),'imag':mp.nstr(mp.im(matrix[i,h]),75)}
                 for h in range(matrix.cols)] for i in range(matrix.rows)]

    def inspect_candidate(interval, bound, det):
        (lo,hi),multiplicity=interval;k0=(lo+hi)/2
        q2=sp.cancel(bound['radicalSquare'].subs(bound['momentum'],k0))
        record={'interval':interval,'midpoint':k0,'radicalSquare':q2,
                'currentComputed':False,'physicalWitnessEstablished':False}
        if q2.is_real is not True or q2==0:
            return {**record,'status':'UNRESOLVED_BRANCH_OR_THRESHOLD'}
        q0=sp.sqrt(q2)  # Positive-frequency principal lift from supplied real-axis germ.
        sub={bound['momentum']:k0,bound['radical']:q0}
        denominator=det['denominator'].subs(sub)
        row_denominators=[v.subs(sub) for v in det['rowDenominators']]
        record.update(radical=q0,denominator=denominator,rowDenominators=row_denominators,
                      bulkAtThisEndMode='propagating' if q2.is_positive is True else 'evanescent')
        if abs(number(denominator))<mp.mpf('1e-25') or any(abs(number(v))<mp.mpf('1e-25') for v in row_denominators):
            return {**record,'status':'UNRESOLVED_DENOMINATOR_NEAR_ZERO'}
        raw=numeric_matrix(bound['matrix'],sub)
        scales=[max(mp.mpf(1),mp.sqrt(sum(abs(raw[i,h])**2 for h in range(5)))) for i in range(5)]
        scaled=mp.matrix([[raw[i,h]/scales[i] for h in range(5)] for i in range(5)])
        _,singular,v=mp.svd(scaled)
        null=[i for i in range(5) if singular[i]<mp.mpf('1e-25')]
        ambiguous=any(mp.mpf('1e-25')<=s<mp.mpf('1e-12') for s in singular)
        right=v.H
        basis=mp.matrix([[right[i,h] for h in null] for i in range(5)]) if null else None
        record.update(matrix=encoded(raw),rowScales=[mp.nstr(x,75) for x in scales],
                      singularValues=[mp.nstr(x,75) for x in singular],nullity=len(null),
                      status='AMBIGUOUS_RANK' if ambiguous else 'CANDIDATE_NULLSPACE' if null else 'REJECTED_ON_SELECTED_LIFT')
        if not null:return record
        residual=scaled*basis
        thickness=mp.sqrt(sum(abs(basis[i,h])**2 for i in (3,4) for h in range(len(null))))
        record.update(rightBasis=encoded(basis),scaledResidual=encoded(residual),
                      thicknessCoordinateNorm=mp.nstr(thickness,75),
                      thicknessClassification='NONZERO_AT_DIAGNOSTIC_PRECISION' if thickness>mp.mpf('1e-12') else 'UNRESOLVED_TINY' if thickness>mp.mpf('1e-25') else 'NO_RESOLVED_THICKNESS_COMPONENT',
                      classificationScope='Entire numerical nullspace in inherited field reference units; not a flux selector or absence theorem.')
        # Addressed one-sided row/column corruption against the actual full basis.
        scores=[sum(abs(basis[h,c])**2 for c in range(len(null))) for h in range(5)]
        column=max(range(5),key=lambda h:scores[h]);mutated=scaled.copy();mutated[0,column]+=1
        changed=mutated*basis
        movement=mp.sqrt(sum(abs(changed[i,h]-residual[i,h])**2 for i in range(5) for h in range(len(null))))
        record.update(mutation={'row':0,'column':column,'addedReferenceCoefficient':1,
                                'baselineMatrix':encoded(scaled),'mutatedMatrix':encoded(mutated),
                                'baselineResidual':encoded(residual),'mutatedResidual':encoded(changed),
                                'movement':mp.nstr(movement,75),'responsive':bool(movement>mp.mpf('1e-12'))})
        return record

    target=sp.Rational(spec['candidate']['omega'])
    unit_packet=op('restore-unit-map',restore,spec['inputs']['unitMap'])
    unit_map=unit_packet['coordinateMap']
    save(out/'physical-unit-map.json',display(unit_map))
    require(len(unit_map['fieldReferenceUnits'])==len(unit_map['equationReferenceUnits'])==5,
            'actual five-field/equation reference-unit map')
    summaries={}
    for label in ('REFERENCE','LEFT','RIGHT'):
        packet=op(label+'-restore',restore,spec['inputs'][label])
        bound=op(label+'-bind',bind,packet,target)
        require(all(bound['checks'].values()),'saved pencil binding/schema')
        det=op(label+'-determinant',determinant,bound)
        require(all(det['checks'].values()),'nonzero pencil determinant')
        elimination=op(label+'-elimination',eliminate,det,bound)
        require(all(elimination['checks'].values()),'supported finite necessary-root problem')
        require(elimination['realDivisionRemainder']==elimination['imagDivisionRemainder']==0,'exact common-factor reconstruction')
        roots=op(label+'-real-root-isolation',isolate,elimination,bound)
        require(all(roots['checks'].values()),'bounded real candidate census')
        candidates=[]
        for i,interval in enumerate(roots['intervals']):
            value=op(label+'-candidate-%02d'%i,inspect_candidate,interval,bound,det)
            if 'mutation' in value:require(value['mutation']['responsive'],'actual candidate basis must detect corruption')
            candidates.append(value)
        summaries[label]={'radicalSquare':str(bound['radicalSquare']),
                          'realCandidates':len(candidates),
                          'candidateStatuses':[v['status'] for v in candidates],
                          'thicknessClassifications':[v.get('thicknessClassification') for v in candidates],
                          'physicalIncidentCurrentEstablished':False,'physicalOutgoingFluxEstablished':False}
    return {'status':'ONE_INPUT_NECESSARY_CHANNEL_DIAGNOSTIC_COMPLETE_NOT_A11_A12_CLEARANCE',
            'candidate':spec['candidate'],'physicalUnitMap':display(unit_map),'ends':summaries,
            'currentNormalization':'UNCOMPUTED: no derivative or nullspace norm is relabelled as physical current.',
            'next':'Adjudicate necessary-channel evidence with c1/c2 source/domain review; stop with a cost/scope decision.',
            'responseOrBulkFluxConstructed':False,'frequencyScan':False,'automaticRetry':False}


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--input-manifest',required=True,type=Path)
    parser.add_argument('--gate-receipt',required=True,type=Path)
    parser.add_argument('--run-directory',required=True,type=Path)
    args=parser.parse_args();spec=json.loads(args.input_manifest.read_text());gate=json.loads(args.gate_receipt.read_text())
    require(gate['status']=='READY_FOR_ONE_GUARDED_RADIATING_FEASIBILITY_DIAGNOSTIC','independent method/instrument gate incomplete')
    require(gate['workerSha256']==sha(Path(__file__)) and gate['inputManifestSha256']==sha(args.input_manifest),'worker/input pins')
    require(gate['goalApproved'] is True and gate['repairedBuildReviewOutstanding'] is False,'scope and independent clearance')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False,'one ordinary bounded job')
    require(gate['sharedGuardSha256']==sha(ROOT/'scripts/s11c_guarded_run.py'),'unchanged guard')
    require(gate['supervisorSha256']==sha(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'unchanged supervisor')
    for record in spec['inputs'].values():require(route(record['path'])==record,'input identity')
    observed=containment();out=args.run_directory.resolve();out.relative_to(ROOT/'_scratch/s11c')
    out.mkdir(parents=True,exist_ok=False);j=Journal(out);started=time.monotonic()
    save(out/'native-containment.json',observed);save(out/'input-manifest.json',spec);save(out/'gate.json',gate)
    try:
        result=construct(spec,out,j)
        post={k:route(v['path']) for k,v in spec['inputs'].items()}
        save(out/'posthashes.json',post);require(post==spec['inputs'],'source posthashes')
        result.update(wallSeconds=time.monotonic()-started,completedOperations=len(j.records),allSourcesUnchanged=True)
        save(out/'operation-index.json',j.records);save(out/'checks.json',result)
        print(json.dumps(result,indent=2,allow_nan=False))
    except BaseException:
        signal.alarm(0)
        save(out/'failure.json',dict(traceback=traceback.format_exc(),incompleteOperation=j.active,
                                    completedOperations=len(j.records),wallSeconds=time.monotonic()-started,automaticRetry=False))
        save(out/'partial-operation-index.json',j.records)
        save(out/'failure-posthashes.json',{k:route(v['path']) for k,v in spec['inputs'].items()})
        raise
    finally:
        save(out/'artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})


if __name__=='__main__':main()
