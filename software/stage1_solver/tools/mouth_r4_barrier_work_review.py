"""Independent reduction of all executed signed banks for the held R4 barrier transition."""
import argparse
import hashlib
import json
from math import fsum, isfinite
from pathlib import Path


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def guard(path, expected=0):
    out=read(path/'outcome.json')
    limits=read(path/'limit-validation.json')
    if not (out['exitCode']==out['childOutcome']['exitCode']==expected and not out['childOutcome']['guardReason']
            and out['limitsVerified'] and limits['verified'] and int(limits['actual']['memory.swap.max'])==0):
        raise ValueError(f'Unverified scientific output at {path}')


def measure(scratch, held_sha, evidence_only=False):
    source=Path(__file__).resolve().parent
    if digest(source/'mouth_r4_barrier_held.py')!=held_sha:
        raise ValueError('Held source changed since its independent CLI review')
    work=scratch/'r4_barrier_work01'
    gate=read(work/'GATE.json')
    guard(work/'collection_guard', 0 if gate['passed'] else 1)
    if gate['absolute_tolerance']!=1.e-7:
        raise ValueError('Unexpected absolute derivative tolerance')
    failed=[]
    for name,sha in gate['sources'].items():
        if digest(source/name)!=sha:
            raise ValueError(f'Source differs from derivative evaluation: {name}')
    modes=['log_depth','scaled_log_exit','mean_width','depth_asymmetry','waist_profile','polar_equatorial']
    rows={}
    for case,eps in [('primary',1.e-8),('secondary',1.e-10)]:
        shape=scratch/'r4_eulerian_stage2_01'/case/'label_stall.npz'
        repair=scratch/'r4_slaved_melt_barrier01'/case/'fixed_active.npz'
        ck=gate['cases'][case]
        if ck['source_sha256']!=digest(shape) or ck['fixed_trace_sha256']!=digest(repair):
            raise ValueError(f'Unmatched physical baseline: {case}')
        max_err=max_act=max_unit=max_weighted=max_res=max_move=0.
        min_chart=1.e300
        detailed={}
        for name in ['active',*modes]:
            directory=work/(case+'_'+name)
            guard(directory.with_name(directory.name+'_guard'))
            doc=read(directory/'readout.json')
            if (doc!=ck['modes'][name] or doc['epsilon']!=eps or doc['test']!=name
                    or doc['case']!=case or doc['fixed_trace_sha256']!=ck['fixed_trace_sha256']
                    or doc['source_sha256']!=ck['source_sha256']):
                raise ValueError(f'Wrong linked bank {directory}')
            analytic=doc['work']['analytic']
            err=[]; gu=[]
            h0=2.e-6 if name=='active' else 1.e-5 if name=='log_depth' else 1.e-4
            probe=doc['work']['steps']
            if len(probe)!=2 or [x['h'] for x in probe]!=[h0,.5*h0]:
                raise ValueError(f'Unexpected signed widths for {directory}')
            for s in probe:
                banks=s['banks']
                if set(banks)!={'-1','1'} or not all(b['valid'] for b in banks.values()):
                    raise ValueError(f'Missing feasible sign for {directory}')
                num=.5*fsum(b['one_sided'] for b in banks.values())
                delta=num-analytic
                if not isfinite(delta) or abs(delta-s['central_minus_analytic'])>1.e-13:
                    raise ValueError(f'Could not reproduce central derivative at {directory}: {delta}')
                if abs(delta)>1.e-7:
                    failed.append([case,name,s['h'],'central check'])
                if name=='active':
                    max_act=max(max_act,abs(delta))
                else:
                    max_err=max(max_err,abs(delta))
                    gd=.5*fsum(b['one_sided_unit_gauge_force'] for b in banks.values())-doc['analytic_unit_gauge']
                    if not isfinite(gd):
                        raise ValueError(f'Nonfinite gauge-only derivative: {directory}')
                    max_unit=max(max_unit,abs(gd))
                    max_weighted=max(max_weighted,eps*abs(gd))
                    gu.append(gd)
                err.append(delta)
                for b in banks.values():
                    if not (b['min_chart']>0 and b['inner_residual']<=3.e-8 and b['iterations']>=0):
                        raise ValueError(f'Infeasible or nonstationary signed endpoint: {directory}')
                    max_res=max(max_res,b['inner_residual'])
                    min_chart=min(min_chart,b['min_chart'])
                    max_move=max(max_move,b['melted_label_change'])
            detailed[name]=dict(signed_central_minus_analytic=err,unit_gauge_central_minus_analytic=gu)
        rows[case]=dict(max_shape_error=max_err,max_active_error=max_act,max_unit_gauge_error=max_unit,
                        max_weighted_gauge_error=max_weighted,min_certified_chart_bound=min_chart,
                        max_unit_inner_force=max_res,max_label_adjustment=max_move,derivatives=detailed)
    if failed!=gate['failed'] or bool(gate['passed'])==bool(failed):
        raise ValueError('Independent failure list differs from the primary reduction')
    if failed and not evidence_only:
        raise ArithmeticError(f'Derivative gate failed: {failed}')
    verdict=('Both saved R4 preparations passed full and half two-sided reduced-work checks for the limiting active radial label and all six shape coordinates. The paired fixed-shape descent is eligible.'
             if not failed else 'The original reduced total-energy differences exceeded their tolerance at the listed widths. The original gate remains closed.')
    return dict(gate_sha256=digest(work/'GATE.json'),held_source_sha256=held_sha,rows=rows,
                failed=failed, passed=not failed,independent_review=verdict+' Every signed trial certified positive orientation and met the unscaled inner stationarity tolerance.',
                limitations='The derivatives hold phase, density and order nodes fixed. They do not establish held-state convergence, gauge independence at equilibrium, or shape-release readiness. Order stationarity after rematching refers to the physical projected residual with finite numerical mask changes accounted for separately.')


def precise_measure(scratch, held_sha, reduction_guard):
    work=scratch/'r4_barrier_work01'
    raw=measure(scratch,held_sha,evidence_only=True)
    if raw['passed']:
        raise ValueError('Independent precision adjudication requires an original failed gate')
    gate=read(work/'PRECISION_GATE.json')
    prec=read(work/'depth_precision.json')
    if (gate['absolute_tolerance']!=1.e-7 or gate['raw_gate_sha256']!=raw['gate_sha256']
            or gate['precision_measurement_sha256']!=digest(work/'depth_precision.json')):
        raise ValueError('Precision collection changed the original derivative tolerance')
    source=Path(__file__).resolve().parent
    for name,sha in gate['sources'].items():
        if digest(source/name)!=sha:
            raise ValueError(f'Changed precision evaluator: {name}')
    if prec['schemas']!=['order_gradient','double_well'] or not 0<prec['longdouble_epsilon']<2e-16:
        raise ValueError('Missing extended-precision components')
    if (max(abs(e) for e in prec['energy_anchor_delta'])>2.e-9
            or max(abs(a-b) for a,b in zip(prec['analytic'],prec['independent_80bit_analytic'],strict=True))>1.e-8):
        raise ValueError('Energy or analytic anchor changed')
    certified={}
    for case in ('primary','secondary'):
        for mode in raw['rows'][case]['derivatives']:
            folder=work/(case+'_'+mode+'_guard')
            certified[str(folder.relative_to(work)/'outcome.json')]=digest(folder/'outcome.json')
    for name,expected_exit in [('collection_guard',1),('depth_precision_readout_guard',0),(reduction_guard,0)]:
        guard(work/name,expected_exit)
        certified[name+'/outcome.json']=digest(work/name/'outcome.json')
    rows=raw['rows']
    failures=[]
    for case,entry in rows.items():
        limit=0.
        for mode,saved in entry['derivatives'].items():
            raw_bank=read(work/(case+'_'+mode)/'readout.json')
            banks=raw_bank['work']['steps']
            accepted=gate['cases'][case]['modes'][mode]
            if (accepted['work']['analytic']!=raw_bank['work']['analytic']
                    or accepted['source_sha256']!=raw_bank['source_sha256']
                    or accepted['fixed_trace_sha256']!=raw_bank['fixed_trace_sha256']):
                raise ValueError('Signed baseline or analytic derivative changed')
            corrections=[]
            residual=[]
            for s,r,claimed in zip(banks,saved['signed_central_minus_analytic'],accepted['work']['steps'],strict=True):
                h=s['h']
                if h!=claimed['h']:
                    raise ValueError('Original signed width changed')
                c=0.
                if mode=='log_depth':
                    part=prec['banks'][str(h)]
                    if set(part['signed'])!={'-1','1'}:
                        raise ValueError('Missing signed order energy')
                    for component in range(2):
                        native,extended=(fsum(v[field][component] for v in part['signed'].values())/2
                                         for field in ['native_work','extended_work'])
                        if max(abs(native-part['native'][component]),abs(extended-part['extended'][component]))>1e-13:
                            raise ValueError('Signed banks disagree with the recorded central component')
                    c=.5*fsum(u-v for bank in part['signed'].values()
                              for u,v in zip(bank['extended_work'],bank['native_work'],strict=True))
                    if abs(c-part['correction'])>1e-13:
                        raise ValueError('Could not reproduce the correction from the four component differences')
                delta=r+c
                for sign,bank in s['banks'].items():
                    expected=bank['one_sided']
                    if mode=='log_depth':
                        p=prec['banks'][str(h)]['signed'][sign]
                        expected+=fsum(v-u for v,u in zip(p['extended_work'],p['native_work'],strict=True))
                    if abs(claimed['banks'][sign]['one_sided']-expected)>1.e-13:
                        raise ValueError(f'Reported signed-bank correction cannot be reproduced in {case}, {mode}')
                if abs(delta-claimed['central_minus_analytic'])>1.e-13:
                    raise ValueError('Independently reconstructed central difference differs from the precision collection')
                if not isfinite(delta):
                    raise ValueError('Nonfinite corrected central error')
                if abs(delta)>1e-7:
                    failures.append([case,mode,h])
                limit=max(limit,abs(delta))
                residual.append(delta)
                corrections.append(c)
            saved.update(precision_central_minus_analytic=residual,signed_correction=corrections)
        entry['max_precision_absolute_error']=limit
    if bool(failures)==bool(gate['passed']):
        raise ValueError('Producer and referee differ on the finite-difference gate')
    if failures:
        raise ArithmeticError(f'Recomputed extended-precision gate failed: {failures}')
    return dict(gate_sha256=digest(work/'PRECISION_GATE.json'),raw_gate_sha256=raw['gate_sha256'],
                precision_readout_sha256=digest(work/'depth_precision.json'),held_source_sha256=held_sha,
                guards_sha256=certified,rows=rows,original_failed=raw['failed'],
                independent_review='All 56 signed perturbations passed the orientation and inner-stationarity checks. The 28 central differences agree with their analytic counterparts at both original widths within 1e-7 when the two extensive order terms use extended precision. The paired fixed-shape descent is eligible.',
                limitations=raw['limitations'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('held_source_sha')
    p.add_argument('output',type=Path)
    p.add_argument('--evidence-only',action='store_true')
    p.add_argument('--precision-reduction-guard')
    a=p.parse_args()
    value=(precise_measure(a.scratch,a.held_source_sha,a.precision_reduction_guard) if a.precision_reduction_guard
           else measure(a.scratch,a.held_source_sha,a.evidence_only))
    with a.output.open('x') as stream:
        json.dump(value,stream,allow_nan=False,indent=2)
        stream.write('\n')
    print(json.dumps(value,allow_nan=False),flush=True)
