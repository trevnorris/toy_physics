#!/usr/bin/env python3
"""Exact arithmetic checks of emitted axis certificates and full-space gates."""
from collections import Counter
import sympy as sp


def check_certificate(certificate, modes):
    operands = certificate['OPERANDS']
    if not operands['REAL_AXIS_CHART_DEFINED']:
        raise ValueError('unresolved reality chart')
    nonzero = {key:value for key,value in certificate['CHECKS'].items() if value!=0}
    if nonzero:
        raise ValueError(('reality reconstruction residuals',nonzero))
    x = operands['RADICAL_COORDINATE']
    q = next(iter(operands['POLYNOMIAL_COORDINATE'].free_symbols))
    factors = sp.sqf_list(sp.Poly(operands['POLYNOMIAL_COORDINATE'],q))[1]
    joined = Counter()
    root_count = 0
    for axis in certificate['AXES']:
        factor = factors[axis['FACTOR_INDEX']][0].monic()
        multiplier = sp.I if axis['AXIS']=='IMAGINARY' else sp.S.One
        real,imag = (sp.Poly(v,x,domain=sp.QQ) for v in
                     sp.expand_complex(factor.as_expr().subs(q,multiplier*x)).as_real_imag())
        common = sp.gcd(real,imag).monic()
        if common.as_expr()!=axis['GCD_COORDINATE'] or common.degree()!=axis['GCD_DEGREE']:
            raise ValueError('axis gcd differs from actual factor')
        if common.count_roots(-sp.oo,sp.oo)!=len(axis['ROOTS']):
            raise ValueError('incomplete real-axis root inventory')
        normal = sp.Poly(axis['NORMAL_SQUARE_COORDINATE'],x)
        for root in axis['ROOTS']:
            lo,hi = root['RADICAL_INTERVAL']
            if common.count_roots(lo,hi)!=1 or not root['UNIQUE_DISK_JOIN']:
                raise ValueError('axis interval/root disk join unresolved')
            disk_index,margins = root['DISK_JOINS'][0]
            disk = certificate['DISKS'][disk_index]
            expected = tuple(sp.expand(disk['RADIUS']**2-sp.Abs(multiplier*v-disk['CENTER'])**2) for v in (lo,hi))
            if margins!=expected or any(v<=0 for v in margins):
                raise ValueError('axis interval not strictly contained in native disk')
            if disk['FACTOR_INDEX']!=axis['FACTOR_INDEX']:
                raise ValueError('axis factor/disk mismatch')
            extrema = [normal.eval(lo),normal.eval(hi)]
            extrema += [normal.eval(v) for v in sp.solve(normal.diff().as_expr(),x) if lo<=v<=hi]
            bound_lo,bound_hi = root['NORMAL_SQUARE_INTERVAL']
            if not bound_lo<=min(extrema)<=max(extrema)<=bound_hi:
                raise ValueError('normal-square interval does not enclose exact extrema')
            threshold = bool(sp.gcd(common,normal).count_roots(lo,hi))
            sign = 0 if threshold else 1 if min(extrema)>0 else -1 if max(extrema)<0 else None
            if sign is None or sign!=root['NORMAL_SQUARE_SIGN']:
                raise ValueError('normal-square sign unresolved or incorrect')
            root_count+=1;joined[disk_index]+=1
    for disk in certificate['DISKS']:
        status = disk['NORMAL_REALITY_STATUS']
        if status=='UNRESOLVED':
            raise ValueError(('unresolved disk reality',disk['ROOT_DISK_INDEX']))
        if not joined[disk['ROOT_DISK_INDEX']] and not (
                disk['REAL_AXIS_CLEARANCE']>0 and disk['IMAGINARY_AXIS_CLEARANCE']>0 and status=='PROVED_NONREAL'):
            raise ValueError('disk has neither an axis certificate nor two axis exclusions')
    normalized = []
    for mode in modes:
        disk = certificate['DISKS'][mode['ROOT_DISK_INDEX']]
        if mode['NORMAL_REALITY_CERTIFICATE']!=disk:
            raise ValueError('candidate reality certificate join')
        if mode['EXACT_REAL_NORMAL']!=(disk['NORMAL_REALITY_STATUS']=='PROVED_REAL'):
            raise ValueError('candidate exact-reality gate differs')
        if mode.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False):
            if not (disk['NORMAL_REALITY_STATUS']=='PROVED_REAL' and disk['NORMAL_NONZERO'] is True and
                    disk['DENOMINATOR_EXCLUDED'] and mode['SHEET_MEMBERSHIP']==sp.true and
                    mode['BULK_DECAY_DISK_CERTIFIED'] and
                    mode['NORMAL_PAIRING_RANK']==mode['FREQUENCY_PAIRING_RANK']==mode['CURRENT_RANK']==mode['NULLITY']):
                raise ValueError('normalization escaped its independent domain gates')
            normalized.append({'index':mode['INDEX'],'nullity':mode['NULLITY'],
                'k':str(mode['K']),'nativeQ':str(mode['Q']),
                'currentEigenvalues':[float(v) for v in mode['FORMS']['CURRENT_EIGENVALUES']]})
    return {'diskStatuses':dict(Counter(v['NORMAL_REALITY_STATUS'] for v in certificate['DISKS'])),
        'candidateStatuses':dict(Counter(v['NORMAL_REALITY_CERTIFICATE']['NORMAL_REALITY_STATUS'] for v in modes)),
        'axisSets':[{key:axis[key] for key in ('FACTOR_INDEX','AXIS','GCD_DEGREE','REAL_ROOT_COUNT')} for axis in certificate['AXES']],
        'checkedAxisRootIntervals':root_count,'literalCertificateResiduals':len(certificate['CHECKS']),
        'normalizedSubspaces':normalized,'normalizedBasisDirections':sum(v['nullity'] for v in normalized)}


def exceptional_fixtures(certificate):
    """Algorithm-only fixtures derived from the wave's two exceptional factors."""
    from copy import deepcopy
    from S11c_d_modal_subspace_check import engine
    operands = certificate['OPERANDS']
    q = next(iter(operands['POLYNOMIAL_COORDINATE'].free_symbols))
    curve = operands['BOUND_WAVE_COORDINATE']
    k = next(symbol for symbol in curve.free_symbols if symbol!=q)
    threshold = sp.Poly(curve.subs(k,0),q)
    branch = sp.Poly(sp.diff(curve,q),q)
    polynomial = threshold*branch
    coverage,_ = engine.EndSpectrumCoverage.isolate(polynomial)
    result = engine.NormalRealityCoverage.construct(polynomial,curve,branch.as_expr(),k,q,coverage)
    axis_summary = check_certificate(result,[])
    origins = sum(root['ORIGIN_ROOT'] for axis in result['AXES'] for root in axis['ROOTS'])
    thresholds = [disk for disk in result['DISKS'] if disk['NORMAL_NONZERO'] is False]
    broken = deepcopy(coverage)
    broken['ROOT_DISKS'][0]['ONE_ROOT_DISK'] = False
    unresolved = engine.NormalRealityCoverage.construct(polynomial,curve,branch.as_expr(),k,q,broken)
    checks = {'originRecordedOnBothAxes':origins==2,
        'thresholdRootsKeptZero':len(thresholds)==len(sp.solve(threshold.as_expr(),q)),
        'thresholdNormalsReal':all(disk['NORMAL_REALITY_STATUS']=='PROVED_REAL' for disk in thresholds),
        'denominatorIntersectionNotExcluded':not any(disk['DENOMINATOR_EXCLUDED'] for disk in result['DISKS']),
        'uncertifiedDiskRemainsUnresolved':unresolved['DISKS'][0]['NORMAL_REALITY_STATUS']=='UNRESOLVED'}
    return {'scope':'ALGORITHM_FIXTURES_NOT_PHYSICAL_SPECTRUM','axisSummary':axis_summary,
            'originAxisRecords':int(origins),'thresholdDiskCount':len(thresholds),'checks':checks}


def main():
    import argparse,hashlib,json,pickle
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('--controls-directory',type=Path,required=True)
    parser.add_argument('--reference-directory',type=Path)
    args = parser.parse_args()
    # Import the pinned engine module before restoring its SymPy/Python payloads.
    from S11c_d_modal_subspace_check import engine
    results = {}
    for path in sorted(args.controls_directory.glob('*/checks.json')):
        report = json.loads(path.read_text())
        payload = path.parent/'objects.pickle'
        if hashlib.sha256(payload.read_bytes()).hexdigest()!=report['objectsSha256']:
            raise ValueError('control payload checksum')
        objects = pickle.loads(payload.read_bytes())
        modes = objects['modes']
        results[path.parent.name] = check_certificate(modes['NORMAL_REALITY_COVERAGE'],modes['RECORDS'])
    if args.reference_directory:
        report = json.loads((args.reference_directory/'checks.json').read_text())
        payload = args.reference_directory/'objects.pickle'
        if hashlib.sha256(payload.read_bytes()).hexdigest()!=report['objectsSha256']:
            raise ValueError('reference payload checksum')
        modes,_ = pickle.loads(payload.read_bytes())
        results['originalReference'] = check_certificate(modes['NORMAL_REALITY_COVERAGE'],modes['RECORDS'])
    fixtures = exceptional_fixtures(modes['NORMAL_REALITY_COVERAGE'])
    print(json.dumps({'physicalPackets':results,'exceptionalAlgorithmFixtures':fixtures},indent=2))
    if not all(fixtures['checks'].values()):raise ValueError('exceptional algorithm fixture')


if __name__=='__main__':
    main()
