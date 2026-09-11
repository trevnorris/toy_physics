#!/usr/bin/env python3
"""Cached one-case implicit-pair coefficients and direct-pencil differences."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import sys
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import _restore


def array_data(value):
    if isinstance(value,np.ndarray):
        return sp.ImmutableMatrix(value)
    if isinstance(value,dict):
        return {k:array_data(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return tuple(array_data(v) for v in value)
    return value


def run():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--sample',type=int,default=0)
    parser.add_argument('--pullback',action='store_true',
                        help='mathematical parameter-coordinate check on the computed pencil')
    parser.add_argument('--production',action='store_true',help='run the production end-mode packet on the cached pencil')
    args=parser.parse_args()
    if args.production and args.pullback:
        parser.error('--production and --pullback select separate checks')
    manifest=json.loads(args.manifest.read_text())
    base=Path(manifest['run_directory'])
    cache=base/'symbols/RIGHT_LAB_HELD_RHO4_CONSTANT.pickle'
    transcript=base/'full.out'
    for path in (cache,transcript):
        if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['artifacts'][str(path.relative_to(base))]['sha256']:
            raise ValueError(('producer artifact mismatch',str(path)))
    for path,expected in manifest['source_hashes_after'].items():
        if path!='scripts/S11c_d_mixing_scattering_sympy_audit.py' and hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected:
            raise ValueError(('producer input changed',path))
    full,curl,units,known=pickle.loads(cache.read_bytes())[:4]
    symbols={s.name:s for s in known if isinstance(s,sp.Symbol)}
    symbols.update({s.name:s for s in full.free_symbols})
    dims=engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dims.known=dict(known)
    dims.unknown,dims.constraints,dims.solution={},set(),{}
    dims.zero=(sp.S.Zero,)*3
    r=SimpleNamespace(symbols=symbols,omega=symbols['omega'],ell=symbols['L_W'],
                      tangents=tuple(symbols['s11cdTangentialMomentum'+str(i)] for i in (1,2)))
    engine.PHYSICAL_METADATA=engine.PhysicalMetadata(dims,r)
    modes=engine.FullPencilModes(SimpleNamespace(r=r,kn=symbols['s11cdSpectralNormalMomentum']),curl,units)
    quotient,qunits=modes.quotient(full,'RECTANGULAR_CHECK')
    algebraic,relation,join=modes.analytic(quotient)
    mapping=modes.sample(algebraic,relation,args.sample)
    origin={modes.eta:sp.S.Zero,modes.sigma:sp.S.Zero}
    family=algebraic.xreplace(mapping)
    bound_relation=relation.xreplace(mapping)
    eta,sigma=modes.eta,modes.sigma
    pullback={}
    if args.pullback:
        eta,sigma=sp.symbols('s11cdJetCheckAlpha s11cdJetCheckBeta',real=True)
        dims.known.update({eta:dims.zero,sigma:dims.zero})
        pullback={modes.eta:eta+sigma,modes.sigma:sigma}
        family=family.xreplace(pullback)
        bound_relation=bound_relation.xreplace(pullback)
        origin={eta:sp.S.Zero,sigma:sp.S.Zero}
    jets=engine.RectangularModeJets(family,bound_relation,modes.k,modes.q,eta,sigma,origin)
    evaluated=sp.lambdify((modes.k,modes.q,eta,sigma),family,'numpy',cse=True)
    square=sp.solve(bound_relation,modes.q**2)[0]
    radical=sp.lambdify((modes.k,eta,sigma),square,'numpy',cse=True)
    prefix='PY_S11CD_FULL_PENCIL_MODE_PIT_RIGHT_LAB_HELD_RHO4_CONSTANT_'+str(args.sample)+': '
    packet=next(_restore(line.partition(': ')[2]) for line in transcript.open() if line.startswith(prefix))
    provenance={'producerManifestSha256':hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                'producerSources':manifest['source_hashes_after'],
                'engineSha256':hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
                'instrumentSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'sample':args.sample,'parameterPullback':pullback}
    engine.emit('RECTANGULAR_CHECK_PROVENANCE',provenance)
    engine.emit('METADATA_RECTANGULAR_CHECK_PROVENANCE',modes.numeric_metadata(engine.cas(provenance),lambda p:dims.zero))
    def finish():
        constraints=tuple(dims.constraints)
        engine.emit('RECTANGULAR_CHECK_DIMENSION_CONSTRAINTS',constraints)
        engine.emit('METADATA_RECTANGULAR_CHECK_DIMENSION_CONSTRAINTS',modes.numeric_metadata(engine.cas(constraints),lambda p:dims.zero))
        resources={'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        engine.emit('RECTANGULAR_CHECK_RESOURCES',resources)
        engine.emit('METADATA_RECTANGULAR_CHECK_RESOURCES',modes.numeric_metadata(engine.cas(resources),lambda p:
            (0,1,0) if p[0]=='wallSeconds' else dims.zero))
    if args.production:
        modes.solve_sample(algebraic,relation,qunits,'RIGHT_LAB_HELD_RHO4_CONSTANT',1,args.sample)
        finish()
        return
    engine.emit('RECTANGULAR_CHECK_PENCIL_GRADE_SUPPORT',tuple(
        (index,any(v!=0 for v in coefficient)) for index,coefficient in jets.pencil_coefficients.items()))
    engine.emit('METADATA_RECTANGULAR_CHECK_PENCIL_GRADE_SUPPORT',modes.numeric_metadata(engine.cas(tuple(
        (index,any(v!=0 for v in coefficient)) for index,coefficient in jets.pencil_coefficients.items())),lambda p:dims.zero))
    summary=[]
    scalars=[]
    for index,item in enumerate(packet):
        stored={str(k):v for k,v in item}
        k,q=complex(stored['K']),complex(stored['Q'])
        matrix=np.asarray(evaluated(k,q,0,0),dtype=complex)
        left,singular,rh=np.linalg.svd(matrix)
        n=int(np.sum(singular<1e-8*max(1.,singular[0])))
        if not n:
            continue
        right,dual=rh.conj().T[:,-n:],left[:,-n:]
        result=jets.construct(k,q,right,dual)
        values={'index':index,'k':k,'q':q,'nullity':n,'defined':result['DEFINED']}
        if result['DEFINED']:
            def at(series,a,b):
                return sum(v*a**i*b**j for (i,j),v in series.items())
            def direct(a,b):
                rp=at(result['RIGHT'],a,b)
                kp=at(result['K'],a,b)
                roots,vectors=np.linalg.eig(kp)
                columns=[]
                for col,root in enumerate(roots):
                    trial=complex(radical(root,a,b))**0.5
                    continued=min((trial,-trial),key=lambda v:abs(v-q))
                    columns.append(np.asarray(evaluated(root,continued,a,b))@rp@vectors[:,col])
                return np.column_stack(columns)@np.linalg.inv(vectors)
            differences=[]
            for step in (1e-3,5e-4,2.5e-4):
                operands=[direct(a,b) for a,b in ((step,step),(step,0),(0,step),(0,0))]
                mixed=(operands[0]-operands[1]-operands[2]+operands[3])/step**2
                differences.append((step,operands,mixed))
            values.update(rightResiduals={g:np.max(np.abs(d['EQUATION'])) for g,d in
                                          result['RIGHT_DIAGNOSTICS']['COEFFICIENT_RESIDUALS'].items()},
                          leftResiduals={g:np.max(np.abs(d['EQUATION'])) for g,d in
                                         result['LEFT_DIAGNOSTICS']['COEFFICIENT_RESIDUALS'].items()},
                          mixedK=result['K'][(1,1)],mixedQ=result['Q'].get((1,1)),
                          directPencilMixedDifferences=differences)
            scalars.append({'index':index,'rightMaximumCoefficientResidual':max(values['rightResiduals'].values()),
                            'leftMaximumCoefficientResidual':max(values['leftResiduals'].values()),
                            'mixedKMaximum':np.max(np.abs(result['K'][(1,1)])),
                            'mixedQMaximum':np.max(np.abs(result['Q'][(1,1)])),
                            'directPencilMixedResiduals':[np.max(np.abs(d[2])) for d in differences]})
        else:
            values.update(rightStatus=result['RIGHT_DIAGNOSTICS'],leftStatus=result['LEFT_DIAGNOSTICS'])
        values=array_data(values)
        engine.emit('RECTANGULAR_CHECK_CANDIDATE_'+str(index),values)
        # Residuals are coefficients in the declared input numerical unit frame;
        # restored row/column units are emitted separately with the quotient.
        engine.emit('METADATA_RECTANGULAR_CHECK_CANDIDATE_'+str(index),modes.numeric_metadata(engine.cas(values),lambda p:
            dims.measure(modes.k) if p[0] in ('k','mixedK') else
            dims.measure(modes.q) if p[0] in ('q','mixedQ') else dims.zero))
        summary.append((index,n,result['DEFINED']))
    engine.emit('RECTANGULAR_CHECK_SUMMARY',summary)
    engine.emit('METADATA_RECTANGULAR_CHECK_SUMMARY',modes.numeric_metadata(engine.cas(summary),lambda p:dims.zero))
    engine.emit('RECTANGULAR_CHECK_SCALARS',scalars)
    engine.emit('METADATA_RECTANGULAR_CHECK_SCALARS',modes.numeric_metadata(engine.cas(scalars),lambda p:
        dims.measure(modes.k) if p[-1]=='mixedKMaximum' else
        dims.measure(modes.q) if p[-1]=='mixedQMaximum' else dims.zero))
    aggregate={'candidateCount':len(summary),'definedCount':sum(bool(v[2]) for v in summary),
               'maximumRightCoefficientResidual':max(v['rightMaximumCoefficientResidual'] for v in scalars),
               'maximumLeftCoefficientResidual':max(v['leftMaximumCoefficientResidual'] for v in scalars),
               'maximumMixedK':max(v['mixedKMaximum'] for v in scalars),
               'maximumMixedQ':max(v['mixedQMaximum'] for v in scalars),
               'maximumDirectMixedResidualByStep':[max(v['directPencilMixedResiduals'][i] for v in scalars) for i in range(3)]}
    engine.emit('RECTANGULAR_CHECK_AGGREGATE',aggregate)
    engine.emit('METADATA_RECTANGULAR_CHECK_AGGREGATE',modes.numeric_metadata(engine.cas(aggregate),lambda p:
        dims.measure(modes.k) if p[0]=='maximumMixedK' else
        dims.measure(modes.q) if p[0]=='maximumMixedQ' else dims.zero))
    finish()


if __name__=='__main__':
    run()
