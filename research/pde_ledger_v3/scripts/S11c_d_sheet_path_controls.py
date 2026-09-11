#!/usr/bin/env python3
"""Reclassify cached computed candidates and test bulk-root path operands."""
import argparse
import hashlib
import pickle
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.core.symbol import Str
import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import _restore


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--symbol-cache',type=Path,required=True)
    parser.add_argument('--transcript',type=Path,required=True)
    args=parser.parse_args()
    full,curl,units,known=pickle.loads(args.symbol_cache.read_bytes())[:4]
    symbols={v.name:v for v in known if isinstance(v,sp.Symbol)}
    symbols.update({v.name:v for v in full.free_symbols})
    dims=engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dims.known=dict(known)
    dims.unknown,dims.constraints,dims.solution={},set(),{}
    dims.zero=(sp.S.Zero,)*3
    r=SimpleNamespace(symbols=symbols,omega=symbols['omega'],ell=symbols['L_W'])
    engine.PHYSICAL_METADATA=engine.PhysicalMetadata(dims,r)
    modes=engine.FullPencilModes(SimpleNamespace(r=r,kn=symbols['s11cdSpectralNormalMomentum']),curl,units)
    _,relation,join=modes.analytic(full)
    engine.physical('SHEET_PATH_CONTROLS_SOURCE_RELATION',relation)
    provenance={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in (
        ('cache',args.symbol_cache),('transcript',args.transcript),('engine',Path(engine.__file__)),
        ('instrument',Path(__file__)))}
    engine.emit('SHEET_PATH_CONTROLS_PROVENANCE',provenance)
    engine.emit('METADATA_SHEET_PATH_CONTROLS_PROVENANCE',modes.numeric_metadata(engine.cas(provenance),lambda p:dims.zero))
    mappings={}
    collected=[]
    all_rows=[]
    all_loci=[]
    for line in args.transcript.open():
        name,_,value=line.partition(': ')
        if name.startswith('PY_S11CD_SPECTRAL_') and '_CARRIER_VALUES_' in name:
            kind,suffix=name[len('PY_S11CD_SPECTRAL_'):].split('_CARRIER_VALUES_',1)
            mappings[kind+'_'+suffix]=dict(_restore(value))
        if name.startswith('PY_S11CD_FULL_PENCIL_MODE_'):
            label=name[len('PY_S11CD_FULL_PENCIL_MODE_'):]
            mapping=mappings[label]
            bound=relation.xreplace(mapping)
            path=engine.BulkSheetPath(bound,modes.k,modes.q)
            slope=sp.cancel(-sp.diff(bound,modes.k)/sp.diff(bound,modes.q))
            derivative=sp.lambdify((modes.k,modes.q),slope,'numpy')
            rows=[]
            for index,item in enumerate(_restore(value)):
                record={str(k):v for k,v in item}
                k,q=complex(record['K']),complex(record['Q'])
                selected,details=path.classify(k,q)
                row={'index':index,'k':k,'q':q,
                     'previousLabel':record.get('PHYSICAL_BULK_SHEET',Str('UNCOMPUTED')),
                     'selectedLabel':selected,'pathDefined':details['PATH_DEFINED'],
                     'pathClearance':details['MINIMUM_BRANCH_POINT_DISTANCE']}
                if details['PATH_DEFINED']:
                    start=details['START_K']
                    delta=k-start
                    sol=solve_ivp(lambda u,v:np.asarray([delta*derivative(start+u*delta,v[0])]),
                                  (0,1),np.asarray([details['SEED_Q']]),rtol=1e-12,atol=1e-14)
                    row.update(odeStatus=sol.status,
                               transportDifference=details['REFINEMENTS'][-1]['END_Q']-sol.y[0,-1],
                               relativeTransportDifference=abs(details['REFINEMENTS'][-1]['END_Q']-sol.y[0,-1])/
                                   max(abs(details['REFINEMENTS'][-1]['END_Q']),abs(sol.y[0,-1])),
                               branchEquationResidual=complex(bound.subs({modes.k:k,modes.q:q})))
                rows.append(row)
            # Explicitly exercise the computed branch points and outward rays.
            loci=[]
            scale=max(abs(v) for v in path.points)
            for point in path.points:
                target=point+(0.5j*scale if point.imag>0 else -0.5j*scale if point.imag<0 else 0j)
                locus=path.transport(target)
                loci.append({'k':target,'pathDefined':locus['PATH_DEFINED'],
                             'pathClearance':locus['MINIMUM_BRANCH_POINT_DISTANCE']})
            packet={'label':label,'candidates':rows,'branchLocusProbes':loci}
            def unit(p):
                key=p[-1]
                if key in ('k','pathClearance'): return dims.measure(modes.k)
                if key in ('q','transportDifference'): return dims.measure(modes.q)
                if key=='branchEquationResidual': return tuple(2*v for v in dims.measure(modes.q))
                return dims.zero
            engine.emit('SHEET_PATH_CONTROLS_PACKET_'+str(len(collected)),packet)
            engine.emit('METADATA_SHEET_PATH_CONTROLS_PACKET_'+str(len(collected)),modes.numeric_metadata(engine.cas(packet),unit))
            collected.append(label)
            all_rows.extend(rows)
            all_loci.extend(loci)
    engine.emit('SHEET_PATH_CONTROLS_PACKET_LABELS',collected)
    engine.emit('METADATA_SHEET_PATH_CONTROLS_PACKET_LABELS',modes.numeric_metadata(engine.cas(collected),lambda p:dims.zero))
    resolved=[r for r in all_rows if 'transportDifference' in r]
    transitions={}
    for row in all_rows:
        pair=(str(row['previousLabel']),str(row['selectedLabel']))
        transitions[pair]=transitions.get(pair,0)+1
    summary={'packetCount':len(collected),'candidateCount':len(all_rows),
             'labelTransitions':[(a,b,n) for (a,b),n in sorted(transitions.items())],
             'odeStatuses':[r['odeStatus'] for r in resolved],
             'maximumTransportDifference':max(abs(r['transportDifference']) for r in resolved),
             'maximumRelativeTransportDifference':max(r['relativeTransportDifference'] for r in resolved),
             'maximumOriginalRadicalResidual':max(abs(r['branchEquationResidual']) for r in resolved),
             'branchLocusProbeCount':len(all_loci),
             'definedBranchLocusPaths':sum(bool(r['pathDefined']) for r in all_loci)}
    engine.emit('SHEET_PATH_CONTROLS_SUMMARY',summary)
    engine.emit('METADATA_SHEET_PATH_CONTROLS_SUMMARY',modes.numeric_metadata(engine.cas(summary),lambda p:
        dims.measure(modes.q) if p[0]=='maximumTransportDifference' else
        tuple(2*v for v in dims.measure(modes.q)) if p[0]=='maximumOriginalRadicalResidual' else dims.zero))
    engine.emit('SHEET_PATH_CONTROLS_DIMENSION_CONSTRAINTS',sorted(dims.constraints,key=str))
    engine.emit('METADATA_SHEET_PATH_CONTROLS_DIMENSION_CONSTRAINTS',modes.numeric_metadata(engine.cas(tuple(dims.constraints)),lambda p:dims.zero))


if __name__=='__main__':
    run()
