#!/usr/bin/env python3
"""Synthetic algebraic-ansatz control for real-normal denominator crossings."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine
sp=engine.sp


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    # Dimensionless synthetic rational-matrix and radical-curve ansatz. These
    # placeholders are not inherited physical carriers or an upstream repair.
    w,k,q=sp.symbols('probeFrequency probeNormal probeRadical')
    physical=sp.ImmutableMatrix([[1/(q-2-sp.I*(w-1))]])
    relation=q**2-k**2-1
    data=engine.BulkExceptionalSlice.analyze(physical,relation,w,k,q)
    old={key:value for key,value in data['REAL_CONDITIONS'].items()
         if not key.startswith('DENOMINATOR_REAL_NORMAL_')}
    prior=engine.BulkExceptionalSlice.nonnegative_targets(old)
    current=engine.BulkExceptionalSlice.nonnegative_targets(data['REAL_CONDITIONS'])
    points=[]
    for frequency,labels in sorted(current.items(),key=lambda v:float(v[0])):
        for radical in sp.solve(data['DENOMINATOR'].subs(w,frequency),q):
            for normal in sp.solve(relation.subs({w:frequency,q:radical}),k):
                point={w:frequency,k:normal,q:radical}
                points.append({'frequency':str(frequency),'normal':str(normal),'radical':str(radical),
                    'normalRealityResidual':str(sp.simplify(sp.im(normal))),
                    'radicalResidual':str(sp.simplify(relation.subs(point))),
                    'denominatorResidual':str(sp.simplify(data['DENOMINATOR'].subs(point))),
                    'incidentConditions':labels})
    record={'inputKind':'SYNTHETIC_ALGEBRAIC_ANSATZ','physicalClaim':False,
        'ansatz':{'rationalMatrix':sp.srepr(physical),'radicalRelation':sp.srepr(relation)},
        'coefficientDimensionsLTM':[0,0,0],'coefficientGradeEpsEtaSigma':[0,0,0],'coefficientLambdaOrder':0,
        'priorConditionTargets':[(str(r),labels) for r,labels in prior.items()],
        'independentRealProjectionTargets':[(str(r),labels) for r,labels in current.items()],
        'newTargetCount':len(set(current)-set(prior)),
        'projectionReconstructionResiduals':[str(v) for v in data['REAL_NORMAL_PROJECTION']['RECONSTRUCTION_RESIDUALS']],
        'points':points,'engineSha256':hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
        'instrumentSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'newTargetCount':record['newTargetCount'],'pointCount':len(points),
        'projectionReconstructionResiduals':record['projectionReconstructionResiduals']}))


if __name__=='__main__':run()
