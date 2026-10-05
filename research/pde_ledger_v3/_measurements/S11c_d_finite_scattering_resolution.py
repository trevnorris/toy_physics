#!/usr/bin/env python3
"""Selected finite-response comparisons through the unchanged pilot constructor."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import S11c_d_finite_scattering as f

CHECKPOINT=f.M/'S11c_d_finite_scattering_checkpoint.json'
PLAN=f.M/'S11c_d_finite_scattering_resolution_plan.md'
CASES=(('source_profile',65,8,2),('collocation',97,8,2),('momentum',97,16,4))


def inspect(directory,source_joins=None):
    checks=json.loads((directory/'checks.json').read_text())
    for name,sha in checks['sourceFiles'].items():
        current=f.digest(f.ROOT/name)
        joined=source_joins is not None and name in source_joins and (sha,current)==tuple(source_joins[name])
        f.require(sha==f.digest(directory/'source'/name) and (current==sha or joined),('case source',name))
    for name,sha in checks['operandHashes'].items():
        f.require(f.digest(Path(name))==sha,('case input',name))
    for name,record in checks['artifacts'].items():
        f.require(f.digest(directory/name)==record['sha256'],('case artifact',name))
    system=f.unpickle(directory/'finite-system.pickle')
    solution=f.unpickle(directory/'finite-solution.pickle')
    matrix,rhs,coefficients=system['matrix'],system['rhs'],solution['coefficients']
    f.require(np.array_equal(matrix@coefficients-rhs,solution['equationResidual']),'saved equation replay')
    f.require(sorted(i for g in solution['groups'] for i in g['rowIndices'])==list(range(80)),'all native rows')
    direct=np.linalg.solve(matrix,rhs)
    size=len(system['nodes']);positions=np.linspace(-48.,48.,129)
    values=f.polynomial_basis(positions,48.,size)
    fields=np.stack([values@coefficients[j*size:(j+1)*size] for j in range(5)])
    scattering=solution['boundaryAnchoredFluxBasisScattering']
    current=solution['outgoingChannelCurrent'];end_currents=[]
    for end in ('LEFT','RIGHT'):
        indices=[i for i,(label,_) in enumerate(solution['outgoingLabels']) if label==end]
        amplitudes=scattering[indices]
        end_currents.append(np.diag(amplitudes.conj().T@current[np.ix_(indices,indices)]@amplitudes).real/solution['incomingFlux'])
    end_currents=np.array(end_currents)
    total_residual=end_currents.sum(axis=0)-solution['outgoingFluxRatio']
    f.require(np.max(abs(total_residual))<1e-12,'full-current end contraction')
    boundary=[];offset=0
    for end,index in [('LEFT',0),('RIGHT',size-1)]:
        channel=system['channels'][end];count=len(channel['incoming'])
        value=np.vstack([system['derivativeMatrices'][0][index]@coefficients[j*size:(j+1)*size] for j in range(5)])
        derivative=np.vstack([system['derivativeMatrices'][1][index]@coefficients[j*size:(j+1)*size] for j in range(5)])
        target=np.zeros_like(value);target[:,offset:offset+count]=channel['incomingBoundaryData'];offset+=count
        boundary.append(float(np.max(abs(derivative-channel['traceMap']@value-target))))
    return {'directory':str(directory),'checks':checks,'scattering':scattering,
        'current':current,'incomingFlux':solution['incomingFlux'],'endCurrentRatios':end_currents,
        'totalCurrentRatio':solution['outgoingFluxRatio'],'currentReconstructionResidual':total_residual,
        'modalAmplitudes':solution['modalAmplitudes'],'positions':positions,'fields':fields,
        'boundaryResidualMaxima':boundary,
        'independentSolveCoefficientDifference':float(np.max(abs(direct-coefficients))),
        'independentSolveResidual':float(np.max(abs(matrix@direct-rhs)))}


def compare(previous,current):
    f.require(np.array_equal(previous['current'],current['current']) and
        np.array_equal(previous['incomingFlux'],current['incomingFlux']),'unchanged channel currents')
    delta=current['scattering']-previous['scattering']
    currents=current['endCurrentRatios']-previous['endCurrentRatios']
    total=current['totalCurrentRatio']-previous['totalCurrentRatio']
    field=current['fields']-previous['fields']
    resolved=np.maximum(abs(current['scattering']),abs(previous['scattering']))>=0.01
    scale=np.maximum(abs(current['scattering']),abs(previous['scattering']))
    return {'scatteringDifference':delta,'endCurrentDifference':currents,
        'totalCurrentDifference':total,'commonGridFieldDifference':field,
        'modalDifferences':{end:current['modalAmplitudes'][end]-previous['modalAmplitudes'][end] for end in ('LEFT','RIGHT')},
        'summary':{'maximumAmplitudeChange':float(np.max(abs(delta))),
            'maximumResolvedAmplitudeRelativeChange':float(np.max(abs(delta[resolved])/scale[resolved])) if np.any(resolved) else None,
            'maximumEndCurrentChange':float(np.max(abs(currents))),
            'maximumTotalCurrentChange':float(np.max(abs(total))),
            'maximumCommonGridFieldChangeByField':np.max(abs(field),axis=(1,2)).tolist()}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE)
    base.mkdir(parents=True,exist_ok=False)
    checkpoint=json.loads(CHECKPOINT.read_text())
    f.require(checkpoint['status']=='VALIDATED_FINITE_INSTRUMENT_PILOT','accepted pilot')
    baseline=Path(checkpoint['runDirectory'])
    f.require(f.digest(baseline/'checks.json')==checkpoint['checksSha256'],'pilot checkpoint join')
    binding=Path(checkpoint['bindingReuse']['directory'])
    pins={str(p.relative_to(f.ROOT)):f.digest(p) for p in (CHECKPOINT,PLAN,Path(__file__),Path(f.__file__),f.ACCEPTANCE)}
    for name in pins:
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((f.ROOT/name).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'cases':CASES,'sourceOrder':256,'profileOrder':512,
        'relativeReportingTarget':0.01,'amplitudeResolution':1e-4,'currentRatioResolution':1e-6,
        'scope':'Selected finite response resolution; fixed physical box, regulator and approximate modal boundaries.'})
    started=time.monotonic();previous=inspect(baseline);records=[]
    environment=os.environ.copy()
    environment.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    for name,size,outer,panel in CASES:
        remaining=int(900-(time.monotonic()-started));f.require(remaining>0,'remaining selected-comparison time budget')
        child=base/name;child.mkdir();directory=child/'complete'
        command=[sys.executable,'-u',str(Path(f.__file__)),'--run-directory',str(directory),
            '--resume-binding',str(binding),'--size',str(size),'--outer',str(outer),'--panel',str(panel),
            '--source-order','256','--profile-order','512','--seconds',str(remaining)]
        f.save(child/'invocation.json',{'command':command,'threadEnvironment':{k:environment[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')}})
        then=time.monotonic()
        with (child/'stdout.txt').open('xb') as out,(child/'stderr.txt').open('xb') as err:
            run=subprocess.run(command,cwd=f.ROOT,env=environment,stdin=subprocess.DEVNULL,stdout=out,stderr=err)
        f.save(child/'outcome.json',{'exitCode':run.returncode,'wallSeconds':time.monotonic()-then,'stderrBytes':(child/'stderr.txt').stat().st_size})
        f.require(run.returncode==0 and (child/'stderr.txt').stat().st_size==0,('selected case completed',name))
        current=inspect(directory)
        f.require(current['checks']==json.loads((child/'stdout.txt').read_text()),'case stdout identity')
        difference=compare(previous,current)
        f.atomic_pickle(child/'comparison.pickle',{'previous':previous,'current':current,'difference':difference})
        records.append({'name':name,'difference':difference['summary'],'checks':current['checks'],
            'boundaryResidualMaxima':current['boundaryResidualMaxima'],
            'independentSolveCoefficientDifference':current['independentSolveCoefficientDifference'],
            'comparisonSha256':f.digest(child/'comparison.pickle')})
        f.save(base/'inventory.json',records);previous=current
    f.require(all(f.digest(f.ROOT/name)==sha for name,sha in pins.items()),'study sources unchanged')
    summary={'sourceFiles':pins,'records':records,'wallSeconds':time.monotonic()-started,
        'scope':'Finite resolution comparisons only; boundary/domain/regulator checks and continuum expansion remain.'}
    f.save(base/'checks.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
