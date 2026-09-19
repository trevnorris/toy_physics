#!/usr/bin/env python3
"""Finish the selected response study from completed cases and native layout sums."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import S11c_d_finite_scattering_resolution as study

f=study.f
REPAIR=f.M/'S11c_d_finite_scattering_checkpoint_repair.json'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--original',type=Path,required=True)
    args=parser.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE)
    original=args.original.resolve();original.relative_to(f.STORE)
    base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    repair=json.loads(REPAIR.read_text())
    helper=str(Path(f.__file__).resolve().relative_to(f.ROOT))
    joins={helper:(repair['oldFiniteSha256'],repair['newFiniteSha256'])}
    f.require(f.digest(Path(f.__file__))==repair['newFiniteSha256'] and
        f.digest(Path(study.__file__))==repair['newStudySha256'],'tested helper revisions')
    pins={str(p.resolve().relative_to(f.ROOT)):f.digest(p) for p in
          (Path(__file__),Path(study.__file__),Path(f.__file__),REPAIR,study.PLAN,study.CHECKPOINT,f.ACCEPTANCE)}
    for path,record in repair['savedArtifacts'].items():
        f.require(f.digest(Path(path))==record['sha256'],('saved original',path))
    for name,sha in json.loads((original/'inputs.json').read_text())['sourceFiles'].items():
        f.require(f.digest(original/'source'/name)==sha,('frozen study source',name))
        current=f.digest(f.ROOT/name)
        if name==helper:f.require((sha,current)==joins[name],'finite helper join')
        elif name==str(Path(study.__file__).resolve().relative_to(f.ROOT)):
            f.require((sha,current)==(repair['oldStudySha256'],repair['newStudySha256']),'study source guard join')
        else:f.require(current==sha,('unchanged study source',name))
    for name in pins:
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((f.ROOT/name).read_bytes())
    records=json.loads((original/'inventory.json').read_text())
    f.require([v['name'] for v in records]==['source_profile','collocation'],'completed original case census')
    for record in records:
        child=original/record['name'];directory=child/'complete'
        outcome=json.loads((child/'outcome.json').read_text())
        f.require(outcome['exitCode']==0 and (child/'stderr.txt').stat().st_size==0,'completed original child')
        current=study.inspect(directory,joins)
        f.require(current['checks']==record['checks']==json.loads((child/'stdout.txt').read_text()),'original case stdout identity')
        f.require(f.digest(child/'comparison.pickle')==record['comparisonSha256'],'saved comparison hash')
        saved=f.unpickle(child/'comparison.pickle')
        f.require(study.compare(saved['previous'],current)['summary']==record['difference'],'original comparison replay')
    previous=current
    checkpoint=json.loads(study.CHECKPOINT.read_text());binding=Path(checkpoint['bindingReuse']['directory'])
    f.save(base/'inputs.json',{'sourceFiles':pins,'sourceJoins':joins,'originalDirectory':str(original),
        'reusedCases':[v['name'] for v in records],'savedMomentumNodes':repair['savedMomentumNodes'],
        'remainingChildBudgetSeconds':800})
    child=base/'momentum';child.mkdir();directory=child/'complete'
    name,size,outer,panel=study.CASES[-1]
    command=[sys.executable,'-u',str(Path(f.__file__)),'--run-directory',str(directory),
        '--resume-binding',str(binding),'--resume-layouts',str(original/'momentum/complete'),
        '--size',str(size),'--outer',str(outer),'--panel',str(panel),
        '--source-order','256','--profile-order','512','--seconds','800']
    environment=os.environ.copy()
    environment.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    f.save(child/'invocation.json',{'command':command,'threadEnvironment':{k:environment[k] for k in
        ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')}})
    then=time.monotonic()
    with (child/'stdout.txt').open('xb') as out,(child/'stderr.txt').open('xb') as err:
        run=subprocess.run(command,cwd=f.ROOT,env=environment,stdin=subprocess.DEVNULL,stdout=out,stderr=err)
    f.save(child/'outcome.json',{'exitCode':run.returncode,'wallSeconds':time.monotonic()-then,'stderrBytes':(child/'stderr.txt').stat().st_size})
    f.require(run.returncode==0 and (child/'stderr.txt').stat().st_size==0,'remaining momentum case completed')
    current=study.inspect(directory)
    f.require(current['checks']==json.loads((child/'stdout.txt').read_text()),'remaining case stdout identity')
    difference=study.compare(previous,current)
    f.atomic_pickle(child/'comparison.pickle',{'previous':previous,'current':current,'difference':difference})
    records.append({'name':name,'difference':difference['summary'],'checks':current['checks'],
        'boundaryResidualMaxima':current['boundaryResidualMaxima'],
        'independentSolveCoefficientDifference':current['independentSolveCoefficientDifference'],
        'comparisonSha256':f.digest(child/'comparison.pickle')})
    f.save(base/'inventory.json',records)
    for path,record in repair['savedArtifacts'].items():
        f.require(f.digest(Path(path))==record['sha256'],('original operand post hash',path))
    f.require(all(f.digest(f.ROOT/name)==sha for name,sha in pins.items()),'recovery sources unchanged')
    summary={'sourceFiles':pins,'sourceJoins':joins,'records':records,'wallSeconds':time.monotonic()-started,
        'reusedMomentumNodes':repair['savedMomentumNodes'],
        'scope':'Finite resolution comparisons only; original completed cases and native partial sums reused; boundary/domain/regulator and continuum work remain.'}
    f.save(base/'checks.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
