#!/usr/bin/env python3
"""Bounded saved-pair regression for the uncompressed certificate route."""
import argparse
import faulthandler
import hashlib
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
STORE=ROOT.parents[1]/'_scratch/s11c'
DIAGNOSTIC=STORE/'s11c-source-quadrature-20260914/binding-diagnostic'
SOURCES=(Path(__file__).resolve(),ROOT/'_measurements/S11c_d_source_fourier_quadrature_check.py',
         ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py')


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def save(p,value):p.write_text(json.dumps(value,indent=2)+'\n')


def worker(base,test,index):
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    faulthandler.dump_traceback_later(150)
    import S11c_d_source_fourier_quadrature_check as c
    started=time.monotonic()
    packet=c.unpickle(DIAGNOSTIC/'pairs.pickle')
    row=next(r for r in packet['records'] if (r['test'],r['index'])==(test,index))
    left,right=row['current'],row['uses'][0]['expected']
    expanded=c.sp.Integral(c.sp.expand_mul(left.function),*left.limits)
    c.atomic_pickle(base/'operands.pickle',{'source':left,'expanded':expanded,'expected':right,
        'provenance':packet['provenance'],'dimensionState':packet['dimensionState'],
        'representations':tuple(c.sp.srepr(v) for v in (left,expanded,right))})
    certificate=c.engine.BoundedSourceFourierAssembly.reconstruction_certificate(
        expanded.function,right.function,shared=False)
    c.atomic_pickle(base/'certificate.pickle',certificate)
    proofs=tuple(certificate['REPLAY_RESIDUALS'])+tuple(v[1] for v in certificate['PHASE_SPLITS'].values())+tuple(
        v[2] for v in certificate['RADICAL_POWERS'].values())
    preliminary={'test':test,'index':index,'representationsDiffer':expanded!=right,
        'limitsEqual':expanded.limits==right.limits,'residualZero':certificate['RESIDUAL']==0,
        'proofCount':len(proofs),'nonzeroProofs':sum(v!=0 for v in proofs),'shared':False}
    save(base/'certificate-checks.json',preliminary)
    mutation=c.engine.BoundedSourceFourierAssembly.reconstruction_certificate(
        expanded.function,2*right.function,shared=False)
    c.atomic_pickle(base/'mutation.pickle',mutation)
    result=preliminary|{'coefficientMutationNonzero':mutation['RESIDUAL']!=0,
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.glob('*.pickle')}}
    save(base/'checks.json',result); print(json.dumps(result,indent=2))
    faulthandler.cancel_dump_traceback_later()
    if not result['limitsEqual'] or not result['residualZero'] or result['nonzeroProofs'] or not result['coefficientMutationNonzero']:
        raise ValueError('uncompressed binding certificate regression requires inspection')


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--test',type=int); parser.add_argument('--index',type=int)
    args=parser.parse_args(); base=args.run_directory.resolve(); base.relative_to(STORE)
    if args.test is not None:
        worker(base,args.test,args.index);return
    base.mkdir(parents=True,exist_ok=False)
    source_check=json.loads((DIAGNOSTIC/'checks.json').read_text())
    if digest(DIAGNOSTIC/'pairs.pickle')!=source_check['pairsSha256']:
        raise ValueError('diagnostic operand packet changed')
    cases=[(r['test'],r['index']) for r in source_check['records'] if not all(r['exact'])]
    if not cases:raise ValueError('no nonidentical live representation cases in diagnostic')
    pins={str(p.relative_to(ROOT)):digest(p) for p in SOURCES}
    for p in SOURCES:
        dest=base/'source'/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    save(base/'preflight.json',{'sourceFiles':pins,'pairsSha256':source_check['pairsSha256'],
        'caseSourceSha256':digest(DIAGNOSTIC/'checks.json'),'cases':cases,'workerWallLimitSeconds':180,
        'workerAddressSpaceLimitBytes':8*1024**3,'shared':False})
    inventory=[]; started=time.monotonic()
    for test,index in cases:
        directory=base/('case-'+str(test)+'-'+str(index));directory.mkdir()
        command=[sys.executable,'-u',str(Path(__file__).resolve()),'--run-directory',str(directory),
                 '--test',str(test),'--index',str(index)]
        item={'test':test,'index':index,'command':command};begin=time.monotonic()
        with (directory/'stdout').open('xb') as out,(directory/'stderr').open('xb') as err:
            try:
                child=subprocess.run(command,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=out,stderr=err,timeout=180)
                item.update(exitCode=child.returncode,status='exited')
            except subprocess.TimeoutExpired:
                item.update(exitCode=None,status='timeout')
        item.update(wallSeconds=time.monotonic()-begin,stderrBytes=(directory/'stderr').stat().st_size)
        if (directory/'checks.json').exists():item['checks']=json.loads((directory/'checks.json').read_text())
        inventory.append(item);save(base/'inventory.json',inventory)
        if item['exitCode']!=0 or item['stderrBytes']:
            raise ValueError(('bounded binding regression worker requires inspection',test,index,item['status']))
    result={'sourceFiles':pins,'pairsSha256':source_check['pairsSha256'],'cases':inventory,
        'wallSeconds':time.monotonic()-started,'shared':False,
        'scope':'Exact certificates on computed expanded representations of the saved finite-input operands and one-sided mutations; no numerical quadrature or upstream physics rerun.'}
    save(base/'checks.json',result)
    if pins!={str(p.relative_to(ROOT)):digest(p) for p in SOURCES}:raise ValueError('regression sources changed during run')
    if not any(r['checks']['representationsDiffer'] for r in inventory):raise ValueError('nonidentical representation path was not exercised')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
