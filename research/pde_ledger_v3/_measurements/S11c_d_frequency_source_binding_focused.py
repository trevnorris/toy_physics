#!/usr/bin/env python3
"""Actual saved/fresh source pairs and a bounded representation regression."""
import argparse,ast,copy,json,resource,signal,time
from pathlib import Path
import sympy as sp
import S11c_d_frequency_source as q

f=q.f
PREVIOUS=f.STORE/'s11c-frequency-source-20260919/retry-01/complete'
ORIGINAL=f.STORE/'s11c-frequency-source-20260919/binding-repair/original.py'


def ast_join():
    before=ast.parse(ORIGINAL.read_text());after=ast.parse(Path(q.__file__).read_text())
    old={v.name:v for v in before.body if isinstance(v,ast.FunctionDef)}
    new={v.name:v for v in after.body if isinstance(v,ast.FunctionDef)}
    unchanged=('load','census','end_sources')
    for name in unchanged:f.require(ast.dump(old[name])==ast.dump(new[name]),('unchanged frequency helper',name))
    # The entire source calculation is retained inside the new cache branch.
    def source_loop(fn):return next(v for v in fn.body if isinstance(v,ast.For) and isinstance(v.target,ast.Tuple) and [getattr(t,'id',None) for t in v.target.elts]==['key','item'])
    old_loop=source_loop(old['sources']);new_loop=source_loop(new['sources'])
    original_calculation=[]
    for statement in old_loop.body:
        if isinstance(statement,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='path' for t in statement.targets):break
        if isinstance(statement,ast.AugAssign) and isinstance(statement.target,ast.Name) and statement.target.id=='frequency_controls':continue
        original_calculation.append(statement)
    f.require(ast.dump(ast.Module(body=original_calculation,type_ignores=[]))==ast.dump(ast.Module(body=new_loop.body[0].orelse,type_ignores=[])),'entire actual source/derivative calculation unchanged')
    # Original output and final validation are retained; only two new tags
    # (certified residuals and explicit representation/domain provenance) enter.
    old_emit=copy.deepcopy(old['emit_result']);new_emit=copy.deepcopy(new['emit_result'])
    loop=next(v for v in new_emit.body if isinstance(v,ast.For))
    loop.body=[v for v in loop.body if not(isinstance(v,ast.For) and isinstance(v.target,ast.Tuple) and [getattr(t,'id',None) for t in v.target.elts]==['name','comparison'])]
    f.require(ast.dump(old_emit)==ast.dump(new_emit),'original physical emitter unchanged')
    # Everything after result construction in main, including the full metadata
    # and hash guard, is unchanged.
    def tail(fn):
        start=next(i for i,v in enumerate(fn.body) if isinstance(v,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='result' for t in v.targets))
        return ast.dump(ast.Module(body=fn.body[start:],type_ignores=[]))
    f.require(tail(old['main'])==tail(new['main']),'entire final result/output/validation tail unchanged')
    return {'unchangedWholeHelpers':unchanged,'sourceCalculationAstEqual':True,'originalEmitterAstEqual':True,'finalValidationTailAstEqual':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('bounded binding regression; saved pairs remain unresolved')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);started=time.monotonic();joins=ast_join()
    data=q.load(base);inventory=json.loads((PREVIOUS/'record-inventory.json').read_text());cuts=data['old']['domain-binding.pickle']['bound']['cutoffBindings'];origin=data['adapter'].input.origin;w=data['frequency'];records={};nonidentical=0;proof_count=0
    for key,value in inventory.items():
        original=data['packet']['records'][key]['record']['ORIGINAL'];saved=f.unpickle(PREVIOUS/value['path']);f.require(f.digest(PREVIOUS/value['path'])==value['sha256'] and original==saved['original'],'saved source hash and original expression')
        live=q.engine.memo_xreplace(data['live'].bind(original),cuts)
        for label,adapter in (('seed',data['adapter']),('control',data['control'])):
            expected=q.engine.memo_xreplace(adapter.bind(original),cuts);actual=live.subs({w:adapter.input.parameters['omega'],**origin},simultaneous=True)
            comparison=q.binding_comparison(base/'fresh-pairs'/key/label,actual,expected,saved['unit'])
            records[key,label]=comparison;nonidentical+=int(not comparison['operands']['exactEqualLive']);proof_count+=len(comparison['proofResiduals'])
    key='local0Row3Column0';saved=f.unpickle(PREVIOUS/inventory[key]['path'])
    expanded=sp.expand(saved['controlActual']);f.require(expanded!=saved['controlExpected'],'actual nonidentical representation coverage')
    comparison=q.binding_comparison(base/'expanded-control',expanded,saved['controlExpected'],saved['unit'])
    f.require(comparison['certificate'] is not None and comparison['mutation']['RESIDUAL']!=0,'actual expanded-pair certificate and mutation')
    f.atomic_pickle(base/'fresh-comparisons.pickle',records)
    wrong=dict(data,pins=dict(data['pins']),operands=dict(data['operands']),manifest=dict(data['manifest'],input={}))
    completed=q.resume_sources(base,data,PREVIOUS)
    f.require(len(completed)==16 and f.digest(base/'baseline-binding.pickle')==f.digest(PREVIOUS/'baseline-binding.pickle'),'actual exact saved binding reuse')
    for key,record in completed.items():
        f.require(record==f.unpickle(PREVIOUS/inventory[key]['path']),'complete restored record identity')
        for label,left,right in (('seed',record['boundAtReference'],record['acceptedBinding']),('control',record['controlActual'],record['controlExpected'])):
            q.binding_comparison(base/'restored-pairs'/key/label,left,right,record['unit'],restored=True)
    # A changed actually consumed input must fail before any copied source can
    # be offered for reuse. The original run and its data are never edited.
    bad=base/'wrong-input';bad.mkdir()
    rejected=False
    try:q.resume_sources(bad,wrong,PREVIOUS)
    except ValueError as error:rejected=str(error)=='unchanged inputs and finite settings'
    f.require(rejected,'changed actual input rejected')
    checks={'status':'PASSED','astJoins':joins,'originalSourceSha256':f.digest(ORIGINAL),'repairedSourceSha256':f.digest(Path(q.__file__)),
      'savedRecords':len(inventory),'freshPairs':len(records),'freshNonidenticalPairs':nonidentical,'freshProofScalars':proof_count,
      'expandedProofScalars':len(comparison['proofResiduals']),'expandedNormalizedResidualZero':comparison['normalizedResidual']==0,'expandedMutationNonzero':comparison['mutation']['RESIDUAL']!=0,
      'restoredPairs':2*len(completed),'reusedBaselineSha256':f.digest(base/'baseline-binding.pickle'),'changedInputRejected':rejected,
      'sourceFiles':data['pins'],'inputPackets':data['operands'],'artifacts':{str(p.relative_to(base)):f.digest(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts},
      'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'scope':'Instrument comparison and exact saved-operand reuse only; no numerical quadrature or frequency pole result.'}
    for name,sha in data['pins'].items():f.require(f.digest(f.ROOT/name)==sha,'post source hash')
    for name,sha in data['operands'].items():f.require(f.digest(Path(name))==sha,'post input hash')
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
