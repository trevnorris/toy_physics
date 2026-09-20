#!/usr/bin/env python3
"""Preserve repeated complete current checkpoints with exact payload joins."""
import argparse,ast,copy,json,marshal,hashlib,shutil,sys,types
from pathlib import Path
import S11c_d_remaining_case_currents_production as p
h,f=p.h,p.f
ORIGINAL=f.STORE/'s11c-remaining-case-currents-20260920/production'
PLAN=f.M/'S11c_d_remaining_case_currents_checkpoint_plan.md'
WRITES={}
FOCUSED=False


def checkpoint(path,value):
    # Every other packet retains the original strict create-once contract.
    if path.name not in ('acoustic.pickle','pairing.pickle') or not path.exists():
        return f.atomic_pickle(path,value)
    saved=f.unpickle(path);before=f.digest(path)
    f.require(h.same(saved,value),('repeated current checkpoint payload',str(path)))
    key=str(path.relative_to(p.STATE['base']))
    record=WRITES.setdefault(key,{'sha256':before,'exactPayloadJoins':0})
    f.require(record['sha256']==before,'unchanged repeated current packet')
    record['exactPayloadJoins']+=1
    evidence=path.with_name(path.stem+'-write-join.pickle')
    if not evidence.exists():f.atomic_pickle(evidence,{'saved':saved,'requested':value,'packetSha256':before})
    else:
        proof=f.unpickle(evidence)
        f.require(proof['packetSha256']==before and h.same(proof['saved'],saved) and h.same(proof['requested'],value),'unchanged checkpoint join operands')
    record['proofSha256']=f.digest(evidence)
    f.require(f.digest(path)==before,'original packet retained byte-for-byte')
    f.save(p.STATE['base']/'checkpoint-write-joins.json',WRITES)


def prepare(base):
    sources,accepted,paired,manifest=p.prepare(base)
    origin=ORIGINAL/'complete';old=json.loads((origin/'inputs.json').read_text())
    outcome=json.loads((ORIGINAL/'currents_construct.invocation.json').read_text())
    f.require(outcome==json.loads((ORIGINAL/'active.json').read_text()) and outcome['exitCode']==1,'actual original failed supervisor outcome')
    f.require(old['sourceFiles']==manifest['sourceFiles'],'unchanged original production/native sources')
    for name,value in old['sourceFiles'].items():f.require(f.digest(origin/'source'/name)==f.digest(f.ROOT/name)==value,'original/frozen/current producer join')
    for name,value in old['inputPackets'].items():f.require(f.digest(Path(name))==value,'original current input unchanged')
    for name,value in old['copiedInputs'].items():f.require(f.digest(origin/name)==f.digest(base/name)==value,'all original focused copies preserved')
    relative='cases/LAB_HELD__RHO4_CONSTANT/reference/acoustic.pickle'
    extra={str(path.relative_to(origin)) for path in origin.rglob('*.pickle') if 'source' not in path.relative_to(origin).parts and str(path.relative_to(origin)) not in old['copiedInputs']}
    f.require(extra=={relative},'complete original extra packet census before new construction')
    packet=f.unpickle(origin/relative)
    f.require(h.same(packet,accepted['LEFT']['acoustic']),'original saved acoustic packet joins accepted full input result')
    dest=base/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(origin/relative,dest)
    digest=f.digest(dest);f.require(digest==f.digest(origin/relative),'original completed current packet byte identity')
    for path in (origin/relative,origin/'inputs.json',ORIGINAL/'active.json',ORIGINAL/'currents_construct.invocation.json',ORIGINAL/'currents_construct.stderr',ORIGINAL/'currents_construct.stdout'):
        manifest['inputPackets'][str(path)]=f.digest(path)
    manifest['copiedInputs'][relative]=digest
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));manifest['sourceFiles'][name]=f.digest(path)
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    manifest['checkpointRecovery']={'originalRun':str(origin),'originalOutcome':outcome,'copiedCompletedPacket':relative,'copiedCompletedSha256':digest,
        'namespaceJoin':coordinator()[1],'writeJoins':WRITES,'scope':'Repeated acoustic/pairing writes require exact payload identity; no original current construction is repeated.'}
    if FOCUSED:manifest['scope']='Focused full LAB_HELD/RHO4 baseline packet/validation replay at its three backgrounds only; zero new acoustic/pairing constructions and no all-case production acceptance.'
    f.save(base/'inputs.json',manifest);p.STATE['manifest']=manifest
    if FOCUSED:sources={**sources,'sources':{'LAB_HELD__RHO4_CONSTANT':sources['sources']['LAB_HELD__RHO4_CONSTANT']}}
    return sources,accepted,paired,manifest


def coordinator():
    main,join=p.coordinator();code=hashlib.sha256(marshal.dumps(main.__code__)).hexdigest()
    proxy=types.SimpleNamespace(**dict(vars(f),atomic_pickle=checkpoint))
    f.require(all(getattr(proxy,key) is value for key,value in vars(f).items() if key!='atomic_pickle'),'all other packet/validation APIs unchanged')
    f.require(main.__globals__['f'] is f and main.__globals__['load'] is p.prepare,'original coordinator namespace bindings')
    main.__globals__['f']=proxy;main.__globals__['load']=prepare
    f.require(hashlib.sha256(marshal.dumps(main.__code__)).hexdigest()==code,'whole coordinator bytecode unchanged')
    return main,{'wholeCoordinatorAstJoin':join,'wholeCoordinatorBytecodeSha256':code,'unchangedBytecode':True,
        'namespaceOverrides':['f.atomic_pickle for complete acoustic/pairing packets only','load for exact saved original packet'],
        'originalHelperSha256':f.digest(Path(p.__file__)),
        'originalAtomicPickleSource':f.atomic_pickle.__code__.co_filename,
        'originalAtomicPickleSourceSha256':f.digest(Path(f.atomic_pickle.__code__.co_filename))}


def main():
    global FOCUSED
    FOCUSED='--check-writes' in sys.argv
    if FOCUSED:sys.argv.remove('--check-writes')
    f.require('--focused' not in sys.argv and '--check-reuse' not in sys.argv,'explicit checkpoint-recovery mode')
    coordinator()[0]()
    base=p.STATE['base'];checks=json.loads((base/'checks.json').read_text())
    if FOCUSED:
        f.require(checks['newAcousticConstructions']==checks['newTwoFrequencyConstructions']==0 and len(checks['caseInventory'])==3,'three complete original baseline validations only')
        mutations=[]
        for kind,key in (('acoustic','SLAB_SURFACE_DENSITY'),('pairing','SLAB_CURRENT_MATRIX')):
            path=base/'cases/LAB_HELD__RHO4_CONSTANT/right'/(kind+'.pickle');value=f.unpickle(path);before=f.digest(path)
            checkpoint(path,value)
            changed=dict(value);changed[key]=2*value[key]
            f.atomic_pickle(base/(kind+'-write-mutation.pickle'),{'original':value,'changed':changed,'key':key})
            try:checkpoint(path,changed)
            except ValueError:mutations.append(kind)
            else:raise ValueError('different current checkpoint accepted')
            f.require(f.digest(path)==before,'mutation leaves complete original packet unchanged')
        f.require(len(mutations)==2,'actual acoustic and pairing coefficient mutations')
        f.save(base/'write-controls.json',{'status':'ACCEPTED_CURRENT_CHECKPOINT_WRITE_REPAIR','caseInventory':checks['caseInventory'],'mutations':mutations,'writeJoins':WRITES,'namespaceJoin':coordinator()[1],
            'scope':'Full original baseline packet validation plus repeated-write controls. No new current or all-case production result.'})
    else:f.require(checks['newAcousticConstructions']==checks['newTwoFrequencyConstructions']==1,'only one genuinely new current/pairing family computed')


if __name__=='__main__':main()
