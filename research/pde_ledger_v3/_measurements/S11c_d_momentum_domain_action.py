#!/usr/bin/env python3
"""Finite momentum-box actions with an explicit matching-rule baseline."""
import argparse,ast,copy,hashlib,inspect,json,resource,shutil,textwrap,time,types
from pathlib import Path
import numpy as np
import S11c_d_momentum_domain_prepare as prepare
position=prepare.position;parallel=prepare.parallel;native=prepare.native;momentum=prepare.momentum;engine=prepare.engine
ROOT,STORE=prepare.ROOT,prepare.STORE
require,digest,save,atomic_pickle,unpickle=prepare.require,prepare.digest,prepare.save,prepare.atomic_pickle,prepare.unpickle
artifact=position.artifact;assemble=position.assemble
M=ROOT/'_measurements';CHECKPOINT=M/'S11c_d_momentum_domain_prepare_checkpoint.json'
PLAN=M/'S11c_d_momentum_domain_action_plan.md';PREFIX='MOMENTUM_DOMAIN_ACTION_LAB_HELD_RHO4_CONSTANT'
CUTOFFS=(2.,2.,3.,4.)


def settings(reference,index,smoke=False):
    rules={}
    for n,value in reference['layoutSettings'].items():
        s=dict(value,momentumBound=CUTOFFS[index],sourceBound=48.,profileBound=14.,sourceNodes=256,profileNodes=512 if index else 256)
        s['innerOrders']=tuple(value.get('innerOrders',(value['panelOrder'],)*(n-1)))
        if smoke and index:s.update(outerOrder=2,panelOrder=1,innerOrders=(1,)*(n-1),sourceNodes=16,profileNodes=16)
        rules[n]=s
    return rules


def load(base):
    accepted,previous=momentum.source.accepted(CHECKPOINT)
    for n,h in accepted['sourceFiles'].items():require(digest(ROOT/n)==h==digest(previous/'source'/n),'accepted source')
    for a in accepted['workerArtifacts']:require(digest(previous/a['path'])==a['sha256'],'accepted worker/transform/probe artifact')
    bp=unpickle(previous/'momentum-domains.pickle');pp=unpickle(previous/'momentum-domain-preparation.pickle');original=unpickle(previous/'accepted-position-domain.pickle');original_bound=unpickle(previous/'accepted-bound-position-domains.pickle')
    require(bp['sourceFiles']==pp['sourceFiles']==accepted['sourceFiles'] and bp['provenance']==pp['provenance']==accepted['provenance'],'preparation packet joins')
    _,sb=momentum.source.accepted(momentum.source.SOURCE);r,dimensions=momentum.source.native.source.restore_context(unpickle(sb/'reduced-action.pickle'));dimensions.__dict__.update(pp['dimensionState'])
    domains={i:bp['domains'][j] for i,j in enumerate((0,0,1,2))};joins={i:bp['joins'][j] for i,j in enumerate((0,0,1,2))};finest={}
    for ti in range(2):
        ref=next(v for v in original['result']['records'] if v['test']==ti and v['index']==2)
        finest[ti]=dict(original_bound['references'][ti],action=ref['action'],groups=ref['groups'],contributions=[{'index':v['index'],'terms':v['values']} for v in ref['terms']],layoutSettings=ref['settings'],setting=ref['settings'][3])
        require(settings(finest[ti],0)==ref['settings'],'exact retained rule baseline')
        action=assemble(r,domains[0],finest[ti],ref['groups']);require(np.array_equal(action['action'],ref['action']) and all(not np.any(v['differences']) for v in action['terms']),'retained action/term reconstruction')
    pins=dict(accepted['sourceFiles'])
    for path in (CHECKPOINT,PLAN,Path(__file__).resolve(),M/'S11c_d_momentum_domain_action_preflight.py'):pins[str(path.relative_to(ROOT))]=digest(path)
    for n in pins:
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,dest)
    for name in ('momentum-domains.pickle','momentum-domain-preparation.pickle','accepted-position-domain.pickle','accepted-bound-position-domains.pickle'):shutil.copyfile(previous/name,base/('accepted-'+name if not name.startswith('accepted-') else name))
    provenance=dict(accepted['provenance'],MOMENTUM_PREPARATION_CHECKPOINT_SHA256=digest(CHECKPOINT),MOMENTUM_DOMAIN_ACTION_INSTRUMENT_SHA256=digest(Path(__file__)),MOMENTUM_DOMAIN_ACTION_PLAN_SHA256=digest(PLAN))
    data={'r':r,'bound':domains[0],'domains':domains,'joins':joins,'finest':finest,'provenance':provenance,'sourceFiles':pins,'smoke':False}
    atomic_pickle(base/'bound-momentum-domains.pickle',{'domains':domains,'joins':joins,'references':finest,'provenance':provenance,'sourceFiles':pins,'dimensionState':dict(vars(dimensions))})
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'cutoffs':CUTOFFS,'rows':80,'sources':70,'profiles':6,'retainedRuleReconstruction':True,'taskWaves':[[(0,1),(1,1)],[(0,2),(0,3),(1,2),(1,3)]],'methodJoin':METHOD_JOIN})
    return data


def task_adapter():
    original=ast.parse(textwrap.dedent(inspect.getsource(position.evaluate_task)));tree=copy.deepcopy(original)
    nodes=[n for n in ast.walk(tree) if isinstance(n,ast.Compare) and ast.unparse(n)== 'index in (1, 2)']
    require(len(nodes)==1,'native task domain census');nodes[0].comparators[0]=ast.parse('(1,2,3)',mode='eval').body
    calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and ast.unparse(n.func)=='settings']
    require(len(calls)==1 and ast.unparse(calls[0].args[1])=='CHOICES[index]','native settings census');calls[0].args[1]=ast.Name(id='index',ctx=ast.Load())
    reverse=copy.deepcopy(tree)
    next(n for n in ast.walk(reverse) if isinstance(n,ast.Compare) and ast.unparse(n)=='index in (1, 2, 3)').comparators[0]=ast.parse('(1,2)',mode='eval').body
    next(n for n in ast.walk(reverse) if isinstance(n,ast.Call) and ast.unparse(n.func)=='settings').args[1]=ast.parse('CHOICES[index]',mode='eval').body
    require(ast.dump(reverse)==ast.dump(original),'entire native worker AST reverse join')
    namespace=dict(vars(position),settings=settings)
    exec(compile(ast.fix_missing_locations(tree),'<momentum-domain-native-worker>','exec'),namespace)
    return namespace['evaluate_task'],{'nativeWorkerAstSha256':hashlib.sha256(ast.dump(original).encode()).hexdigest(),'restoredWorkerAstSha256':hashlib.sha256(ast.dump(reverse).encode()).hexdigest()}

evaluate_task,METHOD_JOIN=task_adapter()


def dispatch(base,data):
    old=parallel.evaluate_task;parallel.evaluate_task=evaluate_task;packets={};manifests=[]
    try:
        for wave,tasks in enumerate(([(0,1),(1,1)],[(0,2),(0,3),(1,2),(1,3)])):
            packets.update(parallel.dispatch(base,data,tasks));manifest=json.loads((base/'workers.json').read_text());save(base/f'workers-wave-{wave}.json',manifest);manifests.append(manifest)
    finally:parallel.evaluate_task=old
    save(base/'workers.json',{'status':'completed','waves':manifests,'outcomes':[v for m in manifests for v in m['outcomes']],'resultFiles':{k:h for m in manifests for k,h in m['resultFiles'].items()},'maximumConcurrentWorkers':4})
    return packets


def combine(base,data,workers):
    records=[];inventory=[];artifacts=[]
    for ti in range(2):
        ref=data['finest'][ti];previous=None
        for index in range(4):
            if index==0:item={'test':ti,'index':0,'method':'retained-profile256','comparisonKind':'retained','groups':ref['groups'],'settings':settings(ref,0),**assemble(data['r'],data['domains'][0],ref,ref['groups'])}
            else:
                w=workers[(ti,index)];directory=base/f'worker-{ti}-{index}'
                require(w['task']==(ti,index) and w['settings']==settings(ref,index,data['smoke']) and w['sourceFiles']==data['sourceFiles'] and w['provenance']==data['provenance'],'worker assignment/rule/provenance')
                artifacts.append(artifact(base,directory/'result.pickle'))
                for a in w['partialArtifacts']:
                    require(digest(directory/a['path'])==a['sha256'],'worker layout/partial hash');artifacts.append(dict(a,path=str((directory/a['path']).relative_to(base))))
                item=dict(w,test=ti,index=index,method='gauss',comparisonKind='profile-rule-change' if index==1 else 'momentum-box-change',precedingActionDifference=w['action']-previous['action'])
                item['integralChanges']=[b['values']-a['values'] for a,b in zip(previous['groups'],w['groups'])];item['termChanges']=[np.asarray(b['values'])-np.asarray(a['values']) for a,b in zip(previous['terms'],w['terms'])]
                for a,b in zip(previous['terms'],w['terms']):require(a['index']==b['index'] and a['rows']==b['rows'],'native term joins')
                if not data['smoke']:
                    for n,s in w['settings'].items():
                        old=previous['settings'][n];changed={k for k in set(s)|set(old) if s.get(k)!=old.get(k)}
                        require(changed==({'profileNodes'} if index==1 else {'momentumBound'}),'one independent changed setting')
            path=base/'records'/f'{ti}-{index}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,item);inventory.append(artifact(base,path));save(base/'record-inventory.json',inventory);records.append(item);previous=item
    result={'records':records,'recordArtifacts':inventory,'workerArtifacts':artifacts,'joins':data['joins']}
    atomic_pickle(base/'momentum-domain-action.pickle',{'result':result,'provenance':data['provenance'],'sourceFiles':data['sourceFiles'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))});return result


def emit(result,r,bound,rows,finest,provenance):
    # Reuse unchanged action/term/cache emission; the domain operand family is
    # momentum limits here, so emit it separately from position-domain joins.
    fn=types.FunctionType(position.emit.__code__,dict(vars(position),PREFIX=PREFIX),position.emit.__name__,position.emit.__defaults__,position.emit.__closure__)
    fn(dict(result,joins={}),r,bound,rows,finest,provenance)
    d=engine.PHYSICAL_METADATA.dimensions;metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=r
    def numeric(name,value,unit):
        body=momentum.number(value);engine.emit(PREFIX+'_'+name,body);engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    for index,j in result['joins'].items():
        numeric('BOUND_'+str(index),j['cutoff'],lambda p:bound['momentumUnit'])
        for row in j['rows']:
            for name in ('old','new','restored'):numeric('LIMIT_'+name.upper()+'_'+str(index)+'_'+str(row['row']),row[name],lambda p:bound['momentumUnit'])
    for item in result['records']:engine.physical(PREFIX+'_COMPARISON_KIND_'+str(item['test'])+'_'+str(item['index']),item['comparisonKind'])


def replay_adapter():
    cache=types.FunctionType(position.cache_for_tag.__code__,dict(vars(position),PREFIX=PREFIX),position.cache_for_tag.__name__)
    factory=types.FunctionType(position.replay_adapter.__code__,dict(vars(position),PREFIX=PREFIX,emit=emit,cache_for_tag=cache),position.replay_adapter.__name__)
    return factory()


def finish(base,data,workers,started):
    result=combine(base,data,workers);before={p.name:digest(p) for p in base.glob('*.pickle')};replay,join=replay_adapter();engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder()
    entries,keys,paths=replay(base,result,data['r'],data['bound'],data['bound']['rows'],data['finest'],data['provenance'])
    norms=[{'test':v['test'],'index':v['index'],'comparisonKind':v['comparisonKind'],'rawIntegralChanges':[float(np.max(abs(a))) for a in v['integralChanges']],'actionChange':float(np.max(abs(v['precedingActionDifference']))),'rawTermChange':max(float(np.max(abs(t))) if len(t) else 0. for t in v['termChanges'])} for v in result['records'] if 'integralChanges' in v]
    summary={'runDirectory':str(base),'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'smoke':data['smoke'],'rows':80,'sources':70,'profiles':6,'records':len(result['records']),'norms':norms,'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,'replayMethodJoin':join,'workerMethodJoin':METHOD_JOIN,'recordArtifacts':result['recordArtifacts'],'workerArtifacts':result['workerArtifacts'],'workerManifest':json.loads((base/'workers.json').read_text()),'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},'wallSeconds':time.monotonic()-started,'coordinatorPeakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'Finite common-rule momentum boxes with explicit profile-rule baseline. Wider-box momentum refinement, infinite tails, uniform/independent-grade coverage, Abel limits, scattering and poles are not established.'}
    for a in (*result['recordArtifacts'],*result['workerArtifacts']):require(digest(base/a['path'])==a['sha256'],'final record/layout/partial identity')
    require(before==summary['packetHashesAfterEmission'] and data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']} and not engine.PHYSICAL_METADATA.dimensions.constraints,'final packet/source/dimension guard')
    save(base/'checks.json',summary);return summary


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);a=parser.parse_args();base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=load(base);workers=dispatch(base,data);summary=finish(base,data,workers,started);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
