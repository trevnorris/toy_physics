#!/usr/bin/env python3
"""Isolated whole-grid workers, exact serial-order resume and deterministic replay."""
import argparse,ast,contextlib,copy,hashlib,inspect,json,multiprocessing as mp,os,resource,sys,textwrap,time,traceback
from multiprocessing.connection import wait
from pathlib import Path
import numpy as np
import S11c_d_three_momentum_source_check as native
from S11c_d_three_momentum_source_check import engine,ROOT,STORE,digest,save,atomic_pickle,unpickle

PLAN=ROOT/'_measurements/S11c_d_parallel_momentum_plan.md'
PACKETS=('accepted-bound-momentum.pickle','accepted-momentum-action.pickle','accepted-single-momentum.pickle','accepted-pair-momentum.pickle','accepted-three-momentum.pickle')


def require(ok,detail):
    if not ok:raise ValueError(detail)


def generated_methods():
    original=ast.parse(textwrap.dedent(inspect.getsource(engine.BoundedSourceFourierQuadrature.FiniteMomentum.group)))
    modified=copy.deepcopy(original);body=modified.body[0].body
    index=next(i for i,n in enumerate(body) if isinstance(n,ast.For))
    injection=ast.parse('if self.resume_state is not None:\n values,mutated,mass,node_count,batch_count=self.seed(test,variables,setting,width,positions,rows)').body[0]
    body.insert(index,injection);proof=copy.deepcopy(modified);proof.body[0].body.pop(index)
    require(ast.dump(proof)==ast.dump(original),'resumed kernel AST join')
    namespace=dict(vars(engine));exec(compile(ast.fix_missing_locations(modified),'<native-group-with-resume>','exec'),namespace)
    source=ast.parse(textwrap.dedent(inspect.getsource(engine.BoundedSourceFourierQuadrature.FiniteMomentum.source_value)))
    trace=copy.deepcopy(source);trace.body[0].name='trace_frequency';statements=trace.body[0].body
    stop=next(i for i,n in enumerate(statements) if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Tuple))
    # The first tuple assignment is symbols,function; retain through census update.
    stop=next(i for i,n in enumerate(statements) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='unique' for t in ast.walk(n.targets[0])))
    trace.body[0].body=statements[:stop]
    require(ast.dump(ast.Module(body=trace.body[0].body,type_ignores=[]))==ast.dump(ast.Module(body=source.body[0].body[:stop],type_ignores=[])),'frequency trace AST join')
    exec(compile(ast.fix_missing_locations(trace),'<native-frequency-prefix>','exec'),namespace)
    return namespace['group'],namespace['trace_frequency'],{'nativeGroupAstSha256':hashlib.sha256(ast.dump(original).encode()).hexdigest(),'nativeFrequencyPrefixAstSha256':hashlib.sha256(ast.dump(ast.Module(body=source.body[0].body[:stop],type_ignores=[])).encode()).hexdigest()}


GROUP,TRACE,METHOD_JOINS=generated_methods()


class ResumableMomentum(engine.BoundedSourceFourierQuadrature.ThreeMomentum):
    group=GROUP
    trace_frequency=TRACE

    def __init__(self,*args,resume_state=None,**kwargs):
        super().__init__(*args,**kwargs);self.resume_state=resume_state;self.resume_prefix_nodes=0;self.active_test=None

    def seed(self,test,variables,setting,width,positions,rows):
        s=self.resume_state;self.active_test=test
        require(s['test']==test and tuple(s['variables'])==tuple(variables) and s['setting']==setting and s['width']==width,'resume domain join')
        require(s['rowIndices']==[r['index'] for r in rows] and np.array_equal(s['positions'],positions),'resume rows/positions')
        require(s['values'].shape==(len(rows),len(positions)) and np.all(np.isfinite(s['values'])) and np.all(np.isfinite(s['measureMutationValues'])),'resume arrays')
        return s['values'].copy(),s['measureMutationValues'].copy(),s['quadratureMass'],s['nodeCount'],s['batchCount']

    def batches(self,variables,setting,pairs,width):
        cutoff=0 if self.resume_state is None else self.resume_state['batchCount']
        for index,(points,weights) in enumerate(super().batches(variables,setting,pairs,width),1):
            if index<=cutoff:
                environment={v:points[:,i] for i,v in enumerate(variables)}
                environment[self.r.regulator]=np.full(len(weights),setting['regulator'])
                rows=[r for r in self.rows if tuple(l[0] for l in r['limits'])==tuple(variables)]
                if index==1:
                    # Actual evaluations initialize observed profile/cache state, even
                    # for an already complete grid. Their census is not counted twice.
                    census=dict(self.source_frequency_census);profiles={};sources={}
                    for row in rows:
                        for factor in row['factors']:
                            self.coefficient_value(factor['coefficient'],environment,np.asarray(self.resume_state['positions']),setting,profiles)
                            self.source_value(self.sources[(self.active_test,factor['sourceIndex'])],environment,setting,sources)
                    self.source_frequency_census=census
                for row in rows:
                    for factor in row['factors']:
                        self.trace_frequency(self.sources[(self.active_test,factor['sourceIndex'])],environment,setting,{})
                self.resume_prefix_nodes+=len(weights)
                continue
            if cutoff:require(self.resume_prefix_nodes==self.resume_state['nodeCount'],'resume generated node count')
            yield points,weights
        if cutoff:require(self.resume_prefix_nodes==self.resume_state['nodeCount'],'complete resume generated node count')


class PrefixComplete(Exception):pass


def cache_state(worker):
    caches=[]
    for (points,order),(nodes,weights) in worker.fixed_rules.items():
        a,b=engine.BoundedSourceFourierQuadrature.rule(points,order)
        caches.append({'points':points,'order':order,'nodesResidual':nodes-a,'weightsResidual':weights-b,'readOnly':not nodes.flags.writeable and not weights.flags.writeable,'bytes':nodes.nbytes+weights.nbytes,'profileRule':(points,order) in worker.profile_rule_keys})
    require(all(not v['nodesResidual'].any() and not v['weightsResidual'].any() and v['readOnly'] for v in caches),'worker cache joins')
    return {'cacheChecks':caches,'sourceFrequencyCensus':worker.source_frequency_census,'evaluatedProfileIntegrals':worker.profile_integrals,'phaseWorkspaceEstimateBytes':worker.peak_phase_workspace_estimate,'batchCacheEstimateBytes':worker.peak_batch_cache_estimate,'workspaceBudgetBytes':worker.workspace_bytes}


def evaluate_task(directory,data,task,*,max_batches=None,resume_state=None):
    r,bound,rows,variables,finest,provenance,pins=data
    ti,index=task;setting=dict(finest[ti]['setting'],sourceNodes=256,profileNodes=128 if index==1 else 256)
    width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
    worker=ResumableMomentum(bound['rows'],bound['sources'],r,resume_state=resume_state)
    artifacts=[];last=[];started=time.monotonic()
    def checkpoint(state):
        if state['batchCount']%64==0 or (max_batches and state['batchCount']==max_batches):
            p=directory/'partials'/f'{state["batchCount"]:08}.pickle';p.parent.mkdir(exist_ok=True);atomic_pickle(p,state)
            artifacts.append({'path':str(p.relative_to(directory)),'bytes':p.stat().st_size,'sha256':digest(p)})
            save(directory/'partial-inventory.json',artifacts)
        if max_batches is not None and state['batchCount']==max_batches:last.append(state);raise PrefixComplete()
    if resume_state is not None:
        atomic_pickle(directory/'resume-state.pickle',resume_state)
        # Preserve the imported prefix as its own operand. No tail is fabricated.
        artifacts.append({'path':'resume-state.pickle','bytes':(directory/'resume-state.pickle').stat().st_size,'sha256':digest(directory/'resume-state.pickle')})
    try:
        group=worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'],checkpoint)
        kind='complete'
    except PrefixComplete:
        group=last[0];kind='prefix'
    if max_batches is not None:require(kind=='prefix','bounded worker exceeded test scope')
    else:
        require(abs(group['volumeResidual'])<1e-10*(1+abs(group['boxVolume'])),'worker finite mass')
        require(not np.any((np.abs(group['values'])>1e-9)&(np.abs(group['measureMutationResidual'])<1e-12)),'worker measure response')
    telemetry=cache_state(worker)
    prefix_workspace=0 if resume_state is None else resume_state['peakWorkspaceEstimateBytes']
    # A saved legacy prefix has a combined workspace estimate, not its split.
    # Use that observed combined value conservatively for each split envelope.
    telemetry['phaseWorkspaceEstimateBytes']=max(prefix_workspace,telemetry['phaseWorkspaceEstimateBytes'])
    telemetry['batchCacheEstimateBytes']=max(prefix_workspace,telemetry['batchCacheEstimateBytes'])
    expected={(ti,f['sourceIndex']) for row in rows for f in row['factors']}
    require(set(telemetry['sourceFrequencyCensus'])==expected,'worker source census')
    profiles=set().union(*(f['coefficient'].atoms(native.sp.Integral) for row in rows for f in row['factors']))
    require(telemetry['evaluatedProfileIntegrals']==profiles,'worker profile census')
    packet={'task':task,'kind':kind,'setting':setting,'group':group,'telemetry':telemetry,'partialArtifacts':artifacts,'provenance':provenance,'sourceFiles':pins,
        'resumePrefixNodes':worker.resume_prefix_nodes,'resumePrefixCombinedWorkspaceEstimateBytes':prefix_workspace,
        'resumeFrequencyScope':'Literal native frequency-prefix replay on every saved-prefix node; numerical prefix sums retained. Workspace envelopes include the saved combined prefix estimate.',
        'methodJoins':METHOD_JOINS,'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    atomic_pickle(directory/'result.pickle',packet)
    require(pins=={n:digest(ROOT/n) for n in pins},'worker source changed')
    save(directory/'checks.json',{'task':task,'kind':kind,'resultSha256':digest(directory/'result.pickle'),'wallSeconds':packet['wallSeconds'],'peakRssKiB':packet['peakRssKiB'],'resumePrefixNodes':worker.resume_prefix_nodes,'pid':os.getpid()})
    return packet


def child(directory,data,task,max_batches,resume_state):
    with (directory/'stdout').open('x') as out,(directory/'stderr').open('x') as err,contextlib.redirect_stdout(out),contextlib.redirect_stderr(err):
        try:
            resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
            evaluate_task(directory,data,task,max_batches=max_batches,resume_state=resume_state)
        except BaseException:traceback.print_exc();sys.exit(1)


def dispatch(base,data,tasks,resumes=None,*,max_batches=None):
    require(all(os.environ.get(k)=='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'one native thread per worker required')
    require(1<=len(tasks)<=4 and len(set(tasks))==len(tasks),'worker task census')
    ctx=mp.get_context('fork');workers={};started=time.monotonic();records=[]
    for task in tasks:
        d=base/f'worker-{task[0]}-{task[1]}';d.mkdir();p=ctx.Process(target=child,args=(d,data,task,max_batches,(resumes or {}).get(task)));p.start();workers[p.sentinel]=(p,d,task)
    save(base/'workers.json',{'status':'running','workers':[{'task':task,'pid':p.pid,'directory':str(d)} for p,d,task in workers.values()]})
    while workers:
        for sentinel in wait(list(workers)):
            p,d,task=workers.pop(sentinel);p.join();record={'task':task,'pid':p.pid,'exitCode':p.exitcode,'stderrBytes':(d/'stderr').stat().st_size};records.append(record)
            if p.exitcode or record['stderrBytes']:
                for other,_,_ in workers.values():other.terminate()
                for other,_,_ in workers.values():other.join()
                save(base/'workers.json',{'status':'failed','outcomes':records})
                raise RuntimeError(('worker failed',str(d),(d/'stderr').read_text()[-3000:]))
    packets={}
    for task in sorted(tasks):
        d=base/f'worker-{task[0]}-{task[1]}';c=json.loads((d/'checks.json').read_text());require(digest(d/'result.pickle')==c['resultSha256'],'worker final packet hash')
        packets[task]=unpickle(d/'result.pickle')
    save(base/'workers.json',{'status':'completed','outcomes':records,'wallSeconds':time.monotonic()-started,'resultFiles':{f'{t[0]}-{t[1]}':digest(base/f'worker-{t[0]}-{t[1]}'/'result.pickle') for t in sorted(tasks)}})
    return packets


def construct_from_workers(base,data,packets,constructor=None):
    r,bound,rows,variables,finest,provenance,pins=data
    parent=engine.BoundedSourceFourierQuadrature.ThreeMomentum
    class Replay(parent):
        def group(self,test,variables,setting,pairs,width,positions,batch_checkpoint=None):
            index=1 if setting['profileNodes']==128 else 2;packet=packets[(test,index)]
            require(packet['kind']=='complete' and packet['provenance']==provenance and packet['sourceFiles']==pins,'replay worker join')
            group=packet['group'];require(tuple(group['variables'])==tuple(variables) and packet['setting']==setting and group['test']==test and group['rowIndices']==[row['index'] for row in rows],'replay native limits/settings')
            d=base/f'worker-{test}-{index}'
            for record in packet['partialArtifacts']:
                p=d/record['path'];require(digest(p)==record['sha256'],'worker partial hash')
                state=unpickle(p);require(state['setting']==setting and state['test']==test and state['width']==width and np.array_equal(state['positions'],positions),'worker setting join')
                if batch_checkpoint:batch_checkpoint(state)
            for cache in packet['telemetry']['cacheChecks']:
                self.fixed_rule(cache['points'],cache['order'])
                if cache['profileRule']:self.profile_rule_keys.add((tuple(cache['points']),cache['order']))
            for key,(lo,hi,count) in packet['telemetry']['sourceFrequencyCensus'].items():
                old=self.source_frequency_census.get(key,(float('inf'),float('-inf'),0));self.source_frequency_census[key]=(min(old[0],lo),max(old[1],hi),old[2]+count)
            self.profile_integrals.update(packet['telemetry']['evaluatedProfileIntegrals'])
            self.peak_phase_workspace_estimate=max(self.peak_phase_workspace_estimate,packet['telemetry']['phaseWorkspaceEstimateBytes'])
            self.peak_batch_cache_estimate=max(self.peak_batch_cache_estimate,packet['telemetry']['batchCacheEstimateBytes'])
            return group
    engine.BoundedSourceFourierQuadrature.ThreeMomentum=Replay
    try:return (constructor or native.construct)(base,r,bound,rows,variables,finest)
    finally:engine.BoundedSourceFourierQuadrature.ThreeMomentum=parent


def prepare(base):
    data=list(native.load(base));pins=data[-1]
    for p in (Path(__file__).resolve(),PLAN,ROOT/'_measurements/S11c_d_parallel_momentum_preflight.py'):
        n=str(p.relative_to(ROOT));pins[n]=digest(p);target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes())
    data[-2]=dict(data[-2],PARALLEL_INSTRUMENT_SHA256=digest(Path(__file__)),PARALLEL_PLAN_SHA256=digest(PLAN))
    save(base/'parallel-preflight.json',{'sourceFiles':pins,'provenance':data[-2],'methodJoins':METHOD_JOINS})
    return tuple(data)


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--handoff',type=Path,required=True);a=p.parse_args()
    base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=prepare(base)
    handoff=json.loads(a.handoff.read_text());require(handoff['status']=='SERIAL_STOPPED_FOR_PARALLEL','handoff not complete')
    resumes={}
    for key,record in handoff['resumes'].items():
        path=Path(record['path']);require(digest(path)==record['sha256'],'handoff resume hash');resumes[tuple(map(int,key.split('-')))]=unpickle(path)
    for n,h in handoff['sourceFiles'].items():require(digest(ROOT/n)==h,'handoff source changed')
    save(base/'handoff.json',handoff);data=list(data);data[-2]=dict(data[-2],SERIAL_HANDOFF_SHA256=digest(a.handoff));data=tuple(data)
    before={n:digest(base/n) for n in PACKETS}
    packets=dispatch(base,data,[(ti,index) for ti in range(2) for index in (1,2)],resumes)
    result=construct_from_workers(base,data,packets);before['three-momentum-source.pickle']=digest(base/'three-momentum-source.pickle')
    r,bound,rows,variables,finest,provenance,pins=data
    entries,keys,metadata_paths=native.emit_and_replay(base,result,r,bound,rows,finest,provenance)
    norms=[]
    for item in result['records']:
        v={'test':item['test'],'index':item['index'],'method':item['method'],'setting':item['setting']}
        for n in ('referenceIntegralDifference','precedingIntegralDifference','differenceFromReference','precedingActionDifference'):
            if n in item:v[n]=float(np.max(abs(item[n])))
        norms.append(v)
    worker_artifacts=[]
    for task,packet in sorted(packets.items()):
        d=base/f'worker-{task[0]}-{task[1]}'
        worker_artifacts.append({'path':str((d/'result.pickle').relative_to(base)),'bytes':(d/'result.pickle').stat().st_size,'sha256':digest(d/'result.pickle')})
        for a in packet['partialArtifacts']:worker_artifacts.append(dict(a,path=str((d/a['path']).relative_to(base))))
    summary={'runDirectory':str(base),'sourceFiles':pins,'provenance':provenance,'methodJoins':METHOD_JOINS,'numericalEngineUnchanged':True,
        'rowCount':len(rows),'profileCount':len(result['evaluatedProfileIntegrals']),'sourceCount':len(result['sourceFrequencyCensus']),'records':len(result['records']),'norms':norms,'ruleCaches':len(result['cacheChecks']),
        'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':metadata_paths,'recordArtifacts':result['recordArtifacts'],'partialArtifacts':result['partialArtifacts'],'workerArtifacts':worker_artifacts,
        'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},'wallSeconds':time.monotonic()-started,
        'coordinatorPeakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'workerPeakRssKiB':{f'{t[0]}-{t[1]}':v['peakRssKiB'] for t,v in packets.items()},
        'resumedNodes':{f'{t[0]}-{t[1]}':v['resumePrefixNodes'] for t,v in packets.items()},'workerManifest':json.loads((base/'workers.json').read_text()),
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'scope':'Four independent numerical workers; unchanged serial summation order within each integral. Legacy prefix sums are retained, frequency census replayed and workspace split conservatively bounded by the saved combined estimate. Finite source/profile refinement only; no physical tails, Abel, scattering or poles.'}
    save(base/'checks.json',summary)
    for record in (*result['recordArtifacts'],*result['partialArtifacts'],*worker_artifacts):require(digest(base/record['path'])==record['sha256'],'final saved operand hash')
    require(pins=={n:digest(ROOT/n) for n in pins} and before=={n:digest(base/n) for n in before} and not engine.PHYSICAL_METADATA.dimensions.constraints,'final source/packet/dimension guard')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
