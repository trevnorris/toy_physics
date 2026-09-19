#!/usr/bin/env python3
"""A small source-derived finite-pencil contour diagnostic."""
import argparse,ast,copy,inspect,json,os,resource,shutil,signal,subprocess,sys,textwrap,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_frequency_matrix as m
end=m.end;f=m.f
PLAN=f.M/'S11c_d_frequency_contour_plan.md';CENTER=1.-.01j;RADIUS=.02;COUNT=16;WORKERS=4


def tree(fn):return ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
def compile_function(node,namespace):
    scope=dict(namespace);exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),__file__,'exec'),scope);return scope[node.name]


def adapters():
    # Keep the accepted corrector, but assign the literal target at the last
    # point, avoiding a one-ulp multiply/divide round trip in its address.
    original=tree(end.continue_pair);node=copy.deepcopy(original);node.name='continue_to'
    assignment=next(v for v in ast.walk(node) if isinstance(v,ast.Assign) and len(v.targets)==1 and isinstance(v.targets[0],ast.Name) and v.targets[0].id=='w')
    value=copy.deepcopy(assignment.value);assignment.value=ast.IfExp(ast.Compare(ast.Name('j',ast.Load()),[ast.Eq()],[ast.Name('count',ast.Load())]),ast.Name('target',ast.Load()),assignment.value)
    reverse=copy.deepcopy(node);reverse.name=original.name;next(v for v in ast.walk(reverse) if isinstance(v,ast.Assign) and len(v.targets)==1 and isinstance(v.targets[0],ast.Name) and v.targets[0].id=='w').value=value
    f.require(ast.dump(reverse)==ast.dump(original),'whole endpoint-address adapter reverse join');continued=compile_function(node,vars(end))

    original_maps=tree(m.end_maps);node=copy.deepcopy(original_maps);node.name='continued_maps'
    select=next(v for v in ast.walk(node) if isinstance(v,ast.If) and ast.unparse(v.test)=='w == 1.0');prior=copy.deepcopy(select.orelse)
    select.orelse=ast.parse("directory=base/(label.lower()+f'-end-{seed[\"index\"]}')\ndirectory.mkdir()\npair=end.Pair(data['chart']['ends'][label],seed)\nstate,paths=continue_to(directory,pair,seed['seedState'],w,.01)\norigin={'method':'continued from fixed seed','points':paths}").body
    reverse=copy.deepcopy(node);reverse.name=original_maps.name;next(v for v in ast.walk(reverse) if isinstance(v,ast.If) and ast.unparse(v.test)=='w == 1.0').orelse=prior
    f.require(ast.dump(reverse)==ast.dump(original_maps),'whole end-map selector reverse join');maps=compile_function(node,dict(vars(m),continue_to=continued))

    original_case=tree(m.complex_case);node=copy.deepcopy(original_case);node.name='frequency_case';frequency_assignment=node.body.pop(0)
    f.require(isinstance(frequency_assignment,ast.Assign) and ast.unparse(frequency_assignment.targets[0])=='frequency','native fixed pilot frequency assignment')
    node.args.args.append(ast.arg(arg='frequency'))
    reverse=copy.deepcopy(node);reverse.name=original_case.name;reverse.args.args.pop();reverse.body.insert(0,frequency_assignment)
    f.require(ast.dump(reverse)==ast.dump(original_case),'whole frequency-case reverse join');case=compile_function(node,dict(vars(m),end_maps=maps))
    return case,continued,{'wholeEndpointAddressJoin':True,'wholeEndMapSelectorJoin':True,'wholeFrequencyCaseJoin':True,'unchangedNativeAccumulator':True,'unchangedCorrectorAndFrechetEquations':True}


def frequencies():return [CENTER+RADIUS*np.exp(2j*np.pi*i/COUNT) for i in range(COUNT)]


def load(base):
    data=m.load(base);p=f.M/'S11c_d_frequency_matrix_checkpoint.json';cp=json.loads(p.read_text());f.require(cp['status']=='VALIDATED_FIRST_FINITE_FREQUENCY_MATRIX','accepted actual finite frequency pencil')
    for name,h in cp['sourceFiles'].items():f.require(data['pins'][name]==h==f.digest(f.ROOT/name),'unchanged actual finite frequency constructor')
    accepted=Path(cp['runDirectory']);seed=accepted/'seed/frequency-system.pickle';f.require(f.digest(seed)==cp['artifacts']['seed/frequency-system.pickle']['sha256'],'actual fixed seed scaling packet')
    system=f.unpickle(seed);data['fixedScales']=(system['fixedRowScale'],system['fixedColumnScale']);data['operands'][str(seed)]=f.digest(seed)
    for path in (Path(__file__),PLAN,p):data['pins'][str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    for name,h in data['pins'].items():
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,target)
    data['manifest'].update(sourceFiles=data['pins'],inputPackets=data['operands'],scope='Sixteen actual finite-pencil samples on omega=1-0.01i+0.02 exp(i theta), nested8/16 winding and inverse/response moments. Numerical contour evidence, not certified spectral completeness or a physical bound-pole set.')
    f.require(abs(CENTER-1)+RADIUS<.25,'entire contour inside accepted real-momentum scalar chart')
    f.save(base/'inputs.json',data['manifest']);return data


def winding(signs):
    signs=np.asarray(signs,complex);increments=np.angle(np.roll(signs,-1)*np.conj(signs))
    return {'increments':increments,'winding':float(np.sum(increments)/(2*np.pi)),'nearestInteger':int(np.rint(np.sum(increments)/(2*np.pi))),'maximumPhaseIncrement':float(np.max(abs(increments)))}


def controls(base):
    records={}
    for n in (8,16):
        theta=2*np.pi*np.arange(n)/n;z=RADIUS*np.exp(1j*theta);weights=z/n
        cases={'simple':(z,1),'semisimpleDouble':(z*z,2),'higherOrderDouble':(z*z,2),'twoDistinct':((z-.004)*(z+.006),2)}
        for name,(det,count) in cases.items():
            w=winding(det/abs(det));f.require(w['nearestInteger']==count,'actual synthetic contour multiplicity control');records[name,n]=w
        inverse=1/z**2;mom0=np.sum(weights*inverse);mom1=np.sum(weights*z*inverse)
        f.require(abs(mom0)<1e-12 and abs(mom1-1)<1e-12,'higher-order inverse survives vanishing residue')
        records['higherOrderMoments',n]={'residueMoment':mom0,'nextMoment':mom1}
        wrong=winding((z/abs(z))[::-1]);f.require(wrong['nearestInteger']==-1,'actual contour orientation control');records['wrongOrientation',n]=wrong
    # An explicit alias control prevents promoting matching sampled counts
    # into a certified argument-principle count without between-node evidence.
    aliases={}
    for n in (8,16):
        z=np.exp(2j*np.pi*np.arange(n)/n);aliases[n]=winding(z**16)
    f.require(all(v['nearestInteger']==0 for v in aliases.values()),'known degree16 alias demonstration')
    records['aliasLimitation']={'actualMultiplicity':16,'sampled':aliases,'conclusion':'Agreement of nested sampled winding alone is not a certified zero count.'}
    f.atomic_pickle(base/'contour-controls.pickle',records);return records


def worker(base,index):
    start=time.monotonic();data=load(base);case,_,joins=adapters();f.save(base/'adapter-joins.json',joins);records=[]
    for i in range(index,COUNT,WORKERS):
        w=frequencies()[i];directory=base/f'point-{i:02d}';directory.mkdir()
        frequency=sp.Float(w.real,17)+sp.I*sp.Float(w.imag,17)
        system,solution,groups=case(directory,data,data['fixedScales'],frequency)
        inverse=np.linalg.solve(system['balanced'],np.eye(645));right=system['balanced']@inverse-np.eye(645);left=inverse@system['balanced']-np.eye(645)
        f.require(end.norm(left)<1e-8 and end.norm(right)<1e-8,'actual full inverse identities in fixed frames')
        sign,logabs=np.linalg.slogdet(system['balanced']);f.require(abs(abs(sign)-1)<1e-10 and np.isfinite(logabs),'finite nonzero sampled determinant')
        result={'index':i,'frequency':w,'balancedInverse':inverse,'leftInverseResidual':left,'rightInverseResidual':right,'determinantPhase':complex(sign),'logAbsDeterminant':float(logabs),'openResponse':solution['openOriginScattering'],'condition':solution['fixedFrameCondition'],'smallestSingularValue':float(system['singularValues'][-1]),'rank':system['rank'],'scaledEquationResidual':end.norm(solution['scaledEquationResidual']),'newMomentumNodes':sum(v['nodes'] for v in groups),'fixedRowScale':system['fixedRowScale'],'fixedColumnScale':system['fixedColumnScale'],'sourceRowUnits':system['sourceRowUnits'],'fieldUnits':system['fieldUnits'],'scope':'Full inverse in fixed coefficient frames; convert with recorded diagonal maps for physical inverse entries.'}
        f.atomic_pickle(directory/'contour-point.pickle',result)
        record={'index':i,'frequency':str(w),'directory':str(directory),'sha256':f.digest(directory/'contour-point.pickle'),'artifacts':{str(p.relative_to(directory)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in directory.rglob('*.pickle')},'wallSeconds':sum(g['wallSeconds'] for g in groups),'newMomentumNodes':result['newMomentumNodes']}
        records.append(record);f.save(base/'point-inventory.json',records)
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'worker source/input post hashes')
    result={'workerIndex':index,'records':records,'sourceFiles':data['pins'],'inputPackets':data['operands'],'joins':joins,'threadEnvironment':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')},'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.require(set(result['threadEnvironment'].values())=={'1'},'one native thread per worker');f.save(base/'worker-checks.json',result);return result


def close_end_paths(base,data,inventory,continued):
    frequencies_=frequencies();records=[]
    for label,selected in data['endSeeds'].items():
        for seed in selected['clusters']:
            pair=end.Pair(data['chart']['ends'][label],seed);state=seed['seedState'];initial=None;path_records=[]
            for step,index in enumerate(list(range(COUNT))+[0]):
                directory=base/(label.lower()+f'-cluster-{seed["index"]}-step-{step:02d}');directory.mkdir()
                state,paths=continued(directory,pair,state,frequencies_[index],.01)
                if initial is None:initial=state
                native_maps=f.unpickle(Path(inventory[index]['directory'])/'frequency-end-maps.pickle')[label]
                radial=next(v['state'] for v in native_maps['clusters'] if v['seed']['index']==seed['index'])
                delta=(state['x']-radial['x'])/pair.unknown_scale
                f.require(end.norm(delta)<1e-8,'contour-versus-radial full end continuation')
                path_records.append({'index':index,'points':paths,'radialDifference':delta})
            closure=(state['x']-initial['x'])/pair.unknown_scale;f.require(end.norm(closure)<1e-8,'complete numerical end loop closure')
            record={'label':label,'index':seed['index'],'paths':path_records,'closureDifference':closure};f.atomic_pickle(base/(label.lower()+f'-closure-{seed["index"]}.pickle'),record);records.append(record)
    f.atomic_pickle(base/'end-loop-closure.pickle',records);return records


def summarize(base,data,inventory):
    rows=sorted(inventory,key=lambda v:v['index']);f.require([v['index'] for v in rows]==list(range(COUNT)),'all actual contour nodes')
    results={};metadata=[]
    for n in (8,16):
        selected=rows[::COUNT//n];signs=[];moments=[np.zeros((645,645),complex) for _ in range(4)];responses=[np.zeros((4,4),complex) for _ in range(4)]
        for record in selected:
            p=Path(record['directory'])/'contour-point.pickle';f.require(f.digest(p)==record['sha256'],'actual completed contour packet');x=f.unpickle(p);i=x['index'];w=frequencies()[i];f.require(x['frequency']==w,'literal contour-node address')
            delta=w-CENTER;weight=RADIUS*np.exp(2j*np.pi*i/COUNT)/n;signs.append(x['determinantPhase'])
            for order in range(4):moments[order]+=weight*delta**order*x['balancedInverse'];responses[order]+=weight*delta**order*x['openResponse']
            if n==16:metadata.append({k:x[k] for k in ('index','frequency','condition','smallestSingularValue','rank','scaledEquationResidual','logAbsDeterminant','determinantPhase','newMomentumNodes')})
        result={'count':n,'winding':winding(signs),'inverseMoments':moments,'responseMoments':responses,'scope':'Numerical contour moments; not assumed projectors, residues of a simple pole, or a certified multiplicity count.'}
        f.atomic_pickle(base/f'contour-{n}.pickle',result);results[n]=result
    comparison={'windingIntegerDifference':results[16]['winding']['nearestInteger']-results[8]['winding']['nearestInteger'],'inverseMomentDifferences':[a-b for a,b in zip(results[16]['inverseMoments'],results[8]['inverseMoments'])],'responseMomentDifferences':[a-b for a,b in zip(results[16]['responseMoments'],results[8]['responseMoments'])]}
    f.atomic_pickle(base/'contour-comparison.pickle',comparison);f.atomic_pickle(base/'point-summary.pickle',metadata);return results,comparison,metadata


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);ap.add_argument('--worker',type=int);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));children=[]
    def timeout(*_):raise TimeoutError('bounded contour budget; preserve completed points and partials')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);start=time.monotonic()
    if args.worker is not None:
        result=worker(base,args.worker);signal.alarm(0);print(json.dumps(result,indent=2));return
    data=load(base);_,continued,joins=adapters();controls(base);f.save(base/'preflight.json',{'center':str(CENTER),'radius':RADIUS,'samples':COUNT,'workers':WORKERS,'joins':joins,'scope':data['manifest']['scope']})
    inventory=[];outcomes=[];logs=[]
    try:
        for i in range(WORKERS):
            directory=base/f'worker-{i}';out=(base/f'worker-{i}.stdout').open('xb');err=(base/f'worker-{i}.stderr').open('xb');logs.extend((out,err));command=[sys.executable,str(Path(__file__).resolve()),'--run-directory',str(directory),'--worker',str(i)]
            child=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=out,stderr=err);children.append((i,child,directory))
        for i,child,directory in children:
            code=child.wait();outcomes.append({'index':i,'pid':child.pid,'exitCode':code,'stderrBytes':(base/f'worker-{i}.stderr').stat().st_size});f.save(base/'workers.json',outcomes)
        f.require(all(v['exitCode']==0 and v['stderrBytes']==0 for v in outcomes),'all four contour workers clean')
        for i,child,directory in children:
            checks=json.loads((directory/'worker-checks.json').read_text());f.require(checks==json.loads((base/f'worker-{i}.stdout').read_text()),'worker checks/stdout identity')
            f.require(checks['sourceFiles']==data['pins'] and checks['inputPackets']==data['operands'],'coordinator/worker source/input identity')
            for record in checks['records']:
                for n,v in record['artifacts'].items():f.require(f.digest(Path(record['directory'])/n)==v['sha256'],'every completed point artifact hash')
            inventory.extend(checks['records'])
        inventory.sort(key=lambda v:v['index']);f.save(base/'point-inventory.json',inventory)
        closure_dir=base/'end-loops';closure_dir.mkdir();closures=close_end_paths(closure_dir,data,inventory,continued)
        results,comparison,metadata=summarize(base,data,inventory)
        f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'coordinator source/input post hashes')
        result={'status':'COMPLETED_NUMERICAL_CONTOUR_DIAGNOSTIC','sourceFiles':data['pins'],'inputPackets':data['operands'],'workers':outcomes,'contourCenter':str(CENTER),'contourRadius':RADIUS,'points':COUNT,'winding':{str(n):{k:v for k,v in results[n]['winding'].items() if k!='increments'} for n in results},'windingIntegerDifference':comparison['windingIntegerDifference'],'inverseMomentChangeNorms':[end.norm(v) for v in comparison['inverseMomentDifferences']],'responseMomentChangeNorms':[end.norm(v) for v in comparison['responseMomentDifferences']],'maximumPointCondition':max(v['condition'] for v in metadata),'minimumPointSingularValue':min(v['smallestSingularValue'] for v in metadata),'maximumEndClosureDifference':max(end.norm(v['closureDifference']) for v in closures),'newMomentumNodes':sum(v['newMomentumNodes'] for v in metadata),'artifacts':{str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-start,'peakCoordinatorRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':data['manifest']['scope'],'limitations':['Matching sampled winding is not a certified count; the saved degree16 alias control demonstrates this.','End loop agreement tests the sampled branch; it does not certify all interior exceptional loci.','No physical pole set, empty spectrum or principal-part classification is supplied by this initial diagnostic.']}
        f.save(base/'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))
    finally:
        for _,child,_ in children:
            if child.poll() is None:child.terminate()
        for _,child,_ in children:
            if child.poll() is None:
                try:child.wait(timeout=5)
                except subprocess.TimeoutExpired:child.kill();child.wait()
        for handle in logs:handle.close()


if __name__=='__main__':main()
