#!/usr/bin/env python3
"""Remaining material1D rows with exact saved full-grid operation reuse."""
import argparse
import ast
import copy
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_material_row_pilot as pilot

saved,p,io,sp,np=pilot.saved,pilot.p,pilot.io,pilot.sp,pilot.np
M,F,require,same=pilot.M,pilot.F,pilot.require,pilot.same
NAME='S11c_d_remaining_case_frequency_material_rows_1d'
PILOT_SHA='f8db96abc7cb38c7133786b52431c8939676d8c314376c6b9d63a6fe5fdd9b9b'
CP_SHA='529f4943c8da52dfc153716e4acadc753ca611fe3e0507135cf54cb91aee7e56'
ROWS=(3,4,5,6,7,8,9,10,11,21,22,23,24,25)
COMPARE=(11,25)


def scope(row,settings):
    return {'case':pilot.CASE,'rowIndex':row,'frequency':{'real':1.,'imag':-.01},'finePanels':1024,'coarsePanels':512 if row in COMPARE else None,
        'nativeSettings':settings,'method':'Saved finite trapezoid arrays, native minus-phase source Fourier and single-measure ordinary transpose. Exact full-grid accepted/new values reuse before new calls.',
        'scope':'Missing material1D rows only. Selected comparisons on rows11/25, no automatic extension. Fixed-row spreads are not scattering error bounds, tiny-effect/pole/domain claims or numerical reuse partitions.'}


def frozen_sources(reader,base):
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        frozen=base.parent/'source'/path.name;target=base/'source'/path.name;target.parent.mkdir(exist_ok=True)
        require(not target.exists() and not target.is_symlink(),'fresh exact frozen source reference')
        target.symlink_to(frozen);reader.retain(target,saved.digest(path))


def loader_adapter(text):
    original=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='main')
    begin=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='(reader, journal, cache)')
    end=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='grid_results')
    fn=copy.deepcopy(original);fn.name='load';fn.args=ast.arguments(posonlyargs=[],args=[ast.arg(arg='base')],kwonlyargs=[],kw_defaults=[],defaults=[]);fn.body=copy.deepcopy(original.body[begin:end]);before=copy.deepcopy(fn);changes=[]
    for i,node in enumerate(fn.body):
        new=None
        if i==0:
            new=copy.deepcopy(node);new.value.elts[2]=ast.Name(id='SHARED_PACKET_CACHE',ctx=ast.Load())
        elif isinstance(node,ast.Expr) and isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=='require' and len(node.value.args)>1 and isinstance(node.value.args[1],ast.Constant) and str(node.value.args[1].value).startswith('literal source frequency is'):
            new=copy.deepcopy(node);new.value.args[0]=ast.parse("si==view['sourceIndex'] and same(source['frequency'],variable)",mode='eval').body
        elif isinstance(node,ast.For) and any(isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='dest' for n in node.body):
            new=ast.parse('frozen_sources(reader,base)').body[0]
        elif isinstance(node,ast.Assign) and ast.unparse(node.targets[0])=='scope':new=ast.parse('scope=batch_scope(ROW,settings)').body[0]
        if new is not None:changes.append((i,copy.deepcopy(node)));fn.body[i]=new
    require(len(changes)==4,'exact metadata/cache/source-selection loader changes')
    reverse=copy.deepcopy(fn)
    for i,node in changes:reverse.body[i]=node
    require(ast.dump(reverse)==ast.dump(before),'whole original saved-input prefix reverses; numerical tail excluded')
    fn.body.append(ast.parse('return locals()').body[0])
    return ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))


def physical(inp):
    return {k:inp[k] for k in ('context','settings','fieldUnits','equationUnits','abel','pairs','profileUnits')}


def coefficient_key(inp,momentum):
    return {'expression':inp['coefficient'],'environment':{inp['row']['limits'][0][0]:momentum[:,None],inp['context']['z']:inp['positions'][None,:],inp['context']['regulator']:inp['settings']['regulator']},
        'unit':inp['row']['factors'][0]['unit'],'physical':physical(inp)}


def fourier_key(inp,momentum,nodes,action):
    rec=inp['sourceAction'];path=Path(rec.get('path',rec.get('logical')))
    require(action.shape==(1024,129) and action.dtype==np.dtype(complex),'full saved action shape/type before reference key')
    return {'frequencies':momentum,'frequencyExpression':inp['source']['frequency'],'sourceNodes':nodes,
        'weightedSourceActionIdentity':{'canonical':str(path.resolve()),'sha256':rec['sha256'],'bytes':rec['bytes']},
        'source':inp['source'],'jet':inp['jet'],'unit':inp['jet']['integralUnit'],'physical':physical(inp)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal=saved.Reader(),pilot.prior.inputs.MetadataJournal(base)
    def packet(rec):return reader.packet(rec.get('logical',rec.get('path')),rec['sha256'])
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    cp=reader.json(M/'S11c_d_remaining_case_frequency_material_row_pilot_checkpoint.json',CP_SHA);require(cp['status']=='ACCEPTED_BOUNDED_MATERIAL_1D_ROW_PILOT','clean accepted pilot')
    pc=meta(cp['artifacts']['manifest']);pa=pc['artifacts'];old=packet(cp['rowInput']);old_source=packet(old['sourceAction']);route=meta(pa['saved-source-and-rule-routes.json']);node_packet=packet(route['sourceNodes']['packet']);nodes=node_packet
    for key in route['sourceNodes']['keys']:nodes=nodes[key]
    fine_rule=packet(route['rules']['1024']);coarse_rule=packet(route['rules']['512']);momentum=fine_rule['nodes'];weights=fine_rule['weights'];require(same(momentum[::2],coarse_rule['nodes']),'actual saved rule nested nodes')
    # The pilot saved full fine arrays through an explicit coarse/new merge. Read
    # that actual completed merge; never pretend it was a standalone old call.
    merge=meta(pa['grids/1024/exact-coarse-reuse-completed.json']);require(merge['coefficient']==pa['grids/1024/full-coefficient.pickle'] and merge['fourier']==pa['grids/1024/full-fourier.pickle'],'actual saved complete fine-grid array owner')
    coefficient_cache=[{'key':coefficient_key(old,momentum),'array':packet(merge['coefficient']),'value':merge['coefficient'],'input':merge['input'],'receipt':pa['grids/1024/exact-coarse-reuse-completed.json'],'origin':'ACCEPTED_MERGED_GRID'}]
    fourier_cache=[{'key':fourier_key(old,momentum,nodes,old_source),'array':packet(merge['fourier']),'value':merge['fourier'],'input':merge['input'],'receipt':pa['grids/1024/exact-coarse-reuse-completed.json'],'origin':'ACCEPTED_MERGED_GRID'}]
    weighted_cache=[]
    for panels in (512,1024):
        folder='grids/'+str(panels);wi=packet(pa[folder+'/weighted-coefficient-input.pickle']);c=packet(wi['coefficient']);w=packet(wi['rule'])['weights'];receipt=meta(pa[folder+'/weighted-coefficient-completed.json'])
        require(receipt=={'input':pa[folder+'/weighted-coefficient-input.pickle'],'value':pa[folder+'/weighted-coefficient-value.pickle']},'actual saved weighted coefficient call')
        weighted_cache.append({'key':{'coefficient':c,'weights':w,'unit':old['row']['factors'][0]['unit'],'physical':physical(old)},'array':packet(receipt['value']),'value':receipt['value'],'input':receipt['input'],'receipt':pa[folder+'/weighted-coefficient-completed.json'],'origin':'ACCEPTED'})
    reader.retain(pilot.__file__,PILOT_SHA)
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        reader.retain(path);dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    env=dict(vars(pilot),NAME=NAME,__file__=str(Path(__file__).resolve()),SHARED_PACKET_CACHE={},frozen_sources=frozen_sources,batch_scope=scope)
    exec(compile(loader_adapter(Path(pilot.__file__).read_text()),'<whole saved-input loader only>','exec'),env)
    journal.json('batch-scope.json',{'rows':list(ROWS),'selectedCoarseComparisons':list(COMPARE),'finePanels':1024,'frequency':{'real':1.,'imag':-.01},'previousGridSeconds':[g['wallSeconds'] for g in cp['grids']],
        'stopping':'Only14missing material1D rows; cost gate3x maximum prior completed whole-row cost(initial10s allowance)+60s within900s. No automatic extension.','precision':{'relative':.01,'absoluteAmplitude':1e-4,'absoluteCurrent':1e-6},'nativeSourceProfileRuleEndCalls':0})
    journal.json('whole-loader-and-prior-result-joins.json',{'pilotHelper':reader.retain(pilot.__file__,PILOT_SHA),'pilotCheckpoint':reader.retain(M/'S11c_d_remaining_case_frequency_material_row_pilot_checkpoint.json',CP_SHA),'wholeOriginalLoaderPrefixReverseAST':True,'oldPilotNumericalTailExecuted':False,'acceptedFullGridMerge':pa['grids/1024/exact-coarse-reuse-completed.json']})
    counts={'newCoefficientBatches':0,'savedCoefficientGridUses':0,'newFourierBatches':0,'savedFourierGridUses':0,'newWeightedCoefficientCalls':0,'savedWeightedCoefficientUses':0,'newRowContractions':0,'newSourceProfileRuleEndCalls':0};results=[];costs=[];all_routes={}
    def cache_request(j,folder,key,entries,calculate,kind):
        matches=[v for v in entries if same(v['key'],key)];arg=j.write(folder+'/input.pickle',{'requested':key,'actualMatches':[{'input':v['input'],'value':v['value'],'receipt':v['receipt'],'origin':v['origin']} for v in matches]})
        if matches:
            owner=matches[0];value,vr=owner['array'],owner['value'];require(all(same(value,v['array']) for v in matches),'all exact full-input returns agree');disposition='REUSED_'+owner['origin'];counts['saved'+kind+'Uses']+=1
        else:
            value=calculate(arg);vr=j.write(folder+'/value.pickle',value);disposition='NEW';counts[{'CoefficientGrid':'newCoefficientBatches','FourierGrid':'newFourierBatches','WeightedCoefficient':'newWeightedCoefficientCalls'}[kind]]+=1
        name=folder+'/completed.json';j.json(name,{'input':arg,'value':vr,'disposition':disposition})
        if not matches:entries.append({'key':key,'array':value,'value':vr,'input':arg,'receipt':j.artifacts[name],'origin':'COMPLETED_NEW'})
        require(value.dtype==np.dtype(complex) and np.isfinite(value).all(),'full finite actual numerical cache return');return value,vr
    for ri in ROWS:
        remaining=900-(time.monotonic()-started);reserve=3*max([10.]+costs)+60;journal.json('cost-decisions/'+str(ri)+'.json',{'remainingSeconds':remaining,'requiredReserveSeconds':reserve,'previousRowSeconds':costs});require(remaining>reserve,'next material row fits measured whole-row cost')
        tick=time.monotonic();rb=base/('row-'+str(ri));rb.mkdir();env['ROW']=ri;s=env['load'](rb);rr,j=s['reader'],s['journal'];inp=s['packet'](s['full']);require(same(s['rules'][1024][0],fine_rule) and same(s['rules'][512][0],coarse_rule),'whole saved numerical rules identical')
        ckey=coefficient_key(inp,momentum)
        def new_coefficient(arg):
            with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):return np.broadcast_to(np.asarray(s['evaluate'](ckey['expression'],ckey['environment']),complex),(1025,129))
        c,cr=cache_request(j,'coefficient-grid',ckey,coefficient_cache,new_coefficient,'CoefficientGrid')
        fkey=fourier_key(inp,momentum,s['nodes'],s['action'])
        def new_fourier(arg):
            with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):phase=np.exp(-1j*momentum[:,None]*s['nodes'][None,:])
            pr=j.write('fourier-grid/phase-value.pickle',phase);ar=j.write('fourier-grid/contraction-input.pickle',{'phase':pr,'sourceAction':inp['sourceAction'],'requestedInput':arg,'operation':'phase @ weightedSourceAction'})
            return phase@s['action']
        sf,sfr=cache_request(j,'fourier-grid',fkey,fourier_cache,new_fourier,'FourierGrid')
        if 'fourier-grid/contraction-input.pickle' in j.artifacts:j.json('fourier-grid/contraction-completed.json',{'input':j.artifacts['fourier-grid/contraction-input.pickle'],'value':sfr})
        require(c.shape==sf.shape==(1025,129),'complete fine coefficient/Fourier grid');grids=[];arrays=[]
        for panels in ((512,1024) if ri in COMPARE else (1024,)):
            grid=coarse_rule if panels==512 else fine_rule;cv=c[::2] if panels==512 else c;fv=sf[::2] if panels==512 else sf;w=grid['weights'];folder='grids/'+str(panels)
            selected=j.write(folder+'/operand-routes.pickle',{'coefficient':cr,'sourceFourier':sfr,'indices':np.arange(0,1025,2) if panels==512 else np.arange(1025),'rule':s['rules'][panels][1],'ownRowInput':s['full']})
            wkey={'coefficient':cv,'weights':w,'unit':inp['row']['factors'][0]['unit'],'physical':physical(inp)}
            weighted,wr=cache_request(j,folder+'/weighted-coefficient',wkey,weighted_cache,lambda arg:cv*w[:,None],'WeightedCoefficient')
            ar=j.write(folder+'/row-input.pickle',{'weightedCoefficient':wr,'sourceFourier':sfr,'sourceIndices':np.arange(0,1025,2) if panels==512 else np.arange(1025),'operandRoutes':selected,'ownRowInput':s['full'],'operation':'weightedCoefficient.T @ sourceFourier','transposeConjugates':False})
            value=weighted.T@fv;vr=j.write(folder+'/row-value.pickle',value);j.json(folder+'/row-completed.json',{'input':ar,'value':vr});counts['newRowContractions']+=1;require(value.shape==(129,129) and np.isfinite(value).all(),'full finite material row')
            require(same(rr.packet(vr['path'],vr['sha256']),value),'exact saved full row readback');grids.append({'panels':panels,'input':ar,'value':vr,'maxNorm':float(np.max(abs(value)))});arrays.append(value)
        comparison=None
        if len(grids)==2:
            ci=j.write('comparison-input.pickle',{'coarse':grids[0]['value'],'fine':grids[1]['value'],'absoluteTarget':1e-4,'relativeTarget':.01});difference=arrays[1]-arrays[0];absolute=float(np.max(abs(difference)));norm=grids[1]['maxNorm'];target=max(1e-4,.01*norm)
            vv=j.write('comparison-value.pickle',{'difference':difference,'absoluteSpread':absolute,'fineMaxNorm':norm,'relativeSpread':absolute/norm if norm else None,'target':target,'withinTarget':absolute<=target});comparison={'input':ci,'value':vv,'absoluteSpread':absolute,'fineMaxNorm':norm,'relativeSpread':absolute/norm if norm else None,'target':target,'withinTarget':absolute<=target};j.json('comparison.json',comparison)
        cost=time.monotonic()-tick;costs.append(cost);result={'case':pilot.CASE,'rowIndex':ri,'sourceIndex':s['si'],'rowInput':s['full'],'rowReturn':grids[-1]['value'],'grids':grids,'comparison':comparison,'wallSeconds':cost,'savedReadback':True,'sourceAction':inp['sourceAction']}
        j.json('result.json',result);rr.postcheck();j.json('inputs.json',{'consumedRoutes':dict(rr.routes)});j.json('checks.json',{'result':result,'consumedRoutes':dict(rr.routes),'artifacts':dict(j.artifacts),'scope':scope(ri,s['settings'])['scope']});all_routes.update(rr.routes)
        for rec in j.artifacts.values():reader.retain(rec['path'],rec['sha256'])
        journal.json('row-results/'+str(ri)+'.json',{'result':result,'manifest':reader.retain(rb/'checks.json')});results.append(result)
    all_routes.update(reader.routes)
    for path,r in all_routes.items():reader.retain(path,r['sha256'])
    reader.postcheck();journal.json('inputs.json',{'consumedRoutes':dict(reader.routes),'acceptedPilot':reader.retain(M/'S11c_d_remaining_case_frequency_material_row_pilot_checkpoint.json',CP_SHA)})
    checks={'status':'COMPLETED_BOUNDED_MATERIAL_1D_ROWS','rows':results,'counts':counts,'consumedRoutes':dict(reader.routes),'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':'Fourteen missing material1D rows with two selected fixed-grid comparisons. Own source/unit/grade/context/Abel/settings/end provenance retained; no native source/profile/rule/end calls. No scattering uncertainty/pole/domain or whole-program completion claim.'}
    journal.json('checks.json',checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
