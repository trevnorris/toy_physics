#!/usr/bin/env python3
"""Only missing native rational end entries with complete operation receipts."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import resource
import signal
import textwrap
import time
from types import SimpleNamespace
import S11c_d_remaining_case_frequency_end_threshold_inputs as previous

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_threshold_inputs_checkpoint.json'
CP_SHA='557215f6c1cb16d75b89ecad85e00a2deea723f9ed33887c7a31ab8cfe5aeb59'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_rational_plan.md'
OWNER=previous.OWNER
same=previous.same
equal=previous.equal
source_body=previous.source_body


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_THRESHOLD_INPUTS','accepted threshold input checkpoint')
    for root,stage,sha in ((origin.parent,'frequency_end_threshold_inputs',cp['checksSha256']),(vr,'validate',cp['validation']['checksSha256'])):
        h.source.receipts.inspect_guard(root,stage)
        checks=origin/'checks.json' if stage!='validate' else vr/'checks.json'
        f.require(f.digest(checks)==sha and checks.read_bytes()==(root/(stage+'.stdout')).read_bytes(),'clean accepted checks/stdout')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedThresholdInputs':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    ref=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():ref(origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        p=origin/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal'] and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'original threshold input reference identity')
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        ref(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-threshold-input-checkpoint.json'),(origin/'checks.json','accepted-threshold-input-checks.json'),
                      (origin/'inputs.json','accepted-threshold-input-manifest.json'),(vr/'checks.json','accepted-threshold-input-validation.json')):ref(path,name,f.digest(path))
    for name,item in cp['validation']['artifacts'].items():ref(vr/name,'accepted-threshold-input-validation/'+name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'fresh rational constructor source')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return cp,manifest


def native_adapter():
    """Observe the entire unchanged native per-entry body and its exact calls."""
    original=ast.parse(textwrap.dedent(inspect.getsource(h.inputs.chart.end_tables))).body[0]
    outer=next(n for n in original.body if isinstance(n,ast.For))
    rows=next(n for n in outer.body if isinstance(n,ast.For))
    columns=next(n for n in rows.body if isinstance(n,ast.For))
    calls={'sp.cancel':'cancel','sp.fraction':'fraction','sp.Poly':'Poly','a.terms':'terms','b.terms':'terms'}
    sites=[]
    class Route(ast.NodeTransformer):
        def visit_Call(self,node):
            native=copy.deepcopy(node);name=ast.unparse(node.func)
            node=self.generic_visit(node)
            if name not in calls:return node
            f.require(not node.keywords,'literal native positional rational calls')
            site=len(sites);sites.append({'site':site,'call':ast.unparse(native),'nativeAST':ast.dump(native)})
            args=[ast.Constant(calls[name]),ast.Constant(site)]
            if name in ('a.terms','b.terms'):args.append(node.func.value)
            else:args.extend(node.args)
            return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='router',ctx=ast.Load()),attr='operation',ctx=ast.Load()),args=args,keywords=[]),node)
    routed=Route().visit(copy.deepcopy(columns))
    observed=[];assignments=[]
    for statement in routed.body:
        observed.append(statement)
        if isinstance(statement,ast.Assign):
            names=[n.id for target in statement.targets for n in ast.walk(target) if isinstance(n,ast.Name)]
            for name in names:
                assignments.append(name)
                observed.extend(ast.parse(f"router.observe('{name}', {name})").body)
        if isinstance(statement,ast.Expr) and isinstance(statement.value,ast.Call) and ast.unparse(statement.value.func)=='rows.append':
            observed.extend(ast.parse("router.observe('entry', rows[-1])").body)
    class Reverse(ast.NodeTransformer):
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=='router.observe':return None
            return self.generic_visit(node)
        def visit_Call(self,node):
            if ast.unparse(node.func)=='router.operation':
                return ast.parse(sites[node.args[1].value]['call'],mode='eval').body
            return self.generic_visit(node)
    whole=copy.deepcopy(original)
    wo=next(n for n in whole.body if isinstance(n,ast.For));wr=next(n for n in wo.body if isinstance(n,ast.For));wc=next(n for n in wr.body if isinstance(n,ast.For))
    wc.body=observed;restored=Reverse().visit(copy.deepcopy(whole))
    f.require(ast.dump(restored)==ast.dump(original),'entire native end_tables reversal: only exact call routing and persistence')
    fn=ast.parse('def native_entry(source,i,j,router):\n pass').body[0]
    fn.body=ast.parse("k,z=source['momentum'],source['radical']\nrows=[]").body+copy.deepcopy(observed)+ast.parse('return rows[0]').body
    module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]));env={'f':f}
    exec(compile(module,'<whole native rational entry with persistence>','exec'),env)
    return env['native_entry'],{'wholeNativeAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'entireOriginalBodyReverseAST':True,'nativeEntryStatements':len(columns.body),'sites':sites,
        'observedAssignments':assignments,'adapterSource':ast.unparse(module)}


def joined_sources(base,manifest):
    result={name:source_body(value) for name,value in {
        'endTables':h.inputs.chart.end_tables,'endSources':h.source.q.end_sources,
        'rationalDeterminant':engine.FullPencilModes.rational_determinant,
        'acceptedThresholdAdapter':previous.accepted_threshold.end_adapter,
        'acceptedRealAxisAnalysis':previous.accepted_threshold.real_axis_analysis,
        'coefficientAdapter':previous.native_adapter,'coefficientConstructor':previous.construct,
        'bindingSeedMatrix':previous.previous.Router.seed_matrix,
        'wholePair':h.end_native.Pair,'continuation':h.end_native.continue_pair}.items()}
    old=json.loads((base/'native-end-threshold-input-joins.json').read_text())
    for name,item in result.items():
        relative=str(Path(item['path']).relative_to(f.ROOT))
        f.require(manifest['sourceFiles'][relative]==item['sha256'],'accepted native/caller current and frozen source')
    for key in ('endTables','endSources','rationalDeterminant','acceptedThresholdAdapter','acceptedRealAxisAnalysis'):same(result[key],old[key])
    native,proof=native_adapter();result['entryAdapter']=proof
    for name,value in (('cancel',sp.cancel),('fraction',sp.fraction),('Poly',sp.Poly),('terms',sp.Poly.terms)):
        result['sympy_'+name]=source_body(value)
    result['sympyVersion']=sp.__version__
    f.save(base/'native-end-rational-joins.json',result)
    return native,result


def forbidden(*args,**kwargs):raise RuntimeError('only unfinished native rational entry calls are enabled')


def prohibit():
    # Keep precisely the original functions used by the unchanged entry body.
    operations={name:getattr(sp,name) for name in ('cancel','fraction','Poly')}
    previous.prohibit()
    for name,value in operations.items():setattr(sp,name,value)
    for name in ('load','construct','native_adapter','main'):setattr(previous,name,forbidden)
    previous.CoefficientRouter.convert=forbidden
    return operations


class OperationRouter:
    def __init__(self,base,own,functions):
        self.base,self.own,self.functions=base,own,functions
        self.calls=[];self.atlas=[];self.folder=None;self.key=None;self.index=None
        # Completed native Matrix.applyfunc(cancel) inputs/results are already
        # saved. Retain their actual raw-matrix positions and unit basis; no
        # old normalization or proof check is repeated to prepare this cache.
        packet=own['fullBindingInput'];raw=f.unpickle(base/'seed-matrix-join/raw-residual.pickle')
        receipt=json.loads((base/'seed-matrix-join/completed.json').read_text())
        seed_atlas=[];units=packet['units']['PENCIL_PLUS']
        for call in receipt['nativeCancelCalls']:
            folder=base/'seed-matrix-join'/str(call['index']);ip=folder/'input.pickle';vp=folder/'value.pickle'
            f.require(f.digest(ip)==call['inputSha256'] and f.digest(vp)==call['valueSha256'],'accepted actual scalar receipt identity')
            argument,value=f.unpickle(ip),f.unpickle(vp)
            indices=[(i,j) for i in range(5) for j in range(5) if equal(raw[i,j],argument)]
            f.require(bool(indices),'saved scalar call occurs in its actual raw matrix')
            for i,j in indices:
                key=self.arguments('cancel',(argument,),units[i][j])
                item={'input':key,'value':value,'owner':{'kind':'accepted-seed-scalar','call':call['index'],'index':(i,j),
                    'inputPath':str(ip),'valuePath':str(vp),'inputSha256':f.digest(ip),'valueSha256':f.digest(vp)},'receipt':call}
                seed_atlas.append(item);self.atlas.append(item)
        f.atomic_pickle(base/'completed-native-scalar-reuse-atlas.pickle',{'rawMatrixPath':str(base/'seed-matrix-join/raw-residual.pickle'),
            'rawMatrixSha256':f.digest(base/'seed-matrix-join/raw-residual.pickle'),'actualMatrixUnits':units,'calls':seed_atlas,
            'nativeProducerSource':source_body(previous.previous.Router.seed_matrix)})

    def arguments(self,name,args,unit):
        p=self.own['source']
        return {'function':name,'args':args,'kwargs':{},'entryUnitContext':unit,
            'coordinates':(p['momentum'],p['radical'],p['frequency']),
            'coordinateUnits':((-1,0,0),(0,-1,0),(0,-1,0))}

    def start(self,folder,key,index):self.folder,self.key,self.index=folder,key,index

    def observe(self,name,value):
        path=self.folder/'operands';path.mkdir(exist_ok=True)
        f.atomic_pickle(path/(name+'.pickle'),value)

    def operation(self,name,site,*args):
        index=len(self.calls);folder=self.base/'rational-operations'/str(index);folder.mkdir(parents=True)
        key=self.arguments(name,args,self.key['unit'])
        matches=[v for v in self.atlas if equal(v['input'],key)]
        request={'requested':key,'nativeSite':site,'entryIndex':self.index,'owner':OWNER,
            'fullEntryInput':self.key,'matches':matches,'nativeCaller':'S11c_d_frequency_chart.end_tables'}
        f.atomic_pickle(folder/'input.pickle',request)
        if matches:
            value=matches[0]['value'];owner=matches[0]['owner']
            for other in matches[1:]:same(value,other['value'])
        else:
            if name=='terms':value=args[0].terms()
            else:value=self.functions[name](*args)
            owner={'kind':'new-rational-operation','call':index,'inputPath':str(folder/'input.pickle'),'valuePath':str(folder/'value.pickle')}
        f.atomic_pickle(folder/'value.pickle',value)
        receipt={'index':index,'function':name,'nativeSite':site,'entryIndex':self.index,'reused':bool(matches),'owner':owner,
            'inputSha256':f.digest(folder/'input.pickle'),'valueSha256':f.digest(folder/'value.pickle')}
        f.save(folder/'completed.json',receipt);self.calls.append(receipt)
        self.atlas.append({'input':key,'value':value,'owner':owner,'receipt':receipt})
        return value


def construct(base,cp,native,functions):
    own=f.unpickle(base/'remaining-case-end-threshold-inputs.pickle')
    same(own['owner'],OWNER)
    source=own['source'];families=own['pendingEntryFamilies'];saved=f.unpickle(base/'saved-rational-entry-atlas.pickle')
    units=own['fullBindingInput']['units']['PENCIL_PLUS']
    router=OperationRouter(base,own,functions);entries=[];routes=[];completed={}
    f.require(len(families)==10,'accepted ten pending full native rational calls')
    for route in own['rationalEntryRoutes']:
        i,j=route['index'];index=(i,j);pairpath=base/'rational-entry-inputs'/f'{i}-{j}.pickle'
        f.require(f.digest(pairpath)==route['inputSha256'],'accepted own entry input hash')
        pair=f.unpickle(pairpath);key=pair['requested']
        same(pair['ownSource'],source);same(pair['owner'],OWNER);same(pair['ownIndex'],index)
        same(key,{'original':source['livePencil'][i,j],'momentum':source['momentum'],'radical':source['radical'],'frequency':source['frequency'],'unit':units[i][j]})
        folder=base/'rational-entries'/f'{i}-{j}';folder.mkdir(parents=True)
        accepted=pair['acceptedMatches'];start=len(router.calls)
        f.atomic_pickle(folder/'requested-input.pickle',{'pair':pair,'path':str(pairpath),'sha256':f.digest(pairpath),'acceptedRoute':route,'wholeOwnInputReference':{'path':str(base/'remaining-case-end-threshold-inputs.pickle'),'sha256':f.digest(base/'remaining-case-end-threshold-inputs.pickle')}})
        if accepted:
            first=accepted[0];same(first['input'],key)
            f.require(any(equal(first,item) for item in saved),'exact accepted native input/full return owner')
            # The accepted computational owner retains its original row/column.
            # Only the new consumer's row/column metadata is placed in the table.
            original=first['value'];entry=dict(original,row=i,column=j)
            f.atomic_pickle(folder/'accepted-return-pair.pickle',{'requested':key,'savedFullReturn':original,'savedOwner':first['owner'],
                'savedPath':first['path'],'savedSha256':first['sha256'],'ownIndex':index,'ownReturn':entry})
            for other in accepted[1:]:
                same(other['input'],key)
                for field in ('original','numeratorTerms','denominatorTerms','reconstructionResidual'):same(original[field],other['value'][field])
            receipt={'index':index,'mode':'accepted-full-entry','owner':first['owner'],'originalPath':first['path']}
        else:
            family=route['pendingFamily'];same(families[family]['input'],key)
            if family in completed:
                first=completed[family];entry=dict(first['entry'],row=i,column=j)
                f.atomic_pickle(folder/'completed-return-pair.pickle',{'requested':key,'saved':first,'ownReturn':entry})
                receipt={'index':index,'mode':'completed-new-full-entry','owner':first['index'],'family':family}
            else:
                same(families[family]['firstIndex'],index);router.start(folder,key,index)
                entry=native(source,i,j,router)
                completed[family]={'entry':entry,'index':index,'path':str(folder/'entry.pickle')}
                receipt={'index':index,'mode':'new-native-entry','owner':index,'family':family}
        f.atomic_pickle(folder/'entry.pickle',entry)
        receipt.update({'requestedInputSha256':f.digest(folder/'requested-input.pickle'),'valueSha256':f.digest(folder/'entry.pickle'),
            'operationIndices':list(range(start,len(router.calls)))})
        f.save(folder/'completed.json',receipt);entries.append(entry);routes.append(receipt)
        f.require(entry['reconstructionResidual']==0,'saved native rational reconstruction residual')
    f.require(len(completed)==len(families) and len(entries)==25,'complete pending family and full own table coverage')
    result={'source':source,'entries':entries,'scope':'Actual full end rational coefficients for later complete invariant-pair continuation; no end mode has been continued here.'}
    target=base/'new-rational-end';target.mkdir();f.atomic_pickle(target/'right-rational-end.pickle',result)
    aliases={}
    for label,case in own['aliases'].items():
        aliases[label]={}
        for end,prior in case.items():
            shared=tuple(prior['owner'])==OWNER
            item={'address':(label,end),'owner':prior['owner'],'fullInputRoute':prior,
                'tablePath':str(target/'right-rational-end.pickle') if shared else prior['tablePath'],
                'mode':'shared-new-rational-table' if shared else 'accepted-baseline-rational-table',
                'endDomainOrNumericalReuseAccepted':False}
            path=base/'rational-cases'/label/end.lower();path.mkdir(parents=True);f.save(path/'route.json',item);aliases[label][end]=item
        f.save(base/'rational-cases'/label/'case-summary.json',aliases[label])
    # Use an actual newly constructed entry for full-input/value routing controls.
    index=families[0]['firstIndex'];key=families[0]['input'];entry=completed[0]['entry']
    wrong_input=dict(key,original=key['original']+1)
    wrong_unit=dict(key,unit=(key['unit'][0]+1,*key['unit'][1:]))
    wrong_value=dict(entry,reconstructionResidual=entry['reconstructionResidual']+1)
    changed_terms=list(entry['numeratorTerms']);powers,coefficient=changed_terms[0];changed_terms[0]=(powers,coefficient+1)
    wrong_terms=dict(entry,numeratorTerms=changed_terms)
    mutations={'actualInput':key,'coefficientInput':wrong_input,'unitInput':wrong_unit,'actualOwner':OWNER,'wrongOwner':(OWNER[0],'LEFT'),
        'actualValue':entry,'residualValue':wrong_value,'coefficientValue':wrong_terms,'ownIndex':index}
    f.atomic_pickle(base/'rational-mutation-operands.pickle',mutations)
    controls={'inputCoefficient':not equal(key,wrong_input),'physicalUnit':not equal(key,wrong_unit),'owner':OWNER!=mutations['wrongOwner'],
        'residual':not equal(entry,wrong_value),'returnedCoefficient':not equal(entry,wrong_terms)}
    f.save(base/'rational-mutation-controls.json',controls);f.require(all(controls.values()),'actual rational route/value mutations')
    counts={name:{'new':sum(c['function']==name and not c['reused'] for c in router.calls),
                  'reused':sum(c['function']==name and c['reused'] for c in router.calls)} for name in ('cancel','fraction','Poly','terms')}
    f.save(base/'rational-operation-inventory.json',{'calls':router.calls,'counts':counts,'entries':routes})
    f.atomic_pickle(base/'remaining-case-end-rational.pickle',{'owner':OWNER,'sourceInput':own,'table':result,'entries':routes,
        'operations':router.calls,'aliases':aliases,'thresholdAndEndDomainPending':True})
    return {'physicalEnds':8,'sharedNewUses':2,'baselineAliases':6,'newEndOwners':1,'cases':aliases,
        'rationalEntries':25,'newEntryCalls':len(completed),'acceptedEntryUses':sum(r['mode']=='accepted-full-entry' for r in routes),
        'operationCounts':counts,'actualMutationControls':len(controls),'newDeterminantThresholdBindingDerivativeModeCurrentCalls':0,
        'endDomainContinuationOrNumericalReuseAccepted':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.previous.previous.previous.protect(base)
    cp,manifest=load(base);native,joins=joined_sources(base,manifest);functions=prohibit();result=construct(base,cp,native,functions)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen source posthash')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'reference postidentity')
    for key,item in joins.items():
        if key.startswith('sympy_'):f.require(f.digest(Path(item['path']))==item['sha256'],'native SymPy operation source posthash')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_RATIONAL_ENTRIES','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
