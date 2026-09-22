#!/usr/bin/env python3
"""Missing native threshold coefficient images and saved rational-call routing."""
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
import S11c_d_remaining_case_frequency_end_bindings as previous
import S11c_d_frequency_source_finish as accepted_threshold

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_bindings_checkpoint.json'
CP_SHA='c556a8ebb82b839558d36946b0317fb671cb2929ef1e7d7af294cfbd700249e6'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_threshold_inputs_plan.md'
OWNER=previous.OWNER
same=previous.same
equal=previous.previous.previous.modes_reader.same
source_body=previous.previous.source_body


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_BINDINGS','accepted binding checkpoint')
    for root,stage,sha in ((origin.parent,'frequency_end_bindings',cp['checksSha256']),(vr,'validate',cp['validation']['checksSha256'])):
        h.source.receipts.inspect_guard(root,stage)
        checks=origin/'checks.json' if stage!='validate' else vr/'checks.json'
        f.require(f.digest(checks)==sha and checks.read_bytes()==(root/(stage+'.stdout')).read_bytes(),'clean accepted checks/stdout')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedBindings':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    ref=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():ref(origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        p=origin/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal'] and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'original binding reference identity')
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        ref(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-binding-checkpoint.json'),(origin/'checks.json','accepted-binding-checks.json'),
                      (origin/'inputs.json','accepted-binding-manifest.json'),(vr/'checks.json','accepted-binding-validation.json')):ref(path,name,f.digest(path))
    for name,item in cp['validation']['artifacts'].items():ref(vr/name,'accepted-binding-validation/'+name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'fresh threshold input source')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return cp,manifest


def native_adapter():
    """Route exactly the two coefficient-image calls in the whole native body."""
    original=ast.parse(textwrap.dedent(inspect.getsource(h.source.q.end_sources))).body[0]
    loop=next(n for n in original.body if isinstance(n,ast.For))
    start=next(i for i,n in enumerate(loop.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='zero' for t in n.targets))
    operations=copy.deepcopy(loop.body[start:start+2])
    routes={
        'live.subs(modes.k, 0).xreplace({frequency: coordinate, modes.q: radical_coordinate})':"router.convert('coefficientMatrix', live, modes.k, frequency, modes.q, coordinate, radical_coordinate)",
        'wave.subs(modes.k, 0).xreplace({frequency: coordinate, modes.q: radical_coordinate})':"router.convert('wave', wave, modes.k, frequency, modes.q, coordinate, radical_coordinate)"}
    counts={s:0 for s in routes}
    class Forward(ast.NodeTransformer):
        def visit(self,node):
            if isinstance(node,ast.expr) and ast.unparse(node) in routes:
                key=ast.unparse(node);counts[key]+=1;return ast.copy_location(ast.parse(routes[key],mode='eval').body,node)
            return super().visit(node)
    adapted=[Forward().visit(n) for n in copy.deepcopy(operations)]
    f.require(set(counts.values())=={1},'two exact native coefficient-image calls')
    reverse={v:k for k,v in routes.items()}
    class Reverse(ast.NodeTransformer):
        def visit(self,node):
            if isinstance(node,ast.expr) and ast.unparse(node) in reverse:return ast.copy_location(ast.parse(reverse[ast.unparse(node)],mode='eval').body,node)
            return super().visit(node)
    whole=copy.deepcopy(original);wl=next(n for n in whole.body if isinstance(n,ast.For))
    wl.body[start:start+2]=[Reverse().visit(n) for n in copy.deepcopy(adapted)]
    f.require(ast.dump(whole)==ast.dump(original),'entire end_sources reverse AST with unchanged completed prefix and pending threshold tail')
    fn=ast.parse('def coefficient_images(live,wave,modes,frequency,coordinate,radical_coordinate,router):\n pass').body[0]
    fn.body=adapted+ast.parse('return zero,curve').body
    module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]));env={}
    exec(compile(module,'<native coefficient-frame input operations>','exec'),env)
    return env['coefficient_images'],{'wholeNativeAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'routes':routes,'counts':counts,'completedPrefixStatements':start,'coefficientStatements':len(operations),
        'pendingNativeStatements':len(loop.body)-start-2,'entireOriginalBodyReverseAST':True,'adapterSource':ast.unparse(module)}


def joined_sources(base,manifest):
    native={'endSources':h.source.q.end_sources,'endTables':h.inputs.chart.end_tables,
        'rationalDeterminant':engine.FullPencilModes.rational_determinant,
        'acceptedThresholdAdapter':accepted_threshold.end_adapter,'acceptedRealAxisAnalysis':accepted_threshold.real_axis_analysis,
        'bindingAdapter':previous.native_adapter,'bindingConstructor':previous.construct,
        'wholePair':h.end_native.Pair,'continuation':h.end_native.continue_pair,'seeds':h.end_native.seeds,'maps':h.end_native.maps}
    result={name:source_body(value) for name,value in native.items()}
    for item in result.values():
        name=str(Path(item['path']).relative_to(f.ROOT));f.require(manifest['sourceFiles'][name]==item['sha256'],'whole native source already accepted/frozen')
    old=json.loads((base/'native-end-binding-joins.json').read_text())
    for name in ('endSources','endTables'):same(result[name],old[name])
    # The accepted producer's full adapter replaces the old real-only gate by
    # its exact complex polynomial / real-common-factor analysis. Never run
    # the obsolete gate or repeat any accepted analysis during this stage.
    adapter_tree=ast.parse(result['acceptedThresholdAdapter']['source'])
    f.require(any(isinstance(n,ast.Constant) and isinstance(n.value,str) and 'analysis=real_axis_analysis(poly,coordinate)' in n.value for n in ast.walk(adapter_tree)),'actual accepted complex-threshold route')
    frequency_cp=json.loads((base/'end-origins/frequency-checkpoint.json').read_text())
    for key in ('acceptedThresholdAdapter','acceptedRealAxisAnalysis'):
        path=Path(result[key]['path']);name=str(path.relative_to(f.ROOT))
        f.require(frequency_cp['sourceFiles'][name]==result[key]['sha256'],'original accepted frequency producer helper')
    fn,proof=native_adapter();result['coefficientAdapter']=proof
    result['scope']='Only missing coefficient images and saved threshold/table call routing; no elimination, polynomial analysis, table construction or end continuation.'
    f.save(base/'native-end-threshold-input-joins.json',result)
    return fn,result


def forbidden(*args,**kwargs):raise RuntimeError('threshold input stage cannot replay accepted science or construct thresholds/tables')


def prohibit():
    previous.prohibit()
    for module,names in ((previous,('load','construct','native_adapter','main')),
                         (accepted_threshold,('real_axis_analysis','end_adapter','reuse','main'))):
        for name in names:setattr(module,name,forbidden)
    for name in ('cancel','expand','simplify','factor','factor_list','resultant','solve','nroots','diff','gcd','gcdex','div','lambdify'):setattr(sp,name,forbidden)
    engine.FullPencilModes.rational_determinant=staticmethod(forbidden)


class CoefficientRouter:
    def __init__(self,base,own,baselines,units):
        self.base,self.own,self.units=base,own,units;self.atlas=[];self.receipts=[]
        for label,item in baselines.items():
            p,t=item['pencil'],item['threshold']
            for name,key,unit in (('coefficientMatrix','livePencil',item['matrixUnits']),('wave','wave',(0,-2,0))):
                call=self.arguments(p[key],p['momentum'],p['frequency'],p['radical'],t['frequencyCoordinate'],t['radicalCoordinate'],unit)
                self.atlas.append((call,t[name],{'end':label,'field':name,'path':item['thresholdPath'],'sourcePath':item['pencilPath']}))
        f.atomic_pickle(base/'saved-coefficient-call-atlas.pickle',self.atlas)

    @staticmethod
    def arguments(expression,momentum,frequency,radical,coordinate,radical_coordinate,unit):
        return {'expression':expression,'momentum':momentum,'normalMomentumCoefficient':0,'frequency':frequency,
            'radical':radical,'frequencyCoordinate':coordinate,'radicalCoordinate':radical_coordinate,'inputUnit':unit,
            'frequencyUnit':(0,-1,0),'radicalUnit':(0,-1,0),'momentumUnit':(-1,0,0),'coefficientCoordinateUnits':((0,0,0),(0,0,0))}

    def convert(self,name,expression,momentum,frequency,radical,coordinate,radical_coordinate):
        unit=self.units if name=='coefficientMatrix' else (0,-2,0)
        requested=self.arguments(expression,momentum,frequency,radical,coordinate,radical_coordinate,unit)
        folder=self.base/'coefficient-operations'/name;folder.mkdir(parents=True)
        matches=[(key,value,owner) for key,value,owner in self.atlas if equal(key,requested)]
        f.atomic_pickle(folder/'input.pickle',{'requested':requested,'fullOwnInput':self.own,'candidates':matches})
        if matches:
            value=matches[0][1];owner=matches[0][2]
            for key,other,_ in matches:same(key,requested);same(value,other)
        else:
            # These are the literal two native operations, each persisted
            # immediately. No native binding or prior coordinate map repeats.
            substituted=expression.subs(momentum,0)
            f.atomic_pickle(folder/'normal-zero-value.pickle',substituted)
            value=substituted.xreplace({frequency:coordinate,radical:radical_coordinate})
            owner={'end':OWNER,'field':name,'path':str(folder/'value.pickle')}
        f.atomic_pickle(folder/'value.pickle',value)
        receipt={'name':name,'reused':bool(matches),'owner':owner,'inputSha256':f.digest(folder/'input.pickle'),
            'valueSha256':f.digest(folder/'value.pickle'),'nativeOperations':[] if matches else ['subs(normalMomentum,0)','xreplace(actual frequency/radical coefficient coordinates)']}
        f.save(folder/'completed.json',receipt);self.receipts.append(receipt)
        self.atlas.append((requested,value,owner))
        f.require(not value.free_symbols-{coordinate,radical_coordinate},'actual coefficient-frame complete variable binding')
        return value


def entry_input(source,row,column,units):
    return {'original':source['livePencil'][row,column],'momentum':source['momentum'],'radical':source['radical'],
        'frequency':source['frequency'],'unit':units[row][column]}


def construct(base,cp,native):
    own=f.unpickle(base/'remaining-case-end-bindings.pickle');packet=own['pencil'];context=own['coordinateInput']['sourceInput']['context']
    same(own['owner'],OWNER);same(packet,f.unpickle(base/'new-end/right-frequency-pencil.pickle'))
    units=previous.matrix_units(context['frame'],(0,0,0));same(units,own['units']['PENCIL_PLUS'])
    baselines={}
    for label in ('REFERENCE','LEFT','RIGHT'):
        pp=base/'end-accepted/frequency'/(label.lower()+'-frequency-pencil.pickle');tp=base/'end-accepted/frequency'/(label.lower()+'-threshold-candidates.pickle')
        p,t=f.unpickle(pp),f.unpickle(tp);cm=t['coordinateMap']
        u=previous.matrix_units({'fieldUnits':cm['fieldReferenceUnits'],'equationUnits':cm['equationReferenceUnits']},(0,0,0))
        baselines[label]={'pencil':p,'threshold':t,'matrixUnits':u,'pencilPath':str(pp),'thresholdPath':str(tp),
            'pencilSha256':f.digest(pp),'thresholdSha256':f.digest(tp)}
    f.atomic_pickle(base/'saved-native-threshold-input-values.pickle',{'own':own,'baselines':baselines})
    for label,item in baselines.items():
        t=item['threshold'];f.require('realAxisAnalysis' in t,'actual accepted complex-polynomial threshold packet')
        f.require(t['frequencyCoordinate'].name=='s11cdFrequencyCoefficientCoordinate' and t['frequencyCoordinate'].is_real is True and t['radicalCoordinate'].name=='s11cdRadicalCoefficientCoordinate' and t['radicalCoordinate'].is_complex is True,'exact original coordinate declarations')
        f.require(t['coordinateMap']['normalMomentumCoefficient']==0,'literal native zero normal momentum')
    t=baselines['RIGHT']['threshold'];coordinate,radical_coordinate=t['frequencyCoordinate'],t['radicalCoordinate']
    router=CoefficientRouter(base,own,baselines,units)
    zero,curve=native(packet['livePencil'],packet['wave'],SimpleNamespace(k=packet['momentum'],q=packet['radical']),packet['frequency'],coordinate,radical_coordinate,router)
    coordinate_map={'frequencyUnit':(0,-1,0),'radicalUnit':(0,-1,0),'normalMomentumUnit':(-1,0,0),'normalMomentumCoefficient':0,
        'unitFrame':context['seedInput']['frame'],'fieldReferenceUnits':context['frame']['fieldUnits'],'equationReferenceUnits':context['frame']['equationUnits'],
        'matrixConvention':t['coordinateMap']['matrixConvention']}
    requested={'coefficientMatrix':zero,'inputUnits':units,'frequencyCoordinate':coordinate,'radicalCoordinate':radical_coordinate}
    det_candidates=[];branch_candidates=[]
    for label,item in baselines.items():
        old=item['threshold'];det_input={'coefficientMatrix':old['coefficientMatrix'],'inputUnits':item['matrixUnits'],
            'frequencyCoordinate':old['frequencyCoordinate'],'radicalCoordinate':old['radicalCoordinate']}
        det_candidates.append({'owner':label,'input':det_input,'value':((old['numerator'],old['denominator']),old['clearedMatrix'],old['rowDenominators']),
            'threshold':old,'path':item['thresholdPath'],'sha256':item['thresholdSha256'],'exactInput':equal(requested,det_input)})
        branch_input={'wave':old['wave'],'radicalCoordinate':old['radicalCoordinate'],'frequencyCoordinate':old['frequencyCoordinate'],'unit':(0,0,0)}
        branch_candidates.append({'owner':label,'input':branch_input,'value':old['bulkBranchFrequencies'],'path':item['thresholdPath']})
    branch_request={'wave':curve,'radicalCoordinate':radical_coordinate,'frequencyCoordinate':coordinate,'unit':(0,0,0)}
    pairs={'owner':OWNER,'requestedDeterminant':requested,'determinantCandidates':det_candidates,'requestedBranches':branch_request,
        'branchCandidates':branch_candidates,'wave':curve,'coordinateMap':coordinate_map,'actualOwnInput':own}
    f.atomic_pickle(base/'threshold-native-call-input-pairs.pickle',pairs)
    det_matches=[d for d in det_candidates if d['exactInput']]
    for match in det_matches[1:]:same(match['value'],det_matches[0]['value'])
    branch_matches=[d for d in branch_candidates if equal(d['input'],branch_request)]
    for match in branch_matches[1:]:same(match['value'],branch_matches[0]['value'])
    elimination_matches=[];analysis_matches=[]
    pending={'determinant':[],'elimination':[],'analysis':[],'branches':[],'rationalEntries':[]}
    if not det_matches:pending['determinant'].append(requested)
    if det_matches:
        numerator=det_matches[0]['value'][0][0]
        elimination_input={'numerator':numerator,'wave':curve,'radicalCoordinate':radical_coordinate,'coefficientUnit':(0,0,0)}
        for label,item in baselines.items():
            old=item['threshold'];key={'numerator':old['numerator'],'wave':old['wave'],'radicalCoordinate':old['radicalCoordinate'],'coefficientUnit':(0,0,0)}
            if equal(elimination_input,key):elimination_matches.append({'owner':label,'input':key,'value':old['elimination'],'path':item['thresholdPath']})
        f.atomic_pickle(base/'elimination-saved-input-pairs.pickle',{'requested':elimination_input,'matches':elimination_matches})
        for item in elimination_matches[1:]:same(item['value'],elimination_matches[0]['value'])
        if not elimination_matches:pending['elimination'].append(elimination_input)
        else:
            analysis_input={'elimination':elimination_matches[0]['value'],'coordinate':coordinate,'analysisSource':source_body_from_join(base,'acceptedRealAxisAnalysis'),'callerSource':source_body_from_join(base,'acceptedThresholdAdapter'),'polynomialConstructor':'sp.Poly(elimination, coordinate)'}
            for label,item in baselines.items():
                old=item['threshold'];key={'elimination':old['elimination'],'coordinate':old['frequencyCoordinate'],'analysisSource':analysis_input['analysisSource'],'callerSource':analysis_input['callerSource'],'polynomialConstructor':analysis_input['polynomialConstructor']}
                if equal(key,analysis_input):analysis_matches.append({'owner':label,'input':key,'value':old['realAxisAnalysis'],'path':item['thresholdPath']})
            f.atomic_pickle(base/'analysis-saved-input-pairs.pickle',{'requested':analysis_input,'matches':analysis_matches})
            for item in analysis_matches[1:]:same(item['value'],analysis_matches[0]['value'])
            if not analysis_matches:pending['analysis'].append(analysis_input)
    if not branch_matches:pending['branches'].append(branch_request)
    threshold_summary={'determinantCandidateOwners':[v['owner'] for v in det_matches],'branchCandidateOwners':[v['owner'] for v in branch_matches],
        'eliminationCandidateOwners':[v['owner'] for v in elimination_matches],'analysisCandidateOwners':[v['owner'] for v in analysis_matches],
        'eliminationDeferredUntilDeterminant':not det_matches,'analysisDeferredUntilElimination':not elimination_matches,
        'newThresholdOrRootOperations':0,'completeEndDomainAccepted':False}
    f.save(base/'threshold-call-routing-summary.json',threshold_summary)
    # Every actual native rational-entry argument and accepted full return is
    # retained. Matching a scalar call is not numerical matrix/row reuse.
    table_atlas=[]
    for label in ('LEFT','RIGHT'):
        path=base/'end-accepted/chart'/(label.lower()+'-rational-end.pickle');table=f.unpickle(path)
        same(table['source'],baselines[label]['pencil']);f.require(len(table['entries'])==25,'accepted complete 5x5 rational table')
        for entry in table['entries']:
            i,j=entry['row'],entry['column'];same(entry['original'],table['source']['livePencil'][i,j])
            table_atlas.append({'input':entry_input(table['source'],i,j,baselines[label]['matrixUnits']),
                'value':entry,'owner':(label,i,j),'path':str(path),'sha256':f.digest(path)})
    f.atomic_pickle(base/'saved-rational-entry-atlas.pickle',table_atlas)
    entry_routes=[];pending_families=[]
    (base/'rational-entry-inputs').mkdir()
    for i in range(5):
        for j in range(5):
            key=entry_input(packet,i,j,units);matches=[v for v in table_atlas if equal(v['input'],key)]
            index=(i,j);path=base/'rational-entry-inputs'/f'{i}-{j}.pickle'
            packet_pair={'owner':OWNER,'ownIndex':index,'requested':key,'acceptedMatches':matches,'ownSource':packet,
                'ownUnitFrame':coordinate_map,'fullOwnInputPath':str(base/'remaining-case-end-bindings.pickle')}
            f.atomic_pickle(path,packet_pair)
            for candidate in matches[1:]:
                for name in ('original','numeratorTerms','denominatorTerms','reconstructionResidual'):same(candidate['value'][name],matches[0]['value'][name])
            route={'index':index,'inputPath':str(path),'inputSha256':f.digest(path),'acceptedOwners':[v['owner'] for v in matches]}
            if not matches:
                family=next((n for n,v in enumerate(pending_families) if equal(v['input'],key)),None)
                if family is None:family=len(pending_families);pending_families.append({'input':key,'firstIndex':index,'uses':[]})
                pending_families[family]['uses'].append(index);route['pendingFamily']=family;pending['rationalEntries'].append({'input':key,'index':index,'family':family})
            entry_routes.append(route);f.save(base/'rational-entry-inputs'/f'{i}-{j}.json',route)
    f.atomic_pickle(base/'pending-end-threshold-table-inputs.pickle',pending)
    f.atomic_pickle(base/'pending-rational-entry-families.pickle',pending_families)
    f.save(base/'rational-entry-routing-summary.json',{'entries':len(entry_routes),'acceptedUses':sum(bool(v['acceptedOwners']) for v in entry_routes),
        'pendingUses':len(pending['rationalEntries']),'pendingFamilies':len(pending_families),'routes':entry_routes,'newTableConstructionCalls':0})
    aliases={}
    for label,case in own['aliases'].items():
        aliases[label]={}
        (base/'threshold-input-cases'/label).mkdir(parents=True)
        for end,route in case.items():
            item={'address':(label,end),'bindingRoute':route,'completeOwnInputPath':str(base/'coordinate-cases'/label/end.lower()/'full-input-route.pickle'),
                'mode':'new-coefficient-input-owner' if tuple(route['owner'])==OWNER else 'accepted-baseline-threshold-table-input-route',
                'owner':route['owner'],'fullEndDomainOrNumericalReuseAccepted':False}
            if tuple(route['owner'])==OWNER:item['coefficientInputPath']=str(base/'threshold-native-call-input-pairs.pickle')
            else:
                item['thresholdPath']=baselines[end]['thresholdPath'];item['tablePath']=str(base/'end-accepted/chart'/(end.lower()+'-rational-end.pickle'))
            (base/'threshold-input-cases'/label/end.lower()).mkdir()
            f.save(base/'threshold-input-cases'/label/end.lower()/'route.json',item);aliases[label][end]=item
        f.save(base/'threshold-input-cases'/label/'case-summary.json',aliases[label])
    # Actual changed inputs, persisted before the exact comparison guards.
    changed=sp.MutableDenseMatrix(zero);changed[0,0]+=1;changed=sp.ImmutableMatrix(changed)
    changed_units=tuple(tuple((v[0]+1,v[1],v[2]) if (i,j)==(0,0) else v for j,v in enumerate(row)) for i,row in enumerate(units))
    original_entry=entry_input(packet,0,0,units);wrong_entry=dict(original_entry,original=original_entry['original']+1)
    mutations={'actualMatrix':zero,'changedMatrix':changed,'actualUnits':units,'changedUnits':changed_units,
        'actualOwner':OWNER,'changedOwner':(OWNER[0],'LEFT'),'actualCoordinate':coordinate,'changedCoordinate':2*coordinate,
        'actualEntryInput':original_entry,'changedEntryInput':wrong_entry}
    f.atomic_pickle(base/'threshold-input-mutation-operands.pickle',mutations)
    controls={'coefficientMatrix':not equal(zero,changed),'physicalEntryUnit':not equal(units,changed_units),'owner':OWNER!=mutations['changedOwner'],
        'coordinateCoefficient':not equal(coordinate,mutations['changedCoordinate']),'rationalEntryCoefficient':not equal(original_entry,wrong_entry)}
    f.save(base/'threshold-input-mutation-controls.json',controls);f.require(all(controls.values()),'actual threshold input routing controls')
    output={'owner':OWNER,'source':packet,'fullBindingInput':own,'coefficientMatrix':zero,'wave':curve,'coordinateMap':coordinate_map,
        'frequencyCoordinate':coordinate,'radicalCoordinate':radical_coordinate,'coefficientReceipts':router.receipts,'thresholdRoutes':threshold_summary,
        'rationalEntryRoutes':entry_routes,'pending':pending,'pendingEntryFamilies':pending_families,'aliases':aliases,
        'thresholdsTablesDomainsAndContinuationUnaccepted':True}
    f.atomic_pickle(base/'remaining-case-end-threshold-inputs.pickle',output)
    return {'cases':aliases,'physicalEnds':8,'newOwnerUses':2,'baselineAliases':6,'coefficientCalls':len(router.receipts),
        'newCoefficientImages':sum(not v['reused'] for v in router.receipts),'reusedCoefficientImages':sum(v['reused'] for v in router.receipts),
        'rationalEntryUses':len(entry_routes),'acceptedRationalEntryUses':sum(bool(v['acceptedOwners']) for v in entry_routes),
        'pendingRationalEntryUses':len(pending['rationalEntries']),'pendingRationalEntryFamilies':len(pending_families),
        'thresholdRoutes':threshold_summary,'actualMutationControls':len(controls),'newBindingDerivativeDeterminantEliminationTableRootCalls':0,
        'thresholdsTablesDomainsAndContinuationUnaccepted':True}


def source_body_from_join(base,name):
    return json.loads((base/'native-end-threshold-input-joins.json').read_text())[name]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.previous.previous.protect(base);cp,manifest=load(base);native,joins=joined_sources(base,manifest);prohibit();result=construct(base,cp,native)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen source posthash')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'reference postidentity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_THRESHOLD_INPUTS','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
