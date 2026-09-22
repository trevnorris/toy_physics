#!/usr/bin/env python3
"""Missing native live-end bindings, using accepted coordinate/tangent values."""
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
import S11c_d_remaining_case_frequency_end_coordinates as previous

h,f,sp,engine=previous.h,previous.f,previous.sp,previous.engine
CP=f.M/'S11c_d_remaining_case_frequency_end_coordinates_checkpoint.json'
CP_SHA='3a7811b4bf3584ea071f338aca105645d6e36786b2cfe0642bec41a10f743141'
PLAN=f.M/'S11c_d_remaining_case_frequency_end_bindings_plan.md'
OWNER=previous.OWNER
same=previous.same


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(f.digest(CP)==CP_SHA and cp['status']=='ACCEPTED_CASE_FREQUENCY_END_COORDINATE_IMAGES','accepted coordinate checkpoint')
    for root,stage,sha in ((origin.parent,'frequency_end_coordinates',cp['checksSha256']),(vr,'validate',cp['validation']['checksSha256'])):
        h.source.receipts.inspect_guard(root,stage)
        checks=origin/'checks.json' if stage!='validate' else vr/'checks.json'
        f.require(f.digest(checks)==sha and checks.read_bytes()==(root/(stage+'.stdout')).read_bytes(),'clean accepted checks/stdout')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],
              'acceptedCoordinates':{'checkpoint':str(CP),'checkpointSha256':CP_SHA,'runDirectory':str(origin),
                  'checksSha256':cp['checksSha256'],'validationDirectory':str(vr),'validatorChecksSha256':cp['validation']['checksSha256']}}
    ref=lambda path,name,sha:h.source.reference(base,manifest,path,name,sha)
    for name,item in cp['artifacts'].items():ref(origin/name,name,item['sha256'])
    for name,item in cp['referencedInputs'].items():
        p=origin/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal'] and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'original coordinate reference identity')
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==value,'accepted current/frozen source')
        ref(origin/'source'/name,'source/'+name,value)
    for path,name in ((CP,'accepted-coordinate-checkpoint.json'),(origin/'checks.json','accepted-coordinate-checks.json'),
                      (origin/'inputs.json','accepted-coordinate-manifest.json'),(vr/'checks.json','accepted-coordinate-validation.json')):ref(path,name,f.digest(path))
    for name,item in cp['validation']['artifacts'].items():ref(vr/name,'accepted-coordinate-validation/'+name,item['sha256'])
    # The original frequency producer consumed these actual saved seed pencils.
    # Read its input inventory; do not derive or re-evaluate a seed curve.
    original=json.loads((base/'end-origins/frequency-checkpoint.json').read_text())
    for end in ('left','right'):
        paths=[Path(p) for p in original['inputPackets'] if p.endswith('/'+end+'-symbolic.pickle')]
        f.require(len(paths)==1,'unique original native frequency seed source')
        path=paths[0];ref(path,'native-frequency-seeds/'+end+'-symbolic.pickle',original['inputPackets'][str(path)])
    for path in (Path(__file__).resolve(),PLAN):
        name=str(path.relative_to(f.ROOT));value=f.digest(path);f.require(name not in manifest['sourceFiles'],'fresh end binding source')
        manifest['sourceFiles'][name]=value;target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return cp,manifest


def native_adapter():
    """Exact prefix of the whole native end loop; threshold suffix stays pending."""
    original=ast.parse(textwrap.dedent(inspect.getsource(h.source.q.end_sources))).body[0]
    loop=next(n for n in original.body if isinstance(n,ast.For))
    split=next(i for i,n in enumerate(loop.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='zero' for t in n.targets))
    prefix=copy.deepcopy(loop.body[:split]);suffix=copy.deepcopy(loop.body[split:])
    routes={
        "modes.analytic(uniform['records'][label]['strong'])":"router.algebraic(uniform['records'][label]['strong'])",
        'inp.mapping(algebraic, relation, (r.omega, modes.k, modes.q, *inp.origin))':'router.mapping(algebraic, relation, (r.omega, modes.k, modes.q, *inp.origin))',
        'algebraic.xreplace(mapping).subs(origin).xreplace({r.omega: frequency})':"router.bind('PENCIL_PLUS', algebraic, mapping, origin, r.omega, frequency)",
        'relation.xreplace(mapping).subs(origin).xreplace({r.omega: frequency})':'router.wave(relation, mapping, origin, r.omega, frequency)',
        "(live.subs(frequency, inp.parameters['omega']) - old['freshPencil']).applyfunc(sp.cancel)":"router.seed_matrix(live, frequency, inp.parameters['omega'], old['freshPencil'])",
        "sp.cancel(wave.subs(frequency, inp.parameters['omega']) - old['curve'])":"router.seed_wave(wave, frequency, inp.parameters['omega'], old['curve'])",
        '-sp.diff(wave, frequency) / sp.diff(wave, modes.q)':'router.transport(wave, frequency, modes.q)',
        'live.diff(frequency) + live.diff(modes.q) * transport':'router.derivative(live, frequency, modes.q, transport)',
        'sp.cancel(wave.diff(frequency) + wave.diff(modes.q) * transport)':'router.tangency(wave, frequency, modes.q, transport)',
    }
    counts={key:0 for key in routes}
    class Forward(ast.NodeTransformer):
        def visit(self,node):
            if isinstance(node,ast.expr):
                key=ast.unparse(node)
                if key in routes:counts[key]+=1;return ast.copy_location(ast.parse(routes[key],mode='eval').body,node)
            return super().visit(node)
    transformed=[Forward().visit(n) for n in copy.deepcopy(prefix)]
    f.require(set(counts.values())=={1},('whole native prefix routing counts',counts))
    reverse={v:k for k,v in routes.items()}
    class Reverse(ast.NodeTransformer):
        def visit(self,node):
            if isinstance(node,ast.expr) and ast.unparse(node) in reverse:return ast.copy_location(ast.parse(reverse[ast.unparse(node)],mode='eval').body,node)
            return super().visit(node)
    restored=[Reverse().visit(n) for n in copy.deepcopy(transformed)]
    f.require(ast.dump(ast.Module(body=restored+suffix,type_ignores=[]))==ast.dump(ast.Module(body=loop.body,type_ignores=[])),'whole native loop reverse AST including untouched threshold suffix')
    whole=copy.deepcopy(original);whole_loop=next(n for n in whole.body if isinstance(n,ast.For));whole_loop.body=restored+suffix
    f.require(ast.dump(whole)==ast.dump(original),'entire original end_sources reverse AST')
    skeleton=ast.parse('def native_prefix(base,data,r,frequency,inp,uniform,modes,label,router):\n    pass\n').body[0]
    skeleton.body=transformed+[ast.Return(ast.Name(id='packet',ctx=ast.Load()))]
    module=ast.fix_missing_locations(ast.Module(body=[skeleton],type_ignores=[]));namespace={'sp':sp,'f':f}
    exec(compile(module,'<native saved-routed end prefix>','exec'),namespace)
    return namespace['native_prefix'],{'wholeNativeAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'routes':routes,'counts':counts,'prefixStatements':split,'pendingThresholdStatements':len(suffix),
        'entireOriginalBodyReverseAST':True,'unchangedNativeSeedBranchFreeSymbolGuards':True,
        'thresholdSuffixExecuted':False,'adapterSource':ast.unparse(module)}


def joined_sources(base):
    native={'endSources':h.source.q.end_sources,'endTables':h.inputs.chart.end_tables,
        'modalPrepare':engine.ModalCurrentSubspaces.prepare,'pairing':engine.ClosedCurrentPairing.construct,
        'mapping':engine.ChannelInput.mapping,'analytic':engine.FullPencilModes.analytic,
        'coordinateConstructor':previous.construct,'coordinateOperation':previous.atomic_operation,
        'originalModeSeed':previous.previous.modes_reader.joined_roots}
    result={name:previous.source_body(value) for name,value in native.items()}
    old=json.loads((base/'native-end-coordinate-joins.json').read_text())
    for name in ('endSources','modalPrepare','pairing','mapping','analytic'):same(result[name],old[name])
    function,adapter=native_adapter();result['adapter']=adapter
    f.save(base/'native-end-binding-joins.json',result)
    return function,result


def forbidden(*args,**kwargs):raise RuntimeError('live end binding may not repeat accepted constructors or derivatives')


def prohibit():
    # Seed cancellation is the only native normalizer needed for a new pair.
    cancel=sp.cancel;previous.prohibit();sp.cancel=cancel
    for module,names in ((previous,('load','construct','atomic_operation','main')),
                         (h.source.q,('load','end_sources','sources','main')),
                         (h.inputs.chart,('load','end_tables','source_chart','main'))):
        for name in names:
            if hasattr(module,name):setattr(module,name,forbidden)
    for name in ('solve','inv','svd','eig','eigh','lstsq','pinv'):setattr(previous.previous.modes_reader.np.linalg,name,forbidden)


def matrix_units(frame,extra):
    rows,fields=frame['equationUnits'],frame['fieldUnits']
    f.require(len(rows)==len(fields)==5 and all(len(v)==3 for v in (*rows,*fields)),'full actual end field/equation unit basis')
    return tuple(tuple(tuple(a-b+c for a,b,c in zip(row,field,extra)) for field in fields) for row in rows)


class Router:
    def __init__(self,base,coordinate,baseline,seed):
        self.base,self.coordinate,self.baseline,self.seed=base,coordinate,baseline,seed
        self.own=coordinate['sourceInput'];self.raw=self.own['rawInput'];self.ctx=self.own['context']
        self.call=self.own['sourceInput']['input'];self.images=coordinate['coordinateImages'];self.map=coordinate['mapping']
        self.omega=self.call['analyticInput']['omega'];self.frequency=self.call['frequency'];self.origin=self.ctx['seedInput']['origin']
        self.units={'PENCIL_PLUS':matrix_units(self.ctx['frame'],(0,0,0)),
            'FREQUENCY_PENCIL_PLUS':matrix_units(self.ctx['frame'],(0,1,0)),
            'NORMAL_PENCIL_PLUS':matrix_units(self.ctx['frame'],(1,0,0)),
            'FREQUENCY_TRANSPORT':(0,0,0),'NORMAL_TRANSPORT':(1,-1,0),'wave':(0,-2,0)}
        self.atlas=[];self.receipts=[];self.bound={}
        for end in ('LEFT','RIGHT'):
            item=f.unpickle(base/'end-call-baseline'/end.lower()/'input-value.pickle');p=item['result']['frequencyPencil']
            actual_units=item['analyticInput']['units']
            f.require(tuple(actual_units['radical'])==(0,-1,0),'actual accepted wave unit declaration')
            for name,expr,key,unit in (('PENCIL_PLUS',p['originalAlgebraic'],'livePencil',matrix_units(actual_units,(0,0,0))),
                                       ('wave',p['originalRelation'],'wave',(0,-2,0))):
                inp=self.arguments(expr,p['mapping'],p['origin'],item['analyticInput']['omega'],p['frequency'],unit)
                self.atlas.append((inp,p[key],{'type':'accepted-original-native-binding','end':end,'field':key,
                    'path':str(base/'end-call-baseline'/end.lower()/'input-value.pickle')}))
        # The complete wave call, seed curve and transport are already saved.
        pair={'requested':self.arguments(coordinate['originalRelation'],self.map,self.origin,self.omega,self.frequency,self.units['wave']),
            'accepted':self.arguments(baseline['originalRelation'],baseline['mapping'],baseline['origin'],self.omega,baseline['frequency'],self.units['wave']),
            'coordinate':coordinate,'baseline':baseline,'acceptedOriginalSeed':seed}
        f.atomic_pickle(base/'wave-and-seed-route-inputs.pickle',pair)
        same(pair['requested'],pair['accepted']);same(coordinate['branchResiduals'],baseline['branchResiduals'])
        same(self.raw['modeInput']['relation'],baseline['originalRelation'])
        same(self.ctx['seedInput']['parameters']['omega'],self.raw['modeInput']['frequency'])
        # The existing acoustic source proof establishes this full curve caller.
        proof=self.own['sourceProofReceipt'];f.require(proof['relationMatchesBaseline'] and proof['branchReceiptMatchesBaseline'],'accepted own acoustic caller proof')
        for atom in baseline['originalRelation'].free_symbols:
            if atom in (self.call['analyticInput']['momentum'],self.call['analyticInput']['radical']):continue
            expected=self.ctx['seedInput']['parameters']['omega'] if atom==self.omega else self.map[atom]
            f.require(self.raw['modeInput']['binding'][atom]==expected,'actual own seed relation binding')
        f.save(base/'wave-and-seed-route-checks.json',{'fullNativeWaveCallJoined':True,'ownSeedRelationBindingJoined':True,
            'acceptedWaveAndSeedCurveReused':True,'newWaveOrAcousticDerivativeCalls':0})

    @staticmethod
    def arguments(expression,mapping,origin,omega,frequency,unit):
        return {'expression':expression,'mapping':mapping,'origin':origin,'omega':omega,'frequency':frequency,'unit':unit}

    def algebraic(self,strong):
        same(strong,self.raw['signature']['strong'])
        return self.images['PENCIL_PLUS'],self.coordinate['originalRelation'],self.coordinate['branchResiduals']

    def mapping(self,algebraic,relation,live):
        saved=f.unpickle(self.base/'saved-native-mapping-input-return.pickle')
        same(algebraic,saved['algebraic']);same(relation,saved['relation']);same(live,saved['live'])
        return self.map

    def bind(self,name,expression,mapping,origin,omega,frequency):
        op=self.arguments(expression,mapping,origin,omega,frequency,self.units[name])
        folder=self.base/'binding-operations'/name;folder.mkdir(parents=True)
        f.atomic_pickle(folder/'input.pickle',op)
        matches=[(v,owner) for key,v,owner in self.atlas if previous.previous.modes_reader.same(key,op)]
        if matches:
            value,owner=matches[0]
            for other,_ in matches[1:]:same(value,other)
            receipt={'name':name,'reused':True,'owner':owner,'completeActualNativeArguments':True}
        else:
            # Literal entire native binding chain, with immediate intermediate
            # checkpoints. Each call exists only for a genuinely missing input.
            mapped=expression.xreplace(mapping);f.atomic_pickle(folder/'mapped.pickle',mapped)
            retained=mapped.subs(origin);f.atomic_pickle(folder/'retained.pickle',retained)
            value=retained.xreplace({omega:frequency});f.atomic_pickle(folder/'value.pickle',value)
            receipt={'name':name,'reused':False,'owner':{'type':'new-native-binding','name':name},'completeActualNativeArguments':True}
        if matches:f.atomic_pickle(folder/'value.pickle',value)
        receipt.update(inputSha256=f.digest(folder/'input.pickle'),valueSha256=f.digest(folder/'value.pickle'))
        f.save(folder/'completed.json',receipt);self.receipts.append(receipt);self.bound[name]=value
        self.atlas.append((op,value,receipt['owner']))
        return value

    def wave(self,relation,mapping,origin,omega,frequency):
        return self.bind('wave',relation,mapping,origin,omega,frequency)

    def saved_binding(self,name):
        return self.bind(name,self.images[name],self.map,self.origin,self.omega,self.frequency)

    def seed_matrix(self,live,frequency,seed_frequency,physical):
        path=self.base/'seed-matrix-join';path.mkdir()
        f.atomic_pickle(path/'input.pickle',{'live':live,'frequency':frequency,'seed':seed_frequency,'ownPhysical':physical,'owner':OWNER,'ownRawInput':self.raw})
        actual=live.subs(frequency,seed_frequency);f.atomic_pickle(path/'seed-value.pickle',actual)
        raw=actual-physical;f.atomic_pickle(path/'raw-residual.pickle',raw)
        # Whole native applyfunc(cancel), checkpoint each actual scalar call.
        calls=[];values=[]
        def cancel(value):
            index=len(calls);folder=path/str(index);folder.mkdir()
            f.atomic_pickle(folder/'input.pickle',value)
            matches=[(i,result) for i,(argument,result) in enumerate(values) if previous.previous.modes_reader.same(argument,value)]
            result=matches[0][1] if matches else sp.cancel(value)
            f.atomic_pickle(folder/'value.pickle',result)
            calls.append({'index':index,'inputSha256':f.digest(folder/'input.pickle'),'valueSha256':f.digest(folder/'value.pickle'),
                'reusedCompletedScalar':bool(matches),'firstOwner':matches[0][0] if matches else index})
            values.append((value,result))
            f.save(folder/'completed.json',calls[-1]);return result
        residual=raw.applyfunc(cancel);f.atomic_pickle(path/'normalized-residual.pickle',residual)
        f.save(path/'completed.json',{'nativeCancelCalls':calls,'ownSeedPhysicalSourceJoined':True})
        return residual

    def seed_wave(self,wave,frequency,seed_frequency,curve):
        same(wave,self.baseline['wave']);same(frequency,self.baseline['frequency']);same(curve,self.seed['curve'])
        same(seed_frequency,self.ctx['seedInput']['parameters']['omega'])
        f.atomic_pickle(self.base/'saved-wave-seed-join.pickle',{'wave':wave,'frequency':frequency,'seedFrequency':seed_frequency,'curve':curve,
            'actualOwnRelation':self.raw['modeInput']['relation'],'ownBinding':self.raw['modeInput']['binding'],
            'acceptedResidual':self.baseline['referenceWaveResidual']})
        return self.baseline['referenceWaveResidual']

    def transport(self,wave,frequency,radical):
        same(wave,self.baseline['wave']);same(frequency,self.baseline['frequency']);same(radical,self.baseline['radical'])
        value=self.saved_binding('FREQUENCY_TRANSPORT')
        f.atomic_pickle(self.base/'saved-transport-native-join.pickle',{'wave':wave,'frequency':frequency,'radical':radical,
            'boundCoordinateTransport':value,'acceptedNativeTransport':self.baseline['radicalFrequencyTransport'],
            'rawCoordinateTransport':self.images['FREQUENCY_TRANSPORT']})
        same(value,self.baseline['radicalFrequencyTransport'])
        return value

    def derivative(self,live,frequency,radical,transport):
        same(live,self.bound['PENCIL_PLUS']);same(frequency,self.frequency);same(radical,self.baseline['radical']);same(transport,self.bound['FREQUENCY_TRANSPORT'])
        value=self.saved_binding('FREQUENCY_PENCIL_PLUS')
        f.atomic_pickle(self.base/'saved-tangent-native-join.pickle',{'livePencil':live,'wave':self.baseline['wave'],'frequency':frequency,'radical':radical,
            'transport':transport,'coordinateTangent':self.images['FREQUENCY_PENCIL_PLUS'],'boundTangent':value,
            'source':self.own,'coordinateProofPath':str(self.base/'operations/FREQUENCY_PENCIL_PLUS/roundtrip.pickle'),
            'justification':'Entire native ModalCurrentSubspaces.prepare tangent, exact accepted constant-scale coordinate conversion and native material/grade binding. No new differentiation.'})
        return value

    def tangency(self,wave,frequency,radical,transport):
        same(wave,self.baseline['wave']);same(frequency,self.baseline['frequency']);same(radical,self.baseline['radical']);same(transport,self.baseline['radicalFrequencyTransport'])
        f.atomic_pickle(self.base/'saved-tangency-native-join.pickle',{'wave':wave,'frequency':frequency,'radical':radical,'transport':transport,'residual':self.baseline['waveTangencyResidual']})
        return self.baseline['waveTangencyResidual']


def construct(base,cp,native):
    coordinate=f.unpickle(base/'remaining-case-end-coordinate-images.pickle');own=coordinate['sourceInput'];raw=own['rawInput'];context=own['context']
    original=f.unpickle(base/'end-call-baseline/right/input-value.pickle');baseline=original['result']['frequencyPencil'];seed=f.unpickle(base/'native-frequency-seeds/right-symbolic.pickle')
    router=Router(base,coordinate,baseline,seed)
    inp=SimpleNamespace(**context['seedInput']);call=own['sourceInput']['input']['analyticInput']
    r=SimpleNamespace(omega=call['omega']);modes=SimpleNamespace(k=call['momentum'],q=call['radical'])
    data={'acceptedEnds':{'RIGHT':{'freshPencil':raw['modeInput']['physical'],'curve':seed['curve']}}}
    uniform={'records':{'RIGHT':{'strong':raw['signature']['strong']}}}
    requested={'owner':OWNER,'coordinate':coordinate,'actualSeed':data['acceptedEnds']['RIGHT'],
        'savedAcousticSeed':seed,'nativeCall':{'omega':r.omega,'frequency':router.frequency,'momentum':modes.k,'radical':modes.q,
            'seedInput':context['seedInput'],'uniform':uniform},'fieldUnits':context['frame']['fieldUnits'],'equationUnits':context['frame']['equationUnits']}
    f.atomic_pickle(base/'native-prefix-input.pickle',requested)
    target=base/'new-end';target.mkdir()
    packet=native(target,data,r,router.frequency,inp,uniform,modes,'RIGHT',router)
    f.save(base/'native-prefix-checks.json',{'branchScalars':len(packet['branchResiduals']),'seedMatrixScalars':len(packet['referencePencilResidual']),
        'allNativeGuardsPassed':True,'newDerivativeCalls':0,'thresholdSuffixPending':True})
    for name in ('NORMAL_TRANSPORT','NORMAL_PENCIL_PLUS'):router.saved_binding(name)
    # The second physical anchoring consumes exactly the accepted shared inputs.
    aliases={}
    for label,case in coordinate['aliases'].items():
        aliases[label]={}
        for end,old in case.items():
            route={'address':(label,end),'coordinateAlias':old,'fullInputPath':str(base/'coordinate-cases'/label/end.lower()/'full-input-route.pickle'),
                'pencilPath':str(target/'right-frequency-pencil.pickle') if old['mode']=='missing-coordinate-image' else str(base/'end-accepted/frequency'/(end.lower()+'-frequency-pencil.pickle')),
                'completeEndDomainOrNumericalReuseAccepted':False}
            if old['mode']=='missing-coordinate-image':route['owner']=OWNER
            else:route['owner']=old['owner']
            folder=base/'binding-cases'/label/end.lower();folder.mkdir(parents=True);f.save(folder/'route.json',route);aliases[label][end]=route
        f.save(base/'binding-cases'/label/'case-summary.json',aliases[label])
    unit=router.units['FREQUENCY_PENCIL_PLUS'];changed_unit=tuple(tuple((v[0]+1,v[1],v[2]) if (i,j)==(0,0) else v for j,v in enumerate(row)) for i,row in enumerate(unit))
    changed=sp.MutableDenseMatrix(packet['livePencil']);changed[0,0]+=1;changed=sp.ImmutableMatrix(changed)
    mapping=dict(router.map);atom=next(iter(mapping));mapping[atom]=mapping[atom]+1
    mutations={'actualPencil':packet['livePencil'],'changedPencil':changed,'actualMapping':router.map,'changedMapping':mapping,
        'actualOwner':OWNER,'changedOwner':(OWNER[0],'LEFT'),'actualDerivativeUnit':unit,'changedDerivativeUnit':changed_unit}
    f.atomic_pickle(base/'binding-mutation-operands.pickle',mutations)
    equal=previous.previous.modes_reader.same
    controls={'pencilCoefficient':not equal(changed,packet['livePencil']),'mappingCoefficient':not equal(mapping,router.map),
        'owner':OWNER!=mutations['changedOwner'],'derivativeUnit':not equal(unit,changed_unit)}
    f.save(base/'binding-mutation-controls.json',controls);f.require(all(controls.values()),'actual own live binding controls')
    output={'owner':OWNER,'pencil':packet,'boundCoordinateImages':router.bound,'coordinateInput':coordinate,
        'bindingReceipts':router.receipts,'aliases':aliases,'units':router.units,'thresholdsTablesAndDomainsPending':True}
    f.atomic_pickle(base/'remaining-case-end-bindings.pickle',output)
    return {'cases':aliases,'physicalEnds':8,'newEndOwners':1,'sharedNewUses':2,'baselineAliases':6,
        'newBindingCalls':sum(not x['reused'] for x in router.receipts),'reusedBindingCalls':sum(x['reused'] for x in router.receipts),
        'newDerivativeAnalyticMappingCurrentModeCalls':0,'actualMutationControls':len(controls),'thresholdsTablesAndDomainsPending':True,
        'completeFrequencyEndOrNumericalReuseAccepted':False}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    previous.previous.protect(base);cp,manifest=load(base);native,joins=joined_sources(base);prohibit();result=construct(base,cp,native)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen source posthash')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal'] and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'reference postidentity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_END_BINDINGS','nativeJoins':joins,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
