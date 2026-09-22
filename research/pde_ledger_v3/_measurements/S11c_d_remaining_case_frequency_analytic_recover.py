#!/usr/bin/env python3
"""Resume unstarted analytic images with exact native derivative-call routing."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil
import textwrap

import S11c_d_remaining_case_frequency_analytic as h

f,n,sp=h.f,h.n,h.sp
PREVIOUS=Path('/var/projects/toy_physics/_scratch/s11c/s11c-remaining-case-frequency-20260921/analytic-sources')
DIAGNOSTIC=PREVIOUS.parent/'analytic-derivative-routing'
PLAN=f.M/'S11c_d_remaining_case_frequency_analytic_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_frequency_analytic_derivative_repair.json'
OriginalRouter=h.ProofRouter
original_protect=h.inputs.protect_references
STATE={}


def protect(base):
    original_protect(base)
    protected_save=f.save
    def save(path,value):
        path=Path(path)
        if path==base/'native-analytic-constructor-joins.json':
            old=json.loads(path.read_text())
            for key in ('source','lift','denominator','nativeSeedComparisonSha256',
                        'nativeBindingComparisonSha256','nativeCertificateSha256','nativeSeedRootSha256'):
                f.require(value[key]==old[key],('unchanged scientific constructor join',key))
            target=base/'recovery-native-analytic-constructor-joins.json'
            f.require(not target.exists(),'fresh recovery native join')
            protected_save(target,value)
            protected_save(base/'recovery-native-join-write-route.json',{
                'originalPath':str(path),'originalSha256':f.digest(path),'newPath':str(target),
                'originalPreserved':True,'onlyChangedJoins':'exact saved derivative route body'})
            return
        return protected_save(path,value)
    f.save=save


def restore_seed_prefix(router,baseline):
    base=router.base
    atlas=base/'saved-seed-root-atlas.pickle';pairs=base/'seed-root-input-pairs.pickle'
    repair=STATE['repair']
    for path in (atlas,pairs):
        name=str(path.relative_to(STATE['base']))
        f.require(f.digest(path)==repair['completedFiles'][name]['sha256'],'completed seed prefix hash')
    # The original completed substitution/equality guards are reused through
    # their complete input/result bytes and the diagnostic's clean hash audit.
    router.seed=f.unpickle(atlas);raw=f.unpickle(pairs)
    accepted_roots=f.unpickle(STATE['base']/'accepted/accepted-chart/root-chart.pickle')
    accepted_baseline=STATE['base']/'accepted/accepted-chart/analytic-sources.pickle'
    f.require(n.same(router.roots,accepted_roots),'same actual root input for completed seed prefix')
    f.require(f.digest(accepted_baseline)==repair['baselineAnalyticSha256'],'same actual baseline analytic records')
    f.require(len(raw)==repair['completedSeedInputCount']==260 and len(router.seed)==6,'complete saved seed prefix')
    for index,item in enumerate(raw):
        saved=f.unpickle(base/'seed-root-inputs'/str(index)/'input.pickle')
        pair=baseline[item['sourceKey']]['bindingComparisons'][item['point']]
        f.require(n.same(saved,{'sourceKey':item['sourceKey'],'point':item['point'],
                              'join':item['join'],'roots':router.roots}), 'literal completed seed input')
        f.require(any(n.same(item['join'],entry) for entry in pair['rootIdentities']),
                  'actual saved producer comparison member')
        entry=router.seed[item['arguments']]
        f.require(n.same(entry['value'],item['join']['rootIdentity']),'completed seed result route')
    f.save(STATE['base']/'completed-seed-prefix-reuse.json',{
        'inputs':len(raw),'identities':len(router.seed),'original':str(PREVIOUS/'complete/analytic-proof-reuse'),
        'rawPairsSha256':f.digest(pairs),'atlasSha256':f.digest(atlas),
        'diagnosticChecksSha256':repair['diagnosticChecksSha256'],
        'substitutionOrProofGuardRepeated':False,'fullSavedInputJoined':True})


def initializer_adapter():
    original=ast.parse(textwrap.dedent(inspect.getsource(OriginalRouter.__init__)))
    changed=copy.deepcopy(original);fn=changed.body[0]
    raw_index=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='raw')
    loop_index=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.For) and ast.unparse(x.iter)=='baseline.items()')
    loop=fn.body[loop_index]
    pair_loop=next(x for x in loop.body if isinstance(x,ast.For) and ast.unparse(x.iter)=="record['bindingComparisons'].items()")
    seed_index=next(i for i,x in enumerate(pair_loop.body) if isinstance(x,ast.For) and ast.unparse(x.iter)=="pair['rootIdentities']")
    seed_loop=copy.deepcopy(pair_loop.body.pop(seed_index))
    tail=copy.deepcopy(fn.body[loop_index+1:]);raw_stmt=copy.deepcopy(fn.body[raw_index])
    fn.body[raw_index]=ast.parse('restore_seed_prefix(self,baseline)').body[0]
    del fn.body[loop_index+1:]
    reverse=copy.deepcopy(changed);rf=reverse.body[0]
    rf.body[raw_index]=raw_stmt;rf.body.extend(tail)
    rp=next(x for x in rf.body[loop_index].body if isinstance(x,ast.For) and ast.unparse(x.iter)=="record['bindingComparisons'].items()")
    rp.body.insert(seed_index,seed_loop)
    f.require(ast.dump(reverse)==ast.dump(original),'whole original router initializer reverse AST: completed seed prefix only')
    ns=dict(vars(h),restore_seed_prefix=restore_seed_prefix)
    exec(compile(ast.fix_missing_locations(changed),'<saved-seed-prefix-initializer>','exec'),ns)
    return ns['__init__'],{'wholeInitializerReverseAST':True,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'removedOnlyCompletedSeedInputProofAndWrites':True}


class CallRouter(OriginalRouter):
    def add_call(self,expression,unit,variables,value,owner):
        unit=tuple(unit);variables=tuple(variables)
        key=(expression,unit,variables);existing=self.derivatives.get(key)
        if existing is not None:
            f.require(n.same(existing['value'],value),'same actual native call has identical saved result')
        else:
            self.derivatives[key]={'expression':expression,'unit':unit,'variables':variables,
                                   'value':value,'owner':owner}

    def add_derivative(self,expression,unit,order,value,owner):
        # This legacy entrypoint is used only by the original analytic initializer
        # and its direct native diff(analytic,w[,2]) calls.
        f.require(owner['kind']=='accepted-analytic-derivative' and order in (1,2),'actual analytic producer call')
        variables=(self.roots['frequency'],) if order==1 else (self.roots['frequency'],2)
        self.add_call(expression,unit,variables,value,dict(owner,nativeCall='sp.diff(analytic,w)' if order==1 else 'sp.diff(analytic,w,2)'))

    def add_source_derivatives(self,record,label,key,packet):
        w=self.roots['frequency'];unit=tuple(record['unit'])
        f.require(packet['frequency']==w and tuple(packet['dimensionState']['known'][w])==(0,-1,0),
                  'actual frequency variable and input unit for native chained derivative')
        owner={'kind':'accepted-live-source-derivative','case':label,'key':key,
               'address':record['address'],'sourceOwner':record['owner'],'ownerAddress':record['ownerAddress']}
        self.add_call(record['liveFrequencyAndGrades'],unit,(w,),record['firstFrequencyDerivative'],
                      dict(owner,field='firstFrequencyDerivative',nativeCall='sp.diff(live,w)'))
        self.add_call(record['firstFrequencyDerivative'],(unit[0],unit[1]+1,unit[2]),(w,),record['secondFrequencyDerivative'],
                      dict(owner,field='secondFrequencyDerivative',nativeCall='sp.diff(first,w)'))

    def derivative(self,expression,*variables,**kwargs):
        f.require(not kwargs and variables in ((self.roots['frequency'],),(self.roots['frequency'],2)),
                  'only exact native first/second frequency derivative call shapes')
        unit=tuple(self.recorder.unit);key=(expression,unit,tuple(variables))
        target=self.base/'derivative-calls'/str(len(self.derivative_calls));target.mkdir(parents=True,exist_ok=False)
        f.atomic_pickle(target/'input.pickle',{'expression':expression,'variables':variables,'unit':unit,
                       'active':str(self.recorder.active),'routing':'literal actual call arguments'})
        entry=self.derivatives.get(key)
        if entry is None:
            value=sp.diff(expression,*variables)
            entry={'expression':expression,'unit':unit,'variables':variables,'value':value,
                   'owner':{'kind':'new-native-analytic-derivative','path':str(target)}}
            f.atomic_pickle(target/'value.pickle',entry)
            self.add_call(expression,unit,variables,value,entry['owner']);reused=False
        else:
            f.require(n.same(key,(entry['expression'],tuple(entry['unit']),tuple(entry['variables']))),
                      'literal native argument/unit/variable tuple result route')
            f.atomic_pickle(target/'value.pickle',entry);reused=True
        self.derivative_calls.append({'path':str(target),'reused':reused,'owner':entry['owner']})
        return entry['value']


def constructor_adapter():
    original=ast.parse(inspect.getsource(h.construct));changed=copy.deepcopy(original)
    mkdir=next(x for x in ast.walk(changed) if isinstance(x,ast.Call) and ast.unparse(x.func)=='proof_base.mkdir')
    mkdir.keywords.append(ast.keyword(arg='exist_ok',value=ast.Constant(True)))
    loop=next(x for x in ast.walk(changed) if isinstance(x,ast.For) and ast.unparse(x.iter)=="packet['records'].items()")
    pos=next(i for i,x in enumerate(loop.body) if isinstance(x,ast.For) and ast.unparse(x.target)=='(order, field)')
    old=copy.deepcopy(loop.body[pos]);loop.body[pos]=ast.parse('router.add_source_derivatives(record,label,key,packet)').body[0]
    reverse=copy.deepcopy(changed)
    rm=next(x for x in ast.walk(reverse) if isinstance(x,ast.Call) and ast.unparse(x.func)=='proof_base.mkdir');rm.keywords=[]
    rl=next(x for x in ast.walk(reverse) if isinstance(x,ast.For) and ast.unparse(x.iter)=="packet['records'].items()")
    rl.body[pos]=old
    f.require(ast.dump(reverse)==ast.dump(original),'whole constructor reverse AST: source actual-call route and saved directory only')
    ns=dict(vars(h),ProofRouter=CallRouter)
    exec(compile(ast.fix_missing_locations(changed),'<native-analytic-call-route-recovery>','exec'),ns)
    return ns['construct'],{'wholeConstructorReverseAST':True,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'sourceDerivativeRegistrationReplacements':1,'existingDirectoryKeywords':1,
        'allNativeScientificBodiesUnchanged':True}


def diagnostic_input_join(base):
    evidence=base/'accepted-derivative-diagnostic'
    checks=json.loads((evidence/'checks.json').read_text())
    f.require(checks['literalConflicts']==0 and checks['coarseConflicts']==checks['coarseConflictsWithDifferentActualCalls']==372,
              'diagnostic resolves only distinct actual call collisions')
    summaries=json.loads((evidence/'coarse-conflicts.json').read_text())
    for item in summaries:
        pair=f.unpickle(evidence/('coarse-conflict-'+str(item['index'])+'.pickle'))
        left,right=pair['existing'],pair['candidate']
        f.require(not pair['sameSavedValues'] and not pair['sameNativeCallArguments'], 'saved diagnostic collision disposition')
        f.require(not n.same((left['argument'],left['variables'],left['argumentUnit']),
                             (right['argument'],right['variables'],right['argumentUnit'])), 'actual raw call tuple differs')
        f.require(left['formalOrder']==right['formalOrder']==2,'actual failed formal-order collision')
    first=f.unpickle(evidence/'coarse-conflict-0.pickle')['existing']
    original=(first['argument'],tuple(first['argumentUnit']),tuple(first['variables']))
    mutations={'argument':(first['argument']+1,original[1],original[2]),
               'unit':(original[0],(original[1][0],original[1][1]+1,original[1][2]),original[2]),
               'variables':(original[0],original[1],(original[2][0],3))}
    f.atomic_pickle(base/'derivative-route-mutation-operands.pickle',{'savedInput':first,'originalCall':original,'mutations':mutations})
    controls={}
    for name,changed in mutations.items():
        try:f.require(n.same(changed,original),'exact saved derivative caller input')
        except ValueError:controls[name]=True
        else:controls[name]=False
    f.save(base/'derivative-route-mutation-controls.json',controls)
    f.require(all(controls.values()),'actual argument/unit/order routing controls reject')
    f.save(base/'derivative-diagnostic-reuse.json',{'pairs':len(summaries),'literalConflicts':0,
        'checksSha256':f.digest(evidence/'checks.json'),'noNewDerivativeOrEquivalenceTest':True})


def load(base):
    repair=json.loads(REPAIR.read_text());STATE.update(base=base,repair=repair)
    f.require(f.digest(Path(__file__))==repair['newHelperSha256'] and f.digest(PLAN)==repair['newPlanSha256'],'reviewed recovery source')
    oldbase=PREVIOUS/'complete';old=json.loads((oldbase/'inputs.json').read_text())
    f.require(f.digest(oldbase/'inputs.json')==repair['originalInputsSha256'] and f.digest(Path(h.__file__))==repair['originalHelperSha256'],
              'immutable original helper and completed input loader')
    f.require(h.body(h.main)==repair['originalMainSha256'],'whole original main unchanged')
    receipt=h.source.receipts.inspect_guard(DIAGNOSTIC,'diagnose')
    f.require(f.digest(DIAGNOSTIC/'checks.json')==repair['diagnosticChecksSha256']
              and (DIAGNOSTIC/'checks.json').read_bytes()==(DIAGNOSTIC/'diagnose.stdout').read_bytes(), 'actual final clean saved routing diagnostic')
    manifest=dict(old,runDirectory=str(base),sourceFiles=dict(old['sourceFiles']),inputPackets=dict(old['inputPackets']),referencedInputs={})
    actual={str(p.relative_to(oldbase)) for p in oldbase.rglob('*') if p.is_file()
            and 'source' not in p.relative_to(oldbase).parts and p!=oldbase/'inputs.json'}
    f.require(actual==set(repair['completedFiles']),'all and only completed failed constructor files')
    for name,item in repair['completedFiles'].items():
        path=oldbase/name
        f.require(f.digest(path)==item['sha256'] and path.stat().st_size==item['bytes']
                  and str(path.resolve())==item['resolved'] and (str(path.readlink()) if path.is_symlink() else None)==item['rawLink'],
                  ('completed original bytes/address',name))
        h.source.reference(base,manifest,path,name,item['sha256'])
    h.source.reference(base,manifest,oldbase/'inputs.json','original-analytic-inputs.json',repair['originalInputsSha256'])
    for root,prefix,inventory in ((PREVIOUS,'original-analytic-logs',repair['originalLogs']),
                                  (DIAGNOSTIC,'accepted-derivative-diagnostic',repair['diagnosticFiles'])):
        for name,item in inventory.items():h.source.reference(base,manifest,root/name,prefix+'/'+name,item['sha256'])
    for path in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(path.relative_to(f.ROOT));value=f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'new distinct source pins')
        manifest['sourceFiles'][name]=value
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,('current source pin',name))
        if name in old['sourceFiles']:f.require(f.digest(oldbase/'source'/name)==value,'unchanged original frozen source')
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest)
        f.require(f.digest(dest)==value,'new frozen source')
    for path,value in manifest['inputPackets'].items():f.require(f.digest(Path(path))==value,'all original/reference input prehashes')
    diagnostic_input_join(base)
    manifest['derivativeCallRecovery']={'original':str(PREVIOUS),'originalInputsSha256':repair['originalInputsSha256'],
        'diagnostic':str(DIAGNOSTIC),'diagnosticChecksSha256':repair['diagnosticChecksSha256'],'diagnosticGuard':receipt,
        'completedFilesReused':len(repair['completedFiles']),'completedSeedInputs':repair['completedSeedInputCount'],
        'actualCallRouting':True,'initializerJoin':STATE['initializerJoin'],'constructorJoin':STATE['constructorJoin'],
        'originalMainSha256':repair['originalMainSha256']}
    f.save(base/'inputs.json',manifest)
    f.save(base/'completed-analytic-reference-reuse.json',manifest['referencedInputs'])
    f.save(base/'derivative-call-recovery-join.json',manifest['derivativeCallRecovery'])
    return manifest,tuple(repair['cases'])


if __name__=='__main__':
    initialize,ij=initializer_adapter();CallRouter.__init__=initialize
    construct,cj=constructor_adapter();STATE.update(initializerJoin=ij,constructorJoin=cj)
    repair=json.loads(REPAIR.read_text())
    f.require(ij==repair['initializerJoin'] and cj==repair['constructorJoin'],'reviewed literal whole-body recovery adapters')
    h.inputs.protect_references=protect
    h.load=load;h.ProofRouter=CallRouter;h.construct=construct
    h.main()
