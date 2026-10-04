#!/usr/bin/env python3
"""Guarded saved-prefix continuation of the bounded contracted J/direct subset.

New algebra/restoration only after containment; no completed scientific function
or integral is replayed. Source/JSON inspection and synthetic tests are separate.
"""
import argparse,ast,hashlib,importlib.util,json,math,os,resource,shutil,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
from S11c_d_defect_packet_contracted_continue_resume import tree_inventory,check_tree
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO={'text':'0','srepr':'Integer(0)'}
sp=None

def require(value,message):
    if value is not True:raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def posthash(path,expected):
    try:
        actual=sha(path)
        return {'expected':expected,'actual':actual,'intact':actual==expected}
    except OSError as error:
        return {'expected':expected,'actual':None,'intact':False,'error':str(error)}

def read(path):return json.loads(Path(path).read_text())

def packed(value):
    if sp is not None and isinstance(value,sp.Basic):return {'text':str(value),'srepr':sp.srepr(value)}
    if type(value) is F:return str(value)
    if hasattr(value,'packed'):return value.packed()
    if isinstance(value,dict):
        require(all(type(k) is str for k in value),'string JSON keys')
        return {k:packed(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [packed(v) for v in value]
    # Old JSON receipts contain elapsed seconds. Preserve finite metadata floats;
    # native unit arithmetic separately requires exact rational exponents.
    if type(value) is float:
        require(math.isfinite(value),'finite inherited JSON metadata');return value
    require(type(value) in (type(None),str,int,bool),'supported JSON values only')
    return value

def save(path,value):
    with Path(path).open('x') as f:
        json.dump(packed(value),f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)

class Journal:
    def __init__(self,out):self.out=out;self.sequence=0;self.previous='0'*64;self.active=None;self.completed=[]
    def emit(self,name,value):
        p=self.out/(name+'.json');save(p,value)
        record={'sequence':self.sequence,'name':name,'sha256':sha(p),'bytes':p.stat().st_size,'previous':self.previous}
        digest=hashlib.sha256(canonical(record).encode()).hexdigest();record['chainSha256']=digest
        with (self.out/'evidence-chain.jsonl').open('a') as f:
            f.write(canonical(record)+'\n');f.flush();os.fsync(f.fileno())
        self.sequence+=1;self.previous=digest
        return record
    def start(self,name,args):
        require(self.active is None,'one active exact operation');self.active=name;self.emit(name+'-input',args)
    def finish(self,value):
        require(self.active is not None,'active exact operation');name=self.active;self.emit(name+'-return',value)
        self.completed.append(name);self.active=None

def copy_inputs(m,J):
    observed=tree_inventory(m['priorRoot']);J.emit('prior-tree-input',{'observed':observed,'expected':m['priorFiles'],'expectedCount':8271,'expectedBytes':158628313})
    check_tree(observed,m['priorFiles'])
    prior={}
    for relative,receipt in m['priorFiles'].items():
        source=Path(m['priorRoot'])/relative;dest=J.out/'prior'/relative;dest.parent.mkdir(parents=True,exist_ok=True)
        require(sha(source)==receipt['sha256'] and source.stat().st_size==receipt['bytes'],'immutable prior file '+relative)
        shutil.copyfile(source,dest);require(sha(dest)==receipt['sha256'],'identical prior copy '+relative)
        prior[relative]={'source':str(source),'path':str(dest.relative_to(J.out)),**receipt}
        J.emit('prior-copy-'+str(len(prior)),{'relative':relative,**prior[relative]})
    after=tree_inventory(m['priorRoot']);copied_tree=tree_inventory(J.out/'prior')
    J.emit('prior-tree-after-copy',{'source':after,'copy':copied_tree})
    check_tree(after,m['priorFiles']);check_tree(copied_tree,m['priorFiles'])
    J.emit('prior-copy-index',prior)
    raw={};copies={}
    for alias,receipt in m['savedInputs'].items():
        source=Path(receipt['path']);dest=J.out/'saved'/alias;dest.parent.mkdir(parents=True,exist_ok=True)
        require(sha(source)==receipt['sha256'] and source.stat().st_size==receipt['bytes'],'original input '+alias)
        shutil.copyfile(source,dest);require(sha(dest)==receipt['sha256'],'copied input '+alias)
        copies[alias]={'source':str(source),'path':str(dest.relative_to(J.out)),'sha256':receipt['sha256'],'bytes':receipt['bytes']}
        # A receipt per file survives even if a later parse or join refuses.
        J.emit('copy-'+str(len(copies)),{'alias':alias,**copies[alias]});raw[alias]=read(dest)
    J.emit('saved-copy-index',copies)
    return raw,copies

def verify_gate(path,manifest_path,m):
    g=read(path)
    require(g['status']=='READY_FOR_ONE_PACKET_CONTRACTED_CONTINUATION' and g['independentBuildClearance'] is True,'actual continuation build clearance')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'actual continuation worker/manifest')
    require(g['sourcePins']==m['sourcePins'],'complete continuation source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','library','authority','buildReviewRecord','continuationMethod','tailMethodRecord','priorResultRecord'):
        require(sha(g[key])==g[key+'Sha256'],'gate document '+key)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual unchanged guard and supervisor')
    require(g['launcher']==m['launcher'] and g['library']==m['librarySource'] and g['authority']==m['executionAuthority'] and g['buildReviewRecord']==m['reviewRecordWillBe'],'actual execution documents')
    require(g['continuationMethod']==m['continuationMethod'] and g['tailMethodRecord']==m['tailMethodRecord'] and g['priorResultRecord']==m['priorResultRecord'],'actual amended method and preserved failure')
    r=read(g['buildReviewRecord'])
    policy=read(m['reviewPolicy'])
    require(policy['status']=='ACTIVE_TEMPORARY_CLAUDE_ONLY_USER_POLICY' and policy['selectedReviewers']==g['reviewers']==m['reviewers']==['claude'] and g['reviewPolicySha256']==sha(m['reviewPolicy']),'actual explicit Claude-only policy')
    require(r['reviewers']==['claude'] and r['reviewPolicySha256']==g['reviewPolicySha256'] and r['sourcePins']==g['sourcePins'],'exact reviewed helper/test/log/policy pins')
    expectedVerdict='CLEAR FOR THIS SAVED-PREFIX CONTRACTED NUMERICAL CONTINUATION BUILD'
    require(g['literalBuildVerdicts']=={e:expectedVerdict for e in ('claude',)} and {e:r['reports'][e]['literalVerdict'] for e in ('claude',)}==g['literalBuildVerdicts'],'actual Claude literal build verdict under temporary policy')
    require(r['allChecksPassed'] is True and r['independentBuildClearance'] is True and r['amendedTailMethodAssessed'] is True,'independent concrete build and amended method assessment')
    for k in ('workerSha256','manifestSha256','librarySha256','launcherSha256','sharedGuardSha256','supervisorSha256','continuationMethodSha256'):require(r[k]==g[k],'exact assessed build '+k)
    baseline=read(m['methodRecord']);require(g['methodRecordSha256']==sha(m['methodRecord']) and baseline['jointIndependentMethodClearance'] is True and baseline['methodSha256']==sha(m['methodPath']),'original numerical method remains fixed')
    prior=read(m['priorResultRecord']);require(prior['status']=='FAILED_PRESERVED_TAIL_ALLOCATION_BEFORE_NUMERICAL_QUADRATURE' and prior['allInspectionChecksPassed'] is True and prior['records']==m['priorFiles'],'exact failed prefix index')
    history=read(m['tailMethodRecord']);require(history['jointIndependentMethodClearance'] is False and history['reports']['claude']['literalVerdict']=='NEEDS REVISION' and history['reports']['grok']['literalVerdict']=='CLEAR FOR THIS SEPARATE ANALYTIC-TAIL ALLOCATION METHOD','literal earlier method findings preserved')
    a=read(g['authority']);require(a['scope']==g['scope']==m['scope'] and a['boundedInstrumentAuthorized'] is True and a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and a['automaticScientificRetry'] is False and a['noDeadline'] is True and g['durationLimits'] is None,'standing bounded continuation authority')
    require(m['resources']['durationLimits'] is None and m['noCompletedPrefixReplay'] is True,'no deadlines or completed replay')
    return g



def load_library(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module

def run(m,J):
    raw,copies=copy_inputs(m,J)
    C=load_library(m['restoreLibrary'],'contracted_restore')
    G=load_library(m['geometryLibrary'],'contracted_geometry')
    N=load_library(m['librarySource'],'contracted_numeric')
    P=load_library(m['prepareLibrary'],'contracted_prepare')
    E=load_library(m['storeLibrary'],'contracted_store')
    I=load_library(m['requestLibrary'],'contracted_request')
    contexts,plans,tails=P.prepare(raw,m,J,sp,C,G,N)
    rules={name:raw['saved/rules/'+name+'.json'] for name in ('A-GL24','A-GL48','B-G7-K15')}
    store=E.EvidenceStore(J.out/'numerical.sqlite');index=I.RequestIndex(store,cache_bytes=0)
    try:
        mathns=store.namespace({'route':'mathematical-inputs','settings':{'scope':m['scope'],'sources':m['sourcePins'],'manifestSha256':m['selfSourceIdentity'],'noOldFunctionCalls':True}})
        mathreceipt=index._put(mathns,'complete-numerical-context',{'contexts':contexts,'ruleOperands':rules,'sourceInputs':{k:v for k,v in m['savedInputs'].items()},'tails':packed(tails),'newGeometryReceipts':contexts['geometryReceipts']})
        J.emit('numeric-storage-admission',{'freeBytes':shutil.disk_usage(J.out).free,'reserveBytes':index.reserve,'maxRecordBytes':index.maximum,'payloadCacheBytes':index.cache_limit,'ruleConstructionCalled':False,'mathematicalReceipt':mathreceipt,'activeLeavesOnDisk':True})
        require(shutil.disk_usage(J.out).free>index.reserve,'initial numerical storage reserve')
        engine=N.Evaluator(store,index,contexts,rules,G,mathreceipt);results={}
        for purpose in ('baseline','wrong-root-mutant','derivative-mutant'):
            ns=(0,2) if purpose=='baseline' else ((0,) if purpose=='wrong-root-mutant' else (2,))
            for key,plan in plans.items():
                for route in ('A24','A48','B50'):
                    for n in ns:
                        name=purpose+'/'+key+'/'+route+'/'+str(n)
                        J.start('numeric-'+name.replace('/','-'),{'purpose':purpose,'K':plan['K'],'T':plan['T'],'carrier':plan['carrier'],'route':route,'n':n,'geometryReceipt':contexts['geometryReceipts'][key]})
                        value=engine.action(route,purpose,n,plan)
                        encoded=N.encode(value);J.emit('numeric-'+name.replace('/','-')+'-full',encoded)
                        results[name]=encoded;J.finish({'complete':True,'sqliteRecords':store.sequence,'sqliteChain':store.previous})
        c=engine.mp
        def addressed(purpose,key,route):
            rows={}
            for name,record in results.items():
                if name.startswith(purpose+'/'+key+'/'+route+'/'):
                    value=N.decode(c,record)
                    require(not (set(rows)&set(value['totals'])),'unique numeric address/primitive participation')
                    rows.update(value['totals'])
            return rows
        def grouped(rows):
            values={key:N.Estimate(v['value'],v['error']) for key,v in rows.items()}
            entries={e['addressId']:e for e in contexts['entries']}
            for key,v in rows.items():
                ident,primitive=key.split('/');entry=entries[int(ident)];part=N.Estimate(v['value'],v['error'])
                for group in ('face/'+entry['face']+'/'+primitive,'grade/11/'+primitive,'component/'+('J' if primitive=='J' else 'D'),'subtotal/J-plus-D'):
                    values[group]=values.get(group,N.Estimate(c.zero,c.zero))+part
            return values
        comparisons=[]
        def compare(label,left,right,reference):
            require(set(left)==set(right)==set(reference),'full comparison vector census')
            rows=[]
            for key in left:
                diff=abs(left[key].value-right[key].value);tau=c.mpf('1e-9')+c.mpf('1e-7')*abs(reference[key].value)
                rows.append({'key':key,'left':left[key],'right':right[key],'difference':diff,'tau':tau,'passed':bool(diff<=tau),'empiricalIndicators':left[key].error+right[key].error})
            J.emit('new-comparison-'+label,N.encode(rows));require(all(r['passed'] for r in rows),'actual independent numerical comparison '+label);comparisons.extend(rows)
        for purpose in ('baseline','wrong-root-mutant','derivative-mutant'):
            for key in plans:
                vals={r:grouped(addressed(purpose,key,r)) for r in ('A24','A48','B50')}
                compare(purpose+'-'+key.replace('/','-')+'-A24-A48',vals['A24'],vals['A48'],vals['A48'])
                compare(purpose+'-'+key.replace('/','-')+'-A48-B50',vals['A48'],vals['B50'],vals['A48'])
            for carrier in ('0','1'):
                for route in ('A24','A48','B50'):
                    base=grouped(addressed(purpose,'27/'+carrier,route));large=grouped(addressed(purpose,'29/'+carrier,route))
                    compare(purpose+'-'+carrier+'-'+route+'-window-enlargement',base,large,grouped(addressed(purpose,'29/'+carrier,'A48')))
        controls=[]
        for purpose,address in [('wrong-root-mutant','8347/Dr'),('derivative-mutant','8350/J')]:
            for carrier in ('0','1'):
                values={}
                for K in (27,29):
                    key=str(K)+'/'+carrier
                    values[K]={route:{p:addressed(p,key,route)[address] for p in ('baseline',purpose)} for route in ('A24','A48','B50')}
                movement=[]
                # All contributions to the finite-window indicator are positive;
                # no inherited baseline tail is assigned to a changed mutant.
                enlargement=sum(abs(values[27][r][p]['value']-values[29][r][p]['value']) for r in values[27] for p in ('baseline',purpose))
                cross=sum(abs(values[K]['A24'][p]['value']-values[K]['A48'][p]['value'])+abs(values[K]['A48'][p]['value']-values[K]['B50'][p]['value']) for K in (27,29) for p in ('baseline',purpose))
                for K in (27,29):
                    for route in values[K]:
                        pair=values[K][route];delta=pair[purpose]['value']-pair['baseline']['value']
                        envelope=pair[purpose]['error']+pair['baseline']['error']+cross+enlargement
                        movement.append({'K':K,'route':route,'baseline':pair['baseline'],'mutant':pair[purpose],'movement':delta,'finiteEmpiricalEnvelope':envelope,'passed':bool(abs(delta)>10*envelope),'allRealMutantClaim':False})
                record={'purpose':purpose,'address':address,'carrier':carrier,'enlargement':enlargement,'crossRoute':cross,'rows':movement}
                J.emit('new-numerical-control-'+purpose+'-'+carrier,N.encode(record));require(all(r['passed'] for r in movement),'responsive actual finite-window control '+purpose);controls.append(N.encode(record))
        J.emit('new-all-address-results',results)
        J.emit('new-baseline-positive-tails',tails)
        return {'status':'BOUNDED_ORDINARY_J_AND_ADDED_DIRECT_SUBTOTALS_COMPLETE','results':results,'eligibleAddresses':20,'allSelectedAddresses':544,'addressPrimitives':40,'controls':controls,'comparisons':len(comparisons),'completePacketAction':False,'leakageFactor':None,'tails':tails,'limits':['The two full routes share an assessed analytic constant-Gaussian identity.','Adaptive and comparison indicators are empirical, not rigorous total quadrature error bounds.','H, flat, height/contact/PV and slope remain pending. J is not full native mixed.','Finite-window controls do not establish all-real mutated responses.','Inferred gamma units and original algebra are inherited. No current, inverse, scattering or loss.']}
    finally:
        J.emit('numerical-journal-final-receipt',{'records':store.sequence,'head':store.previous,'bytes':store.path.stat().st_size,'requestStates':dict(store.db.execute('SELECT state,count(*) FROM request_index GROUP BY state')),'leafState':N.leaf_count(store.db),'freeBytes':shutil.disk_usage(J.out).free,'reserveBytes':index.reserve,'noAutomaticResume':True})
        store.close()


def main():
    global sp
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);args=p.parse_args()
    m=read(args.inputs);g=verify_gate(args.gate,args.inputs,m)
    tail=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(sys.argv==tail and g['command'][-len(tail):]==tail and str(args.out)==g['outputDirectory'],'actual worker argv/output')
    tree=ast.parse(Path(m['helperSource']).read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='containment'];require(len(nodes)==1,'one inert containment helper')
    ns={'require':require,'Path':Path,'os':os,'resource':resource,'THREADS':THREADS};exec(compile(ast.Module(body=nodes,type_ignores=[]),m['helperSource'],'exec'),ns)
    enforced=ns['containment']();args.out.mkdir(exist_ok=False);J=Journal(args.out);start=time.monotonic();error=None;identities={}
    documents={str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate),g['buildReviewRecord']:g['buildReviewRecordSha256'],g['authority']:g['authoritySha256'],m['methodRecord']:g['methodRecordSha256']}
    m['selfSourceIdentity']=sha(args.inputs)
    try:
        for i,(path,digest) in enumerate(documents.items()):
            dest=args.out/('identity-'+str(i)+'-'+Path(path).name);require(sha(path)==digest,'identity source');shutil.copyfile(path,dest);require(sha(dest)==digest,'identity copy');identities[path]={'path':dest.name,'sha256':digest,'bytes':dest.stat().st_size}
        J.emit('additional-identity-copies',identities);J.emit('actual-containment',enforced)
        import sympy as sp
        import mpmath
        require(mpmath.__version__==m['mpmathVersion'],'pinned numerical runtime version')
        require(sp.__version__==m['sympyVersion'],'pinned symbolic runtime version')
        result=run(m,J)
    except BaseException:
        error=traceback.format_exc();result={'status':'FAILED_PRESERVED','failure':error,'activeOperation':J.active,'completedOperations':J.completed};J.emit('failure',result)
    finally:
        post={p:posthash(p,h) for p,h in {**m['sourcePins'],**documents}.items()};identity_post={v['path']:posthash(args.out/v['path'],v['sha256']) for v in identities.values()}
        copied={str(p.relative_to(args.out)):sha(p) for p in (args.out/'saved').rglob('*') if p.is_file()};expected={'saved/'+a:r['sha256'] for a,r in m['savedInputs'].items()}
        J.emit('posthashes',{'sources':post,'identities':identity_post,'copied':copied,'expectedCopied':expected,'copiesIntact':copied==expected})
        result.update(wallMilliseconds=round((time.monotonic()-start)*1000),sourcePosthashesIntact=all(v['intact'] for v in post.values()),copiesIntact=copied==expected,identityCopiesIntact=len(identities)==len(documents) and all(v['intact'] for v in identity_post.values()),completedOperations=J.completed,activeOperation=J.active)
        prior_post={name:{'source':posthash(Path(m['priorRoot'])/name,v['sha256']),'copy':posthash(args.out/'prior'/name,v['sha256'])} for name,v in m['priorFiles'].items()}
        prior_trees={'source':tree_inventory(m['priorRoot']),'copy':tree_inventory(args.out/'prior')}
        J.emit('prior-final-trees',prior_trees)
        tree_error=None
        try:
            for observed in prior_trees.values():check_tree(observed,m['priorFiles'])
        except ValueError as ex:tree_error=str(ex)
        J.emit('prior-copy-posthashes',{'files':prior_post,'treeRefusal':tree_error})
        result['priorCopiesIntact']=tree_error is None and all(v['source']['intact'] and v['copy']['intact'] for v in prior_post.values())
        J.emit('journal-result',result);save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text());sys.stdout.flush()
    if error:sys.stderr.write(error);return 1
    require(result['sourcePosthashesIntact'] and result['copiesIntact'] and result['identityCopiesIntact'] and result['priorCopiesIntact'],'posthash integrity');return 0


if __name__=='__main__':sys.exit(main())
