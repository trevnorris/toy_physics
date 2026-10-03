#!/usr/bin/env python3
"""Saved-evidence RIGHT grazing/control continuation; never replay the prefix."""
import argparse,ast,base64,builtins,hashlib,importlib,io,json,os,pickle,resource,shutil,sys,time,traceback
from collections import OrderedDict
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
G=((0,0),(1,0),(0,1),(1,1))
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
ZERO={'text':'0','srepr':'Integer(0)'}

def require(value,message):
    if value is not True:raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def save(path,value):
    """Stream JSON without building a second complete string; keep partial failures."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n');f.flush();os.fsync(f.fileno())

def definitions(text,names):
    nodes=[n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in nodes}==set(names),'exact inert definition census')
    return ast.Module(body=nodes,type_ignores=[])

def source_parts(text):
    run=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
    ends=[n for n in run.body if isinstance(n,ast.For) and ast.unparse(n.target)=='(end, side)']
    require(len(ends)==1,'one original end loop');end=ends[0]
    hits=[(i,n) for i,n in enumerate(end.body) if isinstance(n,ast.For) and ast.unparse(n.target)=='sign']
    require(len(hits)==1,'one original grazing loop');i,loop=hits[0]
    pairing=[j for j,n in enumerate(loop.body) if "J.join(end + '-native-%d-pairing' % sign" in ast.unparse(n)]
    require(len(pairing)==1 and pairing[0]==5,'exact failed pairing statement')
    coverage=[n for n in end.body if isinstance(n,ast.For) and ast.unparse(n.target)=='(index, c)']
    require(len(coverage)==1,'original unpublished coverage assembly')
    return {'grazingFull':loop.body,'grazingPending':loop.body[pairing[0]:],
            'controls':end.body[i+1:],'coverage':coverage}

def ast_sha(nodes):
    return hashlib.sha256(ast.dump(ast.Module(body=nodes,type_ignores=[]),include_attributes=False).encode()).hexdigest()

def execute_fragment(nodes,context,label):
    exec(compile(ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[])),label,'exec'),context)

def select_path(value,path):
    for key in path:
        require(type(key) in (str,int),'literal selector only');value=value[key]
    return value

class OperandReferences:
    """Exact typed operands already live in pinned immutable raw blobs.

    This stores their full blob receipt and literal selector, not a printed summary.
    Only the exact object selected during contained restoration may use the receipt.
    """
    def __init__(self,blobs,receipts):self.blobs,self.receipts,self.objects=blobs,receipts,{}
    def register(self,member,selector=()):
        value=select_path(self.blobs[member],selector)
        ref={**self.receipts[member],'member':member,'selector':list(selector),'representation':'PINNED_TYPED_OPERAND'}
        self.objects.setdefault(id(value),(value,ref));return value
    def reference(self,value):
        require(id(value) in self.objects,'operand has no original byte/selector provenance')
        original,ref=self.objects[id(value)]
        require(value is original and select_path(self.blobs[ref['member']],ref['selector']) is value,'exact operand object provenance')
        return dict(ref)

def make_evidence(base,exact,references):
    class Evidence(base):
        def __init__(self,out):super().__init__(out);self.join_cache={}
        def emit(self,name,value):
            if name.endswith('-actual-input'):
                return super().emit(name,{'operand':references.reference(value),'printedScientificSummary':False})
            return super().emit(name,value)
        def join(self,name,a,b):
            # Persist exact operand receipts BEFORE a comparison can fail. The
            # original unchanged predicate is used; no serialization comparison.
            previous=self.active;self.active=name
            operands={'actual':references.reference(a),'expected':references.reference(b),'predicate':'unchanged exact_structure'}
            inp=super().emit(name+'-operands',operands);key=(id(a),id(b))
            if key in self.join_cache:
                prior=self.join_cache[key];passed=prior['passed'];reuse=prior['return']
            else:passed=exact(a,b);reuse=None
            out=super().emit(name,{'operands':inp,'passed':passed,'reusedReturn':reuse,'predicateCalled':reuse is None})
            require(passed,name);self.join_cache[key]={'passed':passed,'return':out};self.active=previous
    return Evidence

def verify_gate(path,manifest_path,m):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_SAVED_END_UNIFORM_CONTINUATION','one continuation gate')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest identity')
    require(g['sourcePins']==m['sourcePins'],'source pin census')
    for p,h in m['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'actual containment route')
    for key in ('sharedGuard','supervisor','launcher','baseGate','baseManifest','baseBuildRecord','authority','repairRecord','priorFiles','completionRecord'):
        require(sha(g[key])==g[key+'Sha256']==m['sourcePins'][g[key]],'pinned route '+key)
    old=json.loads(Path(g['baseGate']).read_text());build=json.loads(Path(g['baseBuildRecord']).read_text())
    require(old['workerSha256']==sha(m['baseWorker']) and old['manifestSha256']==sha(g['baseManifest']),'preserved actual prior worker')
    require(old['buildReviewRecordSha256']==sha(g['baseBuildRecord']) and build['independentBuildClearance'] is False,'literal historical build')
    require(g['continuationIndependentBuildClearance'] is False and g['localToolingRepairAuthorized'] is True,'honest tooling authority')
    for key in ('sharedGuard','supervisor'):require(old[key]==g[key] and old[key+'Sha256']==g[key+'Sha256'],'unchanged resource helper')
    authority=json.loads(Path(g['authority']).read_text())
    require(authority['savedEvidenceContinuationAuthorized'] is True and authority['automaticScientificRetry'] is False and authority['scientificRunsAuthorized']==1 and authority['durationLimits'] is None and authority['scope']==m['scope'],'one bounded saved continuation')
    require(g['scope']==m['scope'] and g['scientificRunsAuthorized']==1 and g['durationLimits'] is None and g['launcher']==m['launcher'],'exact gate scope')
    require(m['resources']==json.loads(Path(g['baseManifest']).read_text())['resources'],'unchanged resource envelope')
    return g

def restore_prior(m,out):
    prior=Path(m['priorDirectory']);target=out/'prior';target.mkdir();tree=json.loads(Path(m['priorFiles']).read_text());copies={}
    require(len(tree)==6885 and sum(v['bytes'] for v in tree.values())==540133852,'exact prior file census')
    for name,r in tree.items():
        src=prior/name;dst=target/name
        require(sha(src)==r['sha256'] and src.stat().st_size==r['bytes'],'prior bytes '+name)
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        require(sha(dst)==r['sha256'],'byte-identical prior copy '+name)
        copies[name]={**r,'source':str(src),'path':str(dst.relative_to(out))}
    save(out/'prior-copy-index.json',copies)
    # Read only receipt metadata here. Every published artifact is joined to its
    # original hash-chain entry, including all completed scalar arguments/returns.
    chain=[];previous=None
    for line in (target/'evidence-chain.jsonl').open():
        rec=json.loads(line);pay=rec['payload'];ref=pay['artifact'];name=ref['path']
        require(pay['sequence']==len(chain) and pay['previousSha256']==previous and rec['sha256']==hashlib.sha256(json.dumps(pay,sort_keys=True).encode()).hexdigest(),'original evidence chain')
        require(tree[name]=={'sha256':ref['sha256'],'bytes':ref['bytes']},'original artifact receipt')
        chain.append({**ref,'status':'RESTORED_PUBLISHED_EVIDENCE','functionCalled':False});previous=rec['sha256']
    require(len(chain)==6796,'all completed evidence artifacts')
    save(out/'restored-evidence-index.json',{'records':chain,'originalLastSha256':previous,'recomputed':False})
    failure=json.loads((target/'failure.json').read_text())
    require(failure['incompleteOperation']=='end-uniform-comparison' and failure['traceback'].rstrip().endswith('MemoryError') and 'native_args[\'pairing\'],pairing' in failure['traceback'],'exact preserved storage failure')
    require(not (target/'RIGHT-native--1-pairing.json').exists() and not (target/'RIGHT-comparison.json').exists(),'failed join was never published')
    return target,copies,tree

def run_continuation(m,J,ns,base,blobs,references,prior,tree):
    D=ns['decode'];used=set()
    def raw(name):used.add(name);return json.loads((prior/name).read_text())
    def load(name):return D(raw(name))
    def receipt(name):return {'path':'prior/'+name,**tree[name]}
    bm=json.loads(Path(m['baseManifest']).read_text());operations={r['name']:r for r in bm['uniformReceipts']['completedOperations']}
    require(len(operations)==18,'all old complete operations available')
    def operation(name):
        rec=operations[name];return blobs[rec['input']['member']],blobs[rec['result']['member']]
    J.emit('restored-old-operations',{'records':[{**r,'status':'RESTORED_PRIOR_COMPLETE_RETURN','functionCalled':False} for r in operations.values()]})
    restored={name:operation('restore/'+name)[1] for name in bm['originalPackets']}
    # Pinned copies already supply the input joins from the successful prefix.
    for name,pin in bm['originalPackets'].items():
        inp,_=operation('restore/'+name)
        require(inp['pin']==pin and inp['originalBytes']['sha256']==pin['sha256'] and inp['originalBytes']['bytes']==pin['bytes'],'original saved restore arguments '+name)
    J.emit('restored-prefix-disposition',{'completedEvidence':6796,'literalZeroReturns':2281,'sourceChartUnitsNormalizationReplayed':False,'selectedComparisonsReplayed':False,'leftReplayed':False,'rightPending':'native--1-pairing','originalFailure':receipt('failure.json')})
    context=load('saved/ends/local-context.json');cells=load('saved/ends/symbols.json');depth=load('saved/ends/depth-domain.json')
    p,q,cs=depth['p'],depth['q'],depth['cs'];eta,sigma=context['context']['independentGrades']
    S=load('actual-source-chart.json')['matrixWeakFromUniform'];side_evidence=raw('actual-end-side-map.json')
    side_map={r['end']:r['weakSide'] for r in side_evidence};require(side_map=={'LEFT':'minus','RIGHT':'plus'},'actual saved side map')
    end='RIGHT';side=side_map[end];inp,R=operation(end+'/restriction');freq=restored[end+'Frequency'];pairing,dims=restored[end+'Pairing']
    atoms_R=set().union(*(x.free_symbols for x in (R['P'],R['D'],R['wave'],R['lift'],R['conversion'])))
    def symbol(name):
        found=[s for s in atoms_R if s.name==name];require(len(found)==1,'saved symbol '+name);return found[0]
    w,k,cold,qold=(symbol(n) for n in ('uniformFrequency','uniformNormal','uniformSoundSpeed','uniformPhysicalDepth'))
    common_map={w:sp.Integer(3),k:p,cold:cs,qold:q};old=lambda v:v.xreplace(common_map)
    Dold,lift=map(old,(R['D'],R['lift']))
    require(base['exact_structure'](R['D'],load('RIGHT-invariant-D.json')['actual']) and base['exact_structure'](R['lift'],load('RIGHT-invariant-lift.json')['actual']),'actual saved raw restriction context')
    # Only the unsaved small matrix/argument context is reconstituted from the
    # published cell values. No cell, endpoint, source or grade function runs.
    groups={g:sp.zeros(5) for g in G};closed={g:sp.zeros(5) for g in G};seen=set()
    for c in cells:
        if c['side']!=side:continue
        key=(ROWS.index(c['row']),FIELDS.index(c['field']),tuple(c['grade']));i,j,g=key
        require(key not in seen and c['sumIdentity']['cancelled']==0,'saved right cell uniqueness and proof');seen.add(key)
        groups[g][i,j]=c['symbol'];closed[g][i,j]=c['closedGrazingValue']
    require(len(cells)==200 and len(seen)==100,'all saved end cells')
    Egrades={g:S.T*groups[g]*S for g in G};Ezero={g:S.T*closed[g]*S for g in G}
    origin={s.name:v for s,v in freq['origin'].items()}
    require(origin=={'eta_bg':sp.Rational(1,100),'sigma_W':sp.Rational(1,1000)},'actual physical origin')
    finite_origin={eta:origin['eta_bg'],sigma:origin['sigma_W']}
    Ephys=sum((eta**a*sigma**b*Egrades[a,b] for a,b in G),sp.zeros(5)).subs(finite_origin)
    Eg0=sum((eta**a*sigma**b*Ezero[a,b] for a,b in G),sp.zeros(5)).subs(finite_origin)
    Rnew=sp.MutableDenseMatrix(5,2,lambda i,j:load('RIGHT/selected-R-%d-%d-input.json'%(i,j))['value'])
    selected=[[load('RIGHT/selected-A-%d-%d-return.json'%(i,j)) for j in range(2)] for i in range(5)]
    residual=[[load('RIGHT/selected-R-%d-%d-return.json'%(i,j)) for j in range(2)] for i in range(5)]
    full=[[load('RIGHT/full-Delta-%d-%d-return.json'%(i,j)) for j in range(5)] for i in range(5)]
    grade=[]
    for i in range(5):
        for j in range(2):
            nm='RIGHT/grade-%d-%d'%(i,j);a=load(nm+'-old-return.json');b=load(nm+'-new-return.json')
            jt={g:load(nm+'-J-%d%d-return.json'%g) for g in G};at=load(nm+'-attribution-return.json')
            require(a['status']==b['status']=='EXTRACTED_ON_DECLARED_REGULAR_DOMAIN' and b['remainder']==0 and at['status']=='ZERO','saved successful grade attribution')
            require(all(v['status']=='ZERO' for v in jt.values()) and selected[i][j]['status']=='ZERO','saved selected grade statuses')
            grade.append({'row':i,'column':j,'retained':'AGREEMENT','finite':selected[i][j]['status'],'classification':'AGREEMENT','J':jt,'H':a['remainder'],'regularity':a['originDenominator']})
    require(all(v['status']=='ZERO' for rows in (selected,residual,full) for row in rows for v in row),'actual saved right conditional scalar results')
    J.emit('restored-right-partial',{'fullDifference':full,'selectedFinite':selected,'fiveRowResidual':residual,'grades':grade,'recomputed':False})
    J.emit('reconstituted-small-context',{'symbols':(p,q,cs,eta,sigma,w,k,cold,qold),'map':common_map,'chart':S,'physicalOrigin':finite_origin,'Egrades':Egrades,'closedGrades':Ezero,'Ephys':Ephys,'closedE':Eg0,'Dold':Dold,'lift':lift,'savedRnew':Rnew,'cellSource':receipt('saved/ends/symbols.json'),'rawRestriction':operations['RIGHT/restriction'],'sourceFunctionsCalled':False})
    # Every large operand is a selector into one of the exact copied raw blobs.
    for name,selector in [('restore/RIGHTPairing',(0,)),('restore/RIGHTNativeSource',()),('RIGHT/restriction',()),('RIGHT/native-reconstruction',())]:
        references.register(operations[name]['result']['member'],selector)
    native_member=operations['RIGHT/native-reconstruction']['input']['member']
    for selector in (('pairing',),('native',)):references.register(native_member,selector)
    for sign in (-1,1):
        member=operations['RIGHT/exact-limits/'+str(sign)]['input']['member']
        for selector in ((),('restriction',),('native',)):references.register(member,selector)
    source=Path(m['baseWorker']).read_text();parts=source_parts(source)
    require({k:ast_sha(v) for k,v in parts.items()}==m['originalFragmentHashes'],'unchanged original tail ASTs')
    left=load('LEFT-comparison.json')
    require(len(left['grazing'])==2 and len(left['controls'])==3 and all(v['responsive'] is True for v in left['controls']),'completed LEFT grazing and controls restored')
    results={'LEFT':{'retainedStatuses':sorted({v['retained'] for v in left['grades']}),'finiteStatuses':sorted({v['status'] for row in left['selectedFinite'] for v in row}),'controlCoverage':{c['name']:c['responsive'] for c in left['controls']}}}
    ctx={'sp':sp,'J':J,'require':require,'finite_constant':base['finite_constant'],'operation':operation,'blobs':blobs,'restored':restored,
         'end':end,'side':side,'R':R,'pairing':pairing,'old':old,'w':w,'k':k,'qold':qold,'p':p,'q':q,'cs':cs,'eta':eta,'sigma':sigma,
         'Eg0':Eg0,'Egrades':Egrades,'Ephys':Ephys,'Dold':Dold,'lift':lift,'Rnew':Rnew,'finite_origin':finite_origin,'G':G,'ROWS':ROWS,'FIELDS':FIELDS,
         'cells':cells,'roworder':(1,2,0,3,4),'S':S,'coverage':[],'grazing':[],'full':full,'selected':selected,'residual':residual,'grade':grade,'results':results}
    # The first five statements of the minus target were published. Restore
    # those receipts and actual old operands, then enter at the failed join.
    args,limits=operation('RIGHT/exact-limits/-1');native_args,native_return=operation('RIGHT/native-reconstruction')
    prefix=[]
    for name in ('RIGHT-limit--1-actual-input.json','RIGHT-limit--1-restriction.json','RIGHT-limit--1-native.json'):
        if name.endswith(('restriction.json','native.json')):require(raw(name)['passed'] is True,'saved completed minus-prefix join')
        prefix.append(receipt(name))
    # Metadata byte joins above link the complete old operation to precisely the
    # pending input already published by the prior worker; no native function.
    J.emit('restored-minus-prefix',{'published':prefix,'input':operations['RIGHT/exact-limits/-1']['input'],'nativeOperation':operations['RIGHT/native-reconstruction'],'completedFunctionsCalled':False})
    ctx.update(sign=-1,args=args,limits=limits,native_args=native_args,native_return=native_return)
    execute_fragment(parts['grazingPending'],ctx,'original-worker#unfinished-RIGHT-minus-grazing')
    ctx['sign']=1;execute_fragment(parts['grazingFull'],ctx,'original-worker#unfinished-RIGHT-plus-grazing')
    execute_fragment(parts['controls'],ctx,'original-worker#unfinished-RIGHT-controls')
    # Participation metadata was computed but not published by the failed run.
    # Reconstitute only this small assembly from the same saved cells and lifts.
    for end,side in side_map.items():
        _,rr=operation(end+'/restriction');atoms=set().union(rr['lift'].free_symbols)
        lk=next(s for s in atoms if s.name=='uniformNormal')
        ctx.update(end=end,side=side,lift=rr['lift'].xreplace({lk:p}))
        execute_fragment(parts['coverage'],ctx,'original-worker#unpublished-participation-context')
    J.emit('selected-cell-coverage',ctx['coverage']);require(len(ctx['coverage'])==200,'complete actual saved-cell participation')
    J.emit('consumed-prior-receipts',{n:receipt(n) for n in sorted(used)})
    return {'executionStatus':'COMPLETE_COMPARISON_PENDING_INSPECTION','ends':results,'cells':200,'restoredOperations':18,'restoredLeftComparison':receipt('LEFT-comparison.json'),
            'completedPrefixReplayed':False,'newRoots':False,'fieldOrLoss':False,'frequencyDerivativeTransfer':False,'oldFunctionsCalled':False,
            'originalFailurePreserved':receipt('failure.json'),'continuationIndependentBuildClearance':False}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--gate',type=Path,required=True);args=parser.parse_args()
    m=json.loads(args.inputs.read_text());g=verify_gate(args.gate,args.inputs,m)
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(sys.argv==expected and g['command'][-len(expected):]==expected and args.out.resolve()==Path(g['outputDirectory']).resolve(),'actual worker invocation')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;code=1;result={};start=time.monotonic();copies={}
    pins={**m['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    try:
        bm=json.loads(Path(m['baseManifest']).read_text())
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(bm['helperSource']).read_text(),('require','containment','decode')),'inert-contained-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp,np
        import sympy as sp
        import numpy as np
        from sympy.core.symbol import Str
        from sympy.core.function import FunctionClass
        ns.update(sp=sp,Str=Str)
        codec_ns={'pickle':pickle,'OrderedDict':OrderedDict,'builtins':builtins,'np':np,'sp':sp,'importlib':importlib,'io':io,'require':require}
        exec(compile(definitions(Path(bm['uniformWorker']).read_text(),('SavedCodec','decode')),'unchanged-saved-codec','exec'),codec_ns)
        base={'sp':sp,'np':np,'FunctionClass':FunctionClass,'Path':Path,'base64':base64,'hashlib':hashlib,'json':json,'os':os,'save':save,'sha':sha,'require':require}
        exec(compile(definitions(Path(m['baseWorker']).read_text(),('Evidence','exact_structure','finite_constant')),'unchanged-scientific-predicate-and-tail-helper','exec'),base)
        prior,copies,tree=restore_prior(m,args.out)
        refs=json.loads((prior/'opaque-copy-index.json').read_text());blobs={};receipts={}
        for row in bm['uniformReceipts']['blobReceipts']:
            member=row['member'];ref=refs[member];path=prior/ref['path']
            require(sha(path)==row['sha256']==ref['sha256'] and path.stat().st_size==row['bytes']==ref['bytes'],'actual saved raw operand '+member)
            blobs[member]=codec_ns['decode'](path.read_bytes());receipts[member]={**ref,'path':'prior/'+ref['path']}
        require(len(blobs)==58,'all actual raw blob copies');references=OperandReferences(blobs,receipts)
        Evidence=make_evidence(base['Evidence'],base['exact_structure'],references);J=Evidence(args.out)
        result=J.stage('saved-end-uniform-right-continuation',{'manifest':g['manifestSha256'],'priorCompletion':g['completionRecordSha256'],'baseWorker':sha(m['baseWorker'])},lambda:run_continuation(m,J,ns,base,blobs,references,prior,tree));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        for r in copies.values():pins[r['source']]=r['sha256'];pins[str(args.out/r['path'])]=r['sha256']
        posts={}
        for p,h in pins.items():
            try:posts[p]={'expected':h,'actual':sha(p),'error':None}
            except OSError as e:posts[p]={'expected':h,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',posts)
        if any(v['expected']!=v['actual'] for v in posts.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        if J is not None:result=J.encode(result)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
        if J is not None:save(args.out/'evidence-final-receipt.json',{'records':J.count,'lastSha256':J.previous,'chainSha256':sha(args.out/'evidence-chain.jsonl'),'completeProcess':code==0})
    return code

if __name__=='__main__':sys.exit(main())
