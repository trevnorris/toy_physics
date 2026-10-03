#!/usr/bin/env python3
"""Saved weak-end controls only; no completed prefix replay."""
import argparse
import ast
import hashlib
import itertools
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback
ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','posthash_records','containment','Journal','decode','one_symbol','expanded_sinh_arguments')
G=((0,0),(1,0),(0,1),(1,1))
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
ZERO={'text':'0','srepr':'Integer(0)'}


def require(v, message):
    if v is not True:
        raise ValueError(message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def save(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())

def definitions(text, names):
    nodes = [n for n in ast.parse(text).body
             if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in nodes} == set(names), 'exact helper census')
    return ast.Module(body=nodes, type_ignores=[])

def verify_helper_paths(gate):
    require(gate['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py'), 'actual guard route')
    require(gate['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),
            'actual supervisor route')

def verify_invocation(args, gate, argv):
    require(args.out.resolve()==Path(gate['outputDirectory']).resolve(),'gate output route')
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),
              '--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(list(argv)==expected,'actual worker argv')
    require(gate['command'][-len(expected):]==expected,'gate worker command tail')

class EvidenceLog:
    """Append-only operands/returns without rewriting an ever-growing index."""
    def __init__(self,path,encode):
        self.path=path;self.encode=encode;self.count=0;self.previous=None
        with path.open('x') as f:f.flush();os.fsync(f.fileno())
    def append(self,kind,value):
        payload={'sequence':self.count,'kind':kind,'previousSha256':self.previous,'value':self.encode(value)}
        digest=hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False).encode()).hexdigest()
        with self.path.open('a') as f:
            f.write(json.dumps({'payload':payload,'sha256':digest},allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
        self.count+=1;self.previous=digest
        return {'sequence':self.count-1,'sha256':digest}

def known_true(value,true_atom):
    """Accept only explicit true atoms; never coerce unknowns, numbers or strings."""
    return value is True or value is true_atom


def source_parts(source):
    run=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
    hits=[i for i,n in enumerate(run.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
          and isinstance(n.value.func,ast.Name) and n.value.func.id=='require' and len(n.value.args)==2
          and isinstance(n.value.args[1],ast.Constant) and n.value.args[1].value=='control speed in scoped interval']
    require(len(hits)==1,'unique stopped interval predicate')
    setup=[n for n in run.body[:hits[0]] if isinstance(n,ast.Assign)
           and any(isinstance(x,ast.Name) and x.id in ('controls','point','cs2') for x in n.targets)]
    require(len(setup)==3,'original unpublished point setup')
    nested=[n for n in run.body if isinstance(n,ast.FunctionDef) and n.name in ('zero','constant')]
    require({n.name for n in nested}=={'zero','constant'},'unchanged new-control scalar helpers')
    return setup,run.body[hits[0]+1:],nested


def call_tail(source,context):
    _,tail,_=source_parts(source)
    fn=ast.FunctionDef(name='unfinished_controls',args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[]),body=tail,decorator_list=[])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'original-worker#unmodified-control-tail','exec'),context)
    return context['unfinished_controls']()


def verify_gate(path,manifest_path,manifest):
    g=json.loads(Path(path).read_text());verify_helper_paths(g)
    require(g['status']=='READY_FOR_ONE_SAVED_WEAK_END_CONTROL_CONTINUATION','continuation gate')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'current worker/manifest')
    require(g['sourcePins']==manifest['sourcePins'],'full pin census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','baseGate','baseManifest','baseBuildRecord','authority','repairRecord'):
        require(sha(g[key])==g[key+'Sha256']==g['sourcePins'][g[key]],'pinned helper/authority '+key)
    require(g['launcher']==manifest['launcher'],'actual launcher path')
    old=json.loads(Path(g['baseGate']).read_text());build=json.loads(Path(g['baseBuildRecord']).read_text())
    require(build['independentBuildClearance'] is True and all(r['literalVerdict']=='CLEAR FOR THIS TRANSLATED WEAK-END BUILD' for r in build['reports'].values()),'literal original paired build')
    require(sha(manifest['baseWorker'])==build['workerSha256']==old['workerSha256'],'reviewed original worker')
    require(sha(g['baseManifest'])==build['manifestSha256']==old['manifestSha256'],'reviewed original manifest')
    for key in ('sharedGuard','supervisor'):
        require(g[key]==old[key] and g[key+'Sha256']==old[key+'Sha256']==build[key+'Sha256'],'same reviewed runtime '+key)
    require(g['continuationIndependentBuildClearance'] is False and g['localToolingRepairAuthorized'] is True,'honest local-tooling authority')
    authority=json.loads(Path(g['authority']).read_text())
    require(authority['savedEvidenceContinuationAuthorized'] is True and authority['automaticScientificRetry'] is False
            and authority['scientificRunsAuthorized']==1 and authority['durationLimits'] is None and authority['scope']==manifest['scope'],'one saved-control continuation')
    require(g['scientificRunsAuthorized']==1 and g['scope']==manifest['scope'] and g['durationLimits'] is None and g['priorTree']==manifest['priorTree'],'exact no-deadline bounded scope')
    return g


def run_continuation(manifest,J,helpers):
    prior=Path(manifest['priorDirectory']);copied=J.out/'prior';copied.mkdir();copies={};used=set()
    for n,r in manifest['priorTree'].items():
        src=prior/n;dst=copied/n;require(sha(src)==r['sha256'] and src.stat().st_size==r['bytes'],'prior input '+n)
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);require(sha(dst)==r['sha256'],'prior copy '+n)
        copies[n]={'source':str(src),'path':str(dst.relative_to(J.out)),'sha256':r['sha256'],'bytes':r['bytes']}
    save(J.out/'prior-copy-index.json',copies)
    def raw(n):used.add(n);return json.loads((copied/n).read_text())
    D=helpers['decode'];load=lambda n:D(raw(n));index=raw('artifact-index.json')
    for n,r in index.items():require(r['path']==n and r['sha256']==sha(copied/n) and r['bytes']==(copied/n).stat().st_size,'published artifact '+n)
    failure=raw('failure.json');require(failure['incompleteOperation']=='translated-weak-ends' and 'control speed in scoped interval' in failure['traceback'],'exact prior stop')
    require('translated-weak-ends-return.json' not in index and not any(n.startswith('new-control-') for n in index),'unfinished controls only')
    restored=[];endpoints={};ops=raw('operation-index.json')
    for op in ops:
        require(op['name'].startswith('new-field-endpoints-'),'only prior completed endpoint stage')
        for key in ('input','result'):require(op[key]==index[op[key]['path']],'original completed receipt')
        inp=raw(op['input']['path']);ret=raw(op['result']['path']);require(ret['inherited']==inp,'complete endpoint input/return join')
        decoded_inputs=D(inp);decoded_return=D(ret);fid=inp['fieldId'];endpoints[fid]=decoded_return
        require(decoded_return['inherited']==decoded_inputs,'actual endpoint decode arguments')
        restored.append({**op,'status':'RESTORED_COMPLETE_RETURN','functionCalled':False})
    require(len(restored)==len(endpoints)==34,'all34 complete endpoint returns restored')
    J.emit('restored-complete-operation-index',restored)
    # Published exact prefix returns are restored with the original operands, not recomputed.
    zero_records=[];zero_index=[];byname={};pending=None;previous=None;count=0
    for line in (copied/'exact-evidence.jsonl').open():
        record=json.loads(line);pay=record['payload'];v=pay['value']
        require(pay['sequence']==count and pay['previousSha256']==previous and
                record['sha256']==hashlib.sha256(json.dumps(pay,sort_keys=True,allow_nan=False).encode()).hexdigest(),'prior append-only chain')
        if pay['kind']=='zero-input':require(pending is None,'paired zero input');pending=v
        elif pay['kind']=='zero-return':
            require(pending=={k:v[k] for k in ('name','left','right')} and v['cancelled']==ZERO,'published zero with actual operands')
            decoded=D(v);require(decoded['cancelled']==0,'decoded published zero return')
            zero_records.append(decoded);byname.setdefault(v['name'],[]).append(decoded)
            zero_index.append({'name':v['name'],'sequence':count,'chainSha256':record['sha256'],'functionCalled':False});pending=None
        previous=record['sha256'];count+=1
    require(pending is None and len(zero_records)==333 and count==27796,'all333 completed exact returns and chain')
    J.emit('restored-exact-prefix-index',{'records':zero_index,'originalChainRecords':count,'originalLastSha256':previous,'newlyCalculated':False})
    dispersion=byname['control-outgoing-physical-dispersion'];require(len(dispersion)==1 and dispersion[0]['left']==4 and dispersion[0]['right']==4,'saved physical control-point identity')
    J.emit('restored-control-point-dispersion',{'savedReturn':dispersion[0],'functionCalled':False})
    domain=load('new-closed-depth-domain.json');p=domain['p'];q=domain['q'];cs=domain['cs']
    translation=load('new-translated-response-argument.json');reversed_record=load('new-reversed-native-phase-limits.json')
    H=translation['heightLimits'];reversed_limits=reversed_record['reversedHeightLimits']
    Bh=byname['saved-height-Holder-diagonal-plus'][0]['right']
    physical_context=load('saved/full/extended-binding-context.json');W=physical_context['numeric']['W_0']
    require(W==1 and physical_context['numeric']['omega']==3 and domain['range']==[1,2],'actual saved physical context')
    require(Bh.free_symbols=={q} and all(v.free_symbols==set() for v in H.values()),'restored response/height symbol assumptions')
    symbols=load('new-complete-retained-weak-end-symbols.json')
    require(len(symbols)==200 and len({(s['side'],s['row'],s['field'],tuple(s['grade'])) for s in symbols})==200,'all saved end cells')
    require(all(s['sumIdentity']['cancelled']==0 and s['symbol'].free_symbols<={p,q} for s in symbols),'saved cell identities and actual symbols')
    source=Path(manifest['baseWorker']).read_text();tree=ast.parse(source)
    selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='control_eligible'];require(len(selected)==1,'original pure metadata selector')
    metadata_namespace={'require':require};exec(compile(ast.Module(body=selected,type_ignores=[]),'unchanged-control-metadata','exec'),metadata_namespace)
    possible_controls={k:[] for k in ('omit-height-contact','reverse-translation-phase','omit-lower-normal-sign')}
    join_records=[];seen=set();candidate_records=[]
    for row in ROWS:
        rows=raw(row+'-new-end-addresses.json')
        for r in rows:
            ar=r['inheritedAddress'];i=ar['addressId'];require(i not in seen and i==r['join']['addressId'] and ar['row']==row,'saved actual address route');seen.add(i);join_records.append(r['join'])
            kinds=[k for k in possible_controls if r['endContributions']['plus']!=ZERO and metadata_namespace['control_eligible'](ar,k)]
            if kinds:
                decoded=D(r);a=decoded['inheritedAddress'];s=endpoints[a['sourceTransform']['coefficientId']];c=endpoints[a['consumerTransform']['coefficientId']]
                require(decoded['waveAtP'].free_symbols<={p} and decoded['responseLimit']['normal'].free_symbols<={q},'saved control momentum/depth arguments')
                item=(a,s,c,decoded['waveAtP'],decoded['responseLimit'],decoded['endContributions'])
                for kind in kinds:possible_controls[kind].append(item)
                candidate_records.append({'addressId':i,'rowArtifact':row+'-new-end-addresses.json','kinds':kinds,'sourceField':a['sourceTransform']['coefficientId'],'consumerField':a['consumerTransform']['coefficientId'],'completedCalculationsCalled':False})
    require(len(seen)==13260 and len(join_records)==13260,'all saved address returns reused')
    require({k:len(v) for k,v in possible_controls.items()}=={'omit-height-contact':48,'reverse-translation-phase':48,'omit-lower-normal-sign':24},'actual saved surviving candidates')
    J.emit('restored-control-candidate-routes',candidate_records)
    audit=EvidenceLog(J.out/'exact-evidence.jsonl',J.encode)
    context={'sp':sp,'J':J,'require':require,'audit':audit,'p':p,'q':q,'cs':cs,'W':W,'Bh':Bh,'H':H,'reversed_limits':reversed_limits,
             'symbols':symbols,'possible_controls':possible_controls,'join_records':join_records,'seen':seen,'manifest':manifest,'used':used,'copies':copies}
    setup,tail,nested=source_parts(source)
    exec(compile(ast.Module(body=nested+setup,type_ignores=[]),'original-new-control-helpers-and-small-point-setup','exec'),context)
    cs2=context['cs2'];lower=cs2>1;upper=cs2<4
    J.emit('interval-predicate-operands',{'csSquared':cs2,'lowerEndpoint':1,'upperEndpoint':4,'lowerComparison':lower,'upperComparison':upper,
        'lowerType':type(lower).__module__+'.'+type(lower).__qualname__,'upperType':type(upper).__module__+'.'+type(upper).__qualname__,
        'oldFinalComparisonIsPythonTrue':upper is True,'exactTrueAtom':sp.S.true,
        'lowerAccepted':known_true(lower,sp.S.true),'upperAccepted':known_true(upper,sp.S.true),
        'savedDispersion':dispersion[0],'originalPredicateSource':"require(cs2>1 and cs2<4,'control speed in scoped interval')"})
    require(known_true(lower,sp.S.true) and known_true(upper,sp.S.true),'control speed in scoped interval')
    J.emit('restored-control-context',{'p':p,'q':q,'cs':cs,'W':W,'Bh':Bh,'heights':H,'reversedLimits':reversed_limits,
        'point':list(context['point'].items()),'csSquared':cs2,'savedSymbolCells':len(symbols),'savedAddressReturns':len(seen),
        'originalTailAstSha256':hashlib.sha256(ast.dump(ast.Module(body=tail,type_ignores=[])).encode()).hexdigest(),
        'smallUnpublishedSetup':[ast.unparse(n) for n in setup],'noCompletedFunctionCalled':True})
    result=call_tail(source,context)
    result.update(savedEvidenceContinuation=True,priorCompleteOperationsRestored=34,priorLiteralZeroReturnsRestored=333,priorSavedCellsReused=200,completedPrefixReplayed=False)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest);verify_invocation(args,gate,sys.argv)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;result={};code=1;started=time.monotonic()
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(manifest['helperSource']).read_text(),HELPERS),'unchanged-inert-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('saved-weak-end-controls',{'manifestSha256':gate['manifestSha256'],'baseWorkerSha256':sha(manifest['baseWorker']),'priorDirectory':manifest['priorDirectory']},lambda:run_continuation(manifest,J,ns));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        if (args.out/'exact-evidence.jsonl').exists():
            evidence=args.out/'exact-evidence.jsonl';save(args.out/'exact-evidence-final-receipt.json',{'path':evidence.name,'sha256':sha(evidence),'bytes':evidence.stat().st_size,'mayBeIncomplete':code!=0})
        if (args.out/'prior-copy-index.json').exists():
            for v in json.loads((args.out/'prior-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(v['expected']!=v['actual'] for v in records.values()):result['integrityFailure']=True;result['executionStatus']='INTEGRITY_FAILURE_PRESERVED';code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False);save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':sys.exit(main())
