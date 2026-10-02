#!/usr/bin/env python3
"""Resume only unpublished reference-grazing checks; exact constant-zero tooling."""
import argparse
import ast
from fractions import Fraction
import hashlib
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
BASE_DEFINITIONS=('nonnegative_polynomial','exact_nonzero_number','trace_subtraction_coefficient')
TAIL_MARKER='new-tail-polynomial'


def require(value,message):
    if value is not True:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())


def exact_rational_constant(constructor):
    """Restricted inert constructor syntax, never eval/sympify and never floats/symbols."""
    tree=ast.parse(constructor,mode='eval');steps=[]
    def integer(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return n.value
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub) and isinstance(n.operand,ast.Constant) and type(n.operand.value) is int:return -n.operand.value
        raise ValueError('literal integer required')
    def visit(n):
        require(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and not n.keywords,'restricted constant constructor')
        name=n.func.id
        if name=='Integer':
            require(len(n.args)==1,'integer arity');v=Fraction(integer(n.args[0]))
        elif name=='Rational':
            require(len(n.args)==2,'rational arity');v=Fraction(integer(n.args[0]),integer(n.args[1]))
        elif name in ('Add','Mul'):
            require(len(n.args)>0,'nonempty arithmetic constructor');v=Fraction(0 if name=='Add' else 1)
            for arg in n.args:
                x=visit(arg);v=v+x if name=='Add' else v*x
        elif name=='Pow':
            require(len(n.args)==2,'power arity');base=visit(n.args[0]);power=visit(n.args[1])
            require(power.denominator==1 and -16<=power.numerator<=16,'bounded exact integer power');v=base**power.numerator
        else:raise ValueError('unsupported constant constructor '+name)
        steps.append({'constructor':ast.unparse(n),'numerator':v.numerator,'denominator':v.denominator})
        return v
    value=visit(tree.body)
    return {'constructor':constructor,'steps':steps,'numerator':value.numerator,'denominator':value.denominator,'exactZero':value==0}


def tail_ast(source):
    tree=ast.parse(source);function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
    hits=[i for i,n in enumerate(function.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='nonnegative_polynomial' and len(n.value.args)>1 and isinstance(n.value.args[1],ast.Constant) and n.value.args[1].value==TAIL_MARKER]
    require(len(hits)==1,'unique unfinished tail boundary')
    return function.body[hits[0]:]


def tail_callable(source,namespace):
    function=ast.FunctionDef(name='unfinished_reference_tail',args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[]),body=tail_ast(source),decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[]))
    exec(compile(module,'reviewed-reference-worker#unchanged-unfinished-tail','exec'),namespace)
    return namespace['unfinished_reference_tail']


def verify_gate(path,manifest_path,manifest):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_SAVED_REFERENCE_CONTINUATION','continuation gate')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'current worker/manifest')
    require(g['sourcePins']==manifest['sourcePins'],'complete source pin set')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for k in ('sharedGuard','supervisor','launcher','baseGate','baseManifest','baseBuildRecord','authority'):
        require(sha(g[k])==g[k+'Sha256'],'actual helper/evidence '+k)
    for k in ('sharedGuard','supervisor','launcher'):require(g[k+'Sha256']==g['sourcePins'][g[k]],'executed pinned helper '+k)
    build=json.loads(Path(g['baseBuildRecord']).read_text());old=json.loads(Path(g['baseGate']).read_text())
    require(build['independentBuildClearance'] is True and all(x['literalVerdict']=='CLEAR FOR THIS BOUNDED REFERENCE-GRAZING BUILD' for x in build['reports'].values()),'literal original build')
    require(build['workerSha256']==sha(manifest['baseWorker'])==old['workerSha256'],'original reviewed worker')
    require(build['manifestSha256']==sha(g['baseManifest'])==old['manifestSha256'],'original reviewed manifest')
    require(g['continuationIndependentBuildClearance'] is False and g['localToolingRepairAuthorized'] is True,'honest continuation authority')
    a=json.loads(Path(g['authority']).read_text())
    require(a['savedEvidenceContinuationAuthorized'] is True and a['automaticScientificRetry'] is False,'standing saved-evidence scope')
    require(g['scientificRunsAuthorized']==1 and g['durationLimits'] is None and g['scope']==manifest['scope'],'one bounded no-deadline continuation')
    require(g['priorTree']==manifest['priorTree'],'actual prior result census')
    return g


def run_continuation(manifest,J,helpers,base):
    prior=Path(manifest['priorDirectory']);copied=J.out/'prior';copied.mkdir();copy_index={}
    for name,r in manifest['priorTree'].items():
        source=prior/name;require(sha(source)==r['sha256'] and source.stat().st_size==r['bytes'],'prior input '+name)
        dest=copied/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest)
        require(sha(dest)==r['sha256'],'prior byte copy '+name)
        copy_index[name]={'source':str(source),'path':str(dest.relative_to(J.out)),'sha256':r['sha256'],'bytes':r['bytes']}
    save(J.out/'prior-copy-index.json',copy_index)
    metadata={name:json.loads((copied/name).read_text()) for name in manifest['priorTree'] if name.endswith('.json')}
    require(metadata['checks.json']['incompleteOperation']=='new-tail-root-gap-reconstruction','exact original failure')
    require(metadata['failure.json']['automaticRetry'] is False,'preserved failure')
    require('operation-index.json' not in metadata,'no fabricated completed top-level stage')
    index=metadata['artifact-index.json'];restored=[]
    for name,record in index.items():
        require(record['sha256']==sha(copied/name) and record['bytes']==(copied/name).stat().st_size,'published artifact '+name)
        value=metadata[name]
        if name.endswith('-return.json') and isinstance(value,dict) and 'cancelled' in value and name!='new-tail-root-gap-reconstruction-return.json':
            stem=name[:-len('-return.json')];inputs=metadata[stem+'-input.json'];raw=metadata[stem+'-raw.json']
            require(value['cancelled']=={'text':'0','srepr':'Integer(0)'} and raw['residual']==value['raw'],'prior literal zero receipt '+name)
            # Actual returns and arguments are decoded under containment, never recomputed.
            decodedInputs=helpers['decode'](inputs);decodedReturn=helpers['decode'](value)
            require(decodedReturn['cancelled']==0,'restored zero '+name)
            require(set(decodedInputs)=={'left','right'},'complete restored operands '+name)
            restored.append({'name':stem,'status':'RESTORED_PUBLISHED_ZERO_RETURN','input':index[stem+'-input.json'],'return':record,'functionCalled':False})
    require(len(restored)==129,'all129 published zero returns')
    J.emit('restored-zero-return-index',restored)
    failed=metadata['new-tail-root-gap-reconstruction-return.json'];operands=metadata['new-tail-root-gap-reconstruction-input.json'];certificate=metadata['new-tail-root-gap-certificate.json']
    J.emit('saved-failed-reconstruction',{'input':operands,'return':failed,'polynomialCertificate':certificate,'functionReplayed':False,'residualSource':'actual published cancelled constructor, not a reconstructed or asserted result'})
    exact=exact_rational_constant(failed['cancelled']['srepr']);J.emit('saved-residual-rational-certificate',exact)
    require(exact['exactZero'] is True,'saved residual exact constant zero')
    # The original raw residual, term decomposition and arguments stay byte-identical.
    # Only the failed constant decision is completed. No polynomial extraction replay.
    cache={};used={}
    def restore(name):
        if name not in cache:cache[name]=helpers['decode'](metadata[name])
        used[name]=True;return cache[name]
    def load(name):return restore('saved/'+name)
    native=restore('first-shape-native-transport.json');domain=restore('physical-domain.json')
    trace={label:restore(label+'-new-native-trace.json') for label in ('plus','minus')}
    ref={label:restore(label+'-reference-before-guards.json') for label in ('plus','minus')}
    slot={label:restore(label+'-final-native-slot-routing.json') for label in ('plus','minus')}
    pv=restore('right-height-PV-operands.json');sheet=restore('native-middle-sheet.json');leftpv=restore('left-height-subtracted-PV.json');census=restore('retained-response-census.json')
    beta=restore('restored-direct-beta-input.json')['right'];C=restore('height-left-plus-trace-input.json')['right']
    prior_domain=restore('reused-domain-certificates.json')['domain']
    exprs=[native['proposed'],domain['frequencyContinuation'],pv['density'],sheet['transported'],leftpv['ordinarySubtractedIntegrand'],beta,C,census['heightCoefficient'],*ref['plus']['uncombined'],*[trace[label]['newTrace'] for label in trace],*[slot[label]['savedHeight'] for label in slot]]
    one=helpers['one_symbol'];context={}
    for local,saved in [('omega','reference_unrestricted_frequency'),('qi','reference_qi'),('qm','reference_qm'),('qo','reference_qo'),('k','reference_k'),('l','reference_l'),('m','reference_m'),('t','reference_t'),('cs','reference_cs'),('delta','reference_delta'),('eta','eta_bg'),('sigma','sigma_W'),('s','reference_left_height_transfer')]:context[local]=one(exprs,saved)
    def function(name):
        found={f.func for x in exprs for f in x.atoms(sp.Function) if f.func.__name__==name}
        require(len(found)==1,'unique actual saved function '+name);return found.pop()
    h=function('reference_height_hat');j=function('reference_slope_hat')
    physical=json.loads(Path(manifest['physicalInput']).read_text());values=physical['parameters']
    W=sp.Rational(values['W_0']);L=sp.Rational(values['L_W']);rho=sp.Rational(values['rho_m']);freq=domain['frequency'];edge2=sum(x*x for x in domain['edges'])
    require(W==1 and L==10 and rho==sp.Rational(1,10) and freq==3 and edge2==sp.Rational(1,20),'exact saved physical context')
    require(domain['cs']==[1,2] and domain['k/l']==[-3,3],'actual bounded domain')
    faces={label:{'T':trace[label]['newTrace'],'normal':ref[label]['normal'],'height':slot[label]['savedHeight']} for label in trace}
    source=Path(manifest['baseWorker']).read_text();tree=ast.parse(source)
    phs,psh,tracehs=ref['plus']['uncombined'];z=restore('new-tail-root-gap-input.json')['variables'][0]
    context.update(sp=sp,J=J,helpers=helpers,require=require,load=load,used=used,faces=faces,h=h,j=j,W=W,L=L,freq=freq,edge2=edge2,beta=beta,b0=prior_domain['betaMinimum'],kmin=prior_domain['kappaMinimum'],O=domain['frequencyContinuation'],phs=phs,psh=psh,tracehs=tracehs,C=C,Jdensity=pv['density'],firstref=census['heightCoefficient'],native_middle=sheet['transported'],z=z,c2=Path(manifest['c2Source']).read_text())
    # Explicit small unsaved context only: same elementary scalar/argument bindings
    # and the two ephemeral callables from their actual original assignment ASTs.
    context['mu']=rho*context['omega'];context['Q']=context['l']-context['k'];context['u']=z+12
    context['kd']=sp.sqrt((freq**2-context['delta']**2)/context['cs']**2-edge2)
    run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_science')
    for name in ('A','R'):
        found=[n for n in run.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in n.targets)]
        require(len(found)==1 and isinstance(found[0].value,ast.Lambda),'exact original callable '+name)
        exec(compile(ast.Module(body=found,type_ignores=[]),'original-ephemeral-'+name,'exec'),context)
    for name in BASE_DEFINITIONS:context[name]=base[name]
    J.emit('restored-context-and-arguments',{'savedFilesUsed':dict(used),'symbols':{k:v for k,v in context.items() if k in ['omega','qi','qm','qo','k','l','m','t','cs','delta','eta','sigma','s','z']},'functions':[h.__name__,j.__name__],'physicalParameters':values,'frequency':freq,'edges':domain['edges'],'restoredCoefficients':{'phs':phs,'psh':psh,'tracehs':tracehs,'C':C,'Jdensity':pv['density'],'firstHeight':census['heightCoefficient'],'beta':beta},'smallReconstructedContext':['mu=rho*omega','Q=l-k','u=z+12','kd=sqrt((freq^2-delta^2)/cs^2-edge2)','same A and R lambda ASTs'],'priorCompletedFunctionsCalled':False,'originalUnfinishedTailSha256':hashlib.sha256(ast.dump(ast.Module(body=tail_ast(source),type_ignores=[])).encode()).hexdigest()})
    result=tail_callable(source,context)()
    result.update(savedEvidenceContinuation=True,priorLiteralZeroReturnsRestored=129,priorFailedResidualCertifiedByExactRationalArithmetic=True,completedPrefixReplayed=False)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());g=verify_gate(args.gate,args.inputs,manifest)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False)
    J=None;result={};code=1;started=time.monotonic()
    try:
        helperSource=Path(manifest['helperSource']).read_text();baseSource=Path(manifest['baseWorker']).read_text()
        names=('require','sha','save','replace_json','posthash_records','containment','Journal','decode','one_symbol','expanded_sinh_arguments','named','function_source','assignment_source')
        nodes=[n for n in ast.parse(helperSource).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
        require({n.name for n in nodes}==set(names),'unchanged helper census')
        helpers={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'unchanged-reference-helpers','exec'),helpers)
        save(args.out/'containment.json',helpers['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        helpers.update(sp=sp,Str=Str)
        definitions=[n for n in ast.parse(baseSource).body if isinstance(n,ast.FunctionDef) and n.name in BASE_DEFINITIONS]
        require({n.name for n in definitions}==set(BASE_DEFINITIONS),'original helper definition census')
        base={'sp':sp,'require':require};exec(compile(ast.Module(body=definitions,type_ignores=[]),'unchanged-reference-definitions','exec'),base)
        class ConstantExactJournal(helpers['Journal']):
            def zero(self,name,left,right):
                previous=self.active;self.active=name
                self.emit(name+'-input',{'left':left,'right':right});raw=left-right;self.emit(name+'-raw',{'residual':raw})
                residual=sp.cancel(sp.together(raw));self.emit(name+'-return',{'raw':raw,'cancelled':residual})
                if residual!=0:
                    constructor=sp.srepr(residual);self.emit(name+'-constant-decision-input',{'constructor':constructor})
                    decision=exact_rational_constant(constructor);self.emit(name+'-constant-decision',decision)
                    require(decision['exactZero'] is True,'exact constant zero '+name)
                self.active=previous
        J=ConstantExactJournal(args.out)
        result=J.stage('saved-reference-grazing-continuation',{'manifestSha256':g['manifestSha256'],'baseWorkerSha256':sha(manifest['baseWorker']),'priorDirectory':manifest['priorDirectory']},lambda:run_continuation(manifest,J,helpers,base));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        copyfile=args.out/'prior-copy-index.json'
        if copyfile.exists():
            for v in json.loads(copyfile.read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(v['expected']!=v['actual'] for v in records.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False);save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':sys.exit(main())
