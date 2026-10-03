#!/usr/bin/env python3
"""Bounded Fourier request bank. No response integral or packet action."""
import argparse
import ast
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode','one_symbol')
G=((0,0),(1,0),(0,1),(1,1))
COMPONENTS=('NATIVE_FLAT','NATIVE_HEIGHT','NATIVE_SLOPE','NATIVE_MIXED_ITERATION','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')
ZERO={'text':'0','srepr':'Integer(0)'}


def require(v,message):
    # No bool(value): unknowns and numeric truthiness must never become proofs.
    exact_true=getattr(getattr(globals().get('sp'),'S',None),'true',None)
    if v is not True and not (exact_true is not None and v is exact_true):
        raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())


def definitions(text,names):
    selected=[n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in selected}==set(names),'exact inert helper set')
    return ast.Module(body=selected,type_ignores=[])


def verify_gate(path,manifest_path,m):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_PACKET_FOURIER_BANK','gate status')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==m['sourcePins'],'source census')
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and
        g['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'actual containment route')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source '+p)
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['launcher']==m['launcher'] and g['buildReviewRecord']==m['reviewRecordWillBe'] and g['authority']==m['executionAuthority'],'exact document routes')
    r=json.loads(Path(g['buildReviewRecord']).read_text())
    require(r['independentBuildClearance'] is True and r['allChecksPassed'] is True,'actual substantive build assessment')
    for e in ('claude','grok'):
        require(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS PACKET-ACTION FOURIER BUILD','literal build report')
    for key in ('workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256'):
        require(r[key]==g[key],'review/gate '+key)
    method=json.loads(Path(m['methodRecord']).read_text())
    require(method['jointIndependentMethodClearance'] is True and method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'exact method')
    a=json.loads(Path(g['authority']).read_text())
    require(a['scope']==m['scope']==g['scope'] and a['boundedInstrumentAuthorized'] is True and
        a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and
        a['automaticScientificRetry'] is False and a['noDeadline'] is True and g['durationLimits'] is None,'standing bounded authority')
    return g


def verify_invocation(args,g,argv):
    tail=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(list(argv)==tail and g['command'][-len(tail):]==tail,'worker argv')
    require(args.out.resolve()==Path(g['outputDirectory']).resolve(),'output route')



def rational_pair(value):
    re,im=value.as_real_imag()
    require(re.is_Rational is True and im.is_Rational is True,'finite exact rational pair')
    return [str(re),str(im)]


def family_plan(selected):
    families={}
    for a in selected:
        if a['status']!='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED':continue
        jet=a['jet']
        for role,degree in [('X',0),('Y',0),('Y',1)]:
            spec={'role':role,'argumentDerivative':degree,
                'fieldId':a[('source' if role=='X' else 'consumer')+'Transform']['coefficientId'],
                'timeOrder':jet['timeOrder'] if role=='X' else 0,
                'spatialOrders':jet['spatialOrders'] if role=='X' else [0,0,0],
                'center':'-5/2' if role=='X' else '5/2','width':'8','profileLength':'10'}
            key=json.dumps(spec,sort_keys=True)
            if key not in families:families[key]={'spec':spec,'addressIds':[]}
            families[key]['addressIds'].append(a['addressId'])
    return [families[k] for k in sorted(families)]


def run_science(m,J,ns):
    D=ns['decode'];raw={};copies={}
    for alias,r in m['savedInputs'].items():
        p=Path(r['path']);q=J.out/'saved'/alias;q.parent.mkdir(parents=True,exist_ok=True)
        require(sha(p)==r['sha256'],'saved source '+alias);shutil.copyfile(p,q)
        require(sha(q)==r['sha256'],'saved copy '+alias)
        raw[alias]=json.loads(q.read_text());copies[alias]={'source':str(p),'path':str(q.relative_to(J.out)),'sha256':r['sha256'],'bytes':q.stat().st_size}
        ns['replace_json'](J.out/'saved-copy-index.json',copies)
    prior=raw['preflight/return.json'];require(prior['status']=='PACKET_PREFLIGHT_CERTIFICATES_COMPLETE_NO_ACTION' and prior['quadratureRun'] is False,'saved completed preflight')
    require(prior['K']==27 and prior['U']==75 and prior['T']==122,'actual saved tail radii')
    require(raw['preflight/physical-plan.json']['frequency']==3 and D(raw['preflight/physical-plan.json']['cs'])==sp.sqrt(6)/2 and D(raw['preflight/physical-plan.json']['kappa'])==sp.sqrt(595)/10,'saved physical frequency and selected speed')
    require(D(raw['preflight/physical-plan.json']['centers'])==[-sp.Rational(5,2),sp.Rational(5,2)] and raw['preflight/physical-plan.json']['s']==8,'actual packet centres and width')
    selected=raw['selected/pressure-addresses.json']['selected'];original=raw['inventory/THETA_BALANCE-ordered-addresses.json']
    require(selected==[a for a in original if a['jet']['channel']=='e_W'] and len(selected)==544,'complete selected row')
    require(raw['preflight/selection.json']['addressIds']==[a['addressId'] for a in selected],'preflight actual IDs')
    fields=raw['pressure/fields.json'];certs=raw['pressure/coefficient-certificates.json'];vectors={};vector_inputs={}
    for fid,field in fields.items():
        saved=raw['preflight/field-'+fid+'-adapter-input.json'];bound=raw['preflight/field-'+fid+'-strip-bound.json']
        require(saved['original']==field and saved['certificate']==certs[fid] and saved['savedPolynomial']==raw['field/'+fid+'-polynomial.json'] and saved['inheritedProofInput']==raw['field/'+fid+'-reconstruction-input.json'],'inherited complete vector arguments')
        require(bound['coefficients']==saved['coefficientVector'] and bound['fieldRecalculated'] is False,'completed coefficient-vector join')
        pairs=[]
        for i,v in enumerate(saved['coefficientVector']):
            old=raw['preflight/field-'+fid+'-coefficient-'+str(i)+'-input.json'];ret=raw['preflight/field-'+fid+'-coefficient-'+str(i)+'-return.json']
            require(old['left']==v and ret['cancelled']==ZERO,'inherited vector proof operands/zero')
            value=D(v);pair=rational_pair(value)
            J.emit('coefficient-'+fid+'-'+str(i),{'savedValue':v,'pair':pair,'originalProofInput':old,'originalReturn':ret,'functionCalled':False})
            J.zero('coefficient-'+fid+'-'+str(i)+'-transport',value,sp.Rational(pair[0])+sp.I*sp.Rational(pair[1]));pairs.append(pair)
        vectors[fid]=list(reversed(pairs));vector_inputs[fid]=copies['preflight/field-'+fid+'-adapter-input.json']
    families=family_plan(selected)
    require(len(families)==19 and sum(f['spec']['role']=='X' for f in families)==13,'complete distinct native families')
    live={a['addressId']:a for a in selected if a['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED'}
    for fam in families:
        s=fam['spec'];fid=s['fieldId'];s['coefficientsAscending']=vectors[fid]
        for n in fam['addressIds']:
            a=live[n];route=raw['preflight/address-'+str(n)+'-adapter-input.json']
            require(route['address']==a and fields[fid]['field']==a[('source' if s['role']=='X' else 'consumer')+'Field'],'actual native field/address family join')
            if s['role']=='X':
                require((s['timeOrder'],s['spatialOrders'])==(a['jet']['timeOrder'],a['jet']['spatialOrders']),'native derivative order')
                require(raw['preflight/address-'+str(n)+'-wave-jet-input.json']['left']==a['waveMultiplier'] and raw['preflight/address-'+str(n)+'-wave-jet-return.json']['cancelled']==ZERO,'inherited actual native derivative proof')
        J.emit('family-'+str(families.index(fam)),{'spec':s,'addresses':fam['addressIds'],'sourceVectorReceipt':vector_inputs[fid]})
    unit=raw['preflight/inherited-unit-contract.json'];require(unit['sourceUnit']==[1,-1,0] and unit['forwardFactor']=='1/(2*pi)' and len(unit['actualConsumerJoins'])==4,'actual inherited Fourier and source contract')
    # Numeric reference-unit transforms only. Full response summand dimensions remain future work.
    J.emit('transform-convention',{'unit':unit,'sourceDerivativeBeforeCoefficient':True,'testNoConjugation':True,
        'X':'hat[b D_j u](k)','Y':'2*pi*hat[c v](-l)','argumentDerivative':'d/dk for X; d/dl for Y; Y prime multiplies +i*x',
        'fullSummandUnitsEstablished':False,'perRequestTolerance':'absolute 1e-12 in fixed inherited numeric reference coordinates'})
    import importlib.util
    source=Path(m['librarySource']);require(sha(source)==m['sourcePins'][str(source)],'exact numerical library')
    spec=importlib.util.spec_from_file_location('packet_fourier_runtime',source);lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib)
    import mpmath
    require(str(Path(mpmath.__file__).resolve())==m['runtimeLibrary']['initPath'] and mpmath.__version__==m['runtimeLibrary']['version'],'actual numerical runtime')
    J.emit('numeric-runtime',{'file':mpmath.__file__,'version':mpmath.__version__,'sourcePins':m['runtimeLibrary']['sourcePins']})
    store=lib.DurableStore(J.out/'fourier-evidence.sqlite');summaries=[]
    try:
        evaluator=lib.FourierEvaluator(store)
        # Fixed physical transform arguments: -K,-kappa,0,kappa,K. Y argument is -l.
        arguments=[{'base':'0','offset':'27','sign':-1},{'base':'0','offset':'sqrt(595)/10','sign':-1},
            {'base':'0','offset':'0','sign':1},{'base':'0','offset':'sqrt(595)/10','sign':1},{'base':'0','offset':'27','sign':1}]
        requests=[]
        for p0 in ['sqrt(595)/10','0']:
            for fi,family in enumerate(families):
                f=dict(family['spec']);f['carrier']=p0 if f['role']=='X' or p0=='0' else '-sqrt(595)/10'
                for argument in arguments:requests.append({'family':fi,'packetCarrier':p0,'spec':f,'argument':argument})
        require(len(requests)==190,'fixed bounded request census')
        J.emit('request-plan',{'requests':requests,'futureEveryActualRequestStillCompared':True,'noUniformTransformErrorClaim':True})
        for index,r in enumerate(requests):
            J.active='Fourier-request-'+str(index);value=evaluator.request(r['spec'],r['argument'])
            summary={'requestIndex':index,'family':r['family'],'packetCarrier':r['packetCarrier'],'request':r,'return':value['receipt']}
            summaries.append(summary);J.emit('request-'+str(index)+'-receipt',summary)
        J.active='packet-fourier';store.put('complete-request-receipts',summaries)
    finally:
        store.close()
        J.emit('numerical-journal-receipt',{'path':'fourier-evidence.sqlite','sha256':sha(J.out/'fourier-evidence.sqlite'),'bytes':(J.out/'fourier-evidence.sqlite').stat().st_size})
    return {'status':'BOUNDED_FOURIER_BANK_COMPLETE_NO_PACKET_ACTION','families':len(families),'requests':len(summaries),
        'numericalTransformsComputed':True,'responseIntegralEvaluated':False,'packetValue':None,'currentOrLoss':None,
        'scope':'190 declared source/test/Y-prime transforms only; every future actual argument still requires its own comparison',
        'futureCompleteActionRoutesRequired':True,'scientificAcceptance':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);args=p.parse_args()
    m=json.loads(args.inputs.read_text());g=verify_gate(args.gate,args.inputs,m);verify_invocation(args,g,sys.argv)
    pins={**m['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;code=1;started=time.monotonic();result={}
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(m['helperSource']).read_text(),HELPERS),'pinned-inert-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('packet-fourier',{'manifestSha256':g['manifestSha256'],'buildReviewSha256':g['buildReviewRecordSha256']},lambda:run_science(m,J,ns));code=0
    except BaseException:
        result={'status':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        index=args.out/'saved-copy-index.json'
        if index.exists():
            for r in json.loads(index.read_text()).values():pins[str(args.out/r['path'])]=r['sha256']
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(r['actual']!=r['expected'] for r in records.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result if J is None else J.encode(result));sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
