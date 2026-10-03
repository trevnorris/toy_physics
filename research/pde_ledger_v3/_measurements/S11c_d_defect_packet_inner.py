#!/usr/bin/env python3
"""Bounded collision-aware inner kernel bank. No complete packet action."""
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
    require(g['status']=='READY_FOR_ONE_PACKET_INNER_BANK','gate status')
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
        require(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS PACKET-ACTION INNER BUILD','literal build report')
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



def run_science(m,J,ns):
    D=ns['decode'];raw={};copies={}
    for alias,r in m['savedInputs'].items():
        p=Path(r['path']);q=J.out/'saved'/alias;q.parent.mkdir(parents=True,exist_ok=True)
        require(sha(p)==r['sha256'],'saved source '+alias);shutil.copyfile(p,q)
        require(sha(q)==r['sha256'],'saved copy '+alias)
        raw[alias]=json.loads(q.read_text());copies[alias]={'source':str(p),'path':str(q.relative_to(J.out)),'sha256':r['sha256'],'bytes':q.stat().st_size}
        ns['replace_json'](J.out/'saved-copy-index.json',copies)
    prior=raw['preflight/return.json'];physical=raw['preflight/physical-plan.json']
    require(prior['status']=='PACKET_PREFLIGHT_CERTIFICATES_COMPLETE_NO_ACTION' and prior['quadratureRun'] is False,'complete preflight')
    require((prior['K'],prior['U'],prior['T'])==(27,75,122),'saved tail radii')
    require(physical['frequency']==3 and D(physical['cs'])==sp.sqrt(6)/2 and D(physical['kappa'])==sp.sqrt(595)/10,'actual saved speed/frequency')
    require(raw['fourier/return.json']['status']=='BOUNDED_FOURIER_BANK_COMPLETE_NO_PACKET_ACTION' and raw['fourier/return.json']['requests']==190,'completed Fourier bank')
    for name,receipt in raw['rules/extraction.json']['records'].items():
        require(receipt['record']=='rules/'+name and receipt['sha256']==m['savedInputs']['rules/'+name+'.json']['sha256'],'full rule operand receipt')
    require(raw['rules/extraction.json']['database']==raw['fourier/journal-receipt.json'],'original complete SQLite receipt')
    inherited=[]
    def inherit(prefix):
        op=raw[prefix+'-input.json'];ret=raw[prefix+'-return.json'];J.emit('inherited-'+str(len(inherited)),{'input':op,'return':ret,'functionCalled':False})
        require(ret['cancelled']==ZERO,'literal inherited zero '+prefix);inherited.append(prefix);return op
    selected=raw['selected/pressure-addresses.json']['selected'];original=raw['inventory/THETA_BALANCE-ordered-addresses.json']
    require(selected==[a for a in original if a['jet']['channel']=='e_W'] and len(selected)==544,'entire selected native row')
    adapters=raw['preflight/numeric-factor-adapters.json'];require(len(adapters['definitions'])==20 and len(adapters['addressJoins'])==544,'all face/slot/component templates')
    for a,route in zip(selected,adapters['addressJoins']):
        label='-'.join((a['face'],a['slot'],a['component']));entry=adapters['definitions'][label]
        J.emit('address-'+str(a['addressId']),{'address':a,'route':route,'completeAdapter':entry})
        require(route=={'addressId':a['addressId'],'adapter':label},'actual template address')
        old=raw['factors/'+a['fullFactorProof']['proof']+'-operands.json']
        require(entry['original']==old['mappedAddressFactor'],'inherited exact original complete factor')
        require(old['addressNormalOriginal']==a['normalOriginal'] and old['requiredMap']==a['fullFactorProof']['completeNormalMap'],'original normal argument joins')
    for label,entry in adapters['definitions'].items():
        op=inherit('preflight/numeric-factor-'+label)
        args=raw['preflight/numeric-factor-'+label+'-arguments.json']
        require(op['left']==entry['mapped'] and op['right']==entry['template'] and args['actualCompleteFactor']==entry['original'],'complete native template proof operands')
        require(args['address']['addressId']==entry['firstAddressId'],'first original adapter address')
    for tag,alias in m['wholeDefinitionInputs'].items():
        require(raw['pressure/whole-tags.json'][tag]['savedDefinition']==raw[alias],'saved whole signature '+tag)
        require(m['savedInputs'][alias]['sha256']==raw['pressure/whole-tags.json'][tag]['sha256'],'original whole receipt')
    import importlib.util
    def load(name,path):
        require(sha(path)==m['sourcePins'][path],'pinned library')
        spec=importlib.util.spec_from_file_location(name,path);lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib);return lib
    lib=load('packet_inner_runtime',m['librarySource']);oldlib=load('packet_fourier_inert_storage',m['storageLibrary'])
    # New arithmetic adapter joins; no old source/response function is called.
    k,l,t=sp.symbols('packet_k packet_l packet_t',real=True);qi,qo,qh,qs=sp.symbols('packet_qi packet_qo packet_qh packet_qs')
    mu=sp.Rational(3,10);aa=1/(1-3*sp.I/10);beta=sp.cancel(aa*mu)
    A=lambda z:5*z/(2*sp.sinh(5*sp.pi*z))
    values=lib.kernel_components(k,l,t,qi,qo,qh,qs,A(t),A(l-k-t),aa,mu,beta,sp.Integer(1),sp.Integer(10),sp.I)
    jop=inherit('preflight/new-J-numerical-adapter');dop=inherit('preflight/new-D-numerical-adapter')
    oldJ=D(jop['right']);qm=next(s for s in oldJ.free_symbols if s.name=='packet_qm')
    J.zero('new-runtime-J-arithmetic',values[0],oldJ.xreplace({qm:qh}))
    J.zero('new-runtime-D-arithmetic-sum',sum(values[1:]),D(dop['right']))
    # Actual source parameter bindings precede all numeric defaults.
    source=raw['physical-input.json'];params=source['parameters']
    require(params['W_0']=='1' and params['L_W']=='10' and params['rho_m']=='1/10' and params['tau_A']=='1/10' and params['Lambda_A_0']=='1/100','same physical material/scale')
    J.zero('runtime-mu',mu,3*sp.Rational(params['rho_m']))
    J.zero('runtime-a',aa,sp.Rational(params['Lambda_A_0'])/(sp.Rational(params['rho_m'])**2*(1-3*sp.I*sp.Rational(params['tau_A']))))
    # Missing native h/j adapter. Native w and jet are restored, not regenerated.
    binding=raw['native/binding.json'];require(binding['physicalInput']==source,'actual native profile input')
    require(source['profiles']['w']=='(1+tanh(xi))/2','original tanh profile declaration')
    x=sp.Symbol('composition_x',real=True);w=(1+sp.tanh(x/10))/2
    lemma=D(raw['preflight/Fourier-profile-lemma.json']);hp=lemma['hPhysical'];jp=lemma['jPhysical']
    for face,sign in [('plus',1),('minus',-1)]:
        op=inherit('native/'+face+'-bound-native-height');height=D(op['left']);symbols={z.name:z for z in height.free_symbols}
        require(set(symbols)=={'eta_bg','w1_profile'},'native height free-symbol contract')
        J.zero('new-'+face+'-physical-height-join',height.subs({symbols['eta_bg']:1,symbols['w1_profile']:w},simultaneous=True),sign*hp)
    geom=D(raw['native/lower-geometry.json']);slope=geom['outwardSlopeDefinition'];symbols={z.name:z for z in slope.free_symbols}
    require(set(symbols)=={'sigma_W','w1_profile_d1'},'native slope free-symbol contract')
    # w_xi prime is the declared derivative of the actual dimensionless tanh profile;
    # physical x derivative carries 1/L. This derives only the new physical-product join.
    J.zero('new-physical-slope-join',slope.subs({symbols['sigma_W']:1,symbols['w1_profile_d1']:(1-sp.tanh(x/10)**2)/2},simultaneous=True),jp)
    J.zero('new-physical-product-derivative',hp*jp,jp/4-sp.Rational(10,8)*sp.diff(jp,x))
    hop=inherit('preflight/new-H-contact-adapter');hsop=inherit('preflight/new-H-subtracted-adapter')
    J.zero('new-runtime-H-contact',D(hop['right']),5*A(l-k)/4)
    positive=D(hsop['right']);negative=positive.xreplace({t:-t})
    J.zero('new-paired-H-density',positive+negative,5*A(t)*(A(l-k-t)-A(l-k+t))/(2*sp.I*t))
    # Exact even-cutoff cancellation for the general height action, no test evaluated.
    av,fp,fm,f0,chi=sp.symbols('A_even f_plus f_minus f_zero chi_even')
    J.zero('new-height-chi-cancellation',av*(fp-f0*chi)/t+av*(fm-f0*chi)/(-t),av*(fp-fm)/t)
    contract=raw['native/fourier-contract.json'];require(contract['profileForwardPower']==-3 and contract['invariantEdgeCoordinates']==2,'native 1/(2pi) reduced transform')
    unit=raw['preflight/inherited-unit-contract.json'];require(unit['forwardFactor']=='1/(2*pi)' and unit['sourceUnit']==[1,-1,0],'inherited source unit and measure')
    J.emit('physical-H-and-PV-lemma',{'nativeBinding':binding,'nativeGeometry':raw['native/lower-geometry.json'],'nativeFourierContract':contract,'savedProfileLemma':raw['preflight/Fourier-profile-lemma.json'],'h':hp,'j':jp,'product':hp*jp,'physicalScale':10,'Hanalytic':'jhat(Q)*(1/4-10*i*Q/8)','productToConvolutionFactor':1,'analyticTransformAndDistributionIdentities':'Assessed identities; exact local algebra above is not a CAS proof of distribution theory.','fullSummandUnitsEstablished':False})
    # New fixed-k,l tail constants derived from the published global envelopes.
    bounds=D(raw['preflight/saved-global-bound-inputs.json']);b=bounds['domain']['betaMinimum'];require(b==sp.Rational(3000,11101),'actual beta lower bound')
    require(bounds['profile']['productBound']=='121 exp(-|t|)' and bounds['profile']['L']==10 and bounds['domain']['aMagnitudeUpper']==1,'actual common-sheet profile/domain bound')
    K=sp.Integer(prior['K']);T=sp.Integer(prior['T']);ex=sp.Rational(1,2)**T
    Jtail=sp.Rational(9,40)*121/b**3*2*ex*(T*T+2*T+2+3*K*(T+1)+2*K*K)
    Dref=sp.Rational(3,4)*121/b**2*2*ex*K*(T+1+2*K)
    Dquad=sp.Rational(3,4)*121/b**2*2*ex*(K*K+9)
    tails={'J':Jtail,'Dreflected':Dref,'Dheight':Dref,'Dquadratic':Dquad,'Dsum':2*Dref+Dquad,'H':sp.Rational(55,3)*ex}
    J.emit('new-inner-tail-certificates',{'K':K,'T':T,'b':b,'actualBounds':raw['preflight/saved-global-bound-inputs.json'],'tails':tails,'inequalities':['T>=K+4, kappa<3 imply both internal depth moduli>1 outside |t|<=T','quadrant |qh+qi|>=|qh| and >=|qi|; each q+beta has modulus>=b','J uses (|t|+K)(|t|+2K); direct uses SUM K(|t|+2K),K(|t|+2K),K^2+9','2^-T overestimates exp(-T), with exact polynomial exponential-tail moments'],'evaluatedIntegral':False,'noCancellationPaysTail':True})
    require(T>=K+4 and all(v.is_positive is True and v<sp.Rational(1,10**11) for v in tails.values()),'positive per-component tail certificates')
    plan=lib.point_plan();J.emit('fixed-point-plan',{'points':plan,'kappa':physical['kappa'],'T':122,'enlargedT':124,'noExactExternalGrazing':True,'noOuterPacketIntegral':True})
    # Join the runtime final mixed/direct pressure and normal arithmetic to all 8 native templates.
    hv,jv,dv=sp.symbols('packet_H packet_J packet_D')
    mixed=-sp.I*mu*k*qo*hv/((qo+beta)*(qi+beta))+jv
    for face,sign in [('plus',1),('minus',-1)]:
        for slot in ['pressure','normal']:
            for component,expr in [('NATIVE_MIXED_ITERATION',mixed),('INHERITED_DIRECT_WHOLE_OFF_DIAGONAL',dv)]:
                key='-'.join([face,slot,component]);expected=expr*(sp.I*sign*qo if slot=='normal' else 1)
                J.zero('new-final-runtime-'+key,expected,D(adapters['definitions'][key]['template']))
    import mpmath
    require(str(Path(mpmath.__file__).resolve())==m['runtimeLibrary']['initPath'] and mpmath.__version__==m['runtimeLibrary']['version'],'actual pinned numerical runtime')
    rules={n:raw['rules/'+n+'.json'] for n in ['A-GL24','A-GL48','B-G7-K15']}
    store=oldlib.DurableStore(J.out/'inner-evidence.sqlite')
    try:
        evaluator=lib.InnerEvaluator(store,oldlib,rules);J.active='inner-numerical-bank'
        result=evaluator.run(plan)
    finally:
        store.close();J.emit('numerical-journal-receipt',{'path':'inner-evidence.sqlite','sha256':sha(J.out/'inner-evidence.sqlite'),'bytes':(J.out/'inner-evidence.sqlite').stat().st_size})
    return {'status':'BOUNDED_INNER_KERNEL_BANK_COMPLETE_NO_PACKET_ACTION',**result,'inheritedZeroReturns':len(inherited),'responseIntegralsEvaluated':True,'packetActionEvaluated':False,'currentOrLoss':None,'fullActionAccuracyClaim':False,'scientificAcceptance':False}

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
        result=J.stage('packet-inner',{'manifestSha256':g['manifestSha256'],'buildReviewSha256':g['buildReviewRecordSha256']},lambda:run_science(m,J,ns));code=0
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
