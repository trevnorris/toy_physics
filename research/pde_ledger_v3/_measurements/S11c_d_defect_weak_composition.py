#!/usr/bin/env python3
"""Global weak-composition certificates from saved operands; no convolution/solve."""
import argparse
import ast
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import resource
import shutil
import sys
import time
import traceback
from types import SimpleNamespace

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS = ('require','sha','save','replace_json','posthash_records','containment',
           'Journal','decode','one_symbol','named','function_source','assignment_source',
           'expanded_sinh_arguments')
TEXT_HELPERS = ('require','text_sha','literal_record','tuple_arguments','literal_key',
                'named','selected_case')
G = ((0,0),(1,0),(0,1),(1,1))
SLOTS = ('delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus')
ROWS = ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')


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


def verify_gate(path, manifest_path, manifest):
    gate = json.loads(Path(path).read_text())
    verify_helper_paths(gate)
    require(gate['status']=='READY_FOR_ONE_WEAK_COMPOSITION_INSTRUMENT','gate status')
    require(gate['workerSha256']==sha(__file__) and
            gate['manifestSha256']==sha(manifest_path), 'worker/manifest')
    require(gate['sourcePins']==manifest['sourcePins'], 'source census')
    for p,h in gate['sourcePins'].items():
        require(sha(p)==h,'source '+p)
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(gate[key])==gate[key+'Sha256'],'gate '+key)
    require(gate['launcher']==manifest['launcher'], 'launcher route')
    require(gate['buildReviewRecord']==manifest['reviewRecordWillBe'],'manifest review route')
    review=json.loads(Path(gate['buildReviewRecord']).read_text())
    require(review['independentBuildClearance'] is True and
            review['methodAssessed'] is True and review['allChecksPassed'] is True,
            'substantive implementation and corrected-method assessment')
    for engine in ('claude','grok'):
        require(review['reports'][engine]['literalVerdict']==
                'CLEAR FOR THIS GLOBAL WEAK-COMPOSITION BUILD','literal build verdict')
    for key in ('workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256'):
        require(review[key]==gate[key], 'review/gate '+key)
    require(review['methodSha256']==sha(manifest['methodPath']), 'exact corrected method')
    authority=json.loads(Path(gate['authority']).read_text())
    require(authority['boundedInstrumentAuthorized'] is True and
            authority['automaticScientificRetry'] is False, 'standing bounded authority')
    require(gate['durationLimits'] is None and gate['scientificRunsAuthorized']==1 and
            gate['scope']==manifest['scope'], 'one bounded no-deadline job')
    method_record=json.loads(Path(manifest['methodRecord']).read_text())
    require(method_record['jointIndependentMethodClearance'] is True and
            method_record['methodSha256']==sha(manifest['methodPath']),'method assessment and exact source')
    return gate


def select_control_address(addresses, component, face, slot, source_grade, consumer_grade, predicate):
    candidates=[a for a in addresses if a['component']==component and a['face']==face and a['slot']==slot
        and a['sourceGrade']==list(source_grade) and a['consumerGrade']==list(consumer_grade)
        and a['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED' and predicate(a)]
    require(bool(candidates),'applicable saved control address')
    return min(candidates,key=lambda a:a['addressId'])


def coefficient_polynomial(expr, x, T):
    """New small polynomial certificate, never a full source-domain Poly."""
    substituted=expr.xreplace({sp.tanh(x/10):T})
    require(not substituted.atoms(sp.Function) and substituted.free_symbols<= {T},'polynomial field leaves')
    numerator,denominator=sp.fraction(sp.cancel(substituted))
    require(not denominator.free_symbols and denominator.is_finite is True and denominator.is_zero is False,
            'exact finite constant field denominator')
    p=sp.Poly(numerator,T)
    require(all(not v.free_symbols and v.is_finite is True for v in p.all_coeffs()),'finite exact field coefficients')
    return p,denominator


def run_science(manifest,J,ns,exact_nonzero):
    raw={};copies={};used=set();D=ns['decode'];one=ns['one_symbol']
    for alias,v in manifest['savedInputs'].items():
        source=Path(v['path']);target=J.out/'saved'/alias;target.parent.mkdir(parents=True,exist_ok=True)
        require(sha(source)==v['sha256'],'saved input '+alias);shutil.copyfile(source,target)
        require(sha(target)==v['sha256'],'saved byte copy '+alias)
        copies[alias]={'source':str(source),'path':str(target.relative_to(J.out)),'sha256':v['sha256'],'bytes':target.stat().st_size}
        raw[alias]=json.loads(target.read_text())
    save(J.out/'saved-copy-index.json',copies)
    def take(alias):used.add(alias);return raw[alias]
    def load(alias):return D(take(alias))
    inherited=[]
    def inherit_zero(alias):
        r=take(alias);require(r.get('cancelled')=={'text':'0','srepr':'Integer(0)'},'published literal zero '+alias)
        inherited.append({'alias':alias,'receipt':copies[alias],'functionCalled':False,'status':'RESTORED_PUBLISHED_ZERO_RETURN'})
    # No old source/row/grade reconstruction is called. Exact metadata joins attach new bounds.
    fields=take('inventory/fields.json');coverage=take('inventory/grade-coverage.json')
    require(len(fields)==34 and len(coverage)==320,'completed field/grade census')
    addresses=[]
    for row in ROWS:
        addresses.extend(take('inventory/'+row+'-ordered-addresses.json'))
        partition=take('inventory/'+row+'-full-native-partition.json')
        require(not partition['unknownPressureAtoms'] and partition['localPartUnchanged'] is True,'saved native partition')
        require(all(h['name'] in SLOTS for c in partition['children'] for h in c['hits']),'native pressure slot set')
    require(sorted(a['addressId'] for a in addresses)==list(range(13260)),'complete address ID census')
    require({a['slot'] for a in addresses}=={'pressure','normal'},'no tangential response consumer')
    byid={a['addressId']:a for a in addresses}
    scale=take('inventory/native-profile-scale-join.json');require(scale['savedLength']['srepr']=='Integer(10)' and scale['physicalLength']=='10' and scale['declaredLength']==10,'actual L=10')
    speed=take('inventory/prebinding-native-speed-inventory.json');require(speed['noNativeSpeedSymbols'] and not speed['legacyCsRetuned'],'prebinding speed independence')
    require(all(not v['speedSymbols'] for v in speed['inventory']),'actual empty speed hits')
    for face in ('plus','minus'):
        st=take('inventory/'+face+'-native-source-join-status.json')
        require(st['densityRestoredAndJoined'] and st['chemicalAmplitudeJoined'] and st['inheritedNormalizationJoined'],'saved actual native source joins')
        inherit_zero('inventory/'+face+'-raw-stage2-live-density-join-return.json')
    transforms={k:D(v['field']) for k,v in fields.items()}
    x=one(list(transforms.values()),'composition_x');T=sp.Symbol('weak_tanh_variable',real=True)
    polys={};fieldcert={}
    for fid,expr in transforms.items():
        J.emit('field-'+fid+'-operands',{'fieldId':fid,'original':expr,'x':x,'T':T,'profileScale':scale})
        poly,den=coefficient_polynomial(expr,x,T);pexpr=poly.as_expr()/den
        coeff=poly.all_coeffs();degree=0 if poly.is_zero else int(poly.degree())
        J.emit('field-'+fid+'-polynomial',{'polynomial':poly.as_expr(),'denominator':den,'coefficients':coeff,'degree':degree})
        J.zero('field-'+fid+'-reconstruction',expr,pexpr.subs(T,sp.tanh(x/10)))
        nextpoly=sp.expand((1-T**2)*sp.diff(pexpr,T)/10)
        J.zero('field-'+fid+'-first-derivative',sp.diff(expr,x),nextpoly.subs(T,sp.tanh(x/10)))
        bound=sum(sp.Abs(v/den) for v in coeff)
        J.emit('field-'+fid+'-derivative-class',{'recurrence':'P[n+1]=(1-T**2)*dP[n]/dT/10','P0':pexpr,'P1':nextpoly,
            'zerothDerivativeBound':bound,'argument':'Induction: finite polynomials at each n; |T|<=1 bounds every derivative by a finite coefficient sum.',
            'allOrdersComputed':False,'analyticInductionNotMachineTheorem':True,'schwartzMultiplierClass':'smooth with all derivatives bounded'})
        polys[fid]=pexpr;fieldcert[fid]={'degree':degree,'polynomial':pexpr,'constant':degree==0,'derivativeClass':'bounded at every finite order by assessed induction'}
    J.emit('all-coefficient-certificates',fieldcert)
    # Restore complete factor proofs once; attach the new envelope class at every address.
    proofset=set();route=[];jets={};component_degrees={'NATIVE_FLAT':(0,0),'NATIVE_HEIGHT':('PV_Pk^(3/2)','PV_Pk^(3/2)'),
        'NATIVE_SLOPE':(1,2),'NATIVE_MIXED_ITERATION':(2,3),'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':(2,3)}
    for a in addresses:
        source=fields[a['sourceTransform']['coefficientId']]['field'];consumer=fields[a['consumerTransform']['coefficientId']]['field']
        require(source==a['sourceField'] and consumer==a['consumerField'],'actual field IDs at address')
        require(a['responseMap']['frequency']==3 and a['responseMap']['positiveRegulatorContinuation'] is False,'real source/response address')
        flat=a['responseGrade']==[0,0]
        require(a['responseInputDepth']==('q(l)' if flat else 'q(k)') and a['responseOutputDepth']=='q(l)','depth role')
        require(a['sourceTransform']['transfer']==('l-p' if flat else 'k-p') and a['consumerTransform']['transfer']=='r-l','weak Fourier order')
        proof=a['fullFactorProof']['proof']
        if proof not in proofset:
            for tail in ['-full-mapped-residual-return.json','-normal-source-join-return.json']:inherit_zero('inventory/factors/'+proof+tail)
            proofset.add(proof)
        operands=take('inventory/factors/'+proof+'-operands.json')
        map_keys=('sha256','original','mapped','map','symbolAssumptions','flatSupport','frequency','positiveRegulatorContinuation')
        require(all(operands['actualResponseMap'][key]==a['responseMap'][key] for key in map_keys),'saved factor map actual arguments')
        require(operands['addressNormalOriginal']==a['normalOriginal'] and operands['requiredMap']==a['fullFactorProof']['completeNormalMap'],'actual inherited normal role map')
        require(a['component'] in component_degrees and not a['wholeValueEvaluated'],'whole object inventory')
        if a['component']=='INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':
            require(a['responseGrade']==[1,1] and a['sourceGrade']==a['consumerGrade']==[0,0],'direct grade isolation')
            require(a['responseCoefficient']['text'].startswith('Dwhole_'+a['face']+'('),'actual whole direct tag')
        jets[a['jet']['name']]=a['jet']
        route.append({'addressId':a['addressId'],'face':a['face'],'slot':a['slot'],'component':a['component'],
            'gradeTriple':[a[k] for k in ('consumerGrade','responseGrade','sourceGrade')],
            'sourceFieldId':a['sourceTransform']['coefficientId'],'consumerFieldId':a['consumerTransform']['coefficientId'],
            'sourceJet':a['jet'],'status':a['status'],'inheritedFactorProof':proof,
            'newEnvelopeClass':component_degrees[a['component']][a['slot']=='normal'],'kernelValueComputed':False})
    J.emit('weak-address-coverage',route);J.emit('actual-source-jets',list(jets.values()))
    # Constants and source law are new global-certificate inputs; old compact proofs remain saved.
    omega,qi,qh,qs,qo,qm=sp.symbols('weak_Omega weak_qi weak_qh weak_qs weak_qo weak_qm')
    k,l,t=sp.symbols('weak_k weak_l weak_t',real=True);delta=sp.Symbol('weak_delta',nonnegative=True)
    cs=sp.Symbol('weak_cs',positive=True);mu=omega/10;aa=1/(1-sp.I*omega/10);beta=omega/(10-sp.I*omega)
    A=lambda s:10*s/(4*sp.sinh(5*sp.pi*s));j=lambda s:5*A(s)
    rc=load('saved/reference/retained-response-census.json');jr=load('saved/reference/right-height-PV-operands.json');dr=load('saved/direct/closed-density.json');hr=load('saved/reference/left-height-subtracted-PV.json')
    def lift(expr,names):
        actual={a:names[a.name] for a in expr.free_symbols if a.name in names}
        J.emit('lift-'+str(len(J.artifacts)),{'original':expr,'map':[[a,b] for a,b in actual.items()],'unmapped':[a for a in expr.free_symbols if a not in actual]})
        require(expr.free_symbols<=set(actual),'complete new envelope symbol map')
        return expr.subs(actual,simultaneous=True)
    refmap={'reference_unrestricted_frequency':omega,'reference_qi':qi,'reference_qo':qo,'reference_qm':qm,'reference_k':k,'reference_l':l,'reference_t':t,'reference_left_height_transfer':t}
    dmap={'grazing_unrestricted_frequency':omega,'grazing_qi':qi,'grazing_qo':qo,'grazing_qh':qh,'grazing_qs':qs,'k':k,'grazing_output':l,'grazing_transfer':t}
    actualJ=lift(jr['density'],refmap);actualD=lift(dr['density'],dmap);actualBc=lift(dr['Bc'],dmap)
    Jformula=aa*mu**2*sp.Rational(5,2)*A(t)*A(l-k-t)*(k+t)*(2*k+t)*qi/(qm*(qo+beta)*(qm+beta)*(qi+beta)*(qm+qi))
    Bformula=-sp.I*mu/((qi+beta)*(qo+beta))*(k*(2*l-t)/(qs+qo)+k*(t+2*k)*qi/(qh*(qh+qi))+qi**2/qh)
    Dformula=sp.Rational(5,2)/sp.I*A(t)*A(l-k-t)*Bformula
    J.emit('actual-global-density-operands',{'savedJ':actualJ,'savedD':actualD,'savedBc':actualBc,'Jformula':Jformula,'Dformula':Dformula,'BcFormula':Bformula,'Jroute':'qm=q(k+t)','Droutes':['qh=q(k+t)','qs=q(l-t)']})
    argument_maps=load('inventory/actual-whole-density-argument-maps.json')
    source_depth=argument_maps['depthFunction'];source_p=source_depth.args[0];source_cs=argument_maps['cs'];depth_definition=argument_maps['depthDefinition']
    require(argument_maps['noValueAtGrazing'] is True and isinstance(depth_definition,sp.Piecewise),'saved outgoing depth domain')
    J.emit('actual-global-depth-input',{'sourceFunction':source_depth,'definition':depth_definition,'cs':source_cs,'newGlobalLaw':'sqrt(Omega^2/cs^2 - 1/20 - p^2), first quadrant; real boundary limit', 'positiveDeltaSourceBinding':False})
    for branch in range(2):
        J.zero('actual-real-depth-radicand-'+str(branch),depth_definition.args[branch][0]**2,9/source_cs**2-sp.Rational(1,20)-source_p**2)
    for label,names in [('direct',{'grazing_qi':'input','grazing_qo':'output','grazing_qh':'height','grazing_qs':'reflected'}),('iteration',{'reference_qi':'input','reference_qo':'output','reference_qm':'height'})]:
        entry=argument_maps[label];Kactual=one([entry['mapped']],'composition_k');Lactual=one([entry['mapped']],'composition_l');transfer=entry['boundVariable']
        expected_arguments={'input':Kactual,'output':Lactual,'height':Kactual+transfer,'reflected':Lactual-transfer}
        actual_map={a.name:b for a,b in entry['map']}
        J.emit(label+'-global-root-route-input',{'actualMap':entry['map'],'requiredArguments':expected_arguments,'boundVariable':transfer})
        for name,role in names.items():
            value=actual_map[name];require(value.func==source_depth.func and len(value.args)==1,'native outgoing function route')
            J.zero(label+'-'+role+'-global-root-route',value.args[0],expected_arguments[role])
    J.zero('new-flat-global-coefficient',lift(rc['flat'],refmap),mu/(qi+beta))
    for key,tag,expected in [('heightCoefficient','reference_height_hat',-sp.I*mu*qi/(qo+beta)),('slopeCoefficient','reference_slope_hat',mu*k/((qi+beta)*(qo+beta)))]:
        expr=rc[key];tags=[a for a in expr.atoms(sp.Function) if a.func.__name__==tag];require(len(tags)==1,'actual first-order profile tag')
        J.zero('new-'+key+'-global-coefficient',lift(sp.expand(expr).coeff(tags[0]),refmap),expected)
    inheritedH=one([rc['mixedIteration']],'whole_height_slope_convolution')
    J.zero('new-C-global-coefficient',lift(sp.expand(rc['mixedIteration']).coeff(inheritedH),refmap),-sp.I*mu*k*qo/((qo+beta)*(qi+beta)))
    J.sinh_zero('new-J-envelope-source-join',actualJ,Jformula);J.sinh_zero('new-D-envelope-source-join',actualD,Dformula);J.zero('new-Bc-envelope-source-join',actualBc,Bformula)
    J.sinh_zero('new-H-contact-source-join',lift(hr['contact'],refmap),j(l-k)/4)
    J.sinh_zero('new-H-subtraction-source-join',lift(hr['ordinarySubtractedIntegrand'],refmap),A(t)*(j(l-k-t)-j(l-k))/(2*sp.I*t))
    b=sp.Rational(3000,11101);amin=sp.sqrt(879)/20;Cq=36/sp.sqrt(amin)
    dom=load('saved/direct/domain-bound-certificate.json');oldbeta=lift(dom['beta'],{'grazing_delta':delta})
    J.zero('native-beta-global-law',oldbeta,beta.subs(omega,3+sp.I*delta))
    J.zero('inherited-beta-minimum',dom['betaMinimum'],b);J.zero('inherited-kappa-minimum',dom['kappaMinimum'],amin)
    K,L,Z=sp.symbols('weak_abs_k weak_abs_l weak_abs_t',nonnegative=True);P=1+K+L
    def nonnegative_poly(name,larger,smaller,variables):
        gap=sp.expand(larger-smaller);coeff=sp.Poly(gap,*variables).coeffs()
        J.emit(name,{'larger':larger,'smaller':smaller,'expandedGap':gap,'nonnegativeVariables':list(variables),'coefficients':coeff})
        require(all(c.is_nonnegative is True for c in coeff),'polynomial nonnegative gap '+name)
    nonnegative_poly('global-q-triangle-gap',(K+4)**2,K**2+sp.Rational(453,50),(K,))
    nonnegative_poly('J-numerator-envelope',2*P**2*(1+Z)**2,(K+Z)*(2*K+Z),(K,L,Z))
    nonnegative_poly('D-reflected-numerator-envelope',2*P**2*(1+Z),K*(2*L+Z),(K,L,Z))
    nonnegative_poly('D-height-numerator-envelope',18*P**2*(1+Z),K*(Z+2*K)+16*(1+K)**2,(K,L,Z))
    nonnegative_poly('normal-polynomial-growth',4*P,4*(1+L),(K,L))
    nonnegative_poly('height-root-difference-growth',3*(1+K),2*K+1,(K,))
    J.zero('beta-real-law',sp.re(beta.subs(omega,3+sp.I*delta)).expand(complex=True),30/((10+delta)**2+9))
    J.zero('beta-imag-law',sp.im(beta.subs(omega,3+sp.I*delta)).expand(complex=True),(9+delta*(10+delta))/((10+delta)**2+9))
    J.emit('global-parameter-domain',{'cs':[1,2],'delta':[0,sp.Rational(1,10)],'externalMomenta':'all real','beta':beta,'betaMinimum':b,'kappaMinimum':amin,'kappaMaximum':3,
        'muMagnitudeUpper':sp.Rational(2,5),'aMagnitudeUpper':1,'betaMagnitudeUpper':sp.Rational(2,5),
        'qMagnitudeUpper':'|p|+4','proof':'|Omega|^2<=9.01; |Omega^2/cs^2-1/20|<=9.06; quadrant inequalities inherited with actual native law.',
        'endpointMajorant':'a_*^(-1/2) SUM |p-s*kappa_delta|^(-1/2), s=+-1; moving real comparison endpoints; a.e. only'})
    J.zero('beta-minimum-endpoint',30/((10+sp.Rational(1,10))**2+9),b)
    J.zero('kappa-minimum-endpoint',(9-sp.Rational(1,100))/4-sp.Rational(1,20),amin**2)
    J.emit('global-profile-envelope',{'A':A(t),'removableValue':1/(2*sp.pi),'L':10,'globalBound':'11 exp(-|t|)',
        'smallArgument':'|t|<=1: A<=1/(2pi)<=11/e','largeArgument':'|t|>=1: A<=10|t|exp(-5pi|t|)<=11exp(-|t|)',
        'noQExponentialPrefactor':True,'productBound':'121 exp(-|t|)','analyticInequalitiesAssessed':True})
    w=(1+Z)**2*sp.exp(-Z);primitive=-(Z**2+4*Z+5)*sp.exp(-Z)
    J.zero('weight-antiderivative-identity',sp.diff(primitive,Z),w)
    J.zero('weight-derivative-identity',sp.diff(w,Z),(1+Z)*(1-Z)*sp.exp(-Z))
    J.zero('weight-half-line-mass',-primitive.subs(Z,0),sp.Integer(5))
    J.emit('shift-uniform-square-root-bound',{'weight':w,'supremumUpper':2,'fullMass':10,'nearUnitBallBound':8,'farBound':10,'perRootBound':18,
        'inverseDepthBound':Cq,'routeCenters':['-k +/- kappa_delta','l +/- kappa_delta'],
        'proof':'near root integrate |r|^-1/2 over (-1,1)=4; far factor<=1. Each inverse depth uses TWO roots. Sup w=4/e<2 by e>2.',
        'responseIntegralEvaluated':False})
    J.emit('H-global-bound',{'savedContact':hr['contact'],'savedSubtraction':hr['ordinarySubtractedIntegrand'],'contactBound':5/(8*sp.pi),
        'nearBound':25/(4*sp.pi),'tailBound':10*sp.exp(-5*sp.pi)/sp.pi**2,'looseUpper':100,
        'rationalCoarseUpper':sp.Rational(135,8),'positiveGapTo100':sp.Rational(665,8),'argument':'pi>1 and exp(-5pi)<1; bounds on j,j-prime independent of Q; no |Q|<=6 assumption'})
    envelopes={'J':sp.Rational(4,5)*121*Cq/b**3,'D':18*121*2*Cq/b**2}
    J.emit('global-whole-kernel-envelopes',{'constants':envelopes,'pressureWeight':'(1+|k|+|l|)^2','normalWeight':'(1+|k|+|l|)^3',
        'normalExtraConstant':4,'Hbound':100,'CpressureDegree':1,'CnormalDegree':2,'noComputedResponseIntegral':True,
        'continuityArgument':'On each compact external set: moving-endpoint uniform absolute continuity + compact-external tail + a.e. convergence. Global polynomial envelope then dominates Schwartz pairing.',
        'singleGlobalPointwiseDominantAssumed':False,'machineMeasureTheoryProof':False})
    # NEW global PV coefficients are joined to both saved native signs.
    B=-sp.I*mu*qi/(qo+beta);Bi=-sp.I*mu*qi/(qi+beta)
    J.zero('height-difference-identity',B-Bi,-sp.I*mu*qi*(qi-qo)/((qo+beta)*(qi+beta)))
    for face,sign in [('plus',1),('minus',-1)]:
        actualNormal=lift(rc['normalJet'][face],refmap);J.zero(face+'-actual-normal-for-global-bound',actualNormal,sp.I*sign*qo)
        BN=sp.I*sign*qo*B;BNi=sign*mu*qi**2/(qi+beta)
        J.zero(face+'-normal-height-difference',BN-BNi,sign*mu*qi*beta*(qo-qi)/((qo+beta)*(qi+beta)))
        J.emit(face+'-global-height-PV-certificate',{'pressureCoefficient':B,'normalCoefficient':BN,'beta':beta,
            'pressureSupConstant':8/(5*b),'pressureHolderConstant':8*sp.sqrt(3)/(5*b**2),
            'normalSupConstant':sp.Rational(8,5),'normalHolderConstant':16*sp.sqrt(3)/(25*b**2),
            'weight':'Pk^(3/2) for variation, Pk for magnitude','subtraction':'A(Q)[f_k(Q)-f_k(0) 1_|Q|<=1]/Q plus W/4 contact',
            'normalCoefficientBoundedDirectly':True,'qTimesTestAssumedC1':False,'testSpace':'Schwartz','sourceFrequencyHeldReal':3})
    fourier=take('inventory/fourier-and-unit-provenance.json')
    require(fourier['oneDimensionalForwardNormalization']=='1/(2*pi)' and fourier['measures']==['dl','dl dk'] and not fourier['evaluatedTransform'],'native Fourier convention')
    J.emit('new-weak-duality-and-order',{'inheritedConvention':fourier['contract'],'source':'X(k)=hat[b D_j u](k)','consumer':'Y(l)=2pi hat[c v](-l)',
        'pairing':'integral Y(l)F(l,k)X(k) dl dk; bilinear, not power','constantDeltaTest':'F=delta(l-k) gives integral c v b D_j u dx',
        'flatSupport':'k=l applies to whole flat pole and normal prefactor','sourceDerivativeBeforeMultiplication':True,
        'sourceDerivativeVariable':'original p','normalVariable':'response output l','nativeEpsilonOnce':True,
        'argument':'Fourier inversion and bilinear duality with forward 1/(2pi); coefficient multiplication preserves S by all-field induction.'})
    # Restore identical q(r) controls; no function call or scalar recomputation.
    oldcontrols=take('inventory/responsive-formal-controls.json');restored=[]
    for r in oldcontrols:
        if not r['name'].endswith('normal-q-l-to-r'):continue
        op=take('inventory/'+r['name']+'-control-operands.json');a=byid[op['context']['addressId']]
        J.emit(r['name']+'-reused-arguments',{'actualAddress':a,'savedControlOperands':op,'savedReturn':r})
        require(a['slot']=='normal' and a['component']=='NATIVE_SLOPE' and a['consumerGrade']==[1,0],'reused normal control address')
        require(a['consumerTransform']['coefficientId']==op['context']['consumerTransform']['coefficientId'] and not fieldcert[a['consumerTransform']['coefficientId']]['constant'],'nonconstant actual consumer')
        require(a['jet']==op['context']['sourceJet']['spec'] and a['sourceField']==op['context']['sourceJet']['field'],'actual source arguments for reused normal control')
        require(r['certificate']['finiteNonzero'] and r['formalTagCoefficientOnly'] and op['physicalConvolutionNotEvaluated'],'old control scoped result')
        frac=take('inventory/'+r['name']+'-movement-fraction.json')
        require(all(all(row) for row in frac['finite']) and all(any(row) for row in frac['signedNonzero']),'actual prior finite signed-component certificate')
        for tail in ['-numerator-components-return.json','-denominator-components-return.json','-fraction-reconstruction-return.json']:
            inherit_zero('inventory/'+r['name']+'-movement'+tail)

        restored.append({'name':r['name'],'addressId':a['addressId'],'functionCalled':False,'movement':r['movement'],'status':'RESTORED_PRIOR_CONTROL_RETURN'})
    require(len(restored)==4,'all four identical normal controls reused');J.emit('restored-normal-controls',restored)
    controls=[]
    def control(name,baseline,corrupt,context):
        movement=sp.cancel(corrupt-baseline);J.emit(name+'-operands',{'baseline':baseline,'corrupt':corrupt,'movement':movement,'context':context,'formalCoefficientOnly':True})
        cert=exact_nonzero(J,name+'-movement',movement);controls.append({'name':name,'movement':movement,'certificate':cert,'formalCoefficientOnly':True})
    for face in ('plus','minus'):
        a=select_control_address(addresses,'NATIVE_FLAT',face,'pressure',(1,0),(0,0),lambda v:
            v['jet']['spatialOrders'][0]>0 and fieldcert[v['sourceTransform']['coefficientId']]['degree']==1)
        field=transforms[a['sourceTransform']['coefficientId']];n=a['jet']['spatialOrders'][0];U=sp.Function('weak_trial')(x)
        comm=sp.expand(sp.diff(field*U,x,n)-field*sp.diff(U,x,n));expected=sum(sp.binomial(n,jj)*sp.diff(field,x,jj)*sp.diff(U,x,n-jj) for jj in range(1,n+1))
        J.emit(face+'-Leibniz-control-input',{'address':a,'actualField':field,'actualDerivativeOrder':n,'commutator':comm,'jetFormula':expected,'notFourierValue':True})
        J.zero(face+'-Leibniz-reconstruction',comm,expected)
        local_movement=n*sp.diff(field,x).subs(x,0)
        flat_actual=lift(D(a['responseOriginal']),refmap);flat_point=flat_actual.subs(qo,sp.Rational(3,2))
        actual_consumer=D(a['consumerField']);require(not actual_consumer.free_symbols,'constant addressed pressure consumer')
        movement=local_movement*flat_point*actual_consumer
        control(face+'-Leibniz-interchange',sp.S.Zero,movement,{'addressId':a['addressId'],'coefficient':'independent trial derivative of order n-1; coefficient field at x=0',
            'fieldIsNotFourierTransform':True,'localCoefficientMovement':local_movement,'actualFlatResponse':flat_actual,'responsePoint':flat_point,
            'actualConsumer':actual_consumer,'responsePointDomain':{'l':2,'csSquared':'10/7','qo':'3/2'},'formalCoefficientSensitivityOnly':True})
        haddr=select_control_address(addresses,'NATIVE_MIXED_ITERATION',face,'pressure',(0,0),(0,0),lambda v:v['jet']['name']=='e_W')
        expr=D(haddr['responseOriginal']);Hs=[a for a in expr.atoms(sp.Function) if a.func.__name__=='Hwhole'];require(len(Hs)==1,'actual H tag')
        actualC=sp.expand(expr).coeff(Hs[0]);C=-sp.I*mu*k*qo/((qo+beta)*(qi+beta))
        mappedC=lift(actualC,refmap);J.zero(face+'-H-contact-actual-address-coefficient',mappedC,C.subs(omega,3))
        src=D(haddr['sourceField']);consumer=D(haddr['consumerField']);require(not src.free_symbols and not consumer.free_symbols,'constant control source and consumer')
        cpoint=mappedC.subs({k:sp.Rational(3,2),qo:sp.Rational(3,2),qi:sp.Integer(2)},simultaneous=True)*src*consumer/4
        control(face+'-omit-H-contact',cpoint,sp.S.Zero,{'addressId':haddr['addressId'],'point':{'k':'3/2','l':'2','csSquared':'10/7'},'perUnitFactor':'j(1/2)>0 by A positive on positive real argument','retainedContact':'W j(Q)/4','notComputedH':True})
    # New reflected-route failure certificate acts on the actual saved Bc first term.
    cs2=sp.Rational(10,7);kp=sp.Rational(3,2);lp=sp.Integer(2);tp=sp.Rational(1,2)
    rad=lambda p:9/cs2-sp.Rational(1,20)-p**2
    qip=sp.sqrt(rad(kp));qop=sp.sqrt(rad(lp));qhp=sp.sqrt(rad(kp+tp));qsp=sp.sqrt(rad(lp-tp))
    physical={omega:sp.Integer(3),k:kp,l:lp,t:tp,qi:qip,qo:qop,qh:qhp,qs:qsp}
    baseline=actualBc.subs(physical,simultaneous=True);wrong=actualBc.subs(qs,qh).subs(physical,simultaneous=True)
    J.emit('reflected-root-control-input',{'actualSavedBc':actualBc,'map':[[a,b] for a,b in physical.items()],
        'correctQsMomentum':lp-tp,'wrongQhMomentum':kp+tp,'actualQs':qsp,'wrongQh':qhp,'csSquared':cs2,'actualDomain':rad(lp-tp),
        'corruptedRadicandResidual':qhp**2-rad(lp-tp),'profileFactorsStripped':'A(1/2)A(0), both nonzero removable/positive factors; no convolution value'})
    J.zero('new-control-qi-dispersion',qip**2,rad(kp));J.zero('new-control-qo-dispersion',qop**2,rad(lp));J.zero('new-control-qh-dispersion',qhp**2,rad(kp+tp));J.zero('new-control-qs-dispersion',qsp**2,rad(lp-tp))
    for face in ('plus','minus'):
        ad=select_control_address(addresses,'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL',face,'pressure',(0,0),(0,0),lambda v:v['jet']['name']=='e_W')
        src=D(ad['sourceField']);cons=D(ad['consumerField']);require(not src.free_symbols and not cons.free_symbols,'actual direct control source/consumer constants')
        control(face+'-wrong-reflected-root',baseline*src*cons,wrong*src*cons,{'addressId':ad['addressId'],'savedObject':'direct Bc',
            'actualSource':src,'actualConsumer':cons,'correct':'q(l-t)','corrupt':'q(k+t)','notIterationMiddleControl':True,'perCommonNonzeroProfileFactor':True})
    J.emit('new-formal-controls',controls);J.emit('inherited-zero-returns',inherited)
    J.emit('analytic-conclusion',{'status':'SOURCE_JOINED_GLOBAL_WEAK_PRESSURE_CERTIFICATE','testSpace':'S(R) x S(R), complex bilinear',
        'ordinaryGrowthDegreeAtMost':3,'heightAction':'subtracted PV plus contact; Holder 1/2 bound, weighted by Pk^(3/2)',
        'assessment':'Analytic induction, uniform integrability and dominated convergence independently assessed; not machine measure theory',
        'scope':'real omega3 source/consumer; response-only outgoing regularization; strict rest bulk LAB_HELD/RHO4; cs[1,2]; all real internal momenta',
        'wholeDirectOnce':True,'nativeIterationOnce':True,'hiddenProjection':False,'planeWaveScattering':False,'finiteInverse':False,'loss':False,
        'sourceDimensions':'inherited/required, not reconstructed','localSlabPartIncluded':False})
    J.emit('consumed-source-index',{'files':sorted(used),'copies':copies,'oldFunctionsCalled':False})
    return {'executionStatus':'COMPLETED_GLOBAL_WEAK_COMPOSITION_CERTIFICATES','fields':len(fieldcert),'addressEntries':len(route),'sourceJets':len(jets),
        'newFormalControls':len(controls),'priorControlsRestored':len(restored),'responseIntegralsEvaluated':False,'finiteSolves':0,
        'scienceAcceptance':False,'weakClaim':'bounded continuous bilinear Schwartz pressure contribution, conditional on assessed analytic arguments and saved dependencies',
        'scatteringOrLossClaim':False,'productionChanges':False}


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
        exact={'sp':sp,'require':require};exec(compile(definitions(Path(manifest['exactHelperSource']).read_text(),('exact_nonzero_number',)),'unchanged-exact-nonzero','exec'),exact)
        result=J.stage('global-weak-composition',{'manifestSha256':gate['manifestSha256'],'reviewSha256':gate['buildReviewRecordSha256']},lambda:run_science(manifest,J,ns,exact['exact_nonzero_number']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        if (args.out/'saved-copy-index.json').exists():
            for v in json.loads((args.out/'saved-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(v['expected']!=v['actual'] for v in records.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False);save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
