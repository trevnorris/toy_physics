#!/usr/bin/env python3
"""Packet-action adapter and analytic truncation preflight. No response integral."""
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
    require(g['status']=='READY_FOR_ONE_PACKET_PREFLIGHT','gate status')
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
        require(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS PACKET-ACTION PREFLIGHT BUILD','literal build report')
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


def validate_selection(addresses,cells):
    selected=[a for a in addresses if a['jet']['channel']=='e_W']
    local=[c for c in cells if c['row']=='THETA_BALANCE' and c['field']=='e_W']
    require(len(addresses)==2652 and len(selected)==544 and len({a['addressId'] for a in selected})==544,'complete selected addresses')
    require(all(a['row']=='THETA_BALANCE' for a in addresses),'actual row')
    statuses={s:sum(a['status']==s for a in selected) for s in {a['status'] for a in selected}}
    require(statuses=={'FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED':102,'EXACT_ZERO_SOURCE_JET':106,'EXACT_ZERO_CONSUMER':336},'all statuses')
    require(len(cells)==400 and len(local)==16 and {(c['xOrder'],tuple(c['grade'])) for c in local}==set(itertools.product(range(4),G)),'local derivative/grade rectangle')
    triples={(c,r,s) for c,r,s in itertools.product(G,repeat=3) if tuple(c[i]+r[i]+s[i] for i in range(2)) in G}
    require(len(triples)==16 and {(tuple(a['consumerGrade']),tuple(a['responseGrade']),tuple(a['sourceGrade'])) for a in selected}==triples,'separate pressure triples')
    for a in selected:
        require(tuple(a['targetGrade'])==tuple(sum(a[n][i] for n in ('consumerGrade','responseGrade','sourceGrade')) for i in range(2)),'target grade')
    return selected,local,statuses


def rpair(value):
    value=sp.cancel(value);re,im=map(sp.simplify,value.as_real_imag())
    require(re.is_Rational is True and im.is_Rational is True,'exact complex-rational coefficient')
    return re,im


def l1(value):
    re,im=rpair(value)
    return abs(re)+abs(im)


def moments_upper(n,s):
    # Gaussian absolute moments: M0=s sqrt(2pi)<3s, M1=2s^2,
    # Mn=(n-1)s^2 M(n-2). Rational upper bounds, no quadrature.
    result=[3*s,2*s*s]
    for j in range(2,n+1):result.append((j-1)*s*s*result[j-2])
    return result[:n+1]


def poly_majorant(n,y,c=5):
    # |p0|<3 for both declared carriers. Triangle-rule derivative recurrence.
    value=sp.Integer(1)
    for _ in range(n):value=sp.expand(sp.diff(value,y)+(3+(y+c)/64)*value)
    require(all(v.is_nonnegative is True and v.is_Rational is True for v in sp.Poly(value,y).all_coeffs()),'positive rational derivative envelope')
    return value


def full_moment_bound(poly,y):
    pp=sp.Poly(poly,y);n=max(0,int(pp.degree())) if not pp.is_zero else 0
    moments=moments_upper(n,sp.Integer(8))
    return sum(pp.nth(i)*moments[i] for i in range(n+1))


def exponential_moment(n,start,rate):
    # Exact polynomial prefactor of integral_start^infty x^n exp(-rate*x) dx.
    return sum(sp.Rational(math.factorial(n),math.factorial(j))*start**j/rate**(n-j+1) for j in range(n+1))


def weighted_tail(degree,start,rate=5):
    # For integer rate*start, e^(-rate*start) <= 2^(-rate*start).
    require(isinstance(start,int) and start>=0 and isinstance(rate,int) and rate>0,'integer tail plan')
    rate=sp.Integer(rate)
    return 2*sp.Rational(1,2)**int(rate*start)*sum(math.comb(degree,n)*exponential_moment(n,sp.Integer(start),rate) for n in range(degree+1))


def run_science(m,J,ns):
    D=ns['decode'];one=ns['one_symbol'];raw={};copies={}
    for alias,r in m['savedInputs'].items():
        p=Path(r['path']);q=J.out/'saved'/alias;q.parent.mkdir(parents=True,exist_ok=True)
        require(sha(p)==r['sha256'],'saved source '+alias);shutil.copyfile(p,q);require(sha(q)==r['sha256'],'saved copy '+alias)
        raw[alias]=json.loads(q.read_text());copies[alias]={'source':str(p),'path':str(q.relative_to(J.out)),'sha256':r['sha256'],'bytes':q.stat().st_size}
        ns['replace_json'](J.out/'saved-copy-index.json',copies)
    inherited=[]
    def inherit(alias):
        require(raw[alias].get('cancelled')==ZERO,'literal saved zero '+alias)
        inherited.append({'alias':alias,'receipt':copies[alias],'functionCalled':False})
    selected,local,statuses=validate_selection(raw['inventory/THETA_BALANCE-ordered-addresses.json'],raw['local/all-local-cells.json'])
    require(selected==raw['selected/pressure-addresses.json']['selected'] and local==raw['selected/local-cells.json']['selected'],'actual complete originals versus reviewed projections')
    J.emit('selection',{'addressIds':[a['addressId'] for a in selected],'localCells':[(c['xOrder'],c['grade']) for c in local],'statuses':statuses,'allSelectedOperands':'saved/inventory/THETA_BALANCE-ordered-addresses.json'})
    physical=raw['physical-input.json'];context=D(raw['local/context.json']);params={k:sp.Rational(v) for k,v in physical['parameters'].items()}
    J.emit('physical-source-input',{'context':context,'physicalInput':physical})
    require(context['physical']==physical and context['saved']['physicalInput']==physical and context['frequencyOverride']=={'old':'1','actual':3},'actual same physical source and held-frequency override')
    require(params['W_0']==1 and params['L_W']==10 and params['rho_m']==sp.Rational(1,10) and params['Lambda_A_0']==sp.Rational(1,100) and params['tau_A']==sp.Rational(1,10),'material and profile constants')
    require(params['s11cdTangentialMomentum1']==sp.Rational(1,5) and params['s11cdTangentialMomentum2']==sp.Rational(1,10),'actual saved tangent bindings')
    require(context['numeric']['omega']==3 and context['effectiveSpeedOnlyInPressure'] is True,'actual real held frequency and speed scope')
    scale=raw['inventory/native-profile-scale-join.json'];require(scale['savedLength']['srepr']=='Integer(10)' and scale['physicalLength']=='10' and scale['declaredLength']==10,'native physical length')
    speed=raw['inventory/prebinding-native-speed-inventory.json'];require(speed['noNativeSpeedSymbols'] is True and all(not v['speedSymbols'] for v in speed['inventory']),'unbound source speed independence')
    origin=D(raw['ends/left-source-binding.json'])['origin'];pairs=origin['mappingPairs']
    require({str(k):v for k,v in pairs}=={'eta_bg':sp.Rational(1,100),'sigma_W':sp.Rational(1,1000)},'actual optional origin')
    cs=sp.sqrt(6)/2;kappa=sp.sqrt(595)/10
    match=D(raw['ends/left-match.json']);point=match['point'];point_values={str(k):v for k,v in point['mappingPairs']} if isinstance(point,dict) else {str(k):v for k,v in point}
    require(point_values['weak_end_cs']==cs and point_values['weak_end_p']==kappa and point_values['weak_end_q']==0,'saved selected matching operands')
    J.zero('matching-dispersion',9/cs**2-sp.Rational(1,20),kappa**2)
    J.emit('physical-plan',{'originalPhysicalInput':'saved/physical-input.json','frequency':3,'cs':cs,'kappa':kappa,'carriers':[kappa,sp.S.Zero],'centers':[-sp.Rational(5,2),sp.Rational(5,2)],'s':8,'origin':origin,'noSourceRebind':True})
    fields=raw['pressure/fields.json'];certs=raw['pressure/coefficient-certificates.json'];require(set(fields)==set(certs) and len(fields)==34,'complete saved fields')
    T=sp.Symbol('weak_tanh_variable',real=True);x=sp.Symbol('composition_x',real=True);norms={};vectors={}
    for fid,field in fields.items():
        poly=D(raw['field/'+fid+'-polynomial.json']);proof=raw['field/'+fid+'-reconstruction-input.json'];inherit('field/'+fid+'-reconstruction-return.json')
        require(proof['left']==field['field'],'original field proof operand')
        quotient=D(certs[fid]['polynomial']);coeff=sp.Poly(quotient,T).all_coeffs();original_coeff=[v/poly['denominator'] for v in poly['coefficients']]
        J.emit('field-'+fid+'-adapter-input',{'original':field,'certificate':certs[fid],'savedPolynomial':poly,'inheritedProofInput':proof,'coefficientVector':coeff})
        require(len(coeff)==len(original_coeff) and quotient.free_symbols <= {T},'coefficient variables/count')
        for i,(v,w) in enumerate(zip(coeff,original_coeff)):J.zero('field-'+fid+'-coefficient-'+str(i),v,w)
        require(D(proof['right'])==quotient.subs(T,sp.tanh(x/10)),'actual right operand of inherited proof')
        vectors[fid]=coeff;norms[fid]=sum(l1(v) for v in coeff)
        J.emit('field-'+fid+'-strip-bound',{'coefficients':coeff,'sumAbsRealImag':norms[fid],'strip':[-5,5],'argument':'|tanh((x+iy)/10)|<=1 for |y|<=5<10*pi/4; modulus <= coefficient real/imag L1 sum','fieldRecalculated':False})
    coverage=raw['weak/weak-address-coverage.json'];byid={a['addressId']:a for a in coverage};require(len(coverage)==13260 and len(byid)==13260,'inherited full address routes')
    route=[];proofs=set();live=[]
    for a in selected:
        J.emit('address-'+str(a['addressId'])+'-adapter-input',{'address':a,'savedRoute':byid[a['addressId']]})
        for role in ('source','consumer'):
            fid=a[role+'Transform']['coefficientId'];require(fields[fid]['field']==a[role+'Field'],'actual field ID')
        old=byid[a['addressId']]
        require(old['face']==a['face'] and old['slot']==a['slot'] and old['component']==a['component'] and old['status']==a['status'] and old['sourceJet']==a['jet'] and old['gradeTriple']==[a[k] for k in ('consumerGrade','responseGrade','sourceGrade')],'exact inherited route arguments')
        require(old['sourceFieldId']==a['sourceTransform']['coefficientId'] and old['consumerFieldId']==a['consumerTransform']['coefficientId'],'inherited field map')
        name=a['fullFactorProof']['proof'];op=raw['factors/'+name+'-operands.json']
        # operandSha256 is the old producer's tuple digest, NOT a file digest.
        # Byte receipts come from savedInputs/copy index; the operands join below.
        for key in ('sha256','original','mapped','map','symbolAssumptions','flatSupport','frequency','positiveRegulatorContinuation'):
            require(op['actualResponseMap'][key]==a['responseMap'][key],'factor response map '+key)
        require(op['addressNormalOriginal']==a['normalOriginal'] and op['requiredMap']==a['fullFactorProof']['completeNormalMap'],'normal complete map')
        if name not in proofs:
            full=raw['factors/'+name+'-full-mapped-residual-input.json'];normal=raw['factors/'+name+'-normal-source-join-input.json']
            J.emit(name+'-inherited-inputs',{'operands':op,'full':full,'normal':normal})
            require(full['right']==op['mappedAddressFactor'] and normal['left']==op['addressNormalOriginal'] and normal['right']==op['savedNormal'],'original proof arguments')
            inherit('factors/'+name+'-full-mapped-residual-return.json');inherit('factors/'+name+'-normal-source-join-return.json');proofs.add(name)
        flat=a['component']=='NATIVE_FLAT'
        require(a['responseInputDepth']==('q(l)' if flat else 'q(k)') and a['responseOutputDepth']=='q(l)' and a['sourceTransform']['transfer']==('l-p' if flat else 'k-p') and a['consumerTransform']['transfer']=='r-l','actual depth/Fourier support')
        require(a['responseMap']['frequency']==3 and a['responseMap']['positiveRegulatorContinuation'] is False,'real source family')
        jet=a['jet'];nt=jet['timeOrder'];n1,n2,n3=jet['spatialOrders'];p=sp.Symbol('composition_p',real=True)
        J.zero('address-'+str(a['addressId'])+'-wave-jet',D(a['waveMultiplier']),(-3*sp.I)**nt*(sp.I*p)**n1*(sp.I/5)**n2*(sp.I/10)**n3)
        if a['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED':
            require(a['epsilonCount']==1 and norms[a['sourceTransform']['coefficientId']]>0 and norms[a['consumerTransform']['coefficientId']]>0,'live addressed product')
            live.append(a)
        else:
            require(a['epsilonCount']==0,'saved zero epsilon convention')
            fid=a[('consumer' if a['status']=='EXACT_ZERO_CONSUMER' else 'source')+'Transform']['coefficientId'];require(norms[fid]==0,'actual exact zero field')
        route.append({'addressId':a['addressId'],'status':a['status'],'proof':name,'sourceFieldId':a['sourceTransform']['coefficientId'],'consumerFieldId':a['consumerTransform']['coefficientId'],'jet':jet,'tagCheckInherited':old['wholeTagCheck']})
    J.emit('address-adapter-result',route)
    local_bounds=[]
    for c in local:
        n=c['xOrder'];g=tuple(c['grade']);name='local-'+str(n)+'-'+''.join(map(str,g));v=D(c)
        J.emit(name+'-input',c)
        localT=sp.Symbol('full_weak_tanh_variable',real=True);localx=sp.Symbol('full_weak_x',real=True)
        require([z['name'] for z in v['identities']]==['cell-source-sum','cell-physical-polynomial','cell-derivative'] and
            v['identities'][0]['left']==v['coefficient'] and
            v['identities'][1]['left']==v['coefficient'].subs(localT,sp.tanh(localx/10)) and
            v['identities'][1]['right']==v['polynomial']['polynomial'].subs(localT,sp.tanh(localx/10)) and
            all(z['cancelled']==0 for z in v['identities']),'inherited local proof operands and returns')
        require(v['polynomial']['original']==v['coefficient'],'local original coefficient')
        for r in v['polynomial']['coefficients']:
            J.zero(name+'-numeric-coefficient-'+str(r['order']),sp.Poly(v['polynomial']['polynomial'],localT).nth(r['order']),r['value'])
        vals=[r['value'] for r in v['polynomial']['coefficients']];bound=sum(l1(a) for a in vals);local_bounds.append({'xOrder':n,'grade':list(g),'coefficients':v['polynomial']['coefficients'],'stripL1':bound,'coefficient':v['coefficient']})
    J.emit('local-adapter-result',local_bounds)
    unit=D(raw['inventory/fourier-and-unit-provenance.json']);require(unit['sourceUnit']==[1,-1,0] and unit['oneDimensionalForwardNormalization']=='1/(2*pi)','saved source/Fourier units')
    units=[v for v in unit['consumerUnitJoins'] if v['row']=='THETA_BALANCE'];require(len(units)==4 and all(v['total']==v['expected'] for v in units),'all four native consumer unit joins')
    J.emit('inherited-unit-contract',{'actualConsumerJoins':units,'sourceUnit':unit['sourceUnit'],'forwardFactor':unit['oneDimensionalForwardNormalization'],'dualTestConvention':'fixed dual THETA row amplitude, complex bilinear; no power','boundNumbersUse':'original reference-unit numeric coordinates','noNewPostBindingDimensionProof':True})
    # Missing numerical adapter identities; the old source/response functions are never called.
    omega=sp.Integer(3);mu=omega*params['rho_m'];aa=params['Lambda_A_0']/(params['rho_m']**2*(1-sp.I*omega*params['tau_A']));beta=sp.cancel(aa*mu)
    k,l,t=sp.symbols('packet_k packet_l packet_t',real=True);qi,qm,qh,qs,qo=sp.symbols('packet_qi packet_qm packet_qh packet_qs packet_qo')
    tags=D(raw['pressure/whole-tags.json']);A=lambda z:10*z/(4*sp.sinh(5*sp.pi*z));j=lambda z:5*A(z)
    for tag,alias in m['wholeDefinitionInputs'].items():
        require(raw['pressure/whole-tags.json'][tag]['savedDefinition']==raw[alias],'actual saved whole definition '+tag)
        require(sha(m['savedInputs'][alias]['path'])==raw['pressure/whole-tags.json'][tag]['sha256'],'whole definition byte receipt '+tag)
    # Complete factors at each original scalar address, including both normal signs.
    # Whole tags are independent formal values; this never evaluates a convolution.
    ck,cl,ccs=sp.Symbol('composition_k',real=True),sp.Symbol('composition_l',real=True),sp.Symbol('composition_cs',positive=True)
    hvalue,jvalue,Hvalue,Jvalue,Dvalue=sp.symbols('packet_h packet_j packet_H packet_J packet_D')
    adapters={};adapted=[]
    for a in selected:
        key=(a['face'],a['slot'],a['component']);label='-'.join(key)
        actual=D(a['responseCoefficient'])*D(a['normalMultiplier'])
        op=raw['factors/'+a['fullFactorProof']['proof']+'-operands.json']
        require(actual==D(op['mappedAddressFactor']),'live complete factor joins saved proof operand')
        if label not in adapters:
            suffix=a['face'];qfun=sp.Function('common_outgoing_q')
            expected_calls={qfun(ck):qi,qfun(cl):qo,
                sp.Function('reference_height_hat')(cl-ck):hvalue,
                sp.Function('reference_slope_hat')(cl-ck):jvalue,
                sp.Function('Hwhole')(cl-ck,1,10):Hvalue,
                sp.Function('Jwhole_'+suffix)(cl,ck,3,ccs,sp.Rational(1,5),sp.Rational(1,10),1,10):Jvalue,
                sp.Function('Dwhole_'+suffix)(cl,ck,3,ccs,sp.Rational(1,5),sp.Rational(1,10),1,10):Dvalue}
            calls=actual.atoms(sp.Function)
            J.emit('numeric-factor-'+label+'-arguments',{'address':a,'actualCompleteFactor':actual,'actualCalls':sorted(calls,key=str),'permittedCalls':[[u,v] for u,v in expected_calls.items()]})
            require(calls <= set(expected_calls),'exact whole/depth/profile call signatures')
            mapped=actual.xreplace({v:expected_calls[v] for v in calls}).xreplace({ck:k,cl:l,ccs:cs})
            comp=a['component']
            expected={'NATIVE_FLAT':mu/(qo+beta),'NATIVE_HEIGHT':-sp.I*mu*qi*hvalue/(qo+beta),
                'NATIVE_SLOPE':mu*k*jvalue/((qo+beta)*(qi+beta)),
                'NATIVE_MIXED_ITERATION':-sp.I*mu*k*qo*Hvalue/((qo+beta)*(qi+beta))+Jvalue,
                'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':Dvalue}[comp]
            if a['slot']=='normal':expected*=sp.I*(1 if suffix=='plus' else -1)*qo
            J.zero('numeric-factor-'+label,mapped,expected)
            adapters[label]={'mapped':mapped,'template':expected,'firstAddressId':a['addressId'],'formalWholeValues':True}
        else:require(actual==adapters[label]['original'],'identical complete factor for adapter reuse')
        adapters[label]['original']=actual;adapted.append({'addressId':a['addressId'],'adapter':label})
    require(len(adapters)==20,'all native face/slot/component adapters')
    J.emit('numeric-factor-adapters',{'definitions':adapters,'addressJoins':adapted})
    for label,mapping,expected in [
        ('J',{'reference_k':k,'reference_l':l,'reference_t':t,'reference_qi':qi,'reference_qm':qm,'reference_qo':qo,'reference_unrestricted_frequency':omega},aa*mu**2*10*sp.Rational(1,4)*A(t)*A(l-k-t)*(k+t)*(2*k+t)*qi/(qm*(qo+beta)*(qm+beta)*(qi+beta)*(qm+qi))),
        ('D',{'k':k,'grazing_output':l,'grazing_transfer':t,'grazing_qi':qi,'grazing_qh':qh,'grazing_qs':qs,'grazing_qo':qo,'grazing_unrestricted_frequency':omega},10/(4*sp.I)*A(t)*A(l-k-t)*(-sp.I*mu)/((qi+beta)*(qo+beta))*(k*(2*l-t)/(qs+qo)+k*(t+2*k)*qi/(qh*(qh+qi))+qi**2/qh))]:
        source=tags['Jwhole' if label=='J' else 'Dwhole']['savedDefinition']['density'];mp={s:mapping[s.name] for s in source.free_symbols};require({s.name for s in source.free_symbols}==set(mapping),'complete density variable map')
        J.zero('new-'+label+'-numerical-adapter',source.xreplace(mp),expected)
    hh=tags['H']['savedDefinition'];hm={s:({'reference_k':k,'reference_l':l}[s.name]) for s in hh['contact'].free_symbols}
    J.zero('new-H-contact-adapter',hh['contact'].xreplace(hm),j(l-k)/4)
    hs=hh['ordinarySubtractedIntegrand'];hmap={'reference_k':k,'reference_l':l,'reference_left_height_transfer':t}
    require({s.name for s in hs.free_symbols}==set(hmap),'complete H variable map')
    J.zero('new-H-subtracted-adapter',hs.xreplace({s:hmap[s.name] for s in hs.free_symbols}),A(t)*(j(l-k-t)-j(l-k))/(2*sp.I*t))
    hp=sp.Rational(1,4)*(1+sp.tanh(x/10));jp=sp.Rational(1,4)*(1-sp.tanh(x/10)**2)
    J.zero('new-physical-profile-derivative',sp.diff(hp,x),jp/10)
    J.emit('Fourier-profile-lemma',{'hPhysical':hp,'jPhysical':jp,'hTransformContact':sp.Rational(1,4),'hPV':A(t)/(2*sp.I*t),'jTransform':j(t),'productToConvolutionFactor':1,'assessedAnalyticIdentities':['Fourier derivative rule i*t*hat(h)=hat(h prime)','transform sech^2(x/L) = L^2*t/(2*sinh(pi*L*t/2)) under 1/(2pi)','h even constant part is delta/4; odd tanh part is symmetric PV'],'noIntegralEvaluated':True,'zeroAInherited':'pressure/profile-envelope.json','distributionLemmaNotCASProof':True})
    # Coarse rational envelopes deliberately avoid floating comparisons of a claimed bound.
    y=sp.Symbol('packet_abs_y',nonnegative=True);bounds={}
    for a in live:
        jet=a['jet'];nt=jet['timeOrder'];n1,n2,n3=jet['spatialOrders'];native=sp.Integer(3)**nt/sp.Integer(5)**n2/sp.Integer(10)**n3
        sx=norms[a['sourceTransform']['coefficientId']];cy=norms[a['consumerTransform']['coefficientId']]
        px=poly_majorant(n1,y);py=sp.Integer(1);py1=y+sp.Rational(5,2)+5
        CX=sp.cancel(native*sx*full_moment_bound(px,y)/3);CY=2*cy*full_moment_bound(py,y);CY1=2*cy*full_moment_bound(py1,y)
        require(all(v.is_Rational is True and v>=0 for v in (CX,CY,CY1)),'rational Fourier envelope constants')
        bounds[a['addressId']]={'X':CX,'Y':CY,'Yprime':CY1,'sourcePolynomialMajorant':px,'testDerivativeMajorant':py1,'nativeJetMagnitude':native}
    J.emit('Fourier-envelope-constants',{'bounds':{str(k):v for k,v in bounds.items()},'proof':['|p0|<=3 and |Im z|<=5; Gaussian derivative triangle recurrence','|polynomial(tanh)| <= coefficient L1(real,imag)','exp(25/128)<2; 1/(2*pi)<1/6; sqrt(2*pi)<3','X bound CX exp(-5|k-p0|), Y/Yprime analogous around p0; Y includes 2pi','all carrier/center constants are included; no transform evaluated'],'carrierUniformUpper':3})
    b=sp.Rational(3000,11101);Cq=sp.Integer(36) # a_*=sqrt(879)/20 >1, hence saved 36/sqrt(a_*) <36.
    domain=D(raw['pressure/global-parameter-domain.json']);whole=D(raw['pressure/whole-envelopes.json']);profile=D(raw['pressure/profile-envelope.json'])
    shift=D(raw['pressure/shift-root-bound.json']);height=D(raw['pressure/plus-height-PV.json']);height_minus=D(raw['pressure/minus-height-PV.json'])
    J.emit('saved-global-bound-inputs',{'domain':domain,'whole':whole,'profile':profile,'shift':shift,'heightPlus':height,'heightMinus':height_minus,'actualBeta':beta})
    require(domain['betaMinimum']==b and domain['kappaMinimum']==sp.sqrt(879)/20 and domain['kappaMaximum']==3 and
        domain['muMagnitudeUpper']==sp.Rational(2,5) and domain['aMagnitudeUpper']==1 and domain['betaMagnitudeUpper']==sp.Rational(2,5),'actual global domain constants')
    J.zero('new-beta-real-specialization',domain['beta'].subs(one([domain['beta']],'weak_Omega'),3),beta)
    require(sp.re(beta)>b and kappa**2<9 and domain['kappaMinimum']**2>1,'rational exact bound comparisons')
    require(profile['L']==10 and profile['removableValue']==1/(2*sp.pi) and profile['globalBound']=='11 exp(-|t|)' and profile['productBound']=='121 exp(-|t|)','actual profile envelope')
    J.zero('new-profile-envelope-adapter',profile['A'].subs(one([profile['A']],'weak_t'),t),A(t))
    cq_saved=36/sp.sqrt(domain['kappaMinimum'])
    J.zero('saved-shift-bound-constant-adapter',shift['inverseDepthBound'],cq_saved)
    require(shift['perRootBound']==18 and shift['routeCenters']==['-k +/- kappa_delta','l +/- kappa_delta'],'both actual singular routes')
    J.zero('saved-J-bound-constant-adapter',whole['constants']['J'],sp.Rational(4,5)*121*cq_saved/b**3)
    J.zero('saved-D-bound-constant-adapter',whole['constants']['D'],18*121*2*cq_saved/b**2)
    require(whole['normalExtraConstant']==4 and whole['Hbound']==100 and whole['pressureWeight']=='(1+|k|+|l|)^2' and whole['normalWeight']=='(1+|k|+|l|)^3','actual whole bound weights')
    KJ=sp.Rational(4,5)*121*Cq/b**3;KD=18*121*2*Cq/b**2;Cordinary=1+4*KJ+4*KD+400/b**2
    A0=sp.Rational(8,5)/b;A1=sp.Rational(16,5)/b**2 # sqrt(3)<2
    for face,h in [('plus',height),('minus',height_minus)]:
        J.zero('saved-'+face+'-height-bound-adapter',h['pressureSupConstant'],A0)
        J.zero('saved-'+face+'-height-Holder-adapter',h['pressureHolderConstant'],8*sp.sqrt(3)/(5*b**2))
    require(all(a['slot']=='pressure' for a in live if a['component']=='NATIVE_HEIGHT'),'height bound uses actual pressure-only live addresses')
    F={n:weighted_tail(n,0) for n in range(4)};E15=sp.Integer(3)**15;E30=E15**2
    J.emit('tail-bound-derivation',{'b':b,'CqUpper':Cq,'KJUpper':KJ,'KDUpper':KD,'ordinaryKernelUpper':Cordinary,'heightMagnitude':A0,'heightHolder':A1,'fullExponentialMoments':{str(k):v for k,v in F.items()},'inheritedBounds':D(raw['pressure/whole-envelopes.json']),'proof':['a_*>1, sqrt(3)<2, pi>3, 2<e<3 give rational overestimates','P^3 <= (1+|k|)^3(1+|l|)^3; weighted exponential tails integrated analytically','all-real ordinary outer union bound =2 C CX CY e^30 F3 tail3(K); overlap overcount is safe','height inner bound <=[(A0/4+11 A0)CY+(2 A1 CY+A0 CYprime)/6]*(1+|k|)^2','height k tail uses X only; separate Q>U tail <=11 A0 CX CY e^15 F1 exp(-U), U>=1','T>=K+4, kappa<3 implies all middle/reflected depths have modulus>1 outside |t|<=T','middle J/D tails use polynomial exponential tails of degree2/1, not whole existence constants','H middle tail <=(55/3) exp(-T), coefficient/normal upper (4/b) P^2'],'allAnalyticNotQuadrature':True})
    def contributions(K,U,Tlim):
        data=[]
        for a in live:
            z=bounds[a['addressId']];cx,cy,dy=z['X'],z['Y'],z['Yprime'];comp=a['component'];outer=middle=height=sp.S.Zero
            if comp=='NATIVE_HEIGHT':
                ch=(A0/4+11*A0)*cy+(2*A1*cy+A0*dy)/6
                outer=ch*cx*E15*weighted_tail(2,K);height=11*A0*cx*cy*E15*F[1]*sp.Rational(1,2)**U
            elif comp=='NATIVE_FLAT':
                outer=Cordinary*cx*cy*E30*weighted_tail(3,K) # drop one exponential; conservative diagonal bound
            else:outer=2*Cordinary*cx*cy*E30*F[3]*weighted_tail(3,K)
            if comp=='NATIVE_MIXED_ITERATION':
                middle=4*cx*cy*E30*F[3]**2*(sp.Rational(4,5)*121/b**3)*weighted_tail(2,Tlim,1)
                middle+=(4/b)*cx*cy*E30*F[2]**2*sp.Rational(55,3)*sp.Rational(1,2)**Tlim
            elif comp=='INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':
                middle=4*cx*cy*E30*F[3]**2*(36*121/b**2)*weighted_tail(1,Tlim,1)
            data.append({'addressId':a['addressId'],'component':comp,'grade':a['targetGrade'],'outer':outer,'heightQ':height,'middle':middle,'total':outer+height+middle})
        return data
    K,U,Tlim=4,1,8;budget=sp.Rational(1,10**11);iterations=[]
    while True:
        parts=contributions(K,U,Tlim);totals={k:sum(v[k] for v in parts) for k in ('outer','heightQ','middle')};iterations.append({'K':K,'U':U,'T':Tlim,**totals})
        J.emit('tail-selection-step-'+str(len(iterations)),iterations[-1])
        if all(v<budget/3 for v in totals.values()):break
        if totals['outer']>=budget/3:K+=1
        if totals['heightQ']>=budget/3:U+=1
        if totals['middle']>=budget/3:Tlim+=1
        Tlim=max(Tlim,K+4)
        require(K<=m['preflightCapacity']['maximumOuterRadius'] and U<=m['preflightCapacity']['maximumHeightRadius'] and Tlim<=m['preflightCapacity']['maximumMiddleRadius'],'declared preflight capacity, no quadrature or enlarged scope')
    J.emit('tail-selection-prefix',iterations)
    J.emit('tail-plan',{'K':K,'U':U,'T':Tlim,'allAddresses':parts,'sum':sum(v['total'] for v in parts),'budget':budget,'coversBothCarriers':True,'componentGradeSumsBelowBudgetByPositiveTotal':True,'noCancellationUsed':True,'domain':'ordinary square |k|,|l|<=K; height |k|<=K and 0<Q<=U; H/J/D middle |t|<=T','quadratureReady':False,'remaining':'actual Fourier tails/radii, collision panels, two complete numerical routes and controls have not run'})
    for face in ('plus','minus'):
        require(any(a['face']==face and a['component']=='INHERITED_DIRECT_WHOLE_OFF_DIAGONAL' for a in live),'both face live direct inventory')
    controls={}
    for name,predicate in [('Hcontact',lambda a:a['component']=='NATIVE_MIXED_ITERATION'),('reflectedRoot',lambda a:a['component']=='INHERITED_DIRECT_WHOLE_OFF_DIAGONAL'),('normalDepth',lambda a:a['slot']=='normal' and a['component']=='NATIVE_SLOPE' and not fields[a['consumerTransform']['coefficientId']]['constant']),('Leibniz',lambda a:a['jet']['spatialOrders'][0]>0 and not fields[a['sourceTransform']['coefficientId']]['constant'])]:
        candidates=[a for a in live if predicate(a)];controls[name]={'eligibleAddressIds':[a['addressId'] for a in candidates],'selected':min((a['addressId'] for a in candidates),default=None),'numericalResponse':'NOT_EVALUATED'}
        require(bool(candidates),'predeclared live metadata control '+name)
    J.emit('control-address-preselection',controls);J.emit('inherited-zero-returns',inherited)
    return {'status':'PACKET_PREFLIGHT_CERTIFICATES_COMPLETE_NO_ACTION','addresses':len(selected),'liveAddresses':len(live),'localCells':len(local),'fields':len(fields),'inheritedZeroReturns':len(inherited),'K':K,'U':U,'T':Tlim,'tailBudget':str(budget),'quadratureRun':False,'numericalAction':None,'futureEvaluatorBuildRequired':True,'scientificAcceptance':False}


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
        result=J.stage('packet-preflight',{'manifestSha256':g['manifestSha256'],'buildReviewSha256':g['buildReviewRecordSha256']},lambda:run_science(m,J,ns));code=0
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
