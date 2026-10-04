#!/usr/bin/env python3
"""Guarded exact finite-window packet contraction certificate.

New algebra/restoration only after containment; no old scientific function or
integral is called. Source/JSON inspection and synthetic tests are separate.
"""
import argparse,ast,hashlib,importlib.util,json,math,os,resource,shutil,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
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
    require(g['status']=='READY_FOR_ONE_PACKET_CONTRACTION' and g['independentBuildClearance'] is True,'actual independent build readiness')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==m['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','library','authority','buildReviewRecord'):
        require(sha(g[key])==g[key+'Sha256'],'gate pin '+key)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual guard/supervisor')
    require(g['launcher']==m['launcher'] and g['library']==m['librarySource'] and g['authority']==m['executionAuthority'] and g['buildReviewRecord']==m['reviewRecordWillBe'],'actual document paths')
    r=read(g['buildReviewRecord']);require(r['allChecksPassed'] is True and r['independentBuildClearance'] is True,'build assessment')
    for key in ('workerSha256','manifestSha256','librarySha256','launcherSha256','sharedGuardSha256','supervisorSha256'):require(r[key]==g[key],'actual reviewed '+key)
    require(all(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED FINITE-WINDOW CONTRACTION BUILD' for e in ('claude','grok')),'both literal build reports')
    method=read(m['methodRecord']);require(g['methodRecordSha256']==sha(m['methodRecord']) and method['jointIndependentMethodClearance'] and method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'actual pressure-readiness method')
    a=read(g['authority']);require(a['scope']==g['scope']==m['scope'] and a['boundedInstrumentAuthorized'] and a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and a['automaticScientificRetry'] is False and a['noDeadline'] and g['durationLimits'] is None,'standing bounded authority')
    return g


def load_library(m):
    spec=importlib.util.spec_from_file_location('packet_source_unit_transport',m['librarySource']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def run(m,J,C):
    raw,copies=copy_inputs(m,J)
    def get(alias):return raw[alias]
    def inherited(prefix):
        args,ret=get(prefix+'-input.json'),get(prefix+'-return.json')
        J.emit('inherited-'+prefix.replace('/','--'),{'input':args,'return':ret,'functionCalled':False})
        require(type(ret.get('cancelled')) is dict and ret['cancelled']==ZERO,'published exact cancellation zero '+prefix)
        return args
    def exact(name,left,right,context=None):
        J.start(name,{'left':left,'right':right,'context':context,'newDerivation':True})
        residual=left-right;J.emit(name+'-raw',{'residual':residual})
        reduced=sp.cancel(sp.together(residual));J.emit(name+'-decision',{'cancelledResidual':reduced,'literalZero':reduced is sp.S.Zero})
        require(reduced is sp.S.Zero,'new exact identity '+name);J.finish({'residual':reduced});return reduced
    contracts=get('source-contracts.json')
    for rec in contracts['fragments'].values():
        C.fragment(rec)
        historical=get('accepted-units/posthashes.json')['sources'][rec['source']]
        require(historical['intact'] is True and historical['expected']==historical['actual']==rec['sourceSha256'],'actual source used by accepted unit bridge')
    J.emit('original-source-contracts',contracts)
    accepted=get('accepted-units/result-record.json')
    require(accepted['allChecksPassed'] is True and accepted['status']=='BOUNDED_KERNEL_WAVE_MEASURE_UNITS_ACCEPTED_OUTER_EVALUATOR_PENDING','accepted complete summand unit result')
    for alias,r in m['savedInputs'].items():
        if alias.startswith('accepted-units/') and alias!='accepted-units/result-record.json':
            previous=accepted['records']['complete/'+alias.removeprefix('accepted-units/')]
            require((r['sha256'],r['bytes'])==(previous['sha256'],previous['bytes']),'accepted operand bytes '+alias)
    result=get('accepted-units/journal-result.json')
    require(result['pressureSummandUnitsComplete'] and result['selectedAddresses']==544 and result['templates']==20,'completed unit scope')
    physical=get('preflight/physical-plan.json');origins=get('accepted-units/parameter-quantity-origins.json')
    require(physical==get('saved/preflight/physical-plan.json')==origins['physicalPlan'],'same saved physical plan')
    require(origins['context']==get('saved/native/binding.json')==get('saved/source/transport-context.json')['actual']['restored'],'actual physical/source binding')
    require(origins['context']['physicalInput']==get('saved/physical-input.json'),'same original physical input')
    require(physical['frequency']==3 and physical['cs']['srepr']=='Mul(Rational(1, 2), Pow(Integer(6), Rational(1, 2)))' and physical['kappa']['srepr']=='Mul(Rational(1, 10), Pow(Integer(595), Rational(1, 2)))','fixed real frequency/speed/edge depth')
    require(physical['s']==8 and physical['centers']==[{'text':'-5/2','srepr':'Rational(-5, 2)'},{'text':'5/2','srepr':'Rational(5, 2)'}],'same two packets')
    for name in ('runtime-a','runtime-mu','new-runtime-J-arithmetic','new-runtime-D-arithmetic-sum'):inherited('saved/inner/'+name)
    J.emit('inherited-physical-parameters',{'origins':origins,'actualProofs':{n:get('saved/inner/'+n+'-input.json') for n in ('runtime-a','runtime-mu')},'parameterFunctionsCalled':False})
    aa=C.restore_scalar(sp,get('saved/inner/runtime-a-input.json')['left']);mu=C.restore_scalar(sp,get('saved/inner/runtime-mu-input.json')['left'])
    beta0=aa*mu;kap=C.restore_scalar(sp,physical['kappa'])
    width=C.restore_scalar(sp,origins['context']['numeric']['W_0']);length=C.restore_scalar(sp,origins['context']['numeric']['L_W'])
    require(width is sp.S.One and length==sp.Integer(10),'actual saved physical W and L')
    J.emit('new-physical-domain-decisions',{'beta':beta0,'kappa':kap,'W':width,'L':length,'betaRealPositive':sp.re(beta0).is_positive,'betaImagPositive':sp.im(beta0).is_positive,'kappaPositive':kap.is_positive,'kappaSquareBetweenZeroAndNine':bool(0<sp.Rational(119,20)<9)})
    require(sp.re(beta0).is_positive is True and sp.im(beta0).is_positive is True and kap.is_positive is True,'outgoing physical domain')
    exact('new-bound-beta-binding',beta0,sp.Rational(30,109)+sp.I*sp.Rational(9,109),{'restored_a':aa,'restored_mu':mu})
    exact('new-bound-kappa-square',kap**2,sp.Rational(119,20),physical)
    fields=get('saved/pressure/fields.json');require(len(fields)==34,'all saved coefficient fields')
    for key,field in fields.items():
        proof=inherited('saved/field/'+key+'-reconstruction')
        J.emit('restored-field-'+key,{'field':field,'polynomial':get('saved/field/'+key+'-polynomial.json'),'proofArguments':proof})
        require(proof['left']==field['field'],'original polynomial/field left operand')
    factors={}
    for i in range(17):
        name='address-full-factor-'+str(i);fac=get('saved/factors/'+name+'-operands.json')
        full=inherited('saved/factors/'+name+'-full-mapped-residual');normal=inherited('saved/factors/'+name+'-normal-source-join')
        J.emit('restored-factor-'+str(i),{'actualOperands':fac,'fullProof':full,'normalProof':normal});factors[name]=fac
    # Bind the actual runtime call and its four distinct depth arguments as
    # source ASTs. This inspects definitions, and never invokes the old methods.
    middle=C.fragment(contracts['fragments']['S11c_d_defect_packet_inner_lib.py:middle'])
    assignments={ast.unparse(n.targets[0]):n.value for n in middle.body if isinstance(n,ast.Assign) and len(n.targets)==1}
    expected={'(qi, qo)':'(self.q(c, k), self.q(c, l))','qh':'self.q(c, k + t)','qs':'self.q(c, l - t)','used_qs':'qh if mutate else qs','a':'1 / (1 - 3 * c.j / 10)','mu':'c.mpf(3) / 10','beta':'a * mu','(aa, bb)':'(self.profile(c, t), self.profile(c, l - k - t))','values':'kernel_components(k, l, t, qi, qo, qh, used_qs, aa, bb, a, mu, beta, c.one, c.mpf(10), c.j)'}
    observed={key:ast.unparse(assignments[key]) for key in expected}
    J.emit('original-middle-call-bindings',{'source':contracts['fragments']['S11c_d_defect_packet_inner_lib.py:middle'],'observed':observed,'required':expected,'called':False})
    require(all(ast.dump(assignments[key])==ast.dump(ast.parse(value,mode='eval').body) for key,value in expected.items()),'actual depth/profile/parameter argument routing')
    rootfn=C.fragment(contracts['fragments']['S11c_d_defect_packet_inner_lib.py:q'])
    root_assign=[n for n in rootfn.body if isinstance(n,ast.Assign)];root_return=[n for n in rootfn.body if isinstance(n,ast.Return)]
    require(len(root_assign)==len(root_return)==1,'actual branch definition')
    require(ast.dump(root_assign[0].value)==ast.dump(ast.parse('c.mpf(595)/100-p*p',mode='eval').body) and ast.dump(root_return[0].value)==ast.dump(ast.parse('c.sqrt(d) if d>=0 else c.j*c.sqrt(-d)',mode='eval').body),'same positive outgoing branch and radicand')
    J.emit('original-outgoing-branch',{'source':contracts['fragments']['S11c_d_defect_packet_inner_lib.py:q'],'physicalKappa':physical['kappa'],'actualRadicand':ast.unparse(root_assign[0].value),'actualBranch':ast.unparse(root_return[0].value),'called':False})
    for name in ('H','Jwhole','Dwhole'):
        tag=get('saved/pressure/whole-tags.json')[name];definition=get('saved/whole-origin/'+name+'.json')
        J.emit('inherited-complete-whole-'+name,{'tag':tag,'definition':definition,'called':False})
        require(tag['savedDefinition']==definition,'actual complete whole signature '+name)
    J.emit('inherited-profile-origin',{'source':get('accepted-unit-source-contracts.json')['fragments']['reference-A'],'nativeLemma':get('saved/inner/physical-H-and-PV-lemma.json'),'oldNumericalProfileSource':contracts['fragments']['S11c_d_defect_packet_inner_lib.py:profile'],'scope':'New continuous algebra uses the exact analytic A in saved response definition. The old finite sinhc evaluation and its numerical error remain future evaluator obligations; no Fourier/profile numeric call.'})
    # Interpret only the original finite rational assignment AST on fresh symbols.
    names='k l t qi qo qh qs A1 A2 a mu beta W L'.split();sym={n:sp.Symbol('contract_'+n) for n in names};sym['I']=sp.I
    fn=C.fragment(contracts['fragments']['S11c_d_defect_packet_inner_lib.py:kernel_components'])
    J.start('new-original-density-interpretation',{'source':contracts['fragments']['S11c_d_defect_packet_inner_lib.py:kernel_components'],'arguments':sym,'oldFunctionCalled':False})
    arithmetic=C.Arithmetic(sym,sp.Integer)
    try:densities=arithmetic.statements(fn)
    finally:J.emit('new-original-density-intermediates',arithmetic.events)
    require(len(densities)==4,'all four added primitive returns');J.finish({'orderedDensities':densities})
    k,l,t,qi,qo,qh,qs,A1,A2,a,mmu,beta,W,L=[sym[n] for n in names]
    mid,qm,Amk,Alm,X,Y,N=sp.symbols('contract_m contract_qm contract_Amk contract_Alm contract_X contract_Y contract_N')
    x=X/(qi+beta);y=N*Y/(qo+beta);Cj=a*mmu**2*W*L/4;Cd=(-sp.I*mmu)*W*L/(4*sp.I)
    # Independent contractions are symbols only AFTER their complete integrands
    # have been recorded. No Integral, Fourier request or numerical rule is used.
    integrands={'Y0':Alm*y,'X1':Amk*k*x,'X1T':Amk*k*x,'X2T':Amk*k**2*x,'X02T':Amk*qi**2*x,
        'C0T':Amk*qi*x/(qm+qi),'C1T':Amk*k*qi*x/(qm+qi),'C2T':Amk*k**2*qi*x/(qm+qi),
        'Y0CT':Alm*y/(qm+qo),'Y1CT':Alm*l*y/(qm+qo),'Y0C':Alm*y/(qm+qo),'Y1C':Alm*l*y/(qm+qo)}
    variables={z:('l' if z.startswith('Y') else 'k') for z in integrands}
    windows={z:('I_T(m)' if z.endswith('T') else '[-K,K]') for z in integrands}
    z={name:sp.Symbol('contraction_'+name) for name in integrands}
    final={'J':Cj*z['Y0']*mid/(qm*(qm+beta))*(z['C1T']+mid*z['C0T']),
        'Dr':Cd*z['X1']*(z['Y1CT']+mid*z['Y0CT']),
        'Dh':Cd*z['Y0']/qm*(z['C2T']+mid*z['C1T']),
        'Dq':Cd*z['Y0']/qm*z['X02T'],
        'Dr_wrong_root':Cd*(2*z['X1T']*z['Y1C']+(z['X2T']-mid*z['X1T'])*z['Y0C'])}
    J.emit('new-contraction-definitions',{'integrands':integrands,'variables':variables,'windows':windows,'outerIntegrands':final,'x':x,'y':y,'sourceTestIndependent':True,'normalAtOutput':N,'integralsEvaluated':False})
    mapped={};complete={}
    for i,name in enumerate(('J','Dr','Dh','Dq')):
        original=densities[i]*X*N*Y
        # The unused root is explicitly mapped to the other affine argument,
        # never silently identified with q(m).
        qother=sp.Symbol('contract_q_k_plus_l_minus_m')
        mapping={t:l-mid if name=='Dr' else mid-k,A1:Alm if name=='Dr' else Amk,A2:Amk if name=='Dr' else Alm,
                 qs:qm if name=='Dr' else qother,qh:qother if name=='Dr' else qm}
        left=original.xreplace(mapping);right=final[name].xreplace({z[n]:v for n,v in integrands.items()})
        context={'originalDensity':densities[i],'source':X,'test':Y,'normal':N,'mapping':[[s,v] for s,v in mapping.items()],'clippedVariable':'l' if name=='Dr' else 'k','profileEvennessRequired':False,'allDenominatorsNonzeroPointwise':True,'exceptionalGrazingSetsExcluded':True}
        exact('new-full-factorization-'+name,left,right,context);mapped[name]=left;complete[name]=right
        require(qother not in left.free_symbols,'absent opposite root in separate added primitive')
    # A(z) is even. The reflected substitution sends A1=A(l-m), A2=A(m-k)
    # directly, so no extra sign or phase is introduced by reordering products.
    wrong_original=densities[1].xreplace({qs:qh})*X*N*Y
    wrong=wrong_original.xreplace({t:mid-k,A1:Amk,A2:Alm,qh:qm})
    exact('new-wrong-root-mutant-factorization',wrong,final['Dr_wrong_root'].xreplace({z[n]:v for n,v in integrands.items()}),{'actualOriginalMutant':wrong_original,'clippedVariable':'k','outputContractionUnclipped':True,'notNumericalControl':True})
    # All twenty saved complete native templates retain exactly one injected
    # J/D whole value and the actual normal multiplier. Other slots get zero.
    adapters=get('saved/preflight/numeric-factor-adapters.json');require(len(adapters['definitions'])==20,'all original complete templates')
    injections={}
    for label,adapter in adapters['definitions'].items():
        accepted_template=get('accepted-units/complete-template-'+label+'.json')
        args=get('saved/preflight/numeric-factor-'+label+'-arguments.json')
        proof=inherited('saved/preflight/numeric-factor-'+label)
        require(proof['left']==adapter['mapped'] and proof['right']==adapter['template'],'original full template equality operands')
        J.emit('template-context-'+label,{'accepted':accepted_template,'actualArguments':args,'adapter':adapter})
        require(accepted_template['actualAdapter']==adapter and accepted_template['actualArguments']==args,'original complete template/arguments')
        require(accepted_template['wholeDAdditionalMiddleIntegral'] is False and accepted_template['wholeDAdditionalResolvents'] is False,'whole addend once')
        template=C.restore_scalar(sp,adapter['template']);symbols=C.named_symbols(template)
        normal=sp.Integer(1) if '-pressure-' in label else sp.I*accepted_template['normalSign']*symbols['packet_qo']
        require(accepted_template['normalSign'] in (-1,1) if '-normal-' in label else accepted_template['normalSign'] is None,'actual native normal sign')
        delta=sp.Symbol('contract_injected_whole')
        for tag,component in [('packet_J','NATIVE_MIXED_ITERATION'),('packet_D','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')]:
            active=label.endswith(component);old=symbols.get(tag)
            require((old is not None)==active,'actual whole tag occurrence')
            mutant=template if old is None else template.xreplace({old:old+delta})
            exact('new-native-injection-'+label+'-'+tag,mutant-template,normal*delta if active else sp.S.Zero,{'actualTemplate':adapter,'normalSource':accepted_template['sourceNormalInsertion'],'normalDepth':accepted_template['normalDepth'],'noAdditionalResolvents':True})
        injections[label]={'normal':normal,'template':template,'normalUnit':C.MOMENTUM if '-normal-' in label else C.ZERO}
    selection=get('saved/selected/pressure-addresses.json');selected=selection['selected']
    originalrow=get('saved/inventory/THETA_BALANCE-ordered-addresses.json')
    require(len(originalrow)==2652 and selected==[v for v in originalrow if v['row']=='THETA_BALANCE' and v['jet']['channel']=='e_W'],'complete original row selection')
    require(len(selected)==544 and len({v['addressId'] for v in selected})==544,'complete selected address set')
    oldunits=get('accepted-units/complete-summand-dimensions.json');unitby={r['addressId']:r for r in oldunits}
    require(set(unitby)=={a['addressId'] for a in selected},'same full unit address census')
    address_results=[];env_bound=get('preflight/Fourier-envelope-constants.json')
    for index,address in enumerate(selected):
        ident=address['addressId'];label='-'.join((address['face'],address['slot'],address['component']))
        ui=get('accepted-units/summand-'+str(ident)+'-input.json');ur=get('accepted-units/summand-'+str(ident)+'-return.json')
        J.start('address-'+str(ident),{'address':address,'acceptedUnitArguments':ui,'acceptedUnitReturn':ur,'adapter':adapters['definitions'][label],'originalFactorOperands':factors[address['fullFactorProof']['proof']]})
        require(ui['address']==ui['sourceTransportInput']['address']==address and ur==unitby[ident],'full source/grade/row/address and unit join')
        require(ui['adapter']==adapters['definitions'][label] and adapters['addressJoins'][index]=={'addressId':ident,'adapter':label},'original complete native routing')
        require(ui['normal']==address['normalMultiplier'] and ui['waveProof']['left']==address['waveMultiplier'],'actual normal and wave arguments')
        require(ui['sourceTransportReturn']['nativeSourceAndSavedGradeProfileJoined'] is True and ui['sourceTransportReturn']['addressId']==ident,'inherited native source/grade/profile argument')
        require(factors[address['fullFactorProof']['proof']]['mappedAddressFactor']==ui['adapter']['original'],'complete native factor proof operand')
        for role in ('source','consumer'):
            fieldid=address[role+'Transform']['coefficientId'];require(fields[fieldid]['field']==address[role+'Field'],'actual field identity')
        require(address['epsilonCount']==(0 if address['status'].startswith('EXACT_ZERO') else 1),'same epsilon and explicit zero')
        require(address['responseInputDepth']=='q(k)' and address['responseOutputDepth']=='q(l)' or address['responseMap']['flatSupport']=='k=l','actual off-diagonal depth ordering')
        applicable=address['component'] in ('NATIVE_MIXED_ITERATION','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')
        report={'addressId':ident,'allOriginalJoins':True,'explicitZero':address['status'].startswith('EXACT_ZERO'),'contractedOrdinaryComponent':applicable,'otherComponentsUnchanged':not applicable,'nativeNormalAt':'l','actualUnitsInherited':ur,'numericReuse':False}
        if applicable:
            require(address['responseMap']['flatSupport'] is None and address['measure']=='dl dk','ordinary off-diagonal support')
            nunit=injections[label]['normalUnit'];du={s:C.MOMENTUM for s in (k,l,t,qi,qo,qh,qs,mid,qm,beta)}
            du.update({A1:C.ZERO,A2:C.ZERO,Amk:C.ZERO,Alm:C.ZERO,a:C.unit(get('accepted-units/new-a-units-return.json')['unit']),mmu:C.unit(get('accepted-units/new-mu-units-return.json')['unit']),W:C.unit(origins['registry']['W_0']),L:C.unit(origins['registry']['L_W']),X:C.unit(ur['X']),Y:C.unit(ur['Y']),N:nunit})
            traces=[];cu={name:C.add(C.dimension(sp,v,du,traces),C.MOMENTUM) for name,v in integrands.items()}
            fu={**du,**{z[name]:u for name,u in cu.items()}}
            terms=('J',) if address['component']=='NATIVE_MIXED_ITERATION' else ('Dr','Dh','Dq')
            totals={name:C.add(C.dimension(sp,final[name],fu,traces),C.MOMENTUM) for name in terms}
            J.emit('address-'+str(ident)+'-new-unit-walk',{'events':traces,'contractionUnits':cu,'outerTotals':totals,'originalTotal':ur['total'],'integralsNotEvaluated':True})
            require(all(v==C.unit(ur['total']) for v in totals.values()),'new contraction measure units agree with actual complete summand')
            report.update(contractionUnits=cu,totals=totals)
            if not report['explicitZero']:
                bound=env_bound['bounds'][str(ident)];cx=C.restore_scalar(sp,bound['X']);cy=C.restore_scalar(sp,bound['Y'])
                require(cx.is_Rational is True and cy.is_Rational is True and bool(cx>=0) and bool(cy>=0),'restored rational source envelopes')
                majorants=[]
                for radius,tail in ((27,122),(29,124)):
                    maxm=radius+tail;b=sp.Rational(30,109);bx=cx/b;by=cy*(radius+3 if address['slot']=='normal' else 1)/b
                    cjcap=sp.Rational(130,109)*sp.Rational(3,10)**2*width*length/4;cdcap=sp.Rational(3,10)*width*length/4
                    constants={'J':cjcap*maxm*(maxm+radius)*2*bx*by/b,'Dh':cdcap*radius*(maxm+radius)*2*bx*by,'Dq':cdcap*(radius+3)**2*bx*by,'Dr':cdcap*radius*(radius+maxm)*2*bx*by}
                    require(all(v.is_Rational is True and v>=0 for v in constants.values()),'finite rational compact majorants')
                    majorants.append({'K':radius,'T':tail,'M':maxm,'boundForAbsoluteTripleDensity':{name:constants[name] for name in terms},'denominator':'abs(q(m))','additionalIntegrated_k_l_area':(2*radius)**2,'cutoffDependent':True})
                J.emit('address-'+str(ident)+'-compact-majorants',{'originalEnvelope':bound,'originalProof':env_bound['proof'],'actualParameterOrigins':origins['oldRuntimeProofs'],'sourceUnitReturn':ur,'majorants':majorants,'lemma':'quadrant and simple-root certificates below; sqrt(2)<=2; |A|<=1; |q(k)|<=K+3; |a|<=130/109','notNumericalErrorBound':True})
                report['compactMajorants']=majorants
        J.finish(report);address_results.append(report)
    # General affine/inequality identities certify the domains, with witnesses
    # exercising central intervals, both full wings, endpoints and empty cases.
    KK,TT=sp.symbols('contract_K contract_T',positive=True)
    for name,var,old_t in [('input',k,mid-k),('output',l,l-mid)]:
        exact('new-'+name+'-forward-inverse',(var+old_t if name=='input' else var-old_t),mid)
        if name=='input':slacks=[(TT+old_t,var-(mid-TT)),(TT-old_t,(mid+TT)-var)]
        else:slacks=[(TT+old_t,var-(mid-TT)),(TT-old_t,(mid+TT)-var)]
        # Input dt bound signs are exchanged; output signs are as written.
        if name=='input':slacks=[(TT-old_t,var-(mid-TT)),(TT+old_t,(mid+TT)-var)]
        for i,(left,right) in enumerate(slacks):exact('new-'+name+'-interval-slack-'+str(i),left,right,{'assumptions':'T>K>0, real k/l/m','equivalentBounds':['-K<=variable<=K','m-T<=variable<=m+T'],'intersection':'[max(-K,m-T),min(K,m+T)]'})
        J.emit('new-'+name+'-orientation',{'forward':'m=k+t' if name=='input' else 'm=l-t','inverse':'t=m-k' if name=='input' else 't=l-m','signedJacobian':1 if name=='input' else -1,'reverseMappedLimits':name=='output','absoluteJacobian':1,'clippedVariable':str(var),'outer':'[-K-T,K+T]','central':'abs(m)<=T-K','wings':'T-K<abs(m)<=T+K','inequalityEquivalence':'max lower/min upper is exactly conjunction; nonempty iff -K-T<=m<=K+T','boundaryMeasureZero':True})
    domainrecords=[]
    for radius,tail in ((27,122),(29,124)):
        MM=radius+tail;c=tail-radius
        for mv in (-MM-1,-MM,-MM+1,-c-1,-c,0,c,c+1,MM-1,MM,MM+1):
            window=C.finite_window(radius,tail,mv);domainrecords.append({'K':radius,'T':tail,'m':mv,'window':window})
            require(window['empty']==(abs(mv)>MM),'full retained outer wings')
    J.emit('new-exact-domain-witnesses',domainrecords)
    controls=[]
    for radius,tail in ((27,122),(29,124)):
        mv=tail+radius-F(1,2);kv=-radius;lv=radius
        good=C.mapped_contains(radius,tail,kv,lv,mv,'l');bad=C.mapped_contains(radius,tail,kv,lv,mv,'k')
        witness={'K':radius,'T':tail,'k':kv,'l':lv,'m':mv,'original_t':F(lv)-mv,'correctOutputWindow':good,'wrongInputWindow':bad}
        J.emit('control-window-'+str(radius),witness);require(good and not bad and C.rectangle_contains(radius,tail,kv,lv,F(lv)-mv),'actual wrong reflected window control');controls.append(witness)
        retained=C.mapped_contains(radius,tail,radius,0,mv,'k');pruned=retained and abs(mv)<=radius
        witness={'K':radius,'T':tail,'k':radius,'l':0,'m':mv,'original_t':mv-radius,'correct':retained,'omitWings':pruned}
        J.emit('control-wing-'+str(radius),witness);require(retained and not pruned,'actual full-wing omission control');controls.append(witness)
    for name,component,baseline_poly,mutant_poly in [('reflected-sign','Dr',k*(l+mid),k*(3*l-mid)),('J-polynomial','J',mid*(k+mid),mid*k)]:
        baseline=complete[component];common=sp.cancel(baseline/baseline_poly);mutant=common*mutant_poly
        exact('new-control-baseline-factor-'+name,baseline,common*baseline_poly,{'actualPrimitive':component})
        point={v:sp.Integer(1) for v in baseline.free_symbols|mutant.free_symbols};point.update({k:sp.Integer(1),l:sp.Integer(2),mid:sp.Integer(3)})
        J.start('control-'+name,{'actualFullBaseline':baseline,'actualFullMutant':mutant,'baselinePolynomial':baseline_poly,'mutantPolynomial':mutant_poly,'commonFactor':common,'predeterminedPoint':[[a,b] for a,b in point.items()]})
        movement=sp.cancel(sp.together(mutant-baseline));value=sp.cancel(movement.xreplace(point))
        record={'fullMovement':movement,'pointMovement':value,'formalAlgebraSensitivityOnly':True};J.emit('control-'+name+'-decision',record)
        require(value.is_Rational is True and bool(value!=0),'actual complete algebra control responds');J.finish(record);controls.append(record)
    # Transport every saved union label. We do not construct or approve a new
    # quadrature mesh: branch and profile intersections remain future obligations.
    labelrecords=[]
    for alias in m['geometryPlans']:
        old=get(alias);plan=[]
        for line in ([] if 'height' in alias else old['lines']):
            for primitive in ('J','Dr','Dh','Dq'):
                deps=sorted(str(v) for v in densities[('J','Dr','Dh','Dq').index(primitive)].free_symbols if v in (qi,qo,qh,qs))
                plan.append({'originalLine':line,'primitive':primitive,'map':'t=l-m' if primitive=='Dr' else 't=m-k','external_k_l_constraintUnchanged':True,'actualDepthDependencies':deps,'pairedMiddleRootCollisionPresent':qh in densities[('J','Dr','Dh','Dq').index(primitive)].free_symbols and qs in densities[('J','Dr','Dh','Dq').index(primitive)].free_symbols,'newBranchLabels':['m=-kappa','m=+kappa','k=-kappa','k=+kappa','l=-kappa','l=+kappa'],'newWindowLabels':['m=-K-T','m=-T+K','m=T-K','m=K+T'],'profileArguments':['m-k','l-m'],'numericalMeshReady':False})
        J.emit('new-label-transport-'+Path(alias).stem,{'originalPlan':old,'transport':plan,'heightPVUnchanged':'height' in alias,'noOldGeometryRecalculation':True});labelrecords.extend(plan)
    # Independent nonnegative real parts imply the outgoing quadrant inequalities.
    ux,uy,vx,vy=sp.symbols('contract_ux contract_uy contract_vx contract_vy',nonnegative=True)
    expr=(ux+vx)**2+(uy+vy)**2-(ux**2+uy**2+vx**2+vy**2)
    J.start('new-quadrant-polynomial',{'expression':expr,'variables':[ux,uy,vx,vy]})
    certificate=C.quadrant_polynomial_certificate(sp,sp.expand(expr),(ux,uy,vx,vy));J.emit('new-quadrant-polynomial-decision',certificate);J.finish(certificate)
    rr,ss=sp.symbols('contract_r contract_s',nonnegative=True)
    exact('new-quadrant-root-sum-bound',2*((ux+vx)**2+(uy+vy)**2)-(rr+ss)**2, (rr-ss)**2+4*(ux*vx+uy*vy)+2*(ux**2+uy**2-rr**2)+2*(vx**2+vy**2-ss**2),{'rSquared':'ux^2+uy^2','sSquared':'vx^2+vy^2','remainingTermsNonnegative':'(r-s)^2 and quadrant dot products'})
    radius=sp.Symbol('contract_kappa',positive=True);zz=sp.Symbol('contract_z',real=True)
    exact('new-simple-root-factorization',radius**2-zz**2,(radius-zz)*(radius+zz))
    J.emit('new-compact-integrability-lemma',{'actualBeta':beta0,'betaRealLower':sp.Rational(30,109),'betaImaginary':sp.Rational(9,109),'actualKappa':kap,'kappaDomain':'0<kappa<3; kappa^2=119/20','profile':'real A(z) obeys |A|<=1/(2*pi)<1 from |sinh y|>=|y|; continuous value A(0)=1/(2*pi)','qEnvelope':'|q(p)|=sqrt(abs(kappa^2-p^2))<=abs(p)+3','rootNeighborhoods':'|m-sign*kappa|<kappa/2 implies |m+sign*kappa|>=3*kappa/2','localMajorant':'sqrt(2/(3*kappa))*|m-sign*kappa|^(-1/2)','localIntegralFact':'integral_0^epsilon s^(-1/2) ds=2 sqrt(epsilon); assessed analytic fact, no quadrature','exceptionalSets':['m=+-kappa','simultaneous q(m)=q(k)=0 or q(m)=q(l)=0'],'zeroOverZeroAssigned':False,'fubini':'finite domains, measurable source-bound integrable absolute majorants; analytic application independently assessed, not machine measure theory','positiveRegulatorUsed':False})
    J.emit('complete-address-certificates',address_results)
    return {'status':'BOUNDED_FINITE_WINDOW_CONTRACTION_CERTIFICATE_COMPLETE','selectedAddresses':len(address_results),'templates':len(injections),'fields':len(fields),'controls':controls,'primitiveIdentities':4,'wrongRootAlgebraJoined':True,'originalWindows':[[27,122],[29,124]],'fullWingsRetained':True,'allGeometryLabelsPreserved':True,'labelRecords':len(labelrecords),'numericalEvaluatorReady':False,'integralsEvaluated':False,'completedFunctionsReplayed':False,'limits':['Analytic absolute-integrability and Fubini reasoning are assessed mathematics, not machine measure theory.','Inference-dependent gamma units and accepted original algebra remain dependencies.','No discrete quadrature identity, new nodes, interpolation, uniform error, cost/runtime estimate, integral/action/current/loss or inverse.','Full numerical discretization, request identity, nested error propagation, storage and complete independent routes remain required.']}

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
    try:
        for i,(path,digest) in enumerate(documents.items()):
            dest=args.out/('identity-'+str(i)+'-'+Path(path).name);require(sha(path)==digest,'identity source');shutil.copyfile(path,dest);require(sha(dest)==digest,'identity copy');identities[path]={'path':dest.name,'sha256':digest,'bytes':dest.stat().st_size}
        J.emit('additional-identity-copies',identities);J.emit('actual-containment',enforced)
        import sympy as sp
        require(sp.__version__==m['sympyVersion'],'pinned symbolic runtime version')
        U=load_library(m);result=run(m,J,U)
    except BaseException:
        error=traceback.format_exc();result={'status':'FAILED_PRESERVED','failure':error,'activeOperation':J.active,'completedOperations':J.completed};J.emit('failure',result)
    finally:
        post={p:posthash(p,h) for p,h in {**m['sourcePins'],**documents}.items()};identity_post={v['path']:posthash(args.out/v['path'],v['sha256']) for v in identities.values()}
        copied={str(p.relative_to(args.out)):sha(p) for p in (args.out/'saved').rglob('*') if p.is_file()};expected={'saved/'+a:r['sha256'] for a,r in m['savedInputs'].items()}
        J.emit('posthashes',{'sources':post,'identities':identity_post,'copied':copied,'expectedCopied':expected,'copiesIntact':copied==expected})
        result.update(wallMilliseconds=round((time.monotonic()-start)*1000),sourcePosthashesIntact=all(v['intact'] for v in post.values()),copiesIntact=copied==expected,identityCopiesIntact=len(identities)==len(documents) and all(v['intact'] for v in identity_post.values()),completedOperations=J.completed,activeOperation=J.active)
        J.emit('journal-result',result);save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text());sys.stdout.flush()
    if error:sys.stderr.write(error);return 1
    require(result['sourcePosthashesIntact'] and result['copiesIntact'] and result['identityCopiesIntact'],'posthash integrity');return 0

if __name__=='__main__':sys.exit(main())
