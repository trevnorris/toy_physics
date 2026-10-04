#!/usr/bin/env python3
"""Guarded source-bound kernel, physical wave and measure unit transport.

No CAS restoration, producer, grade/field/integral or geometry replay. The strict
constructor-text dimensional interpreter runs only after shared containment.
"""
import argparse,ast,hashlib,importlib.util,json,math,os,resource,shutil,sys,time,traceback
from fractions import Fraction as F
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO={'text':'0','srepr':'Integer(0)'}

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
    require(g['status']=='READY_FOR_ONE_PACKET_KERNEL_UNITS' and g['independentBuildClearance'] is True,'actual independent build readiness')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==m['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','library','authority','buildReviewRecord'):
        require(sha(g[key])==g[key+'Sha256'],'gate pin '+key)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual guard/supervisor')
    require(g['launcher']==m['launcher'] and g['library']==m['librarySource'] and g['authority']==m['executionAuthority'] and g['buildReviewRecord']==m['reviewRecordWillBe'],'actual document paths')
    r=read(g['buildReviewRecord']);require(r['allChecksPassed'] is True and r['independentBuildClearance'] is True,'build assessment')
    for key in ('workerSha256','manifestSha256','librarySha256','launcherSha256','sharedGuardSha256','supervisorSha256'):require(r[key]==g[key],'actual reviewed '+key)
    require(all(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED KERNEL WAVE/MEASURE UNIT-TRANSPORT BUILD' for e in ('claude','grok')),'both literal build reports')
    method=read(m['methodRecord']);require(g['methodRecordSha256']==sha(m['methodRecord']) and method['jointIndependentMethodClearance'] and method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'actual pressure-readiness method')
    a=read(g['authority']);require(a['scope']==g['scope']==m['scope'] and a['boundedInstrumentAuthorized'] and a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and a['automaticScientificRetry'] is False and a['noDeadline'] and g['durationLimits'] is None,'standing bounded authority')
    return g


def load_library(m):
    spec=importlib.util.spec_from_file_location('packet_source_unit_transport',m['librarySource']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def run(m,J,K):
    U=K.U;raw,copies=copy_inputs(m,J)
    def get(name):return raw[name]
    contracts=get('source-contracts.json')
    for rec in contracts['fragments'].values():K.contract(rec)
    J.emit('actual-source-contracts',contracts)
    for rec in contracts['fragments'].values():
        origins=[]
        for name in ('preflight','inner','fourier','reference'):
            previous=get('origins/'+name+'-posthashes.json')
            if rec['source'] in previous:
                item=previous[rec['source']]
                require(item['expected']==item['actual']==sha(rec['source']) and item.get('error') is None,'source matches executed historical snapshot')
                origins.append({'run':name,'posthash':item})
        require(bool(origins),'actual old source origin '+rec['name'])
        J.emit('source-origin-'+rec['name'],{'source':rec,'originalRunPins':origins})
    source_result=get('source/result-record.json')
    require(source_result['allChecksPassed'] is True and source_result['status']=='BOUNDED_SOURCE_GRADE_PROFILE_UNIT_TRANSPORT_ACCEPTED_KERNEL_MEASURE_PENDING','accepted source transport')
    require(get('source/journal-result.json')['sourceGradeFieldUnitTransportComplete'] is True,'completed source transport')
    # No scientific interpretation of copied bytes was done during preparation.
    for alias,r in m['savedInputs'].items():
        if alias.startswith('source/') and alias!='source/result-record.json':
            old=source_result['records']['complete/'+alias.removeprefix('source/')]
            require((old['sha256'],old['bytes'])==(r['sha256'],r['bytes']),'accepted source file '+alias)
        for run_name in ('preflight','inner','fourier'):
            if alias.startswith(run_name+'/') and alias not in (run_name+'/artifact-index.json',run_name+'/result-record.json'):
                old=get(run_name+'/artifact-index.json')[alias.removeprefix(run_name+'/')]
                require((old['sha256'],old['bytes'])==(r['sha256'],r['bytes']),'original indexed operand '+alias)
    context=get('native/binding.json');registry=get('native/effective-registry.json')['units']
    require(get('source/transport-context.json')['actual']['restored']==context and get('source/transport-context.json')['registry']==registry,'same actual binding and unit registry')
    require(context['frequency']=={'text':'3','srepr':'Integer(3)'} and context['physicalInput']==get('physical-input.json'),'same held omega and original physical input')
    require(context['physicalInput']['parameters']['omega']=='1','original omega1 separate from held3')
    for name,value in context['numeric'].items():
        old=get('source/numeric-binding-'+name+'-return.json')
        J.emit('inherited-binding-'+name,{'return':old,'original':get('source/numeric-binding-'+name+'-input.json'),'functionCalled':False})
        require(old['symbol']==name and old['value']==value and K.unit(old['originalUnit'])==K.unit(registry[name]),'actual numeric value and original unit binding identity')
    def inherited(prefix):
        record=K.inherit_zero(get(prefix+'-input.json'),get(prefix+'-return.json'))
        J.emit('restore-'+prefix.replace('/','--'),record);return record['input']
    def walk(name,n,env,expected,functions=None,params=None,origin=None):
        J.start(name,{'expression':ast.unparse(n) if isinstance(n,ast.AST) else n,'environment':env,'functions':functions or {},'parameters':params or {},'expected':expected,'sourceOrigin':origin,'newDimensionalInterpretationOnly':True})
        w=K.Walk(env,functions,params)
        try:v=w.dim(n)
        finally:J.emit(name+'-walk',w.events)
        r={'unit':v,'expected':expected};J.emit(name+'-decision',r);require(v==K.unit(expected),'homogeneous source formula '+name);J.finish(r);return v
    def source_expr(name):return K.assignment_expr(contracts['fragments'][name])[1]
    base={'omega':K.unit(registry['omega'])};params={n:registry[n] for n in ('rho_m','Lambda_A_0','tau_A')}
    mu=walk('new-mu-units',source_expr('preflight-mu'),base,(-4,-1,1),params=params,origin=contracts['fragments']['preflight-mu'])
    aa=walk('new-a-units',source_expr('preflight-aa'),base,(3,1,-1),params=params,origin=contracts['fragments']['preflight-aa'])
    beta=walk('new-beta-units',source_expr('preflight-beta'),{'aa':aa,'mu':mu},(-1,0,0),origin=contracts['fragments']['preflight-beta'])
    for name in ('runtime-mu','runtime-a'):inherited('inner/'+name)
    physical=get('preflight/physical-plan.json')
    require(physical['frequency']==3 and physical['cs']['srepr']=='Mul(Rational(1, 2), Pow(Integer(6), Rational(1, 2)))','same selected physical speed')
    require(K.unit(registry['W_0'])==K.LENGTH and K.unit(registry['L_W'])==K.LENGTH and K.unit(registry['rho_m'])==(-4,0,1),'actual native width/length/bulk density')
    # Numeric values remain physical magnitudes in their original unit frame.
    J.emit('parameter-quantity-origins',{'context':context,'registry':registry,'physicalPlan':physical,'oldRuntimeProofs':[get('inner/runtime-'+n+'-input.json') for n in ('mu','a')],'boundNumbersAssignedDimensions':False})
    geom=get('native/lower-geometry.json');geom_units={}
    for key,expected in [('labHeight',K.LENGTH),('outwardSlopeDefinition',K.ZERO)]:
        J.start('native-'+key+'-units',{'original':geom[key],'registry':registry,'expected':expected})
        w=U.Units(registry)
        try:v=w.dimension(U.parse(geom[key]['srepr']))
        finally:J.emit('native-'+key+'-walk',w.events)
        J.emit('native-'+key+'-decision',{'unit':v});require(v==expected,'native geometry dimensions');J.finish({'unit':v});geom_units[key]=v
    profile=get('inner/physical-H-and-PV-lemma.json')
    require(profile['nativeBinding']==context and profile['nativeGeometry']==geom and profile['nativeFourierContract']==get('native/fourier-contract.json'),'actual native profile/normal/Fourier join')
    for name in m['inheritedInnerProofs']:inherited('inner/'+name)
    require(get('inner/new-native-slope-scale-operands.json')['savedScale']==get('native/profile-scale.json'),'original derivative L scale')
    first=get('native-input/reference/first-shape-native-transport.json')
    J.start('new-original-first-shape-3D-units',{'original':first,'registry':registry,'nativeMeasure':get('native/fourier-contract.json')})
    uw=U.Units(registry)
    try:firstunit=uw.dimension(U.parse(first['restored']['srepr']))
    finally:J.emit('new-original-first-shape-unit-walk',uw.events)
    J.emit('new-original-first-shape-decision',{'unit':firstunit,'edgeDeltaReduction':[-2,0,0]})
    require(firstunit==(0,-1,1),'actual unbound native three-dimensional first shape')
    J.finish({'nativeUnit':firstunit,'reducedUnit':K.add(firstunit,(-2,0,0))})
    profiles=get('reference/restored-profile-arguments.json');ordered=get('raw-direct/raw-ordered-before-cancel.json')
    require(profiles['original']==ordered,'same actual ordered pre-cancel profile operands')
    for name in ('height-scale-join','profile-phase-join','right-height-PV-density'):inherited('reference/'+name)
    for name in ('height-profile-factor-join','slope-profile-factor-join'):
        original=get('reference/'+name+'-original-input.json');rewrite=get('reference/'+name+'-argument-expansion.json')
        proof=inherited('reference/'+name+'-canonical')
        J.emit('inherited-profile-'+name,{'original':original,'exactArgumentRewrite':rewrite,'canonicalProof':proof,'newProfileCalculation':False})
    J.emit('original-W-L-profile-origins',{'ordered':ordered,'restoredProfileArguments':profiles,'actualGenericScaleSource':contracts['fragments']['raw-direct-scale'],'firstShapeTransport':first,'widthBinding':get('source/numeric-binding-W_0-return.json'),'lengthBinding':get('source/numeric-binding-L_W-return.json'),'inheritedScaleProof':get('reference/height-scale-join-input.json')})
    tags=get('pressure/whole-tags.json')
    for name in ('H','Jwhole','Dwhole'):
        require(tags[name]['savedDefinition']==get('whole-origin/'+name+'.json'),'full actual whole definition '+name)
        J.emit('inherited-whole-'+name,{'tag':tags[name],'completeOriginal':get('whole-origin/'+name+'.json'),'functionsCalled':False})
    # h has the native height unit, j the independent native slope unit. Fourier
    # transformation adds one physical length; the native edge deltas are retained.
    hhat=K.add(geom_units['labHeight'],K.LENGTH);jhat=K.add(geom_units['outwardSlopeDefinition'],K.LENGTH)
    env={n:K.MOMENTUM for n in ('k','l','t','s','Q','qi','qo','qh','qs','qm')}
    env.update(W=K.unit(registry['W_0']),L=K.unit(registry['L_W']),I=K.ZERO,mu=mu,a=aa,aa=aa,beta=beta)
    anode=source_expr('reference-A');require(isinstance(anode,ast.Lambda) and len(anode.args.args)==1,'actual one-variable A profile')
    functions={'A':{'arguments':[anode.args.args[0].arg],'body':ast.unparse(anode.body)}}
    Aunit=walk('new-profile-A-units',anode.body,{**env,anode.args.args[0].arg:K.MOMENTUM},K.ZERO,origin=contracts['fragments']['reference-A'])
    pv=source_expr('reference-originalPV')
    require(isinstance(pv,ast.BinOp) and isinstance(pv.op,ast.Mult) and isinstance(pv.left,ast.BinOp) and isinstance(pv.left.op,ast.Mult),'original ordered height-times-slope syntax')
    hpv=walk('new-height-PV-units',pv.left.right,env,(2,0,0),functions,origin=contracts['fragments']['reference-originalPV'])
    jpv=walk('new-slope-transform-units',pv.right,env,(1,0,0),functions,origin=contracts['fragments']['reference-originalPV'])
    require(hpv==hhat and jpv==jhat,'native height/slope and source PV transform units')
    # The actual generic H source assignments retain W/L before their numerical bindings.
    hsub=walk('new-H-subtracted-units',source_expr('reference-Hsub'),env,(2,0,0),functions,origin=contracts['fragments']['reference-Hsub'])
    hnode=K.contract(contracts['fragments']['reference-H-record'])
    require(isinstance(hnode,ast.Expr) and isinstance(hnode.value,ast.Call),'actual H persistence call')
    hdict=hnode.value.args[1];require(isinstance(hdict,ast.Dict),'full H definition dictionary')
    hd={key.value:value for key,value in zip(hdict.keys,hdict.values) if isinstance(key,ast.Constant)}
    Hcontact=walk('new-H-contact-units',hd['contact'],env,(2,0,0),functions,origin=contracts['fragments']['reference-H-record'])
    Hdensity=walk('new-H-PV-density-units',hd['ordinarySubtractedIntegrand'],{**env,'Hsub':hsub},(3,0,0),functions,origin=contracts['fragments']['reference-H-record'])
    require(K.add(Hdensity,K.MOMENTUM)==Hcontact,'same H contact and integrated PV dimensions')
    # Extract the actual W/4 coefficient from the persisted generic H contact,
    # before multiplication by jhat. Its own unit is L, not hhat's L^2.
    contact_node=hd['contact'];require(isinstance(contact_node,ast.BinOp) and isinstance(contact_node.op,ast.Mult),'actual H contact product')
    heightcontact=walk('new-height-contact-coefficient-units',contact_node.left,env,(1,0,0),functions,origin=contracts['fragments']['reference-H-record'])
    require(K.add(heightcontact,K.LENGTH)==hpv,'delta and PV pieces have the same height-transform unit')
    # Source-level function identity and its complete old proof operands are joined.
    # No invocation of kernel_components or a numerical route occurs here.
    fn=K.contract(contracts['fragments']['inner-kernel-components'])
    require([x.arg for x in fn.args.args]==['k','l','t','qi','qo','qh','qs','A1','A2','a','mu','beta','W','L','I'],'actual generic kernel signature')
    J.start('new-inner-density-units',{'functionSource':contracts['fragments']['inner-kernel-components'],'actualRuntimeCall':contracts['fragments']['inner-runtime-call'],'actualProofCall':contracts['fragments']['inner-proof-call'],'environment':{**env,'A1':Aunit,'A2':Aunit},'oldJProof':get('inner/new-runtime-J-arithmetic-input.json'),'oldDProof':get('inner/new-runtime-D-arithmetic-sum-input.json')})
    dw=K.Walk({**env,'A1':Aunit,'A2':Aunit})
    try:densities,dw=K.statements(fn,{**env,'A1':Aunit,'A2':Aunit},walk=dw)
    finally:J.emit('new-inner-density-walk',dw.events)
    J.emit('new-inner-density-decision',{'units':densities,'components':['J','reflected','height','quadratic']})
    require(len(densities)==4 and all(v==(-1,-1,1) for v in densities),'each added density separately homogeneous')
    J.finish({'units':densities,'middleMeasure':K.MOMENTUM,'wholeUnits':[K.add(v,K.MOMENTUM) for v in densities],'notRootProduct':True})
    # All four q arguments remain distinct in the actual runtime call and saved maps.
    call=K.contract(contracts['fragments']['inner-runtime-call'])
    proofcall=K.contract(contracts['fragments']['inner-proof-call'])
    J.emit('actual-kernel-call-and-definition-joins',{'runtime':ast.unparse(call),'proof':ast.unparse(proofcall),'genericSource':contracts['fragments']['inner-kernel-components'],'restoredNumericalAdapterProofs':[get('inner/new-runtime-'+n+'-arithmetic'+('-sum' if n=='D' else '')+'-input.json') for n in ('J','D')],'wholeDefinitions':get('pressure/whole-definitions.json'),'oldFunctionCalled':False,'sourceAndAcceptedAlgebraDependency':True})
    # Exact dictionary expressions used by the accepted numeric adapter, with
    # physical origins for all previously bound parameters and whole tags.
    templates=K.template_nodes(contracts['fragments']['preflight-templates'])
    tenv={**env,'hvalue':hhat,'jvalue':jhat,'Hvalue':Hcontact,'Jvalue':K.add(densities[0],K.MOMENTUM),'Dvalue':K.add(densities[1],K.MOMENTUM)}
    basic={name:walk('new-template-'+name,node,tenv,(-3,-1,1) if name=='NATIVE_FLAT' else (-2,-1,1),origin=contracts['fragments']['preflight-templates']) for name,node in templates.items()}
    adapters=get('preflight/numeric-factor-adapters.json');selected=get('selected/pressure-addresses.json')['selected']
    require(len(adapters['definitions'])==20 and len(selected)==544 and len(adapters['addressJoins'])==544,'complete template and address counts')
    template_units={}
    for label,entry in adapters['definitions'].items():
        face,slot,component=label.split('-',2);args=get('preflight/numeric-factor-'+label+'-arguments.json');proof=inherited('preflight/numeric-factor-'+label)
        require(entry['firstAddressId']==args['address']['addressId'] and proof['left']==entry['mapped'] and proof['right']==entry['template'] and args['actualCompleteFactor']==entry['original'],'actual full adapter and proof operands')
        dimension=basic[component]
        if slot=='normal':dimension=K.add(dimension,K.MOMENTUM)
        require(slot in ('normal','pressure') and face in ('plus','minus'),'native face/slot')
        template_units[label]=dimension
        J.emit('complete-template-'+label,{'sourceExpression':ast.unparse(templates[component]),'sourceNormalInsertion':contracts['fragments']['preflight-normal'],'face':face,'normalSign':(1 if face=='plus' else -1) if slot=='normal' else None,'actualArguments':args,'actualAdapter':entry,'unit':dimension,'normalDepth':'q(l)' if slot=='normal' else None,'normalFactorCount':1 if slot=='normal' else 0,'wholeDAdditionalMiddleIntegral':False,'wholeDAdditionalResolvents':False})
    Fourier=get('native/fourier-contract.json');convention=get('fourier/transform-convention.json')
    require(Fourier['profileForwardPower']==-3 and Fourier['coordinates']==3 and Fourier['invariantEdgeCoordinates']==2 and Fourier['sourceForwardPower']==0 and Fourier['sourceInversePower']==-3,'actual native edge/2pi normalization')
    require(convention['sourceDerivativeBeforeCoefficient'] is True and convention['testNoConjugation'] is True,'actual Fourier derivative/bilinear convention')
    J.emit('native-Fourier-measure-transport',{'native':Fourier,'actualTransform':convention,'actualProductSource':contracts['fragments']['fourier-product'],'remainingProfileForwardPower':-1,'sourceXFactor':'1/(2*pi)','testYFactor':'1 = 2*pi times normalized forward transform at -l','plainMiddleMeasure':True,'measureRule':'One physical spatial integration adds L; one physical momentum integration subtracts L; delta(momentum) adds L.','dualAmplitude':'fixed formal THETA dual factor carried outside unit triples','noFourierFunctionCalled':True})
    families=[get(alias) for alias in m['fourierFamilies']]
    # Restore the already-responsive source/profile control, never call it again.
    J.emit('restored-profile-control',{'input':get('source/control-missing-profile-L-input.json'),'return':get('source/control-missing-profile-L-return.json'),'functionCalled':False})
    require(get('source/control-missing-profile-L-return.json')['refused'] is True,'completed profile control return')
    results=[]
    for index,a in enumerate(selected):
        ident=a['addressId'];label='-'.join((a['face'],a['slot'],a['component']));old=get('source/address-'+str(ident)+'-input.json');sr=get('source/address-'+str(ident)+'-return.json');wave=inherited('preflight/address-'+str(ident)+'-wave-jet')
        J.start('summand-'+str(ident),{'address':a,'sourceTransportInput':old,'sourceTransportReturn':sr,'waveProof':wave,'actualWaveSource':contracts['fragments']['preflight-wave-proof'],'adapter':adapters['definitions'][label],'normal':a['normalMultiplier'],'requiredDimension':[-2,-1,1]})
        require(old['address']==a and sr['addressId']==ident and sr['nativeSourceAndSavedGradeProfileJoined'] is True,'same accepted source/consumer actual address')
        require(adapters['addressJoins'][index]=={'addressId':ident,'adapter':label},'actual complete adapter routing')
        require(wave['left']==a['waveMultiplier'],'actual inherited wave factor')
        fac=get('factors/'+a['fullFactorProof']['proof']+'-operands.json')
        require(fac['mappedAddressFactor']==adapters['definitions'][label]['original'],'same accepted native full factor')
        require(a['epsilonCount']==(0 if a['status'].startswith('EXACT_ZERO') else 1),'epsilon exactly once or explicit zero')
        jet=a['jet'];require(jet['channel']=='e_W' and len(jet['spatialOrders'])==3 and all(type(n) is int and n>=0 for n in [jet['timeOrder'],*jet['spatialOrders']]),'actual independent physical jet')
        jdim=(-sum(jet['spatialOrders']),-jet['timeOrder'],0)
        route=K.address_route(a)
        r=K.assemble_route(sr['sourceRequiredUnit'],jdim,sr['consumerRequiredUnit'],template_units[label],route)
        # Native numeric wave factor was joined above; its accepted complete
        # p/time/edge expression keeps derivative-before-coefficient order.
        r.update(addressId=ident,status=a['status'],flatSupport=a['responseMap']['flatSupport'],normalMap=K.normal_map(a),sourceUnitProof='inherited actual source address return',zeroHasIntrinsicDimension=False,nonzero=not a['status'].startswith('EXACT_ZERO'))
        J.emit('summand-'+str(ident)+'-decision',r)
        require(K.pressure_dimension(r),'full addressed pressure dimension')
        # Check finite-bank family signatures only as interface evidence. No
        # equality of argument values or cache reuse is claimed by this join.
        matches=[]
        for fam in families:
            sp=fam['spec'] if 'spec' in fam else fam
            if sp.get('role')=='X' and sp.get('fieldId')==a['sourceTransform']['coefficientId'] and sp.get('timeOrder')==jet['timeOrder'] and sp.get('spatialOrders')==jet['spatialOrders']:matches.append(sp)
        r['matchingXFamilyInterfaces']=matches
        ymatches=[f for f in families if f['spec']['role']=='Y' and f['spec']['argumentDerivative']==0 and f['spec']['fieldId']==a['consumerTransform']['coefficientId']]
        r['matchingYFamilyInterfaces']=ymatches
        if r['nonzero']:require(bool(ymatches),'live consumer Fourier family interface')
        if a['component']=='NATIVE_HEIGHT':
            coefficientUnit=K.add(template_units[label],K.scale(hhat,-1))
            contactKernel=K.add(coefficientUnit,heightcontact)
            pvKernel=K.add(coefficientUnit,hpv)
            contactTotal=K.total_unit(sr['sourceRequiredUnit'],jdim,sr['consumerRequiredUnit'],contactKernel,K.MOMENTUM)
            pvTotal=K.total_unit(sr['sourceRequiredUnit'],jdim,sr['consumerRequiredUnit'],pvKernel,K.scale(K.MOMENTUM,2))
            r['heightOwnRoutes']={'contactSingleK':contactTotal,'pairedPV_k_positiveQ':pvTotal,'QJacobianUnit':K.ZERO,'oppositeQBranchesSameUnits':True,'argumentDerivativeYUnit':K.add(r['Y'],K.LENGTH),'actualChiCancellation':get('inner/new-height-chi-cancellation-input.json'),'actualProfileLemma':profile['savedProfileLemma']}
            require(contactTotal['total']==pvTotal['total']==r['total'],'actual separate height contact and paired-PV measure units')
        if r['nonzero']:require(bool(matches),'live source Fourier family interface')
        J.finish(r);results.append(r)
    controls=[]
    for comp,unit0 in zip(('J','reflected','height','quadratic'),densities):
        if comp not in ('J','reflected'):continue
        J.start('control-missing-dt-'+comp,{'actualDensitySource':contracts['fragments']['inner-kernel-components'],'densityUnit':unit0,'baselineMeasure':K.MOMENTUM,'mutatedMeasure':K.ZERO,'mutation':'omit the actual middle dt in this density route'})
        baseline=K.add(unit0,K.MOMENTUM);mutant=K.add(unit0,K.ZERO);r={'baseline':baseline,'mutant':mutant,'responded':baseline!=mutant,'notNumericalResponse':True};J.emit('control-missing-dt-'+comp+'-decision',r);require(r['responded'],'missing measure response');J.finish(r);controls.append(r)
    J.start('control-a-unitless',{'actualSource':contracts['fragments']['inner-kernel-components'],'baselineEnvironment':{**env,'A1':Aunit,'A2':Aunit},'mutation':{'a':K.ZERO}})
    mutated,w=K.statements(fn,{**env,'A1':Aunit,'A2':Aunit,'a':K.ZERO});J.emit('control-a-unitless-walk',w.events);r={'baseline':densities,'mutant':mutated,'responded':mutated[0]!=densities[0],'notNumericalResponse':True};J.emit('control-a-unitless-decision',r);require(r['responded'],'actual J dimensional a control');J.finish(r);controls.append(r)
    candidate=next(v for v in results if v['nonzero'] and v['flatSupport'] is not None)
    a=next(v for v in selected if v['addressId']==candidate['addressId']);sr=get('source/address-'+str(a['addressId'])+'-return.json')
    label='-'.join((a['face'],a['slot'],a['component']));jet=a['jet'];jdim=(-sum(jet['spatialOrders']),-jet['timeOrder'],0)
    J.start('control-flat-delta',{'actualAddress':a,'actualSourceTransport':sr,'actualAdapter':adapters['definitions'][label],'actualSummand':candidate,'jetUnit':jdim,'coefficientKernelUnit':template_units[label],'requiredDimension':[-2,-1,1],'mutation':'Keep actual flat coefficient and both dl dk integrations; remove only delta(k-l) from the supported unreduced route.'})
    r=K.flat_delta_control(a,sr['sourceRequiredUnit'],jdim,sr['consumerRequiredUnit'],template_units[label],candidate)
    J.emit('control-flat-delta-decision',r);require(r['responded'],'actual supported flat delta control');J.finish(r);controls.append(r)
    J.emit('complete-summand-dimensions',results)
    return {'status':'BOUNDED_PACKET_KERNEL_WAVE_MEASURE_UNIT_TRANSPORT_COMPLETE','selectedAddresses':len(results),'formalAddresses':sum(v['nonzero'] for v in results),'templates':len(template_units),'newControls':controls,'restoredProfileControl':True,'pressureSummandUnitsComplete':True,'numericPolynomialIndependentlyDimensional':False,'numericalEvaluatorReady':False,'newIntegralOrAction':False,'priorFunctionsReplayed':False,'limits':['Dimensional homogeneity and original accepted algebra, including inferred gamma units, remain dependencies.','Unit consistency and sensitivity do not prove numerical value or common kernel convention.','No uniform quadrature error, exact cache matches, full outer evaluator, action/current/loss or sweep.']}

def main():
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
        J.emit('additional-identity-copies',identities);J.emit('actual-containment',enforced);U=load_library(m);result=run(m,J,U)
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
