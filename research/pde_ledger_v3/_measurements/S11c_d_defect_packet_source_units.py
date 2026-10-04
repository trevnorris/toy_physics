#!/usr/bin/env python3
"""Guarded saved source-grade/profile unit transport; kernel and measure units pending.

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
    require(g['status']=='READY_FOR_ONE_PACKET_SOURCE_UNITS' and g['independentBuildClearance'] is True,'actual independent build readiness')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==m['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','library','authority','buildReviewRecord'):
        require(sha(g[key])==g[key+'Sha256'],'gate pin '+key)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual guard/supervisor')
    require(g['launcher']==m['launcher'] and g['library']==m['librarySource'] and g['authority']==m['executionAuthority'] and g['buildReviewRecord']==m['reviewRecordWillBe'],'actual document paths')
    r=read(g['buildReviewRecord']);require(r['allChecksPassed'] is True and r['independentBuildClearance'] is True,'build assessment')
    for key in ('workerSha256','manifestSha256','librarySha256','launcherSha256','sharedGuardSha256','supervisorSha256'):require(r[key]==g[key],'actual reviewed '+key)
    require(all(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED SOURCE UNIT-TRANSPORT BUILD' for e in ('claude','grok')),'both literal build reports')
    method=read(m['methodRecord']);require(g['methodRecordSha256']==sha(m['methodRecord']) and method['jointIndependentMethodClearance'] and method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'actual pressure-readiness method')
    a=read(g['authority']);require(a['scope']==g['scope']==m['scope'] and a['boundedInstrumentAuthorized'] and a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and a['automaticScientificRetry'] is False and a['noDeadline'] and g['durationLimits'] is None,'standing bounded authority')
    return g


def load_library(m):
    spec=importlib.util.spec_from_file_location('packet_source_unit_transport',m['librarySource']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def certify(J,U,name,text,registry,expected):
    J.start(name,{'constructor':text,'registry':registry,'expected':expected,'newUnitCheckNotSourceEvaluation':True})
    checker=U.Units(registry)
    try:unit=checker.dimension(U.parse(text))
    except BaseException as error:
        J.emit(name+'-decision-operands',{'status':'UNIT_WALK_REFUSED','expected':expected,'errorType':type(error).__name__,'error':str(error),'partialWalk':True})
        raise
    finally:J.emit(name+'-unit-walk',checker.events)
    result={'unit':unit,'expected':expected,'literalZeroUnitUnspecified':unit is None}
    J.emit(name+'-decision-operands',result)
    require(unit==tuple(map(F,expected)),'native homogeneity '+name)
    J.finish(result);return unit


def join_native_case(J,U,name,original,case,table_tags=()):
    J.start(name,{'originalConstructor':original,'savedCase':case,'tableTags':list(table_tags),'constructorExecuted':False})
    try:
        table=U.parse(original)
        for tag in table_tags:table=U.tagged(table,tag)
        index,payload,census=U.select_case(table,case['case'])
        value=U.tagged(payload,'VALUE')
        result={'originalLabels':census,'selectedIndex':index,'savedIndex':case['caseIndex'],
            'selectedPayload':ast.unparse(payload),'selectedValue':ast.unparse(value),
            'payloadMatches':U.same(payload,U.parse(case['caseConstructorText'])),
            'valueMatches':U.same(value,U.parse(case['valueConstructorText'])),
            'savedPayloadHashMatches':hashlib.sha256(case['caseConstructorText'].encode()).hexdigest()==case['caseConstructorSha256']}
    except BaseException as error:
        J.emit(name+'-decision-operands',{'status':'SOURCE_CASE_REFUSED','errorType':type(error).__name__,'error':str(error)})
        raise
    J.emit(name+'-decision-operands',result)
    require(index==case['caseIndex'] and result['payloadMatches'] and result['valueMatches'] and result['savedPayloadHashMatches'],'actual native case labels/index/payload/value '+name)
    J.finish(result);return value


def run(m,J,U):
    raw,copies=copy_inputs(m,J)
    def get(name):return raw[name]
    def sc(name):return get('composition/'+name+'.json')
    def prior(name,left=None,right=None):
        args=sc(name+'-input');ret=sc(name+'-return')
        J.emit('inherited-'+name,{'input':args,'return':ret,'functionCalled':False})
        require(ret['cancelled']==ZERO,'literal inherited zero '+name)
        if left is not None:require(U.equal(args['left'],left),'left operand join '+name)
        if right is not None:require(U.equal(args['right'],right),'right operand join '+name)
        return args
    native_result=get('native/result-record.json')
    require(native_result['allChecksPassed'] is True and native_result['status']=='BOUNDED_NATIVE_PRESSURE_UNIT_BRIDGE_ACCEPTED_TRANSPORT_AND_EVALUATOR_PENDING','accepted native bridge')
    native_files=native_result['records'];composition_files=get('composition-files.json');composition_result=get('composition-result.json')
    require(composition_result['allInspectionChecksPassed'] is True and composition_result['status']=='SUPPORTED_BOUNDED_SOURCE_CONSUMER_COMPOSITION_INVENTORY','accepted composition inventory')
    identity_joins=[]
    for alias,record in m['savedInputs'].items():
        if alias.startswith('composition/'):old=composition_files[alias.removeprefix('composition/')]
        elif alias.startswith('native/') and alias!='native/result-record.json':old=native_files['complete/'+alias.removeprefix('native/')]
        else:continue
        require(old['sha256']==record['sha256'] and old['bytes']==record['bytes'],'actual accepted operand receipt '+alias)
        identity_joins.append({'alias':alias,'acceptedReceipt':old})
    J.emit('actual-accepted-file-joins',identity_joins)
    context=get('original/binding-context.json');actual=sc('actual-binding-context')
    require(get('native/original-input-context.json')['context']==context,'native certificate and composition binding context')
    require(actual['restored']==context and context['physicalInput']==get('original/physical-input.json'),'same physical bindings')
    require(context['frequency']=={'text':'3','srepr':'Integer(3)'} and actual['noSigmaEtaIdentification'] is True,'real omega and independent grades')
    registry=get('native/effective-registry.json')['units']
    J.emit('transport-context',{'actual':actual,'registry':registry,'nativeResult':native_result['status'],'interpretation':'Parameter numbers are values in the original units, never new dimensionless physical quantities.'})
    # These are missing substitution checks, not a replay of native source homogeneity.
    substitution_results=[]
    for name,value in {**context['profileEqualities'],**context['densityMap']}.items():
        require(name not in ('eta_bg','sigma_W','epsilon_shape'),'no grade/epsilon redefinition')
        require(name in registry,'original substitution unit')
        unit=certify(J,U,'new-substitution-'+name,value['srepr'],registry,registry[name]);substitution_results.append({'name':name,'unit':unit})
    for name in ('eta_bg','sigma_W','epsilon_shape'):
        require(U.unit_tuple(registry[name])==U.ZERO,'dimensionless independent grades and epsilon')
    # A numeric parameter binding preserves its physical unit annotation. No arithmetic on its value is used to infer units.
    numeric=[]
    for name,value in context['numeric'].items():
        require(name in registry,'unaccounted original numeric symbol '+name)
        require(not U.symbols(value['srepr']),'constant saved parameter binding '+name)
        numeric.append({'symbol':name,'originalUnit':registry[name],'value':value,'dimensionlessValueInOriginalUnits':True})
    J.emit('numeric-binding-provenance',numeric)
    scale=sc('native-profile-scale-join');require(scale['savedLength']==context['numeric']['L_W'] and scale['physicalLength']==context['physicalInput']['parameters']['L_W'],'actual original length join')
    U.verify_contracts(m,get('source-contracts.json'))
    J.emit('actual-source-contracts',get('source-contracts.json'))
    # Preserve complete unit walks; select reciprocal bases by their actual AST paths.
    native_units={};denominators={}
    for name in m['inheritedUnitOperations']:
        args=get('native/'+name+'-input.json');ret=get('native/'+name+'-return.json');walk=get('native/'+name+'-unit-walk.json')
        J.start('restore-'+name,{'input':args,'return':ret,'walk':walk})
        require(U.unit_tuple(ret['unit'])==U.unit_tuple(ret['expected']),'inherited native unit return')
        require(U.same(U.parse(args['constructor']),U.parse(walk[-1]['constructor'])) and U.unit_tuple(walk[-1]['unit'])==U.unit_tuple(ret['unit']),'actual inherited unit root')
        denominators[name]=U.reciprocal_bases(args['constructor'],walk)
        native_units[name]={'constructor':args['constructor'],'unit':ret['unit'],'reciprocalBases':denominators[name]}
        J.finish(native_units[name])
    J.emit('unbound-reciprocal-homogeneity',{'operations':denominators,'normalizedBoundDenominatorAssignedUnit':False,'argument':'Homogeneous native reciprocals plus homogeneous substitutions transport a rational source. Cancellation may rescale numerator and denominator together; no unique unit is assigned to the saved numeric normalization.'})
    groups={};source_origins={}
    for face in ('plus','minus'):
        inp=get('original/'+face+'-source-input.json');old=get('original/'+face+'-source-grade-split.json');chem=get('original/native-chemical-amplitude.json')
        ops=sc(face+'-native-source-join-operands');bound=sc(face+'-native-source-join-bound')
        J.emit(face+'-transport-source-operands',{'nativeUnits':native_units[face+'-unbound-source-units'],'original':inp,'chemical':chem,'sourceJoin':ops,'boundJoin':bound,'savedSource':old})
        require(ops['saved']==inp and ops['nativeChemical']==chem and ops['epsilon']==context['epsilon'],'same actual stage2 operands')
        require(U.equal_text(native_units[face+'-unbound-source-units']['constructor'],inp['raw']['srepr']),'same unbound normalized source')
        require(U.equal_text(native_units[face+'-native-velocity-units']['constructor'],inp['nativeVelocity']['srepr']),'same native velocity')
        require(U.equal_text(native_units['new-unbound-chemical-homogeneity']['constructor'],U.tuple_leaf(chem['raw']['srepr'],1)),'same native chemical expression')
        stage={U.symbol_name(U.parse(a['srepr'])):b for a,b in ops['stage2Map']}
        require(stage=={'s11cc1_V_lab_held_'+face:inp['velocityAmplitude'],'s11cc1_mu_theta_lab_held_'+face:chem['amplitude']},'actual simultaneous stage2 map')
        require(bound['target']==inp['combined'],'same source target')
        prior(face+'-native-chemical-raw-join',right=U.encoded_leaf(chem['raw'],1))
        prior(face+'-native-density-map-join',left=context['densityMap']['rho_br_bg_rho4_constant'])
        prior(face+'-native-chemical-epsilon-join',right=chem['amplitude'])
        prior(face+'-own-velocity-normalization',right=inp['velocityAmplitude'])
        prior(face+'-raw-stage2-live-density-join',left=bound['bound'],right=inp['combined'])
        prior(face+'-native-composed-source',left=inp['combined'],right=old['full'])
        # Original normalization and chemical domain are complete inherited evidence.
        prior(face+'-inherited-velocity-normalization-join',left=inp['velocityCoefficient'],right=get('original/inherited-source-normalization.json')['velocityCoefficient'])
        domain=get('original/chemical-amplitude-domain.json')
        prior(face+'-chemical-domain-original-join',left=domain['original'],right=chem['amplitude'])
        prior(face+'-chemical-domain-reduced-join',left=domain['reduced'],right=chem['amplitude'])
        prior(face+'-chemical-domain-fraction-join',left=domain['numerator'])
        prior(face+'-chemical-domain-zero-grade-join',right=domain['denominatorAtZero'])
        U.regularity(raw,'composition/'+face+'-chemical-domain-denominator',domain['denominatorAtZero'],J)
        source_origins[face]=old['full']
    consumer_units=get('native/inherited-consumer-unit-returns.json')
    consumer_input=get('original/THETA_BALANCE-consumer-input.json')
    require(len(consumer_units)==4 and len({v['slot']['text'] for v in consumer_units})==4,'all four inherited consumer slots')
    prior('THETA_BALANCE-pressure-bound-join',right=consumer_input['bound'])
    prior('THETA_BALANCE-pressure-source-join',right=consumer_input['raw'])
    prior('THETA_BALANCE-affine-pressure-reconstruction',left=consumer_input['bound'])
    for name in m['gradeGroups']:
        operands=sc(name+'-operands');split=sc(name+'-split')
        J.start('transport-grade-'+name,{'operands':operands,'split':split,'dimensionlessGrades':context['independentGrades']})
        require(operands['grades']==context['independentGrades'],'actual independent grade arguments')
        if name.endswith('-source'):
            face=name.split('-')[0];unit=U.unit_tuple(native_units[face+'-unbound-source-units']['unit']);origin=source_origins[face]
        else:
            slot=name.removeprefix('THETA_BALANCE-');item=next(u for u in consumer_units if u['slot']['text']==slot)
            unit=U.matrix_unit(item['coefficientDimension']);origin=consumer_input['slotCoefficients'][slot]
            # Save full actual unbound coefficient. Its unit proof is inherited, not recalculated.
            J.emit('consumer-origin-'+slot,item)
            prior(name+'-native-unit-coefficient',right=origin)
            prior(name+'-full-coefficient',left=origin)
        require(U.equal(operands['full'],origin),'same full rational grade source '+name)
        U.regularity(raw,'composition/'+name+'-regular-denominator',split['denominatorAtZero'],J)
        prior(name+'-saved-zero',left=split['retained']['(0, 0)'],right=operands['savedZero'])
        for g in ('00','10','01','11'):prior(name+'-quotient-'+g,right=ZERO)
        require(set(split['retained'])=={'(0, 0)','(1, 0)','(0, 1)','(1, 1)'},'full independent rectangle')
        require(split['nativeShapeCoefficientsCalled'] is False and split['zeroGradeReused'] is True,'saved quotient origin')
        # Assessed homogeneity lemma: regular Taylor coefficients in dimensionless eta/sigma retain the source unit.
        groups[name]={'unit':unit,'split':split,'homogeneityArgument':'Native homogeneous source + unit-preserving binding and dimensionless independent regular quotient; saved recurrence/remainder identities inherited. Not a dimensional proof from numeric polynomial.'}
        J.finish({'sourceUnit':unit,'retainedUnits':{g:None if U.zero(v) else unit for g,v in split['retained'].items()},'zeroRequiredUnit':unit,'normalizedBoundDenominatorUnit':None,'homogeneousUnboundOriginJoined':True,'regularityInherited':True})
    fields=get('inventory/transforms.json');selected=get('inventory/selected.json')['selected'];index=get('transport-index.json')
    require(len(selected)==544 and len({a['addressId'] for a in selected})==544,'complete selected address census')
    source_records={};consumer_records={}
    for face in ('plus','minus'):
        for code in ('00','10','01','11'):
            table=sc(face+'-source-jets-'+code);grade=str(tuple(map(int,code)))
            require(U.equal(table['source'],groups[face+'-source']['split']['retained'][grade]),'actual jet grade source')
            prior(face+'-source-jet-reconstruction-'+code,left=table['source'],right=table['reconstruction'])
    for key,loc in index['sourceLocations'].items():
        face,code,position=key.split('/');table=sc(face+'-source-jets-'+code);item=table['jets'][int(position)]
        candidates=[get(m['indexedAliases'][path]) for path in loc['profileCandidateRecords']]
        J.start('source-field-'+key.replace('/','-'),{'location':loc,'item':item,'profileCandidates':candidates,'sourceGroup':face+'-source','grade':code})
        require(item['spec']['channel']=='e_W','selected actual scalar source channel')
        jet=U.jet_unit(item['spec'],registry,item['atom']);required=U.add(groups[face+'-source']['unit'],U.scale(jet,-1))
        require(tuple(item['jetDimension'])==jet and tuple(item['requiredCoefficientDimension'])==required,'original jet unit metadata joins, not its proof')
        maps=[]
        for p in candidates:
            require(U.equal(p['original'],item['coefficient']) and U.equal(p['oneDimensional'],item['field']),'actual profile operands')
            maps.append(U.profile_transport(p,context,registry))
        require(maps,'at least one exact profile map')
        source_records[key]={'requiredUnit':required,'transportedUnit':None if U.zero(item['field']) else required,'zero':U.zero(item['field']),'jetUnit':jet,'item':item,'maps':maps,'proof':'Homogeneity and actual inherited linear independent-jet reconstruction; no coefficient re-extraction.'}
        J.finish(source_records[key])
    for key,loc in index['consumerLocations'].items():
        slot,grade=key.split('/');g=tuple(int(v.strip()) for v in grade.strip('()').split(','));code=''.join(map(str,g));name='THETA_BALANCE-'+slot
        original=groups[name]['split']['retained'][grade];epsproof=prior(name+'-epsilon-'+code,left=original)
        stripped=U.strip_epsilon(epsproof['right'],context['epsilon'])
        candidates=[get(m['indexedAliases'][path]) for path in loc['profileOutputCandidateRecords']]
        matching=[p for p in candidates if U.equal_text(p['original']['srepr'],stripped)]
        J.start('consumer-field-'+slot+'-'+code,{'location':loc,'original':original,'epsilonProof':epsproof,'strippedConstructor':stripped,'profileCandidates':candidates,'matchingCount':len(matching)})
        require(matching,'actual epsilon-stripped consumer profile input')
        maps=[U.profile_transport(p,context,registry) for p in matching];field=matching[0]['oneDimensional']
        require(all(U.equal(p['oneDimensional'],field) for p in matching),'identical consumer field for all input matches')
        unit=groups[name]['unit'];consumer_records[key]={'requiredUnit':unit,'transportedUnit':None if U.zero(field) else unit,'zero':U.zero(field),'field':field,'maps':maps,'original':original}
        J.finish(consumer_records[key])
    addresses=[]
    for loc in index['addressLocations']:
        a=selected[int(loc['selectedJsonPointer'].rsplit('/',1)[1])];s=source_records[loc['source']];c=consumer_records[loc['consumer']]
        J.start('address-'+str(a['addressId']),{'address':a,'sourceLocation':loc['source'],'consumerLocation':loc['consumer'],'fieldIds':loc['fieldIds']})
        require(a['addressId']==loc['addressId'] and a['row']=='THETA_BALANCE','actual address identity')
        require(a['sourceAtom']==s['item']['atom'] and a['jet']==s['item']['spec'] and U.equal(a['sourceOriginal'],s['item']['coefficient']) and U.equal(a['sourceField'],s['item']['field']),'source argument identity')
        require(U.equal(a['consumerOriginal'],c['original']) and U.equal(a['consumerField'],c['field']),'consumer argument identity')
        for role,field_id in zip(('source','consumer'),loc['fieldIds']):
            require(a[role+'Transform']['coefficientId']==field_id and U.equal(fields[field_id]['field'],a[role+'Field']),'whole coefficient field join')
            require(hashlib.sha256(a[role+'Field']['srepr'].encode()).hexdigest()==field_id,'field ID is extra identity, not dimensional proof')
            poly=get('fields/'+field_id+'-polynomial.json');args=get('fields/'+field_id+'-reconstruction-input.json');ret=get('fields/'+field_id+'-reconstruction-return.json')
            U.field_proof(poly,args,ret,a[role+'Field'])
        require(a['epsilon']==context['epsilon'] and a['epsilonCount']==(0 if a['status'].startswith('EXACT_ZERO') else 1),'actual epsilon ancestry')
        result={'addressId':a['addressId'],'status':a['status'],'sourceUnit':s['transportedUnit'],'sourceRequiredUnit':s['requiredUnit'],'consumerUnit':c['transportedUnit'],'consumerRequiredUnit':c['requiredUnit'],'nativeSourceAndSavedGradeProfileJoined':True,'numericPolynomialIndependentlyDimensional':False,'waveKernelMeasureUnitsChecked':False}
        addresses.append(result);J.finish(result)
    # New responsive unit transport controls use actual selected profile maps and live jet entries.
    controls=[]
    candidates=[(k,s) for k,s in source_records.items() if not s['zero'] and any(v['derivativeOrder']>0 and not v['transverseZero'] for mp in s['maps'] for v in mp['entries'])]
    require(candidates,'applicable live source profile derivative')
    key,s=candidates[0];entry=next(v for mp in s['maps'] for v in mp['entries'] if v['derivativeOrder']>0 and not v['transverseZero'])
    label='missing-profile-L';mutant={'LExponent':0}
    J.start('control-'+label,{'sourceLocation':key,'actualItem':s['item'],'actualProfileEntry':entry,'mutation':mutant})
    expected=U.ZERO;measured=U.add(U.scale(U.unit_tuple(registry['L_W']),0),(-entry['derivativeOrder'],0,0))
    result={'control':label,'baselineUnit':expected,'mutatedUnit':measured,'refused':measured!=expected,'scope':'Actual saved profile derivative and native scale rule, not a field value'}
    J.emit('control-'+label+'-decision-operands',result);require(result['refused'],'responsive missing native scale');J.finish(result);controls.append(result)
    # Mutate the actual original chemical unit calculation, not a displayed required-unit label.
    op=get('native/new-unbound-chemical-homogeneity-input.json');mutated_registry=dict(op['registry']);mutated_registry['eta_bg']=[1,0,0]
    J.start('control-grade-assigned-length',{'actualInheritedInput':op,'baselineReturn':get('native/new-unbound-chemical-homogeneity-return.json'),'mutatedRegistry':mutated_registry})
    checker=U.Units(mutated_registry);refused=False;reason=None
    try:mutated=checker.dimension(U.parse(op['constructor']))
    except ValueError as error:refused=True;reason=str(error);mutated=None
    finally:J.emit('control-grade-assigned-length-unit-walk',checker.events)
    result={'control':'grade-assigned-length','refused':refused,'reason':reason,'mutatedUnit':mutated,'scope':'Actual original chemical expression under changed grade units; no field calculation.'}
    J.emit('control-grade-assigned-length-decision-operands',result);require(refused and 'inhomogeneous' in reason,'actual native grade-unit refusal');J.finish(result);controls.append(result)
    result={'status':'BOUNDED_SOURCE_GRADE_PROFILE_UNIT_TRANSPORT_COMPLETE','selectedAddresses':len(addresses),'sourceLocations':len(source_records),'consumerLocations':len(consumer_records),'gradeGroups':len(groups),'controls':controls,'nativeInferenceDependency':True,'sourceGradeFieldUnitTransportComplete':True,'pressureSummandUnitsComplete':False,'waveKernelMeasureUnitsComplete':False,'numericalEvaluatorReady':False,'newIntegralOrAction':False,'priorFunctionsReplayed':False,'limits':['Homogeneity transported through original accepted algebra, not new independent source/grade calculation.','Zero coefficients carry required-unit expectations only.','Gamma registry units remain inference-dependent.','Wave/profile-transform/kernel/measure/all-summand proof remains required.']}
    J.emit('complete-address-unit-transports',addresses);return result


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
