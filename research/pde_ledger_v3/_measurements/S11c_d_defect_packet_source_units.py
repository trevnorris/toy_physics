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
    def assembly(name,expected,actual):
        J.emit(name+'-assembly-input',{'constructedFromActualSavedOperands':ast.unparse(expected),'savedOperand':actual,'newUnitProvenanceJoin':True,'oldProducerOrQuotientRecurrenceCalled':False})
        left=U.ring(expected);right=U.ring(U.parse(actual['srepr']))
        result={'leftNormalForm':U.dump_ring(left),'rightNormalForm':U.dump_ring(right),'exactlyEqual':left==right}
        J.emit(name+'-assembly-decision',result);require(left==right,'new exact saved-operand assembly '+name)
        return result
    def rational_join(name,left,right):
        J.emit(name+'-rational-input',{'left':ast.unparse(left),'right':ast.unparse(right),'newIdentityOnly':True})
        try:decision=U.rational_identity(left,right)
        except BaseException as error:
            J.emit(name+'-rational-refusal',{'errorType':type(error).__name__,'error':str(error)});raise
        J.emit(name+'-rational-decision',decision);require(decision['matches'],'actual rational source identity '+name)
        for i,item in enumerate(decision['originalInverseDomains']):
            certificate=U.inverse_origin_certificate(item['inverseBase'],context['independentGrades'])
            J.emit(name+'-original-domain-'+str(i),certificate)
            require(certificate['finiteNonzero'],'original rational factor finite/nonzero at grade origin')
        return decision
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
        J.emit('numeric-binding-'+name+'-input',{'symbol':name,'value':value,'unit':registry.get(name),'originalPhysicalParameters':context['physicalInput']['parameters'],'heldFrequency':context['frequency']})
        decision=U.numeric_binding(name,value,context,registry);J.emit('numeric-binding-'+name+'-return',decision);numeric.append(decision)
    J.emit('numeric-binding-provenance',numeric)
    scale=sc('native-profile-scale-join');require(scale['savedLength']==context['numeric']['L_W'] and scale['physicalLength']==context['physicalInput']['parameters']['L_W'],'actual original length join')
    U.verify_contracts(m,get('source-contracts.json'))
    J.emit('actual-source-contracts',get('source-contracts.json'))
    definitions=get('source-contracts.json')['profileDefinitions'];text=Path(definitions['source']).read_text()
    for item in definitions['statements']:
        matches=[n for n in ast.walk(ast.parse(text)) if isinstance(n,ast.Assign) and n.lineno==item['line'] and n.col_offset==item['column']]
        require(len(matches)==1 and ast.get_source_segment(text,matches[0])==item['text'],'actual native profile base definition')
    require(U.same(ast.parse(definitions['statements'][0]['text']).body[0],ast.parse('w=(1+sp.tanh(x/profile_length))/2').body[0]) and U.same(ast.parse(definitions['statements'][1]['text']).body[0],ast.parse('m=(1-sp.tanh(x/profile_length)**2)/3').body[0]),'actual profile base coefficients')
    J.emit('profile-native-base-definitions',definitions)
    all_profile_records={name:value for name,value in raw.items() if name.startswith('composition/profile-map-')}
    catalog,catalog_origins=U.scale_catalog(all_profile_records,context);J.emit('profile-scale-catalog',{'entries':[{'base':k[0],'order':k[1],'mapped':v,'origins':catalog_origins[k]} for k,v in catalog.items()]})
    verifier=U.ScaleVerifier(catalog,context,J)
    for (base,order),mapped in sorted(catalog.items()):verifier.verify(base,order,mapped)

    # Preserve complete unit walks; select reciprocal bases by their actual AST paths.
    native_units={};denominators={}
    for name in m['inheritedUnitOperations']:
        args=get('native/'+name+'-input.json');ret=get('native/'+name+'-return.json');walk=get('native/'+name+'-unit-walk.json')
        J.start('restore-'+name,{'input':args,'return':ret,'walk':walk})
        require(args['registry']==registry,'actual inherited unit registry '+name)
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
        stage_units=[]
        for name,op in [('s11cc1_V_lab_held_'+face,face+'-native-velocity-units'),('s11cc1_mu_theta_lab_held_'+face,'new-unbound-chemical-homogeneity')]:
            item=U.stage_unit_join(name,stage[name],registry,native_units[op]['unit']);stage_units.append(item)
        J.emit(face+'-actual-stage2-unit-joins',stage_units)
        require(all(item['matches'] for item in stage_units),'same stage2 placeholder and replacement units')
        require(bound['target']==inp['combined'],'same source target')
        chemical_leaf=U.encoded_leaf(chem['raw'],1);density_leaf=U.encoded_leaf(context['density'],1)
        prior(face+'-native-chemical-raw-join',left={'srepr':native_units['new-unbound-chemical-homogeneity']['constructor']},right=chemical_leaf)
        prior(face+'-native-density-raw-join',left={'srepr':native_units['new-unbound-live-density-units']['constructor']},right=density_leaf)
        prior(face+'-native-chemical-amplitude-join',left=inp['chemicalAmplitude'],right=chem['amplitude'])
        # New missing unit-provenance joins on actual original operands. These
        # do not call the old source builder, bind, grade or profile functions.
        chemical_expected,chemical_steps=U.binding_operand(ast.unparse(U.times(U.parse(chemical_leaf['srepr']),U.expr_call('Pow',U.parse(context['epsilon']['srepr']),U.integer_node(-1)))),context)
        epsilon_args=sc(face+'-native-chemical-epsilon-join-input')
        J.emit(face+'-chemical-unit-binding-stages',{'stages':chemical_steps,'originalPhysicalUnit':native_units['new-unbound-chemical-homogeneity']['unit'],'dimensionlessEpsilon':context['epsilon']})
        assembly(face+'-chemical-epsilon-left',chemical_expected,epsilon_args['left'])
        velocity_expected,velocity_steps=U.binding_operand(ast.unparse(U.times(U.parse(inp['nativeVelocity']['srepr']),U.expr_call('Pow',U.parse(context['epsilon']['srepr']),U.integer_node(-1)))),context)
        velocity_args=sc(face+'-own-velocity-normalization-input');J.emit(face+'-velocity-unit-binding-stages',velocity_steps)
        assembly(face+'-velocity-normalization-left',velocity_expected,velocity_args['left'])
        identified=U.constructor_substitute(U.parse(inp['raw']['srepr']),stage)
        assembly(face+'-stage2-identified',identified,bound['identified'])
        full_expected,full_steps=U.binding_operand(ast.unparse(identified),context);J.emit(face+'-full-source-unit-binding-stages',full_steps)
        assembly(face+'-full-bound-left',full_expected,bound['bound'])
        prior(face+'-native-density-map-join',left=context['densityMap']['rho_br_bg_rho4_constant'])
        prior(face+'-native-chemical-epsilon-join',left=epsilon_args['left'],right=chem['amplitude'])
        prior(face+'-own-velocity-normalization',left=velocity_args['left'],right=inp['velocityAmplitude'])
        prior(face+'-raw-stage2-live-density-join',left=bound['bound'],right=inp['combined'])
        prior(face+'-native-composed-source',left=inp['combined'],right=old['full'])
        # Original normalization and chemical domain are complete inherited evidence.
        prior(face+'-inherited-velocity-normalization-join',left=inp['velocityCoefficient'],right=get('original/inherited-source-normalization.json')['velocityCoefficient'])
        domain=get('original/chemical-amplitude-domain.json')
        prior(face+'-chemical-domain-original-join',left=domain['original'],right=chem['amplitude'])
        prior(face+'-chemical-domain-reduced-join',left=domain['reduced'],right=chem['amplitude'])
        prior(face+'-chemical-domain-fraction-join',left=domain['numerator'])
        domain_args=sc(face+'-chemical-domain-zero-grade-join-input')
        assembly(face+'-chemical-domain-origin',U.at_grade_zero(U.parse(domain['denominator']['srepr']),context['independentGrades']),domain_args['left'])
        prior(face+'-chemical-domain-zero-grade-join',left=domain_args['left'],right=domain['denominatorAtZero'])
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
            unit_args=sc(name+'-native-unit-coefficient-input')
            unit_left,steps=U.binding_operand(item['coefficient']['srepr'],context);J.emit(name+'-native-consumer-unit-binding-stages',steps)
            rational_join(name+'-native-consumer-unit-left',unit_left,U.parse(unit_args['left']['srepr']))
            prior(name+'-native-unit-coefficient',left=unit_args['left'],right=origin)
            prior(name+'-full-coefficient',left=origin)
        require(U.equal(operands['full'],origin),'same full rational grade source '+name)
        assembly(name+'-denominator-origin',U.at_grade_zero(U.parse(split['denominator']['srepr']),context['independentGrades']),split['denominatorAtZero'])
        U.regularity(raw,'composition/'+name+'-regular-denominator',split['denominatorAtZero'],J)
        rational_join(name+'-full-fraction',U.parse(operands['full']['srepr']),U.times(U.parse(split['numerator']['srepr']),U.expr_call('Pow',U.parse(split['denominator']['srepr']),U.integer_node(-1))))
        prior(name+'-saved-zero',left=split['retained']['(0, 0)'],right=operands['savedZero'])
        for field,expected in U.full_remainder_operands(operands['full'],split,context['independentGrades']).items():assembly(name+'-'+field,expected,split[field])
        J.emit(name+'-full-quotient-proof-provenance',{'numerator':split['numerator'],'denominator':split['denominator'],'retained':split['retained'],'excludedPure':split['excludedPure'],'fullHigherRemainder':split['fullHigherRemainder'],'quotientRingNumeratorRemainder':split['quotientRingNumeratorRemainder'],'completedGradeSource':m['compositionSource'],'acceptedSourceHash':composition_result['workerSha256'],'quotientRecurrenceRerun':False,'note':'The following literal zero observations are inherited outputs of this exact saved grade operation, not independent 0=0 proofs of the full remainder.'})
        require(sha(m['compositionSource'])==composition_result['workerSha256'],'same accepted grade-operation source')
        for g in ('00','10','01','11'):prior(name+'-quotient-'+g,left=ZERO,right=ZERO)
        require(set(split['retained'])=={'(0, 0)','(1, 0)','(0, 1)','(1, 1)'},'full independent rectangle')
        require(split['nativeShapeCoefficientsCalled'] is False and split['zeroGradeReused'] is True,'saved quotient origin')
        # Assessed homogeneity lemma: regular Taylor coefficients in dimensionless eta/sigma retain the source unit.
        groups[name]={'unit':unit,'split':split,'homogeneityArgument':'Native homogeneous source + unit-preserving binding and dimensionless independent regular quotient; saved recurrence/remainder identities inherited. Not a dimensional proof from numeric polynomial.'}
        J.finish({'sourceUnit':unit,'retainedUnits':{g:None if U.zero(v) else unit for g,v in split['retained'].items()},'zeroRequiredUnit':unit,'normalizedBoundDenominatorUnit':None,'homogeneousUnboundOriginJoined':True,'regularityInherited':True,'componentIdentitiesNew':True,'rationalSourceAndOriginJoined':True})
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
            maps.append(U.profile_transport(p,context,registry,verifier))
        require(bool(maps),'at least one exact profile map')
        source_records[key]={'requiredUnit':required,'transportedUnit':None if U.zero(item['field']) else required,'quantitySemantics':{'savedNumber':'numeric magnitude in original fixed base units','physicalCoefficientUnit':None if U.zero(item['field']) else required,'jetPhysicalUnit':jet,'physicalSourceUnit':groups[face+'-source']['unit'],'unitReferences':['L_ref','T_ref','M_ref'],'bareNumericExpressionPassedToUnitEngine':False},'zero':U.zero(item['field']),'jetUnit':jet,'item':item,'maps':maps,'proof':'Homogeneity and actual inherited linear independent-jet reconstruction; no coefficient re-extraction.'}
        J.finish(source_records[key])
    for key,loc in index['consumerLocations'].items():
        slot,grade=key.split('/');g=tuple(int(v.strip()) for v in grade.strip('()').split(','));code=''.join(map(str,g));name='THETA_BALANCE-'+slot
        original=groups[name]['split']['retained'][grade];epsproof=prior(name+'-epsilon-'+code,left=original)
        stripped=U.strip_epsilon(epsproof['right'],context['epsilon'])
        candidates=[get(m['indexedAliases'][path]) for path in loc['profileOutputCandidateRecords']]
        matching=[p for p in candidates if U.equal_text(p['original']['srepr'],stripped)]
        J.start('consumer-field-'+slot+'-'+code,{'location':loc,'original':original,'epsilonProof':epsproof,'strippedConstructor':stripped,'profileCandidates':candidates,'matchingCount':len(matching)})
        require(bool(matching),'actual epsilon-stripped consumer profile input')
        maps=[U.profile_transport(p,context,registry,verifier) for p in matching];field=matching[0]['oneDimensional']
        require(all(U.equal(p['oneDimensional'],field) for p in matching),'identical consumer field for all input matches')
        unit=groups[name]['unit'];consumer_records[key]={'requiredUnit':unit,'transportedUnit':None if U.zero(field) else unit,'quantitySemantics':{'savedNumber':'numeric magnitude in original fixed base units','physicalCoefficientUnit':None if U.zero(field) else unit,'unitReferences':['L_ref','T_ref','M_ref'],'bareNumericExpressionPassedToUnitEngine':False},'zero':U.zero(field),'field':field,'maps':maps,'original':original}
        J.finish(consumer_records[key])
    require(len(source_records)==64 and len(consumer_records)==16,'all original source/consumer locations')
    addresses=[];pointers=[v['selectedJsonPointer'] for v in index['addressLocations']]
    J.emit('selected-address-coverage-input',{'selectedIds':[v['addressId'] for v in selected],'locations':index['addressLocations']})
    require(len(pointers)==544 and len(set(pointers))==544 and set(pointers)=={'/selected/'+str(i) for i in range(544)},'exact complete selected pointer coverage')
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
    require(bool(candidates),'applicable live source profile derivative')
    key,s=candidates[0];entry=next(v for mp in s['maps'] for v in mp['entries'] if v['derivativeOrder']>0 and not v['transverseZero'])
    actual_record=next(p for p in all_profile_records.values() if any(v[0]==entry['original'] and v[1]==entry['savedMapped'] for v in p['map']))
    mutant=json.loads(json.dumps(actual_record));position=next(i for i,v in enumerate(mutant['map']) if v[0]==entry['original'])
    original_mapped=mutant['map'][position][1]
    changed=U.times(U.parse(original_mapped['srepr']),U.expr_call('Pow',U.parse(context['numeric']['L_W']['srepr']),U.integer_node(-entry['derivativeOrder'])))
    mutant['map'][position][1]={'text':ast.unparse(changed),'srepr':ast.unparse(changed)}
    J.start('control-missing-profile-L',{'sourceLocation':key,'actualItem':s['item'],'actualProfileRecord':actual_record,'mutatedRecord':mutant,'mutation':'Divide the actual saved mapped derivative by its original L^r; call the same profile transport path.'})
    refused=False;reason=None
    try:U.profile_transport(mutant,context,registry,verifier)
    except ValueError as error:refused=True;reason=str(error)
    result={'control':'missing-profile-L','refused':refused,'reason':reason,'actualProfileTransportCalled':True,'scope':'Actual mapped-profile scale certificate, not a numerical field/action value.'}
    J.emit('control-missing-profile-L-decision-operands',result);require(refused and 'actual mapped profile' in reason,'responsive actual profile scale refusal');J.finish(result);controls.append(result)
    # Mutate the actual original chemical unit calculation, not a displayed required-unit label.
    op=get('native/new-unbound-chemical-homogeneity-input.json');mutated_registry=dict(op['registry']);mutated_registry['eta_bg']=[1,0,0]
    J.start('control-grade-assigned-length',{'actualInheritedInput':op,'baselineReturn':get('native/new-unbound-chemical-homogeneity-return.json'),'mutatedRegistry':mutated_registry})
    checker=U.Units(mutated_registry);refused=False;reason=None
    try:mutated=checker.dimension(U.parse(op['constructor']))
    except ValueError as error:refused=True;reason=str(error);mutated=None
    finally:J.emit('control-grade-assigned-length-unit-walk',checker.events)
    result={'control':'grade-assigned-length','refused':refused,'reason':reason,'mutatedUnit':mutated,'scope':'Actual original chemical expression under changed grade units; no field calculation.'}
    J.emit('control-grade-assigned-length-decision-operands',result);require(refused and 'inhomogeneous' in reason,'actual native grade-unit refusal');J.finish(result);controls.append(result)
    J.emit('selected-address-coverage-return',{'joinedIds':[v['addressId'] for v in addresses],'selectedIds':[v['addressId'] for v in selected]})
    require(len(addresses)==544 and len({v['addressId'] for v in addresses})==544 and {v['addressId'] for v in addresses}=={v['addressId'] for v in selected},'every selected address exactly once')
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
