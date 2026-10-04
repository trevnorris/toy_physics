#!/usr/bin/env python3
"""Guarded native pressure-source/closure unit bridge; complete summand proof pending.

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
    # the geometry library separately refuses every floating arithmetic input.
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
    require(g['status']=='READY_FOR_ONE_PACKET_NATIVE_UNITS' and g['independentBuildClearance'] is True,'actual independent build readiness')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==m['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','library','authority','buildReviewRecord'):
        require(sha(g[key])==g[key+'Sha256'],'gate pin '+key)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual guard/supervisor')
    require(g['launcher']==m['launcher'] and g['library']==m['librarySource'] and g['authority']==m['executionAuthority'] and g['buildReviewRecord']==m['reviewRecordWillBe'],'actual document paths')
    r=read(g['buildReviewRecord']);require(r['allChecksPassed'] is True and r['independentBuildClearance'] is True,'build assessment')
    for key in ('workerSha256','manifestSha256','librarySha256','launcherSha256','sharedGuardSha256','supervisorSha256'):require(r[key]==g[key],'actual reviewed '+key)
    require(all(r['reports'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED NATIVE PRESSURE-UNIT BUILD' for e in ('claude','grok')),'both literal build reports')
    method=read(m['methodRecord']);require(g['methodRecordSha256']==sha(m['methodRecord']) and method['jointIndependentMethodClearance'] and method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'actual pressure-readiness method')
    a=read(g['authority']);require(a['scope']==g['scope']==m['scope'] and a['boundedInstrumentAuthorized'] and a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and a['automaticScientificRetry'] is False and a['noDeadline'] and g['durationLimits'] is None,'standing bounded authority')
    return g


def load_library(m):
    spec=importlib.util.spec_from_file_location('packet_native_unit_interpreter',m['librarySource']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


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
    native=raw['native/source.json'];merged=raw['units/merged.json'];context=raw['consumer/binding-context.json']
    J.emit('original-input-context',{'context':context,'physical':raw['physical-input.json'],'scope':m['scope']})
    require(context['physicalInput']==raw['physical-input.json'] and context['frequency']=={'text':'3','srepr':'Integer(3)'},'same original input and held omega3')
    # Exact native line and quoted constructor identities, without restoring them.
    def provenance(record):
        source=Path(record['source']);lines=source.read_bytes().splitlines(keepends=True);line=lines[record['valueLine']-1]
        require(hashlib.sha256(line).hexdigest()==record['sourceLineSha256'] and len(line)==record['sourceLineBytes'],'original native source line')
        text=U.quoted_restore_line(line.decode())
        require(hashlib.sha256(text.encode()).hexdigest()==record['constructorSha256'] and len(text.encode())==record['constructorBytes'],'complete original constructor')
        return text
    sources={}
    for record in [native['chemicalSource']['provenance'],native['faceResponseSources']['provenance'],native['geometry']['face_velocity']['bLiteralSource'],native['geometry']['background_density_map']['bLiteralSource']]:
        key=record['source']+':'+str(record['valueLine']);J.emit('native-provenance-'+str(len(sources)),record);sources[key]=provenance(record)
    # Native schema is literal metadata. Do not invoke Inputs or infer_dimensions.
    schema_source=Path(m['nativeC2']).read_text();tree=ast.parse(schema_source)
    schema_nodes=[n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='DIMENSION_SCHEMA' for t in n.targets)]
    require(len(schema_nodes)==1,'unique original static unit schema');schema=ast.literal_eval(schema_nodes[0].value)
    J.emit('native-static-schema',{'source':m['nativeC2'],'sourceSha256':sha(m['nativeC2']),'statement':ast.get_source_segment(schema_source,schema_nodes[0]),'schema':schema})
    registry=dict(schema);unavailable={};origins={}
    for label in ['LEFTNativeSource','RIGHTNativeSource','LEFTPairing','RIGHTPairing']:
        o=raw['units/'+label+'-origin.json'];receipt=o['objectReceipt'];J.emit('registry-origin-'+label,o)
        require(sha(receipt['path'])==receipt['sha256'] and Path(receipt['path']).stat().st_size==receipt['bytes'],'original opaque registry object receipt')
        require(o['sourceReexecuted'] is False and o['registryEntries']>0,'inherited registry origin');origins[label]=o
    for name,records in merged['generatedEntries'].items():
        J.emit('registry-'+name,records)
        require({r['registry'] for r in records}==set(origins) and len(records)==4,'four original registries')
        require(all(r['symbol']['text']==name for r in records),'native generated symbol label')
        try:units=[U.unit_tuple(r['unit']) for r in records]
        except (ValueError,TypeError,ZeroDivisionError,KeyError) as exc:
            unavailable[name]={'records':records,'reason':str(exc)};registry.pop(name,None);continue
        require(all(v==units[0] for v in units),'native registry disagreement')
        require(name not in registry or U.unit_tuple(registry[name])==units[0],'static/generated unit conflict');registry[name]=units[0]
    J.emit('effective-registry',{'units':registry,'unavailable':unavailable,'newInference':False,'independentOfOriginalInference':False})
    for name,expected in [('rho_m',[-4,0,1]),('rho_br',[-3,0,1]),('epsilon_shape',[0,0,0]),('eta_bg',[0,0,0]),('sigma_W',[0,0,0])]:require(U.unit_tuple(registry[name])==tuple(expected),'actual native unit '+name)
    # Chemical constructor is the original tuple's second leaf, not a rebuilt source.
    chemical=native['chemicalSource'];require(chemical['case']==['LAB_HELD','RHO4_CONSTANT'],'declared chemical case')
    chemical_value=join_native_case(J,U,'chemical-original-case',sources[chemical['provenance']['source']+':'+str(chemical['provenance']['valueLine'])],chemical)
    chemical_node=U.tuple_args(chemical_value)[1]
    chemical_old=raw['consumer/native-chemical-amplitude.json']
    require(U.same(chemical_node,U.tuple_args(U.parse(chemical_old['raw']['srepr']))[1]),'actual original chemical operand')
    certify(J,U,'new-unbound-chemical-homogeneity',ast.unparse(chemical_node),registry,[-1,-2,1])
    density=next(v for v in native['geometry']['background_density_map']['cases'] if v['case']==['RHO4_CONSTANT'])
    density_value=join_native_case(J,U,'density-original-case',sources[density['provenance']['source']+':'+str(density['provenance']['valueLine'])],density)
    density_leaf=U.tuple_args(density_value)[1]
    J.emit('native-live-density-arguments',{'case':density,'selectedLeaf':ast.unparse(density_leaf),'savedMap':context['densityMap']})
    require(U.same(density_leaf,U.parse(context['densityMap']['rho_br_bg_rho4_constant']['srepr'])),'actual density map leaf')
    certify(J,U,'new-unbound-live-density-units',ast.unparse(density_leaf),registry,[-3,0,1])
    flat=raw['native/flat-join.json'];require(flat['residual']==ZERO and flat['left']==flat['right'],'inherited native flat identity')
    dtn=raw['native/dtn-operands.json'];dtn_tree=U.parse(dtn['literal'])
    # All native cases supplied; select by literal labels, never a fitted mapping.
    dtn_cases=U.tuple_args(dtn_tree)
    kernel_source=Path(m['nativeC1Exports']).read_text();kernel_ast=ast.parse(kernel_source)
    J.start('original-dtn-export-key',{'source':m['nativeC1Exports'],'sourceSha256':sha(m['nativeC1Exports']),'key':'dtn_kernel','savedLiteral':dtn['literal']})
    keyed_literal,keyed_call=U.export_restore_literal(kernel_ast,'dtn_kernel')
    inputs_class=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='Inputs')
    constructor=next(n for n in inputs_class.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    lookup=[n for n in constructor.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Attribute) and isinstance(t.value,ast.Name) and t.value.id=='self' and t.attr=='kernel' for t in n.targets)]
    expected_lookup=ast.parse("self.kernel = cases(values['dtn_kernel'])").body[0]
    keyed_result={'selectedRestoreCall':ast.unparse(keyed_call),'lookupStatements':[ast.unparse(n) for n in lookup],
        'expectedLookup':ast.unparse(expected_lookup),'literalMatches':keyed_literal==dtn['literal'],
        'lookupMatches':len(lookup)==1 and U.same(lookup[0],expected_lookup)}
    J.emit('original-dtn-export-key-decision-operands',keyed_result)
    require(keyed_result['literalMatches'] and keyed_result['lookupMatches'],'actual original DTN export key and consumer lookup');J.finish(keyed_result)
    contract=raw['source-contracts.json']['kernel_bridge'];native_bridge=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='kernel_bridge')
    require(U.signature(ast.parse(contract).body[0])==U.signature(native_bridge),'actual kernel bridge source contract')
    required_statements=["diagonal = named(raw, 'FLAT_DIAGONAL')", "deltas = diagonal.atoms(sp.DiracDelta)", "z0out = diagonal.xreplace({d: sp.S.One for d in deltas})", "DIMENSION_SCHEMA[z.name] = dimension(z0out)"]
    statements=[U.signature(n) for n in native_bridge.body]
    require(all(U.signature(ast.parse(text).body[0]) in statements for text in required_statements),'actual dynamic dimension override routing')
    J.emit('native-flat-override-source',{'dtnLiteralSha256':hashlib.sha256(dtn['literal'].encode()).hexdigest(),'kernelBridge':contract,'requiredStatements':required_statements})
    face_results={};controls=[]
    for face,sign in [('plus',1),('minus',-1)]:
        response=next(r for r in native['faceResponseSources']['cases'] if r['case']==['LAB_HELD',sign,'RHO4_CONSTANT'])
        velocity=next(r for r in native['geometry']['face_velocity']['cases'] if r['case']==['LAB_HELD',sign,'DELTA_W'])
        for name,case in [('response',response),('velocity',velocity)]:
            J.emit(face+'-'+name+'-native-case',case)
            join_native_case(J,U,face+'-'+name+'-original-case',sources[case['provenance']['source']+':'+str(case['provenance']['valueLine'])],case,('CASES',) if name=='response' else ())
        selected_dtn=[]
        for case in dtn_cases:
            labels,payload=U.tuple_args(case);labels=U.tuple_args(labels)
            if U.kind(labels[0])=='Str' and labels[0].args[0].value=='LAB_HELD' and U.rational(labels[1])==sign:selected_dtn.append(U.tagged(payload,'VALUE'))
        require(len(selected_dtn)==1,'unique native dtn face')
        diagonal=U.tagged(selected_dtn[0],'FLAT_DIAGONAL');require(U.kind(diagonal)=='Mul','native flat diagonal product')
        delta_nodes=[v for v in diagonal.args if U.kind(v)=='DiracDelta'];remaining=[v for v in diagonal.args if U.kind(v)!='DiracDelta']
        delta_removed=ast.Call(func=ast.Name(id='Mul',ctx=ast.Load()),args=remaining,keywords=[])
        J.emit(face+'-actual-flat-delta-removal',{'diagonal':ast.unparse(diagonal),'deltaFactors':[ast.unparse(v) for v in delta_nodes],'deltaRemoved':ast.unparse(delta_removed),'savedFlatJoin':flat})
        require(len(delta_nodes)==3 and U.same(delta_removed,U.parse(flat['left']['srepr'])),'actual three-dimensional flat coefficient join')
        flat_unit=certify(J,U,face+'-delta-removed-native-flat-units',ast.unparse(delta_removed),registry,[-3,-1,1])
        saved=raw['consumer/'+face+'-source-input.json'];require(U.same(U.parse(velocity['valueConstructorText']),U.parse(saved['nativeVelocity']['srepr'])),'native velocity argument')
        certify(J,U,face+'-native-velocity-units',velocity['valueConstructorText'],registry,[1,-1,0])
        resp=U.parse(response['valueConstructorText']);dp=U.tagged(resp,'DELTA_P');require(U.kind(dp)=='Mul','native deltaP product')
        source_add=[v for v in dp.args if U.kind(v)=='Add'];require(len(source_add)==1,'one complete source addend')
        rawsource=U.parse(saved['raw']['srepr']);require(U.kind(rawsource)=='Mul','normalized raw source product')
        eps_inv=[v for v in rawsource.args if U.kind(v)=='Pow' and U.kind(v.args[0])=='Symbol' and U.symbol_name(v.args[0])=='epsilon_shape' and U.rational(v.args[1])==-1]
        remain=[v for v in rawsource.args if v not in eps_inv]
        J.emit(face+'-native-source-selector',{'deltaP':ast.unparse(dp),'selectedSource':ast.unparse(source_add[0]),'savedRaw':saved['raw'],'epsilonDivisors':[ast.unparse(v) for v in eps_inv],'remaining':[ast.unparse(v) for v in remain]})
        require(len(eps_inv)==1 and len(remain)==1 and U.same(remain[0],source_add[0]),'actual normalized native source identity')
        source_unit=certify(J,U,face+'-unbound-source-units',saved['raw']['srepr'],registry,[1,-1,0])
        # Join dynamic native Z unit override to actual flat operand. This is the
        # original kernel_bridge rule, not a correction to its static placeholder.
        zname='s11cc1_dtn_operator_lab_held_'+face;face_registry=dict(registry);face_registry[zname]=flat_unit
        definition=U.tuple_args(U.tagged(resp,'RESOLVENT_DEFINITION'))
        require(len(definition)==3 and U.same(definition[0],U.tagged(resp,'RESOLVENT')),'native resolvent definition')
        require(U.kind(definition[1])=='Add','native inverse operand sum')
        feedback=[v for v in definition[1].args if U.kind(v)=='Mul' and any(U.kind(t)=='Symbol' and U.symbol_name(t)==zname for t in v.args)]
        require(len(feedback)==1,'one original feedback coefficient')
        feedback_factors=[v for v in feedback[0].args if not (U.kind(v)=='Symbol' and U.symbol_name(v)==zname)]
        coefficient=ast.Call(func=ast.Name(id='Mul',ctx=ast.Load()),args=feedback_factors,keywords=[])
        certify(J,U,face+'-unbound-feedback-coefficient-unit',ast.unparse(coefficient),registry,[3,1,-1])
        J.emit(face+'-dynamic-Z-contract',{'flat':flat,'actualFlatConstructor':ast.unparse(delta_removed),'flatUnitReturnOperation':face+'-delta-removed-native-flat-units','staticZUnit':registry[zname],'actualAppliedSymbolUnit':face_registry[zname],'kernelBridge':raw['source-contracts.json']['kernel_bridge'],'definition':ast.unparse(definition[1])})
        certify(J,U,face+'-native-resolvent-homogeneity',ast.unparse(definition[1]),face_registry,[0,0,0])
        certify(J,U,face+'-native-pressure-unit',ast.unparse(dp),face_registry,[-2,-2,1])
        wrong=dict(face_registry);wrong[zname]=[0,0,0]
        for label,expr,reg in [('opaque-Z-placeholder',definition[1],wrong),('brane-density-used-as-bulk',rawsource,dict(registry,rho_m=[-3,0,1]))]:
            J.start(face+'-control-'+label,{'actualConstructor':ast.unparse(expr),'mutatedRegistry':reg,'baselineRegistry':face_registry if label=='opaque-Z-placeholder' else registry})
            checker=U.Units(reg);message=None
            try:checker.dimension(expr)
            except ValueError as exc:message=str(exc)
            finally:J.emit(face+'-control-'+label+'-walk',checker.events)
            result={'name':label,'face':face,'refusal':message,'unitSensitivityOnly':True};J.emit(face+'-control-'+label+'-decision',result)
            require(message is not None and message.startswith('inhomogeneous native Add'),'actual unit mutant refusal');J.finish(result);controls.append(result)
        face_results[face]={'sourceUnit':source_unit,'nativePressureUnit':[-2,-2,1],'actualNativeSourceJoined':True,'dynamicZFromFlat':True}
    # Consumer dimensions were already checked. Restore actual complete coefficient
    # operands/results; do not rerun those dimension calculations.
    consumers=raw['consumer/consumer-unit-joins.json'];selected=[v for v in consumers if v['row']=='THETA_BALANCE']
    J.emit('inherited-consumer-unit-returns',selected)
    require(len(selected)==4 and {v['slot']['text'] for v in selected}=={'delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus'},'all four inherited THETA slots')
    require(all(v['total']==v['expected'] for v in selected),'inherited consumer unit results')
    J.emit('pending-pressure-unit-transport',{'sourceGradeFieldTransport':'REQUIRED; not supplied by this native bridge','all544AddressSummandProof':'REQUIRED','completeHJDProfileAndMeasureUnits':'REQUIRED','nativeConsumers':'inherited complete operands/results, not recomputed','noNumericalActionOrCacheClearance':True})
    return {'status':'NATIVE_PRESSURE_UNIT_BRIDGE_COMPLETE_PENDING_INSPECTION','faceResults':face_results,'controls':controls,'registryUnavailable':unavailable,'sourceGradeFieldUnitTransportComplete':False,'pressureSummandUnitsComplete':False,'numericalEvaluatorReady':False,'sourceConstructorsCalled':False,'integralsEvaluated':0,'scope':m['scope'],'savedCopies':len(copies)}


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
