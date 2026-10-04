#!/usr/bin/env python3
"""Guarded exact outer geometry and static dependency planning; no integration.

This is the geometry subset of pressure-readiness. Missing pressure summand
units remain a separate prerequisite. Saved symbolic JSON is joined as data,
never decoded into a CAS and never used to replay a completed function.
"""
import argparse
import ast
import copy
from fractions import Fraction as F
import hashlib
import importlib.util
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
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO={'text':'0','srepr':'Integer(0)'}
LIVE='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED'


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


def verify_gate(path,manifest_path,m):
    g=read(path)
    require(g['status']=='READY_FOR_ONE_PACKET_GEOMETRY','gate status')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker and manifest pins')
    require(g['sourcePins']==m['sourcePins'],'complete source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority','library'):
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and
        g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual containment path')
    require(g['launcher']==m['launcher'] and g['library']==m['librarySource'] and
        g['buildReviewRecord']==m['reviewRecordWillBe'] and g['authority']==m['executionAuthority'],'actual document routes')
    r=read(g['buildReviewRecord'])
    require(g['independentBuildClearance'] is True and r['independentBuildClearance'] is True and r['allChecksPassed'] is True,'actual paired build assessment')
    for key in ('workerSha256','manifestSha256','librarySha256','launcherSha256','sharedGuardSha256','supervisorSha256'):
        require(r[key]==g[key],'review/gate '+key)
    require(all(r['reports'][reviewer]['literalVerdict']=='CLEAR FOR THIS BOUNDED PACKET-GEOMETRY BUILD' for reviewer in ('claude','grok')),'literal build reports')
    method=read(m['methodRecord'])
    require(g['methodRecordSha256']==sha(m['methodRecord']) and method['jointIndependentMethodClearance'] is True and
        method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'actual assessed method')
    authority=read(g['authority'])
    require(authority['boundedInstrumentAuthorized'] is True and authority['scope']==m['scope']==g['scope'] and
        authority['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and
        authority['automaticScientificRetry'] is False and authority['noDeadline'] is True and g['durationLimits'] is None,'standing bounded authority')
    return g


def load_geometry(m):
    spec=importlib.util.spec_from_file_location('s11c_packet_exact_geometry',m['librarySource'])
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module


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


def dependency_inventory(raw,J,m):
    full=raw['inventory/THETA_BALANCE-ordered-addresses.json'];selection=raw['selected/pressure-addresses.json']
    addresses=[a for a in full if a['row']=='THETA_BALANCE' and a['jet']['channel']=='e_W']
    J.emit('selected-arguments',{'selection':selection,'originalReceipt':m['savedInputs']['inventory/THETA_BALANCE-ordered-addresses.json']})
    require(addresses==selection['selected'] and len(addresses)==544 and len({a['addressId'] for a in addresses})==544,'all actual selected addresses')
    original=m['savedInputs']['inventory/THETA_BALANCE-ordered-addresses.json']
    require(selection['sourceSha256']==original['sha256'] and selection['sourceBytes']==original['bytes'],'selection original byte join')
    old=raw['preflight/address-adapter-result.json'];byid={a['addressId']:a for a in old}
    require(len(old)==len(byid)==544,'complete prior adapter census')
    fields=raw['pressure/fields.json'];require(len(fields)==34,'complete original fields')
    for fid,field in fields.items():
        inp=raw['field/'+fid+'-reconstruction-input.json'];ret=raw['field/'+fid+'-reconstruction-return.json']
        J.emit('field-'+fid+'-inherited',{'field':field,'polynomial':raw['field/'+fid+'-polynomial.json'],'input':inp,'return':ret})
        require(inp['left']==field['field'] and ret['cancelled']==ZERO,'inherited field proof with original operand')
    adapters=raw['preflight/numeric-factor-adapters.json'];require(len(adapters['definitions'])==20,'all complete numeric templates')
    joins={a['addressId']:a for a in adapters['addressJoins']};require(len(joins)==544,'complete original template joins')
    selected_by_id={a['addressId']:a for a in addresses}
    for label,adapter in adapters['definitions'].items():
        arguments=raw['preflight/numeric-factor-'+label+'-arguments.json']
        inp=raw['preflight/numeric-factor-'+label+'-input.json'];ret=raw['preflight/numeric-factor-'+label+'-return.json']
        J.emit('numeric-'+label+'-inherited',{'definition':adapter,'arguments':arguments,'input':inp,'return':ret})
        require(arguments['address']==selected_by_id[adapter['firstAddressId']] and
            arguments['actualCompleteFactor']==adapter['original'],'original complete template arguments')
        require(inp['left']==adapter['mapped'] and inp['right']==adapter['template'] and ret['cancelled']==ZERO,'complete template zero join')
    family_maps={'X':{},'Y':{},'response':{}};entries=[];inherited=set()
    def intern(kind,key,address):
        encoded=canonical(key)
        if encoded not in family_maps[kind]:family_maps[kind][encoded]={'signature':key,'addressIds':[]}
        family_maps[kind][encoded]['addressIds'].append(address)
        return hashlib.sha256(encoded.encode()).hexdigest()
    for a in addresses:
        aid=a['addressId'];o=byid[aid];name=a['fullFactorProof']['proof'];op=raw['factors/'+name+'-operands.json']
        label='-'.join(a[k] for k in ('face','slot','component'));adapter=adapters['definitions'][label]
        J.emit('address-'+str(aid)+'-arguments',{'address':a,'priorAdapter':o,'factorOperands':op,'numericJoin':joins[aid],'numericDefinition':adapter})
        require(o['status']==a['status'] and o['proof']==name and o['jet']==a['jet'],'prior address argument identity')
        for role in ('source','consumer'):
            fid=a[role+'Transform']['coefficientId']
            require(fields[fid]['field']==a[role+'Field'] and o[role+'FieldId']==fid,'actual '+role+' field identity')
        for key in ('sha256','original','mapped','map','symbolAssumptions','flatSupport','frequency','positiveRegulatorContinuation'):
            require(op['actualResponseMap'][key]==a['responseMap'][key],'actual factor argument '+key)
        require(op['addressNormalOriginal']==a['normalOriginal'] and op['requiredMap']==a['fullFactorProof']['completeNormalMap'],'complete native normal map')
        if name not in inherited:
            inp=raw['factors/'+name+'-full-mapped-residual-input.json'];normal=raw['factors/'+name+'-normal-source-join-input.json']
            ret=raw['factors/'+name+'-full-mapped-residual-return.json'];nret=raw['factors/'+name+'-normal-source-join-return.json']
            J.emit(name+'-inherited',{'operands':op,'input':inp,'return':ret,'normalInput':normal,'normalReturn':nret})
            require(inp['right']==op['mappedAddressFactor'] and normal['left']==op['addressNormalOriginal'] and normal['right']==op['savedNormal'],'actual inherited factor operands')
            require(ret['cancelled']==nret['cancelled']==ZERO,'inherited literal factor zeros');inherited.add(name)
        # Numeric definitions are joined through their original, already computed
        # input/return records; do not recalculate the saved expression product.
        numeric=raw['preflight/numeric-factor-'+label+'-input.json'];nret=raw['preflight/numeric-factor-'+label+'-return.json']
        require(numeric['left']==adapter['mapped'] and numeric['right']==adapter['template'] and nret['cancelled']==ZERO,'actual numeric factor proof')
        require(adapter['original']==op['mappedAddressFactor'],'native full factor to original numeric adapter')
        require(joins[aid]['adapter']==label,'numeric adapter address key')
        flat=a['component']=='NATIVE_FLAT'
        require(a['sourceTransform']['transfer']==('l-p' if flat else 'k-p') and a['consumerTransform']['transfer']=='r-l' and
            a['responseInputDepth']==('q(l)' if flat else 'q(k)') and a['responseOutputDepth']=='q(l)','native source/output/depth ordering')
        live=a['status']==LIVE
        if live:require(a['epsilonCount']==1,'one live epsilon')
        else:
            require(a['status'] in ('EXACT_ZERO_SOURCE_JET','EXACT_ZERO_CONSUMER') and a['epsilonCount']==0,'explicit zero status')
            require(a['consumerField' if a['status']=='EXACT_ZERO_CONSUMER' else 'sourceField']==ZERO,'actual exact zero operand')
        source={'field':a['sourceField'],'jet':a['jet'],'waveMultiplier':a['waveMultiplier'],
            'order':'derivatives of Gaussian before source coefficient multiplication','normalization':'1/(2*pi)',
            'units':'PENDING_UNBOUND_PRESSURE_SUMMAND_CERTIFICATE'}
        consumer={'field':a['consumerField'],'order':'coefficient times dual Gaussian; no complex conjugation',
            'normalization':'2*pi*hat[c*v](-l)','units':'PENDING_UNBOUND_PRESSURE_SUMMAND_CERTIFICATE'}
        response={k:a[k] for k in ('face','slot','component','responseGrade','responseOriginal','responseMap','normalOriginal','normalMultiplier','responseInputDepth','responseOutputDepth','measure','freeMomenta')}
        response['completeNumericAdapter']=adapter;response['completeNativeFactor']=op['mappedAddressFactor']
        ids={kind:intern(kind,key,aid) for kind,key in [('X',source),('Y',consumer),('response',response)]} if live else {}
        entries.append({'addressId':aid,'status':a['status'],'face':a['face'],'slot':a['slot'],'component':a['component'],
            'targetGrade':a['targetGrade'],'gradeTriple':[a[k] for k in ('consumerGrade','responseGrade','sourceGrade')],
            'families':ids,'zeroAddressesPreserved':not live,'unitCertificatePending':True})
    require(sum(e['status']==LIVE for e in entries)==102,'all 102 formal products, not asserted nonzero')
    families={kind:[{'id':hashlib.sha256(k.encode()).hexdigest(),**v} for k,v in group.items()] for kind,group in family_maps.items()}
    J.emit('dependency-families',{'families':families,'entries':entries,'scope':'candidate mathematical families only; no numerical cache authorized',
        'rule':'actual carrier/center/width/momentum/rule/route/precision/tail/unit certificate still required',
        'sharedEvaluatedValuesAcrossRoutes':False,'oldBankRequestMatches':0})
    return entries,families


def count_work(entries,families,summary,domain,G):
    counts=G.static_counts(summary);slabs=summary['slabs']
    # Keep both ungrouped address occurrences and exactly identical formal
    # family occurrences. Distinct numerical coordinates are not inferred.
    primitive={'NATIVE_FLAT':['flat'],'NATIVE_HEIGHT':['height-contact','height-paired-PV'],
        'NATIVE_SLOPE':['slope'],'NATIVE_MIXED_ITERATION':['H','J'],
        'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':['D']}
    records=[]
    for comp,names in primitive.items():
        a=[e for e in entries if e['component']==comp and e['status']==LIVE]
        for name in names:
            applicable=(domain=='height') if name=='height-paired-PV' else (domain=='square')
            if not applicable:continue
            per={}
            for route,n in [('A24',24),('A48',48),('B_initial',None)]:
                one=slabs*(2*n if n else 15)
                nodes=one if name in ('flat','height-contact') else counts[route].get('outerNodeOccurrences',counts[route].get('outerNodeOccurrencesLowerBound'))
                per[route]={'nodeOccurrencesPerFamily':nodes,'addressProductOccurrences':nodes*len(a),
                    'candidateResponseFamilyOccurrences':nodes*len({e['families']['response'] for e in a}),
                    'candidateXFamilyOccurrencesBeforeCoordinateReuse':nodes*len({e['families']['X'] for e in a}),
                    'candidateYFamilyOccurrencesBeforeCoordinateReuse':nodes*len({e['families']['Y'] for e in a}),
                    'heightYBranchesPerNode':2 if name=='height-paired-PV' else 1,
                    'lowerBoundOnAdaptiveWork':route=='B_initial'}
            records.append({'primitive':name,'liveAddressIds':[e['addressId'] for e in a],
                'perGrade':[{'grade':list(g),'addressIds':[e['addressId'] for e in a if e['targetGrade']==list(g)]} for g in ((0,0),(1,0),(0,1),(1,1))],
                'routes':per,'addressedScalarReturnOccurrencesA24A48WithoutSharing':sum(v['addressProductOccurrences'] for k,v in per.items() if k!='B_initial'),
                'recordCountConvention':'Conditional unshared format: one scalar address/primitive return record per occurrence. A grouped or vector record format can have fewer records; no serialized-size or required-RSS inference follows.'})
    return {'counts':counts,'primitives':records,'flatContactRule':'Conservative reuse of all square k slabs; flat/contact do not acquire a second outer integral',
        'innerFourierAdaptiveNodeCounts':'UNKNOWN','serializedBytesAndPeakRSS':'UNKNOWN; prior database sizes do not prove native memory failure',
        'wallTimeEstimate':None,'uniformQuadratureErrorClaim':False,'candidateSharingNotEvaluatorReadiness':True}


def run(m,J,G):
    raw,copies=copy_inputs(m,J)
    J.emit('prior-bank-context',{'records':{alias:value for alias,value in raw.items() if alias.startswith('accepted/')},
        'scope':'Historical counters/opaque database receipts only. No numerical bank reopened or matched, no throughput or process/cache inference.'})
    rules={n:raw['rules/'+n+'.json'] for n in ('A-GL24','A-GL48','B-G7-K15')}
    J.emit('immutable-rule-count-operands',rules)
    for n,order in [('A-GL24',24),('A-GL48',48)]:
        require(len(rules[n]['nodes'])==len(rules[n]['weights'])==order and rules[n]['precision']==30,'saved open A rule count')
    require(len(rules['B-G7-K15']['kronrodNodes'])==len(rules['B-G7-K15']['kronrodWeights'])==15 and
        len(rules['B-G7-K15']['gaussNodes'])==7 and rules['B-G7-K15']['precision']==50,'saved physical B rule count')
    physical=raw['preflight/physical-plan.json'];params=raw['physical-input.json']['parameters'];context=raw['local/context.json'];pre=raw['preflight/checks.json']
    J.emit('physical-plan-arguments',{'physical':physical,'original':raw['physical-input.json'],'context':context,'preflightChecks':pre})
    require(context['physical']==context['saved']['physicalInput']==raw['physical-input.json'] and context['frequencyOverride']=={'old':'1','actual':3},'same real-frequency sources')
    require(physical['frequency']==3 and physical['cs']['srepr']=='Mul(Rational(1, 2), Pow(Integer(6), Rational(1, 2)))' and
        physical['kappa']['srepr']=='Mul(Rational(1, 10), Pow(Integer(595), Rational(1, 2)))','saved effective speed/matching operands')
    require(physical['carriers']==[physical['kappa'],ZERO] and physical['s']==8 and
        physical['centers']==[{'text':'-5/2','srepr':'Rational(-5, 2)'},{'text':'5/2','srepr':'Rational(5, 2)'}],'same two packets')
    require(params['L_W']=='10' and params['W_0']=='1' and params['s11cdTangentialMomentum1']=='1/5' and params['s11cdTangentialMomentum2']=='1/10','native scale and tangents')
    require((pre['K'],pre['U'],pre['T'])==(27,75,122),'actual accepted radii, no retuning')
    require(pre['status']=='PACKET_PREFLIGHT_CERTIFICATES_COMPLETE_NO_ACTION' and pre['addresses']==544 and pre['liveAddresses']==102,'accepted preflight scope')
    scale=raw['inventory/native-profile-scale-join.json']
    require(scale['savedLength']['srepr']=='Integer(10)' and scale['physicalLength']==params['L_W'] and scale['declaredLength']==10,'native profile length origin')
    entries,families=dependency_inventory(raw,J,m)
    kappa=G.Quad(0,1,F(119,20));width=F(physical['s']);length=F(params['L_W']);summaries=[];controls=[]
    # This exact constant bridge is new geometry arithmetic, not a replay of
    # matching dispersion, a saved root construction, or a scientific producer.
    require(kappa*kappa==F(595,100),'exact positive geometry basis bridge')
    J.emit('geometry-basis',{'savedKappa':physical['kappa'],'positiveBasis':kappa,'square':kappa*kappa,'width':width,'length':length})
    for carrier_id,carrier in [('matching',kappa),('zero',kappa*0)]:
        for enlargement in (0,2):
            K=pre['K']+enlargement;U=pre['U'];T=pre['T']+enlargement
            specs=G.specifications(kappa,K,U,carrier,width,length)
            for domain,(lines,cuts,box) in specs.items():
                name=carrier_id+'-K'+str(K)+'-'+domain
                args={'carrier':carrier,'K':K,'U':U,'T':T,'lines':lines,'mandatoryCuts':cuts,'box':box}
                J.start(name,args)
                plan=G.arrangement(lines,*box,cuts);J.emit(name+'-full-arrangement',plan)
                result=G.audit(plan,lines,cuts);counts=count_work(entries,families,result,domain,G)
                J.emit(name+'-work-counts',counts);J.finish({'audit':result,'counts':counts,'noNodesOrIntegralsEvaluated':True})
                summaries.append({'name':name,'domain':domain,'carrier':carrier_id,'K':K,'U':U,'T':T,'audit':result,'counts':counts})
                if carrier_id=='matching' and enlargement==0 and domain=='square':
                    mutated=[l for l in lines if 'collision:sum:2' not in l.labels]
                    require(len(mutated)==len(lines)-1,'actual required collision-line removal')
                    J.start('control-missing-collision',{'baseline':name,'removedLine':[l for l in lines if l not in mutated],'mutantLines':mutated,'requiredLines':lines,'requiredCuts':cuts})
                    bad=G.arrangement(mutated,*box,cuts);J.emit('control-missing-collision-arrangement',bad)
                    mutant_coverage=G.audit(bad,mutated,cuts)
                    try:G.audit(bad,lines,cuts)
                    except ValueError as exc:message=str(exc)
                    else:raise ValueError('missing-collision control failed to refuse')
                    require(message=='required collision/resolution incidence','addressed missing-line refusal')
                    controls.append({'name':'missing-collision','refusal':message,'mutantStillCoversBox':mutant_coverage['coverage']});J.finish(controls[-1])
                    bad=copy.deepcopy(plan);cell=bad['slabs'][0]['cells'][0];cell['lowerLine'],cell['upperLine']=cell['upperLine'],cell['lowerLine']
                    J.start('control-reversed-cell',{'baseline':name,'baselineCell':plan['slabs'][0]['cells'][0],'mutantCell':cell,'slabIndex':0,'cellIndex':0})
                    try:G.audit(bad,lines,cuts)
                    except ValueError as exc:message=str(exc)
                    else:raise ValueError('reversed-cell control failed to refuse')
                    require(message=='oriented adjacent cell','addressed orientation refusal')
                    controls.append({'name':'reversed-cell','refusal':message});J.finish(controls[-1])
    J.emit('future-request-identity-contract',{'actualCoordinates':'Exact mathematical descriptors: plan/slab/cell, exact affine endpoints, side, immutable rule and node index. Independently rounded by each route; printed values never keys.',
        'X':['complete field operand and native jet/wave order','carrier/center/width','physical transform argument','normalization','unit certificate','precision/route/contour/radius/tail settings'],
        'response':['face/slot/component','complete native factor and map','physical k/l/depth bindings','omega/cs/edges/materials/W/L','route/order/precision/tail/inner panel settings'],
        'independentRoutes':'A24/A48 and physical B have independent evaluated transforms, inner arrays and outer samples. Immutable source/rule provenance may be shared. No fixed-bank convenience complete assembler called.',
        'oldBankRequestsReused':0,'units':'PENDING_SEPARATE_PRESSURE_CERTIFICATES','numericalEvaluatorReady':False})
    return {'status':'BOUNDED_GEOMETRY_DEPENDENCY_PLAN_COMPLETE_PENDING_INSPECTION','plans':summaries,'addresses':544,'formalAddresses':102,'explicitZeroAddresses':442,
        'candidateFamilyCounts':{k:len(v) for k,v in families.items()},'controls':controls,'pressureSummandUnits':'PENDING_SEPARATE_REQUIRED_CERTIFICATE',
        'numericalEvaluatorReady':False,'integralsEvaluated':0,'completedFunctionsReplayed':False,'numericalCacheValuesReused':0,'fullAdaptiveCost':'UNKNOWN',
        'scope':m['scope'],'savedCopies':len(copies)}


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);args=p.parse_args()
    m=read(args.inputs);gate=verify_gate(args.gate,args.inputs,m)
    document_pins={str(args.inputs):gate['manifestSha256'],str(args.gate):sha(args.gate)}
    for key in ('buildReviewRecord','authority'):document_pins[gate[key]]=gate[key+'Sha256']
    document_pins[m['methodRecord']]=gate['methodRecordSha256']
    tail=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(sys.argv==tail and gate['command'][-len(tail):]==tail and str(args.out)==gate['outputDirectory'],'actual invocation and output')
    # Import only the original inert containment function, never an old main or science function.
    tree=ast.parse(Path(m['helperSource']).read_text());functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='containment']
    require(len(functions)==1,'one inert containment helper')
    ns={'require':require,'Path':Path,'os':os,'resource':resource,'THREADS':THREADS}
    exec(compile(ast.Module(body=functions,type_ignores=[]),m['helperSource'],'exec'),ns)
    enforced=ns['containment']();args.out.mkdir(exist_ok=False);J=Journal(args.out);start=time.monotonic();error=None
    try:
        identities={}
        for index,(path,digest) in enumerate(document_pins.items()):
            target=args.out/('identity-'+str(index)+'-'+Path(path).name)
            require(sha(path)==digest,'identity source');shutil.copyfile(path,target)
            require(sha(target)==digest,'identity copy')
            identities[path]={'path':target.name,'sha256':digest,'bytes':target.stat().st_size}
        J.emit('additional-identity-copies',identities)
        J.emit('actual-containment',enforced);G=load_geometry(m);result=run(m,J,G)
    except BaseException:
        error=traceback.format_exc();result={'status':'FAILED_PRESERVED','failure':error,'activeOperation':J.active,'completeOperations':J.completed}
        J.emit('failure',result)
    finally:
        post={p:posthash(p,h) for p,h in {**m['sourcePins'],**document_pins}.items()}
        identity_post={v['path']:posthash(args.out/v['path'],v['sha256']) for v in identities.values()}
        copied={str(p.relative_to(args.out)):sha(p) for p in (args.out/'saved').rglob('*') if p.is_file()}
        expected={'saved/'+a:r['sha256'] for a,r in m['savedInputs'].items()}
        copy_intact=copied==expected
        J.emit('posthashes',{'sources':post,'identities':identity_post,'copied':copied,'expectedCopied':expected,'copiesIntact':copy_intact})
        result.update(wallMilliseconds=round((time.monotonic()-start)*1000),sourcePosthashesIntact=all(v['intact'] for v in post.values()),
            identityCopiesIntact=all(v['intact'] for v in identity_post.values()) and len(identities)==len(document_pins),copiesIntact=copy_intact)
        result['completedOperations']=J.completed;result['activeOperation']=J.active
        J.emit('journal-result',result)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text());sys.stdout.flush()
    if error:sys.stderr.write(error);return 1
    require(result['sourcePosthashesIntact'] and result['copiesIntact'] and result['identityCopiesIntact'],'posthash integrity')
    return 0


if __name__=='__main__':sys.exit(main())
