#!/usr/bin/env python3
"""Native local coefficients plus inherited weak pressure: no integral or inverse."""
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

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS = ('require','sha','save','replace_json','posthash_records','containment',
           'Journal','decode','one_symbol','expanded_sinh_arguments')
G = ((0,0),(1,0),(0,1),(1,1))
ROWS = ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS = ('u_1','u_2','u_3','theta','e_W')
SLOTS = ('delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus')


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

def memory_domain_certificate(den, bind, constant, terms, zero_record, emit, sp):
    """Exact native time-sign recognition including bound prefactors and powers."""
    rec={'raw':den,'bound':bind(den)}
    emit('memory-domain-input',rec)
    rec['certificate']=constant(rec['bound'],False)
    symbols=list(den.free_symbols);names=[s.name for s in symbols]
    require(len(names)==len(set(names)),'unique native denominator symbols')
    taus=[s for s in symbols if s.name in ('tau_A','tau_V','tau_X')]
    if taus:
        require(len(taus)==1 and set(names)<={'omega','rho_br','W_0','L_W',taus[0].name}
                and 'omega' in names,'native local memory denominator form')
        om=next(s for s in symbols if s.name=='omega');tau=taus[0]
        powers=terms(den,om,sp.Symbol('full_weak_memory_unused'))
        require(powers and all(j==0 and type(i) is int and i>=0 for i,j in powers),
                'native memory polynomial powers')
        degree=max(i for (i,j),v in powers.items() if v!=0)
        require(degree in (1,2),'source-censused native memory power')
        scale=sp.cancel(den.subs(om,0))
        rec.update(omega=om,tau=tau,power=degree,prefactor=scale,
                   coefficients=[{'power':i,'value':v} for (i,j),v in powers.items()],
                   expected=scale*(1-sp.I*om*tau)**degree)
        emit('memory-sign-operands',rec)
        require({s.name for s in scale.free_symbols}<={'rho_br','W_0','L_W'},
                'frequency-independent native memory prefactor')
        rec['prefactorCertificate']=constant(bind(scale),False)
        zero_record(rec,'time-memory-sign',den,rec['expected'])
    emit('memory-domain-return',rec)
    return rec


def verify_build_assessment(gate, review):
    require(review['methodAssessed'] is True and review['allChecksPassed'] is True,
            'actual completed build assessment')
    if review['independentBuildClearance'] is True:
        require(gate['independentBuildClearance'] is True,'honest paired build clearance')
        for engine in ('claude','grok'):
            require(review['reports'][engine]['literalVerdict']==
                    'CLEAR FOR THIS BOUNDED FULL-WEAK BUILD','literal build verdict')
        for key in ('workerSha256','manifestSha256','launcherSha256'):
            require(review[key]==gate[key],'review/gate '+key)
    else:
        require(gate['independentBuildClearance'] is False and
                gate['localToolingRepairAccepted'] is True,'honest tested local repair authority')
        require(review['reports']['claude']['literalVerdict']=='CLEAR FOR THIS BOUNDED FULL-WEAK BUILD'
                and review['reports']['grok']['literalVerdict']=='NEEDS REVISION','preserved literal build pair')
        require(sha(gate['repairRecord'])==gate['repairRecordSha256'],'local repair record pin')
        repair=json.loads(Path(gate['repairRecord']).read_text())
        require(repair['toolingOnly'] is True and repair['testsPassed'] is True
                and repair['noScientificPayloadRestored'] is True,'tested exact-predicate repair')
        require(repair['reviewRecordSha256']==gate['buildReviewRecordSha256'], 'actual assessed baseline')
        for key in ('workerSha256','manifestSha256','launcherSha256'):
            require(repair[key]==gate[key] and repair['reviewed'][key]==review[key],
                    'reviewed/repaired identity '+key)
        require(repair['methodSha256']==review['methodSha256'],'unchanged assessed method')
        for p,h in repair['evidencePins'].items():require(sha(p)==h,'repair evidence '+p)
    for key in ('sharedGuardSha256','supervisorSha256'):
        require(review[key]==gate[key],'review/gate '+key)


def verify_gate(path, manifest_path, manifest):
    gate = json.loads(Path(path).read_text())
    verify_helper_paths(gate)
    require(gate['status']=='READY_FOR_ONE_FULL_WEAK_INSTRUMENT','gate status')
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
    verify_build_assessment(gate,review)
    require(review['methodSha256']==sha(manifest['methodPath']), 'exact corrected method')
    authority=json.loads(Path(gate['authority']).read_text())
    require(gate['authority']==manifest['executionAuthority'] and
            manifest['sourcePins'][gate['authority']]==gate['authoritySha256'], 'actual authority route and pin')
    require(authority['boundedInstrumentAuthorized'] is True and
            authority['automaticScientificRetry'] is False and
            authority['scienceExecutionsAuthorized']==1 and authority['noDeadline'] is True and
            authority['scope']==manifest['scope'], 'standing bounded authority')
    require(gate['durationLimits'] is None and gate['scientificRunsAuthorized']==1 and
            gate['scope']==manifest['scope'], 'one bounded no-deadline job')
    method_record=json.loads(Path(manifest['methodRecord']).read_text())
    require(method_record['jointIndependentMethodClearance'] is True and
            method_record['methodSha256']==sha(manifest['methodPath']),'method assessment and exact source')
    return gate

def validate_partition(record):
    tree=ast.parse(record['fullConstructor'],mode='eval').body
    constructor_bytes=record['fullConstructor'].encode()
    require(b'\n' not in constructor_bytes,'native constructor single-line encoding')
    require(isinstance(tree,ast.Call) and isinstance(tree.func,ast.Name) and tree.func.id=='Add','native row Add')
    children=record['children'];require(len(children)==len(tree.args),'complete child count')
    local=set(record['localChildIndices']);pressure=set(record['pressureChildIndices'])
    require(not local&pressure and local|pressure==set(range(len(children))),'disjoint complete partition')
    for i,(entry,node) in enumerate(zip(children,tree.args)):
        require(node.lineno==node.end_lineno==1,'native scalar source location')
        text=constructor_bytes[node.col_offset:node.end_col_offset].decode()
        require(entry['childIndex']==i and entry['constructorText']==text and
                hashlib.sha256(text.encode()).hexdigest()==entry['sha256'],'actual native child bytes')
        hits=[]
        for call in ast.walk(node):
            if isinstance(call,ast.Call):
                require(isinstance(call.func,ast.Name) and call.func.id in ('Add','Mul','Pow','Symbol','Integer','Rational'),'native scalar constructor')
                if call.func.id=='Symbol':
                    name=ast.literal_eval(call.args[0]);require(isinstance(name,str),'literal native symbol')
                    if 'delta_p' in name or 'd_w_' in name:
                        require(name in SLOTS,'unknown native pressure');hits.append({'constructor':'Symbol','name':name})
        require(hits==entry['hits'] and bool(hits)==(i in pressure) and entry['completeNameCoverage'] is True,'actual pressure complement')
    require(not record['unknownPressureAtoms'],'no unknown pressure names')
    return {'children':len(children),'local':len(local),'pressure':len(pressure),'exactSourcePartition':True}


def numeric_extension(physical, shared, symbolic_names):
    # Values still serialized here: no silently bound eta or speed.
    allowed={k:v for k,v in physical.items() if k not in set(symbolic_names)|{'c_s0'}}
    allowed['omega']='3'
    return allowed


def pressure_child_coverage(partition, addresses, row):
    """New assembly join: actual child hashes per slot, including zero U rows."""
    by_slot={s:[] for s in SLOTS};expected=set()
    for child in partition['children']:
        if child['childIndex'] not in partition['pressureChildIndices']:continue
        require(child['hits'] and all(h['name'] in by_slot for h in child['hits']),'known addressed pressure child')
        expected.add(child['sha256'])
        for slot in {h['name'] for h in child['hits']}:by_slot[slot].append(child['sha256'])
    covered=set()
    for a in addresses:
        require(a['row']==row and a['face'] in ('plus','minus') and a['slot'] in ('pressure','normal'),'native address selection')
        slot=('delta_p_' if a['slot']=='pressure' else 'd_w_delta_p_')+a['face']
        require(a['nativeRowChildHashes']==by_slot[slot],'address to actual native slot children')
        covered.update(a['nativeRowChildHashes'])
    require(covered==expected,'complete native pressure-child union')
    return {'row':row,'bySlot':by_slot,'coveredHashes':sorted(covered),'expectedHashes':sorted(expected),'addresses':len(addresses)}


def inherited_slot_join(native_coefficient, operands, returned):
    # Different sides of a completed cancel-based identity need not be
    # structurally equal. Join the actual input and inherit its published zero.
    require(native_coefficient==operands['left'],'actual inherited slot input')
    require(returned['cancelled']=={'text':'0','srepr':'Integer(0)'},'published slot identity return')
    return {'left':operands['left'],'right':operands['right'],'returned':returned,'functionCalled':False}


def factor_address_join(address, operands, proof_address):
    require(proof_address['addressId']==operands['addressId'] and
            proof_address['fullFactorProof']['proof']==address['fullFactorProof']['proof'],'actual original proof address')
    # This legacy digest is a serialized operand-tuple identity, not a file hash.
    # Keep it attached to its original address; join actual operands below.
    require(proof_address['fullFactorProof']['operandSha256']==address['fullFactorProof']['operandSha256'],'same inherited proof identity')
    for a in (address,proof_address):
        require(a['responseMap']['id']==a['face']+'-'+str(tuple(a['responseGrade']))+'-'+a['component'],'native face/grade/component map label')
    require(proof_address['responseMap']==operands['actualResponseMap'],'original proof response map')
    # The same actual scalar map may serve both native faces. Its descriptive
    # ID differs; all mathematical operands and assumptions must be identical.
    require({k:v for k,v in address['responseMap'].items() if k!='id'}==
            {k:v for k,v in operands['actualResponseMap'].items() if k!='id'},'actual inherited response arguments')
    require(address['normalOriginal']==operands['addressNormalOriginal'],'actual inherited normal operand')
    require(address['fullFactorProof']['completeNormalMap']==operands['requiredMap'],'complete inherited depth/normal arguments')
    require(address['responseOriginal']==operands['actualResponseMap']['original'] and
            address['responseCoefficient']==operands['actualResponseMap']['mapped'],'actual response coefficient arguments')
    return {'proof':address['fullFactorProof']['proof'],'addressId':address['addressId'],'normalOriginal':address['normalOriginal'],
            'normalMultiplier':address['normalMultiplier'],'responseCoefficient':address['responseCoefficient'],'completeNormalMap':address['fullFactorProof']['completeNormalMap']}


def chunks(values,size=16):
    for i in range(0,len(values),size):yield i//size,values[i:i+size]


def checkpoint_batch(J,name,inputs,derive):
    """Preserve every completed scalar even if another scalar in its batch fails."""
    def work():
        done=[]
        try:
            for item in inputs:done.append(derive(item))
        except BaseException:
            J.emit(name+'-partial-returns',done)
            raise
        return done
    return J.stage(name,{'children':inputs,'oldFunctionsCalled':False},work)


class EvidenceLog:
    """Append-only operands/returns without rewriting an ever-growing index."""
    def __init__(self,path,encode):
        self.path=path;self.encode=encode;self.count=0;self.previous=None
        with path.open('x') as f:f.flush();os.fsync(f.fileno())
    def append(self,kind,value):
        payload={'sequence':self.count,'kind':kind,'previousSha256':self.previous,'value':self.encode(value)}
        digest=hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False).encode()).hexdigest()
        with self.path.open('a') as f:
            f.write(json.dumps({'payload':payload,'sha256':digest},allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
        self.count+=1;self.previous=digest
        return {'sequence':self.count-1,'sha256':digest}


def source_rules(c2, composition):
    tree=ast.parse(c2);assign={t.id:n.value for n in tree.body if isinstance(n,ast.Assign)
        for t in n.targets if isinstance(t,ast.Name)}
    waves=ast.literal_eval(assign['WAVE_NAMES'].args[0]);dimensions=ast.literal_eval(assign['DIMENSION_SCHEMA'])
    wave=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='wave_jet')
    source=ast.get_source_segment(c2,wave)
    require("name = 'theta_d' + name.rsplit('_', 1)[1]" in source,'actual grad_theta alias')
    require("if '_tt' in suffix:" in source and 'sp.diff(value, TIME, 2)' in source and 'sp.diff(value, TIME)' in source,'native time derivative rule')
    nodes=[n for n in ast.walk(ast.parse(composition)) if isinstance(n,ast.FunctionDef) and n.name=='wave']
    require(len(nodes)==1,'saved pressure wave mapping')
    phase=next(n.value for n in ast.walk(nodes[0]) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='out' for t in n.targets))
    expected=ast.parse("(-sp.I*3)**spec['timeOrder']",mode='eval').body
    require(ast.dump(phase)==ast.dump(expected),'pressure negative time carrier')
    at_source=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='at_source')
    clauses=[n for n in ast.walk(at_source) if isinstance(n,ast.If) and "w1_profile" in ast.unparse(n.test)]
    expected=ast.parse("if base in ('w1_profile','m1_profile'):\n    value*=self.values['L_W']**len(indices)").body[0]
    require(len(clauses)==1 and ast.dump(clauses[0])==ast.dump(expected),'native L profile factor')
    return {'waveNames':waves,'dimensions':dimensions,'waveSource':source,
            'profileRule':ast.get_source_segment(c2,clauses[0]),'pressureTimeRule':ast.unparse(phase)}


def run_science(manifest,J,ns,unused_exact_helper):
    raw={};copies={};used=set();D=ns['decode']
    audit=EvidenceLog(J.out/'exact-evidence.jsonl',J.encode)
    for alias,v in manifest['savedInputs'].items():
        src=Path(v['path']);dst=J.out/'saved'/alias;dst.parent.mkdir(parents=True,exist_ok=True)
        require(sha(src)==v['sha256'],'saved source '+alias);shutil.copyfile(src,dst)
        require(sha(dst)==v['sha256'],'saved copy '+alias)
        copies[alias]={'source':str(src),'path':str(dst.relative_to(J.out)),'sha256':v['sha256'],'bytes':dst.stat().st_size}
        raw[alias]=json.loads(dst.read_text())
    save(J.out/'saved-copy-index.json',copies)
    def take(alias):used.add(alias);return raw[alias]
    def load(alias):return D(take(alias))
    helpers={'sp':sp,'require':require,'re':re}
    exec(compile(definitions(Path(manifest['compositionSource']).read_text(),('jet_spec','polynomial_terms','quotient_recurrence')),'scalar-algebra-helpers','exec'),helpers)
    jet_spec=helpers['jet_spec'];terms=helpers['polynomial_terms'];quotient=helpers['quotient_recurrence']
    c2=Path(manifest['c2Source']).read_text();rules=source_rules(c2,Path(manifest['compositionSource']).read_text())
    J.emit('native-wave-profile-contract',rules)
    context=load('inventory/actual-binding-context.json')['restored'];physical=json.loads(Path(manifest['physicalInput']).read_text())
    require(context['physicalInput']==physical and context['frequency']==3,'actual saved physical input/frequency')
    eta,sigma=context['independentGrades'];eps=context['epsilon'];L=context['numeric']['L_W']
    require((eta.name,sigma.name,eps.name)==('eta_bg','sigma_W','epsilon_shape') and L==10,'native grades amplitude L')
    numeric={k:sp.Rational(v) for k,v in numeric_extension(physical['parameters'],context['numeric'],(eta.name,sigma.name,eps.name)).items()}
    require(all(numeric[k]==v for k,v in context['numeric'].items()),'all shared numeric bindings identical')
    require(all(s.name not in numeric for s in (eta,sigma,eps)) and 'c_s0' not in numeric,'no grade or speed binding')
    J.emit('extended-binding-context',{'saved':context,'physical':physical,'numeric':numeric,'additionalNames':sorted(set(numeric)-set(context['numeric'])),
        'frequencyOverride':{'old':physical['parameters']['omega'],'actual':3},'independentGrades':[eta,sigma],'epsilon':eps,'effectiveSpeedOnlyInPressure':True})
    T=sp.Symbol('full_weak_tanh_variable',real=True);x=sp.Symbol('full_weak_x',real=True)
    profiles={};profile_records=[];den_cache={};poly_cache={};inherited=[]
    def zero_record(evidence,label,left,right):
        audit.append('zero-input',{'name':label,'left':left,'right':right})
        raw_residual=left-right;residual=sp.cancel(sp.together(raw_residual))
        rec={'name':label,'left':left,'right':right,'raw':raw_residual,'cancelled':residual}
        audit.append('zero-return',rec);evidence.setdefault('identities',[]).append(rec)
        require(residual==0,label)
    def constant(value,allow_zero=True):
        audit.append('constant-input',{'value':value,'allowZero':allow_zero})
        require(not value.free_symbols,'constant coefficient has no free symbols')
        re_part,im_part=sp.expand_complex(value).as_real_imag();re_part=sp.cancel(re_part);im_part=sp.cancel(im_part)
        audit.append('constant-components',{'value':value,'real':re_part,'imaginary':im_part})
        require(re_part.is_Rational is True and im_part.is_Rational is True,'exact rational complex coefficient')
        require(re_part.is_finite is True and im_part.is_finite is True,'finite components')
        require(sp.cancel(value-re_part-sp.I*im_part)==0,'constant component reconstruction')
        require(allow_zero or re_part!=0 or im_part!=0,'nonzero constant domain')
        return {'value':value,'real':re_part,'imaginary':im_part,'finite':True,'nonzero':bool(re_part!=0 or im_part!=0)}
    def polynomial(expr):
        audit.append('polynomial-input',{'expression':expr,'variable':T})
        original=expr;num,den=sp.fraction(sp.cancel(expr));dc=constant(den,False)
        mapping=terms(num,T,sp.Symbol('full_weak_unused'))
        require(all(j==0 for i,j in mapping),'one polynomial variable')
        coefficients={i:sp.cancel(v/den) for (i,j),v in mapping.items()};certs=[constant(v) for v in coefficients.values()]
        rebuilt=sum(v*T**i for i,v in coefficients.items());residual=sp.cancel(original-rebuilt)
        audit.append('polynomial-return',{'original':original,'reconstruction':rebuilt,'residual':residual,'coefficients':[{'order':i,'value':v} for i,v in coefficients.items()]})
        require(residual==0,'polynomial reconstruction')
        return {'original':original,'polynomial':rebuilt,'denominator':dc,'coefficients':[{'order':i,'value':v} for i,v in sorted(coefficients.items())],
            'coefficientCertificates':certs,'degree':max([0]+[i for i,v in coefficients.items() if v!=0]),'zero':all(v==0 for v in coefficients.values()),
            'bound':sum(sp.Abs(v) for v in coefficients.values())}
    def profile_atom(atom):
        name=atom.name;match=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',name)
        require(match is not None,'profile atom syntax');base,suffix=match.groups();directions=re.findall(r'd([123])',suffix)
        if name not in profiles:
            initial=(1+T)/2 if base=='w' else (1-T**2)/3;poly=initial
            for _ in directions:poly=sp.expand((1-T**2)*sp.diff(poly,T)/L)
            unscaled=poly if all(d=='1' for d in directions) else sp.S.Zero;scaled=L**len(directions)*unscaled
            rec={'name':name,'directions':directions,'order':len(directions),'L':L,'unscaled':unscaled,'value':scaled,'initial':initial,
                 'transverseProfileDerivative':any(d!='1' for d in directions)}
            expected=sp.diff(initial.subs(T,sp.tanh(x/L)),x,len(directions)) if not rec['transverseProfileDerivative'] else sp.S.Zero
            zero_record(rec,'actual-profile-derivative',unscaled.subs(T,sp.tanh(x/L)),expected)
            if directions:
                zero_record(rec,'left-profile-jet',scaled.subs(T,-1),sp.S.Zero);zero_record(rec,'right-profile-jet',scaled.subs(T,1),sp.S.Zero)
            profiles[name]=rec;profile_records.append(rec)
        return profiles[name]
    def bind(expr):
        result=expr
        for s in result.free_symbols:
            if s.name in (eta.name,sigma.name,eps.name):require(s=={eta.name:eta,sigma.name:sigma,eps.name:eps}[s.name],'original symbolic grade assumptions')
        result=result.xreplace({s:context['densityMap'][s.name] for s in result.free_symbols if s.name in context['densityMap']})
        for _ in range(len(context['profileEqualities'])+2):
            new=result.xreplace({s:context['profileEqualities'][s.name] for s in result.free_symbols if s.name in context['profileEqualities']})
            if new==result:break
            result=new
        else:raise ValueError('background binding cycle')
        return result.xreplace({s:numeric[s.name] for s in result.free_symbols if s.name in numeric})
    # Restore prior pressure joins by their actual arguments; new local equations do not recompute them.
    for row in ROWS:
        for tag in ['pressure-source-join','pressure-bound-join','affine-pressure-reconstruction']:
            alias='inventory/'+row+'-'+tag
            result=take(alias+'-return.json');require(result['cancelled']=={'text':'0','srepr':'Integer(0)'},'inherited exact pressure identity')
            inherited.append({'name':alias,'input':copies[alias+'-input.json'],'return':copies[alias+'-return.json'],'functionCalled':False})
    # Completed pressure certificate and addresses are restored, not re-derived.
    conclusion=take('weak/analytic-conclusion.json');certs=take('weak/all-coefficient-certificates.json');fields=take('inventory/fields.json')
    prior_manifest=take('weak/input-manifest.json');prior_input=take('weak/global-weak-composition-input.json')
    prior_return=take('weak/global-weak-composition-return.json');prior_copies=take('weak/saved-copy-index.json')
    require(prior_input['manifestSha256']==copies['weak/input-manifest.json']['sha256'],'actual inherited run manifest')
    for k in ('frequency','anchoring','density','drain','edges','effectiveCs','sourceFrequencyContinued','responseRegularization'):
        require(prior_manifest['scope'][k]==manifest['scope'][k],'actual inherited pressure scope '+k)
    for alias in ['inventory/fields.json','inventory/whole-tag-definitions.json','inventory/typed-direct-objects.json','inventory/fourier-and-unit-provenance.json']+['inventory/'+r+'-ordered-addresses.json' for r in ROWS]:
        require(prior_copies[alias]['sha256']==copies[alias]['sha256'] and prior_copies[alias]['source']==copies[alias]['source'],'actual pressure computation operands '+alias)
    definitions_record=take('weak/whole-definition-certificate.json');duality=take('weak/new-weak-duality-and-order.json')
    definition_inputs=take('weak/whole-definition-inputs.json')
    require(definition_inputs['tags']==take('inventory/whole-tag-definitions.json') and
            definition_inputs['typedDirect']==take('inventory/typed-direct-objects.json'),'actual completed whole-definition arguments')
    require(definition_inputs['savedReference']==take('saved/reference/retained-response-census.json'),'actual reference used by completed whole-definition check')
    for tag,alias in [('H','saved/reference/left-height-subtracted-PV.json'),('Jwhole','saved/reference/right-height-PV-operands.json'),('Dwhole','saved/direct/closed-density.json')]:
        require(definition_inputs['tags'][tag]['savedDefinition']==take(alias) and
                definition_inputs['tags'][tag]['sha256']==copies[alias]['sha256']==prior_copies[alias]['sha256'],'actual inherited whole definition '+tag)
    require(definitions_record=={'savedDefinitionsJoined':['H','Jwhole','Dwhole'],'typedFactorsDistinct':True,'wholeDirectOnce':True,'nativeIterationOnce':True},'actual whole-definition certificate')
    require(duality['inheritedConvention']==take('inventory/fourier-and-unit-provenance.json')['contract'] and duality['nativeEpsilonOnce'] is True and duality['sourceDerivativeBeforeMultiplication'] is True,'inherited Fourier/epsilon/order operands')
    require(duality['sourceDerivativeVariable']=='original p' and duality['normalVariable']=='response output l','actual pressure derivative placement')
    for fid,field in fields.items():
        original=take('weak/field-'+fid+'-operands.json');derivative=take('weak/field-'+fid+'-derivative-class.json')
        recon=take('weak/field-'+fid+'-reconstruction-input.json');ret=take('weak/field-'+fid+'-reconstruction-return.json')
        require(original['fieldId']==fid and original['original']==field['field']==recon['left'] and derivative['P0']==certs[fid]['polynomial'],'inherited actual field coefficient and bound')
        field_map={'fieldId':fid,'savedArguments':original,'savedDerivative':derivative,'savedReconstruction':recon}
        zero_record(field_map,'new-field-argument-'+fid,D(recon['right']),D(derivative['P0']).xreplace({D(original['T']):sp.tanh(D(original['x'])/L)}))
        require(ret['cancelled']=={'text':'0','srepr':'Integer(0)'},'restored field reconstruction return')
        inherited.append({'name':'field-'+fid+'-reconstruction','functionCalled':False,'input':copies['weak/field-'+fid+'-reconstruction-input.json'],'return':copies['weak/field-'+fid+'-reconstruction-return.json']})
    require(conclusion['status']=='SOURCE_JOINED_GLOBAL_WEAK_PRESSURE_CERTIFICATE' and conclusion['testSpace']=='S(R) x S(R), complex bilinear' and conclusion['localSlabPartIncluded'] is False,'actual inherited pressure scope')
    require(conclusion['wholeDirectOnce'] and conclusion['nativeIterationOnce'] and not conclusion['hiddenProjection'],'inherited pressure multiplicity')
    require(set(certs)==set(fields) and len(fields)==34,'actual inherited coefficient set')
    weak_routes=take('weak/weak-address-coverage.json');byid={a['addressId']:a for a in weak_routes};require(len(byid)==len(weak_routes)==13260,'all weak pressure addresses')
    native_addresses={a['addressId']:a for row in ROWS for a in take('inventory/'+row+'-ordered-addresses.json')}
    require(len(native_addresses)==13260,'all native pressure argument addresses')
    aggregate=[];allids=[];factor_joins={};wave_joins={};child_joins=[]
    for row in ROWS:
        oldrow=load('consumer/'+row+'-consumer-input.json')
        require(oldrow['raw']==load('inventory/'+row+'-pressure-source-join-input.json')['right'] and oldrow['bound']==load('inventory/'+row+'-affine-pressure-reconstruction-input.json')['left'],'actual saved native consumer source')
        affine=load('inventory/'+row+'-affine-pressure-reconstruction-input.json')
        native_slots={s.name:s for s in affine['right'].free_symbols if s.name in SLOTS}
        require(len(native_slots)==len([s for s in affine['right'].free_symbols if s.name in SLOTS]),'unambiguous native pressure symbols')
        require(all(v==0 or slot in native_slots for slot,v in oldrow['slotCoefficients'].items()),'nonzero coefficient has actual native slot')
        affine_sum=sum(oldrow['slotCoefficients'][slot]*native_slots.get(slot,sp.Symbol(slot)) for slot in SLOTS)
        new_join={'row':row,'savedBound':oldrow['bound'],'savedAffine':affine,'slotCoefficients':oldrow['slotCoefficients'],
            'slotSymbols':native_slots,'newSum':affine_sum,'newAssemblyJoin':True}
        zero_record(new_join,'new-'+row+'-affine-slot-sum',affine['right'],affine_sum)
        zero_record(new_join,'new-'+row+'-bound-slot-sum',oldrow['bound'],affine_sum)
        J.emit(row+'-new-pressure-slot-assembly',new_join)
        slots={}
        for slot in SLOTS:
            split=load('inventory/'+row+'-'+slot+'-split.json');full=load('inventory/'+row+'-'+slot+'-full-coefficient-input.json')
            r=take('inventory/'+row+'-'+slot+'-full-coefficient-return.json')
            inherited_slot_join(oldrow['slotCoefficients'][slot],full,r)
            inherited.append({'name':row+'-'+slot+'-full-coefficient','functionCalled':False,'input':copies['inventory/'+row+'-'+slot+'-full-coefficient-input.json'],'return':copies['inventory/'+row+'-'+slot+'-full-coefficient-return.json']})
            slots[slot]=split['retained']
        addresses=take('inventory/'+row+'-ordered-addresses.json')
        child_join=pressure_child_coverage(take('inventory/'+row+'-full-native-partition.json'),addresses,row)
        J.emit(row+'-addressed-native-pressure-children',child_join);child_joins.append(child_join)
        for a in addresses:
            wr=byid[a['addressId']];slot=('delta_p_' if a['slot']=='pressure' else 'd_w_delta_p_')+a['face']
            require(a['row']==row and wr['face']==a['face'] and wr['slot']==a['slot'] and wr['component']==a['component'],'same native row/face/slot pressure address')
            require(wr['sourceJet']==a['jet'] and wr['gradeTriple']==[a[k] for k in ('consumerGrade','responseGrade','sourceGrade')],'same wave/grade route')
            require(wr['sourceFieldId']==a['sourceTransform']['coefficientId'] and wr['consumerFieldId']==a['consumerTransform']['coefficientId'],'actual field IDs')
            require(fields[wr['sourceFieldId']]['field']==a['sourceField'] and fields[wr['consumerFieldId']]['field']==a['consumerField'],'actual field operands')
            require(D(a['consumerOriginal'])==slots[slot][str(tuple(a['consumerGrade']))],'exact actual native grade slot operand')
            require(wr['status']==a['status'] and wr['inheritedFactorProof']==a['fullFactorProof']['proof'],'same saved scope and factor proof')
            proof=a['fullFactorProof']['proof'];pa='inventory/factors/'+proof
            operands=take(pa+'-operands.json');fj=factor_address_join(a,operands,native_addresses[operands['addressId']])
            if proof not in factor_joins:
                old_factor=load(pa+'-full-mapped-residual-input.json');old_normal=load(pa+'-normal-source-join-input.json')
                for tail in ['full-mapped-residual','normal-source-join']:
                    ret=take(pa+'-'+tail+'-return.json');require(ret['cancelled']=={'text':'0','srepr':'Integer(0)'},'inherited factor/normal zero return')
                    inherited.append({'name':proof+'-'+tail,'functionCalled':False,'input':copies[pa+'-'+tail+'-input.json'],'return':copies[pa+'-'+tail+'-return.json']})
                require(old_normal['left']==D(a['normalOriginal']) and old_normal['right']==D(operands['savedNormal']),'actual saved native normal identity arguments')
                require(old_factor['right']==D(operands['mappedAddressFactor']),'actual saved complete factor argument')
                mapped_normal=D(a['normalOriginal']).xreplace(dict(D(a['fullFactorProof']['completeNormalMap'])))
                zero_record(fj,'new-'+proof+'-normal-argument',D(a['normalMultiplier']),mapped_normal)
                zero_record(fj,'new-'+proof+'-factor-argument',old_factor['right'],D(a['normalMultiplier'])*D(a['responseCoefficient']))
                factor_joins[proof]=fj
            else:
                require(all(fj[k]==factor_joins[proof][k] for k in ['normalOriginal','normalMultiplier','responseCoefficient','completeNormalMap']),'identical factor reuse arguments')
            jet=a['jet'];require(jet==jet_spec(jet['name']) and jet['name'] in rules['waveNames'],'actual native pressure wave jet')
            if jet['name'] not in wave_joins:
                p=sp.Symbol('composition_p',real=True);expected_wave=(-sp.I*3)**jet['timeOrder']
                for n,mom in zip(jet['spatialOrders'],(p,sp.Rational(1,5),sp.Rational(1,10))):expected_wave*=(sp.I*mom)**n
                wj={'jet':jet,'savedMultiplier':D(a['waveMultiplier']),'sourceMomentum':p,'nativeExpected':expected_wave}
                zero_record(wj,'new-pressure-wave-argument-'+jet['name'],wj['savedMultiplier'],expected_wave);wave_joins[jet['name']]=wj
            else:require(D(a['waveMultiplier'])==wave_joins[jet['name']]['savedMultiplier'],'identical wave argument reuse')
            require(a['epsilon']==take('inventory/actual-binding-context.json')['restored']['epsilon'] and a['epsilonCount'] in (0,1),'actual inherited native amplitude')
            require(a['epsilonCount']==(0 if a['status'].startswith('EXACT_ZERO_') else 1),'actual addressed zero/epsilon convention')
            tags=wr['wholeTagCheck'];require(len(tags)==2 and {t['operand'] for t in tags}=={'responseOriginal','responseCoefficient'},'both inherited tag operands')
            # The native mixed response is the sum of H and Jwhole, one each;
            # the separate direct response contains Dwhole only once.
            expected_tag_count={'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL':1,'NATIVE_MIXED_ITERATION':2}.get(a['component'],0)
            require(all(t['tagCount']==expected_tag_count for t in tags),'inherited whole-tag multiplicities')
            require(a['responseMap']['frequency']==3 and a['responseMap']['positiveRegulatorContinuation'] is False,'real source composition')
            allids.append(a['addressId'])
        aggregate.append({'row':row,'pressureChildIndices':take('inventory/'+row+'-full-native-partition.json')['pressureChildIndices'],'addressIds':[a['addressId'] for a in addresses],
            'actualConsumerSlotsJoined':list(slots),'pressureKernelFunctionsCalled':False})
    require(sorted(allids)==list(range(13260)),'once-only pressure addresses')
    require(sum(len(v['coveredHashes']) for v in child_joins)==12,'all twelve pressure child hashes covered')
    J.emit('new-pressure-factor-arguments',factor_joins);J.emit('new-pressure-wave-arguments',wave_joins)
    J.emit('pressure-source-assembly-joins',aggregate);J.emit('restored-pressure-identities',inherited)
    J.emit('inherited-pressure-law',{'conclusion':conclusion,'duality':duality,'wholeDefinitions':take('inventory/whole-tag-definitions.json'),
        'typedObjects':take('inventory/typed-direct-objects.json'),'definitionCertificate':definitions_record,'operationReturn':prior_return,
        'wholeKernelEnvelopes':take('weak/global-whole-kernel-envelopes.json'),'oldFunctionsCalled':False})
    source_origin=take('inventory/'+ROWS[0]+'-full-native-partition.json')['source']
    native_path=Path(source_origin['source']);require(sha(native_path)==manifest['sourcePins'][str(native_path)],'native export unchanged')
    with native_path.open('rb') as f:
        native_line=next(line for i,line in enumerate(f,1) if i==source_origin['valueLine'])
    require(hashlib.sha256(native_line).hexdigest()==source_origin['sourceLineSha256'],'native original source line')
    pressure_joins=[];cells={};child_results=[];new_jets={};raw_names=set();raw_denominators={};counts={}
    for row in ROWS:
        partition=take('inventory/'+row+'-full-native-partition.json');counts[row]=validate_partition(partition)
        require(partition['source']==source_origin,'all native row provenance joins')
        J.emit(row+'-native-partition-join',counts[row])
        selected=[e for e in partition['children'] if e['childIndex'] in partition['localChildIndices']]
        # The pressure operands are restored only for this new combined-row argument identity.
        old_source=load('inventory/'+row+'-pressure-source-join-input.json')
        old_bound=load('inventory/'+row+'-pressure-bound-join-input.json')
        old_affine=load('inventory/'+row+'-affine-pressure-reconstruction-input.json')
        source_pressure=sp.Add(*(sp.sympify(e['constructorText']) for e in partition['children'] if e['childIndex'] in partition['pressureChildIndices']))
        join={'row':row,'pressureIndices':partition['pressureChildIndices'],'savedSource':old_source,'savedBinding':old_bound,'savedAffine':old_affine,'nativeChildSum':source_pressure,'oldZeroFunctionsCalled':False}
        zero_record(join,'new-full-row-pressure-argument',source_pressure,old_source['left'])
        old_consumer=load('consumer/'+row+'-consumer-input.json')
        join['savedConsumer']=old_consumer;audit.append('pressure-assembly-input',join)
        require(old_source['right']==old_consumer['raw'] and old_bound['left']==old_consumer['bound'] and old_bound['right']==old_affine['left'],'old source/binding/affine operand ancestry')
        pressure_joins.append(join);J.emit(row+'-inherited-pressure-arguments',join)
        def derive(entry):
            ev={'row':row,'childIndex':entry['childIndex'],'sourceSha256':entry['sha256'],'sourceConstructor':entry['constructorText']}
            try:
                original=sp.sympify(entry['constructorText']);ev['original']=original
                names={s.name:s for s in original.free_symbols};require(len(names)==len(original.free_symbols),'unambiguous local symbol assumptions')
                raw_names.update(names)
                require(not any(n.lower().startswith(('c_s','cs_')) or n.lower()=='cs' or 'speed' in n.lower() for n in names),'raw local speed independence')
                waves=[s for s in original.free_symbols if s.name in rules['waveNames']]
                require(len(waves)==1,'one linear wave jet in native local child');wave=waves[0];spec=jet_spec(wave.name);require(spec is not None,'actual native wave syntax')
                new_jets[wave.name]=spec
                allowed=set(numeric)|{eta.name,sigma.name,eps.name}|set(rules['waveNames'])|set(context['densityMap'])|set(context['profileEqualities'])
                require(all(n in allowed or re.fullmatch(r'[wm]1_profile((?:_?d[123])*)',n) for n in names),'complete native local symbol classification')
                raw_coefficient=sp.cancel(original/(eps*wave));ev['nativeCoefficient']=raw_coefficient
                require(not raw_coefficient.has(eps,wave) and not any(s.name in rules['waveNames'] for s in raw_coefficient.free_symbols),'native epsilon and wave linearity')
                zero_record(ev,'native-child-reconstruction',original,eps*wave*raw_coefficient)
                for power in original.atoms(sp.Pow):
                    if power.exp.is_negative is True:raw_denominators.setdefault(sp.srepr(power.base),power.base)
                bound=bind(raw_coefficient);ev['boundCoefficient']=bound
                num,den=sp.fraction(sp.cancel(bound));nt=terms(num,eta,sigma);dt=terms(den,eta,sigma)
                den0=dt.get((0,0),sp.S.Zero);dc=constant(den0,False)
                table=quotient(nt,dt,G,sp.cancel);retained=sum(table[g]*eta**g[0]*sigma**g[1] for g in G)
                remainder=num-den*retained;rt=terms(remainder,eta,sigma)
                ev['grade']={'numerator':num,'denominator':den,'zeroDenominator':dc,'table':[{'grade':g,'coefficient':table[g]} for g in G],
                    'unprojected':bound,'excludedRemainder':bound-retained,'quotientNumeratorRemainder':remainder}
                for g in G:zero_record(ev,'quotient-'+str(g),rt.get(g,sp.S.Zero),sp.S.Zero)
                multiplier=(-sp.I*3)**spec['timeOrder']*(sp.I/5)**spec['spatialOrders'][1]*(sp.I/10)**spec['spatialOrders'][2]
                ev['jet']=spec;ev['waveMultiplier']=multiplier;ev['fieldColumn']=FIELDS.index(spec['channel']);ev['xOrder']=spec['spatialOrders'][0]
                ev['nativeJetDimension']=rules['dimensions'][wave.name];ev['dimensionScope']='Native wave/row units inherited; coefficient dimensions after numeric binding required, not independently reconstructed'
                mapped=[];maps=[]
                for g in G:
                    coefficient=table[g];mapping={}
                    for atom in coefficient.free_symbols:
                        pr=profile_atom(atom);mapping[atom]=pr['value'];maps.append({'symbol':atom,'profileRule':pr['name'],'value':pr['value']})
                    value=sp.cancel(coefficient.xreplace(mapping)*multiplier);cert=polynomial(value)
                    mapped.append({'grade':g,'value':value,'polynomial':cert,'profileMap':[[a,b] for a,b in mapping.items()]})
                ev['mappedGrades']=mapped;ev['profileMaps']=maps;ev['completed']=True
                return ev
            except BaseException:
                ev['traceback']=traceback.format_exc();ev['completed']=False;J.emit(row+'-failed-local-child-'+str(entry['childIndex']),ev);raise
        row_results=[]
        for batch,entries in chunks(selected):row_results.extend(checkpoint_batch(J,row+'-local-batch-'+str(batch).zfill(3),entries,derive))
        for ev in row_results:
            child_results.append(ev)
            for record in ev['mappedGrades']:
                key=(row,ev['fieldColumn'],ev['xOrder'],tuple(record['grade']));cells.setdefault(key,[]).append((ev['childIndex'],record['value']))
        # Full native row = new local child sum + restored pressure sum, matched to all native children.
        local_sum=sum(ev['original'] for ev in row_results)
        # Reuse the just-interpreted local children and restored pressure operand.
        # The full constructor's exhaustive child membership was checked as AST,
        # without deserializing any old slab producer or redoing pressure binding.
        full=sp.Add(*(ev['original'] for ev in row_results),source_pressure)
        ev={'row':row,'actualFullRow':full,'newLocalSum':local_sum,'inheritedPressure':old_source['left'],'localIndices':partition['localChildIndices'],'pressureIndices':partition['pressureChildIndices']}
        zero_record(ev,'new-full-native-row-reconstruction',full,local_sum+old_source['left'])
        slotmap={s:sp.S.Zero for s in full.free_symbols if s.name in SLOTS}
        zero_record(ev,'new-local-slot-zero-reconstruction',full.xreplace(slotmap),local_sum)
        J.emit(row+'-full-row-reconstruction',ev)
    require(sum(v['local'] for v in counts.values())==len(child_results)==2908 and sum(v['pressure'] for v in counts.values())==12,'complete actual local coverage')
    J.emit('native-raw-symbol-census',{'names':sorted(raw_names),'speedAbsent':True,'additionalNumeric':sorted(set(numeric)-set(context['numeric'])),'jets':new_jets})
    J.emit('profile-jet-certificates',profile_records)
    # Bind and certify actual original denominator bases. Memory phases retain the native sign.
    domain=[]
    for i,den in enumerate(raw_denominators.values()):
        domain.append(memory_domain_certificate(den,bind,constant,terms,zero_record,audit.append,sp))
    J.emit('native-denominator-and-time-sign-joins',domain)
    cell_records=[];max_order=max(k[2] for k in cells)
    for row,column,n,g in itertools.product(ROWS,range(5),range(max_order+1),G):
        key=(row,column,n,g);entries=cells.get(key,[]);value=sp.cancel(sum((v for _,v in entries),sp.S.Zero));cert=polynomial(value)
        first=sp.expand((1-T**2)*sp.diff(cert['polynomial'],T)/L)
        rec={'row':row,'field':FIELDS[column],'fieldColumn':column,'xOrder':n,'grade':g,'sourceChildren':[i for i,_ in entries],
            'summands':[v for _,v in entries],'coefficient':value,'polynomial':cert,'firstDerivativePolynomial':first,
            'leftEndpoint':sp.cancel(value.subs(T,-1)),'rightEndpoint':sp.cancel(value.subs(T,1)),
            'endpointScope':'Local coefficients only, not full constant-height pressure or a uniform end pencil','epsilonPower':0 if cert['zero'] else 1}
        zero_record(rec,'cell-source-sum',value,sum((v for _,v in entries),sp.S.Zero))
        zero_record(rec,'cell-physical-polynomial',value.subs(T,sp.tanh(x/L)),cert['polynomial'].subs(T,sp.tanh(x/L)))
        zero_record(rec,'cell-derivative',sp.diff(value.subs(T,sp.tanh(x/L)),x),first.subs(T,sp.tanh(x/L)))
        cell_records.append(rec)
    J.emit('all-local-cells',cell_records)
    # Reconstruct every retained physical row with independent test jets.
    formal_jets={(col,n):sp.Symbol('local_field_'+str(col)+'_jet_'+str(n)) for col in range(5) for n in range(max_order+1)}
    for row in ROWS:
        for g in G:
            left=sum(rec['value']*formal_jets[ev['fieldColumn'],ev['xOrder']] for ev in child_results if ev['row']==row for rec in ev['mappedGrades'] if tuple(rec['grade'])==g)
            right=sum(rec['coefficient']*formal_jets[rec['fieldColumn'],rec['xOrder']] for rec in cell_records if rec['row']==row and tuple(rec['grade'])==g)
            J.zero(row+'-new-physical-row-grade-'+str(g[0])+str(g[1]),left,right)
    # New addressed local controls: coefficients, never response or field values.
    controls=[]
    mixed=next((ev for ev in child_results if any(tuple(r['grade'])==(1,1) and not r['polynomial']['zero'] for r in ev['mappedGrades'])),None)
    require(mixed is not None,'applicable actual mixed child')
    mg=next(r for r in mixed['mappedGrades'] if tuple(r['grade'])==(1,1))
    cell=next(c for c in cell_records if (c['row'],c['fieldColumn'],c['xOrder'],tuple(c['grade']))==(mixed['row'],mixed['fieldColumn'],mixed['xOrder'],(1,1)))
    movement=polynomial(-mg['value']);audit.append('mixed-control',{'baseline':cell['coefficient'],'omitted':mg['value'],'movement':movement});require(not movement['zero'],'mixed omission responds')
    controls.append({'name':'omit-native-mixed-child','row':mixed['row'],'childIndex':mixed['childIndex'],'cell':{k:cell[k] for k in ['fieldColumn','xOrder','grade']},
        'baseline':cell['coefficient'],'corrupt':sp.cancel(cell['coefficient']-mg['value']),'movement':movement,'formalCoefficientOnly':True})
    candidates=[];chosen=None
    for ev in child_results:
        for rec in ev['mappedGrades']:
            if chosen is not None:break
            before=next(v['coefficient'] for v in ev['grade']['table'] if tuple(v['grade'])==tuple(rec['grade']))
            for atom,value in rec['profileMap']:
                pr=profiles[atom.name]
                if not pr['order'] or pr['transverseProfileDerivative']:continue
                mapping=dict(rec['profileMap']);mapping[atom]=pr['unscaled'];wrong=sp.cancel(before.xreplace(mapping)*ev['waveMultiplier']);move=polynomial(wrong-rec['value'])
                candidates.append({'row':ev['row'],'childIndex':ev['childIndex'],'grade':rec['grade'],'profile':atom,'zero':move['zero']})
                if not move['zero']:
                    cell=next(c for c in cell_records if (c['row'],c['fieldColumn'],c['xOrder'],tuple(c['grade']))==(ev['row'],ev['fieldColumn'],ev['xOrder'],tuple(rec['grade'])))
                    chosen={'name':'omit-native-profile-L','row':ev['row'],'childIndex':ev['childIndex'],'grade':rec['grade'],'actualProfileRule':pr,
                        'baselineChild':rec['value'],'corruptChild':wrong,'baselineCell':cell['coefficient'],'corruptCell':sp.cancel(cell['coefficient']+wrong-rec['value']),
                        'movement':move,'formalCoefficientOnly':True};break
        if chosen is not None:break
    J.emit('profile-control-candidates',candidates);require(chosen is not None,'applicable nonzero native L omission');controls.append(chosen)
    leibniz=next((c for c in cell_records if c['xOrder']>0 and c['polynomial']['degree']>0 and not c['polynomial']['zero']),None)
    if leibniz is None:J.emit('Leibniz-control-unavailable',{'reason':'No surviving nonconstant coefficient with positive x derivative; no synthetic replacement'})
    else:
        n=leibniz['xOrder'];coeff=leibniz['coefficient'];ds=coeff;movements=[]
        for j in range(1,n+1):
            ds=sp.expand((1-T**2)*sp.diff(ds,T)/L);movements.append({'trialDerivative':n-j,'coefficient':sp.binomial(n,j)*ds,'certificate':polynomial(sp.binomial(n,j)*ds)})
        audit.append('Leibniz-control-coefficients',{'coefficient':coeff,'order':n,'movements':movements})
        require(any(not m['certificate']['zero'] for m in movements),'actual Leibniz movement')
        U=sp.Function('full_weak_trial')(x);a=coeff.subs(T,sp.tanh(x/L));ev={'name':'local-Leibniz-interchange','cell':{k:leibniz[k] for k in ['row','fieldColumn','xOrder','grade']},'coefficient':coeff,'movements':movements,'formalCoefficientOnly':True}
        zero_record(ev,'actual-Leibniz-identity',sp.diff(a*U,x,n)-a*sp.diff(U,x,n),sum(m['coefficient'].subs(T,sp.tanh(x/L))*sp.diff(U,x,m['trialDerivative']) for m in movements))
        controls.append(ev)
    J.emit('new-local-responsive-controls',controls)
    J.emit('complete-weak-assembly',{'status':'SOURCE_JOINED_FULL_RETAINED_WEAK_OPERATOR','testSpace':'S(R)^5 x S(R)^5, complex bilinear','localSourceChildren':2908,
        'pressureSourceChildren':12,'localCells':len(cell_records),'pressureAddresses':len(allids),'maximumLocalXOrder':max_order,'nativeEpsilonOnce':True,'independentGrades':G,
        'pressureConclusionInherited':conclusion,'sourceRows':ROWS,'fields':FIELDS,'wholeDirectOnce':True,'nativeIterationOnce':True,'oldReducedLocalObjectUsed':False,
        'analyticArgument':'Finite sum of local continuous Schwartz pairings (bound 2 C_a p_2,0(v) p_0,n(u)) and the inherited continuous pressure form; not machine measure theory.',
        'derivativeRecurrence':'P[n+1]=(1-T^2)*P[n]\u2032/10','localEndpointsOnly':True,'physicalPositiveRegulatorComposition':False,'planeWaveScattering':False,'finiteInverse':False,'loss':False,'calibration':False,'drain':False})
    J.emit('consumed-source-index',{'files':sorted(used),'copies':copies,'oldFunctionsCalled':False})
    J.emit('exact-evidence-receipt',{'path':audit.path.name,'records':audit.count,'lastSha256':audit.previous,'sha256':sha(audit.path),'bytes':audit.path.stat().st_size})
    return {'executionStatus':'COMPLETED_FULL_RETAINED_WEAK_OPERATOR_CERTIFICATES','localChildren':len(child_results),'localCells':len(cell_records),'localWaveJets':len(new_jets),
        'maximumLocalXOrder':max_order,'pressureAddressesRestored':len(allids),'priorZeroReturnsRestored':len(inherited),'newLocalControls':len(controls),
        'responseIntegralsEvaluated':False,'finiteSolves':0,'weakOperatorOnly':True,'scienceAcceptance':False,'productionChanges':False}


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
        result=J.stage('native-local-full-weak',{'manifestSha256':gate['manifestSha256'],'reviewSha256':gate['buildReviewRecordSha256']},lambda:run_science(manifest,J,ns,exact['exact_nonzero_number']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        if (args.out/'exact-evidence.jsonl').exists():
            evidence=args.out/'exact-evidence.jsonl'
            save(args.out/'exact-evidence-final-receipt.json',{'path':evidence.name,'sha256':sha(evidence),'bytes':evidence.stat().st_size,'mayBeIncomplete':code!=0})
        if (args.out/'saved-copy-index.json').exists():
            for v in json.loads((args.out/'saved-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(v['expected']!=v['actual'] for v in records.values()):result['integrityFailure']=True;result['executionStatus']='INTEGRITY_FAILURE_PRESERVED';code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False);save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':sys.exit(main())
