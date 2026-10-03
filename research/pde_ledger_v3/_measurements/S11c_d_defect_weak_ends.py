#!/usr/bin/env python3
"""Translated Schwartz weak ends; no integral, root solve, matrix or inverse."""
import argparse
import ast
import hashlib
import itertools
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback
ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','posthash_records','containment','Journal','decode','one_symbol','expanded_sinh_arguments')
G=((0,0),(1,0),(0,1),(1,1))
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
ZERO={'text':'0','srepr':'Integer(0)'}


def require(v, message):
    if v is not True:
        raise ValueError(message)

def first_shape_grade_parts(first, eta, sigma, expand):
    """Select independent linear grades; callers save and check reconstruction."""
    expanded=expand(first)
    height_coefficient=expanded.coeff(eta,1).coeff(sigma,0)
    slope_coefficient=expanded.coeff(sigma,1).coeff(eta,0)
    return {'expandedOriginal':expanded,'heightCoefficient':height_coefficient,
            'slopeCoefficient':slope_coefficient,'height':eta*height_coefficient,
            'slope':sigma*slope_coefficient}

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
    require(gate['status']=='READY_FOR_ONE_WEAK_ENDS_INSTRUMENT','gate status')
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

def verify_build_assessment(gate,review):
    require(review['methodAssessed'] is True and review['allChecksPassed'] is True and
            review['independentBuildClearance'] is True and gate['independentBuildClearance'] is True,
            'actual paired build assessment')
    for engine in ('claude','grok'):
        require(review['reports'][engine]['literalVerdict']=='CLEAR FOR THIS TRANSLATED WEAK-END BUILD',
                'literal build verdict')
    for key in ('workerSha256','manifestSha256','launcherSha256','sharedGuardSha256','supervisorSha256'):
        require(review[key]==gate[key],'review/gate '+key)


def inherit_zero(record):
    require(record['cancelled']==ZERO,'published exact zero return')
    return {'returned':record,'functionCalled':False}


def cell_key(cell):
    require(cell['row'] in ROWS and cell['field'] in FIELDS and
            cell['fieldColumn']==FIELDS.index(cell['field']) and
            type(cell['xOrder']) is int and 0<=cell['xOrder']<=3 and
            tuple(cell['grade']) in G,'actual local cell route')
    return cell['row'],cell['field'],tuple(cell['grade']),cell['xOrder']


def address_metadata_join(a,coverage,fields,factors,waves):
    require(a['row'] in ROWS and a['face'] in ('plus','minus') and
            a['slot'] in ('pressure','normal') and a['jet']['channel'] in FIELDS,'address route')
    grades=[tuple(a[k]) for k in ('consumerGrade','responseGrade','sourceGrade')]
    require(all(g in G for g in grades) and tuple(a['targetGrade']) in G and
            tuple(map(sum,zip(*grades)))==tuple(a['targetGrade']),'retained ordered grade sum')
    for key in ('face','slot','component','status'):require(a[key]==coverage[key],'prior addressed '+key)
    require(a['addressId']==coverage['addressId'] and coverage['gradeTriple']==[list(g) for g in grades]
            and a['jet']==coverage['sourceJet'],'actual coverage arguments')
    for role in ('source','consumer'):
        fid=a[role+'Transform']['coefficientId']
        require(fid==coverage[role+'FieldId'] and a[role+'Field']==fields[fid]['field'],'actual '+role+' field')
    require(a['epsilonCount']==(1 if a['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED' else 0),'zero or once epsilon')
    if a['status']=='EXACT_ZERO_CONSUMER':require(a['consumerField']==ZERO,'original zero consumer')
    elif a['status']=='EXACT_ZERO_SOURCE_JET':require(a['sourceField']==ZERO,'original zero source')
    else:require(a['status']=='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED','known formal address status')
    f=factors[a['fullFactorProof']['proof']];w=waves[a['jet']['name']]
    for key in ('normalOriginal','normalMultiplier','responseCoefficient'):
        require(a[key]==f[key],'inherited factor '+key)
    require(a['fullFactorProof']['completeNormalMap']==f['completeNormalMap'],'complete normal map')
    require(a['jet']==w['jet'] and a['waveMultiplier']==w['savedMultiplier'],'actual source wave multiplier')
    for rec in f['identities']+w['identities']:inherit_zero(rec)
    require(a['responseMap']['original']==a['responseOriginal'] and
            a['responseMap']['mapped']==a['responseCoefficient'] and
            a['responseMap']['frequency']==3 and a['responseMap']['positiveRegulatorContinuation'] is False,
            'actual real-frequency response map')
    require(a['responseMap']['id']==a['face']+'-'+str(tuple(a['responseGrade']))+'-'+a['component'],
            'actual native face/grade/component')
    if grades[1]==(0,0):
        require(a['component']=='NATIVE_FLAT' and a['responseMap']['flatSupport']=='k=l' and
                a['responseInputDepth']=='q(l)' and 'k=l' in a['deltaSupport'],'actual flat pole support')
    elif grades[1]==(1,0):require(a['component']=='NATIVE_HEIGHT','actual height component')
    elif grades[1]==(0,1):require(a['component']=='NATIVE_SLOPE','actual slope component')
    else:require(a['component'] in ('NATIVE_MIXED_ITERATION','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL'),
                 'separate whole mixed components')
    return {'addressId':a['addressId'],'factorProof':a['fullFactorProof']['proof'],
            'inheritedZeroReturns':[x['name'] for x in f['identities']+w['identities']],
            'oldFunctionsCalled':False}


def literal_field_joins(fid,fields,certificates,operands,polynomial,derivatives,proof,returned):
    require(operands['fieldId']==fid and operands['original']==fields[fid]['field'], 'actual coefficient field')
    require(certificates[fid]['polynomial']==derivatives['P0'],
            'actual polynomial and bounded-derivative operands')
    require(operands['profileScale']['physicalLength']=='10' and
            operands['profileScale']['declaredLength']==10,'actual pressure profile scale')
    require(proof['left']==operands['original'],'actual field inherited proof argument')
    inherit_zero(returned)
    return {'fieldId':fid,'operands':operands,'polynomial':polynomial,'derivativeClass':derivatives,
            'oldProofInput':proof,'oldProofReturn':returned,'functionCalled':False}


def source_statement_join(source,statement):
    candidates=[n for n in ast.walk(ast.parse(source)) if isinstance(n,(ast.Assign,ast.AnnAssign))]
    node=ast.parse(statement).body[0]
    require(any(ast.dump(n,include_attributes=False)==ast.dump(node,include_attributes=False) for n in candidates),
            'saved native assignment in actual source')
    return {'statement':statement,'actualSourceSha256':hashlib.sha256(source.encode()).hexdigest(),'functionCalled':False}


def native_phase_exponent(expression, bindings, imaginary_unit):
    """Bind only the actual native exp/arithmetic/sum(zip()) grammar; never eval."""
    root=ast.parse(expression,mode='eval').body
    def dotted(node,name):
        return isinstance(node,ast.Attribute) and isinstance(node.value,ast.Name) and node.value.id=='sp' and node.attr==name
    require(isinstance(root,ast.Call) and dotted(root.func,'exp') and len(root.args)==1 and not root.keywords,
            'native scalar exponential grammar')
    def scalar(node,env):
        if isinstance(node,ast.Name):
            require(node.id in env,'bound native phase name');return env[node.id]
        if isinstance(node,ast.Constant):
            require(type(node.value) is int,'integer native phase literal');return node.value
        if dotted(node,'I'):return imaginary_unit
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -scalar(node.operand,env)
        if isinstance(node,ast.BinOp):
            left=scalar(node.left,env);right=scalar(node.right,env)
            if isinstance(node.op,ast.Add):return left+right
            if isinstance(node.op,ast.Sub):return left-right
            if isinstance(node.op,ast.Mult):return left*right
            raise ValueError('unsupported native phase arithmetic')
        require(isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='sum'
                and len(node.args)==1 and not node.keywords,'native phase sum')
        generator=node.args[0]
        require(isinstance(generator,ast.GeneratorExp) and len(generator.generators)==1,'native phase generator')
        clause=generator.generators[0];iterator=clause.iter
        require(isinstance(clause.target,ast.Tuple) and len(clause.target.elts) in (2,3) and
                all(isinstance(n,ast.Name) for n in clause.target.elts) and not clause.ifs and not clause.is_async,
                'native phase pair target')
        require(isinstance(iterator,ast.Call) and isinstance(iterator.func,ast.Name) and iterator.func.id=='zip'
                and len(iterator.args)==len(clause.target.elts) and not iterator.keywords and all(isinstance(n,ast.Name) for n in iterator.args),
                'native phase paired coordinates')
        sequences=[env[n.id] for n in iterator.args]
        require(all(len(s)==3 for s in sequences),'actual native three-dimensional phase')
        names=[n.id for n in clause.target.elts];require(len(set(names))==len(names),'distinct native phase dummy names')
        return sum(scalar(generator.elt,{**env,**dict(zip(names,values))}) for values in zip(*sequences))
    return scalar(root.args[0],bindings)


def control_eligible(address,kind):
    """Eligibility only. Actual endpoint nonzero/finite checks occur under guard."""
    if address['status']!='FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED':return False
    if kind in ('omit-height-contact','reverse-translation-phase'):
        return tuple(address['responseGrade'])==(1,0) and address['component']=='NATIVE_HEIGHT'
    if kind=='omit-lower-normal-sign':
        return (address['face']=='minus' and address['slot']=='normal' and
                tuple(address['responseGrade'])==(0,0) and address['component']=='NATIVE_FLAT')
    raise ValueError('unknown addressed control')


def run_science(manifest,J,ns):
    raw={};copies={};used=set();D=ns['decode'];audit=EvidenceLog(J.out/'exact-evidence.jsonl',J.encode)
    for alias,v in manifest['savedInputs'].items():
        src=Path(v['path']);dst=J.out/'saved'/alias;dst.parent.mkdir(parents=True,exist_ok=True)
        require(sha(src)==v['sha256'],'saved source '+alias);shutil.copyfile(src,dst)
        require(sha(dst)==v['sha256'],'saved copy '+alias)
        copies[alias]={'source':str(src),'path':str(dst.relative_to(J.out)),'sha256':v['sha256'],'bytes':dst.stat().st_size}
        raw[alias]=json.loads(dst.read_text())
    save(J.out/'saved-copy-index.json',copies)
    def take(alias):used.add(alias);return raw[alias]
    def load(alias):return D(take(alias))
    full_receipts=take('full/artifact-index.json')
    for alias,v in manifest['savedInputs'].items():
        if alias.startswith('full/') and alias!='full/artifact-index.json':
            name=alias.split('/',1)[1];receipt=full_receipts[name]
            require(receipt['path']==name and receipt['sha256']==v['sha256'] and receipt['bytes']==v['bytes'],
                    'completed full-weak artifact receipt '+name)
    def zero(label,left,right):
        audit.append('zero-input',{'name':label,'left':left,'right':right})
        residual=sp.cancel(sp.together(left-right));record={'name':label,'left':left,'right':right,'cancelled':residual}
        audit.append('zero-return',record);require(residual==0,label);return record
    def constant(value,nonzero=False):
        audit.append('constant-input',{'value':value,'nonzeroRequired':nonzero})
        require(not value.free_symbols,'constant has no free symbols')
        real,imag=sp.expand_complex(value).as_real_imag();real=sp.cancel(real);imag=sp.cancel(imag)
        rec={'value':value,'real':real,'imaginary':imag};audit.append('constant-components',rec)
        require(real.is_Rational is True and imag.is_Rational is True and
                real.is_finite is True and imag.is_finite is True,'finite rational components')
        rec['reconstruction']=zero('constant-reconstruction',value,real+sp.I*imag)
        require(not nonzero or real!=0 or imag!=0,'exact constant nonzero')
        rec.update(finite=True,nonzero=bool(real!=0 or imag!=0));return rec
    def symbol(expr,name):return ns['one_symbol']([expr],name)
    def named_map(expr,bindings):
        found={s.name:s for s in expr.free_symbols}
        require(len(found)==len(expr.free_symbols),'unique symbol assumptions')
        require(set(found)<=set(bindings),'no omitted symbol mapping '+str(set(found)-set(bindings)))
        return {s:bindings[name] for name,s in found.items()}

    full=take('full/complete-weak-assembly.json');context=take('full/extended-binding-context.json')
    physical=json.loads(Path(manifest['physicalInput']).read_text())
    require(context['physical']==context['saved']['physicalInput']==physical,'actual local and pressure physical input')
    numeric=D(context['numeric']);eta,sigma=D(context['independentGrades']);eps=D(context['epsilon'])
    require(numeric['omega']==3 and numeric['W_0']==1 and numeric['L_W']==10 and
            numeric['rho_m']==sp.Rational(1,10) and numeric['Lambda_A_0']==sp.Rational(1,100)
            and numeric['tau_A']==sp.Rational(1,10),'same actual physical values')
    require('c_s0' not in numeric and full['status']=='SOURCE_JOINED_FULL_RETAINED_WEAK_OPERATOR'
            and full['nativeEpsilonOnce'] is True and full['independentGrades']==[list(g) for g in G],
            'accepted full weak scope')
    speed_census=take('full/native-raw-symbol-census.json')
    require(speed_census['speedAbsent'] is True and 'c_s0' not in speed_census['names'], 'raw local speed absence')
    # Ordinary whole kernels are inherited with their actual hypotheses, not recomputed.
    envelopes=take('weak/global-whole-kernel-envelopes.json');weak=take('weak/analytic-conclusion.json')
    require(weak['status']=='SOURCE_JOINED_GLOBAL_WEAK_PRESSURE_CERTIFICATE' and
            weak['wholeDirectOnce'] is True and weak['nativeIterationOnce'] is True and
            envelopes['pressureWeight']=='(1+|k|+|l|)^2' and envelopes['normalWeight']=='(1+|k|+|l|)^3',
            'inherited whole bounds and continuity')
    assembly=take('full/pressure-source-assembly-joins.json');cells=take('full/all-local-cells.json')
    require(len(cells)==400 and len({cell_key(c) for c in cells})==400,'complete saved local cell census')
    expected={(r,f,g,n) for r in ROWS for f in FIELDS for g in G for n in range(4)}
    require({cell_key(c) for c in cells}==expected,'complete local route grid')
    for c in cells:
        require(c['identities'][0]['left']==c['coefficient'] and
                c['coefficient']==c['polynomial']['original'],'actual local coefficient proof inputs')
        for z in c['identities']:inherit_zero(z)
    J.emit('inherited-local-context',{'physical':physical,'context':context,'fullAssembly':full,
        'localCells':len(cells),'oldEndpointCalculationsCalled':False,'localSourceDerivationsCalled':False})
    fields=take('inventory/fields.json');certs=take('weak/all-coefficient-certificates.json')
    require(set(fields)==set(certs) and len(fields)==34,'all pressure coefficient fields')
    endpoints={}
    for fid in fields:
        op=take('weak/field-'+fid+'-operands.json');poly=take('weak/field-'+fid+'-polynomial.json')
        der=take('weak/field-'+fid+'-derivative-class.json');proof=take('weak/field-'+fid+'-reconstruction-input.json')
        ret=take('weak/field-'+fid+'-reconstruction-return.json')
        joined=literal_field_joins(fid,fields,certs,op,poly,der,proof,ret)
        def endpoint_work():
            polynomial=D(certs[fid]['polynomial']);T=D(op['T'])
            require(polynomial.free_symbols<={T},'saved one-variable polynomial')
            rec={'inherited':joined,'newEndpointArguments':[{'variable':T,'value':-1},{'variable':T,'value':1}]}
            audit.append('field-endpoint-input',rec)
            actual_argument=polynomial.xreplace({T:sp.tanh(D(op['x'])/D(op['profileScale']['savedLength']))})
            rec['restoredProofArgument']={'actual':actual_argument,'saved':D(proof['right'])}
            audit.append('field-restored-proof-argument',rec['restoredProofArgument'])
            require(actual_argument==D(proof['right']),'actual completed polynomial proof arguments')
            rec['minus']=sp.cancel(polynomial.subs(T,-1));rec['plus']=sp.cancel(polynomial.subs(T,1))
            rec['minusCertificate']=constant(rec['minus']);rec['plusCertificate']=constant(rec['plus'])
            rec['analyticConvergence']='Saved bounded-derivative recurrence plus compact convergence and Schwartz tails; not a machine topology proof.'
            return rec
        endpoints[fid]=J.stage('new-field-endpoints-'+fid,joined,endpoint_work)

    # Common real-frequency physical depth remains a formal argument until later tasks.
    p=sp.Symbol('weak_end_p',real=True);q=sp.Symbol('weak_end_q');cs=sp.Symbol('weak_end_cs',positive=True)
    mu=numeric['rho_m']*numeric['omega'];alpha=numeric['Lambda_A_0']/(numeric['rho_m']**2*(1-sp.I*numeric['omega']*numeric['tau_A']))
    beta=sp.cancel(alpha*mu);constant(beta,True);zero('actual-beta-at-three',beta,(30+9*sp.I)/109)
    edge=(numeric['s11cdTangentialMomentum1'],numeric['s11cdTangentialMomentum2'])
    rad=9/cs**2-sum(t*t for t in edge)-p*p;F0=mu/(q+beta);Bh=-sp.I*mu*q/(q+beta)
    W=numeric['W_0'];L=numeric['L_W'];H={'minus':sp.S.Zero,'plus':W/2}
    profile=load('reference/restored-profile-arguments.json');global_profile=load('weak/global-profile-envelope.json')
    A0=profile['removableA0'];zero('saved-profile-A0-join',A0,global_profile['removableValue'])
    zero('profile-A0-normalization',2*sp.pi*A0,1)
    require(profile['noProfileTransformReplayed'] is True and global_profile['L']==10,'actual saved profile family')
    census=load('saved/reference/retained-response-census.json')
    require(census['directMultiplicity']==1 and census['middleIntegrationOfWholeDirect'] is False,'direct once')
    flat=census['flat'];flat_map=named_map(flat,{'reference_unrestricted_frequency':sp.Integer(3),'reference_qi':q})
    zero('saved-flat-to-new-closed-end',flat.xreplace(flat_map),F0)
    closure_return=take('reference/unrestricted-closure-00-return.json')
    inherit_zero(closure_return)
    J.emit('inherited-flat-defining-equation',{'operands':take('reference/unrestricted-closure-00-input.json'),
        'returned':closure_return,'censusFlat':flat,'newDiagonalMap':list(flat_map.items()),'oldFunctionCalled':False})
    J.emit('new-closed-depth-domain',{'p':p,'q':q,'cs':cs,'radicand':rad,'range':[1,2],
        'prescription':'q=sqrt(radicand) for positive radicand, i*sqrt(-radicand) for negative, 0 at grazing; outgoing boundary value',
        'beta':beta,'mu':mu,'alpha':alpha,'savedGlobalEnvelope':envelopes,
        'nonzeroDenominatorArgument':'q in closed first quadrant and Re beta=30/109>0; |q+beta|>=30/109',
        'grazingExtension':'Only closed q+beta expressions evaluated at q=0; no raw 0/0',
        'sourceFrequencyContinued':False,'rootsSolved':False})

    # New translation/Dirichlet algebra; the theorems used below are assessed arguments.
    a,k,l,x,y,Q=sp.symbols('weak_end_a weak_end_k weak_end_l weak_end_x weak_end_y weak_end_Q',real=True)
    pA=profile['heightNumerator'];gA=global_profile['A']
    zero('actual-profile-numerator-argument',pA.xreplace({symbol(pA,'reference_t'):Q}),
         gA.xreplace({symbol(gA,'weak_t'):Q}))
    duality=take('weak/new-weak-duality-and-order.json');native=Path(manifest['c2Source']).read_text()
    Fourier=duality['inheritedConvention'];source_statement_join(native,'phase1 = '+Fourier['sourcePhase'])
    source_statement_join(native,'phase = '+Fourier['profilePhase'])
    require(duality['sourceDerivativeBeforeMultiplication'] is True and duality['normalVariable']=='response output l'
            and Fourier['profileForwardPower']==-3 and Fourier['invariantEdgeCoordinates']==2,'inherited Fourier/jet order')
    x2,x3,y2,y3=sp.symbols('weak_end_x2 weak_end_x3 weak_end_y2 weak_end_y3',real=True)
    coordinate_bindings={'kout':(l,*edge),'kin':(k,*edge),'ko':(l,*edge),'ki':(k,*edge),
                         'X':(x,x2,x3),'Y':(y,y2,y3)}
    native_source_exponent=native_phase_exponent(Fourier['sourcePhase'],coordinate_bindings,sp.I)
    native_profile_exponent=native_phase_exponent(Fourier['profilePhase'],coordinate_bindings,sp.I)
    shifted_source=native_source_exponent.subs({x:x+a,y:y+a},simultaneous=True)
    phase_increment=sp.expand(shifted_source-native_source_exponent)
    phase_coefficient=sp.cancel(phase_increment/(sp.I*(l-k)*a))
    profile_coefficient=sp.cancel(native_profile_exponent/(-sp.I*(l-k)*y))
    audit.append('native-phase-operands',{'native':Fourier,'bindings':coordinate_bindings,'sourceExponent':native_source_exponent,
        'profileExponent':native_profile_exponent,'shiftedSourceExponent':shifted_source,'increment':phase_increment,
        'sourceSignCoefficient':phase_coefficient,'profileSignCoefficient':profile_coefficient})
    phase=zero('native-translation-phase',phase_increment,sp.I*(l-k)*a)
    zero('native-profile-forward-phase',native_profile_exponent,-sp.I*(l-k)*y)
    zero('native-source-profile-sign-join',phase_coefficient,profile_coefficient)
    # The new translated Q phase is derived from the bound native source sign.
    translated_exponent=sp.I*phase_coefficient*Q*a
    flipped_profile=native_phase_exponent(Fourier['profilePhase'],coordinate_bindings,-sp.I)
    flipped_profile_coefficient=sp.cancel(flipped_profile/(-sp.I*(l-k)*y))
    phase_sign_movement=sp.cancel(flipped_profile_coefficient-phase_coefficient)
    phase_sign_certificate=constant(phase_sign_movement,True)
    J.emit('new-native-phase-argument-join',{'sourceExpression':Fourier['sourcePhase'],'profileExpression':Fourier['profilePhase'],
        'bindings':coordinate_bindings,'sourceExponent':native_source_exponent,'profileExponent':native_profile_exponent,
        'shiftedSource':shifted_source,'phaseIncrement':phase_increment,'newTransferExponent':translated_exponent,
        'profileSignMutation':{'actualExpressionExponent':native_profile_exponent,'flippedExponent':flipped_profile,
            'movement':phase_sign_movement,'certificate':phase_sign_certificate},
        'nativeFunctionsCalled':False,'commonEdgesCancelInProfileAndTranslation':True})
    f,fzero,A,chi=sp.symbols('weak_end_f weak_end_fzero weak_end_A weak_end_chi')
    E=sp.exp(translated_exponent);subtracted=E*A*(f-fzero*chi)/Q;diagonal=E*A*fzero*chi/Q
    decomposition=zero('translated-PV-complete-decomposition',subtracted+diagonal,E*A*f/Q)
    limits={side:sp.cancel(W/4+W/(2*sp.I)*(sp.I*sp.pi*A0*s*phase_coefficient)) for side,s in [('minus',-1),('plus',1)]}
    for side in H:zero('native-half-height-'+side,limits[side],H[side])
    reversed_exponent=-translated_exponent
    reversed_sign=sp.cancel(reversed_exponent/(sp.I*Q*a))
    reversed_scalar={side:sp.I*sp.pi*A0*s*reversed_sign for side,s in [('minus',-1),('plus',1)]}
    reversed_limits={side:sp.cancel(W/4+W/(2*sp.I)*value) for side,value in reversed_scalar.items()}
    zero('reversed-phase-left-height-exchange',reversed_limits['minus'],limits['plus'])
    zero('reversed-phase-right-height-exchange',reversed_limits['plus'],limits['minus'])
    J.emit('new-reversed-native-phase-limits',{'originalExponent':translated_exponent,'reversedExponent':reversed_exponent,
        'reversedPVScalarLimits':reversed_scalar,'contact':W/4,'originalHeightLimits':limits,'reversedHeightLimits':reversed_limits,
        'samePVSubtraction':'Reverse E in BOTH subtracted and oscillatory diagonal terms; ordinary L1 terms still vanish weakly',
        'argument':'Reviewed Dirichlet limit is odd under reversal of the actual native phase; no integral evaluated'})
    for face,sign in [('plus',1),('minus',-1)]:
        pv=load('weak/'+face+'-global-height-PV-certificate.json')
        require(pv['normalCoefficientBoundedDirectly'] is True,'normal Holder certificate, no C1(qY) assumption')
        bindings={'weak_Omega':sp.Integer(3),'weak_qi':q,'weak_qo':q}
        mapped=pv['pressureCoefficient'].xreplace(named_map(pv['pressureCoefficient'],bindings))
        mappedn=pv['normalCoefficient'].xreplace(named_map(pv['normalCoefficient'],bindings))
        zero('saved-height-Holder-diagonal-'+face,mapped,Bh)
        zero('saved-normal-Holder-diagonal-'+face,mappedn,sp.I*sign*q*Bh)
    J.emit('new-translated-response-argument',{'phaseIdentity':phase,'phase':E,'dualityInherited':duality,
        'PVDecomposition':decomposition,'contact':W*fzero/4,'subtractedTerm':subtracted,'oscillatoryDiagonalTerm':diagonal,
        'bothPVTermsPrefactor':W/(2*sp.I),
        'profileA0':A0,'DirichletScalarLimits':{'minus':-sp.I*sp.pi*A0,'plus':sp.I*sp.pi*A0},'heightLimits':limits,
        'coefficientTopology':'P(tanh((y+a)/L))*D^n u -> P(+-1)*D^n u in S, using inherited all-derivative bounds, compact convergence and Schwartz tails.',
        'uniformTranslatedBound':'Ordinary kernels: absolute weighted L1 independent of a. Height subtraction: inherited Holder 1/2 and tails. Diagonal PV: split A0 plus integrable (A-A0)/Q; sine integral bound 4, without differentiating exp(iQa).',
        'ordinaryLimits':'Slope, C*H, Jwhole and Dwhole tend weakly to zero by Riemann-Lebesgue on weighted L1(R^2), not pointwise zero.',
        'heightLimitsArgument':'Subtracted L1 term and (A-A0)*chi/Q term vanish by Riemann-Lebesgue. Dirichlet gives i*pi*A0*sign(a); add surviving contact.',
        'compactCsUniformity':'Inherited weighted L1 continuity and common polynomial bounds yield compact L1 image; finite net makes Riemann-Lebesgue uniform. Inherited Holder and beta bounds control the height family.',
        'parameterScope':'Real3 source/local; delta auxiliary in inherited response proof only. New symbols at delta=0.',
        'analyticTheoremsMachineProved':False,'integralsEvaluated':False,'decayRateOrFiniteBoxError':False})

    # Select the independent height grade; slope remains live and is not pointwise zero.
    first=load('reference/native-first-shape-input.json')['right']
    first_parts=first_shape_grade_parts(first,eta,sigma,sp.expand)
    audit.append('native-first-shape-grade-operands',{'savedOriginal':first,'eta':eta,'sigma':sigma,'parts':first_parts})
    first_reconstruction=zero('native-first-shape-independent-grade-reconstruction',
                              first,first_parts['height']+first_parts['slope'])
    require(not (first_parts['heightCoefficient'].free_symbols | first_parts['slopeCoefficient'].free_symbols)
            & {eta,sigma},'independent first-shape coefficient grades')
    first_bindings={'reference_unrestricted_frequency':sp.Integer(3),'reference_qi':q,'reference_qo':q,
                    'reference_k':p,'reference_l':p,eta.name:eta,sigma.name:sigma}
    height_diag=first_parts['height'].xreplace(named_map(first_parts['height'],first_bindings))
    slope_diag=first_parts['slope'].xreplace(named_map(first_parts['slope'],first_bindings))
    J.emit('new-native-first-shape-grade-join',{'savedOriginal':first,'parts':first_parts,
        'reconstruction':first_reconstruction,'bindings':first_bindings,'heightOnCommonDepthSupport':height_diag,
        'slopeOnCommonDepthSupport':slope_diag,'slopePointwiseZeroAsserted':False,
        'slopeWeakLimit':'Inherited ordinary-kernel translated Riemann-Lebesgue argument, not j(0)=0 or sigma=0',
        'domain':'First-height comparison on common nongrazing depth; eta and sigma stay independent'})
    new_native={}
    for face,sign in [('plus',1),('minus',-1)]:
        tr=load('reference/'+face+'-new-native-trace.json');slot=load('reference/'+face+'-final-native-slot-routing.json')
        hc=load('reference/'+face+'-native-height-constant.json')
        for statement in tr['nativeAssignmentSources'].values():source_statement_join(native,statement)
        source_statement_join(native,slot['source'])
        height=tr['restoredTrace']['height'];normal=tr['restoredTrace']['normalJet']
        require(height==hc['savedHeight'] and tr['restoredTrace']['valueCoefficient']==1 and slot['valueCoefficient']==1,
                'actual native height/value inputs')
        oldmap=dict(tr['actualHeightMap']);mappedheight=height.xreplace(oldmap)
        zero('native-height-map-'+face,mappedheight,slot['savedHeight'])
        original_q=symbol(normal,'q_o');normal_q=normal.xreplace({original_q:q})
        zero('native-outward-jet-'+face,normal_q,sp.I*sign*q)
        hsymbol=hc['profileSymbol'];record={}
        for side,endpoint in [('minus',0),('plus',1)]:
            lab=height.subs(hsymbol,endpoint);product=sp.cancel(lab*normal_q)
            zero('constant-height-trace-factor-'+face+'-'+side,product,sp.I*q*eta*H[side])
            dtn=height_diag.replace(lambda z:getattr(z,'is_Function',False) and z.func.__name__=='reference_height_hat',lambda z:H[side])
            zero('constant-height-physical-first-height-grade-'+face+'-'+side,dtn,0)
            trace=1+product;inverse_retained=1-product
            excluded=zero('constant-height-retained-inverse-'+face+'-'+side,trace*inverse_retained,1-product**2)
            candidate=sp.expand(inverse_retained*F0)
            zero('native-end-response-'+face+'-'+side,candidate,F0+eta*H[side]*Bh)
            # Actual affine native reference-location equation; insert flat jet and physical flat closure.
            sol=slot['equationSolution'];named={face+'_affine_height':lab,face+'_jet_slot':normal_q*F0,face+'_physical_target':F0}
            affine=sol.xreplace(named_map(sol,named))
            zero('native-affine-end-slot-'+face+'-'+side,affine,F0+eta*H[side]*Bh)
            record[side]={'labHeight':lab,'normal':normal_q,'trace':trace,'inverseRetained':inverse_retained,
                'excludedEta2':-product**2,'inverseIdentity':excluded,'reference':candidate,'normalResponse':normal_q*candidate,
                'affineReference':affine,'closedGrazingReference':sp.cancel(candidate.subs(q,0)),
                'closedGrazingNormal':sp.cancel((normal_q*candidate).subs(q,0))}
        new_native[face]=record
        J.emit('new-native-constant-height-'+face,{'restoredTrace':tr,'nativeSlot':slot,'heightContext':hc,
            'actualNewConstantSpecializations':record,'domain':'raw first-height grade specialized for q!=0; closed response extended continuously to q=0',
            'noOldTraceOrClosureFunctionCalled':True})

    # Every address is joined to its published source, factor, wave and whole-tag record.
    coverage={r['addressId']:r for r in take('weak/weak-address-coverage.json')}
    factors=take('full/new-pressure-factor-arguments.json');waves=take('full/new-pressure-wave-arguments.json')
    require(len(coverage)==13260 and len(factors)==17 and len(waves)==39,'accepted route census')
    groups={(side,row,field,g):{'local':[],'pressure':[],'localAncestry':[],'addressIds':[]} for side in H for row in ROWS for field in FIELDS for g in G}
    local_records=[]
    for c in cells:
        route=cell_key(c);rec={'savedCell':c,'newTerms':{}}
        for side,endkey in [('minus','leftEndpoint'),('plus','rightEndpoint')]:
            val=D(c[endkey]);require(not val.free_symbols,'saved local endpoint constant')
            term=val*(sp.I*p)**c['xOrder'];rec['newTerms'][side]=term
            group=groups[(side,c['row'],c['field'],tuple(c['grade']))]
            group['local'].append(term);group['localAncestry'].append({'row':c['row'],'field':c['field'],'xOrder':c['xOrder'],'grade':c['grade']})
        local_records.append(rec)
    J.emit('restored-local-endpoints-new-symbol-terms',local_records)
    address_records=[];join_records=[];seen=set();response_cache={};possible_controls={kind:[] for kind in
        ('omit-height-contact','reverse-translation-phase','omit-lower-normal-sign')}
    assembly_by_row={r['row']:r for r in assembly}
    for row in ROWS:
        addresses=take('inventory/'+row+'-ordered-addresses.json')
        require(assembly_by_row[row]['addressIds']==[a['addressId'] for a in addresses],'actual full weak row address ancestry')
        row_results=[]
        try:
            for ar in addresses:
                idx=ar['addressId'];require(idx not in seen,'unique address');seen.add(idx)
                joined=address_metadata_join(ar,coverage[idx],fields,factors,waves);join_records.append(joined)
                rg=tuple(ar['responseGrade']);face=ar['face'];sign=1 if face=='plus' else -1
                source=endpoints[ar['sourceTransform']['coefficientId']];consumer=endpoints[ar['consumerTransform']['coefficientId']]
                wave=D(ar['waveMultiplier']);wp=D(waves[ar['jet']['name']]['sourceMomentum']);wave=wave.xreplace({wp:p})
                require(wave.free_symbols<={p},'source derivative at diagonal original p')
                cache_key=(ar['fullFactorProof']['proof'],face,ar['slot'])
                if cache_key not in response_cache:
                    expr=D(ar['responseCoefficient']);normal=D(ar['normalMultiplier']);original=expr
                    maps={s:p for s in expr.free_symbols|normal.free_symbols if s.name in ('composition_p','composition_k','composition_l','composition_r')}
                    if rg in ((0,0),(1,0)):
                        mapped=expr.xreplace(maps);normal_end=normal.xreplace(maps)
                        mapped=mapped.replace(lambda z:getattr(z,'is_Function',False) and z.func.__name__=='common_outgoing_q',lambda z:q)
                        normal_end=normal_end.replace(lambda z:getattr(z,'is_Function',False) and z.func.__name__=='common_outgoing_q',lambda z:q)
                        zero('addressed-native-normal-'+str(idx),normal_end,1 if ar['slot']=='pressure' else sp.I*sign*q)
                        if rg==(0,0):
                            zero('addressed-flat-end-'+str(idx),mapped,F0);response={side:F0 for side in H}
                        else:
                            hs=list(mapped.atoms(sp.Function));require(len(hs)==1 and hs[0].func.__name__=='reference_height_hat' and hs[0].args==(0,), 'actual height transform after diagonal support')
                            response={side:mapped.xreplace({hs[0]:H[side]}) for side in H}
                            for side in H:zero('addressed-height-end-'+str(idx)+'-'+side,response[side],H[side]*Bh)
                    else:
                        normal_end=None;response={side:sp.S.Zero for side in H}
                    response_cache[cache_key]={'original':original,'responseGrade':rg,'normal':normal_end,'response':response,
                        'reason':'flat delta or contact plus Dirichlet limit' if rg in ((0,0),(1,0)) else 'ordinary-kernel translated weak Riemann-Lebesgue limit; not pointwise zero'}
                rr=response_cache[cache_key];record={'addressId':idx,'inheritedAddress':ar,'join':joined,'waveAtP':wave,'endContributions':{},'responseLimit':rr}
                audit.append('address-end-input',record)
                for side in H:
                    term=sp.S.Zero if rg in ((0,1),(1,1)) else source[side]*consumer[side]*wave*rr['normal']*rr['response'][side]
                    term=sp.cancel(term)
                    require(term.free_symbols<={p,q},'only retained diagonal depth and momentum')
                    record['endContributions'][side]=term
                    group=groups[(side,row,ar['jet']['channel'],tuple(ar['targetGrade']))]
                    group['pressure'].append(term);group['addressIds'].append(idx)
                    if side=='plus' and term!=0:
                        for kind in possible_controls:
                            if control_eligible(ar,kind):
                                possible_controls[kind].append((ar,source,consumer,wave,rr,record['endContributions']))
                audit.append('address-end-return',record);row_results.append(record)
        except BaseException:
            J.emit(row+'-partial-end-address-returns',row_results);raise
        J.emit(row+'-new-end-addresses',row_results);address_records.extend(row_results)
    require(seen==set(coverage),'all pressure addresses retained including zeros')
    require(len({tuple(tuple(a[k]) for k in ('consumerGrade','responseGrade','sourceGrade')) for r in ROWS for a in raw['inventory/'+r+'-ordered-addresses.json']})==16,
            'all sixteen ordered grade triples')
    symbols=[]
    for (side,row,field,grade),v in groups.items():
        rec={'side':side,'row':row,'field':field,'grade':grade,**v};audit.append('end-symbol-input',rec)
        rec['localSum']=sp.cancel(sum(v['local'],sp.S.Zero));rec['pressureSum']=sp.cancel(sum(v['pressure'],sp.S.Zero))
        rec['symbol']=sp.cancel(rec['localSum']+rec['pressureSum'])
        rec['sumIdentity']=zero('new-end-symbol-sum',rec['symbol'],sum(v['local']+v['pressure'],sp.S.Zero))
        rec['closedGrazingValue']=sp.cancel(rec['symbol'].subs(q,0));require(rec['closedGrazingValue'].free_symbols<={p},'closed finite grazing symbol')
        require(not rec['closedGrazingValue'].has(sp.zoo,sp.nan,sp.oo,-sp.oo),'finite closed grazing expression')
        rec['epsilonPower']=0 if rec['symbol']==0 else 1
        rec['epsilonNormalization']='Symbol is the coefficient of the one native epsilon_shape; exact zero has no power'
        rec['weakPairingFactor']=2*sp.pi
        audit.append('end-symbol-return',rec);symbols.append(rec)
    J.emit('new-complete-retained-weak-end-symbols',symbols)

    # New responsive mutations applied through actual surviving address and its full symbol cell.
    controls=[];point={p:sp.Integer(1),q:sp.Integer(2)};cs2=sp.Rational(180,101)
    zero('control-outgoing-physical-dispersion',rad.subs({p:1,cs:sp.sqrt(cs2)}),4)
    require(cs2>1 and cs2<4,'control speed in scoped interval')
    for kind in ('omit-height-contact','reverse-translation-phase','omit-lower-normal-sign'):
        chosen=None
        for ar,source,consumer,wave,rr,end_terms in possible_controls[kind]:
            per_side={};right_nonzero=False
            for side in (('minus','plus') if kind=='reverse-translation-phase' else ('plus',)):
                base=source[side]*consumer[side]*wave;term=end_terms[side]
                if kind=='omit-height-contact':mutated=term-base*rr['normal']*W*Bh/4
                elif kind=='reverse-translation-phase':mutated=base*rr['normal']*reversed_limits[side]*Bh
                else:mutated=-term
                movement=sp.cancel(mutated-term);actual=sp.cancel(movement.subs(point))
                audit.append('control-candidate-input',{'kind':kind,'side':side,'address':ar,'base':base,'response':rr,
                    'original':term,'mutated':mutated,'movement':movement,'point':list(point.items()),'actual':actual,
                    'reversedPhaseLimits':reversed_limits if kind=='reverse-translation-phase' else None})
                cert=constant(actual,False);right_nonzero=right_nonzero or (side=='plus' and cert['nonzero'])
                key=(side,ar['row'],ar['jet']['channel'],tuple(ar['targetGrade']))
                cell=next(s for s in symbols if (s['side'],s['row'],s['field'],tuple(s['grade']))==key)
                mutated_cell=sp.cancel(cell['symbol']-term+mutated)
                identity=zero('control-full-row-movement-'+kind+'-'+side,mutated_cell-cell['symbol'],movement)
                per_side[side]={'originalTerm':term,'mutatedTerm':mutated,'actualOriginalCell':cell['symbol'],
                    'actualMutatedCell':mutated_cell,'cellMovement':movement,'movementCertificate':cert,'rowIdentity':identity}
            if not right_nonzero:continue
            if kind=='reverse-translation-phase' and not per_side['minus']['movementCertificate']['nonzero']:continue
            chosen={'kind':kind,'actualAddressId':ar['addressId'],'actualAddress':ar,'endCells':per_side,
                'point':{'p':1,'q':2,'csSquared':cs2,'outgoing':True},
                'interpretation':'formal coefficient sensitivity, not computed field or power',
                'reversePhaseEvidence':'new-reversed-native-phase-limits.json' if kind=='reverse-translation-phase' else None}
            break
        J.emit('new-control-'+kind,{'selected':chosen,'availableSurvivingAddresses':len(possible_controls[kind])})
        require(chosen is not None,'responsive applicable actual-row control '+kind);controls.append(chosen)
    J.emit('inherited-address-argument-joins',join_records)
    J.emit('consumed-inputs',{'usedAliases':sorted(used),'allCopiedAliases':sorted(copies),'oldFunctionsCalled':False})
    conclusion={'status':'SOURCE_JOINED_TRANSLATED_RETAINED_WEAK_ENDS','localCells':400,'pressureFields':34,'addresses':len(seen),
        'endSymbolCells':len(symbols),'orderedGradeTriples':16,'controls':len(controls),'testSpace':'S(R)^5 x S(R)^5, complex bilinear',
        'heights':H,'pairing':'B_end,g(v,u)=2pi integral hat(v)(-p)^T E_end,g(p) hat(u)(p) dp',
        'scope':manifest['scope'],'analyticArgument':'Assessed Schwartz convergence, uniform translated response bounds, Riemann-Lebesgue and PV/Dirichlet limits; finite checks are not machine functional analysis.',
        'localPressureCoefficientFunctionsReplayed':False,'integralEvaluated':False,'finiteMatrixBuilt':False,'scatteringOrCurrent':False,
        'loss':False,'oldModeAcceptance':False,'decayRateOrFiniteBoxError':False,'primitiveCalibration':False,'drain':False}
    J.emit('translated-weak-end-conclusion',conclusion)
    return {'executionStatus':'COMPLETE_PENDING_INSPECTION','allChecksPassed':True,**conclusion}


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
        result=J.stage('translated-weak-ends',{'manifestSha256':gate['manifestSha256'],'reviewSha256':gate['buildReviewRecordSha256']},lambda:run_science(manifest,J,ns));code=0
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
