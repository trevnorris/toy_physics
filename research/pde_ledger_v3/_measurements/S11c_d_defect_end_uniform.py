#!/usr/bin/env python3
"""Bounded source-joined end/uniform equation comparison; no producer or roots."""
import argparse,ast,base64,builtins,hashlib,importlib,io,json,math,os,pickle,re,resource,shutil,sqlite3,sys,time,traceback
from collections import OrderedDict
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
G=((0,0),(1,0),(0,1),(1,1))
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO={'text':'0','srepr':'Integer(0)'}


def require(value,message):
    if value is not True:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()
    with path.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())


def definitions(text,names):
    nodes=[n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in nodes}==set(names),'inert definition census')
    return ast.Module(body=nodes,type_ignores=[])


def source_method(text,cls,name):
    c=next(n for n in ast.parse(text).body if isinstance(n,ast.ClassDef) and n.name==cls)
    return ast.unparse(next(n for n in c.body if isinstance(n,ast.FunctionDef) and n.name==name))


def exact_structure(a,b):
    if type(a) is not type(b):return False
    if isinstance(a,dict):return a.keys()==b.keys() and all(exact_structure(a[k],b[k]) for k in a)
    if isinstance(a,(list,tuple)):return len(a)==len(b) and all(exact_structure(x,y) for x,y in zip(a,b))
    if isinstance(a,np.ndarray):return a.shape==b.shape and a.dtype==b.dtype and a.tobytes()==b.tobytes()
    r=a==b
    return r is True or r is sp.S.true


def convolution(a,b,rectangle=True):
    out={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            if not rectangle or (i+k<=1 and j+l<=1):out[i+k,j+l]=out.get((i+k,j+l),0)+x*y
    return out


def rotated_name(name,scalar_names):
    """Proper cyclic basis map: weak axes 1,2,3 are uniform axes 3,1,2."""
    cycle={'1':'3','2':'1','3':'2'}
    m=re.fullmatch(r'(u_([123])|theta|e_W|w1_profile|m1_profile)((?:_t|t)*)(.*)',name)
    if m:
        head,component,time_suffix,derivatives=m.groups()
        require(re.fullmatch(r'(?:_?d[123])*',derivatives) is not None,'native jet suffix')
        if component:head='u_'+cycle[component]
        directions=re.findall(r'd([123])',derivatives)
        suffix='' if not directions else '_'+''.join('d'+i for i in sorted(cycle[d] for d in directions))
        return head+time_suffix+suffix
    m=re.fullmatch(r'grad_theta_([123])',name)
    if m:return 'grad_theta_'+cycle[m[1]]
    require(name in scalar_names,'unclassified native tensor/scalar '+name)
    return name


class Evidence:
    def __init__(self,out):self.out=out;self.active=None;self.count=0;self.previous=None;self.operations=[]
    def encode(self,v):
        if isinstance(v,(sp.Basic,sp.MatrixBase,FunctionClass)):return {'text':str(v),'srepr':sp.srepr(v)}
        if isinstance(v,np.ndarray):return {'numpyDtype':v.dtype.str,'shape':list(v.shape),'bytesBase64':base64.b64encode(v.tobytes()).decode()}
        if isinstance(v,dict):
            if all(isinstance(k,str) for k in v):return {k:self.encode(x) for k,x in v.items()}
            return {'mappingPairs':[[self.encode(k),self.encode(x)] for k,x in v.items()]}
        if isinstance(v,(list,tuple)):return [self.encode(x) for x in v]
        if isinstance(v,(set,frozenset)):return {'setValues':[self.encode(x) for x in sorted(v,key=repr)]}
        if isinstance(v,(complex,np.complexfloating)):return {'pythonComplex':[float(v.real),float(v.imag)]}
        if isinstance(v,(np.integer,np.bool_,np.floating)):return v.item()
        require(v is None or type(v) in (str,bool,int,float),'supported evidence type '+str(type(v)))
        return v
    def emit(self,name,value):
        p=self.out/(name+'.json');save(p,self.encode(value))
        ref={'path':str(p.relative_to(self.out)),'sha256':sha(p),'bytes':p.stat().st_size}
        body={'sequence':self.count,'previousSha256':self.previous,'artifact':ref}
        h=hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest()
        with (self.out/'evidence-chain.jsonl').open('a') as f:
            f.write(json.dumps({'payload':body,'sha256':h})+'\n');f.flush();os.fsync(f.fileno())
        self.count+=1;self.previous=h;return ref
    def stage(self,name,arguments,callback):
        previous=self.active;self.active=name;inp=self.emit(name+'-input',arguments)
        value=callback();out=self.emit(name+'-return',value)
        record={'name':name,'status':'COMPLETE','input':inp,'result':out}
        with (self.out/'operation-index.jsonl').open('a') as f:f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
        self.operations.append(record);self.active=previous;return value
    def join(self,name,a,b):
        self.emit(name,{'actual':a,'expected':b,'join':'exact typed structure','passed':exact_structure(a,b)})
        require(exact_structure(a,b),name)
    def zero(self,name,a,b):
        previous=self.active;self.active=name;self.emit(name+'-input',{'left':a,'right':b})
        raw=a-b;v=sp.cancel(sp.together(raw));self.emit(name+'-return',{'raw':raw,'cancelled':v})
        require(v==0,name);self.active=previous;return v


def verify_gate(path,manifest_path,m):
    gate=json.loads(Path(path).read_text())
    require(gate['status']=='READY_FOR_ONE_END_UNIFORM_COMPARISON','actual readiness gate')
    require(gate['workerSha256']==sha(__file__) and gate['manifestSha256']==sha(manifest_path),'worker/manifest')
    require(gate['sourcePins']==m['sourcePins'],'source pins')
    for p,h in m['sourcePins'].items():require(sha(p)==h,'pin '+p)
    require(gate['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and gate['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'actual helper routes')
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(gate[key])==gate[key+'Sha256'],'gate '+key)
    require(gate['launcher']==m['launcher'] and gate['buildReviewRecord']==m['reviewRecordWillBe'] and gate['authority']==m['executionAuthority'],'actual manifest routes')
    review=json.loads(Path(gate['buildReviewRecord']).read_text())
    require(review['independentBuildClearance'] is True and gate['independentBuildClearance'] is True and review['allChecksPassed'] is True,'independent concrete build assessment')
    for name in ('claude','grok'):require(review['reports'][name]['literalVerdict']=='CLEAR FOR THIS BOUNDED END-UNIFORM BUILD','literal build verdict')
    for key in ('workerSha256','manifestSha256','launcherSha256','sharedGuardSha256','supervisorSha256'):require(review[key]==gate[key],'build identity '+key)
    method=json.loads(Path(m['methodRecord']).read_text())
    require(method['jointIndependentMethodClearance'] is True and method['methodSha256']==sha(m['methodPath'])==review['methodSha256'],'actual cleared method')
    authority=json.loads(Path(gate['authority']).read_text())
    require(m['sourcePins'][gate['authority']]==gate['authoritySha256'],'pinned exact authority')
    require(authority['boundedInstrumentAuthorized'] is True and authority['scienceExecutionsAuthorized']==1 and authority['automaticScientificRetry'] is False and authority['noDeadline'] is True and authority['scope']==m['scope']==gate['scope'],'bounded standing authority')
    require(gate['durationLimits'] is None and gate['scientificRunsAuthorized']==1,'no deadline or retry')
    return gate


def polynomial(expr,x,y):
    if not expr.has(x,y):return {(0,0):expr}
    if expr==x:return {(1,0):sp.S.One}
    if expr==y:return {(0,1):sp.S.One}
    if expr.is_Add:
        out={}
        for term in expr.args:
            for key,value in polynomial(term,x,y).items():out[key]=out.get(key,0)+value
        return out
    if expr.is_Mul:
        out={(0,0):sp.S.One}
        for term in expr.args:out=convolution(out,polynomial(term,x,y),False)
        return out
    if expr.is_Pow and expr.exp.is_Integer is True and 0<=expr.exp<=128:
        base=polynomial(expr.base,x,y);out={(0,0):sp.S.One}
        for _ in range(int(expr.exp)):out=convolution(out,base,False)
        return out
    raise ValueError('unsupported polynomial dependence; no EX domain inference')


def coefficients(value,eta,sigma,J,name):
    J.emit(name+'-input',{'value':value,'eta':eta,'sigma':sigma})
    numerator,denominator=sp.fraction(sp.together(value))
    n=polynomial(numerator,eta,sigma);d=polynomial(denominator,eta,sigma);d0=sp.cancel(d.get((0,0),0))
    J.emit(name+'-domain',{'numerator':numerator,'denominator':denominator,'originDenominator':d0,'regularityCondition':sp.Ne(d0,0,evaluate=False)})
    if d0==0:return {'status':'REGULARITY_UNAVAILABLE','originDenominator':d0}
    c={}
    for i,j in G:
        accum=sum(v*c[i-a,j-b] for (a,b),v in d.items() if (a,b)!=(0,0) and a<=i and b<=j and (i-a,j-b) in c)
        c[i,j]=sp.cancel((n.get((i,j),0)-accum)/d0)
    retained=sum(eta**i*sigma**j*c[i,j] for i,j in G);H=sp.cancel(value-retained)
    J.zero(name+'-reconstruction',value,retained+H)
    # Exact coefficient certificate for denominator*(value-retained), not a claim H=0.
    hc=convolution(d,c,False)
    for g in G:J.zero(name+'-excluded-grade-%d%d'%g,n.get(g,0)-hc.get(g,0),sp.S.Zero)
    out={'status':'EXTRACTED_ON_DECLARED_REGULAR_DOMAIN','coefficients':c,'retained':retained,'remainder':H,'originDenominator':d0}
    J.emit(name+'-return',out);return out


def finite_constant(expr):
    value=sp.simplify(expr)
    if value.free_symbols:return {'status':'UNRESOLVED_FREE_SYMBOLS','value':value}
    real,imag=sp.expand_complex(value).as_real_imag();real=sp.simplify(real);imag=sp.simplify(imag)
    reconstruction=sp.simplify(value-real-sp.I*imag)
    finite=real.is_finite is True and imag.is_finite is True and reconstruction==0
    nonzero=finite and (real.is_zero is False or imag.is_zero is False)
    return {'status':'NONZERO_CERTIFIED' if nonzero else ('ZERO' if finite and real==imag==0 else 'UNRESOLVED'),
            'value':value,'real':real,'imaginary':imag,'finite':finite,'reconstruction':reconstruction}


def wave_test(value,p,q,cs,J,name):
    J.emit(name+'-input',{'value':value,'wave':q**2-(9/cs**2-p**2-sp.Rational(1,20))})
    n,d=sp.fraction(sp.together(value));terms=polynomial(n,q,sp.Dummy('unused_q_coefficient'))
    require(all(b==0 for a,b in terms),'only physical depth reduced')
    radical=9/cs**2-p**2-sp.Rational(1,20);a=sp.S.Zero;b=sp.S.Zero
    for (order,_),v in terms.items():
        if order%2:b+=v*radical**(order//2)
        else:a+=v*radical**(order//2)
    a=sp.cancel(a);b=sp.cancel(b);remainder=a+b*q
    # Synthetic division builds an explicit quotient; no full-domain Poly.
    coeff={i:sp.cancel(v) for (i,_),v in terms.items()};quo={}
    for i in range(max(coeff,default=0),1,-1):
        leading=coeff.get(i,0);quo[i-2]=leading;coeff[i]=0;coeff[i-2]=coeff.get(i-2,0)+leading*radical
    quotient=sum(v*q**i for i,v in quo.items());J.zero(name+'-division',n,quotient*(q*q-radical)+remainder)
    domain=sp.Ne(d,0,evaluate=False);point={p:sp.S.One,q:sp.Integer(2),cs:sp.sqrt(sp.Rational(180,101))}
    witness={'point':point,'denominator':finite_constant(d.subs(point)),'value':finite_constant(value.subs(point))}
    zero=a==0 and b==0
    state='ZERO' if zero else ('NONZERO_CERTIFIED' if witness['denominator']['status']=='NONZERO_CERTIFIED' and witness['value']['status']=='NONZERO_CERTIFIED' else 'UNRESOLVED_NONZERO_SYMBOLIC_REMAINDER')
    out={'status':state,'numerator':n,'denominator':d,'domain':domain,'quotient':quotient,'remainder':remainder,
         'coefficient0':a,'coefficient1':b,'positiveSheetWitness':witness,'symbolicCoverage':'on wave, outgoing sheet, original denominator nonzero; no sample-to-interval inference'}
    J.emit(name+'-return',out);return out


def inspect_matrix(matrix,p,q,cs,J,name):
    return [[wave_test(matrix[i,j],p,q,cs,J,name+'-%d-%d'%(i,j)) for j in range(matrix.cols)] for i in range(matrix.rows)]


def run_science(m,J,ns,codec):
    D=ns['decode'];copied={};raw={}
    for alias,pin in m['savedInputs'].items():
        source=Path(pin['path']);require(sha(source)==pin['sha256'],'saved pin '+alias)
        target=J.out/'saved'/alias;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        require(sha(target)==pin['sha256'],'saved copy '+alias)
        copied[alias]={'source':str(source),'path':str(target.relative_to(J.out)),'sha256':pin['sha256'],'bytes':target.stat().st_size}
        raw[alias]=json.loads(target.read_text())
    save(J.out/'saved-copy-index.json',copied)
    # Restore only complete old operation returns and their exact operands.
    receipts=m['uniformReceipts'];database=Path(m['uniformDatabase']);conn=sqlite3.connect('file:'+str(database)+'?mode=ro',uri=True)
    blobs={};opaque_copies={}
    try:
        for rec in receipts['blobReceipts']:
            payload,h,size=conn.execute('select payload,sha256,bytes from blobs where name=?',(rec['member'],)).fetchone()
            require(hashlib.sha256(payload).hexdigest()==h==rec['sha256'] and len(payload)==size==rec['bytes'],'opaque old blob '+rec['member'])
            target=J.out/'uniform-raw'/rec['member'];target.parent.mkdir(parents=True,exist_ok=True)
            with target.open('xb') as f:f.write(payload);f.flush();os.fsync(f.fileno())
            opaque_copies[rec['member']]={'path':str(target.relative_to(J.out)),'sha256':h,'bytes':size}
            value=codec(payload);blobs[rec['member']]=value
            J.emit('restore-opaque-%03d'%len(blobs),{'receipt':rec,'copy':str(target.relative_to(J.out)),'decodedType':str(type(value)),'oldFunctionCalled':False})
    finally:
        conn.close();save(J.out/'opaque-copy-index.json',opaque_copies)
    operations={r['name']:r for r in receipts['completedOperations']}
    def operation(name):
        r=operations[name];return blobs[r['input']['member']],blobs[r['result']['member']]
    restored={}
    for name,pin in m['originalPackets'].items():
        inp,value=operation('restore/'+name);J.join('restore-pin-'+name,inp['pin'],pin)
        require(inp['originalBytes']['sha256']==pin['sha256'] and inp['originalBytes']['bytes']==pin['bytes'],'actual original bytes '+name)
        require(sha(pin['path'])==pin['sha256'],'original source stays pinned '+name);restored[name]=value
    physical=json.loads(Path(m['physicalInput']).read_text());context=D(raw['ends/local-context.json']);cells=D(raw['ends/symbols.json']);depth=D(raw['ends/depth-domain.json'])
    require(context['physical']==physical and context['fullAssembly']['nativeEpsilonOnce'] is True and len(cells)==200,'inherited physical/epsilon/cell context')
    params={k:sp.Rational(v) for k,v in physical['parameters'].items()};params['omega']=sp.Integer(3)
    p,q,cs=depth['p'],depth['q'],depth['cs'];eta,sigma=context['context']['independentGrades']
    require((eta.name,sigma.name)==('eta_bg','sigma_W'),'actual independent grades')
    require(params['eta_bg']==sp.Rational(1,100) and params['W_0']==1 and params['L_W']==10,'actual finite origin and profile')
    require(raw['controls/conclusion.json']['status']=='SOURCE_JOINED_TRANSLATED_RETAINED_WEAK_ENDS' and raw['controls/conclusion.json']['endSymbolCells']==200,'saved concluded cells')
    J.zero('physical-wave-form',depth['radicand'],9/cs**2-p**2-sp.Rational(1,20))
    beta=depth['beta'];br,bi=sp.expand_complex(beta).as_real_imag()
    J.emit('actual-closed-depth-domain',{'source':depth,'betaReal':br,'betaImaginary':bi,'argument':'For outgoing real q>=0, Re(q+beta)>=Re(beta)>0; for q=i*y with y>=0, Re(q+beta)=Re(beta)>0.','newRootOrLimit':False})
    require(br==sp.Rational(30,109) and bi==sp.Rational(9,109) and depth['range']==[1,2],'actual beta/domain arguments')
    # Native coordinate covariance is new missing evidence, never a producer replay.
    tex={'ast':ast,'hashlib':hashlib,'Path':Path}
    exec(compile(definitions(Path(m['textHelperSource']).read_text(),('require','text_sha','literal_record','tuple_arguments','literal_key','named','selected_case')),'inert-source-text-only','exec'),tex)
    original_text,provenance=tex['literal_record'](Path(m['bSource']),'slab_operator')
    _,case_text=tex['selected_case'](original_text,('LAB_HELD','RHO4_CONSTANT'));value_text=tex['named'](case_text,'VALUE')
    actual_rows={}
    for key in ('U_BODY_BALANCE','THETA_BALANCE','E_W_BALANCE'):
        part=tex['named'](tex['named'](value_text,key),'EXPANDED')
        if key=='U_BODY_BALANCE':actual_rows.update({'U'+str(i):v for i,v in enumerate(tex['tuple_arguments'](part))})
        else:actual_rows[key]=part
    for row in ROWS:
        require(raw['native/'+row+'.json']['source']==provenance and raw['native/'+row+'.json']['fullConstructor']==actual_rows[row],'actual native raw source row '+row)
    J.emit('native-source-text-joins',{'provenance':provenance,'rowHashes':{k:hashlib.sha256(v.encode()).hexdigest() for k,v in actual_rows.items()},'sourceConstructorCalled':False})
    contracts=D(raw['native/wave-profile.json'])
    c2_text=Path(m['c2Source']).read_text();c2_tree=ast.parse(c2_text)
    assignments={t.id:n.value for n in c2_tree.body if isinstance(n,ast.Assign) for t in n.targets if isinstance(t,ast.Name)}
    require(contracts['dimensions']==ast.literal_eval(assignments['DIMENSION_SCHEMA']),'actual native unit schema')
    wave_node=next(n for n in c2_tree.body if isinstance(n,ast.FunctionDef) and n.name=='wave_jet')
    require(contracts['waveSource']==ast.get_source_segment(c2_text,wave_node),'actual native wave/derivative convention')
    native={row:D({'text':'native source operand','srepr':raw['native/'+row+'.json']['fullConstructor']}) for row in ROWS}
    atoms=set().union(*(v.free_symbols for v in native.values()));byname={s.name:s for s in atoms};require(len(byname)==len(atoms),'unique native assumptions')
    scalar_names=set(physical['parameters'])|{'epsilon_shape','sigma_W','delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus'}
    rotation={}
    for atom in atoms:
        name=rotated_name(atom.name,scalar_names);require(name in byname,'mapped native atom present '+name)
        target=byname[name];require(atom.assumptions0==target.assumptions0,'rotation assumptions');rotation[atom]=target
    rowmap={'U0':'U2','U1':'U0','U2':'U1','THETA_BALANCE':'THETA_BALANCE','E_W_BALANCE':'E_W_BALANCE'}
    J.emit('source-chart-operands',{'weakMomentum':(p,sp.Rational(1,5),sp.Rational(1,10)),'uniformMomentum':(sp.Rational(1,5),sp.Rational(1,10),p),'rotation':rotation,'rowMap':rowmap,'sourceRows':native,'nativeContract':contracts,'uniformChart':source_method(Path(m['engineSource']).read_text(),'EdgeReduction','__init__')})
    for row in ROWS:J.zero('source-cyclic-covariance-'+row,native[row].xreplace(rotation),native[rowmap[row]])
    S=sp.ImmutableMatrix([[0,0,1,0,0],[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    require(S.det()==1 and S.T*S==sp.eye(5),'proper orthonormal chart map')
    # Dimensions are checked on unbound native operands, not guessed after binding.
    dimensions=contracts['dimensions'];zero_dim=(sp.S.Zero,)*3
    def dimension(expr):
        if expr.is_Number or expr==sp.I:return zero_dim
        if expr.is_Symbol:
            require(expr.name in dimensions,'native dimensions '+expr.name);return tuple(map(sp.Rational,dimensions[expr.name]))
        if expr.is_Add:
            ds=[dimension(v) for v in expr.args if v!=0];require(not ds or all(d==ds[0] for d in ds),'native additive dimensions');return ds[0] if ds else zero_dim
        if expr.is_Mul:
            ds=[dimension(v) for v in expr.args];return tuple(map(sum,zip(*ds)))
        if expr.is_Pow and expr.exp.is_Rational:return tuple(expr.exp*d for d in dimension(expr.base))
        raise ValueError('unsupported native dimension node '+str(type(expr)))
    rowunits=[dimension(native[r]) for r in ROWS];fieldunits=[tuple(dimensions[f]) for f in FIELDS]
    J.emit('source-units',{'nativeRows':rowunits,'nativeFields':fieldunits,'frame':physical['unit_frame'],'unitMap':'proper permutation only, scalar factor one from strong epsilon coefficient; no power map'})
    lift_input,L=operation('selected-lift');U=restored['uniformSource'];common=restored['uniformCommon'];response=restored['uniformResponse']
    native_hashes=[h for path,h in U['sourceFiles'].items() if path.endswith('/S11c_b_exports.py')]
    require(native_hashes==[sha(m['bSource'])],'same original physical slab export across old/new sources')
    J.join('lift-input-curl',lift_input['uniform'],U['curl']);J.join('lift-input-common',lift_input['common'],common);J.join('lift-input-parameters',lift_input['parameters'],params)
    require(lift_input['source']==sha(m['engineSource']),'actual original lift source version')
    L={**L,'fieldUnits':response['fieldUnits'],'currentUnit':response['currentUnit']}
    J.join('lift-source-evidence',L['source'],blobs['lift/source-joins'])
    require(not any(s.name in ('eta_bg','sigma_W') for s in L['lift'].free_symbols),'grade-independent lift')
    require(tuple(L['source']['fieldOrder'])==('u1','u2','u3','theta','eW'),'source field order')
    expected_units=[fieldunits[i] for i in (1,2,0,3,4)];require(list(map(tuple,L['fieldUnits']))==expected_units,'native field unit permutation')
    J.emit('source-scalar-normalization',{'strongMatrixSource':source_method(Path(m['engineSource']).read_text(),'ConstantEndPencil','strong_matrix'),'pairingSource':source_method(Path(m['engineSource']).read_text(),'ClosedCurrentPairing','construct'),'epsilonCount':1,'factor':1,'excludedPowerMap':'PLUS_ROW_POWER_MAP','weakDuality':raw['native/weak-duality.json']})
    require(context['fullAssembly']['sourceRows']==list(ROWS) and context['fullAssembly']['fields']==list(FIELDS),'native weak row field identities')
    # Restore old proof operands and source profile values rather than calculating limits.
    profile=blobs['profiles/constant-end-and-zero-jets']
    for item in profile['profileJoins']:
        require(item['residual']==0 and all(v==0 for v in item['jets']),'inherited profile returns')
        require(item['supplied']==U['profileEndpoints'][item['source'][2]],'actual profile endpoint operand')
    groups={side:{g:sp.zeros(5) for g in G} for side in ('minus','plus')};closed={side:{g:sp.zeros(5) for g in G} for side in groups}
    coverage=[];seen=set()
    for index,c in enumerate(cells):
        side,row,col,g=c['side'],ROWS.index(c['row']),FIELDS.index(c['field']),tuple(c['grade']);key=(side,row,col,g)
        require(key not in seen and side in groups and g in G,'unique end cell');seen.add(key)
        require(c['sumIdentity']['cancelled']==0 and c['sumIdentity']['left']==c['symbol'],'saved sum proof actual argument')
        require(c['weakPairingFactor']==2*sp.pi and c['epsilonPower']==(0 if c['symbol']==0 else 1),'epsilon and Fourier contract')
        groups[side][g][row,col]=c['symbol'];closed[side][g][row,col]=c['closedGrazingValue']
    require(len(seen)==200,'all end cells')
    results={}
    for end,side in (('LEFT','minus'),('RIGHT','plus')):
        inp,R=operation(end+'/restriction');freq=restored[end+'Frequency'];pairing,dims=restored[end+'Pairing'];r=pairing['result']
        J.join(end+'-restriction-source',inp['actualSource'],r['CLOSED_PENCIL_LEGS']);J.join(end+'-restriction-frequency',inp['frequency'],freq)
        J.join(end+'-restriction-lift',inp['lift'],L);J.join(end+'-restriction-parameters',inp['parameters'],params)
        inv=blobs[end+'/restriction/full-five-row-invariance'];algebra=blobs[end+'/restriction/physical-algebraic-joins']
        for key in ('P','D','lift','gram'):J.join(end+'-invariant-'+key,R[key],inv[key])
        J.join(end+'-invariant-raw',R['invariant'],inv['rawResidual']);J.join(end+'-invariant-onwave',R['onWaveInvariant'],inv['onWaveResidual'])
        require(all(v==0 for v in R['onWaveInvariant']),'inherited on-wave zero only')
        for key in ('P','algebraic','relation','wave','conversion'):J.join(end+'-algebraic-'+key,R[key],algebra[key])
        atoms_R=set().union(*(x.free_symbols for x in (R['P'],R['D'],R['wave'],R['lift'],R['conversion'])))
        def symbol(name):
            found=[s for s in atoms_R if s.name==name];require(len(found)==1,'actual old symbol '+name);return found[0]
        w,k,cold,qold=symbol('uniformFrequency'),symbol('uniformNormal'),symbol('uniformSoundSpeed'),symbol('uniformPhysicalDepth')
        common_map={w:sp.Integer(3),k:p,cold:cs,qold:q}
        old=lambda v:v.xreplace(common_map)
        Pold,Dold,lift,Iold=map(old,(R['P'],R['D'],R['lift'],R['invariant']))
        Egrades={g:S.T*groups[side][g]*S for g in G};Ezero={g:S.T*closed[side][g]*S for g in G}
        origin=freq['origin'];require({s.name for s in origin}=={'eta_bg','sigma_W'},'old two independent grades')
        origins={s.name:v for s,v in origin.items()};require(origins=={'eta_bg':sp.Rational(1,100),'sigma_W':sp.Rational(1,1000)},'exact physical origin')
        finite_origin={eta:origins['eta_bg'],sigma:origins['sigma_W']}
        E=sum((eta**a*sigma**b*Egrades[a,b] for a,b in G),sp.zeros(5));Ephys=E.subs(finite_origin)
        Eg0=sum((eta**a*sigma**b*Ezero[a,b] for a,b in G),sp.zeros(5)).subs(finite_origin)
        # Raw binding is a missing new comparison; no old EndBinding method called.
        source=r['CLOSED_PENCIL_LEGS'][0].xreplace(U['profileEndpoints'])
        wl,wr=r['FREQUENCY_LEGS'];kl,kr=r['NORMAL_LEGS'];ql,qr=r['BULK_LEGS']
        mapping={wl:sp.Integer(3),wr:sp.Integer(3),kr:p,qr:q}
        for s in source.free_symbols-set(mapping):
            if s.name=='eta_bg':mapping[s]=eta
            elif s.name=='sigma_W':mapping[s]=sigma
            elif s.name=='c_s0':mapping[s]=cs
            elif s.name in ('omega','s11cdFrequency'):mapping[s]=sp.Integer(3)
            elif s.name in params:mapping[s]=params[s.name]
            else:raise ValueError('raw source binding unavailable '+s.name)
        Praw=source.xreplace(mapping);J.emit(end+'-raw-source-binding',{'nativeSource':r['CLOSED_PENCIL_LEGS'][0],'endpoints':U['profileEndpoints'],'source':source,'map':mapping,'raw':Praw,'origin':finite_origin,'oldPhysical':Pold,'halfHeights':raw['controls/conclusion.json']['heights']})
        require(not Praw.atoms(sp.Integral,sp.Derivative,sp.Limit,sp.Subs),'no unjoined original source operator')
        require(Praw.free_symbols<=set((p,q,cs,eta,sigma)),'complete raw binding')
        originjoin=inspect_matrix(Praw.subs(finite_origin)-Pold,p,q,cs,J,end+'/raw-origin')
        require(all(x['status']=='ZERO' for row in originjoin for x in row),'raw finite origin/source join')
        # Entry units derive from the unbound native rows and the inherited field frame.
        roworder=(1,2,0,3,4);native_entry_units=[[tuple(a-b for a,b in zip(rowunits[i],fieldunits[j])) for j in roworder] for i in roworder]
        require([[tuple(U['units']['strong'][(5*i+j,)]) for j in range(5)] for i in range(5)]==native_entry_units,'strong native row/field units')
        J.emit(end+'-unit-join',{'native':native_entry_units,'old':U['units']['strong'],'restoredRestrictionUnits':R['units'],'factor':1})
        action=Ephys*lift;A=(Ephys-Pold)*lift;Rnew=action-lift*Dold;Delta=Ephys-Pold
        for i in range(5):
            for j in range(2):J.zero(end+'-raw-AR-consistency-%d-%d'%(i,j),Rnew[i,j]-A[i,j],Iold[i,j])
        selected=inspect_matrix(A,p,q,cs,J,end+'/selected-A');residual=inspect_matrix(Rnew,p,q,cs,J,end+'/selected-R');full=inspect_matrix(Delta,p,q,cs,J,end+'/full-Delta')
        raw_action=Praw*lift;new_action=E*lift;grade=[]
        for i in range(5):
            for j in range(2):
                nm=end+'/grade-%d-%d'%(i,j);oldc=coefficients(raw_action[i,j],eta,sigma,J,nm+'-old');newc=coefficients(new_action[i,j],eta,sigma,J,nm+'-new')
                if not all(c['status']=='EXTRACTED_ON_DECLARED_REGULAR_DOMAIN' for c in (oldc,newc)):
                    grade.append({'row':i,'column':j,'status':'ATTRIBUTION_UNAVAILABLE','old':oldc,'new':newc});continue
                joins={g:newc['coefficients'][g]-oldc['coefficients'][g] for g in G}
                jt={g:wave_test(v,p,q,cs,J,nm+'-J-%d%d'%g) for g,v in joins.items()}
                attrib=sum(eta**a*sigma**b*joins[a,b] for a,b in G)-oldc['remainder']
                J.zero(nm+'-attribution',A[i,j],attrib.subs(finite_origin))
                require(newc['remainder']==0,'new retained law has no excluded remainder')
                retained='AGREEMENT' if all(v['status']=='ZERO' for v in jt.values()) else ('RETAINED_MISMATCH' if any(v['status']=='NONZERO_CERTIFIED' for v in jt.values()) else 'UNRESOLVED')
                finite=selected[i][j]['status'];classification='FINITE_TRUNCATION_DIFFERENCE' if retained=='AGREEMENT' and finite=='NONZERO_CERTIFIED' else retained
                grade.append({'row':i,'column':j,'retained':retained,'finite':finite,'classification':classification,'J':jt,'H':oldc['remainder'],'regularity':oldc['originDenominator']})
        for index,c in enumerate(cells):
            if c['side']!=side:continue
            wi=ROWS.index(c['row']);wj=FIELDS.index(c['field']);native_lift=S*lift
            weighted=[c['symbol']*native_lift[wj,j] for j in range(2)]
            coverage.append({'end':end,'cellIndex':index,'row':c['row'],'field':c['field'],'grade':c['grade'],'localAncestry':c['localAncestry'],'pressureAddresses':c['addressIds'],'selectedWeights':[native_lift[wj,j] for j in range(2)],'weightedSymbol':weighted,'zeroWeight':all(native_lift[wj,j]==0 for j in range(2)),'pressureSum':c['pressureSum'],'localSum':c['localSum']})
        # Closed grazing uses saved p and speed targets. Do not call any limit.
        grazing=[]
        for sign in (-1,1):
            args,limits=operation(end+'/exact-limits/'+str(sign));J.emit(end+'-limit-%d-actual-input'%sign,args)
            J.join(end+'-limit-%d-restriction'%sign,args['restriction'],R)
            native_args,native_return=operation(end+'/native-reconstruction')
            J.join(end+'-limit-%d-native'%sign,args['native'],native_return)
            J.join(end+'-native-%d-pairing'%sign,native_args['pairing'],pairing)
            J.join(end+'-native-%d-original'%sign,native_args['native'],restored[end+'NativeSource'])
            require(args['sign']==sign and native_args['frequency']==3,'actual limit/native arguments')
            require(limits['status']=='COMPUTED','saved limits completed')
            normal=old(limits['normal']);speed=old(R['cT']);mp={p:normal,cs:speed,q:sp.S.Zero}
            J.zero(end+'-grazing-wave-%d'%sign,(9/cs**2-p**2-sp.Rational(1,20)).subs(mp),sp.S.Zero)
            require(limits['normal']==sign*sp.sqrt(R['k23']),'actual inherited normal target')
            C=Eg0.subs(mp);LL=lift.subs(mp);DD=Dold.subs(mp);g=C*LL-LL*DD
            entries=[[finite_constant(g[i,j]) for j in range(2)] for i in range(5)]
            regular=all(finite_constant(x).get('finite',False) for obj in (C,LL,DD) for x in obj)
            paths={name:v['selected']['pencil'] for name,v in limits['paths'].items()}
            require(set(paths)=={'RADIATING','EVANESCENT'},'both actual limit paths')
            for route in paths:
                ev=blobs[end+'/limits/sign-'+str(sign)+'/'+route+'/wave-path']
                require(all(x==0 for x in ev['waveResiduals']),'saved actual wave path zero')
                t=ev['parameter'];require(t.is_positive is True,'actual approach parameter')
                require(ev['mapping'][qold]==(t if route=='RADIATING' else sp.I*t),'actual outgoing approach')
                require(ev['mapping'][k]==limits['normal'] and ev['mapping'][w]==3,'actual frequency/normal approach')
                J.emit(end+'-grazing-%d-'%sign+route+'-path',ev)
            require(all(v==0 for v in limits['selectedJoin']['pencil']),'inherited two-sided selected pencil join')
            J.emit(end+'-grazing-%d'%sign,{'point':mp,'closedE':C,'lift':LL,'D':DD,'residual':g,'entryCertificates':entries,'regular':regular,'savedSelectedPencilLimits':paths,'savedSelectedJoin':limits['selectedJoin'],'gradeAttribution':'UNAVAILABLE_NO_NEW_REMAINDER_LIMIT_COMPUTED','oldFunctionsCalled':False})
            require(regular,'closed grazing operands finite')
            grazing.append({'sign':sign,'point':mp,'R0':entries,'grazingGradeAttribution':'UNAVAILABLE'})
        # Three controls use full changed-minus-baseline row movements.
        point={p:sp.S.One,q:sp.Integer(2),cs:sp.sqrt(sp.Rational(180,101))};controls=[]
        choices=[(g,i,j,(eta**g[0]*sigma**g[1]*Egrades[g][i,j]).subs(finite_origin)) for g in G for i in range(5) for j in range(5) if Egrades[g][i,j]!=0 and any(lift[j,c]!=0 for c in range(2))]
        def record_control(name,baseline,changed,ancestry):
            movement=(changed-baseline).subs(point);certs=[finite_constant(x) for x in movement]
            baseline_certs=[finite_constant(x.subs(point)) for x in baseline];changed_certs=[finite_constant(x.subs(point)) for x in changed]
            responsive=any(c['status']=='NONZERO_CERTIFIED' for c in certs) and all(c.get('finite',False) for c in certs+baseline_certs+changed_certs)
            rec={'name':name,'point':point,'baseline':baseline,'changed':changed,'movement':movement,'certificates':certs,'baselineCertificates':baseline_certs,'changedCertificates':changed_certs,'responsive':responsive,'ancestry':ancestry,'meaning':'sensitivity only'};J.emit(end+'/control-'+name,rec);controls.append(rec);return responsive
        chosen=None
        for g,i,j,term in choices:
            move=sp.zeros(5,2)
            for c in range(2):move[i,c]=-term*lift[j,c]
            if any(finite_constant(x.subs(point))['status']=='NONZERO_CERTIFIED' for x in move):chosen=(g,i,j,term,move);break
        if chosen:
            g,i,j,term,move=chosen
            original=next(c for c in cells if c['side']==side and c['row']==ROWS[roworder[i]] and c['field']==FIELDS[roworder[j]] and tuple(c['grade'])==g)
            record_control('omit-cell',Rnew,Rnew+move,{'row':i,'field':j,'grade':g,'removedTerm':term,'originalSavedCell':original})
        else:controls.append({'name':'omit-cell','responsive':False,'status':'INAPPLICABLE'})
        bad=lift.xreplace({p:-p});responsive=record_control('reverse-lift-normal',Rnew,Ephys*bad-bad*Dold,{'normalSignHeldInEquation':True})
        if not responsive:
            for i,j in ((i,j) for i in range(5) for j in range(2) if lift[i,j]!=0):
                bad=sp.MutableDenseMatrix(lift);bad[i,j]=0;movement=Ephys*bad-bad*Dold-Rnew
                if any(finite_constant(x.subs(point))['status']=='NONZERO_CERTIFIED' for x in movement):
                    record_control('curl-entry-ablation',Rnew,Ephys*bad-bad*Dold,{'row':i,'column':j,'directionProof':False});break
        sheet_choice=None
        for j in (3,4):
            baseline=Ephys[:,j];changed=Ephys.xreplace({q:-q})[:,j]
            require(all(not Egrades[g][i,j].has(q) or any(c['side']==side and c['field']==FIELDS[roworder[j]] and c['pressureSum']!=0 for c in cells) for g in G for i in range(5)),'sheet movement has actual pressure ancestry')
            if any(finite_constant(x.subs(point))['status']=='NONZERO_CERTIFIED' for x in changed-baseline):sheet_choice=(j,baseline,changed);break
        if sheet_choice:
            j,baseline,changed=sheet_choice;record_control('opposite-sheet',baseline,changed,{'unitColumn':j,'isMode':False,'pressureAncestry':[c for c in cells if c['side']==side and c['field']==FIELDS[roworder[j]] and c['pressureSum']!=0]})
        else:controls.append({'name':'opposite-sheet','responsive':False,'status':'INAPPLICABLE'})
        result={'fullDifference':full,'selectedFinite':selected,'fiveRowResidual':residual,'grades':grade,'grazing':grazing,'controls':controls,'currentTransfer':False,'comparisonDoesNotDischargeZeroWeightedPressureCells':True}
        J.emit(end+'-comparison',result);results[end]={'retainedStatuses':sorted(set(g.get('retained',g.get('status','UNRESOLVED')) for g in grade)),'finiteStatuses':sorted(set(v['status'] for row in selected for v in row)),'controlCoverage':{c['name']:c['responsive'] for c in controls}}
    J.emit('selected-cell-coverage',coverage)
    return {'executionStatus':'COMPLETE_COMPARISON_PENDING_INSPECTION','ends':results,'cells':len(cells),'restoredOperations':len(operations),'newRoots':False,'fieldOrLoss':False,'frequencyDerivativeTransfer':False,'oldFunctionsCalled':False}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--gate',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    m=json.loads(args.inputs.read_text());g=verify_gate(args.gate,args.inputs,m)
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(sys.argv==expected and g['command'][-len(expected):]==expected and args.out.resolve()==Path(g['outputDirectory']).resolve(),'actual argv/output')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;code=1;result={};start=time.monotonic()
    pins={**m['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(m['helperSource']).read_text(),('require','containment','decode')),'inert-contained-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp,np,FunctionClass
        import sympy as sp
        import numpy as np
        from sympy.core.symbol import Str
        from sympy.core.function import FunctionClass
        ns.update(sp=sp,Str=Str)
        codec_ns={'pickle':pickle,'OrderedDict':OrderedDict,'builtins':builtins,'np':np,'sp':sp,'importlib':importlib,'io':io,'require':require}
        exec(compile(definitions(Path(m['uniformWorker']).read_text(),('SavedCodec','decode')),'saved-codec-only-no-old-guard','exec'),codec_ns)
        J=Evidence(args.out);result=J.stage('end-uniform-comparison',{'manifest':g['manifestSha256'],'build':g['buildReviewRecordSha256']},lambda:run_science(m,J,ns,codec_ns['decode']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        for index_name in ('saved-copy-index.json','opaque-copy-index.json'):
            if (args.out/index_name).exists():
                for row in json.loads((args.out/index_name).read_text()).values():pins[str(args.out/row['path'])]=row['sha256']
        posts={}
        for p,h in pins.items():
            try:posts[p]={'expected':h,'actual':sha(p),'error':None}
            except OSError as e:posts[p]={'expected':h,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',posts)
        if any(r['expected']!=r['actual'] for r in posts.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        if J is not None:result=J.encode(result)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
        if J is not None:save(args.out/'evidence-final-receipt.json',{'records':J.count,'lastSha256':J.previous,'chainSha256':sha(args.out/'evidence-chain.jsonl'),'completeProcess':code==0})
    return code

if __name__=='__main__':sys.exit(main())
