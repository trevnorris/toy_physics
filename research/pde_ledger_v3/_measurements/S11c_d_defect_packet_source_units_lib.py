"""Inert exact source/receipt helpers and bounded unit-transport bookkeeping.

Native calls occur only after containment. No SymPy import or constructor execution.
Synthetic tests may call these helpers on manufactured inputs.
"""
import ast
import hashlib
import importlib.util
import re
from pathlib import Path
from fractions import Fraction as F

_BASE=Path(__file__).with_name('S11c_d_defect_packet_native_units_lib.py')
_spec=importlib.util.spec_from_file_location('inherited_native_unit_engine',_BASE)
_U=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_U)
for _name in ('require','parse','kind','signature','same','rational','symbol_name','tuple_args','unit_tuple','add','scale','Units','ZERO'):
    globals()[_name]=getattr(_U,_name)


def equal_text(a,b):return same(parse(a),parse(b))
def equal(a,b):return type(a) is dict and type(b) is dict and 'srepr' in a and 'srepr' in b and equal_text(a['srepr'],b['srepr'])
def zero(a):return equal(a,{'srepr':'Integer(0)'})
def symbols(text):return {symbol_name(n) for n in ast.walk(parse(text)) if kind(n)=='Symbol'}
def tuple_leaf(text,i):return ast.unparse(tuple_args(parse(text))[i])
def encoded_leaf(value,i):return {'srepr':tuple_leaf(value['srepr'],i)}


def matrix_unit(value):
    node=parse(value['srepr']);require(kind(node) in ('MutableDenseMatrix','ImmutableDenseMatrix') and len(node.args)==1,'saved unit vector')
    rows=node.args[0];require(isinstance(rows,ast.List) and len(rows.elts)==3,'three saved unit rows')
    require(all(isinstance(r,ast.List) and len(r.elts)==1 for r in rows.elts),'one unit column')
    return tuple(rational(r.elts[0]) for r in rows.elts)


def reciprocal_bases(text,events):
    bypath={v['address']:v for v in events};require(len(bypath)==len(events),'unique inherited unit-walk paths')
    result=[]
    for entry in events:
        n=parse(entry['constructor'])
        if kind(n)=='Pow' and rational(n.args[1])<0:
            path=entry['address']+'/base';require(path in bypath,'inherited reciprocal base evidence')
            base=bypath[path];require(same(parse(base['constructor']),n.args[0]) and base['unit'] is not None,'same finite-dimensional unbound reciprocal base')
            result.append({'address':entry['address'],'baseConstructor':base['constructor'],'baseUnit':base['unit'],'power':str(rational(n.args[1])),'unitProofInherited':True,'nonzeroDomainProvedHere':False})
    require(equal_text(text,events[-1]['constructor']),'whole reciprocal census root')
    return result


def regularity(raw,prefix,value,J):
    inp=raw[prefix+'-input.json'];require(equal(inp['value'],value),'actual regularity argument')
    fraction=raw.get(prefix+'-fraction.json')
    require(fraction is not None,'complete inherited fraction/component evidence required; flags alone are insufficient')
    require(len(fraction['finite'])==len(fraction['signedNonzero'])==2 and all(len(v)==2 for v in fraction['finite']+fraction['signedNonzero']),'two real/imaginary component pairs')
    require(all(v is True for row in fraction['finite'] for v in row),'saved exact component finiteness')
    # A complex number is nonzero if at least one real/imaginary component is
    # provably nonzero; the other component may be exactly zero.
    require(all(type(v) is bool for row in fraction['signedNonzero'] for v in row) and all(any(v is True for v in row) for row in fraction['signedNonzero']),'saved signed nonzero numerator and denominator')
    evidence={'input':inp,'fraction':fraction,'joins':{}}
    for suffix in ('fraction-reconstruction','numerator-components','denominator-components'):
        a=raw[prefix+'-'+suffix+'-input.json'];r=raw[prefix+'-'+suffix+'-return.json']
        require(zero(r['cancelled']),'saved exact regularity join')
        expected=value if suffix=='fraction-reconstruction' else fraction['numerator' if suffix=='numerator-components' else 'denominator']
        require(equal(a['left'],expected),'same regularity operand')
        evidence['joins'][suffix]={'input':a,'return':r}
    J.emit('regularity-'+prefix.replace('/','-'),evidence)
    return evidence


def verify_contracts(manifest,contracts):
    for entry in contracts['fragments']:
        text=Path(entry['source']).read_text();tree=ast.parse(text)
        matches=[n for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==entry['name'] and n.lineno==entry['line']]
        require(len(matches)==1 and ast.get_source_segment(text,matches[0])==entry['text'],'actual original source contract '+entry['name'])
        require(hashlib.sha256(entry['text'].encode()).hexdigest()==entry['sha256'],'source contract digest')
    require({'grade_split','quotient_recurrence','profile','join_source_input','bind','wave_jet','at_source'}<=set(v['name'] for v in contracts['fragments']),'complete native/transport algorithm contracts')


def jet_unit(spec,registry,atom):
    name=symbol_name(parse(atom['srepr']));require(spec['name']==name and spec['channel']=='e_W','actual scalar jet name/channel')
    require(type(spec['timeOrder']) is int and spec['timeOrder']>=0 and len(spec['spatialOrders'])==3 and all(type(v) is int and v>=0 for v in spec['spatialOrders']),'native nonnegative derivative orders')
    require(name in registry,'native jet exists in original schema')
    unit=(-sum(spec['spatialOrders']),-spec['timeOrder'],0)
    require(unit_tuple(registry[name])==unit,'actual native scalar-jet registry and derivative rule')
    return unit


def profile_transport(record,context,registry,verifier):
    require(record['L']==context['numeric']['L_W'] and record['independentSigma'] is True and record['profileTransformEvaluated'] is False,'actual original profile context')
    entries=[];seen=set()
    for original,mapped in record['map']:
        name=symbol_name(parse(original['srepr']));match=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',name)
        require(match is not None and name not in seen,'unique original profile atom');seen.add(name)
        directions=re.findall('d([123])',match[2]);order=len(directions);transverse=any(d!='1' for d in directions)
        require(name in registry and unit_tuple(registry[name])==ZERO,'original dimensionless profile jet')
        require(not transverse or zero(mapped),'actual saved transverse profile zero')
        actual_scale=None if transverse else verifier.verify(match[1],order,mapped)
        formal=add(scale(unit_tuple(registry['L_W']),order),(-order,0,0))
        require(formal==ZERO,'native L^r physical-derivative unit')
        entries.append({'original':original,'savedMapped':mapped,'derivativeOrder':order,'directions':directions,'transverseZero':transverse,'LExponent':order,'formalMappedUnit':formal,'numericDerivativeReevaluated':False,'actualScaleCertificate':actual_scale})
    original_profiles={s for s in symbols(record['original']['srepr']) if re.fullmatch(r'[wm]1_profile(?:_?d[123])*',s)}
    require(seen==original_profiles,'all actual coefficient profile atoms mapped')
    return {'original':record['original'],'savedField':record['oneDimensional'],'entries':entries,'actualSourceAlgorithmPinned':True,'profileComputationReplayed':False}


def strip_epsilon(value,epsilon):
    if zero(value):return 'Integer(0)'
    n=parse(value['srepr']);e=parse(epsilon['srepr']);args=list(n.args) if kind(n)=='Mul' else [n]
    hits=[i for i,v in enumerate(args) if same(v,e)];require(len(hits)==1,'one explicit saved epsilon factor')
    del args[hits[0]]
    require(all(symbol_name(x)!='epsilon_shape' for a in args for x in ast.walk(a) if kind(x)=='Symbol'),'epsilon removed exactly once')
    if not args:return 'Integer(1)'
    if len(args)==1:return ast.unparse(args[0])
    return 'Mul('+', '.join(ast.unparse(a) for a in args)+')'


def field_proof(poly,args,ret,field):
    require(equal(args['left'],field) and zero(ret['cancelled']),'actual original field reconstruction input and literal return')
    require(type(poly['degree']) is int and poly['degree']>=0 and len(poly['coefficients'])==poly['degree']+1,'complete original field quotient vector')
    require(not symbols(poly['denominator']['srepr']),'saved constant polynomial denominator')
    require(not zero(poly['denominator']),'nonzero original polynomial denominator')
    require(all(not symbols(v['srepr']) for v in poly['coefficients']),'saved constant field coefficient vector')
    # Coefficient arithmetic is already accepted. The exact original input/return/vector files are pinned together.
    return {'fieldReconstructionInherited':True,'numericPolynomialIndependentlyDimensional':False}

# Exact finite ring normal forms for NEW operand-assembly/covariance joins only.
# No quotient coefficients, old source function, old derivative or CAS constructor
# is evaluated. Negative powers of non-monomial polynomials remain opaque atoms.
def cadd(a,b):return (a[0]+b[0],a[1]+b[1])
def cmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cpow(a,n):
    if n<0:
        den=a[0]*a[0]+a[1]*a[1];require(den!=0,'exact constant nonzero denominator');a=(a[0]/den,-a[1]/den);n=-n
    out=(F(1),F(0))
    while n:
        if n%2:out=cmul(out,a)
        a=cmul(a,a);n//=2
    return out

def radd(*polys):
    out={}
    for p in polys:
        for k,v in p.items():out[k]=cadd(out.get(k,(F(0),F(0))),v)
    return {k:v for k,v in out.items() if v!=(0,0)}

def rmul(a,b):
    out={}
    for ma,ca in a.items():
        for mb,cb in b.items():
            powers={}
            for atom,n in ma+mb:powers[atom]=powers.get(atom,0)+n
            key=tuple(sorted((atom,n) for atom,n in powers.items() if n));out[key]=cadd(out.get(key,(F(0),F(0))),cmul(ca,cb))
    return {k:v for k,v in out.items() if v!=(0,0)}

def rpower(poly,n):
    require(type(n) is int,'integer assembly power')
    if not poly:
        require(n>0,'nonzero assembly denominator');return {}
    if n<0:
        if len(poly)==1:
            mon,c=next(iter(poly.items()));return {tuple((a,p*n) for a,p in mon):cpow(c,n)}
        pivot=sorted(poly.items())[0][1];normalized={k:cmul(c,cpow(pivot,-1)) for k,c in poly.items()}
        key=repr(tuple(sorted(normalized.items())));return {(('inverse-base:'+key,n),):cpow(pivot,n)}
    out={(): (F(1),F(0))}
    while n:
        if n%2:out=rmul(out,poly)
        poly=rmul(poly,poly);n//=2
    return out

def ring(node):
    k=kind(node)
    if k in ('Integer','Rational'):
        v=rational(node);return {} if v==0 else {():(v,F(0))}
    if isinstance(node,ast.Name) and node.id=='I':return {():(F(0),F(1))}
    if k=='Symbol':symbol_name(node);return {((signature(node),1),):(F(1),F(0))}
    if k=='Add':return radd(*(ring(v) for v in node.args))
    if k=='Mul':
        out={(): (F(1),F(0))}
        for v in node.args:out=rmul(out,ring(v))
        return out
    if k=='Pow':
        n=rational(node.args[1]);require(n.denominator==1,'integer assembly exponent');return rpower(ring(node.args[0]),int(n))
    raise ValueError('unsupported exact assembly constructor '+str(k))

def dump_ring(p):return [{'monomial':list(k),'real':str(c[0]),'imaginary':str(c[1])} for k,c in sorted(p.items())]
def expr_call(name,*nodes):return ast.Call(func=ast.Name(id=name,ctx=ast.Load()),args=list(nodes),keywords=[])
def integer_node(n):return parse('Integer('+str(n)+')')
def times(*nodes):return expr_call('Mul',*nodes)
def plus(*nodes):return expr_call('Add',*nodes)
def negative(node):return times(integer_node(-1),node)


def constructor_substitute(node,mapping):
    # A new typed-unit provenance join: annotate/replace original leaves only.
    # The original bind/producer functions are not invoked.
    if kind(node)=='Symbol':return parse(mapping[symbol_name(node)]['srepr']) if symbol_name(node) in mapping else node
    if isinstance(node,ast.Call):return ast.Call(func=node.func,args=[constructor_substitute(v,mapping) for v in node.args],keywords=node.keywords)
    if isinstance(node,ast.UnaryOp):return ast.UnaryOp(op=node.op,operand=constructor_substitute(node.operand,mapping))
    return node


def binding_operand(text,context):
    node=constructor_substitute(parse(text),context['densityMap']);steps=[{'kind':'live-density','constructor':ast.unparse(node)}]
    for _ in range(len(context['profileEqualities'])+2):
        new=constructor_substitute(node,context['profileEqualities'])
        if same(new,node):break
        node=new;steps.append({'kind':'background-substitution','constructor':ast.unparse(node)})
    else:raise ValueError('source binding cycle')
    node=constructor_substitute(node,context['numeric']);steps.append({'kind':'original-parameter-magnitudes','constructor':ast.unparse(node),'physicalUnitsNotErasedFromProvenance':True})
    return node,steps


def full_remainder_operands(full,split,grades):
    e,s=[parse(v['srepr']) for v in grades];terms=[]
    for a,b in ((0,0),(1,0),(0,1),(1,1)):
        terms.append(times(parse(split['retained'][str((a,b))]['srepr']),expr_call('Pow',e,integer_node(a)),expr_call('Pow',s,integer_node(b))))
    retained=plus(*terms)
    return {'fullHigherRemainder':plus(parse(full['srepr']),negative(retained)),
            'quotientRingNumeratorRemainder':plus(parse(split['numerator']['srepr']),negative(times(parse(split['denominator']['srepr']),retained)))}


def at_grade_zero(node,grades):return constructor_substitute(node,{symbol_name(parse(v['srepr'])):{'srepr':'Integer(0)'} for v in grades})


def profile_polynomial(mapped,context):
    # Parse the ACTUAL saved mapped value. tanh's entire physical argument must
    # be x / the original L; no coefficient is identified by its numeric value.
    length=context['numeric']['L_W'];expected=ring(times(parse("Symbol('composition_x',real=True)"),expr_call('Pow',parse(length['srepr']),integer_node(-1))))
    T=parse("Symbol('unit_transport_T',real=True)");argument_joins=[]
    def visit(n):
        if kind(n)=='tanh':
            require(len(n.args)==1,'one tanh argument');actual=ring(n.args[0]);argument_joins.append({'actual':ast.unparse(n.args[0]),'expected':'composition_x / saved L_W','left':dump_ring(actual),'right':dump_ring(expected)})
            require(actual==expected,'actual profile x/L argument');return T
        if kind(n)=='Symbol':raise ValueError('unmatched symbol in saved profile map')
        if isinstance(n,ast.Call):return ast.Call(func=n.func,args=[visit(v) for v in n.args],keywords=n.keywords)
        if isinstance(n,ast.UnaryOp):return ast.UnaryOp(op=n.op,operand=visit(n.operand))
        return n
    poly=ring(visit(parse(mapped['srepr'])));key=signature(T);coeff={}
    for mon,value in poly.items():
        require(not mon or len(mon)==1 and mon[0][0]==key and mon[0][1]>=0,'univariate saved profile polynomial')
        degree=0 if not mon else mon[0][1];require(value[1]==0,'real native profile coefficients');coeff[degree]=value[0]
    return coeff,argument_joins


def profile_recurrence(previous):
    # NEW scale-covariance identity: (1-T^2) dP/dT, checked against the next
    # independently saved profile value. No old physical derivative is called.
    out={}
    for n,c in previous.items():
        if n:out[n-1]=out.get(n-1,F(0))+n*c;out[n+1]=out.get(n+1,F(0))-n*c
    return {n:c for n,c in out.items() if c}

def dump_coeff(p):return [{'power':n,'coefficient':str(c)} for n,c in sorted(p.items())]


def scale_catalog(records,context):
    catalog={};origins={}
    for alias,record in records.items():
        for original,mapped in record['map']:
            name=symbol_name(parse(original['srepr']));match=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',name);require(match is not None,'original profile atom')
            directions=re.findall('d([123])',match[2])
            if any(v!='1' for v in directions):require(zero(mapped),'actual transverse profile zero');continue
            key=match[1],len(directions)
            if key in catalog:require(equal(catalog[key],mapped),'all repeated original map values identical')
            catalog[key]=mapped;origins.setdefault(key,[]).append({'record':alias,'original':original,'mapped':mapped})
    return catalog,origins


def scale_certificate(base,order,mapped,catalog,context):
    # The base definitions are pinned original setup assignments, not guessed
    # normalizations. m's zeroth field need not occur in a selected coefficient.
    initial={'w':{0:F(1,2),1:F(1,2)},'m':{0:F(1,3),2:F(-1,3)}}[base]
    actual,args=profile_polynomial(mapped,context)
    if order==0:expected=initial;previous=None
    else:
        if order==1:previous=initial
        else:
            require((base,order-1) in catalog,'saved lower derivative order');previous,_=profile_polynomial(catalog[base,order-1],context)
        expected=profile_recurrence(previous)
    return {'base':base,'order':order,'actualMapped':mapped,'previousSaved':catalog.get((base,order-1)),'previousCoefficients':None if previous is None else dump_coeff(previous),'actualCoefficients':dump_coeff(actual),'expectedCoefficients':dump_coeff(expected),'argumentJoins':args,'exactCoefficientResidual':dump_coeff({n:actual.get(n,0)-expected.get(n,0) for n in set(actual)|set(expected) if actual.get(n,0)!=expected.get(n,0)}),'matches':actual==expected,'newIdentity':'P_r(T) = (1-T^2) P_(r-1)\u2032(T); chain rule gives L^r physical derivative','priorDerivativeCalled':False}


class ScaleVerifier:
    def __init__(self,catalog,context,J):self.catalog=catalog;self.context=context;self.J=J;self.cache={};self.count=0
    def verify(self,base,order,mapped):
        key=(base,order,signature(parse(mapped['srepr'])))
        if key in self.cache:return self.cache[key]
        name='profile-scale-certificate-'+str(self.count);self.count+=1
        self.J.emit(name+'-input',{'base':base,'order':order,'actualMapped':mapped,'originalLength':self.context['numeric']['L_W'],'previousSaved':self.catalog.get((base,order-1))})
        certificate=scale_certificate(base,order,mapped,self.catalog,self.context)
        self.J.emit(name+'-decision-operands',certificate)
        require(certificate['matches'],'actual mapped profile fails native L-scaled recurrence')
        result={'name':name,'inputKey':list(key),'certificate':certificate};self.cache[key]=result
        self.J.emit(name+'-return',result);return result
