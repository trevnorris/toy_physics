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
    if fraction is None:
        state=raw.get(prefix+'-state.json')
        require(state is not None and state.get('finite') is True and state.get('nonzero') is True,'literal original regularity flags')
        evidence={'input':inp,'state':state}
    else:
        require(all(v is True for row in fraction['finite'] for v in row),'saved exact component finiteness')
        require(all(any(v is True for v in row) for row in fraction['signedNonzero']),'saved signed nonzero numerator and denominator')
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


def profile_transport(record,context,registry):
    require(record['L']==context['numeric']['L_W'] and record['independentSigma'] is True and record['profileTransformEvaluated'] is False,'actual original profile context')
    entries=[];seen=set()
    for original,mapped in record['map']:
        name=symbol_name(parse(original['srepr']));match=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',name)
        require(match is not None and name not in seen,'unique original profile atom');seen.add(name)
        directions=re.findall('d([123])',match[2]);order=len(directions);transverse=any(d!='1' for d in directions)
        require(name in registry and unit_tuple(registry[name])==ZERO,'original dimensionless profile jet')
        require(not transverse or zero(mapped),'actual saved transverse profile zero')
        formal=add(scale(unit_tuple(registry['L_W']),order),(-order,0,0))
        require(formal==ZERO,'native L^r physical-derivative unit')
        entries.append({'original':original,'savedMapped':mapped,'derivativeOrder':order,'directions':directions,'transverseZero':transverse,'LExponent':order,'formalMappedUnit':formal,'numericDerivativeReevaluated':False})
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
