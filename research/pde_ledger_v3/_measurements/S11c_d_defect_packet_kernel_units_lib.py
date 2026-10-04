"""Inert source-AST dimensional interpreter. Native calls require containment.

No numerical function, CAS constructor, transform or integral is evaluated.
Tests use manufactured source and units only.
"""
import ast, hashlib, importlib.util
from pathlib import Path
from fractions import Fraction as F

_base=Path(__file__).with_name('S11c_d_defect_packet_source_units_lib.py')
_spec=importlib.util.spec_from_file_location('inherited_source_unit_helpers',_base)
U=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(U)
require=U.require
ZERO=(F(0),F(0),F(0)); LENGTH=(F(1),F(0),F(0)); MOMENTUM=(-F(1),F(0),F(0))
add=U.add;scale=U.scale;unit=U.unit_tuple


def same(a,b):return ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
def expression(text):return ast.parse(text,mode='eval').body


def rational(n):
    if isinstance(n,ast.Constant) and type(n.value) is int:return F(n.value)
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -rational(n.operand)
    if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Div):return rational(n.left)/rational(n.right)
    if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='sp' and n.func.attr=='Rational' and len(n.args)==2 and not n.keywords:return rational(n.args[0])/rational(n.args[1])
    raise ValueError('only literal exact rational exponent')


class Walk:
    def __init__(self,env,functions=None,params=None):
        self.env={k:None if v is None else unit(v) for k,v in env.items()}
        self.functions=functions or {};self.params=params or {};self.events=[]
    def dim(self,n,path='root'):
        if isinstance(n,str):n=expression(n)
        try:result=self._dim(n,path)
        except BaseException as error:
            self.events.append({'path':path,'source':ast.unparse(n),'refused':True,'reason':str(error)});raise
        self.events.append({'path':path,'source':ast.unparse(n),'unit':result,'literalZero':result is None});return result
    def _dim(self,n,p):
        if isinstance(n,ast.Constant):
            require(type(n.value) is int,'no floats/strings/bools in dimensional arithmetic');return None if n.value==0 else ZERO
        if isinstance(n,ast.Name):
            require(n.id in self.env,'unbound dimensional name '+n.id);return self.env[n.id]
        if isinstance(n,ast.Attribute):
            require(isinstance(n.value,ast.Name) and n.value.id=='sp' and n.attr in ('I','pi'),'unsupported dimensional attribute');return ZERO
        if isinstance(n,ast.Subscript):
            require(isinstance(n.value,ast.Name) and n.value.id=='params' and isinstance(n.slice,ast.Constant) and type(n.slice.value) is str and n.slice.value in self.params,'unbound native parameter');return unit(self.params[n.slice.value])
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):return self.dim(n.operand,p+'.operand')
        if isinstance(n,ast.BinOp):
            a=self.dim(n.left,p+'.left')
            if isinstance(n.op,ast.Pow):
                e=rational(n.right);require(not(a is None and e<=0),'undefined zero power');return None if a is None else scale(a,e)
            b=self.dim(n.right,p+'.right')
            if isinstance(n.op,(ast.Add,ast.Sub)):
                require(a is None or b is None or a==b,'inhomogeneous source addition');return b if a is None else a
            if isinstance(n.op,ast.Mult):return None if a is None or b is None else add(a,b)
            if isinstance(n.op,ast.Div):
                require(b is not None,'literal zero denominator');return None if a is None else add(a,scale(b,-1))
            raise ValueError('unsupported binary operator')
        if isinstance(n,ast.Call):
            require(not n.keywords,'no dimensional keyword calls')
            if isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='sp':
                if n.func.attr in ('cancel','expand'):
                    require(len(n.args)==1,'unary inert algebra wrapper');return self.dim(n.args[0],p+'.argument')
                if n.func.attr in ('sinh','tanh','cosh','exp'):
                    require(len(n.args)==1 and self.dim(n.args[0],p+'.argument') in (ZERO,None),'dimensionless transcendental argument');return ZERO
                if n.func.attr in ('Rational','Integer'):
                    require(len(n.args)==(2 if n.func.attr=='Rational' else 1),'literal constructor arity')
                    value=rational(n) if n.func.attr=='Rational' else rational(n.args[0]);return None if value==0 else ZERO
            require(isinstance(n.func,ast.Name) and n.func.id in self.functions,'unsupported dimensional function')
            fn=self.functions[n.func.id]
            require(len(n.args)==len(fn['arguments']),'complete dimensional function arity')
            dims=[self.dim(x,p+'.argument'+str(i)) for i,x in enumerate(n.args)]
            if 'body' in fn:
                child=Walk({**self.env,**dict(zip(fn['arguments'],dims))},self.functions,self.params)
                try:return child.dim(expression(fn['body']),p+'.lambda')
                finally:self.events.extend(child.events)
            require(dims==[unit(x) for x in fn['inputUnits']],'actual function argument units');return unit(fn['outputUnit'])
        raise ValueError('unsupported source AST '+type(n).__name__)


def contract(record):
    text=Path(record['source']).read_text();tree=ast.parse(text)
    nodes=[n for n in ast.walk(tree) if getattr(n,'lineno',None)==record['line'] and getattr(n,'col_offset',None)==record['column'] and type(n).__name__==record['kind']]
    require(len(nodes)==1 and ast.get_source_segment(text,nodes[0])==record['text'],'exact original source segment '+record['name'])
    require(hashlib.sha256(record['text'].encode()).hexdigest()==record['sha256'],'exact source segment hash')
    return nodes[0]


def assignment_expr(record):
    n=contract(record)
    require(isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name),'one source assignment')
    return n.targets[0].id,n.value


def statements(fn,env,functions=None,params=None,walk=None):
    require(isinstance(fn,ast.FunctionDef),'function AST only')
    walk=walk if walk is not None else Walk(env,functions,params);outputs=None
    for n in fn.body:
        if isinstance(n,ast.Expr) and isinstance(n.value,ast.Constant) and type(n.value.value) is str:continue
        if isinstance(n,ast.Assign):
            require(len(n.targets)==1 and isinstance(n.targets[0],ast.Name),'scalar assignment only')
            walk.env[n.targets[0].id]=walk.dim(n.value,'assignment.'+n.targets[0].id)
        elif isinstance(n,ast.Return):
            require(isinstance(n.value,ast.List),'explicit kernel return list');outputs=[walk.dim(x,'return.'+str(i)) for i,x in enumerate(n.value.elts)]
        else:raise ValueError('unexpected kernel statement '+type(n).__name__)
    require(outputs is not None,'kernel return present');return outputs,walk


def template_nodes(record):
    n=contract(record);require(isinstance(n,ast.Assign) and isinstance(n.value,ast.Subscript) and isinstance(n.value.value,ast.Dict),'actual component dictionary')
    result={}
    for k,v in zip(n.value.value.keys,n.value.value.values):
        require(isinstance(k,ast.Constant) and type(k.value) is str and k.value not in result,'unique component');result[k.value]=v
    require(len(result)==5,'five original response components');return result


def inherit_zero(args,ret):
    require(ret['cancelled']=={'text':'0','srepr':'Integer(0)'},'literal inherited zero')
    for side in ('left','right'):
        require(type(args[side]) is dict and 'srepr' in args[side], 'complete inherited operand '+side)
        if side in ret:require(U.equal(ret[side],args[side]),'original inherited return operand '+side)
    return {'input':args,'return':ret,'functionCalled':False,'interpretation':'Accepted identity on its original domain; opposite sides need not be structurally equal.'}


def total_unit(source,jet,consumer,kernel,measure):
    physical_source=add(unit(source),unit(jet));X=add(physical_source,LENGTH);Y=add(unit(consumer),LENGTH)
    return {'physicalSource':physical_source,'X':X,'Y':Y,'kernel':unit(kernel),'measure':unit(measure),'total':add(add(X,Y),add(unit(kernel),unit(measure)))}
