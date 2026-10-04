"""Strict constructor-text unit interpreter; no symbolic evaluation or producers.

Actual source unit work is scientific and only called within the shared guard.
Tests use synthetic constructor strings and unit registries.
"""
import ast
from fractions import Fraction as F

ZERO=(F(0),)*3

def require(value,message):
    if value is not True:raise ValueError(message)

def parse(text):return ast.parse(text,mode='eval').body

def kind(node):return node.func.id if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) else None

def signature(node):return ast.dump(node,include_attributes=False)

def same(a,b):return signature(a)==signature(b)

def rational(node):
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -rational(node.operand)
    if isinstance(node,ast.Constant) and type(node.value) is int:return F(node.value)
    k=kind(node)
    require(k in ('Integer','Rational') and not node.keywords,'literal exact rational')
    args=[rational(v) for v in node.args]
    require(len(args)==(1 if k=='Integer' else 2) and all(v.denominator==1 for v in args),'rational constructor arity')
    return args[0] if k=='Integer' else args[0]/args[1]

def symbol_name(node):
    require(kind(node)=='Symbol' and len(node.args)==1 and isinstance(node.args[0],ast.Constant) and type(node.args[0].value) is str,'native Symbol constructor')
    # Assumptions are preserved in the structural joins, never executed.
    require(all(v.arg is not None and isinstance(v.value,ast.Constant) and type(v.value.value) in (bool,int) for v in node.keywords),'literal symbol assumptions')
    return node.args[0].value

def tuple_args(node):
    require(kind(node)=='Tuple' and not node.keywords,'native Tuple constructor');return node.args

def tagged(node,name):
    matches=[]
    for item in tuple_args(node):
        args=tuple_args(item)
        if len(args)==2 and kind(args[0])=='Str' and len(args[0].args)==1 and isinstance(args[0].args[0],ast.Constant) and args[0].args[0].value==name:matches.append(args[1])
    require(len(matches)==1,'unique native tag '+name);return matches[0]

def case_label(node):
    if kind(node)=='Str':
        require(len(node.args)==1 and not node.keywords and isinstance(node.args[0],ast.Constant) and type(node.args[0].value) is str,'literal native case string')
        return node.args[0].value
    value=rational(node)
    require(value.denominator==1,'integer native case label')
    return int(value)

def select_case(node,expected):
    require(type(expected) is list and all(type(v) in (str,int) for v in expected),'literal expected case labels')
    census=[];matches=[]
    for index,case in enumerate(tuple_args(node)):
        parts=tuple_args(case);require(len(parts)==2,'native labeled case pair')
        labels=[case_label(v) for v in tuple_args(parts[0])]
        census.append(labels)
        if labels==expected:matches.append((index,parts[1]))
    require(len(matches)==1,'unique original labeled case '+repr(expected))
    return matches[0][0],matches[0][1],census

def literal_dict_item(node,key):
    require(isinstance(node,ast.Dict),'literal source dictionary')
    matches=[v for k,v in zip(node.keys,node.values) if isinstance(k,ast.Constant) and type(k.value) is str and k.value==key]
    require(len(matches)==1,'unique original dictionary key '+key)
    return matches[0]

def export_restore_literal(module,key):
    assignments=[n for n in module.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='_LEDGER' for t in n.targets)]
    require(len(assignments)==1,'unique original _LEDGER assignment')
    entry=literal_dict_item(assignments[0].value,key);value=literal_dict_item(entry,'value')
    require(kind(value)=='_restore' and len(value.args)==1 and not value.keywords and isinstance(value.args[0],ast.Constant) and type(value.args[0].value) is str,'original keyed restore literal')
    return value.args[0].value,value

def unit_tuple(value):
    require(type(value) in (list,tuple) and len(value)==3,'three unit exponents')
    return tuple(F(v) if type(v) in (int,str,F) else rational(parse(v['srepr'])) for v in value)

def add(*units):return tuple(sum(v[i] for v in units) for i in range(3))

def scale(unit,value):return tuple(v*value for v in unit)

class Units:
    def __init__(self,registry):
        self.registry={k:unit_tuple(v) for k,v in registry.items()};self.events=[]
    def dimension(self,node,address='root'):
        k=kind(node)
        if k in ('Integer','Rational'):
            value=rational(node);result=None if value==0 else ZERO
        elif isinstance(node,ast.Name) and node.id in ('I','pi','E'):result=ZERO
        elif k=='Symbol':
            name=symbol_name(node);require(name in self.registry,'unknown native unit '+name);result=self.registry[name]
        elif k in ('Add','Mul'):
            require(bool(node.args) and not node.keywords,'native arithmetic arity')
            ds=[self.dimension(v,address+'/'+str(i)) for i,v in enumerate(node.args)]
            if k=='Add':
                live=[v for v in ds if v is not None]
                require(not live or all(v==live[0] for v in live),'inhomogeneous native Add at '+address)
                result=live[0] if live else None
            else:result=None if any(v is None for v in ds) else add(*ds)
        elif k=='Pow':
            require(len(node.args)==2 and not node.keywords,'native power arity');exponent=rational(node.args[1]);base=self.dimension(node.args[0],address+'/base')
            require(base is not None or exponent>0,'zero denominator/undefined unit');result=None if base is None else scale(base,exponent)
        elif k in ('exp','sin','cos','sinh','cosh','tanh'):
            require(len(node.args)==1 and not node.keywords,'native function arity');d=self.dimension(node.args[0],address+'/argument');require(d in (None,ZERO),'dimensionful transcendental argument');result=ZERO
        elif k=='DiracDelta':
            require(len(node.args)==1 and not node.keywords,'undifferentiated native delta');d=self.dimension(node.args[0],address+'/argument');require(d is not None,'delta of literal zero is unresolved');result=scale(d,-1)
        else:raise ValueError('unsupported native unit constructor '+str(k or type(node).__name__))
        self.events.append({'address':address,'constructor':ast.unparse(node),'unit':None if result is None else [str(v) for v in result],'literalZeroUnitUnspecified':result is None})
        return result


def quoted_restore_line(text):
    outer=ast.parse('{'+text.strip()+'}',mode='eval').body
    require(isinstance(outer,ast.Dict) and len(outer.keys)==1 and isinstance(outer.keys[0],ast.Constant) and outer.keys[0].value=='value','native value-line wrapper')
    value=outer.values[0]
    require(kind(value)=='_restore' and len(value.args)==1 and not value.keywords and isinstance(value.args[0],ast.Constant) and type(value.args[0].value) is str,'quoted native restore argument')
    return value.args[0].value
