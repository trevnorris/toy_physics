"""Strict arithmetic/source helpers for a NEW guarded contraction certificate.

No scientific imports or native execution at import time. A guarded caller
supplies SymPy. Standalone tooling tests use manufactured stdlib operands only.
No old scientific function, numerical integral or rule constructor is called.
"""
import ast
from fractions import Fraction as F
import hashlib
from pathlib import Path

ZERO=(F(0),F(0),F(0))
MOMENTUM=(F(-1),F(0),F(0))
def require(v,message):
    if v is not True:raise ValueError(message)

def unit(v):
    require(type(v) in (list,tuple) and len(v)==3,'unit triple')
    require(all(type(x) in (str,int,F) for x in v),'rational unit exponents only')
    return tuple(F(x) for x in v)
def add(*values):return tuple(sum((unit(v)[i] for v in values),F(0)) for i in range(3))
def scale(v,n):
    require(type(n) in (int,str,F),'rational dimension power')
    return tuple(x*F(n) for x in unit(v))

def fragment(record):
    text=Path(record['source']).read_text()
    require(hashlib.sha256(text.encode()).hexdigest()==record['sourceSha256'],'original source hash')
    nodes=[n for n in ast.walk(ast.parse(text)) if getattr(n,'lineno',None)==record['line'] and getattr(n,'end_lineno',None)==record['endLine'] and ast.get_source_segment(text,n)==record['text']]
    require(len(nodes)==1,'unique exact original source fragment')
    require(hashlib.sha256(record['text'].encode()).hexdigest()==record['fragmentSha256'],'fragment hash')
    return nodes[0]

class Arithmetic:
    """Interpret only the actual polynomial/rational assignment AST, not a call.

    Every intermediate and complete return is retained by the caller. No eval,
    imports, attributes, scientific functions or Python function execution.
    """
    def __init__(self,env,integer):self.env=dict(env);self.integer=integer;self.events=[]
    def expression(self,n):
        if isinstance(n,ast.Name):
            require(n.id in self.env,'bound source operand '+n.id);return self.env[n.id]
        if isinstance(n,ast.Constant):
            require(type(n.value) is int,'integer arithmetic literals only');return self.integer(n.value)
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -self.expression(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=self.expression(n.left),self.expression(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
            if isinstance(n.op,ast.Pow):
                require(isinstance(n.right,ast.Constant) and type(n.right.value) is int and 0<=n.right.value<=4,'source power scope');return a**n.right.value
        raise ValueError('unsupported original arithmetic syntax '+ast.dump(n))
    def statements(self,function):
        require(isinstance(function,ast.FunctionDef),'original function definition')
        require(not function.decorator_list and not function.args.defaults and not function.args.kwonlyargs and function.args.vararg is None and function.args.kwarg is None,'simple original signature')
        require({a.arg for a in function.args.args}==set(self.env),'exact original formal arguments')
        answer=None
        for n in function.body:
            if isinstance(n,ast.Expr) and isinstance(n.value,ast.Constant) and type(n.value.value) is str:continue
            require(answer is None,'nothing after original return')
            if isinstance(n,ast.Assign):
                require(len(n.targets)==1 and isinstance(n.targets[0],ast.Name),'simple assignment')
                value=self.expression(n.value);name=n.targets[0].id;self.env[name]=value
                self.events.append({'assignment':name,'source':ast.unparse(n.value),'value':value})
            elif isinstance(n,ast.Return):
                require(isinstance(n.value,(ast.List,ast.Tuple)),'complete original vector return')
                answer=[self.expression(v) for v in n.value.elts]
                self.events.append({'return':answer})
            else:raise ValueError('unsupported original statement '+ast.dump(n))
        require(answer is not None,'original return present');return answer

def restore_scalar(sp,record):
    """Restore the limited actual numeric-template constructor tree, no eval."""
    require(type(record) is dict and type(record.get('srepr')) is str,'full constructor operand')
    def walk(n):
        if isinstance(n,ast.Constant) and type(n.value) in (str,int,bool):return n.value
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):
            require(isinstance(n.operand,ast.Constant) and type(n.operand.value) is int,'signed integer literal');return -n.operand.value
        if isinstance(n,ast.Name) and n.id=='I':return sp.I
        require(isinstance(n,ast.Call) and isinstance(n.func,ast.Name),'constructor call only')
        name=n.func.id;require(name in ('Integer','Rational','Symbol','Add','Mul','Pow'),'allowed scalar constructor '+name)
        args=[walk(v) for v in n.args]
        if name=='Symbol':
            require(len(args)==1 and type(args[0]) is str,'original symbol name')
            require(all(v.arg in ('real','positive','nonnegative','nonzero','finite','commutative') and isinstance(v.value,ast.Constant) and type(v.value.value) is bool for v in n.keywords),'original boolean assumptions only')
            require(len({v.arg for v in n.keywords})==len(n.keywords),'unique assumptions')
            return sp.Symbol(args[0],**{v.arg:v.value.value for v in n.keywords})
        require(not n.keywords,'no constructor keyword execution')
        if name=='Integer':require(len(args)==1 and type(args[0]) is int,'integer arity/type')
        if name=='Rational':require(len(args)==2 and all(type(v) is int for v in args) and args[1]!=0,'rational arity/domain')
        if name=='Pow':require(len(args)==2,'power arity')
        if name in ('Add','Mul'):require(len(args)>0,'nonempty nary constructor')
        return getattr(sp,name)(*args)
    return walk(ast.parse(record['srepr'],mode='eval').body)

def named_symbols(expr):
    result={}
    for s in expr.free_symbols:
        require(str(s) not in result,'unique symbol name/assumptions');result[str(s)]=s
    return result

def dimension(sp,expr,env,events):
    if expr in env:v=unit(env[expr])
    elif expr==sp.I or expr.is_Rational:v=ZERO
    elif expr.is_Add:
        terms=[dimension(sp,x,env,events) for x in expr.args]
        require(all(t==terms[0] for t in terms),'new expression homogeneous addition');v=terms[0]
    elif expr.is_Mul:v=add(*(dimension(sp,x,env,events) for x in expr.args))
    elif expr.is_Pow:
        a,b=expr.args;require(b.is_Rational is True,'rational dimensional power')
        v=scale(dimension(sp,a,env,events),F(int(b.p),int(b.q)))
    else:raise ValueError('unbound new dimensional operand '+str(expr))
    events.append({'operand':expr,'unit':v});return v

def finite_window(K,T,m):
    require(all(type(v) in (int,str,F) for v in (K,T,m)),'exact rational domain operands')
    K,T,m=F(K),F(T),F(m);require(T>K>0,'declared finite window ordering')
    lo,hi=max(-K,m-T),min(K,m+T)
    return {'lower':lo,'upper':hi,'empty':lo>hi,'degenerate':lo==hi,'central':abs(m)<=T-K,'insideM':abs(m)<=T+K}

def rectangle_contains(K,T,k,l,t):
    K,T,k,l,t=map(F,(K,T,k,l,t));return -K<=k<=K and -K<=l<=K and -T<=t<=T

def mapped_contains(K,T,k,l,m,window_variable):
    require(window_variable in ('k','l'),'explicit actual clipped variable')
    K,T,k,l,m=map(F,(K,T,k,l,m));w=finite_window(K,T,m)
    clipped=k if window_variable=='k' else l;full=l if window_variable=='k' else k
    return w['insideM'] and not w['empty'] and w['lower']<=clipped<=w['upper'] and -K<=full<=K

def quadrant_polynomial_certificate(sp,expression,variables):
    poly=sp.Poly(expression,*variables)
    require(all(v.is_Rational is True and v>=0 for v in poly.coeffs()),'literal nonnegative polynomial coefficients')
    require(all(v.is_nonnegative is True for v in variables),'declared nonnegative real operands')
    return {'expression':expression,'variables':variables,'terms':poly.terms(),'allCoefficientsNonnegative':True,'domain':'all declared variables real and nonnegative','notMeasureTheoryProof':True}
