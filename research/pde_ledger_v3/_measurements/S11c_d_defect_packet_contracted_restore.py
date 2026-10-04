"""Task-local constructor transport for the reviewed saved whole operands.

Planning is stdlib syntax inspection only. Materialization is called only with
the scientific context supplied by the guarded worker. No eval or sympify.
"""
import ast

FUNCTION_ARITIES={'common_outgoing_q':1,'reference_height_hat':1,
                  'reference_slope_hat':1,'Hwhole':3,
                  'Jwhole_plus':8,'Jwhole_minus':8,
                  'Dwhole_plus':8,'Dwhole_minus':8}
ASSUMPTIONS={'real','positive','nonnegative','nonzero','finite','commutative'}


def require(value,message):
    if value is not True:raise ValueError(message)


def constructor_plan(record):
    require(type(record) is dict and type(record.get('srepr')) is str,'full saved constructor operand')
    def integer(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return n.value
        require(isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub) and isinstance(n.operand,ast.Constant) and type(n.operand.value) is int,'integer literal only')
        return -n.operand.value
    def string(n):
        require(isinstance(n,ast.Constant) and type(n.value) is str,'literal constructor name')
        return n.value
    def walk(n):
        if isinstance(n,ast.Name):
            require(n.id in ('I','pi'),'explicit scalar constant only')
            return {'kind':'constant','name':n.id}
        require(isinstance(n,ast.Call),'constructor call only')
        if isinstance(n.func,ast.Call):
            f=n.func
            require(isinstance(f.func,ast.Name) and f.func.id=='Function' and len(f.args)==1 and not f.keywords and not n.keywords,'named applied function only')
            name=string(f.args[0]);require(name in FUNCTION_ARITIES and len(n.args)==FUNCTION_ARITIES[name],'allowed original function and exact arity')
            return {'kind':'function','name':name,'args':[walk(v) for v in n.args]}
        require(isinstance(n.func,ast.Name),'no attributes or dynamic call targets')
        name=n.func.id;require(name in ('Integer','Rational','Symbol','Add','Mul','Pow','sinh'),'explicit original scalar constructor')
        if name=='Symbol':
            require(len(n.args)==1,'symbol arity');symbol=string(n.args[0])
            require(all(k.arg in ASSUMPTIONS and isinstance(k.value,ast.Constant) and type(k.value.value) is bool for k in n.keywords),'original boolean assumptions only')
            require(len({k.arg for k in n.keywords})==len(n.keywords),'unique symbol assumptions')
            return {'kind':'Symbol','name':symbol,'assumptions':{k.arg:k.value.value for k in n.keywords}}
        require(not n.keywords,'no constructor keyword execution')
        if name in ('Integer','Rational'):
            require(len(n.args)==(1 if name=='Integer' else 2),'rational/integer arity')
            args=[integer(v) for v in n.args]
            require(name!='Rational' or args[1]!=0,'nonzero rational denominator')
        else:
            require(len(n.args)>0 if name in ('Add','Mul') else len(n.args)==(2 if name=='Pow' else 1),'exact scalar constructor arity')
            args=[walk(v) for v in n.args]
        return {'kind':name,'args':args}
    return walk(ast.parse(record['srepr'],mode='eval').body)


def restore_scalar(sp,record):
    plan=constructor_plan(record)
    constructors={'Integer':sp.Integer,'Rational':sp.Rational,'Add':sp.Add,
                  'Mul':sp.Mul,'Pow':sp.Pow,'sinh':sp.sinh}
    def materialize(node):
        kind=node['kind']
        if kind=='constant':return sp.I if node['name']=='I' else sp.pi
        if kind=='Symbol':return sp.Symbol(node['name'],**node['assumptions'])
        if kind=='function':return sp.Function(node['name'])(*[materialize(a) for a in node['args']])
        args=node['args'] if kind in ('Integer','Rational') else [materialize(a) for a in node['args']]
        return constructors[kind](*args)
    return materialize(plan)
