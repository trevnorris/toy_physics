"""Translate the two accepted S15 print grammars and compare on their recorded domain."""
from collections import Counter
from functools import lru_cache
from itertools import product
import math
import re

import numpy as np
from scipy.special import roots_legendre, roots_jacobi, spence
import sympy as s
from sympy.core.function import AppliedUndef
from sympy.core.relational import Relational

import S15_cross_engine_comparator as original

Unsupported = original.Unsupported
NAMES = original.Names(original.NAME_TABLE)
# Reading-engine definitions: Python lines 357–415, Wolfram lines 226–286.
for row in original.NAME_TABLE:
    if row.py.startswith(('rho_', 'u_', 'W_', 'Pi_o_', 'f_')) or row.py in ('varrho','eta','zeta'):
        key = row.py + '_change'
        NAMES.maps['py'][key] = key
        NAMES.maps['wl'][f'change[{row.wl}]'] = key
        NAMES.arities[key] = row.arity
for j in range(1, 4):
    key = f'X_1_{j}_velocity_change'
    NAMES.maps['py'][key] = key
    NAMES.maps['wl'][f'centreVelocityChange[{j}]'] = key
    NAMES.arities[key] = 1
for j in (1, 2):
    for head in ('A', 'B'):
        key = f'{head}_{j}'
        NAMES.maps['py'][key] = key
        NAMES.maps['wl'][f'{head*2}[{j}]'] = key
        NAMES.arities[key] = 1

for i, key in enumerate(sorted(NAMES.arities)):
    if NAMES.arities[key] and key not in original.CODES:
        original.CODES[key] = 200+i

EXP = s.symbols('a_1 a_2 b_1 b_2')
DIRECTION = s.symbols('N_1:4')
CHART = {f'directionSeparation[{j}]':f'N_{j}' for j in range(1, 4)}
CHART.update({f'{k*2}[{j}]':f'{k}_{j}' for k in ('a','b') for j in (1,2)})
CHART['dist'] = 'd'
for key, value in CHART.items():
    NAMES.maps['wl'][key] = value
    NAMES.maps['py'][value] = value

def token_pattern(key):
    return re.escape(key).replace(r'\,', ',').replace(',', r',\s*')

# Only source-declared identifiers are translated. Unknown heads remain unsupported.
WL_SCALARS = re.compile(r'(?<![\w$])('+'|'.join(token_pattern(k) for k in CHART)+r')(?![\w$])')
EXTRA = {'log':s.log, 'Log':s.log, 'PolyLog':s.polylog, 'polylog':s.polylog,
         'AppellF1':s.appellf1, 'appellf1':s.appellf1}

def parse(text, engine):
    if engine == 'wl':
        text = WL_SCALARS.sub(lambda m:CHART[m[0].replace(' ','')], text)
        node = original.wl_parse(text)
    else:
        node = original.py_parse(text)
    return translate(node, engine, ())


@lru_cache(maxsize=32768)
def translate(node, engine, scope):
    h, args = node.head, node.args
    def cv(n):
        return translate(n, engine, scope)
    if h == 'Symbol':
        name = node.value
        if name in dict(scope):
            return dict(scope)[name]
        if name in ('Pi','pi'):
            return s.pi
        if name in ('True','False'):
            return name == 'True'
        if name in CHART.values():
            return s.Symbol(name)
        if name not in NAMES.maps[engine]:
            # SymPy's lexical dummies are bound by Subs or Derivative.
            if engine == 'py' and name.startswith('_'):
                return s.Symbol(name)
            raise Unsupported('unknown free scalar '+name)
        return s.Symbol(NAMES.name(name,engine))
    if h == 'Number':
        return s.Rational(node.value)
    if h == 'Boolean':
        return s.true if node.value.lower()=='true' else s.false
    if h == 'Sequence':
        return s.Tuple(*(cv(q) for q in args))
    if h == 'Text':
        return node.value
    if h == 'Apply':
        head = args[0]
        if head.head=='Apply' and original.call(head.args[0], 'Derivative'):
            name = original.key_text(head.args[1])
            key = NAMES.name(name, engine)
            points = tuple(cv(q) for q in args[1:])
            orders = tuple(int(cv(q)) for q in head.args[0].args[1:])
            if len(points)!=len(orders) or NAMES.arities.get(key)!=len(points):
                raise Unsupported('unknown derivative slot or head '+name)
            return original.partial_derivative(s.Function(key), orders, points)
        if head.head=='Prime':
            return original.math_node(node,NAMES,engine,dict(scope))
        name, operands = original.application(node)
        if name is None:
            name = original.key_text(head)
        args = operands
        key = NAMES.maps[engine].get(name)
        if key in NAMES.arities and NAMES.arities[key] > 0:
            if len(args)!=NAMES.arities[key]:
                raise Unsupported('field arity '+name)
            return s.Function(key)(*(cv(q) for q in args))
        if name in EXTRA:
            return EXTRA[name](*(cv(q) for q in args), evaluate=False)
        h = original.BUILTINS.get(name,name)
    if h == 'Subs':
        variables,points = [(q.args if q.head=='Sequence' else (q,)) for q in args[1:]]
        if len(variables)!=len(points) or any(v.head!='Symbol' for v in variables):
            raise Unsupported('Subs binders')
        local = dict(scope)
        slots = tuple(s.Symbol('_reading_slot'+str(len(scope)+i)) for i in range(len(variables)))
        local.update((v.value,slot) for v,slot in zip(variables,slots))
        obj = translate(args[0],engine,tuple(local.items()))
        return original.normalize_partial(s.Subs(obj,slots,tuple(cv(q) for q in points)))
    if h == 'Derivative':
        orders = [(cv(q.args[0]),int(cv(q.args[1]))) if q.head=='Sequence' else cv(q) for q in args[1:]]
        return original.normalize_partial(s.Derivative(cv(args[0]),*orders,evaluate=False))
    if h == 'Integral':
        bounds = [tuple(cv(q)) for q in args[1:]]
        if any(len(spec)!=3 for spec in bounds):
            raise Unsupported('nondefinite integral')
        # Both sources print the rectangular spherical chart. When Wolfram iterators
        # have dependent bounds, its outer-first convention needs reversal.
        if engine=='wl':
            bounds.reverse()
        if all(not (bound.free_symbols & {r[0] for r in bounds})
               for _,lo,hi in bounds for bound in (lo,hi)):
            order = {'s':0,'theta':1,'phi':2}
            bounds.sort(key=lambda r:order[str(r[0])])
        integrand = cv(args[0])
        return s.Integral(integrand,*bounds) if integrand!=0 else s.S.Zero
    if h in ('Add','Mul','Pow','sin','cos','sqrt','Abs','Rational'):
        return getattr(s,h)(*(cv(q) for q in args))
    if h in original.RELATIONS:
        return original.native_relation(h,*(cv(q) for q in args))
    if h in ('And','Or','Not'):
        return getattr(s,h)(*(cv(q) for q in args),evaluate=False)
    if h=='Piecewise':
        pairs = args[0].args if engine=='wl' else args
        rows=[(cv(q.args[0]),cv(q.args[1])) for q in pairs]
        if engine=='wl':
            rows.append((cv(args[1]) if len(args)>1 else s.S.Zero,True))
        return s.Piecewise(*rows,evaluate=False)
    raise Unsupported('unhandled '+engine+' head '+h)


@lru_cache(maxsize=8192)
def outer_polynomial(expr):
    """Expand only ring operations. Free functions and inverse powers are opaque."""
    atoms = {}
    def lift(v):
        if v.is_Add or v.is_Mul:
            return v.func(*(lift(q) for q in v.args))
        if v.is_Pow and v.exp.is_Integer and 0 < v.exp <= 6:
            return lift(v.base)**v.exp
        if v.is_Number or v.is_Symbol:
            return v
        if v not in atoms:
            atoms[v] = s.Dummy()
        return atoms[v]
    lifted = lift(expr)
    return s.expand(lifted).xreplace({token:value for value,token in atoms.items()})


@lru_cache(maxsize=16384)
def normal(expr):
    """Exact algebra in the printed rectangular domain. Only subexpressions used."""
    if not expr.args:
        return expr
    args = tuple(normal(q) for q in expr.args)
    if expr.func in (s.Subs,s.Derivative):
        result = expr if args==expr.args else expr.func(*args)
        return original.normalize_partial(result)
    if isinstance(expr, s.Integral):
        return s.S.Zero if args[0]==0 else s.Integral(*args)
    if expr.is_Pow and args[1].is_Rational and (args[1]<0 or not args[1].is_Integer) and args[0].is_Add:
        # A squared distance may be printed as a sum of squares or expanded.
        # Expand its base without distributing the fractional power.
        return s.Pow(outer_polynomial(args[0]),args[1])
    if isinstance(expr,Relational):
        return expr.func(*args,evaluate=False)
    if expr.func in (s.appellf1,s.polylog,s.log,s.Piecewise):
        return expr.func(*args,evaluate=False)
    return expr if args==expr.args else expr.func(*args)

def exact_residual(a,b):
    if a==b:
        return s.S.Zero
    # Keep native piecewise guards intact. Test each explicit exponent locus
    # separately when algebraic identity is not established in this form.
    if a.has(s.Piecewise,s.appellf1,s.polylog) or b.has(s.Piecewise,s.appellf1,s.polylog):
        return a-b
    a,b=normal(a),normal(b)
    if a==b:
        return s.S.Zero
    # Preserve inner integrals when a Leibniz rate is nested in an angular
    # integral; linearity is valid for the finite, regular domain in the list.
    residual=s.expand(a-b,deep=False)
    groups={}
    rest=[]
    for term in s.Add.make_args(residual):
        factors=s.Mul.make_args(term)
        terms=[v for v in factors if isinstance(v,s.Integral)]
        if len(terms)!=1:
            rest.append(term)
            continue
        integral=terms[0]
        scale=s.Mul(*(v for v in factors if v is not integral))
        if scale.free_symbols.intersection(integral.variables):
            rest.append(term)
            continue
        groups.setdefault(integral.limits,[]).append(scale*integral.function)
    for limits,parts in groups.items():
        v=outer_polynomial(s.Add(*parts))
        if v!=0:
            rest.append(s.Integral(v,*limits))
    return s.Add(*rest)

def branch_loci(text,engine):
    if engine=='py':
        import ast
        clauses=ast.literal_eval(text)
        values=[parse(clause,engine) for clause in clauses]
    else:
        values=list(parse(text,engine))
    loci=set()
    for value in values:
        if value in (True,False):
            continue
        if not isinstance(value,(s.Equality,s.Unequality)):
            raise Unsupported('non-equality exponent branch '+str(value))
        poly=s.Poly(value.lhs-value.rhs,*EXP)
        if poly.total_degree()!=1:
            raise Unsupported('nonlinear exponent branch '+str(value))
        # Store the boundary hyperplane; native polarity and branch order
        # remain in the original operands and in each Piecewise evaluation.
        loci.add(tuple(poly.coeff_monomial(q)/poly.LC() for q in (*EXP,1)))
    return tuple(sorted(loci))

# Numerical fallback. A tensor quadrature retains each native chart, bound, and
# integrand. Samples carry the point, working precision and refinement delta.
FORMS = original.FORMS[:2]
STRATA = original.STRATA
POINTS = original.POINTS
RULES = (10, 16, 24)
GENERIC_EXP = ((1,2,3,4),(2,3,1,2),(3,1,2,3))

def exponent_samples(point,loci):
    base=tuple(map(s.Integer,GENERIC_EXP[point]))
    # Two distinct nonintegral offsets keep the baseline off common exceptional loci.
    base=tuple(v+s.Rational(j+1,8) for j,v in enumerate(base))
    yield 'generic',base
    for index,locus in enumerate(loci):
        slot=next(i for i,v in enumerate(locus[:-1]) if v)
        values=list(base)
        values[slot]=-(locus[-1]+sum(a*x for j,(a,x) in enumerate(zip(locus,base)) if j!=slot))/locus[slot]
        yield f'locus_{index}',tuple(values)


@lru_cache(maxsize=8192)
def jet_function(name,arity,orders,form,stratum,point):
    if name.endswith('_velocity_change'):
        # A variation of instantaneous velocity is a free time function, also on a stationary stratum.
        sample_head = 'response_' + name
        original.CODES[sample_head] = original.CODES[name]
    else:
        sample_head = name
    jet=original.profile_jet(sample_head,arity,orders,form,stratum,point,())
    return s.lambdify(original.FORMAL[:arity],jet,modules='numpy',cse=True,docstring_limit=0)


@lru_cache(maxsize=32768)
def euler_appell(a,b,c,d,x,y,order):
    # Euler's normalized beta integral, valid for these real parameters and bases.
    if not d>a>0 or max(x,y)>=1:
        raise Unsupported('Appell F1 outside the real positive-base Euler branch')
    nodes,weights=roots_jacobi(order,d-a-1,a-1)
    z=(nodes+1)/2
    return float(np.dot(weights,(1-x*z)**(-b)*(1-y*z)**(-c))/weights.sum())


@lru_cache(maxsize=32)
def quadrature(order):
    return roots_legendre(order)


class Sample:
    def __init__(self,form,stratum,point,exponents,order):
        self.form,self.stratum,self.point,self.order=form,stratum,point,order
        self.t=float(POINTS[point][-1])
        self.constants=dict(zip(('x_1','x_2','x_3','t'),map(float,POINTS[point])))
        self.constants.update({str(k):float(v) for k,v in zip(EXP,exponents)})
        def centre(k,j):
            return jet_function(f'X_{k}_{j}',1,(0,),form,stratum,point)(self.t)
        vector=[centre(1,j)-centre(2,j) for j in (1,2,3)]
        distance=math.sqrt(sum(q*q for q in vector))
        self.constants.update({str(k):v/distance for k,v in zip(DIRECTION,vector)})
        self.constants['d']=distance
        radii=[jet_function(k,1,(0,),form,stratum,point)(self.t) for k in ('r_in','R')]
        if not distance>radii[1]>radii[0]>0:
            raise Unsupported('sample outside shell domain')
        self.envs=[{}]
        self.grids={}
        self.cache={}

    def value(self,expr,environment=0):
        key=expr,environment
        if key not in self.cache:
            self.cache[key]=self.compute(expr,environment)
        return self.cache[key]

    def compute(self,v,environment):
        ev=lambda q:self.value(q,environment)
        if v is s.true or v is s.false:
            return v is s.true
        if v.is_Number or v==s.pi:
            return float(v)
        if v.is_Symbol:
            if v in self.envs[environment]:
                return self.envs[environment][v]
            if str(v) in self.constants:
                return self.constants[str(v)]
            raise Unsupported('unbound sample coordinate '+str(v))
        if v.is_Add:
            return sum(ev(q) for q in v.args)
        if v.is_Mul:
            return math.prod(ev(q) for q in v.args)
        if v.is_Pow:
            return np.power(ev(v.base),ev(v.exp))
        if isinstance(v,s.Integral):
            return self.integral(v,environment)
        if isinstance(v,AppliedUndef):
            sig=(v.func,(0,)*len(v.args),v.args)
        elif isinstance(v,(s.Derivative,s.Subs)):
            sig=original.partial_signature(v)
            if sig is None:
                if isinstance(v,s.Subs):
                    expanded=v.doit()
                    if expanded!=v:
                        return ev(expanded)
                raise Unsupported('unevaluated composite derivative or substitution')
        else:
            sig=None
        if sig is not None:
            head,orders,points=sig
            fun=jet_function(head.__name__,len(points),orders,self.form,self.stratum,self.point)
            return fun(*(ev(q) for q in points))
        if v.func in (s.sin,s.cos,s.log):
            return {s.sin:np.sin,s.cos:np.cos,s.log:np.log}[v.func](ev(v.args[0]))
        if v.func==s.appellf1:
            return euler_appell(*(ev(q) for q in v.args),self.order)
        if v.func==s.polylog:
            if v.args[0]!=2: raise Unsupported('polylog order outside Li2')
            return spence(1-ev(v.args[1]))
        if isinstance(v,s.Piecewise):
            for value,condition in v.args:
                if ev(condition): return ev(value)
            raise Unsupported('no selected piecewise branch')
        if isinstance(v,Relational):
            left,right=map(ev,v.args)
            return {'==':lambda:left==right,'!=':lambda:left!=right,'>':lambda:left>right,
                    '<':lambda:left<right,'>=':lambda:left>=right,'<=':lambda:left<=right}[v.rel_op]()
        if v.func in (s.And,s.Or):
            return (all if v.func==s.And else any)(ev(q) for q in v.args)
        raise Unsupported('numerical head '+v.func.__name__)

    def integral(self,v,environment):
        key=environment,v.limits
        if key not in self.grids:
            env=self.envs[environment].copy()
            nodes,weights=quadrature(self.order)
            product_weights=1.
            for symbol,low,high in v.limits:
                if symbol in env:
                    raise Unsupported('shadowed integral variable')
                lo,hi=self.value(low,environment),self.value(high,environment)
                axis={'s':0,'theta':1,'phi':2}[str(symbol)]
                shape=[1,1,1]
                shape[axis]=-1
                env[symbol]=((hi-lo)*nodes/2+(hi+lo)/2).reshape(shape)
                product_weights=product_weights*(weights*(hi-lo)/2).reshape(shape)
            child=len(self.envs)
            self.envs.append(env)
            self.grids[key]=child,product_weights,tuple({'s':0,'theta':1,'phi':2}[str(q[0])] for q in v.limits)
        child,weights,axes=self.grids[key]
        result=np.sum(weights*self.value(v.function,child),axis=axes,keepdims=True)
        return float(result.item()) if result.size==1 else result


def sample_values(a,b,loci):
    rows=[]
    failures=Counter()
    for form,stratum,point in product(FORMS,STRATA,range(len(POINTS))):
        for branch,exponents in exponent_samples(point,loci):
            key=[form,stratum,point,branch]
            try:
                values=[]
                with np.errstate(over='raise',invalid='raise',divide='raise'):
                    for order in RULES:
                        ctx=Sample(form,stratum,point,exponents,order)
                        pair=[float(ctx.value(expr)) for expr in (a,b)]
                        if not all(map(math.isfinite,pair)):
                            raise Unsupported('nonfinite operand')
                        values.append(pair)
                py,wl=values[-1]
                # Empirical error (not an enclosure): two refinements per operand
                # and a floating-point allowance. Never promote it to an identity.
                delta=max(abs(values[j][k]-values[j-1][k]) for j in (1,2) for k in (0,1))
                scale=max(abs(py),abs(wl),1.)
                tolerance=8*delta+1.e-9*scale
                last_change=max(abs(values[2][k]-values[1][k]) for k in (0,1))
                if abs(py-wl)>tolerance:
                    verdict='separated'
                elif last_change > 1.e-6*scale:
                    verdict='poorly_resolved'
                else:
                    verdict='within_estimate'
                rows.append(dict(case=key,exponents=list(map(str,exponents)),values=values,
                                 residual=py-wl,refinement_delta=delta,tolerance=tolerance,status=verdict))
            except (Unsupported,ValueError,ZeroDivisionError,OverflowError,FloatingPointError,TypeError) as exc:
                failures[str(exc)]+=1
                rows.append(dict(case=key,exponents=list(map(str,exponents)),status='unresolved',reason=str(exc)))
    counts=Counter(q['status'] for q in rows)
    if counts['separated']:
        status='sample_difference'
    elif failures or counts['poorly_resolved']:
        status='sample_coverage_incomplete'
    else:
        status='sample_agreement'
    return {'status':status,'sample_counts':dict(counts),'samples':rows}
