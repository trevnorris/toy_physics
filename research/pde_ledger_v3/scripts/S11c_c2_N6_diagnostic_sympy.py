#!/usr/bin/env python3
"""S11c-c2 §5c; construction authority: S11c_c2_N6_route2_spec_astra.md.

1. The script may PRINT computed objects. It may NOT state conclusions. An emit
payload must be a CAS object (an expression, a solved root, a boolean from a
symbolic/modular test), ⛔ never prose describing a result.
2. PRINT the residual; do NOT assert it. `assert residual == 0` is the builder
writing down the expected output. Compute → emit → (never) assert. Here the
residual is a PIT numeric table; residual-zero is ⛔ never an exit.
3. Interpretation belongs to the STEP RECORD. The script does not editorialise.

SUPPLIED: §5c object; §3a close; §3c increment, slot-linearity and extract;
c1 imported face response; S11c-a/b material builders. These are premises.
WITHHELD: any expected value / acceptance criterion for R_N6 or the controls.

Arithmetic DAGs retain rational denominators; formal integrals are never
integrated. Their kernels use alpha-aligned momenta and independent field jets.
The constitutive variation is termwise, before coefficient specialization.
No build_case or build_operator invocation occurs in this module.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from functools import lru_cache
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import resource
import sys
import time

os.environ['S11CB_PROJECTION_WORKERS'] = '1'
import sympy as sp
from sympy.core.function import AppliedUndef
from sympy.core.symbol import Str
import S11c_c2_selfenergy_fold_sympy_audit as c

ROOT = Path(__file__).resolve().parent.parent
GRADES = tuple(itertools.product((0, 1), repeat=2))
ROWS = ('U0', 'U1', 'U2', 'THETA', 'E_W')
BLOCKS = tuple(itertools.product(('TRANSVERSE_TO_THICKNESS', 'THICKNESS_TO_TRANSVERSE'),
                                 ('THETA', 'E_W', 'DIV_U')))
ROW_DIM = {'U0': (-2,-2,1), 'U1': (-2,-2,1), 'U2': (-2,-2,1),
           'THETA': (-3,-1,1), 'E_W': (-1,-2,1)}
WEAK_DIM = {block: ((-3,-1,1) if block == BLOCKS[0] else (-1,-2,1)) for block in BLOCKS}
START = time.monotonic()
MAX_OPS = 180000
CASE = {'anchoring': None, 'density': None}
TEMPLATE_TEST_SPATIAL_DERIVATIVES = []


def encode(x):
    if isinstance(x, dict):
        return {str(k): encode(v) for k,v in x.items()}
    if isinstance(x, (tuple, list, sp.Tuple, sp.MatrixBase)):
        return [encode(v) for v in x]
    if isinstance(x, sp.Integer):
        return int(x)
    if isinstance(x, sp.Basic):
        return str(x)
    return x


def emit(name, payload, **labels):
    print(json.dumps(encode({'object': 'S11CC2_' + name, **CASE,
                      'dimension': [0,0,0], 'N12': [[0,0,0]],
                      **labels, 'data': payload}), sort_keys=True), flush=True)


def progress(stage, **counts):
    emit('N6_PROGRESS', {'stage': stage, 'seconds': round(time.monotonic()-START, 3),
                         'rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, **counts})


def bounded(e, stage):
    size = int(sp.count_ops(e))
    if size > MAX_OPS:
        emit('N6_BLOCK_SIZE', {'operations': size, 'limit': MAX_OPS}, stage=stage)
        raise BlockSize(stage)
    return e


class BlockSize(RuntimeError):
    pass


def sha(e):
    return hashlib.sha256(sp.srepr(c.cas(e)).encode()).hexdigest()


def flatten(rows):
    return dict(zip(ROWS, (*rows['U'], rows['THETA'], rows['E_W'])))


def unflatten(rows):
    return {'U': tuple(rows['U'+str(i)] for i in range(3)),
            'THETA': rows['THETA'], 'E_W': rows['E_W']}


# Rational arithmetic DAG. Degree pairs bound numerator and denominator degrees
# without expanding either polynomial. Exclusion bounds use unique inverse nodes.
@dataclass(frozen=True, eq=False)
class Node:
    op: str
    args: tuple
    degree: tuple
    units: frozenset
    wave_degrees: frozenset = frozenset((0,))


@lru_cache(None)
def number(n):
    return Node('number', (int(n),), (0,0), frozenset(((0,0,0),)) if n else frozenset(),
                frozenset((0,)) if n else frozenset())


@lru_cache(None)
def amplitude():
    return Node('number',(1,),(0,0),frozenset(((0,0,0),)),frozenset((1,)))


@lru_cache(None)
def variable(key,units=(0,0,0)):
    return Node('variable', (key,), (1,0), frozenset((tuple(units),)))


@lru_cache(None)
def measure(signature):
    return Node('number',(1,),(0,0),frozenset((((6-signature) if signature else 0,0,0),)))


@lru_cache(None)
def plus(a,b):
    if a == number(0): return b
    if b == number(0): return a
    return Node('add', (a,b), (max(a.degree[0]+b.degree[1], b.degree[0]+a.degree[1]),
                             a.degree[1]+b.degree[1]), a.units | b.units,a.wave_degrees|b.wave_degrees)


@lru_cache(None)
def times(a,b):
    if a == number(0) or b == number(0): return number(0)
    if a == number(1): return b
    if b == number(1): return a
    return Node('mul', (a,b), tuple(x+y for x,y in zip(a.degree,b.degree)),
                frozenset(tuple(x+y for x,y in zip(u,v)) for u in a.units for v in b.units),
                frozenset(x+y for x in a.wave_degrees for y in b.wave_degrees))


@lru_cache(None)
def power(a,n):
    if n == 0: return number(1)
    if n == 1: return a
    d = a.degree if n > 0 else a.degree[::-1]
    return Node('pow', (a,int(n)), tuple(abs(n)*x for x in d),
                frozenset(tuple(n*x for x in u) for u in a.units),frozenset(n*x for x in a.wave_degrees))


def minus(a,b):
    return plus(a,times(number(-1),b))


def nsum(values):
    out = number(0)
    for v in values: out = plus(out,v)
    return out


def gadd(a,b):
    return {g: plus(a.get(g,number(0)),b.get(g,number(0))) for g in set(a)|set(b)}


def gmul(a,b):
    out = {}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            if i+k <= 1 and j+l <= 1:
                g=(i+k,j+l); out[g]=plus(out.get(g,number(0)),times(v,w))
    return out


def gscale(a,n):
    return {g:times(v,n) for g,v in a.items()}


def gsubtract(a,b):
    return gadd(a,gscale(b,number(-1)))


def gpower(a,n):
    base=a.get((0,0),number(0))
    delta={g:v for g,v in a.items() if g!=(0,0)}
    out={}; product={(0,0):number(1)}
    for order in range(3):
        coefficient=times(number(sp.binomial(n,order)),power(base,n-order))
        if coefficient != number(0): out=gadd(out,gscale(product,coefficient))
        product=gmul(product,delta)
    return out


class Compiler:
    def __init__(self, inputs):
        self.inputs=inputs
        self.cache={}
        self.i=variable(('algebraic_i',))
        self.symbol_bindings={}
        self.unknown=set()

    def scalar(self, e, point='X', binding=None):
        e=sp.sympify(e)
        if binding and e in binding: return binding[e]
        key=(e,point)
        if not binding and key in self.cache: return self.cache[key]
        call=lambda x:self.scalar(x,point,binding)
        if e == self.inputs.eps: out=amplitude()
        elif e == sp.I: out=self.i
        elif e.is_Rational:
            out=times(number(e.p),power(number(e.q),-1))
        elif e in self.symbol_bindings: out=self.symbol_bindings[e]
        elif isinstance(e,sp.Add): out=nsum(call(x) for x in e.args)
        elif isinstance(e,sp.Mul):
            out=number(1)
            for x in e.args: out=times(out,call(x))
        elif isinstance(e,sp.Pow) and e.exp.is_Integer: out=power(call(e.base),int(e.exp))
        elif isinstance(e,sp.Symbol):
            # Profile values at X and Y are distinct; mixed-partial identities
            # are retained by the supplied jet names, never by random renaming.
            local=bool(re.match(r'(W_bg|mu_R_bg|[wm]1_profile)(?:$|_d)',e.name))
            match=re.fullmatch(r'([wm]1_profile)_((?:d[123])+)',e.name)
            canonical=(match[1]+'_'+''.join('d'+i for i in sorted(re.findall(r'd([123])',match[2])))) if match else e.name
            units=tuple(c.DIMENSION_SCHEMA.get(canonical,c.dimension(e)))
            out=variable((canonical,point if local else 'global'),units)
        elif isinstance(e,(AppliedUndef,sp.Derivative)):
            out=variable((sp.srepr(e),'formal_jet'),c.dimension(e))
        elif e == sp.pi: out=variable(('pi','global'))
        else:
            raise TypeError(('unsupported arithmetic circuit',e.func,str(e)[:180]))
        if not binding: self.cache[key]=out
        return out

    def grades(self,e,point='X'):
        e=bounded(sp.sympify(e),'coefficient')
        return {g:self.scalar(v,point) for g,v in c.shape_coefficients(
            e.xreplace(self.inputs.profiles),self.inputs.eta,self.inputs.sigma).items() if g in GRADES}

    def routed(self,e,bindings,point='X'):
        """Evaluate a supplied full source circuit with graded source slots.

        Only arithmetic operations are traversed. This never constructs the
        symbolic closed expression and never substitutes after expansion.
        """
        if e in bindings: return bindings[e]
        if not e.has(*bindings): return self.grades(e,point)
        if isinstance(e,sp.Add):
            out={}
            for arg in e.args: out=gadd(out,self.routed(arg,bindings,point))
            return out
        if isinstance(e,sp.Mul):
            out={(0,0):number(1)}
            for arg in e.args: out=gmul(out,self.routed(arg,bindings,point))
            return out
        if isinstance(e,sp.Pow) and e.exp.is_Integer:
            return gpower(self.routed(e.base,bindings,point),int(e.exp))
        raise TypeError(('graded source binding',e.func))


# Face-only raw-case adapter; source mutation is confined to the normal, before
# all raw builders. The parent's cache and virtual kinematics remain untouched.
def face_factory(a,b,inputs,alpha,rho,route,mu,tilt=1):
    original=a.build_face_source
    def source(*args,**kw):
        result=original(*args,**kw)
        if route == 'EULERIAN' and tilt != 1 and args[1] == 1:
            n=result.normal_exact
            conormal=[sp.cancel(v/n[3]) for v in n]
            jet=a.grad_W[0]
            contribution=sp.diff(conormal[0].subs(result.parameter,0),jet)*jet
            conormal[0]=conormal[0]+(tilt-1)*contribution
            norm=sp.sqrt(sum(v*v for v in conormal))
            result=replace(result,normal_exact=tuple(v/norm for v in conormal))
        return result
    a.build_face_source=source
    try:
        sources={s:source(alpha,s,'DELTA_W',rho,route=route) for s in c.FACES}
        traction={s:a.finalize(a.traction_raw(v)) for s,v in sources.items()}
        closure={s:a.finalize(a.closure_raw(v)) for s,v in sources.items()}
        work={d:a.finalize(a.virtual_work_cases(alpha,'DELTA_W',d,rho,route=route))
              for d in ('DELTA_W','ZETA_C')}
        mass,origins=a.evolution_route(alpha,'DELTA_W',rho,route=route)
        constraint=a.finalize(a.virtual_constraint_route(alpha,'DELTA_W',rho,route=route))
        velocity={s:a.finalize(a.face_velocity_raw(source(alpha,s,'DELTA_W','RHO4_CONSTANT',route=route)))
                  for s in c.FACES}
    finally:
        a.build_face_source=original
    def axes(v): return sp.Tuple(*(sp.Integer(x) if isinstance(x,int) else Str(x) for x in v))
    def entry(name,values): return sp.Tuple(Str(name),sp.Tuple(*(sp.Tuple(axes(k),v) for k,v in values)))
    bundle=sp.Tuple(entry('traction',[((alpha,s,'DELTA_W',rho),traction[s]) for s in c.FACES]),
                    entry('closure_shape_deriv',[((alpha,s,'DELTA_W',rho),closure[s]) for s in c.FACES]),
                    entry('virtual_work_shape_deriv',[((alpha,'DELTA_W',d,rho),work[d]) for d in work]),
                    entry('virtual_constraint',[((alpha,'DELTA_W',rho),constraint)]),
                    entry('evolution_mass_balance',[((alpha,'DELTA_W',rho),mass)]),
                    entry('evolution_term_origins',[((alpha,'DELTA_W',rho),b.casify(origins))]))
    correction=sum(b.bind_mu_theta_operand(v,alpha,mu) for v in closure.values())
    corrected={k:v-correction if k=='TRUE_AREA_FACE_FLUX' else v for k,v in origins.items()}
    folded=b.face_generalized_force_rows(bundle,alpha,rho,
        sp.Tuple(*(sp.Tuple(Str(k),v) for k,v in corrected.items())),mu)
    rows=flatten({'U':folded['U'],'E_W':folded['E_W'],'THETA':mass-correction})
    provenance={'normal':{s:a.finalize(a.face_normal_raw(v)) for s,v in sources.items()},
                'traction':traction,'closure':closure,'work':work,'constraint':constraint,
                'origins':origins,'center':folded['CENTER_FACE_GENERALIZED_ROW'],'velocity':velocity}
    return rows,velocity,provenance


def constitutive(a,b,inputs,alpha,rho):
    progress('energy_source')
    energy=b.construct_energy(alpha,background_depth=3)
    pressure_slots=tuple(inputs.a(prefix+label) for label in ('plus','minus')
                         for prefix in ('delta_p_','d_w_delta_p_'))
    kinetic_sources=(*b.u_tt,b.e_tt,b.mu_W,b.W_bg,b.density_pair(rho)[1])
    emit('N6_BASE_SOURCE_DEPENDENCY',{'bulk_energy_terms':[[term.has(p) for p in pressure_slots] for _,term in energy.terms],
         'kinetic_operands':[[term.has(p) for p in pressure_slots] for term in kinetic_sources]})
    rho4=b.density_pair(rho)[0]
    grad=tuple(b.total_derivative(rho4,i,background_depth=3) for i in range(3))
    adv=b.dot(b.u,grad)/rho4
    t=sp.Symbol('n6AdvectionTag',real=True)
    c.DIMENSION_SCHEMA[t.name]=(0,0,0)
    # Tag enters ONLY the theta substitution of the quadratic source. To use
    # material_pullback verbatim, introduce the complementary source shift
    # before its native theta+advection substitution. The map leaves u fixed.
    shift=(t-1)*adv
    source_tag={b.theta:b.theta+shift,
                **{b.grad_theta[i]:b.grad_theta[i]+b.total_derivative(shift,i,background_depth=3)
                   for i in range(3)}}
    mu_e=[]; mu_m=[]
    for index,(_,term) in enumerate(energy.terms):
        bounded(term,'action_term')
        for monomial in sp.Add.make_args(term):
            pulled=b.material_pullback(monomial.subs(source_tag,simultaneous=True),alpha,rho,background_depth=3)
            # Exact second-return formula from operator_from_density, termwise;
            # pressure-independent U/E variations need not be materialized.
            def variation(density):
                return sp.diff(density,b.theta)-sum(b.total_derivative(
                    sp.diff(density,b.grad_theta[i]),i,background_depth=3) for i in range(3))
            mu_e.append(bounded(variation(monomial),'unpulled_mu_term'))
            mu_m.append(bounded(variation(pulled),'material_mu_term'))
        if index % 5 == 0: progress('constitutive_terms',term=index,total=len(energy.terms))
    # Explicit non-pressure jet boundary, preserving imported symbol identities.
    bridge={a.grad_theta[i]:b.grad_theta[i] for i in range(3)}
    for atom,dim in {**a.SYMBOL_DIMENSIONS,**b.SYMBOL_DIMENSIONS}.items():
        c.DIMENSION_SCHEMA.setdefault(atom.name,tuple(dim))
    c.dimension.cache_clear()
    emit('N6_JET_BRIDGE',{'identities':[(str(x),str(y),c.wave_jet(x)==c.wave_jet(y)) for x,y in bridge.items()]})
    emit('N6_ADVECTION_SOURCE',{'gradient_fingerprint':sha(grad),'advection_operations':sp.count_ops(adv),
                              'tag_derivative_operations':sp.count_ops(sp.diff(shift,t)),
                              'gradient_zero':[v==0 for v in grad]}, probe='N4_ADVECTION')
    return mu_e,mu_m,t,adv


def wave_terms(expressions,inputs):
    result={}
    for expression in expressions:
        expanded=sp.expand(expression)
        for term in sp.Add.make_args(expanded):
            waves=[s for s in term.free_symbols if s.name in c.WAVE_NAMES]
            if len(waves)!=1:
                if term==0: continue
                raise ValueError(('linear constitutive wave census',len(waves),str(term)[:100]))
            w=waves[0]
            coeff=term/w
            result[w]=result.get(w,0)+coeff
    return result


def source_terms(inputs,alpha,rho,face,mu,velocity):
    response=inputs.response[(alpha,face,rho)]
    z=c.atom_named(response,f's11cc1_dtn_operator_{alpha.lower()}_{"plus" if face==1 else "minus"}')
    r=c.named(response,'RESOLVENT')
    src=c.named(response,'DELTA_P').subs({inputs.values['rho_br_bg_rho4_constant']:inputs.density[(rho,)][1]},simultaneous=True)
    src=src.subs({z:1,r:1},simultaneous=True)/inputs.eps
    label='plus' if face==1 else 'minus'
    mus=c.atom_named(src,f's11cc1_mu_theta_{alpha.lower()}_{label}')
    vs=c.atom_named(src,f's11cc1_V_{alpha.lower()}_{label}')
    # The imported response source is differentiated, not retyped.
    mu_factor=sp.diff(src,mus); v_factor=sp.diff(src,vs)
    return wave_terms([mu_factor*x for x in mu]+[v_factor*velocity/inputs.eps],inputs)


def source_value(comp,terms):
    value={}
    for wave,coefficient in terms.items():
        value=gadd(value,gscale(comp.grades(coefficient,'Y'),comp.scalar(c.wave_jet(wave,c.Y))))
    return value


def pressure_table(comp,sources,kernels,normal_jet=False):
    out={}
    for face,terms in sources.items():
        source=source_value(comp,terms)
        for signature,kernel in kernels[face][0].items():
            if normal_jet: kernel=kernel*kernels[face][1]
            for grade,node in gmul(source,comp.grades(kernel)).items():
                out[face,signature,grade]=times(node,measure(signature))
    return out


def pressure_coefficients(rows,slots):
    zero=dict.fromkeys(slots,sp.S.Zero)
    return {(row,p):sp.diff(e,p).xreplace(zero) for row,e in rows.items() for p in slots}


@lru_cache(None)
def template(row,wave,signature,inputs):
    """Single Eulerian extract, on a one-slot/one-wave linear ansatz.

    Coefficients are typed placeholders. Their original circuits are bound to
    this template only after compilation, so no closed physical rows expand.
    """
    cf=sp.Function('n6CarrierCoefficient')(*c.X)
    c.NEW_DIMENSIONS[cf.func]=(0,0,0)
    value=c.wave_jet(wave,c.Y) if wave is not None else sp.S.One
    ko=tuple(inputs.a(f's11cc1_k_output_{j}') for j in range(1,4))
    ki=tuple(inputs.a(f's11cc1_k_input_{j}') for j in range(1,4))
    if signature:
        phase=sp.exp(sp.I*(sum(k*x for k,x in zip(ko,c.X))-sum(k*y for k,y in zip(ko if signature==6 else ki,c.Y))))
        variables=(*ko,*c.Y) if signature==6 else (*ko,*ki,*c.Y)
        if signature==12: variables=(*variables,*c.MIDDLE)
        operand=cf*sp.Integral(phase*value,*((v,-sp.oo,sp.oo) for v in variables))
    else:
        phase=sp.S.One; operand=cf
    rows=dict.fromkeys(ROWS,sp.S.Zero); rows[row]=operand
    weak=c.extract(unflatten(rows),inputs)
    TEMPLATE_TEST_SPATIAL_DERIVATIVES.append(sum(
        isinstance(d.expr,AppliedUndef) and d.expr.func.__name__.startswith('s11cc2Test')
        and any(v in (*c.X,*c.Y) for v in d.variables)
        for branch in weak.values() for expression in branch.values() for d in expression.atoms(sp.Derivative)))
    def strip(e):
        if isinstance(e,sp.Integral): return sp.expand(e.function/phase)
        if isinstance(e,sp.Add): return sp.Add(*(strip(x) for x in e.args))
        if isinstance(e,sp.Mul): return sp.Mul(*(strip(x) for x in e.args))
        return e
    return cf,{block:bounded(strip(weak[block[0]][block[1]]),'weak_template') for block in BLOCKS}


def kernel_coefficients(inputs,alpha,rho,face):
    progress('kernel_bridge',face=face)
    bridge=c.kernel_bridge(inputs,alpha,face,rho,{})
    pmat,ko,ki,qo,second=bridge[6],bridge[7],bridge[8],bridge[9],bridge[10]
    kernels={6:pmat[0,0],9:c.fourier_profiles(inputs,pmat[0,1],ko,ki),12:second}
    # Fourier first jets are derivative transforms of the same profile, with
    # the inherited dimensionless profile-jet scaling, not independent hats.
    for sig,e in kernels.items():
        replacements={}
        for f in e.atoms(AppliedUndef):
            name=f.func.__name__
            if 'ProfileJetHat' in name:
                direction=int(name[-1])-1
                base=sp.Function('s11cc2FourierW1ProfileHatTransfer')(*f.args)
                replacements[f]=sp.I*inputs.values['L_W']*f.args[direction]*base
        kernels[sig]=e.xreplace(replacements)/(2*sp.pi)**3
    reference=face*inputs.values['W_0']/2
    extension=sp.exp(sp.I*face*qo*(c.NORMAL-reference))
    multiplier=sp.diff(extension,c.NORMAL).subs(c.NORMAL,reference)
    return kernels,multiplier


def dimension_record(e,target):
    d=c.dimension(e)
    return {'computed':d,'target':target,'consistent':e==0 or d==target,
            'zero':e==0,'unknown':[str(s) for s in e.free_symbols if s.name not in c.DIMENSION_SCHEMA]}


def build_increment(comp,inputs,coeffs,sources,kernels,slots):
    output={}; dimensions={}; carrier_nodes={}
    for (row,p),coefficient in coeffs.items():
        coefficient=sp.cancel(coefficient/inputs.eps)
        face=1 if p.name.endswith('plus') else -1
        jet=p.name.startswith('d_w_')
        carrier_nodes[(row,p)]=comp.grades(coefficient)
        dims=dimension_record(coefficient,tuple(x-y for x,y in zip(ROW_DIM[row],c.dimension(p))))
        dimensions[(row,str(p))]=dims
        if coefficient == 0:
            continue
        # Differentiate before profiles or numbers; this is the live depth-4
        # weak derivative, including the product rule at the exterior point X.
        dcoeff=[inputs.dx(coefficient,i) for i in range(3)]
        cg=[comp.grades(coefficient)]+[comp.grades(e) for e in dcoeff]
        for signature in (0,6,9,12):
            waves={None:sp.S.One} if signature==0 else sources[face]
            kernel= -p if signature==0 else kernels[face][0][signature]*(kernels[face][1] if jet else 1)
            kg=comp.grades(kernel)
            for wave,source_coefficient in waves.items():
                sg=comp.grades(source_coefficient,'Y')
                cf,weak=template(row,wave,signature,inputs)
                for cg_grade in GRADES:
                    bindings={cf:cg[0].get(cg_grade,number(0))}
                    bindings.update({sp.diff(cf,c.X[i]):cg[i+1].get(cg_grade,number(0)) for i in range(3)})
                    for block,expression in weak.items():
                        w=comp.scalar(expression,binding=bindings)
                        products=gmul({cg_grade:w},gmul(sg,kg))
                        for grade,value in products.items():
                            key=(block,grade,signature,face)
                            output[key]=plus(output.get(key,number(0)),times(times(value,measure(signature)),amplitude()))
        progress('carrier_slot',row=row,slot=str(p))
    # Include both faces independently and their computed sum; include every
    # retained grade/signature/block even for a structurally empty circuit.
    for block,grade,signature in itertools.product(BLOCKS,GRADES,(0,6,9,12)):
        for face in c.FACES: output.setdefault((block,grade,signature,face),number(0))
        output[(block,grade,signature,0)]=plus(output[(block,grade,signature,1)],output[(block,grade,signature,-1)])
    return output,dimensions


def slot_guard(comp,rows,slots,coeffs):
    out={}; native={}; carrier={}
    zero=dict.fromkeys(slots,sp.S.Zero)
    for row,e in rows.items():
        # All four evaluations remain circuits; no symbolic closed slab.
        circuit=comp.grades(e)
        base=comp.grades(e.xreplace(zero))
        linear={}
        for p in slots: linear=gadd(linear,comp.grades(coeffs[row,p]*p))
        for grade in GRADES:
            key=(row,'linear',grade)
            native[key]=minus(circuit.get(grade,number(0)),base.get(grade,number(0)))
            carrier[key]=linear.get(grade,number(0))
            out[key]=minus(native[key],carrier[key])
        for p,q in itertools.combinations_with_replacement(slots,2):
            for grade,node in comp.grades(sp.diff(e,p,q)).items(): out[(row,'cross',str(p),str(q),grade)]=node
        for p in slots:
            denominators=[x.base for x in sp.preorder_traversal(e) if isinstance(x,sp.Pow) and x.exp.is_negative]
            out[(row,'denominator_dependency',str(p))]=number(sum(x.has(p) for x in denominators))
    return native,carrier,out


def closure_guard(comp,inputs,rows,coeffs,slots,sources,kernels,mu_slot,mu_terms):
    """P and chi(P) evaluations of the complete row circuit, plus its carrier.

    Formal integral signatures have independent probe weights in this guard;
    the primary six-block tables retain the signatures separately.
    """
    mu={}
    for term in mu_terms: mu=gadd(mu,comp.grades(term))
    chi={}
    for p in slots:
        face=1 if p.name.endswith('plus') else -1
        source={}
        for w,coef in sources[face].items():
            source=gadd(source,gscale(comp.grades(coef,'Y'),comp.scalar(c.wave_jet(w,c.Y))))
        value={}
        for signature,kernel in kernels[face][0].items():
            if p.name.startswith('d_w_'): kernel=kernel*kernels[face][1]
            weighted=gmul(comp.grades(kernel),source)
            weight=variable(('formal_integral_guard',signature),(6-signature,0,0))
            value=gadd(value,gscale(weighted,weight))
        chi[p]=value
    full={}; carried={}
    for face in (1,-1,0):
        selected=[p for p in slots if face==0 or p.name.endswith('plus' if face==1 else 'minus')]
        for row,e in rows.items():
            opening=comp.routed(e,{mu_slot:mu})
            closing=comp.routed(e,{mu_slot:mu,**{p:chi[p] for p in selected}})
            carrier={}
            for p in selected:
                delta=gsubtract(chi[p],comp.grades(p))
                carrier=gadd(carrier,gmul(comp.grades(coeffs[row,p]),delta))
            change=gsubtract(closing,opening)
            for grade in GRADES:
                key=(row,face,grade)
                full[key]=change.get(grade,number(0)); carried[key]=carrier.get(grade,number(0))
    return full,carried,residual(full,carried)


def residual(a,b):
    return {k:minus(a.get(k,number(0)),b.get(k,number(0))) for k in set(a)|set(b)}


class Sample:
    def __init__(self,prime,seed,cell,inputs,comp):
        self.p=prime; self.seed=seed; self.cell=cell; self.values={}; self.cache={}; self.fractions={}
        self.i=int(sp.sqrt_mod(-1,prime))
        self.values[('algebraic_i',)]=self.i
        omega=self.atom((inputs.values['omega'].name,'global'))
        cs=self.atom((inputs.values['c_s0'].name,'global'))
        h=omega*pow(cs,-1,prime)%prime
        # A real lift has |v_j|<1/16, omega sign cell[0], cs>0.
        # Sphere and hyperboloid stereographic charts respect q^2+k^2=h^2.
        # Reduction uses scale 16*(p+1), invertible mod p.
        self.values[(inputs.values['omega'].name,'global')]=cell[0]*omega%prime
        h=cell[0]*h%prime
        self.momenta={}; self.dispersion=[]
        for leg,branch in zip(('output','input','middle'),cell[1:]):
            v=[self.atom(('chart',leg,j))*pow(16,-1,prime)%prime for j in range(3)]
            if branch==0:
                denominator=(1+sum(x*x for x in v))%prime
                if not denominator: raise ZeroDivisionError('chart')
                k=[2*h*x*pow(denominator,-1,prime)%prime for x in v]
                q=h*(1-sum(x*x for x in v))*pow(denominator,-1,prime)%prime
            else:
                v[2]=-cell[0]*v[2]%prime
                denominator=(1+v[0]**2+v[1]**2-v[2]**2)%prime
                if not denominator: raise ZeroDivisionError('chart')
                lam=-2*h*pow(denominator,-1,prime)%prime
                k=[(h+lam)%prime,lam*v[0]%prime,lam*v[1]%prime]
                q=self.i*lam*v[2]%prime
            self.momenta[leg]=(k,q)
            self.dispersion.append((q*q+sum(x*x for x in k)-h*h)%prime)
            symbols=tuple(inputs.a(f's11cc1_k_{leg}_{j}') for j in range(1,4)) if leg!='middle' else c.MIDDLE
            qs=inputs.a('s11cc1_q_out_'+leg) if leg!='middle' else c.MIDDLE_Q
            for sym,val in zip((*symbols,qs),(*k,q)): self.values[(sym.name,'global')]=val
        if any(self.momenta[x][0]==self.momenta[y][0] for x,y in itertools.combinations(self.momenta,2)):
            raise ZeroDivisionError('coincident_momenta')

    def atom(self,key):
        if key not in self.values:
            counter=0
            while True:
                digest=hashlib.sha256(repr((self.seed,key,counter)).encode()).digest()
                value=int.from_bytes(digest,'big') & ((1<<(self.p-1).bit_length())-1)
                if value<self.p-1: break
                counter+=1
            self.values[key]=1+value
        return self.values[key]

    def evaluate(self,node):
        stack=[node]
        while stack:
            current=stack[-1]
            if current in self.cache:
                stack.pop(); continue
            op=current.op; args=current.args
            missing=[a for a in args if isinstance(a,Node) and a not in self.cache]
            if missing:
                stack.extend(missing); continue
            if op in ('number','variable'):
                numerator=args[0]%self.p if op=='number' else self.atom(args[0])
                denominator=1
            elif op in ('add','mul'):
                an,ad=self.fractions[args[0]]; bn,bd=self.fractions[args[1]]
                numerator=(an*bd+bn*ad if op=='add' else an*bn)%self.p
                denominator=ad*bd%self.p
            elif op=='pow':
                an,ad=self.fractions[args[0]]; exponent=args[1]
                if exponent<0: an,ad=ad,an
                numerator=pow(an,abs(exponent),self.p)
                denominator=pow(ad,abs(exponent),self.p)
            else: raise ValueError(op)
            self.fractions[current]=(numerator,denominator)
            self.cache[current]=numerator*pow(denominator,-1,self.p)%self.p
            stack.pop()
        return self.cache[node]


def census(roots):
    seen=set(); stack=list(roots); excluded=0; variables=set()
    while stack:
        node=stack.pop()
        if node in seen: continue
        seen.add(node)
        if node.op=='variable': variables.add(node.args[0])
        if node.op=='pow' and node.args[1]<0: excluded+=node.args[0].degree[0]
        stack.extend(a for a in node.args if isinstance(a,Node))
    return len(seen),excluded,variables


def pit(objects,inputs,comp,draws,seed):
    primes=(1000000009,998244353,1000000033)
    cells=tuple((sign,*bs) for sign in (-1,1) for bs in itertools.product((0,1),repeat=3))
    roots=[v for values in objects.values() for v in values.values()]
    nodes,E,variables=census(roots)
    # Twelve on-shell coordinates share three quadratic chart denominators
    # and c_s0. Each coordinate numerator has degree <=4. A common product of
    # the twelve possible denominator expressions has degree <=36. Hence 128
    # times (numerator_degree + denominator_degree) is conservative after
    # substitution, including excluded inverses and all three leg charts.
    bounds={name:{str(k):128*sum(v.degree) for k,v in values.items()} for name,values in objects.items()}
    E=128*E+100
    D=max((d for values in bounds.values() for d in values.values()),default=0)
    family=sum(len(values) for values in objects.values())*len(cells)
    if min(primes)-1-E <= 16*D:
        generated=[]; candidate=max(max(primes),16*(D+E+1))
        while len(generated)<3:
            candidate=int(sp.nextprime(candidate))
            if candidate%4==1: generated.append(candidate)
        primes=tuple(generated)
    requested_draws=draws
    rate=max(sp.Rational(D,p-1-E) for p in primes)
    while family*rate**draws > sp.Rational(1,2**64): draws+=1
    per_prime=[min(sp.S.One,sp.Rational(D,p-1-E))**draws if E<p-1 else sp.S.One for p in primes]
    delta=min(sp.S.One,family*max(per_prime))
    emit('N6_PRIMES',{'values':primes,'prime_tests':[sp.isprime(p) for p in primes],
                       'i_roots':[sp.sqrt_mod(-1,p) for p in primes]})
    emit('N6_PIT_PROVENANCE',{'seed':seed,'primes':primes,'draws_per_prime_cell':draws,
        'cells':cells,'D':D,'E':E,'family_size':family,'delta':delta,'nodes':nodes,'requested_draws':requested_draws,
        'conditional_good_prime':True,'formal_kernel_only':True,'hash_sampler_random_oracle_condition':True,
        'chart_scale_multiples':[16*(p+1) for p in primes], 'degree_bounds':bounds,
        'excluded_transition_equations':['omega=0','q=0','chart_denominator=0','response_denominator=0','k_in=k_out']})
    ko=tuple(inputs.a(f's11cc1_k_output_{j}') for j in range(1,4))
    ki=tuple(inputs.a(f's11cc1_k_input_{j}') for j in range(1,4))
    signatures={0:(),6:(*ko,*c.Y),9:(*ko,*ki,*c.Y),12:(*ko,*ki,*c.Y,*c.MIDDLE)}
    emit('N6_FORMAL_KERNEL_SCHEMA',{'bound_variables':{k:[str(v) for v in vs] for k,vs in signatures.items()},
         'measure_dimensions':{k:tuple(sum(c.dimension(v)[j] for v in vs) for j in range(3)) for k,vs in signatures.items()},
         'spatial_test_derivative_count':sum(TEMPLATE_TEST_SPATIAL_DERIVATIVES),
         'temporal_ibp':False,'fourier_phase_evaluated':False,'integrals_evaluated':False,
         'epsilon_homogeneous_columns':sum(len(n.wave_degrees)<=1 for n in roots),'columns':len(roots)})
    schemas={}
    for name,table in objects.items():
        schemas[name]={}
        for key,node in table.items():
            if name[0]=='N6_PRESSURE_AMPLITUDE': target=(-2,-2,1)
            elif name[0]=='N6_NORMAL_JET_AMPLITUDE': target=(-3,-2,1)
            elif name[0]=='N6_FACE_VELOCITY': target=(1,-1,0)
            elif name[0]=='N6_FACE_NORMAL': target=(0,0,0)
            elif name[0]=='N6_MU_AMPLITUDE': target=(-1,-2,1)
            elif name[0]=='N6_PRESSURE_SOURCE_TAG_SENSITIVITY':
                target=tuple(a-b for a,b in zip((1,-1,0),c.dimension(inputs.a(key[1]))))
            elif name[0]=='N6_CARRIER_SOURCE_TAG_CONTRACTION':
                target=tuple(a-b+v-w for a,b,v,w in zip(ROW_DIM[key[0]],c.dimension(inputs.a(key[1])),
                                  (1,-1,0),c.dimension(inputs.a(key[2]))))
            elif isinstance(key[0],tuple) and key[0] in WEAK_DIM:
                target=WEAK_DIM[key[0]]
            elif key[0] in ROW_DIM:
                target=ROW_DIM[key[0]]
                if len(key)>1 and isinstance(key[1],str) and key[1] in inputs.atoms:
                    target=tuple(a-b for a,b in zip(target,c.dimension(inputs.a(key[1]))))
                if len(key)>1 and key[1]=='cross':
                    for slot in key[2:4]: target=tuple(a-b for a,b in zip(target,c.dimension(inputs.a(slot))))
                if len(key)>1 and key[1]=='denominator_dependency': target=(0,0,0)
            elif isinstance(key[0],str) and key[0] in inputs.atoms:
                target=tuple(a-b for a,b in zip((-1,-2,1),c.dimension(inputs.a(key[0]))))
            else: target=None
            schemas[name][str(key)]={'computed':sorted(node.units,key=str),'target':target,
                'consistent':not node.units or (target is not None and node.units==frozenset((target,))),
                'epsilon_support':sorted(node.wave_degrees)}
    rejected=0
    columns={name:sorted(table,key=str) for name,table in objects.items()}
    tables={name:[] for name in objects}
    sample_catalog=[]
    for pi,p in enumerate(primes):
        for ci,cell in enumerate(cells):
            for draw in range(draws):
                attempt=0
                while True:
                    sample_seed=(seed,pi,ci,draw,attempt)
                    try:
                        sample=Sample(p,sample_seed,cell,inputs,comp)
                        # Jointly evaluate before emitting any compared operand.
                        for root in roots: sample.evaluate(root)
                        break
                    except (ZeroDivisionError,ValueError) as exc:
                        if isinstance(exc,ValueError) and 'invertible' not in str(exc): raise
                        rejected+=1; attempt+=1
                        if attempt>200: raise RuntimeError('singular sample exhaustion') from exc
                sample_catalog.append({'prime':p,'cell':cell,'draw':draw,'seed':sample_seed,
                    'dispersion_numerators':sample.dispersion,'rejected_joint':attempt,'variables':len(sample.values)})
                for name in objects:
                    tables[name].append([sample.fractions[objects[name][key]] for key in columns[name]])
            progress('pit_cell',prime=p,cell=cell)
    emit('N6_SAMPLE_PROVENANCE',sample_catalog)
    for (name,probe),table in tables.items():
        keys=columns[name,probe]
        support=set()
        for key in keys:
            grade=next((part for part in reversed(key) if isinstance(part,tuple) and len(part)==2
                        and all(isinstance(x,int) for x in part)),(0,0))
            support.update((e,*grade) for e in objects[name,probe][key].wave_degrees)
        emit(name,{'columns':keys,'numerator_denominator':table,
             'nonzero_modular_numerator':[any(row[i][0]!=0 for row in table) for i in range(len(keys))]},probe=probe,
             dimension=[schemas[name,probe][str(key)] for key in keys],
             N12=sorted(support), declared_retained_rectangle=[[1,*g] for g in GRADES])
    emit('N6_PIT_REJECTIONS',{'rejected_joint':rejected,'accepted':len(primes)*len(cells)*draws})


def run(args):
    global CASE
    progress('imports')
    import S11c_a_interface_geometry_sympy_audit as a
    import S11c_b_brane_operator_sympy_audit as b
    fold,_=c.load_model(str(ROOT/'scripts/S11c_b_exports.py'),str(ROOT/'scripts/S11c_c1_exports.py'))
    inputs=c.bind_inputs(fold)
    emit('N6_PREMISES',{'supplied':['5c','3a','3c','c1_face_response','a_material_sources','b_material_scalar'],
                       'withheld':['residual_expected_value','control_acceptance_criterion'],
                       'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                            [Path(__file__),Path(a.__file__),Path(b.__file__),Path(c.__file__),
                             ROOT/'_measurements/S11c_c2_N6_route2_spec_astra.md']}})
    slots=tuple(inputs.a(prefix+label) for label in ('plus','minus') for prefix in ('delta_p_','d_w_delta_p_'))
    for alpha,rho in itertools.product(c.ANCHORINGS,c.DENSITIES):
        if args.anchoring and alpha!=args.anchoring: continue
        if args.density and rho!=args.density: continue
        CASE={'anchoring':alpha,'density':rho}; progress('case')
        comp=Compiler(inputs)
        emit('N6_PRESSURE_IDENTITY',{'census':[(str(p),p==(a.delta_p_plus if s==1 else a.delta_p_minus) if not p.name.startswith('d_w_') else p==a.dw_delta_p[s],p.assumptions0)
               for s in c.FACES for p in slots if p.name.endswith('plus' if s==1 else 'minus')]})
        imported=flatten(c.expanded_rows(inputs.slab[alpha,rho]))
        e_coeff=pressure_coefficients(imported,slots)
        unpulled,material,t,adv=constitutive(a,b,inputs,alpha,rho)
        mu_e=[inputs.mu[alpha,rho][1]/inputs.eps]
        mu_m=[e.subs(t,1) for e in material]
        progress('native_faces')
        # A typed constitutive circuit slot binds before virtual differentiation;
        # its circuit is specialized at both binding sites in each evaluation.
        mu_slot=sp.Symbol('n6MaterialMu',real=True); c.DIMENSION_SCHEMA[mu_slot.name]=(-1,-2,1)
        face_mu_bindings={mu_slot:mu_m}
        m_rows,m_v,m_prov=face_factory(a,b,inputs,alpha,rho,'MATERIAL',mu_slot)
        e_rows,e_v,e_prov=face_factory(a,b,inputs,alpha,rho,'EULERIAN',mu_slot)
        tilt_rows,_,tilt_prov=face_factory(a,b,inputs,alpha,rho,'EULERIAN',mu_slot,tilt=-1)
        m_coeff=pressure_coefficients(m_rows,slots)
        f_coeff=pressure_coefficients(e_rows,slots)
        tilt_coeff=pressure_coefficients(tilt_rows,slots)
        for label,provenance in [('MATERIAL',m_prov),('EULERIAN_FACTORY',e_prov),('TILT',tilt_prov)]:
            emit('N6_FACE_PROVENANCE',{'fingerprint':sha(provenance),
                'constraint_pressure_dependency':[provenance['constraint'].has(p) for p in slots],
                'origin_pressure_dependency':{k:[v.has(p) for p in slots] for k,v in provenance['origins'].items()},
                'mu_carrier_dependency':[v.has(mu_slot) for v in (m_coeff if label=='MATERIAL' else f_coeff).values()]},route=label)
        objects={}
        def add(name,table,probe=None): objects[name,probe]=table
        imported_mu=wave_terms(mu_e,inputs); reconstructed_mu=wave_terms(unpulled,inputs)
        def mu_table(waves):
            return {(str(w),g):v for w,e in waves.items() for g,v in comp.grades(e,'Y').items()}
        mt=mu_table(imported_mu); rt=mu_table(reconstructed_mu)
        add('N6_MU_RECONSTRUCTION_IMPORTED',mt); add('N6_MU_RECONSTRUCTION_NATIVE',rt)
        add('N6_MU_RECONSTRUCTION_RESIDUAL',residual(mt,rt))
        for route,mu,velocities,provenance in [('EULERIAN',mu_e,e_v,e_prov),('MATERIAL',face_mu_bindings[mu_slot],m_v,m_prov)]:
            add('N6_MU_AMPLITUDE',{('mu',g):v for g,v in source_value(comp,wave_terms(mu,inputs)).items()},route)
            add('N6_FACE_VELOCITY',{(s,g):v for s,e in velocities.items()
                for g,v in source_value(comp,wave_terms([e],inputs)).items()},route)
        for route,provenance in [('EULERIAN',e_prov),('MATERIAL',m_prov),('TILT',tilt_prov)]:
            add('N6_FACE_NORMAL',{(s,order,component,g):v for s,normal in provenance['normal'].items()
                for order in (0,1) for component,e in enumerate(normal[order]) for g,v in comp.grades(e).items()},route)
        def coefficient_table(coeffs):
            return {(row,str(p),g):n for (row,p),e in coeffs.items() for g,n in comp.grades(e).items()}
        ct=coefficient_table(e_coeff); ft=coefficient_table(f_coeff)
        add('N6_CARRIER_RECONSTRUCTION_IMPORTED',ct); add('N6_CARRIER_RECONSTRUCTION_NATIVE',ft)
        add('N6_CARRIER_RECONSTRUCTION_RESIDUAL',residual(ct,ft))
        for route,rows,coeffs in [('EULERIAN',imported,e_coeff),('MATERIAL',m_rows,m_coeff)]:
            for label,table in zip(('NATIVE','CARRIER','RESIDUAL'),slot_guard(comp,rows,slots,coeffs)):
                add('N6_SLOT_GUARD_'+label,table,route)
        kernels={s:kernel_coefficients(inputs,alpha,rho,s) for s in c.FACES}
        ev={s:inputs.geometry['face_velocity'][alpha,s,'DELTA_W'] for s in c.FACES}
        es={s:source_terms(inputs,alpha,rho,s,mu_e,ev[s]) for s in c.FACES}
        ms={s:source_terms(inputs,alpha,rho,s,mu_m,m_v[s]) for s in c.FACES}
        for route,sources in [('EULERIAN',es),('MATERIAL',ms)]:
            add('N6_PRESSURE_AMPLITUDE',pressure_table(comp,sources,kernels),route)
            add('N6_NORMAL_JET_AMPLITUDE',pressure_table(comp,sources,kernels,True),route)
        for route,rows,coeffs,sources,mu in [('EULERIAN',imported,e_coeff,es,mu_e),('MATERIAL',m_rows,m_coeff,ms,mu_m)]:
            for label,table in zip(('NATIVE','CARRIER','RESIDUAL'),closure_guard(comp,inputs,rows,coeffs,slots,sources,kernels,mu_slot,mu)):
                add('N6_CLOSURE_GUARD_'+label,table,route)
        E,edim=build_increment(comp,inputs,e_coeff,es,kernels,slots)
        M,mdim=build_increment(comp,inputs,m_coeff,ms,kernels,slots)
        add('REP_INVARIANCE_EULERIAN_OPERAND',E); add('REP_INVARIANCE_MATERIAL_OPERAND',M)
        add('REP_INVARIANCE_RESIDUAL',residual(E,M))
        F,_=build_increment(comp,inputs,f_coeff,es,kernels,slots)
        T,_=build_increment(comp,inputs,tilt_coeff,es,kernels,slots)
        add('CONTROL_INDEPENDENCE_BASE',F,'TILT'); add('CONTROL_INDEPENDENCE_CORRUPTED',T,'TILT')
        add('CONTROL_INDEPENDENCE_RESIDUAL',residual(T,F),'TILT')
        add('CONTROL_INDEPENDENCE_COMPANION_BASE',M,'TILT')
        add('CONTROL_INDEPENDENCE_COMPANION_CORRUPTED',M,'TILT')
        if rho=='RHOBR_CONSTANT':
            mu_zero=[e.subs(t,0) for e in material]
            mu_zero_slot=sp.Symbol('n6MaterialMuAdvection0',real=True)
            c.DIMENSION_SCHEMA[mu_zero_slot.name]=(-1,-2,1)
            face_mu_bindings[mu_zero_slot]=mu_zero
            zero_rows,zero_v,zero_prov=face_factory(a,b,inputs,alpha,rho,'MATERIAL',mu_zero_slot)
            zero_coeff=pressure_coefficients(zero_rows,slots)
            zero_sources={s:source_terms(inputs,alpha,rho,s,mu_zero,m_v[s]) for s in c.FACES}
            M0,_=build_increment(comp,inputs,zero_coeff,zero_sources,kernels,slots)
            add('CONTROL_INDEPENDENCE_BASE',M,'N4_ADVECTION'); add('CONTROL_INDEPENDENCE_CORRUPTED',M0,'N4_ADVECTION')
            add('CONTROL_INDEPENDENCE_RESIDUAL',residual(M0,M),'N4_ADVECTION')
            add('CONTROL_INDEPENDENCE_COMPANION_BASE',E,'N4_ADVECTION')
            add('CONTROL_INDEPENDENCE_COMPANION_CORRUPTED',E,'N4_ADVECTION')
            add('N6_MU_TAG_SENSITIVITY',mu_table(wave_terms([sp.diff(e,t) for e in material],inputs)),'N4_ADVECTION')
            source_sensitivity={}
            contraction={}
            for s in c.FACES:
                for w in set(ms[s])|set(zero_sources[s]):
                    term=ms[s].get(w,0)-zero_sources[s].get(w,0)
                    sg=comp.grades(term,'Y')
                    for g,n in sg.items(): source_sensitivity[(s,str(w),g)]=n
                    for (row,p),coefficient in m_coeff.items():
                        if not p.name.endswith('plus' if s==1 else 'minus'): continue
                        cg=comp.grades(coefficient/inputs.eps)
                        for g,v in cg.items():
                            for h,z in sg.items():
                                grade=tuple(i+j for i,j in zip(g,h))
                                key=(row,str(p),str(w),grade)
                                contraction[key]=plus(contraction.get(key,number(0)),times(v,z))
            add('N6_PRESSURE_SOURCE_TAG_SENSITIVITY',source_sensitivity,'N4_ADVECTION')
            add('N6_CARRIER_SOURCE_TAG_CONTRACTION',contraction,'N4_ADVECTION')
            emit('N6_ADVECTION_BINDING',{'baseline_mu_fingerprint':sha(mu_m),
                 'corrupted_mu_fingerprint':sha(mu_zero),'baseline_face_fingerprint':sha(m_prov),
                 'corrupted_face_fingerprint':sha(zero_prov),'t_values':[1,0]},probe='N4_ADVECTION')
        else:
            emit('N6_ADVECTION_DEPENDENCY',{'density_gradient_zero':tuple(b.total_derivative(b.density_pair(rho)[0],i)==0 for i in range(3)),
                 'mu_tag_dependency':[e.has(t) for e in material]},probe='N4_ADVECTION')
        emit('N6_DIMENSIONS',{'eulerian':edim,'material':mdim,
             'failure_operand':dimension_record(sp.Add(slots[0],inputs.values['W_0']*slots[0],evaluate=False),c.dimension(slots[0]))})
        pit(objects,inputs,comp,args.draws,args.seed)
        progress('case_finished')
    progress('finished')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--draws',type=int,default=8)
    parser.add_argument('--seed',type=int,default=110602)
    parser.add_argument('--anchoring',choices=c.ANCHORINGS)
    parser.add_argument('--density',choices=c.DENSITIES)
    args=parser.parse_args()
    if args.draws<8: parser.error('at least eight valid draws per prime and cell')
    try: run(args)
    except BlockSize: return 2
    return 0


if __name__=='__main__':
    sys.exit(main())
