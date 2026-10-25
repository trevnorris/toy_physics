#!/usr/bin/env python3
"""S15: direct Cartesian tensor construction and polynomial support grouping.

The supplied equations, radial fields, split, and spherical chart are the only
physical inputs. All contractions, derivatives, partitions and integrals below
are computed. Wolfram code and output are not inputs. No value-based guard.
"""
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import sys
import sympy as S

OUT = Path(__file__).resolve().parent / 'out' / 'S15_two_drain_region.out'


def main(log):
    def emit(name, obj):
        log.write(f'{name}: {obj}\n')
        log.flush()

    def progress(name):
        print(f'WORKING_ON: {name}', flush=True)
        emit('WORKING_ON', name)

    def solo(name, obj):
        emit(name, obj)
        emit(name + '/comparison', 'There is nothing to compare it with in this engine.')

    def pair(name, a, b):
        emit(name + '/operand_a', a)
        emit(name + '/operand_b', b)
        residual = (a-b).applyfunc(S.expand) if isinstance(a, S.MatrixBase) else S.expand(a-b)
        emit(name + '/residual', residual)

    progress('supplied fields and domain')
    t = S.Symbol('t', real=True)
    xyz = S.symbols('x_1:4', real=True)
    point = S.Matrix(xyz)
    theta, phi = S.symbols('theta phi', real=True)
    radius = S.Symbol('s', positive=True)
    rin = S.Function('r_in')(t)
    rout = S.Function('R')(t)
    marks = S.symbols('label_1 label_2 label_r')
    labels = ('drain1', 'drain2', 'remainder')
    subsets = [frozenset(c) for k in range(1, 4) for c in combinations(range(3), k)]
    names = {q: '+'.join(labels[j] for j in sorted(q)) for q in subsets}
    centres = [S.Matrix([S.Function(f'X_{k}_{i}')(t) for i in range(1,4)]) for k in (1,2)]
    distances = [S.sqrt((point-c).dot(point-c)) for c in centres]
    rho_heads = [S.Function(f'rho_{k}') for k in (1,2)]
    u_heads = [S.Function(f'u_{k}') for k in (1,2)]
    ds = [head(r,t) for head,r in zip(rho_heads,distances)]
    vs = [head(r,t)*(point-c)/r for head,r,c in zip(u_heads,distances,centres)]
    args = (*xyz,t)
    ds.append(S.Function('varrho')(*args))
    vs.append(S.Matrix([S.Function(f'W_{i}')(*args) for i in range(1,4)]))
    density = sum(a*d for a,d in zip(marks,ds))
    velocity = sum((a*v for a,v in zip(marks,vs)), S.zeros(3,1))
    pressure = S.Function('p_br')
    eta, zeta = [S.Function(name)(*args) for name in ('eta','zeta')]
    other = S.Matrix(3,3,lambda i,j:S.Function(f'Pi_o_{i+1}_{j+1}')(*args))
    force = S.Matrix([S.Function(f'f_{i}')(*args) for i in range(1,4)])
    gradient = velocity.jacobian(point)
    stress = -pressure(density)*S.eye(3) + eta*(gradient+gradient.T-S.Rational(2,3)*S.eye(3)*S.trace(gradient)) + zeta*S.eye(3)*S.trace(gradient)
    G = density*velocity
    Pi = density*velocity*velocity.T - stress + other
    strip = lambda e: e.xreplace(dict.fromkeys(marks,S.Integer(1)))
    separation = S.sqrt((centres[0]-centres[1]).dot(centres[0]-centres[1]))
    domain = S.And(separation>rout, rout>rin, rin>0)
    emit('DOMAIN_FOR_EVERY_OBJECT', (domain,'both cores excluded; projected slab domain of S14a part 1; drain 1 outside outer surface'))
    emit('PRESSURE_DOMAIN', 'Whole pressure requires p_br defined at total density. Expansion requires p_br defined and differentiable at rho_2+varrho. Subtraction requires both. Conditions hold pointwise on the relevant integration surface.')
    emit('REGULARITY', 'Classical derivatives, integrals and moving-domain transport are used wherever they exist; no limiting regime is chosen.')
    emit('PREMISES', 'Supplied P5 and P6 applied to projected fields; P6 adopted as steady is used as written for moving drains; departures reside in Pi_o. f includes all face/bulk momentum exchange and stays undivided.')
    emit('CHOICES', 'Radial fields about moving centres; labelled split with unrestricted remainders; moving concentric spherical control surfaces; first-order pressure expansion where differentiable.')
    emit('BOOKKEEPING', 'Dimensionless label markers track factors only. Marker coefficient extraction is not an assignment or a freeze of any physical field.')
    solo('c_comp_not_present_in_supplied_equations', S.Function('c_comp')(*args))
    for name,obj in [('distances',distances),('density_profiles_and_remainder',ds),('velocity_profiles_and_remainder',vs),('total_density',strip(density)),('total_velocity',strip(velocity)),('P5_G',strip(G)),('P6_T',strip(stress)),('Pi',strip(Pi))]:
        solo(name,obj)
    divergence = S.Matrix([sum(S.diff(Pi[i,j],xyz[j]) for j in range(3)) for i in range(3)])
    local_lhs = strip(G.diff(t)+divergence)
    pair('supplied_local_balance',local_lhs,force)
    emit('supplied_local_balance/equation', S.Eq(local_lhs,force,evaluate=False))

    # Decompose the supplied tensor through its expression tree. No channel's
    # physical formula is separately supplied, and no live field is zeroed.
    pcall = pressure(density)
    def channel(term):
        if term.has(pressure):
            return 'pressure'
        if term.has(eta,zeta):
            return 'viscous'
        if any(term.has(q) for q in other):
            return 'other'
        return 'carried'
    channels = {k:S.zeros(3,3) for k in ('carried','pressure','viscous','other')}
    for i,j in product(range(3),repeat=2):
        for term in S.Add.make_args(S.expand(Pi[i,j], deep=False)):
            channels[channel(term)][i,j] += term
    pair('Pi_channel_partition',Pi,sum(channels.values(),S.zeros(3,3)))

    # Direct monomial-exponent partition, including all nonempty label sets.
    def groups(expr):
        buckets = {q:S.Integer(0) for q in subsets}
        for powers,coefficient in S.Poly(S.expand(expr),*marks).terms():
            support = frozenset(j for j,power in enumerate(powers) if power)
            if support:
                buckets[support] += coefficient
        emit('partition/unlabelled_coefficient',S.Poly(S.expand(expr),*marks).coeff_monomial((0,0,0)))
        return buckets

    # Differentiate in a formal amplitude; retain only its value and first
    # derivative at the expansion point. No second differentiability assumed.
    eps = S.Symbol('epsilon', real=True)
    amplitude_pressure = pcall.xreplace({marks[0]:eps*marks[0]})
    p0 = amplitude_pressure.subs(eps,0)
    p1 = S.diff(amplitude_pressure,eps).subs(eps,0)
    expansion = p0+p1
    p_groups = {q:S.Integer(0) for q in subsets}
    for term in S.Add.make_args(S.expand(expansion)):
        support = frozenset(j for j,m in enumerate(marks) if term.has(m))
        p_groups[support] += strip(term)
    solo('pressure_whole',strip(pcall))
    solo('pressure_expansion_at_order_one',strip(expansion))
    emit('pressure_expansion/order_and_condition', 'Labelled first-order choice in rho_1 about rho_2+varrho; p_br defined and differentiable there. The local remainder is o(rho_1) as rho_1 tends to zero; no quadratic bound or convergence regime is asserted.')
    solo('pressure_whole_minus_expansion',strip(pcall-expansion))
    pair('pressure_expansion_partition',strip(expansion),sum(p_groups.values()))

    # The spatial chart is supplied; measures and velocities are differentiated.
    e = S.Matrix([S.sin(theta)*S.cos(phi),S.sin(theta)*S.sin(phi),S.cos(theta)])
    def chart(b):
        return centres[1]+b*e
    def measure(mapping,variables):
        jac = mapping.jacobian(S.Matrix(variables))
        return S.trigsimp(S.sqrt(S.trigsimp((jac.T*jac).det()))).replace(S.Abs(S.sin(theta)),S.sin(theta))
    volume_measure = measure(chart(radius),(radius,theta,phi))
    def pull(expr,b):
        return expr.subs(dict(zip(xyz,chart(b))),simultaneous=True)
    def surface_integral(expr,b,area):
        return S.Integral(pull(expr,b)*area,(theta,0,S.pi),(phi,0,2*S.pi))
    def volume_integral(expr):
        return S.Integral(pull(expr,radius)*volume_measure,(radius,rin,rout),(theta,0,S.pi),(phi,0,2*S.pi))

    # Dimensional homogeneity is computed on expression trees. Density retains
    # an abstract independent dimension; no mass or normalization is inferred.
    Z=(S.Integer(0),)*3
    L=(S.Integer(0),S.Integer(1),S.Integer(0))
    T=(S.Integer(0),S.Integer(0),S.Integer(1))
    D=(S.Integer(1),S.Integer(0),S.Integer(0))
    add=lambda a,b:tuple(x+y for x,y in zip(a,b))
    scale=lambda a,n:tuple(x*n for x in a)
    sub=lambda a,b:add(a,scale(b,-1))
    dimensions={**dict.fromkeys(xyz,L),t:T,radius:L,theta:Z,phi:Z,eps:Z,**dict.fromkeys(marks,Z)}
    signatures={}
    def register(head,out,inputs):
        signatures[head]=(out,tuple(inputs))
    for c in centres:
        for q in c:
            register(q.func,L,(T,))
    for head in u_heads: register(head,sub(L,T),(L,T))
    for head in rho_heads: register(head,D,(L,T))
    for q in vs[-1]: register(q.func,sub(L,T),(L,L,L,T))
    register(ds[-1].func,D,(L,L,L,T))
    pd=add(D,scale(sub(L,T),2))
    vd=add(D,sub(scale(L,2),T))
    for q in (eta,zeta): register(q.func,vd,(L,L,L,T))
    for q in other: register(q.func,pd,(L,L,L,T))
    for q in force: register(q.func,add(D,sub(L,scale(T,2))),(L,L,L,T))
    register(pressure,pd,(D,))
    for q in (rin,rout): register(q.func,L,(T,))
    register(S.Function('c_comp'),sub(L,T),(L,L,L,T))
    argument_audit={}
    @lru_cache(None)
    def dim(expr, bindings=()):
        env=dict(bindings)
        if expr == 0: return frozenset()
        if expr.is_Number or expr in (S.pi,): return frozenset((Z,))
        if isinstance(expr,S.Symbol): return frozenset((env.get(expr,dimensions.get(expr)),)) if expr in env or expr in dimensions else unknown(expr)
        if isinstance(expr,S.Subs):
            more={**env}
            for variable,value in zip(expr.variables,expr.point):
                candidates=dim(value,bindings)
                if len(candidates)!=1:
                    emit('dimension_substitution_operands',(expr,candidates))
                    raise ValueError('Substitution dimension is not uniquely specified')
                more[variable]=next(iter(candidates))
            return dim(expr.expr,tuple(more.items()))
        if isinstance(expr,S.Derivative):
            result=dim(expr.expr,bindings)
            for variable,count in expr.variable_count:
                result=frozenset(sub(a,scale(b,count)) for a in result for b in dim(variable,bindings))
            return result
        if isinstance(expr,S.Integral):
            result=dim(expr.function,bindings)
            for variable,*bounds in expr.limits:
                vdims=dim(variable,bindings)
                for bound in bounds:
                    argument_audit[('integral_bound',str(variable),str(bound))]=(vdims,dim(bound,bindings))
                result=frozenset(add(a,b) for a in result for b in vdims)
            return result
        if expr.is_Add:
            return frozenset().union(*(dim(q,bindings) for q in expr.args))
        if expr.is_Mul:
            alternatives=[dim(q,bindings) for q in expr.args]
            return frozenset(tuple(sum(v[j] for v in row) for j in range(3)) for row in product(*alternatives))
        if expr.is_Pow:
            return frozenset(scale(a,expr.exp) for a in dim(expr.base,bindings))
        if expr.func in (S.sin,S.cos,S.Abs):
            if expr.func==S.Abs: return dim(expr.args[0],bindings)
            argument_audit[(str(expr.func),str(expr.args))]=(frozenset((Z,)),dim(expr.args[0],bindings))
            return frozenset((Z,))
        if expr.func in signatures:
            out,inputs=signatures[expr.func]
            for j,(arg,expected) in enumerate(zip(expr.args,inputs)):
                argument_audit[(str(expr.func),j,str(arg))]=(frozenset((expected,)),dim(arg,bindings))
            return frozenset((out,))
        return unknown(expr)
    def unknown(expr):
        emit('dimension_unregistered_operand',expr)
        raise ValueError(f'Unregistered dimension node: {type(expr).__name__}')
    def dimension_audit(name,expr):
        candidates=sorted(dim(expr))
        emit(name+'/dimension_operands',candidates)
        emit(name+'/dimension_pair_residuals',[sub(a,b) for a,b in combinations(candidates,2)])
    emit('dimension_basis',('Density','Length','Time'))
    emit('dimension_field_signatures',signatures)
    emit('dimension_convention','Exponent triples use abstract Density, Length, Time. Literal zero has no preferred dimension (empty list); all distinct additive term dimensions are retained.')
    emit('chart_domain','0 < theta < pi, 0 <= phi < 2*pi, with angular seams/poles interpreted by chart continuation; radii positive on the global domain.')

    surface_parts=[]
    progress('surface crossings')
    for name,b,orientation in (('inner',rin,-1),('outer',rout,1)):
        progress(name)
        n=orientation*e
        v_surface=chart(b).diff(t)
        v_centre=centres[1].diff(t)
        v_radius=v_surface-v_centre
        area=measure(chart(b),(theta,phi)).refine(S.Q.positive(b))
        for key,obj in [('chart',chart(b)),('outward_normal',n),('surface_velocity',v_surface),('centre_velocity',v_centre),('radius_velocity',v_radius),('area_measure',area)]:
            solo(name+'/'+key,obj)
        pair(name+'/surface_velocity_partition',v_surface,v_centre+v_radius)
        crossings={key:-tensor*n for key,tensor in channels.items()}
        crossings['transport_centre']=G*v_centre.dot(n)
        crossings['transport_radius']=G*v_radius.dot(n)
        components=[{} for i in range(3)]
        for key,vector in crossings.items():
            progress(name+'/'+key)
            for i,raw in enumerate(vector):
                tag=f'{name}/{key}/component_{i+1}'
                obj=surface_integral(strip(raw),b,area)
                solo(tag,obj); dimension_audit(tag,obj)
                components[i][key]=obj
                if key not in ('pressure','other'):
                    buckets=groups(raw)
                    for support,part in buckets.items():
                        group_obj=surface_integral(part,b,area)
                        solo(tag+'/labels_'+names[support],group_obj)
                        dimension_audit(tag+'/labels_'+names[support],group_obj)
                    pair(tag+'/group_partition_integrand',strip(raw),sum(buckets.values()))
                if key=='pressure':
                    pressure_atoms={q for q in raw.atoms(S.Function) if q.func==pressure}
                    expanded=raw.xreplace(dict.fromkeys(pressure_atoms,expansion))
                    expanded_obj=surface_integral(strip(expanded),b,area)
                    solo(tag+'/first_order_choice',expanded_obj)
                    dimension_audit(tag+'/first_order_choice',expanded_obj)
                    emit(tag+'/first_order_domain','p_br defined and differentiable at rho_2+varrho pointwise; retained degree one in rho_1; subtraction also requires whole p_br defined at total density.')
                    remainder_obj=surface_integral(strip(raw-expanded),b,area)
                    solo(tag+'/whole_minus_expansion',remainder_obj)
                    dimension_audit(tag+'/whole_minus_expansion',remainder_obj)
                    pressure_slot=S.Dummy("pressure_slot")
                    pressure_as_polynomial=raw.xreplace(dict.fromkeys(pressure_atoms,pressure_slot))
                    coefficient=S.expand(pressure_as_polynomial).coeff(pressure_slot)
                    for support,part in p_groups.items():
                        group_obj=surface_integral(strip(coefficient)*part,b,area)
                        solo(tag+'/expansion_labels_'+names[support],group_obj)
                        dimension_audit(tag+'/expansion_labels_'+names[support],group_obj)
        surface_parts.append(components)

    progress('region momentum and balance')
    solo('region/chart',chart(radius)); solo('region/volume_measure',volume_measure)
    for i,raw in enumerate(G):
        tag=f'region/component_{i+1}'
        momentum=volume_integral(strip(raw))
        # SymPy differentiates the live bounds and the translated integrand.
        rate=S.diff(momentum,t)
        forcing=volume_integral(force[i])
        for key,obj in [('momentum',momentum),('momentum_rate',rate),('integrated_force',forcing)]:
            solo(tag+'/'+key,obj); dimension_audit(tag+'/'+key,obj)
        buckets=groups(raw)
        for support,part in buckets.items():
            group_momentum=volume_integral(part)
            solo(tag+'/momentum_labels_'+names[support],group_momentum)
            solo(tag+'/rate_labels_'+names[support],S.diff(group_momentum,t))
        pair(tag+'/momentum_group_partition_integrand',strip(raw),sum(buckets.values()))
        parts=[surface_parts[0][i],surface_parts[1][i],forcing]
        rhs=sum(obj for surface in surface_parts for obj in surface[i].values())+forcing
        emit(tag+'/balance_parts',parts)
        pair(tag+'/balance',rate,rhs)
        emit(tag+'/balance/reduction', 'No route in this engine reduces this balance residual with the supplied local balance.')
        emit(tag+'/balance_equation',S.Eq(rate,rhs,evaluate=False))
        dimension_audit(tag+'/balance_rhs',rhs)
        a,b=dim(rate),dim(rhs)
        emit(tag+'/balance_dimension_operand_a',sorted(a))
        emit(tag+'/balance_dimension_operand_b',sorted(b))
        emit(tag+'/balance_dimension_residuals',[sub(da,db) for da in sorted(a) for db in sorted(b)])
    progress('dimensional arguments')
    for key,(a,b) in argument_audit.items():
        emit('dimension_argument/name',key)
        emit('dimension_argument/operand_a',sorted(a))
        emit('dimension_argument/operand_b',sorted(b))
        emit('dimension_argument/residuals',[sub(da,db) for da in sorted(a) for db in sorted(b)])
    emit('balance_construction','Printed objects: the moving-domain momentum integral, its time derivative, the surface parts, the integrated force density, and the balance operands, equation and residual.')
    emit('scope','No force law, coefficient, normalization, mass, throat-force identification, or relation to J, Q, Theta_Q or N_infinity_3 is computed.')
    emit('END_OF_TRANSCRIPT','S15 SymPy objects')


if __name__ == '__main__':
    with OUT.open('w') as transcript:
        main(transcript)
