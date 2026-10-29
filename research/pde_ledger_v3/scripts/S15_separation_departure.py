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

OUT = Path(__file__).resolve().parent / 'out' / 'S15_separation_departure.out'


def main(log):
    objects = {}
    def is_object(name):
        parts = name.split('/')
        return ((parts[0] in ('inner','outer') and len(parts)>2 and parts[1] in
                 ('carried','pressure','viscous','other','transport_centre','transport_radius')) or
                (parts[0]=='region' and len(parts)>2 and
                 (parts[2].startswith(('momentum','rate_labels_')) or parts[2]=='integrated_force')))

    def emit(name, obj):
        log.write(f'{name}: {obj}\n')
        log.flush()

    current_object = 'construction'
    def progress(name):
        nonlocal current_object
        current_object = name
        print(f'WORKING_ON: {name}', flush=True)
        emit('WORKING_ON', name)

    def cas_call(step, function, *args, **kwargs):
        step += ' at source line ' + str(sys._getframe(1).f_lineno)
        print(f'WORKING_ON: {current_object} / {step}', flush=True)
        emit('WORKING_ON', current_object + ' / ' + step)
        result = function(*args, **kwargs)
        print(f'FINISHED_STEP: {current_object} / {step}', flush=True)
        return result

    def solo(name, obj):
        if is_object(name): objects[name] = obj
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
    integral_origins={}
    rate_origins={}
    def surface_integral(expr,b,area):
        obj=S.Integral(pull(expr,b)*area,(theta,0,S.pi),(phi,0,2*S.pi))
        integral_origins[obj]=(expr,b,S.cancel(area/S.sin(theta)),False)
        return obj
    def volume_integral(expr):
        obj=S.Integral(pull(expr,radius)*volume_measure,(radius,rin,rout),(theta,0,S.pi),(phi,0,2*S.pi))
        integral_origins[obj]=(expr,radius,S.cancel(volume_measure/S.sin(theta)),True)
        return obj

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
        rate_origins[rate]=momentum
        forcing=volume_integral(force[i])
        for key,obj in [('momentum',momentum),('momentum_rate',rate),('integrated_force',forcing)]:
            solo(tag+'/'+key,obj); dimension_audit(tag+'/'+key,obj)
        buckets=groups(raw)
        for support,part in buckets.items():
            group_momentum=volume_integral(part)
            solo(tag+'/momentum_labels_'+names[support],group_momentum)
            group_rate=S.diff(group_momentum,t)
            rate_origins[group_rate]=group_momentum
            solo(tag+'/rate_labels_'+names[support],group_rate)
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
    # Reading construction: every Q was constructed above from this engine's
    # own tensor. The registry records values, never imports another engine.
    for i in range(3):
        sums = []
        for j, surf in enumerate(('inner','outer')):
            entries = surface_parts[j][i]
            rest = sum(entries[k] for k in ('carried','pressure','viscous','other'))
            transport = sum(entries[k] for k in ('transport_centre','transport_radius'))
            row = dict(rest=rest,transport=transport,crossing=rest+transport)
            sums.append(row)
            for key,value in row.items(): objects[f'{surf}/sum_{key}/component_{i+1}']=value
        for key in sums[0]:
            objects[f'both_surfaces/sum_{key}/component_{i+1}']=sums[0][key]+sums[1][key]

    progress('reading family and derivative operator')
    epsilon = S.Symbol('family_epsilon',real=True)
    dimensions[epsilon]=Z
    changing = [*u_heads,*rho_heads,ds[-1].func,*[q.func for q in vs[-1]],
                eta.func,zeta.func,*[q.func for q in other],*[q.func for q in force]]
    responses={h:S.Function(str(h)+'_change') for h in changing}
    for h,hp in responses.items(): register(hp,*signatures[h])
    # Freeze the instantaneous centre jets during the positional variation.
    # X1-dot has its own live family response; it is not d/dt of epsilon D.
    jets={S.diff(c,t):S.Function(str(c.func)+'_velocity')(t) for cc in centres for c in cc}
    for h in jets.values(): register(h.func,sub(L,T),(T,))
    jet_responses={jets[S.diff(c,t)]:S.Function(str(c.func)+'_velocity_change')(t) for c in centres[0]}
    for h in jet_responses.values(): register(h.func,sub(L,T),(T,))
    back_jets={v:k for k,v in jets.items()}
    shift={a:a+epsilon*(a-b) for a,b in zip(centres[0],centres[1])}
    emit('family/held',(centres[1],centres[1].diff(t),rin,rout,S.diff(rin,t),S.diff(rout,t),pressure))
    emit('family/changing_function_heads',responses)
    emit('family/centre_velocity_changes',jet_responses)
    emit('family/centre_substitution',shift)
    emit('family/jet_notation',back_jets)
    emit('family/choice','Instantaneous family; epsilon dimensionless; D is the unperturbed X_1-X_2. No time differentiation of this family is used.')
    emit('family/response_dimension_signatures',{h:signatures[h] for h in responses.values()})

    @lru_cache(None)
    def variation(expr,with_response):
        if isinstance(expr,S.Integral):
            return S.Integral(variation(expr.function,with_response),*expr.limits)
        if expr.is_Add:
            return S.Add(*(variation(q,with_response) for q in expr.args))
        if expr.is_Mul and expr.has(S.Integral):
            return S.Add(*(variation(f,with_response)*S.Mul(*(g for j,g in enumerate(expr.args) if j!=i)) for i,f in enumerate(expr.args)))
        q=expr.xreplace(jets)
        if with_response:
            substitutions={}
            for h,hp in responses.items():
                aa=tuple(S.Dummy(f'arg_{j}') for j in range(len(signatures[h][1])))
                substitutions[h]=S.Lambda(aa,h(*aa)+epsilon*hp(*aa))
            q=q.replace(lambda node:node.func in substitutions,lambda node:substitutions[node.func](*node.args)).doit()
            q=q.xreplace({h:h+epsilon*hp for h,hp in jet_responses.items()})
        q=q.subs(shift,simultaneous=True)
        return S.diff(q,epsilon).subs(epsilon,0).doit().xreplace(back_jets)

    d=S.Symbol('d',positive=True)
    direction_D=S.symbols('N_1:4',real=True)
    exponents_u=S.symbols('a_1:3',real=True)
    exponents_rho=S.symbols('b_1:3',real=True)
    amplitudes_u=[S.Function(f'A_{k}')(t) for k in (1,2)]
    amplitudes_rho=[S.Function(f'B_{k}')(t) for k in (1,2)]
    radial_arg,time_arg=S.symbols('radial_arg time_arg',positive=True)
    profile_subs={}
    for k in range(2):
        profile_subs[u_heads[k]]=S.Lambda((radial_arg,time_arg),amplitudes_u[k].subs(t,time_arg)*radial_arg**exponents_u[k])
        profile_subs[rho_heads[k]]=S.Lambda((radial_arg,time_arg),amplitudes_rho[k].subs(t,time_arg)*radial_arg**exponents_rho[k])
    d_subs={a:b+d*n for a,b,n in zip(centres[0],centres[1],direction_D)}
    emit('profile_family/substitution',profile_subs)
    emit('profile_family/separation_chart',d_subs)
    emit('profile_family/chart_condition',(d>rout,S.Eq(sum(n*n for n in direction_D),1)))
    emit('profile_family/chart_note','N is an arbitrary real unit vector in the original Cartesian basis; velocities remain independent instantaneous jets.')
    @lru_cache(None)
    def on_family(expr):
        # Apply profile definitions before the separation chart, protecting the
        # independent instantaneous velocities from differentiating that chart.
        q=expr.xreplace(jets)
        q=q.replace(lambda node:node.func in profile_subs,lambda node:profile_subs[node.func](*node.args))
        # Evaluate only derivative/substitution nodes, never the integrals here.
        q=q.replace(lambda n:isinstance(n,(S.Derivative,S.Subs)),lambda n:n.doit())
        return q.subs(d_subs,simultaneous=True).xreplace(back_jets)

    allowed_heads={q.func for cc in centres for q in cc} | {rin.func,rout.func} | {q.func for q in amplitudes_u+amplitudes_rho}
    def has_general_field(expr):
        from sympy.core.function import AppliedUndef
        return any(q.func not in allowed_heads for q in expr.atoms(AppliedUndef))

    # Item 3 uses the pre-chart integrands already constructed from Pi and G.
    # The following kernels contain mathematical dummy variables only.
    mu, az = S.symbols('mu az', real=True)
    axis = S.symbols('axis_1:4', real=True)
    nv = S.symbols('normal_1:4', real=True)
    hh = S.Symbol('profile_distance', positive=True)
    rr = S.Symbol('integration_radius', positive=True)
    dd = S.Symbol('kernel_distance', positive=True)
    ww = S.Symbol('squared_distance', positive=True)
    pp, qq = S.symbols('kernel_power radial_power', real=True)
    low, high = S.symbols('kernel_low kernel_high', positive=True)
    sphere_gb = S.groebner([sum(z*z for z in axis)-1],*axis,domain='EX')
    def unit_reduce(expr):
        numerator,denominator=S.fraction(S.cancel(expr))
        return S.cancel(sphere_gb.reduce(numerator)[1]/sphere_gb.reduce(denominator)[1])
    pole=S.Matrix([0,0,1]); av=S.Matrix(axis)
    hv=pole-av
    rotation=S.eye(3)-hv*hv.T/(1-axis[2])
    rotated=rotation*S.Matrix([S.sqrt(1-mu**2)*S.cos(az),S.sqrt(1-mu**2)*S.sin(az),mu])
    pair('item3_kernel/rotation_metric',rotation.T*rotation,S.eye(3))
    emit('item3_kernel/rotation_metric_unit_reduction',(rotation.T*rotation-S.eye(3)).applyfunc(unit_reduce))
    emit('item3_kernel/rotation_chart_condition','axis is real unit; intermediate chart axis_3 != 1; polynomial moments extend to that pole by continuity.')
    @lru_cache(None)
    def azimuth_moment(powers):
        operand=S.expand(S.prod(v**n for v,n in zip(rotated,powers)))
        answer=0
        for term in S.Add.make_args(operand):
            coefficient,angular=term.as_independent(az,as_Add=False)
            answer+=coefficient*cas_call('integrate', S.integrate, angular,(az,0,2*S.pi))
        answer=S.expand(unit_reduce(answer))
        solo('item3_kernel/azimuth_operand',S.Integral(operand,(az,0,2*S.pi)))
        solo('item3_kernel/azimuth_object',answer)
        return answer
    # The squared distance and Jacobian come from the rotated chart itself.
    squared=S.trigsimp(S.expand((rr*rotated-dd*av).dot(rr*rotated-dd*av)))
    squared=S.trigsimp(unit_reduce(squared))
    solo('item3_kernel/distance_equation_operand',S.Eq(ww,squared))
    distance_solutions=cas_call('solve', S.solve, S.Eq(ww,squared),mu)
    solo('item3_kernel/distance_equation_solutions',distance_solutions)
    mu_inverse=distance_solutions[0]
    pair('item3_kernel/distance_change_relation',squared.subs(mu,mu_inverse),ww)
    jacobian=-S.diff(mu_inverse,ww)
    angular_bounds=[S.factor(squared.subs(mu,z)) for z in (1,-1)]
    endpoint_distance=S.Symbol('endpoint_distance',positive=True)
    endpoint_assumptions=S.Q.positive(dd-rr)&S.Q.negative(rr-dd)&S.Q.positive(dd+rr)
    angular_distance_bounds=[S.refine(S.sqrt(bound),endpoint_assumptions) for bound in angular_bounds]
    solo('item3_kernel/endpoint_domain',dd>rr)
    solo('item3_kernel/distance_endpoints',angular_distance_bounds)
    for bound,root in zip(angular_bounds,angular_distance_bounds):
        pair('item3_kernel/distance_endpoint_relation',root**2,bound)
    solo('item3_kernel/squared_distance',squared)
    solo('item3_kernel/change_of_variable',(mu_inverse,jacobian,angular_bounds))
    primitive_kernel=cas_call('integrate', S.integrate, ww**pp,ww,conds='piecewise')
    solo('item3_kernel/power_primitive_operand',ww**pp)
    solo('item3_kernel/power_primitive',primitive_kernel)
    pair('item3_kernel/power_primitive_derivative',S.diff(primitive_kernel,ww),ww**pp)

    # A finite-interval Appell representation avoids singular lower parameters
    # at negative integer exponents. The coefficient identity derives the
    # integration kernel; it is not a physical ansatz or a truncated expansion.
    m,n=S.symbols('kernel_m kernel_n',integer=True,nonnegative=True)
    aa,bb,xx,yy,zz=S.symbols('kernel_a kernel_b kernel_x kernel_y kernel_z')
    binomial_coefficient=S.rf(-aa,m)*S.rf(-bb,n)/(S.factorial(m)*S.factorial(n))
    integrated_coefficient=binomial_coefficient*cas_call('integrate', S.integrate, zz**(m+n),(zz,0,1))
    appell_coefficient=binomial_coefficient*S.rf(1,m+n)/S.rf(2,m+n)
    pair('item3_kernel/appell_coefficients',integrated_coefficient,appell_coefficient)
    solo('item3_kernel/appell_coefficient_residual',S.simplify(integrated_coefficient-appell_coefficient))
    emit('item3_kernel/appell_continuation','The coefficient identity defines the analytic germ of the integrated binomial product. Analytic continuation along the real integration interval is used only where both affine bases are positive; Appell lower parameter 2 is fixed, so all real powers are covered.')

    @lru_cache(None)
    def affine_integral(power_r, power_affine, sign):
        print(f'WORKING_ON: {current_object} / finite-interval affine-power integration',flush=True)
        emit('WORKING_ON',current_object+' / finite-interval affine-power integration')
        # Substitute the endpoints into an affine unit-interval chart, derive
        # its two binomial parameters and prefactor by evaluation/differentiation.
        mapping=low+(high-low)*zz
        bases=[mapping,dd+sign*mapping]
        constants=[z.subs(zz,0) for z in bases]
        parameters=[-S.diff(z,zz)/c for z,c in zip(bases,constants)]
        return (S.diff(mapping,zz)*constants[0]**power_r*constants[1]**power_affine*
                S.appellf1(1,-power_r,-power_affine,2,*parameters))

    # Logarithmic angular branches are integrated by parts. The remaining
    # rational-power integral uses the same nonsingular Appell kernel.
    log_primitives={}
    vv=S.Symbol('log_argument',positive=True)
    for sign in (-1,1):
        lp=cas_call('integrate', S.integrate, S.log(1+sign*vv)/vv,vv)
        # On 0 < vv < 1 these polar arguments are the principal real branches.
        lp=lp.replace(S.exp_polar,lambda arg:S.exp(arg))
        log_primitives[sign]=lp
        solo('item3_kernel/log_primitive_operand',S.log(1+sign*vv)/vv)
        solo('item3_kernel/log_primitive',lp)
        pair('item3_kernel/log_primitive_derivative',S.diff(lp,vv),S.log(1+sign*vv)/vv)
    @lru_cache(None)
    def radial_log(power_r,sign):
        print(f'WORKING_ON: {current_object} / radial logarithm integration by parts',flush=True)
        emit('WORKING_ON',current_object+' / radial logarithm integration by parts')
        primitive=cas_call('integrate', S.integrate, rr**qq,rr,conds='none').subs(qq,power_r)
        boundary=(primitive*S.log(dd+sign*rr)).subs(rr,high)-(primitive*S.log(dd+sign*rr)).subs(rr,low)
        regular=boundary-sign*affine_integral(power_r+1,-1,sign)/(power_r+1)
        special=S.log(dd)*S.log(high/low)+log_primitives[sign].subs(vv,high/dd)-log_primitives[sign].subs(vv,low/dd)
        return S.Piecewise((regular,S.Ne(power_r,-1)),(special,True))

    @lru_cache(None)
    def radial_integral(expr):
        print(f'WORKING_ON: {current_object} / radial term integration',flush=True)
        emit('WORKING_ON',current_object+' / radial term integration')
        # At most one endpoint power/log occurs in an angular monomial term.
        expr=S.expand(expr,power_base=False,power_exp=False,log=False)
        answer=0
        for term in S.Add.make_args(expr):
            term=S.factor_terms(term)
            numerator,denominator=S.fraction(term)
            numerator_powers=numerator.as_powers_dict()
            denominator_powers=denominator.as_powers_dict()
            power_r=(numerator_powers.pop(rr,S.Integer(0))-
                     denominator_powers.pop(rr,S.Integer(0)))
            # Rebuild the remaining product without expanding an affine power
            # into a polynomial times another symbolic power.
            rest=(S.prod(base**exponent for base,exponent in numerator_powers.items())/
                  S.prod(base**exponent for base,exponent in denominator_powers.items()))
            logs=list(rest.atoms(S.log))
            if logs:
                logarithm=logs[0]
                sign=S.diff(logarithm.args[0],rr)
                coefficient=S.cancel(rest/logarithm)
                pair('item3_kernel/radial_log_factorization',term,coefficient*rr**power_r*logarithm)
                emit('item3_kernel/radial_log_parts',(power_r,coefficient,logarithm,sign))
                if coefficient.has(rr) or sign not in (-1,1):
                    raise ValueError(('radial log factorization',term))
                answer+=coefficient*radial_log(power_r,sign)
            else:
                affine_powers=[z for z in S.Mul.make_args(rest) if z.has(rr)]
                if affine_powers:
                    factor=affine_powers[0]
                    base,exponent=factor.as_base_exp()
                    sign=S.diff(base,rr)
                    coefficient=S.cancel(rest/factor)
                    pair('item3_kernel/radial_affine_factorization',term,coefficient*rr**power_r*base**exponent)
                    emit('item3_kernel/radial_affine_parts',(power_r,coefficient,base,exponent,sign))
                    if coefficient.has(rr) or sign not in (-1,1) or S.expand(base-sign*rr)!=dd:
                        raise ValueError(('radial affine factorization',term))
                    answer+=coefficient*affine_integral(power_r,exponent,sign)
                else:
                    integral=cas_call('integrate', S.integrate, rr**qq,(rr,low,high),conds='piecewise').subs(qq,power_r)
                    answer+=rest*integral
        return answer

    @lru_cache(None)
    def angular_kernel(powers,power,volume):
        print(f'WORKING_ON: {current_object} / angular kernel integration',flush=True)
        emit('WORKING_ON',current_object+' / angular kernel integration')
        moment=azimuth_moment(powers)
        if power.is_number:
            value=cas_call('integrate', S.integrate, moment*squared**(power/2),(mu,-1,1))
            return radial_integral(S.expand(value*rr**qq)) if volume else value
        polynomial=S.Poly(S.expand(moment.subs(mu,mu_inverse)*jacobian),ww)
        result=0
        for (degree,),coefficient in polynomial.terms():
            exponent=pp/2+degree
            # Compute all exceptional powers from the primitive's conditions.
            integral=primitive_kernel.subs(pp,exponent)
            branches=[]
            for primitive,condition in integral.args if isinstance(integral,S.Piecewise) else [(integral,S.true)]:
                bounds=[]
                # Simplify through a positive dummy distance first; only then
                # substitute the positive roots derived on d > r > 0.
                endpoint=S.expand_log(S.powdenest(primitive.subs(ww,endpoint_distance**2),force=True),force=True)
                for root in angular_distance_bounds:
                    bounds.append(endpoint.subs(endpoint_distance,root))
                value=S.expand(coefficient*(bounds[1]-bounds[0]),power_base=False,power_exp=False,log=False)
                if volume:
                    value=radial_integral(S.expand(value*rr**qq,power_base=False,power_exp=False,log=False))
                branches.append((value,condition))
            result+=S.Piecewise(*branches).subs(pp,power)
        return S.simplify(result) if power.is_number else result

    trig_to_normal={S.cos(theta):nv[2],S.sin(theta):S.Symbol('chart_sine'),
                    S.cos(phi):nv[0]/S.Symbol('chart_sine'),S.sin(phi):nv[1]/S.Symbol('chart_sine')}
    @lru_cache(None)
    def integrated_origin(obj):
        print(f'WORKING_ON: {current_object} / profile integral construction',flush=True)
        emit('WORKING_ON',current_object+' / profile integral construction')
        raw,b,measure_without_sine,is_volume=integral_origins[obj]
        # Lift positive radial distances before substituting power profiles.
        distance_bases={S.expand(distances[0]**2):hh**2,S.expand(distances[1]**2):rr**2}
        lifted=raw.replace(lambda z:z.is_Pow and z.base.has(*xyz) and S.expand(z.base) in distance_bases,
                           lambda z:distance_bases[S.expand(z.base)]**z.exp)
        lifted=lifted.replace(lambda z:z.func in profile_subs,lambda z:profile_subs[z.func](*z.args))
        lifted=lifted.xreplace(jets).subs(dict(zip(xyz,centres[1]+rr*S.Matrix(nv))),simultaneous=True)
        lifted=lifted.subs(d_subs,simultaneous=True).xreplace(trig_to_normal)
        lifted=S.expand(S.powdenest(lifted,force=True))
        measure_generic=S.simplify(measure_without_sine.subs(b,rr)) if not is_volume else measure_without_sine.subs(radius,rr)
        lifted=S.expand(lifted*measure_generic)
        polynomial=S.Poly(lifted,*nv)
        answer=0
        for powers,coefficient in polynomial.terms():
            for term in S.Add.make_args(S.expand(coefficient)):
                term=S.powsimp(term,force=True)
                p=term.as_powers_dict().get(hh,S.Integer(0))
                radial_power=term.as_powers_dict().get(rr,S.Integer(0))
                prefactor=S.cancel(term/(hh**p*rr**radial_power))
                pair('item3_kernel/profile_factorization',term,prefactor*hh**p*rr**radial_power)
                emit('item3_kernel/profile_parts',(powers,p,radial_power,prefactor))
                if prefactor.has(hh,rr,theta,phi):
                    raise ValueError(('profile monomial factorization',prefactor))
                kernel=angular_kernel(powers,p,is_volume)
                substitutions={dd:d,**dict(zip(axis,direction_D))}
                if is_volume: substitutions.update({qq:radial_power,low:rin,high:rout})
                else: substitutions[rr]=b
                answer+=prefactor*(1 if is_volume else b**radial_power)*kernel.subs(substitutions,simultaneous=True)
        return answer.xreplace(back_jets)

    # Differentiate the computed closed momentum with the chain rule for the
    # time-dependent separation chart. The velocity coefficients are derived
    # from X1(t), X2(t); they are not assigned by hand.
    position_delta=centres[0]-centres[1]
    distance_time=S.diff(separation,t)
    direction_time=(position_delta/separation).diff(t)
    radial_time=distance_time.xreplace(jets).subs(d_subs,simultaneous=True)
    normal_time=direction_time.xreplace(jets).subs(d_subs,simultaneous=True)
    # Reduce the unit-direction norm before restoring velocity jets.
    norm=sum(z*z for z in direction_D)
    radial_time=S.simplify(radial_time).subs(norm,1).xreplace(back_jets)
    normal_time=normal_time.applyfunc(lambda z:S.simplify(z).subs(norm,1)).xreplace(back_jets)
    solo('item3_kernel/separation_time_derivative',radial_time)
    solo('item3_kernel/direction_time_derivative',normal_time)
    geometry_distance=S.simplify(on_family(variation(separation,False))).subs(norm,1)
    geometry_direction=[S.simplify(on_family(variation(z,False))).subs(norm,1) for z in position_delta/separation]
    solo('item3_kernel/geometric_distance_derivative',geometry_distance)
    solo('item3_kernel/geometric_direction_derivative',geometry_direction)
    def closed_geometric_derivative(expr):
        print(f'WORKING_ON: {current_object} / closed geometric differentiation',flush=True)
        emit('WORKING_ON',current_object+' / closed geometric differentiation')
        return geometry_distance*S.diff(expr,d)+sum(v*S.diff(expr,n) for n,v in zip(direction_D,geometry_direction))
    @lru_cache(None)
    def closed_q(obj):
        if obj in rate_origins:
            momentum_closed=closed_q(rate_origins[obj])
            return S.diff(momentum_closed,t)+S.diff(momentum_closed,d)*radial_time+sum(S.diff(momentum_closed,n)*v for n,v in zip(direction_D,normal_time))
        if isinstance(obj,S.Integral): return integrated_origin(obj)
        if obj.is_Add: return S.Add(*(closed_q(z) for z in obj.args))
        if obj.is_Mul: return S.Mul(*(closed_q(z) for z in obj.args))
        return on_family(obj)

    # Direction coefficients are reduced on the real unit sphere before
    # attempting the d identity. Piecewise branches involve exponents only;
    # distributing coefficients, not whole branches, avoids losing strata.
    direction_zero=(0,0,0)
    direction_basis=S.groebner([sum(n*n for n in direction_D)-1],*direction_D,domain='EX')
    @lru_cache(None)
    def direction_coefficients(expr):
        if not expr.has(*direction_D): return {direction_zero:expr}
        if isinstance(expr,S.Piecewise):
            rows=[(direction_coefficients(value),condition) for value,condition in expr.args]
            monomials=set().union(*(set(row) for row,condition in rows))
            return {powers:S.Piecewise(*[(row.get(powers,0),condition) for row,condition in rows]) for powers in monomials}
        if expr.is_Add:
            result={}
            for value in expr.args:
                for powers,coefficient in direction_coefficients(value).items():
                    result[powers]=result.get(powers,0)+coefficient
            return result
        if expr.is_Mul:
            result={direction_zero:S.Integer(1)}
            for value in expr.args:
                following={}
                for left,a in result.items():
                    for right,b in direction_coefficients(value).items():
                        powers=tuple(x+y for x,y in zip(left,right))
                        following[powers]=following.get(powers,0)+a*b
                result=following
            return result
        return dict(S.Poly(expr,*direction_D,domain='EX').terms())

    def abbreviate_identity(coefficient,parameters):
        # Exact definitions compress analytic kernels, retaining every free
        # argument in the matching problem. They are not independent fields.
        # The solver may fail to classify their identities; no analytic
        # independence or complete solution is inferred from such a return.
        definitions={};names={}
        def descend(value):
            if (isinstance(value,S.Piecewise) or value.func in (S.appellf1,S.polylog)
                    or (value.is_Pow and value.exp.has(*parameters))):
                if value not in names:
                    arguments=tuple(sorted(value.free_symbols,key=str))
                    abbreviation=S.Function('condition_kernel_'+str(len(names)))(*arguments)
                    names[value]=abbreviation
                    definitions[abbreviation]=value
                return names[value]
            if not value.args:return value
            return value.func(*(descend(argument) for argument in value.args))
        abbreviated=descend(coefficient)
        return abbreviated,definitions

    def solve_conditions(name,result,qclosed):
        tag='reading/'+name+'/item4'
        # Parameters are instant values/jets. For a polynomial in D direction,
        # reduction on the unit sphere supplies necessary coefficient equations.
        expr=result.xreplace(jets)
        qjet=qclosed.xreplace(jets)
        params=sorted(expr.atoms(S.Derivative)|expr.atoms(S.core.function.AppliedUndef)|
                      qjet.atoms(S.Derivative)|qjet.atoms(S.core.function.AppliedUndef)|{rin,rout},key=str)
        params+=list(exponents_u+exponents_rho)
        slots=tuple(S.Symbol('condition_parameter_'+str(i),positive=True) if p in (rin,rout)
                    else S.Symbol('condition_parameter_'+str(i),real=True) for i,p in enumerate(params))
        forward=dict(zip(params,slots)); backward=dict(zip(slots,params))
        lifted=expr.xreplace(forward)
        emit(tag+'/parameter_slots',backward)
        emit(tag+'/condition_domain',(rout>rin,rin>0,*[S.Contains(a,S.S.Reals) for a in (*exponents_u,*exponents_rho)]))
        emit(tag+'/identity_operand',lifted)
        coefficients=direction_coefficients(lifted)
        polynomial=sum(coefficient*S.prod(n**p for n,p in zip(direction_D,powers)) for powers,coefficient in coefficients.items())
        emit(tag+'/direction_reduction_operand',polynomial)
        print(f'WORKING_ON: {current_object} / direction polynomial reduction',flush=True)
        emit('WORKING_ON',current_object+' / direction polynomial reduction')
        # Normal-form reduction is linear over the scalar coefficient field.
        # Reduce the monomials in the fixed unit-sphere ideal, then recombine
        # the untouched coefficient objects (including every Piecewise branch).
        quotient=[S.Integer(0) for _ in direction_basis.polys]
        reduced_coefficients={}
        for powers,coefficient in coefficients.items():
            monomial=S.prod(n**p for n,p in zip(direction_D,powers))
            emit(tag+'/monomial_reduction_operand',monomial)
            print(f'WORKING_ON: {current_object} / monomial reduction {powers}',flush=True)
            emit('WORKING_ON',current_object+' / monomial reduction '+str(powers))
            local_quotient,local_remainder=direction_basis.reduce(monomial)
            emit(tag+'/monomial_reduction_return',(local_quotient,local_remainder))
            quotient=[a+coefficient*b for a,b in zip(quotient,local_quotient)]
            for remainder_power,number in S.Poly(local_remainder,*direction_D).terms():
                reduced_coefficients[remainder_power]=reduced_coefficients.get(remainder_power,0)+coefficient*number
        reduced=sum(coefficient*S.prod(n**p for n,p in zip(direction_D,powers))
                    for powers,coefficient in reduced_coefficients.items())
        emit(tag+'/direction_polynomial',polynomial)
        emit(tag+'/unit_sphere_quotient',quotient)
        emit(tag+'/unit_sphere_remainder',reduced)
        print(f'WORKING_ON: {current_object} / unit-sphere reconstruction',flush=True)
        emit('WORKING_ON',current_object+' / unit-sphere reconstruction')
        pair(tag+'/unit_sphere_reconstruction',polynomial,reduced+sum(a*b.as_expr() for a,b in zip(quotient,direction_basis.polys)))
        emit(tag+'/direction_coefficients',list(reduced_coefficients.items()))
        emit(tag+'/identity_domain',(d>rout.xreplace(forward),S.Eq(sum(z*z for z in direction_D),1)))
        attempted=[]
        for powers,coefficient in reduced_coefficients.items():
            emit(tag+'/scalar_identity_operand',(powers,coefficient))
            print(f'WORKING_ON: {current_object} / kernel abbreviation {powers}',flush=True)
            emit('WORKING_ON',current_object+' / kernel abbreviation '+str(powers))
            solver_operand,definitions=abbreviate_identity(coefficient,slots) if coefficient.has(d) else (coefficient,{})
            emit(tag+'/kernel_definitions',definitions)
            emit(tag+'/scalar_solver_operand',(powers,solver_operand))
            pair(tag+'/kernel_definition_reconstruction',coefficient,solver_operand.xreplace(definitions))
            emit(tag+'/kernel_definition_scope','Exact abbreviations with every free argument retained. Coefficient matching does not establish analytic independence of these kernels or classify their exponent strata; an unresolved identity remains a gap.')
            try:
                result_for_coefficient=cas_call('solve_undetermined_coeffs '+str(powers),S.solve_undetermined_coeffs,solver_operand,slots,d,dict=True)
            except (NotImplementedError,S.PolynomialError) as exc:
                result_for_coefficient=type(exc).__name__+': '+str(exc)
            attempted.append((powers,result_for_coefficient))
            emit(tag+'/scalar_identity_solver_return',(powers,result_for_coefficient))
        emit(tag+'/identity_solver_return',attempted)
        complete=(not reduced.has(d,*direction_D,S.Piecewise,S.appellf1,S.polylog))
        if complete:
            # Factor zero sets avoid division by any possibly zero amplitude
            # or velocity. Positive radius powers are nonzero on the domain.
            factored=S.factor_terms(reduced)
            emit(tag+'/factor_operand',factored)
            pair(tag+'/factorization',reduced,factored)
            solutions=[]
            if factored==0:
                solutions=[{}]
                emit(tag+'/factor_zero_set',S.S.UniversalSet)
            else:
                for factor in S.Mul.make_args(factored):
                    if factor.is_nonzero is True:
                        emit(tag+'/factor_zero_set',(factor,S.EmptySet))
                        continue
                    variables=sorted(factor.free_symbols & set(slots),key=str)
                    if len(variables)!=1:
                        emit(tag+'/factor_zero_set',(factor,'Unresolved multivariate factor after the identity-solver attempt.'))
                        complete=False
                        continue
                    variable=variables[0]
                    try:
                        roots=cas_call('solveset', S.solveset, factor,variable,domain=S.S.Reals)
                    except (NotImplementedError,ValueError,TypeError) as exc:
                        roots=S.ConditionSet(variable,S.Eq(factor,0),S.S.Reals)
                        emit(tag+'/factor_solver_exception',str(exc))
                    emit(tag+'/factor_zero_set',(factor,variable,roots))
                    if roots is S.EmptySet: continue
                    if isinstance(roots,S.FiniteSet):
                        solutions.extend({variable:root} for root in roots)
                    else: complete=False
        else: solutions=[]
        if complete:
            emit(tag+'/condition_union',S.Or(*(S.And(*(S.Eq(k.xreplace(backward),v.xreplace(backward),evaluate=False)
                 for k,v in solution.items())) for solution in solutions)).xreplace(back_jets))
        for branch,solution in enumerate(solutions):
            restrictions={k.xreplace(backward):v.xreplace(backward) for k,v in solution.items()}
            condition=S.And(rout>rin,rin>0,*(S.Eq(k,v,evaluate=False) for k,v in restrictions.items()))
            emit(tag+'/condition_'+str(branch),condition.xreplace(back_jets))
            # Restrict independent instantaneous slots. A(t)=0 at this instant
            # must not silently replace its independent time derivative by zero.
            operand=qjet.xreplace(forward)
            substituted=operand.subs(solution,simultaneous=True)
            solo(tag+'/condition_residual_'+str(branch),S.simplify(reduced.subs(solution,simultaneous=True)).xreplace(backward).xreplace(back_jets))
            restricted=S.simplify(substituted).xreplace(backward).xreplace(back_jets)
            solo(tag+'/Q_under_condition_'+str(branch),restricted)
            pair(tag+'/restriction_substitution_'+str(branch),operand,substituted)
        solo(tag+'/Q_on_family',qclosed)
        if complete:
            emit(tag+'/status','Computed factor zero-set conditions and Q after their substitutions are printed; the separation variables do not occur in the unit-sphere remainder.')
        else:
            reason='NOT OBTAINED EXACTLY: Unit-sphere monomial reduction, exact kernel abbreviation retaining all arguments, and solve_undetermined_coeffs on each scalar coefficient were attempted. Kernel identities and exceptional exponent strata remain unclassified; complete conditions and all corresponding restricted Q objects were not obtained.'
            emit(tag+'/status',reason)
            gaps.append((name,'item4',reason))

    existence=('SUFFICIENT existence condition (not necessary or complete): on an open neighborhood of the compact integration shell/surfaces, avoiding both cores, '
               'all field profiles and responses are C^4 in their displayed arguments, centres and radii are C^2 in time, and the epsilon-family admits these derivatives '
               'uniformly with an integrable dominating function for the differentiated integrands. p_br is C^3 on open intervals containing every total/base density used. '
               'The original S15 domain and pointwise pressure conditions also hold. These hypotheses permit differentiation under each integral and mixed partial interchange.')
    gaps=[]
    emit('Q_manifest',list(objects))
    for name,q in objects.items():
        progress('reading/'+name)
        tag='reading/'+name
        solo(tag+'/Q',q)
        emit(tag+'/departure_existence_condition',existence)
        geo_derivative=variation(q,False)
        full_derivative=variation(q,True)
        geo=geo_derivative+2*q
        full=full_derivative+2*q
        response=full-geo
        solo(tag+'/separation_derivative',full_derivative)
        solo(tag+'/E_geo',geo)
        solo(tag+'/E_resp',response)
        pair(tag+'/departure_partition',full,geo+response)
        for objname,obj in (('E_geo',geo),('E_resp',response)):
            dimension_audit(tag+'/'+objname,obj)
            emit(tag+'/'+objname+'/dimension_reference_operand',sorted(dim(q)))
            emit(tag+'/'+objname+'/dimension_reference_residuals',[sub(a,b) for a in dim(obj) for b in dim(q)])
        for k in range(2):
            replacement={S.diff(c,t):S.Integer(0) for c in centres[k]}
            without=q.xreplace(replacement)
            emit(tag+f'/velocity_part_{k+1}/operand_a',q)
            emit(tag+f'/velocity_part_{k+1}/operand_b',without)
            part=q-without
            solo(tag+f'/velocity_part_{k+1}/residual',part)
            dimension_audit(tag+f'/velocity_part_{k+1}',part)
            emit(tag+f'/velocity_part_{k+1}/dimension_reference_residuals',[sub(a,b) for a in dim(part) for b in dim(q)])
        family_geo=on_family(geo)
        family_q=on_family(q)
        solo(tag+'/item3/family_integral',family_geo)
        emit(tag+'/item3/existence_condition',existence)
        if has_general_field(family_geo):
            result=family_geo
            status='Exact family object retained unevaluated: integrand contains a general field, as permitted by item 3.'
            qclosed=family_q
        else:
            qclosed=closed_q(q)
            solo(tag+'/item3/Q_closed_operand',qclosed)
            result=closed_geometric_derivative(qclosed)+2*qclosed
            status='Closed form from exact azimuth moments, squared-distance integration and finite-interval Appell integration, with logarithmic exponent branches.'
        solo(tag+'/item3/object',result)
        if not has_general_field(family_geo) and result.has(S.Integral,S.nan,S.zoo,S.oo,-S.oo):
            status='NOT OBTAINED IN CLOSED FORM: exact azimuth, squared-distance and Appell-kernel integration was attempted; an Integral or undefined/nonfinite expression remains in its return.'
            gaps.append((name,'item3',status))
        emit(tag+'/item3/status',status)
        if result.has(S.Integral,S.nan,S.zoo,S.oo,-S.oo) or has_general_field(result):
            emit(tag+'/item4/status','Not applicable: the printed item-3 form retains a general field or is not a usable closed form.')
            solo(tag+'/item4/Q_on_family',family_q)
        else:
            solve_conditions(name,result,qclosed)
    for name,item,reason in gaps: emit('EXACT_OBJECT_GAP',(name,item,reason))
    emit('reading/dimension_argument_audit',argument_audit)
    emit('END_OF_TRANSCRIPT','S15 separation departure SymPy objects')



if __name__ == '__main__':
    with OUT.open('w') as transcript:
        main(transcript)
