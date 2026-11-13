"""Ground native velocity-part values and tensor response without suppressing other native variations."""
from native_challenges import *
from scipy.special import roots_legendre

def velocity_grounding():
    # V1=x-X1, rho1=|x-X1|^2; the native drain1 pure group already isolates reference.
    u=m.original.FORMAL; instant=m.POINTS[0][-1];tau=u[0]-instant
    ordinary=m.jet_function
    @lru_cache(None)
    def jets(name,arity,orders,form,stratum,point):
        if name not in ('r_in','R','u_1','rho_1') and not m.geometry_head(name):
            return ordinary(name,arity,orders,form,stratum,point)
        src_stop=stratum=='X1_velocity_zero'
        tar_stop=stratum in ('X2_velocity_zero','X2_and_radii_velocity_zero')
        rad_stop=stratum in ('radii_velocity_zero','X2_and_radii_velocity_zero')
        if name=='r_in':value=n.inner+tau*(0 if rad_stop else n.radial_rates[0])+tau**2
        elif name=='R':value=n.outer+tau*(0 if rad_stop else n.radial_rates[1])+2*tau**2
        elif name=='u_1':value=u[0]
        elif name=='rho_1':value=u[0]**2
        else:
            _,k,j=name.split('_');i=int(j)-1
            pos,v,stop=(n.source,n.v1,src_stop) if k=='1' else (n.target,n.v2,tar_stop)
            value=pos[i]+tau*(0 if stop else v[i])+(2*i+int(k))*tau**2
        ds=[(v,k) for v,k in zip(u[:arity],orders) if k]
        return s.lambdify(u[:arity],s.diff(value,*ds) if ds else value,'numpy')
    o=getobj('region/component_1/rate_labels_drain1')
    raw={q:parsed(operands(o,q)[1]) for q in ('Q','velocity_part_1/residual','velocity_part_2/residual','velocity_part_1/operand_b','velocity_part_2/operand_b')}
    with patch.object(m,'jet_function',jets):
        for stratum,reference in zip(m.STRATA,[
            ('-68009/60','-4021/6','-1837/12'),('-4902/5','-4021/6','0'),('-3293/4','-4021/6','-1837/12'),
            ('-4021/6','-4021/6','0'),('-27799/60','0','-1837/12')]):
            C=ctx(stratum=stratum)
            vals={tag:[float(C.value(q))/math.pi for q in pair] for tag,pair in raw.items()}
            total,w1,w2=map(s.Rational,reference);exact=[total,w1,w2,total-w1,total-w2]
            err=max(abs(x-float(r)) for r,ps in zip(exact,vals.values()) for x in ps)
            print('native group rate, velocity residues, each zeroed operand /pi',stratum,vals,'max absolute error',err)
            assert err<1.e-9

def response_grounding():
    # Add an independent, nonsymmetric tensor-valued variation, leaving *every other native term* as in the archive.
    u=m.original.FORMAL
    plus={}
    for i,j in itertools.product(range(3),repeat=2):
        code=11+3*i+j;value=code+sum(code**(d+1)*v for d,v in enumerate(u[:3]))
        plus[f'Pi_o_{i+1}_{j+1}_change']=value
    old=m.jet_function
    @lru_cache(None)
    def aug(name,arity,orders,form,stratum,point):
        base=old(name,arity,orders,form,stratum,point)
        if name not in plus:return base
        ds=[(v,k) for v,k in zip(u[:arity],orders) if k]
        extra=s.lambdify(u[:arity],s.diff(plus[name],*ds) if ds else plus[name],'numpy')
        return lambda *xs:base(*xs)+extra(*xs)
    for j in range(3):
        o=getobj(f'both_surfaces/sum_crossing/component_{j+1}')
        ref,raw=operands(o,'E_resp'); expr=parsed(raw)
        C=ctx(); base=[C.value(v) for v in expr]
        rg=float(old('R',1,(0,),m.FORMS[0],m.STRATA[0],0)(C.t))
        ri=float(old('r_in',1,(0,),m.FORMS[0],m.STRATA[0],0)(C.t))
        exact=-4*math.pi*(rg**3-ri**3)/3*sum((11+3*j+k)**(k+1) for k in range(3))
        with patch.object(m,'jet_function',aug):
            changed=ctx(); new=[changed.value(v) for v in expr]
        delta=[float(y-x) for x,y in zip(base,new)]
        print('change to archived E_resp by independently specified tensor; all else live',j+1,ref,delta,'divergence theorem',exact,'errors',[x-exact for x in delta])
        assert max(abs(x-exact) for x in delta)<1.e-7
if __name__=='__main__':
    velocity_grounding()
    response_grounding()
