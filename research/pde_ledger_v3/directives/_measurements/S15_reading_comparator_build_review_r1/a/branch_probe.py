"""Numerical and physical ground truth for the exceptional intersection in accepted object 019."""
from native_challenges import obj,operands,parsed,m,r,np,mp,s,math,ctx,member,tf
from scipy.special import roots_legendre
import itertools

def reference(exponents,form,point,order):
    # Unit sphere integration uses z=cos(theta) and a periodic azimuth trapezoid.
    t=float(m.POINTS[point][-1]); stratum='generic'
    def at(head):
        return m.jet_function(head,1,(0,),form,stratum,point)(t)
    D=np.array([at(f'X_1_{i}')-at(f'X_2_{i}') for i in (1,2,3)])[:,None,None]
    b=at('r_in')
    z,w=roots_legendre(order);phi=np.arange(2*order)*(math.pi/order)
    z=z[:,None];p=phi[None,:];cp=np.sqrt(1-z*z)
    n=np.stack(np.broadcast_arrays(cp*np.cos(p),cp*np.sin(p),z))
    vec=b*n-D
    rad=np.sqrt(np.sum(vec**2,axis=0))
    a1,a2,b1,b2=map(float,exponents)
    v1=at('A_1')*vec*rad**(a1-1);v2=at('A_2')*b**a2*n
    rho1=at('B_1')*rad**b1;rho2=at('B_2')*b**b2
    K=np.sum(vec*(-D),axis=0)
    dv1=at('A_1')*((a1-1)*rad**(a1-3)*K*vec-rad**(a1-1)*D)
    dr=at('B_1')*b1*rad**(b1-2)*K
    def flux(a,b): return a[1]*np.sum(b*n,axis=0)
    cross=flux(v1,v2)+flux(v2,v1)
    q=(rho1*(cross+flux(v2,v2))+rho2*(cross+flux(v1,v1)))
    dc=flux(dv1,v2)+flux(v2,dv1)
    dq=dr*(cross+flux(v2,v2))+rho1*dc+rho2*(dc+flux(dv1,v1)+flux(v1,dv1))
    # Inner normal is -n, while the rest flux has an extra minus.
    return [float(b**2*math.pi/order*np.sum(v*w[:,None])) for v in (q,dq+2*q)]
def main():
    o=obj(19); fields=('item3/Q_closed_operand','item3/object');ex=(s.Integer(-1),s.Rational(9,4),s.Integer(-2),s.Rational(9,2))
    form=m.FORMS[0];p=0
    allpairs=[parsed(operands(o,name)[1]) for name in fields]
    with mp.workdps(55):
        C=r.PrecisionSample(form,'generic',p,ex)
        print('joint exponents',list(map(str,ex)),'point',list(map(str,m.POINTS[p])),'distance',mp.nstr(C.env[s.Symbol('d')],20))
        for tag,(a,b) in zip(fields,allpairs):
            x,y=C.value(a),C.value(b);print('joint 55 decimal',tag,mp.nstr(x,52),mp.nstr(y,52),mp.nstr(x-y,12))
    for order in (30,60,90):
        print('independent native angular flux and geometric departure',order,reference(ex,form,p,order))
    # Independently realize a 1-unit error active only when native branches meet. Candidate baseline
    # enumerates 14 single hyperplanes plus the generic condition, even for this 14-plane printed object.
    native=operands(o,fields[1])[1][1]
    from native_challenges import m as mm
    mutant=mm.parse('('+native+') + Piecewise[{{1, aa[1] == -1 && bb[1] == -2}}, 0]','wl')
    original=allpairs[1][1]; delta=mutant-original
    mask={}
    for label,expr in [('printed intersection',ex)]+list(mm.exponent_samples(0,[tuple(map(s.Rational,g)) for g in o['exponent_loci']])):
        C=ctx(form,point=0,exp=expr)
        mask[label]=C.value(delta)
    print('native gated form ablation signed synthetic delta WL-mutant minus WL, point 0',mask)
    assert mask['printed intersection']==1 and all(x==0 for name,x in mask.items() if name!='printed intersection')
    old=member(tf,'refinement/019.json.gz')['fields'][0]['rechecks'][0]
    high=r.refine_closed(*allpairs[1],dict(case=old['case'],exponents=old['exponents'],residual=old['old_residual']))
    print('archive flagged native 45/75 check identical',high==old,'status',high['status'],'high residual',high['values_and_residual'][1][-1])
    assert high==old
if __name__=='__main__':main()
