"""Reciprocal, fixed-slip material test across two sections of a planar conversion front.

This diagnostic lets an independently supplied stream change its order profile.
One trial contact area per specified volume can back-load conversion. Both cut
tractions are measured. No mechanical slip is performed on a traveling contact.
"""
import json
import numpy as np
from scipy.integrate import solve_bvp, simpson
from scipy.special import expit,logit

from mouth_order_wave_pulse import emit


def H(c):
    return c*c*(3-2*c)


def run():
    radius,rate,density,width,kappa,mobility,tau=1.,.25,1.,.04,.08,5.,26.0591
    area=4*np.pi*radius**3/3
    alpha=3/radius  # active wall area per volume in a radius-one straight guide
    speed=rate/(area*density)
    href=speed/(6*mobility*density*width)
    high,low=.85,.45
    length=width*(logit(high)-logit(low))
    sample=np.linspace(0,length,1001)
    mesh=np.linspace(0,length,101)
    start=expit(logit(high)-mesh/width)
    def check(beta,head,b):
        def field(z,c):
            return b*np.exp(beta*(H(high)-H(c))-z/(speed*tau))
        def well(c):
            return .5*kappa*(c*(1-c)/width)**2
        def ode(z,y):
            c,dc=y
            g=field(z,c)
            Hp=6*c*(1-c)
            bare=kappa*c*(1-c)*(1-2*c)/width**2
            reaction=alpha*g**2*(.5-beta*c*Hp)
            return [dc,(speed*dc/mobility+bare+reaction+density*head*Hp)/kappa]
        sol=solve_bvp(ode,lambda l,r:[l[0]-high,r[0]-low],mesh,
                      [start,-start*(1-start)/width],tol=5e-10,max_nodes=9000)
        if sol.status:
            raise ArithmeticError(sol.message)
        c,dc=sol.sol(sample)
        rate_c=speed*dc
        rate_y=-6*c*(1-c)*rate_c
        g=field(sample,c)
        tr=c*g
        def integrate(x):
            return area*simpson(x,x=sample)
        supplied=integrate(density*head*rate_y)
        heat=integrate(rate_c**2/mobility)
        soft=integrate(alpha*g**2*rate_c/2)
        stroke=integrate(alpha*beta*tr*g*rate_y)
        renew=integrate(alpha*tr*g/tau)
        outgoing=area*speed*alpha*(low*g[-1]**2-high*g[0]**2)/2
        stress=well(c)-kappa*dc**2/2
        cut_work=area*speed*(stress[-1]-stress[0])
        force=integrate(alpha*tr)
        out=dict(test='fixed_slip_coupon',beta_input=beta,head=head,prestrain=b,
                 supply=supplied,conversion_heat=heat,order_work=soft,stroke_work=stroke,net=soft+stroke,
                 renewal=renew,increase_in_carried_spring_power=outgoing,flow_cut_stress_work=cut_work,
                 force_on_material=force,force_on_fixture=-force,
                 min_order=c.min(),max_order=c.max(),max_order_slope=dc.max(),
                 contact_closure=outgoing-soft-stroke+renew,supply_closure=supplied-heat-cut_work-soft-stroke)
        if b:
            # Rest lengths and parcel transit give beta without using work, force, or the mouth.
            z=sample[len(sample)//2]
            mid=len(sample)//2
            inferred=(np.log(abs(g[mid]/g[0]))+z/(speed*tau))/(H(high)-H(c[mid]))
            out['beta_inferred_at_middle']=inferred
            out['predicted_outlet_extension']=g[mid]*np.exp(inferred*(H(c[mid])-H(low))-(length-z)/(speed*tau))
            out['measured_outlet_extension']=g[-1]
        emit(**out)
        if (abs(out['supply_closure'])>2e-10 or abs(out['contact_closure'])>2e-12 or
            abs(supplied-rate*head*(H(high)-H(low)))>2e-11 or dc.max()>=0 or
            (b and abs(out['predicted_outlet_extension']-g[-1])>2e-12)):
            raise ArithmeticError('Coupon transfer, monotonicity or held-out stroke failed')
    emit(test='coupon_setup',radius=radius,intake=rate,area=area,area_per_volume=alpha,material_speed=speed,
         width=width,kappa=kappa,mobility=mobility,tau=tau,head_reference=href,ordered_side=high,
         less_ordered_side=low,length=length,net_local_beta_lower_bound=9/16)
    for beta in (0.,.5,1.,2.):
        for b in (0.,.001,-.001,.0005,-.0005):
            for ratio in (.99,1.,1.01):
                check(beta,ratio*href,b)


if __name__=='__main__':
    run()
