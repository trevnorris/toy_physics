"""Bend the mouth only after launching its meridional wave, then hold both ends.

The kinetic law follows material displacement in the translated sections.
q denotes the pulled-back one-form. In addition to M(b), its equation
contains -2 bdot (N.T-b M2) qdot and -bddot (N.T-b M2) q. The restoring
reaction is reported positive toward increasing flare displacement.
"""
import argparse

import numpy as np
from scipy.sparse.linalg import splu

from mouth_ringing_bend import digest, emit, setup
from pathlib import Path


def path(t,a):
    x=np.clip((t-a.start)/a.ramp,0.,1.)
    return (a.bend*x**3*(10-15*x+6*x*x),
            a.bend*30*x*x*(1-x)**2/a.ramp,
            a.bend*60*x*(1-x)*(1-2*x)/a.ramp**2)


def run(a):
    (m,n,m2,k,k1,k2,d),lift,good,probes=setup(a)
    kt=k1*.5
    fr=lambda s:s[good][:,good].tocsr()
    m0,n0,ms,k0,ks,kk,drag=map(fr,(m,n,m2,k,k1,k2,d))
    cross=n0+n0.T
    dt=a.dt
    a0=(m0+.5*dt*drag+.25*dt**2*k0).tocsr()
    ab= -cross+.25*dt**2*ks
    abb=ms+.25*dt**2*kk
    solve0=splu(a0.tocsc()).solve
    solve1=splu((a0+a.bend*ab+a.bend**2*abb).tocsc()).solve if a.bend else solve0
    # All boundary lifts and their geometry derivatives are assembled once.
    ml,nl,ntl,sl,kl,k1l,k2l,dl=[mat@lift for mat in (m,n,n.T,m2,k,k1,k2,d)]
    x,y=np.zeros(good.sum()),np.zeros(good.sum())
    u,v=np.zeros_like(lift),np.zeros_like(lift)
    wall_u=wall_v=launch_work=bend_work=lost=0.
    max_iteration=max_residual=error=peak_energy=0.
    trace=[]
    for i in range(round(a.end/dt)):
        t=(i+1)*dt
        b,bd,bdd=path(t-.5*dt,a)
        bend_next,bdnext,_=path(t,a)
        z=np.clip(t/3,0,1)
        next_v=a.stroke*10*z*z*(1-z)**2
        wv=.5*(next_v+wall_v); next_u=wall_u+dt*wv
        wu=.5*(wall_u+next_u)
        dv=next_v-wall_v
        changing=bool(bd or bdd)
        if changing:
            mass=m0-b*cross+b*b*ms
            stiff=k0+b*ks+b*b*kk-bdd*(n0.T-b*ms)
        elif b==0.:
            mass,stiff=m0,k0
        else:
            if i==round((a.start+a.ramp)/dt):
                final_m=m0-b*cross+b*b*ms;final_k=k0+b*ks+b*b*kk
            mass,stiff=final_m,final_k
        inert_l=ml-b*(nl+ntl)+b*b*sl
        damp_l=dl-2*bd*(ntl-b*sl)
        stiff_l=kl+b*k1l+b*b*k2l-bdd*(ntl-b*sl)
        load= -dv*inert_l-dt*(damp_l*wv+stiff_l*wu)
        rhs=mass@y-.5*dt*(stiff@x)+.5*load[good]
        if changing:
            perturb=(b*ab+b*b*abb
                     -(dt*bd+.25*dt**2*bdd)*(n0.T-b*ms))
            current=solve0(rhs)
            # Precondition by the undeformed wave matrix for this small bend.
            for it in range(15):
                correction=solve0(rhs-a0@current-perturb@current)
                current+=correction
                diff=np.linalg.norm(correction)/max(np.linalg.norm(current),1.e-30)
                if diff<2e-12:
                    break
            else:
                raise ValueError('moving-guide midpoint did not converge')
            max_iteration=max(max_iteration,it+1);max_residual=max(max_residual,diff)
        else:
            current=(solve0 if not b else solve1)(rhs)
        xn=x+dt*current; yn=2*current-y
        un=lift*next_u; un[good]=xn
        vn=lift*next_v; vn[good]=yn
        um=.5*(u+un);vel=.5*(vn+v);acc=(vn-v)/dt
        if wv:
            launch_work+=dt*wv*(inert_l@acc+damp_l@vel+stiff_l@um)
        dissip=dt*vel@(d@vel);lost+=dissip
        # Force from virtual acceleration of the physical w component and from
        # the stress energy. This includes the time derivative of wave momentum.
        metric=-vel@(n@vel)-um@(kt@um)
        force=um@(n@acc)-um@(kt@um)
        if b:
            metric+=b*(vel@(m2@vel)-um@(k2@um))
            force-=b*(um@(m2@acc)+um@(k2@um))
        if changing:
            force-=2*bd*(um@(m2@vel))+bdd*(um@(m2@um))
        bend_work-=dt*bd*force
        kin=vn@(m@vn)-2*bend_next*(vn@(n@vn))+bend_next**2*(vn@(m2@vn))
        if bdnext:
            kin+= -2*bdnext*(un@(n@vn)-bend_next*(un@(m2@vn)))+bdnext**2*un@(m2@un)
        energy=.5*(kin+un@(k@un)+bend_next*(un@(k1@un))+bend_next**2*un@(k2@un))
        p_end=un@(n@vn)-bend_next*(un@(m2@vn))
        trace.append([t,bend_next,launch_work,bend_work,lost,energy,metric,force,p_end,*list(probes@un)])
        error=max(error,abs(energy+lost-launch_work-bend_work))
        peak_energy=max(peak_energy,energy)
        u,v,x,y=un,vn,xn,yn
        wall_u,wall_v=next_u,next_v
    data=np.asarray(trace)
    file=a.outdir/'bend_trace.npz'
    np.savez_compressed(file,columns=['t','bend','work','bend_work','loss','energy','metric_force','flare_force','momentum','collar_stress','deep_stress'],readings=data)
    windows={}
    for lo,hi in ((24.,48.),(32.,56.),(40.,64.),(64.,88.)):
        if hi>a.end or lo<a.start+a.ramp:
            continue
        tau=(data[:,0]-.5*dt-lo)/(hi-lo)
        within=(tau>=0)&(tau<=1)
        weights=np.sin(np.pi*tau[within])**2
        values=np.average(data[within],axis=0,weights=weights)
        windows[f'{lo:g}-{hi:g}']=dict(wave_force_flare=values[7],wave_force_entrance=-values[7],
                                        metric_flare=values[6],energy=values[5])
    emit(test='poststroke_bend',bend=a.bend,amplitude=a.stroke,start=a.start,ramp=a.ramp,
         dt=dt,work=launch_work,bend_work=bend_work,energy=energy,max_energy_error=error,
         peak_energy=peak_energy,iteration_count=max_iteration,last_correction_bound=max_residual,
         hann_windows=windows,trace_sha256=digest(file))
    if error>max(2e-6*peak_energy,1e-20):
        raise ValueError('bend work budget is unresolved')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None,feed_spacing=.1)
    for key,default in [('bend',.01),('stroke',.01),('outer',11.),('cut',1.85),
                        ('dt',.025),('end',64.),('start',12.),('ramp',8.)]:
        p.add_argument('--'+key,type=float,default=default)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    if a.start<3 or a.ramp<=0 or a.end<=a.start+a.ramp or not all(np.isclose(z/a.dt,round(z/a.dt)) for z in (a.start,a.ramp,a.end)):
        p.error('Choose grid-aligned preparation, ramp and release intervals')
    run(a)
