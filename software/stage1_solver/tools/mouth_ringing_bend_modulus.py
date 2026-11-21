"""Straight-wave curvature term for a shear modulus sensitive to slow material strain.

mu=mu0 exp(lambda tr[(F.T F-I)^2]/2); slow section shear is b h'(w).
At zero bend the additional flare force per b and lambda is -int mu C² h'².
Use the same monopole edge field, contact history and energy norm as the
geometric force calculation; both tangents and quiet-wave effects subtract.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_ringing_bend import emit, mesh_geometry, prepared, quadratic, setup, taper, tensor


def weight_matrix(a):
    dom,chi,_,_=prepared(a)
    indices=(0,1,2,5)
    grids=[dom.grids[i] for i in indices]
    p,v=grids[2]
    cut=p[np.argmin(abs(p-a.cut))]
    grids[2]=(p[p<=cut],v)
    class Charts:
        @staticmethod
        def chart(i,p,v):
            return dom.chart(indices[i],p,v)
    _,blocks=mesh_geometry(Charts(),grids,a.sub,False)
    cells=np.concatenate([block[2] for block in blocks])
    pairs=cells[:,np.array([[0,2],[1,3],[0,1],[2,3]])]
    edges,ids=np.unique(np.sort(pairs,axis=-1).reshape(-1,2),axis=0,return_inverse=True)
    ids=ids.reshape(-1,4)
    sg=np.sign(pairs[:,:,1]-pairs[:,:,0])*np.array([1.,-1.,-1.,1.])
    # Only conduit elements can have nonzero h'.
    start=sum(len(b[2]) for b in blocks[:2])
    pp,vv,elems=blocks[2]
    ed=ids[start:start+len(elems)]
    sign=sg[start:start+len(elems)]
    i,j=[m.ravel() for m in np.indices((len(pp)-1,len(vv)-1))]
    roots,wg=leggauss(4)
    t,v=[m.ravel() for m in np.meshgrid((1+roots)/2,(1+roots)/2,indexing='ij')]
    dp,dv=np.diff(pp)[i,None],np.diff(vv)[j,None]
    p,z=pp[i,None]+dp*t,vv[j,None]+dv*v
    r,w,rp,rv,wp,wv=dom.chart(2,p,z)
    jac=dp*dv*abs(rp*wv-rv*wp)
    po,vo,old=dom.order_blocks[2]
    oi,oj=i//a.sub,j//a.sub
    ax=quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
    ay=quadratic(2*(z-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
    order=np.einsum('eqi,ei->eq',tensor(ax,ay),chi[old[oi*(len(vo)-1)+oj]])
    per=(4*np.pi*r*r*order*taper(w,.55,1.15)[1]**2/jac)@(np.outer(wg,wg).ravel()/4)
    local=per[:,None,None]*sign[:,:,None]*sign[:,None,:]
    return sparse.coo_matrix((local.ravel(),(np.repeat(ed,4,axis=1).ravel(),
           np.tile(ed,(1,4)).ravel())),shape=(len(edges),len(edges))).tocsr()


def run(a):
    w=weight_matrix(a)
    mats,lift,good,_=setup(a)
    ne=w.shape[0]
    m,k,d=[mats[i][:ne,:ne].tocsr() for i in (0,3,6)]
    good,lift=good[:ne],lift[:ne]
    mf,kf,df=[mat[good][:,good] for mat in (m,k,d)]
    dt=a.dt
    solve=splu((mf+.5*dt*df+.25*dt**2*kf).tocsc()).solve
    ml,kl,dl=[mat@lift for mat in (m,k,d)]
    x=np.zeros(good.sum());y=x.copy();wall_q=wall_v=work=0.
    result=[]
    for i in range(round(a.end/dt)):
        t=(i+1)*dt
        g=np.clip(t/3,0,1)
        nv=.01*10*g*g*(1-g)**2
        vm=.5*(wall_v+nv); nq=wall_q+dt*vm; qm=.5*(wall_q+nq)
        rhs=mf@y-.5*dt*(kf@x)-.5*(ml[good]*(nv-wall_v)+dt*(kl[good]*qm+dl[good]*vm))
        ym=solve(rhs); yn=2*ym-y; xn=x+dt*ym
        u=lift*qm; u[good]=.5*(x+xn)
        if vm:
            vel=lift*vm;vel[good]=ym
            work+=vm*(ml[good]@(yn-y)+(ml@lift)*(nv-wall_v)+dt*(kl@u+dl@vel))
        result.append((t-.5*dt,u@(w@u)))
        x,y,wall_q,wall_v=xn,yn,nq,nv
    times,weights=np.array(result).T
    avgs={}
    for lo,hi in ((24,48),(32,56),(40,64)):
        if hi<=a.end:
            sel=(times>lo)&(times<hi)
            window=np.sin(np.pi*(times[sel]-lo)/(hi-lo))**2
            avgs[f'{lo}-{hi}']=np.average(weights[sel],weights=window)
    emit(test='incremental_modulus',work=work,dt=dt,mesh=a.mesh,sub=a.sub,stiffness_per_lambda=avgs)
    if min(avgs.values())<=0:
        raise ValueError('no positive occupied-mode overlap')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None,feed_spacing=.1,mesh=34,outer=11.,cut=1.85)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--dt',type=float,default=.025)
    p.add_argument('--end',type=float,default=64.)
    run(p.parse_args())
