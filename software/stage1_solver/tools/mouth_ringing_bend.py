"""Two-end conservative force of the formation wave on an attached, bent material guide.

Reference sections translate by d(w) = q_e + (q_f-q_e) h(w). Material
coefficients follow the sections. Pull back the displacement one-form, use
the Euclidean physical curl and vector norm, and retain the spherical sector
and every polar dipole component about the translation axis. In the material
frame, time differentiation of the physical wave displacement also gives
a polarization momentum term when varying the bend velocity.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_inertia import taper
from mouth_formation_wall_work import CENTER, WIDTH, profile
from mouth_order_wave_pulse import digest, emit, mesh_geometry, quadratic, tensor
from mouth_shear_ringdown import prepared


def sparse_local(rows,cols,data,size):
    return sparse.coo_matrix((data.ravel(),(rows,cols)),shape=(size,size)).tocsr()


def setup(a):
    dom,chi,intake,path=prepared(a)
    indices=(0,1,2,5)
    grids=[dom.grids[i] for i in indices]
    pg,vg=grids[2]
    cut=pg[np.argmin(abs(pg-a.cut))]
    grids[2]=(pg[pg<=cut],vg)
    class Charts:
        @staticmethod
        def chart(i,p,v):
            return dom.chart(indices[i],p,v)
    pts,blocks=mesh_geometry(Charts(),grids,a.sub,False)
    # Global edge circulations with orientation from smaller to larger node ID.
    boundaries=np.array([[0,2],[1,3],[0,1],[2,3]])
    all_cells=np.concatenate([block[2] for block in blocks])
    pairs=all_cells[:,boundaries]
    edges,ids=np.unique(np.sort(pairs,axis=-1).reshape(-1,2),axis=0,return_inverse=True)
    edge_ids=ids.reshape(-1,4)
    orient=np.sign(pairs[:,:,1]-pairs[:,:,0])
    ne=len(edges)
    total=2*ne+len(pts)
    # M, N, M2, K, K1, K2, D. M(d)=M-d*(N+N.T)+d² M2.
    mats=[sparse.csr_matrix((total,total)) for _ in range(7)]
    lift=np.zeros(total)
    fixed=np.zeros(total,dtype=bool)
    axis=(pts[:,0]==0.)
    fixed[2*ne:]=axis
    fixed[ne+np.flatnonzero(np.all(axis[edges],axis=1))]=True
    probes=np.zeros((2,total))
    norms=np.zeros(2)
    batch=0
    xx,ww=leggauss(4)
    t,v=[q.ravel() for q in np.meshgrid((xx+1)/2,(xx+1)/2,indexing='ij')]
    weight=np.outer(ww,ww).ravel()/4

    for which,(pp,vv,cells) in zip(indices,blocks,strict=True):
        nrow,ncol=len(pp)-1,len(vv)-1
        ib,jb=[m.ravel() for m in np.indices((nrow,ncol))]
        po,vo,old=dom.order_blocks[which]
        for offset in range(0,len(cells),384):
            sel=slice(offset,offset+384)
            i,j=ib[sel],jb[sel]
            n=len(i)
            raw=np.arange(batch+offset,batch+offset+n)
            ed=edge_ids[raw]
            sg=orient[raw,None,:]
            loc=np.column_stack([ed,ne+ed,2*ne+cells[sel]])
            rows=np.repeat(loc,12,axis=1).ravel()
            cols=np.tile(loc,(1,12)).ravel()
            def operands(t,v):
                dp,dv=np.diff(pp)[i,None],np.diff(vv)[j,None]
                p,z=pp[i,None]+dp*t,vv[j,None]+dv*v
                r,w,rp,rv,wp,wv=dom.chart(which,p,z)
                rt,ry,wt,wy=rp*dp,rv*dv,wp*dp,wv*dv
                jac=rt*wy-ry*wt
                oi,oj=i//a.sub,j//a.sub
                ax=quadratic(2*(p-po[oi,None])/np.diff(po)[oi,None]-1)[0]
                ay=quadratic(2*(z-vo[oj,None])/np.diff(vo)[oj,None]-1)[0]
                order=np.einsum('eqi,ei->eq',tensor(ax,ay),chi[old[oi*(len(vo)-1)+oj]])
                null=np.zeros_like(t)
                at=sg*np.stack((1-v,v,null,null),axis=-1)
                av=sg*np.stack((null,null,1-t,t),axis=-1)
                ur=(wy[:,:,None]*at-wt[:,:,None]*av)/jac[:,:,None]
                uw=(-ry[:,:,None]*at+rt[:,:,None]*av)/jac[:,:,None]
                curl=sg*np.array((1.,-1.,-1.,1.))/jac[:,:,None]
                ns=np.stack(((1-t)*(1-v),(1-t)*v,t*(1-v),t*v),axis=-1)
                nt=np.stack((v-1,-v,1-v,v),axis=-1)
                nv=np.stack((t-1,1-t,-t,t),axis=-1)
                nr=(wy[:,:,None]*nt-wt[:,:,None]*nv)/jac[:,:,None]
                nw=(-ry[:,:,None]*nt+rt[:,:,None]*nv)/jac[:,:,None]
                shape=(n,len(t),12)
                U,V,P,Q,B,C,C1,D,E=(np.zeros(shape) for _ in range(9))
                U[:,:,:4],V[:,:,:4],C[:,:,:4]=ur,uw,curl
                P[:,:,4:8],Q[:,:,4:8],C1[:,:,4:8]=ur,uw,curl
                B[:,:,8:]=ns/r[:,:,None]
                D[:,:,4:8],D[:,:,8:]=-ur/r[:,:,None],nr/r[:,:,None]
                E[:,:,4:8],E[:,:,8:]=uw/r[:,:,None],-nw/r[:,:,None]
                return r,w,jac,rt,ry,wt,wy,order,(U,V,P,Q,B,C,C1,D,E)

            r,w,jac,rt,ry,wt,wy,order,fields=operands(t,v)
            U,V,P,Q,B,C,C1,D,E=fields
            weight0=4*np.pi*r*r*abs(jac)*order*weight
            dh=-taper(w,.55,1.15)[1]
            def gram(f,g,wt):
                return np.einsum('eqi,eqj->eij',wt[:,:,None]*f,g)
            def term(f,wt=weight0):
                return gram(f,f,wt)
            nmat=(gram(U,Q,weight0*dh)+gram(P+2*B,V,weight0*dh))/3
            mx=weight0*dh**2
            # <u_x²>=U²/3+(P+2B)²/9+4(P-B)²/45, with <n_x²>=1/3.
            groups=(term(U)+term(V)+(term(P)+term(Q)+2*term(B))/3,
                    nmat,term(U,mx)/3+term(P+2*B,mx)/9+4*term(P-B,mx)/45,
                    term(C)+(term(C1)+2*term(E)+2*term(D))/3,
                    -2*(gram(C,D,weight0*dh)+gram(D,C,weight0*dh))/3,
                    2*term(D,mx)/3)
            for k,values in enumerate(groups):
                mats[k]+=sparse_local(rows,cols,values,total)
            if which==2:
                wall=j==0
                if np.any(wall):
                    e=ed[wall,0]
                    fixed[e]=fixed[ne+e]=True
                    fixed[2*ne+cells[sel][wall][:,[0,2]]]=True
                    zw=pp[i[wall],None]+np.diff(pp)[i[wall],None]*((xx+1)/2)
                    lift[e]=orient[raw[wall],0]*np.diff(pp)[i[wall]]*(profile(zw)@ww/2)
                    # Read physical round-component shear stress just inside the wall.
                    rr,zz,_,_,_,_,_,ch,wall_fields=operands((xx+1)/2,np.zeros_like(xx))
                    for k,depth in enumerate((CENTER,1.12)):
                        band=np.maximum(0.,1-((zz[wall]-depth)/.045)**2)**2
                        wt0=band*ww*np.diff(pp)[i[wall],None]/2
                        norms[k]+=wt0.sum()
                        data=(wt0[:,:,None]*ch[wall,:,None]*wall_fields[5][wall]).sum(axis=1)
                        np.add.at(probes[k],loc[wall].ravel(),data.ravel())
            if which in (2,5):
                on=i==nrow-1
                if np.any(on):
                    rr,zz,_,_,ry,_,wy,ch,edge_fields=operands(np.ones_like(xx),(xx+1)/2)
                    U,V,P,Q,B,*_=edge_fields
                    wt0=4*np.pi*rr**2*np.hypot(ry,wy)*ch*(ww/2)*on[:,None]
                    mer,dip=(U,P) if which==2 else (V,Q)
                    loss=term(mer,wt0)+(term(dip,wt0)+2*term(B,wt0))/3
                    mats[6]+=sparse_local(rows,cols,loss,total)
        batch+=len(cells)
    probes/=norms[:,None]
    for mat in mats:
        mat.eliminate_zeros()
    good=~fixed
    emit(test='mesh',order_path=str(path),order_sha256=digest(path),intake=intake,cut=cut,
         cells=len(all_cells),dofs=total,free_dofs=int(good.sum()),sub=a.sub,mesh=a.mesh,outer=a.outer)
    return mats,lift,good,probes


def run(a):
    (m,n,m2,k,k1,k2,damp),lift,good,probes=setup(a)
    b=a.bend
    mass=(m-b*(n+n.T)+b*b*m2).tocsr()
    stiff=(k+b*k1+b*b*k2).tocsr()
    inertia=n-b*m2
    half_k=k1*.5+b*k2
    # Use the attached-material kinetic energy of d(F^-T u)/dt. Its extra
    # derivative contributes +d/dt [u.T (N-b*M2) udot] to the flare force.
    dt=a.dt
    free=lambda mat:mat[good][:,good].tocsr()
    mf,kf,df=map(free,(mass,stiff,damp))
    solve=splu((mf+dt/2*df+dt**2/4*kf).tocsc()).solve
    ml,kl,dl=[(mat@lift)[good] for mat in (mass,stiff,damp)]
    q=np.zeros(len(good)); v=np.zeros_like(q)
    x=np.zeros(good.sum()); y=np.zeros_like(x)
    wall_q=wall_v=work=loss=momentum=0.
    bands=((8,12),(20,24),(28,32),(40,44),(44,48))
    records={span:np.zeros(5) for span in bands if a.end>=span[1]}
    trace=[]
    max_error=0.
    for j in range(round(a.end/dt)):
        time=(j+1)*dt
        g=np.clip(time/3,0,1)
        next_w=a.stroke*30*g*g*(1-g)**2/3
        vm=.5*(wall_v+next_w)
        next_q=wall_q+dt*vm
        source= -ml*(next_w-wall_v)-dt*(kl*(wall_q+next_q)*.5+dl*vm)
        ym=solve(mf@y-.5*dt*(kf@x)+.5*source)
        xnew=x+dt*ym; ynew=2*ym-y
        u=q.copy(); u[good]=(x+xnew)*.5; u[~good]=lift[~good]*(wall_q+next_q)*.5
        vel=v.copy();vel[good]=ym;vel[~good]=lift[~good]*vm
        q[good]=xnew;q[~good]=lift[~good]*next_q
        v[good]=ynew;v[~good]=lift[~good]*next_w
        if vm:
            work+=vm*(ml@(ynew-y)+lift@(mass@lift)*(next_w-wall_v)+dt*((stiff@lift)@u+(damp@lift)@vel))
        loss+=dt*vel@(damp@vel)
        metric=-vel@(inertia@vel)-u@(half_k@u)   # wave force on positive-depth end
        next_mom=q@(inertia@v)
        force=metric+(next_mom-momentum)/dt
        energy=.5*(v@(mass@v)+q@(stiff@q))
        err=energy+loss-work
        max_error=max(max_error,abs(err))
        for (lo,hi),s in records.items():
            if lo<time<=hi+1e-10:
                s+=dt*np.array((metric,force,-force,energy,loss))
        trace.append([time,work,loss,energy,metric,force,next_mom,*list(probes@q)])
        x,y=xnew,ynew
        wall_q,wall_v,momentum=next_q,next_w,next_mom
    file=a.outdir/'bend_trace.npz'
    np.savez_compressed(file,readings=trace,columns=['t','work','loss','energy','metric_force','flare_force','momentum','collar_stress','deep_stress'])
    emit(test='bend_response',bend=a.bend,stroke=a.stroke,dt=a.dt,end=a.end,
         work=work,max_balance_error=max_error,
         means={f'{lo}-{hi}':dict(zip(('metric_flare_force','wave_flare_force','wave_entrance_force','energy','port_loss'),
                                  (arr/(hi-lo)).tolist(),strict=True)) for (lo,hi),arr in records.items()},
         trace_sha256=digest(file))
    if max_error>2e-8*max(abs(work),1e-25):
        raise ValueError('wall work, field storage and losses fail the midpoint balance')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None,feed_spacing=.1)
    p.add_argument('--bend',type=float,required=True)
    p.add_argument('--stroke',type=float,default=.01)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--outer',type=float,default=11.)
    p.add_argument('--cut',type=float,default=1.85)
    p.add_argument('--dt',type=float,default=.025)
    p.add_argument('--end',type=float,default=48.)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    if a.sub<1 or not 0<a.dt<=.1 or not np.isclose(1/a.dt,round(1/a.dt)) or a.end<3 or a.outer<=9:
        p.error('Use a positive mesh, integer sampling windows and extended feed')
    run(a)
