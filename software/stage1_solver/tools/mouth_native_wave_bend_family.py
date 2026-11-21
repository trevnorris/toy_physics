"""Read the lateral curvature for four prepared native patterns and two independent bends."""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh,splu

from mouth_native_wave_bend import (Bend,angular_stiffness,assemble,bend_forms,digest,emit,find_pole,
    impedance,momentum_and_mass,pencil,prepared,truncate)


def run(a):
    dom,chi,intake,src=prepared(a)
    dom,chi,cut=truncate(dom,chi,a.cut)
    pts,blocks,_,parts,outer,_=assemble(dom,chi,a.sub,1,a.feed_spacing)
    masses,_=momentum_and_mass(dom,blocks,chi,len(pts),0.)
    ang=sum(angular_stiffness(dom,blocks,chi,len(pts)))
    k0=sum(parts)
    cz,_,_=impedance(dom,blocks,pts,chi,cut,0.)
    mass=masses.sum(axis=0)
    visible=mass>0
    if a.boundary=='waist': visible &= pts[:,1]>.5000000001
    live={l:visible&((pts[:,0]>0)|(l==1)) for l in range(1,6)}
    scale={l:sparse.diags(1/np.sqrt(mass[ix])) for l,ix in live.items()}
    def restrict(mat,i,j):
        return (scale[i]@mat[live[i]][:,live[j]]@scale[j]).tocsr()
    k={l:restrict(k0+(l*(l+1)-2)*ang,l,l) for l in live}
    damp={l:((outer+cz)[ix]/mass[ix] if a.boundary=='open' else np.zeros(ix.sum())) for l,ix in live.items()}
    bends=[Bend(bend_forms(dom,blocks[2],chi,len(pts),s)) for s in ((.55,1.15),(1.2,2.))]
    emit(test='settings',radius=a.radius,mesh=a.mesh,sub=a.sub,cut=cut,outer=a.outer,boundary=a.boundary,
         intake=intake,order_sha256=digest(src),coordinates=['collar offset','deeper offset'])
    for l in range(1,5):
        behind=pts[live[l],1]>.5000000001
        vals,vec=eigsh(k[l][behind][:,behind],k=1,sigma=0,v0=np.ones(behind.sum()),tol=1e-12)
        seed=np.zeros(live[l].sum()); seed[behind]=vec[:,0]
        omega,u,err=find_pole(k[l],damp[l],np.sqrt(vals[0])+1e-4j,seed)
        ed=(abs(omega)**2+np.vdot(u,k[l]@u).real)/4
        # Disjoint depth intervals imply zero explicit (frozen) mixed second derivative.
        direct=[restrict(b.square(l,l),l,l)@u for b in bends]
        conjugated=np.diag([np.vdot(u,v) for v in direct])
        analytic=np.diag([u@v for v in direct])
        for j in (l-1,l+1):
            if j==0: continue
            source=np.column_stack([restrict(b.linear(j,l),j,l)@u for b in bends])
            response=-splu(pencil(k[j],damp[j],omega)).solve(source)
            conjugated+=source.conj().T@response
            analytic+=source.T@response
        work=conjugated.real/(2*ed)
        pole=(2*analytic/(u@((2*omega+1j*damp[l])*u))).real/omega.real
        sym=(work+work.T)/2
        eig,dirs=np.linalg.eigh(sym)
        # Collar fraction t and deeper fraction 1-t give total end offset one.
        direction=np.array([1.,-1.])
        base=np.array([0.,1.])
        denom=direction@sym@direction
        fractions=[0.,1.]
        if abs(denom)>1e-14:
            fractions.append(float(np.clip(-(direction@sym@base)/denom,0,1)))
        values=[(base+t*direction)@sym@(base+t*direction) for t in fractions]
        emit(test='pattern',degree=l,omega=omega.real,loss=-omega.imag,residual=err,
             seed_overlap_squared=abs(seed@u)**2,mean_work_stiffness=work.tolist(),pole_stiffness=pole.tolist(),
             shape_eigenvalues=eig.tolist(),positive_eigenvector_if_any=dirs[:,-1].tolist() if eig[-1]>0 else None,
             fraction_for_maximum=fractions[np.argmax(values)],max_monotone_stiffness=max(values),
             fraction_for_minimum=fractions[np.argmin(values)],min_monotone_stiffness=min(values))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--cut',type=float,default=2.15)
    p.add_argument('--outer',type=float,default=16.)
    p.add_argument('--feed-spacing',type=float,default=.09)
    p.add_argument('--boundary',choices=('waist','open'),default='open')
    run(p.parse_args())
