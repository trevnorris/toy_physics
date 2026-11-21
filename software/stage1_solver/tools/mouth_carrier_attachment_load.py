"""Find the internal load needed to carry the whole shape at one acceleration.

At leading order in axial bias, two independent transverse wall motions have a
positive fluid mass matrix in the resting domain. Integrate the held weak-flow
tractions against the same two displacements. Their sum translates the full
shape. The partition can be shifted along w as a control.
"""
import argparse
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu,spsolve

from mouth_carrier_free_center import basis
from mouth_carrier_lateral_response import quadrature_radii, weighted_gram
from mouth_carrier_lateral_inertia import taper
from mouth_carrier_weak_radius import PhysicalDomain
from mouth_flared_receiver import emit


def ports(dom,field,flow,switch):
    root,weights=leggauss(10)
    mass_source=np.zeros((len(field),2))
    forces=np.zeros((2,2))
    for block in (0,2):
        pg,_,_=dom.blocks[block]
        i=np.arange(len(pg)-1)[:,None]
        hp=np.diff(pg)[:,None]/2
        p=pg[i]+hp*(1+root)
        v=np.zeros_like(p)
        idx,r,wp,_,n,gr,gw=basis(dom,block,p,v,i,np.zeros_like(i))
        w=dom.chart(block,p,v)[1]
        near,_=taper(w,*switch)
        shape=np.stack([near,1-near],axis=-1)
        # Axial carrier phase and two angular stream potentials.
        bfield=field[idx]
        fr=(gr*bfield).sum(axis=-1)
        fw=(gw*bfield).sum(axis=-1)
        u=flow[idx]
        ur=(gr[...,None]*u).sum(axis=-2)
        uw=(gw[...,None]*u).sum(axis=-2)
        darea=4*np.pi*r*r*wp*hp*weights
        sources=np.einsum('eq,eqi,eqj->eij',darea,n,shape)
        np.add.at(mass_source,idx[:,0,:].ravel(),sources.reshape(-1,2))
        forces-=np.einsum('eq,eqi,eqj->ij',darea/3,shape,fr[...,None]*ur+fw[...,None]*uw)
    return mass_source,forces


def run(a):
    dom=PhysicalDomain(1.,1.,a.mesh,a.outer,a.roof,1.)
    k=dom.stiff
    baseline=dom.bulk.astype(float)
    baseline[dom.free]=spsolve(k[dom.free][:,dom.free],-(k@baseline)[dom.free])
    D=baseline@(k@baseline)
    scalar=k+weighted_gram(dom.N,2*dom.weight/quadrature_radii(dom)**2)
    interior=dom.free & (dom.points[:,0]>0)
    lu=splu(scalar[interior][:,interior].tocsc())
    flow=dom.points[:,0,None]*np.c_[dom.feed,dom.bulk]
    flow[interior]=lu.solve(-(scalar@flow)[interior])
    # Normalize the base flux to one. A and B below are per intake, while M is at zero intake.
    for switch in ((.55,1.15),(1.2,2.0)):
        source,f=ports(dom,baseline/D,flow,switch)
        potential=np.zeros_like(source)
        potential[interior]=lu.solve(source[interior])
        mass=potential.T@(scalar@potential)/3
        np.testing.assert_allclose(mass,source.T@potential/3,rtol=3e-10)
        eig=np.linalg.eigvalsh(mass)
        if min(eig)<=0:
            raise ArithmeticError('Independent displacement fields did not yield positive inertia')
        one=np.ones(2)
        m=mass@one
        total=sum(m)
        held=f[:,0]
        acommon=sum(held)/total
        internal=m*acommon-held
        # Positive connector force here acts on the far coordinate.
        # Unconnected accelerations ignore prestress and bending elasticity.
        uncoupled=np.linalg.solve(mass,held)
        reduced=np.linalg.det(mass)/total
        emit(test='attachment',mesh=a.mesh,outer=a.outer,roof=a.roof,switch=switch,
             total_mass=total,mass_matrix=mass.tolist(),common_mass_shares=m.tolist(),
             held_per_intake=held.tolist(),bulk_per_intake=f[:,1].tolist(),
             tow_per_intake=(-f.sum(axis=1)).tolist(),common_a_per_incident_per_intake=acommon,
             connector_on_distal_per_held_force=internal[1]/sum(held),
             internal_forces=internal.tolist(),uncoupled_a_per_incident_per_intake=uncoupled.tolist(),
             relative_mode_mass=reduced)
        if abs(sum(internal))>1e-10 or abs(sum(held)-1)>.005:
            raise ArithmeticError('Unclosed internal transfer or mismatched axial input')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--outer',type=float,default=14.)
    p.add_argument('--roof',type=float,default=15.)
    run(p.parse_args())
