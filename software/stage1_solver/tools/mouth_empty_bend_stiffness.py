"""Wave-empty carrier stiffness of the two native lateral bends, with fixed remote baths."""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_carrier_free_center import basis
from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import Dipole, prepare, trace
from mouth_flared_receiver import emit

INTERVALS = ((.55,1.15),(1.2,2.))


def quadrature(dom):
    roots,_=leggauss(4)
    xi,eta=(z.ravel() for z in np.meshgrid(roots,roots,indexing='ij'))
    radii=[]; depths=[]; radial=[]; axial=[]; columns=[]
    for block,(pg,vg,dofs) in dom.blocks.items():
        i,j=(z.ravel()[:,None] for z in np.indices((len(pg)-1,len(vg)-1)))
        p=pg[i]+np.diff(pg)[i]*(1+xi)/2
        v=vg[j]+np.diff(vg)[j]*(1+eta)/2
        nodes,r,_,_,_,gr,gw=basis(dom,block,p,v,i,j)
        radii.append(r.ravel()); depths.append(dom.chart(block,p,v)[1].ravel())
        radial.append(gr.reshape(-1,9)); axial.append(gw.reshape(-1,9))
        columns.append(np.broadcast_to(nodes,gr.shape).reshape(-1,9))
    cols=np.concatenate(columns).ravel()
    rows=np.arange(0,len(cols)+1,9)
    def matrix(parts):
        return sparse.csr_matrix((np.concatenate(parts).ravel(),cols,rows),shape=dom.N.shape)
    return np.concatenate(radii),np.concatenate(depths),matrix(radial),matrix(axial)


def shapes(r,w,extension,far,roof):
    if extension=='moving-box':
        rgate=np.ones_like(r); zgate=np.ones_like(w)
        dr=np.zeros_like(r); dz=np.zeros_like(w)
    else:
        lo=2.0 if extension=='compact' else 3.4
        rgate,dr=taper(r,lo,far-.35)
        zgate,dz=taper(w,3.3 if extension=='compact' else 4.3,roof-.35)
    hs=[]; sr=[]; sw=[]
    for lo,hi in INTERVALS:
        near,der=taper(w,lo,hi)
        h=1-near
        hs.append(h*rgate*zgate)
        sr.append(h*dr*zgate)
        sw.append(rgate*(h*dz-der*zgate))
    return np.column_stack(hs),np.column_stack(sr),np.column_stack(sw)


def wall_slope(dom,field,tangents):
    """Natural-wall pressure from the carrier stress, independently evaluated on its trace.

    Shifts are constant over each wall sphere. Under the x-only deformation,
    n_x dA is invariant. The reflecting entrance is held fixed in these two tests.
    """
    pg=dom.blocks[2][0]; x,weight=leggauss(9)
    hp=np.diff(pg)[:,None]/2
    p=(pg[:-1,None]+hp*(1+x)).ravel()
    r,w,_,_,_,_,f,fr,fw=trace(dom,2,p,np.zeros_like(p),np.column_stack([field,tangents]))
    h,_,hw=shapes(r,w,'moving-box',dom.far,dom.roof)
    re=lambda a,b: (a.conj()*b).real
    dlag=(re(fr[:,0,None],fr[:,1:])+re(fw[:,0,None],fw[:,1:]-hw*fr[:,0,None])
          +((abs(f[:,0,None])**8-1)/2)*re(f[:,0,None],f[:,1:]))
    measure=4*np.pi/3*r*r*(hp*weight).ravel()
    return h.T@(measure[:,None]*dlag)


def run(a):
    dom,field=prepare(a)
    r,w,gr,gw=quadrature(dom)
    _,sr,sw=shapes(r,w,a.extension,dom.far,dom.roof)
    f,fr,fw=dom.N@field,gr@(field-1),gw@(field-1)
    q=dom.weight
    drive=sr*fr[:,None]+sw*fw[:,None]
    force=.5*(gr.T@(q[:,None]*drive)
            +gw.T@(q[:,None]*(sw*fr[:,None]-sr*fw[:,None]))
            +dom.N.T@(q[:,None]*(2*drive/r[:,None]-sr*((abs(f)**8-1)*f/2)[:,None])))
    d=Dipole(dom,field)
    pairs=[d.solve(np.zeros(len(field),complex),force[:,i]) for i in (0,1)]
    v=np.column_stack([z for z,_ in pairs])
    weight=q*abs(fr)**2/3
    frozen=sr.T@(weight[:,None]*sr)+sw.T@(weight[:,None]*sw)
    relaxation=2*(force.conj().T@v).real/3
    stiffness=frozen-relaxation
    wall=wall_slope(dom,field,v)
    emit(test='empty_carrier',mesh=len(dom.blocks[2][1])-1,extension=a.extension,far=dom.far,roof=dom.roof,
         frozen=frozen.tolist(),relaxation=relaxation.tolist(),holding_stiffness=stiffness.tolist(),
         surface_only=wall.tolist(),eigenvalues=np.linalg.eigvalsh(stiffness).tolist(),
         dipole_residuals=[err for _,err in pairs],reciprocity_error=np.max(abs(stiffness-stiffness.T)))
    if np.max(abs(stiffness-stiffness.T))>1e-9 or not np.isfinite(stiffness).all():
        raise ArithmeticError('Carrier work reciprocity failed')
    if a.save:
        np.savez_compressed(a.save,field=field,tangents=v,stiffness=stiffness,points=dom.points)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(profile='cubic',rate=.25)
    p.add_argument('--source',type=Path)
    p.add_argument('--mesh',type=int,default=20)
    p.add_argument('--far',type=float,default=9.)
    p.add_argument('--roof',type=float,default=10.)
    p.add_argument('--extension',choices=('compact','wide','moving-box'),default='compact')
    p.add_argument('--save',type=Path)
    return p


if __name__=='__main__':
    run(parser().parse_args())
