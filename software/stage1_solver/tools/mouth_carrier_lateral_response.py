"""Dipolar response of the stationary carrier on the four-dimensional flared mouth.

Coordinates are (r,w) with exact integration over brane directions on S².
The round opening has no boundary data. Two independent remote reservoirs
drive the density and phase response through the same curved throat.
"""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_curved_phase import quadratic, tensor
from mouth_curved_phase_stability import weighted_gram
from mouth_flared_receiver import EXIT, ROUNDING, FlaredDomain, emit
from mouth_free_contact import ARCHIVE


def prepare(args):
    if args.source:
        folder = args.source
        record = json.loads((folder/'stdout').read_text().splitlines()[0])
        dom = FlaredDomain(record['n'],record['profile'],record['roof'],record['feed_far'])
        path = folder/'field.npz'
        with path.open('rb') as stream:
            digest = file_digest(stream,'sha256').hexdigest()
        if digest!=json.loads(ARCHIVE.read_text())['output_field_sha256'][str(path)]:
            raise ValueError('The archived carrier field changed')
        with np.load(path) as data:
            np.testing.assert_array_equal(dom.points,data['points'])
            field = data['psi']
        emit(test='preparation',source=str(path),sha256=digest,mesh=record['n'],
             far=dom.far,roof=dom.roof,rate=record['bulk'],profile=dom.profile)
    else:
        dom = FlaredDomain(args.mesh,args.profile,args.roof,args.far)
        if args.rate==0:
            field = np.ones(len(dom.points),complex)
            record = dict(feed=0.,bulk=0.,residual=0.)
        else:
            field,record = dom.solve(args.rate)
        emit(test='preparation',source='new',mesh=args.mesh,far=dom.far,roof=dom.roof,
             profile=dom.profile,rate=args.rate,**record)
    return dom,field


def quadrature_radii(dom):
    roots,_ = leggauss(4)
    u,v = [x.ravel() for x in np.meshgrid(roots,roots,indexing='ij')]
    chunks = []
    for block,(pg,vg,dofs) in dom.blocks.items():
        i,j = np.unravel_index(np.arange(len(dofs)),(len(pg)-1,len(vg)-1))
        p = pg[i,None]+.5*np.diff(pg)[i,None]*(1+u)
        t = vg[j,None]+.5*np.diff(vg)[j,None]*(1+v)
        chunks.append(dom.chart(block,p,t)[0].ravel())
    return np.concatenate(chunks)


class Dipole:
    def __init__(self,dom,field,degree=1,active=None):
        """Linearized carrier in one angular harmonic with fixed complementary data."""
        if degree < 1 or degree != int(degree):
            raise ValueError("Expected a positive angular degree")
        self.dom, self.field = dom,field
        self.count = len(field)
        self.rq = quadrature_radii(dom)
        N = dom.N
        f = N@field
        n = abs(f)**2
        scalar = dom.stiff*.5 + weighted_gram(N,dom.weight*((n**4-1)/4 + degree*(degree+1)/(2*self.rq**2)))
        slope = 2*dom.weight*n**3
        cross = weighted_gram(N,slope*f.real*f.imag)
        self.operator = sparse.bmat([[scalar+weighted_gram(N,slope*f.real**2), cross],
                                     [cross,scalar+weighted_gram(N,slope*f.imag**2)]],format='csr')
        regular = dom.free & (dom.points[:,0]>0)
        if active is not None:
            regular &= active
        self.free = np.tile(regular,2)
        self.factor = splu(self.operator[self.free][:,self.free].tocsc())

    def solve(self,boundary,source=None):
        values = np.r_[boundary.real,boundary.imag]
        rhs = -(self.operator@values)
        if source is not None:
            rhs += np.r_[source.real,source.imag]
        values[self.free] = self.factor.solve(rhs[self.free])
        resid = self.operator@values
        if source is not None:
            resid -= np.r_[source.real,source.imag]
        error = float(max(abs(resid[self.free])))
        if error>3e-8:
            raise ArithmeticError(f'Dipole equation residual {error}')
        return values[:self.count]+1j*values[self.count:],error

    def steady(self,feed,bulk):
        r = self.dom.points[:,0]
        bath = feed*self.dom.feed + bulk*self.dom.bulk
        return self.solve(1j*self.field*r*bath)


def trace(dom,block,p,v,fields):
    pg,vg,dofs = dom.blocks[block]
    i = np.clip(np.searchsorted(pg,p,side='right')-1,0,len(pg)-2)
    j = np.clip(np.searchsorted(vg,v,side='right')-1,0,len(vg)-2)
    hp,hv = np.diff(pg)[i]/2,np.diff(vg)[j]/2
    bp,dp = quadratic((p-pg[i])/hp-1)
    bv,dv = quadratic((v-vg[j])/hv-1)
    nodal = fields[dofs[i*(len(vg)-1)+j]]
    def product(shape):
        return np.einsum('...i,...ic->...c',shape,nodal)
    f = product(tensor(bp,bv))
    fp = product(tensor(dp,bv))/hp[...,None]
    fv = product(tensor(bp,dv))/hv[...,None]
    r,w,rp,rv,wp,wv = dom.chart(block,p,v)
    det = rp*wv-rv*wp
    fr = (wv[...,None]*fp-wp[...,None]*fv)/det[...,None]
    fw = (rp[...,None]*fv-rv[...,None]*fp)/det[...,None]
    return r,w,rp,rv,wp,wv,f,fr,fw


def surface(dom,fields,block,along,value,outward,clip=None,accelerating=False):
    """Integrate through outward meridian normal. First column is the base, next the drive."""
    pg,vg,_ = dom.blocks[block]
    edges = pg if along=='p' else vg
    if clip is not None:
        lo,hi=clip
        edges=np.r_[lo,edges[(edges>lo)&(edges<hi)],hi]
    root,weights = leggauss(8)
    half=np.diff(edges)/2
    q=(edges[:-1,None]+half[:,None]*(root+1)).ravel()
    dq=(half[:,None]*weights).ravel()
    fixed=np.full_like(q,value)
    p,v=(q,fixed) if along=='p' else (fixed,q)
    r,_,rp,rv,wp,wv,f,fr,fw=trace(dom,block,p,v,fields)
    # The normal includes the line Jacobian. Select the sign by its dot
    # product with a vector pointing into the neighboring exterior.
    tr,tw=(rp,wp) if along=='p' else (rv,wv)
    nr,nw=tw,-tr
    projection=nr*outward[0]+nw*outward[1]
    if np.any(abs(projection)<1e-12):
        raise ValueError('Ambiguous surface orientation')
    sign=np.sign(projection)
    nr,nw=nr*sign,nw*sign
    fn=nr[:,None]*fr+nw[:,None]*fw
    n=abs(f[:,0])**2
    dlag=(fr[:,0].conj()*fr[:,1]+fw[:,0].conj()*fw[:,1]).real
    dlag+=(n**4-1)/2*(f[:,0].conj()*f[:,1]).real
    if accelerating:
        # Coaccelerating gauge includes +a*x. Time derivative of the velocity
        # solution enters the isotropic momentum stress in the same frame.
        dlag += n*r+(f[:,0].conj()*f[:,2]).imag
    pressure=-nr*dlag
    derivative=(fr[:,0].conj()*fn[:,1]+(fr[:,1]+2*f[:,1]/r).conj()*fn[:,0]).real
    measure=4*np.pi/3*r**2*dq
    return float(measure@(pressure+derivative)),float(measure@pressure)


def surfaces(dom,fields,accelerating=False):
    def face(block,along,val,normal,clip=None):
        return surface(dom,fields,block,along,val,normal,clip,accelerating)
    wall = {}
    # Every reflecting boundary with a lateral normal is in these three pieces.
    wall['rounded_lip']=face(0,'p',0.,(1.,1.))
    wall['collar']=face(2,'p',0.,(1.,0.),(ROUNDING,1.2))
    wall['flare']=face(2,'p',0.,(1.,0.),(1.2,EXIT))
    # The solid is OUTSIDE the fluid at r > a(w).
    # On the rounded feed lip it lies toward (r=rc,w=.5).
    wall['feed_plane']=face(1,'p',0.,(0.,1.))
    wall['feed_floor']=face(1,'p',1.,(0.,-1.))
    # On part of block zero, v=1 traces the reflecting feed floor until it hits the axis.
    wall['inner_floor']=face(0,'p',1.,(0.,-1.),(-np.pi/2,dom.bottom_turn))
    wall['bulk_plane']=face(4,'v',EXIT,(0.,-1.))
    outer = {'feed':face(1,'v',dom.far,(1.,0.))[0],
             'bulk_side':face(4,'p',1.,(1.,0.))[0],
             'roof_inner':face(3,'v',dom.roof,(0.,1.))[0],
             'roof_outer':face(4,'v',dom.roof,(0.,1.))[0]}
    cuts = {str(z):face(2,'v',z,(0.,1.))[0] for z in (.7,1.4,2.4,3.0)}
    return wall,outer,cuts


def steady_result(dom,field,h,drive,err):
    wall,outer,cuts=surfaces(dom,np.column_stack([field,h]))
    pressure={k:v[1] for k,v in wall.items()}
    # Exact slip trace gives the pressure-only force on the rigid shape.
    force=sum(pressure.values())
    balance=force+sum(outer.values())
    emit(test='steady',drive=drive,dipole_residual=err,wall_force=force,wall_pieces=pressure,
         wall_full_trace={k:v[0] for k,v in wall.items()},outer_flux=outer,
         section_flux=cuts,outward_sum=balance)
    return force


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path)
    p.add_argument('--mesh',type=int,default=22)
    p.add_argument('--far',type=float,default=9.)
    p.add_argument('--roof',type=float,default=10.)
    p.add_argument('--rate',type=float,default=.25)
    p.add_argument('--profile',choices=('cubic','delayed'),default='cubic')
    return p


def main():
    args=parser().parse_args()
    dom,field=prepare(args)
    model=Dipole(dom,field)
    feed,ef=model.steady(1.,0.)
    bulk,eb=model.steady(0.,1.)
    A=steady_result(dom,field,feed,'feed',ef)
    B=steady_result(dom,field,bulk,'bulk',eb)
    tow,et=model.steady(-1.,-1.)
    T=steady_result(dom,field,tow,'tow',et)
    emit(test='superposition',A=A,B=B,tow=T,independent_tow_error=abs(T+A+B),
         field_error=float(max(abs(tow+feed+bulk))))
    if max(abs(tow+feed+bulk))>3e-9 or abs(T+A+B)>1e-9:
        raise ArithmeticError('Independent towing solve failed superposition')


if __name__=='__main__':
    main()
