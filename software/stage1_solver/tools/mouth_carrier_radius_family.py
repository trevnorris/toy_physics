"""Change physical mouth radius at fixed feed thickness and material coefficients.

The entire rounded lip, collar and flare scale with the mouth. The lateral bath
extent and roof are fixed in physical units. A common imposed phase difference
is one candidate preparation; each mouth's material intake is an output.
"""
import argparse
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_carrier_lateral_response import Dipole, weighted_gram
from mouth_carrier_lateral_inertia import WINDOWS, volume
from mouth_flared_receiver import EXIT, ROUNDING, FlaredDomain, emit


class RadiusFamily(FlaredDomain):
    def __init__(self, radius, depth_unit, mesh, law, outer_radius, roof, thickness):
        self.length_unit=1.
        self.depth_stretch=depth_unit/radius
        self.feed_depth=thickness/depth_unit
        super().__init__(mesh,law,round(roof/depth_unit,11),far=round(outer_radius/radius,11))
        self.length_unit=radius
        self.points=self.stretched_points*radius
        del self.stretched_points
        self.stiff=self.stiff*radius**2
        self.weight*=radius**4
        np.testing.assert_allclose(self.points[:,1].min(),-thickness,atol=2e-10)
        np.testing.assert_allclose(self.points[:,0].max(),outer_radius,atol=2e-10)
        radial_edge=np.isclose(self.points[:,0],outer_radius,rtol=0,atol=2e-9)
        heights=self.points[:,1]
        np.testing.assert_array_equal(self.feed,radial_edge & (heights<=1e-10))
        np.testing.assert_array_equal(self.bulk,(radial_edge & (heights>=depth_unit*EXIT-1e-10)) |
                                      np.isclose(heights,roof,atol=2e-9,rtol=0))

    def assemble(self,grids):
        self.bottom_turn=-np.arctan((self.feed_depth+ROUNDING)/self.rc)
        mesh=len(grids[0][1])-1
        breaks=[-np.pi/2,np.arctan(self.slope)-np.pi/2,self.bottom_turn,0.]
        angular=np.r_[np.concatenate([np.linspace(a,b,max(1,int(np.ceil((b-a)*mesh))),
                              endpoint=False) for a,b in zip(breaks[:-1],breaks[1:])]),0.]
        super().assemble([(angular,grids[0][1]),*grids[1:]])
        # Parent identifies the far baths using chart coordinates. Keep the
        # physical stretched points and leave chart coordinates for that step.
        self.stretched_points=self.points
        self.points=np.round(self.points/[1.,self.depth_stretch],11)

    def reach(self,beta):
        m,dm=super().reach(beta)
        low=beta<self.bottom_turn
        sine=np.where(low,np.sin(beta),-1.)
        extent=-(self.feed_depth+ROUNDING)/sine
        return np.where(low,extent,m),np.where(low,-extent*np.cos(beta)/sine,dm)

    def chart(self,block,p,v):
        if block==1:
            one=np.ones_like(p)
            zero=np.zeros_like(p)
            coords=(p,-self.feed_depth*v,one,zero,zero,-self.feed_depth*one)
        else:
            coords=super().chart(block,p,v)
        r,w,rp,rv,wp,wv=coords
        return tuple(self.length_unit*q for q in
                     (r,w*self.depth_stretch,rp,rv,wp*self.depth_stretch,wv*self.depth_stretch))


def fixed_head(dom,bias):
    """Solve native nonlinear complex equation with independently given bath phases."""
    base=dom.bulk.astype(float)
    mask=dom.free
    k=dom.stiff*.5
    kf=k[mask][:,mask]
    base[mask]=spsolve(kf,-(k@base)[mask])
    field=np.exp(1j*bias*base)
    nf=dom.N[:,mask]
    for iteration in range(20):
        f=dom.N@field
        n=abs(f)**2
        scalar=dom.weight*(n**4-1)/4
        residual=k@(field-1)+dom.N.T@(scalar*f)
        error=float(max(abs(residual[mask])))
        if error<3e-12:
            break
        slope=2*dom.weight*n**3
        mixed=weighted_gram(nf,slope*f.real*f.imag)
        operator=sparse.bmat([[kf+weighted_gram(nf,scalar+slope*f.real**2),mixed],
                              [mixed,kf+weighted_gram(nf,scalar+slope*f.imag**2)]],format='csc')
        step=spsolve(operator,-np.r_[residual[mask].real,residual[mask].imag])
        x,y=np.split(step,2)
        field[mask]+=x+1j*y
    else:
        raise ArithmeticError('The fixed-bias carrier did not settle')
    conjugate=(field.conj()*residual).imag
    feed=-2*conjugate[dom.feed].sum()
    bulk=2*conjugate[dom.bulk].sum()
    if abs(feed-bulk)>1e-8 or min(n)<.5:
        raise ArithmeticError(f'Reservoir currents or density failed: {feed=}, {bulk=}, {error=}, {min(n)=}')
    return field,dict(feed=feed,bulk=bulk,min_density=min(n),residual=error,iterations=iteration)


def run(args):
    r=args.radius
    depth=r if args.depth=='scaled' else 1.
    dom=RadiusFamily(r,depth,args.mesh,'cubic',args.outer,args.roof,args.thickness)
    field,base=fixed_head(dom,args.bias)
    emit(test='family_preparation',radius=r,mesh=args.mesh,thickness=args.thickness,depth_law=args.depth,depth_unit=depth,
         outer_radius=args.outer,roof=args.roof,bias=args.bias,**base)
    model=Dipole(dom,field)
    feed,_=model.steady(1.,0.)
    bulk,_=model.steady(0.,1.)
    tow,_=model.steady(-1.,-1.)
    source=dom.N.T@(dom.weight*(1j*(dom.N@tow)-model.rq*(dom.N@field)))
    acceleration,_=model.solve(np.zeros_like(field),source)
    windows=np.array(WINDOWS)*[r,r,depth,depth]
    if windows[:,1].max()>=args.outer or windows[:,3].max()>=args.roof:
        raise ValueError('Control volume intersects driven outer boundary')
    for w in windows:
        forces=[volume(dom,np.column_stack([field,h]),w,False)['virtual_wall_force'] for h in (feed,bulk,tow)]
        mass=-volume(dom,np.column_stack([field,acceleration,tow]),w,True)['virtual_wall_force']
        emit(test='family_response',radius=r,mesh=args.mesh,window=w.tolist(),
             intake=base['feed'],held_feed=forces[0],held_bulk=forces[1],tow=forces[2],
             mass=mass,acceleration_per_incident_speed=forces[0]/mass)
        if abs(sum(forces))>2e-9 or mass<=0:
            raise ArithmeticError('Independent forcing or reactive-mass check failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--depth',choices=('scaled','fixed'),default='scaled')
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--bias',type=float,required=True)
    p.add_argument('--outer',type=float,default=14.)
    p.add_argument('--roof',type=float,default=15.)
    p.add_argument('--thickness',type=float,default=1.)
    run(p.parse_args())
