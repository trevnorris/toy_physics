"""Follow the open degree-four torsional pole on the volume-neutral neck/flare path.

The uniform feeding slab is terminated by its separated outgoing, frequency-
dependent modal traction. The exterior bulk is shear-free. The carrier and
geometric forces are separate observables. Time dependence is exp(-i omega t).
"""
import argparse
import json
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh, splu

from mouth_exchange_geometry import ExchangeShape

DEGREE = 4
STEPS = (0., -.01, -.005, .005, .01)


def uniform_grid(start, end, cells):
    return np.linspace(start, end, cells+1)


def block_matrices(radial, axial, mapping, fade, length, degree):
    # Bilinear quadrilaterals in a meridian with exact S² angular integral.
    xx, zz = np.meshgrid(radial, axial, indexing='ij')
    r, bprime = mapping(xx, zz)
    points = np.c_[r.ravel(),zz.ravel()]
    nodes = np.arange(xx.size).reshape(xx.shape)
    cells = np.stack((nodes[:-1,:-1], nodes[1:,:-1],nodes[:-1,1:],nodes[1:,1:]),axis=2).reshape(-1,4)
    x0,z0 = np.meshgrid(radial[:-1],axial[:-1],indexing='ij')
    hx,hz = np.meshgrid(np.diff(radial),np.diff(axial),indexing='ij')
    x0,z0,hx,hz = [v.ravel() for v in (x0,z0,hx,hz)]
    km = np.zeros((len(cells),4,4)); mm=np.zeros_like(km)
    gauss, weights=leggauss(3)
    for a,wa in zip((gauss+1)/2,weights/2):
        for b,wb in zip((gauss+1)/2,weights/2):
            x=x0+a*hx
            w=z0+b*hz
            r,rw=mapping(x,w)
            r0,_=mapping(x0,w)
            r1,_=mapping(x0+hx,w)
            dr=r1-r0
            if np.any(dr<=0):
                raise ArithmeticError('Folded radial coordinates')
            n=np.array([(1-a)*(1-b),a*(1-b),(1-a)*b,a*b])
            nr=np.array([b-1,1-b,-b,b])/dr[:,None]
            nw=np.array([a-1,-a,1-a,a])/hz[:,None]-rw[:,None]*nr
            # Both shear modulus and transverse inertia vanish in the bulk.
            # If a diffuse outlet is requested, interpolate them together.
            chi=np.ones_like(r)
            if fade:
                t=np.clip((length-w)/fade,0.,1.)
                chi=t*t*(3-2*t)
            weight=wa*wb*dr*hz*chi*r**2
            curl=2*n+r[:,None]*nr
            km += weight[:,None,None]*(degree*(degree+1)*np.outer(n,n)
                   +curl[:,:,None]*curl[:,None,:]+r[:,None,None]**2*nw[:,:,None]*nw[:,None,:])
            mm += (weight*r*r)[:,None,None]*np.outer(n,n)
    return points,cells,km,mm


def assemble(shape, mesh, port_radius, fade):
    neck, span = shape.neck, shape.length
    nx=mesh
    nz=mesh
    ns=int(round(mesh*(port_radius-1.)))
    nt=int(np.ceil(span*mesh))
    lead_w=uniform_grid(-1.,0.,nz)
    throat_w=uniform_grid(0.,span,nt)
    if fade:
        # End of the constant-modulus region is a fixed material coordinate.
        j=np.argmin(abs(throat_w-(span-fade)))
        if j not in (0,len(throat_w)-1):
            throat_w[j]=span-fade
    inner=uniform_grid(0.,neck,nx)
    side=uniform_grid(neck,port_radius,ns)
    def slab(x,w):
        return x,np.zeros_like(w)
    def conduit(x,w):
        b,ratio=shape.radius(w)
        return x*b,x*ratio
    chunks=[block_matrices(inner,lead_w,slab,0.,span,DEGREE),
            block_matrices(side,lead_w,slab,0.,span,DEGREE),
            block_matrices(uniform_grid(0.,1.,nx),throat_w,conduit,fade,span,DEGREE)]
    # Number shared nodes by first appearance in the material chart. Sorting
    # by physical radius instead would reorder unrelated field components when
    # radial grid lines pass each other during a smooth deformation.
    coordinates=np.round(np.vstack([b[0] for b in chunks]),10)
    _,first,index=np.unique(coordinates,axis=0,return_index=True,return_inverse=True)
    ordering=np.argsort(first)
    remap=np.empty_like(ordering);remap[ordering]=np.arange(len(ordering))
    points=coordinates[first[ordering]]
    index=remap[index]
    offset=0; cells=[]; matrices=[]; masses=[]
    for pts,ids,k,m in chunks:
        cells.append(index[offset+ids]); matrices.append(k); masses.append(m)
        offset+=len(pts)
    ids=np.vstack(cells)
    rows,cols=np.repeat(ids,4,axis=1).ravel(),np.tile(ids,(1,4)).ravel()
    live=points[:,0]>0
    def mat(data):
        a=sparse.coo_matrix((np.vstack(data).ravel(),(rows,cols)),shape=(len(points),len(points))).tocsr()
        return a[live][:,live].tocsc()
    # Group mass measures needed to distinguish a throat pole from a slab pole.
    throat_mass=mat([np.zeros_like(masses[0]),np.zeros_like(masses[1]),masses[2]])
    return points[live],mat(matrices),mat(masses),throat_mass


class OutgoingSlab:
    """FE depth eigenchannels; exact continuous radial traction in every channel."""
    def __init__(self, points, radius):
        self.radius=radius
        self.ids=np.flatnonzero(abs(points[:,0]-radius)<1.e-9)
        self.ids=self.ids[np.argsort(points[self.ids,1])]
        z=points[self.ids,1]
        if abs(z[0]+1)>1e-10 or abs(z[-1])>1e-10:
            raise ArithmeticError('Open port does not span entire feeding slab')
        mass=np.zeros((len(z),len(z))); stiff=np.zeros_like(mass)
        for i,h in enumerate(np.diff(z)):
            mass[i:i+2,i:i+2]+=h*np.array([[2.,1.],[1.,2.]])/6
            stiff[i:i+2,i:i+2]+=np.array([[1.,-1.],[-1.,1.]])/h
        lam,V=eigh(stiff,mass)
        self.cutoff=np.sqrt(np.maximum(lam,0.))
        self.project=V.T@mass
        self.coeff=np.array([1j**n*math.factorial(DEGREE+n)/(math.factorial(n)*math.factorial(DEGREE-n)*2**n)
                             for n in range(DEGREE+1)])

    def traction(self, omega):
        k=np.sqrt((omega+0j)**2-self.cutoff**2)
        # Continue the propagating waves outward and evanescent waves in their
        # decaying sheet from positive real omega to its lower-half-plane pole.
        k=np.where((omega.real<self.cutoff)&(k.imag<0),-k,k)
        y=1/(self.radius*k)
        P=np.polynomial.polynomial.polyval(y,self.coeff)
        dP=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(self.coeff))
        return 1j*k-y*dP/(self.radius*P)

    def block(self,omega):
        return self.radius**4*(self.project.T*self.traction(omega))@self.project


def solve(shape, mesh, port_radius, fade, prior=None):
    points,stiff,mass,inside=assemble(shape,mesh,port_radius,fade)
    port=OutgoingSlab(points,port_radius)
    cavity=points[:,1]>0
    val,vec=eigsh(stiff[cavity][:,cavity],k=1,M=mass[cavity][:,cavity],sigma=0.,
                  which='LM',v0=np.ones(cavity.sum()),tol=1e-11)
    seed=np.zeros(len(points));seed[cavity]=vec[:,0]
    omega=np.sqrt(val[0])-1e-4j if prior is None else prior
    v=seed.astype(complex); seed_dual=mass@seed
    port_ids=port.ids
    edge=(np.repeat(port_ids,len(port_ids)),np.tile(port_ids,len(port_ids)))
    for i in range(28):
        B=port.block(omega)
        dB=(port.block(omega+2e-6)-port.block(omega-2e-6))/4e-6
        op=stiff-(omega**2)*mass-sparse.coo_matrix((B.ravel(),edge),shape=mass.shape).tocsc()
        res=op@v
        residual=np.linalg.norm(res)/max(1e-15,np.linalg.norm(stiff@v))
        if residual<3e-10 and abs(seed_dual@v-1)<1e-9:
            break
        deriv=(-2*omega)*(mass@v)
        deriv[port_ids]-=dB@v[port_ids]
        # Bordered Newton avoids inverting an exactly singular resonance operator.
        tangent=sparse.bmat([[op,sparse.csc_matrix(deriv[:,None])],
                             [sparse.csc_matrix(seed_dual[None,:]),None]],format='csc')
        correction=splu(tangent).solve(-np.r_[res,seed_dual@v-1])
        v+=correction[:-1];omega+=correction[-1]
    else:
        raise ArithmeticError(f'Open pole failed, residual {residual}, last frequency {omega}')
    norm=np.vdot(v,mass@v).real
    trapped=np.vdot(v,inside@v).real/norm
    # Flux at the real frequency: power/energy of the mode inside the mouth/slab cut.
    ph=port.project@v[port_ids]
    tr=port.traction(omega.real)
    powers=.5*omega.real*port_radius**4*abs(ph)**2*tr.imag
    kinetic=.25*abs(omega)**2*norm
    strain=.25*np.vdot(v,stiff@v).real
    gamma=-omega.imag
    if not (0<gamma and trapped>.5 and omega.real>0 and np.sum(powers)>0):
        raise ArithmeticError(f'Wrong outgoing or delocalized branch {omega}, {trapped}')
    # Integrate the shape variation of physical curl energy and kinetic energy
    # in comoving trial functions. The real energy work omits exterior storage;
    # the bilinear pole variation additionally accounts for the exact open port.
    delta=1e-4
    form=[]
    for change in (-delta,delta):
        coords, kshift, mshift, _=assemble(type(shape)(shape.step+change),mesh,port_radius,fade)
        if len(coords)!=len(points):
            raise ArithmeticError('Virtual displacement changed topology')
        form.append((kshift,mshift))
    dk=(form[1][0]-form[0][0])/(2*delta)
    dm=(form[1][1]-form[0][1])/(2*delta)
    stored_work = omega.real*(np.vdot(v,dk@v)-abs(omega)**2*np.vdot(v,dm@v)).real/(4*(kinetic+strain))
    db=(port.block(omega+2e-6)-port.block(omega-2e-6))/4e-6
    pole_work=(v@(dk@v)-omega**2*(v@(dm@v)))/(2*omega*(v@(mass@v))+v[port_ids]@(db@v[port_ids]))
    return dict(frequency=omega.real, stored_energy_slope=stored_work,
                complex_pole_slope_real=pole_work.real, complex_pole_slope_imag=pole_work.imag,amplitude_loss=gamma,energy_lifetime=1/(2*gamma),
                period_count_per_energy_efold=omega.real/(4*np.pi*gamma),
                confined_kinetic_fraction=trapped,
                clamped_seed_overlap=abs(seed_dual@v)**2/norm,
                outward_power_over_interior_energy=powers.sum()/(kinetic+strain),
                open_channels=int(sum(port.cutoff<omega.real)),pole_residual=residual,
                newton_iterations=i+1,nodes=len(v)),v/np.sqrt(norm),mass


def run(a):
    rows=[];reference=None
    for step in STEPS:
        out,v,m=solve(ExchangeShape(step),a.mesh,a.port,a.fade,reference)
        if step==0:
            reference=complex(out['frequency'],-out['amplitude_loss'])
            v0=v
        overlap=abs(np.vdot(v0,m@v))**2/np.vdot(v0,m@v0).real
        out.update(step=step,overlap_with_zero_on_material_nodes=overlap)
        if len(v)!=len(v0) or overlap<.96:
            raise ArithmeticError('Mode correspondence failed')
        rows.append(out)
        print(json.dumps(dict(test='mode',mesh=a.mesh,port=a.port,fade=a.fade,**out)),flush=True)
    for s in (.005,.01):
        v={r['step']:r for r in rows}
        w0=v[0]['frequency'];wp=v[s]['frequency'];wm=v[-s]['frequency']
        print(json.dumps(dict(test='signed',mesh=a.mesh,port=a.port,fade=a.fade,step=s,
                              frequency_first=(wp-wm)/(2*s),frequency_second=(wp+wm-2*w0)/s**2,
                              energy_work_curvature=(v[s]['stored_energy_slope']-v[-s]['stored_energy_slope'])/(2*s),
                              differentiated_pole_curvature=(v[s]['complex_pole_slope_real']-v[-s]['complex_pole_slope_real'])/(2*s))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--port',type=float,default=3.)
    p.add_argument('--fade',type=float,default=0.)
    run(p.parse_args())
