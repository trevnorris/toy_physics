"""Static normal phase stiffness of a bent layer with and without an opening.

The two sides have the same neutral carrier. A thin layer gives an excess
resistance to normal phase current on surviving material. A missing layer
enforces phase continuity there. Geometry is prescribed in four spatial
dimensions with three rotations reduced to the radial area 4*pi*r^2.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from mouth_flat_phase_control import saved_bath
from mouth_free_front_wave_head import emit


def link_resistance(strength):
    # Neutral kink width and far carrier stiffness are one. Normal sheet
    # resistance is the excess over an equally thick homogeneous carrier.
    if strength == 1.:
        return np.inf
    a = np.sqrt(strength/(1-strength))
    return 4*a*np.arctan(a)


def gram(a,b):
    return a[:,:,None]*b[:,None,:]


def matrix(ids, local, size):
    return coo_matrix((local.ravel(), (np.repeat(ids,ids.shape[1],axis=1).ravel(),
                      np.tile(ids,(1,ids.shape[1])).ravel())),shape=(size,size)).tocsr()


def solve(mesh, strength, hole, shift, phase, flat=False):
    dr=16./mesh
    radii=np.linspace(0.,16.,mesh+1)
    height=np.zeros_like(radii) if flat else np.sqrt(100.**2-radii**2)-100.
    # A normal-work probe implemented as an axial displacement. Value and
    # slope vanish at r=8 and 16, so this mode does no work on either rim.
    t=np.clip((radii-8.)/8.,0.,1.)
    mode=16*t**2*(1.-t)**2
    height=height+shift*mode
    resistance=link_resistance(strength)
    low_ids=np.arange((mesh+1)**2).reshape(mesh+1,mesh+1)
    up_ids=np.empty_like(low_ids)
    shared=(radii<=8.) if hole else np.zeros_like(radii,dtype=bool)
    if resistance==0.:
        shared[:]=True
    up_ids[:,0]=np.where(shared,low_ids[:,-1],-1)
    up_ids[:,1:]=-1
    missing=up_ids<0
    up_ids[missing]=np.arange(low_ids.size,low_ids.size+np.count_nonzero(missing))
    size=up_ids.max()+1
    x,w=leggauss(3)
    x,w=(x+1)/2,w/2

    def volume(upper, labels):
        eta=np.linspace(0.,1.,mesh+1)
        if upper:
            z=height[:,None]*(1.-eta)+16.*eta
            velocity=mode[:,None]*(1.-eta)
        else:
            z=-16.*(1.-eta)+height[:,None]*eta
            velocity=mode[:,None]*eta
        ids=np.stack([labels[:-1,:-1],labels[1:,:-1],labels[1:,1:],labels[:-1,1:]],axis=-1).reshape(-1,4)
        def corners(field):
            return np.stack([field[:-1,:-1],field[1:,:-1],field[1:,1:],field[:-1,1:]],axis=-1).reshape(-1,4)
        zc,vc=corners(z),corners(velocity)
        radial_start=np.repeat(radii[:-1],mesh)
        local=np.zeros((len(ids),4,4))
        deriv=np.zeros_like(local)
        for i,s in enumerate(x):
            for j,t in enumerate(x):
                ds=np.array([t-1,1-t,t,-t])
                dt=np.array([s-1,-s,s,1-s])
                zs,zt=zc@ds,zc@dt
                vs,vt=vc@ds,vc@dt
                gw=dt[None,:]/zt[:,None]
                gr=(ds[None,:]-zs[:,None]*gw)/dr
                ratio=vt/zt
                dgw=-gw*ratio[:,None]
                dgr=-(vs[:,None]*gw+zs[:,None]*dgw)/dr
                weight=w[i]*w[j]*4*np.pi*(radial_start+s*dr)**2*dr*zt
                base=gram(gr,gr)+gram(gw,gw)
                cross=gram(gr,dgr)+gram(gw,dgw)
                local+=weight[:,None,None]*base
                deriv+=weight[:,None,None]*(cross+cross.transpose(0,2,1)+ratio[:,None,None]*base)
        return matrix(ids,local,size),matrix(ids,deriv,size)

    lower,dlower=volume(False,low_ids)
    upper,dupper=volume(True,up_ids)
    active=radii[:-1]>=8. if hole else np.ones(mesh,dtype=bool)
    lo,up=low_ids[:,-1],up_ids[:,0]
    ids=np.stack([lo[:-1],lo[1:],up[:-1],up[1:]],axis=1)
    g=1/resistance if resistance>0 else 0.
    slope=np.diff(height)/dr
    dslope=np.diff(mode)/dr
    metric=np.hypot(1.,slope)
    local=np.zeros((mesh,4,4))
    deriv=np.zeros_like(local)
    capillary_work=0.
    for s,weight in zip(x,w):
        v=np.array([s-1.,-s,1.-s,s])
        area=4*np.pi*(radii[:-1]+s*dr)**2*dr*weight
        local+=(g*area*metric*active)[:,None,None]*np.outer(v,v)
        deriv+=(g*area*slope*dslope/metric*active)[:,None,None]*np.outer(v,v)
        capillary_work+=np.sum((1/3)*area*slope*dslope/metric*active)
    junction=matrix(ids,local,size)
    A=lower+upper+junction
    dA=dlower+dupper+matrix(ids,deriv,size)
    fixed=np.r_[low_ids[:,0],up_ids[:,-1]]
    free=np.setdiff1d(np.arange(size),fixed)
    psi=np.zeros(size)
    psi[up_ids[:,-1]]=phase
    psi[free]=spsolve(A[free][:,free],-(A[free]@psi))
    port=A@psi
    flux=port[up_ids[:,-1]].sum()
    incoming=-port[low_ids[:,0]].sum()
    energy=.5*psi@(A@psi)
    generalized_work=.5*psi@(dA@psi)
    # For distinct traces read transmitted flux by integrating the link law.
    # The residual of the lower-side weak equation supplies aperture flux.
    seam=lower@psi
    lower_link=np.zeros(size)
    np.add.at(lower_link, ids[:,:2], np.einsum("eij,ej->ei",local[:,:2],psi[ids]))
    if resistance>0:
        leak=sum(np.sum(local[i,2:,:]@psi[ids[i]]) for i in np.flatnonzero(active))
    else:
        # For a transparent junction the exact phase is affine in physical w.
        # Projected area accounts for n_w dA. This is the open-sheet reference.
        leak=phase/32.*(4*np.pi/3)*(16.**3-(8.**3 if hole else 0.))
    aperture=0.
    if hole:
        aperture=np.sum((seam+lower_link)[lo[shared]]) if resistance>0 else flux-leak
    return dict(total_current=flux, intact_material_current=leak, aperture_current=aperture,
                incoming=incoming, aperture_plus_sheet_minus_port=aperture+leak-flux,
                largest_interior_residual=np.max(abs(port[free])), phase_energy=energy,
                twice_energy_minus_port_flux=2*energy-phase*flux,
                material_mode_capillary_hold=capillary_work, material_mode_phase_hold=generalized_work,
                peak_phase_jump=max(abs(psi[up]-psi[lo])), phase=phase)


def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--mesh',type=int,required=True)
    arg.add_argument('--strength',type=float,required=True)
    arg.add_argument('--sheet',choices=('intact','punctured'),required=True)
    arg.add_argument('--flat',action='store_true')
    arg.add_argument('--force-check',action='store_true')
    args=arg.parse_args()
    if not 0<=args.strength<=1 or args.mesh%2 or args.mesh<4:
        raise ValueError('Strength must be in [0,1] and mesh must be even and at least four')
    strength=args.strength
    res=link_resistance(strength)
    phase=saved_bath()/10
    hole=args.sheet=='punctured'
    numerical=None
    if strength<1:
        # Excess resistance from K_perp=1-lambda*4chi(1-chi) on a neutral kink.
        numerical=quad(lambda s: (strength/np.cosh(s/2)**2)/(1.-strength/np.cosh(s/2)**2),
                       -40.,40.,points=[0.],epsabs=2e-10)[0]
    emit(test='law',strength=strength, sheet=args.sheet, mesh=args.mesh,
         flat=args.flat, resistance=None if np.isinf(res) else res,
         zero_link=bool(np.isinf(res)), resistance_by_diffuse_integration=numerical,
         phase=phase, aperture_radius=8., box_radius=16., plate_distance=32., curvature_radius=100.)
    row=solve(args.mesh,strength,hole,0.,phase,args.flat)
    flat_reference=0. if np.isinf(res) else phase*4*np.pi*16**3/(3*(32+res))
    emit(test='transport',flat_full_sheet_reference=flat_reference,**row)
    scale=max(1.,abs(row['total_current']))
    if (abs(row['incoming']-row['total_current'])>scale*1e-9 or
        abs(row['aperture_plus_sheet_minus_port'])>scale*1e-9 or
        abs(row['twice_energy_minus_port_flux'])>2e-9 or row['largest_interior_residual']>3e-10):
        raise ArithmeticError('Current conservation or variational energy check failed')
    if numerical is not None and abs(numerical-res)>5e-8*max(1.,res):
        raise ArithmeticError('Diffuse resistance disagreement')
    if args.flat and not hole and abs(row['total_current']-flat_reference)>1e-8:
        raise ArithmeticError('Flat-sheet series-resistance check failed')
    if args.force_check:
        for step in (.01,.005):
            plus=solve(args.mesh,strength,hole,step,phase,args.flat)['phase_energy']
            minus=solve(args.mesh,strength,hole,-step,phase,args.flat)['phase_energy']
            force=(plus-minus)/(2*step)
            emit(test='held_mode',step=step,work_from_resolved_fields=force,
                 direct_variational_work=row['material_mode_phase_hold'], absolute_difference=abs(force-row['material_mode_phase_hold']))
            if abs(force-row['material_mode_phase_hold'])>3e-7:
                raise ArithmeticError('Force from independent displaced geometries disagrees')


if __name__=='__main__':
    main()
