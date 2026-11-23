"""Resolve the spatially varying steady state underlying the radius-two acoustic gate."""
import argparse
import hashlib
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_exchange_carrier import BIAS, REGIONS, mapped_nodes
from mouth_exit_witness import inputs
from mouth_live_waist_fields import Carrier, ORDER_INERTIA, WIDTH
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def reference(offset=0.):
    trial = Trial(2, 72, bonded_core=True)
    path = ROOT/'_scratch/s11c/mouth_exit_release01/r2n72/endpoint.npz'
    shifts, chi, _, digest, status, *_ = inputs(path, 2, 72, trial)
    archive = ROOT/'research/pde_ledger_v3/_measurements/mouth_breathing_exit_support.json'
    data = json.loads(archive.read_bytes())
    if data['physical_preparations']['r2n72'] != digest or status != 'balanced':
        raise ValueError('Unverified radius-two preparation')
    rows = [json.loads(line) for line in data['runs']['mouth_breathing_matched01/r2n72']['literal_stdout'].splitlines()]
    config = trial.x+shifts[:4]
    config[1] += np.log1p(offset)
    exterior = trial.cr.copy()
    exterior[0] += shifts[5]
    return Meridian(2, config), chi, trial.mesh, rows[0]['original_exit_scale'], exterior, rows, digest


class SteadyMouth:
    """Reconstruct Q1 fields for each side of the sealed wall, in physical (r,w)."""
    def __init__(self, shape, chi, field_mesh):
        n = field_mesh
        self.shape = shape
        self.radial, self.axial = mapped_nodes(n, shape.length, 28., 22.)
        carrier = Carrier(shape, chi.shape[1]-1, chi.shape[0]-1, n)
        energy, _, self.info = carrier.response(chi, BIAS)
        gr, gw = carrier.GradR@chi.ravel(), carrier.GradW@chi.ravel()
        rr, rw, ww = carrier.tensor
        r2 = gr**2+gw**2
        projection = gr[:,None,None]**2*rr + (gr*gw)[:,None,None]*rw + gw[:,None,None]**2*ww
        pieces = []
        self.fields = {}
        for region, ids, K, density_ids in carrier.parts:
            fcell = carrier.f[density_ids]
            rho = fcell.mean(axis=1)**2
            if region == 'core':
                K = K - (ORDER_INERTIA/(rho+ORDER_INERTIA*r2))[:,None,None]*projection
            factor = rho**2 if region == 'cap_layer' else rho
            pieces.append((factor[:,None,None]*K).ravel())
            self.fields[region] = (fcell, ids)
        matrix = sparse.coo_matrix((np.concatenate(pieces),(carrier.rows,carrier.cols)),shape=carrier.stiff.shape).tocsr()
        phase = np.zeros(matrix.shape[0])
        phase[carrier.bulk] = BIAS
        free = carrier.free
        phase[free] = spsolve(matrix[free][:,free], -(matrix@phase)[free])
        self.fields = {key:(f,phase[ids]) for key,(f,ids) in self.fields.items()}
        self.order_r, self.order_w = gr, gw
        self.mesh = n
        self.energy = energy
        current = (matrix@phase)[carrier.bulk].sum()
        if abs(current-self.info['terminal_current']) > 1e-6*abs(current):
            raise ArithmeticError('Reconstructed steady phase changed intake')

    def disk(self, r):
        """Mass current density evaluated using the assembled half-layer jump law."""
        y = np.clip(np.asarray(r)/self.shape.exit,0.,1.)
        edges = self.radial[0]
        i = np.clip(np.searchsorted(edges,y,side='right')-1,0,self.mesh-1)
        a = (y-edges[i])/(edges[i+1]-edges[i])
        f, phase = self.fields['cap_layer']
        lo, hi = (phase[i,2]-phase[i,0]), (phase[i,3]-phase[i,1])
        return f[i].mean(axis=-1)**4*((1-a)*lo+a*hi)/(WIDTH/12)

    def details(self, region, r, w):
        if isinstance(region,(int,np.integer)):
            region = REGIONS[region]
        slab, side = divmod(REGIONS.index(region),2)
        r, w = np.broadcast_arrays(r,w)
        rad, slope = self.shape.radius(w)
        width = 28.-rad if side else rad
        x = (r-rad)/width if side else r/width
        xs, zs = self.radial[side], self.axial[slab]
        i = np.clip(np.searchsorted(xs,x,side='right')-1,0,self.mesh-1)
        k = np.clip(np.searchsorted(zs,w,side='right')-1,0,self.mesh-1)
        a = np.clip((x-xs[i])/(xs[i+1]-xs[i]),0,1)
        b = np.clip((w-zs[k])/(zs[k+1]-zs[k]),0,1)
        cell = k*self.mesh+i
        def interpolant(values):
            v0,v1,v2,v3 = np.moveaxis(values[cell],-1,0)
            value = (1-b)*((1-a)*v0+a*v1)+b*((1-a)*v2+a*v3)
            dr = ((1-b)*(v1-v0)+b*(v3-v2))/(width*(xs[i+1]-xs[i]))
            dw = ((1-a)*(v2-v0)+a*(v3-v1))/(zs[k+1]-zs[k])
            dw -= ((1-x) if side else x)*slope*dr
            return value, dr, dw
        f, df_r, df_w = interpolant(self.fields[region][0])
        theta, u_r, u_w = interpolant(self.fields[region][1])
        density = f*f
        if region=='core':
            gr, gw = self.order_r[cell],self.order_w[cell]
            factor = ORDER_INERTIA*(u_r*gr+u_w*gw)/(density+ORDER_INERTIA*(gr*gr+gw*gw))
            u_r, u_w = u_r-factor*gr,u_w-factor*gw
        return density,u_r,u_w,df_r,df_w,theta

    def base(self,region,r,w):
        return self.details(region,r,w)[:3]


def diagnostics():
    q, weights = np.polynomial.legendre.leggauss(4)
    for offset in (0.,-.012,-.006,.006,.012):
        shape,chi,mesh,_,_,_,digest = reference(offset)
        for n in ((mesh,32,48) if offset==0 else (mesh,)):
            state=SteadyMouth(shape,chi,n)
            edges=state.radial[0]*shape.exit
            r=edges[:-1,None]+np.diff(edges)[:,None]*(1+q)/2
            area=4*np.pi*r*r*np.diff(edges)[:,None]*weights/2
            flux=state.disk(r)
            density=state.details('core',r,np.full_like(r,shape.length))[0]
            mach=flux/density**3   # u = j/n, c^2=n U''(n)=n^4 in the homogeneous law.
            total=np.sum(area*flux)
            emit(kind='steady_profile',offset=offset,field_mesh=n,material_sha256=digest,
                 disk_radius=shape.exit,depth=shape.length,flow=total,mean_current=total/area.sum(),
                 current_min=float(flux.min()),current_max=float(flux.max()),max_normal_mach=float(mach.max()),
                 area_fraction_local_mach_over_one=float(area[mach>1].sum()/area.sum()),
                 throughflow_fraction_local_mach_over_one=float((area*flux)[mach>1].sum()/total),
                 radius_fraction=(r[:,2]/shape.exit).tolist(),current_density=flux[:,2].tolist(),
                 cap_input_residual=abs(total-state.info['terminal_current']))
            if abs(total-state.info['terminal_current'])>1e-6*abs(total) or flux.min()<=0:
                raise ArithmeticError('Disk quadrature and independently stored current disagree')


if __name__=='__main__':
    diagnostics()
