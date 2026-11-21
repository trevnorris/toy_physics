"""Reciprocal anti-plane momentum transport beside two layers of converting sites.

The channel is periodic across y and averaged along x over eight patch lengths.
The only velocity component is u_x(y,z). The mass of a cell is rho Lx dy dz.
Isotropic kinematic viscosity nu and the termination are trial material data.
"""
import argparse

import numpy as np

from mouth_order_wave_pulse import emit
from neutral_patch_sites import COUNT, GAP

WIDTH = 8.
LENGTH = 8.
TOP = 2.5
SUPPORT = .375


def kernel(offset):
    """Compact, even regularization in physical patch length units."""
    return np.where(abs(offset) < SUPPORT, 1 + np.cos(np.pi*offset/SUPPORT), 0.)


class ShearChannel:
    def __init__(self, grid, back_gap, backing, nu, rho, seeds_cases, upper0, warm=True, back_impedance=None):
        if (grid not in (4,8) or backing not in ('rest','free','robin')
                or back_gap <= .5 or not np.isclose(grid*back_gap,round(grid*back_gap))
                or rho <= 0 or nu <= 0):
            raise ValueError('Unsupported channel grid, termination, or transport coefficient')
        self.backing = backing
        self.nu, self.rho = nu, rho
        self.B = int(seeds_cases)
        self.ny = int(WIDTH*grid)
        self.nz = round((TOP+back_gap)*grid)
        self.h = 1/grid
        self.y = -WIDTH/2 + (np.arange(self.ny)+.5)*self.h
        self.z = -back_gap + (np.arange(self.nz)+.5)*self.h
        self.mass = rho*LENGTH*self.h**2
        # Slip length is viscosity divided by backing traction per unit speed.
        self.slip_length = 0.
        if backing=='robin':
            if back_impedance is None or not 0 < back_impedance < np.inf:
                raise ValueError('Robin backing requires positive finite impedance')
            self.slip_length = rho*nu/back_impedance
        self.bottom_factor = 0. if backing=='free' else self.h/(self.h/2+self.slip_length)
        self.traction_scale = rho*nu*LENGTH  # viscosity * Lx * dy/dz
        self.low = round((back_gap-.5)*grid)
        self.high = round((back_gap+.5)*grid)
        if self.low <= 0:
            raise ValueError('Lower CVs must use different lower faces')
        if np.asarray(upper0).shape != (self.B,):
            raise ValueError('Need one upper-wall velocity for every copy')
        self.u = np.zeros((self.B,self.nz,self.ny))
        if warm:
            steady = (np.ones(self.nz) if backing=='free' else
                      (self.z+back_gap+self.slip_length)/(TOP+back_gap+self.slip_length))
            self.u[:] = np.asarray(upper0)[:,None,None]*steady[None,:,None]
        z_weights = np.stack([kernel(self.z),kernel(self.z-GAP)])
        z_weights /= z_weights.sum(axis=1)[:,None]
        self.z_weights = []
        for w in z_weights:
            cells, = np.nonzero(w)
            self.z_weights.append((slice(cells[0],cells[-1]+1),w[cells,None]))
        if (np.any(z_weights[0,:self.low]) or np.any(z_weights[0,self.high:])
                or np.any(z_weights[1,:self.high])):
            raise ValueError('Control region does not separate the regularized site planes')
        self.offsets = np.arange(int(np.ceil(2*SUPPORT*grid))+1)
        self.batch_layers = (2*np.arange(self.B)[:,None]+np.repeat([0,1],COUNT))*self.ny
        self.size = self.B*2*self.ny

    def sample(self, y):
        """Sample local velocity and save the same weights for the reaction."""
        if y.shape != (self.B,2*COUNT):
            raise ValueError('Unexpected material patch layout')
        indices = np.floor((y+WIDTH/2-SUPPORT)/self.h-.5).astype(int)[:,:,None]+1+self.offsets
        displacement = (indices+.5)*self.h-(y[:,:,None]+WIDTH/2)
        weight = kernel(displacement)
        weight /= weight.sum(axis=2,keepdims=True)
        bins = indices%self.ny + self.batch_layers[:,:,None]
        local = np.stack([(self.u[:,cells,:]*w).sum(axis=1) for cells,w in self.z_weights],axis=1)
        return (local.ravel()[bins]*weight).sum(axis=2), (bins,weight)

    def layer_deposits(self, stencil, force_on_fluid):
        bins,weight = stencil
        return np.bincount(bins.ravel(),(weight*force_on_fluid[:,:,None]).ravel(),
                           minlength=self.size).reshape(self.B,2,self.ny)

    def step(self, dt, top, stencil, force_on_fluid):
        """Advance old-time fluxes and read two fluid control volumes independently."""
        if not 0 < dt*self.nu/self.h**2 <= .2:
            raise ValueError('Explicit diffusion step violates monotonicity limit')
        if top.shape != (self.B,) or force_on_fluid.shape != (self.B,2*COUNT):
            raise ValueError('Need one boundary drive and one force per site in every copy')
        u = self.u
        mom_z_old = u.sum(axis=2)*self.mass
        du = u[:,1:]-u[:,:-1]
        lap = np.roll(u,-1,axis=2)+np.roll(u,1,axis=2)-2*u
        lap[:,:-1] += du
        lap[:,1:] -= du
        from_top = 2*(top[:,None]-u[:,-1])
        lap[:,-1] += from_top
        self.top_impulse = dt*self.traction_scale*from_top.sum(axis=1)
        self.bottom_impulse = np.zeros(self.B)
        if self.bottom_factor:
            from_bottom = -self.bottom_factor*u[:,0]
            lap[:,0] += from_bottom
            self.bottom_impulse = dt*self.traction_scale*from_bottom.sum(axis=1)
        upper_flux = dt*self.traction_scale*du[:,self.high-1].sum(axis=1)
        lower_flux = dt*self.traction_scale*du[:,self.low-1].sum(axis=1)
        self.lower_shear_impulse = upper_flux-lower_flux
        u += (dt*self.nu/self.h**2)*lap
        deposits = self.layer_deposits(stencil,force_on_fluid)
        for layer,(cells,w) in enumerate(self.z_weights):
            u[:,cells] += (dt/self.mass)*w*deposits[:,layer,None,:]
        dm = u.sum(axis=2)*self.mass-mom_z_old
        self.lower_storage_change = dm[:,self.low:self.high].sum(axis=1)
        lower = self.lower_shear_impulse-self.lower_storage_change
        alternate = upper_flux+self.bottom_impulse-dm[:,:self.high].sum(axis=1)
        residual = dm.sum(axis=1)-self.top_impulse-self.bottom_impulse-dt*force_on_fluid.sum(axis=1)
        return lower, alternate, residual


def check(grid, backing):
    seed = np.random.default_rng(27184)
    nu, rho, back_gap = .8, 1.7, 1.5
    top = np.array([.4,-.7,0.])
    channel = ShearChannel(grid,back_gap,backing,nu,rho,len(top),top,back_impedance=.9)
    y = seed.uniform(-10,10,(len(top),2*COUNT))
    force = seed.normal(size=y.shape)
    zero = np.zeros_like(force)
    _,s = channel.sample(y)
    prior = channel.u.copy()
    dt = .15*channel.h**2/nu
    channel.step(dt,top,s,zero)
    steady_error = float(np.max(abs(prior-channel.u)))
    channel.u[:] = seed.normal(size=channel.u.shape)
    site_u,s = channel.sample(y)
    fluid_bins = channel.layer_deposits(s,force)
    power_grid = np.zeros(len(top))
    for layer,(cells,weights) in enumerate(channel.z_weights):
        power_grid += (channel.u[:,cells]*weights*fluid_bins[:,layer,None,:]).sum(axis=(1,2))
    power_error = float(np.max(abs(power_grid-(site_u*force).sum(axis=1))))
    low,alt,global_resid = channel.step(dt,top,s,force)
    source = -dt*force[:,:COUNT].sum(axis=1)
    cv_error = max(np.max(abs(val-source)) for val in (low,alt))
    kz = np.pi/(TOP+back_gap)*(1. if backing=='rest' else .5)
    ky = 2*np.pi/WIDTH
    vertical = (np.sin if backing=='rest' else np.cos)(kz*(channel.z+back_gap))
    if backing=='robin':
        from scipy.optimize import brentq
        H = TOP+back_gap
        kz = brentq(lambda k: k*channel.slip_length*np.cos(k*H)+np.sin(k*H),np.pi/(2*H),np.pi/H)
        vertical = np.sin(kz*(TOP-channel.z))
    pattern = vertical[:,None]*np.sin(ky*channel.y)
    channel.u[:] = pattern
    steps = int(np.ceil(.4/dt))
    dt = .4/steps
    for _ in range(steps):
        channel.step(dt,np.zeros_like(top),s,zero)
    decay_error = float(np.max(abs(channel.u-pattern*np.exp(-.4*nu*(kz*kz+ky*ky)))))
    emit(test='shear_channel',grid=grid,backing=backing,steady_error=steady_error,
         adjoint_power_max_error=power_error,lower_source=source.tolist(),
         cv_primary=low.tolist(),cv_alternate=alt.tolist(),global_balance=global_resid.tolist(),
         max_cv_error=cv_error,empty_mode_error=decay_error)
    if max(steady_error,power_error,cv_error,np.max(abs(global_resid))) > 1e-11 or decay_error > .003:
        raise ArithmeticError('Channel reciprocity, momentum, or propagation check failed')
    return decay_error


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    for boundary in ('rest','free','robin'):
        errors = [check(n,boundary) for n in (4,8)]
        if errors[1]>=.4*errors[0]:
            raise ArithmeticError('Empty shear mode did not converge under refinement')
