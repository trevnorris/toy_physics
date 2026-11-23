"""Read the excess coherent work of moving the held R4 inlet and first wall cell.

The positive physical displacement opens the rim by one unit. It scales the
inlet disk radially, moves the immediately adjoining wall and vanishes at the
next axial wall node. Two interior extensions in the homogeneous feed have
identical displacements on all physical interfaces and outer faces. The native
straight-vortex and volume elements are evaluated at their displaced positions.
No crystal, neutral face or fixture force is inferred from this partial work.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_joined_front_audit import Front
from mouth_last_rings_fluid import GRIDS, RimCarrier
from mouth_mobile_joined_fields import JointFields
from mouth_rim_pinned_comparison import digest, read, validate_guard
from mouth_supercore_joined_front import ORDER_INERTIA, emit


def taper(z, mode, wall_width):
    # Both extensions vanish at the closed bottom of the feeding slab (-1).
    upstream = z < 0
    p = np.clip(1 + z, 0, 1)
    s = np.clip(1 - z / wall_width, 0, 1)
    return (np.where(upstream, p**mode, s),
            np.where(upstream, mode * p**(mode-1),
                     np.where(z < wall_width, -1 / wall_width, 0)))


def dot(x, y):
    return np.einsum('ij,ij->i', x, y)


def side_work(m, u, height, step):
    # The perturbed part of the two sealed banks is below the polar bump.
    weight_z = np.diff(m.z[1])[0]
    da = np.diff(m.a)
    out = np.zeros(2)
    for qa in (0.5-0.5/np.sqrt(3), 0.5+0.5/np.sqrt(3)):
        angle = m.a[:-1] + qa*da
        sn = np.sin(angle)
        for qz in (0.5-0.5/np.sqrt(3), 0.5+0.5/np.sqrt(3)):
            z = qz * weight_z
            r = m.shape.radius(z)[0] + step*taper(z, 1, height)[0]
            f, df = m._vortex(r*sn)
            area = np.pi * da * weight_z * sn**2 * r**2 * f * df
            for bank in (0, 1):
                nodes = m.ids[1,bank][0 if bank else -1, :, :2]
                v = (1-qz)*u[0, nodes[:,0]] + qz*u[0, nodes[:,1]]
                v = (1-qa)*v[:-1] + qa*v[1:]
                out[bank] += (1-2*bank)*.5*(area @ (v*v-1))
    return out


class CarrierWork:
    def __init__(self, m, u, height):
        self.m, self.u, self.height = m, u, height
        take = m.D[:,3] < height
        wt, r, a, z, *_ = m.D[take].T
        b, bz = m.shape.radius(z)
        self.bank = np.where(z<0, 0, 2)+(r>b).astype(int)
        self.wt,self.r,self.a,self.z = wt,r,a,z
        self.b,self.bz = b,bz
        self.v = m.V[take] @ u[0]
        if np.any(m.V[take] @ u[2]) or np.any(m.D[take,8]):
            raise ArithmeticError('Entrance lift intersects mobile conversion order')
        gradients = [matrix[take] for matrix in m.G]
        self.gv,self.gp = [np.column_stack([g @ v for g in gradients]) for v in u[:2]]
        if any(np.any(g @ u[2]) for g in gradients):
            raise ArithmeticError('Mobile-order gradient reaches the aperture probe')
        self.gc = m.ref[take,1:]

    def energy(self, mode, step):
        r,a,z,b = self.r,self.a,self.z,self.b
        f, df = taper(z, mode, self.height)
        inside = self.bank%2==0
        scale = np.where(inside, r/b, (28-r)/(28-b))
        scale_r = np.where(inside, 1/b, -1/(28-b))
        moved = r + step*f*scale
        radius_jac = 1 + step*f*scale_r
        shear = step*(df*scale - f*scale*scale_r*self.bz)
        if np.any(radius_jac<=0) or np.any(moved<=0):
            raise ArithmeticError('Signed virtual map reverses a sampled element')
        def transport(gradient):
            radial = gradient[:,0]/radius_jac
            return np.column_stack((radial, gradient[:,1]*r/moved,
                                    gradient[:,2]-shear*radial))
        gv,gp,gc = map(transport, (self.gv,self.gp,self.gc))
        vortex = self.m._vortex(moved*np.sin(a))[0]
        y = self.v**2
        n = vortex**2*y
        pi = ORDER_INERTIA*dot(gp,gc)/(n+ORDER_INERTIA*dot(gc,gc))
        currents = n[:,None]*(gp-pi[:,None]*gc)
        density = (.5*vortex**2*dot(gv,gv)
                   + vortex**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20
                   + .5*dot(gp,currents))
        weights = self.wt*(moved/r)**2*radius_jac
        volume = np.bincount(self.bank, weights=weights*density, minlength=4)
        return np.r_[volume, side_work(self.m,self.u,self.height,step)]


class SignedChart:
    """Move only the spatial mesh of the carrier. Keep the registered R4 label."""
    def __init__(self, center, height, step):
        self.center,self.height,self.step = center,height,step
        self.neck,self.length = center.neck,center.length

    def radius(self, z):
        r, slope = self.center.radius(z)
        f, df = taper(np.asarray(z), 1, self.height)
        return r+self.step*f, slope+self.step*df


def check_native(f, pair, height, step, work):
    moved = SignedChart(f.shape, height, step)
    def profile(x, z):
        return f.profile(x,z)
    def outside(r, z, above):
        b = f.shape.radius(z)[0]
        new_b = moved.radius(z)[0]
        return f.profile.outside(b+(r-new_b)*(28-b)/(28-new_b), z, above)
    profile.outside = outside
    m = RimCarrier(moved, f.order, 1, f.grid, precondition=False,
                   polar_shift=f.polar_shift, terminal_profile=profile)
    if any(not np.array_equal(ids, m.ids[key]) for key,ids in f.model.ids.items()):
        raise ArithmeticError('Signed carrier mesh changed topology')
    # Only the derivative is needed. The remote material coupling is independent
    # of this lift and is zero in both native carrier-only evaluations.
    native = np.array([m.evaluate(u,read=True)[0]-f.model.evaluate(u,read=True)[0]
                       for u in pair])
    local = np.array([(item.energy(1,step)-item.energy(1,0.)).sum() for item in work])
    emit(kind='native_energy_check', step=step, signed_change=native.tolist(),
         local_change=local.tolist(), error=(native-local).tolist())
    if max(abs(native-local))>3e-9:
        raise ArithmeticError('Signed local work disagrees with native assembled energy')


def run(endpoint):
    validate_guard(endpoint)
    command = read(endpoint/'invocation.json')['command']
    script = next(i for i,word in enumerate(command) if Path(word).name=='mouth_whole_r4_relax.py')
    state = Path(command[script+1])
    rec = read(endpoint/'endpoint.json')
    saved = endpoint/rec['file']
    if digest(saved)!=rec['sha256']:
        raise ValueError('Saved field differs from endpoint')
    with np.load(saved) as s, np.load(endpoint/'endpoint_controls.npz') as c:
        if s['state_sha'].item()!=digest(state):
            raise ValueError('Unmatched reference preparation')
        x,pair = s['x'], (s['field'],c['rest'])
    base=Front(state,GRIDS['double'],1)
    f=JointFields(base,state,x,grid='double',precondition=False)
    m=f.model
    for u,bias in zip(pair,(BIAS,0.)):
        if np.any(u[1,m.feed]) or np.any(u[1,m.bulk]!=bias):
            raise ValueError('Incorrect held phase in saved pair')
    height=m.z[1][1]
    if height>=min(.06*m.shape.length,m.shape.length-m.h):
        raise ArithmeticError('First wall cell leaves the unshifted inlet neighborhood')
    emit(kind='inputs', mesh=endpoint.name, saved=digest(saved),
         controls=digest(endpoint/'endpoint_controls.npz'), wall_taper_height=height,
         bank_order=['feed_inner','feed_outer','shaft','receiving_beside_shaft',
                     'inner_sealed_wall','outer_sealed_wall'])
    work=[CarrierWork(m,u,height) for u in pair]
    check_native(f,pair,height,1e-4,work)
    for mode in (1,2):
        centers=np.array([w.energy(mode,0.) for w in work])
        for step in (1e-4,5e-5):
            minus=np.array([w.energy(mode,-step) for w in work])-centers
            plus=np.array([w.energy(mode,step) for w in work])-centers
            slopes=(plus-minus)/(2*step)
            emit(kind='signed_work', feed_extension_power=mode, step=step,
                 minus=minus.sum(axis=1).tolist(), plus=plus.sum(axis=1).tolist(),
                 driven_rest_by_term=slopes.tolist(), driven_rest=slopes.sum(axis=1).tolist(),
                 driven_minus_rest=float((slopes[0]-slopes[1]).sum()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('endpoint',type=Path)
    run(p.parse_args().endpoint)
