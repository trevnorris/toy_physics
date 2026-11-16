"""Brane reaction to prescribed slow axial translation in the fixed lateral reservoir.

Write the laboratory-frame perturbation as p X + u V + h A. For a harmonic
translation at frequency omega, Lp=0, Lu=i(p-omega^2 h), Lh=i u. Only p
has an inhomogeneous mouth condition. This representation subtracts the
static displaced shape before measuring the acceleration response and also
gives the finite-domain zero-frequency derivative without cancellation.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from lateral_reservoir_field import LateralSlab
from mouth_feed_weak_check import dx4


class Translation:
    def __init__(self, slab, base, current):
        self.slab = slab
        self.f = base
        self.n = abs(base)**2
        self.j = current
        self.face = base[-1]
        self.amp = abs(self.face)
        self.rate = np.sum(current*slab.dv)
        x = slab.x[:, None]
        g0 = np.maximum(1-x*x-slab.r**2, 0)
        self.g = g0**3
        self.speed = current/self.n[-1]
        self.normal_factor = -self.g*(1-1/self.amp) + 1j*self.speed
        self.potential = (self.n**4-1)/4

        # Translate the calibrated normal-speed shape g(x-X)/n0(x-X).
        # Interpolation is used only for the isolated density in this shape.
        interp = CubicSpline(slab.x, self.n[-1], axis=0)
        dg = 6*x*g0**2
        shape = self.g/self.n[-1]
        dshape = (dg + shape*interp(slab.x, 1))/self.n[-1]
        norm = np.sum(self.g*slab.dv)
        dv = (self.rate*dshape - self.speed*np.sum(dshape*self.n[-1]*slab.dv))/norm
        self.displaced_normal = (-dg*(1-1/self.amp)+1j*dv)*self.face

        def fixed_trace_normal(offset):
            g = np.maximum(1-(x-offset)**2-slab.r**2, 0)**3
            s = g/interp(slab.x-offset)
            v = self.rate*s/np.sum(s*self.n[-1]*slab.dv)
            return (-g*(1-1/self.amp)+1j*v)*self.face
        eps = 1e-4*slab.dx
        diff = (fixed_trace_normal(eps)-fixed_trace_normal(-eps))/(2*eps)
        if (abs(self.rate-0.25)>1e-12 or
                np.max(abs(diff-self.displaced_normal))>3e-8):
            raise ArithmeticError("Static translation of the fixed-rate mouth law failed")

    def boundary(self, h):
        dh = h[-1]
        dn = 2*(self.face.conj()*dh).real
        mean = np.sum(self.slab.dv*self.speed*dn)/self.rate
        return self.normal_factor*dh + self.face*(-self.g*dn/(2*self.amp**3)-1j*self.speed*mean)

    def potential_change(self, h):
        dn = 2*(self.f.conj()*h).real
        return self.potential*h + self.n**3*dn*self.f

    def step(self, h, rhs, boundary_drive=0):
        b = self.boundary(h) + boundary_drive
        work = rhs-self.potential_change(h)
        return (self.slab.elliptic(2*h.real+work.real,b.real,2)
                + 1j*self.slab.elliptic(work.imag.copy(),b.imag,0))

    def residual(self, h, rhs, boundary_drive=0):
        # Use real-space cylindrical differences and homogeneous bath ghosts
        # independently of the separated elliptic solver.
        s = self.slab
        edge = np.concatenate((-h[:, :1],h,-h[:, -1:]), axis=1)
        lap = np.diff(edge,n=2,axis=1)/s.dx**2
        edge = np.concatenate((h[:, :, :1],h,-h[:, :, -1:]),axis=2)
        radii = np.arange(len(s.r)+1)*s.dx
        lap += np.diff(radii*np.diff(edge,axis=2),axis=2)/(s.r*s.dx**2)
        normal = self.boundary(h)+boundary_drive
        edge = np.concatenate((h[1:2],h,(h[-2]+2*s.dz*normal)[None]),axis=0)
        lap += np.diff(edge,n=2,axis=0)/s.dz**2
        return float(np.max(abs(-lap/2+self.potential_change(h)-rhs)))

    def static_tangent(self):
        h = np.zeros_like(self.f)
        for _ in range(500):
            new = self.step(h,0,self.displaced_normal)
            error = np.max(abs(new-h))
            h = new
            if not np.isfinite(error):
                raise ArithmeticError("Nonfinite static tangent")
            if error < 2e-13:
                residual = self.residual(h,0,self.displaced_normal)
                if residual > 3e-9:
                    raise ArithmeticError(f"Static residual {residual}")
                return h, residual
        raise ArithmeticError("Static translation did not converge")

    def dynamic_tangents(self, p, omega):
        u, h = np.zeros_like(p), np.zeros_like(p)
        for step in range(1,501):
            new_u = self.step(u,1j*(p-omega**2*h))
            new_h = self.step(h,1j*new_u)
            error = max(np.max(abs(new_h-h)),np.max(abs(new_u-u)))
            u, h = new_u, new_h
            if not np.isfinite(error):
                raise ArithmeticError("Nonfinite dynamic tangent")
            if error < 2e-12:
                residual = max(self.residual(u,1j*(p-omega**2*h)),self.residual(h,1j*u))
                if residual > 3e-9:
                    raise ArithmeticError(f"Dynamic residual {residual}")
                return h, residual, step
        raise ArithmeticError(f"Dynamic tangents did not converge at omega={omega}")

    def load(self, h):
        # Positive x is the direction of imposed acceleration. Compare at the
        # same position against the static solution, so geometric cut shifts
        # cancel. A tangentially moving opening has zero normal cut speed.
        s = self.slab
        n, face, dh = self.n[-1], self.face, h[-1]
        dn = 2*(face.conj()*dh).real
        dj = self.j*(dn/n - np.sum(self.speed*dn*s.dv)/self.rate)
        dx_face, dx_h = dx4(face,s.dx), dx4(dh,s.dx)
        vx = (face.conj()*dx_face).imag/n
        dvx = ((dh.conj()*dx_face+face.conj()*dx_h).imag-dn*vx)/n
        carry = np.sum(s.dv*(dj*vx+self.j*dvx))
        slope = -2*self.g*(1-1/self.amp)
        dslope = -self.g*dn/self.amp**3
        stress = np.sum(s.dv*(dslope*dx4(n,s.dx)+slope*dx4(dn,s.dx)))/4
        # Second reading from interior fields, without the prescribed normal
        # derivatives or the integration-by-parts stress identity.
        def normal(f):
            return (11*f[-1]-18*f[-2]+9*f[-3]-2*f[-4])/(6*s.dz)
        dn_layers = 2*(self.f[-4:].conj()*h[-4:]).real
        mixed = ((dx_h.conj()*normal(self.f) + dx_face.conj()*normal(h)).real
                 - dx4(normal(dn_layers),s.dx)/4)
        direct = np.sum((s.dv*mixed)[self.g>0])
        if abs(np.sum(dj*s.dv))>1e-11 or not np.isfinite([carry,stress,direct]).all():
            raise ArithmeticError("Invalid integrated momentum flux or changing intake")
        return carry, stress, direct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet",type=Path)
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(k):v for run in archive["runs"] for k,v in run.get("fields",{}).items()}
    if sha256(args.packet.read_bytes()).hexdigest()!=hashes[args.packet]:
        raise ValueError("Isolated source packet mismatch")
    with np.load(args.packet,allow_pickle=False) as data:
        x,r,dz = data["x"],data["r"],float(data["dz"])
        slab = LateralSlab(round(r[-1]+r[0]),round(1/(2*r[0])),round(1/dz))
        if (int(data["d"])!=0 or not np.array_equal(x,slab.x)
                or not np.array_equal(r,slab.r) or slab.dz!=dz):
            raise ValueError("Geometry differs from the isolated preparation")
        experiment = Translation(slab,data["field"],data["currents"][0])
    p, static_residual = experiment.static_tangent()
    print(f"extent={r[-1]+r[0]:g} dx={slab.dx:.8g} dz={dz:.8g} rate={experiment.rate:.12g} "
          f"static_residual={static_residual:.3g}",flush=True)
    for omega in (0., 0.005, 0.01):
        h, residual, steps = experiment.dynamic_tangents(p,omega)
        carry, stress, direct = experiment.load(h)
        print(f"omega={omega:g} carried_per_A={carry:.12g} stress_per_A={stress:.12g} "
              f"load_per_A={carry+stress:.12g} compensating_per_A={-carry-stress:.12g} "
              f"direct_compensating_per_A={-direct:.12g} residual={residual:.3g} steps={steps}",
              flush=True)


if __name__ == "__main__":
    main()
