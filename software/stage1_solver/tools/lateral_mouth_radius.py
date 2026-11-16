"""Held load and brane-side translation reaction for equal mean intake speed."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

from lateral_mouth_acceleration import Translation
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import check, verified_packet
from mouth_feed_weak_check import dx4


class Target:
    def __init__(self, slab, radius, unit_density, mean_speed):
        self.slab = slab
        self.radius = radius
        self.volume = 4*np.pi*radius**3/3
        self.mean_speed = mean_speed
        # Even continuation across the radial axis for scaled sample positions.
        self.density = RectBivariateSpline(
            slab.x, np.r_[-slab.r[::-1], slab.r], np.column_stack((unit_density[:, ::-1], unit_density)))
        self.g, self.speed, self.dg, self.dspeed = self.profile(0)

    def profile(self, offset):
        s, a = self.slab, self.radius
        x = (s.x[:, None]-offset)/a
        r = s.r[None, :]/a
        cut = np.maximum(1-x*x-r*r, 0)
        g = cut**3
        dg = 6*x*cut**2/a
        # Only the unit ball contributes. Clipping avoids bath extrapolation.
        xx, rr = np.clip(x, -1, 1), np.minimum(r, 1)
        ref = self.density.ev(xx, rr)
        shape = g/ref
        dshape = (dg + shape*self.density.ev(xx, rr, dx=1)/a)/ref
        norm = np.sum(shape*s.dv)
        speed = self.mean_speed*self.volume*shape/norm
        dspeed = (self.mean_speed*self.volume*dshape-speed*np.sum(dshape*s.dv))/norm
        return g, speed, dg, dspeed


def steady(slab, target, distance, unit_density):
    restoring = target.g.copy()
    if distance:
        src = np.maximum(1-(slab.x[:, None]+distance)**2-slab.r**2, 0)**3
        shift = round(-distance/slab.dx)
        if shift*slab.dx != -distance or np.any((src > 0) & (target.g > 0)):
            raise ValueError("Unsupported companion geometry")
        source_shape = src/np.roll(unit_density, shift, axis=0)
        restoring += src

    def boundary(face):
        n = abs(face)**2
        currents = [n*target.speed]
        if distance:
            raw = source_shape*n
            currents.insert(0, 0.25*raw/np.sum(slab.dv*raw))
        normal = (-restoring*(1-1/abs(face)) + 1j*sum(currents)/n)*face
        return np.array(currents), normal

    field = np.ones(slab.eigenvalues.shape, complex)
    for _ in range(500):
        _, normal = boundary(field[-1])
        potential = (abs(field)**8-1)/4
        re = slab.elliptic(2*(field.real-1)-potential*field.real, normal.real, 2)
        im = slab.elliptic(-potential*field.imag, normal.imag, 0)
        updated = 1+re+1j*im
        error = np.max(abs(updated-field))
        if not np.isfinite(error):
            raise ArithmeticError("Nonfinite steady update")
        field = updated
        if error < 5e-13:
            currents, _ = boundary(field[-1])
            return field, currents
    raise ArithmeticError("Steady solve did not converge")


def held_reading(slab, target, field, currents, distance):
    face = field[-1]
    n = abs(face)**2
    j = currents[-1]
    actual_speed = j/n
    avg = np.sum(slab.dv*actual_speed)/target.volume
    if (np.max(abs(actual_speed-target.speed)) > 1e-13 or abs(avg-target.mean_speed) > 1e-13
            or np.any(j[target.g == 0] != 0)):
        raise ArithmeticError("Target intake speed mismatch")
    rate = np.sum(j*slab.dv)
    rates = np.array([0.25,rate] if distance else [rate])
    g = target.g
    if distance:
        g = g + np.maximum(1-(slab.x[:, None]+distance)**2-slab.r**2, 0)**3
    normal = (-g*(1-1/abs(face)) + 1j*currents.sum(axis=0)/n)*face
    residual, bath, _ = check(slab, field, normal, currents, rates)
    dxface = dx4(face,slab.dx)
    carried = -np.sum(slab.dv*j*(face.conj()*dxface).imag/n)
    slope = -2*target.g*(1-1/abs(face))
    stress = -np.sum(slab.dv*slope*dx4(n,slab.dx))/4
    def derivative(f):
        return (11*f[-1]-18*f[-2]+9*f[-3]-2*f[-4])/(6*slab.dz)
    mixed = (dxface.conj()*derivative(field)).real-dx4(derivative(abs(field[-4:])**2),slab.dx)/4
    direct = -np.sum((slab.dv*mixed)[target.g>0])
    parts = np.array([carried,stress,direct])
    if not np.isfinite(parts).all():
        raise ArithmeticError("Nonfinite held load")
    print(f"state_d={distance} target_rate={rate:.12g} mean_speed={avg:.12g} "
          f"bath={bath:.12g} residual={residual:.5g}", flush=True)
    return parts


class SpeedTranslation(Translation):
    def __init__(self, target, base, current):
        self.slab = target.slab
        self.f = base
        self.face = base[-1]
        self.amp = abs(self.face)
        self.n = abs(base)**2
        self.j = current
        self.rate = np.sum(self.slab.dv*current)
        self.g, self.speed = target.g, target.speed
        self.normal_factor = -self.g*(1-1/self.amp) + 1j*self.speed
        self.potential = (self.n**4-1)/4
        self.displaced_normal = (-target.dg*(1-1/self.amp)+1j*target.dspeed)*self.face
        def fixed_trace(offset):
            g, speed, _, _ = target.profile(offset)
            return (-g*(1-1/self.amp)+1j*speed)*self.face
        eps = 1e-4*self.slab.dx
        error = np.max(abs((fixed_trace(eps)-fixed_trace(-eps))/(2*eps)-self.displaced_normal))
        if not np.isfinite(error) or error > 3e-8:
            raise ArithmeticError(f"Translation slope mismatch: {error}")

    def boundary(self, h):
        dn = 2*(self.face.conj()*h[-1]).real
        return self.normal_factor*h[-1]-self.face*self.g*dn/(2*self.amp**3)

    def load(self, h):
        # At fixed speed dj = speed*dn. Restore the integrated-rate term
        # subtracted by the fixed-rate readout. It vanishes for a centered dipole.
        carry, stress, direct = super().load(h)
        dn = 2*(self.face.conj()*h[-1]).real
        dj_sum = np.sum(self.slab.dv*self.speed*dn)
        vx = (self.face.conj()*dx4(self.face,self.slab.dx)).imag/self.n[-1]
        carry += dj_sum/self.rate*np.sum(self.slab.dv*self.j*vx)
        if abs(dj_sum) > 1e-10:
            raise ArithmeticError("Centered acceleration changes the integrated rate at first order")
        return carry, stress, direct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("radius", type=float, choices=(0.75, 1., 1.5))
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(k):v for run in archive["runs"] for k,v in run.get("fields",{}).items()}
    with verified_packet(args.packet, hashes) as packet:
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
        slab = LateralSlab(round(r[-1]+r[0]), round(1/(2*r[0])), round(1/dz))
        if (int(packet["d"]) != 0 or not np.array_equal(x,slab.x)
                or not np.array_equal(r,slab.r) or slab.dz != dz):
            raise ValueError("Calibration grid mismatch")
        unit = packet["field"]
        unit_j = packet["currents"]
    unit_density = abs(unit[-1])**2
    mean = np.sum(slab.dv*unit_j[0]/unit_density)/(4*np.pi/3)
    target = Target(slab, args.radius, unit_density, mean)
    print(f"radius={args.radius:g} extent={r[-1]+r[0]:g} dx={slab.dx:.12g} dz={dz:.12g} "
          f"calibrated_speed={mean:.12g}", flush=True)
    if args.radius == 1:
        base, currents = unit, unit_j
    else:
        del unit
        base, currents = steady(slab, target, 0, unit_density)
    alone = held_reading(slab, target, base, currents, 0)
    if max(abs(alone)) > 1e-11:
        raise ArithmeticError("Isolated mouth has a lateral load")
    translation = SpeedTranslation(target, base, currents[-1])
    p, sr = translation.static_tangent()
    reactions = {}
    for omega in (0., 0.01):
        h, residual, steps = translation.dynamic_tangents(p,omega)
        carry, stress, direct = translation.load(h)
        reactions[omega] = -(carry+stress)
        print(f"omega={omega:g} carry_per_A={carry:.12g} stress_per_A={stress:.12g} "
              f"reaction={reactions[omega]:.12g} direct={-direct:.12g} static_residual={sr:.5g} "
              f"dynamic_residual={residual:.5g} steps={steps}", flush=True)
    del translation, p, h, base
    for distance in (4,8):
        field, currents = steady(slab, target, distance, unit_density)
        carry, stress, direct = held_reading(slab, target, field, currents, distance)-alone
        total = carry+stress
        print(f"d={distance} carry={carry:.12g} stress={stress:.12g} load={total:.12g} "
              f"direct={direct:.12g}", flush=True)
        if reactions[0.] != 0:
            print(f"ratio_d={distance} load_per_slow_reaction={total/reactions[0.]:.12g}", flush=True)
        del field


if __name__ == "__main__":
    main()
