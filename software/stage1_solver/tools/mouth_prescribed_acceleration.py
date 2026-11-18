"""Prescribe X=A*t^2/2 on the saved resting drain and read the feeding-side load."""

import argparse
import json

import numpy as np

from held_crossflow_budget import archived_packet
from mouth_feed_weak_flow import FeedPatch


class MovingFeed:
    def __init__(self):
        path = archived_packet("mouth_feed_weak_run.json", "n32_restoring_slope_u0.npz")
        self.patch = FeedPatch(32, 65, "restoring_slope")
        p = self.patch
        with np.load(path, allow_pickle=False) as data:
            self.base = data["field"]
            if (self.base.shape != p.shape or float(data["speed"]) != 0
                    or (float(data["dx"]), float(data["dz"])) != (p.dx, p.dz)
                    or not np.array_equal(data["current"], p.current)
                    or not np.array_equal(data["mouth"], p.mouth)):
                raise ValueError("Mismatched resting drain")
        self.rate = p.dx**3*p.current.sum()
        nref = np.abs(self.base[-1])**2
        self.reference_spectrum = np.fft.rfft(nref, axis=0)
        self.k = 2*np.pi*np.fft.rfftfreq(p.n, p.dx)[:, None, None]
        self.transverse_squared = p.x.reshape(1,-1,1)**2 + p.x.reshape(1,1,-1)**2

    def boundary(self, face, shift):
        p = self.patch
        # Shift the radius-one opening and its reference density together in
        # the periodic lab coordinates. The lab reservoir remains at rest.
        r2 = (p.x-shift)**2 + self.transverse_squared
        p.mouth = r2 < 1
        p.restoring = np.maximum(1-r2, 0)**3
        shifted_ref = np.fft.irfft(self.reference_spectrum*np.exp(-1j*self.k*shift),
                                   n=p.n, axis=0)
        amp = np.abs(face)
        n = amp**2
        raw = p.restoring*n/shifted_ref
        p.current = self.rate*raw/(p.dx**3*raw.sum())
        return (-p.restoring*(1-1/amp) + 1j*p.current/n)*face

    def rhs(self, field, shift):
        p = self.patch
        normal = self.boundary(field[-1], shift)
        lap = np.zeros_like(field)
        lap[1:-1] = (field[2:] - 2*field[1:-1] + field[:-2])/p.dz**2
        lap[-1] = 2*(field[-2] - field[-1] + p.dz*normal)/p.dz**2
        interior = field[1:]
        for axis in (1,2,3):
            lap[1:] += (np.roll(interior, 1, axis) - 2*interior
                        + np.roll(interior, -1, axis))/p.dx**2
        lap[1:] = 1j*(lap[1:]/2 - p.mu0*(np.abs(interior)**8-1)*interior)
        return lap  # fixed reservoir gives zero time derivative at plane zero

    def reading(self, field, shift):
        self.boundary(field[-1], shift)
        p = self.patch
        read = p.read(field, 0.0)
        imposed_rate = p.dx**3*p.current.sum()
        if not np.isclose(imposed_rate, self.rate, atol=1e-13, rtol=0):
            raise ArithmeticError("Translated intake lost its normalization")
        return dict(carried=read["carried"], stress=read["stress"],
                    feed_load=read["carried"]+read["stress"], inward_rate=imposed_rate)

    def drive(self, acceleration, intervals):
        field = self.base.copy()
        dt = 0.001/intervals
        rows = []
        for i in range(intervals):
            t = i*dt
            start = acceleration*t**2/2
            middle = acceleration*(t+dt/2)**2/2
            end = acceleration*(t+dt)**2/2
            k1 = self.rhs(field, start)
            k2 = self.rhs(field + dt*k1/2, middle)
            k3 = self.rhs(field + dt*k2/2, middle)
            k4 = self.rhs(field + dt*k3, end)
            field += dt*(k1 + 2*k2 + 2*k3 + k4)/6
            if i+1 in (intervals//4, intervals//2, intervals):
                rows.append(dict(t=t+dt, x=end, **self.reading(field, end)))
        return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("acceleration", type=float, choices=(-0.02,-0.01,0,0.01,0.02))
    acceleration = parser.parse_args().acceleration
    model = MovingFeed()
    rest = model.reading(model.base, 0)
    residual = np.max(np.abs(model.rhs(model.base, 0)))
    if residual > 2e-8 or abs(rest["feed_load"]) > 1e-13:
        raise ArithmeticError("Resting feed or stationary residual failed")
    result = dict(acceleration=acceleration, grid=list(model.base.shape),
                  rest=rest, stationary_residual=residual)
    for steps in (16,32):
        result[str(steps)] = model.drive(acceleration, steps)
    print(json.dumps(result, allow_nan=False), flush=True)
    coarse = [r["feed_load"] for r in result["16"]]
    fine = [r["feed_load"] for r in result["32"]]
    np.testing.assert_allclose(coarse, fine, rtol=3e-5, atol=3e-17)


if __name__ == "__main__":
    main()
