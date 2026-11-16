"""One-mouth brane load at the onset of prescribed translation from rest."""

import json
from pathlib import Path
import shlex
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import CubicSpline

from lateral_mouth_acceleration import Translation
from lateral_paired_acceleration import target_load
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


def carrier_rate(motion, perturbation, shift):
    """Lab-frame time derivative from the real-space linearized carrier equation."""
    slab = motion.slab
    z = perturbation
    axial = np.concatenate((-z[:, :1], z, -z[:, -1:]), axis=1)
    lap = np.diff(axial, n=2, axis=1)/slab.dx**2
    edges = np.concatenate((z[:, :, :1], z, -z[:, :, -1:]), axis=2)
    faces = np.arange(len(slab.r)+1)*slab.dx
    lap += np.diff(faces*np.diff(edges, axis=2), axis=2)/(slab.r*slab.dx**2)
    normal = motion.boundary(z) + shift*motion.displaced_normal
    vertical = np.concatenate((z[1:2], z, (z[-2]+2*slab.dz*normal)[None]), axis=0)
    lap += np.diff(vertical, n=2, axis=0)/slab.dz**2
    return -1j*(-lap/2 + motion.potential_change(z))


def frozen_load(motion, density_spline, shift):
    slab = motion.slab
    g = np.maximum(1-(slab.x[:, None]-shift)**2-slab.r**2, 0)**3
    shape = g/density_spline(slab.x-shift)
    speed = motion.rate*shape/np.sum(slab.dv*shape*motion.n[-1])
    return target_load(motion.face, g, speed, slab)


def onset(motion, slope, intervals, duration):
    # All perturbations are per unit imposed acceleration. X=t²/2 at release.
    z = np.zeros_like(motion.f)
    dt = duration/intervals
    readings = []
    for step in range(intervals):
        t = step*dt
        k1 = carrier_rate(motion, z, t*t/2)
        half_shift = (t+dt/2)**2/2
        k2 = carrier_rate(motion, z+dt*k1/2, half_shift)
        k3 = carrier_rate(motion, z+dt*k2/2, half_shift)
        k4 = carrier_rate(motion, z+dt*k3, (t+dt)**2/2)
        z += dt*(k1+2*k2+2*k3+k4)/6
        if (step+1) in (intervals//4, intervals//2, intervals):
            time = (step+1)*dt
            carried, stress, _ = motion.load(z[-4:])
            load = slope*time**2/2 + carried + stress
            readings.append(dict(t=time, incremental_load_per_A=load, twice_load_over_t2=2*load/time**2))
    return readings


def main():
    archive = json.loads(Path("docs/_measurements/lateral_mouth_acceleration.json").read_bytes())
    for run in archive["runs"]:
        path = Path(shlex.split(run["command"])[-1])
        if not run["limits_verified"] or run["guard_exit"] or run["worker_exit"]:
            raise ValueError("Unverified motion reference")
        with verified_packet(path, {path: run["packet_sha256"]}) as packet:
            x, r = packet["x"], packet["r"]
            if int(packet["d"]) != 0 or packet["currents"].shape[0] != 1:
                raise ValueError("Expected one isolated mouth")
            dx = x[1]-x[0]
            slab = SimpleNamespace(x=x, r=r, dx=dx, dz=float(packet["dz"]), dv=2*np.pi*r*dx**2)
            motion = Translation(slab, packet["field"], packet["currents"][0])
        spline = CubicSpline(x, motion.n[-1], axis=0)
        held = frozen_load(motion, spline, 0)
        dg = 6*x[:, None]*np.maximum(1-x[:, None]**2-r**2, 0)**2
        dspeed = (motion.displaced_normal/motion.face).imag
        advective = np.sum(slab.dv*dspeed*(motion.face.conj()*dx4(motion.face, dx)).imag)
        stress = -np.sum(slab.dv*dg*(1-1/motion.amp)*dx4(motion.n[-1], dx))/2
        slope = advective + stress
        shifts = (1e-4*dx, 5e-5*dx)
        fd = [(frozen_load(motion, spline, h)-frozen_load(motion, spline, -h))/(2*h) for h in shifts]
        if abs(held)>1e-11 or not np.allclose(fd, slope, rtol=2e-7, atol=1e-13):
            raise ArithmeticError("Isolated load or translated fixed-trace derivative failed")
        duration = 0.001
        coarse = onset(motion, slope, 4, duration)
        fine = onset(motion, slope, 8, duration)
        values = [v["incremental_load_per_A"] for v in fine]
        if not np.allclose(values, [v["incremental_load_per_A"] for v in coarse], rtol=2e-6, atol=1e-17):
            raise ArithmeticError("Time-step refinement failed on the prescribed-motion test")
        # Cancel the O(t) term in 2 ΔL/t² and compare its onset with the analytic slope.
        extrapolated = 2*fine[0]["twice_load_over_t2"]-fine[1]["twice_load_over_t2"]
        if not np.isclose(extrapolated, slope, rtol=3e-3, atol=1e-12):
            raise ArithmeticError("Carrier initial-value solution disagrees with its onset expansion")
        slow, = (line for line in run["stdout"] if line.startswith("omega=0 "))
        parts = dict(token.split("=") for token in slow.split())
        print(json.dumps(dict(input=str(path), L=r[-1]+r[0], dx=dx, dz=slab.dz, rate=motion.rate,
                              isolated_load=held, frozen_slope=slope, slope_parts=[advective, stress],
                              centered_slopes=fd, short_time=fine, extrapolated_slope=extrapolated,
                              max_step_change=max(abs(a["incremental_load_per_A"]-b["incremental_load_per_A"])
                                                  for a,b in zip(coarse,fine)),
                              archived_slow_reaction=float(parts["compensating_per_A"])),
                         allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
