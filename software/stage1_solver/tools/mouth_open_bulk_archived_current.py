"""Paired current traces as Neumann data for the same ideal 4D receiver.

Reconstruct the archive's calibrated taper times the paired/alone density
ratio on each outlet. Interpolate the density ratio and normalize each
continuous trace to its archived integral. Phi and tangential velocity on
the feed are free on this receiver side.
"""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline

from incident_flow_falloff import check_current, checked_grid
from lateral_source_rate import verified_packet
from mouth_front_pressure_dipole import grid
from mouth_open_bulk_pressure_moment import FEED, PEAK, RATE, exterior_series, interior_rays, unit_gauss


def ball(order):
    r, wr = unit_gauss(order)
    mu, wm = leggauss(order)
    weights = 2 * np.pi * (wr * r**2)[:, None] * wm
    return r[:, None], mu, weights


class Trace:
    def __init__(self, x, a, excess, order):
        # x points away from increasing archive x, towards the companion.
        # Reflect the axisymmetric scalar density through the radial axis.
        self.spline = RectBivariateSpline(-x[::-1], np.r_[-a[::-1], a],
                                          np.c_[excess[::-1, ::-1], excess[::-1]],
                                          kx=3, ky=3, s=0)
        r, mu, weights = ball(order)
        h = self.spline.ev(r * mu, r * np.sqrt(1 - mu**2))
        # Normalize the excess relative to the exact, rate-0.25 radial base.
        self.mean = np.sum(weights * PEAK * (1 - r**2)**3 * h) / RATE

    def perturbation(self, x, a, gradient=False):
        taper = np.maximum(1 - x**2 - a**2, 0)
        h = self.spline.ev(x, a) - self.mean
        peak = PEAK / (1 + self.mean)
        if not gradient:
            return peak * taper**3 * h
        gx = taper * self.spline.ev(x, a, dx=1) - 6 * x * h
        ga = taper * self.spline.ev(x, a, dy=1) - 6 * a * h
        return peak*taper**2*gx, peak*taper**2*ga


def load(cells, d, hashes):
    scratch = Path("_scratch/s11c")
    s = grid(96, cells, 4*cells)
    alone_path = scratch / "incident_flow_falloff01" / f"L96C{cells}_baseline.npz"
    pair_path = scratch / "held_force_incident01" / f"L96C{cells}_d{d}.npz"
    with verified_packet(alone_path, hashes) as packet:
        checked_grid(s, packet)
        assert int(packet["d"]) == 0
        n0 = abs(packet["layers"][-1])**2
        j0 = check_current(s, 0, packet["currents"])
        shape = np.maximum(1 - s.x[:, None]**2 - s.r**2, 0)**3
        np.testing.assert_allclose(j0, RATE*shape/np.sum(s.dv*shape), rtol=1e-12, atol=1e-15)
    with verified_packet(pair_path, hashes) as packet:
        checked_grid(s, packet)
        np.testing.assert_array_equal(packet["centers"], [-d, 0])
        np.testing.assert_array_equal(packet["rates"], [RATE, RATE])
        n = abs(packet["layers"][-1])**2
        currents = packet["currents"]
    traces, differences = [], []
    for index, center in enumerate((-d, 0)):
        current = check_current(s, center, currents[index:index+1])
        offset = round(center/s.dx)
        weight = np.maximum(1 - (s.x[:, None] - center)**2 - s.r**2, 0)**3
        excess = n / np.roll(n0, offset, axis=0) - 1
        raw = weight*(1 + excess)
        reconstructed = RATE * raw / np.sum(s.dv*raw)
        differences.append(float(np.max(abs(reconstructed - current))))
        # Keep extra density samples around the ball for the cubic spline.
        crop = np.arange(-cells-4, cells+4) + len(s.x)//2 + offset
        traces.append(Trace(s.x[crop] - center, s.r[:cells+4], excess[crop, :cells+4], 72))
    print(json.dumps(dict(cells=cells, d=d, archived_current_max_abs_errors=differences)), flush=True)
    assert max(differences) < 2e-14
    return traces


def interior_correction(trace, radius, mu, nangle, nstep, reflect):
    # Integrate -grad(delta j)/(2*pi²*|x-y|²) on rays from the observation
    # point. The 3D ray volume element cancels the kernel before quadrature.
    # Orient the polar axis along the observation radius, with mu its
    # direction cosine relative to the line between the mouths.
    t, wt = unit_gauss(nangle)
    cosine = np.r_[-t, t][:, None, None]
    sine = np.sqrt(1 - cosine**2)
    azimuth = (np.arange(nangle)+0.5)*np.pi/nangle
    az = np.cos(azimuth)[None, :, None]
    bz = np.sin(azimuth)[None, :, None]
    q, wq = unit_gauss(nstep)
    direction = np.sqrt(1 - mu**2)
    vx, va = np.empty_like(mu), np.empty_like(mu)
    root = np.sqrt(1-radius**2 + (radius*cosine)**2)
    length = root - radius*cosine
    step = q*length
    weights = -np.tile(wt, 2)[:, None, None] * length * wq / (nangle*np.pi)
    for i, (c, s) in enumerate(zip(mu, direction, strict=True)):
        nx = cosine*c - sine*az*s
        ny = cosine*s + sine*az*c
        xx = radius*c + step*nx
        yy = radius*s + step*ny
        aa = np.hypot(yy, step*sine*bz)
        gx, ga = trace.perturbation(reflect*xx, aa, gradient=True)
        vx[i] = np.sum(weights * reflect*gx)
        va[i] = np.sum(weights * ga * yy/aa)
    return vx, va


def companion_correction(trace, x, a, center, order, reflect):
    # Integrate delta j times the observation gradient of the Neumann kernel.
    # Azimuth is analytic: int dphi/(B-C*cos phi)^2 = 2*pi*B/(B²-C²)^(3/2),
    # and inserting cos phi in the numerator replaces B by C.
    r, mu, weights = ball(order)
    ys, ring = r*mu, r*np.sqrt(1 - mu**2)
    mass = weights*trace.perturbation(reflect*ys, ring)/np.pi**2
    vx = np.zeros_like(x)
    va = np.zeros_like(a)
    for p, q, m in zip(ys.flat, ring.flat, mass.flat, strict=True):
        u = (x-center) - p
        b, c = u**2 + a**2 + q**2, 2*a*q
        inv = m / ((b-c)*(b+c))**1.5
        vx += inv*u*b
        va += inv*(a*b-q*c)
    return vx, va


def measure(traces, d, order, angular, step, reflect=1):
    source, target = traces
    r, mu, weights = ball(order)
    x, a = r*mu, r*np.sqrt(1 - mu**2)
    center = reflect*d
    dist = np.hypot(x-center, a)
    vr = interior_rays(r[:, 0], 64)[:, None]
    vs = exterior_series(dist, 12)
    ux = vr*mu + vs*(x-center)/dist
    ua = vr*np.sqrt(1 - mu**2) + vs*a/dist
    radial_mass = PEAK*(1 - r**2)**3
    direct_j = target.perturbation(reflect*x, a)
    vtx, vta = np.empty_like(x), np.empty_like(a)
    for i, radius in enumerate(r[:, 0]):
        vtx[i], vta[i] = interior_correction(target, radius, mu, angular, step, reflect)
    vsx, vsa = companion_correction(source, x, a, center, max(order, 16), reflect)
    tx = weights*x
    base = np.sum(tx*(-vr*vs*(r-center*mu)/dist - vs**2/2))
    normal = np.sum(tx*(-radial_mass*direct_j - direct_j**2/2))
    target_tangent = np.sum(tx*(-ux*vtx - ua*vta - (vtx**2 + vta**2)/2))
    source_tangent = np.sum(tx*(-(ux+vtx)*vsx - (ua+vta)*vsa - (vsx**2 + vsa**2)/2))
    return dict(order=[order, angular, step], reflect=reflect,
                target_rate=np.sum(weights*(radial_mass+direct_j)),
                source_rate=np.sum(weights*(radial_mass+source.perturbation(reflect*x,a))),
                radial_moment=base, normal_correction=normal,
                target_tangent_correction=target_tangent, source_tangent_correction=source_tangent,
                corrected_bulk=base+normal+target_tangent+source_tangent)


def main():
    manifest = json.loads(FEED.read_text())
    hashes = {Path(k): v for k,v in manifest["raw_fields_sha256"].items()}
    lines = manifest["runs"]["mouth_front_pressure_dipole02"]["stdout_lines"]
    feeds = {(row["grid"][1], row["d"]): row["moment"]["pressure"]
             for line in lines if line.startswith("{") for row in [json.loads(line)]
             if row["grid"][0] == 96}
    for cells in (4, 6):
        for d in (8, 16):
            traces = load(cells, d, hashes)
            rows = []
            for args in ((12, 8, 8), (20, 14, 10), (30, 20, 16), (12, 8, 8, -1)):
                row = measure(traces, d, *args)
                print(json.dumps(dict(cells=cells, d=d, **row)), flush=True)
                rows.append(row)
            value = rows[2]["corrected_bulk"]
            feed = feeds[cells, d]
            print(json.dumps(dict(cells=cells, d=d, feed=feed, corrected_bulk=value,
                                  correction_to_radial=value-rows[2]["radial_moment"],
                                  bulk_over_feed=value/feed, sigma_s1=(feed-value)/(4*np.pi/3))),
                  flush=True)
            for row in rows:
                np.testing.assert_allclose([row["target_rate"], row["source_rate"]], RATE,
                                           atol=2e-10, rtol=0)
            np.testing.assert_allclose(rows[3]["corrected_bulk"], -rows[0]["corrected_bulk"],
                                       rtol=5e-11, atol=0)
            np.testing.assert_allclose(rows[1]["corrected_bulk"], value, rtol=5e-5, atol=0)


if __name__ == "__main__":
    main()
