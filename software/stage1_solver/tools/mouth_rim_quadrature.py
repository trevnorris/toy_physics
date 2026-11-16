"""Integrate saved sideways-load coefficients over cells clipped by the unit mouth."""

from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline

from mouth_feed_weak_check import dx4


def coefficients(packet):
    x, r = packet["x"], packet["r"]
    f, u = packet["field"], packet["u"]
    face, du = f[-1], u[-1]
    dx, dz = 2*r[0], float(packet["dz"])
    n, dn = abs(face)**2, 2*(f.conj()*u).real
    g = np.maximum(1-x[:, None]**2-r**2, 0)**3
    dv = 2*np.pi*r*dx**2
    current = (packet["normal"]/face).imag*n
    k = np.sum(current*dv)/np.sum(g*dv)
    if np.max(abs(current-k*g)) > 1e-13:
        raise ValueError("Saved radius-one intake differs from k times the taper")

    df, ddu = dx4(face, dx), dx4(du, dx)
    # Remove the known taper analytically. The smooth coefficients exist
    # even in a rim cell whose center lies outside the aperture.
    carried = k*(du.conj()*df+face.conj()*ddu).imag/n
    stress = (-dn[-1]*dx4(n, dx)/(4*abs(face)**3)
              -(1-1/abs(face))*dx4(dn[-1], dx)/2)

    def normal(q):
        return (11*q[-1]-18*q[-2]+9*q[-3]-2*q[-4])/(6*dz)

    mixed = (ddu.conj()*normal(f)+df.conj()*normal(u)).real-dx4(normal(dn), dx)/4
    return k, np.stack((carried, stress, mixed))


def clipped_cell(i, j, cells, nodes, weights):
    """Yield quadrature coordinates and 3D volume weights inside one annular cell."""
    x0, x1 = i/cells, (i+1)/cells
    r0, r1 = j/cells, (j+1)/cells
    rm = (j+0.5)/cells
    # Split at spline knots and wherever the sphere crosses a radial edge.
    for lo, hi in ((r0, rm), (rm, r1)):
        edges = [x0, (i+0.5)/cells, x1]
        for radius in (lo, hi):
            cross = np.sqrt(max(1-radius**2, 0))
            edges.extend(v for v in (-cross, cross) if x0 < v < x1)
        edges = sorted(set(edges))
        for a, b in zip(edges, edges[1:]):
            if ((a+b)/2)**2+lo**2 >= 1:
                continue
            x = (a+b+(b-a)*nodes)/2
            top = np.minimum(hi, np.sqrt(np.maximum(1-x*x, 0)))
            r = (top[:, None]+lo+(top[:, None]-lo)*nodes)/2
            vol = (np.pi/2)*(b-a)*weights[:, None]*(top[:, None]-lo)*weights*r
            yield np.broadcast_to(x[:, None], r.shape).ravel(), r.ravel(), vol.ravel()


def integrate(cells, x, r, coeff, k, order):
    cropped_x, cropped_r = abs(x) < 2, r < 2
    rr = r[cropped_r]
    splines = []
    for field in coeff:
        local = field[cropped_x][:, cropped_r]
        splines.append(RectBivariateSpline(x[cropped_x], np.r_[-rr[::-1], rr],
                                          np.column_stack((local[:, ::-1], local)), s=0))
    nodes, weights = leggauss(order)
    # Columns are carried, stress, brane-layer mixed flux and inward rate.
    modes = ("center", "rim_volume", "rim_constant", "rim_cubic", "all_cubic")
    totals = {name: np.zeros(4) for name in modes}
    moments = np.zeros(3)
    for i in range(-cells, cells):
        ix = len(x)//2+i
        for j in range(cells):
            if min(i*i, (i+1)**2)+j*j >= cells*cells:
                continue
            cut = max(i*i, (i+1)**2)+(j+1)**2 > cells*cells
            center = np.r_[coeff[:, ix, j], k]
            g0 = max(1-x[ix]**2-r[j]**2, 0)**3
            midpoint = center*np.array([g0, g0, float(g0 > 0), g0])*2*np.pi*r[j]/cells**2
            constant, cubic = np.zeros(4), np.zeros(4)
            volume = 0.
            for xp, rp, vol in clipped_cell(i, j, cells, nodes, weights):
                gv = vol*np.maximum(1-xp*xp-rp*rp, 0)**3
                v, w = vol.sum(), gv.sum()
                volume += v
                moments += (v, w, np.dot(gv, xp*xp))
                constant += center*np.array([w, w, v, w])
                cubic += [np.dot(gv, splines[0].ev(xp, rp)),
                          np.dot(gv, splines[1].ev(xp, rp)),
                          np.dot(vol, splines[2].ev(xp, rp)), k*w]
            totals["center"] += midpoint
            totals["rim_volume"] += center*np.array([g0, g0, 1., g0])*volume if cut else midpoint
            totals["rim_constant"] += constant if cut else midpoint
            totals["rim_cubic"] += cubic if cut else midpoint
            totals["all_cubic"] += cubic
    return totals, moments


def main():
    folder = Path("docs/_measurements")
    crossed = json.loads((folder/"mouth_drift_crossed_grid.json").read_text())
    diagonal = json.loads((folder/"mouth_drift_refinement.json").read_text())
    results = {}
    for run in diagonal["runs"]+crossed["runs"]:
        name, digest = next(iter(run["field_sha256"].items()))
        with open(name, "rb") as stream:
            if file_digest(stream, "sha256").hexdigest() != digest:
                raise ValueError(f"Field hash mismatch: {name}")
        with np.load(name, allow_pickle=False) as packet:
            x, r = packet["x"], packet["r"]
            cells, normal = round(1/(2*r[0])), round(1/float(packet["dz"]))
            if float(packet["radius"]) != 1 or cells != run["cells"]:
                raise ValueError("Unexpected mouth or lateral grid")
            k, coeff = coefficients(packet)
        values = dict(item.split("=") for line in run["stdout"] for item in line.split())
        if normal != int(values["normal_cells"]) or (cells, normal) in results:
            raise ValueError("Unexpected or repeated normal grid")
        low, _ = integrate(cells, x, r, coeff, k, 8)
        high, moments = integrate(cells, x, r, coeff, k, 16)
        step = max(np.max(abs(v-low[name])) for name, v in high.items())
        expected = np.array([4*np.pi/3, 64*np.pi/315, 64*np.pi/3465])
        geometry_error = np.max(abs(moments-expected))
        original = [float(values[k]) for k in ("carry","stress","layer_cut","rate")]
        replay_error = np.max(abs(high["center"]-original))
        print(f"cells={cells} normal={normal} original_error={replay_error:.12g} "
              f"geometry_error={geometry_error:.12g} order8_to_16_max={step:.12g}")
        result = {}
        for mode, (carry, stress, mixed, rate) in high.items():
            total = carry+stress
            result[mode] = total
            print(f"cells={cells} normal={normal} mode={mode} carry={carry:.12g} "
                  f"stress={stress:.12g} boundary={total:.12g} layer_ball={mixed:.12g} "
                  f"rate_readout={rate:.12g} boundary_at_matched_rate={carry*0.25/rate+stress:.12g}")
        if (not np.isfinite([*result.values(), step, geometry_error, replay_error]).all()
                or geometry_error > 1e-12 or replay_error > 1e-13 or step > 2e-12):
            raise ArithmeticError("Quadrature, geometry or midpoint replay check failed")
        results[cells, normal] = result
    if set(results) != {(8,32), (8,48), (12,32), (12,48)}:
        raise ValueError("Incomplete four-grid comparison")
    for normal in (32, 48):
        changes = {m: results[12, normal][m]-results[8, normal][m] for m in results[8, normal]}
        for mode, change in changes.items():
            print(f"normal={normal} mode={mode} lateral_change={change:.12g} "
                  f"change_minus_center={change-changes['center']:.12g} "
                  f"percent_removed={100*(1-change/changes['center']):.12g}")


if __name__ == "__main__":
    main()
