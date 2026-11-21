"""Read physical radial current at the circular rim from staggered complex-field links."""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


MANIFEST = Path('research/pde_ledger_v3/_measurements/mouth_front_pressure_dipole.json')


def link_fields(path, digest, distance):
    if sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f'Incorrect packet digest for {path}')
    with np.load(path, allow_pickle=False) as packet:
        if distance:
            np.testing.assert_array_equal(packet['centers'], [-distance, 0])
        else:
            assert int(packet['d']) == 0
        x, r = packet['x'], packet['r']
        dx = r[1] - r[0]
        np.testing.assert_allclose(np.diff(x), dx, rtol=1e-12)
        np.testing.assert_allclose(np.diff(r), dx, rtol=1e-12)
        ix, ir = abs(x) < 1 + 3*dx, r < 1 + 3*dx
        psi = packet['layers'][-1][ix][:, ir]
        x, r = x[ix], r[ir]
    # Mirror the scalar about the cylindrical axis. Its link flux is then
    # exactly zero on that axis, with its physical normal pointing toward +r.
    psi = np.column_stack((psi[:, 0], psi))
    r = np.r_[-r[0], r]
    jx = (psi[:-1].conj()*psi[1:]).imag/dx
    jt = (psi[:, :-1].conj()*psi[:, 1:]).imag/dx
    return (RegularGridInterpolator(((x[1:] + x[:-1])/2, r), jx),
            RegularGridInterpolator((x, (r[1:] + r[:-1])/2), jt))


def product_bounds(a, b, c, d):
    products = np.stack((a*c, a*d, b*c, b*d))
    return products.min(axis=0), products.max(axis=0)


def read_rim(jx, jt):
    # Split also at every interpolation line crossed by (cos theta, sin theta).
    # Both coordinates are monotonic on each arc. Floating-point roundoff at
    # mesh crossings is immaterial at the precision of the reported sign.
    breaks = [np.linspace(0, np.pi, 1001)]
    for interp in (jx, jt):
        x, r = interp.grid
        breaks.append(np.arccos(x[np.abs(x) < 1]))
        angles = np.arcsin(r[(r > 0) & (r < 1)])
        breaks.extend((angles, np.pi - angles))
    angles = np.unique(np.concatenate(breaks))
    x, r = np.cos(angles), np.sin(angles)
    x[[0, -1]], r[[0, -1]] = [1, -1], 0
    points = np.column_stack((x, r))
    sampled = jx(points)*x + jt(points)*r
    # Bound the entire arc, using the extrema of each bilinear component on
    # its enclosing rectangle, followed by interval multiplication with x,r.
    x0, x1 = x[1:], x[:-1]
    r0, r1 = np.minimum(r[:-1], r[1:]), np.maximum(r[:-1], r[1:])
    corners = np.stack([np.column_stack((u, v)) for u in (x0, x1) for v in (r0, r1)])
    bounds = []
    for interp, a, b in ((jx, x0, x1), (jt, r0, r1)):
        values = interp(corners)
        bounds.append(product_bounds(a, b, values.min(axis=0), values.max(axis=0)))
    lower, upper = np.sum(bounds, axis=0)
    return dict(far_pole=sampled[0], near_pole=sampled[-1],
                sample_min=sampled.min(), sample_max=sampled.max(),
                entire_rim_lower=lower.min(), entire_rim_upper=upper.max())


def main():
    manifest = json.loads(MANIFEST.read_text())['runs']['mouth_front_pressure_dipole02']
    assert manifest['guard_exit'] == manifest['child_exit'] == 0 and manifest['limits_verified']
    packets = {}
    for raw in manifest['stdout_lines']:
        row = json.loads(raw)
        if row['grid'] == [96, 6, 24]:
            for name, digest in row['inputs'].items():
                distance = 0 if 'baseline' in Path(name).stem else row['d']
                entry = (Path(name), digest)
                if distance in packets and packets[distance] != entry:
                    raise ValueError('Inconsistent feed packets')
                packets[distance] = entry
    assert set(packets) == {0, 8, 16}
    for distance, (path, digest) in sorted(packets.items()):
        result = read_rim(*link_fields(path, digest, distance))
        print(json.dumps(dict(d=distance, **result), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
