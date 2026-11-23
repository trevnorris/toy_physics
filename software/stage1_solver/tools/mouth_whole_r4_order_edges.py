"""Measure the correction constraints on the two saved R4 whole-mouth endpoints.

Run as a guarded readback.  Supply the original run directory, output JSON
and output figure.  The field/solid definitions are imported from its pin.
"""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def read(path):
    return json.loads(path.read_text())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def annulus(model, base_order, driven, archived):
    """Radial cell containing the median dense, subsonic, positive incident flux."""
    from mouth_joined_front_trace import face
    w, n, mach, j, _ = face(model, base_order, driven, -1, 12)
    pos = w*np.maximum(j, 0.)*((n > .5) & (mach < 1.))
    cells = pos.reshape(12, 12, len(model.rr[0])-1, len(model.a)-1).sum(axis=(0, 1, 3))
    if abs(cells.sum()-archived['dense_subsonic_current']) > 1.e-10 or cells.sum() <= 0:
        raise ArithmeticError('Dense incident annulus differs from the endpoint audit')
    i = np.searchsorted(np.cumsum(cells), cells.sum()/2)
    return model.rr[0][i:i+2].tolist()


def charts(model, base_order, driven, rest):
    base = RegularGridInterpolator(tuple(np.linspace(0., 1., n) for n in base_order.shape), base_order)
    result = {}
    for slab in (1, 2):
        x, _, z = np.meshgrid(model.rr[0], model.a, model.z[slab], indexing='ij')
        ref = (base((z/model.shape.length, x)) if slab == 1 else
               model.terminal_profile(x, z-model.shape.length))
        ids = model.ids[slab, 0]
        fields = {}
        for name, u in (('wound', driven), ('rest', rest)):
            values = np.moveaxis(u[:, ids], 0, -1).copy()
            values[..., 2] += ref
            fields[name] = values
        result[slab] = ref, fields
    return result


def profile(model, chart, x):
    # A meridian at fixed r/b(z) and theta=pi/2. On this curve the Q1 field
    # is piecewise linear in depth, including the native cells outside the band.
    rows = []
    for slab, (ref, fields) in chart.items():
        rx = model.rr[0]
        i = min(np.searchsorted(rx, x, side='right')-1, len(rx)-2)
        f = (x-rx[i])/(rx[i+1]-rx[i])
        interp = lambda v: (1-f)*v[i, -1]+f*v[i+1, -1]
        z = model.z[slab]-model.shape.length
        row = np.array([z, interp(ref), interp(fields['rest'][..., 2]),
                        interp(fields['wound'][..., 2])]).T
        if slab == 2 and not np.allclose(rows[-1][-1], row[0], rtol=0, atol=3.e-12):
            raise ArithmeticError('Order discontinuous at the common disk trace')
        rows.append(row if slab == 1 else row[1:])
    return np.vstack(rows)


def half_sites(profile, h, col):
    """Represent isolated half-order crossings and exactly half-order intervals."""
    band = profile[np.abs(profile[:, 0]) <= h+1.e-13]
    z, y = band[:, 0], band[:, col]-.5
    y[np.abs(y) <= 5.e-14] = 0.
    pieces = []
    for j in range(len(z)-1):
        if y[j] == y[j+1] == 0:
            pieces.append([z[j], z[j+1]])
        elif y[j]*y[j+1] < 0:
            cross = z[j]-(z[j+1]-z[j])*y[j]/(y[j+1]-y[j])
            pieces.append([cross, cross])
        elif y[j] == 0 or y[j+1] == 0:
            site = z[j] if y[j] == 0 else z[j+1]
            pieces.append([site, site])
    merged = []
    for lo, hi in pieces:
        if merged and abs(lo-merged[-1][1]) < 1.e-13:
            merged[-1][1] = hi
        else:
            merged.append([lo, hi])
    return [dict(offset_range=[a, b], distance_to_lower_edge_range=[h+a, h+b],
                 distance_to_upper_edge_range=[h-b, h-a]) for a, b in merged]


def edge_trace(model, chart, slab, x, j, edge):
    """Physical one-sided trace from the specified cell at its axial endpoint."""
    from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA
    from mouth_wound_carrier import straight_vortex
    rx, angles, zs = model.rr[0], model.a, model.z[slab]
    i = min(np.searchsorted(rx, x, side='right')-1, len(rx)-2)
    f = (x-rx[i])/(rx[i+1]-rx[i])
    z = zs[edge]
    radius, slope = model.shape.radius(z)
    r = radius*x
    F, _ = straight_vortex(model.winding, 28.)(r)
    ref, fields = chart[slab]
    output = {}
    for name, values in fields.items():
        def sample(plane):
            return (1-f)*plane[i]+f*plane[i+1]
        line = sample(values[:, -1])
        at = line[edge]
        chart_slope = (line[j+1]-line[j])/(zs[j+1]-zs[j])
        radial = (values[i+1, -1, edge]-values[i, -1, edge])/(radius*(rx[i+1]-rx[i]))
        axial = chart_slope-x*slope*radial
        if x:
            angular = sample(values[:, -1, edge]-values[:, -2, edge])/(r*(angles[-1]-angles[-2]))
            _, phase, order = np.stack((radial, angular, axial), axis=1)
            n = (F*at[0])**2
            c = ORDER_INERTIA*(phase@order)/(n+ORDER_INERTIA*(order@order))
            flux = KAPPA*order-n*c*(phase-c*order)
            flux_vector = flux.tolist()
        else:
            # The vortex factor vanishes at the origin and hence pi_chi=0 in
            # the regular one-sided limit. Only the normal component is used.
            flux_vector = [None, None, float(KAPPA*axial[2])]
        output[name] = dict(chi=float(at[2]), slope_along_chart=float(chart_slope[2]),
                            physical_axial_slope=float(axial[2]), order_flux_r_theta_z=flux_vector)
    baseline = (1-f)*ref[i, -1]+f*ref[i+1, -1]
    bz = (baseline[j+1]-baseline[j])/(zs[j+1]-zs[j])
    br = (ref[i+1, -1, edge]-ref[i, -1, edge])/(radius*(rx[i+1]-rx[i]))
    return dict(native_chi=float(baseline[edge]), native_physical_axial_slope=float(bz-x*slope*br),
                axial_cell_width=float(zs[j+1]-zs[j]), **output)


def weak_reactions(model, driven, rest, load):
    """dE/d(delta chi) at fixed chi0, shape, solid and all other fields."""
    from mouth_breathing_mobile_patch import KAPPA, WELL
    model.material_load = load
    c = model.ref[:, 0]
    ref = model.V.T@(model.hq*2*WELL*c*(1-c)*(1-2*c))
    for i, gradient in enumerate(model.G):
        ref += gradient.T@(model.hq*KAPPA*model.ref[:, i+1])
    g = model.evaluate(driven)[1][2]
    r = model.evaluate(rest)[1][2]
    terms = dict(reference_well_and_gradient=ref, crystal_deposition=load,
                 rest_field_minus_reference=r-ref-load, drive_minus_rest=g-r,
                 wound_total=g, rest_total=r)
    weights = {}
    for slab, z in ((1, model.shape.length-model.h), (2, model.shape.length+model.h)):
        j = int(np.argmin(abs(model.z[slab]-z)))
        plane = model.ids[slab, 0][:, :, j]
        regions = dict(open_disk=plane[:-1], lateral_rim=plane[-1:])
        if np.max(abs(model.z[slab][j]-z)) > 2.e-14 or model.masks[2, plane].any():
            raise ArithmeticError('Plane is not a constrained order edge')
        weights[str(slab)] = {}
        for region, ids in regions.items():
            ids = np.unique(ids)
            weights[str(slab)][region] = {name: dict(signed_sum=float(t[ids].sum()),
                                                      nodal_abs_sum=float(abs(t[ids]).sum()),
                                                      peak_abs_node=float(abs(t[ids]).max()))
                                          for name, t in terms.items()}
    return weights


def measure(run, mesh):
    from mouth_free_shape_solid import Solid
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_mobile_joined_fields import JointFields
    from mouth_rim_pinned_comparison import validate_guard
    folder = run/mesh
    validate_guard(folder)
    command = read(folder/'invocation.json')['command']
    state = Path(command[command.index('-u')+2])
    report = read(folder/'endpoint.json')
    checkpoint = folder/report['file']
    endpoint_hash = digest(folder/'endpoint.json')
    receipt = json.loads((folder/'stdout').read_text().splitlines()[-1])
    if receipt['kind'] != 'mesh_endpoint' or receipt['endpoint'] != endpoint_hash or digest(checkpoint) != report['sha256']:
        raise ValueError('Unmatched guarded endpoint')
    with np.load(checkpoint) as f, np.load(folder/'endpoint_controls.npz') as controls:
        if str(f['state_sha']) != digest(state):
            raise ValueError('Parent order input changed')
        x, q, u, rest = f['x'], f['q'], f['field'], controls['rest']
    base = Front(state, GRIDS['double'], 1)
    fields = JointFields(base, state, x, precondition=False)
    m = fields.model
    if m.shape.length-m.h <= .62*m.shape.length:
        raise ArithmeticError('Polar support intersects the band')
    for name, v in (('wound', u), ('rest', rest)):
        fields._check_field(m, v)
    options = {key: int(command[command.index('--'+key)+1]) for key in ('n', 'k', 'angles', 'degree', 'subdivisions')}
    b = Solid(state, m.shape, m.polar_shift, options['n'], options['k'], options['angles'], m,
              degree=options['degree'], subdivisions=options['subdivisions'])
    _, _, load, _ = b.evaluate(q, gradient=False)
    weak = weak_reactions(m, u, rest, load)
    archived = next(row for row in report['audit']['physical_disk'] if row['side'] == -1 and row['degree'] == 12)
    shell = annulus(m, base.c, u, archived)
    chart = charts(m, base.c, u, rest)
    sites = {}
    for key, frac in (('axis', 0.), ('annulus', np.mean(shell)), ('near_rim', .98), ('bonded_rim', 1.)):
        row = profile(m, chart, frac)
        edge_values = {}
        for slab, side, sign in ((1, 'lower', -1), (2, 'upper', 1)):
            j = int(np.argmin(abs(m.z[slab]-m.shape.length-sign*m.h)))
            inside, outside = (j, j-1) if slab == 1 else (j-1, j)
            edge_values[side] = {bank: edge_trace(m, chart, slab, frac, k, j)
                                 for bank, k in (('inside', inside), ('outside', outside))}
        # Include the neighboring frozen cells beyond each edge. These spans
        # are well inside the polar-support-free neighborhood of the cap.
        keep = (row[:, 0] >= -.25) & (row[:, 0] <= .25)
        sites[key] = dict(radial_fraction=float(frac), disk_radius=float(frac*m.shape.exit),
                          columns=['w_minus_L','reference','rest','wound'], profile=row[keep].tolist(),
                          half_order={k: half_sites(row, m.h, i) for i, k in
                                      enumerate(('reference', 'rest', 'wound'), 1)}, edges=edge_values)
    return dict(endpoint_sha256=endpoint_hash, state_sha256=report['sha256'],
                controls_sha256=digest(folder/'endpoint_controls.npz'), preparation_sha256=digest(state),
                depth=float(m.shape.length), exit_radius=float(m.shape.exit), halfwidth=m.h,
                carrier_grid=list(GRIDS['double']), solid_discretization=options,
                annulus_fraction_cell=shell, sites=sites, weak_reactions=weak)


def draw(results, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.4), sharex=True, sharey=True, constrained_layout=True)
    for axrow, (mesh, record) in zip(axes, results.items()):
        for ax, (name, site) in zip(axrow, record['sites'].items()):
            data = np.asarray(site['profile'])
            ax.axvline(-record['halfwidth'], color='#687783', lw=.9, ls=':')
            ax.axvline(record['halfwidth'], color='#687783', lw=.9, ls=':')
            ax.axhline(.5, color='#b0bac2', lw=.6)
            for col, color, style, label in ((3, '#126482', '-', 'Driven'), (2, '#aa4222', '--', 'Rest'),
                                               (1, '#333333', ':', 'Native reference')):
                ax.plot(data[:, 0], data[:, col], style, color=color, label=label, lw=1.7)
            if mesh == 'coarse':
                ax.set_title(name.replace('_', ' ').capitalize())
            ax.set(xlim=(-.19, .19), ylim=(-.03, 1.03))
        axrow[0].set_ylabel(f'{mesh.capitalize()} solid\nTotal order χ')
    for ax in axes[-1]:
        ax.set_xlabel('Depth w − L')
    axes[1, 0].legend(loc='upper right', frameon=False, fontsize=8)
    fig.savefig(path, dpi=155)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('figure', type=Path)
    args = parser.parse_args()
    pin = args.run/'pin'
    for path, sha in read(args.run/'sources.json').items():
        if digest(pin/path) != sha:
            raise ValueError(f'Pinned input changed: {path}')
    sys.path.insert(0, str(pin/'software/stage1_solver/tools'))
    results = {}
    for mesh in ('coarse', 'fine'):
        results[mesh] = measure(args.run, mesh)
        gc.collect()
    draw(results, args.figure)
    args.output.write_text(json.dumps(dict(command=sys.argv, pin_manifest_sha256=digest(args.run/'sources.json'),
                                          results=results), indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
