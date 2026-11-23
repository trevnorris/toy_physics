"""One-sided Q1 disk traces and independent elemental balances of the held R4 mouth."""
import argparse
import gc
import json
from pathlib import Path
import sys

import numpy as np

from mouth_whole_r4_carrier_refine import digest, endpoint, pin, read


def phase_terms(state):
    """Write j_z=A*p_z+B, with A and B dependent on the order gradient."""
    from mouth_disk_pillbox import dot, ORDER_INERTIA
    _, p, c, n, _, _ = state
    den = n+ORDER_INERTIA*dot(c, c)
    tangential_order = dot(c[:, :2], c[:, :2])
    projection = dot(p[:, :2], c[:, :2])
    A = n*(n+ORDER_INERTIA*tangential_order)/den
    B = -n*ORDER_INERTIA*c[:, 2]*projection/den
    return A, B


def disk(reader, degree):
    from mouth_disk_pillbox import gauss_cells, lag, stresses
    model = reader.m
    length = model.shape.length
    radius, slope = model.shape.radius(length)
    xs, wx = gauss_cells(model.rr[0], degree)
    ang, wa = gauss_cells(model.a, degree)
    x, a = [elem.ravel() for elem in np.meshgrid(xs, ang, indexing='ij')]
    w = 4*np.pi*(radius*x)**2*np.sin(a)*radius*(wx[:, None]*wa).ravel()
    banks = [reader.sample(s, x, a, length) for s in (1, 2)]
    lower, upper = [bank[0] for bank in banks]
    pl, pu = lower[1], upper[1]
    Al, Bl = phase_terms(lower)
    Au, Bu = phase_terms(upper)
    A = .5*(Al+Au)
    # Positive orientation throughout is toward the receiving bulk. The lower chart
    # follows the sloped wall, whereas r=b*x is depth-independent on the upper side.
    map_correction = slope*x*pl[:, 0]
    contributions = dict(phase_chart_normal=A*(pu[:, 2]-pl[:, 2]-map_correction),
                         phase_chart_to_physical=A*map_correction,
                         order_normal_coefficient=.5*(pu[:, 2]+pl[:, 2])*(Au-Al),
                         tangential_projection=Bu-Bl)
    output = dict(degree=degree, disk_area=float(w.sum()), banks={},
                  max_pointwise_jump=dict(density=float(abs(upper[3]-lower[3]).max()),
                      order=float(abs(upper[5]-lower[5]).max()),
                      tangential_phase_gradient=float(abs(pu[:, :2]-pl[:, :2]).max()),
                      tangential_order_gradient=float(abs(upper[2][:, :2]-lower[2][:, :2]).max())))
    currents = []
    for name, (u, rest) in zip(('incident', 'receiving'), banks):
        j = u[3]*lag(u)[1][:, 2]
        currents.append(j)
        T = stresses(u)[1].sum(axis=0)
        T0 = stresses(rest)[1].sum(axis=0)
        output['banks'][name] = dict(signed_current=float(w@j),
            rest_current=float(w@(rest[3]*lag(rest)[1][:, 2])), positive_current=float(w@np.maximum(j,0)),
            axial_phase_term_without_order=float(w@(u[3]*u[1][:, 2])),
            order_momentum_correction=float(w@(j-u[3]*u[1][:, 2])),
            axial_momentum=float(w@T), drive_minus_rest_axial_momentum=float(w@(T-T0)))
    diff = currents[1]-currents[0]
    if np.max(abs(diff-sum(contributions.values()))) > 2.e-11:
        raise ArithmeticError('Common-face phase/order decomposition failed')
    output['receiving_minus_incident'] = float(w@diff)
    output['contributions'] = {key: float(w@term) for key, term in contributions.items()}
    breaks = model.rr[0][[0, *[np.argmin(abs(model.rr[0]-v)) for v in (.9, .98)], -1]]
    output['radial_partition'] = []
    for lo, hi in zip(breaks, breaks[1:]):
        keep = (lo<=x) & (x<hi)
        output['radial_partition'].append(dict(fraction=[float(lo), float(hi)],
            receiving_minus_incident=float(w[keep]@diff[keep]),
            absolute_pointwise_difference_integral=float(w[keep]@abs(diff[keep]))))
    return output


def rim_sides(reader, degree):
    """Raw outward j and T_zi on r=b(w) inside each immediately adjoining cell."""
    from mouth_disk_pillbox import gauss_cells, lag, stresses
    m = reader.m
    ang, wa = gauss_cells(m.a, degree)
    rows = {}
    for slab, tag in ((1, 'sealed_duct_wall'), (2, 'receiving_cylinder')):
        knots = m.z[1][-2:] if slab==1 else m.z[2][:2]
        zs, wz = gauss_cells(knots, degree)
        z, a = [e.ravel() for e in np.meshgrid(zs, ang, indexing='ij')]
        b, slope = m.shape.radius(z)
        weight = 4*np.pi*b*b*np.sin(a)*(wz[:, None]*wa).ravel()
        u, rest = reader.sample(slab, np.ones_like(z), a, z)
        vel = lag(u)[1]
        shear, axial = (v.sum(axis=0) for v in stresses(u))
        sr, ar = (v.sum(axis=0) for v in stresses(rest))
        rows[tag] = dict(axial_width=float(knots[1]-knots[0]),
            signed_outward_current=float(weight@(u[3]*(vel[:, 0]-slope*vel[:, 2]))),
            drive_minus_rest_outward_axial_momentum=float(weight@(shear-sr-slope*(axial-ar))))
    return rows


def weak(model, wound, rest):
    """Integrate test gradients on disjoint core and receiving elements sharing the trace."""
    from mouth_disk_pillbox import lag, physical, stresses
    # Unit nodal disk lift. The exterior sectors touch its rim hat on exactly the
    # same footing as the central receiver cells. The open-interior hat omits the
    # rim nodes and tapers to zero in the last radial cell on both sides.
    plane = model.ids[1, 0][:, :, -1]
    if not np.array_equal(plane, model.ids[2, 0][:, :, 0]):
        raise ArithmeticError('Unshared cap degrees of freedom')
    lifts = np.zeros((model.size, 2))
    lifts[plane, 0] = 1.
    lifts[plane[:-1], 1] = 1.
    grads = np.stack([G@lifts for G in model.G], axis=2)
    blocks = {}
    start = 0
    for (slab, bank), ids in model.ids.items():
        stop = start+int(8*np.prod(np.array(ids.shape)-1))
        rows = np.flatnonzero(np.any(grads[start:stop]!=0, axis=(1,2)))+start
        start = stop
        if not len(rows):
            continue
        D = model.D[rows]
        wt, r, a, z, F, Fp = D[:, :6].T
        if np.any(z <= .62*model.shape.length):
            raise ArithmeticError('Lift enters the polar deformation support')
        vals = model.V[rows]
        gg = [G[rows] for G in model.G]
        current, force = [], []
        for u in (wound, rest):
            values = (vals@u.T).T
            values[2] += model.ref[rows, 0]
            derivatives = np.stack([(G@u.T).T for G in gg], axis=2)
            derivatives[2] += model.ref[rows, 1:]
            state = physical(model.winding, r, a, values, derivatives, F, Fp)
            vel = lag(state)[1]
            j = state[3][:, None]*vel
            gf, _, c, n, _, _ = state
            # T_z,theta has a polar component even under axisymmetry. The common
            # lift is angularly constant; accumulate the full tensor row anyway.
            from mouth_breathing_mobile_patch import KAPPA
            T = gf*gf[:, 2:3]+(n*vel[:, 2])[:, None]*vel+KAPPA*c*c[:, 2:3]
            T[:, 2] = stresses(state)[1].sum(axis=0)
            current.append(np.einsum('q,qti,qi->t',wt,grads[rows],j))
            force.append(np.einsum('q,qti,qi->t',wt,grads[rows],T))
        blocks[slab, bank] = np.stack((current[0], current[1], force[0]-force[1]))
    if start != len(model.D) or set(blocks) != {(1,0),(2,0),(1,1),(2,1)}:
        raise ArithmeticError('Unrecognized test support or quadrature layout')
    rows = {}
    for it, tag in enumerate(('whole_disk', 'open_interior_lift')):
        rows[tag] = {}
        for field, col in (('current',0),('rest_current',1),('drive_minus_rest_fluid_axial_work',2)):
            terms = {f'{slab}_{bank}': float(value[col, it]) for (slab, bank), value in blocks.items()}
            incident = terms['1_0']
            receiving = -(terms['2_0']+terms['1_1']+terms['2_1'])
            rows[tag][field] = dict(elemental_outward_work=terms, incident=incident,
                                    receiving=receiving, receiving_minus_incident=receiving-incident)
    return rows


def measure(base, model, wound, rest, port):
    from mouth_disk_pillbox import FieldTraces
    reader = FieldTraces(model, base.c, [wound, rest])
    traces = [disk(reader, d) for d in (2,6,12,18)]
    local = weak(model, wound, rest)
    flux = local['whole_disk']['current']
    out = dict(depth=float(model.shape.length), radius=float(model.shape.exit),
               incident_radius_slope=float(model.shape.radius(model.shape.length)[1]),
               carrier_sizes=[len(model.rr[0])-1,len(model.a)-1,len(model.z[1])-1,len(model.z[2])-1],
               port_current=port, raw=traces, weak=local, raw_rim_sides=rim_sides(reader,12),
               finest_relative_gap=traces[-1]['receiving_minus_incident']/port,
               checks=dict(weak_ports=max(abs(flux[k]-port) for k in ('incident','receiving'))<3.e-5,
                   common_traces=max(traces[-1]['max_pointwise_jump'].values())<1e-9,
                   face_quadrature=max(abs(traces[-1]['banks'][k]['signed_current']-
                       traces[-2]['banks'][k]['signed_current']) for k in ('incident','receiving'))<1e-6))
    return out


def run(args):
    refs = pin(args.run, args.edges)
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_mobile_joined_fields import JointFields
    from mouth_rim_pinned_comparison import validate_guard
    report = dict(command=sys.argv, script_sha256=digest(Path(__file__)),
                  input_reader_sha256=digest(Path(__file__).with_name('mouth_whole_r4_carrier_refine.py')),
                  source_manifest_sha256=digest(args.run/'sources.json'), edge_readback_sha256=digest(args.edges), results={})
    for mesh in (('fine',) if args.refined else ('coarse','fine')):
        ref = refs[mesh]
        state, x, _, wound, rest, original = endpoint(args.run, mesh, ref)
        base = Front(state, GRIDS['double'], 1)
        port = original['controls']['wound']['feed_current']
        if args.refined:
            folder, log = Path(args.refined[0]), Path(args.refined[1])
            validate_guard(folder.parent)
            r = read(log)
            receipt = json.loads((folder.parent/'stdout').read_text().splitlines()[-1])
            if (receipt['output_sha256']!=digest(log) or r['original_sha256']!=ref['state_sha256']
                    or r['state_sha256']!=ref['preparation_sha256']):
                raise ValueError('Unmatched field refinement')
            pair = []
            for name in ('wound','rest'):
                path = folder/(name+'.npz')
                if digest(path)!=r['comparisons'][name]['sha256']:
                    raise ValueError('Changed rematched field')
                with np.load(path) as f:
                    pair.append(f['field'])
            wound, rest = pair
            port = r['comparisons']['wound']['rematched']['feed_current']
            report['refinement_sha256'] = digest(log)
        f = JointFields(base, state, x, grid='local' if args.refined else 'double', precondition=False)
        for field in (wound, rest):
            f._check_field(f.model, field)
        result = measure(base, f.model, wound, rest, port)
        result['original_endpoint_sha256'] = ref['endpoint_sha256']
        if not args.refined:
            archive = {row['side']:row['signed_current'] for row in original['audit']['physical_disk'] if row['degree']==12}
            recovered = result['raw'][-2]['banks']
            result['checks']['original_face_reader'] = max(abs(recovered[key]['signed_current']-archive[s]) for key,s
                in (('incident',-1),('receiving',1))) < 1e-10
        report['results'][mesh] = result
        del base, f
        gc.collect()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(sha256=digest(args.output), checks={k:v['checks'] for k,v in report['results'].items()})))
    if not all(all(item['checks'].values()) for item in report['results'].values()):
        raise ArithmeticError('A disk diagnostic check failed; see output')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run', 'edges', 'output'):
        p.add_argument(key, type=Path)
    p.add_argument('--refined', nargs=2, metavar=('FIELD_FOLDER','READBACK'))
    run(p.parse_args())
