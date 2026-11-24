"""Read the two saved R4 label maps at the projected held stop, without a relaxation step.

Separate whole-chart orientation, melt-gauge corner determinants and physical
inverse volumes on the production integration grid. The last failed energy
proposal is absent from both endpoint files.
"""
import argparse
import ctypes
import hashlib
import json
import subprocess
import sys
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import mpmath as mp
import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinate(point):
    return np.array([point[k] for k in ('cylindrical_radius', 'axial_brane_coordinate', 'w')])


def split_kernel(pin, work, solid):
    """Change just the first extrema accumulator from global to strictly melted cells."""
    from mouth_r4_eulerian_solid import _kernel

    src = (pin/'mouth_r4_eulerian_solid.cpp').read_text()
    old = 'if(state.Jinv<diagnostics[0])'
    if src.count(old) != 1:
        raise ValueError('Cannot identify native minimum accumulator')
    src = src.replace(old, 'if(!has_order && state.Jinv<diagnostics[0])')
    path = work/'melt_only_min.cpp'
    path.write_text(src)
    library = work/'melt_only_min.so'
    subprocess.run(['g++', '-std=c++17', '-O3', '-shared', '-fPIC', '-I', str(pin),
                    str(path), '-o', str(library)], check=True)
    function = ctypes.CDLL(str(library)).r4_eulerian_material
    function.argtypes = _kernel().argtypes
    function.restype = ctypes.c_int
    solid.function = function
    return sha(path)


def corner_sites(axes, labels, active):
    """Locate one-sided determinant minima using differences of the four defining nodes."""
    rows = {}
    dims = active.shape
    widths = np.broadcast_arrays(*np.ix_(*map(np.diff, axes)))
    for corner in product((0, 1), repeat=3):
        columns = []
        for d in range(3):
            ends = []
            for v in (0, 1):
                bits = list(corner)
                bits[d] = v
                ends.append(labels[tuple(slice(b, b+n) for b,n in zip(bits, dims))])
            columns.append((ends[1] - ends[0])/widths[d][..., None])
        values = np.einsum('...i,...i->...', columns[0], np.cross(columns[1], columns[2]))
        for key, mask in [('active', active), ('melt', ~active)]:
            cell = np.unravel_index(np.argmin(np.where(mask, values, np.inf)), dims)
            if key not in rows or values[cell] < rows[key][0]:
                rows[key] = (float(values[cell]), tuple(map(int, cell)), corner)
    return rows


def local_corner(body, axes, labels, order, cell, bits):
    """High-precision one-sided Q1 gradient, including the angular hoop limit on a pole."""
    from mouth_r4_eulerian_corner_audit import num, profile
    from mouth_r4_refinement_review import _collar_reference
    node = tuple(i+b for i,b in zip(cell,bits))
    x, t, a = [num(ax[i]) for ax,i in zip(axes, node)]
    xx, tt, aa = map(num, labels[node])
    length = num(body.data.length)
    r, rt = profile(body,t)
    y = (t-num(.06))/num(.56)
    band = dband = mp.mpf(0)
    if 0 < y < 1:
        band = 16*y*y*(1-y)**2*num(body.epsilon)
        dband = 32*y*(1-y)*(1-2*y)*num(body.epsilon)/num(.56)
    leg = (3*mp.cos(a)**2-1)/2
    radius = x*r + x**3*band*leg
    physical = dict(radius=float(radius), theta=float(a), cylindrical_radius=float(radius*mp.sin(a)),
                    axial_brane_coordinate=float(radius*mp.cos(a)), w=float(t*length))
    Q = mp.matrix(3,3)
    for d,axis in enumerate(axes):
        low = list(node)
        high = list(node)
        low[d] = cell[d]
        high[d] = cell[d]+1
        h = num(axis[high[d]])-num(axis[low[d]])
        for k in range(3):
            Q[k,d] = (num(labels[tuple(high)][k])-num(labels[tuple(low)][k]))/h
    if not radius or not xx:
        raise ValueError('A radial-axis witness needs a separate radial limit')
    B = mp.matrix([[r+3*x*x*band*leg, x*rt+x**3*dband*leg, -3*x**3*band*mp.sin(a)*mp.cos(a)],
                   [0,length,0], [0,0,radius]])
    if xx <= num(body.cut):
        rr, rrt = profile(body, tt, True)
        ratio = num(body.data.core_ratio)
        C = mp.matrix([[rr*ratio,xx*rrt*ratio,0], [0,num(body.data.ref_length),0], [0,0,xx*rr*ratio]])
        side = 'core'
    else:
        C = _collar_reference(body, xx, tt)
        side = 'collar'
    if a == 0 and aa != 0:
        raise ValueError('Unregistered polar label at a symmetry boundary')
    hoop_inv = (C[2,2]/radius)*(mp.sin(aa)/mp.sin(a) if a else Q[2,2])
    grad = Q*(B**-1)
    j_inv = mp.det(Q)*abs(mp.det(C))*hoop_inv/mp.det(B)
    if j_inv <= 0:
        raise ValueError('Saved corner inverts')
    F = (C*grad)**-1
    J = 1/j_inv
    W = (mp.norm(F)**2+1/hoop_inv**2-4)/2 - mp.log(J) + (J-1)**2/4
    chi = float(order[node])
    cell_order = order[tuple(slice(i, i+2) for i in cell)]
    return dict(cell=list(cell), local=list(bits), node=list(node), physical=physical,
                chart_point=list(map(float,(x,t,a))), label=list(map(float,(xx,tt,aa))), material_side=side,
                chi=chi, positive_order_vertices=int(np.count_nonzero(cell_order)),
                max_vertex_order=float(cell_order.max()), chart_det=float(mp.det(Q)),
                physical_Jinv=float(j_inv), formal_W=float(W), chi_W=float(chi*W),
                chi_Jinv_W=float(chi*j_inv*W), grad_xi_norm=float(mp.norm(grad)))


def gradient_witness(barrier, gradient, order, axes, body):
    from mouth_cp10_label_condition import physical_frame
    rows = []
    for c, cm in enumerate(barrier.components):
        g = cm.sum(gradient[...,c])
        groups = np.flatnonzero(cm.dead)
        group = groups[np.argmax(np.abs(g[groups]))]
        nodes = np.argwhere(cm.groups.reshape(order.shape) == group)
        node = tuple(nodes[0])
        chart = np.array([[ax[n] for ax,n in zip(axes,node)]])
        frame, hoop = physical_frame(body, chart)
        radius = frame[0,2,2]
        rows.append(dict(component=c,group=int(group),node=[int(n) for n in node],equivalent_nodes=nodes.tolist(),
                         signed_unit_gradient=float(g[group]),chi=float(order[node]),
                         physical=dict(cylindrical_radius=float(hoop[0]),
                                       axial_brane_coordinate=float(radius*np.cos(chart[0,2])),
                                       w=float(chart[0,1]*body.data.length))))
    return max(rows,key=lambda r:abs(r['signed_unit_gradient']))


def run(area, previous, output):
    pin = area/'pin/software/stage1_solver/tools'
    sys.path.insert(0,str(pin))
    from mouth_cp10_label_grid import get_cp10
    from mouth_r4_eulerian_corner_audit import point
    from mouth_r4_eulerian_relax_labels import LabelRelaxation, chart_orientation
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_melt_barrier import MeltBarrier
    from mouth_r4_refinement_review import _collar_reference
    from mouth_r4_stage2_outcome import verify_guard
    verify_guard(area/'primary_guard', 1)
    folder = area/'primary'
    body, cp, axes = get_cp10(area.parent)
    docs = {name:read(folder/(name+'.json')) for name in ('held_stop','failed_trial','last_active')}
    states = {}
    for name,doc in docs.items():
        path = folder/doc['file']
        if sha(path) != doc['sha256']:
            raise ValueError(f'State digest differs for {name}')
        with np.load(path,allow_pickle=False) as f:
            if (int(f['radius']), int(f['winding']), int(f['iteration']), float(f['epsilon'])) != (4,1,10,1e-8):
                raise ValueError('Incorrect pinned preparation')
            if not np.array_equal(f['x'], cp['x']):
                raise ValueError('Held shape differs')
            states[name] = (f['labels'], f['field'])
    accepted, fields = states['held_stop']
    inner, trial_fields = states['failed_trial']
    if not np.array_equal(accepted,states['last_active'][0]) or not np.array_equal(fields, trial_fields):
        raise ValueError('Endpoint map/field pairing differs')
    order = np.ascontiguousarray(fields[2,body.ids].transpose(0,2,1))
    opt = LabelRelaxation(body,accepted,order,1e-8)
    old = read(previous)['accepted_corner']
    old_xyz = coordinate(old['physical'])
    old_cell = tuple(old['cell'])
    old_corner = tuple(map(int,old['local']))
    journal_path = folder/'inner_0000034.jsonl'
    journal = [json.loads(line) for line in journal_path.read_bytes().splitlines()]
    status = docs['failed_trial']['inner']
    if (journal[-1]['i'] != status['accepted_steps'] or status != docs['held_stop']['stop']['inner']
            or journal[-1]['melt_residual'] != status['residual']):
        raise ValueError('Endpoint inner journal and stop differ')
    result = dict(command=' '.join(sys.orig_argv),
        state_sha256={name:d['sha256'] for name,d in docs.items()},
        input_sha256={str(p):sha(p) for p in (previous,journal_path,area/'sources.json',
                       *(folder/(n+'.json') for n in docs), area.parent/'r4_melted_survival01/coarse/state_00010.npz',
                       area.parent/'supercore_joined_front01/adopt309/state.npz',area.parent/'r4_whole_mouth01/coarse/endpoint.json')},
        physical_quadrature='Three-point Gauss rules in depth, angle and each native cut radial segment. '
            'Ranges of signed physical Jinv over visited points. Exactly melted cells have eight zero-order nodes.',
        coordinates='Cells, corners and nodes are zero-based (radial fraction, depth fraction, polar angle). '
            'Physical distances are Euclidean in (cylindrical radius, axial brane coordinate, depth) at common azimuth.',
        old_corner={k:old[k] for k in ('cell','local','physical','chart_det','Jinv','chi')},
        active_cells=int(opt.active.sum()), strictly_melted_cells=int((~opt.active).sum()),
        max_active_trace_difference=float(abs((accepted-inner)[opt.stencil]).max()),
        max_label_difference=float(abs(accepted-inner).max()), labels_changed_components=int(np.count_nonzero(accepted!=inner)),
        search=dict(stop=status, last_accepted_inner_step=journal[-1],
            residual_target_unscaled=3e-8, next_fraction_after_last_rejection=status['alpha'],
            last_rejected_fraction=2*status['alpha'], final_rejected_perturbation_saved=False),
        states={})
    with TemporaryDirectory(prefix='r4_projected_orientation_') as tmp:
        work = Path(tmp)
        barrier = MeltBarrier(body,opt,work)
        solid = SpatialCrystal(body)
        result['split_kernel_sha256'] = split_kernel(pin,work,solid)
        for key, labels in [('held_stop',accepted),('failed_trial',inner)]:
            certificate = chart_orientation(axes,labels,opt.active)
            unit, grad, parts, _ = barrier.sample(labels)
            extrema = corner_sites(axes,labels,opt.active)
            elastic,_,_,info = solid.evaluate(labels,order,False)
            with mp.workdps(70):
                corners = {}
                for region, (det,cell,corner) in extrema.items():
                    row = local_corner(body,axes,labels,order,cell,corner)
                    if abs(det-row['chart_det']) > 1e-13:
                        raise ArithmeticError('Corner readback disagrees')
                    row['distance_from_old_corner'] = float(np.linalg.norm(coordinate(row['physical'])-old_xyz))
                    corners[region] = row
                previous_now = local_corner(body,axes,labels,order,old_cell,old_corner)
                minima = {}
                for region,worst in [('strictly_melted',info['worst_global']),('positive_order',info['worst_active'])]:
                    check = point(body,axes,labels,order,worst['cell'],worst['local'],_collar_reference)
                    if (max(abs(worst[a]-check[b])/max(1.,abs(worst[a])) for a,b in (
                        ('Jinv','Jinv'),('chi','chi'),('W','W_formal'),('grad_xi_norm','grad_xi_frobenius'))) > 5e-12
                        or (region=='strictly_melted') != (check['positive_vertices']==0)):
                        raise ArithmeticError('Independent native Gauss-point readback differs')
                    minima[region] = dict(worst,physical=check['physical'],positive_order_vertices=check['positive_vertices'])
            gauge_grad = gradient_witness(barrier,grad,order,axes,body)
            for region,row in corners.items():
                node = tuple(row['node'])
                flat = np.ravel_multi_index(node,order.shape)
                row['joined_unit_gauge_gradient_at_corner_vertex'] = []
                for c,cm in enumerate(barrier.components):
                    group = cm.groups[flat]
                    role = 'registered' if cm.held[group] else ('melt_only_free' if cm.dead[group] else 'active_trace_free')
                    row['joined_unit_gauge_gradient_at_corner_vertex'].append(
                        dict(component=c,role=role,signed=float(cm.sum(grad[...,c])[group])))
                row['distance_from_peak_melted_gauge_gradient'] = float(np.linalg.norm(
                    coordinate(row['physical'])-coordinate(gauge_grad['physical'])))
            result['states'][key] = dict(chart_certificate=certificate,minimum_one_sided_chart_corner=corners,
                original_corner_on_this_map=previous_now, native_sampled_physical_minima=minima, crystal_energy=elastic,
                gauge=dict(unit_energy=float(unit),**parts,peak_melt_only=gauge_grad))
    # Compare the projection diagnostic to the unmodified physical evaluator stored at the retained state.
    retained = result['states']['held_stop']
    profile = docs['held_stop']['physical']
    if (retained['native_sampled_physical_minima']['positive_order']['Jinv'] != profile['min_Jinv_active']
            or min(v['Jinv'] for v in retained['native_sampled_physical_minima'].values()) != profile['min_Jinv_global']
            or retained['crystal_energy'] != docs['held_stop']['elastic_energy']):
        raise ArithmeticError('Signed native measurements disagree with the retained-state record')
    failed = result['states']['failed_trial']
    if (failed['gauge']['corner_min'] != status['corner_min']
            or abs(failed['gauge']['peak_melt_only']['signed_unit_gradient']) != status['residual']):
        raise ArithmeticError('Signed melt measurements disagree with the inner journal')
    for state in result['states'].values():
        active = state['minimum_one_sided_chart_corner']['active']
        cert = state['chart_certificate']
        if cert['chart_lowest_cell'] != active['cell'] or abs(active['chart_det']-cert['chart_det_global_lower'])>5e-13:
            raise ArithmeticError('Minimum chart certificate does not resolve to the identified corner')
        if abs(state['minimum_one_sided_chart_corner']['melt']['chart_det']-state['gauge']['corner_min']) > 1e-14:
            raise ArithmeticError('Melt-gauge corner identification differs')
    modules = {Path(m.__file__).name for m in sys.modules.values() if getattr(m,'__file__',None)
               and Path(m.__file__).parent==pin}
    modules.update(('mouth_r4_eulerian_solid.cpp','mouth_free_shape_solid.cpp','mouth_cp10_spatial_crystal.cpp',
                    'mouth_r4_melt_barrier.cpp','mouth_r4_projected_held.py','mouth_r4_barrier_held.py'))
    frozen = read(area/'sources.json')
    result['pinned_source_sha256'] = {n:sha(pin/n) for n in sorted(modules)}
    if any(frozen[n] != v for n,v in result['pinned_source_sha256'].items()):
        raise ValueError('Pinned measurement source differs')
    result['script_sha256'] = sha(Path(__file__))
    with output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps({name:{'chart':st['chart_certificate'],'melt_corner':st['minimum_one_sided_chart_corner']['melt'],
        'active_corner':st['minimum_one_sided_chart_corner']['active'],'gauss':st['native_sampled_physical_minima'],
        'gauge':st['gauge']} for name,st in result['states'].items()},allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path)
    parser.add_argument('previous',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    run(args.area.resolve(),args.previous.resolve(),args.output.resolve())
