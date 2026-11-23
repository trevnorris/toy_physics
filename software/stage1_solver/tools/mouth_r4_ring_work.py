"""Signed physical ring probes and a fixed-geometry energy budget for held R4 fields."""
import argparse
import ctypes as ct
import json
from math import fsum
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from mouth_whole_r4_order_edges import digest, read


def ring_kernel(include_dir, build_dir):
    from mouth_free_shape_solid import Data, Pointer
    source = Path(__file__).with_suffix('.cpp')
    build_dir.mkdir(exist_ok=True, parents=True)
    library = build_dir / (digest(source) + '-' + digest(include_dir/'mouth_free_shape_solid.cpp') + '.so')
    if not library.exists():
        with tempfile.TemporaryDirectory(dir=build_dir) as temp:
            stage = Path(temp)/'work.so'
            subprocess.run(['g++', '-std=c++17', '-O3', '-shared', '-fPIC',
                            '-I', str(include_dir), str(source), '-o', str(stage)], check=True)
            stage.replace(library)
    work = ct.CDLL(str(library)).radial_ring_work
    work.argtypes = [ct.POINTER(Data), Pointer, ct.c_int, Pointer]
    work.restype = ct.c_int
    return work


def measure_ring(solid, q, correction, include_dir, build_dir):
    """Each free angular site on the last radial/depth ring moves outward by h.

    The terminal cells lie outside the support of the polar chart deformation,
    so the chart radial motion is exactly h/b(z) at fixed depth and polar angle.
    Both absolute chi (with a zero solid base) and reference-plus-correction work.
    """
    from mouth_free_shape_solid import pointer
    work = ring_kernel(include_dir, build_dir)
    offsets = solid.unpack(q)
    times = solid.t[-2]+offsets[-2, -2, :, 1]
    if min(times) <= .62 or not np.all(solid.free[-2, -2, :, 0]):
        raise ValueError('Physical radial probe requires the free terminal ring outside the polar bump')
    b = solid.sh.radius(solid.sh.length*times)[0]
    unit = np.zeros_like(offsets)
    unit[-2, -2, :, 0] = 1./b
    lift = np.ascontiguousarray(unit[..., 0])
    direction = solid.pack(unit)
    at_nodes = np.ascontiguousarray(correction[solid.ids], dtype=float)
    # Full-domain solver gradient on every material coordinate is independent
    # of the local star restriction and the constitutive product split.
    total, native_g, _, info = solid.evaluate(q, correction, load=False)
    full = -float(native_g@direction)
    solver_ring = solid.residual(q, native_g)['last_free_ring_angular_sum']
    def probe(h, differentiate=False):
        delta = solid.unpack(q+h*direction)
        solid.data.offset, solid.data.correction = pointer(delta), pointer(at_nodes)
        solid.data.diagnose = 0
        values = np.zeros(8)
        code = work(ct.byref(solid.data), pointer(lift), int(differentiate), pointer(values))
        if code or not np.isfinite(values).all():
            raise ArithmeticError(f'Invalid solid probe at {h}, error {code}')
        return values

    raw = probe(0., True)
    names = ('stretch_squared', 'minus_log_J', 'quadratic_compression', 'sampling_native_order',
             'sampling_order_argument', 'quadrature_point_and_weight_motion')
    forces = dict(zip(names, (-raw[2:]).tolist()))
    probes = []
    for h in (2.e-5, 1.e-5, 5.e-6, 2.5e-6):
        minus, plus = probe(-h)[:2], probe(h)[:2]
        signed = (minus-plus)/(2*h)
        probes.append(dict(outward_step=h, inward_energy=minus.tolist(), outward_energy=plus.tolist(),
                           radial_force=float(signed[0]), argument_weighted_force=float(signed[1]),
                           native_weighted_force=float(signed[0]-signed[1])))
    checks = dict(full_vs_solver_radial=abs(full-solver_ring[0])<2.e-11,
                  work_terms_vs_full=abs(fsum(forces.values())-full)<2.e-8,
                  smallest_probe_vs_full=abs(probes[-1]['radial_force']-full)<2.e-5,
                  last_halving=abs(probes[-1]['radial_force']-probes[-2]['radial_force'])<2.e-5)
    return dict(physical_sign='positive outward', mesh_reference=[solid.t[-2], solid.x[-2]],
                landed_depth_fraction=[times.min(),times.max()], angular_sites=len(times),
                probe_energy_columns=['chi_W_entire_affected_star','field_argument_times_W_affected_star'],
                full_domain_chi_W=total, full_domain_argument_times_W=info['transfer_energy'],
                full_domain_radial_force=full, full_domain_radial_axial=solver_ring,
                constitutive_radial_forces=forces, constitutive_sum=fsum(forces.values()), signed_probes=probes,
                direct_ring_force=dict(chi_W=full, order_gradient=0., double_well=0., carrier=0., phase=0.),
                direct_field_terms_basis='Exact at fixed Eulerian field coefficients and chart: no solid coordinate enters these terms.',
                checks=checks)


def coherent_energy(model, field):
    """Excess over the same vortex and reference order used by Joined.evaluate."""
    from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA, WELL
    v, p, dc = field
    wt, F, c0 = model.D[:,0], model.D[:,4], model.ref[:,0]
    if np.any(model.D[:,10]):
        raise ValueError('The material deposition would be counted twice')
    chi = c0+model.V@dc
    gv, gp, gc = (np.array([g@c for g in model.G]) for c in (v,p,dc))
    ref_grad = model.ref[:,1:].T
    gc += ref_grad
    def square(a):
        return np.einsum('ij,ij->j', a, a)
    y = (model.V@v)**2
    n = F**2*y
    gp_dot_gc = np.einsum('ij,ij->j',gp,gc)
    return dict(order_gradient=float(.5*KAPPA*(wt@(square(gc)-square(ref_grad)))),
                double_well=float(WELL*(wt@((chi*(1-chi))**2-(c0*(1-c0))**2))),
                density_amplitude_gradient=float(.5*(wt@(F**2*square(gv)))),
                density_potential=float(wt@(F**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20)),
                meridional_phase=float(.5*(wt@(n*square(gp)))),
                phase_order_coupling=float(-.5*ORDER_INERTIA*(wt@(n*gp_dot_gc**2/(n+ORDER_INERTIA*square(gc))))),
                vortex_boundary=float(.5*((v-1)@(model.B@(v+1)))))


def joint_probe(b, q, u, external_load=0.):
    m = b.model
    original_load = m.material_load
    # Exercise the full joint energy, including the change in the material
    # source installed in the field action. The chart and coefficients stay
    # fixed at both physical ring displacements.
    offsets = b.unpack(q)
    z = b.sh.length*(b.t[-2]+offsets[-2,-2,:,1])
    displacement = np.zeros_like(offsets)
    displacement[-2,-2,:,0] = 1./b.sh.radius(z)[0]
    unit = b.pack(displacement)
    joint = []
    h = 1.e-5
    for sign in (-1.,1.):
        reference, _, signed_load, _ = b.evaluate(q+sign*h*unit,gradient=False)
        m.material_load = signed_load+external_load
        native = m.evaluate(u,True)[0]
        signed_terms = coherent_energy(m,u)
        joint.append(dict(physical_displacement=sign*h, native_plus_coupled_field=reference+native,
                          field_without_deposition=native-float(m.material_load@u[2]),
                          explicit_field_terms=signed_terms))
    m.material_load = original_load
    field_forces = {k: (joint[0]['explicit_field_terms'][k]-joint[1]['explicit_field_terms'][k])/(2*h)
                    for k in joint[0]['explicit_field_terms']}
    joint_force = (joint[0]['native_plus_coupled_field']-joint[1]['native_plus_coupled_field'])/(2*h)
    return joint, field_forces, joint_force


def emit(summary, output):
    output.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    for key, row in summary['preparations'].items():
        ring = row['ring']
        print(key, 'radial', ring['full_domain_radial_force'], ring['constitutive_radial_forces'],
              'smallest_probe', ring['signed_probes'][-1]['radial_force'],
              'total_energy', row['total_energy'])
    print('wider_minus_narrow_at_fixed_q', json.dumps(summary['wider_minus_narrow_at_fixed_q']))
    print('checks', summary['accepted'], flush=True)


def run(args, include_dir):
    from mouth_free_shape_solid import Solid
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_mobile_joined_fields import JointFields
    from mouth_rim_pinned_comparison import validate_guard

    folder = args.run/'fine'
    validate_guard(folder)
    validate_guard(args.wider/'..')
    end = read(folder/'endpoint.json')
    measurement = Path('research/pde_ledger_v3/_measurements')
    wide_path = measurement/'mouth_whole_r4_wider_order.json'
    wide = read(wide_path)
    if read(args.wider/'readback.json')['comparisons'] != wide['comparisons'] or not wide['accepted']:
        raise ValueError('Saved wider comparison and accepted measurement differ')
    inv = read(folder/'invocation.json')['command']
    state = Path(inv[inv.index('-u')+2])
    trace = folder/end['file']
    controls = folder/'endpoint_controls.npz'
    if any(digest(p)!=wide[k] for p,k in ((state,'preparation_sha256'),(trace,'state_sha256'),
              (controls,'controls_sha256'),(args.run/'sources.json','source_manifest_sha256'))):
        raise ValueError('Held endpoint and wider inputs differ')
    with np.load(trace) as data, np.load(controls) as rest:
        x,q = data['x'],data['q']
        fields = dict(narrow_wound=data['field'],narrow_rest=rest['rest'])
    wide_states = {name: args.wider/(name+'.npz') for name in ('wound','rest')}
    for name,path in wide_states.items():
        with np.load(path) as data:
            fields['wider_'+name] = data['field']
    f = JointFields(Front(state,GRIDS['double'],1),state,x,precondition=False)
    m = f.model
    args_solid = {k: int(inv[inv.index('--'+k)+1]) for k in ('n','k','angles','degree','subdivisions')}
    b = Solid(state, m.shape, m.polar_shift, args_solid['n'],args_solid['k'],args_solid['angles'],m,
              degree=args_solid['degree'],subdivisions=args_solid['subdivisions'])
    base_energy, _, load, _ = b.evaluate(q,gradient=False)
    m.material_load = load
    summary = dict(command=sys.argv, source_hashes={str(p): digest(p) for p in (
        Path(__file__),Path(__file__).with_suffix('.cpp'),args.run/'sources.json',
        include_dir/'mouth_free_shape_solid.cpp', wide_path, state, trace, controls, *wide_states.values())},
        convention='negative energy derivative per unit uniform physical outward ring motion; held field chart',
        stationary_preparations='saved fields, no fresh minimization', static_geometry_energy=f.static,
        native_order_solid_energy=base_energy, preparations={}, wider_minus_narrow_at_fixed_q={})
    for label, u in fields.items():
        width,name = label.split('_')
        ring = measure_ring(b,q,u[2],include_dir,args.build)
        terms = coherent_energy(m,u)
        material_energy = ring['full_domain_chi_W']
        energy = fsum([f.static,material_energy,*terms.values()])
        stored = wide['narrow_energy'][name] if width=='narrow' else wide['comparisons'][name]['solver']['energy']
        archived_ring = wide['comparisons'][name][width+'_solid']['last_free_ring_angular_sum']
        legacy = m.evaluate(u,True)[0]
        checks = dict(native_field_energy=abs(legacy-stored)<2.e-9,
                      separate_field_energy=abs(fsum([*terms.values(),load@u[2],-legacy]))<2.e-9,
                      coupled_full_energy=abs(energy-(wide['static_energy']+stored))<2.e-8,
                      archived_ring=abs(ring['full_domain_radial_force']-archived_ring[0])<2.e-11,
                      linear_order_deposition=abs(material_energy-base_energy-load@u[2])<2.e-8)
        joint, field_forces, joint_force = joint_probe(b,q,u)
        checks['full_signed_probe'] = abs(joint_force-ring['full_domain_radial_force'])<2.e-5
        checks['signed_fixed_chart_fields'] = all(f==0. for f in field_forces.values())
        summary['preparations'][label]=dict(ring=ring, field_excess_terms=terms, total_energy=energy,
                    joint_signed_energies=joint, full_joint_probe_force=joint_force,
                    measured_direct_field_forces=field_forces, checks={k: bool(v) for k,v in checks.items()})
        (args.build/(label+'.json')).write_text(json.dumps(summary['preparations'][label],allow_nan=False)+'\n')
    for name in ('wound','rest'):
        old, new = [summary['preparations'][width+'_'+name] for width in ('narrow','wider')]
        split = {k: new['field_excess_terms'][k]-v for k,v in old['field_excess_terms'].items()}
        split['chi_W'] = new['ring']['full_domain_chi_W']-old['ring']['full_domain_chi_W']
        dF = {k: new['ring']['constitutive_radial_forces'][k]-v for k,v in old['ring']['constitutive_radial_forces'].items()}
        summary['wider_minus_narrow_at_fixed_q'][name] = dict(energy=split, energy_sum=fsum(split.values()),
             direct_ring_force_change=dF, ring_force_change=fsum(dF.values()))
    summary['accepted'] = all(all(row['checks'].values()) and all(row['ring']['checks'].values())
                              for row in summary['preparations'].values())
    emit(summary, args.output)
    if not summary['accepted']:
        raise ArithmeticError('Readback checks failed; see reported probe and energy differences')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run','wider','build','output'):
        p.add_argument(key,type=Path)
    args = p.parse_args()
    pin = args.run/'pin'
    for key,sha in read(args.run/'sources.json').items():
        if digest(pin/key)!=sha:
            raise ValueError(f'Changed pinned input {key}')
    include_dir = pin/'software/stage1_solver/tools'
    sys.path.insert(0,str(include_dir))
    run(args,include_dir)
