"""Move the CP10 inverse labels with material during a signed grid-depth update."""
import argparse
import ctypes as ct
import json
import subprocess
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import allnode_inverse, prisms
from mouth_cp10_spatial_crystal import Data, Pointer, integrate, load_context, signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def depth_transport(root, base, labels, rates, out, small_steps=False, split_diagonals=False):
    summary = read(base / 'result.json')
    assert summary['labels'] == digest(labels) and summary['checkpoint'] == CHECKPOINT_SHA
    assert digest(rates) == read(rates.parent / 'readback.json')['transport_sha256']
    library = base / 'crystal.so'
    if split_diagonals:
        out.parent.mkdir(exist_ok=True, parents=True)
        library = out.with_suffix('.so')
        if library.exists(): raise FileExistsError(library)
        source = Path(__file__).with_name('mouth_cp10_spatial_crystal.cpp')
        subprocess.run(['g++','-O3','-fPIC','-std=c++17','-shared','-DCP10_SPLIT_REGISTRY_DIAGONALS',str(source),'-o',str(library)],check=True)
    lib = ct.CDLL(str(library)).cp10_spatial_integral
    lib.argtypes = [ct.POINTER(Data), ct.c_int] + [Pointer] * 5 + [ct.c_int] + [Pointer] * 4
    lib.restype = ct.c_int
    ctx = load_context(root, labels)
    body, xi, chi = ctx[2], ctx[-2], ctx[-1]
    with np.load(root / 'r4_melted_survival01/coarse/state_00010.npz') as s:
        mesh = prisms(body, s['solid'])
    with np.load(labels) as f, np.load(rates) as rr:
        active = f['active_stencil']
        t = f['t']; dxi = rr['dxi_dt']; counts = rr['counts']
    native = read(root / 'r4_inverse_material01/cp10/result.json')['lagrangian_crystal']
    data = dict(checkpoint=CHECKPOINT_SHA, labels=digest(labels), rates=digest(rates),
                static_energy_levels=[], depth=[])
    out.parent.mkdir(exist_ok=True, parents=True)
    if out.exists(): raise FileExistsError(out)
    def save(): out.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    hs = (1.e-5, 2.e-6, 4.e-7, 8.e-8) if small_steps else (1.e-4, 5.e-5)
    if split_diagonals:
        hs=(1.e-4,5.e-5,1.e-5,2.e-6,4.e-7,8.e-8)
        bodies=[body]+[signed_geometry(ctx,k,h) for k in range(6) for h in (-1.e-4,1.e-4,-5.e-5,5.e-5)]
        vals, orient, _, _=integrate(lib,bodies,xi,chi,6)
        if orient[13]: raise ArithmeticError('Inverted base crystal')
        pairs=vals[7:].reshape(6,4)
        data['fixed_label_split_registry']=dict(energy=vals[0],delta=vals[0]-native,analytic=vals[1:7].tolist(),
                signed_energies=pairs.tolist(),force_full=((pairs[:,0]-pairs[:,1])/2.e-4).tolist(),force_half=((pairs[:,2]-pairs[:,3])/1.e-4).tolist())
        save()
    for h in (h*sign for h in hs for sign in (-1., 1.)):
        signed = signed_geometry(ctx, 0, h)
        axes = [signed.buffers['field_'+z] for z in ('x','t','a')]
        dt = axes[1]-t
        # Invert the unchanged native deformation at the *new* chart query
        # positions. These labels preserve the normalized material-wall map.
        position = np.stack(np.meshgrid(*axes,indexing='ij'), axis=-1)
        certain, nroot, duplicate, err = allnode_inverse(axes, position, mesh, active)
        if len(duplicate) or np.any(nroot[active] != 1): raise ArithmeticError('Nonunique transported material labels')
        transported = xi.copy(); transported[active] = certain[active]
        # Independent signed transport tangent; compare its *nodal* prediction
        # without using it for the physical crystal energy.
        deriv = np.where(dt[None,:,None,None]>0, dxi[:,:,:,1], dxi[:,:,:,0])
        valid = np.where(dt[None,:,None]>0, counts[:,:,:,1], counts[:,:,:,0])
        if np.any(active & (dt[None,:,None] != 0) & (valid == 0)): raise ArithmeticError('Missing one-sided transport derivative')
        pred = xi + deriv * dt[None,:,None,None]
        vals, orientation, _, _ = integrate(lib, [signed], transported, chi, 6)
        if orientation[13]: raise ArithmeticError('Transport candidate inverted ordered labels')
        row=dict(h=h, energy=vals[0], max_inverse_landing_error=err,
            max_nodal_transport_tangent_error=float(np.max(abs(certain[active]-pred[active]))),
            min_sampled_active_Jinv=orientation[15])
        if split_diagonals:
            if abs(h) in (1.e-4,5.e-5):
                idx=(0 if h<0 else 1)+(0 if abs(h)==1.e-4 else 2)
                row['same_motion_fixed_nodal_energy']=data['fixed_label_split_registry']['signed_energies'][0][idx]
            else:
                row['same_motion_fixed_nodal_energy']=float(integrate(lib,[signed],xi,chi,6)[0][0])
        data['depth'].append(row)
        save()
    es = [d['energy'] for d in data['depth']]
    signed = [(es[2*i]-es[2*i+1])/(2*h) for i,h in enumerate(hs)]
    half = signed[-1]
    projected_base = summary['levels'][1]['analytic'][0]
    if split_diagonals:
        projected_base=data['fixed_label_split_registry']['analytic'][0]
        se=[d['same_motion_fixed_nodal_energy'] for d in data['depth']]
        sf=[(se[2*i]-se[2*i+1])/(2*h) for i,h in enumerate(hs)]
        data['fixed_nodal_depth_by_h']=dict(zip(hs,sf))
        data['paired_transport_delta_by_h']=dict(zip(hs,[s-f for s,f in zip(signed,sf)]))
    data.update(material_matched_depth_force_by_h=dict(zip(hs,signed)),
                fixed_grid_nodal_depth=projected_base, transport_delta=half-projected_base)
    save()
    if small_steps or split_diagonals: return
    for q in (12, 16):
        total, orientation, _, _ = integrate(lib, [body], xi, chi, q)
        if orientation[13]: raise ArithmeticError('Static check found inverted ordered labels')
        data['static_energy_levels'].append(dict(cubature_degree=q, energy=total[0], delta=total[0]-native,
            min_active_Jinv=orientation[15]))
        save()


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    for key in ('root','base','labels','rates','out'): args.add_argument(key,type=Path)
    args.add_argument('--small-steps',action='store_true')
    args.add_argument('--split-diagonals',action='store_true')
    a=args.parse_args(); depth_transport(*(getattr(a,key).resolve() for key in ('root','base','labels','rates','out')),small_steps=a.small_steps,split_diagonals=a.split_diagonals)
