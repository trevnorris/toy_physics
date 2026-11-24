"""Frozen CP10 physical-domain motion and production fixed-nodal shape partials."""
import argparse
from copy import copy
import gc
import inspect
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_cp10_field_shape import Bath, Shape, axes, bump
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, acquired, digest, emit, read


def coordinates(model, zaxes, radius, shift, weight):
    """Physical cylindrical coordinates for each full finite-element block."""
    out = {}
    for key in model.ids:
        slab, side = key
        x, theta, z = np.meshgrid(model.rr[side],model.a,zaxes[slab],indexing='ij')
        b = radius(z)[0]
        rad = b + x*(28-b) if side else x*b
        factor = 1-x if side else x**3
        rad += shift*weight(z)*factor*(3*np.cos(theta)**2-1)/2
        out[key] = np.stack((rad*np.sin(theta),rad*np.cos(theta),z),axis=-1)
    return out


def boundary_groups(blocks):
    def combine(*faces):
        return np.concatenate([face.reshape(-1,3) for face in faces])
    return dict(
        all_six_blocks=combine(*blocks.values()),
        entrance_w0=combine(*(blocks[s,b][:,:,-1 if s==0 else 0] for s in (0,1) for b in (0,1))),
        clamped_entrance_ring=combine(blocks[1,0][-1,:,0]),
        lower_slab_face=combine(*(blocks[0,b][:,:,0] for b in (0,1))),
        exit_disk=combine(blocks[1,0][:,:,-1],blocks[2,0][:,:,0]),
        exit_rim=combine(*(blocks[s,b][0 if b else -1,:,-1 if s==1 else 0] for s in (1,2) for b in (0,1))),
        roof_bulk=combine(*(blocks[2,b][:,:,-1] for b in (0,1))),
        radial_brane_cutoff_and_feed_port=combine(blocks[0,1][-1]),
        radial_bulk_cutoff=combine(blocks[1,1][-1],blocks[2,1][-1]),
        whole_bulk_port=combine(blocks[1,1][-1],blocks[2,1][-1],blocks[2,1][:,:,-1],blocks[2,0][:,:,-1]),
        two_lateral_wall_faces=combine(blocks[1,0][-1],blocks[1,1][0]),
        polar_axis=combine(*(a[:,0] for a in blocks.values())),
        equator=combine(*(a[:,-1] for a in blocks.values())),
        radial_origin=combine(*(blocks[s,0][0] for s in range(3))))


def model_inputs(sc):
    return acquired(sc/'r4_melted_survival01', sc/'r4_exact_support01/patch')


def boundaries(sc, out):
    build,x,*_ = model_inputs(sc)
    from mouth_joined_angular_work import bump as native_bump
    # Shape construction assembles the full legacy FE operator. Retain just
    # its physical geometry so that each signed operator can be released.
    def geometry(state):
        m = build(state).model
        return SimpleNamespace(**{k:getattr(m,k) for k in
                ('ids','rr','a','z','shape','polar_shift','feed','bulk')})
    m0 = geometry(x)
    gc.collect()
    ref = inspect.getclosurevars(build).nonlocals['owner'].reference_depth
    h=1.e-5
    snapshots=[]
    for step in (0.,-h,h):
        off=np.array([step,0.,0.,0.,0.,0.])
        m=m0 if step==0 else geometry(x+off)
        gc.collect()
        if not all(np.array_equal(m.ids[k],m0.ids[k]) for k in m0.ids):
            raise ArithmeticError('Joined topology changed')
        feed=np.unique(m0.ids[0,1][-1]); bulk=np.unique(np.r_[m0.ids[1,1][-1].ravel(),
                          m0.ids[2,1][-1].ravel(),m0.ids[2,1][:,:,-1].ravel(),m0.ids[2,0][:,:,-1].ravel()])
        if not np.array_equal(feed,m.feed) or not np.array_equal(bulk,m.bulk):
            raise ArithmeticError('Reservoir masks changed')
        sh=Shape(m0.shape,m0.polar_shift,off)
        ez=axes(sh.length,ref,m0)
        old=coordinates(m,m.z,m.shape.radius,m.polar_shift,lambda z:native_bump(z,m.shape.length)[0])
        new=coordinates(m0,ez,sh.radius,sh.shift,lambda z:bump(z/sh.length,sh.length)[0])
        snapshots.append((boundary_groups(old),boundary_groups(new)))
    report={}
    for name in snapshots[0][0]:
        old,new=([a[i][name] for a in snapshots] for i in (0,1))
        within=[value[2]-value[1] for value in (old,new)]
        report[name]=dict(samples_per_shape=len(old[0]),
            max_coordinate_delta=float(max(np.max(abs(a-b)) for a,b in zip(old,new))),
            max_difference_of_signed_displacements=float(np.max(abs(within[0]-within[1]))),
            max_difference_of_central_velocity=float(np.max(abs(within[0]-within[1]))/(2*h)),
            max_magnitude_of_legacy_velocity=float(np.max(np.linalg.norm(within[0],axis=1))/(2*h)))
    if max(m['max_coordinate_delta'] for m in report.values())>5.e-14:
        raise ArithmeticError('Physical charts diverged')
    emit(out,dict(checkpoint_sha256=CHECKPOINT_SHA,log_depth_step=h,
             feed_nodes=len(feed),bulk_nodes=len(bulk),groups=report,
             physical_support='Depth varies the exit and wall; nodes can slide tangentially on the fixed outer radial cutoff. Both charts share every displacement. Counts include both traces where specified.'))


def work(sc, out):
    arc=Path('research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json')
    source=read(arc)
    if source['checkpoint_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Unexpected input state')
    commands=source['commands']; archive=commands['split_registry_total']
    reused=Path(archive['result'])
    if digest(reused)!=archive['result_sha256']:
        raise ValueError('Prior force table changed')
    old=read(reused)
    fields=sc/'r4_cp10_two_metric_fields01/field.json'
    registry=sc/'r4_cp10_spatial_crystal01/transport/registry.json'
    projected=sc/'r4_cp10_spatial_crystal01/readout/result.json'
    prior=Path('research/pde_ledger_v3/_measurements/mouth_r4_eulerian_reference.json')
    for path in (fields,registry,projected,prior):
        if digest(path)!=old['source_sha256'][path.name]:
            raise ValueError(f'Changed force input {path}')
    for folder in (sc/'r4_cp10_two_metric_fields01/guard',sc/'r4_cp10_spatial_crystal01/registry_guard',
                   sc/'r4_cp10_two_metric01/guard_registry_total'):
        outcome=read(folder/'outcome.json'); limits=read(folder/'limit-validation.json')
        if outcome['exitCode'] or not outcome['limitsVerified'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
            raise ValueError(f'Invalid previous guard {folder}')
    # Five independent single central differences at h=1e-4 already pass. The
    # depth derivative needs a shorter step to resolve its O(h^2) stream term.
    # Preserve its completed signed crystal energies at exactly that step.
    h=1.e-5
    data=read(registry)
    pair=[v['same_motion_fixed_nodal_energy'] for v in data['depth'] if abs(v['h'])==h]
    if len(pair)!=2 or [v['h'] for v in data['depth'] if abs(v['h'])==h]!=[-h,h]:
        raise ValueError('Missing matched depth central pair')
    build,x,q,u,_,_=model_inputs(sc)
    trial=build(x)
    owner=inspect.getclosurevars(build).nonlocals['owner']
    bath=Bath(trial,owner.base,u,owner.reference_depth)
    # Use the original quadrature and coefficients; promote both sides' coordinates
    # and contractions before subtracting the large far-bath energy.
    m=bath.model=copy(bath.model)
    cast=lambda a: np.asarray(a,dtype=np.longdouble)
    m.a=cast(m.a); m.rr=tuple(map(cast,m.rr)); m.z=tuple(map(cast,m.z))
    bath.u=cast(bath.u)
    unit=np.array([1.,0.,0.,0.,0.,0.],dtype=np.longdouble)
    minus,plus=(bath.energy(sign*np.longdouble(h)*unit) for sign in (-1,1))
    components=(minus-plus)/(2*np.longdouble(h))
    analytic=np.array(old['analytic_total'])
    direct=np.array(old['central_total']['0.0001']);direct[0]=float(components.sum())+(pair[0]-pair[1])/(2*h)
    residual=analytic-direct
    item=dict(checkpoint_sha256=CHECKPOINT_SHA,label_sha256=old['labels'],
         archive={str(p):digest(p) for p in (reused,fields,registry,projected,prior)},
         convention='Nodal label and order coefficients held on the production moving chart, epsilon=0. No re-inversion of the labels.',
         modes=old['modes'],analytic=analytic.tolist(),steps=[h]+[1.e-4]*5,
         depth_bath_components=list(map(float,components)),depth_signed_crystal_energies=pair,
         central=direct.tolist(),analytic_minus_central=residual.tolist(),
         max_abs_direct=float(np.max(abs(residual))),raw_archive=old['central_total'],
         archive_richardson=old['richardson'],max_abs_analytic_minus_richardson=float(np.max(np.abs(old['analytic_minus_richardson']))))
    emit(out,item)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind',choices=('boundaries','work'))
    parser.add_argument('sc',type=Path)
    parser.add_argument('output',type=Path)
    opts=parser.parse_args()
    {'boundaries':boundaries,'work':work}[opts.kind](opts.sc.resolve(),opts.output.resolve())
