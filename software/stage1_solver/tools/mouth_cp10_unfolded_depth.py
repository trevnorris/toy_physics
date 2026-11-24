"""Matched signed depth derivative with an independently unfolded frozen crystal.

The physical mouth, order, reservoirs and phase agree with checkpoint 10. Only
the crystal material map differs. Each signed spatial grid obtains its labels
by direct inversion of that control map. The six disjoint native polar sectors
partition the original composite Gauss quadrature without changing its points.
"""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_label_grid import allnode_inverse, prisms
from mouth_cp10_nested_forces import subdivide
from mouth_cp10_spatial_crystal import signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_free_shape_solid import Data, Pointer, pointer
from mouth_r4_unfolded_control import load_control_context


H = 2e-7
LEVELS = (1,2,4)
DEGREES = (6,8)
TOOLS = Path(__file__).resolve().parent
KERNEL = TOOLS/'mouth_cp10_matched_depth.cpp'


def read(p):
    return json.loads(p.read_bytes())


def write_new(path, value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write('\n')


def control_context(root, area):
    checkpoint = root/'r4_melted_survival01/coarse/state_00010.npz'
    if digest(checkpoint)!=CHECKPOINT_SHA:
        raise ValueError('Changed checkpoint-10 physical state')
    ctx = load_control_context(root,area)
    with np.load(area/'control.npz', allow_pickle=False) as control, np.load(checkpoint, allow_pickle=False) as original:
        if str(control['parent_sha256'])!=CHECKPOINT_SHA:
            raise ValueError('Wrong parent state')
        for name in ('shape_params','shape_extra','x','field','radius','winding'):
            if not np.array_equal(control[name],original[name]):
                raise ValueError(f'Physical data changed: {name}')
        solid = np.ascontiguousarray(control['solid'])
        if (np.array_equal(solid,original['solid']) or solid.shape!=ctx[2].chart_shape or
                not np.array_equal(ctx[2].unpack(control['q']),solid)):
            raise ValueError('Unfolded solid and packed material coordinates disagree')
    return ctx,solid


def prepare(root,area):
    ctx,solid=control_context(root,area)
    state=area/'control.npz'
    dest=area/'depth';dest.mkdir(exist_ok=False)
    mesh=prisms(ctx[2],solid)
    base=ctx[2].buffers['field_t']
    record=dict(source_sha256=digest(state),parent_sha256=CHECKPOINT_SHA,h=H,grids={})
    for level in LEVELS:
        path=area/str(level)/'labels.npz'
        if digest(path)!=read(path.parent/'result.json')['sha256']:
            raise ValueError('Unfolded label input changed')
        with np.load(path, allow_pickle=False) as f:
            if str(f['source_sha256'])!=record['source_sha256']:
                raise ValueError('Labels belong to a different crystal')
            x,t,a=(np.ascontiguousarray(f[k]) for k in 'xta')
            order=np.ascontiguousarray(f['order'])
            active=f['active_stencil']
            xi=f['labels']
            rates=(f['dxi_dt_minus'],f['dxi_dt_plus'])
        for key,axis in zip('xta',(x,t,a)):
            if not np.array_equal(axis,subdivide(ctx[2].buffers['field_'+key],level)):
                raise ValueError(f'Changed physical {key} knots')
        save=dict(x=x,a=a,order=order); signs=[]
        for side,sign in enumerate((-1,1)):
            moved=signed_geometry(ctx,0,sign*H)
            lgrid=moved.buffers['field_t']
            z=subdivide(lgrid,level)
            pos=np.stack(np.meshgrid(x,z,a,indexing='ij'),axis=-1)
            actual,count,ambiguous,landing=allnode_inverse((x,z,a),pos,mesh,active)
            if ambiguous or (count[active]!=1).any() or not np.isfinite(actual[active]).all():
                raise ArithmeticError('Control inverse is missing or multiple')
            dt=subdivide(lgrid-base,level)[None,:,None,None]
            vel=np.where(dt>=0,rates[1],rates[0])
            linear=xi+np.multiply(dt,vel,out=np.zeros_like(xi),where=dt!=0.)
            discrepancy=float(np.abs(linear[active]-actual[active]).max())
            mapped=xi.copy();mapped[active]=actual[active]
            tag='minus' if sign<0 else 'plus'
            save[tag+'_labels']=mapped;save[tag+'_knots']=z
            signs.append(dict(sign=sign,length=moved.data.length,active_nodes=int(active.sum()),
                              maximum_landing_error=landing,max_exact_minus_tangent=discrepancy))
        target=dest/f'signed_{level}.npz'
        np.savez_compressed(target,**save)
        record['grids'][str(level)]=dict(static_labels_sha256=digest(path),signed_sha256=digest(target),signs=signs)
    write_new(dest/'prepared.json',record)
    print(json.dumps(record,allow_nan=False),flush=True)


def sector(root,area,k):
    ctx,solid=control_context(root,area)
    body=ctx[2];body.data.offset=pointer(solid)
    dest=area/'depth'; prep=read(dest/'prepared.json')
    if prep['source_sha256']!=digest(area/'control.npz') or prep['h']!=H:
        raise ValueError('Unexpected control input')
    work=dest/f'polar_{k}';work.mkdir(exist_ok=False)
    lib=work/'kernel.so'
    subprocess.run(['g++','-O3','-fPIC','-shared','-std=c++17',str(KERNEL),'-o',str(lib)],check=True)
    fn=ct.CDLL(str(lib)).matched_depth
    fn.restype=ct.c_int
    fn.argtypes=([ct.POINTER(Data),Pointer,ct.c_int,Pointer,ct.c_int,ct.c_int]+[Pointer]*5+
                 [ct.c_double]*3+[Pointer]*2+[ct.c_int]*3+[Pointer])
    if len(body.a)!=7 or k not in range(6):
        raise ValueError('Expected full coverage by six native polar sectors')
    report=dict(source_sha256=prep['source_sha256'],parent_sha256=CHECKPOINT_SHA,kernel_sha256=digest(KERNEL),
                h=H,sector=k,range=body.a[k:k+2].tolist(),rows={})
    for level in LEVELS:
        path=dest/f'signed_{level}.npz'; spec=prep['grids'][str(level)]
        if digest(path)!=spec['signed_sha256']: raise ValueError('Signed fields have changed')
        with np.load(path,allow_pickle=False) as f:
            x,a=f['x'],f['a']
            tm,tp=f['minus_knots'],f['plus_knots']
            ym,yp=f['minus_labels'],f['plus_labels']; chi=f['order']
        for degree in DEGREES:
            nodes,weights=np.polynomial.legendre.leggauss(degree)
            gp,gw=np.ascontiguousarray((1+nodes)/2),np.ascontiguousarray(weights/2)
            vals=np.zeros(10)
            rc=fn(ct.byref(body.data),pointer(x),len(x),pointer(a),len(a),len(tm),
                  pointer(tm),pointer(tp),pointer(ym),pointer(yp),pointer(chi),
                  *(s['length'] for s in spec['signs']),H,pointer(gp),pointer(gw),degree,k,k+1,pointer(vals))
            if rc or not np.isfinite(vals).all() or vals[9]<=0 or vals[7]<=0:
                raise ArithmeticError(f'Ordered signed kernel failed: {rc}, {vals.tolist()}')
            em,ep=vals[:2],vals[2:4];forces=vals[4:6]
            if max(abs((em-ep)/(2*H)-forces))>1e-8:
                raise ArithmeticError('Energy sum and paired local work differ')
            report['rows'][f'{level}_{degree}']=dict(energy_minus=em.tolist(),energy_plus=ep.tolist(),
                forces=forces.tolist(),shared_points=int(vals[7]),sampled_positive=int(vals[8]),min_ordered_Jinv=vals[9])
    write_new(work/'result.json',report)
    print(json.dumps(report,allow_nan=False),flush=True)


def combine(root,area):
    dest=area/'depth';prep=read(dest/'prepared.json')
    groups=[read(dest/f'polar_{k}/result.json') for k in range(6)]
    if prep['source_sha256']!=digest(area/'control.npz'):raise ValueError('Wrong control state')
    archive=read(TOOLS.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_matched_knot.json')
    common=archive['common_bath_depth_force']
    report=dict(source_sha256=prep['source_sha256'],parent_sha256=CHECKPOINT_SHA,
                prepared_sha256=digest(dest/'prepared.json'),h=H,common_bath_force=common,rows={})
    for k,data in enumerate(groups):
        oc=read(area/f'depth_guard_{k}/outcome.json')
        if (oc['exitCode'] or oc['childOutcome']['exitCode'] or not oc['limitsVerified'] or k!=data['sector'] or
                data['source_sha256']!=prep['source_sha256'] or data['kernel_sha256']!=digest(KERNEL)):
            raise ValueError(f'Incomplete or mismatched sector {k}')
    for level in LEVELS:
        for degree in DEGREES:
            key=f'{level}_{degree}'
            parts=[d['rows'][key] for d in groups]
            sums={name:np.sum([p[name] for p in parts],axis=0,dtype=np.longdouble) for name in ('forces','energy_plus','energy_minus')}
            f=sums['forces'];sigma=(sums['energy_minus']-sums['energy_plus'])/(2*H)
            if max(abs(f-sigma))>1e-8: raise ArithmeticError('Combined energy/work discrepancy')
            report['rows'][key]=dict(native_total=float(f[0])+common,spatial_total=float(f[1])+common,
                paired_crystal_forces=f.astype(float).tolist(),spatial_minus_native=float(f[1]-f[0]),
                reconstructed_fd_minus_paired=(sigma-f).astype(float).tolist(),
                signed_energies=[[float(e) for e in sums[name]] for name in ('energy_minus','energy_plus')],
                min_ordered_Jinv=min(p['min_ordered_Jinv'] for p in parts),points=sum(p['shared_points'] for p in parts))
    write_new(dest/'result.json',report)
    print(json.dumps(report,allow_nan=False),flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('action',choices=('prepare','sector','combine'))
    cli.add_argument('root',type=Path);cli.add_argument('area',type=Path)
    cli.add_argument('sector',nargs='?',type=int)
    args=cli.parse_args(); root=args.root.resolve(); area=args.area.resolve()
    if args.action=='sector':
        sector(root,area,args.sector)
    else:
        {'prepare':prepare,'combine':combine}[args.action](root,area)
