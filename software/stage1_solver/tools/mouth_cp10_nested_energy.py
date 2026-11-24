"""Read-only nested spatial-label energy at the native CP10 physical material points."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, acquired, digest

GROUPS = ('same_material_prism', 'other_prism_same_bond_side', 'opposite_bond_side')
COLUMNS = ('positive_order_points', 'saved', 'spatial', 'saved_label_with_spatial_gradient', 'absolute_pointwise_energy_error')


def summarize(rows):
    count, old, new, hybrid, absolute = np.sum(rows.reshape(-1,len(COLUMNS)),axis=0,dtype=np.longdouble)
    return dict(count=int(count),original=float(old),projected=float(new),difference=float(new-old),
                absolute_pointwise_sum=float(absolute),gradient_first=float(hybrid-old),
                label_values_and_metric_second=float(new-hybrid))


def acquire_labels(path, body, factor):
    with np.load(path, allow_pickle=False) as d:
        if str(d['source_sha256'])!=CHECKPOINT_SHA:
            raise ValueError('Refined labels have the wrong checkpoint source')
        axes = [np.ascontiguousarray(d[k]) for k in ('x','t','a')]
        xi = np.ascontiguousarray(d['labels'])
        stencil = d['active_stencil']
        for k,new in zip(('x','t','a'),axes):
            old = body.buffers['field_'+k]
            if len(new)!=factor*(len(old)-1)+1 or not np.array_equal(new[::factor],old):
                raise ValueError('Nested grid does not retain original '+k+' knots')
        if xi.shape!=tuple(map(len,axes))+(3,) or not np.isfinite(xi[stencil]).all():
            raise ValueError('Invalid labels on original positive-order supports')
    return axes,xi


def calculate(root,labels,factor,out):
    build,state,q,u,_,_ = acquired(root/'r4_melted_survival01', root/'r4_exact_support01/patch')
    trial = build(state)
    body = trial.b
    from mouth_free_shape_solid import Data, Results, Pointer, pointer
    axes,xi = acquire_labels(labels,body,factor)
    out.mkdir(parents=True,exist_ok=False)
    lib = out/'nested_image.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(Path(__file__).with_suffix('.cpp')),'-o',str(lib)],check=True)
    kernel = ct.CDLL(str(lib)).nested_image_energy
    kernel.argtypes = [ct.POINTER(Data),ct.POINTER(Results)]+[Pointer,ct.c_int]*3+[Pointer]*4
    kernel.restype = ct.c_int
    table = np.zeros((len(body.t)-1,len(body.x)-1,len(body.a)-1,2,len(GROUPS),len(COLUMNS)))
    stats = np.zeros(4)
    witness = np.zeros(12)
    knot_args = [value for knots in axes for value in (pointer(knots),len(knots))]
    body.work = lambda data, results: kernel(data,results,*knot_args,pointer(xi),pointer(table),pointer(stats),pointer(witness))
    try:
        native = body.evaluate(q,u[2],load=False,gradient=False)[0]
    except (ValueError, ArithmeticError):
        (out/'invalid.json').write_text(json.dumps(dict(physical=witness[:3].tolist(),original_label=witness[3:6].tolist(),
            projected_label=witness[6:9].tolist(),chi=witness[9],label_chart_det=witness[10],native_forward_chart_det=witness[11]),indent=2)+'\n')
        raise
    reference=json.loads((root/'r4_cp10_two_metric01/image/result.json').read_bytes())
    if reference['checkpoint']!=CHECKPOINT_SHA:
        raise ValueError('Wrong 1x identical-point reference')
    stats_all=summarize(table)
    if stats_all['count']!=5291652 or abs(stats_all['original']-reference['native'])>2e-14 or abs(native-stats_all['original'])>1e-12:
        raise ArithmeticError('Changed native positive-order quadrature control')
    if stats[3]>2e-14 or not stats[0]>0:
        raise ArithmeticError('Reference-volume identity or ordered orientation failed')
    old,new = (table[...,index].sum(axis=-1,dtype=np.longdouble) for index in (1,2))
    diff = new-old
    ids = np.argsort(np.abs(diff),axis=None)[-8:][::-1]
    by_depth=[dict(material_depth_lower=float(body.t[i]),material_depth_upper=float(body.t[i+1]),
                   **summarize(table[i])) for i in range(len(body.t)-1)]
    by_radius=[dict(material_radial_chart_lower=float(body.x[i]),material_radial_chart_upper=float(body.x[i+1]),
                    **summarize(table[:,i])) for i in range(len(body.x)-1)]
    cell_path=out/'element_energies.npy'
    np.save(cell_path,table,allow_pickle=False)
    rows = dict(checkpoint=CHECKPOINT_SHA,label_sha256=digest(labels),factor=factor,gauss_points_unchanged=stats_all['count'],
        native_reference_tolerance=2e-14,min_sampled_Jinv=stats[0],max_sampled_Jinv=stats[1],
        max_roundtrip_abs_point_energy=stats[2],sum_absolute_roundtrip=stats[3],
        total=stats_all,signed_relative_deficit=-stats_all['difference']/stats_all['original'],
        by_correspondence={key:summarize(table[...,i,:]) for i,key in enumerate(GROUPS)},
        same_side=summarize(table[...,:2,:]),by_material_depth=by_depth,by_material_radius=by_radius,
        largest_changed_prisms=[dict(index=list(map(int,ijk)),original=float(old[ijk]),projected=float(new[ijk]),delta=float(diff[ijk]))
                                for ijk in (np.unravel_index(i,diff.shape) for i in ids)],
        sum_absolute_element_differences=float(np.sum(abs(diff),dtype=np.longdouble)),element_sha256=digest(cell_path),
        element_axes=['material_depth','material_radius','material_polar','prism_triangle','correspondence','column'],
        correspondence=GROUPS,columns=COLUMNS,
        decomposition='Exact labels and metric with substituted spatial inverse chart gradient first, then replace label values and reference metric at the same physical point.')
    (out/'result.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'factor':factor,'total':stats_all,'same_side':rows['same_side'],'relative_deficit':rows['signed_relative_deficit']},allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('labels',type=Path)
    p.add_argument('factor',type=int,choices=(2,4));p.add_argument('output',type=Path)
    args=p.parse_args();calculate(args.root.resolve(),args.labels.resolve(),args.factor,args.output.resolve())
