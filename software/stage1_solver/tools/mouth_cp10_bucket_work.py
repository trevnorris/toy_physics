"""Energy and order-shift work by physical cell, original depth band and weld side."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_nested_energy import acquire_labels
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, acquired, digest


COLUMNS = ('E_native','E_spatial','F_native','F_spatial','delta_E','delta_F',
           'grad_only_delta_E','grad_only_delta_F','abs_delta_E','abs_delta_F')
WITNESS = ('weighted_point_delta_F','r','w','theta','chi','minus_dchi','J_L','J_E',
           'W_L','W_E','dV','material_depth','material_side_index','material_x','material_a','chart_x','chart_t')


def calculate(root, degree, output):
    output.mkdir(parents=True, exist_ok=False)
    build, shape, q, fields, _, _ = acquired(root/'r4_melted_survival01', root/'r4_exact_support01/patch')
    from mouth_free_shape_solid import Data, Grid, Pointer, Results, pointer
    solid = build(shape).b
    if np.any(solid.buffers['base']):
        raise ValueError('Expected absolute order')
    src = Path(__file__).with_suffix('.cpp')
    library = output/'bucket_work.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(src),'-o',str(library)],check=True)
    work = ct.CDLL(str(library)).bucket_work
    byte = ct.POINTER(ct.c_ubyte)
    extended = ct.POINTER(ct.c_longdouble)
    work.argtypes = [ct.POINTER(Data),ct.POINTER(Results),ct.POINTER(Grid),ct.POINTER(Pointer),
                     Pointer,byte]+[extended]*3+[ct.POINTER(extended),Pointer,ct.POINTER(ct.c_longlong)]
    work.restype = ct.c_int
    paths = [root/'r4_cp10_labels01/labels.npz',
             *(root/f'r4_nested_labels01/{factor}/labels.npz' for factor in (2,4))]
    tables = [acquire_labels(path,solid,level) for path,level in zip(paths,(1,2,4))]
    axes = (Grid*9)(*(Grid(len(a),pointer(a)) for grid,_ in tables for a in grid))
    labels = (Pointer*3)(*(pointer(xi) for _,xi in tables))
    t = solid.buffers['field_t']
    length = solid.data.length
    tdot = np.where(t*length<length-.12,(t*length/(length-.12))**25,1.)-t
    chi = fields[2][solid.ids]
    last = np.where(np.any(chi>0,axis=(0,1)))[0][-1]
    corner_t = solid.t[:,None,None]+solid.unpack(q)[...,1]
    min_t = np.minimum.reduce([corner_t[j:j+len(solid.t)-1,i:i+len(solid.x)-1,k:k+len(solid.a)-1]
                               for j in (0,1) for i in (0,1) for k in (0,1)])
    active = np.ascontiguousarray(min_t<=t[min(last+1,len(t)-1)],dtype=np.uint8)
    locations, weights = np.polynomial.legendre.leggauss(degree)
    locations, weights = (locations+1.)/2., weights/2.
    solid.data.gauss = Grid(degree,pointer(locations))
    solid.data.weights = pointer(weights)
    total = np.zeros((3,len(COLUMNS)),np.longdouble)
    shape_native = tuple(len(a)-1 for a in tables[0][0])
    shape_fine = tuple(len(a)-1 for a in tables[2][0])
    spatial = np.zeros(shape_native+total.shape,np.longdouble)
    depth_side = np.zeros((len(solid.t)-1,2)+total.shape,np.longdouble)
    per_grid = [np.zeros(tuple(len(v)-1 for v in axes)+(2,),np.longdouble) for axes,_ in tables]
    fine = per_grid[2][...,1]
    witness = np.zeros(shape_fine+(len(WITNESS),))
    count = ct.c_longlong()
    cell_pointers = (extended*3)(*(a.ctypes.data_as(extended) for a in per_grid))
    sums = [a.ctypes.data_as(extended) for a in (total,spatial,depth_side)]
    solid.work = lambda data,result: work(data,result,axes,labels,pointer(tdot),active.ctypes.data_as(byte),
                                          *sums,cell_pointers,pointer(witness),ct.byref(count))
    solid.evaluate(q,fields[2],load=False,gradient=False)

    measurements = src.parents[3]/'research/pde_ledger_v3/_measurements'
    archive = json.loads((measurements/'mouth_r4_matched_knot.json').read_bytes())
    if archive['checkpoint_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Checkpoint provenance mismatch')
    checks = {}
    if degree in (6,8):
        item = archive['order_work'][str(degree)]
        if item['label_sha256']!=[digest(p) for p in paths] or count.value!=item['ordered_points']:
            raise ValueError('Incompatible ordered quadrature')
        baseline = np.array(item['sums'])
        checks['max_abs_archive_order_work'] = float(max(abs(total[:,2]-baseline[0,1]).max(),
                                                        abs(total[:,3]-baseline[1:,1]).max()))
        checks['max_abs_archive_signed_order_difference'] = float(np.max(abs(
            total[:,5]-np.array(item['projected_minus_native'])[:,1])))
        checks['max_abs_archive_pointwise_absolute_order_difference'] = float(np.max(abs(total[:,9]-baseline[1:,2])))
    else:
        old = json.loads((measurements/'mouth_r4_spatial_refinement.json').read_bytes())
        if count.value!=5291652 or old['checkpoint']!=CHECKPOINT_SHA:
            raise ValueError('Incompatible native degree-3 checkpoint')
        checks['max_abs_archive_energy_difference'] = float(np.max(abs(total[:,4]-[
            old['rows'][str(l)]['delta'] for l in (1,2,4)])))
        checks['max_abs_archive_native_energy'] = float(np.max(abs(total[:,0]-old['reference_energy'])))
        checks['max_abs_archive_order_difference'] = float(np.max(abs(total[:,5]-[
            old['checks'][f'fixed_chart_depth_{l}']['projected_minus_native']['order_knot_motion'] for l in (1,2,4)])))
    # View global cancellation against small groups before converting the arrays to float64.
    checks['max_spatial_closure'] = float(abs(spatial.sum(axis=(0,1,2))-total).max())
    checks['max_depth_side_closure'] = float(abs(depth_side.sum(axis=(0,1))-total).max())
    checks['per_grid_cell_closure'] = float(max(abs(g.sum(axis=(0,1,2))-total[i,4:6]).max()
                                                for i,g in enumerate(per_grid)))
    checks['difference_vs_difference_of_raw_sums'] = float(max(abs(total[:,1]-total[:,0]-total[:,4]).max(),
                                                           abs(total[:,3]-total[:,2]-total[:,5]).max()))
    portable = np.asarray(spatial,dtype=float)
    checks['max_saved_spatial_closure'] = float(abs(portable.sum(axis=(0,1,2),dtype=np.longdouble)-total).max())
    if max(checks.values())>6.e-14 or not np.isfinite(total).all():
        raise ArithmeticError(f'Conservation or source discrepancy: {checks}')
    ids = np.argsort(abs(fine),axis=None)[-10:][::-1]
    top = []
    for i in ids:
        key = np.unravel_index(i,shape_fine)
        point = dict(zip(WITNESS,map(float,witness[key])))
        point['material_side'] = ('core','collar')[int(point.pop('material_side_index'))]
        top.append(dict(fine_cell=list(map(int,key)),native_cell=[int(k//4) for k in key],
                        delta_F=float(fine[key]),**point))
    path = output/'cells.npz'
    extra = {f'cells_{level}x':np.asarray(g,dtype=float) for level,g in zip((1,2,4),per_grid)}
    np.savez_compressed(path,**extra,spatial=portable,depth_side=np.asarray(depth_side,dtype=float),
                        material_depth_edges=solid.t,total=np.asarray(total,dtype=float),columns=COLUMNS)
    record = dict(degree=degree,checkpoint_sha256=CHECKPOINT_SHA,label_sha256=[digest(p) for p in paths],
                  source_sha256=digest(src),cells_sha256=digest(path),ordered_points=count.value,
                  native_shape=shape_native,fine_shape=shape_fine,levels=[1,2,4],columns=COLUMNS,
                  totals=np.asarray(total,dtype=float).tolist(),checks=checks,top10=top,
                  witness_definition='At the Gauss point carrying largest absolute weighted signed work difference in each 4x cell; rank cells by absolute net signed order-work difference.')
    with (output/'readout.json').open('x') as dest:
        json.dump(record,dest,indent=2,allow_nan=False)
    print(json.dumps(dict(totals=record['totals'],checks=checks,ordered_points=count.value)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('degree',type=int,choices=(3,6,8))
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    calculate(args.root.resolve(),args.degree,args.output.resolve())
