"""Trace the physical quadrature, reference metric, and cut branches of ranked R4 cells."""
import argparse
import ctypes as ct
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_label_gauge import depth_grid
from mouth_r4_aligned_derivatives import load
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_smooth_label_derivatives import CrystalBanks

LD=np.longdouble
NAMES=['t_index','a_index','segment','x_index','u','v','a','u_left','u_right','r','w','theta',
       'chi','Jinv','W','weighted_energy','core','physical_spline','reference_core_spline','polar_shape',
       'reg_eta','reg_t','upper_triangle','xi_x','xi_t','xi_a','left_x','left_t','left_a','right_x','right_t','right_a',
       'left_nx','left_nt','left_na','right_nx','right_nt','right_na',
       *[f'Finv_{i}{j}' for i in range(3) for j in range(3)],'azimuthal_Finv','chi_times_unit','dxix_du',
       'physical_x','physical_t','xi_minus_weld','native_x','native_t','native_a']


def digest(path):
    with path.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def compile_trace(directory):
    src=Path(__file__).with_suffix('.hpp').resolve()
    unit=directory/'trace.cpp'
    unit.write_text('#include "audit.cpp"\n#include '+json.dumps(str(src))+'\n')
    library=directory/'trace.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-fPIC','-shared',str(unit),'-o',str(library)],check=True)
    f=ct.CDLL(str(library)).depth_cell_trace
    f.argtypes=[ct.POINTER(HPData),P,P,HPGrid,P,ct.c_void_p,P,ct.c_int,P]
    f.restype=ct.c_int
    return f


def sigmas(row):
    mm=np.asarray(row[:,38:47],float).reshape(-1,3,3)
    s=np.linalg.svd(mm,compute_uv=False)
    if not np.allclose(np.abs(np.linalg.det(mm)*row[:,47]),row[:,13],atol=1e-13,rtol=3e-13):
        raise ArithmeticError('Four-dimensional determinant and inverse-map singular frame disagree')
    return np.minimum(s[:,-1],abs(row[:,47]))


def witness(rows,ind):
    return dict(ray=rows[ind,:2].astype(int).tolist(),segment=int(rows[ind,2]),radial_node=int(rows[ind,3]),
                fraction=rows[ind,4:7].astype(float).tolist(),r_w_theta=rows[ind,9:12].astype(float).tolist(),
                order=float(rows[ind,12]))


def inspection(row):
    s=sigmas(row)
    k=np.argmin(row[:,13]);j=np.argmin(s)
    return dict(samples=len(s),min_Jinv=float(row[k,13]),min_sigma_inv=float(s[j]),
                Jinv_witness=witness(row,k),sigma_witness=witness(row,j))


def topology(rows,axes,pressure):
    def cut_name(row,at):
        f=row[7+at]
        if f==0 or f==1:return str(int(f))
        xi=(1-f)*row[26:29]+f*row[29:32]
        n=row[32+3*at:35+3*at]
        if np.count_nonzero(n)>1:
            return f'reg_diagonal(e={int(row[20])},t={int(row[21])})'
        component=np.flatnonzero(n)[0]
        tables=[('x',axes['x']),('reg_e',np.exp(-pressure*axes['reg_e']))] if component==0 else (
                [('t',axes['t']),('reg_t',axes['reg_t'])] if component==1 else [('a',axes['a'])])
        choices=[(abs(y-xi[component]),name,i) for name,values in tables for i,y in enumerate(values)]
        mind=min(k[0] for k in choices)
        return '+'.join(f'{name}[{i}]' for er,name,i in choices if er<=mind+8e-17)
    rays={}
    for t in range(3):
        for a in range(3):
            ray=rows[(rows[:,0]==t)&(rows[:,1]==a)]
            first=ray[::3]
            cuts=[cut_name(r,0) for r in first]+['1']
            # Compare which actual metric branches each positive-volume sample executes.
            branches=[tuple(map(int,r[16:23])) for r in ray]
            def compact(xs):
                return [v for i,v in enumerate(xs) if i==0 or xs[i-1]!=v]
            rays[t,a] = dict(cuts=cuts,regions=compact(branches),sample_regions=branches,
                energy=float(np.sum(ray[:,15],dtype=LD)),range_t=ray[0,[27,30]].astype(float).tolist())
    return rays


def describe(cell,axes,order,groups,shadow):
    i,j,k=cell
    x0,x1=axes['field_x'][i:i+2]
    corner_chi=order[i:i+2,j:j+2,k:k+2]
    vertices=shadow[i:i+2,j:j+2,k:k+2,0]
    weld=axes['x'][(len(axes['x'])-1)*2//3]
    boundaries=[('axis','outer_wall'),('entrance','terminal_disk'),('polar_axis','equator')]
    faces=[n[s] for k,n,ax in zip(cell,boundaries,('x','t','a'))
             for s in range(2) if k==(len(axes['field_'+ax])-2 if s else 0)]
    return dict(cell=list(map(int,cell)),reference_weld_label=float(weld),
                spatial_grid_weld=bool(x0<=weld<=x1),material_weld_contact=bool(vertices.min()<=weld<=vertices.max()),
                material_label_radial_range=[float(vertices.min()),float(vertices.max())],domain_boundary=faces,
                chi_vertex_range=[float(corner_chi.min()),float(corner_chi.max())],
                states={name:inspection(matrix) for name,matrix in groups.items() if name in ('zero','minus_wide','plus_wide')})


def run(study,base,output):
    archive=output/'samples.npz'
    if (output/'dossier.json').exists():
        raise FileExistsError(output/'dossier.json')
    output.mkdir(parents=True,exist_ok=True)
    primary=study/'primary'
    rank_path=primary/'cells.json'
    rank=json.loads(rank_path.read_bytes())['top']
    cells=np.array([r['cell'] for r in rank],dtype=np.int32)
    if len(cells)!=10:raise ValueError('Expected ten cells ranked on the accepted primary map')
    with np.load(primary/'cells.npz',allow_pickle=False) as f:
        target=f['energy']; dims=tuple(f['dims'])
    driver,t,solid,order=load(base,'primary',output)
    if dims!=tuple(np.array(order.shape)-1):raise ValueError('Inconsistent cell grid')
    length=LD(t.shape.length)
    steps=LD(2)**-18*np.array([0,-1,1,-.5,.5,-.25,.25],LD)
    names=['zero','minus_wide','plus_wide','minus_half','plus_half','minus_narrow','plus_narrow']
    if archive.exists():
        with np.load(archive,allow_pickle=False) as saved:
            if not np.array_equal(saved['cells'],cells) or not np.array_equal(saved['steps'],steps) or list(saved['columns'])!=NAMES:
                raise ValueError('Saved quadrature archive does not match selected cells and steps')
            raw={k:saved[k] for k in saved.files if k.startswith('cell') and k!='cells'}
    else:
        bank=CrystalBanks(t.b,order,output)
        sample=compile_trace(output)
        original_t=bank.buffers['field_t'].copy()
        base_t=depth_grid(length,driver.reference_depth)
        raw={}
        for b,(name,h) in enumerate(zip(names,steps)):
            l=length*np.exp(h)
            dt=depth_grid(l,driver.reference_depth)-base_t
            bank.buffers['field_t'][:]=original_t+dt
            raw_l=np.array(l,LD)
            ct.memmove(ct.addressof(bank.data)+HPData.length.offset,raw_l.ctypes.data,raw_l.nbytes)
            lab=driver.labels.astype(LD);lab[...,1]+=dt[None,:,None]
            for rank,cell in enumerate(cells):
                row=np.empty((36000,len(NAMES)),LD)
                coord=np.zeros(3,LD)
                n=sample(ct.byref(bank.data),lab.ctypes.data_as(P),bank.order.ctypes.data_as(P),bank.rule,
                         bank.weights.ctypes.data_as(P),cell.ctypes.data,row.ctypes.data_as(P),len(row),coord.ctypes.data_as(P))
                if n<=0:raise RuntimeError(('Sample collection',cell,n))
                raw[f'cell{rank+1}_{name}']=row[:n].copy()
        with archive.open('xb') as stream:
            np.savez_compressed(stream,**raw,columns=NAMES,cells=cells,steps=steps)
    from mouth_cp10_label_condition import physical_frame
    results=[]; expected_deltas=[]
    for rank,cell in enumerate(cells):
        flat=np.ravel_multi_index(cell,dims)
        maps={}; ray_top={}; state_error=[]
        axes=t.b.buffers
        xyz=np.array([np.mean(axes['field_'+name][ind:ind+2]) for name,ind in zip('xta',cell)])
        frame,_=physical_frame(t.b,xyz[None])
        center=np.array([frame[0,2,2],float(length)*xyz[1],xyz[2]])
        for b,(name,h) in enumerate(zip(names,steps)):
            matrix=maps[name]=raw[f'cell{rank+1}_{name}']
            observed=matrix[:,15].sum(dtype=LD)
            err=abs(observed-target[b,flat])
            if err>1e-18*max(1.,abs(observed)):raise ArithmeticError((list(cell),name,str(observed),str(target[b,flat]),str(err)))
            state_error.append(float(err))
            ray_top[name]=topology(matrix,axes,LD(t.b.data.pressure))
        row=describe(cell,axes,order,maps,driver.labels)
        row.update(rank=rank+1,center_r_w_theta=center.astype(float).tolist(),closure_max_abs=max(state_error))
        original=ray_top['zero']; differences=[]; ray_ladder=[]
        for ray,main in original.items():
            for name in names[1:]:
                test=ray_top[name][ray]
                changed=[x for x in ('cuts','regions','sample_regions') if main[x]!=test[x]]
                if changed:
                    differences.append(dict(ray=list(ray),bank=name,changed=changed,baseline_cuts=main['cuts'],
                        new_cuts=test['cuts'],baseline_regions=main['regions'],new_regions=test['regions'],
                        left_right_reference_depth=test['range_t'],baseline_left_right_depth=main['range_t']))
            f=[]
            for k in range(3):
                minus=maps[names[1+2*k]];plus=maps[names[2+2*k]];zero=maps['zero']
                mm=minus[(minus[:,0]==ray[0])&(minus[:,1]==ray[1]),15].sum(dtype=LD)
                pp=plus[(plus[:,0]==ray[0])&(plus[:,1]==ray[1]),15].sum(dtype=LD)
                zz=zero[(zero[:,0]==ray[0])&(zero[:,1]==ray[1]),15].sum(dtype=LD)
                h=steps[2+2*k]
                f.append(dict(h=float(h),force=float((mm-pp)/(2*h)),minus=float((zz-mm)/(-h)),plus=float((zz-pp)/h)))
            ray_ladder.append(dict(ray=list(ray),levels=f,wide_minus_narrow_force=f[0]['force']-f[-1]['force']))
        row.update(ray_changes=differences,ray_ladder=sorted(ray_ladder,key=lambda v:abs(v['wide_minus_narrow_force']),reverse=True))
        results.append(row);expected_deltas.extend(state_error)
    result=dict(primary_state=digest(base/'primary/rematched.npz'),ranked_cells_sha256=digest(rank_path),
                cell_work_sha256=digest(primary/'cells.npz'),four_dimensional_orthonormal_frame=True,
                source_sha256={f.name:digest(f) for f in [Path(__file__),Path(__file__).with_suffix('.hpp')]},
                quadrature='Identical radial cuts and 3x3x3 integration, reported minima are over all physical integration points in each cell. Coordinates are r,w,theta at its mapped center.',
                singular_value='Smallest singular value of reference gradient in the native 3D orthonormal meridional frame and the azimuthal stretch, selecting the smaller.',
                domain_faces='Axes follow the actual sampled theta endpoints.',
                material_branch_tuple=['core','physical_spline_segment','reference_core_spline_segment','polar_window','registry_eta','registry_depth','upper_triangle'],
                physical_energy_max_absolute_difference=max(expected_deltas),raw_sha256=digest(output/'samples.npz'),ranked=results)
    (output/'dossier.json').write_text(json.dumps(result,allow_nan=False,indent=2)+'\n')
    print(json.dumps(dict(ranked=cells.tolist(),energy_max_error=max(expected_deltas),
                          ray_changes=[len(r['ray_changes']) for r in results],sha256=digest(output/'dossier.json'))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ('study','matched_study','output'):parser.add_argument(arg,type=Path)
    a=parser.parse_args()
    run(a.study.resolve(),a.matched_study.resolve(),a.output.resolve())
