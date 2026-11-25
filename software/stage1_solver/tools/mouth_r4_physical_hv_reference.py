"""Binary128 centered-force reference for the held R4 physical crystal tangent.

The reference promotes the *archived first variation*, including its moving
weld terms, independently of the binary64 tangent. Signed labels, material
gradients and their difference are all evaluated in binary128. The seven
directions and both probe widths are fixed without reading tangent results.
"""
import argparse
import ctypes as ct
from itertools import product
from pathlib import Path
import re
import subprocess

import numpy as np

from mouth_exact_checkpoint_control import digest, read
from mouth_free_shape_solid import Data, Grid
from mouth_r4_aligned_crystal import AlignedCrystal
from mouth_r4_collar166_state import load_held
from mouth_r4_cubic_held import approved_reference
from mouth_r4_eulerian_release import order_grid
from mouth_whole_r4_relax import atomic_json


STEP_SHA = '346fc15d7dc678db51b420bbfde7e07ba8cbffac68aabc9d5cdbd4dd4b149f4c'
VECTORS_SHA = 'b3a633b19af74acfac934330415fe165427f5a319e6cf4c6f3aaa698233d57d0'
BITS = np.array(list(product(range(2), repeat=3)))
HERE = Path(__file__).resolve().parent


def compile_oracle(pinned, out):
    def source(name):
        return (pinned/name).read_text()
    base = source('mouth_free_shape_solid.cpp').split('struct TriSlice {', 1)[0]
    chart = source('mouth_cp10_spatial_crystal.cpp').split('namespace {', 1)[1].split('struct Evaluated {', 1)[0]
    crystal = source('mouth_r4_eulerian_solid.cpp').split('namespace euler {', 1)[1]
    text = base+'\nnamespace {\n'+chart+'\n}\nnamespace euler {\n'+crystal
    loop = 'Q1 q;bool has_order=false;int ids[2][2][2];'
    signature = 'double *diagnostics,double *intact){'
    if text.count(loop)!=1 or text.count(signature)!=1:
        raise ValueError('Unrecognized archived physical gradient')
    text = text.replace(signature, 'double *diagnostics,double *intact,const int *selected){')
    text = text.replace(loop, 'if(!selected[(i*(nt-1)+j)*(na-1)+k])continue;\n'+loop)
    text = re.sub(r'\b(?:long\s+)?double\b','Real',text)
    for name, mapped in {'abs':'fabsq','isfinite':'finiteq','cos':'cosq','sin':'sinq',
                         'exp':'expq','log':'logq','sqrt':'sqrtq','signbit':'signbitq'}.items():
        text = re.sub(r'std::'+name+r'\b',mapped,text)
    text = text.replace('std::numeric_limits<Real>::infinity()','HUGE_VALQ')
    text = text.replace('n<52;n++','n<120;n++')
    text = '#include <quadmath.h>\n#include <stdexcept>\nusing Real=__float128;\n'+text
    src = out/'oracle.cpp'
    src.write_text(text+(HERE/'mouth_r4_physical_hv_quad_bridge.hpp').read_text())
    so = out/'oracle.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-fPIC','-shared',str(src),'-o',str(so),'-lquadmath'],check=True)
    kernel = ct.CDLL(str(so.resolve())).r4_quad_derivative
    kernel.argtypes = [ct.POINTER(Data)]+[ct.c_void_p]*3+[ct.c_double,Grid]+[ct.c_void_p]*3
    kernel.restype = ct.c_int
    return kernel


def frozen(root,folder):
    study=root/'r4_projector_repair01'
    step=study/'primary/step_00382.npz'
    bank=root/'r4_step383_audit01/cg_vectors.npz'
    if digest(step)!=STEP_SHA or digest(bank)!=VECTORS_SHA:
        raise ValueError('Frozen checkpoint or archived CG vector digest mismatch')
    driver,trial,solid,_=load_held(study,folder,1e-8,solve_fields=False)
    approved_reference(trial.b)
    with np.load(step,allow_pickle=False) as rec:
        if int(rec['radius'])!=4 or int(rec['winding'])!=1 or not np.array_equal(rec['x'],driver.x):
            raise ValueError('Held R4 preparation changed')
        labels=rec['labels']
        order=order_grid(trial,rec['field'])
    joint=AlignedCrystal(driver,trial,solid,labels,order)
    return labels,joint,bank


def weld_vectors(joint,labels):
    body=joint.trial.b
    axis=body.buffers['x']
    seam=axis[(len(axis)-1)*2//3]
    dims=joint.order.shape
    eligible=[]
    points=joint.solid.points
    fractions=np.array(list(product(points,repeat=2)))
    bta=np.prod(np.where(BITS[:4,1:][None,:,:],fractions[:,None,:],1-fractions[:,None,:]),axis=2)
    radial=labels[...,0].ravel()
    groups=[cm.groups.reshape(dims) for cm in joint.components]
    offsets=np.r_[0,joint.cuts]
    for index in np.argwhere(joint.opt.active):
        cells=index+BITS
        ids=np.ravel_multi_index(tuple(cells.T),dims)
        chi=joint.order.ravel()[ids]
        if chi.max()==0 or radial[ids].min()>=seam or radial[ids].max()<=seam:
            continue
        # Choose supports on eight entirely free, unaliased nodes. Each has a
        # physical load, and its perturbation never displaces a held label.
        if any(np.any(cm.held[cm.groups[ids]]) or
               np.any(cm.representative[cm.groups[ids]]!=ids)
               for cm in joint.components):
            continue
        left,right=bta@radial[ids[:4]],bta@radial[ids[4:]]
        crossed=((left<seam)&(right>seam))|((left>seam)&(right<seam))
        if np.any(crossed):
            eligible.append((tuple(int(i) for i in index),ids,chi,left,right,crossed))
    eligible.sort(key=lambda item:(item[0][1],item[0][2],item[0][0]))
    if len(eligible)<4:
        raise ValueError(f'Only {len(eligible)} distinct freely perturbable weld cells')
    for a,rank in enumerate(np.linspace(0,len(eligible)-1,4,dtype=int)):
        index,ids,chi,left,right,crossed=eligible[rank]
        packed=np.zeros(len(joint.active))
        for c,cm in enumerate(joint.components):
            ng=groups[c].ravel()[ids]
            indices=offsets[c]+(np.cumsum(~cm.held)-1)[ng]
            coeff=np.array((1., .25*(-1)**a, .4*(-1)**(a//2)))[c]
            profile=(1.+.25*(a+1)*BITS[:,0]+.2*BITS[:,1]-.15*BITS[:,2])
            packed[indices]=coeff*profile*joint.width[indices]
        u=(seam-left[crossed])/(right[crossed]-left[crossed])
        direct=joint.lift(packed)[...,0].ravel()
        dleft,dright=(bta@direct[ids[:4]])[crossed],(bta@direct[ids[4:]])[crossed]
        speed=-((1-u)*dleft+u*dright)/(right[crossed]-left[crossed])
        if np.min(abs(speed))<=0:
            raise ValueError('Target failed to move a witness weld ray')
        yield 'weld'+str(a+1),packed,dict(cell=index,chi_vertex_range=[float(chi.min()),float(chi.max())],
               eligible_cells=len(eligible),weld=float(seam),weld_cut_fractions=u.tolist(),
               signed_weld_cut_speed=speed.tolist())


def norm(v):
    return float(np.sqrt(np.sum(np.square(v,dtype=np.longdouble),dtype=np.longdouble)))


def reference(root,out):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'reference.json').exists():
        raise FileExistsError('The reference bank already exists')
    pinned=root/'r4_projector_repair01/pinned/software/stage1_solver/tools'
    labels,joint,bank=frozen(root,out)
    with np.load(bank,allow_pickle=False) as vectors:
        u,v=vectors['p30'],vectors['p35']
    directions=[('u',u,{}),('v',v,{}),('v_minus_u',v-u,{})]
    directions.extend(weld_vectors(joint,labels))
    kernel=compile_oracle(pinned,out)
    order=np.ascontiguousarray(joint.order)
    gp=joint.solid.rule
    records={}
    for name,packed,info in directions:
        lifted=joint.lift(np.where(joint.active,packed,0.))
        windows=np.lib.stride_tricks.sliding_window_view(np.any(lifted!=0,axis=-1),(2,2,2))
        selected=np.ascontiguousarray(joint.opt.active & windows.any(axis=(-3,-2,-1)),np.int32)
        width_scale=float(np.max(abs(np.where(joint.active,packed,0.))/joint.width))
        hs=[2.**-b/width_scale for b in (29,32)]
        products=[]
        for h in hs:
            derivative=np.zeros_like(labels)
            rc=kernel(ct.byref(joint.trial.b.data),labels.ctypes.data,order.ctypes.data,lifted.ctypes.data,
                     h,gp,joint.solid.weights.ctypes.data,selected.ctypes.data,derivative.ctypes.data)
            if rc or not np.isfinite(derivative).all():
                raise ArithmeticError(f'Binary128 centered oracle {name} failed ({rc})')
            products.append(derivative)
        full=norm(products[1]-products[0])/norm(products[1])
        freed=norm(joint.reduce(products[1]-products[0]))/norm(joint.reduce(products[1]))
        records[name]=dict(info,widths=hs,selected_cells=int(selected.sum()),full_width_spread=full,
                           free_width_spread=freed,full_nodal_norm=norm(products[1]),free_norm=norm(joint.reduce(products[1])))
        np.savez(out/(name+'.npz'),packed=packed,wide=products[0],narrow=products[1])
        records[name]['bank_sha256']=digest(out/(name+'.npz'))
        print(name,records[name],flush=True)
    atomic_json(out/'reference.json',dict(state_sha256=STEP_SHA,vectors_sha256=VECTORS_SHA,records=records,
                scalar_mantissa=113,oracle_source_sha256=digest(out/'oracle.cpp'),accepted_steps=0))


def compare(root,out):
    metadata=read(out/'reference.json')
    if metadata['state_sha256']!=STEP_SHA or metadata['vectors_sha256']!=VECTORS_SHA:
        raise ValueError('Oracle inputs do not match the requested gate')
    labels,joint,_=frozen(root,out)
    rows={}
    archived={}
    for name,record in metadata['records'].items():
        bank=out/(name+'.npz')
        if digest(bank)!=record['bank_sha256']:
            raise ValueError('Directional oracle bank changed')
        with np.load(bank,allow_pickle=False) as data:
            lifted=joint.lift(np.where(joint.active,data['packed'],0.))
            native=joint.solid.hessian_action(labels,joint.order,lifted)
            diff=native-data['narrow']
            rows[name]=dict(full_relative=norm(diff)/norm(data['narrow']),
                            free_relative=norm(joint.reduce(diff))/norm(joint.reduce(data['narrow'])),
                            max_absolute=float(np.max(abs(diff))),free_width_spread=record['free_width_spread'],
                            full_width_spread=record['full_width_spread'])
            if name in ('u','v','v_minus_u'):
                archived[name]=joint.reduce(data['narrow']).astype(np.longdouble)
    passed=all(max(row[k] for k in ('full_relative','free_relative','free_width_spread','full_width_spread'))<=1e-11 for row in rows.values())
    difference=archived['v_minus_u']
    reference_linearity=norm(difference-(archived['v']-archived['u']))/norm(difference)
    atomic_json(out/'comparison.json',dict(rows=rows,passed=passed,accepted_steps=0,state_sha256=STEP_SHA,
                reference_physical_subtractivity=reference_linearity))
    print(rows,flush=True)
    if not passed:
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['reference','compare'])
    parser.add_argument('scratch',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    (reference if args.stage=='reference' else compare)(args.scratch.resolve(),args.output.resolve())
