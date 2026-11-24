"""Measure one accepted R4 material prism under the producing constitutive law."""
import argparse
import ctypes as ct
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import mpmath as mp
import numpy as np


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def assemble_trial(checkpoint, state_root):
    from mouth_continuous_shape_relax import ContinuousRelease
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    run = checkpoint.parent
    invocation = read(run/'invocation.json')['command']
    original, _, proof_file, folder, level = invocation[-5:]
    original = Path(original)/level
    proof = read(Path(proof_file))['results'][level]
    prep = Path(read(original/'invocation.json')['command'][5])
    endpoint = original/'endpoint.json'
    if (digest(prep) != proof['preparation_sha256'] or digest(endpoint) != proof['endpoint_sha256']
            or run != Path(folder) or level != 'coarse'):
        raise ValueError('Original coarse preparation provenance mismatch')
    with np.load(checkpoint) as saved:
        if str(saved['state_sha']) != digest(prep):
            raise ValueError('Checkpoint preparation does not match')
        x, offsets, q, u = (saved[name].copy() for name in ('x','solid','q','field'))
    owner = SimpleNamespace(base=Front(prep,GRIDS['double'],1),args=SimpleNamespace(state=prep),
                            discretization=proof['solid_discretization'],reference_depth=read(endpoint)['depth'])
    trial = ContinuousRelease.build(owner,x,solve_fields=False)
    if not np.array_equal(offsets,trial.b.unpack(q)):
        raise ValueError('Packing changes saved material offsets')
    hashes = {str(p):digest(p) for p in (prep,endpoint,Path(proof_file),run/'invocation.json',state_root/'sources.json')}
    return trial,q,offsets,u,hashes


CORNERS = (((0,0),(1,0),(0,1)),((1,1),(0,1),(1,0)))


class Prism:
    def __init__(self, body, offsets, cell, triangle):
        self.j,self.i,self.k = cell
        self.triangle = triangle
        self.sign = 1-2*triangle
        self.grid = [v[z:z+2] for v,z in zip((body.x,body.t,body.a),(self.i,self.j,self.k))]
        self.indices = [(self.j+v,self.i+u) for u,v in CORNERS[triangle]]
        self.nodes = []
        for u,v in CORNERS[triangle]:
            pair = []
            for h in (0,1):
                off = offsets[self.j+v,self.i+u,self.k+h]
                coord = (self.grid[0][u],self.grid[1][v],self.grid[2][h])
                pair.append((coord,off))
            self.nodes.append(pair)

    def ref_point(self,bary,angle):
        u,v = np.array(bary) @ np.array(CORNERS[self.triangle])
        return [(1-f)*lo+f*hi for (lo,hi),f in zip(self.grid,(u,v,angle))]

    def chart(self,bary,y,cast,rounded):
        # Reconstruct nodal edge differences in exact rationals or high
        # precision. In the solver the unknown is an offset. Its chart value
        # is the mathematical grid + offset; the audit first rounds that sum.
        def number(x):
            return cast(float(x))
        nodes = [[[number(c+o) if rounded else number(c)+number(o)
                   for c,o in zip(coord,off)] for coord,off in pair] for pair in self.nodes]
        dx,dt,da = [number(hi)-number(lo) for lo,hi in self.grid]
        y=number(y)
        w=[number(b) for b in bary]
        # Interpret the input as barycentric coordinates with the first weight
        # fixed by partition of unity, as for the production triangle basis.
        w[0]=1-w[1]-w[2]
        M = []
        for c in range(3):
            row=[(1-y)*pair[0][c]+y*pair[1][c] for pair in nodes]
            M.append([self.sign*(row[1]-row[0])/dx,self.sign*(row[2]-row[0])/dt,
                      sum(ws*(pair[1][c]-pair[0][c])/da for ws,pair in zip(w,nodes))])
        return M


def determinant(m):
    a,b,c=m
    return a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0])


def q1(body,u,landed):
    axes=[body.buffers[name] for name in ('field_x','field_t','field_a')]
    cells=[min(max(np.searchsorted(z,v,side='right')-1,0),len(z)-2) for z,v in zip(axes,landed)]
    widths=[ax[i+1]-ax[i] for ax,i in zip(axes,cells)]
    y=[(v-ax[i])/h for ax,i,v,h in zip(axes,cells,landed,widths)]
    vals=np.zeros(3)
    grads=np.zeros((3,3))
    for inds in product((0,1),repeat=3):
        w=[v if b else 1-v for b,v in zip(inds,y)]
        a,b,c=np.add(cells,inds)
        d=u[:,body.ids[a,c,b]]
        vals+=d*np.prod(w)
        for z in range(3):
            grads[:,z]+=d*(2*inds[z]-1)*w[(z+1)%3]*w[(z+2)%3]/widths[z]
    return vals,grads,list(map(int,cells))


def native_callback(root,build):
    from mouth_free_shape_solid import Data,Pointer
    src=Path(__file__).resolve().with_suffix('.cpp')
    lib=build/'constitutive.so'
    build.mkdir(parents=True,exist_ok=True)
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC','-I',str(root),str(src),'-o',str(lib)],check=True)
    fn=ct.CDLL(str(lib)).fold_point
    fn.argtypes=[ct.POINTER(Data),ct.c_int,ct.c_int,ct.c_int,ct.c_int,ct.c_double,ct.c_double,ct.c_double,Pointer]
    fn.restype=None
    return fn


def point(prism,bary,y,trial,u,callback,*,tensors=False):
    from mouth_free_shape_solid import pointer
    b=trial.b
    x,t,a=prism.ref_point(bary,y)
    data=np.empty(47)
    callback(ct.byref(b.data),prism.i,prism.j,prism.k,prism.triangle,x,t,a,pointer(data))
    landed=data[:3]
    M,C,ref_inv,native_F=data[3:39].reshape(4,3,3)
    hoop,scale,native_det,native_J,native_W,S,code,native_chiW=data[39:]
    vals,grad,cells=q1(b,u,landed)
    v,_,chi=vals
    grad=np.linalg.solve(C.T,grad.T).T
    R=C[2,2]
    s=R*np.sin(landed[2])
    vortex,dvortex=trial.model._vortex(s)
    n=(vortex*v)**2
    grad_n=2*vortex*v*(vortex*grad[0]+v*dvortex*np.array([np.sin(landed[2]),0.,np.cos(landed[2])]))
    matrix=mp.matrix(prism.chart(bary,y,mp.mpf,False))
    F=mp.matrix(C.tolist())*matrix*mp.matrix(ref_inv.tolist())
    stretch=sorted([float(vv) for vv in mp.svd(F,compute_uv=False)]+[float(hoop)])
    det=mp.det(matrix)
    physical_j=mp.mpf(float(scale))*det
    signed_f=mp.det(F)*hoop
    log_term=-mp.log(physical_j) if physical_j>0 else None
    stretch_term=(sum(ff**2 for ff in F)+hoop**2-4)/2
    compression_term=(physical_j-1)**2/4
    W=float(stretch_term+log_term+compression_term) if log_term is not None else None
    domain_err=abs(float(native_chiW)-float(native_W)*chi) if code==0 else None
    row=dict(barycentric=list(map(float,bary)),angular_fraction=float(y),reference_x_t_theta=list(map(float,(x,t,a))),
                landed_x_t_theta=landed.tolist(),physical_r_z_theta=[R,b.sh.length*landed[1],landed[2]],
                cylindrical_radius=s,reference_radius=S,polar_hoop=hoop,chart_to_physical_J_scale=scale,
                det_chart=float(det),J=float(physical_j),signed_det_F=float(signed_f),
                det_ref_inverse=float(np.linalg.det(ref_inv)),principal_stretches=stretch,
                stretch_term=float(stretch_term),negative_log_J=None if log_term is None else float(log_term),
                compression_term=float(compression_term),W=W,chi=chi,chi_W=None if W is None else chi*W,
                n=n,order_gradient_r_z_polar=grad[2].tolist(),density_gradient_r_z_polar=grad_n.tolist(),
                field_cell_x_t_angle=cells,
                native=dict(det_chart=native_det,J=native_J,W=None if not np.isfinite(native_W) else native_W,
                            sample_status=int(code),sample_chi_W=None if code!=0 else native_chiW,
                            pointwise_chiW_difference=domain_err,
                            chart_gradient_max_abs_difference=float(np.max(abs(M-np.array(matrix.tolist(),float)))),
                            F_max_abs_difference=float(np.max(abs(native_F-np.array(F.tolist(),float))))))
    if tensors:
        U,sv,Vh=mp.svd(F)
        nx=[float(Vh[2,c]) for c in range(3)]
        orient=1 if nx[np.argmax(np.abs(nx))]>=0 else -1
        row['meridional_polar_blocks']=dict(chart_gradient=np.array(matrix.tolist(),float).tolist(),
                current_chart_to_physical=C.tolist(),inverse_reference=ref_inv.tolist(),
                physical_F=np.array(F.tolist(),float).tolist(),
                smallest_stretch=float(sv[2]),smallest_right_direction_r_z_polar=[orient*c for c in nx],
                smallest_left_direction_r_z_polar=[orient*float(U[c,2]) for c in range(3)])
    return row


def forces(trial,q,g,prism,vertex):
    from mouth_joined_angular_work import bump
    b=trial.b
    off=b.unpack(q)
    raw=np.zeros(b.chart_shape)
    raw[b.free]=-g
    jj,ii=prism.indices[vertex]
    def ring(j,i):
        X=off[j,i]+np.column_stack([np.full_like(b.a,b.x[i]),np.full_like(b.a,b.t[j]),b.a])
        xp,tp,a=X.T
        z=tp*b.sh.length
        width,slope=b.sh.radius(z)
        h,hz=bump(z,b.sh.length)
        leg=(3*np.cos(a)**2-1)/2
        e=b.epsilon
        R=width*xp+e*h*xp**3*leg
        rx=width+3*e*h*xp**2*leg
        rz=slope*xp+e*hz*xp**3*leg
        ra=-3*e*h*xp**3*np.cos(a)*np.sin(a)
        fr=raw[j,i,:,0]/rx
        fz=raw[j,i,:,1]/b.sh.length-rz*fr
        fa=(raw[j,i,:,2]-ra*fr)/R
        phys=np.column_stack((fr,fz,fa))
        phys[~b.free[j,i]]=0.
        return dict(node_t_x=[int(j),int(i)],landed_r_z_theta=np.column_stack((R,z,a)).tolist(),
                    physical_components_r_z_polar=phys.tolist(),summed_radial_axial=phys[:,:2].sum(axis=0).tolist())
    return dict(fold_vertex_ring=ring(jj,ii),terminal_free_ring=ring(len(b.t)-2,len(b.x)-2))


def main(args):
    from mouth_free_shape_solid import pointer
    from mouth_joint_endpoint_audit import orientation
    dec=read(args.decision)
    cp=Path(dec['checkpoint'])
    if digest(cp)!=dec['sha256'] or cp.name!=dec['latest']['file'] or dec['sha256']!=dec['latest']['sha256']:
        raise ValueError('Gate and accepted checkpoint differ')
    state_root=cp.parent.parent
    pin=state_root/'pin'
    for rel,sha in read(state_root/'sources.json').items():
        if digest(pin/rel)!=sha:
            raise ValueError('Input source pin changed: '+rel)
    trial,q,offsets,u,hashes=assemble_trial(cp,state_root)
    b=trial.b
    audit=orientation(trial.shape,trial.polar_shift,b,offsets)
    if audit!=dec['latest']['orientation']:
        raise ValueError('Accepted orientation audit did not reproduce')
    loc=audit['extremizers']['min']
    selected=Prism(b,offsets,loc['cell_t_x_angle'],loc['triangle'])
    vertex=loc['vertex']
    y0=loc['angular_fraction']
    start=np.eye(3)[vertex]
    middle=np.ones(3)/3
    fn=native_callback(pin/'software/stage1_solver/tools',args.build)
    source=Path(__file__).resolve()
    summary=dict(command=sys.argv,input_sha256={**hashes,**{str(p):digest(p) for p in (cp,source,source.with_suffix('.cpp'))}},
                 iteration=dec['iteration'],selected=loc,checkpoint_sha256=dec['sha256'],gate_det=audit['prism_chart_determinant_range'][0],
                 coordinate_order='chart x,t,theta; physical r,z,polar arc,azimuthal arc; W per reference four-volume',
                 precision='Chart based on nodal grid plus saved binary64 offsets at 75 decimal digits. C, reference inverse and hoop from pinned double-precision kernel. Direct native implementation reported alongside.',
                 discrete_orientation='Signed det(F) includes handedness of initial collar labels. J is the model ratio to the unsigned reference volume with signed current chart orientation.',
                 reference_metrics='At a label-grid boundary the production reference lookup selects its own one-sided tensor. The inward ray reports the limiting branch separately.',
                 paths='fraction 0 is selected corner, fraction 1 is the prism barycenter and angular midpoint',
                 precision_at_extremum={})
    for rounded,label in ((False,'exact_grid_plus_offset'),(True,'binary64_landed_nodes')):
        rat=determinant(selected.chart(start,y0,Fraction,rounded))
        summary['precision_at_extremum'][label]=dict(rational=str(rat),value=float(rat),sign=int(rat>0)-int(rat<0))
    energy,g,_,full=b.evaluate(q,u[2],load=False)
    residual=b.residual(q,g)
    summary['material']=dict(native_integrated_chi_W=energy,residual=residual,
                             force_vectors=forces(trial,q,g,selected,vertex))
    dc=np.ascontiguousarray(u[2,b.ids])
    b.data.offset,b.data.correction=pointer(offsets),pointer(dc)
    b.data.diagnose=0
    path=[]
    for step in (0.,1.e-12,1.e-10,1.e-8,1.e-6,1.e-4,.001,.01,.05,.1,.25,.5,1.):
        row=point(selected,(1-step)*start+step*middle,(1-step)*y0+.5*step,trial,u,fn,tensors=step in (0.,1.))
        path.append(dict(fraction_toward_prism_center=step,**row))
    summary['ray']=path
    parent=[]
    for tr in (0,1):
        prism=Prism(b,offsets,loc['cell_t_x_angle'],tr)
        for y in (0.,.5,1.):
            for bary in (*np.eye(3),middle):
                parent.append(dict(triangle=tr,**point(prism,bary,y,trial,u,fn)))
    summary['parent_prisms_sampled']=parent
    summary['sampled_range_convention']='Two triangles, each at all three meridional vertices and the barycenter on both angular faces and the angular midpoint; finite samples, no physical-J bound asserted.'
    summary['parent_sampled_ranges']={key:[min(vv),max(vv)] for key in ('det_chart','J','chi','n','W','chi_W')
                                      if (vv:=[e[key] for e in parent if e[key] is not None])}
    err=[s['native']['pointwise_chiW_difference'] for s in path+parent if s['native']['pointwise_chiW_difference'] is not None]
    summary['checks']=dict(native_integrated_matches=abs(energy-dec['latest']['field']['solid_energy'])<1.e-10,
                            archived_ring_reproduced=residual==dec['latest']['material'],
                            terminal_force_vector_sum=residual['last_free_ring_angular_sum']==summary['material']['force_vectors']['terminal_free_ring']['summed_radial_axial'],
                            exact_fraction_matches_high_precision=abs(path[0]['det_chart']-summary['precision_at_extremum']['exact_grid_plus_offset']['value'])<1.e-50,
                            native_pointwise_chiW_error_max=max(err,default=0.),
                            regular_chart_at_corner=bool(path[0]['cylindrical_radius']>0 and path[0]['polar_hoop']>0 and path[0]['chart_to_physical_J_scale']>0))
    args.output.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('CHECKPOINT',cp,digest(cp),'CORNER',loc,'DET',summary['precision_at_extremum'],'CHECKS',summary['checks'],flush=True)
    if not (all(v for k,v in summary['checks'].items() if k!='native_pointwise_chiW_error_max')
            and max(err,default=0.)<2.e-9):
        raise ArithmeticError('Readback checks failed; see output')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for label in ('decision','build','output'):
        p.add_argument(label,type=Path)
    args=p.parse_args()
    pin=Path(read(args.decision)['checkpoint']).parent.parent/'pin/software/stage1_solver/tools'
    sys.path.insert(0,str(pin))
    with mp.workdps(75):
        main(args)
