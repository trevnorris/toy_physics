"""Locate the first appearance of the axial registry-plane cut on saved radial quadrature rays."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_cp10_label_condition import physical_frame
from mouth_cp10_label_gauge import depth_grid
from mouth_r4_aligned_derivatives import load
from mouth_r4_depth_cell_geometry import digest

LD=np.longdouble


def frame_at_face(body, xyz, col, dep, upper):
    """Exact two affine-registry triangles at a common physical label (radial coordinate is logarithmic)."""
    x,t,angle=xyz
    eta=-np.log(x)/LD(body.data.pressure)
    egrid,zgrid=(body.buffers['reg_'+k].astype(LD) for k in ('e','t'))
    points=body.buffers['labels'].reshape(len(zgrid),len(egrid),2).astype(LD)
    de=egrid[col+1]-egrid[col]; dt=zgrid[dep+1]-zgrid[dep]
    u,v=(eta-egrid[col])/de,(t-zgrid[dep])/dt
    if upper:
        p,r,z=points[dep+1,col+1],points[dep+1,col],points[dep,col+1]
        u,v,sign=1-u,1-v,-1
    else:
        p,r,z=points[dep,col],points[dep,col+1],points[dep+1,col]
        sign=1
    natural=p+u*(r-p)+v*(z-p)
    natural[0]+=4
    mat=np.zeros((3,3),LD)
    mat[:2,0]=sign*(r-p)/(-LD(body.data.pressure)*x*de)
    mat[:2,1]=sign*(z-p)/dt;mat[2,2]=natural[0]
    return natural,mat


def event(driver,trial,order,cell,ray,rows):
    i,j,k=cell
    a,b=ray
    indices=np.asarray(list(product((0,1),repeat=3)))
    nodes=tuple((indices+cell).T)
    val=rows['zero'][(rows['zero'][:,0]==a)&(rows['zero'][:,1]==b)][0]
    wide=rows['minus_wide'][(rows['minus_wide'][:,0]==a)&(rows['minus_wide'][:,1]==b)]
    body=trial.b; axes=body.buffers
    depth_plane=LD(axes['reg_t'][1]); length=LD(trial.shape.length)
    if depth_plane!=LD(axes['t'][2]):
        raise ValueError('Registry and native subdivision do not share the axial plane')
    crossed=[side for side in (0,1) if (val[27+3*side]-depth_plane)*(wide[0,27+3*side]-depth_plane)<0]
    if len(crossed)!=1:raise ValueError('Expected exactly one shared radial face to cross')
    side=crossed[0]; local=np.array([side,val[5],val[6]],LD)
    mult=np.where(indices,local,1-local)
    basis=np.prod(mult,axis=1)
    bgrad=np.column_stack([(2*indices[:,d]-1)*np.prod(np.delete(mult,d,axis=1),axis=1) for d in range(3)])
    lab=driver.labels[nodes].astype(LD)
    base_t=depth_grid(length,driver.reference_depth)
    def moved(h):
        delta=depth_grid(length*np.exp(LD(h)),driver.reference_depth)-base_t
        n=lab.copy();n[:,1]+=delta[j+indices[:,1]]
        return n,delta
    def f(h):return basis @ moved(h)[0][:,1]-depth_plane
    l,r=-(LD(2)**-18),-(LD(2)**-19)
    if f(l)*f(r)>=0:raise ArithmeticError('Archived widths do not bracket the same crossing')
    for iteration in range(100):
        mid=(l+r)/2
        if mid==l or mid==r:break
        if (f(l)>0)==(f(mid)>0):l=mid
        else:r=mid
    threshold=(l+r)/2
    material,dt=moved(threshold)
    xi=basis @ material; xi[1]=depth_plane
    physical_axes=[axes['field_x'],axes['field_t']+dt,axes['field_a']]
    width=np.array([ax[ii+1]-ax[ii] for ii,ax in zip(cell,physical_axes)],LD)
    pos=np.array([ax[ii]+w*q for ax,ii,w,q in zip(physical_axes,cell,width,local)],LD)
    physical,_=physical_frame(body,np.asarray(pos[None],float))
    physical=physical[0].astype(float);physical[1,1]=float(length*np.exp(threshold))
    B=physical*np.array(width,float);Q=material.T@bgrad
    chi=float(basis@order[nodes]); sides=[]; natural=[]; maps=[]
    for dep,upper in [(0,True),(1,False)]:
        coords,H=frame_at_face(body,xi,6,dep,upper);natural.append(coords);maps.append(H)
        A=H.astype(float)@np.asarray(Q,float)@np.linalg.inv(B)
        hoop=float(coords[0]*np.sin(xi[2])/(physical[2,2]*np.sin(pos[2])))
        ji=abs(np.linalg.det(A))*hoop
        J=1/ji; Finv=np.linalg.inv(A)
        W=.5*(np.sum(Finv**2)+hoop**-2-4)-np.log(J)+.25*(J-1)**2
        measure=4*np.pi*np.linalg.det(B)*physical[2,2]*float(np.sin(pos[2]))
        sides.append(dict(registry=[6,dep,int(upper)],reference_gradient=H.astype(float).tolist(),
                          Jinv=ji,W=float(W),chi_Jinv_W=chi*ji*float(W),
                          weighted_cell_integrand=measure*chi*ji*float(W)))
    gap=max(abs(natural[1]-natural[0]))
    if gap>1e-14:raise ArithmeticError('Shared material plane failed position continuity')
    return dict(cell=list(cell),ray=list(ray),negative_threshold_h=float(threshold),bracket_h=[float(l),float(r)],
                radial_endpoint=side,local=local.astype(float).tolist(),xi=xi.astype(float).tolist(),chi=chi,
                xi_t_minus_plane_baseline=float(f(0)),xi_t_minus_plane_minus_wide=float(f(-LD(2)**-18)),
                xi_t_minus_plane_minus_half=float(f(-LD(2)**-19)),
                entered_segment_width_minus_wide=float(np.max(wide[wide[:,21]==1,8]-wide[wide[:,21]==1,7])),
                physical_r_w_theta=[float(physical[2,2]),float(pos[1]*length*np.exp(threshold)),float(pos[2])],
                natural_position_gap=float(gap),reference_gradient_jump=(maps[1]-maps[0]).astype(float).tolist(),
                sided_point_values=sides,chi_Jinv_W_jump=sides[1]['chi_Jinv_W']-sides[0]['chi_Jinv_W'])


def run(study,base,output):
    source=study/'dossier'; report=json.loads((source/'dossier.json').read_bytes())
    if digest(source/'samples.npz')!=report['raw_sha256']:raise ValueError('Saved trace arrays differ')
    driver,trial,_,order=load(base,'primary',study)
    if output.exists():raise FileExistsError(output)
    top=[]; events=[]
    with np.load(source/'samples.npz',allow_pickle=False) as z:
        for row in report['ranked']:
            rank=row['rank'];cell=row['cell']
            label=dict(rank=rank,cell=cell,center_r_w_theta=row['center_r_w_theta'],material_weld_contact=row['material_weld_contact'],
                        spatial_grid_weld=row['spatial_grid_weld'],domain_boundary=row['domain_boundary'],chi_vertex_range=row['chi_vertex_range'],
                        sampled=row['states'],ray_changes=row['ray_changes'],energy_max_absolute_difference=row['closure_max_abs'])
            if row['ray_changes']:
                if len(row['ray_changes'])!=1 or row['ray_changes'][0]['bank']!='minus_wide':
                    raise ArithmeticError('Unexpected branch pattern for the threshold calculation')
                ray=row['ray_changes'][0]['ray']
                mats={n:z[f'cell{rank}_{n}'] for n in ('zero','minus_wide')}
                ev=event(driver,trial,order,cell,ray,mats);events.append(ev)
                rr=row['ray_ladder'][0]
                if rr['ray']!=ray:raise ArithmeticError('Largest residual ray does not match the branch change')
                ev.update(rank=rank,ray_force_ladder=rr['levels'],wide_minus_narrow_ray_force=rr['wide_minus_narrow_force'])
            top.append(label)
    result=dict(state_sha256=report['primary_state'],source_sha256={f.name:digest(f) for f in
             [Path(__file__),Path(__file__).with_name('mouth_r4_depth_cell_geometry.py'),Path(__file__).with_name('mouth_r4_depth_cell_geometry.hpp')]},
            trace_sha256=report['raw_sha256'],raw_work_sha256=report['cell_work_sha256'],
            footprint=report['quadrature'],coordinates='Cell centers in physical r,w,theta. Jinv and smallest singular value use both material metrics, including azimuthal stretch.',
            capture_guard=dict(invocation=json.loads((study/'dossier_guard/invocation.json').read_bytes())['command'],
                     exit_code=json.loads((study/'dossier_guard/outcome.json').read_bytes())['exitCode'],
                     note='Wrote every sample before emitting a NumPy int64 in the first summary. Readback used samples.npz without replay.'),
            readback_guard=dict(command=json.loads((study/'dossier_read_guard/invocation.json').read_bytes())['command'],
                     outcome=json.loads((study/'dossier_read_guard/outcome.json').read_bytes())['exitCode']),
            branch_tuple=report['material_branch_tuple'],cell_energy_max_abs_error=report['physical_energy_max_absolute_difference'],
            top=top,events=events)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(digest=digest(output),crossings=[(e['cell'],e['negative_threshold_h'],e['chi_Jinv_W_jump']) for e in events]),
                     allow_nan=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('study','base','output'):p.add_argument(key,type=Path)
    args=p.parse_args();run(args.study.resolve(),args.base.resolve(),args.output.resolve())
