"""Read-only Q1 inverse-map probe on the actual physical carrier cell of an R4 crystal witness.

Curvature, four-dimensional hoop stretch and the exact side of an order-grid knot are
retained. Melted-node inverses are diagnostic continuations of the saved labels, not a
solution of a harmonic gauge. The necessary nodal range holds for any Q1 extension.
"""
import argparse
from bisect import bisect_right
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json

import numpy as np

from mouth_active_prism_orientation import _Support
from mouth_free_shape_solid import Solid
from mouth_r4_conforming_prisms import (_geometry, _children, seam_fraction, child_witness,
                                         STAR, TRIANGLES)
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def current(body, p):
    x,t,a=p
    L=body.data.length
    base,slope=body.sh.radius(t*L)
    y=(t-.06)/.56
    band=16*body.epsilon*y*y*(1-y)**2 if 0<y<1 else 0.
    derivative=32*body.epsilon*y*(1-y)*(1-2*y)/.56 if 0<y<1 else 0.
    leg=(3*np.cos(a)**2-1)/2
    r=x*base+x**3*band*leg
    gx=base+3*x*x*band*leg
    gt=x*L*slope+x**3*derivative*leg
    ga=-3*x**3*band*np.sin(a)*np.cos(a)
    sn,co=np.sin(a),np.cos(a)
    B=np.array([[sn*gx,sn*gt,sn*ga+co*r],[co*gx,co*gt,co*ga-sn*r],[0,L,0]])
    return np.array([r*sn,r*co,L*t]),B


def physical_reference(body,q):
    x,t,a=q
    e=-np.log(x)/body.data.pressure
    er,tr=body.buffers['reg_e'],body.buffers['reg_t']
    i=np.clip(np.searchsorted(er,e,side='right')-1,0,len(er)-2)
    j=np.clip(np.searchsorted(tr,t,side='right')-1,0,len(tr)-2)
    u=(e-er[i])/(er[i+1]-er[i]); v=(t-tr[j])/(tr[j+1]-tr[j])
    tri=int(u+v>1)
    uv=np.array([u,v]) if not tri else 1-np.array([u,v])
    weights=np.r_[1-sum(uv),uv]
    offsets=TRIANGLES[tri]
    labels=body.buffers['labels'].reshape(len(tr),len(er),2)
    X=weights@np.asarray([labels[j+dj,i+di] for dj,di in offsets])
    S=4+X[0]
    return np.array([S*np.sin(a),S*np.cos(a),X[1]])


def exact_coordinate(b,solid,surplus,w):
    key=(*w['cell_t_x_angle'],w['triangle'])
    ref,pts,_=_geometry(b,solid,key)
    a=Q(w['angular_fraction'])
    by,bz=map(Q,w['barycentric'][1:]); bary=[1-by-bz,by,bz]
    if key in STAR:
        f=seam_fraction(b); up=int(a>=f)
        ref,pts=_children(key,ref,pts,surplus,f)[up]
        a=(a-f)/(1-f) if up else a/f
    return np.asarray(bary)@((1-a)*ref[0]+a*ref[1]),np.asarray(bary)@((1-a)*pts[0]+a*pts[1])


def pieces(body,solid,surplus):
    # Build a broad-phase bound. A P1-triangle x P1-angle placement is a convex
    # combination of its six vertices irrespective of the sign of its Jacobian.
    Qchart=np.stack(np.meshgrid(body.t,body.x,body.a,indexing='ij'),axis=-1)[...,[1,0,2]]
    Y=Qchart+solid
    keys=[]; refs=[]; landed=[]
    f=float(seam_fraction(body))
    for j,i,k,tri in product(range(24),range(24),range(6),range(2)):
        key=j,i,k,tri
        rc=np.array([[Qchart[j+dj,i+di,k+l] for dj,di in TRIANGLES[tri]] for l in range(2)])
        pc=np.array([[Y[j+dj,i+di,k+l] for dj,di in TRIANGLES[tri]] for l in range(2)])
        if key in STAR:
            rm=(1-f)*rc[0]+f*rc[1]; pm=(1-f)*pc[0]+f*pc[1]
            pm[STAR[key]]+=surplus
            for h in (0,1):
                refs.append([rm,rc[1]] if h else [rc[0],rm]); landed.append([pm,pc[1]] if h else [pc[0],pm]); keys.append(key)
        else:
            refs.append(rc); landed.append(pc); keys.append(key)
    P=np.asarray(landed)
    return keys,np.asarray(refs),P,P.min((1,2)),P.max((1,2))


def roots(mesh,dest):
    keys,refs,placed,lo,hi=mesh
    hits=[]
    tol=1e-12
    for idx in np.flatnonzero(np.all((lo-tol<=dest)&(dest<=hi+tol),axis=1)):
        pc=placed[idx]; low=pc[0]; diff=pc[1]-low
        v,w=low[1]-low[0],low[2]-low[0]
        sv,sw=diff[1]-diff[0],diff[2]-diff[0]
        t,dt=dest-low[0],-diff[0]
        n0,n1,n2=np.cross(v,w),np.cross(sv,w)+np.cross(v,sw),np.cross(sv,sw)
        poly=np.array([dt@n2, t@n2+dt@n1,t@n1+dt@n0, t@n0])
        for z in np.roots(poly):
            if abs(z.imag)>1e-9 or z.real < -1e-9 or z.real > 1+1e-9:continue
            a=z.real
            pts=low+a*diff
            V=np.stack((pts[1]-pts[0],pts[2]-pts[0]),axis=1)
            eta=np.linalg.lstsq(V,dest-pts[0],rcond=None)[0]
            bc=np.r_[1-eta.sum(),eta]
            if min(bc)<-1e-9 or max(abs(bc@pts-dest))>tol:continue
            ref=bc@((1-a)*refs[idx,0]+a*refs[idx,1])
            if not any(np.max(abs(ref-r))<1e-9 for r in hits):
                hits.append(ref)
    return hits


def test(body,solid,surplus,order,w):
    truth=child_witness(body,solid,surplus,order,w)
    reference,p=exact_coordinate(body,solid,surplus,w)
    px=np.array(p,float)
    support=_Support(body,order)
    value=float(support.at(p))
    if value!=truth['chi'] or value<=0:
        raise ValueError('Locate a strictly ordered point')
    phys,B=current(body,px); invB=np.linalg.inv(B)
    orig=[list(map(Q,body.buffers['field_'+a])) for a in ('x','t','a')]
    cells=[bisect_right(v,coord)-1 for v,coord in zip(orig,p)]
    mesh=pieces(body,solid,surplus)
    normI=1/min(truth['principal_stretches'])
    ref_norm_bound = max(np.linalg.norm(physical_reference(body, np.array([float(a),float(t),.7])))
                          for a,t in product(body.x[body.cut_node:],body.t))
    report=[]
    for factor in (1,2,4,8,32,128):
        cp=[]
        for a,i,coord in zip(orig,cells,p):
            h=(a[i+1]-a[i])/factor
            n=min(factor-1,(coord-a[i])//h)
            cp.append((a[i]+n*h,a[i]+(n+1)*h))
        widths=np.array([float(hi-lo) for lo,hi in cp])
        # The largest physical inverse singular value is meridional/polar.
        # Q1 derivatives along each chart dimension are convex averages of
        # nodal edge differences divided by cell width. |I - F^-1| >= |F^-1|-1.
        if normI<=1/truth['hoop']+1e-10:
            raise ValueError('This bound needs a meridional or polar limiting stretch')
        scale=np.sqrt(3)*np.linalg.norm(1/widths)*np.linalg.norm(invB,2)
        req=max(0.,(normI-1)/scale)
        minimum_ranges=dict(max_nodal_span_needed_for_exact_stretch=req,
                            coefficient_bound_required_for_99pct_of_inverse_gradient=max(0.,(.99*normI-1)/scale))
        counts=[]; order_corners=[]; choices=[]; residuals=[]
        for node,bits in zip(product(*cp), product((0,1),repeat=3)):
            rc=roots(mesh,np.array(node,float))
            counts.append(len(rc));order_corners.append(float(support.at(node)))
            # Where material labels have vanished, measure the archived label
            # continuation if it is single-valued. Do not select a root by convenience.
            if len(rc)==1: choices.append(rc[0])
        data=dict(factor=factor,cell_index=cells,chart_widths=widths.tolist(),inverse_sheet_counts=counts,
                  ordered_vertex_count=sum(v>0 for v in order_corners),
                  order_vertex_values=order_corners,**minimum_ranges)
        if len(choices)==8 and all(r[0]>body.cut for r in choices):
            U=np.array([current(body,np.array(node,float))[0]-physical_reference(body,r)
                        for node,r in zip(product(*cp),choices)])
            active=np.array(order_corners)>0
            loc=np.array([float((r-a)/(b-a)) for r,(a,b) in zip(p,cp)])
            weights=np.array([[1-v,v] for v in loc])
            u_interp=np.zeros(3);grad=np.zeros((3,3))
            for u,bits in zip(U,product((0,1),repeat=3)):
                prod=np.array([weights[a,idx] for a,idx in enumerate(bits)])
                u_interp+=np.prod(prod)*u
                for a,side in enumerate(bits):
                    grad[:,a]+=(2*side-1)*u*np.prod(np.delete(prod,a))/widths[a]
            C=np.eye(3)-grad@invB
            circ=(phys[0]-u_interp[0])/phys[0]
            # The saved collar's registered meridional chart reverses orientation.
            # Compare absolute determinant and singular values here.
            jac=1/abs(np.linalg.det(C)*circ)
            svs=[*(1/np.linalg.svd(C,compute_uv=False)), 1/abs(circ)]
            strain=(sum(v*v for v in svs)-4)/2-np.log(jac)+(jac-1)**2/4
            err=np.linalg.norm((phys-u_interp)-physical_reference(body,np.array(reference,float)))
            data.update(Q1_saved_melt_label_continuation=dict(max_nodal_u_component_span=float(np.max(np.ptp(U,axis=0))),
                inverse_spectral_norm=float(max(1/np.array(svs))), min_stretch=min(svs), J=jac,W=strain,
                absolute_reference_coordinate_error=err, signed_cartesian_det_inverse=float(np.linalg.det(C)*circ)))
        report.append(data)
    return dict(truth=truth, spatial_point_float64=px.tolist(), point_order_at_float64=float(support.at(np.array([Q(v) for v in px]))),
                exact_order=value, smallest_forward_stretch=min(truth['principal_stretches']),
                target_inverse_norm=normI, physical_xyz=phys.tolist(),
                largest_norm_at_registered_collar_sample=ref_norm_bound,
                levels=report)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args()
    root=args.scratch
    measurement=root/'r4_local_refinement01/first_step_bounds.json'
    w=json.loads(measurement.read_bytes())['bounds']['best_J_witness']
    prep=root/'supercore_joined_front01/adopt309/state.npz'
    endpoint=root/'r4_whole_mouth01/coarse/endpoint.json'
    states={'child':(root/'r4_local_refinement01/first_step/material.npz','2a970cd6b6ffb4ac263d554666e474493daff79f04eb2ed0641f32a132697ddc'),
            'cp10':(root/'r4_melted_survival01/coarse/state_00010.npz','51727c2d48eea199c29bacbfd4bf07b7b11a3e997a737ebcc8c2f69925cd112b')}
    out={}
    for name,(path,pin) in states.items():
        if digest(path)!=pin:raise ValueError('Unexpected material state')
        with np.load(path) as f:
            if str(f['state_sha'])!=digest(prep) or int(f['winding'])!=1 or int(f['radius'])!=4:
                raise ValueError('Wrong R4 state')
            if 'shape_params' in f:
                pars,extra=f['shape_params'],float(f['shape_extra'])
            else:
                with np.load(root/'r4_feasible_descent01/coarse/stalled_00011.npz') as old:
                    if not np.array_equal(old['x'],f['x']):raise ValueError('Changed shape')
                    pars,extra=old['shape_params'],float(old['shape_extra'])
            shape=Meridian(4.,pars,extra)
            g=axes_model(shape,json.loads(endpoint.read_bytes())['depth'])
            body=Solid(prep,shape,float(f['x'][5]),8,24,6,g,degree=3,subdivisions=2)
            body.buffers['base'].fill(0)
            sub=f['surplus'] if 'surplus' in f else np.zeros(3)
            # If the corresponding material point is melted at CP10, take an
            # interior ordered sample of the same prism, fixed across resolution.
            locations=[w,dict(cell_t_x_angle=[5,18,0],triangle=0,barycentric=[.5,.25,.25],angular_fraction=.2)]
            for probe in locations:
                _,point=exact_coordinate(body,f['solid'],sub,probe)
                if _Support(body,f['field'][2]).at(point)>0: break
            else:raise ValueError('No active comparison in CP10 prism')
            out[name]=test(body,f['solid'],sub,f['field'][2],probe)
            out[name]['sha256']=pin
    with args.output.open('x') as h:json.dump(out,h,indent=2,allow_nan=False)
    print(json.dumps({n:dict(J=v['truth']['J'], chi=v['truth']['chi'],
        inverse_norm=v['target_inverse_norm'],required_span=v['levels'][0]['max_nodal_span_needed_for_exact_stretch'],
        readouts=[(r['factor'],r['inverse_sheet_counts'],r.get('Q1_saved_melt_label_continuation',{}).get('J'))
                  for r in v['levels']]) for n,v in out.items()}),flush=True)


if __name__=='__main__':
    main()
