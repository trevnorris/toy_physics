"""Intersect the saved step-167 Jacobi line with axial planes on the welded collar.

The two coordinates of a ray endpoint are affine in the search fraction. An
interior radial/angular cut meets the depth plane if the determinant of the
two endpoint offsets is zero (quadratic in that fraction). At a fixed depth
a registry triangle diagonal also has a fixed radial coordinate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

LD = np.longdouble
CAP = LD(1)/8192


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def lerp(a,b,f):
    return a+f*(b-a)


def endpoints(arr,t,a):
    return lerp(lerp(arr[:,:-1,:-1],arr[:,1:,:-1],t),
                lerp(arr[:,:-1,1:],arr[:,1:,1:],t),a)


def roots(c0,c1,c2):
    if c2==0:
        return [] if c1==0 else [-c0/c1]
    discr=c1*c1-4*c0*c2
    if discr<=0:
        return []  # tangent contact does not cross the plane
    q=-LD('.5')*(c1+np.copysign(np.sqrt(discr),c1))
    return [q/c2,c0/q]


def radial_levels(body,plane,weld):
    # A diagonal can cross a horizontal reference edge only at a registry
    # vertex. Interior native depth planes can also intersect a diagonal.
    x,e,z=(body.buffers[k] for k in ('x','reg_e','reg_t'))
    out={}
    def add(v,kind):
        out.setdefault(LD(v),[]).append(kind)
    for v in x:
        if v>=weld:
            add(v,'native radial')
    p=LD(body.data.pressure)
    for eta in e:
        v=np.exp(-p*LD(eta))
        if v>=weld:
            add(v,'registry radial')
    if plane in z:
        for key, kinds in out.items():
            if 'registry radial' in kinds:
                kinds.append('registry diagonal junction')
    else:
        j=np.searchsorted(z,plane)-1
        f=(plane-z[j])/(LD(z[j+1])-z[j])
        for a,b in zip(e,e[1:]):
            eta=LD(a)+(1-f)*(LD(b)-a)
            v=np.exp(-p*eta)
            if v>=weld:
                add(v,'registry diagonal')
    return out


def measure(area, ceiling=CAP, floor=LD(0)):
    from mouth_r4_aligned_derivatives import load
    from mouth_r4_eulerian_release import order_grid
    from mouth_cp10_label_condition import physical_frame

    study=area.parent/'r4_segmented_held01'
    start=study/'primary/step_00166.npz'
    end=json.loads(start.with_suffix('.json').read_text())
    direction=area/'directions/directions.npz'
    proof=json.loads((area/'directions/result.json').read_text())
    if not proof['passed'] or digest(start)!=end['sha256']:
        raise ValueError('Wrong step or unverified Jacobi direction')
    dr,tr,solid,_=load(area.parent/'r4_rematch39_01','primary',area/'crossings')
    with np.load(start,allow_pickle=False) as st,np.load(direction,allow_pickle=False) as d:
        if (not np.array_equal(d['labels'],st['labels']) or not np.array_equal(st['x'],dr.x)
                or int(st['radius'])!=4 or int(st['winding'])!=1):
            raise ValueError('Direction, source, or held shape mismatch')
        labels=st['labels'].astype(LD)
        rate=d['pj'].astype(LD)
        chi=order_grid(tr,st['field'])
    dims=np.array(chi.shape)-1
    active=np.max(np.lib.stride_tricks.sliding_window_view(chi,(2,2,2)),axis=(-3,-2,-1))>0
    tx,regt=tr.b.buffers['t'],tr.b.buffers['reg_t']
    planes=np.unique(np.r_[tx,regt])[1:-1].astype(LD)
    weld=LD(tr.b.buffers['x'][(len(tr.b.buffers['x'])-1)*2//3])
    angle=tr.b.buffers['a']
    peak=end['material']['worst_free_active']
    axes=solid.axes
    events={}
    identically_at_endpoint=0
    diagnostics=dict(ray_axial_plane_candidates=0,split_levels_examined=0,roots_at_baseline=0)
    def add(beta,i,j,k,v,a,it,ia,u,x0,t0,angle0,plane,rule,cellactive,c):
        alpha=beta*ceiling
        if alpha<=floor:
            return
        if u==0 or u==1:
            key=('endpoint',int(i+u),j,k,it,ia,float(plane))
        else:
            key=('split',i,j,k,it,ia,float(plane),float(beta),float(u))
        if key in events:
            ev=events[key]
            ev['cells'].append([i,j,k])
            ev['ordered_cell']|=bool(cellactive)
            ev['types']=sorted(set(ev['types']+rule))
            return
        local=np.array([u,v,a],LD); cell=[i,j,k]
        point=np.array([ax[n]+f*(ax[n+1]-ax[n]) for n,f,ax in zip(cell,local,axes)],float)
        # Physical radius is component B_theta,theta including at the axis.
        frame,_=physical_frame(tr.b,point[None]);rad=float(frame[0,2,2])
        theta=float(point[2]);depth=float(point[1]*tr.shape.length)
        distance=np.sqrt((depth-peak['depth'])**2+(rad-peak['radius'])**2
                         +4*rad*peak['radius']*np.sin((theta-peak['theta'])/2)**2)
        events[key]=dict(alpha=float(alpha),cell=cell,cells=[cell],local=local.astype(float).tolist(),
            ray=[it,ia],knot=key[0],types=rule,plane=float(plane),material_x_t_theta=list(map(float,(x0,t0,angle0))),
            r_w_theta=[rad,depth,theta],meridional_angular_distance_to_peak=float(distance),
            absolute_depth_distance_to_peak=abs(depth-peak['depth']),ordered_cell=bool(cellactive),chi=float(c))
    for it,v in enumerate(np.array(solid.points,LD)):
        for ia,a in enumerate(np.array(solid.points,LD)):
            val=endpoints(labels,v,a)
            vel=ceiling*endpoints(rate,v,a)
            cval=endpoints(chi.astype(LD),v,a)
            # Bounding boxes exact in real arithmetic: each ray coordinate is
            # bilinear in (radial fraction, beta=alpha/CAP).
            lo=np.minimum.reduce((val[:-1],val[1:],val[:-1]+vel[:-1],val[1:]+vel[1:]))
            hi=np.maximum.reduce((val[:-1],val[1:],val[:-1]+vel[:-1],val[1:]+vel[1:]))
            scale=LD(16)*np.finfo(LD).eps
            for T in planes:
                tl=endpoints(labels[...,1]-T,v,a)
                # Compute offsets before interpolation to retain a tiny
                # vertex-plane distance without cancellation of absolute depth.
                z0,z1=tl[:-1],tl[1:];vz0,vz1=vel[:-1,...,1],vel[1:,...,1]
                zlo=np.minimum.reduce((z0,z1,z0+vz0,z1+vz1))
                zhi=np.maximum.reduce((z0,z1,z0+vz0,z1+vz1))
                mask=(zlo<=0)&(zhi>=0)&(hi[...,0]>=weld-scale)
                cand=np.argwhere(mask)
                if not len(cand):
                    continue
                diagnostics['ray_axial_plane_candidates']+=len(cand)
                rel=radial_levels(tr.b,T,weld)
                for idx in cand:
                    i,j,k=map(int,idx)
                    left,right=val[i,j,k],val[i+1,j,k]
                    dl,drt=vel[i,j,k],vel[i+1,j,k]
                    P=tl[i,j,k];Q=tl[i+1,j,k]
                    A=dl[1];B=drt[1]
                    for u,P0,D0,X,D in ((0,P,A,left,dl),(1,Q,B,right,drt)):
                        if P0==0:
                            identically_at_endpoint+=int(D0==0)
                            diagnostics['roots_at_baseline']+=int(D0!=0)
                        if D0==0:
                            continue
                        beta=-P0/D0
                        if not (0<beta<=1):
                            continue
                        xx=X+beta*D
                        if xx[0]<weld:
                            continue
                        add(beta,i,j,k,v,a,it,ia,LD(u),*xx,T,['ray endpoint'],active[i,j,k],cval[i+u,j,k])
                    levels=[(0,x,why) for x,why in rel.items() if lo[i,j,k,0]-scale<=x<=hi[i,j,k,0]+scale]
                    levels.extend((2,LD(x),['native angle']) for x in angle
                                  if lo[i,j,k,2]-scale<=x<=hi[i,j,k,2]+scale)
                    for m,target,why in levels:
                        diagnostics['split_levels_examined']+=1
                        C,D=left[m]-target,right[m]-target
                        ca,da=dl[m],drt[m]
                        # P(beta)*D(beta)-Q(beta)*C(beta) = 0. Roots at
                        # endpoints already counted above; coincident traces
                        # that follow the axial plane have no transverse crossing.
                        c0=P*D-Q*C; c1=P*da+A*D-Q*ca-B*C; c2=A*da-B*ca
                        for beta in roots(c0,c1,c2):
                            if not (0<beta<=1):
                                continue
                            f,g=C+beta*ca,D+beta*da
                            if f==g:
                                continue
                            u=-f/(g-f)
                            if not (0<u<1):
                                continue
                            xx=lerp(left+beta*dl,right+beta*drt,u)
                            if xx[0]<weld-scale:
                                continue
                            if abs(lerp(P+beta*A,Q+beta*B,u))>scale:
                                raise ArithmeticError('Unresolved collision residual')
                            c=lerp(cval[i,j,k],cval[i+1,j,k],u)
                            add(beta,i,j,k,v,a,it,ia,u,*xx,T,why,active[i,j,k],c)
    ordered=sorted(events.values(),key=lambda e:e['alpha'])
    loaded=[ev for ev in ordered if ev['ordered_cell']]
    if floor==0 and not loaded:
        raise ValueError('The archived line-search range has no axial reference crossing on active cells')
    first=loaded[0] if floor==0 else json.loads((area/'crossings/result.json').read_text())['first']
    ladder=[]
    for multiple in [.001,.01,.1,.5,.9,1.1,2,4,None]:
        alpha=float(CAP) if multiple is None else first['alpha']*multiple
        if alpha>float(ceiling) or alpha<=float(floor):
            continue
        ladder.append(dict(alpha=alpha,multiple=multiple,crossings=sum(e['alpha']<=alpha for e in loaded),
                           all_geometry_crossings=sum(e['alpha']<=alpha for e in ordered)))
    return dict(state_sha256=digest(start),directions_sha256=digest(direction),tool_sha256=digest(Path(__file__)),
                command=sys.orig_argv,interval_left_open=float(floor),interval_right_closed=float(ceiling),alpha_cap=float(CAP),alpha_cross=first['alpha'],first=first,
                first_all_geometry=ordered[0] if ordered else None,peak=peak,
                planes=planes.astype(float).tolist(),registry_depth_planes=regt.tolist(),weld=float(weld),
                active_cell_count=int(active.sum()),total_cell_count=int(active.size),
                diagnostics=diagnostics,ray_endpoint_incidents_identically_on_plane=identically_at_endpoint,
                definition='Unique physical shared ray endpoints or interior intersections with each distinct radial/angular cut, restricted to the collar side. Events counted on active cells even for zero order at a shared face. Coincident native and material axial planes counted once; tangent contacts and motion along the plane excluded.',
                all_events=ordered,ladder=ladder)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('--extension',action='store_true',help='Only enumerate events beyond the archived line-search fraction through 4 times the first crossing')
    args=p.parse_args();area=args.area.resolve()
    out=area/('crossings_extension' if args.extension else 'crossings')
    out.mkdir(parents=True,exist_ok=True)
    if args.extension:
        original=json.loads((area/'crossings/result.json').read_text())
        record=measure(area,4*LD(original['alpha_cross']),CAP)
        for row in record['ladder']:
            row['crossings']+=sum(ev['ordered_cell'] for ev in original['all_events'])
            row['all_geometry_crossings']+=len(original['all_events'])
        record['prior_sha256']=digest(area/'crossings/result.json')
    else:
        record=measure(area)
    with (out/'result.tmp').open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False);stream.write('\n')
    (out/'result.tmp').rename(out/'result.json')
    print(json.dumps({k:record[k] for k in ('alpha_cross','first','first_all_geometry','ladder','diagnostics')}),flush=True)
