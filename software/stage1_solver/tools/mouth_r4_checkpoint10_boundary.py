"""Locate every folded material prism touched by the rejected spatial Q1 order basis at R4.

Read-only diagnostic on the saved 24 x 24 checkpoint. The production full-prism
certificate can over-enclose a folded region; a positive basis at an actual
negative-determinant witness is a stricter test.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys
from types import SimpleNamespace

import mpmath as mp
import numpy as np


SHA10 = "51727c2d48eea199c29bacbfd4bf07b7b11a3e997a737ebcc8c2f69925cd112b"
OFFSETS = (((0,0),(0,1),(1,0)),((1,1),(1,0),(0,1)))


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(run, build):
    pin = run/"pin/software/stage1_solver/tools"
    sys.path.insert(0,str(pin))
    from mouth_active_prism_orientation import (
        _Support, _certify_wedge, _coefficients, _down, _easy_wedges,
        _exact_array, _quadratic_min, _rational, _up, UNIT,
    )
    from mouth_continuous_order import topology
    from mouth_free_shape_solid import Solid, pointer
    from mouth_r4_fold_order import axes_model
    from mouth_r4_fold_constitutive import Prism, native_callback, point
    from mouth_sealed_shape_relax import Meridian
    from mouth_wound_carrier import straight_vortex

    for name,checksum in read(run/"sources.json").items():
        if digest(run/"pin"/name)!=checksum:
            raise ValueError(f"Changed producer source {name}")
    folder=run/"coarse"
    source=folder/"state_00010.npz"
    if digest(source)!=SHA10:
        raise ValueError("Checkpoint 10 hash differs")
    logs=[json.loads(row) for row in (folder/"stdout").read_text().splitlines() if row.startswith("{")]
    saved={row["iteration"]:row for row in logs if row["kind"]=="accepted"}
    assert saved[10]["sha256"]==SHA10
    residuals={row["iteration"]:row for row in logs if row["kind"]=="residuals"}
    trace=[]
    for i in (8,9,10):
        s,r=saved[i],residuals[i]
        assert digest(folder/s["file"])==s["sha256"]
        trace.append(dict(iteration=i,state_sha256=s["sha256"],forces=r["shape_force"],depth=s["depth"],
                          exit=s["exit_radius"],intake=r["current"]["feed_current"],
                          dense_subsonic=r["current"]["dense_total_subsonic_share"]))
    after=logs[next(n for n,row in enumerate(logs) if row["kind"]=="residuals" and row["iteration"]==10)+1:]
    rejected=[row for row in after if row["kind"]=="rejected_shape_rematch"]
    ids={int(re.search(r"joined node (\d+)",row["detail"]).group(1)) for row in rejected}
    if len(ids)!=1:
        raise ValueError(f"Ambiguous rejected order node {ids}")
    target=ids.pop()
    old=Path(read(folder/"invocation.json")["command"][-3])/"coarse"
    old_cmd=read(old/"invocation.json")["command"]
    original=Path(old_cmd[-5])/"coarse"
    endpoint=original/"endpoint.json"
    prep=Path(read(original/"invocation.json")["command"][5])
    proof=read(Path(old_cmd[-3]))["results"]["coarse"]
    assert digest(prep)==proof["preparation_sha256"] and digest(endpoint)==proof["endpoint_sha256"]
    with np.load(source,allow_pickle=False) as f:
        assert int(f["winding"])==1 and int(f["radius"])==4 and str(f["state_sha"])==digest(prep)
        shape=Meridian(4.,f["shape_params"],float(f["shape_extra"]))
        shift,solid,u,q=float(f["x"][5]),f["solid"],f["field"],f["q"]
    geom=axes_model(shape,read(endpoint)["depth"])
    b=Solid(prep,shape,shift,8,24,6,geom,degree=3,subdivisions=2)
    b.buffers["base"].fill(0.)
    assert np.array_equal(b.unpack(q),solid)
    b.data.offset=pointer(solid)
    correction=np.ascontiguousarray(u[2,b.ids])
    b.data.correction=pointer(correction)
    numbers,_=topology(geom)
    indicator=(numbers==target).astype(float)
    hit=_Support(b,indicator)
    alive=_Support(b,u[2])
    locs=np.argwhere(indicator[b.ids])
    assert len(locs)==1
    rad,ang,ax=locs[0]
    joined=dict(id=target, field_index_x_theta_t=locs[0].tolist(), x=float(geom.rr[0][rad]),
                theta=float(geom.a[ang]), w=float(geom.z[1][ax]), chi=float(u[2,b.ids[rad,ang,ax]]))
    lattice=np.stack(np.meshgrid(b.t,b.x,b.a,indexing="ij"),axis=-1)[...,[1,0,2]]
    coords=lattice+solid
    trial=SimpleNamespace(b=b,model=SimpleNamespace(_vortex=straight_vortex(1,28.)))
    callback=native_callback(pin,build)
    mp.mp.dps=65
    def idx(j,i,k,tri):
        return [(j+dj,i+di,k+h) for h in (0,1) for dj,di in OFFSETS[tri]]
    def planes(indices):
        ix=tuple(np.array(indices).T)
        return (_exact_array(lattice[ix])+_exact_array(solid[ix])).reshape(2,3,3)
    def actual_witness(cell,bary,y):
        j,i,k,tri=cell
        where=point(Prism(b,solid,(j,i,k),tri),bary,y,trial,u,callback)
        return {name:where[name] for name in ("barycentric","angular_fraction","physical_r_z_theta",
                "landed_x_t_theta","J","det_chart","principal_stretches","chi","n")}

    def witness_state(p,scale,support):
        counts=dict(exact_subprisms=0,max_depth_used=0,inactive_fold_enclosures=0)
        cover=set()
        result,witness=_certify_wedge(p,scale,support,cover,counts,16)
        if result=="unresolved":
            raise ValueError("Fold/support intersection unresolved")
        return result,witness,cover

    def incident(vertex):
        t,x,a=vertex
        for j in (t-1,t):
            for i in (x-1,x):
                for k in (a-1,a):
                    if not (0<=j<len(b.t)-1 and 0<=i<len(b.x)-1 and 0<=k<len(b.a)-1):
                        continue
                    for tri in (0,1):
                        keys=idx(j,i,k,tri)
                        if vertex in keys:
                            yield (j,i,k,tri),keys

    def vertex_record(node):
        v=planes([node]*6)[0,0]
        samples=[]
        # Identify active neighboring prisms by lower-bound witnesses. An
        # all-zero enclosing support certifies a prism inactive. Any remaining
        # positive enclosure with no positive witness is explicitly unresolved.
        for cell,keys in incident(node):
            p=planes(keys)
            if not alive.positive[alive.box(p)].any():
                continue
            witnesses=[(alive.at(pos),pos) for pos in p.reshape(6,3)]
            mid=p.sum(axis=(0,1))/6
            witnesses.append((alive.at(mid),mid))
            val,pos=max(witnesses,key=lambda pair:pair[0])
            samples.append(dict(cell=list(cell),chi=float(val),at_x_t_theta=list(map(float,pos))))
        best=max(samples,key=lambda v:v["chi"],default=None)
        return dict(node_t_x_theta=list(node),landed=list(map(float,v)),chi=float(alive.at(v)),
                    ordered_neighbor=best,incident_prisms=len(list(incident(node))),
                    classification="shared_active" if best and best["chi"]>0 else
                    ("unresolved" if samples else "strictly_melted"))

    folded=[]
    cache={}
    for tri,offsets in enumerate(OFFSETS):
        corners=np.stack([coords[j:j+len(b.t)-1,i:i+len(b.x)-1] for j,i in offsets],axis=-2)
        lo,hi=_down(corners),_up(corners)
        easy=_easy_wedges((lo[:,:,:-1],hi[:,:,:-1]),(lo[:,:,1:],hi[:,:,1:]))
        for j,i,k in np.argwhere(~easy):
            j,i,k=map(int,(j,i,k))
            cell=(j,i,k,tri)
            nodes=idx(*cell)
            p=planes(nodes)
            mins=[_quadratic_min(row) for row in _coefficients(p)]
            val,vertex,y=min((v,v_id,y) for v_id,(v,y) in enumerate(mins))
            if val>0:
                continue
            scale=Fraction(1)
            for axis,index in ((b.x,i),(b.t,j),(b.a,k)):
                scale*=_rational(axis[index+1])-_rational(axis[index])
            status,_,cover=witness_state(p,scale,alive)
            assert status=="certified"
            result,at,_=witness_state(p,scale,hit)
            vlist=[]
            for node in nodes:
                if node not in cache:
                    cache[node]=vertex_record(node)
                vlist.append(cache[node])
            item=dict(cell_t_x_angle_triangle=list(cell),certificate_includes_target=target in numbers[list(cover)],
                      vertices=vlist,minimum=actual_witness(cell,UNIT[vertex],y),vertex_of_minimum=vertex)
            if item["certificate_includes_target"]:
                item["ordered_vertex_states"]=[actual_witness(cell,UNIT[h%3],h//3)
                        for h,v in enumerate(vlist) if v["chi"]>0.]
            if result=="active_fold":
                item["target_witness"]=actual_witness(cell,at["barycentric"],at["angular_fraction"])
                item["target_basis_weight"]=at["order"]
            folded.append(item)
    return dict(checkpoint_sha256=SHA10,target=joined,rejections=rejected,trace=trace,folded_prisms=folded,
                coordinate_convention="zero-based [depth, radial, polar, triangle]. The six vertices are the bottom triangle followed by the top triangle.",
                local_stretches="magnitudes of singular values plus hoop, evaluated from saved labels; signed model volume ratio reported separately",
                orientation=saved[10]["orientation"],field=residuals[10]["field"]["raw_residual"])


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run",type=Path)
    ap.add_argument("build",type=Path)
    ap.add_argument("result",type=Path)
    args=ap.parse_args()
    result=audit(args.run.resolve(),args.build.resolve())
    args.result.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("folds, target overlaps, conservative enclosures",len(result["folded_prisms"]),
          sum("target_witness" in k for k in result["folded_prisms"]),
          sum(k["certificate_includes_target"] for k in result["folded_prisms"]))
