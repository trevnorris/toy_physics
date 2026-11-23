"""Direct finite-step energy of the strips that change order element.

Integrate between each original carrier face and the exact preimage of the
corresponding moving material face. At intersections, switch the radial
order element first and the axial element second to avoid double counting.
The carrier field always comes from its own original element.
"""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_closed_rim_surface_review import Reader
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_rim_disk_work_review import Disk
from mouth_rim_mixed_face_work import projection, traces
from mouth_shared_rim_shift import SharedRimChart
from mouth_supercore_joined_front import file_sha


def radius_at(shape,bank,x,z):
    wall = shape.radius(z)[0]
    return wall+x*(28-wall) if bank else x*wall


def radial_coordinate(shape,bank,r,z):
    wall = shape.radius(z)[0]
    return (r-wall)/(28-wall) if bank else r/wall


def material_depth(material,bank,r,target,step):
    z = target.copy()
    for _ in range(5):
        motion,_,dz = material(bank,r,z)
        z -= (z+step*motion-target)/(1+step*dz)
    if np.max(abs(z+step*material(bank,r,z)[0]-target)) > 2e-13:
        raise ArithmeticError("Material map inversion failed in a swept layer")
    return z


def moved_gradient(g,dr,dz,step):
    result = g.copy()
    result[:,2] /= 1+step*dz
    result[:,0] -= step*dr*result[:,2]
    return result


def finite_layer(disk,chart,material,carrier,step):
    shape, model = disk.sh, disk.model
    q = (.5-.5/np.sqrt(3.),.5+.5/np.sqrt(3.))
    work = np.zeros((2,2,2))
    for slab,bank,axis in product((1,2),(0,1),(0,2)):
        grids = (model.rr[bank],model.a,model.z[slab])
        ranges = [np.arange(len(g)-1) for g in grids]
        ranges[axis] = np.arange(1,len(grids[axis])-1)
        j = ranges[2]
        ranges[2] = j[(grids[2][j+int(axis!=2)]>chart.tents[0]) & (grids[2][j]<chart.tents[-1])]
        idx = tuple(t.ravel() for t in np.meshgrid(*ranges,indexing="ij"))
        for parameters in product(q,repeat=2):
            xi = iter(parameters)
            coords, measure = [], np.ones(len(idx[0]))
            for d,(g,i) in enumerate(zip(grids,idx)):
                if axis == d:
                    coords.append(g[i])
                else:
                    h = np.diff(g)[i]
                    coords.append(g[i]+h*next(xi))
                    measure *= h/2
            x0,a0,z0 = coords
            r0 = radius_at(shape,bank,x0,z0)
            relative = carrier(bank,r0,z0)[0]-material(bank,r0,z0)[0]
            take = relative!=0
            if not take.any():
                continue
            x,a,z = [u[take] for u in coords]
            face_indices = [i[take] for i in idx]
            face = x if axis==0 else z
            preimage = face.copy()
            for _ in range(5):
                if axis==0:
                    r = radius_at(shape,bank,preimage,z)
                    target = z+step*carrier(bank,r,z)[0]
                    zm = material_depth(material,bank,r,target,step)
                    preimage = radial_coordinate(shape,bank,radius_at(shape,bank,x,zm),z)
                else:
                    r = radius_at(shape,bank,x,preimage)
                    preimage = z+step*(material(bank,r,z)[0]-carrier(bank,r,preimage)[0])
            sweep = preimage-face
            carrier_idx = list(face_indices)
            carrier_idx[axis] = carrier_idx[axis]-(sweep<0)
            for v in q:
                cx,cz = (face+sweep*v,z) if axis==0 else (x,face+sweep*v)
                r = radius_at(shape,bank,cx,cz)
                V,Vr,Vz = carrier(bank,r,cz)
                zm = material_depth(material,bank,r,cz+step*V,step)
                xm = radial_coordinate(shape,bank,r,zm)
                _,Sr,Sz = material(bank,r,zm)
                actual = list(carrier_idx)
                for d,points in ((0,xm),(2,zm)):
                    actual[d] = np.clip(grids[d].searchsorted(points,side="right")-1,0,len(grids[d])-2)
                before = list(carrier_idx)
                after = list(carrier_idx)
                if axis==2:
                    before[0] = after[0] = actual[0]
                after[axis] = actual[axis]
                if np.any(abs(after[axis]-before[axis])>1):
                    raise ValueError("More than one cell crossed")
                R, slope = shape.radius(cz)
                slope_m = shape.radius(zm)[1]
                measure_on_layer = (measure[take]*abs(sweep)/2 * (1+step*Vz) *
                                    4*np.pi*r*r*np.sin(a)*(28-R if bank else R))
                F = disk.reader.vortex(r*np.sin(a))[0]
                for s, interp in zip((1.,-1.),disk.reader.fields[slab,bank]):
                    value, g = traces(model,interp.values,carrier_idx,(cx,a,cz),bank,r,slope)
                    gp = moved_gradient(g[:,0],Vr,Vz,step)
                    n = (value[:,0]*F)**2
                    def branch(index):
                        _, field_grad = traces(model,interp.values,index,(xm,a,zm),bank,r,slope_m)
                        return projection(n,gp,moved_gradient(field_grad[:,1],Sr,Sz,step))
                    work[int(axis==2),slab-1,bank] -= s*(measure_on_layer @ (branch(after)-branch(before)))/step
    return float(work.sum())


def run(args):
    disk = Disk(args.state,args.wall,args.fields)
    front = Front(args.state,GRIDS["native"],1)
    x,d,_,_ = read_source(args.state,args.wall)
    native = Trial(front,args.state,x,(d.shape[1]-1)//3,d.shape[0]-1,False).model
    chart = SharedRimChart(disk.model,native,args.axes)
    oracle = json.loads(args.faces.read_text())
    if oracle["field"] != file_sha(args.zero) or oracle["wall"] != file_sha(args.wall):
        raise ValueError("Mismatched zero state")
    expected = oracle["readings"][0]["cases"]
    with np.load(args.zero) as f:
        disk.reader = Reader(disk.model,front.c,f["wound"],f["rest"])
    rows = []
    for name,n in (("fixed_wall_difference",0),("rim_four",4),("rim_eight",8)):
        material = chart.physical if n else chart.still
        carrier = (lambda b,r,z: chart.carrier(b,r,z,n)) if n else chart.stationary
        for sign in (-1,1):
            limit = expected[name]["negative" if sign<0 else "positive"]
            for h in (2.e-5,1.e-5):
                force = finite_layer(disk,chart,material,carrier,sign*h)
                rows.append(dict(case=name,step=sign*h,exact_layer=force,linear_limit=limit,difference=force-limit))
    return dict(grid=disk.level,zero=oracle["field"],measurements=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("state","wall","fields","axes","zero","faces"):
        p.add_argument(name,type=Path)
    print(json.dumps(run(p.parse_args()),allow_nan=False))
