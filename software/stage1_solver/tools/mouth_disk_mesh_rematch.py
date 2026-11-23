"""Rematch saved radius-four reservoir endpoints on a new field grid, at fixed solid displacement."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_bonded_front_adjoint import Coupled, get_material_load
from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_joined_front_stationarity import solve
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    return list(map(json.loads,(path.parent/'stdout').read_text().splitlines()))


def inputs(state, joint, fields):
    for path in (state, joint, fields):
        check_guard(path)
    source = rows(joint)
    control = rows(fields)
    if not (source[0]['source_sha256'] == control[0]['source_sha256'] == digest(state)
            and source[-1]['sha256'] == control[0]['joint_sha256'] == digest(joint)
            and control[-1]['sha256'] == digest(fields)):
        raise ValueError('Geometry, bonded displacement and matched control must have common provenance')
    with np.load(joint) as archive:
        displacement = archive['displacement']
    with np.load(fields) as archive:
        pair = [archive[key] for key in ('fields','rest')]
    base = Front(state,tuple(control[0].get('grid',source[0]['grid'])),control[0]['winding'])
    body = Bonded(state,(displacement.shape[1]-1)//3,displacement.shape[0]-1,full_depth=True)
    return base, body, displacement, pair


def remap(old, new, field):
    out = np.zeros((3,new.size))
    for (slab,side), ids in new.ids.items():
        axes = (old.rr[side],old.a,old.z[slab])
        values = field[:,old.ids[slab,side]].transpose(1,2,3,0)
        mesh = np.meshgrid(new.rr[side],new.a,new.z[slab],indexing='ij')
        out[:,ids] = RegularGridInterpolator(axes,values)(np.stack(mesh,axis=-1)).transpose(3,0,1,2)
    return out


def run(args):
    base,body,displacement,pair = inputs(args.state,args.joint,args.fields)
    grid = tuple(args.grid)
    if grid == tuple(base.grid):
        raise ValueError('This tool only solves a new field grid')
    old = base.joined(base.sh,False)
    new = Coupled(base.sh,base.c,base.N,grid,terminal_profile=old.terminal_profile)
    q = displacement[body.mask]*body.scale
    new.material_load = get_material_load(body,q,new)
    new.warm_path = None
    cold = new.seed.copy()
    result = []
    emit(kind='rematch',source_sha256=digest(args.state),joint_sha256=digest(args.joint),
         winding=base.N,grid=grid,seed_grid=base.grid,seed_sha256=digest(args.fields))
    for rest, field in enumerate(pair):
        new.seed = cold.copy()
        new.seed[new.masks] = remap(old,new,field)[new.masks]
        updated, info = solve(new,bool(rest))
        result.append(updated)
        emit(kind='endpoint',rest=bool(rest),**info,**new.measure(updated))
    np.savez_compressed(args.save,fields=result[0],rest=result[1])
    emit(kind='fields',sha256=digest(args.save))


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for key in ('state','joint','fields','save'):
        ap.add_argument(key,type=Path)
    ap.add_argument('--grid',type=int,nargs=3,required=True)
    run(ap.parse_args())
