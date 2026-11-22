"""Integrate saved occupied-mode kinetic energy either side of the internal constriction."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar
from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import Meridian, DELTA

BASE=Path('_scratch/s11c/live_waist_dynamics01')


def share(wave,shape,field,weights,mesh,zlo,zhi):
    """Integrate a truncated Q1 element using its actual four saved nodal amplitudes."""
    n_ax=len(weights)//mesh
    x=np.arange(mesh)[:,None]/mesh
    w0=np.arange(n_ax)[None,:]*shape.length/n_ax
    lower=np.maximum(w0,zlo)
    width=np.maximum(0,np.minimum(w0+shape.length/n_ax,zhi)-lower)
    local=np.r_[field,0j][wave.ids].reshape(mesh,n_ax,4)
    total=np.zeros((mesh,n_ax))
    nodes,quad=leggauss(3)
    for aa,wa in zip((nodes+1)/2,quad/2):
        s=x+aa/mesh
        for b,wb in zip((nodes+1)/2,quad/2):
            w=lower+b*width
            bb=(w-w0)*n_ax/shape.length
            radius=shape.radius(w)[0]
            v=((1-aa)*(1-bb)*local[:,:,0]+aa*(1-bb)*local[:,:,1]
                +(1-aa)*bb*local[:,:,2]+aa*bb*local[:,:,3])
            total+=wa*wb*width*radius/mesh*(s*radius)**4*abs(v)**2
    return np.sum(total*weights.reshape(mesh,n_ax))


def read(case):
    path=BASE/case
    inv=json.loads((path/'invocation.json').read_bytes())
    lines=[json.loads(s) for s in (path/'stdout').read_text().splitlines()]
    prep=next((row for row in lines if row.get('kind')=='preparation'),None)
    if prep is None:
        from mouth_live_waist_order import state
        command=inv['command']
        def option(name):
            return command[command.index(name)+1]
        if '--radius' in command:
            radius=int(option('--radius'))
        elif Path(command[1]).name=='mouth_live_waist_support.py':
            radius=1
        else:
            raise ValueError('Saved field has no radius in its preparation')
        prep=dict(radius=radius, field_mesh=int(option('--field-mesh')), mesh=int(option('--mesh')),
                  coordinates=state(radius,option('--point'))['coordinates'])
    snapshot=Path(inv['command'][inv['command'].index('--save')+1])
    field=np.load(snapshot)
    center=[x for x in lines if (x.get('action')==1 or x.get('kind')=='bounded') and 'wave_frequency' in x]
    if not center or abs(field['frequency']-center[-1]['wave_frequency'])>1e-9:
        raise ArithmeticError(f'Saved center cannot be verified for {case}')
    m,n=prep['field_mesh'],prep['mesh']
    shape=Meridian(prep['radius'],prep['coordinates'])
    obj=OccupiedWave(shape,n,2*n,m,4.2)
    weights=obj.project@field['order'].T.ravel()
    v=field['vector']
    cut=minimize_scalar(lambda z:shape.radius(z)[0],bounds=(0,shape.length),method='bounded',
                        options={'xatol':1e-12}).x
    direct=share(obj,shape,v,weights,m,0,shape.length)
    vec=np.r_[v,0j][obj.ids]
    stored=weights@np.einsum('ci,cij,cj->c',vec.conj(),obj.mc,vec).real
    if abs(direct/stored-1)>2e-10:
        raise ArithmeticError('Independent direct quadrature fails full-pocket mass identity')
    edges=[0., max(0,cut-2*DELTA), cut,min(shape.length,cut+2*DELTA),shape.length]
    fractions=np.array([share(obj,shape,v,weights,m,a,b)/stored for a,b in zip(edges[:-1],edges[1:])])
    if abs(sum(fractions)-1)>1e-5:
        raise ArithmeticError('Kinetic fraction from cut cells does not partition pocket')
    print(json.dumps(dict(case=case,snapshot=str(snapshot),sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
                       radius=prep['radius'],frequency=float(field['frequency']),pocket_mass_norm=stored,
                       waist_depth=cut,waist_radius=float(shape.radius(cut)[0]),band_half_width=2*DELTA,
                       strictly_upstream_of_center=fractions[:2].sum(),strictly_downstream_of_center=fractions[2:].sum(),
                       in_waist_band=fractions[1:3].sum(),beyond_downstream_band=fractions[3],
                       inside_inlet_side_before_band=fractions[0],fraction_partition_error=sum(fractions)-1)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('case',nargs='+')
    for case in p.parse_args().case: read(case)
