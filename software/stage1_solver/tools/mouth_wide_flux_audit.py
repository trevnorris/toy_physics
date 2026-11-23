"""Independent flux partition and signed original-field energy check on widening exits."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_circulation_material import profile
from mouth_poly_shallow_release import emit
from mouth_sealed_shape_relax import Meridian
from mouth_winding_raw_work import amplitudes_at, topology
from mouth_wound_carrier import prepared


def check(folder, raw_work, step):
    verified=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if (not verified['limitsVerified'] or verified['exitCode'] or verified['childOutcome']['exitCode']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Invalid carrier guard {folder}')
    archive=folder/'fields.npz'
    with np.load(archive,allow_pickle=False) as fields:
        radius=int(fields['radius'])
        degree=int(fields['winding'])
        if radius not in (1,2) or abs(degree)!=1:
            raise ValueError('Unexpected preparation')
        config=fields['shape_params']
        shape=Meridian(radius, config)
        b=shape.exit
        disk=dict(zip(map(str,fields['disk_columns']),fields['disk'].T,strict=True))
        area=disk['area']
        n=disk['density']
        j=disk['j_w']
        s=disk['s']
        incoming=area*np.maximum(j,0.)
        if np.min(n)<=0 or min(area)<=0 or not np.isclose(area.sum(),4*np.pi*b**3/3,rtol=1e-10):
            raise ValueError('Invalid densities, quadrature weights, or disk three-volume')
        norm_ma=np.abs(j)/n**3
        total_ma=np.sqrt((j/n)**2+disk['v_r']**2+disk['v_polar']**2+(degree/s)**2)/n**2
        if not np.allclose(disk['v_phi']*s,degree,rtol=1e-12) or np.min(total_ma-norm_ma)<-1e-10:
            raise ValueError('Incorrect azimuthal field or speed')
        def fractions(mask):
            return {'area':float(area[mask].sum()/area.sum()),
                    'intake':float(incoming[mask].sum()/incoming.sum())}
        emit(kind='flux_audit', radius=radius, exit=b, mesh=int(fields['mesh']),
             samples=len(n), min_n=float(n.min()),max_n=float(n.max()),
             min_normal_Mach=norm_ma.min(), min_total_Mach=total_ma.min(),
             net_mass_rate=area@j, backwards_rate=-area@np.minimum(j,0),
             dense=fractions(n>=.5),normal_subsonic=fractions(norm_ma<1),
             dense_normal_subsonic=fractions((n>=.5)&(norm_ma<1)),
             total_subsonic=fractions(total_ma<1),
             dense_total_subsonic=fractions((n>=.5)&(total_ma<1)),
             disk_volume=area.sum(),input_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
        if not raw_work:
            return
        n,nt=int(fields['mesh']),int(fields['polar'])
        outer,roof=float(fields['outer']),float(fields['roof'])
        amplitude,phi=fields['amplitude_ratio'],fields['phase']
        material_sha=str(fields['material_sha256'])
    _,chi,sha=prepared(radius)
    if sha!=material_sha:
        raise ValueError('Unmatched crystal geometry')
    representation=profile(abs(degree),max(80.,2*outer),2.e-9)[0]
    phase_labels,amplitude_labels=topology(shape,n,nt,outer,roof)
    opts=(radius,n,nt,outer,roof,degree,5)
    work=[]
    for sign in (-1,1):
        coords=config.copy()
        coords[1]+=np.log1p(sign*step)
        probe=Meridian(radius,coords)
        work.append(amplitudes_at(phi,probe,amplitude,opts,phase_labels,amplitude_labels,chi,representation))
    emit(kind='independent_carrier_work',radius=radius,exit=b,mesh=n,polar=nt,h=step,
         raw_energy_force=(work[0]['continuous_raw']-work[1]['continuous_raw'])/(2*b*step),
         matched_cell_averages=(work[0]['original_with_native_lumping']-work[1]['original_with_native_lumping'])/(2*b*step))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directories',type=Path,nargs='+')
    parser.add_argument('--raw-work',action='store_true')
    parser.add_argument('--step',type=float,default=.003)
    args=parser.parse_args()
    for directory in args.directories:
        check(directory,args.raw_work,args.step)
