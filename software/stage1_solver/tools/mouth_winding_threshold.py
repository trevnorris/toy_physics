"""Compare wound and unwound exit work on matched spatial grids and reservoir preparations."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_poly_shallow_release import emit
from mouth_winding_speed_gate import preparation


def verified(path):
    outcome=json.loads((path/'outcome.json').read_bytes())
    limits=json.loads((path/'limit-validation.json').read_bytes())
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or
        not outcome['limitsVerified'] or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Unverified {path}')
    raw=(path/'stdout').read_bytes()
    rows=list(map(json.loads,raw.splitlines()))
    return rows,hashlib.sha256(raw).hexdigest()


def compare(paths):
    groups={}
    for path in paths:
        rows,fingerprint=verified(path)
        spec=rows[0]
        if spec['kind']!='input':
            raise ValueError(f'Unexpected source {path}')
        if spec['mode']!='signed':
            raise ValueError('Expected signed exit work')
        shape,_,digest,load=preparation(spec['radius'])
        if digest!=spec['material_sha256'] or abs(shape.exit-spec['reference_exit'])>1e-12:
            raise ValueError(f'Unmatched material {path}')
        by_offset={r['relative_exit']:r for r in rows if r['kind']=='geometry'}
        for reading in (r for r in rows if r['kind']=='signed_work'):
            h=reading['half_offset']
            lo,hi=by_offset[-h],by_offset[h]
            force=(lo['excess_energy']-hi['excess_energy'])/(2*h*shape.exit)
            component_error=sum(reading['components'].values())-force
            if max(abs(force-reading['force_on_exit']),abs(component_error))>1e-8*(1+abs(force)):
                raise ArithmeticError('Signed energy and separated work terms disagree')
            group=(spec['radius'],spec['mesh'],spec['polar_cells'],spec['outer'],spec['roof'],h,load)
            counts=groups.setdefault(group,{})
            key=spec['winding']
            if key in counts:
                raise ValueError(f'Duplicate {key}, {group}')
            counts[key]={'F':force,'nearby_mean_intake':.5*(lo['total_current']+hi['total_current']),
                         'sha256':fingerprint,'path':str(path)}
    for group,counts in sorted(groups.items()):
        if 0 not in counts:
            raise ValueError(f'Unpaired zero-winding state for {group}')
        radius,n,polar,outer,roof,h,load=group
        for m,case in sorted(counts.items()):
            delta=case['F']-counts[0]['F']
            emit(kind='matched_exit',radius=radius,mesh=n,polar=polar,outer=outer,roof=roof,h=h,
                 count=m,circulation=2*np.pi*m,zero_and_wound_sources=[counts[0]['sha256'],case['sha256']],
                 two_sided_carrier_force=case['F'],unwound_carrier_force=counts[0]['F'],
                 added_outward_force=delta,archived_closing_load=load,net_outward_force=delta-load,
                 fraction_closing_load_offset=delta/load,
                 nearby_intake_over_unwound=case['nearby_mean_intake']/counts[0]['nearby_mean_intake'])


def disk_readout(paths):
    for path in paths:
        rows,fingerprint=verified(path)
        case=rows[0]
        if case['mode']!='center':
            raise ValueError('Expected the centered disk field')
        samples=path/'fields.npz'
        with np.load(samples,allow_pickle=False) as archive:
            labels=list(archive['disk_columns'])
            columns={name:values for name,values in zip(labels,archive['disk'].T)}
            if case['material_sha256']!=str(archive['material_sha256']) or case['winding']!=archive['winding'].item():
                raise ValueError('Unmatched field and measurement')
        area=columns['area']
        j=columns['j_w']
        density=columns['density']
        weight=area*j
        current=sum(weight)
        geometry,=[r for r in rows if r['kind']=='geometry']
        if abs(current-geometry['total_current'])>1e-7 or min(density)<=0:
            raise ArithmeticError('Unresolved field current or density')
        normal=abs(j)/density**3
        merid=(columns['v_r']**2+columns['v_polar']**2+(j/density)**2)
        total=np.sqrt(merid+columns['v_phi']**2)/density**2
        edge=np.argmax(columns['s'])
        def measured_fraction(condition,measure):
            return float(measure[condition].sum()/measure.sum())
        angular_flux=area@(j*columns['s']*columns['v_phi'])
        emit(kind='disk_summary',radius=case['radius'],mesh=case['mesh'],polar=case['polar_cells'],
             winding=case['winding'],exit=case['reference_exit'],sha256=fingerprint,
             fields_sha256=hashlib.sha256(samples.read_bytes()).hexdigest(),current=current,
             angular_momentum_advected=angular_flux,angular_current_error=angular_flux-case['winding']*current,
             density_range=[density.min(),density.max()],closest_equatorial_rim_s=columns['s'][edge],
             equatorial_rim_density=density[edge],equatorial_rim_total_speed=np.sqrt(merid[edge]+columns['v_phi'][edge]**2),
             equatorial_rim_sound_speed=density[edge]**2,equatorial_rim_normal_Mach=normal[edge],
             equatorial_rim_total_Mach=total[edge],min_sampled_total_Mach=total.min(),
             total_supersonic_area_fraction=measured_fraction(total>=1.,area),
             normal_supersonic_area_fraction=measured_fraction(normal>=1.,area),
             outward_current_fraction_above_normal_sound=measured_fraction(normal>=1.,area*np.maximum(j,0)),
             negative_mass_flux_fraction=-np.minimum(weight,0).sum()/current)
        if abs(angular_flux-case['winding']*current)>2e-9*(1+abs(current)):
            raise ArithmeticError('Inconsistent angular transport across disk')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--disk',action='store_true')
    p.add_argument('guard_directories',type=Path,nargs='+')
    args=p.parse_args()
    (disk_readout if args.disk else compare)(args.guard_directories)
