"""Join held-exit shape work and disk readouts by the physical exit coordinate."""
import argparse
import json
from collections import defaultdict
from pathlib import Path

from mouth_poly_shallow_release import emit


def read_guarded(path):
    guard = json.loads((path.parent/'outcome.json').read_bytes())
    limits = json.loads((path.parent/'limit-validation.json').read_bytes())
    if (guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Unverified measurement: {path}')
    return [json.loads(line) for line in path.read_bytes().splitlines()]


def collect(root):
    carriers, unwound = {}, {}
    materials = defaultdict(list)
    for f in sorted((root/'carrier').glob('*/stdout')):
        rows = read_guarded(f)
        inp = rows[0]
        key = (inp['radius'], inp['exit'], inp['mesh'])
        disk, = [r for r in rows if r['kind'] == 'disk']
        if inp['winding'] == 0:
            unwound[key] = disk['current']
            continue
        force, = [r for r in rows if r['kind'] == 'force']
        if key in carriers or inp['winding'] != 1 or not (disk['current']>0):
            raise ValueError(f'Unexpected carrier: {key}')
        carriers[key] = (force['carrier_outward'], disk, str(f.relative_to(root)), inp['material_sha256'])
    for f in sorted((root/'material').glob('*/stdout')):
        states = {}
        for r in read_guarded(f):
            if r['kind'] == 'material_input':
                states[r['radius'],r['exit']] = r['order_sha256']
            if r['kind'] == 'material_force':
                key = (r['radius'],r['exit'])
                materials[key].append((r['step'], r['absolute_material_force'],
                                       str(f.relative_to(root)), states[key]))
    for f in sorted((root/'review').glob('*/stdout')):
        rows = read_guarded(f)
        if rows[0]['kind'] != 'input' or 'relative_signed_step' not in rows[0]:
            continue
        for r in rows:
            if r['kind'] == 'audit_force':
                materials[r['radius'],r['exit']].append((r['half_offset'],r['material'],
                                                        str(f.relative_to(root)),rows[0]['material_sha256']))
    results=[]
    for (radius,b,mesh),(carrier,disk,car_source,digest) in sorted(carriers.items()):
        step,material,mat_source,mat_digest = min(materials[radius,b])
        if digest!=mat_digest or abs(disk['angular_momentum_out']-disk['current'])>1e-7:
            raise ValueError('Preparation or transported angular momentum mismatch')
        results.append(dict(radius=radius,exit=b,mesh=mesh,force=carrier+material,
                            carrier_force=carrier,material_force=material,material_half_step=step,
                            material_source=mat_source,carrier_source=car_source,order_sha256=digest,
                            unwound_current=unwound.get((radius,b,mesh)),disk=disk))
    return results


def figure(rows, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,(left,right)=plt.subplots(1,2,figsize=(9.6,4.2),layout='constrained')
    for radius,color in [(1,'#a44825'),(2,'#225d83')]:
        data=[r for r in rows if r['radius']==radius and r['mesh']==13]
        b=[r['exit'] for r in data]
        left.plot(b,[r['force'] for r in data], 'o-', color=color, ms=4, lw=1.,label=f'Entrance {radius}')
        if radius==2:
            for prop,label,style in [('normal_subsonic','Normal speed below sound','o-'),
                         ('dense_and_total_subsonic','Dense, total speed below sound','s--')]:
                right.plot(b,[100*r['disk'][prop]['inward_flux'] for r in data], style,
                           ms=4,lw=1.2,label=label)
            right.axvspan(1.040,1.042,color='.2',alpha=.6)
    left.set(yscale='symlog',ylabel='Net outward force (solver units)',title='Force on a held exit')
    left.axhline(0,color='0.35',lw=.75)
    right.set(ylabel='Share of inward current (%)',ylim=(-4,100),title='Flow through the radius-2 mouth')
    right.text(1.09,71,'Force crossing',fontsize=9)
    for ax in (left,right):
        ax.set_xlabel('Exit radius (solver units)')
        ax.spines[['top','right']].set_visible(False)
        ax.legend(frameon=False,fontsize=9)
    fig.savefig(out,dpi=155)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('data',type=Path)
    p.add_argument('--figure',type=Path)
    args=p.parse_args()
    if args.data.is_dir():
        data=collect(args.data)
        for r in data:
            emit(kind='comparison',**r)
    else:
        data=read_guarded(args.data)
        if not data or any(r.get('kind')!='comparison' for r in data):
            raise ValueError('Expected saved comparison rows')
    if args.figure:
        figure(data,args.figure)
