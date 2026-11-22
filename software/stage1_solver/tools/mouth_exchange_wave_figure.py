"""Plot the tested shape and the measured outgoing wave loss on the five signed shapes."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_exchange_geometry import ExchangeShape, reference

ROOT = Path(__file__).resolve().parents[3]
MEASUREMENTS = ROOT/'research/pde_ledger_v3/_measurements'


def main():
    wave = json.loads((MEASUREMENTS/'mouth_exchange_shear.json').read_text())
    run = wave['runs']['n24_p3_f0']
    raw = run['literal_stdout']
    if (hashlib.sha256(raw.encode()).hexdigest()!=run['stdout_sha256'] or
            run['worker_exit'] or run['guard_exit'] or not run['limits_verified']):
        raise ValueError('Unverified wave readout')
    modes = sorted((r for r in map(json.loads,raw.splitlines()) if r['test']=='mode'),
                   key=lambda r:r['step'])
    zero, = [r for r in modes if r['step']==0.]
    shape = ExchangeShape(0.)
    arc = np.linspace(0.,shape.length,300)
    radius,_ = shape.radius(arc)
    loss_percent = [-100*np.expm1(-2*r['amplitude_loss']*400*np.pi) for r in modes]
    figure,(geometry,leak) = plt.subplots(1,2,figsize=(9.4,3.65),layout='constrained')
    geometry.fill_between(arc,-radius,radius,color='#e5eef3')
    geometry.plot(arc,radius,color='#294a63',lw=2)
    geometry.plot(arc,-radius,color='#294a63',lw=2)
    for w in (.12*shape.length,.88*shape.length):
        b=float(shape.radius(w)[0])
        direction=np.sign(shape.displacement(w))
        for sign in (-1,1):
            geometry.annotate('',xy=(w,sign*(b+.16*direction)),
                              xytext=(w,sign*(b-.16*direction)),
                              arrowprops=dict(arrowstyle='-|>',color='#b04421',lw=1.8))
    geometry.set(xlim=(-.09,shape.length+.09),ylim=(-1.9,1.9),
                 xlabel='Depth into the bulk (waist radii)',ylabel='Radial section (waist radii)',
                 title='Direction of the positive displacement')
    geometry.text(.03,.02,'Arrows show direction only.',transform=geometry.transAxes,fontsize=9)
    leak.plot([100*r['step'] for r in modes],loss_percent,'o-',color='#276892')
    leak.set(xlabel='Nominal neck displacement (% of waist radius)',
             ylabel='Wave energy escaping into the brane (%)',
             title='Leakage over one reference slow cycle',ylim=(0.,1.18*max(loss_percent)))
    leak.set_xticks([-1,-.5,0,.5,1])
    leak.grid(axis='y',alpha=.25)
    for ax in (geometry,leak):
        ax.spines[['right','top']].set_visible(False)
    path = ROOT/'software/stage1_solver/reports/figures/mouth_exchange_wave.png'
    figure.savefig(path,dpi=155)
    bare = reference()[0]
    second = wave['summary']['n24_p3_f0']['second_per_action']
    energy = zero['frequency']*(-bare['area_stiffness'])/second
    print(json.dumps(dict(stored_energy_at_wave_only_threshold_over_sigmaR3=energy,
                          energy_ratio_to_lateral_surface_energy=energy/bare['area'],
                          energy_density_over_sigma_per_radius=energy/bare['volume'],
                          image_sha256=hashlib.sha256(path.read_bytes()).hexdigest())))


if __name__=='__main__':
    main()
