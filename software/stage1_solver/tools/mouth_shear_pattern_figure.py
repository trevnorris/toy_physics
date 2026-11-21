"""Plot measured energy after release and its dependence on the transverse inertia."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_order_wave_pulse import digest, emit
from mouth_shear_ringdown_compare import read


def main(root, out):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.4, 4.1), constrained_layout=True)
    colors = ('#146b68', '#88421f')
    for folder, degree, branch, ls in (('l1b1',1,1,'-'), ('l1b2',1,2,'--'),
                                       ('l4b1',4,1,'-'), ('l2b2',2,2,'--')):
        path = root/f'pattern_{folder}_01'
        final = json.loads((path/'stdout').read_text().splitlines()[-1])
        if json.loads((path/'outcome.json').read_text())['exitCode'] or digest(path/'trace.npz') != final['trace_sha256']:
            raise ValueError('Unverified modal release')
        with np.load(path/'trace.npz') as f:
            data = dict(zip(f['columns'], f['trace'].T, strict=True))
        ax.plot(data['time'], 100*data['active_throat'], ls, color=colors[branch-1],
                alpha=1 if degree==1 else .71, lw=2,
                label=f'Degree {degree}, '+('lowest' if branch==1 else 'first overtone'))
    ax.set(xlabel='Time after release', ylabel='Initial energy in ordered throat (%)',
           title='Which pattern stays near the mouth?', xlim=(0,48), ylim=(-3,104))
    ax.legend(frameon=False, fontsize=8, loc='center right')
    vals = []
    for name in ('drive_fund_fine', 'material_residual16_fine', 'material_constant_fine'):
        s, d = read(root/f'pattern_{name}_01')
        t = d['time']-s['drive_duration']
        window = (t>=44)&(t<48)
        w = d['work'][-1]
        vals.append([np.mean(e[window])/w*100 for e in (
            d['ordered']+d['upper']+d['lower'],d['tail'],d['feed']+d['loss'],d['cut_export']+d['receiver'])])
    base = np.zeros(3)
    for n, (color,label) in enumerate((('#146b68','Ordered throat'),('#d1a344','Weak-shear throat'),
                                       ('#4369a4','Feeding brane'),('#ac6876','Cut or receiver'))):
        v = np.array(vals)[:,n]
        bx.bar(np.arange(3),v,bottom=base,color=color,width=.66,label=label)
        base += v
    np.testing.assert_allclose(base,100.,atol=1e-7)
    bx.set(title='Same localized torque at 3.184', ylabel='Injected work at times 44–48 (%)',
           xlabel='Inertia remaining at zero shear order', ylim=(0,105),
           xticks=[0,1,2],xticklabels=['Zero','1/16','Full'])
    bx.legend(loc='upper center', bbox_to_anchor=(.5,-.20),ncol=2,frameon=False,fontsize=8)
    for panel in (ax,bx):
        panel.spines[['right','top']].set_visible(False)
    fig.savefig(out,dpi=160)
    plt.close(fig)
    emit(figure=str(out),sha256=digest(out),percent_by_material=vals)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('_scratch/s11c'))
    p.add_argument('output',type=Path)
    a=p.parse_args()
    main(a.root,a.output)
