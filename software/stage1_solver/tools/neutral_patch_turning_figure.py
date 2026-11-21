"""Plot internal turning stiffness and the separate whole-patch force trial."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from neutral_patch_resting_order_hold import digest, emit


def records(name):
    folder = Path('_scratch/s11c')/name
    result = json.loads((folder/'outcome.json').read_text())
    limits = json.loads((folder/'limit-validation.json').read_text())
    if result['exitCode'] or result['childOutcome']['exitCode'] or not limits['verified']:
        raise ValueError('Input run incomplete')
    return [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]


def main():
    torques = records('patch_neighbor_torque_v1')
    forces = records('turning_patch_force01')
    fig, (left,right) = plt.subplots(1,2,figsize=(9.2,4.4))
    orders = [.95, .5, .1, .01]
    for ax in (left,right):
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',color='#d8dde2',lw=.7)
    for law, color in [('central','#1e5a88'),('layered','#a64627')]:
        for meter, style, label in [('pooled','-','joint turn'),('within_plane','--','within each layer')]:
            group = {r['chi']:r['stiffness'] for r in torques if r['test']=='neighbor_torque'
                     and r['law']==law and r['inner_per_plane']==4 and r['meter']==meter}
            left.plot(range(4),[group[c]/group[.95] for c in orders],marker='o',ls=style,
                      color=color,label=f'{law.capitalize()}, {label}',lw=2)
    left.set_xticks(range(4), ['95%','50%','10%','1%'])
    left.set_xlabel('Local shear order')
    left.set_yscale('log')
    left.set_ylabel('Stiffness / value at 95% order')
    left.set_title('Turning sites against nearby sites',fontsize=11)
    left.legend(frameon=False,fontsize=8,loc='lower left')
    x = np.arange(2)
    names = ['central','layered']
    entry = {r['patch']:r for r in forces if r['test']=='patch_reference' and r['order']==.5}
    finite = {(r['patch'],r['sheet']):r['material_force'][2] for r in forces if r['test']=='bonded_material'
              and r['order']==.5 and r['ratio']==8}
    values = [[entry[law]['detached_force'] for law in names],
              [finite[law,'curl'] for law in names], [finite[law,'symmetric'] for law in names]]
    for i,(series,tag,color) in enumerate(zip(values,['Free edge','Curl-elastic sheet','Strain-elastic sheet'],
                                               ['#8a94a0','#277491','#c57348'],strict=True)):
        right.bar(x+(i-1)*.24,series,width=.21,label=tag,color=color,zorder=2)
    right.set_xticks(x,['Central bonds','Layered bonds'])
    right.set_ylabel('Interlayer force / imposed strain')
    right.set_title('Half order, after one transit',fontsize=11)
    right.legend(frameon=False,fontsize=8)
    fig.tight_layout()
    target = Path('software/stage1_solver/reports/figures/neutral_patch_neighbor_turn.png')
    fig.savefig(target,dpi=160)
    plt.close(fig)
    emit('figure',path=str(target),sha256=digest(target))


if __name__=='__main__':
    main()
