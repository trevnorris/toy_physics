"""Plot the radius-two joint-radiation gate from its verified measurement record."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).resolve().parents[3]
manifest=root/'research/pde_ledger_v3/_measurements/mouth_curved_exit_gate.json'
runs=json.loads(manifest.read_text())['runs']


def readings(key):
    data=runs[key]
    if data['child_exit'] or data['guard_exit']:
        raise ValueError(f'Unvalidated result {key}')
    return list(map(json.loads,data['literal_stdout'].splitlines()))


def draw():
    profile,=readings('curved_r2_sound01/flux22')
    edges=np.r_[0.,np.cumsum(profile['area'])]
    edges=(edges/edges[-1])**(1/3)
    meshes=[]
    for key in runs:
        if not key.startswith('mouth_curved_exit02/case_'):
            continue
        setup,pole,_=readings(key)
        args=setup['args']
        if (args['front'],args['box'],args['offset'])==('mean',6.,0.):
            meshes.append((args['mesh'],100*pole['pole_cycle_loss']))
    meshes.sort()
    fig,(left,right)=plt.subplots(1,2,figsize=(9,3.9),layout='constrained')
    for field,color,label in [('through_plus','#225e89','Fixed reference disk'),
                              ('relative_plus','#aa4c28','Relative to moving disk')]:
        complex_pairs=np.array(profile['bins'][field])
        magnitude=np.hypot(complex_pairs[:,0],complex_pairs[:,1])
        left.stairs(magnitude,edges,baseline=None,lw=2,color=color,label=label)
    left.set(xlabel='Radial position / exit radius',ylabel='Oscillatory current density (mode units)',
             title='Where the disk current changes',xlim=(0,1),yscale='log')
    left.legend(frameon=False,loc='upper left',fontsize=9)
    n,loss=zip(*meshes)
    right.plot(n,loss,marker='o',lw=2,color='#316987')
    right.axhline(1,color='#a74526',lw=1.6,ls='--',label='1% release threshold')
    right.set(xticks=n,xlim=(10.5,23.5),ylim=(0,17),xlabel='Acoustic cells across disk radius',
              ylabel='Breathing energy lost per cycle (%)',title='Loss under mesh refinement')
    for i,value in zip(n,loss):
        right.text(i,value+.28,f'{value:.2f}%',ha='center',fontsize=9)
    right.legend(frameon=False,loc='lower left',bbox_to_anchor=(.02,.09))
    fig.suptitle('Radius 2  ·  conversion front calibrated at the mean current',fontsize=13)
    for ax in (left,right):
        ax.spines[['right','top']].set_visible(False)
    destination=root/'software/stage1_solver/reports/figures/mouth_curved_exit_gate.png'
    fig.savefig(destination,dpi=160)
    plt.close(fig)


if __name__=='__main__':
    draw()
