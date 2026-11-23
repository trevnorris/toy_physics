"""Plot the crystalline order on the initial and relaxed meridional sections at entrance four."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_sealed_shape_relax import Meridian


def read_state(path):
    guard=json.loads((path.parent/'outcome.json').read_bytes())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('Unverified material endpoint')
    with np.load(path,allow_pickle=False) as z:
        if z['radius']!=4:
            raise ValueError('Expected entrance four')
        extra=float(z['extra']) if 'extra' in z else float(z.get('shape_extra',0.))
        shape=Meridian(4.,z['shape_params'],extra)
        chi=z['order']
    return shape,chi


def plot(initial, endpoint, target):
    fig,axes=plt.subplots(2,1,figsize=(7.7,4.45),sharex=True,sharey=True,layout='constrained')
    for ax,(shape,chi),title in zip(axes,(read_state(initial),read_state(endpoint)),
             ('Initial preparation','After joint shape and order relaxation')):
        t=np.linspace(0,1,chi.shape[0])
        radial=np.linspace(0,1,chi.shape[1])[None,:]*shape.profile(t)[0][:,None]
        depth=np.broadcast_to((shape.length*t)[:,None],chi.shape)
        color=ax.pcolormesh(radial,depth,chi,vmin=0,vmax=1,cmap='viridis',shading='gouraud',rasterized=True)
        ax.plot(radial[:,-1],depth[:,-1],color='black',lw=1.25)
        ax.plot(radial[-1],depth[-1],color='black',lw=1.25)
        ax.set(title=title,ylabel='Depth w',aspect='equal',ylim=(1.08,-.025),xlim=(0,4.18))
        ax.set_title(title,loc='left',fontsize=11)
        ax.set_title('',loc='center')
        ax.spines[['right','top']].set_visible(False)
    axes[-1].set_xlabel('Distance from the mouth center in the brane')
    bar=fig.colorbar(color,ax=axes,shrink=.70,ticks=[0,.5,1])
    bar.set_label('Crystalline order')
    fig.savefig(target,dpi=165)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('initial','endpoint','figure'):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    plot(a.initial,a.endpoint,a.figure)
