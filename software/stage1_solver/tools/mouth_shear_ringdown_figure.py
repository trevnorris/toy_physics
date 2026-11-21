"""Compare energy retained in shear-supporting material after the calibrated burst."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from mouth_order_wave_pulse import ROOT,emit
from mouth_shear_ringdown_compare import read


def series(tag, end=18.):
    setup,d=read(ROOT/f'shear_{tag}_01')
    t=d['time']-setup['drive_duration']
    work=d['work'][-1]
    active=(d['ordered']+d['upper']+d['lower'])/work
    tail=d['tail']/work
    feed=(d['loss']+d['feed'])/work
    exported=d.get('cut_export',np.zeros_like(t))/work
    receiver=d['receiver']/work
    err=max(abs((active+tail+feed+exported+receiver-1)[t>0]))
    if err>2e-8:
        raise ArithmeticError('Drive-off energy partitions do not sum')
    def average(x,mask): return float(np.trapz(x[mask],t[mask])/(t[mask][-1]-t[mask][0]))
    for lo,hi in ((.5,4.5),(4.5,8.5),(14,18)):
        m=(t>=lo)&(t<=hi)
        emit(test='fraction',tag=tag,window=[lo,hi],retained_active=average(active,m),
             weak_layer=average(tail,m),net_to_feed=average(feed,m),net_to_receiver=average(receiver,m),
             absorbed_at_cut=average(exported,m),input_work=work)
    # Downsample both curves to averages of the same 0.4 model-time bins.
    edges=np.arange(0.,min(end,t[-1])+1e-7,.4)
    centers=(edges[:-1]+edges[1:])/2
    envelopes=np.array([np.mean(active[(t>=a)&(t<b)]) for a,b in zip(edges[:-1],edges[1:])])
    return centers,envelopes


def main(out):
    colors={'constant':'#265d83','residual':'#bd662c','ordered':'#317b59'}
    fig,ax=plt.subplots(1,2,figsize=(9.3,3.7),layout='constrained',sharey=True)
    for tag,label,c in (('burst_one_f','Full transverse inertia','constant'),
                        ('burst_one_mobile_f','1/16 remains','residual'),
                        ('vanishing_one_f','Inertia proportional to shear order','ordered')):
        t,y=series(tag)
        ax[0].plot(t,100*y,label=label,color=colors[c],lw=2.)
    for tag,label,c,ls in (('cavity_full','Full transverse inertia',colors['constant'],'-'),
                            ('cavity_shallow','Order-proportional, shallow absorber',colors['ordered'],'-'),
                            ('cavity_deep','Order-proportional, deep absorber','#76558a','--'),
                            ('cavity_free','Order-proportional, free end','#3d5364',':')):
        t,y=series(tag,48.)
        ax[1].plot(t,100*y,label=label,color=c,ls=ls,lw=1.8)
    for tag in ('burst_half','burst_double','vanishing_half','vanishing_double',
                'vanishing_one','vanishing_one_free','vanishing_one_deep','vanishing_one_free_deep'):
        series(tag)
    for a,title in zip(ax,('Drive at the entrance','Drive inside the throat')):
        a.set_yscale('log')
        a.set_ylim(1e-2,105)
        a.set_xlim(0,18 if a is ax[0] else 48)
        a.set_title(title,loc='left',fontsize=11,fontweight='bold')
        a.grid(axis='y',alpha=.2)
        a.spines[['top','right']].set_visible(False)
        a.set_xlabel('Model time after drive stops')
        a.legend(loc='lower left',frameon=True,framealpha=.92,edgecolor='none',fontsize=8)
    ax[0].set_ylabel('Input work in throat with χ ≥ 0.1 (%)')
    fig.suptitle('Where the shear energy remains at radius 1',fontsize=14,fontweight='bold')
    fig.savefig(out,dpi=170)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    main(p.parse_args().output)
