"""Compare the two reference-transport experiments at the same held displacement."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mouth_order_wave_pulse import digest, emit


def rows(folder):
    status=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if (status['exitCode'] or status['childOutcome']['exitCode'] or not status['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Unverified job in {folder}')
    return [json.loads(text) for text in (folder/'stdout').read_text().splitlines()]


def select(data,name,**fields):
    return [r for r in data if r['test']==name and all(r[k]==v for k,v in fields.items())]


def main(args):
    result=rows(args.source)
    initial=select(result,'initial')[0]
    clock=initial['projected_clock_to_deep_end']
    carrier=np.array(initial['carrier'])
    site=select(result,'fixed_sites')
    drifting=select(result,'travelling_hold',lifetime_over_clock=1.)
    end=next(d for d in drifting if d['time_over_clock']==1.)
    for name,t,net in [('initial',0,carrier+initial['frozen']),
                       ('local',clock,carrier+select(result,'fixed_sites',time_over_lifetime=1.)[0]['shear']),
                       ('passing',clock,carrier+end['shear'])]:
        emit(test='deep_bend',model=name,time=t,offset=.01,holder=(net@np.array([0.,.01])).tolist())
    cross=select(result,'site_check')[0]['first_zero_over_lifetime']
    emit(test='local_first_zero',transits=cross,model_time=cross*clock)
    # Same mesh, rates changed in opposite directions. Only the reference
    # advection speed changes; the carrier keeps the same intake.
    full,slow=rows(args.matched),rows(args.half_speed)
    for stress_free in [False,True]:
        def plateau(data,lifetime):
            return np.array(select(data,'inlet_comparison',lifetime_over_clock=lifetime,
                                   stress_free_entry=stress_free)[0]['net'])
        a,b=plateau(full,1.),plateau(slow,2.)
        emit(test='equal_transport_distance',stress_free_entry=stress_free,
             maximum_dc_difference=float(abs(a-b).max()))
        if abs(a-b).max()>1e-9:
            raise ArithmeticError('Equal contact transport lengths changed the steady result')
    # Positive values support every direction in this two-dimensional test.
    fig,ax=plt.subplots(figsize=(7.1,4.2),layout='constrained')
    x1=np.array([v['time_over_lifetime'] for v in site])
    y1=np.array([v['net_eigs'][0] for v in site])
    x2=np.r_[0.,[v['time_over_clock'] for v in drifting]]
    y2=np.r_[initial['net_eigs'][0],[v['net_eigs'][0] for v in drifting]]
    ax.axhline(0,color='#565656',lw=1.1)
    ax.axvline(1,color='#818181',ls=':',lw=1.2)
    for x,y,label,color,mark in [(x1,y1,'Contacts renew in place','#9a4c25','o'),
                                  (x2,y2,'Contacts travel with the material','#186b71','s')]:
        use=x<=4
        ax.plot(x[use],y[use],label=label,color=color,marker=mark,lw=1.65,ms=5)
    ax.scatter([cross],[0],marker='x',s=55,color='#9a4c25',zorder=4)
    ax.set(xlim=(-.07,4.08),ylim=(-.014,.009),
           xlabel='Time held / carrier transit from waist to end of bend',
           ylabel='Least restoring stiffness (model units)',
           title='Holding the bend through material replacement')
    ax.spines[['top','right']].set_visible(False)
    ax.legend(loc='lower left',frameon=False)
    fig.savefig(args.figure,dpi=180)
    plt.close(fig)
    emit(test='figure',path=str(args.figure),sha256=digest(args.figure))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','matched','half-speed','figure']:
        p.add_argument('--'+name,required=True,type=Path)
    main(p.parse_args())
