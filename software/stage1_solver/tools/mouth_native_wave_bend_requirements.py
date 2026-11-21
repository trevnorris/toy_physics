"""Turn the two measured shape curvatures into a wave-empty wall requirement."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh
from scipy.optimize import brentq

from mouth_carrier_attachment_tension import area_stiffness
from mouth_flared_receiver import FlaredDomain
from mouth_order_wave_pulse import digest, emit


def measured(path):
    def load(name): return json.loads((path/name).read_text())
    inv,out,limits=load('invocation.json'),load('outcome.json'),load('limit-validation.json')
    if not (out['limitsVerified'] and limits['verified'] and out['exitCode']==out['childOutcome']['exitCode']==0
            and inv['seconds']==inv['swapMax']==0):
        raise ValueError(f'Unverified run {path}')
    return [json.loads(q) for q in (path/'stdout').read_text().splitlines()]


def run(a):
    dom=FlaredDomain(20,'cubic',15.,far=14.)
    areas=[area_stiffness(dom,interval) for interval in ((.55,1.15),(1.2,2.0))]
    wall=np.diag([k for k,_ in areas])
    emit(test='material_wall',area_stiffness_per_unit_tension=wall.diagonal().tolist(),
         finite_offset_area_checks=[fd for _,fd in areas],radius=1.)
    for folder in a.families:
        rows=measured(folder)
        for row in rows[1:]:
            raw=np.array(row['mean_work_stiffness'])
            sym=(raw+raw.T)/2
            def directional(t):
                v=np.array([t,1-t])
                return v@sym@v
            thresholds=eigvalsh(-sym,wall)
            critical=max(thresholds[-1],0.)
            smallest=np.linalg.eigvalsh(sym+critical*wall)[0]
            if abs(smallest)>2e-12: raise ArithmeticError('Area-tension threshold failed')
            emit(test='wave_requirement',run=folder.name,degree=row['degree'],
                 collar_stiffness=sym[0,0],deep_stiffness=sym[1,1],
                 stiffness_at_five_percent_collar=directional(.05),
                 monotone_zero=brentq(directional,0,1) if directional(0)*directional(1)<0 else None,
                 redistribution_at_fixed_total=np.array([1.,-1.])@sym@np.array([1.,-1.]),
                 marginal_wall_tension_per_stored_energy=critical)
    if a.figure:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,(left,right)=plt.subplots(1,2,figsize=(9.8,4.0),layout='constrained')
        rows=measured(a.offsets)
        cases=[v for v in rows if v['test']=='displaced' and v['interval']==[.55,1.15]]
        coeff=next(v['virtual_work_over_energy'] for v in rows if v['test']=='bend')
        x=np.linspace(-.021,.021,100)
        left.plot(x,-coeff*x,color='#873b21',label='Small-offset prediction')
        left.scatter([v['offset'] for v in cases],[v['flare_force_over_energy'] for v in cases],
                     color='#252525',s=32,zorder=3,label='Independently solved offsets')
        left.set(title='Degree four, bend near entrance',xlabel='Deeper end offset',ylabel='Force on deeper end / stored energy')
        left.ticklabel_format(axis='y',style='sci',scilimits=(-3,-3))
        left.legend(frameon=False,fontsize=9)
        t=np.linspace(0,.05,150)
        vec=np.column_stack((t,1-t))
        for folder,label,color in zip(a.families[:3],('Coarse','Fine','Finest'),('#aaa196','#237773','#202020'),strict=True):
            row=measured(folder)[1]
            h=np.array(row['mean_work_stiffness'])
            right.plot(t*100,-np.einsum('ni,ij,nj->n',vec,h,vec),label=label,color=color)
        right.set(title='Degree one, bend mostly at depth',xlabel='Share of total offset near entrance (%)',
                  ylabel='Force / (stored energy × offset)')
        right.ticklabel_format(axis='y',style='sci',scilimits=(-4,-4))
        right.legend(frameon=False,fontsize=9)
        for panel in (left,right):
            panel.axhline(0,lw=.8,color='#555555')
            panel.spines[['right','top']].set_visible(False)
        fig.savefig(a.figure,dpi=160)
        plt.close(fig)
        emit(test='figure',path=str(a.figure),sha256=digest(a.figure))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('families',nargs='+',type=Path)
    p.add_argument('--figure',type=Path)
    p.add_argument('--offsets',type=Path)
    run(p.parse_args())
