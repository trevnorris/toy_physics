"""Compare the allowed winding sectors on the two held mouths."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[3]
record=root/'research/pde_ledger_v3/_measurements/mouth_static_circulation.json'


def plot():
    raw=json.loads(record.read_text())['runs']['circulation_gate01/force']
    if raw['child_exit'] or raw['guard_exit']:
        raise ValueError('Unverified winding comparison')
    rows=[r for r in map(json.loads,raw['stdout'].splitlines()) if r['mesh']==13 and r['h']==.006]
    fig,(force,flow)=plt.subplots(1,2,figsize=(8.7,3.8),layout='constrained')
    for radius,offset,color,marker in [(1,-.055,'#195c82','o'),(2,.055,'#ac4e25','s')]:
        data=sorted((r for r in rows if r['radius']==radius),key=lambda r:r['count'])
        if [r['count'] for r in data]!=[0,1,2]:
            raise ValueError('Missing winding sector')
        # Horizontal offset keeps the values at each integer visible.
        x=[r['count']+offset for r in data]
        force.scatter(x,[r['net_outward_force'] for r in data],color=color,marker=marker,s=54,label=f'Radius {radius}')
        flow.scatter(x,[100*r['nearby_intake_over_unwound'] for r in data],color=color,marker=marker,s=54)
    force.axhline(0,color='.35',linewidth=1,linestyle='--')
    force.set(title='Force on the held exit',ylabel='Net outward force (solver units)',ylim=(-.46,.48))
    force.legend(loc='upper left',frameon=False)
    flow.set(title='Inward flow at the same phase bias',ylabel='Intake relative to unwound mouth (%)',
             yscale='log',ylim=(.0001,200),yticks=[.0001,.001,.01,.1,1,10,100],
             yticklabels=['0.0001','0.001','0.01','0.1','1','10','100'])
    for panel in (force,flow):
        panel.set(xlim=(-.35,2.35),xticks=[0,1,2],xlabel='Winding quanta')
        panel.spines[['top','right']].set_visible(False)
    flow.grid(axis='y',alpha=.2)
    fig.savefig(root/'software/stage1_solver/reports/figures/mouth_static_circulation.png',dpi=160)
    plt.close(fig)


if __name__=='__main__':
    plot()
