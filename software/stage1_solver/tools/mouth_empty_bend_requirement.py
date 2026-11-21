"""Read the carrier curvature and express the independently missing wall stiffness."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

from mouth_order_wave_pulse import emit,digest


def results(p):
    read=lambda name:json.loads((p/name).read_text())
    inv,guard,out=read('invocation.json'),read('limit-validation.json'),read('outcome.json')
    if not (out['limitsVerified'] and guard['verified'] and out['exitCode']==out['childOutcome']['exitCode']==0
            and inv['swapMax']==inv['seconds']==0):
        raise ValueError('Unverified job')
    return [json.loads(x) for x in (p/'stdout').read_text().splitlines()]


def main(a):
    archive=json.loads(a.wave_record.read_text())['runs']
    def prior(label):
        row=archive['_scratch/s11c/native_bend_'+label+'_01']
        if not row['verified_limits'] or row['worker_exit'] or row['guarded_exit']:
            raise ValueError('Unverified reference')
        return [json.loads(x) for x in row['literal_stdout']]
    wall=np.diag(prior('requirements')[0]['area_stiffness_per_unit_tension'])
    waves=[row for row in prior('family_third') if row['test']=='pattern']
    base_dir=Path('_scratch/s11c')
    for name in ('c26','c34','c48','alt34','larger'):
        carrier=results(base_dir/f'empty_bend_{name}_01')[-1]
        k=np.array(carrier['holding_stiffness'])
        threshold=eigvalsh(-k,wall)[-1]
        n0=np.linalg.eigvalsh(k)
        twice=np.linalg.eigvalsh(k+2*threshold*wall)
        emit(test='empty_threshold',run=name,marginal_sigma=threshold,unloaded_eigenvalues=n0.tolist(),
             twice_threshold_eigenvalues=twice.tolist(),area_matrix=wall.tolist())
        if not (n0[-1]<0 and twice[0]>0 and abs(np.linalg.eigvalsh(k+threshold*wall)[0])<1e-12):
            raise ArithmeticError('Independent tension threshold check failed')
    carrier=results(base_dir/'empty_bend_c48_01')[-1]
    k=np.array(carrier['holding_stiffness'])
    for row in waves:
        l=row['degree']; wave=np.array(row['mean_work_stiffness'])
        sym=(wave+wave.T)/2
        emit(test='occupied_requirement',degree=l,carrier=k.tolist(),wave_per_energy=sym.tolist(),
             empty_sigma=eigvalsh(-k,wall)[-1],wave_only_sigma_per_energy=eigvalsh(-sym,wall)[-1])
    if a.figure:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6.6,4.1),layout='constrained')
        for name,j,label,color in [('fdnear3',0,'Bend near entrance','#276c72'),('fddeep3',1,'Bend deeper in flare','#8f4d27')]:
            data=[v for v in results(base_dir/f'empty_bend_{name}_01') if v['test']=='finite_offset']
            x=np.array([v['offset'][j] for v in data])
            y=-np.array([v['holding_force'][j] for v in data])
            ls=np.linspace(-.021,.021,30)
            ax.plot(ls,-k[j,j]*ls,color=color,lw=1.6)
            ax.scatter(x,y,label=label,color=color,s=28)
        ax.axhline(0,color='#888',lw=.7); ax.axvline(0,color='#888',lw=.7)
        ax.set(xlabel='Signed offset of the bend (model length)',ylabel='Force exerted by flowing material (model units)',
               title='With the wave absent, the flow pushes each bend farther')
        ax.ticklabel_format(style='sci',axis='y',scilimits=(-4,-4))
        ax.legend(frameon=False)
        ax.spines[['top','right']].set_visible(False)
        fig.savefig(a.figure,dpi=160)
        plt.close(fig)
        emit(test='figure',path=str(a.figure),sha256=digest(a.figure))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wave-record',type=Path,default=Path('research/pde_ledger_v3/_measurements/mouth_native_wave_bend.json'))
    p.add_argument('--figure',type=Path)
    main(p.parse_args())
