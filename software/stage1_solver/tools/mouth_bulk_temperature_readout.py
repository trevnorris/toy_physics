"""Reduce independent zero-pressure temperature histories; keep signed and noisy readings visible."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import t as student

from neutral_patch_sites import emit


def guarded(folder):
    outcome,child,limits,duration,inv=[json.loads((folder/name).read_text()) for name in
        ('outcome.json','child-outcome.json','limit-validation.json','duration-validation.json','invocation.json')]
    assert outcome['limitsVerified'] and outcome['childOutcome']==child and child['guardReason'] is None
    assert outcome['exitCode']==child['exitCode']==0 and limits['verified'] and duration['verified']
    assert int(limits['actual']['memory.swap.max'])==0 and inv['seconds']==0 and duration['actual']['RuntimeMaxUSec']=='infinity'
    assert int(limits['actual']['memory.max'])==inv['memoryMax'] and not (folder/'stderr').read_bytes()
    return list(map(json.loads,(folder/'stdout').read_text().splitlines()))


def sample(folder):
    rows=guarded(folder)
    setup=rows[0]; frames=rows[1:]
    assert setup['test']=='setup' and len(frames)==49 and [r['frame'] for r in frames]==list(range(49))
    assert all(r['test']=='frame' for r in frames)
    ref=Path(setup['source'])
    assert hashlib.sha256(ref.read_bytes()).hexdigest()==setup['source_sha256']
    return setup,frames


def mean(values):
    data=np.array(values,dtype=float)
    return float(data.mean())


def interval(values):
    data=np.array(values,dtype=float)
    center=data.mean()
    radius=student.ppf(.975,len(data)-1)*data.std(ddof=1)/np.sqrt(len(data))
    return [float(center-radius),float(center+radius)]


def run(root):
    plan=json.loads((root/'plan.json').read_text())
    assert set(json.loads((root/'DONE.json').read_text())['completed'])==set(plan)
    energy_baseline={}
    for rung in range(5):
        cases=[]
        for history in range(3):
            groups={}
            for kind in ('prepare','undriven','pressure','shear','half_shear'):
                name=f'r{rung}h{history}_{kind}'
                setup,frames=sample(root/name)
                assert setup['rung']==rung and setup['history']==history
                groups[kind]=(setup,frames)
            _,un=groups['undriven']
            _,full=groups['shear']; _,half=groups['half_shear']; _,press=groups['pressure']
            for pair in (full,half):
                amp=groups['shear' if pair is full else 'half_shear'][0]['shears']
                for row in pair:
                    x,y=row['readings']
                    np.testing.assert_allclose((x['shear']-y['shear'])/(amp[0]-amp[1]),row['shear_slope'])
            for row in press:
                x,y=row['readings']
                applied=groups['pressure'][0]['applied_pressures']
                np.testing.assert_allclose(np.log(x['density']/y['density'])/(applied[0]-applied[1]),row['log_volume_compliance'],atol=1e-12)
            late_un=[r['readings'][0] for r in un[-12:]]
            late=mean([r['shear_slope'] for r in full[-12:]])
            halflate=mean([r['shear_slope'] for r in half[-12:]])
            early=full[0]['shear_slope']
            compliance=mean([r['log_volume_compliance'] for r in press[-12:]])
            temp=groups['undriven'][0]['temperature']
            energy=mean([r['energy_per_site'] for r in late_un])
            if rung==0:
                energy_baseline[history]=energy
            result=dict(early_shear=early,late_shear=late,half_late_shear=halflate,
                shear_retained=late/early,half_retained=halflate/half[0]['shear_slope'],
                full_minus_half=late-halflate,compliance=compliance,
                density=mean([r['density'] for r in late_un]),actual_pressure=mean([r['pressure'] for r in late_un]),
                compression_branch_pressure=[mean([r['readings'][i]['pressure'] for r in press[-12:]]) for i in (0,1)],
                initial_neighbors=un[0]['initial'],net_original_lost=un[-1]['original_missing'],
                ever_original_lost=un[-1]['ever_missing'],gained_neighbors=un[-1]['new'],
                neighbor_loss_fraction=un[-1]['original_missing']/un[0]['initial'] if un[0]['initial'] else None,
                largest_cluster=un[-1]['largest_cluster'],min_cluster_during_hold=min(r['largest_cluster'] for r in un),
                rms_nonaffine=un[-1]['rms_nonaffine'],initial_diffraction=un[0]['diffraction'],late_diffraction=un[-1]['diffraction'],
                density_quarters=[mean([r['readings'][0]['density'] for r in un[i+1:i+13]]) for i in (0,12,24,36)],
                potential_energy_per_site=energy,excess_potential_per_site=energy-energy_baseline[history])
            emit(test='history',rung=rung,history=history,temperature=temp,**result)
            cases.append(result)
        ensemble={}
        for name in ['early_shear','late_shear','half_late_shear','shear_retained','half_retained','full_minus_half','compliance','density','actual_pressure','excess_potential_per_site']:
            vals=[case[name] for case in cases]
            ensemble[name]=dict(mean=mean(vals),histories=vals,approximate_95=interval(vals))
        j=ensemble['compliance']
        bounds=j['approximate_95']
        k=1/j['mean'] if j['mean']>0 else None
        emit(test='rung',rung=rung,temperature=temp,**ensemble,dynamic_bulk_modulus=k,
             dynamic_bulk_interval=[1/bounds[1],1/bounds[0]] if bounds[0]>0 else None,
             ideal_density_times_T=temp*ensemble['density']['mean'],
             note='Intervals are approximate Student intervals across three independent heating histories; the finite-time zero-pressure trajectories do not establish an equilibrium phase.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('batch',type=Path)
    run(p.parse_args().batch)
