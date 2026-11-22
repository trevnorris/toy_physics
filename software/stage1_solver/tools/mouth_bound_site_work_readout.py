"""Read the conditional physical work and a lower bound on slip heat at full crossings."""
import argparse
import json
from pathlib import Path

from neutral_patch_sites import emit


def read(run):
    guard=json.loads((run/'outcome.json').read_text())
    child=json.loads((run/'child-outcome.json').read_text())
    if guard['childOutcome']!=child or guard['exitCode'] or child['exitCode'] or child['guardReason'] or not guard['limitsVerified']:
        raise ValueError('Unverified release trajectory')
    rows=[json.loads(line) for line in (run/'stdout').read_text().splitlines()]
    path=next(r for r in rows if r['test']=='path')
    for mid in (r for r in rows if r['test']=='midpoint'):
        for r in rows:
            if r.get('outcome')!='passed_half_order' or r['law']!=mid['law'] or r['gamma']<=0:
                continue
            loss=r['integrated_specific_dissipation']  # gamma * integral v^2 dt
            work=r['gamma']*mid['phase'][-1]-loss
            independent=mid['mu'][-1]+.5*(r['speed']**2-path['entry_speed']**2)
            # Cauchy-Schwarz on sqrt(gamma)*(u-v) and sqrt(gamma)*v.
            # A passive bath losing energy to the site also supplies the slip heat.
            heat_floor=work**2/loss
            emit(test='crossing_work',law=r['law'],gamma=r['gamma'],mesh=path['mesh'],
                 segments=path['segments'],host_side=path['host_side'],site_work=work,
                 endpoint_site_work=independent,slip_heat_lower_bound=heat_floor,
                 net_bath_energy_debit_lower_bound=work+heat_floor)
            if abs(work-independent)>max(5e-10,2*r['energy_ledger_max_error']):
                raise ArithmeticError('Work mismatch')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs',type=Path,nargs='+')
    for folder in p.parse_args().runs:
        read(folder)
