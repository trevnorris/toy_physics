"""Held-end shape test for local reversible flow; use both signed fields."""
import argparse
import json
from mouth_exchange_carrier import assemble
from mouth_exchange_interior import InteriorShape

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--mesh',type=int,default=80)
    p.add_argument('--sheet',type=float,default=.04/6)
    args=p.parse_args()
    step=.005
    rows=[assemble(t,args.mesh,args.sheet,28.,22.,geometry=InteriorShape)
          for t in (-step,step)]
    slope=(rows[1]['theta_energy']-rows[0]['theta_energy'])/(2*step)
    mean_work=sum(r['work'] for r in rows)/2
    print(json.dumps(dict(states=rows,interior_force_slope=slope,average_virtual_work=mean_work)))
    if abs(slope-mean_work)>2e-5:
        raise ArithmeticError('Held-end work disagrees with the energy difference')
