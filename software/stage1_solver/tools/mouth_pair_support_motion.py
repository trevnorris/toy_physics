"""Let the surrounding support translate in the 0PN packet model's steady balance.

One passive scalar impedance opposes its lateral velocity relative to the
stationary distant reservoir. The complete aperture input/export is checked
using the separate packet export channel under a common Galilean shift.
"""
import argparse
from pathlib import Path
import json
import shlex

import numpy as np

from mouth_pair_neighbor_impulse import (MANIFEST, POLICIES, emit, feeding_packet, grid,
                                       samples)


def inspect(folder):
    def read(name):
        return json.loads((folder/name).read_text())
    if (read('child-outcome.json')['exitCode'] or read('outcome.json')['exitCode']
            or not read('limit-validation.json')['verified']):
        raise ValueError('Incomplete packet experiment')
    rows=[json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    records=json.loads(MANIFEST.read_text())['runs']
    for key in POLICIES:
        rays=[r for r in rows if r['test']=='ray' and r['policy']==key]
        labels=np.array([r['rate_label'] for r in rays])
        coefficients=np.array([r['responses'] for r in rays])
        for r in rows:
            if r['test']!='opening' or r['policy']!=key or r['method']!='cell' or r['feed_grid'][:2]!=[96,6]:
                continue
            rec,=[a for a in records if a['name']==f'L96C6_d{r["d"]}']
            args=shlex.split(rec['command'])
            target=Path(args[args.index('--baseline')+1])
            g=grid(*r['feed_grid'])
            data,_=feeding_packet(target,{target:r['input_pins'][str(target)]},g,[0])
            s,mu,vol,fields=samples(g,data,'cell')
            inward=fields[:,0]*vol
            _,idx=np.unique(s.round(12),return_inverse=True)
            group=np.bincount(idx,weights=inward)
            coord=((group.cumsum()-.5*group)/sum(group))[idx]
            weight=inward*(1-mu**2)
            table=np.column_stack([np.interp(coord,labels,col) for col in coefficients.T])
            C,Bu,Eu,minusG,Bv,Ev=weight @ table
            G=-minusG
            # Increasing W is equivalent to subtracting W from both drives
            # and restoring the lab velocity to every outgoing parcel.
            mouth_w=G-C
            brane_w=-(Bu+Bv)
            export_w=weight.sum()-Eu-Ev
            if not np.isclose(C,r['uniform_held_slope']) or not np.isclose(G,r['towing_drag_slope']):
                raise ArithmeticError('Mismatched fixed-support calibration')
            if min(C,G,Bu,Bv)<=0:
                raise ArithmeticError('Measured coefficients do not support positive entrainment bound')
            for mode, impedance in [('free',0.),('matched',Bu+Bv),('fixed',None)]:
                held_w=0. if impedance is None else r['brane']/(Bu+Bv+impedance)
                tow_w=0. if impedance is None else Bv/(Bu+Bv+impedance)
                held_rim=r['mouth']+mouth_w*held_w
                held_brane=r['brane']+brane_w*held_w
                held_out=r['export']+export_w*held_w
                drag=G-mouth_w*tow_w
                tow_brane=Bv+brane_w*tow_w
                tow_out=Ev+export_w*tow_w
                balance=held_rim+held_brane+held_out-r['carrier_total']
                tow_balance=-drag+tow_brane+tow_out
                emit(state=r['state'],policy=key,d=r['d'],mode=mode,
                     input_record=str(folder/'stdout'),C=C,G=G,Bu=Bu,Bv=Bv,
                     impedance=impedance,held_support_speed=held_w,held_mouth=held_rim,
                     held_brane=held_brane,held_export=held_out,tow_support_speed_per_V=tow_w,
                     drag=drag,tow_brane=tow_brane,tow_export=tow_out,held_balance=balance,
                     tow_balance=tow_balance,carrier_total=r['carrier_total'],
                     cancel_drag_W_over_V=G/(G-C) if G!=C else None)
                if (abs(balance)>2e-11 or abs(tow_balance)>3e-9 or
                    (impedance is not None and abs(held_brane-impedance*held_w)>2e-12)
                    or not 0<=tow_w<=1 or drag<min(C,G)-1e-12):
                    raise ArithmeticError('Floating support or outgoing momentum fails balance')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folders',type=Path,nargs='+')
    for f in p.parse_args().folders:
        inspect(f)
