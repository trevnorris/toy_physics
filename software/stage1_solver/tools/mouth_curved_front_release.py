"""Initial velocity after release of a locally held converting spherical front."""
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mouth_curved_free_front_head import solve
from mouth_free_front_wave_head import emit


@dataclass
class LocalHolder:
    center: float
    strength: float
    width: float = 5.

    def affinity(self, radius):
        # An external energy density -lambda*b(r)*chi, zero outside the patch.
        t = np.asarray((radius-self.center)/self.width)
        active = t*t < 1.
        shape = np.zeros_like(t)
        shape[active] = np.exp(-t[active]**2/(1-t[active]**2))
        return -self.strength*shape


def main():
    record = Path(__file__).resolve().parents[3] / 'research/pde_ledger_v3/_measurements/mouth_coupled_conversion.json'
    reference = json.loads(record.read_text())['runs']['curved_pressure']
    if not (reference['limits_verified'] and reference['guard_exit'] == reference['worker_exit'] == 0):
        raise ValueError('Unverified zero-holder branch')
    raw = reference['literal_stdout'].encode()
    sha=hashlib.sha256(raw).hexdigest()
    if sha!='83eb66f14f91cf670cb3159cee91e15e335cb44c0376bb7d408887dbb8d8db40':
        raise ValueError('Changed zero-holder reference')
    rows=[json.loads(line) for line in raw.splitlines()]
    for action, magnitudes in ((.001,(1e-6, .5e-6)), (0.,(1e-6,))):
        base,=[row for row in rows if row.get('test')=='free_spherical_front'
               and row['brane_pressure_excess']==.0045 and row['trial_action_per_inner_area']==action]
        center=base['interface_radius']
        emit(test='reference', action=action, source_sha256=sha, zero_holder_radius=center,
             brane_pressure_excess=.0045, holder_half_width=5.,
             initial_velocity='one-sided t=0+ at chi=0.5 after removing holder')
        for strength in magnitudes:
            pair=[]
            for sign in (-1.,1.):
                holder=LocalHolder(center,sign*strength)
                pair.append(solve(200.,50.,.0045,action,center-200.,.18,8e-10,holder=holder))
            (dm,vm),(dp,vp)=pair
            emit(test='signed_pair',action=action, strength=strength, negative_offset=dm,
                 positive_offset=dp, negative_release_speed=vm, positive_release_speed=vp,
                 odd_offset=(dp-dm)/2, odd_speed=(vp-vm)/2,
                 speed_over_offset=(vp-vm)/(dp-dm))


if __name__=='__main__':
    main()
