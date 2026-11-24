"""Exercise the equilibrium admission condition with inverted shared spatial labels."""
import argparse
import ctypes as ct
import json
from pathlib import Path

import numpy as np

from mouth_cp10_spatial_crystal import Data, Pointer, integrate, load_context
from mouth_exact_checkpoint_control import digest


class OrderedInversion(ArithmeticError):
    def __init__(self, witness):
        super().__init__('Nonpositive inverse determinant on ordered quadrature support')
        self.witness = witness


def admissible_energy(kernel, body, labels, order, degree=3):
    totals, witness, _, _ = integrate(kernel, [body], labels, order, degree)
    if witness[13] != 0:
        raise OrderedInversion(witness)
    return totals[0]


def control(kernel, body, labels, order, core):
    # Pick an actually ordered cell interior; all cells touched by its nodal
    # perturbation must lie wholly on the requested side of the original weld.
    split = body.x[(len(body.x)-1)*2//3]
    radial = labels[...,0]
    cells = np.lib.stride_tricks.sliding_window_view(order,(2,2,2))
    for i,j,k in np.argwhere(cells.min(axis=(-1,-2,-3))>1.e-3):
        if not (1<i<len(radial)-3 and 0<j<order.shape[1]-3 and 1<k<order.shape[2]-3):continue
        near = radial[i-1:i+3,j-1:j+3,k-1:k+3]
        if core and near.max()>=split or not core and near.min()<=split: continue
        sl = np.s_[i:i+2,j:j+2,k:k+2]
        changed = labels.copy()
        rr = radial[sl]
        changed[sl+(0,)] = rr.max()+rr.min()-rr
        try: admissible_energy(kernel,body,changed,order)
        except OrderedInversion as e:
            w=e.witness
            if bool(w[12])!=core or w[4]<=0 or w[5]>=0:
                raise ArithmeticError('Perturbation did not isolate the selected side of the weld')
            return dict(region='core' if core else 'collar',perturbed_cell=[int(i),int(j),int(k)],
                        changed_nodal_components=int(np.count_nonzero(labels!=changed)),
                        accepted=False,reason=str(e),inverted_ordered_samples=int(w[13]),
                        witness_cell=w[9:12].astype(int).tolist(),witness_physical_chart=w[6:9].tolist(),
                        chi=w[4],Jinv=w[5])
        raise ArithmeticError('Inverted shared nodal labels were accepted')
    raise ArithmeticError('No wholly ordered patch on requested side')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','base','labels','out'): p.add_argument(key,type=Path)
    a=p.parse_args(); root=a.root.resolve();base=a.base.resolve();lab=a.labels.resolve()
    ctx=load_context(root,lab)
    k=ct.CDLL(str(base/'crystal.so')).cp10_spatial_integral
    k.argtypes=[ct.POINTER(Data),ct.c_int]+[Pointer]*5+[ct.c_int]+[Pointer]*4
    k.restype=ct.c_int
    out=dict(labels=digest(lab),sides=[control(k,ctx[2],ctx[-2],ctx[-1],side) for side in (True,False)])
    with a.out.open('x') as f:json.dump(out,f,indent=2, allow_nan=False); f.write('\n')
