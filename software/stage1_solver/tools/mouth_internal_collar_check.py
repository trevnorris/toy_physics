"""Check native work on the saved material map and its continuation into melted cells."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_internal_collar import CollarPair, Meridian


def read(path):
    with np.load(path,allow_pickle=False) as state:
        if int(state["radius"])!=4:
            raise ValueError("Expected a radius-four material state")
        shape=Meridian(4.,state["shape_params"],float(state["shape_extra"]))
        labels=state["labels"]
        chart=CollarPair(shape,state["registration"],(labels.shape[1]-1)//2,labels.shape[0]-1)
        order=state["order"]
        inner=(np.linspace(0,1,order.shape[0]),np.linspace(0,1,order.shape[1]),order)
        outer=(state["outer_t"],state["outer_x"],state["outer_order"])
    return chart,labels,inner,outer


def measure(chart,labels,inner,outer):
    chi,_=chart._measure(inner,outer)
    native=chart.evaluate(labels,inner,outer,True)
    density=chart.energy_density(labels,chart.r,chart.L*chart.z)
    # The direct spatial reference uses dV=4 pi r^2 |dr/de| L de dt.
    density_energy=np.dot(chart.qw*chart.r**2*(-chart.re)*chart.L*chi,density)
    cells=len(chart.triangles)
    vertices=labels.reshape(-1,2)[chart.triangles]
    h=(vertices[:,1]-vertices[:,0])*2*chart.na
    v=(vertices[:,2]-vertices[:,0])*chart.nt
    d=h[:,0]*v[:,1]-v[:,0]*h[:,1]
    phi=chi.reshape(cells,-1)
    density_q=(chart.qw*chart.r**2*(-chart.re)*chart.L*chi*density).reshape(cells,-1)
    report={}
    eq=chart.eta.reshape(cells,-1);tq=chart.z.reshape(cells,-1)
    for tag,mask in (("empty",phi.max(axis=1)==0),("weak",(phi.max(axis=1)>0)&(phi.max(axis=1)<1.e-4)),
                     ("present",phi.max(axis=1)>=1.e-4),("ordered",phi.min(axis=1)>=.5)):
        if np.any(mask):
            j=np.flatnonzero(mask)[np.argmin(d[mask])]
            report[tag]=dict(count=int(sum(mask)),det_min=float(d[j]),e=float(eq[j].mean()),t=float(tq[j].mean()),
                             order=[float(phi[j].min()),float(phi[j].max())],energy=float(density_q[mask].sum()))
    return dict(native={k:v for k,v in native.items() if k not in ("gradient","inner_source","outer_source")},
                spatial_energy=float(density_energy),density_minus_native=float(density_energy-native["energy"]),area=report)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("states",nargs="+",type=Path)
    p.add_argument("--gauge",nargs="+",type=float,default=[])
    args=p.parse_args()
    for state in args.states:
        o=json.loads((state.parent/'outcome.json').read_text())
        if not o["limitsVerified"] or o['exitCode'] or o['childOutcome']['exitCode']:
            raise ValueError("Unvalidated saved state")
        chart,labels,inner,outer=read(state)
        ident=dict(source_sha256=hashlib.sha256(state.read_bytes()).hexdigest(),mesh=[chart.na,chart.nt])
        print(json.dumps(dict(**ident,gauge=None,**measure(chart,labels,inner,outer)),allow_nan=False),flush=True)
        for gauge in args.gauge:
            other,solve=chart.solve(inner,outer,initial=labels,gauge=gauge)
            diag=measure(chart,other,inner,outer)
            print(json.dumps(dict(**ident,gauge=gauge,optimizer=solve["optimizer_status"],iterations=solve["iterations"],
                                  **diag),allow_nan=False),flush=True)
