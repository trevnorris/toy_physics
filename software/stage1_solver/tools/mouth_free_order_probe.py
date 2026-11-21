"""Read a local shear-stiffness discriminator and draw the measured half-order fronts."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from mouth_flared_receiver import FlaredDomain, ROUNDING
from mouth_free_order_diffuse import read


def read_fields(folder):
    outcome = json.loads((folder/"outcome.json").read_text())
    assert outcome["exitCode"] == outcome["childOutcome"]["exitCode"] == 0 and outcome["limitsVerified"]
    setup, *rows = map(json.loads,(folder/"stdout").read_text().splitlines())
    source = Path(setup["source"])
    carrier = json.loads((source.parent/"stdout").read_text().splitlines()[0])
    dom = FlaredDomain(carrier["n"],carrier["profile"],carrier["roof"],carrier["feed_far"])
    rows.sort(key=lambda row: row["midpoint"][-1]["w"])
    fronts = []
    for row in rows:
        with np.load(folder/f"order_seed{row['seed']}.npz") as data:
            chi = data["chi"]
        samples = row["midpoint"]
        bounds = max(ROUNDING+1e-6, min(s["w"] for s in samples)-.07), max(s["w"] for s in samples)+.07
        assert bounds[0] > ROUNDING and bounds[1] < 3.05
        radial_fraction=np.linspace(0,1,45)
        def order_at(z, eta):
            radius,_ = dom.radius(z)
            return float(read(dom,chi,eta*radius,z)[0])
        z = np.array([brentq(lambda z: order_at(z,t)-.5,*bounds) for t in radial_fraction])
        r = radial_fraction*dom.radius(z)[0]
        # Check the connected-level interpretation of the loop on many adjacent
        # rays. Each pair of level roots brackets the interval used for work.
        slopes=[]
        for t in radial_fraction:
            end=[brentq(lambda z: order_at(z,t)-level,*bounds) for level in (.55,.45)]
            assert end[1]>end[0]
            for w in np.linspace(*end,9):
                a,ap=dom.radius(w)
                _,cr,cw=read(dom,chi,t*a,w)
                slopes.append(cw+t*ap*cr)
        assert max(slopes)<0
        probe=order_at(1.,0.)
        if not fronts or abs(z[-1]-fronts[-1][1][-1])>1e-7:
            fronts.append((r,z))
        print(json.dumps(dict(folder=str(folder),seed=row["seed"],contact=z[-1],center=z[0],
                              probe_r=0,probe_w=1,chi_at_probe=probe,
                              rays_checked=len(radial_fraction), max_axial_order_slope_in_band=float(max(slopes)),
                              material_fractions_scope="local shear stiffness coefficient divided by ordered value")),flush=True)
    return dom, fronts


def figure(results, dest):
    fig, axes=plt.subplots(1,2,figsize=(8.3,4.45),sharex=True,sharey=True)
    near,deep="#00756e","#b24a25"
    for ax,(dom,fronts) in zip(axes,results,strict=True):
        assert len(fronts)==2
        z=np.linspace(dom.lip,2,501)
        w=np.maximum(z,ROUNDING)
        a,_=dom.radius(w)
        lower=z<ROUNDING
        a[lower]=dom.rc-np.sqrt(np.maximum(0,ROUNDING**2-(z[lower]-ROUNDING)**2))
        ax.fill_betweenx(z,-a,a,color="#f1f6fa")
        for sign in (-1,1):
            ax.plot(sign*a,z,color="#3f4e56",lw=1.8)
        for (r,w),color,label in zip(fronts,(near,deep),("Near equilibrium","Deep equilibrium"),strict=True):
            ax.plot(np.r_[-r[:0:-1],r],np.r_[w[:0:-1],w],lw=2.9,color=color,label=label)
            ax.scatter([-r[-1],r[-1]],[w[-1]]*2,s=19,color=color,zorder=4)
        ax.scatter(0,1.,color="#222222",marker="x",s=30,label="Shear probe")
        ax.text(0,.18,"feeding side",ha="center",fontsize=9,color="#52626b")
        ax.set_title("Earlier flare" if dom.profile=="cubic" else "Delayed flare",fontsize=13)
        ax.set(xlim=(-1.45,1.45),ylim=(0,1.9),aspect="equal",xlabel="Across the throat ($R$)")
        ax.spines[["top","right"]].set_visible(False)
    axes[0].set_ylabel("Depth into the bulk, $w/R$")
    handles,labels=axes[0].get_legend_handles_labels()
    fig.subplots_adjust(top=.81,bottom=.21,wspace=.18)
    fig.legend(handles,labels,loc="lower center",ncol=3,frameon=False)
    fig.suptitle("Where shear order falls to one half",fontsize=15)
    fig.savefig(dest,dpi=170)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("folders",nargs="+",type=Path)
    p.add_argument("--figure",type=Path)
    args=p.parse_args()
    results=[read_fields(d) for d in args.folders]
    if args.figure:
        figure(results,args.figure)
