"""Lower-dimensional control: let the measured material current carry transverse momentum.

Model the round cross-section by its volume, allow a fixed finite transverse
inertia per 4-volume, and use the section-averaged saved shear gate. The
carrier has unit-density plug flow Q/A. The trial has an anechoic feed and
an advecting downstream section with zero shear stiffness.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from mouth_free_order_wave_inertia import ARCHIVE, FEED, packet, sections, verified
from mouth_flared_receiver import EXIT, FlaredDomain


def run(binary, folder, counts):
    rows=list(map(json.loads,(folder/"stdout").read_text().splitlines()))
    setup=rows[0]
    state=json.loads((Path(setup["source"]).parent/"stdout").read_text().splitlines()[0])
    dom=FlaredDomain(state["n"],state["profile"],state["roof"])
    hashes=json.loads(ARCHIVE.read_text())["diffuse_orders_sha256"]
    w=np.linspace(FEED,2.,2801)
    q=.25
    for row in rows[1:]:
        if row["seed"] not in (.44,1.55):
            continue
        order=folder/f"order_seed{row['seed']}.npz"
        verified(order,hashes)
        with np.load(order) as fields:
            _,chi=sections(dom,fields["chi"],w)
        for n in counts:
            # Uniform anechoic entrance lead, downstream end is the actual
            # axial exit. Beyond w=2 the trial carries transverse momentum
            # with the material only. Taper the already weak gate over .04R.
            lo=-.20
            dx=(EXIT-lo)/n
            z=lo+dx*(np.arange(n)+.5)
            zs=np.maximum(z,FEED)
            radius=dom.radius(zs)[0]
            low=zs<dom.rounding
            radius[low]=dom.rc-np.sqrt(dom.rounding**2-(zs[low]-dom.rounding)**2)
            area=4*np.pi*radius**3/3
            tail=np.clip((2.04-z)/.04,0,1)
            gate=np.interp(zs,w,chi)*tail**2*(3-2*tail)
            v=packet(z)*1e-4
            inp=[f"{n} {dx:.17g} {q}"]
            inp.extend(f"{m:.17g} {k:.17g} {speed:.17g}" for m,k,speed in zip(area,area*gate,v))
            out=subprocess.run([str(binary)],input="\n".join(inp),text=True,capture_output=True,check=True)
            values=np.loadtxt(out.stdout.splitlines())
            print(json.dumps(dict(folder=str(folder),seed=row["seed"],n=n,mean_gate_at_taper=chi[-1],
                columns=["t","remaining","feed","advected_across_exit","work_from_background","energy_defect"],
                samples=values.tolist())),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    parser.add_argument("--cells",nargs="+",type=int,default=[400,800])
    args=parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        exe=Path(tmp)/"material_wave"
        source=Path(__file__).with_suffix(".cpp")
        subprocess.run(["g++","-O3","-std=c++17",str(source),"-o",str(exe)],check=True)
        run(exe,args.folder,args.cells)
