"""Check identification of radius-four bank tractions from saved wall work and momentum.

The manufactured flux is a null perturbation of cell balances, without a
constitutive stress law. It is not an alternative physical solution.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def same_trace(rows):
    output = []
    keys = sorted({(row["strip"], row["bank"], row["quadrature"]) for row in rows})
    for strip, bank, order in keys:
        four, eight = sorted(
            (row for row in rows if (row["strip"], row["bank"], row["quadrature"]) ==
             (strip, bank, order)), key=lambda row: row["radial_cells"]
        )
        assert [four["radial_cells"], eight["radial_cells"]] == [4, 8]
        assert four["depth_knots"] == eight["depth_knots"]
        assert abs(four["outward_wall_axial"] - eight["outward_wall_axial"]) < 1e-12
        difference = four["total_native"] - eight["total_native"]
        output.append(dict(strip=strip, bank=bank, quadrature=order, four_minus_eight=difference,
                           face_only_minimax_bound=abs(difference)/2,
                           signed_step_spread=max(abs(r["direct_reaction"][0]-r["direct_reaction"][1])
                                                  for r in (four,eight))))
    return output


def check_null(axes, z_reference, measurements, rows):
    x,z,a = (axes[key] for key in ("inside_x","duct_z","angles"))
    lo, hi = sorted((x[0],x[-1]))
    assert x[0] == lo and hi > lo and np.all(np.diff(z)>0)
    # The two terminal doubled-grid cells contribute -1/2, +3/2. Their
    # mean weights in the wall tent are 3/8, 1/8. Restore the far flux
    # through a compensating bank load well upstream of all measured tents.
    packets = [(z_reference[1],z_reference[2],-1.),
               (z_reference[-3],z_reference[-2],-.5),
               (z_reference[-2],z_reference[-1],1.5)]
    primitive = sum(total*np.clip((z-first)/(last-first),0,1)
                    for first,last,total in packets)
    for first,last,_ in packets:
        assert min(abs(z-first))<1e-12 and min(abs(z-last))<1e-12
    b = np.clip((x-0.9*hi)/(0.1*hi),0,1)
    psi = np.outer(b,primitive)
    radial = np.diff(psi,axis=1)
    axial = -np.diff(psi,axis=0)
    balance = np.diff(radial,axis=0)+np.diff(axial,axis=1)
    # The integrated correction may be distributed on angular cells with
    # positive normalized sin(theta) weights and zero angular-face flux.
    angle_shares = abs(np.diff(np.cos(a)))
    angle_shares /= angle_shares.sum()
    max_balance = float(abs(balance).max()*angle_shares.max())
    end_max = float(max(abs(axial[:,0]).max(),abs(axial[:,-1]).max(),abs(radial[0]).max()))
    assert max_balance<1e-12 and end_max<1e-12
    tents={}
    for row in rows:
        first,middle,last=row["depth_knots"]
        tent = np.maximum(0,np.minimum((z-first)/(middle-first),(last-z)/(last-middle)))
        tents[row["strip"]]=float(radial[-1] @ ((tent[1:]+tent[:-1])/2))
    surfaces=[]
    for measure in measurements:
        low, high = sorted(measure["lower_half_thickness"])
        A=np.clip((z-(z[-1]-high))/(high-low),0,1)
        bank=float(radial[-1] @ ((A[1:]+A[:-1])/2))
        exterior=float(-np.sum(np.diff(A)*(axial[:,1:]+axial[:,:-1])/2))
        surfaces.append(dict(axial_band=measure["axial_band"],bank_duct=bank,bank_receiving=0,
                             enclosing_flux=exterior,net=bank+exterior))
    assert max(abs(x) for x in tents.values())<1e-12
    assert all(abs(r["bank_duct"]-1)<1e-12 and abs(r["net"])<1e-12 for r in surfaces)
    return dict(cell_balance_max=max_balance,end_and_axis_flux_max=end_max,
                depth_packets=packets,wall_tents=tents,closed_surfaces=surfaces)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("fields",type=Path)
    p.add_argument("wall_work",type=Path)
    args=p.parse_args()
    work=json.loads(args.wall_work.read_text())["data"]
    with np.load(args.fields/"double"/"fields.npz") as axis:
        ref=axis["duct_z"]
    result=dict(wall_work_sha256=digest(args.wall_work),levels={})
    for level in ("double","local"):
        directory=args.fields/level
        archived=work[level]
        outcome=json.loads((directory/"outcome.json").read_text())
        field=directory/"fields.npz"
        sha=digest(field)
        assert outcome["exitCode"]==0 and outcome["limitsVerified"] and outcome["childOutcome"]["exitCode"]==0
        assert archived["guard_exit"]==0 and sha==archived["source"]["field"]
        records=[json.loads(line) for line in (directory/"stdout").read_text().splitlines()]
        surfaces=[r for r in records if r["kind"]=="closed_work" and r["lateral_band"]==4]
        assert sorted(r["axial_band"] for r in surfaces)==[2,3]
        with np.load(field) as axes:
            null=check_null(axes,ref,surfaces,archived["results"])
        result["levels"][level]=dict(field_sha256=sha,same_trace=same_trace(archived["results"]),
                                     conservation_null_control=null)
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
