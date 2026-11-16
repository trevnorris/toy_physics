"""Fixed companion, variable target. Measure paired-minus-alone opening loads.

Use the same normal-speed calibration as the source-rate experiment. A
separate source-only solution measures incident velocity at the absent target.
The target stays centered in the same finite reservoir for every subtraction.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from lateral_reservoir_field import LateralSlab
from lateral_reservoir_readout import mouth_load
from lateral_source_rate import check, solve_pair, verified_packet
from mouth_feed_weak_check import dx4


def profile(slab, center, density):
    weight = np.maximum(1-(slab.x[:, None]-center)**2-slab.r**2, 0)**3
    offset = round(center/slab.dx)
    if offset*slab.dx != center:
        raise ValueError("Mouth translation must be grid-aligned")
    return weight, weight/np.roll(density, offset, axis=0)


def single(slab, center, rate, density):
    weight, shape = profile(slab, center, density)

    def boundary(face):
        n = abs(face)**2
        raw = shape*n
        current = rate*raw/np.sum(raw*slab.dv)
        normal = (-weight*(1-1/abs(face))+1j*current/n)*face
        return current[None], normal

    f = np.ones(slab.eigenvalues.shape, complex)
    for step in range(1,501):
        _, normal = boundary(f[-1])
        potential = (abs(f)**8-1)/4
        re = slab.elliptic(2*(f.real-1)-potential*f.real,normal.real,2)
        im = slab.elliptic(-potential*f.imag,normal.imag,0)
        updated = 1+re+1j*im
        error = np.max(abs(updated-f))
        if not np.isfinite(error):
            raise ArithmeticError("Nonfinite isolated field")
        f = updated
        if error < 5e-13:
            currents, normal = boundary(f[-1])
            return f, currents, normal, step
    raise ArithmeticError(f"Isolated field failed to converge: {error}")


def checked_solve(slab, reference, centers, rates, path):
    rates = np.array(rates)
    if len(centers)==1:
        f, currents, normal, steps = single(slab, centers[0], rates[0], reference)
    else:
        f, currents, normal, steps = solve_pair(slab, -centers[0], rates, reference)
    residual, bath, measured = check(slab, f, normal, currents, rates)
    for center, current in zip(centers, currents, strict=True):
        _, shape = profile(slab, center, reference)
        cut = shape>0
        ratio = current[cut]/(shape[cut]*abs(f[-1][cut])**2)
        if (np.any(current[~cut]!=0) or np.ptp(ratio)>2e-12 or not np.isfinite(ratio).all()):
            raise ArithmeticError("Changed intake profile")
    layers = f[-4:].copy()
    np.savez_compressed(path,layers=layers,currents=currents,x=slab.x,r=slab.r,dz=slab.dz,
                        centers=centers,rates=rates)
    print(f"  field={path.stem} steps={steps} residual={residual:.4g} bath={bath:.12g} "
          f"rates={','.join(format(j,'.12g') for j in measured)}",flush=True)
    return layers, currents


def force(slab, layers, current, target_rate):
    face = layers[-1]
    n = abs(face)**2
    weight = np.maximum(1-slab.x[:, None]**2-slab.r**2, 0)**3
    cut = weight>0
    if abs(np.sum(current*slab.dv)-target_rate)>2e-13 or np.any(current[~cut]!=0):
        raise ArithmeticError("Mismatched target for force reading")
    dxf = dx4(face,slab.dx)
    vx = (face.conj()*dxf).imag/n
    slope = -2*weight*(1-1/np.sqrt(n))
    carried = -np.sum(slab.dv*current*vx)
    stress = -np.sum(slab.dv*slope*dx4(n,slab.dx))/4
    # Direct off-diagonal momentum current from the complex field, including
    # the density term. The source lies along the negative x direction.
    def normal(a):
        return (11*a[-1]-18*a[-2]+9*a[-3]-2*a[-4])/(6*slab.dz)
    flux = (dxf.conj()*normal(layers)).real - dx4(normal(abs(layers)**2),slab.dx)/4
    direct = -np.sum((slab.dv*flux)[cut])
    return np.array([carried,stress,direct])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(k):v for run in archive["runs"] for k,v in run.get("fields",{}).items()}
    with verified_packet(args.reference/"d0.npz",hashes) as p:
        x, r, dz = p["x"], p["r"], float(p["dz"])
        slab = LateralSlab(round(r[-1]+r[0]),round(1/(2*r[0])),round(1/dz))
        if int(p["d"])!=0 or not np.array_equal(x,slab.x) or not np.array_equal(r,slab.r) or slab.dz!=dz:
            raise ValueError("Archived isolated grid mismatch")
        field = p["field"]
        reference = abs(field[-1])**2
        baseline = {0.25: (field[-4:].copy(),p["currents"])}
        del field
    args.output.mkdir(exist_ok=False,parents=True)
    rates = (0.125,0.25,0.5)
    print(f"extent={r[-1]+r[0]:g} dx={slab.dx:.12g} dz={dz:.12g}",flush=True)
    alone = {}
    for target in rates:
        if target not in baseline:
            baseline[target] = checked_solve(slab,reference,[0],[target],args.output/f"alone_J{target:g}.npz")
        layers, currents = baseline[target]
        alone[target] = force(slab,layers,currents[0],target)
        if max(abs(alone[target]))>1e-11:
            raise ArithmeticError("Isolated target has a nonzero lateral baseline")
        print(f"  alone_Jt={target:g} carried={alone[target][0]:.12g} stress={alone[target][1]:.12g} "
              f"direct={alone[target][2]:.12g}",flush=True)
    totals = {}
    for d in (4,8):
        incident, _ = checked_solve(slab,reference,[-d],[0.25],args.output/f"source_d{d}.npz")
        face = incident[-1]
        incident_speed = -(face.conj()*dx4(face,slab.dx)).imag/abs(face)**2
        for target in rates:
            if target==0.25:
                with verified_packet(args.reference/f"d{d}.npz",hashes) as p:
                    if int(p["d"])!=d or not np.array_equal(p["x"],x) or not np.array_equal(p["r"],r) or float(p["dz"])!=dz:
                        raise ValueError("Archived pair grid mismatch")
                    layers, currents = p["layers"],p["currents"]
            else:
                layers,currents = checked_solve(slab,reference,[-d,0],[0.25,target],
                                                args.output/f"d{d}_J{target:g}.npz")
            measured = np.sum(currents*slab.dv,axis=(1,2))
            if len(currents)!=2 or max(abs(measured-[0.25,target]))>2e-13:
                raise ArithmeticError("Pair rates changed before readout")
            paired = force(slab,layers,currents[-1],target)
            if target==0.25:
                old = mouth_load(layers,currents,x,r,dz,2)
                if max(abs(paired[:2]-old[:2]))>1e-13:
                    raise ArithmeticError("Boundary-law force disagrees with archived method")
            carried,stress,direct = paired-alone[target]
            total = carried+stress
            # Fixed isolated-target weighting on a companion-only field.
            available = np.sum(slab.dv*baseline[target][1][0]*incident_speed)
            if not np.isfinite([total,direct,available]).all() or available==0:
                raise ArithmeticError("Unusable force or incident reference")
            totals[d,target] = total
            print(f"  d={d} J_t={target:g} carried={carried:.12g} stress={stress:.12g} "
                  f"total={total:.12g} direct={direct:.12g} incident_reference={available:.12g} "
                  f"load_over_incident={total/available:.12g} scaled={d*d*total/(0.25*target):.12g}",flush=True)
        print(f"  d={d} doubling_low={totals[d,0.25]/totals[d,0.125]:.12g} "
              f"doubling_high={totals[d,0.5]/totals[d,0.25]:.12g}",flush=True)
    for target in rates:
        print(f"  J_t={target:g} F8_over_F4={totals[8,target]/totals[4,target]:.12g}")


if __name__=="__main__":
    main()
