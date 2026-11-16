"""Change the companion's rate at fixed mouth geometry and target preparation."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from lateral_reservoir_field import LateralSlab
from lateral_reservoir_readout import mouth_load


def verified_packet(path, hashes):
    if sha256(path.read_bytes()).hexdigest() != hashes[path]:
        raise ValueError(f"Input field checksum mismatch: {path}")
    return np.load(path, allow_pickle=False)


def solve_pair(slab, d, rates, base_density):
    # Use the isolated rate-0.25 calibration for the shape of the normal
    # speed at both mouths. Only its scalar rate normalization is varied.
    shapes = []
    restoring = np.zeros_like(base_density)
    for center in (-d, 0):
        weight = np.maximum(1-(slab.x[:, None]-center)**2-slab.r**2, 0)**3
        restoring += weight
        shift = round(center/slab.dx)
        if shift*slab.dx != center:
            raise ValueError("Shift must coincide with a grid point")
        shapes.append(weight/np.roll(base_density, shift, axis=0))
    shapes = np.array(shapes)

    def boundary(face):
        n = abs(face)**2
        raw = shapes*n
        totals = np.sum(slab.dv*raw, axis=(1, 2))
        currents = raw*(rates/totals)[:, None, None]
        normal = (-restoring*(1-1/abs(face)) + 1j*currents.sum(axis=0)/n)*face
        return currents, normal

    f = np.ones(slab.eigenvalues.shape, complex)
    for step in range(1, 501):
        _, normal = boundary(f[-1])
        potential = (abs(f)**8-1)/4
        re = slab.elliptic(2*(f.real-1)-potential*f.real, normal.real, 2)
        im = slab.elliptic(-potential*f.imag, normal.imag, 0)
        updated = 1+re+1j*im
        error = np.max(abs(updated-f))
        if not np.isfinite(error):
            raise ArithmeticError("Nonfinite update")
        f = updated
        if error < 5e-13:
            currents, normal = boundary(f[-1])
            return f, currents, normal, step
    raise ArithmeticError(f"Steady pair did not converge: {error}")


def check(s, f, normal, currents, rates):
    # Discrete real-space residual including axial, radial and normal faces.
    dx, dz = s.dx, s.dz
    radial_faces = np.arange(len(s.r)+1)*dx
    end = np.concatenate((2-f[:, :1], f, 2-f[:, -1:]), axis=1)
    lap = np.diff(end, n=2, axis=1)/dx**2
    side = np.concatenate((f[:, :, :1], f, 2-f[:, :, -1:]), axis=2)
    lap += np.diff(radial_faces*np.diff(side, axis=2), axis=2)/(s.r*dx**2)
    face = np.concatenate((f[1:2], f, (f[-2]+2*dz*normal)[None]), axis=0)
    lap += np.diff(face, n=2, axis=0)/dz**2
    residual = np.max(abs(-lap/2+(abs(f)**8-1)*f/4))
    ends = np.sum((f[:, 0].imag+f[:, -1].imag)*2*np.pi*s.r*dx, axis=1)
    wall = 2*np.pi*radial_faces[-1]*dx*f[:, :, -1].imag.sum(axis=1)
    inflow = dz*np.sum(s.sqrt_z**2*(ends+wall))*2/dx
    measured = np.sum(currents*s.dv, axis=(1, 2))
    if (not np.isfinite([residual, inflow, *measured]).all() or residual>1e-8
            or max(abs(measured-rates))>3e-14 or abs(inflow-rates.sum())>2e-8):
        raise ArithmeticError(f"Unverified field: {residual=} {inflow=} {measured=}")
    return residual, inflow, measured


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Directory with the matched archived fields")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(k):v for run in archive["runs"] for k,v in run.get("fields",{}).items()}
    with verified_packet(args.reference/"d0.npz", hashes) as p:
        x, r, dz = p["x"], p["r"], float(p["dz"])
        slab = LateralSlab(round(r[-1]+r[0]),round(1/(2*r[0])),round(1/dz))
        if int(p["d"])!=0 or not np.array_equal(slab.x,x) or not np.array_equal(slab.r,r) or slab.dz!=dz:
            raise ValueError("Invalid isolated grid")
        field = p["field"]
        density = abs(field[-1])**2
        alone = mouth_load(field[-4:],p["currents"],x,r,dz,1)
        del field
    args.output.mkdir(exist_ok=False, parents=True)
    print(f"extent={r[-1]+r[0]:g} dx={slab.dx:.12g} dz={dz:.12g} "
          f"alone_force={alone[:2].sum():.12g}", flush=True)
    force = {}
    for d in (4, 8):
        for rate in (0.125, 0.25, 0.5):
            if rate == 0.25:
                with verified_packet(args.reference/f"d{d}.npz", hashes) as p:
                    if (int(p["d"])!=d or not np.array_equal(x,p["x"]) or
                            not np.array_equal(r,p["r"]) or dz!=float(p["dz"])):
                        raise ValueError("Archived pair differs from reference geometry")
                    layers, currents = p["layers"],p["currents"]
                origin = "archive"
            else:
                rates = np.array([rate, 0.25])
                f, currents, normal, step = solve_pair(slab,d,rates,density)
                residual, bath, measured = check(slab,f,normal,currents,rates)
                layers = f[-4:].copy()
                del f
                np.savez_compressed(args.output/f"d{d}_J{rate:g}.npz",layers=layers,currents=currents,
                                    x=x,r=r,dz=dz,d=d,source_rate=rate,target_rate=0.25)
                print(f"  d={d} steps={step} residual={residual:.4g} bath={bath:.12g} "
                      f"source={measured[0]:.12g} target={measured[1]:.12g}", flush=True)
                origin = "new"
            if max(abs(np.sum(currents*slab.dv,axis=(1,2))-[rate,0.25]))>3e-14:
                raise ArithmeticError("Intake rates changed before force readout")
            parts = mouth_load(layers,currents,x,r,dz,2)-alone
            if not np.isfinite(parts).all():
                raise ArithmeticError("Nonfinite mouth load")
            total, direct = parts[:2].sum(), parts[2:].sum()
            force[d,rate] = total
            print(f"  d={d} J_s={rate:g} origin={origin} carried={parts[0]:.12g} stress={parts[1]:.12g} "
                  f"total={total:.12g} direct={direct:.12g} F_over_Js={total/rate:.12g}",flush=True)
        if force[d,0.125]!=0 and force[d,0.25]!=0:
            print(f"  d={d} doubling_low={force[d,0.25]/force[d,0.125]:.12g} "
                  f"doubling_high={force[d,0.5]/force[d,0.25]:.12g}",flush=True)
    for rate in (0.125,0.25,0.5):
        if force[4,rate]!=0:
            print(f"  J_s={rate:g} F8_over_F4={force[8,rate]/force[4,rate]:.12g}")


if __name__ == "__main__":
    main()
