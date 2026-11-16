"""Read momentum storage and outer flux per unit slow mouth velocity in a fixed cylinder."""

import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from scipy.fft import dct, dst, idct, idst

from lateral_mouth_radius import SpeedTranslation, Target
from lateral_radius_momentum import link_budget
from lateral_reservoir_field import LateralSlab


class FaceReconstruction:
    """Solve below a given top trace, with reflecting bottom and fixed lateral bath."""

    def __init__(self, slab):
        self.slab = slab
        levels = len(slab.sqrt_z)-1
        kz = np.pi*(np.arange(levels)+0.5)/levels
        self.eigenvalues = 2*np.sin(kz[:, None, None]/2)**2/slab.dz**2+slab.eigenvalues[0]
        self.root = slab.sqrt_z[:-1, None, None]

    def elliptic(self, rhs, trace, shift):
        s = self.slab
        rhs[-1] += trace/(2*s.dz**2)
        coeff = dst(dct(rhs*self.root, type=3, axis=0, norm="ortho"),
                    type=2, axis=1, norm="ortho")
        coeff = ((coeff*s.sqrt_r) @ s.modes)/(self.eigenvalues+shift)
        coeff = (coeff @ s.modes.T)/s.sqrt_r
        return idct(idst(coeff, type=2, axis=1, norm="ortho"),
                    type=3, axis=0, norm="ortho")/self.root

    def field(self, trace, bath, source):
        field = np.full(self.slab.eigenvalues.shape, bath, dtype=complex)
        field[-1] = trace
        for step in range(1, 501):
            b = source(field)[:-1]
            re = self.elliptic(2*(field[:-1].real-bath)+b.real, trace.real-bath, 2)
            im = self.elliptic(b.imag.copy(), trace.imag, 0)
            update = bath+re+1j*im
            change = np.max(abs(update-field[:-1]))
            field[:-1] = update
            if not np.isfinite(change):
                raise ArithmeticError("Nonfinite interior reconstruction")
            if change < 4e-14:
                return field, step
        raise ArithmeticError(f"Face reconstruction did not converge: {change}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", type=int)
    parser.add_argument("normal_cells", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--field-manifest", type=Path,
                        help="Read the saved mouth preparation from this measurement manifest")
    args = parser.parse_args()
    folder = Path("docs/_measurements")
    manifests = [args.field_manifest] if args.field_manifest is not None else [
        folder/"mouth_drift_refinement.json", folder/"mouth_drift_crossed_grid.json"]
    records = []
    for path in manifests:
        records.extend(json.loads(path.read_text())["runs"])
    record, = [row for row in records if row["cells"] == args.cells
               and row.get("normal_cells", 4*row["cells"]) == args.normal_cells]
    filename, expected = next(iter(record["field_sha256"].items()))
    with open(filename, "rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != expected:
            raise ValueError("Changed input packet")
    s = LateralSlab(24, args.cells, args.normal_cells)
    with np.load(filename, allow_pickle=False) as packet:
        if (not np.array_equal(packet["x"], s.x) or not np.array_equal(packet["r"], s.r)
                or float(packet["radius"]) != 1 or float(packet["dz"]) != s.dz):
            raise ValueError("Grid or mouth differs from archived preparation")
        archived_f, archived_u, normal = packet["field"], packet["u"], packet["normal"]

    reconstruction = FaceReconstruction(s)
    f, nf = reconstruction.field(archived_f[-1], 1., lambda v: -(abs(v)**8-1)*v/4)
    n = abs(f[-1])**2
    current = (normal/f[-1]).imag*n
    s.check_field(f, normal, [current])
    mean_speed = np.sum(s.dv*current/n)/(4*np.pi/3)
    motion = SpeedTranslation(Target(s, 1., n, mean_speed), f, current)
    p, rp = motion.static_tangent()
    u, nu = reconstruction.field(archived_u[-1], 0., lambda v: 1j*p-motion.potential_change(v))
    ru = motion.residual(u, 1j*p)
    f_error = np.max(abs(f[-4:-1]-archived_f[:-1]))
    u_error = np.max(abs(u[-4:-1]-archived_u[:-1]))
    print(f"cells={args.cells} normal_cells={args.normal_cells} f_steps={nf} u_steps={nu} "
          f"withheld_f_error={f_error:.12g} withheld_u_error={u_error:.12g} "
          f"static_residual={rp:.12g} velocity_residual={ru:.12g}", flush=True)
    if (not np.isfinite([rp, ru, f_error, u_error]).all() or max(f_error,u_error) > 3e-11
            or max(rp,ru) > 3e-9):
        raise ArithmeticError("Reconstruction disagrees with held-out layers or mouth equation")

    # At X=0 and constant imposed V, d(psi)/dt=p*V, while psi-f=u*V.
    # Keep enough full-depth data to inspect the fixed cylinders after this run.
    axial = np.flatnonzero(abs(s.x)<9)
    cut = slice(axial[0], axial[-1]+1)
    radial = slice(0, np.count_nonzero(s.r<9))
    dnormal = motion.boundary(u)[cut, radial]
    normal = normal[cut, radial]
    x, r = s.x[cut], s.r[radial]
    f, p, u = (a[:,cut,radial] for a in (f,p,u))
    with args.output.open("xb") as stream:
        np.savez_compressed(stream, x=x, r=r, dz=s.dz, field=f, p=p, u=u,
                            normal=normal, dnormal=dnormal)
    prior = dict(w.split("=") for line in record["stdout"] for w in line.split())
    print(f"archived_boundary={prior['lab']} archived_whole={prior['whole_layer']} "
          f"archived_link4={prior['law_link']}", flush=True)
    for reach in (2,4,8):
        b1, b2 = [link_budget(f,p,u,normal,dnormal,x,r,s.dx,s.dz,reach,m) for m in (1,2)]
        for order, budget in ((2,b1), (4,(4*b1-b2)/3)):
            storage, pressure, gradient, side, mouth, defect = budget
            outer = pressure+gradient+side
            inferred = -storage-outer
            balanced = inferred+defect
            error = balanced-mouth
            print(f"reach={reach} order={order} storage={storage:.12g} "
                  f"cap_pressure={pressure:.12g} cap_gradient={gradient:.12g} "
                  f"side={side:.12g} outer={outer:.12g} exterior={inferred:.12g} "
                  f"lattice_body={defect:.12g} exterior_corrected={balanced:.12g} "
                  f"direct_link={mouth:.12g} closure={error:.12g}", flush=True)
            if (not np.isfinite(budget).all() or abs(error) > 2e-11
                    or (order == 4 and abs(mouth-float(prior["law_link"])) > 2e-12)):
                raise ArithmeticError("Discrete balance or archived link flux check failed")


if __name__ == "__main__":
    main()
