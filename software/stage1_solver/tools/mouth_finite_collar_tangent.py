"""Check rank-one stability of the registered collar and its overlap with the deep cavity wave."""
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss, legval

from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import ROOT, PORT, occupied, read_case


def directional_margin(data):
    # Legendre-Hadamard rank-one test in all material directions, on each saved quadrature point.
    rs, rz, ws, wz, hoop = data["gradient"].T
    jac, chi = data["J"], data["chi"]
    A = np.stack((rs, rz, ws, wz), axis=-1).reshape(-1, 2, 2)
    values = np.linalg.svd(A, compute_uv=False)
    greatest = np.maximum(values[:, 0], hoop)
    least = np.minimum(values[:, -1], hoop)
    trace = np.sum(A*A, axis=(1, 2)) + 2*hoop**2
    k = 3*trace/(8*np.sqrt(jac)) + 1 - np.log(jac)
    margin = k - (greatest-least)**2/(4*np.sqrt(jac))
    index = np.argmin(margin)
    return dict(minimum=float(margin[index]), min_pure_longitudinal=float(k.min()),
                negative_points=int(sum(margin<0)), points=len(jac),
                max_order_on_negative_points=float(max(chi[margin<0], default=0.)),
                order_at_witness=float(chi[index]), J_at_witness=float(jac[index]),
                coordinates_at_witness=[float(data[key][index]) for key in ("a", "t")],
                reference_order_weighted_negative_fraction=float(
                    np.average(margin<0, weights=data["weights_reference_with_order"])))


def cavity_overlap(radius, shape, coefficients, prep, chi, seed, ref, order_n, degree):
    mesh = prep["field_mesh"]
    prior = complex(ref["frequency"], -ref["amplitude_loss"])
    omega, vector, _ = occupied(shape, order_n, mesh, chi, prior, seed)
    if abs(omega.real-prior.real)>1e-8:
        raise ArithmeticError("Archived symmetric pole changed")
    wave = OccupiedWave(shape, order_n, 2*order_n, mesh, PORT)
    weights = np.asarray(wave.project @ chi.T.ravel()).reshape(mesh, -1)
    nz = weights.shape[1]
    values = np.r_[vector,0j][wave.ids].reshape(mesh, nz, 4)
    target = np.vdot(vector, wave.matrix(weights.ravel()[:, None, None]*wave.mc) @ vector).real
    x0 = np.arange(mesh)[:, None]/mesh
    q, qw = leggauss(degree)
    full = collar = 0.
    for zb, wb in zip((q+1)/2, qw/2):
        t = (np.arange(nz)[None, :]+zb)/nz
        radius_at_t = shape.profile(t)[0]
        cut = np.exp(-.5*np.exp(legval(2*t-1, coefficients)))
        for region in ("all", "ordered_collar"):
            left = x0 if region=="all" else np.maximum(x0,cut)
            width = np.maximum(x0+1/mesh-left, 0.)
            for xa, wa in zip((q+1)/2, qw/2):
                x = left + xa*width
                a = (x-x0)*mesh
                v = ((1-a)*(1-zb)*values[:,:,0]+a*(1-zb)*values[:,:,1]
                     +(1-a)*zb*values[:,:,2]+a*zb*values[:,:,3])
                measure = wa*wb*shape.length/nz * radius_at_t * width * (x*radius_at_t)**4
                integral = np.sum(measure*abs(v)**2*weights)
                if region=="all":
                    full += integral
                else:
                    collar += integral
    return dict(overlap_fraction=float(collar/full), full_mass_relative_error=float(full/target-1))


def main():
    archive = json.loads((ROOT/"research/pde_ledger_v3/_measurements/mouth_finite_bore_support.json").read_text())
    for radius in (1, 2):
        prep, chi, seed, ref, _, _, _ = read_case(radius, 72)
        for depth in ("long", "unit"):
            for side in ("ordered", "receiver"):
                field = f"r{radius}_{depth}_g2_{side}.npz"
                file = ROOT/"_scratch/s11c/finite_bore_graph01/fields"/field
                if hashlib.sha256(file.read_bytes()).hexdigest() != archive["central_fields"][field]:
                    raise ValueError("Collar field fingerprint changed")
                with np.load(file) as data:
                    result = directional_margin(data)
                    shape = Meridian(radius,data["coordinates"])
                    if side=="ordered" and depth=="long":
                        result["wave_tests"] = [cavity_overlap(radius, shape, data["coeffs"],
                                                              prep, chi, seed, ref, 72, g)
                                                for g in (3, 5)]
                print(json.dumps(dict(radius=radius,depth=depth,side=side,**result),
                                 allow_nan=False),flush=True)
                tests = result.get("wave_tests", [])
                if tests and (max(abs(t["full_mass_relative_error"]) for t in tests)>2e-5
                              or abs(tests[0]["overlap_fraction"]-tests[1]["overlap_fraction"])>5e-4):
                    raise ArithmeticError("Wave-sector quadrature unresolved")


if __name__ == "__main__":
    main()
