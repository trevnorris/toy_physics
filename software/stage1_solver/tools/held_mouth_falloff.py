"""Check the held-pair distance law against the isolated drain and the slab phase equation."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
from scipy.fft import dctn, idctn

from mouth_feed_weak_check import dx4


def linear_face_phase(current, dx, dz, depth):
    # For a lateral mode with discrete wavenumber k, the depth recurrence
    # has a solution proportional to sinh(i*theta) at depth grid index i.
    # The centred flux condition at the mouth gives phi_N/j =
    # dz*tanh(N*theta)/sinh(theta). The zero lateral mode has phi_N/j = H.
    size = current.shape[0]
    k2 = (2*np.sin(np.pi*np.arange(size)/(2*size))/dx)**2
    k = np.sqrt(k2[:, None, None]+k2[None, :, None]+k2[None, None, :])
    theta = 2*np.arcsinh(dz*k/2)
    kernel = np.full_like(k, depth)
    nonzero = k > 0
    kernel[nonzero] = dz*np.tanh((depth/dz)*theta[nonzero])/np.sinh(theta[nonzero])
    return idctn(kernel*dctn(current, norm="ortho"), norm="ortho")


def compare(path):
    with np.load(path, allow_pickle=False) as packet:
        dx, dz, depth = (float(packet[key]) for key in ("dx", "dz", "depth"))
        if int(packet["separation"]) != 0:
            raise ValueError("Only isolated-mouth fields enter the prediction")
        current = packet["current"]
        face = packet["layers"][-1]
    def reflect_x(f):
        return np.concatenate((f[::-1], f), axis=0)
    phase = linear_face_phase(current, dx, dz, depth)
    face = reflect_x(face)
    isolated_v = (face.conj()*dx4(face, dx)).imag/abs(face)**2
    linear_v = dx4(reflect_x(phase), dx)
    weight = reflect_x(current)*4*dx**3
    if not np.isclose(weight.sum(), 0.25, rtol=0, atol=2e-13):
        raise ValueError("Target weighting has the wrong intake")
    # Translate the isolated target intake without solving any paired state.
    # Positive is toward the source at x=0 for a target centred at x=d.
    loads = {}
    for d in (3, 4, 6):
        target = np.roll(weight, round(d/dx), axis=0)
        loads[d] = [-np.sum(target*v) for v in (linear_v, isolated_v)]
        print(f"  d={d} linear_carried={loads[d][0]:.12g} isolated_carried={loads[d][1]:.12g}")
    q = np.pi/(2*depth)
    first_mode = np.exp(-3*q)*(1+6*q)/(4*(1+3*q))
    print(f"  F6_over_F3 lowest_mode={first_mode:.12g} "
          f"linear={loads[6][0]/loads[3][0]:.12g} "
          f"isolated={loads[6][1]/loads[3][1]:.12g}")


if __name__ == "__main__":
    archive = json.loads(Path("docs/_measurements/held_mouth_pair_field.json").read_text())
    for run in archive["runs"]:
        for name, digest in run.get("fields", {}).items():
            path = Path(name)
            if path.name != "d0.npz":
                continue
            if sha256(path.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Archived isolated field has changed: {path}")
            print(path.parent.parent.name)
            compare(path)
