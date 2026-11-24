"""Read the original all-oriented radius-four state with its pinned source."""
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

if len(sys.argv) != 2:
    raise SystemExit('Pass the directory containing the independent reference and candidate pin')
area = Path(sys.argv[1]).resolve()
root = area.parents[2]
run = root / "_scratch/s11c/r4_continuous_coarse01"
source = run / "pin/software/stage1_solver/tools"
sys.path.insert(0, str(source))
sys.path.append(str(root / "software/stage1_solver/tools"))

from mouth_breathing_mobile_patch import KAPPA, WELL, ORDER_INERTIA
from mouth_r4_fold_constitutive import assemble_trial


def currents_and_order(trial, u, inner):
    m = trial.model
    value = m.V @ u[2]
    grad_c = np.stack([D @ u[2] for D in m.G], axis=1)
    grad_p = np.stack([D @ u[1] for D in m.G], axis=1)
    n = (m.D[:, 4] * (m.V @ u[0]))**2
    den = n + ORDER_INERTIA * np.sum(grad_c**2, axis=1)
    drag = ORDER_INERTIA * np.sum(grad_c * grad_p, axis=1) / den
    velocity = grad_p - drag[:, None] * grad_c
    j = n[:, None] * velocity
    work = KAPPA * grad_c - (n * drag)[:, None] * velocity
    source = 2 * WELL * value * (1 - value) * (1 - 2 * value)
    order = inner + trial.ext + m.V.T @ (m.hq * source)
    for axis, D in enumerate(m.G):
        order += D.T @ (m.hq * work[:, axis])
    flux = {}
    for label, nodes, sign in (
            ("feed_current", m.feed, -1),
            ("bulk_current", m.bulk, 1),
            ("weak_disk_current", m.ids[1, 0][:, :, -1], 1)):
        lifting = np.zeros(m.size)
        lifting[nodes] = sign
        domain = m.D[:, 7] == -1 if label == "weak_disk_current" else slice(None)
        flux[label] = sum(
            np.sum(m.hq[domain] * j[domain, axis] * (D @ lifting)[domain],
                   dtype=np.longdouble)
            for axis, D in enumerate(m.G))
    return trial.pull.T @ order, {k: float(v) for k, v in flux.items()}


def main():
    data = json.loads((root / "_scratch/s11c/r4_fold_gate01/decision.json").read_text())
    cp = run / "coarse/state_00006.npz"
    assert sha256(cp.read_bytes()).hexdigest() == data["sha256"]
    assert cp == Path(data['checkpoint'])
    ref = data["latest"]
    t, q, _, u, _ = assemble_trial(cp, run)
    assert Path(sys.modules["mouth_free_shape_solid"].__file__).parent == source
    b = t.b
    energy, gradient, inner, _ = b.evaluate(q, u[2])
    full_inner = t.pull.T @ inner
    c = u[2, t.first]
    m = t.model
    m.material_load = inner + t.ext
    total_energy, field_grad = m.evaluate(u)
    raw, current = currents_and_order(t, u, inner)
    raw_error = float(np.max(abs(raw - t.pull.T @ field_grad[2])))
    ix, iz, ia = json.loads((root / "research/pde_ledger_v3/_measurements/mouth_r4_fold_constitutive.json").read_text())["ray"][0]["field_cell_x_t_angle"]
    ids = np.unique(t.number[b.ids[ix:ix+2, ia:ia+2, iz:iz+2]])
    assert len(ids) == 8 and np.all(t.free[ids]) and np.all(c[ids] == 0.)
    # An independent admissible order excitation on all corners of the
    # zero-order cell containing the recorded material extremum.
    h = 0.0001
    probe = u[2].copy()
    probe[np.isin(t.number, ids)] += h
    ew = b.evaluate(q, probe, load=False, gradient=False)[0]
    order_secant = (ew - energy) / h
    order_derivative = full_inner[ids].sum()
    # Use the most heavily loaded generalized material coordinate for a
    # signed local displacement. The full native energy is independently
    # re-evaluated on either side and the step is halved.
    iq = np.abs(gradient).argmax()
    slopes = []
    for step in (5e-7, 2.5e-7):
        dq = np.zeros_like(q)
        dq[iq] = step
        ep = b.evaluate(q+dq, u[2], load=False, gradient=False)[0]
        em = b.evaluate(q-dq, u[2], load=False, gradient=False)[0]
        slopes.append((ep-em)/(2*step))
    current_errors = {k: v - ref["current"][k] for k,v in current.items()}
    results = dict(checkpoint_sha256=data["sha256"], solid_energy=energy, total_energy=total_energy,
                   archived_solid_energy_error=energy - ref["field"]["solid_energy"],
                   archived_total_energy_error=total_energy-ref["energy"],
                   independent_raw_order_difference=raw_error,
                   lower_bound_nodes=int(np.sum((c==0)&t.free)),
                   folded_corner_eight_zero_order_nodes=ids.tolist(),
                   folded_corner_inner_order_load=full_inner[ids].tolist(),
                   folded_corner_obstacle_reaction=raw[ids].tolist(),
                   folded_corner_order_excitation=dict(slope=float(order_derivative), secant=float(order_secant)),
                   independent_current=current, archived_current_difference=current_errors,
                   material_probe=dict(index=int(iq), derivative=float(gradient[iq]), secants=slopes))
    (area/"reference.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    print(json.dumps(results, allow_nan=False), flush=True)
    assert max(abs(v) for v in current_errors.values()) < 2e-8
    assert raw_error < 2e-8 and abs(total_energy-ref["energy"]) < 2e-8
    assert abs(energy-ref["field"]["solid_energy"]) < 2e-9
    assert max(abs(s - gradient[iq]) for s in slopes) < 1e-3
    assert abs(order_secant-order_derivative) < 1e-5
    np.savez_compressed(area/"reference.npz", gq=gradient, inner=inner,
                        joined_inner=full_inner, raw=raw, total=total_energy,
                        solid=energy, field_gradient=field_grad,
                        shape_force=ref["shape_force"])


if __name__ == "__main__":
    main()
