"""Reconstruct the R4 fixed-shape stopping witness from accepted labels, without a descent step."""
import argparse
import hashlib
from itertools import product
from pathlib import Path
import json

import mpmath as mp
import numpy as np

from mouth_cp10_label_condition import physical_frame, reference_frame
from mouth_cp10_spatial_crystal import load_context
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_fold_order import all_low_after
from mouth_supercore_shape import preparation


def as_mp(array):
    return mp.matrix(np.asarray(array).astype(float).tolist())


def det_mp(matrix):
    return (matrix[0, 0]*(matrix[1, 1]*matrix[2, 2]-matrix[1, 2]*matrix[2, 1])
            - matrix[0, 1]*(matrix[1, 0]*matrix[2, 2]-matrix[1, 2]*matrix[2, 0])
            + matrix[0, 2]*(matrix[1, 0]*matrix[2, 1]-matrix[1, 1]*matrix[2, 0]))


def run(scratch, case):
    mp.mp.dps = 65
    area = scratch / "r4_eulerian_stage2_01"
    guard = area / (case + "_guard")
    out = read(guard / "outcome.json")
    assert out["exitCode"] == out["childOutcome"]["exitCode"] == 1 and out["limitsVerified"]
    assert not out["childOutcome"]["guardReason"]
    assert int(read(guard / "effective-limits.json")["memory.swap.max"]) == 0
    src = read(area / "sources.json")
    assert all(digest(area/"pin/software/stage1_solver/tools"/key)==value for key, value in src.items())
    stall = read(area/case/"label_stall.json")
    assert digest(area/case/stall["file"]) == stall["sha256"]
    with np.load(area/case/stall["file"]) as f:
        assert str(f["parent_sha256"]) == CHECKPOINT_SHA
        labels, fields, x = (f[key].copy() for key in ("labels","field","x"))
    with np.load(scratch/"r4_melted_survival01/coarse/state_00010.npz") as cp:
        assert np.array_equal(cp["x"], x) and int(cp["radius"]) == 4 and int(cp["winding"]) == 1
    log = [json.loads(row) for row in (guard/"stdout").read_bytes().splitlines()]
    cuts = [row for row in log if row.get("event") == "label_accepted"]
    start = next(row for row in log if row.get("event") == "label_start")
    data = stall["label"]
    assert len(cuts) == data["accepted_steps"] and all(r["energy_delta"] < 0 for r in cuts)
    assert abs((start["elastic_energy"]+start["gauge_energy"]+sum(row["energy_delta"] for row in cuts))
               -(data["elastic_energy"]+data["gauge_energy"])) < 1e-13
    reject = data["last_rejection"]
    if "code" in reject:
        assert reject["code"] == 2 and reject["chart_bound"] < 0
        local = reject["local"]
    else:
        assert reject["reason"] == "Nonpositive or roundoff-unresolved chart determinant"
        assert 0 < reject["chart_det"] <= reject["roundoff_margin"] and reject["chart_lower"] < 0
        local = reject["chart_fractions"]
    cell = np.array(reject["cell"])
    bits = np.array(local, dtype=int)
    assert np.array_equal(local, bits)
    index = cell + bits
    _, _, body, _, _, _, _, _ = load_context(scratch, scratch/"r4_cp10_labels01/labels.npz")
    axes = [body.buffers["field_"+s] for s in "xta"]
    order = np.ascontiguousarray(fields[2,body.ids].transpose(0,2,1))
    node = tuple(index)
    pos = np.array([ax[v] for ax,v in zip(axes,node)])
    xi = labels[node]
    B, lever = physical_frame(body, pos[None])
    C, natural = reference_frame(body, preparation(4)[0], xi[None])
    B, C = as_mp(B[0]), as_mp(C[0])
    w = float(body.sh.length*pos[1])
    r = float(B[2,2])
    incident = []
    target = None
    for offset in product((0,1), repeat=3):
        lower = index-np.array(offset)
        if any(c<0 or c+1>=len(ax) for c,ax in zip(lower,axes)):
            continue
        support = tuple(slice(i,i+2) for i in lower)
        # At a vertex the three columns are the incident edges; use the stored
        # binary-64 endpoints exactly when subtracting to resolve the tiny volume.
        G = mp.matrix(3,3)
        for component, ax in enumerate(axes):
            a = index.copy()
            a[component] = lower[component]
            b = a.copy()
            b[component] += 1
            span = mp.mpf(float(ax[b[component]])) - mp.mpf(float(ax[a[component]]))
            for coord in range(3):
                G[coord,component] = (mp.mpf(float(labels[tuple(b)][coord]))-
                                      mp.mpf(float(labels[tuple(a)][coord])))/span
        inverse = C*G*(B**-1)
        hoop = mp.mpf(float(natural[0])) / mp.mpf(float(lever[0]))
        inv_j = mp.sign(det_mp(G))*abs(det_mp(inverse)*hoop)
        assert inv_j > 0
        F = inverse**-1
        J = 1/inv_j
        energy = (sum(v*v for v in F)+1/hoop**2-4)/2 - mp.log(J) + (J-1)**2/4
        sample = dict(cell=lower.tolist(), max_order=float(order[support].max()),
                      chart_det=float(det_mp(G)), Jinv=float(inv_j), W=float(energy),
                      grad_xi_frobenius=float(mp.mnorm(G*(B**-1),"F")))
        incident.append(sample)
        if np.array_equal(lower, cell):
            target = sample
    assert target is not None and order[node] == 0.0
    if "chi" in reject:
        assert order[node] == reject["chi"]
    angles = np.array([row["alpha"] for row in cuts])
    return dict(case=case, epsilon=data["epsilon"], saved_sha256=stall["sha256"],
                accepted_steps=len(cuts), strict_energy_change=sum(row["energy_delta"] for row in cuts),
                initial_free_force=start["max_free_label_force"], final_free_force=data["max_free_label_force"],
                accepted_last_alpha=angles[-1], accepted_min_alpha=angles.min(),
                accepted_last_energy_change=cuts[-1]["energy_delta"], trial_alpha=data["trial_alpha"],
                # The saved trial_alpha is already halved after this rejection.
                last_rejected_alpha=2*data["trial_alpha"], rejected=reject,
                chart_coordinates=pos.tolist(), material_label=xi.tolist(),
                physical_r_w_theta=[r,w,pos[2]], chi_at_vertex=float(order[node]),
                accepted_corner=target, incident=incident, active_chart_lower=data["min_chart_Jinv_active_bound"],
                sampled_diagnostic=data["worst_global"], harmonic_free_residual=data["max_unit_harmonic_residual"],
                gauge_energy=data["gauge_energy"], energy_elastic=data["elastic_energy"],
                sampled_front=all_low_after(axes[1]*body.sh.length,order.transpose(0,2,1),0.),
                shape_coordinates=x.tolist(), missing_checkpoint11=not (area/case/"checkpoint_00011.json").exists())


def current_budget(scratch):
    """Contract the saved phase current with three independent boundary lifts."""
    from mouth_breathing_mobile_patch import ORDER_INERTIA
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS, RimCarrier

    area = scratch / "r4_eulerian_stage2_01"
    states = [area/case/"label_stall.npz" for case in ("primary", "secondary")]
    for path in states:
        assert digest(path) == read(path.with_suffix(".json"))["sha256"]
    with np.load(states[0]) as a, np.load(states[1]) as b:
        assert np.array_equal(a["field"], b["field"]) and np.array_equal(a["x"], b["x"])
        assert str(a["parent_sha256"]) == str(b["parent_sha256"]) == CHECKPOINT_SHA
        fields = a["field"]
    pin = area/"pin/software/stage1_solver/tools"
    for source, expected in read(area/"sources.json").items():
        assert digest(pin/source) == expected
    prep, reference, body, *_ = load_context(scratch, scratch/"r4_cp10_labels01/labels.npz")
    base = Front(prep, GRIDS["double"], 1)
    m = RimCarrier(body.sh, base.c, 1, base.grid, polar_shift=body.epsilon, precondition=False,
                   axes=axes_for_depth(body.sh.length, reference))
    phase_grad = np.column_stack([g@fields[1] for g in m.G])
    chi_grad = np.column_stack([g@fields[2] for g in m.G])
    wt, _, _, _, f, _, s, cut, _, surface, _ = m.D.T
    n = (f*(m.V@fields[0]))**2
    coupling = ORDER_INERTIA*np.sum(phase_grad*chi_grad, axis=1)/(n+ORDER_INERTIA*np.sum(chi_grad**2, axis=1))
    velocity = phase_grad-coupling[:,None]*chi_grad
    flux = n[:,None]*velocity
    throat = cut == -1
    values = {}
    for name, ids, sign, mask in [("feed", m.feed, -1, slice(None)), ("bulk", m.bulk, 1, slice(None)),
                                  ("disk", m.ids[1, 0][:, :, -1], 1, throat)]:
        lift = np.zeros(m.size)
        lift[ids] = 1.
        values[name] = float(sign * sum((wt[mask]*(g@lift)[mask])@flux[mask, i] for i,g in enumerate(m.G)))
    inward = surface*np.maximum(flux[:, 2], 0.)
    subsonic = throat & (n>.5) & (np.sum(velocity**2, axis=1) + 1/s**2 < n**4)
    values["subsonic_share"] = float(inward[subsonic].sum()/inward[throat].sum())
    values["field_sha256"] = hashlib.sha256(fields.tobytes()).hexdigest()
    values["max_conservation_residual"] = max(abs(values["feed"] - values[k]) for k in ("bulk","disk"))
    for path in states:
        native = read(path.with_suffix(".json"))["held_field_current"]
        for name, key in [("feed","feed_current"), ("bulk","bulk_current"), ("disk","weak_disk_current"),
                          ("subsonic_share","dense_total_subsonic_share")]:
            assert abs(values[name]-native[key]) < 1e-10
    return values


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scratch", type=Path)
    parser.add_argument("case", choices=("primary","secondary","shared"))
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    row = current_budget(args.scratch) if args.case == "shared" else run(args.scratch, args.case)
    with args.out.open("x") as stream:
        json.dump(row, stream, indent=2, allow_nan=False)
        stream.write("\n")
