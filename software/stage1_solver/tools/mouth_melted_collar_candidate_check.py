"""Probe the candidate on the accepted R4 state and two controlled label displacements."""
from fractions import Fraction
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
src = area / "pin/software/stage1_solver/tools"
sys.path[:0] = [str(src),str(run/"pin/software/stage1_solver/tools"),str(root/"software/stage1_solver/tools")]
from mouth_active_prism_orientation import active_orientation
from mouth_free_shape_solid import InvalidPlacement
from mouth_r4_fold_constitutive import assemble_trial, determinant, Prism


def get_raised_field(trial, u, ids):
    out = u[2].copy()
    out[np.isin(trial.number, np.unique(trial.number[ids]))] = .001
    return out


def main():
    gate = json.loads((root/"_scratch/s11c/r4_fold_gate01/decision.json").read_text())
    cp = Path(gate["checkpoint"])
    assert sha256(cp.read_bytes()).hexdigest() == gate["sha256"]
    t, q, offsets, u, _ = assemble_trial(cp, run)
    b = t.b
    b.melted_labels = True
    m = t.model
    e, g, inner, _ = b.evaluate(q,u[2])
    m.material_load = inner+t.ext
    total, grad = m.evaluate(u)
    with np.load(area/"reference.npz") as baseline:
        errors = {k: float(np.max(abs(v-baseline[k]))) for k,v in
                  (("gq",g),("inner",inner),("joined_inner",t.pull.T@inner),("raw",t.pull.T@grad[2]),
                   ("total",total),("solid",e),("field_gradient",grad))}
    print("PINNED BASELINE", errors, flush=True)
    # The production placements can be read exactly as binary nodal offsets.
    # Cross the recorded corner determinant in the depth displacement at
    # that corner. Try both signs and select by the independently calculated
    # rational determinant; the trial is independent of the new certificate.
    ix=(5,17,0)
    loc=gate["latest"]["orientation"]["extremizers"]["min"]
    assert tuple(loc["cell_t_x_angle"]) == ix and loc["triangle"] == 0
    key=(5,18,1,1)
    samples=[]
    for sign in (-1,1):
        near=offsets.copy()
        near[key]+=sign*1e-8
        p=Prism(b,near,ix,0)
        det=determinant(p.chart([0,1,0],1.,Fraction,False))
        samples.append((det,near))
    det, near = min(samples,key=lambda pair: pair[0])
    assert det < 0, float(det)
    tiny = b.pack(near)
    cert=active_orientation(b,near,u[2])
    _,_,_,low=b.lower_bound_load(tiny,u[2])
    ids=np.array(cert["potential_activation_nodes"],int)
    inside=json.loads((root/"research/pde_ledger_v3/_measurements/mouth_r4_fold_constitutive.json").read_text())["ray"][0]["field_cell_x_t_angle"]
    i,j,k=inside
    eight=b.ids[i:i+2,k:k+2,j:j+2].ravel()
    raised = get_raised_field(t,u,eight)
    reactivated = active_orientation(b,near,raised)
    subpoint=dict(exact_corner_Jchart=float(det), original_zero_cert=cert["status"],
                  first_positive_Q1_cert=reactivated["status"],
                  candidate_witness=reactivated["witness"],
                  activation_nodes=len(ids),contains_eight=bool(np.isin(eight,ids).all()),
                  invalid_quadrature=low["nonpositive_J_points"])
    print("SUBQUADRATURE",json.dumps(subpoint,allow_nan=False),flush=True)
    assert cert["all_active_cells_oriented"] and np.isin(eight,ids).all()
    assert reactivated["status"]=="active_fold" and low["nonpositive_J_points"]==0
    # A completely separate interior node, many full cells beyond the last
    # positive order node. Move it by one and a half nominal depth cells.
    deep=offsets.copy()
    deep[18,10,3,1] += .065
    qd=b.pack(deep)
    ed,gd,_,_=b.evaluate(qd,u[2],load=False)
    _,_,lower,info=b.lower_bound_load(qd,u[2])
    deepcert=active_orientation(b,deep,u[2])
    native_ids=np.array(info["undefined_order_nodes"],int)
    geometric_ids=np.array(deepcert["potential_activation_nodes"],int)
    # Cross-check that the geometric certificate encloses every queried
    # invalid native site (the real-space chart is fully retained).
    overlap=bool(np.isin(native_ids,geometric_ids).all())
    m.material_load=lower+t.ext
    new_total,field=m.evaluate(u)
    affected=np.unique(t.number[geometric_ids])
    joined=t.pull.T@field[2]
    free=affected[t.free[affected]]
    measured=m.measure(u)
    reference=json.loads((area/"reference.json").read_text())
    bad=get_raised_field(t,u,native_ids)
    try:
        b.evaluate(qd,bad,load=False)
        raised_native=False
    except InvalidPlacement:
        raised_native=True
    deep_probe=dict(old_material_energy=e, material_delta=ed-e, force_delta=float(np.max(abs(gd-g))),
                    total_delta=new_total-total, native_invalid=info["nonpositive_J_points"],
                    minimum_J=info["min_nonpositive_J"], native_nodes=len(native_ids),
                    geometric_nodes=len(geometric_ids), geometric_coverage=overlap,
                    cert=deepcert["status"], raised_native_rejects=raised_native,
                    raw_reaction_min=float(joined[free].min()), zero_order=bool(np.all(u[2,t.first][affected]==0)),
                    joined_load_change=float(np.max(abs(t.pull.T@(lower-inner)))),
                    current_error={k:float(measured[k]-v) for k,v in reference["independent_current"].items()})
    print("DEEP FOLD",json.dumps(deep_probe,allow_nan=False),flush=True)
    output=dict(candidate_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in src.iterdir() if p.is_file()},
                reference_checkpoint=gate["sha256"], reference_differences=errors, subquadrature=subpoint, deep=deep_probe)
    (area/"candidate.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    assert max(errors.values()) < 1e-11
    assert deep_probe["native_invalid"] > 0 and deep_probe["cert"]=="certified" and overlap and raised_native
    assert deep_probe["zero_order"] and deep_probe["raw_reaction_min"]>=0 and deep_probe["joined_load_change"]>0.
    assert max(abs(deep_probe[k]) for k in ("material_delta","force_delta","total_delta")) < 1e-11
    assert max(abs(v) for v in deep_probe["current_error"].values()) < 2e-8


if __name__ == "__main__":
    main()
