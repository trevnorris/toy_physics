"""Compare the old and candidate production evaluations on the frozen R4 checkpoint 10.

Run the reference and candidate separately with the original simulation pin first
on sys.path. Neither phase rematches fields or advances material/shape coordinates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


CHECKPOINT_SHA = "51727c2d48eea199c29bacbfd4bf07b7b11a3e997a737ebcc8c2f69925cd112b"
PORTS = ("feed_current", "bulk_current", "weak_disk_current", "incident_slice_current",
         "receiving_slice_current", "positive_slice_current", "port_error", "disk_minus_feed_current")


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def emit(path, data):
    text = json.dumps(data, indent=2, allow_nan=False)
    path.with_suffix(".partial").write_text(text+"\n")
    path.with_suffix(".partial").replace(path)


def acquired(run, candidate):
    pin = run/"pin/software/stage1_solver/tools"
    sys.path.insert(0, str(pin))
    if candidate is not None:
        # The immutable overlay contains the guard and its orientation helper.
        # Every other producer module still comes from the verified original pin.
        allowed = {"mouth_active_prism_orientation.py", "mouth_free_shape_solid.py"}
        for file in candidate.glob("*.py"):
            old = pin/file.name
            if old.exists() and file.name not in allowed and digest(old)!=digest(file):
                raise ValueError(f"Unexpected replacement of producer module {old.name}")
        sys.path.insert(0, str(candidate))
    from mouth_continuous_shape_relax import ContinuousRelease
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_melted_collar_survival import MeltTrial
    from mouth_exchange_carrier import BIAS

    for name, checksum in read(run/"sources.json").items():
        if digest(run/"pin"/name) != checksum:
            raise ValueError(f"Producer source changed: {name}")
    case = run/"coarse"
    checkpoint = case/"state_00010.npz"
    record = read(case/"latest.json")
    if digest(checkpoint) != CHECKPOINT_SHA or record["sha256"] != CHECKPOINT_SHA or record["iteration"] != 10:
        raise ValueError("Expected frozen checkpoint 10")
    previous = Path(read(case/"invocation.json")["command"][-3])/"coarse"
    invocation = read(previous/"invocation.json")["command"]
    original = Path(invocation[-5])/"coarse"
    proof = read(Path(invocation[-3]))["results"]["coarse"]
    state = Path(read(original/"invocation.json")["command"][5])
    assert digest(state) == proof["preparation_sha256"]
    assert digest(original/"endpoint.json") == proof["endpoint_sha256"]
    assert BIAS == .4510084282857874
    with np.load(checkpoint, allow_pickle=False) as d:
        assert str(d["state_sha"])==digest(state) and int(d["radius"])==4 and int(d["winding"])==1
        x,q,u,solid = (d[k].copy() for k in ("x","q","field","solid"))

    trial_class = MeltTrial
    if candidate is not None:
        from mouth_exact_melt_trial import ExactMeltTrial
        trial_class = ExactMeltTrial
        if Path(sys.modules[trial_class.__module__].__file__).resolve().parent != candidate:
            raise ValueError("Candidate import path does not identify supplied code")
    if Path(sys.modules[MeltTrial.__module__].__file__).resolve().parent != pin:
        raise ValueError("Producer import bypassed its frozen pin")
    owner = SimpleNamespace(base=Front(state,GRIDS["double"],1), args=SimpleNamespace(state=state),
            discretization=proof["solid_discretization"],reference_depth=read(original/"endpoint.json")["depth"],
            trial_type=trial_class)
    return lambda z: ContinuousRelease.build(owner,z,False), x,q,u,solid,record


def folded_corner_control(t,q,u):
    # Whole-prism readback identified corner (j=5, i=18, k=1). Give every Q1
    # basis nonzero there the same value, preserving joined spatial labels.
    from mouth_active_prism_orientation import _Support, _exact_array
    b = t.b
    offset = b.unpack(q)
    chart = _exact_array([b.x[18],b.t[5],b.a[1]])+_exact_array(offset[5,18,1])
    pos = np.asarray(chart,dtype=float)
    axes = (b.buffers["field_x"], b.buffers["field_t"], b.buffers["field_a"])
    cell = [np.searchsorted(z,p,side="right")-1 for z,p in zip(axes,pos)]
    i,j,k = cell
    nodes = np.unique(t.number[b.ids[i:i+2,k:k+2,j:j+2]])
    test = u.copy()
    test[2,np.isin(t.number,nodes)]=.001
    actual = _Support(b,test[2]).at(chart)
    if float(actual) != .001:
        raise ArithmeticError("Control did not put chi=0.001 at folded corner")
    # The same production feasibility call gates a shape/fixed-work evaluation.
    cert = t.certificate(q,test[2],False)
    if cert["status"] != "active_fold" or cert["all_active_cells_oriented"]:
        raise ArithmeticError(f"Folded-corner negative control failed: {cert}")
    from mouth_free_shape_solid import InvalidPlacement
    try:
        t.fixed(q,test)
    except (InvalidPlacement, ArithmeticError) as exc:
        stopped = {"exception":type(exc).__name__,"message":str(exc)}
    else:
        raise ArithmeticError("Production signed work accepted the actual inverted ordered sample")
    return dict(corner_t_x_polar=[5,18,1],chi=float(actual),landed_x_t_theta=pos.tolist(),
                node_ids=nodes.tolist(),status=cert["status"],witness=cert["witness"],fixed_rejection=stopped)


def evaluate(run, overlay, area):
    from_name = "new" if overlay is not None else "old"
    build,x,q,u,solid,rec = acquired(run,overlay)
    from mouth_whole_r4_relax import STEP
    t = build(x)
    if not np.array_equal(t.b.unpack(q),solid):
        raise ArithmeticError("Labels changed during checkpoint acquisition")
    inner = t.crystal_load_for_match(q,u)
    load = inner+t.ext
    t.model.material_load = load
    energy,g = t.model.evaluate(u)
    stationarity = t.check_settled_fields(q,u,g)
    e,gq,_,_ = t.b.evaluate(q,u[2],load=False)
    flow = t.measure(u,load)
    values = dict(energy=np.array(energy),solid_energy=np.array(e),solid_force=-gq,
                  order_deposition=inner,field_gradient=g,ports=np.array([flow[k] for k in PORTS]))
    control = None if overlay is None else folded_corner_control(t,q,u)
    del t
    work = []
    for mode,unit in enumerate(np.eye(6)):
        energies = [build(x+sign*STEP*unit).fixed(q,u) for sign in (-1,1)]
        work.append(energies)
        print(json.dumps(dict(side=from_name,mode=mode,signed_energies=energies)),flush=True)
    values["signed_energies"]=np.array(work)
    values["shape_force"]=(values["signed_energies"][:,0]-values["signed_energies"][:,1])/(2*STEP)
    archive=dict(energy=rec["energy"],solid_energy=rec["field"]["solid_energy"],
                 shape_force=rec["shape_force"],ports=[rec["current"][k] for k in PORTS])
    archive_error={k:float(np.max(abs(values[k]-v))) for k,v in archive.items()}
    if any(archive_error.values()):
        raise ArithmeticError(f"Evaluations differ from saved readouts: {archive_error}")
    np.savez(area/(from_name+".npz"),**values)
    result=dict(checkpoint=CHECKPOINT_SHA,step=STEP, archived_differences=archive_error,
                producer_manifest_sha=digest(run/"sources.json"),
                producer_source_sha=digest(run/"pin/software/stage1_solver/tools/mouth_melted_collar_survival.py"),
                total_energy=float(energy),solid_energy=float(e),n_material_forces=len(gq),
                signed_shape_forces=values["shape_force"].tolist(), ports={k:flow[k] for k in PORTS},
                stationarity=stationarity["raw_residual"],negative_control=control)
    if overlay is not None:
        with np.load(area/"old.npz",allow_pickle=False) as f:
            if set(f.files)!=set(values):
                raise ValueError("Reference lacks required observables")
            difference={k:float(np.max(abs(values[k]-f[k]))) for k in values}
            bits={k:np.asarray(values[k]).tobytes()==f[k].tobytes() for k in values}
        result.update(exact_differences=difference,bitwise_equal=bits,
                      exact_orientation_sha=digest(overlay/"mouth_active_prism_orientation.py"),
                      native_solid_sha=digest(overlay/"mouth_free_shape_solid.cpp"),
                      native_interface_sha=digest(overlay/"mouth_free_shape_solid.py"),
                      candidate_sha=digest(overlay/"mouth_exact_melt_trial.py"))
        if not all(bits.values()):
            emit(area/"difference.json",result)
            raise ArithmeticError(f"Reference and candidate differ: {difference}")
    emit(area/(from_name+".json"),result)
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("run",type=Path)
    p.add_argument("area",type=Path)
    p.add_argument("--candidate",type=Path)
    args=p.parse_args()
    args.area.mkdir(exist_ok=True,parents=True)
    evaluate(args.run.resolve(),None if args.candidate is None else args.candidate.resolve(),args.area.resolve())
