"""Sum the signed field and spatial-crystal variations of the same six CP10 virtual motions.

The wave-free coherent carrier is common to the legacy and two-metric lifts.
Reconstruct the total central difference componentwise to avoid cancellation
against the order energy of the far bath.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def reduce(fields, projected, legacy, output, registry=None):
    field, spatial, archive = map(read, (fields, projected, legacy))
    old = archive['states']['cp10']['signed_crystal_work']
    if any(k != CHECKPOINT_SHA for k in (field["checkpoint"], spatial["checkpoint"], old["state_sha256"])):
        raise ValueError("Incompatible checkpoint for signed work")
    material, = [row for row in spatial["levels"] if "analytic" in row]
    if material['order_positive_nonpositive_Jinv_samples'] or material['modes'] != old['modes']:
        raise ValueError("Projected material state/modes failed validation")
    inputs = [fields, projected, legacy]
    if registry is not None:
        split = read(registry)
        if split['labels'] != spatial['labels'] or split['checkpoint'] != CHECKPOINT_SHA:
            raise ValueError('Registry-split inputs differ from the projected state')
        fixed = split['fixed_label_split_registry']
        es = np.asarray(fixed['signed_energies'])
        material = dict(material, analytic=fixed['analytic'],
                        forces_h=(es[:,0]-es[:,1])/2e-4,
                        forces_half=(es[:,2]-es[:,3])/1e-4)
        inputs.append(registry)
    euler = np.array(material["analytic"])
    modes = field["modes"]
    assert [r['mode'] for r in modes] == list(range(6))
    bath = np.array([r["analytic"] for r in modes])
    conjugate = euler+bath
    cp10 = np.array(old['saved_total'])
    old_solid = np.array(old['crystal_at_matched_kinematics'])
    fd = {}
    field_fd = {}
    for step, key in ((1.e-4, "forces_h"), (5.e-5, "forces_half")):
        field_fd[str(step)] = np.array([m["central"][str(step)]["force"] for m in modes])
        fd[str(step)] = field_fd[str(step)]+np.array(material[key])
    rich = (4*fd['5e-05']-fd['0.0001'])/3
    full = np.array(material['forces_h'])
    change = full-old_solid
    field_gap = field_fd["0.0001"]-(cp10-old_solid)
    if np.max(abs((change+field_gap)-(fd['0.0001']-cp10))) > 5e-15:
        raise ArithmeticError("Matched virtual-work decomposition did not close")
    result = dict(checkpoint=CHECKPOINT_SHA, labels=spatial['labels'],
        source_sha256={p.name:digest(p) for p in inputs},
        scope="Held projected CP10 spatial field at zero gauge; signed derivatives transport nodal labels and order on the moving grid. Richardson uses both central differences.",
        modes=material['modes'], native_saved=cp10.tolist(), analytic_total=conjugate.tolist(),
        analytic_minus_saved=(conjugate-cp10).tolist(), field_analytic=bath.tolist(),
        inner_crystal_analytic=euler.tolist(), degree=material['degree'],
        central_total={k:v.tolist() for k,v in fd.items()}, richardson=rich.tolist(),
        analytic_minus_richardson=(conjugate-rich).tolist(),
        matched_full_step_crystal_shift=change.tolist(),
        matched_full_step_remaining_fields_delta=field_gap.tolist(),
        absolute_field_energy_delta_against_saved=field['static_delta'])
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("fields","projected","legacy","output"):
        parser.add_argument(arg,type=Path)
    parser.add_argument("--registry",type=Path,help="Use signed work with explicit collar reference-triangle splits")
    args=parser.parse_args()
    reduce(args.fields,args.projected,args.legacy,args.output,args.registry)
