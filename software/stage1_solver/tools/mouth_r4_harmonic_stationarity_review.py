"""Independent first variation on every unforced melted-label degree of freedom."""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from mouth_cp10_label_gauge import gradient_basis, profile
from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read

VERTICES = np.array(list(product((0, 1), repeat=3)))


def field_work(body, axes, chi, label):
    """Assemble work and weak reactions cellwise in the physical metric."""
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    cells = np.argwhere(np.max(np.lib.stride_tricks.sliding_window_view(
        chi, (2,2,2)), axis=(-3,-2,-1)) == 0.)
    idx = cells[:,None,:] + VERTICES[None,:,:]
    nodes = tuple(idx[:,:,k] for k in range(3))
    displacement = (label - grid)[nodes]
    widths = np.column_stack([np.diff(a)[cells[:,k]] for k,a in enumerate(axes)])
    origin = grid[tuple(cells.T)]
    energy = 0.
    local_force = np.zeros_like(displacement)
    for point in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)), repeat=3):
        x,t,angle = (origin + widths*point).T
        radius, radius_t = profile(body.sh, t, np.zeros(6))
        u = (t-.06)/.56
        bump = np.where((u>0)&(u<1), 16*u**2*(1-u)**2, 0.)
        bump_t = np.where((u>0)&(u<1), 32*u*(1-u)*(1-2*u)/.56, 0.)
        p = (3*np.cos(angle)**2 - 1)/2
        b = body.epsilon
        r = x*radius + b*x**3*bump*p
        dr_dx = radius + 3*b*x**2*bump*p
        dr_dt = x*radius_t + b*x**3*bump_t*p
        dr_da = -3*b*x**3*bump*np.sin(angle)*np.cos(angle)
        metric = np.zeros((len(cells),3,3))
        metric[:,0,0] = dr_dx
        metric[:,0,1] = dr_dt
        metric[:,0,2] = dr_da
        metric[:,1,1] = body.sh.length
        metric[:,2,2] = r
        grad = (gradient_basis(point)[None,:,:]/widths[:,None,:]) @ np.linalg.inv(metric)
        grad_u = displacement.transpose(0,2,1) @ grad
        dv = np.pi/2*np.prod(widths,axis=1)*dr_dx*r**2*body.sh.length*np.sin(angle)
        energy += .5*np.einsum('e,ecd,ecd->',dv,grad_u,grad_u)
        local_force += (dv[:,None,None]*grad) @ grad_u.transpose(0,2,1)

    force = np.zeros_like(label)
    np.add.at(force,nodes,local_force)
    # Q1 order is nonnegative: any positive corner marks the whole cell as
    # active for protecting labels. Intersect that stencil with the registry.
    stencil = np.zeros_like(chi,dtype=bool)
    active = np.max(np.lib.stride_tricks.sliding_window_view(
        chi, (2,2,2)),axis=(-3,-2,-1))>0
    for corner in VERTICES:
        stencil[tuple(slice(k,k+n) for k,n in zip(corner,active.shape))] |= active
    residual = []
    for j in range(3):
        ids = np.arange(chi.size).reshape(chi.shape)
        fixed = stencil.copy()
        fixed[:,[0,-1],:] = True
        fixed[-1,:,:] = True
        if j == 0:
            fixed[0,:,:] = True
        elif j == 1:
            ids[0,:,:] = ids[0,:,:1]
        else:
            fixed[:,:,[0,-1]] = True
            axis_free = ~fixed[0]
            ids[0][axis_free] = ids[1][axis_free]
        dof_force = np.bincount(ids.ravel(),weights=force[...,j].ravel(),
                                minlength=chi.size)
        free = np.unique(ids[~fixed])
        free = np.setdiff1d(free,np.unique(ids[fixed]))
        residual.append(float(np.max(np.abs(dof_force[free]))))
    return energy,residual


def review(root):
    body, checkpoint, axes = get_cp10(root)
    run = root/"r4_eulerian_stage2_01"
    output = {}
    for case in ("primary","secondary"):
        guard = run/(case+"_guard")
        endpoint = read(run/case/"label_stall.json")
        path = run/case/endpoint['file']
        if digest(path)!=endpoint["sha256"]:
            raise ValueError(f"Changed saved labels at {case}")
        with np.load(path, allow_pickle=False) as st:
            if str(st["parent_sha256"])!=CHECKPOINT_SHA or not np.array_equal(st["x"],checkpoint["x"]):
                raise ValueError("Requires the held R4 checkpoint-10 chart")
            order = st["field"][2,body.ids].transpose(0,2,1)
            energy, residual = field_work(body,axes,order,st["labels"])
            epsilon = float(st["epsilon"])
        rows = [json.loads(x) for x in (guard/"stdout").read_bytes().splitlines()]
        steps = [r for r in rows if r.get("event")=="label_accepted"]
        normals = np.array([r["max_unit_harmonic_residual"] for r in steps])
        small = [r for r in steps if r["alpha"]<1.e-4]
        first = small[0]["iteration"]
        tail = [r for r in steps if r["iteration"]>=first]
        native = endpoint["label"]
        output[case] = dict(
            saved_state_sha256=digest(path), epsilon=epsilon,
            independent_unit_gauge_energy=energy, unit_gauge_energy_difference=energy-native["gauge_energy"]/epsilon,
            independent_melted_residual=residual,
            maximum_observed_unit_melted_residual_over_accepted_steps=float(normals.max()),
            total_accepted=len(steps), accepted_below_1e_minus4=len(small),
            first_small_alpha=small[0]["alpha"], first_small_iteration=first,
            final_uninterrupted_subthreshold_steps=len(list(take_small(steps))),
            final_alpha=steps[-1]["alpha"], final_active_force=steps[-1]["max_free_label_force"],
            active_force_at_first_small=tail[0]["max_free_label_force"],
            terminal_melted_residual=native["max_unit_harmonic_residual"])
        if abs(energy - native['gauge_energy']/epsilon)>1.e-11 or max(residual)>1.e-10:
            raise ArithmeticError(output[case])
    return output


def take_small(steps):
    for step in reversed(steps):
        if step['alpha']>=1.e-4:
            break
        yield step


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scratch",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    findings = review(args.scratch)
    print(json.dumps(findings,allow_nan=False),flush=True)
    with args.output.open("x") as stream:
        json.dump(findings,stream,indent=2,allow_nan=False)
        stream.write("\n")
