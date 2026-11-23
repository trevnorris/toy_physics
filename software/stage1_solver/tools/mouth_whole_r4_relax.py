"""Release the R4 wall and bonded solid under the rematched coherent flow.

The material and shape blocks take energy-lowering steps in the same objective.
The material need not be equilibrated on an intermediate held shape. Its free
load and the independent shape forces must both pass at the final shape.
"""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from mouth_free_shape_solid import InvalidPlacement, Solid
from mouth_joined_front_audit import Front
from mouth_joint_endpoint_audit import audit_endpoint, geometry, orientation, solid_readout
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_mobile_joined_fields import JointFields
from mouth_rim_pinned_comparison import digest, read, validate_guard


STEP = 1.e-4
SHAPE_TOL = 4.e-4
SOLID_TOL = 3.e-5
FIELD_TOL = 2.e-7
COORDINATES = ("log_depth", "scaled_log_exit", "mean_width", "depth_asymmetry",
               "waist_profile", "polar_equatorial")


def log(kind, **values):
    print(json.dumps(dict(kind=kind, **values), allow_nan=False), flush=True)


def atomic_json(path, value):
    tmp = path.with_suffix(".partial")
    tmp.write_text(json.dumps(value, allow_nan=False, indent=2)+"\n")
    tmp.replace(path)


def save_fields(path, **fields):
    tmp = path.with_suffix(".partial.npz")
    np.savez_compressed(tmp, **fields)
    tmp.replace(path)


class Trial:
    def __init__(self, args, base, x, solve_fields):
        self.f = JointFields(base, args.state, x, grid="double", precondition=solve_fields)
        self.b = Solid(args.state, self.f.shape, self.f.polar_shift, args.n, args.k,
                       args.angles, self.f.model, degree=args.degree, subdivisions=args.subdivisions)

    def fixed(self, q, u):
        reference, _, load, _ = self.b.evaluate(q, gradient=False)
        return fsum((self.f.static, reference, self.f.energy(u, load)))

    def match(self, q, initial):
        reference, _, load, _ = self.b.evaluate(q, gradient=False)
        u, info = self.f.solve(initial, load, tolerance=FIELD_TOL)
        e, g, updated_load, _ = self.b.evaluate(q, u[2])
        if not np.array_equal(load, updated_load):
            raise ArithmeticError("Crystal deposition changed when only order amplitudes changed")
        exchange = fsum((e, -reference, -float(load@u[2])))
        if abs(exchange) > 1.e-8:
            raise ArithmeticError(f"Unbalanced material/field work {exchange}")
        energy = fsum((self.f.static, reference, info["energy"]))
        return u, energy, g, load, dict(field=info, exchange=exchange,
                                        receiver=self.f.receiver["residual"])


def inverse_step(gradient, history, metric):
    v = gradient.copy()
    coeff = []
    for dx, dg, inv in reversed(history):
        c = inv*(dx@v)
        coeff.append(c)
        v -= c*dg
    v /= metric
    for (dx, dg, inv), c in zip(history, reversed(coeff)):
        v += dx*(c-inv*(dg@v))
    return v


def infeasible(error):
    # Other failures invalidate the objective and should reach the run record.
    return isinstance(error, InvalidPlacement) or isinstance(error, ValueError) and any(word in str(error).lower() for word in
        ("placement", "nonpositive", "inverting", "orientation", "outside", "folded", "infeasible"))


def relax_material(trial, q, u, history, steps):
    """Advance the material at fixed fields, retaining its secant preconditioner."""
    b = trial.b
    energy, grad, _, _ = b.evaluate(q, u[2], load=False)
    count = 0
    last_failure = None
    for _ in range(steps):
        info = b.residual(q, grad)
        if material_norm(info) < SOLID_TOL:
            break
        quasi = inverse_step(grad, history, b.metric)
        accepted = False
        for direction in (quasi, grad/b.metric):
            if grad@direction <= 0:
                continue
            scale = min(1., .025/max(1.e-20, np.max(abs(direction))))
            for half in range(38):
                step = -scale*2.**-half*direction
                if np.array_equal(q+step, q):
                    break
                if not orientation(trial.f.shape, trial.f.polar_shift, b,
                                   b.unpack(q+step))["all_open_cells_oriented"]:
                    last_failure = "Trial folds a material element"
                    continue
                try:
                    new_e, new_g, _, _ = b.evaluate(q+step, u[2], load=False)
                except ValueError as exc:
                    if not infeasible(exc):
                        raise
                    last_failure = str(exc)
                    continue
                slope = float(grad@step)
                if new_e-energy < 1.e-4*slope:
                    accepted = True
                    break
                # Native energy can be much larger than last-ring work at the
                # force tolerance. Compare two work integrals before acceptance.
                if abs(slope) < 32*np.spacing(max(1., abs(energy))) and new_e <= energy+4*np.spacing(abs(energy)):
                    _, middle, _, _ = b.evaluate(q+step/2, u[2], load=False)
                    edge = .5*((grad+new_g)@step)
                    simpson = (edge+2*(middle@step))/3
                    if simpson < 1.e-4*slope and abs(edge-simpson) < .03*abs(slope):
                        accepted = True
                        break
            if accepted:
                break
        if not accepted:
            last_failure = last_failure or "No energy-lowering material step on either search direction"
            break
        if np.array_equal(new_g,grad):
            last_failure = "Accepted displacement no longer changes the material gradient at machine precision"
            break
        delta = new_g-grad
        product = delta@step
        if product > 1.e-12*np.linalg.norm(step)*np.linalg.norm(delta):
            history.append((step, delta, 1./product))
            del history[:-16]
        q, energy, grad = q+step, new_e, new_g
        count += 1
    return q, count, last_failure


def material_norm(info):
    return max(info["material_force_max"], info["angular_sum_force_max"])


class Release:
    def __init__(self, args):
        self.args = args
        base = args.state
        prior = read(args.proof)["comparison"]
        meta = next(c["input"] for c in prior if c["grid"] == "double")
        for key, p in (("state_sha",base), ("wall_sha",args.wall), ("zero_sha",args.field)):
            if digest(p) != meta[key]:
                raise ValueError(f"Unverified radius-four input {key}")
        validate_guard(args.field.parent)
        self.base = Front(base, GRIDS["double"], 1)
        x, old_d, _, s = read_source(base, args.wall)
        self.x = x
        self.trial = self.build(x, True)
        self.q = self.trial.b.seed(old_d, s)
        with np.load(args.field) as f:
            self.u, self.rest = f["wound"], f["rest"]
        if self.u.shape != self.trial.f.model.seed.shape:
            raise ValueError("Saved coherent field incompatible with reference mesh")
        self.history = []
        self.iteration = 0
        self.shape_metric = None
        self.last = None
        log("input", radius=4, winding=1, coordinates=COORDINATES, carrier="double",
            solid=[3*args.n, args.k, args.angles], solid_degree=args.degree,
            subdivisions=args.subdivisions, held_sha=meta["wall_sha"],
            field_sha=meta["zero_sha"], state_sha=meta["state_sha"], start=x.tolist())

    def build(self, x, solve_fields=False):
        return Trial(self.args, self.base, x, solve_fields)

    def work(self, x, q, u, e, h):
        pairs = []
        for unit in np.eye(len(x)):
            pairs.append([self.build(x+sign*h*unit).fixed(q,u) for sign in (-1,1)])
        points = np.array(pairs)
        return (points[:,1]-points[:,0])/(2*h), ((points[:,1]-e)+(points[:,0]-e))/h**2

    def write_state(self, trial, q, u, e, internal, meta):
        # Save each accepted physical shape, plus the latest material state
        # even when all residuals have yet to pass.
        path = self.args.folder/f"state_{self.iteration:05d}.npz"
        if path.exists():
            raise FileExistsError(path)
        save_fields(path, x=self.x, solid=trial.b.unpack(q), q=q, field=u,
                    energy=e, shape_params=trial.f.shape.params, shape_extra=trial.f.shape.extra,
                    state_sha=digest(self.args.state), radius=4, winding=1)
        volumes = orientation(trial.f.shape, trial.f.polar_shift, trial.b, trial.b.unpack(q))
        self.last = dict(iteration=self.iteration, file=path.name, sha256=digest(path), orientation=volumes,
                         x=self.x.tolist(), energy=e, material=internal, field=meta,
                         **geometry(trial.f.shape, trial.f.polar_shift))
        atomic_json(self.args.folder/"latest.json", self.last)
        log("accepted", **self.last)
        if not volumes["all_open_cells_oriented"]:
            raise ArithmeticError(f"Candidate fails open-cell orientation check: {volumes}")

    def shape_step(self, grad, diagonal, q, u, energy):
        base_metric = np.diag(1./np.maximum(abs(diagonal), 1.))
        if self.shape_metric is None:
            self.shape_metric = base_metric
        for metric in (self.shape_metric, base_metric):
            direction = -(metric@grad)
            direction *= min(1., .035/max(1.e-20, abs(direction).max()))
            if direction@grad >= 0:
                continue
            for half in range(26):
                dx = direction * 2.**-half
                try:
                    proposed = self.build(self.x+dx)
                    if not orientation(proposed.f.shape, proposed.f.polar_shift, proposed.b,
                                       proposed.b.unpack(q))["all_open_cells_oriented"]:
                        continue
                    fixed = proposed.fixed(q,u)
                    del proposed
                except ValueError as exc:
                    if not infeasible(exc):
                        raise
                    continue
                if fixed < energy + 1.e-4*(grad@dx):
                    return dx, metric
        return None, None

    def run(self):
        q,u,t = self.q, self.u, self.trial
        self.trial = None
        last_x = last_g = None
        u, current_energy, current_gradient, load, info = t.match(q, u)
        initial = t.b.residual(q, current_gradient)
        edge = solid_readout(self.args.state, t.f.shape, t.f.polar_shift, t.b, t.b.unpack(q))
        atomic_json(self.args.folder/"initial_material.json", edge)
        self.write_state(t, q, u, current_energy, initial, info)
        self.iteration += 1
        if not edge["orientation"]["all_open_cells_oriented"]:
            raise ArithmeticError("Initial placement fails orientation check")
        while True:
            # Gauss-Seidel update of the three physical blocks with every
            # internal coordinate free on each shape. Do not reject a shape
            # just because its fixed-shape inner iteration is still moving.
            e0 = current_energy
            q, steps, attempt = relax_material(t, q, u, self.history, 36)
            u, e, gq, load, info = t.match(q, u)
            if e>e0+1.e-9:
                raise ArithmeticError(f"Material and coherent update raised the energy by {e-e0}")
            internal = t.b.residual(q, gq)
            self.write_state(t, q, u, e, internal, info)
            g, diag = self.work(self.x, q, u, e, STEP)
            surface = t.f.measure(u,load)
            values = dict(iteration=self.iteration, shape_force=(-g).tolist(),
                          fixed_material_shape_curvatures=diag.tolist(),
                          solid=internal, field=info, material_steps=steps,
                          material_stop=attempt, current=surface)
            self.last.update(values)
            atomic_json(self.args.folder/"latest.json", self.last)
            log("residuals", **values)
            if last_x is not None:
                dx, dg = self.x-last_x, g-last_g
                prod = dx@dg
                if prod>1.e-10*np.linalg.norm(dx)*np.linalg.norm(dg):
                    b = np.eye(len(dx))-np.outer(dx,dg)/prod
                    self.shape_metric = b@self.shape_metric@b.T+np.outer(dx,dx)/prod
            if abs(g).max()<SHAPE_TOL and material_norm(internal)<SOLID_TOL:
                half, half_curv = self.work(self.x, q,u,e,STEP/2)
                self.last.update(half_shape_force=(-half).tolist(),
                                 half_fixed_material_shape_curvatures=half_curv.tolist(),
                                 half_force_disagreement=float(abs(half-g).max()))
                if abs(half-g).max()>1.e-4 or abs(half).max()>5.e-4:
                    raise ArithmeticError(f"Signed shape gradients disagree: {(-g).tolist()} / {(-half).tolist()}")
                self.finish(t,q,u,load,surface)
                return
            dx, metric = (self.shape_step(g, diag, q,u,e) if abs(g).max()>=SHAPE_TOL
                          else (None,None))
            self.iteration += 1
            if dx is None:
                last_x = last_g = None
                if not steps:
                    raise ArithmeticError(f"No resolved energy descent, shape {(-g).tolist()}, material {internal}")
                current_energy = e
                continue
            last_x, last_g = self.x.copy(), g.copy()
            self.x = self.x+dx
            self.shape_metric = metric
            del t
            t = self.build(self.x, True)
            u, current_energy, _, _, _ = t.match(q, u)
            if current_energy > e+1.e-9:
                raise ArithmeticError(f"Rematching after a shape step raised energy by {current_energy-e}")

    def finish(self,t,q,u,load,measured):
        rest, plain, controls, zero = t.f.controls(u,load,rest_initial=self.rest,return_model=True)
        audit = audit_endpoint(self.base.c,t.f.model,u,rest,zero,plain,
                               state_path=self.args.state,body=t.b,solid=t.b.unpack(q))
        save_fields(self.args.folder/"endpoint_controls.npz",rest=rest,unwound=plain)
        report = dict(**self.last, controls=controls, audit=audit)
        atomic_json(self.args.folder/"endpoint.json", report)
        order_bounds = [controls[k]["nodal_order_range_in_band"] for k in ("wound","rest","unwound")]
        checks = (measured["feed_current"]>0, (measured["dense_total_subsonic_share"] or 0.)>0,
                  measured["port_error"]<3.e-5, abs(measured["disk_minus_feed_current"])<3.e-5,
                  all(lo>=-2.e-12 and hi<=1.+2.e-12 for lo,hi in order_bounds),
                  audit["solid"]["all"]["minimum_J"]>0, audit["solid"]["orientation"]["all_open_cells_oriented"],
                  all(s["placement_error"]<1.e-9 for s in audit["solid"]["physical_neighborhood_sites"]),
                  all(c["stationarity"]<FIELD_TOL and abs(c["imbalance"])<3.e-5 for c in audit["controls"].values()),
                  controls["unwound"]["feed_current"]>0, audit["enclosure"]["resolved"])
        if not all(checks):
            raise ArithmeticError(f"Current, volume or enclosing surface unresolved: {checks}")
        # Comparing the shell forces and physical disk flux on both meshes is
        # part of the independent release decision and requires both outputs.
        log("mesh_endpoint", iteration=self.iteration, endpoint=digest(self.args.folder/"endpoint.json"),
            qualification="Await comparison to the other material mesh and independent audit")

    def failure(self,error):
        record=dict(error=str(error), error_type=type(error).__name__, last=self.last)
        atomic_json(self.args.folder/"issue.json", record)
        if self.last:
            try:
                with np.load(self.args.folder/self.last["file"]) as f:
                    trial=self.build(f["x"])
                    record["saved_material"]=solid_readout(self.args.state,trial.f.shape,trial.f.polar_shift,
                                                          trial.b,f["solid"])
            except (ArithmeticError,ValueError,RuntimeError,OSError) as exc:
                record["saved_material_readout_error"]=repr(exc)
            atomic_json(self.args.folder/"issue.json", record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("state","wall","field","proof","folder"):
        parser.add_argument(key,type=Path)
    for key in ("n","k","angles","degree","subdivisions"):
        parser.add_argument("--"+key,type=int,required=True)
    args=parser.parse_args()
    run=Release(args)
    try:
        run.run()
    except (ArithmeticError,ValueError,RuntimeError) as exc:
        run.failure(exc)
        raise
