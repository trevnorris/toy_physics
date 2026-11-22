"""Release a prestressed collar with a separately prepared, permanently bonded core.

Follow the joint initial downhill displacement through its first directional
minimum or loss of a resolved mode. Test independent shape and order work there.
The entrance and the topology of the material chart are specified throughout.
"""
import argparse
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import expit, logit

from mouth_poly_shallow_release import Trial, emit
from mouth_sealed_shape_relax import Meridian


def run(radius, mesh):
    trial = Trial(radius, mesh, bonded_core=True)
    initial_order = trial.chi
    origin = np.zeros(6)
    energy0, grad, order, fields, metrics, _ = trial.evaluate(origin, initial_order)
    emit(kind="prepared", radius=radius, mesh=mesh, inputs=trial.fingerprints, **metrics)
    seeds = {0.: fields.wave.seed}

    def material_gradient(position, chi, seed, axes, step):
        slopes = []
        for axis in axes:
            probe = np.zeros(6)
            probe[axis] = step
            plus = trial.evaluate(position+probe, chi, seed)[0]
            minus = trial.evaluate(position-probe, chi, seed)[0]
            slopes.append((plus-minus)/(2*step))
        return np.array(slopes)

    gradient = material_gradient(origin, initial_order, seeds[0.], range(6), .0005)
    displacement = -gradient/np.max(abs(gradient))
    mobile = order.free & (initial_order > 0.) & (initial_order < 1.)
    diagonal = (2*order.stiff.diagonal()+1.5*order.lump).reshape(grad.shape)
    order_rate = -grad[mobile]/diagonal[mobile]
    weights = initial_order[mobile]*(1-initial_order[mobile])
    order_rate /= np.max(abs(weights*order_rate))
    initial_logits = logit(initial_order[mobile])
    order_slope = grad[mobile] @ (weights*order_rate)
    emit(kind="initial_slopes", axes=["log_depth", "log_exit_radius", "middle",
         "upper_minus_lower", "inner_log_spacing", "outer_log_spacing"],
         slopes=gradient.tolist(), displacement=displacement.tolist(),
         order_slope=order_slope, total_slope=gradient@displacement+order_slope)
    energies = {0.: energy0}

    def order_at(step):
        chi = initial_order.copy()
        chi[mobile] = expit(initial_logits + step*order_rate)
        return chi

    def evaluate(step):
        if step in energies:
            return energies[step]
        shifts = step*displacement
        shape = Meridian(radius,trial.x+shifts[:4])
        if shape.wall()["minimum_radius"] <= 0.:
            raise ValueError("The open positive-radius material chart has ended")
        anchor = min(seeds, key=lambda prior: abs(step-prior))
        result = trial.evaluate(shifts, order_at(step), seeds[anchor])
        energy, _, _, loads, readout, _ = result
        if readout["mode_overlap"] < .94:
            raise ArithmeticError("Wave continuation unresolved")
        if abs(readout["current"]-readout["outflow"]) > 1e-6*abs(readout["current"]):
            raise ArithmeticError("Terminal mass balance unresolved")
        energies[step] = energy
        seeds[step] = loads.wave.seed
        emit(kind="bonded_path", step=step, from_step=anchor, change=energy-energy0, **readout)
        return energy

    step = .001
    previous = 0.
    while True:
        try:
            energy = evaluate(step)
        except (ArithmeticError, ValueError, RuntimeError) as error:
            # Keep a usable physical bracket even when one distant trial fails.
            # A smaller displacement distinguishes a continuation gap from an
            # upward turn in total energy. Never assign energy to the failed state.
            emit(kind="rejected_probe", step=step, last_valid=previous, reason=str(error))
            step = .5*(previous+step)
            if step-previous < .0001:
                emit(kind="unresolved_endpoint", last_valid=previous, next_trial=step)
                return
            continue
        if energy > energies[previous]:
            break
        previous = step
        step *= 2.
    if previous == 0.:
        raise ArithmeticError("No finite descent on the bonded branch")
    # Use the last lower-energy point and one point on each side. Rejected
    # states are never included among the potential interpolation nodes.
    left = max(s for s in seeds if s < previous)
    minimum = minimize_scalar(evaluate, bracket=(left,previous,step), method="brent",
                              options={"xtol": 2e-5})
    if not minimum.success:
        raise ArithmeticError("Directional energy minimum unresolved")
    s = float(minimum.x)
    evaluate(s)
    shifts, chi, seed = s*displacement, order_at(s), seeds[s]
    end, g, _, _, info, _ = trial.evaluate(shifts, chi, seed)
    # Log depth and middle width can move independently at this minimum.
    checks = [material_gradient(shifts,chi,seed,[0,2],h).tolist() for h in (.001,.0005)]
    local_slope = g[mobile] @ (chi[mobile]*(1-chi[mobile])*order_rate)
    emit(kind="directional_minimum", step=s, change=end-energy0, signed_steps=[.001,.0005],
         energy_slopes_log_depth_and_middle=checks, order_slope=local_slope,
         max_signed_step_difference=float(np.max(abs(np.diff(checks,axis=0)))),
         **info)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=int, choices=(1,2), required=True)
    parser.add_argument("--mesh", type=int, choices=(48,72), required=True)
    args = parser.parse_args()
    run(args.radius,args.mesh)
