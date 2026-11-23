"""Polish a nearly stationary bonded state using integrated directional work.

The comparison uses gradients of the exact same energy on feasible material
maps. A non-smooth step which fails refinement supplies no stationary reading.
"""
import numpy as np


def polish(body, q, target=3.e-5):
    q = q.copy()
    energy, gradient = body.energy(q)
    history = []
    for iteration in range(300):
        residual = abs(body.native[body.mask]).max()
        if residual < target:
            return q, dict(iterations=iteration, residual=residual)
        direction = gradient.copy()
        coeffs = []
        for step, change, inverse in reversed(history):
            a = inverse*(step@direction)
            coeffs.append(a)
            direction -= a*change
        if history:
            step, change, _ = history[-1]
            direction *= (step@change)/(change@change)
        for (step, change, inverse), a in zip(history, reversed(coeffs)):
            direction += step*(a-inverse*(change@direction))
        moved = False
        for search in (direction,gradient):
            if gradient@search <= 0:
                continue
            for half in range(35):
                step = -(2.**-half)*search
                try:
                    after, new = body.energy(q+step)
                    _, mid = body.energy(q+step/2)
                except ValueError as exc:
                    if str(exc)!='Nonpositive trial map or outside crystalline order grid':
                        raise
                    continue
                slope = gradient@step
                edge = ((gradient+new)@step)/2
                work = (edge+2*(mid@step))/3
                if (work < 1.e-4*slope and abs(edge-work)<.03*abs(slope)
                        and after <= energy+1.e-12*max(1.,abs(energy))):
                    moved=True
                    break
            if moved:
                break
        if not moved:
            break
        change = new-gradient
        prod = change@step
        if prod > 1.e-11*np.linalg.norm(change)*np.linalg.norm(step):
            history.append((step,change,1/prod))
            history=history[-16:]
        q += step
        energy, gradient = body.energy(q)
    body.energy(q)
    residual = abs(body.native[body.mask]).max()
    if residual < target:
        return q, dict(iterations=iteration+1, residual=residual)
    raise ArithmeticError(f"No verified smooth material equilibrium at residual {residual}")
