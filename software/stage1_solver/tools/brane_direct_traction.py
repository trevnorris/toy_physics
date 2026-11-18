"""Direct net traction on the homogeneous rotational brane at a spherical rim."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp


def check_admittance():
    # Derive both local traces from the outgoing azimuthal vector potential.
    # At the equator, positive tangent is along the drive, opposite e_theta.
    a, r, s, c, rho = sp.symbols("a r s c rho", positive=True)
    k = s/c
    f = sp.exp(-k*(r-a)) * (1+k*r)/r**2
    p = 2*f/r
    d = sp.diff(r*f, r)/r
    velocity = s*d
    traction = rho*c**2 * (p-sp.diff(r*d, r))/r
    admittance = sp.factor((velocity/traction).subs(r, a))
    print("SymPy equatorial velocity per net traction =", admittance, flush=True)
    tau = a/c
    expected = (1 + 1/(s*tau*(1+s*tau)))/(rho*c)
    if sp.simplify(admittance-expected) != 0:
        raise ArithmeticError("Local traction/velocity normalization failed")


def velocity_history(radius, density, speed, duration, peak, cells):
    """Radial wave solve for h=r*(curl u)_phi/sin(theta), from rest."""
    modulus = density*speed**2
    # No wave reaches the far boundary during the sample window.
    dr = 1.5*speed*duration/cells
    r = radius + np.arange(1, cells)*dr
    angular_term = 2/r**2
    steps = 4 * int(np.ceil(speed*duration/(4*0.5*dr)))
    dt = duration/steps
    count = len(r)
    # Interior h, interior time derivative of h, and equatorial brane velocity.
    state = np.zeros(2*count+1)

    def rhs(t, values):
        h = np.empty(cells+1)
        h[0] = radius * peak * (t/duration)**3 / modulus
        h[1:-1] = values[:count]
        h[-1] = 0.
        out = np.empty_like(values)
        out[:count] = values[count:-1]
        laplacian = (h[2:] - 2*h[1:-1] + h[:-2])/dr**2
        out[count:-1] = speed**2 * (laplacian - angular_term*h[1:-1])
        # Acceleration from the material momentum equation at a+, using only
        # the one-sided spatial derivative of h. No imposed rim velocity.
        gradient = (-3*h[0] + 4*h[1] - h[2])/(2*dr)
        out[-1] = -speed**2 * gradient/radius
        return out

    rows = []
    for i in range(steps):
        t = i*dt
        k1 = rhs(t, state)
        k2 = rhs(t+dt/2, state+dt*k1/2)
        k3 = rhs(t+dt/2, state+dt*k2/2)
        k4 = rhs(t+dt, state+dt*k3)
        state += dt*(k1+2*k2+2*k3+k4)/6
        if i+1 in (steps//4, steps//2, steps):
            t = (i+1)*dt
            rows.append(dict(t=t, traction=peak*(t/duration)**3, velocity=state[-1]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--density", type=float, required=True)
    parser.add_argument("--wave-speed", type=float, required=True)
    parser.add_argument("--duration", type=float, required=True)
    parser.add_argument("--peak-traction", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    params = (args.radius, args.density, args.wave_speed, args.duration)
    if not all(np.isfinite(p) and p > 0 for p in (*params, args.peak_traction)):
        parser.error("All material, scale and load parameters must be positive and finite")
    check_admittance()
    grids = (384, 768, 1536)
    runs = []
    for factor in (0, 1, -1, 0.5, -0.5):
        peak = factor*args.peak_traction
        histories = {str(n): velocity_history(*params, peak, n) for n in grids}
        runs.append(dict(peak=peak, grids=histories))
    print("Homogeneous direct traction; peak, then rim velocities at T/4, T/2, T (finest grid)")
    for run in runs:
        print(f"{run['peak']:+.12g}", *(f"{row['velocity']:+.12g}" for row in run["grids"]["1536"]))
    packet = dict(radius=args.radius, density=args.density, speed=args.wave_speed,
                  duration=args.duration, runs=runs)
    with args.output.open("x") as file:
        json.dump(packet, file, indent=2, allow_nan=False)
        file.write("\n")
    for n in grids:
        traces = np.array([[row["velocity"] for row in run["grids"][str(n)]] for run in runs])
        np.testing.assert_array_equal(traces[0], 0.)
        np.testing.assert_allclose(traces[1:], np.array([1, -1, .5, -.5])[:, None]*traces[1],
                                   rtol=1e-12, atol=0)


if __name__ == "__main__":
    main()
