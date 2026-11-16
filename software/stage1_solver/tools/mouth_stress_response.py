"""Linear carrier load from the mouth-side normal slopes, and local matched controls."""
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import CubicSpline
import sympy as sp

sys.path.insert(0, str(Path("software/stage1_solver/tools").resolve()))
from held_crossflow_budget import archived_packet


def angular_projection():
    r, a = sp.symbols("r a", positive=True)
    mu, U = sp.symbols("mu U", real=True)
    n0, n1, j0, j1 = [sp.Function(k)(r) for k in ("n0", "n1", "J0", "J1")]
    n = n0 + U*n1*mu
    jump = j0 + U*j1*mu  # jump in normal log-density slope, plus minus minus
    dx_jump = mu*sp.diff(jump, r) + (1-mu**2)/r*sp.diff(jump, mu)
    direct = 2*sp.pi*r**2*sp.integrate(sp.diff(n*dx_jump, U).subs(U, 0), (mu, -1, 1))
    weight = 4*sp.pi/3
    bulk = weight*r**2*(n1*sp.diff(j0, r)-j1*sp.diff(n0, r))
    boundary = weight*r**2*n0*j1
    print("angular radial integrand before multiplying by hbar^2/(4m):", sp.factor(direct))
    print("after integration by parts:", bulk, "plus boundary", boundary)
    if sp.simplify(direct-bulk-sp.diff(boundary, r)) != 0:
        raise ArithmeticError("Angular differentiation and radial reduction disagree")
    # A noncompact jump equal to x on a constant density has a pure surface term.
    uniform = sp.simplify(direct.subs({n0:1, n1:0, j0:0, j1:r}).doit())
    sphere = sp.integrate(uniform, (r,0,a))
    print("edge control, integral of d_x x on ball =", sphere)
    assert sp.simplify(sphere-4*sp.pi*a**3/3) == 0

    # Differentiate the translated compact surface energy in Cartesian coordinates.
    x, X, t = sp.symbols("x X t", real=True)  # t = squared transverse radius
    b, s, kappa, c = sp.symbols("b s kappa c", real=True, nonzero=True)
    radius2 = ((x-X)**2+t)/a**2
    ref = 1+b*radius2
    face = 1+b*(x**2+t)/a**2+U*s*x/a
    energy = kappa*(1-radius2)**3*(face-ref)**2/2  # inside the ball; cubic cutoff kills edge work
    force1 = -sp.diff(energy, X, U).subs({X:0,U:0})
    integrand = force1.subs({x:r*mu, t:r**2*(1-mu**2)})
    energy_force = sp.integrate(2*sp.pi*r**2*sp.integrate(integrand, (mu,-1,1)), (r,0,a))
    # c = hbar^2/(2m); real surface energy gives jump = 2 e_n/c.
    jump1 = 2*kappa/c * (1-(r/a)**2)**3 * s*r/a
    stress_integrand = direct.subs({n0:1+b*(r/a)**2, n1:s*r/a,j0:0,j1:jump1}).doit()
    stress_force = sp.factor(c*sp.integrate(stress_integrand,(r,0,a))/2)
    print("change of weak force, translated energy =", sp.factor(energy_force))
    print("change of weak force, two-sided stress =", stress_force)
    assert stress_force != 0 and sp.simplify(energy_force-stress_force) == 0


def trace(speed):
    name = "mouth_feed_weak_run.json" if speed == 0 else "mouth_intake_redistribution.json"
    path = archived_packet(name, f"n32_restoring_slope_u{speed:g}.npz")
    with np.load(path, allow_pickle=False) as packet:
        field = packet["field"]
        density = abs(field[-1])**2
        # Register the face independently by its physical depth in the feed solver.
        depths = np.arange(len(field)) * float(packet["dz"])
        mouth_plane, = np.flatnonzero(depths == 1.0)
        if depths[-1] != 1.0:
            raise ValueError("The feeding slab must end at the mouth")
        registered = field[mouth_plane]
        np.testing.assert_allclose(density, registered.real**2 + registered.imag**2,
                                   rtol=2e-15, atol=0, err_msg="Density is not from the z=1 mouth plane")
        mouth, current = packet["mouth"], packet["current"]
        dx, actual = float(packet["dx"]), float(packet["speed"])
        assert actual == speed and density.shape == current.shape == mouth.shape == (32,32,32)
        assert np.all(current[mouth]>0) and np.all(current[~mouth]==0) and np.all(density>0)
        np.testing.assert_allclose(dx**3*current.sum(), .25, rtol=0, atol=1e-12)
    return density, mouth, dx


def archived_traces():
    ref, mouth, dx = trace(0)
    x = (np.arange(ref.shape[0]) + .5 - ref.shape[0]/2)*dx
    axial = x[:,None,None]
    transverse = x[None,:,None]**2+x[None,None,:]**2
    def envelope(offset):
        return np.maximum(1-(axial-offset)**2-transverse,0)
    env = envelope(0)
    h = env**3
    hx = -6*axial*env**2
    np.testing.assert_array_equal(mouth, h>0)
    np.testing.assert_allclose(ref, ref[::-1], atol=1e-13, rtol=0)
    # Register a mouth-attached reference density through the sampled rest profile.
    # The lateral box is periodic; the reference spline translates with the opening.
    spline = CubicSpline(np.append(x, x[-1]+dx), np.concatenate([ref,ref[:1]]),
                         axis=0, bc_type="periodic")
    ref_x = spline(x, 1)
    def force(density):
        delta = density-ref
        return dx**3*np.sum(.5*hx*delta**2 - h*delta*ref_x)
    def shifted_energy(density, offset):
        return .5*dx**3*np.sum(envelope(offset)**3*(density-spline(x-offset))**2)
    def from_energy(density, step):
        return (shifted_energy(density, -step)-shifted_energy(density, step))/(2*step)
    print("added rest energy, chemical-potential jump, force:",
          shifted_energy(ref, 0), np.max(abs(h*(ref-ref))), force(ref))
    for speed in (.01, .02):
        pos, aperture, spacing = trace(speed)
        np.testing.assert_array_equal(aperture, mouth)
        assert dx == spacing
        neg = pos[::-1]
        odd = (force(pos)-force(neg))/(2*speed)
        dipole = -dx**3*np.sum(h*(pos-neg)/(2*speed)*ref_x)
        estimates = [(from_energy(pos, step)-from_energy(neg, step))/(2*speed)
                     for step in (dx/1000, dx/2000)]
        print(f"|U|={speed:g} per_unit_kappa odd_force/U={odd:.12g} "
              f"rest_gradient_contraction={dipole:.12g} translated_energy={estimates[0]:.12g},{estimates[1]:.12g}")
        np.testing.assert_allclose(force(pos), -force(neg), rtol=1e-10, atol=1e-15)
        np.testing.assert_allclose(estimates, odd, rtol=3e-6, atol=1e-14)
        assert abs(odd) > 1e-12


if __name__ == "__main__":
    angular_projection()
    archived_traces()
