"""Linear sideways translation of a wall and a standing shear envelope on a short segment.

Both displacements and both momenta are independent. Only the geometric rim
position of the wall is driven. A calibrated terminal load replaces the
unresolved continuation. Units are those of the supplied draining state.
"""

import numpy as np
from scipy.integrate import solve_ivp


def fit_terminal_load(train_velocity, train_force, check_velocity, check_force, relative_tolerance):
    """Fit F_out = Z v at the terminal section and test independent readings."""
    def samples(v, f):
        v, f = np.asarray(v, dtype=float), np.asarray(f, dtype=float)
        if v.ndim != 2 or v.shape[1] != 2 or v.shape != f.shape:
            raise ValueError("Terminal readings need two velocities and two forces per sample")
        if not np.all(np.isfinite(v)) or not np.all(np.isfinite(f)) or np.linalg.matrix_rank(v) < 2:
            raise ValueError("Terminal readings must resolve both independent motions")
        return v, f

    if not 0 < relative_tolerance < 1:
        raise ValueError("Provide a relative calibration tolerance between zero and one")
    v, f = samples(train_velocity, train_force)
    vc, fc = samples(check_velocity, check_force)
    # Fit a reciprocal load directly, using both measured traction components.
    design = np.zeros((len(v), 2, 3))
    design[:, 0, :2] = v
    design[:, 1, 1:] = v
    a, b, d = np.linalg.lstsq(design.reshape(-1, 3), f.ravel(), rcond=None)[0]
    z = np.array([[a, b], [b, d]])
    errors = [np.linalg.norm(speed @ z.T-force) /
              max(np.linalg.norm(force), np.finfo(float).tiny)
              for speed, force in ((v,f), (vc,fc))]
    if np.linalg.eigvalsh(z).min() < 0 or max(errors) > relative_tolerance:
        raise ValueError(f"Passive instantaneous load fails calibration: relative errors {errors}")
    return z, errors


def solve_segment(z, cut, line_mass, tension, coupling, termination, acceleration, times,
                  distributed_force=None):
    """From incremental rest, impose X=A t²/2. Input profiles include all mode/area weights.

    z contains the rim, interior measuring cut and downstream end. The wave is
    traction-free at the rim. distributed_force(t) optionally returns a force
    per unit intrinsic length on each field and node, in the lab frame.
    """
    z, times = np.asarray(z, dtype=float), np.asarray(times, dtype=float)
    line_mass, tension, coupling, termination = (
        np.asarray(p, dtype=float) for p in (line_mass, tension, coupling, termination))
    if z.ndim != 1 or len(z) < 3 or z[0] != 0 or np.any(np.diff(z) <= 0):
        raise ValueError("The segment must start at zero and contain increasing depths")
    if not 0 < cut < len(z)-1 or not isinstance(cut, (int, np.integer)):
        raise ValueError("The cut must be an interior node")
    if (times.ndim != 1 or not len(times) or times[0] <= 0 or np.any(np.diff(times) <= 0)
            or any(not np.isfinite(p).all() for p in
                   (z, times, line_mass, tension, coupling, termination, acceleration))):
        raise ValueError("Invalid segment values or sample times")
    if (line_mass.shape != (2,len(z)) or tension.shape != line_mass.shape
            or coupling.shape != z.shape or np.any(line_mass < 0)
            or np.any(tension < 0) or np.any(coupling < 0)):
        raise ValueError("Expected two nonnegative mass/tension profiles and one coupling profile")
    if (termination.shape != (2,2) or not np.allclose(termination,termination.T,atol=0,rtol=1e-12)
            or np.linalg.eigvalsh(termination).min() < 0):
        raise ValueError("Terminal impedance must be reciprocal and passive")

    # Integrate piecewise-linear mass profiles over the nodal control volumes.
    dz = np.diff(z)
    left_mass = (3*line_mass[:,:-1]+line_mass[:,1:])*dz/8
    right_mass = (line_mass[:,:-1]+3*line_mass[:,1:])*dz/8
    mass = np.zeros_like(line_mass)
    mass[:,:-1] += left_mass
    mass[:,1:] += right_mass
    if np.any(mass <= 0):
        raise ValueError("There must be a resolved mode with positive weight in every nodal volume")
    width = np.r_[dz[0], dz[:-1]+dz[1:], dz[-1]]/2
    stiffness = (tension[:,1:]+tension[:,:-1])/(2*dz)
    # Treat each lumped momentum as uniform within its cell. Use one geometric
    # cut for the two fields, so internal wall-wave forces cancel there too.
    share = dz[cut-1]/(dz[cut-1]+dz[cut])
    inside = np.zeros_like(mass)
    inside[:,:cut], inside[:,cut] = 1, share
    size = 2*len(z)

    def forces(t, position, velocity):
        flux = -stiffness*np.diff(position,axis=1)
        end_load = termination @ velocity[:,-1]
        spring = coupling*(position[1]-position[0])
        sources = np.stack([spring,-spring])*width
        if distributed_force is not None:
            applied = np.asarray(distributed_force(t),dtype=float)
            if applied.shape != sources.shape or not np.isfinite(applied).all():
                raise ValueError("Distributed force must provide both fields on this grid")
            sources += applied*width
        net = sources.copy()
        net[:,:-1] -= flux
        net[:,1:] += flux
        net[:,-1] -= end_load
        rim_force = mass[0,0]*acceleration-net[0,0]
        net[0,0] = mass[0,0]*acceleration
        export = (1-share)*flux[:,cut-1]+share*flux[:,cut]
        return net, export, rim_force, np.sum(inside*sources,axis=1)

    def rhs(t, state):
        pos, vel = state[:2*size].reshape(2,2,-1)
        net, export, rim, source = forces(t,pos,vel)
        return np.r_[vel.ravel(), (net/mass).ravel(), export, rim, source]

    initial = np.zeros(2*size+5)
    sol = solve_ivp(rhs,(0,times[-1]),initial,t_eval=times,method="DOP853",
                    rtol=2e-10,atol=5e-17)
    if not sol.success:
        raise ArithmeticError(sol.message)
    pos, vel = sol.y[:2*size].reshape(2,2,len(z),len(times))
    stored = np.einsum("ij,ijt->it",mass*inside,vel)
    cut_flux = np.column_stack([forces(t,pos[:,:,i],vel[:,:,i])[1] for i,t in enumerate(times)])
    impulse = sol.y[2*size:2*size+2]
    rim_impulse = sol.y[2*size+2]
    source_impulse = sol.y[2*size+3:]
    residual = stored + impulse - source_impulse
    residual[0] -= rim_impulse
    return dict(z=z, times=times, cut_depth=z[cut], displacement=pos, velocity=vel,
                cut_flux=cut_flux, stored=stored, cut_impulse=impulse, rim_impulse=rim_impulse,
                source_impulse=source_impulse, balance_residual=residual)
