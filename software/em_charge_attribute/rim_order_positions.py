"""Independent geometric/order markers at exposure onset and a contact control."""

from pathlib import Path
from zipfile import ZipFile

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre, spherical_jn


def emit(label, value):
    print(f"{label} = {value}", flush=True)


# Verify which traces the archived candidates actually contain.
for archive in (
    "_scratch/s11c/lateral_bath_e24_c4_01/fields/d0.npz",
    "_scratch/s11c/signed_mouth_collar01/fields/benchmark_archive_128_positive.npz",
):
    with ZipFile(archive) as saved:
        emit(str(Path(archive)), " ".join(sorted(saved.namelist())))

nu = sp.symbols("nu", real=True)
z = sp.symbols("z", positive=True)
# Half-difference of sin^2((z/2)*nu +/- pi/4) in the strain energy.
odd_strain = sp.trigsimp((sp.sin(z*nu/2 + sp.pi/4)**2 -
                          sp.sin(z*nu/2 - sp.pi/4)**2)/2, method="fu")
projected = sp.factor(sp.Rational(3,2)*sp.integrate(nu*odd_strain, (nu,-1,1)))
emit("dipole contour speed divided by M mu A^2 k^2/(4p)", projected)

# Smooth counterpart of the tested weak jump.  The contour is halfway
# between the two plateau values.  All lengths and slow-time inputs here
# are dimensionless.  The fields below are first Taylor increments only.
eps, width, k, amplitude = 0.02, 0.005, 0.5, 0.1
slope = eps/(4*width)
directions, weights = roots_legendre(28)
def centers(radius, step, phase, amp, boost):
    drift = boost*step
    rim_edges, order_edges = [], []
    for direction in directions:
        def signed_distance(r):
            return np.sqrt(r*r+drift*drift-2*r*drift*direction) - radius

        def order_level(r):
            initial = eps*(expit(signed_distance(r)/width) - 0.5)
            energy = amp**2*k**2*np.sin(k*(r*direction-drift)+phase)**2/4
            return initial-step*energy  # mobility and stiffness are one

        bracket = (radius-0.1, radius+0.1)
        rim_edges.append(brentq(signed_distance,*bracket,xtol=1e-14))
        order_edges.append(brentq(order_level,*bracket,xtol=1e-14))

    def centroid(edges):
        r = np.asarray(edges)
        return 0.75*np.dot(weights, directions*r**4)/np.dot(weights, r**3)
    return np.array([centroid(rim_edges), centroid(order_edges)])


def odd_speed(radius, step, amp, boost):
    baseline = centers(radius,step,0,0,boost)
    plus = centers(radius,step,np.pi/4,amp,boost) - baseline
    minus = centers(radius,step,-np.pi/4,amp,boost) - baseline
    return (plus-minus)/(2*step), np.max(np.abs(plus+minus))/step


def exact_front_speed(radius):
    return 3*amplitude**2*k**2*spherical_jn(1,2*k*radius)/(8*slope)


errors = []
boost_residuals = []
for step in (1/4,1/16,1/64):
    measured = []
    for radius in (1,1.5):
        speed, mirror_error = odd_speed(radius,step,amplitude,0.03)
        at_rest, _ = odd_speed(radius,step,amplitude,0.)
        boost_residuals.extend(abs(speed-at_rest))
        measured.extend(speed)
        assert mirror_error < 1e-9
    error = max(abs(measured[j] - exact_front_speed(r))
                for j,r in ((1,1),(3,1.5)))
    errors.append(error)
    emit(f"tangent step {step:g}: mid rim/front and end wall/front speeds", measured)
emit("max errors from analytic order speeds by step", errors)
emit("largest boost difference in velocities", max(boost_residuals))
ratio = exact_front_speed(1.5)/exact_front_speed(1)
emit("order speed at end divided by midplane", ratio)
half_read, _ = odd_speed(1,1/64,amplitude/2,0.)
full_read, _ = odd_speed(1,1/64,amplitude,0.)
emit("half/full amplitude contour speed ratio at finest step", half_read[1]/full_read[1])

def flare(s):
    return 1+(3*s*s-2*s**3)/2

bow_speed = 2*quad(lambda s: exact_front_speed(flare(s))*np.cos(np.pi*s),0,1)[0]
emit("order cosine-mode onset speed for reflected positive-minus-negative drive", bow_speed)
assert errors[-1] < errors[0] and errors[-1] < 1e-8
assert max(boost_residuals) < 1e-9

# Compare two permitted illustrative attachments. Coordinates are in units
# V*t0 with V the preceding signed cosine-mode speed for positive phase.
# Time is in units t0, with M=gamma=Kwall=Korder=1 in these units.
def trajectory(contact, forcing):
    def rhs(t,state):
        x, velocity, y = state
        transmitted = contact*(y-x)
        return [velocity, transmitted-x, forcing-y-transmitted]
    return solve_ivp(rhs,(0,1),[0.,0.,0.],t_eval=[1.],rtol=1e-11,atol=1e-13,
                     method="DOP853").y[[0,2],-1]


for contact in (0,1):
    quiet = trajectory(contact,0)
    plus = trajectory(contact,1) - quiet
    minus = trajectory(contact,-1) - quiet
    emit(f"contact {contact}: normalized wall, order at t=1", list((plus-minus)/2))
    assert max(abs(plus+minus)) < 1e-12
