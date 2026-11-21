"""A renewing material contact, tested with held twists and translational drives.

Local units: one mean feed residence, one lever length, one constituent mass.
The participating arrival rate j is fixed; store and force outputs are per j.
The global angular test uses only the sealed finite-shell energy certificates.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp


ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_collar_charge_transmission.json")
RATE = 1/256  # Internal angular registration rate in inverse feed-dwell units.
STRENGTH = 2.  # j k ell² / K, in inverse feed-dwell units.


def calibrated_contacts():
    # Test the lever model directly. A balanced pair of springs, each k/2,
    # at opposite unit offsets is free to swivel before registration. After
    # registration the reference direction is supplied by the prepared patch.
    def energy(x, angle):
        offset = np.array([math.cos(angle)-1, math.sin(angle)])
        center = np.array([x, 0.])
        return sum(np.dot(center+sign*offset, center+sign*offset)/4 for sign in (-1,1))
    for x, angle in ((.02,.04),(-.02,-.04),(.01,.02),(.13,math.pi/2)):
        assert math.isclose(energy(x,angle), x*x/2+1-math.cos(angle), abs_tol=3e-16)
    # Recover the turn strength per survivor from a signed cohort probe. Fit
    # at one age and reserve the second age and half-amplitude for prediction.
    def cohort_torque(age, angle):
        registered = quad(lambda t: RATE*math.exp(-RATE*t), 0, age)[0]
        return registered * math.sin(angle)
    angle = .04
    p_fit = (cohort_torque(128,angle)-cohort_torque(128,-angle))/(2*math.sin(angle))
    alpha = -math.log1p(-p_fit)/128
    p_check = cohort_torque(256,angle/2)/math.sin(angle/2)
    print(f"held_cohort surviving_registered_age128={p_fit:.12g} fit_rate={alpha:.12g} "
          f"surviving_registered_age256={p_check:.12g} predicted_age256={-math.expm1(-alpha*256):.12g}")
    assert abs(p_check+math.expm1(-alpha*256))<1e-13
    contacts = {}
    for law in ("delayed","immediate","free_swivel"):
        values = []
        for name, tau in (("feed",1.),("rim",128.)):
            def age_integrand(a):
                # The integration variable is age in local residence units.
                alignment = {"delayed":-math.expm1(-alpha*tau*a),
                             "immediate":1., "free_swivel":0.}[law]
                return math.exp(-a)*alignment
            aligned = quad(age_integrand,0,math.inf,epsabs=1e-13)[0]
            modulus = STRENGTH*tau*aligned
            print(f"twist law={law} patch={name} residence={tau:g} aligned_fraction={aligned:.12g} "
                  f"modulus_over_K={modulus:.12g} torque_at_point04_over_K={modulus*math.sin(angle):.12g}")
            values.append(modulus)
        contacts[law]=values
    return contacts


def classify(contacts):
    rows = json.loads(ARCHIVE.read_text())["runs"]["mouth_collar_attachment_window01"]["stdout_lines"]
    def value(row, key):
        return float(row.split(key+"=")[1].split()[0])
    weak = min(value(row,"transmission_for_g_below") for row in rows)
    strong = max(value(row,"screening_for_g_above") for row in rows)
    assert len(rows)==8 and {int(value(row,"eta")) for row in rows}=={-1,1}
    print(f"archived finite_cases={len(rows)} weak_ceiling={weak:.12g} strong_floor={strong:.12g}")
    # G_i/K = STRENGTH*tau_i²/(registration_time+tau_i).
    print(f"registration_time_in_feed_dwells transmitted_above={STRENGTH/weak-1:.12g} "
          f"rim_certificate_through={STRENGTH*128**2/64-128:.12g} "
          f"screened_below={STRENGTH/strong-1:.12g}")
    for name,(feed,rim) in contacts.items():
        if feed==rim==0:
            # Free field minimum has zero spring and vacuum energy at zero flux.
            result="0,0"
        elif rim>=64 and feed<weak:
            result="+1,-1"
        elif rim>=64 and feed>strong:
            result="0,0"
        else:
            result="unclassified"
        print(f"angular_minimum law={name} mirror_pair={result}")


def mechanics(tau):
    # m=k=ell=1, Ω=1. Relative spring extension is zero at birth, and a
    # captured parcel enters with the incident sideways velocity U.
    decay = 1/tau
    # For constant U,V, a parcel of age a has v=U+(V-U)(1-cos(a)),
    # e=(V-U)sin(a). Average over residence without oscillatory quadrature
    # over a long lifetime: scipy's cosine/sine weighted Laplace integrals.
    def trigonometric_mean(mode, frequency):
        return quad(lambda age: decay*math.exp(-decay*age),0,math.inf,
                    weight=mode,wvar=frequency,epsabs=1e-13)[0]
    cosine = trigonometric_mean("cos",1.)
    sine = trigonometric_mean("sin",1.)
    cosine2 = trigonometric_mean("cos",2.)
    drag = tau*sine                  # F_external/(j m) at unit towing speed
    outgoing = 1-cosine              # net departure momentum / (j m)
    heat = (1-cosine2)/4            # spring energy disposed at release / j
    kinetic = (1-2*cosine+(1+cosine2)/2)/2
    print(f"steady residence={tau:g} mass_in_stock_over_mass_rate={tau:g} "
          f"tow_force_over_mass_rate_V={drag:.12g} held_pull_over_mass_rate_U={drag:.12g} "
          f"exit_px_over_mass_rate_V={outgoing:.12g} "
          f"exit_kinetic_over_mass_rate_V2={kinetic:.12g} release_heat_over_mass_rate_V2={heat:.12g}")
    assert np.allclose([outgoing,kinetic+heat],[drag,drag],atol=2e-12,rtol=0)

    # Force-controlled, zero geometric mass. A prescribed smooth push F/(jm)
    # starts from rest. Its measured drive force equals total spring reaction,
    # so <e>=F/(Nk); solve for mean material momentum, then for port velocity.
    end = 8*tau
    ramp = min(tau,4.)
    def ramp_shape(t):
        z = math.exp(-t/ramp)
        return (1-z)**3, 3*z*(1-z)**2/ramp
    speeds = []
    for amplitude in (.01,-.01,.005,-.005):
        def rhs(t,p):
            return amplitude*ramp_shape(t)[0]-decay*p
        times = np.array([0.,tau/8,tau,end])
        sol = solve_ivp(rhs,(0,end),[0.],t_eval=times,rtol=1e-11,atol=1e-14)
        assert sol.success
        p = sol.y[0]   # stored momentum divided by the mass throughput
        v = []
        for t,measured in zip(times,p):
            force, force_slope = (amplitude*x for x in ramp_shape(t))
            v.append((force_slope+decay*force+measured)/tau)
            history = quad(lambda a: amplitude*ramp_shape(t-a)[0]*math.exp(-decay*a),0,t)[0]
            assert abs(history-measured)<2e-10
        speeds.append(np.array(v))
        print(f"push residence={tau:g} F_plateau_over_mass_rate={amplitude:g} "
              f"speed_at_tau_over_8={v[1]:.12g} speed_at_tau={v[2]:.12g} speed_at_8tau={v[3]:.12g} "
              f"terminal_speed={amplitude/drag:.12g}")
    assert np.allclose(speeds[0],-speeds[1]) and np.allclose(speeds[0],2*speeds[2])
    assert np.allclose(speeds[1],2*speeds[3])


if __name__ == "__main__":
    classify(calibrated_contacts())
    mechanics(1.)
    mechanics(128.)
