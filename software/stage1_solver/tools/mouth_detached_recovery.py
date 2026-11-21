"""Infer detached recovery from a weak-shear return on a flat quarter-order patch.

Predictions integrate exact geometric crossing times across the sites occupied
at exit. A separate clocked force instrument tracks sites and the opposite-face
reaction, without accessing the predicted event factors.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


PRIOR = Path("research/pde_ledger_v3/_measurements/mouth_contact_crossing_force.json")
WIDTH = .25
HOLD = .5
BASE_GAP = 3/16


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def attached_factor(gap, plateau, tau):
    # The initial and final partner is 0, and partner 1 is briefly reached at
    # the turning point. For each exit site both off intervals have length gap.
    nodes, weights = leggauss(32)
    phase = .5*WIDTH*nodes
    speed = (1-WIDTH)/gap
    outside = WIDTH/2-phase
    inside = 1-WIDTH/2-phase
    initial = HOLD + outside/speed
    middle = 2*(1-inside)/speed
    final = HOLD + outside/speed
    pre = (1-plateau)*(-np.expm1(-initial/tau))
    middle_motion = (1-plateau) + plateau*np.exp(-middle/tau)
    post = plateau*(-np.expm1(-final/tau))
    return float(weights @ (pre*middle_motion*post)/2)


def support_reading(gap, off_rate, steps, sites):
    # Model units k=K=zeta_on=1. Varying positive detached drag alone
    # changes off_rate=K/zeta_off but leaves every intact test identical.
    anchors = (np.arange(sites)+.419)/sites
    speed = (1-WIDTH)/gap
    travel = 1/speed
    duration = 2*(HOLD+travel)
    dt = duration/steps
    offsets = np.zeros(sites)
    birth = np.zeros(sites)
    active = abs(anchors - np.rint(anchors)) < WIDTH/2
    partner = np.rint(anchors).astype(int)
    continuously_present = active.copy()
    for j in range(1, steps+1):
        # The applied unit step is a normalization; the tangent equation is linear.
        old = offsets.copy() if j == steps else None
        eq = np.where(active, .5*(1+birth), 0.)
        offsets += (eq-offsets)*np.where(active, -np.expm1(-2*dt), -np.expm1(-off_rate*dt))
        position = anchors + np.interp(j*dt, [0., HOLD, HOLD+travel, HOLD+2*travel, duration],
                                       [0., 0., 1., 0., 0.])
        current_partner = np.floor(position + .5).astype(int)
        engaged = abs(position-current_partner) < WIDTH/2
        new = engaged & (~active | (partner != current_partner))
        birth[new] = offsets[new]-1
        continuously_present &= engaged & (partner==current_partner)
        active, partner = engaged, current_partner
    # There is a final stationary hold. Both faces read the same total force
    # including drag on the still-detached features, whose elastic and viscous
    # contributions cancel in the instantaneous limit.
    drag = np.where(active, 1., 1/off_rate)
    force = np.sum(offsets + drag*(offsets-old)/dt)/active.sum()
    lower = np.sum((1+birth-offsets)[active])/active.sum()
    return float(force), float(lower), int(continuously_present.sum())


def main():
    prior = json.loads(PRIOR.read_text())["runs"]["mouth_contact_crossing_control01"]
    rows = [json.loads(s) for s in prior["stdout_lines"]]
    flat, = [r for r in rows if r['test']=='stationary_patch']
    p, tau = flat['inferred_plateau'], flat['inferred_tau']
    lam = p/tau
    emit(test="intact_input", source_sha256=sha256(PRIOR.read_bytes()).hexdigest(),
         p=p, tau_on=tau, tau_off_same_support=1/lam)
    for r in rows:
        if r["test"] != "flat_return" or not r["mobile"]:
            continue
        measured = r["normalized_gains"][0]
        a = attached_factor(BASE_GAP, p, tau)
        moment = measured/a
        emit(test="archived_return", steps=r["steps"], gap=BASE_GAP, attached_factor=a,
             predicted=a*np.exp(-2*lam*BASE_GAP), previous_reading=measured,
             gap_moment=moment, effective_single_gap_retention=np.sqrt(moment),
             same_support_retention=np.exp(-lam*BASE_GAP),
             opposite_face_discrepancy=r["support_minus_link_gain"])
        if r["steps"] == 2048:
            observed = moment
    # This alternative deliberately fits one archived return by mixing two
    # subpopulations with different detached drags and identical attached laws.
    rates = np.array([.5, 2.])
    at_reference = np.exp(-2*rates*BASE_GAP)
    slow_weight = (observed-at_reference[1])/(at_reference[0]-at_reference[1])
    if not 0 < slow_weight < 1:
        raise ArithmeticError("No positive matching mixture")
    weights = np.array([slow_weight, 1-slow_weight])
    emit(test="alternative_fit", rates=rates.tolist(), weights=weights.tolist())
    for steps, sites in ((2048, 1024), (8192, 4096)):
        first_moment = {}
        for multiple in (1, 2, 4):
            gap = multiple*BASE_GAP
            a = attached_factor(gap,p,tau)
            force, link, retained = zip(*(support_reading(gap, rate, steps, sites)
                                          for rate in (1., *rates)))
            for label, predicted_moment, measured, opposite in (
                    ("same_support", np.exp(-2*lam*gap), force[0], link[0]),
                    ("two_drags", weights @ np.exp(-2*rates*gap),
                     weights @ force[1:], weights @ link[1:])):
                inferred = measured/a
                if multiple == 1:
                    first_moment[label] = inferred
                spread = inferred-first_moment[label]**multiple
                emit(test="new_return", law=label, steps=steps, sites_per_species=sites,
                     gap=gap, gap_multiple=multiple, retained=max(retained),
                     attached_factor=a, predicted=a*predicted_moment, measured=measured,
                     link_force=opposite, measured_moment=inferred,
                     difference_from_single_rate=spread)
                if max(retained) or abs(opposite-measured)>.0005 or abs(a*predicted_moment-measured)>.0007:
                    raise ArithmeticError("Contact-event or traction cross-check failed")


if __name__ == "__main__":
    main()
