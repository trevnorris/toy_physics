"""Follow contact renewal and weak incremental shear on archived mouth streamlines.

Contact identities follow the carrier in this trial. At the inlet old contacts
have zero reference shear, and new contacts are born relaxed. Contact stiffness
and imposed amplitude divide out of the reported response. Frequencies use the
same model clock for every preparation.
"""
from hashlib import file_digest
from pathlib import Path
import json

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from mouth_free_contact import Carrier
from mouth_free_order_diffuse import read


BASE = Path("_scratch/s11c")
PRIOR = Path("research/pde_ledger_v3/_measurements/mouth_mixed_order_ripple.json")
FREQUENCIES = np.array([0., 1., 4., 12.])
SITE_COUNT = 1_000_000
TRIAL_GAMMA = 1.


def emit(**data):
    print(json.dumps(data, allow_nan=False), flush=True)


def trace_path(fields, steps_per_radius):
    carrier = Carrier(Path(fields["carrier_file"]).parent)
    domain = carrier.domain
    order_path = Path(fields["order_file"])
    with order_path.open("rb") as stream:
        if file_digest(stream, "sha256").hexdigest() != fields["order_sha256"]:
            raise ValueError("Changed order field")
    with np.load(order_path) as data:
        order = data["chi"]

    def derivative(w, rt):
        r = rt[0]
        if w <= .5 or not 0 < r < domain.radius(w)[0]:
            raise ValueError("Path left the post-waist chart")
        psi, pr, pw = read(domain, carrier.field, r, w)
        ur, uw = (pr/psi).imag, (pw/psi).imag
        if uw <= 0:
            raise ValueError("No inward material path")
        return np.array([ur/uw, 1/uw])

    def advance(w, rt, step):
        k1 = derivative(w, rt)
        k2 = derivative(w+.5*step, rt+.5*step*k1)
        k3 = derivative(w+.5*step, rt+.5*step*k2)
        k4 = derivative(w+step, rt+step*k3)
        return rt+step*(k1+2*k2+2*k3+k4)/6

    w = fields["w"]
    rt = np.array([fields["r"], 0.])
    chi = read(domain, order, rt[0], w)[0]
    if abs(chi-.25) > 1e-8:
        raise ValueError("Wrong exit probe")
    points = [[rt[1], rt[0], w, chi]]
    step = -1/steps_per_radius
    while chi < .75:
        next_rt = advance(w, rt, step)
        trial = read(domain, order, next_rt[0], w+step)[0]
        if trial <= chi:
            raise ValueError("One-way loss is inadmissible on this path")
        if trial >= .75:
            def inlet(s):
                return read(domain, order, advance(w, rt, s)[0], w+s)[0]-.75
            step = brentq(inlet, step, 0., xtol=5e-14)
            next_rt = advance(w, rt, step)
        w += step
        rt = next_rt
        chi = read(domain, order, rt[0], w)[0]
        points.append([rt[1], rt[0], w, chi])
        if trial >= .75:
            break
    t, r, w, q = np.array(points[::-1]).T
    return t-t[0], r, w, q


def quadrature(t, q, gamma):
    # The contact population is normalized to the occupied fraction at
    # readout. A contact retains its individual reference until destruction.
    rates = gamma*(1-q)
    exposure = cumulative_trapezoid(rates, t, initial=0.)
    tail = np.exp(exposure-exposure[-1])
    drive_at_birth = np.exp(1j*(t[-1]-t[:, None])*FREQUENCIES)
    reset = np.trapz((rates*tail)[:, None]*drive_at_birth, t, axis=0)
    gain = 1.-reset
    return tail[0], gain


def cohort_sample(t, q, gamma, seed):
    # Individual contacts disappear with equal probability across birth ages.
    # Binomial draws count survivors exactly for a bin of identical contacts.
    # Change the shear reference only on newly occupied sites.
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(t), dtype=np.int64)
    counts[0] = round(SITE_COUNT*q[0])
    for j, dt in enumerate(np.diff(t)):
        mid = (q[j]+q[j+1])/2
        born = gamma*mid*(1-mid)*dt
        removed = born+q[j]-q[j+1]
        p_birth = born/(1-q[j])
        p_loss = removed/q[j]
        if not (0 <= p_birth <= 1 and 0 <= p_loss <= 1):
            raise ValueError("Transition probabilities outside unit interval")
        present = counts[:j+1]
        empty = SITE_COUNT-present.sum()
        present -= rng.binomial(present, p_loss)
        counts[j+1] = rng.binomial(empty, p_birth)
    old = counts[0]/counts.sum()
    birth_time = (t[:-1]+t[1:])/2
    # Every old contact has zero reference. Each new contact carries its
    # measured extension since birth. Complex amplitudes combine two weak
    # drives with independent in-phase and quadrature entrance preparations.
    ext = 1.-np.exp(1j*(t[-1]-birth_time[:, None])*FREQUENCIES)
    gain = (counts[0]+counts[1:] @ ext)/counts.sum()
    return old, counts.sum()/SITE_COUNT, gain


def readings(values):
    return [[float(z.real), float(z.imag)] for z in values]


def run():
    prior = json.loads(PRIOR.read_text())["runs"]["mouth_mixed_order_ripple02"]
    patches = [json.loads(s) for s in prior["stdout_lines"]]
    cases = [p for p in patches if p["test"] == "archive_patch"]
    computed = {}
    for p in cases:
        label = p["profile"], p["state"], p["mesh"]
        coarse = trace_path(p, 1000)
        t, r, w, q = trace_path(p, 2000)
        exposure = np.trapz(1-q, t)
        old, gain = quadrature(t, q, TRIAL_GAMMA)
        old_coarse, gain_coarse = quadrature(coarse[0], coarse[-1], TRIAL_GAMMA)
        error = max(abs(gain-gain_coarse))
        emit(test="trace", profile=label[0], state=label[1], mesh=label[2], intake=p["intake"],
             input_carrier=p["carrier_file"], carrier_sha256=p["carrier_sha256"],
             input_order=p["order_file"], order_sha256=p["order_sha256"],
             inlet=[r[0], w[0], q[0]], outlet=[r[-1], w[-1], q[-1]],
             transit=t[-1], integrated_disorder_time=exposure,
             time_refinement_change=abs(t[-1]-coarse[0][-1]), gain_refinement_change=error,
             old_fraction=old, coarse_old_fraction=old_coarse)
        if error > .002 or max(abs(q[[0, -1]]-[.75, .25])) > 1e-8:
            raise ValueError("Path or response integration failed refinement")
        computed[label] = t, q, exposure
    training = computed["cubic", "near", 48]
    training_old, occupancy, _ = cohort_sample(*training[:2], TRIAL_GAMMA, seed=51493)
    fitted = -np.log(training_old)/training[2]
    emit(test="unforced_contact_calibration", sites=SITE_COUNT, profile="cubic", state="near", mesh=48,
         trial_gamma=TRIAL_GAMMA, old_fraction_of_occupied=training_old, occupancy=occupancy,
         fitted_gamma=fitted, integrated_disorder_time=training[2])
    for i, ((profile, state, mesh), (t, q, exposure)) in enumerate(computed.items()):
        if mesh != 48:
            continue
        for gamma in (0., TRIAL_GAMMA):
            predicted_old, predicted = quadrature(t, q, fitted if gamma else 0.)
            sample_old, occupancy, observed = cohort_sample(t, q, gamma, seed=99998+i)
            emit(test="drive", profile=profile, state=state, mesh=mesh, trial_gamma=gamma,
                 omega=FREQUENCIES.tolist(), predicted_gain=readings(predicted), sampled_gain=readings(observed),
                 gain_magnitude=abs(predicted).tolist(), gain_phase_deg=np.angle(predicted, deg=True).tolist(),
                 predicted_old=predicted_old, sampled_old=sample_old, sampled_occupancy=occupancy,
                 max_complex_error=max(abs(predicted-observed)))
            if max(abs(predicted-observed)) > .01 or abs(occupancy-q[-1]) > .004:
                raise ValueError("Particle sample and independent tag prediction differ")


if __name__ == "__main__":
    run()
