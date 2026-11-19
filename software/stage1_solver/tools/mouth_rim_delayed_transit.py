"""Prepare, wait and probe a rim site with its bulk channel open throughout.

Every site including the parcel is hard-core. In an open path the one- and
two-excitation propagators are U and its antisymmetric second power. All
energies are in units of the bulk hopping, times in hbar / hopping.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal


BAND_CENTER = 4.0
PUMP_FRACTION = 0.8
DWELL = 0.33
CONTACT = np.sqrt(2.0)


class Patch:
    def __init__(self, bulk_sites, rim_energy, parcel_energy):
        self.parcel_energy = parcel_energy
        self.diagonal = np.full(bulk_sites + 2, BAND_CENTER)
        self.diagonal[:2] = parcel_energy, rim_energy
        self.free_edges = np.ones(bulk_sites + 1)
        self.free_edges[0] = 0
        self.levels, self.basis = eigh_tridiagonal(self.diagonal, self.free_edges)
        self.parcel = np.eye(1, len(self.diagonal))[0]
        detuning = rim_energy - BAND_CENTER
        self.bound = None
        if abs(detuning) > 1:
            profile = (1 / detuning) ** np.arange(bulk_sites + 1)
            self.bound = np.r_[0.0, profile / np.linalg.norm(profile)]

    def free(self, time, vector):
        return self.basis @ (np.exp(-1j * self.levels * time) * (self.basis.T @ vector))

    def pulse(self, duration, strength=CONTACT):
        edges = self.free_edges.copy()
        edges[0] = strength
        levels, basis = eigh_tridiagonal(self.diagonal, edges)
        active = (basis * np.exp(-1j * levels * duration)) @ basis.T
        undo_free = (self.basis * np.exp(1j * self.levels * duration)) @ self.basis.T
        return active, undo_free @ active

    def energy(self, v):
        return np.sum(self.levels * abs(self.basis.T @ v)**2)

    def ledger(self, u, v, strength):
        out = u @ v
        gain = self.energy(out) - self.energy(v)
        release = 2 * strength * (np.conj(v[0]) * v[1] - np.conj(out[0]) * out[1]).real
        return gain, release


def fringe(echo, state, weight, coherence_squared=0.0):
    """Coherences relative to free flight, then subtract the empty-memory phase."""
    a = echo[0, 0]
    b = echo[0] @ state
    c = np.vdot(state, echo[:, 0])
    d = np.vdot(state, echo @ state)
    ground = 1 + weight * (d - 1)
    excited = (1 - weight) * a + weight * (a * d - b * c)
    determinant = ground * excited - coherence_squared * b * c
    return np.angle(determinant * a.conjugate()) / 2, abs(determinant), ground, excited


def calculate(length, rim_energy, parcel_energy=3.5):
    patch = Patch(length, rim_energy, parcel_energy)
    active, echo = patch.pulse(DWELL)
    _, half = patch.pulse(DWELL / 2)
    departed = active[:, 0].copy()
    departed[0] = 0
    deposited = np.vdot(departed, departed).real
    weight = PUMP_FRACTION * deposited
    state = departed / np.sqrt(deposited)
    actual, release = patch.ledger(active, patch.parcel, CONTACT)
    prefix = f"bulk={length} rim={rim_energy:g} parcel={parcel_energy:g}"
    print(f"prepare {prefix} receiver_quanta={weight:.12g} "
          f"free_energy_gain={PUMP_FRACTION*actual:.12g} balance_error={abs(actual-release):.3g}")
    assert abs(actual - release) < 2e-12
    series = []
    for delay in (0, 1, 4, 16, 40):
        delayed = patch.free(delay, state)
        phase, visibility, _, _ = fringe(echo, delayed, weight)
        coherent, _, _, _ = fringe(echo, delayed, weight, weight*(1-PUMP_FRACTION))
        half_phase = fringe(half, delayed, weight)[0]
        gain, edge = patch.ledger(active, delayed, CONTACT)
        # For an independently mixed probe, energy is additive in occupied one-particle orbitals.
        probe_gain = PUMP_FRACTION * actual + weight * gain
        probe_edges = PUMP_FRACTION * release + weight * edge
        print(f"delay {prefix} wait={delay} local={weight*abs(delayed[1])**2:.12g} "
              f"phase={phase:.12g} pure_pump_phase={coherent:.12g} half_dwell_phase={half_phase:.12g} "
              f"contrast_product={visibility:.12g} probe_energy_gain={probe_gain:.12g} "
              f"balance_error={abs(probe_gain-probe_edges):.3g}")
        assert abs(probe_gain - probe_edges) < 2e-12 and abs(phase) < .2
        series.append((phase, coherent, half_phase, weight*abs(delayed[1])**2))
    n_bound, asymptotic, stationary = 0., 0., 0.
    if patch.bound is not None:
        bound = patch.bound
        n_bound = PUMP_FRACTION * abs(np.vdot(bound, departed))**2
        asymptotic = fringe(echo, bound, n_bound)[0]
        survival = abs(np.vdot(bound, active @ bound))**2
        injection = abs(np.vdot(bound, active[:, 0]))**2
        stationary = n_bound / (1 - survival)
        for passes in (1, 6, 20):
            accumulated = stationary * (1 - survival**passes)
            print(f"separated {prefix} passes={passes} occupation={accumulated:.12g} "
                  f"local={accumulated*abs(bound[1])**2:.12g} phase={fringe(echo,bound,accumulated)[0]:.12g}")
        print(f"train_limit {prefix} bound_survival_per_contact={survival:.12g} "
              f"capture_per_excited_arrival={injection:.12g} bound_population={stationary:.12g} "
              f"phase={fringe(echo,bound,stationary)[0]:.12g}")
        assert 0 <= stationary <= PUMP_FRACTION + 1e-12
    print(f"long_delay {prefix} bound_population={n_bound:.12g} "
          f"local={n_bound*(abs(patch.bound[1])**2 if patch.bound is not None else 0.):.12g} "
          f"phase={asymptotic:.12g}")
    return np.r_[np.ravel(series), n_bound, asymptotic, stationary]


if __name__ == "__main__":
    for rim, parcel in ((4.,3.5), (5.5,3.5), (5.5,7.5)):
        coarse = calculate(144, rim, parcel)
        fine = calculate(224, rim, parcel)
        difference = np.max(abs(coarse - fine))
        print(f"length_check rim={rim:g} parcel={parcel:g} max_difference={difference:.3g}")
        assert difference < 2e-10
