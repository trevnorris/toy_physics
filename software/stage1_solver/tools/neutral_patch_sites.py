"""Zero-charge, zero-intake material patch with mobile sites on two parallel planes.

All pairs on/in between the planes interact. No sites, bonds, edges of a rigid
pad, or registry phases are prescribed in the tangent plane. The planes and
the finite number of constituents on each one are part of the hypothesis.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.csgraph import connected_components


COUNT = 12
GAP = 1.
PITCH_UNIT = 1.  # Isolated same-plane preferred distance; the patch pitch is measured.
WITHIN = (1., 7., PITCH_UNIT)  # well energy, inverse range, rest length
BETWEEN = (.6, 10., .94)
STEP = .0005
SAMPLES = 4000
TEMPERATURE = .008  # k_B T in the within-plane well-energy unit
LAWS = ("central", "layered")


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def morse(r, depth, alpha, length):
    weight = np.exp(-alpha*(r-length))
    return depth*weight*(weight-2), 2*depth*alpha*weight*(1-weight)


class Patch:
    def __init__(self, law):
        if law not in LAWS:
            raise ValueError("Unknown interplane interaction")
        self.law = law
        self.i, self.j = np.triu_indices(2*COUNT, 1)
        self.cross = (self.i < COUNT) & (self.j >= COUNT)
        self.depth, self.alpha, self.length = [
            np.where(self.cross, cross, same) for cross, same in zip(BETWEEN, WITHIN, strict=True)]

    def evaluate(self, points, removed=None):
        delta = points[self.j]-points[self.i]
        tangent_sq = np.einsum("ij,ij->i", delta, delta)
        radius_sq = tangent_sq.copy()
        projection = np.ones_like(radius_sq)
        if self.law == "layered":
            radius_sq[self.cross] = tangent_sq[self.cross]**2/PITCH_UNIT**2
            projection[self.cross] = 2*tangent_sq[self.cross]/PITCH_UNIT**2
        radius_sq[self.cross] += GAP**2
        radius = np.sqrt(radius_sq)
        if np.any(radius == 0):
            raise ArithmeticError("Coincident in-plane particles")
        energy, radial_derivative = morse(radius, self.depth, self.alpha, self.length)
        if removed is not None:
            energy[removed] = radial_derivative[removed] = 0.
        gradients = delta*(projection*radial_derivative/radius)[:,None]
        nodal = np.stack([np.bincount(self.j, weights=g, minlength=2*COUNT) -
                          np.bincount(self.i, weights=g, minlength=2*COUNT)
                          for g in gradients.T], axis=1)
        return energy, nodal, np.sqrt(tangent_sq)

    def relaxed(self, positions, fixed=()):
        free = np.ones(2*COUNT, dtype=bool)
        free[list(fixed)] = False
        state = positions.copy()

        def objective(flat):
            state[free] = flat.reshape(-1, 2)
            energy, gradient, _ = self.evaluate(state)
            return float(sum(energy)), gradient[free].ravel()

        result = minimize(objective, state[free].ravel(), method="L-BFGS-B", jac=True,
                          options={"maxiter": 5000, "gtol": 3e-9, "ftol": 1e-15, "maxls": 70})
        _, gradient = objective(result.x)
        error = float(np.max(abs(gradient)))
        if not result.success or error > 8e-6:
            raise ArithmeticError(f"Patch failed to equilibrate: {result.message}, residual {error}")
        return state.copy(), error

    def neighbors(self, energy, distances, fraction):
        # Every sufficiently deeply bound same-plane pair is in the network.
        # Between planes a registered link also has to be mutual nearest in x,y.
        pair_distance = distances[self.cross].reshape(COUNT, COUNT)
        lower = np.arange(COUNT)
        upper = pair_distance.argmin(axis=1)
        mutual = pair_distance.argmin(axis=0)[upper] == lower
        cross_links = np.zeros((COUNT, COUNT), dtype=bool)
        cross_links[lower[mutual], upper[mutual]] = True
        active = energy < -fraction*self.depth
        active[self.cross] &= cross_links.ravel()
        return active

    def pitch(self, distances):
        same = ~self.cross
        nearest = np.full(2*COUNT, np.inf)
        np.minimum.at(nearest, self.i[same], distances[same])
        np.minimum.at(nearest, self.j[same], distances[same])
        return float(np.median(nearest))

    def component_count(self, edges):
        graph = np.zeros((2*COUNT, 2*COUNT), dtype=bool)
        graph[self.i[edges], self.j[edges]] = True
        return int(connected_components(graph, directed=False, return_labels=False))


def initial(seed):
    rng = np.random.default_rng(seed)
    planes = []
    for _ in range(2):
        points = []
        while len(points) < COUNT:
            trial = rng.uniform(-2.1, 2.1, 2)
            if trial@trial < 2.1**2 and all(np.linalg.norm(trial-p) > .62 for p in points):
                points.append(trial)
        planes.extend(points)
    return np.array(planes)


def prepare(patch, seed):
    raw = initial(seed)
    start, residual = patch.relaxed(raw)
    raw_energy = patch.evaluate(raw)[0].sum()
    energy, gradient, distances = patch.evaluate(start)
    final_energy = energy.sum()
    fractions = (.4, .5, .6)
    survivors = {f: patch.neighbors(energy, distances, f) for f in fractions}
    counts = {f: int(sum(edges & patch.cross)) for f, edges in survivors.items()}
    initial_components = patch.component_count(survivors[.5])
    eligible = np.flatnonzero(survivors[.6] & patch.cross)
    tagged = None
    if len(eligible):
        mid = (start[patch.i[eligible]] + start[patch.j[eligible]])/2
        tagged = int(eligible[np.linalg.norm(mid-start.mean(axis=0), axis=1).argmin()])
    # Brownian dynamics at unit mobility. Record identity at every integration
    # step. The tag is chosen at the first minimum, before the survival test.
    rng = np.random.default_rng(seed+10000)
    state = start.copy()
    noise = np.sqrt(2*TEMPERATURE*STEP)
    for _ in range(SAMPLES):
        state += -STEP*gradient + noise*rng.normal(size=state.shape)
        energy, gradient, distances = patch.evaluate(state)
        for f in fractions:
            survivors[f] &= patch.neighbors(energy, distances, f)
    rested, final_residual = patch.relaxed(state)
    e_end, _, d_end = patch.evaluate(rested)
    for f in fractions:
        survivors[f] &= patch.neighbors(e_end, d_end, f)
    summary = dict(seed=seed, law=patch.law, constituents_per_layer=COUNT,
                   initial_energy=float(raw_energy), minimum_energy=float(final_energy),
                   measured_pitch=patch.pitch(d_end),
                   tag=None if tagged is None else [int(patch.i[tagged]), int(patch.j[tagged])],
                   tag_offset_before_clamp=None if tagged is None else float(d_end[tagged]),
                   observed_model_time=SAMPLES*STEP, dt=STEP, kBT=TEMPERATURE,
                   initial_connected_components=initial_components,
                   persistent_components=patch.component_count(survivors[.5]),
                   interplane_contacts={str(f): [counts[f], int(sum(edges & patch.cross))]
                                        for f, edges in survivors.items()},
                   tag_survived={str(f): None if tagged is None else bool(edges[tagged])
                                 for f, edges in survivors.items()},
                   max_free_gradient=max(residual, final_residual))
    emit(test="material_patch", **summary)
    return rested, tagged, summary


def radial_calibration():
    # Independent isolated neutral pair, normal extension only. The virtual
    # work instrument reads potential energies on opposite sides of each gap.
    grid = np.linspace(.90, 1.40, 251)
    h = 1e-5
    high = morse(grid+h, *BETWEEN)[0]
    low = morse(grid-h, *BETWEEN)[0]
    return grid, (high-low)/(2*h)


def calibrate(path):
    grid, tension = radial_calibration()
    data = {"normal_separations": grid.tolist(), "normal_tensions": tension.tolist(),
            "force_unit_measured": float(tension.max()), "gap": GAP,
            "layered_crossover_length": PITCH_UNIT, "patches": {}}
    for law in LAWS:
        _, _, result = prepare(Patch(law), 191)
        data["patches"][law] = result
    path.write_text(json.dumps(data, allow_nan=False)+"\n")
    emit(test="independent_tensile_calibration", samples=len(grid),
         normal_tension_at_gap=float(tension[np.argmin(abs(grid-GAP))]),
         sampled_maximum_tension=data["force_unit_measured"], gap=GAP)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calibration", required=True, type=Path)
    calibrate(parser.parse_args().calibration)
