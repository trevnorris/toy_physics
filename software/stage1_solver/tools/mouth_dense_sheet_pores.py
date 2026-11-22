"""Bound the static test-site energy on a frozen 4D D4 coupon's ordered sheet.

Distinct from a finite-temperature or fully mobile passage free energy.
"""
import argparse
import heapq
import itertools
import json
import math

import numpy as np


def measure(density=1.65, tolerance=0.1):
    spacing = (0.5 / density) ** 0.25
    cutoff = 2.5
    # Tangential positions are in lattice steps, with all normal positions at
    # multiples of six steps. Two-step translations and independent reflections
    # map the full sheet into [0,1]^3 while preserving FCC parity.
    span = math.ceil(cutoff / spacing)
    nodes = np.array([p for p in itertools.product(range(-span, 2 + span), repeat=3)
                      if sum(p) % 2 == 0], dtype=float)

    def potential(r):
        e = np.exp(-7 * (r - 1))
        t = np.clip(4 * (cutoff - r), 0., 1.)
        return e * (e - 2) * (t**3 * (10 + t * (-15 + 6*t)))

    def plane_energy(point, w=0.):
        tangential2 = ((nodes - point)**2).sum(axis=1)
        return float(sum(np.sum(potential(spacing * np.sqrt(tangential2 + (w-k*6)**2)))
                         for k in (-1, 0, 1)))

    def lower_bound(low, high):
        # Each held site's distance varies between two componentwise bounds.
        # The switched Morse energy decreases to its minimum at r=1 and then
        # increases to zero, so clamp r=1 to that interval. On the sheet all
        # normal replicas lie beyond cutoff.
        near_vec = np.maximum(np.maximum(low-nodes, nodes-high), 0)
        far_vec = np.maximum(np.abs(nodes-low), np.abs(nodes-high))
        rnear = spacing * np.sqrt((near_vec**2).sum(axis=1))
        rfar = spacing * np.sqrt((far_vec**2).sum(axis=1))
        return float(potential(np.clip(1., rnear, rfar)).sum())

    octa = np.array([1., 0., 0.])
    tetra = np.full(3, .5)
    best = plane_energy(octa)
    witness = octa.copy()
    low = np.zeros(3)
    high = np.ones(3)
    boxes = [(lower_bound(low, high), 0, low, high)]
    next_id = 1
    while best - boxes[0][0] > tolerance:
        _, _, low, high = heapq.heappop(boxes)
        mid = (low + high) / 2
        for corner in itertools.product((0, 1), repeat=3):
            left = np.where(corner, mid, low)
            right = np.where(corner, high, mid)
            bound = lower_bound(left, right)
            center = (left + right) / 2
            value = plane_energy(center)
            if value < best:
                best = value
                witness = center.copy()
            heapq.heappush(boxes, (bound, next_id, left, right))
            next_id += 1

    # In lattice-step units the FCC Voronoi cell has |x|+|y|<=1,
    # |y|+|z|<=1, |x|+|z|<=1. Its vertices are the signed permutations
    # of (1,0,0) and (1/2,1/2,1/2); the largest vertex radius is 1.
    return {"density": density, "spacing": spacing, "temperature": 1.024,
            "sheet_test_particle_nearest_site_max_distance": spacing,
            "single_pair_energy_at_max_distance": float(potential(spacing)),
            "octahedral_energy": plane_energy(octa),
            "tetrahedral_energy": plane_energy(tetra),
            "sheet_min_energy_lower": boxes[0][0],
            "sheet_min_energy_upper": best, "witness_lattice_units": witness.tolist(),
            "energy_tolerance": tolerance, "boxes_evaluated": next_id,
            "sheet_only_path_at_octahedral_xyztangent": [
                {"w_in_lattice_steps": w, "energy": plane_energy(octa, w)}
                for w in [0., .25, .5, .75, 1., 1.5, 2., 3.]]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--density", type=float, default=1.65)
    parser.add_argument("--tolerance", type=float, default=.1)
    args = parser.parse_args()
    if args.density <= 0 or args.tolerance <= 0:
        parser.error("density and tolerance must be positive")
    # The fixed replica and on-sheet bound in this test require 6a>cutoff.
    if 6*(0.5/args.density)**.25 <= 2.5:
        parser.error("the surface-normal period must exceed the potential's cutoff")
    print(json.dumps(measure(args.density, args.tolerance), indent=2, allow_nan=False))
