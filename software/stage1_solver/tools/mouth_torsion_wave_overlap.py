"""Reciprocal virtual work between a prepared mouth and imposed acoustic modes."""

import numpy as np

from mouth_torsion_contact import cut_graph, mouth_response
from mouth_torsion_power import spectrum


def main():
    for radius in (5, 10):
        sites, edges, _, _ = cut_graph(radius)
        origin, end = sites[edges[:, 0]], sites[edges[:, 1]]
        axes = end - origin
        tau, _ = mouth_response(radius, 1)
        directed = tau[:, None] * axes
        moment = .5 * np.cross((origin + end) / 2, directed).sum(axis=0)
        print(f"R={radius} circulating_torque_moment={np.round(moment, 9).tolist()}")
        for ray in ((1, 0, 0), (1, 2, 3)):
            direction = np.array(ray) / np.linalg.norm(ray)
            leading = np.linalg.norm(np.cross(direction, moment))
            for kr in (.2, .1, .05):
                k = kr / radius
                # The bulk eigenvectors store each outgoing edge at its base.
                source = np.exp(-1j * k * (origin @ direction)) @ directed
                _, modes = spectrum(k * direction)
                work = np.linalg.norm(modes[:3, 4:6].conj().T @ source)
                print(f"R={radius} ray={ray} kR={kr:g} "
                      f"acoustic_work_norm={work:.12g} work_over_k={work/k:.9f} "
                      f"moment_prediction={leading/2:.9f}")
                # I+2J = 4. An arbitrary basis in the degenerate subspace
                # changes the two individual overlaps, but preserves the norm.
                if kr == .05:
                    assert leading > 0 and abs(work/k-leading/2) < .005*leading/2


if __name__ == "__main__":
    main()
