"""Read the local incident speed in the target-absent periodic feed reference."""

import numpy as np

from mouth_feed_weak_flow import FeedPatch


def main():
    # Same periodic box, depth mesh and feed as the signed held-mouth budget.
    # Remove both aperture boundary terms. Keep its footprint only for sampling.
    patch = FeedPatch(32, 65, "restoring_slope")
    patch.current.fill(0)
    patch.restoring = np.zeros_like(patch.restoring)
    footprint = patch.mouth[1:-1]
    if np.any(patch.mouth[[0, -1]]) or not footprint.any():
        raise ValueError("Footprint must lie inside the centered phase stencil")

    for speed in (0.0, 0.01, -0.01, 0.02, -0.02):
        envelope, steps, _ = patch.solve(speed)
        # Restore the physical phase. The centered phase increment measures the
        # speed at the sampling plane independently of the solver's ux readout.
        face = envelope[-1] * np.exp(1j * speed * patch.x)
        phase_step = np.angle(face[2:] * face[:-2].conj())
        ux = (phase_step / (2 * patch.dx))[footprint]
        density_error = np.max(np.abs(np.abs(envelope)**2 - 1))
        normal_link = np.imag(envelope[-2].conj() * envelope[-1]) / patch.dz
        max_normal = np.max(np.abs(normal_link))
        print(f"feed={speed:+.4f} local_min={ux.min():+.12g} local_max={ux.max():+.12g} "
              f"max_speed_error={np.max(np.abs(ux-speed)):.3g} "
              f"max_density_error={density_error:.3g} max_normal_link={max_normal:.3g} "
              f"steps={steps}", flush=True)
        np.testing.assert_allclose(ux, speed, rtol=0, atol=1e-12)
        np.testing.assert_allclose([density_error, max_normal], 0, rtol=0, atol=1e-11)


if __name__ == "__main__":
    main()
