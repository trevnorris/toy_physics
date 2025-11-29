import numpy as np


class OnePNOrbitIntegrator:
    """
    Fast 2-body orbit integrator using the scalar 1/2 toy 1PN force law:
        a(r) = -[ μ/r^2 + 3 μ^2/(c_s^2 r^3) ] r_hat
    """

    def __init__(self, mu, c_s):
        self.mu = mu
        self.c_s = c_s

    def acceleration(self, pos, vel):
        """
        pos, vel: arrays of shape (N, 3) in simulation units.
        Returns acceleration array of shape (N, 3).
        """
        r_vec = pos
        r2 = np.sum(r_vec**2, axis=1)
        r = np.sqrt(r2)
        eps = 1e-12
        r_safe = np.maximum(r, eps)
        r_hat = r_vec / r_safe[:, None]

        a_mag = -(self.mu / r_safe**2 + 3 * (self.mu**2) / (self.c_s**2 * r_safe**3))
        return a_mag[:, None] * r_hat

    def step_leapfrog(self, pos, vel, dt):
        """
        Advance one leapfrog step.
        """
        acc = self.acceleration(pos, vel)
        vel_half = vel + 0.5 * dt * acc
        pos_new = pos + dt * vel_half
        acc_new = self.acceleration(pos_new, vel_half)
        vel_new = vel_half + 0.5 * dt * acc_new
        return pos_new, vel_new
