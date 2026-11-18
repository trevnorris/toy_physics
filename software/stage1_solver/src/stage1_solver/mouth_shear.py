"""Outgoing rotational displacement around a spherical rim in homogeneous brane.

The contact used here is an explicit diagnostic hypothesis. Its tangential
traction per area is slip_drag times the relative tangential velocity. The
scalar carrier and the longitudinal material sector do not drive this field.
"""

from math import isfinite, pi

from scipy.integrate import solve_ivp


class MouthShearField:
    """Linear translation dipole from rest; positive force acts on the exterior."""

    def __init__(self, radius, density, wave_speed, slip_drag, mouth_velocity, end_time):
        params = (radius, density, wave_speed, end_time)
        if not all(isfinite(x) and x > 0 for x in params):
            raise ValueError("Radius, density, wave speed and duration must be positive")
        if not isfinite(slip_drag) or slip_drag < 0:
            raise ValueError("Slip drag must be finite and nonnegative")
        self.radius = radius
        self.wave_speed = wave_speed
        self.modulus = density * wave_speed**2
        self.tau = radius / wave_speed
        self.mass_scale = 8 * pi * density * radius**3 / 3
        self.drag = 8 * pi * radius**2 * slip_drag / 3
        self.mouth_velocity = mouth_velocity
        self.end_time = end_time
        # B, B', B'' are the retarded outgoing dipole amplitude and its
        # derivatives. Solving at the rim determines the entire exterior.
        result = solve_ivp(self._rhs, (0, end_time), (0., 0., 0.), method="DOP853",
                           dense_output=True, rtol=1e-11, atol=1e-24)
        if not result.success:
            raise ArithmeticError(result.message)
        self._amplitude = result.sol

    def _rhs(self, t, state):
        _, first, second = state
        tau = self.tau
        scale = self.mass_scale + self.drag * tau
        third = (self.drag * (self.mouth_velocity(t) - first) / scale - second) / tau
        return first, second, third

    def displacement(self, t, gap):
        """Return p,d where u_r=p cos(theta), u_theta=-d sin(theta) at r=a+gap."""
        if not (0 <= t <= self.end_time and isfinite(gap) and gap >= 0):
            raise ValueError("Readout outside the solved time interval or exterior")
        past = t - gap / self.wave_speed
        if past <= 0:
            return 0., 0.
        radius = self.radius + gap
        delay_scale = radius / self.wave_speed
        value, first, second = self._amplitude(past)
        common = value + delay_scale * first
        factor = (self.radius / radius)**3
        return -2 * factor * common, factor * (common + delay_scale**2 * second)

    def rim_velocity(self, t):
        if not 0 <= t <= self.end_time:
            raise ValueError("Time outside the solved interval")
        state = self._amplitude(t)
        first, second, third = self._rhs(t, state)
        return first + self.tau * second + self.tau**2 * third

    def shell_force(self, t, gap, spacing):
        """Axial curl traction integrated on a cut, from sampled material only."""
        if not 0 < spacing <= gap:
            raise ValueError("The centered stencil must lie in the brane exterior")
        radius = self.radius + gap
        radial, _ = self.displacement(t, gap)
        _, left = self.displacement(t, gap - spacing)
        _, right = self.displacement(t, gap + spacing)
        gradient = ((radius + spacing) * right - (radius - spacing) * left) / (2 * spacing)
        return 8 * pi * self.modulus * radius * (radial - gradient) / 3
