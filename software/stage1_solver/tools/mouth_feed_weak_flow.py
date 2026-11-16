"""Steady neutral carrier on the feeding side of one periodic, flat aperture.

Radius and inward number rate are prescribed. The two mouth density laws are
test inputs. The nonlinear solves share the same feed density and phase drive.
All three brane coordinates are periodic. Only 0 <= z <= 1 is modeled, with
z increasing from the feed to the mouth. m = hbar = 1, P = n^5/5.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.fft import irfftn, rfftn


class FeedPatch:
    def __init__(self, n, depth_points, law):
        self.dx = 4 / n
        self.dz = 1 / (depth_points - 1)
        self.n = n
        self.law = law
        coord = (np.arange(n) - n / 2 + 0.5) * self.dx
        self.x, y, z = np.meshgrid(coord, coord, coord, indexing="ij", sparse=True)
        r2 = self.x**2 + y**2 + z**2
        self.mouth = r2 < 1
        self.weight = np.maximum(1 - r2, 0)**3
        self.current = 0.25 * self.weight / (self.dx**3 * self.weight.sum())
        self.restoring = self.weight if law == "restoring_slope" else np.zeros_like(self.weight)
        self.shape = (depth_points, n, n, n)
        # Eigenvalues of the positive, periodic second difference.
        k = 4 * np.sin(np.pi * np.arange(n) / n)**2 / self.dx**2
        self.lateral_diag = k[:, None, None] + k[None, :, None] + k[None, None, :n//2+1]
        self.mu0 = 0.25
        self.shift = 2.0

    def x_derivative(self, f):
        return (np.roll(f, -1, axis=-3) - np.roll(f, 1, axis=-3)) / (2 * self.dx)

    def elliptic_step(self, rhs, derivative):
        # At the feed, the complex envelope is 1. At the mouth, the
        # centred ghost formula enforces its supplied normal derivative.
        h = 0.5 / self.dz**2
        system = rhs[1:].copy()
        system[0] += h
        system[-1] += derivative / self.dz
        # The RHS is complex. Transform real and imaginary parts separately
        # using the real-FFT path for all three periodic brane coordinates.
        def solve_real(part):
            work = rfftn(part, axes=(1, 2, 3))
            diag = self.shift + self.lateral_diag / 2 + 2 * h
            upper = []
            for j in range(len(work)):
                lower = -2 * h if j == len(work)-1 else -h
                pivot = diag if j == 0 else diag - lower * upper[-1]
                if j:
                    work[j] -= lower * work[j-1]
                work[j] /= pivot
                upper.append(-h / pivot)
            for j in range(len(work)-2, -1, -1):
                work[j] -= upper[j] * work[j+1]
            return irfftn(work, s=(self.n,)*3, axes=(1, 2, 3))
        solution = solve_real(system.real) + 1j*solve_real(system.imag)
        result = np.empty_like(rhs)
        result[0] = 1
        result[1:] = solution
        return result

    def solve(self, speed, reference_density=None):
        # Gauge out exp(i U x). The feed then has phase gradient U exactly;
        # choose mu(U)=mu0+U^2/2 to keep its density and thermodynamic state.
        f = np.ones(self.shape, complex)
        for iteration in range(1, 251):
            n = abs(f)**2
            rhs = (self.shift + self.mu0 - self.mu0*n**4)*f + 1j*speed*self.x_derivative(f)
            # The phase slope j/n enforces the given inward rate. The real
            # amplitude slope is zero or relaxes departures from feed density.
            wall = f[-1]
            amp = abs(wall)
            log_amp_z = -self.restoring*(amp-1)/amp
            if reference_density is not None:
                # Keep the zero-drive normal-speed profile. A single scalar
                # rescales it to preserve the total number rate as n changes.
                raw_current = self.weight * n[-1] / reference_density
                self.current = 0.25 * raw_current / (self.dx**3 * raw_current.sum())
            derivative = (log_amp_z + 1j*self.current/n[-1])*wall
            updated = self.elliptic_step(rhs, derivative)
            error = np.max(abs(updated - f))
            if not np.isfinite(error):
                raise ArithmeticError("Nonfinite field iteration")
            f = updated
            if error < 2e-11:
                return f, iteration, error
        raise ArithmeticError(f"Steady nonlinear iteration failed, final update {error}")

    def read(self, f, speed):
        # One-sided derivative from the feed. Reconstruct the actual normal
        # flux and mixed stress from the field, independently of input j.
        face = f[-1]
        dzf = (3*face - 4*f[-2] + f[-3]) / (2*self.dz)
        n = abs(f)**2
        j = (face.conj()*dzf).imag
        ux = speed + (face.conj()*self.x_derivative(face)).imag/n[-1]
        logn_z = -2*self.restoring*(1-1/abs(face))
        stress = -0.25*n[-1]*self.x_derivative(logn_z)
        logn = np.log(n)
        slope_fd = (11*logn[-1] - 18*logn[-2] + 9*logn[-3] - 2*logn[-4])/(6*self.dz)
        stress_fd = -0.25*n[-1]*self.x_derivative(slope_fd)
        def integral(field):
            return self.dx**3 * field[self.mouth].sum()
        # The discrete probability flux on a link follows directly from the
        # imaginary part of the stationary nearest-neighbour field equation.
        feed_rate = self.dx**3 * (f[0].conj()*f[1]).imag.sum() / self.dz
        return dict(rate=integral(j), carried=integral(self.current*ux), stress=integral(stress),
                    carried_fd=integral(j*ux), stress_fd=integral(stress_fd), feed_rate=feed_rate,
                    density_dipole=integral(n[-1]*np.broadcast_to(self.x, self.mouth.shape)),
                    peak_current_error=np.max(abs(j-self.current)), feed_speed=speed,
                    min_density=np.min(n))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", type=int, nargs="+", default=[16, 24, 32])
    parser.add_argument("--depth-cells", type=int, help="Override normal resolution")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for n in args.levels:
        for law in ("flat_density_slope", "restoring_slope"):
            depth_cells = args.depth_cells if args.depth_cells is not None else n//2
            if n < 4 or depth_cells < 3:
                raise ValueError("At least four lateral and three normal intervals are required")
            patch = FeedPatch(n, depth_cells+1, law)
            control = None
            for speed in (0.0, 0.01, 0.02):
                f, steps, error = patch.solve(speed)
                data = patch.read(f, speed)
                np.savez_compressed(args.output/f"n{n}_{law}_u{speed:g}.npz",
                                    field=f, current=patch.current, dx=patch.dx, dz=patch.dz,
                                    speed=speed, mouth=patch.mouth)
                if speed == 0:
                    control = data
                slopes = [(data[key]-control[key])/speed if speed else 0.0
                          for key in ("carried", "stress", "carried_fd", "stress_fd", "density_dipole")]
                print(f"n={n} law={law} U={speed:g} steps={steps} update={error:.3g} "
                      f"feed_rate={data['feed_rate']:.10g} mouth_rate_fd={data['rate']:.10g} "
                      f"max_j_error={data['peak_current_error']:.3g} "
                      f"min_n={data['min_density']:.9g} carried={data['carried']:.10g} "
                      f"stress={data['stress']:.10g} slopes(flow,stress,flow_fd,stress_fd,density_moment)="
                      + ",".join(f"{x:.10g}" for x in slopes), flush=True)
                if (data["rate"] <= 0 or data["min_density"] <= 0
                        or abs(data["feed_rate"] - 0.25) > 2e-7):
                    raise ArithmeticError("Lost steady drainage, flux balance or positive density")


if __name__ == "__main__":
    main()
