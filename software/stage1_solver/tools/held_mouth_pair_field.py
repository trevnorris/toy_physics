"""Steady neutral carrier feeding one or two held spherical openings in a periodic slab.

The transverse and pair-midpoint reflection planes bisect grid cells. A DCT-II
on each half period implements the even subspace of the periodic difference
operator. The feed is fixed at psi=1. No incident lateral velocity is imposed.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.fft import dctn, idctn


class Slab:
    def __init__(self, width, depth, cells_per_radius, cells_per_depth_unit):
        self.dx = 1 / cells_per_radius
        self.dz = 1 / cells_per_depth_unit
        size = width * cells_per_radius // 2
        if 2*size != width*cells_per_radius:
            raise ValueError("Half period must contain a whole number of cells")
        self.coord = (np.arange(size) + 0.5)*self.dx
        self.intervals = depth*cells_per_depth_unit
        self.shape = (self.intervals+1, size, size, size)
        eig = 2*np.sin(np.pi*np.arange(size)/(2*size))**2 / self.dx**2
        lateral = eig[:, None, None] + eig[None, :, None] + eig[None, None, :]
        self.h = 0.5/self.dz**2
        self.inverse_pivots = {}
        for shift in (0, 2):
            diagonal = shift + lateral + 2*self.h
            pivots = np.empty((self.intervals, *lateral.shape))
            pivots[0] = 1/diagonal
            for i in range(1, self.intervals):
                factor = 2 if i == self.intervals-1 else 1
                pivots[i] = 1 / (diagonal - factor*self.h**2*pivots[i-1])
            self.inverse_pivots[shift] = pivots

    def elliptic(self, rhs, normal, shift):
        # Zero variation at the feed. A centred ghost enforces the specified
        # derivative at the mouth, giving a doubled last lower coefficient.
        rhs[-1] += normal/self.dz
        axes = (1, 2, 3)
        work = dctn(rhs, axes=axes, norm="ortho", overwrite_x=True)
        pivots = self.inverse_pivots[shift]
        work[0] *= pivots[0]
        for i in range(1, self.intervals):
            lower = 2*self.h if i == self.intervals-1 else self.h
            work[i] += lower*work[i-1]
            work[i] *= pivots[i]
        for i in range(self.intervals-2, -1, -1):
            work[i] += self.h*pivots[i]*work[i+1]
        return idctn(work, axes=axes, norm="ortho", overwrite_x=True)

    def preparation(self, separation, alone_density):
        x = self.coord
        offset = separation/2
        r2 = (x[:, None, None]-offset)**2 + x[None, :, None]**2 + x[None, None, :]**2
        weight = np.maximum(1-r2, 0)**3
        if separation == 0:
            return weight, weight
        # Translate the isolated normal-speed shape onto the right mouth.
        # Both mouths have equal rates by midpoint reflection.
        shift = round(offset/self.dx)
        if abs(shift*self.dx-offset) > 1e-12 or offset+1 >= x[-1]:
            raise ValueError("Use a grid-aligned midpoint and a nonoverlapping periodic image")
        i = np.arange(len(x)) - shift
        return weight, weight/alone_density[np.where(i < 0, -i-1, i)]

    def boundary(self, face, weight, speed_shape, separation):
        n = abs(face)**2
        raw = speed_shape*n if separation else speed_shape
        multiplicity = 4 if separation else 8
        current = raw * (0.25/(multiplicity*self.dx**3*raw.sum()))
        normal = (-weight*(1-1/abs(face)) + 1j*current/n)*face
        return current, normal

    def solve(self, separation, alone_density):
        weight, speed_shape = self.preparation(separation, alone_density)
        field = np.ones(self.shape, complex)
        for iteration in range(1, 501):
            _, normal = self.boundary(field[-1], weight, speed_shape, separation)
            bulk = field[1:]
            potential = (abs(bulk)**8-1)/4
            real = self.elliptic(2*(bulk.real-1) - potential*bulk.real, normal.real, 2)
            imag = self.elliptic(-potential*bulk.imag, normal.imag, 0)
            updated = (1+real) + 1j*imag
            update = np.max(abs(updated-bulk))
            if not np.isfinite(update):
                raise ArithmeticError("Nonfinite carrier update")
            field[1:] = updated
            if update < 4e-13:
                current, normal = self.boundary(field[-1], weight, speed_shape, separation)
                return field, weight, current, normal, iteration
        raise ArithmeticError(f"Nonlinear iteration failed, final update {update}")

    def validate(self, f, normal, current, separation):
        # Independent nearest-neighbour real-space equation, including the
        # even lateral ghosts and the mouth ghost. Read every feed-depth link.
        link_error = 0.0
        residual = 0.0
        for k in range(1, len(f)):
            layer = f[k]
            hi = f[k+1] if k+1 < len(f) else f[k-1] + 2*self.dz*normal
            lap = (hi - 2*layer + f[k-1])/self.dz**2
            for axis in range(3):
                pad = [(0, 0)]*3
                pad[axis] = (1, 1)
                lap += np.diff(np.pad(layer, pad, mode="symmetric"), n=2, axis=axis)/self.dx**2
            eq = -lap/2 + (abs(layer)**8-1)*layer/4
            residual = max(residual, float(np.max(abs(eq))))
            flux = 8*self.dx**3*(f[k-1].conj()*layer).imag.sum()/self.dz
            link_error = max(link_error, abs(flux-(0.5 if separation else 0.25)))
        rate = (4 if separation else 8)*self.dx**3*current.sum()
        if (not np.isfinite([residual, link_error, rate]).all() or residual > 2e-8
                or link_error > 2e-7 or abs(rate-0.25) > 2e-14 or np.any(f[0] != 1)):
            raise ArithmeticError(f"Carrier balance failed: {residual=} {link_error=} {rate=}")
        return residual, link_error, rate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, required=True)
    parser.add_argument("--depth", type=int, required=True)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--normal-cells", type=int, required=True, help="Intervals per depth unit")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    slab = Slab(args.width, args.depth, args.cells, args.normal_cells)
    alone_density = None
    for d in (0, 3, 4, 6):
        f, weight, current, normal, steps = slab.solve(d, alone_density)
        residual, link_error, rate = slab.validate(f, normal, current, d)
        if d == 0:
            alone_density = abs(f[-1])**2
        np.savez_compressed(args.output/f"d{d}.npz", layers=f[-4:], current=current,
                            weight=weight, dx=slab.dx, dz=slab.dz, width=args.width,
                            depth=args.depth, separation=d)
        print(f"d={d} rate_per_mouth={rate:.12g} steps={steps} "
              f"max_equation_residual={residual:.3g} max_link_rate_error={link_error:.3g}",
              flush=True)
        del f


if __name__ == "__main__":
    main()
