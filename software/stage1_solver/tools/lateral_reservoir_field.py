"""Axisymmetric steady carrier with a lateral bath, one held target and optional companion.

The brane is three-dimensional. x runs along the pair and r is radius in the other
two brane directions. Both normal faces reflect except inside the open mouths.
The cylindrical outer wall and endcaps fix psi=1 at their faces. No lateral
current or inverse-distance profile is prescribed.
"""

import argparse
from pathlib import Path

import numpy as np
from scipy.fft import dct, dst, idct, idst
from scipy.linalg import eigh_tridiagonal


class LateralSlab:
    def __init__(self, extent, cells, normal_cells):
        self.dx = h = 1/cells
        self.dz = 1/normal_cells
        nr = extent*cells
        self.x = (np.arange(2*nr)+0.5-nr)*h
        self.r = r = (np.arange(nr)+0.5)*h
        self.dv = 2*np.pi*r*h*h
        # Radial finite volumes use annular face areas. Symmetrize the
        # operator in the r-weighted measure before diagonalization.
        radial_diag = np.full(nr, 1/h**2)
        radial_diag[-1] += extent/(2*r[-1]*h**2)
        radial_off = -(r[:-1]+h/2)/(2*h**2*np.sqrt(r[:-1]*r[1:]))
        er, self.modes = eigh_tridiagonal(radial_diag, radial_off)
        self.sqrt_r = np.sqrt(r)
        self.sqrt_z = np.ones(normal_cells+1)
        self.sqrt_z[[0,-1]] = np.sqrt(0.5)
        ex = 2*np.sin(np.pi*np.arange(1, 2*nr+1)/(4*nr))**2/h**2
        ez = 2*np.sin(np.pi*np.arange(normal_cells+1)/(2*normal_cells))**2/self.dz**2
        self.eigenvalues = ez[:, None, None]+ex[None, :, None]+er[None, None, :]

    def elliptic(self, rhs, normal, shift):
        rhs[-1] += normal/self.dz
        rhs *= self.sqrt_z[:, None, None]
        spec = dst(dct(rhs, type=1, axis=0, norm="ortho"), type=2, axis=1, norm="ortho")
        spec = ((spec*self.sqrt_r) @ self.modes)/(self.eigenvalues+shift)
        spec = (spec @ self.modes.T)/self.sqrt_r
        return idct(idst(spec, type=2, axis=1, norm="ortho"),
                    type=1, axis=0, norm="ortho")/self.sqrt_z[:, None, None]

    def solve(self, d, reference_density=None):
        centers = [0] if d == 0 else [-d]
        if reference_density is not None:
            centers.append(0)
        weights, shapes = [], []
        for center in centers:
            r2 = (self.x[:, None]-center)**2 + self.r**2
            weight = np.maximum(1-r2, 0)**3
            weights.append(weight)
            if reference_density is None:
                shapes.append(weight)
            else:
                offset = round(center/self.dx)
                if not np.isclose(offset*self.dx, center):
                    raise ValueError("Mouth translations must be grid-aligned")
                shapes.append(weight/np.roll(reference_density, offset, axis=0))
        restoring = sum(weights)

        def boundary(face):
            n = abs(face)**2
            currents = []
            for shape in shapes:
                raw = shape if reference_density is None else shape*n
                currents.append(0.25*raw/np.sum(raw*self.dv))
            normal = (-restoring*(1-1/abs(face)) + 1j*sum(currents)/n)*face
            return currents, normal

        f = np.ones(self.eigenvalues.shape, complex)
        for step in range(1, 501):
            _, normal = boundary(f[-1])
            potential = (abs(f)**8-1)/4
            re = self.elliptic(2*(f.real-1)-potential*f.real, normal.real, 2)
            im = self.elliptic(-potential*f.imag, normal.imag, 0)
            update = 1+re+1j*im
            delta = np.max(abs(update-f))
            if not np.isfinite(delta):
                raise ArithmeticError("Nonfinite carrier update")
            f = update
            if delta < 5e-13:
                currents, normal = boundary(f[-1])
                self.check_field(f, normal, currents)
                return f, np.array(currents), step
        raise ArithmeticError(f"Steady nonlinear solve did not converge: {delta}")

    def check_field(self, f, normal, currents):
        # Apply the real-space Laplacian, independently of the spectral solve.
        # The external half-cell Dirichlet faces have ghosts 2-f.
        axial = np.concatenate((2-f[:, :1], f, 2-f[:, -1:]), axis=1)
        lap = np.diff(axial, n=2, axis=1)/self.dx**2
        radial = np.concatenate((f[:, :, :1], f, 2-f[:, :, -1:]), axis=2)
        radii = np.arange(len(self.r)+1)*self.dx
        lap += np.diff(radii*np.diff(radial, axis=2), axis=2)/(self.r*self.dx**2)
        vertical = np.concatenate((f[1:2], f, (f[-2]+2*self.dz*normal)[None]), axis=0)
        lap += np.diff(vertical, n=2, axis=0)/self.dz**2
        residual = np.max(abs(-lap/2 + (abs(f)**8-1)*f/4))
        # The outer face link samples the bath half a lateral cell away.
        endcap_area = 2*np.pi*self.r*self.dx
        ends = ((f[:, 0].imag+f[:, -1].imag)*endcap_area).sum(axis=1)
        wall = 2*np.pi*radii[-1]*self.dx*f[:, :, -1].imag.sum(axis=1)
        inflow = self.dz*np.sum(self.sqrt_z**2*(ends+wall))*2/self.dx
        rates = np.array([np.sum(j*self.dv) for j in currents])
        if (not np.isfinite([residual,inflow]).all() or residual > 1e-8
                or max(abs(rates-0.25)) > 3e-14 or abs(inflow-rates.sum()) > 2e-8):
            raise ArithmeticError(f"Equation or current mismatch: {residual=} {inflow=} {rates=}")
        print(f"  mouths={len(currents)} residual={residual:.4g} bath_inward={inflow:.12g} "
              f"per_mouth={','.join(format(j,'.12g') for j in rates)}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extent", type=int, required=True)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--normal-cells", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.extent <= 10 or args.cells < 2 or args.normal_cells < 4:
        raise ValueError("Insufficient room or resolution for the mouths and readouts")
    args.output.mkdir(parents=True, exist_ok=False)
    slab = LateralSlab(args.extent, args.cells, args.normal_cells)
    isolated, currents, steps = slab.solve(0)
    print(f"isolated steps={steps}", flush=True)
    reference = abs(isolated[-1])**2
    np.savez_compressed(args.output/"d0.npz", field=isolated, currents=currents,
                        x=slab.x, r=slab.r, dz=slab.dz, d=0)
    for d in (3, 4, 6, 8):
        f, currents, steps = slab.solve(d, reference)
        np.savez_compressed(args.output/f"d{d}.npz", layers=f[-4:], currents=currents,
                            x=slab.x, r=slab.r, dz=slab.dz, d=d)
        print(f"d={d} steps={steps}", flush=True)


if __name__ == "__main__":
    main()
