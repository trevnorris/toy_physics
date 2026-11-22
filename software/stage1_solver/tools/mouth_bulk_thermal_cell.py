"""Four-periodic skew-cell Brownian dynamics using the switched WITHIN pair energy."""
import argparse
import ctypes as ct
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from time import perf_counter

import numpy as np

from neutral_patch_sites import TEMPERATURE, WITHIN, emit


def library(work: Path):
    """Build the native kernel in work and return its typed ctypes interface."""
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_suffix('.cpp')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    shared = work/f'thermal_cell_{digest}.so'
    if not shared.exists():
        # Atomic publication also permits different jobs to use the same build directory.
        with tempfile.TemporaryDirectory(dir=work) as build:
            compiled = Path(build)/'thermal.so'
            subprocess.run(['g++', '-std=c++17', '-O3', '-march=native', '-fPIC', '-shared',
                            str(source), '-o', str(compiled)], check=True)
            os.replace(compiled, shared)
    lib = ct.CDLL(str(shared))
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    lib.thermal_make.argtypes = [ct.c_int, array, array, ct.c_double, array]
    lib.thermal_make.restype = ct.c_void_p
    lib.thermal_free.argtypes = [ct.c_void_p]
    lib.thermal_free.restype = None
    lib.thermal_read.argtypes = [ct.c_void_p, array, array, array]
    lib.thermal_read.restype = ct.c_double
    lib.thermal_advance.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, array, array]
    lib.thermal_advance.restype = ct.c_int
    return lib


class ThermalCell:
    """Row-vector lattice. The Cartesian strain derivative returned by read is extensive."""
    def __init__(self, lib, lattice, n, cutoff):
        a = np.array(lattice, dtype=np.float64, order='C', copy=True)
        if a.shape != (4, 4) or not np.all(np.isfinite(a)):
            raise ValueError('Expected a finite 4 by 4 lattice')
        if not (0 < n < np.iinfo(np.int32).max//4 and int(n) == n):
            raise ValueError('Expected a positive site count')
        if not (np.isfinite(cutoff) and cutoff > .25):
            raise ValueError('Expected a cutoff larger than the switching width')
        inverse = np.ascontiguousarray(np.linalg.inv(a))
        self.lib, self.lattice, self.n, self.rc = lib, a, int(n), float(cutoff)
        self.lattice.setflags(write=False)
        self.volume = abs(np.linalg.det(a))
        self.handle = lib.thermal_make(self.n, a, inverse, self.rc, np.array(WITHIN))
        if not self.handle:
            raise MemoryError('Could not initialize the thermal cell')

    def __del__(self):
        if getattr(self, 'handle', None):
            self.lib.thermal_free(self.handle)

    def read(self, pos):
        """Return potential energy, Cartesian site gradients, and dE/d(strain) at fixed fractions."""
        x = np.ascontiguousarray(pos, dtype=np.float64)
        if x.shape != (self.n, 4) or not np.isfinite(x).all():
            raise ValueError('Expected finite Cartesian positions with shape (n,4)')
        grad, virial = np.empty_like(x), np.empty((4, 4))
        energy = self.lib.thermal_read(self.handle, x, grad, virial)
        if not np.isfinite(energy):
            raise ArithmeticError('Could not evaluate finite periodic site forces')
        return energy, grad, virial

    def advance(self, pos, noise, dt):
        """Update Cartesian pos in place using pre-scaled increments of shape (steps,n,4)."""
        if (not isinstance(pos, np.ndarray) or pos.dtype != np.float64 or
            pos.shape != (self.n, 4) or not pos.flags.c_contiguous or
            not pos.flags.writeable or not np.isfinite(pos).all()):
            raise ValueError('Expected writable, finite, contiguous float64 positions with shape (n,4)')
        noise = np.ascontiguousarray(noise, dtype=np.float64)
        if noise.ndim != 3 or noise.shape[1:] != (self.n, 4) or len(noise) > np.iinfo(np.int32).max:
            raise ValueError('Expected noise of shape (steps,n,4)')
        if not (np.isfinite(dt) and dt > 0):
            raise ValueError('Expected a positive finite time step')
        failure = self.lib.thermal_advance(self.handle, len(noise), dt, noise, pos)
        if failure:
            raise ArithmeticError(f'Thermal integrator failed at batch step {failure}')


def compare(lib, fractional, lattice, rc, name):
    """Reference force and strain from the independent all-image Python implementation."""
    from mouth_bulk_relaxed_shape import PeriodicMaterial

    n = len(fractional)
    cart = np.ascontiguousarray(fractional@lattice)
    native = ThermalCell(lib, lattice, n, rc)
    oracle = PeriodicMaterial(n, rc)
    e, g, w = native.read(cart)
    ref_e, ref_g, ref_w = oracle.read(fractional, lattice)
    ref_virial = lattice.T@ref_w
    # A nonsymmetric test deformation includes rotation and mixes Cartesian axes.
    transform = np.array([[.4, .25, 0., -.1], [-.34, -.5, .29, .2],
                          [.07, .16, .13, .31], [-.5, -.16, -.27, .33]])
    # Simultaneous Cartesian and cell perturbations check the deformation convention.
    travel = np.random.default_rng(20261121).normal(size=cart.shape)
    dx = 1.e-5
    samples = []
    for sign in (1., -1.):
        a = lattice@(np.eye(4)+sign*dx*transform)
        y = fractional+sign*dx*travel@np.linalg.inv(a)
        samples.append(oracle.read(y, a)[0])
    fd = (samples[0]-samples[1])/(2*dx)
    strain_and_force = np.sum(w*transform)+np.sum(g*travel)
    stats = dict(test='thermal_kernel_comparison', state=name, n=n, energy=e,
                 reference_energy=ref_e, max_reference_force=float(np.max(abs(ref_g))),
                 energy_error=abs(ref_e-e), max_gradient_error=float(np.max(abs(ref_g-g))),
                 max_virial_error=float(np.max(abs(ref_virial-w))),
                 finite_difference=fd, work=float(strain_and_force),
                 finite_difference_error=abs(float(strain_and_force)-fd))
    emit(**stats)
    np.testing.assert_allclose(e, ref_e, rtol=5e-13, atol=2e-11)
    np.testing.assert_allclose(g, ref_g, rtol=2e-11, atol=2e-9)
    np.testing.assert_allclose(w, ref_virial, rtol=2e-11, atol=1e-8)
    np.testing.assert_allclose(strain_and_force, fd, rtol=5e-7, atol=3e-4)
    return native


def check(cold: Path, work: Path):
    from mouth_bulk_relaxed_shape import PeriodicMaterial

    outcome = json.loads((cold.parent.parent/'outcome.json').read_text())
    if not outcome['limitsVerified'] or outcome['exitCode'] or outcome['childOutcome']['exitCode']:
        raise ValueError(f'Unverified input at {cold}')
    with np.load(cold, allow_pickle=False) as data:
        y, a, rc = data['y'], data['cell'], float(data['rc'])
    lib = library(work)
    cell = compare(lib, y, a, rc, 'free_shape_cold')
    x = np.ascontiguousarray(y@a)
    dt = .0002
    rng = np.random.default_rng(48021)
    noise = rng.normal(size=(1024, len(y), 4))*np.sqrt(2*TEMPERATURE*dt)
    ref_g = PeriodicMaterial(len(y),rc).read(y,a)[1]
    one_step = x+noise[0]-dt*ref_g
    cell.advance(x, noise[:1], dt)
    step_error = np.max(abs(x-one_step))
    begin = perf_counter()
    cell.advance(x, noise[1:], dt)
    elapsed = perf_counter()-begin
    emit(test='thermal_kernel_integration', n=len(y), dt=dt, force_steps=1023,
         seconds=elapsed, steps_per_second=1023/elapsed, first_step_max_error=float(step_error))
    np.testing.assert_allclose(step_error, 0., atol=5e-13)
    compare(lib, x@np.linalg.inv(a), a, rc, 'free_shape_warm')
    # Exercise the retained Verlet list across arbitrary site-by-site windings.
    old = cell.read(x)
    rewound = cell.read(x+rng.integers(-3,4,size=x.shape)@a)
    errors = [float(np.max(abs(np.asarray(p)-q))) for p,q in zip(old,rewound)]
    emit(test='thermal_kernel_winding', errors_energy_gradient_virial=errors)
    np.testing.assert_allclose(errors, 0., atol=1e-8)
    # Sample 87's rows exceed the cutoff. A small oblique one-site cell is an
    # independent test of the self images' finite energy, strain and zero force.
    small = np.array([[1.5,.26,0.,0.],[.1,1.7,.12,0.],[.15,0.,1.9,.12],[0.,0.,.33,1.8]])
    compare(lib, np.array([[.21,.17,.42,.55]]), small, rc, 'one_site_self_images')
    emit(test='thermal_kernel_guard', passed=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cold', required=True, type=Path)
    parser.add_argument('--work', required=True, type=Path)
    args = parser.parse_args()
    check(args.cold, args.work)
