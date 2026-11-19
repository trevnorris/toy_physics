"""Stationary neutral carrier through a transparent aperture from a slab into a bath."""

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from stage1_solver.config import TensorGridSpec
from stage1_solver.open_sleeve import OpenSleeve


class JoinedMouth(OpenSleeve):
    """Reuse the nonlinear carrier solve with two baths and an open internal face."""

    def __post_init__(self):
        n = self.mouth_cells
        self.h = h = 1/n
        nr, nw = self.grid.nr, self.grid.nw
        ids = np.arange(nr*nw).reshape(nr, nw)
        edges = np.arange(nr+1)*h
        self.shell = 4*np.pi/3*np.diff(edges**3)
        areas = 4*np.pi*edges**2
        across = np.ones((nr, nw-1), bool)
        across[n:, n-1] = False  # insulating screen outside the radius-one opening
        zlinks = np.broadcast_to(self.shell[:, None]/h, across.shape)
        self.left = np.r_[ids[:-1].ravel(), ids[:, :-1][across]]
        self.right = np.r_[ids[1:].ravel(), ids[:, 1:][across]]
        self.conductance = np.r_[np.repeat(areas[1:-1], nw), zlinks[across]]
        self.scale = areas[-1]  # one common row scaling, including the prescribed rate
        self.volume = np.repeat(self.shell*h, nw)/self.scale
        inlet, outlet = np.zeros_like(ids, float), np.zeros_like(ids, float)
        inlet[-1, :n] = 2*areas[-1]    # side of feeding slab
        outlet[-1, n:] = 2*areas[-1]   # side and top of receiving bath
        outlet[:, -1] += 2*self.shell/h
        self.inlet, self.outlet = inlet.ravel()/self.scale, outlet.ravel()/self.scale
        off = coo_matrix((-self.conductance/self.scale, (self.left, self.right)),
                         shape=(ids.size, ids.size)).tocsr()
        off = off + off.T
        diag = -np.asarray(off.sum(axis=1)).ravel() + self.inlet + self.outlet
        self.kinetic = 0.5
        self.operator = self.kinetic*(off + diags(diag))


def normal_export(f, shell, h, n, radius):
    """Flux out of the fluid on the artificial top and radial side of a cylinder."""
    i = round(radius/h)
    top = n+i
    density = abs(f)**2
    dw = (f[:i, top]-f[:i, top-1])/h
    nww = (density[:i, top-2]-density[:i, top-1]-density[:i, top]
           +density[:i, top+1])/(2*h*h)
    # Interpolate density itself so that the pressure does not pick up a
    # spurious density depression from averaging two different phases.
    ntop = (density[:i, top-1]+density[:i, top])/2
    top_stress = abs(dw)**2 + np.expm1(5*np.log(ntop))/5 - nww/4
    top_flux = np.dot(shell[:i], top_stress)

    dr = (f[i]-f[i-1])/h
    mid = (f[i]+f[i-1])/2
    dn = (density[i]-density[i-1])/h
    # Mixed stress on the cylindrical side. Reflect each domain separately
    # at the solid boundaries. Integrate the mixed density derivative by endpoints.
    def neumann_face(first, second):
        return (9*first-second)/8
    bottom_density_slope = neumann_face(dn[0], dn[1])
    below_screen = neumann_face(dn[n-1], dn[n-2])
    above_screen = neumann_face(dn[n], dn[n+1])
    cap_density_slope = (dn[top-1]+dn[top])/2
    dmid = np.empty(top, complex)
    for a, b in ((0, n), (n, top)):
        after = mid[b] if b == top else mid[b-1]
        padded = np.r_[mid[a], mid[a:b], after]
        dmid[a:b] = (padded[2:]-padded[:-2])/(2*h)
    mixed = h*np.vdot(dr[:top], dmid).real
    mixed -= (below_screen-bottom_density_slope+cap_density_slope-above_screen)/4
    side_flux = 4*np.pi*radius**2*mixed
    if not np.all(np.isfinite([side_flux, top_flux])):
        raise ArithmeticError("Unresolved exterior flux")
    return top_flux, side_flux


def measure(s, rate, cuts):
    n, h = s.mouth_cells, s.h
    flat, phase, residual = s.solve(rate/s.scale)
    f = flat.reshape(s.grid.nr, s.grid.nw)
    radial = (np.arange(s.grid.nr)+0.5)*h
    normal = (f[:n, n-1].conj()*f[:n, n]).imag/h
    inward = s.shell[:n]*normal
    bath_in = s.scale*np.dot(s.inlet, (np.exp(-1j*phase)*flat).imag)
    bath_out = -s.scale*np.dot(s.outlet, flat.imag)
    # Independent one-sided extrapolations inside the opening avoid the singular rim.
    below = (3*f[:n, n-1]-f[:n, n-2])/2
    above = (3*f[:n, n]-f[:n, n+1])/2
    def contrast(trace):
        angle = np.angle(trace)
        axis = (9*angle[0]-angle[1])/8
        return axis-np.interp(0.5, radial[:n], angle)
    print(f"cells_per_radius={n} extent={s.grid.r_max:g} rate={rate:g} "
          f"feed={bath_in:.12g} mouth={inward.sum():.12g} out={bath_out:.12g} "
          f"far_feed_phase={phase:+.12g} scaled_residual={residual:.3g}", flush=True)
    print(f"  inner_half_rate={inward[:n//2].sum():.12g} min_local_current={normal.min():.12g} "
          f"feed_contrast={contrast(below):+.12g} receiver_contrast={contrast(above):+.12g}", flush=True)
    if not np.allclose([bath_in, inward.sum(), bath_out], rate, rtol=0, atol=1e-8):
        raise ArithmeticError("Number current failed a separate cut")
    for radius in cuts:
        top, side = normal_export(f, s.shell, h, n, radius)
        print(f"  radius={radius} top_export={top:+.12g} side_export={side:+.12g} "
              f"fluid_on_slab={-(top+side):+.12g}", flush=True)
    # Opposite pressure loads on the two sides of the screen, on an annulus
    # clear of the singular lip. At a reflecting wall the normal phase slope vanishes.
    annulus = (radial >= 1.5) & (radial < 3)
    def wall_load(region):
        density = abs(region[annulus, :3])**2
        # Even extension at a cell face. A quartic through 1/2,3/2,5/2
        # reproduces its value and normal curvature.
        surface = density @ np.array([150, -25, 3])/128
        curvature = density @ np.array([-34, 39, -5])/(24*h*h)
        traction = np.expm1(5*np.log(surface))/5 - curvature/4
        return np.dot(s.shell[annulus], traction)
    feed = wall_load(f[:, n-1::-1])
    recv = -wall_load(f[:, n:])
    print(f"  screen_annulus_1.5_to_3 feed={feed:+.12g} receiver={recv:+.12g} "
          f"net={feed+recv:+.12g}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cells", type=int)
    parser.add_argument("extent", type=int)
    parser.add_argument("rates", type=float, nargs="+")
    args = parser.parse_args()
    n, extent = args.cells, args.extent
    if n < 4 or n % 2 or extent < 6:
        parser.error("Use an even resolution >=4 and extent >=6")
    mesh = TensorGridSpec(r_max=extent, nr=extent*n, w_min=-1,
                          w_max=extent, nw=(extent+1)*n)
    solver = JoinedMouth(mesh, mouth_cells=n, sleeve_cells=n, eos_k=0.2)
    cuts = [r for r in (2, 4, 8) if r+1 < extent]
    for throughput in args.rates:
        measure(solver, throughput, cuts)
