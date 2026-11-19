"""Standing-wave and rim-load tests on a locally level conversion patch."""
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import legder, leggauss, legval, legvander
from scipy.interpolate import BSpline
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, kron
from scipy.special import j0, jn_zeros


def strip():
    # Constant-width strip, intrinsic distance measured back from the free edge.
    z, k = sp.symbols("z k", positive=True)
    mode = sp.besselj(0, k*z)
    residual = sp.simplify(sp.expand_func((z*mode.diff(z)).diff(z) + k*k*z*mode))
    edge = sp.limit(z*mode.diff(z), z, 0)
    moment = z**2*(mode**2 + sp.besselj(1, k*z)**2)/2
    norm_residual = sp.simplify(sp.expand_func(moment.diff(z) - z*mode**2))
    root = jn_zeros(0, 1)[0]
    print(f"strip ODE={residual} free_edge_traction={edge} norm_identity={norm_residual} "
          f"kL={root:.12f} midpoint_U={j0(root/2):.12f} node_U={j0(root):.3e}", flush=True)
    if (residual, edge, norm_residual) != (0, 0, 0) or abs(j0(root)) > 1e-12:
        raise ArithmeticError("Strip equation, boundary or norm failed")


def transition(x):
    return x**3*(10-15*x+6*x**2)


class Collar:
    """Axisymmetric scalar on a 3D graph. Cubic radial splines and angular Legendre modes."""

    def __init__(self, cells, lmax):
        edges = np.linspace(1, 3, cells+1)
        knots = np.r_[np.repeat(1., 3), edges, np.repeat(3., 3)]
        basis = BSpline(knots, np.eye(len(knots)-4), 3)
        q, wr = leggauss(7)
        radius = ((edges[1:, None]+edges[:-1, None])/2 + q/cells).ravel()
        mu, wm = leggauss(2*lmax+12)
        polys = np.eye(lmax+1)
        angular = legvander(mu, lmax)
        angular_der = np.column_stack([legval(mu, legder(p)) for p in polys])
        # Only the endpoint value is fixed at the inner node. The final spline
        # is retained at the free outer edge. Regularity at the angular axes
        # is implicit in the expansion in cos(theta).
        rad = csr_matrix(basis(radius)[:, 1:])
        self.value = kron(rad, angular, format="csr")
        self.dr = kron(csr_matrix(basis(radius, nu=1)[:, 1:]), angular, format="csr")
        self.dm = kron(rad, angular_der, format="csr")
        self.r, self.mu = [a.ravel() for a in np.meshgrid(radius, mu, indexing="ij")]
        self.rim_dr = np.kron(basis(1., nu=1)[None, 1:], angular)
        self.rim_weights = 2*np.pi*wm
        self.epsilon = 0.2
        self.chi = transition(np.clip(3-self.r, 0, 1))
        self.weight = 2*np.pi*np.kron(np.tile(wr/cells, cells), wm)*self.r**2*self.chi
        self.basis = basis
        self.lmax = lmax

    def solve(self, beta):
        r, mu, e = self.r, self.mu, self.epsilon
        x = np.clip(r-1, 0, 1)
        gate = transition(x)
        gate_dr = 30*x**2*(1-x)**2
        hr = e*(r-beta*mu*(gate+r*gate_dr))
        hm = -e*beta*r*gate
        angular_metric = (1-mu**2)/r**2
        volume = np.sqrt(1+hr**2+angular_metric*hm**2)
        scale = self.weight/volume

        def gram(a, coeff, b):
            return a.T @ b.multiply(coeff[:, None])

        mixed = gram(self.dr, -scale*hr*angular_metric*hm, self.dm)
        stiff = (gram(self.dr, scale*(1+angular_metric*hm**2), self.dr)
                 + mixed + mixed.T + gram(self.dm, scale*angular_metric*(1+hr**2), self.dm))
        mass = gram(self.value, self.weight*volume, self.value)
        vals, vecs = eigh(stiff.toarray(), mass.toarray(), subset_by_index=(0, 1))
        u = vecs[:, 0]
        omega = np.sqrt(vals[0])
        defect = np.linalg.norm(stiff@u - vals[0]*(mass@u))/np.linalg.norm(stiff@u)
        # For unit mass norm, energy=omega^2/2. Mean outward pressure on
        # the reflecting inner boundary is |normal derivative|^2/4.
        rim_per_energy = np.dot(self.rim_weights, (self.rim_dr@u)**2)/((1+e*e)*2*vals[0])
        midpoint = abs(np.kron(self.basis(2.5)[1:], np.ones(self.lmax+1)) @ u)
        if defect > 1e-8 or vals[0] <= 0 or not np.isfinite(rim_per_energy) or rim_per_energy <= 0:
            raise ArithmeticError(f"Mode or rim load unresolved: {vals[0]}, {defect}, {rim_per_energy}")
        return omega, rim_per_energy, midpoint, defect, np.sqrt(vals[1])


if __name__ == "__main__":
    strip()
    for cells, angular in [(12, 6), (24, 10), (36, 14)]:
        collar = Collar(cells, angular)
        frequencies = {}
        for tilt in [0, 2, 2.49, 2.5, 2.51, 3, 4]:
            omega, load, value, residual, second = collar.solve(tilt)
            frequencies[tilt] = omega
            print(f"cells={cells} angular={angular} tilt={tilt:g} omega={omega:.10f} "
                  f"rim_load_per_energy={load:.10f} abs_U_at_r2.5_north={value:.10f} "
                  f"relative_residual={residual:.3e} next_axisymmetric_omega={second:.8f}", flush=True)
        force = -(frequencies[2.51]-frequencies[2.49])/0.02
        print(f"cells={cells} generalized_tilt_force_per_action_at_2.5={force:+.10f}", flush=True)
