"""Independent pressure-driven stationary front and flat-patch relaxation checks."""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import simpson, solve_bvp
from scipy.linalg import eigh_tridiagonal
from scipy.special import expit


@dataclass(frozen=True)
class Material:
    a: float = 2.0
    kappa: float = 1.0
    mobility: float = 0.4
    pin: float = 1.5
    density: float = 2.0
    area: float = 3.0

    @property
    def delta(self):
        return (self.kappa/(2*self.a))**0.5

    def reference(self, z):
        chi = expit(-z/self.delta)
        c = chi*(1-chi)
        return chi, -c/self.delta, 4*c, 4*c*(2*chi-1)/self.delta


def branch(material, head, extent):
    """Apply a pressure drop using endpoint data, solving for the crossing speed."""
    m = material
    z = np.linspace(-extent*m.delta, extent*m.delta, 501)
    chi, deriv, _, _ = m.reference(z)
    pressure = head*m.density

    def equations(z, state, parameter):
        order, slope, _ = state
        reference, _, weight, _ = m.reference(z)
        speed, = parameter
        dw = 2*m.a*order*(1-order)*(1-2*order)
        pinning = m.pin*weight*(order-reference)
        return np.array([slope, (dw+pinning+speed*slope/m.mobility)/m.kappa,
                         speed*slope**2/m.mobility])

    def boundary(left, right, _parameter):
        return np.array([left[0]-chi[0], right[0]-chi[-1], left[2], right[2]-pressure])

    seed = np.array([chi, deriv, pressure*(1-chi)])
    sol = solve_bvp(equations, boundary, z, seed, p=[pressure],
                    tol=2e-9, bc_tol=2e-10, max_nodes=6000)
    if not sol.success:
        raise RuntimeError(sol.message)
    z = np.linspace(z[0], z[-1], 2001)
    order, slope, _ = sol.sol(z)
    ref, ref_prime, weight, weight_prime = m.reference(z)
    speed = sol.p[0]
    I = simpson(slope**2, x=z)
    # Generalized force acting ON the fixed anchor: -partial_X F per area.
    force = m.pin*simpson(weight_prime*(order-ref)**2/2-weight*(order-ref)*ref_prime, x=z)
    dissipation = m.area*speed**2*I/m.mobility
    rate = m.area*m.density*speed
    error = dissipation-rate*head
    print(f"extent_delta={extent} head={head:.6g} rate={rate:.12g} rate/head={rate/head:.12g} "
          f"anchor_dp={force:.12g} imposed_dp={pressure:.12g} power_error={error:.3g} "
          f"half_level={np.interp(0.5, order[::-1], z[::-1]):.12g}")
    assert abs(force-pressure)<2e-8 and abs(error)<1e-9
    assert np.all(np.diff(order) <= 1e-9) and pressure <= m.pin/2
    return rate


def relaxation(material):
    """Measure the rate of the same material-depth eigenmode at two lateral k."""
    m = material
    z = np.linspace(-20*m.delta, 20*m.delta, 1601)
    dz = z[1]-z[0]
    ref, _, w, _ = m.reference(z[1:-1])
    local = 2*m.a*(1-6*ref*(1-ref))+m.pin*w
    diag = 2*m.kappa/dz**2 + local
    offdiag = np.full(len(diag)-1, -m.kappa/dz**2)
    wave_numbers = np.array([.4, .9])
    rates = []
    for wave in wave_numbers:
        eigen = eigh_tridiagonal(diag+m.kappa*wave**2, offdiag,
                                 select="i", select_range=(0, 0), eigvals_only=True)
        rates.append(m.mobility*eigen[0])
    diffusivity = (rates[1]-rates[0])/(wave_numbers[1]**2-wave_numbers[0]**2)
    power = (-1+(25-8*m.pin/m.a)**.5)/2
    isolated_eigenvalue = m.a*(4-power**2)/2
    observed = rates[0]/m.mobility - m.kappa*wave_numbers[0]**2
    sigma = m.kappa/(6*m.delta)
    conductance = m.area*m.density**2*diffusivity/sigma
    print(f"wave_rates={rates[0]:.12g},{rates[1]:.12g} diffusivity={diffusivity:.12g} "
          f"normal_ground_eigen={observed:.12g} analytic={isolated_eigenvalue:.12g} "
          f"zero_head_conductance={conductance:.12g}")
    assert abs(diffusivity-m.mobility*m.kappa)<1e-9
    assert abs(isolated_eigenvalue-observed)<4e-5
    return conductance


if __name__ == "__main__":
    material = Material()
    predicted = relaxation(material)
    for head in (.002, .001):
        short = branch(material, head, 16)
        long = branch(material, head, 22)
        assert abs(short-long)<1e-8 and abs(long/head-predicted)<.005
