"""One proposed neighbor Hamiltonian tested in a held lane and a static mouth.

For the finite-annulus certificate, retain the earlier product-direction state
space on the two material faces. Add six copies of the cut angular mesh on
successive sections, sharing directions along each radial lane in the trial.
The dynamical test isolates a radial lane, with inward carrier replacement.
"""
import numpy as np
import sympy as sp
from scipy.linalg import expm

from mouth_collar_orientation_pulse import pulse, sensors, push_inward, swap_matrix
from mouth_torsion_exterior import spring_energy


LAYERS, EXCHANGE = 6, .002  # energies in units of the angular spring K
HBAR = 1.  # K times the model time unit


def radial_link():
    eps, imaginary = sp.symbols("epsilon b", real=True)
    # For rays (1,epsilon*z_i) with real normalization, the overlap has
    # imaginary part O(epsilon^2). Expanding its positive real denominator
    # around 1 does not change the leading coefficient of |v|-Re(v).
    realpart = 1+sp.Symbol("a",real=True)*eps**2
    imagpart = imaginary*eps**2
    old = sp.sqrt(realpart**2+imagpart**2)-realpart
    leading = sp.simplify(sp.series(old,eps,0,5).removeO())
    swap = swap_matrix(2,0)
    energies = np.linalg.eigvalsh(EXCHANGE*(np.eye(4)-swap))
    print(f"radial_old_link_first_term={leading} exchange_symmetric_energy={energies[0]:.10g} "
          f"antisymmetric_energy={energies[-1]:.10g}")
    assert sp.simplify(leading-imaginary**2*eps**4/2)==0 and energies[0]>=0
    return (energies[-1]-energies[0])/2


def trial_edges(R):
    # Enumerate each physical angular edge once. Use the same radial directions
    # on scaled copies of this mesh at the other five radial layers.
    vals = np.arange(-R,R+1)
    grid = np.stack(np.meshgrid(vals,vals,vals,indexing="ij"),axis=-1)
    on = np.max(np.abs(grid),axis=-1)==R
    count = 0
    costs = []
    for side in (1,-1):
        cost = 0.
        for axis in range(3):
            a, b = [slice(None)]*3, [slice(None)]*3
            a[axis], b[axis] = slice(None,-1), slice(1,None)
            use = on[tuple(a)] & on[tuple(b)]
            x = side*grid[tuple(a)][use]
            y = side*grid[tuple(b)][use]
            dot = np.sum(x*y,axis=1)/(np.linalg.norm(x,axis=1)*np.linalg.norm(y,axis=1))
            overlap2_one_face = (1+dot)/2
            cost += np.sum(1-overlap2_one_face**2)
            if side==1:
                count += len(x)
        costs.append(cost)
    assert count==48*R**2 and max(costs)<=24 and costs[0]==costs[1]
    return count,costs[0]


def certificate(R):
    # Earlier lower bounds for the old energy only. New exchange is positive
    # on every competitor; an upper bound is needed only on the radial trial.
    g, pin, angle = 4., 400., .1
    h = np.cos(1/(2*R))
    error = 4*np.sin(angle/4)
    slip = (np.pi-1/(2*R**2))/4-np.arcsin(error/h)
    wrong_rim = 2*g*(h-error)*(1-np.cos(slip))
    floors = [wrong_rim, pin*(1-np.cos(angle)),12*spring_energy(np.pi/3),2*spring_energy(np.pi)]
    upper = (72+48/R+12/R**2)/(16*(R-1))
    count,geometry = trial_edges(R)
    bound = upper+24*LAYERS*EXCHANGE
    print(f"static R={R} edges={count} measured_exchange_cost_six_layers={EXCHANGE*LAYERS*geometry:.10g} "
          f"old_upper={upper:.12g} combined_upper={bound:.12g} "
          f"wrong_floor={min(floors):.12g} margin={min(floors)-bound:.10g} "
          f"sufficient_J_over_K_below={(min(floors)-upper)/(24*LAYERS):.10g}")
    assert min(floors)>bound and geometry*LAYERS*EXCHANGE<=bound-upper+1e-14


def held_lane(J, transit):
    size = LAYERS
    dt = transit/size
    # Isolate the resonant transverse transition of identical complete registers.
    # The full swap leaves this two-level local subspace invariant.
    H = sum(J*(np.eye(2**size)-swap_matrix(size,i)) for i in range(size-1))
    propagation = expm(-1j*dt*H/HBAR)
    rho = pulse(size,np.ones(size))
    edges = np.diff(np.eye(size),axis=0)
    analytic = np.eye(size,k=-1)@expm(-1j*dt*J/HBAR*(edges.T@edges))
    predicted = np.ones(size,dtype=complex)
    for n in range(size):
        rho = push_inward(propagation@rho@propagation.conj().T)
        predicted = analytic@predicted
    direct = sensors(rho,size)
    maximum = max(abs(direct-predicted))
    signal = np.mean(direct)
    print(f"always_on J={J:.8g} tau={transit:.10g} J_tau_over_hbar={J*transit/HBAR:.10g} "
          f"one_transit_I={signal.real:.10g} Q={signal.imag:.10g} amplitude={abs(signal):.10g} "
          f"max_full_spin_error={maximum:.3g}")
    assert maximum<1e-12
    return signal


def main():
    measured_J = radial_link()
    for R in (5,10):
        certificate(R)
    # Separate static (cubic cut) and dynamic (spherical feeding lane) tests.
    # Radius-five example; Q is a supplied volume/time for this one benchmark.
    Q = .25
    transit = 4*np.pi*(6**3-5**3)/(3*Q)
    assert held_lane(0.,transit)==0
    assert abs(held_lane(measured_J,transit))>1e-2
    held_lane(measured_J,transit/2)


if __name__=="__main__":
    main()
