"""Does an orientation pulse survive complete carrier replacement at rest?

Equal-travel-time parcels on a single radial lane. A synchronous inward step
replaces the outer parcel and discharges the inner parcel every transit/N.
Exchange acts only on the two-state orientation, leaving material IDs attached
to their carriers. The readout is the linear complex tilt in a frame rotating
at the independently measured common transition frequency. Full tensor-product
evolution is compared with a prediction using only a held pair's two meters.
"""
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from mouth_marked_collar import TRANSIT


def swap_matrix(size, left):
    states = np.arange(2**size)
    flip = ((states >> (size-1-left)) ^ (states >> (size-2-left))) & 1
    partners = states ^ (flip << (size-1-left)) ^ (flip << (size-2-left))
    return np.eye(2**size)[partners]


def pulse(size, weights):
    # Derivative of a product of states |ground> + epsilon*w_i |excited>.
    rho = np.zeros((2**size,2**size), complex)
    excited = 1 << np.arange(size-1,-1,-1)
    rho[excited,0] = weights
    rho[0,excited] = np.conj(weights)
    return rho


def sensors(rho, size):
    rows = np.arange(len(rho))
    values = []
    for site in range(size):
        bit = 1 << (size-1-site)
        low = rows[(rows & bit)==0]
        values.append(np.sum(rho[low+bit,low]))
    return np.array(values)


def master(size, rate, angles):
    # Independent Poisson encounters. Work with the complete tensor-product
    # density matrix, including sectors not excited by the weak pulse.
    from scipy.sparse import eye, kron
    dimension = 2**size
    generator = csr_matrix((dimension**2,dimension**2),dtype=complex)
    identity = eye(dimension,format="csr")
    for left in range(size-1):
        swap = csr_matrix(swap_matrix(size,left))
        for angle in angles:
            unitary = np.cos(angle)*identity-1j*np.sin(angle)*swap
            generator += (rate/len(angles))*(kron(unitary,unitary.conj(),format="csr")-eye(dimension**2,format="csr"))
    return generator


def push_inward(rho):
    # Read and trace out the innermost site. Inject an independent ground site
    # at the outer radius. A state handoff cannot duplicate its departing copy.
    other = len(rho)//2
    remaining = np.trace(rho.reshape(other,2,other,2),axis1=1,axis2=3)
    fresh = np.zeros_like(rho)
    fresh[:other,:other] = remaining
    return fresh


def measure_pair(rate, angles):
    rho = pulse(2, [1.,0.])
    evolution = master(2,rate,angles)
    new_coherence = sensors((evolution@rho.ravel()).reshape(4,4),2)[1]
    return new_coherence


def run(size, rate, angles, name):
    interval = TRANSIT/size
    measured = measure_pair(rate,angles)
    # Predict from the initial held-neighbor transfer rate and known clock only.
    edges = np.diff(np.eye(size),axis=0)
    collisions = expm(-interval*measured*(edges.T@edges))
    inward = np.eye(size,k=-1)
    step = inward@collisions
    spectral_radius = max(abs(np.linalg.eigvals(step)))
    ones = np.ones(size)
    injection = np.eye(size)[0]
    # Time integral of the response to a unit stamp entering just after a shift.
    # For real mixing this is the pulse's mean exit time. For complex mixing
    # it is the integrated rotating-frame amplitude, not a lifetime statistic.
    entry_integral = interval*np.sum(np.linalg.solve(np.eye(size)-step,injection))
    model = ones.astype(complex)
    rho = pulse(size,ones)
    full = interval*master(size,rate,angles)
    budget, max_error = 0j,0.
    for n in range(1,2*size+1):
        between = expm_multiply(full,rho.ravel()).reshape(rho.shape)
        budget += sensors(between,size)[-1]
        rho = push_inward(between)
        direct = sensors(rho,size)
        model = step@model
        max_error = max(max_error,max(abs(direct-model)))
        if n in (size//2,size,2*size):
            survivor = max(0.,1-n/size)
            signal = np.mean(direct)
            print(f"channel={name} steps={n} elapsed_over_transit={n/size:g} old_ID_fraction={survivor:g} "
                  f"tilt_I={signal.real:.10g} tilt_Q={signal.imag:.10g} explained_fraction="
                  f"{((budget+np.sum(direct))/size).real:.10g}")
        assert abs(budget+sum(direct)-size)<2e-10 and abs(np.trace(rho))<2e-11
    print(f"channel={name} held_neighbor_rate=({measured.real:.10g},{measured.imag:.10g}) "
          f"largest_mode_modulus={spectral_radius:.10g} entry_integral_over_transit="
          f"{(entry_integral/TRANSIT).real:.10g} max_full_spin_error={max_error:.3g}")
    assert abs(entry_integral/TRANSIT-1)<1e-10 and spectral_radius<1 and max_error<1e-10
    return step, measured


def main():
    size, rate = 6, .75
    print(f"radial_sites={size} transit={TRANSIT:.10g} collision_rate={rate:g}")
    pure, _ = run(size,0.,[0.],"carriage")
    random, real_rate = run(size,rate,[np.pi/4,-np.pi/4],"mixed_phase")
    coherent, complex_rate = run(size,rate,[np.pi/4],"fixed_phase")
    def remaining(step):
        return np.mean(np.linalg.matrix_power(step,size)@np.ones(size))
    print(f"after_one_transit_I_difference={float((remaining(coherent)-remaining(random)).real):.10g} "
          f"Q_difference={float((remaining(coherent)-remaining(random)).imag):.10g}")
    assert abs(remaining(pure))<1e-12 and remaining(random).real>.01
    assert abs(real_rate.real-complex_rate.real)<1e-13 and abs(remaining(coherent).imag)>.01


if __name__=="__main__":
    main()
