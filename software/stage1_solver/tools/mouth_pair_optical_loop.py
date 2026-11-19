"""Read the same physical constituent on a closed square with true circulation."""
import numpy as np


def test(hoppings, binding):
    size = len(hoppings)
    eye = np.eye(size)
    motion = np.zeros((size,size))
    phases = []
    for site, hop in enumerate(hoppings):
        neighbor = (site+1)%size
        motion[site,neighbor] = motion[neighbor,site] = -hop
        current = np.zeros((size,size),complex)
        current[neighbor,site] = -1j*hop
        current[site,neighbor] = 1j*hop
        phases.append(current)
    # Uniformly shared directed phase F/4 winds once around the square.
    loop = sum(phases)/size
    def odd(op):
        return np.kron(op,eye)-np.kron(eye,op)
    def common(op):
        return np.kron(op,eye)+np.kron(eye,op)
    base = common(motion) + np.diag((binding*(1-eye)).ravel())
    energies, states = np.linalg.eigh(base)
    ground = states[:,0]
    gaps = energies-energies[0]
    low = (energies<=1e-9) & (gaps>1e-9)
    # Also redistribute the same time-varying flux by a site rephasing.
    # A physically equivalent experiment includes its on-site scalar term.
    charges = np.diag([.2,-.1,.35,-.45])
    gen = odd(charges)
    difference = 1j*(gen@base-base@gen)
    raw = states.conj().T @ (odd(loop)@ground)
    moved = states.conj().T @ ((odd(loop)+difference)@ground)
    fixed = moved + 1j*gaps*(states.conj().T@(gen@ground))
    local_real = np.zeros_like(motion)
    local_real[0,1] = local_real[1,0] = -hoppings[0]
    common_lines = states.conj().T @ (common(local_real)@ground)
    single = np.linalg.eigvalsh(motion)
    odd_edge = binding+single[:2].sum()
    # Antisymmetric products have zero double occupancy, so this is the
    # exact odd-sector spectrum even while the two particles attract.
    antisymmetric = []
    for i in range(size):
        for j in range(i+1,size):
            a = np.zeros((size,size))
            a[i,j], a[j,i] = 1/np.sqrt(2), -1/np.sqrt(2)
            antisymmetric.append(a.ravel())
    project = np.array(antisymmetric)
    odd_measured = np.linalg.eigvalsh(project@base@project.T)
    odd_predicted = sorted(binding+single[i]+single[j] for i in range(size) for j in range(i+1,size))
    bright = (abs(raw)>1e-7) & (gaps>1e-9)
    print(f"hoppings={hoppings} Delta={binding:g} low_states_including_ground={int(low.sum())+1} "
          f"ground={energies[0]:.9f} odd_gap_floor={odd_edge-energies[0]:.9f} "
          f"lowest_loop_active_gap={gaps[bright].min():.9f}")
    print(f"loop_strength_into_low={sum(abs(raw[low])**2):.3g} "
          f"loop_strength_into_high={sum(abs(raw[bright])**2):.9f} "
          f"common_local_strain_strength_into_low={sum(abs(common_lines[low])**2):.9f} "
          f"frame_residual={max(abs(fixed-raw)):.3g}")
    assert np.allclose(odd_measured,odd_predicted)
    assert low.sum()==size-1 and abs(raw[low]).max()<1e-11
    assert max(abs(fixed-raw))<1e-12 and sum(abs(common_lines[low])**2)>0.001


if __name__ == "__main__":
    for bonds in ([1.,1.,1.,1.],[1.,.85,1.12,.93]):
        test(bonds,16)
