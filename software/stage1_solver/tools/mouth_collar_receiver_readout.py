"""Read the prepared material, the collective receivers, and one receiver at a time.

A zero-temperature interaction lowers one transition. The joint state contains
vacuum and at most one emitted quantum over all receiver slots. This basis
represents the pure state exactly for the prescribed fixed pin on each lane.
"""
import numpy as np
from scipy.spatial.transform import Rotation

from mouth_collar_defect import face_angle
from mouth_collar_order_preparation import COUPLING,DURATION,PAULI,full_block,pin
from mouth_material_links import section


def exact_state(axis, slots):
    ground, excited = pin(axis)
    incoming = np.array([1.,0.])
    survivors, lowering = np.cos(np.sqrt(2)*COUPLING*DURATION), -1j*np.sin(np.sqrt(2)*COUPLING*DURATION)/np.sqrt(2)
    joint = np.empty((2,1+3*slots),complex)
    joint[:,0] = (ground+survivors**slots*excited)@incoming
    emissions = np.column_stack([lowering*ground@sigma@excited@incoming for sigma in PAULI])
    for j in range(slots):
        joint[:,1+3*j:4+3*j] = survivors**j*emissions
    return joint


def dense_check():
    axis = np.array([.8,0.,-.6])
    block=full_block(axis,DURATION).reshape(2,4,2,4)
    dense = np.array([1.,0.],complex)
    for j in range(3):
        dense = np.tensordot(block[:,:,:,0],dense,axes=(2,0))
    dense = dense.transpose(0,3,2,1)  # spin, earliest slot, second, latest slot
    joint = exact_state(axis,3)
    sparse = np.zeros_like(dense)
    sparse[:,0,0,0] = joint[:,0]
    for j in range(3):
        index = [slice(None),0,0,0]
        index[1+j] = slice(1,4)
        sparse[tuple(index)] = joint[:,1+3*j:4+3*j]
    err = np.linalg.norm(dense-sparse)
    print(f"three_events_full_tensor_error={err:.3g}")
    assert err<1e-12


def degree(rays, triangles):
    return np.sum(face_angle(np.asarray(rays)[triangles]))/(2*np.pi)


def scan(level,sign,slots):
    normals,tri = section(level)
    normals = Rotation.from_rotvec([.17,-.24,.07]).apply(normals)
    full_rays,spin_rays,receiver_rays = [],[],[]
    individual = [[] for _ in range(slots)]
    worst_vacuum=1.
    for normal in normals:
        joint = exact_state(sign*normal,slots)
        # Joint row/column indices refer to the spin and receiver ket bases.
        spin = joint@joint.conj().T
        probabilities,vectors = np.linalg.eigh(spin)
        preferred=vectors[:,-1]
        partner=preferred.conj()@joint/np.sqrt(probabilities[-1])
        spin_rays.append(preferred)
        receiver_rays.append(partner)
        full_rays.append(joint.ravel())
        for j in range(slots):
            group = joint[:,1+3*j:4+3*j]
            state = np.column_stack((joint[:,0],group))
            marginal = state.T@state.conj()
            marginal[0,0] += 1-np.trace(marginal).real  # emissions in all other slots
            worst_vacuum=min(worst_vacuum,marginal[0,0].real)
            individual[j].append(np.linalg.eigh(marginal)[1][:,-1])
    counts = [degree(rays,tri) for rays in (spin_rays,receiver_rays,full_rays)]
    separate = max(abs(degree(rays,tri)) for rays in individual)
    print(f"level={level} eta={sign:+d} contacts={slots} spin_collective_joint="
          f"{counts[0]:.8g},{counts[1]:.8g},{counts[2]:.8g} "
          f"max_individual_abs_degree={separate:.3g} min_individual_vacuum={worst_vacuum:.8g}")
    expected = sign if np.cos(np.sqrt(2)*COUPLING*DURATION)**(2*slots)<.5 else 0
    assert max(abs(np.array(counts)-[expected,-expected,0]))<1e-9
    assert separate<1e-9 and worst_vacuum>.5


if __name__=="__main__":
    dense_check()
    for level in (2,4):
        for face in (1,-1):
            for count in (1,3,4,6):
                scan(level,face,count)
