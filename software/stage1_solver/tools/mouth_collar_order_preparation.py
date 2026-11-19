"""Can an initially uniform streaming orientation acquire a signed texture?

The existing local radial pin fixes the resonant transition axis. An added
isotropic energy-loss channel couples that transition to three degenerate
receiver states. Each encounter starts the receiver in the same scalar
vacuum, and keeps its outgoing state in the energy account. Calculations
are in the resonant interaction picture, for one orienting face of a collar.
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from mouth_material_links import section, state, preferred_bonds, fluxes, count
from mouth_pair_constituent_response import PAULI


IDENTITY = np.eye(2)
VACUUM = np.diag([1.,0.,0.,0.])
INPUT = np.diag([1.,0.])
COUPLING = 1.
DURATION = .33


def pin(axis):
    spin = np.einsum("a,aij->ij",axis,PAULI)
    return (IDENTITY+spin)/2, (IDENTITY-spin)/2


def full_block(axis, duration):
    ground, excited = pin(axis)
    coupling = np.zeros((8,8),complex)
    for a, sigma in enumerate(PAULI):
        emission = np.zeros((4,4))
        emission[a+1,0] = 1
        term = np.kron(ground@sigma@excited, emission)
        coupling += COUPLING*(term+term.conj().T)
    return expm(-1j*duration*coupling)


def event(rho, block):
    joint = block@np.kron(rho,VACUUM)@block.conj().T
    tensor = joint.reshape(2,4,2,4)
    return np.trace(tensor,axis1=1,axis2=3), np.trace(tensor,axis1=0,axis2=2)


def bloch(rho):
    return np.array([np.trace(rho@sigma).real for sigma in PAULI])


def predicted(normals, eta, transfer, events):
    axis = eta*normals
    longitudinal = axis[:,2]
    remaining = transfer**events
    transverse = np.array([0.,0.,1.])-longitudinal[:,None]*axis
    return np.sqrt(remaining)*transverse+(1-(1-longitudinal)*remaining)[:,None]*axis


def calibration():
    axis = np.array([0.,0.,1.])
    _, excited = pin(axis)
    rho, bath = event(excited, full_block(axis,DURATION))
    loss = np.trace(bath[1:,1:]).real
    retained = np.trace(rho@excited).real
    print(f"held g={COUPLING:g} duration={DURATION:g} excited_retained={retained:.12g} "
          f"receiver_quanta={loss:.12g} energy_residual={abs(1-retained-loss):.3g}")
    assert abs(retained+loss-1)<1e-12
    return retained


def finite_checks(retained):
    normals = np.array([[0.,0.,1.],[.8,0.,-.6],[0.,.6,.8],[0.,0.,-1.]])
    for eta in (1,-1):
        for n in (1,3,4,6):
            target = predicted(normals,eta,retained,n)
            errors, residuals = [], []
            for direction, expected in zip(normals,target):
                axis = eta*direction
                _, excited = pin(axis)
                rho = INPUT
                block = full_block(axis,DURATION)
                quanta = 0.
                for _ in range(n):
                    rho, bath = event(rho,block)
                    quanta += np.trace(bath[1:,1:]).real
                errors.append(np.linalg.norm(bloch(rho)-expected))
                residuals.append(abs(np.trace((INPUT-rho)@excited)-quanta))
            print(f"full eta={eta:+d} encounters={n} polarization_error={max(errors):.3g} "
                  f"emitted_energy_residual={max(residuals):.3g}")
            assert max(errors+residuals)<2e-12


def measure_mesh(retained):
    for level in (2,4):
        normals, triangles = section(level)
        for eta in (1,-1):
            for n in (0,1,3,4,6):
                p = predicted(normals,eta,retained,n)
                lengths = np.linalg.norm(p,axis=1)
                assert min(lengths)>1e-6
                rays = [state(direction) for direction in p/lengths[:,None]]
                phases = fluxes(preferred_bonds(rays,triangles),triangles)
                measured = count(phases)
                quanta_per_site = (1-retained**n)/2  # uniform sphere mean initial excitation = 1/2
                print(f"mesh level={level} eta={eta:+d} encounters={n} retained_power={retained**n:.9g} "
                      f"min_Bloch_length={min(lengths):.9g} eigenray_degree={measured} "
                      f"mean_emitted_quanta_per_site={quanta_per_site:.9g}")
                assert measured==(0 if retained**n>.5 else eta)
    # Locate the spectral degeneracy independently using a full event at the
    # opposed pole, after three completed encounters and a fractional fourth.
    opposed = np.array([0.,0.,-1.])
    exact = full_block(opposed,DURATION)
    rho = INPUT
    for _ in range(3):
        rho, _ = event(rho,exact)
    def gap_signed(duration):
        pole, _ = event(rho,full_block(opposed,duration))
        return bloch(pole)[2]
    crossing = brentq(gap_signed,0.,DURATION)
    print(f"critical fourth_contact_fraction={crossing/DURATION:.10g} "
          f"pole_Bloch_z={gap_signed(crossing):.3g} eigenray_at_crossing=undefined")
    # With no energy receiver, local precession is a smooth unitary map.
    normals, triangles = section(3)
    for eta in (-1,1):
        rays = []
        for direction in normals:
            ground, excited = pin(eta*direction)
            rays.append((ground+np.exp(-2.3j)*excited)@np.array([1.,0.]))
        degree=count(fluxes(preferred_bonds(rays,triangles),triangles))
        print(f"pin_precession_only eta={eta:+d} degree={degree}")
        assert degree==0


def main():
    retained = calibration()
    finite_checks(retained)
    measure_mesh(retained)


if __name__=="__main__":
    main()
