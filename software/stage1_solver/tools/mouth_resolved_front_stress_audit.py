"""Test the resolved front traction by virtual work in the unreduced Hamiltonian."""
import numpy as np
from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA, WELL
from mouth_poly_shallow_release import emit


def physical_energy(density, order, grad_f, grad_phase, grad_order):
    # Minimize over order momentum analytically in the original Hamiltonian.
    gradient_sq=grad_order@grad_order
    pi=density*(grad_phase@grad_order)/(density/ORDER_INERTIA+gradient_sq)
    velocity=grad_phase-pi/density*grad_order
    energy=(density*(velocity@velocity)/2 +pi**2/(2*ORDER_INERTIA)+grad_f@grad_f/2
            +KAPPA*gradient_sq/2+WELL*(order*(1-order))**2+(density**5-5*density+4)/20)
    return energy,velocity


def main():
    rng=np.random.default_rng(114303)
    max_error=0.
    total=0
    for density in (.012,.15,.6,.92):
        for _ in range(4):
            f,phase,q=rng.normal(size=(3,4))
            f[2],q[2]=0.,0.
            s=rng.uniform(.15,2.45)
            phase[2]=1/s
            order=rng.uniform(.05,.95)
            q*=8.
            energy,velocity=physical_energy(density,order,f,phase,q)
            for component in (0,1,3):
                variation=np.zeros((4,4))
                variation[3,component]=1
                def transformed(h):
                    F=np.eye(4)+h*variation
                    def push(g):
                        return np.linalg.solve(F.T,g)
                    return np.linalg.det(F)*physical_energy(density,order,push(f),push(phase),push(q))[0]
                # Axial outward momentum differs in sign from energy under the embedding shift.
                h=1.e-5
                measured=(transformed(-h)-transformed(h))/(2*h)
                stress=density*velocity[3]*velocity[component]+f[3]*f[component]+KAPPA*q[3]*q[component]
                if component==3:stress-=energy
                max_error=max(max_error,abs(measured-stress)/(1+abs(stress)))
                total+=1
    emit(kind='stress_virtual_work',virtual_deformations=total,max_relative_to_one_plus_abs_stress=max_error)
    if max_error>1e-8:
        raise ValueError('Hamiltonian work and tested normal stress disagree')


if __name__=='__main__':
    main()
