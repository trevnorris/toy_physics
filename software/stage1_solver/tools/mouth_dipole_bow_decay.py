"""Can two passive responses of the same bow give different dipole lifetimes?

Linearize the smooth-gate material about its sustained coherent monopole,
using the measured mixed gate derivative. Prescribe identical bare mass and
stiffness for the bow in both cases. Only its damping changes. Neither is
an estimate of the actual material damping.
"""
import argparse
from pathlib import Path
from itertools import permutations

import numpy as np
from scipy.integrate import solve_ivp

from mouth_wave_linewidth_budget import row


def generator(time,ws,wd,bending,mass,stiffness,damping):
    driver=bending*np.cos(ws*time)
    return np.array([[0.,1.,0.,0.],[-wd**2,0.,-driver,0.],
                     [0.,0.,0.,1.],[-driver/mass,0.,-stiffness/mass,-damping/mass]])


def period_map(ws,wd,coupling,mass,stiffness,damping):
    period=2*np.pi/ws
    def rhs(time,flat):
        return (generator(time,ws,wd,coupling,mass,stiffness,damping)@flat.reshape(4,4)).ravel()
    solution=solve_ivp(rhs,(0,period),np.eye(4).ravel(),method="DOP853",rtol=2e-12,atol=2e-14)
    assert solution.success
    return np.linalg.eig(solution.y[:,-1].reshape(4,4))


def compare(ws,wd,overlap,damping,gap):
    mass,stiffness=12.,6.
    action=1.
    coupling=overlap*np.sqrt(2*action/ws)
    period=2*np.pi/ws
    old,vec=period_map(ws,wd,0,mass,stiffness,damping)
    marker=np.argmin(abs(old-np.exp(1j*wd*period)))
    for fraction in (.1,.25,.5,.75,1.):
        values,vectors=period_map(ws,wd,fraction*coupling,mass,stiffness,damping)
        # Follow all four branches to avoid tagging a lost mechanical mode as a dipole.
        perm=min(permutations(range(4)),key=lambda p:np.sum(abs(values[list(p)]-old)**2))
        old=values[list(perm)]
        vec=vectors[:,list(perm)]
    multiplier=old[marker]
    gamma=-2*np.log(abs(multiplier))/period
    state=vec[:,marker]
    def energy(v):
        q,qd,b,bd=v
        return (.5*(abs(qd)**2+wd**2*abs(q)**2+mass*abs(bd)**2+stiffness*abs(b)**2)
                +coupling*(q.conj()*b).real)
    def ledger(time,aux):
        q,qd,b,bd=aux[:4]
        force_rate=-ws*coupling*np.sin(ws*time)
        flow=generator(time,ws,wd,coupling,mass,stiffness,damping)@aux[:4]
        return np.r_[flow,force_rate*(q.conj()*b).real,damping*abs(bd)**2]
    evolved=solve_ivp(ledger,(0,period),np.r_[state,0.,0.],method="DOP853",rtol=2e-12,atol=2e-14)
    assert evolved.success
    final=evolved.y[:,-1]
    energy_error=energy(final[:4])-energy(state)-final[4].real+final[5].real
    def loss_response(freq):
        return (1/(stiffness-mass*freq**2-1j*damping*freq)).imag
    # Coherent classical support calculation. Both combination tones scale
    # with I_s; there is no additional spontaneous '+1' in the lower channel.
    perturbation=coupling**2/(4*wd)*(loss_response(wd-ws)+loss_response(wd+ws))
    print(f"bow damping={damping:g} bare_mass={mass:g} bare_stiffness={stiffness:g} "
          f"ImChi_diff={loss_response(wd-ws):.12g} ImChi_sum={loss_response(wd+ws):.12g} "
          f"perturbative_rate={perturbation:.12g} exact_linear_rate={gamma:.12g} "
          f"action_time={1/gamma:.12g} width_over_initial_split={gamma/gap:.12g} "
          f"max_multiplier={max(abs(old)):.12g} energy_error={abs(energy_error):.3g}")
    assert max(abs(old))<1 and abs(energy_error)<1e-9 and gamma>0


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("spectrum",type=Path)
    parser.add_argument("bend_projection",type=Path)
    args=parser.parse_args()
    p=row(args.spectrum,"point n=1 ")
    s=row(args.spectrum,"levels ")
    bend=row(args.bend_projection,"bow h=0.015 ")
    for damping in (.1,6.5):
        compare(float(p["omega_s"]),float(p["omega_d"]),float(bend["overlap"]),
                damping,abs(float(s["U_over_action"])))
