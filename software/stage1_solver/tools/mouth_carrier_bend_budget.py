"""Receiving-geometry controls and a separate coherent-support work ledger."""
import argparse
import json
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from mouth_carrier_bend_work import OpenSleeve, TensorGridSpec, response


def settings():
    path=Path("research/pde_ledger_v3/_measurements/mouth_trapped_wave_spectrum.json")
    runs=json.loads(path.read_text())["runs"]
    def row(run,prefix):
        line,=[v for v in runs[run]["stdout_lines"] if v.startswith(prefix)]
        return dict(field.split("=",1) for field in line.split() if "=" in field)
    point=row("smooth_dipole_l015_01","point n=1 ")
    energy=row("smooth_dipole_l015_01","levels ")
    bend=row("mouth_wave_linewidth_budget01","bow h=0.015 ")
    flow=row("mouth_wave_depth_response01","archive ")
    archive=yaml.safe_load(Path("software/stage1_solver/reports/pathA_26_derrick_results.yaml").read_text())
    return point,energy,bend,flow,archive


def control(cells,ratio,lengths,static,zero):
    point,_,bend,flow,archive=settings()
    p=archive["parameters"]
    length=float(point["L"])
    radius=archive["phase_A"]["q_star"][0]
    grid=TensorGridSpec(ratio*radius,ratio*cells,0.,length*(1+lengths)/2,(1+lengths)*cells)
    s=OpenSleeve(grid,cells,cells,hbar=p["hbar"],constituent_mass=p["m"],eos_k=p["K_eos"],bath_density=p["rho0"]/p["m"])
    rate=0. if zero else float(flow["reference_inward_current"])/p["m"]
    psi,phi,error=s.solve(rate)
    flux=s.budget(psi,phi)
    assert max(abs(np.array(flux[:3])-rate))<1e-10 and flux[3]<1e-8
    frequencies=[0.] if static else [float(bend["difference"]),float(bend["sum"])]
    for ports in ("feed+receiver","receiver"):
        for omega in frequencies:
            mean,powers,cut,force=response(s,psi,phi,rate,length,omega,ports,True,2)
            print(f"control cells={cells} receiver_radii={ratio} receiver_lengths={lengths} achieved_intake={flux[0]:.12g} "
                  f"ports={ports} omega={omega:.12g} power_per_b2={mean:.12g} feed_per_b2={powers['feed']:.12g} "
                  f"cut_per_b2={cut:.12g} induced_reaction_real={force.real:.12g}",flush=True)


def pump():
    point,_,bend,_,archive=settings()
    ws,wd=float(point["omega_s"]),float(point["omega_d"])
    g=float(bend["overlap"])*np.sqrt(2*archive["parameters"]["I_wave"]/ws)
    for omega in (wd-ws,wd+ws):
        # b=cos(omega*t), q initially empty; results per b_peak squared.
        def rhs(t,y):
            q,v,_,_=y
            b=np.cos(omega*t)
            return (v,-wd**2*q-g*b*np.cos(ws*t),
                    -g*omega*np.sin(omega*t)*np.cos(ws*t)*q,
                    -g*ws*np.sin(ws*t)*b*q)
        period=2*np.pi/omega
        times=period*np.array([1,4,16])
        sol=solve_ivp(rhs,(0,times[-1]),np.zeros(4),t_eval=times,
                      method="DOP853",rtol=5e-12,atol=1e-13)
        assert sol.success
        # Analytic Duhamel solution: one resonant cosine and one off-resonant.
        other=2*ws-wd if omega==wd-ws else wd+2*ws
        for cycles,t,y in zip((1,4,16),times,sol.y.T):
            q,v,wall_work,pump_work=y
            q_check=-g*t*np.sin(wd*t)/(4*wd)-.5*g*(np.cos(other*t)-np.cos(wd*t))/(wd**2-other**2)
            e=.5*(v*v+wd**2*q*q)
            interaction=g*np.cos(ws*t)*np.cos(omega*t)*q
            residual=e+interaction-pump_work-wall_work
            assert abs(q-q_check)<2e-9*max(1,abs(q_check)) and abs(residual)<2e-9*max(1,e)
            print(f"support omega={omega:.12g} cycles={cycles} external_wall_work_per_b2={wall_work:.12g} "
                  f"pump_work_per_b2={pump_work:.12g} dipole_energy_per_b2={e:.12g} "
                  f"interaction_energy_per_b2={interaction:.12g} balance_error={residual:.3g}")


def interpret(dynamic_path,static_path):
    point,energy,bend,_,archive=settings()
    ws,wd=float(point["omega_s"]),float(point["omega_d"])
    g=float(bend["overlap"])*np.sqrt(2*archive["parameters"]["I_wave"]/ws)
    gap=abs(float(energy["U_over_action"]))
    def load(path, tag):
        return [dict(v.split("=",1) for v in line.split()[1:]) for line in path.read_text().splitlines() if line.startswith(tag)]
    static=load(static_path,"control ")
    work=load(dynamic_path,"work ")
    for ports in ("feed+receiver","receiver"):
        s0,=[r for r in static if r["ports"]==ports and float(r["omega"])==0.]
        table=[r for r in work if float(r["current"])>0 and r["ports"]==ports and r["drift"]=="1" and r["trace"]=="2"]
        assert len(table)==2
        for arms in (1,2):
            response_values=[]
            for r in table:
                omega=float(r["omega"])
                power=float(r["work_per_b2"])
                # Take previously illustrative K to be measured STATIC stiffness,
                # and subtract the induced zero-frequency carrier reaction.
                load_change=complex(float(r["reactive_F_per_b"])-float(s0["induced_reaction_real"]),2*power/omega)
                chi=1/(6-12*omega**2-arms*load_change)
                response_values.append(chi)
                if arms==1:
                    print(f"cycle ports={ports} omega={omega:.12g} mean_power_per_b2={power:.12g} "
                          f"energy_per_cycle_b1e5={(2*np.pi/omega)*power*1e-10:.12g} "
                          f"feed_fraction={float(r['feed_per_b2'])/power:.12g} effective_C={2*power/omega**2:.12g}")
            gamma=g*g/(4*wd)*sum(v.imag for v in response_values)
            shift=-g*g/(8*wd)*sum(v.real for v in response_values)
            # Equal area, separated preparation and selective readout of 01/12.
            x=np.linspace(-gap-8*gamma,8*gamma,20001)
            lines=(gamma/2)/(x*x+(gamma/2)**2)+(3*gamma/2)/((x+gap)**2+(3*gamma/2)**2)
            peaks,_=find_peaks(lines)
            print(f"transfer ports={ports} identical_arms={arms} first_order_energy_rate={gamma:.12g} "
                  f"dipole_shift={shift:.12g} loss_over_conservative_gap={gamma/gap:.12g} conditional_peaks={len(peaks)}")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(required=True,dest="mode")
    case=sub.add_parser("control")
    case.add_argument("--cells",required=True,type=int)
    case.add_argument("--radii",type=int,default=2)
    case.add_argument("--lengths",type=int,default=1)
    case.add_argument("--static",action="store_true")
    case.add_argument("--rest",action="store_true")
    sub.add_parser("pump")
    data=sub.add_parser("interpret")
    data.add_argument("dynamic",type=Path)
    data.add_argument("static",type=Path)
    args=parser.parse_args()
    if args.mode=="control":
        control(args.cells,args.radii,args.lengths,args.static,args.rest)
    elif args.mode=="pump":
        pump()
    else:
        interpret(args.dynamic,args.static)
