"""Local spectral check with the carrier impedance evaluated across the line.

Keep the central dipole and its two directly coupled bend sidebands. Other
support harmonics, sources of dephasing and intrinsic loss are omitted.
"""
import argparse
import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, minimize_scalar

from mouth_carrier_bend_work import OpenSleeve, TensorGridSpec, response
from mouth_carrier_bend_budget import settings


def dipole_widths(detunings, forces, ws, wd, g, gap, arms):
    # The choice of frequency samples covers both half-height crossings.
    def read(narrow):
        sample = np.arange(0,len(detunings),2) if narrow else np.arange(len(detunings))
        kernels=[(PchipInterpolator(detunings[sample],f[sample].real,extrapolate=False),
                  PchipInterpolator(detunings[sample],f[sample].imag,extrapolate=False)) for f in forces]
        def absorptive(d):
            z=wd+d
            self_energy=0j
            for bi, shift in enumerate((-ws,ws)):
                re,im=kernels[bi]
                D=6-12*(z+shift)**2-arms*(re(d)+1j*im(d))
                self_energy+=g**2/(4*D)
            return (1/(wd**2-z**2-self_energy)).imag
        lo,hi=detunings[sample[0]],detunings[sample[-1]]
        center=minimize_scalar(lambda d:-absorptive(d),bounds=(lo,hi),method="bounded",
                               options={"xatol":1e-13}).x
        half=absorptive(center)/2
        left=brentq(lambda d:absorptive(d)-half,lo,center)
        right=brentq(lambda d:absorptive(d)-half,center,hi)
        return center,right-left,(right-center)/(center-left)
    center,width,ratio=read(False)
    coarse=read(True)
    assert max(abs(center-coarse[0]),abs(width-coarse[1]))<2e-5
    x=np.linspace(-gap-8*width,8*width,20001)
    peaks,_=find_peaks(.5*width/(x*x+(.5*width)**2) + 1.5*width/((x+gap)**2+(1.5*width)**2))
    return center,width,ratio,len(peaks),abs(coarse[1]-width)


def calculate(cells):
    point,energy,bend,flow,archive=settings()
    p=archive["parameters"]
    ws,wd,length=[float(point[k]) for k in ("omega_s","omega_d","L")]
    g=float(bend["overlap"])*np.sqrt(2*p["I_wave"]/ws)
    gap=abs(float(energy["U_over_action"]))
    radius=archive["phase_A"]["q_star"][0]
    grid=TensorGridSpec(2*radius,2*cells,0.,length,2*cells)
    solver=OpenSleeve(grid,cells,cells,hbar=p["hbar"],constituent_mass=p["m"],
                      eos_k=p["K_eos"],bath_density=p["rho0"]/p["m"])
    rate=float(flow["reference_inward_current"])/p["m"]
    psi,phase,_=solver.solve(rate)
    assert max(abs(np.array(solver.budget(psi,phase)[:3])-rate))<1e-10
    offsets=np.array([-.016,-.012,-.009,-.0075,-.006,-.004,-.002,0.,.004])
    for ports in ("feed+receiver","receiver"):
        def force(omega):
            return response(solver,psi,phase,rate,length,omega,ports,True,2)[3]
        static=force(0.).real
        forces=[]
        for sign in (-1,1):
            row=np.array([force(wd+sign*ws+d)-static for d in offsets])
            assert min(row.imag)>0
            forces.append(row)
        for arms in (1,2):
            shift,width,ratio,peaks,error=dipole_widths(offsets,forces,ws,wd,g,gap,arms)
            print(f"lineshape cells={cells} ports={ports} arms={arms} three_tone_shift={shift:.12g} "
                  f"dipole_FWHM={width:.12g} width_over_gap={width/gap:.12g} right_over_left={ratio:.9g} "
                  f"frequency_sampling_change={error:.3g} conditional_peak_count={peaks}",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cells",type=int,required=True)
    calculate(p.parse_args().cells)
