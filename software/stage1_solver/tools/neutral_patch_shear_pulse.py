"""One held draining mouth: measure shear pulses at two interior fluid planes.

The pulse and absolute-time calibration drives are specified independently of
companion flow or force data. The release supplies only stored initial shear
velocity. Sensor data omit the backing and use identical zero-drive histories.
"""
import argparse
import math
from pathlib import Path

import numpy as np

from mouth_order_wave_pulse import digest, emit
from neutral_patch_reciprocal_flow import OBSERVABLES, experiment

PLANES = np.array([-.5,-1.])
PULSE_CASES = ('zero','pulse','minus_pulse','half_pulse','minus_half_pulse')
AMP = .001


class Pulse:
    cases = (*PULSE_CASES,'b0','b1','b2','b3','release')

    def __init__(self, output):
        self.output = output

    def prepare(self, time):
        duration = time[-1]-time[0]
        p = (time-time[0])/duration
        bump = np.where(p < .45,np.sin(np.pi*p/.45)**2,0.)
        drives = [np.zeros_like(p),bump,-bump,bump/2,-bump/2]
        drives.extend(math.comb(3,k)*p**k*(1-p)**(3-k) for k in range(4))
        drives.append(np.zeros_like(p))
        top = AMP*np.column_stack(drives)
        initial = np.zeros(len(self.cases))
        initial[-1] = AMP
        self.speed = np.empty((len(time),len(PULSE_CASES),len(PLANES)))
        self.stress = np.empty_like(self.speed)
        emit(test='probe_definition',amplitude=AMP,pulse_ends=.45*duration,
             pulse_peaks=.225*duration,z=PLANES.tolist(),basis='cubic Bernstein on absolute transit time',
             release='unit empty-channel steady velocity profile scaled by amplitude, then zero top drive')
        return top,initial

    def sample(self, step, fluid):
        # Adjacent-cell averages at interior faces, before updating velocities.
        # tau is the positive-z shear force per unit area on the material below.
        right = np.searchsorted(fluid.z,PLANES)
        if step==0:
            np.testing.assert_allclose((fluid.z[right-1]+fluid.z[right])/2,PLANES,atol=1e-14)
            self.right = right
            self.h,self.nu,self.rho = fluid.h,fluid.nu,fluid.rho
        left_u = fluid.u[:,self.right-1,:].mean(axis=2)
        right_u = fluid.u[:,self.right,:].mean(axis=2)
        for values,store in (((right_u+left_u)/2,self.speed),
                             (fluid.rho*fluid.nu*(right_u-left_u)/fluid.h,self.stress)):
            views = values.reshape(-1,len(self.cases),len(PLANES))[:,:len(PULSE_CASES)]
            store[step] = (views-views[:,:1]).mean(axis=0)

    def finish(self,time,response):
        if not np.isfinite(self.speed).all() or not np.isfinite(self.stress).all():
            raise ValueError('Invalid internal shear gauge')
        odd = (self.speed[:,1]-self.speed[:,2])/2
        half = self.speed[:,3]-self.speed[:,4]
        odd_tau = (self.stress[:,1]-self.stress[:,2])/2
        half_tau = self.stress[:,3]-self.stress[:,4]
        peak = np.argmax(odd,axis=0)
        emit(test='probe_gauge',u_peak=odd.max(axis=0).tolist(),t_peak=time[peak].tolist(),
             tau_at_u_peak=odd_tau[peak,np.arange(len(PLANES))].tolist(),
             max_half_u_discrepancy=(np.max(abs(odd-half),axis=0)/np.max(abs(odd),axis=0)).tolist(),
             max_half_tau_discrepancy=(np.max(abs(odd_tau-half_tau),axis=0)/np.max(abs(odd_tau),axis=0)).tolist())
        temp = self.output.with_suffix('.partial.npz')
        np.savez_compressed(temp,t=time,u=self.speed,tau=self.stress,z=PLANES,
                            probe_cases=PULSE_CASES,case_names=self.cases,
                            amplitude=AMP,h=self.h,rho=self.rho,nu=self.nu,obs=OBSERVABLES,
                            response=response.mean(axis=0),response_by_seed=response)
        temp.replace(self.output)
        emit(test='probe_saved',path=str(self.output),sha256=digest(self.output))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--law',choices=('central','layered'),required=True)
    parser.add_argument('--backing',choices=('rest','free','robin'),required=True)
    parser.add_argument('--back-impedance',type=float)
    parser.add_argument('--grid',type=int,choices=(4,8),default=8)
    parser.add_argument('--refine',type=int,choices=(1,2),default=2)
    parser.add_argument('--pieces',type=int,default=160)
    parser.set_defaults(histories=4,back_gap=1.5,nu=1.,cold=False)
    args=parser.parse_args()
    experiment(args,Pulse(args.output))
