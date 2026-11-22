"""Count sites formed in a thin capture slice under opposite weak incident flows.

An undepleted incoming population samples the near-wall band. Each entrant
is captured with constant hazard k while traversing a short fixed w span L.
Bound residents keep their arrival angle and have the same mean lifetime at
all angles. The product k L has units of velocity. This is a capture model,
not a measurement of the chance to bind to a physical mouth.
"""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from mouth_carrier_lateral_response import Dipole
from mouth_curved_attachment_supply import section
from mouth_curved_companion_order import PRIOR as FRONT_RECORD
from mouth_curved_site_front import MEASURE, front_read, prepare, sha
from mouth_curved_site_pair import validated
from mouth_flared_receiver import emit
from mouth_segment_contact_hold import PITCH


RECORD = MEASURE/'mouth_curved_companion_order.json'


def recorded(name):
    run=json.loads(RECORD.read_text())['guarded_runs'][name]
    assert run['limits_verified'] and run['child_exit']==run['guard_exit']==0
    return list(map(json.loads,run['literal_stdout'].splitlines()))


class Gate:
    def __init__(self, dom, fields, depth, width):
        self.f,self.fw,self.weights,_=section(dom,fields,depth,width,12)
        psi,dpsi=self.f.T
        v,dv=self.fw.T
        self.n=abs(psi)**2
        self.j=(psi.conj()*v).imag
        self.dn=2*(psi.conj()*dpsi).real
        self.dj=(psi.conj()*dv+dpsi.conj()*v).imag
        self.throughput=self.weights@self.j
        self.mean_speed=self.throughput/(self.weights@self.n)
        if min(self.j)<=0:
            raise ValueError('The local passage law requires an inward stream')

    def analytic(self, optical):
        q=optical*self.mean_speed*self.n/self.j
        p=-np.expm1(-q)
        count=self.weights@(self.j*p)
        change=self.dj*p + self.j*np.exp(-q)*q*(self.dn/self.n-self.dj/self.j)
        return count/self.throughput, (self.weights@change)/count

    def signed(self,optical,drive):
        x,weights=leggauss(10)
        f=self.f[:,0,None]+drive*self.f[:,1,None]*x
        fw=self.fw[:,0,None]+drive*self.fw[:,1,None]*x
        j=(f.conj()*fw).imag
        n=abs(f)**2
        if j.min()<=0:
            raise ValueError('Outward crossing in reconstructed tangent field')
        counts=self.weights@(j*-np.expm1(-optical*self.mean_speed*n/j))
        return (weights@counts)/2, ((weights*x)@counts)/(weights@counts)


def run(mesh, folder):
    data=recorded('mesh'+str(mesh))
    supply=validated(folder)
    dom,psi,chi=prepare(mesh)
    dipole,residual=Dipole(dom,psi).steady(1.,0.)
    fields=np.column_stack([psi,dipole])
    emit(test='preparation_checks',mesh=mesh,prior_sha256=sha(RECORD),
         prior_isolated_sha256=sha(FRONT_RECORD),supply_sha256=sha(folder/'stdout'),
         carrier_free_residual=residual)
    for loc,depth in [('entry',1.175),('half',front_read(dom,chi,[],1.)[0])]:
        material=[r for r in data if r['test']=='binding' and r['location']==loc]
        if len(material)!=2:
            raise ValueError('Expected the two relaxed adsorption laws')
        for band,width in [('half_pitch_unit',PITCH/2),('one_pitch_unit',PITCH),
                           ('two_pitch_units',2*PITCH)]:
            gate=Gate(dom,fields,depth,width)
            supply_row,=[r for r in supply if r.get('test')=='band' and
                         r['mesh']==mesh and r['location']==loc and r['band']==band]
            key=dict(mesh=mesh,location=loc,band=band)
            slow=gate.analytic(1e-8)[1]
            fast=gate.analytic(80.)[1]
            np.testing.assert_allclose([slow,fast],[supply_row['inventory_beta'],supply_row['rate_beta']],
                                       atol=2e-8,rtol=0.)
            for row in material:
                critical=row['fractional_count_dipole_per_U_for_zero_load']
                initial=row['fixed_count_medium_force_per_mass_U']
                multiplier=row['count_dipole_force_coefficient']
                grid=np.linspace(-16.,4.,121)
                diff=np.array([gate.analytic(np.exp(v))[1]-critical for v in grid])
                brackets=np.flatnonzero(diff[1:]*diff[:-1]<0)
                emit(test='counts',**key,law=row['law'],width=width,critical_count=critical,
                     force_per_mass_U_slow=initial+multiplier*slow,
                     force_per_mass_U_fast=initial+multiplier*fast,
                     beta_slow=slow,beta_fast=fast,reversals_in_sampled_range=len(brackets))
                for i in brackets:
                    logroot=brentq(lambda v:gate.analytic(np.exp(v))[1]-critical,
                                   grid[i],grid[i+1],xtol=1e-13)
                    optical=np.exp(logroot)
                    fraction,beta=gate.analytic(optical)
                    signed=[]
                    for drive in (.002,.001):
                        gain,left=gate.signed(optical,drive)
                        loss,right=gate.signed(optical,-drive)
                        np.testing.assert_allclose(gain,loss,atol=1e-13)
                        signed.append(3*(left-right)/(2*drive))
                    emit(test='capture_threshold',**key,law=row['law'],width=width,
                         signed_beta=signed,analytic_beta=beta,mean_speed=gate.mean_speed,
                         k_times_L=gate.mean_speed*optical,mean_speed_optical_depth=optical,
                         zero_flow_fraction_of_band_crossings_captured=fraction,
                         incident_fraction_of_whole_section=fraction*supply_row['band_fraction'])
                    np.testing.assert_allclose(signed,beta,atol=2e-6,rtol=1e-5)


if __name__=='__main__':
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(26,34,48),required=True)
    parser.add_argument('--supply',type=Path,required=True)
    args=parser.parse_args()
    run(args.mesh,args.supply)
