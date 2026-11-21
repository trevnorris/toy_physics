"""Periodic work on a local wave probe as fresh material crosses the planar order front."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_bvp
from scipy.special import expit, logit

from mouth_conversion_stroke_coupon import H
from mouth_order_wave_pulse import emit


LEDGER = Path('research/pde_ledger_v3/_measurements')
KAPPA, MOBILITY = .08, 5.


def settings():
    record = json.loads((LEDGER/'mouth_wave_drain_power.json').read_text())
    run = record['runs']['drain_wave_power_n34_t192_01']
    assert run['limits_verified'] and run['guarded_exit']==run['worker_exit']==0
    wave, = [r for r in map(json.loads, run['literal_stdout'].splitlines()) if r['test']=='scalar']
    # The planar material patch, including its strain and force units, is the
    # preparation of the earlier zero-slip contact calculation.
    folder = Path('_scratch/s11c/contact_supply_coupon_01')
    outcome = json.loads((folder/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    setup = json.loads((folder/'stdout').read_text().splitlines()[0])
    assert setup['kappa']==KAPPA and setup['mobility']==MOBILITY and setup['intake']==.25
    return setup, wave['omega']


def run(args):
    setup, omega = settings()
    speed, width, tau = (setup[k] for k in ('material_speed','width','tau'))
    high, low = setup['ordered_side'], setup['less_ordered_side']
    area, alpha, length = (setup[k] for k in ('area','area_per_volume','length'))
    h = args.head_ratio*setup['head_reference']
    mesh = np.linspace(0, length, 65)
    w0 = KAPPA/width**2
    def bare(c):
        return .5*w0*(c*(1-c))**2
    def unforced(z, y):
        c, cz = y
        return np.array([cz,(speed*cz/MOBILITY+c*(1-c)*(w0*(1-2*c)+6*h))/KAPPA])
    guess = expit(logit(high)-mesh/width)
    base = solve_bvp(unforced, lambda l,r: [l[0]-high,r[0]-low], mesh,
                     np.array([guess,-guess*(1-guess)/width]), tol=1e-9)
    if not base.success:
        raise ArithmeticError(base.message)
    z = np.linspace(0, length, 1025)
    c0, cz0 = base.sol(z)
    integral = lambda values: area*simpson(values, x=z, axis=-1)
    supply = setup['intake']*h*(H(high)-H(low))
    base_heat = integral((speed*cz0)**2/MOBILITY)
    base_cut = area*speed*(bare(low)-bare(high)-.5*KAPPA*(cz0[-1]**2-cz0[0]**2))
    emit(test='unforced', head=h, head_ratio=args.head_ratio, omega=omega, intake=setup['intake'],
         bare_pure_phase_energy_drop=bare(1.)-bare(0.), supplied_per_entering=h*(H(high)-H(low)),
         supply=supply, heat=base_heat, outward_order_work=base_cut, transit=length/speed)
    np.testing.assert_allclose(supply,base_heat+base_cut,atol=5e-12)
    count = args.phases
    theta = 2*np.pi*np.arange(count)/count
    slip = np.cos(theta)[:,None]
    velocity = -omega*np.sin(theta)[:,None]
    freq = np.fft.fftfreq(count)*count
    time_derivative = omega*np.fft.ifft(1j*freq[:,None]*np.fft.fft(np.eye(count),axis=0),axis=0).real
    contact = args.law=='contact'
    period = 2*np.pi/omega
    for amplitude in args.amplitudes:
        a2 = amplitude**2
        beta = args.beta
        def dynamics(position, state):
            background, slope = base.sol(position)
            b, bz = state[:count], state[count:2*count]
            delta = a2*b
            c = background+delta
            hp = 6*c*(1-c)
            db = time_derivative @ b + speed*bz
            dc = speed*slope+a2*db
            # Divided difference of the derivative of the bare + head energy.
            du = b*(w0*(1-6*background+6*background**2+(-3+6*background)*delta+2*delta**2)
                    +6*h*(1-2*background-delta))
            if contact:
                extension = state[2*count:]
                reaction = alpha*extension**2*(.5-beta*c*hp)
                de = (velocity-time_derivative @ extension-(1/tau+beta*hp*dc)*extension)/speed
            else:
                reaction = .5*slip**2
            order = [bz, (db/MOBILITY+du+reaction)/KAPPA]
            return np.vstack(order+[de] if contact else order)

        def boundary(left, right):
            return np.r_[left[:count], right[:count], left[2*count:]]

        sol = solve_bvp(dynamics,boundary,mesh,np.zeros((count*(2+int(contact)),len(mesh))),
                        tol=3e-7,max_nodes=1800)
        if not sol.success:
            raise ArithmeticError(sol.message)
        state = sol.sol(z)
        b, bz = state[:count], state[count:2*count]
        c = c0+a2*b
        cz = cz0+a2*bz
        dc = speed*cz+a2*(time_derivative @ b)
        hp = 6*c*(1-c)
        mean = lambda values: np.mean(values,axis=0)
        if contact:
            extension = state[2*count:]
            work = -integral(alpha*a2*mean(c*extension*velocity))
            elastic = alpha*a2*mean(c*extension**2)/2
            renew = integral(alpha*a2*mean(c*extension**2)/tau)
            from_order = integral(alpha*a2*mean(extension**2*(.5-beta*c*hp)*dc))
        else:
            # The stationary background force is exactly in quadrature.
            work = -integral(a2**2*mean(b*slip*velocity))
            elastic = a2*mean(c*slip**2)/2
            renew = 0.
            from_order = integral(a2*mean(slip**2*dc)/2)
        export = area*speed*(elastic[-1]-elastic[0])
        heat = integral(mean(dc**2)/MOBILITY)
        cut = area*speed*(bare(low)-bare(high)-.5*KAPPA*(mean(cz**2)[-1]-mean(cz**2)[0]))
        emit(test='periodic', law=args.law,beta=beta,head_ratio=args.head_ratio,
             strain_amplitude=amplitude,phases=count,axial_nodes=len(sol.x),period=period,
             wave_mean_power=work,wave_cycle_work=period*work,work_per_passed_constituent=work/setup['intake'],
             head_power=supply,order_to_elastic=from_order,order_heat=heat,renewal_heat=renew,
             outward_order_work=cut,carried_elastic_out_minus_in=export,
             contact_balance=from_order-work-renew-export,order_balance=supply-from_order-heat-cut,
             min_order=c.min(),max_order=c.max(),maximum_material_order_rate=dc.max(),
             zero_drive_tagged_extension_gain=np.exp(beta*(H(high)-H(low))-length/(speed*tau)) if contact else None)
        np.testing.assert_allclose(from_order,work+renew+export,rtol=2e-7,atol=1e-13)
        np.testing.assert_allclose(supply,from_order+heat+cut,rtol=0,atol=2e-11)
        if c.min() < -1e-7 or c.max()>1+1e-7 or dc.max()>=0:
            raise ValueError('Material left the converting branch')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--law',choices=('shear','contact'),required=True)
    parser.add_argument('--beta',type=float,default=0.)
    parser.add_argument('--head-ratio',type=float,required=True)
    parser.add_argument('--phases',type=int,required=True)
    parser.add_argument('--amplitudes',type=float,nargs='+',required=True)
    args = parser.parse_args()
    if args.phases<9 or args.phases%2==0 or min(args.amplitudes)<=0:
        parser.error('Use an odd number of phases >=9 and positive strain amplitudes.')
    run(args)
