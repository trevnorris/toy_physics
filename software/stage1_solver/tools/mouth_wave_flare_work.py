"""Axial work and lateral following of a trial bound wave on the carrier flare.

Use the archived massive scalar, with a smooth well between the waist and the
exit, in three brane dimensions and w. The neutral carrier is unchanged when
this probe is occupied. All reported forces are extra forces per wave action.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import eigsh, spsolve
from scipy.special import expit

from mouth_carrier_lateral_inertia import taper
from mouth_curved_phase import RIM_SLOPE
from mouth_flared_receiver import BORE, EXIT, FLARE_START, ROUNDING, FlaredDomain, emit
from mouth_trapped_wave_ports import ARCHIVE


RECORD = Path(__file__).resolve().parents[3] / 'research/pde_ledger_v3/_measurements/mouth_carrier_0pn.json'
SWITCHES = ((.55, 1.15), (1.2, 2.))
PROFILE = SimpleNamespace(profile='cubic', neck=BORE + ROUNDING*(np.sin(np.arctan(RIM_SLOPE))-1))


def requirements():
    runs = json.loads(RECORD.read_text())['runs']
    inputs = {}
    for name in ('one', 'inverse_half', 'inverse_two', 'inverse_four'):
        run = runs['carrier_depth_' + name + '_01']
        if run['exit_code'] or run['worker_exit'] or not run['limits_verified']:
            raise ValueError('Invalid carrier reference')
        row = json.loads(run['stdout_lines'][-1])
        if name == 'one':
            force = -row['gradient']
            inputs[1.] = (row['depth'], force, force, force)
        else:
            inputs[row['radius']] = (row['depth'], row['required_restoration'],
                                     row['pressure_only'], row['tension_only'])
    return inputs


def radial_profile(s):
    return FlaredDomain.radius(PROFILE, np.clip(s, ROUNDING, EXIT))


def binding_bound(radius, depth, p):
    # Inside the radial surface bound g_r^2 by one. Outside, use its squared
    # exponential tail. The two axial tails have constant endpoint radii.
    br, bw = p['delta_r']/2, p['delta_parallel']/2
    def radial_integral(s):
        a = radius*radial_profile(s)[0]
        return 4*np.pi*(a**3/3 + br*a*a/2 + br**2*a/2 + br**3/4)
    nodes, weights = leggauss(16)
    volume = 0.
    for lo, hi in ((ROUNDING, FLARE_START), (FLARE_START, EXIT)):
        s = lo+(hi-lo)*(1+nodes)/2
        volume += depth*(hi-lo)/2 * (weights @ radial_integral(s))
    volume += bw*(radial_integral(ROUNDING)+radial_integral(EXIT))/2
    sharp_constant = 4*np.pi*np.sqrt(6)/3
    strength = p['mu_out']**2-p['mu_in']**2
    return strength*np.sqrt(volume)/(p['c_w']**2*sharp_constant)


def laplacian(n):
    return sparse.diags((-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)), (-1, 0, 1))


class Wave:
    def __init__(self, radius, depth, step, pad, p):
        self.radius, self.depth, self.p = radius, depth, p
        self.contrast = p['mu_out']**2-p['mu_in']**2
        self.r = (step*np.arange(1, int(np.ceil((1.5*radius+pad)/step))))[:, None]
        self.w = (step*np.arange(-int(np.ceil(pad/step))+1,
                                  int(np.ceil((EXIT*depth+pad)/step))))[None, :]
        self.free = (p['c_w']/step)**2*sparse.kronsum(
            laplacian(self.w.size), laplacian(self.r.size), format='csc')
        self.seed = (self.r*np.exp(-self.r**2/(1+radius**2)-(self.w-2*depth)**2)).ravel()

    def gate(self, h):
        s = self.w/h
        b, slope = radial_profile(s)
        br, bw = self.p['delta_r']/2, self.p['delta_parallel']/2
        radial = expit((self.radius*b-self.r)/br)
        left, right = expit((self.w-ROUNDING*h)/bw), expit((EXIT*h-self.w)/bw)
        value = radial*left*right
        dh = value*((1-radial)*(-self.radius*s*slope/h)/br
                    + (EXIT*(1-right)-ROUNDING*(1-left))/bw)
        return value, dh, radial

    def ground(self, h):
        g, dh, gr = self.gate(h)
        op = self.free + sparse.diags((self.p['mu_out']**2-self.contrast*g).ravel(), format='csc')
        val, vec = eigsh(op, k=1, sigma=0., which='LM', v0=self.seed, ncv=14, tol=2e-11)
        u = vec[:, 0]
        return float(np.sqrt(val[0])), u, op, dh, gr

    def read(self, sideways):
        h = self.depth
        omega, u, op, dh, gr = self.ground(h)
        force = self.contrast*np.dot(u*u, dh.ravel())/(2*omega)
        delta = .002*h
        fd = (self.ground(h-delta)[0]-self.ground(h+delta)[0])/(2*delta)
        residual = np.linalg.norm(op@u-omega**2*u)
        emit(test='occupied_minus_empty', radius=self.radius, depth=h, omega=omega,
             bulk_threshold=self.p['mu_out'], force_per_action=force, frequency_difference_force=fd,
             eigen_residual=residual, complementary_gate_norm_weight=1.-np.dot(u*u,self.gate(h)[0].ravel()))
        if residual>1e-7 or abs(force-fd)>max(.004*abs(force),2e-5):
            raise ValueError('Eigenvalue or axial work check failed')
        if sideways:
            # The ground-state derivative occupies the l=1 channel. Compute
            # the wave-action centroid for shifting a depth-dependent center.
            dipole = op + sparse.diags(np.broadcast_to(2*self.p['c_w']**2/self.r**2,
                                      (self.r.size,self.w.size)).ravel()-omega**2,format='csc')
            dvdr = self.contrast*self.gate(h)[0]*(1-gr)/(self.p['delta_r']/2)
            def centroid(weight):
                rhs = (weight*dvdr).ravel()*u/np.sqrt(3)
                change = spsolve(dipole,rhs)
                return 2*np.sum(self.r.ravel()[:,None]*(u*change).reshape(dvdr.shape))/np.sqrt(3)
            rigid = centroid(1.)
            for lo, hi in SWITCHES:
                entrance = taper(self.w/h,lo,hi)[0]
                near, far = centroid(entrance), centroid(1-entrance)
                emit(test='action_center', switch=[lo,hi], rigid_translation=rigid,
                     entrance_only=near, flare_only=far, additivity_error=near+far-rigid)
                np.testing.assert_allclose(near+far,rigid,rtol=1e-9)
        return omega, force


def main(args):
    p = yaml.safe_load(ARCHIVE.read_text())['parameters']
    inputs = requirements()
    modes = {}
    emit(test='settings', step=args.step, padding=args.pad, wave_speed=p['c_w'],
         mu_in=p['mu_in'], mu_out=p['mu_out'], delta_r=p['delta_r'], delta_w=p['delta_parallel'])
    for r in sorted(inputs):
        h = inputs[r][0]
        ratio = binding_bound(r,h,p)
        emit(test='continuum_binding_bound', radius=r, strength_over_sobolev=ratio,
             necessary_contrast=(p['mu_out']**2-p['mu_in']**2)/ratio)
        if ratio <= 1.:
            continue
        omega, force = Wave(r,h,args.step,args.pad,p).read(sideways=r==1.)
        if omega < p['mu_out']:
            modes[r] = (omega, force)
    if 1. in modes:
        reference = inputs[1.][1]
        for r in sorted(modes.keys()-{1.}):
            _, required, pressure, tension = inputs[r]
            restore = min(pressure,tension)/reference
            growth = modes[r][1]/modes[1.][1]
            emit(test='reference_recalibration', radius=r, min_restoration_growth=restore,
                 force_growth_at_common_action=growth, empty_deficit=min(pressure,tension)-required,
                 minimum_action_ratio_at_large_occupation=restore/growth)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=float, required=True)
    parser.add_argument('--pad', type=float, required=True)
    main(parser.parse_args())
