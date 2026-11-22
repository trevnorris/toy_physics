"""Fixed-action wave force for an axial shift of the order front in a held mouth.

Open torsional resonance on the verified fixed-intake carrier. Both squared curl
energy and transverse inertia use the same displaced order field. Only material
order in the conduit moves. The rim, geometric walls and distant cuts stay put.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh

from mouth_curved_phase import divided
from mouth_native_wave_bend import find_pole
from mouth_shear_pattern_release import (angular_stiffness, assemble, digest, emit,
    impedance, linear, momentum_and_mass, prepared, product, quadratic, tensor, truncate)


def smooth(t):
    s = np.clip(t, 0.0, 1.0)
    return s**3 * (10.0 - 15.0*s + 6.0*s*s)


class Front:
    def __init__(self, dom, order):
        p, v, cells = dom.order_blocks[2]
        self.w = divided(p)
        i, j = np.indices((len(p)-1, len(v)-1))
        self.ids = np.empty((2*len(p)-1, 2*len(v)-1), dtype=int)
        for a in range(3):
            for b in range(3):
                self.ids[2*i+a, 2*j+b] = cells[:,3*a+b].reshape(i.shape)
        self.base = order
        self.profile = CubicSpline(self.w, order[self.ids], axis=0)
        self.envelope = smooth((self.w-.6)/.25) * smooth((1.79-self.w)/.20)

    def fields(self, shift):
        argument = self.w - shift*self.envelope
        vals = self.profile(argument)
        deriv = -self.profile(argument, 1)*self.envelope[:, None]
        field = self.base.copy()
        field[self.ids] = vals
        derivative = np.zeros_like(field)
        derivative[self.ids] = deriv
        return field, derivative

    def half_depths(self, field):
        # Reconstruct half-level positions using the actual nodal quadratic
        # interpolation along the wall and axis used by the shear solver.
        def depth(column):
            vals = field[self.ids[:, column]]
            i = np.flatnonzero((vals[::2][:-1]>.5)&(vals[::2][1:]<=.5)).item()
            def level(w):
                ref = 2*(w-self.w[2*i])/(self.w[2*i+2]-self.w[2*i])-1
                return quadratic(np.array([ref]))[0][0] @ vals[2*i:2*i+3] - .5
            return brentq(level,self.w[2*i],self.w[2*i+2])
        return [depth(0), depth(-1)]


def mean_force(dom, block, order, derivative, wave, pole, degree, energy):
    """Physical virtual work: independently integrate strain and kinetic work densities."""
    p0, v0, refs = dom.order_blocks[2]
    p, v, ids = block
    roots, weights = leggauss(3)
    xi, eta = (a.ravel() for a in np.meshgrid(roots,roots,indexing='ij'))
    shape = product(linear(xi),linear(eta))
    dlin = np.tile([-.5,.5], (len(xi),1))
    sp, sv = product(dlin,linear(eta)), product(linear(xi),dlin)
    i, j = (a.ravel() for a in np.indices((len(p)-1,len(v)-1)))
    moving = (p[i]<1.82)&(p[i+1]>.57)
    i, j, ids = i[moving], j[moving], ids[moving]
    hp, hv = np.diff(p)[i,None]/2, np.diff(v)[j,None]/2
    z, b = p[i,None]+hp*(1+xi), v[j,None]+hv*(1+eta)
    old_i = np.searchsorted(p0,(p[i]+p[i+1])/2)-1
    old_j = np.searchsorted(v0,(v[j]+v[j+1])/2)-1
    a = quadratic(2*(z-p0[old_i,None])/np.diff(p0)[old_i,None]-1)[0]
    c = quadratic(2*(b-v0[old_j,None])/np.diff(v0)[old_j,None]-1)[0]
    basis = tensor(a,c)
    old_ids = refs[old_i*(len(v0)-1)+old_j]
    raw = np.einsum('eqi,ei->eq',basis,order[old_ids])
    dchi = np.einsum('eqi,ei->eq',basis,derivative[old_ids])*((raw>0)&(raw<1))
    r, _, rp, rv, wp, wv = dom.chart(2,z,b)
    det = rp*wv-rv*wp
    local = wave[ids]
    f = local@shape.T
    fp, fv = (local@sp.T)/hp, (local@sv.T)/hv
    fr, fw = (wv*fp-wp*fv)/det,(rp*fv-rv*fp)/det
    wt = (8*np.pi/3)*r**2*abs(det)*hp*hv*np.outer(weights,weights).ravel()*dchi
    strain = wt*(degree*(degree+1)*abs(f)**2+abs(2*f+r*fr)**2+r*r*abs(fw)**2)
    kinetic = wt*r**2*(abs(local)**2 @ shape.T)*abs(pole)**2
    scale = pole.real/(4*energy)
    return -scale*strain.sum(), scale*kinetic.sum()


def run(a):
    raw, base, intake, source = prepared(a)
    front = Front(raw, base)
    emit(test='preparation', intake=intake, radius=a.radius, mesh=a.mesh, subdivisions=a.sub,
         cut_request=a.cut, outer=a.outer, feed_spacing=a.feed_spacing, degree=a.degree,
         stored_action_units=1, state=str(source), state_sha256=digest(source),
         pure_translation_depths=[.85,1.59], virtual_shift_vanishes_outside=[.6,1.79])
    rows = {}
    for q in (0.0, -.02, -.01, .01, .02):
        shifted, dchi = front.fields(q)
        wall, axis = front.half_depths(shifted)
        dom, chi, cut = truncate(raw, shifted, a.cut)
        pts, blocks, _, parts, feed, _ = assemble(dom,chi,a.sub,1,a.feed_spacing)
        masses, _ = momentum_and_mass(dom,blocks,chi,len(pts),0.)
        extra = angular_stiffness(dom,blocks,chi,len(pts))
        k = sum(parts)+(a.degree*(a.degree+1)-2)*sum(extra)
        cz, _, _ = impedance(dom,blocks,pts,chi,cut,0.)
        m = masses.sum(axis=0)
        live = (m>0)&((pts[:,0]>0)|(a.degree==1))
        mass = m[live]
        scale = sparse.diags(1/np.sqrt(mass))
        stiff = (scale@k[live][:,live]@scale).tocsr()
        damping = (feed[live]+cz[live])/mass
        behind = pts[live,1]>.5000000001
        val, vec = eigsh(stiff[behind][:,behind],k=1,sigma=0,tol=1e-12,
                         v0=np.ones(behind.sum()))
        seed = np.zeros(live.sum())
        seed[behind] = vec[:,0]
        pole, normed, resid = find_pole(stiff,damping,np.sqrt(val[0])+1e-4j,seed)
        energy = (abs(pole)**2+np.vdot(normed,stiff@normed).real)/4
        phys = np.zeros(len(pts),dtype=complex)
        phys[live] = normed/np.sqrt(mass)
        elastic, kinetic = mean_force(dom,blocks[2],chi,np.r_[dchi,0.],phys,pole,
                                      a.degree,energy)
        if q == 0:
            base_normed = normed
            base_live = live
        if not np.array_equal(live,base_live):
            raise ArithmeticError('Active mesh changed')
        overlap = float(abs(np.vdot(normed,base_normed))**2)
        record = dict(q=q, waist_to_wall_half=wall-.5, waist_to_axis_half=axis-.5,
                      actual_cut=cut, frequency=pole.real, amplitude_loss=-pole.imag,
                      wave_force_per_action=elastic+kinetic, shear_part=elastic,
                      kinetic_part=kinetic, empty_wave_force=0.0,
                      mass_norm_overlap_squared=overlap, closed_seed_overlap_squared=abs(seed@normed)**2,
                      kinetic_fraction_behind_waist=float(np.sum(abs(normed[behind])**2)),
                      residual=resid)
        emit(test='front',**record)
        if (not .90<overlap<1.000001 or resid>2e-8 or
                not (0.0 <= -pole.imag < .03*pole.real)):
            raise ArithmeticError('Mode identification or residual failed')
        rows[q] = record
    compare(rows)


def compare(rows):
    for h in (.02,.01):
        minus, zero, plus = [rows[q] for q in (-h,0.,h)]
        omega = [r['frequency'] for r in (minus,zero,plus)]
        fdforce = -(omega[2]-omega[0])/(2*h)
        stiffness = (omega[2]+omega[0]-2*omega[1])/h**2
        stress_k = -(plus['wave_force_per_action']-minus['wave_force_per_action'])/(2*h)
        emit(test='signed_comparison', excursion=h, fixed_action_force_secant=fdforce,
             frequency_curvature_per_action=stiffness, direct_stress_stiffness_per_action=stress_k,
             zero_shift_direct_force=zero['wave_force_per_action'],
             endpoint_forces_have_opposite_sign=bool(minus['wave_force_per_action']*plus['wave_force_per_action']<0))
        if not np.isclose(zero['wave_force_per_action'],fdforce,rtol=.025,atol=1e-7):
            raise ArithmeticError('Virtual work and fixed-action energy derivative disagree')


def recover(directory):
    # The original summary printer rejected a NumPy boolean after all five
    # eigensolves had emitted complete rows. Do not repeat those solves.
    source = directory.parent/'solves_original.py'
    if digest(source) != 'b273775b4e993a4dd2b1ee8249d1d7aa08a70c5c6f0409086865464333db5455':
        raise ValueError('Unexpected source for recovered rows')
    guard = json.loads((directory/'outcome.json').read_text())
    failure = (directory/'stderr').read_text()
    if not (guard['limitsVerified'] and guard['exitCode']==guard['childOutcome']['exitCode']==1
            and guard['childOutcome']['guardReason'] is None and
            "Object of type bool_ is not JSON serializable" in failure and
            "emit(test='signed_comparison'" in failure):
        raise ValueError('The previous run did not reach the identified summary error')
    header, *fields = map(json.loads,(directory/'stdout').read_text().splitlines())
    if header['test']!='preparation' or digest(Path(header['state']))!=header['state_sha256']:
        raise ValueError('Order source mismatch')
    if [line['q'] for line in fields] != [0,-.02,-.01,.01,.02]:
        raise ValueError('Incomplete signed rows')
    for line in fields:
        if not (line['test']=='front' and .9<line['mass_norm_overlap_squared']<1.000001
                and line['residual']<2e-8 and 0<=line['amplitude_loss']<.03*line['frequency']):
            raise ValueError('Invalid saved branch')
    emit(test='recovered', original_source=digest(source), path=str(directory),
         raw_stdout_sha256=digest(directory/'stdout'), verified_signed_rows=len(fields))
    compare({line['q']:line for line in fields})


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.set_defaults(radius=1.,state=None)
    parser.add_argument('--mesh',type=int,default=34)
    parser.add_argument('--sub',type=int,default=1)
    parser.add_argument('--cut',type=float,default=2.15)
    parser.add_argument('--outer',type=float,default=12.)
    parser.add_argument('--feed-spacing',type=float,default=.09)
    parser.add_argument('--degree',type=int,choices=(1,4),default=4)
    parser.add_argument('--complete-rows',type=Path)
    args = parser.parse_args()
    if args.complete_rows:
        recover(args.complete_rows)
    else:
        run(args)
