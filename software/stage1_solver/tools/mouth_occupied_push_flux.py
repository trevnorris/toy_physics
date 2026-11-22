"""Lab-frame momentum budget of a weak harmonic entrance push at fixed remote baths.

The field tangent is in coordinates carried by a lateral displacement. The
reservoir boundaries are stationary and the material wall is impermeable.
Output is complex amplitude per unit amplitude of the external force.
"""
import argparse
from pathlib import Path

import numpy as np
import yaml
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import prepare, surface, trace
from mouth_empty_bend_stiffness import quadrature, shapes
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response
from mouth_occupied_companion import OccupiedMode, parts
from mouth_two_depth_load import DEPTHS, cut_flux
from mouth_wave_flare_work import ARCHIVE, SWITCHES, Wave


def components(real_field, imag_field):
    """Cosine and sine-direction coefficients of a complex-valued physical field."""
    return real_field.real + 1j*imag_field.real, real_field.imag + 1j*imag_field.imag


def current(dom, field, tangent, volume):
    # Under x'=x+d(r,w), the Jacobian 1+partial_x d cancels the
    # inverse-Jacobian change in partial_x psi at first order.
    r, _, gr, _ = volume
    f, fr = dom.N@field, gr@(field-1)
    h, hr = dom.N@tangent, gr@tangent
    return dom.weight@((h.conj()*fr+f.conj()*(hr+2*h/r)).imag)/3


def face_points(dom, block, along, fixed):
    pg, vg, _ = dom.blocks[block]
    grid = pg if along == 'p' else vg
    nodes, weights = leggauss(8)
    half = np.diff(grid)[:, None]/2
    x = (grid[:-1, None]+half*(1+nodes)).ravel()
    other = np.full_like(x, fixed)
    p, v = (x, other) if along == 'p' else (other, x)
    return p, v, (half*weights).ravel()


def reservoir(dom, field, tangent, block, along, fixed, normal):
    p, v, weight = face_points(dom, block, along, fixed)
    fields = np.column_stack([field, tangent])
    r, _, rp, rv, wp, wv, values, radial, axial = trace(dom, block, p, v, fields)
    np.testing.assert_allclose(values[:, 1], 0., atol=1e-11)
    # Dirichlet tangent and its time derivative vanish on a remote boundary.
    # The flow part of the gradient stress is n v_x v_n. Keep the
    # isotropic action stress in the complementary, nonadvective part.
    f, h = values.T
    fr, hr = radial.T
    fw, hw = axial.T
    tr, tw = (rp, wp) if along == 'p' else (rv, wv)
    nr, nw = tw, -tr
    sign = np.sign(nr*normal[0]+nw*normal[1])
    nr, nw = nr*sign, nw*sign
    a, ar, aw = abs(f), abs(f)*(fr/f).real, abs(f)*(fw/f).real
    ur, uw = (fr/f).imag, (fw/f).imag
    g = h/f
    gr, gw = (hr-g*fr)/f, (hw-g*fw)/f
    un = nr*ur+nw*uw
    dux, dun = gr.imag+2*g.imag/r, nr*gr.imag+nw*gw.imag
    adv = a*a*((dux+2*g.real*ur)*un+ur*dun)
    br, bw = ar*g.real+a*gr.real, aw*g.real+a*gw.real
    grad = (br+2*a*g.real/r)*(nr*ar+nw*aw)+ar*(nr*br+nw*bw)
    measure = 4*np.pi/3*r*r*weight
    total, isotropic = surface(dom, fields, block, along, fixed, normal)
    result = np.array([total, measure@adv, measure@grad, isotropic])
    np.testing.assert_allclose(result[0], result[1:].sum(), rtol=1e-8, atol=2e-11)
    return result


def wall_pressure(dom, field, tangent, tangent_rate, q, qdot, block, extension):
    p, v, weight = face_points(dom, block, 'p', 0.)
    r, w, _, _, wp, _, values, radial, axial = trace(
        dom, block, p, v, np.column_stack([field, tangent, tangent_rate]))
    f, h, ht = values.T
    fr, hr = radial[:, :2].T
    fw, hw = axial[:, :2].T
    shape, sr, sw = shapes(r, w, extension, dom.far, dom.roof)
    rg, dr = taper(r, 2. if extension == 'compact' else 3.4, dom.far-.35)
    zg, dz = taper(w, 3.3 if extension == 'compact' else 4.3, dom.roof-.35)
    displacement_r = np.column_stack([dr*zg, sr])@q
    displacement_w = np.column_stack([rg*dz, sw])@q
    velocity = np.column_stack([rg*zg, shape])@qdot
    # Material Neumann condition gives Re(psi_x^* psi_n)=V_n j_x
    # at the moving wall. The net outward x momentum is isotropic stress;
    # n_x dA is unchanged by a deformation along x.
    lag = (fr.conj()*(hr-displacement_r*fr) + fw.conj()*(hw-displacement_w*fr)
           + (abs(f)**8-1)*f.conj()*h/2).real + (f.conj()*(ht-velocity*fr)).imag
    return -(4*np.pi/3*r*r*wp*weight)@lag


def wave_velocity_probe(model, step, pad, nu):
    """Direct -<phi_t grad_x phi> from separately resolved optical sidebands."""
    p = yaml.safe_load(ARCHIVE.read_text())['parameters']
    well = Wave(1., 1., step, pad, p)
    omega, u, op, _, radial = well.ground(1.)
    # Ground-state phases are arbitrary; our dipole source and u below are
    # built from the same eigenvector, independently of the impedance solve.
    gate = well.gate(1.)[0]
    radius = np.broadcast_to(well.r, gate.shape).ravel()
    d = op+sparse.diags(2*p['c_w']**2/radius**2-omega**2, format='csc')
    vr = well.contrast*gate*(1-radial)/(p['delta_r']/2)
    shape = np.column_stack([np.ones(len(u))]+[
        np.broadcast_to(1-taper(well.w, *span)[0], gate.shape).ravel() for span in SWITCHES])
    source = (vr.ravel()*u/np.sqrt(3))[:, None]*shape
    def sideband(frequency):
        return splu(d+sparse.diags(np.full(len(u), omega**2-frequency**2),
                                   format='csc')).solve(source)
    lower, upper = sideband(omega-nu), sideband(omega+nu)
    # u is r times the physical mode, normalized in the transformed L2
    # norm. Both radial endpoints have Dirichlet values for u.
    padded = np.pad(u.reshape(gate.shape), ((1, 1), (0, 0)))
    dx_u = ((padded[2:]-padded[:-2]).ravel()/(2*step)-u/radius)/np.sqrt(3)
    coeff = (omega*(lower-upper)/1j+1j*nu*(lower+upper)/2)
    np.testing.assert_allclose(omega, model.omega, atol=1e-10)
    return -(dx_u@coeff)/omega, [omega-nu, omega+nu, p['mu_out']]


def run(args):
    if args.nested_source:
        from mouth_occupied_push_baths import load_bath
        dom, field = load_bath(args.nested_source)
        emit(test='preparation', source=str(args.nested_source), far=dom.far, roof=dom.roof)
    else:
        dom, field = prepare(args)
    extension = args.extension
    wall = Response(dom, field, extension)
    wave = OccupiedMode(args.step, args.pad)
    volume = quadrature(dom)
    gauge = np.array([1., 0., 0.])
    surfaces = {
        'feed': (1, 'v', dom.far, (1., 0.)),
        'bulk_side': (4, 'p', 1., (1., 0.)),
        'roof_inner': (3, 'v', dom.roof, (0., 1.)),
        'roof_outer': (4, 'v', dom.roof, (0., 1.)),
    }
    for nu in args.omega:
        rate = 1j*nu
        op = wall.operator-rate*wall.skew
        rhs = wall.bend+rate*wall.motion
        field_per_q = splu(op.tocsc()).solve(rhs)
        zc = wall.frozen-(2/3)*(wall.bend-rate*wall.motion).T@field_per_q
        zw = wave.force(rate)
        wave_p, bands = wave_velocity_probe(wave, args.step, args.pad, nu)
        emit(test='sidebands', frequency=nu, frequencies_and_exterior_threshold=bands,
             wave_inertia_from_momentum=parts(wave_p/rate),
             wave_inertia_from_common_force=parts(zw[0]/rate**2))
        for tension in (.02, .04):
            zt = tension*wall.area
            q = np.linalg.solve(zc+wave.action*zw+zt, gauge)
            pair = np.zeros(2*len(field), complex)
            pair[wall.dipole.free] = field_per_q@q
            tan = components(*np.split(pair, 2))
            derivative = components(*np.split(rate*pair, 2))
            readings = []
            for i, h in enumerate(tan):
                qi, qdi = (q.real, (rate*q).real) if i == 0 else (q.imag, (rate*q).imag)
                flux = [reservoir(dom, field, h, *spec) for spec in surfaces.values()]
                wall_x = [wall_pressure(dom, field, h, derivative[i], qi, qdi, block, extension)
                          for block in (0, 2)]
                cuts = [cut_flux(dom, field, h, z) for z in DEPTHS]
                readings.append((np.array(flux), np.array(wall_x), np.array(cuts),
                                 current(dom, field, h, volume)))
            flux, wall_x, cuts, pcar = [a+1j*b for a, b in zip(*readings)]
            outward = flux[:, 0].sum()
            dpc, dpw = rate*pcar, rate*wave.action*(wave_p@q)
            fc, fw, ft = zc@q, wave.action*(zw@q), zt@q
            emit(test='momentum', frequency=nu, tension=tension, mesh=len(dom.blocks[2][1])-1,
                 force_reference='unit entrance push; exp(+i nu t); positive flux out of carrier',
                 coordinates=parts(q), physical_station_acceleration=parts(rate**2*np.cumsum(q)),
                 carrier_applied_generalized=parts(fc), wave_applied_generalized=parts(fw),
                 tension_applied_generalized=parts(ft), generalized_residual=parts(fc+fw+ft-gauge),
                 reservoir_names=list(surfaces),
                 reservoir_components=['total', 'particle_carriage', 'amplitude_gradient', 'isotropic'],
                 reservoir_flux=parts(flux), feed_out=parts(flux[0, 0]),
                 bulk_out=parts(sum(flux[1:, 0])), total_reservoir_out=parts(outward),
                 carrier_dPdt=parts(dpc), wave_dPdt=parts(dpw),
                 direct_momentum_sum=parts(dpc+dpw+outward),
                 carrier_minus_weak=parts(dpc+outward-fc[0]), wave_minus_weak=parts(dpw-fw[0]),
                 outward_through_moving_lip_and_conduit=parts(wall_x),
                 direct_carrier_balance=parts(dpc+outward+sum(wall_x)),
                 fixed_reflecting_planes_lateral_traction=0.,
                 internal_cut_depths=DEPTHS.tolist(),
                 internal_cut_components=['total', 'particle_carriage', 'amplitude_gradient'],
                 internal_cut_toward_bulk=parts(cuts))
            if max(abs(fc+fw+ft-gauge)) > 1e-8 or not np.isfinite([dpc, dpw, outward]).all():
                raise ArithmeticError('Applied-force or finite-momentum check failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument('--source', type=Path)
    source.add_argument('--nested-source', type=Path)
    parser.add_argument('--mesh', type=int, default=34)
    parser.add_argument('--far', type=float, default=9.)
    parser.add_argument('--roof', type=float, default=10.)
    parser.add_argument('--profile', default='cubic', choices=['cubic'])
    parser.add_argument('--rate', type=float, default=.25)
    parser.add_argument('--extension', choices=['compact', 'wide'], default='wide')
    parser.add_argument('--step', type=float, default=.03)
    parser.add_argument('--pad', type=float, default=5.)
    parser.add_argument('--omega', type=float, nargs='+', default=[.005])
    run(parser.parse_args())
