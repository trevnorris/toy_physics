"""A supported radius-one mouth in a weak feeding-brane stream, at fixed intake."""
import argparse
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from scipy.linalg import eigvalsh
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare, weighted_gram
from mouth_free_bend_drive import Response, real_pair
from mouth_free_bend_held_flow import volume_reaction
from mouth_free_bend_incident import boundary_momentum
from mouth_free_bend_release import pencil
from mouth_wave_flare_work import ARCHIVE, SWITCHES, Wave, emit, taper


def parts(a):
    return dict(real=np.real(a).tolist(), imag=np.imag(a).tolist())


class OccupiedMode:
    """Ground-mode susceptibility per action, with both bends acting on one well."""
    def __init__(self, step, pad):
        p = yaml.safe_load(ARCHIVE.read_text())['parameters']
        self.action = p['I_wave']
        well = Wave(1., 1., step, pad, p)
        omega, u, op, _, radial = well.ground(1.)
        self.omega = omega
        if omega >= p['mu_out']:
            raise ArithmeticError('The prepared wave lies outside the bound band')
        gate = well.gate(1.)[0]
        width = p['delta_r']/2
        vr = well.contrast*gate*(1-radial)/width
        radial_source = vr.ravel()*u/np.sqrt(3)
        radius = np.broadcast_to(well.r, gate.shape).ravel()
        self.dipole = op + sparse.diags(2*p['c_w']**2/radius**2-omega**2, format='csc')
        shape = np.column_stack([np.ones(len(u))]+[
            np.broadcast_to(1-taper(well.w,*span)[0], gate.shape).ravel() for span in SWITCHES])
        self.q = radial_source[:, None]*shape
        inverse = splu(self.dipole).solve
        self.first = inverse(self.q)
        # For a radially translated well the contact term for coordinates i,j
        # uses s_i*s_j. Calibrate its common-translation identity at each depth;
        # check the uncorrected radial second derivative separately below.
        local = shape.T@((radial_source*self.first[:,0])[:,None]*shape)
        virtual = self.q.T@self.first
        self.static = (local-virtual)/omega
        vxx = (vr*((2*radial-1)/width+2/well.r)/3).ravel()
        direct = (shape.T@((.5*u*u*vxx)[:,None]*shape)-virtual)/omega
        self.mass = (self.first.T@self.first+4*omega**2*self.first.T@inverse(self.first))/omega
        emit(test='wave', step=step, padding=pad, omega=omega, action=self.action,
             static=self.static.tolist(), direct_static=direct.tolist(), mass=self.mass.tolist())
        np.testing.assert_allclose(self.static[:,0], 0., atol=2e-12)
        np.testing.assert_allclose(self.static, self.static.T, atol=2e-12)
        if min(eigvalsh(self.mass)) <= 0 or min(eigvalsh(self.static[1:,1:])) <= 0:
            raise ArithmeticError('Unresolved occupied-mode inertia or alignment')

    def force(self, rate):
        if rate == 0:
            return self.static
        result = self.static.copy()
        # Pair the positive and negative frequency sidebands. For a real
        # growth rate they are conjugate. The resolvent identity avoids
        # subtracting two nearly equal full solutions.
        shifts = [-rate**2+sign*2j*self.omega*rate for sign in (-1,1)]
        if np.isreal(rate):
            shifts = shifts[:1]
        for shift in shifts:
            shift = np.real_if_close(shift).item()
            inverse = splu(self.dipole-sparse.diags(np.full(self.dipole.shape[0],shift),
                                                    format='csc')).solve
            correction = shift*self.q.T@inverse(self.first)/(2*self.omega)
            result -= (2 if len(shifts)==1 else 1)*correction.real
        return result


class IncidentCarrier:
    def __init__(self, source):
        dom, field = prepare(argparse.Namespace(source=source))
        self.wall = Response(dom, field, 'wide')
        self.model = self.wall.dipole
        self.area = self.wall.area
        mass = weighted_gram(dom.N,dom.weight)
        skew = sparse.bmat([[None,-mass],[mass,None]], format='csr')
        boundary_field = 1j*dom.points[:,0]*field*dom.feed
        boundary = real_pair(boundary_field)
        self.field_load = -(self.model.operator@boundary)[self.model.free]
        self.field_time = (skew@boundary)[self.model.free]
        self.traction = volume_reaction(dom,field,boundary_field,'wide')
        self.momentum = boundary_momentum(dom,field,boundary_field,'wide')

    def condensed(self, rate):
        wall = self.wall
        op = wall.operator-rate*wall.skew
        solve = self.model.factor.solve if rate==0 else splu(op.tocsc()).solve
        rhs = np.column_stack((wall.bend+rate*wall.motion, self.field_load+rate*self.field_time))
        field = solve(rhs)
        if np.max(abs(op@field-rhs)) > 4e-8:
            raise ArithmeticError('Carrier field solve failed')
        coupling = (2/3)*(wall.bend-rate*wall.motion).T
        z = wall.frozen-coupling@field[:,:3]
        load = self.traction-rate*self.momentum+coupling@field[:,3]
        return z, load

    def full_check(self, rate, extra, predicted):
        energy, time = pencil(self.wall,0.)
        system = energy-rate*time+sparse.block_diag(
            (sparse.csc_matrix(self.wall.operator.shape), sparse.csc_matrix(extra)), format='csc')
        rhs = np.r_[(2/3)*(self.field_load+rate*self.field_time),self.traction-rate*self.momentum]
        independent = splu(system.tocsc()).solve(rhs)
        residual = np.max(abs(system@independent-rhs))
        emit(test='full_field', residual=residual,
             relative_coordinate_error=np.linalg.norm(independent[-3:]-predicted)/np.linalg.norm(predicted))
        np.testing.assert_allclose(independent[-3:], predicted, rtol=3e-6, atol=1e-7)
        if residual > 1e-7:
            raise ArithmeticError('Full fluid and occupied-coordinate balance failed')


def run(args):
    wave = OccupiedMode(args.step,args.pad)
    carrier = IncidentCarrier(args.source)
    sigma_area = carrier.area
    static, load = carrier.condensed(0.)
    combined = static+wave.action*wave.static
    # Relax the common translation in this finite bath before asking how
    # much independent wall area tension would remove a negative bend.
    def bends(matrix):
        return matrix[1:,1:]-np.outer(matrix[1:,0],matrix[0,1:])/matrix[0,0]
    tension_threshold = max(eigvalsh(-bends(combined),sigma_area[1:,1:]))
    action_threshold = max(eigvalsh(-bends(static),wave.static[1:,1:]))
    emit(test='static', carrier=static.tolist(), incident=load.tolist(), area=sigma_area.tolist(),
         occupied_eigen=eigvalsh(combined).tolist(), critical_tension=tension_threshold,
         critical_action_at_zero_tension=action_threshold)
    checks = [(wave.action,0.),(wave.action,.02),(wave.action,.04),(0.,.02)]
    for nu in (.01,.005,.0025):
        rate = 1j*nu
        zw = wave.force(rate)
        zc, drive = carrier.condensed(rate)
        np.testing.assert_allclose(zc,zc.T.conj(),atol=2e-9)
        emit(test='operators', frequency=nu, wave=zw.tolist(), carrier=parts(zc), drive=parts(drive))
        for action,tension in checks:
            stat = static+action*wave.static+tension*sigma_area
            if min(eigvalsh(stat)) <= 0:
                continue
            extra = action*zw+tension*sigma_area
            z = zc+extra
            x = np.linalg.solve(z,drive)
            accel = -nu**2*x[0]
            rigid = -nu**2*drive[0]/z[0,0]
            stations = np.cumsum(x)
            emit(test='driven', frequency=nu, action=action, tension=tension,
                 acceleration_per_speed=parts(accel), rigid_per_speed=parts(rigid),
                 rigid_relative_error=abs(accel-rigid)/abs(rigid),
                 in_phase_relative_change=(accel.real-rigid.real)/rigid.real,
                 relative_to_entrance=parts((stations-x[0])/x[0]),
                 flare_mismatch=abs((stations[-1]-x[0])/x[0]),
                 absolute_entrance_excursion_at_1e_minus_6=1e-6*abs(x[0]),
                 residual=np.max(abs(z@x-drive)))
            if nu==.005 and action==wave.action and tension==.02:
                carrier.full_check(rate,extra,x)

    if min(eigvalsh(combined)) < 0:
        # Look for a genuinely growing homogeneous solution. A negative
        # static stiffness alone would not settle gyroscopic stability.
        def growth_matrix(s):
            if s==0:
                return combined
            return carrier.condensed(s)[0]+wave.action*wave.force(s)
        def determinant(s):
            return np.linalg.det(growth_matrix(s))
        left, dleft = 0., np.linalg.det(combined)
        for right in (.025,.05,.10,.20,.40):
            dright = determinant(right)
            if np.sign(dleft) != np.sign(dright):
                root = brentq(determinant,left,right,xtol=1e-10)
                z = growth_matrix(root)
                _,singular,vt = np.linalg.svd(z)
                mode = vt[-1]/np.max(abs(vt[-1]))
                emit(test='unforced_growth', rate=root, bracket=[left,right],
                     determinants=[dleft,dright], mode=mode.tolist(), residual=np.max(abs(z@mode)))
                if singular[-1] > 2e-8:
                    raise ArithmeticError('Unresolved growth-mode balance')
                break
            left,dleft = right,dright
        else:
            emit(test='growth_unresolved', last_rate=left)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--step',type=float,required=True)
    parser.add_argument('--pad',type=float,required=True)
    run(parser.parse_args())
