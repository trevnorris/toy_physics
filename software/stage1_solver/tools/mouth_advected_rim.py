"""Stationary transport of a thin converting layer over a three-ball, R=1."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.interpolate import BSpline, RegularGridInterpolator
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from incident_flow_falloff import checked_grid, check_current
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4
from mouth_front_pressure_dipole import grid, VOLUME


EVIDENCE = Path('research/pde_ledger_v3/_measurements/mouth_front_pressure_dipole.json')
RATE = .25
TENSION = 1/3  # Trial units; also obtained from the independent flat order energy below.


def face_trace(path, pins, slab, centers):
    with verified_packet(path, pins) as packet:
        checked_grid(slab, packet)
        if 'centers' in packet:
            np.testing.assert_array_equal(packet['centers'], centers)
        else:
            assert int(packet['d']) == -centers[0]
        currents = packet['currents']
        assert len(currents) == len(centers)
        for i, center in enumerate(centers):
            check_current(slab, center, currents[i:i+1])
        ix = abs(slab.x)<1+5*slab.dx
        ir = slab.r<1+5*slab.dx
        face = packet['layers' if 'layers' in packet else 'field'][-1][ix][:,ir]
        current = currents[-1][ix][:,ir]
    x, r = slab.x[ix], slab.r[ir]
    density = abs(face)**2
    j_x = (face.conj()*dx4(face,slab.dx)).imag
    # Scalars reflect evenly about r=0; the transverse radial current is odd.
    padded = np.concatenate((face[:,1::-1],face),axis=1)
    dr_face = dx4(padded.T,slab.dx).T[:,2:]
    j_r = (face.conj()*dr_face).imag
    raw = np.stack((density,-j_x,j_r),axis=-1)[2:-2,:-2]
    axis = raw[:,:1].copy()
    axis[:,:,2] = 0
    interp = RegularGridInterpolator((x[2:-2],np.r_[0,r[:-2]]),np.concatenate((axis,raw),axis=1))
    return SimpleNamespace(interp=interp, density=density, current=current, x=x, r=r)


def archive_drive(record):
    slab = grid(*record['grid'])
    pins = {Path(name):digest for name,digest in record['inputs'].items()}
    paths = list(pins)
    assert len(paths) == 2 and 'baseline' in paths[0].stem or (len(paths)==2 and paths[0].stem=='d0')
    reference = face_trace(paths[0],pins,slab,[0])
    pair = face_trace(paths[1],pins,slab,[-record['d'],0])
    weight = np.maximum(1-reference.x[:,None]**2-reference.r**2,0)**3
    def prepare(trace):
        profile = weight*trace.density/reference.density
        factor = trace.current.max()/profile.max()
        np.testing.assert_allclose(trace.current, factor*profile, atol=2e-15, rtol=1e-12)
        def sample(radius, mu):
            points = np.column_stack((-radius*mu,radius*np.sqrt(1-mu**2)))
            density, j_axis, j_transverse = trace.interp(points).T
            q = factor*np.maximum(1-radius**2,0)**3*density/reference.interp(points)[:,0]
            return np.array([q,j_axis,j_transverse])
        return sample
    return prepare(reference), prepare(pair)


class Sphere:
    def __init__(self, elements, modes, angular_nodes):
        self.modes = modes
        breaks = np.linspace(0,1,elements+1)
        knots = np.r_[np.zeros(3),breaks,np.ones(3)]
        radial_count = len(knots)-4
        self.basis = BSpline(knots,np.eye(radial_count),3)
        self.radial_count = radial_count
        self.mu, self.angular_weights = leggauss(angular_nodes)
        self.p = np.stack([Legendre.basis(l)(self.mu) for l in range(modes)],axis=1)
        self.dp = np.stack([Legendre.basis(l).deriv()(self.mu) for l in range(modes)],axis=1)
        self.sin = np.sqrt(1-self.mu**2)
        roots,weights = leggauss(5)
        self.r = (breaks[:-1,None]+(roots+1)/(2*elements)).ravel()
        wr = np.tile(weights/(2*elements),elements)
        self.w = 2*np.pi*(wr*self.r**2)[:,None]*self.angular_weights
        self.br, self.db = self.basis(self.r), self.basis(self.r,1)
        self.indices = np.arange(radial_count*modes).reshape(radial_count,modes)
        self.free = np.ones_like(self.indices,dtype=bool)
        self.free[-1] = False
        self.free[0,1:] = False

    def prepare(self, sample, normalize=True):
        rr,mm = np.broadcast_arrays(self.r[:,None],self.mu[None,:])
        q,ja,jt = sample(rr.ravel(),mm.ravel()).reshape(3,*rr.shape)
        projected_rate = np.sum(self.w*q)
        if normalize:
            q = q*(RATE/projected_rate)
        matrix_rows, matrix_cols, diffusion, drift = [],[],[],[]
        rhs = np.zeros(self.indices.size)
        for i,r in enumerate(self.r):
            active = np.flatnonzero(self.br[i])
            ids = self.indices[active].ravel()
            def outer(a,b):
                return (a[None,:,None]*b[:,None,:]).reshape(len(self.mu),-1)
            basis = outer(self.br[i,active],self.p)
            radial = outer(self.db[i,active],self.p)
            theta = outer(self.br[i,active],-self.sin[:,None]*self.dp)/r
            jr = ja[i]*self.mu+jt[i]*self.sin
            jtheta = -ja[i]*self.sin+jt[i]*self.mu
            w = self.w[i,:,None]
            matrix_rows.extend(np.repeat(ids,len(ids)))
            matrix_cols.extend(np.tile(ids,len(ids)))
            diffusion.extend((radial.T@(w*radial)+theta.T@(w*theta)).ravel())
            drift.extend((-(jr[:,None]*radial+jtheta[:,None]*theta).T@(w*basis)).ravel())
            rhs[ids] += basis.T@(self.w[i]*q[i])
        coords = (matrix_rows,matrix_cols)
        size = (self.indices.size,)*2
        def sparse(values):
            return coo_matrix((values,coords),shape=size).tocsr()
        q_edge, ja_edge, jt_edge = sample(np.ones_like(self.mu),self.mu)
        return SimpleNamespace(diffusion=sparse(diffusion),drift=sparse(drift),rhs=rhs,q=q,ja=ja,
            projected_rate=projected_rate,edge_jr=ja_edge*self.mu+jt_edge*self.sin,edge_q=q_edge)

    def solve(self, data, K, transport=True):
        matrix = K*data.diffusion
        if transport:
            matrix = matrix+data.drift
        active = self.free.ravel()
        solution = np.zeros(self.indices.size)
        solution[active] = spsolve(matrix[active][:,active],data.rhs[active])
        coefficients = solution.reshape(self.indices.shape)
        h = (self.br@coefficients)@self.p.T
        s = -(self.basis(1,1)@coefficients)
        edge_slope = self.p@s
        edge_q = data.edge_q + edge_slope*data.edge_jr*int(transport)
        radial = (self.db@coefficients)@self.p.T
        angular = ((self.br@coefficients)@self.dp.T)*self.sin/self.r[:,None]
        lateral = np.sum(self.w*h*data.ja)*int(transport)
        normal = np.sum(self.w*self.r[:,None]*self.mu*data.q)
        return dict(s1=float(s[1]),mean_s=float(s[0]),moment=float(normal+lateral),
                    flat_current_moment=float(normal),collar_moment=float(lateral),
                    edge_projection_error=float(K*VOLUME*s[1]-normal-lateral),
                    edge_q_range=[float(edge_q.min()),float(edge_q.max())],
                    edge_jr_range=[float(data.edge_jr.min()),float(data.edge_jr.max())],
                    edge_slope_range=[float(edge_slope.min()),float(edge_slope.max())],
                    max_slope=float(np.hypot(radial,angular).max()),height=float(h.max()),
                    flat_interpolated_rate=float(data.projected_rate))


def sign_control():
    # An exactly known radial parabola and a separately manufactured dipole.
    sphere = Sphere(12,4,24)
    def sample(r,mu):
        return np.array([1+.15*r*mu,np.zeros_like(r),np.zeros_like(r)])
    result = sphere.solve(sphere.prepare(sample,normalize=False),.5)
    print('polynomial slopes and analytic target',result['mean_s'],result['s1'],2/3,.06,flush=True)
    np.testing.assert_allclose([result['mean_s'],result['s1']],[2/3,.06],atol=2e-10)
    # Constant tangential stream, only uniform q; test drift sign at zero stream.
    def driven(sign):
        def sample(r,mu):
            return np.array([np.ones_like(r),sign*.0001*np.ones_like(r),np.zeros_like(r)])
        return sphere.solve(sphere.prepare(sample,normalize=False),.5)['s1']
    slope = (driven(1)-driven(-1))/.0002
    print('stream derivative, analytic target',slope,1/(15*.5**2),flush=True)
    np.testing.assert_allclose(slope,1/(15*.5**2),rtol=2e-7)
    # Independent 1-D equilibrium tension in trial units; no companion data.
    kappa,delta=.04,.02
    z, w = leggauss(220)
    gradient = 1/(4*delta*np.cosh(25*z/2)**2)
    sigma = kappa*25*delta*np.dot(w,gradient**2)
    print('flat order energy tension, analytic target',sigma,kappa/(6*delta),flush=True)
    np.testing.assert_allclose(sigma,TENSION,rtol=1e-11)


def measure():
    archive=json.loads(EVIDENCE.read_text())['runs']['mouth_front_pressure_dipole02']
    assert archive['guard_exit']==archive['child_exit']==0 and archive['limits_verified']
    for raw in archive['stdout_lines']:
        record=json.loads(raw)
        base, pair=archive_drive(record)
        for elements, modes, angles in ((16,5,48),(32,8,96)):
            sphere = Sphere(elements,modes,angles)
            single, double = (sphere.prepare(f) for f in (base,pair))
            for K,transport in ((.4,True),(.8,True),(1.6,True),(.4,False)):
                a=sphere.solve(single,K,transport)
                b=sphere.solve(double,K,transport)
                ds = b['s1']-a['s1']
                moment = b['moment']-a['moment']
                pressure = TENSION*VOLUME*ds
                print(json.dumps(dict(grid=record['grid'],d=record['d'],elements=elements,modes=modes,
                    K=K,transport=transport,alone=a,pair=b,delta_slope=ds,
                    moment_slope=moment/(K*VOLUME),pressure_jump_moment=pressure,
                    feed_pressure_moment=record['moment']['pressure'],
                    required_bulk_moment=record['moment']['pressure']-pressure),allow_nan=False),flush=True)
                assert abs(a['s1'])<1e-12 and b['height']>0
                np.testing.assert_allclose(ds,moment/(K*VOLUME),rtol=.008,atol=5e-10)
                np.testing.assert_allclose(K*4*np.pi*b['mean_s'],RATE,rtol=.004)


if __name__=='__main__':
    sign_control()
    measure()
