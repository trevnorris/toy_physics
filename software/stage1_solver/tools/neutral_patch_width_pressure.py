"""Squeeze a stripe of a settled neutral bilayer and measure lateral exchange.

The pressure is conjugate to the actual volume between the two triangulated
faces in a fixed laboratory window. A second, independent bias is conjugate to
the weighted number of sites in the window. All constituents can move.
"""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cho_factor, cho_solve, eigh

from neutral_patch_sites import emit
from neutral_patch_wrinkle import WrinkleCell

ROOT = Path('_scratch/s11c')
STATES = {1: 'wrinkle24_p1_01', 2: 'wrinkle24_p2_01'}
SIZES = (1., 2., 4.)
STEP = .001


def triangles(side):
    i,j = np.indices((side,side)).reshape(2,-1)
    def index(x,y):
        return (x%side)*side + y%side
    a = index(i,j)
    b = index(i+1,j)
    c = index(i+j%2,j+1)
    d = index(i+j%2-1,j+1)
    plane = np.vstack((np.column_stack([a,b,c]), np.column_stack([a,c,d])))
    return np.vstack([plane,plane+side**2])


class Stripe:
    def __init__(self, cell, size, center, degree=5):
        self.cell, self.size, self.center = cell,size,center
        self.ids = triangles(cell.side)
        self.sign = np.repeat([-1.,1.], len(self.ids)//2)
        x, wg = leggauss(degree)
        u,v = (z.ravel() for z in np.meshgrid((1+x)/2,(1+x)/2,indexing='ij'))
        self.basis = np.column_stack([(1-u)*(1-v),u,(1-u)*v])
        self.weights = np.outer(wg,wg).ravel()*(1-u)/4
        self.area = cell.area if size is None else (
            cell.box[0]*np.sqrt(2*np.pi)*size/self.peak())
    def peak(self):
        period = self.cell.box[1]
        return np.exp(-.5*(np.arange(-2,3)*period/self.size)**2).sum()
    def weight(self, y):
        if self.size is None:
            return np.ones_like(y), np.zeros_like(y)
        period = self.cell.box[1]
        d = y-self.center
        d = d-period*np.rint(d/period)
        weight, slope = np.zeros_like(d),np.zeros_like(d)
        for image in range(-2,3):
            distance=d+image*period
            term=np.exp(-.5*(distance/self.size)**2)
            weight+=term
            slope-=term*distance/self.size**2
        return weight/self.peak(),slope/self.peak()

    def read(self, pos):
        # Unwrap each material triangle and integrate both faces at the same xy.
        tri = pos[self.ids].copy()
        for k in (1,2):
            dx = tri[:,k,:2]-tri[:,0,:2]
            tri[:,k,:2] -= self.cell.box*np.rint(dx/self.cell.box)
        a,b = tri[:,1,:2]-tri[:,0,:2],tri[:,2,:2]-tri[:,0,:2]
        det = a[:,0]*b[:,1]-a[:,1]*b[:,0]
        if np.min(det)<=0:
            raise ValueError('An original material face no longer defines a height graph')
        y, z = tri[:,:,1]@self.basis.T, tri[:,:,2]@self.basis.T
        g, dg = self.weight(y)
        wg = g*self.weights
        val = (wg*z).sum(axis=1)
        faces = np.sum((det*val).reshape(2,-1),axis=1)/self.area
        areas = np.sum((det*wg.sum(axis=1)).reshape(2,-1),axis=1)
        # Derivative of the signed, pressure-weighted volume, including slopes.
        ddet = np.empty((len(det),3,2))
        ddet[:,1,0],ddet[:,1,1] = b[:,1],-b[:,0]
        ddet[:,2,0],ddet[:,2,1] = -a[:,1],a[:,0]
        ddet[:,0] = -(ddet[:,1]+ddet[:,2])
        local = np.empty_like(tri)
        local[:,:,:2] = ddet*val[:,None,None]
        local[:,:,1] += det[:,None]*((dg*z*self.weights)@self.basis)
        local[:,:,2] = det[:,None]*(wg@self.basis)
        vgrad = np.stack([np.bincount(self.ids.ravel(),weights=(self.sign[:,None]*local[:,:,d]).ravel(),
                                      minlength=len(pos)) for d in range(3)],axis=1)
        w, derivative = self.weight(pos[:,1])
        ngrad = np.zeros_like(pos)
        ngrad[:,1]=derivative
        volume = self.area*(faces[1]-faces[0])
        values = dict(volume=volume, width=volume/self.area, midplane=faces.mean(),
                      number=w.sum(), density=w.sum()/volume)
        if volume<=0 or not np.allclose(areas,self.area,atol=2e-8,rtol=1e-10):
            raise ArithmeticError(f'Invalid width or quadrature area: {volume}, {areas-self.area}')
        return values, vgrad, ngrad


def tangent_solver(cell, pos):
    h, _, _ = cell.evaluate(pos,tangent=True)
    for k in range(3):
        h[k::3,k::3] += 2000./cell.n
    eig,vec=eigh(h,subset_by_index=(0,1))
    if eig[0]<=0:
        raise ValueError(f'Unstable preparation: {eig}')
    fac=cho_factor(h,overwrite_a=True)
    def solve(g):
        centered = g-g.mean(axis=0)
        return cho_solve(fac,centered.ravel()).reshape(pos.shape)
    return solve,eig,vec[:,0]


def experiment(power, center_fraction):
    cell = WrinkleCell(power)
    path = ROOT/STATES[power]/'state.npz'
    with np.load(path) as state:
        pos=state['position'].copy()
        np.testing.assert_allclose(cell.box,state['box'])
        assert int(state['power'])==power
    base_force=cell.evaluate(pos)[1]
    solve,eig,phase=tangent_solver(cell,pos)
    center=cell.box[1]*center_fraction
    sensor=Stripe(cell,1.,center)
    ref,_,_=sensor.read(pos)
    emit(test='preparation', power=power, center_fraction=center_fraction, sites=cell.n,
         periods=cell.box.tolist(), reference=ref, residual=abs(base_force).max(),
         state=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
         least_stiffness=eig.tolist())
    for size in (*SIZES,None):
        stripe=Stripe(cell,size,center)
        zero,vg,ng=stripe.read(pos)
        dv,dn=solve(vg),solve(ng)
        cvv,cvn,cnn=np.sum(vg*dv),np.sum(ng*dv),np.sum(ng*dn)
        assert np.allclose(cvn,np.sum(vg*dn),atol=1e-10)
        check=Stripe(cell,size,center,degree=7)
        alt,av,an=check.read(pos)
        emit(test='response', size=size, loaded_area=stripe.area, base=zero,
             volume_pressure=cvv, number_pressure=cvn, number_affinity=cnn,
             width_pressure=cvv/stripe.area,
             width_pressure_at_fixed_weighted_number=(cvv-cvn**2/cnn)/stripe.area if cnn>0 else cvv/stripe.area,
             soft_fraction_volume=(phase@vg.ravel())**2/(eig[0]*cvv),
             net_force=vg.sum(axis=0).tolist(),force_quadrature_difference=abs(av-vg).max())
        np.testing.assert_allclose(av,vg,atol=2e-8)
        np.testing.assert_allclose(alt['width'],zero['width'],atol=1e-10)
        np.testing.assert_allclose(vg.sum(axis=0),0,atol=3e-8)
        if size is None: # In this control the whole periodic cell is closed to lateral exchange.
            np.testing.assert_allclose(ng,0)
        loads=[('pressure',vg,dv)]
        if size is not None:
            loads.append(('affinity',ng,dn))
        for drive,force,dx in loads:
            susceptibility=np.sum(force*dx)
            phase_fraction=(phase@force.ravel())**2/(eig[0]*susceptibility)
            emit(test='drive_check', size=size, drive=drive, soft_compliance_fraction=phase_fraction,
                 peak_linear_displacement=np.max(abs(STEP*dx)))
            if phase_fraction>.01 or np.max(abs(STEP*dx))>.01:
                raise ArithmeticError('Pressure/particle drive strongly excites a soft translation')
            for load in (STEP,-STEP,STEP/2,-STEP/2):
                moved=pos+load*dx
                # Quasistatic finite-load solutions, independently checked against
                # full nonlinear forces. Preserve only the mean translation gauge.
                for _ in range(12):
                    vals,v,n=stripe.read(moved)
                    energy,grad,_=cell.evaluate(moved)
                    residual=grad-load*(v if drive=='pressure' else n)
                    if np.max(abs(residual))<1.5e-9:
                        break
                    moved-=solve(residual)
                current,_,_=sensor.read(moved)
                changes={key:np.log(current[key]/ref[key]) for key in ('width','number','density')}
                vals,v,n=stripe.read(moved)
                full_check=cell.evaluate(moved)[1]-load*(v if drive=='pressure' else n)
                emit(test='load', power=power,size=size,drive=drive,load=load,center=current,
                     center_log_changes=changes, loaded_region=vals,
                     center_midplane_change=current['midplane']-ref['midplane'],
                     peak_displacement=abs(moved-pos).max(),residual=abs(full_check).max())
                if abs(full_check).max()>5e-9:
                    raise ArithmeticError('Nonlinear response unresolved')
                np.testing.assert_allclose(changes['density'],changes['number']-changes['width'],atol=4e-15)
                exact=(vals['volume']-zero['volume']) if drive=='pressure' else (vals['number']-zero['number'])
                expected=cvv if drive=='pressure' else cnn
                np.testing.assert_allclose(exact/load,expected,rtol=.009,atol=2e-6)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('power',type=int,choices=(1,2))
    p.add_argument('--center',type=float,default=.5)
    args=p.parse_args()
    experiment(args.power,args.center)
