"""The saved radius-one volume-preserving neck/flare mode in physical coordinates."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

from mouth_conversion_shape_mode import profile

MANIFEST = Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_conversion_shape.json'


@lru_cache(maxsize=1)
def reference():
    run = json.loads(MANIFEST.read_text())['primary']['exchange']
    if (run['worker_exit'] or run['guard_exit'] or not run['limits_verified']
            or hashlib.sha256(run['literal_stdout'].encode()).hexdigest() != run['stdout_sha256']):
        raise ValueError('The reference meridian measurement did not validate')
    row, = (r for r in json.loads(run['literal_stdout']) if
            (r['waist'], r['exit_to_waist'], r['cosines']) == (1., 1.5, 10))
    pressure, _, length, _, solution = profile(1.5)
    x, weights = leggauss(row['quadrature'])
    z, dz = (1+x)*length/2, weights*length/2
    b, slope = solution(z)
    darea = 4*np.pi*b**2*np.hypot(1.,slope)*dz
    frequency = np.pi*np.arange(1,11)/length
    mean = np.cos(np.outer(frequency,z)) @ darea / sum(darea)
    return row, pressure, length, solution, z, dz, frequency, mean


class ExchangeShape:
    """Radius b and db/dw on the displaced wall; values outside are constant endpoint extensions."""

    def __init__(self, step):
        if abs(step) > .02:
            raise ValueError('Displacement outside the tested small-amplitude mode neighborhood')
        self.step = float(step)
        row, pressure, self.length, solution, z, dz, frequency, mean = reference()
        coeff = np.asarray(row['coefficients'])

        def base(t):
            b, slope = solution(t)
            arc = np.hypot(1.,slope)
            f = coeff @ (np.cos(np.outer(frequency,t))-mean[:,None])
            df = -(coeff*frequency) @ np.sin(np.outer(frequency,t))
            return b, slope, arc, pressure-2/(b*arc), f, df

        b, slope, arc, curvature, f, df = base(z)

        def volume_error(offset):
            h = self.step*f+offset
            jac = 1+h*curvature-self.step*df*slope/arc
            return (4*np.pi/3)*dz @ ((b+h/arc)**3*jac) - row['volume']

        self.offset = brentq(volume_error,-.04,.04,xtol=1e-15) if step else 0.
        h = self.step*f+self.offset
        weight = dz*(b+h/arc)**2*arc*(1+h*curvature)
        offset_rate = -(weight@f)/sum(weight)
        self.volume = row['volume']+volume_error(self.offset)
        t = np.linspace(0.,self.length,2049)
        b, slope, arc, curvature, f, df = base(t)
        h = self.step*f+self.offset
        w = t-h*slope/arc
        w[0],w[-1] = 0.,self.length
        jac = 1+h*curvature-self.step*df*slope/arc
        if np.min(jac)<=0:
            raise ArithmeticError('Displaced meridian folded')
        radii = b+h/arc
        tangent = (slope*(1+h*curvature)+self.step*df/arc)/jac
        tangent[0],tangent[-1] = 0.,0.
        self.neck,self.exit = float(radii[0]),float(radii[-1])
        self._radius = CubicHermiteSpline(w,radii,tangent)
        # Differentiate the normal graph at fixed *physical* w, including
        # the inverse map and the volume-correcting offset.
        self._grid = w
        self._speed = (f+offset_rate)*(1+tangent*slope)/arc

    def radius(self,w):
        clipped = np.clip(w,0.,self.length)
        return self._radius(clipped), self._radius(clipped,1)

    def displacement(self,w):
        return np.interp(w,self._grid,self._speed)


if __name__ == '__main__':
    row, _, _, _, _, _, _, _ = reference()
    x, weights = leggauss(200)
    for step in (-.01,-.005,0.,.005,.01):
        shape = ExchangeShape(step)
        w, dw = shape.length*(1+x)/2, shape.length*weights/2
        b,slope=shape.radius(w)
        vol = (4*np.pi/3)*dw @ b**3
        area = 4*np.pi*dw @ (b*b*np.hypot(1.,slope))
        h = 1e-6
        finite = (ExchangeShape(step+h).radius(w)[0]-ExchangeShape(step-h).radius(w)[0])/(2*h)
        error = np.max(abs(finite-shape.displacement(w)))
        print(json.dumps(dict(step=step,neck=shape.neck,exit=shape.exit,length=shape.length,
                              volume=vol,volume_relative_error=vol/row['volume']-1.,
                              area=area,area_change=area-row['area'],speed_error=error)))
        np.testing.assert_allclose(vol,row['volume'],rtol=2e-12)
        if error>2e-6:
            raise ArithmeticError('Fixed-coordinate shape derivative failed')
