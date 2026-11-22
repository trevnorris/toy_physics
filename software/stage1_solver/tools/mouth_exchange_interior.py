"""Independent redistribution between the two halves of the side wall.

The feed aperture, terminal radius, axial length and enclosed volume are
fixed. With those constraints this measures distributed work on the flare.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from mouth_exchange_geometry import reference


class InteriorShape:
    def __init__(self, step):
        self.step = float(step)
        if abs(step) > .02:
            raise ValueError('Interior wall test is a local variation')
        row, _, self.length, self._solution, *_ = reference()
        self.neck, self.exit = 1., 1.5
        x, q = leggauss(80)
        w = self.length*np.r_[x+1,x+3]/4
        dw = self.length*np.tile(q,2)/4
        b = self._solution(w)[0]
        low, high, _, _ = self._bumps(w)
        def change(second):
            db = step*low+second*high
            return dw @ (db*(3*b*b+db*(3*b+db)))
        self.second = brentq(change,-.06,.06,xtol=1e-15) if step else 0.
        b = b+step*low+self.second*high
        self._slope = -(dw@(b*b*low))/(dw@(b*b*high))
        self.volume = (4*np.pi/3)*dw@b**3
        np.testing.assert_allclose(self.volume,row['volume'],rtol=2e-12)

    def _bumps(self,w):
        x = np.clip(np.asarray(w)/self.length,0.,1.)
        theta = 2*np.pi*x
        y = np.sin(theta)
        b = y**4
        db = (8*np.pi/self.length)*y**3*np.cos(theta)
        first = x < .5
        return b*first,b*(~first),db*first,db*(~first)

    def radius(self,w):
        x = np.clip(np.asarray(w,dtype=float),0.,self.length)
        b, slope = self._solution(x.ravel()).reshape((2,)+x.shape)
        low, high, dl, dh = self._bumps(x)
        return b+self.step*low+self.second*high, slope+self.step*dl+self.second*dh

    def displacement(self,w):
        low, high, _, _ = self._bumps(w)
        return low+self._slope*high


if __name__=='__main__':
    x, q = leggauss(100)
    for step in (-.005,0.,.005):
        shape = InteriorShape(step)
        w = shape.length*np.r_[x+1,x+3]/4
        dw = shape.length*np.tile(q,2)/4
        b,slope = shape.radius(w)
        area = 4*np.pi*dw@(b*b*np.hypot(1.,slope))
        speed = (InteriorShape(step+1e-6).radius(w)[0]-InteriorShape(step-1e-6).radius(w)[0])/2e-6
        error = np.max(abs(speed-shape.displacement(w)))
        print(json.dumps(dict(step=step,second=shape.second,area=area,volume=shape.volume,
                              speed_error=error,neck=float(shape.radius(0.)[0]),
                              exit=float(shape.radius(shape.length)[0]))))
        np.testing.assert_allclose(shape.radius([0.,shape.length])[0],[1.,1.5],atol=1e-12)
        if error>2e-7:
            raise ArithmeticError('Wrong shape velocity on independently held-end deformation')
