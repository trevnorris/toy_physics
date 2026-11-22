"""Variational linear-element check with direct energy differences in physical units."""
import json

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import expit


class Strip:
    def __init__(self, n, ambient):
        self.a = 75.
        self.delta = .04
        half = .5
        edge = [-half-ambient*self.delta, -half, half, half+ambient*self.delta]
        pieces = [np.linspace(lo,hi, 1+round((hi-lo)*n/self.delta))
                  for lo,hi in zip(edge,edge[1:])]
        self.x = np.concatenate([part[:-1] for part in pieces]+[pieces[-1][-1:]])
        self.h = np.diff(self.x)
        nodes, weights = np.polynomial.legendre.leggauss(3)
        self.left = (1.-nodes)/2.
        self.right = (1.+nodes)/2.
        self.s = self.x[:-1,None]*self.left + self.x[1:,None]*self.right
        self.quad = .5*self.h[:,None]*weights
        self.fixed = np.searchsorted(self.x,edge)
        self.free = np.ones(len(self.x),dtype=bool)
        self.free[self.fixed] = False
        self.seed = expit((half-self.x)/self.delta)*expit((half+self.x)/self.delta)
        self.seed[self.fixed] = [0., .5, .5, 0.]

    def values(self, c, k, Y):
        l, r = self.left, self.right
        mid = c[:-1,None]*l+c[1:,None]*r
        d = np.diff(c)[:,None]/self.h[:,None]
        jac = 1.+k*self.s
        f = self.a*(mid**2*(1.-mid)**2 + self.delta**2*d**2)
        stress = Y*mid**2*k*self.s
        addition = .5*stress*k*self.s
        work = np.sum(self.quad*(f*self.s+stress*self.s))
        energy = np.sum(self.quad*(jac*f+addition))
        C = jac*self.a*2.*mid*(1.-mid)*(1.-2.*mid)+Y*mid*(k*self.s)**2
        G = jac*2.*self.a*self.delta**2*d
        grad = np.zeros(len(c))
        grad[:-1] += np.sum(self.quad*(l*C - G/self.h[:,None]),axis=1)
        grad[1:] += np.sum(self.quad*(r*C + G/self.h[:,None]),axis=1)
        U = self.quad*(jac*self.a*(2.-12.*mid+12.*mid**2)+Y*(k*self.s)**2)
        V = self.quad*jac*2.*self.a*self.delta**2/self.h[:,None]**2
        diag = np.zeros(len(c))
        diag[:-1] += np.sum(U*l*l+V,axis=1)
        diag[1:] += np.sum(U*r*r+V,axis=1)
        off = np.sum(U*l*r-V,axis=1)
        off[~(self.free[1:]&self.free[:-1])] = 0.
        diag[self.fixed] = 1.
        grad[self.fixed] = 0.
        return energy, work, grad, np.stack((np.r_[0.,off],diag,np.r_[off,0.]))

    def equilibrate(self, k, Y):
        c = self.seed.copy()
        for attempt in range(15):
            E, M, grad, L = self.values(c,k,Y)
            if max(abs(grad)) < 1.e-12:
                return E, M, c
            c -= solve_banded((1,1),L,grad)
        raise ArithmeticError(max(abs(grad)))


if __name__ == '__main__':
    for n in [5, 10, 20]:
        strip = Strip(n,20.)
        for Y in [0., 1.]:
            energies = {}
            virtual = {}
            direct = {}
            kstep = .0001
            for k in [0., -.08, .08, -.04, .04]:
                E, M, c = strip.equilibrate(k,Y)
                energies[k] = E
                virtual[k] = M
                ep, _, _ = strip.equilibrate(k+kstep,Y)
                em, _, _ = strip.equilibrate(k-kstep,Y)
                direct[k] = (ep-em)/(2.*kstep)
            print(json.dumps(dict(elements_per_width=n, Y=Y, energy=energies,
                                  energy_difference_moment=direct, material_moment=virtual)),flush=True)
