"""Fixed-preparation shape work in a sealed duct with a reversible terminal front.

Each trial shape has a new density, phase and spatially resolved open wave.
Both physical reservoirs, slab thickness, front thickness, action and wall
energy coefficients are held in reference-radius units.
"""
import argparse
import json
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import norm
from scipy.interpolate import CubicHermiteSpline

from mouth_exchange_carrier import BIAS
from mouth_exchange_carrier_density import solve as carrier_solve
from mouth_exchange_geometry import reference
from mouth_exchange_shear import solve as shear_solve


DELTA = .04
SIGMA = 1.
ACTION = 1.
BEND = -(5./6) * SIGMA * DELTA**2
GAUSS = (np.pi**2/3-2) * SIGMA * DELTA**2
COORDINATES = ('log_depth', 'log_exit_radius', 'middle', 'upper_minus_lower')
DIFFERENCE = .003


def emit(**items):
    print(json.dumps(items, allow_nan=False), flush=True)


@lru_cache(maxsize=1)
def reference_spline():
    _, pressure, length, solution, *_ = reference()
    w = np.linspace(0, length, 2001)
    r, slope = solution(w)
    return pressure, length, CubicHermiteSpline(w, r, slope)


class Meridian:
    """Free depth, exit opening and two interior widths at a prescribed aperture."""
    def __init__(self, entrance, amplitudes, extra=0.):
        self.neck = entrance
        self.params = np.asarray(amplitudes)
        self.extra = extra
        self.pressure, self.unit_length, self.reference = reference_spline()
        self.length = entrance*self.unit_length*np.exp(self.params[0])
        self.exit = 1.5*entrance*np.exp(self.params[1])

    def profile(self, t):
        b = self.reference(t*self.unit_length)
        slope = self.reference(t*self.unit_length, 1)
        first = slope*self.unit_length
        second = (2*(1+slope*slope)/b-self.pressure*(1+slope*slope)**1.5)*self.unit_length**2
        e = self.exit/self.neck-1.5
        m, a = self.params[2:]
        # All added modes vanish at the entrance. The last three also
        # vanish at the exit. The unused quartic probes generalization.
        t1 = t*(1-t)
        b = b + e*t + 4*m*t1 + 4*a*t1*(2*t-1) + 16*self.extra*t1**2
        first = first + e + 4*m*(1-2*t) + 4*a*(-6*t*t+6*t-1) + 32*self.extra*t1*(1-2*t)
        second = second - 8*m + 4*a*(6-12*t) + 32*self.extra*(1-6*t+6*t*t)
        return self.neck*b, self.neck*first/self.length, self.neck*second/self.length**2

    def radius(self, w):
        r, slope, _ = self.profile(np.clip(w/self.length, 0, 1))
        return r, np.where((w<0)|(w>self.length),0.,slope)

    def wall(self):
        x, weights = leggauss(60)
        wgt = weights*self.length/2
        r, slope, second = self.profile((x+1)/2)
        arc = np.hypot(1,slope)
        turn = -second/arc**3
        hoop = 1/(r*arc)
        trace = turn+2*hoop
        second_invariant = 2*turn*hoop + hoop**2
        darea = 4*np.pi*r**2*arc*wgt
        area = sum(darea)
        curvature_work = darea@(BEND*trace**2/2 + GAUSS*second_invariant)
        r_all, s_all, s2_all = self.profile(np.linspace(0,1,201))
        max_curv = max(np.max(abs(s2_all)/(1+s_all**2)**1.5),
                       np.max(1/(r_all*np.hypot(1,s_all))))
        cap = 4*np.pi*self.exit**3/3
        return dict(area=area, cap=cap, curvature_energy=curvature_work,
                    minimum_radius=min(r_all), maximum_abs_principal_delta=max_curv*DELTA,
                    end_slopes=[s_all[0],s_all[-1]])


class Experiment:
    def __init__(self, radius, mesh, port):
        self.radius, self.mesh, self.port = radius, mesh, port
        self.previous_mode = (None, None)
        self.axial_cells = int(np.ceil(2.77*radius*mesh))

    @lru_cache(maxsize=25)
    def fields(self, coordinates, extra=0.):
        shape = Meridian(self.radius, coordinates, extra)
        wall = shape.wall()
        if wall['minimum_radius'] <= 0 or shape.length >= 22:
            raise ValueError('Shape outside the resolved positive-radius finite reservoir')
        carrier = carrier_solve(0, self.mesh, BIAS, lambda step:shape, sealed_side=True)
        prior, seed = self.previous_mode
        wave, vec, matrix = shear_solve(shape, self.mesh, self.port, 0., prior=prior,
                                        axial_cells=self.axial_cells, shape_work=False,
                                        incoming_vector=seed)
        if not np.isclose(carrier['current_by_cut']['cap_layer'],carrier['current'],rtol=2e-7):
            raise ArithmeticError('Sealed duct port current differs from terminal disk current')
        total = SIGMA*(wall['area']+wall['cap']) + wall['curvature_energy'] + carrier['energy'] + ACTION*wave['frequency']
        metrics = dict(depth=shape.length,exit_radius=shape.exit,**wall,carrier_energy=carrier['energy'],
                       current=carrier['current'], wave_frequency=wave['frequency'], wave_rate=2*wave['amplitude_loss'],
                       feed_energy_flux_over_stored=wave['outward_power_over_interior_energy'],
                       wave_pocket_fraction=wave['confined_kinetic_fraction'],
                       clamped_overlap=wave['clamped_seed_overlap'], total=total)
        return metrics,vec,matrix

    def signed(self, x, h, extra=False):
        base,v,_ = self.fields(tuple(x))
        history = []
        indices = range(5) if extra else range(4)
        for i in indices:
            d = np.zeros(4)
            if i < 4:
                d[i] = h
            samples=[]
            for sign in (-1.,1.):
                row,vs,mass = self.fields(tuple(x+sign*d),sign*h if i==4 else 0.)
                overlap=abs(np.vdot(v,mass@vs))**2/np.vdot(v,mass@v).real
                if overlap<.9:
                    raise ArithmeticError(f'No consistent wave under signed variation, overlap {overlap}')
                samples.append(row)
            minus,plus=samples
            quantities=('total','area','cap','curvature_energy','carrier_energy','wave_frequency')
            slopes={k:(plus[k]-minus[k])/(2*h) for k in quantities}
            history.append(slopes)
        return history


def run(a):
    exp=Experiment(a.radius,a.mesh,a.port)
    x=np.array(a.start)
    emit(test='preparation', radius=a.radius, mesh=a.mesh, port=a.port, axial_cells=exp.axial_cells,
         coordinates=COORDINATES, physical_action=ACTION, sigma=SIGMA,delta=DELTA,B=BEND,m2=GAUSS,
         pressure_jump=0.,phase_bias=BIAS, side_normal_current=0.,terminal_sheet_length=DELTA/6)
    if a.mode=='point':
        full=exp.signed(x,DIFFERENCE,extra=True)
        half=exp.signed(x,DIFFERENCE/2,extra=True)
        emit(test='point', coordinates=x.tolist(), measurements=exp.fields(tuple(x))[0],
             signed=full, half_signed=half)
        return
    for step in range(50):
        base,v,_=exp.fields(tuple(x))
        emit(test='state', index=step, coordinates=x.tolist(),**base)
        if a.mode == 'continuation':
            exp.previous_mode = complex(base['wave_frequency'], -base['wave_rate']/2), v
        # A pocket of less than half the kinetic norm has ceased to be the
        # localized branch targeted by the archived fixed-action preparation.
        # Resolve the boundary from both sides instead of calling it equilibrium.
        stop_reason = None
        if base['wave_pocket_fraction']<.5:
            stop_reason='delocalized'
        elif a.mode=='relax' and base['maximum_abs_principal_delta']>.10:
            stop_reason='thin_front'
        elif a.mode in ('extension','continuation') and base['minimum_radius'] <= DELTA*np.log(9):
            # The planar order profile requires delta*log(9) between the
            # half-order dividing surface and material with 90 percent order.
            stop_reason='waist comparable to planar half-to-90-percent width'
        if stop_reason:
            emit(test='limit', reason=stop_reason,
                 coordinates=x.tolist(), measurements=base)
            return
        slopes=exp.signed(x,DIFFERENCE)
        g=np.array([row['total'] for row in slopes])
        emit(test='force', index=step, derivatives=slopes)
        if norm(g)<.01:
            emit(test='stationary_candidate', coordinates=x.tolist(), measurements=base)
            return
        descent=-g/norm(g);distance=.13
        while distance >= .13/128:
            candidate=x+distance*descent
            trial_meridian=Meridian(a.radius,candidate)
            if min(trial_meridian.profile(np.linspace(0,1,301))[0])<=0:
                distance/=2
                continue
            cm,cv,m=exp.fields(tuple(candidate))
            overlap=abs(np.vdot(v,m@cv))**2/np.vdot(v,m@v).real
            if overlap>.90 and cm['total']<=base['total']-.05*distance*norm(g):
                emit(test='accepted', distance=distance, overlap=overlap)
                x=candidate
                break
            distance/=2
        else:
            emit(test='limit', reason='local descent or mode correspondence unresolved', coordinates=x.tolist())
            return
    emit(test='limit', reason='finite descent inconclusive', coordinates=x.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--port',type=float,default=4.2)
    p.add_argument('--mode',choices=('relax','point','extension','continuation'),default='relax')
    p.add_argument('--start',type=float,nargs=4,default=(0.,0.,0.,0.))
    run(p.parse_args())
