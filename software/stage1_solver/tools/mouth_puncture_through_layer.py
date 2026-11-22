"""Steady phase transfer across a smooth, curved order transition on a stream path.

The projected inner opening is empty in the sheet FEM. Here a path feeding a
drain has to change order. Compare a barrier sensitive only to local order to
one sensitive to an independently stipulated spanning crystalline structure.
Area is per unit midsurface area, n0=delta_B=K0=1, no energy tilt between phases.
"""
import numpy as np
from scipy.integrate import quad

from mouth_flat_phase_control import saved_bath
from mouth_free_front_wave_head import emit


def measure(strength, spanning):
    radius, length = 100.,32.
    phase = saved_bath()/10
    def density(s):
        area = (1+s/radius)**3
        k = 1. - strength*spanning/np.cosh(s/2)**2
        return 1./(area*k)
    def integrate(gap):
        return sum(quad(density,lo,hi,epsabs=2e-10,epsrel=1e-11)[0]
                   for lo,hi in [(-length/2,-gap),(gap,length/2)])
    if strength*spanning==1.:
        for gap in (.2,.1,.05,.025):
            resistance = integrate(gap)
            emit(test='critical_layer_cutoff', omitted_half_width=gap, retained_resistance=resistance,
                 upper_bound_on_positive_current_per_midarea=phase/resistance,
                 leading_resistance_times_cutoff=resistance*gap)
        emit(test='smooth_separating_layer',strength=strength,spanning=spanning,
             phase=phase, finite_energy_stationary_current=0., limit='divergent normal resistance',
             continuous_classical_phase_solution_for_nonzero_drop=False)
        return
    resistance=integrate(0.)
    current=phase/resistance
    # The independent phase-energy quadrature uses the flux solution.
    energy=quad(lambda s: .5*current**2*density(s),-length/2,length/2,
                points=[0.],epsabs=2e-14)[0]
    emit(test='smooth_separating_layer',strength=strength,spanning=spanning,
         phase=phase, finite_energy_stationary_current=current, resistance=resistance,
         phase_energy=energy, twice_energy_minus_flux_work=2*energy-current*phase)
    if abs(2*energy-current*phase)>2e-10:
        raise ArithmeticError('Phase energy and fixed port phases disagree')


if __name__=='__main__':
    for strength,spanning in ((0.,1),(.9,1),(.9999,1),(1.,1),(1.,0)):
        measure(strength,spanning)
