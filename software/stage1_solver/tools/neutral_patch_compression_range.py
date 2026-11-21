"""Check the size of the EOS mismatch over densities present on the archived deep wall."""
import numpy as np
from numpy.polynomial import Polynomial

from mouth_middle_tag_grip import cubic_trace
from neutral_patch_carrier_compression import CALIBRATION, Cell, carrier_pressure
from neutral_patch_sites import emit


def wall_density_range(trace, start, end):
    samples = []
    # The wall field is complex Q2 on each element; density is its squared modulus.
    for lo,hi,values in zip(trace.edges[:-1],trace.edges[1:],trace.psi):
        a,b = max(lo,start),min(hi,end)
        if b<a:
            continue
        left,right = (2*(np.array([a,b])-lo)/(hi-lo)-1)
        f = Polynomial([values[1],(values[2]-values[0])/2,(values[2]+values[0])/2-values[1]])
        density = Polynomial(f.coef.real)**2 + Polynomial(f.coef.imag)**2
        roots = density.deriv().roots()
        for x in [left,right,*[r.real for r in roots if abs(r.imag)<1e-10 and left<r.real<right]]:
            w = lo+(1+x)*(hi-lo)/2
            check = abs(trace.sample(trace.psi,np.array([w]))[0][0])**2
            np.testing.assert_allclose(density(x),check,atol=2e-13,rtol=0)
            samples.append((float(density(x)),float(w)))
    return min(samples),max(samples)


def run():
    cell = Cell()
    for mesh in (26,34):
        trace,intake = cubic_trace(mesh)
        low,high = wall_density_range(trace,1.175,2.)
        emit(test='wall_density',mesh=mesh,intake=intake,depth_interval=[1.175,2.],
             minimum_density_and_depth=low,maximum_density_and_depth=high)
        for gap_in_brane in (False,True):
            def raw(n):
                return cell.measure(n,gap_in_brane,True)[0]['raw_pressure']
            reference = raw(1.)
            scale = carrier_pressure(CALIBRATION)/(raw(CALIBRATION)-reference)
            def derivative(h):
                return (raw(np.exp(h))-raw(np.exp(-h)))/(2*h)
            # Separate, linear-only reading about the undriven resting state.
            slope = (4*derivative(.0005)-derivative(.001))/3
            rows = []
            for n in np.linspace(low[0],high[0],9):
                target = carrier_pressure(n)
                change = raw(n)-reference
                rows.append([n,target,scale*change/target-1,change/(slope*target)-1])
            emit(test='operating_range',mesh=mesh,stack='brane' if gap_in_brane else 'w',
                 zero_strain_scale=1/slope,one_percent_scale=scale,
                 columns=['density','EOS_pressure_increment','one_percent_fit_fractional_error','zero_strain_fit_fractional_error'],
                 rows=rows)


if __name__=='__main__':
    run()
