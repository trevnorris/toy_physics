"""Carry the inlet through the same prescribed deeper bend before reading the deep rows' force."""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator

from mouth_carrier_lateral_inertia import taper
from mouth_empty_bend_stiffness import INTERVALS
from mouth_lining_symmetric_passage import PITCH, ROW_SPACING, SUPPORT, cubic_trace, read_inputs, row_potential
from neutral_patch_sites import emit


class Strip:
    def __init__(self, mesh, stride, grid):
        trace, intake = cubic_trace(mesh)
        patch, _, _ = read_inputs(mesh)
        assert abs(intake-.25) < 1.e-8
        depth = np.unique(np.r_[np.linspace(1.175,1.65,8193),1.55,
                                trace.edges[(trace.edges > 1.175) & (trace.edges < 1.65)]])
        clock = cumulative_trapezoid(1/trace.speed(depth), depth, initial=0.)
        self.pitch = patch['pitch']*PITCH
        arc = cumulative_trapezoid(np.hypot(1.,trace.geometry.radius(depth)[1]),depth,initial=0.)/self.pitch
        self.xrows = ROW_SPACING*np.arange(int(np.floor(arc[-1]/ROW_SPACING))+1)
        self.depths = PchipInterpolator(arc, depth)(self.xrows)
        self.deep = self.depths >= 1.55
        self.shape = 1-taper(self.depths, *INTERVALS[1])[0]
        self.time = np.linspace(0,clock[-1],len(self.xrows)*stride+1)
        self.dt = self.time[1]
        self.gate = int(np.argmin(abs(self.time-np.interp(1.55,depth,clock))))
        self.x = PchipInterpolator(clock,arc)(self.time)
        self.chi = trace.sample(trace.chi,PchipInterpolator(clock,depth)(self.time))[0]
        self.mass_scale = 2*patch['shared_scale']/self.pitch**2
        self.ny, self.nx = grid, grid//2
        u = np.arange(self.ny+1)/self.ny
        x = np.arange(SUPPORT*self.nx+1)[:,None]/self.nx
        g = row_potential(x,u[None,:],patch['pitch'],True)[1]
        # Periodize at the image-truncation roundoff scale and window the axial end.
        np.testing.assert_allclose(g,-g[:,::-1],atol=5e-13)
        g = .5*(g-g[:,::-1])
        g[:,-1] = g[:,0]
        window = np.clip(SUPPORT-x,0.,1.)
        self.spline = CubicSpline(u,g*window**2*(3-2*window),axis=1,bc_type='periodic')
        emit(test='fixture', intake=intake, native_pitch=patch['pitch'], pitch=self.pitch,
             rows=len(self.xrows), deep_rows=int(self.deep.sum()), first_deep_depth=float(self.depths[self.deep][0]),
             mouth_depths=[1.175,1.55,1.65], transit=self.time[-1], gauge_time=self.time[self.gate],
             step=self.dt, shape_range=[self.shape[0],self.shape[-1]], grid=grid,stride=stride,mesh=mesh)

    def tables(self, q):
        size = self.ny
        result = np.zeros((len(self.time),size+1))
        deep = np.zeros_like(result)
        row_modes = np.arange(size)/size
        for j,center in enumerate(self.xrows):
            relative = (row_modes-(j%2)/2-q*self.shape[j])%1
            row_table = self.spline(relative)
            distance = self.nx*abs(self.x-center)
            active = distance < SUPPORT*self.nx
            index = distance[active].astype(int)
            frac = (distance[active]-index)[:,None]
            g = -self.mass_scale*self.chi[active,None]*(
                (1-frac)*row_table[index]+frac*row_table[index+1])
            result[active,:size] += g
            if self.deep[j]:
                deep[active,:size] += g
        result[:,-1]=result[:,0]
        deep[:,-1]=deep[:,0]
        return result, deep


def interpolate(row, index, fraction):
    return row[index]+fraction*(row[index+1]-row[index])


def transit(strip, q, phases):
    total, deep = strip.tables(q)
    y = (np.arange(phases)+.5)/phases
    velocity = np.zeros_like(y)
    def sample(k):
        u = strip.ny*(y-np.floor(y))
        j = np.minimum(u.astype(np.intp),strip.ny-1)
        return interpolate(total[k],j,u-j), interpolate(deep[k],j,u-j)
    a, row_a = sample(0)
    dt = strip.dt
    after, deep_from_entry, deep_after = 0.,0.,0.
    for k in range(1,len(strip.time)):
        y += dt*velocity + .5*dt*dt*a
        new, new_row = sample(k)
        velocity += .5*dt*(a+new)
        dp = .5*dt*(row_a.mean()+new_row.mean())
        deep_from_entry += dp
        if k == strip.gate:
            gauge_v = velocity.mean()
            gauge_rms = np.sqrt(np.mean(velocity**2))
        if k > strip.gate:
            after += .5*dt*(a.mean()+new.mean())
            deep_after += dp
        a, row_a = new, new_row
    np.testing.assert_allclose(after,velocity.mean()-gauge_v,atol=3.e-12,rtol=1.e-8)
    return dict(gauge_mean_v=gauge_v,gauge_rms_v=gauge_rms,exit_mean_v=velocity.mean(),
                all_rows_after_gauge=after,deep_rows_full_path=deep_from_entry,
                deep_rows_after_gauge=deep_after)


def run(args):
    strip=Strip(args.mesh,args.stride,args.grid)
    result={}
    for q in (0.,args.step,-args.step,args.step/2,-args.step/2):
        result[q]=transit(strip,q,args.phases)
        emit(test='passage',q_pitches=q,phases=args.phases,**result[q])
    for q in (args.step,args.step/2):
        for key in ('gauge_mean_v','exit_mean_v','all_rows_after_gauge','deep_rows_full_path','deep_rows_after_gauge'):
            plus,minus,base=result[q][key],result[-q][key],result[0.][key]
            emit(test='signed_slope',field=key,step_pitches=q,
                 odd=(plus-minus)/(2*q),even=((plus-base)+(minus-base))/(2*q))
    # Both momentum and displacement use the same measured pitch. Slopes are
    # holder force per incoming tag mass flux per *bend-amplitude* in inverse time.


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,required=True,choices=(26,34))
    parser.add_argument('--stride',type=int,required=True)
    parser.add_argument('--grid',type=int,required=True)
    parser.add_argument('--phases',type=int,required=True)
    parser.add_argument('--step',type=float,default=.001)
    args=parser.parse_args()
    if min(args.phases,args.grid,args.stride)<2 or args.phases%2 or args.grid%2 or not 0<args.step<.1:
        parser.error('Expected positive resolution, even phase and table counts, and a weak bend step')
    run(args)
