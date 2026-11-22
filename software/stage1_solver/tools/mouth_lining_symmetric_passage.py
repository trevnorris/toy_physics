"""Resolve the reflection-symmetric cold inlet of the mobile transverse row strip.

Opposite relative phases have opposite accelerations on *each* material row.
One representative from each reflected pair therefore evolves in the unshifted
staggered row pattern even with every material row free. Uniform inlet phase is
integrated by an equal-flux midpoint rule. Speeds below are in pitches per time.
"""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_mobile_lining import ROW_SPACING, SUPPORT, row_potential
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import emit


def observe(y, v, x, n_rows, entry, depth, step, reference='interacting'):
    j = min(n_rows-1, int(np.rint(x/ROW_SPACING)))
    offset = y - (j%2)/2
    # Signed cosine moments of the full distribution; all sine moments vanish
    # by reflection about the actual row. The two members have opposite speeds.
    moments = [np.mean(np.cos(2*np.pi*k*offset)) for k in (1,2,4)]
    emit(test='crossing_reading', reference=reference, depth=depth, step=step, row=j,
         relative_phase_cosines=moments, fraction_within_eighth_pitch=np.mean(abs(offset-np.rint(offset))<=.125),
         rms_relative_speed=np.sqrt(np.mean(v**2)), rms_unwrapped_shift_pitches=np.sqrt(np.mean((y-entry)**2)))


def run(args):
    tr, intake = cubic_trace(args.mesh)
    patch, _, _ = read_inputs(args.mesh)
    pitch_native=patch['pitch']
    pitch=pitch_native*PITCH
    w=np.unique(np.r_[np.linspace(1.175,1.65,8193),tr.edges[(tr.edges>1.175)&(tr.edges<1.65)]])
    clock=cumulative_trapezoid(1/tr.speed(w),w,initial=0.)
    arc=cumulative_trapezoid(np.hypot(1,tr.geometry.radius(w)[1]),w,initial=0.)/pitch
    n=int(np.floor(arc[-1]/ROW_SPACING))+1
    steps=n*args.stride
    t=np.linspace(0,clock[-1],steps+1)
    depth=PchipInterpolator(clock,w)(t)
    x=PchipInterpolator(clock,arc)(t)
    order=tr.sample(tr.chi,depth)[0]
    np.testing.assert_allclose(clock[-1],tr.transit(w[0],w[-1],12),rtol=1e-6)
    ny, nx = args.grid, args.grid//2
    ys = np.arange(ny+1)[None,:]/ny
    xs = np.arange(SUPPORT*nx+1)[:,None]/nx
    row_force=row_potential(xs,ys,pitch_native,True)[1]
    np.testing.assert_allclose(row_force,-row_force[:,::-1],atol=3e-13)
    window=np.clip(SUPPORT-xs,0,1)
    row_force*=window**2*(3-2*window)
    # Sum the contributions from the two staggered row types before any tag
    # interpolation. For the invariant branch this equals applying each row force.
    force=np.zeros((len(t),ny))
    for j in range(n):
        dx=abs(x-j*ROW_SPACING)*nx
        mask=dx<SUPPORT*nx
        i=dx[mask].astype(int)
        f=(dx[mask]-i)[:,None]
        row=(1-f)*row_force[i,:ny]+f*row_force[i+1,:ny]
        force[mask]+=np.roll(row,ny//2 if j%2 else 0,axis=1)
    force*=-order[:,None]*2*patch['shared_scale']/pitch**2
    force=np.concatenate((force,force[:,:1]),axis=1)
    np.testing.assert_allclose(force,-force[:,::-1],atol=4e-12)
    if args.phases%2 or ny%2: raise ValueError('Use a reflection-compatible grid')
    entry=(np.arange(args.phases//2)+.5)/args.phases
    y=entry.copy()
    v=np.zeros_like(y)
    probes={int(np.argmin(abs(depth-d))) for d in (1.175,1.35,1.50,1.54,1.55,1.56,1.65)}
    gate=int(np.argmin(abs(depth-1.55)))
    dt=t[1]
    emit(test='symmetric_setup', mesh=args.mesh, intake=intake, full_phase_count=args.phases,
         stride=args.stride, grid=ny, integrated_rows=n, transit=t[-1], pitch=pitch,
         gauge_depth=depth[gate], dt=dt)
    def acceleration(k):
        u=(y-np.floor(y))*ny
        i=np.minimum(u.astype(np.intp),ny-1)
        return force[k,i]+(u-i)*(force[k,i+1]-force[k,i])
    a=acceleration(0)
    for k in range(steps+1):
        if k in probes: observe(y,v,x[k],n,entry,depth[k],k)
        if k==gate: initial_y, initial_v = y.copy(),v.copy()
        if k==steps: break
        y+=dt*v+.5*dt**2*a
        a_next=acceleration(k+1)
        v+=.5*dt*(a+a_next)
        a=a_next
    ballistic=initial_y+(t[-1]-t[gate])*initial_v
    observe(ballistic,initial_v,x[-1],n,entry,depth[-1],steps,'free_after_gauge')
    emit(test='downstream_difference', rms_speed_change=np.sqrt(np.mean((v-initial_v)**2)),
         rms_separation_from_ballistic_pitches=np.sqrt(np.mean((y-ballistic)**2)),
         incoming_rms_speed=np.sqrt(np.mean(initial_v**2)), time_after_gauge=t[-1]-t[gate])
    if not np.isfinite(y).all() or not np.isfinite(v).all(): raise ArithmeticError('Nonfinite transit')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('mesh','stride','phases','grid'): parser.add_argument('--'+name,required=True,type=int)
    run(parser.parse_args())
