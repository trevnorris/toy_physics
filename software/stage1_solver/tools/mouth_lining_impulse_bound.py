"""Phase-independent upper bound on a fresh site's transverse impulse after the flare gauge."""
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_mobile_lining import ROW_SPACING, SUPPORT, row_potential
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, GAP, emit


def continuum_gradient_cap(pitch):
    # Pick the axially nearest row (on the infinite lattice), then the nearest
    # periodic site on each row separately. Image distances obey these bounds for
    # *any* row phases. The actual strip/window only deletes or weakens contacts.
    depth, alpha, rest = BETWEEN
    bound=0.
    images=12
    count=0
    for j in range(-int(np.ceil(SUPPORT/ROW_SPACING))-1, int(np.ceil(SUPPORT/ROW_SPACING))+2):
        dx=max(0.,(abs(j)-.5)*ROW_SPACING)
        if dx>=SUPPORT: continue
        count+=1
        k=np.arange(-images,images+1)
        dy=np.maximum(0.,abs(k)-.5)
        r=np.sqrt(GAP**2+pitch**2*(dx**2+dy**2))
        # r>=gap>rest, |V'|<=2 D alpha exp[-alpha(r-rest)].
        # dr/dy=pitch^2*y/r has magnitude at most pitch.
        bound+=np.exp(-alpha*(r-rest)).sum()
    # For both image tails r>=pitch*(|k|-.5).
    tail=2*count*np.exp(-alpha*(pitch*(images+.5)-rest))/(-np.expm1(-alpha*pitch))
    return float(2*depth*alpha*pitch*(bound+tail))


def continuous_phase_cap(pitch, grid):
    """Upper bound for the exact periodic Morse force on the entire finite strip.

    Each row takes its own worst transverse phase. The sum of row-force
    magnitudes also bounds a finite subset, including the axial taper.
    """
    depth, alpha, rest = BETWEEN
    outer = int(np.floor(SUPPORT/ROW_SPACING+.5))
    offsets = np.linspace(-ROW_SPACING/2,ROW_SPACING/2,grid+1)[:,None]
    phases = np.linspace(-.5,.5,grid+1)[None,:]
    upper_on_grid = np.zeros(len(offsets))
    for j in range(-outer,outer+1):
        grad = row_potential(offsets-j*ROW_SPACING,phases,pitch,True)[1]
        upper_on_grid += np.max(np.abs(grad),axis=1)
    # Omitted periodic images in the grid evaluation, all rows, both tails.
    tail = 2*(2*outer+1)*2*depth*alpha*pitch*np.exp(-alpha*(pitch*10.5-rest))/(-np.expm1(-alpha*pitch))
    # The Hessian of V(r) has radial/tangential eigenvalues V'' and V'/r.
    # Since 0<exp[-alpha(r-rest)]<1, both are bounded by 2D alpha
    # times max(alpha,1/GAP) times this exponential. Both coordinate
    # derivatives of the transverse force are bounded by the Hessian norm.
    hessian_sum = pitch*max(alpha,1/GAP)*continuum_gradient_cap(pitch)
    spacing_error = hessian_sum*(ROW_SPACING+1)/(2*grid)
    return upper_on_grid.max()+tail+spacing_error, upper_on_grid.max(), spacing_error


def run(args):
    trace, intake = cubic_trace(args.mesh)
    patch, _, _ = read_inputs(args.mesh)
    p=patch['pitch']
    pitch=p*PITCH
    w=np.unique(np.r_[np.linspace(1.175,1.65,8193),trace.edges[(trace.edges>1.175)&(trace.edges<1.65)]])
    if args.native_path: w=np.unique(np.r_[w,1.55])
    arc=cumulative_trapezoid(np.hypot(1.,trace.geometry.radius(w)[1]),w,initial=0.)/pitch
    clock=cumulative_trapezoid(1/trace.speed(w),w,initial=0.)
    n=int(np.floor(arc[-1]/ROW_SPACING))+1
    t=np.linspace(0,clock[-1],n*args.stride+1)
    depth=PchipInterpolator(clock,w)(t)
    x=PchipInterpolator(clock,arc)(t)
    gate=int(np.argmin(abs(depth-1.55)))
    t, depth, x=t[gate:],depth[gate:],x[gate:]
    order=trace.sample(trace.chi,depth)[0]
    if order.min()<0: raise ValueError('Negative order')
    ny,nx=args.grid,args.axial_grid
    ys=np.arange(ny+1)[None,:]/ny
    xs=np.arange(SUPPORT*nx+1)[:,None]/nx
    table=row_potential(xs,ys,p,True)[1]
    window=np.clip(SUPPORT-xs,0,1)
    table*=window**2*(3-2*window)
    cap=np.zeros_like(x)
    for j in range(n):
        dx=abs(x-j*ROW_SPACING)*nx
        mask=dx<SUPPORT*nx
        i=dx[mask].astype(int)
        frac=(dx[mask]-i)[:,None]
        # At this axial separation the force is piecewise linear in transverse
        # phase, so its absolute maximum lies on a phase-grid node.
        row=(1-frac)*table[i]+frac*table[i+1]
        cap[mask]+=np.max(abs(row),axis=1)
    scale=2*patch['shared_scale']/pitch**2
    dt=np.diff(t)
    acceleration=scale*order*cap
    bound=np.sum(.5*dt*(acceleration[1:]+acceleration[:-1]))
    envelope=continuum_gradient_cap(p)
    # For the exact quadratic wall-order interpolant use its maximum throughout
    # each actual intersected element, not just the endpoints or sampled max.
    maxima=[]
    for j,(lo,hi) in enumerate(zip(trace.edges[:-1],trace.edges[1:],strict=True)):
        left=max(lo,depth[0])
        right=min(hi,depth[-1])
        if right<left: continue
        nodes=trace.chi[j]
        A=2*(nodes[0]-2*nodes[1]+nodes[2]); B=-3*nodes[0]+4*nodes[1]-nodes[2]
        us=[(left-lo)/(hi-lo),(right-lo)/(hi-lo)]
        if A and us[0]<-B/(2*A)<us[1]: us.append(-B/(2*A))
        maxima.append(max(abs((A*u+B)*u+nodes[0]) for u in us))
    order_max=max(maxima)
    if order_max < order.max()-1e-12: raise ValueError('Order supremum failed')
    # Duration is from the archived carrier integration. Bound conditioning on
    # that duration avoids any numerical assertion about a true continuum PDE.
    smooth_cap, grid_value, grid_margin = continuous_phase_cap(p,args.envelope_grid)
    emit(test='absolute_downstream_impulse_cap',mesh=args.mesh,intake=intake,phase_grid=ny,axial_grid=nx,native_path=args.native_path,
         steps_per_packet=args.stride,gauge_depth=depth[0],duration=t[-1]-t[0],
         discrete_impulse_cap_pitches_per_time=bound,unwindowed_row_gradient_envelope=envelope,
         abs_order_supremum=order_max, exact_force_gradient_cap=smooth_cap,
         sampled_gradient_peak=grid_value, derivative_margin=grid_margin, envelope_grid=args.envelope_grid,
         sharper_continuous_impulse_cap=scale*order_max*smooth_cap*(t[-1]-t[0]),continuous_force_cap_pitches_per_time_squared=scale*order_max*envelope,
         continuum_cap_at_sampled_duration=scale*order_max*envelope*(t[-1]-t[0]))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('mesh','stride','grid','axial-grid','envelope-grid'): parser.add_argument('--'+name,required=True,type=int)
    parser.add_argument('--native-path',action='store_true')
    run(parser.parse_args())
