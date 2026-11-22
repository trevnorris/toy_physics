"""Local signed force, front velocity and edge support in the matched-curvature controls.

A neutral signed-distance order kink is prepared on a spherical cap in 4D.
The phased carrier initially has n=1 and a uniform velocity along w. The
inverse local order equation gives the instantaneous force required to hold
this prepared profile. It does not prescribe an unmeasured rim constitutive
law or claim the two-sided baths maintain that profile at all later times.
"""
import json

import numpy as np
from scipy.integrate import simpson
from scipy.special import expit

from mouth_flat_phase_control import saved_bath


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def local(radius, sample_r, span, count, speed):
    kappa, mobility = 2., .5
    s = np.linspace(-span, span, count)
    q = expit(-s)
    dq = -q*(1-q)
    ddq = dq*(2*q-1)
    radial = radius+s
    weight = (radial/radius)**3
    normal_w = np.sqrt(1.-(sample_r/radius)**2)
    u_n = normal_w*speed
    well = (q*(1-q))**2
    bare_mu = 2*q*(1-q)*(1-2*q)-kappa*(ddq+3*dq/radial)
    holding_mu = -u_n*dq/mobility-bare_mu
    hold = simpson(weight*holding_mu*dq,x=s)
    preload = -simpson(weight*bare_mu*dq,x=s)
    def energy(R):
        # Fixed solid angle. Differentiate, then normalize by reference area.
        return simpson((R+s)**3*(well+.5*kappa*dq**2),x=s)/radius**3
    step = .001
    energy_derivative = (energy(radius+step)-energy(radius-step))/(2*step)
    center=count//2
    dot_order = -u_n*dq-mobility*bare_mu
    released_speed=-dot_order[center]/dq[center]
    zero_held = np.max(abs(u_n*dq+mobility*(bare_mu+holding_mu)))
    shape_norm = simpson(weight*dq**2,x=s)
    emit(test='local_order_hold',radius=radius,sample_r=sample_r,normal_w=normal_w,
         physical_halfwidth=span,native_integration_samples=count,carrier_speed_w=speed,
         current_normal=u_n,normal_initial_level_speed=released_speed,
         initial_relative_material_current=u_n-released_speed,
         normal_external_order_force=hold,neutral_curvature_preload=preload,
         area_energy_derivative=energy_derivative,flow_part=hold-preload,
         direct_shape_norm=shape_norm,linear_force_per_w_speed=-shape_norm*normal_w/mobility,
         held_order_rate_residual=zero_held)
    if (abs(preload-energy_derivative)>2e-10 or zero_held>1e-13
            or abs(released_speed-(u_n-3*mobility*kappa/radius))>1e-12):
        raise ArithmeticError('Local curvature, force or order-rate check failed')


def edge(radius, inner, outer):
    # Virtual common w translation of an ideal thin, constant-tension patch.
    sigma = 1/3
    proj_area = 4*np.pi*(outer**3-inner**3)/3
    flange = 4*np.pi*sigma/radius
    distributed = 3*sigma*proj_area/radius
    inner_w = flange*inner**3
    outer_w = -flange*outer**3
    emit(test='ideal_edge_translation',radius=radius,inner=inner,outer=outer,
         tension=sigma,distributed_external_w=distributed,inner_conormal_w=inner_w,
         outer_conormal_w=outer_w,total_virtual_work_per_unit_translation=distributed+inner_w+outer_w)
    if abs(distributed+inner_w+outer_w) > 1e-12:
        raise ArithmeticError('Unforced common translation does capillary work')


if __name__=='__main__':
    radius, full = 100., saved_bath()/32.
    emit(test='preparation',order='logistic signed distance, ordered inside',
         span=32.,radius=radius,opening_radius=8.,outer_radius=16.,equal_density=1.,
         mobility=.5,kappa=2.,trial_action=0.,force_units='per three-area at the sampled ray',
         free_time='instant immediately after order holder removed without impulse')
    for points, width in ((1601,12.),(3201,14.)):
        for velocity in (0., full, -full, .5*full, -.5*full):
            local(radius,12.,width,points,velocity)
    for inner in (0.,8.):
        edge(radius,inner,16.)
