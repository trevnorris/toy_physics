"""Holding force of an independently compression-calibrated trial mouth coating.

Cold, fixed-gap two-layer cells cover a three-dimensional wall with a specified
reference thickness. Registered sites follow the deformed wall locally. Shear
order is frozen at the stationary intake during this static constitutive test.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_inertia import taper
from mouth_empty_bend_stiffness import INTERVALS
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_periodic_compression import Bilayer
from neutral_patch_sites import GAP, emit, morse

DATA = Path('research/pde_ledger_v3/_measurements')


def read_inputs(mesh):
    calibration = json.loads((DATA/'neutral_patch_zero_load.json').read_text())['production']
    record = json.loads((DATA/'mouth_occupied_companion.json').read_text())
    mouth = record['guarded_runs'][f'occupied_companion_n{mesh}_01']
    if calibration['exit_code'] or mouth['guard_exit'] or mouth['child_exit'] or not mouth['limits_verified']:
        raise ValueError('Unverified calibration or occupied carrier record')
    patch = json.loads(calibration['literal_stdout'].splitlines()[0])
    rows = {r['test']: r for r in map(json.loads, mouth['literal_stdout'].splitlines())}
    if patch['fit_density'] != 1.01 or rows['wave']['action'] != 1.:
        raise ValueError('This test uses the neutral 1% calibration and one unit of wave action')
    preparation = rows['preparation']
    if abs(preparation['rate']-.25)>1.e-9 or (preparation['far'],preparation['roof'])!=(9.,10.):
        raise ValueError('Unmatched bath or intake')
    baseline = np.array(rows['static']['carrier']) + np.array(rows['wave']['static'])
    return patch, baseline, np.diag(rows['static']['area'])


def surface_rule(trace, span):
    low, high = span
    edges = np.r_[low, trace.edges[(trace.edges>low)&(trace.edges<high)], high]
    points, weights = leggauss(6)
    half = np.diff(edges)[:,None]/2
    depths = (edges[:-1,None] + half*(points+1)).ravel()
    radius, slope = trace.geometry.radius(depths)
    arc = np.hypot(1., slope)
    h_s = -taper(depths,low,high)[1]/arc
    area = (half*weights).ravel()*4*np.pi*radius**2*arc
    order = trace.sample(trace.chi,depths)[0]
    # The coverage law is the physical [0,1] part of the interpolated field.
    if order.min() < -1.e-5 or order.max() > 1.00001:
        raise ValueError('Order outside its physical range')
    return area, h_s, slope/arc, np.clip(order,0.,1.)


def bond_forces(cell, pitch, q, h_s, meridian_x, angular_x, angular_weights):
    """Local derivative of the finite 4D squared bond lengths for a meridional sheet.

    A bond's in-sheet vector is (a,b). Its shifted vector gains q*a*h_s e_x.
    The normal layer gap stays perpendicular to that shifted plane.
    """
    a, b = (cell.vectors*pitch).T
    projection = meridian_x[:,None]*a + angular_x[:,None]*b
    velocity = h_s*a
    squared = a*a + b*b + cell.cross*GAP**2
    r = np.sqrt(squared+2*q*velocity*projection+(q*velocity)**2)
    _, derivative = morse(r,*cell.params)
    return angular_weights @ np.sum(derivative*velocity*(projection+q*velocity)/r,axis=1)


def calculate(mesh):
    patch, baseline, area_per_tension = read_inputs(mesh)
    trace, intake = cubic_trace(mesh)
    if abs(intake-.25)>1.e-9:
        raise ValueError('Order trace has a different intake')
    # Multiplying energy per reference four-volume by normal depth converts
    # the calibration into energy per reference three-area of represented wall.
    thickness = 2*GAP*PITCH
    surface_scale = thickness*patch['shared_scale']
    modes = patch['homogeneous_modes']
    mu = surface_scale*modes['area_preserving_shear']/4
    c11 = surface_scale*(modes['dilation']+modes['area_preserving_shear'])/4
    emit(test='calibration', mesh=mesh, intake=intake, normal_thickness=thickness,
         pitch=patch['pitch'], surface_scale=surface_scale, planar_shear=mu,
         longitudinal=c11, rest_stress=patch['in_plane_stress'])
    cell = Bilayer(1,1,True,patch['shell'])
    u, wu = leggauss(6)
    phi = np.arange(8)*2*np.pi/8
    normal_x = np.repeat(u,len(phi))
    angular_x = (np.sqrt(1-u[:,None]**2)*np.cos(phi)).ravel()
    angular_weights = np.repeat(wu/(2*len(phi)),len(phi))
    signed_steps = np.array([1.e-4,-1.e-4,5.e-5,-5.e-5,0.])
    for i, span in enumerate(INTERVALS,start=1):
        area, h_s, meridian_radial, order = surface_rule(trace,span)
        # For a uniformly oriented spherical section and transverse tangent,
        # <n_x²> = <t_x²> = 1/3; <n_x t_x> = 0.
        kernel = area*h_s**2*(c11*meridian_radial**2+mu)/3
        per_depth = np.array([
            [bond_forces(cell,patch['pitch'],q,hs,radial*normal_x,angular_x,angular_weights)
             for q in signed_steps] for hs,radial in zip(h_s,meridian_radial,strict=True)])
        for label, coverage in [('persistent',np.ones_like(order)),('ordered_fraction',order)]:
            result = kernel@coverage
            forces = surface_scale*(area*coverage)@per_depth
            signed_slopes = (forces[:-1]-forces[-1])/signed_steps[:-1]
            np.testing.assert_allclose(signed_slopes, result, rtol=3.e-5,atol=2.e-9)
            emit(test='clamp', bend=i, engagement=label, material_slope=result,
                 material_force_zero=forces[-1], steps=signed_steps[:-1].tolist(),
                 signed_force_per_step=signed_slopes.tolist(), carrier_plus_wave=baseline[i,i],
                 total=baseline[i,i]+result, equivalent_area_tension=result/area_per_tension[i],
                 candidate_totals=(baseline[i,i]+np.array([.02,.04])*area_per_tension[i]).tolist())
        # Each bond of a sheet lying within one spherical section translates
        # by the same amount at both ends. Its exact length is unchanged.
        emit(test='circumferential_control',bend=i,material_slope=0.,total=baseline[i,i])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(26,34),required=True)
    calculate(parser.parse_args().mesh)
