"""Put two equal-mass, dilute sites on opposite sides of one held converting mouth.

E_site = m mu_ads(int f chi), with the same normalized material-side Gaussian as
the one-site calculation. All responses are first derivatives at m=0. The force
reported here is the medium's order-sector force on the attached wall and sites.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_curved_site_front import prepare, sha, source
from mouth_surface_site_transfer import BASIS, Surface
from neutral_patch_sites import emit


RECORD = Path('research/pde_ledger_v3/_measurements/mouth_bound_site_transfer.json')


def validated(path):
    outcome = json.loads((path/'outcome.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']:
        raise ValueError(f'Incomplete source {path}')
    return list(map(json.loads,(path/'stdout').read_text().splitlines()))


def material_calibration(side):
    record=json.loads(RECORD.read_text())['guarded_runs']['primary']
    assert record['limits_verified'] and record['child_exit']==record['guard_exit']==0
    rows=list(map(json.loads,record['literal_stdout'].splitlines()))
    sample,=[r for r in rows if r['test']=='sample']
    return Surface(side,sample['shell'],sample['pitch'],False),sample['scale']


def run(args):
    dom,psi,chi=prepare(args.mesh)
    front=validated(args.front)
    stress=validated(args.stress)
    centers = [r for r in front if r['test']=='site']
    holding = [r for r in stress if r['test']=='linear_bare_order_wall_load']
    assert [r['label'] for r in centers]==[r['location'] for r in holding]==['entry','half']
    assert all(r['mesh']==args.mesh and r['fields_sha256']==sha(args.fields) for r in holding)
    with np.load(args.fields) as fields:
        np.testing.assert_array_equal(fields['chi'],chi)
        for site in centers:
            np.testing.assert_equal(fields['depth_'+site['label']],site['depth'])
            np.testing.assert_equal(fields['source_dipole_grad_'+site['label']],site['source_dipole_grad'])
    # Using a real unit carrier here sets the integrating factor to one.
    # The overlap in a site binding energy is an ordinary physical-volume average.
    averages=[]
    for site in centers:
        loads,n,_,_=source(dom,np.ones_like(psi),chi,site['depth'],site['width'],12)
        np.testing.assert_allclose(n,site['physical_gaussian_fraction_inside'],rtol=1e-13)
        np.testing.assert_allclose(loads[0].sum(),1.,rtol=1e-13)
        averages.append(loads[0]@chi)
    direct=np.array([r['direct_external_holding']-r['bulk_source'] for r in holding])
    stress_sum=np.array([r['bulk_advection']-sum(r['reservoir_outward_flux'].values()) for r in holding])
    surface,scale=material_calibration(args.side)
    emit(test='input',mesh=args.mesh,side=args.side, width=centers[0]['width'], empty_order_overlaps=averages,
         direct_order_holding_per_h=direct.tolist(),volume_order_holding_per_h=stress_sum.tolist(),
         threshold_half_to_entry_abs_slope=direct[0]/direct[1],
         terminal_sources_sha256=[sha(args.front/'stdout'),sha(args.stress/'stdout')])
    for mean,law in [(False,'weaken_all_contacts'),(True,'retain_mean_attraction')]:
        x=np.zeros(surface.dim)
        x[:2]=surface.pitch*BASIS.sum(axis=0)/3
        x[2]=.75
        slopes=[]
        energies=[]
        for z in averages:
            x,e,slope,_=surface.equilibrate(z,mean,x)
            assert slope<0 and surface.stability(x,z,mean)>0
            slopes.append(scale*slope)
            energies.append(scale*e)
        _,e,_,_=surface.equilibrate(averages[0],mean,x)
        np.testing.assert_allclose(scale*e,energies[0],atol=1e-8,rtol=0.)
        # At first order independent site energies add. The holding load from
        # each site points toward it; the medium-on-assembly force is opposite.
        coeff=direct*np.abs(slopes)
        bulk=stress_sum*np.abs(slopes)
        normals=np.array([1.,-1.])
        net=-normals@coeff
        emit(test='opposite_pair',law=law,specific_binding_derivatives=slopes,
             half_to_entry_abs_slope=abs(slopes[1]/slopes[0]),
             absolute_force_per_single_site_mass=coeff.tolist(),
             medium_force_entry_positive_half_negative_per_site_mass=net,
             bulk_reconstruction_per_site_mass=-normals@bulk,
             exchanged_sites_force_per_site_mass=normals@coeff,
             opposed_entry_sites_force_per_site_mass=-normals@np.repeat(coeff[0],2),
             opposed_half_sites_force_per_site_mass=-normals@np.repeat(coeff[1],2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(26,34,48),required=True)
    parser.add_argument('--side',type=int,choices=(5,8),default=5)
    for name in ('front','stress','fields'):
        parser.add_argument('--'+name,type=Path,required=True)
    run(parser.parse_args())
