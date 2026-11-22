"""Map inward GNLS crossings at the half-order cut onto a held local patch.

cosine is n.x and clearance_fraction is (R-r)/PITCH on the fixed w cut.
The microscopic x tangent points downstream on the wall, y is the meridian
toward the positive pole, and z points into the fluid. A third wall tangent
orthogonal to the meridian has zero carrier velocity in this axisymmetric
preparation. The patch's orientation, height origin and thermal spread are
additional material assumptions. The two tangent components can be rotated
together into a chosen triangular lattice orientation by the collision engine.
"""
import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import Dipole, trace
from mouth_curved_attachment_supply import gauss
from mouth_curved_site_front import MEASURE, front_read, prepare, sha
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import TEMPERATURE, emit


def record_rows(name, key):
    path = MEASURE / (name+'.json')
    run = json.loads(path.read_text())['guarded_runs'][key]
    if not run['limits_verified'] or run['guard_exit'] or run['child_exit']:
        raise ValueError(f'Unverified input {path}/{key}')
    emit(test='input', record=str(path), run=key, sha256=sha(path))
    return list(map(json.loads, run['literal_stdout'].splitlines()))


@dataclass
class HalfOrderBand:
    depth: float
    radius: float
    wall_slope: float
    scale: float
    pitch: float
    scalar_order: float
    order_dipole: float
    edges: np.ndarray
    coefficients: np.ndarray
    undriven_rate: float

    def save(self, path):
        np.savez(path, **asdict(self))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            return cls(**{k: v.item() if v.ndim == 0 else v for k,v in data.items()})

    def launch(self, drive, cosine, clearance_fraction):
        """Evaluate at common proposal nodes, uniform in cosine and clearance.

        The importance weight is the local axial crossing rate relative to the
        undriven band's mean. At U=0 its joint mean is one. All drives share
        that denominator, retaining the reconstructed (even) throughput change.
        tangent/polar/normal are in microscopic velocity units, with positive
        normal away from the sheet. v_w and j_w retain original mouth units.
        z_offset is normal distance from the smooth tangent plane in native
        length units. It must be added to a separately chosen surface height.
        """
        mu, s = np.broadcast_arrays(cosine, clearance_fraction)
        if (not np.isfinite(drive) or not np.all(np.isfinite(mu)) or not np.all(np.isfinite(s))
                or np.any(abs(mu)>1) or np.any((s<0) | (s>1))):
            raise ValueError('Use finite drive, cosine in [-1,1], clearance in [0,1]')
        index = np.clip(np.searchsorted(self.edges, s, 'right')-1, 0, len(self.edges)-2)
        t = 2*(s-self.edges[index]) / np.diff(self.edges)[index]-1
        c = self.coefficients[index]
        t = t[..., None,None]
        f, fr, fw = np.moveaxis(c[...,0,:,:]+t*(c[...,1,:,:]+t*c[...,2,:,:]), -2, 0)
        field = f[...,0]+drive*mu*f[...,1]
        r = self.radius-PITCH*s
        vw = np.imag((fw[...,0]+drive*mu*fw[...,1])/field)
        vr = np.imag((fr[...,0]+drive*mu*fr[...,1])/field)
        vp = drive*np.sqrt(1-mu*mu)*np.imag(f[...,1]/field)/r
        current = abs(field)**2*vw
        if np.any(vw <= 0):
            raise ValueError('Outward crossing in reconstructed axial stream')
        h = np.hypot(1., self.wall_slope)
        velocity_unit = np.sqrt(self.scale)
        return dict(weight=4*np.pi*PITCH*r*r*current/self.undriven_rate,
                    tangent=(vw+self.wall_slope*vr)/(h*velocity_unit),
                    polar=vp/velocity_unit, normal=(self.wall_slope*vw-vr)/(h*velocity_unit),
                    z_offset=s/h, v_w=vw, j_w=current)


def prepare_half():
    """Resolve the immutable archived carrier, order overlap and calibration."""
    dom, psi, chi = prepare(48)
    w = float(front_read(dom, chi, [], 1.)[0])
    radius, slope = map(float, dom.radius(w))
    dpsi, residual = Dipole(dom,psi).steady(1.,0.)
    fields = np.column_stack([psi,dpsi])
    material = record_rows('mouth_bound_site_transfer', 'primary')
    sample, = [row for row in material if row['test']=='sample']
    binding = record_rows('mouth_curved_companion_order', 'mesh48')
    half = [row for row in binding if row['test']=='binding' and row['location']=='half'
            and row['width']==PITCH]
    if len(half)!=2:
        raise ValueError('Missing half-order contact calibration')
    order, dipole = half[0]['order_average'], half[0]['order_dipole_overlap_per_U']
    np.testing.assert_equal([order,dipole], [half[1]['order_average'],half[1]['order_dipole_overlap_per_U']])

    # At fixed w, f, f_r and f_w are at most quadratic in s *within each
    # element*. The native physical f_w includes the chart slope. Interior
    # nodes retain both one-sided derivative traces at an element interface.
    vertices = radius*dom.blocks[2][1]/PITCH
    edges = np.r_[0., vertices[(vertices>0) & (vertices<1)], 1.]
    t = np.array([-.75,0.,.75])
    points = edges[:-1,None]+(1+t)*np.diff(edges)[:,None]/2
    def radial(x):
        values = trace(dom,2,np.full(x.shape,w),x*PITCH/radius,fields)
        return np.stack(values[-3:],axis=-2)
    exact = radial(points)
    coefficients = np.stack((exact[:,1], (exact[:,2]-exact[:,0])/1.5,
                             (exact[:,2]+exact[:,0]-2*exact[:,1])/(2*.75**2)), axis=1)
    s, weights = gauss(edges, 10)
    f, _, fw = np.moveaxis(radial(s),-2,0)
    current = (f[:,0].conj()*fw[:,0]).imag
    undriven = weights @ (4*np.pi*PITCH*(radius-PITCH*s)**2*current)
    band = HalfOrderBand(w,radius,slope,sample['scale'],sample['pitch'],order,dipole,
                         edges,coefficients,float(undriven))
    reconstructed = band.launch(0.,0.,s)
    err = np.max(abs(reconstructed['j_w']-current))
    emit(test='half_band', depth=w, radius=radius, wall_slope=slope, edges=edges.tolist(),
         scale=band.scale, native_length_in_mouth_units=PITCH,
         native_velocity_in_mouth_units=np.sqrt(band.scale),
         native_time_in_mouth_units=PITCH/np.sqrt(band.scale),
         pitch=band.pitch, scalar_order=order, order_dipole=dipole,
         carrier_residual=residual, undriven_rate=undriven, section_reconstruction_error=err)
    np.testing.assert_allclose(reconstructed['j_w'],current,rtol=0.,atol=2e-13)
    return band


def thermal_at_fixed_axial(band, launch, delta_normal):
    """Add a specified microscopic normal velocity while preserving axial flux.

    This is a conditional radial dispersion at fixed v_w. A genuine isotropic
    thermal reservoir instead requires a new positive-v_w flux weighting.
    """
    return dict(launch, normal=launch['normal']+delta_normal,
                tangent=launch['tangent']-band.wall_slope*delta_normal)


def incident_normal_speed(q, temperature=TEMPERATURE):
    """Positive incoming speed for an independent Maxwellian wall collision.

    q is a uniform quantile. The density of speeds incident through the plane
    is u exp(-u²/2T)/T for a unit-mass thermal constituent. This distribution
    is conditional on hitting the wall, separate from the axial crossing law.
    """
    q = np.asarray(q)
    if np.any(~np.isfinite(q)) or np.any((q<=0)|(q>=1)) or not 0<temperature<np.inf:
        raise ValueError('Use quantiles strictly between 0 and 1 and positive finite temperature')
    return np.sqrt(-2*temperature*np.log1p(-q))


def diagnostics(band):
    s, ds = gauss(band.edges,32)
    mu, dmu = leggauss(22)
    proposal = ds[:,None]*dmu[None,:]/2
    for drive in (0., .002,-.002,.001,-.001):
        row = band.launch(drive, mu, s[:,None])
        measure = proposal*row['weight']
        total = measure.sum()
        p = measure/total
        # Ranges refer to the listed quadrature nodes. The field extrapolated
        # to the surface can switch the sign of a tiny numerical normal speed.
        vals = dict(cosine=mu, **{k:row[k] for k in ('tangent','normal','polar','z_offset','v_w')})
        stats = {}
        for k, a in vals.items():
            v = np.broadcast_to(a, p.shape).ravel()
            ix = np.argsort(v)
            # These are quantiles of discrete quadrature samples. For the
            # physical clearance and grazing angle use the resolved grid below.
            quantiles = np.interp([.1,.5,.9], np.cumsum(p.ravel()[ix]),v[ix])
            stats[k] = dict(mean=p.ravel()@v, sampled_range=[v.min(),v.max()],
                            quantiles_10_50_90=quantiles.tolist())
        # The held wall's continuum normal velocity is exactly zero at its
        # surface. Its Q2 trace has a small numerical boundary residual.
        face = band.launch(drive, np.array([-1.,0.,1.]), 0.)
        emit(test='axial_distribution', drive=drive, total_over_undriven=total,
             observables=stats, increasing_z_normal_fraction=np.sum(p*(row['normal']>0)),
             wall_trace_normal_native_at_mu_minus_zero_plus=face['normal'].tolist())
        if drive==0:
            np.testing.assert_allclose(total,1.,rtol=0.,atol=3e-13)
    emit(test='incident_wall_beam', native_temperature=TEMPERATURE,
         incoming_speeds_at_quantiles_10_50_90=incident_normal_speed([.1,.5,.9]).tolist(),
         mean_incoming_speed=np.sqrt(np.pi*TEMPERATURE/2), rms_incoming_speed=np.sqrt(2*TEMPERATURE),
         native_normal_sigma_for_axial_gaussian=np.sqrt(TEMPERATURE))


def verify_geometry(band):
    """Resolve incidence angles and check the previous gate on this sampler."""
    rates = record_rows('mouth_curved_attachment_counts', 'supply48')
    archived, = [r for r in rates if r['test']=='band' and r['location']=='half'
                 and r['band']=='one_pitch_unit']
    gate_rows = record_rows('mouth_curved_attachment_counts', 'gate48')
    gate, = [r for r in gate_rows if r['test']=='capture_threshold' and r['location']=='half'
             and r['band']=='one_pitch_unit' and r['law']=='retain_mean_attraction']
    s, ds = gauss(band.edges, 12)
    mu, dmu = leggauss(10)
    proposal = ds[:,None]*dmu[None,:]/2
    fraction, cosine, selected = [], [], []
    for u in (0.,.002,-.002,.001,-.001):
        row = band.launch(u,mu,s[:,None])
        p = proposal*row['weight']
        captured = p*-np.expm1(-gate['k_times_L']/row['v_w'])
        fraction.append(captured.sum()/p.sum())
        cosine.append(np.sum(p*mu)/p.sum())
        selected.append(np.sum(captured*mu)/captured.sum())
    measured = 3*(np.array(cosine)[1::2]-cosine[2::2])/np.array([.004,.002])
    captured = 3*(np.array(selected)[1::2]-selected[2::2])/np.array([.004,.002])
    emit(test='previous_gate_crosscheck', rate=band.undriven_rate, archived_rate=archived['rate'],
         measured_incident_beta=measured.tolist(), archived_incident_beta=archived['rate_beta'],
         measured_old_law_fraction=fraction[0], archived_old_law_fraction=gate['zero_flow_fraction_of_band_crossings_captured'],
         measured_old_law_beta=captured.tolist(), archived_old_law_beta=gate['signed_beta'])
    np.testing.assert_allclose(band.undriven_rate, archived['rate'], rtol=1e-11)
    np.testing.assert_allclose(measured, archived['rate_beta'], atol=1e-10,rtol=0)
    np.testing.assert_allclose(fraction[0], gate['zero_flow_fraction_of_band_crossings_captured'], rtol=1e-12)
    np.testing.assert_allclose(captured, gate['signed_beta'], atol=2e-11,rtol=0.)

    # Uniform subcells in physical radial clearance resolve the tiny numerical
    # slip near the wall. Mid-weight quantiles on unique samples center atoms
    # of the discrete quadrature distribution; the zero-flow pole is uniform.
    grid = (np.arange(8192)+.5)/8192
    row = band.launch(0.,0.,grid)
    p = row['weight']/np.sum(row['weight'])
    angle = np.degrees(np.arctan2(-row['normal'], row['tangent']))
    def summary(v):
        ix = np.argsort(v)
        quantiles = np.interp([.1,.5,.9], np.cumsum(p[ix])-.5*p[ix], v[ix])
        return dict(mean=float(p@v), sample_range=[v.min(),v.max()], p10_p50_p90=quantiles.tolist())
    emit(test='grazing_incidence', drive=0., radial_midpoint_cells=len(grid),
         angle_toward_wall_from_tangent_degrees=summary(angle), native_normal_clearance=summary(row['z_offset']),
         native_tangent_speed=summary(row['tangent']), native_into_fluid_normal_speed=summary(row['normal']),
         sampled_outward_from_sheet_fraction=p[row['normal']>0].sum(),
         native_max_normal_displacement_ballistic_in_one_pitch_gate=np.max(abs(row['normal'])*np.sqrt(band.scale)/row['v_w']))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--save', type=Path, help='Build the archived mesh-48 section, save and summarize it')
    group.add_argument('--verify', type=Path, help='Check incidence geometry and the prior hazard using a saved section')
    args = parser.parse_args()
    if args.verify:
        verify_geometry(HalfOrderBand.load(args.verify))
    else:
        band = prepare_half()
        band.save(args.save)
        diagnostics(HalfOrderBand.load(args.save))
