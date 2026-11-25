"""Finite R4 shape trials on the step-395 cubic, two-metric internal functional."""
import argparse
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
import gc
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth
from mouth_cp10_label_gauge import depth_grid
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_free_shape_solid import Solid
from mouth_joint_endpoint_audit import geometry
from mouth_r4_aligned_crystal import AlignedCrystal
from mouth_r4_aligned_derivatives import shape_bank
from mouth_r4_aligned_rematch import FIELDS, front
from mouth_r4_barrier_shape import BarrierShapeWork
from mouth_r4_cubic_held import R4_COEFFICIENT_SHA, approved_reference
from mouth_r4_cubic_physical_orientation import physical_certificate
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, chart_orientation
from mouth_r4_eulerian_release import Driver, block_axes, order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_harmonic_shape import ReducedShapeWork
from mouth_r4_slaved_melt import FLOW
from mouth_r4_smooth_fields import GaugeGeometry, read_fields
from mouth_r4_smooth_static import Gauge
from mouth_whole_r4_relax import atomic_json

ROOT = Path(__file__).resolve().parents[3]
SCRATCH = ROOT/'_scratch/s11c'
ARCHIVE = ROOT/'research/pde_ledger_v3/_measurements/mouth_r4_held_pair395.json'
ARCHIVE_SHA = 'e3434eeae512e882573ebe2212b294ab22a24a0ea0ec2d1bcd4b73c5f542c7e9'
PIN = SCRATCH/'r4_adaptive382_01/pinned/manifest.json'
PIN_SHA = 'bfe9c1e8897cd679f082ff930b4e9299d2adbcbc27a3e815db7fc1db294aebe0'


@contextmanager
def current_functional_only():
    def forbidden(*args, **kwargs):
        raise AssertionError('Legacy solid, shape objective or vertex-only gauge reached the R4 release')
    with ExitStack() as stack:
        for cls, method in ((Driver, 'parts_at'), (Solid, 'evaluate'),
                            (BarrierShapeWork, '__call__'), (ReducedShapeWork, '__call__'),
                            (Gauge, 'cell_barrier')):
            stack.enter_context(patch.object(cls, method, forbidden))
        yield


def verify_sources():
    if digest(ARCHIVE) != ARCHIVE_SHA or digest(PIN) != PIN_SHA:
        raise ValueError('Step-395 measurement or solver pin changed')
    for path, sha in read(PIN).items():
        if digest(ROOT/path) != sha:
            raise ValueError(f'Prior solver source changed: {path}')


@dataclass
class Snapshot:
    q: np.ndarray
    labels: np.ndarray
    fields: np.ndarray
    t_axis: np.ndarray


def start_state(model):
    rec = read(ARCHIVE)['endpoints']['primary']
    file = ROOT/rec['state_file']
    if digest(file) != rec['state_sha256'] or not rec['passed']:
        raise ValueError('Unverified step-395 primary state')
    with np.load(file, allow_pickle=False) as data:
        if (int(data['outer']) != 395 or int(data['radius']) != 4 or int(data['winding']) != 1
                or float(data['epsilon']) != 1e-8 or str(data['parent_sha256']) != CHECKPOINT_SHA
                or not np.array_equal(data['x'], rec['shape'])):
            raise ValueError('Wrong primary preparation')
        # Use the original binary64 production knots, including roundoff.
        length = model.shape(data['x']).length
        axis = axes_for_depth(length, model.reference_depth)[1][1]/length
        return Snapshot(data['x'], data['labels'], data['field'], axis), rec


class GeometryFailure(ValueError):
    """A proposed geometry or its transported material chart is infeasible."""


class ReleaseModel:
    # Only geometry construction and the order-knot diagnostic are imported
    # from the former outer driver. Its energy/step methods are inaccessible.
    shape = Driver.shape
    trial = Driver.trial
    forces_at_order_knots = Driver.forces_at_order_knots

    def __init__(self, folder, epsilon):
        Driver.__init__(self, SCRATCH, folder, epsilon)
        self.x = self.fields = self.labels = None

    def open(self, saved, q):
        self.x = np.asarray(q, dtype=float).copy()
        if self.x.shape != (6,) or not np.isfinite(self.x).all():
            raise GeometryFailure('Nonfinite six-coordinate shape')
        shape = self.shape(self.x)
        if not .35 < shape.length < 21.5:
            raise GeometryFailure(f'Axial guard at {shape.length}')
        _, zs = axes_for_depth(shape.length, self.reference_depth)
        if any(not np.all(np.diff(z)>0) for z in zs):
            raise GeometryFailure('Depth mesh is not strictly increasing')
        trial = self.trial()
        approved_reference(trial.b)
        trial.smooth_order_geometry = GaugeGeometry(trial, aligned=True)
        solid = SpatialCrystal(trial.b)
        axis = block_axes(trial.b)[1]
        self.labels = saved.labels.copy()
        self.labels[...,1] += (axis-saved.t_axis)[None,:,None]
        self.fields = trial.transported(saved.fields)
        self.fields[0,~trial.model.masks[0]] = 1.
        self.fields[1,trial.model.feed] = 0.
        self.fields[1,trial.model.bulk] = .4510084282857874
        self.current_trial = trial
        certificate = self.orientation(trial)
        return trial, solid, certificate

    def orientation(self, trial):
        axes = block_axes(trial.b)
        order = order_grid(trial, self.fields)
        active = np.lib.stride_tricks.sliding_window_view(order,(2,2,2)).max(axis=(-3,-2,-1))>0
        try:
            chart = chart_orientation(axes,self.labels,active)
            cert = physical_certificate(trial.b,self.labels,chart['chart_det_global_lower'])
        except (ArithmeticError,ValueError) as exc:
            raise GeometryFailure(f'Whole-cell certificate: {exc}') from exc
        lo, hi = cert['physical_profile_bounds']
        if (cert['cell_count']!=28080 or not cert['physical_global_lower']>0
                or not lo>3*abs(trial.polar_shift)+.05 or not hi*np.exp(.125)<27.5):
            raise GeometryFailure(f'R4 guard band or orientation: {cert}')
        return cert

    def snapshot(self, trial):
        return Snapshot(self.x.copy(), self.labels.copy(), self.fields.copy(), block_axes(trial.b)[1].copy())

    def observe(self, trial, solid, cert):
        order = order_grid(trial,self.fields)
        joint = AlignedCrystal(self,trial,solid,self.labels,order)
        _, _, material = joint.evaluate(self.labels)
        physical,load,field = read_fields(self,trial,solid,self.labels,self.fields)
        if abs(field['gauge_energy']-material['gauge_energy'])>1e-15:
            raise ArithmeticError('Order and crystal gauge work disagree')
        flow = trial.measure(self.fields,load)
        residuals = dict(label=material['max_free_label_force'],
                         melted_label=self.epsilon*material['max_unit_melt_residual'],
                         unit_melt=material['max_unit_melt_residual'], **{k:field[k] for k in FIELDS})
        geo = geometry(trial.shape,trial.polar_shift)['extrema']
        low = min(geo.values(),key=lambda v:v['minimum'])
        w = front(trial,self.fields)
        return dict(shape=self.x.tolist(),depth=trial.shape.length,exit_radius=trial.shape.exit,
                    min_radius=low['minimum'],min_radius_t=low['narrowest_depth']/trial.shape.length,
                    max_radius=max(v['maximum'] for v in geo.values()), total_energy=physical+material['gauge_energy'],
                    crystal_energy=material['elastic_energy'],gauge_energy=material['gauge_energy'],
                    native_residuals=residuals,max_internal=max(v for k,v in residuals.items() if k!='unit_melt'),
                    certified_Jinv=cert['physical_global_lower'],certificate=cert,
                    current={k:float(flow[k]) for k in FLOW},sampled_melted_front=w,
                    front_fraction=w/trial.shape.length if w is not None else None,
                    free_zero=int(np.sum(trial.free & (self.fields[2,trial.first]==0))))


def forces(model, trial, solid):
    bank = shape_bank(model,trial,solid,order_grid(trial,model.fields))
    rows = bank['rows']
    if not bank['passed'] or len(rows)!=6:
        raise ArithmeticError(f'Moved shape signed derivative failure: {rows}')
    curvature = []
    for row in rows:
        wide, narrow = row['central']
        if (wide['h'],narrow['h']) != (2.**-18,2.**-20):
            raise ValueError('Unexpected signed bank widths')
        curvature.append([(a['negative']-a['positive'])/a['h'] for a in (wide,narrow)])
    return np.array([r['analytic'] for r in rows]), np.asarray(curvature), bank


def audit(folder):
    verify_sources()
    model = ReleaseModel(folder,1e-8)
    state, prior = start_state(model)
    trial, solid, cert = model.open(state,state.q)
    base = model.observe(trial,solid,cert)
    f, curvature, bank = forces(model,trial,solid)
    deltas = {k:base[k]-prior[k] for k in ('total_energy','crystal_energy','gauge_energy')}
    for k in ('native_residuals','current'):
        if set(base[k])!=set(prior[k]):
            raise ValueError(f'Missing baseline channels: {k}')
        deltas.update({f'{k}.{s}':x-prior[k][s] for s,x in base[k].items()})
    deltas['certified_Jinv'] = base['certified_Jinv']-prior['orientation']['physical_certificate']['physical_global_lower']
    df = f-np.array(prior['shape_forces'])
    result = dict(source_sha256=prior['state_sha256'],code_sha256=digest(Path(__file__)),
                  accepted_shape_steps=0,baseline=base,shape_force=f.tolist(),
                  held_deltas=deltas,force_deltas=df.tolist(),reference_sha256=R4_COEFFICIENT_SHA,
                  curvatures_wide_narrow=curvature.tolist(),curvature_convention='transported, negative minus positive one-sided force divided by h',
                  shape_audit=bank,probes=[])
    atomic_json(folder/'audit.json',result)
    if any(v!=0 for v in deltas.values()) or np.any(df!=0):
        raise ArithmeticError('Held-state invariance failed; see audit.json')
    delta = 1.e-4*f/max(abs(f))
    model.current_trial = None
    del trial, solid
    gc.collect()
    for scale in (1., .5):
        trial, solid, cert = model.open(state,state.q+scale*delta)
        moved = model.observe(trial,solid,cert)
        work = moved['total_energy']-base['total_energy']
        linear = -float(f@(model.x-state.q))
        result['probes'].append(dict(scale=scale,observables=moved,delta_energy=work,
                                     linear=linear,remainder=work-linear))
        atomic_json(folder/'audit.json',result)
        model.current_trial = None
        del trial, solid
        gc.collect()
    full,half = result['probes']
    noise = 16*abs(np.spacing(base['total_energy']))
    if not (full['delta_energy']<0 and abs(full['remainder'])<=.02*abs(full['linear'])+noise
            and abs(half['remainder'])<=.35*abs(full['remainder'])+noise):
        raise ArithmeticError('Transported directional work failed its first-order/half-width check')
    result['passed'] = True
    atomic_json(folder/'audit.json',result)
    print('Held-state and transported-work checks passed.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    args=parser.parse_args()
    output=args.folder.resolve()
    output.mkdir(parents=True,exist_ok=False)
    try:
        with current_functional_only():
            audit(output)
    except Exception as exc:
        atomic_json(output/'ISSUE.json',dict(type=type(exc).__name__,reason=str(exc)))
        raise
