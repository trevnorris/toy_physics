"""Verify R4 step-108 log-depth work on each smooth collar registry branch."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mouth_cp10_field_shape import Bath
from mouth_cp10_label_gauge import depth_grid
from mouth_r4_aligned_crystal import AlignedGauge
from mouth_r4_aligned_derivatives import load, SignedGauge
from mouth_r4_depth_cell_work import SOURCE
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_smooth_endpoint import promoted_bath, SHAPE_WORK_TOL
from mouth_r4_smooth_label_derivatives import CrystalBanks
from mouth_r4_smooth_static import digest, LD


def tangent_kernel(directory):
    header = Path(__file__).with_name('mouth_r4_depth_cell_work.hpp').resolve()
    source = directory/'segment_tangent.cpp'
    source.write_text('#include "audit.cpp"\n#include "'+str(header)+'"\n')
    lib = directory/'segment_tangent.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-shared','-fPIC',str(source),'-o',str(lib)],check=True)
    call = ct.CDLL(str(lib)).depth_cell_gradient
    call.argtypes = [ct.POINTER(HPData),P,P,HPGrid,P,P,P]
    call.restype = ct.c_int
    return call


class FrozenDepth:
    def __init__(self, driver, trial, order, directory):
        self.driver, self.trial, self.order = d, t, order = driver, trial, order
        self.banks = CrystalBanks(t.b, order, directory)
        self.kernel = tangent_kernel(directory)
        dims = tuple(n-1 for n in order.shape)
        self.active = np.maximum.reduce([order[i:i+dims[0],j:j+dims[1],k:k+dims[2]]
                                         for i,j,k in np.ndindex(2,2,2)]).ravel()>0
        self.labels = np.asarray(d.labels, LD)
        self.length = LD(t.shape.length)
        self.origin = depth_grid(self.length,d.reference_depth)
        self.tbase = self.banks.buffers['field_t'].copy()
        self.bath = Bath(t,d.base,d.fields,d.reference_depth)
        self.gauge = AlignedGauge(t,d.reference_depth,order)
        self.signed = SignedGauge(t,d.reference_depth,order)
        self.direction = np.eye(6,dtype=LD)[0]

    def transport(self, q):
        length = self.length*np.exp(q)
        knots = depth_grid(length, self.driver.reference_depth)
        motion = knots-self.origin
        self.banks.buffers['field_t'][:] = self.tbase+motion
        val = np.array(length,LD)
        ct.memmove(ct.addressof(self.banks.data)+HPData.length.offset,val.ctypes.data,val.nbytes)
        moved = self.labels.copy()
        moved[...,1] += motion[None,:,None]
        rates = np.where(knots*length < length-.12,(knots*length/(length-.12))**25,1.)-knots
        return moved, np.ascontiguousarray(rates,LD)

    def energy(self, q):
        moved,_ = self.transport(q)
        dx = q*self.direction
        return (self.banks.energy(moved,self.active),
                promoted_bath(self.bath,self.trial,self.driver.reference_depth,dx),
                LD(self.driver.epsilon)*self.signed.evaluate(self.labels,dx))

    def analytic(self, q):
        moved, rates = self.transport(q)
        b = self.banks
        terms = np.zeros((self.active.size,3),LD)
        rc = self.kernel(ct.byref(b.data),moved.ctypes.data_as(P),b.order.ctypes.data_as(P),
                         b.rule,b.weights.ctypes.data_as(P),rates.ctypes.data_as(P),terms.ctypes.data_as(P))
        if rc or not np.isfinite(terms).all():
            raise ArithmeticError(f'Invalid analytic crystal derivative at {q}')
        # The two smooth contributions use their production complex directional derivative.
        eta = LD('1e-24')
        dx = (q+1j*eta)*self.direction
        return (terms.sum(axis=1,dtype=LD),
                promoted_bath(self.bath,self.trial,self.driver.reference_depth,dx).imag/eta,
                LD(self.driver.epsilon)*self.gauge.evaluate(self.labels,dx).imag/eta)


def compact(parts, cell_indices):
    components = np.array([np.sum(v,dtype=LD) for v in parts],LD)
    return np.r_[components,components.sum(dtype=LD),parts[0][cell_indices]]


def difference(right, left, cell_indices):
    return compact([r-l for r,l in zip(right,left)],cell_indices)


def measure(root, case, out, events):
    out.mkdir(parents=True,exist_ok=False)
    study = root/'r4_rematch39_01'
    if digest(study/case/'rematched.npz') != SOURCE[case]:
        raise ValueError('Unexpected rematch')
    driver, trial, _, order = load(study, case, out)
    field = FrozenDepth(driver,trial,order,out)
    dims = tuple(n-1 for n in field.order.shape)
    cells = [e['cell'] for e in events]
    indices = np.ravel_multi_index(np.asarray(cells).T,dims)
    h = LD(2)**-18
    faces = sorted(set(LD(e['negative_threshold_h']) for e in events))
    if len(faces)!=2 or not (-h < faces[0] < faces[1] < 0):
        raise ValueError('Expected two distinct negative registry-face crossings')
    limits = [-h,*faces,h]
    energy = [field.energy(q) for q in limits]
    parts = []
    work = np.zeros(4+len(cells),LD)
    for left,right,em,ep in zip(limits,limits[1:],energy,energy[1:]):
        width = right-left
        m = (right+left)/2
        # Open two-point Gauss integration never samples the ambiguous face trace.
        gpts = [m+sign*width/(2*np.sqrt(LD(3))) for sign in (-1,1)]
        slopes = [compact(field.analytic(q),indices) for q in gpts]
        integrated = width*(slopes[0]+slopes[1])/2
        observed = difference(ep,em,indices)
        work += integrated
        tangent = compact(field.analytic(m),indices)
        # A distinct centered energy probe fully inside each smooth branch.
        # Keep the shortest probe wide enough to resolve extensive bath energy.
        step = LD('.4')*width
        em_inner, ep_inner = (field.energy(m+sgn*step) for sgn in (-1,1))
        fd = difference(ep_inner,em_inner,indices)/(2*step)
        parts.append(dict(interval=[float(left),float(right)],midpoint=float(m),
                          analytic_at_midpoint=tangent.astype(float).tolist(),
                          centered_width=float(step),centered_slope=fd.astype(float).tolist(),
                          centered_minus_analytic=(fd-tangent).astype(float).tolist(),
                          quadrature_points=list(map(float,gpts)),analytic_work=integrated.astype(float).tolist(),
                          endpoint_energy_change=observed.astype(float).tolist(),
                          work_minus_energy_per_log_depth=((integrated-observed)/width).astype(float).tolist()))
    delta = difference(energy[-1],energy[0],indices)
    wide_error = (work-delta)/(2*h)
    archived = json.loads((root/'r4_rematch39_01/derivatives'/(case+'_shape')/'record.json').read_bytes())['rows'][0]
    # Reconcile to the independent archived total-energy finite difference,
    # which did not use the branch-specific analytic kernel.
    wide_archive_error = -delta[3]/(2*h)-LD(archived['central'][0]['central'])
    worst_inside = max(abs(err) for row in parts for err in row['centered_minus_analytic'])
    worst_segment = max(abs(err) for row in parts for err in row['work_minus_energy_per_log_depth'])
    ok = (max(worst_inside,worst_segment,np.max(abs(wide_error)),abs(wide_archive_error))<=SHAPE_WORK_TOL)
    data = dict(case=case,epsilon=field.driver.epsilon,state_sha256=SOURCE[case],command=sys.orig_argv,
                source_sha256=digest(Path(__file__)),header_sha256=digest(Path(__file__).with_name('mouth_r4_depth_cell_work.hpp')),
                definition='Positive sign is dE/d(log depth). Last four entries are individual physical crystal cells.',
                columns=['crystal','bath','weighted gauge','total']+[str(tuple(c)) for c in cells],
                faces=[dict(cell=e['cell'],ray=e['ray'],log_depth=e['negative_threshold_h']) for e in events],
                active_crystal_cells=int(np.count_nonzero(field.active)),intervals=parts,
                summed_integrated_analytic_work=work.astype(float).tolist(),wide_energy_change=delta.astype(float).tolist(),
                wide_work_minus_energy_per_log_depth=wide_error.astype(float).tolist(),
                archived_wide_force_residual=float(wide_archive_error),max_centered_error=worst_inside,
                max_segment_force_equivalent_error=worst_segment,tolerance=SHAPE_WORK_TOL,passed=bool(ok))
    (out/'result.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:data[k] for k in ('case','state_sha256','max_centered_error','max_segment_force_equivalent_error',
                                         'wide_work_minus_energy_per_log_depth','archived_wide_force_residual','passed')}),flush=True)
    return ok


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('case',choices=SOURCE)
    parser.add_argument('out',type=Path)
    args = parser.parse_args()
    dossier = Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_depth_cell_geometry.json'
    record = json.loads(dossier.read_bytes())
    if record['state_sha256'] != SOURCE['primary']:
        raise ValueError('Unexpected face-trace source')
    sys.exit(0 if measure(args.scratch.resolve(),args.case,args.out.resolve(),record['events']) else 1)
