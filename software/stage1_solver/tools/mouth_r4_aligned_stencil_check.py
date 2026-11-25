"""Check the 39-point determinant rule on a Q1 edge inversion and saved R4 labels."""
from hashlib import sha256
import json
import sys
from pathlib import Path

import numpy as np

from mouth_r4_aligned_stencil import (EDGES, EDGE_WEIGHT, GRID, MIDCOUNT, VERTICES,
                                      WEIGHTS, barrier_means, gradient_basis, samples)


def direct(points, widths, q):
    jac = (points-points[:1]).T@gradient_basis(q)/widths
    return jac[:,0]@np.cross(jac[:,1],jac[:,2])


def response(points, widths, direction, span):
    """Energy slopes and two centered derivatives of the complex-step slope."""
    tangent = np.clongdouble('1.e-25j')
    def slope(s):
        shifted = points.astype(np.clongdouble)+(s+tangent)*direction
        return np.array(barrier_means(shifted[None],widths)).ravel().imag/tangent.imag
    return slope(0), [((slope(h)-slope(-h))/(2*h)).tolist() for h in (span,span/2)]


def audit():
    scale = np.longdouble
    widths = np.ones((1,3),scale)
    shift = scale('0.375')
    x,y,z = VERTICES.T.astype(scale)
    controls = []
    for margin in map(scale,('0.00001','0','-0.00001')):
        # All vertical edge polynomials are (z-3/8)^2+margin.
        labels = np.array([(z-shift)*x-y, margin*x+(z-shift)*y, z]).T
        det,t,adapt = samples(labels[None],widths)
        explicit = np.array([direct(labels,widths[0],GRID[l]/2+(GRID[r]-GRID[l])*v/2)
                             for (l,m,r),v in zip(EDGES,t[0])])
        controls.append(dict(analytic_min=float(margin),fixed_min=float(det.min()),
                             moving_min=float(adapt.min()),max_direct_error=float(max(abs(explicit-adapt[0])))))
    synthetic_curvature = []
    direction = np.array([0*x,x,0*x]).T
    for margin in map(scale,('0.0001','0.000001','0.00000001')):
        labels = np.array([(z-shift)*x-y,margin*x+(z-shift)*y,z]).T
        measured,pairs = response(labels,widths,direction,margin/1024)
        fixed_z = GRID[:,2].astype(scale)/2
        fixed_slope = 1-1/((fixed_z-shift)**2+margin)
        old_slope = fixed_slope[MIDCOUNT==0].mean()
        aligned_slope = WEIGHTS@fixed_slope
        for left,mid,right in EDGES:
            if GRID[left,2]==GRID[right,2]:
                aligned_slope += EDGE_WEIGHT*fixed_slope[mid]
            else:
                a, c = shift**2+margin, (1-shift)**2+margin
                t = np.sqrt(a)/(np.sqrt(a)+np.sqrt(c))
                dt = t*(1-t)*(1/a-1/c)/2
                j = (t-shift)**2+margin
                aligned_slope += EDGE_WEIGHT*(1-1/j)*(1+2*(t-shift)*dt)
        aligned_slope /= 216
        synthetic_curvature.append(dict(margin=float(margin),
            old_and_aligned_slopes=measured.astype(float).tolist(),
            analytic_slopes=[float(old_slope),float(aligned_slope)],
            analytic_minus_complex_relative=float(max(abs(np.array([old_slope,aligned_slope])/measured-1))),
            signed_old_and_aligned=[[float(x) for x in row] for row in pairs],
            margin_squared_times_aligned=float(margin**2*pairs[-1][1])))
    # A constant affine chart has the original value under the normalized rule.
    chart = np.array([2*x+y,-x+3*y,z]).T[None]
    old,new = barrier_means(chart,widths)
    weight = int(WEIGHTS.sum()+len(EDGES)*EDGE_WEIGHT)
    # Use the separately verified cell volume. Scalar J depends only on the
    # product of chart widths, so one normalized column suffices for this test.
    study = Path('_scratch/s11c/r4_newton_held01/primary/last_active.npz')
    edge_path = Path('research/pde_ledger_v3/_measurements/mouth_r4_newton_edge.json')
    oracle = json.loads(edge_path.read_bytes())
    state_sha = sha256(study.read_bytes()).hexdigest()
    if state_sha != 'a2210361aa5b38712de7db618b5220985e12f39f3a089351dfd3767bef0a7178':
        raise ValueError('Different step-108 state')
    with np.load(study,allow_pickle=False) as saved:
        points = saved['labels'][tuple((VERTICES+[25,11,10]).T)].astype(scale)
    widths[0,0] = oracle['normalizing_cell_volume']
    det,t,adapt = samples(points[None],widths)
    e = EDGES.index((0,1,2))
    check = direct(points,widths[0],np.array([0,0,t[0,e]],dtype=scale))
    row = dict(step=108,cell=[25,11,10],old_corner_min=float(det[0,MIDCOUNT==0].min()),
               fixed_lobatto_min=float(det.min()),blocking_edge_midpoint=float(det[0,1]),
               moving_fraction=float(t[0,e]),moving_determinant=float(adapt[0,e]),
               direct_determinant=float(check),stationary_minimum=oracle['min_on_blocking_edge'],
               frozen_corner_mean=float(barrier_means(points[None],widths)[0][0]),
               frozen_aligned_mean=float(barrier_means(points[None],widths)[1][0]))
    # Normalize a single nodal coordinate to unit change in the edge minimum.
    # At an interior stationary point its displacement makes no linear
    # contribution. No independent coordinate is changed in the saved state.
    q = np.array([0,0,oracle['minimizing_fraction_on_blocking_edge']],dtype=scale)
    unit = np.eye(points.size,dtype=scale).reshape(points.size,8,3)
    height = np.clongdouble('1.e-25j')
    slopes = np.array([direct(points+height*v,widths[0],q).imag/height.imag for v in unit])
    j = np.argmax(abs(slopes))
    span = scale(oracle['min_on_blocking_edge'])/1024
    frozen_slope,frozen_curvature = response(points,widths,unit[j]/slopes[j],span)
    row['local_stiffness'] = dict(vertex=VERTICES[j//3].tolist(),component=int(j%3),
        unit_change_to_stationary_edge_det=float(slopes[j]),
        signed_old_and_aligned_opening_slopes=frozen_slope.astype(float).tolist(),
        widths=[float(span),float(span/2)],
        signed_old_and_aligned=[[float(x) for x in r] for r in frozen_curvature])
    result=dict(command=sys.orig_argv,script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                stencil_sha256=sha256(Path(__file__).with_name('mouth_r4_aligned_stencil.py').read_bytes()).hexdigest(),
                edge_oracle_sha256=sha256(edge_path.read_bytes()).hexdigest(),
                synthetic_curvature=synthetic_curvature,state_sha256=state_sha,determinant_weight_sum=weight,
                affine_mean_difference=float((new-old)[0]),synthetic=controls,saved=row)
    print(json.dumps(result,allow_nan=False,indent=2),flush=True)
    assert all(r['analytic_minus_complex_relative']<1e-10 for r in synthetic_curvature)
    assert abs(synthetic_curvature[-1]['margin_squared_times_aligned']-1/27)<1e-7
    assert frozen_curvature[-1][1]>frozen_curvature[-1][0]>0
    assert max(abs(np.asarray(frozen_curvature[0])/frozen_curvature[1]-1))<5e-5
    assert weight==216 and abs((new-old)[0])<1e-17
    assert min(c['fixed_min'] for c in controls)>0
    assert controls[0]['moving_min']>0 and abs(controls[1]['moving_min'])<1e-16 and controls[2]['moving_min']<0
    assert all(c['max_direct_error']<1e-16 for c in controls)
    assert abs(check-adapt[0,e])<1e-17
    assert 0<oracle['min_on_blocking_edge']<=row['moving_determinant']<row['fixed_lobatto_min']


if __name__=='__main__':
    audit()
