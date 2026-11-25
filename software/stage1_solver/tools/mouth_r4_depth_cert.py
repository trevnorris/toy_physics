"""Read the step-108 close-approach map and resolve its two orientation enclosures."""
import argparse
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_cp10_label_condition import physical_frame, reference_frame
from mouth_supercore_shape import preparation
from mouth_r4_edge_descent import CELL, NODES, STATE, edge_data
from mouth_r4_eulerian_corner_audit import point
from mouth_r4_eulerian_relax_labels import chart_orientation, cell_order
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_refinement_review import _collar_reference

BITS = tuple(product(range(2),repeat=3))
BOX = tuple(product(range(3),repeat=3))
HESSIAN_SHA = '5c622b400a9d9cd91df008927976d079040cc143459d1dc0084332556d2b0423'


def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def chart_matrix(vertices, widths, local):
    """Exact derivative of the stored Q1 label with respect to radius fraction, depth fraction, angle."""
    D = np.zeros((3,3),dtype=object)
    for j,bits in enumerate(BITS):
        b = [v if k else 1-v for k,v in zip(bits,local)]
        for axis in range(3):
            w = (2*bits[axis]-1)*b[(axis+1)%3]*b[(axis+2)%3]/widths[axis]
            D[:,axis] += vertices[j]*w
    return D


def det(m):
    return sum(m[0,c]*(m[1,(c+1)%3]*m[2,(c+2)%3]-m[1,(c+2)%3]*m[2,(c+1)%3])
               for c in range(3))


def bernstein(nodes):
    vals = nodes.copy()
    for axis in range(3):
        m = np.moveaxis(vals,axis,-1)
        m[...,1] = 2*m[...,1]-(m[...,0]+m[...,2])/2
    return vals


def without_floats(bern, axis):
    """de Casteljau subdivision at exactly one half, preserving rational coefficients."""
    b = np.moveaxis(bern,axis,-1)
    a = (b[...,0]+b[...,1])/2
    z = (b[...,1]+b[...,2])/2
    m = (a+z)/2
    return [np.moveaxis(q,-1,axis) for q in
            (np.stack((b[...,0],a,m),axis=-1),np.stack((m,z,b[...,2]),axis=-1))]


def children(bern):
    for i,a in enumerate(without_floats(bern,0)):
        for j,b in enumerate(without_floats(a,1)):
            for k,c in enumerate(without_floats(b,2)):
                yield (i,j,k), c


def certify(vertices,widths,address):
    root = np.empty((3,3,3),object)
    for ijk in BOX:
        root[ijk] = det(chart_matrix(vertices,widths,[Q(k,2) for k in ijk]))
    root = bernstein(root)
    pending = [((0,0,0),root)]
    log = []
    completed_min = None
    native_leaf = None
    for depth in range(13):
        next_level = []
        uncertain = []
        leaves = []
        for index,b in pending:
            bound = b.min()
            if bound>0:
                leaves.append(bound)
                completed_min = min(completed_min,bound) if completed_min is not None else bound
                continue
            uncertain.append((index,bound))
            if depth==8 and index==address:
                coords = tuple(Q(k,2**depth) for k in index)
                samples = [det(chart_matrix(vertices,widths,[v+Q(j,2**(depth+1)) for v,j in zip(coords,s)])) for s in BOX]
                native_leaf = dict(depth=depth,address=list(index),lo=list(map(float,coords)),
                    hi=[float(v+Q(1,2**depth)) for v in coords],
                    coefficient_index=list(map(int,np.unravel_index(b.argmin(),b.shape))),
                    exact_minimum_coefficient=float(bound),min_of_27_actual_samples=float(min(samples)),
                    minimum_coefficient_after_one_bisection=float(min(c.min() for _,c in children(b))))
            for bits,c in children(b):
                next_level.append((tuple(2*n+j for n,j in zip(index,bits)),c))
        log.append(dict(depth=depth,visited=len(pending),resolved=len(leaves),unresolved=len(uncertain),
                        minimum_open_coefficient=float(min(b for _,b in uncertain)) if uncertain else None))
        if not next_level:
            break
        pending = next_level
    if next_level or native_leaf is None or not completed_min>0:
        raise ValueError('Exact rational Bernstein cover did not close on the native limiting cell')
    return dict(levels=log,certified_lower=float(completed_min),native_leaf=native_leaf)


def metrics(body,axes,labels,order,vertices,widths,coord):
    mat = chart_matrix(vertices,widths,coord)
    def high(v): return mp.mpf(v.numerator)/v.denominator
    ms = mp.matrix([[high(mat[i,j]) for j in range(3)] for i in range(3)])
    chart_s = list(reversed(list(mp.svd(ms,compute_uv=False))))
    # The existing two-metric evaluator takes binary64 coordinates. pm identifies
    # that rounding explicitly for the exact quadratic-edge minimizer.
    pm = tuple(map(float,coord))
    physical = point(body,axes,labels,order,CELL,pm,_collar_reference)
    Ji = physical['Jinv']
    C = dict(min_eigen=1/physical['max_singular_Finv']**2,
             max_eigen=1/physical['min_singular_Finv']**2,det=1/Ji**2,
             inverse_gram_min_eigen=physical['min_singular_Finv']**2)
    return dict(chart_singular_values_ascending=list(map(float,chart_s)),exact_chart_det=float(det(mat)),
                rounded_coordinate_physical=physical,C_equals_FtF=C)


def run(archive):
    tools = Path(__file__).resolve().parent
    signed = tools.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_aligned_hessian.json'
    if digest(signed)!=HESSIAN_SHA:
        raise ValueError('Source edge audit changed')
    pin = archive/'hessian_pin'
    inputs = json.loads(signed.read_bytes())
    manifest = json.loads((pin/'manifest.json').read_bytes())
    if digest(pin/'manifest.json')!=inputs['evaluator_manifest_sha256']['hessian']:
        raise ValueError('Edge evaluator inventory differs')
    sources = {n:digest(tools/n) for n in
               ('mouth_r4_eulerian_solid.cpp','mouth_cp10_spatial_crystal.cpp','mouth_free_shape_solid.cpp',
                'mouth_r4_eulerian_solid.py','mouth_r4_eulerian_relax_labels.py','mouth_cp10_label_grid.py',
                'mouth_r4_eulerian_corner_audit.py','mouth_r4_refinement_review.py','mouth_r4_edge_descent.py')}
    if any(manifest[n]!=v for n,v in sources.items()):
        raise ValueError('The diagnostic must use the frozen geometry and certificate implementations')
    body, cp, axes = get_cp10(archive.parent)
    arm = None
    source = None
    files = {}
    order = None
    for key in ('primary','secondary'):
        here = archive/key
        f = here/'rematched.npz'
        s = inputs['strengths'][key]['hessian']['measurement']
        if digest(f)!=s['state_sha256'] or json.loads((here/'hessian.json').read_bytes())!=s:
            raise ValueError('Saved state or edge probe changed')
        if arm is not None and arm!=s['orientation']:
            if (arm['samples'][2]['physical_rejection']!=s['orientation']['samples'][2]['physical_rejection']
                    or arm['perturbed_radial_node']!=s['orientation']['perturbed_radial_node']):
                raise ValueError('Gauge cases do not reconstruct a common rejected map')
        arm = s['orientation']
        with np.load(f,allow_pickle=False) as data:
            if (str(data['source_sha256'])!=STATE or not np.array_equal(cp['x'],data['x']) or
                    int(data['radius'])!=4 or int(data['winding'])!=1):
                raise ValueError('Unexpected state')
            if source is None:
                source=data['labels']
                order=np.ascontiguousarray(data['field'][2,body.ids].transpose(0,2,1))
            elif not np.array_equal(source,data['labels']):
                raise ValueError('Fixed labels differ between the gauges')
        files[key] = digest(f)
    cp.close()
    w = [Q(float(a[i+1]))-Q(float(a[i])) for a,i in zip(axes,CELL)]
    _,b0,normal,_ = edge_data(source,w)
    p = arm['perturbed_radial_node']
    node = tuple(p['node'])
    step = float(b0/normal[0][0])
    probe = arm['samples'][2]
    if step!=p['direction'] or probe['parameter']!=-.95 or node!=NODES[0]:
        raise ValueError('Original edge displacement changed')
    labels = source.copy()
    labels[node+(0,)] += probe['parameter']*step
    active = cell_order(order)
    cert = chart_orientation(axes,labels,active)
    if cert['chart_det_global_lower']!=probe['whole_chart_lower']:
        raise ValueError('Whole-domain lower bound differs from source')
    try:
        SpatialCrystal(body).evaluate(labels,order,gradient=False)
    except InvalidSpatialCrystal as exc:
        diagnostic=dict(reason=str(exc),**exc.diagnostics)
    else:
        raise ValueError('The original certificate unexpectedly admitted the close approach')
    if probe['physical_rejection']!=diagnostic:
        raise ValueError(f'Native refusal changed: {diagnostic}')
    local = [Q(v) for v in diagnostic['local']]
    address = tuple(int(k*256) for k in local)
    vertices = np.array([[Q(float(a)) for a in labels[n]] for n in NODES],dtype=object)
    proof = certify(vertices,w,address)
    t,j,_,_ = edge_data(labels,w)
    with mp.workdps(80):
        report = {k:metrics(body,axes,labels,order,vertices,w,u) for k,u in
                  [('diagnostic_lo',local),('quadratic_edge_minimum',(Q(0),Q(0),t))]}
    delta = diagnostic['chart_bound'] - proof['native_leaf']['exact_minimum_coefficient']
    margin = 128*np.finfo(float).eps*np.prod([max(np.linalg.norm(
            (labels[tuple(np.array(CELL)+np.array(bits)+np.eye(3,dtype=int)[d])]-
             labels[tuple(np.array(CELL)+np.array(bits))])/float(w[d])) for bits in BITS if bits[d]==0)
            for d in range(3)])
    if not (proof['native_leaf']['exact_minimum_coefficient'] < 0 < proof['native_leaf']['min_of_27_actual_samples']):
        raise ArithmeticError('Exact subbox coefficients and actual samples do not reproduce the native refusal')
    return dict(source_checkpoint='cc11a9f3',inputs=files,source_measurement_sha256=digest(signed),source_digests=sources,
        reconstruction={'node':list(node),'radial_direction':step,'parameter':probe['parameter'],
                        'unchanged_input_value':float(source[node+(0,)]),'trial_value':float(labels[node+(0,)]),
                        'labels_sha256_after':hashlib.sha256(labels.tobytes(order='C')).hexdigest(),
                        'identical_between_gauges':True,'all_eight_cell_orders':order[tuple(np.array(NODES).T)].tolist()},
        predicates={'native':'orientation samples positive, Bernstein minimum <= 0 at depth == 8; returns 3 with lower bound and subbox lower corner',
                    'other':'chart_orientation uses a roundoff margin and de Casteljau subdivision up to depth 12',
                    'native_extra_positive_floor':0.,'chart_margin_at_native_cell':float(margin),
                    'native_location_role':'lower corner of inconclusive subbox; the negative number is a polynomial coefficient enclosure, not a determinant sampled here',
                    'metric_or_cholesky_used_before_refusal':False},
        native=diagnostic,chart=cert,target_cell_exact_proof=proof,native_minus_exact_coefficient=delta,
        exact_cell_bound_minus_roundoff_buffer_minus_reported_global_bound=float(proof['certified_lower']-margin-cert['chart_det_global_lower']),
        original_edge_minimum=float(b0),trial_edge_minimum_fraction=float(t),trial_edge_minimum=float(j),
        points=report, float64={'epsilon':np.finfo(float).eps,'min_positive_normal':np.finfo(float).tiny},
        interpretation='The exact binary-input determinant stays positive on the whole limiting cell. Its depth-eight Bernstein control net straddles zero. The exact cover records the first dyadic subdivision level that certifies the whole cell. The physical strain tensor is never formed on this rejection path.')


def complete_spectra(archive, report):
    """Read orthonormal volume frames on two points of the verified map, without evaluating a solver."""
    body, cp, axes = get_cp10(archive.parent)
    cp.close()
    source = archive/'primary/rematched.npz'
    if digest(source)!=report['inputs']['primary']:
        raise ValueError('Changed source')
    with np.load(source,allow_pickle=False) as f:
        labels=f['labels']
    r=report['reconstruction']
    labels[tuple(r['node'])+(0,)]+=r['radial_direction']*r['parameter']
    if hashlib.sha256(labels.tobytes(order='C')).hexdigest()!=r['labels_sha256_after']:
        raise ArithmeticError('Mismatched trial map')
    vertices = np.array([[Q(float(a)) for a in labels[n]] for n in NODES],dtype=object)
    widths = [Q(float(a[i+1]))-Q(float(a[i])) for a,i in zip(axes,CELL)]
    reference = preparation(4)[0]
    for data in report['points'].values():
        witness = data['rounded_coordinate_physical']
        chart=[Q(v) for v in witness['local']]
        xyz=np.array([witness['physical_chart']])
        xi=np.array([witness['label']])
        current, c_hoop=physical_frame(body,xyz)
        natural, n_hoop=reference_frame(body,reference,xi)
        gradient=chart_matrix(vertices,widths,chart)
        with mp.workdps(80):
            g=mp.matrix([[mp.mpf(v.numerator)/v.denominator for v in row] for row in gradient])
            d=mp.matrix(natural[0])*g*(mp.matrix(current[0])**-1)
            inv_sigma=sorted([*mp.svd(d,compute_uv=False),abs(mp.mpf(float(n_hoop[0]))/mp.mpf(float(c_hoop[0])))])
            C=sorted([1/s**2 for s in inv_sigma])
            relative=abs(float(mp.fprod(inv_sigma))/witness['Jinv']-1)
            if (relative>1e-12 or abs(float(inv_sigma[0])/witness['min_singular_Finv']-1)>1e-12
                    or abs(float(inv_sigma[-1])/witness['max_singular_Finv']-1)>1e-12):
                raise ArithmeticError(f'Orthonormal geometry differs from the high-precision physical frame: {relative}')
            data['physical_singular_values_Finv_ascending']=list(map(float,inv_sigma))
            data['C_equals_FtF']['eigenvalues_ascending']=list(map(float,C))
            data['spectrum_product_relative_error_vs_high_precision_Jinv']=relative
    report['coordinates']='Chart singular values use d(xi)/d(radius fraction, depth fraction, polar angle). Physical singular values use both stress-free and current orthonormal squared-length metrics and the azimuthal hoop. Physical positions r and w are in solver length units and theta is in radians.'
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--extend-metrics',type=Path,help='Existing exact cover; evaluate only the two local frame spectra')
    args=parser.parse_args()
    result=run(args.archive.resolve()) if args.extend_metrics is None else json.loads(args.extend_metrics.read_bytes())
    result=complete_spectra(args.archive.resolve(),result)
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(result,allow_nan=False),flush=True)
