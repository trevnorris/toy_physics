"""Read-only cell dossier in the registered cubic R4 reference at step 287."""
import argparse
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from mouth_cp10_label_condition import physical_frame, reference_frame
from mouth_cp10_label_gauge import VERTICES
from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import digest
from mouth_r4_cubic_held import approved_reference
from mouth_r4_eulerian_solid import chart_bounds
from mouth_supercore_shape import preparation


STEP287 = '5683615a04966ef35ba992026c8990965005c13bfec8e6d440534ade1e4b46c0'


class Cell:
    """Physical chart, stored metric and exact Q1 shape functions on this cell."""
    def __init__(self, body, labels, order, key, reference):
        self.body = body
        self.reference = reference
        self.axes = [body.buffers['field_' + name][i:i+2] for name, i in zip('xta', key)]
        self.lo = np.array([axis[0] for axis in self.axes])
        self.width = np.array([axis[1] - axis[0] for axis in self.axes])
        self.vertices = VERTICES + np.asarray(key)
        self.xi = labels[tuple(self.vertices.T)]
        self.order = order[tuple(self.vertices.T)]
        self.chart_labels = np.ascontiguousarray(self.xi.reshape(2,2,2,3))

    def values(self, local):
        q = np.atleast_2d(local)
        weights = np.where(VERTICES, q[:,None,:], 1 - q[:,None,:])
        shape = np.prod(weights, axis=2)
        gradients = np.stack([(2*VERTICES[:,j]-1)*weights[:,:,(j+1)%3]*weights[:,:,(j+2)%3]
                              / self.width[j] for j in range(3)], axis=2)
        # Same translation cancellation as the native whole-cell orientation kernel.
        dx = self.xi - self.xi[:1]
        jac = np.einsum('pvd,vc->pcd', gradients, dx)
        chart = np.linalg.det(jac)
        xi = self.xi[0] + shape @ dx
        pos = self.lo + q*self.width
        spatial, hoop = physical_frame(self.body,pos)
        natural, hoop0 = reference_frame(self.body,self.reference,xi)
        if np.any(hoop <= 0) or np.any(hoop0 <= 0):
            raise ValueError('The implicated cells must be away from the coordinate axes')
        physical = chart * (np.abs(np.linalg.det(natural))*hoop0 / (np.linalg.det(spatial)*hoop))
        return chart, physical, shape@self.order, spatial[:,2,2], pos

    def witness(self, local):
        chart, physical, chi, r, pos = self.values(local)
        return dict(local=np.asarray(local).tolist(), chart_det=float(chart[0]),
                    physical_Jinv=float(physical[0]), chi=float(chi[0]), r=float(r[0]),
                    w=float(pos[0,1]*self.body.data.length), theta=float(pos[0,2]))

    def measured_minimum(self, points):
        grid = np.array(list(product(np.linspace(0,1,9), repeat=3)))
        starts = np.vstack((grid, np.array(points).reshape(-1,3)))
        chart, phys, *_ = self.values(starts)
        def measure(col, initial):
            best = starts[np.argmin(initial)]
            low = float(np.min(initial))
            # Start from several low boundary/edge points as well as the witnesses.
            for x in starts[np.argsort(initial)[:8]]:
                result = minimize(lambda q: self.values(q)[col][0], x, method='L-BFGS-B',
                                  bounds=[(0.,1.)]*3,
                                  options=dict(ftol=1e-15, gtol=1e-11, maxiter=120))
                if np.isfinite(result.fun) and result.fun<low:
                    best, low = result.x, float(result.fun)
            return self.witness(best)
        return measure(0, chart), measure(1, phys)


def site_values(body, labels, order, sites):
    reference = preparation(4)[0]
    return [dict(cell=s['cell'], **Cell(body,labels,order,s['cell'],reference).witness(s['local']))
            for s in sites]


def dossier(body, labels, order, sites, certificate):
    """Give cell bounds and labelled witness values; certificate must apply to these labels."""
    reference = preparation(4)[0]
    keys = sorted({tuple(s['cell']) for s in sites})
    rows = []
    vertex_sets = []
    for key in keys:
        cell = Cell(body,labels,order,key,reference)
        selected = [s for s in sites if tuple(s['cell'])==key]
        pts = [s['local'] for s in selected]
        chart_point, phys_point = cell.measured_minimum(pts)
        chart = chart_bounds(cell.axes,cell.chart_labels,np.ones((1,1,1),np.uint8))
        lower = chart['chart_det_global_lower']
        # Fall back to the previously enclosed global factor on the core; three
        # downward roundings protect the two multiplications at binary scale.
        core = np.nextafter(certificate['core_metric_lower']*certificate['global_xi_r_over_x_lower'],-np.inf)
        core = np.nextafter(core*certificate['global_xi_r_over_x_lower'],-np.inf)
        factor = min(([core] if min(cell.xi[:,0])<=body.cut else [])
                     + ([certificate['collar_metric_lower']] if max(cell.xi[:,0])>body.cut else []))
        phys_lower = np.nextafter(lower*factor,-np.inf)
        phys_lower = np.nextafter(phys_lower*certificate['global_sine_ratio_lower'],-np.inf)
        phys_lower = max(phys_lower,certificate['physical_global_lower'])
        if not 0 < phys_lower <= phys_point['physical_Jinv'] or not lower<=chart_point['chart_det']:
            raise ArithmeticError('Certified lower and measured upper bounds do not enclose the minimum')
        vertices = []
        for node,c in zip(cell.vertices,cell.order):
            neighborhood = tuple(slice(max(0,i-1),i+2) for i in node)
            vertices.append(dict(index=node.tolist(),chi=float(c),
                                 incident_to_ordered_element=bool(np.any(order[neighborhood]>0))))
        vertex_sets.append(set(tuple(v['index']) for v in vertices))
        rows.append(dict(cell=list(key),vertices=vertices,chi_range=[float(cell.order.min()),float(cell.order.max())],
                         all_vertices_melted=bool(np.all(cell.order==0)),
                         certified_chart_lower=lower,certified_physical_lower=float(phys_lower),
                         estimated_chart_minimum=chart_point, estimated_physical_minimum=phys_point,
                         witnesses=[{**s,**cell.witness(s['local'])} for s in selected]))
    shared = [dict(cells=[list(keys[i]),list(keys[j])],shared_vertices=sorted(vertex_sets[i]&vertex_sets[j]))
              for i,j in combinations(range(len(keys)),2)]
    return dict(cells=rows,shared=shared,
                minima='Certified positive lower bounds per cell. Numerical minima use a 9x9x9 closed-cell grid, given witnesses and bounded multistart minimization. A measured value upper-bounds the actual cell minimum up to floating-point roundoff.')


def supplement(data, body, labels, order, failure):
    joint_path = failure.with_name('result.json')
    capture = json.loads(joint_path.read_bytes())
    if capture['accepted_steps']!=0 or data['rejected_trial']['failure_sha256']!=digest(failure):
        raise ValueError('Unmatched joint-probe or failure record')
    retained = data['rejected_trial']['sites']
    extra = [dict(step=288,role='joint_probe_added',**s) for s in capture['joint']['constraints']
             if not any(v['cell']==s['cell'] and v['local']==s['local'] for v in retained)]
    if not extra:
        raise ValueError('No additional joint-probe sites to examine')
    raw_path = failure.with_name('rejected_labels.npz')
    if digest(raw_path)!=data['rejected_trial']['label_sha256']:
        raise ValueError('Original rejected trial changed')
    with np.load(raw_path,allow_pickle=False) as raw:
        at_raw = site_values(body,raw['labels'],order,extra)
    for s,accepted,rejected in zip(extra,site_values(body,labels,order,extra),at_raw):
        cell = next(c for c in data['cells'] if c['cell']==s['cell'])
        cell['witnesses'].append({**s,**accepted})
        data['rejected_trial']['sites'].append({**s,**rejected})
    data.update(joint_result_sha256=digest(joint_path),supplemental_sites=extra)
    return data


def run(root, output, failure, saved_dossier):
    path = root/'r4_cubic_held01/primary/step_00287.npz'
    record = json.loads(path.with_suffix('.json').read_text())
    if digest(path)!=STEP287 or record['sha256']!=STEP287:
        raise ValueError('Wrong accepted R4 state')
    boundary = root/'r4_cubic_held01/stop_readout/result.json'
    proof = json.loads(boundary.read_text())
    if proof['source_sha256']!=STEP287:
        raise ValueError('Physical certificate uses a different state')
    body,cp,_ = get_cp10(root)
    with cp, np.load(path,allow_pickle=False) as z:
        if int(z['radius'])!=4 or not np.array_equal(z['x'],cp['x']):
            raise ValueError('Wrong R4 geometry')
        labels = z['labels']
        order = np.ascontiguousarray(z['field'][2,body.ids].transpose(0,2,1))
    approved_reference(body)
    sites = [dict(step=287,alpha=r['alpha'],cell=r['site']['cell'],local=r['site']['chart_fractions'])
             for r in record['search']['rejections'] if 'chart_fractions' in r['site']]
    capture = None
    if failure is not None:
        capture = json.loads(failure.read_bytes())
        sites += [dict(step=288,role='retained',**s) for s in capture['retained_sites']]
        sites += [dict(step=288,role=k,**capture[k]) for k in ('newly_witnessed','failed_row_site')]
    output.parent.mkdir(parents=True,exist_ok=True)
    if saved_dossier is not None:
        data = json.loads(saved_dossier.read_bytes())
        if failure is None or data['accepted_state_sha256']!=STEP287:
            raise ValueError('A checked step-287 dossier and a failure record are required')
        data = supplement(data,body,labels,order,failure)
        data['original_dossier_sha256'] = digest(saved_dossier)
    else:
        data = dossier(body,labels,order,sites,proof['physical_certificate'])
    data.update(accepted_state_sha256=STEP287,script_sha256=digest(Path(__file__)),case='accepted step 287')
    if capture is not None and saved_dossier is None:
        rejected = failure.with_name('rejected_labels.npz')
        reference = preparation(4)[0]
        with np.load(rejected,allow_pickle=False) as trial:
            raw = trial['labels']
        observations = site_values(body,raw,order,sites)
        witness = site_values(body,raw,order,[capture['newly_witnessed']])[0]
        hit = capture['rejections'][-1]['site']
        if hit['cell']!=witness['cell'] or hit['chart_fractions']!=witness['local']:
            raise ValueError('Raw rejected label witness does not match the stored rejection')
        if not np.isclose(witness['chart_det'],hit['chart_det'],rtol=1e-8,atol=1e-13):
            raise ArithmeticError('Independent rejected-label point evaluation disagrees with production chart')
        samples = []
        for c in data['cells']:
            cell = Cell(body,raw,order,c['cell'],reference)
            chart_min,physical_min = cell.measured_minimum([w['local'] for w in c['witnesses']])
            samples.append(dict(cell=c['cell'], estimated_chart_minimum=chart_min,
                                estimated_physical_minimum=physical_min))
        data['rejected_trial']=dict(label_sha256=digest(rejected),failure_sha256=digest(failure),
                                    alpha=capture['trial_fraction'],new_witness=witness,
                                    production_chart_det=hit['chart_det'],
                                    sites=[{**s,**obs} for s,obs in zip(sites,observations)],
                                    sampled_cell_minima=samples)
    with output.open('x') as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'cells':[(r['cell'],r['certified_physical_lower'],r['estimated_physical_minimum']['physical_Jinv'],
                                r['chi_range']) for r in data['cells']], 'shared':data['shared']},allow_nan=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--failure',type=Path,help='Step-288 failure.json next to the rejected labels')
    parser.add_argument('--supplement',type=Path,help='Augment existing signed dossier with new joint-probe witnesses')
    args=parser.parse_args();run(args.root,args.output,args.failure,args.supplement)
