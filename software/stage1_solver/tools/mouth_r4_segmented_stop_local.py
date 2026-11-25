"""Read the last held R4 state at the step-167 trial witnesses and the entrance force peak."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def read_state(folder, name):
    path = folder/(name+'.npz')
    record = json.loads(path.with_suffix('.json').read_bytes())
    if sha(path)!=record['sha256']:
        raise ValueError('Changed accepted archive '+name)
    with np.load(path,allow_pickle=False) as data:
        if (int(data['radius']),int(data['winding']),float(data['epsilon']))!=(4,1,1e-8):
            raise ValueError('Unexpected preparation')
        return record,{k:data[k] for k in ('field','labels','x')}


def main(study, output):
    pin = study/'held_pin'
    src = pin/'software/stage1_solver/tools'
    if any(sha(src/n)!=s for n,s in json.loads((pin/'manifest.json').read_bytes()).items()):
        raise ValueError('Held source pin mismatch')
    sys.path.insert(0,str(src))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_order import KAPPA, WELL
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_r4_aligned_rematch import front
    from mouth_r4_eulerian_corner_audit import point
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_refinement_review import _collar_reference
    from mouth_r4_slaved_melt import FLOW
    from mouth_r4_smooth_fields import GaugeGeometry, SmoothOrder, material_load
    from mouth_r4_smooth_held import SmoothHeld
    folder = study/'primary'
    fresh,prev = read_state(folder,'matched_00148')
    last,st = read_state(folder,'step_00166')
    if not np.array_equal(st['field'],prev['field']) or not np.array_equal(st['x'],prev['x']):
        raise ValueError('Final burst did not retain the saved matched fields and fixed shape')
    issue = json.loads((folder/'ISSUE.json').read_bytes())['witness']
    if issue['material']!=last['material']:
        raise ValueError('Issue does not arise at last accepted material state')
    # Explicit static assembly: no sparse preconditioner or relaxation step is needed.
    driver = SmoothHeld(study.parent,output.parent,1e-8)
    if not np.array_equal(driver.x,st['x']):
        raise ValueError('Held shape differs')
    trial = ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
        axes=axes_for_depth(driver.shape(driver.x).length,driver.reference_depth))
    trial.model.extended_energy_sum = True
    trial.diag = (sum(g.power(2).T@(KAPPA*trial.model.hq) for g in trial.G)
                 +trial.V.power(2).T@(2*WELL*trial.model.hq))
    trial.smooth_order_geometry = GaugeGeometry(trial,aligned=True)
    solid = SpatialCrystal(trial.b)
    u,labels = st['field'],st['labels']
    order = order_grid(trial,u)
    load,err = material_load(trial,solid,labels,u)
    trial.model.material_load = load
    phys,g = trial.model.evaluate(u)
    chi = u[2,trial.first]
    unit,greg = SmoothOrder(trial,labels,driver.epsilon).gauge(chi)
    joined = trial.pull.T@g[2]+driver.epsilon*greg
    kkt = np.where(chi==0,np.minimum(joined,0),np.where(chi==1,np.maximum(joined,0),joined))
    node = np.flatnonzero(trial.free)[np.argmax(abs(kkt[trial.free]))]
    current = trial.measure(u,load)
    field = dict(order_kkt=float(abs(kkt[node])),order_kkt_node=int(node),order_at_peak=float(chi[node]),
                 signed_order_derivative=float(joined[node]),order_kkt_node_depth=float(trial.depth[node]),
                 order_projected=float(np.max(abs(trial.projected(chi,joined)))),
                 amplitude_raw=float(np.max(abs(g[0,trial.model.masks[0]]))),
                 phase_raw=float(np.max(abs(g[1,trial.model.masks[1]]))))
    archive = [json.loads(line) for line in (study/'primary_guard/stdout').read_text().splitlines()]
    fresh_live = next(row for row in archive if row['kind']=='held_fresh_fields' and row['step']==148)
    if (any(float(current[k])!=fresh_live['current'][k] for k in FLOW)
            or any(field[k]!=fresh['field'][k] for k in ('amplitude_raw','phase_raw'))
            or abs(float(phys+driver.epsilon*unit)-last['total_energy'])>64*np.spacing(float(phys))):
        raise ArithmeticError('Saved-state current, carrier forces or energy failed held-state readback')
    axis = solid.axes
    def local(cell,f):
        with mp.workdps(64):
            return point(trial.b,axis,labels,order,cell,f,_collar_reference)
    inversion = [r for r in issue['tests'] if r['site'].get('code')==2]
    peak = last['material']['worst_free_active']
    by_vertex = np.array(peak['node'])
    one_sided=[]
    for bits in np.ndindex(2,2,2):
        cell=by_vertex-np.array(bits)
        if np.all(cell>=0) and np.all(cell<np.array(order.shape)-1):
            one_sided.append(local(cell,bits))
    sites = [local(row['site']['cell'],row['site']['chart_fractions']) for row in inversion]
    report = dict(command=sys.orig_argv,script_sha256=sha(Path(__file__)),pin_sha256=sha(pin/'manifest.json'),
         saved_states={'matched_148':fresh['sha256'],'accepted_166':last['sha256']},
         source_issue_sha256=sha(folder/'ISSUE.json'),fixed_shape=st['x'].tolist(),field_arrays_identical=True,
         field_at_last_match=fresh_live,unrematched_fields_166=field,total_energy=float(phys+driver.epsilon*unit),
         gauge_energy=float(driver.epsilon*unit),physical_energy=float(phys),shared_work_error=float(err),
         flow_at_166={k:float(current[k]) for k in FLOW},front=front(trial,u),
         accepted_peak=peak,accepted_peak_incident_elements=one_sided,
         failed_inversion_trials=inversion,accepted_state_at_inversion_witness=sites,
         terminal_trials=issue['tests'],global_orientation=last['material']['chart'],
         energy_convention='W_formal is the crystal law at the label map; chi_Jinv_W is its physical-volume energy density.')
    with output.open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(dict(step166_fields=field,accepted_witness=sites,peak_force=peak),allow_nan=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('study',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    main(args.study.resolve(),args.output.resolve())
