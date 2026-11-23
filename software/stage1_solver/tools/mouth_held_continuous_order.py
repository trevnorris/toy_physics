"""Continuously matched R4 preparation on the saved fine lattice under rest and drive."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def read(p):
    return json.loads(p.read_text())


def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(folder, name, **fields):
    p = folder/(name+'.partial.npz')
    np.savez_compressed(p, **fields)
    p.replace(folder/(name+'.npz'))


def write(p, data):
    tmp = p.with_suffix('.partial.json')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    tmp.replace(p)


def run(args):
    from mouth_continuous_order import ContinuousTrial, seed_order, norm
    from mouth_joined_front_audit import Front
    from mouth_last_rings_fluid import GRIDS
    from mouth_poly_shallow_release import CUT
    from mouth_receiver_order_equilibrium import receiver
    from mouth_receiver_order_control import side_source
    from mouth_rim_pinned_comparison import validate_guard
    src = args.original
    folder = src/'fine'
    validate_guard(folder)
    endpoint = read(folder/'endpoint.json')
    audit = read(args.edges)['results']['fine']
    state = Path(read(folder/'invocation.json')['command'][5])
    cp = folder/endpoint['file']
    controls = folder/'endpoint_controls.npz'
    for path,key in ((state,'preparation_sha256'),(cp,'state_sha256'),(controls,'controls_sha256'),
                     (folder/'endpoint.json','endpoint_sha256')):
        if digest(path)!=audit[key]:
            raise ValueError(f'Changed fine input {path}')
    with np.load(cp) as f, np.load(controls) as c:
        x, q, narrow = f['x'], f['q'], dict(wound=f['field'], rest=c['rest'])
    base = Front(state, GRIDS['double'], 1)
    opt = audit['solid_discretization']
    t = ContinuousTrial(base, state, x, **opt)
    source = side_source(base,t.shape)
    profile, _ = receiver(t.shape, base.order_mesh, CUT, lambda w,r: source,
                          polar_shift=t.polar_shift, material_path=state)
    m = t.model
    def init(name):
        return t.transported(seed_order(m, base.c, profile, narrow[name]))
    initial = init('rest')
    args.saved.mkdir(parents=True,exist_ok=True)
    results = dict(state=digest(state), fine=digest(cp),controls=digest(controls),
                   saved_original_pin=digest(src/'sources.json'), code=digest(Path(__file__)),
                   model=digest(Path(__file__).with_name('mouth_continuous_order.py')),
                   generated_by=sys.argv, solid=opt, coordinates=x.tolist(), field_grid=list(base.grid),
                   carrier_nodes=m.size, absolute_order_nodes=len(t.free),free_order_nodes=int(sum(t.free)),
                   far_dirichlet_nodes=len(t.fixed_order), exterior_crystal_unit_order=float(t.ext.sum()),
                   solutions={})
    for name in ('rest','wound'):
        if name=='wound':
            initial[:2]=narrow['wound'][:2]
        iteration_offset = [0]
        def progress(iteration, e, resid, c, scale, damping):
            iteration_offset[0] += 1
            total = iteration_offset[0]
            if total == 1 or total%10 == 0:
                save(args.saved,'intermediate_'+name, x=x,q=q,absolute_order=c,iteration=total)
                print(json.dumps(dict(field=name,order_iteration=total,energy=e,
                                      previous_projected_residual=resid,step=scale,regularization=damping)),flush=True)
        u, energy, gq, load, meta = t.match(q,initial,rest=name=='rest',log=progress)
        save(args.saved,name, x=x,q=q,field=u,load=load)
        current = m.measure(u)
        c = u[2,t.first]
        old_inner = m.ids[1,0]
        old_d = np.broadcast_to(m.z[1][None,None,:],old_inner.shape)
        old_r = np.broadcast_to(m.rr[0][:,None,None],old_inner.shape)
        original = RegularInner(base.c,old_d/t.shape.length,old_r)+narrow[name][2,old_inner]
        result = dict(sha256=digest(args.saved/(name+'.npz')), meta=meta,energy=energy,
                      solid=t.b.residual(q,gq), current=current, order_range=[float(c.min()),float(c.max())],
                      in_shaft_axis=list(map(float,u[2,old_inner[0,0,[0,len(m.z[1])//2,-1]]])),
                      prior_high_order_nodes=int(np.count_nonzero(original>.9)),
                      now_low_among_prior_high=int(np.count_nonzero((original>.9)&(u[2,old_inner]<.1))),
                      absolute_trace_max=norm(u[2]-c[t.number]), lower_obstacle_count=int(sum(c[t.free]==0)),
                      upper_obstacle_count=int(sum(c[t.free]==1)))
        results['solutions'][name] = result
        write(args.saved/'record.json',results)
        initial = u.copy()
        print(json.dumps(dict(finished=name,energy=energy,meta=meta,ring=result['solid']['last_free_ring_angular_sum'])),flush=True)
        if (c.min()<0 or c.max()>1 or result['absolute_trace_max']>1.e-13
                or meta['field']['stationarity']>=2.e-7 or current['port_error']>=3.e-5):
            raise ArithmeticError(f'{name} field check failed')


def RegularInner(values,t,x):
    from scipy.interpolate import RegularGridInterpolator
    axes = tuple(np.linspace(0,1,n) for n in values.shape)
    return RegularGridInterpolator(axes,values)((t,x))


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('original','edges','saved'):
        p.add_argument(key,type=Path)
    args = p.parse_args()
    pin = args.original/'pin'
    for path,sha in read(args.original/'sources.json').items():
        if digest(pin/path)!=sha:
            raise ValueError(f'Changed original pinned source {path}')
    sys.path.insert(1,str(pin/'software/stage1_solver/tools'))
    run(args)
