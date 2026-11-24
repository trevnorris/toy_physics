"""Fixed-field, unit-gauge precision probes on the last R4 melt-only iterate.

The old L-BFGS history was not checkpointed. Use the actual archived iterate
and the solver's specified history-free preconditioner. Compare the old and
extended-precision energies on the very same rounded proposals, with a separate
unrounded path to distinguish coordinate precision from energy arithmetic.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import mpmath as mp
import numpy as np


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def triple(a, b, c):
    return np.sum(a * np.cross(b, c), axis=-1)


class CellWork:
    def __init__(self, barrier, labels):
        self.b = barrier
        self.edges = self.columns(labels)
        self.dets = triple(*self.edges)
        if self.dets.min() <= 0:
            raise ArithmeticError('Archived state has a nonpositive corner')
        self.K = barrier.opt.K.astype(np.longdouble)
        self.x = labels.reshape(-1, 3).astype(np.longdouble) - barrier.opt.base.reshape(-1, 3)
        self.linear = (self.K @ self.x + self.K.T @ self.x) / 2

    def columns(self, field):
        v = field.reshape(-1, 3).astype(np.longdouble)[self.b.ids]
        ids = np.arange(8)
        return tuple(((v[:, ids ^ bit] - v) * np.where(ids & bit, -1, 1)[None,:,None]
                      / self.b.width[:,axis,None,None]) for axis, bit in enumerate((4,2,1)))

    def path(self, direction):
        a, b, c = self.edges
        u, v, w = self.columns(direction)
        d1 = triple(u,b,c) + triple(a,v,c) + triple(a,b,w)
        d2 = triple(u,v,c) + triple(u,b,w) + triple(a,v,w)
        d3 = triple(u,v,w)
        p = direction.reshape(-1,3).astype(np.longdouble)
        hl = np.sum(p*self.linear)
        hq = np.sum(p*(self.K @ p))/2
        def delta(alpha):
            t = np.longdouble(alpha)
            dj = t*(d1+t*(d2+t*d3))
            r = dj/self.dets
            if np.min(r) <= -1:
                raise ArithmeticError('Precision probe left feasible domain')
            # This is the scalar difference, evaluated without subtracting two
            # nearly identical determinants or gauge energies.
            correction = np.empty_like(r)
            small = abs(r) < 1e-3
            z = r[small]
            correction[small] = z*z*(1/np.longdouble(2)+z*(-1/np.longdouble(3)+z*(
                1/np.longdouble(4)+z*(-1/np.longdouble(5)+z*(1/np.longdouble(6)
                - z/np.longdouble(7))))))
            correction[~small] = r[~small]-np.log1p(r[~small])
            return t*(hl+t*hq)+np.sum(self.b.weight[:,None]*(dj*(1-1/self.dets)+correction))
        return delta


def coordinate_fd(barrier, labels, c, group, step):
    """Independent 75-decimal scalar energy evaluated directly on affected cells."""
    cm = barrier.components[c]
    changed = cm.groups == group
    rows = np.flatnonzero(changed)
    ci = np.flatnonzero(changed[barrier.ids].any(axis=1))
    K = barrier.opt.K.tocsr()
    x = labels.reshape(-1,3)
    base = barrier.opt.base.reshape(-1,3)
    m = lambda z: mp.mpf(float(z))
    # Harmonic central difference is confined to entries incident to rows;
    # coalesce a pair with both indices in the set exactly once.
    entries = {}
    for i in rows:
        for j,v in zip(K.getrow(i).indices, K.getrow(i).data):
            entries[i,int(j)] = m(v)
        col = K.getcol(i).tocoo()
        for j,v in zip(col.row,col.data):
            entries[int(j),i] = m(v)
    dx = {i:m(x[i,c])-m(base[i,c]) for ij in entries for i in ij}
    def scalar(h):
        energy = mp.fsum(v*(dx[i]+int(changed[i])*h)*(dx[j]+int(changed[j])*h)/2
                         for (i,j),v in entries.items())
        for ix in ci:
            ids = barrier.ids[ix]
            w = [m(v) for v in barrier.width[ix]]
            points = [[m(v) for v in x[i]] for i in ids]
            for k,i in enumerate(ids):
                if changed[i]:
                    points[k][c] += h
            for k in range(8):
                mat = mp.matrix(3,3)
                for axis, bit in enumerate((4,2,1)):
                    sign = -1 if k & bit else 1
                    for j in range(3):
                        mat[j,axis] = sign*(points[k^bit][j]-points[k][j])/w[axis]
                det = mp.det(mat)
                energy += m(barrier.weight[ix])*(det-1-mp.log(det))
        return energy
    return (scalar(step)-scalar(-step))/(2*step)


def sweep(barrier, extended, labels, step, slope, old_det):
    """Central derivative and work, first along the ideal ray and then at the exact rounded trial."""
    if not slope < 0:
        raise ArithmeticError('Supplied direction is not a strict descent at the archived iterate')
    opt = barrier.opt
    axes = opt.axes
    from mouth_r4_eulerian_relax_labels import chart_orientation
    ideal = extended.path(step)
    fd = []
    for k in (0,4,8,12):
        h = 2.**-k
        estimate = (ideal(h)-ideal(-h))/(2*h)
        fd.append(dict(step=h,central=float(estimate),error=float(estimate-slope)))
    profile = []
    for k in range(32):
        alpha = 2.**-k
        trial = labels + alpha*step
        # Reproduce the recorded legacy Armijo expression, including double
        # corner determinants and fp64 log1p before long-double summation.
        _, _, _, det = barrier.sample(trial)
        dr = trial-labels
        dv = np.sum(barrier.weight[:,None]*((det-old_det)-np.log1p((det-old_det)/old_det)),
                    dtype=np.longdouble)
        mid = ((labels-opt.base)+.5*dr).reshape(-1,3)
        dv += np.sum(dr.reshape(-1,3)*(opt.K @ mid),dtype=np.longdouble)
        actual = extended.path(dr)(1)
        line = alpha*slope
        item = dict(alpha=alpha,linear_slope=line,legacy=float(dv),extended_rounded=float(actual),
                    extended_ideal=float(ideal(alpha)),changed_components=int(np.count_nonzero(dr)),
                    legacy_armijo=bool(dv < 1e-4*line),extended_armijo=bool(actual < 1e-4*line))
        if item['extended_armijo'] and not any(r.get('chart_global_lower') is not None for r in profile):
            item['chart_global_lower'] = chart_orientation(axes,trial,opt.active)['chart_det_global_lower']
        profile.append(item)
    return dict(slope=slope, central_differences=fd, alpha_profile=profile,
                first_extended_descent=next((r for r in profile if r['extended_armijo']),None),
                old_false_rejections=sum(not r['legacy_armijo'] and r['extended_armijo'] for r in profile),
                old_false_accepts=sum(r['legacy_armijo'] and not r['extended_armijo'] for r in profile))


def run(area, out, historical_direction):
    pin = area/'pin/software/stage1_solver/tools'
    sys.path.insert(0,str(pin.resolve()))
    from mouth_cp10_label_grid import get_cp10
    from mouth_r4_eulerian_relax_labels import LabelRelaxation
    from mouth_r4_melt_barrier import MeltBarrier
    from mouth_r4_stage2_outcome import verify_guard

    verify_guard(area/'primary_guard',1)
    folder = area/'primary'
    metadata = json.loads((folder/'failed_trial.json').read_bytes())
    state = folder/metadata['file']
    journal = folder/'inner_0000034.jsonl'
    last = json.loads(journal.read_bytes().splitlines()[-1])
    if digest(state) != metadata['sha256'] or metadata['inner']['residual'] != last['melt_residual']:
        raise ValueError('Wrong inner endpoint')
    b, cp, axes = get_cp10(area.parent)
    with np.load(state, allow_pickle=False) as f:
        if int(f['radius']) != 4 or int(f['winding']) != 1 or not np.array_equal(f['x'],cp['x']):
            raise ValueError('Wrong R4 held geometry')
        labels = f['labels']
        order = np.ascontiguousarray(f['field'][2,b.ids].transpose(0,2,1))
        archive_keys = f.files
    opt = LabelRelaxation(b,labels,order,1e-8)
    with TemporaryDirectory() as tmp:
        barrier = MeltBarrier(b,opt,Path(tmp))
        native_energy, gr, parts, old_det = barrier.sample(labels)
        g = barrier.extract(gr,'dead')
        if abs(g).max() != last['melt_residual']:
            raise ArithmeticError('Production inner gradient does not reproduce the stall')
        extended = CellWork(barrier,labels)
        fresh = -barrier.pre(g)
        step = barrier.melt_direction(fresh)
        scale = max(1.,*(abs(step[...,k]).max()/(.25*np.median(np.diff(ax))) for k,ax in enumerate(axes)))
        fresh /= scale
        step /= scale
        slope = float(g @ fresh)
        if slope >= 0:
            raise ArithmeticError('History-free direction fails descent')
        coords = []
        with mp.workdps(75):
            for c,cm in enumerate(barrier.components):
                nodal = cm.sum(gr[...,c])
                groups = np.flatnonzero(cm.dead)
                group = groups[np.argmax(abs(nodal[groups]))]
                derivative = coordinate_fd(barrier, labels,c,group,mp.mpf('1e-22'))
                coords.append(dict(component=c,group=int(group),
                    nodes=np.argwhere(cm.groups.reshape(order.shape)==group).tolist(),
                    signed_analytic=float(nodal[group]), fd_75_digit=mp.nstr(derivative,34),
                    signed_error=float(mp.mpf(float(nodal[group]))-derivative)))
        fresh_result = sweep(barrier, extended, labels, step, slope, old_det)
        historic_result = None
        if historical_direction is not None:
            reconstruction_path = historical_direction.with_name('reconstruction.json')
            reconstruction = json.loads(reconstruction_path.read_bytes())
            if (reconstruction['direction_sha256'] != digest(historical_direction)
                    or reconstruction['inner_journal_sha256'] != digest(journal)
                    or reconstruction['journal_rows_all_scalars_equal'] != 106
                    or not reconstruction['final_labels_bit_exact']):
                raise ValueError('Uncertified secant history')
            verify_guard(historical_direction.parent.parent/'direction_guard',0)
            with np.load(historical_direction,allow_pickle=False) as vectors:
                expected = sha256(np.ascontiguousarray(labels).tobytes()).hexdigest()
                if str(vectors['labels_sha256']) != expected:
                    raise ValueError('Historical direction belongs to a different map')
                packed = vectors['packed']
                historic_step = vectors['step']
                if packed.shape != g.shape or not np.array_equal(barrier.melt_direction(packed),historic_step):
                    raise ValueError('Historical direction violates the registered melt-only reduction')
            historic_result = sweep(barrier, extended, labels, historic_step, float(g @ packed), old_det)
            historic_result['direction_sha256'] = digest(historical_direction)
            if historic_result['slope'] != reconstruction['slope']:
                raise ValueError('Reconstructed analytic slope differs')
            historic_result['fixed_path_certification'] = reconstruction
    frozen = json.loads((area/'sources.json').read_bytes())
    modules = sorted({Path(m.__file__).name for m in sys.modules.values()
              if getattr(m,'__file__',None) and Path(m.__file__).parent == pin.resolve()} | {'mouth_r4_melt_barrier.cpp'})
    hashes = {k:digest(pin/k) for k in modules}
    if any(frozen[k]!=v for k,v in hashes.items()):
        raise ValueError('Pinned input source changed')
    result = dict(command=sys.orig_argv,script_sha256=digest(Path(__file__)),
                  inner_state_sha256=digest(state),journal_sha256=digest(journal),
                  pinned_source_sha256=hashes,keys_in_inner_archive=archive_keys,
                  historical_direction=('Certified reconstruction supplied' if historic_result is not None else
                    'Not saved or reconstructed. Secant vectors are absent. Saved direction files are outer proposals.'),
                  previous_inner=metadata['inner'], final_accepted=last, all_active_and_order_fixed=True,
                  normalized_gauge_energy=float(native_energy),native_components=parts,
                  extended_mantissa_bits=int(np.finfo(np.longdouble).nmant)+1,coords=coords,
                  history_free=fresh_result,historical=historic_result)
    with out.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    brief = lambda r: None if r is None else {k:v for k,v in r.items() if k != 'alpha_profile'}
    print(json.dumps(dict(coords=coords,fresh=brief(fresh_result),historical=brief(historic_result))),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--historical-direction',type=Path)
    args = parser.parse_args()
    run(args.area.resolve(),args.output,args.historical_direction)
