"""Read-only reconstruction of a failed R4 crystal projection and a joint-direction test."""
import argparse
import itertools
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import nnls

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_aligned_crystal import AlignedCrystal
from mouth_r4_collar166_state import load_held
from mouth_r4_cubic_held import approved_reference
from mouth_r4_cubic_physical_orientation import physical_certificate
from mouth_r4_eulerian_release import order_grid
from mouth_r4_eulerian_relax_labels import InfeasibleLabels
from mouth_r4_eulerian_solid import InvalidSpatialCrystal
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


def saved(folder, stem):
    record = read(folder / (stem + '.json'))
    path = folder / record['file']
    if digest(path) != record['sha256']:
        raise ValueError(f'State checksum differs for {stem}')
    with np.load(path, allow_pickle=False) as state:
        if int(state['radius']) != 4 or int(state['winding']) != 1:
            raise ValueError('Incorrect preparation')
        return record, state['labels'], state['field'], state['x']


def key(cell, local):
    return tuple(int(x) for x in cell), tuple(float(x) for x in local)


def site_record(site):
    return dict(cell=list(site[0]), local=list(site[1]))


def witness(error):
    return getattr(error, 'witness', getattr(error, 'diagnostics', {}))


def checked_context(study, out):
    verify_guard(study / 'primary_guard', 1)
    before, l286, fields, shape = saved(study / 'primary', 'matched_00286')
    after, labels, same_fields, same_shape = saved(study / 'primary', 'step_00287')
    driver, trial, solid, _ = load_held(study, out, 1e-8, solve_fields=False)
    approved_reference(trial.b)
    if not (np.array_equal(shape, driver.x) and np.array_equal(shape, same_shape)
            and np.array_equal(fields, same_fields)):
        raise ValueError('Saved order, carriers or held shape differ')
    order = order_grid(trial, fields)
    joint = AlignedCrystal(driver, trial, solid, l286, order)
    value, gradient, material = joint.evaluate(l286)
    if material['max_free_label_force'] != before['material']['max_free_label_force']:
        raise ValueError('Step-286 crystal residual failed reconstruction')
    # Recompute the archived step only in memory to recover the exact cached
    # diagonal, CG forcing history, and retained cell witnesses for step 288.
    candidate, _, _, m287, search = joint.step(l286, value, gradient, material)
    check = dict(max_label_difference=float(np.max(np.abs(candidate-labels))),
                 step_fraction=search['alpha'], iterations=search['cg_iterations'],
                 work=search['energy_change'], backtracks=search['backtracks'],
                 peak_force=m287['max_free_label_force'],
                 matching_rejections=search['rejections']==after['search']['rejections'],
                 source_hashes=[before['sha256'],after['sha256']])
    for field in ('alpha','cg_iterations','energy_change','backtracks','slope','forcing_eta'):
        if search[field] != after['search'][field]:
            raise ValueError(f'Archived step-287 {field} differs: {check}')
    if check['max_label_difference'] != 0 or not check['matching_rejections']:
        raise ValueError(f'Archived step-287 direction differs: {check}')
    return joint, labels, check


def recover(joint, labels, args):
    value, gradient, material = joint.evaluate(labels)
    try:
        joint.step(labels, value, gradient, material)
    except ArithmeticError as exc:
        if str(exc) != 'Unresolved corner projection':
            raise
        levels = {}
        tb = exc.__traceback__
        while tb:
            levels[tb.tb_frame.f_code.co_name] = tb.tb_frame.f_locals.copy()
            tb = tb.tb_next
    else:
        raise ArithmeticError('The recorded projection failure did not reproduce')
    step, project, failure = (levels[k] for k in ('step','project','project_halfspaces'))
    raw = project['raw']
    sites = [key(c,t) for c,t in joint.tight]
    raw_normals = [joint.corner_normal(labels,*s)[0] for s in sites]
    participating = [next(i for i,n in enumerate(raw_normals) if np.array_equal(row,n))
                     for row in failure['normals']]
    normalized, result = failure['N'], failure['result']
    slack = normalized @ result
    # Extended precision measures primal slack independently of Gram truncation.
    nz = np.any(normalized != 0, axis=0)
    n = normalized[:,nz].astype(np.longdouble)
    direct = n @ result[nz].astype(np.longdouble)
    vals = failure['values']
    eig, vec, keep = failure['eig'], failure['vec'], failure['keep']
    old_lambda = failure['multipliers']
    ns = len(normalized)
    dropped_rhs = vec[:,~keep] @ (vec[:,~keep].T @ vals)
    grams = failure['grams']
    new_witness = step['tests'][-1]['site']
    new_site = key(new_witness['cell'],new_witness['chart_fractions'])
    if new_site != sites[-1]:
        raise ValueError('The newly detected corner differs from the projector history')
    failed_site = sites[participating[int(np.argmin(slack))]]
    capture = dict(reproduction=args, direction='Newton' if raw is step['p'] else 'preconditioned gradient',
                   retained_sites=[site_record(s) for s in sites],rows=participating,
                   trial_fraction=step['alpha'],newly_witnessed=site_record(new_site),
                   failed_row_site=site_record(failed_site),rejections=step['tests'],
                   min_normalized_slack=float(min(slack)),normalized_slacks=slack.tolist(),
                   extended_slack=direct.astype(float).tolist(),
                   permitted_negative_slack=-1e-10*max(1.,float(np.linalg.norm(raw))),
                   raw_slacks=vals.tolist(),multipliers=old_lambda.tolist(),
                   eigenvalues=eig.tolist(),threshold=float(eig.max()*ns*np.finfo(float).eps),
                   retained_rank=int(keep.sum()),dropped_rhs_component=dropped_rhs.tolist(),
                   gram_residual=(vals+grams@old_lambda).tolist(), cg_iterations=step['cg'],
                   cg_reason=step['reason'],gradient_dot_raw=float(gradient@raw))
    return capture, raw, gradient, material, sites, step['trial']


def joint_direction(raw, diagonal, normals):
    """Solve the diagonal-metric projection using the local unsquared row operator.

    The passive set from NNLS is polished on the original data in 70-digit arithmetic.
    A small inward margin resolves roundoff when applying the projected solution.
    """
    involved = np.any(normals!=0,axis=0)
    n = normals[:,involved]
    diag = diagonal[involved]
    root = np.sqrt(diag)
    ns = np.linalg.norm(n*root,axis=1)
    a = n*root / ns[:,None]
    y = raw[involved]/root
    guess, _ = nnls(a.T, -y)
    with mp.workdps(70):
        N = mp.matrix(n.tolist())
        D = mp.diag(diag.tolist())
        p = mp.matrix(raw[involved].tolist())
        gram = N*D*N.T
        sizes = [mp.sqrt(gram[i,i]) for i in range(N.rows)]
        A = mp.diag([1/t for t in sizes])*N
        G = A*D*A.T
        b = A*p
        # A common positive margin in the whitened row norm, three orders
        # below our displayed feasibility tests even for the strongest rows.
        margin = mp.mpf('3e-14')*max(1,max(map(abs,b)))
        b = b - mp.matrix([margin]*N.rows)
        seed = tuple(int(i) for i in np.where(guess>0)[0])
        def subsets():
            yield seed
            for count in range(N.rows+1):
                for bank in itertools.combinations(range(N.rows),count):
                    if bank!=seed:
                        yield bank
        for bank in subsets():
            lam = mp.matrix(N.rows,1)
            if bank:
                block = mp.matrix([[G[i,j] for j in bank] for i in bank])
                rhs = mp.matrix([-b[i] for i in bank])
                try:
                    active = mp.lu_solve(block,rhs)
                except (ValueError,ZeroDivisionError):
                    continue
                if min(active)<0:
                    continue
                for i,l in zip(bank,active):
                    lam[i] = l
            slack = b+G*lam
            if min(slack)<-mp.mpf('1e-45')*max(1,max(map(abs,b))):
                continue
            direction = raw.copy()
            direction[involved] = np.array([float(v) for v in p+D*A.T*lam])
            dots = n@direction[involved]
            normal_slacks = dots/np.linalg.norm(n,axis=1)
            if np.any(dots<0):
                raise ArithmeticError(f'Unresolved binary64 joint direction: {normal_slacks}')
            return direction,dict(active_rows=list(bank),slacks=normal_slacks.tolist(),
                                  margin_in_scaled_rows=float(margin),precision_digits=70,
                                  scaled_condition=float(np.linalg.cond(a)), support_dofs=int(involved.sum()))
    raise ArithmeticError('No resolved joint projection under the inward roundoff margin')


def test_joint(joint, labels, raw, gradient, material, sites):
    trials = []
    while True:
        normals = np.array([joint.corner_normal(labels,*s)[0] for s in sites])
        direction, certificate = joint_direction(raw,joint.jacobi,normals)
        slope = float(gradient@direction)
        if slope>=0:
            return dict(feasible_descent=False, reason='Non-descending joint projection',slope=slope,
                        constraints=[site_record(s) for s in sites],projection=certificate,trials=trials)
        alpha = 1.
        new_site = None
        while alpha>1e-4:
            test = labels+alpha*joint.lift(direction)
            try:
                _,_,after = joint.evaluate(test)
                de = after['elastic_energy']-material['elastic_energy']
                dg = float(joint.epsilon*joint.weighted.work(labels,test))
                work = de+dg
                ok = (work<0 and work<=1e-4*alpha*slope and
                      (de!=0 or np.array_equal(test[joint.opt.stencil],labels[joint.opt.stencil])))
                trials.append(dict(alpha=alpha,valid=True,signed_work=work,strict_armijo=bool(ok)))
                if ok:
                    orientation = physical_certificate(joint.trial.b,test,after['chart']['chart_det_global_lower'])
                    return dict(feasible_descent=True,alpha=alpha,slope=slope,physical_work=de,gauge_work=dg,
                                total_work=work,armijo_rhs=1e-4*alpha*slope,peak_force=after['max_free_label_force'],
                                chart=after['chart'],sampled_physical_min=after['min_Jinv_global'],
                                physical_certificate=orientation, constraints=[site_record(s) for s in sites],
                                projection=certificate,trials=trials)
            except (InfeasibleLabels,InvalidSpatialCrystal) as error:
                info = witness(error)
                trials.append(dict(alpha=alpha,valid=False,reason=str(error),site=info))
                where = info.get('chart_fractions',info.get('local'))
                if where is not None and 'cell' in info:
                    pair = key(info['cell'],where)
                    if pair not in sites:
                        new_site = pair
                        break
            alpha *= .5
        if new_site is None:
            return dict(feasible_descent=False,reason='No verified descent at production fractions',
                        slope=slope,projection=certificate,constraints=[site_record(s) for s in sites],trials=trials)
        sites.append(new_site)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output',type=Path)
    options=parser.parse_args()
    study,out=options.study.resolve(),options.output.resolve()
    out.mkdir(exist_ok=False)
    joint,labels,control = checked_context(study,out)
    np.seterr(invalid='raise',over='raise')
    capture,raw,g,material,sites,rejected = recover(joint,labels,control)
    with (out/'rejected_labels.npz').open('xb') as stream:
        np.savez_compressed(stream,labels=rejected)
    atomic_json(out/'failure.json',capture)
    normals = np.stack([joint.corner_normal(labels,*s)[0] for s in sites])
    used = np.any(normals!=0,axis=0)
    with (out/'local_projection.npz').open('xb') as f:
        np.savez_compressed(f,indices=np.flatnonzero(used),normals=normals[:,used],raw=raw[used],
                            diagonal=joint.jacobi[used],active_rows=capture['rows'])
    answer = test_joint(joint,labels,raw,g,material,sites.copy())
    row=dict(accepted_steps=0,source_commit='f28f283f',epsilon=1e-8,
             step_287_peak_force=material['max_free_label_force'],failure=capture,joint=answer)
    atomic_json(out/'result.json',row)
    print('rejected row',capture['failed_row_site'],'new witness',capture['newly_witnessed'],
          'joint descent',answer['feasible_descent'],'accepted steps',0,flush=True)


if __name__=='__main__':
    main()
