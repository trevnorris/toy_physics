"""Local material and coherent-field response on the last saved, rejected R4 wall.

Hold the exterior shape fixed. Probe nine original material sites around the
peak and take the Schur complement of the field and neighboring displacement
coordinates. This tests the tangent on the native quadrature, including its
one-sided interpolation checks. It does not solve the whole-mouth equilibrium.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize._numdiff import group_columns
from scipy.sparse.linalg import splu

from mouth_bonded_front_adjoint import set_front
from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_ring_solid import source
from mouth_supercore_joined_front import file_sha


def field_tangent(model, field):
    # The three fields share hex element connectivity. Columns sharing no
    # affected equation can be differentiated together by complex step.
    V = model.V.astype(bool)
    adjacency = (V.T @ V + model.B.astype(bool)).tocsr()
    mask = model.masks
    pattern = sparse.bmat([[adjacency[a][:, b] for b in mask] for a in mask], format='csc')
    groups = group_columns(pattern)
    mat = pattern.astype(float)
    zero = np.zeros_like(field)
    def action(direction):
        du = model.unpack(direction, zero)
        return model.vars(model.evaluate(field + 1j*1e-25*du)[1].imag)/1e-25
    for color in np.unique(groups):
        cols = np.flatnonzero(groups == color)
        grad = action((groups == color).astype(float))
        for c in cols:
            loc = slice(mat.indptr[c], mat.indptr[c+1])
            mat.data[loc] = grad[mat.indices[loc]]
    probe = np.random.default_rng(182).normal(size=mat.shape[1])
    direct = action(probe)
    check = np.max(abs(mat @ probe - direct))/np.max(abs(direct))
    skew = sparse.linalg.norm(mat-mat.T)/sparse.linalg.norm(mat)
    if max(check, skew) > 1e-9:
        raise ArithmeticError(f'Independent field tangent checks failed {check} {skew}')
    return mat, dict(columns=len(groups), colors=len(np.unique(groups)),
                     action_relative=check, skew_relative=skew)


def relax_neighbors(matrix, target):
    other = np.setdiff1d(np.arange(len(matrix)), target)
    block = matrix[np.ix_(target, other)]
    return (matrix[np.ix_(target,target)] -
            block @ np.linalg.solve(matrix[np.ix_(other,other)],block.T))


def run(state, checkpoint):
    _, d, u, x = source(state, checkpoint)
    trial = Trial(Front(state,(13,10,14),1), state, x, 11, 32, True)
    b, m = trial.body, trial.model
    set_front(b,m,u)
    sites = [(j,i,c) for j in (29,30,31) for i in (30,31,32) for c in (0,1)]
    index = tuple(np.array(sites).T)
    targets = np.array([sites.index((31,32,c)) for c in (0,1)])
    if not b.mask[index].all():
        raise ValueError('Held material site in probe')
    selected = np.isin(b.id, np.ravel_multi_index(index[:2], d.shape[:2])).any(axis=1)
    disp = np.einsum('ij,ij->i',b.basis,d[:,:,1].ravel()[b.id])
    band = float(((b.w+disp)[selected]/b.L).min())
    if band <= .62 or b.tt[selected].min() <= .62:
        raise ValueError('Angular deformation could modify the local tangent')
    def read(dd):
        q = dd[b.mask]*b.scale
        Bonded.energy(b,q)
        return b.native[index].copy(), b.get_material_load(q,m)
    g0, m.material_load = read(d)
    stationarity = float(np.max(abs(m.norm(m.vars(m.evaluate(u)[1])))))
    if stationarity > 2e-7:
        raise ValueError(f'Field preparation residual {stationarity}')
    Hf, metrics = field_tangent(m,u)
    factor = splu(Hf)
    rows=[]
    for h in (2e-7, 1e-7):
        K = np.empty((len(sites),len(sites)))
        B = np.empty((Hf.shape[0],len(sites)))
        sides = []
        for k, site in enumerate(sites):
            dd=d.copy()
            dd[site]+=h
            gp,lp=read(dd)
            dd[site]-=2*h
            gm,lm=read(dd)
            K[:,k]=(gp-gm)/(2*h)
            load=np.zeros_like(u)
            load[2]=(lp-lm)/(2*h)
            B[:,k]=m.vars(load)
            sides.append(float(np.max(abs((gp-2*g0+gm)/h))))
        W=factor.solve(B)
        resid=np.max(abs(Hf@W-B))/max(np.max(abs(B)),1e-30)
        correction=B.T@W
        raw=float(np.max(abs(K-K.T)))
        K=(K+K.T)/2
        relaxed=K-(correction+correction.T)/2
        data={}
        for name,mat in [('held_field',K),('responsive_field',relaxed)]:
            for kind,thin in [('held_neighbors',mat[np.ix_(targets,targets)]),
                              ('responsive_neighbors',relax_neighbors(mat,targets))]:
                data[f'{name}_{kind}']={'matrix':thin.tolist(),'eigenvalues':np.linalg.eigvalsh(thin).tolist()}
        rows.append(dict(step=h,field_solve_relative=resid,material_skew_absolute=raw,
                         one_sided_max_difference=max(sides),patch_spectrum_min=np.linalg.eigvalsh(relaxed)[0],
                         field_correction_max_abs=float(np.max(abs(correction))),readings=data))
    return dict(source_sha=file_sha(state),checkpoint_sha=file_sha(checkpoint),grid=[33,32],
                local_nodes=[list(s) for s in sites[::2]],active_min_actual_depth=band,
                field_preconditioned_residual=stationarity, signed_peak_force=(-g0[targets]).tolist(),
                field_tangent_checks=metrics,steps=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('checkpoint',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    result=run(args.state,args.checkpoint)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
