"""Independent affine-triangle check of static stress work on the attached annulus."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_sealed_shape_relax import Meridian


def triangle_response(radius, saved_mesh, stride, manifest):
    name = f"r{radius}_broad_{saved_mesh}"
    prep = json.loads(manifest["runs"][name]["literal_stdout"].splitlines()[0])
    path = Path("_scratch/s11c/live_waist_dynamics01/fields") / (name+".npz")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != manifest["saved_field_sha256"][str(path)]:
        raise ValueError("Archived snapshots differ")
    with np.load(path) as archive:
        order = archive["order"][::stride, saved_mesh//2::stride]
    shape = Meridian(radius, prep["coordinates"])
    length = shape.length
    t = np.linspace(0., 1., order.shape[0])
    x = np.linspace(.5, 1., order.shape[1])
    t, x = np.meshgrid(t, x, indexing="ij")
    r = x * shape.radius(length*t)[0]
    xy = np.c_[r.ravel(), (length*t).ravel()]
    index = np.arange(t.size).reshape(t.shape)
    lo, up = index[:-1], index[1:]
    triangles = np.concatenate([np.stack([lo[:, :-1], lo[:, 1:], up[:, 1:]], axis=-1),
                                np.stack([lo[:, :-1], up[:, 1:], up[:, :-1]], axis=-1)]).reshape(-1, 3)
    points = xy[triangles]
    affine = np.concatenate([np.ones((*points.shape[:2], 1)), points], axis=2)
    basis_gradient = np.linalg.inv(affine)[:, 1:, :]
    deriv_r, deriv_w = basis_gradient[:, 0], basis_gradient[:, 1]
    areas = .5*np.linalg.det(affine)
    if np.min(areas) <= 0.:
        raise ValueError("Inverted physical triangle")
    # Four-point cubic rule. The 4D integration weight contains r^2
    # and each hoop strain contributes 1/r, so cubic is exact in each triangle.
    nodes = np.array([[1/3., 1/3., 1/3.], [.6,.2,.2], [.2,.6,.2], [.2,.2,.6]])
    weights = np.array([-27.,25.,25.,25.])/48.
    local = np.zeros((len(points), 6, 6))
    components = []
    for bary, wt in zip(nodes, weights):
        rq, wq = np.einsum("i,cij->jc", bary, points)
        order_q = order.ravel()[triangles] @ bary
        measure = 4.*np.pi*areas*wt*rq**2*order_q
        eps = np.zeros((len(points), 5, 6))
        eps[:, 0, 0::2] = deriv_r
        eps[:, 1, 1::2] = deriv_w
        eps[:, 2, 0::2] = eps[:, 3, 0::2] = bary / rq[:, None]
        eps[:, 4, 0::2], eps[:, 4, 1::2] = .5*deriv_w, .5*deriv_r
        trace = eps[:, :4].sum(axis=1)
        contraction = np.einsum("cij,i,cik->cjk", eps, [2.,2.,2.,2.,4.], eps)
        local += measure[:, None, None]*(contraction+.5*trace[:, :, None]*trace[:, None, :])
        components.append((rq, wq, measure, eps))
    dofs = (2*triangles[:, :, None] + np.arange(2)).reshape(-1, 6)
    rows = np.broadcast_to(dofs[:, :, None], local.shape)
    cols = np.broadcast_to(dofs[:, None, :], local.shape)
    stiffness = sparse.coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                                  shape=(2*t.size, 2*t.size)).tocsr()
    free = np.repeat(((x>.5) & (x<1.) & (t>0.) & (t<1.)).ravel(), 2)
    # Positive 0.001 offset of the natural placement. All four specimen edges
    # register to the displayed meridian, whose present displacement is zero.
    u = np.zeros((t.size, 2))
    u[:, 0] = (-.001 * .25**2 * 16*(t*(1.-t))**2/r**2).ravel()
    u = u.ravel()
    u[free] = spsolve(stiffness[free][:, free], -(stiffness[free][:, ~free] @ u[~free]))
    forces = np.zeros(2)
    for rq, wq, measure, eps in components:
        e = np.einsum("cij,cj->ci", eps, u[dofs])
        isotropic = .5*e[:, :4].sum(axis=1)
        rr, ww, hh, hp = 2.*e[:, :4].T + isotropic
        rz = 2.*e[:, 4]
        time = wq/length
        width, wall_slope = shape.radius(wq)
        widen = 4.*time*(1.-time)
        radial_rate = widen/width
        axial_rate = rq/width*(4.*(1.-2.*time)/length - widen*wall_slope/width)
        # Continuous virtual displacements for width and length are applied
        # directly to the stress; boundary approximation converges on refinement.
        forces -= [measure @ ((rr + hh + hp)*radial_rate + rz*axial_rate), measure @ ww/length]
    return dict(radius=radius, input_snapshot=str(path), input_sha256=digest,
                input_grid=saved_mesh, subsampling=stride, triangles=len(points), natural_offset=.001,
                elastic_opening_force=float(forces[0]), elastic_lengthening_force=float(forces[1]))


if __name__ == "__main__":
    manifest = json.loads(Path("research/pde_ledger_v3/_measurements/mouth_live_waist_order.json").read_bytes())
    rows = [triangle_response(r, 72, stride, manifest) for r in (1, 2) for stride in (2, 1)]
    print(json.dumps(rows, allow_nan=False), flush=True)
