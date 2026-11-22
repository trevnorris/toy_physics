"""Crystal virtual work across both sides of a conditionally welded held waist."""
import json

import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import splu

from mouth_unified_puncture import Annulus, INSET, OFFSET, RECORD, ROOT, archived
from mouth_waist_outer_order import exterior

BALANCE = ROOT / "research/pde_ledger_v3/_measurements/mouth_coupled_waist_balance.json"
CUT = 4.6


class ReceiverBand(Annulus):
    """An equal number of cells outside the half-order surface, from x=1 to 1.5."""

    def __init__(self, shape, chart, wa, xout, field):
        super().__init__(shape, chart)
        self.x += 1-INSET
        width = shape.radius(self.t*shape.length)[0]
        self.r = self.x*width
        outer_coordinate = (self.r-width)/(CUT-width)
        order = RegularGridInterpolator((wa, xout), field)
        self.chi = order(np.stack((self.t*shape.length, outer_coordinate), -1))


def dof_map(ids):
    return np.stack((2*ids, 2*ids+1), -1).ravel()


def joined(radius, mesh, receiver_mesh, manifest, balance):
    shape, order, source = archived(radius, mesh, manifest)
    info, (wa, xout, field) = exterior(shape, receiver_mesh, CUT, return_side_field=True)
    archived_row, = [r for r in balance["readings"] if r["radius"] == radius and
                     r["state"] == "broad" and r["order_mesh"] == mesh and
                     r["exterior_mesh"] == receiver_mesh]
    center, = [r for r in archived_row["energies"]["radial"] if r["displacement"] == 0]
    if abs(info["energy"]-center["exterior_order"]) > 1.e-9:
        raise ValueError("Exterior field does not reproduce the archived central energy")
    inner = Annulus(shape, order)
    outer = ReceiverBand(shape, order, wa, xout, field)
    jump = np.max(abs(inner.chi[:, -1]-outer.chi[:, 0]))
    if jump > 1.e-12 or np.min(outer.chi) < 0:
        raise ValueError("Order discontinuity or negative material modulus on the specimen")
    # Share both displacement components at b. Only the two remote cuts
    # and the ends have prescribed displacement. No traction is imposed at b.
    ids = np.arange((2*mesh+1)*(mesh+1)).reshape(2*mesh+1, mesh+1)
    half = mesh//2
    numbers = [dof_map(part) for part in (ids[:, :half+1], ids[:, half:])]
    bands = (inner, outer)
    matrices, volumes = zip(*(band.stiffness() for band in bands))
    entries, rows, cols = [], [], []
    n = 2*ids.size
    modes = np.zeros((n, 3))
    for band, mat, dofs in zip(bands, matrices, numbers):
        local = mat.tocoo()
        entries.append(local.data)
        rows.append(dofs[local.row])
        cols.append(dofs[local.col])
        modes[dofs] = band.boundary_modes()
    joint = sparse.coo_matrix((np.concatenate(entries), (np.concatenate(rows), np.concatenate(cols))),
                               shape=(n, n)).tocsc()
    free = np.ones(n, dtype=bool)
    cuts = {name: dof_map(nodes) for name, nodes in (
        ("inner_radial_cut", ids[1:-1, 0]), ("outer_radial_cut", ids[1:-1, -1]),
        ("inlet_end", ids[0]), ("outlet_end", ids[-1]))}
    for nodes in cuts.values():
        free[nodes] = False
    relaxed = modes.copy()
    relaxed[free] = splu(joint[free][:, free]).solve(-(joint[free][:, ~free] @ modes[~free]))
    local_reactions = [mat @ relaxed[dofs] for mat, dofs in zip(matrices, numbers)]
    residual = joint @ relaxed
    inside_join = dof_map(np.arange(inner.r.size).reshape(inner.r.shape)[1:-1,-1])
    outside_join = dof_map(np.arange(outer.r.size).reshape(outer.r.shape)[1:-1,0])
    interface = [r[ix] for r, ix in zip(local_reactions, (inside_join, outside_join))]
    jump_work = np.max(abs(interface[0]+interface[1]))
    join_dofs = dof_map(ids[:, half])
    join_motion = relaxed[join_dofs].reshape(-1, 2, 3)
    join_test = modes[join_dofs].reshape(-1, 2, 3)
    # At fixed baseline coordinates, the normal is (1,-b')/sqrt(1+b'^2).
    # The natural-map column below gives relative, not absolute, displacement.
    # Absolute displacement per +d at the held cuts is test-minus-relaxed.
    widths, slopes = shape.radius(inner.t[:, -1]*shape.length)
    def normal(motion):
        return (motion[:, 0]-slopes[:, None]*motion[:, 1])/np.hypot(1., slopes)[:, None]
    join_normal = normal(join_motion)
    test_normal = normal(join_test)
    waist = np.argmin(widths)
    data = []
    for d in (-OFFSET, 0., OFFSET):
        # Relative displacement in natural coordinates, to first order in d.
        relative = -d*relaxed[:, 0]
        force = d*residual[~free, 0] @ modes[~free, 1:]
        boundary = {name: (d*residual[nodes, 0] @ modes[nodes, 1:]).tolist()
                    for name, nodes in cuts.items()}
        traces = modes[dof_map(ids[1:-1, half]), 1:]
        interface_work = [(d*r[:, 0] @ traces).tolist() for r in interface]
        # The stress of the outer band is integrated only on b<r<1.5b.
        # Its interfacial nodal traction on the common join is opposite the
        # inner reaction. Neither adds a surface or scalar-front energy.
        stored, volume_work = [], []
        for band, mat, dofs in zip(bands, matrices, numbers):
            local = relative[dofs]
            stored.append(float(.5*local @ (mat @ local)))
            volume_work.append((-band.stress_work(local, modes[dofs, 1:])).tolist())
        signed = []
        for m in (1, 2):
            def elastic_energy(q):
                u = relative + q*relaxed[:, m]
                return .5*u @ (joint @ u)
            h = .0002
            signed.append((elastic_energy(-h)-elastic_energy(h))/(2*h))
        if max(abs(force-signed)) > 2.e-10 or max(abs(force-np.sum(volume_work, axis=0))) > 2.e-10:
            raise ArithmeticError("Joined stress and signed shape energies disagree")
        data.append(dict(natural_offset=d, elastic_force=force.tolist(), from_signed_energy=signed,
                         force_at_external_cuts=boundary, stored_elastic_by_side=stored,
                         internal_join_reaction_work_by_side=interface_work,
                         stress_volume_work_by_side=volume_work))
    return dict(radius=radius, order_mesh=mesh, receiver_mesh=receiver_mesh, **source,
                outer_order=info, archived_outer_energy=center["exterior_order"],
                volume_by_side=list(volumes), specimen_outer_min_order=float(np.min(outer.chi)),
                interface_order_jump=float(jump), free_residual=float(np.max(abs(residual[free]))),
                max_join_reaction_sum_for_unit_modes=float(jump_work), sides=["inner","outer"],
                waist_normal_displacement_per_unit_cut_motion=join_normal[waist, 1:].tolist(),
                waist_normal_displacement_of_affine_chart=test_normal[waist, 1:].tolist(),
                max_join_normal_shift_under_positive_offset=float(
                    OFFSET*np.max(abs((test_normal-join_normal)[:, 0]))),
                force_axes=["opening","physical lengthening"], data=data)


def main():
    manifest = json.loads(RECORD.read_bytes())
    balance = json.loads(BALANCE.read_bytes())
    data = [joined(r, i, o, manifest, balance) for r in (1, 2) for i, o in ((48, 40), (72, 70))]
    print(json.dumps(data, allow_nan=False, indent=2))
    if any(max(row["free_residual"], row["max_join_reaction_sum_for_unit_modes"]) > 2.e-11
           for row in data):
        raise ArithmeticError("Joined specimen has an unresolved material reaction")


if __name__ == "__main__":
    main()
