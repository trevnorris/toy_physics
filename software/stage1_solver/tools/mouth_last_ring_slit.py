"""Locate the carrier-grid change in axial traction along a fixed held rim.

Use the saved wound driven/resting pairs, interpolate each on its own carrier
nodes, and integrate both faces over the same physical strips on all grids.
No new trajectory, wave, or material response is solved here.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_endpoint import wall_reaction
from mouth_closed_rim_surface_review import Reader
from mouth_joint_polar_shape import Trial
from mouth_joined_front_audit import Front
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_supercore_joined_front import Joined, check_guard, emit, file_sha


class Topology(Joined):
    def matrix(self, order, elastic_density):
        """Reader rebuilds the physical gradients directly from saved field nodes."""


def integrate(state, checkpoint, fields):
    check_guard(fields)
    rows = [json.loads(line) for line in (fields.parent/'stdout').read_text().splitlines()]
    if not any(row.get('kind') == 'fields' and row['sha256'] == file_sha(fields) for row in rows):
        raise ValueError('The carrier field pair has no verified readout')
    x, _, _, _ = read_source(state, checkpoint)
    base = Front(state, GRIDS['native'], 1)
    trial = Trial(base, state, x, 11, 32, False)
    with np.load(fields) as arrays:
        if arrays['checkpoint_sha'].item() != file_sha(checkpoint):
            raise ValueError('Field pair belongs to a different solid')
        mesh = Topology(trial.shape, base.c, 1, arrays['grid'],
                        terminal_profile=trial.model.terminal_profile, precondition=False)
        for field_name, shape in (('wound', (3,mesh.size)), ('rest', (3,mesh.size))):
            if arrays[field_name].shape != shape:
                raise ValueError('Incompatible saved field coordinates')
        mesh.rr = (arrays['inside_x'], arrays['outside_x'])
        mesh.z = (arrays['feed_z'], arrays['duct_z'], arrays['receiver_z'])
        mesh.a = arrays['angles']
        reader = Reader(mesh, base.c, arrays['wound'], arrays['rest'])
        level = str(arrays['level'])
    length = mesh.shape.length
    emit(kind='source', level=level, checkpoint=file_sha(checkpoint), fields=file_sha(fields),
         halfwidth=mesh.h, exit_radius=mesh.shape.exit)
    # Four halving windows *inside* the same weak integration plateau.
    # Above this z every wall point lies outside the polar shape deformation.
    previous = np.zeros((2,2))
    previous_depth = 0.
    for fraction in (1/16, 1/8, 1/4, 1/2, 1., 1.25):
        depth = mesh.h*fraction
        def flat(z):
            return np.ones_like(z)
        coarse = wall_reaction(reader, flat, length-depth, 6)
        total = wall_reaction(reader, flat, length-depth, 10)
        if np.max(abs(coarse-total)) > 2e-6:
            raise ArithmeticError('Need membrane-trace quadrature refinement')
        strip = total-previous
        emit(kind='strip', level=level, depth_from_exit=[previous_depth,depth],
             cumulative_bank_mass_axial=total.tolist(), strip_bank_mass_axial=strip.tolist(),
             strip_axial_sum=float(strip[:,1].sum()), cumulative_axial_sum=float(total[:,1].sum()),
             mean_axial_force_per_depth=float(sum(strip[:,1])/(depth-previous_depth)),
             bank_quadrature_change=float(np.max(abs(coarse-total))))
        previous_depth, previous = depth, total


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('state','checkpoint','fields'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    integrate(args.state,args.checkpoint,args.fields)
