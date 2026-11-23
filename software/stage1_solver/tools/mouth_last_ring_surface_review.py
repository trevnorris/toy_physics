"""Sharp nested cylinder and slit fluxes on the independently rematched rejected R4 wall."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_surface_review import Reader
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_supercore_joined_front import Joined, check_guard, emit, file_sha


class FieldChart(Joined):
    """Rebuild the node numbering only. Every read uses the archived physical axes."""
    def matrix(self, order, elastic_density):
        pass


def run(a):
    x, displacement, _, _ = read_source(a.state, a.wall)
    check_guard(a.fields)
    ref = Front(a.state, GRIDS['native'], 1)
    trial = Trial(ref, a.state, x, (displacement.shape[1]-1)//3, displacement.shape[0]-1, False)
    digest = file_sha(a.fields)
    with np.load(a.fields) as f:
        level = f['level'].item()
        if tuple(f['grid']) != GRIDS[level] or f['checkpoint_sha'].item() != file_sha(a.wall):
            raise ValueError('Wrong carrier state or mesh')
        report = [json.loads(row) for row in (a.fields.parent/'stdout').read_text().splitlines()]
        if not any(r.get('kind') == 'fields' and r['sha256'] == digest for r in report):
            raise ValueError('Paired archive differs from guarded result')
        if level == 'native':
            model = trial.model
        else:
            rr = tuple(f[k] for k in ('inside_x','outside_x'))
            z = tuple(f[k] for k in ('feed_z','duct_z','receiver_z'))
            model = FieldChart(trial.shape, ref.c, 1, GRIDS[level], precondition=False,
                               terminal_profile=trial.model.terminal_profile)
            if (any(len(a) != len(b) for a,b in zip((*rr,*z), (*model.rr,*model.z)))
                    or not np.array_equal(model.a, f['angles']) or model.size != f['wound'].shape[1]):
                raise ValueError('Physical axes do not match the field topology')
            model.rr, model.z = rr, z
        reader = Reader(model,ref.c,f['wound'],f['rest'])
    emit(kind='surface_review',level=level,paired_field=digest,rejected_wall=file_sha(a.wall))
    for h in (.16,.235):
        slit=reader.sealed_increment(0.,h,8)
        for radius in (4.55,5.25):
            values={}
            for q in (4,8):
                _,net=reader.integrate(radius,h,q)
                values[q]=net
            emit(kind='surface',radius=radius,physical_half_thickness=h,sealed_slit=slit.tolist(),
                 closed=values[8].tolist(),with_slit=(values[8]+slit).tolist(),
                 quadrature_change=(values[8]-values[4]).tolist())


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','wall','fields'):
        p.add_argument(name,type=Path)
    run(p.parse_args())
