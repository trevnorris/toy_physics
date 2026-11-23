"""Hold the historical rejected fine geometry with the new interior coordinate."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_fixed_wall_angular import solve_fixed
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import check_guard, emit, file_sha


def check(state, snapshot, target, folder):
    check_guard(target)
    prep, _, x, d, u = accepted(state,snapshot)
    data = json.loads(target.read_text())
    if (data['accepted_sha256'] != file_sha(snapshot)
            or not np.array_equal(x,data['original_x']) or prep['material_cells'] != [33,32]):
        raise ValueError('The historical candidate and material state differ')
    emit(kind='preparation', state=file_sha(state),snapshot=file_sha(snapshot),
         coordinates=data['x'], target=file_sha(target), old_residual=data['prior_attempt_residual'])
    solve_fixed.input_sha = file_sha(snapshot)
    solve_fixed(state,np.asarray(data['x']),d,u,folder)
    emit(kind='saved',sha256=file_sha(folder/'equilibrium.npz'))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('state','snapshot','target','folder'):
        p.add_argument(arg, type=Path)
    check(**vars(p.parse_args()))
