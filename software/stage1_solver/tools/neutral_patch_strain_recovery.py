"""Check the signed impulse after release separately from the driven peak."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from neutral_patch_strain_forecast import LAWS, PULSE_FRACTION, read_material, rms, sha


def main(folder, frozen, expected):
    if sha(frozen) != expected:
        raise ValueError('Frozen hold-based prediction changed')
    with np.load(frozen) as data:
        for law in LAWS:
            path = folder/f'pulse_{law}.npz'
            t, orders, _, _, pulse, _, _ = read_material(path, 'pulse')
            np.testing.assert_array_equal(data['time'], t)
            on = t <= PULSE_FRACTION*t[-1]
            off = t >= PULSE_FRACTION*t[-1]
            for i, order in enumerate(orders):
                direct = pulse[-1, i].mean(axis=0)
                pred = data[law+'_predicted'][-1,i].mean(axis=0)
                load = trapezoid(direct[on], t[on])
                recovery = trapezoid(direct[off], t[off])
                predicted = trapezoid(pred[off], t[off])
                print(json.dumps(dict(test='release_impulse',law=law,order=order,
                                      impulse_unit='patch force times carrier time per unit peak strain',
                                      signed_loading=load,signed_recovery=recovery,predicted_recovery=predicted,
                                      recovery_over_signed_loading=recovery/load,
                                      error_over_recovery_rms=rms(pred[off]-direct[off])/rms(direct[off]),
                                      pulse_sha256=sha(path),forecast_sha256=expected),allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('frozen',type=Path)
    parser.add_argument('expected')
    args = parser.parse_args()
    main(args.folder,args.frozen,args.expected)
