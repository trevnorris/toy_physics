"""Shape work on a flowing, initially wave-free order field just after wave preparation."""
import argparse
import json
import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_live_waist_order import Order, load_for, state
from mouth_sealed_shape_relax import Meridian


def measure(radius, tag, mesh, field_mesh):
    sample = state(radius, tag)
    coords = np.array(sample["coordinates"])
    geometry = Meridian(radius, coords)
    order = Order(geometry, mesh)
    fields = load_for(geometry, mesh, field_mesh, sample["wave_frequency"])
    chi, flow = order.condition(order.seed, fields, 0., BIAS)
    ef, gf, info = fields.evaluate(chi, 1., BIAS)
    eo, go = order.energy_gradient(chi)
    print(json.dumps(dict(kind="center", radius=radius, tag=tag, mesh=mesh,
                          field_mesh=field_mesh, chi_axis=order.axis(chi), minimum=float(chi.min()),
                          maximum=float(chi.max()), order_energy=eo, coupled_energy=eo+ef,
                          flow_stationarity=flow["order_residual"],
                          new_order_residual=float(np.linalg.norm((gf+go)[order.free])),
                          **info), allow_nan=False), flush=True)
    for coordinate, index, scale in (("radial", 2, 1/radius), ("axial", 0, 1.)):
        for offset in (-.002, -.001, .001, .002):
            x = coords.copy()
            x[index] += scale*offset
            shifted = Meridian(radius, x)
            grid = Order(shifted, mesh)
            load = load_for(shifted, mesh, field_mesh, info["wave_frequency"], fields.wave.seed)
            ef, _, record = load.evaluate(chi, 1., BIAS)
            eo, _ = grid.energy_gradient(chi)
            print(json.dumps(dict(kind="signed", axis=coordinate, offset=offset, order_energy=eo,
                                  coupled_energy=eo+ef, **record), allow_nan=False), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--radius", type=int, choices=(1, 2), required=True)
    p.add_argument("--point", choices=("broad", "mid", "tight"), required=True)
    p.add_argument("--mesh", type=int, required=True)
    p.add_argument("--field-mesh", type=int, required=True)
    a = p.parse_args()
    measure(a.radius, a.point, a.mesh, a.field_mesh)
