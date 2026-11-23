"""Measure held-exit virtual work and flux using the same wound phase in both reservoirs."""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_exchange_carrier import BIAS
from mouth_poly_shallow_release import emit
from mouth_sealed_shape_relax import Meridian
from mouth_wound_carrier import Wound, prepared


def exit_shape(radius, original, exit_radius):
    config = original.copy()
    config[1] = np.log(exit_radius / (1.5 * radius))
    shape = Meridian(radius, config)
    if not np.isclose(shape.exit, exit_radius, rtol=1.e-13, atol=0.):
        raise ValueError('Exit map did not preserve the physical radius')
    return shape


def cap_phase_traces(model, phase):
    """Use the disk quadrature and the two independent meridional phase traces."""
    axis = (leggauss(3)[0] + 1) / 2
    minus, plus = phase[model.minus], phase[model.plus]
    pairs = []
    for a in axis:
        for c in axis:
            basis = np.kron([1-a, a], [1-c, c])
            pairs.append(np.column_stack((minus @ basis, plus @ basis)))
    return np.concatenate(pairs)


def summarize_disk(fields, model, info):
    area, _, _, s, j, density, ur, ut, spin = fields['disk'].T
    normal = np.abs(j) / density**3
    total = np.sqrt((j/density)**2 + ur**2 + ut**2 + spin**2) / density**2
    positive = area * np.maximum(j, 0)
    angular = (area*j) @ (s*spin)
    if abs(angular - model.count*info['total_current']) > 1.e-7 * (1+abs(angular)):
        raise ArithmeticError('Terminal cut did not conserve angular momentum')
    if not np.isfinite(total).all() or density.min() <= 0 or positive.sum() <= 0:
        raise ArithmeticError('Invalid disk measurements')

    def fractions(mask):
        return dict(area=area[mask].sum()/area.sum(), inward_flux=positive[mask].sum()/positive.sum())

    return dict(current=area@j, negative_flux_fraction=np.sum(np.maximum(-area*j, 0))/positive.sum(),
                min_normal_mach=normal.min(), min_total_mach=total.min(),
                rim_density=density[s.argmax()], max_density=density.max(),
                rim_normal_mach=normal[s.argmax()], rim_total_mach=total[s.argmax()],
                normal_subsonic=fractions(normal < 1), total_subsonic=fractions(total < 1),
                density_half_bulk=fractions(density >= .5),
                dense_and_normal_subsonic=fractions((normal < 1) & (density >= .5)),
                dense_and_total_subsonic=fractions((total < 1) & (density >= .5)),
                outer_half_swirl_radius=fractions(s >= .5 * model.shape.exit),
                angular_momentum_out=angular)


def run(args):
    original, order, digest = prepared(args.radius)
    b = args.exit
    reference = Meridian(args.radius, original).exit
    emit(kind='input', radius=args.radius, exit=b, unwound_reference_exit=reference,
         material_sha256=digest, mesh=args.mesh, polar=args.polar, half_step=args.h,
         meridional_bias=BIAS, outer=args.outer, roof=args.roof, winding=args.winding,
         circulation=2*np.pi*args.winding, steps=args.steps)
    results = {}
    for relative_step in args.steps:
        shape = exit_shape(args.radius, original, b*(1+relative_step))
        model = Wound(shape, order, args.winding, args.mesh, args.polar, args.outer, args.roof)
        amplitude, phase, info = model.solve(BIAS)
        results[relative_step] = (info['excess_energy'], info['terms'])
        if relative_step == 0:
            fields = model.center_fields(amplitude, phase, info)
            phase_traces = cap_phase_traces(model, phase)
            if len(phase_traces) != len(fields['disk']):
                raise ArithmeticError('Disk trace and quadrature lengths differ')
            emit(kind='disk', radius=args.radius, exit=b, winding=args.winding, **summarize_disk(fields,model,info))
            args.save.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(args.save, **fields, cap_phase_traces=phase_traces,
                                material_sha256=digest,winding=args.winding,radius=args.radius,exit=b,
                                shape_params=shape.params, outlet_phase_bias=BIAS,
                                disk_columns=['area','r/b','polar_angle','s','j_w','density','v_r','v_polar','v_phi'])
        emit(kind='state', radius=args.radius, exit=shape.exit, relative_step=relative_step,
             excess_energy=info['excess_energy'], current=info['total_current'], terms=info['terms'],
             newton_error=info['residual'], iterations=info['iterations'])
    if -args.h in results and args.h in results:
        lo, lower = results[-args.h]
        hi, upper = results[args.h]
        scale = 2 * args.h * b
        force = (lo-hi) / scale
        components = {k:(lower[k]-upper[k])/scale for k in lower}
        if abs(sum(components.values())-force) > 2e-8*(1+abs(force)):
            raise ArithmeticError('Terms fail to recover total held work')
        emit(kind='force',radius=args.radius,exit=b,winding=args.winding,carrier_outward=force,parts=components)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius',type=int,choices=(1,2))
    parser.add_argument('exit',type=float)
    parser.add_argument('save',type=Path)
    parser.add_argument('--mesh',type=int,default=13)
    parser.add_argument('--polar',type=int,default=10)
    parser.add_argument('--winding',type=int,default=1,choices=(0,1))
    parser.add_argument('--h',type=float,default=.006)
    parser.add_argument('--outer',type=float,default=28.)
    parser.add_argument('--roof',type=float,default=22.)
    parser.add_argument('--center-only',action='store_true')
    args=parser.parse_args()
    if args.exit <= 0 or not 0 < args.h < .1:
        parser.error('Expected a positive exit and a small positive force step')
    args.steps = [0.] if args.center_only else [0.,-args.h,args.h]
    run(args)
