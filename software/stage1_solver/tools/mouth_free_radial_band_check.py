"""Check radial equilibrium of each sampled shape independently in eight depth bands."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def pressure_bound(design, pressure):
    """Minimum achievable worst pressure error at nonnegative tension and wave action."""
    scales = np.linalg.norm(design, axis=0)
    scaled = design/scales
    count = len(pressure)
    constraints = np.column_stack((np.vstack([scaled,-scaled]), -np.ones(2*count)))
    fit = linprog([0.,0.,0.,0.,1.], A_ub=constraints, b_ub=np.r_[pressure,-pressure],
                  bounds=[(None,None), (None,None), (0,None), (0,None), (0,None)], method='highs',
                  options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
    if not fit.success:
        raise ArithmeticError(fit.message)
    params = fit.x[:-1]/scales
    hi, lo = np.split(fit.ineqlin.marginals,2)
    witness = hi-lo
    residual = pressure-design@params
    return dict(optimum=fit.fun, parameters=params.tolist(), residuals=residual.tolist(),
                signed_band_weights=witness.tolist(), projected_columns=(witness@design).tolist(),
                norm_one=np.abs(witness).sum(), lower_bound=witness@pressure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', type=Path, nargs='+')
    args = parser.parse_args()
    for folder in args.runs:
        outcome = json.loads((folder/'outcome.json').read_text())
        if not outcome['limitsVerified'] or outcome['exitCode'] or outcome['childOutcome']['exitCode']:
            raise ValueError(f'Incomplete guarded shape run at {folder}')
        source = (folder/'stdout').read_bytes()
        rows = [json.loads(line) for line in source.splitlines()]
        forces = [row for row in rows if row['test']=='force']
        if len(forces)!=5 or len({row['q'] for row in forces})!=5:
            raise ValueError(f'Incomplete signed and half-strength shape rows at {folder}')
        for row in sorted(forces,key=lambda row: row['q']):
            bands = row['bands']
            if len(bands) != 8:
                raise ValueError('Expected eight bands')
            design = np.array([[1., b['mean_depth'], b['curvature'], -b['wave_pressure_per_action']]
                               for b in bands])
            surface, volume = [[b[key] for b in bands] for key in ('carrier_surface','carrier_volume')]
            max_readout_difference = max(abs(np.subtract(surface, volume)))
            for label, pressure in [('volume', volume), ('surface', surface)]:
                result = pressure_bound(design, np.array(pressure))
                print(json.dumps(dict(source_sha256=sha256(source).hexdigest(), mesh=row['mesh'],
                                      q=row['q'], carrier=label,
                                      max_band_surface_volume_difference=max_readout_difference,
                                      **result),allow_nan=False),flush=True)
                # Each signed witness bounds every permitted wall law at this profile.
                np.testing.assert_allclose(result['projected_columns'][:2],0.,atol=3e-10)
                if max(result['projected_columns'][2:])>3e-10 or result['norm_one']>1+1e-10:
                    raise ArithmeticError('Dual witness violates a material sign or normalization')
                np.testing.assert_allclose(result['lower_bound'],result['optimum'],atol=2e-10,rtol=0)
                np.testing.assert_allclose(max(abs(np.array(result['residuals']))),result['optimum'],
                                           atol=2e-10,rtol=0)


if __name__=='__main__':
    main()
