"""Independently check the conditional translating order-front bend on archived inputs."""

import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import linprog
import yaml


ROOT = Path(__file__).resolve().parents[3]
COEFFICIENTS = ROOT / 'research/pde_ledger_v3/_measurements/mouth_bend_coefficients.json'
DERRICK = ROOT / 'software/stage1_solver/reports/pathA_26_derrick_results.yaml'


def archived_row(runs, name, start):
    line, = [row for row in runs[name]['stdout_lines'] if row.startswith(start)]
    return dict(token.split('=', 1) for token in line.split() if '=' in token)


def bend(u):
    return math.sin(math.pi*u)**2


def modal_primitive(u):
    return 3*u/8 - np.sin(2*np.pi*u)/(4*np.pi) + np.sin(4*np.pi*u)/(32*np.pi)


def cellwise_minimum_contrast(target, n):
    # Scale the unknown positive densities by their minimum. Their overall
    # normalization drops out of the modal-to-radial mass ratio. The last
    # linear-program variable is the upper density bound in these units.
    means = n*np.diff(modal_primitive(np.linspace(0., 1., n+1)))
    ceiling = sparse.hstack([sparse.eye(n), -np.ones((n, 1))], format='csr')
    balance = np.r_[means-target, 0.][None, :]
    result = linprog(np.r_[np.zeros(n), 1.], A_ub=ceiling, b_ub=np.zeros(n),
                     A_eq=balance, b_eq=[0.], bounds=(1., None), method='highs')
    assert result.success, result.message
    return result.fun


def main():
    runs = json.loads(COEFFICIENTS.read_text())['runs']
    settings = yaml.safe_load(DERRICK.read_text())['parameters']
    geom = archived_row(runs, 'bend_projection_wave015', 'geometry ')
    fold = archived_row(runs, 'bend_projection_fold64', 'threshold ')
    uniform = archived_row(runs, 'bend_projection_fit64',
                           'peaks cells=64 sheet=uniform_sheet loss=carrier_and_scalar ')
    radius, half = float(geom['radius']), float(geom['half'])
    radial, tension = float(geom['breathing']), settings['sigma']
    surface = 8*math.pi*radius**2*half
    density = radial/surface
    speed = math.sqrt(tension/density)
    angular = quad(lambda z: z*z/2, -1., 1.)[0]
    unweighted = radial*angular*quad(bend, 0., 1.)[0]
    modal = radial*angular*quad(lambda u: bend(u)**2, 0., 1.)[0]
    print(f'homogeneous density={density:.12g} c_order={speed:.12g} '
          f'wall_Px_over_bdot_both={unweighted:.12g} modal_mass_both={modal:.12g} '
          f'modal_mass_half={modal/2:.12g} archived_fine_peaks={uniform["selective_count"]}')
    assert math.isclose(modal, float(uniform['M_wall']), rel_tol=1e-10)
    assert math.isclose(density, settings['rho0']*radius/3, rel_tol=1e-10)

    target = 3*float(fold['wall_mass_at_fold'])/radial
    edge = math.asin(target**.25)/math.pi
    width = 1-2*edge
    positive = quad(lambda u: bend(u)**2-target, edge, 1-edge)[0]
    contrast = 1+(target-3/8)/positive
    low = 1/(1+(contrast-1)*width)
    high = contrast*low
    areas = [quad(lambda u: bend(u)**2, a, b)[0]
             for a, b in [(0., edge), (edge, 1-edge), (1-edge, 1.)]]
    modal_witness = radial*(low*(areas[0]+areas[2])+high*areas[1])/3
    print(f'fine_fold={fold["wall_mass_at_fold"]} target_mean_f_squared={target:.12g} '
          f'high_interval_start={edge:.12g} high_interval_end={1-edge:.12g} '
          f'positive_excess={positive:.12g} density_contrast={contrast:.12g} '
          f'speed_contrast={math.sqrt(contrast):.12g}')
    print(f'fold_witness low_over_mean={low:.12g} high_over_mean={high:.12g} '
          f'c_slow={speed/math.sqrt(high):.12g} c_fast={speed/math.sqrt(low):.12g} '
          f'normalized_total={low*(1-width)+high*width:.12g} modal_mass={modal_witness:.12g}')
    assert math.isclose(modal_witness, float(fold['wall_mass_at_fold']), rel_tol=1e-12)
    for n in (32, 256, 2048):
        discrete = cellwise_minimum_contrast(target, n)
        print(f'cellwise_density_LP n={n} contrast={discrete:.12g} '
              f'excess_above_continuum={discrete-contrast:.5g}')
        assert -1e-10 < discrete-contrast < 1000/n**2


if __name__ == '__main__':
    main()
