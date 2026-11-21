"""Independent quadrature and three-amplitude spectrum for two trial wall profiles."""

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.signal import find_peaks
import yaml


ROOT = Path(__file__).resolve().parents[3]
MEASUREMENTS = ROOT / 'research/pde_ledger_v3/_measurements'


def archived_runs(name):
    return json.loads((MEASUREMENTS / f'{name}.json').read_text())['runs']


def rows(lines, prefix):
    return [dict(field.split('=', 1) for field in line.split() if '=' in field)
            for line in lines if line.startswith(prefix)]


def interpolate(table, x, column):
    offsets = [float(s['offset']) for s in table]
    assert np.all(np.diff(offsets) > 0) and offsets[0] <= x.min() <= x.max() <= offsets[-1]
    return np.interp(x, offsets, [float(s[column]) for s in table])


def absorption_width(shift, matrix):
    drive = np.zeros((len(shift), 3, 1))
    drive[:, 0, 0] = 1
    response = np.linalg.solve(matrix, drive)[:, 0, 0].imag
    peak, = find_peaks(response)[0]
    half = response[peak] / 2
    def crossing(index):
        return brentq(lambda xx: np.interp(xx, shift, response) - half,
                      shift[index], shift[index + 1], xtol=1e-15)
    lower = crossing(np.flatnonzero(response[:peak] < half)[-1])
    upper = crossing(peak + np.flatnonzero(response[peak:] < half)[0] - 1)
    return upper - lower


def selective_readout(cells, mass, geometry, direct, coeff, trap, carrier, action):
    preparation = trap['smooth_dipole_l015_01']['stdout_lines']
    point, = rows(preparation, 'point n=1 ')
    energy, = rows(preparation, 'levels ')
    bend, = rows(trap['mouth_wave_linewidth_budget01']['stdout_lines'], 'bow h=0.015 ')
    ws, wd = float(point['omega_s']), float(point['omega_d'])
    gap = abs(float(energy['U_over_action']))
    coupling = float(bend['overlap']) * math.sqrt(2 * action / ws)
    offset = np.linspace(-.016, .08, 48001)
    z = wd + offset
    matrix = np.zeros((len(z), 3, 3), complex)
    matrix[:, 0, 0] = wd**2 - z**2
    matrix[:, 0, 1:] = matrix[:, 1:, 0] = coupling / 2
    run = {48: 'open_reservoir_lines48', 64: 'final_open_lines64'}[cells]
    lines = carrier[run]['stdout_lines'] + coeff[f'bend_projection_rad{cells}']['stdout_lines']
    waves = rows(coeff['bend_projection_wresp03']['stdout_lines'], 'wave ')
    common = interpolate(waves, offset, 'central')
    matrix[:, 1, 2] = matrix[:, 2, 1] = direct / 2 - common
    for j, tone, nu in ((1, 'difference', z-ws), (2, 'sum', z+ws)):
        table = sorted(rows(lines, f'tone={tone} '), key=lambda s: float(s['offset']))
        for s in table:
            s['imag'] = 2 * float(s['power_per_b2']) / float(s['omega'])
        fluid = interpolate(table, offset, 'induced_real') + 1j*interpolate(table, offset, 'imag')
        wave = interpolate(waves, offset, 'lower') if j == 1 else (
            interpolate(waves, offset, 'upper_real') + 1j*interpolate(waves, offset, 'upper_imag'))
        matrix[:, j, j] = (float(geometry['tensile_both']) + direct - mass*nu**2
                           - common - wave - fluid)
    width = absorption_width(offset, matrix)
    # A, 3A are the half widths, with equal areas. Count real roots of
    # minus the derivative with positive slope, using the line spacing as unit.
    a = width / (2*gap)
    x = np.polynomial.Polynomial([0., 1.])
    numerator = x*((x+1)**2+9*a**2)**2 + 3*(x+1)*(x**2+a**2)**2
    peaks = sum(abs(r.imag) < 1e-9 and numerator.deriv()(r.real) > 0 for r in numerator.roots())
    return width, width/gap, peaks


def main():
    coeff = archived_runs('mouth_bend_coefficients')
    trap = archived_runs('mouth_trapped_wave_spectrum')
    carrier = archived_runs('mouth_open_reservoirs')
    archive = yaml.safe_load((ROOT/'software/stage1_solver/reports/pathA_26_derrick_results.yaml').read_text())
    geo, = rows(coeff['bend_projection_wave015']['stdout_lines'], 'geometry ')
    support, = rows(coeff['bend_projection_wave015']['stdout_lines'], 'support ')
    r, h, radial = (float(geo[k]) for k in ('radius', 'half', 'breathing'))
    parameters = archive['parameters']
    uniform_speed = math.sqrt(parameters['sigma'] * (8*math.pi*r**2*h) / radial)
    avg = 1+8*35/128
    all_times = []
    for s in (0., .15):
        def density(u):
            return (1+8*math.sin(math.pi*(u-s))**8)/avg
        def local_speed(u):
            return uniform_speed / math.sqrt(density(u))
        quarters = [h*quad(lambda u: 1/local_speed(u), j/4, (j+1)/4,
                           epsabs=1e-13, epsrel=1e-13)[0] for j in range(4)]
        all_times.append(sum(quarters))
        weighted = quad(lambda u: math.sin(math.pi*u)**4*density(u), 0., 1.)[0]
        # Only the constant, first and second cosine modes of the shifted
        # density contribute to the bend projection.
        cross = 105/1024 + 7/64*math.cos(2*math.pi*s) + 7/512*math.cos(4*math.pi*s)
        mass = radial/3 * (3/8+8*cross)/avg
        numeric = radial/3 * weighted
        print(f'shift={s:g} M_exact={mass:.12g} M_integral={numeric:.12g} '
              f'c_min={local_speed((s+.5)%1):.12g} c_max={local_speed(s):.12g} '
              f'total_T={sum(quarters):.12g} quarter_times='+','.join(f'{v:.12g}' for v in quarters))
        assert math.isclose(numeric, mass, rel_tol=1e-13)
        assert math.isclose(quad(density, 0., 1.)[0], 1., rel_tol=1e-13)
        for cells in (48, 64):
            width, ratio, count = selective_readout(cells, mass, geo, float(support['direct_K']),
                                                   coeff, trap, carrier, parameters['I_wave'])
            print(f'shift={s:g} cells={cells} linear_table_FWHM={width:.12g} '
                  f'width_over_gap={ratio:.12g} raw_peaks=1 selective_peaks={count}')
            fold, = rows(coeff[f'bend_projection_fold{cells}']['stdout_lines'], 'threshold ')
            assert count == 1 + (mass > float(fold['wall_mass_at_fold']))
    assert math.isclose(*all_times, rel_tol=1e-13)


if __name__ == '__main__':
    main()
