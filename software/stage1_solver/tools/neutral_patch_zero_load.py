"""Predict compression from a stress-free in-plane pitch at fixed layer gap."""
import numpy as np
from scipy.optimize import brentq

from neutral_patch_carrier_compression import CALIBRATION, DENSITIES, carrier_pressure
from neutral_patch_periodic_compression import Bilayer
from neutral_patch_sites import GAP, emit, morse


def stress_free_pitch(cell):
    # Only the planar lattice changes during preparation. The old fixed-gap
    # virial vanishes at dE/d(log pitch)=0, independent of volume normalization.
    def pressure_at_pitch(pitch):
        return cell.equilibrium(pitch**-3, False)[1]

    bracket = (0.95, 1.02)
    pitch = brentq(pressure_at_pitch, *bracket, xtol=3e-15)
    return pitch, [pressure_at_pitch(a) for a in bracket]


def homogeneous_checks(cell, relaxed):
    planar = cell.vectors + relaxed[cell.j] - relaxed[cell.i]
    vertical_squared = cell.cross * GAP**2
    radius = np.sqrt(np.sum(planar**2, axis=1) + vertical_squared)
    _, first = morse(radius, *cell.params)
    # Derivatives of cell energy per two sites with respect to planar strain.
    stress = planar.T @ ((first/radius)[:, None] * planar) / cell.count
    np.testing.assert_allclose(stress, 0., rtol=0, atol=2e-10)

    def energy(diagonal_strain):
        deformed = planar * np.exp(diagonal_strain)
        radii = np.sqrt(np.sum(deformed**2, axis=1) + vertical_squared)
        return morse(radii, *cell.params)[0].sum()/cell.count

    # Sixfold symmetry leaves one dilation and two equivalent planar shears.
    h = 1.e-4
    curvatures = {}
    for name, direction in [('dilation', [1, 1]), ('area_preserving_shear', [1, -1])]:
        step = h*np.array(direction)
        curve = (energy(step) + energy(-step) - 2*energy(0.))/h**2
        if curve <= 0:
            raise ArithmeticError(f'Unstable homogeneous {name} mode')
        curvatures[name] = float(curve)
    return stress.tolist(), curvatures


def measure(nx, ny, shell):
    cell = Bilayer(nx, ny, True, shell)
    pitch, bracket_pressures = stress_free_pitch(cell)
    cell.reference *= pitch
    cell.vectors *= pitch
    energy, rest, offset = cell.equilibrium(1., False)
    stress, modes = homogeneous_checks(cell, offset)
    rng = np.random.default_rng(124)
    shaken = cell.equilibrium(1., False, rng.normal(0., .02, offset.shape))
    np.testing.assert_allclose(shaken[2], offset, rtol=0, atol=2e-7)
    np.testing.assert_allclose(shaken[:2], [energy, rest], rtol=0, atol=2e-9)
    rows = []
    for density in DENSITIES:
        _, pressure, state = cell.equilibrium(density, False, offset)
        least_mode = cell.stability(density, False, state)
        volume_ratio = 1 / density
        dv = 1e-5
        def relaxed_energy(volume):
            return cell.equilibrium(1/volume, False, state)[0]
        from_energy = (relaxed_energy(volume_ratio-dv)-relaxed_energy(volume_ratio+dv))/(2*dv)
        np.testing.assert_allclose(from_energy, pressure, rtol=1e-7, atol=5e-8)
        rows.append(dict(n=density, change=pressure-rest, min_stiffness=least_mode,
                         work_check_error=abs(from_energy-pressure)))
    curvature = []
    for h in (0.001, 0.0005):
        plus, minus = [cell.equilibrium(np.exp(s*h), False, offset)[1]-rest for s in (1, -1)]
        curvature.append((plus+minus)/(h*(plus-minus)/2))
    return dict(sites=2*cell.count, shell=shell, pitch=pitch, bracket_pressures=bracket_pressures,
                rest_energy=energy, rest_pressure=rest, in_plane_stress=stress, homogeneous_modes=modes,
                curvature_log_steps=[0.001, 0.0005], curvature=curvature, rows=rows)


def main():
    cases = [(4, 3, 6), (4, 3, 9), (10, 10, 6)]
    baseline = measure(*cases[0])
    fit = next(row for row in baseline['rows'] if row['n'] == CALIBRATION)
    if fit['change'] <= 0:
        raise ArithmeticError('The 1% density calibration requires a positive pressure increase')
    scale = carrier_pressure(CALIBRATION)/fit['change']
    for index, case in enumerate(cases):
        result = baseline if index == 0 else measure(*case)
        np.testing.assert_allclose(result['pitch'], baseline['pitch'], rtol=0, atol=1e-10)
        for row, reference in zip(result['rows'], baseline['rows'], strict=True):
            np.testing.assert_allclose(row['change'], reference['change'], rtol=0, atol=3e-9)
            target = carrier_pressure(row['n'])
            prediction = scale*row['change']
            row.update(target=target, predicted=prediction,
                       fractional_error=prediction/target-1 if row['n'] != 1. else None)
        emit(test='zero_in_plane_load', fit_density=CALIBRATION, shared_scale=scale, **result)


if __name__ == '__main__':
    main()
