"""Independent energy-current check on four continuum translations of saved feed fields."""
import math
from pathlib import Path

import numpy as np

from mouth_carrier_momentum import read_state
from mouth_release_input_check import manifest_fields


def aperture_data(n, base, driven):
    records = Path("docs/_measurements")
    scratch = Path("_scratch/s11c")
    packets = {}
    for speed, folder in ((0., base), (.01, driven), (.02, driven)):
        manifest = "mouth_feed_weak_run.json" if speed == 0 else "mouth_intake_redistribution.json"
        f, current, dx, dz = read_state(scratch/folder/"fields", n, "restoring_slope",
                                         speed, manifest_fields(records/manifest))
        x = (np.arange(n) + .5 - n/2)*dx
        x2, y2, z2 = np.meshgrid(x*x, x*x, x*x, indexing="ij", sparse=True)
        radius2 = x2 + y2 + z2
        weight = np.maximum(1-radius2, 0.)**3
        face = f[-1]
        wall = (-weight*(1-1/abs(face)) + 1j*current/abs(face)**2)*face
        layers = (11*f[-1]-18*f[-2]+9*f[-3]-2*f[-4])/(6*dz)
        if np.any(wall[radius2>=1]) or np.max(abs(f[0]-1.)) > 1e-12:
            raise ValueError("The support or the feed reference differs from the calculation")
        packets[round(100*speed)] = face, wall, layers
    return packets, radius2 < 1, dx


def readings(packet, mask, dx, stream, translation):
    relative = stream-translation
    face, wall, layers = packet
    if relative < 0:
        face, wall, layers = (arr[::-1] for arr in (face, wall, layers))
    frequency = .25 + stream**2/2
    k = 2*math.pi*np.fft.fftfreq(face.shape[0], d=dx)[:, None, None]
    # An even periodic grid has an ambiguous Nyquist derivative. Use the
    # symmetric trigonometric interpolant, shared by both signed translations.
    spectrum = np.fft.fft(face, axis=0)
    slope_multiplier = 1j*k.copy()
    slope_multiplier[face.shape[0]//2] = 0
    dx_face = np.fft.ifft(slope_multiplier*spectrum, axis=0)

    def field_at(time):
        shift = np.exp(-1j*k*translation*time)
        mid = face.shape[0]//2
        shift[mid] = np.cos(k[mid]*translation*time)
        return np.exp(-1j*frequency*time)*np.fft.ifft(spectrum*shift, axis=0)

    # Direct lab-time differentiation on the moving envelope at fixed spatial
    # coordinate. The e^(i U x) factor cancels against the normal derivative.
    dt = .025
    time_derivative = (field_at(-2*dt)-8*field_at(-dt)
                       +8*field_at(dt)-field_at(2*dt))/(12*dt)
    def integral(values):
        return dx**3*values[mask].sum()

    outputs = []
    for normal in (wall, layers):
        number = integral((face.conj()*normal).imag)
        momentum = integral((dx_face.conj()*normal).real)
        direct = -integral((time_derivative.conj()*normal).real)
        product = frequency*number + translation*momentum
        if abs(direct-product)>2e-11:
            raise ValueError("Boosted current and time derivative disagree")
        outputs.append((number, momentum, product, direct))
    return np.array(outputs)


def compare(n, base, driven):
    packets, mask, dx = aperture_data(n, base, driven)
    for strength in (.01, .005):
        drives = ((1, 1), (1, -1), (-1, 1), (-1, -1))
        values = [readings(packets[round(100*abs(i-j)*strength)], mask, dx,
                           i*strength, j*strength) for i, j in drives]
        signed = (values[0]-values[1]-values[2]+values[3])/(4*strength**2)
        for i, name in enumerate(("wall_law", "layers")):
            transfer = values[1][i, 1]/(2*strength)
            response = signed[i, 2]
            mass_artifact = (.25+strength**2/2)*signed[i, 0]
            print(f"n={n} speed={strength:g} derivative={name} "
                  f"inward_rate_min={min(v[i,0] for v in values):.12g} "
                  f"mixed_current_over_UV={response:.12g} "
                  f"mixed_direct_over_UV={signed[i,3]:.12g} "
                  f"envelope_momentum_over_relative={transfer:.12g} "
                  f"mixed_rate_term_over_UV={mass_artifact:.12g}")
            if abs(response-transfer-mass_artifact)>1e-11:
                raise ValueError("Four-sign energy accounting failed")


if __name__ == "__main__":
    compare(24, "mouth_feed_weak01", "mouth_redistribution_coarse01")
    compare(32, "mouth_feed_depth01", "mouth_redistribution_depth01")
