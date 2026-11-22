"""Normal acoustic pressure transfer through an immobile phase link.

The two neutral half-spaces share n0=c=delta=1. The held plane carries the
same quadratic phase contact as the static test, with no independent normal
mechanical transmission. We use nondispersive, small-amplitude acoustics; the
sigma=4 Gaussian pressure shape matches the earlier neutral carrier pulse.
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import erfcx

from mouth_free_front_wave_head import emit
from mouth_puncture_phase_barrier import link_resistance


def fixed_transmission(k, resistance):
    if np.isinf(resistance):
        return 0.
    return 1./(1.+(.5*k*resistance)**2)


def scattering(k, resistance):
    # Unit incident pressure with exp(-i omega t). Pressure=i omega theta,
    # normal current=(p_inc-p_ref)/Z=p_tra/Z, Z=n0 c=1 here.
    # theta_plus-theta_minus = R * j_n.
    if np.isinf(resistance):
        return 1.+0j, 0.+0j
    matrix = np.array([[1.,1.],[-1.,1.-1j*k*resistance]])
    return np.linalg.solve(matrix, [1.,1.])


def gaussian_packet(strength, width):
    resistance = link_resistance(strength)
    def integrand(k):
        return np.exp(-(width*k)**2) * fixed_transmission(k,resistance)
    direct = quad(integrand,0.,np.inf,epsabs=4e-13)[0] * 2*width/np.sqrt(np.pi)
    if resistance == 0.:
        closed=1.
    elif np.isinf(resistance):
        closed=0.
    else:
        x = 2*width/resistance
        closed=np.sqrt(np.pi)*x*erfcx(x)
    r,t = scattering(1./width,resistance)
    emit(test='fixed_flat_acoustic_link',strength=strength,width=width,
         resistance=None if np.isinf(resistance) else resistance,
         fixed_support=True, packet_transmitted_fraction=direct,
         analytic_packet_fraction=closed, probe_k=1./width,
         transmitted_at_probe=abs(t)**2, reflected_at_probe=abs(r)**2,
         probe_impedance_prediction=fixed_transmission(1./width,resistance))
    if (abs(direct-closed)>2e-11 or abs(abs(r)**2+abs(t)**2-1)>1e-12 or
            abs(abs(t)**2-fixed_transmission(1./width,resistance))>1e-12):
        raise ArithmeticError('Acoustic scattering or packet energy closure failed')


def impermeable_sheet(k, areal_inertia):
    # Separate freely moving nonporous planar sheet. Equal velocities at the
    # two faces and jump in pressure = areal inertia times sheet acceleration.
    # No tension restoring force in the uniform plane-translation mode.
    matrix=np.array([[1.,1.,-0.],[-1.,0.,-1.],[1.,-1.,1j*k*areal_inertia]])
    reflected,transmitted,velocity=np.linalg.solve(matrix,[1.,-1.,-1.])
    emit(test='impermeable_moving_sheet', wavenumber=k, example_areal_inertia=areal_inertia,
         incident_pressure=1., lab_sheet_speed_magnitude=abs(velocity),
         relative_normal_flux_magnitude=abs(transmitted-velocity),
         reflected_fraction=abs(reflected)**2,transmitted_fraction=abs(transmitted)**2)
    if abs(transmitted-velocity)>1e-12 or abs(abs(reflected)**2+abs(transmitted)**2-1)>1e-12:
        raise ArithmeticError('Moving-sheet mass or energy balance failed')


if __name__=='__main__':
    for strength in (0.,.9,.9999,1.):
        gaussian_packet(strength,4.)
    for areal_inertia in (0.,1.):
        impermeable_sheet(.25,areal_inertia)
