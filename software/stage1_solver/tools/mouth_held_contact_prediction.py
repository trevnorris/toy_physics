"""Infer a scalar contact stiffness from one held mouth's complex near field.

Synthetic test with local frequency-dependent contact. Each held-mouth solve
keeps both springs and the physical slip coordinate. The pair predictor uses
only the inferred contact torque response and the known static background.
"""
import numpy as np
from scipy.special import spherical_jn

from mouth_pair_translating_cuts import AMPLITUDE, MODULUS, SPEED, Pair, hankel, polynomial
from mouth_pair_contact_motion import experiment


def held_measurement(k, omega_reference, loss, l):
    # Input a regular harmonic normalized to one at R=1. Sample the total
    # azimuthal angle on a stationary sphere outside the held mouth.
    omega = k*SPEED
    outer = 1.55
    g1 = h = 16*MODULUS
    mass, drag = h/omega_reference**2, loss*h/omega_reference
    def shape(incident,r,derivative=False):
        radial = spherical_jn if incident else hankel
        factor = k if derivative else 1.
        return factor*radial(l,k*r,derivative)/radial(l,k)
    di = 1+shape(True,1,True)
    do = 1+shape(False,1,True)
    retained = np.array([[g1-MODULUS*do,-g1],
                         [-g1,g1+h-mass*omega**2-1j*omega*drag]])
    source = np.array([MODULUS*di-g1,g1])
    scatter, _ = np.linalg.solve(retained,source)
    readout = shape(True,outer) + scatter*shape(False,outer)
    observed = (readout-shape(True,outer))/shape(False,outer)
    # Torque balance gives G/C in units R=1, independent of the mechanism
    # used to generate the sample. No source displacement in this calibration.
    inferred = (di+do*observed)/(1+observed)
    return observed, inferred


def predict(pair,k,g):
    ell = np.arange(1,pair.lmax+1)
    tangential = np.stack([polynomial(l,pair.u)[1] for l in ell])
    tangential *= ((2*ell+1)/(2*ell*(ell+1)))[:,None]*pair.weights
    project = lambda field: tangential@field[1]+(pair.pr@field[0])*(g/(ell*(ell+1)))[:,None]
    mat,rhs=[],[]
    for center,x in zip(pair.centers,(AMPLITUDE,-AMPLITUDE)):
        mat.append(project(pair.field(center,1,k)))
        step=.0001
        gradient = (pair.background(center,step)-pair.background(center,-step))/(2*step)
        rhs.append(-x*project(gradient[:,:,None])[:,0])
    coeff=np.linalg.solve(np.vstack(mat),np.concatenate(rhs))
    return pair.outgoing(coeff,k)


def main():
    pair = Pair(np.array([-2.,2.]),10,8.)
    for ratio in (.5, .95, np.sqrt(18/17), 1.15):
        k=.06*ratio
        for loss in (0.,.02):
            scatter,g = held_measurement(k,.06*SPEED,loss,1)
            _,g2 = held_measurement(k,.06*SPEED,loss,2)
            predicted = predict(pair,k,g)
            full = experiment(pair.centers,[AMPLITUDE,-AMPLITUDE],k,ratio,loss)[0]
            relative = abs(predicted/full-1)
            print(f"ratio={ratio:.9f} loss={loss:g} scatter_real={scatter.real:.9g} "
                  f"scatter_imag={scatter.imag:.9g} g_real={g.real:.9g} g_imag={g.imag:.9g} "
                  f"dipole_quadrupole_difference={abs(g-g2):.3g} pair_power={predicted:.9g} "
                  f"relative_prediction_error={relative:.3g}")
            assert abs(g-g2)<1e-10 and relative<1e-5


if __name__ == "__main__":
    main()
