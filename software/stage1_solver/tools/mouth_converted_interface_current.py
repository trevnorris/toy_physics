"""Normal-current test of a flat neutral brane-to-bulk interface with one continuous phase.

The tested ambiguity is the density derivative of the order-gradient energy.
Both laws retain the native carrier current j=n grad(theta), EOS, and bare
order tension at homogeneous reference density. A half-order sheet is centered
by the fixture. Reversible order relaxation and the specified mean external
sheet potential apply throughout. All currents refer to unit interface area.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.special import expit

from mouth_curved_site_front import sha
from mouth_free_order_diffuse import SIGMA
from mouth_surface_carrier_reflection import PITCH, sheet_mean, stationary
from neutral_patch_sites import emit


def records(path, run):
    r = json.loads(path.read_text())['guarded_runs'][run]
    assert r['limits_verified'] and r['child_exit']==r['guard_exit']==0
    return [json.loads(line) for line in r['literal_stdout'].splitlines()]


def run(args):
    setups = records(args.registration,'collective')
    setup, = [r for r in setups if r['test']=='setup']
    band_record = records(args.impact,'incidence')
    source, = [r for r in band_record if r['test']=='preparation']
    width = source['delta']
    nb, height = setup['background_density'], args.height
    kap0 = 6*SIGMA*width
    def potential(x):
        return setup['energy_scale']*sheet_mean(x/PITCH,setup['pitch'])[0]

    mesh = np.unique(np.r_[np.linspace(0,.06,181),np.linspace(.06,1.,141),
                           np.linspace(1.,height,1+int(20*(height-1)))])
    base = stationary(mesh,potential,nb,args.tol)
    x = base.x
    c = expit(-x/width)
    init = np.vstack((base.y, c, -kap0*c*(1-c)/width,np.zeros_like(x)))
    emit(test='setup', half_extent=height, tolerance=args.tol, bulk_density=nb, width=width,
         kappa_reference=kap0, sigma_reference=SIGMA, specific_sheet_energy=setup['energy_scale'],
         pitch=setup['pitch'], bare_compression_modulus=nb**5, homogeneous_static_shear_bulk=0.,
         position_reference='fixture fixes sheet and half order at z=0; opposite reservoirs are in opposite order states',
         registration_record=str(args.registration), registration_sha256=sha(args.registration),
         impact_record=str(args.impact), impact_sha256=sha(args.impact))
    reference = None
    for eta in (0.,.5):
        dk = kap0*eta/nb
        def ode(z,y,j):
            f,df,c,q,theta = y
            kappa = kap0*(1-eta)+dk*f*f
            dc = q/kappa
            phase_grad = j[0]/f**2
            chemical = potential(z)+(f**8-nb**4)/4 + .5*dk*dc**2
            return np.array([df,(2*chemical+phase_grad**2)*f,dc,
                             kap0/width**2*c*(1-c)*(1-2*c),phase_grad])

        zero = None
        for bias in (0.,.002,-.002,.001,-.001):
            def bc(left,right,j):
                return np.array([left[1],left[2]-.5,left[4],
                                 right[0]-np.sqrt(nb),right[2],right[4]-bias/2])
            if zero is None:
                guess_x = x if reference is None else reference.x
                guess_y = init if reference is None else reference.y
                starting_current = 0.
            else:
                guess_x, guess_y = zero.x, zero.y.copy()
                guess_y[4] = bias*guess_x/(2*height)
                starting_current = bias/total_length
            sol = solve_bvp(ode,bc,guess_x,guess_y,p=[starting_current],
                            tol=args.tol,bc_tol=1e-11,max_nodes=35000)
            if not sol.success:
                raise ArithmeticError(sol.message)
            z = np.r_[np.linspace(0,.05,251),np.linspace(.06,1.,100),np.linspace(1.1,height,100)]
            f, df, c, q, theta = sol.sol(z)
            j = float(sol.p[0])
            # Check the independent stationary normal-stress identity.
            # Hamiltonian per unit area includes the conserved phase term;
            # its derivative is minus n V' for an imposed external sheet.
            kappa = kap0*(1-eta)+dk*f*f
            hamiltonian = (.5*df**2+.5*j*j/f**2+.5*q*q/kappa
                -(f**10-5*nb**4*f**2+4*nb**5)/20
                -kap0/(2*width**2)*(c*(1-c))**2-potential(z)*f*f)
            def integrand(s):
                grad = setup['energy_scale']/PITCH*sheet_mean(s/PITCH,setup['pitch'])[1]
                return sol.sol(s)[0]**2 * grad
            force = sum(quad(integrand,a,b,epsabs=1e-10)[0] for a,b in ((0,.01),(.01,.05),(.05,.2)))
            stress_error = hamiltonian[-1]-hamiltonian[0]+force
            count = np.min(f)
            emit(test='signed', eta=eta,bias=bias,current=j,mid_density=f[0]**2,
                 measured_mid_chi_gradient=q[0]/kappa[0],
                 normal_stress_residual=float(stress_error), nodes=len(sol.x))
            assert count>0 and abs(stress_error)<5e-6 and np.max(np.diff(c))<1e-9
            if bias==0:
                zero = sol
                if eta==0:
                    reference=sol
                integrations=((0,.01),(.01,.05),(.05,1.),(1.,height))
                total_length = 2*sum(quad(lambda s:1/sol.sol(s)[0]**2,a,b,epsabs=3e-9)[0]
                                     for a,b in integrations)
                intrinsic = total_length-2*height/nb
                emit(test='linear', eta=eta, chamber_slope=1/total_length, intrinsic_slope=1/intrinsic,
                     excess_phase_length=intrinsic*nb, density_ratio=f[0]**2/nb,
                     fixed_uniform_density_reference_tension=SIGMA, width_from_mid_gradient=-.25*kappa[0]/q[0])
                assert abs(j)<1e-12
            else:
                assert np.isclose(j/bias,1/total_length,rtol=2e-6)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',type=Path,default=Path('research/pde_ledger_v3/_measurements/mouth_surface_registration.json'))
    parser.add_argument('--impact',type=Path,default=Path('research/pde_ledger_v3/_measurements/mouth_attachment_impact.json'))
    parser.add_argument('--height',type=float,default=18.)
    parser.add_argument('--tol',type=float,default=2e-7)
    args=parser.parse_args()
    if args.height<10 or args.tol<=0:
        parser.error('Use half-extent >= 10 and positive residual tolerance')
    run(args)
