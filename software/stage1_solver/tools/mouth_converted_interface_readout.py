"""Separate contact force from transverse carriage on a smooth patch; calibrate its order cost."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit

from mouth_curved_site_front import sha
from neutral_patch_sites import emit


def run(folder):
    outcome=json.loads((folder/'outcome.json').read_bytes())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    rows=[json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    setup,=[r for r in rows if r['test']=='setup']
    width,kap=setup['width'],setup['kappa_reference']
    emit(test='input',run=str(folder),stdout_sha256=sha(folder/'stdout'))
    for eta in (0.,.5):
        # Independent, clamped homogeneous-density test with no imposed particle
        # sheet potential. In a free one-dimensional order kink the two order
        # energy densities are equal. Infer its width from the measured energy.
        sigmas=[]
        for density_ratio in (.98,.99,1.,1.01,1.02):
            stiffness=kap*(1-eta+eta*density_ratio)
            well=kap/(2*width**2)
            def terms(scale):
                def gradient(z):
                    c=expit(-z/scale)
                    return .5*stiffness*(c*(1-c)/scale)**2
                def local(z):
                    c=expit(-z/scale)
                    return well*(c*(1-c))**2
                return [2*quad(term,0,60*width,epsabs=1e-12)[0] for term in (gradient,local)]
            def dilation(scale):
                grad,local=terms(scale)
                return local-grad
            scale=brentq(dilation,.8*width,1.2*width,xtol=2e-14)
            tension=sum(terms(scale))
            sigmas.append(tension)
            emit(test='clamped_density',eta=eta,density_over_reference=density_ratio,
                 positive_gradient_cost=stiffness,tension=tension,width=scale)
        emit(test='clamped_log_slope',eta=eta,
             full=(sigmas[-1]-sigmas[0])/(.04*sigmas[2]),half=(sigmas[3]-sigmas[1])/(.02*sigmas[2]))
        for row in rows:
            if row.get('test')!='signed' or row['eta']!=eta:
                continue
            # At the symmetry plane f'=0; phase_z is fixed by current. A uniform
            # tangential boost relative to the held sheet leaves n fixed when the
            # remote chemical potential takes up the common U²/2 kinetic term.
            f=np.sqrt(row['mid_density'])
            dzpsi=1j*row['current']/f
            for slip in (0.,.002,-.002,.001,-.001):
                dxpsi=1j*slip*f
                mixed_stress=(dzpsi.conjugate()*dxpsi).real
                emit(test='coherent_slip',eta=eta,bias=row['bias'],slip=slip,
                     signed_flux_per_area=mixed_stress, sheet_tangential_force_per_area=0.)
                assert abs(mixed_stress-row['current']*slip)<1e-18


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('normal_run',type=Path)
    run(p.parse_args().normal_run)
