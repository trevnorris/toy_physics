#!/usr/bin/env python3
"""Small noncommuting arithmetic controls; no physical solve or quadrature."""
import json
from pathlib import Path
import numpy as np
import S11c_d_continuum_response as c


def main():
    rng=np.random.default_rng(11919)
    matrix={g:rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)) for g in c.G}
    matrix[(0,0)]+=8*np.eye(4)
    expected={g:rng.normal(size=(4,3))+1j*rng.normal(size=(4,3)) for g in c.G}
    # Independent untruncated convolution, then retain the requested rectangle.
    rhs={g:sum((matrix[a]@expected[b] for a in c.G for b in c.G
                if tuple(x+y for x,y in zip(a,b))==g),np.zeros((4,3),complex)) for g in c.G}
    result=c.solve(matrix,rhs)
    difference={g:result['coefficients'][g]-expected[g] for g in c.G}
    c.f.require(c.boundary.norm(difference)<1e-12,'noncommuting full mixed recursion')
    seed={g:rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)) for g in c.G}
    root={g:a+a.conj().T for g,a in seed.items()};root[(0,0)]+=12*np.eye(3)
    metric={g:sum((root[a]@root[b] for a in c.G for b in c.G
                  if tuple(x+y for x,y in zip(a,b))==g),np.zeros((3,3),complex)) for g in c.G}
    square,residual=c.sqrt_series(metric)
    square_difference=c.boundary.subtract(square,root)
    c.f.require(c.boundary.norm(square_difference)<1e-11,'noncommuting current square root')
    reversed_metric=dict(metric)
    reversed_metric[(1,1)]=metric[(1,1)]+root[(1,0)]@root[(0,1)]-root[(0,1)]@root[(1,0)]
    mutation=c.boundary.norm(c.boundary.subtract(c.J.multiply(square,square),reversed_metric))
    c.f.require(mutation>1e-6,'mixed-order matrix mutation response')
    # Field and phase maps generally do not commute. The field-coordinate phase
    # is conjugated by the same basis map, never applied in the old coordinates.
    field=np.array([[2,1],[0,1]],complex);phase=np.array([[1,2j],[3,1]],complex);s=np.array([[1,2],[3,5]],complex)
    joined=field@phase@s;separate=(field@phase@np.linalg.inv(field))@(field@s)
    phase_residual=c.boundary.norm(joined-separate);wrong=c.boundary.norm(phase@field@s-joined)
    c.f.require(phase_residual<1e-12 and wrong>0,'field-coordinate phase and sign/order control')
    report={'scriptSha256':c.f.digest(Path(c.__file__)),'focusedSha256':c.f.digest(Path(__file__)),
        'coefficientDifference':c.boundary.norm(difference),'independentDifference':c.boundary.norm(result['independentDifference']),
        'mixedForcingMutation':c.boundary.norm(result['mixedForcingMutation']),
        'squareRootDifference':c.boundary.norm(square_difference),'squareResidual':c.boundary.norm(residual),
        'mixedOrderingMutation':mutation,'phaseCoordinateResidual':phase_residual,'wrongPhaseCoordinateResponse':wrong,
        'scope':'Noncommuting finite-matrix arithmetic controls; no physical response.'}
    c.f.save(c.f.M/'S11c_d_continuum_response_focused.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':main()
