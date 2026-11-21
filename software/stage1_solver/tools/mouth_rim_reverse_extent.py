"""Locate reversed crossing in the trial front interior from its solved curvature."""
import json
import numpy as np
from scipy.sparse.linalg import spsolve

from mouth_advected_rim import Sphere, archive_drive, EVIDENCE, RATE


def readout(sphere, drive, K):
    data = sphere.prepare(drive)
    matrix = K*data.diffusion+data.drift
    active=sphere.free.ravel()
    c=np.zeros(sphere.indices.size)
    c[active]=spsolve(matrix[active][:,active],data.rhs[active])
    c=c.reshape(sphere.indices.shape)
    ell = np.arange(sphere.modes)
    r=sphere.r[:,None]
    curvature=(sphere.basis(sphere.r,2)@c + 2*(sphere.db@c)/r
               -(sphere.br@c)*(ell*(ell+1))/r**2)@sphere.p.T
    q=-K*curvature
    signed=np.sum(q*sphere.w)
    reversed_rate=-np.sum(np.minimum(q,0)*sphere.w)
    r_end = []
    for column in q.T:
        positive=np.flatnonzero(column>0)
        i=positive[-1]
        if i==len(column)-1 or sphere.r[i]<.7:
            raise ValueError('Edge change of sign is unresolved')
        r_end.append(sphere.r[i]+(sphere.r[i+1]-sphere.r[i])*column[i]/(column[i]-column[i+1]))
    # Check the sign by the conservative projected current directly inside
    # the computed negative band, independent of differentiating h twice.
    def test_radius(value):
        mu=sphere.mu
        def current(r,theta):
            m=np.cos(theta)
            _,ja,jt=drive(r*np.ones_like(m),m)
            return ja*m+jt*np.sin(theta), -ja*np.sin(theta)+jt*m
        theta=np.arccos(mu)
        eps=1e-5
        jr,jth=current(value,theta)
        div=(((value+eps)**2*current(value+eps,theta)[0]-(value-eps)**2*current(value-eps,theta)[0])/(2*eps*value**2)
             +(np.sin(theta+eps)*current(value,theta+eps)[1]-np.sin(theta-eps)*current(value,theta-eps)[1])/(2*eps*value*np.sin(theta)))
        h=sphere.basis(value)@c@sphere.p.T
        hr=sphere.basis(value,1)@c@sphere.p.T
        hth=-(sphere.basis(value)@c@sphere.dp.T)*sphere.sin/value
        q0=drive(np.full_like(mu,value),mu)[0]*RATE/data.projected_rate
        q_transport=q0-h*div-jr*hr-jth*hth
        return [float(q_transport.min()),float(q_transport.max())]
    interior=test_radius(.99)
    edge_slope=-(sphere.basis(1,1)@c)@sphere.p.T
    needed=-edge_slope*data.edge_jr
    result=dict(signed_rate=signed, reverse_rate=reversed_rate, reverse_over_net=reversed_rate/signed,
                outer_zero_radius=[min(r_end),max(r_end)], q_at_099_from_continuity=interior,
                required_edge_w_current=[needed.min(),needed.max()])
    print(json.dumps(result,allow_nan=False),flush=True)
    assert abs(signed-RATE)<.001 and interior[-1]<0 and reversed_rate>0
    return result


if __name__=='__main__':
    run=json.loads(EVIDENCE.read_text())['runs']['mouth_front_pressure_dipole02']
    assert run['guard_exit']==run['child_exit']==0 and run['limits_verified']
    for line in run['stdout_lines']:
        row=json.loads(line)
        if row['grid'] != [96,6,24]:
            continue
        alone,pair=archive_drive(row)
        states=[('paired',pair)] if row['d']==16 else [('alone',alone),('paired',pair)]
        for state,drive in states:
            for K in (.4,.8,1.6):
                previous=None
                for n,m,a in ((48,9,120),(96,9,120)):
                    sphere=Sphere(n,m,a)
                    print(f'state={state} d={row["d"]} K={K} elements={n} modes={m}',flush=True)
                    result=readout(sphere,drive,K)
                    if previous:
                        np.testing.assert_allclose(result['reverse_rate'],previous['reverse_rate'],rtol=.02)
                    previous=result
