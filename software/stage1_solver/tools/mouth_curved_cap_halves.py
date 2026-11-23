"""Two independent half-interface ports sharing the half-order displacement.

At the center one uses the actual perturbed trace. At the outer face the traces
are extrapolated back to the center through a homogeneous lead. Therefore the
right half can be added directly only to a solver retaining the full, mobile
left layer. Joining both half ports and freeing their common phase/density
provides an independent check of the whole diffuse front.
"""
from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import splu

from mouth_curved_cap_ports import (
    background, evolution_oblique, order_impedance, uniform,
)
from mouth_poly_shallow_release import emit


@lru_cache(maxsize=128)
def half_port(omega,j,kt=0.,side=1,extent=.48,spacing=.0012):
    if side not in (-1,1):
        raise ValueError('Use -1 upstream or +1 downstream')
    fields = background(j)
    g0 = fields(0.)[-1].item()
    area = max(81,round(extent/spacing)+1)
    points = np.linspace(-extent,0,area) if side<0 else np.linspace(0,extent,area)
    homogeneous = evolution_oblique(np.array([0.]),uniform,omega,j,kt)[0,:4,:4]
    eye = np.eye(6,dtype=complex)
    maps = [eye[[2,0,3,1]].copy(),eye[[2,0,3,1]].copy()]
    outer = (side+1)//2
    center = 1-outer
    maps[outer][:,:4] = expm(-side*extent*homogeneous)[[2,0,3,1]]
    left = np.zeros((6,6),complex)
    right = np.zeros_like(left)
    bc = [left,right]
    left[[0,2]] = maps[0][:2]
    right[[1,3]] = maps[1][:2]
    bc[outer][4] = eye[5]-order_impedance(omega,j,kt,side)*eye[4]
    bc[center][5] = eye[4]
    rhs = np.zeros((6,5),complex)
    rhs[:4,:4] = np.eye(4)
    rhs[5,4] = -g0
    # One linear banded collocation, including the six physical boundary
    # equations. Keep true face coordinates at the joined material interface.
    h=np.diff(points)[:,None,None]
    at=evolution_oblique(points,fields,omega,j,kt)
    mid=evolution_oblique((points[1:]+points[:-1])/2,fields,omega,j,kt)
    lo=-eye-h/6*(at[:-1]+4*mid@(.5*eye+h/8*at[:-1]))
    hi=eye-h/6*(at[1:]+4*mid@(.5*eye-h/8*at[1:]))
    mat=sparse.hstack((sparse.block_diag(lo),sparse.csr_matrix((6*(area-1),6))))
    mat+=sparse.hstack((sparse.csr_matrix((6*(area-1),6)),sparse.block_diag(hi)))
    mat=sparse.vstack((mat,sparse.hstack((left,sparse.csr_matrix((6,6*(area-2))),right))),format='csc')
    loads=np.zeros((6*area,5),complex)
    loads[-6:]=rhs
    data=splu(mat).solve(loads)
    traces=[data[:6],data[-6:]]
    fl=maps[0][2:]@traces[0]
    fr=maps[1][2:]@traces[1]
    return np.stack((-fl[0],fr[0],-fl[1],fr[1],side*g0*traces[center][5]))


def run():
    import json
    from mouth_poly_shallow_release import ROOT
    folder=ROOT/'_scratch/s11c/mouth_curved_cap01/ports'
    guard=json.loads((folder/'outcome.json').read_bytes())
    if guard['exitCode'] or not guard['limitsVerified'] or guard['childOutcome']['exitCode']:
        raise ArithmeticError('Previous full-interface guard did not complete')
    cases=[json.loads(line) for line in (folder/'stdout').read_bytes().splitlines()]
    omega,j=cases[0]['omega'],cases[0]['j']
    for kt in (0.,1.8,4.33):
        lhs=half_port(omega,j,kt,-1)
        rhs=half_port(omega,j,kt,1)
        matrix=np.zeros((7,7),complex)
        for block,indices in [(lhs,[0,5,2,6,4]),(rhs,[5,1,6,3,4])]:
            matrix[np.ix_(indices,indices)]+=block
        x=np.linalg.solve(matrix[5:,5:],matrix[5:,:5])
        combined=matrix[:5,:5]-matrix[:5,5:]@x
        # Conserves the same as-run direct whole-interface solution.
        row,=[r for r in cases if r['kind']=='sheet' and r['kt']==kt and r['extent']==.48 and r['n']==4]
        saved=np.asarray(row['external_port'])
        full=saved[...,0]+1j*saved[...,1]
        relative=np.linalg.norm(combined-full)/np.linalg.norm(full)
        x_err=np.linalg.norm(combined[:,-1]-full[:,-1])/np.linalg.norm(full[:,-1])
        emit(kind='joined_half_layers',kt=kt,omega=omega,j=j,relative_error=relative,
             driver_column_relative_error=x_err,
             left_hermitian=np.linalg.norm(lhs-lhs.conj().T)/np.linalg.norm(lhs),
             right_hermitian=np.linalg.norm(rhs-rhs.conj().T)/np.linalg.norm(rhs))
        if max(relative,x_err)>1.e-4:
            raise ArithmeticError('The joined half layers disagree with the direct sheet')


if __name__=='__main__':
    run()
