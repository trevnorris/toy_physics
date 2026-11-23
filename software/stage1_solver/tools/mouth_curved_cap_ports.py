"""Intrinsic moving conversion sheet for coherent sound, with uniform leads removed.

Return jumps (+ receiver minus − feed) in phase, sqrt-density, normal current and
density slope, and the external force on the half-order level, as a linear map of
their four mean traces and its prescribed displacement. A homogeneous carrier
fills the acoustic leads on either side. The sheet retains no lead radiation
impedance; the surrounding domain supplies the two incoming states.
"""
import argparse
from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.linalg import splu

from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA, WIDTH, WELL, evolution, stationary_front
from mouth_poly_shallow_release import emit


@lru_cache(maxsize=96)
def background(j):
    return stationary_front(j, 16.)


def evolution_oblique(z, fields, omega, j, kt):
    mat=evolution(z,fields,omega,j)
    if kt:
        f, _, _, g=fields(z)
        n=f**2
        pi=ORDER_INERTIA*j*g/n
        k2=kt**2
        mat[:,1,0]+=k2
        mat[:,3,2]+=k2*n
        mat[:,3,4]-=k2*pi
        mat[:,5,2]-=k2*pi
        mat[:,5,4]+=k2*(KAPPA+pi*pi/n)
    return mat


def uniform(z):
    one=np.ones_like(z)
    zero=np.zeros_like(z)
    return one,zero,zero,zero


def order_impedance(omega,j,kt,side):
    roots=np.roots([KAPPA-ORDER_INERTIA*j*j, 2*ORDER_INERTIA*j*omega,
                    2*WELL+KAPPA*kt**2-ORDER_INERTIA*omega**2])
    ko,=[k for k in roots if side*k.imag>0]
    return 1j*((KAPPA-ORDER_INERTIA*j*j)*ko+ORDER_INERTIA*j*omega)


@lru_cache(maxsize=320)
def sheet(omega, j, kt=0., extent=.48, refinement=2):
    """Mean-to-jump scattering form, without propagating or evanescent lead loads.

    All acoustic field traces refer to z=0, extended from ±extent using the uniform
    carrier. For finite j the physical stationary density relaxes on a healing
    length, so extent is to be tested, in particular for growing evanescent inputs.
    """
    fields=background(j)
    half=np.linspace(0,extent,round(extent/WIDTH*6)*refinement+1)
    z=np.r_[-half[::-1],half]
    mid=len(half)-1
    starts=np.delete(np.arange(z.size-1),mid)
    step=np.diff(z)[starts,None,None]
    mat=evolution_oblique(z,fields,omega,j,kt)
    mhalf=evolution_oblique((z[starts]+z[starts+1])/2,fields,omega,j,kt)
    eye=np.eye(6)
    left=-eye-step/6*(mat[starts]+4*mhalf@(.5*eye+step/8*mat[starts]))
    right=eye-step/6*(mat[starts+1]+4*mhalf@(.5*eye-step/8*mat[starts+1]))
    ids=6*np.arange(z.size)[:,None]+np.arange(6)
    rows=[]; cols=[]; values=[]
    def block(row,col,matrix):
        rows.append(np.broadcast_to(row[..., :,None],matrix.shape).ravel())
        cols.append(np.broadcast_to(col[...,None,:],matrix.shape).ravel())
        values.append(matrix.ravel())
    block(ids[:len(starts)],ids[starts],left)
    block(ids[:len(starts)],ids[starts+1],right)
    row=6*len(starts)
    rhs=np.zeros((len(z)*6,5),complex)
    Au=evolution_oblique(np.array([0.]),uniform,omega,j,kt)[0,:4,:4]
    maps=[]
    for side,node in [(-1,0),(1,-1)]:
        extrapolate=expm(-side*extent*Au)[[2,0,3,1]]
        virtual=np.zeros((4,6),complex)
        virtual[:,:4]=extrapolate
        maps.append(virtual)
        block(np.arange(row,row+4),ids[node],.5*virtual)
    rhs[row:row+4,:4]=np.eye(4)
    row+=4
    for side,node in [(-1,0),(1,-1)]:
        line=eye[5:6]-order_impedance(omega,j,kt,side)*eye[4:5]
        block(np.array([row]),ids[node],line)
        row+=1
    block(np.arange(row,row+5),ids[mid],eye[:5])
    block(np.arange(row,row+5),ids[mid+1],-eye[:5])
    row+=5
    block(np.array([row]),ids[mid],eye[4:5])
    g0=fields(0.)[-1].item()
    rhs[row,4]=-g0
    A=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),
                        shape=(len(rhs),len(rhs))).tocsc()
    profiles=splu(A).solve(rhs).reshape(len(z),6,5)
    traces=[maps[0]@profiles[0],maps[1]@profiles[-1]]
    out=np.empty((5,5),complex)
    out[:4]=traces[1]-traces[0]
    out[4]=g0*(profiles[mid+1,5]-profiles[mid,5])
    return out


def port(omega,j,kt=0.,extent=.48,refinement=2):
    """External dynamic stiffness, trace order (θ−,θ+,δsqrt(n)−,δsqrt(n)+,X).

    The two-field traces use the same coordinate, with outward forces (-j−,+j+,
    -a'−,+a'+,F). Near singular density-jump charts use sheet() with a multiplier.
    """
    G=sheet(omega,j,kt,extent,refinement)
    eye=np.eye(5)
    mean=.5*np.stack((eye[0]+eye[1],eye[2]+eye[3]))
    jump=np.stack((eye[1]-eye[0],eye[3]-eye[2]))
    known=G[:2,:2]@mean+G[:2,4:5]@eye[4:5]
    mean_flux=np.linalg.solve(G[:2,2:4],jump-known)
    outputs=G@np.vstack((mean,mean_flux,eye[4]))
    minus=-mean_flux+outputs[2:4]/2
    plus=mean_flux+outputs[2:4]/2
    return np.stack((minus[0],plus[0],minus[1],plus[1],outputs[4]))


def plane_modes(omega,j,kt,side):
    roots=np.roots([.25,0,1+kt**2/2-j*j,2*omega*j,kt**2*(1+kt**2/4)-omega**2])
    modes=[]
    for k in roots:
        k2=k*k+kt*kt
        intrinsic=omega-j*k
        theta=2*(1+k2/4)/(1j*intrinsic)
        state=np.array([theta,1,2*j+1j*k*theta,1j*k])
        if abs(k.imag)<1.e-8:
            group=(j+k*(1+k2/2)/intrinsic).real
            if group*side>0:
                power=side*2*group*omega*intrinsic.real/k2.real
                modes.append((0,k,state,power))
        elif k.imag*side>0:
            modes.append((1,k,state,0.))
    return sorted(modes,key=lambda x:(x[0],abs(x[1].imag)))


def scatter(G,omega,j,kt):
    """Restore homogeneous outgoing conditions and scatter both incident modes plus a wall drive."""
    modes=[plane_modes(omega,j,kt,side) for side in (-1,1)]
    if any(len(seq)!=2 or seq[0][0]!=0 for seq in modes):
        return None
    channels=[np.stack([v[2] for v in seq],axis=1) for seq in modes]
    L=np.zeros((4,4),complex); R=np.zeros_like(L)
    L[:,:2]=channels[0]; R[:,2:]=channels[1]
    ins=[modes[1][0],modes[0][0]]
    left=np.zeros((4,3),complex);right=np.zeros_like(left)
    left[:,0]=ins[0][2];right[:,1]=ins[1][2]
    X=np.array([0.,0.,1.])
    coeff=np.linalg.solve(R-L-.5*G[:4,:4]@(L+R),
                          -right+left+.5*G[:4,:4]@(right+left)+G[:4,4,None]*X)
    mean=(L@coeff+R@coeff+left+right)/2
    force=G[4,:4]@mean+G[4,4]*X
    mout=np.sqrt([modes[0][0][3],modes[1][0][3]])
    incoming=np.sqrt([v[3] for v in ins])
    outgoing=mout[:,None]*coeff[[0,2]]
    S=outgoing[:,:2]/incoming
    Qin=np.diag(np.r_[incoming*incoming,0.])
    work=.25*(np.outer(1j*omega*X,force)+np.outer(force.conj(),-1j*omega*X))
    balance=outgoing.conj().T@outgoing-Qin-work
    return dict(S=S,core_force=force[2],force_from_incoming=force[:2], unitary_error=np.linalg.norm(S.conj().T@S-np.eye(2)),
                power_error=np.max(abs(balance)), receiver_power=abs(outgoing[1,2])**2,
                feed_power=abs(outgoing[0,2])**2)


def cplx(v):
    v=np.asarray(v)
    return np.stack((v.real,v.imag),axis=-1).tolist()


def run():
    import hashlib
    import json
    from mouth_poly_shallow_release import ROOT
    root=ROOT/'_scratch/s11c/mouth_breathing_mobile_patch01/matched_r2'
    raw=(root/'stdout').read_bytes()
    guard=json.loads((root/'outcome.json').read_bytes())
    if guard['exitCode'] or not guard['limitsVerified'] or guard['childOutcome']['exitCode']:
        raise ArithmeticError('Unvalidated acoustic reference')
    rows=list(map(json.loads,raw.splitlines()))
    setup=rows[0]
    omega=setup['frequency']; j=setup['measured_mean_flux']
    emit(kind='prepared',omega=omega,j=j,reference_sha=hashlib.sha256(raw).hexdigest(),
         shape_sha=setup['shape_hash'])
    qdiag=np.diag([1,-1,1,-1,1])
    for kt in (0.,1.,1.8,4.33):
        reference=None
        for extent,n in [(.32,3),(.48,2),(.48,4),(.72,3),(1.2,3)]:
            g=sheet(omega,j,kt,extent,n)
            stiffness=port(omega,j,kt,extent,n)
            test=scatter(g,omega,j,kt)
            herm=np.linalg.norm(stiffness-stiffness.conj().T)/np.linalg.norm(stiffness)
            row=dict(kind='sheet',kt=kt,extent=extent,n=n,mean_to_jump=cplx(g),external_port=cplx(stiffness),
                     density_phase_jump_condition=np.linalg.cond(g[:2,2:4]),hermitian_defect=herm)
            if test is not None:
                row.update(S=cplx(test['S']),force_from_incoming=cplx(test['force_from_incoming']),
                           translation_force=cplx(test['core_force']),unitarity=test['unitary_error'],
                           power_error=test['power_error'],feed=test['feed_power'],receiver=test['receiver_power'])
            if extent==.48 and n==4:
                neg=port(omega,-j,kt,extent,n)
                # Phase is odd, sqrt-density and the half-order displacement even under time reversal.
                parity=np.diag([-1,-1,1,1,1])
                row['reversal_relative_defect']=np.linalg.norm(stiffness-parity@neg.T@parity)/np.linalg.norm(stiffness)
            if kt==0 and test is not None:
                previous,=[r for r in rows if r['kind']=='two_sided_patch' and r['normal_current']==j and r['mesh_refinement']==6]
                row['diffuse_force_mismatch']=abs(test['core_force']-complex(previous['force_re'],previous['force_im']))
            emit(**row)
            if herm>2e-4 or test is not None and max(test['unitary_error'],test['power_error'])>2e-4:
                raise ArithmeticError('Flux or lossless intrinsic stiffness failed')
    # Held, zero-current neutral front at normal incidence gives an independent reciprocal test.
    g=sheet(omega,0.,0.,.48,4)
    res=scatter(g,omega,0.,0.)
    emit(kind='neutral',transmit_power=abs(res['S'][1,0])**2,
         reflection_power=abs(res['S'][0,0])**2, equal_side_emission=abs(res['feed_power']/res['receiver_power']-1),
         S=cplx(res['S']),unitarity=res['unitary_error'])
    if abs(res['S'][1,0]-res['S'][0,1])>2e-6 or abs(res['feed_power']/res['receiver_power']-1)>2e-6:
        raise ArithmeticError('Neutral two-port symmetry failed')


if __name__=='__main__':
    run()
