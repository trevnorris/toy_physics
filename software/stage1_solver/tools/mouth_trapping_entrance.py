"""Locally planar, real-frequency shear refraction from the bonded mouth into its flat feed."""
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legder, legval

from mouth_exit_witness import inputs
from mouth_finite_bore_support import graph_values
from mouth_poly_bonded_core import core_stretches
from mouth_poly_shallow_release import Trial, emit
from mouth_sealed_shape_relax import Meridian


def strata(outer,cut,q,p,angular_length):
    # At pure angular tangential slowness with p_r=0, |p_h| = angular_length/r.
    # The transmitted threshold is r=angular_length in current coordinates.
    intervals = []
    lo,hi = q*angular_length,min(cut,angular_length)
    if lo<hi:
        intervals.append(['core',lo,hi])
    limiting_label = max(0.,angular_length-outer,-np.log(angular_length/outer)/p)
    if limiting_label<.5:
        intervals.append(['collar',cut,outer*np.exp(-p*limiting_label)])
    return intervals


def run(radius,n):
    root=Path('_scratch/s11c')
    tr=Trial(radius,n,bonded_core=True)
    sh,chi,_,digest,status,*_=inputs(root/f'mouth_exit_release01/r{radius}n{n}/endpoint.npz',radius,n,tr)
    shape=Meridian(radius,tr.x+sh[:4])
    c=tr.co.copy()
    c[0]+=sh[4]
    lp=legval(-1.,c)
    p=np.exp(lp)
    cut=radius*np.exp(-p/2)
    q=core_stretches(shape,c,*tr.core_reference,np.array([0.]),np.array([cut]))[0].item()
    folder=root/f'mouth_trapping_tail01/r{radius}n{n}'
    guard=json.loads((folder/'outcome.json').read_bytes())
    raw=(folder/'stdout').read_bytes()
    rows=[json.loads(x) for x in raw.splitlines()]
    if (status!='balanced' or guard['exitCode'] or guard['childOutcome']['exitCode'] or
            not guard['limitsVerified'] or rows[0]['endpoint_sha256']!=digest or np.min(chi[0])<=0):
        raise ValueError('Material, order or spectrum preparation mismatch')

    # Independent spatial point test directly on the meridional tensor.
    rr=radius*(np.arange(10000)+.5)/10000
    collar=rr>cut
    t=np.zeros_like(rr)
    rad,mix,ax,jac=core_stretches(shape,c,*tr.core_reference,t,rr)
    axial=np.full_like(rr,ax)
    hoop=rad.copy()
    a=-np.log(rr[collar]/radius)/p
    _,rad[collar],mix[collar],hoop[collar],_,_=graph_values(
        a,t[collar],shape,lp,2*legval(-1.,legder(c)))
    axial[collar]=-shape.length
    brr,brw,bww=rad**2+mix**2,mix*axial,axial**2
    schur=brr-brw**2/bww
    if np.min(schur)<=0 or np.max(abs(schur-rad**2))>1e-11:
        raise ArithmeticError('Invalid normal dispersion discriminant')
    # At one material tag, test a signed incoming normal group velocity, the
    # continued flat evanescent root, and both dispersion polynomials explicitly.
    rmid=radius*np.exp(-p*.25)
    h=rmid/(radius+.25)
    ph=(1+1/h)/2
    pw=-np.sqrt(1-h*h*ph*ph)/shape.length
    flat_pw=-1j*np.sqrt(ph*ph-1)  # exp(i omega p_w w) decays as w -> -infinity
    err=max(abs(shape.length**2*pw**2+h*h*ph*ph-1.),abs(flat_pw**2+ph**2-1.))
    if err>1e-12 or pw>=0 or flat_pw.imag>=0:
        raise ArithmeticError('Local transmitted root failed')
    emit(kind='entrance',radius=radius,order_mesh=n,endpoint_sha256=digest,
         frequency_scan_sha256=hashlib.sha256(raw).hexdigest(),cut_radius=cut,
         bonded_core_hoop_stretch=q,collar_min_hoop_limit=cut/(radius+.5),
         inward_collar_tag=.25,tag_radius=rmid,tag_hoop_speed=h,
         tag_pure_hoop_window=[1.,1/h],example_p_h=ph,
         incident_normal_slowness=pw,incident_normal_group_velocity=shape.length**2*pw,
         flat_transmission_normal_decay_per_omega=-flat_pw.imag,dispersion_error=err,
         collar_fraction_of_open_area=1-(cut/radius)**3,
         sampled_incident_order=[float(min(chi[0])),float(max(chi[0]))])
    for degree in (6,12):
        row,=[s for s in rows if s['kind'] in ('start','angular_mode') and s['degree']==degree]
        omega=row['frequency']
        for correction in (-2.,0.):
            angular_length=np.sqrt(degree*(degree+1)+correction)/omega
            intervals=strata(radius,cut,q,p,angular_length)
            coverage=sum(hi**3-lo**3 for _,lo,hi in intervals)/radius**3
            angular=angular_length/rr
            # Both valid radicands are tested separately; p_r is fixed to zero.
            hit=(1-hoop**2*angular**2>0)&(1-angular**2<0)
            sampled=np.sum(rr[hit]**2)/np.sum(rr**2)
            if abs(sampled-coverage)>3.01e-4:
                raise ArithmeticError('Material-label and physical-coordinate windows disagree')
            emit(kind='angular_ray_estimate',radius=radius,order_mesh=n,degree=degree,
                 omega_real=omega,angular_eigenvalue_offset=correction,
                 pure_hoop_angular_length=angular_length,subluminal_window_radii=intervals,
                 open_area_fraction=coverage,independent_area_samples=sampled,
                 tangential_radial_slowness_assumed=0.)


if __name__=='__main__':
    for radius in (1,2):
        for mesh in (48,72):
            run(radius,mesh)
