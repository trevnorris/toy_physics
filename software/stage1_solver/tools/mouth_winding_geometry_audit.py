"""Independent two-sided winding geometry and continuum Mach bounds for two held exits."""
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def signed_work_control(shape):
    """Partition a common, smooth azimuthal energy between two contiguous materials.

    n=s²/(s²+a²) is a manufactured healing profile. For unit winding the
    angular energy density n/(2s²) has no singularity. This field is held fixed
    in physical space. It need not satisfy a reflecting wall's density condition.
    """
    nodes, weights = leggauss(96)
    t=(nodes+1)/2
    b,slope,_=shape.profile(t)
    # Differentiate with respect to the dimensional terminal radius; db(t)=t db.
    weight=weights*shape.length/2
    for a in (.4,1.):
        def polar_integral(radius):
            # Integral over the physical B³ using cylindrical s and both axial hemispheres.
            return 2*np.pi*quad(lambda x:
                radius**3*np.sin(x)*np.cos(x)**2/(a*a+(radius*np.sin(x))**2),
                0,np.pi/2,epsabs=1.e-12,epsrel=2.e-13)[0]

        def surface_integral(radius):
            # Surface stress by Reynolds transport. This integrates the actual S².
            cos=nodes
            s2=radius**2*(1-cos**2)
            return np.pi*radius**2*(weights@(1/(s2+a*a)))

        direct=weight@(t*np.array([surface_integral(r) for r in b]))
        outer=shape.neck*5
        for step in (.001,.0005):
            # Integrate the OUTER region independently in spherical coordinates,
            # with its fixed far boundary, rather than setting its force to -F_inner.
            def outside(lo):
                r=lo+(outer-lo)*(nodes+1)/2
                s2=r[:,None]**2*(1-nodes[None,:]**2)
                sphere=np.pi*r*r*((1/(s2+a*a))@weights)
                return .5*(outer-lo)*(weights@sphere)
            E0=[]; E1=[]
            for sign in (-1,1):
                local=b+sign*step*t
                E0.append(weight@np.array([polar_integral(q) for q in local]))
                E1.append(weight@np.array([outside(q) for q in local]))
            Fin=(E0[0]-E0[1])/(2*step)
            Fout=(E1[0]-E1[1])/(2*step)
            emit(kind='common_field_partition',radius=shape.neck,core_width_for_control=a,
                 signed_step=step,inside_force=Fin,outside_force=Fout,net=Fin+Fout,
                 independent_inner_surface_force=-direct)
            if max(abs(Fin+Fout),abs(Fin+direct))>3.e-6*abs(direct):
                raise ArithmeticError('Partition and surface work disagree')


def run():
    record=json.loads((ROOT/'research/pde_ledger_v3/_measurements/mouth_breathing_exit_support.json').read_bytes())
    for radius in (1,2):
        trial=Trial(radius,72,bonded_core=True)
        path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz'
        shifts,_,_,sha,status,*_=inputs(path,radius,72,trial)
        if sha!=record['physical_preparations'][f'r{radius}n72'] or status!='balanced':
            raise ValueError('Reference mismatch')
        shape=Meridian(radius,trial.x+shifts[:4]); b=shape.exit
        run=record['runs']['mouth_breathing_matched01/'+('r1n72_direction' if radius==1 else 'r2n72')]
        raw=run['literal_stdout'].encode()
        assert hashlib.sha256(raw).hexdigest()==run['stdout_sha256'] and run['guard_exit']==run['child_exit']==0
        lines=[row for row in map(json.loads,raw.splitlines())
               if row.get('index')==0 and row.get('exit_step')==.006 and 'current_closing_load' in row]
        assert len(lines)==1
        for degree in (1,2):
            cyl=min(degree/b,1.)
            vol=1-(1-cyl**2)**1.5
            rim=1-np.sqrt(1-cyl**2)
            direct=3*quad(lambda x:x*np.sqrt(1-x*x),0,cyl,epsabs=1.e-13)[0]
            emit(kind='kinematic_bound',radius=radius,endpoint_sha256=sha,exit_radius=b,
                 closing_load=lines[0]['current_closing_load'],degree=degree,
                 circulation=2*np.pi*degree,minimum_total_mach_lower_bound=degree/b,
                 volume_fraction_certified=vol,rim_area_fraction_certified=rim,
                 cap_volume_crosscheck=abs(vol-direct))
            if abs(vol-direct)>1.e-11:
                raise ArithmeticError('Independent cylindrical volume mismatch')
        signed_work_control(shape)


if __name__=='__main__':
    run()
