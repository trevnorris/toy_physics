"""Compare layer slip and wave curl using four admissible fixed-gap motions.

Both planes and their normal can rotate together in this additional test.
All cases preserve the gap exactly. The site energy is evaluated from the
rotated physical bonds, using the same radial and layered microscopic laws.
"""
import argparse
import numpy as np
from neutral_patch_ringing_modulus import Crystal, DIRECTIONS, collar_orders
from neutral_patch_sites import BETWEEN, GAP, PITCH_UNIT, WITHIN, emit, morse


def physical_energy(bonds, cross, law, normal):
    z = bonds@normal
    tangent = bonds-z[:,None]*normal
    q2 = np.sum(tangent*tangent,axis=1)
    if cross and law=='layered':
        q2 = q2*q2/PITCH_UNIT**2
    return morse(np.sqrt(z*z+q2), *(BETWEEN if cross else WITHIN))[0].sum()


def energy_density(crystal, origin, slip, angle, order):
    c,s = np.cos(angle),np.sin(angle)
    rotation = np.array([[c,0.,-s],[0.,1.,0.],[s,0.,c]])
    xy = crystal.xy[np.any(crystal.ij!=0,axis=1)]
    intra = np.column_stack([xy,np.zeros(len(xy))])@rotation.T
    inter = np.column_stack([crystal.xy+origin+GAP*slip*DIRECTIONS[0],
                             np.full(len(crystal.xy),GAP)])@rotation.T
    normal = rotation[:,2]
    # Two planes, each with one half of the intra-plane bond sum per cell.
    e = physical_energy(intra,False,crystal.law,normal)
    e += order/.95*physical_energy(inter,True,crystal.law,normal)
    return e/(crystal.area*GAP)


def run(law):
    order = collar_orders()[1]
    crystal = Crystal(law,8)
    origin = crystal.register()
    e0 = energy_density(crystal,origin,0.,0.,order)
    k0 = crystal.reading(origin)[2][0,0]*GAP*order/(crystal.area*.95)
    for amp in (.001,.0005,.00025):
        for name, slip, turn in [('tangent_by_depth',1,0),('depth_by_tangent',-1,-1),
                                  ('symmetric',2,1),('rotation',0,1)]:
            plus,minus = [energy_density(crystal,origin,i*amp*slip,i*amp*turn,order)
                          for i in (1,-1)]
            measured = (plus+minus-2*e0)/(amp*amp*k0)
            # Derivative at the straight reference of R(angle) (I + slip e_x ⊗ e_w).
            gradient = np.array([[0.,slip-turn],[turn,0.]])
            curl = gradient[1,0]-gradient[0,1]
            emit(test='motions',law=law,order=order,probe=amp,prestrain=0.,motion=name,
                 gradient_xw_wx=[gradient[0,1],gradient[1,0]],unit_slip_stiffness=k0,
                 site_energy_curvature_ratio=measured,unit_curl_energy_curvature_ratio=curl*curl)
            np.testing.assert_allclose(measured,slip*slip,atol=8e-5,rtol=2e-5)
    for s in (1/32,1/64,1/128):
        for direction in (-1.,1.):
            shift = origin+direction*s*GAP*DIRECTIONS[0]
            expected = crystal.reading(shift)[2]
            for h in (.0002,.0001):
                columns = []
                for v in DIRECTIONS:
                    gp = crystal.reading(shift+GAP*h*v)[1]
                    gm = crystal.reading(shift-GAP*h*v)[1]
                    columns.append((gp-gm)/(2*GAP*h))
                fd = np.column_stack(columns)
                emit(test='vector_probe',law=law,prestrain=direction*s,probe=h,
                     stress_jacobian=(fd*GAP*order/(crystal.area*.95)).tolist(),
                     max_stiffness_error_over_zero=float(abs(fd-expected).max()/crystal.reading(origin)[2][0,0]))
                np.testing.assert_allclose(fd,expected,atol=3e-5,rtol=3e-6)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('law',choices=('central','layered'))
    run(p.parse_args().law)
