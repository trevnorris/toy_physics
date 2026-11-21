"""Independent area-energy estimate of axial force transmission by a tense wall.

Treat the existing weak-intake two-coordinate inertia and held loads as data.
Translate each circular cross section without changing radius or base shape.
Test the quadratic area stiffness by directly integrating both tilted surfaces.
"""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_inertia import taper
from mouth_flared_receiver import FlaredDomain, emit


def area_stiffness(dom,switch):
    pg=dom.blocks[2][0]
    lo,hi=switch
    edges=np.r_[lo,pg[(pg>lo)&(pg<hi)],hi]
    nodes,weights=leggauss(12)
    half=np.diff(edges)[:,None]/2
    w=(edges[:-1,None]+half*(1+nodes)).ravel()
    weight=(half*weights).ravel()
    r,slope=dom.radius(w)
    _,deriv=taper(w,lo,hi)
    k=4*np.pi/3*np.dot(weight,r*r*deriv**2/(1+slope**2)**1.5)
    def area(d):
        metric=np.sqrt(1+(slope[:,None]+d*deriv[:,None]*nodes)**2)
        return 2*np.pi*np.dot(weight,r*r*(metric@weights))
    estimates=[(area(d)+area(-d)-2*area(0.))/d**2 for d in (.002,.001)]
    if abs(estimates[-1]/k-1)>1e-5 or abs(estimates[-1]-k)>abs(estimates[0]-k)*.32:
        raise ArithmeticError('Area stiffness did not reproduce the direct deformations')
    return k,estimates


if __name__=='__main__':
    reference=Path('_scratch/s11c/carrier_attachment_40_01')
    if json.loads((reference/'outcome.json').read_text())['exitCode'] or not json.loads((reference/'limit-validation.json').read_text())['verified']:
        raise ValueError('Invalid reference')
    # Only the common radius-one wall profile enters the area integral.
    dom=FlaredDomain(20,'cubic',15.,far=14.)
    for line in (reference/'stdout').read_text().splitlines():
        row=json.loads(line)
        switch=row['switch']
        k,fd=area_stiffness(dom,switch)
        mass=np.array(row['mass_matrix'])
        held=np.array(row['held_per_intake'])
        total=mass.sum()
        h=mass.sum(axis=0)/total
        mu=np.linalg.det(mass)/total
        fraction=abs(row['connector_on_distal_per_held_force'])
        direction=np.array([-1.,1.])
        for omega in (.02,.01):
            # A 5% differential excursion is an illustrative resolution criterion.
            tolerance=.05
            threshold=omega**2*(mu+total*fraction/tolerance)/k
            readings=[]
            for tension in (.5*threshold,threshold,2*threshold):
                x=np.linalg.solve(tension*k*np.outer(direction,direction)-omega**2*mass,held)
                relative=abs(direction@x/(h@x))
                readings.append(relative)
            emit(test='tension',switch=switch,omega=omega,mass=total,relative_mass=mu,
                 stiffness_per_unit_tension=k,area_second_differences=fd,
                 fractional_excursion=tolerance,threshold_tension=threshold,
                 tension_multipliers=[.5,1,2],relative_excursions=readings)
            np.testing.assert_allclose(readings[1],tolerance,rtol=5e-11)
