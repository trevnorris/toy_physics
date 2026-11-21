"""Material order-energy crossing and stress work on cuts of completed 4D order fields."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss

from mouth_bend_drive_control import line_rule
from mouth_carrier_lateral_response import trace
from mouth_empty_bend_stiffness import shapes
from mouth_flared_receiver import FlaredDomain, emit
from mouth_order_wave_pulse import digest


def boundary(dom, fields, offset, meta, block, fixed_p, fixed, normal, interval=None):
    pg,vg,_=dom.blocks[block]
    edges=vg if fixed_p else pg
    if interval:
        lo,hi=interval
        edges=np.r_[lo,edges[(edges>lo)&(edges<hi)],hi]
    q,weights=line_rule(edges)
    p,v=(np.full_like(q,fixed),q) if fixed_p else (q,np.full_like(q,fixed))
    r,w,rp,rv,wp,wv,f,fr,fw=trace(dom,block,p,v,fields)
    tangent_r,tangent_w=(rv,wv) if fixed_p else (rp,wp)
    nr,nw=tangent_w,-tangent_r
    sign=np.sign(nr*normal[0]+nw*normal[1])
    if np.any(sign==0):
        raise ValueError('Surface orientation ambiguous')
    nr,nw=nr*sign,nw*sign
    # All chosen surfaces lie within the uncut bend profile or outside it.
    _,sr,sw=shapes(r,w,'compact',dom.far,dom.roof)
    br,bw=sr@offset,sw@offset
    if np.any(abs(br)>1e-13):
        raise ValueError('Unexpected shape taper on a measured surface')
    modes=f.shape[-1]//2
    angular,ang_wt=leggauss(max(12,4*modes))
    laws=[Legendre.basis(l) for l in range(modes)]
    result=np.zeros(4)
    kappa=meta['nu']/meta['mobility']
    for cos,aw in zip(angular,ang_wt/2,strict=True):
        y=np.array([law(cos) for law in laws])
        dy=np.array([law.deriv()(cos) for law in laws])
        def evaluate(s):
            phi=f[:,s]@y
            dr=fr[:,s]@y
            da=(f[:,s]@dy)/r
            return phi,dr,da,fw[:,s]@y-bw*(cos*dr+(1-cos*cos)*da)
        psi,pr,pa,pw=evaluate(slice(0,modes))
        c,cr,ca,cw=(x.real for x in evaluate(slice(modes,None)))
        ur,ua,uw=(np.imag(x/psi) for x in (pr,pa,pw))
        adv=ur*cr+(1-cos*cos)*ua*ca+uw*cw
        e=kappa*(cr*cr+(1-cos*cos)*ca*ca+cw*cw+(c*(1-c)/meta['delta'])**2)/2
        normal_w=nw-nr*bw*cos
        vn=nr*ur+normal_w*uw
        gn=nr*cr+normal_w*cw
        area=4*np.pi*aw*weights*r*r
        transport=area@(e*vn)
        stress=area@(-kappa*adv*gn)
        result+=[transport,stress,transport+stress,area@(abs(psi)**2*vn)]
    return dict(zip(('advected_order_energy','gradient_stress_work','sum','number'),result,strict=True))


def run(folder):
    def read(name):
        return json.loads((folder/name).read_text())
    status,limits=read('outcome.json'),read('limit-validation.json')
    if (status['exitCode'] or status['childOutcome']['exitCode'] or not status['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError('Source field run has not passed guard verification')
    lines=[json.loads(q) for q in (folder/'stdout').read_text().splitlines()]
    prep,setup=lines[:2]
    dom=FlaredDomain(setup['mesh'],prep['profile'],prep['roof'],prep['far'])
    for info in (row for row in lines if row['test']=='field'):
        path=Path(info['path'])
        if digest(path)!=info['sha256']:
            raise ValueError('Field readback mismatch')
        with np.load(path) as field:
            data=np.column_stack([field['psi'],field['chi']])
        offset=np.array([0.,info['bend']])
        def cut(block,along_p,fixed,normal,interval=None):
            return boundary(dom,data,offset,setup,block,along_p,fixed,normal,interval)
        interior={}
        for lo,hi in ((1.173,2.067),(1.36,1.52)):
            interior[str((lo,hi))]=dict(enter=cut(2,True,lo,(0,1)),leave=cut(2,True,hi,(0,1)),
                                         wall=cut(2,False,0.,(1,0),(lo,hi)))
        # The bath normals here all point out of the computational domain.
        ends=dict(feed=cut(1,True,dom.far,(1,0)),bulk_side=cut(4,False,1.,(1,0)),
                  roof_inner=cut(3,True,dom.roof,(0,1)),roof_outer=cut(4,True,dom.roof,(0,1)))
        emit(test='order_energy_cuts',mesh=setup['mesh'],bend=info['bend'],source_sha256=info['sha256'],
             segments=interior,baths=ends)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folders',nargs='+',type=Path)
    for folder in p.parse_args().folders:
        run(folder)
