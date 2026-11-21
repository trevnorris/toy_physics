"""Conditional 0PN momentum transfer from the paired feeding ball into a held rim.

Reuse the isolated curved carrier to clock material packets, and the separate
flat-face paired preparations to supply the incident tangential carriage.
Cumulative isolated intake registers shells between the two charts. This
matches transported rate labels only. All fresh packets start unstretched.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator

from incident_flow_falloff import check_current, checked_grid
from lateral_source_rate import verified_packet
from lateral_target_rate import force
from mouth_contact_age import PRIOR, emit
from mouth_feed_weak_check import dx4
from mouth_free_contact import Carrier
from mouth_free_order_diffuse import read
from mouth_front_pressure_dipole import MANIFEST, grid


PATCH = Path('research/pde_ledger_v3/_measurements/mouth_neighbor_patch.json')
POLICIES = ('intact','reconnect','increased_drag')
COLUMNS = ('mouth','brane','out','tow_mouth','tow_brane','tow_out')


def constitutive_data():
    record = json.loads(PATCH.read_text())['runs']['mouth_neighbor_patch48_ramp_01']
    if record['exit_code'] or not record['limits_verified']:
        raise ValueError('Invalid patch calibration')
    row = json.loads(record['stdout_lines'][0])
    return np.array(row['spring']), np.array(row['drag']), row['contact_spring']


def transit(carrier, order, eta, section):
    """Start at an upstream circular section, then follow one r,w streamline."""
    dom = carrier.domain
    def quantities(w, r):
        psi,pr,pw = read(dom,carrier.field,r,w)
        return (pr/psi).imag, (pw/psi).imag, abs(psi)**2
    def rhs(w, y):
        ur,uw,_ = quantities(w,y[0])
        if not (0<y[0]<dom.radius(w)[0] and uw>0):
            raise ValueError('Streamline left inward chart')
        return [ur/uw,1/uw]
    def level(value):
        def cross(w,y):
            return read(dom,order,y[0],w)[0]-value
        cross.direction=-1
        cross.terminal=value==.25
        return cross
    values=(.75,.60,.40,.25)
    sol=solve_ivp(rhs,(section,section+.45),[eta*dom.radius(section)[0],0.],
                  events=[level(v) for v in values],dense_output=True,
                  rtol=1e-9,atol=1e-11,max_step=.003)
    if not sol.success or any(len(v)!=1 for v in sol.t_events):
        raise ValueError('Missing an order crossing')
    times=np.array([a[0,1] for a in sol.y_events])
    times-=times[0]
    w=np.linspace(sol.t_events[0][0],sol.t_events[-1][0],129)
    r,t=sol.sol(w)
    t-=sol.y_events[0][0,1]
    n=np.array([quantities(wi,ri)[2] for wi,ri in zip(w,r)])
    q=np.array([read(dom,order,ri,wi)[0] for ri,wi in zip(r,w)])
    if np.any(np.diff(q)>=0) or np.min(n)<=0:
        raise ValueError('Nonmonotone passage or nonpositive material density')
    # Per unit entry volume, with unit mass per constituent. Spring and
    # viscous coefficients per current volume are carried into that volume.
    return times, n[0], PchipInterpolator(t,np.c_[q/.25*n[0]/n,n[0]/n],axis=0)


def throughput_coordinate(carrier, eta, w):
    # Include every radial element boundary of the Q2 chart.
    vg=carrier.domain.blocks[2][1]
    cuts=np.unique(np.r_[0.,(1-vg)[(1-vg>0)&(1-vg<eta)],eta])
    nodes,weights=leggauss(6)
    rhat=((cuts[1:]+cuts[:-1])[:,None]+np.diff(cuts)[:,None]*nodes)/2
    a,_,psi,_,uw=carrier.eval(rhat.ravel(),np.full(rhat.size,w))
    flux=(4*np.pi*a**3*rhat.ravel()**2*abs(psi)**2*uw).reshape(rhat.shape)
    if np.min(uw)<=0:
        raise ValueError('Rate label is not monotone')
    return .5*np.diff(cuts)@(flux@weights)


def packet(law, path, policy, U, V, W=0.):
    """Lab impulses per entry mass. Q tracks viscous work; S tracks gate work."""
    k,d,kc=law
    clocks,m,scale=path
    ds=scale.derivative()
    # v, three neighbor extensions, contact extension, mouth/brane impulses,
    # dashpot dissipation, work from evolving elastic coefficients.
    y=np.zeros(9)
    y[0]=U
    rupture=0.
    intervals=[(0.,clocks[-1],True)] if policy=='intact' else [
        (clocks[0],clocks[1],True),(clocks[1],clocks[2],False),(clocks[2],clocks[3],True)]
    for t0,t1,on in intervals:
        if not on:
            rupture+=.5*scale(t0)[0]*kc*y[4]**2
            y[4]=0.
        drags=d.copy()
        if not on and policy=='increased_drag':
            drags[0]*=4
        def rhs(t,z):
            s,vol=scale(t)
            v,e,link=z[0],z[1:4],z[4]
            visc=vol*drags
            far=(v-W+np.sum(s*k*e/visc))/np.sum(1/visc)
            de=(far-s*k*e)/visc
            rim=s*kc*link if on else 0.
            return np.r_[-(rim+far)/m,de,(v-V if on else 0.),rim,far,
                         visc@(de*de),.5*ds(t)[0]*(k@(e*e)+on*kc*link**2)]
        sol=solve_ivp(rhs,(t0,t1),y,rtol=1e-9,atol=3e-12,method='DOP853',
                      max_step=clocks[-1]/192)
        if not sol.success:
            raise ArithmeticError(sol.message)
        y=sol.y[:,-1].copy()
    v,e,link=y[0],y[1:4],y[4]
    rim,far,heat,work=y[5:]
    remaining=.5*scale(clocks[-1])[0]*(k@(e*e)+kc*link**2)
    dp=m*(U-v)-rim-far
    de=.5*m*(U**2-v**2)+work-rim*V-far*W-heat-rupture-remaining
    if abs(dp)>2e-8*m or abs(de)>4e-7*m or min(heat,rupture,remaining)<-1e-10:
        raise ArithmeticError(f'Packet balance failed: {dp/m=}, {de/m=}, {heat=}, {rupture=}, {remaining=}')
    return np.array([rim/m,far/m,v]), np.array([dp/m,de/m,work/m,heat/m,rupture/m,remaining/m])


def transfer(mesh, state, count):
    records=map(json.loads,json.loads(PRIOR.read_text())['runs']['mouth_mixed_order_ripple02']['stdout_lines'])
    p,=[a for a in records if a['test']=='archive_patch' and (a['mesh'],a['profile'],a['state'])==(mesh,'cubic',state)]
    carrier=Carrier(Path(p['carrier_file']).parent)
    path=Path(p['order_file'])
    if sha256(path.read_bytes()).hexdigest()!=p['order_sha256']:
        raise ValueError('Order pin mismatch')
    with np.load(path) as f:
        order=f['chi']
    section=.501 if state=='near' else 1.15
    total=throughput_coordinate(carrier,1.,section)
    if abs(total-.25)>.0005:
        raise ValueError('Flared-section intake mismatch')
    law=constitutive_data()
    labels=[]
    table={key:[] for key in POLICIES}
    for i,eta in enumerate(np.linspace(.005,.995,count)):
        coord=throughput_coordinate(carrier,eta,section)/total
        path=transit(carrier,order,eta,section)
        labels.append(coord)
        for key in POLICIES:
            held,hbudget=packet(law,path,key,1.,0.)
            tow,tbudget=packet(law,path,key,0.,1.)
            coeff=np.r_[held,tow]
            if i==count//2:
                mixed,_=packet(law,path,key,-.5,.25,.125)
                predicted=-.625*held+.125*tow+[0.,0.,.125]
                if max(abs(mixed-predicted))>2e-7:
                    raise ArithmeticError('Independent mixed-drive or boost check failed')
            table[key].append(coeff)
            emit(test='ray',state=state,mesh=mesh,eta=eta,rate_label=coord,policy=key,
                 clocks=path[0].tolist(),entry_mass=path[1],responses=coeff.tolist(),
                 held_budget=hbudget.tolist(),tow_budget=tbudget.tolist(),section_rate=total)
    if min(np.diff(labels))<=0 or labels[0]>.002 or labels[-1]<.975:
        raise ValueError('Insufficient coverage of cumulative flow')
    return np.array(labels), {k:np.array(a) for k,a in table.items()}


def feeding_packet(path,pins,g,centers):
    with verified_packet(path,pins) as f:
        checked_grid(g,f)
        if 'centers' in f:
            np.testing.assert_array_equal(f['centers'],centers)
        elif int(f['d'])!=-centers[0]:
            raise ValueError('Feeding preparation mismatch')
        currents=f['currents']
        if len(currents)!=len(centers):
            raise ValueError('Wrong mouth count')
        for c,j in zip(centers,currents):
            check_current(g,c,j[None])
        layers=f['layers'] if 'layers' in f else f['field'][-4:]
    j=currents[-1]
    carried,stress,direct=force(g,layers,j,.25)
    face=layers[-1]
    n=abs(face)**2
    ux=(face.conj()*dx4(face,g.dx)).imag/n
    even=np.concatenate([face[:,1::-1],face,face[:,-1:-3:-1]],axis=1)
    radial=(even[:,:-4]-8*even[:,1:-3]+8*even[:,3:-1]-even[:,4:])/(12*g.dx)
    ur=(face.conj()*radial).imag/n
    ix=abs(g.x)<1+2*g.dx
    ir=g.r<1+2*g.dx
    data=np.stack([j,j*ux,j*ur],axis=-1)[ix][:,ir]
    return data, np.array([carried,stress,direct])


def samples(g,data,method):
    x=g.x[abs(g.x)<1+2*g.dx]; r=g.r[g.r<1+2*g.dx]
    if method=='cell':
        xx,rr=np.meshgrid(x,r,indexing='ij')
        radius=np.hypot(xx,rr)
        inside=radius<1.
        return radius[inside],(xx/radius)[inside],(2*np.pi*rr*g.dx**2)[inside],data[inside]
    # Integrate the cylindrical interpolant in independent spherical coordinates.
    # The radial vector component is odd and scalar/axial components even at the axis.
    imag=data[:,:1,:].copy()
    imag[:,:,2]*=-1
    interp=RegularGridInterpolator((x,np.r_[-r[0],r]),np.concatenate([imag,data],axis=1))
    a,b=leggauss(128)
    radius=(a+1)/2
    s,mu=np.meshgrid(radius,a,indexing='ij')
    weights=np.pi*s**2*b[:,None]*b[None,:]
    points=np.c_[(s*mu).ravel(),(s*np.sqrt(1-mu**2)).ravel()]
    return s.ravel(),mu.ravel(),weights.ravel(),interp(points)


def paired_case(record, labels, table, mesh, state, count):
    if record['guard_exit'] or record['worker_exit'] or not record['limits_verified']:
        raise ValueError('Feeding archive guard failed')
    ref=json.loads(record['stdout'][-1]); d=ref['d']
    g=grid(*ref['grid'])
    pins={Path(k):v for k,v in (record['input_sha256']|record['output_fields_sha256']).items()}
    args=shlex.split(record['command'])
    src=[Path(args[args.index(k)+1]) for k in ('--baseline','--pair' if '--pair' in args else '--output')]
    (alone,A),(pair,P)=[feeding_packet(p,pins,g,c) for p,c in zip(src,([0],[-d,0]))]
    flux=P-A
    np.testing.assert_allclose(flux,[ref['carried'],ref['stress'],ref['direct']],rtol=5e-12,atol=5e-15)
    for method in ('cell','spherical'):
        s,mu,volume,base=samples(g,alone,method)
        diff=samples(g,pair,method)[-1]-base
        j=volume*base[:,0]
        _,groups=np.unique(s.round(12),return_inverse=True)
        rates=np.bincount(groups,weights=j)
        rank=((rates.cumsum()-.5*rates)/rates.sum())[groups]
        # Axis x points away from companion. Both cylindrical components of
        # the actual momentum vector enter the in-shell (angular) projection.
        tangent=volume*(-(1-mu**2)*diff[:,1]+mu*np.sqrt(1-mu**2)*diff[:,2])
        axial=-volume*diff[:,1]
        radial=axial-tangent
        for key,coeff in table.items():
            c=np.stack([np.interp(rank,labels,col) for col in coeff.T])
            mouth,brane,exit_t=c[:3]@tangent
            bypass=radial.sum()+flux[1]
            out=exit_t+bypass
            projected_rate=j*(1-mu**2)
            pull_slope=c[0]@projected_rate
            tow_drag=-c[3]@projected_rate
            if abs(mouth+brane+out-axial.sum()-flux[1])>1e-10 or (method=='cell' and abs(axial.sum()-flux[0])>3e-14):
                raise ArithmeticError('Opening momentum ledger mismatch')
            emit(test='opening',state=state,clock_mesh=mesh,rays=count,feed_grid=ref['grid'],d=d,
                 policy=key,method=method,input_pins={str(p):pins[p] for p in src},
                 integrated_mass=j.sum(),tangent=tangent.sum(),radial=radial.sum(),
                 feed_stress_passed=flux[1],carrier_total=axial.sum()+flux[1],
                 archived_direct=flux[2],mouth=mouth,brane=brane,export=out,
                 export_tangent=exit_t,uniform_held_slope=pull_slope,towing_drag_slope=tow_drag,
                 source_incident_speed=ref['incident_speed'],uniform_predicts=pull_slope*ref['incident_speed'],
                 force_balance_error=mouth+brane+out-axial.sum()-flux[1])
    return


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(34,48),required=True)
    parser.add_argument('--state',choices=('near','deep'),required=True)
    parser.add_argument('--rays',type=int,required=True)
    args=parser.parse_args()
    labels,table=transfer(args.mesh,args.state,args.rays)
    for record in json.loads(MANIFEST.read_text())['runs']:
        if record['name']!='compare':
            paired_case(record,labels,table,args.mesh,args.state,args.rays)


if __name__=='__main__':
    main()
