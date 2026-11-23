"""Assemble the axial configurational residual of each FE node on the held R4 lip.

Compare the reactions on the actual boundary nodes with held-out, signed wall
work. Interior residuals are measured before reading any displacement benchmark.
Both banks meet at one shared meridional vertex at the disk; report its separate
contribution once, including the cells on the receiving side of the disk.
"""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_disk_pillbox import dot, physical, stresses
from mouth_last_ring_solid import interpolation
from mouth_rim_disk_work_review import Disk
from mouth_rim_shared_crystal import TerminalQuadrature
from mouth_supercore_joined_front import file_sha

SECTORS=('amplitude','phase','order','vortex_volume','crystal','vortex_wall')


def weights(at):
    return np.where(np.array(list(product((0,1),repeat=3))),at,1-at)


def volume(disk,slab,bank,radial,depth,order):
    m=disk.model
    xs,ang,zs=m.rr[bank],m.a,m.z[slab]
    # The displacement is axisymmetric. Sum the local angle contributions
    # before assembling the four vertices of each meridional element.
    ii,aa,jj=np.meshgrid(radial,np.arange(len(ang)-1),depth,indexing='ij')
    ix,ia,iz=ii.ravel(),aa.ravel(),jj.ravel()
    bits=np.array(list(product((0,1),repeat=3)))
    corner=[f.values[ix[:,None]+bits[:,0],ia[:,None]+bits[:,1],iz[:,None]+bits[:,2]]
            for f in disk.reader.fields[slab,bank]]
    dx,da,dz=np.diff(xs)[ix],np.diff(ang)[ia],np.diff(zs)[iz]
    nodes=np.zeros((len(xs),len(zs),len(SECTORS)))
    q,wg=np.polynomial.legendre.leggauss(order)
    q=(q+1)/2;wg=wg/2
    for i,j,k in product(range(order),repeat=3):
        f=weights(np.array([q[i],q[j],q[k]]))
        basis=f.prod(axis=1)
        dN=np.stack([(2*bits[:,dim]-1)*f[:,(dim+1)%3]*f[:,(dim+2)%3] for dim in range(3)])
        x=xs[ix]+q[i]*dx
        a=ang[ia]+q[j]*da
        z=zs[iz]+q[k]*dz
        R,s=disk.sh.radius(z)
        width=28.-R if bank else R
        lever=1-x if bank else x
        r=R+width*x if bank else R*x
        wt=4*np.pi*r*r*np.sin(a)*width*dx*da*dz*wg[i]*wg[j]*wg[k]
        F,Fp=disk.reader.vortex(r*np.sin(a))
        tr=np.zeros((4,len(x)));tz=np.zeros_like(tr)
        for sign,n in zip((1.,-1.),corner):
            val=np.einsum('i,qij->qj',basis,n)
            grad=np.einsum('ki,qij->qjk',dN,n)
            grad[:,:,0]/=(width*dx)[:,None]
            grad[:,:,1]/=(r*da)[:,None]
            grad[:,:,2]=grad[:,:,2]/dz[:,None]-(s*lever)[:,None]*grad[:,:,0]
            state=physical(1,r,a,val.T,grad.transpose(1,0,2),F,Fp)
            Tr,Tz=stresses(state)
            tr[:3]+=sign*Tr;tz[:3]+=sign*Tz
            gf,_,_,density,swirl,_=state
            gv=grad[:,0]
            y=val[:,0]**2
            u=F**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20
            delta=.5*(F**2*dot(gv,gv)-dot(gf,gf)-density*swirl)+u-(density**5-5*density+4)/20
            tr[3]+=sign*(F**2*gv[:,0]*gv[:,2]-gf[:,0]*gf[:,2])
            tz[3]+=sign*(F**2*gv[:,2]**2-gf[:,2]**2-delta)
        for di,dk in product((0,1),repeat=2):
            sx=(q[i] if di else 1-q[i]);sz=(q[k] if dk else 1-q[k])
            gr=(2*di-1)*sz/(dx*width)
            gz=(2*dk-1)*sx/dz-s*lever*gr
            force=((tr*gr+tz*gz)*wt).T.reshape(*ii.shape,4).sum(axis=1)
            nodes[np.ix_(radial+di,depth+dk,np.arange(4))]+=force
    return nodes


def add_wall(disk,bank,depth,order,nodes):
    m=disk.model
    q,w=np.polynomial.legendre.leggauss(order);q=(q+1)/2;w=w/2
    zz=m.z[1]
    a,az=np.meshgrid(np.arange(len(m.a)-1),depth,indexing='ij')
    da=np.diff(m.a)[a]
    trace=m.rr[bank][0 if bank else -1]
    for i,j in product(range(order),repeat=2):
        angle=m.a[a]+q[i]*da
        z=zz[az]+q[j]*np.diff(zz)[az]
        r=disk.sh.radius(z)[0]
        F,Fp=disk.reader.vortex(r*np.sin(angle))
        point=np.column_stack([np.full(z.size,trace),angle.ravel(),z.ravel()])
        amps=[f(point)[:,0].reshape(z.shape) for f in disk.reader.fields[1,bank]]
        weight=(1-2*bank)*.5*4*np.pi*r*r*np.sin(angle)**2*F*Fp*da*w[i]*w[j]
        # The oriented surface energy contributes minus its axial variation.
        # Integration dz cancels the 1/dz of the two endpoint derivatives.
        work=(weight*(amps[0]**2-amps[1]**2)).sum(axis=0)
        face=0 if bank else -1
        nodes[face,depth,5]+=work
        nodes[face,depth+1,5]-=work


def add_solid(disk,state,wall,fields,nodes):
    m=disk.model
    with np.load(fields) as f:
        diff=(f['wound'][2]-f['rest'][2])[m.ids[1,0]]
    delta=np.einsum('raj,a->jr',diff,angular_weights(m.a))
    quad=TerminalQuadrature(state,wall,m.rr[0],m.z[1],64)
    xs,zs=quad.ax,quad.z
    for v in quad.blocks():
        x,z=v['x'],v['z']
        ii=np.clip(xs.searchsorted(x)-1,0,len(xs)-2)
        jj=np.clip(zs.searchsorted(z)-1,0,len(zs)-2)
        dx=np.diff(xs)[ii];dz=np.diff(zs)[jj]
        p=(x-xs[ii])/dx;t=(z-zs[jj])/dz
        chi=interpolation(zs,xs,delta,z,x)[0]
        R,s=disk.sh.radius(z)
        K,J=v['K'],v['J']
        rr=-chi*v['measure']*K[:,0,1]
        rz=-chi*v['measure']*(K[:,1,1]-1+.5*J*(J-1))
        for di,dk in product((0,1),repeat=2):
            gr=(2*di-1)*(t if dk else 1-t)/(dx*R)
            gz=(2*dk-1)*(p if di else 1-p)/dz-x*s*gr
            np.add.at(nodes[:,:,4],(ii+di,jj+dk),rr*gr+rz*gz)


def run(a):
    disk=Disk(a.state,a.wall,a.fields)
    m=disk.model
    with np.load(a.reference) as f:
        grid0=[f['inside_x'],f['outside_x']]
        zd=f['duct_z']
    # The exact original native taper nodes occur in both archived meshes.
    # Reconstruct the native exterior axis rather than infer a cylindrical radius.
    from mouth_exchange_carrier import mapped_nodes
    native_x=mapped_nodes(13,disk.sh.length,28.,22.)[0][1]
    mid=np.argmin(abs(zd-(disk.sh.length-m.h)))
    strips=dict(near=zd[[-9,-5,-1]],middle=zd[[-25,-17,-9]],upstream=zd[[mid-4,mid-3,mid-1]])
    regions={}
    for slab,bank in product((1,2),(0,1)):
        z=m.z[slab]
        lo=min(strips['upstream'][0],disk.windows[-1][0])
        hi=disk.windows[-1][-1]
        depth=np.flatnonzero((z[1:]>lo+1e-13)&(z[:-1]<hi-1e-13))
        radial=np.arange(len(m.rr[bank])-1) if bank==0 else np.flatnonzero(m.rr[bank][:-1]<native_x[6]-1e-13)
        nodes=volume(disk,slab,bank,radial,depth,a.quad)
        if slab==1:add_wall(disk,bank,depth,a.quad,nodes)
        regions[slab,bank]=nodes
    add_solid(disk,a.state,a.wall,a.fields,regions[1,0])
    # No input from the measured wall work or pointwise slit tractions entered
    # in the above calculation. Extract physical bank test functionals first.
    probes=[]
    for bank in (0,1):
        face=0 if bank else -1
        x,z=m.rr[bank],m.z[1]
        nodes=regions[1,bank]
        for strip,(lo,peak,hi) in strips.items():
            tent=np.maximum(0,np.minimum((z-lo)/(peak-lo),(hi-z)/(hi-peak)))
            boundary=np.einsum('z,zs->s',tent,nodes[face])
            readings=[]
            for size in (4,8):
                edge=grid0[bank][size if bank else -1-size]
                taper=np.clip(1-x/edge if bank else (x-edge)/(1-edge),0,1)
                total=np.einsum('x,z,xzs->s',taper,tent,nodes)
                readings.append(dict(width=size,total=total.tolist(),interior=(total-boundary).tolist()))
            probes.append(dict(bank=('duct','outside')[bank],strip=strip,face=boundary.tolist(),readings=readings))
    enclosures=[]
    for k,win in enumerate(disk.windows,2):
        l0,l1,h0,h1=win
        for j in (4,5):
            body=np.zeros(len(SECTORS));slit=np.zeros_like(body);joint=np.zeros_like(body)
            incident_tip=np.zeros_like(body);contrib=[]
            for (slab,bank),nodes in regions.items():
                z=m.z[slab]
                az=np.clip((z-l0)/(l1-l0),0,1)*np.clip((h1-z)/(h1-h0),0,1)
                x=m.rr[bank]
                rad=np.ones_like(x) if not bank else np.clip((native_x[j+1]-x)/(native_x[j+1]-native_x[j]),0,1)
                all_nodes=np.einsum('x,z,xzs->s',rad,az,nodes)
                face=0 if bank else -1
                tip=nodes[face,-1 if slab==1 else 0]
                bank_load=az[:-1]@nodes[face,:-1] if slab==1 else np.zeros_like(tip)
                joint+=tip;slit+=bank_load
                if slab==1:incident_tip+=tip
                body+=all_nodes-tip-bank_load
                contrib.append(dict(slab=slab,side=bank,tip=tip.tolist(),bank_without_tip=bank_load.tolist()))
            enclosures.append(dict(axial_band=k,lateral_band=j,bulk_nodes=body.tolist(),
                                   component_rows=contrib,slit_separate=slit.tolist(),shared_tip=joint.tolist(),
                                   incident_only_tip=incident_tip.tolist(),native_translation_reaction=(slit+joint+body).sum(),
                                   recovered_banks_with_shared_tip=(slit+joint).sum()))
    # Now open the held-out control records and compare the predictions.
    ledger=json.loads(a.bank_measure.read_text())['data'][disk.level]
    if ledger['source']['field']!=disk.source_hash['field']:raise ValueError('Benchmark on different fields')
    for row in probes:
        for r in row['readings']:
            bench=next(b for b in ledger['results'] if b['strip']==row['strip'] and b['bank']==row['bank']
                       and b['radial_cells']==r['width'] and b['quadrature']==a.quad)
            actual=bench['total_native']
            r.update(signed_energy=actual,pointwise=bench['outward_wall_axial'],
                     boundary_minus_signed=sum(row['face'])-actual,volume_check=sum(r['total'])-actual)
            if abs(r['volume_check'])>3e-5:raise ArithmeticError('Independent signed-work validation failed')
    closed=[json.loads(s) for s in (a.fields.parent/'stdout').read_text().splitlines() if s.startswith('{')]
    for v in enclosures:
        old=next(b for b in closed if b['kind']=='closed_work' and b['axial_band']==v['axial_band']
                 and b['lateral_band']==v['lateral_band'])
        ext=old['outward_momentum']
        v.update(external=ext,strong_closed=old['closed_plus_slit'],
                 conditional_closed=ext+v['recovered_banks_with_shared_tip'],
                 closure_check=ext+v['native_translation_reaction'])
        if abs(v['closure_check'])>(2e-6 if a.quad==2 else .02):raise ArithmeticError('Common-translation audit')
    return dict(level=disk.level,source=disk.source_hash,gauss_per_axis=a.quad,sectors=SECTORS,
                sign='Reaction = minus the driven-minus-resting signed energy derivative. Integrated nodal moments use unit displacement and include crystal and vortex surface terms.',
                extent='Full interior radius; exterior through both fixed taper radii. Same physical heights and tapers across meshes.',
                probes=probes,enclosures=enclosures)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','fields','reference','bank_measure','output'):p.add_argument(key,type=Path)
    p.add_argument('--quad',type=int,choices=(2,4),default=2)
    a=p.parse_args()
    a.output.write_text(json.dumps(run(a),indent=2,allow_nan=False)+'\n')
