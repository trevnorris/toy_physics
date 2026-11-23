"""Independent measurements on a supplied live R4 shape, lattice and three fields.

The axial test uses two physical cylinders of the same radius with different
half-heights. The volume between them contains a sloping, two-sided slit.
A separate smooth virtual translation, held constant in physical units on
every grid, measures the work between nested supports.
"""
from itertools import product

import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import legval
from scipy.interpolate import PPoly
from scipy.optimize import least_squares

from mouth_closed_rim_surface_review import Reader, quadrature
from mouth_disk_pillbox import dot, physical, stresses
from mouth_exchange_carrier import BIAS
from mouth_internal_collar_check import read
from mouth_joined_angular_work import bump
from mouth_joined_front_trace import face
from mouth_last_ring_solid import tri
from mouth_polar_shape_endpoint_review import extrema
from mouth_supercore_joined_front import ORDER_INERTIA
from mouth_supercore_shape import preparation


def geometry(shape, epsilon):
    windows = dict(pole=extrema(shape, epsilon, 1.), equator=extrema(shape, epsilon, -.5))
    at = np.array([0., .25, .34, .5, .75, 1.])*shape.length
    r = shape.radius(at)[0]
    h = bump(at, shape.length)[0]
    return dict(depth=shape.length, exit_radius=shape.exit, wall_polar_amplitude=epsilon,
                extrema=windows, depth_samples=at.tolist(), pole_radius=(r+epsilon*h).tolist(),
                equator_radius=(r-.5*epsilon*h).tolist(),
                minimum_radius=min(v['minimum'] for v in windows.values()))


def material_map(shape, epsilon, body, solid, t, x, angle):
    """Map interpolated labels (x,t,angle) to (r,z,theta) in physical space."""
    ids, bary, _, _ = tri(body.x, body.t, len(body.x), len(body.t), x, t)
    j = np.clip(body.a.searchsorted(angle)-1, 0, len(body.a)-2)
    frac = (angle-body.a[j])/np.diff(body.a)[j]
    nodes = np.asarray(solid).reshape(-1,len(body.a),3)
    lo = np.einsum('ij,ijk->ik', bary, nodes[ids,j[:,None]])
    hi = np.einsum('ij,ijk->ik', bary, nodes[ids,j[:,None]+1])
    dx, dt, dtheta = ((1-frac[:,None])*lo+frac[:,None]*hi).T
    v, z, a = x+dx, (t+dt)*shape.length, angle+dtheta
    r = v*shape.radius(z)[0]+epsilon*v**3*bump(z,shape.length)[0]*(3*np.cos(a)**2-1)/2
    return np.column_stack((r,z,a))



def orientation(shape, epsilon, body, solid):
    """Extremize the material chart Jacobian over every entire triangular prism."""
    grid=np.stack(np.meshgrid(body.t,body.x,body.a,indexing='ij'),axis=-1)[...,[1,0,2]]
    coords=grid+solid
    ranges=[[float(coords[...,k].min()),float(coords[...,k].max())] for k in range(3)]
    if not np.all(np.isfinite(coords)):
        raise ValueError('Nonfinite lattice')
    extremes=[np.inf,-np.inf]
    location={}
    for triangle,upper in enumerate((False,True)):
        corners=((1,1),(0,1),(1,0)) if upper else ((0,0),(1,0),(0,1))
        vertex=np.stack([coords[j:j+len(body.t)-1,i:i+len(body.x)-1] for i,j in corners],axis=-2)
        sign=-1. if upper else 1.
        # Columns 0,1 are constant within a triangle at any fixed angle.
        gx=sign*(vertex[...,1,:]-vertex[...,0,:])/np.diff(body.x)[None,:,None,None]
        gt=sign*(vertex[...,2,:]-vertex[...,0,:])/np.diff(body.t)[:,None,None,None]
        ga=np.diff(vertex,axis=2)/np.diff(body.a)[None,None,:,None,None]
        def det(at):
            xx=(1-at)*gx[:,:,:-1]+at*gx[:,:,1:]
            tt=(1-at)*gt[:,:,:-1]+at*gt[:,:,1:]
            return np.einsum('...m,...jm->...j',np.cross(xx,tt),ga)
        start,middle,end=(det(k) for k in (0.,.5,1.))
        aa=2*(start+end-2*middle)
        bb=end-start-aa
        critical=np.clip(np.divide(-bb,2*aa,out=np.zeros_like(aa),where=aa!=0.),0.,1.)
        for at,value in ((0.,start),(1.,end),(critical,start+critical*(bb+critical*aa))):
            for idx,find in ((0,np.argmin),(1,np.argmax)):
                key=np.unravel_index(find(value),value.shape)
                val=float(value[key])
                if (idx==0 and val<extremes[idx]) or (idx==1 and val>extremes[idx]):
                    extremes[idx]=val
                    location['min' if idx==0 else 'max']=dict(cell_t_x_angle=list(map(int,key[:3])),
                        triangle=triangle,vertex=int(key[3]),angular_fraction=float(at if np.isscalar(at) else at[key]))
    # Exact piecewise-polynomial meridional extrema, including its quartic waist mode.
    s=shape.reference
    knots=s.x/shape.unit_length
    poly=Polynomial([0.,1.])
    u=poly*(1-poly)
    pert=(shape.exit/shape.neck-1.5)*poly+4*shape.params[2]*u+4*shape.params[3]*u*(2*poly-1)+16*shape.extra*u*u
    coeff=np.vstack((np.zeros(len(knots)-1),s.c*(shape.unit_length**np.arange(3,-1,-1))[:,None]))
    from math import factorial
    for k in range(5):
        coeff[4-k]+=pert.deriv(k)(knots[:-1])/factorial(k)
    pp=PPoly(shape.neck*coeff,knots,extrapolate=False)
    ids=pp.derivative().roots(extrapolate=False)
    bmin=float(np.min(pp(np.r_[0.,1.,ids[np.isfinite(ids)]])))
    slope_bound=bmin+min(0.,3.*epsilon,-1.5*epsilon)
    bounds_error=max(0.,-ranges[0][0],ranges[0][1]-1.,-ranges[1][0],ranges[1][1]-1.,-ranges[2][0],ranges[2][1]-np.pi/2)
    return dict(prism_chart_determinant_range=extremes,extremizers=location,placement_node_ranges=ranges,
                meridian_minimum=bmin,physical_radial_chart_derivative_lower_bound=slope_bound,
                max_parameter_box_overshoot=float(bounds_error),all_open_cells_oriented=bool(
                    bounds_error<=2e-13 and slope_bound>0 and extremes[0]>0))


def rim_neighborhood(shape,epsilon,body,solid):
    """Find material preimages of a physical target set identical on both grids."""
    targets=[]
    for delta,depth,angle in product((.05,.0722,.09),(.997,.9989,.9995),(.09,.48,1.12,1.53)):
        radius=shape.exit-delta
        z=depth*shape.length
        def err(y):
            x,t,a=np.array(y).reshape(3,1)
            r=material_map(shape,epsilon,body,solid,t,x,a)[0]
            return r-[radius,z,angle]
        sol=least_squares(err,[np.clip(radius/shape.radius(z)[0],0.,1.),depth,angle],
                          bounds=([0.,0.,0.],[1.,1.,np.pi/2]),xtol=1e-14,ftol=1e-14,gtol=1e-14)
        targets.append((delta,depth,angle,sol.x,float(abs(sol.fun).max())))
    return targets


def solid_readout(state_path, shape, epsilon, body, solid):
    """Sample both native triangles per cell and cluster extra sites inside the disk rim.

    Jacobians below differentiate the complete physical and relaxed positions,
    using centered differences of those maps. SVD of their quotient includes
    meridional, polar and azimuthal stretches. Sites cover both materials,
    with separate reports for their shared interface and the terminal edge.
    """
    chart, labels, _, _ = read(state_path)
    init_shape, _, init_registration = preparation(4)
    if len(chart.co) != 1 or np.shape(solid) != (len(body.t),len(body.x),len(body.a),3):
        raise ValueError('Invalid R4 registered lattice')
    p = np.exp(chart.co[0])
    cut = np.exp(-p/2)

    def reference(t,x,angle):
        core = x < cut
        out = np.empty((len(x),3))
        out[:,2] = angle
        b0 = init_shape.profile(t[core])[0]
        p0 = np.exp(legval(2*t[core]-1,init_registration))
        out[core,0] = x[core]*b0*np.exp(-p0/2)/cut
        out[core,1] = t[core]*init_shape.length
        if np.any(~core):
            e = -np.log(x[~core])/p
            ids, bary, _, _ = tri(chart.a,chart.t,len(chart.a),len(chart.t),e,t[~core])
            out[~core,:2] = np.einsum('ij,ijk->ik',bary,labels.reshape(-1,2)[ids])
            out[~core,0] += shape.neck
        return out

    ii,jj,kk = np.indices((len(body.x)-1,len(body.t)-1,len(body.a)-1))
    cells = (ii.ravel(),jj.ravel(),kk.ravel())
    edge = (ii == ii.max()) & (jj == jj.max())
    sampling = [(cells,a,b,c) for a,b in ((.21,.27),(.76,.81)) for c in (.22,.78)]
    # Approach the corner from inside both one-sided neighboring triangles.
    sampling += [((ii[edge],jj[edge],kk[edge]),a,b,c)
                 for a,b in ((.975,.99),(.99,.975),(.985,.009),(.009,.985))
                 for c in (.02,.5,.98)]
    targets=rim_neighborhood(shape,epsilon,body,solid)
    material,errors=np.stack([v[3] for v in targets]).T,np.array([v[4] for v in targets])
    x,t,theta=material
    axes=body.x,body.t,body.a
    i,j,k=[np.clip(ax.searchsorted(m)-1,0,len(ax)-2) for ax,m in zip(axes,material)]
    inverse_fractions=tuple((m-ax[idx])/np.diff(ax)[idx] for m,ax,idx in zip(material,axes,(i,j,k)))
    sampling.append(((i,j,k),*inverse_fractions))
    totals = {}
    target_values=[]
    for batch,((i,j,k),a,b,c) in enumerate(sampling):
        x = body.x[i]+a*np.diff(body.x)[i]
        t = body.t[j]+b*np.diff(body.t)[j]
        theta = body.a[k]+c*np.diff(body.a)[k]
        actual = lambda t,x,w: material_map(shape,epsilon,body,solid,t,x,w)
        def jacobian(fn):
            center = fn(t,x,theta)
            columns = []
            for it, step in enumerate((np.diff(body.x)[i],np.diff(body.t)[j],np.diff(body.a)[k])):
                dx,dt,da = (int(q==it)*8e-7*step for q in range(3))
                deriv=(fn(t+dt,x+dx,theta+da)-fn(t-dt,x-dx,theta-da))/(1.6e-6*step[:,None])
                deriv[:,2] *= center[:,0]
                columns.append(deriv)
            return center,np.stack(columns,axis=-1)
        current,deriv=jacobian(actual)
        zero,basis=jacobian(reference)
        F=deriv@np.linalg.inv(basis)
        hoop=current[:,0]*np.sin(current[:,2])/(zero[:,0]*np.sin(theta))
        stretches=np.linalg.svd(F,compute_uv=False)
        J=np.linalg.det(deriv)/np.abs(np.linalg.det(basis))*hoop
        r,z,th=current.T
        regions=dict(all=np.ones(len(x),bool), inside_disk_edge=(i==len(body.x)-2)&(j==len(body.t)-2),
                     bonded_join=(i==body.cut_node-1)|(i==body.cut_node))
        if batch==len(sampling)-1:
            good=errors<1e-9
            regions={key:value&good for key,value in regions.items()}
            regions['fixed_physical_neighborhood']=good
            for ktarget,(dr,depth,polar,_,mismatch) in enumerate(targets):
                item=dict(edge_minus_radius=dr,z_over_depth=depth,angle=polar,placement_error=mismatch)
                if good[ktarget]:
                    item.update(min_stretch=float(min(stretches[ktarget,-1],abs(hoop[ktarget]))),
                                J=float(J[ktarget]))
                target_values.append(item)
        for key, mask in regions.items():
            ix=np.flatnonzero(mask)
            if not len(ix):continue
            summary=totals.setdefault(key,dict(sites=0,minimum_stretch=None,maximum_stretch=None,
                                               minimum_J=None,maximum_J=None))
            summary['sites']+=len(ix)
            for label,values,fn,select in (
                    ('minimum_stretch',np.minimum(stretches[:,-1],abs(hoop)),np.argmin,min),
                    ('maximum_stretch',np.maximum(stretches[:,0],abs(hoop)),np.argmax,max),
                    ('minimum_J',J,np.argmin,min),('maximum_J',J,np.argmax,max)):
                n=ix[fn(values[ix])]
                val=float(values[n])
                if summary[label] is None or select(summary[label],val)==val:
                    summary[label]=val
                    summary[label+'_site']=dict(material_t_x_angle=[t[n],x[n],theta[n]],
                                                physical_r_z_angle=[r[n],z[n],th[n]])
    totals['physical_neighborhood_sites']=target_values
    totals['orientation']=orientation(shape,epsilon,body,solid)
    return totals


def fluxes(model, u):
    """Unassembled variational port currents from the physical phase gradient."""
    n=(model.D[:,4]*(model.V@u[0]))**2
    p=np.column_stack([m@u[1] for m in model.G])
    c=model.ref[:,1:]+np.column_stack([m@u[2] for m in model.G])
    v=p-c*(ORDER_INERTIA*dot(p,c)/(n+ORDER_INERTIA*dot(c,c)))[:,None]
    current=n[:,None]*v
    wt=model.D[:,0]
    currents={}
    for key,ids,sign in (('feed',model.feed,-1),('bulk',model.bulk,1)):
        test=np.zeros(model.size)
        test[ids]=1.
        currents[key]=float(sign*sum((wt*(m@test))@current[:,i] for i,m in enumerate(model.G)))
    currents['imbalance']=currents['feed']-currents['bulk']
    currents['phase_at_ports']=[float(u[1,model.feed].mean()),float(u[1,model.bulk].mean())]
    gradient = model.evaluate(u)[1]
    currents['stationarity'] = float(np.max(abs(model.norm(model.vars(gradient)))))
    currents['raw_free_gradient_max'] = {
        name: float(np.max(abs(row[free])))
        for name, row, free in zip(('amplitude', 'phase', 'order'), gradient, model.masks)
    }
    return currents


def physical_patch(reader, slab, bank, x, a, z):
    """Independent Q1 gradient in cylindrical coordinates for the volume-work test."""
    m=reader.model
    axes=(m.rr[bank],m.a,m.z[slab])
    point=(x,a,np.broadcast_to(z,x.shape))
    base=[np.clip(v.searchsorted(y,side='right')-1,0,len(v)-2) for v,y in zip(axes,point)]
    gap=[np.diff(v)[j] for v,j in zip(axes,base)]
    frac=[(y-v[j])/d for y,v,j,d in zip(point,axes,base,gap)]
    R,slope=m.shape.radius(z)
    width=28.-R if bank else R
    r=R+width*x if bank else width*x
    F,Fp=reader.vortex(r*np.sin(a))
    tr=np.zeros(len(x));tz=np.zeros(len(x))
    factor_r=np.zeros(len(x));factor_z=np.zeros(len(x))
    for sign,interp in zip((1.,-1.),reader.fields[slab,bank]):
        val=np.zeros((len(x),3))
        grad=np.zeros((len(x),3,3))
        for bits in product((0,1),repeat=3):
            weights=[f if side else 1-f for f,side in zip(frac,bits)]
            N=interp.values[tuple(i+v for i,v in zip(base,bits))]
            val+=np.prod(weights,axis=0)[:,None]*N
            for k in range(3):
                grad[:,:,k]+=(((2*bits[k]-1)/gap[k])*weights[(k+1)%3]*weights[(k+2)%3])[:,None]*N
        grad[:,:,0]/=width
        grad[:,:,1]/=r[:,None]
        grad[:,:,2]-=(slope*(1-x if bank else x))[:,None]*grad[:,:,0]
        state=physical(m.winding,r,a,val.T,grad.transpose(1,0,2),F,Fp)
        vr,vz=stresses(state)
        tr+=sign*vr.sum(axis=0)
        tz+=sign*vz.sum(axis=0)
        gf,_,_,density,swirl,_=state
        gv=grad[:,0]
        y=val[:,0]**2
        interaction=F**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20
        lag_change=.5*(F**2*dot(gv,gv)-dot(gf,gf)-density*swirl)+interaction-(density**5-5*density+4)/20
        factor_r+=sign*(F**2*gv[:,0]*gv[:,2]-gf[:,0]*gf[:,2])
        factor_z+=sign*(F**2*gv[:,2]**2-gf[:,2]**2-lag_change)
    return tr,tz,factor_r,factor_z


def segment(axis, lo, hi, degree, extra=()):
    knots=np.r_[axis,extra]
    return quadrature(np.r_[lo,knots[(lo<knots)&(knots<hi)],hi],degree)


def shell_work(reader, near, far, radial, degree):
    """Virtual work of the physical stress for near minus far axial displacement."""
    m=reader.model
    L=m.shape.length
    lo,hi=radial
    ang,wa=quadrature(m.a,degree)
    def support(s,edges):
        a,b=edges
        return np.clip((b-s)/(b-a),0,1), np.where((s>a)&(s<b),-1/(b-a),0.)
    # Split all axial taper knots. Split radial knots independently at each
    # depth because a fixed physical cylinder cuts moving chart columns.
    terms=np.zeros(4)
    for slab,side in ((1,-1.),(2,1.)):
        edges=L+side*np.array([*near,*far])
        zq,wz=segment(m.z[slab],edges.min(),edges.max(),degree,edges)
        for z,dz in zip(zq,wz):
            az,daz=support(abs(z-L),near)
            af,daf=support(abs(z-L),far)
            A=az-af
            Az=side*(daz-daf)
            R=m.shape.radius(z)[0]
            for bank in (0,1):
                width=28.-R if bank else R
                offset=R if bank else 0.
                boundaries=(np.array(radial)-offset)/width
                xs,wx=segment(m.rr[bank],0.,min(1.,boundaries[1]),degree,boundaries)
                xx,aa=np.meshgrid(xs,ang,indexing='ij')
                x,a=xx.ravel(),aa.ravel()
                r=offset+width*x
                rad,dr=support(r,radial)
                wt=(4*np.pi*width*dz*(wx[:,None]*wa).ravel()*r**2*np.sin(a))
                tr,tz,factor_r,factor_z=physical_patch(reader,slab,bank,x,a,z)
                terms[:3] += [wt@(tr*dr*A),wt@(tz*rad*Az),
                              wt@(factor_r*dr*A+factor_z*rad*Az)]
            if slab==1:
                # Both duct banks have the same radius and separate amplitude traces.
                # The radial taper equals one on their entire common wall.
                F,Fp=reader.vortex(R*np.sin(ang))
                for bank,sign in ((0,1.),(1,-1.)):
                    coords=np.column_stack((np.full(len(ang),1-bank),ang,np.full(len(ang),z)))
                    u,v=[q(coords)[:,0] for q in reader.fields[1,bank]]
                    terms[3]-=sign*2*np.pi*R**2*dz*Az*(wa@(np.sin(ang)**2*F*Fp*(u*u-v*v)))
    return terms


def face_work(reader, windows, radial, degree):
    """Direct, separately differentiated physical flux on two caps and a cylinder."""
    m=reader.model
    L=m.shape.length
    lo,hi=radial
    a,wa=quadrature(m.a,degree)
    faces={}
    for slab,side in ((1,-1.),(2,1.)):
        a0,b0=windows
        lim=np.sort(L+side*np.array(windows))
        zs,wz=segment(m.z[slab],*lim,degree)
        cap=np.zeros(2)
        for z,w in zip(zs,wz):
            b=m.shape.radius(z)[0]
            for bank in (0,1):
                mesh=b+(28.-b)*m.rr[1] if bank else b*m.rr[0]
                r,wr=segment(mesh,mesh[0],min(mesh[-1],hi),degree,[lo])
                rr,aa=np.meshgrid(r,a,indexing='ij')
                weight=4*np.pi*rr**2*np.sin(aa)*(wr*np.clip((hi-r)/(hi-lo),0,1))[:,None]*wa
                cap+=w/(b0-a0)*(weight.ravel()@reader.sample(slab,bank,rr,aa,z,np.array([0,0,side])))
        faces['lower_cap' if slab==1 else 'upper_cap']=cap
        edge=L+side*b0
        zs,wz=segment(m.z[slab],min(edge,L),max(edge,L),degree,[L+side*a0])
        cylinder=np.zeros(2)
        for z,w in zip(zs,wz):
            b=m.shape.radius(z)[0]
            mesh=b+(28.-b)*m.rr[1]
            r,wr=segment(mesh,lo,hi,degree)
            rr,aa=np.meshgrid(r,a,indexing='ij')
            weight=4*np.pi*rr**2*np.sin(aa)*wr[:,None]*wa
            A=np.clip((b0-abs(z-L))/(b0-a0),0,1)
            cylinder+=w*A/(hi-lo)*(weight.ravel()@reader.sample(slab,1,rr,aa,z,np.array([1,0,0])))
        faces['lower_side' if slab==1 else 'upper_side']=cylinder
    return {k:v.tolist() for k,v in faces.items()},sum(faces.values())


def enclosed(model, base_order, wound, rest):
    """Record current matched-geometry support or the actual reason a cut is invalid."""
    near,far,radial=(.145,.17),(.195,.22),(4.45,5.05)
    L=model.shape.length
    info=dict(near_half_heights=near,far_half_heights=far,physical_radial_taper=radial,
              order_halfwidth=model.h,wall_polar_support_end=.62*L)
    lower=L-far[1]
    if lower<=.62*L or near[0]<=model.h:
        return dict(**info,resolved=False,reason='Physical shell intersects the polar chart or mobile order')
    if extrema(model.shape,0.,1.)['maximum']>=radial[0]:
        # A wide entrance away from this shell is allowed.
        if model.shape.radius(np.linspace(lower,L,1025))[0].max()>=radial[0]:
            return dict(**info,resolved=False,reason='Exterior cylinder crosses the current wall')
    dc=wound[2]-rest[2]
    outside=list(np.ravel(dc[model.ids[1,0][:,:,model.z[1]<L-model.h-1e-12]]))
    outside.extend(np.ravel(dc[model.ids[2,0][:,:,model.z[2]>L+model.h+1e-12]]))
    max_off_band=max(np.max(abs(np.asarray(outside))),
                     *(np.max(abs(dc[model.ids[slab,1]])) for slab in (1,2)))
    if max_off_band!=0.:
        return dict(**info,resolved=False,reason='Differential crystal order reaches the outer shell')
    reader=Reader(model,base_order,wound,rest)
    readings=[]
    for degree in (2,4):
        fluxes={}
        for label,bounds in (('near',near),('far',far)):
            faces,out=face_work(reader,bounds,radial,degree)
            left,right=bounds
            # Integrate only on the two interior banks. Splitting at the end
            # of the ramp preserves physical taper width on nonaligned meshes.
            slit_flat=reader.sealed_increment(0.,left,degree)
            slit_taper=reader.sealed_increment(left,right,degree,L-left)
            slit=slit_flat+slit_taper
            fluxes[label]=dict(faces=faces,outward_mass=float(out[0]),negative_outward_axial=float(-out[1]),
                               slit_mass_axial=slit.tolist(),outward_plus_slit_mass_axial=(out+slit).tolist())
        shell=shell_work(reader,near,far,radial,degree)
        near_force=fluxes['near']['negative_outward_axial']
        far_force=fluxes['far']['negative_outward_axial']
        readings.append(dict(degree=degree,**fluxes,
                             work_radial_axial=shell[:2].tolist(),far_plus_shell=float(far_force+shell[:2].sum()),
                             closure=float(near_force-far_force-shell[:2].sum()),
                             native_vortex_volume_wall=shell[2:].tolist(),differential_crystal_shell_work=0.,
                             far_plus_native_shell=float(far_force+shell.sum()),
                             near_minus_far_plus_native_shell=float(near_force-far_force-shell.sum())))
    return dict(**info,resolved=True,maximum_off_band_order_change=max_off_band,measurements=readings,
                convention='wound drive minus wound rest; force is minus outward axial momentum')


def audit_endpoint(base_order, model, wound, rest, unwound_model, unwound, *, state_path, body, solid):
    if model.winding!=1 or unwound_model.winding!=0 or not np.array_equal(model.shape.params,unwound_model.shape.params):
        raise ValueError('Winding or shape mismatch in same-material controls')
    if (model.polar_shift!=unwound_model.polar_shift or model.shape.extra!=unwound_model.shape.extra or
        any(not np.array_equal(a,b) for a,b in zip((*model.rr,*model.z,model.a),(*unwound_model.rr,*unwound_model.z,unwound_model.a))) or
        not np.array_equal(model.material_load,unwound_model.material_load)):
        raise ValueError('The two windings require the same geometry, axes and crystal deposition')
    controls={}
    for name,m,u,bias in (('wound',model,wound,BIAS),('rest',model,rest,0.),('unwound',unwound_model,unwound,BIAS)):
        if (not np.all(u[1,m.bulk]==bias) or np.any(u[1,m.feed]) or np.any(u[0,np.r_[m.feed,m.bulk]]!=1.)):
            raise ValueError(f'Inconsistent prescribed ports for {name}')
        controls[name]=fluxes(m,u)
    disks=[]
    for side in (-1,1):
        for degree in (6,12):
            w,n,mach,j,_=face(model,base_order,wound,side,degree)
            pos=w*np.maximum(j,0.)
            dense=(n>.5)&(mach<1.)
            disks.append(dict(side=side,degree=degree,signed_current=float(w@j),
                              positive_current=float(pos.sum()),dense_subsonic_current=float(pos[dense].sum()),
                              dense_subsonic_fraction=float(pos[dense].sum()/pos.sum()) if pos.sum()>0 else None))
    plain=controls['unwound']['feed']
    return dict(geometry=geometry(model.shape,model.polar_shift),controls=controls,physical_disk=disks,
                wound_fraction_of_unwound_intake=controls['wound']['feed']/plain if plain>0. else None,
                solid=solid_readout(state_path,model.shape,model.polar_shift,body,solid),
                enclosure=enclosed(model,base_order,wound,rest))
