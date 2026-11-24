"""One conforming angular edge split of the held R4 material map.

All interpolation is of registered material placements, including their
original reference labels. The sole independent added vertex has three
displacement components. Transition cells partition their interface faces
without changing their volume interpolation.
"""
from collections import Counter
from fractions import Fraction as Q

import numpy as np

from mouth_active_prism_orientation import (
    _Support, _certify_wedge, _coefficients, _down, _easy_wedges,
    _exact_array, _possible_fold, _quadratic_min, _radial_bound, _rational, _up)
from mouth_prism_support_decision import decide


CENTER = (5, 18)
TRIANGLES = (((0, 0), (0, 1), (1, 0)), ((1, 1), (1, 0), (0, 1)))
STAR = {(j, i, 0, triangle): offsets.index((5-j, 18-i))
        for j in (4, 5) for i in (17, 18) for triangle, offsets in enumerate(TRIANGLES)
        if (5-j, 18-i) in offsets}


def _geometry(body, solid, key):
    j, i, k, tri = key
    offsets = TRIANGLES[tri]
    ref = _exact_array([[(body.x[i+x], body.t[j+t], body.a[k+y]) for t, x in offsets] for y in (0, 1)])
    points = ref + _exact_array([[solid[j+t, i+x, k+y] for t, x in offsets] for y in (0, 1)])
    scale = (_rational(body.x[i+1])-_rational(body.x[i])) * (_rational(body.t[j+1])-_rational(body.t[j])) * (_rational(body.a[k+1])-_rational(body.a[k]))
    return ref, points, scale


def seam(body):
    """The native crystal angular integration boundary nearest the layer midpoint."""
    edges = np.asarray(body.buffers['angle_edges'])
    internal = edges[(edges > body.a[0]) & (edges < body.a[1])]
    if not len(internal):
        raise ValueError('No internal production quadrature panel in angular strip')
    return float(internal[np.argmin(abs(internal - (body.a[0]+body.a[1])/2))])


def seam_fraction(body):
    return (_rational(seam(body))-_rational(body.a[0]))/(_rational(body.a[1])-_rational(body.a[0]))


def _children(key, ref, points, surplus, fraction):
    rm = (1-fraction)*ref[0]+fraction*ref[1]
    pm = (1-fraction)*points[0]+fraction*points[1]
    pm[STAR[key]] += _exact_array(surplus)
    return ((np.stack((ref[0], rm)), np.stack((points[0], pm))),
            (np.stack((rm, ref[1])), np.stack((pm, points[1]))))


def _cert(points, scale, support, activation, counts, depth):
    status, witness = _certify_wedge(points, scale, support, activation, counts, depth)
    if status == 'unresolved':
        status, witness = decide(points, scale, support)
        if status == 'certified' and activation is not None:
            minima = [_quadratic_min(row)[0] for row in _coefficients(points)]
            activation.update(support.ids[support.box(_possible_fold(minima) @ points)].ravel().tolist())
    return status, witness


def refined_orientation(body, solid, surplus, order, *, collect_activation=False, max_depth=12):
    """Exact positive-order support on unchanged prisms and both refined halves."""
    axes = [np.asarray(a) for a in (body.x, body.t, body.a)]
    grid = np.stack(np.broadcast_arrays(axes[0][None, :, None], axes[1][:, None, None], axes[2][None, None, :]), axis=-1)
    solid, surplus = np.asarray(solid), np.asarray(surplus)
    if (solid.shape != grid.shape or surplus.shape != (3,) or max_depth < 0 or np.ndim(order) != 1
            or not np.all(np.isfinite(solid)) or not np.all(np.isfinite(surplus)) or np.any(body.buffers['base'])):
        raise ValueError('Invalid R4 material displacement or order input')
    support = _Support(body, order)
    if (any((v[0], v[-1]) != (w[0], w[-1]) or np.any(np.diff(v) <= 0) for v, w in zip(axes, support.axes))
            or [a[0] for a in axes] != [0, 0, 0] or [a[-1] for a in axes] != [1, 1, np.pi/2]):
        raise ValueError('Material and order domains disagree')
    radial = _radial_bound(body)
    activation = set() if collect_activation else None
    counts = dict(total_prisms=int(2*np.prod(np.array(grid.shape[:3])-1))+len(STAR),
                  positive_bernstein_prisms=0, exact_prisms=0, exact_subprisms=0,
                  inactive_fold_enclosures=0, max_depth_used=0, active_fold_prisms=0,
                  unresolved_prisms=0, refined_children_checked=0)
    def result(status, witness=None):
        return dict(status=status, all_active_cells_oriented=status == 'certified', witness=witness,
                    physical_radial_chart_derivative_lower_bound=radial, counts=counts,
                    potential_activation_nodes=sorted(activation) if activation is not None else None)
    if not np.isfinite(radial) or not np.isfinite(body.shape.length) or min(radial, body.shape.length) <= 0:
        return result('unresolved', {'reason': 'Physical radial chart multiplier unverified'})
    coords = grid+solid
    for n, axis in enumerate(axes):
        x = coords[..., n]
        outside = (x < axis[0]) | (x > axis[-1])
        for pos in np.argwhere(((x == axis[0]) | (x == axis[-1])) & (solid[..., n] != 0)):
            key = (*pos, n)
            val = _rational(grid[key])+_rational(solid[key])
            if not support.axes[n][0] <= val <= support.axes[n][-1]:
                outside[tuple(pos)] = True
        if outside.any():
            return result('outside_chart', {'node_t_x_angle':np.argwhere(outside)[0].tolist(),'component':n})
    def check(key, points, scale, half=None):
        counts['exact_prisms'] += 1
        status, witness = _cert(points, scale, support, activation, counts, max_depth)
        if status == 'certified':
            return None
        counts[status+'_prisms'] += 1
        w = dict(cell_t_x_angle=list(key[:3]), triangle=key[3], **witness)
        if half is not None:
            w['child_half'] = half
            if 'angular_fraction' in w:
                w['child_angular_fraction'] = w['angular_fraction']
                w['angular_fraction'] = float((0 if half == 0 else fraction) + (fraction if half == 0 else 1-fraction)*Q(w['child_angular_fraction']))
        return result(status, w)
    fraction = seam_fraction(body)
    for key in STAR:
        ref, pts, scale = _geometry(body, solid, key)
        for half, (_, pc) in enumerate(_children(key, ref, pts, surplus, fraction)):
            if any(pc[..., n].min() < axis[0] or pc[..., n].max() > axis[-1] for n, axis in enumerate(support.axes)):
                return result('outside_chart', {'cell_t_x_angle':list(key[:3]),'triangle':key[3],'child_half':half})
            counts['refined_children_checked'] += 1
            error = check(key, pc, scale*(fraction if half == 0 else 1-fraction), half)
            if error:
                return error
    for tri, offsets in enumerate(TRIANGLES):
        v = np.stack([coords[j:j+len(body.t)-1, i:i+len(body.x)-1] for j, i in offsets], axis=-2)
        lo, hi = _down(v), _up(v)
        easy = _easy_wedges((lo[:, :, :-1], hi[:, :, :-1]), (lo[:, :, 1:], hi[:, :, 1:]))
        for j, i, k, triangle in STAR:
            if tri == triangle:
                easy[j,i,k] = True
        counts['positive_bernstein_prisms'] += int(easy.sum())-sum(tri == s[3] for s in STAR)
        for ind in np.argwhere(~easy):
            key = (*map(int, ind), tri)
            _, points, scale = _geometry(body, solid, key)
            error = check(key, points, scale)
            if error:
                return error
    return result('certified')


def connectivity_audit(body):
    """Verify closed cells and reverse-oriented paired facets on the entire mesh.

    The neighbors of the patch carry identical midpoint edge incidences.
    A transition face with one extra collinear vertex is a pentagon. With
    two such vertical edges it is represented by two quadrilateral faces.
    """
    marked = {(j+t,i+x) for j,i,k,tri in STAR for t,x in TRIANGLES[tri]}
    maxt, maxx, maxk = len(body.t)-1, len(body.x)-1, len(body.a)-1
    facets = {}
    cells = subdivided_faces = transition_cells = 0
    def add_cell(faces):
        nonlocal cells
        cells += 1
        edges = Counter()
        for face in faces:
            for a,b in zip(face, face[1:]+face[:1]):
                edges[a,b] += 1
            key = tuple(sorted(face))
            facets.setdefault(key, []).append(face)
        if any(v != 1 or edges[b,a] != 1 for (a,b),v in edges.items()):
            raise ArithmeticError('An oriented cell boundary has an open edge')
    for j,i,k,tri in np.ndindex(maxt,maxx,maxk,2):
        verts = [(j+t,i+x) for t,x in TRIANGLES[tri]]
        split = (j,i,k,tri) in STAR
        slabs = [(2*k,2*k+1),(2*k+1,2*k+2)] if split else [(2*k,2*k+2)]
        if k==0 and not split and any(v in marked for v in verts):
            transition_cells += 1
        for start,end in slabs:
            bot, top = tuple((*v,start) for v in verts),tuple((*v,end) for v in verts)
            faces = [bot[::-1],top]
            for v in range(3):
                w=(v+1)%3
                sb,st,tb,tt=bot[v],top[v],bot[w],top[w]
                sp=(*verts[v],1) if k==0 and not split and verts[v] in marked else None
                tp=(*verts[w],1) if k==0 and not split and verts[w] in marked else None
                if sp and tp:
                    faces += [(sb,tb,tp,sp),(sp,tp,tt,st)]
                    subdivided_faces += 1
                else:
                    faces.append((sb,tb)+((tp,) if tp else ())+(tt,st)+((sp,) if sp else ()))
            add_cell(faces)
    exterior = interior = 0
    for key,uses in facets.items():
        if len(uses)==1 and any(all(v[d]==bound for v in key) for d,hi in enumerate((maxt,maxx,2*maxk)) for bound in (0,hi)):
            exterior += 1
        elif len(uses)==2:
            v,w=uses[0],uses[1][::-1]
            if not any(v==w[n:]+w[:n] for n in range(len(w))):
                raise ArithmeticError('Opposing facets have inconsistent orientation')
            interior += 1
        else:
            raise ArithmeticError(f'Nonconforming facet {key}')
    return dict(parent_prisms_split=len(STAR), children=2*len(STAR), material_prisms=cells,
                angular_midpoint_vertices=len(marked), independent_midpoint_vertices=1,
                transition_cells=transition_cells, transition_quadrilateral_partitions=subdivided_faces,
                paired_interior_facets=interior, exterior_facets=exterior)


def child_witness(body, solid, surplus, order, witness):
    """Exact chart and certified physical determinant at a parent-chart coordinate."""
    from mouth_ordered_j_bounds import IV, iv, profile, reference_factors
    key=(*witness['cell_t_x_angle'],witness['triangle'])
    ref,pts,scale = _geometry(body,solid,key)
    b=witness['barycentric']
    bary=np.array([1-Q(b[1])-Q(b[2]),Q(b[1]),Q(b[2])],object)
    y=Q(witness['angular_fraction'])
    fraction=seam_fraction(body)
    half=int(y>=fraction) if key in STAR else None
    if half is not None:
        ref,pts=_children(key,ref,pts,surplus,fraction)[half]
        lo,hi=(Q(0),fraction) if half == 0 else (fraction,Q(1))
        y,scale=(y-lo)/(hi-lo),scale*(hi-lo)
    X=bary@((1-y)*pts[0]+y*pts[1])
    x,t,a=bary@((1-y)*ref[0]+y*ref[1])
    det=bary@_coefficients(pts)@np.array([(1-y)**2,2*y*(1-y),y*y],object)/scale
    chi=_Support(body,order).at(X)
    xx,tt,aa=X
    h=(tt-Q(.06))/Q(.56)
    f=16*h*h*(1-h)**2*Q(body.epsilon) if 0<h<1 else Q(0)
    p=(3*IV.cos(iv(aa))**2-1)/2
    b0=profile(body,tt,False)
    R=iv(xx*b0)+iv(xx**3*f)*p
    rx=iv(b0)+iv(3*xx**2*f)*p
    angular=(iv(bary@(pts[1,:,2]-pts[0,:,2])/(ref[1,0,2]-ref[0,0,2]))
             if a==0 else IV.sin(iv(aa))/IV.sin(iv(a)))
    vals=[rx*iv(body.data.length)*R**2*angular*iv(det)/(S**2*d) for S,d in reference_factors(body,x,t,key[1]<body.cut_node)]
    if not vals:
        raise ValueError('Missing reference metric in physical determinant witness')
    j_bounds=[float(_down(min(float(v.a) for v in vals))),float(_up(max(float(v.b) for v in vals)))]
    j_value=j_bounds[0]/2+j_bounds[1]/2
    dx=_rational(body.x[key[1]+1])-_rational(body.x[key[1]])
    dt=_rational(body.t[key[0]+1])-_rational(body.t[key[0]])
    sign=1-2*key[3]
    basis=(1-y)*pts[0]+y*pts[1]
    M=np.asarray([sign*(basis[1]-basis[0])/dx, sign*(basis[2]-basis[0])/dt,
                  bary@(pts[1]-pts[0])/(ref[1,0,2]-ref[0,0,2])],float).T
    W,hoop,stretches=_energy_at_witness(body, key, x,t,a,xx,tt,aa,M,j_value)
    return dict(cell_t_x_angle=list(key[:3]), triangle=key[3], barycentric=list(map(float,bary)),
                angular_fraction=float(witness['angular_fraction']), child_half=half, child_angular_fraction=float(y),
                placement_x_t_angle=list(map(float,X)), order=float(chi), chi=float(chi), chart_J=float(det),
                J_interval=j_bounds, J=j_value, W=W, chi_W=None if W is None else float(chi)*W,
                hoop=hoop, principal_stretches=stretches)


def _energy_at_witness(body,key,x,t,a,xx,tt,aa,M,J):
    """Evaluate the unchanged four-stretch energy; use the exact-chart J to take log."""
    from mouth_ordered_j_bounds import profile
    def profile_derivative(z, reference):
        spline=body.buffers['spline']
        k=min(len(spline)-1,max(0,int(len(spline)*z)))
        h=z-Q(k,len(spline))
        c,b,d,_=map(Q,spline[k])
        e,m,u,v=map(Q,body.buffers['shape'][4:] if reference else body.buffers['shape'][:4])
        return 4*((3*c*h+2*b)*h+d+e+4*m*(1-2*z)+4*u*(-6*z*z+6*z-1)+32*v*z*(1-z)*(1-2*z))
    xi,ti,ai=map(float,(xx,tt,aa))
    b=float(profile(body,tt,False))
    bt=float(profile_derivative(tt,False))
    h=(ti-.06)/.56
    f=16*h*h*(1-h)**2*body.epsilon if 0<h<1 else 0.
    ft=32*h*(1-h)*(1-2*h)*body.epsilon/.56 if 0<h<1 else 0.
    c,s=np.cos(ai),np.sin(ai)
    p=(3*c*c-1)/2
    radius=xi*b+xi**3*f*p
    C=np.array([[b+3*xi**2*f*p,xi*bt+xi**3*ft*p,-3*xi**3*f*c*s],
                [0.,body.data.length,0.],[0.,0.,radius]])
    if key[1]<body.cut_node:
        sx=float(profile(body,t,True))*body.data.core_ratio
        st=float(x*profile_derivative(t,True))*body.data.core_ratio
        tx,ttref=0.,body.data.ref_length
        S=float(x)*sx
    else:
        er,tr=body.buffers['reg_e'],body.buffers['reg_t']
        lab=body.buffers['labels'].reshape(len(tr),len(er),2)
        xi_ref,ti_ref=float(x),float(t)
        e=-np.log(xi_ref)/body.data.pressure
        i=min(len(er)-2,max(0,int(np.searchsorted(er,e,side='right'))-1))
        j=min(len(tr)-2,max(0,int(np.searchsorted(tr,ti_ref,side='right'))-1))
        de,dt=er[i+1]-er[i],tr[j+1]-tr[j]
        u,v=(e-er[i])/de,(ti_ref-tr[j])/dt
        up=int(u+v>1)
        z=lab[j+up,i+up]
        le,lt=lab[j+up,i+1-up]-z,lab[j+1-up,i+up]-z
        S=4+z[0]+(1-u if up else u)*le[0]+(1-v if up else v)*lt[0]
        sx,tx=(1-2*up)*le/(-body.data.pressure*xi_ref*de)
        st,ttref=(1-2*up)*lt/dt
    inv_ref=np.linalg.inv([[sx,st,0.],[tx,ttref,0.],[0.,0.,S]])
    F=C@M@inv_ref
    hoop=radius*(M[2,2] if a==0 else s/np.sin(float(a)))/S
    stretches=sorted([*map(float,np.linalg.svd(F,compute_uv=False)),float(abs(hoop))])
    W=float((np.sum(F*F)+hoop*hoop-4)/2-np.log(J)+(J-1)**2/4) if J>0 else None
    return W,float(hoop),stretches


refined_witness=child_witness


def _exact_trace_audit(body, solid, fraction):
    """Compare every mapped face trace adjoining the edge, including the new basis.

    Each restricted side face is bilinear, hence its four control values
    determine the complete trace. Caps are affine and have three controls.
    """
    by_edge={}
    for j in range(CENTER[0]-2,CENTER[0]+2):
        for i in range(CENTER[1]-2,CENTER[1]+2):
            for tri,offsets in enumerate(TRIANGLES):
                key=(j,i,0,tri)
                verts=[(j+t,i+x) for t,x in offsets]
                for idx in range(3):
                    a,b=verts[idx],verts[(idx+1)%3]
                    by_edge.setdefault(tuple(sorted((a,b))),[]).append((key,verts))
    shared=exterior=0
    # One independent scalar hat multiplies each displacement component.
    # Setting every component to one therefore exercises all three traces.
    probe=np.ones(3)
    for edge,faces in by_edge.items():
        if not any(key in STAR for key,_ in faces):
            continue
        if len(faces)!=2:
            raise ArithmeticError('Local patch misses a face neighbor')
        for half,(lo,hi) in enumerate(((0,fraction),(fraction,1))):
            baseline=derivative=None
            for key,verts in faces:
                ref,pts,_=_geometry(body,solid,key)
                if key in STAR:
                    _,pc=_children(key,ref,pts,probe,fraction)[half]
                    _,p0=_children(key,ref,pts,np.zeros(3),fraction)[half]
                else:
                    pc=p0=np.stack(((1-lo)*pts[0]+lo*pts[1],(1-hi)*pts[0]+hi*pts[1]))
                select=[verts.index(v) for v in edge]
                controls=p0[:,select]
                tangent=pc[:,select]-controls
                if baseline is not None and (np.any(baseline!=controls) or np.any(derivative!=tangent)):
                    raise ArithmeticError('Mapped shared-face trace or new displacement mode is discontinuous')
                baseline,derivative=controls,tangent
            if all(key in STAR for key,_ in faces):
                shared+=1
            else:
                exterior+=1
    for key in STAR:
        ref,pts,_=_geometry(body,solid,key)
        (_,below),(_,above)=_children(key,ref,pts,probe,fraction)
        if np.any(below[1]!=above[0]) or np.any(below[0]!=pts[0]) or np.any(above[1]!=pts[1]):
            raise ArithmeticError('Refined cap trace discontinuity')
    return dict(star_interface_half_faces=shared, transition_half_faces=exterior,
                in_prism_seams=len(STAR),parent_endcaps=2*len(STAR),
                maximum_exact_placement_trace_difference=0.0,
                maximum_exact_hierarchical_basis_trace_difference=0.0)


def geometry_audit(body, solid, order):
    """Unmoved mapped labels and gradients agree as exact parent polynomials."""
    from bisect import bisect_right
    from itertools import product
    topology=connectivity_audit(body)
    zero=np.zeros(3)
    peak_map=peak_grad=Q(0)
    worst=None
    fraction=seam_fraction(body)
    for key in STAR:
        ref,pts,s=_geometry(body,solid,key)
        for half,(rc,pc) in enumerate(_children(key,ref,pts,zero,fraction)):
            lo,hi=(Q(0),fraction) if half == 0 else (fraction,Q(1))
            for original,child in ((ref,rc),(pts,pc)):
                for k in (0,1):
                    y=lo+k*(hi-lo)
                    expect=(1-y)*original[0]+y*original[1]
                    peak_map=max(peak_map, max(abs(child[k]-expect).flat))
                    peak_grad=max(peak_grad,max(abs((child[k,1:]-child[k,0])-(expect[1:]-expect[0])).flat))
                peak_grad=max(peak_grad,max(abs((child[1]-child[0])/(hi-lo)-(original[1]-original[0])).flat))
            for i,row in enumerate(_coefficients(pc)):
                d,u=_quadratic_min(row)
                d/=s*(hi-lo)
                if worst is None or d<worst[0]:
                    worst=(d,key,half,i,u,pc,s*(hi-lo))
    if peak_map or peak_grad:
        raise ArithmeticError('Hierarchical split changed the parent map')
    d,key,half,i,u,pc,s=worst
    if d>=0:
        raise ValueError('No folded target for the required negative control')
    point=(1-u)*pc[0,i]+u*pc[1,i]
    support=_Support(body,order)
    if support.at(point)!=0:
        raise ValueError('Negative control point was already ordered')
    positions=[]
    for x,axis in zip(point,support.axes):
        lo=min(len(axis)-2,bisect_right(axis,x)-1)
        positions.append([j for j in (lo,lo+1) if (x>axis[lo] if j==lo+1 else x<axis[lo+1])])
    perturbed=np.array(order,copy=True)
    for pos in product(*positions):
        perturbed[support.ids[pos]]=0.001
    neg_order=_Support(body,perturbed)
    status,item=_cert(pc,s,neg_order,None,{'exact_subprisms':0,'max_depth_used':0,'inactive_fold_enclosures':0},12)
    negative=refined_orientation(body,solid,zero,perturbed)
    if neg_order.at(point)!=Q(.001) or status!='active_fold' or negative['status']!='active_fold':
        raise ArithmeticError('Required child-level negative control failed')
    return dict(connectivity=topology, shared_traces=_exact_trace_audit(body, solid, fraction), split_angle=seam(body), split_parent_fraction=float(fraction), exact_max_reference_or_placement_delta=float(peak_map),
                exact_max_reference_or_placement_derivative_delta=float(peak_grad),
                negative_control=dict(cell_t_x_angle=list(key[:3]),triangle=key[3],child_half=half,
                                      child_angular_fraction=float(u),barycentric_vertex=i,
                                      placement_x_t_angle=list(map(float,point)),chart_J=float(d),order=float(neg_order.at(point)),
                                      child_status=status, global_status=negative['status']))
