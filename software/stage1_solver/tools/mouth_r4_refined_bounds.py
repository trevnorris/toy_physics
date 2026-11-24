"""One-sided bounds on ordered-region minima in the six refined R4 parent prisms."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from fractions import Fraction as Q
from math import log

import numpy as np

from mouth_active_prism_orientation import _coefficients, _Support
from mouth_r4_conforming_prisms import STAR, _geometry, _children, seam_fraction, child_witness


ZERO, ONE, HALF = Q(0), Q(1), Q(1,2)
UNIT = tuple(tuple(Q(i==j) for i in range(3)) for j in range(3))
MID = tuple(tuple((v+w)/2 for v,w in zip(UNIT[i],UNIT[(i+1)%3])) for i in range(3))
CENTROID = (Q(1,3),)*3


def _solve_knots(p,q,axes):
    """All positions at spatial Q1 knots along one mapped straight segment."""
    cuts = {ZERO,ONE}
    for i,axis in enumerate(axes):
        if p[i]==q[i]:
            continue
        low,high=sorted((p[i],q[i]))
        for value in axis[bisect_right(axis,low):bisect_left(axis,high)]:
            cuts.add((value-p[i])/(q[i]-p[i]))
    return sorted(cuts)


def _parameters(t,p,q):
    return tuple((1-t)*a+t*b for a,b in zip(p,q))


def _from_witness(item):
    return (*item['cell_t_x_angle'],item['triangle'])


def refined_bounds(body, solid, surplus, order, witnesses=()):
    """Bound infima by attainable ordered-point witnesses after the exact guard.

    The caller must first certify its complete refined state. Any nonpositive
    J found at a strictly ordered evaluation point contradicts that guard and
    aborts this readout. Candidate sites use serialized parent coordinates.
    """
    support=_Support(body,order)
    fraction=seam_fraction(body)
    trials={}
    explicit={}
    boundary_faces=0
    finite_energy_boundary=None
    sheet_nodes = len(body.t)*len(body.x)
    if sheet_nodes!=625:
        raise ValueError('Expected the held coarse R4 grid')
    geometry={}
    if isinstance(witnesses,dict):
        witnesses=[witnesses] if 'triangle' in witnesses else list(witnesses.values())
    for key in sorted(STAR):
        ref,pts,scale=_geometry(body,solid,key)
        geometry[key]=[]
        for half,(_,pc) in enumerate(_children(key,ref,pts,surplus,fraction)):
            lo,hi=(ZERO,fraction) if half==0 else (fraction,ONE)
            geometry[key].append((pc,_coefficients(pc)/(scale*(hi-lo)),lo,hi))
    def record(key,p,explicitly=False):
        b1,b2,y=map(float,p)
        if not (0<=b1<=1 and 0<=b2<=1 and 0<=y<=1):
            return
        tag=(*key,b1,b2,y)
        if tag in trials:
            if explicitly:
                explicit[tag]=True
            return
        b1,b2,y=map(Q,(b1,b2,y))
        bary=np.asarray((1-b1-b2,b1,b2),object)
        if bary[0]<0:
            return
        half=int(y>=fraction)
        pc,coef,lo,hi=geometry[key][half]
        u=(y-lo)/(hi-lo)
        point=bary@((1-u)*pc[0]+u*pc[1])
        chi=support.at(point)
        if chi<=0:
            return
        det=bary@coef@np.array([(1-u)**2,2*u*(1-u),u*u],object)
        if det<=0:
            raise ArithmeticError(f'Positive-order refined sample with nonpositive determinant at {tag}')
        trials[tag]=(float(chi),float(det))
        if explicitly:
            explicit[tag]=True
    def sample_segment(key,geom,p,q):
        nonlocal boundary_faces, finite_energy_boundary
        pc,coef,lo,hi=geom
        def spatial(param):
            b1,b2,y=param
            u=(y-lo)/(hi-lo)
            return np.array((1-b1-b2,b1,b2),object)@((1-u)*pc[0]+u*pc[1])
        a,b=spatial(p),spatial(q)
        cuts=_solve_knots(a,b,support.axes)
        for c in cuts:
            record(key,_parameters(c,p,q))
        # At a Q1-cell boundary the physically ordered side may have a tiny
        # tail. Interior positions along this segment have a single support
        # pattern, so a zero endpoint with positive midpoint is a witnessed
        # continuous crossing. Approach it with representable coordinates.
        for v,w in zip(cuts,cuts[1:]):
            m=(v+w)/2
            pos=lambda t: (1-t)*a+t*b
            if support.at(pos(m))<=0:
                continue
            record(key,_parameters(m,p,q))
            for edge in (v,w):
                if support.at(pos(edge))==0:
                    boundary_faces+=1
                    # Prove the weighted-energy infimum if this ordered
                    # segment approaches a regular, zero-order endpoint.
                    # J stays positive by continuity, so W stays bounded.
                    # For any J>0, AM-GM of four squared stretches gives
                    # W >= 2*(sqrt(J)-1-log(sqrt(J))) + (J-1)^2/4 >= 0.
                    if finite_energy_boundary is None:
                        b1,b2,y=_parameters(edge,p,q)
                        bary=np.array((1-b1-b2,b1,b2),object)
                        u=(y-lo)/(hi-lo)
                        chart_det=bary@coef@np.array([(1-u)**2,2*u*(1-u),u*u],object)
                        point=pos(edge)
                        if chart_det>0 and min(point[0],point[2],y)>0:
                            ref,_,_=_geometry(body,solid,key)
                            x,t,a=bary@((1-y)*ref[0]+y*ref[1])
                            from mouth_ordered_j_bounds import reference_factors
                            metric=reference_factors(body,x,t,key[1]<body.cut_node)
                            if metric and all(S>0 and d>0 for S,d in metric):
                                finite_energy_boundary=dict(cell_t_x_angle=list(key[:3]),triangle=key[3],
                                    barycentric=list(map(float,bary)),angular_fraction=float(y),
                                    exact_order=0., chart_J=float(chart_det),
                                    positive_midpoint_order=float(support.at(pos(m))),
                                    strict_J_inside=bool(chart_det>0),
                                    ordered_approach='One Q1 support interval with positive order in its interior')
                    for shift in (8,24,42,51):
                        near=edge+(m-edge)/2**shift
                        record(key,_parameters(near,p,q))
        if p[:2]==q[:2]:
            b1,b2=p[:2]
            c=np.array((1-b1-b2,b1,b2),object)@coef
            square=c[0]-2*c[1]+c[2]
            if square>0:
                y=(c[0]-c[1])/square
                if 0<y<1:
                    record(key,_parameters(y,p,q))
    for key,halves in geometry.items():
        for geom in halves:
            _,_,lo,hi=geom
            for bary in (*UNIT,*MID,CENTROID):
                b1,b2=bary[1:]
                sample_segment(key,geom,(b1,b2,lo),(b1,b2,hi))
            for y in (lo,(lo+hi)/2,hi):
                for i in range(3):
                    p,q=UNIT[i],UNIT[(i+1)%3]
                    sample_segment(key,geom,(*p[1:],y),(*q[1:],y))
                    sample_segment(key,geom,(*p[1:],y),(*CENTROID[1:],y))
    for w in witnesses:
        if not w or 'cell_t_x_angle' not in w:
            continue
        key=_from_witness(w)
        if key in STAR and 'barycentric' in w and 'angular_fraction' in w:
            record(key,(*w['barycentric'][1:],w['angular_fraction']),explicitly=True)
    # Rank cheaply using exact chart determinant and exact native Q1 order.
    # Each shortlist is then evaluated with the physical metric and full
    # constitutive energy. Rank every child separately, avoiding domination
    # of all candidate lists by a single nearly melted cell.
    per_child=defaultdict(list)
    for tag in trials:
        per_child[(*tag[:4],int(Q(tag[-1])>=fraction))].append(tag)
    selected=set(explicit)
    for tags in per_child.values():
        for rank in (lambda t:trials[t][0], lambda t:trials[t][1],
                     lambda t:trials[t][0]*(1+abs(log(trials[t][1])))):
            selected.update(sorted(tags,key=rank)[:3])
    readouts=[]
    for tag in sorted(selected):
        j,i,k,tri,b1,b2,y=tag
        w=child_witness(body,solid,surplus,order,dict(cell_t_x_angle=[j,i,k],triangle=tri,
                                barycentric=[1.-b1-b2,b1,b2], angular_fraction=y))
        if w['chi']<=0 or w['J_interval'][0]<=0 or w['W'] is None:
            raise ArithmeticError('Positive exact chart witness failed physical constitutive readback')
        readouts.append(w)
    if not readouts:
        return dict(sampled_positive_points=0, conclusion='No strictly ordered sample in the refined region',
                    certified_region_has_order=None, upper_bounds=None)
    jw=min(readouts,key=lambda w:w['J_interval'][1])
    cw=min(readouts,key=lambda w:w['chi'])
    ew=min(readouts,key=lambda w:w['chi_W'])
    return dict(scope='Six parents and both angular children. Lower J bound assumes the full orientation guard passed.',
                sampled_positive_points=len(trials), full_constitutive_points=len(readouts),
                observed_ordered_children=len(per_child), positive_to_exact_zero_segment_endpoints=boundary_faces,
                infimum_J_bounds=[0.,jw['J_interval'][1]], best_J_witness=jw,
                infimum_chi_bounds=[0.,0. if boundary_faces else float(np.nextafter(cw['chi'],np.inf))], smallest_positive_chi_witness=cw,
                infimum_chi_W_bounds=[0.,0. if finite_energy_boundary else None],
                regular_zero_order_boundary=finite_energy_boundary, chi_W_witness_estimate=ew['chi_W'], smallest_chi_W_witness=ew,
                chi_W_qualification='Upper witness evaluated with the constitutive law to float precision; no certified rounding bound.',
                zero_order_infimum_proved_on_segment=bool(boundary_faces),
                order_values='Distinct positive witnesses. An unresolved lower infimum need not be attained at a point.')
