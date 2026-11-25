"""Signed physical/gauge derivatives and an off-centre melted-edge control at held R4 step 166."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_cp10_spatial_crystal import MODES
from mouth_r4_aligned_crystal import AlignedCrystal, cell_call
from mouth_r4_aligned_derivatives import EPS, label_bank, order_bank, shape_bank
from mouth_r4_aligned_hessian import ReferenceGauge, audit_site, edge_directions
from mouth_r4_aligned_stencil import BASIS, GRID, VERTICES, samples
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, chart_orientation
from mouth_r4_eulerian_solid import InvalidSpatialCrystal
from mouth_r4_newton_audit import (FULL_ERROR_TOL, HALVING_TOL, SYMMETRY_TOL,
                                  WIDTHS, PromotedGradient, directions)
from mouth_r4_smooth_static import LD
from mouth_r4_stage2_outcome import digest


def strict_controls(joint,labels):
    """Isolate an off-centre negative determinant in one still-melted cell on step 166.

    Follow a free radial label incident on a shallow edge dip. At fixed sample
    coordinates this rank-one perturbation is affine in its size, so the
    smallest positive 27-point root bounds the search before any vertex can
    become negative.
    """
    order=joint.order
    axes=joint.opt.axes
    cert=chart_orientation(axes,labels,joint.opt.active)
    m=joint.solid.evaluate(labels,order,gradient=False)[3]
    pos=dict(whole_min=cert['chart_det_global_lower'],crystal_min=m['min_chart_Jinv_bound'],
             whole_active_min=cert['chart_det_active_lower'],crystal_active_min=m['min_chart_Jinv_active_bound'])
    pos['passed']=(pos['whole_min']==pos['crystal_min']>0 and
                   pos['whole_active_min']==pos['crystal_active_min']>0)
    ids=joint.weighted.ids
    dims=tuple(np.array(order.shape)-1)
    width=joint.weighted.width
    points=labels.reshape(-1,3)[ids].astype(LD)
    fixed,_,edges=samples(points,width)
    fully_melted=np.max(order.ravel()[ids],axis=1)==0
    # Rank by depth of the off-centre dip relative to all fixed samples.
    ranks=np.argsort(edges.min(axis=1)/fixed.min(axis=1))
    cm=joint.components[0]
    for idx in ranks:
        if not fully_melted[idx]:continue
        xyz=np.asarray(np.unravel_index(idx,dims))
        for node_index,base_id in enumerate(ids[idx]):
            group=cm.groups[base_id]
            if not cm.dead[group] or np.count_nonzero(cm.groups==group)!=1:continue
            dy=np.zeros((8,3),LD)
            dy[node_index,0]=1
            d=np.einsum('vi,qvj->qij',dy,BASIS)/width[idx,None,None,:]
            matrix=(points[idx]-points[idx,:1]).T@np.asarray(BASIS)/width[idx,None,None,:]
            # Only the radial row changes. Its exact determinant derivative is
            # the corresponding cofactor contraction.
            sens=np.sum(d[:,0,:]*np.cross(matrix[:,1,:],matrix[:,2,:]),axis=1)
            for sign in (-1,1):
                delta=sign*sens
                falling=delta<0
                if not np.any(falling):continue
                wall=np.min(-fixed[idx][falling]/delta[falling])
                for frac in (.2,.55,.85,.97,.999,.99999):
                    shift=float(sign*frac*wall)
                    y=labels.copy()
                    node=tuple(map(int,xyz+VERTICES[node_index]))
                    y[node+(0,)]+=shift
                    p=y.reshape(-1,3)[ids[idx:idx+1]].astype(LD)
                    f,t,j=samples(p,width[idx:idx+1])
                    if f.min()<=0 or j.min()>=-1e-9:continue
                    edge=int(np.argmin(j))
                    # A perturbation on a shared vertex may invert a different
                    # adjacent element first. Compare the two complete traversals.
                    chart=physical=barrier=None
                    try:chart_orientation(axes,y,joint.opt.active)
                    except InfeasibleLabels as exc:chart=exc.witness
                    try:joint.solid.evaluate(y,order,gradient=False)
                    except InvalidSpatialCrystal as exc:physical=exc.diagnostics
                    try:cell_call('values',p,width[idx:idx+1])
                    except InfeasibleLabels as exc:barrier=exc.witness
                    if not (chart and physical and barrier):continue
                    if not(chart.get('code')==physical.get('code')==2 and
                           chart['cell']==physical['cell'] and chart['chart_fractions']==physical['local']
                           and chart['chart_lower']==physical['chart_bound'] and chart['chart_det']<0
                           and not joint.opt.active[tuple(chart['cell'])]
                           and barrier.get('stencil_index',-1)>=27):continue
                    negative=dict(cell=xyz.tolist(),node=node,shift=shift,
                        edges_affected=int((j<0).sum()),minimum_edge_index=edge,
                        fraction=float(t[0,edge]),edge_det=float(j[0,edge]),
                        fixed_min=float(f.min()),vertex_min=float(f[0,np.all(GRID!=1,axis=1)].min()),
                        local_order_max=float(order.ravel()[ids[idx]].max()),
                        chart=chart,physical=physical,stencil=barrier,passed=True)
                    return dict(passed=pos['passed'],positive=pos,negative=negative)
    raise ArithmeticError('Could not exhibit an inverted melted edge with all 27 fixed samples positive')


def hessian(driver,t,solid,order,out):
    joint=AlignedCrystal(driver,t,solid,driver.labels,order)
    oracle=ReferenceGauge(t,driver.reference_depth,order)
    promoted=PromotedGradient(t.b,out)
    edge,eds=edge_directions(joint,driver.labels)
    cases=[(name,spec,dirs,WIDTHS) for name,spec,dirs in directions(t,joint.opt,order)]
    cases.append(('melted_edge',edge,eds,(2.**-9,2.**-10,2.**-11)))
    sites={name:audit_site(joint,driver,promoted,oracle,spec,dirs,widths)
           for name,spec,dirs,widths in cases}
    return dict(sites=sites,thresholds=dict(relative_Hv=FULL_ERROR_TOL,
                reference_halving=HALVING_TOL,symmetry=SYMMETRY_TOL),
                passed=all(row['passed'] for row in sites.values()))


def run(root,case,sector,out):
    from mouth_r4_collar166_state import load_held
    out.mkdir(exist_ok=False,parents=True)
    d,t,solid,order=load_held(root,out,EPS[case],smooth=True)
    if sector.startswith('label'):
        result=label_bank(d,t,solid,order,out,int(sector[-1]))
    elif sector=='order':
        result=order_bank(d,t,solid,order,out)
    elif sector=='shape':
        result=shape_bank(d,t,solid,order)
        if [r['mode'] for r in result['rows']]!=MODES:
            raise ArithmeticError('Incomplete six-force bank')
        wide=next(b for b in result['rows'][0]['central'] if b['h']==2.**-18)
        result.update(wide_depth=wide)
        result['passed']=bool(result['passed'] and abs(wide['minus_analytic'])<=result['threshold'])
    elif sector=='orientation':
        result=strict_controls(AlignedCrystal(d,t,solid,d.labels,order),d.labels)
    else:
        result=hessian(d,t,solid,order,out)
    result.update(epsilon=EPS[case],sector=sector,case=case,step=166,state_sha256=(
        '4b84a0d2ad53889a1d619e68960df70233ceb589fc4fe6c336c6e59e832a2694'),
        source_sha256=digest(Path(__file__)),held_shape=d.x.tolist(),mantissa_bits=np.finfo(LD).nmant+1)
    (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(epsilon=EPS[case],sector=sector,passed=result['passed'])),flush=True)
    if not result['passed']:raise ArithmeticError('Frozen smooth-reference audit failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('study',type=Path)
    p.add_argument('case',choices=EPS)
    p.add_argument('sector',choices=('label0','label1','label2','order','shape','hessian','orientation'))
    p.add_argument('output',type=Path)
    a=p.parse_args()
    run(a.study.resolve(),a.case,a.sector,a.output.resolve())
