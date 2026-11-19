"""Bound quasistatic work by a converting stream on the compact mouth links.

The original incident ensemble has three equal nonzero eigenvalues 1/3.
Each marked twist derivative has eigenvalues t,t,0,0,-t,-t. In the gate
whose site-scalar energies agree, each impulse has zero full-turn integral.
A version of the first-slip bound does not need the zero-period symmetry.
"""
import argparse
from pathlib import Path

from mpmath import iv

from mouth_collar_link_bounds import R,PI,TAU,mesh,dot,cross,spring,rounded_bound
from mouth_standing_wave_collar_cost import rows


def area_weights(level):
    verts,faces=mesh(level)
    w=[]
    for i,j,k in faces:
        a,b,c=verts[i],verts[j],verts[k]
        w.append(2*iv.atan2(dot(a,cross(b,c)),1+dot(a,b)+dot(b,c)+dot(c,a))/(4*PI))
    assert sum(w).a<=1<=sum(w).b
    return w


def certify(path):
    entries=list(rows(path,"forest"))
    level=int(entries[0]["level"])
    w=area_weights(level)
    m2=sum(x**2 for x in w)
    wmax=max(x.b for x in w)
    nf=len(w)
    # Unit working scale: every arrival visits one local loop; mean dwell 1,
    # hopping t=1, total number rate 1/4. No free intake fit across radii.
    A=R(2)/3*R(".25")
    c=spring(PI,4)/PI**2
    print(f"geometry level={level} area_weight_sq_sum_upper={rounded_bound(m2,10,True):.10f} "
          f"largest_weight_upper={rounded_bound(wmax,9,True):.9f} "
          f"max_force_sum_upper={rounded_bound(A,8,True):.8f}")
    for d in entries:
        assert d["level"]==str(level)
        scale=R(d["R/R0"])
        K=1/scale
        exterior=PI/(32*scale)
        Cmin=K*min(spring(PI,4).a,(5*spring(TAU/5,4)).a,(7*nf*spring(TAU/nf,4)).a)
        B=A**2*m2/(4*K*c)
        assert B<Cmin
        modal=R(d["all_occupations_gain_upper"])
        # Allow all principal surface angles to vary. Since |Wparcel|<=
        # 2 sqrt(B E_collar), minimize E_collar-2 sqrt(B E_collar) on [Cmin,∞).
        net=exterior+Cmin-2*iv.sqrt(B*Cmin)-modal
        # Before the first half-turn, all face lifts are principal. Pay for
        # that face; use Young on every other driven face and drop exterior.
        # The previous modal table covers any field at fixed sorted actions.
        before_slip=K*spring(PI,4)-PI*A*wmax-B-modal
        critical_A=(exterior+Cmin-modal)*iv.sqrt(K*c/(m2*Cmin))
        # This local threshold follows from the square-root bound only while
        # B<Cmin; the first obstacle is whichever comparison is lost first.
        threshold=min(critical_A.a,(2*iv.sqrt(K*c*Cmin/m2)).a)
        print(f"combined level={level} R/R0={d['R/R0']} wave_upper={d['all_occupations_gain_upper']} "
              f"cost_after_converter_lower={rounded_bound(net+modal,8):.8f} "
              f"unpaid_energy_lower={rounded_bound(net,8):.8f} "
              f"first_slip_shortfall_lower={rounded_bound(before_slip,8):.8f} "
              f"exclusion_if_D_times_tT_below={rounded_bound(3*threshold/2,6):.6f}",flush=True)
        assert net>0 and before_slip>0


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("forests",nargs="+",type=Path)
    for evidence in parser.parse_args().forests:
        certify(evidence)
