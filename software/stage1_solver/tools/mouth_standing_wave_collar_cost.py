"""Combine the rounded modal upper bounds with the separately known field cost."""
import argparse
from pathlib import Path

import sympy as sp
from mpmath import iv

from mouth_collar_link_bounds import R, PI, TAU, spring, rounded_bound


def rows(path,prefix):
    for line in path.read_text().splitlines():
        if line.startswith(prefix+" "):
            yield dict(item.split("=",1) for item in line.split()[1:])


def main(ground,forests):
    # Ground and excited fractions of the axial wave load are additive.
    # Differentiation is at fixed modal action and independent L and R.
    L,nu,I,k=sp.symbols("L nu I k",positive=True)
    omega=sp.sqrt(nu**2+(k/L)**2)
    axial=-sp.diff(I*omega,L)
    print("symbolic axial force times L^3 omega / action =",sp.factor(axial*L**3*omega/I))
    assert sp.simplify(axial*L**3*omega/I-k**2)==0
    ground_rows={(d["level"],d["R/R0"]):d for d in rows(ground,"budget") if d["eta"]=="+1"}
    for path in forests:
        for d in rows(path,"forest"):
            lev=int(d["level"])
            ref=d["R/R0"]
            nv_faces=8*4**lev
            gap=min(spring(PI,4).a,(5*spring(TAU/5,4)).a,
                    (7*nv_faces*spring(TAU/nv_faces,4)).a)
            # Values read from earlier stdout are conservative in the printed
            # direction. The compact cost is added only when the exterior has
            # nonzero net flux. Rb=16*(R/R0), C=1, K0=1 in the scaled example.
            upper=R(d["all_occupations_gain_upper"])
            ext=R(d["exterior_cost_lower"])
            g=R(ground_rows[(str(lev),ref)]["ground_saving_upper"])
            s=R(ref)
            k_req=max(R(0).a,((upper-ext)/gap).b)
            f_req=max(R(0).a,((ext-g)/(upper-g)).a)
            available=ext+gap/s
            print(f"cost level={lev} R/R0={ref} total_release_upper={rounded_bound(upper,8,True):.8f} "
                  f"scaled_collar_plus_exterior_lower={rounded_bound(available,8):.8f} "
                  f"sufficient_Kmin_above={rounded_bound(k_req,6,True):.6f} "
                  f"necessary_excited_load_fraction_exterior_only={rounded_bound(f_req,6):.6f}")
            assert upper < available
            if upper>=ext:
                assert 0<f_req<1


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("ground",type=Path)
    parser.add_argument("forests",nargs="+",type=Path)
    args=parser.parse_args()
    main(args.ground,args.forests)
