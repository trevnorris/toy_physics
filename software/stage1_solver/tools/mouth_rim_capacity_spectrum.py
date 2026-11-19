"""Predict the mean site energy from independent one- and two-quantum data.

One isolated receiver mode, with an onsite energy increment U for storing
the second quantum. The passing parcel has two internal states. The
exchange matrix elements are G at the first rung and sqrt(2) G at the next.
"""
import numpy as np
import sympy as sp


def weak_result():
    detuning, second_increment, g = sp.symbols("Delta U G", positive=True)
    first = (sp.sqrt(detuning**2 + 4*g**2)-detuning)/2
    second = (sp.sqrt((detuning+second_increment)**2 + 8*g**2)-detuning-second_increment)/2
    mean = first - second/2
    weak = sp.factor(mean.diff(g,2).subs(g,0)/2)
    harmonic = sp.factor(mean.subs(second_increment,0).diff(g,4).subs(g,0)/24)
    print(f"weak_centroid_over_G2={weak} harmonic_centroid_over_G4={harmonic}")
    assert sp.simplify(weak - second_increment / (detuning*(detuning+second_increment))) == 0


def measured_split(detuning, g):
    matrix = np.array([[0.,g],[g,detuning]])
    levels = np.linalg.eigvalsh(matrix)
    return levels[1] - detuning if detuning > 0 else levels[0] - detuning


if __name__ == "__main__":
    weak_result()
    for delta, u in ((2,0), (2,1), (2,2), (2,18), (-2,1), (-2,4)):
        # This perturbative test excludes exact crossing at Delta+U=0.
        for g in (.1, .05):
            first = measured_split(delta, g)
            second = measured_split(delta+u, np.sqrt(2)*g)
            spin_e, spin_g = first - second, first
            mean = (spin_e+spin_g)/2
            weak = u/(delta*(delta+u))
            print(f"levels delta={delta} U={u} G={g} shift_g={spin_g:.12g} shift_e={spin_e:.12g} "
                  f"common={mean:.12g} common_over_G2={mean/g**2:.12g} "
                  f"predicted_weak={weak:.12g} ratio_to_single_capacity={mean/first:.12g}")
            assert abs(mean / g**2 - weak) < 0.06
