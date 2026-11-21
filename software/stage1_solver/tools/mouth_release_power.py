"""Exact release-force and signed power identities for the two-support comparison."""
import sympy as sp


def main():
    incoming, rim, support, share = sp.symbols("U V W h", real=True)
    angle = sp.symbols("theta", real=True)
    cosine, sine = sp.cos(angle), sp.sin(angle)
    decision_v = rim + (incoming-rim)*cosine
    decision_stretch = (rim-incoming)*sine
    anchor = (1-share)*rim + share*support
    # Finish the portion of the spring orbit remaining to a full period.
    final_v = sp.trigsimp(anchor + (decision_v-anchor)*cosine - decision_stretch*sine)
    final_stretch = sp.trigsimp(decision_stretch*cosine + (decision_v-anchor)*sine)
    recovered = sp.trigsimp((final_v-incoming)/(1-cosine))
    second = final_v-decision_v
    rim_impulse = decision_v-incoming + (1-share)*second
    other_impulse = share*second
    energy_error = sp.trigsimp(rim*rim_impulse + support*other_impulse
                               - (final_v**2-incoming**2+final_stretch**2)/2)
    A, B = sp.symbols("A B", nonnegative=True)
    rim_mean = sp.expand((A-B)*(rim-incoming) + B*sp.simplify(rim_impulse/(1-cosine)))
    c = -sp.diff(rim_mean, incoming)
    g = sp.diff(rim_mean, rim)
    gap = sp.factor(A*g-c**2)
    print("returned_cohort_exit_change_over_1_minus_cos =", sp.factor(recovered))
    print("returned_cohort_energy_residual =", energy_error)
    print("held_coefficient =", c, "tow_coefficient =", g, "A*g-c^2 =", gap)
    assert sp.simplify(recovered-share*(support-rim)) == 0
    assert energy_error == 0 and sp.expand(gap-B*(A-B)*(1-share)**2) == 0
    assert sp.expand(rim_mean - c*(support-incoming) - g*(rim-support)) == 0

    # Independent of the spring model. One cold incident sideways stream, a
    # stationary support, no incremental drain power and no depleting storage.
    held, drag = sp.symbols("c g", real=True)
    available = incoming**2/2 + rim*(drag*rim-held*incoming)
    equivalent = (incoming-held*rim)**2/2 + (drag-held**2/2)*rim**2
    contrast = (available - available.subs(rim,-rim) - available.subs(incoming,-incoming)
                + available.subs({rim:-rim,incoming:-incoming}))/4
    ideal_deficit = -available.subs({held:1, drag:0, rim:2*incoming})
    print("minimum_available_at_U_equals_cV =", sp.factor(available.subs(incoming,held*rim)))
    print("four_signed_runs_cross_power_per_mass_rate =", sp.factor(contrast))
    print("ideal_branch_deficit_at_V_equals_2U_per_mass_rate =", ideal_deficit)
    assert sp.expand(available-equivalent) == 0
    assert sp.expand(contrast + held*incoming*rim) == 0
    assert ideal_deficit == 3*incoming**2/2


if __name__ == "__main__":
    main()
