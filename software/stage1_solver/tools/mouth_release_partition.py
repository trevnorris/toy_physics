"""Determine which impulse measurements resolve a single 0PN mouth release."""

import sympy as sp


def main():
    # All impulses are signed in the lab and relative to the undriven control.
    # J is the applied impulse plus the feeding-side momentum input.
    stored, through, other = sp.symbols("dP I_through I_other")
    known, through_obs, other_obs, stored_obs = sp.symbols("J Q B P")
    unknowns = (stored, through, other)
    balance = stored + through + other - known

    # Even with the change in geometric velocity measured, an exchange
    # proportional to acceleration is invisible in the integrated sum.
    k, dv = sp.symbols("k delta_V")
    shifted = balance.subs({stored: stored + k*dv, through: through - k*dv},
                           simultaneous=True)
    print(f"budget change for opposite shifts at fixed delta_V: {sp.expand(shifted-balance)}")

    cases = [
        ("force and incoming", []),
        ("also other impulse", [other - other_obs]),
        ("also through impulse only", [through - through_obs]),
        ("other and through", [other - other_obs, through - through_obs]),
        ("other and assigned inventory", [other - other_obs, stored - stored_obs]),
    ]
    for label, readings in cases:
        matrix, rhs = sp.linear_eq_to_matrix([balance, *readings], unknowns)
        solution, = sp.linsolve((matrix, rhs), unknowns)
        free_directions = [tuple(vector) for vector in matrix.nullspace()]
        print(f"{label}: (dP,I_through,I_other)={solution} free_directions={free_directions}",
              flush=True)
        residual = matrix * sp.Matrix(solution) - rhs
        if residual.applyfunc(sp.expand) != sp.zeros(matrix.rows, 1):
            raise ArithmeticError(f"Impulse reconstruction failed for {label}")


if __name__ == "__main__":
    main()
