#!/usr/bin/env python3
"""C7 faithfulness mutator for the stage030 -> stage031 ``zero_mode`` edge.

The C7 harness (``composite_build.py:check_c7``) invokes this script with
``env[C7_FACET] = exported_facet`` and reads the LAST stdout line as a JSON
object ``{consumer_stage: tooth_or_"PASS"}``. Its job is to prove that
stage031's zero-mode-consuming tooth is a GENUINE dependency, not a decorative
one: if we perturb stage030's transverse zero mode ``f0`` so that it is no
longer the true kernel of the transverse operator ``O_perp``, then stage031's
``PASS_F0_MOUTH_VALUE_EVALUATED`` tooth must fire.

HONESTY CONTRACT (non-negotiable). For the zero-mode facet this script does NOT
look up an answer. It ACTUALLY perturbs the kernel exponent (2 -> 1,
``f0 = 1/(ell*cosh(w/ell))``) and COMPUTES the transverse-operator residual

    O_perp f0 = -d^2 f0 / dw^2 + V_H * f0 ,   V_H = (4 - 6*sech(w/ell)^2)/ell^2

with sympy (V_H and the residual form are copied verbatim from the stage030
audit ``ledger_stage030_electric_scalar_localized_h_closure_sympy_audit.py``,
locus ``PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE``). The consumer tooth is emitted
ONLY when that computed residual simplifies to something nonzero. The genuine
kernel (power 2) gives residual 0 -> PASS; the perturbation (power 1) gives a
nonzero residual -> the zero-mode property is genuinely broken -> the tooth
fires. The tooth string is therefore derived from the computed outcome, never
printed unconditionally.

``--decorative`` reproduces the dishonest contrast: it IGNORES the facet and the
computation and always prints PASS, which the checker must catch as
DECORATIVE_DEPENDENCY.
"""

from __future__ import annotations

import argparse
import json
import os

import sympy as sp

# The exact facet string this producer export declares (== exported_facet ==
# c7_expect.facet_used). Any other facet is not ours -> honest mutators ignore
# unknown facets and stay silent (PASS).
ZERO_MODE_FACET = "zero_mode_kernel_power"

# The stage031 tooth whose assertion consumes stage030's zero mode. It re-
# evaluates f0(0)=1/ell and protects the sech^2 shape (curvature at w=0); it is
# an existing tooth in stage031.verification.teeth (mutation "f0_power",
# claim f0_mouth_value). Emitted ONLY when the perturbation is computably real.
CONSUMER_STAGE = "stage031"
CONSUMER_TOOTH = "PASS_F0_MOUTH_VALUE_EVALUATED"

# The genuine kernel exponent; the perturbation drops it to 1.
CANONICAL_POWER = 2
PERTURBED_POWER = 1


def transverse_operator_residual(power: int) -> sp.Expr:
    """Return simplify(O_perp f0) with f0 = 1/(ell*cosh(w/ell)**power).

    O_perp = -d^2/dw^2 + V_H,  V_H = (4 - 6*sech(w/ell)^2)/ell^2  (stage030).
    Zero for the true kernel (power 2); nonzero otherwise.
    """
    w, ell = sp.symbols("w ell", positive=True, real=True)
    V_H = (4 - 6 / sp.cosh(w / ell) ** 2) / ell**2
    f0 = 1 / (ell * sp.cosh(w / ell) ** power)
    return sp.simplify(-sp.diff(f0, w, 2) + V_H * f0)


def compute_zero_mode_outcome() -> tuple[str, sp.Expr]:
    """Honestly compute the consumer outcome for the zero-mode facet.

    Perturb the kernel exponent 2 -> 1 and test whether the perturbed f0 is
    still annihilated by O_perp. Returns (outcome, residual) where outcome is
    the tooth name iff the residual is computably nonzero, else "PASS".
    """
    residual = transverse_operator_residual(PERTURBED_POWER)
    if residual != 0:
        # The perturbed profile is NOT the zero mode -> f0(0)=1/ell no longer
        # follows from the true kernel; stage031's tooth genuinely fails.
        return CONSUMER_TOOTH, residual
    return "PASS", residual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--decorative",
        action="store_true",
        help="Dishonest contrast: ignore the facet and always print PASS.",
    )
    parser.add_argument(
        "--show-residual",
        action="store_true",
        help="Also print the computed residual (both powers) to stderr for audit.",
    )
    args = parser.parse_args()

    facet = os.environ.get("C7_FACET", "")

    if args.decorative:
        # Decorative dependency: no computation, no facet sensitivity.
        print(json.dumps({CONSUMER_STAGE: "PASS"}, sort_keys=True))
        return

    if facet == ZERO_MODE_FACET:
        outcome, residual = compute_zero_mode_outcome()
        if args.show_residual:
            canonical = transverse_operator_residual(CANONICAL_POWER)
            print(
                f"[audit] V_H = (4 - 6*sech(w/ell)**2)/ell**2\n"
                f"[audit] O_perp f0 (power {CANONICAL_POWER}, true kernel) = {canonical}\n"
                f"[audit] O_perp f0 (power {PERTURBED_POWER}, perturbed)   = {sp.simplify(residual)}\n"
                f"[audit] residual nonzero -> {residual != 0} -> outcome {outcome!r}",
                flush=True,
            )
        out = {CONSUMER_STAGE: outcome}
    else:
        # Unknown / non-owned facet: an honest mutator does not warn on facets
        # it does not model.
        out = {CONSUMER_STAGE: "PASS"}

    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()
