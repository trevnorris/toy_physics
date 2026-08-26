#!/usr/bin/env python3
"""FORM ABLATION 1 — freeze traced-perturbation evaluation at flat face.

Change: in physical_trace_fields, pass h_eval = face*W0/2 into compose_physical_trace
instead of h0 = face*W_bg/2. Geometry (normals, dh, measure) unchanged.

Re-emit FACE_SHIFT pressure + velocity_w and RELATIVE_FLUX for one case.
Report literal before/after diff.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import sympy as sp

ROOT = Path("/tmp/s11c_a_review_grok")
BASE = ROOT / "engine_copy.py"
ABLATED = ROOT / "engine_ablate_eval_face.py"


def patch_eval_face(src: str) -> str:
    """Replace physical_trace_fields so composition uses flat face, not h0=s W_bg/2."""
    old = '''def physical_trace_fields(
    face: int,
    h0: sp.Expr,
    dh: sp.Expr,
    parameter: sp.Symbol,
    rho4_bg_exact: sp.Expr,
) -> tuple[sp.Expr, tuple[sp.Expr, ...], sp.Expr, tuple[sp.Expr, ...]]:
    """Build all physical traces from the common §3c field composition."""
    pressure_reference = delta_p_plus if face == 1 else delta_p_minus
    pressure_perturbation = affine_bulk_perturbation(
        pressure_reference, dw_delta_p[face], face,
    )
    velocity_perturbation = tuple(
        affine_bulk_perturbation(
            delta_v_bulk[face][i], dw_delta_v_bulk[face][i], face,
        )
        for i in range(4)
    )
    density_perturbation = affine_bulk_perturbation(
        delta_rho4_face[face], dw_delta_rho4_face[face], face,
    )
    current_perturbation = tuple(
        affine_bulk_perturbation(j_bulk[i], dw_delta_j_bulk[face][i], face)
        for i in range(4)
    )

    # These are the physical members of B0 supplied in §§2d, 3b, and 3c.
    # Keeping their construction here makes every normal jet enter through an
    # actual d/dw performed by the trace composition, including the current
    # background rho_4D^0 v_bulk^0 and each density representative.
    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    density_background = rho4_bg_exact.subs(parameter, 0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )

    pressure_trace = compose_physical_trace(
        pressure_background, pressure_perturbation, h0, dh, parameter,
    )
    velocity_trace = tuple(
        compose_physical_trace(
            velocity_background[i], velocity_perturbation[i], h0, dh, parameter,
        )
        for i in range(4)
    )
    density_trace = compose_physical_trace(
        density_background, density_perturbation, h0, dh, parameter,
    )
    current_trace = tuple(
        compose_physical_trace(
            current_background[i], current_perturbation[i], h0, dh, parameter,
        )
        for i in range(4)
    )
    return pressure_trace, velocity_trace, density_trace, current_trace'''

    new = '''def physical_trace_fields(
    face: int,
    h0: sp.Expr,
    dh: sp.Expr,
    parameter: sp.Symbol,
    rho4_bg_exact: sp.Expr,
) -> tuple[sp.Expr, tuple[sp.Expr, ...], sp.Expr, tuple[sp.Expr, ...]]:
    """Build all physical traces from the common §3c field composition."""
    # FORM ABLATION: freeze evaluation face at flat s*W0/2 instead of h0=s*W_bg/2.
    h_eval = sp.Rational(face, 2) * W0
    pressure_reference = delta_p_plus if face == 1 else delta_p_minus
    pressure_perturbation = affine_bulk_perturbation(
        pressure_reference, dw_delta_p[face], face,
    )
    velocity_perturbation = tuple(
        affine_bulk_perturbation(
            delta_v_bulk[face][i], dw_delta_v_bulk[face][i], face,
        )
        for i in range(4)
    )
    density_perturbation = affine_bulk_perturbation(
        delta_rho4_face[face], dw_delta_rho4_face[face], face,
    )
    current_perturbation = tuple(
        affine_bulk_perturbation(j_bulk[i], dw_delta_j_bulk[face][i], face)
        for i in range(4)
    )

    # These are the physical members of B0 supplied in §§2d, 3b, and 3c.
    # Keeping their construction here makes every normal jet enter through an
    # actual d/dw performed by the trace composition, including the current
    # background rho_4D^0 v_bulk^0 and each density representative.
    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    density_background = rho4_bg_exact.subs(parameter, 0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )

    pressure_trace = compose_physical_trace(
        pressure_background, pressure_perturbation, h_eval, dh, parameter,
    )
    velocity_trace = tuple(
        compose_physical_trace(
            velocity_background[i], velocity_perturbation[i], h_eval, dh, parameter,
        )
        for i in range(4)
    )
    density_trace = compose_physical_trace(
        density_background, density_perturbation, h_eval, dh, parameter,
    )
    current_trace = tuple(
        compose_physical_trace(
            current_background[i], current_perturbation[i], h_eval, dh, parameter,
        )
        for i in range(4)
    )
    return pressure_trace, velocity_trace, density_trace, current_trace'''

    if old not in src:
        raise SystemExit("ABORT: target physical_trace_fields block not found exactly")
    return src.replace(old, new, 1)


def emit_objects(mod):
    # Clear face cache
    mod._FACE_CACHE.clear()
    mod._FINAL_CACHE.clear()
    src = mod.build_face_source("LAB_HELD", 1, "ZETA_C", "RHO4_CONSTANT")
    pressure = mod.finalize(mod.pressure_trace_raw(src))
    velocity = mod.finalize(
        mod.epsilon * sp.ImmutableMatrix(
            mod.shape(src.bulk_velocity_trace_exact, src.parameter)
        )
    )
    density = mod.finalize(mod.epsilon * mod.shape(src.density_trace_exact, src.parameter))
    flux = mod.finalize(mod.relative_flux_raw(src))
    return {
        "FACE_SHIFT.pressure": pressure,
        "FACE_SHIFT.velocity_w": velocity[3],
        "FACE_SHIFT.density": density,
        "RELATIVE_FLUX": flux,
        "h0_geometry": src.h0,
    }


def main():
    src = BASE.read_text()
    ABLATED.write_text(patch_eval_face(src))
    print("Wrote", ABLATED)
    print("FORM CHANGE: compose_physical_trace evaluation face h0=s*W_bg/2 -> h_eval=s*W0/2")
    print()

    sys.path.insert(0, str(ROOT))
    # Import pristine under a clean name by loading from file via runpy-style
    import importlib

    # Load baseline
    if "engine_copy" in sys.modules:
        del sys.modules["engine_copy"]
    base = importlib.import_module("engine_copy")
    before = emit_objects(base)

    # Load ablated under distinct module name
    import importlib.util
    spec = importlib.util.spec_from_file_location("engine_ablate_eval_face", ABLATED)
    abl = importlib.util.module_from_spec(spec)
    sys.modules["engine_ablate_eval_face"] = abl
    spec.loader.exec_module(abl)
    after = emit_objects(abl)

    print("=" * 72)
    print("LITERAL BEFORE / AFTER / DIFF")
    print("=" * 72)
    any_moved = False
    for key in before:
        b, a = before[key], after[key]
        if isinstance(b, sp.Basic) and isinstance(a, sp.Basic):
            diff = sp.simplify(sp.expand(a - b))
            moved = diff != 0
        else:
            diff = f"(non-expr) before={b!r} after={a!r}"
            moved = b != a
        any_moved = any_moved or moved
        print(f"\n--- {key} ---")
        print(f"BEFORE: {b}")
        print(f"AFTER:  {a}")
        print(f"DIFF:   {diff}")
        print(f"MOVED:  {moved}")

    print("\n" + "=" * 72)
    print(f"ABLATION_MOVED_ANY_OBJECT = {any_moved}")
    print(
        "JUDGEMENT: evaluation face of traced perturbation is LOAD-BEARING iff MOVED.\n"
        "If AFTER loses the eta_bg*(W0/2)*w1*d_w_delta_* piece while geometry h0 still\n"
        "carries W_bg, the physics of §3c (evaluate df at h_s0, not freeze at flat) is\n"
        "confirmed present in the pristine engine."
    )


if __name__ == "__main__":
    main()
