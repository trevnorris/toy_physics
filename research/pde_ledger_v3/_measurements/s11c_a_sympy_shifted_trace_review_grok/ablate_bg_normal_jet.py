#!/usr/bin/env python3
"""FORM ABLATION 2 — give a supplied background a nonzero w-dependence.

Change: density_background from the supplied w-independent rho4_bg to
  rho4_bg * (1 + (w - h0)/W0)
so ∂_w f⁰ |_{h0} = rho4_bg/W0 ≠ 0. That feeds the §3c jet term δh·∂_w f⁰.

Re-emit FACE_SHIFT.density (and pressure as a negative control: bg still 0).
Report literal before/after diff.
"""
from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path

import sympy as sp

ROOT = Path("/tmp/s11c_a_review_grok")
BASE = ROOT / "engine_copy.py"
ABLATED = ROOT / "engine_ablate_bg_jet.py"


def patch_bg_jet(src: str) -> str:
    old = '''    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    density_background = rho4_bg_exact.subs(parameter, 0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )'''

    new = '''    velocity_background = (sp.Integer(0),) * 4
    pressure_background = sp.Integer(0)
    # FORM ABLATION: replace supplied w-independent density background by a
    # w-dependent form so ∂_w f⁰(x,h_s⁰) is nonzero and feeds δh·∂_w f⁰.
    rho4_bg0 = rho4_bg_exact.subs(parameter, 0)
    density_background = rho4_bg0 * (1 + (w - h0) / W0)
    current_background = tuple(
        density_background * component for component in velocity_background
    )'''

    if old not in src:
        raise SystemExit("ABORT: target background block not found exactly")
    return src.replace(old, new, 1)


def emit_objects(mod):
    mod._FACE_CACHE.clear()
    mod._FINAL_CACHE.clear()
    src = mod.build_face_source("LAB_HELD", 1, "ZETA_C", "RHO4_CONSTANT")
    pressure = mod.finalize(mod.pressure_trace_raw(src))
    density = mod.finalize(mod.epsilon * mod.shape(src.density_trace_exact, src.parameter))
    # Also show the raw (pre-finalize) density linearisation to expose dh*dw_f0
    rho_lin = sp.diff(src.density_trace_exact, src.parameter).subs(src.parameter, 0)
    # Background jet that compose should have used
    # dens_bg = rho4*(1+(w-h0)/W0); d/dw at h0 = rho4/W0
    return {
        "FACE_SHIFT.pressure": pressure,
        "FACE_SHIFT.density": density,
        "density_lin_raw": rho_lin,
        "dh": src.dh,
        "h0": src.h0,
    }


def main():
    ABLATED.write_text(patch_bg_jet(BASE.read_text()))
    print("Wrote", ABLATED)
    print("FORM CHANGE: density_background := rho4_bg0 * (1 + (w-h0)/W0)")
    print("  => ∂_w density⁰ |_{h0} = rho4_bg0/W0  (was 0)")
    print()

    sys.path.insert(0, str(ROOT))
    if "engine_copy" in sys.modules:
        # re-import fresh
        pass
    base = importlib.import_module("engine_copy")
    # force reload in case prior ablation left state — engine_copy file is pristine
    base = importlib.reload(base)
    before = emit_objects(base)

    spec = importlib.util.spec_from_file_location("engine_ablate_bg_jet", ABLATED)
    abl = importlib.util.module_from_spec(spec)
    sys.modules["engine_ablate_bg_jet"] = abl
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
            diff = f"(non-expr compare)"
            moved = sp.simplify(sp.expand(a - b)) != 0 if isinstance(b, sp.Basic) else (b != a)
        any_moved = any_moved or moved
        print(f"\n--- {key} ---")
        print(f"BEFORE: {b}")
        print(f"AFTER:  {a}")
        print(f"DIFF:   {diff}")
        print(f"MOVED:  {moved}")

    # Expected jet contribution: dh * (rho4_bg0 / W0)
    rho4_bg0 = base.rho_br / base.W0
    expected_jet = before["dh"] * (rho4_bg0 / base.W0)
    print("\n--- Expected new jet piece dh * ∂_w ρ⁰ = dh * (ρ4_bg0/W0) ---")
    print(f"  expected_jet (raw) = {expected_jet}")
    print(f"  finalize(eps*expected_jet) = {base.finalize(base.epsilon * expected_jet)}")
    dens_diff = sp.simplify(sp.expand(after["FACE_SHIFT.density"] - before["FACE_SHIFT.density"]))
    print(f"  observed density DIFF = {dens_diff}")
    print(
        f"  DIFF matches finalize(eps*jet)? "
        f"{sp.simplify(dens_diff - base.finalize(base.epsilon * expected_jet)) == 0}"
    )

    print("\n" + "=" * 72)
    print(f"ABLATION_MOVED_ANY_OBJECT = {any_moved}")
    print(
        "JUDGEMENT: background normal jet pathway is LOAD-BEARING iff density MOVED\n"
        "and pressure (still bg=0) did NOT. Matching dh*(∂_w ρ⁰) confirms the engine\n"
        "implements δ[f]=δf(h_s⁰)+δh·∂_w f⁰ via actual d/dw of the supplied background,\n"
        "not a fabricated free-premise symbol — and that under the true supplied\n"
        "w-independent density the jet correctly vanishes."
    )


if __name__ == "__main__":
    main()
