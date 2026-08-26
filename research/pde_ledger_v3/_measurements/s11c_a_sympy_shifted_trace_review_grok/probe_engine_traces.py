#!/usr/bin/env python3
"""Probe S11c-a SymPy engine against independent §3c derivation.

Imports the /tmp engine COPY. Asks of each traced physical object:
  which line computed this? Does it contain a fabricated BACKGROUND normal-jet?
  Is the perturbation evaluated at h_s0 = s W_bg/2?
  Does CONORMAL keep the shifted-evaluation term on a generic probe?
Reports any assert preceding a guarded value.
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path("/tmp/s11c_a_review_grok")
sys.path.insert(0, str(ROOT))
import engine_copy as eng  # noqa: E402


def free_symbols_named(expr, substrings):
    names = []
    if not isinstance(expr, sp.Basic):
        return names
    for sym in expr.free_symbols:
        for sub in substrings:
            if sub in str(sym):
                names.append(str(sym))
                break
    return sorted(set(names))


def main():
    print("=" * 72)
    print("PROBE — engine traced-bulk objects vs §3c")
    print("=" * 72)

    # --- asserts ---
    src = (ROOT / "engine_copy.py").read_text()
    assert_lines = [
        (i + 1, line.rstrip())
        for i, line in enumerate(src.splitlines())
        if "assert " in line and not line.strip().startswith("#")
    ]
    print("\nASSERTS in engine (any assert preceding a guarded value):")
    if not assert_lines:
        print("  (none)")
    else:
        for n, line in assert_lines:
            print(f"  L{n}: {line}")

    # --- build one face source and inspect raw traces before finalize ---
    source = eng.build_face_source(
        "LAB_HELD", 1, "DELTA_W", "RHO4_CONSTANT", route="EULERIAN",
    )
    print("\n--- FaceSource structural fields (LAB_HELD, UPPER, DELTA_W, RHO4) ---")
    print(f"h0 (exact, pre-finalize) = {source.h0}")
    print(f"dh = {source.dh}")
    print(f"h0 == face*W_bg/2 ? {sp.simplify(source.h0 - eng.W_bg / 2) == 0}")

    # Linearised physical traces (the §3c objects)
    p_lin = sp.diff(source.pressure_trace_exact, source.parameter).subs(source.parameter, 0)
    v_lin = [
        sp.diff(component, source.parameter).subs(source.parameter, 0)
        for component in source.bulk_velocity_trace_exact
    ]
    rho_lin = sp.diff(source.density_trace_exact, source.parameter).subs(source.parameter, 0)
    j_lin = [
        sp.diff(component, source.parameter).subs(source.parameter, 0)
        for component in source.current_trace_exact
    ]

    print("\n--- Linearised traces (pre-finalize, shape derivative) ---")
    print(f"pressure_lin = {p_lin}")
    print(f"velocity_lin[0..3] = {v_lin}")
    print(f"density_lin = {rho_lin}")
    print(f"current_lin[0..3] = {j_lin}")

    # Fabricated BACKGROUND normal-jet symbols? Look for anything like d_w_* of a
    # background (not the wave-field d_w_delta_* coordinates).
    bg_jet_needles = (
        "d_w_v_bulk_0", "d_w_f0", "dw_f0", "d_w_pressure_0", "d_w_rho_0",
        "background_normal", "v_bulk_normal_jet", "dw_bg",
    )
    all_lin = [p_lin, rho_lin, *v_lin, *j_lin]
    print("\n--- Fabricated BACKGROUND normal-jet symbols in linearised traces? ---")
    found_fabricated = False
    for label, expr in (
        ("pressure", p_lin),
        ("velocity", v_lin),
        ("density", rho_lin),
        ("current", j_lin),
    ):
        flat = expr if isinstance(expr, sp.Basic) else sp.Tuple(*expr)
        hits = free_symbols_named(flat, bg_jet_needles)
        # Also check whether v_dr (drain) appears — §3c says inert, not a jet premise
        has_v_dr = flat.has(eng.v_dr) if isinstance(flat, sp.Basic) else any(e.has(eng.v_dr) for e in expr)
        print(f"  {label}: fabricated_bg_jet_hits={hits or 'none'}; has_v_dr={has_v_dr}")
        if hits or has_v_dr:
            found_fabricated = True
    print(f"  FABRICATED_BACKGROUND_JET = {found_fabricated}")

    # Wave-field normal jets (coordinates) SHOULD appear via the affine ansatz
    # when evaluated off the flat face: (h0 - s W0/2)*dw_delta_*
    print("\n--- Wave-field normal jets (COORDINATE, not premises) present? ---")
    for label, expr in (
        ("pressure", p_lin),
        ("velocity_w_comp", v_lin[3]),
        ("density", rho_lin),
        ("current_w_comp", j_lin[3]),
    ):
        hits = free_symbols_named(expr, ("d_w_delta_",))
        print(f"  {label}: {hits or 'NONE — check whether h0-shift of perturbation is absent'}")

    # Expand with profile subs: does eta appear through (h0-flat)*dw?
    print("\n--- After PROFILE_SUBS + O(eta): does shifted eval of perturbation carry eta? ---")
    for label, expr in (
        ("pressure", p_lin),
        ("velocity_w", v_lin[3]),
        ("density", rho_lin),
    ):
        subbed = expr.xreplace(eng.PROFILE_SUBS)
        # first order in eta_bg
        eta0 = subbed.subs(eng.eta_bg, 0)
        eta1 = sp.diff(subbed, eng.eta_bg).subs(eng.eta_bg, 0)
        print(f"  {label} O(eta^0) = {eta0}")
        print(f"  {label} O(eta^1) = {eta1}")
        print(f"  {label} O(eta^1) is zero? {sp.simplify(eta1) == 0}")

    # Independent derivation expectation for pressure (bg=0):
    #   delta_p_flat + (h0 - W0/2)*d_w_delta_p
    # with h0 = W_bg/2
    expected_p = eng.delta_p_plus + (source.h0 - eng.W0 / 2) * eng.dw_delta_p[1]
    print("\n--- Pressure linearisation vs independent formula ---")
    print(f"  engine  = {sp.expand(p_lin)}")
    print(f"  expect  = {sp.expand(expected_p)}")
    print(f"  residual= {sp.simplify(sp.expand(p_lin - expected_p))}")

    expected_vw = eng.delta_v_bulk[1][3] + (source.h0 - eng.W0 / 2) * eng.dw_delta_v_bulk[1][3]
    print("\n--- Velocity_w linearisation vs independent formula ---")
    print(f"  engine  = {sp.expand(v_lin[3])}")
    print(f"  expect  = {sp.expand(expected_vw)}")
    print(f"  residual= {sp.simplify(sp.expand(v_lin[3] - expected_vw))}")

    # Density: bg = rho4_ref (w-independent) for RHO4_CONSTANT; jet term from bg vanishes;
    # perturbation still gets (h0-flat)*dw_delta_rho
    rho4_ref = eng.rho_br / eng.W0
    expected_rho = eng.delta_rho4_face[1] + (source.h0 - eng.W0 / 2) * eng.dw_delta_rho4_face[1]
    # Plus possible dh * dw_rho0 = 0
    print("\n--- Density linearisation vs independent formula (dw_rho0=0) ---")
    print(f"  density_background used = {source.rho4_bg_exact.subs(source.parameter, 0)}")
    print(f"  engine  = {sp.expand(rho_lin)}")
    print(f"  expect  = {sp.expand(expected_rho)}")
    print(f"  residual= {sp.simplify(sp.expand(rho_lin - expected_rho))}")

    # Frozen-at-flat wrong formula residual
    frozen_p = eng.delta_p_plus  # no (h0-flat)*dw
    print("\n--- If engine froze perturbation at flat face, pressure would equal delta_p_plus ---")
    print(f"  engine - frozen = {sp.simplify(sp.expand(p_lin - frozen_p))}")
    print(f"  (nonzero means engine DID evaluate at h_s0, not freeze at flat)")

    # --- CONORMAL: shifted-evaluation term on generic probe ---
    print("\n--- CONORMAL_DERIV generic probe shifted-evaluation term ---")
    con = eng.conormal_raw(source)
    # con = Tuple(background, eps*derivative, full)
    deriv = con[1] / eng.epsilon
    print(f"  conormal shape-deriv (raw) = {deriv}")
    dw_syms = free_symbols_named(deriv, ("d_w_trace_grad_f",))
    print(f"  contains d_w_trace_grad_f_* (shifted-eval / face-motion jet)? {dw_syms or 'NO'}")
    # Also check that dh multiplies those jets
    has_dh_times_jet = False
    expanded = sp.expand(deriv)
    for term in sp.Add.make_args(expanded):
        if term.has(source.dh) and any(term.has(s) for s in eng.conormal_raw.__code__.co_consts if False):
            pass
        # simpler: any term containing both a d_w_trace_grad_f and dependence on dh
        if any(str(s).startswith("d_w_trace_grad_f") for s in term.free_symbols):
            if source.dh in term.free_symbols or term.has(source.dh):
                has_dh_times_jet = True
                print(f"  term with dh*jet: {term}")
    # dh may have been substituted already as zeta expression
    # source.dh for DELTA_W LAB = face * W0 * e_W / 2 = W0*e_W/2
    print(f"  source.dh = {source.dh}")
    # Re-check: expand conormal construction directly
    grad_f = tuple(sp.Symbol(f"trace_grad_f_{i}") for i in range(1, 5))
    dw_grad_f = tuple(sp.Symbol(f"d_w_trace_grad_f_{i}") for i in range(1, 5))
    # From source lines 947-953:
    traced = tuple(grad_f[i] + source.parameter * source.dh * dw_grad_f[i] for i in range(4))
    exact = eng.dot(source.normal_exact, traced)
    cderiv = sp.diff(exact, source.parameter).subs(source.parameter, 0)
    print(f"  recomputed conormal deriv = {cderiv}")
    print(f"  dh*dw_grad_f present in recomputed? "
          f"{any(cderiv.has(dw_grad_f[i]) for i in range(4))}")

    # --- Which lines compute which objects ---
    print("\n--- Which lines compute traced-bulk objects ---")
    print("  compose_physical_trace:     L588-601  (δ[f]= compose f0+q*df at h0+q*dh)")
    print("  physical_trace_fields:      L629-685  (bg from B0; affine pert about flat)")
    print("  affine_bulk_perturbation:   L578-585  (pert jet about s*W0/2)")
    print("  build_face_source h0:       L844      (h0 = face*W_bg/2)  [EULERIAN]")
    print("  build_material_face h0:     L776-777  (h0 = face_height|q=0)")
    print("  relative_flux_raw:          L894-897  (uses bulk_velocity_trace_exact)")
    print("  traction_raw:               L919-924  (uses pressure_trace_exact)")
    print("  kinematic_raw:              L931-943  (uses bulk_velocity_trace_exact)")
    print("  conormal_raw:               L946-953  (generic probe keeps dh*dw_grad_f)")
    print("  build_face_shift_raw:       L1062-1088 (emits all four physical traces)")
    print("  velocity_background:        L659      (=0; no free premise jet)")
    print("  pressure_background:        L660      (=0)")
    print("  density_background:         L661      (rho4_bg_exact; w-indep via density_pair)")
    print("  current_background:         L662-664  (rho0 * v0 = 0)")

    # Emit FACE_SHIFT one case via finalize to show eta content
    print("\n--- FACE_SHIFT single case after finalize (pressure, velocity_w, density) ---")
    # Clear caches that might interfere? build_face_shift uses build_face_source
    fs = eng.build_face_source("LAB_HELD", 1, "ZETA_C", "RHO4_CONSTANT")
    p = eng.finalize(eng.pressure_trace_raw(fs))
    v = eng.finalize(eng.epsilon * sp.ImmutableMatrix(
        eng.shape(fs.bulk_velocity_trace_exact, fs.parameter)
    ))
    d = eng.finalize(eng.epsilon * eng.shape(fs.density_trace_exact, fs.parameter))
    print(f"  FACE pressure finalized = {p}")
    print(f"  FACE velocity finalized = {v}")
    print(f"  FACE density finalized  = {d}")
    print(f"  pressure has eta_bg? {p.has(eng.eta_bg)}")
    print(f"  velocity[3] has eta_bg? {v[3].has(eng.eta_bg)}")
    print(f"  density has eta_bg? {d.has(eng.eta_bg)}")


if __name__ == "__main__":
    main()
