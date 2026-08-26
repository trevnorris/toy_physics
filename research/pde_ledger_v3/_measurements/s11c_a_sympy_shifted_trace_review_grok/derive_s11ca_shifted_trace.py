#!/usr/bin/env python3
"""Independent derivation of S11c-a §3c shifted-trace linearisation.

Compose f = f0 + eps*df at the moving face h_s = s*W_bg/2 + eps*dh,
with W_bg = W0*(1 + eta*w1). Expand to first order in eps and first
order in eta. Print the linearised TRACE for:
  (a) background identically zero
  (b) background with genuinely nonzero normal derivative

Source of truth only: S11c_a_SHARED_PHYSICS.md §2a, §2d, §3c.
Does NOT import or read the engine.
"""
from __future__ import annotations

import sympy as sp

eps, eta = sp.symbols("epsilon eta", real=True)
s, W0, w1, dh = sp.symbols("s W0 w1 delta_h", real=True)
x, w, t = sp.symbols("x w t", real=True)

# Bookkeeping: W_bg and background face (§2a, §3c)
W_bg = W0 * (1 + eta * w1)
h_s0 = s * W_bg / 2  # background face
h_s = h_s0 + eps * dh  # moving face

# Generic bulk field as Taylor jet in w about a reference height w_ref:
#   f0(w) = f0_ref + (w - w_ref)*dw_f0_ref + ...
#   df(w) = df_ref  + (w - w_ref)*dw_df_ref  + ...
# Symbols below are the jets at the FLAT reference face w_flat = s*W0/2.
w_flat = s * W0 / 2
f0_flat, dw_f0_flat = sp.symbols("f0_at_flat  dw_f0_at_flat", real=True)
df_flat, dw_df_flat = sp.symbols("df_at_flat   dw_df_at_flat", real=True)


def field_at(w_eval, f0_flat_val, dw_f0_flat_val, df_flat_val, dw_df_flat_val):
    """First-order jet of f = f0 + eps*df about the flat face."""
    f0 = f0_flat_val + (w_eval - w_flat) * dw_f0_flat_val
    df = df_flat_val + (w_eval - w_flat) * dw_df_flat_val
    return f0 + eps * df


def expand_trace(label, f0_flat_val, dw_f0_flat_val, df_flat_val, dw_df_flat_val):
    """Compose at moving face, expand O(eps^1)*O(eta^0..1), extract delta[trace]."""
    f_face = field_at(h_s, f0_flat_val, dw_f0_flat_val, df_flat_val, dw_df_flat_val)
    # Exact composition minus background value at background face, then /eps,
    # truncate to O(eps^0) in the quotient and O(eta^1) in shape.
    f_bg_at_hs0 = field_at(h_s0, f0_flat_val, dw_f0_flat_val, 0, 0)
    # Background face value is the eps=0 part of f at h_s0
    f_bg_at_hs0 = sp.simplify(f_bg_at_hs0.subs(eps, 0))

    delta_trace_exact = sp.simplify((f_face - f_bg_at_hs0) / eps)
    # Drop O(eps) remainders inside the already-divided object
    delta_trace = sp.series(delta_trace_exact, eps, 0, 1).removeO()
    # First shape order in eta
    delta_trace = sp.series(delta_trace, eta, 0, 2).removeO()
    delta_trace = sp.expand(delta_trace)

    # Textbook §3c formula evaluated carefully:
    #   delta[f(x,h_s)] = df(x, h_s0) + dh * dw_f0(x, h_s0)
    # with jets of df and dw_f0 themselves expanded about the flat face to O(eta).
    dw_f0_at_hs0 = dw_f0_flat_val  # first-order jet: dw_f0 constant in w for this probe
    df_at_hs0 = df_flat_val + (h_s0 - w_flat) * dw_df_flat_val
    law = sp.expand(sp.series(df_at_hs0 + dh * dw_f0_at_hs0, eta, 0, 2).removeO())

    # Frozen-wrong alternative: evaluate perturbation at flat face, drop shift of face in df
    frozen_wrong = sp.expand(df_flat_val + dh * dw_f0_flat_val)

    print("=" * 72)
    print(f"CASE: {label}")
    print("-" * 72)
    print(f"W_bg               = {W_bg}")
    print(f"h_s0 = s*W_bg/2    = {sp.expand(h_s0)}")
    print(f"h_s  = h_s0+eps*dh = {sp.expand(h_s)}")
    print(f"w_flat = s*W0/2    = {w_flat}")
    print(f"f0_flat            = {f0_flat_val}")
    print(f"dw_f0_flat         = {dw_f0_flat_val}")
    print(f"df_flat            = {df_flat_val}")
    print(f"dw_df_flat         = {dw_df_flat_val}")
    print()
    print(f"COMPOSITION delta[f(x,h_s)]  (O(eps^1), O(eta^0..1)):")
    print(f"  {delta_trace}")
    print()
    print(f"SECTION 3c LAW  df(x,h_s0) + dh*dw_f0(x,h_s0)  (O(eta^0..1)):")
    print(f"  {law}")
    print()
    residual = sp.simplify(sp.expand(delta_trace - law))
    print(f"RESIDUAL (composition - law): {residual}")
    print()
    print(f"FROZEN-AT-FLAT (wrong) df(flat)+dh*dw_f0(flat):")
    print(f"  {frozen_wrong}")
    print(f"RESIDUAL (composition - frozen): {sp.simplify(sp.expand(delta_trace - frozen_wrong))}")
    print()

    # Grade census
    terms = sp.Add.make_args(sp.expand(delta_trace)) if delta_trace != 0 else ()
    print("TERM CENSUS (each monomial and whether it carries eta):")
    if not terms:
        print("  (zero)")
    for term in terms:
        has_eta = term.has(eta)
        print(f"  term={term}   has_eta={has_eta}")
    print()
    return delta_trace, law, residual


def main():
    print("INDEPENDENT DERIVATION — S11c-a §3c shifted-trace law")
    print("Truncation: first order in epsilon (wave), first order in eta (shape).")
    print()

    # (a) background identically zero (velocity / pressure / current in this scope)
    a_comp, a_law, a_res = expand_trace(
        "(a) background identically zero  [f0=0, dw_f0=0]",
        f0_flat_val=0,
        dw_f0_flat_val=0,
        df_flat_val=df_flat,
        dw_df_flat_val=dw_df_flat,
    )

    # (b) background with genuinely nonzero normal derivative
    b_comp, b_law, b_res = expand_trace(
        "(b) background nonzero WITH nonzero normal derivative",
        f0_flat_val=f0_flat,
        dw_f0_flat_val=dw_f0_flat,
        df_flat_val=df_flat,
        dw_df_flat_val=dw_df_flat,
    )

    # (c) density-like: nonzero background but w-independent (supplied ρ_4D,bg⁰)
    c_comp, c_law, c_res = expand_trace(
        "(c) density-like: nonzero bg, w-independent  [dw_f0=0]  (supplied §3c density)",
        f0_flat_val=f0_flat,
        dw_f0_flat_val=0,
        df_flat_val=df_flat,
        dw_df_flat_val=dw_df_flat,
    )

    print("=" * 72)
    print("SUMMARY — what the law requires for physical fields in this scope")
    print("=" * 72)
    print(
        """
From composition:
  delta[f(x,h_s)] = df(x, h_s0) + dh * dw_f0(x, h_s0)
  with h_s0 = s*W_bg/2 = s*W0/2 + eta*(s*W0*w1/2)

Expanding df(x,h_s0) about the flat face to O(eta):
  df(x,h_s0) = df(flat) + eta*(s*W0*w1/2)*dw_df(flat) + O(eta^2)

So two structural pieces, distinct in role:
  (1) PERTURBATION EVALUATION FACE: df is evaluated at h_s0 (carries eta through W_bg),
      NOT frozen at w = s*W0/2. The eta*dw_df piece is O(eps)*O(eta) and is retained
      at first shape order (§2a: do not drop mixed bookkeepers).
  (2) BACKGROUND NORMAL JET: dh * dw_f0(x,h_s0). This is present only when the
      supplied background has nonzero dw_f0. Per §3c / §2d in this scope:
        v_bulk^0 = 0,  delta_p^0 = 0,  j^0 = 0  =>  dw_f0 = 0 for those
        rho_4D,bg^0 depends on in-plane anchor, not on w  =>  dw_rho0 = 0
      Therefore NO free premise jet may be invented; the jet term vanishes because
      the supplied background makes it vanish, not because the law was truncated.

Conormal on an ARBITRARY probe field is not one of the supplied physical fields:
its shifted-evaluation term (the eta*dw_df piece for a generic probe) must be left
intact — it is an operator identity, not a physical-field specialisation.
"""
    )
    print(f"case (a) residual composition-vs-law: {a_res}")
    print(f"case (b) residual composition-vs-law: {b_res}")
    print(f"case (c) residual composition-vs-law: {c_res}")


if __name__ == "__main__":
    main()
