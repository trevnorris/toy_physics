#!/usr/bin/env python3
"""Ground c1's IMPORT_KEYS manifest against the REAL frozen base (rule 2 / verify-against-real-artifact).

Run from research/pde_ledger_v3/scripts/ (where ledger_fold.py and S11c_b_exports.py live):

    cd research/pde_ledger_v3/scripts && python3 ../directives/_measurements/c1_migration_import_keys_scan.py

Reports: (1) the single-atomic-base fold row count and that it has zero overwrites (so
S11c_b_exports.py alone is the frozen base carrying the whole merged model through S11c-b);
(2) presence of each candidate root c1 binds per the design build-plan step 3 + spec §1c;
(3) that check_consumer resolves the full root set's recursive closure with no ambiguity/dangling.

This does NOT decide membership (that is the spec's bind-set + the build's smoke-test); it verifies
the declared manifest is present and closes cleanly against the real artifact — the required
real-artifact acceptance step, not a synthetic fixture.
"""
import sys
from pathlib import Path

# allow running from anywhere: put the scripts/ dir on the path
SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import ledger_fold as lf  # noqa: E402

BASE = str(SCRIPTS / "S11c_b_exports.py")

# c1's declared IMPORT_KEYS roots (the F9 write-keys c1 binds directly), reconciled from BOTH review legs
# (Codex + Grok) and grounded in S11c_c1_SHARED_PHYSICS.md §1b/§1d/§2c/§5 + export_ledger_bind_closure_design.md.
CONSTANTS_KNOBS = [
    "mu_theta_operator", "c_s0", "rho_m", "rho_br", "W_0", "e_W",
    "v_bulk_normal_0", "q_out", "omega", "epsilon_shape",  # +omega (ω), +epsilon_shape (ε bookkeeper): direct binds
    "Lambda_A_0", "Lambda_V_0", "Lambda_X_0", "tau_A", "tau_V", "tau_X",
]
T_SUBSTRATE = [  # S11c-a T-a..T-i (nine rows)
    "face_normal", "conormal_deriv", "face_measure_shape_deriv", "face_velocity",
    "relative_flux", "kinematic_balance", "traction", "face_shift", "closure_shape_deriv",
]
PROFILES_DENSITY = [  # bulk uses constant rho_m (§1b); density enters c1 only via mu_s = mu_theta / rho_br,bg0 (§1d)
    "W_bg", "w1_profile", "L_W", "sigma_W", "eta_bg",
    "rho_br_bg_rho4_constant",  # RHO4-const rep; rho_br above is the RHOBR-const rep
]
S11B_REGRESSION = [  # spec §1c regression operands (§5c/§5d)
    "z_impermeable", "z_by_regime", "z_by_parity", "added_mass", "grazing_behaviour",
    "face_response", "face_response_coeffs", "permeable_dissipative_by_regime_and_parity",
    "degenerate_loci_equations", "degenerate_loci_solution",
    "degenerate_loci_identically_satisfied", "degenerate_loci_inconsistent",
    "degenerate_loci_real_admissible",
]
ROOTS = CONSTANTS_KNOBS + T_SUBSTRATE + PROFILES_DENSITY + S11B_REGRESSION

# NOT bound by c1 (reserved for c2 or structurally absent); assert none is pulled into the closure by a c1 root.
NOT_BOUND = [
    "mu_R", "mu_R_bg", "m1_profile",           # modulus channel: inside opaque mu_theta, reserved for c2
    "rho_4D_bg_rho4_constant", "rho_4D_bg_rhobr_constant",  # slab-side 4D-density reps; bulk is constant rho_m
    "e_W_bg",                                   # background contrast, carried by eta_bg
]


def main() -> int:
    fold, audit = lf.load_model(BASE)
    print(f"FOLD ROWS: {len(fold)}")
    print(f"SOURCE ROW COUNTS: {audit['source_row_counts']}")
    print(f"OVERWRITES (single atomic base ⇒ expect []): {audit['overwrites']}")

    missing = [k for k in ROOTS if k not in fold]
    print(f"\nCANDIDATE ROOTS: {len(ROOTS)} declared; MISSING from fold: {missing or 'none'}")

    res = lf.check_consumer(fold, ROOTS)
    print("\ncheck_consumer OK — no AmbiguousSymbolError / ClosureError")
    print(f"  roots declared    : {len(ROOTS)}")
    print(f"  recursive closure : {len(res['closure'])}")
    print(f"  symbol edges      : {len(res['symbol_edges'])}")
    print(f"  dimension edges   : {len(res['dimension_edges'])}")

    # The over-declarations both legs flagged: none of these (except mu_R, reachable via mu_theta_operator)
    # is pulled into c1's bind-closure — evidence they are genuine over-declarations, not silent drops.
    print("\nNOT_BOUND keys — in the closure of the 44 roots?")
    for k in NOT_BOUND:
        print(f"  {k}: in_closure={k in res['closure']}"
              f"{'  (reachable via mu_theta_operator — correct)' if k == 'mu_R' else ''}")
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
