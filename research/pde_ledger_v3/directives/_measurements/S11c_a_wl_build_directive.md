# Measurements twin — S11c_a_wl_build_directive.md

Every claim carries the command that produced it (CLAUDE.md rule 2).

## Physics authority + tag grammar + blindness are all SPEC-native
CMD: sed -n '528,551p' research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
  ⇒ §7: grammar `<ENGINE>_S11CA_<QUANTITY>` (PY|WL); "The SymPy engine imports and carries the S11b LEDGER;
     the Wolfram engine RE-DERIVES the supplied §§1–3 inputs WITHOUT an import." ⇒ blindness + WL_ tags are spec-mandated.
CMD: grep -n 'tag = f"PY_S11CA' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py ⇒ line 359 (sibling emits PY_S11CA_<QUANTITY>; comparator joins on S11CA_<QUANTITY> stem).
CMD: sed -n '452,481p' research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md ⇒ §5a route-2 = material-coordinate w'=[w-zeta_c]/[W_bg+delta_W], mapped back to Eulerian; one-sided mutates ONLY the direct route.

## Blindness control anchor
CMD: sed -n '12,22p' research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md
  ⇒ "Each step's Wolfram engine imports nothing and re-derives independently — that is the only blindness control in this design."

## Structural template (mechanical precedent, not authority)
CMD: cat research/pde_ledger_v3/directives/S11b_wl_build_directive.md ⇒ thin WL build directive shape (authority/boundary, blind engine, mechanical precedent, run discipline 2 seats + kill criteria, executable acceptance blindness/flush/no-stray-write, conflicts).
CMD: ls research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl ⇒ named as this engine's mechanical precedent (another step's physics; do NOT import/transcribe).

## T7 comparator contract (frozen; downstream, NOT built here)
CMD: sed -n '53,60p' research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md
  ⇒ T7: comparator frozen before it sees either output; REJECTS a native boolean as a residual operand ⇒ WL engine must emit booleans as the CAS object, never host-native.

## §5a — the guard SECTION WAS REMOVED in the fold (2 directive legs, 2026-08-23)
  Codex + Grok both caught a BLINDNESS LEAK: my §5a guard grouped the objects as "T-c, T-d, T-i, T-g (and T-h)"
  = the SIBLING's defect partition (spec §5a order is T-g,T-h,T-c,T-d,T-i, all five equal), and collapsed the two
  w' anchoring maps (W_bg(x(X,t)) vs W_bg(X)) into one displayed formula. FIX: removed the whole §5a guard section.
  The directive now relies on the AUTHORITY BLOCK's §5 pointing (spec-faithful, no re-render) + the WL engine's OWN
  2 legs to enforce §5a route independence (rule 12: catch in review, ⛔ never warn the blind builder). BLIND-SAFE.
  Also folded: mu_R≡mu_R_bg conflation (F3), §4/§5/§6/§7 tag-list narrowing→§§2-6 (F4/Codex), boolean→relational
  retaining operands (Codex F3), NOT_ESTABLISHED(§0)→§8 builder report (F5, spec §0 has 0 occurrences).
