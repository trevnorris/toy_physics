#!/usr/bin/env python3
"""
Band 181-253 SCRIPT/OUTPUT numbering fixes — explicit, orchestrator-vetted edit list.

Like the band 091-180 applier, this band's actionable content-stale cross-refs are few; each
was individually verified by the orchestrator against the cited stage's filename-stem + the
deliverable's content + corroborating siblings/cross-engine twin. Source-only (.py/.wl); .txt
outputs are refreshed by re-running afterward. Every edit is digit-only (strip-all-digits guard),
applied on its exact line with a negative look-ahead so `Stage 24` never matches inside `Stage 249`.
Idempotent. Run with --apply to write.

Vetted edits (see redteam/script_output_maps/BAND_181-253_RESULT.md for the full evidence):
  stage182 x2   Stage 249 -> Stage 181  (zeta-cancellation / support-blindness "lives upstream";
                                         stage181=coherent_tracking_defect owns "zeta drops out of
                                         T^2 and R_target"; stage249=helicity export owns eta_h
                                         cancellation, not zeta; 249>182 so "upstream" was impossible)
  stage206 x4   Stage 239 -> Stage 205  (tau_log2 "quadratic log predictor" is DEFINED at stage205
                                         =directional_hessian_and_quadratic_root_refinement [19 hits];
                                         stage239 has 0; incl. the compound "Stage 206/239"->"206/205"
                                         miss on the expect_zero label)
  stage218 x1   Stage 249 -> Stage 215  (the "600 + 5*2*54" = 1140 full support<=4 budget is owned by
                                         stage215=full_primitive_quadruple_ranking_theorem [line 208];
                                         249>218 "upstream" impossible; 232=249-17 has no 1140)
  stage227 x1   Stage-225 -> Stage-223  (the K compatibility value 24.4737548792909... is the stage223
                                         =..._primitive_branch_compatibility... deliverable; the
                                         cross-engine .py twin labels the same computation "Stage 223
                                         branch"; -2 content drift, decided by the twin)
  stage233 x5   Stage 239 -> Stage 188  (the branch-observable compiler B_*/dln_Rtr/Xi/R packet is
                                         OWNED by stage188 "THE EXACT FIRST-ORDER OBSERVABLE COMPILER"
                                         with identical symbols+formulas; stage239=rigid_mouth normal
                                         form owns none of it; stage222=239-17 is pole-census, no match)
                Stage 240 -> Stage 223  (Pbar/compat point P0_target_compat originates at stage223)
                Stage 241 -> Stage 224  (ceilings dict + budgets owned by stage224; stage226 calls them
                                         "Transported Stage 224 budgets"; 240/241-17=223/224 corroborates)
"""
import os, re, sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DRY = "--apply" not in sys.argv

S182_PY = "scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py"
S182_WL = "mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl"
S206_PY = "scripts/moving_throat_pde_stage206_certified_ray_ranking_and_local_bracketing_theorem_sympy_audit.py"
S218_PY = "scripts/moving_throat_pde_stage218_full_support_cardinality_5_completion_and_local_mixed_ray_search_closure_sympy_audit.py"
S227_WL = "mathematica/moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_mathematica_audit.wl"
S233_PY = "scripts/moving_throat_pde_stage233_rigid_mouth_orbit_lock_compiler_and_the_static_turbulence_gate_sympy_audit.py"

# (relpath, line, old_token, new_token)
EDITS = [
    (S182_PY, 140, "Stage 249", "Stage 181"),
    (S182_WL, 135, "Stage 249", "Stage 181"),
    (S206_PY, 131, "Stage 239", "Stage 205"),
    (S206_PY, 133, "Stage 239", "Stage 205"),
    (S206_PY, 142, "Stage-239", "Stage-205"),
    (S206_PY, 145, "Stage 206/239", "Stage 206/205"),
    (S218_PY, 347, "Stage 249", "Stage 215"),
    (S227_WL, 154, "Stage-225", "Stage-223"),
    (S233_PY, 20, "Stage 239", "Stage 188"),
    (S233_PY, 32, "Stage 239", "Stage 188"),
    (S233_PY, 113, "Stage 241", "Stage 224"),
    (S233_PY, 129, "Stage 240", "Stage 223"),
    (S233_PY, 129, "Stage 241", "Stage 224"),
]

def strip_digits(s):
    return re.sub(r"\d", "", s)

byfile = {}
for relpath, line, old, new in EDITS:
    byfile.setdefault(relpath, []).append((line, old, new))

applied = skipped = files_changed = 0
skips = []
for relpath in sorted(byfile):
    path = os.path.join(ROOT, relpath)
    if not os.path.exists(path):
        skips.append(f"{relpath}: FILE NOT FOUND"); continue
    lines = open(path, encoding="utf-8").read().split("\n")
    changed = False
    for line, old, new in sorted(byfile[relpath]):
        if line < 1 or line > len(lines):
            skips.append(f"{relpath}:{line} `{old}` out of range"); skipped += 1; continue
        if strip_digits(old) != strip_digits(new):
            skips.append(f"{relpath}:{line} NON-DIGIT-ONLY `{old}`->`{new}`"); skipped += 1; continue
        cur = lines[line-1]
        pat = re.escape(old) + r"(?!\d)"
        cur2, n = re.subn(pat, new, cur, count=1)
        if n == 0:
            if old not in cur and new in cur:
                skipped += 1  # idempotent (already applied)
            else:
                skips.append(f"{relpath}:{line} `{old}` NOT FOUND (got: {cur.strip()[:70]!r})"); skipped += 1
            continue
        lines[line-1] = cur2
        applied += 1
        changed = True
    if changed:
        files_changed += 1
        if not DRY:
            open(path, "w", encoding="utf-8").write("\n".join(lines))

print("DRY-RUN (pass --apply to write)" if DRY else "APPLIED")
print(f"  edits applied: {applied}   skipped: {skipped}   files changed: {files_changed}   total: {len(EDITS)}")
for s in skips:
    print("    -", s)
