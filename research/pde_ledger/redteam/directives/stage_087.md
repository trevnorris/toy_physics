---
unit_id: 087
batch: III.5
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
verification_status: scripts_pass_pending_verifier
needs_user_resolution: true
---

# Codex directive — unit 087

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim (with banner-level target_mismatch)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_087.tex:14-20` quote:
  > `\boxed{\rho_\alpha=\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}\quad\text{is the only remaining support-side input from the outgoing branch.}}`
  > `Given \(\rho_\alpha\), the Family--1 support/source verdict follows from Stage~086.  No separate dependence on \(N_Q^{\rm target}\), \(\beta_0\), \(s_-\), \(\lambda_-\), or \(\widehat m_-\) remains.`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md:22-38` quote (paraphrased): the variables `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req` all drop out, leaving only `rho_alpha`.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:19` quote: `banner("STAGE 70 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:32` quote: `banner["STAGE 070 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE"];`
- Neither script introduces any of `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req` or checks their cancellation.

## Resolve before fix_loop

The paper card says stage 087's `Output` is the *cancellation* claim — that the explicit Family-1 support theorem has lost dependence on `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req` and now depends only on `rho_alpha`. The scripts do not check that. They check properties of an auxiliary `zeta_req(rho_alpha, eps_blk)` at `eps_blk = 0`, with banners reading "STAGE 70". Is stage 087 (a) a status/checkpoint consolidation whose cancellations were verified in upstream stages (068-069, 085, 086), or (b) a load-bearing claim that should have explicit cancellation checks in its own scripts?

Possible directions (the user picks one):

- (a) **Stage 087 is a status/checkpoint consolidation.** No new math. Direction: update both script banners from "STAGE 70" / "STAGE 070" to "STAGE 087", update the script docstrings/comments to explicitly say "carry-forward from stages 085-086; cancellations verified upstream", and keep the existing `zeta_req` window-value sanity checks as downstream consistency. F2 must still be resolved (the existing checks are tautological even as sanity probes — see F2 question).
- (b) **Stage 087 is load-bearing.** Direction: rewrite both scripts to introduce `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req` as free symbols, exhibit the explicit support-theorem expression in those symbols, and assert that under the closure `mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5)` the dependence collapses to `rho_alpha = alpha_req/alpha_mix` only. This is a substantial new derivation.
- (c) **Mixed.** Keep the scripts as carry-forward (a) but add a minimal symbolic check that `zeta_req` is the right shape for the upstream window (Stage 069/086) — clarify exactly which upstream identity is being repeated here.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — tautological_check

**Subtype:** tautological_check (linked to F1's resolution)

**Paper side / upstream source:** the threshold values `rho_suff = 3.46622291347846`, `rho_fail = 3.46752913273870`, `rho_max = 3.46752922945601` appear in `notes/stages/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md:60-63` but their derivation is in upstream stages (068-069).

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:27` — `expect_zero("unblocked zeta_req", zeta_req.subs(eps_blk, 0) - (rho_alpha - 1))` — algebraically `(rho_alpha - 1) - (rho_alpha - 1) == 0`.
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:41-43` — `expect_zero("zeta_suff - 2.46622291347846", sp.N(zeta_suff - sp.Float("2.46622291347846"), 18))` and two siblings, each comparing `rho_X - 1` to a hand-typed literal that equals `rho_X - 1`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:44` and `:58-60` — identical pattern.

## Resolve before fix_loop

If F1 direction = (a) status/checkpoint, are these `zeta_suff/zeta_fail/zeta_max` numeric checks meant to test (i) the literal threshold values (in which case they should compare against an upstream-derived source like a result quoted from Stage 069's script output, not against `rho_X - 1`); or (ii) just to print the window values for visual confirmation (in which case they should be `print` statements, not `expect_zero`/`expectApprox`)? Codex needs to know whether to delete or replace them.

If F1 direction = (b), F2 is resolved by the new derivation work and the existing tautologies are simply replaced.

The orchestrator will not invoke Codex on this finding until F1 is resolved.

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:32-63`

**Issue:** The Mathematica script mirrors the SymPy script's variable choreography (same `zeta_req` definition, same derivative print, same eps_blk -> 0 substitution, same three numeric ratio substitutions, same tautological residuals). Only the `d zeta_req exact formula` check is independent. This violates the two-engine independence policy.

**Required change:**
HOLD until F1 and F2 resolve. Once the user has chosen a direction for F1/F2, re-author the Mathematica script (if it remains in scope) so that any kept assertion is derived independently from the physical premises rather than mirroring the sympy step list. If F1 direction = (a) (status/checkpoint), F3 can be closed as won't-fix.

**Verification:**
After F1/F2 land, the verifier confirms the mathematica script no longer matches the sympy script line-for-line in its body.
