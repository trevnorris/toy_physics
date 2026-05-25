---
unit_id: 010
batch: I.1
created_at: 2026-05-25T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 010

This directive contains only `paper_misalignment` findings. The orchestrator is holding for user resolution. Do NOT edit `paper/`, `notes/`, or scripts to "fix" any paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

For each `## Resolve before fix_loop` block below, the orchestrator surfaces the question to the user; once the user picks a direction, a new directive is written specifying the concrete edit (and may then be routed to Codex if a script-side change is the chosen direction).

---

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim (C2, C3 not exercised) and target_mismatch (script's `δP_2`, `δP_4` are adjacent-but-distinct objects from the paper's `δu_2`, `δu_4`)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex:25-33` quotes:
  ```
  \begin{equation}\label{eq:stage010-du2}
  \delta u_2=\frac{D_0z_2-D_2z_0}{D_0^2},
  \end{equation}
  \begin{equation}\label{eq:stage010-du4}
  \delta u_4=
  \frac{D_0^2z_4-D_0(2D_2z_2+D_4z_0)+2D_2^2z_0}{D_0^3},
  \end{equation}
  ```
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex:40-42` `\stagefield{Output}` quote: "Stage~010 exports the transport map \eqref{eq:stage006-projected-shifts}--\eqref{eq:stage010-dp0}."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:50-52` quote:
  ```python
  P0p = N0p / D0p
  P2p = (D0p * N2p - 2 * D2p * N0p) / D0p**2
  P4p = (D0p**2 * N4p - 2 * D0p * (D2p * N2p + D4p * N0p) + 3 * D2p**2 * N0p) / D0p**3
  ```
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:58-83` asserts closed forms for `dP_0`, `dP_2`, `dP_4`. Only `dP_0` matches paper deliverable C4. No assertion on `δu_2` or `δu_4` anywhere in either engine.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl:31-63` mirrors the SymPy structure with the same gap.

Confirming algebra: the paper's `δu_2 = (D_0 z_2 - D_2 z_0)/D_0^2` is the first variation of `u_2 = -D_2/D_0` under `D_n -> D_n - \varepsilon z_n` (pure D-side). The script's `δP_2` is the first variation of `(D_0 N_2 - 2 D_2 N_0)/D_0^2` under the same shifts plus `N_n -> N_n + \varepsilon n_n` — a full ratio, and contains `n_0, n_2, N_0, N_2` contributions that `δu_2` does not. Setting `n_0 = n_2 = 0` in `δP_2` still leaves `N_2 z_0/D_0^2 + 2 N_0 z_2/D_0^2 - 4 D_2 N_0 z_0/D_0^3`, so `δP_2 \ne δu_2` even in the n-zero limit. Same situation for `δu_4` vs `δP_4`.

## Resolve before fix_loop

The paper card asserts three first-variation identities as the stage's `\stagefield{Output}`: `δu_2 = (D_0 z_2 - D_2 z_0)/D_0^2`, `δu_4 = (D_0^2 z_4 - D_0 (2 D_2 z_2 + D_4 z_0) + 2 D_2^2 z_0)/D_0^3`, and `δP_0 = (D_0 n_0 + N_0 z_0)/D_0^2`. The scripts verify `δP_0` exactly but do not verify `δu_2` or `δu_4`; instead they verify `δP_2` and `δP_4` (full numerator/denominator ratio variations including `n_n` shifts), which are distinct objects.

Which is the correct stage 010 deliverable?

Possible directions (the user picks one):

- (a) **Paper is correct** — `δu_2` and `δu_4` are the intended export. Stage 010 must add explicit `assert_zero` checks for those identities to both engines. Suggested form (sympy):
  ```python
  u2_eps = -(D2 - eps * z2) / (D0 - eps * z0)
  u4_eps = ((D2 - eps * z2)**2 - (D0 - eps * z0) * (D4 - eps * z4)) / (D0 - eps * z0)**2
  du2 = sp.diff(u2_eps, eps).subs(eps, 0)
  du4 = sp.diff(u4_eps, eps).subs(eps, 0)
  assert_zero("delta u_2", du2 - (D0*z2 - D2*z0) / D0**2)
  assert_zero("delta u_4", du4 - (D0**2 * z4 - D0*(2*D2*z2 + D4*z0) + 2*D2**2*z0) / D0**3)
  ```
  with the analogous Mathematica check using `Series[..., {eps, 0, 1}]` and `Coefficient[..., eps, 1]`. The existing `δP_2`, `δP_4` blocks may stay (as cross-checks for the bundle structure) or be removed at the user's discretion — but they do not by themselves satisfy `\stagefield{Output}`.

- (b) **Script is correct** — `δP_2` and `δP_4` (the full ratio variations) are the actual stage 010 exports. The paper card is wrong: equations `eq:stage010-du2` and `eq:stage010-du4` should be rewritten to display the explicit `δP_2` and `δP_4` closed forms (the eight-term and twelve-term expressions in sympy `:60-67` and `:71-82`), the symbol `u_n` should be renamed `P_n`, and any downstream paper text that references `\eqref{eq:stage010-du2}` or `\eqref{eq:stage010-du4}` needs the same rename.

- (c) **Both u_n (pure D-side) and P_n (full ratio) are physically distinct objects the stage should export**, but the paper card only shows the u_n side and the script only shows the P_n side. In that case both sides need extension: paper gets a `δP_n` paragraph, scripts get `δu_n` assertions. This is the only direction where both the paper card and the scripts grow.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. Codex must NOT silently rewrite the paper or the scripts to make them agree.

---

## F2 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex` (entire card): no mention of S, T, P_target, Z_0_slot, D_0_target, real Y_{20} overlap, Gaunt coefficients, real-harmonic lane multipliers, grouped trace/anomaly decomposition, the b = 3a line, primitive static Xi, Delta, Q, q_1, p_1, d_1, or any K-surface compatibility relation.
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:42` row for stage 010 reads only: "Projected Maxwell push into bundle slots / \StatusExactClosure{} / Projected Z_n, N_n slot transport for the grouped response bundle." — no further structure.
- Notes file `notes/em_projected/step_08_projected_maxwell_push_bundle_master_notes.md` (the source the script docstring references at `:2`) is not present in the committed repo — this is one of the EM-projected stages whose per-stage notes were not committed.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:85-169` verifies seven additional clusters: (1) one-pole K-surface from the pole relation `(K - B_0 - Z_{0,\mathrm{slot}} - \varepsilon z_0)(T + \varepsilon z_4) = 3(S + \varepsilon z_2)^2`; (2) normalization K-surface from `(N_0 + \varepsilon n_0)/(K - B_0 - Z_{0,\mathrm{slot}} - \varepsilon z_0) = P_{\mathrm{target}}`; (3) compatibility surface and its first variation `dcompat = n_0/P_{\mathrm{target}} - 6 S z_2/T + 3 S^2 z_4/T^2`; (4) transported-target compatibility analog with `P_{\mathrm{target,transport}} = (N_0 + \varepsilon n_0)/D_{0,\mathrm{target}}`; (5) two negative-control `assert_nonzero` mutation guards; (6) real-Y_{20}-square Gaunt overlap ratios `(1, 1/2, -1)` and a trace/anomaly decomposition along `b = 3a`; (7) primitive static Xi formula under `N_{0,\mathrm{sym}} = P^2/\Delta^2`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl:65-220` mirrors the same seven clusters (M4-M17).

## Resolve before fix_loop

The scripts certify seven content clusters (K-surface compatibility, transported-target compatibility, Gaunt lane multipliers, weak-axisymmetric trace anomaly, primitive static Xi, plus the two sign-flip mutation guards) that the paper card never mentions and that have no anchoring notes file in the repo.

Which is the correct disposition?

Possible directions (the user picks one):

- (a) **Paper card is incomplete** — these seven clusters are load-bearing for stage 010's role in the larger argument (likely needed by stages 011 P_2 bridge, 012 primitive bridge, 013-014 mouth-Taylor). Stage_010.tex needs new paragraphs describing each cluster with explicit equation labels, so that the printed paper exposes what the script verifies. Optionally the part appendix row gets a richer Verification-output column.

- (b) **Scripts over-verify, trim them** — these clusters are scaffolding from an earlier draft (notice the SymPy docstring at `:2` references `step_08_projected_maxwell_push_bundle_master_notes.md`, an EM-projected step number that was renumbered to stage 010 and has no surviving notes file). The script blocks at sympy `:85-169` and mathematica `:65-220` are removed, leaving only the `δP_0` check (and any `δu_n` content from F1) as the audit body.

- (c) **File the missing notes** — author publishes `notes/em_projected/step_08_projected_maxwell_push_bundle_master_notes.md` to anchor the seven clusters as legitimate stage 010 derivations, and the stage card may optionally cross-reference them via `\StageFile`. This is the lightest-touch option (no paper-body change, no script change) but only valid if the notes content is in fact load-bearing for stage 010 specifically (not stage 008 only).

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. Codex must NOT silently delete script content or silently extend the paper card.

---

Note: F1 and F2 are independent — the user may pick different directions for each. F1's resolution affects the `δu_n` vs `δP_n` symbol convention; F2's resolution affects the seven extra clusters. After the user has chosen, a follow-up directive will be written with the concrete edits, routed either to Codex (script-side changes) or to the paper-side workflow.
