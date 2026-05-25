---
unit_id: 010
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 2
findings_total: 2
material_change: false
paper_alignment: aligned
---

# Verification — unit 010 (v2 paper-alignment iteration)

This verification covers the v2 paper-grounded audit pass, which produced two
`paper_misalignment` findings (F1 = `δu_n` vs `δP_n` mismatch; F2 = seven script
clusters with no paper anchor). The user picked direction (c) for F1 ("both
sides grow") and direction (a) for F2 ("paper card is incomplete"); resolutions
are recorded in `redteam/resolutions/batch_I1_paper_alignment.md` (Q5, Q6).
An earlier verification of the v1 (scripts-only) iteration is preserved by git
history; this file replaces it with the v2 outcome.

## Per-finding outcomes

### F1 — paper_misalignment (δu_n vs δP_n)

**Classification:** resolved

**What changed:**

Paper side (`paper/stages/stage_010.tex`):
- New `\paragraph{Bundle shifts.}` block (`:13-58`) keeps the existing
  `\eqref{eq:stage010-du2}`/`\eqref{eq:stage010-du4}`/`\eqref{eq:stage010-dp0}`
  identities, then adds explicit display equations for the full prefactor
  variations under labels `eq:stage010-dP2` (lines 42-48) and `eq:stage010-dP4`
  (lines 49-58), with the bridge text "The denominator-only quantities
  \(u_2,u_4\) and the full prefactor quantities \(P_n\) are distinct."
- `\stagefield{Output}` paragraph (lines 137-144) cites
  `\eqref{eq:stage010-projected-shifts}--\eqref{eq:stage010-dP4}`, picking up
  the new dP4 anchor as the end of the transport-map range.
- Local label rename `eq:stage006-projected-shifts → eq:stage010-projected-shifts`
  (line 16). `rg` over the source tree shows the new label is referenced in
  `paper/stages/stage_010.tex:139` only; remaining `eq:stage006-projected-shifts`
  hits are in stale `.aux` build artifacts (refreshed at next LaTeX build) and
  in historical audit prose under `redteam/`.

Script side — SymPy (`scripts/moving_throat_pde_stage010_*_sympy_audit.py`):
- Lines 50-51 add the denominator-only inversion coefficients
  `u2p = -D2p/D0p`, `u4p = D2p**2/D0p**2 - D4p/D0p` built from the perturbed
  `D_n^p = D_n - eps z_n` quantities.
- Lines 56-57 compute `du2 = sp.diff(u2p, eps).subs(eps, 0)` and
  `du4 = sp.diff(u4p, eps).subs(eps, 0)`.
- Lines 62-66 assert the paper's closed forms:
  `assert_zero("delta u2", du2 - (D0*z2 - D2*z0)/D0**2)` and
  `assert_zero("delta u4", du4 - (D0**2*z4 - D0*(2*D2*z2 + D4*z0) + 2*D2**2*z0)/D0**3)`.

Script side — Mathematica (`mathematica/moving_throat_pde_stage010_*_mathematica_audit.wl`):
- Lines 31-32 add `u2slot[e_] := -den2[e]/den0[e]` and
  `u4slot[e_] := den2[e]^2/den0[e]^2 - den4[e]/den0[e]`.
- Lines 41-42 extract the first-order coefficient via
  `Coefficient[Normal[Series[u2slot[eps], {eps, 0, 1}]], eps, 1]` (and the same
  for `u4slot`).
- Lines 47-49 check `m0aResidual = FullSimplify[u2Linear - (D0 z2 - D2 z0)/D0^2]`
  with `Exit[1]` on `=!= 0`; lines 51-57 mirror for `m0bResidual`.

**Assessment:**

The new assertions are non-tautological. `u2p` and `u4p` are defined purely in
terms of the perturbed `D_n^p` and then differentiated by sympy; the RHS in the
`assert_zero` is the paper's symbolic closed form in `(D0, D2, D4, z0, z2, z4)`.
A wrong sign or coefficient on the paper-side RHS would leave a nonzero residual
after `sp.simplify`. The Mathematica mirror uses a different idiom
(`Series`+`Coefficient` rather than `sp.diff(...).subs(eps, 0)`), so the two
engines reach the same RHS through structurally distinct algebra paths.

The existing `δP_0`, `δP_2`, `δP_4` assertions (sympy `:67-92`,
mathematica `:59-79`) are untouched — they continue to verify the same closed
forms they did in the v1 iteration. The paper additions for `eq:stage010-dP2`
and `eq:stage010-dP4` reproduce those same RHS values exactly, so paper and
scripts now display the same `δP_n` closed forms.

Hand check on a couple of the `δP_4` cross terms against the paper's display:
paper line 54-55 shows `(N_4 z_0 + 2 N_2 z_2 - 2 D_2 n_2 + 2 N_0 z_4 - 2 D_4 n_0)/D_0^2`
matching script `:81-86` term-for-term; paper line 55 `(3 D_2^2 n_0 - 4(D_2 N_2 + D_4 N_0) z_0 - 6 D_2 N_0 z_2)/D_0^3` matches script `:87-89`; paper line 56
`9 D_2^2 N_0 z_0/D_0^4` matches script `:90`. Engines agree.

### F2 — paper_misalignment (7 unanchored script clusters)

**Classification:** resolved

**What changed:**

Paper side (`paper/stages/stage_010.tex`):
- `\paragraph{Compatibility transport.}` (lines 60-94) introduces `S`, `T`,
  `P_{0,target}`, the one-pole and fixed-target normalization K-surfaces
  (label `eq:stage010-k-surfaces`, lines 63-73), the fixed-target compatibility
  variation `δ𝒞_fixed` with explicit `z_0`-cancellation note
  (label `eq:stage010-compat-fixed`, lines 75-83), and the transported-target
  variant `𝒞_tr`/`δ𝒞_tr` (label `eq:stage010-compat-transport`, lines 86-94).
- `\paragraph{Grouped weak-axisymmetric lane.}` (lines 96-112) introduces
  the real-`Y_{20}` square-overlap lane multipliers `(λ_{20}, λ_{21}, λ_{22})
  = (1, 1/2, -1)` (label `eq:stage010-y20-lambdas`) and the trace/anomaly
  decomposition `b_x = 3 a_x` (label `eq:stage010-weak-axisym-trace`).
- `\paragraph{Primitive static prefactor.}` (lines 114-126) introduces the
  primitive mouth data `Q, Δ, P, q_1, d_1, p_1`, the perturbation slopes
  `z_0^{prim}, n_0^{prim}`, and the `Ξ_static` closed form
  (label `eq:stage010-primitive-xi`).
- Sign-flip mutation guards equation (label `eq:stage010-mutation-guards`,
  lines 128-133) records the nonzero residual `6 S^2 z_4/T^2` left by flipping
  the sign of the `3 S^2 z_4/T^2` term in either compatibility variation.
- `\stagefield{Output}` paragraph (lines 137-144) enumerates the new anchors:
  transport map, compatibility transports, weak-axisymmetric lane signature,
  primitive static prefactor anchor, and sign-flip guards.

Paper appendix row updated (`paper/appendices/stage_appendix_part01.tex:42`):
"Projected \(Z_n,N_n\) slot transport, denominator and prefactor variations,
compatibility transports, and weak-axisymmetric lane anchors." (previously:
"Projected \(Z_n,N_n\) slot transport for the grouped response bundle.")

Script side (housekeeping only — no algebra changed):
- SymPy docstring at `:2` now reads "Stage 010 audit for projected-Maxwell
  transport into grouped bundle slots." (previously referenced
  `step_08_projected_maxwell_push_bundle_master_notes.md`, an unsurvived
  EM-projected notes file). Banner print at `:180` updated to
  "STAGE 010 PROJECTED MAXWELL PUSH MASTER AUDIT" (previously "STEP 08").

**Assessment:**

The seven previously-unanchored clusters now have explicit equation labels in
the stage card, and the Output paragraph cites them. The paper text reproduces
the script's identities verbatim where they are displayed, and labels them so
downstream stages (011 P_2 bridge, 012 primitive bridge, 013-014 mouth-Taylor)
can cite them by reference.

No script-side algebra changed for F2: the seven clusters were already
script-verified; the resolution was to grow paper anchors so the audit body is
no longer orphan. The docstring/banner cleanup removes the dangling reference
to the missing `step_08_*_notes.md` file but does not change any assertion.

The label rename `eq:stage006-projected-shifts → eq:stage010-projected-shifts`
is correct: the apply notes flagged that `rg` found no external refs to the
old name, and my own search reproduces this (only stale `.aux` artifacts and
historical audit prose under `redteam/`).

## Exec log assessment

**SymPy:** exit=n/a. Per-stage exec log `redteam/exec_logs/stage_010_sympy.log`
is **not present** for this iteration; only `stage_010_diff.patch` (mtime
2026-05-21) survives from v1. Saved transcript at
`scripts/output/moving_throat_pde_stage010_*_sympy_audit.txt` is also stale
(mtime 2026-05-21 11:51, ~4 days older than the current script mtime
2026-05-25 02:19). The verifier sandbox blocks script execution, so I cannot
refresh either artifact. The user-supplied directive states "Codex apply notes:
... scripts/.../stage010_*_sympy_audit.py: added u2/u4 definitions and
assertions (lines 50-65)" and the resolution-file `sanity_check` line reads
"sympy + mathematica exit 0"; treating that user-confirmed orchestrator
report as authoritative for the missing log.

**Mathematica:** exit=n/a. Same situation —
`redteam/exec_logs/stage_010_mathematica.log` is not present; saved transcript
at `mathematica/output/moving_throat_pde_stage010_*_mathematica_audit.txt`
also predates the current script. User-supplied directive confirms both
engines exit 0 post-edit.

**Output freshness:** the saved `.txt` transcripts are NOT post-edit fresh
(script mtimes 2026-05-25 02:15/02:19, transcript mtimes 2026-05-21 11:51).
This is documented as a side observation; it does not block verification
because (a) the orchestrator's user-relayed report confirms exit 0 from both
engines after the v2 edits, and (b) the script-side edits I can read are
non-tautological by construction (see F1 assessment). If the orchestrator
wants to backfill the per-stage exec logs for the v2 iteration, that is a
log-collection concern, not a substance concern.

## Material-change assessment

`material_change`: **false**.

The δu_n forms are *new* closed forms not previously asserted in either engine,
but they do not redefine or update any value that prior verification or any
downstream stage already depended on. The pre-existing δP_0, δP_2, δP_4
closed forms are byte-identical to the v1 iteration. The seven extra clusters
(K-surfaces, transported-target compatibility, Y20 lane multipliers, weak-
axisymmetric trace anomaly, primitive static Ξ, two mutation guards) were
already script-verified in v1; the v2 edits only added paper-side anchors and
display equations, not new derived values. The label rename is internal to
stage 010 (no external refs).

No downstream unit needs to be marked `upstream_stale` on substance grounds.

## Side observations (non-blocking)

- Per-stage exec logs `redteam/exec_logs/stage_010_sympy.log` and
  `redteam/exec_logs/stage_010_mathematica.log` are absent for this v2
  iteration. Saved `output/*.txt` transcripts also predate the v2 script
  edits (2026-05-21 vs 2026-05-25). Orchestrator's log-collection step
  appears to have skipped this stage for the v2 paper-alignment pass.
  User-supplied report of "exit 0" is being accepted as authoritative.
- Stale `eq:stage006-projected-shifts` references remain in
  `paper/pde_ledger.aux` (build artifact, will refresh on next LaTeX build)
  and in `redteam/{directives,reports}/stage_010.md` (historical audit prose
  that quotes the old label; no need to backfill).
- The previous verification file content for this stage was for the v1
  (scripts-only) iteration; this file replaces it. Git history preserves
  the v1 verification at the prior blob.

## Verdict justification

Both v2 `paper_misalignment` findings are resolved by user-approved Codex
edits. F1 grew both the paper (new `eq:stage010-dP2`, `eq:stage010-dP4`
display equations and explanatory bridge text) and the scripts (new `δu_2`
and `δu_4` assertions in both engines via independent algebra paths —
`sp.diff` for SymPy, `Series`+`Coefficient` for Mathematica). F2 grew the
paper card with explicit equation labels for all seven previously-orphan
clusters and updated the Output paragraph to cite them; the script-side
docstring/banner cleanup removes the dangling reference to the missing
EM-projected notes file. Paper alignment status moves from `partial` to
`aligned`. No previously-verified value changed; `material_change: false`.
Exec logs are not present for this iteration (documented as a side
observation); accepting user-supplied "exit 0" report as authoritative.

stage 010: verified
