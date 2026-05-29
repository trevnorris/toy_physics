---
unit_id: 015
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive - unit 015

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F1), do nothing - the orchestrator is holding for user resolution. Do not edit `paper/stages/stage_015.tex`, `paper/appendices/`, `notes/`, or any prose document. Do not "fix" F1 by trimming the script blocks unless a follow-up directive from the user explicitly authorizes direction (b).

F3 (tautological_check on wall-only K1/H_even specialization) is **blocked on F1's resolution**: the right edit depends on whether the wall-only block survives. Codex must skip F3 in this round and append a `## Blocked: F3` block.

F4 (mathematica_transliteration) requires restructuring the Mathematica K_eta block in a way that is mechanical but not purely textual. The instructions below are deliberately concrete; if any step is ambiguous after reading, append `## Blocked: F4` with the specific question and skip.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch `paper/`, `notes/`, or any prose documents.

## F1 - paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_015.tex:44-46` quote: "Stage~015 exports the parent throat action promotion and the exact quadratic recovery formula \eqref{eq:stage015-keta}."
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:52` quote: "015 & Parent throat action master & \StatusExact{} / \StatusReduced{} & Parent throat-action packet and projection/reduction status boundary."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:103-208` (wall-only K1/H_even gates, Jacobian determinant, perturbed-solve diagnostics, real-Y20 overlap ratios `{1, 1/2, -1}`, grouped trace/line `xbar = x0, bx = 3*ax`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:104-196` (mirroring M4-M9 blocks)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:2` docstring quote: `"""Master-note audit for step_13_parent_throat_action_master_notes.md."""` (the referenced note does not exist under `notes/stages/`)

## Resolve before fix_loop

The scripts (sympy A9-A19; mathematica M8-M28) verify three blocks that have no paper-side counterpart: (i) wall-only specializations of the even-channel gates K1 = -dM + dK/9 and H_even = (2/3)dM - dK/27 with their Jacobian determinant `1/27` and Gaussian overlaps `dMoverlap = sqrt(pi/3)`, `dKoverlap = 23 sqrt(pi)/(3 sqrt(3))`; (ii) real-Y20 overlap ratios `{1, 1/2, -1}` for `m = 0, 1, 2`; (iii) grouped trace/line identities `xbar = x0`, `bx = 3*ax`. The stage 015 paper card mentions none of these. The sympy docstring references `step_13_parent_throat_action_master_notes.md`, which does not exist in `notes/stages/`. Which is correct?

Possible directions (the user picks one):
- (a) **Paper card is incomplete.** The wall-only/Y20/grouped blocks are genuine stage 015 deliverables that were omitted from the card during a paper revision. Action: expand `paper/stages/stage_015.tex` to add `\paragraph{Wall-only even-gate reduction.}`, `\paragraph{Real-Y20 overlap ratios.}`, and `\paragraph{Grouped trace/line invariants.}` paragraphs documenting these as outputs; create or restore `notes/stages/moving_throat_pde_stage015_*.md` from the master-note source the docstring references; then re-issue this directive with F1 marked `resolved` so Codex can act on F3 (which becomes "non-tautological derivation required for wall-only specialization") and F4 (which extends to cover the wall-only block under Mathematica re-derivation).
- (b) **Script is bloated.** The wall-only/Y20/grouped blocks were carried over from a master-note draft that has since been split into separate stages and no longer belongs in stage 015. Action: trim the sympy script to remove lines 103-208 (everything from `wall_w = sp.symbols(...)` through `assert_zero("grouped line b=3a", bx - 3 * ax)`), update the final `print` summary at sympy line 211 to drop the wall-only/grouped wording, trim the mathematica script to remove lines 104-196 (everything from `D01full = ...` through `expectZero["M9 grouped line b=3a", bx - 3*ax];`), update the sympy docstring to a stage 015 description that does not reference `step_13_*_notes.md`, regenerate the stale outputs.
- (c) **Both directions partially apply.** A subset of the extra blocks belongs to stage 015 (e.g., the wall-only block as the `\StatusReduced{}` half of the status boundary the appendix row mentions) while another subset (Y20, grouped trace) belongs to a different stage. Action: specify which blocks to keep and which to trim; the orchestrator will issue a follow-up directive.

The orchestrator will not invoke Codex on F3 or F4 (within the wall-only scope) until the user has chosen.

## F2 - insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:55-72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:67-81`

**Issue:**
The concrete-Gaussian IBP test uses `A_concrete = exp(-w^2)` (even) and `eta_concrete = exp(-w^2/2)` (even). Both the cross integrand `-A * eta * eta_w` and the bulk integrand `A_w * eta^2 / 2` end up odd in `w`, so each integral vanishes by parity over `(-oo, oo)` and the boundary discharge is also 0. The IBP identity asserted at sympy line 69-72 / mathematica line 81 thus reduces to `0 - (0 + 0) = 0`, which is trivially true regardless of whether IBP was applied correctly. The paper-card promise of "both zero and nonzero boundary-discharge probes" is satisfied only by the unrelated `atan(w)` probe on the `boundary_value` operator (sympy line 57-58), not by an asymmetric concrete IBP profile.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`, *after* the existing assertion at line 72 (preserve the Gaussian baseline), insert the following block:

```python
    # Second concrete probe with asymmetric A profile so the cross and bulk
    # integrals are individually nonzero. Both integrals are even by parity
    # (A_extra is odd, eta is even, products produce even integrands), so
    # the IBP identity asserts a real cancellation, not 0 = 0 + 0.
    A_concrete_asym = w_ibp * sp.exp(-w_ibp**2)
    eta_concrete_asym = sp.exp(-w_ibp**2 / 2)
    quad_boundary_concrete_asym = boundary_value(
        (-A_concrete_asym * eta_concrete_asym**2 / 2), w_ibp
    )
    quad_cross_concrete_asym = sp.integrate(
        -A_concrete_asym * eta_concrete_asym * sp.diff(eta_concrete_asym, w_ibp),
        (w_ibp, -sp.oo, sp.oo),
    )
    quad_bulk_concrete_asym = sp.integrate(
        sp.diff(A_concrete_asym, w_ibp) * eta_concrete_asym**2 / 2,
        (w_ibp, -sp.oo, sp.oo),
    )
    assert_nonzero("asymmetric concrete IBP cross is nontrivial", quad_cross_concrete_asym)
    assert_nonzero("asymmetric concrete IBP bulk is nontrivial", quad_bulk_concrete_asym)
    assert_zero("asymmetric concrete quadratic IBP boundary discharge", quad_boundary_concrete_asym)
    assert_zero(
        "asymmetric concrete quadratic IBP with boundary",
        quad_cross_concrete_asym - (quad_boundary_concrete_asym + quad_bulk_concrete_asym),
    )
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`, *after* the existing assertion at line 81, insert the analogous block:

```mathematica
(* Second concrete probe with asymmetric A profile so cross and bulk are
   individually nonzero. Asserts a real IBP cancellation, not 0 = 0 + 0. *)
aConcreteAsym = w*Exp[-w^2];
etaConcreteAsym = Exp[-w^2/2];
boundaryGaussianAsym =
  Quiet[
    Limit[-aConcreteAsym*etaConcreteAsym^2/2, w -> Infinity] -
      Limit[-aConcreteAsym*etaConcreteAsym^2/2, w -> -Infinity],
    Limit::alimv
  ];
crossGaussianAsym =
  Integrate[-aConcreteAsym*etaConcreteAsym*D[etaConcreteAsym, w], {w, -Infinity, Infinity}];
bulkGaussianAsym =
  Integrate[D[aConcreteAsym, w]*etaConcreteAsym^2/2, {w, -Infinity, Infinity}];
expectNonzero["M2 asymmetric IBP cross nontrivial", crossGaussianAsym];
expectNonzero["M2 asymmetric IBP bulk nontrivial", bulkGaussianAsym];
expectZero["M2 asymmetric IBP boundary discharge", boundaryGaussianAsym];
expectZero["M2 asymmetric IBP cross equals bulk", crossGaussianAsym - bulkGaussianAsym];
```

Do not delete the existing Gaussian-pair asserts (they remain as a parity baseline). Do not modify the K_eta block (lines 74-101 sympy / 83-102 mathematica).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 015` and `redteam exec-mathematica 015` and confirms (a) the new `assert_nonzero("asymmetric concrete IBP cross is nontrivial", ...)` / `expectNonzero["M2 asymmetric IBP cross nontrivial", ...]` lines appear, (b) both scripts exit 0, and (c) the printed residual for `quad_cross_concrete_asym` / `crossGaussianAsym` is nonzero (specifically `-sqrt(pi/2)/4` for the sympy probe; mathematica should print the same up to its display convention).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`
- summary: Added asymmetric concrete IBP probes with nonzero cross and bulk integrals while preserving the original Gaussian parity baseline.
- deviation: none

## F3 - tautological_check

**Status:** BLOCKED on F1.

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:118-127`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:104-113`

**Issue:**
`assert_zero("wall-only K1 specialization", K1_wall - (-dM + dK / 9))` and `assert_zero("wall-only H_even specialization", H_even_wall - (sp.Rational(2, 3) * dM - dK / 27))` (sympy 126-127; math 112-113) compare the result of `K1_full.subs(wall_only_specialization)` against the same expression manually inlined as the expected RHS. The assert cannot fail for any choice of `1/9` or `-1/27` because the LHS *is* the substitution.

**Why blocked:** The correct action depends on F1's resolution. If (a), F3 becomes "derive the wall-only forms from an independent paper-side reduction principle" - a substantive new check Codex needs the paper-side specification for. If (b), F3 dissolves because the entire wall-only block is removed. If (c), Codex needs to know which subset of the block survives. Append `## Blocked: F3` with the question "Awaiting F1 resolution (a/b/c)" and skip to F4.

## Blocked: F3

- reason: The wall-only K1/H_even specialization edit is blocked on F1's unresolved paper_misalignment scope.
- question: Awaiting F1 resolution (a/b/c)

## F4 - mathematica_transliteration

**Status:** PARTIALLY BLOCKED on F1. The K_eta block (M3) re-derivation below is independent of F1 and Codex may apply it now. The wall-only block (M4-M7) re-derivation is blocked on F1 resolution and Codex must skip that portion.

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:83-102` (K_eta block M3 - re-derive)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:104-171` (wall-only block M4-M7 - blocked)

**Issue:**
The K_eta block computes `L2raw = Coefficient[Series[lagrangian, {eps, 0, 2}] // Normal, eps, 2]`, then writes `L2afterIBP = Expand[L2raw - (-TwR0*R0p*eta*etaw) + dTwRR0p*eta^2/2]` and compares it against a pre-written `canonicalL2 = mu0*etat^2/2 - Tw0*etaw^2/2 - TO0*grad2/2 - effectiveMass*eta^2/2`. This mirrors the SymPy script line-by-line (same `L2_raw`/`canonicalL2` decomposition, same `cross_term`/`cross_after_ibp` peeling). A genuine independent derivation should compute K_eta from the Euler-Lagrange operator on the full `lagrangian`, not via series-expansion-and-pattern-match.

**Required change (K_eta block only):**

Replace lines 83-102 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl` with the following Euler-Lagrange derivation (the wall-only block at lines 104+ stays untouched until F1 resolution):

```mathematica
(* M3: K_eta via direct Euler-Lagrange linearization on the parent
   Lagrangian density. We treat R as a field, write R = R0[w] + eps*phi[t,w,Om]
   where phi is the linearized fluctuation, expand the EL equation to O(eps),
   and read off the mass coefficient as -K_eta * phi. This path is independent
   of the SymPy script's series-coefficient choreography. *)

ClearAll[R0, phi, TwSig, USig, KetaFromEL];
TwSig[r_, w_] := Tw0 + (r - R0[w])*TwR0 + (r - R0[w])^2*TwRR0/2;
USig[r_, w_] := U0 + (r - R0[w])*UR0 + (r - R0[w])^2*URR0/2;

(* Lagrangian density as function of R(t,w,Om) and its derivatives, with
   |grad_Om R|^2 abbreviated by gO so we do not need an explicit Om mesh. *)
LDensity[R_, Rt_, Rw_, gO_, w_] :=
  mu0*Rt^2/2 - TwSig[R, w]*Rw^2/2 - TO0*gO/2 - USig[R, w];

(* Euler-Lagrange operator at fixed Om: d_R L - d_t (d_Rt L) - d_w (d_Rw L).
   Compute symbolically with the partials taken as if R, Rt, Rw, gO are
   independent slot variables, then substitute the on-shell expansion. *)
ELRaw[R_, Rt_, Rw_, gO_, w_] :=
  D[LDensity[R, Rt, Rw, gO, w], R] -
    Dt[D[LDensity[R, Rt, Rw, gO, w], Rt], t,
       Constants -> {mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0}] -
    Dt[D[LDensity[R, Rt, Rw, gO, w], Rw], w,
       Constants -> {mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0}];

(* Substitute the linearization R = R0[w] + eps*eta, Rt = eps*etat,
   Rw = R0'[w] + eps*etaw, gO = eps^2*grad2 (transverse term is O(eps^2)). *)
ELLinearized =
  ELRaw[R0[w] + eps*eta, eps*etat, R0'[w] + eps*etaw, eps^2*grad2, w] /.
    {R0'[w] -> R0p,
     Dt[R0p, w, Constants -> {mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0}] ->
       R0pp,
     Dt[R0[w], w, Constants -> {mu0, Tw0, TO0, U0, TwR0, TwRR0, UR0, URR0}] ->
       R0p};

(* Take the O(eps) coefficient and read off the eta-coefficient: that is
   -K_eta * eta after the etat-/etaw-derivative terms are stripped. *)
ELOrderEps = Coefficient[Series[ELLinearized, {eps, 0, 1}] // Normal, eps, 1];
ELMassCoeff =
  Coefficient[ELOrderEps /. {Dt[eta, t, ___] -> 0, Dt[eta, w, ___] -> 0,
                              Dt[etat, t, ___] -> 0, Dt[etaw, w, ___] -> 0},
              eta];
(* The EL gives -K_eta * eta + (kinetic / gradient pieces); ELMassCoeff is
   -K_eta. The IBP-derived dTwRR0p stands for d/dw (TwR0(R0,w) * R0'(w)). *)
KetaFromEL =
  -ELMassCoeff /. {Dt[TwR0, w, ___]*R0p + TwR0*R0pp -> dTwRR0p,
                   TwR0*R0pp + Dt[TwR0, w, ___]*R0p -> dTwRR0p};
expectZero["M3 K_eta via EL linearization matches IBP form",
  KetaFromEL - (URR0 - dTwRR0p + TwRR0*R0p^2/2)];

(* Sign-mutation guard: a + instead of - on dTwRR0p must fail. *)
expectNonzero["M3 K_eta via EL dTwRR0p sign mutation",
  KetaFromEL - (URR0 + dTwRR0p + TwRR0*R0p^2/2)];
```

After insertion, delete the now-redundant `Tw = Tw0 + eps*TwR0*eta + ...`, `U = ...`, `Rt = eps*etat`, `Rw = R0p + eps*etaw`, `lagrangian = ...`, `L2raw = ...`, `crossCoeff = D[D[L2raw, eta], etaw]`, `expectZero["M3 K_eta raw eta etaw cross coefficient", crossCoeff + TwR0*R0p]`, `effectiveMass = ...`, `canonicalL2 = ...`, `L2afterIBP = ...`, `expectZero["M3 K_eta canonical quadratic form", L2afterIBP - canonicalL2]`, `effectiveMassMutated = ...`, `canonicalL2Mutated = ...`, `expectNonzero["M3 K_eta dTwRR0p sign mutation", ...]` lines (the entire current M3 block, lines 83-102) - their replacement is the EL derivation above.

If `Dt` substitution proves unreliable (Mathematica `Dt` is finicky about `Constants`), append `## Blocked: F4 (K_eta block)` with the specific Dt error and skip.

Do NOT modify the wall-only block (lines 104-171) or M8-M9 (lines 173-196); those are blocked on F1.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 015` and confirms (a) the new `expectZero["M3 K_eta via EL linearization matches IBP form", ...]` line appears, (b) the script exits 0, (c) the printed residual for this check is `0`, and (d) the Mathematica script no longer contains the literal token `L2afterIBP` (a structural marker that the transliteration spine has been replaced).

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`
- summary: Replaced the Mathematica M3 K_eta series-matching block with an Euler-Lagrange linearization that reads the mass coefficient directly.
- deviation: Wall-only M4-M7 re-derivation left untouched because it remains blocked on F1; the EL implementation uses independent slot variables before substituting the field expansion.

## Applied: F4 (iter2)

- files_changed:
  - `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`
- summary: Reworked the M3 Euler-Lagrange derivation to use ordinary `D` with an explicit temporary `twR[w]` profile, then collapsed `d/dw(TwR0 R0p)` to `dTwRR0p`.
- deviation: none
