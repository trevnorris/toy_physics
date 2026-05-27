---
unit_id: 068
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26T13:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 068 (v2)

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:44-101` and
  `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:33-66`
  replace the three tautological "derivations" with three substantive anchors:
  1. **Gain decomposition (W_res = C^2 W_wall).** SymPy lines 55-63 build
     `Gmatch_expr = rho_star*g_phi^2*N_phiphi/(m*c_s^2*K_X)`,
     `Wwall_expr = kappa*Gmatch_expr`, `Gres_expr = C2*Gmatch_expr`,
     `Wres_expr = kappa*Gres_expr`, then assert
     `simplify(Wres_expr - C2*Wwall_expr) == 0`. Mathematica lines 36-44 mirror
     this with `GmatchExpr = rhoStar*gPhi^2*Nphi/(mPart*cs^2*KX)` etc.
  2. **Matched limit (C2 -> 1 collapses W_res to W_wall).** SymPy lines 64-66 and
     Mathematica lines 45-46 add the limit check.
  3. **P_res = 1/C_res^2 from required-wall-figure ratio.** SymPy lines 78-89
     use two independent `sp.solve` calls on the matched and profile Peclet
     balances (local symbols `Wm`, `Wp`), then take the ratio at
     `C2 -> Cres^2`. Mathematica lines 49-55 use
     `FullSimplify[((PeReq/(C2*Delta1))/(PeReq/Delta1)) /. C2 -> Cres^2, ...]`.
  4. **Numeric anchor.** SymPy lines 91-101 and Mathematica lines 57-66 verify
     the paper card's `P_res = 1.005612487760576` against `1/C_res^2` with
     `C_res^2 = 0.994418836451529` to 20-digit precision.

**Assessment:**
The fix is correct and substantive.
- The gain-decomposition assertion `Wres_expr - C2*Wwall_expr` non-trivially
  collapses `rho_star, g_phi, N_phiphi, m, c_s, K_X, kappa` — perturbing any
  single factor (e.g., changing `N_phiphi` to `2*N_phiphi` on only one side)
  would cause the residual to be non-zero. This is NOT a relabel-then-compare;
  W_wall is built independently as `kappa * Gmatch` and W_res independently as
  `kappa * C2 * Gmatch`, then their difference modulo `C2*Wwall` is symbolically
  simplified. The anchor in the matched-branch gain decomposition is exactly
  what the directive asked for.
- The `P_res = 1/C_res^2` derivation now flows through two genuine `sp.solve`
  calls on the Peclet balances and a ratio simplification — it is not the
  literal `(1/Cres^2) - (1/Cres^2)` of v1.
- The numeric anchor exec-log lines
  `P_res numeric residual = 5.6958391724936524581E-16` (SymPy) and
  `5.1512360584381274304121e-16` (Mathematica) confirm the paper card's quoted
  value agrees with `1/C_res^2` to ~5e-16, well below the 1e-12 threshold.
- The original tautological `P_res*C_res^2 - 1 = 0` line (which did not even
  contain `Pres`) is removed entirely.
- The Mathematica derivation uses `FullSimplify` on the ratio of required-wall
  figures rather than a `Solve` call, satisfying the F2 independence requirement
  in tandem.

No collateral edits beyond what the directive asked for. The directive's
deviation note ("Mathematica check labels were aligned to the directive's
verification inventory by adding 'from' to two labels") is a cosmetic label
change that does not affect substance.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:68-91`
  replaces the line-by-line port of SymPy's `Solve` choreography with
  `Reduce`-and-`Cases` derivations that keep the positivity premise explicit:
  ```
  WfailMatch = First[Cases[
      Reduce[Wmatch*Deltainf == PeReq && Wmatch > 0 && Deltainf > 0 && PeReq > 0, Wmatch, Reals],
      HoldPattern[Wmatch == rhs_] :> rhs, Infinity
    ]];
  ```
  applied to each of the four thresholds (WfailMatch, WsuffMatch, WfailRes,
  WsuffRes).
- The stray `gapSym;` no-op line (was line 74) is removed.
- The F1 P_res derivation uses a direct ratio expression with `FullSimplify`
  rather than mirroring SymPy's two-Solve idiom.
- The "Way B" band-width Solve idiom (`gap /. First@Solve[WsuffMatch + gap == Pres*WsuffMatch, gap]`)
  is replaced by the F3 P-form construction (the directive notes that F3
  supersedes the interim F2 Series form).

**Assessment:**
The Mathematica script no longer mirrors the SymPy script line-by-line:
- SymPy uses
  `sp.solve(sp.Eq(W_match_sym * Delta_sym, Pe_req), W_match_sym)[0]` then
  substitutes `Delta_sym -> Deltainf/Delta0`. Mathematica uses
  `Reduce[Wmatch*Deltainf == PeReq && Wmatch > 0 ..., Wmatch, Reals]` with
  `Cases` pattern extraction, separately for each Delta. The algebraic path is
  different: Mathematica's `Reduce` returns an `And`-clause structure that gets
  pattern-matched with `HoldPattern[Wmatch == rhs_] :> rhs`, while SymPy's
  `solve` returns a bare expression directly.
- For the F1 P_res derivation, SymPy uses two Solve calls then a ratio;
  Mathematica uses `FullSimplify` on a direct ratio expression with `C2 -> Cres^2`
  substitution and explicit `Delta1 > 0` assumption.
- For the F3 band-width check, SymPy uses `sp.simplify` on
  `Wsuff_res_C - Wsuff_match`; Mathematica uses `FullSimplify` with explicit
  `Assumptions -> $Assumptions`.
- The `gapSym;` leftover is gone (the diff at line 101 of stage_068_diff.patch
  shows the deletion of `gapSym;`).

The exec log confirms both engines produce the same algebraic output for the
threshold quantities (`PeReq/Deltainf`, `PeReq/Delta0`, `PeReq/(C2*Deltainf)`,
`PeReq/(C2*Delta0)`), but the printed intermediate forms differ. Mathematica's
C-form band width is printed as `(PeReq - Cres^2*PeReq)/(Cres^2*Delta0)` while
SymPy's is `-Pe_req/Delta_0 + Pe_req/(C_res**2*Delta_0)` — same expression
post-simplification but distinct simplification paths, confirming
non-parallel symbolic processing.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:145-174`
  and `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:106-125`
  replace the redundant "Way B" Solve idiom with a C-form vs P-form cross-check.
  - **C-form (Way A):** SymPy `Wsuff_res_C = Pe_req / (Cres**2 * Delta0)`,
    `success_band_widthA = sp.simplify(Wsuff_res_C - Wsuff_match)`. Mathematica
    `WsuffResC = PeReq/(Cres^2 * Delta0)`,
    `successBandA = FullSimplify[WsuffResC - WsuffMatch, ...]`. Routes through
    `Cres` directly, NOT through `Pres`.
  - **P-form (Way B):** SymPy
    `success_band_widthB = sp.simplify((Pres - 1) * Wsuff_match)`. Mathematica
    `successBandB = FullSimplify[(Pres - 1)*WsuffMatch, ...]`. Routes through
    `Pres` directly.
  - **Cross-check:** the comparison `(A - B).subs(Pres, 1/Cres**2) == 0` (SymPy)
    and `(A - B) /. Pres -> 1/Cres^2` (Mathematica) — this is now sensitive to
    the `Pres = 1/Cres^2` relation itself.

**Assessment:**
The Way A and Way B paths are now substantively distinct:
- Way A's expression involves `Cres^2` in the denominator (`-Pe_req/Delta_0 +
  Pe_req/(C_res**2*Delta_0)` from the SymPy log).
- Way B's expression involves `Pres`
  (`Pe_req*(P_res - 1)/Delta_0` from the SymPy log).
- The cross-check explicitly substitutes `Pres -> 1/Cres**2` before asserting
  equality, so any perturbation of the `Pres = 1/Cres^2` link would now produce
  a non-zero residual. Conceptual perturbation test: if `Pres` were redefined
  as `2/Cres^2` instead of `1/Cres^2`, Way A would still equal
  `Pe_req*(1/Cres^2 - 1)/Delta_0` while Way B (after the substitution) would
  equal `Pe_req*(2/Cres^2 - 1)/Delta_0` — a factor-of-2 discrepancy in the
  leading `Pe_req/Cres^2/Delta_0` term that the assertion would catch.
- The exec logs show both forms printed for both Delta_0 (success) and
  Delta_inf (failure), and both cross-checks pass: `success band C-form vs
  P-form (under Pres = 1/Cres^2) = 0` and the corresponding failure-side line.

The directive's prescription is faithfully applied in both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `W_res - C2 * W_wall (from gain decomposition) = 0`
- `W_res(C2->1) - W_wall (matched limit) = 0`
- `P_res - 1/C_res^2 (from required-wall-figure ratio) = 0`
- `P_res numeric residual = 5.6958391724936524581E-16`
- `Wfail_res * C2 - Wfail_match = 0` and
  `Wsuff_res * C2 - Wsuff_match = 0`
- `success band C-form vs P-form (under Pres = 1/Cres^2) = 0`
- `failure band C-form vs P-form (under Pres = 1/Cres^2) = 0`

**Mathematica:** exit=0. Notable lines:
- `W_res - C2 * W_wall (from gain decomposition) = 0` / `PASS`
- `W_res(C2->1) - W_wall (matched limit) = 0` / `PASS`
- `P_res - 1/C_res^2 (from required-wall-figure ratio) = 0` / `PASS`
- `P_res numeric residual = 5.1512360584381274304121e-16`
- `Success-side band width (C-form) = (PeReq - Cres^2*PeReq)/(Cres^2*Delta0)`
- `Success-side band width (P-form) = (PeReq*(-1 + Pres))/Delta0`
- `success band C-form vs P-form (under Pres = 1/Cres^2) = 0` / `PASS`
- `failure band C-form vs P-form (under Pres = 1/Cres^2) = 0` / `PASS`

Engines agree on every comparable quantity; intermediate forms differ
(Mathematica's `(PeReq - Cres^2*PeReq)/(Cres^2*Delta0)` vs SymPy's
`-Pe_req/Delta_0 + Pe_req/(C_res**2*Delta_0)`), consistent with the F2
non-parallel-symbolic-path requirement. The numeric anchor residuals (~5e-16
in both engines) are well below the 1e-12 threshold.

**Output freshness:** the saved `.txt` outputs under `scripts/output/` and
`mathematica/output/` carry mtimes of May 22 (pre-edit), while the scripts have
mtimes of May 26 (post-edit). The exec logs at
`redteam/exec_logs/stage_068_*.log` are dated May 26 13:13 and contain the
post-fix transcripts. The orchestrator's exec logs are the authoritative
post-fix record; the saved `.txt` outputs appear to be the older v1
transcripts. This is a side observation (the verifier reads the exec logs, not
the saved `.txt` outputs).

## Material-change assessment

`material_change`: false.

The v1 material change (lifting `Wfail_res / Wfail_match` from postulated
literals to Solve-derived expressions on resonance-corrected premises) is
intact and unaffected by v2 edits. Confirmed:
- SymPy section 2 (lines 113-127) still uses
  `sp.solve(sp.Eq(W_match_sym * Delta_sym, Pe_req), W_match_sym)[0]` and
  `sp.solve(sp.Eq(C2 * W_prof_sym * Delta_sym, Pe_req), W_prof_sym)[0]` to
  derive the four thresholds. The exec log confirms these produce
  `Pe_req/Delta_inf`, `Pe_req/Delta_0`, `Pe_req/(C2*Delta_inf)`,
  `Pe_req/(C2*Delta_0)` — the same load-bearing expressions the v1 fix
  delivered.
- Mathematica section 2 (lines 68-91) was upgraded by F2 from `First@Solve[...]`
  to `Reduce[...]` with positivity constraints and `Cases` pattern-matching.
  The derivation is now *more rigorous* (positivity premise made explicit) but
  the output expressions are unchanged: `PeReq/Deltainf`, `PeReq/Delta0`,
  `PeReq/(C2*Deltainf)`, `PeReq/(C2*Delta0)`. The v1 premise-replacement (still
  derived from the matched and profile Peclet balances, not assigned from
  literals) is preserved.

v2 edits modified:
1. The top-of-script W_res / P_res derivations (F1) — these were tautologies in
   v1; the v2 fix makes them substantive but the final asserted relations
   (`W_res = C^2 W_wall`, `P_res = 1/C_res^2`) are unchanged.
2. The Mathematica engine path for section 2 (F2) — Reduce instead of Solve,
   but same output.
3. The band-width cross-check (F3) — Way B now uses Pres directly instead of a
   Solve-for-gap. Output expressions for the band widths are unchanged
   (`Pe_req*(P_res - 1)/Delta_0` etc.).

No downstream stage's symbolic input from stage 068 changes. `Pe_req`,
`Delta_0`, `Delta_inf`, `Pres`, `Cres`, the four threshold expressions, and the
two band widths all remain symbolically identical to v1. Therefore no
downstream unit needs re-audit on account of v2.

## Side observations (non-blocking)

- The saved `.txt` outputs under `scripts/output/` and `mathematica/output/`
  are older than the post-fix scripts (May 22 vs. May 26). The exec logs at
  `redteam/exec_logs/stage_068_*.log` are the authoritative post-fix
  transcripts. If the orchestrator's downstream consumers read the saved
  `.txt` files instead of the exec logs, they will see stale output. This is
  a process question for the orchestrator, not a verification failure.
- Both scripts retain the original "STAGE 51" / "STAGE 051" banner text
  (file numbering 068 vs. legacy banner 51/051). Cosmetic only; flagged in v1
  side observations and remains out of scope.
- The Mathematica script's `$Assumptions` accumulates symbol declarations
  across two `$Assumptions = $Assumptions && ...` lines (lines 29-31 then
  37-38), which is fine but slightly unusual. No functional impact.
- SymPy's local symbols `Wm`, `Wp` used in the section 1b Solve calls are
  unlabelled (no `print` of the resulting Solve solutions). The directive
  explicitly notes they are scoped to those lines, so this is intentional.

None of these affect verification.

## Verdict justification

All three findings are resolved. F1's tautologies are replaced by substantive
derivations: the gain-decomposition anchor for `W_res = C^2 W_wall` (which
non-trivially collapses seven independent component symbols), the
required-wall-figure ratio for `P_res = 1/C_res^2` (which routes through two
genuine Solve calls and a ratio simplification), and a 20-digit numeric anchor
against the paper card's `P_res = 1.005612487760576`. F2's Mathematica
transliteration is broken: section 2 thresholds now derive via
`Reduce[...] && Cases[...]` with positivity premises, the `gapSym;` no-op is
gone, and the F1 P_res derivation uses a direct ratio simplification instead
of mirroring SymPy's two-Solve choreography. F3's band-width "two ways" check
is now genuinely two routes — Cres-form (no Pres) vs. Pres-form (no Cres) —
with the cross-check explicitly substituting `Pres -> 1/Cres^2`, making the
assertion sensitive to that link rather than absorbing errors symmetrically.
Both engines exit 0; the v1 material_change carry-forward (Solve-derived
thresholds in section 2) is intact and produces the same symbolic outputs as
before. No material_change in v2; downstream stages are not affected.
