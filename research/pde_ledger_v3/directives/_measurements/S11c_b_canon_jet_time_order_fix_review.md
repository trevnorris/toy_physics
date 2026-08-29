# canon_jet_name time-order fix directive — decision-list review (two legs, HOLD → folded → build)

**Legs (orchestrator-written ⇒ Codex + Grok):** prompt `_legs/S11c_b_canon_jet_time_order_fix_review.md`.
Raw: Codex `~/.s11_build/S11c_b_jetfix_directive_codex.txt`; Grok `~/.s11_build/S11c_b_jetfix_directive_grok.txt`.
**Both legs: diagnosis + named canonical form (`_t`/`_tt`) SOUND; HOLD on two directive defects (folded).**

## Both legs confirmed (independent hand-traces)
`canon_jet_name` (~L798-808) reduces N time tokens to a Boolean `_t`; `jet_suffix_from` correctly emits N `t`
codes; PY spells order-1 `_t`, order-2 single-token `_tt` (`S11c_b_..._sympy_audit.py:214-217,248-252`); WL
`D[…,{time,2}]` → `…_t_t`. So WL `u_1_t_t`→`u_1_t` while PY `u_1_tt` stays — asymmetric. Fix = count time tokens,
emit PY's `_tt`, recognize PY's `tt` token. S11c-a has NO multi-time jets (0 PY `_tt`, 0 WL time-order-2), so its
comparator diff is a clean byte-identical regression check.

## Blocking findings (folded)
1. **COUPLING_KERNEL is NOT tt-free** (Codex; the leg legs DISAGREED, settled by orchestrator computation).
   Grok read the post-collapse residual (0 `_tt`) and called it tt-free; Codex traced the RAW WL transcript.
   Orchestrator verification (rule 1/13):
   ```
   $ grep -oE 'Derivative\[[0-9], [0-9], [0-9], 2\]\[transverseTrialPotential[A-Za-z]+\]' \
       research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out | sort | uniq -c
     864 Derivative[0, 1, 1, 2][transverseTrialPotentialOne]
     864 Derivative[1, 0, 1, 2][transverseTrialPotentialTwo]
     864 Derivative[1, 1, 0, 2][transverseTrialPotentialThree]     # + Derivative[0,0,0,2] of eWave/uOne/Two/Three
   ```
   The WL coupling kernel DOES carry ∂²ₜ of the transverse trial; the DEFECT collapses them to `_t`, so the
   residual only LOOKS tt-free (it shows `A_T_s11cb_*_t_*`). The fix will change COUPLING_KERNEL and EXPOSE a
   genuine cross-engine asymmetry — WL keeps ∂²ₜ A_T, PY zeros `u_tt`/`e_tt` in its sector restriction
   (`S11c_b_..._sympy_audit.py:696`). My directive's "COUPLING_KERNEL tt-free/unchanged" was FALSE and would have
   told the builder to SUPPRESS a real finding. **Fix:** removed; the builder reports the pre/post diff + the
   newly-visible order-2 jets and records the asymmetry for the adjudication (rule 6). ⚠ This asymmetry is now a
   key candidate finding for Phase 2 (representational-via-dispersion `u_tt→−ω²` vs a real disagreement).
2. **Fixture must pin EXACT canonical strings** (both legs): pairwise equality passes a boolean+`tt` wrong-form
   that collapses both to `u_1_t`. **Fix:** pin `== "u_1_tt"`, `== "e_W_tt"`, `== "u_1_tt_d1"`, the coupling case
   `canon_jet_name("A_T_s11cb_1_t_t_d2_d3") == "A_T_s11cb_1_tt_d2d3"`, order-1 unchanged, spatial equal to
   PRE-FIX form (`mu_R_bg_d1_d2 → mu_R_bg_d1d2`).

## Cleared
Fix mechanism (counter-based canon_jet_name) sound; scope (time branch only) leaves spatial/order-≤1 identical;
S11c-a byte-identical regression check adequate; `t_hold_*` names protected pre-canon; no leaked inertia residual
value.

## Disposition
Both concrete corrections folded into the fix directive (one pass). Fix mechanism confirmed by both legs ⇒
release the builder (rule 7 fold-and-go); the build legs verify the actual fix + S11c-a regression + the exposed
coupling asymmetry.
