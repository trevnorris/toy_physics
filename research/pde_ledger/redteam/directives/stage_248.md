---
unit_id: 248
batch: VIII.1
created_at: 2026-06-03T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-03T10:20:15-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 248

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch `paper/stages/stage_248.tex` or any prose document, EXCEPT the single notes line authorized in `## F0` below (user-approved 2026-06-03). Otherwise the red-team only modifies scripts.

## F0 — notes correction (USER-APPROVED 2026-06-03)

**Target:** `notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md:506`

**Change (apply verbatim):** `\times 168\%` → `\times 100\%` (NOTES-ONLY typo; the published card is UNAFFECTED). The notes' own next line states `\approx 23.3128\%`, and the SymPy script computes `improve_pct = 100.0*(ratio_num - 1.0)` = `23.3128` (py:238, asserted py:286) — both confirm the multiplier is `100`, not `168` (the recurring stale "168" previously corrected at stages 148/232). Edit ONLY notes:506; do NOT change the script or the card. Codex applies; Claude reviews.

## Applied: F0

- files_changed:
  - `notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md`
- summary: Corrected the notes-only transmission increase multiplier from `168\%` to `100\%`.
- deviation: none

## F1 — insufficient_verification (dynamic diagnostics identity untested)

**Target (SymPy):** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
**Target (Mathematica):** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`

**Issue:** The card boxes the diagnostic identity (eq stage248-diagnostics) `lambda_th(E) = |V(r_+(E))/V'(r_+(E))| = |E/V'(r_+(E))|` and `Xi_turn(E) = Xi_1(r_+(E))`. The scripts only assert these hardcoded benchmark numbers are positive (`assert Xi_turn_num > 0` / `assert lambda_th_num > 0` at py:290-291; `Sign[XiTurnNum]-1` / `Sign[lambdaThNum]-1` at wl:243-244). A positivity check on a constant exercises no Stage-248 relation. The load-bearing equality — that the two boxed forms of `lambda_th` agree because `V(r_+)=E` at the outer turning point — is never verified in either engine.

**Required change:**

In the SymPy script, in Section 3 (turning-point / WKB compiler, after the existing transport-law asserts at line 179, before Section 4 begins at line 184), add a non-tautological symbolic diagnostic check. Introduce an undefined function for the barrier scalar and reuse the existing outer-turning-point symbolic machinery:

```python
    # Dynamic diagnostics: lambda_th(E) = |V(r_+)/V'(r_+)| = |E/V'(r_+)|,
    # the second equality holding because V(r_+(E)) = E at the outer turning point.
    Vp_rp = sp.diff(V(r_plus_E), r_plus_E)             # V'(r_+(E_turn)), nonzero symbol
    lambda_th_def = V(r_plus_E) / Vp_rp                # |V/V'| form (drop abs for the identity)
    lambda_th_E = E_turn / Vp_rp                        # |E/V'| form
    # Pre-substitution: the two forms are NOT identically equal (non-vacuous guard).
    lambda_th_gap_raw = sp.simplify(lambda_th_def - lambda_th_E)
    # Post turning-condition V(r_+(E)) = E: they coincide.
    lambda_th_gap = sp.simplify(lambda_th_gap_raw.subs(V(r_plus_E), E_turn))
    # Xi_turn is the front-end barrier scalar sampled at the outer turning point.
    Xi1 = sp.Function("Xi1")
    Xi_turn_sym = Xi1(r_plus_E)

    print("lambda_th gap (raw)   =", lambda_th_gap_raw)
    print("lambda_th gap (V=E)   =", lambda_th_gap)
    print("Xi_turn symbolic      =", Xi_turn_sym)

    assert lambda_th_gap_raw != 0          # non-vacuous: identity is not free
    assert lambda_th_gap == 0              # identity holds under the turning condition
    assert Xi_turn_sym == Xi1(r_plus_E)    # carried tag sampled at r_+(E)
```

Keep the existing Section-5 numeric readbacks (`Xi_turn_num`, `lambda_th_num`, the `> 0` asserts) unchanged — they are the benchmark layer.

In the Mathematica script, add the matching symbolic check in Section III (after the existing transport-law `expectZero` at wl:176, before Section IV at wl:178). Use a native route, not a transliteration of the Python above:

```mathematica
(* Dynamic diagnostics: lambda_th = |V(r_+)/V'(r_+)| = |E/V'(r_+)|, the second
   equality from the outer turning condition V(r_+(E)) = E. *)
VpRp = D[Vfun[rp], rp];
lambdaThDef = Vfun[rp]/VpRp;
lambdaThE = Eturn/VpRp;
lambdaThGapRaw = FullSimplify[lambdaThDef - lambdaThE, Assumptions -> $Assumptions];
lambdaThGap = FullSimplify[(lambdaThDef - lambdaThE) /. Vfun[rp] -> Eturn, Assumptions -> $Assumptions];
Print["lambda_th gap (raw) = ", fmt[lambdaThGapRaw]];
If[TrueQ[lambdaThGapRaw == 0], fail["lambda_th identity non-vacuous guard", lambdaThGapRaw], pass["lambda_th identity non-vacuous guard"]];
expectZero["lambda_th identity under V(r_+)=E", lambdaThGap];
expectZero["Xi_turn sampled at r_+", Xi1fun[rp] - Xi1fun[rp]];
```

Note `rp = rpFun[Eturn]` is already defined at wl:153, so reuse it; `Xi1fun` is an undefined head (introduce it as a fresh symbolic head). The last `expectZero` is the definitional anchor for the carried tag (it is genuinely trivial — keep it only as a presence marker; the substantive new check is the `lambda_th` identity-under-turning-condition pair).

**Self-test (performed by auditor):** `lambda_th_def - lambda_th_E = (V(r_+) - E)/V'(r_+)`, which is nonzero in general (so `lambda_th_gap_raw != 0` holds) and reduces to 0 exactly when `V(r_+) = E` is substituted (so `lambda_th_gap == 0` holds). The check is non-vacuous and anchored to the card's boxed identity. `V` is an undefined function, so no variable-independence/identically-zero trap.

**Verification command:** verifier runs `redteam exec-sympy 248` and `redteam exec-mathematica 248`; confirms the new `lambda_th gap (V=E)` check appears and is 0, the raw-gap guard is nonzero, and both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- summary: Added symbolic dynamic-diagnostic checks for the non-vacuous `lambda_th` identity under the turning condition and the `Xi_turn` sampling marker.
- deviation: none

## F2 — mathematica_transliteration (Section II / V echo SymPy)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl:95-118`

**Issue:** Section II of the `.wl` is a line-by-line port of py:68-107: same variable choreography (`ElaunchNew`/`E_launch_new`, `vcritNew`/`vcrit_new`, `vcritNewSolved`/`vcrit_new_solved` via `Solve[...] /. ConditionalExpression` mirroring SymPy's `solve(...) if not could_extract_minus_sign`), same assertion order. This gives little independent confidence beyond re-running SymPy's algebra in Mathematica syntax.

**Required change:** Re-derive the two threshold speeds via a different native primitive so the second engine is an independent route, not an echo. Replace the `vcritNewSolved` / `vcontactCoulSolved` derivations (wl:97-100 and wl:107-110) so they use `Reduce` over the positive-root domain rather than `First[Solve[...]] /. ConditionalExpression`:

```mathematica
vcritNewSolved = FullSimplify[
  v0 /. ToRules[Reduce[ElaunchNew == Vpeak && v0 > 0 && Vpeak > V0, v0]],
  Assumptions -> $Assumptions
];
...
vcontactCoulSolved = FullSimplify[
  v0 /. ToRules[Reduce[ElaunchCoul == 1/rContact && v0 > 0 && r0 > rContact, v0]],
  Assumptions -> $Assumptions
];
```

The downstream `expectZero["solve(v_crit,new) - compiler", vcritNewSolved - vcritNew]` (wl:115) and `expectZero["solve(v_contact,coul) - compiler", vcontactCoulSolved - vcontactCoul]` (wl:117) must still pass (the `Reduce` root must equal the hand-written `Sqrt[...]` compiler form). Leave Sections I, III, IV unchanged (their hard steps are already engine-native). Leave Section V benchmark literals unchanged (a shared numeric readback is acceptable).

If `Reduce`/`ToRules` returns a form that `FullSimplify` cannot collapse against the existing `vcritNew` within the script's tolerance, append `## Blocked: F2` with the residual rather than forcing it.

**Verification command:** verifier runs `redteam exec-mathematica 248`; confirms Section II no longer uses the `First[Solve[...]] /. ConditionalExpression` pattern, the `solve(...) - compiler` checks still PASS, and the script exits 0.

## Blocked: F2

- reason: The prescribed `v0 /. ToRules[Reduce[...]]` form does not produce usable replacement rules in this Wolfram version; `Reduce` returns parameter-domain predicates conjoined with the root, `ToRules` rejects that conjunction, and after simplifying under `$Assumptions` the equation is oriented as `Sqrt[...] -> v0`, leaving the downstream residual nonzero.
- question: Should F2 allow extracting and reorienting the `v0 == root` equality from the `Reduce` result, or should a different exact `Reduce`/domain form be specified?

## F2 — ITER 2 REFRAME (orchestrator, 2026-06-03; supersedes the Blocked F2 above)

The `Reduce`/`ToRules` route is WITHDRAWN — it is a Wolfram-version idiom dead-end (confirmed by the Block above), not a math error. F2's GOAL is unchanged and still required (248 is a checkpoint → higher bar): Section II's two threshold speeds (`vcritNew`, `vcontactCoul`) must be confirmed by a route that is genuinely INDEPENDENT of the SymPy `Solve(...)`-then-compare that the `.wl` currently mirrors (the `First[Solve[...]] /. ConditionalExpression` + minus-sign branch logic).

**Requirement (Codex designs the implementation — this is acceptance criteria, NOT a hand-coded route):** Re-confirm both threshold speeds WITHOUT re-running `Solve`/`Reduce` to re-derive them. Use the independent **satisfaction** primitive: show each hand-written compiler closed form SATISFIES its own defining energy equality under the script's `$Assumptions`, and is the physical positive root. (If you prefer a different genuinely-native route — e.g. a correctly-extracted `Reduce` root via `Cases[..., HoldPattern[v0 == r_] :> r]` — that is also acceptable, as long as it is NOT the `Solve`/`ConditionalExpression` mirror.)

**Acceptance criteria (the verifier will check all five):**
1. The `vcritNewSolved` / `vcontactCoulSolved` derivations using `First[Solve[...]] /. ConditionalExpression` (and the withdrawn `Reduce`/`ToRules` attempt) are GONE from Section II.
2. In their place, a native substitution check that each compiler form is a root of its defining equation — i.e. `FullSimplify[(<launch-energy expression in v0> /. v0 -> vcritNew) - <peak/contact target>] === 0` under `$Assumptions`, and likewise `vcontactCoul` against `1/rContact`.
3. A NON-VACUITY guard: the pre-substitution gap `(<launch-energy in v0>) - <target>` is NOT identically zero (it genuinely depends on `v0`), so the satisfaction check would fail for a wrong closed form.
4. A positive-branch guard: `vcritNew > 0` and `vcontactCoul > 0` hold under `$Assumptions` (the verified root is the physical one, matching SymPy's branch selection).
5. The new check is engine-native (substitution + `FullSimplify`), NOT a re-run of `Solve`; Sections I/III/IV/V stay unchanged; both engines exit 0.

Apply on iter 2, then append `## Applied: F2` (or `## Blocked: F2` with the residual if even the satisfaction route resists `FullSimplify`).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- summary: Replaced the solver-derived threshold-speed witnesses with native substitution, non-vacuity, and positive-branch checks for the hand-written compiler forms.
- deviation: Used the same local physical domain inequalities from the withdrawn positive-root request (`Vpeak > V0`, `r0 > rContact`) because the global `$Assumptions` do not imply the square-root branches are positive.

## F3 — stale banner label (STAGE 231 → STAGE 248)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl:65`

**Issue:** The banner prints `"STAGE 231 — DYNAMIC EVENT-CHAIN COMPILER FROM THE RELAXED STATIONARY BARRIER FRONT END"` for a stage-248 file (leftover from a clone). The captured output header also reads "STAGE 231".

**Required change:** Change the banner string at wl:65 from `"STAGE 231 — DYNAMIC EVENT-CHAIN COMPILER FROM THE RELAXED STATIONARY BARRIER FRONT END"` to `"STAGE 248 — DYNAMIC EVENT-CHAIN COMPILER FROM THE RELAXED STATIONARY BARRIER FRONT END"`.

**Verification command:** verifier runs `redteam exec-mathematica 248`; confirms the output header reads "STAGE 248" and the script exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- summary: Updated the Mathematica audit banner from Stage 231 to Stage 248.
- deviation: none
