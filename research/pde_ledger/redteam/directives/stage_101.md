---
unit_id: 101
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 101

Apply each non-`paper_misalignment` finding below (F1, F2, F3) in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For F4 (`paper_misalignment`), DO NOTHING. The orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or the scripts to "fix" F4 unless a follow-up directive explicitly authorizes a direction.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

---

## F1 — tautological_check (Mathematica)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:40`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:41`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:47`

**Issue:** Each `expectZero` checks `(nQExact /. {...}) - <expected>`, where `<expected>` is exactly what `nQExact` literally reduces to after the substitution. The residual is `0` by construction of `nQExact` from `Solve`; no physics is exercised.

**Required change:**

Replace the three `expectZero` calls so the residual is anchored to the **input factorized equation** `mHat0^2*chiQ*nQ - 1`, not to `nQExact` substitutions. Concretely:

1. At line 40, replace the existing call with:
   ```
   expectZero["point-particle natural branch reduction", (mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, nQ -> 1/chiQ}];
   ```
2. At line 41, replace the existing call with:
   ```
   expectZero["canonical compact outgoing branch gives NQ=1", (mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, chiQ -> 1, nQ -> 1}];
   ```
3. At line 47, replace the existing call with:
   ```
   expectZero["exact replacement chiQ=1+DeltaQ", FullSimplify[(mHat0^2*chiQ*nQ - 1) /. {mHat0 -> 1, chiQ -> 1 + deltaQ, nQ -> 1/(1 + deltaQ)}, Assumptions -> $Assumptions]];
   ```

Do **not** delete the lines defining `nQSol`/`nQExact` (lines 33-34) or the three diagnostic `Print[...]` lines (36-38) — those remain useful as transcripts of the solved form.

**Why this is non-tautological:** the residual on the LHS is now the *input* equation `mhat0^2 * chiQ * NQ - 1`, and the substituted RHS uses a *proposed* candidate value of `NQ`. If the proposed candidate is wrong (typo, sign flip, factor error), the residual no longer simplifies to `0`. By contrast, the previous form substituted `nQExact` (the solved value) into itself, which trivially gives `0`.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 101` and confirms (a) the `.wl` still exits 0 with three `PASS` lines, and (b) each `expectZero` argument now mentions `mHat0^2*chiQ*nQ - 1` (not `nQExact`).

---

## F2 — missing_verification_script (subtype: script_doesnt_cover_claim — SymPy)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py:1-18`

**Issue:** SymPy script contains zero `assert` statements. Every numerical/symbolic step is `print(...)`-only, so the script's "PASS" exit code is decoupled from the correctness of the math. A wrong proposed value of `NQ` would still produce exit 0.

**Required change:**

Insert four `assert` statements. After line 11 (`print('canonical compact outgoing branch chiQ=1 gives = ...')`), insert:

```python
# Anchor the reductions to the INPUT factorization mhat0^2 * chiQ * NQ = 1.
# Substituting the proposed solved NQ on the natural source-map branch must zero the residual.
assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, NQ: 1/chiQ})) == 0, \
    'point-particle natural branch reduction NQ = 1/chiQ failed against input factorization'
assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, chiQ: 1, NQ: 1})) == 0, \
    'canonical compact outgoing branch NQ = 1 failed against input factorization'
```

After line 16 (`print('small-DeltaQ series = ...')`), insert:

```python
# Anchor exact NQ = 1/(1 + DeltaQ) on the natural source-map branch to the input factorization,
# with chiQ = 1 + DeltaQ.
assert sp.simplify((mhat0**2 * chiQ * NQ - 1).subs({mhat0: 1, chiQ: 1 + DeltaQ, NQ: 1/(1 + DeltaQ)})) == 0, \
    'exact replacement NQ = 1/(1+DeltaQ) failed against input factorization'
# Confirm the small-DeltaQ linearization matches the paper's stated form NQ - 1 = -DeltaQ + O(DeltaQ^2).
# The script's series_delta is the order-2 truncation of -DeltaQ/(1+DeltaQ).
assert sp.expand(series_delta - (-DeltaQ + DeltaQ**2)) == 0, \
    'small-DeltaQ series does not match -DeltaQ + DeltaQ**2'
```

The variable `NQ` is already declared positive/real on line 5; `assert` statements use the existing symbol, no new declarations needed.

**Trivial-case pre-check (self-test from auditor):**
- `(mhat0**2*chiQ*NQ - 1).subs({mhat0:1, NQ:1/chiQ})` → `chiQ*(1/chiQ) - 1 = 0`. ✓
- `(mhat0**2*chiQ*NQ - 1).subs({mhat0:1, chiQ:1, NQ:1})` → `1*1 - 1 = 0`. ✓
- `(mhat0**2*chiQ*NQ - 1).subs({mhat0:1, chiQ:1+DeltaQ, NQ:1/(1+DeltaQ)})` → `(1+DeltaQ)/(1+DeltaQ) - 1 = 0`. ✓
- `series_delta - (-DeltaQ + DeltaQ**2)` where `series_delta = DeltaQ**2 - DeltaQ` → `0`. ✓
- Falsifiability: mutating `1/chiQ` to `1/chiQ**2` makes the first residual `chiQ/chiQ**2 - 1 = 1/chiQ - 1 ≠ 0`. ✓

**Verification command:**
After Codex applies, verifier runs `redteam exec-sympy 101` and confirms the script exits 0 with the four new `assert` statements visible in the file (grep `assert` should return at least 4 lines).

---

## F3 — insufficient_verification (small-`Delta_Q` linearization not anchored)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py:14-16`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:44-46`

**Issue:** The series `N_Q - 1 = -Delta_Q + O(Delta_Q^2)` (notes line 77) is computed and printed but never compared to its paper-stated form. A sign flip or order mismatch would not fail the audit.

**Required change:**

The SymPy side is handled by the `series_delta` assertion already inserted under F2. Do not add a duplicate.

On the Mathematica side, insert a new `expectZero` immediately after line 46 (the existing `Print["small-DeltaQ series = ...]`) and before line 47:

```
expectZero["small-DeltaQ series matches paper", Expand[seriesDelta - (-deltaQ + deltaQ^2)]];
```

**Trivial-case pre-check (self-test from auditor):**
The transcript shows `seriesDelta = -deltaQ + deltaQ^2`, so the residual is `0`. ✓ Mutating to `deltaQ + deltaQ^2` makes the residual `2*deltaQ ≠ 0`. ✓

**Verification command:**
After Codex applies, verifier confirms a fourth `PASS` line appears in the Mathematica transcript: `PASS: small-DeltaQ series matches paper`.

---

## F4 — paper_misalignment (subtype: script_missing_paper_claim)

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_101.tex:21-25` — the `\stagefield{Checks}` block lists three checks:
  > "Check the product `\widehat m_0^{\,2}\chi_QN_Q` keeps source, conservative, and outgoing factors separate. Check that higher odd terms begin beyond the point-particle 2.5PN coefficient. Check the outgoing `l=2` DtN fingerprint against the normalized `z=\omega a/c_s` expansion."
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:295` quote: "Higher odd denominator terms beginning at `O(\omega^7)` are invisible to the point-particle 2.5PN coefficient." (Check 2 substance)
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:305-326` — derivation of `\Lambda_2^{\rm out}(z) = -3 + z^2/3 + z^4/9 + i z^5/9 + O(z^6)` and `\widehat Y_2^{\rm out}(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)`, with `\chi_Q = 1` fixed by comparison (Check 3 substance).

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py:1-18` — no `omega` symbol, no `z`, no Hankel-function expansion. Neither Check 2 nor Check 3 is exercised.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl:1-52` — same.

## Resolve before fix_loop

The stage card promises three `Checks`, but the scripts only honor the first (implicitly, via symbol separation). Two questions for the user:

1. **Does Stage 101 own Check 2** (higher odd terms begin beyond point-particle 2.5PN), or is it carried forward from an upstream stage? The notes for Stage 101 do not derive it; the appendix asserts it in a single sentence (line 295) without scripted backing in this unit.
2. **Does Stage 101 own Check 3** (outgoing `l=2` DtN fingerprint vs `z = omega a/c_s` expansion), or is it Stage 80's deliverable? The notes file (line 41-51) explicitly attributes the canonical compact branch's `chi_Q = 1` identification to Stage 80, suggesting Check 3 belongs there, not here.

Possible directions (the user picks one for each check):

- **(a) Stage 101 owns Check 2 and Check 3.** Then the scripts must be extended:
  - Check 2: introduce `omega`, expand the denominator `1 - omega^2/Omega_Q^2 - i*chi_Q*sigma_Q_can*omega^5` to `O(omega^7)`, and assert the first non-`(omega^2, omega^5)` term is at `omega^7` or higher.
  - Check 3: introduce `z = omega*a/c_s`, compute `Lambda_2_out(z) = z * d/dz log(SphericalHankelH1(2, z))` (or the closed-form `-3 + z^2/3 + z^4/9 + i z^5/9 + O(z^6)`), assert series equals the appendix-stated form, then assert `chi_Q = 1` falls out of matching `\widehat Y_2^{\rm out}` to the `\widehat Y_Q^{\rm ret}` one-pole form.
  Then re-run sympy+mathematica.

- **(b) Stage 101 does NOT own these checks; they belong to upstream stages.** Then the stage_101.tex `\stagefield{Checks}` block (lines 21-25) should be trimmed: Check 2 cited to whatever stage actually verifies it (likely the upstream `2.5PN audit` unit referenced in the notes), Check 3 cited to Stage 80. No script change in unit 101.

- **(c) Hybrid: Stage 101 owns Check 1 only, with Check 2 and Check 3 explicitly cross-referenced** — paper card adds "(see Stage NNN)" pointers after each promised check, no script change.

The orchestrator will not invoke Codex on F4 in unit 101 until the user has chosen a direction. F1, F2, F3 may proceed independently.

---
## Applied: F4 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster B direction (c))

- files_changed: scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py
- summary: SymPy docstring annotates Check 2 (higher odd terms beyond 2.5PN) and Check 3 (DtN l=2 fingerprint pinning chi_Q = 1) as upstream carry-ins from stage 102 and stage 097 respectively. Stage 101 owns Check 1 (factorization mhat_0^2 chi_Q N_Q = 1) which is verified via F1/F2 input-equation anchors.
- deviation: none

## Applied: F1
- files_changed: mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl
- summary: Replaced three tautological expectZero calls that substituted nQExact into its own definition with non-tautological anchors against the INPUT factorization mhat0^2 * chiQ * nQ - 1 substituted with candidate NQ values. Plus banner sweep.
- deviation: none

## Applied: F2
- files_changed: scripts/moving_throat_pde_stage101_natural_source_map_reduction_sympy_audit.py
- summary: Added 4 substantive asserts anchored to the input factorization mhat0^2*chiQ*NQ - 1 with candidate NQ values, plus a small-DeltaQ linearization check.
- deviation: none

## Applied: F3
- files_changed: mathematica/moving_throat_pde_stage101_natural_source_map_reduction_mathematica_audit.wl
- summary: Added expectZero "small-DeltaQ series matches paper" asserting seriesDelta - (-deltaQ + deltaQ^2) = 0.
- deviation: none
