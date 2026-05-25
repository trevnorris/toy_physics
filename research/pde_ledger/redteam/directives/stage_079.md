---
unit_id: 079
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-23T05:34:45Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 079

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:67-68`

**Issue:** The two `expectApprox` calls at lines 67-68 compare `zeta0` and `zetaInf` against decimal constants (`1.00005192880219492404747131934` and `2.46752922945601123498982913352`) that were copied byte-for-byte from the SymPy output file. This is not an independent verification — Mathematica is round-tripping SymPy's own numbers back through itself. The fix in F1 and F2 is the same edit; perform it once and mark both findings applied.

**Required change:**
Replace lines 67-68:

Before:
```
expectApprox["zeta_F1(0+) numeric check", zeta0, ToExpression["1.00005192880219492404747131934`30"], 10^-14];
expectApprox["zeta_max^(F1) numeric check", zetaInf, ToExpression["2.46752922945601123498982913352`30"], 10^-14];
```

After:
```
zeta0Sym = FullSimplify[Limit[aF1*omega^2, pe -> 0, Direction -> "FromAbove"], Assumptions -> $Assumptions];
zetaInfSym = FullSimplify[Limit[aF1*omega^2, pe -> Infinity], Assumptions -> $Assumptions];
expectZero["zeta_F1(0+) - A_F1", zeta0Sym - aF1];
expectZero["zeta_F1(inf) - A_F1 Pi^2/4", zetaInfSym - aF1*Pi^2/4];
```

These new checks take the symbolic limit of the *product* `aF1*omega^2` (not of `omega` alone followed by multiplication), exercising the same derivation path SymPy uses in its A3/A4 checks on script lines 62-63.

If `FullSimplify` of the residual does not symbolically reduce to exact `0` (because `aF1` is a high-precision numeric atom rather than a symbol), fall back to using `expectApprox` with target `0` and tolerance `10^-40`:
```
expectApprox["zeta_F1(0+) - A_F1", N[zeta0Sym - aF1, 50], 0, 10^-40];
expectApprox["zeta_F1(inf) - A_F1 Pi^2/4", N[zetaInfSym - aF1*Pi^2/4, 50], 0, 10^-40];
```
Choose `expectZero` first; only fall back to `expectApprox 1e-40` if a manual symbolic reduction is not guaranteed for a numeric `aF1`. Either way, no decimal target near `1.00005...` or `2.46752...` may remain in the file.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 079` and confirm:
1. The strings `1.00005192880219492404747131934` and `2.46752922945601123498982913352` no longer appear in the `.wl`.
2. The script's output contains `PASS: zeta_F1(0+) - A_F1` and `PASS: zeta_F1(inf) - A_F1 Pi^2/4`.
3. Exit code 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`
- summary: Replaced the hardcoded zeta numeric assertions with symbolic product-limit residual checks using the permitted approximate-zero fallback.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:67-68`

**Issue:** Because `omega0` (line 52) and `omegaInf` (line 53) are already evaluated to exact symbolic values (`1` and `Pi/2`) before lines 61-62 form `zeta0` and `zetaInf`, the two `expectApprox` checks on lines 67-68 reduce to `aF1 == hardcoded` and `aF1*Pi^2/4 == hardcoded` — they cannot detect an error in how Mathematica computes `Limit[aF1*omega^2, ...]`, because that limit is never taken on the *product*; the substitution happens before the assertion. The substitution also makes the comparison value tautological in `omega0` / `omegaInf`.

**Required change:**
Apply the same edit as F1 (the two findings share lines 67-68). The replacement code in F1 introduces `zeta0Sym = Limit[aF1*omega^2, pe -> 0, ...]` and `zetaInfSym = Limit[aF1*omega^2, pe -> Infinity]` so that Mathematica's `Limit` is exercised on the product, not on pre-substituted scalars. Lines 61-62 (`zeta0 = N[aF1*omega0^2, 50]; zetaInf = N[aF1*omegaInf^2, 50];`) and the `Print` calls on lines 65-66 may remain untouched (they are diagnostic prints, not assertions). Do not delete them.

After F1's edit lands, this finding is resolved as a byproduct — mark `## Applied: F2` referencing the same `files_changed` and note `"resolved by F1 edit"` in `summary`.

**Verification command:**
Same as F1.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`
- summary: resolved by F1 edit.
- deviation: none

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:74` (insert new check immediately after line 74; do not modify lines 32-74's existing checks)

**Issue:** The Mathematica script is a syntactic transliteration of the SymPy script: same closed-form for `Omega(Pe)`, same `A_F1` formula, same series-comparison pattern. Both engines compute the same limits and series of the same closed expression. If the closed form for `Omega(Pe)` is wrong, both engines pass identically. The script needs one substantive independent check that exercises a different code path than `sp.series`. Symbolic differentiation `D[omega, pe]` followed by `Limit[..., pe -> 0]` lands on the slope `(4-Pi)/(2*Pi)` via a different route than SymPy's `series` method, and forces Mathematica's symbolic-differentiation engine to corroborate the closed form.

**Required change:**
Insert the following lines immediately after the existing `small-Pe linear coefficient check` at line 74, and before the blank `Print[""]` at line 76:

```
omegaPrime0 = FullSimplify[Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"], Assumptions -> $Assumptions];
Print["Omega'(0+) = ", fmt[omegaPrime0]];
expectApprox["Omega'(0+) - (4-Pi)/(2 Pi)", N[omegaPrime0 - (4 - Pi)/(2*Pi), 50], 0, 10^-40];
```

Rationale for `expectApprox` (rather than `expectZero`): `Limit[D[omega, pe], pe -> 0]` may not symbolically reduce to exact `(4-Pi)/(2*Pi)` after `FullSimplify` because the derivative of `omega` is a non-trivial rational expression in `Exp[pe]` and `pe`; numeric comparison at 50-digit precision against tolerance `10^-40` is the robust form. The check is substantive because:
- `omega` truly depends on `pe` (so `D[omega, pe]` is non-trivial, not identically zero — confirmed by the small-Pe series check at line 74 producing a nonzero linear coefficient).
- The target `(4 - Pi)/(2 Pi)` is independently derivable from the closed form by Taylor expansion: `numerator ~ pi^2*Pe + 2*pi*Pe^2 + ...`, `denominator ~ pi^2*Pe*(1 + Pe/2) + ...`, ratio `~ 1 + (2/Pi - 1/2)*Pe + O(Pe^2) = 1 + ((4-Pi)/(2*Pi))*Pe + O(Pe^2)`. So `Omega'(0+) = (4-Pi)/(2*Pi) ≈ 0.13662`.
- It uses Mathematica's `D` (symbolic differentiation) rather than `Series`, so a bug in either's small-Pe expansion would be caught by disagreement with `(4-Pi)/(2*Pi)`.

Do not change any other lines. Do not modify line 76 (blank `Print[""]`), line 77 (passed banner), or line 79 (`Exit[0]`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 079` and confirm:
1. The script contains a `D[omega, pe]` expression (it did not before).
2. The output contains `PASS: Omega'(0+) - (4-Pi)/(2 Pi)`.
3. Exit code 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`
- summary: Added the independent symbolic-differentiation slope check for `Omega'(0+)`.
- deviation: none
