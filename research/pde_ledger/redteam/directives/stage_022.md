---
unit_id: 022
batch: I.2
created_at: 2026-05-21T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-21T19:45:45Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 022

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl:33-149`

**Issue:**

The Mathematica audit script is a section-by-section transliteration of the SymPy audit script. Each of Sections I, II, IV, V uses the identical "build factored form → take `Series` in omega around 0 → read `Coefficient` of each power" recipe, and Section V uses the same explicit `j2`, `y2`, `h2 = j2 + i y2` polynomial-rational form as SymPy rather than Mathematica's built-in `SphericalHankelH1[2, z]`. This violates the two-engine independence policy: both engines walk the same algebraic path, so a bug in the shared idea cannot be caught.

**Required change:**

Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl` as follows. Keep every `expectZero` target literal on the RHS unchanged — these are the claim. Change only the *route* by which the LHS of each assertion is computed, so it is not "build the same intermediate and read the same Series coefficient" as the SymPy script.

1. **Section I (currently lines 33-41) — inverse-relation route for `u2, u4`.**

   Replace lines 33-41 with:

   ```
   dCons = d0 + d2*omega^2 + d4*omega^4;
   (* Inverse-relation route: assume Y_resp = 1 + u2Sym omega^2 + u4Sym omega^4 + O(omega^6)
      and impose Y_resp * dCons - d0 = 0 modulo omega^6. *)
   Clear[u2Sym, u4Sym];
   yRespCand = 1 + u2Sym*omega^2 + u4Sym*omega^4;
   prod = Expand[yRespCand*dCons - d0];
   coeffEqs = Table[Coefficient[prod, omega, k] == 0, {k, 0, 4}];
   sol = First[Solve[coeffEqs, {u2Sym, u4Sym}]];
   u2 = FullSimplify[u2Sym /. sol, Assumptions -> $Assumptions];
   u4 = FullSimplify[u4Sym /. sol, Assumptions -> $Assumptions];

   banner["SECTION I — CONSERVATIVE OPERATOR TO RESPONSE"];
   Print["u2 (inverse-relation route) = ", fmt[u2]];
   Print["u4 (inverse-relation route) = ", fmt[u4]];
   expectZero["u2 formula", u2 + d2/d0];
   expectZero["u4 formula", u4 - (d2^2 - d0*d4)/d0^2];
   ```

   (Removes the `Series[d0/dCons, ...]` route entirely. The new route solves coefficient equations on the *product* `Y_resp * dCons - d0`, which is algebraically distinct from "Taylor-expand the quotient".)

2. **Section II (currently lines 43-66) — inverse-relation route for `P0, P2, P4`; direct polynomial expansion for `K0..K4, Gamma5`.**

   Replace lines 43-54 with:

   ```
   nFac = n0 + n2*omega^2 + n4*omega^4;
   (* Inverse-relation route: assume pref = p0Sym + p2Sym omega^2 + p4Sym omega^4 + O(omega^6)
      and impose pref * dCons^2 - d0 * nFac = 0 modulo omega^6. *)
   Clear[p0Sym, p2Sym, p4Sym];
   prefCand = p0Sym + p2Sym*omega^2 + p4Sym*omega^4;
   prodPref = Expand[prefCand*dCons^2 - d0*nFac];
   coeffEqsPref = Table[Coefficient[prodPref, omega, k] == 0, {k, 0, 4}];
   solPref = First[Solve[coeffEqsPref, {p0Sym, p2Sym, p4Sym}]];
   p0 = FullSimplify[p0Sym /. solPref, Assumptions -> $Assumptions];
   p2 = FullSimplify[p2Sym /. solPref, Assumptions -> $Assumptions];
   p4 = FullSimplify[p4Sym /. solPref, Assumptions -> $Assumptions];

   (* Branch coefficients: direct polynomial expansion (Expand only, no Series). *)
   y2Out = 1 + aa*omega^2 + bb*omega^4 + I*g5*omega^5;
   prefTrunc = p0 + p2*omega^2 + p4*omega^4;
   yBranch = Expand[prefTrunc*y2Out];
   k0 = FullSimplify[Coefficient[yBranch, omega, 0], Assumptions -> $Assumptions];
   k2 = FullSimplify[Coefficient[yBranch, omega, 2], Assumptions -> $Assumptions];
   k4 = FullSimplify[Coefficient[yBranch, omega, 4], Assumptions -> $Assumptions];
   gamma5 = FullSimplify[Coefficient[yBranch, omega, 5]/I, Assumptions -> $Assumptions];
   ```

   Leave lines 56-66 (the `banner` and the seven `expectZero` calls for `P0..Gamma5`) unchanged.

3. **Section III (currently lines 68-80) — leave as is.**

   No edit. The grouped real P2 inverse map is a 3-coefficient algebraic identity; both engines do the same algebra because there is only one route. Add a comment after line 73 (`banner["SECTION III — GROUPED REAL P2 INVERSE MAP"];`):

   ```
   (* Intentional parallel check: the 3x3 inverse-map identity admits no
      engine-independent route. Both engines verify the same algebra as a
      sanity cross-check. *)
   ```

4. **Section IV (currently lines 82-99) — inverse-relation route for `N0, N2, N4`; drop the IV.2 round-trip (handled by F2).**

   Replace lines 87-91 with:

   ```
   dProto = delta0 - s2*omega^2 + omega^4;
   pProto = p0proto - gW*omega^2;
   (* Inverse-relation route: assume nSeries = n0Cand + n2Cand omega^2 + n4Cand omega^4
      and impose nSeries * dProto^2 - pProto^2 = 0 modulo omega^6. *)
   Clear[n0Cand, n2Cand, n4Cand];
   nCand = n0Cand + n2Cand*omega^2 + n4Cand*omega^4;
   prodN = Expand[nCand*dProto^2 - pProto^2];
   coeffEqsN = Table[Coefficient[prodN, omega, k] == 0, {k, 0, 4}];
   solN = First[Solve[coeffEqsN, {n0Cand, n2Cand, n4Cand}]];
   n0Proto = FullSimplify[n0Cand /. solN, Assumptions -> $Assumptions];
   n2Proto = FullSimplify[n2Cand /. solN, Assumptions -> $Assumptions];
   n4Proto = FullSimplify[n4Cand /. solN, Assumptions -> $Assumptions];
   ```

   Leave lines 93-99 (the `banner` and the three `expectZero` calls for `N0, N2, N4`) unchanged.

5. **Section V (currently lines 111-136) — use `SphericalHankelH1[2, z]` instead of explicit `j2, y2, h2 = j2 + i y2`.**

   Replace lines 121-123 with:

   ```
   (* Use Mathematica's built-in spherical Hankel function instead of the
      explicit polynomial-rational form, so the derivation of A, B, G5 is
      independent of the SymPy script's choice of j2, y2 expressions. *)
   h2 = SphericalHankelH1[2, z];
   ```

   (Delete the explicit `j2 = ...`, `y2 = ...`, and `h2 = FullSimplify[j2 + I y2, ...]` lines; they are replaced by the single `SphericalHankelH1` definition.)

   Leave lines 124-136 unchanged. The `lambda2 = (omega D[h2, z]/h2) /. z -> omega a/cS` derivation, the series in omega, the `y2Resp`/`y2Static`/`y2Hat` construction, and the three `expectZero` assertions on `aStage4, bStage4, g5Stage4` against the target literals all stay as written. (Mathematica's `SphericalHankelH1[2, z]` is mathematically identical to `j_2(z) + i y_2(z)`, so the assertions still pass — but the LHS is now reached via Mathematica's built-in, not via a hand-typed polynomial that mirrors SymPy.)

6. **Section VI (currently lines 138-149) — leave as is.**

   No edit. This section just plugs the already-derived `g5Stage4 = a^5/(27 c_s^5)` and `aStage4, bStage4` into `Solve[mhat^2*p0Static*g5Stage4 == gammaGR]` and asserts the constant arithmetic. The route is identical to SymPy's, but the inputs (`g5Stage4` etc.) now arrive via the independent route established in step 5, so the chain is independent end-to-end.

**Self-test trace (the LHS variables and how they reach the assertions):**

- `u2 = u2Sym /. Solve[...]` returns `u2Sym -> -d2/d0`, so `u2 + d2/d0 → 0`. ✓
- `u4 = u4Sym /. Solve[...]` returns `u4Sym -> (d2^2 - d0 d4)/d0^2`, so `u4 - (d2^2 - d0 d4)/d0^2 → 0`. ✓
- `p0Sym -> n0/d0`, `p2Sym -> (d0 n2 - 2 d2 n0)/d0^2`, `p4Sym -> (d0^2 n4 - 2 d0(d2 n2 + d4 n0) + 3 d2^2 n0)/d0^3`. All three `expectZero` calls reduce to `0`. ✓
- `n0Cand -> p0proto^2/delta0^2`, `n2Cand -> 2 p0proto(p0proto s2 - delta0 gW)/delta0^3`, `n4Cand -> (delta0^2 gW^2 - 2 delta0 p0proto^2 - 4 delta0 p0proto s2 gW + 3 p0proto^2 s2^2)/delta0^4`. ✓
- `SphericalHankelH1[2, z]` returns `-((3 - 3 I z - z^2) E^(I z))/z^3` in closed form. `omega D[h, z]/h /. z -> omega a/cS` followed by Series gives the same `Lambda2` as the explicit form (since `j_2 + i y_2 = -((3 - 3 i z - z^2) e^(i z))/z^3` is a standard identity), so the subsequent `aStage4, bStage4, g5Stage4` extraction matches and the three target assertions still pass.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 022` and confirm:

(a) the script exits 0 and every `expectZero` check passes;

(b) Sections I, II, IV use `Solve[coeffEqs, ...]` on coefficient-equation lists (the inverse-relation route);

(c) Section V uses `SphericalHankelH1[2, z]` and no longer defines `j2 = ...`, `y2 = ...` explicitly;

(d) all original `expectZero` target literals (`u2 + d2/d0`, `u4 - (d2^2 - d0 d4)/d0^2`, `p0 - n0/d0`, etc., through `mhat=1 K4 target`) remain on the RHS unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- summary: Replaced Mathematica quotient-series routes with inverse-relation coefficient solves and switched the Stage-4 Hankel derivation to `SphericalHankelH1[2, z]`.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py:245-261`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl:104-109`

**Issue:**

The IV.2 "round-trip into Stage-4 symbols" check tests that `Series[Nproto.subs(dict_back)].coeff(omega, k) == Series[Nproto].coeff(omega, k).subs(dict_back)`. Since `dict_back` (which sends `Delta0, S2, P0_proto` to `Omega_A^2 Omega_W^2 - R^2`, `Omega_A^2 + Omega_W^2`, `Omega_A^2 gW + R g_A`) introduces no `omega`-dependence, Taylor expansion in `omega` commutes with the substitution. The residual is identically zero by construction — the assertion cannot fail and contributes no verification value.

**Required change:**

1. **SymPy.** In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`, delete lines 245-261 (the entire block beginning with `Nproto_stage4 = (` through the closing parenthesis of the third `expect_zero` call). The lines to delete in full are:

   ```
   Nproto_stage4 = (
       (dict_back[P0_proto] - gW * omega**2) ** 2
       / (dict_back[Delta0] - dict_back[S2] * omega**2 + omega**4) ** 2
   )
   Nseries_stage4 = sp.expand(sp.series(Nproto_stage4, omega, 0, 6).removeO())
   expect_zero(
       "N0 round-trip into Stage-4 symbols",
       Nseries_stage4.coeff(omega, 0) - sp.simplify(N0.subs(dict_back)),
   )
   expect_zero(
       "N2 round-trip into Stage-4 symbols",
       Nseries_stage4.coeff(omega, 2) - sp.simplify(N2.subs(dict_back)),
   )
   expect_zero(
       "N4 round-trip into Stage-4 symbols",
       Nseries_stage4.coeff(omega, 4) - sp.simplify(N4.subs(dict_back)),
   )
   ```

   Leave the `dict_back = {...}` construction and the three `print(...)` / `sp.pprint(...)` calls at lines 232-243 in place — those are documentation, not assertions.

2. **Mathematica.** In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`, delete lines 104-109. The lines to delete in full are:

   ```
   nStage4 = Expand[
     Normal[Series[((p0Back - gW*omega^2)^2)/(deltaBack - s2Back*omega^2 + omega^4)^2, {omega, 0, 4}]]
   ];
   expectZero["N0 round-trip", Coefficient[nStage4, omega, 0] - (n0Proto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back})];
   expectZero["N2 round-trip", Coefficient[nStage4, omega, 2] - (n2Proto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back})];
   expectZero["N4 round-trip", Coefficient[nStage4, omega, 4] - (n4Proto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back})];
   ```

   Leave the dictionary definitions (`deltaBack = ...`, `s2Back = ...`, `p0Back = ...`) at lines 101-103 in place; they are referenced by no other code after deletion, but removing them is a stylistic refactor and out of scope. (If linting complains about unused variables, that is a separate concern.)

**Self-test trace:**

- After deletion, `Nproto_stage4`, `Nseries_stage4`, and `nStage4` are no longer referenced anywhere else in either script. Verified by scanning the remaining sections (Section V in both scripts uses only `A_stage4, B_stage4, G5_stage4, NQ_target, K2_target, K4_target`, none of which depend on the deleted block).
- The IV.1 prototype-formula assertions (`N0 prototype`, `N2 prototype`, `N4 prototype`) at SymPy lines 217-222 and Mathematica lines 94-99 remain and continue to assert the non-tautological closed-form `N_k` expressions.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 022` and `redteam exec-mathematica 022`. Both must exit 0. The output transcripts must no longer contain the lines `N0 round-trip into Stage-4 symbols`, `N2 round-trip into Stage-4 symbols`, `N4 round-trip into Stage-4 symbols` (SymPy) or `N0 round-trip`, `N2 round-trip`, `N4 round-trip` (Mathematica). The IV.1 prototype assertions (`N0 prototype`, `N2 prototype`, `N4 prototype`) must still appear and still print residual `0`.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- summary: Removed the tautological IV.2 round-trip assertions from both audit scripts while preserving the dictionary printouts.
- deviation: none
