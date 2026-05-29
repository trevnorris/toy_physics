---
unit_id: 117
batch: IV.3
created_at: 2026-05-29
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-29T13:20:57-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 117

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## RESOLVED — Claude+Codex consult 019e74f7 (2026-05-29)

All four findings were resolved by a Claude+Codex consult; both parties CONCUR'd on every direction. Summary of the agreed dispositions:

- **F1** = accepted policy-mirror. The `.wl` is a sanctioned near-line-by-line cross-engine port of a pure rational series-coefficient classification; a Padé/residue reformulation adds little adversarial value. NO script change — record only. (Per the Mathematica mirror policy.)
- **F2** = upstream-import de-tautologization. The hardcoded `kappa0 = (1+r_c)/3` is replaced by `kappa0` built from the stage-116 **forward** tube-length closed form `L_W = pi*a*sqrt((1+r_c)/3)/2` via `kappa0 = 4*L_W**2/(pi**2*a**2)`, so the `kappa_c = 1/3` / core-residual check now EXERCISES the tube-length law rather than restating its target.
- **F3** = comment-only provenance fix. `gamma0 = (1+r_c)/9` is NOT derived anywhere — it is a pure-scale **ansatz** (modeling choice) of the canonical `l=2` branch, postulated in the stage-116 note §"Bare outgoing normalization". The script's false `# (Stage 119)` / "derived upstream" attributions for `gamma0` are corrected; the `gamma_c = 1/9` assertion is unchanged (acknowledged consistency-of-assumption check). Logged to the planned ansatz catalog.
- **F4** = reconcile-already-applied. The capstone survivability flags are already wired to the computed section-1–5 residuals (no longer hardcoded `True`); confirm + record.

Neither F2 nor F3 is conceptual: `kappa0` is genuinely forward-derived upstream (stage 116, D/N half-wave eigenvalue → tube-length law), while `gamma0` is explicitly postulated as an ansatz in the stage-116 note. The fixes only (F2) reroute an already-true substitution through its genuine forward source and (F3) correct prose comments; no assertion target value changes.

## F1 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl` (whole file)

**Issue:**
The `.wl` is a near-line-by-line port of the `.py`. The payload is pure rational series-coefficient classification: which outlet-deformation classes preserve `Ŷ₂ = 1 + z²/9 + 4z⁴/81 + i z⁵/27` to O(z⁵). An independent Padé/residue reformulation of this rational-series bookkeeping adds little adversarial value over the SymPy run.

**Required change:**
NONE. Per the Claude+Codex consult 019e74f7 and the Mathematica mirror policy, this `.wl` is accepted as a sanctioned cross-engine mirror, not as an independent-route second engine. Do not rewrite it for F1. (You will still edit the `.wl` under F2 and F3 below for the mirrored substitution and comment fixes; that is unrelated to this F1 disposition.)

**Verification command:**
No verification needed for F1 beyond the F2/F3/F4 runs. Record the disposition only.

When recording, append:

```
## Applied: F1

- files_changed: none
- summary: accepted policy-mirror per Claude+Codex consult 019e74f7; no independent-route rewrite
- deviation: no change (mirror accepted)
```

## Applied: F1

- files_changed: none
- summary: accepted policy-mirror per Claude+Codex consult 019e74f7; no independent-route rewrite
- deviation: no change (mirror accepted)

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:164` and `:174` (and the comment block 165-168)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:118` and `:126` (and the comment block 119-122)

**Issue:**
The `kappa_c = 1/3` chain in Section 5 is currently tautological. Two facts establish this:

1. `Lw_required` (SymPy L164 / `lWRequired` `.wl` L118) is computed by INVERTING the very relation under test: `sp.solve(sp.Eq(4*Lw**2/(pi**2*a**2), (1+r_c)/3), Lw)[0]`. Feeding that back would be tautological — it solves for `L_W` from `kappa0 = (1+r_c)/3`, then would reconstruct `kappa0 = (1+r_c)/3` by definition. (`Lw_required` is also presently dead: it is never consumed by `delta_core`.)
2. The load-bearing `delta_core` substitution (SymPy L174 / `.wl` L126) hardcodes `kappa0 -> (1+r_c)/3`. With `kappa0` hardcoded to its target, the resulting `kappa_c = kappa0/(1+r_c) = 1/3` and the core-residual check cannot fail on a wrong tube-length-law coefficient.

The genuine upstream forward result is at stage 116: from the D/N half-wave eigenvalue `k_W = pi/(2 L_W)` one derives `kappa0 = 4 L_W**2/(pi**2 a**2)` (stage116 SymPy L43–47), and the compensation requirement then fixes the tube length to the closed form `L_W = pi*a*sqrt((1+r_c)/3)/2` (stage116 SymPy L60–62; stage-116 note §"Compensation-selected tube length"). The fix routes the stage-117 substitution through that forward chain.

**De-tautologization route chosen (anti-tautology / X−X pitfall):**
Because `Lw_required`/`lWRequired` is obtained by inverting the target, we do NOT feed it back. Instead we build `kappa0` from the stage-116 **forward** tube-length closed form written explicitly:

```
L_W_forward = pi*a*sqrt((1+r_c)/3)/2          # stage-116 boxed result (forward)
kappa0_from_tube = 4*L_W_forward**2/(pi**2*a**2)   # = (1+r_c)/3 IFF the law is correct
```

then substitute `kappa0 -> kappa0_from_tube` into `delta_core`. This is genuinely non-tautological: if the tube-length-law coefficient were wrong — e.g. `sqrt((1+r_c)/3)` were instead `sqrt((1+r_c)/2)`, or the prefactor `pi*a/2` were `pi*a`, or the exponent on `L_W` differed — then `kappa0_from_tube != (1+r_c)/3`, so `kappa_c != 1/3`, and the `delta_core - delta_core_expected` residual at O(z⁶) (SymPy L178 / `.wl` L128-131) would be NONZERO and the check would FAIL. The check now verifies that the stage-116 closed form actually yields `kappa_c = 1/3`.

**Required change (SymPy):**

Before (lines 164-168):
```python
Lw_required = sp.solve(sp.Eq(4 * Lw**2 / (sp.pi**2 * a**2), (1 + r_c) / 3), Lw)[0]
# Stage 116 fixes kappa_0_bare = (1+r_c)/3 via the D/N tube; carrying forward
# to kappa_c = kappa_0/(1+r_c) = 1/3 is then arithmetic, not an independent check.
print("carrying forward (Stage 116): kappa_0_bare = (1+r_c)/3 -> kappa_c = 1/3")
print("carrying forward (Stage 119): gamma_0_bare = (1+r_c)/9 -> gamma_c = 1/9")
```

After (the `Lw_required` inversion is replaced by the FORWARD tube-length closed form, and `kappa0_from_tube` is derived from it; F3 also rewrites the two print lines — see F3 for those two lines):
```python
# De-tautologized (F2): build kappa0 from the stage-116 FORWARD tube-length law,
# not by inverting kappa0 = (1+r_c)/3. The stage-116 boxed result is
#   L_W = pi a sqrt((1+r_c)/3) / 2   (forward, from the D/N half-wave eigenvalue),
# and kappa0 = 4 L_W^2 / (pi^2 a^2). The kappa_c = 1/3 / core-residual check below
# therefore EXERCISES that closed form: a wrong tube-length coefficient would make
# kappa0_from_tube != (1+r_c)/3, so delta_core would no longer collapse and the
# O(z^6) residual check would FAIL.
L_W_forward = sp.pi * a * sp.sqrt((1 + r_c) / 3) / 2
kappa0_from_tube = sp.simplify(4 * L_W_forward**2 / (sp.pi**2 * a**2))
print("carrying forward (Stage 116): L_W = pi a sqrt((1+r_c)/3)/2 -> kappa_0_bare = 4 L_W^2/(pi^2 a^2) -> kappa_c = 1/3")
print("gamma_0_bare = (1+r_c)/9 is a pure-scale ANSATZ of the canonical l=2 branch (stage-116 note), not derived; gamma_c = 1/9 is a consistency-of-assumption check")
```

Then, in the `delta_core` substitution dictionary (lines 170-176), replace the hardcoded `kappa0: (1 + r_c) / 3` with the tube-derived value. Keep `gamma0: (1 + r_c) / 9` unchanged (it is the ansatz; see F3).

Before (lines 170-176):
```python
delta_core = sp.simplify(
    rho_c - sigma_c / (1 - kappa_c * z**2 - I * gamma_c * z**5)
).subs({
    gq: gq_solutions[0],
    kappa0: (1 + r_c) / 3,
    gamma0: (1 + r_c) / 9,
})
```

After:
```python
delta_core = sp.simplify(
    rho_c - sigma_c / (1 - kappa_c * z**2 - I * gamma_c * z**5)
).subs({
    gq: gq_solutions[0],
    kappa0: kappa0_from_tube,
    gamma0: (1 + r_c) / 9,
})
```

Note: `kappa0_from_tube` simplifies to `(1+r_c)/3`, so the existing assertion LABELS and the `delta_core` collapse residual (L178) are unchanged — output strings still match. The difference is that the substituted value now flows from the explicit tube-length closed form.

Also update the section-5 comment block (lines 141-149) for the `kappa0` provenance phrasing — see F3 for the exact comment edit (F3 owns all comment-text changes; F2 owns only the code substitution above). Apply F2's code edits and F3's comment edits both; they touch adjacent lines but do not overlap.

**Required change (Mathematica):**

Mirror the same substitution in the `.wl`.

Before (lines 118-122):
```mathematica
lWRequired = FullSimplify[lW /. First[Solve[4 lW^2/(Pi^2 a^2) == (1 + rC)/3, lW, Reals]]];
(* Stage 116 fixes kappa_0_bare = (1+r_c)/3 via the D/N tube; carrying forward
   to kappa_c = kappa_0/(1+r_c) = 1/3 is then arithmetic, not an independent check. *)
Print["carrying forward (Stage 116): kappa_0_bare = (1+r_c)/3 -> kappa_c = 1/3"];
Print["carrying forward (Stage 119): gamma_0_bare = (1+r_c)/9 -> gamma_c = 1/9"];
```

After:
```mathematica
(* De-tautologized (F2): build kappa0 from the stage-116 FORWARD tube-length law,
   not by inverting kappa0 = (1+r_c)/3. Stage-116 boxed result:
   L_W = Pi a Sqrt[(1+r_c)/3]/2, and kappa0 = 4 L_W^2/(Pi^2 a^2). A wrong
   tube-length coefficient would make kappa0FromTube != (1+r_c)/3, so deltaCore
   would no longer collapse and the O(z^5) residual check would FAIL. *)
lWForward = Pi a Sqrt[(1 + rC)/3]/2;
kappa0FromTube = FullSimplify[4 lWForward^2/(Pi^2 a^2)];
Print["carrying forward (Stage 116): L_W = Pi a Sqrt[(1+r_c)/3]/2 -> kappa_0_bare = 4 L_W^2/(Pi^2 a^2) -> kappa_c = 1/3"];
Print["gamma_0_bare = (1+r_c)/9 is a pure-scale ANSATZ of the canonical l=2 branch (stage-116 note), not derived; gamma_c = 1/9 is a consistency-of-assumption check"];
```

Then, in the `deltaCore` replacement-rule list (line 126), replace the hardcoded `kappa0 -> (1 + rC)/3` with the tube-derived value, keeping `gamma0 -> (1 + rC)/9` (the ansatz):

Before (lines 124-126):
```mathematica
deltaCore = FullSimplify[
  rhoC - sigmaC/(1 - kappaC z^2 - I gammaC z^5)
] /. First[gqSolutions] /. {kappa0 -> (1 + rC)/3, gamma0 -> (1 + rC)/9};
```

After:
```mathematica
deltaCore = FullSimplify[
  rhoC - sigmaC/(1 - kappaC z^2 - I gammaC z^5)
] /. First[gqSolutions] /. {kappa0 -> kappa0FromTube, gamma0 -> (1 + rC)/9};
```

`kappa0FromTube` simplifies to `(1+rC)/3`, so the `deltaCore` collapse residual (line 128-131) and its `PASS` label are unchanged.

Also update the section-5 comment block (lines 98-103) per F3 (comment text only).

**Claim manifest:**
- The `kappa_c = 1/3` consequence is no longer asserted by hardcoding `kappa0 = (1+r_c)/3`; it now follows from the stage-116 forward tube-length closed form `L_W = pi a sqrt((1+r_c)/3)/2`. The load-bearing assertion (`concrete core collapses to the compensated hybrid class`, SymPy L178 / `.wl` L128-131) would FAIL if that closed form's coefficient were wrong.

**Verification command:**
After Codex applies, run `python3 scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py` and `math -script mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`. Both must exit 0. The existing check line `concrete core collapses to the compensated hybrid class = 0` (SymPy) and `PASS: concrete core collapses to the compensated hybrid class` (Mathematica) must still appear; no assertion label changes.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`
- summary: replaced the inverted/hardcoded kappa0 substitution with kappa0 derived from the stage-116 forward tube-length closed form in both audit scripts
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py` — comment block at lines 141-149 (the "derived upstream" phrasing that lumps `kappa_0` and `gamma_0`), and the second `print(...)` line (originally line 168, now the F2-rewritten `gamma_0` print). These are COMMENTS / printed-status strings only.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl` — comment block at lines 98-103, and the second `Print[...]` (originally line 122, now the F2-rewritten `gamma_0` print).

**Issue:**
`gamma0 = (1+r_c)/9` is NOT derived anywhere. It is a pure-scale **ansatz** — a modeling choice — postulated in `notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md:57-75` (§"Bare outgoing normalization": "A simple concrete realization is to take the bare mixed outlet to be a pure-scale deformation of the canonical compact outgoing `l=2` branch"). It is carried as a hardcoded input at stage 115 (`...stage115_core_balance_compensation_sympy_audit.py:47`, `gamma0_can = (1+r_c)/9`) and stage 116 (`...stage116_..._sympy_audit.py:73`, whose comment correctly says it is an "upstream-carried input ... not derived in-stage ... an expect_zero here would be tautological").

The stage-117 script currently mis-attributes `gamma0`'s provenance:
1. The section-5 comment (SymPy L141-149 / `.wl` L98-103) says the forward expressions for BOTH `kappa_0_bare = (1+r_c)/3` AND `gamma_0_bare = (1+r_c)/9` are "derived upstream" — true for `kappa_0`, FALSE for `gamma_0`.
2. The status print (original SymPy L168 / `.wl` L122) read `carrying forward (Stage 119): gamma_0_bare = (1+r_c)/9 ...`. Stage 119 derives only `kappa_0` (its parent-balance/tube-length algebra); it does NOT derive `gamma_0`. (Note: F2 above already rewrites this exact print line to the corrected ansatz wording — so for the print line, F2's "After" text IS the F3 fix; do not double-edit it. F3 here additionally fixes the comment BLOCK.)

The `gamma_c = 1/9` assertion itself is NOT changed: it is an acknowledged consistency-of-assumption check and cannot be de-tautologized, because `gamma0` is postulated rather than forward-derived. (Logged to the planned ansatz catalog.)

This is a COMMENT-ONLY change (red-team may edit script comments). The `kappa_0` "derived upstream" attribution is CORRECT and stays, pointing to stage 116 (D/N eigenvalue → tube-length law).

**Required change (SymPy) — comment block lines 141-149:**

Before:
```python
# This stage is a status-consolidation card. The forward expressions for
# kappa_0_bare = (1+r_c)/3 and gamma_0_bare = (1+r_c)/9 are derived upstream:
#   - Stage 115: parent-overlap reparametrization and Schur reduction;
#   - Stage 116: D/N half-wave eigenvalue k_W = pi/(2 L_W) on q(0)=0, q'(L_W)=0,
#                yielding kappa_0_bare from the tube-length law.
# Here we substitute those forward expressions and read off the *arithmetic
# consequences* kappa_c = kappa_0/(1+r_c) = 1/3 and gamma_c = gamma_0/(1+r_c) = 1/9.
# Those substitutions are not independent derivations; the load-bearing check in
# this stage is the residual delta_core - delta_core_expected at order z^6.
```

After:
```python
# This stage is a status-consolidation card. Provenance of the two bare
# mixed-channel coefficients DIFFERS and must not be lumped together:
#   - kappa_0_bare = (1+r_c)/3 is FORWARD-DERIVED upstream at Stage 116 from the
#     D/N half-wave eigenvalue k_W = pi/(2 L_W) on q(0)=0, q'(L_W)=0, giving
#     kappa_0 = 4 L_W^2/(pi^2 a^2) and the tube-length law L_W = pi a sqrt((1+r_c)/3)/2.
#     F2 routes the substitution below through that forward closed form, so the
#     kappa_c = 1/3 / core-residual check exercises (and could falsify) the law.
#   - gamma_0_bare = (1+r_c)/9 is NOT derived: it is a pure-scale ANSATZ (modeling
#     choice) of the canonical compact outgoing l=2 branch, postulated in the
#     stage-116 note (sec. "Bare outgoing normalization") and carried as a hardcoded
#     input at Stages 115/116. The gamma_c = gamma_0/(1+r_c) = 1/9 line is therefore
#     an acknowledged consistency-of-assumption check, not an independent derivation,
#     and cannot be de-tautologized while gamma_0 stays postulated.
# The load-bearing check in this stage is the residual delta_core - delta_core_expected
# at order z^6.
```

(The `kappa_c`/`gamma_c` and `r_c` etc. assignment lines that follow are unchanged. The `print` lines previously at 165-168 are handled by F2's rewrite above — the corrected `gamma_0` ansatz print is part of F2's "After" block.)

**Required change (Mathematica) — comment block lines 98-103:**

Before:
```mathematica
(* Status-consolidation card. Forward expressions for kappa_0_bare = (1+r_c)/3
   and gamma_0_bare = (1+r_c)/9 come from stages 115 (Schur reduction) and 116
   (D/N half-wave eigenvalue: k_W = Pi/(2 L_W) on q(0)=0, q'(L_W)=0). Here we
   substitute those forward expressions and read off kappa_c = 1/3, gamma_c = 1/9
   as arithmetic consequences. The load-bearing check in this stage is the
   residual deltaCore - deltaCoreExpected at order z^5. *)
```

After:
```mathematica
(* Status-consolidation card. Provenance of the two bare coefficients DIFFERS:
   - kappa_0_bare = (1+r_c)/3 is FORWARD-DERIVED at Stage 116 from the D/N
     half-wave eigenvalue k_W = Pi/(2 L_W) on q(0)=0, q'(L_W)=0, giving
     kappa_0 = 4 L_W^2/(Pi^2 a^2) and the tube-length law
     L_W = Pi a Sqrt[(1+r_c)/3]/2. F2 routes the substitution below through that
     forward closed form, so the kappa_c = 1/3 / core-residual check exercises it.
   - gamma_0_bare = (1+r_c)/9 is NOT derived: it is a pure-scale ANSATZ of the
     canonical compact outgoing l=2 branch, postulated in the stage-116 note
     ("Bare outgoing normalization") and carried as a hardcoded input at
     Stages 115/116. gamma_c = 1/9 is thus a consistency-of-assumption check.
   The load-bearing check is the residual deltaCore - deltaCoreExpected at z^5. *)
```

(The `gamma_0` status `Print` line is corrected as part of F2's "After" block above.)

**Claim manifest:**
- `gamma_0 = (1+r_c)/9` is documented as a postulated ansatz (stage-116 note §"Bare outgoing normalization"), not a derived quantity; the `gamma_c = 1/9` assertion is an explicit consistency-of-assumption check. No assertion logic changes.

**Verification command:**
After Codex applies, run both scripts. Both must exit 0; all existing `PASS`/`= 0` lines unchanged. Confirm the printed status string no longer references "Stage 119" for `gamma_0` and instead labels it an ansatz.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`
- summary: corrected the section-5 provenance comments and gamma0 status print to document gamma0 as a pure-scale ansatz rather than a derived quantity
- deviation: none

## F4 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:184-204` (the capstone flag computations)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:136-156`

**Issue / current state:**
ALREADY SUBSTANTIALLY APPLIED by a prior pass. The capstone classification survivability flags are NO LONGER hardcoded `True`. Current state (confirmed by reading the live scripts):

SymPy L184-204:
- `even_ok_*` / `odd_ok_*` are wired to the section-1–4 residuals (e.g. `even_ok_scale = ({sol[beta] for sol in beta_solutions} == {-1, 1})`, `odd_ok_scale = (sp.simplify(chi_arg.subs(beta, 1) - 1) == 0)`, similarly for robin/standalone/hyb_cancel/compensated).
- `nontrivial_scale`, `nontrivial_robin`, `nontrivial_hyb_cancel` are literal `False` — but these are CORRECT physics (those classes are trivial; comments explain "scale class is trivially canonical with rescale", "rho_R = 0 only", "collapses at gamma = 0"), not hardcoded survivability passes.
- `nontrivial_standalone = (sp.simplify(sigma_match) != 0)` — computed.
- `nontrivial_compensated` (L202-204) is wired to the section-5 series residual `delta_core - delta_core_expected == 0` — the load-bearing entry.

Mathematica L136-156 mirrors this exactly (`evenOk*`/`oddOk*` from residuals; `nontrivialScale/Robin/HybCancel = False`; `nontrivialStandalone = !TrueQ[... sigmaMatch === 0]`; `nontrivialCompensated` anchored to the L154-156 series residual).

So the original F4 intent (wire flags to computed booleans, not hardcoded `True`) is already met.

**Required change:**
NONE to the flag logic — it is already wired. Codex must only RECONCILE and run to exit 0.

- Verify the flags are computed (not hardcoded `True`) at the exact ranges above. They are.
- Note that F2 changes `delta_core` (via `kappa0_from_tube`), which feeds `nontrivial_compensated` (SymPy L202-204) / `nontrivialCompensated` (`.wl` L154-156). Since `kappa0_from_tube` simplifies to `(1+r_c)/3`, the residual still reduces to 0 and the flag stays `True` — confirm this after the F2 edit by running the scripts.

When recording, append:

```
## Applied: F4

- files_changed: none (flags already wired by a prior orchestrator-direct pass; F2's delta_core change flows through unchanged)
- summary: confirmed capstone survivability flags are computed from section-1–5 residuals (not hardcoded True); ran both scripts to exit 0
- deviation: already applied by prior orchestrator-direct pass — reconciled
```

If — contrary to the above — any survivability flag is found still hardcoded `True` (it should not be), wire it to the corresponding section residual per the original intent and record that under `## Applied: F4` instead with the edited lines listed.

**Verification command:**
After Codex applies F1-F3 and reconciles F4, run `python3 scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py` and `math -script mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`. Both must exit 0; the capstone must still print exactly one nontrivial survivor: `compensated Robin-mixed core realization` (SymPy `nontrivial_survivors` check L216-218 / `.wl` `nontrivialSurvivors` check L166-169).

## Applied: F4

- files_changed: none (flags already wired by a prior orchestrator-direct pass; F2's delta_core change flows through unchanged)
- summary: confirmed capstone survivability flags are computed from section-1–5 residuals (not hardcoded True); ran both scripts to exit 0
- deviation: already applied by prior orchestrator-direct pass — reconciled
