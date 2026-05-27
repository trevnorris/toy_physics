---
unit_id: 089
batch: III.5
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
verification_status: scripts_pass_pending_verifier
needs_user_resolution: true
---

# Codex directive — unit 089

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit `paper/`, `notes/`, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch `paper/`, `notes/`, or any prose documents.

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_089.tex:13-17` quote: "At zero transport bias, \(\Omega_{\mathrm{Pe}=0}=1\), so Family--1 supplies \(\zeta_{\rm F1}(0)=A_{\rm F1}\simeq1.00005192880220.\)"
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_089.tex:26-29` quote: "\(\boxed{\mathrm{Pe}_{\rm req}=0}\)"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage089_family1_minimal_isotropic_verdict.md:86-109` quote: "Stage 62 proved that the explicit Family-1 transport map obeys `zeta_F1(Pe) = A_F1 Omega_Pe^2` … For the minimal isotropic passive/outgoing branch, `zeta_req^(min) = 1/3 < A_F1 ≈ 1.00005192880220`. So the required transport bias is exactly `Pe_req = 0`."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:91-92` quote: `if not (zeta_min < A_F1): raise AssertionError("Minimal isotropic branch no longer succeeds at zero transport bias.")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:83` quote: `expectTrue["minimal isotropic branch succeeds at zero transport bias", zetaMin < aF1];`

## Resolve before fix_loop

The paper card's `\stagefield{Output}` is the boxed `Pe_req = 0` (eq. app-stage089-Pe-zero). Its derivation hinges on `\Omega_{\rm Pe=0} = 1`, which gives `\zeta_{\rm F1}(0) = A_{\rm F1}`. The scripts verify only the precondition `zeta_min < A_F1` and never construct/check `Pe_req` or the limit `\Omega(Pe \to 0) = 1`. The Omega expression in both scripts is `0/0` at `Pe = 0`, so the limit is non-trivial. The script-side proof of the boxed equation has an unverified link.

Possible directions (the user picks one):
- (a) Strengthen scripts → in both `.py` and `.wl`, add an explicit symbolic check that `sp.limit(Omega, Pe, 0) == 1` / `Limit[omegaPe[pe], pe -> 0] == 1`, then assert `zeta_F1(0) - A_F1 == 0` (i.e., evaluate `zetaF1` at the limit, or equivalently `A_F1 * Omega(0)^2`), closing the chain to `Pe_req = 0`. Add an explicit `Pe_req = 0` print/assert wired to that derivation.
- (b) Accept the carry-forward — declare that `Omega(Pe=0) = 1` is established upstream (in Stage 062 audit scripts) and add an in-script comment naming the exact upstream verification file/line that proves it, so the link is anchored rather than asserted in-stage.
- (c) Update the paper card to weaken the `Output` line from `Pe_req = 0` to `zeta_req^min < A_F1` and demote the `Pe_req = 0` claim to a corollary stated in prose.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. Findings F2-F4 below are independent of this resolution and can be applied first if desired.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:34-89`

**Issue:**
The Mathematica script is a syntactic port of the SymPy script: identical constants (`kappaF1 = 12321/5`, `etaF1 = 37`), identical Pe literals (`peSuffChi = SetPrecision[96.5285247264386, 40]`, `peFailChi = SetPrecision[11220.5441626259, 40]`), identical Q definition, identical intermediate variable names (`zetaSuff = zetaF1[peSuffChi]`, `rhoSuff = q[zetaSuff, 0]`), and identical anchor checks. This violates the two-engine independence policy for checkpoint stages.

**Required change:**
Rewrite the Mathematica derivation so it does not echo SymPy's algebra. Specifically, in `wl:47-48`, remove the hardcoded `peSuffChi` and `peFailChi` literals and instead derive them by `FindRoot[zetaF1[pe] == zetaSuffTarget, {pe, 100}]` and `FindRoot[zetaF1[pe] == zetaFailTarget, {pe, 10000}]`, where `zetaSuffTarget` and `zetaFailTarget` come from the notes-quoted `rho_suff^(chi) ≈ 3.46622291347846` and `rho_fail^(chi) ≈ 3.46752913273870` via `zeta = rho - 1`. Then `rhoSuff` is computed from the rederived Pe. Also: change the `q[zeta, eps]` parameterization on `wl:49` to a different equivalent form, e.g., `q2[zeta_, eps_] := (1 - eps zeta + (1 - 2 eps) zeta)/(1 - eps zeta)` simplified to verify it matches numerically (an independent algebraic route to the same identity), and assert at `wl:57-59` against the algebraic form, not the trivial substitution.

Concretely, replace lines 47-48 with:
```
zetaSuffTarget = SetPrecision[3.46622291347846 - 1, 40];
zetaFailTarget = SetPrecision[3.46752913273870 - 1, 40];
peSuffChi = pe /. FindRoot[zetaF1[pe] == zetaSuffTarget, {pe, 100}, WorkingPrecision -> 40];
peFailChi = pe /. FindRoot[zetaF1[pe] == zetaFailTarget, {pe, 10000}, WorkingPrecision -> 40];
```
Then rebuild `zetaSuff`, `zetaFail`, `rhoSuff`, `rhoFail` from those rederived `peSuffChi`, `peFailChi`. The numeric values of `rhoSuff` and `rhoFail` should match SymPy's to ~10^-12.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 089` and confirm: (a) the literals `96.5285247264386` and `11220.5441626259` no longer appear in the `.wl` file (other than possibly as commented provenance); (b) script exits 0; (c) `rhoSuff` and `rhoFail` numerics still match the SymPy outputs to ~10^-12.

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:61-64`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:57-59`

**Issue:**
Lines `sympy:61-64` and `wl:57-59` perform anchor checks of the form `rho_X - (1 + zeta_X) == 0`. Because `rho_X` is defined as `Q.subs(zeta, zeta_X)` with `Q` already evaluating to `1 + zeta` at `eps_blk = 0` (sympy:53; wl:49), the residual is algebraically zero by construction. The sympy output transcript confirms all four such residuals print as `0` literally. These do not exercise Stage 69's content — they exercise sympy's substitution.

**Required change:**

In `sympy:61-64`, replace the four `expect_zero(...)` anchor lines with explicit numeric cross-checks against the notes' Stage-69 quoted values. Specifically:

```python
# Stage-69 numeric anchors (from notes/stages/moving_throat_pde_stage089_*.md §1).
expect_close = lambda name, val, target, tol=1e-12: (
    print(f"{name} = {sp.N(val, 25)} vs target {target}"),
    (_ for _ in ()).throw(AssertionError(f"{name} off target")) if abs(complex(sp.N(val - target, 50))) > tol else None
)
expect_close("rho_suff vs Stage-69 quote", rho_suff, sp.Float("3.46622291347846"), tol=1e-12)
expect_close("rho_fail vs Stage-69 quote", rho_fail, sp.Float("3.46752913273870"), tol=1e-12)
expect_close("rho_max  vs Stage-69 quote", rho_max,  sp.Float("3.46752922945601"), tol=1e-12)
```

(If introducing `expect_close` as a lambda is awkward, write a proper helper at module top similar to `expect_zero`. Keep the `Stage-62 zeta_max = A_F1 pi^2/4` line at `sympy:60` — that one is a legitimate Stage-62 algebraic identity.)

In `wl:57-59`, replace with the analogous `expectApprox` calls against the same notes-quoted targets:
```
expectApprox["rho_suff vs Stage-69 quote", rhoSuff, SetPrecision[3.46622291347846, 25], 10^-12];
expectApprox["rho_fail vs Stage-69 quote", rhoFail, SetPrecision[3.46752913273870, 25], 10^-12];
expectApprox["rho_max  vs Stage-69 quote", rhoMax,  SetPrecision[3.46752922945601, 25], 10^-12];
```

**Verification command:**
After Codex applies, verifier runs `redteam exec-sympy 089` and `redteam exec-mathematica 089`; both should still PASS with the new lines printing residual `~10^-15` rather than the trivial `0`.

## F4 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:49-50`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_mathematica_audit.wl:47-48`

**Issue:**
`Pe_suff_chi = sp.Float("96.5285247264386")` and `Pe_fail_chi = sp.Float("11220.5441626259")` are load-bearing literals without in-script provenance. The stage card lists `Inputs: zeta_req = 1/3, A_F1` but not these Pe values. The notes file quotes the derived `rho_suff^(chi) ≈ 3.46622291347846` but not the upstream Pe.

**Required change:**
If F2 is applied (which removes the literals from the `.wl` side by rederiving them), apply the parallel change to the SymPy script: replace `sympy:49-50` with

```python
# Stage-69 numeric anchors (rederived from notes §1 rho values).
zeta_suff_target = sp.Float("3.46622291347846") - 1
zeta_fail_target = sp.Float("3.46752913273870") - 1
Pe_suff_chi = sp.nsolve(zeta_F1 - zeta_suff_target, Pe, sp.Float("100"), tol=1e-30, maxsteps=200, prec=80)
Pe_fail_chi = sp.nsolve(zeta_F1 - zeta_fail_target, Pe, sp.Float("10000"), tol=1e-30, maxsteps=200, prec=80)
```

If F2 is NOT applied yet, add provenance comments above lines 49-50 in `.py` and 47-48 in `.wl`:

```python
# CARRY-FORWARD: Pe_suff_chi and Pe_fail_chi from Stage-069 verification.
# Source: scripts/output/moving_throat_pde_stage069_*_sympy_audit.txt
```

Pick the derivation path if Codex can confirm the FindRoot/nsolve converges; otherwise pick the provenance-comment path.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 089` and confirms: (a) either no bare literals at sympy:49-50, or a comment block naming the upstream source; (b) script exits 0; (c) `rho_suff`, `rho_fail` numerics still match notes-quoted Stage-69 values.
