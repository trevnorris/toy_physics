---
unit_id: 014
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 014 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_014.tex`
- notes: (none)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_014.tex`) describes a "mouth-Taylor gate bridge" whose deliverables are three identities. First, the even-gate definitions `K_1 = -z_2 - z_0/9` (eq:stage014-k1) and `H_{\rm even} = -z_4 + (2/3) z_2 - z_0/27` (eq:stage014-heven). Second, the determinantal claim that the even-gate compensation denominator is proportional to `Q S_2 - \Delta H_{\rm port}` (eq:stage007-compensation-denominator). The card's `Output` paragraph then states verbatim: "Stage~014 exports the mechanism sieve: projected EM mouth data can tune the even gates only away from the degeneracy \eqref{eq:stage007-compensation-denominator}." The Part I appendix row corroborates this with the one-line summary "Gate conditions for carrying mouth-local projected data into the grouped bundle." There is no `notes/stages/moving_throat_pde_stage014_*.md` file under the project tree.

## What the script claims to verify

The SymPy script (and its Mathematica twin) verifies several things. (1) Structural linearity of the Taylor lift: every bundle slip (`z2d, z4d, n2d, n4d`) is linear in each primitive slip (sympy lines 102-104). (2) Bundle pull-back consistency: directly-substituted `Xi_load`, `K_1`, `H_{\rm even}` equal the symbolically lifted bundle forms (sympy lines 121, 123, 125; mathematica lines 113-115). (3) Sector independence: `Xi_load` is independent of `Sx, Hx, Gx`, and `K_1, H_{\rm even}` are independent of `Px, Gx` (sympy lines 112-116; mathematica lines 128-136). (4) The mechanism sieve: solving `K_1 = H_{\rm even} = 0` over the source-only `(Q', \Delta')` or spectral-only `(S_2', H_{\rm port}')` sectors yields only the trivial solution (sympy lines 135-138, mathematica M5/M6), and the mixed-sector compensation surface has denominator proportional to `\Delta H_{\rm port} - Q S_2`, equivalently `-\Delta^2 \cdot Z_2` (sympy lines 131-132, 144-146; mathematica M9). (5) Constant-prefactor transport ledger: derivative coefficients of `\delta P_2`, `\delta P_4` with respect to `G_W'` (sympy lines 127-128; mathematica M3, M4). (6) Sign-flip mutation test on the compensation denominator (sympy line 133; mathematica M10).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check(s) | Status |
|---|---|---|
| `K_1 = -z_2 - z_0/9` (eq:stage014-k1) | `K1_bundle` consistency line 123 sympy, line 114 mathematica | match |
| `H_{\rm even} = -z_4 + (2/3) z_2 - z_0/27` (eq:stage014-heven) | `He_bundle` consistency line 125 sympy, line 115 mathematica | match |
| Compensation denominator proportional to `Q S_2 - \Delta H_{\rm port}` (eq:stage007-compensation-denominator) | `Hx_den - 9*Delta**2*(Delta*Hport - Q*S2)` line 131 sympy + Z2-slot equivalence line 145; M9 mathematica | match (sign-conjugate factor: `\Delta H_{\rm port} - Q S_2 = -(Q S_2 - \Delta H_{\rm port})`; same vanishing surface) |
| "Mechanism sieve: even gates only tunable away from the degeneracy" | Pure-sector solves trivial lines 135-138 sympy / M5, M6 mathematica; Jacobian determinants nonzero (M7, M8); compensation surface explicit (sympy line 80, mathematica line 165) | match |
| (no paper claim) | `Xi_load` definition line 67, `Xi_bundle` consistency line 121, `d Xi_load / dPx` line 126, M2 mathematica | extra |
| (no paper claim) | `\delta P_2`, `\delta P_4` definitions lines 74-75, `d\delta P_2/dGx` line 127, `d\delta P_4/dGx` line 128, M3 + M4 mathematica | extra |
| (no paper claim) | `Compat` line 72 (printed only, no assertion) | extra (informational) |

`paper_alignment` is set to `partial` because every paper-side deliverable has a matching, well-anchored script-side check, BUT the scripts exercise additional load-bearing content (`Xi_load`, `\delta P_2`, `\delta P_4`) that the paper card and appendix row do not name.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 102-104 | `assert_zero(sp.diff(_expr, _slip, 2))` for {z2d, z4d, n2d, n4d} × six slips | structural lift linearity (paper-internal pre-requisite) | yes |
| A2 | sympy | 112-113 | `assert_zero(diff(Xi, S2x/Hx/Gx))` | extra (Xi_load not on paper card) | yes (against script docstring) |
| A3 | sympy | 114-116 | `assert_zero(diff(K1/He, Px/Gx))` | claim 1 (K1, H_even sector structure) — auxiliary | yes |
| A4 | sympy | 121 | `assert_zero(Xi - Xi_bundle)` | extra | yes |
| A5 | sympy | 123 | `assert_zero(K1 - K1_bundle)` | claim 1 (eq:stage014-k1) | yes |
| A6 | sympy | 125 | `assert_zero(He - He_bundle)` | claim 1 (eq:stage014-heven) | yes |
| A7 | sympy | 126 | `assert_zero(coef_Xi_Px - 2/P)` | extra | yes |
| A8 | sympy | 127 | `assert_zero(coef_dP2_Gx + 2*P/(D0*Delta**2))` | extra | yes |
| A9 | sympy | 128 | `assert_nonzero(coef_dP4_Gx)` | extra | yes |
| A10 | sympy | 129-130 | `assert_nonzero(qd/sh Jacobian det)` | claim 2 (mechanism sieve) | yes |
| A11 | sympy | 131-132 | `assert_zero(Hx_den - 9*Delta**2*(Delta*Hport - Q*S2))` and S2 analogue | claim 2 (compensation denominator) | yes |
| A12 | sympy | 133 | `assert_nonzero(Hx_den - sign-flip mutation)` | claim 2 (sign anchoring) | yes |
| A13 | sympy | 135-138 | bare `if qd_only != [...]: raise` and likewise for sh | claim 2 (mechanism-sieve negative-existence) | yes |
| A14 | sympy | 142-143 | `assert_zero(K1.subs(comp_surface))` and H_even analogue | claim 2 (compensation closes the gates) | **no — self-flagged tautology** |
| A15 | sympy | 145-146 | `assert_zero(Hx_den - (-9*Delta**4*Z2_slot))` and Sx analogue | claim 2 (factor identification with Z2 slot) | yes |
| M1 | mathematica | 128-136 | `expectZero` independence checks (mirror of A2/A3) | claim 1 / extra | yes |
| M2 | mathematica | 138 | `expectZero[D[XiLoad, Px] - 2/P]` | extra | yes |
| M3 | mathematica | 139 | `expectZero[D[deltaP2, Gx] - (-2*P/(D0*Delta^2))]` | extra | yes |
| M4 | mathematica | 140 | `expectNonzero[D[deltaP4, Gx]]` | extra | yes |
| M5 | mathematica | 142-147 | `expectEqual[qdSolve, {{Qx -> 0, Deltax -> 0}}]` | claim 2 | yes |
| M6 | mathematica | 149-154 | `expectEqual[shSolve, {{Sx -> 0, Hx -> 0}}]` | claim 2 | yes |
| M7 | mathematica | 158, 161 | `expectNonzero[Det[jacobianQD]]` | claim 2 | yes |
| M8 | mathematica | 159, 162 | `expectNonzero[Det[jacobianSH]]` | claim 2 | yes |
| M9 | mathematica | 168-175 | `expectZero` on compensation denominators | claim 2 (mirror of A11) | yes |
| M10 | mathematica | 176-179 | `expectNonzero` on sign-flip mutation | claim 2 (mirror of A12) | yes |

Notes: A14 is self-flagged on sympy lines 139-141 ("tautological by sp.solve's contract... kept for visual symmetry"). A2/A3 are self-flagged on line 106 ("construction-level, not physical") with the bundle pull-back checks A4-A6 serving as the substantive anchor. A4-A6 are genuine: they equate the directly-substituted `Xi/K1/He` (built from closed-form primitive z2/z4/n0 expressions) to the symbolically-lifted `Xi_bundle/K1_bundle/He_bundle` (built from `z0d, z2d, z4d, n0d`). Both engines confirm zero residue.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_014.tex:30-33` (the entire `Output` paragraph)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:67,71-75,93-95,126-128`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_mathematica_audit.wl:113,117-126,138-140`

**What's wrong:**

The paper card's `Output` paragraph (stage_014.tex:30-33) states only:
> "Stage~014 exports the mechanism sieve: projected EM mouth data can tune the even gates only away from the degeneracy \eqref{eq:stage007-compensation-denominator}."

The card body exhibits exactly two equations (`K_1`, `H_{\rm even}`) and one denominator identity. The card does **not** mention `Xi_load`, `\delta P_2`, `\delta P_4`, the isotropic compatibility (`Compat`), or any "5PN bottleneck" / "constant-prefactor transport" content.

The scripts, by contrast, devote a substantial fraction of their load-bearing assertions to these extra quantities:

- Sympy line 67: `Xi = sp.simplify((2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)).subs(subs_der)/mu1)`
- Sympy line 126: `assert_zero("d Xi_load / d Pprime", coef_Xi_Px - 2/P)`
- Sympy lines 74-75: definitions of `deltaP2`, `deltaP4` as bundle-level transport expressions
- Sympy lines 127-128: assertions on `d \delta P_2 / dGx` and `d \delta P_4 / dGx`
- Mathematica line 113: `XiLoad = reduce[n0d/N0[P, Delta] + z0d/D0]`
- Mathematica lines 117-126: definitions of `deltaP2`, `deltaP4`
- Mathematica lines 138-140: assertions M2, M3, M4

The SymPy docstring frames this as "moving-throat 5PN bottlenecks", but that framing is absent from the paper card. The appendix one-liner ("Gate conditions for carrying mouth-local projected data into the grouped bundle") is broad enough to *plausibly* admit the transport ledger, but does not actually name `Xi_load`, `\delta P_2`, or `\delta P_4`.

**Why this matters:**

If the paper card is the canonical definition of what Stage 014 proves, then a significant fraction of the script's load-bearing assertions test content the paper does not claim. A reader who compares paper-to-script will not be able to verify the Stage 014 claim from the paper card alone, because the card omits half of what the script tests. Conversely, if the script is the canonical record, the paper card understates the stage's deliverables. Either way the two are out of sync.

**Required change:**

(Routed to user — see directive.) Codex must NOT silently edit either side. The user chooses whether to (a) expand the paper card's `Output` to enumerate `Xi_load`, `\delta P_2`, `\delta P_4` as further deliverables of Stage 014, (b) trim the scripts to test only the even-gate sieve content the paper actually states (and migrate the transport-ledger assertions to a separate stage), or (c) add one or two sentences to the paper card explicitly acknowledging that the audit script also exercises the constant-prefactor transport coefficients used downstream.

**Verification:**

After user resolution and any follow-up edits, the audit re-run should produce a Paper ↔ script cross-check table with no `extra` rows (or with `extra` rows explicitly justified by a paper-side acknowledgement).

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py:139-143`

**What's wrong:**

Lines 142-143 read:
```
assert_zero("compensation K1", K1.subs(comp_surface))
assert_zero("compensation H_even", He.subs(comp_surface))
```
where `comp_surface` is defined on line 80 as `sp.solve([sp.Eq(K1, 0), sp.Eq(He, 0)], [Hx, Sx], dict=True)[0]`. By `sp.solve`'s contract, substituting the returned solution back into the original system yields zero by construction; the assertions cannot fail no matter what `K1` and `H_even` are. The author has self-flagged this on lines 139-141 ("Note: the next two assertions are tautological by sp.solve's contract... kept here for visual symmetry with the explicit denominator factorization assertions above").

The Mathematica script does not contain the analogous tautological pair (it stops after `compSurface` is computed and goes straight to denominator factorization checks).

**Why this matters:**

Self-flagged tautological assertions clutter the assertion inventory and obscure which checks do real work. A future reader who greps for `assert_zero(... compensation ...)` and sees these pass will believe the compensation surface has been verified, when in fact the substantive verification is the denominator factorization on lines 131-132 and the Z2-slot identification on lines 145-146, plus the sign-anchor mutation on line 133. Keeping tautologies "for visual symmetry" is exactly the trap the audit framework is supposed to surface.

**Required change:**

In `scripts/moving_throat_pde_stage014_projected_maxwell_mouth_taylor_gate_bridge_sympy_audit.py`, delete lines 139-143 inclusive (the two `assert_zero` calls for `compensation K1` and `compensation H_even` plus the preceding three-line self-flag comment). Keep lines 144-146 (the `Z2_slot` denominator identification) intact — those are the substantive checks. The output transcript will lose two `compensation K1 ...` / `compensation H_even ...` status print lines that the current script does not emit anyway (`assert_zero` does not print on success; the assertions raise only on failure).

**Verification:**

After Codex applies, the sympy script should be shorter by exactly 5 lines (lines 139-143 removed). `redteam exec-sympy 014` should still exit 0 with output transcript identical to the current one (no print lines change since `assert_zero` is silent on success). No other assertion's status changes.

## Independent-derivation check (Mathematica)

The Mathematica script is **NOT** a transliteration of the SymPy script. Its derivation strategy is genuinely independent:

- SymPy (lines 48-54) defines `z0, z2, z4, n0, n2, n4` as already-product-rule-expanded closed-form expressions:
  ```python
  z2 = (-Delta**2*h1 + Delta*(Hport*d1 + Q*s1 + S2*q1) - 2*Q*S2*d1)/Delta**3
  ```
- Mathematica (lines 52-58) defines the *primitive* `Z0, Z2, Z4, N0, N2, N4` as functions of their underlying source/spectral arguments, then obtains the Taylor lift via symbolic differentiation `D[Z2[Q + mu1*Qx*ell, ...], ell] /. ell -> 0`:
  ```mathematica
  Z2[Q_, S2_, Hport_, Delta_] := (Q*S2 - Hport*Delta)/Delta^2;
  z2d = reduce[(D[Z2[Q + mu1*Qx*ell, S2 + mu1*Sx*ell, Hport + mu1*Hx*ell, Delta + mu1*Deltax*ell], ell] /. ell -> 0)/mu1];
  ```

These are mathematically equivalent (both yield the linearized z2 in the primitive slips), but they take *different* algebraic routes (closed-form pre-expanded product rule vs. symbolic chain rule applied to an unfactored primitive). A line-by-line transliteration would have copied `(-Delta**2*h1 + Delta*(Hport*d1 + Q*s1 + S2*q1) - 2*Q*S2*d1)/Delta**3` directly into Mathematica syntax; this one does not.

The same pattern holds for the bottlenecks: sympy defines `K1, He` by substituting into the closed-form z2/z4 (sympy lines 68-69), while Mathematica builds K1, He from the symbolically-lifted `z2d, z4d` (mathematica lines 114-115). The independent-derivation requirement is satisfied.

## Engine cross-check

Both engines agree at all checkpoints:

| Item | SymPy | Mathematica |
|---|---|---|
| `d(Xi_load)/dP'` | `2/P` (sympy output line 56) | M2 residual = 0 against `2/P` |
| `d(\delta P_2)/dG_W'` | `-2*P/(D0*Delta**2)` (sympy output line 58) | M3 residual = 0 against `-2*P/(D0*Delta^2)` |
| `d(\delta P_4)/dG_W'` | `2*(D0*Delta*Gw - 2*D0*P*S2 + 2*D2*Delta*P)/(D0**2*Delta**3)` (sympy output line 60) | M4 residual = `(2*D0*Delta*Gw + 4*D2*Delta*P - 4*D0*P*S2)/(D0^2*Delta^3)` — identical up to factor-of-2 distribution |
| `Hx` compensation denominator | `9*Delta**2*(Delta*Hport - Q*S2)` (sympy line 131) | M9 residual = 0 against `9*Delta^2*(Delta*Hport - Q*S2)` |
| `Sx` compensation denominator | `9*Delta*(Delta*Hport - Q*S2)` (sympy line 132) | M9 residual = 0 against `9*Delta*(Delta*Hport - Q*S2)` |
| Mechanism sieve pure-sector solves | `[{Deltax: 0, Qx: 0}]`, `[{Hx: 0, S2x: 0}]` (sympy output lines 67, 70) | M5 = `{{Qx -> 0, Deltax -> 0}}`, M6 = `{{Sx -> 0, Hx -> 0}}` |
| Sign-flip mutation | nonzero (asserted) | M10 residual = `-18*Delta^2*Q*S2` (nonzero) |
| Spectral Jacobian determinant | nonzero (asserted) | M8 = `(-(Delta*Hport) + Q*S2)/Delta^4` = `-(\Delta H_{port} - Q S_2)/\Delta^4` |

No `engine_disagreement` finding.

Output freshness: sympy script mtime 21 May 12:53, output mtime 21 May 15:00 (fresh). Mathematica script mtime 21 May 12:54, output mtime 21 May 15:01 (fresh). No `stale_output` finding.

## Verdict justification

Within the scope of what the paper card actually claims, the scripts hold up under adversarial attack. The even-gate definitions match the paper exactly via the bundle-pull-back consistency assertions A5/A6 (sympy lines 123, 125) and their Mathematica counterparts. The compensation denominator is correctly identified as the `Q S_2 - \Delta H_{\rm port}` surface (with the sign convention explicitly anchored against a sign-flip mutation in A12). The mechanism-sieve claim is exercised both by direct pure-sector solves and by the Jacobian-determinant non-vanishing checks. Mathematica is a genuine second derivation, not a transliteration. The only script-side defect on the substantive content is the self-flagged tautological pair A14 (low severity, F2).

The principal concern is **F1**: the scripts test additional quantities (`Xi_load`, `\delta P_2`, `\delta P_4`, and a printed-only `Compat`) that the paper card does not mention. Per the v2 audit contract this is `paper_misalignment` (subtype `paper_missing_script_claim`) and requires user resolution; Codex is not authorized to edit the paper or trim the scripts unilaterally. The verdict is `findings` (not `clean`) because of F1+F2, but not `stop_cold` because F1 is low-severity (the scripts over-cover the paper, they do not contradict it) and F2 is a minor scaffolding cleanup. Downstream stages presumably consume `Xi_load, \delta P_2, \delta P_4` already, so this does not invalidate downstream content.

## Self-test notes

I checked: (1) variable-independence trap on the proposed F2 deletion — removing the two tautological asserts and their preceding comment block has no side effects (the variables `comp_surface, K1, He` remain defined and are used nowhere else after that point). (2) F1 does not propose a Codex-applied edit, so symmetry/parity/trivial-case traps are not load-bearing for it; the F2 edit is purely a deletion so its trivial-case substitution is "no new check needed". (3) Path specifications are not load-bearing since no `missing_verification_script` finding is filed (both engines are present). (4) Paper round-trip: the F2 deletion does not introduce any new paper-misalignment since it removes (rather than adds) script content; the paper card is silent on whether `K1.subs(comp_surface)` should be checked, so removing the check leaves the paper-script relationship unchanged. (5) I verified the apparent sign discrepancy between paper's `Q S_2 - \Delta H_{\rm port}` and script's `\Delta H_{\rm port} - Q S_2` is a sign-conjugate factor (both define the same vanishing surface), and the script explicitly bridges the two via the `Z2_slot` identification on lines 144-146.
