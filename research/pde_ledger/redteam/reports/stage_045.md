---
unit_id: 045
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md"]
  paper_appendix: present
---

# Audit unit 045 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_045.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 68 plus surrounding context rows)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.txt`

## What the paper claims

Paper card `\stagefield{Output}` (verbatim): "The exact tracking condition \eqref{eq:app-stage045-tracking} and coherent tracking factor \eqref{eq:app-stage045-Rtr}." The boxed body equations enumerate three explicit deliverables: (i) the coherent kernel identity `g_B g_R = g_W g_S` (eq:app-stage045-coherent-condition); (ii) the tracking surface `R_phi = R_U = R_tr` (eq:app-stage045-tracking); (iii) the closed-form tracking factor `R_tr = (1 + chi_0/(1+delta_U))/(1+chi_0)` with bounds `1/(1+delta_U) < R_tr < 1` (eq:app-stage045-Rtr). The body also asserts the rank-2 branch collapses to tracking laws `M_tr = G_tr(xi, delta; R_tr)` and `R_target = F_tr(xi, delta; R_tr)` (eq:app-stage045-tracking-laws). The notes add explicit closed forms for `G_tr` (at lam0=2/9), `F_tr` (at lam0=2/9), the total baseline `M_tr = M_mix + M_supp`, and the assertion that the Stage-27 normalization residual collapses to the Stage-23 tracking form `F_tr`. The part-03 appendix row at line 68 paraphrases: "Coherent condition forcing exact tracking, R_phi=R_U=R_tr." Paper alignment is OK; the script covers all paper-side deliverables and adds none beyond what notes mention.

## What the script claims to verify

The SymPy docstring lists four checks: (1) the coherent kernel identity `g_B g_R = g_W g_S`; (2) `rho_0 = sigma_0`; (3) the R_tr range identities; (4) collapse of the Stage-27 quadratic branch equation to the tracking law. The actual assertions cover those four plus a check that `M_tr = M_mix + M_supp` matches a channel-summed expression, plus a check that the script's hand-written generic-lam0 expression for `F_track` specializes to the notes' `F_tr` at lam0=2/9. The Mathematica script asserts the same set of identities in the same order, using identical algebraic constructions. Banners and the final print line label this as "Stage 28" / "Stage 028" while the filename and paper card use "stage 045".

## Paper ↔ script cross-check

| Paper / notes deliverable | Script-side check | Status |
|---|---|---|
| `g_B g_R = g_W g_S` from coherent kernel (eq-coherent-condition) | sympy line 68 `expect_zero("g_B g_R - g_W g_S", g_B_ext * g_R_ext - g_W_ext * g_S_ext)`; wl line 57 | match |
| `R_phi = R_U = R_tr` (eq-tracking) | indirectly verified via `rho_0 - sigma_0 = 0` (sympy 74, wl 60) — same direction factor follows | match |
| Closed-form `R_tr = (1+chi_0/(1+delta_U))/(1+chi_0)` + range bounds (eq-Rtr) | range identities, sympy 94 / 98, wl 63 / 67 | match |
| `M_tr = G_tr(xi, delta; R_tr)` (eq-tracking-laws + notes sec 5) | tracking quadratic collapse + `G_tr formula` + D/N specialization (sympy 160, 167, 176; wl 113, 118, 123) | match |
| `R_target = F_tr(xi, delta; R_tr)` (eq-tracking-laws + notes sec 6: "Stage-27 normalization residual collapses to the Stage-23 tracking form") | sympy 188 / wl 135 only verify that a hand-written generic-lam0 F_track specializes algebraically to the notes' lam0=2/9 form; no derivation of F_tr from the Stage-27 residual is performed | partial |
| Total baseline `M_tr = M_mix + M_supp` (notes sec 4) | sympy 129 / wl 91 check `M_tr - channel_sum`, but both sides are built from identical formulas via the same `prefactor * Z / (1-eps)` template — tautological | partial |

Set `paper_alignment: aligned` — script covers every paper deliverable and does not introduce claims beyond what notes/appendix mention. Two deliverables are only partially exercised (insufficient_verification rather than misalignment).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63 | `g_W extracted - reference == 0` | g_W coefficient form (kernel definition) | yes (kernel coefficient check) |
| A2 | sympy | 64 | `g_R extracted - reference == 0` | g_R coefficient form | yes |
| A3 | sympy | 65 | `g_B extracted - reference == 0` | g_B coefficient form | yes |
| A4 | sympy | 66 | `g_S extracted - reference == 0` | g_S coefficient form | yes |
| A5 | sympy | 68 | `g_B g_R - g_W g_S == 0` | claim (i) coherent kernel identity | yes |
| A6 | sympy | 74 | `rho_0 - sigma_0 == 0` | claim (ii) R_phi = R_U | yes (the simplified equality is non-trivial: the two ratios use different sqrt(mu) products) |
| A7 | sympy | 94 | `(1 - R_tr) - chi_0 delta_U / ((1+chi_0)(1+delta_U)) == 0` | claim (iii) range identity 1 | yes |
| A8 | sympy | 98 | `R_tr - 1/(1+delta_U) - delta_U / ((1+chi_0)(1+delta_U)) == 0` | claim (iii) range identity 2 | yes |
| A9 | sympy | 129 | `M_tr - channel_sum == 0` | total baseline M_tr = M_mix + M_supp (notes sec 4) | no — both sides built from identical template |
| A10 | sympy | 160 | tracking quadratic collapse `numTrack + collapsed_num == 0` | claim (iv) Stage-27 quadratic reduces to tracking quadratic | yes |
| A11 | sympy | 167 | `M_tr_req - xi(delta+xi)/(delta + (1+lam0 R_U^2) xi) == 0` | G_tr generic formula | yes |
| A12 | sympy | 176 | `G_tr D/N specialization == 0` | G_tr at lam0=2/9 matches notes | yes (substitution check at the load-bearing physical lam0) |
| A13 | sympy | 188 | `F_track at lam0=2/9 - F_tr_expected == 0` | F_tr formula from Stage-27 normalization residual | partial — both sides hand-written by the script; the residual-collapse claim is not exercised |
| A14 | wl | 48–51 | g_* extracted vs reference | mirrors A1–A4 | yes |
| A15 | wl | 57 | `g_B g_R - g_W g_S == 0` | mirrors A5 | yes |
| A16 | wl | 60 | `rho_0 - sigma_0 == 0` | mirrors A6 | yes |
| A17 | wl | 63 | range identity 1 | mirrors A7 | yes |
| A18 | wl | 67 | range identity 2 | mirrors A8 | yes |
| A19 | wl | 91 | M_tr - channel_sum | mirrors A9 | no (same template tautology) |
| A20 | wl | 113 | tracking quadratic collapse | mirrors A10 | yes |
| A21 | wl | 118 | G_tr generic | mirrors A11 | yes |
| A22 | wl | 123 | G_tr D/N | mirrors A12 | yes |
| A23 | wl | 135 | F_tr law | mirrors A13 | partial |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:34-69`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:46-101`

**What's wrong:**
The Mathematica script is a line-by-line port of the SymPy script's algebra, not an independent re-derivation of the stage's claims from the physical premises. Corresponding sections match step for step:

- SymPy (line 47): `coupling_density = sp.expand((lam_W*W_sym + lam_phi*phi_sym) * (eta_sym - gamma*U_sym))`
- Mathematica (line 35): `couplingDensity = Expand[(lamW*Wsym + lamPhi*phisym)*(etasym - gamma*Usym)]`

- SymPy (lines 49–56): four `coupling_density.coeff(...)` extractions producing `g_W_ext`, `g_R_ext`, `g_B_ext`, `g_S_ext` with identical sign conventions.
- Mathematica (lines 36–43): four `Coefficient[Coefficient[...]]` extractions producing `gWext`, `gRext`, `gBext`, `gSext` with the same sign conventions.

- SymPy (lines 114–129): `channels = [("W", Z_W, eps_W_split), ("phi", Z_phi, eps_phi_split)]` and `M_tr_channel_sum = sum(prefactor * Z_i / (1 - eps_i) for ...)`.
- Mathematica (lines 78–83): `channels = {{ZW, epsWSplit}, {ZPhi, epsPhiSplit}}; mTrChannelSum = FullSimplify[Total[prefactor*#[[1]]/(1 - #[[2]]) & /@ channels], ...]`.

- SymPy (lines 178–188): F_track defined explicitly as `(delta + (1 + lam0*R_U^2)*xi)^2 * (delta + (1 + lam0*R_U)*xi)^2 / ((1 - xi) * ((delta + xi)^2 + lam0*R_U^2*xi^2)^2)`.
- Mathematica (lines 125–135): identical expression `fTrack = FullSimplify[(delta + (1 + lambda0*rU^2)*xi)^2 * (delta + (1 + lambda0*rU)*xi)^2 / ((1 - xi) * ((delta + xi)^2 + lambda0*rU^2*xi^2)^2), ...]`.

The only structural divergence is the Mathematica `Series[branchNumRaw, {rPhi, rU, 0}]` route in lines 104–105 (a marginal independence), but every other check uses identical variable choreography and identical hand-coded expressions for F_track, channel structure, branch equation, and G_tr.

**Why this matters:**
The second-engine policy exists so that two engines independently confirm the same physics. A transliteration only confirms that two CAS tools agree on textbook-level algebra — it does not provide a second route to the physical claim. If the SymPy algebra encodes a subtle convention error (e.g., the F_track ansatz at generic lam0 is wrong, or the branch equation has a sign error inherited from Stage-27), the Mathematica script reproduces the same error rather than catching it.

**Required change:**
Restructure the Mathematica script so it derives the same identities from a different starting point. Concrete options that preserve the existing check set:
1. For the F_tr check, instead of hand-writing `fTrack` and verifying its specialization, derive F_tr at lam0=2/9 directly from the Stage-27 normalization residual (notes section 6) by substituting `rPhi -> rU` in the residual and confirming the result equals the closed-form F_tr quoted in the notes — this is the verification the paper actually claims (cf. F3) AND it makes the wl path genuinely different from the py path because the SymPy F_track is hand-written and the Mathematica route would be derivation-based.
2. For the branch equation: instead of mirroring the SymPy `Series` workaround, substitute `rPhi -> rU` in the simplified `branchEq` numerator and verify via independent factoring (`Factor[Numerator[Together[branchEq /. rPhi -> rU]]]`); compare against the collapsed quadratic.
3. For the coherent identity `g_B g_R = g_W g_S`: instead of mirroring the manual coefficient extraction, exploit the factored form directly: `expectZero["g_B g_R - g_W g_S", (lamPhi/Sqrt[muEta*muPhi]) * (gamma*lamW/Sqrt[muU*muW]) - (lamW/Sqrt[muEta*muW]) * (gamma*lamPhi/Sqrt[muU*muPhi])]` so the wl version uses reference forms while the py version uses extracted forms — establishing independence at that check.

Treat option (1) as the primary fix because it also resolves F3.

**Verification:**
After the rewrite, a structural diff between `.py` and `.wl` should not show a one-to-one correspondence of intermediate variable names and hand-written expressions. The `wl` output should still pass the same set of `expectZero` assertions but with at least one assertion derived by a route the SymPy script does not use.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:114-129`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:78-91`

**What's wrong:**
The `M_tr - channel_sum` assertion compares two expressions built from the same template. SymPy:
```
channels = [("W", Z_W, eps_W_split), ("phi", Z_phi, eps_phi_split)]
prefactor = 8 * (1 + chi_0)**2 / (sp.pi**2 * (1 - eps_eta))
M_tr_channel_sum = sum(prefactor * Z_i / (1 - eps_i) for (_, Z_i, eps_i) in channels)
M_mix = prefactor * Z_W / (1 - eps_W_split)
M_supp = prefactor * Z_phi / (1 - eps_phi_split)
M_tr = M_mix + M_supp
expect_zero("M_tr - channel_sum", M_tr - M_tr_channel_sum)
```
Both sides expand to `prefactor * (Z_W/(1-eps_W_split) + Z_phi/(1-eps_phi_split))` by construction. The assertion cannot fail. The notes claim being checked here is `M_tr := M_mix + M_supp` with the explicit prefactor form (notes sec 4) — but neither `M_mix` nor `M_supp` is derived from `Z_W (1+chi_0)^2 / (pi^2 (1-eps_eta) (1-eps_W^split))` anywhere upstream; both are written directly. So the check confirms the algebra of sum-over-list vs sum-of-two-terms, not the physical content (that the M_mix/M_supp formulas have the prefactor structure the notes state).

**Why this matters:**
A reviewer reading the transcript sees "M_tr - channel_sum = 0" and concludes the total baseline formula is verified. It isn't: it's verified that two ways of writing the same sum agree. If the prefactor were wrong (e.g., `8` replaced by `4`, or `(1+chi_0)^2` replaced by `(1+chi_0)`), the check would still pass with the wrong value because both sides change identically.

**Required change:**
Remove the tautological assertion (since it adds no information). Concretely: delete the `M_tr_channel_sum` definition (sympy line 119–121; wl line 80–83) and the `expect_zero("M_tr - channel_sum", ...)` line (sympy line 129; wl line 91). Keep the printed `M_mix`, `M_supp`, `M_tr` lines for documentation. No replacement is required by the paper — the substantive total-baseline claim is downstream of Stage-022 and Stage-026, not Stage-045.

**Verification:**
After the fix, the transcripts no longer contain the tautological `M_tr - channel_sum = 0` line; `M_mix`, `M_supp`, `M_tr` are still printed. Section 3 of the script still runs but no longer claims a non-existent verification.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:172-189`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:125-136`

**What's wrong:**
Notes section 6 claims: "Because the coherent local kernel lands on the tracking surface, the Stage-27 normalization residual also collapses to the Stage-23 tracking form: `R_target = F_tr(xi, delta; R_tr)`" with the explicit closed-form `F_tr` at lam0=2/9. The script's check is:

```python
F_track = sp.simplify(
    (delta + (1 + lam0 * R_U ** 2) * xi) ** 2
    * (delta + (1 + lam0 * R_U) * xi) ** 2
    / ((1 - xi) * ((delta + xi) ** 2 + lam0 * R_U ** 2 * xi ** 2) ** 2)
)
F_tr_expected = sp.simplify(
    (9 * delta + (9 + 2 * R_U ** 2) * xi) ** 2
    * (9 * delta + (9 + 2 * R_U) * xi) ** 2
    / (81 * (1 - xi) * (9 * delta ** 2 + 18 * delta * xi + (9 + 2 * R_U ** 2) * xi ** 2) ** 2)
)
expect_zero("F_tr normalization law", sp.simplify(F_track.subs(lam0, lam0_dn)) - F_tr_expected)
```

`F_track` is a generic-lam0 expression hand-written into the script with no upstream derivation. `F_tr_expected` is the lam0=2/9 form from the notes. The assertion verifies that algebraic substitution `lam0 → 2/9` in the script's hand-written F_track gives the notes' F_tr — that is a pure algebraic identity by construction (the generic form was written so it would specialize this way). The `coherent normalization residual = R_target - F_tr_expected` line is only printed (sympy 189; wl 136), never asserted to vanish.

Nothing in the script ties `F_track` to the Stage-27 normalization residual on the tracking branch (which is what the notes claim collapses to F_tr). The script does not import Stage-27's residual, substitute `R_phi -> R_U`, and confirm the result equals F_tr. So the deliverable "`R_target = F_tr` on the tracking branch" is not exercised; the script only verifies that two equivalent forms of the same hand-written expression agree.

**Why this matters:**
The notes treat the normalization-residual collapse as a substantive Stage-045 result ("the first exact 'concrete-kernel' form of the normalization test"). The script as written cannot detect an error in either (a) the F_tr formula or (b) the claim that Stage-27's normalization residual reduces to F_tr on the tracking branch. A user could change F_track to a wrong generic-lam0 expression and the F_tr_expected to the matching wrong specialization, and the check would still pass.

**Required change:**
Replace the F_tr check with a derivation-style verification rooted in the Stage-27 normalization residual. Concrete plan:

(i) The Stage-27 normalization residual on the rank-2 closure is the upstream object that this stage claims collapses to F_tr. The Stage-27 / Stage-044 audit scripts contain this residual; the directive should ask the user to point at the canonical Stage-27 expression rather than have Codex guess. The substantive check is:

```python
# pseudocode
R_target_residual_stage27 = <load or restate from Stage-27 / Stage-044>
collapsed_residual = sp.simplify(R_target_residual_stage27
                                 .subs(R_phi, R_U)
                                 .subs(lam0, sp.Rational(2, 9)))
F_tr_expected = sp.simplify(
    (9*delta + (9 + 2*R_U**2)*xi)**2 * (9*delta + (9 + 2*R_U)*xi)**2
    / (81 * (1 - xi) * (9*delta**2 + 18*delta*xi + (9 + 2*R_U**2)*xi**2)**2)
)
expect_zero("F_tr collapse from Stage-27 residual",
            collapsed_residual - (R_target - F_tr_expected))
```

The Mathematica fix should derive F_tr by an independent route (e.g., by directly substituting tracking conditions into the explicit branch problem; or by verifying that F_tr_expected as a function of `rU` satisfies the normalization equation the residual encodes). See F1 — the F_tr fix is also the primary independence fix.

(ii) Because reconstructing the Stage-27 residual without help may be risky and the user's prior audit work on Stages 027 / 044 should be the source of truth, the directive's F3 entry includes a `## Resolve before fix_loop` block asking the user to identify the canonical Stage-27 residual expression to import. The orchestrator routes this to the user before Codex applies the fix.

**Verification:**
After the fix, the script's output contains an assertion line of the form `F_tr collapse from Stage-27 residual = 0` whose two sides are derived from genuinely different expressions (one from the Stage-27 residual on the tracking branch; one from the notes' closed form). Removing the tracking substitution `R_phi -> R_U` (or substituting a wrong tracking surface) should cause the assertion to fail.

### F4 — paper_misalignment

**Subtype:** notes_contradicts_script

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_045.tex:1` quote: "Stage 045: Coherent Local D/N Support Kernel and Exact Tracking Reduction"
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex:68` quote: "045 & Coherent local D/N support kernel"
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:3,31,191` quotes: `"""Stage 28 SymPy audit."""`, `banner("STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT")`, `print("\nAll Stage-28 symbolic checks passed.")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:26,139` quotes: `banner["STAGE 028 — COHERENT LOCAL TRACKING"]`, `Print["Stage 045 Mathematica audit passed."]`

**What's wrong:**
The paper card and part-03 appendix label this stage `045`. The SymPy script's docstring, banner, and terminal print line still label it "Stage 28" (the pre-renumbering identifier). The Mathematica banner labels it "STAGE 028" while the final print says "Stage 045 Mathematica audit passed." — internally inconsistent within the same file.

This is a labeling drift, not a math error. The mathematical content is aligned. But because it persists in user-facing transcript output (the `.txt` files in `output/` quote "STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT" / "STAGE 028 — COHERENT LOCAL TRACKING"), it can confuse readers tracing transcripts back to paper stages.

**Why this matters:**
A reviewer or future auditor reading `scripts/output/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.txt` sees "STAGE 28" and may believe they have the wrong transcript. The paper card's `\stagefield{Verification}` block points to this transcript file. The fix direction is unambiguous (script side should match the paper), so this paper_misalignment does not require a user resolution gate — see directive F4.

**Required change:**
Pure cosmetic fix: update banner / docstring / final print strings to say "Stage 045" instead of "Stage 28" / "Stage 028".
- sympy line 3: replace `Stage 28 SymPy audit.` with `Stage 045 SymPy audit.`
- sympy line 31: replace `STAGE 28 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT` with `STAGE 045 — COHERENT LOCAL D/N KERNEL TRACKING AUDIT`
- sympy line 191: replace `All Stage-28 symbolic checks passed.` with `All Stage-045 symbolic checks passed.`
- wl line 26: replace `STAGE 028 — COHERENT LOCAL TRACKING` with `STAGE 045 — COHERENT LOCAL TRACKING`

(The wl line 139 footer `Stage 045 Mathematica audit passed.` is already correct and should remain as is.)

No assertion depends on the label strings, so the fix is mechanically safe.

**Verification:**
After the fix, both `.txt` outputs banner with "STAGE 045" and the final lines say "Stage 045 ... passed." The numeric content of assertions is unchanged.

## Independent-derivation check (Mathematica)

Transliteration. See F1 for the detailed quote pairs. The variable-by-variable correspondence between the two scripts is one-to-one, the channel-list construction is identical, the F_track hand-written expression is identical, and the branch equation is identical. The single divergence — Mathematica's `Series[branchNumRaw, {rPhi, rU, 0}]` vs SymPy's `branch_eq.subs(R_phi, R_U)` — is a minor algorithmic variation in how the substitution is performed, not an independent derivation route.

## Engine cross-check

Both engines produce identical residuals across all corresponding checks:

| Check | SymPy output | Mathematica output |
|---|---|---|
| g_* extracted - reference | 0 (lines 9–12) | 0 (lines 5–12) |
| g_B g_R - g_W g_S | 0 (line 13) | 0 (line 14) |
| rho_0 - sigma_0 | 0 (line 16) | 0 (line 18) |
| range identity 1, 2 | 0 (lines 24–25) | 0 (lines 20–23) |
| M_tr - channel_sum | 0 (line 34) | 0 (line 29) |
| tracking quadratic collapse | 0 (line 41) | 0 (line 33) |
| G_tr formula | 0 (line 43) | 0 (line 36) |
| G_tr D/N specialization | 0 (line 44) | 0 (line 38) |
| F_tr normalization law | 0 (line 45) | 0 (line 40) |
| coherent normalization residual | nonzero printed expression (line 46) | nonzero printed expression (line 41) |

Engines agree. (The "coherent normalization residual" lines are not asserted to vanish — they are status prints that show `R_target - F_tr_expected`. This is consistent across engines but does not constitute a verification.)

## Verdict justification

Verdict: **findings** (4). Paper alignment is OK: every paper-side deliverable maps to at least one script-side check. The script's math content is consistent with the paper's claims, and both engines agree on every assertion. The findings are:

- **F1 (medium)** — the Mathematica script is a line-by-line transliteration of the SymPy script's algebra rather than an independent derivation, violating the second-engine policy.
- **F2 (low)** — the `M_tr - channel_sum` assertion is tautological (both sides built from the same template), so it does not actually verify the total-baseline formula.
- **F3 (medium)** — the F_tr normalization-law check is insufficient: F_track is hand-written into the script and the assertion only verifies its algebraic specialization at lam0=2/9, not the deliverable claim that the Stage-27 normalization residual collapses to F_tr on the tracking branch.
- **F4 (low, paper_misalignment / notes_contradicts_script)** — script banners and the SymPy docstring still say "Stage 28" / "Stage 028" while the paper card and filenames say "Stage 045"; the Mathematica file is internally inconsistent (banner "STAGE 028", footer "Stage 045 ... passed.").

Attacks tried that failed: (i) the coherent kernel identity check — extraction from the factored coupling density does produce identifiable g_*, and the resulting identity `g_B g_R = g_W g_S` follows substantively from the factored form (the factorization is the physical input, but the algebra extracting and re-combining the coefficients is non-trivial in that it would surface a sign error in the extraction); (ii) the `rho_0 - sigma_0` check — at face value this looks tautological after `g_B g_R = g_W g_S`, but the two ratios use different `sqrt(mu)` structures so the simplification is non-trivial; (iii) the range identities for R_tr — direct algebraic identities that are not pre-baked anywhere in the script; (iv) the G_tr D/N specialization — substitutes lam0=2/9 into a generically-derived M_tr_req solved from the collapsed quadratic, so the specialization checks a real derivation chain rather than self-referential algebra.

No `stop_cold` flag. F4 is a labeling drift (low severity) and the fix direction is unambiguous (script side updates to match the paper card); the directive applies the cosmetic fix without requiring user resolution. F3 does require user resolution to identify the canonical Stage-27 normalization residual expression to import; the directive includes a `## Resolve before fix_loop` block for F3 only. None of the findings would invalidate downstream units' results (the math is right; the framing isn't perfect).

## Self-test notes

I checked the following traps. (a) Variable independence: no `sp.diff` calls in the script; not applicable. (b) Symmetry/parity: no integrals over unbounded domains in the script; not applicable. (c) Trivial-case pre-check for the proposed F3 fix: substituting `R_phi -> R_U` in Stage-27's residual at lam0=2/9 should yield F_tr per notes sec 6 — this is exactly the deliverable the script currently fails to exercise; if the Stage-27 residual or the F_tr formula were wrong, the new assertion would fail. (d) Path specifications: all paths in this report are absolute and match the existing file layout. (e) Paper round-trip: F1, F2, F4 do not change paper-side claims. F3's required change uses notes section 6 verbatim and the Stage-27 residual from upstream scripts, so it does not introduce a new paper_misalignment.
