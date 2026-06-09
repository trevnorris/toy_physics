---
unit_id: 164
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage164_microscopic_log_channels.md]
  paper_appendix: present
---

# Audit unit 164 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_164.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage164_microscopic_log_channels.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 59, 265, 272-360 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage164_microscopic_log_channels_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage164_microscopic_log_channels_mathematica_audit.txt`

## What the paper claims

The stage card's `\stagefield{Output}` says the stage "Rewrites the Stage 163 off-family scalar in explicit branch-variable drifts of \((\mathcal Z_q,\rho_w,c_{s,w},c_s,\mathcal T_m,v_{w0},a,L_W)\)." The notes make this concrete in six steps: (1) the two Stage-163 log channels equal `δln(g/r)` and `-δln r_c` exactly via the parent ratios `r=λ/√(KsKq)`, `g=gq√Ks/(gs√Kq)`, `rc=λ²/(KsKq)`; (2) the exact Stage-118 parent-action products; (3) the uniform-overlap + D/N simplification; (4) the healing-locked products `gqKs/(gsλ) = -(27π mψ²/40ħμ₀q*)·Zq c_{s,w}³/(ρw Tm vw0 a² L_W²)` and `KsKq/λ² = (27π³mψ²/320ħμ₀q*²)·Zq cs² c_{s,w}³/(ρw vw0² a² L_W³)`; (5) the explicit `δ⊥` formula with coefficients `A_*=g_*+1/(4√(1+r_*²))`, `B_*=1/(2√(1+r_*²))`, `C_*=2g_*+3/(4√(1+r_*²))` and the Family-1 numerics (`A_*≈0.880589`, `B_*≈0.245108`, etc.); (6) the symbolic + numeric tangency law for `δln Tm`. The appendix (Eqs. app-part05-rg-defs, log-channels-rg, delta-perp-explicit, ABC-defs) mirrors steps 1, 4-grouped, and 5. The card is `\StatusExactClosure`, checkpoint: False — both engines required, exact alignment expected.

## What the script claims to verify

Both scripts verify the same six-step chain. SymPy asserts (5 `expect_zero`): the parent-ratio identities `g/r = gqKs/(gsλ)` and `1/rc = KsKq/λ²`; the two healing-locked product literals (after full general→uniform→healing substitution chain); and the compressed `δ⊥` grouped form. It prints the general and uniform products and all Family-1 numeric coefficients. Mathematica asserts the same five PLUS pins the uniform-overlap products against explicit literals (2 extra `expectZero`), PLUS an independent `Series`/`Coefficient` route deriving the log-channel coefficient vectors and `δ⊥` (3 `expectZero`), PLUS five `expectApprox` numeric checks on `A_*`, `B_*`, and the c_sw/v_w0/L_W coefficients (tol 1e-14).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Channels = δln(g/r), -δln rc (notes §1, app log-channels-rg) | py 41-42 / wl 57-58 expect_zero on parent-ratio identities | match |
| General Stage-118 products (notes §2) | py 64-65 print / wl 70-71 print (no assert, but feed asserted chain) | match (intermediate) |
| Uniform-overlap products (notes §3) | py 80-81 print / wl 93-94 expectZero against literals | match (wl stronger) |
| Healing-locked products (notes §4, app) | py 101-102 / wl 115-116 expect_zero vs literal | match |
| δ⊥ compressed A_*/B_*/C_* form (notes §5, app delta-perp-explicit) | py 130 / wl 183 expect_zero; wl 147-163 series route | match |
| Family-1 numerics A_*≈.880589, B_*≈.245108, coeff list (notes §5) | py 141-158 prints / wl 212-216 expectApprox | match |
| Tangency law δln Tm symbolic (notes §6, app §tangency) | py 133-134 solve+print / wl 186-190 Solve+print | match (symbolic) |
| Tangency law δln Tm NUMERICS (notes §6, 1.16167…, 3.48501…, etc.) | not emitted by either script | not a card deliverable — notes-only detail; not required |

`paper_alignment: aligned`. Every card/appendix deliverable has a script-side check; the notes §6 numeric tangency breakdown is notes-internal exposition (the card's stated Output is the symbolic branch-variable rewrite, not the numeric tangency table), so its absence from the scripts is not a `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 41 | `expect_zero(g/r - gqKs/(gsλ))` | channels=δln(g/r) | yes |
| A2 | sympy | 42 | `expect_zero(1/rc - KsKq/λ²)` | channels=-δln rc | yes |
| A3 | sympy | 101 | `expect_zero(first_prod_heal - literal)` | healing 1st product | yes |
| A4 | sympy | 102 | `expect_zero(second_prod_heal - literal)` | healing 2nd product | yes |
| A5 | sympy | 130 | `expect_zero(delta_perp - grouped)` | δ⊥ A_*/B_*/C_* form | yes |
| B1 | math | 57-58 | `expectZero` parent-ratio identities | channels | yes |
| B2 | math | 93-94 | `expectZero` uniform products vs literal | uniform products | yes |
| B3 | math | 115-116 | `expectZero` healing products vs literal | healing products | yes |
| B4 | math | 147-148 | `expectZero` series-route channel vectors | channel coeff vectors (indep route) | yes |
| B5 | math | 163 | `expectZero` series-route δ⊥ | δ⊥ form (indep route) | yes |
| B6 | math | 183 | `expectZero` δ⊥ compressed form | δ⊥ A_*/B_*/C_* form | yes |
| B7 | math | 212-216 | `expectApprox` 5 Family-1 numerics | Family-1 numeric coeffs | yes |

No tautological rows. A1/A2 verify a non-trivial √(KsKq) cancellation (definitional but load-bearing — the whole point of §1). A3-A4/B3 compare an independently-typed closed-form literal against the full substitution chain, so a typo in the literal would fail. B4-B5 use a genuinely distinct mechanism (multiplicative perturbation + Series + Coefficient).

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent re-derivation, not a transliteration** (high confidence). Three pieces of evidence:

1. **The .wl adds an entire derivation route the .py lacks.** Lines 118-166 ("Independent series route") derive the log-channel coefficient vectors by perturbing each microscopic variable multiplicatively (`p -> p (1 + eps dlnP)`), forming the perturbed/unperturbed product *ratio*, and reading the O(eps) `Coefficient[Normal[Series[...]]]`. The .py instead hand-writes `first_heal = dlnZ + 3*dlncsw - dlnrho - ...` (py 117-118) as a literal vector. These are mathematically equivalent but procedurally unrelated; a port would have copied the literal vector, not synthesized it from `Series`. I verified the route is sound: for `firstHealExpected` the variable exponents (zq:+1, csw:+3, rhoW:-1, tm:-1, vw0:-1, a:-2, lw:-2) reproduce `firstHealHand` exactly.

2. **The .wl pins intermediate forms the .py only prints.** The .py merely `print`s the uniform-overlap products (py 80-81); the .wl asserts them against explicit literals `firstUniformExpected = -9 ks zq/(64 lw² tm a⁴ ell² mu0 qstar vw0)` and reconciles via `expectZero` (wl 88-94). A line-by-line port would inherit the .py's print-only treatment, not add new assertions.

3. **The .wl adds numeric `expectApprox` gates (wl 212-216) absent from the .py**, which only prints the numbers.

Corresponding-section comparison (channel-vector derivation):
- .py 117-118: `first_heal = dlnZ + 3*dlncsw - dlnrho - dlnTm - dlnv - 2*dlna - 2*dlnLw` (literal)
- .wl 138-147: `firstRatio = (firstHealExpected /. pertRule)/firstHealExpected; firstHealSeries = Coefficient[Normal[Series[firstRatio,{eps,0,1}]],eps]; expectZero["first channel via series route", firstHealSeries - firstHealHand]` (synthesized, then cross-checked against the hand form)

The shared scaffolding (banner/expect helpers, same physical premises, same final A_*/B_*/C_* target) is expected and required — both engines must verify the *same claim*. What matters is the derivation path, and the .wl's is genuinely independent.

## Engine cross-check

Both engines agree at the level claimed:
- Parent-ratio identities → both `0` (py txt 5-6; wl txt 5-8).
- General products: py `-sqrt(2)*pi*K_s*Z_q/(2*I_sq*J_s*L_W**(3/2)*T_m*mu_0*qstar*v_w0)` ≡ wl `-(ks*Pi*zq)/(Sqrt[2]*isq*js*lw^(3/2)*mu0*qstar*tm*vw0)` (√2/2 = 1/√2; same form).
- Uniform products: py `-9*K_s*Z_q/(64*L_W**2*T_m*a**4*ell**2*...)` ≡ wl `(-9*ks*zq)/(64*a^4*ell^2*lw^2*...)`; second products both `9π².../512...`. Identical.
- Healing products: py `-27*pi*Z_q*c_sw**3*mpsi**2/(40*L_W**2*T_m*a**2*hbar*...)` ≡ wl `(-27*csw^3*mpsi^2*Pi*zq)/(40*a^2*hbar*lw^2*...)`; both healing exact-formula residuals `0`.
- δ⊥ compressed form: both residuals `0`.
- Family-1 numerics: py A_*=0.88058909004156951300, wl A_*=0.8805890900415696 (agree to printed precision; wl `expectApprox` diff 1.1e-16). B_*, c_sw, v_w0, L_W coefficients all agree (max diff 4.4e-16). The .py's printed coeff[ln(a)] = -1.7611781800831390260 vs wl -1.7611781800831392 agree to 1e-16.

No engine disagreement.

## Verdict justification

**clean.** I attacked: (a) the parent-ratio identities — they verify a real √(KsKq) cancellation, not a tautology; (b) the healing-product literals — they are independently-typed closed forms compared against a 3-step substitution chain, so they can fail; (c) the .wl as a possible port — it is demonstrably independent (adds a Series/Coefficient route and extra assertions the .py lacks); (d) the numerics — I re-derived A_*=g_*+1/(4√(1+r_*²))≈0.758035+0.122554=0.880589, B_*=2b≈0.245108, C_*=2g_*+3b≈1.88373, all matching both engines, the notes §5, and the appendix ABC-defs; (e) staleness — both outputs are newer than their scripts; (f) the stale-168/100π²/4107-radius class — not present in any stage-164 file (this stage doesn't carry that radius). The paper card Output, the notes §1-§5, and the appendix Eqs. all match the script-verified claim. The notes §6 numeric tangency table is notes-internal exposition beyond the card's stated symbolic Output and is not a required script deliverable. Both engines required (checkpoint:False, ExactClosure) and both present, substantive, and agreeing.

## Self-test notes

Checked the traps: (1) Variable independence / derivative-zero — N/A, no `diff`/`D` calls; the only "perturbation" is the .wl's multiplicative `Series` route, whose exponent-counting I verified by hand against both expected monomials. (2) Symmetry/parity — N/A, no integrals (all algebra is rational-function/log-linear). (3) Trivial-case pre-check — substituted the variable exponents into the series route and confirmed `firstHealSeries`/`secondHealSeries` reduce to the hand vectors; confirmed A_*/B_*/C_* numerics by direct arithmetic. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

Outputs are fresh (both `.txt` mtimes newer than their scripts), so the committed outputs are the authoritative record. All emitted deliverable values reconcile against the notes `.md` (the natural carrier for the full numeric/symbolic detail; the terse `.tex` card legitimately states only the Output sentence, and the appendix carries the grouped symbolic forms).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g/r = gqKs/(gsλ)`, `1/rc = KsKq/λ²` (parent-ratio identities) | py 41-42 / wl 57-58; py txt 5-6, wl txt 5-8 | notes §1 (lines 54-66); app log-channels-rg (293-301) | MATCH |
| general 1st product `-√2 π Ks Zq/(2 Isq Js Lw^{3/2} Tm μ₀ q* vw0)` | py 64 / wl 70; py txt 14, wl txt 13 | notes §2 (109-119, exact form before closure) | MATCH |
| general 2nd product `π² Ks Zq cs²/(4 Isq² Lw² μ₀ q*² vw0²)` | py 65 / wl 71; py txt 15, wl txt 14 | notes §2 (122-131) | MATCH |
| uniform 1st product `-9 Ks Zq/(64 Lw² Tm a⁴ ell² μ₀ q* vw0)` | py 80 / wl 88-91; py txt 20, wl txt 19 | notes §3 (172-181, exponent form) | MATCH |
| uniform 2nd product `9π² Ks Zq cs²/(512 Lw³ a⁴ ell² μ₀ q*² vw0²)` | py 81 / wl 89-92; py txt 21, wl txt 20 | notes §3 (184-193, exponent form) | MATCH |
| healing 1st product `-(27π mψ²/40ħμ₀q*) Zq c_sw³/(ρw Tm vw0 a² Lw²)` | py 98/101 / wl 110/115; py txt 26, wl txt 29 | notes §4 (223-230, box) | MATCH |
| healing 2nd product `(27π³mψ²/320ħμ₀q*²) Zq cs² c_sw³/(ρw vw0² a² Lw³)` | py 99/102 / wl 111/116; py txt 27, wl txt 30 | notes §4 (232-239, box) | MATCH |
| δ⊥ compressed A_*/B_*/C_* form | py 121-131 / wl 174-184; py txt 35, wl txt 51 | notes §5 (304-315, box); app delta-perp-explicit (306-327) | MATCH |
| tangency law δln Tm (symbolic) | py 133-134 / wl 186-190; py txt 36, wl txt 52 | notes §6 (372-400, box); app §tangency (328-332) | MATCH |
| g_* = 0.758035078944663 | py 141/148 / wl 194/200; py txt 41, wl txt 57 | notes 319, 344 | MATCH |
| 1/(4√(1+r_*²)) = 0.12255401109690651 | py 144/149 / wl 196/201; py txt 42 | notes 323-325 (≈0.122554011096907) | MATCH |
| A_* = 0.880589090041570 | py 145/150 / wl 197/212; py txt 43 | notes 328, 338 | MATCH |
| B_* = 0.245108022193814 | py 146/151 / wl 198/213; py txt 44 | notes 330, 342 | MATCH |
| coeff[ln c_sw] = 2.64176727012471 | py 153 / wl 205/214; py txt 46 | notes 340 | MATCH |
| coeff[ln Tm] = -0.758035078944663 | py 155 / wl 207; py txt 48 | notes 344 (=g_*) | MATCH |
| coeff[ln v_w0] = -1.00314310113848 | py 156 / wl 208/215; py txt 49 | notes 351 | MATCH |
| coeff[ln a] = -1.76117818008314 | py 157 / wl 209; py txt 50 | notes 353 | MATCH |
| coeff[ln L_W] = -1.88373219118005 (C_*) | py 158 / wl 210/216; py txt 51 | notes 355; app ABC-defs C_* (326) | MATCH |

INTERNAL scaffolding (no finding expected): `r_* = 1.77799353547498` (an imported Family-1 input, not a stage-164-derived deliverable — appears in notes 320/322 as the canonical point); the eps/pertRule perturbation symbols; series-route intermediate `firstRatio`/`secondRatio`; expectApprox diffs and PASS flags; banner strings.

Note: the notes §6 tangency-law NUMERIC breakdown (1.16167271174229, 3.48501813522686, 0.323345423484574, 1.32334542348457, 2.32334542348457, 2.48499999999999) is notes-internal — neither script emits these numbers, and the card's stated Output is the symbolic branch-variable rewrite, not the numeric tangency table. Not a deliverable absent from docs; it is present in the notes and simply not script-emitted. No finding.

reconciliation: complete; 19 deliverable values checked, 0 misaligned.
