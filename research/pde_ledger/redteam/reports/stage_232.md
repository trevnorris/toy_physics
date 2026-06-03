---
unit_id: 232
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 232 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_232.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card `\stagefield{Output}` states: "Current branch verdict: the support/source side is safe by a large margin in the injected data; the active unresolved gate is static orbit-lock/coherent placement rather than support starvation." The stage is a numerical injection step (`\StatusNumerical{}`, `\StatusOpen{}`): it takes the numerically located Family-1 5PN support/source packets and inserts them into the same-charge selected-branch inequalities. The notes (Sections 1-7) enumerate the concrete deliverables the audit script is said to verify (notes Section 5): (1) the refreshed `\Lambda_{\rm EM}` geometry formulas `Lambda_ell = 20√2 π/x01`, `kappa = (9/5) Lambda_ell²`; (2) the Robin support equation `y tan y = eta` and ceiling `zeta_max = A_K π²/4`; (3) the support-drop endpoints `Delta_0`, `Delta_inf`; (4) the fixed-point law `Pe = Xi·Delta(Pe;kappa,eta)` for both wall-depth extractions; (5) the roots `Pe_*^(chi)`, `Pe_*^(J)`; (6) the transported support/source values `zeta_phys`, `rho_alpha,max`; (7) the demand constants `zeta_req = 1/3`, `rho_alpha^req = 4/3`; (8) proximity of both points to the ceiling. The appendix row (line 76) and narrative (lines 867-874) confirm `zeta_req = 1/3`, `rho_alpha^req = 4/3`, `zeta_phys ≈ 2.4675`, large-margin overshoot. The notes also write the figure-of-merit prefactor explicitly: `Xi_chi = 168 Theta_w^(chi) Lambda_ell²` and `Xi_J = 168 Theta_w^(J) Lambda_ell²` (notes lines 153, 157).

## What the script claims to verify

The script reconstructs the geometry constants, solves the Robin equation by bisection, computes `A_K` and `zeta_max`, evaluates the support-drop endpoints in closed form and cross-checks a closed-form `Delta_closed(Pe)` against direct `mp.quad` quadrature at Pe=1,10,100, solves the two fixed-point roots by bisection inside the `[Xi·Delta_0, Xi·Delta_inf]` bracket, computes the overlap boost `Omega_Pe` and the physical support ratios, then asserts the injected safety margins, ratios, and ceiling gaps are positive and equal to the quoted decimals. Each numeric quantity is checked with `assert_close(actual, <independent decimal literal>, tol)`, plus three structural asserts (`margin>0`, `gap>0`, root-in-bracket and residual≈0). The carried wall-depth inputs `Theta_w_chi`, `Theta_w_J` and all `assert_close` targets are hardcoded numerical injections — legitimate for an explicitly "injection" stage. The single substantive disagreement is the figure-of-merit prefactor: script line 149-150 uses `100 * Theta_w * Lambda_ell**2`, not the `168` written in the notes.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| (1) `Lambda_ell = 20√2π/x01`, `kappa = (9/5)Lambda_ell²` | lines 49-64 + asserts 62-64 | match |
| (2) Robin `y tan y = eta`, `zeta_max = A_K π²/4` | lines 69-83 | match |
| (3) `Delta_0`, `Delta_inf` | lines 106-114 + quad cross-check 131-141 | match |
| (4) fixed-point `Pe = Xi·Delta` (both extractions) | lines 162-184 | match (but Xi prefactor disputed, see F1) |
| (5) roots `Pe_*^(chi)`, `Pe_*^(J)` | asserts 186-187 | match numerically (depend on F1) |
| (6) `zeta_phys`, `rho_alpha,max` | asserts 209-212 | match |
| (7) `zeta_req = 1/3`, `rho_alpha^req = 4/3` | lines 217-218 + asserts | match |
| (8) ceiling proximity | lines 230-231, asserts 249-258 | match |
| Notes formula `Xi = 168 Theta_w Lambda_ell²` | line 149-150 uses `100` | mismatch |

Dominant pattern is `match`, with one load-bearing constant `mismatch` between the notes' written formula and the script → `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62-64 | `assert_close(Lambda_ell/chi_s/kappa, literal)` | claim 1 | yes |
| A2 | sympy | 81-83 | `assert_close(y/A_K/zeta_max, literal)` | claim 2 | yes |
| A3 | sympy | 113-114 | `assert_close(Delta_0/Delta_inf, literal)` | claim 3 | yes |
| A4 | sympy | 131-135 | `Delta_quad vs Delta_closed at Pe=1,10,100` | claim 3/4 (closed-form validity) | yes (independent quadrature) |
| A5 | sympy | 138-139 | `Delta_closed(1e-8)≈Delta_0`, `Delta_0<Delta_closed(1000)<Delta_inf` | claim 3 (endpoint theorem) | yes |
| A6 | sympy | 156-157 | `assert_close(Xi_chi/Xi_J, literal)` | claim 4 | partial (matches notes *decimals*, not notes *formula* — F1) |
| A7 | sympy | 181-184 | root-in-bracket + `|Pe - Xi·Delta(Pe)|<1e-30` | claim 4/5 | yes (non-tautological fixed-point residual) |
| A8 | sympy | 186-187 | `assert_close(Pe_chi/Pe_J, literal)` | claim 5 | yes (depends on F1) |
| A9 | sympy | 209-212 | `assert_close(zeta_phys/rho, literal)` | claim 6 | yes |
| A10 | sympy | 249-251 | `margin>0`, `gap_to_ceiling>0` | claim 7/8 (verdict: safe margin) | yes |
| A11 | sympy | 253-258 | `assert_close(margins/ratios, literal)` | claim 7 | yes |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md:153,157`
- `/var/projects/toy_projects/.../scripts/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.py:149-150`

**What's wrong:**
The notes write the wall/source figure-of-merit with prefactor **168**:
- notes:153 `\Xi_\chi = 168\,\Theta_w^{(\chi)}\Lambda_\ell^2 \approx 5.5548332017764099\times 10^5`
- notes:157 `\Xi_J = 168\,\Theta_w^{(J)}\Lambda_\ell^2 \approx 1.2663707072528143\times 10^5`

The script uses prefactor **100**:
- script:149 `Xi_chi = 100 * Theta_w_chi * Lambda_ell**2`
- script:150 `Xi_J = 100 * Theta_w_J * Lambda_ell**2`

The notes are internally inconsistent: their *written formula* (168) does not produce their *quoted decimal*. With `Lambda_ell² ≈ 1365.2826610556` and `Theta_w^(chi) ≈ 4.06863235008162`:
- `100 · 4.06863235008162 · 1365.2826610556 ≈ 555483.320178` → matches the quoted `≈ 5.5548332017764099×10^5` and the script assert at line 156.
- `168 · 4.06863235008162 · 1365.2826610556 ≈ 933211.98` → contradicts the quoted decimal.
- For the J branch: `100·...≈ 126637.070725` matches the quoted `1.2663707072528143×10^5`; `168·...≈ 212750.28` does not.

So the script's prefactor (100) is consistent with the notes' own *numbers*, and the notes' *symbolic formula* (168) is the outlier. This is a genuine notes↔script (and notes-internal) disagreement on a load-bearing constant that feeds the fixed-point roots `Pe_*` and downstream margins; it must be resolved by the user (Codex does not edit notes, and the auditor cannot unilaterally pick which side is "the answer").

**Why this matters:**
`Xi` is the coefficient of the fixed-point law `Pe = Xi·Delta(Pe)`. A wrong prefactor changes `Pe_*` (and the asserted decimals at lines 186-187), and in principle every downstream margin/ratio/ceiling-gap literal. The qualitative verdict is robust to this choice because `Omega_Pe → π/2` as `Pe → ∞`, so `zeta_phys → A_K·π²/4 = zeta_max` regardless of whether `Pe_* ≈ 11155` or larger — the support side saturates near the ceiling either way. But the published *numbers* (`Pe_*`, ceiling gaps `~1e-7`/`~1e-6`) are only correct for one prefactor. Leaving the notes' formula reading 168 while the verified pipeline uses 100 means the paper-side derivation chain cannot be reproduced from the notes as written.

**Required change:**
None applied by Codex. Routed to user — see directive `## Resolve before fix_loop`.

**Verification:**
After user resolution: if 100 is correct, the notes formula at lines 153/157 is corrected to `100` (paper-side edit, post-resolution) and the script is unchanged. If 168 is correct, the script lines 149-150 change to `168 *` and every asserted decimal at lines 156-157, 186-187, 209-212, 253-258, and the ceiling-gap prints must be regenerated; the verifier re-runs `redteam exec-sympy 232` and confirms exit 0 with the new literals.

### F2 — missing_verification_script (subtype missing_mathematica)

**Severity:** medium
**Files:**
- (no `.wl` exists for unit 232)

**What's wrong:**
This unit has a SymPy script and no Mathematica `.wl`. The stage card itself records "Mathematica audit: none yet" (`stage_232.tex:11`). Every claim the SymPy script verifies is independently reproducible with native Mathematica primitives: `BesselJZero[0,1]` for `x01`; closed-form geometry; `FindRoot`/`Reduce` for the Robin root `y tan y = eta` on `(0,π/2)`; `Integrate`/`NIntegrate` of `K(x)·Sigma(x)` over `[0,1]` for `Delta` (an *independent* route to the script's hand-derived `Delta_closed`); `FindRoot` for the fixed-point `Pe = Xi·Delta(Pe)`; and exact-arithmetic margin checks. The dual-engine rule's test is "is it possible," not "is it necessary" — it is clearly possible. A `mathematica_audit.wl` should be written that derives `Delta(Pe)` by *symbolic integration of the kernel against the source profile* (rather than transliterating the `Delta_closed` algebra) and independently solves the fixed-point and margin claims.

**Why this matters:**
Second-engine coverage is the project's standing requirement wherever Mathematica can independently verify. The SymPy `Delta_closed` is a hand-derived closed form; an independent Mathematica `Integrate` of the kernel is exactly the kind of cross-check that would catch an algebra error in that closed form (the SymPy script only cross-checks `Delta_closed` against its own `mp.quad` of the same integrand, which would not catch a wrong analytic primitive if the integrand itself were mis-stated).

**Required change:**
Codex writes a NEW independent-route Mathematica script (see directive F2 claim manifest). Note: do NOT apply F2 until F1 is resolved, because the `Xi` prefactor (100 vs 168) determines the `Pe_*` and margin literals the `.wl` must reproduce. If F1 resolves to "script is correct (100)," the `.wl` targets the current literals; if "168," both engines target the regenerated literals.

**Verification:**
After Codex applies, the verifier runs `math -script mathematica/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.wl` and confirms it exits 0 with all in-file checks passing, and that the final `Delta`, `Pe_*`, `zeta_phys`, and margin values agree with the (resolved) SymPy values.

## Independent-derivation check (Mathematica)

No `.wl` exists, so no transliteration check applies. See F2 for the requirement that a new `.wl` derive `Delta(Pe)` via native symbolic integration rather than porting the SymPy `Delta_closed` algebra.

## Engine cross-check

Only one engine present; no cross-check possible. This is the basis for F2.

## Verdict justification

The SymPy script is structurally sound: the geometry, Robin root, endpoint, fixed-point, overlap, and margin computations are all genuine numeric reconstructions, the `assert_close` targets are independent decimal literals (a wrong computation would fail them — non-tautological), and the `Delta_quad` vs `Delta_closed` and fixed-point-residual checks are real cross-checks. The carried `Theta_w` inputs and asserted decimals being hardcoded is legitimate for an explicitly numerical "injection" stage. The verdict is `findings` for two reasons: (F1) the notes' written figure-of-merit prefactor (168) disagrees with both the script (100) and the notes' own quoted decimals (which 100 reproduces) — a load-bearing constant disagreement routed to the user; and (F2) the dual-engine rule requires a Mathematica audit, which is clearly possible here and absent. Not `CRITICAL_DOWNSTREAM`: this stage emits a qualitative verdict (support side non-bottlenecked), not a numeric constant consumed by later stages, and the verdict is robust to the prefactor choice because `zeta_phys` saturates at `zeta_max` regardless. Attacks tried that failed: probed every `assert_close` for tautology (none — all are independent targets); checked for division-by-zero at the `p±alpha` poles in `Delta_closed` over the actually-evaluated `Pe` range (none hit); verified the Robin bracket sign convention `(f_lo<0, f_hi>0)` is satisfied and the root lands in `(0,π/2)`; verified the `Sigma_Pe` and `Omega_Pe` numerically-stable forms are algebraically identical to the notes' forms (they are).

## Self-test notes

Checked the four relevant traps. (1) No `sp.diff`/`D` calls in this script, so the variable-independence trap does not apply. (2) Symmetry/parity: the `Delta` integral is over `[0,1]` (bounded, not symmetric), so the odd/even-vanishing trap does not apply; instead I confirmed the closed form's `p±alpha` denominators never vanish in the evaluated `Pe` range. (3) Trivial-case pre-check: I hand-verified the prefactor arithmetic (100 vs 168) against the quoted decimals, which is the crux of F1, and confirmed 100 reproduces both `Xi_chi` and `Xi_J`. (4) Path: the F2 `.wl` target is placed under `mathematica/` with the mandatory `_mathematica_audit.wl` suffix; I did not prescribe the script route (Codex designs it), only the claim manifest.
