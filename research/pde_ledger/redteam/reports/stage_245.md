---
unit_id: 245
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 245 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_245.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (row at line 88)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

Stage 245 compiles the abstract Stage-243 non-rigid mouth/dressing `(U,V)` lane into the actual finite physical observables. The carried chart is `U = ln(T^2/T_ref^2)`, `V = ln(eps_eta/eps_eta_ref)`, with the selected-branch identity `R_target * T^2 = Lambda_0 (1 - eps_eta)` (eq:app-part08-stage245-selected-identity) producing the target-ratio compiler `R_target/R_ref = [(1 - eps_ref e^V)/(1 - eps_ref)] e^{-U}`. The reduced free energy `F_nr = (1/2)a_U U^2 + (1/2)a_V V^2 - chi_UV U V - f_U U` has the exact stationary packet `U = a_V f_U/Delta_UV`, `V = chi_UV f_U/Delta_UV` with `Delta_UV = a_U a_V - chi_UV^2`, admissible when `Delta_UV > 0`. Distinct deliverables: (1) the exact stationary `(U,V)` packet, response ratio `V/U = chi_UV/a_V`, and determinant/Hessian admissibility; (2) the finite physical compiler for `T^2`, `eps_eta`, and `R_target` — the `R_target` leg being derived from the selected-branch identity; (3) the positive `U/V` drain `D_UV = chi_UV U V = chi_UV^2 a_V f_U^2/Delta_UV^2 >= 0` (nonnegativity is the stated claim, eq:app-part08-stage245... "it is nonnegative even when U and V have opposite signs"); (4) the dependent microscopic correction `y_nr^dep = (0, -chi_UV f_U/Delta_UV, (a_V - chi_UV) f_U/Delta_UV)`; (5) the weak-axisymmetric first-order packet `Xi_1^nr = u_1`, `R_1^nr + Xi_1^nr = -eps_ref/(1-eps_ref) v_1`, `R_1^nr = -u_1 - eps_ref/(1-eps_ref) v_1`; (6) the orbit-side/support-blind split `partial_Lambda U = partial_varrho U = partial_Lambda V = partial_varrho V = 0`, explicitly conditioned on `f_U` and `chi_UV` being carried as orbit-side data "rather than support-placement coordinates." The appendix row (line 88) summarizes the same deliverables and the `StatusExactClosure` flag. The notes additionally record a Session-I readback (`U ~ 0.14313458`, `V ~ -0.03619791`, `eps_ref = 0.3`, `eps_eta ~ 0.28933482`, `R_target/R_ref ~ 0.87984149`, `R_1^nr ~ -0.12762119`).

## What the script claims to verify

The docstring enumerates 8 deliverables matching the paper. Section 1 solves the stationarity equations with `sp.solve` and compares to the closed forms (`U_expected`, `V_expected`), checks `V/U`, the Hessian determinant, and both limiting cases (`f_U=0`, `chi_UV=0`) — these are genuine, independent derivations. Section 2 computes the drain closed form and matches it to `D_expected`. Section 3 defines `T^2`, `eps_eta`, `R_ratio` by the boxed formulas and checks they reproduce `exp(U_sol)`, `exp(V_sol)`, and a multiplicative round-trip. Section 4 checks the dependent correction matrix against its closed form. Section 5 linearizes `U = eps u_1`, `V = eps v_1` and extracts first-order coefficients of `ln T^2`, `ln(1-eps_eta)`, `ln R_target`. Section 6 differentiates each compiled object w.r.t. support symbols `Lam`, `varrho` and asserts zero. Section 7 reconstructs `(chi_UV, Delta, f_U)` from the observed `(U,V)` then "rebuilds" `(U,V)` and checks the numeric `eps_eta` against the session value.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) stationary packet `U_sol,V_sol`, `V/U`, det, admissibility | §1 solve + closed-form compare + limit cases | match |
| (2a) `T^2 = T_ref^2 e^U`, `eps_eta = eps_ref e^V` | §3 `T2/T2_ref - exp(U_sol)`, `eps_eta/eps_ref - exp(V_sol)` | partial (tautological — see F2) |
| (2b) `R_target/R_ref` FROM selected-branch identity `R_target T^2 = Lambda_0(1-eps_eta)` | §3 round-trip of self-defined `R_ratio` | mismatch (selected-branch identity never used — see F2) |
| (3) drain `D_UV >= 0` (nonnegativity) | §2 closed-form match only; no `>=0` assertion | partial (see F4) |
| (4) dependent correction `y_nr^dep` | §4 matrix vs closed form | match |
| (5) first-order `Xi_1^nr`, `R_1^nr`, `R_1+Xi_1` | §5 series coefficients | match |
| (6) support-blind `partial_{Lam,varrho} = 0` | §6 `diff(obj, Lam/varrho)` of expressions never containing Lam/varrho | mismatch (vacuous self-test trap — see F1) |
| Session-I readback `eps_eta, R_target/R_ref, R_1` | §7: only `eps_eta` asserted; `R_target/R_ref`, `R_1` printed only; `(U,V)` "rebuild" is an inverse round-trip | partial (see F3) |

`paper_alignment` set to `partial`: no value disagreements, but two deliverables (support-blindness and the `R_target` selected-branch leg) are not actually exercised by the script's assertions.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 69 | `simplify(U_sol - U_expected) == 0` | claim 1 (U) | yes |
| A2 | sympy | 70 | `simplify(V_sol - V_expected) == 0` | claim 1 (V) | yes |
| A3 | sympy | 71 | `simplify(ratio_VU - chi_UV/a_V) == 0` | claim 1 (V/U) | yes |
| A4 | sympy | 72 | `simplify(det_H - Delta_UV) == 0` | claim 1 (admissibility) | yes |
| A5 | sympy | 73 | `U_sol.subs(f_U,0) == 0` | claim 1 (limit 2.1) | yes |
| A6 | sympy | 74 | `V_sol.subs(f_U,0) == 0` | claim 1 (limit 2.1) | yes |
| A7 | sympy | 75 | `simplify(V_sol.subs(chi_UV,0)) == 0` | claim 1 (limit 2.2) | yes |
| A8 | sympy | 76 | `simplify(U_sol.subs(chi_UV,0) - f_U/a_U) == 0` | claim 1 (limit 2.2) | yes |
| A9 | sympy | 88 | `simplify(D_UV - D_expected) == 0` | claim 3 (drain form) | partial (form only, not `>=0`) |
| A10 | sympy | 107 | `simplify(T2/T2_ref - exp(U_sol)) == 0` | claim 2a | no (tautological by definition) |
| A11 | sympy | 108 | `simplify(eps_eta/eps_eta_ref - exp(V_sol)) == 0` | claim 2a | no (tautological by definition) |
| A12 | sympy | 109 | `simplify(R_exact_check) == 0` | claim 2b | no (round-trip of self-defined R_ratio) |
| A13 | sympy | 125 | `simplify(y_dep - y_expected) == Matrix([0,0,0])` | claim 4 | yes |
| A14 | sympy | 151 | `simplify(Xi1_nr - u1) == 0` | claim 5 | partial (log∘exp trivial, but matches claim) |
| A15 | sympy | 152 | `simplify(RpXi1 - RpXi1_expected) == 0` | claim 5 | yes |
| A16 | sympy | 153 | `simplify(R1_nr - R1_expected) == 0` | claim 5 | yes |
| A17 | sympy | 154 | `simplify(R1_nr + Xi1_nr - RpXi1) == 0` | claim 5 (consistency) | yes |
| A18 | sympy | 177-178 | `diff(obj, Lam)==0`, `diff(obj, varrho)==0` (×7 objects) | claim 6 | no (vacuous — VAR not in EXPR) |
| A19 | sympy | 216 | `abs(eps_eta_obs - 0.28933482) < 5e-9` | Session-I readback | yes (numeric anchor) |
| A20 | sympy | 217 | `abs(U_rebuilt - U_obs) < 1e-12` | Session-I readback | no (inverse round-trip) |
| A21 | sympy | 218 | `abs(V_rebuilt - V_obs) < 1e-12` | Session-I readback | no (inverse round-trip) |

## Findings

### F1 — insufficient_verification (variable-independence self-test trap)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:159-178`

**What's wrong:**
Section 6 claims to verify the paper's support-blindness deliverable (eq:app-part08-stage245-support-blind, lines 99-107: `partial_Lambda U = partial_varrho U = partial_Lambda V = partial_varrho V = 0`). The script introduces fresh symbols `Lam, varrho = sp.symbols("Lam varrho", ...)` at line 160, then for each compiled object (`U_sol`, `V_sol`, `eps_eta`, `R_ratio`, `D_UV`, `Delta_Keta`, `Delta_mu`) asserts `sp.diff(obj, Lam) == 0` and `sp.diff(obj, varrho) == 0` (lines 177-178). But every one of these objects is built only from `a_U, a_V, chi_UV, f_U, eps_eta_ref` (and constants). None of them ever contains `Lam` or `varrho`. Therefore `sp.diff(obj, Lam)` is identically zero **by construction** — the derivative of an expression w.r.t. a symbol it does not contain is always 0. The assertion cannot fail no matter what the physics is. This is exactly the variable-independence self-test trap the auditing process flags as vacuous: independence "verified" by differentiating w.r.t. a variable the expression never contains.

The paper's claim is *conditional* and substantive: support-blindness holds "provided `f_U` and `chi_{UV}` are treated as orbit-side data rather than support-placement coordinates" (line 107). The physical content is precisely whether `f_U` and `chi_UV` carry hidden `(Lambda, varrho)` dependence. The script side-steps this entirely by declaring `f_U` and `chi_UV` as plain symbols with no `(Lambda, varrho)` dependence and then differentiating — guaranteeing zero.

**Why this matters:**
The support-blind/orbit-side split is one of the stage's six enumerated deliverables and the central structural claim of Sections 7 and 9 (and the complement to Stage 244). The script "passes" it without exercising any physics. If the closure were *wrong* — e.g., if `f_U` legitimately depended on a support coordinate — this check would still pass. The verification is vacuous.

**Required change:**
Make the independence check substantive by exercising the *conditional* the paper states. Construct the compiled objects with `f_U` and `chi_UV` carried as orbit-side data (functions of orbit variables only, e.g. `f_U = F(r)`, `chi_UV = X(r)` with `r` an orbit-side symbol that is explicitly independent of `Lam, varrho`), then differentiate and assert zero — and, to make the assertion non-vacuous, add a positive control: define a deliberately support-contaminated variant `f_U_bad = F(r) + Lam` (or `chi_UV_bad = X(r)*varrho`), rebuild one object with it, and assert that its `Lam`/`varrho` derivative is **nonzero**. The pair (real check passes, contaminated control fails to be zero) demonstrates the check has discriminating power. See directive for exact construction.

**Verification:**
New positive-control assertion(s) appear in Section 6 asserting `sp.diff(obj_contaminated, Lam) != 0` (a nonzero literal/expression), alongside the existing zero checks now built from orbit-carried `f_U`, `chi_UV`. Script still exits 0.

### F2 — paper_misalignment / insufficient_verification (R_target leg never derived from selected-branch identity; T^2 and eps_eta checks tautological)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:96-109`

**What's wrong:**
The paper's target-ratio compiler (eq:app-part08-stage245-rtarget-finite, lines 20-25, and eq:app-part08-stage245-finite-rtarget, lines 48-54) is *derived from* the selected-branch identity `R_target T^2 = Lambda_0 (1 - eps_eta)` (eq:app-part08-stage245-selected-identity, lines 15-18). The script never uses that identity. Instead it writes `R_ratio` directly as the answer (line 99):
```
R_ratio = sp.simplify(((1 - eps_eta_ref * sp.exp(V_sol)) / (1 - eps_eta_ref)) * sp.exp(-U_sol))
```
and then checks `R_exact_check` (line 100):
```
R_exact_check = sp.simplify(R_ratio * sp.exp(U_sol) * (1 - eps_eta_ref) / (1 - eps_eta_ref * sp.exp(V_sol)) - 1)
```
Substituting `R_ratio` into `R_exact_check` gives `1 - 1 = 0` identically — it is a round-trip through the same algebra that defines `R_ratio`, so A12 (line 109) cannot fail. Likewise A10/A11 (lines 107-108) check `T2/T2_ref - exp(U_sol)` and `eps_eta/eps_eta_ref - exp(V_sol)` where `T2` and `eps_eta` were *defined* as `T2_ref*exp(U_sol)` and `eps_eta_ref*exp(V_sol)` on lines 97-98 — tautological by construction.

The mismatch vs. the paper: the paper's deliverable (2b) is that the `R_target` compiler *follows from* the selected-branch identity. The script verifies a self-consistent rearrangement of the already-assumed answer, not the implication from the identity. This is a `target_mismatch`-flavored gap, but because it does not involve a numeric/value disagreement and the fix is purely script-side (derive `R_target` from the identity rather than asserting it), it is repairable without user resolution — I classify the actionable defect as `insufficient_verification` and note the alignment gap. No paper value is overwritten.

**Why this matters:**
The whole point of Section 1 of the card and notes §1 is that `R_target` is the *commuting product* of the two finite legs *via the selected-branch identity*. A reader trusting the SymPy audit would believe the implication `R_target T^2 = Lambda_0(1-eps_eta) ⟹ R_target/R_ref = [(1-eps_ref e^V)/(1-eps_ref)]e^{-U}` has been machine-checked. It has not; only an algebraic rearrangement of the conclusion has.

**Required change:**
Derive `R_ratio` from the selected-branch identity rather than asserting it. Introduce `Lambda_0` (positive symbol), define `R_target = Lambda_0*(1 - eps_eta)/T2` and `R_target_ref = Lambda_0*(1 - eps_eta_ref)/T2_ref` from the identity `R_target*T^2 = Lambda_0(1-eps_eta)`, then form `R_target/R_target_ref` and assert it `simplify`s to `((1 - eps_eta_ref*exp(V_sol))/(1 - eps_eta_ref))*exp(-U_sol)` (the paper's boxed form). This makes the check a genuine implication from the identity. Keep `T2`, `eps_eta` definitions but ensure the `R_target` derivation does not reuse `R_ratio` as an input. See directive.

**Verification:**
New assertion `assert sp.simplify(R_target_over_ref - paper_boxed_form) == 0` where `R_target_over_ref` is built from `Lambda_0*(1-eps_eta)/T2` divided by its ref, and `paper_boxed_form = ((1 - eps_eta_ref*exp(V_sol))/(1 - eps_eta_ref))*exp(-U_sol)`. Script exits 0.

### F3 — tautological_check (Session-I "rebuild" is an inverse round-trip)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:195-218`

**What's wrong:**
Section 7 infers `chi_UV_obs = aV_obs*V_obs/U_obs` (line 195, i.e. invert `V/U = chi_UV/a_V`) and `fU_obs = U_obs*Delta_obs/aV_obs` (line 198, i.e. invert `U = a_V f_U/Delta`). It then "rebuilds" `U_rebuilt = aV_obs*fU_obs/Delta_obs` (line 200) and `V_rebuilt = chi_UV_obs*fU_obs/Delta_obs` (line 201) and asserts they equal `U_obs`, `V_obs` (lines 217-218). But `U_rebuilt = a_V*(U_obs*Delta/a_V)/Delta = U_obs` exactly, and `V_rebuilt = (a_V V_obs/U_obs)*(U_obs Delta/a_V)/Delta = V_obs` exactly — both are round-trips through the inverses just applied. A20/A21 cannot fail; they test the algebra `x = g(f^{-1}(x))`, not the physics. Meanwhile the genuine session anchors `R_target/R_ref ~ 0.87984149` (paper line 522) and `R_1^nr ~ -0.12762119` (paper line 548) are computed and printed (lines 205, 214) but **never asserted**.

**Why this matters:**
The readback is meant to confirm the compiler reproduces the independently-recorded Session-I numbers. Two of the three asserted readback checks are vacuous; the two non-trivial session numbers (`R_target/R_ref`, `R_1`) are print-only. The readback gives a false sense of cross-validation.

**Required change:**
Replace the two tautological `_rebuilt` asserts with assertions against the session-recorded values: `assert abs(float(R_ratio_obs) - 0.87984149) < 5e-9` and `assert abs(float(R1_obs) - (-0.12762119)) < 5e-9` (the session numbers quoted in notes §8 / paper lines 522, 548). The `eps_eta_obs` assert (A19) already does this correctly and should stay. Do not remove the inference lines (they document the implied `chi_UV`, `f_U`); just stop asserting the self-evident round-trips.

**Verification:**
Two new numeric asserts against `0.87984149` and `-0.12762119` appear; the `< 1e-12` round-trip asserts on `U_rebuilt`/`V_rebuilt` are removed or downgraded to prints. Script exits 0.

### F4 — insufficient_verification (drain nonnegativity claim not asserted)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py:81-88`

**What's wrong:**
The paper's deliverable (3) is explicitly that the drain is *nonnegative* — eq:app-part08-stage245... and card line 73 ("it is nonnegative even when U and V have opposite signs"); notes §5 boxes `D_UV = ... >= 0`. The script verifies only the *form* `D_UV == chi_UV^2 a_V f_U^2/Delta_UV^2` (A9, line 88) and never asserts `>= 0`. The closed form makes nonnegativity manifest given `a_V > 0` (declared) and `Delta_UV^2 >= 0` and `chi_UV^2, f_U^2 >= 0`, but the load-bearing "even when chi_UV < 0 / V opposite sign to U" claim is not exercised at all.

**Why this matters:**
"Positive drain even with opposite-sign legs" is the exact-version translation of the Session-I energy-drain statement (notes §5, card line 73). A reader relying on the audit would assume the sign claim was checked; it was not.

**Required change:**
Add an assertion that the drain closed form is nonnegative. Since `D_expected = chi_UV**2 * a_V * f_U**2 / Delta_UV**2`, assert `sp.simplify(D_expected) == sp.Abs(chi_UV)**2 * a_V * f_U**2 / Delta_UV**2` is insufficient; instead use the sign structure: assert `sp.ask(sp.Q.nonnegative(D_expected))` is True under the declared assumptions, OR (more robust) assert `sp.simplify(D_expected - chi_UV**2 * a_V * f_U**2 / Delta_UV**2) == 0` AND add an explicit symbolic check that the numerator `chi_UV**2 * a_V * f_U**2 >= 0` and `Delta_UV**2 >= 0` by structure — concretely, sample the opposite-sign branch: with `chi_UV` negative, substitute a concrete admissible point (e.g. `a_U=2.5, a_V=3.0, chi_UV=-0.76, f_U=0.33`) and assert `D_expected.subs(...) > 0`, mirroring the Session-I branch where `chi_UV < 0`. See directive.

**Verification:**
A new assert demonstrating `D_UV > 0` on a concrete opposite-sign (chi_UV < 0) admissible point appears after line 88. Script exits 0.

### F5 — missing_verification_script (missing_mathematica)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 245 is checkpoint `False` and `is_status_only_candidate: False`; it is a substantive derivation. Every claim is independently checkable in Mathematica: the linear stationarity solve (`Solve`), the Hessian determinant (`Det`), the drain closed form, the `exp` finite compiler, the selected-branch implication, the first-order series (`Series`/`Normal`/`Coefficient`), the support-blindness, and the numeric Session-I readback (`N`). There is no obstruction (no special functions, no SymPy-only routine) — Mathematica can fully and independently verify this stage. Per the dual-engine rule (test is "is it possible," not "is it necessary"), a `.wl` is required and is absent.

**Why this matters:**
The stage currently has zero cross-engine confirmation. Given that the SymPy script has a vacuous support-blindness check (F1), a self-referential `R_target` check (F2), and tautological readback rebuilds (F3), a genuinely independent second engine is exactly what would catch these. The `.wl` must be an independent re-derivation, not a transliteration — in particular it should derive `R_target` *from* the selected-branch identity (the route SymPy skips) and build the support-blindness check with a real positive control.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.wl` using native Mathematica primitives and a different decomposition than the `.py`. See directive claim manifest (M1–M8).

**Verification:**
The verifier runs `math -script` on the new `.wl`; it must exit 0 with `expectZero`/`expectTrue`/`expectApprox` covering M1–M8.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed. The directive prescribes a NEW `.wl` and explicitly requires an independent decomposition (selected-branch-identity route for `R_target`; `Series` rather than a custom `coeff_linear`; positive-control support-blindness; structural nonnegativity for the drain) to avoid producing a line-by-line port of the `.py`.

## Engine cross-check

N/A — only one engine present.

## Verdict justification

`findings` (not clean, not stop_cold). The stage's core stationary-packet algebra (Section 1), dependent correction (Section 4), and first-order packet (Section 5) hold up under attack and faithfully exercise the paper claims. But four script-side defects keep it from clean: (F1) the support-blindness check is the classic vacuous variable-independence self-test — differentiating expressions w.r.t. `Lam`/`varrho` they never contain, guaranteeing zero and verifying no physics; (F2) the `R_target` compiler is asserted by definition and round-tripped rather than derived from the selected-branch identity the paper says it follows from, and the `T^2`/`eps_eta` checks are tautological by construction; (F3) the Session-I readback "rebuilds" `(U,V)` through the very inverses just applied (tautological) while leaving the two genuinely independent session numbers print-only; (F4) the drain's stated nonnegativity is never asserted. No `paper_misalignment` requiring user resolution: there are no numeric/value disagreements between paper and script — every quoted constant matches (`eps_eta ~ 0.28933482`, `R_target/R_ref ~ 0.87984149`, `R_1 ~ -0.12762119` all agree). All defects are repairable script-side without changing any paper value, so `needs_user_resolution: false` and `stop_cold: null`. Separately, the stage is single-engine and a Mathematica route is plainly possible (F5, `missing_mathematica`). Output `.txt` is fresh (mtime 1778525533 > script 1778522332).

## Self-test notes

I checked the variable-independence trap (F1: confirmed §6 differentiates objects w.r.t. `Lam`/`varrho` that never appear in them → identically zero → vacuous; this is the load-bearing finding). I checked the round-trip/tautology trap (F2: `R_exact_check` reduces to `1-1=0`; A10/A11 compare definitions to themselves; F3: `U_rebuilt`/`V_rebuilt` are `x=g(g^{-1}(x))`). For the proposed F1 positive control I confirmed `D[F(r)+Lam, Lam] = 1 != 0`, so the contaminated control genuinely fails the zero test. For F4 I confirmed `chi_UV^2 a_V f_U^2/Delta_UV^2` is `>=0` for `a_V>0` and is strictly `>0` at the opposite-sign sample point (`chi_UV<0`, `f_U!=0`). No value mismatches found between paper and script, so no user-gated `paper_misalignment`.
