---
unit_id: 117
batch: IV.3
auditor_model: claude-opus-4-7
audit_date: 2026-05-27
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage117_outlet_core_status.md"]
  paper_appendix: present
---

# Audit unit 117 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_117.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage117_outlet_core_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 1268: `\input{stages/stage_117}` between stage_116 and stage_118)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage117_outlet_core_status_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_117.tex`) has no explicit `\stagefield{Output}` equation; it carries `\claimstatus{\StatusExactClosure{} / \StatusOpen{}}`, identifies this as the "core outlet realization ledger step", and lists three diagnostic Checks (Schur complement signs in the two-channel core model; D/N mixed-tube length normalization against `L/a`; parent overlap ratios). The notes file `moving_throat_pde_stage117_outlet_core_status.md` (titled "Stage 219" internally, presumably a stale heading) is the load-bearing source: it states the concrete core preserves the canonical outgoing branch iff `g_s^2(K_s K_q + lam^2) = 4(K_s g_q - lam g_s)^2`, `L_W = (pi a/2) sqrt((1+r_c)/3)` with `r_c = lam^2/(K_s K_q)`, and that on this surface the effective outlet reduces to `Lambda_eff(z) = Lambda_2^out(z) + 4 sigma_* - sigma_*/(1 - z^2/3 - i z^5/9)` with `sigma_* = g_s^2/(4 K_s)`, preserving the canonical low-frequency expansion of `Y_2(z)`. The remaining open question is microscopic: does the actual moving-throat core realize this concrete surface?

Note: the paper card's `\section[Stage~134]` heading conflicts with `\label{stage:117}` — a numbering display quirk; not a math-side finding.

## What the script claims to verify

The script (mirrored across SymPy and Mathematica) executes a six-section classification of outlet deformations: (1) the scale/argument class has only `beta = ±1` solutions and the positive branch preserves `chi_Q = 1`; (2) pure Robin loading preserves the canonical branch only at `rho = 0`; (3) a standalone mixed pole forces a formal `kappa = -1/9` and disappearing `sigma`; (4) the hybrid outlet has two even-matching branches (cancellation `kappa=0, rho=sigma` and compensated `kappa=1/3, rho=4 sigma`), with the compensated branch at `gamma = 1/9` collapsing to `(1-sigma) Lambda_out`; (5) the concrete two-channel core obeys the coupling-balance surface `rho_c = 4 sigma_c` with `sigma_c|_balance = sigma_* = g_s^2/(4 K_s)`, and via choosing `L_W` and the bare mixed normalization realizes `kappa_c = 1/3, gamma_c = 1/9`, making `delta Lambda_core` match the notes' compensated form; (6) a capstone classification confirms only the compensated Robin-mixed core realization is nontrivial.

## Paper ↔ script cross-check

| Paper / notes deliverable | Script-side check | Status |
|---|---|---|
| Coupling-balance: `g_s^2(K_s K_q+lam^2)=4(K_s g_q-lam g_s)^2` (notes lines 33-37) | Section 5: `rho_c - 4 sigma_c = 0` solved for `g_q`; equivalent after expansion | match |
| `L_W = (pi a/2) sqrt((1+r_c)/3)` (notes lines 36-39) | Section 5: derives `Lw_required` then re-substitutes; assertion `kappa_c - 1/3 = 0` is tautological by construction | mismatch (assertion shape) |
| `Lambda_eff(z) = Lambda_2^out + 4 sigma_* - sigma_*/(1-z^2/3 - i z^5/9)` (notes lines 42-51) | Section 5 final: `delta_core - delta_core_expected = 0` matches notes formula | match |
| `sigma_* = g_s^2/(4 K_s)` (notes line 51) | Section 5: asserted equal to `sigma_c` on balance branches | match |
| Schur-complement signs (paper Checks item 1) | Implicit via squared form `(K_s g_q - lam g_s)^2`; no explicit positivity test | partial |
| D/N mixed-tube length normalization vs L/a (paper Checks item 2) | Section 5 covers it, but tautologically (see above) | partial |
| Parent overlap ratios (paper Checks item 3) | Not tested in script | missing |
| Status `\StatusExactClosure / \StatusOpen` (capstone) | Section 6 capstone classification; but the survivability flags are hardcoded into the row table rather than computed from sections 1-5 | mismatch (assertion shape) |

Paper alignment: **partial** — the math content the notes anchor IS exercised by the script in most places, but several of the load-bearing "verifications" are circular substitutions, and the capstone "only the compensated branch survives" claim is asserted against a hardcoded table rather than derived. Two paper Checks items (Schur signs, parent overlaps) are not directly tested.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 64-65 | `if {sol[beta] ...} != {-1, 1}: raise` | scale class structure | yes |
| A2 | sympy | 66 | `expect_zero(... chi_arg.subs(beta,1) - 1)` | scale class odd norm | yes |
| A3 | sympy | 80-81 | `if robin_even_solutions != [{rho:0}]: raise` | pure Robin uniqueness | yes |
| A4 | sympy | 82 | `expect_zero(... chi_R.subs(rho,0) - 1)` | pure Robin odd norm | yes |
| A5 | sympy | 98 | `expect_zero("formal even-match forces kappa=-1/9", kappa_match + 1/9)` | standalone mixed pole structure | yes (non-trivial via algebra) |
| A6 | sympy | 99 | `expect_zero("standalone mixed pole disappears", sigma_match)` | standalone mixed pole vanishes | yes |
| A7 | sympy | 100 | `expect_zero("odd norm is then trivial", chi_mix.subs(sigma,0) - 1)` | trivial at sigma=0 | yes |
| A8 | sympy | 122 | `expect_zero("hybrid cancellation branch odd norm", chi_cancel - (1 - 9 sigma gamma))` | hybrid cancellation branch odd-norm shape | yes |
| A9 | sympy | 123-126 | `expect_zero("hybrid cancellation trivial when gamma=0", ...)` | cancellation branch is trivial at gamma=0 | yes |
| A10 | sympy | 127 | `expect_zero("compensated branch odd norm", chi_comp.subs(gamma,1/9) - 1)` | compensated branch odd norm | yes |
| A11 | sympy | 128-138 | `expect_zero("compensated branch collapses to pure scale", ...)` | compensated branch == (1-sigma) Lambda_out | yes |
| A12 | sympy | 149-152 | `expect_zero("both core-balance branches give same sigma_*", ...)` | both g_q roots → same sigma | yes |
| A13 | sympy | 153 | `expect_zero("core-balance sigma_* value", sigma_c - sigma_*)` | sigma_c on balance = g_s^2/(4 K_s) | yes |
| A14 | sympy | 156 | `expect_zero("D/N tube fixes kappa_c=1/3", kappa_c.subs(kappa0, 4 Lw_required^2/(pi^2 a^2)) - 1/3)` | D/N normalization → kappa_c=1/3 | **no — tautological** (Lw_required is solved from exactly that equation) |
| A15 | sympy | 157 | `expect_zero("bare mixed fixes gamma_c=1/9", gamma_c.subs(gamma0, (1+r_c)/9) - 1/9)` | gamma_c normalization | **no — tautological** (gamma_c = gamma0/(1+r_c); substitute gamma0=(1+r_c)/9) |
| A16 | sympy | 167 | `expect_zero("concrete core collapses to compensated hybrid class", ...)` | final Lambda_eff form match notes | yes |
| A17 | sympy | 181-182 | `if nontrivial_survivors != ["compensated Robin-mixed core realization"]: raise` | capstone "only compensated survives" | **no — hardcoded table** |
| Bn  | mathematica | 38-137 | mirrors A1-A17 1:1 | same set | mirror |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:34-137`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:50-182`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script, not an independent re-derivation. Compare:

SymPy (`.py:54`):
```
Y_arg = sp.expand(sp.series((-3 * S) / (S * Lambda_out.subs(z, beta * z)), z, 0, 6).removeO())
```
Mathematica (`.wl:37`):
```
yArg = Normal[Series[(-3 s)/(s (lambdaOut /. z -> beta z)), {z, 0, 5}]];
```

SymPy (`.py:103`):
```
Lambda_hyb = sp.expand(sp.series(Lambda_out + rho - sigma / (1 - kappa * z**2 - I * gamma * z**5), z, 0, 6).removeO())
```
Mathematica (`.wl:75`):
```
lambdaHyb = Normal[Series[lambdaOut + rho - sigma/(1 - kappa z^2 - I gamma z^5), {z, 0, 5}]];
```

SymPy (`.py:141-145`):
```
r_c = sp.simplify(lam**2 / (Ks * Kq))
rho_c = sp.simplify(gs**2 / Ks)
sigma_c = sp.simplify((Ks * gq - lam * gs) ** 2 / (Ks**2 * Kq * (1 + r_c)))
kappa_c = sp.simplify(kappa0 / (1 + r_c))
gamma_c = sp.simplify(gamma0 / (1 + r_c))
```
Mathematica (`.wl:98-102`):
```
rC = FullSimplify[lam^2/(ks kq)];
rhoC = FullSimplify[gs^2/ks];
sigmaC = FullSimplify[(ks gq - lam gs)^2/(ks^2 kq (1 + rC))];
kappaC = FullSimplify[kappa0/(1 + rC)];
gammaC = FullSimplify[gamma0/(1 + rC)];
```

The two scripts share section banners, intermediate variable names (with naming convention conversions), check labels (the verbatim strings inside `expect_zero`/`expectZero`), and the entire choreography of sections 1-6. The Mathematica script does not derive `r_c`, `sigma_c`, `Lambda_eff`, or the balance equation from independent premises — it re-uses the SymPy script's pre-computed algebraic forms. This violates the second-engine policy that both engines must derive the result independently.

**Why this matters:**
A second engine that merely mirrors the first cannot catch transcription errors, sign conventions, or algebraic missteps in the first; it only catches Mathematica-vs-SymPy CAS-implementation disagreement, which is not the intent. The cross-check value here is near zero.

**Required change:**
Rewrite the Mathematica script so it derives the concrete-core balance equation, the D/N tube normalization, and the `Lambda_eff` formula independently. A defensible independent derivation would: (i) introduce the Schur-complement-style two-channel coefficient construction from first principles in Mathematica (rather than carrying over `sigma_c = (K_s g_q - lam g_s)^2 / (K_s^2 K_q (1+r_c))` as a typed-in formula), (ii) compute the series of `Lambda_eff(z)` and compare against `Lambda_out` independently, and (iii) verify the coupling-balance and D/N-tube relations by independent algebraic paths (e.g., solve the original 2x2 Schur system in Mathematica from the physical setup, rather than re-typing SymPy's reduced formula).

Because (a) this is structural rather than a single-line bug, and (b) re-doing the second-engine work is the auditor's job to *flag*, not Codex's job to invent physics, this finding is reported but the directive marks it `## Blocked` — the user must direct how to re-derive in Mathematica (or accept the transliteration as a known limitation).

**Verification:**
After rewrite, the Mathematica script's intermediate variable names, banner strings, and section ordering should differ from the SymPy script's; the same final identities should be reached by different intermediate steps.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:155-156`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:112-113`

**What's wrong:**
The script defines `Lw_required` as the solution of the equation `4 Lw^2/(pi^2 a^2) = (1+r_c)/3`, then asserts:
```python
Lw_required = sp.solve(sp.Eq(4 * Lw**2 / (sp.pi**2 * a**2), (1 + r_c) / 3), Lw)[0]
expect_zero("D/N tube fixes kappa_c = 1/3", kappa_c.subs(kappa0, 4 * Lw_required**2 / (sp.pi**2 * a**2)) - sp.Rational(1, 3))
```
Since `kappa_c = kappa0/(1+r_c)` and `Lw_required` is by construction the value that makes `4 Lw^2/(pi^2 a^2) = (1+r_c)/3`, the substitution `kappa0 → 4 Lw_required^2/(pi^2 a^2)` gives `kappa0 = (1+r_c)/3`, so `kappa_c = ((1+r_c)/3)/(1+r_c) = 1/3` is forced algebraically by the very equation that defined `Lw_required`. The assertion cannot fail.

Notes claim: `L_W = (pi a / 2) sqrt((1+r_c)/3)` (notes line 36-39) — this is the D/N tube normalization the paper card refers to. The script reaches the same `L_W` value but doesn't *derive* it from a D/N tube boundary condition: it solves backward from the desired `kappa_c = 1/3`.

**Why this matters:**
The "D/N mixed-tube length normalization" check is one of the three Checks the paper card explicitly enumerates. Verifying it by inverting the conclusion provides no physical content — any number `c` in place of `1/3` could be substituted at the start and the assertion would still pass.

**Required change:**
Replace the inversion-based check with a forward derivation: start from a D/N tube of length `L_W` with the mixed boundary conditions documented upstream (in the notes or upstream stages) and compute `kappa0` (or whatever bare quantity) from that boundary value problem. Then evaluate `kappa_c = kappa0/(1+r_c)` and assert `kappa_c - 1/3 = 0` using `L_W = (pi a/2) sqrt((1+r_c)/3)` substituted into the *forward* expression for `kappa0`, not into the equation that defined it.

Concrete script-side fix shape (sympy): define `kappa0_DN = f(L_W, a, r_c)` where `f` comes from the D/N tube boundary problem (not from inverting the kappa_c=1/3 requirement), then substitute `L_W = pi a sqrt((1+r_c)/3) / 2` and assert the resulting `kappa_c = 1/3`. If the forward `f` is not derived in this stage and is supposed to be carried forward from an upstream stage, then this check should be marked as a carry-forward and the docstring should cite the upstream unit; the current script does neither.

**Verification:**
The new check's expression must not contain `Lw_required` (or any quantity solved from the conclusion). It must contain a forward expression for `kappa0` in terms of `L_W` (and other physical inputs), with `L_W` then substituted to its stated value.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:157`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:114`

**What's wrong:**
```python
expect_zero("bare mixed normalization fixes gamma_c = 1/9", gamma_c.subs(gamma0, (1 + r_c) / 9) - sp.Rational(1, 9))
```
With `gamma_c = gamma0/(1+r_c)`, substituting `gamma0 → (1+r_c)/9` yields `((1+r_c)/9)/(1+r_c) = 1/9` by pure algebraic cancellation. The "bare mixed normalization" is never derived — the script just *declares* `gamma0 = (1+r_c)/9` is the normalization and then asserts the algebraic consequence.

**Why this matters:**
Same issue as F2: an assertion that cannot fail because the input was chosen to make it succeed. The "bare mixed normalization" should be derived from the physical mixed-channel setup (e.g., from how `gamma0` arises in the mixed pole's bare construction), not assumed.

**Required change:**
Either (a) replace the assertion with a forward derivation of `gamma0` from the bare mixed-channel construction (cite the source notes lines defining what `gamma0` is physically), then substitute and verify `gamma_c = 1/9`, or (b) re-frame the line as a definitional substitution rather than a "fix": rename the print label to "set bare mixed normalization gamma0 := (1+r_c)/9" and remove the misleading `expect_zero` since it's testing nothing. Option (a) is preferred if the upstream derivation is available in this stage's scope; option (b) is acceptable if the user accepts that `gamma0 = (1+r_c)/9` is a carry-forward definition.

**Verification:**
After fix, either the assertion compares a *derived* expression (containing physical inputs like channel widths, residues, etc.) against `1/9`, or the line is no longer an `expect_zero` call but a labeled definition.

### F4 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage117_outlet_core_status_sympy_audit.py:170-182`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl:126-137`

**What's wrong:**
The capstone classification table is hardcoded:
```python
classification_rows = [
    ("scale/argument", True, True, False, "harmless beta = 1 pure-scale branch"),
    ("pure Robin", False, False, False, "rho_R = 0 only"),
    ("standalone mixed pole", False, False, False, "sigma_W = 0 only (formal kappa = -1/9)"),
    ("hybrid cancellation", True, True, False, "gamma_W = 0 reduces to exact cancellation"),
    ("compensated Robin-mixed core realization", True, True, True, "balance surface + D/N tube normalization"),
]
...
nontrivial_survivors = [name for name, even_ok, odd_ok, nontrivial, _ in classification_rows if even_ok and odd_ok and nontrivial]
if nontrivial_survivors != ["compensated Robin-mixed core realization"]:
    raise AssertionError(...)
```
The `True/False` flags in column positions 2-4 are not computed from sections 1-5; they are typed in. The "uniqueness of the surviving branch" is asserted by filtering the typed-in flags. The check cannot fail unless a developer edits the table by hand.

**Why this matters:**
The capstone is supposed to be the "tag the closure" line — `\StatusExactClosure` in the paper card. If the section-1-5 sub-checks fail, the capstone should reflect that. As written, the capstone is independent of the prior section outcomes.

**Required change:**
Wire the flags to computed booleans from sections 1-5. E.g., set `even_ok_scale = (beta_solutions matches {±1})`, `odd_ok_scale = (chi_arg.subs(beta,1) == 1)`, `nontrivial_scale = False` (the scale class is trivially the canonical itself with rescaling), and similar for the other four classes. Then assert `nontrivial_survivors == ["compensated Robin-mixed core realization"]` against the *computed* flags. If the boolean for `nontrivial_compensated` is derived from "the section-5 final delta_core matches delta_core_expected" residual being zero, then the capstone assertion has real content.

**Verification:**
After fix, the `classification_rows` should be a list of named tuples whose 2nd-4th entries are expressions (or variables holding `True`/`False`) computed from earlier section results — not boolean literals. If any section's assertion fails, the corresponding row's `even_ok`/`odd_ok`/`nontrivial` flag should reflect that, so the capstone catches the regression.

## Independent-derivation check (Mathematica)

The Mathematica script is a literal transliteration of the SymPy script. Comparing corresponding sections:

- `.py:50` ↔ `.wl:34` — identical `Lambda_out = -3 + z^2/3 + z^4/9 + I z^5/9` definition.
- `.py:54` ↔ `.wl:37` — identical `Y_arg` series with the same dummy variable `S`/`s` to be canceled.
- `.py:58-62` ↔ `.wl:41` — same simultaneous-equation solve for `beta`.
- `.py:103-115` ↔ `.wl:75-80` — same hybrid series, same solve for `{rho, kappa}` with identical even-match equations.
- `.py:141-145` ↔ `.wl:98-102` — same five concrete-core symbolic definitions.
- `.py:155-156` ↔ `.wl:112-113` — same inversion to get `Lw_required` and same tautological re-substitution.

All section banners are character-identical ("STAGE 100 — CONCRETE OUTLET-CORE STATUS" — note both use the stale "STAGE 100" label, see below), and every `expect_zero` / `expectZero` label string matches verbatim. This is finding F1.

Aside: the banner says "STAGE 100" in both scripts when the unit ID is 117 (and the paper section heading is "Stage 134"). This is purely cosmetic/informational and does not constitute a math finding, but it reinforces that the scripts are derived from a template copied from another unit.

## Engine cross-check

Both engines produce the same output structure with `EXIT_CODE: 0`. All asserted residuals print as `0`. Engine outputs are consistent (engine_disagreement: n/a since the scripts are mirrors). Outputs are fresh: sympy `.txt` mtime 2026-05-11 12:45 vs `.py` mtime 11:56 (fresh); mathematica `.txt` mtime 2026-05-11 13:09 vs `.wl` mtime 11:56 (fresh).

## Verdict justification

**Verdict: findings.** The script does verify the load-bearing identities from the notes (coupling-balance equation, `sigma_* = g_s^2/(4 K_s)`, `Lambda_eff` form match), and those checks (A12, A13, A16) are substantive and non-trivial. However: (i) the Mathematica script is a transliteration, providing no independent cross-engine check (F1); (ii) two of the four "D/N tube + bare mixed normalization" checks central to the paper's Checks list are tautological substitutions of values defined to satisfy them (F2, F3); (iii) the capstone classification asserts against hardcoded boolean flags rather than computed ones (F4). The paper card is also vague (no `\stagefield{Output}` equation; only narrative + status tag), so the only authoritative anchor is the notes file — which IS faithfully verified by the substantive checks, just with the cited circular substitutions for the D/N and gamma normalizations.

Attacks attempted that did not produce findings: (a) verified `chi_comp` algebra on the compensated branch — `(9 sigma gamma - 1)/(sigma - 1)` at `gamma = 1/9` evaluates to 1 ✓; (b) verified `sigma_c|_balance = g_s^2/(4 K_s)` algebraically — both `g_q` roots give `(K_s g_q - lam g_s)^2 = g_s^2 (K_s K_q + lam^2)/4`, times `1/(K_s^2 K_q (1+r_c)) = 1/(K_s(K_s K_q + lam^2))` ✓; (c) verified the coupling-balance form is equivalent to the notes' equation after multiplying through ✓; (d) checked symbol assumption domains — `Ks, Kq, lam, kappa0, gamma0, a, Lw` all positive, `rho, sigma, kappa, gamma, gs, gq` real-only; consistent with the physical setup; (e) checked that the standalone-mixed-pole sigma_match = 0 algebra is correct (it is, after expanding `(1/3 + sigma/9)^2/(3+sigma)^2 - (1/9 - sigma/81)/(-3-sigma) = 4/81` yields `sigma = 0`).

`paper_alignment: partial` — the script's `Output`-equivalent claims (the notes' equations) are reached, but two of the three Checks items the paper enumerates (D/N tube length normalization, parent overlap ratios) are not properly exercised: D/N is asserted tautologically, parent overlaps are not tested at all. Schur sign checks are absorbed implicitly via the squared form.

No `stop_cold` — the findings are bounded, do not invalidate the substantive checks A12/A13/A16, and do not propagate downstream in a way that would alter results carried forward. The compensated-branch identification, `sigma_* = g_s^2/(4 K_s)`, and the `Lambda_eff` form are all real conclusions of the substantive checks; only their D/N and gamma normalizations are wrapped in circular substitutions, and the values themselves (`kappa_c = 1/3`, `gamma_c = 1/9`) are confirmed elsewhere by the section-5 final `delta_core - delta_core_expected = 0` check (A16), which is substantive.

## Self-test notes

- **Variable independence**: No new `sp.diff` or `D[expr, var]` is proposed; F2 and F3 fixes are about restructuring existing algebraic chains, not introducing derivatives.
- **Symmetry/parity**: No new integrals proposed.
- **Trivial-case pre-check**: For F2 — a forward `kappa0 = f(L_W, a, ...)` substituted with `L_W = pi a sqrt((1+r_c)/3)/2` must give `kappa0 = (1+r_c)/3` (not trivially because of the function `f`), then `kappa_c = 1/3` follows non-trivially. Verified mentally. For F3 — analogous. For F4 — the proposed wiring uses booleans computed from existing residuals (e.g., `nontrivial_compensated = (delta_core - delta_core_expected).simplify() == 0`); these were already verified zero in the existing run.
- **Path specifications**: F1, F2, F3, F4 all target existing files; paths verified.
- **Paper round-trip**: The fixes do not introduce new constants; they restructure existing checks. The values `1/3`, `1/9`, `g_s^2/(4 K_s)` remain unchanged and match the notes.
- **F1 is blocked, not auto-applied**: Re-deriving the second engine from independent premises requires physics judgment about the Schur-complement construction and the D/N tube BVP that Codex cannot safely invent. The directive marks F1 with `## Blocked` requesting user guidance.
