---
unit_id: 068
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage068_resonance_thresholds.md]
  paper_appendix: present
---

# Audit unit 068 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_068.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage068_resonance_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 114; `\input{stages/stage_068}` at line 254)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.txt`

## What the paper claims

The paper card states the stage delivers a single load-bearing scalar: the threshold correction factor `P_res = 1.005612487760576` that multiplies the matched-branch wall figure of merit (carried in from stage 066) when the independent sech-Gaussian profile family is used in place of the matched profile. The body text quotes `P_res = 1.005612...` verbatim and asserts that the independent benchmark raises the required wall figure by about 0.56%. The `\stagefield{Output}` line reads "Threshold correction factor `P_res`." The notes elaborate three deliverables: (i) the resonance-corrected wall figure `W_res(r) = C^2(r) W_wall`; (ii) the threshold translation `W_wall <= Pe_req/[C^2 Delta_inf]` (fail) and `W_wall >= Pe_req/[C^2 Delta_0]` (succeed), specialising at the resonance to `P_res * Pe_req / Delta_{inf,0}`; (iii) the profile-sensitivity band of width `(P_res - 1) * Pe_req / Delta_{0,inf}` of about 0.56% on either side. `P_res` itself is defined as `1/C_res^2`, with `C_res^2 = 0.994418...` inherited from stage 067.

## What the script claims to verify

Both scripts (SymPy and the Mathematica mirror) claim to (1) derive `W_res = C^2 W_wall` from "amplitude-squared = power" reasoning; (2) derive `P_res = 1/C_res^2` as the inverse amplification; (3) compute matched-branch thresholds `Wfail_match = Pe_req/Delta_inf` and `Wsuff_match = Pe_req/Delta_0` from `W_match * Delta = Pe_req` via Solve; (4) compute profile-family thresholds `Wfail_res = Pe_req/(C2*Delta_inf)` and `Wsuff_res = Pe_req/(C2*Delta_0)` from `C2 * W_prof * Delta = Pe_req` via Solve; (5) check the algebraic cross-relations `W_res * C2 - W_match = 0` and the resonance substitution `C2 -> 1/P_res` collapsing to `P_res * W_match`; (6) compute the band widths two algebraic ways and assert agreement. No numeric anchor for `P_res ~ 1.005612...` is exercised in-script; that is carried symbolically as `1/C_res^2`.

## Paper-script cross-check

| Paper-side deliverable | Script-side coverage | Status |
|---|---|---|
| `W_res(r) = C^2(r) W_wall` (notes section 1) | SymPy line 60 / Mathematica line 38 - but constructed by relabel-then-compare (see F1) | match (tautological) |
| `P_res = 1/C_res^2` (notes section 2) | SymPy line 66 / Mathematica line 43 - Solve produces exact RHS, then compared to RHS (see F1) | match (tautological) |
| Profile-family thresholds `Pe_req/[C^2 Delta_{0,inf}]` (notes section 2) | SymPy lines 93-94 / Mathematica lines 52-53 via independent Solve on `C2*W_prof*Delta = Pe_req` | match |
| Matched-branch thresholds `Pe_req/Delta_{0,inf}` (notes section 2) | SymPy lines 85-86 / Mathematica lines 47-48 via Solve on `W_match*Delta = Pe_req` | match |
| Resonance specialisation `W^(res,*) = P_res * W_match` | SymPy lines 107-110 / Mathematica lines 65-66 | match |
| Band width `(P_res-1) Pe_req / Delta_{0,inf}` (notes section 4) | SymPy lines 118-125 / Mathematica lines 70-76 (two algebraic paths, both end at same expression) | partial (see F3) |
| Numeric value `P_res ~ 1.005612487760576` (card eq. and notes section 2) | (not exercised) - `P_res` and `C_res` left fully symbolic; numeric value inherited from stage 067 | extra-to-script / acceptable carry-forward |
| `\stagefield{Output}` = `P_res` factor | Combined coverage of (i)-(iii) above implies the factor is correctly defined and applied | match |

Paper alignment is `aligned` in substance: the load-bearing equations exported by the script match the paper's load-bearing equations exactly. The v1 material change (lifting `Wfail_res / Wfail_match` from postulated literals to Solve-derived expressions on the matched and profile Peclet balances) is intact: both pairs are derived via `sp.solve` / `Solve`, and the resulting expressions agree term-by-term with the notes section 2 forms.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 60 | `W_res_derived - C2 * W_wall == 0` | (i) `W_res = C^2 W_wall` | no (tautological - see F1) |
| A2 | sympy | 66 | `Pres_derived - 1/Cres**2 == 0` | (ii) `P_res = 1/C_res^2` | no (tautological - see F1) |
| A3 | sympy | 68 | `(1/Cres**2)*Cres**2 - 1 == 0` | "consistency of `P_res*C_res^2 = 1`" | no (P_res not even in the LHS - see F1) |
| A4 | sympy | 103 | `Wfail_res * C2 - Wfail_match == 0` | thresholds scale by `1/C^2` | yes - links the two independent Solves |
| A5 | sympy | 104 | `Wsuff_res * C2 - Wsuff_match == 0` | thresholds scale by `1/C^2` | yes |
| A6 | sympy | 107-108 | `Wfail_res(C2->1/Pres) - Pres*Wfail_match == 0` | resonance-point specialisation | yes (light - algebraic) |
| A7 | sympy | 109-110 | `Wsuff_res(C2->1/Pres) - Pres*Wsuff_match == 0` | resonance-point specialisation | yes (light - algebraic) |
| A8 | sympy | 133 | success band Way A = Way B | band width identity | partial (Way B doesn't use `Wsuff_res` - see F3) |
| A9 | sympy | 134 | failure band Way A = Way B | band width identity | partial (see F3) |
| A10 | mathematica | 38 | `WresDerived - C2*Wwall == 0` | (i) `W_res = C^2 W_wall` | no (tautological - see F1) |
| A11 | mathematica | 43 | `PresDerived - 1/Cres^2 == 0` | (ii) `P_res = 1/C_res^2` | no (tautological - see F1) |
| A12 | mathematica | 61-62 | `WfailRes*C2 - WfailMatch == 0` etc. | thresholds scale by `1/C^2` | yes |
| A13 | mathematica | 65-66 | resonance substitution checks | resonance-point specialisation | yes (light) |
| A14 | mathematica | 83-84 | band widths Way A = Way B | band width identity | partial (see F3) |

## Findings

### F1 - tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:53-68`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:36-43`

**What's wrong:**
Three load-bearing "derivations" near the top of both scripts are algebraic identities by construction, not derivations.

(a) `W_res = C^2 W_wall`. SymPy lines 53-60:
```
A_in, C_sym = sp.symbols("A_in C", positive=True, real=True)
A_trans = C_sym * A_in
W_wall_def = A_in**2
W_res_derived = sp.simplify((A_trans**2).subs(A_in**2, W_wall) * 1)
W_res_derived = sp.simplify(W_res_derived.subs(C_sym**2, C2))
expect_zero("W_res - C2 * W_wall", W_res_derived - C2 * W_wall)
```
The construction is: `A_trans**2 = C**2 * A_in**2`, then literal substitutions `A_in**2 -> W_wall` and `C**2 -> C2` produce `C2 * W_wall` by relabeling. Subtracting `C2 * W_wall` cannot fail; nothing about "power = amplitude squared" is verified because the only step that uses that premise is the line `W_wall_def = A_in**2`, which is itself an unverified definition.

Mathematica lines 36-38 are worse:
```
WresRule = First@Solve[Wres == (C2)*Wwall, Wres];
WresDerived = Wres /. WresRule;
expectZero["W_res - C2 * W_wall", WresDerived - C2*Wwall];
```
Solve `Wres == C2*Wwall` for `Wres` returns `Wres -> C2*Wwall`. Comparing to `C2*Wwall` is `C2*Wwall - C2*Wwall = 0` by definition.

(b) `P_res = 1/C_res^2`. SymPy lines 64-66:
```
Cres = sp.symbols("C_res", positive=True, real=True)
Pres_derived = sp.simplify(1 / Cres**2)
expect_zero("P_res - 1/C_res^2", Pres_derived - 1/Cres**2)
```
`(1/Cres**2) - (1/Cres**2) = 0` by construction. Mathematica lines 41-43 use Solve to extract the same exact form: tautological.

(c) `P_res*C_res^2 - 1 = 0`. SymPy line 68:
```
expect_zero("P_res*C_res^2 - 1", (1/Cres**2)*Cres**2 - 1)
```
The LHS does not contain `P_res` at all - only `Cres`. It is literally `1 - 1 = 0`. The label is misleading; the symbol `Pres` declared at line 42 is never substituted in.

**Why this matters:**
Three of the script's stated bottom-line claims (i)-(ii) and the consistency check (c) cannot fail under any value of the physics. The notes locate the substantive physical content in section 1: "Power = amplitude^2" is the load-bearing premise that connects `C` (transmission coefficient) to `C^2` (the power-level coherence factor used elsewhere in the ledger). The script writes this premise as a definition (`W_wall_def = A_in**2`) and never tests it against any other characterisation of `W_wall` - there is no `W_wall` carried in from upstream that the script could equate against `A_in^2`. So the chain "transmission -> amplitude -> squared -> power" is a sequence of labels, not a verification. If `W_res` had been defined in stage 066 as `(something not equal to A_in^2)`, the script would not catch it.

**Required change:**
Replace the three tautological checks with substantive anchors:

- **(a) `W_res = C^2 W_wall`**: anchor `W_wall` in the matched-branch gain expression carried from stage 066/notes section 1: `G_match = rho_star * g_phi^2 * N_phiphi / (m * c_s^2 * K_X)` and `W_wall = kappa * G_match`. The profile-family gain is `G_res = C^2 * G_match`, so `W_res = kappa * G_res = C^2 * (kappa * G_match) = C^2 * W_wall`. Encode this as: define `Gmatch_expr` in terms of `rho_star, g_phi, N_phiphi, m, c_s, K_X`; define `Wwall_expr = kappa * Gmatch_expr`; define `Gres_expr = C2 * Gmatch_expr`; define `Wres_expr = kappa * Gres_expr`; assert `simplify(Wres_expr - C2 * Wwall_expr) == 0`. This is substantive because it would fail if any of the four substitutions (`G_match`, `kappa`, `G_res`, `W_res`) were perturbed independently.

- **(b) `P_res = 1/C_res^2`**: anchor `P_res` in the gain-amplification interpretation from notes section 1. The threshold balance at the matched branch is `W_wall * Delta = Pe_req`; replacing `W_wall` with `W_res = C^2 W_wall` gives `C^2 * W_wall * Delta = Pe_req`, i.e. the required wall figure is now `Pe_req/(C^2 Delta)`. The amplification of the required wall figure from matched to profile-family is then `1/C^2`. Define `Pres_expr` as the ratio `Wfail_res / Wfail_match` (or `Wsuff_res / Wsuff_match`) at the resonance `C2 -> Cres^2`, and assert it equals `1/Cres^2`. Both sides will now flow from `sp.solve` on the matched and profile Peclet balances, so the assertion exercises both Solves and the substitution.

- **(c) Replace line 68 entirely** with a numeric anchor: set `Cres_sq_numeric = sp.Float("0.994418836451529", 20)` and `Pres_numeric = 1 / Cres_sq_numeric`, then assert `abs(Pres_numeric - sp.Float("1.005612487760576", 20)) < sp.Float("1e-12", 20)`. The paper card states the numeric value to 15 digits; the audit must check it lines up with stage 067's output. (If the carry-forward value lives in a constants file, import from there; otherwise the paper card's quoted constants are the in-script anchor.)

Apply the analogous changes to the Mathematica script (`mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`), but write them as an independent derivation route (see F2).

**Verification:**
After the fix, `Wres_expr` is built from `kappa, rho_star, g_phi, N_phiphi, m, c_s, K_X` and the assertion `Wres_expr - C2 * Wwall_expr == 0` non-trivially collapses these symbols; perturbing any single one in the script should cause the assertion to fail. The numeric `P_res` check should produce a residual below 1e-12 in the printed output.

### F2 - mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:36-84`

**What's wrong:**
The `.wl` script is structurally a line-by-line port of the `.py` script, not an independent re-derivation. Three corresponding pairs:

SymPy lines 84-94 vs Mathematica lines 46-53 - same Solve choreography:
```
# SymPy:
W_match_sol = sp.solve(sp.Eq(W_match_sym * Delta_sym, Pe_req), W_match_sym)[0]
Wfail_match = sp.simplify(W_match_sol.subs(Delta_sym, Deltainf))
W_prof_sol = sp.solve(sp.Eq(C2 * W_prof_sym * Delta_sym, Pe_req), W_prof_sym)[0]
Wfail_res = sp.simplify(W_prof_sol.subs(Delta_sym, Deltainf))
```
```
(* Mathematica: *)
WmatchSol = First@Solve[Wmatch*DeltaSym == PeReq, Wmatch];
WfailMatch = FullSimplify[(Wmatch /. WmatchSol) /. DeltaSym -> Deltainf];
WprofSol = First@Solve[C2*Wprof*DeltaSym == PeReq, Wprof];
WfailRes = FullSimplify[(Wprof /. WprofSol) /. DeltaSym -> Deltainf];
```
Same symbol names (capitalised), same number and order of Solve calls, same substitution sequence, same intermediate variable names.

SymPy lines 122-125 vs Mathematica lines 73-76 - same "Way B" Solve idiom:
```
gap_sym = sp.symbols("gap", real=True)
success_band_widthB = sp.solve(sp.Eq(Wsuff_match + gap_sym, Pres * Wsuff_match), gap_sym)[0]
```
```
gapSym;
successBandB = gap /. First@Solve[WsuffMatch + gap == Pres*WsuffMatch, gap];
```
The leftover statement `gapSym;` on line 74 of the `.wl` does nothing useful - it's the literal port of `gap_sym = sp.symbols("gap", real=True)` with no Mathematica equivalent function (Mathematica does not need symbol declarations), so it was kept as a no-op. This is a strong tell that the script was generated by translation, not co-design.

SymPy lines 56-59 (relabel chain through `A_trans = C * A_in`) vs Mathematica lines 36-38 (single-line Solve): the Mathematica side is even shorter, suggesting the author saw that the "derivation" was vacuous and compressed it, but did not replace it with an independent derivation.

**Why this matters:**
The dual-engine policy exists so that two independent symbolic engines must each derive the result from the physical premises. If the Mathematica script's algebra mirrors the SymPy choreography, any single error in the algebraic approach (e.g., a wrong sign in the premise, a swapped numerator/denominator, an inverted ratio) will be present in both - agreement between engines no longer constitutes verification.

**Required change:**
Rewrite the Mathematica script so it derives the same final claims from the matched-branch physical premises along a different algebraic path. Concrete suggestions:

- Define the matched threshold and profile threshold not by Solve on a Peclet balance, but by `Reduce[]` on the explicit boundary-condition form (e.g., set up `Reduce[W*C2*Delta == PeReq && W > 0, W]` for each `Delta`).
- Derive the band width by `Series[]` expansion of `Pres` around 1: `Series[Pres*WsuffMatch, {Pres, 1, 1}]` and read off the first-order coefficient as the band width, plus a `FullSimplify` cross-check against the closed form.
- Verify `Wres = C^2 Wwall` from the matched-gain side (see F1 reform), using `FullSimplify[...]` over the four-symbol gain `rho_star * g_phi^2 * N_phiphi / (m * c_s^2 * K_X)` rather than the amplitude-squared relabel.

The two engines must reach `Pe_req/Delta_inf`, `Pe_req/Delta_0`, `Pe_req/(C2 Delta_inf)`, `Pe_req/(C2 Delta_0)`, and band-width `Pe_req (Pres - 1)/Delta_{0,inf}` along non-parallel symbolic paths. The leftover `gapSym;` line should be removed regardless.

**Verification:**
The verifier compares the final printed expressions from both engines: same values are required, but the intermediate `Print` statements (which the script already emits) should now show different algebraic intermediates between `.py` and `.wl`. A spot check: line ~75 of the new `.wl` should no longer be a Solve on `WsuffMatch + gap == Pres*WsuffMatch`; that's the SymPy choreography.

### F3 - insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:118-134`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:70-84`

**What's wrong:**
The "band widths computed two independent ways" check is not actually two independent computations.

Way A (line 118):
```
success_band_widthA = sp.simplify((Wsuff_res - Wsuff_match).subs(C2, 1 / Pres))
```
Substituting `C2 -> 1/Pres` in `Wsuff_res = Pe_req/(C2*Delta_0)` gives `Pres * Pe_req / Delta_0 = Pres * Wsuff_match`. So Way A is `Pres * Wsuff_match - Wsuff_match = (Pres - 1) * Wsuff_match`.

Way B (line 124):
```
success_band_widthB = sp.solve(sp.Eq(Wsuff_match + gap_sym, Pres * Wsuff_match), gap_sym)[0]
```
Solve `Wsuff_match + gap = Pres * Wsuff_match` for `gap`: `gap = (Pres - 1) * Wsuff_match`.

Both ways compute `(Pres - 1) * Wsuff_match` - Way B never references `Wsuff_res`. It is an algebraic restatement of the same identity using a different syntactic path. If hypothetically `Wsuff_res` were defined incorrectly (e.g., `Pe_req/(Pres*Delta_0)` accidentally), Way A would catch it but Way B would not - and the comparison `A - B = 0` would fail. Conversely, if `Pres` were defined incorrectly as something other than `1/Cres^2`, both Ways would absorb the error identically and the comparison would pass.

So the check verifies that `(Wsuff_res at C2=1/Pres) = Pres * Wsuff_match`, which is exactly the same identity already checked at lines 109-110. It is redundant rather than reinforcing.

**Why this matters:**
The script claims (docstring item 4) "Exact width of the profile-sensitive threshold band" is checked. The current check only confirms what was already checked, not a new derivation. The paper notes section 4 give the band width as `Pe_req/Delta - Pe_req/(C^2 Delta)` evaluated at `C^2 = Cres^2`, equal to `(Pres - 1) Pe_req / Delta`. A genuine independent check would derive the band width from the paper's section 4 framing (matched success threshold vs. profile-family success threshold, both expressed in `Cres`) without going through `Pres` substitution.

**Required change:**
Replace Way B with a derivation that does not pre-bake the answer. For example:

- Define `Wsuff_res_at_resonance = Pe_req / (Cres**2 * Delta_0)` (using `Cres` directly, not `Pres`).
- Define `band_width_C_form = Wsuff_res_at_resonance - Wsuff_match`.
- Define `band_width_P_form = (Pres - 1) * Wsuff_match` using `Pres = 1/Cres**2`.
- Assert `simplify(band_width_C_form.subs(Pres, 1/Cres**2) - band_width_P_form.subs(Pres, 1/Cres**2)) == 0`.

This routes the verification through `Cres` (the upstream-anchored symbol) rather than `Pres` (the derived symbol), so the assertion fails if anyone perturbs the `Pres = 1/Cres^2` link.

Apply the analogous reformulation to the Mathematica script.

**Verification:**
After the fix, the printed band-width expressions should reference `Cres^2` somewhere in the chain (not just `Pres`), and a deliberate `Cres -> 2*Cres` perturbation in a side test should cause the check to fail.

## Independent-derivation check (Mathematica)

See F2: the Mathematica script is a line-by-line port of the SymPy script, including a leftover no-op declaration (`gapSym;` on line 74). The Solve choreography, variable choices, and assertion order are parallel. This is the `mathematica_transliteration` finding.

## Engine cross-check

The two saved transcripts agree on every printed expression at the symbol level (after accounting for naming conventions `Delta_inf` <-> `Deltainf`, `Pe_req` <-> `PeReq`, etc.):

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `W_res - C2 * W_wall` | `0` | `0` |
| `Wfail_match` | `Pe_req/Delta_inf` | `PeReq/Deltainf` |
| `Wsuff_match` | `Pe_req/Delta_0` | `PeReq/Delta0` |
| `Wfail_res` | `Pe_req/(C2*Delta_inf)` | `PeReq/(C2*Deltainf)` |
| `Wsuff_res` | `Pe_req/(C2*Delta_0)` | `PeReq/(C2*Delta0)` |
| success band A | `Pe_req*(P_res - 1)/Delta_0` | `(PeReq*(-1 + Pres))/Delta0` |
| failure band A | `Pe_req*(P_res - 1)/Delta_inf` | `(PeReq*(-1 + Pres))/Deltainf` |

The engines are in algebraic agreement. Whether that agreement is meaningful is undercut by F2 (transliteration).

## Verdict justification

The script's load-bearing chain - matched and profile-family threshold derivations via Solve on Peclet balances, followed by cross-relations and resonance specialisation - is faithful to the paper card and notes. The v1 material change (lifting `Wfail_res / Wfail_match` from postulated values to Solve-derived expressions) is intact and produces exactly the exported quantities the notes require: `Pe_req/Delta_inf`, `Pe_req/Delta_0`, `Pe_req/(C^2 Delta_inf)`, `Pe_req/(C^2 Delta_0)`. Paper alignment is therefore `aligned`.

However, three of the script's load-bearing "derivations" near the top - `W_res = C^2 W_wall`, `P_res = 1/C_res^2`, and `P_res*C_res^2 = 1` - are tautological label substitutions, not substantive checks. The Mathematica script is a line-by-line port of the SymPy choreography (with a stray no-op `gapSym;`), violating the dual-engine independence policy. And the "two ways" band-width check is redundant rather than independently reinforcing.

The numeric value `P_res ~ 1.005612487760576` is not exercised in-script. This is acceptable as a carry-forward from stage 067 (whose script verifies `C_res^2`) but the absence is worth noting in the resolution of F1(c), where I prescribe an in-script numeric anchor.

Attacks I tried that failed: I checked whether the `positive=True` assumptions hid a sign-branch issue (they don't - every quantity here is genuinely positive by physical setup and the paper). I checked whether the symbol `C2` correctly represents `C^2` and not `C` (it does - every multiplication and substitution treats it consistently). I checked whether `Pe_req` was being conflated with anything from stage 066 in a way that would change its meaning (no - it's a free symbol throughout, as is appropriate for stage 068). I checked whether the resonance substitution `C2 -> 1/Pres` matches the notes (it does - notes section 2 sets `P_res = 1/C_res^2` at the resonance point `r = sqrt(pi)`).

Verdict: `findings` (3 findings). No stop-cold flag - none of the findings would invalidate downstream stages. Stage 069 uses the threshold values `P_res * Pe_req / Delta_{0,inf}` symbolically; the script's correct delivery of these is independent of the tautological top-of-script issues.

## v1 material_change carry-forward

The user's note flags that stage 068 carries a v1 `material_change: true` flag because v1 replaced the postulated `Wfail_res / Wfail_match` with Solve-derived expressions on the resonance-corrected premises. The v2 re-audit specifically validates this premise-replacement:

- SymPy lines 84 (`W_match_sol = sp.solve(...)`) and 92 (`W_prof_sol = sp.solve(...)`) confirm the two thresholds are now Solve-derived from the matched and profile-family Peclet balances respectively, not assigned from a literal.
- Mathematica lines 46-48 and 51-53 mirror this (each via `First@Solve[...]`).
- The resulting expressions (`Pe_req/Delta_inf`, `Pe_req/Delta_0`, `Pe_req/(C^2 Delta_inf)`, `Pe_req/(C^2 Delta_0)`) match notes section 2 verbatim and match the paper card's `\stagefield{Output}` requirement that `P_res` multiplies the matched window.
- No regression introduced by the v1 fix; the threshold-translation portion of the audit is substantively solid.

The remaining v2 findings (F1, F2, F3) are pre-existing weaknesses in the parts of the script that were not touched by the v1 material change - the `W_res = C^2 W_wall` and `P_res = 1/C_res^2` "derivations" and the band-width "two ways" check. None affect the v1 fix; all three predate it.

## Self-test notes

- **Variable independence**: there are no `sp.diff` calls in this script, so the derivative-tautology trap does not apply.
- **Symmetry/parity**: no integrals over unbounded domains in this script; the parity trap does not apply.
- **Trivial-case pre-check**: For the F1 prescribed `Wres_expr = kappa * Gres_expr` reformation, substituting `kappa, rho_star, g_phi, N_phiphi, m, c_s, K_X` all positive reals and `C2 = 1` should give `Wres_expr - Wwall_expr = 0` (matched limit). For the F3 prescribed `Cres -> 2*Cres` perturbation, `Pres -> 1/(4*Cres^2)` and the band-width assertion would fail by a factor of 4 - confirming the assertion is now sensitive to `Cres`.
- **Path specifications**: the directive does not introduce any new scripts; it only modifies existing files at `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py` and `mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`.
- **Paper round-trip**: the F1 fix introduces `kappa, rho_star, g_phi, N_phiphi, m, c_s, K_X` from notes section 1; these are matched-branch gain components and are already named in the notes - no new paper_misalignment risk. The F1(c) numeric anchor uses the paper card's quoted constant `1.005612487760576` verbatim. The F3 fix routes through `Cres` which is the upstream-anchored symbol from stage 067.
