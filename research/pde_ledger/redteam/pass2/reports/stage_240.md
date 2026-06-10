---
unit_id: 240
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 240 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_240.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows 92, 1179-1235, claimstatus 1381)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.txt`

## What the paper claims

Stage 240 fixes the support-loading ratio selected by the minimal isotropic quadrupole precursor. `\stagefield{Output}`: "Exact selected loading ratio: $\rho_\alpha=4/3$, $\zeta_{\rm req}=1/3$, and $\Pi_{\rm tr}=4C_{\rm mix}/3$, placing the branch on the symmetric lowest-twin support slice." The notes add the full derivation: (1) selected-branch product identities cancel the normalization amplitudes so $\Pi_{\rm tr}/C_{\rm mix}=\alpha_{\rm req}/\alpha_{\rm mix}=:\rho_\alpha$; (2) the contact-plus-pole inverse compiler $c_0=1/\rho_\alpha$, $c_1=(\rho_\alpha-1)/\rho_\alpha$ with inverses $\rho_\alpha=1/c_0=1/(1-c_1)$, $\zeta_{\rm req}=c_1/c_0$; (3) the minimal isotropic module $c_0=3/4,c_1=1/4$ giving $\rho_\alpha=4/3,\zeta_{\rm req}=1/3$; (4) $\Pi_{\rm tr}=(4/3)C_{\rm mix}$, $S_{\rm req}=4/3$; (5) the support-selector reduction $\varrho=2(1-\epsilon_*)/3$; (6) the regime classification $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$ (symmetric lowest twin). Notes §2.1 explicitly states the pole location $\Omega_Q$ does NOT enter the static loading-ratio extraction. The card's `\stagefield{Verification}` says "Mathematica audit: none yet."

## What the script claims to verify

The SymPy and Mathematica scripts assert all six deliverables: the product-ratio identity (plus a spectral-substitution variant), the contact-plus-pole compiler and its inverse formulas, an explicit extraction of $(c_0,c_1)$ from the $\Omega_Q$-bearing precursor $Y_{\rm support}$ together with an $\Omega_Q$-independence check on the extracted weights, the minimal-isotropic specialization to $4/3$ and $1/3$, $\Pi_{\rm tr}=(4/3)C_{\rm mix}$ and $S_{\rm req}=4/3$, the selector reduction $\varrho=2(1-\epsilon_*)/3$, and the strict regime inequality $1<4/3<2$ ⇔ $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| $\Pi_{\rm tr}/C_{\rm mix}=\alpha_{\rm req}/\alpha_{\rm mix}=\rho_\alpha$ | py L64-76 / wl M1 | match |
| compiler $c_0=1/\rho_\alpha$, $c_1=(\rho_\alpha-1)/\rho_\alpha$ + inverses | py L86-132 / wl M2 | match |
| $\Omega_Q$ does not enter static extraction (notes §2.1) | py L113-130 / wl L76-105 | match |
| minimal $c_0=3/4,c_1=1/4 \Rightarrow \rho_\alpha=4/3,\zeta_{\rm req}=1/3$ | py L137-146 / wl M3 | match |
| $\Pi_{\rm tr}=(4/3)C_{\rm mix}$, $S_{\rm req}=4/3$ | py L151-158 / wl M4 | match |
| $\varrho=2(1-\epsilon_*)/3$ | py L163-175 / wl M5 | match |
| regime $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$ | py L180-188 / wl M6 | match |
| card Verification "Mathematica audit: none yet" | `.wl` is present and passing | mismatch (card stale) |

`paper_alignment: partial` — every Output value reconciles exactly; the only discrepancy is the card's stale "Mathematica audit: none yet" claim while a passing `.wl` exists (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 64-67 | `simplify(Pi/C - rho)==0` | product ratio | yes |
| A2 | sympy | 73-76 | spectral ratio ==0 | product ratio (spectral) | yes |
| A3 | sympy | 86 | `Y_support - Y_rho ==0` | compiler form | yes |
| A4 | sympy | 91-94 | inverse formulas ==0 | compiler inverses | yes |
| A5 | sympy | 118-132 | limit-extraction + `diff(.,Omega_Q)==0` + match | $\Omega_Q$-independence (F1 fix) | yes |
| A6 | sympy | 144-146 | minimal $\Rightarrow 4/3,1/3$ | minimal specialization | yes |
| A7 | sympy | 151-158 | $\Pi=(4/3)C$, $S=4/3$ | demand product | yes |
| A8 | sympy | 163-175 | $\varrho=2(1-\eps)/3$ | selector reduction | yes |
| A9 | sympy | 180-188 | regime $1<4/3<2$ | regime classification | yes (literal) |
| M1-M6 | mathematica | 58-142 | `expectZero`/`expectTrue` mirror | all six | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim (card understates verification coverage)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_240.tex:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.wl`

**What's wrong:**
The card states (line 11): `\stagefield{Verification}{SymPy audit: ... Mathematica audit: none yet.}` But a Mathematica audit `.wl` exists, is committed, and its saved output passes all six modules (M1-M6). The card's verification line is stale relative to the dual-engine reality.

**Why this matters:**
Dual-engine rule is in force; the card should reflect that both engines now verify this stage. This is documentation-side only — no result value is affected.

**Required change:**
None for Codex (Codex must not edit paper/). Route to user: update the card's `\stagefield{Verification}` line to cite the existing `.wl`.

**Verification:**
Card line 11 references the `.wl` filename once updated.

## Independent-derivation check (Mathematica)

The `.wl` is parallel in structure (M1-M6 mirror the six SymPy sections) but is NOT a pure transliteration. The weight extraction uses a genuinely different algebraic route: the `.wl` first computes `yApart = Apart[ySupport, omega]` (visible in output line 10 as an explicit partial-fraction decomposition into `omega-omegaQ` and `omega+omegaQ` simple poles), then extracts `c1Extracted = Limit[poleFactor*yApart, omega->omegaQ]` and `c0Extracted = Limit[yApart - c1Extracted/poleFactor, omega->0]`. The SymPy script instead multiplies the un-decomposed `Y_support` by the pole denominator and takes the limit directly (no `apart`). Different decomposition path, same physics. The compiler/inverse algebra and the regime test use independent Mathematica primitives (`FullSimplify`, `Resolve`). Not a transliteration finding.

## Engine cross-check

Both outputs agree at the claimed level. SymPy emits 23 `[ok]` lines + "All Stage 240 symbolic checks passed"; Mathematica emits PASS for M1-M6 and "Stage 240 Mathematica audit passed." The `.wl` output explicitly shows `c0Extracted = alphaMix/alphaReq`, `c1Extracted = 1 - alphaMix/alphaReq` (output L11-12), matching the SymPy extraction `c0_static = alpha_mix/alpha_req`, `c1_static = (alpha_req-alpha_mix)/alpha_req`. No sign or factor disagreement.

## Verdict justification

Every paper Output value ($\rho_\alpha=4/3$, $\zeta_{\rm req}=1/3$, $\Pi_{\rm tr}=(4/3)C_{\rm mix}$, $\varrho=2(1-\epsilon_*)/3$, $S_{\rm req}=4/3$, regime $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$) is exercised by a non-tautological, well-anchored check in BOTH engines, and all reconcile against the card, notes, and appendix. The pass-1 F1 variable-independence self-test trap is genuinely FIXED: the extracted weights are now pulled from `Y_support` (which carries $\Omega_Q$ via its pole), the script guards (py L115-116, 119-120, 124-125; wl L77, L84) confirm the extraction path actually sees $\Omega_Q$ before the limit, and the subsequent `diff(c0_static, Omega_Q)` / `diff(c1_static, Omega_Q)` exercise the real notes-§2.1 claim that $\Omega_Q$ drops out of the static extraction (it could have survived had the precursor weighted the pole by $\Omega_Q$). The only finding is a low-severity stale card Verification line. Attacks tried that failed: hardcoded-result (the $3/4,1/4 \to 4/3,1/3$ chain is derived through the inverse compiler, not asserted as a bare literal); tautology on the extraction (the limit path traverses $\Omega_Q$-bearing objects); symbol-assumption (positivity / $\alpha_{\rm req}>\alpha_{\rm mix}$ in the `.wl` matches the physical setup, the `.py` keeps weaker assumptions and still passes); engine disagreement (none). Verdict: findings (1, low, paper-side documentation only).

## Self-test notes

Checked variable-independence trap (step 1): every `sp.diff`/`D[]` in the F1 block (py L129-130, wl L104-105) differentiates `c0_static`/`c1_static` w.r.t. `Omega_Q`; these are extracted from `Y_support` which DOES carry `Omega_Q`, and the post-limit Omega_Q-free result is the genuine physics claim, not a trap — guards force the pre-limit object to contain `Omega_Q`. Checked parity (step 2): no unbounded integrals, n/a. Checked trivial-case (step 3): substituting the minimal module gives the literal $4/3$/$1/3$ as the outputs show. No directive of script edits is warranted; the lone finding is paper-side, routed to the user.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| $\rho_\alpha=4/3$ | py L144 / wl M3 (out L29-31) | stage_240.tex:15; notes L17,265; appendix:1208 | MATCH |
| $\zeta_{\rm req}=1/3$ | py L146 / wl M3 (out L33) | stage_240.tex:15; notes L19,273; appendix:1210 | MATCH |
| $\Pi_{\rm tr}=(4/3)C_{\rm mix}$ | py L153 / wl M4 (out L35) | stage_240.tex:15; notes L21,285; appendix:1215 | MATCH |
| $S_{\rm req}=4/3$ | py L158 / wl M4 (out L37) | notes L32,321; appendix:1229 | MATCH |
| $\varrho=2(1-\epsilon_*)/3$ | py L166 / wl M5 (out L39) | notes L31,311; appendix:1227 | MATCH |
| $c_0=1/\rho_\alpha$, $c_1=(\rho_\alpha-1)/\rho_\alpha$ | py L88-89 / wl M2 (out L11-12) | notes L214-216; appendix:1206-1208 | MATCH |
| minimal $c_0=3/4,c_1=1/4$ | py L137-138 / wl L108-109 | notes L253-256; appendix:1199-1201 | MATCH |
| regime $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$ | py L186 / wl M6 (out L43) | stage_240.tex:15; notes L361; appendix:1234 | MATCH |

Non-checkpoint constants per audit note all reconcile: $\rho_\alpha=4/3$ (MATCH), $\zeta_{\rm req}=1/3$ (MATCH), $\Pi_{\rm tr}=(4/3)C_{\rm mix}$ (MATCH).

INTERNAL scaffolding (no finding): `pole`/`poleFactor`, `pole_denominator`, `c1_probe`/`poleProbe`, `static_sum_probe`/`contactProbe`, `Y_rho`, `NQ_target_spec`/`nqSpectral`, `C_mix_Lambda`/`cMixSelector`, pass/fail flags, the `.has(Omega_Q)`/`FreeQ` guard residuals.

reconciliation: complete; 8 deliverable values checked, 0 misaligned (the one paper-side finding is a stale Verification-coverage line, not a result-value mismatch).
