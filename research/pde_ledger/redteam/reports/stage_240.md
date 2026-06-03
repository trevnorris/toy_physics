---
unit_id: 240
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 240 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_240.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 92, 200-205, 1179-1235, 1378-1400)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_projects... -> /var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The card's `\stagefield{Output}` states verbatim: "Exact selected loading ratio: $\rho_\alpha=4/3$, $\zeta_{\rm req}=1/3$, and $\Pi_{\rm tr}=4C_{\rm mix}/3$, placing the branch on the symmetric lowest-twin support slice." The notes and the Part VII appendix (eqs. `app-part07-selected-product-identities` through `app-part07-twin-window`) enumerate the deliverables: (1) the product-ratio identity $\Pi_{\rm tr}/C_{\rm mix}=\alpha_{\rm req}/\alpha_{\rm mix}=:\rho_\alpha$ (also in a spectral form via $N_Q^{(\rm target)}=\hat m_-^2\beta_0 s_-/\lambda_-$); (2) the contact-plus-pole compiler $c_0=1/\rho_\alpha,\ c_1=(\rho_\alpha-1)/\rho_\alpha$ with inverses $\rho_\alpha=1/c_0=1/(1-c_1),\ \zeta_{\rm req}=c_1/c_0$; (3) the minimal isotropic specialization $c_0=3/4,\ c_1=1/4 \Rightarrow \rho_\alpha=4/3,\ \zeta_{\rm req}=1/3$; (4) $\Pi_{\rm tr}=(4/3)C_{\rm mix}$, $S_{\rm req}=4/3$; (5) the support-selector reduction $\varrho=\pi^2\Pi_{\rm tr}/(16\Lambda)=2(1-\epsilon_*)/3$ using $C_{\rm mix}=8\Lambda(1-\epsilon_*)/\pi^2$; (6) the regime classification $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$ (symmetric lowest twin). The notes §2.1 additionally make a specific physics assertion: the pole location $\Omega_Q$ controls the dynamical shape but does **not** enter the static loading-ratio extraction; only the normalized weights $(c_0,c_1)$ matter.

## What the script claims to verify

The SymPy docstring (lines 2-26) and assertions cover all six paper deliverables. It builds `rho_alpha = alpha_req/alpha_mix` symbolically, checks the product-ratio cancellation (lines 64-76, including the spectral substitution), the contact/pole compiler rewrite `Y_support == Y_rho` (line 86), the normalization `c0+c1=1` and the inverse formulas (lines 91-111), the `Omega_Q`-independence of `c0,c1` (lines 114-115), the minimal-module specialization to `4/3` and `1/3` (lines 127-129), `Pi_tr=(4/3)C_mix` and `S_req=4/3` (lines 134-141), the `varrho=2(1-eps_*)/3` reduction (lines 146-158), and a numeric regime inequality `1<4/3<2` plus `C_mix<Pi_tr<2C_mix` at `C_mix=1` (lines 163-171). All 21 checks PASS (output exit 0).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) product-ratio identity (+ spectral form) | lines 64-76 | match |
| (2) contact/pole compiler + inverses | lines 86, 91-111 | match |
| (2, §2.1) $\Omega_Q$ drops out of the static extraction | lines 114-115 | mismatch-in-substance (tautological; see F1) |
| (3) minimal module $c_0=3/4,c_1=1/4\Rightarrow 4/3,1/3$ | lines 127-129 | match |
| (4) $\Pi_{\rm tr}=(4/3)C_{\rm mix}$, $S_{\rm req}=4/3$ | lines 134-141 | match (but constructed-then-confirmed; weak, see Verdict) |
| (5) $\varrho=2(1-\epsilon_*)/3$ | lines 146-158 | match |
| (6) regime $C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}$ | lines 163-171 | match |
| second engine (Mathematica) | none | missing (F2) |

`paper_alignment: aligned` — every constant in the script (`3/4, 1/4, 4/3, 1/3, 8 Lambda(1-eps_*)/pi^2, pi^2/(16 Lambda)`) matches the card Output, the notes, and the appendix equations exactly. No value or target mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 64-67 | `simplify(Pi/C - rho)==0` | claim 1 | yes |
| A2 | sympy | 73-76 | `simplify(Pi_spec/C_spec - rho)==0` | claim 1 (spectral) | yes |
| A3 | sympy | 86 | `simplify(Y_support - Y_rho)==0` | claim 2 | yes |
| A4 | sympy | 91 | `c0+c1-1==0` | claim 2 | yes |
| A5 | sympy | 92-94 | inverse formulas | claim 2 | yes |
| A6 | sympy | 100-111 | symbolic-`c0/c1` specializations | claim 2 | yes (redundant w/ A5) |
| A7 | sympy | 114-115 | `diff(c0_expr, Omega_Q)==0` | claim 2 / §2.1 | **no — tautological** |
| A8 | sympy | 127-129 | `1/(3/4)-4/3==0`, `(1/4)/(3/4)-1/3==0` | claim 3 | yes |
| A9 | sympy | 134-141 | `(4/3)C - (4/3)C==0`, `S_req-4/3==0` | claim 4 | partial — constructed-then-confirmed |
| A10 | sympy | 146-158 | `varrho - 2(1-eps_*)/3 == 0` | claim 5 | yes |
| A11 | sympy | 163-171 | `1<4/3<2`, `C<Pi<2C` at `C=1` | claim 6 | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py:113-115`

**What's wrong:**
The notes §2.1 make a concrete physics claim (lines 238-243): "The pole location $\Omega_Q$ controls the dynamical conservative shape of the precursor, but it does **not** enter the static loading-ratio extraction. At this stage only the normalized contact and pole weights $(c_0,c_1)$ matter." The script tries to verify this at lines 113-115:

```python
# Omega_Q does not affect the static loading-ratio extraction.
assert_zero(sp.diff(c0_expr, Omega_Q), "Omega_Q independence of c0")
assert_zero(sp.diff(c1_expr, Omega_Q), "Omega_Q independence of c1")
```

But `c0_expr = sp.simplify(1/rho_alpha) = alpha_mix/alpha_req` (line 88) and `c1_expr = (alpha_req-alpha_mix)/alpha_req` (line 89). Neither expression contains the symbol `Omega_Q` at all — they were *constructed* from `alpha_req, alpha_mix` only. Therefore `sp.diff(c0_expr, Omega_Q)` is identically `0` by construction, and the `assert_zero` cannot fail no matter what the physics is. This is the textbook variable-independence trap: differentiating with respect to a variable the expression does not depend on. The assertion provides zero evidence for the §2.1 claim.

The genuine content of §2.1 is that even though the *full* precursor `Y_support` / `Y_rho` (line 83-84) carries `Omega_Q` through `pole = 1/(1 - omega**2/Omega_Q**2)`, the *static* weights extracted from it (the `omega`-independent contact part and the residue/pole weight) are `Omega_Q`-free. The script already has the `Omega_Q`-dependent objects in hand (`Y_support`, `pole`); the check should extract `c0,c1` *from* `Y_support` (e.g. the `omega=0` contact split vs. the pole coefficient) and show the extracted weights carry no `Omega_Q`, rather than differentiating an expression that was written without `Omega_Q` to begin with.

**Why this matters:**
A reader (or downstream auditor) sees `[ok] Omega_Q independence of c0` in the transcript and concludes the §2.1 physics claim is machine-verified. It is not — it is a no-op. If the static-extraction algebra were actually `Omega_Q`-dependent (a real error mode this check is meant to guard against), this assertion would still pass. The check is decorative.

**Required change:**
Replace the two `sp.diff(..., Omega_Q)` assertions with a check that genuinely exercises the static extraction from the `Omega_Q`-bearing precursor. Concretely: extract the contact fraction and the pole weight *from* `Y_support` (the line-83 expression that does contain `Omega_Q`) — e.g. the `omega -> 0`-vs-pole decomposition or the partial-fraction split of `Y_support` in `omega` — obtain candidate `c0_static, c1_static`, then assert (a) `c0_static` and `c1_static` are each free of `Omega_Q` (`assert_zero(sp.diff(c0_static, Omega_Q))` where `c0_static` was extracted from an `Omega_Q`-bearing object), and (b) `c0_static - c0_expr == 0`, `c1_static - c1_expr == 0`. The directive states the requirement; Codex designs the extraction.

**Verification:**
The verifier confirms (i) the new check derives `c0_static`/`c1_static` from an expression that syntactically contains `Omega_Q` (so the derivative is non-trivially zero), and (ii) the labels `Omega_Q independence ...` still appear and the script exits 0.

### F2 — missing_verification_script (missing_mathematica)

**Severity:** medium
**Files:**
- `(missing)` — no `.wl` for unit 240

**What's wrong:**
Unit 240 has a SymPy script and no Mathematica script. The stage is `is_checkpoint: False`, `is_status_only_candidate: False`. Every claim in this stage is a closed-form algebraic / rational-function identity: a ratio cancellation, a partial-fraction (contact-plus-pole) rewrite, rational substitutions, and a numeric interval check. Mathematica can independently verify all of these via native primitives (`Simplify`/`FullSimplify`, `Apart`/`Together`, `Series`/`Limit` for the static extraction, `Reduce`/numeric comparison for the regime interval). The test is "is it possible" — it clearly is — so a second engine is required.

**Why this matters:**
With only one engine, there is no cross-check on the algebra; a SymPy-specific simplification quirk (or the F1 no-op) would go undetected. The dual-engine policy requires an independent route wherever Mathematica *can* verify.

**Required change:**
Codex writes a NEW independent-route Mathematica script (see directive F2 for the path, claim manifest, and the independence constraint). It must derive the results via a *different decomposition* than the `.py`, not transliterate it.

**Verification:**
`redteam exec-mathematica 240` produces the new `.wl`'s output, the script exits 0, and each manifest claim M1-M6 has a corresponding in-file check (`expectZero[...]` / `If[... Exit[1]]`).

## Independent-derivation check (Mathematica)

No `.wl` exists, so `mathematica_transliteration` does not apply yet. The directive's F2 explicitly constrains the new script to an independent route (different decomposition) and forbids a line-by-line port of the SymPy choreography.

## Engine cross-check

Only one engine present; `engines_agree: n/a`.

## Verdict justification

`findings`. Paper alignment is exact: all six deliverables map to script checks and every constant matches the card/notes/appendix — I attacked the constants (3/4, 1/4, 4/3, 1/3, the `8 Lambda(1-eps_*)/pi^2` and `pi^2/(16 Lambda)` factors) and they hold, so no `paper_misalignment`. Two real script-side defects survive: (F1) the `Omega_Q`-independence assertions are a variable-independence no-op that gives false credit to the §2.1 physics claim, and (F2) the stage lacks a required second engine even though every claim is Mathematica-verifiable. Note also a non-blocking weakness in A9 (lines 134-141): `Pi_tr_selected` is defined as `(4/3)*C_mix` and then asserted equal to `(4/3)*C_mix`, so those two checks are tautological re-statements of the already-verified `rho=4/3` (line 127); I did not raise a separate finding because they faithfully restate the paper's claim 4 and the substantive content is genuinely tested upstream at line 127 — but the new Mathematica engine (F2) should test claim 4 substantively rather than mirror this redundancy.

## Self-test notes

I walked the variable-independence trap on every `sp.diff`: the only derivatives are lines 114-115 w.r.t. `Omega_Q`, and `c0_expr`/`c1_expr` do not contain `Omega_Q` — confirmed identically-zero by construction (this is F1). No integrals exist, so the symmetry/parity trap is N/A. Trivial-case pre-check: the regime inequality at lines 163-171 substitutes concrete `4/3` and `C_mix=1`, giving `1<4/3<2` and `1<4/3<2` — genuinely a numeric pass that could fail if the ratio left (1,2), so it is sound. Path spec for F2 confirmed: `.wl` lives in `mathematica/` (224 siblings, `_mathematica_audit.wl` suffix). Paper round-trip: the F1 fix introduces no new constant (it re-derives the same `c0=alpha_mix/alpha_req`), so no new `paper_misalignment`.
