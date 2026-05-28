---
unit_id: 165
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T15:20:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage165_exact_branch_drifts.md]
  paper_appendix: present
---

# Audit unit 165 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_165.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage165_exact_branch_drifts.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 61, 334-338)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.txt`

## What the paper claims

The stage card `\stagefield{Output}` states verbatim: "Shows \(\delta\ln L_W=\delta\ln a\) and reduces the lower compensated branch to four irreducible drifts \((\delta\ln\mathcal Z_q,\delta\ln\rho_w,\delta\ln c_s,\delta\ln a)\)." The notes expand this into the full derivation: (D1) the D/N selection law \(L_W=\tfrac{\pi a}2\sqrt{(1+\mathfrak r^2)/3}\) gives \(\delta\ln L_W=\delta\ln a\) at fixed \(\mathfrak r\); (D2) the fixed-\(\mathfrak r\) and fixed-\(\mathfrak g\) conditions [\(\delta\ln(K_sK_q/\lambda^2)=0\) and \(\delta\ln(g_qK_s/(g_s\lambda))=0\)] solve for \(\delta\ln v_{w0}\) and \(\delta\ln\mathcal T_m\); (D3) the product/ratio factorization \(\delta\ln(v_{w0}/\mathcal T_m)=2\delta\ln c_s-\delta\ln a\) and \(\delta\ln(v_{w0}\mathcal T_m)=\delta\ln\mathcal Z_q+3\delta\ln c_{s,w}-\delta\ln\rho_w-4\delta\ln a\); (D4) the frozen \(n=5\) wall EOS \(\delta\ln c_{s,w}=2\delta\ln\rho_w\) collapsing the product law to \(\delta\ln\mathcal Z_q+5\delta\ln\rho_w-4\delta\ln a\); (D5) the Stage-249 off-family normal coordinate vanishes identically, \(\delta_\perp=0\), because the off-family channels ARE the fixed-\(\mathfrak g\)/\(\mathfrak r\) conditions. The appendix row (line 61) and lines 334-338 repeat the headline co-transport law.

## What the script claims to verify

The docstring lists five checks matching D1-D5. The load-bearing assertions are: (a) Mathematica `a*D[Log[lwLaw],a]-1 == 0` (D1, Mathematica only); (b) the SymPy/Mathematica solve of the 2x2 linear system {eq_r, eq_g} for (dv, dT), printed but not asserted; (c) `expectZero` on `ratio_law - (2*dcs - da)` and `prod_law - (dZ + 3*dcsw - drho - 4*da)` (D3); (d) `expectZero` on the n=5 product law minus `(dZ + 5*drho - 4*da)` (D4); (e) `expectZero` on `channel_g` and `channel_r` after substituting the solved (dv, dT) (D5). The SymPy script's D1 step is a deliberate no-op (`if False else 1`).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: δlnL_W = δln a | Mathematica line 33 (genuine); SymPy line 41 disabled (`if False`), prints text only | partial (SymPy missing) |
| D2: solve dv, dT from fixed-r/fixed-g | both engines solve [eq_r, eq_g] for (dv, dT); forms match notes boxes | match |
| D3: ratio law / product law | `expectZero` against independent target forms (both engines) | match |
| D4: n=5 EOS product collapse | `expectZero` against `(dZ + 5 drho - 4 da)` (both engines) | match |
| D5: δ⊥ = 0 (off-family channels vanish) | `channel_g`/`channel_r` = subs(solution into eq_g/eq_r) → 0 by construction | match-but-tautological |

`paper_alignment: aligned` — every deliverable maps to a check, every constant and sign in eq_r/eq_g and the target forms matches the notes exactly. The two findings are script-side rigor gaps, not paper↔script disagreements.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41 | `... if False else 1` (no-op, prints text) | D1 | no (disabled) |
| A2 | mathematica | 33 | `expectZero[a*D[Log[lwLaw],a]-1]` | D1 | yes |
| A3 | sympy | 66 | `expect_zero(ratio_law - (2*dcs - da))` | D3 ratio | yes |
| A4 | sympy | 67 | `expect_zero(prod_law - (dZ+3*dcsw-drho-4*da))` | D3 product | yes |
| A5 | sympy | 75-78 | `expect_zero(n=5 product - (dZ+5*drho-4*da))` | D4 | yes |
| A6 | sympy | 83 | `expect_zero(channel_g)` | D5 | no (tautological) |
| A7 | sympy | 84 | `expect_zero(channel_r)` | D5 | no (tautological) |
| A8 | mathematica | 54-55 | `expectZero[ratioLaw...]`, `expectZero[prodLaw...]` | D3 | yes |
| A9 | mathematica | 62-65 | `expectZero[n=5 product...]` | D4 | yes |
| A10 | mathematica | 70-71 | `expectZero[channelG]`, `expectZero[channelR]` | D5 | no (tautological) |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:41`

**What's wrong:**
Deliverable D1 — the stage's headline `\stagefield{Output}` claim \(\delta\ln L_W=\delta\ln a\) — is the most-cited result of this stage (appendix eq. `eq:app-part05-LW-co-transport`, line 337). The SymPy script does not verify it. Line 41 reads:
```python
dlog_LW = sp.simplify(sp.diff(sp.log(LW_law), sp.log(a)) if False else 1)  # documentary only
```
The `if False` guard means SymPy always assigns the literal `1` and never differentiates `LW_law`. Line 43 then prints the claim "At fixed r_* : d ln L_W = d ln a" as plain text. A reader of the SymPy transcript would believe SymPy confirmed D1; it did not. The Mathematica script does verify it correctly at line 33 (`expectZero["d ln L_W - d ln a at fixed r_*", a*D[Log[lwLaw], a] - 1]`), so the *paper* claim is exercised by one engine, but the SymPy coverage of the headline deliverable is a no-op print, not an assertion. (Note: `sp.diff(log(LW_law), sp.log(a))` as originally written would itself be wrong — SymPy `diff` w.r.t. a composite `log(a)` does not compute a logarithmic derivative; the correct form is `a*sp.diff(sp.log(LW_law), a)`. That is presumably why the author short-circuited it, but short-circuiting to a hardcoded `1` rather than fixing the derivative leaves the check absent.)

**Why this matters:**
The single load-bearing output of the stage has no SymPy assertion. If `LW_law` were ever edited to a form with a different scaling exponent in `a` (e.g. a stray `a^2`), the SymPy script would still print "d ln L_W = d ln a" and pass. The check that would catch such a regression exists only in the second engine.

**Required change:**
Replace the no-op assignment with a genuine logarithmic-derivative check matching the Mathematica form. Add an `expect_zero` (not just a print) on `a*sp.diff(sp.log(LW_law), a) - 1`, so SymPy confirms `LW_law` has scaling degree exactly 1 in `a`. Concretely, replace line 41 with:
```python
dlog_LW = sp.simplify(a*sp.diff(sp.log(LW_law), a))
print("D/N law: L_W =", LW_law)
print("At fixed r_* : d ln L_W = d ln a")
expect_zero("d ln L_W - d ln a at fixed r_*", dlog_LW - 1)
```
(Remove the duplicate `print("D/N law: L_W =", LW_law)` at line 42 if both would appear.)

**Verification:**
After the fix, the SymPy transcript should contain a new line `d ln L_W - d ln a at fixed r_* = 0` (the `expect_zero` print), and the script must still exit 0. The residual `a*d/da log(πa√(1+r²)/(2√3)) - 1 = a·(1/a) - 1 = 0`, so the assertion passes for the correct law and would fail if the `a`-exponent of `LW_law` were anything other than 1.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:81-84`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:68-71`

**What's wrong:**
The "Stage 164 channel closure" check (D5, δ⊥=0) is algebraically guaranteed by construction and cannot fail. `dv_sol`/`dT_sol` are obtained from `sp.solve([eq_r, eq_g], [dv, dT])`. The two channel expressions checked,
```python
channel_g = (dZ + 3*dcsw - drho - dT - dv - 2*da - 2*dLW).subs({dv: dv_sol, dT: dT_sol})
channel_r = (dZ + 2*dcs + 3*dcsw - drho - 2*dv - 2*da - 3*dLW).subs({dv: dv_sol, dT: dT_sol})
```
are *exactly the left-hand sides of `eq_g` and `eq_r`*. Substituting the solution of a linear system back into the same system's equations yields 0 by definition of "solution." No physics is exercised: the check verifies that `sp.solve` returns a correct solution, not that δ⊥=0 holds for any independent reason. The Mathematica mirror (lines 68-71, `/. sol`) is identically tautological.

**Why this matters:**
The notes (section 8, lines 386-405) make a genuine claim — that the Stage-249 off-family normal coordinate vanishes *because* the off-family channels coincide with the branch conditions. The script's check does not test that coincidence; it tests `solve()`. The notes' own framing ("the branch drift constraints themselves already compute the cancellation," line 404-405) means the result is true essentially by definition once the channels are identified with eq_r/eq_g — so this is a low-severity rigor gap, not a wrong result. It is flagged so the reader does not over-credit the transcript's `PASS: fixed-g channel` / `PASS: fixed-r channel` lines as independent verification.

**Required change:**
Make the closure check substantive by exercising the *definition* of δ⊥ rather than re-substituting the solution into its own defining system. Construct δ⊥ as a stated linear combination of the two off-family channels (per Stage 249's normal-coordinate definition, e.g. a fixed-coefficient combination of `δln(g_qK_s/(g_sλ))` and `δln(K_sK_q/λ²)`) and assert that this combination evaluates to 0 *only when both branch conditions hold* — i.e. solve the branch conditions and confirm δ⊥ vanishes, while also confirming δ⊥ is NOT identically zero off the branch (an `assert_nonzero` on δ⊥ evaluated with `dv`, `dT` left free). If the precise δ⊥ coefficients are not available from this stage's notes, downgrade the two channel lines to clearly-labelled `print` statements (drop the `expect_zero`/`expectZero`) so the transcript does not assert a tautology as a PASS.

**Verification:**
After the fix, either (i) a new `expect_zero` on the constructed δ⊥(branch) AND a new `assert`/nonzero-print on δ⊥(free) appear in both transcripts, or (ii) the two `fixed-g channel`/`fixed-r channel` PASS lines are replaced by plain prints. Script must still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is a line-by-line transliteration of the `.py`: identical symbol set (`dZ, drho, dcsw, dcs, dT, dv, da, dLW, r, a`), identical `eqR`/`eqG` definitions (`.wl:35-36` vs `.py:46-47`), identical `Solve[{eqR,eqG},{dv,dT}]` vs `sp.solve([eq_r,eq_g],[dv,dT])`, identical `/. dLW -> da` vs `.subs(dLW, da)`, identical `/. dcsw -> 2*drho` vs `.subs(dcsw, 2*drho)`, the same four `expectZero` targets in the same order, and even the same copy-paste artifact banner `"STAGE 148 — EXACT LOWER-BRANCH DRIFT LAWS"` (`.wl:26` and `.py:30`) plus an identical carry-forward `Print` block. This is `mathematica_transliteration` in the strict sense. However, the entire stage is the solution of one 2x2 linear system in log-drifts followed by linear substitutions; there is essentially a unique procedure, so the "independent re-derivation" bar carries little physical content here. I am NOT raising it as a separate numbered finding because the transliteration does not weaken any *correctness* claim (linear algebra has no branch/sign ambiguity for the second engine to independently catch), and fixing F1/F2 is the higher-value action. Noted for the record.

## Engine cross-check

Both engines agree. SymPy output line 17 `d ln v_w0 = -3*dLW/2 + dZ/2 - da + dcs + 3*dcsw/2 - drho/2` equals Mathematica line 20 `(-2*da + 2*dcs + 3*dcsw - 3*dLW - drho + dZ)/2`; line 18 equals Mathematica line 21. All `expectZero`/`expect_zero` residuals are 0 in both transcripts (ratio, product, n=5 product, both channels). The n=5 forms match (`dZ/2 - 5*da/2 + dcs + 5*drho/2` etc.). No `engine_disagreement`.

I independently re-derived the 2x2 solve: from eq_r, `dv = (dZ + 2dcs + 3dcsw − drho − 2da − 3dLW)/2`; back-substituting into eq_g gives `dT = (dZ + 3dcsw − drho − 2dcs − 2da − dLW)/2`. Both match the script outputs and the notes' boxed forms (sections 3-4) exactly, including signs. The ratio `dv_DN − dT_DN = 2dcs − da` and product `dv_DN + dT_DN = dZ + 3dcsw − drho − 4da` reproduce the notes' section-6 boxes; the n=5 substitution `dcsw → 2drho` reproduces the section-7 boxes.

## Verdict justification

`findings`. The math is correct: every constant, sign, and coefficient in eq_r, eq_g, the solved drift laws, the ratio/product factorization, and the n=5 collapse matches the notes' boxed equations and the appendix exactly, and both engines agree. I attacked the eq_r/eq_g coefficients (matched notes lines 99-110 and 154-163), the solve (re-derived by hand, matched), the ratio/product targets (independent forms, non-tautological, matched section 6), and the n=5 substitution (matched section 7) — all held. Two rigor gaps remain: the SymPy script does not actually verify the headline deliverable D1 (`if False` no-op at line 41; F1, medium), and the D5 channel-closure check is tautological in both engines (substituting a linear solution back into its own system; F2, low). Neither alters any result or downstream constant, so no `stop_cold`.

## Self-test notes

Trap 1 (variable independence): F1's proposed `a*sp.diff(sp.log(LW_law), a)` — `LW_law = πa√(1+r²)/(2√3)` does depend on `a`, so the derivative is `1/a` (nonzero), and `a·(1/a) − 1 = 0`; the `expect_zero` passes for the correct law and fails if the `a`-exponent changes, so it is not trivially-true. Trap 3 (trivial-case pre-check): substituting `r=1` into `LW_law` still gives degree-1 scaling in `a`, residual 0. Trap 5 (paper round-trip): the F1 fix introduces no new constant — it only differentiates the same `LW_law` already present and already matching the notes' L_W form (line 67); no new `paper_misalignment` created. F2's required change is offered as either a substantive δ⊥ construction or a downgrade-to-print, both safe.
