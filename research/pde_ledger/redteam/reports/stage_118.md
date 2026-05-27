---
unit_id: 118
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage118_parent_core_extraction.md]
  paper_appendix: present
---

# Audit unit 118 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_118.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage118_parent_core_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (line 1270 input)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage118_parent_core_sympy_audit.txt` (present, exit 0)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage118_parent_core_mathematica_audit.txt` (present, exit 0)

## What the paper claims

The stage card (line 16) states the stage "gives overlap formulas for `K_s, K_q, lambda, g_s, g_q` from a GNLS plus localized-Maxwell throat-core ansatz." The notes (authoritative, since the .tex is terse) supply five boxed deliverables:

1. Shell stiffness `K_s = 4πa²(H_w ℓ/3 + ℏ²/(15 m_ψ ρ_w ℓ))`; on the GNLS healing lock `ℓ = ℏ/(2 m_ψ c_sw)`, this reduces to `K_s = 3πa²ℏ²/(5 m_ψ ρ_w ℓ)`.
2. Mixed side-channel stiffness `K_q = (Z_q/μ₀)(π² c_s²/(4 L_W²))`.
3. Hybridization `λ = − q_* v_w0 I_sq` (note the **minus sign**, notes line 175-181). Uniform-core closure: `λ = −(8√2/3) q_* v_w0 a² ℓ √L_W` (notes line 195-198, explicitly boxed with a leading minus).
4. Shell mouth coupling `g_s = T_m · 4πa²ℓ/3`.
5. Mixed mouth coupling `g_q = (Z_q/μ₀) · π/(√2 L_W^(3/2))`.

Supporting integrals: `I_f = 1/3, I_g = 4/15` (tanh-wall moments); the D/N half-wave norm/derivative/endpoint data; and the bilinear `δ²H_sq = −q_* v_w0 ∫ varrho_s A_q · s q` from Madelung kinetic expansion (notes section 3, eqs (4)).

## What the script claims to verify

Both scripts share the same five-section structure:

- **I.** Tanh-wall moments via the substitution `t = tanh(y)`: assert `I_f = 1/3, I_g = 4/15`.
- **II.** D/N half-wave: assert `∫ χ² = 1`, `∫(χ')² = π²/(4 L_W²)`; print `∫ χ` and `χ'(0)`.
- **III.** Shell stiffness: assert `K_s` (formula built from `I_f, I_g`) equals the boxed closed form, and on healing lock matches `3πa²ℏ²/(5 m_ψ ρ_w ℓ)|_{ℓ=ℏ/(2 m_ψ c_sw)}`.
- **IV.** Bilinear `sq` coefficient: build `(1/2)m(ρ₀+s varrho_s)(v₀ − (q_*/m) q A_q)²`, extract `s q` coefficient, and assert it equals `−q_* v_w0 varrho_s A_q`.
- **V.** Effective mixed stiffness, mouth couplings, and **uniform-core closure of λ**: define `lam_uniform := qstar*v0*I_sq_uniform` (positive!) and assert it equals `+(8√2/3) q_* v_w0 a² ℓ √L_W` (positive!).

Section IV verifies the *correct* sign convention (`−q_* v_w0 varrho_s A_q`). Section V then drops the minus sign when assembling the uniform-core closure of λ.

## Paper ↔ script cross-check

| Paper-side deliverable | Script check | Status |
|---|---|---|
| Tanh-wall `I_f = 1/3, I_g = 4/15` | I.1, I.2 | match |
| D/N `∫ χ² = 1, ∫(χ')² = π²/(4L_W²)` | II.1, II.2 | match |
| D/N `χ'(0) = π/(√2 L_W^(3/2))` | implicitly via `g_q closed form` assertion (V) | match |
| `I_q = 2√(2 L_W)/π` | `I_q closed form` (V) | match |
| `K_s` general closed form | III.1 (`K_s closed form`) | match |
| `K_s` healing-lock reduction | III.2 (`healing-lock K_s`) | match |
| `K_q = (Z_q/μ₀) · π² c_s²/(4 L_W²)` | V (definition), no explicit `expect_zero` | partial (defined, not asserted) |
| `g_s = T_m · 4πa²ℓ/3` (via `J_s = 4πa²ℓ/3`) | V (`J_s closed form` and definition `g_s = T_m·J_s`) | partial (J_s asserted; `g_s = T_m·J_s` is by construction) |
| `g_q = (Z_q/μ₀) · π/(√2 L_W^(3/2))` | V (`g_q closed form`) | match |
| `δ²H_sq` bilinear coefficient `= −q_* v_w0 varrho_s A_q` | IV | match |
| `λ = −q_* v_w0 I_sq` (boxed eq, notes line 175-181) | (no explicit assertion of this from IV → V transition) | missing |
| `λ_uniform = −(8√2/3) q_* v_w0 a² ℓ √L_W` (boxed eq, notes line 195-198) | V (`lambda uniform closure`) but **with WRONG SIGN** (script asserts `+(8√2/3)…`) | mismatch |

Paper alignment is `misaligned` because the script's load-bearing closure assertion for λ contradicts the paper's boxed result by a sign.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 34 | `I_f == 1/3 and I_g == 4/15` (raise on mismatch) | tanh moments | yes |
| A2 | sympy | 49 | `expect_zero("D/N norm check", chi_norm - 1)` | D/N normalization | yes |
| A3 | sympy | 50 | `expect_zero("D/N stiffness check", chi_grad - π²/(4 L_W²))` | D/N stiffness | yes |
| A4 | sympy | 58 | `expect_zero("K_s closed form", K_s − K_s_expected)` | `K_s` general form | yes |
| A5 | sympy | 64 | `expect_zero("healing-lock K_s", K_s_heal − 3πa²ℏ²/(5 m_ψ ρ_w · (ℏ/(2 m_ψ c_sw))))` | `K_s` healing lock | yes (but RHS has ℓ substituted to a literal, not symbolic — see verdict note) |
| A6 | sympy | 77 | `expect_zero("bilinear sq coefficient", sq_coeff + qstar*varrho_s*v0*A_q)` | bilinear `−q_* v_w0 varrho_s A_q` | yes |
| A7 | sympy | 97 | `expect_zero("g_q closed form", g_q − (Z_q/μ₀)·π/(√2 L_W^(3/2)))` | `g_q` boxed form | yes |
| A8 | sympy | 98 | `expect_zero("J_s closed form", J_s − 4πa²ℓ/3)` | `J_s` (feeds `g_s` and `I_sq`) | yes |
| A9 | sympy | 99 | `expect_zero("I_q closed form", I_q − 2√(2 L_W)/π)` | `I_q` (feeds `λ`) | yes |
| A10 | sympy | 100-101 | `expect_zero("lambda uniform closure", lam_uniform − 8√2·qstar·v0·a²·ℓ·√L_W/3)` | `λ` uniform closure — but **wrong sign** | NO (verifies +λ; paper says −λ) |
| A1' … A10' | mathematica | 35-36, 53-54, 66, 71, 81, 105-108 | Mirror of A1–A10 | same | same (A10' has same sign error) |

Row A10/A10' fails alignment with the paper. Row A6/A6' verifies the **correct** sign for the bilinear — so the script is internally inconsistent: section IV finds `sq_coeff = −q_* v_w0 varrho_s A_q`, but section V drops the minus when integrating to closed form.

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** `target_mismatch`
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:88` and `:100-101`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:96` and `:108`

**What's wrong:**

The paper notes (`notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:174-182`) state:
> `λ = − q_* v_w0 𝓘_sq,  𝓘_sq := ∫ d⁴X varrho_s 𝒜_q.`

and the uniform-core closure (`notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:194-199`):
> `λ = −(8√2/3) q_* v_w0 a² ℓ √L_W.`

This minus sign is derived in the notes from the Madelung kinetic expansion:
> `δ²H_sq = − q_* v_w0 ∫ d⁴X varrho_s 𝒜_q · sq.`

The script's section IV correctly verifies the bilinear coefficient with the minus:
- SymPy line 74-77: `sq_coeff = −A_q qstar v0 varrho_s`, and the assertion is `sq_coeff + qstar*varrho_s*v0*A_q == 0` (i.e. confirms the minus).
- Mathematica line 78-81: identical structure.

But the section-V uniform-core closure drops the minus sign:
- SymPy line 88: `lam_uniform = sp.simplify(qstar * v0 * I_sq_uniform)` — **no minus**.
- SymPy line 100-101: `expect_zero("lambda uniform closure", lam_uniform − 8*sqrt(2)*qstar*v0*a²*ell*sqrt(L_W)/3)` — verifies `lam_uniform = +(8√2/3) qstar v0 a² ℓ √L_W`, contradicting the paper's `−(8√2/3) …`.
- Mathematica line 96: `lamUniform = FullSimplify[qStar*v0*iSqUniform, ...]` — same omission.
- Mathematica line 108: `expectZero["lambda uniform closure", lamUniform − (8*Sqrt[2]*qStar*v0*a²*ell*Sqrt[lW])/3]` — same wrong sign.

So the script verifies `λ = +(8√2/3) q_* v_w0 a² ℓ √L_W` while the paper claims `λ = −(8√2/3) q_* v_w0 a² ℓ √L_W`. Both engines share the error identically (consistent transliteration), so they cross-check each other on the wrong claim.

The script also never connects section IV (the `−q_* v_w0 varrho_s A_q` finding) to section V (the closure formula) via an explicit assertion that `λ = −q_* v_w0 I_sq`. The minus sign is simply dropped in the transition.

**Why this matters:**

Stage 118 feeds Stages 125-139 and the positive mouth-source branch (per stage card line 27). The hybridization parameter λ appears in the Stage-97/114 Schur complement only through `λ²` (so the squared quantity is sign-invariant — `K_s K_q + λ²` in denominators, `(K_s g_q − λ g_s)²` in numerators). The numerator term `(K_s g_q − λ g_s)²` does pick up the sign of λ via the cross term `−2 K_s g_q · λ g_s`. Downstream stages that read `λ` as a *signed* coupling — particularly anywhere the bilinear `λ s q` term itself enters a Hamiltonian or force law — would propagate the wrong sign. The notes' explicit derivation `δ²H_sq = − q_* v_w0 …` makes the sign physical (a static background flow `v_w0 > 0` and positive charge `q_* > 0` produce a *negative* coupling), so flipping it inverts the direction of the hybridization shift.

**Required change:**

This is `paper_misalignment` (subtype `target_mismatch`). Codex must **not** silently apply either direction. The directive's body for this finding is a `## Resolve before fix_loop` block — the user must decide whether (a) the paper is correct and the script should be flipped (the dominant reading, since section IV of the script itself derives the minus), or (b) the script's section V is correct and the notes are wrong, or (c) there is a convention disagreement (e.g., `v_w0` vs `−v_w0`, or `q_*` vs `−q_*`) elsewhere that makes both internally consistent.

The most likely correct resolution is (a): the script's section IV result `sq_coeff = −q_* v_w0 varrho_s A_q` is correct (matches the notes), and section V is a sign error introduced when going from coefficient to closure. In that case, `lam_uniform` should be redefined with a leading minus and the `expect_zero` target also given a leading minus.

**Verification:**

After user picks direction (a): `lam_uniform` should equal `−(8√2/3) qstar v0 a² ℓ √L_W`; the saved output's `lambda (uniform-core closure) = …` line should change sign, and the `lambda uniform closure` line should still print `0`.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:79-101`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:83-108`

**What's wrong:**

Of the five paper-side boxed deliverables in the notes ((1) `K_s` general & healing, (2) `K_q`, (3) `λ` general & uniform, (4) `g_s`, (5) `g_q`), two have no explicit `expect_zero`-style assertion against their boxed closed form:

- `K_q`: defined at line 82 as `(Zq/mu0) * (sp.pi**2*c_s**2/(4*L_W**2))` — this *is* the boxed closed form, so there is nothing to verify; the assertion is by construction. The corresponding paper claim is anchored in `chi_grad` (A3), so the `K_q` formula is supported indirectly. Acceptable, but worth a one-line `expect_zero("K_q closed form", K_q - (Zq/mu0) * chi_grad)` to demonstrate that the asserted D/N stiffness moment **is** what makes `K_q` take the boxed form.
- `g_s = T_m · 4πa²ℓ/3`: only the `J_s` factor is asserted (A8); `g_s = T_m * J_s` is a one-line definition, so by construction. Add `expect_zero("g_s closed form", g_s - Tm * 4*pi*a**2*ell/3)` to make this explicit.
- The general (non-uniform) hybridization `λ = − q_* v_w0 𝓘_sq` is never asserted in either form (general or via section IV's bilinear). The section IV check confirms the local bilinear coefficient; the section V definition of `lam_uniform` uses `I_sq_uniform`. The chain "bilinear coefficient → integrated → uniform closure" is implicit. An `expect_zero("lambda integrated from bilinear", lam_from_section_IV + qstar*v0*I_sq_uniform)` (using whatever sign the user picks in F1) would directly verify the connection.

**Why this matters:**

The script's load-bearing assertion list (10 `expect_zero` calls) covers the supporting integrals and `K_s` but leaves `K_q`, `g_s`, and the general-form `λ` as "by construction." A reader who trusts only the asserted lines cannot confirm from the script alone that the paper's boxed `K_q` and `g_s` formulas are the ones produced. This is fixable by three short additional assertions; the underlying derivations are all already present in the script.

**Required change:**

Add three explicit assertions:
1. `expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4 * L_W**2))` after line 90 in SymPy; mirror in Mathematica after line 98.
2. `expect_zero("g_s closed form", g_s - Tm * 4*sp.pi*a**2*ell/3)` after line 93 in SymPy; mirror in Mathematica after line 101.
3. `expect_zero("lambda general form", lam_uniform - <sign> * qstar * v0 * J_s * I_q)` after line 101 in SymPy (`<sign>` per F1 resolution); mirror in Mathematica after line 108.

**Verification:**

After Codex applies, the verifier confirms three new `expect_zero` lines appear in each script and three new `… = 0` lines in the saved output; script exits 0.

## Independent-derivation check (Mathematica)

The `.wl` script is a near-line-by-line transliteration of the `.py` script:

- Both use the identical `t = tanh(y)` substitution trick for `I_f, I_g`, computing `(1/4)(1−t²)` for `I_f`. Mathematica can integrate `sech^4(y)` natively (`Integrate[(1/2 Sech[y]^2)^2, {y, -Infinity, Infinity}]` yields `2/3 → … /3` after dividing by 2), but the script reuses SymPy's avoidance-substitution. (SymPy lines 30-31 vs. Mathematica lines 31-32.)
- The D/N section is variable-renamed but otherwise identical: same `χ = √(2/L_W) sin(πz/(2L_W))`, same four integrals/derivatives in the same order, same two assertions in the same order. (SymPy 39-50 vs. Mathematica 40-54.)
- Section III: identical formula assembly, identical `K_s_expected`, identical healing substitution. (SymPy 54-67 vs. Mathematica 58-71.)
- Section IV: identical `(1/2) m (ρ₀ + s varrho_s)(v₀ − (q_*/m) q A_q)²` expansion, identical `Coefficient[Coefficient[…, sAmp, 1], qAmp, 1]` extraction. (SymPy 71-77 vs. Mathematica 75-81.)
- Section V: identical assembly order and identical four assertions, including the identical sign error on `lam_uniform`. (SymPy 81-101 vs. Mathematica 85-108.)
- `expectZero` in Mathematica is a structural mirror of `expect_zero` in SymPy.

This is borderline `mathematica_transliteration`. However, the actual *computational work* (the integrals and the symbolic simplifications) is performed independently by each engine's kernel — the bug, if anything, is that both scripts share the **same algebraic choice path** (substitute `t = tanh y` even when Mathematica wouldn't need to). Since the integrals are evaluated by each engine's symbolic kernel against the substituted integrand (not by re-quoting a literal answer), I judge this stage *not* to merit a `mathematica_transliteration` finding on its own. The shared sign error in F1 is the more serious symptom of the lack of independence; once F1 is resolved on both sides, the engine cross-check still provides useful redundancy.

## Engine cross-check

Both engines reach identical final values (modulo trivial commutative/print-form differences):

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `I_f, I_g` | `1/3, 4/15` | `1/3, 4/15` |
| `∫ χ² dz` | `1` | `1` |
| `∫ (χ')² dz` | `π²/(4 L_W²)` | `Pi^2/(4 lW^2)` |
| `∫ χ dz` | `2√2 √L_W / π` | `(2 Sqrt[2] Sqrt[lW])/Pi` |
| `χ'(0)` | `√2 π / (2 L_W^(3/2))` | `Pi/(Sqrt[2] lW^(3/2))` (algebraically equal: `√2/2 = 1/√2`) |
| `K_s` | `4πa²(5 H_w ℓ² m_ψ ρ_w + ℏ²)/(15 ℓ m_ψ ρ_w)` | `(4 a² π (ℏ² + 5 ℓ² hWall mPsi rhoW))/(15 ℓ mPsi rhoW)` |
| `K_s` healing | `6πa² c_sw ℏ/(5 ρ_w)` | `(6 a² cSw ℏ π)/(5 rhoW)` |
| `sq coefficient` | `−A_q q_* v₀ varrho_s` | `−(aQ qStar v0 varrhoS)` |
| `λ (uniform)` | `(8√2 √L_W a² ℓ q_* v₀)/3` | `(8 Sqrt[2] a² ell Sqrt[lW] qStar v0)/3` |

All 10 `expect_zero` lines pass in both engines. `engines_agree: true`. Both engines have the **same sign error** on `λ`, so they cross-confirm the wrong claim.

Output mtimes: sympy script `1778522333` → output `1778525128` (output newer, fresh). Mathematica script `1778522213` → output `1778526571` (output newer, fresh). `outputs_fresh: true`.

## Verdict justification

The script's supporting machinery is correct: tanh-wall moments (`I_f = 1/3, I_g = 4/15`), D/N half-wave normalization and stiffness integral, `K_s` general and healing-lock closed forms, `g_q` and `J_s` closed forms, and the bilinear `sq` coefficient `−q_* v_w0 varrho_s A_q` all check out and are correctly anchored to the paper notes. Both engines agree numerically on every printed value.

However, the section-V uniform-core closure of `λ` drops the minus sign that section IV of the same script correctly derived, asserting `λ_uniform = +(8√2/3) q_* v_w0 a² ℓ √L_W` while the paper's boxed result is `λ_uniform = −(8√2/3) q_* v_w0 a² ℓ √L_W`. Both engines share this error identically. This is `paper_misalignment` (subtype `target_mismatch`, F1, high severity). Independently, three paper-side deliverables (`K_q` boxed form, `g_s` boxed form, `λ` general form) are present only as definitions, not as assertions, so the script's verification surface is narrower than the paper's claim surface (F2, medium, `insufficient_verification`).

Verdict: `findings`. `stop_cold: null` (the disagreement is repairable; the downstream impact via `λ²`-only entries in the Schur complement is sign-invariant for `(K_s K_q + λ²)` denominators, and the cross-term `−2 K_s g_q λ g_s` in numerators is the only signed appearance — see F1 "Why this matters"). The sign error is real but the orchestrator should halt for user resolution before either flipping the script or flipping the notes; this is the canonical paper_misalignment workflow.

Attacks tried: (a) re-checked the section-IV bilinear by hand: `(1/2) m (ρ₀ + s varrho_s)(v₀ − (q_*/m) q A_q)²` expanded gives `sq` term `s varrho_s · (1/2) m · 2 v₀ · (−q_*/m) q A_q = −q_* v₀ varrho_s A_q · s q`. Script's section IV assertion is correct, section V is inconsistent with section IV; (b) checked whether `v_w0` could be implicitly negative in the script's convention (e.g., infall velocity sign) — no such convention is documented in either the script or the notes; v_w0 is the static background mixed flow and is implicitly positive (the notes say `v_w0 ≠ 0` but don't restrict sign); (c) checked whether `I_sq_uniform = J_s · I_q` could be the source of a sign — both `J_s = 4πa²ℓ/3 > 0` and `I_q = 2√(2 L_W)/π > 0` (positive geometric integrals), so `I_sq_uniform > 0` and the asserted closure with a leading `+` is genuinely wrong; (d) checked whether the healing-lock check is tautological — the LHS substitutes `H_w, ℓ` into the K_s formula, and the RHS substitutes `ℓ` into `3πa²ℏ²/(5 m_ψ ρ_w ℓ)`. Both reduce to `6πa² c_sw ℏ/(5 ρ_w)` independently of K_s_expected, so the check verifies a real reduction. Not tautological; (e) checked whether `K_s_expected` on line 56 is hard-coded from outside derivation — no, `K_s_expected = 4πa²(H_w ℓ/3 + ℏ²/(15 m_ψ ρ_w ℓ))` is the boxed form from the notes; the actual K_s formula on line 55 is built from `I_f, I_g` which are asserted independently, so this is not `hardcoded_result`; (f) checked whether the symbol assumptions invalidate any step — `a, ell, H_w, hbar, mpsi, rho_w, c_sw, L_W, mu0, Zq, c_s, Tm` are all `positive=True, real=True`, which is consistent with their physical interpretation; `s, q, m, v0, qstar, varrho_s, A_q` are `real=True` (no positivity), correct for bilinear coefficient extraction. No `symbol_assumption_error`.

## Self-test notes

I checked: (1) variable independence is not relevant — no `sp.diff` calls in the load-bearing assertions of section V (the only derivative is `chi'(0)` in section II, which depends on `z` and is non-zero by direct computation); (2) symmetry/parity: the tanh-wall moments use `t ∈ [−1, 1]`; `(1/4)(1−t²)` is even (integral is `2 · (1/4) · 2/3 = 1/3` ✓), `t²(1−t²)` is even (integral is `2 · (1/3 − 1/5) = 4/15` ✓); (3) trivial-case pre-check: in the limit `q_* → 0` or `v_w0 → 0`, the bilinear vanishes and `λ → 0` in any sign convention — both wrong-sign and right-sign closure reduce to `0`, so the limit doesn't distinguish; the sign-distinguishing check would be a finite-`q_*, v_w0` numerical evaluation, which neither engine performs; (4) path specs not relevant (no missing scripts); (5) paper round-trip for F2: adding `expect_zero("K_q closed form", K_q - (Zq/mu0)*pi^2 c_s^2/(4 L_W^2))` uses only constants the paper card states; no new paper_misalignment introduced.
