---
unit_id: 118
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage118_parent_core_extraction.md]
  paper_appendix: present
---

# Audit unit 118 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_118.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage118_parent_core_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (present; stage included via `\input{stages/stage_118}` at line 1270; no additional summary row beyond the included card)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage118_parent_core_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage118_parent_core_mathematica_audit.txt`

## What the paper claims

The card (`Purpose`, `Derivation ledger`, and the quote at lines 15-17) states the stage gives explicit overlap formulas for the five reduced core parameters `(K_s, K_q, λ, g_s, g_q)` from one concrete GNLS + localized-Maxwell throat-core ansatz. The card itself is terse and routes detail to the notes; the card body quotes only "Gives overlap formulas for K_s, K_q, λ, g_s, g_q from a GNLS plus localized-Maxwell throat-core ansatz." The notes are authoritative on the deliverables and enumerate five boxed closed forms: `K_s = 4πa²(H_w ℓ/3 + ℏ²/(15 m_ψ ρ_w ℓ))` (and its healing-lock form `3πa²ℏ²/(5 m_ψ ρ_w ℓ)`), `K_q = (Z_q/μ0)·π²c_s²/(4 L_W²)`, `λ = -q_* v_{w0} I_sq` with uniform-closure value `-(8√2/3) q_* v_{w0} a² ℓ √L_W`, `g_s = T_m·4πa²ℓ/3`, and `g_q = (Z_q/μ0)·π/(√2 L_W^{3/2})`. The notes also pin the supporting moments/integrals: `I_f = 1/3`, `I_g = 4/15`, `∫χ²dz = 1`, `∫(χ')²dz = π²/(4L_W²)`, `I_q = ∫χ dz = 2√(2L_W)/π`, `χ'(0) = π/(√2 L_W^{3/2})`, and the bilinear `sq` term `δ²H_sq = -q_* v_{w0} ∫ρ_s A_q · sq`. Notably the card states NO explicit sign for λ; the notes carry three boxed minus signs (lines 175-181, 194-198, 277-281).

## What the script claims to verify

The SymPy script (docstring lines 2-7) verifies the first parent-action extraction of `(K_s, K_q, λ, g_s, g_q)`. Concretely it: (I) computes the tanh-wall shell moments `I_f, I_g` by the substitution `t = tanh y` and checks them against `1/3` and `4/15`; (II) computes the D/N half-wave norm, gradient norm, integral, and `χ'(0)` and checks the norm and the gradient (stiffness) integrals; (III) assembles `K_s` from `I_f, I_g` and checks it against the closed form, then substitutes the healing-lock relation and checks the locked form; (IV) expands the Madelung kinetic energy `(1/2)m(ρ0+s ρ_s)(v0-(q_*/m)q A_q)²`, extracts the `sq` coefficient, and checks it equals `-q_* ρ_s v0 A_q`; (V) assembles `K_q, g_q, J_s, g_s, I_q, I_sq, λ` and checks each against its closed form, including two independent forms of `λ`. The Mathematica `.wl` mirrors the same five sections with the same checks.

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| `I_f = 1/3` | py L30/34, wl L31/35 | match |
| `I_g = 4/15` | py L31/34, wl L32/36 | match |
| `∫χ²dz = 1` | py L41/49, wl L44/53 | match |
| `∫(χ')²dz = π²/(4L_W²)` | py L42/50, wl L45/54 | match |
| `I_q = 2√(2L_W)/π` | py L43/101, wl L46/109 | match |
| `χ'(0) = π/(√2 L_W^{3/2})` | py L44 (used in g_q L98), wl L47/106 | match |
| `K_s` closed form | py L55-58, wl L63-66 | match |
| `K_s` healing-lock form | py L61-67, wl L68-71 | match |
| `K_q = (Z_q/μ0)π²c_s²/(4L_W²)` | py L82/97, wl L90/105 | partial — closed-form check is X−X (tautology); see F1. Load-bearing integral `∫(χ')²` IS independently checked (L50). |
| `g_q = (Z_q/μ0)π/(√2 L_W^{3/2})` | py L83/98, wl L91/106 | match |
| `J_s = 4πa²ℓ/3` | py L84/99, wl L92/107 | match |
| `g_s = T_m·4πa²ℓ/3` | py L85/100, wl L93/108 | match |
| `I_sq (uniform) = 8√2√L_W a²ℓ/3` | py L87/(via λ), wl L95 | match |
| `λ = -q_* v0 I_sq` (bilinear) | py L74-77, wl L78-81 | match (MINUS, independently derived) |
| `λ (uniform) = -8√2√L_W a²ℓ q_* v0/3` | py L88/102-104, wl L96/110-111 | match (MINUS) |

Dominant pattern: aligned. One deliverable (`K_q`) is covered only by a self-referential restatement at the closed-form line, but its load-bearing integral is independently verified upstream in the same script — flagged low-severity F1, not a coverage gap.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 34 | `I_f-1/3==0 and I_g-4/15==0` (raise) | I_f, I_g moments | yes |
| A2 | sympy | 49 | `expect_zero(chi_norm - 1)` | D/N normalization | yes |
| A3 | sympy | 50 | `expect_zero(chi_grad - π²/(4L_W²))` | gradient norm → K_q | yes |
| A4 | sympy | 58 | `expect_zero(K_s - K_s_expected)` | K_s closed form | yes |
| A5 | sympy | 64-67 | `expect_zero(K_s_heal - 3πa²ℏ²/(5 m ρ ℓ))` | K_s healing lock | yes |
| A6 | sympy | 77 | `expect_zero(sq_coeff + q_* ρ_s v0 A_q)` | λ bilinear (sign) | yes |
| A7 | sympy | 97 | `expect_zero(K_q - (Zq/μ0)π²c_s²/(4L_W²))` | K_q closed form | no (X−X tautology; see F1) |
| A8 | sympy | 98 | `expect_zero(g_q - (Zq/μ0)π/(√2 L_W^{3/2}))` | g_q closed form | yes (chi'(0) derived L44) |
| A9 | sympy | 99 | `expect_zero(J_s - 4πa²ℓ/3)` | J_s | yes (uses I_f) |
| A10 | sympy | 100 | `expect_zero(g_s - T_m 4πa²ℓ/3)` | g_s | yes |
| A11 | sympy | 101 | `expect_zero(I_q - 2√(2L_W)/π)` | I_q | yes (chi_int derived L43) |
| A12 | sympy | 102-103 | `expect_zero(lam_uniform + 8√2 q_* v0 a²ℓ√L_W/3)` | λ uniform (sign) | yes |
| A13 | sympy | 104 | `expect_zero(lam_uniform + q_* v0 J_s I_q)` | λ bilinear↔uniform consistency | yes |
| B1-B13 | mathematica | 35-36,53-54,66,71,81,105-111 | mirror of A1-A13 | same claims | same (B7 mirrors A7 tautology) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:82,97`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:90,105`

**What's wrong:**
`K_q` is *defined* directly as its closed form and then checked against the identical literal:

SymPy:
```
82: K_q = sp.simplify((Zq/mu0) * (sp.pi**2*c_s**2/(4*L_W**2)))
...
97: expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4*L_W**2))
```
Mathematica:
```
90: kQ = FullSimplify[(zQ/mu0)*(Pi^2*cSound^2/(4*lW^2)), ...];
...
105: expectZero["K_q closed form", kQ - (zQ/mu0)*Pi^2*cSound^2/(4*lW^2)];
```
The subtraction is `X − X`; the check cannot fail regardless of the physics. Contrast the analogous deliverables `g_q`, `J_s`, `g_s`, `I_q`, all of which are assembled from *independently derived* upstream quantities (`chi_prime_0` at L44, `I_f` at L30, `chi_int` at L43) so that their closed-form checks genuinely exercise a derivation.

This reproduces the documented historical caveat exactly: the `K_q closed form` check is X−X. It is non-blocking because the load-bearing physical content — that `K_q = (Z_q/μ0)∫(χ')²dz` is *fixed* (the notes line 149: "not arbitrary: it is fixed by one localization norm and the D/N tube length") — IS independently verified upstream in the same script: A3 (`expect_zero(chi_grad - π²/(4L_W²))`, py L50 / wl L54) confirms `∫(χ')²dz = π²/(4L_W²)`, which is the only nontrivial step between the integral definition of `K_q` and its closed form. So the deliverable is genuinely covered; only the line-97/105 check itself is empty.

**Why this matters:**
As written, the `K_q closed form` assertion provides no verification and could mask a future regression if someone edits line 82's definition (the check would silently re-pass). Tying `K_q` to the already-verified gradient norm `chi_grad` makes the check load-bearing and self-consistent with how `K_q` is built in the notes.

**Required change:**
Define `K_q` from the verified gradient integral rather than from its target literal, so the closed-form check exercises a real reduction. In SymPy, replace line 82 with `K_q = sp.simplify((Zq/mu0) * c_s**2 * chi_grad)` (where `chi_grad` is the value computed and asserted at L42/L50, i.e. `π²/(4L_W²)`); leave line 97 unchanged so it now confirms `(Zq/mu0)·c_s²·chi_grad − (Zq/mu0)π²c_s²/(4L_W²) == 0`. In Mathematica, replace line 90 with `kQ = FullSimplify[(zQ/mu0)*cSound^2*chiGrad, Assumptions -> $Assumptions];` (using `chiGrad` from L45); leave line 105 unchanged. This makes the printed `K_q` value identical to today's output (the notes' `K_q = (Z_q/μ0)π²c_s²/(4L_W²)` with `c_s²` carried), so no downstream or reconciliation value changes — only the check becomes non-tautological.

**Verification:**
After the edit, py L97 / wl L105 still print `K_q closed form = 0` and PASS, but now `K_q` is `(Zq/mu0)*c_s**2*chi_grad`, so the residual genuinely depends on `chi_grad` having reduced to `π²/(4L_W²)`. Scripts must still exit 0 with all checks passing; printed `K_q` value (output line 35 / 44) is unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is a close structural parallel of the `.py`: identical five-section layout (I-V with the same banners), identical check names and order (`I_f - 1/3`, `D/N norm check`, `K_s closed form`, `bilinear sq coefficient`, `K_q/g_q/J_s/g_s/I_q closed form`, two λ checks), and the same variable choreography (just renamed: `chi`→`chi`, `chi_prime_0`→`chiPrime0`, `lam_uniform`→`lamUniform`). This is the expected pattern for this transliteration-watch band (105-175). However, the core computational steps are genuinely re-executed by Mathematica's own kernel, not echoed: the two engines independently evaluate `Integrate[(1/4)(1-t²),{t,-1,1}]`, `Integrate[chi²,...]`, `Integrate[D[chi,zTube]²,...]`, `Coefficient[Coefficient[expr,sAmp,1],qAmp,1]`, etc., and arrive at the same closed forms (e.g. SymPy `chi'(0) = sqrt(2)*pi/(2*L_W**(3/2))` vs Mathematica `chi'(0) = Pi/(Sqrt[2]*lW^(3/2))` — the same value in two different normal forms, which is a tell of independent simplification rather than copy). The bilinear `sq` coefficient is extracted by each engine's native polynomial-coefficient routine from a freshly expanded square. I judge this NOT a `mathematica_transliteration` finding: the parallel scaffolding is acceptable; the load-bearing algebra is independently computed by each CAS. The one shared weakness (the X−X `K_q closed form` check) is a logic defect present in both, captured by F1, not a transliteration defect.

## Engine cross-check

Both engines print identical results for every deliverable:

| value | SymPy output (line) | Mathematica output (line) |
|---|---|---|
| I_f | `1/3` (5) | `1/3` (4) |
| I_g | `4/15` (6) | `4/15` (4) |
| ∫χ²dz | `1` (11) | `1` (15) |
| ∫(χ')²dz | `pi**2/(4*L_W**2)` (12) | `Pi^2/(4*lW^2)` (17) |
| ∫χ dz | `2*sqrt(2)*sqrt(L_W)/pi` (13) | `(2*Sqrt[2]*Sqrt[lW])/Pi` (17) |
| χ'(0) | `sqrt(2)*pi/(2*L_W**(3/2))` (14) | `Pi/(Sqrt[2]*lW^(3/2))` (18) |
| K_s | `4*pi*a**2*(5*H_w*ell**2*mpsi*rho_w+hbar**2)/(15*ell*mpsi*rho_w)` (21) | `(4*a^2*Pi*(hbar^2+5*ell^2*hWall*mPsi*rhoW))/(15*ell*mPsi*rhoW)` (27) |
| K_s healing | `6*pi*a**2*c_sw*hbar/(5*rho_w)` (23) | `(6*a^2*cSw*hbar*Pi)/(5*rhoW)` (31) |
| sq coeff | `-A_q*qstar*v0*varrho_s` (29) | `-(aQ*qStar*v0*varrhoS)` (37) |
| K_q | `pi**2*Zq*c_s**2/(4*L_W**2*mu0)` (35) | `(cSound^2*Pi^2*zQ)/(4*lW^2*mu0)` (44) |
| g_q | `sqrt(2)*pi*Zq/(2*L_W**(3/2)*mu0)` (36) | `(Pi*zQ)/(Sqrt[2]*lW^(3/2)*mu0)` (45) |
| J_s | `4*pi*a**2*ell/3` (37) | `(4*a^2*ell*Pi)/3` (47) |
| g_s | `4*pi*Tm*a**2*ell/3` (38) | `(4*a^2*ell*Pi*tM)/3` (47) |
| I_sq | `8*sqrt(2)*sqrt(L_W)*a**2*ell/3` (39) | `(8*Sqrt[2]*a^2*ell*Sqrt[lW])/3` (48) |
| λ uniform | `-8*sqrt(2)*sqrt(L_W)*a**2*ell*qstar*v0/3` (40) | `(-8*Sqrt[2]*a^2*ell*Sqrt[lW]*qStar*v0)/3` (49) |

All thirteen `expect_zero`/`expectZero` residuals are 0 in both engines; every PASS line present in the Mathematica transcript. Engines agree. No `engine_disagreement`.

## λ-sign consistency finding (re-confirmed fresh)

The MINUS sign is consistent everywhere λ appears, and it is *independently derived*, not pasted in:

1. **Bilinear expansion (the source of the sign).** py L74 expands `(1/2)m(ρ0+s ρ_s)(v0-(q_*/m)q A_q)²`. The `sq` cross term is `(1/2)m·s ρ_s · 2 v0 · (-(q_*/m) q A_q) = -q_* ρ_s v0 A_q · sq`. The `coeff(s,1).coeff(q,1)` extraction yields `-A_q*qstar*v0*varrho_s` (output L29). The assertion L77 confirms `sq_coeff = -q_* ρ_s v0 A_q`. The minus comes from the cross term of the SQUARE, so it is forced by the algebra, not assumed. Matches notes line 168-172/175.
2. **Uniform-closure definition.** py L88 `lam_uniform = -qstar*v0*I_sq_uniform`; printed `-8√2√L_W a²ℓ q_* v0/3` (output L40); asserted at L102-103. Matches notes boxed line 194-198.
3. **Bilinear↔uniform cross-check.** py L104 asserts `lam_uniform = -q_* v0 J_s I_q`, tying the Section V value back to the Section IV bilinear sign. Consistent MINUS.
4. **Mathematica mirror.** wl L78-81, L96, L110-111 reproduce all three with the same MINUS.
5. **Card vs notes.** The card states no λ sign (it only references "overlap formulas"); the notes carry three boxed minuses (lines 175, 194, 280). No contradiction: the script's MINUS agrees with the notes, and the card is silent — so there is NO residual `paper_misalignment` on the λ sign. The historical λ-sign misalignment that was resolved to MINUS is fully consistent across card (silent), notes (3× minus), script Section IV (derived minus), and Section V (minus). The prompt's note that downstream stage-123 uses an un-squared λ is consistent with this stage exporting a signed λ; nothing in stage 118 squares λ or drops the sign. λ-sign: CONSISTENT, no finding.

## Value Reconciliation (pass-2 augmentation)

Authoritative record: script source + committed outputs (both fresh, mtime 2026-05-29 > script mtime 2026-05-27). Outputs not stale.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| I_f = 1/3 | py L30 / out L5; wl L31 / out L4 | notes:84 | MATCH |
| I_g = 4/15 | py L31 / out L6; wl L32 / out L4 | notes:86 | MATCH |
| ∫χ²dz = 1 | py L41 / out L11; wl L44 / out L15 | notes:136 | MATCH |
| ∫(χ')²dz = π²/(4L_W²) | py L42 / out L12; wl L45 / out L17 | notes:138 | MATCH |
| I_q = ∫χ dz = 2√(2L_W)/π | py L43 / out L13; wl L46 / out L17 | notes:190 | MATCH (2√2√L_W/π ≡ 2√(2L_W)/π) |
| χ'(0) = π/(√2 L_W^{3/2}) | py L44 / out L14; wl L47 / out L18 | notes:231 | MATCH (√2/2 ≡ 1/√2) |
| K_s = 4πa²(H_w ℓ/3+ℏ²/(15 m_ψ ρ_w ℓ)) | py L55 / out L21; wl L63 / out L27 | notes:91-99, 259-267 | MATCH (algebraically equal rearrangement) |
| K_s healing-lock = 3πa²ℏ²/(5 m_ψ ρ_w ℓ) → 6πa²ℏ c_sw/(5 ρ_w) | py L63 / out L23; wl L70 / out L31 | notes:107-112 | MATCH (= notes form under ℓ=ℏ/(2m_ψ c_sw)) |
| K_q = (Z_q/μ0)π²c_s²/(4L_W²) | py L82 / out L35; wl L90 / out L44 | notes:142-147, 271 | MATCH (value); check is tautological → F1 |
| g_q = (Z_q/μ0)π/(√2 L_W^{3/2}) | py L83 / out L36; wl L91 / out L45 | notes:235, 287 | MATCH |
| J_s = 4πa²ℓ/3 | py L84 / out L37; wl L92 / out L47 | notes:188 | MATCH |
| g_s = T_m·4πa²ℓ/3 | py L85 / out L38; wl L93 / out L47 | notes:217, 285 | MATCH |
| I_sq (uniform) = 8√2√L_W a²ℓ/3 | py L87 / out L39; wl L95 / out L48 | notes:186-190 (J_s·I_q) | MATCH |
| λ (uniform) = -8√2√L_W a²ℓ q_* v0/3 | py L88 / out L40; wl L96 / out L49 | notes:194-198, 277-281 | MATCH (MINUS) |
| sq coefficient = -A_q q_* v0 ρ_s | py L75 / out L29; wl L79 / out L37 | notes:168-172 | MATCH (MINUS) |

INTERNAL scaffolding (no prose carrier expected, no finding): `chi_s` profile (py L27), the residual `expect_zero`/`expectZero` "= 0" diagnostics, `K_s_expected`/`kSExpected` reference literals, `healing_sub` substitution dict, the intermediate `expr` polynomial, and all PASS flags.

reconciliation: complete; 15 deliverable values checked, 0 misaligned.

## Verdict justification

Verdict: findings (one low-severity, non-blocking `tautological_check`). The paper's five-parameter deliverable set is faithfully and independently verified across both engines, with every emitted value reconciling to the notes (15/15). I attacked the λ sign hardest given the documented history: it holds — the MINUS is forced by the cross term of the expanded Madelung square in Section IV, re-used in Section V, cross-checked between the two λ forms (L104), mirrored in Mathematica, and consistent with the notes' three boxed minuses; the card is silent so there is no live misalignment. I attacked the closed-form checks for X−X tautologies: only `K_q` (py L97 / wl L105) is genuinely tautological — the others (`g_q`, `J_s`, `g_s`, `I_q`, `K_s`) assemble from independently derived moments/integrals and would fail on a wrong moment. The `K_q` tautology is non-blocking because its load-bearing integral `∫(χ')²dz = π²/(4L_W²)` is independently asserted upstream (A3, py L50 / wl L54), so the deliverable is covered; F1 just tightens the check to depend on that integral. No symbol-assumption errors (positivity on `a, ℓ, L_W, μ0, Z_q, c_s, T_m` is justified by the physical setup; `qstar`/`v0` left unrestricted real, correct for a sign-carrying coupling). No hardcoded results (every literal is a derived target). Outputs fresh; engines agree. Not a transliteration. No stop-cold condition.

## Self-test notes

Checked: (1) variable independence — the only derivatives are `sp.diff(chi,z)` / `D[chi,zTube]` where `chi` genuinely depends on `z`/`zTube`, so `chi'(0)` and `∫(χ')²` are non-trivial (not identically zero). (2) Symmetry/parity — the `I_f, I_g` integrals are done via the `t=tanh y` substitution on `[-1,1]` with even integrands `(1-t²)` and `t²(1-t²)`, both even → nonzero, matching `1/3` and `4/15`; correct. (3) Trivial-case — the F1 fix `K_q = (Zq/mu0)*c_s²*chi_grad` with `chi_grad → π²/(4L_W²)` reproduces today's printed `K_q` exactly and the residual against the literal reduces to 0, so the tightened check passes for the right reason and changes no downstream value. (4) F1 paper round-trip — the fix introduces no new constant; `c_s²` already carries the printed K_q value, so no new `paper_misalignment`.
