---
unit_id: 118
batch: IV.3
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 118

Apply each non-`paper_misalignment` finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F1), do nothing — the orchestrator is holding for user resolution. Do not edit `paper/`, `notes/`, or scripts to "fix" the paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-`paper_misalignment` finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch `paper.tex`, `notes/`, or any prose documents.

## F1 — paper_misalignment

**Subtype:** `target_mismatch`

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:174-182` quote:
  > "λ = − q_* v_{w0} 𝓘_{sq}, 𝓘_{sq} := ∫ d⁴X varrho_s 𝒜_q."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:194-199` quote:
  > "λ = − (8√2 / 3) q_* v_{w0} a² ℓ √L_W."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:168-173` quote (the bilinear derivation that fixes the sign):
  > "δ²H_{sq} = − q_* v_{w0} ∫ d⁴X varrho_s 𝒜_q · s q."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:88` quote:
  > `lam_uniform = sp.simplify(qstar * v0 * I_sq_uniform)`  (no minus sign)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:100-101` quote:
  > `expect_zero("lambda uniform closure", lam_uniform - 8*sp.sqrt(2)*qstar*v0*a**2*ell*sp.sqrt(L_W)/3)`  (verifies `λ_uniform = + (8√2/3) q_* v₀ a² ℓ √L_W`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:96` quote:
  > `lamUniform = FullSimplify[qStar*v0*iSqUniform, Assumptions -> $Assumptions];`  (no minus)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:108` quote:
  > `expectZero["lambda uniform closure", lamUniform - (8*Sqrt[2]*qStar*v0*a^2*ell*Sqrt[lW])/3];`  (verifies `+ (8√2/3) …`)

Note: the script's own section IV (SymPy line 77; Mathematica line 81) **correctly** finds the bilinear `sq` coefficient as `− q_* v₀ varrho_s A_q`. So section V is internally inconsistent with section IV in addition to disagreeing with the notes.

## Resolve before fix_loop

The paper notes and the script's own section-IV bilinear coefficient both have `λ = − q_* v_{w0} 𝓘_{sq}` (minus sign); the script's section-V closure asserts `λ_uniform = + (8√2/3) q_* v_{w0} a² ℓ √L_W` (plus sign). Which sign is correct?

Possible directions (the user picks one):
- (a) **Notes are correct** (most likely, since the script's own section IV agrees with the notes): change SymPy line 88 to `lam_uniform = sp.simplify(-qstar * v0 * I_sq_uniform)` and SymPy line 100-101 target to `lam_uniform + 8*sp.sqrt(2)*qstar*v0*a**2*ell*sp.sqrt(L_W)/3`; mirror in Mathematica lines 96 and 108. Re-run sympy+mathematica.
- (b) **Script section V is correct** and the notes (`notes/stages/moving_throat_pde_stage118_parent_core_extraction.md:174-182, 194-199`) and the script's section IV (lines 71-77) are wrong — flip the notes' minus signs and re-derive the bilinear. (Unlikely: the Madelung kinetic expansion derivation in the notes is explicit and matches the script's section IV.)
- (c) **Convention disagreement**: somewhere in the chain `v_w0` or `q_*` is implicitly negated, making both consistent under a hidden sign convention. (Unlikely: neither file documents such a convention.)

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:90-102` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:98-109`

**Issue:**

Three of the five paper-side boxed deliverables (`K_q`, `g_s`, and `λ` in its general bilinear form `λ = − q_* v_w0 𝓘_sq`) appear in the script only as definitions or are derived from a separate section without an `expect_zero` connecting the definition to the boxed closed form. The script's coverage of the paper claim is narrower than it should be. Adding three `expect_zero` lines makes the verification surface match the paper's claim surface.

**This finding is DEFERRED until F1 is resolved.** The "lambda general form" assertion (item 3 below) depends on the sign chosen in F1. Codex should **not** apply F2 until the orchestrator confirms F1's resolution. If the user picks F1 direction (a), apply F2 below as written. If the user picks (b) or (c), the F2 instructions are updated accordingly in a follow-up directive.

**Required change (assuming F1 direction (a) — minus sign):**

1. **K_q closed-form assertion.** After SymPy line 95 (after `print("lambda (uniform-core closure) = ", lam_uniform)`), insert before line 97:
   ```python
   expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4*L_W**2))
   ```
   Mirror in Mathematica after line 103 (before line 105):
   ```mathematica
   expectZero["K_q closed form", kQ - (zQ/mu0)*Pi^2*cSound^2/(4*lW^2)];
   ```

2. **g_s closed-form assertion.** Insert immediately after the K_q assertion in SymPy:
   ```python
   expect_zero("g_s closed form", g_s - Tm * 4*sp.pi*a**2*ell/3)
   ```
   Mirror in Mathematica immediately after the K_q assertion:
   ```mathematica
   expectZero["g_s closed form", gS - tM*4*Pi*a^2*ell/3];
   ```

3. **lambda general-form assertion** (sign per F1 direction (a)). After SymPy line 101 (after the existing `lam_uniform - …` assertion, which under F1(a) will have been flipped to `lam_uniform + …`):
   ```python
   expect_zero("lambda from bilinear", lam_uniform + qstar * v0 * J_s * I_q)
   ```
   Mirror in Mathematica after the (flipped) line-108 assertion:
   ```mathematica
   expectZero["lambda from bilinear", lamUniform + qStar*v0*jS*iQ];
   ```
   (Under F1(a), `lam_uniform = − qstar v0 I_sq_uniform = − qstar v0 J_s I_q`, so `lam_uniform + qstar v0 J_s I_q = 0`.)

**Claim manifest:**

M1 — `K_q = (Z_q / μ₀) · π² c_s² / (4 L_W²)`  (notes line 142-147)
M2 — `g_s = T_m · 4πa²ℓ/3`  (notes line 216-219)
M3 — `λ = − q_* v_{w0} · J_s · I_q` in uniform-core closure  (notes line 175-198; sign per F1 resolution)

Each new assertion exercises exactly the boxed identity from the notes it cites.

**Verification command:**

After Codex applies (post-F1 resolution), the verifier will run `redteam exec-sympy 118` and `redteam exec-mathematica 118` and confirm the three new `expect_zero` / `expectZero` lines appear in each script and three new `… = 0` and `PASS:` lines appear in each saved output; both scripts exit 0.
