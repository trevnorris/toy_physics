---
unit_id: 094
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: false
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage094_isotropic_geometry_decoupling.md]
  paper_appendix: present
---

# Audit unit 094 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_094.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage094_isotropic_geometry_decoupling.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_094}` line at L1222 — no separate appendix row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.txt`

## What the paper claims

Stage 094 is the "Isotropic Geometry-Decoupling Theorem." The stage card body quotes (verbatim): "Orthogonality of l=0 geometry and grouped l=2 wall modes gives K_{g,2}=K_{g,4}=0 on the isotropic branch." The notes elaborate: on an isotropic moving-throat branch, the scalar/geometry l=0 lane and the grouped real l=2 lanes are exactly orthogonal in the quadratic wall theory (notes §1–2), so the dynamic l=0 geometry lane contributes no bilinear cross term to the grouped P2 conservative quadrupole module at linear order. Consequently, ε_2 = Ω_Q² K_{g,2}/K_pole = 0 and ε_4 = Ω_Q⁴ K_{g,4}/K_pole = 0, and the 3/4 + 1/4 split is recovered with c_pole = 1/4, c_geom = 3/4 (notes §3). The card's `\stagefield{Checks}` enumerates three items, the first of which is: "Check the static limit ε_2=ε_4=0 returns c_pole=1/4."

## What the script claims to verify

Both scripts compute on real spherical harmonics Y00, Y20, Y21c, Y21s, Y22c, Y22s on S². They (i) verify Y00 and each Y2A are unit-normalized; (ii) verify (-Δ_S²)Y2A = 6 Y2A; (iii) compute and zero-assert the three angular cross integrals ⟨Y00|Y2A⟩, ⟨∇Y00·∇Y2A⟩, ⟨Y00|(-Δ_S²)Y2A⟩; and (iv) form a "Generic isotropic cross coefficient C_{0,A}" out of symbolic coefficients (mu, Tw/tw, TOm/tOm, K/kPot) and the angular integrals, and assert it is identically zero for each A. There is no script-side computation of ε_2, ε_4, K_pole, K_{g,2}, K_{g,4}, c_pole, or the static-limit reduction to c_pole = 1/4.

## Paper ↔ script cross-check

| Paper-side deliverable | Script coverage |
|---|---|
| l=0 vs each l=2 wall-mode orthogonality (mass, gradient, Laplacian) | `match` — sympy L40–42 and `.wl` L63–66 zero-assert all three integrals for every Y2A |
| Block-diagonal bilinear: cross coefficient between g(w,t) Y00 and q_A(w,t) Y2A vanishes for the action μ(∂_t η)² + T_w(∂_w η)² + T_Ω η(-Δ)η + K η² | `mismatch` — SymPy L52 forms `mu*I_mass - Tw*I_mass - TOm*I_lap - K*I_pot` (physically correct: T_w-cross involves angular ∫Y00·Y2A because ∂_w factors out); Mathematica L60 forms `mu*overlap - tw*gradCross - tOm*lapCross - kPot*overlap` (uses `gradCross` instead of `overlap` for the tw term, which represents a different physical coupling). Both zero only because all three angular integrals vanish, but the engines are not computing the same expression. |
| K_{g,2} = 0, K_{g,4} = 0 named explicitly | `partial` — script verifies the underlying orthogonality but never instantiates K_{g,2}, K_{g,4} as symbolic objects nor labels its zero-result as those two scalars. |
| ε_2 = ε_4 = 0; static limit gives c_pole = 1/4 (paper Check #1) | `missing` — script has no ε_2, ε_4, K_pole, c_pole symbols and no static-limit check. |

`paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `assert domega(simplify(Y00*Y)) == 0` for each Y2A | l=0/l=2 mass-overlap orthogonality | yes |
| A2 | sympy | 41 | `assert domega(grad_cross) == 0` for each Y2A | angular gradient orthogonality | yes |
| A3 | sympy | 42 | `assert domega(simplify(Y00*(-lap_s2(Y)))) == 0` for each Y2A | Laplacian cross orthogonality | yes |
| A4 | sympy | 54 | `assert simplify(Ccross) == 0` for each Y2A | bilinear cross coefficient (using `mu*I_mass - Tw*I_mass - TOm*I_lap - K*I_pot`) | yes (correctly tracks the wall action terms) |
| A5 | mathematica | 50 | `expectZero["Y00 normalization - 1", dOmega[y00^2] - 1]` | Y00 normalization | yes |
| A6 | mathematica | 62 | `expectZero["Y" <> label <> " normalization - 1", norm - 1]` | Y2A normalization | yes |
| A7 | mathematica | 63 | `expectZero["<Y00|Y" <> label <> ">", overlap]` | mass-overlap orthogonality | yes |
| A8 | mathematica | 64 | `expectZero["(-Delta)Y" ... " - 6 Y" ...]` | Laplacian eigenvalue | yes |
| A9 | mathematica | 65 | `expectZero["<grad Y00 . grad Y" ...]` | gradient-cross orthogonality | yes |
| A10 | mathematica | 66 | `expectZero["<Y00|(-Delta)Y" ...]` | Laplacian cross orthogonality | yes |
| A11 | mathematica | 67 | `expectZero["Generic isotropic cross coefficient C_0," ..., cCross]` (uses `mu*overlap - tw*gradCross - tOm*lapCross - kPot*overlap`) | bilinear cross coefficient | partial — does not match the wall-action structure used in SymPy and in notes §1 |

## Findings

### F1 — engine_disagreement

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py:52`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:60`

**What's wrong:**
The two engines purport to compute the same object ("Generic isotropic cross coefficient C_{0,A}") but actually compute different symbolic expressions.

SymPy (line 52):
```
Ccross = sp.simplify(mu*I_mass - Tw*I_mass - TOm*I_lap - K*I_pot)
```
Here `Tw` (coefficient of (∂_w η)² in the wall action of notes §1) multiplies `I_mass = ⟨Y00|Y2A⟩`, which is physically correct: in the bilinear cross term `2 T_w (∂_w g)(∂_w q_A) ∫ Y00 Y2A dΩ`, the angular integral is `I_mass`.

Mathematica (line 60):
```
cCross = FullSimplify[mu*overlap - tw*gradCross - tOm*lapCross - kPot*overlap, ...]
```
Here `tw` multiplies `gradCross = ⟨∇_S² Y00 · ∇_S² Y2A⟩` rather than `overlap`. That is a different physical coupling (it represents an angular gradient cross term, which would arise from an operator like `T_Ω' (∇_S² η)·(∇_S² η)` but NOT from `T_w (∂_w η)²`).

Both expressions evaluate to 0 because each individual angular integral (overlap, gradCross, lapCross) vanishes on its own. So the assertions both pass — but the agreement is coincidental at the action-structure level, not enforced by the algebra. Were any one of the three orthogonalities to fail, the two engines would disagree on which `T*` coefficient should multiply the surviving residual, and the cross-engine cross-check would be illusory.

**Why this matters:**
The point of having two independent engines is to detect transliteration errors and physical-setup mistakes. Here the SymPy expression is faithful to the wall action stated in notes §1; the Mathematica expression is not. The mismatch is masked by the fact that the orthogonality theorem holds component-wise, but if a future revision modifies the action (e.g., adds a small anisotropy) the SymPy script will correctly track which coupling survives while the Mathematica script will not. It also undermines the claim (implicit in having two engines) that two independent re-derivations agree.

**Required change:**
Change the Mathematica `cCross` definition on line 60 to mirror the action structure: the `tw` (T_w) coefficient should multiply `overlap`, not `gradCross`. Replace
```
cCross = FullSimplify[mu*overlap - tw*gradCross - tOm*lapCross - kPot*overlap, Assumptions -> $Assumptions];
```
with
```
cCross = FullSimplify[mu*overlap - tw*overlap - tOm*lapCross - kPot*overlap, Assumptions -> $Assumptions];
```
(SymPy line 52 is already in the correct form; do not change it. The `gradCross` value is still independently asserted to vanish on line 65, so removing it from `cCross` does not lose coverage of the gradient-orthogonality check.)

**Verification:**
After the edit, line 60 of the .wl reads `cCross = FullSimplify[mu*overlap - tw*overlap - tOm*lapCross - kPot*overlap, ...]`. Mathematica re-run still PASSes all six `expectZero` lines per Y2A (the change is from `tw*0` to `tw*0`, still zero). The transcript line "Generic isotropic cross coefficient C_0,*" still prints "= 0" / PASS.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py:1-57`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:1-75`

**What's wrong:**
The paper card's `\stagefield{Checks}` enumerates as Check #1: "Check the static limit ε_2=ε_4=0 returns c_pole=1/4." Neither script instantiates ε_2, ε_4, K_pole, K_{g,2}, K_{g,4}, or c_pole as symbolic objects, and neither script performs the static-limit substitution showing c_pole reduces to 1/4. The notes §3 give the algebraic chain:

> `eps_2 = Omega_Q^2 K_(g,2) / K_pole = 0`, `eps_4 = Omega_Q^4 K_(g,4) / K_pole = 0` ... `c_pole = 1/4`, `c_geom = 3/4`.

The scripts verify the upstream premise (l=0/l=2 orthogonality, and therefore the bilinear cross coefficient vanishes) but do not name K_{g,2}, K_{g,4}, ε_2, ε_4, or c_pole, nor do they exercise the static-limit reduction the card explicitly lists as a check.

The script's coverage is a strict subset of what the card and notes describe: the orthogonality theorem (necessary condition) is proved, but the named consequences (sufficient form, with the specific labels and the c_pole = 1/4 reduction) are not.

**Why this matters:**
The auditor cannot, from these scripts alone, verify the load-bearing constant `c_pole = 1/4` that the card claims is the static-limit output. If a later stage cites c_pole = 1/4 as a derived constant from stage 094, the citation traces back to scripts that do not in fact derive that number — they derive only the orthogonality that the notes argue implies it.

**Required change:**
Add a short symbolic block at the end of each script that:
1. Defines `K_g2`, `K_g4`, `K_pole`, `Omega_Q` as symbols (or assigns `K_g2 = 0`, `K_g4 = 0` as the orthogonality result).
2. Defines `eps_2 = Omega_Q**2 * K_g2 / K_pole` and `eps_4 = Omega_Q**4 * K_g4 / K_pole`.
3. Asserts each is 0 under the orthogonality result.
4. Defines `c_pole = sp.Rational(1, 4)` and `c_geom = sp.Rational(3, 4)` and prints them (or, more strongly, plugs ε_2 = ε_4 = 0 into whatever symbolic obstruction formula is referenced from upstream and confirms `c_pole = 1/4`). At minimum, the script should *name* the constants the card claims it produces.

SymPy patch sketch (append after line 56):
```python
# Static-limit check: with K_(g,2) = K_(g,4) = 0 from orthogonality,
# the contamination numbers vanish and the 3/4 + 1/4 split is recovered.
Omega_Q, K_pole = sp.symbols('Omega_Q K_pole', positive=True)
K_g2 = sp.Integer(0)   # established by orthogonality above
K_g4 = sp.Integer(0)   # established by orthogonality above
eps_2 = Omega_Q**2 * K_g2 / K_pole
eps_4 = Omega_Q**4 * K_g4 / K_pole
assert sp.simplify(eps_2) == 0
assert sp.simplify(eps_4) == 0
c_pole = sp.Rational(1, 4)
c_geom = sp.Rational(3, 4)
assert c_pole + c_geom == 1
print('eps_2 =', eps_2, '; eps_4 =', eps_4, '; c_pole =', c_pole, '; c_geom =', c_geom)
```

Mirror an analogous block in the `.wl` (using `0` for `Kg2`, `Kg4`, defining `eps2`, `eps4`, `cPole`, `cGeom` and `expectZero`-ing the contamination numbers; `expectZero["c_pole + c_geom - 1", cPole + cGeom - 1]`).

**Verification:**
After the edit, sympy transcript shows two new `assert` lines pass (no AssertionError) and prints a line containing `eps_2 = 0 ; eps_4 = 0 ; c_pole = 1/4 ; c_geom = 3/4`. Mathematica transcript shows matching new `PASS:` lines for the eps2/eps4 zero check and for `c_pole + c_geom - 1`.

## Independent-derivation check (Mathematica)

Strong transliteration suspicion. The `.wl` script mirrors the `.py` script line-by-line:
- Same harmonic definitions (sympy L24–29 vs. wl L41–46).
- Same loop over the same five Y2A keys (sympy L33 vs. wl L52).
- Same six per-mode quantities computed in the same order: normalization, overlap, lap-residual, gradCross, lapCross, cCross (sympy L34–53 vs. wl L55–67).

However, the two engines do NOT compute the same `cCross` (F1 above); the SymPy form matches notes §1 and the Mathematica form does not. This is unusual: transliteration normally produces identical algebra. Here the algebra is *similar but inequivalent* — the Mathematica author seems to have transliterated the *list of available angular integrals* but recombined them slightly differently when forming `cCross`. I record this as a transliteration-shaped script with a substantive algebraic divergence (F1). I do not file a separate `mathematica_transliteration` finding because the F1 divergence shows the two scripts are not identical at the load-bearing line; flagging both would double-count.

## Engine cross-check

All scalar orthogonality assertions agree (both engines: overlap = 0, gradCross = 0, lapCross = 0, normalization = 1, eigenvalue 6 confirmed for all Y2A). The `cCross` final value is 0 in both engines, but the symbolic expression being asserted differs (see F1). Outputs are fresh (sympy .py mtime 2026-04-01; sympy .txt mtime 2026-05-11; .wl mtime 2026-05-11 11:56; .txt mtime 2026-05-11 13:03 — script ≤ output for both).

## Verdict justification

The unit's stated theorem (l=0 ⊥ l=2 in the isotropic quadratic wall theory) is solidly verified by both scripts via explicit S² integrals over real spherical harmonics, with non-tautological assertions exercising every angular integral the notes invoke. The verdict is `findings`, not `clean`, because (i) the two engines compute different symbolic forms for what they both label "Generic isotropic cross coefficient C_{0,A}" (F1), and (ii) the paper card's explicit Check #1 ("static limit ε_2 = ε_4 = 0 returns c_pole = 1/4") is not exercised by either script (F2). Neither finding compromises the underlying mathematics — the orthogonality theorem is genuinely proved — but both are real gaps between paper-side claim and script-side coverage. `stop_cold: null` because the fixes are local, mechanical, and do not propagate. Not a paper_misalignment because the notes and card do not assert anything the scripts contradict; the issue is missing coverage and inconsistent cCross definition.

## Self-test notes

I checked: (1) the substituted Mathematica `cCross` (with `tw*overlap` replacing `tw*gradCross`) still simplifies to 0 because `overlap = 0` for every Y2A — confirmed by tracing the SymPy run, which uses the same form. (2) The new SymPy block in F2 uses `K_g2 = sp.Integer(0)` as an assignment, so `eps_2 = Omega_Q**2 * 0 / K_pole = 0` simplifies cleanly without division-by-zero risk (K_pole is positive symbol, not zero). (3) The static-limit constants `c_pole = 1/4`, `c_geom = 3/4` are quoted verbatim from notes §3 — the paper round-trip lands on the same values the card claims. (4) Path specifications: F2 edits target the existing scripts in `scripts/` and `mathematica/`; no new files. (5) Symbol parity / domain: `Omega_Q`, `K_pole` declared positive; ε_n are well-defined real numbers.
