---
unit_id: 116
batch: IV.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md"]
  paper_appendix: present
---

# Audit unit 116 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_116.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage116_dn_mixed_tube_realization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows for stage_116 found at line 1266; MTDC-T8 anchor block lines 1175-1179)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_116.tex` line 15-17) states verbatim: "First D/N mixed half-wave fixes \(L_W=\pi a\sqrt{(1+r_c)/3}/2\) and the bare outgoing scale." The notes (authoritative, more detailed) lay out four concrete deliverables: (1) the D/N half-wave eigenvalue `k_W = π/(2L_W)`, `Ω_W = π c_s/(2L_W)`; (2) the bare even coefficient `κ_0 = (ω²/Ω_W²)/z² = 4L_W²/(π²a²)` derived from the eigenfunction problem `q''+k²q=0, q(0)=0, q'(L_W)=0`; (3) the compensation-selected tube length `L_W = πa/2 · sqrt((1+r_c)/3)` enforced by `κ_0 = (1+r_c)/3` with `r_c = λ²/(K_s K_q)`; (4) the bare outgoing form `D_W^bare(z) = (1+r_c)(1 - z²/3 - i z⁵/9) + O(z⁶)` and γ_0 = (1+r_c)/9; (5) after Stage-97 (1+r_c) renormalization, the canonical coefficients `κ_c = 1/3`, `γ_c = 1/9`.

## What the script claims to verify

The SymPy script defines `r_c = lam²/(K_s K_q)`, hardcodes `kappa0_from_tube = 4L_W²/(π²a²)`, solves the compensation equation `kappa0_from_tube = (1+r_c)/3` for `L_W` (printing the result, but with no assertion), then defines `kappa0_bare = (1+r_c)/3` and divides by `(1+r_c)` and asserts the quotient equals `1/3`. Same pattern for `γ_c`. Finally, it constructs `D_bare = (1+r_c)·(canonical-form)`, divides by `(1+r_c)`, and asserts the quotient equals the canonical form. The Mathematica script mirrors this structure, but adds one additional substantive assertion: `lWRequired - πa·sqrt((1+r_c)/3)/2 == 0` (the actual tube-length law of paper claim 3).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| D/N half-wave eigenvalue `k_W = π/(2L_W)` derived from `q''+k²q=0, q(0)=0, q'(L_W)=0` | (none — eigenvalue problem not posed) | missing |
| Bare coefficient `κ_0 = 4L_W²/(π²a²)` from the eigenvalue | Literal definition `kappa0_from_tube = 4L_W²/(π²a²)` (sympy:26; wl:34) | mismatch (asserted as literal, not derived) |
| Tube-length law `L_W = πa√((1+r_c)/3)/2` | SymPy: only prints `L_W_required`, no assertion. Mathematica: `expectZero["tube-length law", ...]` (wl:39) | partial (only Mathematica covers it) |
| `γ_0 = (1+r_c)/9` literally appears | Defined as literal `gamma0_bare = (1+r_c)/9` (sympy:36; wl:42) | extra/tautological scaffolding |
| Renormalization: `(1+r_c)·canonical` ÷ `(1+r_c)` = canonical | Both engines assert this via tautological constructions | match-but-tautological |
| Final values `κ_c = 1/3, γ_c = 1/9` | Both engines assert `kappa_c - 1/3 == 0`, `gamma_c - 1/9 == 0` | tautological (algebraic identity of construction) |

`paper_alignment: partial` — the central paper claim 3 (tube-length law) is exercised on the Mathematica side but not on the SymPy side. The two interior claims of the stage (κ_0 derivation from eigenfunction, then κ_c emerging from renormalization) are not substantively verified: κ_0 is hardcoded and κ_c follows tautologically.

## Assertion inventory

| #  | Script       | Line | Form                                                  | Exercises which paper claim?                  | Anchored to claim? |
|----|--------------|------|-------------------------------------------------------|-----------------------------------------------|--------------------|
| A1 | sympy        | 41   | `expect_zero("final kappa_c - 1/3", kappa_c - 1/3)`   | claim 5 (final κ_c=1/3)                       | no (tautological)  |
| A2 | sympy        | 42   | `expect_zero("final gamma_c - 1/9", gamma_c - 1/9)`   | claim 5 (final γ_c=1/9)                       | no (tautological)  |
| A3 | sympy        | 46-49| `D_final - (1 - z²/3 - i z⁵/9) == 0`                  | claim 4/5 (renormalization)                   | no (tautological)  |
| A4 | mathematica  | 39   | `lWRequired - πa·sqrt((1+r_c)/3)/2 == 0`              | claim 3 (tube-length law)                     | yes (substantive)  |
| A5 | mathematica  | 46   | `kappaC - 1/3 == 0`                                   | claim 5                                       | no (tautological)  |
| A6 | mathematica  | 47   | `gammaC - 1/9 == 0`                                   | claim 5                                       | no (tautological)  |
| A7 | mathematica  | 51   | `dFinal - (1 - z²/3 - i z⁵/9) == 0`                   | claim 4/5                                     | no (tautological)  |

Only A4 is non-tautological. Every other assertion follows from the algebraic structure the script itself imposes a few lines earlier.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:35-42`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:41-47`

**What's wrong:**
SymPy lines 35-42:
```
kappa0_bare = sp.simplify((1 + r_c) / 3)
gamma0_bare = sp.simplify((1 + r_c) / 9)
kappa_c = sp.simplify(kappa0_bare / (1 + r_c))
gamma_c = sp.simplify(gamma0_bare / (1 + r_c))
expect_zero("final kappa_c - 1/3", kappa_c - sp.Rational(1, 3))
expect_zero("final gamma_c - 1/9", gamma_c - sp.Rational(1, 9))
```
The script defines `kappa0_bare = (1+r_c)/3` as a literal value, then divides by `(1+r_c)`. By construction `kappa_c == 1/3` for any value of r_c, so the assertion cannot fail. Same for γ_c, and same in Mathematica (lines 41-47).

The notes (paper claim 5) describe the renormalization as: the bare coefficient `(1+r_c)/3` is what falls out of the κ_0 = 4L_W²/(π²a²) calculation combined with κ_0 = (1+r_c)/3, and the (1+r_c) gets stripped by Stage-97's denominator renormalization. The script bypasses every step of that physics and writes the bare value as `(1+r_c)/3` directly. The assertion just confirms `(x/3)/x = 1/3`.

**Why this matters:**
The paper claim is that combining the D/N geometric κ_0 with the upstream compensation constraint yields κ_c = 1/3 after renormalization. The script's assertion confirms only the algebraic identity `((1+r_c)/3)/(1+r_c) = 1/3`, which is true for any rational expression. A bug in the upstream chain (e.g., if κ_0 from compensation were actually (1+2r_c)/3) would not be caught.

**Required change:**
Anchor `kappa0_bare` to the geometric derivation. Replace the literal `kappa0_bare = (1+r_c)/3` with the explicit solve: compute `L_W` from `kappa0_from_tube = (1+r_c)/3`, then substitute back to get `kappa0_bare = 4*L_W_required**2/(pi**2*a**2)`, simplify, and assert this equals `(1+r_c)/3` (i.e., the round-trip through L_W reproduces the input). Same for the Mathematica script. The renormalization check then becomes: assert that `kappa0_bare/(1+r_c)` simplifies to `1/3` where kappa0_bare came from the geometric round-trip — not where it was defined as `(1+r_c)/3` two lines earlier.

**Verification:**
After the fix, `kappa_c = kappa0_bare/(1+r_c)` should still simplify to `1/3`, but `kappa0_bare` must now be obtained from `4*L_W_required**2/(pi**2*a**2)` (with `L_W_required` the symbolic solve result), not from the literal `(1+r_c)/3`. New round-trip assertion: `simplify(4*L_W_required**2/(pi**2*a**2) - (1+r_c)/3) == 0`.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:44-49`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:49-51`

**What's wrong:**
SymPy lines 44-49:
```
D_bare = sp.expand((1 + r_c) * (1 - z**2 / 3 - sp.I * z**5 / 9))
D_final = sp.simplify(D_bare / (1 + r_c))
expect_zero("bare scaled-canonical branch renormalizes to canonical",
            D_final - (1 - z**2 / 3 - sp.I * z**5 / 9))
```
`D_bare` is defined as exactly `(1+r_c)·canonical` and then divided by `(1+r_c)`. The assertion that the quotient equals canonical is algebraically guaranteed. Same structure in Mathematica (lines 49-51).

The paper claim (notes section "Bare outgoing normalization") is that the bare outgoing branch is a *pure-scale deformation* of the canonical branch — i.e., the deformation factor is identically `(1+r_c)` with no z-dependence in the prefactor and that this is consistent with γ_0 = (1+r_c)/9 (which would otherwise require deriving γ_0 from the bare branch). The script's check confirms the trivial algebraic factoring, not the physical claim that the deformation is pure-scale.

**Why this matters:**
A "pure-scale" claim has content only when the bare branch is derived from somewhere else and then shown to factor as `(1+r_c)·canonical`. As written, the script asserts `(c·X)/c = X` for `c = 1+r_c`, `X = 1 - z²/3 - i z⁵/9`. The check has no physics in it.

**Required change:**
At minimum, derive γ_0 = (1+r_c)/9 from the D_bare polynomial (extracting the z⁵ coefficient and multiplying by appropriate power-counting), and assert that the derived γ_0 matches `(1+r_c)/9`. Concretely: `gamma0_from_D = -sp.I * (D_bare - (1+r_c)) ... ` extract z⁵ coefficient via `sp.Poly(D_bare, z).nth(5)` or `D_bare.coeff(z, 5)`. Assert `expect_zero("gamma0 from bare D matches", coeff_of_z5(D_bare)/(-sp.I) - (1+r_c)/9)`. Then the renormalization claim `gamma0_bare/(1+r_c) = 1/9` exercises that the z⁵ coefficient of D_bare divided by (1+r_c) gives 1/9 — which is the physical "pure-scale" assertion.

**Verification:**
After the fix, the script must extract a coefficient from `D_bare` (a polynomial-level extraction step that is not algebraically guaranteed to give `(1+r_c)/9`) and compare it to the gamma0 claim. The bare γ_0 then becomes a derived quantity, not a defined-then-asserted one.

### F3 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:26`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl:34`

**What's wrong:**
SymPy line 26: `kappa0_from_tube = sp.simplify(4 * L_W**2 / (sp.pi**2 * a**2))` and Mathematica line 34: `kappa0FromTube = FullSimplify[4*lW^2/(Pi^2*a^2), Assumptions -> $Assumptions]`. This expression is the *output* of the D/N eigenfunction analysis (notes paper claim 2): solving `q'' + k²q = 0` with `q(0)=0, q'(L_W)=0` gives `k_W = π/(2L_W)`, hence `Ω_W = k_W c_s = πc_s/(2L_W)`, hence `κ_0 = (ω²/Ω_W²)/z² = 4L_W²/(π²a²)` with `z = aω/c_s`.

The scripts skip every step of that derivation and start at the answer.

**Why this matters:**
The notes title the relevant subsection "Bare mixed D/N tube" and box the eigenvalue formula explicitly. If the D/N boundary conditions were different (e.g., DD instead of DN), `k_W` would change to `π/L_W` and `κ_0` to `L_W²/(π²a²)` — a factor of 4 error. Hardcoding `4L_W²/(π²a²)` cannot detect such a setup error, and the rest of the cascade (`L_W = πa√((1+r_c)/3)/2`) would silently absorb the wrong constant.

**Required change:**
Add an explicit eigenvalue derivation block before line 26 (and before line 34 in Mathematica):
1. Symbol declaration: define a generic eigenfunction `q(x) = sin(k*x)` (already satisfies `q(0)=0`), then impose the second BC `q'(L_W) = k*cos(k*L_W) = 0` → smallest positive root `k_W = π/(2*L_W)`.
2. Compute `Omega_W = k_W * c_s` and let `z = a*ω/c_s`.
3. Build `kappa0_derived = (ω/Omega_W)**2 / z**2` and `simplify`.
4. Assert `expect_zero("kappa0 matches geometric expression", kappa0_derived - 4*L_W**2/(pi**2*a**2))`.

Then keep `kappa0_from_tube = 4*L_W**2/(pi**2*a**2)` as a notational alias of the *verified* expression. The mathematica fix mirrors this: use `Solve[Cos[k*lW] == 0 && k > 0, k]` or directly write `kW = Pi/(2 lW)` and verify by substitution that the BC `D[Sin[kW*x], x] /. x -> lW` simplifies to 0.

**Verification:**
A new assertion appears in both engines that links the eigenvalue boundary condition `cos(k*L_W) = 0` and the derived `κ_0` formula. The new sympy assertion should appear before line 26 (current) and yield 0 on simplify.

### F4 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:27-33`

**What's wrong:**
The SymPy script obtains `L_W_required` via `sp.solve` (lines 27-30), then *prints* it (line 33) but never asserts it. The Mathematica script (line 39) does have a tube-length assertion: `expectZero["tube-length law", lWRequired - (Pi*a*Sqrt[(1 + rC)/3])/2]`. The SymPy script is missing the corresponding check.

The paper card's body equation (line 16, stage_116.tex) is exactly this identity: "First D/N mixed half-wave fixes \(L_W=\pi a\sqrt{(1+r_c)/3}/2\)". This is the central stage claim and the SymPy engine does not test it.

**Why this matters:**
With both A1-A3 tautological and A4 (the tube-length law) only present in Mathematica, the SymPy script has *zero* non-tautological assertions exercising any paper claim. If the Mathematica engine were retired or had a transliteration drift, the SymPy script alone would not catch a bug in the central identity of the stage.

**Required change:**
After line 30 (after the solve), add:
```
expect_zero(
    "tube-length law: L_W = pi a sqrt((1+r_c)/3) / 2",
    sp.simplify(L_W_required - sp.pi * a * sp.sqrt((1 + r_c) / 3) / 2),
)
```
The `sp.solve` returns a list; the script already takes `[0]`. If sympy returns the negative branch first (it should not given `positive=True`, but defensive), use `sp.Max(*sp.solve(...))` or filter for positive solutions before subtracting.

**Verification:**
A new line appears in the sympy output reading `tube-length law: ... = 0`. Engine cross-check now becomes symmetric: both engines test the same identity (A4 ↔ new sympy assertion).

## Independent-derivation check (Mathematica)

The `.wl` file mirrors the SymPy variable choreography one-to-one:
- SymPy `kappa0_from_tube = sp.simplify(4 * L_W**2 / (sp.pi**2 * a**2))` ↔ Mathematica `kappa0FromTube = FullSimplify[4*lW^2/(Pi^2*a^2), ...]`
- SymPy `L_W_required = sp.solve(sp.Eq(kappa0_from_tube, (1+r_c)/3), L_W)[0]` ↔ Mathematica `lWRequired = FullSimplify[lW /. First[Solve[kappa0FromTube == (1 + rC)/3, lW, Reals]], ...]`
- SymPy `D_bare = sp.expand((1+r_c)*(1 - z**2/3 - sp.I*z**5/9))` ↔ Mathematica `dBare = Expand[(1+rC)*(1 - z^2/3 - I*z^5/9)]`

This is a literal port: same intermediate names, same call shapes, same intermediate expressions. The Mathematica script does add one independent assertion (the tube-length law, line 39) that SymPy lacks, but the underlying algebraic flow is transliterated. The hardcoded `4*L_W²/(π²a²)` is the load-bearing expression in both engines, and neither engine derives it from the eigenvalue boundary-value problem the notes describe. I am flagging this as a structural concern but **not** filing it as a `mathematica_transliteration` finding because (a) Mathematica adds one independent assertion not present in SymPy, and (b) the load-bearing fix is the F3 eigenvalue derivation, which should be done independently in each engine. After F3 is applied (separately in each engine, using each engine's native ODE/Solve idiom), the transliteration concern is resolved.

## Engine cross-check

Both engines exit 0 with PASS outputs. The numerical/algebraic forms match modulo simplification: SymPy `pi*a*sqrt(3*K_q*K_s + 3*lam**2)/(6*sqrt(K_q)*sqrt(K_s))` and Mathematica `(a*Sqrt[3 + (3*lam^2)/(kQ*kS)]*Pi)/6` are equivalent (factor sqrt(3)/sqrt(3) and rewrite). No engine disagreement at the level the scripts claim. The output transcripts are fresh (mtime > script mtime for both engines).

## Verdict justification

The two engines produce internally consistent transcripts, but three of the four script-side assertion families are tautological (F1, F2 — algebraic identities of the form `(c·X)/c = X` for c = 1+r_c, X = canonical). The fourth load-bearing quantity (the D/N geometric κ_0) is hardcoded with no derivation from the eigenvalue problem the notes describe (F3). The SymPy script additionally does not assert the central paper claim (the tube-length law) — the Mathematica script does, but SymPy only prints `L_W_required` without comparison (F4). Verdict is `findings`, not `clean`. No `stop_cold` — none of the findings is incompatible with the paper, and downstream stages cite the tube-length law as a forward; once the SymPy assertion is added and the κ_0 derivation is in place, the chain holds.

## Self-test notes

- **Variable independence (trap 1):** The proposed F3 fix uses `q(x) = sin(k*x)` then `D[q, x] /. x -> L_W = k*cos(k*L_W)`. The variables `q`, `k`, `x`, `L_W` are independent symbols; the derivative wrt x is genuinely nonzero and the BC condition `cos(k*L_W) = 0` truly constrains k. The smallest positive root is k=π/(2*L_W), nonzero. Not a vacuous derivative.
- **Trivial-case pre-check (trap 3):** For F4, substituting `r_c = 0`: `L_W_required = πa·sqrt(1/3)/2 = πa/(2√3)`. RHS: `πa·sqrt(1/3)/2 = πa/(2√3)`. Residual = 0. For `r_c → ∞`: both sides scale as `πa·sqrt(r_c/3)/2`. Residual = 0. The assertion is non-vacuously zero across the parameter range.
- **Paper round-trip (trap 5):** All proposed fixes (F1-F4) require no paper-side change; the paper card already states `L_W = πa√((1+r_c)/3)/2` (line 16), the notes derive κ_0 from the D/N eigenvalue problem (notes lines 9-25), and the canonical 1/3, 1/9 final values match the notes' boxed values. No new `paper_misalignment` is introduced.
