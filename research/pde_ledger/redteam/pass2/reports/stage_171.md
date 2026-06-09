---
unit_id: 171
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage171_microscopic_grouped_obstructions.md]
  paper_appendix: present
---

# Audit unit 171 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_171.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage171_microscopic_grouped_obstructions.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 73 status table; lines 416–459 §"Direct outlet map" with eqs `eq:app-part05-K-decomp`, `eq:app-part05-WBZ-defs`, `eq:app-part05-G-decomp`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage171_microscopic_grouped_obstructions_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}{Decomposes \(\mathcal K_A\) and \(\mathcal G_A\) into wall, BdG, Maxwell/mixed, and outgoing-transfer defect bundles.}` The card is terse; the notes are authoritative. Stage 171 takes the two Stage-170 linear grouped outlet obstructions `K_A = δD_{A,2} + (1/9)δD_{A,0}` and `G_A = δN_{A,0} − P0·δD_{A,0}` and makes their microscopic content explicit. The distinct deliverables are: (1) the even split `K_A = W_A − B_A − Z_A` with `W_A=(1/9)δK_A−δM_A`, `B_A=δB_{A,2}+(1/9)δB_{A,0}`, `Z_A=δZ_{A,2}+(1/9)δZ_{A,0}`; (2) the odd split `G_A = −P0·δK_A + P0·δB_{A,0} + N_A` with `N_A=δN_{A,0}+P0·δZ_{A,0}`; (3) the exact BdG linearizations of `B_{A,0}=c²/ϖ²`, `B_{A,2}=c²/ϖ⁴` and the collected `B_A` bundle; (4) the portwise Maxwell/mixed variations `δΔ, δS, δG, δP, δQ` of the invariants `Δ=Ω_U²Ω_W²−R²`, `S=Ω_U²+Ω_W²`, `Q=g_U²Ω_W²+2g_Ug_W R+g_W²Ω_U²`, `G=g_U²+g_W²`, `P=Ω_U²g_W+Rg_U`; (5) the conservative moment variations `δZ_0, δZ_2` of `Z_0=Q/Δ`, `Z_2=(QS−GΔ)/Δ²` and the collected `Z_A` bundle; (6) the outgoing-transfer variation `δN_0` of `N_0=P²/Δ²` and the collected `N_A` bundle; (7) the weak-axisymmetric `(λ_{20},λ_{21},λ_{22})=(1,1/2,−1)` collapse to the scalar pair `(𝔎₁,𝔊₁)`.

## What the script claims to verify

Both scripts assert that the paper's collected closed forms equal the term-by-term first variation. Section 1 checks the `K_A`/`G_A` split identities. Section 2 builds `δB_0, δB_2` via partial derivatives of `c²/w²`, `c²/w⁴` and checks them against the boxed formulas, then checks the collected `B_A` bundle. Section 3 builds primitive variations `δΔ, δS, δG, δP, δQ` from `D[]`/`sp.diff` of the invariant definitions and checks each against the notes boxes; then builds `δZ_0, δZ_2, δN_0` from generic symbols `Δ,S,Q,G,P` and checks them, then checks the collected `Z_A` and `N_A` bundles. Section 4 loops over the three `λ` values and checks that the microscopic-expanded `K`/`G` amplitudes equal the rebuilt Stage-170 forms. The Mathematica script additionally verifies the `Z_A` and `N_A` bundles via an independent `Series`-in-`eps2` route. Every check is an `expect_zero`/`expectZero` of a residual; all pass.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| (1) `K_A = W_A − B_A − Z_A` even split | py:39–41 / wl:31–33 `K_A split` | match |
| (2) `G_A = −P0 δK_A + P0 δB_{A,0} + N_A` odd split | py:43–45 / wl:35–37 `G_A split` | match |
| (3) `δB_0, δB_2`, collected `B_A` | py:56–69 / wl:45–58 | match |
| (4) `δΔ, δS, δG, δP, δQ` Maxwell primitives | py:87–107 / wl:70–89 | match |
| (5) `δZ_0, δZ_2`, collected `Z_A` | py:115–141 / wl:99–130 | match |
| (6) `δN_0`, collected `N_A` | py:120,132,143–148 / wl:101,111,132–143 | match |
| (7) weak-axisym `(1,1/2,−1)` collapse to `(𝔎₁,𝔊₁)` | py:155–164 / wl:148–163 | match |

`paper_alignment: aligned`. Every boxed result in the notes and every appendix decomposition equation has a corresponding non-tautological script check, in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 41 | `expect_zero(K_exact − K_split)` | (1) | yes |
| A2 | sympy | 45 | `expect_zero(G_exact − G_split)` | (2) | yes |
| A3 | sympy | 61–62 | `δB_0/δB_2` D-built vs boxed | (3) | yes |
| A4 | sympy | 69 | `B_A` bundle | (3) | yes |
| A5 | sympy | 96–107 | `δΔ,δS,δG,δP,δQ` D-built vs boxed | (4) | yes |
| A6 | sympy | 122–132 | `δZ_0,δZ_2,δN_0` vs boxed | (5),(6) | yes |
| A7 | sympy | 141 | `Z_A` bundle | (5) | yes |
| A8 | sympy | 148 | `N_A` bundle | (6) | yes |
| A9 | sympy | 163–164 | weak-axisym K/G (×3 λ) | (7) | yes |
| M1 | math | 33,37 | `K_A`/`G_A` split | (1),(2) | yes |
| M2 | math | 50–58 | BdG | (3) | yes |
| M3 | math | 78–89 | Maxwell primitives | (4) | yes |
| M4 | math | 103–111 | `δZ_0,δZ_2,δN_0` | (5),(6) | yes |
| M5 | math | 123 | `Z_A` bundle (D-route) | (5) | yes |
| M6 | math | 130 | `Z_A` bundle (Series route, independent) | (5) | yes |
| M7 | math | 138 | `N_A` bundle (D-route) | (6) | yes |
| M8 | math | 143 | `N_A` bundle (Series route, independent) | (6) | yes |
| M9 | math | 154–161 | weak-axisym K/G (×3 λ) | (7) | yes |

All assertions are anchored. None is tautological in the `x=expr; assert x==expr` sense: in every case the LHS is constructed by an independent mechanism (term-by-term `sp.diff`/`D[]` total differential, or a `Series` expansion) and the RHS is the paper's separately-typed collected closed form, so a wrong coefficient in either the engine derivation or the transcribed paper formula leaves a nonzero residual. The `δB_0` etc. checks compare a D-summed total differential against a hand-typed boxed formula — genuinely falsifiable.

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage171_microscopic_grouped_obstructions_mathematica_audit.wl:31-89,99-111,148-163`

**What's wrong:**
The bulk of the `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The variable choreography, intermediate expressions, formula targets, summation order, and even the literal print strings correspond one-to-one. Three representative pairs:

BdG (py:56–60 vs wl:45–49):
```
# py
dB0_exact = sp.diff(B0, c) * dc + sp.diff(B0, w) * dw
dB2_exact = sp.diff(B2, c) * dc + sp.diff(B2, w) * dw
dB0_formula = 2 * c * dc / w**2 - 2 * c**2 * dw / w**3
dB2_formula = 2 * c * dc / w**4 - 4 * c**2 * dw / w**5
```
```
(* wl *)
dB0Exact = D[b0, c]*dc + D[b0, w]*dw;
dB2Exact = D[b2, c]*dc + D[b2, w]*dw;
dB0Formula = 2*c*dc/w^2 - 2*c^2*dw/w^3;
dB2Formula = 2*c*dc/w^4 - 4*c^2*dw/w^5;
```

Maxwell primitives (py:87–94 vs wl:70–76): identical `diff(expr,var)*dvar` summation in identical order.

Weak-axisymmetric loop (py:158–164 vs wl:148–163): same `eps*lam*(...)` micro/rebuilt construction and the same residual checks.

The `.wl` does add two genuinely independent mechanisms — the `Series`-in-`eps2` routes for the `Z` and `N` bundles (wl:127–130, 140–143) — which compute the first-order variation by Taylor expansion rather than by summing partial-derivative terms, and compare against the same paper-typed collected target. Those two checks are real second-engine value, so this is a *partial* transliteration rather than a pure echo. But ~17 of the 19 shared checks are structurally identical ports.

**Why this matters:**
The second-engine policy is that both engines derive the result independently from the physical premises. Where the `.wl` merely re-runs the `.py`'s algebra with Mathematica syntax, a transcription error common to both engines (e.g., a mistyped boxed formula copied identically into both `*_formula`/`*Formula` lines) would not be caught — both engines would agree on the same wrong target. The two `Series` routes mitigate this only for the `Z`/`N` bundles, not for the BdG, primitive-variation, split, or weak-axisymmetric checks.

**Required change:**
Strengthen the `.wl` so the independent-mechanism coverage extends to the remaining shared checks, mirroring the `Series`-route pattern already present at wl:127–130 and wl:140–143. Concretely, for the BdG bundle (wl:53–58) and the weak-axisymmetric K/G checks (wl:148–163), add a parallel route that builds the first variation by a distinct mechanism (e.g. `Series[... /. {c -> c + eps2 dc, w -> w + eps2 dw}, {eps2,0,1}]` for BdG; and for the weak-axisymmetric collapse, derive the `(𝔎₁,𝔊₁)` amplitudes by substituting the section-1 split into the `eps*lam` ansatz rather than re-typing the collected forms) and compare against the paper's collected target. Keep the existing D[]-route checks; add the second routes alongside.

**Verification:**
After Codex applies, the Mathematica transcript should show new PASS lines such as `BdG obstruction bundle (series route)` and a weak-axisymmetric cross-route, and the script must still `Exit[0]`. The new checks must use a mechanism distinct from the `.py`'s `diff`-summation so the `.wl` is no longer a pure port of the shared sections.

## Independent-derivation check (Mathematica)

Partial. Sections 1 (split), 2 (BdG), 3 (Maxwell primitives + `δZ`/`δN`), and 4 (weak-axisymmetric) are line-by-line ports of the `.py` (see F1 quotes). The `.wl` does include two genuinely independent re-derivations the `.py` lacks: the `Series`-based first-variation of `(z2 + z0/9)` (wl:127–130) and of `(n0Expr + p0·z0)` (wl:140–143), each compared against the paper's collected closed form. Those two are not transliterations — they exercise a distinct expansion mechanism on the two most substantive bundles. The overall call is a low-severity `mathematica_transliteration` (F1): the engine adds independent value on `Z`/`N` but echoes the `.py` algebra everywhere else. Confidence: high that the non-Series sections are ports; medium that this warrants a finding rather than acceptance (the underlying math is fully correct and the two Series routes are a real independent cross-check — the finding is recorded primarily because this is the close of the 105–175 transliteration-watch band and the port is structurally undeniable).

## Engine cross-check

Both engines compute identical residual sets, all `= 0`. The SymPy transcript shows the 21 base checks at 0; the Mathematica transcript shows those same checks plus the two extra `(series route)` checks, all `= 0` / `PASS`, ending `Stage 171 Mathematica audit passed.` No disagreement. The Mathematica `Series`-route results (`Z obstruction bundle (series route) = 0`, `N obstruction bundle (series route) = 0`) agree with the SymPy `D`-route results for the same bundles, which is a meaningful agreement across two distinct mechanisms.

## Verdict justification

The scripts faithfully and non-tautologically verify every deliverable the notes and appendix state for Stage 171; paper alignment is exact (all signs, the `1/9` weights, the `(1,1/2,−1)` axisymmetric signature, and the four-bundle decomposition match). Attacks tried and failed: (a) checked the `δB_0`/`δB_2` and Maxwell-primitive checks for the `x=expr; assert x==expr` tautology pattern — they instead compare an independently summed total differential against a hand-typed boxed formula, so they are falsifiable; (b) checked the `Δ=Ω_U²Ω_W²−R²` definition against the notes given the `U≡Ω_U²` convention — consistent; (c) checked the weak-axisymmetric sign/weight collapse — matches `𝔎₁`/`𝔊₁` exactly; (d) searched for the stale `168π²`/`100π²`/`4107` constant class — none present (this stage is purely symbolic). Verdict is `findings` solely on the low-severity `mathematica_transliteration` (F1): correct math, but the `.wl` ports most of the `.py`'s algebra rather than re-deriving it, with independent value present only on the two `Series` routes. Outputs are fresh (sympy .txt 16:10 > .py 15:55; mathematica .txt 16:26 = .wl 16:26).

## Value Reconciliation (pass-2 augmentation)

This stage emits no numeric constants — it is a purely symbolic identity-verification stage. The "values" are the boxed symbolic closed forms the scripts assert and that the notes/appendix carry. Each is reconciled below.

| value (symbolic deliverable) | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `K_A = (1/9 δK_A − δM_A) − (δB_{A,2}+1/9 δB_{A,0}) − (δZ_{A,2}+1/9 δZ_{A,0})` | py:40 / wl:32; sympy.txt:5, math.txt:5 | notes:116–126 (boxed); appendix:435–449 (`eq:app-part05-K-decomp`,`-WBZ-defs`) | MATCH |
| `G_A = −P0 δK_A + P0 δB_{A,0} + (δN_{A,0}+P0 δZ_{A,0})` | py:44 / wl:36; sympy.txt:7, math.txt:7 | notes:130–137 (boxed); appendix:451–458 (`eq:app-part05-G-decomp`) | MATCH |
| `δB_0 = 2c/ϖ²·δc − 2c²/ϖ³·δϖ` | py:59 / wl:48; sympy.txt:9 | notes:166–174 (boxed) | MATCH |
| `δB_2 = 2c/ϖ⁴·δc − 4c²/ϖ⁵·δϖ` | py:60 / wl:49; sympy.txt:11 | notes:177–186 (boxed) | MATCH |
| `B_A = Σ[2c(1/ϖ⁴+1/9ϖ²)δc − 2c²(2/ϖ⁵+1/9ϖ³)δϖ]` | py:65–68 / wl:54–57; sympy.txt:13 | notes:190–200 (boxed) | MATCH |
| `δΔ = Ω_W²δΩ_U² + Ω_U²δΩ_W² − 2R δR` | py:96 / wl:78; sympy.txt:15 | notes:354–362 (boxed) | MATCH |
| `δS = δΩ_U² + δΩ_W²` | py:97 / wl:79; sympy.txt:17 | notes:365–369 (boxed) | MATCH |
| `δG = 2g_U δg_U + 2g_W δg_W` | py:98 / wl:80; sympy.txt:19 | notes:372–376 (boxed) | MATCH |
| `δP = g_W δΩ_U² + g_U δR + R δg_U + Ω_U² δg_W` | py:99 / wl:81; sympy.txt:21 | notes:379–389 (boxed) | MATCH |
| `δQ = g_W²δΩ_U² + g_U²δΩ_W² + 2g_Ug_W δR + 2(g_UΩ_W²+g_WR)δg_U + 2(g_WΩ_U²+g_UR)δg_W` | py:102–106 / wl:84–88; sympy.txt:23 | notes:392–409 (boxed) | MATCH |
| `δZ_0 = (Δ δQ − Q δΔ)/Δ²` | py:122 / wl:103; sympy.txt:25 | notes:237–241 (boxed) | MATCH |
| `δZ_2 = S/Δ²·δQ + Q/Δ²·δS − δG/Δ + (G/Δ² − 2QS/Δ³)δΔ` | py:125–130 / wl:106–109; sympy.txt:27 | notes:244–254 (boxed) | MATCH |
| `δN_0 = 2P/Δ²·δP − 2P²/Δ³·δΔ` | py:132 / wl:111; sympy.txt:29 | notes:295–301 (boxed) | MATCH |
| `Z_A = (S/Δ²+1/9Δ)δQ + Q/Δ²·δS − δG/Δ + (G/Δ² − Q/9Δ² − 2QS/Δ³)δΔ` | py:135–140 / wl:117–122; sympy.txt:31 | notes:263–273 (boxed) | MATCH |
| `N_A = P0/Δ·δQ + 2P/Δ²·δP − (P0Q/Δ² + 2P²/Δ³)δΔ` | py:144–147 / wl:134–137; sympy.txt:33 | notes:310–322 (boxed) | MATCH |
| weak-axisym `𝔎₁ = 1/9 K⁽¹⁾ − M⁽¹⁾ − (B₂⁽¹⁾+1/9 B₀⁽¹⁾) − (Z₂⁽¹⁾+1/9 Z₀⁽¹⁾)` | py:159,161 / wl:150,152; sympy.txt:39–47 | notes:426–434 (boxed) | MATCH |
| weak-axisym `𝔊₁ = N₀⁽¹⁾ − P0 K⁽¹⁾ + P0 B₀⁽¹⁾ + P0 Z₀⁽¹⁾` | py:160,162 / wl:151,153; sympy.txt:41–47 | notes:437–442 (boxed) | MATCH |
| `(λ_{20},λ_{21},λ_{22})=(1,1/2,−1)` | py:155 / wl:162; sympy.txt:39–47 | notes:419–423; card line 21 `(1,1/2,-1)` | MATCH |

INTERNAL (scaffolding, no finding): `eps`, `eps2` perturbation parameters; generic symbols `Δ,S,Q,G,P` and their variations used as standalone placeholders; banner/PASS/FAIL strings; the three printed carry-forward summary lines.

reconciliation: complete; 18 values checked, 0 misaligned

## Self-test notes

Variable-independence: every `sp.diff(EXPR,VAR)`/`D[expr,var]` is taken w.r.t. a symbol that genuinely appears in EXPR (e.g. `Q=g_U²W+2g_Ug_W R+g_W²U` depends on all of `U,W,R,g_U,g_W` — no identically-zero derivative; `S=U+W` correctly has no `R`/`g` derivative terms). Symmetry/parity: no integrals in this stage. Trivial-case: substituting any single primitive variation (e.g. only `δc≠0`) into the BdG bundle reproduces the boxed `δc` coefficient, confirming the residual genuinely vanishes only for the correct formula. Paper round-trip: the F1 fix adds independent second-route checks against the *same* paper-typed targets, introducing no new constant and no new `paper_misalignment`.
