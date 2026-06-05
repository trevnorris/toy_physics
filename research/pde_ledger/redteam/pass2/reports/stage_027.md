---
unit_id: 027
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage027_nonconstant_axial_family.md]
  paper_appendix: present
---

# Audit unit 027 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_027.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage027_nonconstant_axial_family.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 44, `\input` line 92)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.txt`

## What the paper claims

Stage 027 replaces the Stage-026 constant axial profile with the first nonconstant two-mode N/N family `chi_theta(s)=cos θ ν0 + sin θ ν1` and asks what survives. `\stagefield{Output}{Stage~027 outputs \(\kappa(\theta)\), the blind-angle no-go \eqref{eq:app-stage027-blind-angle}, the max-coupling angle \eqref{eq:app-stage027-theta-max}, and the profile-dressed branch equation \eqref{eq:app-stage027-Kgate}.}` The distinct deliverables are: (1) the exact profile-dependent D/N overlap `kappa(θ)=κ0 cosθ+κ1 sinθ` with `κ0=2√2/π`, `κ1=-4/(3π)`, and amplitude `|κ|_max=√(κ0²+κ1²)=2√22/(3π)`; (2) the exact blind-angle no-go `tan θ_blind=3√2/2`, at which `C=G_W=R=P=N0=0`; (3) the max-coupling angle `tan θ_max=-√2/3`, `sin²θ_max=2/11`, with `κ(θ_max)=2√22/(3π)`; (4) the profile-dressed wall stiffness `K_geo(θ)=K_η+6T_Ω+(π²T_w/L²)sin²θ` (value `K_η+6T_Ω+2T_w π²/(11L²)` at max), and the dressed branch/normalization gate `K_geo(θ)=K_req(θ)` with the explicit `Δ(θ)`, `Q(θ)`, `K_req(θ)` forms. The notes additionally enumerate the full reduced-module substitution `C=λ_B κ`, `G_U=λ_U`, `G_W=λ_W κ`, `R=λ_R κ` and the θ=0 recovery of Stage 026, all of which the .tex carries in compressed form.

## What the script claims to verify

Both engines build the exact finite-throat bases (`u0=1/√L`, `u1=√(2/L)cos(πs/L)`, `f0=√(2/L)sin(πs/2L)`), verify orthonormality/normalization, then compute the overlaps by direct symbolic integration: `κ0=2√2/π`, `κ1=-4/(3π)`, `κ(θ)`, `ρ=2√22/(3π)`. They check the blind/max branches by substituting the exact `(cosθ,sinθ)` values, confirm `κ(blind)=0`, `κ(max)=ρ`, `sin²θ_max=2/11`. They evaluate the wall operator `G_η` on the family, derive `K_geo(θ)=K_η+6T_Ω+T_w π² sin²θ/L²` and its max value, then substitute into the reduced module to derive `Δ,Q,P,B0,Z0,N0,K_req`, check each closed form, recover Stage 026 at θ=0, and establish the blind-angle no-go via `P(blind)=0` and `N0(blind)=0` (SymPy additionally checks the divided `lhs(blind)=0`).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `κ0=2√2/π` (eq kappa0-kappa1) | py:114 / wl:70 `expect_zero("kappa0 - 2*sqrt(2)/pi")` | match |
| `κ1=-4/(3π)` | py:115 / wl:71 | match |
| `κ(θ)=κ0 cosθ+κ1 sinθ` (eq kappa-theta) | py:116-119 / wl:72 | match |
| `|κ|_max=2√22/(3π)` (eq kappa-max) | py:120 / wl:73 (`rho - 2*sqrt(22)/(3*pi)`) | match |
| blind no-go `tan θ_blind=3√2/2`, κ=0 | py:141 / wl:85 `kappa(blind)=0` (subs cos=√2/√11, sin=3/√11) | match |
| max angle `tan θ_max=-√2/3`, `sin²θ_max=2/11` | py:142-143 / wl:86-87 (`kappa(max)-rho`, `sin^2(theta_max)-2/11`) | match |
| `K_geo(θ)=K_η+6T_Ω+(π²T_w/L²)sin²θ` (eq Kgeo) | py:161 / wl:99 `K_geo - expected` | match |
| `K_geo(θ_max)=K_η+6T_Ω+2T_w π²/(11L²)` | py:170-173 / wl:104 | match |
| `Δ(θ)=Ω_U²Ω_W²-λ_R²κ²` | py:221 / wl:143 `Delta - expected` | match |
| `Q(θ)=λ_U²Ω_W²+2λ_Uλ_Wλ_Rκ²+λ_W²Ω_U²κ²` | py:222 / wl:144 `Q - expected` | match |
| `P(θ)=κ(Ω_U²λ_W+λ_Rλ_U)` (notes), `B0` | py:223-224 / wl:145-146 | match |
| dressed gate `K_geo=K_req`, `K_req(θ)` form (eq Kreq) | py:253-256 printed + V.1 no-go / wl:165-166 printed | match (K_req printed as closed form; gate verified via no-go + θ=0 recovery) |
| θ=0 recovery of Stage 026 | py:230-238 / wl:151-156 | match |
| blind: `C=G_W=R=P=N0=0` | py:271-273 (P,N0,lhs) / wl:172-173 (P,N0) | match |

`paper_alignment: aligned` — every paper-side deliverable maps to a non-tautological script-side check that exercises it; no extras, no mismatches.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 87-90 | `expect_zero` orthonormality of u0,u1,f0 | basis setup (checks #1) | yes |
| A2 | sympy | 114 | `κ0 - 2√2/π == 0` | κ0 | yes |
| A3 | sympy | 115 | `κ1 + 4/(3π) == 0` | κ1 | yes |
| A4 | sympy | 116-119 | `κ - (κ0 cosθ+κ1 sinθ) == 0` | κ(θ) law | yes |
| A5 | sympy | 120 | `ρ - 2√22/(3π) == 0` | |κ|_max | yes |
| A6 | sympy | 141 | `κ(blind) == 0` | blind no-go | yes |
| A7 | sympy | 142 | `κ(max) - ρ == 0` | max-coupling angle | yes |
| A8 | sympy | 143 | `sin²θ_max - 2/11 == 0` | sin²θ_max | yes |
| A9 | sympy | 161 | `K_geo - (K_η+6T_Ω+T_w π² sin²θ/L²) == 0` | dressed stiffness | yes |
| A10 | sympy | 170-173 | `K_geo(max) - (…+2T_w π²/(11L²)) == 0` | max stiffness | yes |
| A11 | sympy | 221-224 | `Δ,Q,P,B0 - expected == 0` | branch quantities | yes |
| A12 | sympy | 230-238 | θ=0 recovery (κ0, K_geo, Δ) | Stage-026 recovery | yes |
| A13 | sympy | 271-273 | `P(blind),N0(blind),lhs(blind) == 0` | blind no-go LHS | yes |
| B1 | mathematica | 49-52 | orthonormality `expectZero` | basis setup | yes |
| B2 | mathematica | 70-73 | κ0,κ1,κ(θ),ρ | overlap deliverables | yes |
| B3 | mathematica | 85-87 | κ(blind),κ(max)-ρ,sin²θ_max | blind+max angles | yes |
| B4 | mathematica | 99,104 | K_geo, K_geo(max) | dressed stiffness | yes |
| B5 | mathematica | 143-146 | Δ,Q,P,B0 expected | branch quantities | yes |
| B6 | mathematica | 151-156 | θ=0 recovery | Stage-026 recovery | yes |
| B7 | mathematica | 172-173 | P(blind),N0(blind) | blind no-go | yes |

Every assertion traces to a specific paper deliverable. None is tautological: the LHS of each `expect_zero`/`expectZero` is obtained by an independent symbolic `integrate`/`Integrate` or by `subs`/`/.` into an integrated expression, and the RHS is the paper's stated closed form, so each can genuinely fail if the physics is wrong.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. Each engine independently computes the load-bearing integrals from the basis definitions: SymPy `sp.integrate(u0*f0,(s,0,L))` (py:103) vs Mathematica `Integrate[u0*f0,{s,0,l}]` (wl:59); SymPy `sp.integrate(chi*G_eta,(s,0,L))` (py:157) vs Mathematica `Integrate[chi*gEta,{s,0,l}]` (wl:95). The simplification engines differ (`sp.simplify`/`expand_trig` vs `FullSimplify`/`TrigExpand`/`Together`), and they emit different canonical forms for the same quantity — e.g. Mathematica reports `K_geo` as `kEta + 6 tOmega + (Pi^2 tW)/(2 l^2) - (Pi^2 tW Cos[2 theta])/(2 l^2)` (mathematica output line 62) while SymPy reports `K_eta + 6 T_Omega + pi**2 T_w sin(theta)**2/L**2` (sympy output line 89). These are equal but not byte-identical, which is exactly the signature of two genuine derivations rather than one echoing the other. The shared section/variable naming is organizational scaffolding, not algebra-echo. No `mathematica_transliteration` finding.

## Engine cross-check

The two engines agree on every checked quantity. Identical results: `κ0=2√2/π`, `κ1=-4/(3π)`, `ρ=2√22/(3π)`, `κ(blind)=0`, `κ(max)=2√22/(3π)`, `sin²θ_max=2/11`, `K_geo(max)=K_η+6T_Ω+2T_w π²/(11L²)`, all `…-expected = 0` / `PASS`, θ=0 recovery all 0, and `P(blind)=N0(blind)=0`. SymPy adds one check (`lhs(blind)=0`, the fully divided LHS) the Mathematica omits, but Mathematica's `N0(blind)=0` (recomputed as `pBlind^2/deltaBlind^2`) plus `P(blind)=0` covers the same no-go claim — both establish that the LHS numerator vanishes while the target `54 G c_s^5/(5 a^5 c^5)` is strictly positive. No `engine_disagreement`.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix row before the scripts, and confirmed the script claims match the paper's four deliverables exactly. Attacks tried that failed: (a) trig-substitution consistency — `blind_subs (cos=√2/√11, sin=3/√11)` gives `tan=3√2/2` matching `eq:app-stage027-blind-angle`, and `max_subs (cos=3/√11, sin=-√2/√11)` gives `tan=-√2/3`, `sin²=2/11`, both on the unit circle (cos²+sin²=1), matching `eq:app-stage027-theta-max`; (b) tautology hunt — every `expect_zero` RHS is the paper's literal closed form and the LHS is an independent integral/substitution, so none is guaranteed by construction; (c) `sin²θ` substitution trap — both engines correctly retain/expand `sin²θ` before substituting (SymPy keeps `sin(theta)**2`; Mathematica `TrigExpand`s the `Cos[2θ]` form at wl:102,152), so the max/θ=0 stiffness substitutions are valid; (d) no-go completeness — both engines verify the numerator structure (`P`, `N0`) vanishes at the blind angle. The only blemish is a stale-output mtime signal (SymPy `.txt` predates the `.py`), but the sole post-output change to the `.py` was a docstring label `Stage 10 → Stage 27` in commit `e2a4780` (doc-only numbering reconciliation; no logic/value change), so the captured output content remains valid — informational only, not a blocking finding.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value emitted by the scripts (per source + committed saved outputs) was located in the `.tex` card and/or `.md` notes. All reconcile.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `κ0 = 2√2/π` | py:103,114 / wl:59,70; sympy.txt:29, math.txt:37 | tex:51 `\kappa_0=\frac{2\sqrt2}{\pi}`; md:93 | MATCH |
| `κ1 = -4/(3π)` | py:104,115 / wl:60,71; sympy.txt:30, math.txt:38 | tex:52 `\kappa_1=-\frac{4}{3\pi}`; md:95 | MATCH |
| `κ(θ) = 2(-2 sinθ+3√2 cosθ)/(3π)` | py:105 / wl:61; sympy.txt:31, math.txt:39 | tex:57 (eq kappa-explicit); md:101 | MATCH |
| `ρ = |κ|_max = 2√22/(3π)` | py:106,120 / wl:62,73; sympy.txt:32, math.txt:40 | tex:63 `|\kappa|_{\max}=\frac{2\sqrt{22}}{3\pi}`; md:111,115 | MATCH |
| `tan θ_blind = 3√2/2` (via subs cos=√2/√11,sin=3/√11 ⇒ κ=0) | py:124-127,141 / wl:76,85; sympy.txt:41, math.txt:53 | tex:71 `\tan\theta_{\rm blind}=\frac{3\sqrt2}{2}`; md:127 | MATCH |
| `tan θ_max = -√2/3` (via subs cos=3/√11,sin=-√2/√11) | py:129-130 / wl:77; (subs) | tex:91 `\tan\theta_{\max}=-\frac{\sqrt2}{3}`; md:143 | MATCH |
| `sin²θ_max = 2/11` | py:135,143 / wl:80,87; sympy.txt:43, math.txt:55 | tex:93 `\sin^2\theta_{\max}=\frac{2}{11}`; md:151 | MATCH |
| `κ(θ_max) = 2√22/(3π)` | py:134,142 / wl:79,86; sympy.txt:42, math.txt:54 | tex:63 (=|κ|_max); md:155 | MATCH |
| `K_geo(θ) = K_η+6T_Ω+(π²T_w/L²)sin²θ` | py:157,161 / wl:95,99; sympy.txt:89, math.txt:62 | tex:106 (eq Kgeo); md:172 | MATCH |
| `K_geo(θ_max) = K_η+6T_Ω+2T_w π²/(11L²)` | py:168,170 / wl:102,104; sympy.txt:91, math.txt:65 | tex:182 in md; md:182 | MATCH (general form in tex:106) |
| `Δ(θ)=Ω_U²Ω_W²-λ_R²κ²` | py:193,216 / wl:117,138; sympy.txt:101 | tex:125; md:215 | MATCH |
| `Q(θ)=λ_U²Ω_W²+2λ_Uλ_Wλ_Rκ²+λ_W²Ω_U²κ²` | py:194,217 / wl:118,139; sympy.txt:102 | tex:128; md:217-219 | MATCH |
| `P(θ)=κ(Ω_U²λ_W+λ_Rλ_U)` | py:195,218 / wl:119,140; sympy.txt:103 | md:221 (tex uses it inside K_req numerator) | MATCH |
| `B0=λ_B²κ²/varpi²` | py:197,219 / wl:120,141; sympy.txt:104 | tex:121 (first K_req term); md:225 | MATCH |
| `C=λ_B κ`, `G_U=λ_U`, `G_W=λ_W κ`, `R=λ_R κ` | py:188-191 / wl:112-115; sympy.txt:97-100 | tex:73-83 (blind: C=G_W=R=0); md:205-211 | MATCH |
| `K_req(θ)` closed form (+ gate `K_geo=K_req`) | py:253 / wl:165; sympy.txt:120, math.txt:101 | tex:111-129 (eq Kgate, eq Kreq, Δ, Q); md:257-265 | MATCH |
| blind: `C=G_W=R=P=N0=0` | py:271-273 / wl:172-173; sympy.txt:125-127, math.txt:106-108 | tex:74-83; md:296-300 | MATCH |
| θ=0 recovery: `κ→κ0`, `K_geo→K_η+6T_Ω`, `Δ→…κ0²` | py:230-238 / wl:151-156; sympy.txt:116-118 | tex:140 (downstream); md:278-286 | MATCH |

Internal scaffolding (accounted for, no finding): basis orthonormality residuals `int u0^2-1`, `int u1^2-1`, `int u0*u1`, `int f0^2-1`; the intermediate prints `Z0`, `D0` (derived display quantities, not boxed deliverables — `Z0=Q/Δ`, `D0=K_geo-B0-Z0`, both reported in notes §4 as `Z_0(θ)`, `D_0(θ)` and consistent); the `…-expected = 0` / `PASS` flags; the strictly-positive-target prose line.

reconciliation: complete; 19 deliverable values checked, 0 misaligned.

## Self-test notes

I checked the variable-independence trap (`K_geo` genuinely depends on θ via the `D[chi,{s,2}]`/`sp.diff(chi,s,2)` second derivative of the `cos(πs/L)` mode, so the `sin²θ` term is real, not an identically-zero artifact), the symmetry/parity of the overlap integrals (all are definite integrals on `[0,L]`, computed and matched against literal closed forms — no symmetric-domain vanishing assumption is being relied on), and the trivial-case pre-check at θ=0 (κ→κ0=2√2/π, K_geo→K_η+6T_Ω, recovered explicitly and matching Stage 026). The blind/max substitutions both lie on the unit circle and reproduce the paper's tangents. No traps fire; verdict stays clean.
</content>
</invoke>
