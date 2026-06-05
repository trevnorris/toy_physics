---
unit_id: 001
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md]
  paper_appendix: present
---

# Audit unit 001 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_001.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 24: stage 001)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage001_geometry_lift_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage001_geometry_lift_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage001_geometry_lift_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage001_geometry_lift_mathematica_audit.txt`

## What the paper claims

Stage 001 is the foundational geometry-lift stage. It replaces the finite-dimensional collective variables `a(t), L(t)` by a distributed throat shape field `Σ(X,t)=r-R(Ω,w,t)` (eq:app-stage001-sigma), recovering `a,L` as lowest moments, and sets up the first linearized matter/gauge/geometry system. Its `\stagefield{Output}` (tex:223) verbatim: *"The stage outputs the shape-field closure \eqref{eq:app-stage001-sigma}, the confinement variation \eqref{eq:app-stage001-delta-v}, the harmonic split \eqref{eq:app-stage001-harmonic-expansion}, the modal wall PDE \eqref{eq:app-stage001-modal-wall-pde}, and the response-operator target \eqref{eq:app-stage001-response-operator}. These are the definitions without which later grouped response coefficients have no PDE-side interpretation."* The card's `\stagefield{Verification}` (tex:11) narrows what the audits should check to three load-bearing identities: *"the monopole/P2 harmonic bookkeeping, the confinement chain-rule sign, and the densitized-versus-weighted wall-action convention."* The `\stagefield{Checks}` checklist (tex:225-230) lists: the moment relation `q_{00}=2√π δa`, the sign of `δV_conf` from `δΣ=-η`, harmonic-orthogonality selection (no l=0↔l=2 mixing), and the ontology check (keep `A_w,J^w,F_{μw}` alive). The notes add the boxed linearized Maxwell law `∂_M(Z(w)δF^{MN})+(1/ξ)∂^N(∂·δA)=μ0 δJ^N` (§8) and the eigenvalue specialization `−Δ_{S²}Y_{lm}=l(l+1)Y_{lm}` driving the modal split (§5). It is a checkpoint stage and is `is_status_only_candidate: False`, so both engines are required.

## What the script claims to verify

The SymPy script (docstring lines 8-22) and the assertions verify: (I) normalized real-harmonic bookkeeping — `Y00=1/(2√π)`, unit norms, zero averages for the five real P2 modes, orthogonality `cross(Y00,Y20)=0`, the mouth-average extraction `q00=2√π δa`, and the spherical-Laplacian eigenvalue `−Δ_{S²}=l(l+1)` (i.e. `6` for l=2); (II) the confinement chain-rule sign `δV_conf=−(V'_wall(Σ0/ℓc)/ℓc)η` via a genuine ε-variation; (III) the modal wall Euler–Lagrange equation in both the densitized and weighted-surface conventions, the `K_l` restoring-shift specialization at l=0,2, and the sourced RHS; (IV) a representative two-component (x,w) localized-Maxwell linearization with the `Z(w)` weight, gauge-fixing term, and source current. The Mathematica `.wl` mirrors the same checks and adds an independent `SphericalHarmonicY`-based eigenvalue confirmation (I.3b) and uses `VariationalD` for the Maxwell sector. Every check is an `expect_zero`/`expectZero` against a hand-written target; the saved outputs show all residuals `= 0` and all `PASS`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Shape-field closure Σ=r−R (eq:sigma, Output item 1) | definitional; no checkable identity, but the harmonic structure on the mouth S² it enables is exercised in Section I | match (definition) |
| Confinement variation δV_conf=−(V'_wall/ℓc)η (eq:delta-v, boxed; Output item 2; Check 2) | Section II ε-variation, py:130-142 / wl:141-153 | match |
| Harmonic split + P2 set + no l=0↔l=2 mixing (eq:harmonic-expansion, Output item 3; Check 3) | Section I norms/averages/`cross(Y00,Y20)`, py:82-106 / wl:89-103 | match |
| q_{00}=2√π δa (Check 1) | Section I.2 mouth-avg, py:108-115 / wl:105-112 | match |
| Modal wall PDE (eq:modal-wall-pde, boxed; Output item 4) incl. l(l+1) eigenvalue and densitized/weighted forms (tex:149) | Section I.3 eigenvalues + Section III EL densitized & weighted & K_l shift & source, py:117-127,145-192 / wl:114-139,155-192 | match |
| Response-operator target j_A=ΣZ_eff,AB u_B (eq:response-operator, Output item 5) | no direct check; it is a *definitional target* whose substantive content (block-diagonalization on isotropic throat) rests on the Section I orthogonality, which IS checked | match (definition; substantive content covered) |
| Linearized Maxwell ∂_M(Z δF)+gauge=μ0 δJ (eq:linear-maxwell, boxed body; notes §8; ontology Check 4) | Section IV representative (x,w) EL, py:195-229 / wl:194-213 | match (representative, as paper states) |

`paper_alignment: aligned`. Every Output/Check deliverable maps to a faithful script-side check or is a definition whose substantive content is exercised. No mismatches, no orphaned (extra) script checks — the Maxwell section corresponds to a boxed body equation + notes §8 + Check 4.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 90-93 | `expect_zero(average(Y2m))` ×5 | harmonic split / zero sphere-average of P2 | yes |
| A2 | sympy | 95-98 | `norm(Y00)-1==0`, `norm(Y20)-1==0` | harmonic normalization | yes |
| A3 | sympy | 103-106 | `cross(Y00,Y20)==0` | no l=0↔l=2 mixing (Check 3) | yes |
| A4 | sympy | 114 | `mouth_avg - q00/(2√π)==0` | q00=2√π δa (Check 1) | yes |
| A5 | sympy | 125-127 | `−Δ_{S²}Y - l(l+1)Y==0` | eigenvalue → modal split | yes |
| A6 | sympy | 141 | `first_var - target==0` | δV_conf sign (Check 2, eq:delta-v) | yes |
| A7 | sympy | 167 | densitized EL `==0` | modal wall PDE (eq:modal-wall-pde) | yes |
| A8 | sympy | 177 | weighted EL `==0` | weighted-vs-densitized convention (tex:149) | yes |
| A9 | sympy | 181-182 | `K_l(0)−K_η`, `K_l(2)−(K_η+6T_Ω)` | l(l+1) restoring shift | yes |
| A10 | sympy | 192 | sourced EL `==0` | sourced RHS of modal PDE | yes |
| A11 | sympy | 228-229 | Maxwell x,w EL `==0` | eq:linear-maxwell (representative) | yes |
| B1 | mathematica | 96-103 | mirror of A1-A3 | same | yes |
| B2 | mathematica | 112 | mirror of A4 | q00=2√π δa | yes |
| B3 | mathematica | 121-126 | mirror of A5 (real-harmonic route) | eigenvalue | yes |
| B4 | mathematica | 138-139 | `lapEig[0]`, `lapEig[2]` via `SphericalHarmonicY` | eigenvalue (independent route) | yes |
| B5 | mathematica | 153 | mirror of A6 | δV_conf sign | yes (parallel; see note) |
| B6 | mathematica | 173-182 | mirror of A7-A9 | modal wall PDE | yes |
| B7 | mathematica | 192 | mirror of A10 | sourced RHS | yes |
| B8 | mathematica | 212-213 | Maxwell EL via `VariationalD` | eq:linear-maxwell | yes |

No tautological rows: every assertion compares an *independently computed* quantity (a closed-form sphere integral, an `euler_equations`/`EulerEquations`/`VariationalD` result, or a differentiation-based variation) against a hand-written target. None of the checks restate a definition against itself.

## Findings

None.

### Adversarial attacks attempted (all failed to break the scripts)

1. **Sign of δV_conf.** Verified by hand: `d/dε Vwall((Σ0−εη)/ℓc)|_{ε=0} = −(η/ℓc)V'_wall(Σ0/ℓc)`, and the script's target `−η·∂_{Σ0}Vwall(Σ0/ℓc) = −(η/ℓc)V'_wall`. Both carry the 1/ℓc and the minus sign. Matches boxed eq:app-stage001-delta-v (tex:91) and notes §3. No sign flip.
2. **Euler–Lagrange sign convention (Section III).** SymPy `euler_equations` returns `∂L/∂q − Σ_i d/dx_i(∂L/∂q_{x_i}) = 0`. With `ldens = ½μq_t²−½T_w q_w²−½k_l q²` this gives `−μq_tt + d/dw(T_w q_w) − k_l q`, exactly `target_dens` (py:166). Matches modal wall PDE. The weighted form correctly produces the `d/dw(g·T_w q_w)` term that, on division by g, yields the extra first-derivative term the paper (tex:149) and notes (§5) describe.
3. **Maxwell EL signs (Section IV).** Hand-derived both components. `el_Ax = d/dw(Z F_wx) + d/dx(divA)/ξ − μ0 J_x` and `el_Aw = −d/dx(Z F_wx) + d/dw(divA)/ξ − μ0 J_w`; both equal the script targets (py:225-226). The opposite-sign relation between SymPy's manual EL operator and Mathematica's `VariationalD` is correctly handled by the leading `−` on `VariationalD` (wl:208-209, comment wl:206-207).
4. **Eigenvalue constant "6".** `l(l+1)` at l=2 is 6; the script proves this is not a magic literal both by computing `−Δ_{S²}Y_{2m}=6Y_{2m}` directly (A5/B3) and independently via `SphericalHarmonicY` (B4), and by the `K_l(2)=K_η+6T_Ω` specialization (A9).
5. **q00=2√π δa.** `(1/4π)∫η dΩ = q00·(1/4π)·2√π = q00/(2√π)`, so δa=q00/(2√π) ⇒ q00=2√π δa. Matches tex:127 and notes §4. The P2 terms drop out because their averages vanish (A1) — so the script's inclusion of q20,q21c,q22c in `eta` (py:110) genuinely tests that they don't contaminate the monopole average; this is non-trivial, not scaffolding.
6. **Symbol assumptions.** SymPy: `nonzero=True` on Σ0/η/ℓc/eps (Section II) is harmless and does not over-strengthen any simplify; `0<θ<π` is implicit via the integration domain; Mathematica uses `0<θ<π`, `ℓc>0`, `ξ>0`, `μ0>0`. None of these assumptions are strong enough to force a false `==0` (the targets are subtracted, so an over-aggressive simplify would still need the algebra to actually cancel). No `symbol_assumption_error`.
7. **Stale output.** Script mtimes May 25 02:13; output mtimes May 25 17:24 (sympy) / 17:28 (mathematica) — outputs are NEWER than scripts. Not stale.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a pure transliteration. It shares the same physical premises and the same basis definitions (which is unavoidable — the harmonics are what they are), but it adds genuinely independent derivation routes the SymPy script does not have:
- **I.3b (wl:128-139):** confirms the angular-Laplacian eigenvalue `−l(l+1)` from Mathematica's built-in `SphericalHarmonicY[l,0,θ,φ]`, an entirely different object from the hand-typed real harmonics. This is a real second derivation of the eigenvalue.
- **Section IV (wl:208-209):** uses the built-in `VariationalD` variational derivative, whereas SymPy assembles the Euler–Lagrange operator term-by-term via explicit `sp.diff` of the Lagrangian w.r.t. each field-gradient (py:214-223). Two structurally different machineries producing the same Maxwell residual.
- **Section III:** Mathematica uses `EulerEquations` (its own CAS builtin) vs SymPy's `euler_equations`; the integrals in Section I are evaluated independently by each engine's `Integrate`/`integrate`.

The one section explicitly flagged as a non-independent parallel check is the confinement chain rule (wl:143-145 comment: *"chain rule on a single composite function admits no engine-independent route; this section is an intentional parallel check rather than an independent derivation"*). This honest disclosure is acceptable — a single-composite chain rule genuinely has only one route, and both engines arriving at the same `−(η/ℓc)V'_wall` is a legitimate cross-engine sanity check. No `mathematica_transliteration` finding.

## Engine cross-check

Both saved outputs report identical results: every residual `= 0`, Mathematica prints `PASS` on each. Section by section: harmonic averages/norms (`2√π`, `1`, `0`), `q00` extraction, eigenvalues (incl. the extra `SphericalHarmonicY` route in Mathematica), confinement sign, densitized + weighted EL, `K_l` shifts, sourced RHS, and both Maxwell components all agree. `engines_agree: true`.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix row first, built the claim model, then attacked the scripts. Every load-bearing identity the paper's `\stagefield{Output}`, `\stagefield{Verification}`, and `\stagefield{Checks}` require is exercised by a non-tautological assertion in both engines, and the values/signs all match the paper (boxed eqs delta-v and modal-wall-pde, the `q00=2√π δa` relation, the `l(l+1)` eigenvalue, the densitized-vs-weighted convention, and the representative localized-Maxwell law). I tried to break the sign of `δV_conf`, the Euler–Lagrange sign conventions in both Sections III and IV, the eigenvalue constant, the monopole-average derivation, and the symbol assumptions; all held. The Mathematica script is an independent second engine (with genuinely distinct routes via `SphericalHarmonicY` and `VariationalD`), not a transliteration. The two Output items without a direct numeric/symbolic assertion (the shape-field closure and the response-operator target) are definitions/targets whose substantive content (mouth-S² harmonic structure; isotropic block-diagonalization) rests on the orthogonality that IS checked. Outputs are fresh.

## Self-test notes

I checked the variable-independence trap (every `sp.diff`/`D` differentiates a field that genuinely depends on the variable — e.g. `D[expr,eps]` where `expr` contains `eps`, and the EL `diff` w.r.t. field gradients that the Lagrangian actually contains; no identically-zero derivatives masquerading as passes). I checked parity/symmetry on the S² integrals (the P2 modes are odd-or-mixed in cos θ / in φ against the even measure, so their averages genuinely vanish; Y00 is constant so its average is the nonzero `2√π` — consistent with output). Trivial-case substitution on each `expect_zero` confirms the residual is the difference of two algebraically-equal expressions (not a vacuously-true `0==0`). No directive is written (zero findings).

## Value Reconciliation (pass-2 augmentation)

This is a foundational/definitional stage: it emits symbolic identities and structural results, not numeric figures-of-merit. Every emitted *result* value is a symbolic relation, and each is reflected in the `.tex` card and/or the `.md` notes.

| value | source (py/wl + output line) | .tex / .md location | status |
|---|---|---|---|
| `Y00 = 1/(2√π)` (and ∫Y00 dΩ = 2√π) | py:66 / wl:68; out line 9 (`average = 2*sqrt(pi)`) | tex:116; md:163 | MATCH |
| Five real P2 modes have zero sphere-average | py:90-93 / wl:96-100; out lines 15-19 | tex:109-113 (P2 set), md:165 ("zero sphere average") | MATCH |
| Unit norms `‖Y00‖=‖Y20‖=1`, `⟨Y00,Y20⟩=0` | py:95-106 / wl:101-103; out lines 20-22 | tex:118 (`∫Y00²=1`); orthogonality implied by Check 3 (tex:228), md:147,165 | MATCH |
| `q00 = 2√π δa` | py:114 / wl:112; out line 27 | tex:127 (boxed-relation eq:q00-da), tex:226 (Check 1); md:171 | MATCH |
| `−Δ_{S²}Y_{lm} = l(l+1)Y_{lm}` (=6 at l=2) | py:125-127 / wl:121-139; out lines 33-38 | tex:153 + tex:161; md:235 + md:239 | MATCH |
| `δV_conf = −(V'_wall(Σ0/ℓc)/ℓc) η` | py:141 / wl:153; out line 50/65 | tex:89-91 (boxed eq:delta-v), tex:227; md:135, md:303, md:351 | MATCH |
| Densitized modal-wall EL `−μ q_tt + ∂_w(T_w q_w) − [K_η+l(l+1)T_Ω] q = 0` | py:166-167 / wl:172-173; out line 60 | tex:158-162 (boxed eq:modal-wall-pde); md:239 | MATCH |
| Weighted-surface EL form (extra `g'` first-derivative term) | py:170-178 / wl:175-179; out lines 65-66 | tex:149; md:224-225 | MATCH |
| `K_l(l=0)=K_η`, `K_l(l=2)=K_η+6T_Ω` | py:181-182 / wl:181-182; out lines 71-73 | tex:161 (`K_η+l(l+1)T_Ω`); md:239 | MATCH |
| Sourced RHS `S_lm^{(ψ,A)} + f_lm^{ext}` on modal PDE | py:185-192 / wl:184-192; out line 78 | tex:162 (eq RHS); md:361 | MATCH |
| Linearized Maxwell `∂_M(Z δF^{MN}) + (1/ξ)∂^N(∂·δA) = μ0 δJ^N` (representative x,w) | py:225-229 / wl:210-213; out lines 83-84 | tex:194-197 (boxed eq:linear-maxwell); md:317 | MATCH |

INTERNAL items (verification scaffolding, no prose expected, no finding): the `expect_zero`/`expectZero` residuals themselves (all `= 0`); the `PASS:` flags; the abstract symbol/function placeholders (`μ_eta, T_w, T_Ω, K_η, g(w), Z(w), A_x, A_w, J_x, J_w, Vwall, q(t,w)`); the auxiliary mode amplitudes `q20,q21c,q22c` used only to probe the monopole average; the `eps` variation parameter.

reconciliation: complete; 11 values checked, 0 misaligned.
