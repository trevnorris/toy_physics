---
unit_id: 245
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 245 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_245.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage245_nonrigid_mouth_dressing_packet_and_uv_drain_compiler_mathematica_audit.txt`

## What the paper claims

Stage 245 compiles the abstract Stage-243 non-rigid `(U,V)` free-energy lane into the actual finite physical observables. From the reduced free energy
`F_nr = (1/2)a_U U^2 + (1/2)a_V V^2 - chi_UV U V - f_U U`,
the stationary packet is `U = a_V f_U/Delta_UV`, `V = chi_UV f_U/Delta_UV`, with `Delta_UV = a_U a_V - chi_UV^2` (`>0` for admissibility), and `V/U = chi_UV/a_V`. The deliverables (notes §1-8, card eqs. stage245-direct-observables … support-blind) are: (1) the exact finite compiler `T^2 = T_ref^2 exp(a_V f_U/Delta_UV)`, `eps_eta = eps_eta_ref exp(chi_UV f_U/Delta_UV)`, and the target-ratio `R_target/R_ref = [(1-eps_ref exp(V))/(1-eps_ref)] exp(-U)` derived from the selected-branch identity `R_target T^2 = Lambda_0(1-eps_eta)`; (2) the dependent microscopic correction `y_nr^dep = (0, -chi_UV f_U/Delta_UV, (a_V-chi_UV)f_U/Delta_UV)`; (3) the nonnegative drain `D_UV = chi_UV^2 a_V f_U^2/Delta_UV^2 >= 0`; (4) the weak-axisymmetric first-order packet `Xi_1^nr = u_1`, `R_1^nr = -u_1 - [eps_ref/(1-eps_ref)] v_1`, `R_1+Xi_1 = -[eps_ref/(1-eps_ref)] v_1`; (5) orbit-side / support-blind split `∂_Λ U = ∂_ρ U = ∂_Λ V = ∂_ρ V = 0`; (6) a Session-I numeric readback (`U≈0.14313458, V≈-0.03619791, eps_ref=0.3` → `eps_eta≈0.28933482, R_ratio≈0.87984149, R_1≈-0.12762119`). The appendix row 88 quotes the same: "`\StatusExactClosure{}` … Exact finite compiler for `T^2`, `eps_eta`, `R_target`; non-rigid determinant `Delta_UV`; dependent microscopic correction; positive `U/V` drain; first-order packet `(Xi_1, R_1)`."

## What the script claims to verify

The docstring enumerates eight checks matching the eight deliverables. The sympy assertions: derive the stationary `(U,V)` via `sp.solve` and confirm the closed forms, the `V/U` ratio, the Hessian determinant `= Delta_UV`, and the two limiting slices (`f_U=0 → 0,0`; `chi_UV=0 → f_U/a_U, 0`); confirm `D_UV` closed form and its positivity on the opposite-sign branch; confirm the three finite compilers (with `R_ratio` independently re-derived from the selected-branch identity and equated to the paper form); confirm `y_dep`; confirm the first-order `Xi_1/R_1` coefficients via series expansion; confirm support-blindness (zero derivatives) plus a positive control; confirm the Session-I readback against independently-recorded numbers. The Mathematica script (M1-M8) verifies the same set with `Solve`/`Series`/`D` and a hand-written Hessian.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Stationary `U,V` + `Delta_UV` + `V/U` ratio | sympy L69-72 / M1 L109-112 | match |
| Limiting slices (`f_U=0`, `chi_UV=0`) | sympy L73-76 / M1 L113-116 | match |
| Drain `D_UV ≥ 0` | sympy L88,95 / M2-M3 L123,132 | match |
| Finite compiler `T^2`, `eps_eta`, `R_target` (from selected-branch identity) | sympy L121-124 / M4 L147-149 | match |
| Dependent correction `y_nr^dep` | sympy L140 / M5 L160 | match |
| First-order packet `Xi_1`, `R_1`, `R_1+Xi_1` | sympy L166-169 / M6 L187-190 | match |
| Support-blind split `∂_Λ=∂_ρ=0` | sympy L192-193 / M7 L207-208 (+ positive control L202-203 / L214-215) | match |
| Session-I readback | sympy L241-245 / M8 L230-232 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 69-72 | `simplify(...)==0` (U,V,ratio,detH from solve) | stationary packet | yes |
| A2 | sympy | 73-76 | subs-limit checks | limiting slices | yes |
| A3 | sympy | 88,95 | drain closed form + `float>0` | positive drain | yes |
| A4 | sympy | 121-124 | finite compilers + identity re-derivation | finite compiler | yes |
| A5 | sympy | 140 | matrix residual==0 | dependent correction | yes |
| A6 | sympy | 166-169 | series-coeff checks | first-order packet | yes |
| A7 | sympy | 192-193 | `dLam==0`, `dvarrho==0` (trivial: vars absent) | support-blind | yes (claim is absence) |
| A8 | sympy | 202-203 | `dU_bad_dLam != 0` + exact value (Lam present) | positive control for A7 | yes (load-bearing) |
| A9 | sympy | 241-245 | abs-diff vs recorded Session-I numbers | Session-I readback | yes |
| M1-M8 | math | 109-232 | `expectZero/expectTrue/expectApprox` | same eight | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. Evidence:
1. **Hessian/determinant route differs.** SymPy builds the Hessian by autodiff `sp.hessian(F_nr,(U,V))` then `.det()` (L55-56); Mathematica hand-writes `hessian = {{aU,-chiUV},{-chiUV,aV}}` and `Det[hessian]` (L100,112). Different construction of the same admissibility object.
2. **Stationary-solve form differs.** Mathematica's `Solve` returned `uVar -> (fU - chiUV^2 fU/(-(aU aV)+chiUV^2))/aU` (output L10) — a structurally different closed form than sympy's `a_V f_U/(a_U a_V - chi_UV^2)` (output L9) — then proven equal to the canonical form via `expectZero` (L109). A line-by-line port would have produced the same form.
3. **Series machinery differs.** SymPy uses a custom `coeff_linear` (`series().removeO().coeff(eps,1)`, L28-29); Mathematica uses `Coefficient[Normal[Series[...]],eps]` (L164-177).
The two scripts share the natural M1-M8 / section-1-8 outline because they audit one stage, but each derives the results with its own primitives. Not `mathematica_transliteration`.

## Engine cross-check

Both engines agree at the level claimed. Symbolic forms match (sympy output L9-11 vs math output L11-14: `U_sol = a_V f_U/Delta`, `V_sol = chi_UV f_U/Delta`, `V/U = chi_UV/a_V`). Drain at the opposite-sign point agrees to 16+ digits: sympy `0.003937884170811954` (output L23) vs Mathematica `0.003937884170811952585…` (output L42). Session-I: sympy `eps_eta=0.28933482012…, R_ratio=0.87984149193…, R1=-0.12762119` (output L79,80,93) vs Mathematica `0.28933482012…, 0.87984149193…, -0.12762119` (output L122,123,124). All `PASS`. `engines_agree: true`.

## Verdict justification

`clean`. Attacks tried that failed: (a) **F1 self-test trap** — I checked every `sp.diff`/`D[]`: the M7/Section-6 support-blindness derivatives are indeed trivially zero because `Lam`/`varrho` are absent from the differentiated expressions, BUT this is the intended claim (blindness = absence) and it is guarded by the pass-1 positive control (`f_U_bad = f_U + Lam`, sympy L198-203 / Mathematica L212-215) which differentiates an expression where `Lam` genuinely appears, asserts the result is nonzero AND equals the exact value `a_V/Delta_UV`; the fix is present and load-bearing in both engines. (b) **Tautology check on R_target** — the `R_ratio` is not asserted against itself: it is independently re-derived from the selected-branch identity `R_target T^2 = Lambda_0(1-eps_eta)` (sympy L108-110 / Mathematica L138-140) and equated to the paper closed form (`expectZero` passes), so the identity is genuinely exercised. (c) **Session-I round-trip tautology** — the script explicitly avoids asserting the algebraically-exact inverse round-trip (comment L242-243) and instead checks against independently-recorded numbers `0.28933482/0.87984149/-0.12762119` with `5e-9` tolerance — non-tautological. (d) **Symbol domains** — `a_U,a_V>0`, `0<eps_ref<1` in both engines match the physical setup; no positivity is smuggled into a step that needs the opposite. (e) **Drain positivity** is checked on the harder opposite-sign branch (`chi_UV<0`), the Session-I branch, not just the easy same-sign case. Paper card, notes, and appendix row 88 all match the script's verified claims; both engines present and independent; outputs fresh.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `U = a_V f_U/Delta_UV` | py L52 / out L9; wl L109 / out L12 | tex L213 (uv-stationary); md L220 | MATCH |
| `V = chi_UV f_U/Delta_UV` | py L53 / out L10; wl L110 / out L13 | tex L215; md L221 | MATCH |
| `Delta_UV = a_U a_V - chi_UV^2` | py L51 / out L8; wl L86 / out L11 | tex L208; md L212 | MATCH |
| `V/U = chi_UV/a_V` | py L54 / out L11; wl L111 | md L230 | MATCH |
| `D_UV = chi_UV^2 a_V f_U^2/Delta_UV^2` | py L83 / out L21; wl L121 / out L35 | tex L222; md L356 | MATCH |
| `T^2 = T_ref^2 exp(a_V f_U/Delta)` | py L104 / out L28; wl L136 / out L49 | tex L37-41; md L266-272 | MATCH |
| `eps_eta = eps_ref exp(chi_UV f_U/Delta)` | py L105 / out L29; wl L137 / out L50 | tex L43-47; md L275-281 | MATCH |
| `R_target/R_ref = [(1-eps_ref e^V)/(1-eps_ref)] e^{-U}` | py L111 / out L31; wl L141 / out L51 | tex L49-54; md L285-294 | MATCH |
| `y_nr^dep = (0, -chi_UV f_U/Delta, (a_V-chi_UV)f_U/Delta)` | py L131 / out L37-48; wl L154 / out L62 | tex L64-72; md L333-342 | MATCH |
| `Xi_1^nr = u_1` | py L156 / out L53; wl L179 / out L70 | tex L84; md L409 | MATCH |
| `R_1^nr = -u_1 - eps_ref v_1/(1-eps_ref)` | py L159 / out L55; wl L181 / out L72 | tex L91-95; md L422-427 | MATCH |
| `R_1+Xi_1 = -eps_ref v_1/(1-eps_ref)` | py L160 / out L54; wl L180 / out L71 | tex L84-88; md L415-417 | MATCH |
| Session-I `eps_eta = 0.28933482` | py L241 / out L79; wl L230 / out L122 | md L513 (`≈0.28933482`) | MATCH |
| Session-I `R_ratio = 0.87984149` | py L244 / out L80; wl L231 / out L123 | md L523 (`≈0.87984149`) | MATCH |
| Session-I `R_1 = -0.12762119` | py L245 / out L93; wl L232 / out L124 | md L548 (`≈-0.12762119`) | MATCH |
| Session-I inputs `U=0.14313458, V=-0.03619791, eps_ref=0.3` | py L209-211; wl L219-221 | md L501-504 | MATCH |
| Session-I `y_dep = (0, 0.03619791, 0.17933249)` | py L218 / out L81-86 | md L528-535 | MATCH |

INTERNAL (scaffolding, no finding): drain-positivity probe params `a_U=2.5,a_V=3.0,chi_UV=-0.76,f_U=0.33`; inferred-only Session-I quantities `chi_UV≈-0.7586827, chi_lambda≈-0.7986134, Delta_UV≈6.9244006, f_U≈0.3303737` (printed as round-trip diagnostics, not stage deliverables — the notes §8 deliberately report only U,V,eps_eta,R_ratio,y_dep,Xi_1,R_1); `U_rebuilt`, `V_rebuilt` (round-trip identity checks); positive-control derivative `a_V/Delta_UV`; pass/fail flags, tolerances.

reconciliation: complete; 18 deliverable values checked, 0 misaligned

## Self-test notes

Checked variable-independence (trap): the M7/Section-6 `D/diff` w.r.t. `Lam,varrho` are correctly zero only because those symbols are absent — that is the intended support-blindness claim, and the pass-1 positive control (`f_U+Lam`) is present and load-bearing in both engines, asserting a genuinely nonzero derivative `= a_V/Delta_UV`. Checked symmetry/parity: no unbounded-domain integrals in this stage. Checked trivial-case pre-checks: `f_U=0 → U=V=0` and `chi_UV=0 → U=f_U/a_U,V=0` both reduce as asserted; drain at the concrete opposite-sign point is a positive literal. No directive written (zero findings).
