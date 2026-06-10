---
unit_id: 223
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 223 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_223.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at L58, L100, L574-597 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage223_5pn_isotropic_target_surface_primitive_branch_compatibility_and_dynamic_survival_window_mathematica_audit.txt`

## What the paper claims

Stage 223 asks whether the same primitive finite-throat one-port branch can simultaneously satisfy the conservative isotropic one-pole condition (`u4=4u2^2`) and the outgoing-normalization target (`P0=P0_target`). The card's `\stagefield{Output}` reads verbatim: *"Exact primitive compatibility surface $N_0/P_{0,{\rm target}}=3(M+B_2+Z_2)^2/(B_4+Z_4)$ and finite dynamic survival windows on the compatible branch."* The derivation eliminates the static stiffness `K` through the one-pole condition (`K_pole`) and through the normalization condition (`K_norm`), then equates them. The notes enumerate nine deliverables (notes §9): the overlap constant `κ=2√2/π`, the one-pole numerator identity, the exact compatibility identity, its primitive specialization, the concrete sample values (`P0_target_compat≈0.00207`, `K_compat≈24.4738`, `D0_compat≈24.2373`), the quartic pole census (4 roots), the wall-like residue/linewidth `R_Q` figures, the monotone `λ_W` compatibility-family scan, and the four finite survival-window thresholds in `P0_target_compat`. The appendix (L574-597) reproduces the isotropic D-coefficients and the boxed compatibility equation. Card status is "Mixed: ExactClosure + Numerical."

## What the script claims to verify

The SymPy script verifies, in order: the symbolic overlap integral equals `2√2/π` (L92-93); the one-pole numerator identity `u4-4u2^2 = (D0(B4+Z4)-3(M+B2+Z2)^2)/D0^2` (L149-152); the compatibility surface by SOLVING the stiffness-equality `K_norm==K_pole` for `P0_target` and matching it to the boxed `N0(B4+Z4)/(3(M+B2+Z2)^2)` (L154-162); the primitive specialization (L164-167); the concrete sample numerics via `assert_close` against pinned constants (L240-255); the quartic pole census roots and `R_Q` values (L264-287); the monotone `λ_W` scan with both value and monotonicity asserts (L371-385); and the four bisected survival-window ceilings (L390-418). The Mathematica `.wl` mirrors the same nine claims (M1-M9) with native primitives. Every block ends in an `assert`/`expectZero`/`expectClose`/`fail-Exit[1]`; no block is print-only.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| overlap constant `κ=2√2/π` | py L92-93 / wl M1 (Integrate) | match |
| one-pole numerator identity | py L149-152 / wl M2 | match |
| exact compatibility surface (Output box) | py L154-162 (solve) / wl M3 | match |
| primitive specialization of `P0_target_compat` | py L164-167 / wl M4 | match |
| sample `P0_target_compat`,`K_compat`,`D0_compat` | py L253-255 / wl M6 | match |
| quartic pole census (4 roots) | py L264-271 / wl M5,M7 | match |
| wall-like `R_Q` figures | py L280-287 / wl M8 | match |
| monotone `λ_W` scan | py L371-385 / wl M9 | match |
| 4 finite survival windows | py L415-418 / wl M9 | match |

`paper_alignment: aligned`. Every stated deliverable has a load-bearing dual-engine check; nothing extra, nothing missing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 93 | `simplify(κ - 2√2/π)==0` | overlap constant | yes |
| A2 | sympy | 152 | `cancel(u4-4u2^2 - (...)/D0^2)==0` | one-pole identity | yes |
| A3 | sympy | 156,160-161 | `solve(Eq(K_norm,K_pole),P0_target)[0] - boxed ==0` | compatibility surface (Output) | yes |
| A4 | sympy | 162 | `cancel(D0_compat - N0/P0_target_compat)==0` | compatibility surface | yes |
| A5 | sympy | 167 | `cancel(P0_target_compat - primitive_specialization)==0` | primitive specialization | yes |
| A6 | sympy | 182 | `Poly(F_y,y).degree()==4` | quartic pole structure | yes |
| A7 | sympy | 240-255 | `assert_close` on 17 sample constants | sample numerics | yes |
| A8 | sympy | 270-271 | `assert_close` on 4 pole roots | pole census | yes |
| A9 | sympy | 286-287 | `assert_close` on 4 `R_Q` | wall residues | yes |
| A10 | sympy | 378-385 | `assert_close` scan + monotonicity | `λ_W` scan | yes |
| A11 | sympy | 415-418 | `assert_close` 4 ceilings | survival windows | yes |
| M1 | mathematica | 137 | `expectZero[κ - 2√2/π]` | overlap constant | yes |
| M2 | mathematica | 151 | `expectZero[onePoleResidual]` | one-pole identity | yes |
| M3a | mathematica | 160 | `expectZero[solvedTarget - boxedTarget]` | compatibility surface | yes |
| M3b | mathematica | 161-164 | `expectZero[nZero/boxedTarget - 3(...)^2/(b4+z4)]` | compatibility surface | yes |
| M4 | mathematica | 202 | `expectZero[pCompatPrimitive - primitiveClaim]` | primitive specialization | yes |
| M5 | mathematica | 217-218 | `Exponent==4` + `Last[coeffs]-mass==0` | quartic structure | yes |
| M6 | mathematica | 247-249 | `expectClose` 3 sample values | sample numerics | yes |
| M7 | mathematica | 262-268 | `expectVectorClose` roots + ordering | pole census | yes |
| M8 | mathematica | 280-285 | `expectVectorClose` 4 `R_Q` | wall residues | yes |
| M9 | mathematica | 306-352 | `expectVectorClose` scan + 4 ceilings + thresholds | scan + survival windows | yes |

No tautological or orphaned rows.

## Findings

None. Zero findings.

## Independent-derivation check (Mathematica)

**INDEPENDENCE CALL: independent.**

The `.wl` mirrors the SymPy *claim sequence* (M1-M9 follow the same nine deliverables), which is expected and acceptable — but the derivation route uses native Mathematica primitives throughout, and in several places does strictly more than the SymPy script:

1. **Overlap integral (M1 vs py L92).** SymPy: `sp.integrate(u0*f0, (s,0,L))` then `simplify(... - 2√2/π)`. Mathematica: `FullSimplify[Integrate[uConst fProfile, {s,0,len}], Assumptions -> len>0]` then `expectZero[... - 2√2/π]`. Both perform the genuine symbolic integral from the wavefunction profiles, not a transcribed closed form. Independent.

2. **Compatibility surface (M3 vs py L156-161, the F2-fix block).** SymPy: `P0_target_compat_solved = sp.solve(sp.Eq(K_norm, K_pole), P0_target)[0]`, then asserts `cancel(solved - boxed)==0`. Mathematica: `solvedTarget = stripConditional[p0Target /. First[Solve[kNormIso == kPoleIso, p0Target]]]`, then `expectZero[solvedTarget - boxedTarget]`. Both engines independently SOLVE the stiffness-equality equation rather than echo a pre-computed form. Mathematica ADDS a second independent check M3b (`nZero/boxedTarget - 3(...)^2/(b4+z4)`) absent from SymPy.

3. **Quartic structure (M5 vs py L182).** SymPy checks only `Poly(F_y,y).degree()==4`. Mathematica uses a different decomposition — `CoefficientList[fyExpanded, y]` with `Exponent[...]==4` AND additionally `Last[fyCoefficients] - mass == 0` (leading coefficient equals the mass), a check SymPy does not perform.

4. **Root finding (M7 vs py L262).** SymPy: NumPy `np.roots` on float coefficients. Mathematica: `NSolve[poly==0, y, WorkingPrecision->60]` at 60-digit working precision with an independent positive-real selection. Genuinely different numeric kernels; the agreement is meaningful corroboration, not a port.

The single shared "operation" is `Solve`/`solve` applied to the same stiffness-equality equation and the same symbolic moment definitions (`B0..Z4`, `N0`). Those moment definitions are the *physical premises of the stage*, not borrowed intermediate algebra; both engines must start from the same premises by definition. No SymPy-specific choreography (intermediate substitution order, helper-function naming, lambdify scaffolding) is echoed. **Not a transliteration.**

## Engine cross-check

Both engines pass all checks (SymPy "audit completed successfully", Mathematica "All M1-M9 completed successfully"). Numeric agreement at the level claimed:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `P0_target_compat` (sample) | 0.002069792318062883 | 0.00206979231806288269... |
| `K_compat` (sample) | 24.473754879290976 | 24.47375487929097685925... |
| `D0_compat` (sample) | 24.23730998862225 | 24.23730998862225050219... |
| pole roots | 0.97157, 1.41651, 1.99754, 4.94905 | 0.97157, 1.41651, 1.99754, 4.94905 |
| `R_Q` (4 wall/internal) | 0.15989, 0.000806, 30.1999, 36.1712 | 0.15989, 0.000806, 30.1999, 36.1712 |
| 4 survival ceilings | 0.0028313, 0.0035965, 0.0081734, 0.0116634 | 0.0028313, 0.0035965, 0.0081734, 0.0116634 |

All residuals reported by the `.wl` are ≤ ~1e-13, inside the declared tolerances. Engines agree.

### ⚠️ R_Q reconfirmation (λ_W=0.2 wall values)

- **SymPy** emits (output L55): `(0.2, 0.0005769708798430451, 29.31584648723137, 138.81413694208146, 137.50254660071312)` — i.e. lower/upper wall `R_Q = 138.814136942081 / 137.502546600713`. **CONFIRMED.**
- **Mathematica** `.wl`: the M9 scan begins at `λ_W=2/5` (`scanLambdas = {2/5,3/5,4/5,1}`, wl L301), so it does NOT emit a λ_W=0.2 scan row. The `.wl` therefore does not independently reproduce 138.81/137.50; it corroborates the λ_W=0.4-1.0 rows (and the λ_W=0.4 row's wall `R_Q` 30.1999/36.1712 matches the main M8 census). Cross-engine corroboration of the corrected values is via the matching scan formula at λ_W=0.4 and the shared `rqOmega`/`rq_func` definition, not a literal λ_W=0.2 row.
- **Notes** line 351 reads: `| 0.2 | 0.000576970879843 | 29.3158464872314 | 138.814136942081 | 137.502546600713 |` — matches the SymPy-emitted corrected values. **CONFIRMED.**
- **Surviving 206/205?** A whole-repo grep for `206.814136942081` / `205.502546600713` returns hits ONLY in red-team scaffolding (the first-pass `redteam/directives/stage_223.md`, `redteam/reports/stage_223.md`, `redteam/codex_logs/*`, `redteam/batches/*`, and historical tracker snapshots in `notes/MATHEMATICA_MIRROR_POLICY.md`/`CHECKPOINT_CONSTANT_PROVENANCE.md` recording the resolved typo). **NONE in the live stage notes, the `.tex` card, the appendix, or either script/output.** No surviving 206/205 in any deliverable-carrying document.

### F2 compat-surface de-tautologization

**Holds.** The block (py L154-162) computes `K_pole = 3(M+B2+Z2)^2/(B4+Z4)+B0+Z0` and `K_norm = N0/P0_target+B0+Z0` from independent premises (one-pole side vs normalization side), then `sp.solve(sp.Eq(K_norm,K_pole), P0_target)[0]` solves the equality for `P0_target` and asserts it equals the independently-constructed boxed target `N0(B4+Z4)/(3(M+B2+Z2)^2)`. This is load-bearing: `K_pole` and `K_norm` are built from disjoint constructions (the one-pole numerator structure vs the outgoing-normalization ratio), so their equality solving back to the boxed form is a genuine consistency theorem, not `x==x`. It would fail if either the one-pole `K` elimination or the normalization `K` elimination were wrong. The `.wl` M3 reproduces this with `Solve` and adds the redundant M3b form. The compat-surface check is non-tautological and confirmed in both engines.

## Verdict justification

`clean`. I attacked: (1) the F2 de-taut fix — confirmed `K_pole` and `K_norm` are independently constructed so the `solve`-then-match is a real theorem, not circular; (2) the `.wl` for transliteration — found native primitives (Integrate, Solve, NSolve@60-digit, CoefficientList) plus two checks (M3b, leading-coeff=mass) SymPy lacks, so it is independent; (3) the λ_W=0.2 R_Q reconfirmation — SymPy and the notes both emit the corrected 138.814.../137.502..., no 206/205 survives in any deliverable doc; (4) symbol domains — all `positive=True` declarations (κ, λ's, Ω's, ϖ, K, M, P0_target) match the physical setup and the `.wl` `$Assumptions` mirror them, with the added non-degeneracy guards `kWall-b0-z0≠0`, `b4+z4≠0`, `mass+b2+z2≠0` that protect the divisions; (5) the survival-window bisections — brackets straddle (the `assert_close` ceilings match between engines to ~1e-13). Both engines exit 0, outputs are fresh (16:37 > scripts 16:23/16:26). I read the card, notes, and appendix row, and the script's verified claim matches the paper's `\stagefield{Output}` exactly. The only non-blocking item is the known-deferred numbering class noted below.

NOTE (deferred numbering class, NOT a finding): the system reminder flagged that the first pass left 223's notes TITLE reading "Stage 240". On this read the notes title (L2) reads "Stage 223" correctly — the title appears to have been corrected, and a grep of the notes file for "240" returns no hits. No surviving "Stage 240" reference was found in the stage's notes. Either way this falls under the deferred project-wide stem-keyed renumber pass (P4-53 residual class), not a pass-2 blocker.

## Self-test notes

Checked: (1) variable independence — no `sp.diff`/`D[]` derivative-style assert in this stage (the checks are algebraic identities and numeric round-trips), so the zero-derivative trap is not applicable. (2) Symmetry/parity — the only integral is the overlap `∫_0^L (1/√L)·√(2/L)·sin(πs/2L) ds` over a bounded [0,L] domain with a single-signed positive integrand, so it is legitimately nonzero (=2√2/π); no symmetric-domain cancellation risk. (3) Trivial-case pre-check — the sample slice reduces `P0_target_compat`, `K_compat`, `D0_compat` to the pinned numerics that both engines reproduce to ~1e-16/1e-21, and the F2 `solve` returns a single positive root (P0_target declared positive), so the `[0]`/`First` selection is well-defined. No traps fired.

## Value Reconciliation (pass-2 augmentation)

Both engines and both saved outputs were used as the authoritative record (outputs fresh; nothing run). Every emitted deliverable value reconciles against the `.tex` card and/or the `.md` notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `κ = 2√2/π` | py out L4 / wl out L11 | notes L47 `κ=2√2/π`; appendix narrative | MATCH |
| `P0_target_compat ≈ 0.00206979231806289` (sample) | py out L39 / wl out L46 | notes L254 `≈0.00206979231806289` | MATCH |
| `K_compat ≈ 24.4737548792910` (sample) | py out L40 / wl out L47 | notes L260 `≈24.4737548792910` | MATCH |
| `D0_compat ≈ 24.2373099886222` (sample) | py out L41 / wl out L48 | notes L268 `≈24.2373099886222` | MATCH |
| pole roots 0.971575315129468, 1.41651290122561, 1.99753567893361, 4.94905432364313 | py out L44-47 / wl out L59 | notes L289-298 (all four) | MATCH |
| `R_Q` 0.159888393135835, 0.000806281535937178, 30.1999075602499, 36.1711864832695 | py out L44-47 / wl out L70 | notes L309-318 (all four) | MATCH |
| threshold η=0.1 ≈ 21.8545662963584 | py out L50 / wl out L90 | notes L340 `≈21.8545662963584` | MATCH |
| threshold η=0.3 ≈ 7.86187368416853 | py out L51 / wl out L92 | notes L343 `≈7.86187368416853` | MATCH |
| scan λ_W=0.2 row: 0.000576970879843, 29.3158464872314, **138.814136942081**, **137.502546600713** | py out L55 (wl scan starts at 0.4) | notes L351 (matches, corrected values) | MATCH |
| scan λ_W=0.4 row: 0.002069792318063, 24.4737548792910, 30.1999075602499, 36.1711864832695 | py out L56 / wl out L78 | notes L352 | MATCH |
| scan λ_W=0.6 row: 0.004865681200486, 21.1544287401845, 12.8348600273988, 16.7575510327116 | py out L57 / wl out L79 | notes L353 | MATCH |
| scan λ_W=0.8 row: 0.009169913681573, 19.0298300900561, 7.06074242207991, 9.69035785242054 | py out L58 / wl out L80 | notes L354 | MATCH |
| scan λ_W=1.0 row: 0.014981190324091, 17.7824591822917, 4.45922850098679, 6.30111094469551 | py out L59 / wl out L81 | notes L355 | MATCH |
| 10% lower-wall ceiling ≈ 0.00283133168555932 | py out L64 / wl out L94 | notes L373 `≲0.00283133168555932` | MATCH |
| 10% upper-wall ceiling ≈ 0.00359651058968466 | py out L66 / wl out L95 | notes L378 `≲0.00359651058968466` | MATCH |
| 30% lower-wall ceiling ≈ 0.00817339430971383 | py out L69 / wl out L96 | notes L386 `≲0.00817339430971383` | MATCH |
| 30% upper-wall ceiling ≈ 0.0116633929790174 | py out L71 / wl out L97 | notes L390 `≲0.0116633929790174` | MATCH |
| compatibility surface `N0/P0_target = 3(M+B2+Z2)^2/(B4+Z4)` (symbolic) | py L160-161 / wl M3 | `.tex` Output L15; appendix eq L585-589; notes L155-159 | MATCH |
| compatible target `P0_target_compat = N0(B4+Z4)/(3(M+B2+Z2)^2)` (symbolic) | py L157 / wl L158 | appendix eq L591-596; notes L163-167 | MATCH |

INTERNAL scaffolding (accounted-for, no finding): sample intermediate constants `C≈0.450158, G_U=0.3, G_W≈0.360127, R≈0.225079, Δ≈1.909339, Q≈0.354725, P≈0.427650, B0,B2,B4,Z0,Z2,Z4,N0` (these also appear in notes §4 and reconcile, but are intermediate inputs to the deliverables, not bottom-line outputs); the symbolic `F(y)` quartic; the bisection root λ_W values; pass/fail flags; residuals; tolerances; degree/coefficient-count prints.

reconciliation: complete; 19 deliverable values checked, 0 misaligned.
