---
unit_id: 015
batch: I.2
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
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 015 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_015.tex`
- notes: `(none)` (no files matching `notes/stages/moving_throat_pde_stage015_*.md`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 52)
- sympy: `/var/projects/toy_physics/.../scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
  (full path `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.txt`

## What the paper claims

The card (`\paragraph{Output.}`) states verbatim: "Stage~015 exports the parent throat action promotion and the exact quadratic recovery formula \eqref{eq:stage015-keta}." Two distinct deliverables follow:

1. **Parent-action promotion.** The throat is promoted to an autonomous graph field `R(Ω,w,t)` with action `S_Σ[R] = ∫ dt dw dΩ L_Σ` and density (eq:stage015-parent-density) `L_Σ = ½μ_Σ R_t² − ½T_{w,Σ} R_w² − ½T_{Ω,Σ}|∇_Ω R|² − U_Σ`, embedded in `S_total = S_ψ + S_EM + S_Σ`.
2. **Exact quadratic recovery formula** (eq:stage015-keta): `K_η = U_{Σ,RR}(R0,w) − ∂_w(T_{w,Σ,R}(R0,w) R0′) + ½ T_{w,Σ,RR}(R0,w)(R0′)²`, the quadratic stiffness around the background `R0(w)`.

The card also states the audit "carries the boundary densities explicitly and tests both zero and nonzero boundary-discharge probes." The appendix row (line 52) summarizes: "Parent throat-action packet and projection/reduction status boundary," status `\StatusExact{}/\StatusReduced{}`. No `notes/stages/` file exists for this stage; the `.tex` is therefore authoritative.

## What the script claims to verify

Both engines verify the same physics. (1) A generic quadratic integration-by-parts (IBP) product-rule identity `−A η η′ = ∂_w(−Aη²/2) + (A′η²/2)`, plus a sign-mutation that must fail. (2) Boundary-discharge probes: a Gaussian profile where the boundary term and cross/bulk integrals vanish symmetrically, and an asymmetric profile (`A = w e^{−w²}`) where cross and bulk are each individually nonzero (`Sqrt[Pi/2]/4`) yet the IBP cancellation still holds — so the identity is a real cancellation, not `0 = 0 + 0`. A nonzero endpoint probe (`atan(w)`) confirms the boundary operator can actually detect discharge. (3) The `K_η` recovery: SymPy expands the parent Lagrangian `L`, extracts the `eps²` coefficient (`L2_raw`), applies the IBP on the `−T_{w,R0} R0′ η η_w` cross term, and asserts the resulting quadratic action equals the canonical form carrying `K_η = URR0 − d_TwR_R0p + TwRR0 R0p²/2`; a sign-flipped `K_η` must fail. Mathematica derives `K_η` by an **independent Euler–Lagrange linearization** of the same density and asserts the same closed form, plus its own sign mutation.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Parent-action density `L_Σ` (eq:stage015-parent-density) | SymPy `L = mu0*Rt²/2 − Tw*Rw²/2 − TO0*eps²*grad2/2 − U` (py:93); Mathematica `LDensity` (wl:111-112) | match (structure exercised: it is the source of the quadratic limit) |
| `K_η` recovery formula (eq:stage015-keta) | SymPy `K_eta = URR0 − d_TwR_R0p + TwRR0*R0p²/2` asserted equal to IBP'd quadratic action (py:95-104); Mathematica EL route asserts same form (wl:135-141) | match |
| Boundary densities carried; zero + nonzero discharge probes tested | SymPy `nonzero_boundary_probe` (py:40-41), `quad_boundary_concrete` zero (py:51), asymmetric probe (py:60-79); Mathematica M1/M2 (wl:60-100) | match |
| `S_total = S_ψ + S_EM + S_Σ` (eq:stage015-total-action) | Not separately checked | extra/none — bookkeeping statement, no independent identity to verify; not a deficiency |

`paper_alignment: aligned`. The only un-exercised paper line is the additive total-action bookkeeping (eq:stage015-total-action), which states no algebraic identity for a script to test — it is a definitional sum, not a claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 30-33 | `assert_zero(−Aηη′ − [∂_w(−Aη²/2)+A′η²/2])` | boundary/IBP machinery (deliverable 2 sub-step) | yes |
| A2 | sympy | 34-37 | `assert_nonzero(mutated boundary sign)` | guards A1 | yes |
| A3 | sympy | 41 | `assert_nonzero(boundary_value(atan))` | nonzero-discharge probe | yes |
| A4 | sympy | 51 | `assert_zero(Gaussian boundary discharge)` | zero-discharge probe | yes |
| A5 | sympy | 52-55 | `assert_zero(cross − (boundary+bulk))` | concrete IBP w/ boundary | yes |
| A6 | sympy | 73-74 | `assert_nonzero(asym cross / bulk)` | non-tautology guard for A7 | yes |
| A7 | sympy | 76-79 | `assert_zero(asym cross − (boundary+bulk))` | real-cancellation IBP | yes |
| A8 | sympy | 100 | `assert_zero(raw_cross_coeff + TwR0*R0p)` | isolates the cross term IBP'd | yes |
| A9 | sympy | 104 | `assert_zero(L2_after_ibp − canonical_L2)` | **K_η recovery (deliverable 2)** | yes |
| A10 | sympy | 108 | `assert_nonzero(K_η sign-mutation fails)` | guards A9 | yes |
| M1 | mathematica | 60 | `expectZero(IBP product-rule)` | mirror of A1 | yes |
| M2 | mathematica | 65 | `expectNonzero(mutated sign)` | mirror of A2 | yes |
| M3 | mathematica | 75-100 | Gaussian + asymmetric discharge probes | mirror of A4-A7 | yes |
| M4 | mathematica | 140-141 | `expectZero(K_η via EL − (URR0−dTwRR0p+TwRR0 R0p²/2))` | **K_η recovery, independent EL route** | yes |
| M5 | mathematica | 144-145 | `expectNonzero(K_η EL sign-mutation)` | guards M4 | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration. The boundary/IBP probes (M1/M2) do mirror the SymPy probes line-for-line — appropriate, since these are concrete-integral cross-checks where echoing is benign and is the point of a second engine. The load-bearing `K_η` derivation, however, is genuinely independent:

- SymPy (py:94-104) extracts `L2_raw = ∂²_eps L / 2 |_{eps=0}`, then performs the IBP **by hand** (`L2_raw − cross_term + cross_after_ibp`) to reach the canonical quadratic action — a *second-variation-then-IBP* route.
- Mathematica (wl:114-141) instead forms the **Euler–Lagrange operator** `∂L/∂R − D_t(∂L/∂R_t) − D_w(∂L/∂R_w)` on slot variables, substitutes the perturbed field, takes the `eps¹` coefficient, reads off the mass term, and resolves the product derivative via `R0p*twR′[w] → dTwRR0p − TwR0*R0pp` (wl:137) — an *EL-linearization* route that surfaces and cancels the `R0pp` (`R0″`) terms that the product rule `∂_w(T_{w,R}R0′)` generates. Two different formalisms reaching `URR0 − dTwRR0p + TwRR0 R0p²/2`. This is a true second-engine confirmation, not algebra echo.

## Engine cross-check

Both engines exit `STATUS: PASS`. The `K_η` closed form agrees exactly (SymPy `URR0 − d_TwR_R0p + TwRR0*R0p²/2` ≡ Mathematica `URR0 − dTwRR0p + TwRR0*R0p²/2`). The asymmetric IBP cross and bulk integrals agree numerically/symbolically across engines: Mathematica reports `Sqrt[Pi/2]/4` for both cross and bulk (output lines 10, 12) and their difference is 0; SymPy asserts the same nonzero-then-cancel structure (A6/A7). Mutation residuals are nonzero in both (`−2*dTwRR0p` Mathematica, sign-flip in SymPy), confirming the guards bite. No disagreement.

## Verdict justification

`clean`. I read the card, confirmed there is no `notes/stages/` file (the `.tex` is authoritative), and read the appendix row. Both engines are present, fresh (outputs dated 22:02/22:06 vs scripts 22:00/22:05), and agree. Attacks tried that failed to break it: (1) the `K_η` check is **not** tautological — `canonical_L2` is built from the IBP-transformed `L2_raw`, and the sign-mutation guard (A10/M5) returns nonzero, so the assertion genuinely constrains the formula; (2) the IBP probes are **not** `0=0+0` — the asymmetric probe forces cross and bulk individually nonzero (`Sqrt[Pi/2]/4`) before they cancel, exactly as the card requires ("both zero and nonzero boundary-discharge probes"); (3) parity/derivative traps checked (see self-test) — the Mathematica `D[dLdRt, t]` term legitimately drops (it lands on the `η_t` monomial, not the `η` mass coefficient that defines `K_η`), so it does not corrupt the extracted result; (4) the Mathematica route is independent (EL linearization vs SymPy second-variation+IBP), satisfying the second-engine policy. The script's verified claim matches the paper's `Output` deliverables exactly.

Non-finding note (informational): the SymPy docstring (py:2) and print labels (py:110, "STEP 13 PARENT THROAT ACTION MASTER AUDIT") use the stale pre-renumber label "step 13 / STEP 13" rather than "stage 015." This is the documented incomplete-renumber from the EM-extension realignment; the canonical filename and card are stage 015. It is a cosmetic prose label, carries no numeric/mathematical deliverable, and mechanical label-only numbering work is documented as exhausted/deferred — so it is not raised as a finding here.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `K_η = URR0 − d_TwR_R0p + TwRR0·R0p²/2` (closed-form quadratic stiffness) | py:95 + asserted py:104; wl:141 (`URR0 − dTwRR0p + TwRR0*R0p²/2`); wl output line 18-19 PASS | `stage_015.tex:35-40` eq:stage015-keta: `K_eta = U_{Σ,RR}(R0,w) − ∂_w(T_{w,Σ,R}(R0,w)R0′) + ½ T_{w,Σ,RR}(R0,w)(R0′)²` | MATCH |
| Canonical quadratic Lagrangian `½μ0 η_t² − ½Tw0 η_w² − ½TO0 grad2 − ½K_η η²` (`canonical_L2`) | py:96; asserted py:104 | `stage_015.tex:18-25` eq:stage015-parent-density (`L_Σ = ½μ_Σ R_t² − ½T_{w,Σ}R_w² − ½T_{Ω,Σ}|∇_Ω R|² − U_Σ`); quadratic limit reported via eq:stage015-keta | MATCH (the quadratic limit of the carded `L_Σ`) |
| Parent action `S_Σ[R] = ∫ dt dw dΩ L_Σ` | py:93 (`L` density form); wl:111-112 (`LDensity`) | `stage_015.tex:14-16` eq:stage008-parent-action + density eq:stage015-parent-density | MATCH |
| Raw cross-term coefficient `−TwR0·R0p` | py:100 (`raw_cross_coeff + TwR0*R0p` → 0) | (implicit in eq:stage015-keta product-derivative term `∂_w(T_{w,Σ,R}R0′)`) | MATCH (intermediate of the carded K_η) |

INTERNAL scaffolding (accounted-for, no prose expected, no finding): generic IBP product-rule residual `0`; mutated-boundary-sign residual (`−η(ηA′+2Aη′)` Mathematica) — guard; `atan` boundary probe value `π` — sanity probe; Gaussian boundary discharge `0` and cross−bulk `0` — symmetric probe; asymmetric cross/bulk integral value `Sqrt[Pi/2]/4` — non-tautology demonstrator (intentionally nonzero, then cancels); asymmetric boundary discharge `0`; `K_η` sign-mutation residuals (`−2*dTwRR0p` Mathematica) — guards. None of these are stated deliverables; all are verification machinery.

reconciliation: complete; 4 deliverable values checked, 0 misaligned.

## Self-test notes

Checked the standard traps. (1) Variable-independence / derivative-zero trap: the Mathematica `D[dLdRtSlot, t]` term (wl:128) is identically zero because `etat` is a plain symbol, but this is harmless — it would contribute only to the `η_t` monomial, whereas `K_η` is read off `Coefficient[ELOrderEps, eta]` (the `η` mass term), so the dropped term cannot make any `K_η` assertion pass trivially; the sign-mutation guards (A10/M5) independently confirm the `K_η` assertions are live. (2) Parity/symmetry: the asymmetric probe `A = w e^{−w²}` (odd) × `η² = e^{−w²}` (even) gives an even integrand, so the cross/bulk integrals can be — and are (`Sqrt[Pi/2]/4`) — nonzero on the symmetric domain, exactly as the script's own comment (py:56-59, wl:83-84) reasons; the IBP cancellation is therefore a genuine cancellation. (3) Trivial-case: substituting the Gaussian (even `A`) makes cross and bulk vanish symmetrically (output 0), and the boundary terms vanish at ±∞, all consistent. No directive written (zero findings).
