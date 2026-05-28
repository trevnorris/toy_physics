---
unit_id: 167
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
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
  notes_stage_files: [moving_throat_pde_stage167_bundle_transport_tangent_compensation.md]
  paper_appendix: present
---

# Audit unit 167 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_167.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage167_bundle_transport_tangent_compensation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 65, 334-361, 1463)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Proves arbitrary first-order isotropic bundle drift is tangent to the exact compensated parent family." The appendix row (line 65) restates this and eq:app-part05-tangent-parent-family (lines 351-360) makes the bottom-line concrete: for pure bundle drift `δln r_c = 0`, `δln 𝔯 = 0`, `δln 𝔤 = 0`. The notes derive this from a closed first-order substitution chain: the four bundle observables `(Θ_w, K_s, K_q, P_0)` map to all remaining mouth/background drifts via explicit linear carry-forward laws (notes §1-3, collected in §7), and §5 strengthens the theorem to `δ_⊥ = 0` (the Stage-248/163 off-family normal coordinate vanishes), with the two log-channels `δln(g_q K_s/(g_s λ)) = 0` and `δln(K_s K_q/λ²) = 0`. Distinct deliverables: (D1) the three parent invariants vanish; (D2) the two off-family channels vanish and `δ_⊥ = 0`; (D3) the explicit carry-forward laws for `v_{w0}, T_m, g_s, g_q, λ` (notes §2-3, §7), which are the supporting derivation rather than the headline Output.

## What the script claims to verify

The SymPy docstring (lines 7-11) lists four checks: (1) bundle transport formulas for `c_{s,w}, ℓ, L_W, v_{w0}, T_m`; (2) transport for `g_s, g_q, λ`; (3) first-order invariance of `r_c, 𝔯, 𝔤`; (4) vanishing of the off-family channels and `δ_⊥`. In practice SymPy *asserts* only items (3) and (4) — six `expect_zero` calls at lines 77-79 (invariants) and 84-85, 89 (channels + `δ_⊥`); items (1) and (2) are only printed (lines 44-70, 91-98). The Mathematica script asserts the same invariants/channels AND additionally asserts the individual carry-forward laws for `v_{w0}, T_m, g_s, g_q, λ` (lines 53-54, 68-70) plus two `v/T` sum/difference cross-checks (lines 55-56). Both build all derived quantities from the same upstream chain (`drho, da, dcs, dZ` → `dcsw, dell, dLW` → `dv, dT` → `dgq, dgs, dIsq, dlam`), so the invariant checks genuinely exercise the whole substitution rather than testing definitions against themselves.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1: `δln r_c = δln 𝔯 = δln 𝔤 = 0` | sympy L77-79, wl L76-78 | match |
| D2: two channels `= 0` and `δ_⊥ = 0` | sympy L84-85,89; wl L83-84,87 | match |
| D3: carry-forward laws `v_{w0}, T_m, g_s, g_q, λ` | wl L53-54,68-70 assert; sympy prints only | match (mathematica), partial (sympy) |

The paper's headline Output (D1) and the strengthened theorem (D2) are fully asserted in both engines. The supporting carry-forward laws (D3) are asserted in Mathematica and printed-but-not-asserted in SymPy. Because the SymPy invariant checks are built on `dlam`, `dgs`, `dgq` (which are themselves built on `dv`, `dT`), an error in those carry-forward expressions that did not cancel inside the invariant combinations would still go undetected by SymPy alone — but Mathematica's explicit assertions close that gap. Net `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 77 | `expect_zero(2*dlam - dKs - dKq)` | D1 (r_c) | yes |
| A2 | sympy | 78 | `expect_zero(dlam - (dKs+dKq)/2)` | D1 (𝔯) | yes |
| A3 | sympy | 79 | `expect_zero(dgq + dKs/2 - dgs - dKq/2)` | D1 (𝔤) | yes |
| A4 | sympy | 84 | `expect_zero(dgq + dKs - dgs - dlam)` | D2 (chan1) | yes |
| A5 | sympy | 85 | `expect_zero(dKs + dKq - 2*dlam)` | D2 (chan2) | yes |
| A6 | sympy | 89 | `expect_zero(gstar*chan1 + chan2/(4*sqrt(1+rstar^2)))` | D2 (δ_⊥) | yes |
| A7 | math | 53 | `expectZero(dv - (-3dKs/4+dKq/2+13dTheta/8))` | D3 (v_{w0}) | yes |
| A8 | math | 54 | `expectZero(dT - (-5dKs/4+dKq/2+15dTheta/8-2dP/5))` | D3 (T_m) | yes |
| A9 | math | 55 | `expectZero((dv-dT) - (2dcs - da))` | D3 (v/T) | yes |
| A10 | math | 56 | `expectZero((dv+dT) - (dZ+5drho-4da))` | D3 (vT) | yes |
| A11 | math | 68 | `expectZero(dgs - (-dKs/4+dKq/2+3dTheta/8-2dP/5))` | D3 (g_s) | yes |
| A12 | math | 69 | `expectZero(dgq - (-3dKs/4+dKq+3dTheta/8-2dP/5))` | D3 (g_q) | yes |
| A13 | math | 70 | `expectZero(dlam - (dKs+dKq)/2)` | D3 (λ) | yes |
| A14 | math | 76-78 | `expectZero(drc), (dr), (dg)` | D1 | yes |
| A15 | math | 83-84,87 | `expectZero(chan1),(chan2),(deltaPerp)` | D2 | yes |

All asserted rows are non-tautological: every quantity flows through the full upstream substitution; none is "define `x = e` then assert `x == e`."

## Findings

### F1 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py:29`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl:26`

**What's wrong:**
Both scripts print a banner reading `"STAGE 150 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION"` (sympy line 29, wl line 26), and this propagates verbatim into both saved transcripts (sympy output line 11, mathematica output line 11). This is audit unit 167, not 150. The label is a factual misstatement of which stage the transcript belongs to. (The transcript headers correctly name the stage-167 script files, so the body banner directly contradicts the header.) This is not a math defect — the title string is non-load-bearing — but it is a provenance error that would mislabel the captured evidence for this stage. Note: it is categorized `hardcoded_result` only because it is a baked-in literal that the check confirms nothing about; it is closest in spirit to a stale label.

Separately (informational, not a finding): the SymPy inline comments reference "Stage 166 bundle inversion" (L33), "Stage 165 exact branch drifts" (L52), "Stage 163 off-family channels" (L81), while the notes cite Stages 251/250/248 for the same identities and the appendix prose cites "Stages 165--167". This is renumbering churn between note drafts and is cosmetic; it does not affect the asserted math and is left unflagged.

**Why this matters:**
A future reader auditing the stage-167 transcript would see "STAGE 150" and either distrust the file or mis-file the evidence. It is the kind of copy-paste artifact that erodes trust in the audit corpus even though the underlying algebra is correct.

**Required change:**
Change the banner string at sympy L29 and wl L26 from `STAGE 150` to `STAGE 167`. No other edit. Re-run both engines so the transcripts capture the corrected banner.

**Verification:**
After the fix, sympy output line 11 and mathematica output line 11 both read `STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION`, both scripts still exit 0, and all 6 (sympy) / 15 (mathematica) checks still PASS.

## Independent-derivation check (Mathematica)

The `.wl` is structurally a near port of the `.py`: identical symbol names (`drho, da, dcs, dZ, dcsw, dell, dLW, dv, dT, dgq, dgs, dIsq, dlam, drc, dr, dg, gstar, rstar`), identical definitions in the same order (sympy L34-65 ↔ wl L31-61), and the same `expect_zero`/`expectZero` harness. This borders on `mathematica_transliteration`. However, two factors keep it below the threshold for a finding: (1) the underlying physics is a finite set of linear log-derivative substitutions — there is essentially no degree of freedom for an "independent re-derivation" of a linear identity beyond re-stating the same coefficients; and (2) the `.wl` is strictly *more* thorough than the `.py`, adding seven assertions the SymPy script never makes (the explicit carry-forward laws A7-A13), which is the opposite of a passive echo. The second engine therefore adds genuine verification coverage rather than merely re-running the first engine's algebra, so I do not file a transliteration finding.

## Engine cross-check

Both engines produce identical derived expressions and all checks pass. Side-by-side of the load-bearing outputs:

- `δln v_{w0}`: sympy `dKq/2 - 3*dKs/4 + 13*dTheta/8` (out L21) = math `(4*dKq - 6*dKs + 13*dTheta)/8` (out L24) — identical.
- `δln T_m`: sympy `dKq/2 - 5*dKs/4 - 2*dP/5 + 15*dTheta/8` (out L22) = math `(20*dKq - 50*dKs - 16*dP + 75*dTheta)/40` (out L25) — identical.
- `δln g_q, g_s, λ`, `r_c, 𝔯, 𝔤 = 0`, both channels `= 0`, `δ_⊥ = 0`: identical in both transcripts (sympy out L26-43, math out L38-67).

No disagreement. `engines_agree: true`.

## Verdict justification

I hand-verified every substitution in the chain against the notes: the four carry-forward identities (§1), the frozen `n=5` EOS / healing-lock relations (`δln c_{s,w} = δlnΘ_w`, `δln ℓ = -δlnΘ_w`), the Stage-250 branch laws for `v_{w0}, T_m`, the parent-coupling exponents (`g_q ∝ Z_q L_W^{-3/2}`, `g_s ∝ T_m a² ℓ`, `λ ∝ v_{w0} a² ℓ L_W^{1/2}`), and the three invariant combinations. Every coefficient (the `13/8`, `15/8`, `-2/5 P`, the `P` cancellation in `v_{w0}`, `δln λ = ½(K_s+K_q)`) reproduces the boxed notes results exactly, and all three invariants plus both channels plus `δ_⊥` reduce to `0`. The `expect_zero`/`expectZero` calls are non-tautological (they traverse the full upstream substitution), the `sqrt(1+rstar²)` in `δ_⊥` poses no false-pass hazard because `chan1 = chan2 = 0` independently, and the symbol domains (`real=True`; `gstar, rstar > 0`) are harmless for linear log-derivative algebra. Both engines agree and the transcripts are fresh (outputs mtime 12:47/13:19 > script mtime 11:58). The paper's headline Output (the three invariants) and the strengthened `δ_⊥ = 0` theorem are fully asserted in both engines, so `paper_alignment: aligned`. The single finding is a low-severity mislabeled banner ("STAGE 150" in a stage-167 audit); the math itself holds up under attack.

## Self-test notes

I checked: (1) coefficient-by-coefficient reproduction of every boxed notes formula — all match, including the `P_0` cancellation in `δln v_{w0}` and the `½(K_s+K_q)` in `δln λ`. (2) Tautology trap — confirmed the invariant/channel checks flow through the full upstream substitution, not define-then-assert. (3) False-pass via `sqrt` in `δ_⊥` — confirmed `chan1` and `chan2` vanish independently, so the radical never matters. (4) SymPy coverage gap — SymPy asserts only D1/D2 (the paper Output), prints D3; Mathematica asserts D3 too, so the paper claim is fully covered across the engine pair. The proposed banner fix is a pure string edit with no math impact and cannot introduce a new paper_misalignment.
