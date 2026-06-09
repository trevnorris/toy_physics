---
unit_id: 008
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 008 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_008.tex`
- notes: `(none)` (no file matching `notes/stages/moving_throat_pde_stage008_*.md`)
- part appendix: `/var/projects/toy_physics/.../stage_appendix_part01.tex` → `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 008 at line 38; `\input` at line 95)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.txt`

## What the paper claims

Stage 008 (Part I, anchor MTDC-T4) extends the projected-Maxwell construction to a
weighted gauge-fixing profile H(w) and identifies the **clean matching channel**
under which projection-first reproduces the reduction-first brane coupling. The
card states two matching conditions verbatim:

- Clean gauge channel (`eq:stage005-hz`): `H(w)=Z(w)`, `\xi_{\rm eff}^{\rm proj}=\xi`.
- Source match (`eq:stage005-source-match`): `S(w)=Z(w)/Z_{\rm int}`, `Z_{\rm int}=\int Z(w)\,dw`, giving `\mu_{\rm eff}^{\rm proj}=\mu_0/Z_{\rm int}` ("the reduction-first coupling").

`\stagefield{Output}`: *"Stage~008 exports the matching conditions
\eqref{eq:stage005-hz}--\eqref{eq:stage005-source-match}. These conditions are the
firewall between exact projection and reduced brane coupling."* The card body also
states the audit "checks the Gaussian case and the regulator limit, including
mutation guards distinguishing H=1 from H=Z." The appendix row (line 38) summarizes
the stage as the "Slot-level distinction between source-flux projection and matched
reduction," status `\StatusExact{} / \StatusReduced{}`. Distinct deliverables:
(D1) H=Z ⇒ ξ_eff^proj = ξ for any normalized W; (D2) S=Z/Z_int ⇒ μ0_eff^proj =
μ0/Z_int for any normalized W; (D3) H=1 produces a W-dependent ξ_eff (the
mutation/distinction that motivates the matching); (D4) the reduction-first H=1
regulator limit is degenerate (ξ4(R)→0), which is why H=Z is the clean choice.

## What the script claims to verify

The SymPy script (docstring lines 2-18) states it (1) writes the exact projected
brane law for a general weight H(w), (2) derives the zero-mode effective coupling
`mu0_eff = mu0·I_WS/I_WZ` and gauge parameter `xi_eff = xi·I_WZ/I_WH`, and
(3) compares H=1 vs H=Z with Gaussian matched-kernel examples, asserting that for
H=Z, ξ_eff = ξ for any normalized W. Sections 5 and 5b carry the substantive
verification via real `sp.integrate` over the line: a matched kernel W=Z/Z_int and
an *independent* Gaussian kernel of width σ≠λ. The Mathematica script verifies the
same structural identities (M1 reciprocal, M2 Gaussian and Lorentzian H=Z, M3/M4
matched-source coupling and concrete Gaussian values, M5 regulator limit, M6 H=Z
reduction identity, M7 Lorentzian numeric independent-profile checks). Both scripts
print STATUS: PASS.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: H=Z ⇒ ξ_eff^proj=ξ, any normalized W | py L205 `xi_eff_HZ_match - xi`; py L247 `xi_eff_HZ_indep - xi` (σ≠λ); wl M2a L48, M2b L73, M4 L134, M7 L183 | match |
| D2: S=Z/Z_int ⇒ μ0_eff^proj=μ0/Z_int, any normalized W | py L207 `mu0_eff_source_match - mu0/Z_int`; py L265 independent W; wl M3 L87, M4 L143, M7 L186 | match |
| D3: H=1 ⇒ W-dependent ξ_eff (distinction/mutation guard) | py L206 `xi_eff_H1_match - xi/√2`; py L256-259 σ≠λ guard asserting ξ_eff(H=1)≠ξ; wl M4 L137 `xi/√2` | match |
| D4: reduction-first H=1 regulator limit ξ4(R)→0 | py L288 `xi4_unweighted_limit==0`; wl M5 L153 | match |
| Z_int=∫Z dw (Gaussian = √π λ) | py L201; wl M4 L125 | match |
| μ_eff^proj=μ0/Z_int boxed result | py L156/L207; wl M3/M4 | match |

`paper_alignment: aligned`. Every card-stated matching condition has a
non-tautological script-side check (the integrals are actually computed, in two
engines, including kernels not proportional to Z).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 74 | `inv_xi·xi - 1 == 0` | none (reciprocal bookkeeping) | no (self-labeled tautology) |
| A2 | sympy | 130-133 | `I_WH[H→Z] - I_WZ == 0` | D1 (anchor only) | no (self-labeled tautology) |
| A3 | sympy | 137-140 | gauge-driver factoring | D1 (anchor only) | no (self-labeled tautology) |
| A4 | sympy | 141-144 | `xi·I_WZ/I_WH[H→Z] - xi == 0` | D1 (anchor only) | no (self-labeled tautology) |
| A5 | sympy | 156 | `mu0·(I_WZ/Z_int)/I_WZ - mu0/Z_int` | D2 (anchor only) | no (self-labeled cancellation tautology) |
| A6 | sympy | 201 | `Z_int_gauss - √π·λ == 0` | Z_int value | yes (integrate) |
| A7 | sympy | 202 | `Z2_int - √(2π)λ/2 == 0` | I_WZ building block | yes (integrate) |
| A8 | sympy | 203 | `I_WZ_match - √2/2 == 0` | D1/D2 (matched W) | yes (integrate) |
| A9 | sympy | 204 | `I_WH_HZ_match - I_WZ_match == 0` | D1 | yes (integrate, distinct integrand) |
| A10 | sympy | 205 | `xi_eff_HZ_match - xi == 0` | D1 | yes |
| A11 | sympy | 206 | `xi_eff_H1_match - xi/√2 == 0` | D3 | yes |
| A12 | sympy | 207 | `mu0_eff_source_match - mu0/Z_int == 0` | D2 | yes |
| A13 | sympy | 208 | `mu0_eff_delta_match/(mu0/Z_int) - √2 == 0` | D4-adjacent (delta mismatch) | yes |
| A14 | sympy | 239 | `W_indep_norm - 1 == 0` | normalization of independent W | yes |
| A15 | sympy | 245 | `I_WH_HZ_indep - I_WZ_indep == 0` | D1 (σ≠λ) | yes |
| A16 | sympy | 247 | `xi_eff_HZ_indep - xi == 0` | D1 (σ≠λ, strongest) | yes |
| A17 | sympy | 256-259 | guard: ξ_eff(H=1)≠ξ at σ=1/2,λ=1 | D3 (mutation guard) | yes |
| A18 | sympy | 265-268 | `mu0_eff_source_indep - mu0/Z_int == 0` | D2 (σ≠λ) | yes |
| A19 | sympy | 288 | `xi4_unweighted_limit == 0` | D4 | yes (limit) |
| M1 | mathematica | 19 | reciprocal `xiEff·invXiEff-1` | none | no (symbolic identity) |
| M2a | mathematica | 38,45,48 | Gaussian norm + H=Z gauge | D1 | yes (integrate) |
| M2b | mathematica | 63,70,73 | Lorentzian norm + H=Z gauge | D1 (independent profile) | yes (integrate) |
| M3 | mathematica | 87 | matched source μ0/Z_int (σ-Gaussian W) | D2 | yes |
| M4 | mathematica | 111,125-145 | 7 concrete matched-kernel residuals | D1/D2/D3/D4-adj | yes |
| M5 | mathematica | 153 | regulator limit →0 | D4 | yes (Limit) |
| M6 | mathematica | 162 | `xi·zInt/zInt - xi` | D1 (reduction-first) | no (symbolic identity) |
| M7 | mathematica | 183,186 | Lorentzian numeric H=Z + matched source | D1/D2 (independent profile) | yes (numeric at 2 points) |

## Findings

None. The five tautological `assert_zero` calls on the SymPy side (A1-A5) are
**explicitly self-labeled** as tautologies (lines 73, 120-123, 154-155) and are
retained as readability anchors; the substantive verification is carried by the
real integrations in sections 5 and 5b (A6-A19) and independently by the
Mathematica script. Per the audit rubric a self-disclosed anchor whose substantive
counterpart exists is not a defect — flagging A1-A5 as `tautological_check` would be
a false positive because the load-bearing claim is exercised elsewhere
(non-tautologically) for the same identity. Likewise M1 and M6 are pure symbolic
identities but every physical claim they bookkeep is also checked against actual
computed integrals (M2-M4, M7).

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration of the `.py`. Evidence:
- The Mathematica script introduces a profile the SymPy script never uses: a
  Lorentzian weight `lorentzWeight = (sigma/Pi)/(w^2 + sigma^2)` (wl L26), used for
  M2 Pair B (L53-76) and the M7 numeric checks (L167-189). The SymPy script's
  independent kernel (py L236) is a *Gaussian* of width σ, not a Lorentzian — so the
  two engines test H=Z and the matched-source coupling against **different**
  non-matched profiles. That is genuine independent corroboration, not echoing.
- M7 uses numeric substitution at concrete rationals (`lambda->1, sigma->1/2|2,
  xi->7/5, mu0->11/10`, wl L176-180, evaluated to 30 digits) — a numeric route the
  SymPy script does not take.
- Variable/function choreography differs: the SymPy script builds symbolic Integral
  objects with `subs(H→Z)` substitution gymnastics (py L124-144), which the
  Mathematica script does not mirror; the `.wl` computes closed-form `Integrate[...]`
  directly and compares residuals.

Conclusion: independent derivation confirmed.

## Engine cross-check

Both engines agree at the level they claim. Concrete Gaussian matched-kernel values
agree: SymPy emits `Z_int = sqrt(pi)*lambda`, `Z2_int = sqrt(2)*sqrt(pi)*lambda/2`,
`I_WZ = sqrt(2)/2`, `xi_eff(H=1) = sqrt(2)*xi/2`, delta-source ratio `sqrt(2)`,
matched-source ratio `1` (sympy output L80-96); Mathematica's M4 residual vector is
`{0,0,0,0,0,0,0}` against the identical targets `Sqrt[Pi]*lambda`,
`Sqrt[Pi/2]*lambda`, `Sqrt[2]/2`, `xi`, `xi/Sqrt[2]`, `Sqrt[2]`, `mu0/zArea`
(wl L116-122, output L15). The independent-profile result `I_WZ = lambda/sqrt(lambda^2+sigma^2)`
(sympy output L103) is consistent with the Mathematica Gaussian-Gaussian and
Lorentzian overlaps (M2a/M2b/M7 residuals all 0). No sign, factor-of-2, or branch
disagreement found. `engines_agree: true`.

## Verdict justification

`clean`. Attacks attempted and failed: (1) I checked whether the headline H=Z⇒ξ_eff=ξ
claim was only verified via the trivial substitution `I_WH[H→Z]=I_WZ` — it is also
verified by actually integrating two *distinct* non-matched kernels in each engine
(Gaussian σ≠λ in SymPy A15-A16; Lorentzian in Mathematica M2b/M7), so the "for any
normalized W" universality is genuinely exercised, not assumed. (2) I checked the
delta-source ratio √2 and the H=1 matched ξ_eff=xi/√2 by hand from the closed-form
integrals; both are correct. (3) I checked symbol domains: all profiles use
`positive`/`Reals` assumptions consistent with the physical setup (λ,σ,ξ,μ0,R>0),
and the integrals converge under these. (4) I checked the σ≠λ H=1 guard (A17): it
correctly raises if ξ_eff(H=1) accidentally equals ξ, so the "mutation guard
distinguishing H=1 from H=Z" the card promises is real. (5) Output mtimes
(17:13) are newer than script mtimes (17:11) — outputs fresh; content matches the
current scripts. (6) The Value Reconciliation below confirms every emitted
deliverable value is reflected in the card. I read the card, the (absent) notes, and
the appendix row; the script's verified claims match the card's stated matching
conditions exactly.

Minor cosmetic (non-finding): the equation labels in stage_008.tex are named
`eq:stage005-hz` / `eq:stage005-source-match` (lines 18, 27) — a stale `005` prefix
on labels that are defined and used within stage 008. They resolve correctly
(`pde_ledger.aux` shows them as B.75/B.76 inside the stage-008 block B.8), so this is
purely a naming inconsistency, not a broken reference or a math defect. Noted for the
numbering-reconciliation track, not raised as a red-team finding.

## Self-test notes

I checked the variable-independence trap (no `sp.diff`/`D[]` derivatives appear; the
script is integral-based, so the "derivative w.r.t. an absent symbol" failure mode
cannot occur here). I checked symmetry/parity: all integrands (Gaussian·Gaussian,
Lorentzian·Gaussian) are even in w over the symmetric line, so the integrals are
legitimately nonzero — consistent with the asserted closed forms. I checked
trivial-case substitution: setting σ=λ in the independent-W result λ/√(λ²+σ²) gives
1/√2 = the matched I_WZ, internally consistent; and the delta-source ratio reduces to
√2 by hand. No traps triggered; verdict stands as clean.

## Value Reconciliation (pass-2 augmentation)

Per the augmentation doc, I enumerated every RESULT/deliverable value the two
scripts emit (from source + saved outputs) and located each in the docs. The card is
deliberately terse (it states the matching *conditions* as its Output, not the
intermediate Gaussian numerics), and there is **no `.md` notes file** for this stage,
so the natural carrier for everything is the `.tex` card. Per the augmentation's
guard, intermediate Gaussian/Lorentzian benchmark numbers (Z_int, Z2_int, √2/2,
xi/√2, the √2 delta-source ratio, the σ≠λ closed forms) are **INTERNAL** scaffolding
that demonstrate the deliverable conditions on concrete profiles; they are not
themselves card-level deliverables, so their absence from the terse card is not a
MISSING-DELIVERABLE (it is the legitimate terse-card omission the guard describes).

Deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `H(w)=Z(w)` (clean gauge condition) | py docstring L17, §3/§5; wl M2/M4 | stage_008.tex:19 `H(w)=Z(w)` | MATCH |
| `xi_eff_proj = xi` (H=Z, any normalized W) | py L205/L247; wl M2a/M2b/M4/M7; sympy out L53,L87,L104; wl out L7,L11,L16 | stage_008.tex:20 `\xi_{\rm eff}^{\rm proj}=\xi` | MATCH |
| `S(w)=Z(w)/Z_int` (source match) | py L158/§4; wl M3 | stage_008.tex:28 `S(w)=\frac{Z(w)}{Z_{\rm int}}` | MATCH |
| `Z_int = ∫ Z(w) dw` (definition) | py L262/L181; wl L78 zArea | stage_008.tex:30 `Z_{\rm int}=\int Z(w)\,dw` | MATCH |
| `mu_eff_proj = mu0/Z_int` (reduction-first coupling) | py L156/L207/L265; wl M3 L85, M4 L143; sympy out L63,L95,L106; wl out L13,L16 | stage_008.tex:34 `\mu_{\rm eff}^{\rm proj}=\frac{\mu_0}{Z_{\rm int}}` | MATCH |
| `xi_eff_proj = xi (∫WZ)/(∫WH)` (general form) | py L71/§2/§7; sympy out L36 | implicit via `eq:stage007-effective-parameters` (stage_007.tex:24, cited by stage_008.tex:15); card states the H=Z specialization | MATCH (general form carried by parent stage 007; stage 008 states the specialization the card is responsible for) |

INTERNAL items (genuine scaffolding, no finding): `Z_int_gauss = √π·λ`,
`Z2_int = √(2π)λ/2`, `I_WZ_match = √2/2`, `xi_eff(H=1)=xi/√2` (and the σ≠λ closed
forms `lambda/√(λ²+σ²)`), delta-source ratio `√2`, matched-source ratio `1`,
`xi4_unweighted(R)=√π·λ·xi/(2R)` and its R→∞ limit `0`, the reciprocal residual `0`,
all M-residuals (`0`), M7 numeric residuals (`{0,0}`), and the Lorentzian benchmark
profile. These are concrete-profile demonstrations / verification scaffolding, not
card-level deliverables.

reconciliation: complete; 6 deliverable values checked, 0 misaligned.
