---
unit_id: 162
batch: IV.6
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
  notes_stage_files: [moving_throat_pde_stage162_parent_compensation_rigidity.md]
  paper_appendix: present
---

# Audit unit 162 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_162.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage162_parent_compensation_rigidity.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_162}` at line 1358; no separate narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage162_parent_compensation_rigidity_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage162_parent_compensation_rigidity_mathematica_audit.txt`

## What the paper claims

The stage card's bottom-line claim (the quoted box) is: "Staying on the exact parent compensation family gives automatic D/N similarity preservation and \(\Delta_Q=0\) at first order." The notes sharpen this into concrete identities, all built on the Stage 119 exact parent compensation family (`1+𝔯²=4(𝔤−𝔯)²`, `𝔤±=𝔯±½√(1+𝔯²)`, `L_W/a=(π/2)√((1+𝔯²)/3)`) and the Stages 115–116 odd-normalization `γ₀=(1+𝔯²)/9`. The deliverables are: (1) the automatic similarity identity `δln γ₀ − 2 δln(L_W/a)=0`, i.e. `Ξ_slip=0` along the whole family; (2) the lower-branch differential law `δ𝔤 = (1 − 𝔯/(2√(1+𝔯²))) δ𝔯` with the exact positive decomposition `(4+3𝔯²)/(2√(1+𝔯²)(2√(1+𝔯²)+𝔯))`, hence first-order rigidity `δ𝔤=0 ⟹ δ𝔯=0`; (3) the numerical rigidity factor on Family-1 `𝔯_F1≈1.77799353547498`, `(d𝔤/d𝔯)|_F1≈0.564199521046343`, `δ𝔯≈1.77242…δ𝔤`. The card's `Checks` list (deviations about the renormalized canonical point; even-preservation gate before the odd defect; tangent motion gives `δ⊥=0`) is conceptual framing; the audit-target row points at the four scripted checks.

## What the script claims to verify

Both scripts (identical 4-check structure) verify, symbolically over real `r`: (1) `d ln γ₀ − 2 d ln(L_W/a) = 0` from `γ₀=(1+r²)/9` and `L_W/a=(π/2)√((1+r²)/3)`; (2) `d(g_lower)/dr − (1 − r/(2√(1+r²))) = 0` from `g_lower=r−½√(1+r²)`; (3) `slope − (4+3r²)/(2√(1+r²)(2√(1+r²)+r)) = 0`, establishing manifest positivity; (4) a numerical evaluation at `r_F1=1.77799353547498` giving `(dg/dr)|_F1≈0.5642…` and its reciprocal `≈1.7724…`. The docstring states exactly these four checks. The carry-forward prints restate `Ξ_slip=0` and the `δ𝔤=0 ⟹ δ𝔯=0` rigidity conclusion.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `δln γ₀ − 2 δln(L_W/a)=0` ⟹ `Ξ_slip=0` (notes §2) | sympy L53 / wl L43 `expect_zero("similarity identity", dlog_gamma − 2 dlog_L)` | match |
| Lower-branch differential law `δ𝔤=(1−𝔯/(2√(1+𝔯²)))δ𝔯` (notes §3) | sympy L58 / wl L47 | match |
| Exact positive decomposition `(4+3𝔯²)/(…)` (notes §3 box) | sympy L67 / wl L50 | match |
| Rigidity `δ𝔤=0 ⟹ δ𝔯=0` (notes §3, §5) | not asserted; established by combining checks 2+3 (slope ≠ 0 ⟹ injective) and printed as carry-forward | partial (positivity proven; implication argued in prose) |
| Numerical `𝔯_F1`, `(d𝔤/d𝔯)|_F1`, `d𝔯/d𝔤` (notes §3) | sympy L71–76 / wl L53–58 (printed, no assert) | match (illustrative print) |
| Downstream `Δ_Q=0`, `N_Q−1=0`, `δ𝔅_W=0`, `δγ_W=0` (notes §4) | not independently scripted; carry-forward prints only | partial (consequences of `δ𝔯=0`, stated in prose) |

Dominant pattern: aligned. The three load-bearing symbolic identities each have a faithful non-tautological check in both engines; the downstream collapse (§4) is a prose consequence of `δ𝔯=0`, not a separate scripted target, which is consistent with the card's scope (the card's box claims only "automatic similarity preservation and `Δ_Q=0` at first order," delivered by the `Ξ_slip=0` + rigidity chain).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero(dlog_gamma − 2·dlog_L)` | deliverable 1 (`Ξ_slip=0`) | yes |
| A2 | sympy | 58–61 | `expect_zero(slope − (1 − r/(2√(1+r²))))` | deliverable 2 (lower-branch law) | yes |
| A3 | sympy | 64–67 | `expect_zero(slope − (4+3r²)/(…))` | deliverable 2 (positivity decomp) | yes |
| A4 | sympy | 71–76 | numeric prints of `r_F1`, slope, 1/slope | deliverable 3 (numeric rigidity) | partial (print, no assert) |
| B1 | mathematica | 43 | `expectZero[dlogGamma − 2·dlogL]` | deliverable 1 | yes |
| B2 | mathematica | 47 | `expectZero[slope − (1 − r/(2√(1+r²)))]` | deliverable 2 | yes |
| B3 | mathematica | 49–50 | `expectZero[slope − (4+3r²)/(…)]` | deliverable 2 | yes |
| B4 | mathematica | 53–58 | numeric prints | deliverable 3 | partial (print) |

A1–A3 / B1–B3 are non-tautological: each defines an expression independently (`slope = diff(g_lower, r)`, the candidate RHS `1 − r/(2√(1+r²))`, the decomposition `(4+3r²)/(…)`) and asserts their *difference* vanishes. The RHS is NOT the construction of the LHS — a wrong differentiation or a wrong decomposition would leave a nonzero residual. So the checks genuinely exercise the claimed identities. A4/B4 are illustrative numeric prints with no assertion (consistent with the notes treating them as "numerically substantial, not marginal" color, not a closure gate).

## Findings

### F1 — stale numbering reference (script comment)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py:39`

**What's wrong:**
The SymPy script's source comment reads:
> `# Exact parent family formulas from Stages 99 and 102.`
This cites the wrong source stages. The notes for this stage attribute these exact formulas to **Stage 119** (parent compensation family + balance law + `L_W/a=(π/2)√((1+𝔯²)/3)`) and **Stages 115–116** (`γ₀=(1+𝔯²)/9`):
> notes:40 "Stage 119 rewrote the compensated throat-core condition in terms of the normalized parent ratios"
> notes:58–63 "And Stages 115–116 fix the bare odd normalization by … `γ₀=(1+𝔯²)/9`."
This is the same stale-numbering class that the Step-0 fix already corrected in the notes (`Stages 99 and 170` → `Stage 119`). The garble survived in this script comment as "99 and 102". The Mathematica script carries no such comment, so it is not affected.

**Why this matters:**
It is a non-load-bearing comment (it does not feed any assertion; the formulas are typed directly and the checks still pass), so the math is correct. But it misattributes provenance and propagates the known +17-era numbering drift; a reader auditing provenance would be sent to the wrong stages.

**Required change:**
Edit `scripts/...stage162..._sympy_audit.py:39` from
`# Exact parent family formulas from Stages 99 and 102.`
to
`# Exact parent family formulas from Stage 119 (with gamma0 from Stages 115-116).`

**Verification:**
After the edit, `grep -n "Stages 99 and 102" scripts/...stage162..._sympy_audit.py` returns nothing; the new comment cites Stage 119 (and 115–116 for gamma0); the script still runs to exit 0 with all three `expect_zero` checks at 0 (the comment change cannot affect execution).

## Independent-derivation check (Mathematica)

The two scripts share the same four-check choreography, the same variable names (`r`, `dr`), the same print labels, and the same numeric literal `1.77799353547498`. This parallel structure is, however, forced by the physics: `gamma0`, `Lratio`/`lRatio`, and `g_lower`/`gLower` are **paper-supplied givens** (Stage 119 / 115–116 closed forms), not quantities either engine derives — so identical typing of the starting expressions is expected and legitimate, not an echo. From those givens each engine independently computes the derivatives and simplifications with its own native machinery:

- sympy L49 `sp.diff(sp.log(gamma0), r)` + `sp.simplify`; wl L39 `D[Log[gamma0], r]` + `FullSimplify` — distinct symbolic engines recomputing the same derivative, not the `.wl` re-typing the `.py`'s *result*.
- sympy L28–32 `expect_zero` (`simplify(expand(...))`, raise on `!=0`) vs wl L20–24 `expectZero` (`FullSimplify[Together[Expand[...]]]`, `Exit[1]` on `=!= 0`) — different normalization pipelines (`Together` added on the Mathematica side).
- The outputs differ in surface form in exactly the way two independent CAS engines would: sympy emits `pi*sqrt(r**2/3 + 1/3)/2` (out L6) while Mathematica emits `(Pi*Sqrt[1 + r^2])/(2*Sqrt[3])` (out L6) for the SAME `L_W/a` — a genuinely independent simplification, not a transliteration of one literal form.

Conclusion: this is the borderline "both engines type the same paper givens and differentiate" case, but because (a) the starting expressions are non-derivable givens and (b) each engine independently recomputes the derivatives/simplifications and produces independently-normalized output, it does **not** rise to a `mathematica_transliteration` finding. No finding raised.

## Engine cross-check

The two engines agree on every emitted value:
- `similarity identity = 0` (sympy out L10; wl out L10–11 PASS)
- `lower-branch differential law = 0` (sympy out L12; wl out L13–14 PASS)
- `positive slope decomposition = 0` (sympy out L13; wl out L15–16 PASS)
- `(dg/dr)|_F1`: sympy `0.564199521046342514…` (out L16) vs wl `0.564199521046342510…` (out L19) — agree to ~16 sig figs (the tail differs only at machine-precision noise from the `SetPrecision`/`Float` of an inexact literal).
- `dr/dg|_F1`: sympy `1.77242263188284677…` (out L17) vs wl `1.77242263188284678…` (out L20) — agree to ~16 sig figs.

No engine disagreement.

## Verdict justification

I attacked: (1) tautology — the three `expect_zero`/`expectZero` checks each compare an independently-differentiated LHS against an independently-typed candidate RHS, so a wrong derivative or decomposition would leave a nonzero residual; not tautological. (2) hardcoded result — the only literal is `r_F1=1.77799353547498`, which I confirmed equals the canonical Family-1 radius `√(4107−100π²)/(10π)` to full precision and matches the notes' stated `𝔯_F1`; it drives only an illustrative print, not a closure assertion. (3) symbol-assumption — both engines declare `r, dr` real; positivity of `1+r²` and the decomposition denominator `2√(1+r²)+r` holds for all real `r` (the notes' "> 0 for every real `r`"), so no hidden assumption smuggles in the result. (4) transliteration — judged borderline-but-clean (see above). (5) stale `168π²`/`100π²`/`168%` — none of those constants appear here; the canonical `(1+r²)/9` and `(π/2)√((1+r²)/3)` are used exactly. Paper alignment is exact on the three load-bearing identities; the downstream §4 collapse is a prose consequence of `δ𝔯=0`, consistent with the card's scope. The numbering Step-0 fix is confirmed landed in the notes (Stage 119 at lines 40 and 221) and the card carries no `+17` self-label in `\stagefield{Purpose}` (179 not present). The sole defect is the non-load-bearing stale comment `Stages 99 and 102` in the `.py`, hence verdict `findings` with one low-severity item.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

Outputs are FRESH (sympy .txt 2026-05-28 11:30 > .py 2026-05-27 23:11; wl .txt 2026-05-28 11:32 > .wl 2026-05-27 23:11), so this reconciliation rests on the committed transcripts plus the source.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `γ₀=(1+𝔯²)/9` | py L40 / wl L31; out L5 | notes:61–63 (boxed `γ₀=(1+𝔯²)/9`) | MATCH |
| `L_W/a=(π/2)√((1+𝔯²)/3)` | py L41 / wl L32; out L6 | notes:53–56 (boxed) | MATCH |
| `g_lower=𝔯−½√(1+𝔯²)` | py L42 / wl L33; out L7 | notes:121–124 (`𝔤_−(𝔯)`) | MATCH |
| similarity identity `δln γ₀ − 2 δln(L_W/a) = 0` | py L53 / wl L43; out L10 | notes:92–105 (boxed `=0`, `Ξ_slip=0`) | MATCH |
| lower-branch slope `1 − 𝔯/(2√(1+𝔯²))` & decomposition `(4+3𝔯²)/(2√(1+𝔯²)(2√(1+𝔯²)+𝔯))` | py L60,64–66 / wl L47,49; out L11,15 | notes:127–142 (boxed prefactor) | MATCH |
| `𝔯_F1=1.77799353547498`, `(d𝔤/d𝔯)|_F1≈0.564199521046343`, `d𝔯/d𝔤≈1.77242263188285` | py L71–76 / wl L53–58; out L15–17 | notes:156–164 (all three quoted) | MATCH |

INTERNAL scaffolding (no finding): `banner`/`expect_zero`/`expectZero`/`pass`/`fail`/`fmt` helpers; `dlog_gamma`/`dlog_L` intermediates; pass-flag prints; the four carry-forward prose lines.

All six emitted deliverable values reconcile against the notes. (The `.tex` card is deliberately terse — it states only the box "automatic D/N similarity preservation and `Δ_Q=0`" — so per the augmentation guards the numeric/symbolic deliverables living correctly in the `.md` notes count as MATCH and the card's omission is not a MISSING.) No MISMATCH, no MISSING-DELIVERABLE. The only finding remains F1 (stale comment), which is a `stale-numbering` issue, not a value reconciliation issue.

## Self-test notes

I checked: (1) variable independence — every `diff/D` is w.r.t. `r`, and `gamma0`, `Lratio`, `g_lower` all genuinely depend on `r`, so no derivative is identically zero (the `expect_zero` residuals are nontrivially zero, the printed intermediate derivatives `2 dr r/(r²+1)` etc. are nonzero). (2) trivial-case — at `r=0`: `slope = 1`, RHS `1 − 0 = 1`, decomposition `4/(2·1·2)=1`; all three residuals reduce to 0, and the similarity residual `2r/(r²+1) − 2·r/(r²+1) = 0` identically. (3) No integrals/parity to check. F1 is comment-only, cannot affect execution; the self-test confirms the single prescribed edit cannot introduce a new `paper_misalignment`.
