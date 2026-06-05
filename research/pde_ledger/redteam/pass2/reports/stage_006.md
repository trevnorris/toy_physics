---
unit_id: 006
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
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

# Audit unit 006 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_006.tex`
- notes: `(none)` (no files match `notes/stages/moving_throat_pde_stage006_*.md`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 34 references stage 006; `\input` at line 91)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.txt`

## What the paper claims

Stage 006 exports the **vector form** of the projected Maxwell system. Two field
layers are introduced: the *measured* fields `E_meas = ∫ W E dw`,
`B_meas = ∫ W B dw` (eq:stage006-measured-fields) and the *source-coupled flux*
fields `D_flux = ∫ W Z E dw`, `H_flux = ∫ W Z B dw` (eq:stage006-flux-fields),
with the explicit note that identifying the two pairs is an *additional* closure
assumption, not part of the parent projection. With the stated convention
`F^{i0} = E_i`, the card fixes the homogeneous (Bianchi) signs as
`∂_t B_meas + ∇×E_meas = 0` (and implicitly `div B_meas = 0`), and the
inhomogeneous Ampère sign as
`∇×H_flux − ∂_t D_flux = μ₀ J_proj + L_mix` (eq:stage006-ampere, lines 39–40),
where `L_mix` is the transverse-leakage vector inherited from stage 005.
The `\stagefield{Output}` paragraph states: "Stage 006 exports the field split
(eqs measured–flux) and the vector-law sign convention used by the rest of the
projected Maxwell block." The appendix row (line 34) summarizes the deliverables
as "Measured-field homogeneous equations, source-flux inhomogeneous equations,
and the open-system charge balance" (the Gauss-like law `div D + Leak0 = μ₀ ρ`).

Distinct deliverables: (D1) the measured/flux field-split definitions; (D2) the
homogeneous laws `div B_meas = 0` and `∂_t B_meas + ∇×E_meas = 0`; (D3) the
Gauss-like / charge-balance law `div D_flux + Leak0 = μ₀ ρ_proj`; (D4) the
Ampère-like law with the displayed sign `∇×H_flux − ∂_t D_flux = μ₀ J + L_mix`;
(D5) the F^{i0}=E_i sign convention. No numeric constants are stated as
deliverables anywhere in the card or appendix row.

## What the script claims to verify

Both engines independently rebuild the projected vector system from antisymmetric
tensor data and check four things: (a) the homogeneous Bianchi rearrangement
yields `div B = 0` and `∂_t B + ∇×E = 0` (sympy lines 78–107; wl M1, lines
72–98); (b) the inhomogeneous tensor divergence `∂_μ G^{μν} + Leak` rearranges to
a Gauss-like law and an Ampère-like law (sympy 124–159; wl M2, 100–124); (c) a
concrete Gaussian-projected bulk potential satisfies all four projected laws with
the source `J` *defined from* the bulk tensor divergence (sympy 198–305; wl M4,
135–173); (d) a leakage-normalization value and several adversarial
sign-mutation / parity guards (sympy 252–369; wl M3/M5). The docstring
(sympy 1–12) states the goal as separating measured fields from source-coupled
flux fields and confirming the projected theory distinguishes (E,B) from (D,H).
The script's load-bearing Ampère assertion (sympy `amp1_target`, line 146;
wl `ampereRearranged`, line 119–122) is `−∇×H_flux − ∂_t D_flux + Leak = μ₀ J`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 field split `E_meas=∫WE`, `D_flux=∫WZE`, etc. | sympy 42–55 defs + concrete 218–219, 236–237/267–268; wl 73–74, 101–102, 152–155 | match |
| D2a `div B_meas = 0` | sympy `divB` 78 + concrete 238; wl M1 `divCycleResidual` 93–98, M4 163 | match |
| D2b `∂_t B_meas + ∇×E_meas = 0` | sympy `checks_faraday` 94–107 + concrete 239–250; wl M1 82–97, M4 164 | match (output 13–15 = standard sign) |
| D3 `div D_flux + Leak0 = μ₀ ρ` | sympy `lhs0/rhs0` 131–135 + concrete 290–293; wl M2 `gaussRearranged` 118/123, M4 167 | match |
| D4 Ampère `∇×H_flux − ∂_t D_flux = μ₀ J + L_mix` | sympy `amp*_target` 146–159 (= `−∇×H − ∂_t D + L = μ₀ J`); wl `ampereRearranged` 119–124 (same) | **mismatch** (curl-H sign and leak sign differ from card) |
| D5 sign convention `F^{i0}=E_i` | sympy uses lower-index `F_{i0}=E_i` (comment 72, defs 73–75); wl `fieldF[i,0]=E_i` (78) | partial (index placement up vs down; homogeneous output signs nevertheless match the card, so not load-bearing) |

`paper_alignment: partial` — D1, D2, D3 match exactly; D4 (Ampère sign) is a
genuine `target_mismatch`; D5 is a minor notational up/down-index difference whose
observable consequences (the displayed Faraday/divB signs) still match.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 106–107 | `simplify(far_i − (∂_t B_i + curlE_i)) == 0` | D2b Faraday | yes |
| A2 | sympy | 158–159 | `simplify(lhs_i − amp_i_target) == 0`, with `amp_i_target = −∇×H − ∂_t D + L` | D4 Ampère (opposite curl/leak sign) | partial (verifies a sign-flipped identity vs card) |
| A3 | sympy | 238 | concrete `div B_proj == 0` | D2a divB | yes |
| A4 | sympy | 239–250 | concrete projected Faraday == 0 | D2b | yes |
| A5 | sympy | 290–293 | concrete Gauss `div D + leak0 − μ₀ Pg(J0) == 0` (J0 defined from tensor div) | D3 | yes (sign fixed by construction of J0) |
| A6 | sympy | 294–305 | concrete Ampère with J_i defined from tensor div | D4 | partial (J_i defined from same tensor → cannot independently pin curl-H sign vs card) |
| A7 | sympy | 220 | `E_meas_example − D_flux_example != 0` | D1 (E≠D for nontrivial Z) | yes |
| A8 | sympy | 257–259 | boundary==0, IBP identity, leak==transverse derivative | D3/D4 leakage machinery | yes |
| A9 | sympy | 260 | `leak1 − sqrt(2)/4 == 0` (numeric normalization) | none (internal, test-profile value) | n/a |
| A10 | sympy | 261, 366–369 | sign-mutation guards `!= 0` | D2b/IBP robustness | yes (guards fire) |
| A11 | mathematica | 97–98 | M1 Faraday & divB rearrangement == 0 | D2 | yes |
| A12 | mathematica | 123–124 | M2 Gauss & Ampère rearrangement == 0 (`ampereCurl3 = −curl3`) | D3/D4 | partial (D4 same sign-flip as A2) |
| A13 | mathematica | 163–173 | M4 concrete Bianchi/Gauss/Ampère == 0 | D2/D3/D4 | yes (D4 sign by construction) |
| A14 | mathematica | 131–133 | M3 boundary==0, IBP, `vectorLeakMoment − 1/(2√2) == 0` | none (internal) / IBP machinery | n/a / yes |
| A15 | mathematica | 178–179 | M5 sign-mutation guards `!= 0` | robustness | yes (residuals 1/√2, −3) |

## Findings

### F1 — paper_misalignment

**Subtype:** target_mismatch
**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_006.tex:38-41` (eq:stage006-ampere)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py:146-159`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage006_projected_maxwell_vector_mathematica_audit.wl:119-124,171`

**What's wrong:**
The paper card displays the projected Ampère law with the standard textbook sign:
> `\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux} =\mu_0\mathbf J_{\rm proj}+\mathbf L_{\rm mix}.` (stage_006.tex:38–41)

i.e. rearranged, `∇×H − ∂_t D − L_mix = μ₀ J`.

The script's load-bearing symbolic assertion verifies the **opposite curl sign and
opposite leak sign**. In SymPy the target the tensor divergence is checked equal
to is
```
amp1_target = (∂_z H2 − ∂_y H3) − ∂_t D1 + L1   # = −(∇×H)_1 − ∂_t D1 + L1
```
(`amp*_target`, lines 146–148; `check_amp` asserts the tensor divergence equals
this, lines 150–159; the printed law sets it `= μ₀ J`, lines 153–155, output
lines 37–39). So the script's projected Ampère law is
`−∇×H_flux − ∂_t D_flux + L = μ₀ J`. The Mathematica side derives the *same*
sign-flipped form independently: `ampereCurl3[v][i] = Σ eps3[k,j,i] ∂_j v[k] =
−(curl3 v)[i]` (definition lines 65–68), and `ampereRearranged` asserts the
tensor divergence equals `ampereCurl3[Hflux] − ∂_t Dflux + leak = −∇×H − ∂_t D + L`
(lines 119–124, M2 residual 0 in output line 8). Both engines agree with each
other and disagree with the card. The two forms are not algebraically equal for
general fields (equality would require `∇×H ≡ L` identically), so this is a true
identity-level discrepancy, not a notation rewrite.

Note the concrete-profile Ampère checks (sympy A6 lines 294–305; wl A13 lines
169–173) cannot adjudicate this, because there the current `J` is *defined from*
the same bulk tensor divergence (sympy `J1_bulk` etc. lines 278–289; wl
`bulkCurrent` lines 158–161), so those checks pass for either overall sign by
construction. The symbolic M2/A2 rearrangement is the only check that pins the
curl-H sign, and it pins it opposite to the card.

(The related D5 item — card writes `F^{i0}=E_i` (upper indices) while the script
uses `F_{i0}=E_i` (lower indices, sympy comment line 72, wl line 78) — is a pure
index-placement difference whose observable Faraday/divB signs still match the
card, so it is not raised as a separate finding; it is mentioned only because it
is the convention statement adjacent to the mismatched Ampère sign and may be the
root of the discrepancy.)

**Why this matters:**
The card's `\stagefield{Output}` explicitly states the stage exports "the
vector-law sign convention used by the rest of the projected Maxwell block." If
the exported Ampère sign printed in the paper is opposite to the one both audit
engines verify, downstream stages that cite stage 006's Ampère sign convention
could inherit an inconsistent sign. The displayed equation is a stated deliverable
(D4), so the audit cannot certify the card's Ampère line as verified.

**Required change:**
Route to user — see `## Resolve before fix_loop` in the directive. Codex must not
auto-edit either side. (Either the card's `∇×H` should be `−∇×H` and `+L_mix`
moved/sign-adjusted to match the script, or the script's `amp*_target` /
`ampereCurl3` sign should be flipped to match the textbook card — but which is
"correct" depends on the intended `F^{μν}` raising convention from stages 004/005,
which is the user's call.)

**Verification:**
After the user picks a direction: if paper-side, the card line 39–40 reads the
form matching `amp*_target`; if script-side, `amp*_target` (lines 146–148) and
`ampereCurl3` (lines 65–68) are flipped, M2/M4 still exit 0, and the printed law
(output lines 37–39) matches the card. The two displayed Ampère forms must become
identical.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. The SymPy script hand-writes
every tensor component and every curl component explicitly (e.g. `far1 =
diff(F23,t)+diff(F30,y)+diff(F02,z)`, lines 86–88; `curlE1 = diff(E3,y)−diff(E2,z)`,
line 90). The Mathematica script instead builds the field tensors abstractly via
`LeviCivitaTensor[3]` contractions (`fieldF[i,j] := Σ eps3[k,i,j] Bvec[k]`,
lines 79–80; `fluxG`, lines 107–111) and derives the rearrangements by summing
`projectedInhom[nu] := Σ_μ D[fluxG[μ,ν], coord]` (lines 113–116) and the
epsilon-defined `curl3`/`ampereCurl3` (lines 61–68). The choreography, naming, and
intermediate objects differ substantially; these are two genuinely independent
constructions that happen to converge on the same (sign-flipped-vs-card) Ampère
form — which actually strengthens confidence that the discrepancy is paper-vs-script,
not an artifact of one engine. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at every checked level (SymPy raises no AssertionError, all
residues 0/expected; Mathematica output residuals all 0 except the deliberate M5
mutations 1/√2 and −3). The leakage normalization agrees numerically: SymPy
`leak1 = √2/4` (line 260, output implicit via PASS) and Mathematica
`vectorLeakMoment = 1/(2√2) = √2/4` (line 133, output line 14). No
`engine_disagreement`. Crucially, both engines independently produce the same
`−∇×H − ∂_t D + L` Ampère form, so the F1 discrepancy is robust across engines.

## Verdict justification

`verdict: findings`. The homogeneous block (div B = 0, Faraday) and the Gauss-like
charge-balance law match the card exactly and survive the sign-mutation guards;
the field-split definitions match; both engines are independent and agree. The one
real defect is F1: both engines verify a projected Ampère law whose `∇×H` and leak
signs are opposite to the card's displayed `eq:stage006-ampere`, and the concrete
checks cannot adjudicate it because `J` is defined from the same tensor divergence.
This is a stated deliverable (D4), so it is a `paper_misalignment / target_mismatch`
routed to the user; it is not auto-`CRITICAL_DOWNSTREAM` (downstream impact is
assessed after the user fixes the direction). Attacks tried that failed:
tautology check on Faraday/Gauss (non-tautological — they assert a tensor
rearrangement equals an independently written vector form); parity check on the
leak integrals (correct: `Wgp·Z·w` even → nonzero, `Z=w` makes it odd → zero, as
asserted); variable-independence on the `assert_nonzero` mutation guards (they
fire, residuals 1/√2 and −3); symbol-domain check (t,x,y,z,w real, mu0 nonzero —
appropriate, no positivity over-assumption). I confirm I read the paper card, the
(absent) notes, and the appendix row before auditing the scripts, and that every
deliverable except D4 matches the script.

## Value Reconciliation (pass-2 augmentation)

Outputs are FRESH (both `.txt` mtimes 2026-05-25 17:24/17:29 are newer than both
script mtimes 2026-05-25 02:13), so the reconciliation rests on script source +
committed outputs without running anything.

The stage's stated deliverables (card + appendix row) are entirely **symbolic**
(field-split definitions, homogeneous laws, Gauss-like law, Ampère law, sign
convention). The scripts emit exactly **one** non-scaffolding numeric value, the
leakage normalization for the specific hand-chosen test profile (Z=Gaussian,
F_{w1}=w):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| leakage normalization `√2/4 = 1/(2√2) ≈ 0.35355` | sympy `leak1 − sp.sqrt(2)/4` line 260 (PASS, output implicit); wl `vectorLeakMoment − 1/(2 Sqrt[2])` line 133, output line 14 (`residual = 0`) | absent from both `.tex` (card states no constants) and `.md` (notes absent) | INTERNAL — not a stated deliverable; it is the value of a leakage integral for one illustrative test profile, used only to anchor the IBP/leak machinery. Per augmentation guard, NOT flagged MISSING. |
| projected Ampère law `−∇×H − ∂_t D + L = μ₀ J` (symbolic deliverable) | sympy `amp*_target` 146–148; wl `ampereRearranged` 119–124 | card eq:stage006-ampere lines 39–40 displays `∇×H − ∂_t D = μ₀ J + L_mix` | MISMATCH → folded into F1 (`value_mismatch`/`target_mismatch`, see Findings). |

Symbolic deliverables that reconcile (MATCH): field split E_meas/B_meas/D_flux/H_flux
(card eqs 15–27 ↔ sympy 42–55/218–219/236–237/267–268, wl 73–74/152–155);
`div B_meas = 0` (card ↔ sympy 78/238, wl 93–98/163); Faraday `∂_t B + ∇×E = 0`
(card eq 35 ↔ sympy 94–107/239–250 output 13–15, wl 82–97/164); Gauss-like
`div D + Leak0 = μ₀ ρ` (appendix "open-system charge balance" ↔ sympy 131–135/290–293,
wl 118/167).

INTERNAL scaffolding (accounted-for, no finding): pass/fail residues `[0,0,0]`/`0`;
M5 sign-mutation residuals `1/√2` and `−3` (deliberately nonzero guards);
`E_meas_example − D_flux_example` (E≠D demonstration); bogus-projection / asymmetric-weight
/ antisymmetric-Z demonstration residuals; all boundary-term `=0` checks; intermediate
`Pg`/`Pgp`/`boundary` quantities.

`reconciliation: complete; 1 numeric value + 5 symbolic deliverables checked,
1 misaligned` (the symbolic Ampère deliverable D4, folded into F1; the single
numeric value √2/4 is internal and reconciles as non-deliverable).

## Self-test notes

Checked the traps: (1) variable independence — the `assert_nonzero` mutation
guards (sympy 261, 366–369; wl 178–179) genuinely depend on the differentiated
symbols and the outputs show nonzero residuals (1/√2, −3), so they fire rather
than passing trivially. (2) Parity/symmetry on the leak integrals — `Wgp·Z·w`
with Z=Gaussian is even (nonzero integral, value √2/4), and Z=w makes the
integrand odd (zero), exactly as the asserts claim; the symmetric vs asymmetric
Gaussian weight checks are parity-correct. (3) Trivial-case — the concrete Gauss/
Ampère checks define `J` from the tensor divergence so they pass for either
overall sign, which is precisely why they cannot adjudicate F1; the discriminating
check is the symbolic M2/A2 rearrangement, which I verified pins the curl-H sign
opposite to the card. No directive self-test for a Codex edit is needed because
the sole finding is a `paper_misalignment` routed to the user, not a Codex-applied
script edit.
