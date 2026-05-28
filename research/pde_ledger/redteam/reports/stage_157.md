---
unit_id: 157
batch: IV.6
auditor_model: claude-opus-4-7
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md
  paper_appendix: present
---

# Audit unit 157 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_157.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage157_core_mouth_coevolution_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_157}` row at line 1348; no inline appendix prose for this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.txt`

## What the paper claims

The card (Stage 157, v59 Status, `\StatusNumerical / \StatusOpen`) is a co-evolving core-mouth transport ledger entry. The block-narrative summary (quoted from the card) is: "Identifies the reduced co-evolving target and leaves the deviation-to-normalization map as the next task." The `\stagefield{Checks}` block enumerates three deliverables: (1) deviations are taken about the renormalized co-evolving canonical point; (2) even-preservation constraints are imposed before reading the remaining odd defect; (3) tangent motion on the parent compensation family gives `delta_perp = 0`. The supporting notes pin the carry-forward numerics: `r_{F1} ≈ 1.77799353547498`, `Sigma_0^can ≈ 4.651033550168867`, `T_hat_can ≈ 1.4467083664567613`, `Pi_can ≈ 3.871564377479002`, and identify `R = 1/4` as the lower-compensated branch label and `r_{F1} - sqrt(1+r_{F1}^2)/2` as the corresponding `g_*` value. No additional row text in the part-IV appendix beyond the `\input` line.

## What the script claims to verify

The SymPy docstring enumerates seven checks: (1) `g = g_* <=> R = 1/4` exactly; (2) the self-matched traction law `Sigma_0 = 20 T_hat^2 / 9` (applied to the renormalized canonical tuple); (3) the carried Stage 156 renormalized canonical tuple satisfies the branch identities (`R_can = 1/4`, `Pi_can = Sigma_0^can (1 - S_can/4)`); (4) the renormalized point sits strictly above the original canonical point in all three coordinates; (5) tangent motion on `g = g_-(r)` keeps `delta R = 0`; (6) tangent motion plus `dE2 = dE4 = 0` forces `deltaC = delta_kappa_W = 0`; (7) the Stage 158 expansion point is the renormalized canonical tuple, with the linearised `Pi`-packet `dPi = (1 - S/4) dSigma_0 - Sigma_0 dS/4`. The Mathematica script mirrors this set in identical order with the same expressions and constants.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| Card check (1): deviations about renormalized co-evolving canonical point | sympy lines 76-93 (carry-forward numerics, `R_can = 1/4`, `Pi_can` identity, self-matched traction law, strict-above ordering); `.wl` lines 63-81 | match |
| Card check (2): even-preservation enforces zero defect before odd reading | sympy lines 103-110 (Solve dE2 = dE4 = 0 → `{deltaC: 0, dkappa: 0}`); `.wl` lines 92-100 | match |
| Card check (3): tangent motion on parent compensation family gives `delta_perp = 0` | sympy line 101 (`tangent motion keeps delta R = 0`); `.wl` line 90 | match (with delta_perp identified as delta R on the family) |
| Card "Identifies the reduced co-evolving target" (block-narrative line) | sympy lines 56-62 (`R(g_*) - 1/4 = 0`) plus the carry-forward tuple identities | match |
| Card "leaves the deviation-to-normalization map as the next task" | sympy lines 115-122 / `.wl` lines 105-113 (Stage 158 expansion handoff) | extra-but-aligned (handoff scaffolding rather than verification of a card claim) |
| Notes: `r_{F1} = sqrt(4107 - 100 pi^2)/(10 pi)` closed form | sympy line 64 (`rF1_exact`); `.wl` line 51 | match |
| Notes: numeric values `Sigma_0^can ≈ 4.6510…`, `T_hat_can ≈ 1.4467…`, `Pi_can ≈ 3.8715…` | constants loaded from `scripts/numerical/stage155_156_fixedpoint_samples.json` (matches notes to 12+ digits) | match |

Overall alignment is good. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 61 | `expect_zero("R(g_*) - 1/4", R.subs(g, g_star) - 1/4)` | block-narrative: "reduced co-evolving target" (R = 1/4 branch label) | yes |
| A2 | sympy | 76 | `expect_close("r_F1 radical check", rF1_exact, rF1)` | notes carry-forward (closed-form r_F1 vs JSON value) | yes |
| A3 | sympy | 77 | `expect_close("g_* lower-branch check", ...)` | notes carry-forward (closed form for g_*) | yes |
| A4 | sympy | 78 | `expect_close("self-matched traction law", Sigma0_can, 20/9 T_hat_can^2)` | card check (1): renormalized canonical point identity | yes |
| A5 | sympy | 81 | `expect_close("R_can = 1/4", R_can, 1/4)` | card check (1) | yes |
| A6 | sympy | 82 | `expect_close("Pi_can identity", Pi_can, Sigma0_can*(1 - S_can/4))` | card check (1) | yes |
| A7 | sympy | 92 | `assert Sigma0_can > Sigma0_star and T_hat_can > T_hat_star and Pi_can > Pi_star` | notes "renormalized point sits above original canonical point" | yes |
| A8 | sympy | 101 | `expect_zero("tangent motion keeps delta R = 0", dR.subs(dg, gp*dr))` | card check (3) | yes |
| A9 | sympy | 109-110 | `assert even_preservation == [{deltaC:0, dkappa:0}]` | card check (2) | yes |
| A10 | sympy | 113 | `expect_zero("tangent motion kills delta C", -16 sigma_star * (dR.subs(dg, gp*dr)))` | (claims card check (2) but algebraically reduces to A8) | partial (trivial multiple of A8) |
| A11 | sympy | 121 | `expect_zero("Stage 158 tangent expansion packet", ...)` | Stage 158 handoff (not a card check) | yes (anchored to dPi formula) |
| B1-B11 | mathematica | 48-112 | same checks in same order with identical algebraic content | same as A1-A11 | same |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:1-124`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:1-128`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script: same helper trio (`banner`, `expectZero`, `expectApprox`, `loadConstants`), same four-section structure with identical banner titles, same constants pulled from the same JSON file, same Section 3 hard-coded coefficient expressions (the `dE2 = (deltaC - 9 sigmaStar dKappa)/(27 (1-sigmaStar))` / `dE4 = (5 deltaC - 72 sigmaStar dKappa)/(243 (1-sigmaStar))` system reproduces the Python literally), and the same `deltaCFromTangent = -16 sigmaStar (dR /. dg -> gp dr)` derivation in the same order. Examples of paired transliteration:

- Python (line 60): `g_star = sp.simplify(r - sp.sqrt(1 + r**2) / 2)` ↔ Mathematica (line 47): `gStar = FullSimplify[r - Sqrt[1 + r^2]/2]`.
- Python (lines 106-107): `dE2 = (deltaC - 9 sigma_star * dkappa)/(27*(1 - sigma_star)); dE4 = (5 deltaC - 72 sigma_star * dkappa)/(243*(1 - sigma_star))` ↔ Mathematica (lines 94-95): `dE2 = (deltaC - 9 sigmaStar dKappa)/(27 (1 - sigmaStar)); dE4 = (5 deltaC - 72 sigmaStar dKappa)/(243 (1 - sigmaStar))`.
- Python (line 112): `deltaC_from_tangent = sp.simplify(-16 * sigma_star * dR.subs(dg, gp * dr))` ↔ Mathematica (line 102): `deltaCFromTangent = FullSimplify[-16 sigmaStar (dR /. dg -> gp dr)]`.

Neither engine independently derives the dE2/dE4 expressions, the `-16 sigma_star` projector, or the `Pi = Sigma_0 (1 - S/4)` linearisation from upstream physical premises; both consume them as given inputs and re-test the same identity. The Mathematica check therefore does not add independent algebraic evidence — it confirms only that the Mathematica simplifier reaches the same syntactic conclusion as SymPy on the same expression.

**Why this matters:**
The second-engine policy is meant to catch CAS-specific simplification bugs and human transcription errors by re-deriving the result through a different code path. A line-by-line port cannot catch a coefficient error in the dE2/dE4 system or in the projector `-16 sigma_star`, because both engines see the same literals. If the underlying coefficients are wrong upstream, both engines pass together.

**Required change:**
At least one of the two engines should re-derive the load-bearing pieces from upstream premises rather than retyping the closed-form. Specifically: (a) the dE2/dE4 expressions should be computed in the Mathematica script from the actual canonical-even Galerkin coefficients (or wherever the 27/243 normalisations originate) rather than written as literals; or (b) the `-16 sigma_star` projector that relates `delta R` to `delta C` should be derived in one engine and only the consumed value used in the other, with the Mathematica script independently checking the implication chain `dR = 0 ⇒ dC = 0` via the canonical-even constraint pair rather than via the same `-16 sigma_star` multiplication. The directive flags this for restructuring the `.wl` script in Section 3.

**Verification:**
The fixed `.wl` script should contain at least one Section 3 derivation step that does not appear verbatim (modulo syntax) in the `.py` script — e.g., a Galerkin projection that produces the dE2/dE4 expressions, or an independent computation of `deltaC = -16 sigma_star * dR` from the underlying linearised potential.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:112-113`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:102-103`

**What's wrong:**
The check `tangent motion kills delta C` is `expect_zero("tangent motion kills delta C", -16 sigma_star * dR.subs(dg, gp*dr))`. But the immediately preceding assertion (line 101) already proved `dR.subs(dg, gp*dr) == 0`. Multiplying a verified zero by `-16 sigma_star` is guaranteed to give zero — this assertion cannot fail in any state where the preceding assertion passes, and so it does not independently exercise the comment claim "tangent motion kills delta C." The substantive content (that the canonical-even constraint pair maps a vanishing dR onto a vanishing dC) lives in the dE2/dE4 Solve at lines 103-110, not in this line.

**Why this matters:**
The output banner "tangent motion kills delta C = 0" reads, on its own, as if the script verified a second physical implication. It does not. Anyone reading the transcript may believe both `delta R = 0` and `delta C = 0` were independently verified along the family; in fact only the first is.

**Required change:**
Either (a) replace the body with a substantive check — derive `delta C` from `delta R` via the canonical-even constraints (`dE2 = 0` ⇒ relation between deltaC and dkappa; pair with `dE4 = 0`), and verify that under `dg = gp dr` the resulting deltaC vanishes; or (b) clarify the comment / output label to indicate this line is a trivial consequence of A8 ("delta C = -16 sigma_star * delta R, hence vanishes by A8") rather than presenting it as an independent assertion. Mirror the change in the `.wl` file.

**Verification:**
After fix, the relevant assertion either references the dE2/dE4 system explicitly (option a), or carries an inline comment / banner label explaining the trivial-multiple relationship (option b). Output transcript wording should reflect the choice.

### F3 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py` (mtime 2026-05-11 12:48)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.txt` (mtime 2026-05-11 12:47)

**What's wrong:**
The SymPy saved output is one minute older than the SymPy script. The output content otherwise matches what the current script would produce (same banner titles, same assertion names — including the "Stages 138-139" banner that still appears in the current script). The Mathematica output (13:28) is fresher than its script (13:17), so only the sympy output is stale.

**Why this matters:**
Informational. The verifier will trigger a fresh run; the content currently captured appears consistent with the current script, so this is not blocking.

**Required change:**
A fresh `redteam exec-sympy 157` after F1/F2 land will refresh the output.

**Verification:**
Output mtime > script mtime after the fix iteration.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script, not an independent derivation. See F1 for paired excerpts. Section 3 in particular re-uses the SymPy dE2/dE4 expressions and the `-16 sigma_star` projector verbatim. Section 4's linearised `Pi = Sigma_0 (1 - S/4)` packet is also a direct port. Section 1's `R(g_*) - 1/4 = 0` is an algebraic identity that any engine handles trivially, so transliteration there is unavoidable, but Sections 2-4 are where independence would actually buy redundancy.

## Engine cross-check

Both engines complete with `EXIT_CODE: 0` and report identical assertion names. Residuals agree to working precision:
- `r_F1 radical check`: sympy 1.88e-15, mathematica 1.95e-15.
- `g_* lower-branch check`: sympy 1.73e-16, mathematica 2.16e-16.
- `self-matched traction law`: sympy 1.32e-15, mathematica 1.20e-15.
- `R_can = 1/4`: sympy 4.36e-16, mathematica 4.33e-16.
- `Pi_can identity`: sympy 1.69e-15, mathematica 1.51e-15.
- `dPi_tan at the renormalized canonical point`: both engines yield `-1.162758387542217 dS + 0.832409471081634 dSigma0` (to 30 digits).

Numerical agreement is well within tolerances. The transliteration (F1) is the reason agreement is this tight — both engines are evaluating the same closed-form on the same inputs — but agreement itself is not in dispute.

## Verdict justification

`findings`, not `clean`: the Mathematica script is a port of the SymPy script (F1 medium), and the "tangent motion kills delta C" line is a trivial restatement of the preceding `delta R = 0` check (F2 low). The saved SymPy output is one minute stale (F3 low, informational). No paper_misalignment: the script faithfully exercises the three card checks and the carry-forward numerics match the notes to 12+ digits. No `stop_cold` is warranted — the findings are about engine independence and assertion phrasing, not about the math being wrong. Attacks I tried that failed: I attempted to derive a paper↔script claim mismatch (the card's `delta_perp` vs the script's `delta R`) but the notes use the family parametrisation in which they coincide; I attempted to flag the dE2/dE4 numerators as `hardcoded_result` but they feed a Solve whose outcome (deltaC = dkappa = 0) is the meaningful claim and the values come from the upstream canonical-even decomposition; I attempted to flag the strict-above ordering as too weak a check but it correctly exercises the notes' "renormalized canonical sits above the analyzed canonical" statement.

## Self-test notes

I walked through the dE2/dE4 system mentally: the 2×2 determinant is `1·(-72 sigma_star) - 5·(-9 sigma_star) = -27 sigma_star`, full rank for any nonzero sigma_star, so Solve correctly returns the unique `{deltaC:0, dkappa:0}`. I also confirmed the algebraic identity `R(g_-) = 1/4` reduces to `(1+r^2)/4 / (1+r^2) = 1/4` symbolically — non-tautological because the formula `g_- = r - sqrt(1+r^2)/2` is what's being asserted (substituting any other guess would fail). I confirmed `dR.subs(dg, gp*dr) = 0` is a genuine tangent-vanishing check (the projection `dg = gp dr` lands on the curve `g = g_-(r)`, where `R` is constant `= 1/4`, so `dR = 0` is the real claim, not a tautology). The F2 finding is about the *next* line trivially inheriting that zero, not about A8 itself.
