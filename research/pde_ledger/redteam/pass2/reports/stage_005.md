---
unit_id: 005
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

# Audit unit 005 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_005.tex`
- notes: `(none)` (no files match `notes/stages/moving_throat_pde_stage005_*.md`; the prompt's reading list also declares notes as `(none)`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 32 references stage 005)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.txt`

## What the paper claims

Stage 005 is the covariant projection of the localized parent Maxwell equation onto the brane. The parent law is `∂_N(Z(w)F^{NM}) = μ₀J^M` (eq:stage005-parent-maxwell). For brane-facing components `M=μ`, one applies a normalized observation kernel `W(w)` and integrates in `w`. The `\stagefield{Output}` paragraph states: "Stage~005 exports the projected inhomogeneous law (\ref{eq:stage005-projected-maxwell})--(\ref{eq:stage005-projected-maxwell-expanded}) and projected continuity (\ref{eq:stage005-projected-continuity})." The distinct deliverables are: (1) the projected inhomogeneous Maxwell law `∂_ν∫WZF^{νμ}dw = μ₀∫WJ^μ dw − ∫W∂_w(ZF^{wμ})dw` and its IBP-expanded form `= μ₀∫WJ^μ dw − [WZF^{wμ}]_∂ + ∫(∂_wW)ZF^{wμ}dw` (the explicit boundary-flux + kernel-gradient leakage structure); and (2) the projected charge continuity law `∂_μ∫WJ^μ dw = −[WJ^w]_∂ + ∫(∂_wW)J^w dw`, framed as an exact open-system balance term. The appendix row (line 32) summarizes: "Exact projected inhomogeneous Maxwell equation, boundary flux, kernel-gradient leakage, and projected continuity law." The card carries no numeric constants; the only quantitative content is the structure of the identities. The card also relies on the stage-004 IBP identity (eq:stage004-projection-ibp), which I confirmed exists in `paper/stages/stage_004.tex:28`.

## What the script claims to verify

Both scripts verify, for an explicit Gaussian projection kernel `W(w) = e^{-w²}/√π`, the same family of projected identities the paper exports. The SymPy script (a) prints the generic symbolic derivation (projection commutes with t,x,y,z but not with `∂_w`, where IBP produces a boundary term + a `W'`-weighted leakage term), then (b) runs a concrete-profile audit block that asserts: brane-derivative commutation; the decaying-kernel boundary term vanishes; the projected IBP identity `Pg(∂_w Q) + Pgp(Q) = 0`; a sign-mutant of the IBP that must NOT vanish; the full projected inhomogeneous law `project_bulk_lhs − projected_lhs == 0` with and without the boundary term; a leakage sign-mutant that must NOT vanish; the exact leakage value `−Pgp(w) = 1`; and the projected continuity law `project_bulk_cont − projected_cont == 0`. The Mathematica script independently re-derives the same five claims (M1–M5) with *different* concrete profiles (`Sin[t]`, `Cos[t]`, higher even powers of `w`, `Q = w·e^{−w²/4}`) plus a sign-mutant guard for each. The docstrings on both engines match the paper's stated deliverables verbatim in intent (projected inhomogeneous law, boundary discharge, kernel-gradient leakage, projected continuity).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1a) Projected inhomogeneous law, leakage form `−∫∂_w(ZF^{wμ})` | sympy generic `projected_eq` / `projected_eq_decay` (lines 120–143); concrete `assert_zero` "concrete projected inhomogeneous law" (line 262); math M3 (lines 89–96) | match |
| (1b) IBP-expanded boundary + kernel-gradient form `−[WZF^{wμ}]_∂ + ∫(∂_wW)ZF^{wμ}` | sympy `proj_dw` (line 96), `projected_eq` keeps `boundary(W·Gw) − Pwprime(Gw)` (lines 125–126); concrete boundary term asserted zero (line 260) and "...with boundary" (line 261); math M2 IBP + M3 boundary (lines 53–73, 83–88) | match |
| (1c) Kernel-gradient leakage is nonzero / explicit value | sympy "concrete transverse leakage value" `leakage − 1 == 0` and "concrete transverse leakage" nonzero (lines 277–279); math M4 `leakage == 1` + nonzero guard (lines 107–127) | match |
| (2) Projected continuity (open-system charge law) | sympy `proj_cont` / `proj_cont_decay` (lines 171–191); concrete "concrete projected continuity law" (line 303); math M5 (lines 129–158) | match |
| Brane-derivative commutation (premise of (1)/(2)) | sympy "projection commutes with brane derivative" (line 226); math M1 (lines 37–51) | match |

No script-side check is orphaned; no paper deliverable is unmatched. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 226 | `Pg(∂_t Φ) − ∂_t Pg(Φ) == 0` | commutation premise | yes |
| A2 | sympy | 230 | `boundary_value(Wg·Q) == 0` | (1b) boundary discharge | yes |
| A3 | sympy | 231 | `Pg(∂_w Q) − (boundary − Pgp(Q)) == 0` | (1b) IBP-with-boundary | yes |
| A4 | sympy | 232 | `Pg(∂_w Q) + Pgp(Q) == 0` | (1b) decay-form IBP | yes |
| A5 | sympy | 233 | `Pg(∂_w Q) − (boundary + Pgp(Q)) ≠ 0` (mutant) | (1b) sign guard | yes |
| A6 | sympy | 260 | `boundary_value(Wg·Gwc) == 0` | (1b) boundary discharge | yes |
| A7 | sympy | 261 | `project_bulk_lhs − projected_lhs_with_boundary == 0` | (1a)/(1b) full law w/ boundary | yes |
| A8 | sympy | 262 | `project_bulk_lhs − projected_lhs == 0` | (1a) full inhomogeneous law | yes |
| A9 | sympy | 263–275 | leakage-sign mutant ≠ 0 | (1a) sign guard | yes |
| A10 | sympy | 278 | `(−Pgp(Gwc)) − 1 == 0` | (1c) exact leakage value | yes |
| A11 | sympy | 279 | `(−Pgp(Gwc)) ≠ 0` | (1c) nonzero leakage | yes |
| A12 | sympy | 301 | `boundary_value(Wg·Jwc) == 0` | (2) boundary discharge | yes |
| A13 | sympy | 302 | `project_bulk_cont − (projected_cont + boundary) == 0` | (2) continuity w/ boundary | yes |
| A14 | sympy | 303 | `project_bulk_cont − projected_cont == 0` | (2) projected continuity | yes |
| M1 | math | 42–51 | `projW[D[Φ,t]] − D[projW[Φ],t] == 0` + mutant | commutation premise | yes |
| M2 | math | 57–73 | boundary==0; `projW[D[Q,w]]+projWPrime[Q]==0`; mutant | (1b) IBP | yes |
| M3 | math | 85–105 | boundary==0; full inhomogeneous law ==0; mutant | (1a) full law | yes |
| M4 | math | 113–127 | `leakage==1`; mutant; nonzero guard | (1c) exact leakage | yes |
| M5 | math | 138–158 | boundary==0; projected continuity ==0; mutant | (2) continuity | yes |

Every row traces to a specific paper deliverable. Every `assert_zero`/`==0` of substance is paired with either a sign-mutant `assert_nonzero`/`=!=0` guard or a distinct nonzero-value check, so none is tautological.

## Findings

None. Each attack below failed to break the scripts.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. Three points of evidence:
- **Different concrete profiles.** SymPy uses `Phi = (t²+x)(w²+1)`, `Q = w³+w`, `Gwc = w`, `G0c = t(w²+1)`, etc. Mathematica uses `Phi = (Sin[t]+x²)(w²+3)`, `Q = w·Exp[−w²/4]`, `Gw = w + w³/3`, `G0 = Cos[t](w²+2)`, `Gy = y(w⁴+1)`, etc. These are independently chosen witnesses, not the same numbers rewritten.
- **Different boundary-test design.** SymPy tests boundary decay with a polynomial profile `w³+w` multiplied by the Gaussian kernel; Mathematica deliberately tests with a non-Gaussian decaying factor `w·Exp[−w²/4]` whose product with the kernel decays at a different rate — a genuinely separate witness for the same boundary-discharge claim.
- **Different leakage witness vs. value.** Both arrive at the same analytic constant `−Pgp(w) = 1`, but M4 isolates `leakageGw = w` and computes it from the kernel directly; this is a derived numeric outcome (`1`), not a hardcoded echo of the SymPy block. The agreement on the value `1` is a true cross-engine confirmation, not a copy.

The Mathematica script also adds mutant residuals printed as concrete numbers (`7 Cos[t]`, `16/(5√5)`, `3`, `−2`, `3`) which differ from anything in the SymPy script, confirming independent algebra.

## Engine cross-check

Both engines emit `STATUS: PASS`. The load-bearing identities agree:
- Commutation residual = 0 (both).
- Projected IBP identity residual = 0 (both); both fire a nonzero sign-mutant (SymPy A5 `assert_nonzero`; Mathematica M2 mutant = `16/(5√5)` ≠ 0).
- Projected inhomogeneous law residual = 0 (both); both fire a leakage sign-mutant (SymPy A9; Mathematica M3 mutant = `3` ≠ 0).
- Exact transverse leakage value = `1` (SymPy A10 asserts `leakage − 1 == 0`; Mathematica M4 output `M4 leakage cross-check = 1`). Verified by hand: `Pgp(w) = ∫ W'(w)·w dw = ∫ (−2w/√π)e^{−w²}·w dw = (−2/√π)·(√π/2) = −1`, so `−Pgp(w) = +1`. Both engines agree with the hand computation.
- Projected continuity residual = 0 (both); both fire a sign-mutant (SymPy A-mutant at line 263 family applies to inhomogeneous law; Mathematica M5 mutant = `3` ≠ 0).

`engines_agree: true`.

## Verdict justification

`clean`. I read the paper card, the appendix row, and confirmed there are no stage-005 notes; the card's only cross-stage dependency (eq:stage004-projection-ibp) is real. The two deliverables the card exports — the projected inhomogeneous Maxwell law (with its boundary-flux + kernel-gradient leakage structure and exact leakage constant) and the projected continuity law — are each exercised by both engines with independently chosen concrete profiles and matching sign-mutant guards. Attacks tried and failed: (a) **tautology** — every substantive `==0` is paired with a sign-mutant `≠0` guard or a nonzero-value check, so none is constructed to pass unconditionally; (b) **boundary-term cheat** — I confirmed the boundary terms genuinely vanish for the chosen decaying profiles (`w³+w`, `w·e^{−w²/4}`, `w`, `w+w³/3` times the Gaussian all → 0 at ±∞), so the "with boundary" and "decay form" assertions are not silently identical; (c) **leakage value** — hand computation of `−Pgp(w)=1` matches both engines, so M4/A10 is not a hardcoded number checked against itself (the value is derived from the Gaussian kernel in-script); (d) **derivative-zero trap** — the inhomogeneous-law and continuity checks differentiate `Pg(G0c)` w.r.t. `t` etc., and each projected field genuinely depends on its differentiation variable (`G0c = t(w²+1)` depends on `t`, `Gxc = x·w²` depends on `x`, etc.), so no `assert_nonzero` is trivially satisfied and no `assert_zero` passes by a vanishing derivative; (e) **symbol-domain** — `t,x,y,z` real, `mu0,xi` nonzero; the Gaussian integrals are over real `w` and converge; no assumption masks a branch error. Outputs are fresh (sympy/math `.txt` mtimes May 21 ≥ script mtimes May 4 / May 21). The script verifies exactly the paper's claim, no more, no less.

## Self-test notes

- **Variable independence:** Confirmed every `sp.diff(Pg(...), VAR)` / `D[projW[...], var]` differentiates a projected field that actually depends on `VAR` (`G0c`/`G0` on `t`, `Gxc`/`Gx` on `x`, `Gyc`/`Gy` on `y`, `Gzc`/`Gz` on `z`; continuity J-fields likewise). No identically-zero derivative — the derivative-zero trap from earlier units is absent here.
- **Symmetry/parity:** The leakage integral `∫ W'(w)·w dw` is (odd)·(odd)=even integrand → nonzero (=−1), consistent with the asserted nonzero leakage; boundary integrands like `w·e^{−w²}` are odd but the boundary *value* is a limit difference at ±∞ (both 0), correctly asserted zero. Parity of M3/M5 mixed even/odd `w`-powers (e.g. `w+w³/3`) does not make any projected term spuriously vanish.
- **Trivial-case pre-check:** Mentally substituting the concrete profiles reproduces the saved outputs (M4=1, mutant residuals 7Cos[t], 16/(5√5), 3, −2, 3), all nonzero where required and 0 where required. No directive is written because `findings_count = 0`.

## Value Reconciliation (pass-2 augmentation)

The deliverables of stage 005 are **structural identities**, not numeric constants. The only genuine numeric result value the scripts emit is the exact Gaussian transverse-leakage constant `−Pgp(w) = 1`. All other emitted quantities are either (i) symbolic boxed identities that the paper carries verbatim as equations, or (ii) verification scaffolding (residuals, mutant residuals, PASS flags, boundary values asserted to 0).

### Deliverable-level reconciliation

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Projected inhomogeneous Maxwell law `∂_μ Proj_W[ZF^{μν}] + Boundary[WZF^{wν}] − Proj_{W'}[ZF^{wν}] + Proj_W[Γ^ν] = μ₀ Proj_W[J^ν]` | py lines 150–154 / out lines 37–41; wl M3 lines 89–96 / out lines 12–17 | `stage_005.tex:20–35` (eq:stage005-projected-maxwell, eq:stage005-projected-maxwell-expanded) | MATCH |
| Boundary-decay form `∂_μ Proj_W[ZF^{μν}] − Proj_{W'}[ZF^{wν}] + Proj_W[Γ^ν] = μ₀ Proj_W[J^ν]` | py lines 199–201 / out line 66; wl M3 decay structure | `stage_005.tex:20–26` (eq:stage005-projected-maxwell) | MATCH |
| Projected charge continuity `∂_t Proj_W[J^0] + ∂_a Proj_W[J^a] = Proj_{W'}[J^w]` | py lines 203–204 / out line 70; wl M5 lines 142–149 / out lines 24–29 | `stage_005.tex:39–43` (eq:stage005-projected-continuity) | MATCH |
| Non-commutation / IBP identity `Proj_W[∂_w Q] = Boundary[WQ] − Proj_{W'}[Q]` | py line 96 / out line 12; wl M2 lines 61–67 / out lines 6–11 | `stage_005.tex:27–35` (via eq:stage004-projection-ibp reference + expanded form) | MATCH |
| Exact transverse leakage constant `−Proj_{W'}[w] = 1` (for Gaussian kernel) | py line 278 / (asserted, value implicit); wl M4 line 110 / out line 18 (`M4 leakage cross-check = 1`) | not in `.tex`; no `.md` notes exist | INTERNAL (see note) |

**Note on the leakage constant `1`:** This `1` is a *witness-specific* numeric outcome of choosing the concrete Gaussian kernel `W=e^{−w²}/√π` and the concrete profile `G^{wν}=w` for the audit. It is not one of the stage's stated deliverables — the paper deliverable is the *structural* leakage term `∫(∂_wW)ZF^{wμ}dw` (which the card states symbolically), not the numeric value `1` for one arbitrary test profile. The card is a terse symbolic card and there are no stage-005 notes by design (prompt declares notes `(none)`). Per the augmentation Guards ("a terse `.tex` card legitimately omits intermediate quantities"; mark MISSING only for *deliverables* absent from both `.tex` and `.md`), this is INTERNAL scaffolding, not a MISSING-DELIVERABLE. No finding.

### INTERNAL items (no finding expected)

Pass/fail flags (`STATUS: PASS`); residual-near-zero values (`M1..M5 residual = 0`, all `assert_zero` residues); mutant residuals (`7 Cos[t]`, `16/(5√5)`, `3`, `−2`, `3`); boundary values asserted to 0 (`M2/M3/M5 boundary value = 0`, sympy `boundary_Q`, `boundary_Gw`, `boundary_Jw`); the witness-specific leakage constant `1` (justified above); intermediate concrete profiles (`Phi`, `Q`, `G0c`/`G0`, etc.).

reconciliation: complete; 4 deliverable values checked, 0 misaligned (1 numeric witness value classified INTERNAL with justification).
