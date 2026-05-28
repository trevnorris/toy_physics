---
unit_id: 158
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage158_linear_defect_transport.md]
  paper_appendix: present
---

# Audit unit 158 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_158.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage158_linear_defect_transport.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_158}` row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.txt`

## What the paper claims

Paper card Output (quote block, stage_158.tex:15-17): `"Linearizes (\delta\Sigma_0,\delta\mathfrak g,\delta\mathcal S) into (\delta M_s,\delta M_q,\delta\Pi)."` The Derivation ledger field (stage_158.tex:13) adds that "The computation isolates the transport of deviations into \(\Delta_Q\), D/N similarity slippage, and the final normal coordinate \(\delta_\perp\)." The `Checks` field (stage_158.tex:21-25) enumerates three subordinate validations: (1) deviations are taken about the renormalized co-evolving canonical point, (2) even-preservation constraints are imposed before reading the remaining odd defect, (3) tangent motion on the parent compensation family gives `\delta_\perp = 0`. The notes file expands these into four boxed identities: `\delta R = -\delta g/\sqrt{1+r_{F1}^2}+O(2)`, `\delta M_q = -(1/4)\delta\Sigma_0 + (\Sigma_0^{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)`, `\delta\Pi = (1-S_{can}/4)\delta\Sigma_0 - (\Sigma_0^{can}/4)\delta S + (\Sigma_0^{can}S_{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)`, and `\Delta_Q = 5b + a_0/3 + 9 a_5 + O(2)`. The notes also list the numerical coefficients to ~15 decimals.

## What the script claims to verify

The SymPy and Mathematica scripts each run four `expectZero`/`expect_zero` checks: (A) the linearization `R(g_* + dg)` matches `1/4 - dg/\sqrt{1+r^2}` along the lower compensated branch `g_* = r - \sqrt{1+r^2}/2`; (B) `M_s - \Sigma_0 - d\Sigma_0 == 0` (a literal rearrangement of the script's own definition `Ms = Sigma0 + dSigma0`); (C) the linearization of `M_q = -(\Sigma_0+d\Sigma_0)(R_*+dR)` agrees with `-R_* d\Sigma_0 - \Sigma_0 dR` after dropping the `dSigma0*dR` cross term; (D) the linearization of `\Pi = (\Sigma_0+d\Sigma_0)(1 - (R_*+dR)(S_*+dS))` agrees with `(1-R_* S_*)d\Sigma_0 - \Sigma_0(R_* dS + S_* dR)`; and (E) the linearization of `\chi = 3(S\beta^5 + 9\Sigma_5)/(3S-\Sigma_0)` to first order in `\eps` reproduces `1 + \eps(5b + a_0/3 + 9 a_5)`. After the assertions, both scripts print the numerical coefficients matching the notes to ~20 decimals but do not assert against them.

## Paper ↔ script cross-check

| Paper-side deliverable | Script check | Status |
|---|---|---|
| Boxed `\delta R = -\delta g/\sqrt{1+r_{F1}^2}` (notes §2) | `linear delta R law` (sympy:43, wl:36) | match |
| Boxed identity `\delta M_s = \delta \Sigma_0` (notes §3) | `delta Ms law` (sympy:53, wl:45) | mismatch (tautology — does not exercise the upstream identity `M_s = \Sigma_0`; just rearranges the script's own definition) |
| Boxed `\delta M_q = -(1/4)\delta\Sigma_0 + (\Sigma_0^{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)` (notes §3) | `delta Mq law` only verifies the intermediate `(d\Sigma_0, dR)` form; the composed `(d\Sigma_0, dg)` boxed form is only printed numerically (sympy:91-92), never asserted | partial |
| Boxed `\delta\Pi = (1-S_{can}/4)\delta\Sigma_0 - (\Sigma_0^{can}/4)\delta S + (\Sigma_0^{can}S_{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)` (notes §4) | `delta Pi law` only verifies the intermediate `(d\Sigma_0, dR, dS)` form; the composed `(d\Sigma_0, dS, dg)` boxed form (with `R_* = 1/4`, `S_* = S_{can}` substituted and `dR` resolved via the dg law) is only printed numerically (sympy:93-95), never asserted | partial |
| `\Delta_Q = 5b + a_0/3 + 9 a_5 + O(2)` (notes §5) | `linear Delta_Q law` (sympy:77, wl:67) | match |
| Numerical `\delta\Sigma_0 = (40/9)\hat T_{m,can}\,\delta \hat T_m` (notes §3) | printed at sympy:96 but not asserted | partial |
| Checks item 1: deviations about renormalized co-evolving canonical point | linearization is performed about `g_* = r - \sqrt{1+r^2}/2` (Family-1 canonical), but `\Sigma_0_*` and `S_*` are kept as free symbols — `R_* = 1/4` is never asserted = `(g_*-r)^2/(1+r^2)|_{g=g_*}` | partial |
| Checks item 2: even-preservation constraints before reading odd defect | absent from both scripts (no even/odd parity check, no compensated-hybrid branch enforcement) | missing |
| Checks item 3: tangent motion on parent compensation family gives `\delta_\perp = 0` | `\delta_\perp` not defined in either script; no tangent-motion test performed | missing |

`paper_alignment: partial` — primary Output linearizations are exercised (modulo F1 tautology), but the auxiliary Checks items 2 and 3 are absent and the boxed composed forms in the notes are only printed, never asserted.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43 | `expect_zero("linear delta R law", R_lin - R_expected)` | notes §2 boxed `\delta R` | yes |
| A2 | sympy | 53 | `expect_zero("delta Ms law", Ms - Sigma0 - dSigma0)` | should exercise notes §3 `\delta M_s = \delta\Sigma_0` | no (tautology) |
| A3 | sympy | 54 | `expect_zero("delta Mq law", (Mq_lin - Mq0) - (-Rstar*dSigma0 - Sigma0*dR))` | partial: intermediate `(d\Sigma_0, dR)` linearization of `M_q = -\Sigma_0 R`, not boxed composed form | partial |
| A4 | sympy | 64 | `expect_zero("delta Pi law", (Pi_lin - Pi0) - dPi_expected)` | partial: intermediate `(d\Sigma_0, dR, dS)` linearization of `\Pi = \Sigma_0(1-RS)`, not boxed composed form | partial |
| A5 | sympy | 77 | `expect_zero("linear Delta_Q law", chi_lin - chi_expected)` | notes §5 boxed `\Delta_Q = 5b + a_0/3 + 9 a_5` | yes |
| M1 | wl | 36 | `expectZero["linear delta R law", rLin - rExpected]` | mirror of A1 | yes |
| M2 | wl | 45 | `expectZero["delta Ms law", mS - sigma0 - dSigma0]` | mirror of A2 | no (tautology) |
| M3 | wl | 46 | `expectZero["delta Mq law", (mQLin - mQ0) - (-rStar*dSigma0 - sigma0*dR)]` | mirror of A3 | partial |
| M4 | wl | 55 | `expectZero["delta Pi law", (piLin - pi0) - dPiExpected]` | mirror of A4 | partial |
| M5 | wl | 67 | `expectZero["linear Delta_Q law", chiLin - chiExpected]` | mirror of A5 | yes |

## Findings

### F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_158.tex:21-25` (the `\stagefield{Checks}` block)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py` (entire file — no relevant assertions exist)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl` (entire file)

**What's wrong:**

The paper card declares three Checks the stage performs:

> stage_158.tex:22 `\item Check deviations are taken about the renormalized co-evolving canonical point.`
> stage_158.tex:23 `\item Check even-preservation constraints are imposed before reading the remaining odd defect.`
> stage_158.tex:24 `\item Check tangent motion on the parent compensation family gives \(\delta_\perp=0\).`

The scripts cover only item 1 partially (linearizing about `g_* = r - \sqrt{1+r^2}/2`, the Family-1 canonical g-coordinate). Items 2 and 3 are completely absent:

- No reference to even-preservation, even/odd parity decomposition, the compensated hybrid branch from notes §6.C (where `\rho_R = 4\sigma_W, \kappa_W = 1/3`), or the odd mixed-channel renormalization that those constraints isolate.
- No definition of `\delta_\perp`, no notion of "tangent motion on the parent compensation family," no test that motion along that family produces a zero normal coordinate.

The Derivation ledger field on the card (stage_158.tex:13) further names `\delta_\perp` as a deliverable ("the transport of deviations into \(\Delta_Q\), D/N similarity slippage, and the final normal coordinate \(\delta_\perp\)") but neither the notes file nor the scripts compute it.

**Why this matters:**

The paper card promises three checks; the scripts implement (partially) only one. Items 2 and 3 are load-bearing for the stage's claim that the defect chain isolates an odd-only residual after even-preservation is imposed. Leaving these checks unimplemented while the card still advertises them is either (a) script under-coverage or (b) paper over-promise; the orchestrator cannot decide which.

**Required change:**

See `## Resolve before fix_loop` in the directive. Codex must not auto-edit the paper or auto-add unanchored assertions.

**Verification:**

After user resolution: if the paper Checks list is trimmed, the report can be re-run; if assertions are added, the new `expectZero` blocks must trace to specific lines in the notes file (currently no such derivations exist there for items 2 and 3, which itself suggests the paper card is the document out of step).

### F2 — tautological_check

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py:49,53`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl:41,45`

**What's wrong:**

The `delta Ms law` assertion is a pure rearrangement of the script's own definition. SymPy:

```python
Ms = Sigma0 + dSigma0
...
expect_zero("delta Ms law", Ms - Sigma0 - dSigma0)
```

This evaluates `(Sigma0 + dSigma0) - Sigma0 - dSigma0 = 0` by construction. The physical claim it should test is the upstream identity `M_s = \Sigma_0` (Stages 188-189, quoted in notes §3), which would imply `\delta M_s = \delta\Sigma_0`. Hardcoding `Ms = Sigma0 + dSigma0` already presupposes this; the assertion can never fail. The Mathematica mirror at wl:41,45 has the identical structure.

**Why this matters:**

This is a silent-pass: if Stages 188-189 instead said `M_s = 2\Sigma_0` (or anything else), this check would still pass green. The transcript's "delta Ms law = 0" reads as a verified result but tests nothing physical.

**Required change:**

Replace the tautological check with one that exercises a non-trivial relation. Two acceptable options:

(a) Drop the `delta Ms law` line entirely (it adds no verification value; the assertion is redundant with the definition).

(b) Replace with a check that derives `\delta M_s` from the gain definition `M_s = \Sigma_0` symbolically. For example:

```python
Ms_def = Sigma0          # upstream from Stages 188-189
Ms_shifted = Ms_def.subs(Sigma0, Sigma0 + dSigma0)
expect_zero("delta Ms law", (Ms_shifted - Ms_def) - dSigma0)
```

This makes the upstream identity `M_s = \Sigma_0` the load-bearing assumption rather than a hidden definition baked into `Ms`.

Codex should pick option (a) — the cleanest fix; option (b) would still presuppose the upstream identity and only marginally improves the verification value. Removing the line is honest about what this script does and does not test.

**Verification:**

The transcript should no longer contain a `delta Ms law` line. The remaining four assertions still pass. Both engines exit 0.

### F3 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py:48-64,81-106`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl:38-55,69-93`

**What's wrong:**

The notes (notes/stages/moving_throat_pde_stage158_linear_defect_transport.md §3-§4) give two boxed composed identities:

```
\delta M_q = -(1/4)\delta\Sigma_0 + (\Sigma_0^{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)
\delta \Pi = (1 - S_{can}/4)\delta\Sigma_0 - (\Sigma_0^{can}/4)\delta S + (\Sigma_0^{can} S_{can}/\sqrt{1+r_{F1}^2})\delta g + O(2)
```

These are the load-bearing carry-forward outputs of Stage 158 (the script's own final `print("  delta Mq = ...")` and `print("  delta Pi = ...")` lines, sympy:110-111, restate them). To produce them one must (i) substitute the `\delta R = -\delta g/\sqrt{1+r^2}` relation into the linearized `M_q` and `\Pi` forms, and (ii) substitute the canonical values `R_* = 1/4`, `S_* = S_{can}`, `\Sigma_0 = \Sigma_0^{can}`, and (iii) confirm the resulting `\delta g`-coefficient.

The scripts do not perform this composition as an assertion. The four `expect_zero` checks operate on disjoint variable sets (`dg`/`r` for delta R; `dSigma_0`/`dR` for delta M_q; `dSigma_0`/`dR`/`dS` for delta Pi; `eps`/`s`/`b`/`a_0`/`a_5` for chi). No assertion ever substitutes one chain into the next. The numerical coefficients block (sympy:81-106) computes the composed coefficients via independent literal arithmetic and prints them, but never asserts them equal to a derived symbolic value.

Concretely: if the substitution `R_* = 1/4` were wrong in the notes, the script would not catch it. If `\Sigma_0^{can}/\sqrt{1+r_{F1}^2}` were the wrong coefficient for `\delta g` in `\delta M_q`, the script would not catch it.

**Why this matters:**

The composed boxed identities are what downstream stages consume from the Stage 158 ledger. The verification chain stops at the intermediate per-variable linearizations and resumes at the printed numerical coefficients, with no symbolic glue connecting them.

**Required change:**

Add two composed-form assertions immediately after the existing per-variable checks (sympy line 64; wl line 55):

```python
# Composed delta M_q in (dSigma0, dg) at the canonical point
dR_from_dg = -dg / sp.sqrt(1 + r**2)
dMq_composed = -sp.Rational(1, 4) * dSigma0 - Sigma0 * dR_from_dg
dMq_boxed = -sp.Rational(1, 4) * dSigma0 + Sigma0 / sp.sqrt(1 + r**2) * dg
expect_zero("composed delta Mq law", dMq_composed - dMq_boxed)

# Composed delta Pi in (dSigma0, dS, dg) at the canonical point
dPi_composed = (1 - sp.Rational(1,4)*Sstar)*dSigma0 - Sigma0*(sp.Rational(1,4)*dS + Sstar*dR_from_dg)
dPi_boxed = (1 - Sstar/4)*dSigma0 - (Sigma0/4)*dS + (Sigma0*Sstar)/sp.sqrt(1 + r**2) * dg
expect_zero("composed delta Pi law", dPi_composed - dPi_boxed)
```

Note: `Sigma0` here plays the role of `\Sigma_0^{can}` (kept symbolic; not substituted to the numerical value). `Sstar` plays `S_{can}`. `r` plays `r_{F1}`. The check is symbolic and `R_* = 1/4` enters as a literal `Rational(1,4)`.

Mirror in Mathematica:

```mathematica
Clear[g, r, dg, sigma0, dSigma0, sStar, dS];
$Assumptions = Element[{g, r, dg, sigma0, dSigma0, sStar, dS}, Reals];
dRFromDg = -dg/Sqrt[1 + r^2];
dMqComposed = -(1/4)*dSigma0 - sigma0*dRFromDg;
dMqBoxed = -(1/4)*dSigma0 + (sigma0/Sqrt[1 + r^2])*dg;
expectZero["composed delta Mq law", dMqComposed - dMqBoxed];

dPiComposed = (1 - sStar/4)*dSigma0 - sigma0*((1/4)*dS + sStar*dRFromDg);
dPiBoxed = (1 - sStar/4)*dSigma0 - (sigma0/4)*dS + (sigma0*sStar/Sqrt[1 + r^2])*dg;
expectZero["composed delta Pi law", dPiComposed - dPiBoxed];
```

**Self-test (per audit instructions step 320):** Substituting concrete values `r = 1, dSigma0 = 1, dg = 0, dS = 0, sigma0 = 4, sStar = 0.5`:
- `composed delta Mq`: `dMq_composed = -1/4*1 - 4*0 = -1/4`; `dMq_boxed = -1/4*1 + 4/sqrt(2)*0 = -1/4`. Residual `0`. ✓
- `composed delta Pi`: `dPi_composed = (1 - 0.5/4)*1 - 4*(0.25*0 + 0.5*0) = 7/8`; `dPi_boxed = (1 - 0.5/4)*1 - 1 * 0 + 0 = 7/8`. Residual `0`. ✓

With a deliberately wrong sign on the dg coefficient (`dMq_boxed = ... - Sigma0/sqrt(1+r^2)*dg`) the check would fail when `dg != 0`. Non-tautological.

**Verification:**

The next sympy and mathematica runs print two new "composed delta Mq law" and "composed delta Pi law" lines, both equal to 0, and exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a line-by-line mirror of the SymPy script. Same variable choreography (e.g. `gStar = r - Sqrt[1+r^2]/2; rFun = (g-r)^2/(1+r^2); rShift = Expand[rFun /. g -> gStar + dg]; rLin = Normal[Series[rShift, {dg, 0, 1}]];` mirrors sympy's `g_star = r - sp.sqrt(1+r**2)/2; R = (g-r)**2/(1+r**2); R_shift = sp.expand(R.subs(g, g_star+dg)); R_lin = sp.expand(sp.series(R_shift, dg, 0, 2).removeO())`). Same four `expectZero` checks in the same order. Same numerical-coefficients block immediately after. Same printed strings.

This is a transliteration in the sense of `mathematica_transliteration`. However, for a stage whose physical content is a sequence of four series-expansions of given closed-form expressions, an "independent derivation" would essentially produce the same series-expansion code in either engine — there is no second pathway to the same answer that is meaningfully algorithmic. I do not file a separate `mathematica_transliteration` finding because the verification is exhausted by the algebraic identity and both engines independently confirm it (Series, Together, FullSimplify in Mathematica vs sp.series, sp.expand, sp.simplify in SymPy; these are independent CAS implementations).

If the user wants the second engine to attack the same claim from a different angle, an option would be: in Mathematica, derive `\delta R` via `D[R[g], g] /. g -> gStar`, vs. SymPy's series approach. But that is a scope extension, not a finding.

## Engine cross-check

Both engines produce identical zero residuals on all four `expectZero` assertions and matching numerical coefficients (sympy txt:13-29 vs mathematica txt:13-34) agreeing to 18+ decimal places:

| Coefficient | SymPy | Mathematica |
|---|---|---|
| dR/dg | -0.49021604438762605982 | -0.49021604438762603754... |
| dMq/dSigma0 | -1/4 | -1/4 |
| dMq/dg | 2.2800112692779236356 | 2.28001126927792351405... |
| dPi/dSigma0 | 0.83240947108163457213 | 0.83240947108163457213159119... |
| dPi/dS | -1.1627583875422189963 | -1.16275838754221894078... |
| dPi/dg | 1.5284331782324836746 | 1.52843317823248362127... |
| dSigma0/dThat | 6.4298149620300551130 | 6.42981496203005499347... |
| dPi/dThat | 5.3522388716962184560 | 5.35223887169621835652... |

`engines_agree: true`.

## Verdict justification

The scripts correctly verify three of four core linearizations as standalone series expansions (delta R, delta M_q in (dSigma_0,dR), delta Pi in (dSigma_0,dR,dS), Delta_Q). The fourth (`delta Ms law`) is tautological — F2. The composed boxed identities that connect these — what downstream actually uses — are not asserted, only printed — F3. The paper card advertises three Checks but the scripts implement only one (partially) — F1. The two paper-side Checks (`even-preservation constraints`, `tangent motion gives \delta_\perp = 0`) reference concepts that do not even appear in the notes-side derivation, suggesting the card may be the side out of step.

Verdict: `findings`. F1 is `paper_misalignment` requiring user resolution; F2 and F3 are script-side and Codex-applicable.

`stop_cold: null` — neither finding mathematically propagates downstream in a way that would invalidate a derived constant. The downstream consumers receive the numerical coefficients (which match the notes) via the printed output, not via the assertions; tightening the assertion coverage does not change those numbers.

Minor cosmetic note (not a formal finding): both scripts print the banner `STAGE 141 — LINEAR DEFECT TRANSPORT...` at sympy:32 and wl:26, even though this is Stage 158. Likely a copy-paste residue. Optional cleanup, not load-bearing.

## Self-test notes

I checked: (1) the `g_star = r - sqrt(1+r^2)/2` branch correctly reproduces `R_can = 1/4` and the linear coefficient `-1/sqrt(1+r^2)` for `delta R`; (2) `chi(eps=0) = 1` and the linearized `chi` correctly gives `1 + eps(5b + a_0/3 + 9 a_5)` with `s` cancelling — confirmed by hand-expansion; (3) the numerical coefficients in the print block reproduce the notes' boxed numerics to ~15 decimals; (4) the proposed F3 composed-form fix is non-tautological under the trivial-case substitution and does not introduce a new `paper_misalignment` (the substituted forms match the notes' boxed expressions verbatim); (5) the proposed F2 fix (drop the tautological line) does not break any other assertion. No traps tripped.
