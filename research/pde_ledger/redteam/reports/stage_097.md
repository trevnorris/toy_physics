---
unit_id: 097
batch: IV.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage097_single_normalization_defect.md]
  paper_appendix: present
---

# Audit unit 097 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_097.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage097_single_normalization_defect.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only includes `\input{stages/stage_097}` at line 1228 — no narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.txt`

## What the paper claims

The stage card's body equation (quote block, `stage_097.tex:15-17`) states the bottom-line claim verbatim:

> "The isotropic passive/outgoing one-pole branch reduces to \(N_Q=\overline K_0/\overline K_0^{\rm target}\)."

The notes (sections 1–4) expand this into four sub-identities: (i) `Kbar_2 = Kbar_0/(4 Omega_Q^2)` and `Kbar_4 = Kbar_0/(4 Omega_Q^4)` from the canonical-invariant form; (ii) `Gammabar_5 = 9 Kbar_2^(5/2)/Kbar_0^(3/2) = 9 Kbar_0/(32 Omega_Q^5)`; (iii) `Kbar_0^target = 64 G Omega_Q^5/(45 c^5) = 54 G c_s^5/(5 a^5 c^5)` after `Omega_Q = 3 c_s/(2 a)`, with `Gammabar_5^target = 2G/(5 c^5)`; (iv) `R_0 = R_2 = R_4 = R_5 = N_Q - 1`. In addition, the card's `\stagefield{Checks}` block (lines 21–25) demands three checks: a static limit `eps_2=eps_4=0` returning `c_pole=1/4`; `l=0` and `l=2` orthogonality "before applying the geometry firewall"; and the preservation of the minimal-module hypothesis on any support/source success claim.

## What the script claims to verify

Both engines verify four identities: (A1) `Gammabar_5 = 9 K2^(5/2)/K0^(3/2)` simplifies to `9 K0/(32 Omega^5)` under `K2 := K0/(4 Omega^2)` and `Omega, K0` positive; (A2) `K0_target |_{Omega = 3 c_s/(2 a)} = 54 G c_s^5/(5 a^5 c^5)`; (A3) `Gammabar_5_target = 2G/(5 c^5)`; (A4) substituting `K0 = N_Q * K0_target` makes all four `R_i := X / X_target - 1` reduce to `N_Q - 1`. Neither engine tests anything labelled `eps_2`, `eps_4`, `c_pole`, `l=0`/`l=2` orthogonality, or a "minimal-module hypothesis" carry-forward.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage | Status |
|---|---|---|
| `Kbar_2 = Kbar_0/(4 Omega^2)`, `Kbar_4 = Kbar_0/(4 Omega^4)` | sympy lines 13–14, mathematica lines 33–34 (definitions) | match (definitional carry-in) |
| `Gammabar_5 = 9 Kbar_0/(32 Omega^5)` | sympy line 17, mathematica line 41 | match |
| `K0_target = 64 G Omega^5/(45 c^5) → 54 G c_s^5/(5 a^5 c^5)` | sympy line 23, mathematica line 49 | match |
| `Gammabar_5_target = 2G/(5 c^5)` | sympy line 28, mathematica line 54 | match |
| `R_0 = R_2 = R_4 = R_5 = N_Q - 1` | sympy lines 37–40, mathematica lines 66–69 | match |
| Card "Check" #1: static limit `eps_2 = eps_4 = 0` → `c_pole = 1/4` | (none) | missing |
| Card "Check" #2: `l=0`, `l=2` orthogonality before geometry firewall | (none) | missing |
| Card "Check" #3: support/source success carries minimal-module hypothesis | (none) | missing |

Dominant pattern: the body-equation claim is matched faithfully; the three line-item "Checks" enumerated in the card are not exercised. Setting `paper_alignment: partial`.

## Assertion inventory

| #  | Script        | Line | Form                                                              | Exercises which paper claim?        | Anchored to claim? |
|----|---------------|------|-------------------------------------------------------------------|-------------------------------------|--------------------|
| A1 | sympy         | 17   | `assert sp.simplify(Gamma5 - Gamma5_expected) == 0`               | notes §2 (Gamma5 closed form)       | yes |
| A2 | sympy         | 23   | `assert sp.simplify(K0_target_geom - 54 G c_s^5/(5 a^5 c^5)) == 0`| notes §3 (geometric target)         | yes |
| A3 | sympy         | 28   | `assert sp.simplify(Gamma5_target - 2 G/(5 c^5)) == 0`            | notes §3 (Gamma5_target)            | yes |
| A4 | sympy         | 37   | `assert sp.simplify(R0 - (NQ - 1)) == 0`                          | notes §4 (R_0)                      | partial (R0 is by-construction tautological after the substitution) |
| A5 | sympy         | 38   | `assert sp.simplify(R2 - (NQ - 1)) == 0`                          | notes §4 (R_2)                      | yes |
| A6 | sympy         | 39   | `assert sp.simplify(R4 - (NQ - 1)) == 0`                          | notes §4 (R_4)                      | yes |
| A7 | sympy         | 40   | `assert sp.simplify(R5 - (NQ - 1)) == 0`                          | notes §4 (R_5, fractional-power)    | yes |
| B1 | mathematica   | 41   | `expectZero[..., gamma5 - gamma5Expected]`                        | notes §2                            | yes |
| B2 | mathematica   | 49   | `expectZero["geometric target reduction", ...]`                   | notes §3                            | yes |
| B3 | mathematica   | 54   | `expectZero["Gamma5_target - 2G/(5c^5)", ...]`                    | notes §3                            | yes |
| B4 | mathematica   | 66   | `expectZero["R0 - (N_Q - 1)", r0 - (nQ - 1)]`                     | notes §4 (R_0)                      | partial (same tautological structure as A4) |
| B5 | mathematica   | 67   | `expectZero["R2 - (N_Q - 1)", ...]`                               | notes §4                            | yes |
| B6 | mathematica   | 68   | `expectZero["R4 - (N_Q - 1)", ...]`                               | notes §4                            | yes |
| B7 | mathematica   | 69   | `expectZero["R5 - (N_Q - 1)", ...]`                               | notes §4                            | yes |

Note on A4/B4: `R0` is built as `(N_Q * K0_target) / K0_target - 1`, which is algebraically `N_Q - 1` independent of any physics. Listing it is harmless (it documents the convention), but it is not load-bearing. The genuine content lives in A5–A7 / B5–B7, which test that `K_2`, `K_4`, and the fractional-power combination `Gamma_5` all scale linearly with `K_0` so that all four `R_i` collapse to the same scalar.

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_097.tex:21-25`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py:1-53`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl:1-75`

**What's wrong:**
The stage card's `\stagefield{Checks}` block lists three concrete checks (paper side, `stage_097.tex:21-25`):

> "Check the static limit \(\epsilon_2=\epsilon_4=0\) returns \(c_{\rm pole}=1/4\).
> Check \(l=0\) and \(l=2\) orthogonality before applying the geometry firewall.
> Check that any support/source success statement still carries the minimal-module hypothesis."

Neither the SymPy nor the Mathematica audit script touches `eps_2`, `eps_4`, a `c_pole = 1/4` static limit, `l=0`/`l=2` orthogonality, or the minimal-module hypothesis. The scripts only exercise the body-equation claim (`N_Q` defect identity and the four `R_i = N_Q - 1` identities). Notes §§1–5 likewise do not derive `c_pole`, the eps_i obstruction variables, or the orthogonality conditions — those references in the card appear to be carry-ins from Part III (the `\stagefield{Inputs}` field cites "the Part III minimal isotropic module" and "the static/dynamic geometry split").

Subtype: `script_missing_paper_claim` — the paper enumerates checks the script does not perform. Direction of resolution belongs to the user: either the card's `Checks` block is referencing upstream-verified machinery that this stage does not need to re-verify (in which case the card should be rewritten to say "carries forward" rather than "Check ..."), or those three checks really are part of this stage's verification scope and the scripts must add them.

**Why this matters:**
A reader who trusts the `Checks` block expects to find these three items exercised by the audit script. Right now they cannot be located in either engine. If the items are genuine verification obligations, the audit is incomplete; if they are upstream commitments, the card's wording mis-advertises the audit. Either way, the paper↔script contract is loose.

**Required change:**
See `## Resolve before fix_loop` block in the directive. Codex must not silently edit either side.

**Verification:**
Once the user picks a direction: (a) if scripts must add the checks, the verifier confirms three new `assert` / `expectZero` blocks appear and that the saved transcripts list passes for each; (b) if the card must be rewritten, the verifier confirms the `Checks` block in `stage_097.tex` matches the carry-forward language and the cross-check table in this report flips from `missing` to `carry-forward`.

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl:33-69`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage097_single_normalization_defect_sympy_audit.py:13-40`

**What's wrong:**
The `.wl` file is a line-by-line transliteration of the `.py` file. Compare the algebraic choreography:

Sympy (lines 13–17):
```
K2 = sp.simplify(K0 / (4 * Omega**2))
K4 = sp.simplify(K0 / (4 * Omega**4))
Gamma5 = sp.simplify(9 * K2**sp.Rational(5, 2) / K0**sp.Rational(3, 2))
Gamma5_expected = sp.simplify(9 * K0 / (32 * Omega**5))
assert sp.simplify(Gamma5 - Gamma5_expected) == 0
```

Mathematica (lines 33–41):
```
k2 = FullSimplify[k0/(4*omegaQ^2), ...];
k4 = FullSimplify[k0/(4*omegaQ^4), ...];
gamma5 = FullSimplify[9*k2^(5/2)/k0^(3/2), ...];
gamma5Expected = FullSimplify[9*k0/(32*omegaQ^5), ...];
expectZero["Gamma5 - 9 K0/(32 Omega^5)", gamma5 - gamma5Expected];
```

Same pattern repeats for `K0_target`, `K0_target_geom`, `Gamma5_target`, and the four `R_i` checks (sympy 32–40 vs mathematica 56–69): identical intermediate names, identical substitution order (`K0 -> nQ*K0_target`), identical residual structures `X_target - <literal>` and `R_i - (N_Q - 1)`. The Mathematica script never independently re-derives any quantity from a different algebraic route (e.g., it does not start from the conservative module `Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega^2)` and Taylor-expand to get `K2`, `K4`; it does not compute `Gamma5_target` by an independent power-of-omega expansion of the GR quadrupole formula). Variable-name differences (`K0` → `k0`, `Omega` → `omegaQ`, `G` → `gConst`, etc.) are cosmetic.

**Why this matters:**
The two-engine policy is in place precisely to catch algebra mistakes a single engine might miss (e.g., a sign on a fractional-power branch). A line-by-line port shares the human's algebraic decisions and therefore shares the failure modes; if SymPy's `simplify` mishandles `(K0/(4 Omega^2))^(5/2)` under positivity assumptions, Mathematica's `FullSimplify` performing exactly the same chain of substitutions will reproduce the same wrong answer. The independent-derivation guarantee is lost.

**Required change:**
Rework the Mathematica script so it derives the four identities along a different algebraic path. Concrete suggestion (Codex may adopt or substitute an equivalent independent route):

1. Start from the conservative module expression `YhatCons[w_] := 3/4 + (1/4)/(1 - w^2/omegaQ^2)`, define `KbarCons[w_] := k0 * YhatCons[w]`, and obtain `k2`, `k4` via `SeriesCoefficient[KbarCons[w], {w, 0, 2}]` and `SeriesCoefficient[KbarCons[w], {w, 0, 4}]`, then compare with `k0/(4 omegaQ^2)` and `k0/(4 omegaQ^4)`.
2. Derive `Gamma5` by independently writing the outgoing 2.5PN audit identity `Gammabar_5 = 9 (Kbar_2)^(5/2) / Kbar_0^(3/2)`, but substitute the *series-derived* `k2` from step 1 rather than reusing the `k0/(4 omegaQ^2)` literal.
3. For the geometric target, expand `k0Target /. omegaQ -> 3*cSound/(2*aRad)` by `PowerExpand` + `FullSimplify` (the existing path) but additionally cross-check `Gamma5_target` by first substituting `Omega -> 3 c_s/(2 a)` into the full `9 Kbar_0/(32 Omega^5)` and then simplifying, rather than reusing the `2G/(5 c^5)` literal as the target.
4. For the `R_i` reductions, instead of `subs K0 -> nQ k0Target` then dividing, build the ratios symbolically as `Kbar_2 / Kbar_2^target` etc., where the actual-branch `Kbar_i` are produced by the series route in step 1.

Codex may choose any independent route as long as the algebra is not a renamed copy of the SymPy script. The four bottom-line `expectZero` calls should remain (paper-side claims are unchanged).

**Verification:**
After Codex applies, the verifier (a) confirms the four `expectZero` calls in the `.wl` still pass; (b) inspects the new file diff and confirms the intermediate derivation does not align line-by-line with the SymPy script; (c) re-runs the saved transcript to confirm PASS.

## Independent-derivation check (Mathematica)

See finding F2. The `.wl` is structurally a transliteration of the `.py`: same intermediate variable names (modulo `K0→k0`, `Omega→omegaQ`, `G→gConst`, etc.), same algebraic choreography, same substitution order (`omegaQ -> omegaGeom`, `k0 -> nQ*k0Target`), same residual forms in each `expectZero`. No second algebraic route is taken (no series expansion of the conservative module, no independent derivation of `Gamma5_target` via the GR formula).

## Engine cross-check

Both engines produce the same final residuals: each prints `0` for the four identities and PASSes. SymPy transcript lines 9–18 and Mathematica transcript lines 13–35 list the same symbolic forms (`K2 = K0/(4 Omega**2)`, `Gamma5 = 9 K0/(32 Omega^5)`, `K0_target = 64 G Omega^5/(45 c^5)`, `R_i = N_Q - 1`). Engines agree on the numbers; the disagreement is only on the methodology (F2). No `engine_disagreement` finding.

## Verdict justification

The body-equation claim of the stage — that the actual isotropic passive/outgoing one-pole branch reduces to a single scalar defect `N_Q = Kbar_0/Kbar_0^target` with `R_0 = R_2 = R_4 = R_5 = N_Q - 1` — is faithfully and non-tautologically verified by both engines (A5–A7 / B5–B7 test the linear scaling of the three independent low-frequency coefficients and the fractional-power Gamma_5 combination; the geometric-target reduction A2/B2 is a real algebraic check that `(3 c_s/(2a))^5` propagates correctly). I attempted to break the assertions by: tracing whether `K2 = K0/(4 Omega^2)` is silently re-used in both sides of A1 (it is — but the load-bearing test is the `9 K2^(5/2)/K0^(3/2)` → `9 K0/(32 Omega^5)` algebraic identity, which requires the positivity assumption to handle the fractional power and is therefore non-trivial); checking that `R0` is tautological (it is, but A4 is a documentation artifact, not the load-bearing claim); and verifying that the `Omega → 3 c_s/(2a)` substitution gives `54/5` via `64/(45) * (3/2)^5 = 64*243/(45*32) = 15552/1440 = 10.8 = 54/5` (correct). The two findings are: (F1) the card's `\stagefield{Checks}` block enumerates three checks (static limit eps_2=eps_4=0 → c_pole=1/4; l=0/l=2 orthogonality; minimal-module hypothesis) that neither script touches, which is a `paper_misalignment` (script_missing_paper_claim) the user must resolve; and (F2) the Mathematica script is a line-by-line transliteration of the SymPy script, violating the second-engine policy. Neither finding warrants a stop-cold verdict; the body-equation claim is mathematically sound and downstream stages 100–106 depend on the `N_Q` definition, which is intact. No CRITICAL_DOWNSTREAM. Cosmetic note (not a finding): banner text and docstrings refer to "Stage 80" / "STAGE 080" while the file path / `\label` use `097`. This is a pre-existing renumbering artifact; not load-bearing.

## Self-test notes

- Variable independence: no `sp.diff`/`D[...]` is used; F2's suggested series-coefficient route uses `SeriesCoefficient[KbarCons[w], {w, 0, n}]` which is well-posed (KbarCons explicitly depends on `w`). No trivial-derivative trap.
- Symmetry/parity: no integration in this unit; n/a.
- Trivial-case pre-check: substituted `NQ = 1` mentally — all four `R_i` reduce to `0`, matching `N_Q - 1 = 0`. Substituted `Omega = 1`, `K0 = 1` into `9 K0/(32 Omega^5)` = `9/32` and into `9 (1/4)^(5/2)/1^(3/2) = 9/32` ✓.
- Paper round-trip: F1 routes to the user (no script edit prescribed), so it cannot introduce a new paper_misalignment. F2 prescribes a re-derivation that still verifies the same four paper-side identities (`Gamma5 closed form`, `geometric target reduction`, `Gamma5_target`, four `R_i`); no constant is changed.
- Path specifications: F2 edits the existing `.wl` file in place at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage097_single_normalization_defect_mathematica_audit.wl`; no new file needed.
