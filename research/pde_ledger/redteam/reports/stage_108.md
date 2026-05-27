---
unit_id: 108
batch: IV.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage108_robustness_classes.md
  paper_appendix: present
---

# Audit unit 108 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_108.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage108_robustness_classes.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only a single `\input{stages/stage_108}` line at L1250; no extra row content for this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage108_robustness_classes_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage108_robustness_classes_mathematica_audit.txt`

## What the paper claims

The stage classifies which isotropic DtN deformations of the canonical outgoing kernel `Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9` preserve the odd-normalization fingerprint `chi_Q = 1`. The notes enumerate three classes: (A) pure scale `Lambda^def = S Lambda^out` leaves `chi_Q = 1` exactly; (B) pure scale+argument `Lambda^def(z) = S Lambda^out(beta z)` only preserves canonical even moments at `beta = 1` (yielding `chi_Q = 1`); (C) additive throat-core channel with even-moment matching `Sigma_2 = -Sigma_0/9`, `Sigma_4 = -Sigma_0/27` gives `chi_Q = 3(S + 9 Sigma_5)/(3S - Sigma_0)`, with the special case `Sigma_5 = 0` giving `chi_Q = 3S/(3S - Sigma_0)`. The notes also boxed a general exact preservation submanifold `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27` covering arbitrary `beta`, not just `beta = 1`. The card additionally calls for verification checks against Robin loading and mixed-pole limits, but the notes restrict to scale/argument/additive classes only — so the Robin/mixed-pole language is only a card-side annotation about imported context, while the notes are the authoritative content for what to verify.

## What the script claims to verify

The SymPy script defines `Lambda_out = -3 + z^2/3 + z^4/9 + i z^5/9`, then exercises three classes. Class A: it verifies `Y_scale - Y_can = 0` (pure-scale invariance of the normalized outgoing series). Class B: it series-expands `Y_arg(beta z)`, prints `m2_arg = beta^2/9`, `m4_arg = 4 beta^4/81`, `chi_arg = beta^5`, solves for `beta in {-1, 1}`, and asserts `chi_arg(beta=1) - 1 = 0`. Class C (`beta = 1`): it solves the even-match constraints for `Sigma_2, Sigma_4`, derives `chi_add`, asserts `chi_add = 3(S + 9 Sigma_5)/(3S - Sigma_0)`, then solves `chi_add = 1` for `Sigma_5` and confirms the `beta=1` preservation locus `Sigma_5 = -Sigma_0/27`. The Mathematica script mirrors the same pipeline (with extra explicit residuals on `m2_arg`, `m4_arg`, `chi_arg`) but contains an operator-precedence defect in the `chi_arg(beta=1)` assertion.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Class A: `chi_Q = 1` under pure scale | sympy L32, wl L37 `expectZero["pure scale invariance", ...]` | match |
| Class B: `beta^2 = 1, beta^4 = 1`, `beta = 1`, `chi_Q = 1` | sympy L44-48; wl L48-55 | partial (see F2: wl substitution defect) |
| Class C even-moment match: `Sigma_2 = -Sigma_0/9, Sigma_4 = -Sigma_0/27` | sympy L58-63; wl L65-72 | match |
| Class C odd normalization: `chi_Q = 3(S+9 Sigma_5)/(3S-Sigma_0)` | sympy L66; wl L76 | match |
| Class C special case `Sigma_5 = 0` => `chi_Q = 3S/(3S - Sigma_0)` | (no assertion explicitly substitutes `Sigma_5 = 0`) | partial / missing (covered implicitly by general formula) |
| Exact preservation submanifold `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27` (general beta) | sympy L68-70; wl L78-81 only verify the `beta = 1` reduction `Sigma_5 = -Sigma_0/27` | partial (F1: missing general beta-dependent locus) |
| Card "Checks" item: Robin and standalone mixed-pole limits | (no script-side check) | missing (informational; notes do not require these — see verdict justification) |

`paper_alignment: partial` — the three classes and the general odd-normalization formula are correctly tested, but the general (β-dependent) exact preservation submanifold from the notes is not exercised, and the Mathematica β=1 check is malformed.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 32 | `expect_zero('pure scale invariance', Y_scale - Y_can)` | Class A: chi_Q=1 under scale | yes |
| A2 | sympy | 46-47 | `assert {sol[beta]} == {-1, 1}` | Class B: even-moment roots beta=+/-1 | yes |
| A3 | sympy | 48 | `expect_zero('chi_arg(beta=1) - 1', chi_arg.subs(beta, 1) - 1)` | Class B: chi_Q=1 at beta=1 | yes |
| A4 | sympy | 59-60 | `if len(sol) != 1: raise` | Class C: even-match uniqueness | yes |
| A5 | sympy | 66 | `expect_zero('chi_add - 3(S+9Sigma5)/(3S-Sigma0)', ...)` | Class C odd normalization | yes |
| A6 | sympy | 70 | `expect_zero('preservation locus check', chi_add.subs(Sigma5, chi_pres) - 1)` | beta=1 preservation special case | partial (only beta=1, not general) |
| B1 | wl | 37 | `expectZero["pure scale invariance", yScale - yCan]` | Class A | yes |
| B2 | wl | 48 | `expectZero["m2_arg - beta^2/9", m2Arg - beta^2/9]` | Class B m2 | yes |
| B3 | wl | 49 | `expectZero["m4_arg - 4 beta^4/81", ...]` | Class B m4 | yes |
| B4 | wl | 50 | `expectZero["chi_arg - beta^5", ...]` | Class B chi | yes |
| B5 | wl | 54 | `If[betaRoots =!= {-1, 1}, fail[...]]` | Class B even-moment roots | yes |
| B6 | wl | 55 | `expectZero["chi_arg(beta=1) - 1", chiArg /. beta -> 1 - 1]` | Class B chi(beta=1)=1 | NO (operator-precedence bug; see F2) |
| B7 | wl | 66 | `If[Length[sol] != 1, fail[...]]` | Class C even-match uniqueness | yes |
| B8 | wl | 71-72 | `expectZero["Sigma2... + sigma0/9"], ["Sigma4... + sigma0/27"]` | Class C even-moment values | yes |
| B9 | wl | 76 | `expectZero["chi_add - 3(...)..."]` | Class C odd normalization | yes |
| B10 | wl | 80-81 | `expectZero["Sigma5 preservation locus + sigma0/27"], ["preservation locus check"]` | beta=1 preservation special case | partial (only beta=1, not general) |

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage108_robustness_classes.md:90-92`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:68-70`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:78-81`

**Subtype:** script_missing_paper_claim

**What's wrong:**
The notes box a *general* exact preservation submanifold valid for arbitrary `beta`:

> "The exact condition is `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27`. This includes as special cases: pure scale deformation, pure scale+argument deformation with beta = 1, additive core deformations whose odd slot is locked to the static shift."

The scripts, however, fix `beta = 1` for the Class C analysis (sympy L51: `Lambda_add = sp.expand(S*Lambda_out + ...)` with no `beta z` substitution; wl L57 similarly). Consequently the preservation-locus check (sympy L68-70, wl L78-81) only verifies `Sigma_5 = -Sigma_0/27`, which is the `beta = 1` reduction of the boxed condition. The general `beta`-dependent identity `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27` — explicitly enumerated in the notes as a unifier of all three classes — is never tested.

**Why this matters:**
The notes treat the general preservation submanifold as the load-bearing closure of the classification ("This includes as special cases ..."). Without a script-side check for the `beta`-dependent form, the unification claim that A, B (beta=1), and C all live on the same submanifold is unverified.

**Required change:**
In both scripts, after the existing `beta = 1` Class C block, add a generalized Class C+beta block that combines scale+argument with additive core, computes `chi_Q(S, beta, Sigma_0, Sigma_5)` (with `Sigma_2, Sigma_4` re-solved against the even-moment match in this generalized setup), and verifies that `chi_Q = 1` iff `Sigma_5 = S*(1 - beta^5)/9 - Sigma_0/27`.

**Verification:**
After Codex applies, the verifier should see a new assertion in each script verifying the general preservation submanifold formula. The existing `beta = 1` check should still pass as a sanity reduction.

### F2 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:55`

**What's wrong:**
The line reads:

```
expectZero["chi_arg(beta=1) - 1", chiArg /. beta -> 1 - 1];
```

In Mathematica, `Plus`/`Subtract` bind tighter than `Rule` and `ReplaceAll`. So this parses as `chiArg /. (beta -> (1 - 1))` = `chiArg /. (beta -> 0)`. Since `chiArg = beta^5`, substituting `beta -> 0` gives `0^5 = 0`, which the `expectZero` happily accepts. The assertion is therefore NOT checking "chi at beta=1 minus 1 equals 0"; it is checking "chi at beta=0 equals 0", which holds trivially for `beta^5` and tells us nothing about the beta=1 deformation closing back to canonical. The SymPy counterpart at L48 uses `chi_arg.subs(beta, 1) - 1` and is correct.

This is the load-bearing Class-B `chi_Q = 1` check in Mathematica; right now it would still "pass" if `chiArg` were `beta^5 + beta^3` or `beta^7` or any expression that vanishes at the origin — including expressions that do NOT equal 1 at beta=1.

**Why this matters:**
A second-engine check that passes by precedence accident rather than by computing the intended substitution is not a check. It is also exactly the kind of silent-pass bug the audit framework exists to catch.

**Required change:**
Edit `mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:55` from:

```
expectZero["chi_arg(beta=1) - 1", chiArg /. beta -> 1 - 1];
```

to:

```
expectZero["chi_arg(beta=1) - 1", (chiArg /. beta -> 1) - 1];
```

**Verification:**
After the fix, the printed line `chi_arg(beta=1) - 1 = 0` should still appear (it remains true because `1^5 - 1 = 0`), but the residual now genuinely substitutes beta=1 before subtracting 1. To make the test discriminating, the verifier can also check (mentally) that perturbing `chiArg = beta^5` to `beta^5 + (beta - 1)` would now fail the assertion (residual becomes 0 still — bad test); perturbing to `2*beta^5` would now correctly fail (residual = 1). The fix is the local one-line parenthesization; the broader semantic check is enforced by F1's added general-locus assertion.

### F3 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage108_robustness_classes_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage108_robustness_classes_mathematica_audit.txt`

**What's wrong:**
The SymPy and Mathematica banners in both scripts read `STAGE 91` (sympy L25) and `STAGE 091` (wl L26), but this unit is stage 108 (label `stage:108`, section heading "Stage 125"). The saved transcripts echo "STAGE 91 / 091". Output mtimes are newer than script mtimes, so the captured output is fresh; the banner mismatch is in the *script source* (`banner('STAGE 91 — ROBUSTNESS CLASSES FOR chi_Q')` and `banner["STAGE 091 — ..."]`).

This is not strictly "stale_output" (mtimes confirm freshness), but the printed banner labelling the wrong stage number is a maintenance hazard if outputs are grepped by stage. Flagged at low severity for traceability cleanup.

**Why this matters:**
Stage banners are used by the orchestrator and humans to associate transcripts with the correct unit. A mismatched banner can mislead someone grepping output by stage number into thinking this unit's transcript belongs to stage 91.

**Required change:**
Edit `scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:25` from `banner('STAGE 91 — ROBUSTNESS CLASSES FOR chi_Q')` to `banner('STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q')`. Edit `mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:26` from `banner["STAGE 091 — ROBUSTNESS CLASSES FOR chi_Q"];` to `banner["STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q"];`. (Stage 108 is the canonical unit number used by the audit framework; the card heading "Stage 125" is the published label, but every other artifact in the framework refers to this as unit 108.)

**Verification:**
After re-running, the transcript header should read `STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q` in both engines.

## Independent-derivation check (Mathematica)

The Mathematica script and the SymPy script share the same overall pipeline (define `lambdaOut`/`Lambda_out`, build `Y_scale`, `Y_arg`, do even-moment match, then Class-C additive core). Variable names are independent (`sNorm` vs `S`, lowercase `sigma{0,2,4,5}` vs `Sigma{0,2,4,5}`); the Mathematica script adds explicit intermediate residuals (`m2_arg - beta^2/9`, `m4_arg - 4 beta^4/81`, `chi_arg - beta^5`, `Sigma2 + sigma0/9`, `Sigma4 + sigma0/27`, `Sigma5 + sigma0/27`) that the SymPy script does not check separately. This is more thorough than the SymPy side, not less, so it is not a transliteration. Both engines independently apply the series expansion / coefficient extraction approach mandated by the paper's classification. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce:
- `pure scale invariance = 0`
- `m2_arg = beta^2/9`, `m4_arg = 4 beta^4/81`, `chi_arg = beta^5`
- `beta` solutions = `{-1, 1}`
- `Sigma2 = -Sigma_0/9`, `Sigma4 = -Sigma_0/27`
- `chi_add = 3(S + 9 Sigma_5)/(3S - Sigma_0)` (Mathematica prints the factored form `(-3*(9*sigma5 + sNorm))/(sigma0 - 3*sNorm)` which is identical up to sign of numerator and denominator)
- `Sigma_5` preservation locus at `beta = 1` = `-Sigma_0/27`

Engines agree. The only engine-side defect is F2 (Mathematica chi_arg(beta=1) check parses wrong; passes for the wrong reason); SymPy is correct.

## Verdict justification

`verdict: findings`. The paper's three-class structure (A scale, B scale+argument, C additive core) and the Class-C odd normalization formula `chi_Q = 3(S + 9 Sigma_5)/(3S - Sigma_0)` are correctly verified by both engines. Two substantive defects: (F1) the general `beta`-dependent exact preservation submanifold `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27` from the notes is not exercised — only the `beta=1` reduction is; (F2) the Mathematica `chi_arg(beta=1) - 1` assertion has an operator-precedence bug that makes it substitute `beta=0` instead of `beta=1`, so it passes by accident. F3 is a cosmetic banner mismatch. No `stop_cold` is warranted: F1 is a missing-check rather than a wrong-result, F2 is an isolated typo fixed by parenthesization, and F3 is informational. F1 is `paper_misalignment` (subtype `script_missing_paper_claim`) — the orchestrator will halt for user resolution since the unit's notes assert a more general result than the scripts verify; the user must decide whether the scripts should be extended to the general locus or the notes should be downgraded to the `beta=1` special case.

Attacks tried that failed: (a) checking whether `chi_add`'s `Sigma_2, Sigma_4` substitution might double-count even-moment constraints — no, the `sol = sp.solve(...)` is on `(Sigma_2, Sigma_4)` only, and `chi_add` is then evaluated at that unique solution; sound. (b) checking whether the series truncation order 5 is sufficient — yes, the odd normalization lives in the `z^5` coefficient which is fully captured by `series(..., 6).removeO()`. (c) checking `nonzero=True` assumption on `S`, `beta` — appropriate since `S` divides into `Y_scale` and `beta` appears in series of `Lambda_out(beta z)`; would be ill-defined at 0. (d) attempting to find a tautology in `chi_add - 3(S+9 Sigma_5)/(3S-Sigma_0)` — `chi_add` is built from `(-L5/L0).subs(sol) / (1/27)` where `L0 = S*Lambda_out|_{z=0} + Sigma_0 = -3S + Sigma_0` and `L5 = Sigma_5 / I * I = Sigma_5` (after `im()`); so `chi_add = 27 * (-Sigma_5 / (-3S + Sigma_0)) = 27 Sigma_5 / (3S - Sigma_0)`... wait — this only equals `3(S + 9 Sigma_5)/(3S - Sigma_0)` if `S` also enters `L5`. Re-reading: `Lambda_add = S*Lambda_out + Sigma_0 + Sigma_2 z^2 + Sigma_4 z^4 + I Sigma_5 z^5`. The `z^5` coefficient is `S*(1/9) + I Sigma_5` (from `Lambda_out`'s `I z^5/9` term plus the additive `I Sigma_5 z^5`). Im of that is `S/9 + Sigma_5`. So `L5 = S/9 + Sigma_5`. Then `-L5/L0 = -(S/9 + Sigma_5)/(-3S + Sigma_0) = (S/9 + Sigma_5)/(3S - Sigma_0)`. Divided by `1/27` gives `27(S/9 + Sigma_5)/(3S - Sigma_0) = (3S + 27 Sigma_5)/(3S - Sigma_0) = 3(S + 9 Sigma_5)/(3S - Sigma_0)`. Matches the paper. Sound. The assertion is non-tautological — it would fail if e.g. the `I z^5/9` term in `Lambda_out` were wrong.

## Self-test notes

Checked: (1) F2's Mathematica precedence — verified by checking the Wolfram precedence table (Plus > Rule > ReplaceAll); the unparenthesized form does substitute `beta -> 0`. (2) F1's required change is non-trivial; the directive must instruct Codex to combine `Lambda_out(beta z)` with the additive core, re-solve `Sigma_2, Sigma_4` (which now depend on `beta`), then verify the `Sigma_5 = S(1-beta^5)/9 - Sigma_0/27` locus. I have flagged this in the directive as needing user resolution (since it overlaps with `paper_misalignment`), so Codex will not auto-apply F1 — only F2 and F3 are mechanical. (3) Confirmed that fixing F2's parenthesization does not introduce any new misalignment with paper (`chi_arg(beta=1) = 1` is the correct Class B closure per the notes).
