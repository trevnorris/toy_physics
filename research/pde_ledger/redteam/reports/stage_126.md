---
unit_id: 126
batch: IV.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - notes/stages/moving_throat_pde_stage126_positive_source_families.md
  paper_appendix: present
---

# Audit unit 126 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_126.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage126_positive_source_families.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only a single `\input{stages/stage_126}` at L1286; no extra row content)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage126_positive_source_families_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage126_positive_source_families_mathematica_audit.txt`

## What the paper claims

The stage card (line 16) states: "Self-matched derivative profile gives `g = pi/4`; a positive convex family hits the lower branch." The notes (authoritative since the .tex body is one sentence) give three concrete deliverables:

1. The self-matched normalized source `sigma_match(z) = k cos(kz)` with `k = pi/(2L)` integrates to 1 over `[0,L]` and yields the exact mouth-bias `g_match = pi/4`.
2. The convex positive family `sigma_xi(z) = (1 - xi) k cos(kz) + xi/L` for `0 <= xi <= 1` is normalized and nonnegative on `[0,L]`, with bias `g_xi = (1 - xi)(pi/4) + xi (2/pi)`.
3. Because `2/pi < g_-^F1 < pi/4`, there is a unique `xi_* in (0,1)` such that `g_{xi_*} = g_-^F1`, with closed form `xi_* = (pi/4 - g_-^F1)/(pi/4 - 2/pi)`. The notes also give the explicit surd form (containing a transcription typo discussed in F1) and the numerical value `xi_* ~= 0.183918405511538`.

The card's `Checks` field additionally requests: (i) positivity of the mouth source, (ii) zero-flux and boundary-layer normalizations in the GNLS/localized-Maxwell reduction, (iii) Family-1 compensation against the lower (not equal-normalized) branch. `Inputs` declares that `g_-^F1` is imported from "the lower compensated core branch" (upstream Family-1 branch derivation), so the script may treat its symbolic form as given.

## What the script claims to verify

Both the SymPy and Mathematica scripts implement the same five-step verification:

1. Compute `int_0^L k cos(kz) dz` and assert it equals 1 (`self-matched normalization`).
2. Compute `int_0^L k cos(kz) * cos(kz) dz` and assert it equals `pi/4` (`self-matched bias`).
3. Define `g_-^F1 = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3))/(20 pi)` as a symbolic literal carried in from the Family-1 branch; print numeric values of `delta_g_match`, `T_-/T_match`.
4. Compute `int_0^L sigma_xi cos(kz) dz` for `sigma_xi = (1-xi) k cos(kz) + xi/L`; assert convex-family normalization equals 1; solve `g_xi = g_-^F1` for `xi`; assert `g_xi(xi_*) - g_-^F1 = 0`.
5. Verify numerically that `2/pi < g_-^F1 < pi/4` (Mathematica calls `fail[]` on violation; SymPy only prints `True`/`False`).

The scripts do NOT explicitly verify pointwise nonnegativity of `sigma_match` or `sigma_xi`. Both scripts print a banner that mislabels the stage as "STAGE 109" rather than 126.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `sigma_match` integrates to 1 on `[0,L]` | sympy L17,22; wl L33,38 `expectZero["self-matched normalization", normMatch - 1]` | match |
| `g_match = pi/4` (self-matched bias) | sympy L18,22; wl L34,39 `expectZero["self-matched bias", gMatch - Pi/4]` | match |
| Imported `g_-^F1 = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3))/(20 pi)` | sympy L25-26; wl L41-42 (literal symbolic form, no derivation) | match (Inputs declares this as imported; informational, see F4) |
| Convex family normalized to 1 | sympy L40,45; wl L54,59 | match |
| Convex family pointwise nonnegative on `[0,L]` for `0 <= xi <= 1` | (no script-side check) | missing (see F2) |
| Unique `xi_* in (0,1)` with `g_{xi_*} = g_-^F1` | sympy L48 `sp.solve(...)`, L52-55 assert; wl L61 `Solve[..., Reals]`, L64 assert | match (existence + value via Solve; uniqueness implicit via linear `g_xi(xi)`) |
| `xi_* ~= 0.183918405511538` numeric | sympy/wl print `xi_* numeric` matching to 20 digits | match |
| Interval `2/pi < g_-^F1 < pi/4` | sympy L60 prints bool; wl L68-70 asserts via `fail[]` | match (sympy side is print-only, not asserted; minor weakness — see verdict) |
| Card `Checks` (i) positivity of mouth source | (no explicit assertion) | missing (F2) |
| Card `Checks` (ii) zero-flux / boundary-layer normalizations | the `int = 1` integrations partially cover normalization; zero-flux not addressed | partial (informational; boundary-layer reduction is per Inputs from upstream) |
| Card `Checks` (iii) Family-1 lower branch (not equal-normalized) | scripts hit the explicit `g_-^F1` formula whose numerical value is `~0.758` (lower branch) | match |
| Notes' explicit surd for `xi_*` containing `4107 - 168 pi^2` | scripts produce `4107 - 100 pi^2` (script value is correct; notes have a typo) | mismatch (F1, paper_misalignment) |

`paper_alignment: partial` — the four primary identities (normalizations, `g_match = pi/4`, convex-family compensation) are all faithfully verified, but the notes' written surd form of `xi_*` contradicts the script (F1) and two `Checks` items (positivity, zero-flux/boundary-layer) are not explicitly asserted (F2).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 22 | `if simplify(norm_match - 1) != 0 or simplify(g_match - pi/4) != 0: raise` | self-matched normalization + bias | yes |
| A2 | sympy | 45-46 | `if simplify(norm_xi - 1) != 0: raise` | convex-family normalization | yes |
| A3 | sympy | 54-55 | `if simplify(g_xi_star) != 0: raise` (where `g_xi_star = g_xi(xi_*) - g_-^F1`) | Family-1 compensation `g_{xi_*} = g_-^F1` | yes |
| A4 | sympy | 60 | `print(bool(2/pi < gminus < pi/4))` | interval pre-condition | partial (printed, not asserted) |
| B1 | wl | 38 | `expectZero["self-matched normalization", normMatch - 1]` | self-matched normalization | yes |
| B2 | wl | 39 | `expectZero["self-matched bias", gMatch - Pi/4]` | self-matched bias | yes |
| B3 | wl | 59 | `expectZero["convex-family normalization", normXi - 1]` | convex-family normalization | yes |
| B4 | wl | 64 | `expectZero["g_xi(xi_*) - g_-", (gXi /. xi -> xiStar) - gMinus]` | Family-1 compensation | yes |
| B5 | wl | 68-70 | `If[!TrueQ[intervalCheck], fail[...]]` | interval pre-condition | yes |

No `sigma >= 0` assertion in either engine. No explicit assertion that `xi_*` lies in `(0,1)` (only printed numerically). The Family-1 surd `(2 sqrt(4107 - 100 pi^2) - 37 sqrt(3))/(20 pi)` is hand-typed in both scripts — informational issue (F4).

## Findings

### F1 — paper_misalignment

**Severity:** medium
**Subtype:** notes_contradicts_script
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage126_positive_source_families.md:100`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:25,48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:41,61`

**What's wrong:**
The notes give the boxed exact surd form for the convex-family compensation point as

> `xi_* = (-37*sqrt(3) - 5*pi^2 + 2*sqrt(4107 - 168*pi^2)) / (5*(8 - pi^2)) ~= 0.183918405511538`

(notes line 100, inside the box). The scripts compute (via `sp.solve`/`Solve`)

```
xi_* = (-37*sqrt(3) - 5*pi^2 + 2*sqrt(4107 - 100*pi^2)) / (5*(8 - pi^2))
```

(sympy L48 output L25; wl L61 output L29 in normalized form).

The two forms differ only in `4107 - 168*pi^2` (notes) vs `4107 - 100*pi^2` (scripts). Numerically:

- Script form: `2*sqrt(4107 - 100*pi^2) ~= 2*sqrt(3120.04) ~= 111.715`, so numerator `~= -64.086 - 49.348 + 111.715 = -1.719`, divided by `5*(8 - 9.870) ~= -9.348` gives `xi_* ~= 0.18392`. Matches the script printout `0.18391840551153962831` and the notes' stated numerical value `0.183918405511538`.
- Notes form: `2*sqrt(4107 - 168*pi^2) ~= 2*sqrt(2448.5) ~= 98.965`. Numerator `~= -64.086 - 49.348 + 98.965 = -14.469`, divided by `-9.348` gives `xi_* ~= 1.548`. This is outside `(0,1)` and does NOT equal the notes' own quoted numerical value.

So the notes contain a transcription typo: `168` should be `100`. The notes' numerical answer and the scripts agree; only the algebraic display in the notes is wrong.

Sanity check by re-deriving from the closed form `xi_* = (pi/4 - g_-^F1)/(pi/4 - 2/pi)`:
- Numerator with `g_-^F1 = (2 sqrt(4107 - 100 pi^2) - 37 sqrt(3))/(20 pi)`:
  `pi/4 - g_-^F1 = (5 pi^2 - 2 sqrt(4107 - 100 pi^2) + 37 sqrt(3))/(20 pi)`.
- Denominator: `pi/4 - 2/pi = (pi^2 - 8)/(4 pi) = (pi^2 - 8)/(4 pi)`.
- Ratio: `xi_* = ((5 pi^2 - 2 sqrt(4107 - 100 pi^2) + 37 sqrt(3))/(20 pi)) * (4 pi/(pi^2 - 8)) = (5 pi^2 - 2 sqrt(4107 - 100 pi^2) + 37 sqrt(3))/(5(pi^2 - 8)) = -(5 pi^2 - 2 sqrt(4107 - 100 pi^2) + 37 sqrt(3))/(5(8 - pi^2)) = (2 sqrt(4107 - 100 pi^2) - 5 pi^2 - 37 sqrt(3))/(5(8 - pi^2))`.

This matches the script form exactly with `100*pi^2`, confirming the notes typo.

**Why this matters:**
A reader cross-checking the script against the notes will see a symbolic mismatch in a load-bearing boxed equation. The numerical value happens to agree, but the algebraic display in the notes is wrong, and silently relying on numerical agreement obscures the typo. Future re-derivations or re-formulations that take the notes' symbolic form as truth will produce wrong results.

**Required change:**
**This is a paper_misalignment requiring user resolution.** See the directive's `## Resolve before fix_loop` block. Codex must not auto-edit the notes or the scripts.

**Verification:**
After user resolution: the notes' boxed `xi_*` formula and the scripts' computed `xi_*` symbolic form should be algebraically identical (up to sign normalization), and the numerical value `0.183918405511538...` should remain unchanged.

### F2 — script_missing_paper_claim

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_126.tex:22` (first item of `\stagefield{Checks}`)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage126_positive_source_families.md:78` ("This is normalized and nonnegative on [0,L]")
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:39` (definition of `sigma_xi`; no positivity assertion follows)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:53` (definition of `sigmaXi`; no positivity assertion follows)

**What's wrong:**
The paper card's `Checks` list begins with "Check positivity of the mouth source; do not use signed fitting to reach the upper branch." The notes also boxes "This is normalized and nonnegative on `[0,L]`." The scripts define `sigma_match = k cos(kz)` and `sigma_xi = (1 - xi) k cos(kz) + xi/L`, then check normalization and the Family-1 compensation, but never assert nonnegativity of either source on `z in [0, L]` for `xi in [0, 1]`.

It is mathematically obvious that on `z in [0, L]` we have `kz in [0, pi/2]` and `cos(kz) >= 0`, so `sigma_match >= 0` and `sigma_xi >= 0` for `xi in [0, 1]`. But "obvious" is not a script-side check, and the paper explicitly enumerates positivity as a required check (it is the distinguishing feature versus signed fitting to the upper branch).

**Why this matters:**
The whole purpose of Stage 126 versus other compensation-point analyses is reaching the lower branch with a *positive* source — the card explicitly forbids "signed fitting to reach the upper branch." Verifying that the convex family really is nonnegative for the full parameter range `xi in [0,1]` is the load-bearing positivity gate. Without an assertion, the scripts cannot detect a future regression where someone changes the family to one that goes negative.

**Required change:**
Add to both scripts a nonnegativity check for `sigma_xi(z)` over `z in [0, L]` and `xi in [0, 1]`. Since `sigma_xi = (1-xi) k cos(kz) + xi/L` and both terms are individually nonnegative on the stated domain, a clean assertion is: verify that the minimum of `sigma_xi` over the box `(z, xi) in [0,L] x [0,1]` is `>= 0`. Concretely, one can assert each summand is nonnegative:

- `(1 - xi) >= 0` for `xi in [0, 1]` (trivial bound).
- `k cos(kz) >= 0` for `z in [0, L]` (since `kz in [0, pi/2]`).
- `xi/L >= 0` for `xi in [0, 1]` and `L > 0`.

In SymPy: e.g.,
```python
# Positivity: sigma_xi >= 0 on z in [0, L], xi in [0, 1]
assert sp.simplify(sp.cos(k*z).rewrite(sp.cos).subs(z, 0) - 1) == 0  # boundary check
assert bool(sp.cos(k*sp.Symbol('z0', positive=True)).subs('z0', L).equals(0))  # cos(pi/2) = 0
# Better: directly verify endpoints and use sp.solveset for interior minimum
min_match = sp.minimum(k*sp.cos(k*z), z, sp.Interval(0, L))
assert sp.simplify(min_match) == 0  # exact: cos(pi/2) = 0
# At xi extremes:
assert sp.minimum(sigma_xi.subs(xi, 0), z, sp.Interval(0, L)) == 0
assert sp.simplify(sigma_xi.subs(xi, 1) - 1/L) == 0  # constant 1/L > 0
```

(See directive F2 for the precise patch.)

In Mathematica: use `Minimize[sigmaXi, {z, xi} \[Element] Rectangle[{0, 0}, {lM, 1}]]` and assert the minimum is `>= 0`.

**Verification:**
After Codex applies, the verifier should see a new `expectZero` (or assert) for positivity in both engines, and the saved transcript should print the verified minimum as `>= 0`.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:60`

**What's wrong:**
The SymPy script prints

```python
print("Check 2/pi < g_- < pi/4 ->", bool(sp.N(2/sp.pi) < sp.N(gminus) < sp.N(sp.pi/4)))
```

but does NOT raise on `False`. The Mathematica counterpart (wl L68-70) does correctly halt with `fail[...]` if the interval check fails. The notes treat the inequality `2/pi < g_-^F1 < pi/4` as the load-bearing existence guarantee for `xi_* in (0,1)`; if this inequality were ever to fail (e.g. due to a future revision of the imported `g_-^F1`), the SymPy script would still exit 0 and print a banner-passing transcript.

**Why this matters:**
The script's existence claim for `xi_* in (0,1)` depends on the bracketing inequality. A silent `False` print here would not halt the audit, defeating the purpose of the check.

**Required change:**
Convert the print to a raising assertion in the SymPy script. E.g.,
```python
interval_check = bool(sp.N(2/sp.pi) < sp.N(gminus) < sp.N(sp.pi/4))
print("Check 2/pi < g_- < pi/4 ->", interval_check)
if not interval_check:
    raise AssertionError("g_-^F1 does not lie strictly between 2/pi and pi/4.")
```

**Verification:**
After Codex applies, the SymPy script should still print `True` (output unchanged) but a hypothetical regression would now raise. The Mathematica side already does this.

### F4 — stale_output

**Severity:** low (informational, banner mislabel)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:13`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:26`
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage126_positive_source_families_sympy_audit.txt:11`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage126_positive_source_families_mathematica_audit.txt:11`

**What's wrong:**
Both scripts emit `banner("STAGE 109 — EXPLICIT POSITIVE SOURCE FAMILIES")` (sympy L13, wl L26), and the saved transcripts (output L11 in both) echo this. The unit is stage 126, not 109. Output mtimes (sympy 2026-05-11 12:45, wl 2026-05-11 13:10) are both newer than the script mtimes (sympy Apr 1 21:36, wl May 11 11:56), so the transcripts are fresh — the misnomer is in the script source.

This is not strictly `stale_output` (mtimes are fine), but the banner text labels the wrong stage in the captured transcript, which is a traceability hazard if outputs are grepped by stage number. Flagged at low severity per the audit framework's standing convention.

**Why this matters:**
A grep for `STAGE 126` in `scripts/output/` and `mathematica/output/` will miss this unit's transcripts. A grep for `STAGE 109` will return both the real stage 109 transcripts and this unit's mislabeled transcripts, causing cross-stage confusion.

**Required change:**
- In `scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:13`, change `banner("STAGE 109 — EXPLICIT POSITIVE SOURCE FAMILIES")` to `banner("STAGE 126 — EXPLICIT POSITIVE SOURCE FAMILIES")`.
- In `mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:26`, change `banner["STAGE 109 — EXPLICIT POSITIVE SOURCE FAMILIES"];` to `banner["STAGE 126 — EXPLICIT POSITIVE SOURCE FAMILIES"];`.

**Verification:**
After re-running, the saved transcripts should print `STAGE 126 — EXPLICIT POSITIVE SOURCE FAMILIES` instead of `STAGE 109 — ...`.

## Independent-derivation check (Mathematica)

The Mathematica and SymPy scripts share a near-identical pipeline:

- Same `k = pi/(2 L)`.
- Same `sigma_match` and `sigma_xi` definitions.
- Same integration calls (`Integrate[sigmaMatch, {z, 0, lM}]` <-> `sp.integrate(sigma_match, (z, 0, L))`).
- Same imported `g_-^F1` symbolic literal `(2*sqrt(4107 - 100*pi^2) - 37*sqrt(3))/(20*pi)`.
- Same `Solve[gXi == gMinus, xi, Reals]` <-> `sp.solve(sp.Eq(g_xi, gminus), xi)[0]`.

Variable-name divergences are cosmetic (`L` vs `lM`, `R` vs `rDisc`, `gminus` vs `gMinus`). The Mathematica script does add an explicit `If[!TrueQ[intervalCheck], fail[...]]` halt that SymPy omits (see F3).

The underlying computation is just a definite integral and a `Solve`. The "independent derivation" bar is hard to clear when the claim is structurally that simple — both engines necessarily must compute the same integral. The mirroring is not gratuitous transliteration; both engines genuinely re-derive `g_match = pi/4` from the integral definition. I do not flag `mathematica_transliteration` here, but note that the .wl is structurally a tight mirror of the .py.

## Engine cross-check

Both engines produce (after normalization):

- `sigma_match normalization = 1` (both)
- `g_match = pi/4` (both)
- `g_-^F1 = (-37*sqrt(3) + 2*sqrt(4107 - 100*pi^2))/(20*pi)` (both; same form modulo Mathematica's `Sqrt[]` vs SymPy's `sqrt`)
- Numerical `g_match = 0.78539816339744830962` (SymPy) / `0.78539816339744830961566084581987572105` (Mathematica) — agree to all printed digits
- Numerical `g_-^F1 = 0.75803507894466282692` (SymPy) / `0.75803507894466282691968089041411045778` (Mathematica) — agree
- Numerical `xi_* = 0.18391840551153962831` (SymPy) / `0.18391840551153962830834428246014671147` (Mathematica) — agree
- `g_xi(xi_*) - g_- = 0` (both `PASS`)
- Interval check `True` (both)

The Mathematica simplification gives `g_xi` in the form `-1/4*(Pi*(-1 + xi)) + (2*xi)/Pi` and `xi_* = 1 + (40 + 37*Sqrt[3] - 2*Sqrt[4107 - 100*Pi^2])/(5*(-8 + Pi^2))`, which are algebraically identical to SymPy's `(8*xi + pi**2*(1 - xi))/(4*pi)` and `(-37*sqrt(3) - 5*pi^2 + 2*sqrt(4107 - 100*pi^2))/(5*(8 - pi^2))` respectively (verified by sign-tracking).

Engines agree. `engines_agree: true`.

## Verdict justification

`verdict: findings`. The four primary identities the stage exists to prove — `sigma_match` normalization, `g_match = pi/4`, convex-family normalization, and Family-1 compensation `g_{xi_*} = g_-^F1` — are correctly verified in both engines, and the engines agree numerically and symbolically. The findings are:

- F1 (`paper_misalignment`, `notes_contradicts_script`): the notes' boxed surd form for `xi_*` contains a transcription typo (`168*pi^2` should be `100*pi^2`); the script's symbolic form and the notes' own numerical value agree. Routed to user via `## Resolve before fix_loop` (almost certainly: fix the notes typo).
- F2 (`script_missing_paper_claim`): the paper card's first `Checks` item — "Check positivity of the mouth source" — is not explicitly asserted in either engine, even though mathematical nonnegativity is straightforward. The whole point of this stage is reaching the lower branch with a positive source, so an explicit positivity assertion is required.
- F3 (`insufficient_verification`): SymPy interval check prints but does not raise; Mathematica side correctly raises. Small but real hardening.
- F4 (banner mislabel): both scripts print "STAGE 109" instead of "STAGE 126" — traceability hazard, mechanical fix.

No `stop_cold`. F1 is `paper_misalignment` (notes typo, user resolves direction). F2/F3/F4 are mechanical and do not propagate downstream: positivity is structurally guaranteed by the definitions, the interval check passes numerically, and the banner is cosmetic. The card's Inputs declare `g_-^F1` as an imported result (from the lower-compensated-core-branch upstream Family-1 derivation), so `hardcoded_result` does not apply — the script is honoring the card's declared import boundary.

Attacks tried that failed: (a) checking whether `2/pi < g_-^F1 < pi/4` actually holds — numerically `0.63662 < 0.75804 < 0.78540`, yes. (b) checking whether `g_xi` is monotone in `xi` so `xi_*` is unique — yes, `g_xi = (1-xi)(pi/4) + xi(2/pi)` is affine in `xi`, slope `(2/pi - pi/4) ~= -0.149 < 0`, monotone decreasing, unique solution in `(0,1)` given the bracketing. (c) checking whether `sp.solve(sp.Eq(g_xi, gminus), xi)` could return spurious branches — `g_xi` is linear in `xi`, exactly one solution. (d) checking whether `simplify` under `xi: real=True` (no positivity) could hide a sign issue — no, `xi_*` is computed by direct linear solve, not branch-sensitive. (e) checking whether the assertion `simplify(g_xi_star) != 0` is tautological — `g_xi_star = g_xi(xi_*) - g_-^F1`, where `g_xi` is computed by integration and `xi_*` is solved against `g_-^F1`; the assertion non-tautologically tests that the closed-form `Solve` solution actually substitutes back to zero (would fail if `simplify` couldn't reduce the surds, which is a real risk for nested radicals). (f) checking whether the imported `g_-^F1` literal in the script matches the upstream Family-1 derivation — the symbolic form `(2 sqrt(4107 - 100 pi^2) - 37 sqrt(3))/(20 pi)` produces `~0.75804`, matching the notes' stated `g_-^F1 ~= 0.758035078944663`; concordant.

## Self-test notes

Checked: (1) F1's numerical re-derivation from `xi_* = (pi/4 - g_-^F1)/(pi/4 - 2/pi)` — confirmed the script's `100*pi^2` form is correct and the notes' `168*pi^2` is the typo. (2) F2's positivity claim — verified `cos(kz) >= 0` on `kz in [0, pi/2]` (i.e., `z in [0, L]`) and that `sigma_xi` is a convex combination of nonnegative pieces; the proposed assertion will pass and is non-tautological (would fail if `sigma_xi` were ever changed to include a sign-changing term). (3) F3's interval-check hardening — confirmed the SymPy script currently does not raise on `False` while the Mathematica side does; the proposed assertion change keeps the printed output identical when the check passes. (4) F4's banner edit is a string-replace; no math impact. (5) Confirmed no `paper_misalignment` is auto-introduced by F2/F3/F4 fixes: positivity is asserted in both the card and the notes, so adding the check brings the script *closer* to the paper; F3 is a hardening only; F4 is cosmetic.
