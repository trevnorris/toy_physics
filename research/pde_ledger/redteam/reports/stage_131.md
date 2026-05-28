---
unit_id: 131
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md
  paper_appendix: present
---

# Audit unit 131 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_131.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_131}` row at line 1296; no separate narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.txt`

## What the paper claims

Stage 131 establishes the parent micro-threshold for canonical mouth compensation. The card body quotes verbatim:

> "Canonical branch is the parent threshold \(T_m - q_* A'_0 = \Pi_*\Theta_\sigma/L\)."

The notes file makes four distinct deliverables explicit:
1. The compensation point value `\Pi_* \approx 1.50882951349316` defined as the root of `gPi(\Pi) = g_-^{F1}` for the Family-1 lower branch.
2. The boxed parent threshold identity `T_m - q_* A_0' = \Pi_* \Theta_\sigma / L \approx 1.509 \Theta_\sigma / L`.
3. The linearized deviation slope `g'(\Pi_*) \approx 0.0714453558083195`, giving `\mathfrak g_\Pi - \mathfrak g_-^{F1} \approx 0.0714... (\Pi_m - \Pi_*)`.
4. The collapsed structure: only the single linear combination `T_m - q_* A_0'` matters for the mouth-source bias on this boundary-layer branch.

The `Checks` field also enumerates three "must-check" items: positivity of the mouth source; zero-flux / boundary-layer normalizations in the GNLS/localized-Maxwell reduction; and the Family-1 compensation point versus the lower (not the singular equal-normalized) branch.

## What the script claims to verify

Both scripts compute `\Pi_*` by numerically solving `gPi(\Pi) - g_-^{F1} = 0` (SymPy `nsolve`, Mathematica `FindRoot`), then print `\Pi_*`, `V_{1,*} = \Pi_*\Theta_\sigma/L`, `g'(\Pi_*)`, and a symbolic "Parent bias mismatch formula" given by `(T_m - q_* A_0') - \Pi \Theta_\sigma / L`. The SymPy script contains **no assertions whatsoever** — only `print` statements. The Mathematica script contains exactly **one** asserting check (`expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20]`), which simply re-evaluates `gPi` at the root that `FindRoot` just returned. The "Parent bias mismatch formula" is printed as a free-symbol residual, never asserted against any value.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side coverage | Status |
|---|---|---|
| `\Pi_* \approx 1.50882951349316` (numeric value of compensation point) | Computed and printed by both scripts; **never asserted** against the paper-stated literal | partial |
| `T_m - q_* A_0' = \Pi_* \Theta_\sigma / L` (parent threshold identity) | The expression `(T_m - q_* A_0') - \Pi \Theta_\sigma / L` is printed as a symbolic residual; never substituted at `\Pi = \Pi_*` and never asserted | partial |
| `g'(\Pi_*) \approx 0.0714453558083195` (linearized slope) | Printed by both scripts; **never asserted** against the paper-stated literal | partial |
| Check: positivity of mouth source | `Pi, Theta_sigma, L` declared positive in symbol assumptions; no positivity test on a mouth-source expression | partial |
| Check: zero-flux / boundary-layer normalizations (GNLS / localized-Maxwell) | not exercised | missing |
| Check: F1 compensation point against lower branch (not singular branch) | `gMinus` hardcoded as `0.758035...` (sympy) or `(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi)` (mathematica) without anchor to an upstream-stage derivation; lower-branch-vs-singular-branch distinction not enforced | partial |

`paper_alignment` set to `partial`: every paper-side deliverable is touched by the script in some way (values are computed and printed), but none are converted into a failing-on-violation assertion, so the audit cannot certify any of them.

## Assertion inventory

| #  | Script       | Line | Form                                                                                  | Exercises which paper claim?      | Anchored to claim? |
|----|--------------|------|---------------------------------------------------------------------------------------|-----------------------------------|--------------------|
| A1 | sympy        | —    | (none — only `print` statements)                                                       | none                              | no                 |
| A2 | mathematica  | 46   | `expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20]`         | claim 1 (Pi_*) — tautologically  | no                 |

A2 is tautological because `piStar` is *defined* at line 35 as the root of `FindRoot[gPi == gMinus, ...]`; re-evaluating `gPi` there and comparing to `gMinus` cannot fail short of FindRoot returning a value outside its accuracy tolerance, and that is a check on FindRoot, not on the physics.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Subtype:** `script_doesnt_cover_claim`
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:1-27`

**What's wrong:**
The SymPy script contains zero assertions. `grep -n "assert\|expect\|raise" ...` returns nothing. The entire body is symbol declarations plus four `print` statements:

```python
print("Pi_* =", sp.N(Pi_star, 30))
print("V1_* =", V1_star)
print("g'(Pi_*) =", gprime_star)
print("Parent bias mismatch formula =", threshold_residual)
```

A SymPy "audit" with no assertions cannot fail on any physics violation. The output transcript shows `Status: PASS` and `EXIT_CODE: 0` purely because the script never tested anything.

**Why this matters:**
The paper card states this stage is a verification step ("Stage~131 is a positive mouth-source and boundary layer ledger step. Its audit target is the verification output quoted below."). The notes file lists at least three numerical/symbolic deliverables (Pi_* value, parent threshold identity, g'(Pi_*) value). None of those are converted into checks that can fail. The script currently exits 0 even if every printed number were wrong.

**Required change:**
Add `assert` lines for each paper-stated deliverable. At minimum:
- Assert `abs(sp.N(Pi_star, 50) - sp.Float("1.50882951349316", 50)) < sp.Float("1e-14", 50)` (anchors Pi_* to paper's quoted value).
- Assert that `gPi` evaluated at `Pi = Pi_star` is consistent with `g_minus` within tolerance — but only AFTER recomputing `Pi_star` from an independent route (see F4) or by re-checking via residual; otherwise this is the same tautology as F2.
- Assert `abs(gprime_star - sp.Float("0.0714453558083195", 50)) < sp.Float("1e-14", 50)`.
- Assert `sp.simplify(threshold_residual.subs(Pi, Pi_star) - (Tm - qstar*A0p - sp.Float(Pi_star,50)*Theta_sigma/L)) == 0` (binds the printed formula to the paper's identity at the compensation point).

Each assertion must be `assert <cond>, "<descriptive message>"` and the messages must name the paper anchor (e.g., `"Pi_* mismatch vs notes stage_131 Sec.1"`).

**Verification:**
After the fix, `redteam exec-sympy 131` runs to exit code 0 only if all four assertions pass. The output transcript must contain four `PASS:` lines (or equivalent) tied to the four deliverables.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:35-46`

**What's wrong:**
Line 35 defines `piStar` as the numerical root of `gPi == gMinus`:

```mathematica
piStar = piM /. FindRoot[gPi == gMinus, {piM, 1.5}, WorkingPrecision -> 80, ...];
```

Line 46 then asserts:

```mathematica
expectApprox["Pi_* compensation point", gPi /. piM -> piStar, gMinus, 10^-20];
```

This is the only `expectApprox` / `assert`-equivalent in the file. By construction, `gPi` evaluated at `piStar` equals `gMinus` to whatever accuracy `FindRoot` reports — the comparison cannot fail short of a numerical solver malfunction. The check tests `FindRoot`'s convergence, not any physics claim of the stage. The output transcript confirms: residual `≈ 8.92e-58` — vacuously below the `10^-20` tolerance because `piStar` was defined to make it so.

Secondary cosmetic issue (not separately filed): line 26 banner reads `"STAGE 114 — PARENT MICRO-THRESHOLD..."` but this is stage 131. The output transcript echoes the wrong stage label.

**Why this matters:**
A tautological assertion gives false confidence: the transcript labels this `PASS: Pi_* compensation point`, but nothing about the physics has been verified. None of the four paper-side deliverables (Pi_* numeric value, parent threshold identity, g'(Pi_*), linearization formula) have a non-tautological check.

**Required change:**
Replace the tautological `expectApprox` with substantive checks corresponding to the paper-side deliverables (parallel to F1's SymPy assertions, but written independently of the SymPy choreography to satisfy the second-engine independence policy):
- `expectApprox["Pi_* paper value", N[piStar, 50], 1.50882951349316, 10^-14]` (anchors to the paper's quoted Pi_*).
- `expectApprox["g'(Pi_*) paper value", N[D[gPi, piM] /. piM -> piStar, 50], 0.0714453558083195, 10^-14]`.
- A check that, with `piM -> piStar` substituted into `thresholdResidual`, the symbolic form reduces to `tM - qStar*a0Prime - piStar*thetaSigma/lM` (which is the paper-stated identity at the compensation point).

Also fix the line 26 banner string from `"STAGE 114"` to `"STAGE 131"`.

**Verification:**
After the fix, the `.wl` transcript must show: (a) `PASS: Pi_* paper value` with residual at numerical-precision level, NOT trivially zero by construction; (b) `PASS: g'(Pi_*) paper value`; (c) a `PASS` on the threshold-identity reduction; (d) the banner reads `STAGE 131`.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:25-26`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:40-45`

**What's wrong:**
Both scripts print the "Parent bias mismatch formula" as a free-symbol expression:

- SymPy line 25: `threshold_residual = sp.simplify((Tm - qstar*A0p) - Pi*Theta_sigma/L)`
- Mathematica line 40: `thresholdResidual = FullSimplify[(tM - qStar*a0Prime) - piM*thetaSigma/lM, ...]`

Neither substitutes `Pi -> Pi_*` and neither asserts anything. The resulting printed object `-A0p*q_* + T_m - Pi*Theta_sigma/L` is identically the LHS minus RHS of the canonical-branch identity from notes Sec.2 by construction; it cannot vanish for arbitrary `Pi`, and the script never picks the `Pi = Pi_*` slice that the paper claims is the canonical branch. The paper's `Checks` list also calls for (i) zero-flux / boundary-layer normalization tests and (iii) lower-branch vs singular-branch enforcement — neither is exercised in either script.

**Why this matters:**
The paper's bottom-line claim is the threshold *at* `Pi = Pi_*`. Without that substitution and without the lower-vs-singular branch enforcement, the scripts only restate the *definition* of the bias parameter, not the claim that compensation occurs *at* the specific Pi_* value.

**Required change:**
- SymPy line 26 should be replaced (or augmented) by an assertion of the form:
  `assert sp.simplify(threshold_residual.subs(Pi, Pi_star) - (Tm - qstar*A0p - sp.N(Pi_star, 50)*Theta_sigma/L)) == 0, "parent threshold identity at Pi_*"`. (Note: the equation is structurally an identity once Pi is substituted; the value of the assertion is to anchor `Pi_*` numerically.)
- Add a lower-branch sanity check: assert that `gPi` evaluated at `Pi_star` matches the F1 lower-branch value `g_minus` AND that `gPi` evaluated at a clearly different Pi (e.g., `Pi = 2*Pi_star`) does NOT match `g_minus` — this distinguishes the chosen root from the trivial / singular branches mentioned in the paper's Checks.
- Mirror both in the Mathematica file.

**Verification:**
After the fix, both `.txt` transcripts must show a PASS line referencing the parent threshold identity at Pi_* AND a PASS line confirming that a counter-example point (e.g., 2*Pi_star) gives a clearly nonzero residual against g_minus.

### F4 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:12`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:33`

**What's wrong:**
- SymPy line 12: `g_minus = sp.Float("0.758035078944663")` — bare numeric literal with no citation to a source stage or derivation.
- Mathematica line 33: `gMinus = N[(2*Sqrt[4107 - 100*Pi^2] - 37*Sqrt[3])/(20*Pi), 80]` — exact-form expression with no citation.

The notes file refers to `g_-^{F1}` but does not state its numeric value or closed form. The closed form in the `.wl` (involving `Sqrt[4107 - 100*Pi^2]` and `37*Sqrt[3]`) is a Family-1 lower-branch construction that must come from an upstream stage (e.g., Stages 128–130) which neither script cites. Without anchor, the load-bearing `Pi_*` value depends on a hardcoded import that the audit cannot trace.

Additionally, the two engines independently use *different forms* (a float literal vs an exact closed form). They agree numerically (`0.758035...`), but the SymPy float was apparently transcribed rather than independently constructed, weakening the second-engine independence claim.

**Why this matters:**
`Pi_*` is the root of `gPi(Pi) = g_minus`. If `g_minus` is wrong (typo, stale upstream, sign error), `Pi_*` is wrong, and every downstream stage (133–145 per the card's Downstream Use field) that quotes `Pi_* ≈ 1.50882951349316` is compromised. The script must anchor `g_minus` to a verifiable source (either re-derive from the F1 lower-branch closed form OR cite the upstream-stage file by name in a comment AND assert the value matches the closed form within tolerance).

**Required change:**
- In both scripts, add a comment immediately above the `g_minus` / `gMinus` line citing the upstream source stage (e.g., `# F1 lower-branch value carried from Stage NNN`). The auditor cannot identify the exact upstream stage without reading other units; the implementer must fill this in.
- In the SymPy script, replace the bare float by constructing it from the same closed form as the Mathematica side, then assert numerical agreement:
  `g_minus_exact = (2*sp.sqrt(4107 - 100*sp.pi**2) - 37*sp.sqrt(3)) / (20*sp.pi)`
  `g_minus = sp.N(g_minus_exact, 50)`
  `assert abs(g_minus - sp.Float("0.758035078944663", 50)) < sp.Float("1e-14", 50), "g_minus closed form vs printed literal"`.
- In the Mathematica script, add an `expectApprox` confirming the closed-form numeric value matches the literal `0.758035078944663` printed in the prior stage transcript.

**Verification:**
After the fix, both transcripts contain an explicit `g_minus` consistency check. A grep for `4107` and `100*Pi^2` (or `100*sp.pi**2`) appears in both engines, demonstrating independent (not transliterated) construction from the same closed form.

## Independent-derivation check (Mathematica)

Both scripts share the same numerical strategy (define `gPi`, solve for its intersection with `gMinus`, print `Pi_*`, `V_{1,*}`, `g'(Pi_*)`, and a symbolic "parent bias mismatch formula"). The intermediate variable names map line-for-line:

- SymPy line 14: `gPi = sp.simplify(2*Pi*(2*Pi*sp.exp(Pi) + sp.pi)/((4*Pi**2 + sp.pi**2)*(sp.exp(Pi) - 1)))`
- Mathematica line 34: `gPi = FullSimplify[2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1)), ...]`

- SymPy line 17: `V1 = Pi*Theta_sigma/L`
- Mathematica line 37: `v1 = piM*thetaSigma/lM`

- SymPy line 25: `threshold_residual = sp.simplify((Tm - qstar*A0p) - Pi*Theta_sigma/L)`
- Mathematica line 40: `thresholdResidual = FullSimplify[(tM - qStar*a0Prime) - piM*thetaSigma/lM, ...]`

The choreography is the same: variable-by-variable transliteration in Mathematica syntax. The only structural divergence is that the Mathematica script uses the exact closed form for `gMinus` while the SymPy script uses a float literal — but this is a difference in import-data form, not in derivation.

This is a borderline `mathematica_transliteration` case. However, the stage's verification target is shallow (compute one root and print four quantities), so the "two independent derivations" requirement has limited room to manifest differently. I am NOT filing a separate `mathematica_transliteration` finding because (a) the fix for F4 will introduce some independent construction in the SymPy side from the closed form, and (b) the more pressing issue is that both engines are doing essentially nothing assertion-wise (F1, F2), making the transliteration discussion moot until the assertions exist.

## Engine cross-check

Both engines agree numerically:

- `Pi_*` — SymPy: `1.50882951349315861144664988336`; Mathematica: `1.50882951349315558300555075595`. Difference at the ~30th digit, consistent with `nsolve` tol=1e-30 and Mathematica's `WorkingPrecision -> 80`.
- `g'(Pi_*)` — SymPy: `0.0714453558083194772127035664826`; Mathematica: `0.071445355808319521189460301988136`. Same magnitude, differing only in trailing digits.
- "Parent bias mismatch formula" — SymPy prints `-A0p*q_* + T_m - Pi*Theta_sigma/L`; Mathematica prints `-(a0Prime*qStar) - (piM*thetaSigma)/lM + tM`. Symbolically identical up to symbol naming.

No `engine_disagreement` finding.

## Verdict justification

This is a `findings` verdict with four substantive issues. The SymPy script (F1) contains zero assertions — it is a verification script in name only. The Mathematica script's sole assertion is tautological (F2): it re-checks the defining equation of `Pi_*` against a value `Pi_*` was just constructed to satisfy. Two of the paper's load-bearing numeric outputs (`Pi_* ≈ 1.50882951349316`, `g'(Pi_*) ≈ 0.0714453558083195`) are printed but not asserted against the paper's quoted literals. The lower-branch input `g_minus` is hardcoded (literal float on the SymPy side, exact closed form with no anchor on the Mathematica side) in a way that cannot be traced to an upstream stage from these scripts alone (F4). The "Parent bias mismatch formula" is printed as a symbolic identity but never evaluated at the compensation point `Pi = Pi_*` (F3).

The paper's claims are consistent and internally well-defined; the math is not broken. The defect is verification coverage, not physics. No `stop_cold` flag fires: nothing here propagates a wrong value downstream (the printed numerics agree with the notes' quoted Pi_* to 30 digits, so downstream consumers reading the transcript get the right number — it's just that the audit didn't prove it).

`paper_alignment: partial` because every paper-side deliverable is touched by some script-side computation, but no deliverable is converted into a failing-on-violation assertion.

## Self-test notes

- **Variable independence**: For the proposed `sp.diff(gPi, Pi)` in F1's verification (already present in the original script at line 22), `gPi` is a function of `Pi` only — the differentiation is well-formed. No diff trap.
- **Symmetry/parity**: No integrals over unbounded domains in this stage; not applicable.
- **Trivial-case pre-check**: For F3's counter-example assertion (`gPi[2*Pi_star] != g_minus`), mental substitution at large Pi gives `gPi ~ 2*Pi^2 / (4*Pi^2) = 1/2 = 0.5`, clearly different from `0.758`, so the proposed `assert nonzero residual` is correct in form.
- **Paths**: F1, F2, F3, F4 all target existing files in their correct directories (`scripts/` for `.py`, `mathematica/` for `.wl`). No new files required.
- **Paper round-trip**: The required numeric anchors (`Pi_* ≈ 1.50882951349316`, `g'(Pi_*) ≈ 0.0714453558083195`) are directly quoted in `notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md` lines 8 and 93 — no new constants introduced.
