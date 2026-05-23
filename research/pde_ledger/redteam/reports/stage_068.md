---
unit_id: 068
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 068 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.txt`

## What the script claims to verify

The docstring claims the audit covers four results: (1) the resonance-corrected wall figure `W_res(r) = C^2(r) W_wall`; (2) "exact threshold translation from the matched branch to the profile family"; (3) the explicit resonance-point penalty `P_res = 1/C_res^2`; (4) the "exact width of the profile-sensitive threshold band". The actual assertions exercise none of these as physical content. They reduce to algebraic identities between symbols the script itself just defined three lines earlier: e.g. defining `Wfail_res := Pe_req/(C2*Delta_inf)` and then asserting that substituting `C2 = 1/P_res` yields `P_res*Pe_req/Delta_inf`, which is a one-line algebraic step. No PDE, no resonance condition, no profile family, and no derivation of `P_res = 1/C_res^2` is exercised; the relationship is hard-coded into the substitution `C2 -> 1/P_res`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 67 | `expect_zero(Wfail_res.subs(C2, 1/Pres) - Pres*Wfail_match)` | no (tautological — `Pe_req/((1/Pres)*Delta_inf) - Pres*Pe_req/Delta_inf` is `0` by one algebraic step) |
| A2 | sympy | 68 | `expect_zero(Wsuff_res.subs(C2, 1/Pres) - Pres*Wsuff_match)` | no (same tautology with `Delta_0`) |
| A3 | sympy | 80 | `expect_zero(success_band_width/Wsuff_match - (Pres - 1))` | no (`success_band_width := Pres*Wsuff_match - Wsuff_match`, so the ratio is identically `Pres-1` by construction) |
| A4 | sympy | 81 | `expect_zero(failure_band_width/Wfail_match - (Pres - 1))` | no (same tautology) |
| A5 | mathematica | 46 | `expectZero[(WfailRes /. C2 -> 1/Pres) - Pres*WfailMatch]` | no (mirror of A1) |
| A6 | mathematica | 47 | `expectZero[(WsuffRes /. C2 -> 1/Pres) - Pres*WsuffMatch]` | no (mirror of A2) |
| A7 | mathematica | 57 | `expectZero[successBandRes/WsuffMatch - (Pres - 1)]` | no (mirror of A3; the substitution `C2 -> 1/Pres` collapses `successBand` to `(Pres-1)*WsuffMatch`, so the ratio is `Pres-1` by construction) |
| A8 | mathematica | 58 | `expectZero[failureBandRes/WfailMatch - (Pres - 1)]` | no (mirror of A4) |

Every assertion is tautological. There are zero non-trivial assertions in the file.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:55-68`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:36-47`

**What's wrong:**
The "threshold translation" checks are algebraic identities built directly from the script's own definitions.

In SymPy:
```
Wfail_match = sp.simplify(Pe_req / Deltainf)            # line 55
Wfail_res   = sp.simplify(Pe_req / (C2 * Deltainf))     # line 58
expect_zero("Wfail_res - P_res*Wfail_match",
            Wfail_res.subs(C2, 1 / Pres) - Pres * Wfail_match)   # line 67
```

Direct substitution gives `Wfail_res.subs(C2, 1/Pres) = Pe_req/((1/Pres)*Delta_inf) = Pres*Pe_req/Delta_inf`, which equals `Pres*Wfail_match` by one algebraic step. The assertion can never fail because `Wfail_res` was defined as `Pe_req/(C2*Delta_inf)` two lines above the assertion. The Mathematica script (`mathematica/...stage068...wl:36-47`) repeats the same identity with renamed symbols.

The script's docstring (line 11-12) advertises this as "Exact threshold translation from the matched branch to the profile family", but no matched-branch threshold is computed from any physical premise — both `Wfail_match` and `Wfail_res` are postulated algebraic expressions, and `P_res = 1/C^2` is hard-coded into the substitution rule.

**Why this matters:**
Stage 068 is non-checkpoint and `is_status_only_candidate: false`, so the manifest requires substantive verification. Right now the audit script provides zero evidence that the resonance correction `W_res = C^2 W_wall`, the penalty `P_res = 1/C_res^2`, or the "translation" between matched and profile-family thresholds is physically derivable. Anyone changing `C^2`-vs-`C` conventions, a sign, or the dependence on `r` upstream would still see this script pass.

**Required change:**
Replace the substitution-identity checks with a derivation that ties the resonance correction to a non-trivial premise. The minimum acceptable upgrade has two parts (both engines):
1. Derive `W_res = C^2 W_wall` from a stated physical premise, e.g. a transmission/standing-wave amplification factor `C(r)` acting on a power balance `W_wall * |something|^2 = ...`. The premise should be written symbolically and the algebraic step should be `expect_zero(W_res_derived - C2 * W_wall)`.
2. Derive `P_res = 1/C_res^2` (the resonance-point inverse-square scaling) from the inverse of the amplification factor. The check should look like `expect_zero(Pres_derived - 1/Cres**2)` where `Pres_derived` is computed from the symbolic transmission setup, not declared.

Failing that, mark the unit as `is_status_only_candidate: true` in the manifest. But that is a structural change outside the script's scope — Codex should instead add new derivation steps and replace A1, A2, A5, A6 with their non-tautological versions.

**Verification:**
After Codex applies, A1/A2 (SymPy lines 67-68) and A5/A6 (Mathematica lines 46-47) should each be replaced by an `expect_zero` whose left-hand side is computed from a symbolic premise rather than substituted into the very definition the assertion references. The substituted expression and the comparison expression must not both trace to the same single algebraic definition.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:74-81`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:50-58`

**What's wrong:**
The "profile-sensitive band width" assertions are trivially zero. In SymPy:
```
success_band_width = sp.simplify((Pres * Wsuff_match) - Wsuff_match)        # line 74
failure_band_width = sp.simplify((Pres * Wfail_match) - Wfail_match)        # line 75
expect_zero("success width / matched threshold",
            sp.simplify(success_band_width / Wsuff_match - (Pres - 1)))     # line 80
```

`success_band_width / Wsuff_match = (Pres*Wsuff_match - Wsuff_match)/Wsuff_match = Pres - 1`, so the residual is `(Pres-1) - (Pres-1) = 0`. The assertion is `(Pres-1) == (Pres-1)`, which is true for every symbol.

Mathematica's analog (lines 50-58) takes a slightly different route — defines `successBand = WsuffRes - WsuffMatch` in terms of `C2`, then applies the substitution `C2 -> 1/Pres` before forming the ratio. After substitution: `successBandRes = PeReq*(1/Pres - 1)/((1/Pres)*Delta0) * (-1)` correction — more simply, the ratio reduces to `Pres - 1` by direct algebra. The assertion `expectZero[successBandRes/WsuffMatch - (Pres - 1)]` therefore evaluates `0 - 0 = 0` for the same reason.

**Why this matters:**
Stage 068's docstring (line 14) lists "Exact width of the profile-sensitive threshold band" as one of four things being checked. The current assertions don't exercise any width formula — they only verify that `(P-1)x/x = P-1`. The author's intended claim (presumably `W_suff^(res,*) - W_suff^(match) = (P_res - 1) * W_suff^(match)`) is true only when both endpoints are computed from a common physical premise; the audit script simply asserts the identity between two postulated symbols.

**Required change:**
Either:
(a) Anchor the band-width expressions to the same physical premises that should now drive F1's fix. Define `Wsuff_match` from a matched-branch power balance (e.g. `Pe_req / Delta_0` where `Delta_0` is derived from a Pe→0 limit of the matched solution) and `Wsuff_res` from the resonance-corrected version, so that the assertion `Wsuff_res - Pres * Wsuff_match` (at C2 = 1/Pres) tests a non-trivial relation between two independently derived expressions.
(b) Remove the band-width block entirely if no independent premise is available, and update the docstring accordingly. But removing claims requires touching the docstring of the script (lines 14, 84-95), which Codex may do as part of its script-only edits.

The minimum substantive upgrade: replace lines 74-75 of the SymPy script and lines 50-56 of the Mathematica script so that `success_band_width` and `failure_band_width` are computed from `Wsuff_res - Wsuff_match` directly (which gives `Pe_req*(1-C2)/(C2*Delta_0)` in raw form), and the assertion checks that the simplified raw form equals `(Pres - 1)*Wsuff_match` only after a one-step substitution `C2 -> 1/Pres` chained through the raw form, not through a definition that already contains `Pres - 1` as a factor. Even this is borderline tautological if `Wsuff_res` was itself postulated rather than derived; F1's fix should be in place first.

**Verification:**
After Codex applies, A3/A4 (SymPy lines 80-81) and A7/A8 (Mathematica lines 57-58) should each compare a symbolic expression computed from the F1-derived `Wsuff_res`/`Wfail_res` against an independently-derived band-width formula. The two sides must not collapse to the same algebraic skeleton in one substitution step.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage068_resonance_thresholds_mathematica_audit.wl:26-58` (compared against `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:37-81`)

**What's wrong:**
The Mathematica script is a near line-by-line transliteration of the SymPy script rather than an independent re-derivation. Concrete evidence:

1. Identical helper choreography:
   - SymPy lines 30-34 define `expect_zero(name, expr)` that prints and raises on nonzero.
   - Mathematica lines 20-24 define `expectZero[name_String, expr_]` with the same signature, same printed format, and same fail-on-nonzero behavior.

2. Identical symbol set with mirrored names: SymPy uses `C2, W_wall, Pe_req, Delta_0, Delta_inf, P_res`; Mathematica uses `C2, Wwall, PeReq, Delta0, Deltainf, Pres`. Mathematica's `Wres` symbol is introduced and then immediately overwritten by `FullSimplify[C2*Wwall]` (line 33), matching SymPy's line 48 `W_res = sp.simplify(C2 * W_wall)`.

3. Identical sequence and labels: both scripts compute the four thresholds (matched fail/succeed, profile fail/succeed) in the same order with the same printed labels (`"Matched fail threshold"`, `"Matched succeed threshold"`, `"Profile-family fail thresh"`, `"Profile-family succ thresh"` — verbatim strings appearing in SymPy lines 61-64 and Mathematica lines 41-44).

4. Identical assertion sequence: both scripts then call expect_zero on `Wfail_res(C2 -> 1/Pres) - Pres*Wfail_match`, then on the `Wsuff` analog, then on the success/failure band widths. The order, names, and expression structure match one-for-one.

5. Identical "FINAL LEDGER" prose: SymPy lines 84-95 and Mathematica lines 60-72 print the same multi-line interpretation text verbatim.

A genuine second-engine derivation would solve a transmission/wave problem in Mathematica idiom (`DSolve`, `Reduce`, complex-amplitude algebra) and then compare its symbolic result to the SymPy script's `W_res = C^2 W_wall`. The current .wl does no derivation at all; it restates the same algebraic identities.

**Why this matters:**
The second-engine policy exists to catch single-engine bugs (CAS-specific simplification artifacts, branch-cut handling differences). A transliteration cannot catch any such bug because both engines walk the same algebraic chain. Together with F1/F2, this means stage 068 currently has zero independent verification.

**Required change:**
After F1 lands (which forces the SymPy script to derive `W_res = C^2 W_wall` and `P_res = 1/C_res^2` from a stated transmission/amplification premise), rewrite the Mathematica script so it re-derives the same two results using a different algebraic path — for example, by solving a 1-D Helmholtz-with-impedance-jump problem symbolically with `Solve`/`Reduce` to obtain the transmission coefficient `C(r)`, then asserting that `|C|^2 * W_wall` equals the SymPy-side `W_res` and `1/|C|^2` equals `P_res`. The structural rule: the .wl must not import or restate the symbolic formulas from the .py; it must produce them.

Codex can implement this change concurrently with F1's fix, since both engines need new derivation content.

**Verification:**
After Codex applies, the verifier should be able to point to (a) at least one Mathematica-only construct (e.g. `Solve[...]`, `DSolve[...]`, `Eliminate[...]`) that derives the transmission factor from a physical premise, and (b) a comparison between the Mathematica-derived `W_res` (or `P_res`) and the literal symbolic form `C2*W_wall` (or `1/C2`). If the .wl still consists only of `FullSimplify` on postulated expressions, the finding is not resolved.

## Independent-derivation check (Mathematica)

The .wl is a transliteration, as detailed in F3. The variable choreography, the assertion sequence, the printed-label strings, and the "FINAL LEDGER" prose all match the .py line-for-line. There is no Mathematica-side derivation — `FullSimplify` is applied to expressions that were declared, not computed.

## Engine cross-check

Both engines pass with zero residuals on identically-structured tautologies. The SymPy output (`scripts/output/...stage068...txt:18-19, 22-23`) prints:
```
Wfail_res - P_res*Wfail_match = 0
Wsuff_res - P_res*Wsuff_match = 0
success width / matched threshold = 0
failure width / matched threshold = 0
```
The Mathematica output (`mathematica/output/...stage068...txt:18-21, 28-31`) prints the same four `0` residuals with `PASS:` lines. The agreement is real but uninformative — both scripts are checking the same tautologies in the same order.

## Verdict justification

`findings`. The script's assertions are uniformly tautological — each one is an algebraic identity between two expressions the script itself defined within a few lines. No physical premise is exercised: there is no transmission/amplification setup behind `W_res = C^2 W_wall`, no derivation of `P_res = 1/C_res^2` (the relation is hard-coded into the substitution `C2 -> 1/Pres`), and no independently-derived band width. The Mathematica script compounds the problem by transliterating the SymPy algebra rather than performing an independent derivation. Verdict is not `stop_cold` because the fix is local (re-anchor the assertions to a stated symbolic premise) and there is no evidence the current script's wrong claims propagate to downstream units with a sign or factor flip — every downstream check that consumes `P_res = 1/C_res^2` would face the same redress.

Attacks attempted that the script did survive: (i) symbol-domain attack — the positivity assumptions on `C2, Pe_req, Delta_0, Delta_inf, P_res` are consistent and don't hide branch issues here; (ii) missing-branch attack — the script doesn't claim coverage over `C2 < 1` vs `C2 > 1` differently, so there's no branch to miss; (iii) stale-output attack — both output files post-date the corresponding source files. The script fails the only attack that mattered: assertion non-triviality.

## Self-test notes

I checked: (a) variable independence — not applicable here, no `diff` calls; (b) symmetry/parity — not applicable, no integrals; (c) trivial-case pre-check — I substituted `Pe_req = Delta_0 = Delta_inf = 1, C2 = 1/Pres` into each `expect_zero` call and confirmed each reduces to `0` purely by algebra independent of `Pres`, which is the diagnostic for tautology; (d) path specifications — F1/F2/F3 target existing files in `scripts/` and `mathematica/` and no new file is requested, so the script-target rule is satisfied. F1's required-change describes new derivation steps but does not propose new files.
