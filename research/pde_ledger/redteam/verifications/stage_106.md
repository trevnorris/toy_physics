---
unit_id: 106
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 106

Re-verification against the current directive, which now carries
`## Applied: F1..F4` blocks. This supersedes the 2026-05-27 first pass, whose F1
write-up cited the pre-correction attribution ("Stage 104 is the source of
chi_Q = 1"); that line has since been corrected to Stage 105 and is re-checked
below.

## Per-finding outcomes

### F1 — paper_misalignment (resolved, doc-only)

**Classification:** resolved

**What changed:**
The directive RESOLVED F1 via direction (b)-style: chi_Q=1 is a legitimate
carry-in, no Hankel/z-expansion is added to stage 106, and the only script
action was a citation correction.
- `scripts/...sympy_audit.py:5-13` — docstring now attributes item (ii) higher
  odd → **Stage 102**; item (iii) outgoing l=2 DtN fingerprint → **Stage 104**;
  **chi_Q=1 fixed → Stage 105** (`chiQ_fix_from_outgoing_dtn`). The diff
  (patch lines 101-105) shows the stale "Stage 104 ... is the source of
  chi_Q = 1" wording removed and replaced. No residual "Stage 88" reference.
- `mathematica/...mathematica_audit.wl:29-32` — header comment matches: Check
  (ii) → stage 102; Stage 104 proves the DtN fingerprint; Stage 105 fixes
  chi_Q=1. (Diff patch lines 13-16 show the prior "stages 102 and 104" lumped
  wording was corrected to split 104=fingerprint / 105=chi_Q.)

**Assessment:**
Correct and complete, matching the resolved-directive requirement exactly
(ii→102, iii DtN→104, chi_Q=1→105). NO Hankel/z-expansion verification was
added to stage 106 (no `Lambda_2`, no `z = omega a/c_s` expansion in either
file). The diff touches only the two target scripts — no paper.tex or notes
edits appear in the patch. chi_Q=1 remains a carry-in. Properly scoped doc-only
fix.

### F2 — tautological_check (resolved)

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:65-72` — replaced the old `K4 - 4*K2**2/K0` and
  `Gamma5 - 9*sqrt(K2**5/K0**3)` checks (tautological because K2, K4 were
  defined as `K0/(4Ω²)`, `K0/(4Ω⁴)`) with checks on the hardcoded literals:
  `K0_target*K4_target - 4*K2_target**2` and
  `Gamma5_target - 9*sqrt(K2_target**5/K0_target**3)`.
- `mathematica/...audit.wl:86-93` — equivalent target-literal checks.

**Assessment — LOAD-BEARING target-literal independence judgment:**
Each `*_target` literal is set INDEPENDENTLY (sympy lines 50-53; wl lines
70-73); none is derived from another:
- `K0_target = 54·G·c_s⁵/(5·a⁵·c⁵)`
- `K2_target = 6·G·c_s³/(5·a³·c⁵)`
- `K4_target = 8·G·c_s/(15·a·c⁵)`
- `Gamma5_target = 2·G/(5·c⁵)`

`K0_target·K4_target - 4·K2_target² = 0` is a GENUINE (non-constructed) numerical
agreement: `K0·K4 = (54/5)(8/15) = 432/75 = 144/25` and `4·K2² = 4·(6/5)² =
144/25` (both times the common factor `G²c_s⁶/(a⁶c¹⁰)`). The two
independently-written rational prefactors happen to satisfy K0·K4 = 4K2². The
old check collapsed to identity for ANY K0; the new check would FAIL under a
perturbation of any one literal (e.g. K2_target → 7/5·...). Non-tautological.

### F3 — mathematica_transliteration (resolved)

**Classification:** resolved

**What changed:**
The `.wl` derives the bottom line via an independent one-pole `Yret`
omega-series path:
`Yret[om,chi] = 3/4 + (1/4)/(1 - om²/OmegaQ² - I·chi·sigmaQcan·om⁵)`,
series-expanded in omega; reads off the ω⁵ coefficient; then solves
`Solve[oddScaleFromSeries == gamma5Target, chiQ]` to FIX chi_Q = 1
(wl lines 39-81). The F2 target-literal checks and F4 sensitivity ride on top.

**Assessment — independence judgment:**
Genuinely distinct, not a renamed mirror. The case-sensitive grep returns 0
occurrences of the `.py` intermediate names `nqGeneral`, `nqCanonical`, `k0`,
`k2`, `k4`, `gamma5` in the `.wl` (the only `…Target` symbols are hardcoded
literals, not constraint-derived intermediates). The `.wl` uses `Yret`,
`seriesY`, `OmegaQ`, `sigmaQcan`, `omegaSym`, `oddScaleFromSeries`,
`chiFromOmega5Match`, `sourceNormalizer` — different starting object (retarded
one-pole), different operation (omega-series), different derivation order. It
never calls `Solve[constraint, NQ]` nor performs a branch pick, so a SymPy
branch/sign bug could not be reproduced verbatim. It actually does MORE than the
`.py`: it DERIVES chi_Q=1 from the ω⁵ match rather than importing it, which also
answers the codex_reviews R3 concern. Assertions are non-tautological — the
ω⁵-form match, the chi_Q=1 solve, and the ω⁷-onset check
(`omega7Coeff = (8I/243)a⁷/cS⁷` vs the independently-written closed form
`I·sigmaQcan/(2·OmegaQ²)`, which I confirmed evaluates to `(8I/243)a⁷/cS⁷`) each
test a real relation between the series and an independent target.

### F4 — insufficient_verification (resolved)

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:80-97` — Delta_Q sensitivity: expand
  `gamma_eff_off = (m0hat²·Gamma5).subs(m0hat→1, chi_Q→1+Delta_Q)` to first
  order; assert `linear_coeff == -2G/(5c⁵)` and `zeroth_coeff == Gamma5_target`.
- `mathematica/...audit.wl:109-124` — equivalent Delta_Q first-order block.

**Assessment — non-triviality:**
Non-trivial. On the natural branch gamma_eff = Gamma5_target/(chi_Q·m0hat²); at
m0hat=1, chi_Q=1+Delta_Q this is Gamma5_target·(1 - Delta_Q + O(Delta_Q²)), so
the first-order slope = -Gamma5_target = -2G/(5c⁵) and the zeroth =
Gamma5_target. The check exercises the actual series expansion of 1/(1+Delta_Q);
a sign flip or factor error in N_Q would change the linear coefficient and break
the assertion. Both engines print the two residuals = 0 (sympy log lines 21-22;
mathematica log lines 33-36, both PASS). Not a 0==0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `target identity K0_target K4_target - 4 K2_target^2 = 0`
- `target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3) = 0`
- `Delta_Q first-order sensitivity coefficient = 0`
- `Delta_Q zeroth-order coefficient equals Gamma5_target = 0`

**Mathematica:** exit=0, 11 PASS / 0 FAIL. All 11 are content-bearing (distinct
named residuals, no 0==0):
`sigma_Q^can - 4 a^5/(27 c_s^5)`; `omega^5 coefficient form`;
`first next-odd omega^7 coefficient`; `omega^5 coefficient gives chi_Q
Gamma5_target`; `chi_Q fixed by canonical omega^5 match`;
`target identity K0_target K4_target - 4 K2_target^2`;
`target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3)`;
`N_Q on natural branch at m0hat=1, chi_Q=1`; `canonical gamma_eff - target`;
`Delta_Q first-order sensitivity coefficient`;
`Delta_Q zeroth-order coefficient equals Gamma5_target`.

Engine cross-check: SymPy reaches `Gamma5 = 2G/(5c⁵·chi_Q·m0hat²)` →
gamma_eff = 2G/(5c⁵) on the canonical branch; Mathematica reaches the same via
the retarded one-pole series, and both report the Delta_Q slope -2G/(5c⁵). The
two engines now reach the bottom line via structurally different paths, so the
agreement is substantive.

**Output freshness:** confirmed. Source mtimes (py 1780073275, wl 1780073348)
predate the saved-output mtimes (py .txt 1780073489, wl .txt 1780073508). The
saved sympy `.txt` is byte-identical to the exec-log body (diff returned no
difference).

## Material-change assessment

`material_change`: false.

All four fixes add or correct verification depth (target-literal consistency,
independent `.wl` derivation, Delta_Q sensitivity, citation accuracy). None
changes a derived numerical result: the bottom line is still N_Q = 1 on the
canonical branch and gamma_quad^eff = 2G/(5c⁵). chi_Q=1 remains an upstream
carry-in (Stage 105). No downstream unit is affected.

## Side observations (non-blocking)

- The `.wl` prints (but does not assert on) the ω⁶ coefficient `(16a⁶)/(729cS⁶)`.
  Print-only; the load-bearing odd checks are on ω⁵ and ω⁷. Acceptable.
- codex_reviews/stage_106.md raised R1 (paper card still assigning ii/iii to
  stage 106). R1 is out of scripts-only scope and is logged to
  PAPER_CLEANUP_TRACKER for the manual paper pass per the resolved F1 direction.
  Its R2 (ω⁷ print-only) and R3 (`.wl` imports chi_Q=1) are ALREADY addressed by
  the applied edits: the `.wl` now has an explicit ω⁷-coefficient assertion and
  derives chi_Q=1 from the ω⁵ match. Non-blocking.
- Both scripts now banner `STAGE 106 — REDUCED 2.5PN CLOSURE ON CANONICAL
  OUTGOING DtN BRANCH`; the auditor's stale-banner observation is moot.

## Verdict justification

All four findings are resolved. F1 is a correctly-scoped doc-only citation fix
(ii→102, iii DtN→104, chi_Q=1→105) with no Hankel/z-expansion and no paper/notes
edits. F2's new checks test four independently-hardcoded literals whose
K0·K4 = 4K2² (144/25 = 144/25) is a genuine non-constructed agreement — not a
renamed tautology. F3's `.wl` is a structurally independent one-pole `Yret`
omega-series derivation (0 occurrences of the `.py` lowercase names; it even
derives chi_Q=1 rather than importing it). F4's Delta_Q first-order slope
(-2G/(5c⁵)) exercises the series algebra non-trivially in both engines. SymPy
exits 0; Mathematica exits 0 with 11 content-bearing PASS / 0 FAIL; outputs are
fresh; the diff shows no collateral edits beyond the two target scripts.
Verdict: verified.
