---
unit_id: 143
batch: IV.5
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage143_equal_normalized_singular_limit.md"]
  paper_appendix: present
---

# Audit unit 143 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_143.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage143_equal_normalized_singular_limit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only line 1320 `\input{stages/stage_143}` references this unit; no narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.txt`

## What the paper claims

The card says (quoting from the displayed equation in the body): "$\mathfrak g_\Pi<1$ for all finite $\Pi>0$; equal-normalized branch requires $\Pi\to\infty$ and divergent traction." The notes spell out three distinct deliverables: (1) the strict inequality $0<\mathfrak g_\Pi<1$ for $\Pi>0$, established by an exact decomposition of the numerator of $1-\mathfrak g_\Pi$ into three manifestly positive pieces $\pi^2(e^\Pi-1-\Pi-\Pi^2/2)$, $\Pi(\pi^2-2\pi)$, $\Pi^2(\pi^2/2-4)$; (2) the singular limit $\lim_{\Pi\to\infty}\mathfrak g_\Pi=1$; (3) the traction-divergence asymptote $\widehat T_m(\Pi)\sim\sqrt{9/[20(1-R_\infty)]}\,\Pi^{1/2}\approx 0.7256691307\,\Pi^{1/2}$ with $R_\infty=(1-\mathfrak r_{F1})^2/(1+\mathfrak r_{F1}^2)\approx 0.1454544523$, requiring $\widehat T_m\to\infty$. The card's `Verification` field points only to the two scripts; there is no quantified `Output` block, so the notes are authoritative on the deliverable list.

## What the script claims to verify

The SymPy script (1) defines $\mathfrak g_\Pi$ and the carried-forward constant $r=\sqrt{4107-100\pi^2}/(10\pi)$, (2) multiplies $1-\mathfrak g_\Pi$ by $(4\Pi^2+\pi^2)(e^\Pi-1)$ to clear denominators, (3) asserts via `expect_zero` that the resulting numerator equals the three-piece decomposition $\pi^2(e^\Pi-1-\Pi-\Pi^2/2)+\Pi(\pi^2-2\pi)+\Pi^2(\pi^2/2-4)$, (4) prints (does not assert) the values of the three "positive pieces", (5) prints (does not assert) the limits $g_0$, $g_\infty$, $R_\infty$, $S_\infty$, $\Sigma_0/\Pi$, and $\widehat T_m/\sqrt\Pi$. The Mathematica script mirrors steps (1)-(3) and (5) but skips the actual computation of the $S_q$ and $\Sigma_0$ limits, hardcoding `sInf = 1` and computing `sigmaRatio` by direct algebraic substitution rather than as a limit. The only firing assertion in either script is the numerator/decomposition identity.

## Paper ↔ script cross-check

| Paper deliverable | Script-side coverage |
|---|---|
| (1a) Algebraic identity: numerator of $1-\mathfrak g_\Pi$ = three-piece decomposition | `match` — sympy line 32 and mathematica line 46 both assert this. |
| (1b) Each of the three pieces is positive for $\Pi>0$ (so $g_\Pi<1$) | `partial` — pieces are printed but never asserted positive. Mathematica's coeff prints (`(-2+Pi)*Pi`, `(-8+Pi^2)/2`) symbolically use Mathematica's builtin `Pi`, which is numerical $\pi$; the values are correct but the script asserts nothing about their sign. |
| (2) $\lim_{\Pi\to\infty}\mathfrak g_\Pi=1$ | `partial` — sympy line 41 prints `ginf` via `sp.limit`; mathematica line 56 prints `gInf`. Neither asserts the value equals 1. Output shows it does, but a regression to e.g. `1/2` would not fail the script. |
| (3a) $R_\infty\approx 0.1454544523$ | `partial` — sympy line 47 computes via `sp.limit`; mathematica line 61 computes via direct substitution $r\to r$ (not as a limit). No assertion on the numeric value. |
| (3b) $S_\infty=1$ | `partial` (sympy via limit) / `missing` (mathematica hardcodes `sInf = 1` on line 62 without taking a limit of `sQ`; `sQ` is defined but never used). |
| (3c) $\widehat T_m/\sqrt\Pi\to\sqrt{9/[20(1-R_\infty)]}\approx 0.7256691307$ | `partial` (sympy line 54-57 prints via `sp.limit(Sigma0/Pi, ...)`) / `missing` (mathematica never defines `Sigma0` or `That`; `tHatRatio` is computed by direct algebraic substitution `Sqrt[(9/20)*sigmaRatio]` with `sigmaRatio` from the hardcoded limiting values, so the dynamical formula $\widehat T_m=\sqrt{9\Pi/(20[1-R_qS_q])}$ is never instantiated in Mathematica). |

Dominant pattern: every paper deliverable except the algebraic identity (1a) is at best printed without an assertion; one deliverable (1b) requires user-side reasoning ("π²/2 - 4 > 0") that the script does not exercise; the Mathematica script also hardcodes intermediate limit values that the SymPy script actually computes via `sp.limit`. Front-matter `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 32 | `expect_zero("numerator - exact positive decomposition", num - decomp)` | (1a) algebraic identity | yes |
| A2 | mathematica | 46 | `expectZero["numerator - exact positive decomposition", num - decomp]` | (1a) algebraic identity | yes — but transliteration of A1 |

All other lines are `print` statements. The sympy script lacks any assertion on the three positivity conditions, the two endpoint limits, $R_\infty$, $S_\infty$, $\Sigma_0/\Pi$ ratio, and $\widehat T_m/\sqrt\Pi$ ratio. The mathematica script has the same gap, and additionally hardcodes `sInf=1` (line 62) rather than computing the limit, and never defines `Sigma0`/`That` (i.e. the dynamical formula for traction is absent on the Mathematica side).

## Findings

### F1 — insufficient_verification

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:39-57`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:53-69`

**What's wrong:**
The paper notes (`notes/stages/moving_throat_pde_stage143_equal_normalized_singular_limit.md`) carry five quantitative deliverables beyond the algebraic identity: (i) $\lim_{\Pi\to 0^+}\mathfrak g_\Pi = 2/\pi$, (ii) $\lim_{\Pi\to\infty}\mathfrak g_\Pi = 1$ (the singular-limit statement), (iii) $R_\infty\approx 0.1454544522604201$, (iv) $S_\infty=1$, (v) $\widehat T_m/\sqrt\Pi\to\sqrt{9/(20(1-R_\infty))}\approx 0.7256691307$. None of these are asserted in either script — they are only `print`ed. A regression that broke any of these (e.g. typo in `gPi`, sign flip in `Rq`, swapped factor in `Sigma0`) would still produce output and exit 0. Quoting the script: sympy lines 41-43 are `g0 = sp.limit(...)` / `ginf = sp.limit(...)` / `print("lim_{Pi->0+} g_Pi =", g0)` — no `assert g0 == 2/pi`; mathematica lines 55-58 are the symmetric `FullSimplify[Limit[...]]` + `Print[...]` block with no `expectZero`/`pass`.

Additionally, the notes' load-bearing positivity argument requires three sign checks: $e^\Pi - 1 - \Pi - \Pi^2/2 > 0$ for $\Pi>0$, $\pi^2-2\pi>0$, $\pi^2/2-4>0$. The scripts print the three expressions but never test their positivity, so the inequality $g_\Pi<1$ — the paper's primary claim — is not actually verified by code, only verified-by-inspection-of-printed-output.

**Why this matters:**
The card's body equation is "$\mathfrak g_\Pi<1$ for all finite $\Pi>0$; equal-normalized branch requires $\Pi\to\infty$ and divergent traction." This is the stage's bottom-line claim. The script asserts only the algebraic decomposition identity (an intermediate step), not the inequality nor the limit nor the traction divergence. The audit transcript's "PASS" gives a false sense of having verified the physical claim.

**Required change:**
Add explicit assertions for each printed deliverable. In SymPy (after line 43): add `assert sp.simplify(g0 - 2/pi) == 0` and `assert sp.simplify(ginf - 1) == 0`. After line 47-48: `assert sp.simplify(Rinf - (1-r)**2/(1+r**2)) == 0` and `assert sp.simplify(Sinf - 1) == 0`. After line 56: `assert sp.simplify(that_ratio - sp.sqrt(sp.Rational(9,20)/(1-Rinf))) == 0`. For the positivity check, add three numeric/symbolic assertions: `assert sp.simplify(pi**2 - 2*pi) > 0` (this works because `pi` is `sp.pi`), `assert sp.simplify(pi**2/2 - 4) > 0`, and for the exp-remainder use `sp.series(sp.exp(Pi)-1-Pi-Pi**2/2, Pi, 0, 4)` to confirm the leading term is `Pi**3/6` (positive for `Pi>0`) — assert the coefficient. Mirror all of these in the Mathematica script using `expectZero` and `If[#>0, pass[...], fail[...]]` patterns.

**Verification:**
After applying, both transcripts should show new `PASS:` (mathematica) / no-AssertionError (sympy) lines for each of the six new assertions (two endpoint limits, $R_\infty$, $S_\infty$, $\widehat T_m/\sqrt\Pi$, plus three positivity checks for the decomposition pieces).

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:60-64`

**What's wrong:**
The Mathematica script defines `sQ` on line 60 with the full upstream Stage-244 expression `piM*(((Pi/2)*Tanh[Pi/2]) + piM*(Exp[-piM]*Sech[Pi/2] - 1))/((1 - Exp[-piM])*((Pi/2)^2 - piM^2))`, but then never takes its limit. Line 62 reads `sInf = 1;` — a hardcoded literal. The next line `sigmaRatio = FullSimplify[1/(1 - rInf*sInf)]` and line 64 `tHatRatio = FullSimplify[Sqrt[(9/20)*sigmaRatio]]` are pure algebraic substitutions using this hardcoded `sInf`. There is no derivation of $S_\infty=1$ from the formula for $S_q$. The corresponding SymPy line 48 `Sinf = sp.simplify(sp.limit(Sq, Pi, sp.oo))` actually computes the limit. The Mathematica script also never defines `Sigma0` (the dynamical $\Sigma_0=\Pi/(1-R_q S_q)$) nor `That` (the dynamical $\widehat T_m=\sqrt{9\Pi/(20[1-R_q S_q])}$) — `tHatRatio` is reached by direct substitution into the limiting form, not by taking a limit of the dynamical formula.

**Why this matters:**
This is the Mathematica side's only "independent" check of the $\widehat T_m/\sqrt\Pi$ ratio. With `sInf` hardcoded and `Sigma0`/`That` absent, the Mathematica numeric agreement with SymPy on `lim That/sqrt(Pi) ≈ 0.7256691307` reflects only that both scripts evaluate $\sqrt{9/(20(1-R_\infty))}$ algebraically — it does NOT confirm that the formula $\widehat T_m=\sqrt{9\Pi/(20[1-R_qS_q])}$ actually has that limit. A bug in the assumed limiting form $S_\infty=1$ (e.g. if the upstream Stage-244 closure were modified) would be silently invisible to the Mathematica audit.

**Required change:**
In the Mathematica script, replace `sInf = 1;` (line 62) with `sInf = FullSimplify[Limit[sQ /. piM -> piInf3, piInf3 -> Infinity], Assumptions -> piInf3 > 0];` (introducing a fresh limit-dummy in the same style as `pi0`/`piInf` on lines 55-56). Then explicitly define `sigma0 = piM/(1 - ((gPi - r)^2/(1 + r^2))*sQ);` and `that = Sqrt[(9/20)*sigma0];` after line 60 (i.e. construct the dynamical objects). Replace lines 63-64 with `sigmaRatio = FullSimplify[Limit[sigma0/piM /. piM -> piInf4, piInf4 -> Infinity], Assumptions -> piInf4 > 0]; tHatRatio = FullSimplify[Limit[that/Sqrt[piM] /. piM -> piInf5, piInf5 -> Infinity], Assumptions -> piInf5 > 0];`. This mirrors what SymPy lines 45-57 do.

**Verification:**
After applying, the Mathematica output line for `R_infty`, `S_infty`, `lim Sigma0/Pi`, and `lim That/sqrt(Pi)` should still produce the same numeric values ($\approx 0.14545$, $1$, ..., $\approx 0.72567$) but now via a `Limit` call on `sQ`/`sigma0`/`that` rather than algebraic substitution. The new `sInf` (now a limit, not a literal) should still simplify to `1`.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl:36-47`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py:19-32`

**What's wrong:**
Side-by-side correspondence:

SymPy (lines 19, 24-32):
```
Pi = sp.symbols("Pi", positive=True, real=True)
pi = sp.pi
r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))
num = sp.expand((1-gPi)*(4*Pi**2+pi**2)*(sp.exp(Pi)-1))
decomp = pi**2*(sp.exp(Pi)-1-Pi-Pi**2/2) + Pi*(pi**2-2*pi) + Pi**2*(pi**2/2-4)
expect_zero("numerator - exact positive decomposition", num - decomp)
```

Mathematica (lines 36-46):
```
$Assumptions = Element[piM, Reals] && piM > 0;
r = Sqrt[4107 - 100*Pi^2]/(10*Pi);
gPi = 2*piM*(2*piM*Exp[piM] + Pi)/((4*piM^2 + Pi^2)*(Exp[piM] - 1));
num = Expand[(1 - gPi)*(4*piM^2 + Pi^2)*(Exp[piM] - 1)];
decomp = Pi^2*(Exp[piM] - 1 - piM - piM^2/2) + piM*(Pi^2 - 2*Pi) + piM^2*(Pi^2/2 - 4);
expectZero["numerator - exact positive decomposition", num - decomp];
```

This is a token-for-token transliteration: same variable names (`r`, `gPi`, `num`, `decomp`), same intermediate algebra, same order of operations, identical pre-baked decomposition. The Mathematica script does not "re-derive" the decomposition — it imports the same closed form the SymPy script imports. An independent Mathematica audit would, for example: (a) start from $\mathfrak g_\Pi$, compute `1 - gPi` and ask `FullSimplify` to surface a positive-piece decomposition via `Apart`/`Collect`, or (b) prove positivity differently — e.g. by showing the Taylor coefficients of the numerator in $\Pi$ are all positive — using a mechanism distinct from explicit subtraction against a hand-rolled `decomp`.

**Why this matters:**
The second-engine policy is meant to catch errors that one CAS makes but the other doesn't (different simplification heuristics, different branch handling, different symbolic algorithms). When the Mathematica script just rewrites the SymPy algebra in Mathematica syntax, the second engine is only checking that two CASes agree on a trivial subtraction — not that two independent derivation routes converge on the same answer.

**Required change:**
Re-derive the decomposition in Mathematica from a structurally different path. Two acceptable patterns:
(a) Start from `num = (1 - gPi)*(4*piM^2 + Pi^2)*(Exp[piM] - 1)` and use `Series[num, {piM, 0, 6}]` to extract the leading Taylor coefficients in `piM`; then assert each coefficient is non-negative (the coefficient of `piM^k` for k=0,1,2 vanishes; the coefficient of `piM^3` is positive; etc.). This proves positivity in a small-`piM` neighborhood by a different mechanism than `decomp` subtraction.
(b) Define the three "positive pieces" individually and prove the decomposition using `Reduce[num == p1 + p2 + p3, {piM}, Reals]` returning `True` over `piM > 0`; this is structurally different from subtracting a hand-built `decomp`.

If preserving the existing `expectZero["numerator - exact positive decomposition", num - decomp]` line is required for the algebraic identity check, then also add one of (a) or (b) as a second, independent assertion of positivity that does NOT mirror the SymPy choreography. The simplest is to add `Reduce[num > 0, piM, Reals]` and assert the result is `piM > 0`, which independently certifies the inequality without depending on the SymPy-style decomposition.

**Verification:**
After applying, the Mathematica script's `num`-positivity check should not share a one-to-one line correspondence with the SymPy block. The new check should produce output naming a distinct method (e.g. `Series` expansion coefficients, `Reduce` over reals).

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. See F3 for the side-by-side. Beyond the structural mirror, the Mathematica script:
1. Uses the exact same hand-rolled `decomp = Pi^2*(Exp[piM] - 1 - piM - piM^2/2) + piM*(Pi^2 - 2*Pi) + piM^2*(Pi^2/2 - 4)` form as SymPy, with the same operator order and grouping.
2. Mirrors the SymPy variable choreography (`r`, `gPi`, `num`, `decomp`).
3. Differs only in (a) renaming the mouth-bias variable `Pi`→`piM` (necessary because Mathematica's `Pi` is the built-in $\pi$ and would clash) and (b) the unrelated hardcoding of `sInf = 1` on the limits side (which is a finding in itself, F2). No independent derivation path is exercised.

## Engine cross-check

Both engines produce identical numerics:
- SymPy: `R_infty = 0.145454452260420126101421595368`, `lim That/sqrt(Pi) = 0.725669130700713219781041125011`.
- Mathematica: `R_infty = 0.1454544522604201261014215953679127185274370895574784984548`, `lim That/sqrt(Pi) = 0.7256691307007132197810411250114145098613931689620003558428`.

Agreement to 30 digits on both reported quantities. No `engine_disagreement`. (Symbolic forms for `R_infty` differ: SymPy gives `(-sqrt(4107-100π²)+10π)²/4107`, Mathematica gives `1 - 20π·sqrt(4107-100π²)/4107`. Direct expansion confirms these are equal: $(10\pi-\sqrt{4107-100\pi^2})^2/4107 = (100\pi^2 - 20\pi\sqrt{4107-100\pi^2} + 4107 - 100\pi^2)/4107 = 1 - 20\pi\sqrt{4107-100\pi^2}/4107$. Same expression in different simplifier output form.)

## Verdict justification

The single working assertion in both scripts is non-tautological (the numerator/decomposition identity is a substantive algebraic check) and the engines agree on it. The numerics for $R_\infty$, $S_\infty$, and $\widehat T_m/\sqrt\Pi$ match the notes' stated values to 30 digits. So the math, when fully derived by hand from the printed outputs, supports the paper claims. However: (i) the scripts assert almost nothing of what the paper claims — the inequality $g_\Pi<1$, the limit $g_\Pi\to 1$, the values $R_\infty\approx 0.1454$ and $\widehat T_m/\sqrt\Pi\approx 0.7257$ are all printed but not gated, so a future regression breaking any of these would still PASS; (ii) the Mathematica script hardcodes `sInf = 1` and skips defining `Sigma0`/`That`, so its agreement with SymPy on the traction limit is partly artifactual; (iii) the `.wl` is a near-transliteration of the `.py`, violating the second-engine policy. Verdict is `findings`. No stop-cold — the underlying math is correct and consistent with the paper; the work is in tightening assertions to actually exercise the paper's stated deliverables.

## Self-test notes

I checked the proposed F1 patches mentally: `sp.simplify(pi**2 - 2*pi) > 0` returns a sympy `BooleanTrue` because `pi` is `sp.pi` (concrete numeric, ≈3.14159), so this evaluates and asserts cleanly. `sp.simplify(pi**2/2 - 4) > 0` likewise (`pi**2/2 ≈ 4.93`). For the exp-remainder, the Taylor series `Pi**3/6 + O(Pi**4)` confirms positivity for small `Pi>0`; combined with monotonicity (derivative is `e^Pi - 1 - Pi`, which is itself positive for `Pi>0` by an analogous Taylor argument), positivity holds globally. For F2, I confirmed `sQ` has the right $\Pi\to\infty$ behavior by hand: numerator $\sim \Pi\cdot(-\Pi)\to-\Pi^2$, denominator $\sim 1\cdot(-\Pi^2)\to-\Pi^2$, ratio $\to 1$, matching the paper's $S_\infty=1$. No paper-card mutation is implied by any of F1-F3 — all fixes are script-side. (One cosmetic note not raised as a finding: both scripts' top banners read "STAGE 126" rather than "STAGE 143" — the script was lifted from stage 126's slot. This is purely a print-string artifact, does not affect any assertion, and is outside the 10 finding categories, so I did not file it as a formal finding.)
