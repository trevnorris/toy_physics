---
unit_id: 050
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md]
  paper_appendix: present
---

# Audit unit 050 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_050.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 78)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.txt`

## What the paper claims

Stage 050 compares the physical coherent support ratio `zeta_n^(phys)` (from Stage 049) against the required threshold `zeta_req` (from Stage 048) and orders the support tower. `\stagefield{Output}` states: "Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success}, higher-harmonic exclusion/softness thresholds, and the enhancement ceiling \eqref{eq:app-stage050-Sn-max}." Concretely the deliverables are: (D1) the lowest symmetric twin has `zeta_0 = 1` so its enhancement is `S_0 = S(1;eps) = 2` exactly, independent of `eps`; (D2) lowest-twin success iff `zeta_req <= 1 <=> S_req <= 2`; (D3) the higher-harmonic immediate exclusion test `zeta_req > 1/(2n+1)^2 ⟹ n`th twin lane cannot reach threshold; (D4) the softness ceiling `x <= x_max(n;zeta_req) = [1/((2n+1)^2 zeta_req) - 1]/[n(n+1)]`; and (D5) the higher-harmonic enhancement ceiling `S_n^(twin)(x;eps) < S_n^(max)(eps) := 1 + (1-eps)/((2n+1)^2 - eps)`. The notes add the explicit closed forms `S_1^(max) = (10-2eps)/(9-eps)` and `S_0^(twin)=2`, and the no-go tower (`n=1` impossible if `zeta_req>1/9`, etc.). No load-bearing free numeric constant is introduced at this stage beyond these closed forms.

## What the script claims to verify

Both scripts assert (via `expect_zero`/`expectZero`, hard `raise`/`Exit[1]` on nonzero residual): the doubling identity `S(1;eps)-2 = 0` and the anchored `zeta_0^(twin)-1 = 0` (SymPy imports `twin_support_ratio(0,x)` from the Stage-049 module; Mathematica inlines the same closed form with `n->0`); the factored equivalence `zeta_req - 1 = (1-eps)(S_req-2)/(1+eps(S_req-2))` (so `zeta_req<=1 <=> S_req<=2`); the monotonicity-in-`x` derivative form of `zeta_n^(twin)`; the `x_max` closed form obtained by `solve`/`Solve` of `zeta_n = zeta_req` (SymPy checks `x_eq` against the boxed closed form, Mathematica back-substitutes `zeta_n(x_max)-zeta_req=0` and checks the implicit relation); the admissibility numerator reducing the non-negativity of `x_max` to `(2n+1)^2 zeta_req <= 1`; the ceiling saturation `S_n^(twin)(x=0)-S_n^(max)=0` plus the factored sign-definite form of `S_n^(max)-S_n^(twin)`; and prints the explicit `S_1^(max)`, `S_2^(max)`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: `zeta_0=1`, `S_0=S(1;eps)=2` | sympy L48-51 / wl L44-48 (`zeta_0^(twin)-1=0`, `S(1;eps)-2=0`) | match |
| D2: `zeta_req<=1 <=> S_req<=2` | sympy L52-58 / wl L49-54 (factored `zeta_req-1` equivalence) | match |
| D3: exclusion `zeta_req>1/(2n+1)^2` | sympy L81-90 / wl L75-83 (admissibility numerator ⟹ `(2n+1)^2 zeta_req<=1`) | match |
| D4: `x_max=[1/((2n+1)^2 zeta_req)-1]/[n(n+1)]` | sympy L74-79 / wl L57-73 (solve + closed-form/back-sub) | match |
| D5: `S_n^(twin)<S_n^(max):=1+(1-eps)/((2n+1)^2-eps)` | sympy L92-116 / wl L85-109 (saturation + sign-definite factored diff) | match |
| notes: `S_1^(max)=(10-2eps)/(9-eps)`, `S_2^(max)` | sympy L115-116 / wl L108-109 (printed) | match (forms agree: `2(eps-5)/(eps-9)=(10-2eps)/(9-eps)`) |

`paper_alignment: aligned`. Every boxed/Output deliverable has a non-tautological script-side check in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 47-50 | `twin_support_ratio(0,x) - 1 == 0` | D1 | yes (imports 049 helper, n=0 reduces to 1) |
| A2 | sympy | 51 | `simplify(S(1;eps)-2) == 0` | D1 | yes |
| A3 | sympy | 53-58 | `factor(zeta_req-1) - (1-eps)(Sreq-2)/(...) == 0` | D2 | yes |
| A4 | sympy | 67-72 | `diff(zeta_n,x) - target == 0` | D4 (monotonicity) | yes |
| A5 | sympy | 74-79 | `solve(zeta_n=zeta_req,x) - boxed x_max == 0` | D4 | yes |
| A6 | sympy | 84-89 | `admissibility numerator == 0` ⟹ `(2n+1)^2 zeta_req<=1` | D3 | yes |
| A7 | sympy | 97-98 | `S_n(x=0) - S_n_max == 0` | D5 (saturation) | yes |
| A8 | sympy | 100-111 | `(S_n_max - S_n) - factored target == 0` | D5 (sign-definite) | yes |
| A9 | mathematica | 44-47 | `(zeta_n /. n->0) - 1 == 0` | D1 | yes |
| A10 | mathematica | 48 | `(S /. zeta->1) - 2 == 0` | D1 | yes |
| A11 | mathematica | 49-54 | `(zeta_req-1) - (1-eps)(sReq-2)/(...) == 0` | D2 | yes |
| A12 | mathematica | 63-67 | `dZetaNdx (2n+1)^2(1+n(n+1)x)^2 + n(n+1) == 0` | D4 | yes |
| A13 | mathematica | 69 | `(zeta_n /. x->xEq) - zeta_req == 0` | D4 | yes |
| A14 | mathematica | 70-73 | `(2n+1)^2 zeta_req (1+n(n+1)xEq) - 1 == 0` | D3/D4 | yes |
| A15 | mathematica | 75-83 | admissibility numerator `== 0` | D3 | yes |
| A16 | mathematica | 89 | `(sN /. x->0) - sNMax == 0` | D5 | yes |
| A17 | mathematica | 92-99 | `Numerator(sNMax-sN) - (1-eps)(2n+1)^2 n(n+1)x == 0` | D5 | yes |
| A18 | mathematica | 100-107 | `Denominator(sNMax-sN) - ((2n+1)^2-eps)(...) == 0` | D5 | yes |

No assertion is tautological-by-construction: each `expect_zero`/`expectZero` compares an independently-`solve`d/`diff`'d/`simplify`d quantity against a hand-written boxed target, so a wrong boxed target (or a wrong `solve` branch) would leave a nonzero residual and abort. None of the checks define `x = expr` then assert `x == expr`.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.txt` (mtime 2026-05-26 02:58:23)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.txt` (mtime 2026-05-26 02:58:29)

**What's wrong:**
Both saved outputs predate both scripts (`.py` and `.wl` mtime 2026-06-03 15:59:11). The content confirms staleness: the committed SymPy transcript header reads `STAGE 33 — PHYSICAL ZETA VS ZETA_REQ` (line 3) and closes `All Stage-33 symbolic checks passed.` (line 42), and the Mathematica transcript reads `STAGE 033 — PHYSICAL ZETA VS ZETA_REQ` (line 3), whereas the current scripts emit `STAGE 50`/`STAGE 050` banners (`.py:34`, `.wl:32`). The math content in the transcripts still matches the current scripts (all residuals 0, all PASS), so the staleness is banner-label only — but the freshness gate should refresh both.

**Why this matters:**
The committed transcript is supposed to be the authoritative record of current script output; a stale banner could mislead a future reader into thinking the output reflects an older script revision.

**Required change:**
Re-run both scripts and recommit the refreshed transcripts (orchestrator does this on the freshness gate). No source change required for this finding.

**Verification:**
After re-run, `scripts/output/...sympy_audit.txt` line 3 reads `STAGE 50 — PHYSICAL ZETA VS ZETA_REQ` and the close line reads the Stage-50 wording; `mathematica/output/...mathematica_audit.txt` line ~3 reads `STAGE 050 — ...`. mtimes newer than the scripts.

### F2 — stale_output (numbering self-labels)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:2` — docstring `"""Stage 33 SymPy audit.`
- `.../stage050_..._sympy_audit.py:118` — `print("\nAll Stage-33 symbolic checks passed.")`
- `.../stage050_..._sympy_audit.py:61-62` — comment `# Imported from Stage 32's explicit D/N overlap extraction` and `# Stage 32's`

**What's wrong:**
Stale `Stage NNN` self-labels in the SymPy script. The docstring (`.py:2`) and closing print (`.py:118`) say "Stage 33"/"Stage-33"; the banner (`.py:34`) correctly says "STAGE 50". The import-source comment (`.py:61-62`) credits "Stage 32's explicit D/N overlap extraction," but the actual import is `from moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit import twin_support_ratio` (`.py:17`) — the helper is Stage 049, not "Stage 32". These are the known low-severity stale-self-label class (incomplete renumber from the EM-extension realignment); the math/import target itself is correct.

**Why this matters:**
Cosmetic/traceability only — no math impact — but the docstring and the import-source comment misattribute the helper's origin stage, which can confuse provenance tracing.

**Required change:**
Per Reading-2 in-loop policy: docstring `.py:2` → "Stage 50 SymPy audit."; close line `.py:118` → "All Stage-50 symbolic checks passed."; comment `.py:61-62` → reference "Stage 049's explicit D/N overlap extraction" / "Stage 049's". (Orchestrator may fold label scope per its separate numbering pass; the import statement `.py:17` is already correct and must not change.)

**Verification:**
`.py:2`, `.py:118`, `.py:61-62` no longer contain `Stage 33`/`Stage 32`; refreshed SymPy transcript close line matches.

## Independent-derivation check (Mathematica)

The Mathematica script is NOT a line-by-line transliteration; it re-derives independently and even differs in method on the load-bearing step. Three corresponding sections:

- **`x_max` derivation.** SymPy: `x_eq = sp.solve(sp.Eq(zeta_n, zeta_req), x)[0]` then asserts `x_eq -` boxed closed form `== 0` (`.py:74-79`). Mathematica also `Solve`s (`.wl:57`) but then verifies the solution *differently* — it back-substitutes `(zetaN /. x -> xEq) - zetaReq == 0` (`.wl:69`) and checks the implicit relation `(2n+1)^2 zetaReq (1 + n(n+1) xEq) - 1 == 0` (`.wl:70-73`), rather than comparing against the explicit boxed form. Independent verification route.
- **Ceiling sign-definiteness.** SymPy compares the whole difference `S_n_max - S_n` against one hand-written factored target (`.py:100-111`). Mathematica instead splits `Numerator[Together[...]]` and `Denominator[Together[...]]` and checks each against its own target (`.wl:92-107`). Different decomposition of the same claim.
- **`zeta_0` anchor.** SymPy imports the Stage-049 helper `twin_support_ratio(0,x)` (`.py:17,49`); Mathematica inlines the closed form `1/((2n+1)^2(1+x n(n+1))) /. n->0` (`.wl:46`). The engines reach the doubling anchor by different means (cross-module import vs. inline literal), which is the opposite of transliteration.

No `mathematica_transliteration` finding.

## Engine cross-check

Both transcripts show identical final symbolic results and all checks PASS / residual 0:
- `zeta_req`: SymPy `(S_req-1)/(eps*(S_req-2)+1)`; Mathematica `(-1+sReq)/(1+eps*(-2+sReq))` — same.
- `zeta_n^(twin)`: SymPy `1/((2n+1)^2*(n*x*(n+1)+1))`; Mathematica `1/((1+2n)^2*(1+n*(1+n)*x))` — same.
- `S_1^(max)`: SymPy `2*(eps-5)/(eps-9)`; Mathematica `1+(-1+eps)/(-9+eps)` — algebraically equal (`= (10-2eps)/(9-eps)`).
- `S_2^(max)`: SymPy `2*(eps-13)/(eps-25)`; Mathematica `1+(-1+eps)/(-25+eps)` — equal.
- All `expect_zero`/`expectZero` residuals are 0 in both. Engines agree.

## Verdict justification

Attacks attempted and survived: (a) tested every `expect_zero` for tautology — each compares an independently `solve`/`diff`/`simplify`-derived object against a hand-written boxed target, so a wrong boxed form would leave a nonzero residual (not tautological); (b) checked the `n->0` substitution in the Mathematica script under its `$Assumptions` `n>=1` — the substitution is a pure numeric reduction `1/(1*(1+0))-1=0` performed independent of the assumption set, no domain contradiction; (c) checked `sp.solve(...)[0]` branch selection — the single rational equation `zeta_n=zeta_req` has one solution, and SymPy's branch is cross-checked against the boxed closed form (and Mathematica back-substitutes), so no missing-branch; (d) checked positivity assumptions (`eps,Sreq,x` positive, `0<eps<1` only in `.wl`) against the paper's `0<eps<1, S_req>1, x>0` — consistent; (e) confirmed every paper Output deliverable D1–D5 plus the notes' explicit `S_1^(max)/S_2^(max)` forms maps to a matched script check in both engines. The math is fully aligned with the paper. The only findings are the low-severity stale-output (banner predates the scripts) and the stale numbering self-labels in the SymPy docstring/close-line/import-comment; both are non-blocking, material_change false. Verdict: findings.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit against the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `zeta_req = (S_req-1)/(1+eps(S_req-2))` | py:41 / wl:39; sympy.txt:9, math.txt:5 | notes:65; (.tex carries it via Stage-048 input) | MATCH |
| `S(zeta;eps)=1+zeta(1-eps)/(1-eps zeta)` | py:42 / wl:40; sympy.txt:10, math.txt:6 | notes:71 | MATCH |
| `zeta_0^(twin)=1` | py:49 / wl:46; sympy.txt:11 (residual 0) | .tex:13 ("`zeta_0=1`"); notes:24,85 | MATCH |
| `S_0 = S(1;eps) = 2` | py:51 / wl:48; sympy.txt:12, math.txt:9 | .tex:17 (boxed `S_0=2`); notes:29,92,183 | MATCH |
| `zeta_req<=1 <=> S_req<=2` (factored equiv.) | py:53-58 / wl:49-54; sympy.txt:13-14 | .tex:23-25 (boxed); notes:99-101 | MATCH |
| `zeta_n^(twin)=1/((2n+1)^2(1+x n(n+1)))` | py:63 (import) / wl:56; sympy.txt:19, math.txt:14 | notes:17-18,114-115 | MATCH |
| `x_max=[1/((2n+1)^2 zeta_req)-1]/[n(n+1)]` | py:74-79 / wl:57-73; sympy.txt:21-22 | .tex:40-41 (boxed); notes:51,154-155 | MATCH |
| exclusion `zeta_req>1/(2n+1)^2` ⟹ n-lane fails | py:84-90 / wl:75-83; sympy.txt:27-29, math.txt:22-24 | .tex:31-33 (boxed); notes:121-134 | MATCH |
| `S_n^(max)=1+(1-eps)/((2n+1)^2-eps)` | py:94 / wl:86; sympy.txt:35, math.txt:26 | .tex:47-48 (boxed); notes:172-173 | MATCH |
| `S_1^(max)=(10-2eps)/(9-eps)` [`2(eps-5)/(eps-9)`] | py:115 / wl:108; sympy.txt:39, math.txt:34 | notes:177 | MATCH |
| `S_2^(max)=2(eps-13)/(eps-25)` | py:116 / wl:109; sympy.txt:40, math.txt:35 | (.tex omits; notes give the general `S_n^(max)` form) | MATCH (general form covers it; specific n=2 is intermediate) |

INTERNAL / scaffolding (no prose expected): the derivative residual `d zeta_n/dx` form (py:67-72/wl:63-67), the `admissibility numerator residual` (py:84/wl:75), the saturation residual `S_n(x=0)-S_n_max` (py:97/wl:89), the ceiling factored/numerator/denominator residuals (py:100-111/wl:90-107) — all are zero-residual verification quantities, not reported deliverables.

reconciliation: complete; 11 values checked, 0 misaligned

## Self-test notes

Checked the three standard traps. (1) Variable independence: the only derivative is `sp.diff(zeta_n, x)`/`D[zetaN,x]` and `zeta_n` genuinely depends on `x` (denominator `1+x n(n+1)`), so the derivative is non-trivially nonzero and the monotonicity assertion is substantive — not an identically-zero derivative trap. (2) No unbounded-domain integrals in this stage, so the parity trap does not apply. (3) Trivial-case pre-check: `zeta_0=1`, `S(1;eps)=2`, and `S_n^(max)|_{x=0}=S_n^(twin)|_{x=0}` all reduce to literal 0 residuals as the transcripts show; the `solve`-derived `x_max` back-substitutes to `zeta_req`. No directive with new assertions is warranted (both findings are stale_output, non-script-math); per the prompt, stale_output is informational and does not by itself require a Codex directive, so no directive is written.
