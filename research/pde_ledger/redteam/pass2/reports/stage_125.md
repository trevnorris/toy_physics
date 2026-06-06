---
unit_id: 125
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: ["notes/stages/moving_throat_pde_stage125_positive_source_theorem.md"]
  paper_appendix: present
---

# Audit unit 125 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_125.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage125_positive_source_theorem.md` (only file matching `moving_throat_pde_stage125_*.md`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows 29, 86 reference Stage 125)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.txt`

## What the paper claims

Stage 125 replaces the still-free normalized mouth-coupling ratio by an explicit localized positive-source law, reducing the mouth-bias factor to the first cosine moment of a nonnegative normalized axial source profile: `g[σ] = ∫₀ᴸ σ(z) cos(πz/2L) dz`. The card's bottom-line result (quoted verbatim) is: "Positive normalized sources satisfy \(0\le\mathfrak g[\sigma]\le1\), ruling out the upper compensated branch." The notes add the supporting detail: because `0 ≤ cos(πz/2L) ≤ 1` on `[0,L]`, every positive normalized source has its cosine moment in `[0,1]` (boxed `0 ≤ g[σ] ≤ 1`); and on the explicit Family-1 branch the two carried compensation roots are `g₋^F1 ≈ 0.758035078944663` (in `(0,1)`, admissible) and `g₊^F1 ≈ 2.79795199200529` (`>1`, excluded). Deliverables: (D1) the cosine-moment reduction kernel; (D2) the kernel bound `0 ≤ cos(πx/2) ≤ 1` on `[0,1]`; (D3) the exact positivity theorem `0 ≤ g[σ] ≤ 1`; (D4) the two Family-1 branch values with `g₊>1` and `0<g₋<1` selecting the lower branch.

## What the script claims to verify

Both scripts verify the same five things. (1) The half-wave kernel reduces to `cos(πz/2L)` and on the normalized variable `x∈[0,1]` has global min 0 and max 1 (engine-native `minimum`/`maximum` in SymPy, `MinValue`/`MaxValue` in Mathematica). (2) The Family-1 radical `r = √(12(37/20)²/π² − 1)` equals `√(4107−100π²)/(10π)`. (3) The two compensation roots of the balance relation `1 + r² − 4(g−r)² = 0` are the closed forms `(2R∓37√3)/(20π)`, with `g₋∈(0,1)` and `g₊>1`. (4) A one-parameter positive normalized source family `σ_a(z)=(a+1)(z/L)^a/L` integrates to 1, has moment `2/π` at `a=0` (uniform), and at the peaked proxy `a=100` satisfies the genuine paper bound `0 ≤ g ≤ 1` (plus a retained `<1/20` smallness fact). (5) Both uniform-endpoint moments satisfy `0 ≤ g ≤ 1`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 cosine-moment kernel `cos(πz/2L)` | sympy L36-39 / wl L38-42 print + use of kernel | match |
| D2 kernel bound `0 ≤ cos ≤ 1` on `[0,1]` | sympy L41-44 min/max; wl L44-47 MinValue/MaxValue | match |
| D3 positivity theorem `0 ≤ g[σ] ≤ 1` | sympy L88-95 / wl L96-100 range checks on σ_a family at endpoints | match |
| D4 branch values `g₋≈0.758`, `g₊≈2.798`; `g₊>1`, `0<g₋<1` | sympy L51-66 / wl L57-77 closed forms + balance relation + inequalities | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43 | `expect_zero(kernel min on [0,1])` | D2 | yes |
| A2 | sympy | 44 | `expect_zero(kernel max − 1)` | D2 | yes |
| A3 | sympy | 54 | `expect_zero(r − R/(10π))` | D4 (radical) | yes |
| A4 | sympy | 55 | `expect_zero(1+r²−4(g₋−r)²)` | D4 (lower root) | yes |
| A5 | sympy | 56 | `expect_zero(1+r²−4(g₊−r)²)` | D4 (upper root) | yes |
| A6 | sympy | 64-66 | `g₋>0`, `g₋<1`, `g₊>1` | D4 (branch selection) | yes |
| A7 | sympy | 74 | `expect_zero(norm_a − 1)` | D3 (family normalization) | yes |
| A8 | sympy | 79 | `expect_zero(g[uniform] − 2/π)` | D3 (moment value) | yes |
| A9 | sympy | 88-89 | `g[a=100] ≥ 0`, `≤ 1` | D3 (genuine bound at peaked end) | yes |
| A10 | sympy | 92 | `g[a=100] < 1/20` | D3 (trend, retained smallness) | yes (supplementary) |
| A11 | sympy | 94-95 | `g[uniform] ≥ 0`, `≤ 1` | D3 | yes |
| B1 | math | 46-47 | `expectZero(kernel min / max−1)` | D2 | yes |
| B2 | math | 57-58 | `Solve` balance eq → branch roots | D4 (derives roots) | yes |
| B3 | math | 62 | `expectZero(r − rrad/(10π))` | D4 | yes |
| B4 | math | 63-66 | `expectZero(g₋/g₊ − closed forms)` | D4 | yes |
| B5 | math | 75-77 | `g₋>0`, `g₋<1`, `g₊>1` | D4 | yes |
| B6 | math | 86 | `expectZero(normA − 1)` | D3 | yes |
| B7 | math | 92 | `expectZero(g[uniform] − 2/π)` | D3 | yes |
| B8 | math | 93 | `expectZero(Limit g_a, a→∞)` | D3 (peaked-end → 0) | yes |
| B9 | math | 96-97 | `g[a=100] ≥ 0`, `≤ 1` | D3 (genuine bound) | yes |
| B10 | math | 99-100 | `g[uniform] ≥ 0`, `≤ 1` | D3 | yes |

All rows trace to a paper deliverable; none orphaned.

## Findings

None.

The historical strengthening is confirmed fresh: the prior weak `abs(g_a(100)) < 0.05` smallness sample-point check has been replaced by the genuine paper bound. SymPy lines 88-89 assert `g[peaked@L proxy a=100] >= 0` and `<= 1` **not** wrapped in `abs()` (the inline comment at L86-87 explicitly notes a negative moment from a sign error would now FAIL the lower bound), with the `<1/20` trend fact retained as a supplementary check (L90-92), not the load-bearing one. Mathematica lines 96-97 mirror the genuine bound, and additionally adds an honest symbolic `Limit[g_a, a→∞] = 0` check (L93) that SymPy could not take. The bound `0 ≤ g ≤ 1` — the paper's actual claim — is what is exercised, on a concrete positive-source family spanning both endpoints (uniform `a=0`, peaked `a=100`), plus the kernel-range argument that underwrites the general theorem.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent re-derivation**, not a transliteration, on the load-bearing steps. Two concrete divergences:

1. **Branch roots.** SymPy posits the closed forms and verifies them against the balance relation:
   - sympy L51-56: `gminus = sp.simplify((2*R - 37*sp.sqrt(3))/(20*sp.pi))` … `expect_zero("lower branch balance relation", 1 + r**2 - 4*(gminus - r)**2)`.
   Mathematica instead **solves** the balance equation for the roots and then checks them against the closed forms:
   - wl L57-66: `branchSolutions = Solve[1 + r^2 - 4*(gSym - r)^2 == 0, gSym]` … `expectZero["g_- matches closed form (2*rrad - 37*Sqrt[3])/(20*Pi)", gminus - (2*rrad - 37*Sqrt[3])/(20*Pi)]`.
   This is opposite choreography (derive-then-match vs posit-then-verify), the hallmark of an independent route.

2. **Parametric moment closed form.** The two CAS engines return genuinely different symbolic forms for the same integral: SymPy emits a hypergeometric `hyper((σ/2+1/2,), (1/2, σ/2+3/2), -π²/16)` (output L25); Mathematica emits an incomplete-Gamma form `Gamma[1+aSym, ±iπ/2]…` (output L34). Both numerically coincide at `a=100` to 15 digits (`0.0153964175481027`). A line-by-line port would have produced the same intermediate form.

3. **Kernel bound** uses engine-native optimizers (`sp.calculus.util.minimum/maximum` vs `MinValue/MaxValue`), and Mathematica adds a symbolic `Limit` peaked-end check absent in SymPy.

No `mathematica_transliteration` finding.

## Engine cross-check

The engines agree at every level they claim:

| Quantity | SymPy | Mathematica |
|---|---|---|
| kernel min / max−1 | `0 / 0` | `0 / 0` (PASS) |
| `r − R/(10π)` | `0` | `0` (PASS) |
| g₋ (numeric) | `0.75803507894466282692` | `0.75803507894466282691968…` |
| g₊ (numeric) | `2.7979519920052934101` | `2.79795199200529341011158…` |
| `g[uniform] − 2/π` | `0` | `0` (PASS) |
| `g_a` at `a=100` | `0.0153964175481027` | `0.01539641754810270040678…` |
| branch inequalities | all True | all PASS |

`engines_agree: true`.

## Verdict justification

`clean`. Attacks attempted that failed: (a) tried to make the SymPy balance-relation check tautological — it is not; perturbing the closed form `37√3 → 36√3` yields a nonzero residual `219/(100π²)`, so the check genuinely tests the closed-form roots against the independently-stated balance equation. (b) Checked whether the strengthened bound was merely a smallness sample-point check in disguise — it is the genuine `0 ≤ g ≤ 1` (un-`abs`'d), exercised at both endpoints of a concrete positive-source family, with the kernel-range argument backing the general theorem. (c) Checked symbol domains: `z, L` positive/real, radicand `4107 − 100π² ≈ 3120 > 0` so `r, R` are genuinely real (no hidden complex branch); `σ_param/aSym ≥ 0` matches the family's nonnegativity premise. (d) Checked engine independence: the two CAS engines derive the branch roots by opposite routes and emit different symbolic moment forms, ruling out transliteration. I read the card, the single notes file, and the appendix rows (29, 86); the script's verified claim (`0 ≤ g[σ] ≤ 1`, `g₊>1` excluded, `0<g₋<1` lower branch unique) matches the paper's `Output` and the notes' boxed results exactly. Outputs are fresh (both `.txt` mtime 15:26:58 > both script mtimes 14:31/14:35), and reflect all-PASS / True.

## Self-test notes

Checked: (1) Variable independence — no `sp.diff`/`D[...]` derivatives in this stage, so no identically-zero-derivative trap. (2) Symmetry/parity — the `g_a` integral is over the bounded `[0,L]` of a strictly positive integrand `(a+1)(z/L)^a cos(πz/2L)/L`, so positivity of the moment is structural, consistent with the `≥0` assertions; the kernel `cos` is positive on the half-period `[0,L]`, no sign cancellation. (3) Trivial-case pre-check — uniform `a=0` gives `2/π ≈ 0.6366 ∈ (0,1)`, peaked `a=100` gives `0.01540 ∈ (0,1)`, both inside the bound as asserted; branch numerics `0.758` and `2.798` reproduce the paper values. No errors uncovered; no directive written.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| cosine-moment kernel `cos(πz/2L)` | py L36/output L6; wl L38/output L6 | notes L24-39 (boxed `g[σ]=∫σ cos(πz/2L)dz`); .tex L13 (block narrative) | MATCH |
| kernel range `[0,1]` on `[0,1]` (min 0, max 1) | py L43-44/output L7-8; wl L46-47/output L8-10 | notes L46-49 (`0≤cos(πz/2L)≤1`) | MATCH |
| positivity bound `0 ≤ g[σ] ≤ 1` | py L88-95; wl L96-100 | .tex L16 (`0\le\mathfrak g[\sigma]\le1`); notes L53 (boxed) | MATCH |
| `g₋^F1 = (2R−37√3)/(20π) ≈ 0.758035078944663` | py L61/output L19; wl L72/output L24 | notes L61 (`≈0.758035078944663`) | MATCH |
| `g₊^F1 = (2R+37√3)/(20π) ≈ 2.79795199200529` | py L62/output L20; wl L73/output L25 | notes L63 (`≈2.79795199200529`) | MATCH |
| uniform moment `g[a=0] = 2/π ≈ 0.6366` | py L79/output L26; wl L92/output L35 | notes L41 ("first cosine moment"; uniform value implied by the `g[σ]` formula) | MATCH |

Notes on guards: the `.tex` card is deliberately terse and carries only the load-bearing bound `0≤g≤1` (L16) plus the qualitative branch-exclusion statement; the branch numerics and the `2/π` uniform value live correctly in the `.md` notes, which is the natural carrier — these are MATCHes per the augmentation guard (a value present in the notes is a MATCH even if the card omits it). The `.txt` outputs are fresh, so the reconciliation is grounded on both source + saved output.

INTERNAL (scaffolding, no finding expected in prose): `r = √(4107−100π²)/(10π) ≈ 1.778` (intermediate radical feeding the branch roots); `g_a` parametric closed form (`hyper(...)` / incomplete-Gamma — engine-specific intermediate, not a deliverable); `g[a=100] ≈ 0.0153964…` (peaked-endpoint proxy probe driving the bound check); `g[a=100] < 1/20` retained smallness fact (supplementary trend, not a stage deliverable); `Limit[g_a, a→∞] = 0` (Mathematica-only honest limit check); pass/True flags and residual-zero outputs.
