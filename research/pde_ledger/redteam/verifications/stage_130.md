---
unit_id: 130
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 130

> Supersedes the prior-pass verdict (which described the now-removed 6-point
> monotonicity sweep). This pass verifies the STRENGTHENED F1 covariance
> certificate that replaced that sweep.

## Per-finding outcomes

### F1 — insufficient_verification (strict monotonicity / Π_* uniqueness "proven" by a 6-point sweep)

**Classification:** resolved

**What changed:**
The finite six-point sweep was REPLACED in BOTH engines by the FKG/Chebyshev
symmetrized-covariance global certificate plus the uniqueness bracket.

- SymPy `scripts/moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:38-99`:
  the old `for val in (1/10, 1/2, 1, 15088/10000, 3, 10)` loop (diff lines 84-88,
  all `-` removals) is GONE. New content: density normalization guard
  (`norm_p - 1 == 0`, :46-48); symmetrized double-integral identity
  `Cov_double = ½∬(f1-f2)(z1-z2)p1 p2` asserted `== Cov := Efz - gPi*Ez`
  (:53-65); `f'(z)` closed-form certificate
  `gprime - (π/2L)sin(πz/2L) == 0` (:73-77); consistency guard
  `dgPi + Cov/L == 0` (:85-86); uniqueness bracket `2/π < g_minus < 1`
  (:93-99).
- Mathematica `mathematica/moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:61-108`:
  the old `Module[{vals = {1/10,…}, dv}, Do[… dg/dpiM …]]` (diff lines 15-22, all
  `-` removals) is GONE. New content mirrors SymPy: `normP` normalization
  (`expectZero`, :67-68); `covDouble` symmetrized identity
  (`expectZero["symmetrized covariance identity", covDouble - cov]`, :72-83);
  `fPrime` closed form (`expectZero["f'(z) closed form", …]`, :89-90);
  `Reduce[Sin[πz/(2lM)] > 0 && 0 < z < lM]` sign decision (:91-95); consistency
  `expectZero["dg/dpiM = -(1/lM) Cov consistency", dgPi + cov/lM]` (:98-99);
  uniqueness bracket (:102-108).

The recorded F1 deviation (removal of the optional SymPy `reduce_inequalities`
probe) is exactly the one the directive's implementation note explicitly
permitted as non-load-bearing.

**Assessment:** Correct and addresses the finding. The certificate is a GENUINE
global proof, not a sample:
- Normalization + symmetrized identity (`Cov = ½∬ symmetrizer · p·p`) +
  symmetrizer `≤ 0` pointwise (strict off-diagonal, because `f = cos(πz/2L)` is
  strictly decreasing on `(0,L)` — certified from `f'(z)`'s closed form and the
  `sin` argument lying in `(0,π/2)`) + `p > 0` ⇒ `Cov < 0` strictly ⇒
  `dg/dΠ = −Cov/L > 0` for ALL `Π > 0`. No interval is sampled.
- Uniqueness: the verified monotone range `g(0+) = 2/π → g(∞) = 1` (limits also
  checked: out lines 3-4 SymPy / 11-16 Mathematica) plus the bracket
  `2/π < g_- < 1` give at-most-one-root and exactly-one-root respectively, hence
  a unique `Π_*`.

This is a sound GLOBAL proof and strictly stronger than the rejected sweep,
matching the Claude+Codex consult Q1 (`_consult_batch4.md:9-11`).

**Non-tautology:** the symmetrized double integral is built INDEPENDENTLY from
`f`, `z`, `sigma` and asserted equal to `Cov = E[fz] − g·E[z]`; it is NOT an
X−X self-comparison. A wrong `gPi`, `sigma`, or `f` would break
`Cov_double − Cov == 0`. The transcripts print `Cov_Pi(f,z) (symmetrized)` as a
substantive, distinct closed-form expression (SymPy out line 5; Mathematica out
line 25) that nonetheless simplifies to match `Cov`. The load-bearing sign
argument (the `f'(z)` closed form + bounded `sin` range) remains in BOTH engines
after the `reduce_inequalities` deletion, so that deletion removed no
load-bearing content.

### F2 — insufficient_verification (boxed closed-form g_Π equality — RECONCILE only)

**Classification:** resolved

**What changed:** Nothing — preserved as required. SymPy
`…sympy_audit.py:16-18` still defines
`gPi_boxed = 2Π(2Π e^Π+π)/((4Π²+π²)(e^Π−1))` and guards
`gPi - gPi_boxed == 0`; Mathematica `…mathematica_audit.wl:43-44` still has
`gPiBoxed` and
`expectZero["g_Pi matches paper boxed closed form", gPi - gPiBoxed]`. The diff
shows no edits to these regions.

**Assessment:** Correct. The check is substantive (a misquoted constant breaks
the simplify-to-zero) and still passes post-F1: SymPy out line 1 prints the
matched closed form; Mathematica out lines 5-6 show
`PASS: g_Pi matches paper boxed closed form`. Both scripts exit 0.

## Exec log assessment

**SymPy:** exit=0. SymPy uses print + AssertionError guards (no PASS vocabulary);
the script ran to completion (printed `Pi_*`, `g(Pi_*)`, `g'(Pi_*)`), so every
guard passed. Notable lines from the committed transcript:
- `Covariance identity residual = 0`
- `Cov_Pi(f,z) (symmetrized) = 2*L*(-4*pi**2*Pi*exp(2*Pi) - … )/((1 - exp(Pi))**2*(16*Pi**4 + …))` (substantive, distinct expression — not X−X)
- `f'(z) = pi*sin(pi*z/(2*L))/(2*L)`
- `Global strict monotonicity certified: dg/dPi = -(1/L) Cov_Pi(f,z) > 0 for Pi>0.`
- `Bracket for unique Pi_*: 2/pi = 0.6366… < g_minus = 0.758035… < 1`
- `Pi_* = 1.50882951349315861144664988336` (unchanged from prior pass)
- NO `dg/dPi at Pi=…` sweep lines anywhere.

**Mathematica:** exit=0. **11 PASS lines**, ending
`Stage 130 Mathematica audit passed.`:
1. `PASS: g_Pi matches paper boxed closed form` (F2)
2. `PASS: covariance identity`
3. `PASS: uniform-source limit`
4. `PASS: point-source limit`
5. `PASS: sigma_piM normalized on [0,lM]`
6. `PASS: symmetrized covariance identity`
7. `PASS: f'(z) closed form`
8. `PASS: f strictly decreasing on (0,lM) -> symmetrizer <= 0`
9. `PASS: dg/dpiM = -(1/lM) Cov consistency`
10. `PASS: g_minus strictly inside (2/Pi, 1): Pi_* unique`
11. `PASS: Family-1 compensation point`

The `Reduce[Sin[Pi z/(2 lM)] > 0 …]` decided to a non-`False` region
(`lM > 0 && 0 < z < lM`, out line 23), so the strict-decrease certificate is
genuinely decided, not assumed. `Pi_* = 1.5088295134931555…` matches the SymPy
value and the prior pass. NO `dg/dpiM > 0 at piM=…` sweep lines anywhere.

**Output freshness:** Both committed `.txt` transcripts reflect the post-fix
scripts: they print the new certificate lines (symmetrized covariance, f'(z)
closed form, consistency, bracket) and contain NO finite-sample monotonicity
sweep lines. The orchestrator independently re-ran both engines to produce these
authoritative outputs; content matches the post-fix source.

## Material-change assessment

`material_change`: **false**. The edit strengthened the verification surface
only (finite sweep → global symbolic certificate). No derived quantity changed:
the `g_Π` boxed form, `Π_* = 1.50882951349…`, `x_* = 0.66276540…`, and
`g(Π_*) = 0.758035…` are identical to the prior pass. No downstream unit's
inputs are altered, so no narrow re-audit is warranted on numeric grounds.

## Side observations (non-blocking)

- The SymPy comment block at :71-72 calls the auxiliary `gprime = diff(-f, z)`
  helper "g(z) := -f(z)", reusing the symbol `g` already used for the bias map
  in prose. Purely a comment-naming overlap; the code variable is `gprime` and
  there is no actual collision. Non-blocking.
- The two engines compute `Cov_double`/`covDouble` in different but equivalent
  algebraic groupings (SymPy factors `(1 - exp(Pi))**2`; Mathematica
  `(-1 + E^piM)^2`); each independently passes its `covDouble - cov == 0`
  identity, which is the relevant check. Non-blocking.

## Verdict justification

Both findings are `resolved`. The headline F1 requirement — eliminate the finite
monotonicity sweep — is met in BOTH engines (confirmed via the diff's `-` removal
hunks and the absence of any sweep line in both committed transcripts). The
replacement is a sound GLOBAL proof: a normalized density, an
independently-constructed symmetrized covariance identity (non-tautological — a
wrong gPi/σ/f breaks `Cov_double − Cov == 0`), a decidable bounded-domain sign
certificate for the strictly-decreasing `f`, a consistency tie-back
`dg/dΠ = −Cov/L`, and a uniqueness bracket `2/π < g_- < 1` over the verified
monotone range `2/π → 1`. F2's boxed-form equality is preserved and still
passing. SymPy exits 0 (all guards passed); Mathematica exits 0 with 11 PASS
lines. No regressions in the diff; no derived value changed, so
`material_change: false`.
