---
unit_id: 246
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 246

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage246_..._mathematica_audit.wl` (not in the
captured diff — it is a new add — but present on disk, 180 lines). Covers M1–M9 with
`expectZero`/`expectTrue`/`expectApprox`, exits 0, all 16 PASS lines present in the log.

**Assessment:**
Genuinely independent route, not a port of the `.py`:
- M1 (`wl:56`): `Integrate[source,{x,0,1}]` directly.
- M2/M3 (`wl:57-67`): native `Integrate` of `source·Cos[πx/2]` and `source·K_q` against the
  closed forms — no hand product-to-sum (the `.py`'s choreography is not reused).
- M4 (`wl:90-97`): the σ_min claim is DERIVED via `MinValue[{source/.rules, 0<=x<=1}, x]`
  and `expectApprox`'d against the paper's piecewise targets 3/10, −5/4, 7/16 at tol 10^-20.
  This is exactly the check the SymPy script omits — the minimum is computed from the source,
  NOT re-encoded from the piecewise. This independently closes the F2/F3 tautologies from the
  second engine.
- M5 (`wl:101-115`): `Det[Msrc]` and `Inverse[Msrc].{gt,St}` — native primitives.
- M6 (`wl:119-124`): quarter-ratio identity, both branches.
- M7 (`wl:126-127`): compensation line via `Solve`.
- M8 (`wl:129-156`): boundary candidate `1+b[r]−Abs[a[r]]` under signed assumptions, threshold
  via `Reduce` (the directive explicitly permits the boundary-candidate form here).
- M9 (`wl:160-176`): numeric readback pinning g, S, R, σ_min at tol 5e-9, plus g(0)>1, rThr>reval.
No false positivity assumptions added for the symbolic integrals (a,b,gt,St declared Reals only).

### F2 — tautological_check (quadratic rewrite + σ_min self-test)

**Classification:** resolved

**What changed:**
- `sympy:56`: new `assert simplify(sigma - (1 - b + a·cos(πx) + 2b·cos(πx)²)) == 0`.
- `sympy:88-95`: new independent true-minimum check inside the test loop.

**Assessment:**
- The y=cos(πx) substitution is now genuinely tested: line 56 connects the trig σ(x) to the
  quadratic via `cos(2πx)=2cos²(πx)−1`. Non-tautological (the old `sigma_y` self-equality at
  line 54 is retained only as documentation; the substantive assert is line 56).
- σ_min is now checked against an INDEPENDENT minimum: `sigma_min_true = Min(*cand)` where
  `cand` are σ_y at y=±1 plus the vertex only when 2b>0 and y*∈[−1,1] — computed directly from
  the quadratic, NOT from the piecewise branch logic. Verified by hand the three points: (1/2,−1/5)→
  min(1.3,0.3)=3/10; (5/2,1/4)→vertex y*=−5/2 out of range→min(3.75,−0.75)=−5/4; (1/2,1/2)→
  vertex y*=−1/4 in range→7/16. A wrong vertex coefficient (e.g. a²/(7b)) would now fail. The old
  `sigma_min_expected` self-test at line 96 is retained but harmless.

### F3 — tautological_check (transported σ_min branch selection)

**Classification:** resolved

**What changed:**
`sympy:212-220`: introduces signed orientation symbols `a0p` (positive), `b0n` (negative),
substitutes them into `sigma_min_piece` from Section 2, `piecewise_fold`s, and asserts the
result equals `1 − (a0p − b0n)·s_r`.

**Assessment:**
The transported minimum is now tied to the Section-2 piecewise (`sigma_min_piece`), not a
self-rearrangement. The signed assumptions force `b_r_or<0` → the `Le(b,0)` branch is selected,
and `Abs(a_r_or)=a_r_or` (a0p>0, s_r>0), so the fold reduces to the boundary form. A wrong
branch boundary `4b` or a sign error in the boundary piece would surface. The original line-209
tautological assert is retained as documentation per directive but the new assert is the
substantive one.

### F4 — insufficient_verification (S readback)

**Classification:** resolved

**What changed:**
`sympy:257`: `assert abs(float(S_eval) - 0.67584771) < 5e-9`.

**Assessment:**
The 0.67584771 Session-I readback is now ASSERTED (was print-only). Matches the g/R/σ_min
pattern. Log shows S_eval = 0.6758477114656324, within tolerance.

## Exec log assessment

**SymPy:** exit=0. All symbolic checks passed; printed values unchanged from pre-fix
(e.g. `S[sigma](r_eval) = 0.6758477114656324`, `sigma_min(r_eval) = -0.08979545...`). The four
new asserts raise nothing.

**Mathematica:** exit=0. 16 PASS lines, all M1–M9 covered:
- `PASS: M4 minimum at a=1/2,b=-1/5/...,b=1/4/...,b=1/2` — minimum derived via MinValue, delta=0.
- `PASS: M9 S(r_eval)` delta 1.47e-9; `PASS: M9 g(r_eval)` delta 4.08e-9; all within 5e-9.

**Output freshness:** confirmed. SymPy output mtime 1780502817 > script mtime 1780502394;
Mathematica output mtime 1780502817 > script mtime 1780502586. Both `.txt` outputs regenerated post-fix.

## Material-change assessment

`material_change`: false. All four edits are verification-strength additions (new asserts and a
new independent engine). No derived constant, formula, or printed numeric value changed — the
diff inserts only assertions and a comment block; every Section 1–9 printed value is identical to
the pre-fix output. Downstream Stage 247 consumers of σ_min / quarter-ratio see no value change.

## Side observations (non-blocking)

- F2 line 54 and F3 line 209 (the original tautological asserts) remain in place. Per the directive
  they were explicitly kept as documentation alongside the substantive new asserts; this is by
  design, not a defect.
- M8's Mathematica check uses the boundary candidate `1+b[r]−Abs[a[r]]` rather than `Minimize`;
  the directive explicitly permitted the boundary-candidate form for M8, and the threshold is still
  derived via `Reduce`. Acceptable.

## Verdict justification

All four findings are genuinely resolved. F1 adds a real second-engine route whose M4 derives the
minimum via `MinValue` (the precise check the SymPy script omitted) and computes moments via native
`Integrate`, not a transliteration. F2 now tests the y=cos(πx) substitution and checks σ_min against
an independent finite-candidate minimum. F3 ties the transported σ_min to the Section-2 piecewise via
signed-orientation branch selection. F4 asserts the previously print-only S readback. Both engines
exit 0, outputs are fresh, the diff shows only the four targeted insertions with no collateral edits,
and no derived value changed → material_change false.
