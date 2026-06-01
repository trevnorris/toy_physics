---
unit_id: 197
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 197

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine gap)

**Classification:** resolved

**What changed:**
Codex created the new independent-route Mathematica audit
`mathematica/moving_throat_pde_stage197_conditional_packetA_closure_theorem_mathematica_audit.wl`
(no edit to the SymPy reference, as directed). It re-verifies every load-bearing
assertion of the SymPy script across the same seven sections (I isotropic-P0
projection, II residual collapse, III native z^5 extraction + source-map reduction,
IV finish-line equivalences, V deformation-algebra gate, VI linearization, VII
higher-odd irrelevance) and ends `Exit[0]` after all PASS lines.

**Assessment — genuine independent route (directive's central requirement):**
The `.wl` is NOT a transliteration of the `.py`. The decompositions differ at the
load-bearing core:

- **DtN window source (the decisive difference).** The SymPy script *hardcodes* the
  pre-computed l=2 DtN Taylor coefficients as closed forms
  (`L0=-3S+Σ0`, `L2=Sβ²/3+Σ2`, `L4=Sβ⁴/9+Σ4`, `L5=Sβ⁵/9+Σ5`; py lines 81–87) and
  forms the rational response by hand. The `.wl` instead *derives* the window from the
  actual special function:
  `lambdaOut = FunctionExpand[z*D[SphericalHankelH1[2,z],z]/SphericalHankelH1[2,z]]`
  (wl line 105), then `Series` to order 5 (wl line 106) to obtain
  `lambdaWindow5 = -3 + z²/3 + z⁴/9 + (I/9) z⁵` natively, applies the β-stretch
  `z -> betaStretch*z`, and adds the deformation. This is a different route to the same
  kernel via a native Mathematica primitive (the Hankel log-derivative), not a re-typed
  closed form.
- **Even-moment matching.** SymPy *hardcodes* the matched moments
  `Σ2_match`, `Σ4_match` (py lines 90–91). The `.wl` instead `Solve`s for σ2, σ4 from
  the canonical normalized even targets `c2==1/9, c4==4/81` (wl lines 117–124), checks
  the solution is unique, and only THEN confirms the solved values equal the SymPy forms
  (`Sigma_2/Sigma_4 agreement with SymPy result == 0`, output lines 31–34). Derive-then-
  compare, not assume.
- **χ_Q extraction is double-routed.** The `.wl` extracts χ_Q both by `Coefficient[...,z,5]`
  (wl line 138) AND independently by the fifth derivative `D[...,{z,5}]/5!` (wl line 142),
  and cross-checks the two against each other (wl line 149) before comparing to the SymPy
  closed form (wl line 150). Both internal routes and the cross-engine compare pass.
- **Zero-set equivalence uses `Reduce`.** Where SymPy asserts the reparametrization
  identities algebraically, the `.wl` additionally solves the zero set with
  `Reduce[deltaNormGeneric==0, chiFree, Reals]` and `Reduce[nGeneric==1, chiFree, Reals]`
  (wl lines 177–184), then `Equivalent[zeroSet, chiFree==1]` (wl lines 199–206). The
  equivalence direction is *computed*, not asserted.

**Substantive / non-vacuous (directive + central-question b):** The χ_Q=z⁵ derivation
is genuine. I independently reproduced the `.wl` Section III in a Mathematica kernel:
`lambdaWindow5 = -3 + z²/3 + z⁴/9 + (I/9)z⁵`, the even-match `Solve` returns the unique
`Σ2 = (-3(β²-1)S - Σ0)/9`, `Σ4 = (-3(β⁴-1)S - Σ0)/27` (= SymPy), and the z⁵ extraction
yields `χ_Q = 3(Sβ⁵+9Σ5)/(3S-Σ0)` with `χ - χ_SymPy = 0` — matching the committed
output verbatim. (My first hand-reconstruction used the wrong Hankel sign convention
`-i/9`; Mathematica's `SphericalHankelH1[2,z]` gives `+i/9`, and the script's output is
the correct one — confirming the output is genuine, not stale or hand-edited.) The
off-canonical guard is real: `Delta_norm at chi_Q=6/5 = -9 bigG·soundSpeed⁵/(5 c⁵ a⁵)`,
`expectTrue[... != 0]` PASS — closure genuinely FAILS off the canonical point, so the
`⟸` direction is secured and the equivalence is not vacuous. I confirmed the `Reduce`
equivalence is non-tautological: `Reduce[Delta_norm==0,chiFree]` independently returns
`chiFree==1`, and `Equivalent[chiFree==1, chiFree==2]` does NOT collapse to True, so the
`expectTrue` would fail against a wrong target. No X−X or hardcoded-vs-itself in the
load-bearing checks. (The Section II "first seven slots" and the `Last - deltaNormSlot`
checks are definitional restatements of the carried front end, mirroring the SymPy
decorative checks A3/A5; harmless and non-defect, same as the auditor noted.)

## Exec log assessment

**SymPy:** exit=0 (unchanged reference engine). Notable lines (committed `.txt`):
`chi_Q extractor - deformation algebra formula = 0`;
`Delta_norm at chi_Q = 6/5 = -9*G*c_s**5/(5*a**5*c**5)`;
`(3S - Sigma0)(chi_Q - 1) - closure numerator = 0`; `d chi_Q/dL7 = 0`.

**Mathematica:** exit=0. Notable lines (committed `.txt`):
`chi_Q from Series/Coefficient z^5 extraction = (3*(betaStretch^5*scaleS + 9*sigma5))/(3*scaleS - sigma0)`;
`Series coefficient agrees with fifth-derivative extraction` PASS;
`chi_Q extractor - SymPy deformation algebra formula` PASS;
`Reduce[Delta_norm == 0, chi_Q] = chiFree == 1`; `Delta_norm == 0 iff chi_Q == 1 = True` PASS;
`Delta_norm at chi_Q = 6/5 = (-9*bigG*soundSpeed^5)/(5*lightSpeed^5*radius^5)` nonzero PASS.
All sections PASS; script reaches `Exit[0]`.

**Output freshness:** confirmed. `.wl` mtime 2026-06-01 12:02:15; its output `.txt`
mtime 12:03:34 (newer). The SymPy output `.txt` was also refreshed 12:03:34 (re-run,
unchanged). Outputs are post-fix and consistent with the scripts.

## Material-change assessment

`material_change`: false. No SymPy/derived result changed; this only ADDS a second
engine that reproduces the existing SymPy conclusions exactly. No downstream unit's
derived value moves. (Cross-engine agreement is complete: χ_Q closed form, matched
Σ2/Σ4, Δ_norm, the linearization coefficients 5, 1/(3S), 9/S, and the L7-independence
all agree.)

## Side observations (non-blocking)

- The `.wl` adds a value-add beyond the SymPy script: an independent
  `SphericalHankelH1`-derived DtN window and a `Reduce`-based zero-set proof, both of
  which strengthen the dual-engine story rather than merely mirroring it.
- The `P0_target from Gamma_5 normalization` check (wl lines 99–103) re-derives
  `P0_target` from `Gamma_5 = 2G/(5c⁵)` and `27 c_s⁵ Gamma_5/a⁵`, agreeing with the
  literal `54 G c_s⁵/(5 a⁵ c⁵)` — a small extra independent constant cross-check not in
  the SymPy script. Harmless and correct.

## Verdict justification

`verified`. The single finding (missing second engine) is resolved by a genuinely
independent Mathematica route, not a transliteration: it derives the DtN window from
`SphericalHankelH1[2,z]` natively, solves the even moments rather than hardcoding them,
extracts χ_Q from the z⁵ coefficient by two independent internal methods, and proves the
finish-line equivalence with `Reduce`. I reproduced the load-bearing Section III
extraction in a Mathematica kernel and it matches the committed output exactly (χ_Q
closed form, matched Σ2/Σ4, χ−χ_SymPy=0), confirming the output is genuine. The
off-canonical guard (χ_Q=6/5) genuinely shows closure fails, and the `Reduce`/`Equivalent`
zero-set proof is non-vacuous. Both engines exit 0, outputs are fresh, no regressions.
