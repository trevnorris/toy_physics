# Decision review — the blind Wolfram N6 build directive

## Artifact
`research/pde_ledger_v3/directives/S11c_c2_N6_wl_build_directive.md` (orchestrator-written; a physics-bearing build
directive for a NEW blind Mathematica engine). Review it **until clear**: report every issue that would change what the
WL engine computes or what may be claimed from it. This is a DECISION review of a directive — the `.wl` engine does not
exist yet, so ⛔ do not attempt a fictional-script ablation; executable script-control tests are deferred to the build.

## What this directive is for
The per-engine (SymPy) c2 N6 (representation-invariance) question is RESOLVED: N6 passes as **operator covariance**
(Reading B), established by (1) the **carrier reconcile** `C_E − C_M` (geometric carrier is representation-independent)
and (2) the **source-naturality** residual `R_cov = ms − ms_pred` (the commuting square). This directive commissions the
**blind Wolfram sibling** that must INDEPENDENTLY reproduce both checks, importing nothing. ⛔ **You are reviewing the
DIRECTIVE, not re-adjudicating Reading B** — the covariance verdict is settled and is NOT to be re-litigated; check
whether the directive faithfully and completely commissions the blind WL reproduction.

## Source-of-truth to read FIRST (form your own view before judging the directive), and quote both sides
- `research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` §5c (the N6 object, the fixed-`(α,ρ)` rule, the
  category-error prohibitions), §§1–2, §6, §7 (incl. the "blind Wolfram engine re-derives … importing nothing" clause).
- `research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_directive.md` (cleared) — the carrier/source decomposition,
  the exact three-way affine split, the carrier-bridge direct-coefficient test, `SPLIT_CHECK`.
- `research/pde_ledger_v3/directives/S11c_c2_N6_covariance_directive.md` (cleared) — `R_cov`, the Φ-prolongation crux,
  non-circularity, the two knives, the `R_COV_INCREMENT = closed_response` pin.
- `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py` and `.../S11c_c2_N6_covariance_sympy.py` — the CLEARED
  SymPy instruments (the source-of-truth for exactly what the two checks compute; the WL engine must reproduce the same
  OBJECTS, blind — ⛔ never the same representation).
- `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md` (cleared) — the native material carrier `S_{P,M}`.
- `research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md` — the verdict FORM + the three carry-forward caveats.
- The blind-WL precedent + shared-symbol conventions: `research/pde_ledger_v3/directives/S11c_c1_wl_build_directive.md`.

## Required method (DOCUMENT review — derive your own view from the sources, then judge the directive)
For each finding, quote the directive AND the source it contradicts or omits, with `file:line`. A physics claim in your
report is worth nothing without a `file:line` in a source — ⛔ a prose re-derivation is discarded. Report a finding only
if it catches a way the WL engine, built to this directive, could **compute the wrong thing** or let a **wrong claim** be
made. Check specifically:

1. **Faithful engine-neutral translation.** Do the directive's engine-neutral object statements (§"The two checks,
   stated engine-neutrally") match what the cleared SymPy instruments actually compute? In particular: (a) the carrier
   `C` = the pressure-slot coefficients `C_{r,p}` with the P-projection; `C_E` from the Eulerian slab rows, `C_M` from
   the material face fold; (b) the source amplitude `b_{r,s}` and the `es`/`ms` binding (Eulerian μ_E+V_E vs material
   μ_M+V_M); (c) the increment `I(C,S) = −C·p + B(C,S)` affine, and the exact three-way split `R_N6 = I(ΔC,ms) +
   B(C_M,ΔS) + B(ΔC,ΔS)` with CARRIER using `ms` (not `es`) and SOURCE/CROSS closed-response-only (sig 6/9/12); (d) `Φ`
   prolonged to the jet order PRESENT in `μ_E` (rank-2), NOT a fields+first-jets dict; (e) `R_cov = source(μ_M,V_M) −
   source(μ_E∘Φ,V_E)` and its NON-circularity (prediction from imported/re-derived μ_E + supplied Φ, ⛔ never the
   material pullback); (f) `R_COV_INCREMENT = closed_response(C_M,R_cov)` sig 6/9/12, ⛔ not `build_increment`.
2. **Anything MISSING for the two checks.** Does the directive commission everything the carrier bridge + `R_cov` need
   (both face routes of the S11c-a substrate; μ_E and μ_M constructed, not opaque; the c1 closed response opaque-to-
   coordinate-map; the four pressure slots; V_E and V_M both computed)? Is the "N6-scoped, not the full c2 operator"
   scope actually SUFFICIENT — can the carrier bridge and `R_cov` be built without the full §3a close / §3b re-extract,
   or does the split need machinery the directive omits?
3. **Blindness by ABSENCE, not denylist.** Is blindness enforced by the engine importing nothing (re-derive from
   specs), ⛔ never by a do-not-read sentence? Does anything push the WL engine to reproduce the SymPy
   representation/tag values/shape (designed-to-agree, rule 16)? Is any SymPy-API name left as an instruction a WL
   builder cannot act on?
4. **Rule-17 freezes.** Does the directive anywhere freeze a varying field/coefficient (background `W_bg`, density,
   jets, the advection `a_ρ`, the thickness shift `h_α`) where it must stay live? (`a_ρ=u·∇ρ₄/ρ₄`, `h_α=u·∇W_bg/W_bg` at
   LAB_HELD else 0; `RHO4_CONSTANT ⇒ g_i=0 ⇒ a_ρ=0` is a computed absence, ⛔ not a freeze.)
5. **Controls able-to-fail and FORM.** Are the four controls (carrier covector/normal FORM knife; Φ-coefficient knife
   2·a_ρ; θ-independent junk knife κ_j·J_μ·e_W; source `a_ρ→0` bite) genuinely one-sided, able-to-fail, and FORM (not a
   coefficient rescale, not an A−A)? Would each move the object it targets and leave the others fixed?
6. **No leakage / no pre-judgement.** Any expected value, residual-zero exit/assert, or VERDICT? Is the disposition of
   `R_N6`/`C_E−C_M`/`R_cov` left unjudged (raw nonzero ≠ disagreement)? Is the "⛔ this engine must not subtract its own
   discrepancy and call it covariance" instruction present and correct?
7. **The three caveats.** Are the three carry-forward premise caveats (Φ physical-correctness; V transform; extracted-
   block leakage) carried as NOT closed by this engine, routed downstream — ⛔ not silently pre-cleared?

## What you are handed
The directive, the source-of-truth files above, and read access to the repo. You are NOT handed an expected answer —
there is none; the directive withholds every computed value. ⛔ Do not propose making the WL engine agree with SymPy;
propose only fixes that make the blind reproduction correct and complete.

## Output
A list of findings, each: `file:line` in the directive, the contradicting/omitted source `file:line`, why it changes
what is computed or may be claimed, and the minimal fix. If a section is sound, say so briefly. End with an overall
disposition: CLEAR-TO-BUILD, or FOLD-REQUIRED with the blocking items.
