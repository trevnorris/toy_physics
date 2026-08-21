# Codex consult — adjudicate my S11b WL-engine review disposition (be adversarial)

I am the orchestrator. I ran two independent review legs on a blind Wolfram engine, adjudicated their
findings, and wrote a disposition. Your job is to INDEPENDENTLY check my adjudication and the path I
propose — and to try to REFUTE each finding, not confirm it. A consult that rubber-stamps is worthless;
the SymPy round-2 disposition on the sibling engine was WRONG and a Codex consult corrected it, so assume
this one may be wrong too. Where you agree, say why in your own computed terms; where you disagree, show it.

## Read these (all under /var/projects/toy_physics)
- My disposition: `research/pde_ledger_v3/steps/S11b_wl_engine_review_disposition.md`
- The blind WL engine under review: `research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl`
  (key lines: `equationZeroForm` L367; `ZPERM_SLICE_MAP` L855-867; `solveFace` face law ~L293-300;
  `PRESSURE_WORK_SIGN_CHECK` L1168-1181; `energyEquationRules` L1116-1135)
- The sole physics authority (spec): `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`
  (§4 face closure + affinity ~L208-254; B2c acceptance obligation ~L676-704)
- The sibling SymPy engine's map computation: `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py:1329-1334`
- The Grok leg raw output: `/home/trevnorris/.s11_build/s11b_wl_grok_review.txt`
- My independent derivation: `research/pde_ledger_v3/steps/_measurements/s11b_zperm_slicemap_sign.py`

## Questions — DEMAND COMPUTATION of yourself where a sign or value is at stake
1. **F-WL-1 (the sign).** Before reading my verdict, derive it yourself: from spec §4
   (`J = Λ_A·𝒜 + Λ_V·V`, `𝒜 = μ_s − δp/ρ_m`), on the `μ_s=0` slice, what coefficient does the raw-pressure
   closure `J = Λ_p·δp + Λ_V·V` force for `Λ_p⁰` in terms of `Λ_A⁰`? Write and RUN a script; paste stdout.
   Then read `ZPERM_SLICE_MAP` (L855-867) with `equationZeroForm` (L367) and decide: is the engine's emitted
   `+Λ_A/ρ_m` a genuine bug, or is it defensible under some reading of the spec's wording ("raw-pressure-
   driven WITH COEFFICIENT Λ_p⁰")? Is there any convention under which `+` is correct? Is my "residual-
   coefficient extraction" diagnosis right? ⛔ Do not take my word — show your algebra and the code.
2. **F-WL-2 / F-WL-3 (dead checks).** Is `PRESSURE_WORK_SIGN_CHECK` truly tautological (residual ≡ 0 by
   construction)? Are the energy/two-port/causality/kernel-orientation/grazing checks actually decoupled
   from the assembled system, or did the legs (and I) misread which objects they consume? Point to the lines.
   Could any of these be made genuine, or must they be re-emitted honestly ("no independent second route")?
3. **X-1 (basis count 11 vs 10).** Is `(∇·u)² ≡ tr((∇u)²)` modulo a total divergence (so the two counts
   are a mod-divergence convention difference), or is one engine simply wrong? Which count SHOULD the ledger
   adopt, and is that a comparator finding or a step-record decision?
4. **The path.** I propose: build+freeze the T7 comparator FIRST, run it against the current engines to
   MEASURE F-WL-1 + X-1 + any dimension-tag mismatch, THEN fold all findings into ONE WL repair directive
   (2 legs before the builder), rather than repairing the WL engine now. Is comparator-first right, or is
   there a concrete reason to repair first? Does repairing F-WL-1 before the comparator records it amount to
   "designing away the disagreement," or is it just fixing a bug? Give a recommendation.
5. **Anything I mis-adjudicated or missed** — a finding I dismissed that is real, a "correct" value that is
   actually wrong, or a repair-scope item not on my list.

## Output
Under ~30 lines. For each of 1-5: your verdict (agree / disagree / uncertain), the computation or line
cite that supports it, and — for anything you'd change — the concrete change. Save any script + stdout to
named absolute paths and report them. This is a consult, not a directive; you are checking my reasoning.
