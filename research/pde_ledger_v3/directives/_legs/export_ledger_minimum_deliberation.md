# Design deliberation — what is the MINIMUM that must be EXPORTED (plain-git LEDGER) vs EMIT-ONLY (.out)?

You are one of two independent engines deliberating a design question. Reason from first principles and from
the actual files; do **not** try to agree with the other engine — an independent view (and any disagreement) is
the point. Cite exact rows / spec lines / rule text for every claim; ⛔ a prose "should be fine" is worthless.

## The question
The export chain has **two channels**:
- **EXPORT** — the plain-git `*_exports.py` **LEDGER**, imported by the next chain step (`N1`), which binds
  objects from it **by key** to build its own objects and never re-derives an inherited object. It must stay
  under GitHub's 100 MB plain-git file cap (an `*_exports.py` is ⛔ never annexed).
- **EMIT** — the stdout tag stream captured to the annexed `out/*.out` (git-annex + GIN, **no size cap**), read
  by the cross-engine comparator and the review legs. It carries **everything** the engine prints.

**What is the MINIMUM the LEDGER must contain for the chain to function** (i.e. for every downstream step to be
able to build), and what can therefore live **only** in the `.out`? We are **clearly exporting too much** — I
want a principled answer, not a case-by-case trim.

## Measured problem (the current inherited `S11c_b_exports.py` LEDGER — examine it)
- Total **58.8 MB**. **29% (17 MB) is `*_term_origins` rows** — provenance-trace DIAGNOSTICS.
- `slab_operator` is **18 MB** across **8 cases** (~2.25 MB/case, ~75K SymPy nodes/case). It is fully
  **expanded**: **271 distinct symbols but 365,820 symbol occurrences** (~1,349× — a flat sum of thousands of
  monomials). `coupling_kernel` 7.6 MB. `mu_theta_operator` similarly large.
- Row schema (each row): `display` (rendered str) + `value` (`_restore("<srepr>")`) + `class` + `step` +
  `f9_operands` + `corroborated_steps` + optional `dimension_key`. So each expression is stored **twice**
  (display + srepr), and srepr is ~4–5× the rendered form.

## What to read before answering
- The chain rules: `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` — **F1–F10**
  (esp. the new **F10**: "LEDGER carries the model-level register; step-level diagnostics are emit-only"), and
  `research/pde_ledger_v3/S9_REWRITE_PLAN.md` **D1** ("export every MAIN object; no per-object judgement; errs
  toward exporting more") and **D3** (the round-trip reviver).
- What a downstream step actually binds: `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` §0/§1/§3b
  (c1 imports S11c-b's `slab_operator`, `coupling_kernel`, `mu_theta_operator`, the T-substrate, and the
  constants — and NAMES its own consume-set for c2), and `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`.
- The LEDGER itself: `research/pde_ledger_v3/scripts/S11c_b_exports.py` (grep row keys, `class`, `step`; look at
  a `*_term_origins` row, `slab_operator`, a knob row like `c_s0`, `Lambda_A_0`).
- The comparator (what reads which channel): `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py`
  (line ~1206 `exact_rational_residual`; the header "S11c_b_exports.py is not used as a tag stream").

## Deliver (be concrete and cite)
1. **The PRINCIPLE.** State the minimal rule for what the LEDGER must carry for a downstream step to build.
   Frame it in terms of what a later step BINDS (composes / substitutes / evaluates), not by object category
   name. What is the smallest closed set?
2. **CATEGORIZE the current LEDGER content against that principle** — which of these must be exported vs are
   emit-only, and WHY (cite what binds them, or that nothing does): (a) knob/constant rows (`c_s0`, `rho_m`,
   `Lambda_*`, `tau_*`, `W_0`, …); (b) varying-profile rows (`W_bg`, `w1_profile`, `sigma_W`, …); (c) the
   model-definition operators/kernels/responses a later step composes (`slab_operator`, `coupling_kernel`,
   `mu_theta_operator`, the T-substrate, the DtN/face-response c1 will export); (d) `*_term_origins`; (e)
   dissipation/energy/loci/noninvertibility/Hermitian-reactive/regime-parity structural views; (f) §5/§6 control
   operands; (g) intermediate results named but not bound downstream. For each: EXPORT or EMIT-ONLY + the
   binding (or its absence) that decides it.
3. **The MINIMUM set**, and the delta vs what is currently exported (what should be CUT). Quantify against the
   58.8 MB where you can (e.g. the 29% term_origins).
4. **Under-export risk + how it's bounded.** The `.out` holds everything (recoverable); each step's spec
   declares its consume-set; worst case a step re-derives. Is `D1`'s "errs toward more" default still needed, or
   does F10 + the `.out` safety-net replace it? Where is the real danger, and what guard removes it?
5. **Representation vs set (keep these separate).** Even a minimal SET still includes big operators a later step
   binds (e.g. `slab_operator`, which c2 folds into). Note — but do NOT solve here — whether those bound objects
   must be stored in the **expanded** form, given the comparator canonicalizes the difference itself
   (`sp.cancel(sp.together(sp.expand(A−B)))`). Flag it as a separate lever; the question here is the SET.
6. **Does your answer revise F10 or D1 as written?** If so, how.

## Output
A principled minimal-export rule, the categorization, the cut-list with rough MB, the under-export risk
analysis, and an explicit note of any disagreement you'd expect with the current F10/D1. End with the single
sentence that would go in the export-chain decision list as the rule.
