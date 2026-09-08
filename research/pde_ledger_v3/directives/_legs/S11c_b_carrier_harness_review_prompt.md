# Independent review — S11c-b carrier ablation HARNESSES (astra-written scripts)

You are one of two independent review legs (fresh Claude agent + Grok) for two **committed ablation harness
scripts** written by a code engine (gpt-6-astra) against a fixed directive. Working dir `/var/projects/toy_physics`.
Repo read access; write/run check scripts under `/tmp` ONLY; ⛔ do NOT modify the working tree (copy to /tmp and
ablate the copy). **Derive/verify independently — prose is discounted unless grounded in a cited line or a check
script whose path + literal stdout you report.**

## Artifacts under review (the harnesses)
- SymPy harness: `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_carrier_ablation_harness_sympy.py`
- WL harness:    `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_carrier_ablation_harness.wl`
  (+ its driver `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_carrier_ablation_harness_wl.py`)
- Their run transcripts (evidence): `.../\_measurements/S11c_b_carrier_ablation_harness_{sympy,wl}.md`

## The spec they must satisfy (read it)
- Directive (the knife-list they implement): `.../directives/S11c_b_carrier_ablation_harness_directive.md`
- Orchestrator adjudication key (the expected cones — for your cross-check only):
  `.../\_measurements/S11c_b_carrier_ablation_harness_adjudication_key.md`
- The engines the harnesses WRAP (must not be re-derived):
  `.../scripts/S11c_b_brane_operator_sympy_audit.py`, `.../mathematica/S11c_b_brane_operator_mathematica_audit.wl`

## What to check — report a finding only if it makes the harness certify the wrong thing or fail
1. **Wrap the live engine, don't re-derive.** SymPy must `import` the engine and call
   `build_operator("MATERIAL_ADVECTED","RHO4_CONSTANT","EULERIAN")` (accessing rows via `named_tuple_row`, since it
   returns a nested `casify` Tuple); WL must load definitions only + call `evaluatedModel[...]["OPERATOR"]` (NOT the
   emit `Do`/`extractCouplingData`). ⛔ A hand-typed carrier / re-implemented `face_factory` / re-derived rows is a
   top defect. Confirm the harness computes on the engine's own emitted rows.
2. **Each knife = the directive's ONE site + ONE FORM.** K_A = delete the **expanded** `Lambda_A_0`-bearing closure
   addends (⛔ not a `Lambda_A_0→0` substitution, ⛔ not the `:408` bind), leaving Λ_V; K_T = remove the whole
   traction channel (all four `face_u`/`face_e` additions / WL `virtualWork=0`); K_W = native-pressure face collapse
   `delta_p_minus→delta_p_plus` (+ d_w) via `substrate_substitutions` / WL `pressureField[-1]:=pressureUpper`.
   Verify each patch is on the object's **construction**, applied by copy-and-patch of the production source, not a
   reimplementation.
3. **Observation = complete carrier.** `∂(row)/∂(native pressure atom)` FIRST then `→0`; every row × every slot
   (order `(δp+, ∂_wδp+, δp−, ∂_wδp−)`); printed for baseline + each knife + diff. ⛔ No selection, ⛔ no cone
   labels, ⛔ no expected-outcome annotation baked into the harness.
4. **PRINT-not-PASS.** ⛔ No `PASS`/`FAIL`/verdict/"bites" payload; ⛔ no assertion that a diff is zero/nonzero.
5. **⛔⛔ ABLATE THE HARNESS ITSELF (the core check).** Do NOT trust the harness's self-report. Copy the harness to
   /tmp and mutate it: (a) turn a knife's FORM into a **no-op** (patch nothing) — the printed diff for that knife
   MUST go identically zero (else the harness fabricates a bite); (b) turn a knife's FORM into a **coefficient
   rescale** — confirm the harness does not misreport that as a structural bite; (c) confirm the harness's own three
   self-tests (extractor-order `→0`-before-`∂`; dead-path at a pressure-free site `:2970`/`:1347`; ×2 live rescale)
   are genuine — re-run them and report the literal diffs. Report the ablation script paths + literal stdout.
6. **Drift guard.** The baseline carrier matches a fresh canonical run (imported `build_operator` /
   definitions-loaded `evaluatedModel`) vs an unablated temp copy — coherent and actually compared, not asserted.

## ⛔⛔ EXECUTION BUDGET — hard-bound EVERY run (a prior leg WEDGED on an unbounded WL kernel)
⛔ **Wrap EVERY execution in `timeout`** — `timeout 420 python3 …` for SymPy, `timeout 600 wolframscript …` for WL.
A timeout hit is **not** a reason to wait or retry: **report it as a finding ("<engine> <step> exceeded budget") and
MOVE ON.** ⛔ Never raise a timeout, ⛔ never block waiting on a run, ⛔ never run more than one WL kernel at a time.
⭐ **Minimize expensive runs — the canonical S11c-b build is ~4–5 min each:**
- Do the wrap-engine / site-FORM / print-not-PASS / complete-carrier checks by **READING** (items 1–4, 6) — no
  execution needed.
- For the ablate-the-harness check (item 5), you MAY take astra's committed transcript
  (`_measurements/S11c_b_carrier_ablation_harness_{sympy,wl}.md`) as the **baseline** and run only the **mutated**
  copy — so each knife costs ONE build, not two. Do **one** no-op mutation per knife (confirm its diff → 0) and
  **one** rescale (confirm it is not misreported as a structural bite), plus re-run the harness's own three
  self-tests once. That is sufficient; ⛔ do not re-run every variant from scratch.
- If the WL kernel ablation times out, say so as a finding and fall back to the **static** WL read + astra's WL
  transcript — ⛔ do NOT wedge waiting for it.
⛔ Copy to /tmp; ablate the copy; ⛔ never modify the working tree. ⭐ Save every ablation script + its literal stdout
to named absolute paths and report them.

## Output
Per-item verdict with cited lines / check-script stdout, then a final line: **SOUND** (harnesses correctly certify
the carrier and cannot self-report) or **NOT-SOUND** (exact file + line + fix for every blocking issue).
