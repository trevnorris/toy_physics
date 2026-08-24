# Measurements twin — S11c_a_wl_engine_eval_fix_directive.md

## The measured errors (the claim carries the command)
CMD: wolframscript -file .../S11c_a_interface_geometry_mathematica_audit.wl > run.out 2> run.err   (bg b69732hsc, exit 0)
  run.out: 13705734 bytes, 40 distinct WL_S11CA_ tags (complete parallel set).
CMD: grep -av 'WolframScript.conf' ~/.s11_build/S11c_a_wl_engine_run.err
  ⇒ Part::pkspec1 (expression index cannot be used as a part specification) — appears FIRST;
    AssociateTo::argrx (called with 3 arguments; 2 expected). Both suppressed after a few by General::stop.
CMD: for pat in AssociateTo 'Part[' '$Failed' Missing pkspec; do grep -acE "$pat" run.out; done ⇒ all 0
  ⇒ no error residue in payloads, BUT AssociateTo::argrx returns unevaluated ⇒ DROPS the entry silently.

## Bug 1 root cause (orchestrator-located)
CMD: sed -n '1256,1264p; 1386,1394p' .../S11c_a_interface_geometry_mathematica_audit.wl
  ⇒ AssociateTo[formBaseOperand, virtualKey -> objectRecord[...], evolutionKey -> objectRecord[...]]  = 3 args (2 rules).
    Same at uniformS11ca/uniformS11b (1386-1394). FIX: AssociateTo[assoc, {r1, r2}] (list) or 2 calls.
  Single-rule AssociateTo (rep-invariance 1006-1055, independence 1080-1127) are FINE (2 args) ⇒ §5a controls INTACT.
  objectRecord[expr_, dim_, epsOrder_:1] (line 79) — 2/3-arg optional, not the culprit.
  ⇒ affected: CONTROL_FORM (§5b) + UNIFORM_LIMIT (§5c) missing VIRTUAL_WORK + EVOLUTION keyed entries.

## Bug 2 (Part::pkspec1) — NOT localized by orchestrator; Codex to trace (appears first ⇒ possibly a primary task).
