# Independent physics review — S11c-a SymPy engine, traced-bulk shifted trace (§3c)

## Artifact (a SCRIPT)
`research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py` — the SymPy interface-geometry
engine. It emits interface shape-derivative tags `S11CA_*`.

## What to check
Every object that TRACES a bulk field — `RELATIVE_FLUX`, `TRACTION`, `KINEMATIC_BALANCE`,
`EVOLUTION_MASS_BALANCE`, `EVOLUTION_TERM_ORIGINS`, `VIRTUAL_WORK_SHAPE_DERIV`, `CLOSURE_SHAPE_DERIV`,
`FACE_SHIFT`, and the representation / control-form / control-independence / uniform-limit packages built
from them — must implement the §3c trace-linearisation law
`δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s·∂_w f⁰(x,h_s⁰)` for the traced bulk velocity, pressure, current, and
density, consistently with the supplied background state. Decide whether it does.

## What you are handed
- The engine (above).
- Source of truth: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — read §1, §2a, §2d, §3b, §3c.

## Required method — SCRIPT branch (a prose derivation is discarded)
1. Derive it YOURSELF first. Write a sympy script that composes `f = f⁰ + ε·δf` at the moving face
   `h_s = sW_bg/2 + ε·δh`, `W_bg = W₀(1+η·w₁)`, expands to first order in `ε` and first order in `η`, and
   prints the correct linearised trace for (a) a field whose background is zero and (b) a field whose
   background is genuinely nonzero (nonzero normal derivative). Save script + literal stdout to named
   absolute paths and report them.
2. Then check the engine against YOUR derivation, asking of each traced physical object: which line computed
   this? Specifically —
   - Does any traced physical object (velocity/pressure/current/density) contain a background NORMAL-JET
     symbol that is not `∂_w` of a supplied background — i.e. a fabricated free premise?
   - Is the traced perturbation evaluated at the background face `h_s⁰ = sW_bg/2` (so its dependence on the
     background face position is present), or frozen at the flat reference `w = sW₀/2`?
   - Is the conormal derivative's generic probe field (an operator on an ARBITRARY field) left with its
     shifted-evaluation term? (It should be; it is not one of the supplied physical fields.)
   Report any `assert` that precedes the value it guards.

## ⛔ MANDATORY FORM ABLATION (rule 14 — the only thing that has caught the worst defect)
Ablate load-bearing STRUCTURE in a /tmp COPY and report the LITERAL before/after diff. A coefficient rescale
tests arithmetic; only a FORM change tests physics. Targets, at minimum:
- The evaluation face of the traced perturbation (change `h_s⁰ = sW_bg/2` to the flat `sW₀/2`).
- A supplied background value that feeds `∂_w f⁰` (change one physical background from its supplied value /
  `w`-independence to a nonzero `w`-dependent form).
Change the structure, re-emit the affected traced-bulk object, and report the literal diff for each. Judge
from the diff whether the object is load-bearing and whether the physics is right. You need not re-run all 47
tags — import the engine and re-emit a single traced-bulk object per ablation.

## Physics filter
Report a finding only if it catches a way the physics could be wrong, not "wrong on a different input."

## Ablation sandbox
Copy the engine to /tmp and ablate the COPY. ⛔ Never modify the working tree. A run may take a few minutes
and emits large output; background it and watch for progress. Save every ablation script + its literal stdout
to named absolute paths and report them.
