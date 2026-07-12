# Build + run the native-`Pᵃ` constraint-class gate (dual-engine), iterate to green

**Implement and RUN the gate specified in `software/em_charge_attribute/directive_native_p_constraint_gate.md` (v5, `DIRECTIVE_CLEAN`).** You are the builder: design the scripts, run them, iterate until they run clean and produce the report. Reasoning effort: **xhigh**. Follow the directive exactly; **choose no physics that affects the verdict** — freeze + DOCUMENT the flagged build-implementation choices (explicit operator basis + completeness, spatial topology, boundary/test-function class, puncture treatment, the control gauge-fixing pairs) in the report.

## What to build (dual engine — REQUIRED)
- **SymPy** (`software/em_charge_attribute/native_p_gate_sympy.py`): the full functional Dirac–Bergmann analysis for **THEORY-A and THEORY-C** — Legendre map, primary/secondary constraint chain, class of each constraint via weak Poisson brackets, the search for an ADDITIONAL first-class primary→Gauss generator `G` as a function of the frozen couplings `{g_a}`, the source-free-then-one-linear-`j`-coupling test, the reduced-physical-mode shear-duplicate compare, and generic/tuned/symmetry-protected classification. Plus the SIX controls with their runtime able-to-fail guards (computed, not source-grep).
- **Mathematica** (`software/em_charge_attribute/native_p_gate_dual.wl`): INDEPENDENTLY re-derive the constraint chain, classes, and the `G`-existence verdict (per theory) — a genuine second engine, not a port of the Python.
- **Comparator + runner** (`native_p_gate_compare.py`, `run_native_p_gate.sh`): assert `ENGINE_AGREE` per theory on {constraint classes, `G` existence, verdict}; runner self-resolves its dir, each step `timeout 600`.

## Acceptance (iterate until ALL hold)
- Every script runs clean (exit 0); no step > 10 min (a timeout = reformulate the math, never raise the cap).
- **All SIX controls pass their able-to-fail assertions** (Maxwell → first-class Gauss detected; gauged-hard-unit → MIXED detected; bare O(4) sigma → second-class radial, no extra Gauss; non-conserved current → inconsistent preservation; gauge-fixed Maxwell → second-class/no-local-gauge; global-only → no local Gauss). **Per-tooth ablation:** confirm each assert fires at its own point.
- `ENGINE_AGREE` per theory.
- The report `software/em_charge_attribute/reports/native_p_constraint_gate.md` is produced with, per theory (A, C): the frozen setup + exhibited basis + documented implementation choices; the full Dirac chain (constraints + classes + matrix); the `G`-search result; the six controls; dual-engine logs; and the one-line **VERDICT** from the directive's decision table.

## Notes
- Mathematica: this session has `--sandbox danger-full-access` so `wolframscript`/`math -script` can run; do not exceed ONE concurrent Mathematica invocation.
- Do NOT weaken the directive to make a script pass. If the math forces a genuine `GATE_ILL_POSED` or a control cannot be made to pass honestly, REPORT that — a computed no-go / ill-posed is a first-class result. Do not fabricate `ENGINE_AGREE`.
