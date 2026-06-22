# STATUS — where the Path-A program is (single front door)

**This file is the canonical "you are here."** It is a thin pointer, not a copy — the detail lives in the linked docs.
Updated at every milestone (same moment `software/stage1_solver/decisions/13` §0 is updated). Last update: **2026-06-22.**

---

## Current state (one paragraph)

The toy model targets the GR-quadrupole verdict `P0·χ_Q·g_mhat²·λγ⁵/g_G = 54/5`. Under the **calibrate-predict** discipline
(every value DERIVED or a declared calibration gap — no silent knobs), the factors now stand: `P0` extracted, **`χ_Q ≈ 0.712`
COMPUTED** (Gate 3), `λγ` a genuine model gap (`BETA_GENUINE_GAP`), and — **Gate 4 (2026-06-22) → `GENUINE_BLOCKED`** — the
gravity ratio `g_mhat²/g_G` is **not derivable** from the current action (no target-blind source-map kernel; `α_J` doesn't
cancel; all 22 `m̂` sites in `pde.tex` are target-facing, dual-reviewed). So the model does not derive its own gravity coupling:
`g_G` is calibrated on Newton's `G`, and the **`54/5` quadrupole is an ABSORBED calibration anchor, not a prediction**. The
verdict closes only with the **EM-sector anchor** (which pins `λγ`) — now load-bearing. The falsifiable payoff is the **held-out
surplus** (g−2, 5PN, ringdown, multi-defect), riding the shared *derived* `χ_Q` + `P0/D0` bundle + `c_s`.

## Next step

1. **The EM-sector anchor** — pin `λγ = c_γ/c_s` (and `μ0, q_*, Z_int`) against a charge/fine-structure-sector calibration; this
   closes the verdict count.
2. **The held-out surplus** — the actual predictions (g−2, 5PN, ringdown, multi-defect), no refit, on the shared bundle.

(Recovering `54/5` as a *prediction* would require ADDING a source-map/mass-bridge postulate to the action — a modeling choice,
not a computation. Out of scope unless elected.)

## Map — what you want → which doc holds it

| You want… | Look here |
|---|---|
| Full current state + resume-after-`/compact` pointer | `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0 |
| Every value classified DERIVED / INPUT / gap / benchmark | `software/stage1_solver/decisions/14_value_provenance_and_calibration_map.md` (Gate-4 closeout = §6) |
| The defect-regime + held-out-surplus roadmap | `docs/defect_interaction_map.md` (CURRENT STATUS banner at top) |
| Per-gate derivation detail (Gates 0–4) | `software/stage1_solver/reports/pathA_22b_minimal_combination_xi.md` |
| The directive that ran Gate 4 | `software/stage1_solver/directives/pathA_22b_gate4_ratio_or_blocked.md` |
| The calibrate-predict methodology | `software/stage1_solver/decisions/09_calibrate_predict_methodology.md`; `docs/pathA_preregistration.md` |

## The verdict equation (reference)

```
P0 · χ_Q · g_mhat² · λγ⁵ / g_G  =  54/5
 ✓     ✓     gap1     gap2  cal-on-G        (✓ = derived; gap1 g_mhat absorbs 54/5; gap2 λγ ← EM anchor)
G = (a·c_s²/m_GNLS)·g_G ,  m̂0 = (c_s/(a²·√m_GNLS))·g_mhat ,  c = λγ·c_s
```
