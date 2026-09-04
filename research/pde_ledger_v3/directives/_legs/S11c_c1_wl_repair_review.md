# Independent re-review — S11c-c1 Wolfram engine REPAIR (build leg)

## Artifact
`research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` — the **repaired** blind Wolfram
engine (repair directive `directives/S11c_c1_wl_repair_directive.md`, R1–R4; baseline before repair `e139bc61`).
Your job: confirm each repaired control now genuinely **BITES**, that the fixes are real (not relabelled), and
that the 2-leg-sound core is unchanged. A prose claim is worth nothing — **show your own ablation script and its
literal stdout** or the claim is discarded (rule 2).

## What you are handed
- The artifact above; the repair directive `directives/S11c_c1_wl_repair_directive.md`; the build-review record
  `directives/_measurements/S11c_c1_wl_build_review.md` (the three original defects, code-verified).
- The sole physics authority `directives/S11c_c1_SHARED_PHYSICS.md` (§3a `:247-261`, §3b `:321-330`, §6) + the
  siblings `S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`.
- ⛔ NOT the SymPy engine/export — do not seek them (this is a per-engine re-review; cross-engine is T7's job).

## What each repaired control must now do (ablate it — literal before/after)
**R1 — `DTN_OPERATOR` carries BOTH legs.** The composition `N_0 ∘ M_h ∘ N_0` (and the `N_0^{-1}` sandwich) must
have the **rightmost** factor at `momentumInput` and the **leftmost** at `momentumOutput`. Probe the emitted
`DTN_OPERATOR`: `COMPOSITION_HAS_Q_INPUT` and `COMPOSITION_HAS_Q_OUTPUT` must both be `True`; kernel-expand the
composition on-shell and difference against `DTN_KERNEL` — the correct labelling gives 0, and **re-freezing the
rightmost factor to the output leg (a control) makes it nonzero**. Confirm `DTN_KERNEL` itself is unchanged from
`e139bc61` (still two-leg).

**R2 — energy BULK operand is a genuine far-field Poynting flux of `φ`.** Confirm `farFieldPhase` is **gone** and
the bulk operand is `½Re[δp·v_bulk,n*]` of the half-space **outgoing** solution `φ` on a surface at `|w|→∞`, ⛔ not
the face `½Re(δp·V*)`. **One-sided corruption (the decisive test):** flip ONLY the bulk route (a `q_out` branch
flip in the far-field outgoing solution) — the bulk operand and `ENERGY_RESIDUAL` must move, the **face** operand
must NOT. And there must be **no** substitution (`farFieldPhase→1` or otherwise) that collapses the residual to a
structural zero. If flipping the bulk route also moves the face operand, or the two operands are byte-identical at
any parameter, the routes are still not independent — a finding.

**R3 — energy FACE operand binds the response `t_s`.** Confirm the face operand is built from the closed
response's `t_s = −(δp_s + Λ_X 𝒜_s)n̂_s` (`FACE_HAS_TS_FROM_RESPONSE=True`), ⛔ not a fresh impermeable solve.
**One-sided corruption:** flip the sign of the **response** `t_s` at its source — the face operand and
`ENERGY_RESIDUAL` must move, the bulk operand must NOT. If flipping the response `t_s` leaves the residual
unchanged, the audit still cannot catch a response `t_s` error — a finding.

**R4 NITs:** `SOURCE_EQUATIONS` payload holds the **real** re-parseable equations (no `flatBulk$NNNNN` gensyms);
`PERMEABLE_PORT_HERMITIAN` `DELTA_W`/`ZETA_C` carry actually-computed blocks (or a computed equality if they
coincide), not a duplicated matrix; the `§5a` layer-potential operand is the computed second-route kernel (no
decorative `RadiationPreservingLayerPotential[…]` head with no data dependence); locus `REAL_ADMISSIBLE` comes
from a genuine admissibility test (an undecidable branch is `UNDECIDED`, ⛔ not auto-`ADMISSIBLE`).

## Regression — the sound core must be byte-identical to `e139bc61`
`git diff e139bc61 -- research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` should touch
ONLY the composition, the energy operands, and the four NIT emits. Confirm `DTN_KERNEL`, `DTN_FLAT_SYMBOL`,
`DTN_RIGID_SHIFT_*`, the permeable `FACE_RESPONSE` operator inverse, the T-a..T-i re-derivation, the §5b/§5d/§5e
controls, reserved-name spellings, and `μ_θ` opacity are unchanged and still bite (spot-check one). Blindness
stands: 0 `Get`/`Import`/`<<`/`ReadString`/`OpenRead`/abs-path.

## Method (SCRIPT — derive, then ablate on a /tmp COPY)
Write your own check of the load-bearing objects first (save script + literal stdout to named absolute paths).
Then ablate on a **/tmp copy** and report literal before/after. A **one-sided corruption** (break ONE route, show
only its operand moves) is the only test that the two routes are independent — a zero residual alone proves
nothing.
```
⛔ Wrap EVERY kernel run in `timeout 600` (a 600s hit is a FAILED ablation — report and move on).
⛔ NEVER raise the timeout; run AT MOST ONE kernel at a time (2-seat licence — another leg may run after you).
⛔ Copy the .wl to /tmp and ablate the COPY; ⛔ never modify the working tree.
⭐ Save every ablation script AND its literal stdout to named absolute paths, and report them.
⚠ If a background job is killed with a healthy log, check `free -h` and `ps … --sort=-rss | head` FIRST (orphan kernel).
```

## Output
For each of R1/R2/R3/R4: state whether it now BITES, with the literal one-sided-corruption before/after that
establishes it. Report any control that still does NOT bite as a MUST-FIX with its literal output. State whether
the core regression is clean. If everything bites and the core is unchanged, say so and name the two or three you
ablated most closely. Report only findings that catch a way the physics or a control could be wrong.
