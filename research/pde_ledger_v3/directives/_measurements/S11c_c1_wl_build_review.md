# _measurements — S11c_c1 WL engine build review (2 legs) + orchestrator adjudication

The WL engine `mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` is **Codex-written**, so its build legs are
**a fresh Claude Agent + Grok** (authorship table). Both ablate Mathematica kernels, so they were **serialized**
(2-seat licence): fresh-Claude first (in-process Agent), then Grok (detached). Leg prompt:
`directives/_legs/S11c_c1_wl_build_review.md`. Raw Grok report: `…/scratchpad/grok_wl_build_review.log`; leg-1
report returned inline. Leg saved-artifact paths are listed in each leg's report.

## Commands (literal)
```bash
# Build (Codex, detached; setsid+marker+Monitor; danger-full-access for Mathematica)
codex exec -c model_reasoning_effort=xhigh --sandbox danger-full-access "$(<…/S11c_c1_wl_build_directive.md)" \
  > …/scratchpad/codex_wl_build.log 2>&1            # 603,683 build tokens
# Leg 1 — fresh Claude Agent (in-process, general-purpose; no shell command), serialized FIRST
# Leg 2 — Grok, serialized SECOND (both ablate Mathematica kernels; 2-seat licence)
grok --prompt-file …/_legs/S11c_c1_wl_build_review.md --cwd /var/projects/toy_physics \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain \
  > …/scratchpad/grok_wl_build_review.log 2>&1
```
> ⚠ Note (repair-2 gate, 2026-09-04): Grok's NIT below — "`PERMEABLE_PORT_HERMITIAN` parity keys … both get the
> same per-face matrix" — correctly identified the **fake-parity-axis** pattern; repair-1 fixed it in
> `PERMEABLE_PORT_HERMITIAN` (to the correct congruence, distinct blocks), but the **identical pattern in
> `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`** was outside that NIT's scope and is what repair-2 R2 fixes.

## Deliverable verification (orchestrator)
- `.wl` = 79,824 bytes / 1,862 lines; builder §8 report: 51 tags (50 shared + 1 local), all §3–§6 tasks ran,
  none skipped/deferred; `muTheta` one opaque operand; 603,683 build tokens.
- Blindness: structural scan clean (0 `Get`/`Import`/`<<`/`ReadString`/`OpenRead`/abs-path); builder's isolated
  vs in-repo runs byte-identical (`bebe8f6f…`, 59.4 MB, exactly 51 tag lines); `mathematica/out/` clean.

## Verdicts (the two legs DISAGREED — the gate working, rule 6)
- **Leg 1 (fresh Claude):** "SOUND — 0 MUST-FIX, 2 NITs." Independent SymPy derivation + Mathematica ablation of
  the engine's own functions. Confirmed the two-momentum kernel (both legs live under all 4 corruptions),
  rigid-shift cancellation, operator inverse, T-a re-derivation, every §5 control biting, blindness, reserved
  spellings, `μ_θ` opacity. Its 2 NITs (far-field phase modulus; formal operator heads) — the FIRST it
  under-weighted.
- **Grok:** "NOT sound — 3 MUST-FIX." Grounded every finding in literal ablation stdout
  (`/tmp/s11cc1_wl_grok_review/`). Confirmed the SAME sound core, but found three real defects in the
  operator-composition emit and the energy dissipation audit.

## Orchestrator rule-13 verification — I read the engine code myself and CONFIRMED all 3 Grok MUST-FIX
- **MF1 — `DTN_OPERATOR` composition freezes the input leg** (`operatorCompositionFromDerivation`, `.wl:542-560`).
  `gZero = FourierMultiplier[I qOut[omega, momentumOutput]]` is used in BOTH positions of
  `OperatorComposition[gZero, multiplication, gZero]`; `momentumInput = {kPrime…}` appears nowhere in the
  composition. The rightmost `N_0` (which acts on the input field) is frozen to the OUTPUT momentum — a rule-17
  WKB/left-quantized freeze in the operator emit, internally inconsistent with the verified two-leg `DTN_KERNEL`
  (`.wl:478-479` `qOutputLive`/`qInputLive` both present). Spec §3a (`SHARED_PHYSICS.md:247-261`): the operator
  must carry BOTH legs.
- **MF2 — energy "bulk far-field" operand is the face quantity × a free `farFieldPhase`, not φ at `|w|→∞`**
  (`deriveEnergyOperands`, `.wl:994-1028`). `bulkAmplitude` is solved from the FACE Neumann law
  (`I qOut·A == energyVelocity`, `:1008-1010`); `farFieldPressure = pressureAtFace·farFieldPhase` (`:1015`),
  `farFieldVelocity = energyVelocity·farFieldPhase` (`:1017`), so
  `outgoingFlux = measure·Re[pressureAtFace·Conj[energyVelocity]·|farFieldPhase|²]/2`. At `farFieldPhase=1` the
  bulk operand = the face pairing byte-identical → residual 0 (Grok ablation A12). That is exactly the
  `½Re(δp·V*)` spec §3b (`:321-330`) forbids — the two "independent routes" are one. Same defect class as the
  SymPy engine's build-review finding #1 (repaired R1: real far-field Poynting from φ at `|w|→∞`).
- **MF3 — energy face operand uses a reconstructed impermeable traction, not the response's `t_s`** (`.wl:1011-1012`).
  `tractionVector = tractionOrientation·pressureAtFace·normal` (fresh flat solve, NO `Λ_X`); it never binds
  `FACE_RESPONSE`'s `t_s` (which carries `Λ_X`, `RESPONSE_TRACTION_HAS_LAMBDAX=True`). So a `t_s`-sign or
  `Λ_X`-placement error in the closed response is invisible to the audit — contradicting spec §3b
  (`:329-330`: "a `t_s`-sign … error is caught here, in c1, by the traction-vs-far-field comparison"). Flipping
  the reconstruction's own orientation moves it (what leg 1 saw), but flipping the EMITTED `t_s` cannot.

## What both legs confirm SOUND (2-leg)
Two-momentum `DTN_KERNEL` (both `qOut(k)`,`qOut(k′)`, `Ŵ_bg(k−k′)`; all 4 corruptions move it; matches independent
derivation `iωρ_m(q_k q_{k′}+k·k′−κ²)/(q_k q_{k′})`); rigid-shift cancellation (drop shifted-trace → nonzero);
operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹` with both Fredholm legs; T-a re-derived from the level set + `s(n̂·ŵ)>0`
(tilt ablation moves it); §5b/§5d/§5e controls bite; §5d operand B = bare `Z_0` (no `W_0→W̄₀(1+η)`); blindness;
reserved spellings (`WBg`/`w1Profile`/`etaBg`/`LambdaA0`/`rhoBr`≠`rhoBrBgRho4Constant`/…); `μ_θ` opaque; §4 emit
complete; no `VERDICT`/native-boolean residual.

## NITs to fold into the repair (both legs)
- Leg1-NIT-1 far-field phase modulus — subsumed by MF2.
- Grok: `SOURCE_EQUATIONS` emit `HoldForm[{flatBulk$NNNNN,…}]` Module-local gensyms (unevaluable; the real Solve
  equations were dropped from the payload). `PERMEABLE_PORT_HERMITIAN` parity keys `DELTA_W`/`ZETA_C` looped but
  unused (both get the same per-face matrix — the δW/ζ_c coupling isn't computed). `§5a` layer-potential emits a
  typed head `RadiationPreservingLayerPotential[…]` with no data dependence (the biting second route is computed
  separately — cosmetic). Locus `REAL_ADMISSIBLE` labels any non-`True`/`False` `Reduce` result `ADMISSIBLE`
  (`Head=!=Reduce`) — a weak admissibility test. Leg1-NIT-2 formal operator heads for `DTN_HERMITIAN/REACTIVE`
  are defensible (true-area adjoint, reduces at retained order) — NOT repaired.

## Decision → REPAIR round (mirrors the SymPy engine's repair)
Core physics 2-leg-sound; 3 MUST-FIX all in the operator-composition emit + energy audit. Repair directive
`directives/S11c_c1_wl_repair_directive.md` (fix MF1 momentum labels; MF2 real far-field Poynting from φ at
`|w|→∞`; MF3 bind the response `t_s`; fold the biting NITs), then a repair build + 2 fresh re-review legs that
ablate each repaired control to confirm it now BITES. Reviewed baseline committed before the repair runs.
