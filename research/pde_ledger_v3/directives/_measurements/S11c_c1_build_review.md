# S11c-c1 engine — build-review record (2 legs, Codex-written ⇒ fresh Claude Agent + Grok)

Artifact: `scripts/S11c_c1_bulk_closure_sympy_audit.py` (Codex-written per the migrated build directive; runs
exit 0, writes the 906,467-byte / 44-row own-rows-delta `scripts/S11c_c1_exports.py`). Review prompt (identical
to both legs): `directives/_legs/S11c_c1_build_review.md`. Both legs wrote independent derivations BEFORE opening
the engine and did source-copy form ablations; raw reports saved outside the repo (tree hygiene). Both verdicts:
**CHANGES-REQUESTED.** The orchestrator verified every finding (rule 13); the ONE inter-leg disagreement is
adjudicated below with a script.

## Confirmed SOUND (both legs' independent derivations matched the engine)
- The **two-momentum DtN kernel** carries both branch legs `q_out(ω,k)`, `q_out(ω,k′)` LIVE (form ablation
  `q_input=q_output` MOVES it in both legs; a byte-identical output would have meant a freeze).
- The **permeable face response** is the operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹` (not a scalar division); `Λ_X`
  appears ONLY in the traction (`FLUX_V0_HAS_X False`, `TRACTION_V0_HAS_X True`); `μ_θ` is carried OPAQUE
  (no `theta`/`e_W`/`u` substitution).
- The **own-rows delta**: 44 rows (5 model-level + 39 new `s11cc1_` symbols, NOT the 2441-row model); fold adds
  it with ZERO overwrites; the bare `face_response` stays the S11b row, distinct from c1's F9c
  `s11c_c1_face_response` (predecessor preserved); `IMPORT_KEYS`=44 exactly the set `build_model` binds; the
  `check_consumer` + `assert_lookups_equal_manifest` lookup smoke-test **bites both directions** (drop a bound
  key → `ManifestError: undeclared lookup`; add a bogus key → `declared-but-unused`).
- The §5a-tilt one-sided independence, §5b per-direction `W_bg` form ablation, §5e branch sign-flip/momentum-
  freeze, and §6 homogeneity (source-level dimension corruption gives a nonzero residual) are GENUINE controls.

## ADJUDICATED disagreement — Grok "tangential kernel freeze" is a FALSE POSITIVE (Claude right)
Grok read `KERNEL_HAS_K_OUT False` on the first-shape term (`k_input_sq` only) as a rule-17 freeze of the output
tangential leg. Claude showed `k_input²` (η-order) + the tilt cross-term (σ_W-order) reconstruct the two-leg
`k·k'` form of `Div(h∇)` under the gradient relation `w1_jet_hat = i(k_out−k_in)·w1_hat`. **Verified**
(`_measurements/c1_tangential_kernel_adjudication.py`): `remainder − target ≡ 0` — the engine reconstructs the
correct two-momentum bilinear form; `k_out` enters via the profile hats' transfer argument (`w1_profile_hat_*` /
`w1_profile_jet_hat_*`), exactly the spec-mandated independent η/σ_W bookkeeping (§2c). The DtN kernel is CORRECT.

## Confirmed DEFECTS — all IMPLEMENTATION bugs in emit-only controls (the spec is correct); repair pending
1. **[both legs] The "independent" energy-balance route (§3b obj. 3) is A−A / toothless.** Lines ~1154–1169:
   `outgoing_flux = ½Re(iρω·[z_energy/(iρω)]·|V|²)` — the `(iρω)` cancels, so the "far-field flux" is just
   `z_energy` (= `δp/v` AT THE FACE, which §3b explicitly forbids: "⛔ not δp at the face"), and
   `bulk_comparison = −outgoing_flux = face_power` ⇒ `ENERGY_RESIDUAL[BASELINE] ≡ 0` for correct,
   branch-corrupted, AND arbitrary-garbage `z_energy`. A branch/impedance/`q_out`-sign error is invisible.
   **Fix:** compute the bulk operand independently as the outgoing plane-wave **Poynting flux from φ at |w|→∞**.
2. **[both legs] The representation-invariance control (§5a "the genuine control") is A−A.** The EULERIAN and
   HANZAWA routes are the SAME rational function (`route=` only regroups the algebra); `EULERIAN − HANZAWA ≡ 0`
   identically; grep finds NO coordinate change / layer potential / transformed radiation condition. **Fix:**
   build a genuine radiation-preserving Hanzawa / layer-potential second construction with its own outgoing kernel.
3. **[Claude] The on-shell `xreplace` silently misses** (lines ~911, ~1365): keyed on `Add(k_out_i²)`, but
   `sp.factor` distributes `c_s0²` so the bare `Add` never matches ⇒ the `DTN_RIGID_SHIFT_RESIDUAL` is emitted
   OFF-shell (nonzero-as-emitted) and the §5d ZERO_JET control carries the same off-shell term (reads as a
   spurious thickness-dependence *physics* finding). The kernel is physically correct on-shell. **Fix:**
   substitute the dispersion on `q_out²` directly (or scale/simplify before the match).
4. **[Grok] The Hermitian/two-port dissipation objects are `Z^(1)`-only.** `port_matrix` (line ~777) and
   `DTN_HERMITIAN_PART` (line ~970 `hermitian_kernel(direct_kernels)`) are fed the first-shape kernel only, so
   both VANISH at `(η,σ_W)=0` and cannot audit the leading bulk radiation `H_a[Z_0]` (§3b obj. 1) or the closed
   port map at leading order (§3b obj. 2). **Fix:** feed the full bare DtN `Z_0+Z_1` to the Hermitian/port objects.
5. **[Grok] `assert_delta_is_minimal` is theater.** `own_closure = resolved_keys` (the delta itself), so it
   asserts `delta == delta` and an accidental re-accumulation would pass. (The lookup smoke-test still validates
   the export, so this is a weakened, not absent, guard.) **Fix:** wire `own_closure` to the intended export
   key-set (the 5 model-level roots + their new-symbol closure), NOT the delta's actual keys.

## Disposition
Core physics 2-leg-confirmed sound; the exported model-level rows are unaffected by all 5 fixes (the defects are
in emit-only controls/guard). Fold ONCE into a repair directive (rule 7); the repaired engine gets 2 fresh legs.
