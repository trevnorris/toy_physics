# Measurements twin — S11c_a_wl_engine_repair2_directive.md

## Legs (Grok found both; orchestrator verified at src)
- Grok leg2: log ~/.s11_build/S11c_a_wl_engine_repaired_grok.log; scripts /tmp/s11ca_wl/{derive_s11ca,sliced_tf,sliced_baseline,...}.*
- Agent leg1 (ad986d4578bc1f2ca): CLEARED §5a/T-h but MISSED T-f + T-0-grad.

## F1 (T-f) verified at src (baseline a15bc69c)
CMD: sed -n '656,685p' .../S11c_a_interface_geometry_mathematica_audit.wl ⇒ projectionTermsSource keeps window INSIDE Inactive[Integrate][window·field,{w,-∞,∞}] (correct base).
  shapeDeriv 933-972 + Expand/Series pulls window'' OUT → Inactive[Integrate][1,{w,-∞,∞}] × window''(fixed) = divergent ∫1 dw. Grok stdout /tmp/s11ca_wl/sliced_tf.stdout.txt.
CMD: SymPy projection_terms (1002+) uses sp.Integral(window_bg·rho_shape_t + rho0·delta_window_t, (w,-oo,oo)) with WINDOW jets (flat_window_jets 979-999) as w-functions INSIDE ⇒ convergent. SymPy CLEAN.

## F2 (T-0 grad) verified at src
CMD: sed -n '705,720p' .../wl ⇒ `gradient = D[sigmaZero]/. etaBg W0/LW -> sigmaW` — no-op for RHO4_CONSTANT (∇Σ=η rhoBr/LW ∂w1, pattern etaBg·W0/LW absent) ⇒ graded {0,1,0} not {0,0,1}.
CMD: SymPy build_background_density_raw (897) ⇒ gradient=[dx(sigma_e,i)] via grad_W (σ_W jet)+finalize/multigrade → (ρbr/W0)σ_W w1_grad, σ_W-graded. SymPy CLEAN (round-2 legs confirmed).

## BLINDNESS: directive quotes Grok's finding (spec-derived) + spec §2a/§3c. NO SymPy construction referenced. Fix = "keep window inside integral" (math fact) + "grade via §2a map" (spec). Blind-safe.

## Round-3 directive review — Grok leg (rule 7 second leg; Codex was leg 1 bks19z3r1)
CMD: grok --prompt-file .../_legs/S11c_a_wl_repair2_review_prompt.md ... ⇒ log ~/.s11_build/grok_repair2_directive_review.txt (7656 bytes, exit 0; 2 prior launches spurious-killed at 509/715B pre-CAS).
3 CAS-backed edits, ALL verified by orchestrator & folded:
- A (Defect-2 recipe was a NO-OP as written): "retain W_bg" insufficient — ρ4,ref=ρbr/W0 and W_bg=W0(1+ηw1) self-cancel W0 BEFORE D, so no ∂W_bg unit for the §2a map. Fix = keep W_bg an UNEXPANDED/HELD factor through D, THEN map ∂_{y_i}W_bg → σ_W ∂_ξ w1.
- B (rule-5 leak I'd left in a verbatim quote): grade numerics {0,1,0}/{0,0,1}/{0,0,0} handed the builder the target ⇒ iterate-to-green pastes it. Stripped; kept qualitative "η zero-jet, not σ_W first-jet".
- C (confirm too weak): "no Inactive[Integrate][1,…]" passes (a) field left outside a non-1 integrand, (b) dropping mixed-η terms. Broadened to: no normalCoordinate-dep factor outside any held ∫dw + dynamic operand must still depend on etaBg inside the integrand.

## MY independent CAS verification of edit A (rule 13 / rule 2)
CMD: wolframscript -file .../scratchpad/s11ca_probe/t0_grade_probe.wl ⇒ literal stdout:
  BAD_GRAD1_RHO4: (etaBg*rhoBr*Derivative[1,0,0][w1p][...])/LW ; BAD_GRADE_RHO4 (etaBg,sigmaW): {1,0}  ← current recipe = η zero-jet (defect reproduced)
  GOOD_GRAD1_RHO4: (rhoBr*sigmaW*dw1[1])/W0 ; GOOD_GRADE_RHO4 (etaBg,sigmaW): {0,1}  ← held-factor recipe = σ_W first-jet ✓
  RHOBR_GRAD: {0,0,0} ✓
⇒ Grok's mechanism confirmed exactly; the held-factor fix produces the σ_W first-jet, the "Together/retain" fix does not.

## ⚠ POST-BUILD VERIFICATION TRAP (orchestrator-side, do NOT put in builder directive)
RHOBR gradient is the ZERO VECTOR; this engine's gradeTerms[0] = Exponent[0,·] = {{0,-∞,-∞}}, NOT {{0,0,0}}.
When verifying item (4): check the RHO4 gradient EXPRESSION is σ_W-graded (a first-jet in σ_W) and RHOBR's EXPRESSION is 0 — do NOT reject on a {{0,-∞,-∞}} multigrade for the zero vector. (Grok, run3.out baseline emit.)
