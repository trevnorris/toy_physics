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
