# Projected Maxwell Mouth-Taylor Extension — Master Note

This bundle continues the projected-Maxwell / moving-throat bridge by replacing the
abstract primitive slippages with a concrete near-mouth Taylor ansatz.

## Files
- `step_12_projected_maxwell_mouth_taylor_gate_bridge_notes.md`
- `step_11_projected_maxwell_mouth_taylor_master_sympy.py`
- `step_12_projected_maxwell_mouth_taylor_gate_bridge_sympy.py`

## What is new
1. A concrete one-sided mouth-kernel ansatz
   \[
   X_{\rm proj}(\ell)=X(0)+\ell \mu_1 X'(0)+O(\ell^2)
   \]
   for the primitive one-port packet
   \[
   X\in\{Q,S_2,H_{\rm port},\Delta,P,G_W\}.
   \]

2. Exact derivative-level formulas for
   \[
   z_0,z_2,z_4,n_0,n_2,n_4.
   \]

3. Exact projected-Maxwell contributions to the live weak-axisymmetric / 5PN
   bottlenecks
   \[
   \Xi_{\rm load}=\frac{P_1}{P_0},\qquad
   K_1,\qquad
   H_{\rm even}.
   \]

4. An exact mechanism sieve showing:
   - pure \((Q',\Delta')\) cannot close the even gates nontrivially,
   - pure \((S_2',H'_{\rm port})\) cannot close them nontrivially,
   - a nontrivial repair requires a mixed compensation surface.

5. A clean transport-factorization:
   - \((Q',\Delta',S_2',H'_{\rm port})\) repair the conservative even gates,
   - \(P'\) tunes the prefactor slope \(\Xi_1\),
   - \(G_W'\) first tunes \(\delta P_2,\delta P_4\).

## Why this matters
This is the first concrete near-throat projected-Maxwell ansatz that feeds
directly into the same bottlenecks the moving-throat program now treats as live.
It strongly suggests that the missing near-throat EM piece, if real, should be a
multi-channel mouth correction rather than a single effective scalar.
