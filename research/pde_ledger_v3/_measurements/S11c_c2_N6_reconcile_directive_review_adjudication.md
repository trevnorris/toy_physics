# S11c-c2 N6 RECONCILE directive — decision-review adjudication (2026-09-06)

Orchestrator-written physics-bearing pre-builder directive → 2 decision legs (Codex-sol xhigh + Grok, identical
prompt `directives/_legs/S11c_c2_N6_reconcile_directive_review_prompt.md`). Both EXIT=0. Reports:
`scratchpad/{codex_N6rc_decision.log, grok_N6rc_decision.log}`. **Both legs: FOLD REQUIRED (neither CLEAR TO BUILD).**
Findings CONVERGENT (both independently caught F1); no conflict. My verification (G4) below; all VERIFIED, folded once
into `directives/S11c_c2_N6_reconcile_directive.md` (v2). Because the directive is physics-bearing, the folded v2 gets
a re-review pass (review-until-clear) before the astra build.

## Findings — verified + disposition

**A1 — the `a_ρ` corruption cannot move the CARRIER bridge (BOTH legs; F1).** VERIFIED. `a_ρ` is the constitutive
θ-advection inside material μ (`constitutive` :328-342); the open pressure-slot coefficients carry neither θ nor μ
(route-2 §4 :130-150), and `face_factory` never sees `a_ρ`. So dropping `a_ρ` moves the SOURCE (via μ_M), ⛔ not the
carrier `C_M`. My directive's requirement that it move `CARRIER_BRIDGE_RESIDUAL` too is fabrication-forcing (astra
would have to illicitly couple `a_ρ` into `C_M`), and on `RHO4_CONSTANT` (`g_i=0⇒a_ρ=0`) it is A−A. **Fold:** separate
one-sided corruptions — carrier bite = covector/normal-map corruption with `b_M` fixed; source bite = `a_ρ→0`,
`RHOBR_CONSTANT`-only, moving source/channels and leaving carrier + Eulerian byte-identical; **OBSERVE** the
end-to-end `R_N6` response, ⛔ do not require it to move (route-2 :209).

**A2 — the two-term split does not localize the cross-term (Codex F2).** VERIFIED (algebra reproduced):
`I(ΔC,S_E)=I(ΔC,S_M)+B(ΔC,ΔS)`, so my "carrier channel" `I(ΔC,S_E)` hides the carrier–source cross-term. **Fold:**
adopt the exact three-way split `R_N6 = I(ΔC,S_M) + B(C_M,ΔS) + B(ΔC,ΔS)` — `CARRIER_CHANNEL=I(ΔC,S_M)` (build_increment
on ΔC with the **material** source ms; carries the bare −ΔC·p), `SOURCE_CHANNEL=B(C_M,ΔS)` and `CROSS_CHANNEL=B(ΔC,ΔS)`
(both signatures-{6,9,12}-only contractions, no bare term). Localizes carrier-only / source-only / cross.

**A3 — pin `b` to the diagnostic's actual sources + in-process single `pit()` (Grok F2).** VERIFIED. The increment's
source is the diagnostic `es`/`ms` circuits (`source_terms` :378-389, :841-843), ⛔ not a re-coded `Λ` expression
(which would double-count ε via `V̄=V/ε` :389 or miss imported slot factors, and violates corollary 1). The
`b_{r,s}` formula is the slot IDENTIFICATION only. **Fold:** `SOURCE_BRIDGE_RESIDUAL := es − ms` (source_terms
circuits); `ΔS := es − ms`; recompute `E,M,R_N6` in-process; one `pit()` over the whole reconcile object dict
(including `SPLIT_CHECK`); ⛔ never join stored diagnostic PIT tables (a cross-run join makes `SPLIT_CHECK` vacuous).

**A4 — `SPLIT_CHECK` "structurally zero node" is unsatisfiable (Codex F3).** VERIFIED. `plus`/`minus` (:137,:163)
drop literal-zero operands only; they do not reduce `x−x`, so a faithful independent `SPLIT_CHECK` stays an `add`
node. **Fold:** the `SPLIT_CHECK` guard = **exact samplewise-zero numerators on the shared PIT samples** (every draw:
numerator ≡ 0), ⛔ not a structural-node-identity; ⛔ prohibit satisfying it by emitting `number(0)` directly.

**A5 — add the sanctioned jet-vocabulary bridge (Codex F4).** VERIFIED. Route-2 :128 mandates the explicit
`a.grad_theta[i]` (`theta_d1`) ↔ `b.grad_theta[i]` (`grad_theta_1`) boundary table; the diagnostic exposes it
(`N6_JET_BRIDGE` :351-356). **Fold:** add this exact identity to the predeclared/emitted frozen-relation set; ⛔ forbid
broader symbol renaming.

## Held sound (both legs, not folds)
The corrected zero-modulo-definitions question; absence of a nonzero offset `J`; the fixed-anchoring restriction; the
coefficient-level carrier test (`Increment(ΔC,S_E)=0` alone is weaker); the combined-source formula kept OUTSIDE
opaque `Z[b]`; the affine-split algebra + term-2 excluding signature 0; the one-sided PIT interpretation
(nonzero = certificate, all-zero = "no nonzero found" at conditional δ); no expected value / no residual-zero exit.

## Round 2 (re-review of v2) — 2 more folds (both legs; convergent; sharply converging 5→2)
Both re-review legs (Codex-sol + Grok, `scratchpad/{codex,grok}_N6rc_rereview.log`, EXIT=0) independently re-verified
all five round-1 folds landed faithfully, and **both** found the same 2 residual items → v3. Not defect-breeding
(R2 items are: an incomplete application of A4 I missed, and a carrier-bite target error round 1's A1 fix exposed).

**R2-F1 — stale `SPLIT_CHECK` "structurally zero" at directive:154 (BOTH legs).** VERIFIED: A4 fixed the contract at
the affine-split section but the verification bullet still said "structurally zero"; `plus`/`minus` (:137,:163) don't
reduce `x−x`, so a faithful `SPLIT_CHECK` stays an `add` node and the only way to satisfy "structurally zero" is
`number(0)` — which the directive forbids and which can't detect a wrong split. **Fold:** verify bullet = independent
residual node with samplewise-zero shared-PIT numerators (print node + per-sample numerators; ⛔ not `number(0)`).

**R2-F2 — carrier-bite target error at directive:117-120 (BOTH legs; Grok sharper).** VERIFIED: the carrier bridge is
`C_E − C_M` = imported `e_coeff` (:797) − material `m_coeff` (:809). The Eulerian-factory tilt (:808) writes
`tilt_coeff` and touches NEITHER — so it cannot move the bridge (my example was unsatisfiable, same A1 class). **Fold:**
carrier bite = corrupt the **material** covector/normal-map feeding `C_M` (`material_inverse_transpose`/material
normal, route-2 :113), `ms` held uncorrupted ⇒ moves `CARRIER_BRIDGE`/`CARRIER_CHANNEL`, leaves `C_E`/Eulerian/`SOURCE`
unchanged. The Eulerian-factory tilt is the SEPARATE reconstruction/independence probe, itself gated on
unmodified-factory = imported `C_E` (baseline nonzero = reconstruction drift; route-2 :195,:209).

**Round-2 held sound (both legs):** the corrected question, three-way split, `es`/`ms` pinning + in-process `pit()`,
jet bridge, coefficient-level carrier test, source outside `Z`, one-sided PIT reading, script clauses + corollaries.
⇒ v3 (both R2 folds applied) → one scoped round-3 confirmation before the astra build.
