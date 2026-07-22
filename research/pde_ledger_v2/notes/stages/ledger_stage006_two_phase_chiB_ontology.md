# ledger_stage006_two_phase_chiB_ontology

## ⚠ Decision-16 amendment (2026-07-21) — α_aniso retired

`α_aniso` (the easy-plane anisotropy `α_aniso χ_B (P·ŵ)²`, pin P7) was retired with the polar field `P` (Decision 16,
`software/stage1_solver/decisions/16_retire_brane_polar_field.md`; → failures-paper backlog): removed from the live A1
dimensional surface (its absence asserted, with a reinjection tooth) and from the drift, so **operative `DRIFT(5)`** =
`{χ_B, a_B, κ_B, Γ_B, gating}`. Pin P7 is reframed: `χ_B` is simply **the** postulated wall order parameter (route
(a)), NOT built as `|P_∥|²`; `χ_B = |P_∥|²` (route (c)) stays a named, high-risk, Part-VII-adjacent future gate needing
a new T0 freeze (obsolete-as-carried, not foreclosed). Everything else stands (the recovery reduction, the
θ-as-Maxwell-φ no-go, C1–C4, all other teeth). Dual-engine after the amendment: **SymPy 121 PASS / Mathematica 119
PASS**, both exit 0.

## Status

`ACTION_SPECIFIED_CLASSIFIED` (structure) — the two-phase material-state ontology
made explicit as a **POSTULATED** action, dimensionally classified, with two
EARNED-within-closure pieces and one carried able-to-fail no-go.

**The honest one-line landing: the wall (brane) is made explicit as a postulated
field — it is NOT derived from the one medium.** The stage's job is labeling and
classification: every term of the χ_B free-energy functional and order-state
balance PDE carries a printed provenance class (POSTULATED / CONVENTION / CITED /
EARNED-relative-to-postulates / DEFERRED). Nothing postulated is presented as
derived, and cross-engine agreement is NOT evidence of derivation (both engines
encode the same postulated `F`; agreement proves consistent encoding only).

- **EARNED (relative to the imposed split):** (1) the structural-closure
  identities of the two-phase split (advective form, disorder complement,
  Γ_B-cancelling total-conservation sum); (2) the wall-admission kink check
  (single kink, see the slab caveat below); (3) the **recovery reduction** — at
  `chi_B = 1, Gamma_B = 0, J_chi = 0` the projected two-source law reduces
  exactly (assert-zero) to the frozen OLD-ledger projected-leakage law
  `S_leak` of `research/pde_ledger/paper/stages/stage_243.tex` /
  `stage_244.tex`, including the stage-244 Gaussian one-mode closed form
  re-derived by direct integration (`RECOVERY_REDUCTION_VERIFIED`).
- **Carried no-go (able-to-fail):** θ-as-Maxwell-φ is a `FATAL_FLAW` — a stable
  order-parameter phase has the WRONG-SIGN gradient stiffness for Maxwell's
  electric term; the Maxwell locus is reachable `BY_TUNING` only
  (`C5_RESOLVED_MAXWELL_BY_TUNING`), never `WITH_PROVENANCE`; the only
  provenance sign-flip is a Lifshitz instability = the already-killed pathA_25
  wall. Tokens: `FAIL_CAUCHY_STRAY_LONGITUDINAL` (finite `B_eff`),
  `FAIL_C5_LONGITUDINAL_ZERO_MODE` (`B_eff = 0`).
- **Does NOT earn light** (the no-go is precisely why), does NOT solve any
  wall/throat dynamics, does NOT assert the global returns locally.

Third stage of **Part I (The medium)** (I-3). **Consumes** `ledger_stage004`
(I-1: `{L,T,M}` dictionary + single-well `U`), `ledger_stage005` (I-2:
`c_s0^2 = 5 K rho0^4 / m_GNLS`), and `ledger_stage003` (Part III:
`c_gamma^2 = mu_R/rho_br`, `C_J = -J rho_B0`, `B_eff = rho_B0^2/chi_c`, and the
second-class-pair ⇒ stray-DOF classification rule) as cited inputs with
exact-value citation-integrity checks.

**FRESH-AUTHORED stage:** no prior script pair existed; the audit was authored
from prose sources (provenance below), which is why the directive SPECIFIES the
action rather than extracting a harness.

## Purpose

Formalize the program's brane/bulk ontology — ONE conserved medium, TWO material
states — as an executable, dual-engine-audited object: brane = the medium's
ordered, shear-supporting phase (`chi_B = 1`); bulk = the SAME medium
de-structured, shear-free (`chi_B = 0`); throat = a phase-conversion site
(admittance/outlet, not a suction force); gravity's drain = steady brane-order
loss through an open throat with the return as the separate global process.

## The postulated action and balance laws (all POSTULATED, labels printed)

Free energy (4D spatial `d^4X`, time separate; every integrand term audited to
the 4D free-energy density `M L^-2 T^-2`):

```
F = ∫ d^4X [ (1/2) m_GNLS n |u|^2 + U(n) + f_B(chi_B)
             + (kappa_B/2)|grad_4 chi_B|^2 + chi_B f_shear
             + f_throat + f_mix ]
```

Balance laws:

```
∂_t n + div_4(n u) = 0                                    (total, exact)
∂_t(chi_B n) + div_4(chi_B n u + J_chi) = n Gamma_B ,      [Gamma_B] = T^-1
Gamma_B = Gamma_return - Gamma_drain
```

`Gamma_B < 0` = brane→bulk de-structuring (drain); `Gamma_B > 0` = bulk→brane
re-ordering (return). Total material is conserved; brane-ORDER is not.

### The 13 pinned modeling choices (P1–P13)

| Pin | Choice |
|---|---|
| P1 | `n` = total conserved 4D constituent **number** density, `[n] = L^-4` (I-1 row); kinetic term carries the constituent mass: `(1/2) m_GNLS n |u|^2`. **`KINETIC_MASS_FACTOR_PINNED`** — the handoff's `(1/2) n|u|^2` restatement lacks the mass factor (see Corrections below). |
| P2 | `chi_B` dimensionless, range `[0,1]`; `n_B = chi_B n`. |
| P3 | Double-well `f_B(chi_B) = a_B chi_B^2 (1-chi_B)^2`, `a_B > 0`, minima `{0,1}` (recovery limit sits at a minimum), `[a_B] = M L^-2 T^-2`. POSTULATED — the cited I-1 `U(rho) = (K/4) rho^5` is SINGLE-well, so two-phase structure cannot come from the parent potential. |
| P4 | **`SAME_DENSITY_DEGENERACY_POSTULATED`** — `f_B` is n-independent, `f_B(n,0) = f_B(n,1) = 0` (asserted). A density-split (liquid–vapor) interface is rejected as the brane: no shear ⇒ no light (the pathA_25 lesson). |
| P5 | Interface gradient `(kappa_B/2)|grad_4 chi_B|^2`, `[kappa_B] = M T^-2`. |
| P6 | Shear gate `chi_B f_shear`, `f_shear = (1/2) mu_R^(4) (curl_4 u_d)^2` with **displacement** `u_d` (`[u_d] = L`, distinct from the velocity `u`), `[mu_R^(4)] = M L^-2 T^-2`. Rotational (MacCullagh) form, NOT `(∂u)^2` Cauchy. Brane-projected `mu_R` (`M L^-1 T^-2`, stage003) related by `∫ chi_B mu_R^(4) dw` — dim-consistency asserted, projection POSTULATED/PENDING. |
| P7 | **`chi_B` is THE postulated wall order parameter** — an independent scalar field (route (a)), NOT built as `|P_par|^2`. The anisotropy `alpha_aniso chi_B (P·w_hat)^2` (`POSTULATED_ANISOTROPY`, `[alpha_aniso] = M L^-2 T^-2`) is **RETIRED with `P` (Decision 16 → failures-paper backlog)**, absent from the live A1 surface. `FUTURE_GATE_CHI_B_EQ_ABS_P_PARALLEL_SQ` (route (c) — `chi_B = |P_par|^2`) remains a named high-risk Part-VII-adjacent gate needing a NEW T0 freeze; obsolete-as-carried, not foreclosed. |
| P8 | Dynamics = labeled adjunct: `D_t chi_B = -M_chi mu_chi + Gamma_B`, `mu_chi = δF/δchi_B` (`M L^-2 T^-2`), `[M_chi] = L^2 T M^-1`. Energy ledger pinned to the variationally consistent `P_order = ∫ mu_chi D_t chi_B d^4X` (= `M L^2 T^-3`, power). **`HANDOFF_P_ORDER_N_PLACEMENT_CORRECTED`** (see Corrections). |
| P9 | `J_chi = 0` default (simplest advective closure); `J_chi ≠ 0` deferred, `[J_chi] = L^-3 T^-1`. |
| P10 | Recovery target = the **frozen canonical OLD-ledger** `S_leak` (stage_243/244.tex), with `j^w ≡ n u^w` and unit-normalized `W(w)`: `[W] = L^-1`, `∫W dw = 1`. NOT the handoff restatement (which would be circular). |
| P11 | Global returns `R_0 = -M_0`, `R_1 = -D_1` (pathA_28/29) printed as labeled postulates, NOT asserted locally. |
| P12 | Throat = phase-conversion **ontology, not a solve**: the site of `Gamma_drain > 0` (`chi_B: 1→0`), an admittance/outlet driven by existing stress/μ gradients (`J_repair ~ -M_n grad mu`, `[M_n] = L^-4 T M^-1`) — NOT a new pairwise suction force. Gate-L wall-translation-Goldstone hazard (`δw = u_w`) DEFERRED. |
| P13 | Sign convention pinned: free-energy `+(1/2) kappa_phase (grad theta)^2` (stability ⇒ `kappa_phase > 0`) enters the Lagrangian as `K_theta = -kappa_phase < 0`; Maxwell's electric square needs `K_theta = C_J^2/rho_br > 0`. |

## Leg A — dimensional classification + structural closure

**A1 (dims).** Every OPERATIVE `F`-integrand term = `M L^-2 T^-2` exactly (exponent-triple
bookkeeping, both engines): kinetic `M·L^-4·L^2T^-2`; `K n^5 =
M L^18 T^-2 · L^-20`; `a_B` well; `kappa_B` gradient (`M T^-2 · L^-2`); shear
gate; placeholders `f_throat`/`f_mix` (`DEFERRED_PLACEHOLDER`). The `α_aniso`
anisotropy is **excluded from the operative surface (Decision 16)** — its absence
is asserted (reinjection tooth), and it is dim-audited only as a retired-historical
term (it *was* homogeneous, gone by decision not defect).
Balance rows all `L^-4 T^-1` ⇒ `[Gamma_B] = T^-1`. Adjunct rows: `mu_chi`,
`M_chi`, `P_order` (power), `M_n`, `J_repair`. Projection rows: `[W] = L^-1`,
`[rho_B] = L^-4`, `[S_leak] = L^-4 T^-1`.

**A2 (closure identities, EARNED, computed via real derivative expansion).**
From the two balance laws:

```
D_t chi_B ≡ ∂_t chi_B + u·grad_4 chi_B = Gamma_B - (1/n) div_4 J_chi
```

(advective form); with `n_D = (1-chi_B) n` the disorder complement
`∂_t n_D + div_4(n_D u - J_chi) = -n Gamma_B`; and the SUM of order + disorder
balances reproduces total conservation IDENTICALLY — `Gamma_B` cancels:
conversion moves order between phases, never creates medium. Order
non-conservation is genuine (`n Gamma_B` source nonzero for `Gamma_B ≠ 0`).

**A3 (wall admission, EARNED relative to P3/P5).** The static 1D EL equation
`kappa_B chi_B'' = f_B'(chi_B)` admits the kink between the two minima:

```
chi_B(w) = 1/(1 + e^{-w/delta}) = (1/2)(1 + tanh(w/(2 delta))),
delta = sqrt(kappa_B/(2 a_B))          (width SOLVED from the EL equation)
sigma_wall = ∫ kappa_B (chi_B')^2 dw = sqrt(2 a_B kappa_B)/6,
[sigma_wall] = M L^-1 T^-2             (energy per 3-area in 4D)
```

Both engines verify the EL residual by substitution and integrate `sigma_wall`
in closed form (the `.wl` independently enters via `DSolve`'s logistic).

**⚠ The slab caveat (carried limitation — wall admission ≠ slab stability).**
"Bulk on both sides" (required for `±w` charge, conceptual_foundation §2) makes
the brane a finite `chi_B = 1` **slab** bounded by `chi_B = 0` on both sides —
a kink–antikink pair — and the double-well alone provides NO mechanism selecting
the slab width (kink–antikink pairs generically attract). A3 is a single-kink
admission check only. The slab width `W_slab` is an un-earned input — the old
`Z/W/B_ell` profile scale in new clothes — and maps onto the OLD ledger's known
open item (`L/a` self-selection "requires dynamics", Gate-6 / sim-deferred).
Registered as `FREE-UNREDUCED` in the parameter register.

## Leg B — the recovery reduction (EARNED; the frozen-firewall assert-zero)

Projecting the order balance with `W(w)` (`rho_B = ∫W chi_B n dw`,
`J_B = ∫W(chi_B n u_xyz + J_chi,xyz) dw`) and integrating the w-divergence by
parts gives the two-source law:

```
∂_t rho_B + div_3 J_B = S_flux + S_convert
S_flux    = -[W(chi_B n u^w + J_chi^w)]_∂ + ∫ W'(w)(chi_B n u^w + J_chi^w) dw
S_convert = ∫ W n Gamma_B dw
```

The brane observer loses ordered material two ways: transport through `w`
(`S_flux`) and loss of the order STATE in place (`S_convert`).

**The reduction (assert-zero):** the audits assemble the general projected
object with `chi_B(w)`, `Gamma_B(w)`, `J_chi^w(w)` computationally LIVE on
explicit convergent profile families (SymPy: Gaussian + rational; Mathematica:
its own hyperbolic family), substitute `chi_B → 1, Gamma_B → 0, J_chi → 0`
IN-ENGINE, and assert-zero against an INDEPENDENTLY constructed frozen target
`S_leak = -[W j^w]_∂ + ∫W' j^w dw` (`j^w = n u^w`; stage_243.tex
eq. `stage243-leakage-identity`), pinned to an independently transcribed
closed-form value so identical textual corruption of both sides still fails.
EARNED **relative to the imposed chi_B split + declared `W(w)`**.

**Conditionality (the limit does real work — computed residuals):** `chi_B = c ≠ 1`
leaves `(c-1)·S_leak`; `Gamma_B ≠ 0` leaves the live `S_convert` integral;
`J_chi^w ≠ 0` leaves its flux terms. Each `expect_nonzero` on a residual computed
from the same general object.

**The stage-244 Gaussian anchor:** with `W_lambda = e^{-w^2/lambda^2}/(lambda
sqrt(pi))`, `phi_lambda = (2w/(sqrt(pi) lambda^3)) e^{-w^2/lambda^2}`,
`E_w = -E0 phi_lambda`, `j^w = mu_w rho0 E_w`, direct exact integration gives

```
S_leak = sqrt(2) mu_w rho0 E0 / (2 sqrt(pi) lambda^3)
```

matching `stage_244.tex` eq. `stage244-sleak-e0` — anchoring the reduction to the
OLD ledger's committed number, not to this stage's own algebra.

## Leg C — the θ-as-Maxwell-φ no-go (carried FATAL_FLAW, able-to-fail)

**Claim killed:** complex `chi_B = |chi_B| e^{i theta}` with the phase θ serving
as the MacCullagh/Maxwell scalar potential φ (the C5 resolution). Do NOT carry
into any rung-φ.

**Division of ownership:** `ledger_stage003` owns the full Dirac–Bergmann
machinery and the classification RULE (bracket ≠ 0 ⇒ second-class pair ⇒ 1 stray
longitudinal DOF ⇒ `FAIL_CAUCHY_STRAY_LONGITUDINAL`; `B_eff = 0` variant ⇒
`FAIL_C5_LONGITUDINAL_ZERO_MODE`; bracket = 0 ⇒ first-class ⇒ Maxwell). I-3
CONSUMES that rule and OWNS what the χ_B ontology itself fixes:

- **C1 — the sign lock.** Phase stability of the ordered state requires
  `+(1/2) kappa_phase (grad theta)^2` in `F` with `kappa_phase > 0` (the
  Josephson structure couples `theta_dot` to the conjugate order density, NOT
  `grad theta` to `∂_t u`). Under P13: `K_theta = -kappa_phase < 0` in `L`.
  Maxwell requires `K_theta = C_J^2/rho_br > 0` (consumed `C_J = -J rho_B0`).
  The two required signs are OPPOSITE — asserted as an exact sign predicate.
- **C2 — the discriminator bracket.** For the longitudinal `(u_L, theta)` probe
  (pathA_36 primitive, with `m_theta^2 theta^2` algebraically present in both
  engines): primary constraint `Phi_1 = pi_theta - J k rho_B0 u_L`, bracket

  ```
  {Phi_1, Phi_2} = k^2 (J^2 rho_B0^2 + kappa_phase rho_br) / rho_br
  ```

  Its zero locus, SOLVED in-engine (not hardcoded), is exactly
  `K_theta = +J^2 rho_B0^2 / rho_br` — the only route to first-class is the
  tuned electric sign. Tokens are computed from wired discriminators
  (bracket-zero?, `B_eff`, `m_theta^2`, provenance flags): provenance branch →
  `FAIL_CAUCHY_STRAY_LONGITUDINAL`; `B_eff = 0` →
  `FAIL_C5_LONGITUDINAL_ZERO_MODE`; tuned fixture (`K_theta = J^2 rho_B0^2/
  rho_br`, `B_eff = 0`, `m_theta^2 = 0`) → bracket COMPUTES to zero →
  `C5_RESOLVED_MAXWELL_BY_TUNING`, flagged **BY_TUNING NOT WITH_PROVENANCE**
  (the frozen χ_B definitions force none of the three tunings; the flag is
  computed from provenance booleans with a live counterfactual). Detunings
  (K_theta off-locus / `B_eff ≠ 0` / `m_theta^2 ≠ 0`) each re-fire the
  corresponding pathA_36 control token (`FAIL_C5_LONGITUDINAL_ZERO_MODE`,
  `FAIL_SECOND_CLASS_NOT_MAXWELL`).
- **C3 — the Lifshitz identification.** The only provenance route to
  `K_theta > 0` is `kappa_phase < 0`: with stabilizing `kappa_4 > 0`
  (`[kappa_4] = M L^3 T^-2`) the phase sector `f(k) = (1/2)(kappa_phase k^2 +
  kappa_4 k^4)|theta_k|^2` has its minimum at finite
  `k*^2 = -kappa_phase/(2 kappa_4) > 0` ⟺ `kappa_phase < 0` (asserted both
  directions) — a modulated/Lifshitz instability of the uniform phase, i.e. the
  pathA_25 density-smectic wall (CITED, already falsified
  `FAIL_LIGHT_STARVED`; its `k_Rstar^2 = 40 K m rho0^4/hbar^2` dim-checked
  `L^-2`). The sign-flip escape re-enters a killed wall.
- **C4 — transverse control.** The transverse dispersion is derived in-engine
  from the ε-parametrized `L_T = (1/2) eps u_T_dot^2 - (1/2) mu_R k^2 u_T^2`:
  `omega^2 = (mu_R/eps) k^2`. Baseline `eps = rho_br` reproduces the consumed
  `c_gamma^2 = mu_R/rho_br` → `PASS_TRANSVERSE_UNDISTURBED`; the `eps = 2
  rho_br` fixture shifts the speed → `FAIL_TRANSVERSE_DISTURBED`. Tokens wired
  from the computed comparison (they genuinely flip with ε).

## Drift, consumed inputs, carried tokens

**Operative `DRIFT(5)`, computed from the enumerated list:** `{chi_B (field); a_B; kappa_B;
Gamma_B (conversion law); the chi_B-gating structural choice}` (`alpha_aniso` retired with `P`, Decision 16 →
failures-paper backlog; teeth: re-inject `alpha_aniso` → computed n=6 ≠ 5 fires; corrupt n → token equality fires).
Reconciliation vs `rung_W_reframe.md:140` printed in-script (rung_W counts a
2-constant generic well + no `Gamma_B`; the P3 one-constant `[0,1]` well + the
live conversion law reconcile to the operative `DRIFT(5)`). `M_chi`/`J_chi` are deferred adjuncts,
not live knobs. Dead θ-branch symbols (`theta, J, rho_B0, K_theta/kappa_phase,
chi_c, B`) are NOT admitted as live knobs (`THETA_BRANCH_DEAD_NOT_ADMITTED`).
Cross-reference: `rho_B0, chi_c` already appear in pathA_41's Part-VI drift trio
— the register must not double-count. This drift is the honest price of route
(a) — "add a postulated order field" — and is exactly the `SECOND_MEDIUM_DRIFT`
lineage the program tracks.

**Consumed (cited, exact-value citation-integrity checked, never re-derived):**
- `ledger_stage004`: `[n] = [rho0] = L^-4`, `[m_GNLS] = M`, `[K] = M L^18 T^-2`,
  `[hbar] = M L^2 T^-1`; `U(rho) = (K/4) rho^5` single-well (the P3 justification).
- `ledger_stage005`: `c_s0^2 = 5 K rho0^4 / m_GNLS`.
- `ledger_stage003`: `c_gamma^2 = mu_R/rho_br`; `C_J = -J rho_B0`;
  `B_eff = rho_B0^2/chi_c`; dims of `{mu_R, rho_br, rho_B0, chi_c, B, J,
  K_theta}`; the constraint-classification rule.

**Carried / deferred (none repaired):** the three no-go tokens; the pathA_25
`FAIL_LIGHT_STARVED` lineage; the `SECOND_MEDIUM_DRIFT` lineage note; the slab
caveat (above); Gate-L `δw = u_w` translation Goldstone; `J_chi ≠ 0`;
`f_throat`/`f_mix` placeholders; the dynamics/energy-ledger adjunct status;
global returns P11.

## Corrections to the source prose (pinned + ablation-documented, not silent)

Two dimensional defects in `notes/brane_bulk_handoff.md` were caught during
authoring (confirmed independently by the Codex design review):

1. **Kinetic mass factor (P1):** the handoff's `(1/2) n|u|^2` is inhomogeneous
   for a number density `[n] = L^-4`; the pinned form is `(1/2) m_GNLS n|u|^2`.
   The drop-`m_GNLS` firewall ablation fires on the handoff's form.
2. **`P_order` n-placement (P8):** the handoff's `P_order = ∫ mu_chi n Gamma_B
   d^4X` is inhomogeneous under its own `mu_chi = δF/δchi_B` (extra `L^-4`);
   the variationally consistent `P_order = ∫ mu_chi D_t chi_B d^4X` is power.
   The handoff's form is a firewall ablation that must fire.

An erratum pointer is placed in the handoff note; the corrections live here and
in the audits, the source is not silently rewritten.

## Verification

- **⭐ Decision-16 amendment (2026-07-21):** the `α_aniso` retirement applied per
  `_scratch/decision16_amendment_directive.md` (directive cleared the Codex→Grok→Codex bookend). Dual-engine after
  the amendment: **SymPy 121 PASS / Mathematica 119 PASS**, both exit 0; transcripts regenerated. ⏳ Fresh-agent
  tri-review of the amended scripts + docs is the next gate. The original-build record below is retained (tallies
  110/108 superseded).
- SymPy audit: `scripts/ledger_stage006_two_phase_chiB_ontology_sympy_audit.py`
  — exit 0, **121 PASS / 0 FAIL** (was 110; transcript in `scripts/output/`).
- Mathematica audit:
  `mathematica/ledger_stage006_two_phase_chiB_ontology_mathematica_audit.wl` —
  exit 0, **119 PASS / 0 FAIL** (was 108), genuinely independent route (own exponent
  association, formal-function closure expansion, `DSolve` kink entry, own
  hyperbolic B1 family, Lagrangian/determinant leg-C route vs the `.py`'s
  Hamiltonian/bracket route; transcript in `mathematica/output/`).
- Full tri-review: orchestrator arbiter re-run (runners, unchanged scripts);
  fidelity leg caught the v1 B2 recovery check as an X≡X tautology
  (`FIDELITY_ISSUES`); adversarial-with-ablation leg PROVED it (identical
  corruption of both sides survived at exit 0) plus stamped C4 tokens
  (`ADVERSARIAL_CONCERNS`, 13-mutation matrix, 11 genuine teeth). Remediated
  (general projected law assembled live + limit in-engine + independent frozen
  target; leg-C discriminators wired from algebra + `m_theta^2` made real in
  the `.py`; C4 ε-dispersion derived in-engine; `.py` A2 real derivative
  expansion; `.wl` B1 self-computed profiles) → fresh-agent re-verify
  `REVERIFY_CLEAN` with a 6-run corruption matrix (the previously-passing
  dual-corruption class now fails in both engines) + 2 spot-ablations.

## Provenance

Sources (READ-ONLY): `notes/brane_bulk_handoff.md` (the χ_B closure; Item A /
Tasks A–E; Tests 1–3), `docs/conceptual_foundation.md` §2–§4 (the physical
picture; single-well obstacle; MacCullagh-vs-Cauchy), `notes/rung_W_reframe.md`
(double-well forms; wall stability; the DRIFT tally; the θ FATAL_FLAW
structural argument), `software/stage1_solver/reports/pathA_36_c5_phase_potential.md`
(the computed no-go: probe Lagrangian, constraint pair, bracket, tuned locus,
controls — reshaped as `ledger_stage003`), and the frozen recovery targets
`research/pde_ledger/paper/stages/stage_243.tex` + `stage_244.tex`. Directive +
review + remediation artifacts: `_scratch/ledger_stage006_*`.
