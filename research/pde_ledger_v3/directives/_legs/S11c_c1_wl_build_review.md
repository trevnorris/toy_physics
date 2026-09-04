# Independent physics review — S11c-c1 Wolfram engine (build leg)

## Artifact
`research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` — the **blind Wolfram engine** for
sub-step S11c-c1 (curved-interface bulk closure), written by Codex. It imports nothing and re-derives every object
from the shared-physics specs. Its product is a flushed stdout tag stream (`WL_S11CC1_<QUANTITY>`).

You are one of two independent legs. Your job is to find every way this engine's PHYSICS or CONTROLS could be
wrong — a conclusion emitted without a computation behind it, a control that cannot fail, a frozen operator, a
mis-derived object. A prose re-derivation is worth nothing; **show your own script and its literal stdout, or the
claim is discarded** (`CLAUDE.md` rule 2).

## What you are handed
- The artifact above.
- The **sole physics authority**: `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` (§§0–8), and the
  sibling specs it re-derives from: `directives/S11c_a_SHARED_PHYSICS.md` (§§1–4, the T-a..T-i substrate) and
  `directives/S11b_SHARED_PHYSICS.md` (§1b/§2 bulk acoustics + branch `q_out`; §9 B0 the flat reference).
- The build directive for context: `directives/S11c_c1_wl_build_directive.md`.
- ⛔ You are NOT handed the SymPy engine or its export, and you must not seek them: this leg checks the WL engine's
  own correctness + spec-fidelity + whether its controls BITE — the cross-engine comparison is T7's separate job.

## What to check (derive each yourself from the spec FIRST, then test the engine)
1. **The DtN is a two-momentum OPERATOR carrying BOTH legs** `q_out(ω,k)`, `q_out(ω,k′)` — ⛔ not a single-`k`
   multiplier, ⛔ not a one-leg left-quantized `a(x,k)=W_bg(x)σ(q(k))` (spec §3a). Both diagonalize the operator
   and delete the mode mixing (the rule-17 freeze). **FORM ABLATION (mandatory):** collapse the kernel's two legs
   (`k′→k`, i.e. `Ŵ_bg(0)`, or `q_out(k′)→q_out(k)`) and re-run — the kernel MUST move. A byte-identical
   `DTN_KERNEL`/`DTN_OPERATOR` under that collapse means the second leg was never live.
2. **The rigid-shift residual** (`k=k′`, `Ŵ_bg(0)`) must cancel on-shell (spec §3a) — check
   `DTN_RIGID_SHIFT_RESIDUAL` is a computed object, and that a WRONG kernel does NOT cancel (corrupt the kernel,
   confirm the residual moves).
3. **The permeable response is an operator inverse** `[I+(Λ_A/ρ_m²)Z]⁻¹` (spec §3b), ⛔ not a scalar division.
4. **The three dissipation objects are DISTINct** (spec §3b): (i) `H_a[Z]` bulk-radiation Hermitian part; (ii) the
   two-port permeable-port Hermitian form; (iii) the INDEPENDENT energy balance whose **face operand is the
   true-area traction pairing** (⛔ NOT `½Re(δp·V*)`, which equals the bulk flux and never sees `t_s`) and whose
   **bulk operand is the far-field Poynting flux at `|w|→∞`** (⛔ NOT `δp` at the face). **One-sided corruption:**
   flip the `t_s` sign in the face operand only — the energy residual must move (and the bulk operand must NOT),
   proving the two routes are independent. If flipping `t_s` moves both, the route is decoration.
5. **Blind re-derivation (the only cross-engine control):** the engine imports nothing — structurally scan for
   `Get`/`Import`/`<<`/`ReadString`/`OpenRead` and any absolute repo path (there must be none), and run the
   copy-to-empty-dir check. The T-a..T-i substrate must be **re-derived from the S11c-a maps + level set + the
   orientation law `s(n̂·ŵ)>0`**, ⛔ not transcribed. **FORM ABLATION:** the §5b `W_bg`-tilt ablation must MOVE
   the computed T-a and its DtN descendants — a byte-identical T-a under the tilt ablation is a transcribed
   constant, not a re-derivation.
6. **Every §5 control BITES.** For each of §5a rep-invariance (Eulerian vs Hanzawa/layer-potential — ⛔ not the
   secular global scaling), §5a one-sided tilt independence, §5b per-direction `W_bg` form ablation, §5c
   uniform-limit, §5d zero-jet (operand B = the BARE half-space `Z`, ⛔ no `W_0→W̄₀(1+η)`, ⛔ no two-face
   re-solve), §5e branch/momentum liveness (both one-leg sign-flips AND both one-leg freezes): confirm the control
   is a genuine second route, not an `A−A≡0` re-run. **Ablate each on a /tmp copy and report the literal residual
   before and after** — a control whose baseline and corrupted operands are byte-identical calls is not a control.
7. **`μ_θ` carried opaque** (spec §0/§3b): ⛔ not expanded into slab DOFs, ⛔ not reconstructed from S11c-b, ⛔ not
   copying the SymPy engine's per-case locals. It is one opaque operand `muTheta`.
8. **Reserved-name spellings** (for the downstream T7 join): the shared reserved symbols must use the
   `mechanical_lower_camel` spelling of their reserved snake_case names — `W_bg`→`WBg`, `w1_profile`→`w1Profile`,
   `eta_bg`→`etaBg`, `epsilon_shape`→`epsilonShape`, `sigma_W`→`sigmaW`, `Lambda_A_0`→`LambdaA0` (capital L),
   `tau_A`→`tauA`, `q_out`→`qOut`, `rho_m`→`rhoM`, `c_s0`→`cS0`, `mu_theta`→`muTheta`, the two density reps
   `rho_br`→`rhoBr` and `rho_br_bg_rho4_constant`→`rhoBrBgRho4Constant` (kept distinct). A wrong spelling
   (`w1` for `w1Profile`, one merged density) would be a trivial rename residual at T7 — report any.
9. **Emit-list completeness** vs spec §4 (every §3a/§3b primary + every §5/§6 control operand+residual, one tag
   per named object, cases as a keyed Association).
10. **The three script clauses** (spec §6): prints computed objects (⛔ never a prose/typed-answer conclusion, ⛔
    never a hand-typed CAS object with no data dependence — corollary 1); emits operand A, B, `A−B` before any
    guard (⛔ no tautological residual `A≔B/C` then `A·C−B`); no `VERDICT`/`PASS`/`FAIL`; a boolean is a typed CAS
    object (unevaluated relational / `STATUS_TOKEN`), ⛔ never a native `True`/`False`. ⚠ Report any `assert`/`Quit`
    that precedes the value it guards. Ask of every emitted physics claim: **WHICH LINE COMPUTED THIS?** — give the
    line or report it as uncomputed. Cross-check with `reduction/derived_or_declared.py` (triage, not a verdict).

## Required method (this is a SCRIPT — derive, then ablate)
Write your OWN derivation script (a small SymPy or Wolfram check of the load-bearing objects — the flat symbol
`Z_0=ρ_m ω/q_out`, the two-momentum kernel structure, the operator inverse, the far-field intensity) **before**
opening the artifact, and save **both the script and its literal stdout to named absolute paths**; without them
your derivation claims are discarded. Then ablate every load-bearing check on a **/tmp copy** and report the
literal diff. ⭐ A FORM ablation (change the STRUCTURE of a load-bearing object — collapse the two legs, flip a
sign and an off-diagonal) is MANDATORY and is the only thing that has ever caught the worst defect; a COEFFICIENT
rescale tests only arithmetic.

⛔ **Mathematica hygiene (this artifact spawns kernels):**
```
⛔ Wrap EVERY kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on.
⛔ NEVER raise the timeout, and ⛔ never run more than one kernel at a time (the licence has TWO seats).
⛔ Copy the artifact to /tmp and ablate the COPY. ⛔ Never modify the working tree.
⭐ Save every ablation script AND its literal stdout to named absolute paths, and report those paths.
⚠ If a background job is killed with a healthy-looking log, check `free -h` BEFORE anything else (an orphaned
   kernel leaks memory and the OOM killer takes unrelated jobs); `ps -eo pid,rss,pcpu,etime,comm --sort=-rss | head`.
```

## Physics filter
Report a finding only if it catches a way the physics could be wrong or a control that cannot fail; do not report
"the script would be wrong on a different input". State severity (MUST-FIX vs NIT) and, for each, the literal
ablation output that establishes it. If the engine is sound, say so and name the two or three objects you ablated
most closely with their literal before/after.
