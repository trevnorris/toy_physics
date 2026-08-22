# Independent review — S11b WL engine AFTER repair (F-WL-1/2/3)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl`
(a Mathematica script just repaired in place; ~2272 lines). It re-derives the S11b interface coupling law
from the spec alone and emits one tag per named object; ⛔ no VERDICT.

## What to check
The repair (a) fixed `ZPERM_SLICE_MAP` (F-WL-1: extraction + static level), and (b) made genuine a layer of
previously-tautological/decoupled checks (F-WL-2 energy; F-WL-3a kernel-orientation, 3b causality, 3c
grazing). Verify each repaired check is now GENUINE — its two operands come from INDEPENDENT routes and it
can FAIL under a one-sided form change — and that the repair opened NO new tautology and did NOT change the
protected physics. ⛔ Do not diff against any expected value (the orchestrator holds the withheld references).

## What you are handed
- The artifact above and its repair directive
  `research/pde_ledger_v3/directives/S11b_wl_engine_repair_directive.md` (F-WL-1/2/3 + the acceptance).
- The physics authority `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` (§4, §6, §8, B2c).
- ⛔ You are NOT handed the sibling SymPy engine or its output; derive independently where needed.

## Required method — SCRIPT; ablation mandatory
Derive independently where a claim is at stake; write+run your own script, save script+stdout to named
absolute paths, report them. ⛔⛔ FORM-ablate each repaired check on a /tmp COPY and report the LITERAL diff:
- **F-WL-1**: is `S11B_ZPERM_SLICE_MAP` now extracted from the flux RHS (not the residual zero-form) and
  ω-independent (static)? Test representation-invariance: corrupt `equationZeroForm` / swap the closure
  LHS↔RHS — the map must NOT move (a residual-zero-form extractor WOULD). Is the emitted map free of ω, τ?
  ⛔ Do not state the value; report whether it is representation-invariant and static.
- **F-WL-2**: are the slab-power and bulk-power operands from INDEPENDENT routes (off-shell EOM vs supplied
  bulk acoustic power), NOT both from the solved `fullSystem`? One-sided corruption: break ONLY the slab
  route (a sign/form change to the EOM) — the residual must move; break an unrelated object — it must not.
  Is `PRESSURE_WORK_SIGN_CHECK` still ±the same X (still tautological)?
- **F-WL-3a/b/c**: 3a — does `KERNEL_ORIENTATION_IDENTITIES` now extract from the assembled pre-elimination
  closure (bites when a kernel is mis-oriented)? 3b — does `CAUSALITY_CHECK` now actually inspect records (a
  missing/renamed record FAILS)? 3c — does `GRAZING_MODE_CLASSIFICATION` classify from a non-degenerate
  leading-order block (not `q·MATRIX`=0)? For each, one-sided-corrupt its own load-bearing object and show
  the residual/classification moves; corrupt an unrelated object and show it does not.
- **No regression** (Scope): impedance/regimes, added mass, `S11B_ZPERM_SLICE`, transverse, breathing slice,
  longitudinal dispersion, roots must be behaviourally unchanged (except `KERNEL_PROPAGATION_RESIDUALS`,
  which F-WL-3a may legitimately move). Confirm.
Report any `assert`/`Abort` before an emit, and any check whose residual is still zero by construction.

## Ablation sandbox — MATHEMATICA, TWO-SEAT LICENCE (binds both legs identically)
⛔ Copy the .wl to /tmp and ablate the COPY. ⛔ Never modify the working tree.
⛔ Wrap EVERY kernel run in `timeout 600`; ⛔ NEVER raise it; ⛔ NEVER run more than one kernel at a time.
⛔ Kill on RSS > 6 GB; if a kernel is killed, run `free -h` FIRST, then hunt orphans with
`ps -eo pid,rss,pcpu,etime,comm --sort=-rss | head`. Save every ablation script + literal stdout to named
absolute paths and report them.

## Physics filter / report — under ~30 lines
Report a finding only if it is a way the physics could be wrong or a check that still cannot fail. Numbered
findings (tag, line, ablation+literal stdout path, why); end with a verdict: sound to commit, or repair again?
