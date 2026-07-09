# II-G2b (ledger_stage014) source map — pathA_31 truncation consistency (generalized eig + β_L0 sweep + N-convergence)

> Running-start prep captured 2026-07-08 (post stage013, this session) so the reshape directive can be authored without
> re-discovery. **All line refs below VERIFIED against the current sources** (recon 2026-07-08: `_sympy.py` 1359 lines,
> `.wl` 428, report 112, directive 289 — unchanged from the stage013 map's recon).
> Companion: `part2_gravity_atomic_split.md` (rows 013/014/015 + the Cluster-A reshape-cost map + the pathA_31 trip-ups);
> `stage013_pathA31_breathing_source_map.md` (the SIBLING map — 013 built the M/K this stage consumes).
> Build-order id **014**, Part II. **pathA_31 SPLITS 3-way: 013 (M/K projection, DONE+committed) + 014 (THIS: truncation
> consistency) + 015 (legacy-Hessian + Hellmann–Feynman force, NEXT after this).**
> **Source top-line: `BREATHING_CALIBRATED`** — the joint 3-stage verdict; **014 carries the truncation-consistency
> component (2/3)**: does the 2-mode `{α_a, α_L}` collective truncation faithfully capture the low spectrum of the full
> combined-basis generalized eigenproblem, across the β_L0 sweep and as N→16?

## ⭐ The THREE headline differences from 013 (READ FIRST)
1. **⚠ 014 is NUMERIC / float-BEARING — 013 was symbolic/float-free.** This is the single biggest design shift. 014 runs
   `scipy.linalg.eigh` (generalized eig `K v = ω² M v`), `scipy.integrate.quad` overlaps, and a β_L0 sweep of computed
   `o_1`/`o_2`/`ω²` floats. The anti-gaming teeth must be designed **around numeric outputs** (the overlap-floor predicate
   and the sweep window must be genuinely COMPUTED from the eig, and each per-tooth ablation must perturb a real float and
   confirm the fail fires). The dual-engine agreement is therefore a **transcript-level structural match** (both engines
   independently reach the same pass-pattern over the sweep + the same window boundary + the selected `o_k` to a loose
   numeric tol), NOT a symbolic `FullSimplify[...==0]`. Keep numpy/scipy in the reshape.
2. **⭐ The pathA_31 v1 REJECTION was a GAMED truncation threshold — THIS is the stage that carries that scar.** v1 typed
   the overlap floor `o_k≥0.9` and hand-picked the β-sweep to pass. The directive MUST require: the floor is a genuine
   predicate applied to the COMPUTED `o_1,o_2`; the β_L0 sweep window is COMPUTED (pass at low β, **fail at high β** — see
   §2); a mutant that hardcodes `pass=True` everywhere, or that never reads the computed overlap, must be caught by the
   per-tooth ablation. **This is the reshape's central able-to-fail obligation.**
3. **⭐ The overlap does NOT guard profile-correctness — carry the caveat, it is COMPUTED in the data (see §2, §5).** The
   `constant_one` counterfactual (α_a=1, a WRONG profile) **passes the overlap floor** (`o_1=1.0, o_2=0.974 ≥ 0.9`). Only
   the `degenerate_zero` profile (α_a=0) **fails** it (`o_2=0.223`, rank-deficient basis). So 014's overlap certifies
   *truncation-consistency* (the 2-mode span captures the low spectrum), NOT that the profiles are the correct harmonic
   lifts — that is 013's `𝓛₀[α]=0` residual + **015's HF mismatch** (the real profile guard). 014 must NOT claim the
   constant profile is *rejected* (that is 015); it carries the constant-passes-overlap as the honest scope limit.

## The `.wl` is a HYBRID (same family as 013's, but the 014 block is NUMERIC)
The pathA_31 `.wl` computes its own numeric Galerkin (`galerkinRow` = native `NIntegrate` + generalized `Eigensystem`,
L229–289) — a genuine independent engine — **but then diffs its results against the imported SymPy floats** (`sympyBetaSweep`
etc. via `Get[sympyExprFile]`; the `maxSweepDelta = Max[Abs[mmaSweep − sympyBetaSweep]]` cross-check L296–305). So the
014 reshape is the 013 pattern applied to the numeric block: **sever the sympy-float import; the `.wl` must apply its OWN
overlap floor to its OWN `galerkinRow` output and emit its OWN pass-pattern / window / selected `o_k`** — dual-engine
agreement then lives at the arbiter transcript level (both engines reach the same window + pass-pattern), not in a
`.wl`-imports-`.py`-floats diff.

## §1 The 014 slice (`.py` line ranges) — the CLEAN CUTS (all VERIFIED)
- **Module constants (014-owned; L50–55):** `EPS_TRUNC=0.1` (L50), `OVERLAP_FLOOR=1.0−EPS_TRUNC=0.9` (L51),
  `N_FINAL=16` (L52), `N_CONVERGENCE=[4,8,12,16]` (L53), `BETA_L0_SWEEP=[0.1,0.2,0.5,1.0,1.85,2.0,3.0,5.0,10.0,18.0,30.0,50.0]`
  (L54), `BETA_L0_FROM_R0=37.0/20.0` (L55). ⚠ `BETA_L0_FROM_R0` is the **frozen-packet anchor** consumed from 011/012/013
  (`β·L0=37/20`; §4) — predeclared, NOT best-fit over the sweep. The sweep grid + floor are the 014 anti-gaming surface.
- **`_baseline_functions(beta_l0, profile)` L447–473:** the numeric profile callables. `baseline` = `α_a=sinh(b(1−x))/sinh b`,
  `α_L=sinh(bx)/sinh b` (x=w/L0 normalized; L453–463). The two counterfactual profiles live here: `degenerate_zero` →
  `α_a≡0` (L468); `constant_one` → `α_a≡1` (L470). ⚠ These are 014's counterfactual generators — the guard (below) reads
  their overlaps.
- **`galerkin_overlap(beta_l0, n_modes, profile)` L476–557 — the 014 CRUX:** builds the `(2+N)×(2+N)` mass/stiff matrices
  (`{α_a, α_L, g_1..g_N}`, the `g_n=sin((n−½)πx)` lane modes L494–511) by `quad`; drops rank-deficient rows (`active`,
  L513–516); solves the generalized eig `eigh(stiff_active, mass_active)` L517; sorts L518–520; projects each of the two
  lowest eigenvectors onto the `{α_a,α_L}` subspace via the mass-Gram (`selector`/`gram`/`gram_pinv` L522–527) → the
  overlaps `o_1,o_2` L529–538; `min_omega12_squared` L540; `gap` L541; **`pass_row = o_1≥FLOOR and o_2≥FLOOR and
  min_ω²>0`** L542 (the genuine computed predicate — the v1 anti-gaming fix). Returns the per-(β,N) row L543–557
  (incl. `rank_deficient_basis` L555, `mass_condition` L556).
- **`build_truncation(profile)` L560–607 — the 014 assembly:** the sweep `[galerkin_overlap(β,N_FINAL) for β in
  BETA_L0_SWEEP]` L561; the `selected` row at `BETA_L0_FROM_R0` L562; the N-convergence `[galerkin_overlap(37/20,n) for n
  in N_CONVERGENCE]` L563; the **computed window** `beta_window` from the passing β's L565–572 (`min/max` of the passing
  set + the `predicate` string); `status=selected["pass"]` L574; the `beta_from_R0` provenance dict L588–595
  (`"R0 geometry alone does not derive Tw or K_eta"` — the calibration-honesty line, ties to §6); the `g_lane_eigenproblem`
  descriptor L599–605; the `combined_basis` string L606.
- **The `build_truncation` calls + numeric-WL append (inside `symbolic_engine`) L760–809:** `truncation=build_truncation("baseline")`
  L760; `wrong_zero_trunc=build_truncation("degenerate_zero")` L761; `wrong_const_trunc=build_truncation("constant_one")`
  L762; **the `numeric_wl` scratch-export block L763–809** (writes `sympyBetaSweep`/`sympySelectedOverlap`/
  `sympyWrongZeroOverlap`/`sympyWrongConstOverlap` to `SYM_EXPR_WL` via `open("a")` L808–809 — **THIS IS THE 014 BRIDGE
  TO SEVER**, §3).
- **The `wrong_o_k` slices of the counterfactual guard (014-owned) L819–823 + L833 + L837:** in the `degenerate` block,
  `wrong_o_k={o_1, o_2, rank_deficient_basis}` L819–823 (the overlaps of the α_a=0 profile — `o_2=0.223` FAILS floor);
  in the `nontrivial` block, `wrong_o_k={o_1, o_2}` L833 + `overlap_passes=bool(wrong_const_trunc["status"])` L837 (the
  α_a=1 profile — `o_1=1.0, o_2=0.974` PASSES floor → the caveat, headline diff #3). ⚠ **The rest of the guard is NOT
  014:** `M_det`/`M_posdef` (L817–818, L831–832) = 013; `F_a_dist`/`F_a_legacy_frozen`/`hf_mismatch`/`fails` (L824–827,
  L834–838) = 015. 014 reads ONLY the `wrong_o_k`/`overlap_passes` fields.
- **The 014 verdict rung L1050–1052:** `trunc=payload["truncation"]`; `if not trunc["status"] or
  trunc["min_omega12_squared"]<=0: return "BREATHING_FAIL_TRUNCATION_INCONSISTENT"`. This is 014's own fail token in the
  joint chain (013→structure/dynamical, 014→truncation, 015→HF).
- **⭐ CLEAN CUT — 014 owns the numeric-Galerkin machinery + its constants + the overlap slices of the guard ONLY. Do NOT
  recompute (013/015 territory):**
  - 013: `symbolic_engine` harmonic profiles + `M_AB`/`K_AB` ∫dw projection L630–679; the dimensional block L326–444;
    `M_det`/`M_posdef` guard fields. **013 is DONE + committed — CITE its `M_AB`/`K_AB`, do NOT re-derive.**
  - 015: `build_structure_gate` L217–311; `E_geom`/`H_legacy` + HF force `F_dist`/`F_legacy` + `hf_*_ok` L681–703; the
    `hf_mismatch`/`F_a_*` guard fields L824–838; the static-dynamic limit.

## §1b The `.wl` 014 slice (VERIFIED) — L209–306
- `profileFunctions[b, profile]` L209–227 (numeric profiles; the two counterfactual profiles);
  `galerkinRow[b, nModes, profile]` L229–289 (native `NIntegrate` mass/stiff L240–253 + generalized `Eigensystem[{ka,ma}]`
  L268 + subspace overlaps L277–288 + `min12`/`gap`; returns `{b, o_1, o_2, min12, gap}` L288). This is the genuine
  independent numeric engine — KEEP it.
- The calls L291–294: `mmaSweep` L291 (⚠ currently reads the β grid from `sympyBetaSweep[[All,1]]` — after severing, the
  `.wl` must declare its OWN sweep grid); `mmaSelected` L292; `mmaWrongZero` L293; `mmaWrongConst` L294.
- **⚠ THE BRIDGE TO SEVER L296–305:** `maxSweepDelta = Max[Abs[Flatten[mmaSweep − sympyBetaSweep]]]` L296 +
  `selectedDelta`/`wrongZeroDelta`/`wrongConstDelta` L297–299 + the `numericChecks` assoc L302–305 (all four are
  `.wl`-diffs-`.py`-floats). SEVER: replace with the `.wl` applying its OWN `OVERLAP_FLOOR` to its OWN `galerkinRow`
  output → its OWN pass-pattern + window + selected `o_k` + the two counterfactual overlaps. (013 `.wl` block = L29–207 —
  the symbolic M/K, already reshaped/committed; EXCLUDE. 015 `.wl` = the HF/structure block — EXCLUDE.)

## §2 The 014 claim-set (compute + assert; report quotes)
The physical numbers (report `## Truncation Consistency` :37–59; the honest anchors a fresh fidelity leg re-derives):
- **Selected point (β_L0=1.85=37/20, N=16):** `o_1=0.99311, o_2=0.98776, min(ω_1²,ω_2²)=3.42252, gap=2.22787` (report
  :40). `pass=True` — the 2-mode truncation is clean at the physical wall packet.
- **The COMPUTED β_L0 sweep window (report :44–56):** `pass=True` for β_L0 ∈ {0.1, 0.2, 0.5, 1.0, 1.85, 2.0, 3.0};
  **`pass=False` for β_L0 ∈ {5.0, 10.0, 18.0, 30.0, 50.0}** (at β_L0=5, `o_1=0.8598 < 0.9`). `beta_window` computed =
  `{min_in_sweep: 0.1, max_in_sweep: 3.0}` (report :42). ⭐ **This is the honest validity caveat: clean 2-mode truncation
  only for order-unity wall stiffness** — the STATUS/source phrasing is `K_η/T_w ≲ 2.6` (β²=K_η/T_w; the sweep boundary
  sits between β_L0=3 pass and β_L0=5 fail). ⚠ Treat `≲2.6` as an EMERGENT number from the computed sweep, a fidelity
  anchor — do NOT let the directive type it as a pass criterion. The falsification content = the window has a **real upper
  edge** (sharp walls genuinely fail), not an everywhere-pass.
- **N-convergence (report :58):** at β_L0=1.85, `o_1` is stable across N=4/8/12/16 (0.99301→0.99310→0.99311→0.99311);
  `ω_1²`→3.4225, `pass=True` at every N. The truncation certificate is N-converged (not an artifact of a lucky small N).
  ⚠ `mass_condition` grows with N (3.1e3→1.9e5) — a benign conditioning growth, not a failure; note it so a reviewer does
  not mistake it for one.
- **The two counterfactual overlaps (report :79–82; §5):** `degenerate_zero` → `o_2=0.2227 < FLOOR`, `rank_deficient_basis
  =True` → FAILS truncation (a genuine 014 tooth). `constant_one` → `o_1=1.0, o_2=0.9738 ≥ FLOOR`, `overlap_passes=True`
  → PASSES the overlap (the caveat: overlap ≠ profile-correctness; rejection is 015's HF).
- **What 014 EARNS (structure) vs CALIBRATES:** EARNED = the truncation-consistency certificate is genuinely computed (real
  generalized eig, real floor predicate, real window, N-converged, degenerate genuinely fails). CALIBRATED = the wall
  packet `{T_w, K_η, β}` values feeding the profiles are 013's already-counted calibration knobs (§6) — 014 adds none.
- **The 014-scoped verdict + joint composition:** compute the 014 component (the truncation-consistency certificate) and
  print the joint `BREATHING_CALIBRATED = (013: M/K operator-projection + (a,L) closure) ∧ (014: truncation consistency /
  validity window) ∧ (015: legacy-structure + HF force)`. Do NOT type `BREATHING_CALIBRATED` as 014-earned alone
  (the 009/010 + 011/012 + 013 pattern).

## §3 Reshape cost (the bridge to sever) — hybrid `.wl`, NUMERIC block
Same scratch-yaml payload-mirror family as 011/012/013, applied to the numeric Galerkin block. Strip:
- **`.py` → the numeric-WL scratch append:** `SYM_EXPR_WL` def L46; the `numeric_wl` list build + `SYM_EXPR_WL.open("a")`
  append L763–809. (The symbolic `SYM_EXPR_WL.write_text(...)` L755 is 013's bridge — already handled in the 013 reshape;
  confirm the standalone 014 script does zero file I/O.)
- **`.wl` READS + diffs (SEVER — make the `.wl` self-certify):** `Get[sympyExprFile]` L20 (the numeric symbols
  `sympyBetaSweep`/`sympySelectedOverlap`/`sympyWrongZeroOverlap`/`sympyWrongConstOverlap` come from here); the sweep grid
  read `sympyBetaSweep[[All,1]]` L291; the four `*Delta` diffs L296–299; the `numericChecks` assoc L302–305. Re-target: the
  `.wl` declares its OWN sweep grid + floor + N-grid (predeclared constants, same values as `.py` — NOT readback), applies
  the floor to its OWN `galerkinRow` output, and asserts its OWN pass-pattern/window/selected-`o_k`/counterfactual-overlaps.
- **⭐ The `.wl` 014 KEEP (make it fully independent):** `profileFunctions` L209–227 + `galerkinRow` L229–289 (native
  `NIntegrate`+`Eigensystem`) — re-target its asserts from `sympy*`-diffing to its OWN floor predicate. **Zero file I/O.**
- **Dual-engine agreement (numeric, transcript-level):** both engines independently emit the sweep pass-pattern (12 booleans),
  the window `[0.1, 3.0]`, the selected `o_1,o_2` (agree to a loose tol, e.g. `<1e-6`), the N-convergence stability, and the
  two counterfactual overlaps. The ARBITER compares the two transcripts (the 013 model) — the scripts do NOT cross-read.
- **Arity discipline (standing, stage007–012 lesson):** def/call arity scan + unevaluated-leakage transcript scan on the
  reshaped `.wl` (a mismatched call silently skips its section at exit 0).

## §4 Consumed / exported
- **Consumes (cite, dual-site integrity, don't re-derive — the 013→014 hand-off):**
  1. **013's `M_AB`/`K_AB`** — the combined-basis generalized eigenproblem `K v = ω² M v` is over the `{α_a, α_L, g_1..g_N}`
     basis whose 2×2 `{α_a,α_L}` block IS 013's operator-projected `M_AB`/`K_AB` (report :39, `combined_basis`). ⚠ 014
     must CITE 013's M/K (dual-site integrity), NOT re-derive the profiles / M / K / EOM LHS (that is 013, done+committed).
     Design the dual-site: e.g. the `{α_a,α_L}` sub-block of 014's numeric mass/stiff Gram must MATCH 013's symbolic
     `M_AB`/`K_AB` evaluated at the frozen packet (a genuine cross-check that 014's numeric basis is the same object 013
     projected) — that is the honest 013→014 seam, and it fires if the profiles drift.
  2. **013's frozen packet `{L0=37/20, T_w=1, K_η=1, β=1}`** → `BETA_L0_FROM_R0=37/20` (L55), the sweep anchor. `β·L0=37/20`
     is the branch-determinable `L/a=37/20` (`[[project-ansatz-value-catalog]]` route (b)) — cite, don't refit.
  3. ⚠ **`c_S` is NOT consumed** — the matter-sector `c_s`/BdG `k⁴` stays DEFERRED (`kξ≪1`, `phonon_limit_caveat` report
     :101). 014 is a static-spectrum truncation test; no speed enters.
- **Exports:** the validity window + the truncation certificate → **015** (015's HF/structure runs on the same 2-mode
  closure 014 certifies) + the ℓ=0 map (022/023). 014's share = the generalized-eig truncation consistency, the β_L0
  window, the N-convergence, and the honest overlap≠profile-correctness caveat.

## §5 Teeth candidates (014-specific — ALL float-bearing; per-tooth ablation MANDATORY)
⭐ Because 014 is numeric, EVERY tooth must be shown to fire by perturbing a REAL computed float (not a symbolic
`==0`), per [[feedback-per-tooth-ablation]]. Keep/assign to 014:
1. **The overlap-floor predicate is genuinely computed** (`pass_row = o_1≥0.9 and o_2≥0.9 and min_ω²>0`, L542) — ablation:
   feed the `degenerate_zero` profile → `o_2=0.223 < 0.9` → `pass=False` (the floor genuinely reads the computed overlap;
   a mutant that hardwires `pass=True` or ignores `o_2` must be caught).
2. **The β_L0 sweep window has a real upper edge** — ablation: a high-β point (β_L0=5) genuinely computes `o_1=0.86 < 0.9`
   → `pass=False`; the window `[0.1, 3.0]` is COMPUTED from the passing set (L565–572), not typed. A mutant that pins every
   sweep row to pass must flip the window / be caught.
3. **The `degenerate_zero` counterfactual fails truncation** (rank-deficient basis + `o_2<floor`, L819–823) — the 014
   able-to-fail companion. Ablation: the null α_a genuinely collapses the `{α_a,α_L}` span → the overlap drops.
4. **The `constant_one` counterfactual PASSES the overlap** (`overlap_passes=True`, L837) — ⚠ this is NOT a fail tooth; it
   is the **honest-caveat tooth**: 014 must record that a wrong profile passes the overlap (so the overlap does not certify
   profile-correctness). Per-tooth ablation here confirms the guard genuinely reads `wrong_const_trunc["status"]` (a mutant
   that fakes `overlap_passes=False` would falsely claim the overlap DOES guard profiles — must be caught). The rejection
   of the constant profile is 015's HF mismatch — do NOT fold that into 014.
5. **N-convergence is genuine** — ablation: the o_k are stable across N=4→16 (report :58); a mutant that returns a single
   N or a non-converging series must be visible.
6. **The 013→014 M/K dual-site seam** (§4.1) — the `{α_a,α_L}` sub-block of 014's numeric Gram matches 013's symbolic
   `M_AB`/`K_AB` at the frozen packet; a profile drift fires it on both engines.
7. **`.wl` def/call arity scan + unevaluated-leakage scan** (standing).
⚠ **NOT 014 (do not pull in):** the harmonic-lift residual `𝓛₀[α]=0` (013); the `forbidden_fit_flags`/`M_det→0` structure
guard (013); the HF `x−x` two-route enforcement + `hf_mismatch` (015); the corrupt-`[T_w]` dimensional probe (013).

## §6 Register expectation — likely ZERO new counted knobs (CONFIRM at registration)
014 is the **truncation-consistency of 013's already-counted closure** — it introduces no new physics constants. Expected
register outcome: **zero new counted knobs.** The knobs the sweep uses (`T_w`, `K_η=T_w β²`, `β`, `L0`) are all 013's
(`{μ_η, T_w, β}` counted in 013 at count 13; `K_η` = R29 manifestation; `L0`/`β·L0=37/20` = ACTION-geometry). The 014
numeric controls — `EPS_TRUNC=0.1` / `OVERLAP_FLOOR=0.9` (the truncation tolerance), `N_FINAL=16` / `N_CONVERGENCE`
(basis size), the `BETA_L0_SWEEP` grid — are **method/tolerance parameters, tracked-not-counted** (like a solver rtol; a
convergence knob, not a physics DOF). A structural edge may record the validity-window caveat (`clean 2-mode truncation
only for K_η/T_w ≲ order-unity`) if useful — it discharges nothing. ⚠ CONFIRM at registration + Codex-verify per
[[feedback-parameter-register-every-stage]]: verify none of the sweep/floor/N controls is mis-classed as a derived physics
knob, and that the frozen packet is a CITE of 013's rows (no double-count). This mirrors 010/011/012's "zero new counted
knobs" registrations.

## Verdict tokens + honest scope
014 carries the **truncation-consistency component (2/3) of `BREATHING_CALIBRATED`**: the combined-basis generalized eig
+ the COMPUTED β_L0 window + the N-convergence certify the 2-mode `{α_a,α_L}` truncation is clean at the physical wall
packet (EARNED-structure); the underlying `{T_w,K_η,β}` are 013's CALIBRATION inputs (→ `CALIBRATED`, not `PASS`). Caveats
to CARRY honestly: **(a)** clean 2-mode truncation ONLY for order-unity wall stiffness `K_η/T_w ≲ 2.6` — sharp walls
genuinely fail (a real validity edge, computed from the sweep, not a rescue); **(b)** the modal overlap does NOT guard
profile-correctness (`constant_one` passes it) — the profile guard is 013's residual + 015's HF; **(c)** the phonon-limit /
BdG `k⁴` stays deferred (`kξ≪1`). CITED: 013's `M_AB`/`K_AB` + the frozen packet `L0=37/20` (branch-determinable). Nothing
here re-derives the profiles/M/K (013) or the HF force (015).

## Process (unchanged, calibrated — the UPDATED per-stage pipeline)
Author the II-G2b reshape directive (§1 the clean 014 slice + KEEP-014-only numeric block + §2 faithful cover + §3
bridge-strip incl. the hybrid-`.wl` numeric re-independence + §5 the float-bearing 014 teeth with per-tooth ablation + §6
the zero-new-knobs register note) → **Codex xhigh design-review** → fold to `DIRECTIVE_CLEAN` → **⭐ final Grok-4.5 headless
compute-verify pass** on the same prompt (assess + independently verify each catch — Grok is noisy; it caught a
kernel-preserving residual-tooth defect on stage013) → fold → **Codex confirm-pass on the Grok-folded directive** →
**pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`, background, stdin-from-file,
absolute paths, xhigh) → dual-engine both exit 0 (repo root + foreign CWD) → arbiter re-run via the runners → full
tri-review (fidelity re-derives the o_k/window/N-convergence + adversarial-with-**mandatory per-tooth ablation** over
EVERY check incl. the numeric floor/window/counterfactual/dual-site/arity) → remediate → fresh-agent re-verify → bump
counts 13→14 → parameter register (confirm **zero new counted knobs**) + Codex-verify → note/card/`\input{stages/stage_014}`
+ registration → rebuild PDF → commit + docs/memory sync. Orchestrator authors notes/cards/LaTeX/registration; Codex codes.
Target stem: `ledger_stage014_breathing_truncation_consistency` (confirm slug at authoring).
