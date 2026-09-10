# Independent physics review — S11c-d SymPy BUILD DIRECTIVE (decision review, the G2 TRIGGER)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_sympy_build_directive.md`

## What to check
This is an **orchestrator-written SymPy BUILD DIRECTIVE** for step **S11c-d** (profile-conditioned transverse↔thickness
scattering / mixing / leakage). It is the **decision review that gates the builder** (`CLAUDE.md` G2 TRIGGER: no
builder launches until its build directive has had two legs). It governs the **SymPy engine only** (the blind Wolfram
engine + the T7 comparator are separate downstream artifacts).

The directive is a **THIN** directive: the physics authority is the **CLEARED** spec
`directives/S11c_d_SHARED_PHYSICS.md` (v10, committed `399a8516`, round-10 dual-engine gate both legs SOUND). The
directive must **POINT at** the spec, ⛔ never restate/re-derive/drift the physics. Its job is to fix only the
**build-mechanical** layer + the deferred §1c Fourier-reduction **element census against the real rows** + leak
discipline. It carries **physics-bearing content** (the §1 HELD-PHYSICS pointers, the §3 census, the debt rendering),
so beyond the one-pass decision checks it also meets the spec-review bar: **any finding that would change what is
computed or what may be claimed is a must-fix**, not a nit.

Your job: adjudicate the directive's build-mechanical **decisions** and its physics-bearing **pointers/census** by
**evidence**. This is a decision review (one two-leg pass), but a physics-bearing defect routes to the spec bar.

## What you are handed (there is no do-not-read list — what a leg must not use, it is not given)
- The artifact above.
- The physics authority: `directives/S11c_d_SHARED_PHYSICS.md` (the cleared spec — this is the source of truth; the
  directive must faithfully point at it and must not contradict, drop, or drift it).
- The **real export files** the directive pins facts against: `scripts/S11c_b_exports.py` (base),
  `scripts/S11c_c1_exports.py` (c1 delta), `scripts/S11c_c2_exports.py` (c2 delta), and the loader
  `scripts/ledger_fold.py`. Python has sympy; work from `/var/projects/toy_physics/research/pde_ledger_v3`.
- The c1/c2 SymPy build-directive precedents: `directives/S11c_c2_sympy_build_directive.md`,
  `directives/S11c_c1_sympy_build_directive.md`.
- The build + review skills: `.claude/skills/build/SKILL.md`, `.claude/skills/review-legs/SKILL.md`; and `CLAUDE.md`.
- The orchestrator's own grounding for the directive's §2/§3 facts:
  `directives/_measurements/S11c_d_sympy_build_directive_census.md`. ⚠ This is the orchestrator's fact-lookup, NOT a
  source of truth — a §2/§3 factual claim is confirmed **only** by your own re-run against the real rows, ⛔ never by
  reading this file and agreeing.

## Required method
**This is a DOCUMENT (a directive), and it PINS FACTS against the real ledger rows — so use BOTH:**

1. **Document branch.** Read the **source of truth first** — the cleared spec `S11c_d_SHARED_PHYSICS.md` (and the
   relevant c2 exports/step-record pointers it cites) — form your own view of what the spec establishes, what it
   defers to the build directive, and what it withholds, and **only then** read the directive. Quote both sides
   (spec vs directive) for every finding. Check the directive **does not restate/drift the physics** and **does not
   contradict, drop, or silently narrow** any cleared-spec control.

2. **Computational verification of the pinned facts (⛔ a prose claim is discarded — run it, save the script + its
   literal stdout to a named absolute path, and report those paths).** Settle these contested questions **by
   computation**, from the real files:

   - **(Q1 · import wiring §2)** Run `load_model("scripts/S11c_b_exports.py", "scripts/S11c_c1_exports.py",
     "scripts/S11c_c2_exports.py")`. Verify: base **2441** + c1 **44** + c2 **70** → fold **2555**, strictly
     additive, `overwrites == []`, all pairwise key-intersections empty. Verify `check_consumer` with the two closed
     rows as roots resolves (closure ≈ 213, no ambiguity). Is the **IMPORT_KEYS rule** correct (= the build's actual
     `fold[key]` direct-lookup set; `assert_lookups_equal_manifest` fails on an undeclared lookup AND a declared-but-
     unused key)? Is the **provenance hazard real** — do BOTH the open S11c-b `slab_operator`/`coupling_kernel` AND
     the c2 `s11cc2ClosedSlabOperator`/`s11cc2ClosedCouplingKernel` exist in the fold, and does the guard pass on
     either (so binding the open rows is a SILENT wrong-physics swap)? Is the exact carrier casing right
     (`s11cc2Fieldtheta` lowercase)? Do the two `Closed*` rows lack a `*Dimension` companion?

   - **(Q2 · the Fourier census §3 — the LOAD-BEARING completeness question)** Introspect
     `s11cc2ClosedSlabOperator` and `s11cc2ClosedCouplingKernel` (per `(α,ρ)` case). Verify the directive's census is
     **correct AND COMPLETE**: (a) all four c2 hats (`s11cc2FourierW1ProfileHatTransfer`,
     `s11cc2FourierW1ProfileJetHat{1,2,3}`) present as **applied functions** at the **transfer** argument
     `(k_out−k_in)` AND the **middle-leg** arguments `(k_out−k_mid)`, `(k_mid−k_in)` with `k_mid=s11cc2MiddleMomentum
     {1,2,3}`; (b) the c1-level snake_case hats ABSENT; (c) explicit `sp.Integral` measures over `d³y`
     (`s11cc2Y{1,2,3}`) AND `d³k` (`s11cc1_k_output/input_{1,2,3}`, `s11cc2MiddleMomentum{1,2,3}`); (d) **zero**
     `DiracDelta` in either closed row, and 3 `DiracDelta` on `dtn_kernel` at `(k_output−k_input)`. **Is any
     convention-bearing 3-D element PRESENT in the real rows that the census MISSES** (a missed element left in a 3-D
     convention is a must-fix)? Is `s11cc2OutgoingNormalMomentum` handled? ⛔ Is the census over-claiming an element
     that is not actually there?

   - **(Q3 · recipe-creep guard §3)** Does §3 pin **WHAT** to reduce (the census of elements present in the real
     rows) **WITHOUT** supplying **HOW** (⛔ no typed `[L_W/(2π)]`/`(2π)²L_W`/`δ²(Q_∥)` map; ⛔ `dtn_kernel` not bound
     as the reduction convention; the `(2π)`/dimensional arithmetic is each engine's own **computed** reduction)?
     Does it correctly require the **both-operand** reduction record (3-D operand + reduced operand) and the
     reconstruction round-trip, and correctly forbid the `A−A` tautology? Or does it **drift into specifying the
     reduction result** (which the spec deferred here precisely to avoid — the Fourier recipe forced 5 spec rounds)?

   - **(Q4 · leak discipline §6)** Is the S11c-d **expected-value acceptance criterion** properly WITHHELD (the
     falsification numeric bound / the `O(1)` grating reductio — orchestrator-side, `R1`-blocked)? Is the c2
     cross-engine operand **DEBT rendered QUALITATIVELY** — ⛔ no c2-import status counts, ⛔ no expected value/sign/
     order/parity/grade/baseline handed to the builder as an S11c-d target? Does any **prohibition name the step's
     real expected shape** (a prohibition leaks as surely as an assertion)? Grep the directive for co-occurring
     step symbols in proximity and read every hit.

   - **(Q5 · physics-bearing pointers §1, the spec bar)** Does §1's HELD-PHYSICS faithfully carry the cleared spec's
     load-bearing controls — the **reduced-representation rule** (§2 governing §§3–6: every downstream object built
     from the §1c-reduced rows, ⛔ never unreduced 3-D content); the **withdrawn-F ⇒ compute `K_0`/`A_0`, never type
     `=0`**; the **FULL** off-diagonal vertex (tilt+modulus-gradient+advection, ⛔ not a single subchannel); the
     **reduced-block-vs-kernel residual as a surfaced finding**; **no global `ω(k)`**; the **two DISTINCT** photon-
     kill channels with the profile-functional conditional bound pole (⛔ no 1-D weak-well theorem, `N13`); the
     conditional `N12` labels; the c2 **operand DEBT as a live premise** (name+propagate, ⛔ don't resurface/pre-
     adjudicate)? Does the directive **misrepresent, drop, or silently narrow** any of these (a spec-bar must-fix)?

   - **(Q6 · builder-lane + script obligations)** Are the three script clauses + the structural rule + the four
     corollaries present and correct (input-driven; tag-name-is-output; **no tautological residual**; the **§5a
     routes correctly NOT dressed as a kernel-N6 closure**; emission never conditional on value; FORM ablation via
     §5c)? Does the directive keep the builder in its lane (build → verify-own-deliverable → report → **stop**; ⛔ no
     launching legs/comparator/WL/downstream; ⛔ no reading /build or /review-legs)?

   - **(Q7 · EMIT vs EXPORT §5)** Is the export membership sound — grounded in S11c-e's DECLARED scope (the
     flux-normalized conversion observable + leakage + confinement whose weak limit reproduces S11c-d), with the
     no-S11c-e-manifest caveat? Is anything exported that S11c-e won't bind, or emit-only that it will need?

## Physics filter
Report a finding only if it catches a way the **physics or the build could go wrong** — an incorrect or **incomplete**
census, a wrong/silent import binding, a leaked expected value, a misrepresentation/drop of a cleared-spec control,
recipe-creep in §3, or a broken builder-lane/script-clause obligation. ⛔ Do not report "the directive would be wrong
on a different input." For each finding, state whether it changes **what is computed or what may be claimed** (a
must-fix routed to the spec bar) or is a one-pass decision nit, and quote spec-vs-directive (or command+stdout).

## Ablation sandbox
Copy any file you introspect to `/tmp` and work on the copy; ⛔ never modify the working tree. Save every verification
script AND its literal stdout to named absolute paths and report those paths.

## Bounds
Write your report and exit. ⛔ Do not spawn agents or build `run_all`/`watch`/supervisor orchestration — a leg twice
deadlocked doing exactly that; iterating to clearance is the orchestrator's job, not the leg's. Finish IN-TURN,
foreground blocking only.
