# S11c-b Wolfram engine — build directive

## Authority and boundary

Write `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` in full. Its product is
the flushed stdout tag stream; that is its only write. No other `.wl`, `.py`, `.out` or exports file is a
build premise or an input.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (HEAD `7023420d`) is the sole
physics authority, together with the two sibling SHARED PHYSICS specs it inherits by reference
(`S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md` — these are shared physics authorities, ⛔ not engine
outputs). Every equation, inherited object, ansatz, premise, supplied law, derivation route, control,
dimension obligation and tag obligation enters from those specs alone, and the S11c-b spec wins every conflict
with this directive. Implement its §§0–8 directly. Point to rather than duplicate: §1 (the inherited setup,
DOF, local curl/div sector labels, stored energy and symmetry inputs, uniform decoupling context); §2 (the
background ansatz, which quantities vary, the activated admissibility computation); §3 (the variable-
coefficient basis and new invariants, the divergence-form operator, the off-diagonal kernel, the §3d
background-order balance); §4 (the objects to compute); §5 (the routes and controls); §6 (method, dimensions,
the three script clauses, no VERDICT); §7 (the parallel tag grammar); §8 (supplied vs computed, the builder
report). Add no expected value or acceptance criterion (`CLAUDE.md` rule 5).

⛔⛔ **This is the blind engine.** The build's inputs besides the specs carry no S11c-b computed physics: this
directive, `CLAUDE.md`, and the mechanical precedent named below. **The engine this build writes reads no
file at run time — it imports nothing and re-derives from the specs' equations alone**, including
re-deriving the S11c-a shape-derivative substrate objects the operator consumes (from the S11c-a spec's
§§1–4), ⛔ never importing `S11c_a_exports.py` or any transcript.
`S9_export_chain_rebuild_directive.md:16-18`: this independence is the design's **one** blindness control,
and nothing else may be built pretending to be one. ⛔ This engine imports no `LEDGER` and writes none. Per
spec §7, declare the inherited constants of §1/§2 (`c_s0`, `mu_R`, `rho_br`, `W_0`, `e_W`, `rho_m`, the drain
`v_bulk_normal_0`) and the reserved operand `mu_theta` as the engine's own symbols exactly as the specs
supply them, and emit **no** import tag (that absence is by construction, ⛔ not a denylist). ⛔ The varying
background profiles of §2/§7 (`W_bg`, `mu_R_bg`, `w1_profile`, `m1_profile`, `sigma_W`, the density
representatives, `e_W_bg`, …) are **not** inherited constants — they originate from §2a's ansätze, are
re-derived here, and ⛔ must not reuse a constant key (§7: `mu_R_bg` must not reuse `mu_R`, nor `W_bg` reuse
`W_0`). ⛔⛔ **There is no do-not-read list, and none may be added** (rule 12): blindness is enforced by the
engine importing nothing, ⛔ never by a sentence forbidding a read.

Do not add parallel machinery for the spec's contracts: no per-cell completion registry beyond §4/§5's task
structure, and no directive-local exit policy. §6's emission-not-conditional-on-value rule applies as a
property, without named exceptions. ⛔ **No `VERDICT`, `PASS` or `FAIL`** (§6): a physics disagreement exits
0, and a boolean-valued test is emitted **as the typed CAS object retaining its operands** (the unevaluated
relational, e.g. `a == b`), ⛔ never a host-language native or evaluated boolean (the frozen T7 comparator
rejects a native boolean as a residual operand).

## The conventions shared with the sibling engine

Per spec §7 the emitted grammar is `WL_S11CB_<QUANTITY>`: replace the leading `S11CB` of every base tag named
in §§2–6 with `WL_S11CB`, one tag per named object, a single object's branch/face/DOF/density/direction/
source cases carried as a keyed CAS map in that object's payload, ⛔ not as separately invented tag names. Any
unavoidable engine-local tag carries `_LOCAL_` immediately after `S11CB`, with one local-tag inventory.
⚠ Wolfram symbols cannot contain `_`; that constrains only the engine's **internal** symbol spellings —
⭐ the **emitted-name strings** must be exactly the spec's `S11CB_<QUANTITY>` bases under the `WL_` prefix,
because the downstream comparator joins the two engines on those strings. Take the tag list from **every base
name in §§2–6**; ⛔ do not narrow it to a subset, rename, abbreviate, or merge a base tag.

The four conventions the spec pins, anchored identically to the sibling engine:
- the sector labels are the §1a **local** `∇×`/`∇·` operators, ⛔ not a global projector (`N5`);
- the variable-coefficient energy basis is a **non-unique quotient — construct your OWN representative**, ⛔
  do not pin one, and ⛔ do not obtain the energy by `W_0→W_bg` substitution into the uniform `U` (it omits
  the new gradient-of-background invariants, spec §3a/`N15`);
- carry `mu_theta` as the reserved variable-coefficient operand `δU/δθ`, branch-anchored at `chi` on the
  material-advected branch, ⛔ do not construct it from the uniform energy (§3b/§4);
- admissibility is the §3d **background-order (ε⁰)** balance (generalized body force + per-face traction in
  the `𝒮_hold⁰` pairing), ⛔ not the `ε→0` limit of the §3b wave operator, and ⛔ no `W_bg−W_0` insertion.

## Mechanical precedent

`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` is **mechanical
precedent, not authority.** It shows the working shape under this design: an internal tag mapped to a
standard emitted-name string at the emit site, a duplicate emitted name a hard stop, every payload printed as
`NAME: value` in `InputForm` as computed, flushing established on a redirected stream, poles handled
symbolically (e.g. `1/x == 0` rather than a limit that returns `ConditionalExpression`), and
`ConditionalExpression` stripped before emit. Its physics is another step's. ⛔ Do not import or transcribe
any object, symbol value, or result from it — re-derive S11c-b (and the S11c-a substrate it needs) from the
specs alone.

## Run discipline — binds the build's own demonstration runs

- The full derivation, and with it the writing of
  `mathematica/out/S11c_b_brane_operator_mathematica_audit.out`, belongs to the **orchestrator after
  review**. The build's demonstration runs go under the build's own scratch paths, ⛔ never under
  `mathematica/out/`.
- The licence has two seats: run **at most one kernel** at any time.
- Kill criteria for a demonstration kernel, both by PID, both recorded with the task they interrupted and
  reported, ⛔ never answered by narrowing or cheapening the requested object: **600 seconds with no new
  output** (the runtime rule is observable progress, not elapsed time — a visibly streaming task may run
  long), or resident memory past **6 GB**. Run with `--sandbox danger-full-access` (Mathematica).

## Acceptance — executable, no expected values

1. **Blindness / no-hidden-dependency and no-stray-write.** Copy the finished `.wl` alone into an empty
   scratch directory and run it there and in-repo: both runs exit 0, the streams are non-empty and
   byte-identical, and after each run its scratch directory contains exactly the one file that was copied in.
   Any run-time dependence on a repository file, and any write other than stdout, breaks this — and would
   break the blindness control.
2. **Flush under redirected stdout.** Run with stdout redirected to a file; on at least one task the capture
   file must be observed **growing while the kernel is alive**. `Print` alone does not establish flushing on
   a redirected stream. If every task completes too quickly to observe growth, record that in the report —
   ⛔ do not claim the flush obligation was demonstrated.
3. After all demonstrations: `git status --porcelain research/pde_ledger_v3/mathematica/out/` prints nothing.

## Conflicts

Use the spec's §7/§8 (and the §8 report) for any conflict, ambiguity, unavailable construction, or object you
cannot emit under one-tag-per-object. ⛔ Do not fill such a gap with new physics, an expected result, or a
self-review mechanism. Report any ambiguity or non-computable requested object in the §8 builder report
(spec §8), ⛔ not as an invented tag-stream token. ⛔ Do not run the task set as a completeness loop that
hides a per-task failure; §6's emit-and-continue rule (a physics disagreement emits and continues at exit 0;
nonzero exit is operational failure only) is exact.
