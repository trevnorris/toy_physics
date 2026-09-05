# Decision-leg review — S11c-c1 T7 cross-engine comparator BUILD DIRECTIVE (before any build)

You are an independent, adversarial reviewer of an **orchestrator-written build directive**, BEFORE the
comparator is built. This is the rule-7 decision gate: your job is to find every way the directive would
produce a wrong or misleading instrument, or would violate the frozen contract, so it is fixed in one fold
before a builder ever runs. A clean pass with citations is equally useful; a productive review is not the goal,
finding real defects is.

⛔ **A prose judgement is worth nothing here. Ground every finding in a file/line or a literal payload byte you
actually read in the repo.** Where you claim the directive misdescribes a payload shape, a family, or a symbol,
show the `grep`/`awk` command and its literal output. Where you claim it contradicts the contract or the
sibling, quote both sides with line numbers.

Working dir: `/var/projects/toy_physics` (git repo, branch `ledger-v3-rebuild`).

⚠ **OPERATIONAL — the two `.out` payloads are large (~82 MB WL / ~91 MB PY); several single tags hold
multi-megabyte srepr.** ⛔ Do NOT attempt a full CAS parse of them and do NOT read them whole — you will stall.
Use light `grep`/`awk`/`sed`/`cut` to sample: tag inventories (`grep -oE '^PY_S11CC1_[A-Z0-9_]+:'`), payload
heads (`grep -E '^PY_S11CC1_<FAM>:' … | head -1 | cut -c1-400`; WL via `awk '/^WL_S11CC1_<FAM>:/{p=1} p{print}
p&&/… /{…}'` then `cut`). Sampling the key/container structure is enough to check the directive — you do not
need the full expression trees. Keep each shell command bounded (a `head`/`cut` cap) so nothing hangs.

## The artifact under review
`research/pde_ledger_v3/directives/S11c_c1_comparator_build_directive.md`

## Read these FIRST (source of truth), form your own view, THEN read the directive
- ⭐ The FROZEN T7 contract: `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:580-587` (N8) and the
  §7 tag grammar around `:589-602`. This is the authority — the directive must not override it.
- The reconciliation schema it inherits: `research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md:233-253`
  (axis-typed keys; **never blanket-collapse**; integral linearity; ground backgrounds; closed folds are
  reviewed REGISTRY entries, not comparator name-surgery).
- The MECHANICAL BASE the directive says to reuse: `research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py`
  and the closest sibling `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py` +
  `research/pde_ledger_v3/scripts/test_S11c_b_cross_engine_comparator.py`; the transliteration
  `research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py:147` (`mechanical_lower_camel`).
- The REAL committed payloads (⛔ `datalad get` both first; they are git-annex pointers):
  `research/pde_ledger_v3/scripts/out/S11c_c1_bulk_closure_sympy_audit.out` (PY, 63 tags) and
  `research/pde_ledger_v3/mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out` (WL, 51 tags). Sample the
  tag inventory (`grep -oE '^PY_S11CC1_[A-Z0-9_]+:'` / `^WL_S11CC1_`) and a few family payload heads yourself.

## What to check (report each as CLEAN / GAP / INACCURATE / MISSING with exact evidence)
1. **Contract fidelity.** Does the directive faithfully carry N8 (`:580-587`): join by object name + S11c-a
   axis-typed keys; paired residual operands; reject a native boolean as a residual operand; repoint ablation;
   "computes and prints, deciding nothing"; ⛔ no pre-registered representational fold; a branch/regime keyed
   differently is a computed residual, not a fold? Flag anything it omits, weakens, or contradicts.
2. **Rule 5 — no leaked value / no target.** Does any line of the directive tell the builder what a residual
   equals, is expected to be, or should reduce to (zero or otherwise)? Does any "definition of done" or
   acceptance criterion reference an expected value? A builder iterating to exit 0 converges on any target it
   can see — a leaked outcome is a defect even when phrased as "do not canonicalize to it." Quote any hit.
3. **Rule 6 / rule 12 — SEALS intact.** The directive SEALS four physics-bearing reconciliations for post-run
   adjudication: the two-momentum legs (PY opaque `s11cc1_q_out_output/_input` vs WL live `qOut[...]`); the ω
   real-assumption; the DtN operator algebra (noncommutative Symbols vs `OperatorSum`/`OperatorComposition`);
   any parity/regime/face key the engines shape differently. Verify each seal is (a) genuinely present, (b) not
   secretly reconciled elsewhere in the directive, (c) not turned into a prohibition that itself leaks the
   answer. Is any OTHER physics-bearing reconciliation (a fold that could manufacture agreement) smuggled in as
   a "spelling" fold? In particular audit fold-1's head/symbol spelling list: is every entry a genuine
   `mechanical_lower_camel` inverse of a real PY reserved name, or is any of them a physics identification?
4. **Axis typing — no forced merge.** The directive forbids assuming `FACE_n ≡ DIRECTION_n` and forbids merging
   PY's positional `Integer(1)` onto an axis it does not mean. Check against the REAL keys: PY `DTN_OPERATOR`
   key `(LAB_HELD, 1)` vs WL `LAB_HELD`; WL `DTN_KERNEL` `LAB_HELD|FACE_1`; WL `CONTROL_FORM_RESIDUAL`
   `...|DIRECTION_1`; the PERMEABLE parity keys (PY `THICKNESS` vs WL `PARITY_DELTA_W`). Would the directive as
   written cause a builder to (a) silently merge distinct axes to force a join — manufacturing agreement — or
   conversely (b) fail to join cases that are genuinely the same object, drowning the comparison in unpaired
   accounting? Either is a finding.
5. **Census + family coverage.** Verify the "identical 50 non-`LOCAL_` family" census against the streams. Are
   the named per-family VALUE shapes (DtN core, face response, loci, ports, energy, controls, dimensions)
   accurate to the real payloads, or does the directive name a container that is not there / miss one that is?
   Is any family unassigned to an extraction group?
6. **False-agreement / blanket-collapse.** Could following this directive yield a comparator that prints a zero
   residual where the engines actually disagree — via a blanket `X(args)→X`, an over-eager spelling map, an
   AppliedUndef→Symbol arg-strip, or a `raw_control_case` outer-signature used where a real scalar leaf should
   have been reached? The `raw_control_case` pattern is legitimate ONLY where the two outer containers genuinely
   differ; flag any place the directive would let it hide a comparable leaf.
7. **Emission contract + tests.** Is the no-verdict / operands-before-guard / accounting-mandatory / exit-0-on-
   disagreement contract complete and consistent with the sibling? Are the required synthetic tests (injective
   spelling map, μ_θ never globally collapsed, repoint ablation, native-boolean rejected, disagreement-is-exit-
   0) the right battery, and is the repoint ablation specified so that it actually bites (a broken pointer must
   move a residual)?
8. **Anything that would cost a build round.** Ambiguities a builder would resolve wrongly; a level-above miss;
   an instruction that contradicts the mechanical base's actual API.

## Output
A numbered verdict (CLEAN / GAP / INACCURATE / MISSING) per item with exact evidence (file:line or literal
payload bytes + the command that produced them), then a one-line overall: is this directive sound to hand a
builder as-is, or what specific folds must be applied first. ⛔ Do NOT propose expected residual values or
physics reconciliations to bake in — if a reconciliation is physics-bearing, the correct fix is to keep it
SEALED for post-run adjudication, not to pre-register it.
