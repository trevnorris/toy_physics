# Independent review — S11c-c2 N6 ablation-harness directive, OUTPUT-CONTRACT + WL-BUDGET revision

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_ablation_harness_directive.md`
(a Codex-authored **build directive** = a spec the builder will implement). This is a review **by reading +
cheap mechanical checks**, NOT by execution: the four ablation harnesses this directive specifies **do not
exist yet** — executable ablation of the built harnesses is a separate later build review. Do not try to run
a harness and do not treat their absence as a finding.

## What changed (this is a FOCUSED re-review of a targeted revision)
This directive's knife list was already cleared by a prior two-round review (commit 096d27a9). A build
against it then FAILED on two OUTPUT/BUDGET gaps only, and the directive was revised to fix exactly those:
1. **Output contract → compact.** The harness must print, per variant per certified/DEAD object, only the
   object's **compact PIT fingerprint** (the engine's own `n.pit` GF(p) nonzero-tally / PROBE table) **+ a
   SHA-256 digest** (`n.sha(obj)` for SymPy / `emitted_object_sha256` for WL) — **never** a full symbolic
   object, symbolic difference, arithmetic DAG, PIT sample matrix, or raw engine transcript (including on
   error paths). The prior build dumped full symbolic covariance objects → a ~285 MiB uncommittable
   transcript.
2. **WL budget.** The WL engine has no case/draw knob, so every WL run applies two declared build-scoped,
   **NON-KNIFE** patches to its `/tmp` copy: the driver case iterator (audit.wl:1014) → single case
   `{"LAB_HELD","RHOBR_CONSTANT"}` (matching the three SymPy harnesses' `--anchoring LAB_HELD --density
   RHOBR_CONSTANT` pin), and the adaptive `drawCount` (audit.wl:760) → `drawCount = 4`, each under
   `timeout --kill-after=5 600`, serialized. The prior build's unablated 4-case WL run blew `timeout 600`.

The **knife design is claimed frozen** (every knife scope/site, old fragment, FORM replacement, NO-OP,
COEFFICIENT ×2, DEAD set left byte-identical). Your job is to verify the revision is sound AND did not
disturb the frozen design.

## What you are handed
- The revised directive (path above).
- The exact revision diff: `/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_c2_N6_rereview.diff`
  (PRE → POST; use it to confirm what did and did not change).
- The four live engines the harnesses wrap (read these to check coverage + the budget-patch sites):
  - `research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl`
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py`
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py` (defines `pit` at line 668, `sha` at line 91)
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py`
- The governing standard: `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/feedback_carrier_first_numeric_pit_for_heavy_cas.md`
  (a heavy-CAS ablation harness must print compact PIT fingerprints + digests, never full symbolic).

## Required method
Read the revised directive and the four engines. Verify every claim below **against the engine source**, not
against the directive's own prose. Where a claim is checkable by a `grep`/line read, do it and cite the line.
This is read-only: if you run any shell check, do not modify the working tree.

## Settle these specific questions — give reasoning + any mechanical check for each
**Q1 — Frozen knives untouched.** From the diff, confirm that NO knife block (scope/site, old fragment, FORM
replacement, NO-OP/IDENTITY, COEFFICIENT ×2, DEAD set) for any of the ten knives (WL K_carrier/K_junk/
K_split_route; covariance K_junk/K_circular/K_rank; diagnostic K_EW_rowdrop/K_slotdrop [+K_ewsign];
reconcile K_normal/K_source_route/K_operand_swap) was altered. Report any knife-bearing hunk.

**Q2 — Coverage is real.** The revision claims every certified + DEAD object in every harness is covered by
an engine-emitted compact PIT fingerprint OR a digest, with the ONLY non-PIT objects being `FROZEN_PHI` and
`PHI_DOMAIN_CENSUS` (covered by `n.sha`/`emitted_object_sha256`). Verify against the engines: (a) do the
SymPy engines actually route their certified/DEAD arithmetic objects through `n.pit`, and does `n.pit` emit
a per-object fingerprint with a nonzero tally that a harness can select? (b) does `n.sha` exist and apply to
`FROZEN_PHI`/`PHI_DOMAIN_CENSUS`? (c) does the WL engine emit a per-object PROBE fingerprint + an
`emitted_object_sha256` for every WL primary/DEAD object? Name any certified/DEAD object that is NOT in fact
covered by a compact fingerprint or digest.

**Q3 — Faithful witness (the load-bearing question).** Could a genuine FORM knife bite a certified/DEAD
object WITHOUT moving its compact PIT fingerprint AND without changing its digest — i.e., is there a blind
spot where the compact record fails to witness a real structural change the full-symbolic triple would have
shown? Consider especially: an object whose only witness is the digest (does a FORM change necessarily change
the canonicalized `n.sha` input?); and a PIT nonzero-tally that could stay constant while the underlying
object changes. Conversely, could a NO-OP/IDENTITY or a pure COEFFICIENT change spuriously move a fingerprint
(e.g. seed instability) and be misread as FORM bite?

**Q4 — WL budget fix.** Verify in `audit.wl`: (a) the case-iterator fragment at line 1014 and the
`drawCount = If[...]` fragment at line 760 are each unique and correctly the driver/adaptive-draw sites; (b)
replacing the iterator by the single case and `drawCount` by 4 is genuinely NON-KNIFE — it does not alter the
retained case's per-case symbolic objects, the prime list, branch cells, or any knife site; (c) is
`drawCount = 4` adequate for this harness's COMPARATIVE purpose (baseline vs corrupted fingerprint), given it
overrides the engine's adaptive 2^-80 false-negative draw count — or does 4 draws risk a knife's bite being
invisible? (d) the directive notes `caseOrdinal`-derived PIT seeds shift when the case schedule is restricted
(the retained case moves from ordinal 2 to 1); does this break the baseline-vs-corrupted comparison or any
cross-reference (each variant uses the same restricted schedule)?

**Q5 — No new leak / no new defect.** Does the revised contract actually prevent a symbolic / hundreds-of-MB
blowup, including on failure (the error-path clause: bounded diagnostics only)? Is there any residual path in
the spec where a full symbolic object, symbolic difference, arithmetic DAG, PIT sample matrix, or raw engine
transcript could still reach a committed transcript? Did the revision introduce any ambiguity a builder could
implement wrongly?

## Physics filter
Report a finding only if it catches a way the ablation cert could be wrong or misleading — an uncovered
object, a blind-spot where a knife's bite is invisible, a budget patch that silently changes the physics, a
residual output-volume leak, or a disturbed frozen knife. Do not report "the harness would be wrong on a
different engine" or style preferences.

## Output
A short report: per-question verdict (SOUND / finding), each finding with the directive line + the engine
evidence, and an overall SOUND / NOT-SOUND on whether the revised output contract + WL budget fix are correct
and complete and left the frozen knife design intact.
