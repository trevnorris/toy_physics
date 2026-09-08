# Independent re-review — S11c-c2 N6 ablation-harness directive (OUTPUT-CONTRACT + WL-BUDGET revision)

Leg: fresh Claude agent (independent). Method: read-only inspection of the revised directive, the PRE→POST
diff, the four live engines, `c.cas`/`n.sha`, and the two surviving (failed-build) harness deliverables.
Every claim below is checked against engine/harness source, not the directive's prose. The four harnesses do
not exist in final form; no engine was run.

**OVERALL: SOUND.** The output-contract compaction and the WL budget patch are correct and complete for their
stated purpose, the digest+fingerprint witnesses are faithful on both engines, and the frozen ten-knife design
is untouched by the diff. One genuine second-order caveat (seed-schedule false-positive on one-sided DEAD
controls) and two non-blocking clarity nits are recorded; none changes what is computed or may be claimed.

Note on process: I initially suspected a must-fix — the directive names `emitted_object_sha256` from WL
`RUN_PROVENANCE`, and neither string exists in the WL engine (`grep` = 0 hits in `mathematica/`). On checking
the harness deliverable I found this is a HARNESS-computed record, not an engine field, and it is faithful.
The finding is therefore RETRACTED. Details under Q2/Q3.

---

## Q1 — Frozen knives untouched — SOUND

The PRE→POST diff (`_legs/S11c_c2_N6_rereview.diff`) contains exactly six hunks, none knife-bearing:

1. `@@ -12,29` — Global contract clauses 1,4,5,7,8 rewritten (compact-payload contract, temporary-spool
   clause, WL-budget process clause, compact drift-print, compact-guard + error-path clause). Contract text only.
2. `@@ -42,7` — clause-9 tail: "print these named objects as triples" → "as the compact fingerprint/digest
   triples defined in clause 1". Text only.
3. `@@ -59,7` — Deliverables transcript sentence → "KB-scale compact record … never … full symbolic object …".
4. `@@ -67,9` — Harness-1 header paragraph replaced by the two-edit WL budget-patch spec. This whole section
   sits BEFORE the "Primary print set" and the `K_carrier`/`K_junk`/`K_split_route` blocks; no knife touched.
5. `@@ -265,7` — one phrase in the diagnostic control-extraction preamble: "variant triples" → "compact
   fingerprint/digest variant triples". This is the shared control paragraph immediately above the
   `### K_EW_rowdrop` header; the header and the knife block are unchanged.
6. `@@ -483,3` — four CHANGE-LOG bullets appended (output-contract, coverage, WL-budget). Log only.

No hunk overlaps any of the ten knives' scope/site, Old fragment, FORM replacement, NO-OP/IDENTITY,
COEFFICIENT ×2, or DEAD-set lines (WL K_carrier/K_junk/K_split_route; cov K_junk/K_circular/K_rank; diag
K_EW_rowdrop/K_slotdrop[+K_ewsign]; rec K_normal/K_source_route/K_operand_swap). The frozen design is intact.

## Q2 — Coverage is real — SOUND (all four engines)

(a) SymPy certified/DEAD arithmetic objects route through `n.pit`, which emits a per-object selectable
fingerprint:
- diagnostic: every certified/DEAD object is `add(...)`-ed into `objects` (lines 823–890) → `pit(objects,…)`
  at `diagnostic:899`.
- reconcile: `add(...)` → `objects['N6RC_'+name,None]` (247–273) → `n.pit(objects,…)` at `reconcile:276`.
- covariance: `tables` → `objects = {('N6COV_'+name,None):…}` at `covariance:250` → `n.pit` at `252`; the five
  named arithmetic objects (R_COV, SOURCE_ACTUAL, SOURCE_PREDICTED, R_COV_INCREMENT, R_COV_CONTROL_DELTA) are
  all in `tables` (201–212), hence in `objects`.
- `pit` emits per object `emit(name,{'columns':keys,'numerator_denominator':table,
  'nonzero_modular_numerator':[any(row[i][0]!=0 …)]},probe=probe,…)` (`diagnostic:769–770`). The
  `nonzero_modular_numerator` tally is the compact fingerprint the harness selects; the `numerator_denominator`
  sample matrix (the volume) is what clause 1 forbids in the transcript. Correct.

(b) `n.sha` exists at `diagnostic:91` (`hashlib.sha256(sp.srepr(c.cas(e)).encode())`). `c.cas`
(`S11c_c2_selfenergy_fold_sympy_audit.py:63`) handles str/bool/Mapping/list/tuple/set and otherwise
`sp.sympify`, so it canonicalizes nested payloads. FROZEN_PHI/PHI_DOMAIN_CENSUS: `covariance` PHI_DOMAIN_CENSUS
emits `n.sha(mu_e)`/`n.sha(b.DERIVATIVE_MAP)` directly (`covariance:106,114`); FROZEN_PHI (116–125) carries no
engine sha but is coverable by the clause-1 worker fallback (`compute only n.sha(obj)`), which succeeds on its
emitted payload. Both are emitted separately from `objects`, so correctly digest-only (not PIT).

(c) WL: every primary/DEAD numeric object is inserted by `addNumeric` (`audit.wl:875–900,931`) →
`numericObjects` → `probeCase`, which emits per object `{case,"PROBE",cellIndex,prime} -> {valid,rejected,
Length[leafExpressions]}` (`audit.wl:800`) plus a per-object payload with `PROBE_NUMERATORS` (808),
`PROBE_DENOMINATORS` (809), `SAMPLE_INDEX` (810), and `LOCAL:PROBE` carries the sample index (814). The
directive's WL PIT fingerprint (valid/rejected + nonzero_count from PROBE_NUMERATORS + circuit_leaves) is
constructible from exactly these. The two non-PIT objects FROZEN_PHI (`putMeta` at 975) and PHI_DOMAIN_CENSUS
(`putMeta` at 850) are metadata — but the final emission loop `Do[… If[KeyExistsQ[outputMetadata,id],
emit[{family,name},outputMetadata[id]]] …, {name, covNames}]` (`audit.wl:1017–1022`, `covNames` includes both,
lines 32) emits them on the normal path as `WL_S11CC2_N6COV_FROZEN_PHI` / `…_PHI_DOMAIN_CENSUS`
(`standardEmissionName`, 37). The WL digest is the harness-computed `emitted_object_sha256` =
`sha(payload.strip().encode())` over each emitted `WL_S11CC2_*` line (see the surviving harness
`scripts/S11c_c2_N6_ablation_harness_wl.py:187–194`) — i.e. a SHA-256 of the object's full emitted text,
stored in the harness's own `RUN_PROVENANCE`. So both metadata objects are captured (structural delta +
digest), and every WL PIT object additionally gets that digest.

No certified/DEAD object is uncovered by a compact fingerprint or digest. `emitted_object_sha256`/
`RUN_PROVENANCE` are harness artifacts, not engine fields; the directive's wording ("use the corresponding
`emitted_object_sha256` entry from WL `RUN_PROVENANCE`") is terse but resolvable via the harness shape (nit 1
below).

## Q3 — Faithful witness — SOUND, with one seed-schedule caveat

The compact record is **fingerprint + digest**, and the digest is a hash of the object's full payload taken
BEFORE the transcript projection, so it is a faithful value-level witness:
- SymPy: `n.sha(obj)` over the emitted payload (which contains the `numerator_denominator` sample matrix). A
  FORM change from one nonzero value to another — invisible to the per-column `nonzero_modular_numerator`
  tally — still changes the sample values and hence the digest.
- WL: `emitted_object_sha256` hashes the full emitted object text, including the `PROBE_NUMERATORS` sample
  values. A nonzero→nonzero same-sparsity FORM change (invisible to `nonzero_count`) changes the emitted text
  and hence the digest.
So there is no blind spot of the kind Q3 asks about: a genuine FORM bite that leaves the tally/count fixed is
still witnessed by the digest.

Spurious-move direction:
- NO-OP/IDENTITY replaces the exact fragment by itself → byte-identical re-parsed source → identical emitted
  payload → identical fingerprint AND digest. No spurious move.
- COEFFICIENT ×2 is a real arithmetic change and is MEANT to move the digest; that is the companion's purpose
  (arithmetic vs FORM/physics), not a false positive.

**Caveat (genuine, second-order; not blocking).** Both engines evaluate all objects JOINTLY per sample and
reject a draw if ANY root is singular there (`diagnostic:750–755`; WL fixed-stream acceptance,
`audit.wl:774–792`). A FORM knife that changes an object's singularity structure can shift which sample points
are accepted, which changes the sampled values — and therefore the digest — of even a DEAD/unchanged object.
This is a possible FALSE POSITIVE on a one-sided DEAD control (never a false negative), it is inherent to the
joint-rejection sampler rather than introduced by this revision, and the seed-robust `nonzero` tally/count is
unaffected. Recommendation for the later executable build-review: interpret one-sided DEAD controls primarily
via the seed-robust nonzero fingerprint and read any digest move jointly with the emitted rejection counts
(`N6_PIT_REJECTIONS`; WL `REJECTIONS`). Proper venue is the build review, since it needs the harnesses run.

## Q4 — WL budget fix — SOUND

(a) Sites unique and correct. The driver `Do[caseOrdinal++; Check[buildCase[case],Quit[99]]; … , {case,
Tuples[…]}]` is unique (`buildCase[case]` occurs only at `audit.wl:1011`); the iterator literal `{case,
Tuples[{{"LAB_HELD","MATERIAL_ADVECTED"},{"RHO4_CONSTANT","RHOBR_CONSTANT"}}]}` occurs only at 1014. The
SETUP occurrence at 582 is `"CASES" -> Tuples[…]` (no `{case,` prefix) and the directive correctly excludes
it. The adaptive `drawCount = If[0<bound<1, …, 8];` assignment is unique at 760; the other `drawCount` at 726
is the Module-local default `drawCount = 8` (a different fragment), and 771/774/792 are usages. Matching by
the driver-`Do` scope plus the literal `{case,…}`/`If[…]` fragment (directive clause 3) selects exactly one
each.

(b) Genuinely NON-KNIFE. The case edit only changes the argument list handed to the unchanged `buildCase`, so
the retained case `{LAB_HELD,RHOBR_CONSTANT}` runs identically → its per-case symbolic objects are unchanged.
The drawCount edit is inside `probeCase`, AFTER `numericObjects`/`leafExpressions` are compiled
(`audit.wl:728–730`) and AFTER `excluded`/`numeratorDegree`/`bound` are computed (757–759); it changes only
the count of valid PIT samples drawn (used at 771/774/792), not `numericObjects`, `primeList` (fixed at 733),
`branchCells` (fixed at 627), or any knife site. The emitted bounds/`FAMILY_BOUND` remain production-computed
(759, 797–799).

(c) `drawCount = 4` is adequate for the harness's COMPARATIVE (baseline-vs-corrupted) purpose. The engine's
adaptive count targets an ABSOLUTE 2^-80 false-negative for zero-certification; a relative
different-vs-identical decision needs far less. Per-draw Schwartz–Zippel agreement rate ≈ D/(p−1−E) ≈ 1e-6 for
these ~1e9 primes; across 3 primes × branchCells × 4 draws the chance a genuinely different object matches
baseline at every sample is astronomically small, and the faithful text/payload digest (Q3) records the
difference at even a single differing sample. Not a risk to bite-visibility.

(d) The `caseOrdinal`-derived seed shift is handled. `seed = 731921 + 10000·caseOrdinal + 100·cellIndex + …`
(`audit.wl:763`); restricting the schedule moves the retained case from ordinal 2 to 1, so absolute seeds
differ from the old four-case `.out`. But every compared variant (CANONICAL, NO-OP, FORM, ×2) uses the SAME
restricted schedule → the same ordinal 1 → mutually consistent seeds, so the baseline-vs-corrupted comparison
is valid. The ablation harness compares only within its restricted schedule and does not cross-reference the
full-schedule `.out`; the c2 name-join comparator is a separate step. The directive requires each variant to
print its actual seeds (`audit.wl:793–794` supplies `ATOM_SEEDS`). No break.

## Q5 — No new leak / no new defect — SOUND, two clarity nits

The blowup is prevented. Clause 1 excludes the `numerator_denominator` sample matrix and any full symbolic
object/difference/DAG; clause 4 keeps the raw spool in `/tmp` and discards it with the tree; clause 8's failure
clause forbids raw object-bearing stdout on error paths; a transcript byte-count guard is an added backstop.
These target exactly the surviving failed-build leaks the rebuild must remove: the covariance harness's
`print(raw.read_text())` on subprocess failure (`covariance_ablation_harness.py:158`), the WL harness's
`SUBPROCESS_OUTPUT`/`print(raw.read_text())` (`ablation_harness_wl.py:181,196`), and the WL harness's
full-payload `baseline=a,corrupted=b` triple print (`ablation_harness_wl.py:167`) that produced the ~285 MiB
transcript. I found no residual spec path that would route a full symbolic object, symbolic difference,
arithmetic DAG, PIT sample matrix, or raw transcript into a committed transcript once clause 1/4/8 are honored.

Non-blocking clarity nits (neither changes what is computed or may be claimed):
1. Clause 1's "use the corresponding `emitted_object_sha256` entry from WL `RUN_PROVENANCE`" reads as if it
   were an engine-emitted field; it is actually the harness's own record, `emitted_object_sha256[name] =
   sha256(full emitted `WL_S11CC2_*` text)` (`ablation_harness_wl.py:190,193–194`). One explicit sentence
   defining it (hash of the full emitted object text, taken before projection) would remove any ambiguity and
   make the digest's faithfulness self-evident to a from-scratch rebuild.
2. Clause 1's "Hash the exact selected object payload … with `n.sha(obj)`" for SymPy: the raw arithmetic
   objects are DAG `Node`s (not `sympify`-able), so `obj` must be the PARSED emitted payload (dict/list/int,
   which `c.cas` accepts) — or equivalently a plain sha over the emitted JSON. Stating which prevents a builder
   from calling `n.sha` on a `Node` (fails) or reaching for an unavailable symbolic form. Faithful either way.

---

### Verdict
Per-question: Q1 SOUND · Q2 SOUND · Q3 SOUND (seed-schedule DEAD-control caveat for the build review) ·
Q4 SOUND · Q5 SOUND (two clarity nits). Overall **SOUND** — the revised output contract and WL budget fix are
correct and complete and left the frozen knife design intact.
