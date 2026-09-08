# Independent build review — S11c-c2 N6 ablation harnesses (executable, ablate the harness)

## Artifacts (four astra-written ablation harnesses + their transcripts)
- WL:        `research/pde_ledger_v3/scripts/S11c_c2_N6_ablation_harness_wl.py` (+ `mathematica/S11c_c2_N6_ablation_harness.wl`)
             → transcript `_measurements/S11c_c2_N6_ablation_harness_wl.md`
- covariance: `research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_ablation_harness.py`
             → `_measurements/S11c_c2_N6_covariance_ablation_harness.md`
- diagnostic: `research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_ablation_harness.py`
             → `_measurements/S11c_c2_N6_diagnostic_ablation_harness.md`
- reconcile:  `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_ablation_harness.py`
             → `_measurements/S11c_c2_N6_reconcile_ablation_harness.md`

## Role
Each harness is a committed ablation CERT (the S11c-b harness STANDARD): it WRAPS the live N6 engine, applies
the fixed FORM knives from the cleared directive `directives/S11c_c2_N6_ablation_harness_directive.md`, and
PRINTs, per variant per certified/DEAD object, a compact `{baseline, corrupted, diff}` = **PIT fingerprint +
harness-computed full-payload SHA-256 digest** — ⛔ never a full symbolic object / sample matrix / arithmetic
DAG / raw residual table / raw engine stdout. The directive (knife design) is already review-cleared; THIS
review checks the BUILT harness faithfully IMPLEMENTS it and that the cert is not fakeable.

## The engines the harnesses wrap (read to derive independently)
`research/pde_ledger_v3/scripts/S11c_c2_N6_{covariance,diagnostic,reconcile}_sympy.py` and
`research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl`. The cleared knife list (10 knives) is in
the directive; the WL harness is pinned to case `{MATERIAL_ADVECTED, RHO4_CONSTANT}` (= actualJunkCase, the
unique all-WL-knives-live case), the three SymPy harnesses to `--anchoring LAB_HELD --density RHOBR_CONSTANT`.

## Required method — ABLATE THE HARNESS, ⛔ do NOT trust its self-report
⭐⭐ **The harness's own transcript is its self-report; verify it by ABLATING the harness and re-running, never
by reading it.** For EACH harness, copy it to `/tmp` and ablate the COPY (⛔ never the working tree):

1. **MANDATORY FORM ablation of the harness itself.** Break a load-bearing step — e.g. make the harness NOT
   apply a knife's FORM replacement (turn the FORM variant into a NO-OP), or corrupt the site the knife
   targets — then re-run and report the LITERAL transcript diff. If the harness's FORM-variant fingerprint or
   digest is BYTE-IDENTICAL to CANONICAL/NO-OP when it should move, the knife is not actually implemented (a
   dead cert). Conversely a NO-OP/IDENTITY variant MUST reproduce the CANONICAL fingerprint+digest exactly.
2. **Does the knife actually BITE?** For each of the 10 knives, confirm from a re-run that the FORM variant
   MOVES the certified object's PIT fingerprint and/or full-payload digest vs CANONICAL, and that the
   COEFFICIENT ×2 companion moves it too (arithmetic), while NO-OP/IDENTITY leaves it identical. **Name any
   knife that is INERT on its harness's pinned case** (this is exactly the defect a prior round caught: WL
   K_junk was silently disabled by the case pin). Check WL K_junk specifically bites now on
   `{MATERIAL_ADVECTED, RHO4_CONSTANT}`.
3. **Does the harness WRAP the LIVE engine (not fake / hardcode)?** One-sided corruption: break the ENGINE's
   knife site (in a /tmp engine copy) vs break the HARNESS's application — confirm the harness reads the real
   engine output (if the harness output does not move when the engine's construction changes, it is not
   wrapping the live engine). Confirm no certified fingerprint/digest is a hand-typed constant (e.g. an
   R_N6→18/288 value baked into the harness rather than computed from the engine).
4. **Compact-contract / no leak (incl. error paths).** Confirm NO transcript contains a full symbolic object,
   symbolic difference, arithmetic DAG, PIT sample matrix (`numerator_denominator`), `ARITHMETIC` circuit,
   raw `SparseArray`, or raw engine stdout — and that a FORCED error path (e.g. make a worker exit nonzero)
   prints only bounded diagnostics, never a raw object-bearing dump. The WL transcript is ~1 MB: judge whether
   that is appropriately-compact numeric fingerprints (many objects × cells × primes × variants) or
   over-verbose (a fingerprint that should be summarized further) — report which.
5. **Digest faithfulness.** Confirm the per-object digest is the harness SHA over the FULL emitted payload
   before compaction (so a same-support FORM change with unchanged nonzero tally still moves the digest), and
   is NOT just a re-hash of the compact fingerprint. Exhibit a knife whose bite is witnessed ONLY by the
   digest (same-support), and confirm the digest moves.
6. **Sampler-awareness on DEAD controls (BOTH engines).** The joint-rejection PIT sampler (WL audit.wl:774–791
   AND SymPy diagnostic_sympy.py:744–755) can shift the accepted sample set when a knife changes singularity
   structure, so a DEAD/one-sided control's DIGEST can move even when its object is unchanged — a possible
   false POSITIVE on DEAD controls, ⛔ never a false negative. When you see a moved DEAD digest, distinguish a
   real leak of the knife into a DEAD object from this sampler artifact (check whether the nonzero fingerprint
   — the robust channel — moved, and whether the accepted-sample set changed).

## Operational constraints (BOTH legs get these — they are method, not wrapper)
⛔ The WL harness spawns Mathematica kernels. Wrap EVERY kernel run in `timeout 600`; a 600 s hit is a FAILED
ablation — report it and move on; ⛔ never raise the timeout; ⛔ never run more than one kernel at a time (the
licence has TWO seats). ⛔ Copy each harness (and any engine you corrupt) to `/tmp` and ablate the COPY; ⛔
never modify the working tree. ⭐ Save every ablation script AND its literal stdout to named absolute paths and
report those paths — ⛔ a prose "I re-ran it and it moved" is discarded; show the script and the literal diff.

## Physics filter
Report a finding only if it catches a way the CERT could be wrong or misleading: an inert/unimplemented knife,
a hardcoded/faked value, a harness that does not wrap the live engine, a same-support bite the digest misses,
an output leak (incl. error paths), or a DEAD false-positive misread as a bite. ⛔ No style preferences.

## Output
Per-harness, per-knife verdict (bites / inert / faked) with the literal ablation diff + script path; the
compact-contract verdict per transcript; and an overall SOUND / NOT-SOUND on whether the four harnesses are
faithful, non-fakeable ablation certs implementing the cleared knife list.
