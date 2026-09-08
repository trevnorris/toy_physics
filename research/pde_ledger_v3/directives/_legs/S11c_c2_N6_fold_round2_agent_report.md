# Independent review leg — S11c-c2 N6 ablation-harness directive, FOLD round 2

**Artifact:** `directives/S11c_c2_N6_ablation_harness_directive.md` (Codex-authored build directive = spec)
**Method:** read directive + fold diff + four live engines; read-only mechanical checks; no execution.
**Overall verdict: SOUND.** The fold resolves F1, F2, and F3, introduces no new defect, and leaves the
ten-knife design frozen and byte-untouched. Two forward-looking NITs for the executable build review (neither
blocks the directive; neither changes what is computed or may be claimed at the directive level).

---

## Q1 — F1 fully resolved, nothing else broken by the repin — SOUND

**WL side.** `actualJunkCase = {"MATERIAL_ADVECTED", "RHO4_CONSTANT"}` at `audit.wl:21` — matches the repin
exactly (directive clause 5 lines 53-54; Harness-1 edit-1 lines 122-124). Under that case all three WL knives bite:

- **K_junk** (`actualJunkCoefficient` 0→1, site `:20`): the *only* `case===` gate is
  `activeJunk = If[case === actualJunkCase, actualJunkCoefficient, 0]` at `:845`. Under the OLD pin
  `{"LAB_HELD","RHOBR_CONSTANT"}` that predicate was FALSE ⇒ `activeJunk=0` regardless of the knife ⇒ K_junk was
  inert (the F1 bug, confirmed). Under the repin `case === actualJunkCase` is TRUE ⇒ `activeJunk = actualJunkCoefficient`
  ⇒ knife 0→1 flows into `materialAmplitude`'s `"MU" -> el[pulled] + jk junkMu eW` at `:247`. **Now live.**
- **K_carrier** (`materialNormalKnife`, site `:17`): consumed at `:860`
  `materialGeometry[anchor, s, materialNormalKnife]` — unconditional in `buildCase`, NOT case-gated; feeds `cm` (`:868`)
  → RC carrier/operand certified objects. Source velocity uses a fresh knife-free builder `materialGeometry[anchor,s,0]`
  at `:861` (as the directive notes). **Live.**
- **K_split_route** (`sourceChannel`, site `:892`): `sourceChannel = joinFaceMaps[buildContraction[cm[#], es[#]-ms[#], response[#], #, False]&]`
  — unconditional in `buildCase`, NOT case-gated; feeds `SPLIT_SUM`/`SPLIT_CHECK`/`SOURCE_CHANNEL` (`:894-899`). **Live.**

So the repin is exactly the unique all-WL-knives-live case, as claimed.

**SymPy side.** The three SymPy pins are correctly left `LAB_HELD/RHOBR_CONSTANT` (H2 `:234`, H3 `:331`, H4 `:434`;
none touched by the diff). No SymPy knife is case-gated to a different case the way WL K_junk was — every SymPy knife
targets a module constant or `run`/factory logic, gated only by the `--anchoring/--density` filter that *selects*
the pinned case:

- covariance **K_junk**: `ACTUAL_JUNK` (`:35`) → passed as `kappa_j` (signature `actual_amplitudes(...,kappa_a,kappa_j)` `:137`; call `:187-188`) → appended amplitude `kappa_j*junk*b.e_W` at `:145`. Live.
- covariance **K_circular** (`:186-188`) and **K_rank** (`image_of` `images[atom]=...` `:92-97`): `run`/`prolonged_phi` logic. Live.
- diagnostic **K_EW_rowdrop** (`face_factory` row build `:311`, bound at `m_rows` `:806`) and **K_slotdrop** (`slots=...` `:788`): `run`/`face_factory`; the `for alpha,rho` loop `:789-791` just filters to the pin. Live.
- reconcile **K_normal** (`n.face_factory` `:91`, confirmed exact old fragment), **K_source_route** (`source=closed_response(...)` `:260`), **K_operand_swap** (`ms[s].get(w,sp.S.Zero)` in `ds` `:211-212`): `run`/`build_material_carrier`. Live.

**No inert/case-gated knife found on any harness's pinned case.** F1 resolved.

## Q2 — F2 fully resolved, no residual leak — SOUND

The REV1 carve-out ("Small scalar residuals (`R_N6`, `R_cov`, `SPLIT_CHECK`, slot/closure guard residuals)…remain
direct payloads") is **deleted** (diff old lines 12-14). REV2 clause 1 (lines 31-35) makes those four object classes
"certified objects like every other: print their compact PIT fingerprint…plus the harness-computed digest, and never
print the engine's raw heavy tag or table," and narrows the surviving direct-payload allowance to "a value that is
already a genuine scalar integer/rational with no `ARITHMETIC` or table payload."

That narrowed allowance cannot be satisfied by any heavy object: `R_N6`, `R_COV`, `SPLIT_CHECK`, and the
slot/closure guard residuals are all `addNumeric`-registered (`audit.wl:897-900`) and carry an `ARITHMETIC` circuit
payload (`:804`), so they are excluded from the scalar allowance by construction. The only named non-PIT objects are
`FROZEN_PHI`/`PHI_DOMAIN_CENSUS` (metadata maps), routed through digest + genuinely-scalar fields only.

No remaining leak path in clauses 1/4/5/8, the Deliverables sentence, or the certified/DEAD lists: clause 4 (48-49)
forbids embedding a full symbolic object/difference/DAG/raw transcript "including on an error path"; clause 8 (70-72)
bounds failure output to "exit code, stderr tail/size, spool SHA-256, and missing tags, never the raw object-bearing
stdout or a symbolic payload"; the Deliverables sentence (99-101) restates the forbidden-content list (now adding
"raw residual table"). F2 resolved.

## Q3 — F3 faithful + implementable — SOUND

**Well-defined on both engines.**
- **WL:** the driver SHAs the full emitted `WL_S11CC2_* = payload` text before compaction. The mechanism already
  exists in the failed harness at `ablation_harness_wl.py:190` (`emitted_hashes[name]=sha(payload.strip()...)`, a plain
  text hash over the RHS), and the payload text includes `ARITHMETIC` (`audit.wl:804`) and `PROBE_NUMERATORS`
  (`:808`). Well-defined.
- **SymPy:** `n.sha` = diagnostic `sha(e)=hashlib.sha256(sp.srepr(c.cas(e)).encode())` (`:91`). `c.cas`
  (selfenergy base `S11c_c2_selfenergy_fold_sympy_audit.py:63-72`) is a *recursive* converter over str/bool/Mapping/
  list-tuple/else→`sp.sympify`, so it cleanly handles a **parsed emitted PIT payload** (nested dict/list of
  ints/rationals/strings). It raises on a raw arithmetic-DAG `Node` (hits the `else: sp.sympify` branch; a custom
  `Node` is not sympify-able) — exactly the failure the directive warns against. So the F3 instruction ("apply `n.sha`
  to the parsed emitted PIT payload, NOT the raw `Node`") is both faithful and implementable.

**Genuine same-support witness (no true blind spot for a live knife).** The emitted `n.pit` payload per object is
`{'columns', 'numerator_denominator': table, 'nonzero_modular_numerator': bitmap}` (`diagnostic:769-770`). The
fingerprint keeps `columns`+bitmap and drops the `numerator_denominator` sample matrix; the digest is over the FULL
parsed payload including that matrix. A FORM change that alters an object's value but preserves its nonzero support
still changes the sampled numerators → the digest moves. WL is analogous (ARITHMETIC/PROBE_NUMERATORS change). For a
live knife acting on a genuinely-affected certified object, digest-unchanged AND fingerprint-unchanged would require a
byte-identical payload, i.e. zero computational effect — which is the inert-knife case Q1 already excludes. No true
blind spot.

The recorded Opus-Q3 caveat is correctly directional: a moved DEAD digest is a possible **false positive** (the WL
joint-rejection sampler's accepted set can shift when a knife changes singularity structure), **never a false
negative**; it is scoped to the executable build review (CHANGE-LOG 606-607). F3 resolved.

## Q4 — no new defect; knives frozen — SOUND

The diff's five hunks touch only: clause 1 output contract (F2/F3), clause 5 WL budget case-pin (F1), the
Deliverables sentence, Harness-1 edit-1 case literal (F1), and the CHANGE-LOG (3 new FOLD-round-2 entries + 1
build-review note). A grep of the diff for every knife name and for `FORM replacement`/`COEFFICIENT ×2`/`DEAD print
set` finds knife names **only in clause-5/CHANGE-LOG prose** (diff lines 83, 132) and no knife-definition block —
the ten-knife scope/site/FORM/NO-OP/COEFFICIENT/DEAD blocks (directive lines 167-535) are byte-untouched. Frozen.

Consistency: clause 1 ↔ coverage note (573-582) ↔ residual-compaction CHANGE-LOG entry (596-598) all agree that
R_N6/R_cov/SPLIT_CHECK/guard residuals route through PIT+digest. Clause 5 ↔ Harness-1 edit-1 ↔ WL-case-liveness
CHANGE-LOG entry (593-595) all agree on the repin and on the three SymPy pins staying put. No budget/output
regression: case cardinality stays 1 (single pinned case, was 1, still 1), `drawCount=4` unchanged, one kernel under
`timeout --kill-after=5 600` unchanged. The independent-certificate framing (WL need not share the SymPy case) is
sound — this ablation harness certifies each engine's own knives on its own pinned case; the cross-engine comparison
is the separate c2 comparator step, not this artifact. The Opus-Q3 sampler caveat is correctly scoped to the build
review and framed false-positive-only. No new ambiguity a builder could implement wrongly.

---

## NITs for the executable BUILD review (do not block the directive)

1. **SymPy sampler has the same false-positive-on-DEAD behavior the caveat names only for WL.** The Opus-Q3 caveat
   cites `audit.wl:774-791`, but the SymPy `pit` sampler (`diagnostic:744-755`) also carries a per-sample `attempt`
   in its seed (`sample_seed=(seed,pi,ci,draw,attempt)`) and retries on singular samples, so a knife that changes
   singularity structure can likewise shift the accepted sample set and move a DEAD object's digest. The safety claim
   ("never a false negative") holds for BOTH engines, and the directive states the principle generally, so no
   directive change is required — the build review should simply apply the same sampler-awareness to SymPy DEAD
   digests.
2. **The referenced failed WL harness dumps full raw stdout on a subprocess error** (`ablation_harness_wl.py:181`
   `out('SUBPROCESS_OUTPUT', ...stdout=raw.read_text()...)` and `:196` `print(raw.read_text())`). This is exactly the
   error-path leak the directive's clauses 4/8 now forbid. The directive is correct; this is only a reminder that the
   rebuild must NOT reproduce that error-path dump when it "reproduces the digest mechanism" from that file (the
   directive means only the `:190` digest, not the `:181/:196` stdout dump). Build-review watch item.

## Verdict

- Q1 SOUND · Q2 SOUND · Q3 SOUND · Q4 SOUND.
- **Overall: SOUND.** F1/F2/F3 resolved; no new defect; ten-knife design frozen and byte-intact.
