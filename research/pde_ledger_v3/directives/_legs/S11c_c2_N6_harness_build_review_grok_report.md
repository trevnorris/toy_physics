I'll treat this as an independent build review of the four N6 ablation harnesses: copy them to `/tmp`, ablate the copies, and check that each knife actually bites the live engine. First I'll load the review method and the artifacts.Next I’ll locate the four harnesses, the cleared knife directive, and the engines they wrap so the ablation copies can be set up without touching the working tree.I'll read the knife directive and the four harness implementations next so the `/tmp` ablations can target live knife sites rather than the transcripts.Harnesses and transcripts are in place and no Wolfram kernel is running. Next I’ll inspect the engines’ knife sites and estimate run times so the `/tmp` ablations stay serialized and inside the 600s budget.No hardcoded digests in the harness sources. I’ll copy everything to `/tmp`, add a label filter on the copies, and start the three SymPy re-runs while the uniqueness and compact-contract checks run.Fragment uniqueness checks out. Next I’ll start the three independent SymPy re-runs and tighten the transcript summarizer so DEAD-object sampler shifts aren’t mistaken for knife bites.The three SymPy re-runs are in flight. I’ll set up the error-path test and inspect the WL transcript’s long lines while those jobs run.A Wolfram kernel is still up after the error-path run. I’ll kill it now so it doesn’t hold a license seat or leak memory.Pristine runs are still on CANONICAL. I’ll stage the engine-corruption copies and a tighter fingerprint summarizer so those ablations can start as soon as the re-runs finish.Form-noop and engine-corrupt copies are ready. I’ll wait for the three pristine SymPy harnesses to finish, then run those ablations. The watcher will fire on completion or failure.CANONICAL, unablated, and the first FORM labels already finished with exit 0. I’ll parse those completed triples now instead of waiting for the full harnesses.Independent re-run SHAs match the committed covariance and reconcile transcripts byte-for-byte on the labels that have finished. I’ll wait for the remaining variants, then run the FORM-noop and engine-wrap ablations.IDENTITY variants are reproducing CANONICAL exactly, and FORM/X2 are biting on the completed knives. I’ll keep waiting for the remaining labels.Covariance and reconcile finished with exit 0. I’ll start their FORM-noop ablations now while diagnostic finishes the last two labels.Diag and both FORM-noop runs finished with exit 0. Next I’ll check that turning FORM into a no-op left fingerprints identical to CANONICAL, then start the remaining ablations.FORM-noop made those variants identical to CANONICAL, so the harness is actually applying the `form` field. Next I’ll run the diagnostic FORM-noop and the live-engine corruption tests.Reconcile engine-corrupt exited 2; the other wrap runs succeeded. I’ll inspect that failure and compare the successful wrap SHAs next.K_split_route/IDENTITY was SIGKILL’d and left an orphan kernel. I’ll kill that kernel, then rerun IDENTITY and X2 with the label filter.K_split_route/X2 matches CANONICAL on every certified object — the coefficient companion is inert. I’ll kill the leftover kernel, then run the WL FORM-noop and engine-wrap checks.Kernels are cleared. The independent re-runs are complete; the report below is from those ablations, not from the harnesses’ self-report.

---

# Independent build review — S11c-c2 N6 ablation harnesses

**Overall: SOUND**, with one named inert coefficient companion.

The four harnesses wrap the live engines, apply the cleared knife list, emit compact PIT+digest triples, and are not fakeable under FORM-noop, one-sided engine corruption, or error-path dumps. Every FORM knife **bites** on its pinned case, including WL `K_junk` on `{MATERIAL_ADVECTED, RHO4_CONSTANT}`.

**Finding (does not make the cert fake):** WL **`K_split_route/X2` is INERT** on that pinned case. Independently re-running `2 (es[#] - ms[#])` produced byte-identical certified digests to CANONICAL. `K_split_route/FORM` (affine `False→True`) still moves `SPLIT_CHECK` via `LOCAL_BARE`. The patch is implemented; the doubled operand is a no-op here (source channel at `affine=False` does not change).

All work was on `/tmp` copies. Working tree untouched.

---

## Method (what was actually re-run)

Copies under `/tmp/n6_grok_review/{pristine,form_noop,engine_corrupt,error_path,wl}/`. Ablation scripts and literal stdout:

| Role | Script | Stdout |
|---|---|---|
| fragment uniqueness | `/tmp/n6_grok_review/scripts/check_fragment_uniqueness.py` | `/tmp/n6_grok_review/out/fragment_uniqueness.txt` |
| full SymPy re-runs | harness copies in `pristine/` | `/tmp/n6_grok_review/out/pristine_{cov,diag,rec}.log` (rc 0) |
| full WL re-run | `wl/` driver | `/tmp/n6_grok_review/out/pristine_wl.log` (rc 2: IDENTITY SIGKILL after FORM succeeded; see below) |
| FORM-noop (form→old) | `form_noop/` | `/tmp/n6_grok_review/out/form_noop_{cov,diag,rec}.log`, `wl_form_noop_junk.log` (all rc 0) |
| engine-site corruption | `engine_corrupt/` | `/tmp/n6_grok_review/out/engcorr_{cov,diag,rec}.log`, `wl_engcorr_junk.log` |
| error-path | stub engine / stub `.wl` | `/tmp/n6_grok_review/out/error_{cov,wl}.log` |
| WL X2-only | `S11CC2_LABELS=K_split_route/X2` | `/tmp/n6_grok_review/out/wl_split_x2.log` (rc 0) |
| bite table | — | `/tmp/n6_grok_review/out/evidence/sympy_bite_table.txt` |

Every WL kernel was under the harness’s own `timeout --kill-after=5 600`. One kernel at a time. A 600 s hit on a later CANONICAL retry was treated as a failed ablation and not raised.

---

## 1. Harness FORM ablation (mandatory)

Turn the knife’s `form` field into the old fragment on a `/tmp` copy, re-run, compare to an independent CANONICAL.

| Harness | Ablated knife | Result |
|---|---|---|
| covariance | `K_junk` form=`ACTUAL_JUNK = sp.Integer(0)` | All 7 certified + DEAD: `digest_equal=1`, empty bitmap/fp. Script: `compare_form_noop.py` on `form_noop_cov.log` |
| diagnostic | `K_EW_rowdrop` form=`ROW` | All certified/native/DEAD: `digest_equal=1`. `form_noop_diag.log` |
| reconcile | `K_normal` form=`NORMAL` | All 13 certified + DEAD: `digest_equal=1`. `form_noop_rec.log` |
| WL | `K_junk` form=`actualJunkCoefficient = 0;` | All 5 certified SHAs **byte-identical to independent CANONICAL**, different from real FORM on the 4 objects that bite. `wl_form_noop_junk.log` |

Combined with the real FORM variants **moving** (below): the harness actually applies the `form` field. It is not a hardcoded triple.

---

## 2. Per-knife bite (independent re-run)

Pinned cases: SymPy `--anchoring LAB_HELD --density RHOBR_CONSTANT`; WL `{MATERIAL_ADVECTED, RHO4_CONSTANT}`.

### Harness 1 — WL (`S11c_c2_N6_ablation_harness_wl.py` + `.wl`)

Independent CANONICAL SHAs **byte-matched** the committed transcript (`R_N6=15d0e5f8…`, `SPLIT_CHECK=7f02d074…`, `R_COV=9156f6fa…`).

| Knife | FORM | IDENTITY | X2 | Verdict |
|---|---|---|---|---|
| **K_carrier** | BITES. `R_N6` digest+`nonzero_count` move; SLOT/CLOSURE digest move. DEAD digest_equal=1; only `circuit_leaves` in fp | exact (all de=1, empty fp) | BITES same pattern as FORM | **BITES** |
| **K_junk** | **BITES** on the pinned junk case. `R_N6` and `R_COV` digest+`nonzero_count`; SPLIT/CLOSURE digest. DEAD de=1, `circuit_leaves` only | exact | BITES same as FORM | **BITES** (prior silent-disable is gone) |
| **K_split_route** | BITES `SPLIT_CHECK` only (`nonzero_count`; R_N6 unmoved). DEAD all de=1, empty fp | **FAILED ablation** (worker `exit_code=-9` / SIGKILL, ~91 kB expected vs 10 633 B stdout). Not a physics reading of IDENTITY | **INERT**: X2-only re-run SHAs **identical to CANONICAL** on all 5 certified objects; FORM’s `SPLIT_CHECK` differs | FORM **BITES**; X2 **INERT** |

`K_junk` FORM vs CANONICAL (from `pristine_wl.log`):

- `R_N6` `15d0e5f8…` → `7f1a477d…`, fp fields `{circuit_leaves, nonzero_count}`
- `R_COV` `9156f6fa…` → `5ea2f43e…`, same
- DEAD carriers/predicted/Φ: digest unchanged

### Harness 2 — covariance

Independent CANONICAL/UNABLATED/`K_junk/FORM` SHAs **byte-matched** committed.

| Knife | FORM | IDENTITY | X2 | Verdict |
|---|---|---|---|---|
| **K_junk** | BITES `R_COV`, `SOURCE_ACTUAL` (**digest-only**), `R_COV_CONTROL_DELTA`. DEAD predicted/Φ unmoved | exact | BITES same | **BITES** |
| **K_circular** | BITES R_COV/actual/predicted/increment/control (all **digest-only**). DEAD `SOURCE_ACTUAL` digest-only (sampler; bitmap unmoved) | exact | BITES predicted/R_COV/increment; DEAD actual unmoved (one-sided) | **BITES** |
| **K_rank** | BITES including `FROZEN_PHI`; `PHI_DOMAIN_CENSUS` unmoved (paths preserved). DEAD actual digest-only | exact | BITES R_COV/predicted/increment/Φ | **BITES** |

Digest-only exhibit: `K_junk/FORM` `SOURCE_ACTUAL` `1e29a292…` → `3e32b52d…`, `nonzero_count=0`, empty bitmap. Full-payload digest moved; compact tally did not.

### Harness 3 — diagnostic

| Knife | FORM | IDENTITY | X2 / extra | Verdict |
|---|---|---|---|---|
| **K_EW_rowdrop** | BITES all certified (MATERIAL+RESIDUAL fp; Eulerian/guards **digest-only**). DEAD μ-recon digest-only (sampler) | exact | BITES MATERIAL/RESIDUAL; Eulerian digest_same; DEAD unmoved | **BITES** |
| **K_slotdrop** | BITES all certified. DEAD μ/face-velocity digest-only (sampler; nonzero fp unmoved) | exact | BITES all certified; DEAD unmoved | **BITES** |
| **K_ewsign** (COEFFICIENT extra) | — | — | BITES MATERIAL+RESIDUAL (MATERIAL **digest-only**); Eulerian unmoved | **BITES** |

### Harness 4 — reconcile

Independent CANONICAL and `K_normal/FORM` SHAs **byte-matched** committed (13/13).

| Knife | FORM | IDENTITY | X2 | Verdict |
|---|---|---|---|---|
| **K_normal** | BITES R_N6/carrier-material/MATERIAL_OPERAND/channels; Eulerian/source DEAD unmoved. SPLIT/SOURCE/CROSS **digest-only** | exact | BITES same | **BITES** (scoped wrapper) |
| **K_source_route** | BITES `SPLIT_CHECK` + `SOURCE_CHANNEL` only; DEAD carriers unmoved | exact | BITES same; SOURCE_CHANNEL digest-only | **BITES** |
| **K_operand_swap** | BITES `SPLIT_CHECK` + `SOURCE_CHANNEL`; DEAD carriers unmoved | exact | BITES same | **BITES** |

---

## 3. Wrap live engine (one-sided engine corruption)

Corrupt the engine knife site in `/tmp`, run harness CANONICAL (no knife patch).

| Engine site | Result |
|---|---|
| cov `ACTUAL_JUNK=1` | 4/7 certified SHAs moved vs clean CANONICAL; **6/7 match** harness `K_junk/FORM` (increment extra-move: same tree also had a reconcile `ms→es` sibling patch — two-sided tree, not a fake wrap) | `engcorr_cov.log` |
| diag `rows.pop('E_W')` | **7/7 match** harness `K_EW_rowdrop/FORM`; all moved vs clean CANONICAL | `engcorr_diag.log` |
| rec `ms→es` | Worker `TypeError`/`NameError` (37 B hashed stderr), all objects MISSING, no leak. Output **moved** (crash vs success) | `engcorr_rec.log` |
| WL `actualJunkCoefficient=1` | **5/5 certified SHAs match** independent `K_junk/FORM`; 4/5 moved vs CANONICAL (SLOT unmoved, same as FORM) | `wl_engcorr_junk.log` |

No certified digest is a hand-typed constant in the harness sources (no 64-hex literals; `R_N6` is only a tag name). WL hashes full `WL_S11CC2_* = payload` assignment bytes; SymPy uses `n.sha` on parsed PIT data before dropping `numerator_denominator`.

---

## 4. Compact-contract / error paths

Forbidden tokens (`numerator_denominator`, `SparseArray[`, `ARITHMETIC`, raw DAG) **count 0** in committed transcripts and in every independent `/tmp` transcript.

| Transcript | Size | Compact? |
|---|---|---|
| WL committed / independent | ~1.02 MB / 999 693 B log | **Appropriately compact.** Bulk is interned `PIT_FINGERPRINT` (38 schemas × 51 `{valid,rejected,nonzero_count,circuit_leaves,numerator_count}` cells) + interned `COMPONENT_KEYS` (one 117 137 B axis catalog). Numeric fingerprints, not raw `SparseArray`/circuits. Could intern COMPONENT_KEYS more aggressively; not an object leak. Markdown wrapper (`Working directory` / invocation) sits outside ````text`. |
| covariance | 195–200 kB | compact |
| diagnostic | 252–258 kB | compact |
| reconcile | 268–274 kB | compact |

**Forced error paths** (stub engine prints `LEAK_OBJECT numerator_denominator SparseArray…` then `sys.exit(3)`; stub `.wl` prints the same and `Exit[3]`):

- Transcripts **0** hits on those strings.
- Only hashed `engine.stdout` (16 150 B / similar) + `exit_code` 2 + `MISSING` tags.
- Paths: `/tmp/n6_grok_review/out/error_cov.log`, `error_wl.log`.

---

## 5. Digest faithfulness (same-support)

Digest is **not** a re-hash of the compact tally. Independent exhibits where digest moves and nonzero/bitmap do not:

- **covariance `K_junk/FORM` `SOURCE_ACTUAL`**: digest `1e29a292…`→`3e32b52d…`, `nonzero_count=0`, empty `bitmap_delta`.
- **covariance `K_circular/FORM`**: five certified arithmetic objects digest-only.
- **reconcile `K_normal/FORM`**: `SPLIT_CHECK`, `SOURCE_CHANNEL`, `CROSS_CHANNEL` digest-only.

WL same-support: `K_carrier/FORM` `SPLIT_CHECK` digest moves with fp field `circuit_leaves` only (not `nonzero_count`).

---

## 6. Sampler-awareness on DEAD

Joint-rejection PIT (WL `audit.wl` 774–791; SymPy diagnostic 744–755) can move a DEAD **digest** or shared fingerprint metadata when singularity structure changes. **Never a false negative.**

Observed DEAD moves, all **digest-only or `circuit_leaves` only** — nonzero tally unmoved:

- WL `K_carrier`/`K_junk` FORM: DEAD digest_equal=1; fp only `circuit_leaves` (shared circuit catalog, not object support).
- covariance `K_circular`/`K_rank` FORM: DEAD `SOURCE_ACTUAL` digest-only.
- diagnostic `K_EW_rowdrop`/`K_slotdrop` FORM: DEAD μ-recon / face-velocity digest-only.

**Not a leak into DEAD objects.** Nonzero fingerprint stayed put. X2 variants generally left DEAD fully unmoved (coefficient does not change singularity the same way).

---

## Operational notes (not physics findings)

- WL `K_split_route/IDENTITY` in the full 11-label run: worker `exit_code=-9`, 10 633 B stdout (held-guards only). A CANONICAL retry later hit the **600 s kernel timeout** (silent emit-at-end). Both reported as **failed ablations**; timeout was not raised. IDENTITY is mechanically `old→old` (`ident_eq=True` in uniqueness check); every other IDENTITY variant that completed was exact.
- Reconcile engine-corrupt crashed (`TypeError`/`NameError`); bounded diagnostics only.
- Leftover kernels after `-9` were SIGTERM/SIGKILL’d; none left at end of this review.

---

## Verdict table (10 knives)

| # | Knife | Harness | FORM | X2 | IDENTITY | Fake? |
|---|---|---|---|---|---|---|
| 1 | K_carrier | WL | BITES | BITES | exact | no |
| 2 | K_junk | WL | **BITES** (pinned junk case) | BITES | exact | no |
| 3 | K_split_route | WL | BITES `SPLIT_CHECK` | **INERT** | ablation failed (SIGKILL) | patch real; X2 operand dead |
| 4 | K_junk | covariance | BITES | BITES | exact | no |
| 5 | K_circular | covariance | BITES | BITES | exact | no |
| 6 | K_rank | covariance | BITES | BITES | exact | no |
| 7 | K_EW_rowdrop | diagnostic | BITES | BITES | exact | no |
| 8 | K_slotdrop | diagnostic | BITES | BITES | exact | no |
| 9 | K_normal | reconcile | BITES | BITES | exact | no |
| 10 | K_source_route + K_operand_swap | reconcile | both BITES | both BITES | both exact | no |

Plus diagnostic `K_ewsign` COEFFICIENT extra: BITES.

---

## SOUND / NOT-SOUND

**SOUND.** The four harnesses are faithful, non-fakeable ablation certs of the cleared knife list: they wrap the live engines, FORM knives bite (including WL `K_junk` on `actualJunkCase`), NO-OP/IDENTITY reproduces CANONICAL where the ablation completed, digests witness same-support changes, transcripts stay compact, error paths do not dump objects.

**Named inert companion:** WL `K_split_route/X2` — coefficient ×2 of `es-ms` is byte-identical to CANONICAL on this pin; FORM still bites via the affine route selector. That is case liveness of the coefficient channel, not an unimplemented or hardcoded cert.
