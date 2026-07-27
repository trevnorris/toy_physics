# Dimension inventory + quantity provenance — what we know, where it lives, what is open

**Purpose.** A durable pointer map so the dimension work survives a compact. This is a *thin pointer,
not a copy* — every number below is sourced, and the cited file is authoritative.

Written 2026-07-26 on `ledger-v2-rebuild`. Anchor for the rewrite = `manifests/REWRITE_HANDOFF.md`.

> ⚠ **The failure this doc exists to prevent.** A session resumed on the mechanics of the rewrite
> without reading `docs/model_map.md`, `STATUS.md`, or the plan-of-record commits, and began
> re-deriving conclusions the corpus already recorded (the `K_eta` measure split was re-discovered
> from scratch; it is documented at `measure_register_sufficiency.md` §6.5). **Read the maintained
> docs before touching the scripts** — see the READ ORDER below.

---

## 0. READ ORDER

1. `docs/model_map.md` — the master map: throughline, per-sector atlas, earned/calibrated/R1/departure
   ledger, glossary. Read this to hold the model in your head.
2. `STATUS.md` — the canonical "you are here" (its top three blocks; the rest is dated history).
3. `git log -1 --format=%B aae5d389` and `b5527062` — **the plan of record**: why the parallel fanout
   was cancelled, why a shared dimension module was chosen over an escape hatch, and the standing rule
   *"never adjust the process because the corpus is inconvenient."*
4. `manifests/REWRITE_HANDOFF.md` → `notes/rewrite_reference_table.md` — the per-stage rewrite loop.
5. **This file** — the inventory and the open provenance questions.

## 1. ⭐ WHY THE DIMENSION WORK EXISTS (not refactoring hygiene)

The 44-stage manifest fanout is **blocked**: dimension recovery covers only ~16 of 43 scripts while
the schema requires it for all. The user's decision (`b5527062`) was to **fix the corpus, not weaken
the check** — thirteen dimension idioms across 43 scripts *is* drift and the checker was reporting it
accurately. The shared module lets the checker delete `BARE_TUPLE_DIM_ORDER_BY_SHA256`, the AST
bare-tuple recovery, the live module-execution path and the per-script order registry: **one import
means one recovery path.**

Three consumers, not one:
1. the per-stage cross-engine `.py`-vs-`.wl` comparison;
2. the composite checker's single dimension-recovery path;
3. **Part VII's whole-system dimensional firewall** — a named ledger deliverable (`model_map.md` §3.7).

## 2. ⭐ THE INVENTORY ALREADY EXISTS — DO NOT RE-RUN IT

**`notes/measure_register_sufficiency.md`** (144 lines, measured 2026-07-26). Headline facts:

| | |
|---|---|
| Scripts manipulate | **401 `(stage, quantity)` pairs = 226 distinct quantities** |
| Over | **30 of 43** scripts; **13 have no dimension machinery**: 001, 014, 015, 017, 019, 020, 022, 024, 025, 026, 028, 033, 043 |
| Register defines | **145** quantities across 31 stages |
| Register covers | **105 / 226 = 46.5 %** |
| Bases | **6 conventions over 4 axis sets** — full table with loci at §3 |
| Symbolic exponents | **none anywhere**; fractional in **8** scripts (not an 042-only need) |

⚠ The 226/401 counts are **exact** under the doc's stated criterion. The *bucketing* of the 121
misses is a judgement call at the margins (±5). Do not quote the buckets as precise.

**The 121 misses:** 74 intermediate/derived composites never intended as knobs · 23
fields/coordinates/operators the register deliberately omits · **12 "a registered knob under a
different key, *or new*"** ← §4 below · 8 convention-variants · 4 external constants (`G`, `c`).

Companion: `notes/measure_rewrite_feasibility.md`. Verdict of both: the register is a **parameter**
register, not a **quantity** register — a correctness oracle for the 105 it covers, never a seed.

## 3. ⭐ THE UNIFYING INSIGHT — "different dimension" is usually a REDUCTION LEVEL

The model lives in a **4D bulk** and on a **3D brane**. Integrating out a coordinate shifts L-powers.
The same quantity therefore legitimately carries different dimensions at different reduction levels,
and this is **not** drift:

| Level | Measure | Example — the wall packet |
|---|---|---|
| volume density | `a²dw dΩ` ~ L³ | stage016 `K_eta = M L⁻³T⁻²` |
| line density | `dw` ~ L¹ | stage013 `K_eta = M L⁻¹T⁻²` |
| reduced scalar | L⁰ | stage023 `K_eta = M T⁻²` |

All three integrate to `[K] = M T⁻²`. **Verified across four quantities independently** — `K_eta`,
`T_w` (`M L T⁻²` vs `M L⁻¹T⁻²`), `μ_η` (`M L⁻¹` vs `M L⁻³`), `T_Omega` (`M L⁻³T⁻²` vs `M T⁻²`) — each
differing by exactly the L-power its measure predicts. Four independent confirmations of one
explanation. The corpus states it outright at `scripts/…stage016…py:837` / `mathematica/…016….wl:615`
("stage013's line-density relation does not transfer"), credited to a Grok compute-catch, and
`measure_register_sufficiency.md` §6.5 already recorded it.

⛔ **Therefore: never resolve a dimension conflict by renaming the variants apart.** Giving three
names to three reduction levels makes the conflict vanish *and destroys the consistency check*. Worse,
it **hides reduction debt and silently shrinks the irreducible knob count** — stage043's `[40,49]`
headline. The register already encodes this correctly:
- `:179` `K_η` `M L⁻¹T⁻²` **DERIVED**, `= T_w β²`, *"collapses into {T_w, β} … not independently counted"*
- `:170` 023's `{K_eta, T_Omega}` `M T⁻²` **FREE-UNREDUCED, PENDING scalar-reduction — COUNTED as
  reduction debt**, *"NOT identified with the raw 013/017 densities"*

The differing dimension is the **symptom**; the unperformed reduction is the **thing**.

## 4. ⭐ THE OPEN QUESTION CLASS — "registered under a different key, OR NEW"

This is the question the whole naming issue reduces to, and it is **archaeological, not a matter of
opinion**: Codex derived this math, reviewed by Claude/Grok/GLM, so the answer is in the artifacts —
not in anyone's memory. Nobody can answer it from recollection.

Named instances from `measure_register_sufficiency.md` §2/§4/§5:
- **new-or-renamed (12):** `lambda_T`, `eta_a`, `lambda_Sigma`, `E_g`, `qL`/`qh`/`qM`/`mh` (042);
  `kappa_chi`/`lambda_chi` (044) = register `κ_B`/`a_B`; **`A_E` (036/037/038) has no row at all.**
- **one quantity, several names:** `D` (:202) = `sigma` (040:822) · `η_i` (:209) = `eta_a` = `eta`
  (044) · `s_i` (:226) = `s` = `s1`/`s2`.
- **one name, several quantities:** `κ` · `T_Ω`/`T_Omega`/`T̃_Ω` · `A`/`A_eff`/`A_X`/`A_E` ·
  `A_V` colliding across 032 and 036/037 · **`c_E` (`L T⁻¹`) vs `C_E` (`M⁻¹L⁻⁴T²`) — differing only
  by letter case.**
- **`ε0/ε1` (:169) vs `Z0_ret/Z1_ret` (023):** the register calls them *the same two dofs* yet they
  carry different dimensions (`1` vs `M T⁻²`).

### 4.1 ⛔ THE GENUINE ERROR (not a convention) — `r_BA` in two unit systems

> *"Two stages computing the same `r_BA = q_T²c_γ²/(μ_R A_E)` use incompatible unit systems; both land
> dimensionless, so nothing catches it."* — `measure_register_sufficiency.md` §3

- stage038 (`…py:715–719`, four-axis `(M,L,T,E)`): `μ_R = M²L T⁻⁴E⁻¹`, `A_E = L·E`
- stages 003/007/036/037/040/044: `μ_R = M L⁻¹T⁻²`; stage036 `:289` / 037 `:614`: `A_E = M L³T⁻²`

`r_BA` is the open ratio at the heart of the magnetism↔electricity boost-consistency result
(`model_map.md` §3.5). The dimensional gate passes on both sides **because the answer is dimensionless
either way** — a real cross-stage inconsistency wearing a green check. This one is not a reduction
level; it is two unit systems. **Adjudicate it; do not fold it into the convention bucket.**

### 4.2 Register defects found (separate from the conflicts)

- **Locus mis-attribution.** Rows `:182–:185` say "stage 017, dual-engine" — but **stage017 has no
  dimension machinery in either engine** (`…017…py:62` only cites `CITED_016_DIMENSIONAL_OK = True`).
  The values live at `…stage016…py:355–366`. *Values agree, the pointer does not.* Same shape at
  `:174` (`Vp0/ℓ_c`, no engine computes it) and `:181` (`{κ,χ,σ_a,σ_L}`, no exponent triple given).
  ⇒ so stage016's `T_Omega` **is** registered (mislabelled); its `K_eta`/`μ_η`/`T_w` are **not**.
- **No machine key.** The register's key is a prose table cell mixing symbol, gloss, strikethrough
  (`~~λ_Pu~~`), escaped pipes and set-braces. No `quantity_id` field exists.
- **Retired rows** (`λ_Pu` :139, `α_aniso` :159) carry historical dims with no live/retired flag.

## 5. WHERE THE PROVENANCE ANSWERS COULD LIVE

The EM sector was **rethought** for this ledger, so a quantity's current meaning may not match its
first appearance. Search order for any "is this new or renamed?" question:

- `research/pde_ledger_v2/notes/stages/` — per-stage notes; `notes/*_source_map.md` — per-stage
  source-to-claim maps (the most direct provenance evidence that exists)
- `notes/parameter_register.md` (888 lines) + `notes/midway_knob_audit.md`
- `software/em_charge_attribute/` and `software/stage1_solver/reports/` — the `pathA_*` reports the
  ledger stages were folded FROM (⚠ `docs/model_map.md` §5.1 maps stage → source report)
- `software/stage1_solver/decisions/` — numbered decisions (esp. **16**, the `P`-retirement)
- **git**: 517 commits on this branch, 69 mentioning "dimension". Commit messages on this project are
  unusually substantive — treat them as primary sources.
- `_scratch/` dirs (all gitignored, all readable): `./_scratch`, `docs/_scratch`,
  `research/pde_ledger_v2/_scratch`, `software/stage1_solver/_scratch`,
  `software/em_charge_attribute/_scratch`
- `research/` — 19 dirs of prior papers feeding this ledger, incl. `4d_em_fields/`, `em_fields/`,
  `brane_bulk_ontology/`, and the **superseded** 253-stage quarry `research/pde_ledger/`
- `notes/ledger_exclusions_failures_paper_backlog.md` — where retired approaches went

## 6. ⛔ PROTOCOL FOR ANSWERING A PROVENANCE QUESTION

A question of the form *"which existing knob is `lambda_T`?"* is **fabrication-forcing** — it presumes
an answer exists and its only honest failure mode is an invented identification. Ask instead:

1. Find **every** occurrence of the symbol across the search order in §5. Report loci verbatim.
2. Report whether **any source ASSERTS an identification** — quoting the sentence and its file:line.
3. If no source asserts one, the answer is **`NONE_FOUND`**. That is a *correct, expected, useful*
   answer, not a failure. Never propose an identification the corpus does not make.
4. Every claim must be a **quote with file:line**, so the gauntlet can verify without redoing the
   search. Answers are quotes, not conclusions.

Then gauntlet the answers (Codex → Grok → Codex per `docs/development_pipeline.md`).

⚠ Codex adjudicating math it originally wrote is acceptable **only** in this quote-backed form: it is
doing archaeology on committed artifacts, and every claim is independently checkable. It is not
re-deriving, and it must never be asked to.

## 7. WHAT IS SETTLED (2026-07-26)

- Module is **not structurally blocked** by stage038 (4 axes `(M,L,T,E)`) or stage042
  (`(stiffness,L,T)`, `Fraction` exponents, `set()` homogeneity ⇒ needs hashability, present at
  `ledger_dimensions.py:156`). Proven by execution — `scripts/probe_ledger_dimensions_extremes.py`.
- The `.py` emits an axis-labelled dimension **sidecar**; stdout stays byte-identical to the committed
  transcript (a behaviour-preservation gate that holds 2/2 and is worth keeping).
- ⭐ **A `render="tuple"` transcript cannot detect a basis transposition at all** — the tuple branch
  emits exponents in the declared axis order, so relabelling is unobservable. stage011 transposed is
  byte-identical on stdout AND passes all its own assertions. `render="symbolic"` detects it but has
  blind spots (4 of 9 lines at stage004). ⇒ **the axis-labelled cross-engine comparison is the only
  universal gate.** Proven able-to-fail on both stages.
- ⚠ **Cross-engine coverage varies wildly per stage**: stage004 compares 20/20 quantities; **stage011
  compares 2 of 12**, and its entire able-to-fail proof rests on one quantity. The comparator must
  print coverage every time — a 2-of-12 stage must never read as a clean pass.
