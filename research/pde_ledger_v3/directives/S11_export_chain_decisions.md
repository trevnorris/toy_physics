# Decision list — the export chain, before S11's engines are built

## ⛔⛔⛔ BLOCKED BY BOTH LEGS, 2026-08-10. ⛔ THIS IS A RECORD, ⛔ NOT A DIRECTIVE.

⭐ Codex: **8** findings · Grok: **4** · both verdicts **BLOCKED**. ⭐ Every finding below was verified by
the orchestrator at source before being accepted (rule 13). ⛔ No builder may read `E1`–`E7` as decisions.

| # | finding | legs | verified how |
|---|---|---|---|
| **B1** | ⛔⛔ **`E2` + `E4` DESTROY the chain's only cross-step contention.** Under `D5` two steps deriving one object write **one key** and the writer compares them — ⭐ that is what the 3 rows ARE. Producer-scope them and they write **different keys, and nothing compares them.** | both | construction from my own wording; `S9_REWRITE_PLAN:235,253` |
| **B2** | ⛔ **`E4` then never fires.** `S10::root_ordering_d3` vs `S11::root_ordering_d3` are distinct full keys ⇒ ⛔ the `C20` case `E4` exists for is **renamed past, not detected** — and `E4` *does* fire on a legitimate same-object re-derivation. | Codex | construction |
| **B3** | ⛔⛔ **`DEFECT_REGISTER.md:675` ALREADY PRESCRIBED THE FIX** — *"an imported key it is about to write must be compared as an **object**, and a write over a differing object is a finding that fails loudly."* ⚠ `E2`/`E4` **replaced** it with a mechanism proposed inside a review. | orchestrator | opened `:675` |
| **B4** | ⛔ **`M3` is wrong: a dereferenced pointer IS a binding.** `S11:522` performs `LEDGER[name]['dimension_key']` → `LEDGER[that key]`. ⚠ `E3`'s **conclusion** survives (S10 need not be re-keyed first — both legs agree); ⛔ its **justification** does not. | Codex | `S11:522-523` |
| **B5** | ⛔ **`E2`'s byte-identical carry CONTRADICTS any future re-key.** Moving S10's 455 keys must rewrite `rho_br`/`mu_R`'s `dimension_key` — ⛔ which are carried S9-origin rows. | Codex | `M13`'s two `S9→S10` edges |
| **B6** | ⛔ **The 3 rows are BOTH sides, and `E2` picks neither.** A scoped retrofit holds **620** records (165 carried + 455 observed); ⛔ forcing one side erases an operand. | both | 617 = 162 + 452 + 3 |
| **B7** | ⛔ **`E6` is FALSE.** `M10` measured `§8`'s grammar window and generalised to the document. ⭐ `S10_SHARED_PHYSICS.md:197` orders **"Emit `M_A`, `M_B`, and `M_A − M_B`"**; both engines emit `Q2_MATRIX_A/B`; ⛔ **zero** tags carry `M_A`. ⇒ `C19` is a real deviation and the S10 tag-name question is **load-bearing**, ⛔ not deferrable. | both | ⭐ orchestrator re-read `:193-197` and grepped the `.out` |
| **B8** | ⛔⛔ **`M1` LEAKS A COMPUTED S11 RESULT.** Stating *"2 resolve: `rho_br`, `mu_R`; 6 absent"* hands the builder the exact membership of `Q6R_RESOLVED_COEFFICIENTS` / `Q6R_UNRESOLVED_COEFFICIENTS` ⇒ ⛔ hard-codable, and a later upstream change stays invisible. ⭐ `E7` must state the **generic rule only**. | Codex | ⭐ the measurement belonged in the review prompt, ⛔ not the list |
| **B9** | ⛔ **`corroborated_steps` is an agreement claim with no operands in the file that carries it** — and `S11:527-529` **specifies `Q6r` to propagate it.** ⇒ a consumer forwards a claim it cannot check. ⚠ `E5`'s narrow reversal is right for the **run**; ⛔ it does not extend to the merged export. | both | row fields = `class, corroborated_steps, display, step, value, value_kind` |
| **B10** | ⛔ **`E1` has no publication semantics for a partial `MAIN` sweep** (`S11:895`) ⇒ a skipped cell can leave predecessor rows in an export that looks valid. ⚠ Same defect as blocked-list finding 4. | Codex | `S11:893-898` |

⭐⭐ **WHAT SURVIVED:** `E3`'s conclusion (⛔ S10 is **not** re-keyed before S11 — both legs, by measurement),
`E1`'s obligation (⭐ from `D5`, ⛔ not from S11's spec), and `E7`'s necessity (⭐ restated generically).
⇒ ⭐ **Plan items 2 and 4 are refuted and both legs confirmed it. ⛔ Plan item 1 — option E — is refuted too,
by the party that proposed it.** ⚠ That is exactly what putting `U1` in scope was for.

---

**Orchestrator-written, 2026-08-10. ⛔ Nothing here is applied. Two legs before any builder launches.**

⚠ This list **supersedes items 2 and 4** of `S11_naming_and_chain_plan.md`. Those two were written from
documents; the measurements below were taken from the artifacts and refute them. ⭐ The refutations are
themselves in scope for this review — ⛔ they are mine and no leg has seen them.

---

## What was measured — ⭐ every row carries the command that produced it

| # | measurement | result | command |
|---|---|---|---|
| M1 | keys S11's spec tells its PY engine to look up in `S10_exports.LEDGER` | 8 names (`S11:518`); ⛔ **6 are ABSENT** from the ledger; **2 resolve**: `rho_br`, `mu_R` | `python3 -c "…LEDGER.get(n)…"` over the `:518` map |
| M2 | `class` / `step` of those 2 | **`KNOB` / `S9`** — ⛔ neither is an S10 product | same |
| M3 | how S11 reaches any S10-produced row | ⭐ only through the row's stored `dimension_key` pointer (`S11:523`) — ⛔ it spells no S10 key | `S11_SHARED_PHYSICS.md:522-525` |
| M4 | S10-produced keys literally spelled in S11's spec | ⛔ **0** (one `\bphase\b` hit is English prose) | regex of all 455 `step=='S10'` keys against the spec text |
| M5 | rows S10 overwrites | **3**, all `corroborated_steps ('S9','S10')` | `[k for k,v in LEDGER.items() if 'corroborated_steps' in v]` |
| M6 | what the engine emits for those 3 | ⭐ **name · upstream value · downstream value · residual · both classes · both steps · flag**, then `assert` | `S10…sympy_audit.py:2089-2111`; payload at `out/S10….out:4215` |
| M7 | where the predecessor operand lives | ⭐ committed `S9_exports.py`; `S9_REWRITE_PLAN:243` calls cross-step comparison *"a diff between two committed files"* | `sed -n '235,247p' S9_REWRITE_PLAN.md` |
| M8 | export-key authoring surface in S10's engine | **6** `ExportRecord(` sites + `generated_ledger_key():1949` | `grep -c 'ExportRecord(' …` |
| M9 | emitted-tag surface | **106** `emit(` py · **137** `emit[` wl (S10); 99 · 110 (S11) | `grep -c` on the four engines |
| M10 | what S10's spec supplies as tag names | ⭐ a **grammar** and a parallel-sets requirement (`S10:484-501`); ⛔ **no `<QUANTITY>` vocabulary at all** | `sed -n '484,510p' directives/S10_SHARED_PHYSICS.md` |
| M11 | what S11's spec supplies as tag names | ⭐ explicit quantity names throughout `§6`/`§8` | `S11_SHARED_PHYSICS.md` |
| M12 | whether S11's spec requires S11 to WRITE an export | ⛔ **it does not mention one** — only the `Q6r` import | `grep -n "export\|LEDGER" directives/S11_SHARED_PHYSICS.md` |
| M13 | `dimension_key` targets, by (owner step → target step) | S10→S10 **47** · S9→S10 **2** · S9→S9 **1** | counter over `LEDGER` |

⇒ ⛔ **The plan's item 4 fails M1–M4.** *"S11's PY engine binds to S10's export keys, so renaming them later
forces an S11 rebuild"* is false: the binding surface is two S9-origin knob rows plus a pointer the
generated file rewrites atomically.
⇒ ⛔ **The plan's item 2 fails M5–M7.** *"`S9_REWRITE_PLAN#4` instructs the defect — the overwrite destroys
operand A and leaves agreement asserted"* is false in both clauses: `D10` requires both operands and the
residual, the engine emits them, and operand A is committed in the predecessor file.

---

## The decisions

**E1 · S11's SymPy engine writes `scripts/S11_exports.py`.**
⚠ S11's spec is silent (M12); the chain pattern requires it. ⭐ If this is decided NO, `E2`–`E4` are void
and the collision question moves to whichever step first writes again.

**E2 · Every key an export writer CREATES, from S11 onward, is producer-scoped.**
⭐ Imported rows the step does not touch are carried through **byte-identical, under their existing keys**.
⛔ The scope syntax is not yet chosen and must be — it becomes a permanent chain convention.

**E3 · `S10_exports.py` and `S9_exports.py` are NOT re-keyed.**
⭐ Basis: M1–M4, M8. ⚠ Consequence to state plainly: S9's and S10's keys stay flat and first-come, so
`E2` is asymmetric across the chain forever.

**E4 · An export writer FAILS if a key it creates already exists in the imported `LEDGER`.**
⭐ Required even under `E2` — ⛔ naming discipline cannot prevent an accident. ⚠ `D5`'s *"chain integrity is
checkable for free"* cannot catch this: it exempts rows the step touched, and the colliding write **is** a
touch.

**E5 · The three overwritten rows keep their present shape; ⛔ no `(previous, current, residual)` triple is
added to the export.**
⭐ Basis: M5–M7. ⚠ `S11:527-529` requires `Q6r` to emit `corroborated_steps` as provenance, so the field
stays and is read.

**E6 · S10's step record gains the tag-name traversal it currently lacks.**
⚠ The record says `M_A`; the comparator joins `Q2_MATRIX_A`; `Q2_MATRIX` appears in neither build directive.
⭐ Under M10 that is **not** an engine deviating from its spec — the spec never named quantities. ⛔ The
defect is that a reader cannot get from one to the other.

**E7 · `Q6r`'s live lookup must not raise on the six absent names.**
⚠ M1: they are absent from the `LEDGER` entirely, ⛔ not merely missing a `dimension_key`.
⭐ `S11:524-525`'s resolved/unresolved split must be well-defined for that case.

---

## ⛔ What this list does not decide

- `C17`, `C18`, and S10's requirements registers — ⛔ repaired by nothing here.
- Whether S10's **emitted tag names** ever change. ⭐ M9–M11: S11's spec already fixes its own quantity
  names, so the forward discipline holds without touching S10's 106 + 137 call sites.
