# Independent review — the two S11b per-engine BUILD DIRECTIVES

You are one of two independent legs (the other is a different engine; you will not see its output). These
are **orchestrator-written** documents. Review them **read-only** — ⛔ modify nothing in the working tree.

## Artifacts under review (both, at HEAD of branch `ledger-v3-rebuild`)

1. `research/pde_ledger_v3/directives/S11b_sympy_build_directive.md`
2. `research/pde_ledger_v3/directives/S11b_wl_build_directive.md`

They are the two per-engine build directives for the unified S11b step ("the interface coupling law"). Each
is meant to be a **thin wrapper**: it points the engine at the shared physics spec (which is already
reviewed and committed) and adds **only** the per-engine chain wiring. The shared spec's physics is **not**
under review here and ⛔ must not be re-opened — your job is the wiring, the provenance, the pointers, and
whether anything leaks or misstates.

## What to check — the question

Does either directive state anything **false, leaked, mis-pointed, weaker-than-the-obligation, or
inconsistent with the other directive or with the authorities below?** Report a finding only if it catches a
way the **build** could go wrong (a wrong import, a leaked value, a blindness break, a mis-scoped task, an
inherited obligation silently weakened). ⛔ Do not report "a directive could be worded more nicely."

## What you are handed (read the SOURCES OF TRUTH first, then the two directives)

Form your own view of what each authority requires **before** reading the directives, then quote **both
sides** (directive `file:line` and source `file:line`) for every finding.

- **The shared spec** `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` — the sole physics
  authority both engines wrap. Its §11 (lines ~947–988) already states the intended downstream wiring and
  explicitly **defers the exact per-engine wiring to these directives**; its §10 is the tag grammar; §8 is
  the three clauses / no-verdict / locus protocol; §9 is the task list (B0–B9); §12 is supplied-vs-tested.
- **The decision list** `research/pde_ledger_v3/directives/S11b_unified_decisions.md` — G1–G14. The wiring
  authorities are **G7** (chain mechanism inherited: F1–F9, the blind-Wolfram control, imports, the
  `c_s0`/`μ_R`/`ρ_br⁰` bindings, the digest pin), **G8** (three named deviations: (a) comparator to the
  frozen T7 contract, net-new and separate; (b) restore D3; (c) include `_RELATIONALS`), **G9** (the
  denylist is CUT — the only blindness control is the blind Wolfram engine), **G10** (print objects +
  residuals, no VERDICT), **G12** (required fixes, incl. G12a: does anything presuppose the FORM of an
  answer?).
- **The chain rules** `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` — **F1** (flat
  keys), **F3** (a re-derived row carries its own evidence), **F6** (publish only if every declared primary
  cell completed), **F9** (the three-valued collision: F9a/F9b/F9c; and note line ~203 flags the
  `s11_`→`<step>_` prefix generalization as *not decided* by F9), and the **OWED TO THE BUILD REVIEW**
  checklist (line ~210, which belongs to the review legs, ⛔ not to engine guards).
- **The export-shape rules** `research/pde_ledger_v3/S9_REWRITE_PLAN.md` — **D1** (line ~210: export what
  the primary derivations emit; controls are ablation evidence), **D3** (line ~217: reconstruct-and-compare
  round-trip in the run), **D5** (flat).
- **The frozen comparator contract** `research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md:53-60`
  — the **T7** contract the (separate, later) comparator must meet.
- **The upstream export the SymPy engine imports** `research/pde_ledger_v3/scripts/S11_exports.py` — confirm
  it carries rows keyed `c_s0`, `mu_R`, `rho_br` (and their shapes: `display`/`value`/`value_kind`/`class`/
  `step`/optional `dimension_key`/`f9_operands`/`corroborated_steps`), the `_restore`/`_RELATIONALS`
  reviver, the `BUILD_INPUT_DIGESTS` self-pin, and the `MappingProxyType` + `del _LEDGER` freeze footer.
- **The precedent directives** `research/pde_ledger_v3/directives/S11_sympy_build_directive.md` and
  `research/pde_ledger_v3/directives/S11_wl_build_directive.md` — the pattern these two extend. The S11
  SymPy engine **dropped D3** and shipped **no comparator**; the S11b directives are supposed to restore D3
  and defer the comparator to a separate artifact.

## Required method — DOCUMENT review

Read each authority, form your own view of what it requires, then read the two directives and quote both
sides. For **every** citation a directive makes (`file:line`, a rule name, a row key, a section number),
**open the cited source and confirm it says what the directive claims.** A citation you did not verify is
not reviewed. Specifically decide, per directive:

**SymPy directive:**
- Are the inherited-object bindings right and faithful to spec §11 — `c_s0`→`LEDGER['c_s0']`,
  `μ_R`→`LEDGER['mu_R']`, `ρ_br⁰`→`LEDGER['rho_br']` (settled identity, not engine-decided), and do
  `ρ_m`/`v_dr` correctly originate here with no upstream row? Is anything re-declared under a new identity?
- Is the candidate population correct under **D1** (primary derivations B0–B7,B9 export; **B8 controls are
  ablation evidence, not exports**; `_LOCAL_` tags engine-local, not exports)? Is any primary/control/local
  boundary drawn wrong against §9/§10?
- Are **F1/F9/F3** *pointed at* (not restated, no tag-to-key transform specified)? Is the **F9c prefix =
  `s11b_`** resolution legitimate, or an over-reach beyond what F1/F9 license? Is the three-valued F9
  rendered correctly (F9b carries **no** residual — the operand pair is the object)?
- Does `BUILD_INPUT_DIGESTS` pin **exactly** `{this engine's own source, S11_exports.py,
  S11b_SHARED_PHYSICS.md}` — the three the spec §11 names — and is the freeze pattern right?
- Are the **three G8 deviations** rendered correctly and completely: comparator **separate/not built here**
  (T7, net-new); **D3 restored** (S11 dropped it); **`_RELATIONALS` included** (relations/inequalities
  round-trip)?
- Does it wrongly encode the **OWED TO THE BUILD REVIEW** checklist as engine guards, or make a **check
  audit its own input**? Does it add any expected value or acceptance criterion (rule 5 violation)?

**Wolfram directive:**
- Is blindness enforced **by absence** (imports nothing, re-derives from §1–§6, emits no import tag) and
  **not** by any do-not-read list (rule 12 / G9)? Flag ANY sentence that forbids a read, or any residual
  denylist.
- Is the no-VERDICT rule correct (a boolean test emitted as the CAS object, never a native boolean)?
- Is the run discipline sound (one kernel / two seats; kill on 600 s no-new-output or >6 GB; `out/`
  replacement is the orchestrator's post-review) and are the acceptance items **able to fail** (esp. the
  scratch-dir isolation as a genuine blindness/no-stray-write probe, and the flush observation)?

**Both:**
- Do they wrap the **same** spec §§0–13 consistently, with consistent, non-colliding engine/exports/`out`
  filenames (`S11b_interface_coupling_law_*`)?
- **Trap 1 — false-equal by name reuse:** anywhere a name is reused across a semantic boundary such that
  F9's object comparison could prove a **false equal** (the `v_0` hazard, decision list G3 — S11's `v_0` is
  a different quantity; S11b's drain is `v_bulk_normal_0`). Does either directive risk reintroducing it?
- **Trap 2 — restated, not pointed:** anywhere a directive **restates** an inherited obligation (an F-rule,
  a D-rule, a §8 clause, the blindness control) in **weaker words** instead of pointing at it — a re-wording
  comes out weaker and drifts.
- **Trap 3 — presupposed form (G12a):** does any wiring instruction presuppose the FORM or value of an
  answer the engines are supposed to compute?

## Physics filter

Report a finding only if it catches a way the build could be wrong. ⛔ Do not report "this could be wrong on
some other step" or style.

## Independent-derivation note

This is a document review, so there is no CAS ablation to run and no kernel to spawn — parallel is fine, no
timeout wrapper, ⛔ modify nothing. But **every cited `file:line` / rule / row key must be opened and
confirmed**; a finding resting on an unverified citation is discarded, and a claim that a directive is
"clear" without having checked its citations against the sources is weak evidence.

## Output

Per directive (SymPy, then Wolfram): **CLEAR**, or findings — each as `defect · directive file:line ·
source file:line · the fix`. Then the three cross-cutting traps: CLEAR or finding. End with a one-line
verdict per directive: is it safe to build from, or is there something that must be fixed first?
