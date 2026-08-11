# Decision list — repairing `C17` and `C18` in S11's shared spec

**Orchestrator-written, 2026-08-10. ⛔ Nothing applied. Two legs before any builder launches (rule 7).**

⚠ Both are **OPEN spec defects**, and a shared spec is **physics-bearing**: an error there makes both
engines agree on the same wrong thing, or **disagree for a reason that is not physics**.
⛔ Both engines implement both defects **faithfully**. ⇒ ⛔ **Do not repair an engine for either.**

| | the defect, in one line | locus |
|---|---|---|
| `C17` | ⛔ A stratum's `Q3`/`Q4` rerun is evaluated **at one point**, but emitted under the **component's** scope | `S11:620`, `:627`, `:633` |
| `C18` | ⛔ `_SOLUTION` is specified as *"exactly as your CAS returns it"* ⇒ ⚠ an **engine artifact**, not an object, and the two lists are **not joinable** | `S11:230-251` |

---

## The decisions

**S1 · A stratum's `Q3`/`Q4` objects are computed ON THE STRATUM — ⛔ not at a point of it.**
⭐ They are the objects obtained **under that stratum's defining equations**, by whatever route the CAS
supports (rule 3: ⛔ this pins no derivation path).
⭐ `STRATUM<s>_POINT` remains, ⛔ but as an emitted **witness that the component is non-empty and
§3-allowed** — ⚠ it is no longer the evaluation site.

**S2 · If an engine cannot build `S1`'s object and falls back to evaluating at its point, the tag SAYS SO.**
⭐ Every stratum-scoped `Q3`/`Q4` tag carries an explicit evaluation-scope token distinguishing
*on-the-component* from *at-the-emitted-point*. ⛔ A point-local result may never be emitted under a scope
that reads as the component's.

**S3 · Every locus carries a CANONICAL characterisation both engines compute.**
⚠ Two `STRATUM_ORDERING` lists are today joinable only through the CAS's own presentation.
⭐ The spec pins one engine-independent object over `_EQUATIONS` — ⛔ with its variable order and monomial
order fixed **in the spec**, or it is not canonical. ⭐ `_SOLUTION` stays, as provenance.

**S4 · ⛔ A BRANCH IS NEVER SILENTLY DROPPED.**
⭐ If `_REAL_ADMISSIBLE` excludes a branch, emit **which branch** and the test's operands.
⭐ If the engine cannot decide, emit an explicit **undecided** token — ⛔ never a bare `False`, which is
what makes an omission indistinguishable from a refutation.

**S5 · An undecided admissibility is NOT COMPARABLE.**
⭐ The frozen comparator contract treats it as **not compared** — ⛔ never as agreement, ⛔ never as
disagreement. ⚠ Physics-bearing: ⛔ scoring it either way reports a result no engine established.

---

## ⛔ What this list does not decide

- ⛔ **Any admissibility ALGORITHM.** ⚠ `S11:232-237` already measures that the two CASes have unequal
  real-domain capability ⇒ ⛔ requiring one symmetrically asks an engine for what it cannot do.
- ⛔ **Canonicalising strata by recursive stratification.** ⚠ The register rules it out: cylindrical-
  decomposition-grade work ⛔ no two independently built engines implement identically.
- ⭐ `C19`, `C20`, the export chain ⇒ `S11_export_chain_decisions_v2.md`. ⛔ Not this list.
