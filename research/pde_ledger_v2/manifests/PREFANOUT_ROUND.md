# Pre-fanout consolidated round (#27) — DIRECTIVE

The last prep before the 44-stage fanout. Two phases: (1) fix the two checker gaps the
pilot surfaced; (2) rework the pilot manifests (export-all + `.out` evidence + cleanups)
against the fixed checker, then re-verify. Project root `/var/projects/toy_physics`;
work in `research/pde_ledger_v2/`. HEAD at authoring = `739ada37`.

Process: **CODEX IS THE CODER** — the standing calibrated contract (`docs/development_pipeline.md`
§Roles). Codex designs + writes + runs every script and iterates to exit 0; Claude reviews and
owns directive prose only; agents are REVIEW instruments and never write code. Every review leg =
a FRESH agent. Independently re-run every script (never trust a self-report). The checker is
FROZEN except during Phase 1 (which changes it deliberately).

> ⚠ CORRECTION (2026-07-24): this file previously declared "agent-as-coder is the default,"
> citing two Codex CLI stalls. That claim did NOT survive verification — all 27 Codex runs from
> that day end `___CODEX_BUILD_DONE___(exit=0)` (five of them stage031 extractions), and a fresh
> smoke run returned a correct, self-counted answer in 20s on codex-cli 0.145.0 (current). There
> was no Codex fault. Do not re-derive the pivot from this file's history.

---
## PHASE 1 — fix the two checker gaps (dual-fixture discipline, then Grok)
Checker = `manifests/composite_build.py` (currently 49 self-test fixtures, frozen SHA
`faa7e8f1de1d3755`). Both fixes MUST keep ALL 49 existing fixtures green AND keep real
stage030/031/032 passing — pin BOTH with fixtures so it can't oscillate (this is the
convergence discipline that worked in #25/#26).

**Gap A — `pi`/NumberSymbol dimensionless.** `dimension_of` doesn't treat sympy `pi`
(NumberSymbol, `pi.is_Number` is False) as dimensionless → π in a dimensional expression
hits `DIMENSION_RULE_UNSUPPORTED`. FIX: treat sympy NumberSymbols (`pi`, `E`, …) as
dimensionless in `dimension_of`. Fixtures: (i) a dim expr containing `4*pi*R**2` certifies
(π dimensionless) → PASS; (ii) π does NOT mask a real dim break (a genuinely inhomogeneous
π-bearing expr still → `DIMENSIONAL_INHOMOGENEITY` FAIL).

**Gap B — bare-tuple dim definitions.** The dim-recovery requires a live `Dim` class, but
some scripts (e.g. stage032) define dims as bare tuples `dim_E=(2,-2,1)` — so 032's C4
certificate had to BORROW stage031's script (mis-anchored). FIX: extend recovery to also
read module-level bare-tuple dim assignments (`dim_X = (l,t,m)` literals, order per the
script's own convention/docstring) so a symbol's `dim_source` can anchor to its OWN script.
Fixtures: (i) a symbol dim_source'd to a bare-tuple def certifies → PASS; (ii) a wrong
declared tuple vs the bare-tuple source → `DIM_SOURCE_TUPLE_MISMATCH` FAIL.

**Gap C (schema, small) — `.out` evidence field.** To cite saved Mathematica output, the
`verification` block (or per-claim mathematica evidence) needs an optional field for the
`.out` path + sha256 (today it only has `mathematica_script` = the `.wl`). Add the minimal
optional field to the schema + a check that, when present, the `.out` digest is live
(recompute from `mathematica/out/…`). Fixture: a stale `.out` digest → FAIL.

Acceptance P1: `python3 manifests/composite_build.py --self-test` → all fixtures PASS incl.
the new ones; a smoke shows stage032's dims certify against its OWN script; real
{030,031,032} still clean. THEN **one Grok adversarial pass** on the checker delta (holes /
can't-fail fixtures), fold, re-green. Re-freeze; record the new SHA.

---
## PHASE 2 — rework the pilot manifests against the fixed checker, then re-verify
**Protocol change first:** in `manifests/EXTRACTION_PROTOCOL.md` rule 7, replace "export
operative claim ids" with **"export EVERY operative claim AND every ownership (`declare_*`)
claim; only retired/superseded/departed claims stay unexported"** (user-approved 2026-07-24;
kills the under-export class deterministically). This becomes the fanout default too.

Then rework stage030/031/032 (Codex is the coder; the checker from P1 is frozen again):
1. **Export-all:** re-export each of 030/031/032 to include every operative + ownership
   claim. This makes 031 export `S_gg`'s binding claim (`declare_sgg`) and everything else
   downstream might cite.
2. **Un-workaround 032:** now that `S_gg` is exported, replace 032's opaque `C_V` fold with
   the real cited form `A_V = m_gg·φ²/S_gg²` (drop the opaque-numerator representation).
   Re-anchor 032's every `dim_source` to its OWN script (bare-tuple, now recoverable in P1)
   — no more borrowing stage031's `Dim` class.
3. **Drop the pi workaround:** remove stage031's `declare_pi` dimensionless-local symbol +
   claim (the checker handles π now); ensure no double-definition.
4. **Fix the det_m narrative:** stage031's `extraction.report.prose_script_mismatches[0]`
   falsely says "v1 det_m was M²T⁻⁴; corrected" — v1 was `[L,M,T]` axis order (already
   M⁻²T⁴), a re-order not an error. Correct that prose (value is already right).
5. **Cite the `.out` evidence:** for each of 030/031/032 (and, if cheap, all extracted),
   add the `.out` path + recomputed sha256 to the Mathematica evidence (Phase-1 Gap-C field).

Re-verify (independent): `python3 manifests/composite_build.py` on {030,031,032} →
clean (PARTIAL only on honest edges into unextracted stages; zero FAIL); the C7 `zero_mode`
edge STILL PASS (c7_edges=1/17); a FRESH fidelity agent spot-checks each reworked manifest
(esp. 032's un-opaqued `A_V` + re-anchored dims). Commit.

---
## AFTER #27 → the remaining work
- **043** (record_range path — the one distinct checker path still untested) + optionally
  006/044 (or roll them into the fanout).
- **The 44-stage FANOUT** — user opt-in, via the Workflow tool, per `FANOUT_PLAN.md`
  (Codex as coder, ~5 concurrent, checker frozen, per-stage fidelity agent, export-all +
  `.out` citation by default). Then the ledger-wide composite build + C8 + wiring into
  per-commit.
