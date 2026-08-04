# S9 TeX card — repair: an over-read tool result, and a missing headline

⛔ **Edit in place. Four fixes, nothing else.** ⛔ Do not commit.

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S9_light_requires_shear.tex`
**Source of truth:** `/var/projects/toy_physics/research/pde_ledger_v3/steps/S9_light_requires_shear.md`

⚠ Verification content must stay in `\paragraph{Verification.}`. ⛔ **Never** move it into
`\stagefield{Verification}` — `paper/macros.tex:18-26` suppresses that field in the default build.

---

## T1 — ⛔⛔ BLOCKING: the card presents a FAILING, self-disclaiming triage run as a clean result

The card quotes the automated consumer as *"12 cross-engine quantities agreeing and 0 disagreeing, with
all 1219/1219 dimensionful expressions homogeneous."*

⛔ **It keeps the two flattering numbers and drops everything qualifying.** What the tool actually prints:

```
OPERATIONAL_FAILURE: 1 value(s) were UNPARSED
CONTROL_RESPONSE: compared=170 responsive=150 invariant=20 unparsed=1 unpaired=8
DIMENSIONS: total=1559 checked=1219 homogeneous=1219 non_homogeneous=0
CROSS_ENGINE: agree=12
TAG_PARITY: packages=16 gaps=16
LIMITATION: triage only; this run does not establish physical correctness, completeness, or
            derivation coverage.
```
**and the run exits non-zero (2).**

⚠⚠ **This is the failure mode this ledger has had before: a tool's output quoted as coverage it never
had**, in the one artifact a reader actually trusts, under a heading reading *"Verification"*.

⭐ **Fix — report it honestly, in the same paragraph:**
- add **`150 of 170` compared tags respond to some control**, i.e. ⭐ **20 respond to none** — that is the
  only figure measuring able-to-fail *reach*, and dropping it while keeping the pass numbers is the
  over-read;
- add that **tag-set parity gaps are reported** (they are legitimate, but they are reported);
- state plainly that the **run terminates non-zero** on the unparsed `Piecewise` tag;
- carry the tool's **own** disclaimer: it is **triage**, and ⛔ it does not establish physical
  correctness, completeness, or derivation coverage.

⛔ Do not soften any of these and ⛔ do not move them to a different paragraph from the pass numbers.

## T2 — ⛔ the card never states the MAIN action's computed result

⚠ The card's only quantitative mode count and speed belong to a **control** (gradient elasticity).
⛔ **A reader of this card alone cannot tell the curl-only action was solved at all**, let alone what it
gave. The conditional framing appears **twice** while the positive result appears **nowhere**.

⭐ **Fix — state the main action's computed result plainly, before the conditional:**
```
det M = ω²ρ_br (ω²ρ_br − μ_R k²)²
roots {0, μ_R k²/ρ_br}
at ω² = μ_R k²/ρ_br : two transverse modes (E2 = 2), not longitudinal
at ω² = 0           : E2 = 0, longitudinal
c² = μ_R/ρ_br, wavevector-independent (scaling residual 0)
E2 drops to 1 when the inertia is made anisotropic
```
⭐ **Then** the conditional. ⚠ The step exists to establish this; ⛔ it must not be inferable only from a
control.

## T3 — attribute the one-sided-corruption claim to where it happened

The card states that one-sided corruption proved the two matrix routes independent. ⛔ **That ablation is
not reproducible from the engines' output** — it was performed by review legs on modified copies.
⭐ Attribute it as an **independent review result**, ⛔ not as something a reader can rerun from the
committed scripts.

## T4 — the register citation does not resolve

`\StageFile{notes/parameter\_register.md}` has no such path: v3 has no `notes/` tree and the file is
**v2's** (`research/pde_ledger_v2/notes/parameter_register.md`, lines 137–138; `R4` at line 271). Every
other `\StageFile` in the card is a full repo-relative path. ⭐ Fix the path and ⭐ say it is **v2's**
register.

---

## Report back — under 15 lines

1. The file changed; whether `pdflatex` still completes cleanly.
2. One line per fix T1–T4, quoting the sentence each now reads as.
3. ⭐ Anything else in the card you believe over-reads a tool result or states a conclusion the engines do
   not support. ⭐ This is wanted.

⚠ Note for the record, ⛔ not something to fix here: the committed `pde_ledger_v3.pdf` is **stale** — it
still contains pre-rebuild text. It needs a rebuild before anyone reads the PDF as current.
