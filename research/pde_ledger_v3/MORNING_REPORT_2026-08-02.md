# Morning report — where the derivations settled

**Overnight, 2026-08-02.** Two things were owed: the `n_pol` registry decision, and S10.
Both are done. **S10's physics is banked and both review legs came back CLEAN.**

---

## 1. The registry decision — `D_brane` in, `n_pol` out

⭐ **Settled by test, not by argument**, which turned out to matter because the machinery had an opinion.

**What I tried first:** register both `D_brane = 3` and `n_pol = 2`, with a relation `n_pol = D_brane − 1`
to encode *why* the count is what it is.

**What happened:** `constraint_dimension` **raises** —
`constraint uses symbols outside the dimension variables: ['D_brane', 'n_pol']`. That is a deliberate
guard keeping discrete facts out of the continuous parameter count, not a bug.

**And the schema has no vocabulary for it.** `additionalProperties: false`, no field expressing *"this
discrete value derives from that one."* So a bare `n_pol` row would make a **derived** count look like a
**free choice** — worse than omitting it.

⇒ **Decision: register `D_brane` only.** Verified clean — it lands on the discrete axis beside `n_eos`,
the continuous count stays 10, residue stays 6, acceptance payload unchanged.

⭐ **And this closes a hole that was already open.** S9's dimensions *depend* on the brane being
three-dimensional — an ablation showed both `[ρ_br]` and `[μ_R]` change when it changes — and that input
was **declared nowhere in the registry.** S10 registering it retroactively grounds S9.

⚠ **The count lives in the scripts and the step record, not the registry.** If you want it registered,
that needs a schema change, and I did not make one unasked.

---

## 2. S10 — two transverse photons

### The result

| `D` | `ω² = 0` | `ω² = (μ_R/ρ_br)k²` | sum |
|---|---|---|---|
| 2 | nullity 1, ∥ `k` | nullity 1, ⊥ `k` | 2 |
| **3** | nullity 1, ∥ `k` | **nullity 2, ⊥ `k`** | 3 |
| 4 | nullity 1, ∥ `k` | nullity 3, ⊥ `k` | 4 |
| 5 | nullity 1, ∥ `k` | nullity 4, ⊥ `k` | 5 |

⭐ **The count is COMPUTED** — from the nullity of the dynamical matrix at each root, in both engines.
⛔ **Not** read off a multiplicity in a factorised determinant, which is how I had been getting it by
hand. Both reviewers independently confirmed this with `file:line`, and grepped for a literal count in
each source: none.

### ⭐⭐ The claim this settles, in a stronger form than I first put it

I had said the count is `3 − 1`, not the codimension `4 − 3`. The real statement is better:

> ⛔ **The bulk never enters the computation at all.**

There is no codimension anywhere in the calculation — the nullity depends only on the brane's own
dimension. **Codimension is not merely the wrong number; it is an absent quantity.** The `D`-table is the
evidence, computed case by case with nothing extrapolated.

⇒ ⭐ **Light having exactly two polarisations is a statement that our space is three-dimensional.**

### Why the longitudinal is free

A review leg contributed the clean reason, and it is better than the case-split I walked:

```
(1/2) Σ_{i,j} (k_i a_j − k_j a_i)²  =  |k|²|a|² − (k·a)²
```

⇒ the stiffness form **penalises only the perpendicular part**. The longitudinal mode is not forbidden —
it costs **zero energy**, so there is no restoring force. That is exactly what *curl-only* means: the
medium resists twisting and is indifferent to compression along the wave.

The two transverse modes are **degenerate by symmetry**: nothing in the equation refers to *which*
perpendicular direction, so both carry the same `ω`.

### ⭐ A loose end from S9 closed itself

The SymPy side derived the dimensions as a closed function of `D`:

```
[ρ_br] = (−D,  0, 1)          [μ_R] = (2−D, −2, 1)
```

Recall the S9 ablation: changing the assumed brane dimension corrupted both constituent dimensions while
the audit's dimension check stayed green. I filed that as a blind spot.

⇒ ⭐⭐ **It is not a blind spot, it is an identity.** `[μ_R] − [ρ_br] = (2−D) − (−D) = 2` in the length
slot, **for every `D`**. The speed's dimension *cannot* see the brane dimension. What looked like a weak
check is a property of the physics.

---

## 3. How it was verified — and why the agreement is not vacuous

**Two engines, written independently.** The Mathematica audit is **blind**: no arguments, no imports, no
file reads. It was built **first**. The SymPy audit — which *does* import the shared registry — was built
while the Mathematica file was **moved out of the tree**, so it could not transcribe from it.
✅ **The quarantine held**: no git-recovery attempt in the builder's log, and the restored file is
byte-identical to its committed blob.

⭐ **The routes genuinely differ**, which is what makes the agreement mean something: the `.wl` takes the
Hessian of the plane-wave quadratic action with a **fully symbolic** wavevector; the `.py` derives its
matrix from **Euler–Lagrange** with an explicit numeric unit direction. A reviewer noted the `.py`'s
single-direction specialisation was the one place a direction-dependent degeneracy could hide — and the
`.wl`'s symbolic `k` closes exactly that.

**Predictions were committed before either script ran** (`bc276485`). Every value matched, including both
uncertainties I had flagged: the antisymmetric form **does** reduce to `|∇×u|²` at `D=3` (computed,
residual 0), and `D=2` is **not** degenerate in kind.

**The form control did the work a coefficient control cannot.** Replacing the antisymmetric-derivative
stiffness with full gradient-squared collapses the two roots into **one**, nullity `D`, and the `ω²=0`
root disappears — the longitudinal propagates too. Rescaling `μ_R` instead moves the roots and leaves the
nullities at 1 and 2, untouched. ⇒ **A coefficient control tests the arithmetic; only a form control
tests a claim about a count.** That is the S9 review's lesson applied rather than restated.

**Both review legs: CLEAN.** No blocking findings from either. Five nits between them, none passing the
physics filter; the most substantive is that the `.wl` computes its equations of motion and its dynamical
matrix by two independent routes and never compares them — a reviewer verified by hand that they agree,
so nothing is currently wrong, and cross-engine agreement already covers the physics. **Recorded, not
acted on.**

---

## 4. Two things I got wrong, and one process fix

⚠ **My hand algebra was unverified and you were right to reject it.** Of what I wrote during the walk,
the equation of motion and the dynamical matrix came from the S9 script; the decomposition, the case
split, the count and the degeneracy were **mine and unexecuted**. The scripts have now confirmed all of
it — but it was not evidence when I presented it.

⚠ **I created a leak channel while trying to be rigorous.** I committed my pre-registered predictions to
a **tracked** path so they would prove they pre-dated the results — which put my predicted table where a
reviewer could read it. Grok shows clear independent work (it derived the stiffness identity, which is
nowhere in my file), but I cannot prove it was unanchored. ⇒ **Fix: commit the pre-registration, then
move it out of the tree for the review**, exactly as I did with the `.wl`. I had the technique and did
not apply it to my own artifact.

⛔⛔ **A SILENT NO-OP THAT PASSED MY OWN VERIFICATION RULE — the most useful failure of the night.**
The first attempt to build S10's TeX card **did nothing**, and reported success. My shell's working
directory had persisted from an earlier `cd`, so the heredoc wrote the directive to a **relative** path
that did not resolve, while `codex exec` read it back by **absolute** path — getting an empty string.
Codex started, was handed nothing, asked *"What would you like me to work on?"*, and exited 0.

⚠ **It passed the rule I have on record.** *"Verify a background job by its ARTIFACT — `hook: Stop`,
`tokens used` — never by a process check."* Both markers were present. But they prove **the session
finished**, ⛔ not that **the deliverable exists**. A completed session that did nothing satisfies the
check exactly.

⇒ ⭐ **The rule needs a second clause: verify the DELIVERABLE, not the session.** The tell was in plain
sight and I nearly missed it — **3,130 tokens** against 36,979+ for the real builds. A build that costs
almost nothing did almost nothing.
⇒ ⭐ **And: use absolute paths for anything a background job will read.** The shell cwd persists between
calls; a relative write followed by an absolute read is a silent mismatch.

⚠ **A directive-wording artifact worth not repeating.** I wrote *"an integer that appears anywhere in the
source as the answer is a failed build"*, and the builder over-generalised into avoiding integer literals
at all — writing the control's scale factor `2` as `Length[ordinaryCurlVector] − Boole[...]`. Harmless
(it is a control parameter, never a reported count; I checked). But it is the denylist pattern: **I
banned a spelling when I meant to ban a behaviour.**

---

## 5. Where the ledger stands

**Half one, Part I (Light):** S9 closed · **S10 closed**.

| commit | what |
|---|---|
| `bc276485` | S10 predictions, pre-registered |
| `a942c950` | blind Mathematica audit — the count, computed |
| `90c30ed8` | SymPy audit + the `D_brane` registry row |
| `2769362f` | S10 step record |

All registry gates green, payload unchanged: acceptance `MATCH` (10→6, 10→6, 10→5), dimensional
homogeneity `PASS`, able-to-fail `PASS`, 11 tests.

**In flight when I stopped:** Codex building S10's TeX card (per your rule that it builds the `.tex`).

### ▶ Next is S11 — and it is the interesting one

The stray longitudinal. What S10 leaves: `ω² = 0` means the longitudinal is **non-propagating**, ⛔ **not
absent**. Once the medium is allowed to be compressible, that zero lifts and the slot becomes a real
propagating mode.

⛔ **And per your correction, that is not a defect awaiting repair.** Exact Maxwell would be the
**failure** — it puts charge in by hand, so a model matching it exactly would have no way to physically
anchor charge. The extra mode **is** the anchor, and it is what made the drum-head picture click. That
framing is now recorded in `CHARTER.md#falsification-standard`, `V3_STEP_PLAN.md` § S11, and the charge
memory (which now leads with your drum-head description instead of burying it under build verdicts —
the reason it kept not sticking).

### Still open, unchanged

- **`ħ`'s class** — `postulated` or `calibrated`. The one decision genuinely owed by you.
- **O-02** — is `K` + the EOS exponent one entry or two.
- The requirements-first step **ordering** in the plan (substrate still listed first).
