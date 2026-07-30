# stage023 step-(f) sealed prediction — and its (h2) adjudication

The prediction below was written 2026-07-28 ~05:15 **outside the repo**, before any `.py` work existed and
with the `.wl` reference frozen at 29 records, per `manifests/DIMENSION_REWRITE.md` §4-(e). It was opened
for the first time at step (h2), after both review legs and both remediation rounds had landed. What can
be stated: it was kept outside the repository until (h2), and none of the five directives or prompts
issued to a leg before it was opened — the two review prompts, the two remediation directives, and the
first prose decision list — names it or any prediction file. ⚠ **`UNDETERMINED`:** whether an agent or
leg could nonetheless have read it. Nothing here establishes that. No filesystem-access record was kept,
and the same limit as the custody caveat below applies — an artifact cannot attest to who did not read
it.

⚠ **Custody caveat, unchanged from the sealed text:** committing it afterwards means git cannot prove it
predates the build. Its value is as a *working* pre-registration — it disciplines expectations and records
the ones that turn out wrong — not as cryptographic evidence.

---

## The sealed text, verbatim

> **P1 — the comparator's headline.** `py=29|wl=29|shared=29|py_only=0|wl_only=0|mismatches=0`, exit 0,
> and `ARTIFACT_NAME_WAIVERS` stays `{}`.
>
> **P2 — the name-mapping step is where a defect will land, if one does.** The `.wl` record names are
> snake_case join keys (`sourced_dims.K_eta`, `Z0_ret`, `Omega_U`, `R_mix`, `g_U`, `c_s`, `q_free`,
> `eta_null`, `P0_physical`), but the `.py`'s **variable identifiers** are CamelCase (`Keta`, `Z0ret`,
> `OmegaU`, `Rmix`, `gU`, `c_s`, `q_free`, `eta_null`) while its **SymPy symbol names** already carry the
> snake_case forms. ⇒ I predict the `.py` emitter keys off the symbol `.name`, not the variable, and that
> if it keys off the variable instead, the comparator fires `py_only`/`wl_only` on roughly nine names at
> once — a loud, unmissable failure rather than a silent one.
>
> **P3 — the operator inversion.** Like stage016, the local machinery carries exponent **vectors**
> (`dim_add`, `dim_scale`) while the shared module carries **quantities** (`*`, `/`, `**`). I predict the
> conversion inverts every operator, and that this is the one place a faithful-looking but wrong refactor
> could pass — so the differential check against the pre-rewrite values is the load-bearing gate, not the
> comparator.
>
> **P4 — what will NOT be caught, and I want this on record.** The five intermediates (`P_port`,
> `Delta_port`, `N0_from_port`, `P0_raw`, `K1`) are absent from **both** artifacts. That symmetric
> omission is invisible to the comparator by construction. It is disclosed in the note, and I predict
> nothing in step (f) surfaces it independently.
>
> **P5 — the fractional pair.** `g_U`/`g_W` must serialise as `-1/2, 1/2, -2` in the sidecar. ⚠ Given the
> measured parser gap, a float rendering would pass `load_dimensions` silently, so I predict the sidecar
> must be checked on its **raw token**, and that the comparator alone would not catch a decimal.
>
> **P6 — the risky one.** I predict the `.py` conversion does **not** change any declared value, and
> therefore that the canonical table gains a third `K_eta` member and a second `T_Omega` member, both
> landing as `NEEDS_ADJUDICATION` rather than `AGREE`. ⛔ If instead they land `AGREE`, something has
> silently normalised a reduction level away and the run must stop.
>
> **P7 — the one I most expect to be wrong.** That step (f) needs only one build round. Every build this
> session needed two or more.

---

## Adjudication

**P1 — CONFIRMED.** Measured repeatedly, most recently by the orchestrator re-running the named command
after the final remediation: `RESULT|stage=stage023|status=PASS|mismatches=0` with
`py=29|wl=29|shared=29|py_only=0|wl_only=0`, exit 0, and `ARTIFACT_NAME_WAIVERS|…|py_only=(none)|wl_only=(none)`.

**P2 — ⛔ FALSIFIED IN MECHANISM, confirmed only in outcome, and the falsification is the more useful
half.** The emitted spellings are the snake_case symbol names, as predicted. But the emitter does **not**
key off `symbol.name`: `dimension_records()` is a hand-typed dict literal in which each join key is a
string constant written beside an independently written source expression
(`"sourced_dims.K_eta": SOURCED_DIMS[Keta]`). Name and value are therefore two independent acts of typing,
which is precisely what makes a silent mis-binding constructible — the class recorded as **U13**, where
five same-class rebindings produced a byte-identical payload and a green comparator (artifact retained:
`notes/stage023_step_h_evidence/u13_u14_rescued/`). The prediction's
premise — that a name defect would be *loud* — is wrong for exactly the defect class that matters: had
the emitter keyed off `symbol.name` as predicted, that class could not exist.

**P3 — ⚖ SPLIT: the mechanism half CONFIRMED, the exclusivity half FALSIFIED.** *Confirmed:* the
conversion did invert every operator (`dim_add`/`dim_scale` → `*`/`**`), and the
transliteration-fidelity leg checked the inversion term-by-term and found it value-equivalent and
strictly stronger than what it replaced — the old `zip` truncated a length mismatch silently, while the
module's `_combine` raises. The site was correctly identified as a risk site, and the differential check
found no defect there. ⛔ *Falsified:* the claim that this is **"the one place"** a faithful-looking but
wrong refactor could pass is false as written, because this same session measured a second such place —
the hand-typed join keys in `dimension_records()`, where five same-class rebindings leave the stage
green, the emitted payload byte-identical and the comparator passing, recorded as **U13**. ⭐ The finding
under the falsification: the exclusivity was never scoped to a corpus or to a mechanism set, which is
the defect class `manifests/DIMENSION_REWRITE.md` §4-(c1)(5) exists to catch — appearing here inside an
adjudication of a prediction about detection blindness.

**P4 — CONFIRMED.** The five intermediates remain absent from both artifacts, the symmetric omission
remains invisible to the comparator by construction, and nothing in step (f) surfaced it. What surfaced
it was the step-(c1) physics leg (§1.7(1), and **U11** for why two of the five matter more than the other
three) — a different instrument entirely, which is the point.

**P5 — CONFIRMED, both halves.** The sidecar carries the exact tokens
`sourced_dims.g_U|exponents={-1/2, 1/2, -2}` and the same for `g_W`; all 87 exponent tokens across the
artifact are signed integers or exact `n/d`. The second half — that the comparator alone would not catch a
decimal — was independently measured by the adversarial leg, which observed float exponents surviving the
comparator once the module pin is re-accepted. The raw-token check was therefore load-bearing, as predicted.

**P6 — CONFIRMED, including the stop condition that did not fire.** No declared value changed. The
canonical table now carries `KEta` with three members (stage013 `M L⁻¹T⁻²` line · stage016 `M L⁻³T⁻²`
volume · stage023 `M T⁻²` reduced) and `TOmega` as a new two-member group (stage016 `M L⁻³T⁻²` ·
stage023 `M T⁻²`), **both `NEEDS_ADJUDICATION`**, taking the table from 2 such groups to 3. Neither landed
`AGREE`, so the ⛔ stop condition did not fire. This is correct output, not drift: they are reduction
levels of one object, each off by exactly its measure's L-power (§7).

**P7 — CONFIRMED for what it claimed, and the meta-expectation attached to it was wrong.** Step (f) — the
`.py` build itself — reached green in one round and needed no rebuild. ⚠ Scope, stated precisely: the
author expected this prediction to fail, and it did not; but two further code rounds *were* needed at
step (h), driven by review findings rather than by the build. So "one build round" held while "one round
total for the stage" did not — and the prediction only ever claimed the former. ⚠ The round count for (f)
is taken from the build session's own record, not re-verified here.

**Score: 5 fully confirmed · 1 falsified (P2) · 1 split (P3: mechanism confirmed, exclusivity
falsified).** Both failures changed how the stage is documented: P2's is why **U13** exists, and P3's
exclusivity half is a claim scoped to no corpus and no mechanism set — the defect class
`manifests/DIMENSION_REWRITE.md` §4-(c1)(5) exists to catch, met here inside this adjudication.
