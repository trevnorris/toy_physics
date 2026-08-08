# Scope — two cross-step dimension checks, and the R4 decision behind them

**Status: SCOPE, 2026-08-07. ⛔ Nothing built. ⛔ Not reviewed.**
Replaces the S9 pilot as the proposed next action. Both review legs on `HARNESS_S9_PILOT_PLAN.md`
independently recommended these instead of that pilot.

⚠ **One claim from those legs did NOT survive checking, and it changes the cost** — see W0.

---

## W0 · ⛔ The R4 repair is a DECISION, not a fix — ⛔ it is not "a few lines of YAML"

Both legs said the S9/S10 R4 bindings could be repaired cheaply. ⛔ **Checked against the outputs: they
cannot.**

**What R4 asks for** (`reduction/relations.yaml`):
```
residual = Q.brane.c_gamma − Sqrt( Q.brane.mu_R / Q.brane.rho_br )
designated_output: Q.brane.c_gamma          declared dimension [1,-1,0]
```

**What the bindings supply** — a **squared** speed, dimension `[2,-2,0]`:
```
checks_S9.yaml :194   Q.brane.c_gamma ← WL_S9_CANDIDATE_SPEED_SQUARED1
checks_S10.yaml       Q.brane.c_gamma ← WL_S10_MAIN_D3_Q3_DISTINCT_ROOTS (select: sequence_second)
```

⛔⛔ **NEITHER ENGINE EMITS A NON-SQUARED SPEED AT ALL.** Every S9 speed tag is a `*_SPEED_SQUARED*`
family member; the SymPy side likewise emits only `*_SPEED_SQUARED_CANDIDATES`. ⇒ **there is no tag to
re-point the binding at**, so this is not a re-point.

### ⇒ Three routes, and the choice is the user's

| route | what changes | cost | ⚠ what it commits us to |
|---|---|---|---|
| **A** engines emit `c_γ` itself | both engines gain one tag | engine edit ×2 + rerun + legs ×2 per engine | ⚠ taking a **square root** is a branch choice; R4's own assumption block already declares `positive`, so the positive root is defensible — ⛔ but it must be **emitted from the computation**, not asserted |
| **B** add a squared-form relation to the registry | `relations.yaml` gains `c_gamma² − mu_R/rho_br` | small | ⛔ two relations for one fact; which is canonical, and does the duplicate ever get checked against the original? |
| **C** re-designate the registry quantity | `quantities.yaml` gains `Q.brane.c_gamma_sq` | small | ⚠ grows the registry; ⭐ but the registry **must** grow anyway (5 relations / 14 quantities against ~25 steps) |

⭐ **Recommendation: A.** It is the only route where the engine *computes* the object the registry names,
and W2 below then checks it. ⛔ B and C make the residual green by changing the question.

⚠ ⛔ **Do not do W0 first.** W2 diagnoses it for free — see below.

---

## W1 · Dimension type-check on registry substitutions

**The defect.** `reduction/engine_output_checks.py:2519` accepts any scalar:
```python
if isinstance(selected, bool) or not isinstance(selected, (int, sp.Basic)) or _is_boolean_value(selected):
    status, detail = "INVALID_SHAPE", f"{qid}:{tag}: substitution must be scalar"
```
⇒ a **squared speed** is substituted where a **speed** is designated, and the only symptom is a mysterious
`NONZERO` residual that has been carried as a note at S9 since the step closed — ⛔ and is **also** red at
S10 with the same category error, which the S10 record does not say.

**The change.** Before substituting, compare the substituted value's dimension against the qid's declared
`exponents` in `quantities.yaml`. On mismatch emit a **distinct** status — `DIMENSION_MISMATCH` — ⛔ never
folded into `NONZERO`.

⭐ **What it buys, stated at its real width:** it converts a check that **cannot fail informatively** into
one that names its own failure. ⛔ It does not make R4 true; W0 does that.

**Able-to-fail.** ⭐ It already has a live positive: today's S9 **and** S10 bindings must both report
`DIMENSION_MISMATCH` the moment it exists. ⛔ If either reports clean, the check is wrong.
⚠ Negative control: a correct binding (after W0-A) must report clean, and perturbing the declared exponent
by one must flip it.

**Cost:** one function, inside `check_registry_residuals`. No engine change, no rerun. Hours.

---

## W2 · ⭐⭐ Compare EMITTED dimension vectors against the REGISTRY's declared vectors

**The gap, verified.** Nothing anywhere compares an engine's emitted dimension vector to
`quantities.yaml`. Confirmed by a review leg and re-checked here: **zero `registry` references** in
`engine_output_checks.py:1897-2410`, the whole dimension section. The comparator compares engine-to-engine
dimensions; the homogeneity gate walks registry relations. ⛔ **Neither closes the loop between them.**

**Why this is the highest-value item in this document.** ⭐⭐ **S9's one measured defect that survived two
review legs AND a full ablation suite was a wrong dimension** (`steps/S9_light_requires_shear.md:187`) —
caught only because a second engine happened to disagree. ⇒ this check catches that class **from one
engine**, against an oracle outside every step's shared specification.

**Both sides already exist.** The engines emit dimension vectors and the registry declares them:
```
engine:    WL_S9_SPEED_SQUARED_IMPLIED_DIMENSION: {{2, -2, 0}}
           PY_S9_MAIN_DIM_SPEED_FROM_EXPRESSION: [(mu_R/rho_br, Matrix([[2],[-2],[0]]))]
registry:  Q.brane.c_gamma  exponents [1,-1,0]
```
⭐ **Note what that already shows: `[2,-2,0]` against a declared `[1,-1,0]`.** ⇒ **W2 diagnoses W0 for
free**, which is why W2 comes first.

### ⛔⛔ Two conditions, and skipping either makes the check worthless

1. ⛔ **The axis convention must be CHECKED, ⛔ never assumed.** `quantities.yaml` declares
   `LTM-exponent-vector-v1`; the engines emit bare triples. ⚠ A wrong axis order makes the comparison
   trivially agree or trivially disagree — and this repository has **already measured** an `M,L,T`
   vs `(L,T,M)` disagreement between an engine pair (`STATUS.md`, stage037). ⇒ the convention tag must be
   carried on **both** sides and compared before any exponent is.
2. ⛔ **The `D`-dependence must be handled explicitly.** The engines emit dimensions **symbolic in the
   brane dimension `D`** — `rhoBr -> {-D, 0, 1}`, `muR -> {2-D, -2, 1}` — while the registry declares the
   `D = 3` specialisation `[-3,0,1]` and `[-1,-2,1]`. ⚠ They agree **at D = 3 and only there**. ⇒ either
   specialise explicitly and **record that the check is D-specialised**, or teach the registry the
   `D`-dependence. ⛔ Do not let a `D`-symbolic vector silently compare equal to a numeric one.

**Able-to-fail.** ⭐ Perturb one declared exponent in `quantities.yaml` and confirm the row flips; ⭐ perturb
one **emitted** exponent and confirm it flips independently. ⛔ Both directions, or it is a demo.
⚠ And run it on the **existing** S9/S10 outputs before changing anything — whatever it finds there is a
real finding about closed steps.

**Cost:** a comparison function plus the two conditions above. No engine change. Days, not hours, because
condition 2 is a modelling question, not plumbing.

---

## Sequencing, and what each needs

| # | order | author | legs |
|---|---|---|---|
| **W2** | ⭐ first — it diagnoses W0 and needs no decision | Codex | fresh Claude + Grok |
| **W1** | second — small, and W2 will have shown what it must catch | Codex | fresh Claude + Grok |
| **W0** | ⛔ **last, and only after the user picks a route** | route A ⇒ one engine at a time | fresh Claude + Grok, **per engine** |

⛔ **No `HARNESS_STANDARD.md`, no adjudication artifact, no decision-list rounds, no pilot.** ⭐ These are
edits to instruments that already exist, and each one either catches a way the physics could be wrong or
it does not ship.

---

## ⚠ What this does NOT do — stated so it is not read in

- ⛔ It does **not** address the naming/negotiation cost. That question is **unresolved**: a leg measured
  that deleting the naming layer moves **4 of 12** S9 rows, and the plan conceded twice that deleting it
  relocates the cost into the symbol→value map. ⇒ ⛔ **the days-per-step saving remains unestablished.**
- ⛔ It does **not** close the shared-spec blindness class in general — ⭐ though a leg measured that a
  **repaired** R4 fires on exactly the exhibited mutant and only on the mutant, so W0-A plus W1 closes
  that *instance*.
- ⛔ It does **not** touch S9's physics debts: P2 uncomputed, the longitudinal absence **assumed rather
  than derived**, and the curl-only conditionality.
- ⚠ Coverage stays thin: **5 relations over 14 quantities in 2 sectors.** ⇒ ⭐ growing the registry is part
  of every future step's cost.
