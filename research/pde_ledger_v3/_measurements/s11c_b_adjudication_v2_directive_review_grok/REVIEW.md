# Independent review — S11c-b ADJUDICATION v2 build directive (decision list)

**Artifact:** `directives/S11c_b_adjudication_v2_build_directive.md`
**Verdict:** **NOT SOUND for build launch** — four blocking findings. Do not launch the builder until Bridge D is rewritten to the engines' `PROFILE_GRADE` map (not naive chain-through-`eta_bg`), `mu_R_bg` is in scope or explicitly deferred, divergence fixtures are bound to the same `DEPENDENT_BASES` Euler path used in production, and the preamble's expected residual taxonomy / false ADMISSIBILITY-agree claim are cut.

**Method:** engine sources + committed adjudication/comparator read first; Bridge D chain rule computed in CAS (`bridge_d_chain_rule_vs_profile_grade.py` + `.stdout`); baseline route headers from `/tmp/S11c_b_adjudicated_run.out`. No Mathematica spawn.

---

## What was verified (source-grounded)

### Bridge D zero-jet (correct as far as it goes)

PY `profile_definitions` / `PROFILE_GRADE_SUBS`:

```649:663:research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
        sp.Eq(W_bg, W0 * (1 + eta_bg * w1_profile), evaluate=False),
...
PROFILE_GRADE_SUBS = {
    W_bg: W0 * (1 + eta_bg * w1_profile),
```

WL keeps thickness as `anchoredWidth = backgroundAnchor[...]` (WL:448,535) with zero-jet value `widthValue = WZero (1 + etaBg w1ProfileZero)` (WL:294). Zero-jet substitution `W_bg ↦ W_0·(1+eta_bg·w1_profile)` matches both engines.

### Bridge D jets (directive disagrees with both engines)

PY first-jet definitional map is **`sigma_W · w1_profile_d{i}`**, not `W_0·eta_bg·w1_profile_d{i}`:

```652:665:research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
        sp.Tuple(*(sp.Eq(grad_W[i], sigma_W * w1_grad[i], evaluate=False) for i in DIRECTIONS)),
...
    **{grad_W[i]: sigma_W * w1_grad[i] for i in DIRECTIONS},
```

PY second jet of `W_bg` (operator / background dx paths):

```1864:1866:research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
            derivative_map[grad_W[jet_direction]] = (
                sigma_W * w_profile_second / L_W
            )
```

WL `profileRulesRetainingGeneratedJets` (WL:347–350): `Derivative[…][widthBase] → sigmaW * profileJetSymbol["W_BG", orders] / LWidth^(i+j+k-1)`.

CAS residual (literal stdout):

```
FIRST_JET_RESIDUAL_directive_minus_engine = W_0*eta_bg*w1_profile_d1 - sigma_W*w1_profile_d1
SECOND_JET_RESIDUAL_directive_minus_engine = W_0*eta_bg*w1_profile_d1d1 - sigma_W*w1_profile_d1d1/L_W
```

`eta_bg` and `sigma_W` are independent grade bookkeepers (PY:156–157; prior repair text forbids collapsing them). Differentiating the zero-jet formula through `eta_bg` is continuum chain rule on a different object than `PROFILE_GRADE_SUBS`.

### Baseline FLAG shape (why the wrong jet map matters)

From `/tmp/S11c_b_adjudicated_run.out` (also `baseline_route_headers.txt`):

| route | count |
|---|---|
| ALGEBRAIC_MATCH | 38 |
| ALGEBRAIC_FLAG | 36 |
| CONTAINER_FLAG | 4 |
| STRUCTURE_INCOMPLETE | 69 (includes ENERGY_BASIS_VARIABLE/NEW/OMISSIONS) |

Sample A/B background atoms:

- **SLAB MASS:** A=`W_bg_d*`; B=`eta_bg, sigma_W, w1_profile*, m1_profile_d*`
- **MU_THETA:** A=`W_bg, W_bg_d*, mu_R_bg_d*`; B=`eta_bg, sigma_W, w1_*, m1_*`
- **ADMISSIBILITY BODY_FORCE/THETA:** A=already profile-expanded nonzero (`eta_bg, sigma_W, L_W, w1_profile, w1_profile_d{i}d{i}`); **B=`Integer(0)`** — not a W_bg spelling bridge

### Comparator divergence machinery (usable, with caveats)

`modulo_total_divergence` (`S11c_b_cross_engine_comparator.py:555–582`) returns the spatial Euler signature on `DEPENDENT_BASES` only; it annihilates total in-plane divergences for those fields and does **not** reconstruct `V`. If no dependent-base symbols appear, it returns the expression unchanged (does not blanket-zero). Preferring this path matches the directive's "reuse comparator" preference; fixtures must hit this path.

### Protected ENERGY_BASIS

Directive exclusion of divergence classification on ENERGY_BASIS is stated and matches v1's quotient warning. Baseline keeps ENERGY_BASIS_VARIABLE/NEW/OMISSIONS as `STRUCTURE_INCOMPLETE`; `07/10` and gamma-DivGrad live inside that quotient content. Exclusion as written is complete **if** the builder gates by family name and never runs the classifier on STRUCTURE_INCOMPLETE ENERGY rows.

---

## Findings

### B1 — BLOCKING — Bridge D jet map is not the engines' definitional expansion

**Where:** directive L31–32.

**Fact:** Directive mandates
`W_bg_d{i} ↦ W_0·eta_bg·w1_profile_d{i}`,
`W_bg_d{i}d{j} ↦ W_0·eta_bg·w1_profile_d{i}d{j}`.
Engines define
`W_bg_d{i} ↦ sigma_W·w1_profile_d{i}` (PY:652,665; WL scale `sigmaW` at 324–325),
`∂∂ W_bg ↦ sigma_W·w1_profile_d{i}d{j}/L_W` (PY:1864–1866; WL:347–350).
Measured residual nonzero (see `.stdout`).

**Why it changes the build:** A builder implementing the written chain rule builds the wrong bridge. On SLAB/MU/COUPLING FLAGs (A carries `W_bg_d*`, B carries `sigma_W·w1_d*`), the wrong map leaves `(W_0·eta_bg − sigma_W)·w1_d*` instead of cancelling — and/or invents a bookkeeper collapse §2a treats as illegal. DoD's "jet-consistent definitional expansion" would green-light the wrong map.

**Correction:** Replace the chain-rule bullet with the source object `PROFILE_GRADE_SUBS` / WL `profileRulesRetainingGeneratedJets`, explicitly:

- `W_bg ↦ W_0·(1+eta_bg·w1_profile)` (unchanged)
- `W_bg_d{i} ↦ sigma_W·w1_profile_d{i}`
- `W_bg_d{i}d{j} ↦ sigma_W·w1_profile_d{i}d{j}/L_W` (and higher orders per WL `LWidth^(n-1)` if those tokens appear)
- Cite PY:649,652,663–665,1864–1866 and WL:294,347–350,448 — **not** "chain rule on eta_bg"
- State that this retains `w1_profile` jets (still true) and does **not** identify `sigma_W` with `W_0·eta_bg`

### B2 — BLOCKING — Bridge D omits the same-class `mu_R_bg` definitional map present in FLAG residuals

**Where:** directive L25–37 (W_bg only).

**Fact:** `PROFILE_GRADE_SUBS` also maps
`mu_R_bg ↦ mu_R·(1+eta_bg·m1_profile)`,
`mu_R_bg_d{i} ↦ mu_R·sigma_W·m1_profile_d{i}/W_0` (PY:650,655–656,664,666).
Baseline MU_THETA / COUPLING / SLAB B-sides carry `m1_profile_d*`; A-sides carry `mu_R_bg_d*`.

**Why it changes claims:** Leaving `mu_R_bg` unbridged while calling Bridge D "the" missed representational bridge causes same-class spelling residuals to survive as `RESIDUAL_BULK` and be claimable as genuine findings — or tempts an unstated second bridge.

**Correction:** Either (preferred) include the full thickness+modulus profile-grade map under Bridge D with citations, or explicitly defer `mu_R_bg`/`m1_profile` mismatches as out-of-scope for v2 and forbid treating them as divergence-resolved agreement.

### B3 — BLOCKING — Divergence fixtures are not bound to the production Euler/`DEPENDENT_BASES` path

**Where:** directive L49–55.

**Fact:** Committed `modulo_total_divergence` Euler-reduces only symbols in `DEPENDENT_BASES` (`S11c_b_cross_engine_comparator.py:494–507,555–582`). A placeholder `φ` / `f,g,h` **outside** that set never enters the Euler loop (empty `symbols_by_base` → return expression unchanged). Then:

- bulk fixture can "pass" as `RESIDUAL_BULK` without testing the Euler classifier;
- divergence fixture fails unless the builder adds a **second** recognition path (explicit `∂_i` pattern), which production residuals (product-rule expanded) may not share.

**Why it changes the build / false-agreement risk:** Fixtures can be satisfied by a path that is not the path that classifies operand residuals — the exact hole that lets a bulk term be folded on real cases while fixtures stay green. Reconstructed `V` is optional in the text ("V **or** reduced-to-zero witness"); Euler-zero is fine **if** it is the same function used on cases.

**Correction:**

1. Mandate reuse of `C.modulo_total_divergence` (or an equivalent Euler signature on the same `DEPENDENT_BASES`).
2. Fixtures must use fields from `DEPENDENT_BASES` (e.g. `θ_probe` / `e_W_probe` / `u_T_*`), e.g. bulk `|∇θ_probe|² → RESIDUAL_BULK`, divergence `Σ_i ∂_i(θ_probe·e_W_probe·δ_{1i})` (or similar) → `REPRESENTATIONAL_DIVERGENCE` with printed Euler-zero witness.
3. Require production and fixtures to call one function; `--drop-divergence` bypasses that function.
4. Keep ENERGY_BASIS family gate as a hard skip before that function.

(Definition of "total in-plane divergence" via IBP/Euler in L41–43 is otherwise adequate and is **not** blanket `simplify→0`.)

### B4 — BLOCKING — Preamble leaks expected residual taxonomy and misstates ADMISSIBILITY

**Where:** directive L4–9.

**Fact:** Preamble asserts FLAG residuals "decompose into" (missed W_bg bridge + protected quotient + total-divergence coupling structure), and that engines "AGREE on ADMISSIBILITY operator/support/count (38)". Baseline has **4** `ALGEBRAIC FLAG ADMISSIBILITY_OPERATOR_OPERAND` BODY_FORCE/THETA cases with WL operand `0` and PY nonzero profile-expanded content; `38` is the global `ALGEBRAIC_MATCH` count, not an ADMISSIBILITY-only agree.

**Why it changes the build:** Rule 5 — the builder must not be handed expected classifications. The ADMISSIBILITY THETA residual is exactly the false-agreement danger surface for the divergence step (PY density vs WL 0). Telling the builder that FLAGs are "representational / divergence" is an acceptance leak.

**Correction:** Cut the decomposition forecast and the ADMISSIBILITY-agree claim from the directive. Keep only: apply enumerated bridges; classify surviving residuals with the ablatable divergence predicate; print routes; protected sets untouched. Interpretation stays in the step record.

---

## Non-blocking notes

### N1 — Witness `V` vs Euler signature

Comparator has no `V` reconstructor. "Reduced-to-zero witness" via printed Euler Association is enough if B3's single-path rule holds. Optional: drop "reconstructed V" or mark it as non-required.

### N2 — Ablation / DoD adequacy (modulo B1–B3)

`--drop-bridge-d`, `--drop-divergence`, plus v1 hooks, are the right shape and can fail a jet-lowering or blanket-simplify implementation **once** Bridge D and fixtures are corrected. DoD's "resurfaces the W_bg residual" is an ablation load-bearing check, not a per-case expected value — acceptable after B1 fixes the map.

### N3 — Protected sets wording

L57–61 are sound. Positive family allow-list for the divergence step (operator/coupling algebraic FLAG survivors only) would make the ENERGY_BASIS trap harder to miss.

---

## Check-by-check summary

| # | Check | Result |
|---|---|---|
| 1 | Bridge D + jet-consistency | **Fail** — zero-jet OK; jets wrong (B1); mu omitted (B2); not a jet collapse, but wrong expansion |
| 2 | IBP false-agreement | **Partial** — Euler definition OK; fixtures not path-bound (B3); `--drop-divergence` OK; witness OK if Euler |
| 3 | ENERGY_BASIS / 07/10 / gamma-DivGrad | **Pass** as written, if family-gated |
| 4 | Leak / rule 5 | **Fail** — preamble taxonomy + false ADMISSIBILITY agree (B4); fixtures use placeholders (OK); no numeric residual targets in DoD |
| 5 | Ablation + DoD | **Pass shape**, blocked until B1–B3 so DoD cannot green-light a wrong Bridge D |

**Launch recommendation:** revise directive for B1–B4; re-run decision-list legs; then build.
