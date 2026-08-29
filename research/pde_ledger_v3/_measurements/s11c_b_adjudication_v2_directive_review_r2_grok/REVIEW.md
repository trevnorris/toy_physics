# Independent review — S11c-b ADJUDICATION v2 build directive (decision list, round 2)

**Artifact:** `research/pde_ledger_v3/directives/S11c_b_adjudication_v2_build_directive.md`
**Verdict:** **NOT SOUND for build launch** — one blocking ambiguity on protected-atom routing vs Bridge D, plus one blocking underspec on strong-row `HeldDiv` expansion. Bridge D map, strong/weak whitelist, V-mandatory divergence classifier, fixtures, and rule-5 leak fixes from round 1 are otherwise source-aligned.

**Method:** engine sources + comparator divergence machinery read first; Bridge D checked line-by-line against `PROFILE_GRADE_SUBS`; strong/weak split against `variationalSource` / COUPLING weak-pairing construction; baseline FLAG atoms from `/tmp/S11c_b_adjudicated_run.out`; CAS measurements in `bridge_d_and_protected_atoms.py` + `.stdout`. No Mathematica spawn. Did not read the parallel Codex round-2 leg.

---

## What was read (source-grounded)

### Bridge D object = `PROFILE_GRADE_SUBS` (PY:662–671)

```
W_bg      → W0*(1 + eta_bg*w1_profile)
mu_R_bg   → mu_R*(1 + eta_bg*m1_profile)
W_bg_d{i} → sigma_W*w1_profile_d{i}
mu_R_bg_d{i} → mu_R*sigma_W*m1_profile_d{i}/W0
rho_4D_bg_rho4_constant     → rho_br/W0
rho_br_bg_rho4_constant     → rho_br*(1 + eta_bg*w1_profile)
rho_4D_bg_rhobr_constant    → rho_br/(W0*(1 + eta_bg*w1_profile))
rho_br_bg_rhobr_constant    → rho_br
```

`sigma_W` and `eta_bg` are independent symbols (PY:156–157). Measured: `SIGMA_W_IDENTICAL_W0_ETA = False`.

Second jets are **not** dict entries. They appear only in `operator_dx` (PY:1864–1869):
`∂∂W_bg → sigma_W·w1_second/L_W`, `∂∂mu → mu_R·sigma_W·m1_second/(W0·L_W)`.
Measured: `SECOND_JET_W_bg_d1d1_IN_SUBS = False`. Baseline FLAG residuals carry `W_bg_d{i}` / `mu_R_bg_d{i}` but **no** free `W_bg_d{i}d{j}` atoms — second-jet content is already on the profile side as `w1_profile_d{i}d{j}`.

WL keeps anchored symbols (`anchoredWidth`, WL:448,535); profile jet scale is `sigmaW` / `LWidth^(n-1)` (WL:347–355). Directive’s first-jet / second-jet formulae match the engines; the load-bearing instruction “import/reference `PROFILE_GRADE_SUBS`; do not re-derive” is the correct object.

### Strong vs weak

- Strong: `variationalSource = ∂L/∂q − Inactive[Div](...)` (WL:274–277), used for slab/mu construction (WL:793–795). Comparator extract for `SLAB_OPERATOR` / `MU_THETA_OPERATOR` does **not** set `reduce_divergence=True` (comparator:760–845).
- Weak density: `extractCoupling` builds `weakPairingRecord` / `PAIRING_DENSITY_MODULO_COMPACT_SUPPORT_IBP` (WL:1075–1124). Sector density cases are extracted **without** `reduce_divergence`; only ADJOINTNESS_* objects set it (comparator:866–921).
- Baseline FLAG HeldDiv: all 16 `SLAB_OPERATOR` + all 4 `MU_THETA_OPERATOR` residuals carry `HeldDiv` (WL side); 4/8 `COUPLING_KERNEL` too. Drop helper `_drop_held_divergences` replaces `HeldDiv` with `0` (comparator:547–552) — measured on a toy.

### Divergence classifier / fixtures

Directive requires printed raw `R`, explicit 3-component `V`, and `R − Σ∂_i V_i == 0` — Euler may screen only. Measured with fixture-local formal jet derivative: `a_d1·φ + a·φ_d1 − ∂_1(a·φ) = 0`; bulk `a·φ_d1 − ∂_1(a·φ) = −a_d1·φ`. Comparator Euler on placeholders returns the expression unchanged (`EULER_PLACEHOLDER_BULK = a*phi_d1`) — so V-verification (not bare Euler on `DEPENDENT_BASES`) is the right fixture path; fixture-local registry is necessary. `--drop-divergence` is a real ablation (bypass classifier + witness).

### Protected atoms in baseline FLAGs (not ENERGY_BASIS-only)

| family | HAS_PROTECTED | NO_PROTECTED |
|---|---|---|
| SLAB_OPERATOR | 16 | 0 |
| MU_THETA_OPERATOR | 4 | 0 |
| COUPLING_KERNEL | 4 | 4 |
| ADMISSIBILITY_OPERATOR_OPERAND | 0 | 4 |
| SLAB_OPERATOR_TERM_ORIGINS (alg+container) | 0 | 8 |

Protected atoms present: `gamma_s11cb_{w_bg,mu_r_bg}_{07,10}`, `gamma{Width,Modulus}DivGrad{Theta,Ew}`.
Toy: Bridge D remaps `W_bg_d1` and leaves `gammaWidthDivGradTheta` coeff unchanged (`PROTECTED_COEFF_UNCHANGED = True`; no protected key in `PROFILE_GRADE_SUBS`).

### Leak / DoD

Preamble states measured v1 counts only (38 MATCH / 40 FLAG); no per-case expected classification; no ADMISSIBILITY-agree claim. Fixtures use placeholders. Case-ID multiset + drop ablations value-free. Accounting DoD can fail.

---

## Findings

### B1 — BLOCKING — Protected-atom route vs Bridge D is ambiguous (whole-residual freeze vs term hold)

**Where:** directive L61–66.

**Fact:** Text says any residual *containing* a protected atom “is kept RAW under … `PROTECTED_UNREDUCED`” and parenthetically “no Bridge D fold **of that term**, no divergence step.”

Two incompatible builds:

1. **Whole-residual:** presence of any protected atom → skip Bridge D and divergence → emit `PROTECTED_UNREDUCED`. Then **all 16 SLAB + all 4 MU_THETA** FLAG cases never receive Bridge D, even though their residuals also carry the `W_bg`/`mu_R_bg` spelling Bridge D exists to clear. Measured coexistence: every SLAB/MU FLAG has both protected gammas **and** profile-grade atoms.
2. **Term-hold:** apply Bridge D to the residual (safe: `PROFILE_GRADE_SUBS` does not touch protected symbols); if protected atoms remain, forbid divergence and account as `PROTECTED_UNREDUCED` / keep those monomials raw; only protection-free COUPLING residuals enter the V-classifier.

**Why it changes the build / claims:** Reading (1) leaves representational spelling as claimable differences on the main FLAG surface and largely idles Bridge D. Reading (2) is the one consistent with round-1’s atom-gate intent and with measured Bridge D safety. Leaving both open, a builder can implement either; DoD does not discriminate.

**False-agreement surface:** allowing the **divergence** step on a residual that still contains `07/10` or DivGrad gammas (COUPLING THICKNESS_TO_TRANSVERSE has them) can stamp `REPRESENTATIONAL_DIVERGENCE` on a residual that still carries one-engine quotient physics — that *is* manufactured agreement. Bridge D alone does not create that failure.

**Correction:** Mandate an explicit order, e.g.

1. renames + Bridge A + Bridge D (`PROFILE_GRADE_SUBS` import);
2. if residual free symbols intersect the protected atom set → `PROTECTED_UNREDUCED` (print raw; **no** divergence classifier); Bridge D already applied;
3. else if family is `COUPLING_KERNEL` density → V-classifier;
4. else strong exact → `MATCH`/`FLAG`.

State that ENERGY_BASIS_* never enters Bridge D or the classifier (family gate, already present).

### B2 — BLOCKING — Strong-row `HeldDiv` “EXPANDED” is not an executable replacement rule

**Where:** directive L43–44 (“A HeldDiv in a strong row is EXPANDED and compared, never discarded”).

**Fact:** WL slab/mu operands emit `Inactive[Div]` → comparator `HeldDiv`; PY emits expanded derivative content. Baseline: 16/16 SLAB and 4/4 MU FLAG residuals contain `HeldDiv` on the B side only. The only named HeldDiv helper in-tree is `_drop_held_divergences` → **Zero** (comparator:547–552; measured `AFTER_DROP_HELDDIV = bulk`).

**Why it changes the build:** “Expanded” without a replacement invites reuse of the drop helper or `reduce_divergence=True` on strong extract paths — the exact false-agreement mechanism round 1 removed from the density classifier’s scope. Conversely, leaving `HeldDiv` opaque forever makes strong exact compare a permanent form mismatch even after Bridge D.

**Correction:** Define expand as a printed, ablatable map on strong rows only, e.g.  
`HeldDiv(V) ↦ Σ_i formal_∂_i(V_i)` in the declared jet algebra (same formal derivative used for V-checks), **never** `_drop_held_divergences` / Euler quotient. Require that strong families keep `reduce_divergence=False`. Optionally add a fixture: HeldDiv of a known vector in a strong family expands to the formal divergence and is compared; drop ablation must not be used on that path.

---

## Non-blocking notes

### N1 — Second jets listed in Bridge D prose but absent from `PROFILE_GRADE_SUBS`

Directive L25–28 and DoD L86 mention second jets via `sigma_W/L_W`. The committed dict has only zero- and first-jets + four density bookkeepers (12 entries). Second jets live in `operator_dx` (PY:1864–1869). For residual bridging this is harmless (no free `W_bg_d{i}d{j}` in baseline FLAGs). Keep “Bridge D **is** the imported `PROFILE_GRADE_SUBS` object”; treat second-jet lines as engine documentation, not extra hand-written subs. Build-leg “matches `PROFILE_GRADE_SUBS` exactly” already says the right check.

### N2 — Density bookkeeper parenthetical names two of four

L29 cites `rho4_rho4, rhobr_rho4`; the dict also has `rho4_rhobr, rhobr_rhobr` (actual symbol names `rho_*_bg_*_constant`). Importing the dict covers all four. Spell “all entries of `PROFILE_GRADE_SUBS`” to avoid a partial hand copy.

### N3 — Divergence whitelist is the right positive allow-list

Only `COUPLING_KERNEL` is divergence-eligible; strong list matches variational operators. `COUPLING_KERNEL_TERM_ORIGINS` is not a formed weak density (term-origin containers; baseline STRUCTURE / container). `ADMISSIBILITY_SUPPORT_OPERAND` omitted from the strong bullet is fine under a positive divergence allow-list. No other CORE family is a `weakPairingRecord` density in the WL construction cited.

### N4 — Fixtures + `--drop-divergence`

Product-rule pair + strong-ineligible route + one shared function + fixture-local registry exclude the round-1 Euler/`DEPENDENT_BASES` fixture hole. V-mandatory verdict blocks Euler-zero-as-agreement. Sound as written once B1 forbids divergence on protected-bearing residuals.

### N5 — Rule 5

Round-1 preamble taxonomy / false ADMISSIBILITY agree are gone. Remaining purpose sentence (“representational identity or genuine difference”) is route vocabulary, not a per-case forecast. Acceptable.

---

## Check-by-check summary

| # | Check | Result |
|---|---|---|
| 1 | Bridge D = `PROFILE_GRADE_SUBS`, complete, `sigma_W ≢ W_0·eta_bg`, jets retained | **Pass** (import object; N1/N2 prose nits) |
| 2 | Strong-vs-weak whitelist; no HeldDiv drop on strong | **Fail open** — whitelist right; HeldDiv *expand* underspecified (B2); drop still a misuse hazard |
| 3 | Divergence classifier rigor / fixtures / `--drop-divergence` | **Pass** shape; depends on B1 so protected COUPLING rows cannot enter |
| 4 | Protected atom-gating; ENERGY_BASIS barred | **Fail** — atom list right; whole-vs-term route ambiguity (B1) |
| 5 | Leak + accounting DoD | **Pass** |

**Launch recommendation:** fold B1 (ordered Bridge D → protected stop → COUPLING V-classify) and B2 (HeldDiv → formal divergence, never drop) into the directive; re-leg; then build.

**Measurements:**  
`research/pde_ledger_v3/_measurements/s11c_b_adjudication_v2_directive_review_r2_grok/bridge_d_and_protected_atoms.py`  
`research/pde_ledger_v3/_measurements/s11c_b_adjudication_v2_directive_review_r2_grok/bridge_d_and_protected_atoms.stdout`
