# ledger_stage040_cone_lock_readjudication

## Status

**Part VI — The knit. VI-1 (build-order second; the FIRST stage of the 3-stage
Part-VI split by stage id, built SECOND in the 041→040→042 order, user decision
2026-07-23) — the cone-lock re-adjudication.** All four far-field sectors
(gravity, light, charge, magnetism) live on ONE shared candidate field set. The
gate asks the honest question: are the two calibrated cone locks that would make
those sectors share a single characteristic cone **EARNED facts** of the
committed model, or merely **AVAILABLE calibration choices**? The two locks are

- **Lock A** — `λ_γ=1 ⟺ light cone = gravity-phonon cone ⟺ L_A: m*mu_R − 5*K*rho^4*rho_br = 0`.
- **Lock B** — `c_E=c_gamma ⟺ electric cone = light cone ⟺ L_B: c_E^2*rho_br − mu_R = 0 ⟺ r_cone=1`.

The honest answer is that **neither lock is derived on the earned ledger**: both
are calibration INPUTS with codimension **`Δr=2`**, and the committed model
GENUINELY ADMITS solutions where each lock FAILS (each lock is a `WITNESSED`
non-entailment, witness values `L_A=5`, `L_B=7`; equivalently `r_cone=9/2≠1` at
the lock-B witness). Verdict token (both engines, both exit 0):

- **`CONE_LOCK_CALIBRATED`**, with `delta_r=2` (`dim_before=10 → dim_after=8`),
  both locks `provenance=WITNESSED` (`lock_value_at_witness`: `L_A=5`, `L_B=7`),
  `atomic_riders=()`, Route A `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT`
  (`missing_objects={R1,R2,R3,R4,R5}`), Route B `ROUTE_B_CLOSED_CHECKED_NEGATIVE`,
  freedom `FREEDOM_UNCONSTRAINED{C_hu, rho_br}`.

**⭐ RE-ADJUDICATION (interpretation-level, carried verbatim in intent, NEVER
softened, NEVER inflated).** The exact token `CONE_LOCK_CALIBRATED` STANDS as
computed, but its INTERPRETATION is re-adjudicated to the honest **"NO committed
cone lock"**: the cone locks are **AVAILABLE CALIBRATION CHOICES, NOT earned
facts** — the G0 card + the magnetism build do NOT establish `λ_γ=1` or
`c_E=c_gamma`; `c_s/c_E=2` was CHOSEN expressly to avoid forcing one; and Part
V's `r_cone` left `c_E=c_gamma?` OPEN. This is NOT a token replacement and NOT
"vacuous." Scope class: **CALIBRATED / UNCOMMITTED** (both locks are inputs, not
earned). The card uses `\StatusOpen`, NOT `\StatusExactClosure`.

**⚠ The gate is GENUINELY ABLE-TO-FAIL, not near-vacuous.** The production
verdict reads CALIBRATED **because the committed model admits lock-violating
witnesses** (the two `WITNESSED` non-entailment witnesses `5` and `7`), NOT
because any tooth is a tautology. The controls prove it is able-to-fail from BOTH
directions: the `forced_lock` synthetic-throat-bridge control flips the verdict
to **`CONE_LOCK_DERIVED`** (`Δr=0`, both locks `ENTAILED`); `over_constrained`
fires **`NO_GO(sector-ledger)`**; `freedom_tie` fires **`NO_GO(cone-lock)`**; and
the two single-lock forces give **`CONE_LOCK_PARTIAL(...)`** (`Δr=1`). The verdict
is **load-bearing for Part VII's codimension count** (the `Δr=2` technique).

**⚠ `λ_γ = c_gamma/c_s` is a DERIVED ratio, NOT a free knob** — the free content
is `c_gamma` (ultimately `{mu_R, rho_br}`) and `c_E`. **NEVER assert
`c_E=c_gamma` as settled; `r_cone` (R71, Part V) is the open handle**
(`r_cone = c_E^2/c_gamma^2 = 9/2 ≠ 1` at the witness). Lock B is the Part-V-reopened
lock (wired to `r_cone`, R71); lock A (light ↔ gravity-phonon) is an untouched
separate calibration.

## Purpose

Record, as a first-class CALIBRATED / UNCOMMITTED result, that the two cone locks
are AVAILABLE calibration choices rather than earned facts, with a codimension of
`Δr=2` (two independent locks, NOT one), and that the gate that decides this is
genuinely able-to-fail in both the derived and no-go directions. The decisive
objects are (i) the two per-lock `WITNESSED` provenance facets (each lock its OWN
tooth, wired to its own witness value + — for lock B — Part V's `r_cone`); (ii)
the `Δr=2` codimension (Krull grevlex dimension in `.py`, CAD `RegionDimension`
in `.wl`); (iii) the coupled scalar–vector field overlay
(`det M|cone = −C_hu^2*k^4`, the two mixed poles, and the inherited `OPEN_110`);
(iv) the `{C_hu, rho_br}`-freedom conditionality (citing built stage-041's `C_hu`
certification, `rho_br` still sim-deferred); (v) the earned-vs-calibrated
partition that repairs the source's `earned_equalities` bug; and (vi) the eight
falsifiability controls, each an individually-asserted computed tooth.

Both engines build the source-fact inventory, run the algebraic adjudication, and
reach `CONE_LOCK_CALIBRATED` independently — through the identical pipeline

```text
source-fact inventory  →  Route-A grade + Route-B status + freedom status
                       →  sector SAT + locks SAT (exact positive witnesses)
                       →  per-lock provenance (WITNESSED via non-entailment)  →  lock values 5, 7
                       →  codimension  Δr = dim_before − dim_after = 10 − 8 = 2
                       →  field overlay  det M|cone = −C_hu^2*k^4 + the two poles (+ OPEN_110)
                       →  decide(prepass, algebra)  →  (CONE_LOCK_CALIBRATED, ()).
```

Consumes **nothing new** — it CITES the base medium field set + GNLS substrate
`{rho, K, m, a}` (Part I), the gravity phonon speed `c_s^2 = 5*K*rho^4/m` (Part
II), the light/shear speed `c_gamma^2 = mu_R/rho_br` (Part III / stage003), the
shear-surface inertia/modulus `{rho_br, mu_R}` (Part III / pathA_35, postulated,
dims cited), the Route-A shared nonlinear THROAT solve
(`ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT`, missing `{R1..R5}`,
`REGISTERED_DEFERRED` — the SAME shared solve as gravity `{mu_R,rho_br}` R10/R30,
electric `A_E`/sign R63, magnetism `q_T` R67), built stage-041's `C_hu`-freedom
certification, and Part V's open `r_cone` (R71). It performs none of those solves.

## 1. The two cone locks + the algebraic adjudication (`ROUTE_A_GRADE`, `ROUTE_B_STATUS`, `FREEDOM_STATUS`, `SECTOR_SAT`, `LOCK_SAT`, `PRODUCTION_VERDICT`)

The gate enumerates the source-fact inventory and computes the three prepass
grades:

```text
Route A  →  ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT   (candidate bridge present, {R1..R5} absent)
Route B  →  ROUTE_B_CLOSED_CHECKED_NEGATIVE                    (h and u_T distinct Stage-4 DOF; thin-plate over-import rejected)
freedom  →  FREEDOM_UNCONSTRAINED{C_hu, rho_br}                (the two free-parameter facts)
```

The sector ledger `E` (base equations WITHOUT the locks) is `SAT`, certified by an
exact positive witness against the polynomial equalities + strict stability slack
`sigma>0` (`SECTOR_SAT`); `E` PLUS the two locks is also `SAT`, a common positive
witness satisfying both locks + stability (`LOCK_SAT`). `PRODUCTION_VERDICT` is
the first-match `decide(prepass, algebra)` ASSEMBLED from the computed inputs
`(provA=WITNESSED, provB=WITNESSED, delta_r=2)` per the source branch table — NOT
a stored literal — and the whole `decide`-output PAIR is asserted equal to the
fixed ratified pair `(CONE_LOCK_CALIBRATED, ())` (a spurious production rider makes
the pair differ and fires; this is `computed == fixed`, not a `[]==[]`
self-compare). The nonempty-rider path is exercised inside `VERDICT_REDERIVATION`
by the A-provenance-inconclusive case, which computes `prov_A=INCONCLUSIVE` ⇒
`CONE_LOCK_PROVENANCE_INCONCLUSIVE` with
`atomic_riders=["ENTAILMENT_INCONCLUSIVE(L_A)"]`.

## 2. The two per-lock `WITNESSED` provenance facets + the `r_cone` wire (`PROVENANCE_LOCK_A`, `LOCK_A_WITNESS_VALUE`, `PROVENANCE_LOCK_B`, `LOCK_B_WITNESS_VALUE`)

Each lock has its OWN provenance facet — the load-bearing requirement that neither
lock is earned:

- **Lock A** (`λ_γ=1 ⟺ L_A`) is `WITNESSED`: the Groebner remainder of `L_A`
  modulo the sector ideal is NONZERO **and** an exact real non-entailment witness
  satisfies `E` with `L_A ≠ 0`. At the witness (`m=5, mu_R=2, K=1, rho=1,
  rho_br=1`), `L_A = m*mu_R − 5*K*rho^4*rho_br = 5` (COMPUTED, asserted `== 5`).
  Lock A is an untouched separate calibration (light ↔ gravity-phonon), NOT the
  Part-V one.
- **Lock B** (`c_E=c_gamma ⟺ L_B ⟺ r_cone=1`) is `WITNESSED`: `L_B` is NOT
  entailed by the sector ideal (nonzero remainder + a real witness with
  `L_B ≠ 0`). **Wired to Part V's `r_cone` (R71):**
  `L_B = c_E^2*rho_br − mu_R = mu_R*(r_cone − 1)` with
  `r_cone = c_E^2/c_gamma^2 = c_E^2*rho_br/mu_R`; at the witness
  (`c_E=3, rho_br=1, mu_R=2`), `L_B = 7` and `r_cone = 9/2 ≠ 1` (COMPUTED,
  asserted `L_B == 7` and `r_cone == 9/2`). The electric cone is NOT pinned to the
  light cone — this is the concrete evidence that `c_E=c_gamma` is uncommitted.
  R71 is CITED as provenance (printed, no cross-read to Part V).

A witness value of `0` does NOT establish entailment — it merely ceases to witness
non-entailment; entailment is decided ONLY by the Groebner remainder-zero (`.py`)
/ universal `Resolve[ForAll]` (`.wl`) route, never by a single witness.

## 3. The codimension `Δr=2` (`DELTA_R`)

The cone-lock codimension is `delta_r = dim_before − dim_after = 10 − 8 = 2`,
computed by the SAME algebra the source uses:

- **`.py`** — Groebner grevlex Krull dimension of the sector variety before/after
  adjoining the two locks, real-locus guarded by exact positive witnesses.
- **`.wl`** — CAD `RegionDimension[ImplicitRegion[...]]` before/after.

The assert is `computed delta_r == 2` (with separate computed guards
`dim_before == 10`, `dim_after == 8`) — NOT `delta_r == dim_before − dim_after`
(a banned `X≡X` self-check). Dropping one lock from the "after" ideal (or forcing
one lock `ENTAILED`) makes the drop `≠ 2` and fires. This is the **hidden
multiplicity** result: the two locks are INDEPENDENT (two locks, not one) — the
register's R9 `CODIM-PROVEN` edge and the sobering half of the falsification-first
discipline.

## 4. The coupled scalar–vector field overlay + `OPEN_110` (`FIELD_OVERLAY_DET_ON_CONE`, `FIELD_OVERLAY_POLES`, `OPEN_110_CARRY`)

The coupled scalar–vector dispersion determinant is

```text
det M = (rho_br*omega^2 − B_eff*k^2)*(M_h*omega^2 − K_h*k^2) − C_hu^2*k^4.
```

On the light cone (substitute `omega^2 → (mu_R/rho_br)*k^2` and
`K_h → M_h*mu_R/rho_br`) it COMPUTES to exactly `−C_hu^2*k^4`
(`FIELD_OVERLAY_DET_ON_CONE`; asserted `factor(det_on_cone + C_hu^2*k^4) == 0`).
The two mixed scalar–vector speeds (roots of the determinant) are computed
symbolically to their closed forms (`FIELD_OVERLAY_POLES`), the discriminant
carrying the `4*C_hu^2*M_h*rho_br` coupling term. **The `.wl`'s three source
payload literals (the on-cone determinant + the two `Delta_pole_*_under_B` strings)
are DELETED and replaced by the DERIVED symbolic forms** — the req-(ii) fix.

**⚠ The `c_gamma` guard is DIRECT, not merely the det-residual.** With
`K_h → M_h*mu_R/rho_br` and `q = omega^2/k^2`, the on-cone residual `−C_hu^2*k^4`
vanishes-to-`−C_hu^2` at BOTH diagonal roots `q = mu_R/rho_br` (the correct shear
cone) AND `q = B_eff/rho_br` (a competing uncoupled diagonal cone) — so the
det-residual alone cannot enforce `c_gamma^2 = mu_R/rho_br`. The cone value
actually supplied to the reduction is therefore independently CERTIFIED against the
DIRECT shear relation `rho_br*q = mu_R` (a sub-assert within
`FIELD_OVERLAY_DET_ON_CONE`, the genuine computed replacement for the source's
stripped `c_gamma^2 = mu_R/rho_br` grep). The own-mutation substitutes the
competing diagonal cone `q = B_eff/rho_br`: the `−C_hu^2*k^4` det-residual STILL
passes (both roots give it), while the shear-relation certification FAILS
(`rho_br*(B_eff/rho_br) − mu_R = B_eff − mu_R ≠ 0`) and fires.

**`OPEN_110` (inherited, carried, NOT resolved).** The on-cone determinant
`−C_hu^2*k^4` VANISHES iff `C_hu=0`: the mixed scalar–vector poles sit OFF the cone
with residual `∝ C_hu^2`, and are cone-coincident ONLY IF `C_hu=0`. The computed
coincidence set is exactly `{C_hu=0}` (`Solve[det_on_cone == 0, C_hu]` over `k≠0`
gives `C_hu=0`), status token
`OFF_CONE_under_AB_proportional_to_C_hu_squared_OPEN_110`. This is an INHERITED
open item, carried; it is NOT resolved here.

## 5. The `{C_hu, rho_br}`-freedom conditionality (`CONDITIONALITY_FREEDOM`)

040's non-`NO_GO` verdict is CONDITIONAL on `{C_hu, rho_br}` being FREE. The tooth
asserts the two-branch dependence (a positive/negative pair): under
`FREEDOM_UNCONSTRAINED{C_hu, rho_br}` the verdict is non-`NO_GO`
(`CONE_LOCK_CALIBRATED`); under a `freedom_tie` (C_hu tied to `q_h_sq`) it flips to
`NO_GO(cone-lock)`. The verdict is therefore a genuine FUNCTION of the freedom
status — that dependence is WHY 040 stays conditional.

**⭐ Cite built stage-041 (PRINTED provenance, NOT a cross-read).** Built stage-041
certifies `C_hu` is genuinely free on the current ledger
(`FREEDOM_CERTIFIED_CURRENT_LEDGER{C_hu}`, commit `c663d4a3`), which BACKS the
`C_hu` leg of the freedom assumption. **⚠ `rho_br` remains
`FREEDOM_SIM_DEFERRED{Route-A}` per 041 — so 041 does NOT fully discharge 040's
condition; 040 stays EXPLICITLY CONDITIONAL on the sim-deferred `rho_br`.** The
engine computes its OWN `FREEDOM_UNCONSTRAINED{C_hu, rho_br}` status inline; the
041 `C_hu` cert and the `rho_br` sim-deferred status are CITED (printed) as the
reason the freedom assumption holds for `C_hu` but not yet for `rho_br`. The engine
NEVER opens 041's output. `CONDITIONALITY_FREEDOM` is kept DISTINCT from
`CTRL_FREEDOM_TIE` (the single-control neutralization tooth): this tooth asserts the
two-branch dependence + prints the split provenance.

## 6. The earned-vs-calibrated partition — the source bug fix (`EARNED_VS_CALIBRATED_PARTITION`)

**⚠ Source bug (req v).** The source `pathA_40` `ledger.earned_equalities`
(`.py:813-821`, yaml `:1932-1939`) WRONGLY listed the two UNEARNED locks
`L_A`/`L_B` among the earned equalities. The reshape MOVES them out. The
earned-equalities set is

```text
earned = { c_s^2=5*K*rho^4/m,  c_gamma^2=mu_R/rho_br,  B_eff=rho_B0^2/chi_c,
           K_h=M_h*c_E^2,      B_eff*K_h − C_hu^2 = sigma (sigma>0) },
```

and the calibrated/tested set is `{L_A, L_B}` — the two DISJOINT. The tie to
provenance is the **BICONDITIONAL, computed per-lock in BOTH directions**: for each
lock, `(L_i ∈ earned) ⟺ (prov_i == ENTAILED)` AND
`(L_i ∈ calibrated) ⟺ (prov_i == WITNESSED)`. On production both locks compute
`WITNESSED` ⇒ both live in `calibrated`, neither in `earned`. The converse is also
able-to-fail: applying the SAME rule to the `forced_lock` `ENTAILED` case, both
locks MUST appear in `earned` (`forced_earned = BASE_EARNED ∪ {L_A, L_B}`,
`forced_calibrated = ∅`). Injecting `L_A` (or `L_B`) into `earned` while it is
`WITNESSED` (restoring the source bug) fires the biconditional. **Load-bearing:**
leaving a `WITNESSED` lock in `earned_equalities` would silently re-assert the very
locks the gate proves uncommitted.

## 7. The eight falsifiability controls (`CTRL_*`, `VERDICT_REDERIVATION`)

The eight source controls are the positive/negative-direction falsifiability of the
gate — they are WHAT MAKES IT NOT VACUOUS. Each is an individually-asserted computed
tooth (compute the case's verdict via `decide`, assert `== the named build-native
verdict`), neutralized by its own-mutation. For the two controls whose named verdict
EQUALS the production verdict (`CTRL_ABSENT`, `CTRL_PARTIAL_INVENTORY`), the tooth
asserts the case's distinctive COMPUTED TUPLE `(Route-A grade, missing-set, verdict,
Δr)`, not the verdict alone (a verdict-only assert could not fire).

| Control | Named build-native verdict |
|---|---|
| `CTRL_WELL_POSED` (all `R1..R5` present) | `HALT_ROUTE_A_WELL_POSED` (algebra not run) |
| `CTRL_ABSENT` (pathA_35/36-only, no bridge) | `CONE_LOCK_CALIBRATED` — assert TUPLE `(ROUTE_A_ABSENT, {R1..R5}, CALIBRATED, Δr=2)` |
| `CTRL_PARTIAL_INVENTORY` (R1/R2 present) | `CONE_LOCK_CALIBRATED` — assert TUPLE `(UNDERDETERMINED, {R3,R4,R5}, CALIBRATED, Δr=2)` |
| **`CTRL_FORCED_LOCK`** (⭐ synthetic throat bridge ⇒ both locks `ENTAILED`) | **`CONE_LOCK_DERIVED`** (`Δr=0`, `dim_before=8=dim_after`) |
| `CTRL_A_ONLY_PARTIAL` (force lock A only) | `CONE_LOCK_PARTIAL(derived=A, calibrated=B)` (`Δr=1`) |
| `CTRL_B_ONLY_PARTIAL` (force lock B only) | `CONE_LOCK_PARTIAL(derived=B, calibrated=A)` (`Δr=1`) |
| `CTRL_OVER_CONSTRAINED` (`C_hu^2 = B_eff*K_h + eta_over`) | `NO_GO(sector-ledger)` (sector `UNSAT`) |
| `CTRL_FREEDOM_TIE` (`C_hu^2=q_h_sq`, `q_h_sq*rho_br=2*B_eff*M_h*mu_R`) | `NO_GO(cone-lock)` (`E+locks` `UNSAT`) + `FREEDOM_TIED` |

**`CTRL_FORCED_LOCK` is the KEY not-vacuous control** — a target-blind synthetic
throat bridge (`mu_R=rho_br*tau`, `c_E^2=tau`, `m*tau=5*K*rho^4`) makes BOTH locks
`ENTAILED` and flips the verdict to `CONE_LOCK_DERIVED`. This is what proves the
gate genuinely CAN derive the locks; the production CALIBRATED verdict is therefore
a genuine first-match output, not a foregone conclusion.

`VERDICT_REDERIVATION` mutates a COMPUTED input (a constructed source-fact set /
synthetic relation fed to `decide`), NOT the final token, and asserts it re-derives
the NAMED build-native verdict per the source's own branch table — so all
alternates (`CONE_LOCK_DERIVED`, the two `CONE_LOCK_PARTIAL`, the two `NO_GO`, the
A-provenance-inconclusive `CONE_LOCK_PROVENANCE_INCONCLUSIVE` with the nonempty
rider) are build-native. 040 needs NO authored stage-local alternate token. The
`CONE_LOCK_DERIVED` and the two `NO_GO` rows are the negative-direction
falsifiability; the tooth's own ablation corrupts ONE constructed case's computed
input so it fails to re-derive its NAMED verdict, firing at
`VERDICT_REDERIVATION`'s own assert.

## 8. Ledger accounting — what 040 does NOT do

- **Discharges NO knob.** The two cone locks are recorded **CALIBRATED /
  UNCOMMITTED** (`Δr=2`, both `WITNESSED`) — NOT `DERIVED`, NOT a codim reduction.
  Dressing an `IMPOSED`/uncommitted lock as `DERIVED` would falsely SHRINK the
  irreducible count. The `Δr=2` is a **codimension DIAGNOSTIC** for Part VII's
  count, NOT a discharged reduction edge.
- **The actual nonlinear THROAT solve** that would supply `{R1..R5}` and
  potentially DERIVE the locks is **SIM-DEFERRED** (Part VII's central reduction
  debt) — 040 NAMES the deferred solve as the reason neither lock is earned; it
  does not perform it.
- **The Part-V magnetism `r_cone` closure** (`r_cone=1 ⟺ c_E=c_gamma ⟺ lock B`)
  is OPEN — 040 CITES R71 as the concrete evidence lock B is uncommitted; it does
  not close the cone (that needs the throat solve + `delta_BA=0`).
- **`OPEN_110`** (the mixed scalar–vector cone-coincidence, `∝ C_hu^2`, coincident
  only if `C_hu=0`) is CARRIED as an inherited open item, NOT resolved.
- **Whether the locks are ultimately derivable** is handed to the deferred throat
  solve (Part VII names it). 040 records the CALIBRATED/UNCOMMITTED status + the
  `Δr=2` codimension only.
- **The NG5 one-medium reducibility** (`SECOND_MEDIUM_DRIFT`, the trio
  `{rho_B0,chi_c,C_hu}`) is stage 041 (BUILT, `c663d4a3`); 040 CITES only 041's
  `C_hu`-freedom certification and does not re-adjudicate NG5. **The charge-coupled
  scalar consistency** is stage 042 (built after 040). Not touched here.

## Source-to-stage predicate manifest

Completeness certificate: **no silently-dropped source claim.** The source pathA_40
predicate universe (**22** source predicates: the source's 14 checks + its 8
controls) partitions, **same partition both engines**, into a mutually-disjoint,
engine-qualified four-way map:

```text
{ preserved-folded: 18,  replaced-by-stronger: 1,  scoped-out: 3 }  == 22 source predicates
+ newly-added: 4  (STAGE-ONLY, no source predicate)
```

with the exact source-predicate total (`22`) AND the executable-tooth registry size
(`28`) each asserted computationally (`computed == fixed`, both engines; NOT a
`len==len` self-compare), so a silently-dropped tooth fires it.

- **Preserved-folded (18 source predicates → 23 teeth):** `route_a_inventory` →
  `ROUTE_A_GRADE`; `route_b_check` → `ROUTE_B_STATUS`; `freedom_check` →
  `FREEDOM_STATUS`; `sat_decision(sector)` → `SECTOR_SAT`; `sat_decision(locks)` →
  `LOCK_SAT`; `entailment_status(A)` → `PROVENANCE_LOCK_A` + `LOCK_A_WITNESS_VALUE`;
  `entailment_status(B)` → `PROVENANCE_LOCK_B` + `LOCK_B_WITNESS_VALUE`;
  `dimension_delta` → `DELTA_R`; `decide` → `PRODUCTION_VERDICT` +
  `VERDICT_REDERIVATION`; `field_overlay` → `FIELD_OVERLAY_DET_ON_CONE` +
  `FIELD_OVERLAY_POLES` + `OPEN_110_CARRY` (the `.wl`'s three hardcoded payload
  literals `:268-270` DELETED, replaced by the derived symbolic forms — the
  req-(ii) fix, NOT a separate predicate); the 8 controls → the 8 `CTRL_*` teeth.
- **Replaced-by-stronger (1):** `ledger.earned_equalities` (the source bug listing
  `L_A`/`L_B`) → `EARNED_VS_CALIBRATED_PARTITION` (moves them into the
  calibrated/tested set, tied to `WITNESSED` provenance — req v).
- **Newly-added (4, STAGE-ONLY, disjoint from the 22):** `CONDITIONALITY_FREEDOM`
  (the two-branch freedom dependence + the built-041 `C_hu`-cert /
  `rho_br`-sim-deferred printed provenance — the source carried only a prose
  caveat, no computed tooth, and 041 did not exist); `DIMENSION_HOMOGENEITY`
  (units-restored dim tooth — the source had none); `DUAL_ENGINE_TERMS` (local
  per-engine inventory); `SOURCE_TO_STAGE_MANIFEST` (the partition tooth).
- **Scoped-out (3):** the `require_token`/`import_source_objects` grep + the `.wl`
  `assertContains` filesystem text-scans (`:102-140` / `:34-48`) — dodgeable
  source-greps, barred by the standing rule; provenance CITED in this note instead,
  and the load-bearing `c_gamma^2=mu_R/rho_br` relation re-enforced by the computed
  `FIELD_OVERLAY_DET_ON_CONE` DIRECT shear-relation certification; the
  `argparse`/`--compare` harness + all file-writing; the `comparison_payload`/
  `compare_with_mathematica`/`sha256_json`/`count_agreements` cross-engine payload +
  the `.wl` `Import`/`canon` cross-read (`:275-284`) — replaced by the orchestrator
  arbiter re-run.

This `22` is the SOURCE-MANIFEST arithmetic and is DISTINCT from the **28
EXECUTABLE** stage teeth (23 preserved-folded + 1 replaced-by-stronger + 4
newly-added).

## Verification

- **Dual-engine, both exit 0, 28 executable teeth each, genuinely independent
  implementations.**
  `scripts/ledger_stage040_cone_lock_readjudication_sympy_audit.py` — **SymPy 28
  teeth**.
  `mathematica/ledger_stage040_cone_lock_readjudication_mathematica_audit.wl` —
  **Mathematica 28 teeth** (`ROUTE_A_GRADE`, `ROUTE_B_STATUS`, `FREEDOM_STATUS`,
  `SECTOR_SAT`, `LOCK_SAT`, `PROVENANCE_LOCK_A`, `LOCK_A_WITNESS_VALUE`,
  `PROVENANCE_LOCK_B`, `LOCK_B_WITNESS_VALUE`, `DELTA_R`, `PRODUCTION_VERDICT`,
  `FIELD_OVERLAY_DET_ON_CONE`, `FIELD_OVERLAY_POLES`, `OPEN_110_CARRY`,
  `EARNED_VS_CALIBRATED_PARTITION`, `CONDITIONALITY_FREEDOM`, `CTRL_WELL_POSED`,
  `CTRL_ABSENT`, `CTRL_PARTIAL_INVENTORY`, `CTRL_FORCED_LOCK`, `CTRL_A_ONLY_PARTIAL`,
  `CTRL_B_ONLY_PARTIAL`, `CTRL_OVER_CONSTRAINED`, `CTRL_FREEDOM_TIE`,
  `VERDICT_REDERIVATION`, `DIMENSION_HOMOGENEITY`, `DUAL_ENGINE_TERMS`,
  `SOURCE_TO_STAGE_MANIFEST`). Standalone, print-only, assert-zero
  (`raise SystemExit(1)` / `Exit[1]`), no argparse harness, no JSON/YAML payload,
  **zero file-I/O between engines**. Each engine builds its own source-fact
  inventory, computes its own verdict / `delta_r` / lock provenances / lock values /
  freedom / field overlay, and reaches `CONE_LOCK_CALIBRATED` on its own —
  cross-engine agreement is arbiter-confirmed by the orchestrator re-running BOTH
  engines (NOT an in-script compare, NOT one engine reading the other). **No
  dual-engine disagreement.**
- **The `.wl` is a genuinely INDEPENDENT implementation.** 040 is a
  GENUINELY-SYMBOLIC stage; the two engines use materially distinct CAS primitives.
  The `.py` computes the codimension via **Groebner grevlex Krull dimension + exact
  positive witnesses** and provenance via **Groebner remainder-zero ideal
  membership**; the `.wl` computes the SAME facts via CAD-backed
  **`RegionDimension[ImplicitRegion[...]]`** and
  **`Resolve[Exists/ForAll]` / `FindInstance`** (+ `FullSimplify`/`Solve` for the
  field overlay). These are materially distinct routes to `delta_r=2` and to the
  `WITNESSED` provenances — KEPT distinct, NOT a line-by-line port of either the
  stage `.py` or the source `.wl`. The `.wl`'s three source field-overlay payload
  literals are DELETED and re-derived; no engine hardcodes the verdict / `delta_r` /
  lock provenances / lock values / field overlay, and neither reads the other's
  output. (This "genuinely-symbolic independence" scoping was FLAGGED for + CONFIRMED
  at the Codex→Grok→Codex bookend, with a side-by-side transliteration screen of the
  `.wl` against BOTH the stage `.py` and the source `.wl`.)
- **Tri-review CLEAN, ZERO remediation.** The orchestrator arbiter re-ran both
  engines (both exit 0, both produce the identical verdict `CONE_LOCK_CALIBRATED` /
  `delta_r=2` (`10→8`) / both locks `WITNESSED` with values `5`/`7` / `r_cone=9/2` /
  `FREEDOM_UNCONSTRAINED{C_hu,rho_br}` / field overlay `−C_hu^2*k^4`; each tooth's
  mutation fires); **FIDELITY_CLEAN** (all 10 faithfulness items PASS — the
  provenances genuinely compute `WITNESSED`, the codimension is a real Krull/CAD
  count, the field overlay is derived not hardcoded, the DIRECT shear-relation
  certification is genuine, the earned/calibrated biconditional tracks provenance
  both ways, and the `earned_equalities` source bug is repaired);
  **ABLATION_CLEAN** (28/28 teeth fire at their OWN assert in BOTH engines, verified
  via `FIRED_AT_OWN_ASSERT` under `LEDGER_STAGE040_MUTATION=<PRED>`, 56 mutation runs
  across the two engines; the shear guard is genuine — the competing-cone mutation
  passes the det-residual but fails the shear relation; the `.wl` is independent, no
  cross-read).
- **NO can't-fail conjuncts.** Every headline object — the verdict, `delta_r`, the
  two lock provenances (`WITNESSED` vs `ENTAILED`), the two lock witness-values
  (`5`/`7`), the freedom status, the field overlay, the earned-vs-calibrated
  partition, and every control verdict — is **COMPUTED from the real algebraic
  objects** and corruptible by a mutation; the required able-to-fail shape is
  `computed_actual == fixed_ratified_expected` (e.g. the assembled verdict against
  the fixed ratified token, `delta_r == 2`, `lock_value_A == 5`,
  `lock_value_B == 7`), never `stored == stored`. Unfalsifiable provenance (source
  line refs, the `r_cone`/041-cert citations, the honest caveat prose) is PRINTED,
  not asserted. The `.wl` retains no source hardcode and performs no `.py`
  cross-read.

## Downstream consumers

- **Part VI continues** after 040: **042 (VI-3 — the charge-coupled scalar
  consistency, `SCALAR_DEPARTURE_MAPPED_MAGNITUDE_SIM_GATED`)**. The Part-VI knit is
  041 (VI-2, the reducibility gate, built first) → 040 (VI-1, this stage) → 042
  (VI-3).
- **Parameter register:** the cone locks `λ_γ=1` (lock A) / `c_E=c_gamma` (lock B,
  ⟺ `r_cone=1`) recorded **CALIBRATED / UNCOMMITTED** (`Δr=2`, both `WITNESSED`,
  values `5`/`7`), the inherited `OPEN_110` (`C_hu^2` off-cone), and the
  `{C_hu,rho_br}`-freedom conditionality (041 certifies `C_hu`
  `FREEDOM_CERTIFIED_CURRENT_LEDGER`, `rho_br` `FREEDOM_SIM_DEFERRED`) — edges
  **R78–R80**, refining the seeded cone-lock edges R7–R11 with the built 040
  re-adjudication. Each **discharges NO knob**; the `Δr=2` is a codimension
  DIAGNOSTIC (not a discharged edge), load-bearing for Part VII's codimension count.
  ⚠ The seeded R7/R8 class the locks `IMPOSED`; 040's re-adjudicated finding is "NO
  committed cone lock" (AVAILABLE / UNCOMMITTED — the model admits lock-violating
  witnesses), a reconciliation nuance recorded in the new rows.
- **`docs/model_map.md` §3.6 / §4:** the cone-lock re-adjudication bullet — the
  light/gravity-phonon lock and the electric/light lock are AVAILABLE calibration
  choices, not earned facts (`Δr=2`, `r_cone` open); the dynamical throat solve that
  would derive them is deferred.
- **Part VII:** the `Δr=2` codimension enters the honest codimension count (the
  scaled pathA_40 technique); `c_E=c_gamma` stays OPEN (`r_cone`, R71) until the
  shared throat solve runs; and the whether-the-locks-are-derivable question is
  upgraded (or not) by that solve.

## Provenance

- **Physics source:** `software/stage1_solver/tools/pathA_40_cone_lock_sympy.py`
  (the lean algebraic cone-lock decision: `LOCKS` dict, `import_source_objects`,
  `route_a_inventory`, `route_b_check`, `freedom_check`, `sat_decision` /
  `unsat_certificate` / `witness_for`, `entailment_status` [Groebner remainder-zero
  `ENTAILED` / non-entailment `WITNESSED` + `lock_value_at_witness`],
  `groebner_dimension` / `dimension_delta` [`delta_r = 10 − 8 = 2`], `decide`,
  `field_overlay`, `build_results` [⚠ the `earned_equalities` source bug listing
  `L_A`/`L_B`], `build_case_objects` + the 8 controls) +
  `software/stage1_solver/tools/pathA_40_cone_lock.wl` (the already-independent
  re-derivation: `Resolve[Exists/ForAll]`, `FindInstance`,
  `RegionDimension[ImplicitRegion[...]]`, the field-overlay symbolic derivation —
  ⚠ but with three hardcoded payload string literals and an `Import`/`canon`
  cross-read) + reports `pathA_40_cone_lock.md` / `..._results.yaml`. **The
  `.py`/`.wl` are AUTHORITATIVE over the report prose**, which carries the
  `earned_equalities` source bug. The reshape strips the cross-engine payload
  machinery, the filesystem grep-guards, and the file-writing; computes (does not
  hardcode) the `.wl` field-overlay forms; gives each lock its own provenance facet;
  wires lock B to `r_cone` (R71); carries `OPEN_110` + the freedom conditionality;
  and moves `L_A`/`L_B` out of `earned_equalities`.
- **Consumes:** nothing new — cites the Part-I substrate `{rho,K,m,a}`, the
  gravity phonon speed `c_s^2 = 5*K*rho^4/m` (Part II), the light/shear speed
  `c_gamma^2 = mu_R/rho_br` + `{rho_br, mu_R}` (Part III / stage003, pathA_35, dims
  cited), the Route-A shared throat/cone-lock solve (`REGISTERED_DEFERRED`), built
  stage-041's `C_hu`-freedom certification (`c663d4a3`, printed provenance), and
  Part V's open `r_cone` (R71).
- **Cites (provenance, NOT re-derived, NOT re-adjudicated, NOT re-counted; de-dup
  deferred to Part VII):** `{rho_br,mu_R}` (Part III / pathA_35, R10/R30-adjacent),
  `c_gamma` (stage003), `c_s` (Part II), `c_E` (stage030, R7/R8), `C_hu` (stage030),
  `B_eff` (stage003), the Route-A shared throat solve, and Part V's `r_cone` (R71).
- **Governing:** `notes/ledger_v2_blueprint.md` §5 (reshape spec) + §6 (per-tooth
  ablation); `notes/part6_knit_atomic_split.md` (the RATIFIED 3-stage split, 040 =
  VI-1 = the cone-lock re-adjudication; build order 041→040→042; the § FRESH-READ
  FINDINGS reshape-REQUIREMENTS (040) block (i)–(v));
  `research/pde_ledger/notes/MATHEMATICA_MIRROR_POLICY.md`; `docs/model_map.md` §3.6 + §4 (the departure
  ledger). Reshape directive + review trail:
  `research/pde_ledger_v2/_scratch/stage040_reshape_directive.md`.
