# Moving-Throat PDE — Open-Item Reconciliation

**Purpose.** Before any compute on the nonlinear closure, settle a framing question that has been
*asserted* but never *verified*: are the two "open items" the same thing?

- **Object A — the ledger's open item: "actual branch realization."**
  (`moving_throat_pde_handoff_full.md` §13.2; translation dictionary §15.3.)
- **Object B — `pathA_29`'s "open-item #9 / R0=−M0, R1=−D1."**
  (`software/stage1_solver/reports/pathA_28_cancellation_condition.yaml`; `pathA_29_brane_bulk_return.md`.)

`pathA_29`'s directive claims it addresses "pde_ledger open-item #9." **That mapping is asserted in the
stage1_solver directive, not present in the ledger** — `R0=−M0` / "open-item #9" return **zero** matches in
the 253-stage ledger (grep-verified 2026-06-25; the ledger's `S_leak` appears in §4.1 but the ℓ=0/1
*radiation-cancellation* condition does not). So "#9" is `stage1_solver`'s own numbering. This note maps the
two precisely.

> **⚠️ Symbol firewall.** The ledger *does* use the glyph `R0`, but it means the **mixed-port coupling**
> `R_{A,r}` in `Δ_{A,r}=Ω_{U}²Ω_{W}²−R_{A,r}²` (handoff §10.2; stage_251) — **not** pathA_29's ℓ=0 *return
> moment* `R0=−M0`. Same three characters, unrelated objects. This collision is precisely how the two
> framings get falsely conflated; keep them distinct.

---

## 1. What each object actually is

### Object A — ledger "branch realization" = the **ℓ=2 (quadrupole) required channel**
The ledger has finished **all algebraic compression**. Its single open item is the *branch dynamics of the
solved PDE*: produce, on the true moving-throat solution, the grouped real **P2** (ℓ=2) data —

- conservative/transfer moments `(K_A, M_A, B_{A,n}, Z_{A,n}, N_{A,n})` for lanes `A∈{20,21,22}`
  (handoff §10),
- the isotropy gate `a₂=b₂=a₄=b₄=0` (§11),
- the outgoing prefactor `P₀ = N₀/D₀` (§12.2),

and hit the **single sharp normalization target**

```
m̂₀² · P₀ = 54 G c_s⁵ / (5 a⁵ c⁵).            (handoff §12.3)
```

Physically this is the **GR-quadrupole / gravitational-radiation-emission** end — the channel GR *requires*
(leading mass-quadrupole radiation, 2.5PN/4PN). The job is to **match** it.

### Object B — pathA_29 "#9 / R0=−M0, R1=−D1" = the **ℓ=0,1 forbidden channels**
`pathA_28` computed the raw outgoing amplitudes of the *same* multipole expansion and isolated, for the two
**GR-forbidden** sectors, the moments the brane↔bulk return must cancel
(`pathA_28_cancellation_condition.yaml`):

- **ℓ=0 (monopole):** `M0(ω)=∫_brane S_leak d³x`; return `R0=−M0` must cancel the `O(ω⁰)` source moment
  feeding the raw `O(ω¹)` outgoing coefficient. Residual `A0 ∝ i·a·(ω/c_s)·(M0+R0)`.
- **ℓ=1 (dipole):** `D1_i=∫x_i S_leak + ∫j_i` (incl. carried odd wake); `R1=−D1` cancels the `O(ω⁰)`
  dipole/momentum-rate moment feeding the raw `O(ω³)` coefficient. Residual `A1 ∝ i·a³·(ω/c_s)³·(D1+R1)`.

Physically this is the **forbidden-radiation** end — GR forbids monopole GW (Birkhoff) and dipole GW
(momentum conservation). The job is to **cancel** it. `pathA_29` then showed (RETURN_RESIDUAL_PREDICTION)
that the drain *cannot* fully cancel: a bounded residual `∝ ε₀ = 1−𝒯₀(0)` survives → the falsifiable
departure from GR.

---

## 2. The mapping — same expansion, division of labor by ℓ

A and B are **the multipole sectors of one and the same object**, sliced by angular momentum ℓ:

| ℓ | GR status | role of the return | tracked by |
|---|-----------|--------------------|------------|
| 0 | forbidden (Birkhoff) | cancel `M0` → residual `∝ε₀` | pathA_28/29 ("#9") |
| 1 | forbidden (momentum) | cancel `D1` → residual `∝ε₀` | pathA_28/29 ("#9") |
| 2 | **required** (quadrupole) | produce `P₀` → hit `54/5` | **ledger** ("branch realization") |
| ≥3 | finite-size demoted | (roadmap Phase 5) | — |

The shared root structure is explicit in the handoff:

1. **One source.** pathA_28's `M0/D1` are *moments of the very same* `S_leak = −[Wj^w] + ∫W′j^w` that the
   ledger's projected continuity carries (handoff §4.1). Identical operator, projected onto different ℓ.
2. **One PDE, projected by harmonic.** "Projecting [the linearized geometry PDE] onto harmonics gives
   separate modal equations for `l=0`, grouped real `l=2`, and higher multipoles" (handoff §9.4). A and B
   are two of those projections.
3. **One outgoing operator.** pathA_28 says it computed "the raw ell=0/1/2 outgoing amplitudes." ℓ=0,1 are
   the ones that must vanish (B); ℓ=2 is the one that must *survive and be normalized* (A). Same low-ω
   outgoing fingerprint `Ŷ(ω)` (handoff §12.1), read at different ℓ.

**Behind all three: the single brane↔bulk return law `S_return` (and the nonlinear moving-throat PDE that
determines it).** Projected onto ℓ=0 it gives `R0`; onto ℓ=1, `R1`; onto ℓ=2 it feeds the outgoing-transfer
moments `N_{2m,n}` and hence `P₀`.

---

## 3. So: same, or independent? — **Neither. Complementary.**

- **Not identical.** pathA_29's `R0/R1` are the ℓ=0/1 return moments; the ledger's `N_{A,n}` are the ℓ=2
  outgoing-transfer moments. Different angular-momentum channels, different physical demand
  (*cancel* vs *match*). `R0/R1 ≠ N_{A,n}`.
- **Not independent.** They are driven by the *one* return law `S_return` / the *one* nonlinear PDE. You
  cannot solve one end while leaving the other free: the same `S_return` simultaneously sets the ℓ=0/1
  residual *and* the ℓ=2 prefactor.
- **Therefore complementary sectors of one object.** "#9" is an **additional physical requirement** that
  `stage1_solver` surfaced (the forbidden-channel cancellation), sitting *alongside* — not inside — the
  ledger's ℓ=2 open item. The ledger never tracked it because the ledger's radiation focus was always the
  ℓ=2 quadrupole.

### One honest gap this does *not* close
A and B were computed in **different reductions**. pathA_29's `R0/R1` residual came from a *postulated flat
finite-slab* model (brane at `w=0`, absorber at `w=d`), via projected continuity + a solved Helmholtz
transport phase. The ledger's `N_{A,n}` come from the *full mixed-sector port kernel* (handoff §10.3:
`Ω_U, Ω_W, R_{A,r}, g`'s). They agree in **physical role** but are **not yet the same computation**. The
nonlinear closure's actual job is to **unify them**: one solved PDE that reproduces *both* pathA_29's ℓ=0/1
slab residual *and* the ledger's ℓ=2 `N_{A,n}`/`P₀`.

---

## 4. Consequences for the push (what this buys us)

1. **Don't conflate.** "Close pde_ledger #9" overstates pathA_29. The ledger's open item (ℓ=2 branch
   realization → `54/5`) is **still fully open**. pathA_29 closed the *linear-response, flat-slab* version of
   the ℓ=0/1 forbidden-channel question, with a residual-prediction verdict.

2. **The first compute gate is correctly shared infrastructure.** The Phase-1 **frozen-wall D/N unit test**
   `Z₀₀(ω)=−(ω/c_s)·tan(ωL/c_s)` (handoff §7.2, roadmap Phase 1–2) validates the geometry lift / mode
   structure that **every** ℓ-sector rides on. It is upstream of both A and B → right place to start.

3. **A *new* falsifiable cross-ℓ consistency condition.** Because one `S_return` sets all ℓ, the unified PDE
   must *simultaneously*: (a) leave the bounded ℓ=0/1 residual `∝ε₀` that pathA_29 found (don't secretly
   over-cancel → that would erase the GR-departure prediction), and (b) deliver the ℓ=2 `P₀` that hits
   `54/5`. A return law that "cleanly cancels at all ℓ" would be **suspicious**, not a success
   (falsification ethos). This cross-ℓ check is itself an able-to-fail gate the nonlinear closure must pass —
   and it did not exist as a stated test before this reconciliation.

4. **Calibration boundary (unchanged).** `54/5` is a calibration anchor (Gate 4 = GENUINE_BLOCKED): the PDE
   must deliver the *form/branch* of `P₀` and the *isotropy gate*; the overall normalization of `G` stays a
   labeled knob. The `ε₀` residual amplitude is likewise currently uncalibrated. "Thorough, calibrations
   allowed, nothing asserted" applies to both ends.

---

## 5. One-line answer

> pathA_29's "#9 / R0=−M0, R1=−D1" and the ledger's "actual branch realization" are **the ℓ=0/1 and the ℓ=2
> sectors of the same linearized moving-throat outgoing response**, both driven by the one brane↔bulk return
> law / one nonlinear PDE. Same object, **different multipole channels, opposite demands** (cancel the
> forbidden ℓ=0/1; match the required ℓ=2). They are complementary, not identical and not independent — and
> unifying them in a single solved PDE (so it reproduces *both* the pathA_29 ℓ=0/1 residual *and* the ledger
> ℓ=2 `P₀`) **is** the nonlinear closure, with a new cross-ℓ consistency gate falling out for free.

*Cross-refs:* `moving_throat_pde_handoff_full.md` (§4.1, §9.4, §10.3, §12); `moving_throat_pde_roadmap.md`
(Phase 1–6); `software/stage1_solver/reports/pathA_28_cancellation_condition.yaml`;
`software/stage1_solver/reports/pathA_29_brane_bulk_return.md`; memory `project-moving-throat-pde-push`.
