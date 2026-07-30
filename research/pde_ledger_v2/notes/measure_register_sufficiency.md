# Is `notes/parameter_register.md` sufficient as the seed for a shared dimensions module?

Measured 2026-07-26 on `ledger-v2-rebuild`. Read-only. Paths relative to `research/pde_ledger_v2/`.

## 1. HEADLINE

| | |
|---|---|
| **Coverage fraction** | **105 / 226 = 46.5 %** of the scripts' distinct dimension-bearing model quantities have a register row (parameter-like subset: 104/202 = 51.5 %) |
| **`pending` count** | **2** (`W_slab` :163, `κ_4` :165) — plus 9 rows asserting *no* `{L,T,M}` at all |
| **Distinct bases in the 43 scripts** | **6** conventions over **4** distinct axis sets |

**Verdict: INSUFFICIENT as the seed.** It is a *parameter* register (knobs), not a *quantity* register. It is accurate where it
speaks — 67 semantic agreements against 4 genuine disagreements, sampled across 21 stages — but it never mentions 53 % of what the
scripts dimension, records only one of six axis conventions, has no entry in the native basis of either extreme stage (042/038), and
its key is a prose table cell, not an identifier. Use it as a **correctness oracle** for the 105 quantities it covers, not as the
source a module is generated from.

## 2. COVERAGE (§1)

**Method.** *Register:* master table `:125–240` (116 rows) + EM-track table `:878–885` (8 rows); the 16 bundling rows expanded
(`C_E, C_B` :136; `{κ,χ,σ_a,σ_L}` :181; `{B̃0..}`+`{Z̃0..}` :185; `{φ,q,j,λ}` :228; …); 9 rows whose dim cell is `—`/"**structural,
dimensionless**" subtracted (:154,:161,:229,:230,:232,:234,:236,:239,:240). *Scripts:* criterion = **a name bound to a dimension**
(posted exponent vector, dim-map entry, or computed by dim algebra) **denoting a model symbol/field/coefficient**; basis constants
(`LENGTH`/`MASS`/`ZERO_DIM`…) and check-labels naming an expression are excluded.

* Register defines **137** distinct quantities in the master table (2 retired, 2 `pending`) + **8** EM-track `indexed_*` rows =
  **145**, spanning **31** of 43 stages in its `Enters` column.
* Scripts manipulate **401** (stage, quantity) pairs = **226** distinct quantities, over the **30 / 43** scripts with dimension
  machinery. 13 have none: 001, 014, 015, 017, 019, 020, 022, 024, 025, 026, 028, 033, 043.

* **105 / 226 covered (46.5 %).** The 121 misses: **74** intermediate/derived composites that were never knobs (`A0/A1/T0/T1` 023,
`F_B0…F_Br`/`U_B`/`U_A2` 037, `P0_raw`/`Yhat`/`mu_hat0` 021, `S_leak`/`W`/`sigma_wall` 006, `M_AB`/`K_AB` 013, `M2`/`K2` 016);
**23** fields/coordinates/operators the register deliberately never lists (`u,v,w,k,omega,r,R,theta,n,psi,b_T,curl_u,r_B`); **12** a
registered knob under a *different key*, or new (`kappa_chi`/`lambda_chi` 044 = register `κ_B`/`a_B`; `A_E` 036/037/038 has no row
at all; `lambda_T`,`eta_a`,`lambda_Sigma`,`E_g`,`qL/qh/qM/mh` 042); **8** convention-variants of a registered knob
(`mu_eta_density`,`T_w_density`,`K_eta_density`,`T_Omega_density` 016; `Z0_ret`/`Z1_ret` 023); **4** benchmark/external constants
(`G` 004:240/021:378/027:205, `c` 021/027).

The gap is **(a) plus a category mismatch**, not (b): only **3** stages have dimension machinery and zero register parameter rows
(**008, 021, 042**); conversely **4** stages carry register rows whose `.py` has no dimension machinery (**015, 017, 024, 033**).

## 3. BASIS (§2) — decisive

**stage042 `(stiffness, length, time)`, fractional — ABSENT in that basis; the register has a gap.** `…stage042_…py:808–848` posts
`stiffness_dim=(1,0,0)`, `length_dim=(0,1,0)`, `frequency_dim=(0,0,-1)`, `charge_dim=(Fraction(1,2),Fraction(3,2),Fraction(-1))`,
keying `B,C,K,qL,qh,qM,mh,cE,cg,Mh,QE,omega,d,r`. `qL/qh/qM/mh` have **no register row**; the rest have `{L,T,M}` rows for
*different* objects (`B` = `M L⁻¹T⁻²` :146 vs `stiffness_dim`; `M_h` = `M L⁻¹` :198 vs `zero_dim` 042:844; `Q_E` :143 asserts **no
dimension** — "R1-deferred `magnitude`" — while the script assigns `charge_dim`). **No normalisation between the bases is stated
anywhere**, neither here nor in `manifests/DIM_ORDER_DECISION.md`; 042 reaches the register only via edges R81–R84 (:347–350).

**stage038 four axes `(M,L,T,E-charge)` — present, normalised into `{L,T,M}` with different values, convention unstated.**
`…stage038_…py:692`: `Dimension = tuple[int,int,int,int]  # M, L, T, E-charge`; honest branch :715–719 posts `A_E=(0,1,0,1)`,
`q_T=(1,0,-1,0)`, `c_γ²=c_E²=(0,2,-2,0)`, `mu_R=(2,1,-4,-1)`. `q_T = M T⁻¹` matches row :231 and the speeds match, but `mu_R = M²L
T⁻⁴E⁻¹` **contradicts** row :137 `μ_R = M L⁻¹T⁻²`, and `A_E = L·E` **contradicts** stage037's `AE_DIM=(3,-2,1)` = `M L³T⁻²`
(`…stage037_…py:614`); `A_E` has **no row of its own**, only the formula in row :238. Two stages computing the same `r_BA =
q_T²c_γ²/(μ_R A_E)` use incompatible unit systems; both land dimensionless, so nothing catches it.

**Full basis enumeration** (axis names in the order the script *stores* them):

| # | axes | order | scripts | fractional exponents |
|---|---|---|---|---|
| 1 | 3 | **`(L,T,M)`** | 002 :69, 004 :107, 005 :184, 006 :116, 007 :147, 009 :444, 030 :173, 031 :141, 032 :1046, 034 :329, 035 :340, 036 :285, 037 :604, 044 :416 — **14** | 004/005/006 via `**Rational(1,2)` |
| 2 | 3 | **`(L,M,T)`** | 011 :149, 012 :193, 013 :176, 016 :158, 018 :163, 021 :162, 023 :868, 027 :793 — **8** | 012 (3/2 corrupt walk), 021 (`mu_hat0` −1/2 :671), 023 (`gU`/`gW` ±1/2 :473–474), 027 (`I25` = 5/2 :208) |
| 3 | 3 | **`(M,L,T)`** | 003 :87, 010 :571, 039 :433, 040 :799, 041 :1084 — **5** | no |
| 4 | **2** | **`(L,T)`** — no mass axis | 008 :462–464 — **1** | no |
| 5 | 3 | **`(stiffness, length, time)`** | 042 :808–829 — **1** | **yes** (1/2, 3/2) |
| 6 | **4** | **`(M,L,T,E-charge)`** | 038 :692 — **1** | no |
| — | 0 | none | 001,014,015,017,019,020,022,024,025,026,028,033,043 — **13** | — |

No symbolic exponents anywhere. Fractional exponents occur in **8** scripts, four of them in `{L,T,M}` — fractions are not an
042-only requirement.

## 4. CONFLICTS (§3)

**Same name, two stages, different dimension** (semantic, order-corrected):

| quantity | locus A | locus B | note |
|---|---|---|---|
| `μ_R` | `stage003…py:875`, `007:437`, `037:674`, `040:814`, `044:449` → `M L⁻¹T⁻²` | `stage038…py:719` → `M²L T⁻⁴E⁻¹` | register has only the first |
| `A_E` | `stage036…py:289`, `037:614` → `M L³T⁻²` | `stage038…py:715` → `L·E` | no register row at all |
| `A_V` | register :227 (032) → `M L³T⁻²` | `stage036…py:395`, `037:680` → `L²T⁻²` | different object, same key |
| `K_eta` | `stage013…py:412` + R29 → `M L⁻¹T⁻²` | `stage016…py:362` → `M L⁻³T⁻²`; `stage023…py:466` → `M T⁻²` | **three** values; register documents 013-vs-023 (:170), **not** 016 |
| `μ_η`, `T_w` | `stage013…py:411,412` → `M L⁻¹`, `M L T⁻²` | `stage016…py:360,361` → `M L⁻³`, `M L⁻¹T⁻²` | line- vs volume-density; register records only 013 |
| `T_Omega` | `stage016…py:363` → `M L⁻³T⁻²` | `stage023…py:467` → `M T⁻²` | register **does** carry both (:182, :170) |
| `M0` | `stage009…py:467`, `010:556` → `T⁻¹` | `stage023…py:460` → `M T⁻¹` | no register dim row to arbitrate |
| `B`,`C`,`K`,`M_h`,`d`,`r`,`Q_E` | `{L,T,M}` rows | `stage042…py:831–848` stiffness basis | see §3 |

**Register vs its own script's locus — 40+ entries sampled across 21 stages** (003,004,006,007,009,010,011,012,013,016,018,021,023,027,030,031,032,034,035,036,037,040,044): **67 agreements**, all semantic/order-corrected — the four
primitives + `c_s0`/`a` (004:179–202, 006:302–303, 011:367, 044:427); the eight stage-003 light-sector constants (003:874–881 ≡
006:314–321); `Ω_w`/`g_ℓ`/`ℓ_g` (007:427–439); `d` (009:477); `k_warp` (010:563); `L0` (011:369); `α` (012:565); the wall packet
(013:409–413); 016:358–366; 023:465–488; the nine port scalars (027:205–216,:385); the eleven electric-closure rows (030:208–222);
the seventeen puncture-deflection rows (031:555–575); 034:375–377/:720; `J_T` (035:840); 036:390–398; the four boost ratios
(037:688–691); `Z_χ` (044:431). **Disagreements — 4:**
1. **`Q_E`** — register :143 asserts *no* dimension; `stage042…py:845` assigns `charge_dim`.
2. **`μ_R` / `A_E` at 038** (above).
3. **`μ_η`/`T_w`/`K_η` at 016** — the volume-density values exist in code and nowhere in the register.
4. **Locus mis-attribution.** Rows :182–:185 (`T_Ω`,`β₂`,`M̃/K̃/T̃_Ω`,`{B̃,Z̃}`) say "stage 017,
dual-engine", but **stage017 has no dimension machinery in either engine** — `…stage017_…py:62` only cites `CITED_016_DIMENSIONAL_OK
= True` (mirrored `…stage017…wl:34`); the values live at `…stage016_…py:355–366`. Values agree, the pointer does not. Same shape at
:174 (`Vp0/ℓ_c = M L⁴T⁻²`, honestly labelled "Codex-verified", but **no engine computes it** — stage015 has zero dim machinery in
`.py` *and* `.wl`) and :181 (`{κ,χ,σ_a,σ_L}` claims "stage 015 … dual-engine" and gives no exponent triple).

## 5. IDENTITY (§4)

The key is **neither stable nor unique**, and is not a machine key: a prose cell mixing symbol, gloss, strikethrough (`~~λ_Pu~~`
:139), escaped pipes (`h_A = P₀H\|_A` :204) and set-braces (`{κ, χ, σ_a, σ_L}` :181). No `quantity_id` field exists.

*One name, several quantities:* `κ` — legacy Hessian constant (:181) vs `κ = D/B_eff` `M L T⁻²` (:215), plus `κ_B`/`κ_4`/`κ_phase`.
`T_Ω` (:182 `M L⁻³T⁻²`) vs `T_Omega` (:170 `M T⁻²`) vs `T̃_Ω` (:184). `K_η` (:179) vs `K_eta` (:170), plus a third at
`stage016…py:362`. `A` (:224) vs `A_eff` (:201) vs `A_X` (:227) vs `A_E` (no row); `A_V` collides across 032 and 036/037. `C` (:225
`M²L⁴T⁻⁴`) vs `C_E`/`C_B`/`C_J`/`C_hu` (:136/:149/:151) vs stiffness `C` (042:832). `c_E` (`L T⁻¹` :135) vs `C_E` (`M⁻¹L⁻⁴T²` :136)
— distinguished only by **case**. `q`/`j`/`λ`/`φ` (:228) vs `q_T`/`J`/`λ_c`/`λγ`.

*One quantity, several names:* `D` (:202) = `sigma` (040:822); `a_B`/`κ_B` (:156/:157) = `lambda_chi`/`kappa_chi` (044:432–433);
`η_i` (:209) = `eta_a` (034/035) = `eta` (044:444); `s_i` (:226) = `s` (034/035) = `s1`/`s2` (036/037). And `ε0/ε1` (:169) vs
`Z0_ret/Z1_ret` (023): the register calls them *the same two dofs* yet they carry **different dimensions** (`1` vs `M T⁻²`).

## 6. WHAT A MODULE SEEDED FROM THIS WOULD BE MISSING

1. **121 of 226 quantities** — every intermediate, every field/coordinate, `A_E`, `G`, `c`, and stage042's four charge couplings.
2. **Five of the six axis conventions.** The header (:70) declares `{L,T,M}` only; nothing represents `(L,T)` (008),
   `(stiffness,length,time)` (042) or `(M,L,T,E-charge)` (038), and no normalisation between them is stated.
3. **A machine key** — no `quantity_id`; names collide and one quantity has several spellings.
4. **`dim_source`/`dim_source_order` per symbol** — the register cites a *stage*, and in ≥4 cases (:174, :181, :182–:185) the cited
   stage is not where the value is computed.
5. **Per-stage scoping** — one flat namespace over a corpus that is not: `μ_η`,`T_w`,`K_η` legitimately differ between the 013 line
   and the 016 volume convention.
6. **The ~15 quantities the register carries that no engine computes** — the 2 `pending` rows, the 9 "structural, dimensionless"
   rows, `Vp0/ℓ_c`, `{κ,χ,σ_a,σ_L}`: nothing to emit, no oracle.
7. **Retired rows** (`λ_Pu` :139, `α_aniso` :159) carry historical dims with no live/retired flag.

## Commands run

```
awk 'NR>=123 && NR<=241 && /^\|/' notes/parameter_register.md     # master table (118 pipe-lines)
grep -n 'pending' notes/parameter_register.md                     # 2 live pending dim rows
for f in scripts/ledger_stage*.py; do grep -cEi '\bdim\b|dimension|Dim\(' "$f"; done  # machinery census
grep -ni 'dim' mathematica/ledger_stage{015,017,024}_*.wl         # 015/024 = 0; 017 = cited boolean
python3 _scratch/pass1_dim_survey/measure/{register_rows,coverage}.py    # the two counts — ⛔ SCRIPTS NOT RETAINED
```

⛔ **The two counting scripts are not retained** — they lived in gitignored
`_scratch/pass1_dim_survey/measure/` and no copy survives, so that line is a record of what was run, not a
re-runnable command. The other four commands above are re-runnable as written.
Per-script inventories came from five read-only agents over the 43 `*_sympy_audit.py` files under the §1 criterion; their raw output
backs §2–§4. The 226/401 counts are exact under that criterion; the *bucketing* of the 121 misses is a judgement call at the margins
(±5 between "intermediate" and "renamed knob").
