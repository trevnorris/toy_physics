# Census rows — stage016 (ℓ=2 SO(3) covariance), both engines

**Status:** PILOT census rows, built to `CENSUS_SCHEMA.md`. Two artifacts (§7.2):

- `research/pde_ledger_v2/scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py` — below, **PY**
- `research/pde_ledger_v2/mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl` — below, **WL**

⛔ This file records provenance (§3's three axes). It does **not** rule on dimensional correctness (§11).
The stage's own `21 of 21 CORRECT` verdict and its `12 typed / 9 computed` split are **dimensional**
verdicts on a **different row universe** — 20 of those 21 records are out of this census's universe
(§6 below). Nothing here was steered toward them.

**Method.** Every locus quoted below was opened and read in this pass (§9.1 rule 1). No attribution was
taken from the stage note, the parameter register, or the schema without opening the file it points at —
and doing so found two stale pointers (§5, D4).

---

## 1. The reported-result set for this stage (§5.8 output 17)

Taken **unchanged** from `notes/census/REPORTED_RESULTS.md` §1, which is the fixed starting set (§7.1.1
step 1). ⛔ Not amended by this pass.

| id | claim (abbreviated) | loci |
|---|---|---|
| `RES:016:L2-IRREP` | the ℓ=2 sector is one 5-dim SO(3) irrep — all five real ℓ=2 channels share one `−Δ_S²` eigenvalue `λ_m = 6`, so the ℓ=2 response is m-degenerate | note `.../notes/stages/ledger_stage016_l2_so3_covariance.md:43-53`; PY `:287`, `:290`, `:299`, asserted `:679-682`; WL `:180-186`, `:187-193`, `:202`, asserted `:466-472` |
| `RES:016:K2-FORM` | `K₂ = K̃ + λ_m·T̃_Ω`, `M₂ = M̃`, with the `T̃_Ω` coefficient being the computed eigenvalue rather than an independently chosen constant | note `:55-60`; PY `:304-305`, `:499`, `:500-504`, asserted `:693-697`; WL `:216`, `:219`, `:220-225`, asserted `:481-482` |

---

## 2. The universe, and the builder decisions taken to execute §7.1.1 / §7.2

The universe is the transitive input closure of those two results (§7.1.1). Walking them reaches
**16 quantities** at **38 binding sites** across the two engines. Every row below carries its
reachability witness.

⭐ **Builder decisions (§13 `BUILDER-DECISION`), stated so they are contestable and applied uniformly:**

- **B-D1 — alias stores are not separate occurrences.** A pure re-bind of an already-computed value
  inside the same expression block (PY `:288-289` `neg_laps[name] = …` / `lambdas[name] = rayleigh`) is
  folded into the computing binding site. Without this, PY's two-step store scores two occurrences where
  WL's single `AssociationThread` (`:179-186`) scores one, inflating one engine by transliteration style
  alone. A re-bind that is **separately consumed** downstream (PY `:506` `lambda_ref`, WL `:227`
  `lambdaRef`) **is** counted, per §7.2's plain text.
- **B-D2 — a loop that assigns one quantity across five channels at one line is one occurrence.** PY
  `:286`/`:287` and WL `:179`/`:180` bind `(−Δ_S²)Y_m` and `λ_m` for all five m at a single site with
  identical code; the reported result is precisely that the value is the **same** across m. The five
  harmonics, by contrast, sit at five distinct literal lines and are five occurrences.
- **B-D3 — a free-symbol declaration is a binding site.** PY `:219-221` and WL `:21` declare `M̃`, `K̃`,
  `T̃_Ω` and assign them **no value**. §7.2 defines a binding site as where an artifact "assigns or
  asserts a value", which does not cover this; but those quantities are in `RES:016:K2-FORM`'s closure,
  and dropping them would let an in-universe quantity leave the count silently, which §5.3 forbids
  outright. So they get rows. See §5 C2 for the axis-B consequence.
- **B-D4 — scope is the two engines only.** The stage **note** is also an artifact under §7.2 (its
  markdown table rows are binding sites), and it is **not censused here**. See §8.

---

## 3. The rows

Field order per row: quantity · binding site · reachability witness · axis A + §9 evidence · axis B +
§9 evidence · axis C closure/leaves/`PHYSICS-FED`/`SELF-REFERENTIAL` · `CONVENTION-LADEN` + evidence ·
`is_tier` · `should_be_tier` / `should_be_basis` / `delta` · LIVE|RETIRED.

⚠ **Read every `is_tier` below together with §5 C1.** A schema conflict between §3.3+§5.7 and §5's
framing decides `is_tier` for **22 of 38** rows; the primary reading is used throughout and the
alternative is counted in §4.

---

### 3.1 `QID:Y20`

#### `QID:Y20@scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py#270`
- **binding site** — `.../scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:270`, the literal
  `sp.sqrt(sp.Rational(5,16)/sp.pi)*(3*sp.cos(theta)**2 - 1)` inside `harmonics_l2()` (`:268-275`).
- **reachability witness** — `RES:016:L2-IRREP`, hop 1: `:270` → `y_expr` (`:285`) → `rayleigh` (`:287`),
  the Rayleigh quotient the result's `λ_m` **is**. Also `RES:016:K2-FORM` at hop 2 via `λ_m` → `:499`.
- **axis A** — `A-UNADJUDICATED`. No reduction is performed (typed literal at `:270`); **no route is
  recorded** for the harmonics anywhere in the two engines or the note (note `:27-41` writes the five
  expressions down and asserts orthonormality; it names no derivation of them); no benchmark; no
  postulate is stated as a defining property of the medium; and no framework property is named as
  foreclosing a route, so `A-IRREDUCIBLE-STRUCTURAL`'s §9 evidence is absent. §9.0: axis-A evidence
  missing ⇒ `A-UNADJUDICATED`.
- **axis B** — `B-DECLARED-LITERAL`. Code locus of the declaration `:270`.
- **axis C** — closure `{C-SELF@:270}`. Walk terminates immediately at the occurrence's own literal.
  `PHYSICS-FED = false`. `SELF-REFERENTIAL = false` (§3.3: immediate termination is a plain literal).
- **CONVENTION-LADEN** — `UNADJUDICATED`. Clause (a) **is** met: a transformation group is documented
  for this occurrence — SO(3) acting on the ℓ=2 irrep, which is the stage's own subject (note `:24-25`,
  `:52-53`; PY `:825`). Clause (b) is **not demonstrated**: no admissible alternative basis (a rotated
  orthonormal frame of the same irrep) is ever run. The only basis perturbations executed are the
  **inadmissible** corruptions at `:774` and `:801`, which are controls, not invariance demonstrations.
  ⇒ `convention-unadjudicated` bucket; buys no exclusion (§3.4).
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `no-tier:unadjudicated` / `none` / **no**.
- **LIVE**.

#### `QID:Y20@mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl#166`
- **binding site** — `.../mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl:166`,
  `"20" -> Sqrt[5/(16 Pi)] (3 Cos[theta]^2 - 1)` in the `harmonics` association (`:165-171`).
- **reachability witness** — `RES:016:L2-IRREP`, hop 1: `:166` → `ys` (`:172`) → `localLambdas` (`:180-186`).
- **axis A** — `A-UNADJUDICATED`, same evidence position as the PY twin: literal at `:166`, no route, no
  benchmark, no stated postulate, no named foreclosing property.
- **axis B** — `B-DECLARED-LITERAL`, locus `:166`.
- **axis C** — closure `{C-SELF@:166}`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `UNADJUDICATED` (clause (a) met per WL `:636`; clause (b) undemonstrated —
  the only basis perturbations are the controls at `:582`, `:600`).
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `no-tier:unadjudicated` / `none` / **no**.
- **LIVE**.

### 3.2 `QID:Y21c` — `PY#271`, `WL#167`

Both occurrences carry **exactly** the field values of §3.1 (same axis A with the same missing
evidence, `B-DECLARED-LITERAL` at the cited line, closure `{C-SELF}` ⇒ `PHYSICS-FED = false`,
`SELF-REFERENTIAL = false`, `CONVENTION-LADEN = UNADJUDICATED` on clause (b),
`is_tier = should_be_tier = no-tier:unadjudicated`, basis `none`, delta **no**, **LIVE**), with these
loci substituted:

- `QID:Y21c@scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py#271` — literal
  `-sp.sqrt(sp.Rational(15,4)/sp.pi)*sp.sin(theta)*sp.cos(theta)*sp.cos(phi)`. Witness: `RES:016:L2-IRREP`
  hop 1 via `:285`→`:287`.
- `QID:Y21c@mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl#167` — literal
  `-Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Cos[phi]`. Witness: hop 1 via `:172`→`:180-186`.

⚠ The leading minus sign is the Condon–Shortley phase — a second documented convention on the same
occurrence; it does not change the `UNADJUDICATED` flag, which is already blocked at clause (b).

### 3.3 `QID:Y21s` — `PY#272`, `WL#168`

All fields as §3.2, loci `:272` (`-sp.sqrt(sp.Rational(15,4)/sp.pi)*sp.sin(theta)*sp.cos(theta)*sp.sin(phi)`)
and `:168` (`-Sqrt[15/(4 Pi)] Sin[theta] Cos[theta] Sin[phi]`).

### 3.4 `QID:Y22c` — `PY#273`, `WL#169`

All fields as §3.2, loci `:273` (`sp.sqrt(sp.Rational(15,16)/sp.pi)*sp.sin(theta)**2*sp.cos(2*phi)`)
and `:169` (`Sqrt[15/(16 Pi)] Sin[theta]^2 Cos[2 phi]`). No Condon–Shortley sign on this pair.

### 3.5 `QID:Y22s` — `PY#274`, `WL#170`

All fields as §3.2, loci `:274` (`sp.sqrt(sp.Rational(15,16)/sp.pi)*sp.sin(theta)**2*sp.sin(2*phi)`)
and `:170` (`Sqrt[15/(16 Pi)] Sin[theta]^2 Sin[2 phi]`).

---

### 3.6 `QID:S2-measure` — the S² integration measure

#### `QID:S2-measure@scripts/…_sympy_audit.py#230`
- **binding site** — `.../scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:230`, `integrate_s2`;
  the measure itself is the literal `sp.sin(theta)` factor plus the limits `(phi, 0, 2*sp.pi)`,
  `(theta, 0, sp.pi)` at `:233-235`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1: the Rayleigh quotient at `:287` is a ratio of two
  `integrate_s2` calls; `λ_m` has no value without this measure.
- **axis A** — `A-UNADJUDICATED`. Typed at `:233-235`; no reduction; no route recorded; no benchmark. A
  `A-IRREDUCIBLE-POSTULATE` reading ("the throat's angular section is a round S²") was tested against §9
  and **fails**: no locus in either engine, the note, or `parameter_register.md` states that roundness as
  a defining property of the medium. Note `:180` mentions "an isotropic wall" only inside an explicitly
  `UNDETERMINED` remark about `T_Ω`.
- **axis B** — `B-DECLARED-LITERAL`, locus `:233-235`. (`sp.integrate` executes here, but what it
  produces is `λ_m`/`Gram`, not the measure.)
- **axis C** — closure `{C-SELF@:233-235}`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
  ⚠ See §5 C3: this leaf is simultaneously a same-artifact literal and a pure-mathematical object
  (`C-MATH`); the tie is broken toward `C-SELF` by a stated rule and neither tag confers `PHYSICS-FED`.
- **CONVENTION-LADEN** — `false` (**not a candidate convention**, §3.4 third option). No transformation
  group is *documented* for the measure anywhere in this corpus; the note `:36` presents it as "the
  genuine measure", i.e. as fact, not as a choice.
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` /
  **DELTA**. Basis, kept legible as a hunch: in the finished model the wall's angular section and its
  measure follow from the throat geometry rather than being typed; ⛔ no route to it is named at any
  locus, which is exactly why the basis is not `named-route`.
- **LIVE**.

#### `QID:S2-measure@mathematica/…_audit.wl#107`
- **binding site** — `.../mathematica/ledger_stage016_l2_so3_covariance_mathematica_audit.wl:107`,
  `intS2`; measure literal `Sin[theta]` and limits `{phi, 0, 2 Pi}`, `{theta, 0, Pi}` at `:109-111`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1 via the Rayleigh ratio at `:183`.
- **axis A** — `A-UNADJUDICATED` (same evidence position as the PY twin).
- **axis B** — `B-DECLARED-LITERAL`, locus `:109-111`.
- **axis C** — closure `{C-SELF@:109-111}`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false` (not a candidate convention).
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

---

### 3.7 `QID:neg-Delta-S2` — the `−Δ_S²` Laplace–Beltrami operator

#### `QID:neg-Delta-S2@scripts/…_sympy_audit.py#239`
- **binding site** — `:239`, `laplacian_s2`; the operator is typed at `:241-242` as
  `(1/sin θ)∂_θ(sin θ ∂_θ·) + (1/sin²θ)∂_φ²·`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1: `:239` → `neg_lap` (`:286`) → `rayleigh` (`:287`).
  Also `RES:016:K2-FORM` hop 3 (via `λ_m` → `build_K2` at `:499`).
- **axis A** — `A-UNADJUDICATED`. Typed at `:241-242`; no reduction performed; **no model-equation locus
  is cited for it** anywhere in either engine (the PY docstring `:2-10` names `pathA_32` as the *source
  script* it was reshaped from, and note `:16-18` states the derivation is "self-contained"; neither is a
  citation of an equation the operator follows from). `A-IRREDUCIBLE-STRUCTURAL` was tested and fails
  §9: no locus names the framework property that forecloses deriving the wall's angular operator (the
  deferrals recorded at note `:311-316` are about the **radial** profile and the Gate-4/5/6 sim, not this
  operator), so naming one here would be manufacturing evidence.
- **axis B** — `B-DECLARED-LITERAL`, locus `:241-242`.
- **axis C** — closure `{C-SELF@:241-242}`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
  ⚠ This is the **decisive** axis-C call of the stage: by §3.3.1(1), a routine's operator is a
  `C-FIELDEQ` leaf **only if the artifact cites the locus of the equation being solved**. It does not.
  ⇒ no `C-FIELDEQ` leaf anywhere downstream of it. §3.3.1(1)'s own test-that-this-can-fail — "a solve
  assembled entirely from the artifact's own constants, with no cited equation, is **not**
  `PHYSICS-FED`" — is the case here.
- **CONVENTION-LADEN** — `false` (not a candidate convention: no transformation group is claimed for the
  operator).
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` /
  **DELTA**. Basis: the model's picture has the ℓ=2 angular operator as the angular part of the wall's
  own equation of motion, which would make it (and everything below it) physics-fed. ⛔ Hunch — no route
  is recorded.
- **LIVE**.

#### `QID:neg-Delta-S2@mathematica/…_audit.wl#115`
- **binding site** — `:115`, `lapS2`; operator typed at `:116-118`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1 via `:179` → `:183`.
- **axis A** — `A-UNADJUDICATED` (same evidence position; WL cites no equation locus either — see the
  header comment `:1-8` and the print at `:460`, which names the operator but cites nothing).
- **axis B** — `B-DECLARED-LITERAL`, locus `:116-118`.
- **axis C** — closure `{C-SELF@:116-118}`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

---

### 3.8 `QID:neg-Delta-S2-applied` — the value `(−Δ_S²)Y_m`

#### `QID:neg-Delta-S2-applied@scripts/…_sympy_audit.py#286`
- **binding site** — `:286`, `neg_lap = compact(-laplacian_s2(y_expr))` (loop over the five channels; one
  site by B-D2).
- **reachability witness** — `RES:016:L2-IRREP`, hop 1: feeds the Rayleigh numerator at `:287` and the
  eigenfunction residual at `:290`, both listed as loci of the result.
- **axis A** — `A-REDUCED`. §9 evidence: **the reduction is performed at `:286`**, and it reduces to
  `Y_m` (`:270-274`) and the operator (`:239`, body `:241-242`) — both loci given.
- **axis B** — `B-EXECUTED`. Code locus `:286` (calls `sp.diff` at `:241-242`); input leaves `:270-274`
  and `:241-242`; consumed by the assert at `:680`.
- **axis C** — closure `{C-SELF@:270-274, C-SELF@:241-242}`. `PHYSICS-FED = false` — the closure is fully
  determined and contains no `C-FIELDEQ` and no `C-EXTERNAL` leaf (§3.7 above). `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false` (not a candidate convention).
- **is_tier** — `no-tier:unclassified-nonfed` (§5.1: `A-REDUCED ∧ ¬PHYSICS-FED`). Also §4 near-miss
  bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` /
  **DELTA** (basis as in §3.7 — if the operator were physics-fed this row would be tier 3).
- **LIVE**.

#### `QID:neg-Delta-S2-applied@mathematica/…_audit.wl#179`
- **binding site** — `:179`, `localNegLaps = AssociationThread[localOrder, clean[-lapS2[#]] & /@ localYs]`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1: feeds `:183` and the residual at `:190`.
- **axis A** — `A-REDUCED`; reduction performed at `:179`, reducing to `:166-170` and `:116-118`.
- **axis B** — `B-EXECUTED`; code locus `:179` (native `D` at `:116-118`); leaves `:166-170`, `:116-118`;
  asserted `:467`.
- **axis C** — closure `{C-SELF@:166-170, C-SELF@:116-118}`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

---

### 3.9 `QID:lambda-m` — the `−Δ_S²` eigenvalue (4 occurrences)

#### `QID:lambda-m@scripts/…_sympy_audit.py#287`
- **binding site** — `:287`,
  `rayleigh = compact(integrate_s2(y_expr*neg_lap)/integrate_s2(y_expr**2))`.
- **reachability witness** — `RES:016:L2-IRREP`, **hop 0**: this is the result's own quantity, asserted
  `:679` and `:681`. Also `RES:016:K2-FORM` hop 1 (`build_K2(lambdas[name])`, `:499`).
- **axis A** — `A-REDUCED`. §9 evidence: reduction performed at `:287`; reduces to `Y_m` (`:270-274`),
  the operator (`:241-242`) via `neg_lap` (`:286`), and the measure (`:233-235`).
- **axis B** — `B-EXECUTED`. Code locus `:287`; input leaves `:270-274`, `:241-242`, `:233-235`;
  asserted `:679-682`.
- **axis C** — closure `{C-SELF@:270-274, C-SELF@:241-242, C-SELF@:233-235}`. `PHYSICS-FED = false`
  (determined; no `C-FIELDEQ`/`C-EXTERNAL` leaf — §3.7). `SELF-REFERENTIAL = false` (the walk never
  returns to `:287`).
- **CONVENTION-LADEN** — `false` (not a candidate convention: no transformation group is claimed for the
  eigenvalue; the SO(3) group documented for the *basis* leaves this value invariant, which is the
  result, not a convention claim about it).
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` /
  **DELTA**. ⚠ Counter-reading recorded, not adopted: if `−Δ_S²` is held to be pure mathematics
  (`REPORTED_RESULTS.md:78-80` — "The number 6 is `ℓ(ℓ+1)` and is mathematics"), then no physical input
  will ever feed it and `should_be_tier` collapses to `is_tier` with basis `none`. The adopted reading
  follows the same file's test-(b) reasoning (`:44-47`): the eigenvalue **becomes a coefficient of the
  medium's ℓ=2 equation of motion**, so the picture does expect it to be physics-fed.
- **LIVE**.

#### `QID:lambda-m@scripts/…_sympy_audit.py#506`
- **binding site** — `:506`, `lambda_ref = lambdas["20"]`. Counted separately under B-D1 because it is
  separately consumed (`:508`, `:509`).
- **reachability witness** — `RES:016:K2-FORM`, hop 1: `:506` → `k2_ref = build_K2(lambda_ref)` (`:508`).
- **axis A** — `A-REDUCED`; the reduction is performed at `:287` (cited), reducing to `:270-274`,
  `:241-242`, `:233-235`.
- **axis B** — `B-EXECUTED`; the value was produced by executed code at `:287`; leaves as above.
- **axis C** — closure as `#287`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

#### `QID:lambda-m@mathematica/…_audit.wl#180`
- **binding site** — `:180-186`, `localLambdas = AssociationThread[…]`; the Rayleigh quotient itself at
  `:183`.
- **reachability witness** — `RES:016:L2-IRREP` hop 0 (asserted `:466`, `:471`); `RES:016:K2-FORM` hop 1
  (`:219`).
- **axis A** — `A-REDUCED`; reduction performed at `:183`; reduces to `:166-170`, `:116-118`, `:109-111`.
- **axis B** — `B-EXECUTED`; code locus `:183`; leaves `:166-170`, `:116-118`, `:109-111`; asserted `:466-472`.
- **axis C** — closure `{C-SELF@:166-170, C-SELF@:116-118, C-SELF@:109-111}`; `PHYSICS-FED = false`;
  `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

#### `QID:lambda-m@mathematica/…_audit.wl#227`
- **binding site** — `:227`, `lambdaRef = lambdas["20"]`; separately consumed at `:229`, `:316`.
- **reachability witness** — `RES:016:K2-FORM`, hop 1: `:227` → `k2Ref = buildK2[lambdaRef]` (`:229`).
- **axis A** — `A-REDUCED`; reduction performed at `:183`; reduces to `:166-170`, `:116-118`, `:109-111`.
- **axis B** — `B-EXECUTED`; produced by executed code at `:183`.
- **axis C** — closure as `#180`; `PHYSICS-FED = false`; `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

---

### 3.10 `QID:Gram-l2` — the 5×5 orthonormality matrix

⚠ **Universe call, flagged (see §5 C4).** `REPORTED_RESULTS.md:111` demotes `Gram = I₅` from being a
*result*. It is admitted here as an **input** to the "**5-dimensional** irrep" clause of
`RES:016:L2-IRREP` — the linear independence of the five channels is what makes the shared eigenvalue a
5-dim degeneracy rather than a coincidence among possibly-dependent functions. Being demoted as a result
does not exclude a quantity from another result's closure.

#### `QID:Gram-l2@scripts/…_sympy_audit.py#281`
- **binding site** — `:281`, `gram = sp.Matrix([[integrate_s2(ys[i]*ys[j]) …]])`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1 **via the claim's "5-dimensional" clause**;
  asserted `:676-677`.
- **axis A** — `A-REDUCED`. Reduction performed at `:281`; reduces to `Y_m` (`:270-274`) and the measure
  (`:233-235`).
- **axis B** — `B-EXECUTED`. Code locus `:281`; leaves `:270-274`, `:233-235`; asserted `:676-677`.
- **axis C** — closure `{C-SELF@:270-274, C-SELF@:233-235}`; `PHYSICS-FED = false`;
  `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`. No convention claim is made or implied about `Gram` itself: it is
  asserted as a computed fact (`:676-677`). (The *normalisation* convention question sits on the harmonic
  rows, §3.1-3.5, where it is recorded.)
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `no-tier:unclassified-nonfed` / `none` / **no**.
  ⭐ The asymmetry against §3.9's `DELTA` is deliberate and follows `REPORTED_RESULTS.md:44-47`:
  orthonormality on S² is a property of the sphere and of the written basis and does **not** enter the
  medium's equation of motion, so no physical input is ever expected to feed it.
- **LIVE**.

#### `QID:Gram-l2@mathematica/…_audit.wl#178`
- **binding site** — `:178`, `localGram = Table[intS2[localYs[[i]] localYs[[j]]], …]`.
- **reachability witness** — `RES:016:L2-IRREP`, hop 1 via the "5-dimensional" clause; asserted `:462-463`.
- **axis A** — `A-REDUCED`; reduction at `:178`; reduces to `:166-170`, `:109-111`.
- **axis B** — `B-EXECUTED`; code locus `:178`; leaves `:166-170`, `:109-111`; asserted `:462-463`.
- **axis C** — closure `{C-SELF@:166-170, C-SELF@:109-111}`; `PHYSICS-FED = false`;
  `SELF-REFERENTIAL = false`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unclassified-nonfed`; §4 bucket **executed-but-not-physics-fed**.
- **should_be_tier / should_be_basis / delta** — `no-tier:unclassified-nonfed` / `none` / **no**.
- **LIVE**.

---

### 3.11 `QID:M-tilde` — `M̃`, the ℓ=2 reduced mass scalar

#### `QID:M-tilde@scripts/…_sympy_audit.py#219`
- **binding site** — `:219`, `Mtilde = sp.Symbol("Mtilde", positive=True, real=True)`. A declaration
  that assigns **no value** (B-D3).
- **reachability witness** — `RES:016:K2-FORM`, hop 1: `:219` → `m2_core = Mtilde` (`:507`), i.e. the
  result's `M₂ = M̃`.
- **axis A** — `A-REDUCIBLE-UNDERIVED`. §9 evidence — **the named route and where it is recorded**:
  `M̃ = ∫μ_η β₂² dV`, recorded as edge **R35** at `notes/parameter_register.md:185` (the
  `M̃, K̃, T̃_Ω` row) and `:302` (the R35 edge row), and stated in the stage note at `:212` (H1) and
  `:76-80` (§1.4's display `M₂ = μ_η β₂² dV`). **Not executed**: H1 (note `:210-213`) states outright
  that "Neither engine anywhere writes `M̃=∫μ_ηβ₂²dV`; it lives only in print strings", and
  stage043's manifest records the same fact independently as `moment_integral_executed=False`
  (`scripts/ledger_stage043_irreducible_count_range_sympy_audit.py:382-393`, `shorthand="R35"`).
  ⛔ `A-CALIBRATED` was tested and **fails §9**: the note calls `M̃` a "frozen calibration input"
  (`:311-314`) and stage043 carries `frozen_calibration_input=True` (`:390`), but **no benchmark is
  named at any locus** — the register's own row classes it `DERIVED`, not `CALIB`. A calibration claim
  without a benchmark is not `A-CALIBRATED`.
- **axis B** — `B-UNADJUDICATED`. ⚠ §9.0 default, taken because **§3.2 has no value for this case** —
  see §5 C2. Nothing computes `M̃` and nothing types a value for it; `B-DECLARED-LITERAL` would assert a
  typed value that does not exist, and `B-DERIVED-IN-FORM-UNEXECUTED` would assert a symbolic derivation
  that this artifact does not contain (H1).
- **axis C** — closure `{C-FREE@:219}`: the artifact never assigns it, and no value it stands for is
  fixed at a cited locus (note `:326` cites "011/012/013" as provenance, which is a stage list, not a
  locus with a value, so `C-PEER`'s citation condition fails — §3.3). No leaf establishes
  physics-feeding ⇒ `PHYSICS-FED = C-UNRESOLVED`. `SELF-REFERENTIAL` = **undefined** (§3.3).
- **CONVENTION-LADEN** — `false` (not a candidate convention: no transformation group is claimed for `M̃`).
- **is_tier** — `no-tier:unadjudicated` (§5.7, `PHYSICS-FED = C-UNRESOLVED` trigger). ⚠ **§5 C1
  conflict:** under §5's "projection of axis A" framing this row is `tier1-debt`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / **`named-route`** / **DELTA**. Route
  recorded at `notes/parameter_register.md:185`, `:302` (R35).
- **LIVE**.

#### `QID:M-tilde@mathematica/…_audit.wl#21`
- **binding site** — `:21`, `$Assumptions = … && Mtilde > 0 && …`. The only site at which WL declares the
  symbol; it asserts positivity, not a value. Used at `:228`.
- **reachability witness** — `RES:016:K2-FORM`, hop 1: `:21` → `m2Core = Mtilde` (`:228`).
- **axis A** — `A-REDUCIBLE-UNDERIVED`, same route and same recording loci as the PY twin; not executed
  in this engine either (H1 covers both engines).
- **axis B** — `B-UNADJUDICATED` (§5 C2).
- **axis C** — closure `{C-FREE@:21}`; `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (§5 C1 alternative: `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

### 3.12 `QID:K-tilde` — `K̃`, the ℓ=2 radial stiffness scalar

- `QID:K-tilde@scripts/…_sympy_audit.py#220` — `Ktilde = sp.Symbol("Ktilde", positive=True, real=True)`.
- `QID:K-tilde@mathematica/…_audit.wl#21` — `$Assumptions … Ktilde > 0 …`; used at `:216`.

Both occurrences carry the field values of §3.11 with one substitution: the named route is
`K̃ = ∫[T_w β₂'² + K_η β₂²] dV`, recorded at `notes/parameter_register.md:185` and `:302` (R35), and
obtainable by matching the note's two displays at `:57-58` (`K₂ = K̃ + λ_m·T̃_Ω`) and `:79`
(`K₂ = (T_w β₂'² + K_η β₂² + λ·T_Ω β₂²) dV`) — the matching step is mine and is stated so it can be
contested. Not executed (H1, note `:210-213`; stage043 `:390` `moment_integral_executed=False`).
Reachability witness: `RES:016:K2-FORM` hop 1 via `build_K2` (PY `:305`, WL `:216`).
⇒ axis A `A-REDUCIBLE-UNDERIVED`; axis B `B-UNADJUDICATED`; closure `{C-FREE}` ⇒
`PHYSICS-FED = C-UNRESOLVED`, `SELF-REFERENTIAL` undefined; `CONVENTION-LADEN false`;
`is_tier no-tier:unadjudicated` (alt. `tier1-debt`); `should_be_tier tier3-emergent` / `named-route` /
**DELTA**; **LIVE**.

### 3.13 `QID:T-Omega-tilde` — `T̃_Ω`, the ℓ=2 reduced angular-stiffness scalar

- `QID:T-Omega-tilde@scripts/…_sympy_audit.py#221` — `TomegaTilde = sp.Symbol("TomegaTilde", positive=True, real=True)`.
- `QID:T-Omega-tilde@mathematica/…_audit.wl#21` — `$Assumptions … TomegaTilde > 0`; used at `:216`.

Fields as §3.12, with the named route `T̃_Ω = ∫T_Ω β₂² dV` (`notes/parameter_register.md:185`, `:302`).
Not executed (H1). ⛔ `A-CALIBRATED` tested and **fails**: `parameter_register.md:182` classes the
**density** `T_Ω` as `CALIB` — a **different quantity** from `T̃_Ω`, and in any case no benchmark number
is named there either ("independent calibration input"). ⚠ Note `:180-182` records that `T_Ω`'s
independence is **UNDETERMINED** (an isotropic wall would give `T_Ω = T_w/a²`); that is recorded, not
resolved, and it does not supply a benchmark.
⇒ axis A `A-REDUCIBLE-UNDERIVED`; axis B `B-UNADJUDICATED`; `PHYSICS-FED = C-UNRESOLVED`;
`CONVENTION-LADEN false`; `is_tier no-tier:unadjudicated` (alt. `tier1-debt`);
`should_be_tier tier3-emergent` / `named-route` / **DELTA**; **LIVE**.

⚠ **Identity hazard carried, not resolved (§8.3 — the builder does not adjudicate merges):** note
`:187-192` records that `T̃_Ω`(016) and stage023's `T_Ω` carry the same dimension and the same
`ℓ(ℓ+1)`-shaped reduction one ℓ apart, and are "different quantity until R42 exists". This census keeps
them as **two** QIDs and flags the merge as an open identity question for the physics review leg.

---

### 3.14 `QID:K2` — the ℓ=2 angular stiffness (6 occurrences)

#### `QID:K2@scripts/…_sympy_audit.py#305` — the FORM
- **binding site** — `:305`, `return compact(Ktilde + coeff*TomegaTilde)` inside `build_K2` (`:304-305`).
  This is the occurrence at which the artifact asserts the **form** the result claims.
- **reachability witness** — `RES:016:K2-FORM`, **hop 0**.
- **axis A** — `A-REDUCIBLE-UNDERIVED`. Named route: the ℓ=2 reduction of the wall energy integral
  `K₂ = ∫[T_w β₂'² + (K_η + λ·T_Ω) β₂²] a²dw dΩ`, recorded at note `:76-80` and, in that exact form, at
  `notes/parameter_register.md:182` (the `T_Ω` row). Not executed: the integrand is assembled at `:337-340`
  **for the dimension walk only** and never evaluated, and note `:259-270` (H8) records that this
  assembly is a *parallel reconstruction* whose terms can be dropped without failing anything.
- **axis B** — `B-DECLARED-LITERAL`. The additive form is typed at `:305`.
- **axis C** — closure `{C-FREE@:220 (Ktilde), C-FREE@:221 (TomegaTilde), C-FREE (the unbound `coeff`
  parameter at the definition site)}`. No leaf establishes physics-feeding ⇒
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false` (not a candidate convention).
- **is_tier** — `no-tier:unadjudicated` (§5.7 C-UNRESOLVED trigger; §5 C1 alternative `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**
  (route at note `:76-80`, `parameter_register.md:182`, R35 at `:302`).
- **LIVE**.

#### `QID:K2@scripts/…_sympy_audit.py#499` — assembled per channel
- **binding site** — `:499`, `k2_core = {name: build_K2(lambdas[name]) for name in order}`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0 (this is the assembly the result's "uses the live
  computed λ" clause is about; asserted via `:500-504`, `:693-697`).
- **axis A** — `A-REDUCIBLE-UNDERIVED`, same route/recording loci as `#305`.
- **axis B** — `B-EXECUTED`. Code locus `:499`; input leaves `:220`, `:221` and the live `lambdas`
  computed at `:287` (whose own leaves are `:270-274`, `:241-242`, `:233-235`); asserted `:693-697`.
- **axis C** — closure `{C-FREE@:220, C-FREE@:221, C-SELF@:270-274, C-SELF@:241-242, C-SELF@:233-235}`.
  Contains a `C-FREE` leaf and **no** leaf establishing physics-feeding ⇒ `PHYSICS-FED = C-UNRESOLVED`
  (§3.3, explicitly *not* `false`); `SELF-REFERENTIAL` undefined.
  ⚠ Consequence per §3.3: this row is ⛔ **not** admitted to §4's `executed-but-not-physics-fed`, which
  is a positive claim, despite being `B-EXECUTED`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

#### `QID:K2@scripts/…_sympy_audit.py#508` — the reference assembly
- **binding site** — `:508`, `k2_ref = build_K2(lambda_ref)`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0; a second binding site of the result's quantity
  (§7.2). ⚠ Its downstream consumer (`:509`, the dimensional block) is out of universe; the row is in
  universe because the **quantity** is (see §5 C5).
- **axis A** — `A-REDUCIBLE-UNDERIVED`, same route.
- **axis B** — `B-EXECUTED`; code locus `:508`; leaves `:220`, `:221`, `:506`→`:287`.
- **axis C** — closure as `#499`; `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

#### `QID:K2@mathematica/…_audit.wl#216` — the FORM
- **binding site** — `:216`, `buildK2[coeff_] := clean[Ktilde + coeff TomegaTilde]`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0.
- **axis A** — `A-REDUCIBLE-UNDERIVED`, route as PY `#305`; not executed in this engine (the integrand at
  `:231-236` is dimension-walked only, `:259-264`).
- **axis B** — `B-DECLARED-LITERAL` (form typed at `:216`).
- **axis C** — closure `{C-FREE@:21 (Ktilde), C-FREE@:21 (TomegaTilde), C-FREE (unbound `coeff`)}`;
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

#### `QID:K2@mathematica/…_audit.wl#219` — assembled per channel
- **binding site** — `:219`, `k2Core = AssociationThread[order, buildK2[lambdas[#]] & /@ order]`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0; asserted `:481-482`.
- **axis A** — `A-REDUCIBLE-UNDERIVED`, same route.
- **axis B** — `B-EXECUTED`; code locus `:219`; leaves `:21`, `:21`, and live `lambdas` from `:183`
  (own leaves `:166-170`, `:116-118`, `:109-111`); asserted `:481-482`.
- **axis C** — closure `{C-FREE@:21, C-FREE@:21, C-SELF@:166-170, C-SELF@:116-118, C-SELF@:109-111}`;
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined; not admitted to
  `executed-but-not-physics-fed`.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

#### `QID:K2@mathematica/…_audit.wl#229` — the reference assembly
- **binding site** — `:229`, `k2Ref = buildK2[lambdaRef]`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0 (second binding site; consumer `:316` is out of
  universe — §5 C5).
- **axis A** — `A-REDUCIBLE-UNDERIVED`, same route.
- **axis B** — `B-EXECUTED`; code locus `:229`; leaves `:21`, `:21`, `:227`→`:183`.
- **axis C** — closure as `#219`; `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

---

### 3.15 `QID:M2` — the ℓ=2 reduced mass

#### `QID:M2@scripts/…_sympy_audit.py#507`
- **binding site** — `:507`, `m2_core = Mtilde`.
- **reachability witness** — `RES:016:K2-FORM`, **hop 0** (the result's `M₂ = M̃` half).
- **axis A** — `A-REDUCIBLE-UNDERIVED`. Named route: `M₂ = ∫μ_η β₂² dV` (note `:76-78`;
  `parameter_register.md:185`, `:302` R35). Not executed (H1, note `:210-213`).
- **axis B** — `B-DECLARED-LITERAL` (a typed identification with the free symbol `Mtilde`).
- **axis C** — closure `{C-FREE@:219}`; `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

#### `QID:M2@mathematica/…_audit.wl#228`
- **binding site** — `:228`, `m2Core = Mtilde`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0.
- **axis A** — `A-REDUCIBLE-UNDERIVED`, same route/recording loci.
- **axis B** — `B-DECLARED-LITERAL`.
- **axis C** — closure `{C-FREE@:21}`; `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated` (alt. `tier1-debt`).
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `named-route` / **DELTA**.
- **LIVE**.

---

### 3.16 `QID:K2-TOmegaTilde-coefficient` — the coefficient extracted back out of `K₂`

#### `QID:K2-TOmegaTilde-coefficient@scripts/…_sympy_audit.py#500`
- **binding site** — `:500`, `k2_coeff = {name: extract_k2_coeff(k2_core[name]) …}`; the extractor
  `sp.diff(k2_expr, TomegaTilde)` is at `:309`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0: it is the object of the result's clause "the
  `T̃_Ω` coefficient **being** the computed eigenvalue", asserted through the residual at `:502` and
  `:693-697`.
- **axis A** — `A-REDUCED`. Reduction performed at `:500` (extractor `:309`); reduces to `k2_core`
  (`:499`) and through it to `λ_m` (`:287`), `Y_m` (`:270-274`), the operator (`:241-242`) and the
  measure (`:233-235`).
- **axis B** — `B-EXECUTED`. Code locus `:500`; leaves `:499`, `:220`, `:221`; asserted `:693-697`.
- **axis C** — closure `{C-FREE@:220, C-FREE@:221, C-SELF@:270-274, C-SELF@:241-242, C-SELF@:233-235}`
  under the **syntactic** walk (the expression differentiated at `:309` contains `Ktilde` even though the
  derivative eliminates it — see §5 C6). ⇒ `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined.
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` /
  **DELTA** (same basis as §3.9 — it is `λ_m` by construction).
- **LIVE**. ⚠ See §5 C7: under the *purpose* reading of `SELF-REFERENTIAL` this row is
  `SELF-REFERENTIAL = true` (the artifact hands `λ` in at `:499` and reads it back out at `:500`), which
  §3.3 would then bar from `A-REDUCED`. Both readings leave the row outside every tier; the alternative
  is counted in §4.

#### `QID:K2-TOmegaTilde-coefficient@mathematica/…_audit.wl#220`
- **binding site** — `:220`, `k2Coeff = AssociationThread[order, extractK2Coeff[k2Core[#]] & /@ order]`;
  extractor `D[k2Expr, TomegaTilde]` at `:217`.
- **reachability witness** — `RES:016:K2-FORM`, hop 0; asserted `:481-482`.
- **axis A** — `A-REDUCED`; reduction at `:220` (extractor `:217`); reduces to `:219` → `:183` →
  `:166-170`, `:116-118`, `:109-111`.
- **axis B** — `B-EXECUTED`; code locus `:220`; leaves `:219`, `:21`, `:21`; asserted `:481-482`.
- **axis C** — closure `{C-FREE@:21, C-FREE@:21, C-SELF@:166-170, C-SELF@:116-118, C-SELF@:109-111}`;
  `PHYSICS-FED = C-UNRESOLVED`; `SELF-REFERENTIAL` undefined (alt. reading: `true`, §5 C7).
- **CONVENTION-LADEN** — `false`.
- **is_tier** — `no-tier:unadjudicated`.
- **should_be_tier / should_be_basis / delta** — `tier3-emergent` / `physical-picture-expectation` / **DELTA**.
- **LIVE**.

---

## 4. The required output set (§5.8) — all 17, `is_tier` counts only (§6.1)

**Denominators (§10.1):** 38 occurrences / 16 QIDs, after the §6 exclusions. ⛔ No bare fraction is
quoted below; both levels are carried. The **rank** denominator (§10.1) is not computed here —
`OPEN-METHOD`.

| # | output | occurrence level | QID level |
|---|---|---:|---:|
| 1 | **TIER 1**, as a range (§5.6) | **[0, 30]** | **[0, 13]** |
|   | ├ `tier1-debt` | 0 | 0 |
|   | ├ `tier1-structural` | 0 | 0 |
|   | └ `tier1-postulate` | 0 | 0 |
| 2 | **TIER 2 / calibrated** (§10.2.1) | 0 | 0 |
| 3 | **TIER 3** | 0 | 0 |
|   | ├ split by axis B (§5.4): executed / unexecuted | 0 / 0 | — |
|   | └ split by propagation (§5.5): `tier3-calibration-propagated` / `tier3-held-out` | 0 / 0 | — |
| 4 | **`DERIVED`** (§4) | **0** | 0 |
| 5 | near-miss **`executed-but-not-physics-fed`** | **8** | 3 |
| 6 | near-miss **`derived-in-form-but-unexecuted`** | 0 | 0 |
| 7 | near-miss **`physics-fed-but-declared-literal`** | 0 | 0 |
| 8 | **`unclassified-nonfed`** | **8** | 3 |
| 9 | **convention** bucket | 0 | 0 |
| 10 | **`unadjudicated`** | **30** | 13 |
| 11 | **`convention-unadjudicated`** | **10** | 5 |
| 12 | **`self-referential`** | **0** | 0 |
| 13 | **conflict set** (§10.3) | — | **0** |
| 14 | **stage043 delta**, both directions | see §7 | see §7 |
| 15 | **out-of-scope list**, each row with its reason | see §6 | — |
| 16 | **reached-by-no-reported-result** | see §6 | — |
| 17 | **reported-result set** with loci | §1 (2 results) | — |

⛔ Do not sum this table (§5.8): buckets 5-7 overlap the tiers, and 8/10 are disjoint here only by
accident of this stage's distribution.

**Tier-1 lower bound = 0.** Every in-universe row lands outside the tiers, so no tier-1 assignment is
*established from evidence* and the whole tier-1 account sits in the upper bound. **The span is the
entire denominator.**

**Alternative counts, from the two flagged schema conflicts** (§5 C1, C7) — reported so neither is
invisible, and ⛔ neither replaces the primary:

- **C1 alternative** (axis A governs the projection where axis A is evidenced; `C-UNRESOLVED` guards
  tier 3 only): `tier1-debt` = **16** occurrences / **6** QIDs, `unadjudicated` = 14 / 7,
  TIER 1 = [16, 30] occurrences / [6, 13] QIDs. Tier 3, tier 2 and `DERIVED` are **0** in both readings.
- **C7 alternative** (`SELF-REFERENTIAL` read by purpose rather than by letter): `self-referential` = 2
  occurrences / 1 QID; those two rows lose `A-REDUCED` and are already `no-tier:unadjudicated`, so no
  tier count moves.

**Per-QID aggregation (§10.2):** no QID has a tiered occurrence, so §10.2.2 rule 2 does not fire —
`mixed-adjudication` = 0. All 16 QIDs are wholly tier-less; each QID's occurrences agree **in kind**
across the two engines, so the conflict set is **empty**. ⭐ That agreement is itself a reading: the two
engines are transliterations of one another, so cross-engine agreement carries no independent
provenance information (the stage's own §1.5b/H9 says the same about its dimension records).

---

## 5. Where the schema could not decide — blocked rows and missing rules

⭐ Every item below was hit while classifying. None was resolved by picking the better-reading answer.

**C1 — ⛔ BLOCKING for 22 of 38 rows: `C-UNRESOLVED` overrides an evidenced tier-1 axis-A value, and the
schema says both that it does and that it does not.**
- §3.3 ("`C-UNRESOLVED` rows go to `unadjudicated` (§5.3)") and §5.7 (`no-tier:unadjudicated` ←
  "`PHYSICS-FED = C-UNRESOLVED`") say the tier is `no-tier:unadjudicated`.
- §5's own framing says the tiers are "a projection of axis A (**with one stated axis-C guard on tier 3**,
  §5.1)", and §5/§5.7 say `tier1-debt` ← `A-REDUCIBLE-UNDERIVED`. §9.0 adds that an axis-C failure
  "does **NOT** touch axis A".
- **Missing rule:** *precedence between an evidenced tier-1 axis-A value and a `C-UNRESOLVED` closure.*
  The stakes are the whole tier-1 **lower** bound: 16 occurrences / 6 QIDs move between [0,·] and [16,·].
- Rows affected: §3.11-§3.16 (all `M̃`/`K̃`/`T̃_Ω`/`K₂`/`M₂`/coefficient occurrences), 22 rows.
- Primary reading taken: §3.3 + §5.7 (two explicit rules beat one framing sentence). Alternative counted
  in §4.

**C2 — axis B has no value for a declared free symbol that carries no value.** `M̃`, `K̃`, `T̃_Ω` are
declared (PY `:219-221`, WL `:21`) and never assigned. `B-DECLARED-LITERAL` asserts a typed value that
does not exist; `B-EXECUTED` and `B-DERIVED-IN-FORM-UNEXECUTED` assert a computation this artifact does
not contain; `B-ASSERTED-TARGET` is wrong (nothing is checked against them here). Recorded
`B-UNADJUDICATED` under §9.0 on 6 rows, which mislabels a **positively established absence** as a
failure to establish. **Missing value:** something like `B-NO-VALUE-DECLARED`.

**C3 — a leaf can be both `C-SELF` and `C-MATH`, and §3.3 gives no precedence.** The S² measure
(PY `:233-235`) and the `−Δ_S²` operator (PY `:241-242`) are literals declared in the artifact **and**
pure-mathematical objects. Rule adopted and applied uniformly: **tag `C-SELF` where the artifact declares
it literally, and record the `C-MATH` reading.** ⚠ Not outcome-bearing here — neither tag confers
`PHYSICS-FED` — but it would be if a census ever wanted to separate "the artifact made this up" from
"this is mathematics".

**C4 — §3.3's closure walk assumes an expression; a reported result stated as a prose claim has clauses,
and the schema does not say which binding sites its clauses reach.** `RES:016:L2-IRREP` claims a
"**5-dimensional** SO(3) irrep". `λ_m`'s value does not depend on `Gram`, but the 5-dimensionality
clause does. Decision taken: `Gram` is **in** universe as an input to that clause (§3.10), flagged.
Excluding it would remove 2 rows (both `no-tier:unclassified-nonfed`) and would not move any tier count.

**C5 — §7.1.1 defines universe membership for *quantities*; §7.2 defines occurrences for *binding
sites*; nothing says whether a binding site of an in-universe quantity that is itself in no closure is
in universe.** PY `:508` / WL `:229` (`k2_ref`) bind `K₂` but feed only the out-of-universe dimensional
block. Decision: **in** (a quantity's binding sites are its occurrences), which is the direction that
keeps rows in the account. Affects 2 rows.

**C6 — the closure walk has no rule for inputs the producing operation eliminates.** `∂(K̃ + λT̃_Ω)/∂T̃_Ω`
(PY `:500`) syntactically consumes `K̃` — a `C-FREE` leaf — while the value depends only on `λ`. The
syntactic walk gives `C-UNRESOLVED`; a value-dependence walk gives `PHYSICS-FED = false`. Primary
reading: syntactic (§3.3 says "walk back **the expression** that produced its value"). Both leave the row
outside every tier; under the value-dependence reading the 2 rows move from `unadjudicated` to
`unclassified-nonfed`.

**C7 — `SELF-REFERENTIAL` is defined as reaching *the occurrence's own declaration*, which misses a
round trip through a different occurrence.** PY `:499` puts `λ` into `K₂`; PY `:500` differentiates it
back out; the pair is asserted at `:502`. By the letter, `SELF-REFERENTIAL = false`; by §3.3's stated
purpose ("the closure returned the value it was handed") it is `true`. The stage's own de-count note
(`:68-71`) treats exactly this shape as vacuous when done one step more directly. Letter taken;
alternative counted in §4.

**C8 — one occurrence carrying two competing substrate claims has no adjudication rule.** For `M̃`/`K̃`/`T̃_Ω`
the note says "frozen calibration inputs" (`:311-314`) and stage043 carries
`frozen_calibration_input=True` (`:390`), while `parameter_register.md:185` classes the same quantities
`DERIVED` with route R35. §10.3's conflict machinery is defined over a QID's **occurrences**, not over
two claims about **one** occurrence. Resolved here by §9 evidence (no benchmark ⇒ `A-CALIBRATED` fails;
route recorded ⇒ `A-REDUCIBLE-UNDERIVED` holds) rather than by a schema rule; recorded so it is visible.

**⛔ Rows blocked outright: 0.** Every row carries a value on all three axes, the flag, both tiers and a
basis. C1 and C7 are recorded as *contested* `is_tier` values with both counts reported, not as blanks.

### Attribution defects found by opening the loci (§9.1 rule 2)

**D4 — two documents point at the wrong lines of the PY engine, both by exactly 41 lines.**
- `notes/parameter_register.md:182`, `:183`, `:185` each attribute the stage016 dimension declarations to
  `scripts/ledger_stage016_l2_so3_covariance_sympy_audit.py:355-366`. **Opened:** `:355-366` is the
  `ok = bool(…)` conjunction and the head of `dimension_eval`'s return dict. The declarations are in
  `make_dim_rules` at `:312-326` (the twelve entries at `:314-325`) — exactly **41 lines earlier**.
- The stage note repeats the same stale range at `:194` (`sympy:355-366` ↔ `.wl:239-250`) and at `:212`
  (`sympy:364` for the `M̃`/`K̃`/`T̃_Ω` declarations; `:364` is `"k2_terms": {`, and `364 − 41 = 323` is
  `Mtilde: EXPECTED_M`). The **`.wl` half of both citations is correct** (`makeDimRules` entries at
  `.wl:239-250`), and the note's `sympy:723` is correct.
- ⇒ The values agree; the pointers do not — the same defect shape §9.1 was written from, at a different
  pair of files. ⛔ Not repaired here (this census does not edit substrate); recorded with both loci.

---

## 6. Out of scope (§7.3), with reasons and counts

⛔ Nothing was silently dropped. Exclusions split into **individually enumerated** rows and
**class-level** records; the second kind is marked as such rather than being presented as a row count.

### 6.1 Individually enumerated — reason `reached-by-no-reported-result` (§5.8 output 16)

**64 binding sites.** Neither reported result's closure reaches any of them: the dimensional block
**consumes** `K₂` and `M₂` (PY `:348-349`, WL `:265-266`) and produces nothing either result's value
depends on.

| set | PY loci | WL loci | count |
|---|---|---|---:|
| the 12 `dim_rules.*` exponent-triple declarations | `:314-325` | `:239-250` | 24 |
| the 9 `baseline_dims.*` walked results | `:342-350` (emitted `:478-488`) | `:259-267` (emitted `:509-517`) | 18 |
| dimension helper/target constants (`ZERO_DIM`, `EXPECTED_M`, `EXPECTED_K`, `EXPECTED_RATIO`, the inline `L³` target) | `:160-163`, `:353` | `:23-26`, `:270` | 10 |
| the 6 dimension-walk integrand expressions (`measure`, `m2_integral`, the 3 `K₂` terms, `k2_integral`) | `:335-340` | `:231-236` | 12 |

⭐ Those first two rows are **21 records per engine** — the row universe of the stage's `21 of 21 CORRECT`
dimensional verdict. **20 of the 21 are out of this census's universe**; the one overlap is not a record
at all but the quantities `M̃`/`K̃`/`T̃_Ω`/`K₂`/`M₂` that the records *describe*, which enter this census
through `RES:016:K2-FORM` and not through the dimension machinery. ⇒ The two universes are near-disjoint
by construction, which is why no resemblance between the two distributions should be expected (§11).

### 6.2 Recorded at class level, not per binding site

Each class carries its §7.3 reason and its loci; ⛔ the individual binding sites inside them were **not**
enumerated, and no count is asserted for them beyond "order 10² across the two engines".

| class | §7.3 reason | loci |
|---|---|---|
| the three able-to-fail probes, their self-ablations, the 4 corrupted rule maps and the aggregate battery | controls and deliberate negatives | PY `:406-422`, `:431-458`, `:516-627`; WL `:300-313`, `:318-340`, `:348-439` |
| the per-tooth ablations on copies | controls and deliberate negatives | PY `:770-816`; WL `:576-612` |
| assertion residuals and hash guards (`residuals`, `k2_coeff_residuals`, `self_overlaps`, `input_hashes`) | test scaffolding — outputs of checks, never inputs to a result | PY `:290`, `:502`, `:511`, `:513`; WL `:187-193`, `:221-224`, `:343`, `:345` |
| gate booleans and verdict tokens (`gram_is_identity`, `lambda_all_six`, `residuals_zero`, `covariant_ok`, `gates`, `verdict`, `ISOTROPY_CALIBRATED`…) | gate pass/fail flags (§7.1.1 exclusion list) | PY `:30-34`, `:295-300`, `:504`, `:512-514`, `:628-640`; WL `:15-19`, `:198-203`, `:225`, `:344-346`, `:441-448` |
| run tallies and print/label strings | display and formatting constants; row counts | PY `:27-28`, `:819-867`, `:894-906`; WL `:12-13`, `:630-669`, `:694-702` |
| the WL arity self-check's transient re-run of `evalDimensional` | test scaffolding (and a duplicate of an already out-of-scope object) | WL `:614-628` |

**RETIRED rows (§7.4): 0.** No occurrence in either engine is marked retired, and neither engine carries
a retired-approach value. The stage's own de-counted construct — `k_coeff_equal` (note `:68-71`) —
**has no binding in either engine**: verified by grep, the name survives only inside print strings
(PY `:689`, `:816`, `:839`; WL `:478`, `:611`, `:650`), so it produces no row. ⇒ The note's claim that
it "survives only as a documentation string, not a counted check" (`:71`) is **confirmed by
measurement**, not inherited.

---

## 7. Reconciliation with stage043 (§8.4) — the stage016 slice only

⛔ The 152-ID manifest was **not** re-derived; it was read at
`scripts/ledger_stage043_irreducible_count_range_sympy_audit.py`.

**Direction 1 — `REG:` IDs mapping to ≥ 1 census occurrence: 7.**

| `REG:` ID | stage043 `ratified_category` (locus) | census QID |
|---|---|---|
| `REG:C1:Mtilde` | `extension-convention-open`, `source_status=DERIVED-IN-FORM-UNEXECUTED`, `moment_integral_executed=False` (`:382-393`) | `QID:M-tilde` |
| `REG:C1:Ktilde` | same (`:382-393`) | `QID:K-tilde` |
| `REG:C1:Ttilde_Omega` | same (`:382-393`) | `QID:T-Omega-tilde` |
| `REG:derived:Mtilde_definition` | `structural-no-knob` (`:288-311`) | `QID:M-tilde` |
| `REG:derived:Ktilde_definition` | `structural-no-knob` (`:288-311`) | `QID:K-tilde` |
| `REG:derived:Ttilde_definition` | `structural-no-knob` (`:288-311`) | `QID:T-Omega-tilde` |
| `REG:derived:lambda_m_SO3` | `structural-no-knob` (`:288-311`) | `QID:lambda-m` |

⚠ **Two `REG:` IDs map to one census QID** in three cases (the C1 row and the definition row are the same
quantity under two register roles). §8.4 does not say whether that is a defect; recorded, not resolved.
⚠ `REG:derived:lambda_m_SO3` names edge **R34**, whose text (`parameter_register.md:301`) bundles
`Gram=I₅`, `λ_m=6` **and** the `K₂` form under one ID. `QID:Gram-l2` and `QID:K2` are therefore
**candidate** mappings to that same ID — ⛔ recorded as an **open identity question**, not as a merge
(§8.3: the builder does not adjudicate merges).

**Direction 2 — census QIDs mapping to no `REG:` ID: 12 of 16** — `QID:Y20`, `QID:Y21c`, `QID:Y21s`,
`QID:Y22c`, `QID:Y22s`, `QID:S2-measure`, `QID:neg-Delta-S2`, `QID:neg-Delta-S2-applied`, `QID:K2`,
`QID:M2`, `QID:K2-TOmegaTilde-coefficient`, `QID:Gram-l2`. This is the census's genuine extension for
this stage: **the register has no entry for the basis, the measure, the operator, or the assembled
`K₂`/`M₂` themselves.**

**Not reconciled: 145 of 152 `REG:` IDs.** They belong to stages outside this pilot and were ⛔ **not**
recorded as out-of-scope, because §8.4's option 2 requires a reason referred to §7.1.1 or §7.3, and that
reason is not decidable from stage016 (§7.1.1: membership needs one corpus-wide pass). **This is an open
obligation of the census, not a discharged one.**

---

## 8. Coverage statement

**Covered in full, row by row:** the transitive input closure of both stage016 reported physical results,
across both engines — **16 quantities, 38 occurrences**, every one carrying all three axes, the flag,
the closure with leaf tags, `PHYSICS-FED`, `SELF-REFERENTIAL`, `is_tier`, `should_be_tier`,
`should_be_basis`, `delta`, LIVE/RETIRED, IDs and §9 evidence. **No row was truncated and none is
blocked.**

**Selection rule where I did not enumerate:** §6.2's six exclusion classes are recorded **as classes with
their loci and reasons**, not as individual out-of-scope rows. The rule is: *a binding site excluded by a
§7.3 class other than `reached-by-no-reported-result` is recorded at class level; every
`reached-by-no-reported-result` exclusion is enumerated individually.* That is why §6.1 lists 64 exact
loci and §6.2 asserts no count.

**Not covered, with counts:**

1. **The stage note as an artifact — 0 of its binding sites censused.** `notes/stages/ledger_stage016_l2_so3_covariance.md`
   is an artifact under §7.2 and its markdown rows are binding sites; the two engines' rows above
   therefore under-count every quantity the note also asserts. At minimum the note's `:57-58` (`K₂`/`M₂`),
   `:76-80` (the density relations) and `:29-35` (the five harmonics) are uncensused binding sites of
   in-universe quantities — **≥ 8 occurrences not emitted**. Scoped out by B-D4 (the assignment named the
   two engines).
2. **145 of 152 stage043 `REG:` IDs unreconciled** (§7) — an open §8.4 obligation.
3. **`QID_REGISTRY.md` not written.** §8.3 requires one; this pilot mints the 16 QIDs above but does not
   create the registry file, and no merge is adjudicated (§8.3 reserves that for the physics review leg).
4. **The rank denominator (§10.1) not computed** — `OPEN-METHOD`, and out of a single stage's reach.
5. **Cross-stage duplication not resolved.** `REPORTED_RESULTS.md:298-301` records that stage017 rides
   this covariance; if 017 presents it as *its* result, some rows above may belong to a joint closure.
   Nothing here presumes an answer.
