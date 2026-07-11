# ledger_stage029 — the Post-Newtonian corpus DOI-cite (Check II-P6)

**Part / anchor.** Part II — Gravity (Cluster C, the CITE-ONLY provenance cap). The audited external post-Newtonian
ladder (1PN→4PN + 2.5PN) — seven separately-published, Zenodo-DOI'd papers under `research/` — recorded here by exact
DOI so the by-DOI imports in Stages 020 / 027 / 028 have visible in-Part provenance. Source: the seven corpus dirs
`research/{1pn_orbital_dynamics, 4d_1pn_bridge, 4d_1pn_full, 4d_2pn, 4d_2_5pn, 4d_3pn, 4d_4pn}` + the top-level
`README.md` (the authoritative DOI list). **The LAST Part-II gravity stage** — after it the sector CLOSES and the
scheduled MIDWAY KNOB AUDIT runs.

**Landing.** `PN_CORPUS_CITED` — a **documentation landing**, NOT a verdict token. ⚠ **There is NO `FAIL_*`/`PASS`
physics verdict, NO computed object, and NO dual-engine audit** — this is a CITE-ONLY provenance stage. Per the
blueprint's completeness standard #3 the PN corpus is **cited by DOI, NOT inlined** (the `4d_*pn*` papers are
separately audited, DOI'd, and already in `research/`).

**Status.** Cite-only (no `scripts/`, no `mathematica/`, no `output/` transcripts). Adds **zero** new counted knobs and
**zero** new dims — the PN results are imported external facts, not ledger parameters. Coverage class: **No executable
audit** (the first such stage in the rebuilt ledger).

**⚠⚠ Scope invariant (load-bearing honesty): Stage 029 does NOT re-derive, re-audit, or re-verify any PN result.** The
PN physics audit lives in the seven corpus papers. 029 is provenance only. The ledger's OWN gravity earning is Stages
001–028; 029 records the external audited corpus those stages calibrate-predict against.

---

## 1. What this stage records — the seven-entry PN corpus

The rebuilt ledger's gravity sector is a **calibrate-predict** construction: the PN corpus is the audited, GR-matched
reference against which the ledger's earned form/branch/fingerprint results are checked. All seven titles match their
`README.md` DOI entry exactly (zero mismatches). Scope class for all seven = **CITE-only / external audited corpus**
(imported, not re-derived); two are CONDITIONAL (`4d_2_5pn`, `4d_4pn`), flagged as such.

| # | corpus dir | Zenodo DOI | paper title | PN order / imported headline result |
|---|---|---|---|---|
| 1 | `research/1pn_orbital_dynamics` | [`10.5281/zenodo.19449058`](https://doi.org/10.5281/zenodo.19449058) | *Newtonian and 1PN Orbital Dynamics from a Superfluid Defect Toy Model* | **0PN Newtonian + 1PN** perihelion precession; the Static-Limit Theorem; `β = κ_ρ+κ_add+κ_PV = 3` (`κ_ρ=1, κ_add=1/2, κ_PV=3/2`, with `κ_PV=3/2` fixed by the 1PN GR-matching condition, not derived in that paper) |
| 2 | `research/4d_1pn_bridge` | [`10.5281/zenodo.19449653`](https://doi.org/10.5281/zenodo.19449653) | *Deriving Key Post-Newtonian Coefficients from the Unified 4D Toy Model* | **1PN coefficient bridge**: `κ_add=1/2` (topology discriminator vs `1/3`), EOS exponent `n=5`, `α²=3/4`, `K_vec=2/π²`, universal `v⁴` coefficient `3/8` |
| 3 | `research/4d_1pn_full` | [`10.5281/zenodo.19450102`](https://doi.org/10.5281/zenodo.19450102) | *A Controlled Derivation of the Full First Post-Newtonian Sector from the Unified 4D Toy Model* | **Full conservative 1PN** two-body EIH Lagrangian; `β₁PN=3`, EIH cross-coefficients `C∥=−7/2, C_L=−1/2`; matches the standard EIH Lagrangian + perihelion shift |
| 4 | `research/4d_2pn` | [`10.5281/zenodo.19450284`](https://doi.org/10.5281/zenodo.19450284) | *A Controlled Derivation of the Full Conservative Second Post-Newtonian Sector from the Unified 4D Toy Model* | **Full conservative 2PN** (order `c⁻⁴`); Legendre-transforms EXACTLY to the standard generic-frame conservative ADM Hamiltonian through 2PN; fixes `λ_ρ=1/2` |
| 5 | `research/4d_2_5pn` | [`10.5281/zenodo.19492270`](https://doi.org/10.5281/zenodo.19492270) | *A Conditional Derivation of the Point-Particle 2.5PN Sector from the Unified 4D Toy Model* | **Conditional 2.5PN** dissipative gate; recovers the **Burke–Thorne / Iyer–Will** quadrupole radiation-reaction; ⭐ its single open normalization `m̂₀²Γ₅ = 2G/(5c⁵)` is the item Stage 028 (via pathA_43) MEETS at reduced-closure |
| 6 | `research/4d_3pn` | [`10.5281/zenodo.19501724`](https://doi.org/10.5281/zenodo.19501724) | *A Full Conservative Derivation of the Two-Body 3PN Sector from the Unified 4D Toy Model* | **Full conservative 3PN** (order `c⁻⁶`) in a fixed ADM chart; repair-ledger `(μ_{ρ,3}, d₃, s₂₄) = (1/4, −45/4, −1/16)`; full ledger assigned, `ΔH₃ = −ΔL₃` |
| 7 | `research/4d_4pn` | [`10.5281/zenodo.19561056`](https://doi.org/10.5281/zenodo.19561056) | *A Conditional Full Conservative Derivation of the Two-Body 4PN Sector from the Unified 4D Toy Model* | **Conditional 4PN** (local + tail); the tail coefficient `C_tail = (GM/2c³) γ_quad^eff` ties to the **same** 2.5PN quadrupole-normalization gap (`γ_quad^eff = 2G/(5c⁵) ⇒ C_tail = G²M/(5c⁸)`) |

**Shared foundational root.** The seven papers descend from a common ancestor, the Action paper *4D Toy Model — Action,
Projections, and Controlled Brane Limits* ([`10.5281/zenodo.19449589`](https://doi.org/10.5281/zenodo.19449589)); the
two 4D 1PN papers (`4d_1pn_bridge`, `4d_1pn_full`) cite it explicitly. It is the derivation-chain root of the 4D
toy-model series; it is not one of the seven PN-ladder dirs and so is not tabulated above, only named here as the
shared root.

## 2. Provenance caveat — the DOIs are README-authoritative, not source-cross-validated

The top-level `README.md` "Papers List" (16 DOI-bearing bullets) is the **sole authoritative source** for these DOIs.
A read of all seven `.tex` sources confirms that **none of the seven papers declares its own Zenodo DOI in source** —
the two 4D 1PN papers cite only the foundational Action DOI `19449589`, and the other five carry no Zenodo/DOI
citation. Consequently:
- there is **no self-declared DOI to conflict** with the README → zero DOI mismatches;
- the available cross-checks are (a) exact **title-match** against the README (all seven AGREE) and (b) both 4D 1PN
  papers correctly citing the shared `19449589` root.

The DOIs are Zenodo-registration facts, **not re-derivable from the paper sources**. This note cites them from the
README as authoritative; it does not represent them as source-verified.

## 3. The 2.5PN line-map — backing Stage 028's INV3 import

Stage 028's 2.5PN Burke–Thorne match-back imports the corpus's OWN 2.5PN form (its INV3) and the Burke–Thorne
normalization. In `research/4d_2_5pn/paper/4d_2_5pn.tex` (the display block spans L467–474):

| Line(s) | Content |
|---|---|
| **L469** | the moment invariant `K̄₄ K̄₀ = 4 K̄₂²` |
| **L471–473** | the corpus form `Γ̄₅ = 9 K̄₂^{5/2}/K̄₀^{3/2}` (Stage 028's INV3) |
| **L60**, **L822** (boxed) | the open normalization `m̂₀² Γ₅ = 2G/(5c⁵)` |
| **L797** | `γ_GR = 2G/(5c⁵)` |
| **L793**, **L961** | the Burke–Thorne force with coefficient `−2G/(5c⁵)` |

⚠ Stage 028's card + source map cite a blanket `4d_2_5pn.tex:469` for the corpus form; the corpus form is precisely at
L471–473 (L469 is the invariant). This note tightens that citation. It does **not** edit Stage 028 (committed).

## 4. Why Stage 029 exists — the derivation-chain dependency

The in-Part imports that lean on this corpus by DOI:
- **Stage 028** (the 2.5PN match-back) imports the corpus form `9K̄₂^{5/2}/K̄₀^{3/2}` (INV3) and `Γ̄₅=2G/(5c⁵)` from
  `4d_2_5pn` — the single open normalization item pathA_43 meets at reduced-closure.
- **Stage 020** (source gate pathA_33; `54/5`, `G=GENUINE_BLOCKED`) classifies `2/5·G` as `external_bridge_input` —
  the "external bridge" is this PN/GR corpus (the `2G/5c⁵` radiation normalization).
- **`4d_4pn`** ties its tail coefficient to the same 2.5PN quadrupole-normalization gap (corpus-internal consistency
  with the ledger).
- The 1PN→3PN conservative ladder is the GR-matched backbone of the gravity sector's calibrate-predict validation.

Per the settled recommendation (Part-II atomic split, "PN-cite as a thin Part-II stage rather than a register item"),
this dependency is given a visible in-Part stage rather than a bare Part-VII register pointer, so the provenance reads
where the gravity sector reads.

## 5. Honest CITE-only scope + physical framing

**Scope.** Stage 029 records the corpus; it does not re-derive or re-verify it. The `54/5`/`Γ̄₅`/`2G·5⁻¹c⁻⁵`
normalizations imported into 020/027/028 are `external_bridge_input` with `G=GENUINE_BLOCKED` — the PDE delivers the
form/branch, not Newton's `G` or the Burke–Thorne magnitude from first principles (the full 1PN→4PN from-throat
re-derivation remains SIM-DEFERRED, Gate 6).

**Physical framing (toy-analog).** In the model's picture, gravity is the throat drain (a one-way inflow whose
*change* propagates at the medium ripple speed `c_s`), and the gravity sector's post-Newtonian dynamics are validated
by this GR-matched corpus. Framed strictly as a toy analog: a superfluid-defect toy model reproducing GR-like PN
two-body dynamics through 4PN + 2.5PN radiation reaction — an external audited corpus, cited by DOI.

## 6. Consumed / exported + sector close

- **Consumes (external PROVENANCE):** the seven DOI'd corpus papers + the foundational Action paper `19449589`
  (shared root). Zero file I/O; zero re-derivation.
- **Exports:** the in-Part provenance that Stages 020 / 027 / 028 and the Part-VII calibration map cite.
- ⭐ **Stage 029 CLOSES Cluster C (024–029) AND the entire Part-II gravity sector.** Next: the scheduled
  **MIDWAY KNOB AUDIT** (`notes/parameter_register.md` §"MIDWAY KNOB AUDIT" — the pathA_40 `Δr=2` codimension dry-run
  over Parts I–II + the held-out vs irreducible-route-less tally + the (a) declared-universal-constants / (b)
  reduction-debt split). ⚠ Relevant audit datum: the PN results here are IMPORTED external facts, and the corpus's OWN
  internal coefficients (`κ_ρ, κ_add, κ_PV, λ_ρ, μ_{ρ,3}, …`) are the corpus papers' calibrations — they are NOT ledger
  knobs and must not be counted as such.

## 7. Registration

- **Manifest:** `canonical_stage_count` 28 → 29; new `policy` entry keyed `build_order_stage: "029"`, slug
  `pn_corpus_doi_cite`, CITE-only (no script pair).
- **Coverage:** `canonical`/`reviewed` → 29; `sympy`/`mathematica`/`verified` → 28 (⚠ this corrects a stage028
  registration lag — those 28 stages are dual-engine; at stage028 only `canonical`/`reviewed` had been bumped);
  `verified` stops at 28 (029 is documentary-reviewed, not executable-verified); a new "No executable audit | 1 | 029"
  coverage-class row; Part II → 24 total / 23 SymPy / 23 Mathematica / 0 numerical / 24 review; the provenance-index
  scope line → "001-029".
- **Parameter register:** zero new knobs, zero new dims; a provenance edge **R48** (the audited external PN corpus is
  cited by DOI — imported external facts, discharges nothing).
- **Verification (this stage):** doc-fidelity, not executable — the DOIs, titles, headline-result attributions, and
  line-map were cross-checked for document fidelity (no PN physics re-verified — that audit lives in the corpus papers)
  against the corpus sources + README via a Codex→Grok→Codex directive bookend and the post-authoring doc-fidelity
  tri-review (`FIDELITY_CLEAN` + `ADVERSARIAL_CLEAN`; register-verify `REGISTER_CLEAN` after 2 wording nits).
  There is no `scripts/`/`mathematica/` audit for Stage 029.
