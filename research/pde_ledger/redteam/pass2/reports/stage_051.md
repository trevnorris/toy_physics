---
unit_id: 051
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage051_lowest_twin_criterion.md]
  paper_appendix: present
---

# Audit unit 051 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_051.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage051_lowest_twin_criterion.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 80, 315, 348)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.txt`

## What the paper claims

Stage 051 is an `\StatusExactClosure{}` checkpoint (anchor MTDC-T6) that re-expresses the lowest-twin
support sufficiency test in the branch-product variables of the continuum placement map. The bottom-line
deliverable is the boxed criterion, `\stagefield{Output}{The lowest-twin branch-product criterion
\eqref{eq:app-stage051-twin-criterion}}`, namely `\boxed{\Pi_{\rm tr}\leq2C_{\rm mix}=\frac{16\Lambda(1-\epsilon)}{\pi^2}}`
with `C_{\rm mix}=8\Lambda(1-\epsilon)/\pi^2` and `\Pi_{\rm tr}=F_{\rm tr}G_{\rm tr}`. The card also states the
equivalent threshold scales `\Lambda_{\rm twin,req}=\frac{\pi^2}{16(1-\epsilon)}\Pi_{\rm tr}` and
`M_{\rm mix}^{({\rm twin,req})}=G_{\rm tr}/2`. The notes add three further exact deliverables: the closed product
form `Pi_tr` (notes §1), the two endpoint facts `Pi_tr(0)=0` and `Pi_tr(xi->1-)=+infinity`, the wall/mixed overlap
threshold `Z_W^(twin,req)=pi^2(1-eps_eta)(1-eps)G_tr/[16(1+chi0)^2]` (notes §3), and the exact twin-saturation
depth `xi_(2x)`, the positive quadratic root of `G_tr=2 M_mix` (notes §4). The appendix row 80 summarizes it as
"Branch-product criterion `\Pi_{\rm tr}\leq2C_{\rm mix}`."

## What the script claims to verify

The SymPy script (a) builds `Pi_tr` from the independent product `F_tr*G_tr` and asserts it equals the claimed
closed form (`expect_zero`, l.75); (b) checks the two endpoint limits (l.78–85); (c) builds `zeta_req` from the
Stage-048 map and asserts it vanishes at `Pi=C_mix`, equals 1 at `Pi=2 C_mix`, hence pins the threshold
equivalence `zeta_req<=1 <=> Pi<=2 C_mix` (l.99–107); (d) forms the three threshold scales and asserts the
`Z_W` threshold equals the Stage-047 forward map applied to `M_mix=G_tr/2` (l.130); (e) writes the closed
`xi_(2x)` root and asserts `G_tr(xi_(2x))=2 M_mix` (l.147). The Mathematica script mirrors the claim set but, for
`xi_(2x)`, independently `Solve`s the quadratic and compares the derived root to the closed-form claim before the
substitution check (l.100–111).

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `Pi_tr = F_tr G_tr` closed form (notes §1) | sympy l.75 / wl l.58 `Pi_tr - expected == 0` | match |
| `Pi_tr(0)=0`, `Pi_tr(xi->1-)=+inf` (notes §1) | sympy l.78–85 / wl l.60–69 | match |
| boxed `Pi_tr <= 2 C_mix = 16Λ(1-eps)/π²` (card eq. twin-criterion) | sympy l.99–107 / wl l.80–82 threshold equivalence at boundary | match |
| `Λ_twin,req = π² Π_tr/[16(1-eps)]` (card eq. Lambda-threshold) | sympy l.115 (printed form) / wl l.84 | match (printed, see F-note) |
| `M_mix^(twin,req)=G_tr/2` (card eq. Lambda-threshold) | sympy l.116 / wl l.85; used in Z_W map check l.130 | match |
| `Z_W^(twin,req)` (notes §3) | sympy l.130 / wl l.94 forward-map consistency | match |
| `xi_(2x)` root of `G_tr=2 M_mix` (notes §4) | sympy l.147 / wl l.100–111 (wl independently Solves) | match |

`paper_alignment: aligned`. Every paper/notes deliverable has a faithful, non-tautological script-side counterpart,
with values matching the card and notes exactly. The single finding is the stale committed transcripts.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 75 | `expect_zero(Pi_tr - Pi_expected)` | Pi_tr closed form | yes |
| A2 | sympy | 80–81 | `Pi0 == 0` | endpoint xi->0 | yes |
| A3 | sympy | 82–85 | `Pi1 is oo` | endpoint xi->1- | yes |
| A4 | sympy | 99 | `zeta_req(Pi=C_mix)==0` | threshold equivalence | yes |
| A5 | sympy | 100 | `zeta_req(Pi=2C_mix)-1==0` | boxed criterion boundary | yes |
| A6 | sympy | 104–107 | `(zeta_req-1)(Pi=2C_mix)==0` | same as A5 (redundant) | yes (dup) |
| A7 | sympy | 130 | `ZW_twin_req - map(M_mix=G_tr/2)==0` | Z_W threshold consistency | yes (mild, see note) |
| A8 | sympy | 147 | `G_tr(xi_2x) - 2 M_mix==0` | xi_(2x) root | yes |
| A9 | wl | 58 | `Pi_tr - expected (Factor/Together) == 0` | Pi_tr closed form | yes |
| A10 | wl | 64,68 | `pi0===0`, `1/pi1==0` | endpoints | yes |
| A11 | wl | 80–82 | three zeta_req boundary checks | boxed criterion boundary | yes |
| A12 | wl | 94 | `zWTwinReq - map == 0` | Z_W threshold consistency | yes |
| A13 | wl | 110 | `xi2xDerived (Solve) - xi2xClaim == 0` | xi_(2x) derived independently | yes (strong) |
| A14 | wl | 111 | `G_tr(xi2xDerived) - 2 mMix == 0` | xi_(2x) root | yes |

Notes on individual rows:
- A5/A6 are the same check written twice (`zeta_req.subs(Pi,2Cmix)-1` vs `factor((zeta_req-1).subs(Pi,2Cmix))`).
  Harmless redundancy, not a finding.
- A7/A12 (`Z_W` consistency) is mild: it confirms `1/16 = 1/(8·2)` by composing the Stage-047 forward map (`/8`)
  with `M_mix=G_tr/2`. The `/8` map constant is carried from Stage 047 (notes l.133 states the inverse map), so it is
  anchored and the check is a legitimate cross-form consistency test, not tautological. Both the direct form and the
  via-map form start from genuinely different expressions; the assertion can fail if the `/8`↔`/16` bookkeeping were
  wrong. Acceptable.
- Λ_twin,req (card eq.) is printed (sympy l.115–120) but not independently asserted against a second derivation; it
  is, however, just `π² Π_tr/[16(1-eps)]` = `Π_tr / (2 C_mix · ...)` rearranged, and the boxed-criterion boundary
  (A5) plus the C_mix definition pin the same content. Not a separate finding given A5 covers the load-bearing
  inequality and `Π_tr` itself is asserted (A1).

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.txt:1` (mtime 2026-05-22 17:18)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.txt:1` (mtime 2026-05-22 17:31)

**What's wrong:**
Both committed transcripts predate their scripts (`.py` mtime 2026-06-03 15:59; `.wl` mtime 2026-06-03 15:59).
The captured content also disagrees with the current script's banners, confirming the staleness is material to the
labels: the SymPy transcript header reads `STAGE 34 — EXACT TRACKING-BRANCH PRODUCT` (output l.3) and
`STAGE 34 THEOREM LEDGER` (output l.123), while the current script emits `banner("STAGE 51 — EXACT TRACKING-BRANCH
PRODUCT")` (sympy l.63) and `banner("STAGE 51 THEOREM LEDGER")` (sympy l.149). The Mathematica transcript header
reads `STAGE 034 — LOWEST TWIN CRITERION` (output l.3) while the current script emits
`banner["STAGE 051 — LOWEST TWIN CRITERION"]` (wl l.32). The numeric/algebraic results in the transcripts are
otherwise correct and consistent with the current scripts (all checks `= 0`, all `PASS`); only the stage-number
banner text is out of date.

**Why this matters:**
A checkpoint stage's committed transcript should reflect the current script. Here the body math is unchanged, so the
staleness is label-only, but the "STAGE 34"/"STAGE 034" headers on a Stage-051 transcript are misleading on their
face. The orchestrator's independent re-run will refresh both.

**Required change:**
Re-run both scripts to regenerate the committed transcripts (`python3 .../stage051_..._sympy_audit.py >
.../output/...sympy_audit.txt` and `math -script .../stage051_..._mathematica_audit.wl >
.../output/...mathematica_audit.txt`). No script-math edit is required for the math to pass; this finding is the
freshness signal. See also the stale self-labels noted below (script-side, low severity).

**Verification:**
After regeneration, the SymPy transcript header should read `STAGE 51 — EXACT TRACKING-BRANCH PRODUCT` and
`STAGE 51 THEOREM LEDGER`; the Mathematica transcript header should read `STAGE 051 — LOWEST TWIN CRITERION`; all
`= 0` / `PASS` lines unchanged.

## Stale self-labels (low-severity numbering class — informational, non-blocking)

Per the orchestrator's known `stale_output` numbering class, these self-labels in the script source carry a stale
stage number from the EM-extension renumber and should be corrected when the orchestrator does the label pass (the
canonical topic is stage 051):
- `scripts/.../stage051_..._sympy_audit.py:3` — docstring filename line reads
  `moving_throat_pde_stage34_lowest_twin_criterion_sympy_audit.py` (should be `stage051`).
- `scripts/.../stage051_..._sympy_audit.py:21` — provenance note "Stage 050/034 product law" — the `/034` is a
  stale cross-reference (the Stage-048→050 zeta map / Stage-047 product law is what is actually carried; `034` is
  the pre-renumber index of this very stage).
- `scripts/.../stage051_..._sympy_audit.py:63,149` — banners say `STAGE 51` (2-digit; cosmetic, current topic).

The Mathematica banners (`STAGE 051`, l.32; `Stage 051`, l.114) are already correct. These label items are flagged
per policy but are not blocking and do not change the math; the orchestrator resolves label scope separately.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. The strongest evidence is `xi_(2x)`: where SymPy
only *substitutes* the hand-written root and checks `G_tr(xi_2x)-2 M_mix==0` (sympy l.147), Mathematica
independently solves the quadratic via `Solve[{gTr == 2 mMix, xi > 0}, xi, Reals]` (wl l.100), strips the
`ConditionalExpression` wrapper, and compares the *derived* root to the closed-form claim (`xi2xDerived - xi2xClaim`,
wl l.110) before the substitution check (wl l.111). The Solve-derived form printed in the transcript
(`(-9*delta + 2*mMix*(9 + 2*r^2) + Sqrt[81*(delta + 2*mMix)^2 - ...])/18`, output l.33) is a genuinely different
surface form from the claim (output l.34), and the `==0` check (output l.35–36) proves their equality. The product
`Pi_tr` is canonicalized via an independent `Factor[Together[...]]` route (wl l.48) rather than SymPy's
`simplify(factor(...))` choreography. Endpoint robustness uses the `1/pi1==0` idiom (wl l.68) rather than SymPy's
`is sp.oo` identity test. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at the level claimed. SymPy `Pi_tr - expected closed form = 0` (output l.27) ↔ Mathematica
`Pi_tr - expected closed form = 0 / PASS` (output l.10–11). Endpoints: both give `Pi_tr(0+)=0` and `Pi_tr(1-)=oo`
(sympy output l.32–33; wl output l.16–17, the latter as `Infinity`). The three `zeta_req` boundary checks are `0` in
both (sympy output l.56–58; wl output l.22–27). `Z_W^(twin,req) - forward-map = 0` in both (sympy output l.108; wl
output l.31–32). `G_tr(xi_(2x)) - 2 M_mix = 0` in both (sympy output l.120; wl output l.37–38), and Mathematica
additionally confirms `xi_(2x): Solve vs claim = 0` (output l.35–36). No `engine_disagreement`.

## Verdict justification

`findings` — one low-severity `stale_output` finding (both committed transcripts predate the scripts and carry
"STAGE 34"/"STAGE 034" headers that the current scripts no longer emit). The math itself holds up under attack:
`Pi_tr` is built independently from `F_tr*G_tr` and matched to the closed form; the threshold equivalence
`zeta_req<=1 <=> Pi<=2 C_mix` is pinned at the boundary in both engines; the `Z_W` consistency uses the anchored
Stage-047 `/8` forward map and is non-tautological; and `xi_(2x)` is independently Solved in Mathematica and
verified by substitution in both. Attacks tried and failed: (1) tautology on `Pi_tr` — defeated, `Pi_expected` is
written independently of the `F_tr*G_tr` product; (2) tautology on the `Z_W` map check — defeated, the via-map form
starts from the inverse Stage-047 map and the `8↔16` bookkeeping is exercised; (3) symbol-domain error on
`(1+chi0)^2` — no impact, both sides carry the identical factor; (4) `xi_(2x)` hand-baked — defeated by the
Mathematica `Solve` derivation; (5) paper↔script value drift — none, the boxed criterion, `C_mix`, `Λ_twin,req`,
and `M_mix^(twin,req)=G_tr/2` all match the card verbatim and `Z_W`/`xi_(2x)` match the notes. I read the card,
the notes, and the appendix rows; the script's claim set matches the paper's claim set.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `C_mix = 8Λ(1-eps)/π²` | py l.88; out l.38–42 (`8·Λ·(1-eps)/π²`) | tex l.16 (`eq:app-stage051-Cmix`); md l.183 | MATCH |
| boxed criterion `Pi_tr <= 2 C_mix = 16Λ(1-eps)/π²` | py l.151; out l.126 | tex l.22–23 (boxed); md l.97, l.186 | MATCH |
| `Pi_tr` closed form `xi(xi+δ)(9δ+(9+2R)xi)²(9δ+(9+2R²)xi)/[9(1-xi)(...)²]` | py l.71–74; out l.20–26 | md l.49–53 | MATCH |
| `Pi_tr(xi->0+) = 0` | py l.78; out l.32 | md l.59 | MATCH |
| `Pi_tr(xi->1-) = +infinity` | py l.82; out l.33 | md l.61 | MATCH |
| `Λ_twin,req = π² Π_tr/[16(1-eps)]` | py l.115; out l.88–94 | tex l.28 (`eq:app-stage051-Lambda-threshold`); md l.112, l.189 | MATCH |
| `M_mix^(twin,req) = G_tr/2` | py l.116; out l.96–100 | tex l.30; md l.125, l.192 | MATCH |
| `Z_W^(twin,req) = π²(1-eps_eta)(1-eps)G_tr/[16(1+chi0)²]` | py l.117; out l.102–107 | md l.141–143 | MATCH (notes carrier) |
| `xi_(2x) = [2M_mix(9+2R²)-9δ+sqrt((...)²+648 M_mix δ)]/18` | py l.136–143; out l.113–119; wl l.105 | md l.161–163 | MATCH (notes carrier) |

INTERNAL (scaffolding, no finding expected in prose): `G_tr`, `F_tr` (intermediate transport factors, both also in
md §1 anyway), `S_req = Pi/C_mix`, `zeta_req` symbolic form, `ZW_from_Mmix` forward map, `zeta_req_branch`
(physical-branch substitution), all `expect_zero`/`PASS` residual flags.

reconciliation: complete; 9 values checked, 0 misaligned
