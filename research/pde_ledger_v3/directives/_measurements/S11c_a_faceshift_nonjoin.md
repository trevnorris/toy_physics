# FACE_SHIFT non-join — measured (2026-08-27)

Resolving task #69: the `FACE_SHIFT` (T-e shifted-trace) family had `join=0` in the reviewed
comparator run and was **uncompared**. This file records the literal operands and the computation
behind the diagnosis. Everything below is a printed CAS object + the command that produced it; no
conclusion is asserted by a script.

## 1. The accounting line (reviewed comparator run, `~/.s11_build/comparator_run.out:1796`)

```
ACCOUNTING FACE_SHIFT {join=0, py_only=160, wl_only=80, duplicate_key=0, parse_failed=0, axis_set_mismatch=240} reasons=WL missing DENSITY
```

PY keys each case `(BRANCH, DENSITY, FACE, DOF, FIELD)`; WL keys `(BRANCH, FACE, DOF, FIELD)` — **no
DENSITY axis**. The key tuples never match ⇒ `join=0`. PY has 160 cases (2 branch × **2 density
representatives** × 2 face × 2 dof × 10 field-components), WL has 80 (no representative split).

## 2. The keying gap is degenerate — PY's two representatives are byte-identical

`research/pde_ledger_v3/directives/_measurements/S11c_a_faceshift_nonjoin_verify.py`
(reads the two committed transcripts via the comparator loader):

```
PY cases: 160  WL cases: 80
PY BULK_DENSITY cases: 16
PY reps identical across DENSITY axis for every (branch,face,dof)? True
```

⇒ For `FACE_SHIFT`, PY's `RHO4_CONSTANT` and `RHOBR_CONSTANT` operands are **identical** in every
`(branch,face,dof)`. The DENSITY axis carries no information here; the non-join is a pure keying
artifact (PY carries a degenerate axis, WL does not).

## 3. The only cross-engine content is a WL-only density-background term

Literal operands (LAB_HELD, MINUS, DELTA_W, FIELD=BULK_DENSITY), from `…_verify.py`:

```
PY  value: -W_0*delta_rho_4D_face_minus_dw*eta_bg*w1_profile/2 + delta_rho_4D_face_minus
WL  value:  W_0**2*e_W*eta_bg*rhoBulkBackgroundXJETXdwXdw*w1_profile/4
            - W_0*e_W*rhoBulkBackgroundXJETXdw/2
            - W_0*eta_bg*rhoBulkPerturbationXJETXdw*w1_profile/2
            + rhoBulkPerturbation(x1, x2, x3, (-W_0/2, time))
```

- WL's perturbation terms (`rhoBulkPerturbation` face value + `rhoBulkPerturbationXJETXdw`
  re-centering jet) match PY's `delta_rho_4D_face_minus` + `…_dw` term exactly under the field-name
  bridge (`rhoBulkPerturbation ↔ delta_rho_4D_face`, `e_W ≡ δW/W_0 ↔ eta_bg·w1_profile`, spec §45).
- The **only** WL terms without a PY counterpart carry `rhoBulkBackgroundXJETXdw` (`∂_wρ⁰`) and
  `rhoBulkBackgroundXJETXdwXdw` (`∂²_wρ⁰`) — the shifted-trace background term `δh_s·∂_wρ⁰` and its
  higher-order re-expansion.

## 4. Grounding the density background ⇒ WL − PY = 0 for all 8 cases

`research/pde_ledger_v3/directives/_measurements/S11c_a_faceshift_nonjoin_residual.py` applies the
perturbation field-name bridge and sets `∂_wρ⁰ = ∂²_wρ⁰ = 0` (§2b: the supplied density background
depends on the in-plane anchor, not `w`). Literal stdout:

```
('LAB_HELD','MINUS','DELTA_W',  grounded_resid=0, wl_only_term = W_0**2*eta_bg**2*rhoBulkBackgroundXJETXdwXdw*w1_profile**2/4 - W_0*eta_bg*rhoBulkBackgroundXJETXdw*w1_profile/2)
('MATERIAL_ADVECTED','MINUS','DELTA_W', grounded_resid=0, wl_only_term = … rhoBulkBackgroundXJETXdw / …XdwXdw only …)
('LAB_HELD','MINUS','ZETA_C',  grounded_resid=0, …)
('MATERIAL_ADVECTED','MINUS','ZETA_C', grounded_resid=0, …)
('LAB_HELD','PLUS','DELTA_W',  grounded_resid=0, …)
('MATERIAL_ADVECTED','PLUS','DELTA_W', grounded_resid=0, …)
('LAB_HELD','PLUS','ZETA_C',  grounded_resid=0, …)
('MATERIAL_ADVECTED','PLUS','ZETA_C', grounded_resid=0, …)
```

Every `wl_only_term` is a linear combination of `∂_wρ⁰` and `∂²_wρ⁰` — **nothing else survives**.
Set them to zero ⇒ residual is exactly 0 in all 8 cases.

## 5. Spec grounding (the physics the residual rides on)

- **§3c**: `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰)`; *"Every background face value or normal
  derivative appearing in this law is obtained by differentiating a member of the supplied background
  state 𝔅⁰ (§2d); none may be introduced as a free premise"*; *"the supplied density background
  depends on the in-plane anchor, not on `w`."*
- **§2b**: both representatives (`ρ_4D,bg⁰ = ρ_4D,ref⁰` or `rho_br/W_bg(y)`) are functions of
  `W_bg(y)` — in-plane `y`, not `w`.
- **§2d**: `ρ_4D,bg⁰`, `ρ_br,bg⁰` are members of the supplied state `𝔅⁰`.

⇒ `∂_wρ⁰` must be **obtained by differentiating the supplied (w-independent) background** = 0.

## 6. WL source — the density background is the ONE ungrounded field

`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl:433-435`:

```
pressureZero[coordinates_List, normal_] := 0;                                    (* grounded to 0 *)
rhoBulkZero[coordinates_List, normal_]  := rhoBulkBackground @@ Append[coordinates, {normal, time}]; (* ABSTRACT free premise, w-dependent *)
rhoBulkWave[coordinates_List, normal_]  := rhoBulkPerturbation @@ Append[coordinates, {normal, time}];
```

`rhoBulkBackground` has **no definition anywhere** in the WL source (single usage, L435). WL
explicitly zeroed the pressure/velocity/current backgrounds but left the density background as a bare
abstract function of `{w, time}`, so `∂_w rhoBulkBackground` is a live nonzero symbol.

## 7. Verdict (orchestrator, pending 2 from-spec adjudication legs)

Same class as the WL background-current finding (`6fae82b8`): a background quantity carried as a
**free premise** instead of grounded in the supplied state, which §3c forbids and which the supplied
premise forces to zero. PY (grounded, w-independent) is §3c-faithful; WL is not. The residual is
physically zero but the object WL *emits* differs from PY's. Disposition (fix WL to ground the density
background vs. record the representational difference) goes to two independent legs + the user.

## 8. Two independent from-spec adjudication legs — both CONFIRM (2026-08-27)

The WL engine is Codex-written ⇒ legs = fresh Claude Agent + Grok. Prompt (neutral, no leaked
verdict): `directives/_legs/S11c_a_faceshift_density_adjudication.md`. Both derived the correct §3c
density trace from the spec with runnable SymPy + literal stdout and reached the same verdict.

- **Grok** (`~/.s11_build/faceshift_adj/grok.txt`; scripts under
  `research/pde_ledger_v3/_measurements/s11c_a_density_trace_adjudication_grok/`): PY §3c-faithful;
  WL's `rhoBulkBackground` is an undefined w-dependent free premise §3c bans; `∂_wρ_4D,bg⁰=0` for both
  representatives (differentiated the 𝔅⁰ member); WL−PY is a **spec-forced zero**, not a missing term;
  no spec ambiguity. Noted WL already has the grounded in-plane `rho4Profile` used elsewhere.
- **Fresh Agent** (scratchpad `faceshift_density_adjudication.py`): same verdict; `∂_wρ⁰=0` for all
  four maps + the material-advected branch (stdout); PY does an actual `sp.diff` on the grounded
  `density_pair` → 0; WL keeps `∂_w`/`∂²_w rhoBulkBackground` live. Flagged an honest **wording seam**:
  §3c:375 scopes the law to a "nonzero background face value **or derivative**"; the density's *value*
  is nonzero, but the shift term multiplies the *derivative* (=0), so the density qualifies for the
  shifted-perturbation evaluation while its background shift term vanishes — a seam, not a license.

## 9. Verified: WL is internally inconsistent (grounded density everywhere but T-e)

`mathematica/S11c_a_interface_geometry_mathematica_audit.wl` DOES define grounded, in-plane-only,
representative-aware density profiles and uses them everywhere except the shifted trace:

```
305: rho4Profile["RHO4_CONSTANT", coordinates_List] := rhoBr/W0;
309: rho4Profile["RHOBR_CONSTANT", coordinates_List] := rhoBr/(W0*(1 + ...));   (* in-plane, w-independent *)
472-529: affinity/traction use rhoBrProfile[density, coordinates]   (* grounded *)
653,684-739: objects use rho4Profile[density, ...]                  (* grounded *)
880/901: BACKGROUND_DENSITY_MAP built from the grounded profiles
434-435: rhoBulkZero := rhoBulkBackground[...,{normal,time}]        (* T-e ONLY: ungrounded free premise *)
```

⇒ The fix is trivial and well-scoped: point `rhoBulkZero` at the grounded `rho4Profile[density,coords]`
(w-independent, representative-aware) instead of `rhoBulkBackground`. Then `∂_wρ⁰=0` emerges by
differentiation (§3c-compliant), WL matches PY, AND WL gains the DENSITY axis ⇒ the family JOINS with
PY's two (byte-identical) representative cases, residual → 0. This mirrors the bg-current fix `6fae82b8`.

Confirmed by the comparator accounting: PROJECTION/EVOLUTION families are CLEAN (join=8, no mismatch) —
those use the grounded `rho4Profile`/density map, so the ungrounded `rhoBulkBackground` is confined to
T-e and does NOT propagate to any physics-bearing S11c-a result.

## 10. WL grounding fix built + reviewed; MATERIAL-transport concern adjudicated benign (2026-08-27)

**Fix (Codex build, directive `S11c_a_wl_density_grounding_fix_directive.md`):** `traceSource` now computes
the background face value + its normal derivative by a real `D[·,w]`; the density field is grounded via
`rhoBulkRepresentativeZero[density] := rho4Profile[density, coordinates]` (in-plane §2b representative);
all three shifted-trace sites (primary L1120, form-control L1883, uniform-limit L2096) loop `{density}`
and key by `densityFaceCaseKey`. Result: FACE_SHIFT gains the DENSITY axis, `rhoBulkBackground` gone.

**Cross-engine check (`…/scratchpad/faceshift_postfix_check.py`):** WL FACE_SHIFT cases 80→160; key axes
now `[BRANCH,DENSITY,DOF,FACE,FIELD]`; density cases **joined 16/16**, **all residuals 0**; WL density
operand = perturbation re-centering only (`rhoBulkBackground` absent).

**2 build legs (fresh Agent + Grok):** grounding confirmed §3c-correct; **form ablation bites** (swap the
grounded w-independent background for a w-dependent test fn ⇒ the shift term reappears ∝ `deltaWidth` ⇒
genuinely `D[·,w]`-wired, not a hardcode); collateral clean (only FACE_SHIFT + the 3 UNIFORM_LIMIT slices
changed, zero non-FACE_SHIFT lines; form-control's FACE_SHIFT slice rides inside `CONTROL_FORM_*`);
blindness preserved. Leg 2 (Grok) raised a MATERIAL_ADVECTED×RHOBR concern → adjudicated below.

**MATERIAL-transport adjudication (2 from-spec legs, Agent + Grok, each with SymPy + stdout):** §3a
mandates `ρ_4D,bg^{0,M} ≡ ρ(χ)`; both engines emit the MATERIAL×RHOBR density trace with NO
material-transport term (identical to LAB), so the cross-engine residual (0) can't adjudicate it — a
rule-7 shared-blind-spot candidate. Verdict: **the omission is correct.** §3c's trace law is normal-shift
only (`δh_s·∂_w f⁰`) and `∂_wρ⁰=0`, so the in-plane transport `−u·∇ρ⁰` is not a §3c face-shift term; it
legitimately enters `δh_s` (the MATERIAL face displacement, where it multiplies `∂_wρ⁰=0` ⇒ contributes
nothing to the density trace) and the density representative consumed by T-f/T-g/T-h (projection/evolution
carry `advective = σ_W·u·∇w1`). Injecting it into T-e would **double-count** (both legs showed the residual
= the advection counted twice). NOT a finding — both engines correct, not a shared blind spot.

**One non-physics residual (recorded, not fixed):** the WL `EXACT_TRACE_SOURCE` sub-field displays the
MATERIAL background at lab `ρ(x)` where §3a §3a:316 says `ρ(χ)`. Both adjudication legs confirmed this
changes **no computed object** (at background `u=0 ⇒ χ=x`; `ρ(χ)` also has `∂_w=0`; the shape-derivative
operand is frame-independent). It is a display-fidelity point in an un-differentiated intermediate, not
cross-compared and not a physics defect. Left as-is; noted for the step record.
