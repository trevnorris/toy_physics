# _measurements — S11c_c1_repair_directive.md (the repair fold + its 2-leg re-review)

The repair directive folds the five findings of the c1 build review (`_measurements/S11c_c1_build_review.md`;
2 legs, Codex-written engine ⇒ fresh Claude Agent + Grok). Each finding's provenance is that record; the repair
was itself 2-leg re-reviewed (prompt `directives/_legs/S11c_c1_repair_review.md`). Both re-review legs derived
independently and ablated each control on /tmp copies (working tree untouched). **Both verdicts: CLEAR / SOUND —
no findings.** Raw reports outside the repo (tree hygiene). Literal ablation results below (each control now bites).

## Deliverable verification (orchestrator, rule 13; vs committed baseline `65afa1cd`)
- Engine compiles; runs exit 0 (63 tags, no skipped tasks, ~517 s); regenerates the 906,467-byte export.
- `git diff 65afa1cd -- scripts/S11c_c1_exports.py` = **exactly 1 line** — the engine's own self-referential
  digest; all 44 LEDGER rows (incl. the 5 model-level rows) byte-identical (both legs re-confirmed srepr-equal).
- `dtn_first_kernel`, `closed_coefficients`, `dtn_flat_symbol`, `response_operator_case`, `hermitian_kernel`, and
  the 44-key `IMPORT_KEYS` are byte-identical to `65afa1cd`; `ledger_fold.py` unchanged. Diff scope = only the 5
  controls + new helpers (`outgoing_farfield_poynting`, `hanzawa_first_kernel`, `reduce_on_shell`,
  `first_shape_true_area_weight`, `intended_own_bind_closure`, `simplified_object_difference`, and the
  `INTENDED_EXPORT_CANDIDATE_ROOTS`/`INTENDED_EXPORT_WRITE_ROOTS` constants; `port_matrix`/`port_hermitian`/
  `build_model`/`build_delta`).

## The five controls — each now BITES (concurrent literal ablation, both legs)
- **R1 energy (§3b obj.3).** Emitted `ENERGY_RESIDUAL` eval at a propagating point: `BASELINE 0`; `t_s`-sign-flip
  (face only) `+4`; `q_out`-branch-flip in `outgoing_farfield_poynting` (bulk only) `−4`. The two operands share
  no `z_energy` factor (a shared factor would keep the residual 0 under a φ-only branch flip). Both legs derived
  the far-field acoustic intensity independently and it MATCHES the engine (Claude: FLAT 14.3125, GENERAL −48.66
  vs bare ½Re(z0|v|²)=1.72 — carries curvature cross-terms, the outgoing intensity at |w|→∞, not δp at the face).
- **R2 rep-invariance (§5a).** `hanzawa_first_kernel` = `OUTGOING_NEUMANN_LAYER_POTENTIAL`, kernel
  `exp(iq_out n)/(iq_out)`, radiation-selected branch, ⛔ NOT the secular `w′=[w−ζ_c]/[W_bg+δW]`. `srepr(direct)
  ≠ srepr(hanzawa)` (37 vs 58 ops — not A−A); `EULERIAN−HANZAWA=0` for the correct kernel; one-sided Eulerian
  corruption drives `REP_INVARIANCE_RESIDUAL` nonzero (`−ω ρ_m k_in₁ w1_jet₁ σ_W/(q_in q_out)`). Independent
  Eulerian N-to-D matches the engine Eulerian exactly.
- **R3 on-shell (§3a/§5d).** old `Add(k_i²)` `xreplace` reproduced-missing after `sp.factor`; `reduce_on_shell`
  (`q²→ω²/c_s0²−|k|²` via together/cancel/factor) fires including on the factored form; a wrong kernel does not
  cancel (`I ω ρ_m h ≠ 0`). `DTN_RIGID_SHIFT_RESIDUAL = 0` on all 4 (anchoring,face); `ZERO_JET` operand now `Z_0`
  with no off-shell height term.
- **R4 Hermitian/port Z_0+Z_1 (§3b obj.1–2).** at (η,σ_W)=0 `Z_1=0` but `DTN_HERMITIAN_PART` = `Re(ρ_m ω/q) ≠ 0`
  (= `H_a[Z_0]`, both legs' independent object); old Z_1-only feed vanishes; all port entries nonzero at flat;
  `FULL_BARE_DTN` tags; `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` on 32 non-propagating cases, `Ge(re(·),0)` sign
  tests on 4 propagating cases, no sign in prose.
- **R5 minimum-mode (§D3).** `intended_own_bind_closure` computed from the 5 intended roots' values (5 + 39
  new-symbol closure = 44), `FUNCTION_READS_DELTA_KEYS False`; stray-row injection RAISES `LedgerFoldError: extra
  row(s): …`; a drop RAISES `missing row(s): …`; the old `own_closure=delta.keys()` does not raise.

## No-regression (both legs)
Two-momentum kernel both legs live (form ablation `q_in=q_out` MOVES it); resolvent `[I+(Λ_A/ρ_m²)Z]⁻¹` not a
scalar division; `Λ_X` only in traction; opaque `μ_θ` (no θ/e_W atoms); 44-row own-rows delta, fold 2441+44 zero
overwrites; `check_consumer` + lookup smoke-test bite both ways.

## Non-blocking observations (both legs; ⛔ NOT findings — carried to the c1 step record)
1. The old `dtn_first_kernel(..., route="HANZAWA")` regrouping is unused dead code (the live route is
   `hanzawa_first_kernel`). Inert; optional tidy-up, ⛔ not worth a further build round (rule 15).
2. `ENERGY_RESIDUAL[BASELINE]` prints as an unevaluated `re(...)` tree, not `Integer(0)` — a printer/simplify
   residue that evaluates to 0 and moves under a real error.
3. ⚠ `ZERO_JET_RESIDUAL` is an **omega-assumption identity** — S11b `z_impermeable`'s plain `omega` vs c1's
   `omega(real=True)` are distinct SymPy symbols, so the raw residual prints `ρω/q − ρω/q` (numerically 0), ⛔ NOT
   a thickness dependence. **The T7 comparator must canonicalize omega assumptions** (relevant when c1's σ_W→0
   limit is compared to the S11b flat object, and cross-engine to the WL engine) — [[project_dual_engine_symbol_transliteration]].
