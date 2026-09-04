# _measurements — S11c_c1_wl_build_directive.md (the 2-leg decision gate + the folds)

The S11c-c1 Wolfram build directive is **orchestrator-written**, so its rule-7 decision gate is **Codex + Grok**
(the authorship table: whatever writes does not review). Both legs read the sole-authority c1 spec + the sibling
S11c-a/S11b specs + the SymPy build directive + `S9…:14-18` (the one blindness control) + `N8`/`N9`, formed their
own view first, then reviewed the directive. Leg prompt: `directives/_legs/S11c_c1_wl_build_directive_review.md`.
Raw reports outside the repo (tree hygiene): `…/scratchpad/{codex,grok}_wl_dirreview.log`.

## Launch (detached — escapes the ~87 s `run_in_background` reap)
```
setsid nohup bash -c "cd /var/projects/toy_physics && codex exec -c model_reasoning_effort=xhigh \
  --sandbox workspace-write \"$(cat …/_legs/S11c_c1_wl_build_directive_review.md)\" > …/codex_wl_dirreview.log 2>&1; \
  echo EXIT=$? > …/codex_wl_dirreview.marker" < /dev/null & disown
setsid nohup bash -c "cd /var/projects/toy_physics && grok --prompt-file …/_legs/S11c_c1_wl_build_directive_review.md \
  --cwd /var/projects/toy_physics --model grok-4.6 --effort high --permission-mode bypassPermissions \
  --output-format plain > …/grok_wl_dirreview.log 2>&1; echo EXIT=$? > …/grok_wl_dirreview.marker" < /dev/null & disown
```
Both `EXIT=0`. Codex log 1.17 MB (162 k tokens); Grok log 20 KB. Neither reviewed the other's output.

## Verdicts
- **Codex:** "Do not release this directive to the builder" — **5 MUST-FIX**, no NITs. Core DtN/response/
  dissipation/emit-list/WL-run obligations otherwise wired correctly.
- **Grok:** "Not safe to build as-is" — **3 MUST-FIX** (= Codex's 1,2+3,5). §3a/§3b/Λ_X wiring, emit list, three
  clauses, and run-time "imports nothing" otherwise sound. Grok **explicitly judged the blindness framing sound
  and the copy-test correctly-scoped** (disagreeing with Codex's finding 4 as a full defect).

## Orchestrator rule-13 verification of the load-bearing findings
- **`mu_theta_operator` is a composite, not an atom** (Codex #3 / Grok Finding 2b) — CONFIRMED. Command:
  `python3 -c "import importlib.util; …; print(str(LEDGER['mu_theta_operator'])[:240])"` →
  `{'display': '(((LAB_HELD, RHO4_CONSTANT), ((VALUE, (mu_theta, B_rho_3*epsilon_shape*eta_bg*theta*w1_profile + …`
  — an anchoring/density-keyed `Tuple` carrying the full expanded `μ_θ`, ~216 KB, NOT `Symbol('mu_theta')`. Both
  legs independently confirmed the committed SymPy c1 engine represents `μ_θ` as per-anchoring/face opaque locals
  `s11cc1_mu_theta_{anchoring}_{face}` (engine `:261-265`) + a composite `FACE_RESPONSE` sidecar
  `OPAQUE_MU_THETA_OPERATOR` (`:1266`). My directive's "same single opaque atom" claim was false.
- **`mechanical_lower_camel` spellings** (Codex #5 / Grok Finding 3) — verified against the real function
  (`scripts/S11b_cross_engine_comparator.py:147`: `pieces=name.split("_"); pieces[0]+"".join(p[:1].upper()+p[1:]
  for p in pieces[1:])` — the FIRST piece keeps its case). Ground-truth map: `W_bg→WBg`, `w1_profile→w1Profile`,
  `eta_bg→etaBg`, `epsilon_shape→epsilonShape`, `sigma_W→sigmaW`, `rho_br→rhoBr`,
  `rho_br_bg_rho4_constant→rhoBrBgRho4Constant`, `Lambda_A_0→LambdaA0` (⚠ capital L), `tau_A→tauA`, `q_out→qOut`,
  `rho_m→rhoM`, `c_s0→cS0`, `W_0→W0`, `e_W→eW`, `v_bulk_normal_0→vBulkNormal0`, `mu_theta→muTheta`.

## The five folds applied (one pass — rule 7: fold and go, ⛔ not re-legged)
1. **T-a "check against" (Codex #1 / Grok F1) — cut.** The directive no longer quotes the c1 spec §2a expanded
   normal as "an input to check against"; it now instructs re-deriving T-a..T-i from the S11c-a maps + level set
   `F_s^α` + orientation law `s(n̂·ŵ)>0` (⛔ not `sign(∇₄F)`), states §2a's expanded normal is a computed T-a
   **result** (NOT a derivation input or target), and requires the §5b `W_bg`-tilt ablation to MOVE the computed
   T-a (byte-identical = transcribed constant = a finding).
2. **`μ_θ` "cancels out of every residual" (Codex #2 / Grok F2) — deleted.** No pre-registered agreement (rule 16
   / spec §7 "no representational fold pre-registered").
3. **`μ_θ` atom-identity (Codex #3 / Grok F2b) — corrected.** The directive now states the row is a composite and
   the committed SymPy engine's per-case-local + composite-sidecar representation; instructs WL to build its
   **own** opaque `muTheta`, ⛔ never copying PY's `s11cc1_*` locals/composite (blindness + designed-to-agree),
   with the PY↔WL `μ_θ` representational residual **adjudicated after the run** (like the omega artifact).
4. **Blindness/copy-test (Codex #4 — the valid precision part only; Grok judged the framing sound).** Acceptance
   check #1 rescoped to **run-time** import-freedom, with a structural scan for
   `Get`/`Import`/`<<`/`ReadString`/`OpenRead` + absolute repo path, run outside the repo; noted that
   construction-time blindness comes from the builder being handed only the specs (⛔ the "denylist" reframing was
   NOT taken — it matches the committed master WL directive + N9, per Grok).
5. **Shared-symbol table (Codex #5 / Grok F3) — replaced** with the exact reserved names in verified
   `mechanical_lower_camel` spelling; the two density reps kept distinct (⛔ not one `ρ_br,bg⁰`); the frequency
   kernel `Λ_I(ω)` marked a constructed object (from `LambdaI0`,`tauI`,`omega`), not a bound symbol; "only emitted
   names are standardized" clarified to govern TAG STRINGS, not payload-symbol identity; the `k`/`k′` legs named
   WL-legally (`k`,`kPrime`), ⛔ not PY's `s11cc1_k_*`.

## Noted latent observation (⛔ NOT folded — carried to the T7 join / step record)
- **The c1 spec §2a (committed `f90e7630`) already exposes the first-order outward normal
  `n̂_s^α=(−½∇_x W_bg^α,s)+O(|∇W_bg|²)` (a T-a result) as if supplied.** Both legs flagged it; Codex asked to
  repair the spec too. Adjudication: the directive **neutralizes** it (fold 1) rather than reopening the
  committed, 2-leg-reviewed spec (whose content sha the committed SymPy engine's `BUILD_INPUT_DIGESTS` pins — a
  spec edit would cascade a SymPy re-verify). The method's actual defense holds: the WL T-a is **computed** from
  the setup, the §5b form ablation must move it, and the WL-recomputed-T-a vs SymPy-imported-S11c-a-T-a residual
  is a real cross-engine control. If the exposure bites at T7 (a suspiciously-zero substrate residual), revisit
  the spec then.
