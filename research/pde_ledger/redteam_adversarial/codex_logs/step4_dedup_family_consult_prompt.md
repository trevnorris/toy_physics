# Step 4 consult — dedup/family rules produced unsound proposals (Claude+Codex, read-only)

You are the Codex side of a structured consult per the signed-off adjudication shape (execution plan §5; operating-contract item 4 makes dedup/family/batch granularity a Claude+Codex call, NOT a user gate — this is not conceptual). Read-only sandbox: verify against the files, propose, don't edit. End with a `## CONSULT VERDICT` section, per-question AGREE / DISAGREE-with-alternative, grounded in files you open.

## What happened

Step 4 built Phase B tooling. The CODE passed a clean review (both load-bearing checks — `phase-b-ingest` rejects any origin/constraint claim without a `notes/` citation; `benchmark-ingest` rejects any entry without a real source citation — are correctly enforced). But the two read-only PROPOSALS the tooling generated are unsound. The agreed D1 process placed a clean-agent review BEFORE applying the alias map precisely to catch this; it fired.

Files (all under `/var/projects/toy_physics/`, read-only):
- `research/pde_ledger/redteam_adversarial/provenance/_dedup_proposal.yaml` — 140 alias groups, 291 aliases, 631 canonicals, 342 ambiguous.
- `research/pde_ledger/redteam_adversarial/provenance/_family_map.yaml` — 1358 families, 1329 singletons, 0 unmapped (marked `pre_alias_map`).
- `research/pde_ledger/redteam_adversarial/fit_insertion_points.yaml` — human-readable candidate union.
- Tooling: `.claude/skills/adversarial-audit/lib/core.py` (dedup `build_dedup_proposal` ~:635 + helpers `records_value_compatible` ~:306, `citations_adjacent` ~:284; family `build_family_map` ~:857).

## Finding 1 — dedup over-merges (orchestrator-VERIFIED on stage073)

The alias group `fit_stage073_Lambda_ell_fixed_37` (8 members, `_dedup_proposal.yaml` ~line 2245) merges SEVEN distinct parameters at different lines: `epsilon_r=0.05=1/20` (line 30), aspect ratio (line 48), `Lambda_ell` (line 72), `eta` (line 94), plus `k_m`, `t_x`, `Lambda_star` — span lines 6–94. These are NOT the same physical insertion point. The group's own `adjacency_rule` text says "line distance <= 3", but the operative merge is the `records_value_compatible` path: same `anchor_stage` + overlapping `parameter_family` union + no *conflicting* extracted literal value — and most members carry `literal_values: []`, so nothing conflicts and everything in the stage's parameter soup transitively merges. A clean review agent found 127/140 groups (91%) span >3 lines; ≥4-member groups (stage073=8, stage098=7, stage111=4, stage247=12, stage026=4) systematically collapse distinct fitted values. The sound merges are the small 2-member ones (paper-vs-script-vs-pass2 of the *same* value, e.g. stage026 κ_n=√2/((n+½)π); path-prefix dupes; mechanical-key vs agent-key for the same (stage,param,value)).

**A false merge is the dangerous error**: aliases inherit the canonical's Phase C verdict, so e.g. `epsilon_r=1/20`, the `eta`-pin, and the 37/20 aspect ratio would all get ONE verdict — a real fit could hide behind another candidate's PASS.

## Finding 2 — alias/ambiguous self-contradiction

~153 candidate ids appear in BOTH the alias-merge and the ambiguous-keep-separate sections (e.g. stage026 `g_u=1/2` and `k_req=54/5` are declared aliases into `fit_stage026_K_req_solved_exactly_from_target` AND members of a stage-026 ambiguous group flagged "conflicting literal values → stay separate"). The ambiguous bucket's call is the correct one; the alias merge is wrong. No apply-time tiebreak exists.

## Finding 3 — family map fragments every headline chain into singletons

Family key = normalized `parameter_names` set + EACH extracted literal value. Consequences: (a) the same conceptual parameter at a different stage/value forms a NEW family, so the load-bearing cross-stage chains fragment — 54/5 quadrupole, 37/20 aspect cascade, and Sigma0_can/T_can are all scattered across singleton families; the Σ0 **…867 vs …876 digit-transposition** carriers land in SEPARATE families (and the two purpose-built digit-variant candidates have empty value-keys), so the exact drift the audit targets is invisible at family level; (b) noisy multi-name candidates explode — `fit_stage_105_chi_q` maps to 18 families, `fit_stage_099_k0_sym`→13. Only `chi_Q=1` (`fam_0332`, 12 members) clustered well. Root cause is upstream: Phase A `parameter_names` are whole-sentence symbol captures (one candidate carries ~28 names), so the param-set+value key is both too sparse (values empty) and too noisy (param unions).

## Claude's proposed fix (verify and refine)

**Q1 — Confirm Finding 1.** Open `_dedup_proposal.yaml`, verify stage073 (and one more ≥4-member group, e.g. stage247 or stage098) are false merges, and confirm the operative path is `records_value_compatible` over a `parameter_family` union, not the advertised 3-line adjacency. AGREE the dedup map is unsafe to apply as-is?

**Q2 — Dedup fix = strict, value-anchored, conservative.** Re-spec `dedup-propose`: alias two candidates ONLY when same `anchor_stage` AND same single normalized target parameter AND same extracted literal value AND (identical `(path,line)` OR same normalized path within ≤3 lines OR a pure path-prefix variant). DROP the `records_value_compatible` no-conflicting-value fallback entirely. Any candidate with empty `literal_values`, or with a multi-parameter (polluted) name set that can't be reduced to one target, is NOT merged → it stays a separate canonical (and the alias/ambiguous self-contradiction disappears because nothing is both). Expect the canonical count to rise substantially (toward ~800+); accept that — over-merge hides fits, under-merge only costs redundant Phase C audits. AGREE?

**Q3 — Family fix = conceptual-parameter identity, value as attribute.** Re-spec `family-build`: derive ONE primary conceptual parameter per candidate from the `candidate_key` / `reason` (e.g. `fit_stage073_eta_pinned_37`→`eta`; `fit_stage_105_chi_q`→`chi_Q`) rather than the polluted `parameter_names` union; key the family on that normalized conceptual parameter ALONE (NOT param+value); carry each occurrence's literal value as a per-member ATTRIBUTE so cross-stage occurrences LINK and value divergences (Σ0 867 vs 876) surface as a within-family `value_divergence` finding. One primary family per candidate (kill the 18-way explosion); a small number of genuinely multi-target candidates may map to >1. Remaining singletons in the long tail are acceptable (genuinely distinct single insertions). AGREE?

**Q4 — Seed the known chains?** The D3 consult already enumerated the headline chains (54/5+gamma_GR; 37/20→F1 incl. 4107; chi_Q/m0_hat/N_Q; Sigma0_can/T_can incl. 867/876; barrier 222–224 incl. V_known; calibration 245–253 incl. lambda_L+CODATA; wall-action+m̂/P0). Should `family-build` seed these as explicit named families (a safety net so they cannot fragment), with mechanical membership filled by the conceptual-parameter rule and any candidate the rule misses flagged for review — or do you prefer to trust the re-keyed rule alone and verify coverage post-hoc? Recommend one.

**Q5 — Anything missed.** Is purely-mechanical grouping the right tool here given the upstream `parameter_names` noise, or should one step be agent-assisted? Should we instead fix the Phase A `parameter_names` pollution at the source (re-derive a single target parameter per candidate) before re-running dedup/family — and if so, mechanically or via agents? Flag anything that should change before I write the rebuild directive.
