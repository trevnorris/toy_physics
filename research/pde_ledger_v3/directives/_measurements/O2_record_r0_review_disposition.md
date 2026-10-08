# O2 record r0: review dispositions (orchestrator)

**Artifact:** the O2 step record, written by Codex (gpt-6.1-sol, session `01a119a8…`) to the record directive
`directives/O2_record_directive.md` (`c9db665d`). Its four files are frozen in
`_scratch/s9b_build/o2_record_review_baseline_r0.sha256`:

| File | sha |
|---|---|
| `steps/O2_steady_brane_balance.md` | `411b8259…` |
| `steps/_measurements/O2_record_measurements.md` | `9fa98a00…` |
| `scripts/O2_record_measurements.py` | `2f89e772…` |
| `SUBSTRATE_REQUIREMENTS.md` (O2 pass) | `e5df1ef9…` |

The record is Codex-written physics-bearing prose, so it gets fresh Claude + Grok, source-first, reviewed until clear
on what may be claimed (G1, G4).

**Legs.** Both used the identical prompt `_scratch/s9b_build/o2_record_review_prompt.md`, with no build directive in
the packet. Both reported before adjudication.
- **Grok:** NEEDS REVISION, 1 finding (`_scratch/s9b_build/o2_record_review_r0_grok.txt`; evidence in
  `…_r0_grok_evidence/`).
- **Fresh Claude (opus):** NEEDS REVISION, 6 findings (report in this session's transcript; evidence in
  `_scratch/s9b_build/o2_record_review_r0_claude_evidence/`).

**What both legs established independently:**
- **Row classes.** The record's row classes and counts match each leg's own enumeration of the filed comparator
  output: 10 kinds, 160 joined, 204 unjoined, 0 unaccounted, and no nonzero residual.
- **Measurements.** The generator reruns byte-identical to the filed measurements file, by retrieval only.
- **Statuses and premises.** Statuses match spec §§1–3 and 7–8. Premises 1–4 carry their adopted label, and none
  is presented as derived.
- **Owners.** O1 and O3–O7 owners match scoping §1, and none is closed.
- **No leakage.** No expected value, sign or light-bending outcome is stated, and `ρ_br V` is not supplied.
- **Routed items.** The routed items and the comparator's acceptance limits are carried.

Each verification below is a mechanical lookup. The legs' filed stdout is quoted verbatim. Commands and literal
output are in `O2_record_r0_review_disposition_lookups.md`.

| # | Finding (legs) | Disposition | What must be true after repair |
|---|---|---|---|
| F1 | `R-S1-03` (time-reversibility; target S1) had its target qualifier rewritten to "O2 explicitly names S8". O2 has no time-reversibility content, and the entry is not one of the six the pass says gained an O2 source. (Claude 1, Grok 1) | **ACCEPT.** Lookups F1: committed `c9db665d` L373 and L390 carried "register inference; no owner named by the records". The working tree carries the new string at L394 (`R-S1-03`) and L411 (`R-S8-06`). The record has 0 matches for `revers`/`onsager`. | Only entries that O2 sources bear on change. Each change states only what O2's sources support. `R-S1-03` carries no O2 claim. |
| F2 | §4 row 2 calls unbound roles "representational … to that stated extent" on the strength of the comparator's label. The label is printed for every action key missing from the comparator's table. (Claude 2) | **ACCEPT.** Lookups F2: record L173. The label comes from comparator L1108/L1125–1126, a missing table binding. It covers PY rotational work and the WL native normal and hold actions. Face geometry is among the items the engines' acceptance routed as unsettled. | An unbound role is reported as an OPEN cross-engine difference for which no comparison was formed, with the label's provenance cited. It is called representational only where a cited physics source settles it. |
| F3 | The record leaves out comparator results that bear on what may be claimed. (a) In all six balance components, both engines have identical role/orientation/OPEN-free entry multisets. (b) No paired OPEN inventory has an empty difference (234 nonempty, 0 empty). (c) No formed residual is nonzero (136 exact zero, 212 not formed), and several zero rows and the "nested sibling absent" reason are missing from the listing. (Claude 3) | **ACCEPT.** Lookups F3: leg stdout `identical_role_orientation_multisets=True` for every balance component; `paired role entries with nonempty differences: 234 with all-empty differences: 0`; `NONZERO exact residuals: 0`. These are printed comparator fields, which the record may report by retrieval. | The record states, by retrieval and within the printed limits, what agrees and what does not: the balance role/orientation structure, the OPEN-inventory difference counts, and the residual outcome totals. Its zero-row and not-formed-reason listings are complete. |
| F4 | The WL-only derivative content is narrowed to "in-plane". The same eight WL-only live keys appear in all six balance components. (Claude 4) | **ACCEPT.** Lookups F4: record L179–181 says "in-plane". Leg stdout lists `ProfileDerivative` of `V_r, delta, f, h, j_n, mu_perp, o2_rho_br_live` (order 1) and `xi_w` (order 2), each with `balances=[energy_balance, hold_bulk, hold_inplane 0/1/2, hold_normal]`. | Every printed cross-engine live-object difference is reported with every balance component, and every key, it occurs in. |
| F5 | The Part D handoff does not say which content of the named objects is usable. "Accepted emitted geometry" includes the native face geometry, which is never compared. Every OPEN role's inventory differs. Adopting either engine's narrower inventory would drop dependences. The mass-law qualification is not carried. (Claude 5) | **ACCEPT.** Lookups F5: record L302–309. Spec L124–125: "…gradients and material history remain admissible dependences". Spec L330–336 gives the induced-metric relative-`O(ε)` qualification. | The handoff names exactly what Part D may rely on, distinguishing compared from uncompared content. No dependence printed by either engine is dropped, neither engine's OPEN inventory is taken as the object's, and the qualifications attached to supplied inputs travel with them. |
| F6 | `R-S12-01` (target S12) absorbs the O2 net-supplier/budget obligation. The sources keep that obligation separate and name no S12 owner. (Claude 6) | **ACCEPT.** Lookups F6: register L438. Spec L282 makes `𝒮_E,net`/`𝒫_E,supply` an "O2 requirement". Spec L312–314: "This conditional obligation is separate from S12's non-variational source partners…". Contract L459–461 keeps the non-passive-interface condition separate. | Each register entry records one object with the owner its sources name, or "no owner named". An obligation the sources keep separate is not merged. |

**Note: measurements size.** The filed measurements file is 10.9 MB of plain-git text (`O2_record_measurements.md`).
That is acceptable once in the preservation commit. The repair should keep each retrieval's output to what the
cited claim rests on, so that each revision does not re-commit a file of that size. This is not a finding.

**Pattern (G4).** This is the record's first review round and its first repair. No defect was bred by a repair. The
same author is kept.

**Next.** Preserve r0 as the reviewed baseline (not accepted). Resume the same Codex session with a repair brief
stating F1–F6 as what must be true. Then run review round 2 with fresh legs.
