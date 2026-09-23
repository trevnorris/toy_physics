# S11c-d scope and review inventory — decision draft

Prepared 2026-09-22 against `0d77af53561968d4593b49261589e2ee259eadbb`.
**Draft for scope/review planning; not a revised specification, physics clearance,
or S11c-d completion record.** Numerical work remains paused. This inventory was
made from source text, saved reports and JSON metadata; no scientific operation,
pickle deserialization, production validator or export job was run. Subsequent
fresh Claude/Grok document reviews both requested revision; their exact draft-1
packet and reports are preserved at `2e8ec5a3`. Draft 2 and round-2 reports are
preserved at `75086f8a`: one amendment CLEAR, one NEEDS REVISION, both inventories
NEEDS REVISION. Draft 3 and both round-3 reports are preserved at `134d033d`,
again with one amendment CLEAR and both inventories needing revision. The user
requested that Codex retain authorship; the proposed Claude authoring task never
launched. Draft 4 and round-4 reports are preserved at `116f6b36`: Grok cleared
both documents; Claude requested the bulk-reconstruction/source/check correction
and inventory fixes. The user has now resumed work after the requested pause.
Codex retains authorship under that preference; this fourth substantive fold
records an explicit G4 author-change exception and remains subject to fresh
Claude/Grok review. It grants no scientific clearance or production authority.

**Adopted direction, pending amendment review:** following the consolidated
proposal, the user said on 2026-09-22, “if we've solidified the plan, then let's
continue.” The decision is **Option B**, with A9 first. The [consolidated
amendment draft][amendment] proposes pole deferral, practical acceptance and the
kept pole claims together. It also addresses non-vacuous open-thickness FORM
coverage. User adoption authorizes the two non-author reviews; it is not their
clearance and production remains paused. The tables retain the original
obligations and explicitly identify the amendment's added coverage item A11,
the clarified retained bulk-depth flux obligation A12, the signed-balance versus
additional loss-attribution disposition A13, and the pole deferral. None is
silently marked complete.

## What the decision actually concerns

**The original real-frequency scattering calculation is substantially further
along than the recent row-by-row updates suggested.** At the approved development
input `omega=1`, all four anchoring/density cases already have complete finite
responses, separate retained-grade continuum responses, current bookkeeping and
the selected uniform/profile/coordinate controls. These are scoped computed
results, with published component transcripts. They are not yet the finished,
independently reviewed four-case engine and bindable export.

The unfinished eight material row families belong to the newer **complex-frequency
calculation at `omega=1-0.01i`**. They support the two material frequency responses
and subsequent bounded pole work. They do not block the already computed
real-frequency responses. My previous shorthand that “material responses and
continuum remain” was too broad; it mixed these two stages.

Three other distinctions determine the remaining work:

- `exploratoryAcceptanceV1` already permits practical numerical scattering and a
  separate bounded numerical pole search. It makes general Fredholm/Riesz
  certification, rigorous tail/Abel estimates and global exceptional-locus
  coverage optional unless the corresponding stronger claim is made. It does
  **not** remove the required pole-search output or authorize an empty answer.
- The governing `nonlinearPoleV2` correction does not have established two-leg
  physical-spec clearance in the records inspected here. However, its bounded
  **Lean NP1–NP4 mathematical core has since received CLEAR reviews from Claude
  and Grok**. Those reviews explicitly exclude analytic existence theory and
  the physical S11c application. Neither “the pole mathematics is all unreviewed”
  nor “the physical pole contract is cleared” describes the current state.
- The Mathematica audit of the upstream repairs **already ran**, followed by
  the approved Mathematica pressure-trace repair. Repeating that audit from
  scratch is not the missing stabilization step. Current-version non-author
  review, and eventual full d cross-engine comparison, are separate obligations.

**Updated recommendation:** use the consolidated Option B amendment to defer
the pole branch and make the general weak-order FORM/export the first
construction priority. Package the repaired upstream inputs and existing
real-frequency evidence for their required reviews. Retain the full two-ended
scattering calculation; do not replace it with an unjustified elementary Born
matrix element. Keep production paused until the adopted amendment
clears its two non-author reviews. The original inventory recommended retaining
the numerical pole branch; this follow-up deliberately changes that recommendation.

## Authority, scope and review status

| Authority | What it requires or permits | Review state / consequence |
|---|---|---|
| [Decisions N2, N4–N6][decisions] | Named localized profile; retained weak-order form; separate profile anchoring from coordinate representation. N2 permits the later spec to refine the initial step boundary. | Starting decisions are not a replacement for the later spec. A general dispersion relation and order-unity physical leakage magnitude remain out of scope. |
| [Shared physics v10 §§2–5][s2] | Complete two-ended distorted scattering; separate continuum re-expansion; conversion/survival; profile-conditional pole data and controls. The shortcut to a simple uniform-background matrix element needs computed premises. | Round-10 clearance is recorded at `399a8516`, with [Opus][spec-opus] and [Grok][spec-grok] reports; the pinned v10 includes their post-clearance wording folds. This historical status does not transfer automatically to later corrections or code. |
| [nonlinearPoleV2 §§1–7][pole] | Correctly typed principal parts, multiplicity, source/observation maps, conditional projections and scoped pole claims. | User approved (`7c98b8ee`); the document explicitly disclaims independent review. Whole governing physical-contract clearance is **not established** by the records inspected. Review its interaction with the later acceptance policy before further pole production. |
| [Exploratory acceptance, numerical sufficiency and optional follow-up][acceptance] | Roughly 1% stability of resolved observables, declared absolute resolution, selected independent comparisons and honest unresolved results. General certification is optional follow-up. | User approved 2026-09-18. No separate two-leg clearance of this change in completion criteria was located in the inspected evidence. Include it in the scope reconciliation; do not imply that user approval is a review report. |
| [Retained solver/export contract, items 3–7][retained] | Real-frequency continuum bookkeeping stays separate from complex-frequency pole data. Transparent symbolic exports remain bindable; numerical records retain their input/domain/status. An unsuccessful search is unresolved. | Preserve the retained text and its SHA256 `f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2`. A fixed-frequency table does not fulfill the generic symbolic handoff. |
| [NP1–NP4 fidelity closure][lean-review] | Conditional algebra, finite Laurent contour identities and discriminating examples. | Two independent CLEAR reports for packet `9c7e8e41dd1820647c3dc787ad171043944038c70d37b873699d414f49396e57`. This bounded mathematical work is complete. Physical application hypotheses and the whole pole directive are outside its clearance. |
| [Project review discipline][review-policy] and [spec §7][s7] | Non-author review of governing physics, builds/instruments and records; separate blind Wolfram engine and T7 comparison. | Saved-result validation, two numerical algorithms, two coordinate implementations and two non-author reviewers are four different kinds of evidence. One does not silently discharge another. |

## How to read the inventory and costs

This is a scientific-workstream inventory, not a line-by-line audit of roughly
700 historical commits. A workstream includes its listed sources, helpers,
recovery files, checks and output records. Historical “next” paragraphs in reports
are superseded by later completion checkpoints where available.

**D** = required deliverable. **S** = necessary support for the indicated
deliverable/claim. **O** = optional follow-up or an implementation choice, not a
separate required physical output. “Computed” means the cited saved result exists
within its stated domain; it does not mean independently cleared or complete for
all parameters.

**Review status:** **V/P** means construction/validation evidence exists, but
two current-version non-author physics clearances have not been established in
this inventory. **P** means review remains pending. **CLEAR-scoped** means the
specific identified packet is cleared, within its recorded scope. The first
independent document-review round requested revision; it did not validate the
underlying scientific calculations or clear the present revised document.

Cost bands are rough **additional focused implementation/packaging effort**:
**S: 0.5–2 hours; M: 2–8 hours; L: 1–3 working days; XL: more than 3 working days,
with no reliable ceiling yet.** They are judgment estimates, not measured
forecasts or elapsed-time promises. Reviewer/provider waiting and revisions are
additional; related rows share work and must not be summed mechanically.
Measured worker times below give scale only. Costs assume the present physical
inputs remain valid. A substantive review finding can invalidate that assumption.

### A. Real-frequency deliverables and their supporting construction

| ID / item | Exact clause served | Class | Finished / remaining | Independent review | Rough cost to finish |
|---|---|---|---|---|---|
| A1. Stabilize the four upstream repairs: conservative inertia, mechanical loading, pressure trace, thickness coordinate | [§1b inherited-input status][s1b]; [§6 method][s6]; project G1/G4 | S for every downstream claim | Source repairs, regenerations and focused checks exist. The [Mathematica audit][wl-audit] (`946acc84`) found no matching defect on its audited domains for three issues, and found the pressure-trace defect in the native Mathematica c2 N6 instrument. Its [approved repair][wl-repair] (`5acfdf30`) includes focused and native four-case evidence. The [thickness-coordinate repair report][thickness-repair] records b/c1/c2/d regeneration/export evidence but does not claim a fresh result for its retained two-frequency endpoint pairing. That diagnostic has unestablished completion in this report; locate any later exact-result disposition during A2/A4/A11 input inspection. It is not an automatic second open-thickness frequency mandate. | V/P; neither an author checkpoint nor the audit is full two-leg clearance of every repaired producer/instrument. Original cross-engine operand debt remains explicit. | **L** to assemble exact version/claim packets and a coverage map; review/revision time additional. No blanket audit/regeneration rerun. Actual Wolfram repair production was about 53 minutes, its validation about 10 minutes. |
| A2. Full reduced operator/kernel, named profile functions, dimensions, independent grades, current and end-mode construction | [§1c][s1c], [§2][s2], [§3a][s3a], [§6][s6] | S; derived operators/modes/currents also D | Computed reduction, end modes, native current-normalization records, source/row/boundary operands and all four own-case routes exist. A1's unestablished post-repair endpoint-pairing/full current-adjoint-normalization coverage also applies here; existing records do not establish that diagnostic's completion. Locate its actual later disposition before declaring normalization coverage closed. Explicitly map `K_0`, `K_-`, `K_+` and the surfaced reduced-operator-block versus reduced-kernel off-diagonal residual to their actual emitted operands; blanket reduction status is not their individual coverage. Their conditions and omitted parent orders remain. Final integration/export must consume these objects, not substitute a baseline or fixed-frequency matrix for a pencil. | V/P. Prior spec clearance is not implementation clearance. | **M** to map the consumed dependency/version closure for review, shared with A1 and C1. No new spectrum/current construction justified merely for packaging. |
| A3. Complete two-ended finite response at the approved real input | [§2, especially 379–407][s2], [§3a][s3a], acceptance numerical item 1 | D | **Computed for all four cases at `omega=1`**: each 645 unknowns/four incident columns; baseline reused, three new cases. [Response report/checkpoint][real-response] retains all actual inputs, own ends, residuals and outputs. Positive regulator 0.1 and approximate modal boundaries remain. This is the finite retained-model response, not the separately required continuum expansion. | V/P; separate review and final engine integration remain. | **S** for a compact consumer/result index, shared with C1; no remaining original-input material solve. Original three-case producer took about 687 seconds. |
| A4. Retained continuum grades, amplitude decomposition, conversion/survival/current bookkeeping | [§2 continuum requirement][s2], [§3a][s3a], [§3c][s3c], numerical [§5d bookkeeping][s5d] | D | **Computed for all four cases at `omega=1`** in the [response][real-response] and [current][flux] records: independent grades, finite/continuum responses, baseline/interference/quadratic terms, incident denominators and closed/cross-current contributions. Four open transverse directions and **no open thickness direction** occur at this input. The reported total open-current ratio equals transverse survival here because all open outputs are transverse. Separately, the saved real-momentum bulk-depth domain is empty: `q_depth^2=-k_normal^2-1/25`; depth-integrated interface-normal bulk current is not an escape flux. Closed thickness matching remains present. `C_{T→H}` here means thickness end-channel conversion only; it excludes separately attributed bulk-depth escape and closure dissipation. The distinct generic bulk-flux obligation remains in A12. Required signed exchange/balance evidence and the unavailable separate loss-attribution question are distinguished in A13. A1's unestablished endpoint-pairing coverage also applies to these current denominators and normalization claims; locate the actual disposition before declaring it closed. This establishes neither generic leakage absence nor complete N13 confinement. | V/P; generic representation/export, nonempty-channel coverage A11 and interpretation review remain. | **M** to assemble the weak-order output/claim map, shared with A9/C1. No second continuum solve at the same input. The three-case current producer took about 278 seconds. |
| A5. Eulerian/material-coordinate covariance within each fixed physical anchoring | [§5a coordinate routes][s5a] | D control | **All four selected coordinate controls computed and published**; [coordinate report][coordinate] records finite differences about 2.22e-13 (amplitudes), 2.97e-13 (currents), and forcing/covector omission controls. These do not establish the distinct density-advection mutation in A8. LAB_HELD and MATERIAL_ADVECTED are physically different anchorings. The selected chart regression does not close upstream kernel N6/N3/N4 debt. | V/P; numerical consistency is not an independent derivation of the imported kernel. | **S** to index existing coordinate evidence, shared with review packaging. No new chart sweep. A8's remaining operands are separate. |
| A6. Three uniform-background controls per case | [§5b][s5b] | D control | **Twelve own-background response controls complete**, [publication `44de056f`][uniform], including reflection/identity and incident-sign comparisons. Evidence for the separate reduced-coupling `K_uniform,end−K_uniform,reference` both-operand triplets is **not established by that report**; locating/validating those exact existing operands, or constructing genuinely missing ones, remains. Response identity is not a substitute. Neither establishes exact variable-profile transparent boundaries. | V/P for response controls; coupling-triplet coverage P. | **S** response index; **M–L** conditional on what the coupling operand inspection finds, shared with A2/C1. No repeat homogeneous solve; producer about 53 seconds. |
| A7. Profile FORM, thickness edge/bump and independent modulus discrimination | [§5c][s5c] | D control | [Four-case profile FORM responses][profile] are computed, with reported channel changes below the declared absolute resolutions. That report references a thickness-bump moment; it does **not establish the exact separate thickness and modulus operand/residual payloads** required by §5c: `(w1-prime)_red(0)`, `Delta w1`, their residual, and the analogous modulus pair/residual. Locate the actual saved payloads before constructing anything missing. A Fourier identity as one payload cannot substitute for either pair. A bump moment is not a bump-scattering calculation. | V/P for reported responses; exact discriminator-pair coverage P until joined to actual records. No profile universality follows from sub-resolution changes. | **M** for locate-first operand/output mapping; **L/conditional** only if required controls are genuinely missing, shared with C1. No broad profile sweep or automatic replay. Original three-case response production about 432 seconds. |
| A8. First-jet sensitivity and separate density-advection operands | [§5a explicit mutations and absence operands][s5a], [§3c][s3c] | D control | [Four-case first-jet responses/output][first-jet] are computed, including four finite and three new continuum controls; amplitude/current changes at most 8.243e-7/5.482e-7 are below the stated reporting goals. Their discriminating power for a resolved channel-level claim is not established. The one-sided closed-operator probe is not a consistent new profile or an isolated advection/c2 closure. **Unestablished exact coverage:** the one-coordinate-route `w1-prime`/`grad w1` baseline/mutated/residual triplet; the explicit one-route density-advection omission/reversal in RHOBR; and the RHO4 structural-absence source/density-gradient operands. Locate saved records first. A valid shape-sensitivity probe need not define a new consistent profile; the gap is the exact required operand/route evidence, not that limitation alone. Existing material-covector forcing controls do not discharge the density operands. | V/P for the numerical probe; exact §5a triplet and density/absence coverage P until source operands and records are joined. | **M–L** to inspect/complete genuinely missing controls, shared with A2/A5/C1. Do not repeat completed output phases (about 71 minutes historically) merely to locate evidence. |
| A9. Generic weak conversion FORM, transparent symbolic weak coefficients and downstream handoff | [§3c–§3d][s3c], [§4][s4], [§5d][s5d], retained contract item 6; amendment §2 | D, first construction priority after gates | Numerical coefficients/bookkeeping, full symbolic reduced sources and profile operands exist. **A complete bindable, casewise downstream export is not finished.** The real-frequency arrays do not replace the required differentiable operator/continuum/weak-coefficient expressions. Exact representation and clause coverage need integration, including non-vacuous end-channel coverage A11 and the separately typed bulk-depth escape FORM A12. | P for final artifact; existing components V/P. | **L**, shared with C1, plus conditional A11/A12 effort; overlap must be assessed before summing estimates. Missing generic grade/FORM construction may require new work. Numeric tables alone do not close this item. |
| A10. Practical numerical reliability of the observables actually reported | [acceptance numerical items 3–5 and stopping rule][acceptance] | S | Baseline [resolution/domain/regulator comparisons][domain] support dominant finite real-frequency observables at the stated goals. They do **not** establish resolution/regulator stability of the retained-grade weak coefficients prioritized in A9. That claim-specific assessment remains, as do relevant other-case comparisons. Small claim-carrying effects need tighter checks. Positive regulator, finite boundaries, unresolved tiny signals and omitted pure-second-order terms stay explicit. | V/P for scoped existing comparisons; weak-coefficient assessment P. Residuals and row spreads are not observable error bounds. | **M–L**, after choosing claims, shared with A9/A11; no fixed sweep count. Baseline domain/regulator set took 170 seconds. New batches need a question, measured pilot and stopping rule. |
| A11. Non-vacuous thickness end-channel FORM coverage | [amendment §3.3][amendment], explicitly added acceptance requirement, not a pre-existing second-frequency mandate | D under the proposed amendment | **Unfinished.** Check `J_H,out/J_T,in` on an admissible physical domain with incoming transverse and nonempty outgoing thickness-like full-end channel space, independently symbolic/analytic or numerical. The `omega=1` zero-rank selector does not pay; bulk depth is also closed at that saved point. Bulk-only radiation at some other point would still not discharge this end-channel check. Bounded inspection must distinguish actual root-closure mechanisms and any regulator dependence; a proven structural absence on a stated domain is recorded separately from an unfound witness, then triggers a scope/domain decision, not an indefinite search or automatic waiver. Cross-use B1's saved end relations/threshold metadata only as bounded admissibility inputs. A mathematical continuation outside a physical domain does not establish this coverage. | P; neither a physical witness nor the required independent generic-domain check is established here. | **M** for bounded operand/domain inspection; **L–XL/uncertain** if new generic construction or a physical frequency response is needed. All 80 baseline nonlocal rows depend on frequency; a new point can require new rows, modes and currents. Full input matches determine reuse. One point is not a cheap-call guarantee or permission for a sweep. |
| A12. Separately typed bulk-depth escape flux and weak FORM | [v10 §3b(i)/N13][s3b], [S11b energy accounting][bulk-accounting], [S11b standing limit][rest-source], [c1 §2b][c1-validity]; amendment §2 | D retained within A9; not part of pole deferral | **Unfinished / actual construction coverage not established.** Source the half-space solution map from c1 §§1b/1d/2a/3 acoustics, radiation and physical-face closure data on the solved slab state, with repaired shifted/tilted-face and reference/physical-pressure conventions. Reduce only in-plane content; keep the depth coordinate. Reconcile new construction/import roots in the reviewed consume set. Emit literal physical-face trace residuals and an independent acoustic-face-power versus bulk-control-surface/far-field flux check with matched measures, limits and any lateral contribution. These are explicit added detailed d construction/check/coverage duties. Existing c1 restricted energy operands may supply genuine reuse only after full input/domain joins; no completed solve is replayed. Retain incident normalization, baseline/interference and independent grades. End `J_H`, depth-integrated normal current and a deficit do not substitute. At `omega=1` bulk radiating support is closed. The [frequency-source report][frequency-source] lists algebraic bulk branches `±sqrt(5)` but does not establish an admissible radiating witness. Reconcile the source-context validity wording under amendment §2 and establish nonempty physical support for an independent check. Bounded saved-operand/domain inspection comes first. If support or validity remains unestablished, stop for a scope/domain decision; record an actual structural absence separately from no witness found. Neither pays nonempty-support coverage or permits a zero export. | P for generic flux, sourced half-space reconstruction/consume-set joins, trace and power checks, validity/support and export. Strict-rest-frame grazing is separately labelled; off-grazing smallness is not a uniform moving-background justification. No physical total-loss/confinement answer is claimed. | **M** bounded inspection; **L–XL/uncertain** if admissible generic construction/checks are missing, shared with A2/A9/C1. No automatic extension, extra frequency or global survey; unresolved admissibility stops rather than accumulating production. Scientific costs remain unmeasured. |
| A13. Signed closure/interface balance and separate loss-attribution boundary | [§3a current identity][s3a], acceptance item 2, [S11b energy-accounting discriminators][bulk-accounting] as method, [§6][s6], amendment §§2/3.1/5 | D/S for required signed-balance/control evidence; O/new scope for a separately attributed absorption observable | Required signed exchanges, their transport meanings and applicable independent balance/sign checks remain current duties under A2/A4/A12. Their complete d-level operand/coverage map is **not established by the cited response summaries**. Separately attributed closure/interface absorption is unconstructed; a balance remainder is not assumed positive or assigned an exclusive mechanism. Additional observable design belongs to **“Closure/interface loss attribution — future scope decision”**, not current production or an S11c-e assignment. This does not defer any already required balance operand/check. The retained end conversion, bulk escape and survival do not bound total N13 loss. | P for the full balance/control coverage map and consumer interpretation; no separate absorption result or clearance claimed. | **M** locate-first balance/control mapping, shared with A2/A4/A12/C1–C3; missing required work is costed there. **0 new separate-attribution production** under this amendment; any new observable requires its own scope decision and estimate. No absent loss is exported as zero. |

### B. Complex-frequency and pole branch — historical work and explicit Option-B disposition

| ID / item | Exact clause served | Class | Finished / remaining | Independent review | Rough cost to finish |
|---|---|---|---|---|---|
| B1. Complex-frequency source, scalar, end and numerical-row preparation | [§3b][s3b]; [pole §1][pole]; retained contract item 5; A2/A11 for shared end-relation metadata | S for bounded frequency/pole work; shared local scattering support remains retained | Own complex input families, coefficients, end continuation and saved-input routes exist. Of **32 initially missing owners at `1-0.01i`, 24 are computed**: nine **LAB_HELD/RHOBR_CONSTANT** owners ([remainder-row report][lab-remainder] and [checkpoint][lab-remainder-cp]) plus fifteen **MATERIAL_ADVECTED/RHO4_CONSTANT 1D** owners ([latest material 1D checkpoint][material1d], accepted at `0d77af53`). Actual full inputs govern aliases; those counts do not establish other-density coverage or numerical reuse. End relations and algebraic threshold metadata are also A11 inspection inputs, not proof of an open physical channel. | V/P. | **0 new pole-support production under B**; historical estimate **M** for remaining eight owners if that branch were retained. Shared A11 inspection is costed there. Existing LAB batch timings do not guarantee material costs. |
| B2. Eight material 2D row owners and two material complex responses | Same as B1; full source/end obligations in [§2–§3a][s2] | S under the unamended pole scope; **deferred under B** | **Deferred to Localized-interface pole study — deferred from S11c-d:** MATERIAL/RHO4 rows 60–63 and MATERIAL/RHOBR rows 50–53, preserving the finite profile remainder, positive-regulator Abel term and all mixed coefficients. Source actions already exist. Then two full own-case 645-by-645/four-incident responses at `1-0.01i`. The [routing draft][pause] is preserved and unlaunched. | P for new results; methods/components V/P. | **0 new production under B**. Historical unamended estimate **M–L** for rows, assembly and focused saved-result checks together, overlapping B1. Benchmark: LAB four 2D rows took 12.88 seconds; LAB response production 31.29 seconds (solve 0.56 seconds). These are not a guaranteed material runtime. |
| B3. Two LAB complex-frequency responses and selected row-to-observable sensitivity | [§3b][s3b]; [acceptance stopping rule][acceptance] | Historical S; preserve under B | [LAB response][lab] complete at `1-0.01i`: baseline whole response reused, own RHOBR response rank 645, condition about 5655, residual 9.38e-16. [Row46 sensitivity][sensitivity] changes all 16 amplitudes by at most 4.67e-14; stop that momentum refinement for this response. These complex amplitudes are not Hermitian flux probabilities. | V/P; independent saved-result checking is not two-leg physical review. | **S** reference packaging under B only; **0 new production**; no new LAB solve or row refinement called for. |
| B4. Bounded numerical frequency-pole searches | [§3b][s3b], [pole §1][pole], [acceptance pole paragraph][acceptance] | D under unamended v10; **deferred under B** | Baseline [16/32-point search][contour] on `abs(omega-(1-0.01i))=0.02` is **closed without a resolved candidate**. It is not certified empty; do not double it automatically. Other case searches are deferred to **Localized-interface pole study — deferred from S11c-d**; no new search is current S11c-d production. A completed response at one complex point is not a search. Chart/sheet, seed path, count, practical precision and stop decision must be specified before new batches. | V/P for existing baseline diagnostic; governing correction review pending. | **0 new search under B**. Historical unamended estimate **L–XL** for remaining scoped diagnostics after B2 and the scope decision; candidate-dependent. Baseline added 16 midpoints in 639.84 seconds. Old multiworker timings do not forecast current serialized 2 GiB runs. No promise of a resolved yes/no. |
| B5. Actual candidate principal parts, overlaps and physical-bound classification | [pole §§2–5][pole-data]; [§3b][s3b] | Conditional D under unamended v10; **deferred under B** | Any future candidate work belongs to **Localized-interface pole study — deferred from S11c-d**. No physical profile-dependent pole is established by the current baseline search. Therefore no actual candidate residue, full principal part or physical-bound disposition is complete. Synthetic examples, bulk determinants and zero-normal-momentum algebraic end candidates are not substitutes. Retain full source/observation maps, multiplicity, sheet/decay/width/all-channel tests whenever such a claim is made. | P for the physical application. Bounded NP1–NP4 core is CLEAR-scoped. | **0 production under B**, even if historical data later suggest a candidate; separate authorization/plan required. Historical cost **XL/unknown if a candidate is found**, especially at a threshold or defective pole. Do not prebuild an exhaustive exceptional-locus campaign. |
| B6. General Fredholm/Riesz realization, rigorous contour/tail/Abel bounds and promotion to the parent operator | [pole §§1, 4, 6][pole]; [acceptance optional follow-up][acceptance] | O for numerical results; S only for the corresponding theorem/certified or promoted claim | Conditional mathematical handoffs exist; physical hypotheses are not thereby established. General certification, global exceptional-locus coverage and parent-model pole promotion are not completed. **Defer as a blanket campaign.** Local validity checks needed for an actually claimed numerical candidate remain in B4/B5. | NP1–NP4 CLEAR-scoped; full physical certification P. | **0 in the recommended numerical scope; XL/unbounded if elected separately.** No theorem-strength language without its hypotheses. |

### C. Integration, independent review and final closure

| ID / item | Exact clause served | Class | Finished / remaining | Independent review | Rough cost to finish |
|---|---|---|---|---|---|
| C1. Complete SymPy engine, own-rows delta, fresh keys, recursive bind closure and compact/expanded equivalence | [§7][s7], retained contract items 6–8; subsequent reviewed build reconciliation | D | Many component transcripts and checkpoints are published. The prototype engine exists; **`scripts/S11c_d_exports.py` does not exist at this checkpoint**. Final integration must map every required root to an actual computed expression or honestly scoped numerical record, preserve sources/units/grades and avoid replay of completed scientific work merely for output plumbing. | P for final engine/export; component V/P is not final build clearance. | **L–XL**, including A9 rather than additional to it. Output encoding and legacy helper integration have caused substantial historical overhead; a new all-ancestry audit is not the solution. |
| C2. Non-author review of physical authorities, upstream repairs, d builds/instruments and final record | [project G1/G4 and artifact table][review-policy] | S, required review process | Original v10 has identified historical clearance. The actual build directive/program brief need reconciliation to the amended scope under their own review gate. The consolidated amendment, repaired upstream versions, accumulated d implementation/instruments and final physics-bearing prose need exact-version review accounting. Two non-author reviews of the bounded Lean core are complete and must not be scheduled again as if missing. | Mixed: identified packets CLEAR-scoped; outstanding physical versions P. Round 3 gave one amendment CLEAR and one NEEDS REVISION; both inventories NEEDS REVISION. Draft 4 remains P, with Codex authorship retained at the user's request and fresh Claude/Grok reviews required. | **L** to prepare dependency-scoped packets; external review/folds **XL/uncertain**. Applicable artifact/version clearance remains required; grouping does not waive checks. |
| C3. Blind Wolfram d engine and T7 cross-engine comparison | [§7, comparator and blind engine][s7]; amended retained/deferred coverage | D in whole-step process; outside current SymPy builder lane | The upstream Mathematica audit/repair exists. A complete blind **d** engine/result and completed d T7 comparison were not identified. Retained v10 coverage includes the both-operand reduced-operator/kernel joins and reconstruction from those rows, modal currents, signed balance and end-channel scattering/weak FORM/controls. The amendment adds the c1-sourced half-space solution map, physical-face trace residuals, independent acoustic-power/far-field check and bulk-depth flux/weak FORM to independent build and T7 coverage; these are not already inherited/completed payloads. Keep actual dimensions/grades and propagate c2 operand debt. Under B, pole families are explicitly outside comparison membership, never fabricated zero residuals. | P. A same-code saved-result validator or upstream audit does not close it. | **XL**, no reliable detailed estimate before a scoped independent build plan. A fresh engine is more than replaying an existing audit. |
| C4. Reviewed interpretation, downstream FORM handoff and S11c-d completion record | [§3d][s3d], [§5d][s5d], [§8][s8]; record review policy | D | Concise constituent reports exist. Final generic-vs-numerical scope, unresolved search records, practical resolution, inherited debts and S11c-e handoff still need consolidation after C1–C3 or an explicitly revised contract. No order-unity slit-edge leakage magnitude or physical falsification bound is available from this toy-model step. | P for final record. | **M–L**, plus its own review/folds. Physical magnitude remains R1-blocked; it is not another numerical job to finish here. |
| C5. Saved-input journals, completion hooks, guards, recovery adapters, validation and output inventories | [§6 script obligations][s6], [§7 provenance/output][s7], [local safeguards][agents] | S where needed for trustworthy execution/reuse; exhaustive repetition O | Existing evidence/failure histories are preserved. They support real outputs but are not extra physics deliverables. All future scientific jobs retain the one-worker guard; no old science is replayed for plumbing. Direct references replace recursive copying. Repeated exhaustive unchanged-source/ancestry validation is not the forward plan. | Operational checks exist; a physics-bearing instrument still belongs in C2. Hash identity alone is not physics review. | **S–M** to document the current consumer boundaries. Further machinery only for a concrete missing input/result/check; no open-ended journal/refactor project. |

## What to defer, and what cannot silently be cut

The following require no new S11c-d production under the adopted Option-B direction (effective only after amendment clearance):

- B2's eight material complex rows and two material complex responses, B4's new
  casewise pole searches, and B5's candidate principal-part/overlap/classification
  production are deferred to **Localized-interface pole study — deferred from
  S11c-d**. B3 is preserved/reference-packaged only. Shared local end metadata
  remains available to retained A2/A11; no new pole calculation is implied.

- Another baseline contour doubling, broad profile/frequency grids, global
  exceptional-locus surveys, or intermediate refinement after its effect on the
  reported response is demonstrably below the chosen resolution.
- General infinite-dimensional/Fredholm, tail, Abel-limit and parent-pole
  promotion certificates absent a claim needing them. Existing bounded Lean
  results stay available; no additional theorem campaign is implied.
- Separately attributed closure/interface absorption as an additional observable:
  **“Closure/interface loss attribution — future scope decision.”** This is not
  authorization to omit existing signed exchange/balance operands or checks;
  those remain retained in A13 and its supporting rows. No total N13 loss bound
  is inferred from the incomplete attribution.
- Reconstructing unsaved old symbolic internals, repeating completed science or
  copying/revalidating the entire ancestral inventory to create a new handoff.
- Order-unity strong-edge dynamics, a generic global dispersion relation, the
  R1-dependent physical leakage magnitude, and the explicitly deferred giant
  upstream cross-engine operand families. These are outside this step's current
  boundary, not hidden remaining rows.

The complete two-ended response, retained continuum bookkeeping, specified
controls, bindable export, applicable independent review and d cross-engine
process **cannot** be dropped while calling the result full compliance with the
unamended v10. Under that unamended authority, the casewise bounded pole-search
output is required even though its general-certification layer is optional.
Option B changes that membership explicitly through the reviewed amendment;
its deferred pole work is not current S11c-d production. The end-channel and
bulk-depth continuum FORM obligations remain retained and separately typed.

## The decision options

| Option | Deliverable and next work | Cost / contract consequence |
|---|---|---|
| **A. Retain the practical full S11c-d contract** | Preserve the already computed four-case real response/continuum/controls; resolve and independently review the pole/acceptance wording and upstream repair packets; package the existing outputs. Then finish only the eight material complex rows, two material responses and selected remaining bounded searches. Candidate work is conditional. Complete final engine/export and separate reviews/Wolfram/T7. | Numerical gaps appear modest relative to integration and review, but whole-step completion is still substantial. **Several focused working days plus independent review/engine cycles**, with no defensible fixed finish date. No new global certification campaign. |
| **B. Scattering/FORM handoff with pole deferral — adopted direction, amendment pending review** | Finish the retained scattering/continuum/FORM deliverable; defer remaining casewise pole work to its named later package. Preserve existing results. Complete the locate-first A6–A8 control coverage, A9, new A11 nonempty-thickness-current coverage, clarified retained A12 bulk-depth escape FORM and required A13 balance evidence/claim disposition. A11/A12 unresolved admissibility stops for a scoped decision; no forced witness or automatically added separate absorption calculation. | Saves unfinished production whose purpose is B1–B5 pole work, not A9–A12/C1–C4 or reviews. A11 can require frequency-dependent rows, modes and currents and can cost comparably to a response stage; A12 also has unresolved construction cost. Net savings are not guaranteed. Requires the **reviewed amendment** before production. No full original-scope completion claim. |
| **C. Replace full scattering by a simple closed-form Born element** | Only legitimate under the computed shortcut premises in §2, or through an explicitly reviewed change of the physical problem. Existing full-scattering results should inform that assessment. | Not an automatic cost reduction or currently justified substitute. A new premise audit/spec revision could cost more than packaging existing results. Not recommended as the default. |

The user adopted Option B and three document-review rounds have completed.
After the third round the user requested that Codex retain drafting. The
external authoring proposal never launched. This third substantive Codex fold
uses fresh Claude/Grok non-author reviews; the inventory and amendment remain
pending clearance and do not authorize production or restart the preserved
material-row draft.
Their adopted direction, review clearance and later build authorization remain
distinct records.

## Evidence boundary and preservation

The accompanying [evidence index][evidence] records selected authority/report/
checkpoint paths, current hashes and JSON status fields, review entry points,
the paused draft and the absent final export/engine paths. It is a navigation
index, not a new scientific validation or a claim that every historical commit
has been reviewed. In particular, formal clearance must come from the actual
versioned review packet and reports, not an `ACCEPTED_*` checkpoint label.

Source files, accepted outputs, failure histories and the retained builder suffix
were not changed. The working tree was clean at `0d77af53` before this draft;
no relevant scientific worker was running in the inspected process snapshot.
The exact reviewed first draft was committed for preservation at `2e8ec5a3`,
explicitly without clearance. Draft 2 and its reports are preserved at `75086f8a`,
and draft 3 and its reports at `134d033d`. The present revision and evidence
index remain pending document review; neither reruns nor validates the science.

[decisions]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_decisions.md:47
[s1b]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:132
[s1c]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:183
[s2]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:332
[s3a]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:421
[s3b]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:479
[s3c]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:545
[s3d]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:631
[s4]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:681
[s5a]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:711
[s5b]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:766
[s5c]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:791
[s5d]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:813
[s6]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:827
[s7]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:857
[s8]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SHARED_PHYSICS.md:909
[pole]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_NONLINEAR_POLE_CONTRACT.md:3
[pole-data]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_NONLINEAR_POLE_CONTRACT.md:44
[acceptance]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_EXPLORATORY_ACCEPTANCE.md:8
[retained]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_sympy_builder_report.md:488
[review-policy]: /var/projects/toy_physics/CLAUDE.md:140
[spec-opus]: /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_d_shared_physics_review_r10_opus.md
[spec-grok]: /var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_d_shared_physics_review_r10_grok.md
[lean-review]: /var/projects/toy_physics/research/pde_ledger_v3/lean/s11/POLE_FIDELITY_REVIEW.md
[wl-audit]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_wolfram_repair_audit_report.md
[wl-repair]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_wolfram_pressure_trace_repair_report.md
[real-response]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_response_report.md
[flux]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_flux_report.md
[coordinate]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_coordinate_response_report.md:274
[uniform]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_uniform_report.md
[profile]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_profile_response_report.md
[first-jet]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_first_jet_response_report.md
[domain]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_finite_scattering_domain_report.md
[material1d]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_frequency_material_rows_1d_checkpoint.json
[pause]: /var/projects/toy_physics/_scratch/s11c/s11c-remaining-case-frequency-20260921/material-2d-routes/paused-for-scope-assessment.json
[lab]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_frequency_lab_operator_recover_report.md
[sensitivity]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_frequency_lab_row_sensitivity_report.md
[contour]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_frequency_contour_refine_report.md
[agents]: /var/projects/toy_physics/AGENTS.md
[evidence]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_scope_inventory_20260922_evidence.json
[amendment]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md
[thickness-repair]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_thickness_coordinate_repair_report.md:387
[bulk-accounting]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md:452
[rest-source]: /var/projects/toy_physics/research/pde_ledger_v3/steps/S11b_interface_coupling_law.md:155
[c1-validity]: /var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:204
[frequency-source]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_frequency_source_report.md
[lab-remainder]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_frequency_remainder_rows_report.md
[lab-remainder-cp]: /var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_remaining_case_frequency_remainder_rows_checkpoint.json
