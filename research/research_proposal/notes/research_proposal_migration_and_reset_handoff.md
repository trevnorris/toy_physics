# Migration and Reset Handoff for the Lean Toy-Model Research Proposal

## Purpose

This handoff resets the research proposal after the original skeleton expanded to forty chapters and the first ten completed chapters grew into a book-length combination of proposal, ontology restatement, research-governance manual, status report, derivation-ledger specification, verification SOP, and project-management system.

The earlier work remains valuable. It contains carefully developed distinctions, failure conditions, cross-sector risks, and useful language. It is now treated as a **reference library**, not as a manuscript that must be preserved or assembled wholesale.

The new controlling structure is:

- `research_proposal.tex`

This handoff controls how the old material is used:

- `research_proposal_migration_and_reset_handoff.md`

Together, these two files supersede the old forty-chapter proposal skeleton as the structural authority for the proposal.

---

# 1. Reset decision

The proposal will answer five primary questions:

1. What is the model, in plain English?
2. What must light and each far-field force sector reproduce?
3. What calculations will determine the requirements imposed by each sector on the shared brane--bulk medium?
4. Can one frozen parent medium satisfy the full intersection of those requirements?
5. What throat, particle, interaction, and reservoir work follows if the common medium survives?

Everything else is supporting infrastructure.

The final proposal will therefore contain:

- front matter;
- eight numbered chapters;
- a short companion-document list;
- no default technical appendices.

The target is approximately **23,000--34,000 words**, excluding equations, references, and short tables. Depending on mathematical density and figures, this will likely produce roughly **55--85 manuscript pages**. That is still substantial, but it is a research proposal rather than a technical monograph.

The proposal is not intended to absorb:

- the full canonical ontology and closure ledger;
- a dated audit of current derivation status;
- the complete light L0--L11 program;
- the complete opposite-orientation contact derivation plan;
- the full cross-sector reviewer memorandum;
- the complete derivation ledger;
- detailed CAS or numerical SOPs;
- reusable templates, registries, or revision-control systems.

Those artifacts remain beside the proposal.

---

# 2. Recommended project setup

Loading the full source set and all earlier proposal chapters into one project is appropriate. It gives later drafting sessions access to useful formulations without requiring the files to be repeatedly attached.

The project should distinguish **controlling files**, **scientific source files**, **status/evidence files**, and **reference-only drafting material**.

## 2.1 Controlling proposal files

These govern the new manuscript:

1. `research_proposal.tex`
2. `research_proposal_migration_and_reset_handoff.md`

When either conflicts with the old forty-chapter skeleton, these files control the proposal structure.

## 2.2 Canonical scientific source

The canonical ontology and closure ledger controls:

- ontology;
- definitions;
- symbol meanings;
- epistemic status;
- conceptual freeze rule;
- failure conditions;
- dependency and closure order;
- provenance claims recorded in that ledger.

Current source file:

- `toy_model_ontology_summary(20260819-194322).md`

A later controlled ontology version supersedes it.

## 2.3 Detailed companion research plans

These remain the technical companions for their domains:

- `light_guided_photon_soliton_research_plan(3).md`
- `opposite_orientation_throat_coincidence(6).md`
- `cross_sector_research_burdens_and_compatibility_gates(20260819).md` or the exact filename used in the project

The proposal summarizes and organizes these documents. It does not reproduce them in full.

## 2.4 Current status and executable evidence

The current state of derivations is controlled by:

- the derivation ledger;
- repository artifacts;
- SymPy and Mathematica scripts;
- numerical outputs;
- an independent Codex repository audit.

The proposal will not contain a dated `Current State of the Program` chapter. A separate audit may be maintained as, for example:

- `toy_model_program_status_audit.md`

That report may be updated without rewriting the proposal.

## 2.5 Reference-only proposal material

The following files remain useful sources but are **not controlling assembly inputs**:

- `toy_model_research_proposal_skeleton.tex`
- `chapter_00_document_control.tex`
- `chapter_00a_provisional_abstract.tex`
- `chapter_01_executive_summary.tex`
- `chapter_02_plain_english_ontology.tex`
- `chapter_03_central_research_question_hypothesis_falsifiability.tex`
- `chapter_04_scope_epistemic_status_and_claims_not_being_made.tex`
- `chapter_05_current_state_of_the_program.tex`
- `chapter_06_capability_first_multi_pass_research_strategy.tex`
- `chapter_07_the_two_orders_of_work.tex`
- `chapter_08_frozen_parameters_provenance_and_independent_verification.tex`
- `chapter_09_standard_template_for_every_sector_chapter.tex`
- `chapter_10_acceptance_budgets_and_known_architecture_gates.tex`

Any additional section fragments produced from the old skeleton are also reference-only unless deliberately migrated into one of the new eight chapters.

## 2.6 Recommended project instruction

Use the following as a persistent project instruction or as the opening instruction in the reset drafting thread:

> The lean skeleton and migration handoff are the controlling proposal structure. The canonical ontology controls definitions and epistemic status. The companion plans control detailed light, contact, and cross-sector research programs. The old forty-chapter skeleton and Chapters 0--10 are reference-only sources to distill, not prose or structure to preserve. Current calculation status comes from the derivation ledger and independent Codex audit, not from the proposal. Do not add chapters, large appendices, or new governance systems without explicit approval.

---

# 3. Authority and conflict rules

Use the following precedence order.

## 3.1 Scientific definitions and status

1. Current canonical ontology and closure ledger
2. Later explicit user corrections
3. Detailed companion documents within their stated scope
4. Research proposal
5. Earlier proposal chapters and skeleton notes

The proposal may simplify language, but it may not silently alter the ontology.

## 3.2 Proposal structure and length

1. Lean skeleton
2. This migration handoff
3. Explicit later user instruction
4. Old proposal skeleton only as a source of ideas

## 3.3 Current derivation status

1. Independent Codex repository audit and current derivation ledger
2. Current executable artifacts
3. Canonical provenance table, within its dated limits
4. Old `Current State of the Program` chapter only as historical reference

## 3.4 Conflict treatment

When a source conflict is found:

- do not reconcile it silently;
- preserve the canonical definition;
- record the conflict for the user or Codex;
- treat the old passage as superseded until reviewed;
- do not use the conflict as an excuse to expand the proposal.

---

# 4. Lean proposal structure

## Front matter

### Document Control

**Form:** one-page table embedded in the master document.

**Contains:**

- document role;
- canonical source hierarchy;
- version and date;
- statement that old proposal chapters are reference-only;
- statement that current status comes from the derivation ledger and Codex audit;
- eight-chapter structure control.

**Does not contain:**

- a chapter-length governance system;
- requirement-ID rules;
- discrepancy workflows;
- release-control procedures;
- full source registers;
- detailed reopen protocols.

The full old Document Control chapter may remain as a separate internal project-governance artifact.

### Provisional Abstract

**Target:** 300--450 words.

Summarize:

- the one-medium ontology;
- the central common-medium question;
- the far-field-first program;
- the Candidate Parent System phase;
- the later throat and reservoir phase;
- negative results as legitimate outcomes.

Do not list every chapter or every verification artifact.

---

## Chapter 1: Executive Summary

**Target:** 1,000--1,500 words.

**Purpose:** Give the full proposal in compact form.

**Required content:**

- one-medium brane--bulk ontology;
- light and particles in the ontology;
- central nonempty-intersection test;
- far-field discovery order;
- common-medium freeze and rederivation;
- later supported-throat and global-closure work;
- meaning of full, partial, and negative outcomes.

**Primary old source material:**

- `chapter_01_executive_summary.tex`
- `chapter_00a_provisional_abstract.tex`
- compact mental picture from the canonical ontology

**Migration rule:** Rewrite from scratch. Preserve the clear conceptual sequence, but remove detailed governance, evidence taxonomy, and long lists of open tasks.

---

## Chapter 2: Plain-English Ontology

**Target:** 3,000--4,500 words.

**Purpose:** Give a readable and technically faithful mental picture of the model.

**Required sections:**

1. Arena and one fundamental medium
2. Ordered brane and de-structured bulk
3. Light in the brane
4. Particles as finite supported throats
5. Charge, gravity, magnetism, and radiation
6. Inertia and the complete dressed particle
7. Global conservation, return, and reservoir closure
8. What is postulated and what remains to be shown

**Primary source material:**

- canonical ontology and closure ledger
- `chapter_02_plain_english_ontology.tex`

**Keep and condense:**

- one substance in two material states;
- finite-thickness shear-supporting brane;
- light as transverse brane shear;
- particles as finite oriented throats;
- support mode as structural rather than automatically observed mass;
- charge orientation versus species;
- gravity as source-side transport/stress/conversion/return response;
- magnetism and radiation from moving and accelerating oriented throats;
- active, passive, and inertial roles as distinct;
- global drain--return and reservoir closure.

**Leave outside the proposal:**

- full variable dictionary;
- extensive flux notation;
- complete operational definitions;
- long equation chains;
- full failure-condition list;
- provenance table;
- detailed closure checklist;
- detailed pair-contact set definitions.

Use only the equations necessary to orient the reader. The canonical ontology remains the place for formal definitions.

---

## Chapter 3: Research Question and Strategy

**Target:** 2,000--3,000 words.

**Purpose:** Explain how the program tests the model without turning the proposal into a governance manual.

**Required sections:**

1. Central research question and scope
2. Capability-first, far-field-first reasoning
3. Requirement-discovery order and mathematical closure order
4. Compact sector-analysis pattern
5. Frozen parameters and proposal--ledger cross-check
6. Falsifiability and levels of outcome

**Primary old source material:**

- `chapter_03_central_research_question_hypothesis_falsifiability.tex`
- selected material from `chapter_04_scope_epistemic_status_and_claims_not_being_made.tex`
- `chapter_06_capability_first_multi_pass_research_strategy.tex`
- `chapter_07_the_two_orders_of_work.tex`
- selected material from `chapter_08_frozen_parameters_provenance_and_independent_verification.tex`
- the compact pattern hidden inside `chapter_09_standard_template_for_every_sector_chapter.tex`
- a concise statement of gates from `chapter_10_acceptance_budgets_and_known_architecture_gates.tex`

**Keep:**

- existential common-medium question;
- separate requirement-discovery and mathematical closure orders;
- far-field capability contracts before Candidate V1;
- provisional effective sources before solved throats;
- one frozen parent system and shared parameter point;
- proposal as top-down contract and ledger as bottom-up evidence;
- independent SymPy and Mathematica checks;
- distinction between sector failure, candidate failure, and architecture failure.

**Condense sharply:**

- six epistemic labels: one short paragraph or compact table;
- claims not being made: a short scope paragraph;
- sector template: one compact sequence, not fourteen modules;
- freeze rule: one concise statement plus change propagation;
- evidence expectations: summary only.

**Leave outside:**

- extensive outcome taxonomies;
- sentence-level drafting rubric;
- claim-envelope schema;
- requirement-ID governance;
- discrepancy-classification system;
- complete CAS SOP;
- numerical method-specific SOP;
- parameter-registry design;
- dependency-graph specification;
- reusable work-package checklist;
- acceptance-budget database.

### Compact sector pattern to retain

Every scientific sector should answer:

1. What behavior must a brane observer see?
2. What material mechanism is proposed?
3. What source and response calculation will test it?
4. What static and dynamic behavior is required?
5. What unwanted channels must be calculated?
6. What does the result require or forbid in the common medium?
7. What outcome passes, conditionally survives, or rejects the branch?

That is enough. The template does not need its own chapter.

---

## Chapter 4: Light Research Program

**Target:** 3,000--4,500 words.

**Purpose:** State why light leads and summarize the two light tracks.

**Required sections:**

1. Research goals and material picture
2. Track A: ordinary guided transverse light
3. Track B: localized classical photon-packet candidate
4. Calculations, unwanted channels, and early failure tests
5. Requirements handed to the common-medium synthesis

**Primary scientific source:**

- `light_guided_photon_soliton_research_plan(3).md`

**Supporting sources:**

- light and longitudinal sections of the canonical ontology;
- concise architecture-gate discussion from the cross-sector review;
- any verified derivation-ledger outputs supplied by Codex.

**Summarize rather than reproduce:**

- two transverse polarizations;
- guided normal profile and confinement/leakage question;
- complete longitudinal and bulk continuum audit;
- finite-thickness dispersion;
- even helper-field hypothesis;
- three-dimensional packet localization;
- nonlinearity versus ordinary vacuum linearity;
- collision and radiation tests;
- handoff of shared coefficients and profiles.

**Do not reproduce:**

- the full L0--L11 sequence;
- all spectral diagnostics;
- all candidate envelope equations;
- the full frozen-light ledger;
- complete artifact lists;
- every individual acceptance criterion and failure condition.

A concise table may list the main outputs required from light. Detailed tests remain in the companion plan.

---

## Chapter 5: Far-Field Force Program

**Target:** 6,000--8,500 words total.

**Purpose:** Walk through the far-field force sectors in the selected discovery order and turn each into requirements for the common medium.

**Required major sections:**

1. Common method for the far-field sectors
2. Gravity
3. Static electricity and charge
4. Magnetism
5. Accelerating charge and electromagnetic radiation
6. Passive response, inertia, extra modes, and drag
7. Far-field requirement summary

### Common section format

Each sector should contain only:

- behavior to reproduce;
- proposed material mechanism;
- calculation path;
- medium requirements;
- principal unwanted modes or departures;
- major failure conditions;
- handoff to Chapter 6.

### Gravity

Use the canonical ontology to summarize:

- source-side material flow, pressure, stress, conversion, return, and reservoir response;
- probe-independent environmental field before passive response;
- inverse-square local regime and attractive sign;
- orientation-even gravity under reflected conditions;
- return screening, overshoot, and matching scale;
- dynamic response poles and propagation;
- one frozen calibration.

Do not insert a full active-mass formalism or all return-profile equations unless essential for the research direction.

### Static electricity and charge

Summarize:

- charge orientation as a finite-thickness boundary condition;
- complete odd/even static and retarded blocks;
- need for an odd eigenbranch with leading `k^2` static stiffness;
- physical interaction sign from the correct branch and force formalism;
- factorization and unit-charge universality;
- oriented count/current and conservation;
- dynamical classification of the electric carrier.

Do not reproduce the full two-interface derivation or contact geometry here.

### Magnetism

Summarize:

- moving oriented throat source;
- sign-dependent transverse response;
- tensor structure and falloff;
- separation from sign-independent inertial wake and velocity-dependent gravity;
- requirement that magnetic normalization follow from the electric source and shared transverse response.

### Accelerating charge and radiation

Summarize:

- one current and coupling across static electricity, magnetism, and transverse radiation;
- energy and momentum transfer to both light polarizations;
- required calculation of radiation into odd, even, longitudinal, bulk, and conversion/reservoir branches;
- electromagnetic closure as a cross-sector test rather than three analogies.

### Passive response, inertia, extra modes, and drag

Summarize:

- distinct active, passive, and inertial roles;
- later requirement for positive ordinary response and low-energy universality;
- full spectrum and characteristic speeds;
- Landau-like thresholds in reversible branches;
- low-frequency friction in dissipative or mixed branches;
- extra radiation and preferred-medium departures.

Detailed calculations remain in the ledger and later throat program.

---

## Chapter 6: Cross-Sector Integration and the Common Medium

**Target:** 3,000--4,500 words.

**Purpose:** Tie the sector programs together and define the transition from requirements to one actual candidate medium.

**Required sections:**

1. Medium Capability Matrix
2. Shared fields, coefficients, speeds, sources, and constraints
3. Major cross-sector compatibility gates
4. Common viable-region test
5. Freezing Candidate Parent System V1
6. Stable brane, complete spectrum, and far-field rederivation
7. Decision outcomes

**Primary sources:**

- `chapter_10_acceptance_budgets_and_known_architecture_gates.tex`, selectively;
- `cross_sector_research_burdens_and_compatibility_gates(20260819).md`;
- old skeleton chapters on far-field synthesis, Candidate V1, stable background, and rederivation;
- canonical closure sequence.

**Keep:**

- one common parameter point or region;
- no weighted average compensating for a failed mandatory capability;
- candidate branch, field inventory, operator basis, and reservoir choice must be explicit;
- stable finite brane and complete constrained spectrum;
- every far-field sector rederived from Candidate V1;
- empty, fragile, and robust common regions as distinct outcomes.

### Seven architecture gates to summarize

Use one concise subsection or table for:

1. broad architecture versus one frozen candidate;
2. dynamic completion of the odd Coulomb sector;
3. characteristic-speed compatibility with light guidance and continuum thresholds;
4. universality, conservation, and discreteness of throat families;
5. photon binding versus ordinary vacuum linearity;
6. persistent drainage, longevity, and reservoir closure;
7. projected coincidence, contact, and annihilation topology.

Do not reproduce the full reviewer memorandum or the entire acceptance-budget framework. Detailed metrics and thresholds belong in the capability matrix and derivation ledger.

### Candidate Parent System V1 declaration

The proposal should say what must be frozen:

- field content;
- conservative, relaxational, or mixed branch;
- equations and constraints;
- operator basis and derivative order;
- symmetries;
- boundary and reservoir data;
- independent parameters;
- source and observable conventions;
- forbidden later additions.

It need not include the full configuration-control schema.

---

## Chapter 7: Beyond the Far Field: Throat and Particle Program

**Target:** 3,000--4,500 words.

**Purpose:** Outline the nonlinear particle work that begins after a common background medium survives.

**Required sections:**

1. One-throat existence, support, and stability
2. Species structure, orientation partners, and charge universality
3. Moving throats, inertia, passive response, and radiation
4. Two-throat interactions and finite-thickness crossover
5. Projected coincidence, core contact, and annihilation
6. Local--global drain--return and reservoir fixed point

**Primary sources:**

- throat and closure sections of the canonical ontology;
- `opposite_orientation_throat_coincidence(6).md`;
- old skeleton chapters on one throat, antiparticle/species, moving throat, two-throat interactions, and global closure.

**Keep:**

- support-mode seed versus self-consistent nonlinear supported throat;
- stability or lifetime;
- orientation reversal versus full antiparticle map;
- species structure and charge universality/discreteness;
- moving-throat inertia, passive response, magnetism, radiation, and drag;
- far-field-to-finite-thickness crossover;
- projected overlap versus full core contact;
- contact versus annihilation;
- same-species antiparticle pair versus different-species neutral pair;
- local--global reservoir fixed point.

**Do not reproduce:**

- the full core-set formalism;
- every graph-regime equation;
- the full Green-matrix program;
- all pair-force ensemble distinctions;
- detailed collision parameter sweeps;
- complete global reservoir ledger.

Those remain in the canonical and focused companion documents.

---

## Chapter 8: Verification, Milestones, and Possible Outcomes

**Target:** 1,800--2,800 words.

**Purpose:** Close the proposal with the evidence process and decision sequence, without duplicating the derivation-ledger handbook.

**Required sections:**

1. Proposal-to-ledger verification
2. Required symbolic and numerical evidence
3. Milestones and decision gates
4. Revision and downstream reopening
5. Interpreting success, partial survival, and failure

**Primary old source material:**

- selected concise material from `chapter_08_frozen_parameters_provenance_and_independent_verification.tex`;
- old skeleton chapters on milestones, reproducibility, and outcomes;
- selected front-matter controls from `chapter_00_document_control.tex`.

**Keep:**

- top-down proposal versus bottom-up ledger;
- independent SymPy and Mathematica implementations;
- numerical evidence where symbolic closure is unavailable;
- assumptions and frozen inputs attached to claims;
- upstream changes reopen real dependencies;
- major gates from sector contracts through global closure;
- outcome hierarchy: sector result, common far-field medium, stable throat, moving/interacting particle sector, global closure.

**Leave outside:**

- detailed artifact directory schemas;
- immutable-output procedures;
- complete discrepancy taxonomy;
- method-specific numerical SOPs;
- full parameter-search governance;
- revision-history appendix;
- work-package review checklist.

---

# 5. Migration map from the old proposal

| Old item | Action | New destination | Notes |
|---|---|---|---|
| Chapter 0: Document Control | Move outside and reduce | Front matter; Chapter 8 only where needed | Keep a one-page table. Preserve the long version as an internal governance artifact. |
| Provisional Abstract | Rewrite | Front matter | Reduce to 300--450 words after body stabilizes. |
| Chapter 1: Executive Summary | Rewrite and condense | New Chapter 1 | Keep conceptual sequence; remove detailed governance and long evidence lists. |
| Chapter 2: Plain-English Ontology | Rewrite and condense | New Chapter 2 | Use canonical ontology as authority; retain mental picture, not full ledger. |
| Chapter 3: Central Question and Falsifiability | Merge | New Chapter 3 | Keep existential test, levels of failure, and outcome meaning. |
| Chapter 4: Scope and Epistemic Status | Compress heavily | New Chapter 3 | One scope subsection and compact status explanation. |
| Chapter 5: Current State of the Program | Remove | Separate Codex audit | No dated status chapter in the proposal. |
| Chapter 6: Capability-First Strategy | Merge | New Chapter 3 | Keep three-pass logic and capability contracts. |
| Chapter 7: Two Orders of Work | Merge | New Chapter 3 | Keep the distinction and concise sequences. |
| Chapter 8: Frozen Parameters and Verification | Split and condense | New Chapters 3 and 8 | Freeze rule in Chapter 3; evidence process in Chapter 8. |
| Chapter 9: Standard Template | Remove as chapter | Half-page pattern in Chapter 3 | Retain detailed template separately if useful for internal drafting. |
| Chapter 10: Acceptance Budgets and Gates | Condense | New Chapter 6 | Summarize gate categories and seven architecture gates. Keep detailed budgets outside. |
| Old Work Package 1: Light | Merge and rewrite | New Chapter 4 | Detailed stages remain in the light companion. |
| Old Work Packages 2--7 | Merge | New Chapter 5 | Use sections rather than separate chapters. |
| Old Far-Field Compatibility Synthesis | Merge | New Chapter 6 | Central capability-matrix section. |
| Old Candidate Parent System V1 | Merge | New Chapter 6 | Keep only freeze contents and decision role. |
| Old Stable Brane and Spectrum | Merge | New Chapter 6 | State required background solve and spectrum. |
| Old Far-Field Rederivation | Merge | New Chapter 6 | One decisive validation stage. |
| Old One-Throat chapter | Merge | New Chapter 7 | Summarize existence/support/stability. |
| Old Antiparticle and Species chapter | Merge | New Chapter 7 | Summarize full map and universality/discreteness. |
| Old Moving-Throat chapter | Merge | New Chapter 7 | Summarize mass response, radiation, and drag. |
| Old Two-Throat chapter | Merge | New Chapter 7 | Preserve contact distinctions; details stay in companion plan. |
| Old Global Closure chapter | Merge | New Chapter 7 | Focus on local--global fixed point; cosmology remains a later extension. |
| Old Milestones chapter | Merge | New Chapter 8 | Use a small gate sequence. |
| Old Reproducibility chapter | Condense | New Chapter 8 | Evidence summary only. |
| Old Outcome Interpretation chapter | Merge | New Chapter 8 | Keep layered outcome meaning. |
| Old Appendices Overview | Remove | Companion list only | No appendix program by default. |
| Glossary and Symbol Dictionary | Separate project artifact | Optional later small appendix | Add only if the final proposal genuinely needs it. |
| Epistemic-Status Rubric | Separate project artifact | None by default | Canonical ontology already controls status. |
| Capability Matrix Template | Separate working artifact | Summarized in Chapter 6 | The filled matrix belongs beside the proposal. |
| Acceptance-Budget Template | Separate working artifact | Summarized in Chapter 6 | Numerical bounds remain in the ledger/matrix. |
| Proposal-to-Ledger Crosswalk | Separate working artifact | Summarized in Chapter 8 | Do not publish the entire schema in the proposal. |
| Frozen-Parameter Table | Separate ledger artifact | Mentioned in Chapters 3 and 8 | Not a proposal appendix. |
| Failure Register | Canonical/ledger artifact | Major failures summarized by sector | Do not repeat the full register. |
| Provenance Table | Canonical/ledger artifact | Mentioned in Chapter 8 | Current status comes from audit. |
| Architecture-Gate Summary | Companion artifact | Concise Chapter 6 summary | Seven gates, a few pages maximum. |
| Work-Package Review Checklist | Internal authoring artifact | None | Do not include in proposal. |
| Revision History | Repository/version-control artifact | Front-matter version only | Add a short revision table only if required for release. |

---

# 6. Hard anti-bloat rules

These rules apply to every rewritten chapter.

## 6.1 Structural limits

- Eight numbered chapters only.
- No new chapter without explicit approval and a revision to the lean skeleton.
- No default appendix chapters.
- A subsection should exist because it advances the scientific proposal, not because a prior document had a section on the topic.

## 6.2 Word limits

Treat the target ranges in the lean skeleton as hard drafting budgets. A chapter may exceed its upper bound only after the user explicitly approves the reason.

Before delivery, report the approximate word count and remove avoidable repetition.

## 6.3 Source delegation

Prefer a short summary and companion reference over reproducing:

- a full derivation plan;
- a complete failure ledger;
- a detailed verification SOP;
- a long formal definition already controlled by the ontology;
- a current-status inventory controlled by Codex.

## 6.4 No duplicated ledgers

The proposal must not become a second:

- ontology ledger;
- provenance ledger;
- parameter registry;
- status audit;
- failure register;
- acceptance-budget database;
- artifact manifest.

## 6.5 Readability rules

- Lead with physical meaning and research direction.
- Use equations only when they clarify the required behavior or decisive test.
- Use one compact table where a table replaces several pages of prose.
- Keep requirement IDs mainly in the capability matrix and derivation ledger, not in every paragraph.
- Summarize major failure conditions; do not reproduce every edge case.
- Distinguish static, dynamic, and nonlinear questions, but do not restate the distinction repeatedly.

## 6.6 Rewrite from scratch

Do not reduce an old chapter by deleting paragraphs sequentially. That tends to preserve the old chapter's oversized logic.

Instead:

1. start from the new chapter purpose and word budget;
2. outline only the questions that belong in the proposal;
3. consult the canonical and companion sources;
4. retrieve useful old formulations selectively;
5. write a clean new chapter;
6. remove anything that belongs in another artifact.

---

# 7. Drafting workflow for each new chapter

## Step 1: Load controlling context

Read:

- the lean skeleton's chapter plan;
- this handoff's chapter-specific migration notes;
- the relevant canonical ontology sections.

## Step 2: Load only relevant reference material

Retrieve only the old chapters and companion sections needed for the current chapter. The fact that all files are in the project does not require all of them to be considered equally in every drafting turn.

## Step 3: Create a word allocation

Before drafting, assign a rough word budget to each section. This prevents one methodological issue from consuming the chapter.

## Step 4: Draft from the scientific question outward

For sector material, use:

1. behavior to reproduce;
2. proposed material interpretation;
3. calculation path;
4. requirements for the medium;
5. major failure channels;
6. downstream handoff.

## Step 5: Perform a delegation edit

For each paragraph, ask:

- Is this already defined canonically elsewhere?
- Is this a detailed derivation better left to a companion plan?
- Is this current status better left to Codex?
- Is this a reusable process rule better left to an internal artifact?

If yes, summarize or remove it.

## Step 6: Perform a repetition edit

Remove repeated explanations of:

- the one-medium test;
- frozen parameters;
- static versus dynamic response;
- current open status;
- the difference between proposal and ledger.

State each centrally, then refer back briefly.

## Step 7: Deliver an input-ready fragment

Each rewritten chapter should:

- use its new `proposal_chapter_...tex` filename;
- contain no preamble;
- contain no `\begin{document}` or `\end{document}`;
- be compatible with the lean master preamble;
- avoid PDF generation unless explicitly requested;
- receive a syntax-only LaTeX check where practical.

## Step 8: Record what was intentionally omitted

The delivery note should identify any major old material deliberately left in a companion or internal artifact. This helps prevent it from being reintroduced later through uncertainty.

---

# 8. New fragment filenames

Use these names so the rewritten proposal cannot accidentally include the old chapters:

- `proposal_provisional_abstract.tex`
- `proposal_chapter_01_executive_summary.tex`
- `proposal_chapter_02_plain_english_ontology.tex`
- `proposal_chapter_03_research_question_and_strategy.tex`
- `proposal_chapter_04_light_research_program.tex`
- `proposal_chapter_05_far_field_force_program.tex`
- `proposal_chapter_06_cross_sector_integration_and_common_medium.tex`
- `proposal_chapter_07_beyond_far_field_throat_particle_program.tex`
- `proposal_chapter_08_verification_milestones_and_outcomes.tex`

The old chapter filenames should remain unchanged as historical/reference artifacts.

---

# 9. Reusable prompt for drafting a chapter in the project

Use or adapt this prompt for each new chapter:

```text
Draft the next chapter of the lean toy-model research proposal as an input-ready
LaTeX fragment.

Controlling structure:
- research_proposal.tex
- research_proposal_migration_and_reset_handoff.md

Scientific authority:
- the current ontology and closure ledger controls definitions and epistemic
  status;
- the relevant companion plan controls detailed research content within its
  scope;
- the derivation ledger and Codex audit control current calculation status.

The old forty-chapter skeleton and old Chapters 0--10 are reference-only. Use
useful ideas or language from them, but do not preserve their structure or
length. Rewrite the chapter from scratch.

Chapter to draft: [INSERT CHAPTER NUMBER AND TITLE]
Target length: [INSERT SKELETON RANGE]
Output filename: [INSERT NEW proposal_chapter_... FILENAME]

Hard constraints:
- do not add chapters or appendices;
- do not create a dated current-status section;
- do not duplicate the canonical ontology, provenance ledger, failure register,
  acceptance-budget database, or detailed verification SOP;
- refer to companion documents instead of reproducing their full derivation
  programs;
- every substantial paragraph must explain the model, define a sector test,
  state a calculation, extract a common-medium requirement, connect sectors, or
  explain the later research path;
- contain no preamble or document environment;
- do not build a PDF;
- perform a syntax-only check where practical.

Before delivery, edit to the word budget and state which detailed material was
left in companion or internal project documents.
```

---

# 10. First rewrite sequence

The cleanest order is:

1. Chapter 2: Plain-English Ontology
2. Chapter 3: Research Question and Strategy
3. Chapter 4: Light Research Program
4. Chapter 5: Far-Field Force Program
5. Chapter 6: Cross-Sector Integration and the Common Medium
6. Chapter 7: Beyond the Far Field
7. Chapter 8: Verification, Milestones, and Possible Outcomes
8. Chapter 1: Executive Summary
9. Provisional Abstract

Writing the body before the Executive Summary and Abstract will reduce later rewrites and prevent those front sections from promising structure that the final body does not contain.

The existing long Executive Summary and Abstract remain useful source material, but the final versions should be written after Chapters 2--8 stabilize.

---

# 11. Completion test for the proposal

The proposal is ready for assembly when a reader can answer the following without consulting the internal governance files:

1. What substance and geometry does the model propose?
2. What are light, particles, charge, gravity, magnetism, radiation, and inertia within that picture?
3. Why does the research begin with far-field capabilities?
4. What will be calculated for each sector?
5. What does each sector require from the common medium?
6. How are conflicting requirements identified?
7. What must Candidate Parent System V1 specify and survive?
8. What throat and reservoir work follows?
9. How will the independent derivation ledger verify the proposal?
10. What results would count as partial success, candidate failure, or architecture failure?

The proposal is still too large if a reader must work through long chapters on document control, status taxonomies, templates, provenance schemas, or acceptance-budget administration before reaching the actual light and force research program.

---

# 12. Reset declaration

The previous expansion was not wasted work. It produced:

- a strong source library;
- carefully articulated failure conditions;
- a clear separation of static, dynamic, periodic, and dissipative problems;
- a mature freeze and verification discipline;
- a set of major cross-sector architecture gates;
- detailed companion plans that can remain outside the main proposal.

The reset changes the **document architecture**, not the scientific standard.

The lean proposal will present the model and research direction clearly. The detailed ledgers, companion derivations, verification procedures, and historical chapters will remain available in the project to support and audit it without turning the proposal itself into a book.
