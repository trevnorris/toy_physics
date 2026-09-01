# Script-catalog task

Apply `00-snapshot-contract.md`. This task refreshes the declared `entries`
region of one mixed, domain-grouped script catalog.

## Inputs

The transaction task supplies the domain, staged base page, entry-point unit
records, prepared semantic/excerpt content, identity-only metadata, related
capsules/topics, and target commit. Catalog only the meaningful entry points
selected by the task. Do not enumerate every helper, test, generated module,
cache, or output file.

## Entry extraction

For each entry point or coordinated family, write:

- repository path and language;
- configured source-unit ID and role: builder, validator, comparator, runner,
  reporter, or visualizer;
- named objects it computes/builds/compares/emits, without a scientific
  conclusion;
- inputs, arguments, fixtures, premises, or imported exports actually present;
- stdout tags and tracked result/output paths visible in the packet;
- guards, invariants, controls, and their exact coverage;
- exact recorded invocation, or `not recorded`;
- separate interpretive paper/step/result source and anchor, or `none`;
- related topic/capsule links supplied by the task;
- target commit and member blob identities.

Never state that a script proves, establishes, validates, or shows a physical
claim merely because it asserts, exits successfully, or prints `PASS`. Use
`measured` only when the packet carries the exact recorded invocation, literal
stdout/transcript with the relevant operands and residual, and separate
interpretation. Otherwise use `provisional` for the catalog description and
name missing invocation, output, or interpretation. For identity-only
artifacts, record identity and path only.

When a runner coordinates engines, comparator, and report, one family entry is
preferred over repetitive helper entries. Still identify which member supplies
each role and whether engines share premises or exports. A visualizer's render
is not a verification unless a supplied source explicitly records that role.

## Output contract

Start from the staged base page and replace only the exact
`BEGIN GENERATED:entries` region. Preserve markers and all unmarked body
bytes. Generated statements and the page roll-up remain
`memory_review: ai_draft`. Update only delegated machine-owned frontmatter;
sort and deduplicate direct original source paths.

Use one `## \`path\`` heading per entry/family and the field order in schema
version 2. End each entry with member-specific `Sources` citations. Do not cite
transaction paths or memory capsules as evidence.
Obey the catalog budget (default ceilings: 250 words per entry, 10 entries, and
2,500 words for the generated region). Every substantive behavior or result
statement uses the shared status block; neutral path/language/member metadata
need not.
