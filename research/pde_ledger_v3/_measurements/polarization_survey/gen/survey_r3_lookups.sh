#!/bin/bash
# Mechanical fact-lookups (sha256sum/grep/sed/cat only) behind the polarization survey review, round 3.
cd /var/projects/toy_physics || exit 1
run() { echo '```'; echo "\$ $*"; eval "$@" 2>&1 | cut -c1-600; echo '```'; echo; }
V=_scratch/polarization/POLARIZATION_SURVEY.md
P=_scratch/polarization
L=research/pde_ledger_v3
echo "# Lookups — polarization survey review, round 3 (generated $(date '+%Y-%m-%d %H:%M'))"
echo
echo "Generator: \`_scratch/polarization/gen/survey_r3_lookups.sh\` (sha256sum/grep/sed/cat only)."
echo
run "(cd $P && sha256sum -c survey_review_baseline_r3.sha256)"
run "grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*\|\*\*Verdict: [A-Z][A-Z ]*' $P/survey_review_r3_claude.txt $P/survey_review_r3_grok.txt"
echo "## Finding 1: E01's class"
run "grep -n '| E01 ' $V | sed -n '\$p'"
run "grep -n 'R-S1-01' $V | head -4"
run "sed -n 134,136p $L/SUBSTRATE_REQUIREMENTS.md"
run "sed -n 374,376p $L/V3_STEP_PLAN.md"
run "sed -n 182,183p $L/steps/S10_two_transverse_photons.md"
run "sed -n 73,76p $L/steps/S10_two_transverse_photons.md"
run "grep -n 'C1' $P/survey_r0_review_disposition.md | head -3"
run "grep -n 'Reproduced\|supplied input that states it\|E01 and E22' $P/survey_repair3_fresh_author_prompt.md"
echo "## Finding 2: where the parity condition is located"
run "grep -n -c 'S11bB_interface_assembly.md:53' $V"
run "grep -n -o 'S11bB_interface_assembly.md:53[^;]*' $V | head -4"
run "sed -n '17p;53,55p;76p' $L/steps/S11bB_interface_assembly.md"
run "sed -n '284p;288p' $L/directives/S11b_SHARED_PHYSICS.md"
run "sed -n 358p $L/directives/S11bB_SHARED_PHYSICS.md"
run "grep -c 'S11bB\?_SHARED_PHYSICS' $V"
run "sed -n 110,111p $L/steps/O2_steady_brane_balance.md"
run "cat $P/survey_review_r3_claude_evidence/chiral_passive_check_stdout.txt"
echo "## Cause of finding 2: the round-1 disposition and repair-2 brief listed S11bB:53–54 among the parity sources"
run "grep -n -o 'S11bB:53–54[^;|]*' $P/survey_r1_review_disposition.md"
run "grep -n 'S11bB_interface_assembly.md:53' $P/survey_repair2_prompt.md"
