#!/bin/bash
# Mechanical fact-lookups (sha256sum/grep/sed only) behind the polarization survey review, round 1.
cd /var/projects/toy_physics || exit 1
run() { echo '```'; echo "\$ $*"; eval "$@" 2>&1 | cut -c1-300; echo '```'; echo; }
V=_scratch/polarization/POLARIZATION_SURVEY.md
L=research/pde_ledger_v3
echo "# Lookups — polarization survey review, round 1 (generated $(date '+%Y-%m-%d %H:%M'))"
echo
echo "Generator: \`_scratch/polarization/gen/survey_r1_lookups.sh\` (sha256sum/grep/sed only)."
echo
run "(cd _scratch/polarization && sha256sum -c survey_review_baseline_r1.sha256)"
run "grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*' _scratch/polarization/survey_review_r1_claude.txt _scratch/polarization/survey_review_r1_grok.txt"
echo "## The survey's class rule and the rows it is applied to (findings 1, 2)"
run "sed -n 192p $V"
run "grep -n '| E01 \|| E1[0-7] \|| E2[5-9] \|| E3[2-4] ' $V | awk -F'|' 'NF>3 {print \$1\"|\"\$2\"|\"\$3}'"
run "grep -n 'required but with no mechanism\|not addressed by the model\|testable statement' _scratch/polarization/survey_prompt_r0.md"
echo "## Missed v3 sources (finding 3)"
run "sed -n 52,64p $L/steps/S11_stray_longitudinal.md"
run "ls $L/steps | grep -i '^O2'"
run "sed -n 89,92p $L/steps/O2_steady_brane_balance.md"
run "sed -n 110,111p $L/steps/O2_steady_brane_balance.md"
run "sed -n 115,117p $L/steps/O2_steady_brane_balance.md"
run "sed -n 53,54p $L/steps/S11bB_interface_assembly.md"
run "git show ede8aa21:$L/directives/S9b_SHARED_PHYSICS.md | sed -n 117,119p"
run "grep -c 'S11_stray_longitudinal.md:5[2-9]\|S11:5[2-9]\|S11:6[0-4]' $V"
run "grep -c 'S11bB' $V"
echo "## Photon-mass bounds (finding 4)"
run "grep -n -i 'ryutov\|bonetti\|fast radio burst\|FRB' $V | head"
