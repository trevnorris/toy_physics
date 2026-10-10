#!/bin/bash
# Mechanical fact-lookups (sha256sum/grep/sed only) behind the polarization survey review, round 2.
cd /var/projects/toy_physics || exit 1
run() { echo '```'; echo "\$ $*"; eval "$@" 2>&1 | cut -c1-420; echo '```'; echo; }
V=_scratch/polarization/POLARIZATION_SURVEY.md
L=research/pde_ledger_v3
echo "# Lookups — polarization survey review, round 2 (generated $(date '+%Y-%m-%d %H:%M'))"
echo
echo "Generator: \`_scratch/polarization/gen/survey_r2_lookups.sh\` (sha256sum/grep/sed only)."
echo
run "(cd _scratch/polarization && sha256sum -c survey_review_baseline_r2.sha256)"
run "grep -n -o 'Verdict:\*\* [A-Z][A-Z ]*\|Verdict: [A-Z][A-Z ]*' _scratch/polarization/survey_review_r2_claude.txt _scratch/polarization/survey_review_r2_grok.txt"
echo "## Finding 1: E12 and the rule-5 wording"
run "grep -n 'no model source bears on it' _scratch/polarization/survey_repair2_prompt.md"
run "grep -n -o 'no source bears on that response' $V"
run "grep -n '| E1[0-2] ' $V | sed -n '4,6p' | cut -c1-260"
run "sed -n 1183p $L/V3_STEP_PLAN.md"
echo "## Finding 2: how the survey reads S11:52–64"
run "grep -n -o 'S11_stray_longitudinal.md:52–64[^.]*' $V | head -3"
run "sed -n 52,55p $L/steps/S11_stray_longitudinal.md"
run "grep -c 'a defect splits them' $V"
echo "## Finding 3: N_eff"
run "grep -c -i 'N_eff\|N_{\\\\rm eff}\|N_{eff}\|effective number of' $V"
echo "## Finding 4: rules 3 and 4, E01 and E22"
run "grep -n 'A fact that follows only from a supplied input' _scratch/polarization/survey_repair2_prompt.md"
run "grep -n -o 'not a computed tomography/physical-state-selection object for rules 1 or 3' $V"
run "sed -n 73,76p $L/steps/S10_two_transverse_photons.md"
echo "## Finding 5: missed sources"
run "sed -n 74,76p $L/steps/S11bB_interface_assembly.md"
run "sed -n 425p $L/SUBSTRATE_REQUIREMENTS.md"
run "grep -c 'S11bB_interface_assembly.md:74\|S11bB:74' $V"
run "grep -c 'SUBSTRATE_REQUIREMENTS.md:425\|R-S8-06' $V"
