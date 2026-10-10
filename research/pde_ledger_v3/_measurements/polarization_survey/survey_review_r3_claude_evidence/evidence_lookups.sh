#!/bin/bash
cd /var/projects/toy_physics
S=_scratch/polarization/POLARIZATION_SURVEY.md
echo '### 1. survey E01 rationale';  grep -n '^| E01 |' $S | grep -oE 'Its action, field content and dimensional input do \*\*not themselves state the two-direction result\*\*; supplying D=3 is not supplying D−1=2'
echo '### 2. survey Part 3 R-S1-01 row'; sed -n '157p' $S | grep -oE 'R-S1-01: [^|]*'
echo '### 3. register R-S1-01'; nl -ba research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | sed -n '134,136p'
echo '### 4. plan S10 entry'; nl -ba research/pde_ledger_v3/V3_STEP_PLAN.md | sed -n '374,376p;387,388p'
echo '### 5. S10 supplied field content'; nl -ba research/pde_ledger_v3/steps/S10_two_transverse_photons.md | sed -n '182,183p'
echo '### 6. survey cites of plan 374-391 in body'; sed -n '1,254p' $S | grep -cE 'V3_STEP_PLAN\.md:3(7[4-9]|8[0-9]|9[01])'
echo '### 7. S11bB odd couplings and headline'; nl -ba research/pde_ledger_v3/steps/S11bB_interface_assembly.md | sed -n '16,17p;53,55p;59p;76p'
echo '### 8. survey E10 parity sentence'; grep -n '^| E10 |' $S | grep -oE 'The unresolved parity/chirality condition lies \*\*outside that class\*\*: [^[]*'
echo '### 9. survey Part 3 row title'; grep -n 'Odd/chiral interface couplings' $S | cut -c1-90
echo '### 10. supplied symmetry group'; nl -ba research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md | sed -n '280,284p;288p'; nl -ba research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md | sed -n '358p;362p'
echo '### 11. survey cites of S11b/S11bB specs'; grep -cE 'S11bB?_SHARED_PHYSICS' $S
