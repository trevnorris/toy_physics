#!/usr/bin/env bash
# Leg script 07: FORM ablation of the v4 author's potential_trace_pullback_audit.py.
# The ablated copy reads a corrupted copy of the governing spec whose bulk law is made
# rotational (v_bulk = grad4 phi + curl4 A).  Compare the two stdouts.
cd /tmp/s11cd_clean_review_r4_claude/author_ablation
python3 base.py > base.out 2>&1; echo "BASE_EXIT $?"
python3 ablated.py > ablated.out 2>&1; echo "ABLATED_EXIT $?"
echo "DIFF_BEGIN"; diff base.out ablated.out; echo "DIFF_END"
for tag in SUPPLIED_TRACE_VELOCITY SUPPLIED_NORMAL_JET_VELOCITY CLASS_P_ODD_BULK_TRACE_COMPONENT CLASS_P_ODD_BULK_NORMAL_JET_COMPONENT D_VIRTUAL_X3_D_PHYSICAL_U3 VIRTUAL_MAP_IS_PHYSICAL_FRECHET_ROW; do
  echo "TAG $tag BASE=[$(grep "^$tag " base.out)] ABLATED=[$(grep "^$tag " ablated.out)]"
done
grep -n "print(\"CLASS_P_ODD\|print(\"SUPPLIED\|print(\"D_VIRTUAL\|print(\"VIRTUAL_MAP_IS" base.py
