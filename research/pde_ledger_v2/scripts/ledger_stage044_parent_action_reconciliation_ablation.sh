#!/usr/bin/env bash
# Run both clean stage-044 engines, byte-compare their canonical JSON, then
# demonstrate every registered primitive mutation failing at its own assertion.

set -uo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd -- "${script_dir}/../../.." && pwd)"
py_script="${project_root}/research/pde_ledger_v2/scripts/ledger_stage044_parent_action_reconciliation_sympy_audit.py"
wl_script="${project_root}/research/pde_ledger_v2/mathematica/ledger_stage044_parent_action_reconciliation_mathematica_audit.wl"
scratch_dir="${project_root}/research/pde_ledger_v2/_scratch/stage044"
ablation_log="${scratch_dir}/OUT_stage044_ablation.txt"
py_verdict="${scratch_dir}/verdict_py.json"
wl_verdict="${scratch_dir}/verdict_wl.json"
run_tmp="$(mktemp -d)"

cleanup() {
  rm -rf -- "${run_tmp}"
}
trap cleanup EXIT

mkdir -p -- "${scratch_dir}"
exec > >(tee "${ablation_log}") 2>&1

echo "STAGE044_ABLATION_BEGIN"
echo "PROJECT_ROOT=${project_root}"
echo "PY_SCRIPT=${py_script}"
echo "WL_SCRIPT=${wl_script}"

if ! timeout 600 python3 "${py_script}" > "${run_tmp}/base_py.out" 2>&1; then
  echo "BASE_PY=FAIL"
  sed -n '1,240p' "${run_tmp}/base_py.out"
  exit 1
fi
echo "BASE_PY=PASS"

if ! math -script "${wl_script}" > "${run_tmp}/base_wl.out" 2>&1; then
  echo "BASE_WL=FAIL"
  sed -n '1,240p' "${run_tmp}/base_wl.out"
  exit 1
fi
echo "BASE_WL=PASS"

if ! cmp -s -- "${py_verdict}" "${wl_verdict}"; then
  echo "BYTE_IDENTITY=FAIL"
  cmp -l -- "${py_verdict}" "${wl_verdict}" | sed -n '1,40p'
  exit 1
fi
echo "BYTE_IDENTITY=PASS bytes=$(wc -c < "${py_verdict}")"

timeout 600 python3 "${py_script}" --mutation-map > "${run_tmp}/mutation_map.json"
mutation_count="$(jq 'length' "${run_tmp}/mutation_map.json")"
echo "MUTATION_COUNT=${mutation_count}"
if [[ ${mutation_count} -ne 43 ]]; then
  echo "MUTATION_COUNT=FAIL expected=43 actual=${mutation_count}"
  exit 1
fi

passed=0
while IFS=$'\t' read -r mutation expected_tooth; do
  py_out="${run_tmp}/py_${passed}.out"
  wl_out="${run_tmp}/wl_${passed}.out"

  LEDGER_STAGE044_MUTATION="${mutation}" \
    timeout 600 python3 "${py_script}" > "${py_out}" 2>&1
  py_rc=$?
  py_marker="FIRED_AT_OWN_ASSERT=${expected_tooth}"
  if [[ ${py_rc} -eq 0 ]] || ! grep -Fqx -- "${py_marker}" "${py_out}"; then
    echo "ABLATION mutation=${mutation} tooth=${expected_tooth} engine=py rc=${py_rc} fired=NO"
    sed -n '1,240p' "${py_out}"
    exit 1
  fi

  LEDGER_STAGE044_MUTATION="${mutation}" \
    math -script "${wl_script}" > "${wl_out}" 2>&1
  wl_rc=$?
  wl_marker="FIRED_AT_OWN_ASSERT=${expected_tooth}"
  if [[ ${wl_rc} -eq 0 ]] || ! grep -Fqx -- "${wl_marker}" "${wl_out}"; then
    echo "ABLATION mutation=${mutation} tooth=${expected_tooth} engine=wl rc=${wl_rc} fired=NO"
    sed -n '1,240p' "${wl_out}"
    exit 1
  fi

  passed=$((passed + 1))
  echo "ABLATION mutation=${mutation} tooth=${expected_tooth} py_rc=${py_rc} py_fired=YES wl_rc=${wl_rc} wl_fired=YES"
done < <(jq -r 'to_entries[] | [.key,.value] | @tsv' "${run_tmp}/mutation_map.json")

if [[ ${passed} -ne ${mutation_count} ]]; then
  echo "ABLATION_COVERAGE=FAIL passed=${passed} expected=${mutation_count}"
  exit 1
fi

echo "ABLATION_COVERAGE=PASS fired=${passed}/${mutation_count} engines=2"
echo "BYTE_IDENTITY_FINAL=PASS"
echo "STAGE044_ABLATION_END status=GREEN"
