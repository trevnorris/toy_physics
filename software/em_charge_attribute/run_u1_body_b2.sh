#!/usr/bin/env bash
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
stage="$here/reports/u1_body_dynamics_artifacts/stage_b2_0_intake_radiative_contract"
scratch="$here/_scratch/u1_b2_stage0_run"
input="$here/u1_body_b2_stage0_inputs.yaml"
certificate="$here/reports/u1_body_dynamics_artifacts/b1_orchestrator_approval.yaml"
commit="${STARTUP_CONTRACT_COMMIT:-}"
trust_mode="${TRUST_MODE:-}"

usage() {
  echo "usage: STARTUP_CONTRACT_COMMIT=<40hex> TRUST_MODE=traced $0 --stage0" >&2
  echo "   or: ... $0 --resume-initial-full-b1 (after an interrupted initial full-B1 replay)" >&2
  echo "   or: ... $0 --repair-initial-full-b1-liveness2 (after stale partial liveness output is detected)" >&2
  echo "   or: ... $0 --finish-initial-stage0 (only after a fresh replay closure-only halt)" >&2
  echo "   or: ... STAGE0_CONTRACT_DIGEST=<sha256> ENVIRONMENT_IDENTITY_DIGEST=<sha256> $0 --resume" >&2
  echo "   or: ... STAGE0_CONTRACT_DIGEST=<sha256> ENVIRONMENT_IDENTITY_DIGEST=<sha256> $0 --continue-stage0-after-full-replay" >&2
  exit 64
}

[[ "$commit" =~ ^[0-9a-f]{40}$ ]] || usage
[[ "$trust_mode" == traced ]] || usage

shim_hash="$(sha256sum "$here/u1_body_b2_no_network.c" | awk '{print $1}')"
shim="/var/tmp/u1_body_b2_no_network_${shim_hash}.so"
gcc -shared -fPIC -O2 -Wall -Wextra -o "$shim" "$here/u1_body_b2_no_network.c" -ldl
audit_hash="$(sha256sum "$here/u1_body_b2_first_use_audit.c" | awk '{print $1}')"
audit="/var/tmp/u1_body_b2_first_use_audit_${audit_hash}.so"
gcc -shared -fPIC -O2 -Wall -Wextra -o "$audit" "$here/u1_body_b2_first_use_audit.c" -ldl
real_wolframscript="$(command -v wolframscript)"
real_wolfram_kernel="/usr/local/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel"
test -x "$real_wolfram_kernel"
serial_bin="/var/tmp/u1_body_b2_serial_bin"
mkdir -p "$serial_bin"
cp "$here/u1_body_b2_wolframscript_serial.sh" "$serial_bin/wolframscript"
chmod 755 "$serial_bin/wolframscript"
pycache="/var/tmp/u1_b2_pycache"
mkdir -p "$pycache"

bwrap_exec() {
  local -a sandbox=(--unshare-net --ro-bind / / --proc /proc --dev-bind /dev /dev --chdir "$root" --tmpfs /tmp)
  local prod="$here/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production"
  for writable in "$here/_scratch" "$stage" "$prod" "$pycache"; do
    if [[ -d "$writable" ]]; then sandbox+=(--bind "$writable" "$writable"); fi
  done
  for aggregate in "$here/reports/u1_body_dynamics_results.yaml" "$here/reports/u1_body_dynamics.md"; do
    if [[ -f "$aggregate" ]]; then sandbox+=(--bind "$aggregate" "$aggregate"); fi
  done
  # The certificate-pinned B1 replay mutation harness creates and unlinks its
  # three source copies beside the protected engines.  This replay-only nested
  # mount does not change the read-only root or any B2 production invocation.
  if [[ "${B2_B1_REPLAY_MUTATION_WRITES:-0}" == 1 ]]; then sandbox+=(--bind "$here" "$here"); fi
  timeout 600 bwrap "${sandbox[@]}" env B2_READ_ONLY_ROOT_SANDBOX=1 B2_WRITABLE_MOUNT_POLICY="generated stage/scratch outputs, ephemeral tmpfs, merger-owned aggregate files" B2_FIRST_USE_LOG="${B2_FIRST_USE_LOG_OVERRIDE:-$scratch/traces/stage0/isolated.firstuse}" PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX="$pycache" LD_PRELOAD="$audit:$shim" PATH="$serial_bin:$PATH" WOLFRAMSCRIPT_REAL="$real_wolframscript" WOLFRAM_KERNEL_REAL="$real_wolfram_kernel" "$@"
}

snapshot_initial_existence() {
  local trace_dir="$1" snapshot="$1/initial_existence.paths"
  mkdir -p "$trace_dir"
  find "$root" -path "$root/.git" -prune -o -type f -print | sort > "$snapshot"
}

isolated() {
  local first_use_log="${B2_FIRST_USE_LOG_OVERRIDE:-$scratch/traces/stage0/isolated.firstuse}"
  mkdir -p "$(dirname "$first_use_log")"
  B2_FIRST_USE_LOG_OVERRIDE="$first_use_log" bwrap_exec "$@"
}

trace_run() {
  local trace_dir="$1" name="$2" log="$3" first_use writable aggregate; shift 3
  local -a sandbox=(--unshare-net --ro-bind / / --proc /proc --dev-bind /dev /dev --chdir "$root" --tmpfs /tmp)
  mkdir -p "$trace_dir"
  first_use="$trace_dir/$name.firstuse"
  rm -f "$first_use"
  local prod="$here/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production"
  for writable in "$here/_scratch" "$stage" "$prod" "$pycache"; do
    if [[ -d "$writable" ]]; then sandbox+=(--bind "$writable" "$writable"); fi
  done
  for aggregate in "$here/reports/u1_body_dynamics_results.yaml" "$here/reports/u1_body_dynamics.md"; do
    if [[ -f "$aggregate" ]]; then sandbox+=(--bind "$aggregate" "$aggregate"); fi
  done
  if [[ "${B2_B1_REPLAY_MUTATION_WRITES:-0}" == 1 ]]; then sandbox+=(--bind "$here" "$here"); fi
  echo "[B2 trace] $name"
  timeout 600 strace -ff -qq -yy -ttt -s 4096 \
    -e trace=openat,openat2,execve,chdir,fchdir,close,dup,dup2,dup3,clone,clone3,fork,vfork,connect,sendto,sendmsg,sendmmsg,socket,bind,listen,unlink,unlinkat,rename,renameat,renameat2,truncate \
    -o "$trace_dir/$name.strace" \
    bwrap "${sandbox[@]}" env B2_READ_ONLY_ROOT_SANDBOX=1 B2_WRITABLE_MOUNT_POLICY="generated stage/scratch outputs, ephemeral tmpfs, merger-owned aggregate files" B2_FIRST_USE_LOG="$first_use" PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX="$pycache" LD_PRELOAD="$audit:$shim" PATH="$serial_bin:$PATH" WOLFRAMSCRIPT_REAL="$real_wolframscript" WOLFRAM_KERNEL_REAL="$real_wolfram_kernel" "$@" >"$log" 2>&1
}

environment_assert() {
  local consumer="$1" output_root="$2" stage0_digest="${STAGE0_CONTRACT_DIGEST:-}" environment_digest="${ENVIRONMENT_IDENTITY_DIGEST:-}"
  [[ "$stage0_digest" =~ ^[0-9a-f]{64}$ && "$environment_digest" =~ ^[0-9a-f]{64}$ ]] || usage
  mkdir -p "$output_root/environment_assertions"
  isolated python3 "$here/u1_body_b2_environment_gate.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "$stage0_digest" --environment-identity-digest "$environment_digest" --no-network-shim "$shim" --consumer "$consumer" --output "$output_root/environment_assertions/${consumer}.yaml"
}

trace_authenticated() {
  local trace_dir="$1" name="$2" log="$3" target="$4" evidence_root="$5"; shift 5
  environment_assert "$name" "$evidence_root"
  mkdir -p "$evidence_root/first_use_authentication"
  trace_run "$trace_dir" "$name" "$log" python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer "$name" --target "$target" --evidence "$evidence_root/first_use_authentication/${name}.yaml" -- "$@"
}

trace_authenticated_bound() {
  local trace_dir="$1" name="$2" log="$3" target="$4" evidence_root="$5" producer_record="$6"; shift 6
  environment_assert "$name" "$evidence_root"
  mkdir -p "$evidence_root/first_use_authentication"
  trace_run "$trace_dir" "$name" "$log" python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer "$name" --target "$target" --evidence "$evidence_root/first_use_authentication/${name}.yaml" --producer-record "$producer_record" -- "$@"
}

full_b1_replay() {
  local out="$scratch/full_b1" trace="$scratch/traces/full_b1" logs="$scratch/logs/full_b1"
  snapshot_initial_existence "$trace"
  mkdir -p "$out" "$logs"
  trace_run "$trace" phase_a_sympy "$logs/phase_a_sympy.log" python3 "$here/u1_body_dynamics_sympy.py" --output "$out/sympy_phase_a.json"
  trace_run "$trace" phase_a_wolfram "$logs/phase_a_wolfram.log" wolframscript -file "$here/u1_body_dynamics_dual.wl" --output "$out/mathematica_phase_a.json"
  trace_run "$trace" b1_sympy "$logs/b1_sympy.log" python3 "$here/u1_body_mechanics_sympy.py" --phase-artifact "$out/sympy_phase_a.json" --output "$out/sympy_phase_b1.yaml"
  trace_run "$trace" b1_wolfram "$logs/b1_wolfram.log" wolframscript -file "$here/u1_body_mechanics_dual.wl" --phase-artifact "$out/mathematica_phase_a.json" --output "$out/mathematica_phase_b1.yaml"
  trace_run "$trace" b1_compare_stage1 "$logs/b1_compare_stage1.log" python3 "$here/u1_body_mechanics_compare.py" --input "$here/u1_body_mechanics_inputs.yaml" --sympy-artifact "$out/sympy_phase_b1.yaml" --math-artifact "$out/mathematica_phase_b1.yaml" --phase-a-artifact "$out/sympy_phase_a.json" --verify-only --stage1-complete "$out/stage1_complete.yaml"
  local paths=() core_batches=12
  for batch in $(seq 0 $((core_batches-1))); do
    paths+=("$out/b1_mutation_core_${batch}.yaml")
    B2_B1_REPLAY_MUTATION_WRITES=1 B2_CAPTURE_WOLFRAM_ASSERT_FROM_ARTIFACT=1 trace_run "$trace" "mutation_core_${batch}" "$logs/mutation_core_${batch}.log" python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$out" --core-only --core-batch "$batch/$core_batches" --output "$out/b1_mutation_core_${batch}.yaml"
  done
  trace_run "$trace" mutation_core_aggregate "$logs/mutation_core_aggregate.log" python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$out" --core-only --core-batch-results "${paths[@]}" --output "$out/b1_mutation_core.yaml"
  for batch in 0 1 2; do
    trace_run "$trace" "mutation_liveness_${batch}" "$logs/mutation_liveness_${batch}.log" python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$out" --core-results "$out/b1_mutation_core.yaml" --liveness-batch "$batch/3" --output "$out/b1_liveness_${batch}.yaml"
  done
  trace_run "$trace" mutation_final "$logs/mutation_final.log" python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$out" --core-results "$out/b1_mutation_core.yaml" --liveness-results "$out/b1_liveness_0.yaml" "$out/b1_liveness_1.yaml" "$out/b1_liveness_2.yaml" --output "$out/b1_mutation_results.yaml"
  trace_run "$trace" seed_results "$logs/seed_results.log" /bin/cp "$here/reports/u1_body_dynamics_results.yaml" "$out/replay_results.yaml"
  trace_run "$trace" seed_report "$logs/seed_report.log" /bin/cp "$here/reports/u1_body_dynamics.md" "$out/replay_report.md"
  trace_run "$trace" b1_compare_stage2 "$logs/b1_compare_stage2.log" python3 "$here/u1_body_mechanics_compare.py" --artifacts "$out" --sympy-artifact "$out/sympy_phase_b1.yaml" --math-artifact "$out/mathematica_phase_b1.yaml" --phase-a-artifact "$out/sympy_phase_a.json" --mutation-results "$out/b1_mutation_results.yaml" --report "$out/replay_report.md" --results "$out/replay_results.yaml"
  local generated_args=(--generated-prefix "$out" --generated-prefix "$here/_scratch/u1_b1_mutations")
  for id in F1_wrong_angular_gram F1_wrong_normalization F1_perturbed_support_bound; do
    generated_args+=(--generated-prefix "$here/.mutation_${id}_u1_body_mechanics_sympy.py" --generated-prefix "$here/.mutation_${id}_u1_body_mechanics_dual.wl")
  done
  isolated python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$trace" --scope full_b1 --certificate "$certificate" --initial-existence "$trace/initial_existence.paths" "${generated_args[@]}" --output "$stage/b1_full_closure.yaml"
}

targeted_b1_replay() {
  local out="$scratch/targeted_b1" trace="$scratch/traces/targeted_b1" logs="$scratch/logs/targeted_b1"
  local approved_reemit=0
  if [[ "${STAGE0_CONTRACT_DIGEST:-}" =~ ^[0-9a-f]{64}$ && "${ENVIRONMENT_IDENTITY_DIGEST:-}" =~ ^[0-9a-f]{64}$ ]]; then approved_reemit=1; fi
  snapshot_initial_existence "$trace"
  mkdir -p "$out" "$logs"
  if (( approved_reemit )); then
    trace_authenticated "$trace" phase_a_sympy "$logs/phase_a_sympy.log" "$here/u1_body_dynamics_sympy.py" "$stage" python3 "$here/u1_body_dynamics_sympy.py" --output "$out/sympy_phase_a.json"
    trace_authenticated "$trace" phase_a_wolfram "$logs/phase_a_wolfram.log" "$here/u1_body_dynamics_dual.wl" "$stage" wolframscript -file "$here/u1_body_dynamics_dual.wl" --output "$out/mathematica_phase_a.json"
    trace_authenticated "$trace" targeted_phase_a_producer_seal "$logs/targeted_phase_a_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$stage" python3 "$here/u1_body_b2_producer_seal.py" --producer approved_reemit_targeted_phase_a --path "$out/sympy_phase_a.json" --path "$out/mathematica_phase_a.json" --output "$stage/b1_targeted_phase_a_producer_digests.yaml"
    trace_authenticated_bound "$trace" b1_sympy "$logs/b1_sympy.log" "$here/u1_body_mechanics_sympy.py" "$stage" "$stage/b1_targeted_phase_a_producer_digests.yaml" python3 "$here/u1_body_mechanics_sympy.py" --phase-artifact "$out/sympy_phase_a.json" --output "$out/sympy_phase_b1.yaml"
    trace_authenticated_bound "$trace" b1_wolfram "$logs/b1_wolfram.log" "$here/u1_body_mechanics_dual.wl" "$stage" "$stage/b1_targeted_phase_a_producer_digests.yaml" wolframscript -file "$here/u1_body_mechanics_dual.wl" --phase-artifact "$out/mathematica_phase_a.json" --output "$out/mathematica_phase_b1.yaml"
    trace_authenticated "$trace" targeted_producer_seal "$logs/targeted_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$stage" python3 "$here/u1_body_b2_producer_seal.py" --producer approved_reemit_targeted_engines --path "$out/sympy_phase_a.json" --path "$out/mathematica_phase_a.json" --path "$out/sympy_phase_b1.yaml" --path "$out/mathematica_phase_b1.yaml" --output "$stage/b1_targeted_producer_digests.yaml"
    trace_authenticated "$trace" targeted_compare "$logs/targeted_compare.log" "$here/u1_body_b2_targeted_replay.py" "$stage" python3 "$here/u1_body_b2_targeted_replay.py" --certificate "$certificate" --producer-record "$stage/b1_targeted_producer_digests.yaml" --protected-sympy "$here/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_b1.yaml" --protected-mathematica "$here/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_b1.yaml" --replay-sympy "$out/sympy_phase_b1.yaml" --replay-mathematica "$out/mathematica_phase_b1.yaml" --protected-phase-sympy "$here/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_a.json" --protected-phase-mathematica "$here/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_a.json" --replay-phase-sympy "$out/sympy_phase_a.json" --replay-phase-mathematica "$out/mathematica_phase_a.json" --output "$stage/b1_targeted_replay.yaml"
    environment_assert targeted_trace_audit "$stage"
  else
    trace_run "$trace" phase_a_sympy "$logs/phase_a_sympy.log" python3 "$here/u1_body_dynamics_sympy.py" --output "$out/sympy_phase_a.json"
    trace_run "$trace" phase_a_wolfram "$logs/phase_a_wolfram.log" wolframscript -file "$here/u1_body_dynamics_dual.wl" --output "$out/mathematica_phase_a.json"
    trace_run "$trace" b1_sympy "$logs/b1_sympy.log" python3 "$here/u1_body_mechanics_sympy.py" --phase-artifact "$out/sympy_phase_a.json" --output "$out/sympy_phase_b1.yaml"
    trace_run "$trace" b1_wolfram "$logs/b1_wolfram.log" wolframscript -file "$here/u1_body_mechanics_dual.wl" --phase-artifact "$out/mathematica_phase_a.json" --output "$out/mathematica_phase_b1.yaml"
    trace_run "$scratch/traces/stage0" targeted_producer_seal "$logs/targeted_producer_seal.log" python3 "$here/u1_body_b2_producer_seal.py" --producer stage0_targeted_replay_engines --path "$out/sympy_phase_a.json" --path "$out/mathematica_phase_a.json" --path "$out/sympy_phase_b1.yaml" --path "$out/mathematica_phase_b1.yaml" --output "$stage/b1_targeted_producer_digests.yaml"
    trace_run "$scratch/traces/stage0" targeted_compare "$logs/targeted_compare.log" python3 "$here/u1_body_b2_targeted_replay.py" --certificate "$certificate" --producer-record "$stage/b1_targeted_producer_digests.yaml" --protected-sympy "$here/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_b1.yaml" --protected-mathematica "$here/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_b1.yaml" --replay-sympy "$out/sympy_phase_b1.yaml" --replay-mathematica "$out/mathematica_phase_b1.yaml" --protected-phase-sympy "$here/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_a.json" --protected-phase-mathematica "$here/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_a.json" --replay-phase-sympy "$out/sympy_phase_a.json" --replay-phase-mathematica "$out/mathematica_phase_a.json" --output "$stage/b1_targeted_replay.yaml"
  fi
  if (( approved_reemit )); then
    isolated python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer targeted_trace_audit --target "$here/u1_body_b2_trace_audit.py" --evidence "$stage/first_use_authentication/targeted_trace_audit.yaml" -- python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$trace" --scope targeted_b1 --certificate "$certificate" --initial-existence "$trace/initial_existence.paths" --generated-prefix "$out" --generated-prefix "$stage" --output "$stage/b1_targeted_closure.yaml"
  else
    isolated python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$trace" --scope targeted_b1 --certificate "$certificate" --initial-existence "$trace/initial_existence.paths" --generated-prefix "$out" --generated-prefix "$stage" --output "$stage/b1_targeted_closure.yaml"
  fi
}

finalize_initial_stage0() {
  isolated python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$scratch/traces/stage0" --scope stage0 --certificate "$certificate" --initial-existence "$scratch/traces/stage0/initial_existence.paths" --generated-prefix "$stage" --generated-prefix "$scratch/stage0_mutations" --output "$stage/stage0_trace_closure.yaml"
  trace_run "$scratch/traces/stage0" contract_probe "$scratch/logs/stage0/contract_probe.log" python3 "$here/u1_body_b2_stage0_contract.py" --input "$input" --startup "$stage/startup_gate.yaml" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --agreement "$stage/stage0_engine_agreement.yaml" --schema-integration "$stage/stage0_schema_integration.yaml" --mutations "$stage/stage0_mutations.yaml" --b1-full-closure "$stage/b1_full_closure.yaml" --b1-targeted-replay "$stage/b1_targeted_replay.yaml" --b1-targeted-closure "$stage/b1_targeted_closure.yaml" --stage0-trace-closure "$stage/stage0_trace_closure.yaml" --output "$stage/stage0_contract.yaml"
  trace_run "$scratch/traces/stage0" summary_probe "$scratch/logs/stage0/summary_probe.log" python3 "$here/u1_body_b2_stage0_summary.py" --contract "$stage/stage0_contract.yaml" --output "$here/_scratch/_u1_phaseB2_build_stage0_summary.md"
  trace_run "$scratch/traces/stage0" trace_audit_probe "$scratch/logs/stage0/trace_audit_probe.log" python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$scratch/traces/stage0" --scope stage0 --certificate "$certificate" --initial-existence "$scratch/traces/stage0/initial_existence.paths" --generated-prefix "$stage" --generated-prefix "$scratch/stage0_mutations" --output "$scratch/stage0_trace_probe.yaml"
  isolated python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$scratch/traces/stage0" --scope stage0 --certificate "$certificate" --initial-existence "$scratch/traces/stage0/initial_existence.paths" --generated-prefix "$stage" --generated-prefix "$scratch/stage0_mutations" --output "$stage/stage0_trace_closure.yaml"
  isolated python3 "$here/u1_body_b2_stage0_contract.py" --input "$input" --startup "$stage/startup_gate.yaml" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --agreement "$stage/stage0_engine_agreement.yaml" --schema-integration "$stage/stage0_schema_integration.yaml" --mutations "$stage/stage0_mutations.yaml" --b1-full-closure "$stage/b1_full_closure.yaml" --b1-targeted-replay "$stage/b1_targeted_replay.yaml" --b1-targeted-closure "$stage/b1_targeted_closure.yaml" --stage0-trace-closure "$stage/stage0_trace_closure.yaml" --output "$stage/stage0_contract.yaml"
  isolated python3 "$here/u1_body_b2_stage0_summary.py" --contract "$stage/stage0_contract.yaml" --output "$here/_scratch/_u1_phaseB2_build_stage0_summary.md"
  test "$(wc -l < "$here/_scratch/_u1_phaseB2_build_stage0_summary.md")" -le 20
  echo "B2_STAGE0_HALT: artifact=$stage/stage0_contract.yaml sha256=$(sha256sum "$stage/stage0_contract.yaml" | awk '{print $1}')"
  exit 42
}

stage0() {
  rm -rf "$stage" "$scratch"
  rm -rf "$here/_scratch/u1_b1_mutations"
  mkdir -p "$stage" "$scratch/logs/stage0" "$scratch/traces/stage0"
  snapshot_initial_existence "$scratch/traces/stage0"
  trace_run "$scratch/traces/stage0" startup "$scratch/logs/stage0/startup.log" python3 "$here/u1_body_b2_startup_gate.py" --input "$input" --startup-contract-commit "$commit" --trust-mode "$trust_mode" --no-network-shim "$shim" --output "$stage/startup_gate.yaml"
  trace_run "$scratch/traces/stage0" sympy "$scratch/logs/stage0/sympy.log" python3 "$here/u1_body_b2_stage0_sympy.py" --input "$input" --output "$stage/sympy_stage0.yaml"
  trace_run "$scratch/traces/stage0" wolfram "$scratch/logs/stage0/wolfram.log" wolframscript -file "$here/u1_body_b2_stage0_dual.wl" --input "$input" --output "$stage/mathematica_stage0.yaml"
  trace_run "$scratch/traces/stage0" compare "$scratch/logs/stage0/compare.log" python3 "$here/u1_body_b2_stage0_compare.py" --input "$input" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --output "$stage/stage0_engine_agreement.yaml"
  trace_run "$scratch/traces/stage0" schema_integration "$scratch/logs/stage0/schema_integration.log" python3 "$here/u1_body_b2_schema_integration.py" --agreement "$stage/stage0_engine_agreement.yaml" --input "$input" --output "$stage/stage0_schema_integration.yaml"
  local mutation_batches=()
  for batch in 0 1 2 3; do
    mutation_batches+=("$stage/stage0_mutations_source_${batch}.yaml")
    trace_run "$scratch/traces/stage0" "mutations_source_${batch}" "$scratch/logs/stage0/mutations_source_${batch}.log" python3 "$here/u1_body_b2_stage0_mutations.py" --input "$input" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --work "$scratch/stage0_mutations" --source-batch "$batch/4" --output "$stage/stage0_mutations_source_${batch}.yaml"
  done
  for batch in 0 1 2 3 4 5; do
    mutation_batches+=("$stage/stage0_mutations_comparator_${batch}.yaml")
    trace_run "$scratch/traces/stage0" "mutations_comparator_${batch}" "$scratch/logs/stage0/mutations_comparator_${batch}.log" python3 "$here/u1_body_b2_stage0_mutations.py" --input "$input" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --work "$scratch/stage0_mutations" --comparator-batch "$batch/6" --output "$stage/stage0_mutations_comparator_${batch}.yaml"
  done
  mutation_batches+=("$stage/stage0_mutations_trace.yaml")
  trace_run "$scratch/traces/stage0" mutations_trace "$scratch/logs/stage0/mutations_trace.log" python3 "$here/u1_body_b2_stage0_mutations.py" --input "$input" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --work "$scratch/stage0_mutations" --trace-only --output "$stage/stage0_mutations_trace.yaml"
  trace_run "$scratch/traces/stage0" mutations_aggregate "$scratch/logs/stage0/mutations_aggregate.log" python3 "$here/u1_body_b2_stage0_mutations.py" --aggregate-results "${mutation_batches[@]}" --output "$stage/stage0_mutations.yaml"
  full_b1_replay
  targeted_b1_replay
  finalize_initial_stage0
}

finish_initial_stage0() {
  test ! -e "$stage/stage0_contract.yaml"
  test -s "$stage/startup_gate.yaml" -a -s "$stage/stage0_engine_agreement.yaml" -a -s "$stage/stage0_mutations.yaml"
  test -s "$scratch/full_b1/b1_mutation_results.yaml" -a -d "$scratch/traces/full_b1"
  local generated_args=(--generated-prefix "$scratch/full_b1" --generated-prefix "$here/_scratch/u1_b1_mutations")
  for id in F1_wrong_angular_gram F1_wrong_normalization F1_perturbed_support_bound; do
    generated_args+=(--generated-prefix "$here/.mutation_${id}_u1_body_mechanics_sympy.py" --generated-prefix "$here/.mutation_${id}_u1_body_mechanics_dual.wl")
  done
  isolated python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$scratch/traces/full_b1" --scope full_b1 --certificate "$certificate" --initial-existence "$scratch/traces/full_b1/initial_existence.paths" "${generated_args[@]}" --output "$stage/b1_full_closure.yaml"
  rm -rf "$scratch/targeted_b1" "$scratch/traces/targeted_b1" "$scratch/logs/targeted_b1"
  targeted_b1_replay
  finalize_initial_stage0
}

refresh_stage0_comparator_mutations() {
  local trace="$scratch/traces/stage0" logs="$scratch/logs/stage0"
  local mutation_batches=()
  for batch in 0 1 2 3; do
    test -s "$stage/stage0_mutations_source_${batch}.yaml"
    mutation_batches+=("$stage/stage0_mutations_source_${batch}.yaml")
  done
  for batch in 0 1 2 3 4 5; do
    rm -f "$trace"/mutations_comparator_${batch}.strace* "$trace/mutations_comparator_${batch}.firstuse" \
      "$logs/mutations_comparator_${batch}.log" "$stage/stage0_mutations_comparator_${batch}.yaml"
    mutation_batches+=("$stage/stage0_mutations_comparator_${batch}.yaml")
    trace_run "$trace" "mutations_comparator_${batch}" "$logs/mutations_comparator_${batch}.log" python3 "$here/u1_body_b2_stage0_mutations.py" --input "$input" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --work "$scratch/stage0_mutations" --comparator-batch "$batch/6" --output "$stage/stage0_mutations_comparator_${batch}.yaml"
  done
  test -s "$stage/stage0_mutations_trace.yaml"
  mutation_batches+=("$stage/stage0_mutations_trace.yaml")
  rm -f "$trace"/mutations_aggregate.strace* "$trace/mutations_aggregate.firstuse" \
    "$logs/mutations_aggregate.log" "$stage/stage0_mutations.yaml"
  trace_run "$trace" mutations_aggregate "$logs/mutations_aggregate.log" python3 "$here/u1_body_b2_stage0_mutations.py" --aggregate-results "${mutation_batches[@]}" --output "$stage/stage0_mutations.yaml"
}

clear_initial_liveness2_state() {
  local out="$scratch/full_b1" trace="$scratch/traces/full_b1" logs="$scratch/logs/full_b1"
  local case_dir case_id
  for case_dir in "$here"/_scratch/u1_b1_mutations/input_liveness/[0-9][0-9][0-9][0-9]; do
    [[ -d "$case_dir" ]] || continue
    case_id="$(basename "$case_dir")"
    if (( 10#$case_id % 3 == 2 )); then rm -rf "$case_dir"; fi
  done
  for name in mutation_liveness_2 mutation_final seed_results seed_report b1_compare_stage2; do
    rm -f "$trace"/${name}.strace* "$trace/${name}.firstuse" "$logs/${name}.log"
  done
  rm -f "$out/b1_liveness_2.yaml" "$out/b1_mutation_results.yaml" \
    "$out/replay_results.yaml" "$out/replay_report.md" "$stage/b1_full_closure.yaml"
}

complete_initial_full_b1_from_liveness2() {
  local out="$scratch/full_b1" trace="$scratch/traces/full_b1" logs="$scratch/logs/full_b1"
  trace_run "$trace" mutation_liveness_2 "$logs/mutation_liveness_2.log" python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$out" --core-results "$out/b1_mutation_core.yaml" --liveness-batch 2/3 --output "$out/b1_liveness_2.yaml"
  trace_run "$trace" mutation_final "$logs/mutation_final.log" python3 "$here/u1_body_mechanics_mutations.py" --artifacts "$out" --core-results "$out/b1_mutation_core.yaml" --liveness-results "$out/b1_liveness_0.yaml" "$out/b1_liveness_1.yaml" "$out/b1_liveness_2.yaml" --output "$out/b1_mutation_results.yaml"
  trace_run "$trace" seed_results "$logs/seed_results.log" /bin/cp "$here/reports/u1_body_dynamics_results.yaml" "$out/replay_results.yaml"
  trace_run "$trace" seed_report "$logs/seed_report.log" /bin/cp "$here/reports/u1_body_dynamics.md" "$out/replay_report.md"
  trace_run "$trace" b1_compare_stage2 "$logs/b1_compare_stage2.log" python3 "$here/u1_body_mechanics_compare.py" --artifacts "$out" --sympy-artifact "$out/sympy_phase_b1.yaml" --math-artifact "$out/mathematica_phase_b1.yaml" --phase-a-artifact "$out/sympy_phase_a.json" --mutation-results "$out/b1_mutation_results.yaml" --report "$out/replay_report.md" --results "$out/replay_results.yaml"
  finish_initial_stage0
}

resume_initial_full_b1() {
  local out="$scratch/full_b1" trace="$scratch/traces/full_b1" logs="$scratch/logs/full_b1"
  test ! -e "$stage/stage0_contract.yaml"
  test -s "$stage/startup_gate.yaml" -a -s "$stage/stage0_engine_agreement.yaml" -a -s "$stage/stage0_mutations.yaml"
  test -s "$trace/initial_existence.paths" -a -s "$out/b1_mutation_core.yaml"
  test -s "$out/b1_liveness_0.yaml" -a -s "$out/b1_liveness_1.yaml"
  test ! -s "$out/b1_liveness_2.yaml"
  test ! -s "$out/b1_mutation_results.yaml"

  # Refresh the comparator-only teeth after the first-time-construction fix.
  # The expensive source-ablation batches are unchanged and remain valid.
  refresh_stage0_comparator_mutations

  # The interrupted liveness batch can leave generated case outputs even when
  # its trace/output artifact was cut off.  Clear exactly the batch-2/3 case
  # directories as well as that step's partial evidence before replaying it.
  clear_initial_liveness2_state
  complete_initial_full_b1_from_liveness2
}

repair_initial_full_b1_liveness2() {
  test ! -e "$stage/stage0_contract.yaml"
  test -s "$stage/startup_gate.yaml" -a -s "$stage/stage0_engine_agreement.yaml" -a -s "$stage/stage0_mutations.yaml"
  test -s "$scratch/traces/full_b1/initial_existence.paths" -a -s "$scratch/full_b1/b1_mutation_core.yaml"
  test -s "$scratch/full_b1/b1_liveness_0.yaml" -a -s "$scratch/full_b1/b1_liveness_1.yaml"
  clear_initial_liveness2_state
  complete_initial_full_b1_from_liveness2
}

resume() {
  local stage0_digest="${STAGE0_CONTRACT_DIGEST:-}" environment_digest="${ENVIRONMENT_IDENTITY_DIGEST:-}" contract="$stage/stage0_contract.yaml"
  [[ "$stage0_digest" =~ ^[0-9a-f]{64}$ && "$environment_digest" =~ ^[0-9a-f]{64}$ ]] || usage
  local prod="$here/reports/u1_body_dynamics_artifacts/stage_b2_1_intake_radiative_production" trace="$scratch/traces/production" logs="$scratch/logs/production"
  rm -rf "$prod" "$trace" "$logs"; mkdir -p "$prod" "$logs"
  snapshot_initial_existence "$trace"
  environment_assert resume_entry "$prod"
  trace_authenticated "$trace" manifest_gate "$logs/manifest_gate.log" "$here/u1_body_b2_manifest_gate.py" "$prod" python3 "$here/u1_body_b2_manifest_gate.py" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --output "$prod/manifest_gate.yaml"
  trace_authenticated "$trace" startup_recheck "$logs/startup_recheck.log" "$here/u1_body_b2_startup_gate.py" "$prod" python3 "$here/u1_body_b2_startup_gate.py" --input "$input" --startup-contract-commit "$commit" --trust-mode "$trust_mode" --no-network-shim "$shim" --environment-identity-digest "$environment_digest" --output "$prod/startup_recheck.yaml"
  local tout="$prod/targeted_b1_outputs"
  mkdir -p "$tout"
  trace_authenticated "$trace" targeted_phase_a_sympy "$logs/targeted_phase_a_sympy.log" "$here/u1_body_dynamics_sympy.py" "$prod" python3 "$here/u1_body_dynamics_sympy.py" --output "$tout/sympy_phase_a.json"
  trace_authenticated "$trace" targeted_phase_a_wolfram "$logs/targeted_phase_a_wolfram.log" "$here/u1_body_dynamics_dual.wl" "$prod" wolframscript -file "$here/u1_body_dynamics_dual.wl" --output "$tout/mathematica_phase_a.json"
  trace_authenticated "$trace" targeted_phase_a_producer_seal "$logs/targeted_phase_a_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$prod" python3 "$here/u1_body_b2_producer_seal.py" --producer resume_targeted_phase_a_engines --path "$tout/sympy_phase_a.json" --path "$tout/mathematica_phase_a.json" --output "$prod/targeted_phase_a_producer_digests.yaml"
  trace_authenticated_bound "$trace" targeted_b1_sympy "$logs/targeted_b1_sympy.log" "$here/u1_body_mechanics_sympy.py" "$prod" "$prod/targeted_phase_a_producer_digests.yaml" python3 "$here/u1_body_mechanics_sympy.py" --phase-artifact "$tout/sympy_phase_a.json" --output "$tout/sympy_phase_b1.yaml"
  trace_authenticated_bound "$trace" targeted_b1_wolfram "$logs/targeted_b1_wolfram.log" "$here/u1_body_mechanics_dual.wl" "$prod" "$prod/targeted_phase_a_producer_digests.yaml" wolframscript -file "$here/u1_body_mechanics_dual.wl" --phase-artifact "$tout/mathematica_phase_a.json" --output "$tout/mathematica_phase_b1.yaml"
  trace_authenticated "$trace" targeted_producer_seal "$logs/targeted_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$prod" python3 "$here/u1_body_b2_producer_seal.py" --producer resume_targeted_replay_engines --path "$tout/sympy_phase_a.json" --path "$tout/mathematica_phase_a.json" --path "$tout/sympy_phase_b1.yaml" --path "$tout/mathematica_phase_b1.yaml" --output "$prod/targeted_producer_digests.yaml"
  trace_authenticated "$trace" targeted_replay "$logs/targeted_replay.log" "$here/u1_body_b2_targeted_replay.py" "$prod" python3 "$here/u1_body_b2_targeted_replay.py" --certificate "$certificate" --producer-record "$prod/targeted_producer_digests.yaml" --protected-sympy "$here/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_b1.yaml" --protected-mathematica "$here/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_b1.yaml" --replay-sympy "$tout/sympy_phase_b1.yaml" --replay-mathematica "$tout/mathematica_phase_b1.yaml" --protected-phase-sympy "$here/reports/u1_body_dynamics_artifacts/stage1/sympy_phase_a.json" --protected-phase-mathematica "$here/reports/u1_body_dynamics_artifacts/stage1/mathematica_phase_a.json" --replay-phase-sympy "$tout/sympy_phase_a.json" --replay-phase-mathematica "$tout/mathematica_phase_a.json" --output "$prod/targeted_replay.yaml"
  trace_authenticated "$trace" sympy "$logs/sympy.log" "$here/u1_body_b2_sympy.py" "$prod" python3 "$here/u1_body_b2_sympy.py" --input "$input" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --output "$prod/sympy_b2.yaml"
  trace_authenticated "$trace" wolfram "$logs/wolfram.log" "$here/u1_body_b2_dual.wl" "$prod" wolframscript -file "$here/u1_body_b2_dual.wl" --input "$input" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --output "$prod/mathematica_b2.yaml"
  trace_authenticated "$trace" production_producer_seal "$logs/production_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$prod" python3 "$here/u1_body_b2_producer_seal.py" --producer production_dual_engines --path "$prod/sympy_b2.yaml" --path "$prod/mathematica_b2.yaml" --output "$prod/production_producer_digests.yaml"
  trace_authenticated "$trace" compare "$logs/compare.log" "$here/u1_body_b2_compare.py" "$prod" python3 "$here/u1_body_b2_compare.py" --sympy "$prod/sympy_b2.yaml" --mathematica "$prod/mathematica_b2.yaml" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --producer-record "$prod/production_producer_digests.yaml" --output "$prod/engine_agreement.yaml"
  trace_authenticated "$trace" agreement_producer_seal "$logs/agreement_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$prod" python3 "$here/u1_body_b2_producer_seal.py" --producer production_comparator --path "$prod/engine_agreement.yaml" --output "$prod/agreement_producer_digests.yaml"
  trace_authenticated "$trace" mutations "$logs/mutations.log" "$here/u1_body_b2_mutations.py" "$prod" python3 "$here/u1_body_b2_mutations.py" --sympy "$prod/sympy_b2.yaml" --mathematica "$prod/mathematica_b2.yaml" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --work "$scratch/production_mutations" --output "$prod/mutations.yaml"
  trace_authenticated_bound "$trace" merger "$logs/merger.log" "$here/u1_body_b2_merge.py" "$prod" "$prod/agreement_producer_digests.yaml" python3 "$here/u1_body_b2_merge.py" --agreement "$prod/engine_agreement.yaml" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --b1-results-snapshot "$here/reports/u1_body_dynamics_artifacts/b1_final_results_snapshot.yaml" --b1-report-snapshot "$here/reports/u1_body_dynamics_artifacts/b1_final_report_snapshot.md" --results "$here/reports/u1_body_dynamics_results.yaml" --report "$here/reports/u1_body_dynamics.md" --provenance "$prod/merger_provenance.yaml"
  environment_assert production_trace_audit "$prod"
  isolated python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --consumer production_trace_audit --target "$here/u1_body_b2_trace_audit.py" --evidence "$prod/first_use_authentication/production_trace_audit.yaml" -- python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$trace" --scope stage0 --certificate "$certificate" --initial-existence "$trace/initial_existence.paths" --generated-prefix "$prod" --generated-prefix "$scratch/production_mutations" --output "$prod/production_trace_closure.yaml"
  trace_authenticated "$trace" completion_producer_seal "$logs/completion_producer_seal.log" "$here/u1_body_b2_producer_seal.py" "$prod" python3 "$here/u1_body_b2_producer_seal.py" --producer production_trace_and_merger --path "$prod/production_trace_closure.yaml" --path "$prod/merger_provenance.yaml" --output "$prod/completion_producer_digests.yaml"
  trace_authenticated_bound "$trace" completion_gate "$logs/completion_gate.log" "$here/u1_body_b2_completion_gate.py" "$prod" "$prod/completion_producer_digests.yaml" python3 "$here/u1_body_b2_completion_gate.py" --stage0 "$contract" --stage0-contract-digest "$stage0_digest" --environment-identity-digest "$environment_digest" --trace-closure "$prod/production_trace_closure.yaml" --merger-provenance "$prod/merger_provenance.yaml" --authentication-root "$prod" --output "$prod/completion_gate.yaml"
  echo "B2_RESUME: complete"
}

continue_stage0_after_full_replay() {
  [[ "${STAGE0_CONTRACT_DIGEST:-}" =~ ^[0-9a-f]{64}$ && "${ENVIRONMENT_IDENTITY_DIGEST:-}" =~ ^[0-9a-f]{64}$ ]] || usage
  environment_assert stage0_reemit "$stage"
  test -s "$stage/startup_gate.yaml"
  test -s "$stage/stage0_engine_agreement.yaml"
  test -s "$scratch/full_b1/b1_mutation_results.yaml"
  local generated_args=(--generated-prefix "$scratch/full_b1" --generated-prefix "$here/_scratch/u1_b1_mutations")
  for id in F1_wrong_angular_gram F1_wrong_normalization F1_perturbed_support_bound; do
    generated_args+=(--generated-prefix "$here/.mutation_${id}_u1_body_mechanics_sympy.py" --generated-prefix "$here/.mutation_${id}_u1_body_mechanics_dual.wl")
  done
  environment_assert stage0_reemit_full_closure "$stage"
  isolated python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$scratch/traces/full_b1" --scope full_b1 --certificate "$certificate" --initial-existence "$scratch/traces/full_b1/initial_existence.paths" "${generated_args[@]}" --output "$stage/b1_full_closure.yaml"
  rm -rf "$scratch/targeted_b1" "$scratch/traces/targeted_b1" "$scratch/logs/targeted_b1"
  targeted_b1_replay
  environment_assert stage0_reemit_trace_closure "$stage"
  isolated python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer stage0_reemit_trace_closure --target "$here/u1_body_b2_trace_audit.py" --evidence "$stage/first_use_authentication/stage0_reemit_trace_closure.yaml" -- python3 "$here/u1_body_b2_trace_audit.py" --trace-dir "$scratch/traces/stage0" --scope stage0 --certificate "$certificate" --initial-existence "$scratch/traces/stage0/initial_existence.paths" --generated-prefix "$stage" --generated-prefix "$scratch/stage0_mutations" --output "$stage/stage0_trace_closure.yaml"
  local candidate="$scratch/stage0_contract_reemit_candidate.yaml" candidate_seal="$scratch/stage0_contract_reemit_producer.yaml"
  environment_assert stage0_reemit_contract "$stage"
  isolated python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer stage0_reemit_contract --target "$here/u1_body_b2_stage0_contract.py" --evidence "$stage/first_use_authentication/stage0_reemit_contract.yaml" -- python3 "$here/u1_body_b2_stage0_contract.py" --input "$input" --startup "$stage/startup_gate.yaml" --sympy "$stage/sympy_stage0.yaml" --mathematica "$stage/mathematica_stage0.yaml" --agreement "$stage/stage0_engine_agreement.yaml" --schema-integration "$stage/stage0_schema_integration.yaml" --mutations "$stage/stage0_mutations.yaml" --b1-full-closure "$stage/b1_full_closure.yaml" --b1-targeted-replay "$stage/b1_targeted_replay.yaml" --b1-targeted-closure "$stage/b1_targeted_closure.yaml" --stage0-trace-closure "$stage/stage0_trace_closure.yaml" --container-path "$stage/stage0_contract.yaml" --output "$candidate"
  environment_assert stage0_reemit_candidate_seal "$stage"
  isolated python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer stage0_reemit_candidate_seal --target "$here/u1_body_b2_producer_seal.py" --evidence "$stage/first_use_authentication/stage0_reemit_candidate_seal.yaml" -- python3 "$here/u1_body_b2_producer_seal.py" --producer approved_stage0_reemit_contract --path "$candidate" --output "$candidate_seal"
  environment_assert stage0_reemit_summary "$stage"
  isolated python3 "$here/u1_body_b2_authenticated_exec.py" --stage0 "$stage/stage0_contract.yaml" --stage0-contract-digest "${STAGE0_CONTRACT_DIGEST}" --consumer stage0_reemit_summary --target "$here/u1_body_b2_stage0_summary.py" --evidence "$stage/first_use_authentication/stage0_reemit_summary.yaml" --producer-record "$candidate_seal" -- python3 "$here/u1_body_b2_stage0_summary.py" --contract "$candidate" --logical-contract "$stage/stage0_contract.yaml" --output "$here/_scratch/_u1_phaseB2_build_stage0_summary.md"
  test "$(wc -l < "$here/_scratch/_u1_phaseB2_build_stage0_summary.md")" -le 20
  environment_assert stage0_reemit_install "$stage"
  test "$(awk '/^[[:space:]]+sha256:/{print $2; exit}' "$candidate_seal")" = "$(sha256sum "$candidate" | awk '{print $1}')"
  mv "$candidate" "$stage/stage0_contract.yaml"
  echo "B2_STAGE0_HALT: artifact=$stage/stage0_contract.yaml sha256=$(sha256sum "$stage/stage0_contract.yaml" | awk '{print $1}')"
  exit 42
}

case "${1:-}" in --stage0) stage0 ;; --resume-initial-full-b1) resume_initial_full_b1 ;; --repair-initial-full-b1-liveness2) repair_initial_full_b1_liveness2 ;; --finish-initial-stage0) finish_initial_stage0 ;; --continue-stage0-after-full-replay) continue_stage0_after_full_replay ;; --resume) resume ;; *) usage ;; esac
