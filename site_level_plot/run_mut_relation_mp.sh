#!/usr/bin/env bash
# parallel_plot_mutation_relation.sh
set -euo pipefail

# 同時跑幾個任務，可用環境變數覆蓋：MAX_JOBS=6 ./parallel_plot_mutation_relation.sh
MAX_JOBS="${MAX_JOBS:-4}"

wait_group() {
    local ok=0 fail=0
    for pid in $(jobs -rp); do
    if wait "$pid"; then ok=$((ok+1)); else fail=$((fail+1)); fi
    done
    echo "[$(date '+%F %T')] GROUP DONE: ok=${ok}, fail=${fail}"
    # 想讓失敗直接終止整體流程的話，取消下一行註解
    # (( fail == 0 )) || exit 1
}

# 你的三組 (species, csv, outdir)
# groups=(
#   "ALG-1,/media/disk1/shangyi/MutaCLASH_additional_analysis/Output/ALG-1_pair_level_result/ALG-1_pair_level_result.csv,Output/ALG-1_mixed_relation"
#   "PRG-1,/media/disk1/shangyi/MutaCLASH_additional_analysis/Output/PRG-1_pair_level_result/PRG-1_pair_level_result.csv,Output/PRG-1_mixed_relation"
#   "PRG-1,/media/disk1/shangyi/MutaCLASH_additional_analysis/Output/PRG-1_pair_level_result_WAGO_target/PRG-1_WAGO_target_pair_level_result.csv,Output/PRG-1_wago_target_mixed_relation"
# )

ts_md() { date +%Y_%m_%d_%H_%M_%S; }

ts="$(ts_md)"

# site level group
mkdir -p Output/ALG-1_site_level_mixed_relation_$ts \
            Output/PRG-1_site_level_mixed_relation_$ts \
            Output/Homo_site_level_mixed_relation_$ts \
            Output/Fly_site_level_mixed_relation_$ts

groups=(
    "ALG-1,Output/ALG-1_site_level_2025_09_26_17_41_22/ALG-1_site_level_result.csv,Output/ALG-1_site_level_mixed_relation_$ts"
    "PRG-1,Output/PRG-1_site_level_2025_09_26_17_41_22/PRG-1_site_level_result.csv,Output/PRG-1_site_level_mixed_relation_$ts"
    "Homo,Output/Homo_site_level_2025_09_26_17_41_22/Homo_site_level_result.csv,Output/Homo_site_level_mixed_relation_$ts"
    "Fly,Output/Fly_site_level_2025_09_27_04_37_36/Fly_site_level_result.csv,Output/Fly_site_level_mixed_relation_$ts"
)

# # hybrid level groups
# groups=(
#     "ALG-1,/media/disk1/shangyi/MutaCLASH_additional_analysis/Input/ALG-1_rep3.csv,Output/ALG-1_mixed_relation_hybrid"
#     "PRG-1,/media/disk1/shangyi/MutaCLASH_additional_analysis/Input/PRG-1_piRNA_neworder_cutadapt.csv,Output/PRG-1_new_setting_mixed_relation_hybrid"
#     "Homo,HOMO_CLIP1-4_2025-09-18_19-57-03_MutaCLASH_right_preprocess/HOMO_CLIP1-4.csv,Output/Homo_mixed_relation_hybrid"
#     "Fly,/media/disk1/shangyi/test_MutaCLASH/MutaCLASH0908/data/output/SRR29206632_2025-09-19_03-48-13_MutaCLASH/SRR29206632.csv,Output/Fly_mixed_relation_hybrid"
#     # "Fly,/media/disk1/shangyi/test_MutaCLASH/MutaCLASH0908/data/output/SRR29206632_2025-09-20_00-37-02_MutaCLASH_with_a/SRR29206632.csv,Output/Fly_mixed_relation_hybrid_with_a"
# )

modes=(score readcount CIMS)
# muts=(mis del all)
muts=(mis del)

mkdir -p logs
start_time=$(date +%s)

run() {
  local name="$1"; shift
  "$@" >"logs/${name}.log" 2>&1 &
  # 節流：若達到並行上限，等任一個結束
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n || true
  done
}



for g in "${groups[@]}"; do
  IFS=',' read -r species csv outdir <<<"$g"
  mkdir -p "$outdir"

  for mode in "${modes[@]}"; do
    for mut in "${muts[@]}"; do
      # 組一個好辨識的任務名稱當 log 檔名
      name="$(basename "$outdir")_${species}_${mode}_${mut}"
      run "$name" \
        python plot_mixed_relation.py \
          --csv "$csv" \
          --species "$species" \
          --mut "$mut" \
          --output "$outdir" \
          --mode "$mode" \
          --pair
    done
  done
done

# 等最後一批
wait_group

elapsed=$(( $(date +%s) - start_time ))
echo "ALL DONE in ${elapsed}s. Logs under ./logs/"
