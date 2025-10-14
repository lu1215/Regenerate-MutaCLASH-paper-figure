#!/usr/bin/env bash
# parallel_abu_sgl_stat.sh
set -euo pipefail

start_time=$(date +%s)
mkdir -p logs

# 同時跑幾個任務；可用環境變數覆蓋：MAX_JOBS=6 ./parallel_abu_sgl_stat.sh
MAX_JOBS="${MAX_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

# 小工具：啟動背景任務 + 節流
run() {
  local name="$1"; shift
  echo "[$(date '+%F %T')] START $name"
  "$@" >"logs/${name}.log" 2>&1 &
  # 達上限就等任一任務完成
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n || true
  done
}

# 等待目前這批全部完成
wait_group() {
  local ok=0 fail=0
  for pid in $(jobs -rp); do
    if wait "$pid"; then ok=$((ok+1)); else fail=$((fail+1)); fi
  done
  echo "[$(date '+%F %T')] GROUP DONE: ok=${ok}, fail=${fail}"
  # 想讓失敗直接終止整體流程的話，取消下一行註解
  # (( fail == 0 )) || exit 1
}

ts_md() { date +%Y_%m_%d_%H_%M_%S; }

ts="$(ts_md)"

################
## site level ##
################

mkdir -p Output/PRG-1_abu_score_sgl_stat_site_$ts/ \
         Output/PRG-1_abu_readcount_sgl_stat_site_$ts \
         Output/PRG-1_abu_sgl_stat_combined_site_$ts \
         Output/ALG-1_abu_score_sgl_stat_site_$ts \
         Output/ALG-1_abu_readcount_sgl_stat_site_$ts \
         Output/ALG-1_abu_sgl_stat_combined_site_$ts 

run PRG1_sgl_score_site \
  python 22GAbuGroupbysglstatsig.py \
    --csv /media/disk1/shangyi/single_read_mut_stat/output_2025_09_28_15_41_03/PRG-1_site_level_result_with_sgl_mut_stat.csv \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode score \
    --out_dir Output/PRG-1_abu_score_sgl_stat_site_$ts/ \
    --pair true


run PRG1_sgl_readcount_site \
  python 22GAbuGroupbysglstatsig.py \
    --csv /media/disk1/shangyi/single_read_mut_stat/output_2025_09_28_15_41_03/PRG-1_site_level_result_with_sgl_mut_stat.csv \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode readcount \
    --out_dir Output/PRG-1_abu_readcount_sgl_stat_site_$ts/ \
    --boundaries 10000/4/1/0 \
    --pair true

run PRG1_sgl_combined_site \
  python 22GAbuGroupbysglstatsig.py \
    --csv /media/disk1/shangyi/single_read_mut_stat/output_2025_09_28_15_41_03/PRG-1_site_level_result_with_sgl_mut_stat.csv  \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode combined \
    --out_dir Output/PRG-1_abu_sgl_stat_combined_site_$ts \
    --pair true

#（如要啟用 ALG-1 版，照樣加 run … 指令即可；你原本那兩行我保留為註解）
run ALG1_sgl_score_site python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/output_2025_09_28_15_41_03/ALG-1_site_level_result_with_sgl_mut_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/ALG-1_abu_score_sgl_stat_site_$ts --boundaries 200/140/100/60 --RNA miRNA  --pair true
run ALG1_sgl_readcount_site python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/output_2025_09_28_15_41_03/ALG-1_site_level_result_with_sgl_mut_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_abu_readcount_sgl_stat_site_$ts --boundaries 10000/4/1/0 --RNA miRNA --pair true 
run ALG1_sgl_combined_site python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/output_2025_09_28_15_41_03/ALG-1_site_level_result_with_sgl_mut_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/ALG-1_abu_sgl_stat_combined_site_$ts --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0 --pair true

wait_group   # Stage A 全部完成後再進 Stage B

#======== Stage B：吃 combined 輸出的統計（可並行；此處只有一個）========#
mkdir -p Output/PRG-1_abu_sgl_stat_combined_site_$ts
run PRG1_sgl_proc_site \
  python process_stat_csv.py \
    -i Output/PRG-1_abu_sgl_stat_combined_site_$ts/all_stats_summary.csv \
    -p PRG-1_wago -r piRNA \
    -o Output/PRG-1_abu_sgl_stat_combined_site_$ts \
    -t WAGO-1 --suff_rev

run ALG1_sgl_proc_site python process_stat_csv.py -i Output/ALG-1_abu_sgl_stat_combined_site_$ts/all_stats_summary.csv -p ALG-1 -r miRNA -o Output/ALG-1_abu_sgl_stat_combined_site_$ts --suff_rev

wait_group   # Stage B 全部完成


##################################
## site level without site CIMS ##
##################################

mkdir -p Output/PRG-1_abu_score_sgl_stat_site_wo_mut_$ts/ \
         Output/PRG-1_abu_readcount_sgl_stat_site_wo_mut_$ts \
         Output/PRG-1_abu_sgl_stat_combined_site_wo_mut_$ts \
         Output/ALG-1_abu_score_sgl_stat_site_wo_mut_$ts \
         Output/ALG-1_abu_readcount_sgl_stat_site_wo_mut_$ts \
         Output/ALG-1_abu_sgl_stat_combined_site_wo_mut_$ts 

run PRG1_sgl_score_site_wo_mut \
  python 22GAbuGroupbysglstatsig.py \
    --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/PRG-1_site_level_result_with_sgl_mut_stat_site_wo_mut_.csv \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode score \
    --out_dir Output/PRG-1_abu_score_sgl_stat_site_wo_mut_$ts/ \
    --pair true


run PRG1_sgl_readcount_site_wo_mut \
  python 22GAbuGroupbysglstatsig.py \
    --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/PRG-1_site_level_result_with_sgl_mut_stat_site_wo_mut_.csv \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode readcount \
    --out_dir Output/PRG-1_abu_readcount_sgl_stat_site_wo_mut_$ts/ \
    --boundaries 10000/4/1/0 \
    --pair true

run PRG1_sgl_combined_site_wo_mut \
  python 22GAbuGroupbysglstatsig.py \
    --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/PRG-1_site_level_result_with_sgl_mut_stat_site_wo_mut_.csv  \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode combined \
    --out_dir Output/PRG-1_abu_sgl_stat_combined_site_wo_mut_$ts \
    --pair true

#（如要啟用 ALG-1 版，照樣加 run … 指令即可；你原本那兩行我保留為註解）
run ALG1_sgl_score_site_wo_mut python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/ALG-1_site_level_result_with_sgl_mut_stat_site_wo_mut_.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/ALG-1_abu_score_sgl_stat_site_wo_mut_$ts --boundaries 200/140/100/60 --RNA miRNA  --pair true
run ALG1_sgl_readcount_site_wo_mut python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/ALG-1_site_level_result_with_sgl_mut_stat_site_wo_mut_.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_abu_readcount_sgl_stat_site_wo_mut_$ts --boundaries 10000/4/1/0 --RNA miRNA --pair true 
run ALG1_sgl_combined_site_wo_mut python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/ALG-1_site_level_result_with_sgl_mut_stat_site_wo_mut_.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/ALG-1_abu_sgl_stat_combined_site_wo_mut_$ts --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0 --pair true

wait_group   # Stage A 全部完成後再進 Stage B

#======== Stage B：吃 combined 輸出的統計（可並行；此處只有一個）========#
mkdir -p Output/PRG-1_abu_sgl_stat_combined_site_wo_mut_$ts
run PRG1_sgl_proc_site_wo_mut \
  python process_stat_csv.py \
    -i Output/PRG-1_abu_sgl_stat_combined_site_wo_mut_$ts/all_stats_summary.csv \
    -p PRG-1_wago -r piRNA \
    -o Output/PRG-1_abu_sgl_stat_combined_site_wo_mut_$ts \
    -t WAGO-1 --suff_rev

run ALG1_sgl_proc_site_wo_mut python process_stat_csv.py -i Output/ALG-1_abu_sgl_stat_combined_site_wo_mut_$ts/all_stats_summary.csv -p ALG-1 -r miRNA -o Output/ALG-1_abu_sgl_stat_combined_site_wo_mut_$ts --suff_rev

wait_group   # Stage B 全部完成
# ############################
# # ## without hybrid mutation
# ############################
# #======== Stage A：產出各 out_dir（可並行）========#
# mkdir -p Output/PRG-1_abu_score_sgl_stat_wo_mut_$ts/ \
#          Output/PRG-1_abu_readcount_sgl_stat_wo_mut_$ts \
#          Output/PRG-1_abu_sgl_stat_combined_wo_mut_$ts \
#          Output/ALG-1_abu_score_sgl_stat_wo_mut_$ts \
#          Output/ALG-1_abu_readcount_sgl_stat_wo_mut_$ts \
#          Output/ALG-1_abu_sgl_stat_combined_wo_mut_$ts 

# run PRG1_sgl_score_wo_mut \
#   python 22GAbuGroupbysglstatsig.py \
#     --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/PRG-1_piRNA_neworder_cutadapt_hybrid_wo_mut_.csv \
#     --target Input/add_two_HCLee.RNAseq.master.xlsx \
#     --id_map Input/mRNA_WS275_IDtoName.csv \
#     --mode score \
#     --out_dir Output/PRG-1_abu_score_sgl_stat_wo_mut_$ts/

# run PRG1_sgl_readcount_wo_mut \
#   python 22GAbuGroupbysglstatsig.py \
#     --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/PRG-1_piRNA_neworder_cutadapt_hybrid_wo_mut_.csv \
#     --target Input/add_two_HCLee.RNAseq.master.xlsx \
#     --id_map Input/mRNA_WS275_IDtoName.csv \
#     --mode readcount \
#     --out_dir Output/PRG-1_abu_readcount_sgl_stat_wo_mut_$ts/ \
#     --boundaries 10000/4/1/0

# run PRG1_sgl_combined_wo_mut \
#   python 22GAbuGroupbysglstatsig.py \
#     --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/PRG-1_piRNA_neworder_cutadapt_hybrid_wo_mut_.csv  \
#     --target Input/add_two_HCLee.RNAseq.master.xlsx \
#     --id_map Input/mRNA_WS275_IDtoName.csv \
#     --mode combined \
#     --out_dir Output/PRG-1_abu_sgl_stat_combined_wo_mut_$ts

# #（如要啟用 ALG-1 版，照樣加 run … 指令即可；你原本那兩行我保留為註解）
# run ALG1_sgl_score_wo_mut python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/ALG-1_rep3_with_sgl_stat_hybrid_wo_mut_.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/ALG-1_abu_score_sgl_stat_wo_mut_$ts --boundaries 200/140/100/60 --RNA miRNA 
# run ALG1_sgl_readcount_wo_mut python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/ALG-1_rep3_with_sgl_stat_hybrid_wo_mut_.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_abu_readcount_sgl_stat_wo_mut_$ts --boundaries 10000/4/1/0 --RNA miRNA
# run ALG1_sgl_combined_wo_mut python 22GAbuGroupbysglstatsig.py --csv /media/disk1/shangyi/single_read_mut_stat/wo_mut_res/ALG-1_rep3_with_sgl_stat_hybrid_wo_mut_.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/ALG-1_abu_sgl_stat_combined_wo_mut_$ts --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0

# wait_group   # Stage A 全部完成後再進 Stage B

# #======== Stage B：吃 combined 輸出的統計（可並行；此處只有一個）========#
# mkdir -p Output/PRG-1_abu_sgl_stat_combined_wo_mut_$ts
# run PRG1_sgl_proc_wo_mut \
#   python process_stat_csv.py \
#     -i Output/PRG-1_abu_sgl_stat_combined_wo_mut_$ts/all_stats_summary.csv \
#     -p PRG-1_wago -r piRNA \
#     -o Output/PRG-1_abu_sgl_stat_combined_wo_mut_$ts \
#     -t WAGO-1 --suff_rev

# run ALG1_sgl_proc_wo_mut python process_stat_csv.py -i Output/ALG-1_abu_sgl_stat_combined_wo_mut_$ts/all_stats_summary.csv -p ALG-1 -r miRNA -o Output/ALG-1_abu_sgl_stat_combined_wo_mut_$ts --suff_rev

# wait_group   # Stage B 全部完成

elapsed=$(( $(date +%s) - start_time ))
echo "ALL DONE in ${elapsed}s. See logs/ for details."
