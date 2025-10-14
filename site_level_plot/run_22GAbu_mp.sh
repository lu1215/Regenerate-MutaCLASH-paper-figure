#!/usr/bin/env bash
# Bash 4.3.48(1)

set -euo pipefail

start_time=$(date +%s)
mkdir -p logs

# 並行度（同時跑幾個），可依機器核心數調整或用環境變數覆蓋：MAX_JOBS=8 ./run.sh
MAX_JOBS="${MAX_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

#======== 共用小工具 ========#
# 以背景行程啟動一個任務，輸出寫到 logs/<name>.log
run() {
  local name="$1"; shift
  echo "[`date '+%F %T'`] START $name"
  "$@" >"logs/${name}.log" 2>&1 &
  # 節流：若達到並行上限，就等任一任務結束再繼續塞
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n || true
  done
}

# 等待目前所有背景任務並報告成功/失敗
wait_group() {
  local ok=0 fail=0
  # jobs -rp：只列出仍在跑的背景行程 PID
  for pid in $(jobs -rp); do
    if wait "$pid"; then
      ok=$((ok+1))
    else
      fail=$((fail+1))
    fi
  done
  echo "[`date '+%F %T'`] GROUP DONE: ok=${ok}, fail=${fail}"
  # 若需要讓整體腳本在子任務失敗時也失敗，取消下一行註解：
  # (( fail == 0 )) || exit 1
}

ts_md() { date +%Y_%m_%d_%H_%M_%S; }

ts="$(ts_md)"

## processing site data
run PRG1_site_rc      python 22GAbuAnalysis.py --csv Output/PRG-1_site_level_2025_09_26_17_41_22/PRG-1_site_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/PRG-1_site_level_22G_readcount_$ts --pairing true --boundaries 10000/4/1/0
run PRG1_site_score   python 22GAbuAnalysis.py --csv Output/PRG-1_site_level_2025_09_26_17_41_22/PRG-1_site_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score     --out_dir Output/PRG-1_site_level_22G_score_$ts     --pairing true
run PRG1_site_comb    python 22GAbuAnalysis.py --csv Output/PRG-1_site_level_2025_09_26_17_41_22/PRG-1_site_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined  --out_dir Output/PRG-1_site_level_combined_$ts      --pairing true

run ALG1_site_rc      python 22GAbuAnalysis.py --csv Output/ALG-1_site_level_2025_09_26_17_41_22/ALG-1_site_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_site_level_abu_readcount_$ts --pairing true --RNA miRNA --boundaries 10000/4/1/0
run ALG1_site_score   python 22GAbuAnalysis.py --csv Output/ALG-1_site_level_2025_09_26_17_41_22/ALG-1_site_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score     --out_dir Output/ALG-1_site_level_abu_score_$ts     --pairing true --RNA miRNA --boundaries 200/140/100/60
run ALG1_site_comb    python 22GAbuAnalysis.py --csv Output/ALG-1_site_level_2025_09_26_17_41_22/ALG-1_site_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined  --out_dir Output/ALG-1_site_level_combined_$ts      --pairing true --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0


# Ctrl-C 時把孩子們一起關掉，避免殭屍行程
trap 'echo "Interrupted, stopping children..."; jobs -pr | xargs -r kill; wait' INT TERM

#======== Stage A：22GAbuAnalysis 同步啟動（彼此獨立可並跑）========#
# run PRG1_score        python 22GAbuAnalysis.py --csv Input/PRG-1_piRNA_neworder_cutadapt.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score     --out_dir Output/PRG-1_22G_test_score_cutadapt_neworder
# run PRG1_readcount    python 22GAbuAnalysis.py --csv Input/PRG-1_piRNA_neworder_cutadapt.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/PRG-1_22G_test_readcount_cutadapt_neworder --boundaries 10000/4/1/0
# run PRG1_pair_rc      python 22GAbuAnalysis.py --csv Output/PRG-1_pair_level_result_cutadapt_neworder/PRG-1_piRNA_neworder_cutadapt_pair_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/PRG-1_pairing_level_22G_readcount_cutadapt_neworder --pairing true --boundaries 10000/4/1/0
# run PRG1_pair_score   python 22GAbuAnalysis.py --csv Output/PRG-1_pair_level_result_cutadapt_neworder/PRG-1_piRNA_neworder_cutadapt_pair_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score     --out_dir Output/PRG-1_pairing_level_22G_score_cutadapt_neworder    --pairing true
# run PRG1_pair_comb    python 22GAbuAnalysis.py --csv Output/PRG-1_pair_level_result_cutadapt_neworder/PRG-1_piRNA_neworder_cutadapt_pair_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined  --out_dir Output/PRG-1_pair_combined_cutadapt_neworder            --pairing true
# run PRG1_hybrid_comb  python 22GAbuAnalysis.py --csv Input/PRG-1_piRNA_neworder_cutadapt.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined  --out_dir Output/PRG-1_hybrid_combined_cutadapt_neworder

# run ALG1_score python 22GAbuAnalysis.py --csv Input/ALG-1_rep3.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/ALG-1_abu_test_score_violin_and_datapoint --boundaries 200/140/100/60 --RNA miRNA
# run ALG1_readcount python 22GAbuAnalysis.py --csv Input/ALG-1_rep3.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_abu_test_readcount_violin_and_datapoint --boundaries 10000/4/1/0 --RNA miRNA
# run ALG1_pair_rc python 22GAbuAnalysis.py --csv Output/ALG-1_pair_level_result_new_merged_method/ALG-1_pair_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/ALG-1_pairing_level_abu_score_merged_pos --pairing true --RNA miRNA --boundaries 200/140/100/60
# run ALG1_pair_score python 22GAbuAnalysis.py --csv Output/ALG-1_pair_level_result_new_merged_method/ALG-1_pair_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_pairing_level_abu_readcount_merged_pos --pairing true --RNA miRNA --boundaries 10000/4/1/0
# run ALG1_pair_comb python 22GAbuAnalysis.py --csv Output/ALG-1_pair_level_result_new_merged_method/ALG-1_pair_level_result.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/ALG-1_pair_combined_merged_pos --pairing true --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0
# run ALG1_hybrid_comb python 22GAbuAnalysis.py --csv Input/ALG-1_rep3.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/ALG-1_hybrid_combined_violin_and_datapoint --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0


# 若還要把 ALG-1 那組也並跑，照樣加 run <name> <command>（你原本是註解掉的）

wait_group   # Stage A 全部完成後再進入下一階段

#======== Stage B：依賴 Stage A 產物的統計任務，並跑 ========#
# run ALG1_pair_stats    python process_stat_csv.py -i Output/ALG-1_pair_combined_merged_pos/all_stats_summary.csv   -p ALG-1      -r miRNA -o Output/ALG-1_pair_combined_merged_pos
# run PRG1_pair_stats    python process_stat_csv.py -i Output/PRG-1_pair_combined_cutadapt_neworder/all_stats_summary.csv   -p PRG-1_wago -r piRNA -o Output/PRG-1_pair_combined_cutadapt_neworder/   -t WAGO-1
# run ALG1_hybrid_stats  python process_stat_csv.py -i Output/ALG-1_hybrid_combined_violin_and_datapoint/all_stats_summary.csv -p ALG-1      -r miRNA -o Output/ALG-1_hybrid_combined_violin_and_datapoint/
# run PRG1_hybrid_stats  python process_stat_csv.py -i Output/PRG-1_hybrid_combined_cutadapt_neworder/all_stats_summary.csv -p PRG-1_wago -r piRNA -o Output/PRG-1_hybrid_combined_cutadapt_neworder/ -t WAGO-1
run PRG1_site_stats python process_stat_csv.py -i Output/PRG-1_site_level_combined_$ts/all_stats_summary.csv   -p PRG-1_wago -r piRNA -o Output/PRG-1_site_level_combined_$ts/   -t WAGO-1
# run PRG1_pair_stats python process_stat_csv.py -i Output/PRG-1_site_level_combined_2025_09_27_13_01_16/all_stats_summary.csv   -p PRG-1_wago -r piRNA -o Output/PRG-1_site_level_combined_$ts/   -t WAGO-1
run ALG1_site_stats python process_stat_csv.py -i Output/ALG-1_site_level_combined_$ts/all_stats_summary.csv   -p ALG-1 -r miRNA -o Output/ALG-1_site_level_combined_$ts/

wait_group   # Stage B 全部完成

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "Execution time: ${elapsed} seconds" | tee -a time.log
