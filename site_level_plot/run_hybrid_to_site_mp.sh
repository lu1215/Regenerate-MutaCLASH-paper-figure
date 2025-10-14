#!/usr/bin/env bash
set -euo pipefail
start_time=$(date +%s)
mkdir -p logs

MAX_JOBS="${MAX_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

run() {
    local name="$1"; shift
    echo "[$(date '+%F %T')] START $name"
    "$@" >"logs/${name}.log" 2>&1 &
    # 達上限就等任一任務完成
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n || true
    done
}

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

mkdir -p Output/PRG-1_site_level_$ts/ \
         Output/ALG-1_site_level_$ts \
         Output/Homo_site_level_$ts \
         Output/Fly_site_level_$ts \
        
mkdir -p Output/PRG-1_wago_site_level_$ts/


run gen_PRG-1_site python HybridReadTransToPairLevel_pos_aware.py --csv_file Input/PRG-1_piRNA_neworder_cutadapt.csv --sample_name PRG-1 --output Output/PRG-1_site_level_$ts
run gen_ALG-1_site python HybridReadTransToPairLevel_pos_aware.py --csv_file Input/ALG-1_rep3.csv --sample_name ALG-1 --output Output/ALG-1_site_level_$ts
run gen_Fly_site python HybridReadTransToPairLevel_pos_aware.py --csv_file Input/SRR29206632_with_cage_abu.csv --sample_name Fly --output Output/Fly_site_level_$ts
run gen_Homo_site python HybridReadTransToPairLevel_pos_aware.py --csv_file Input/HOMO_CLIP1-4.csv --sample_name Homo --output Output/Homo_site_level_$ts
run gen_PRG-1_wago_site python HybridReadTransToPairLevel_pos_aware.py --csv_file Input/PRG-1_piRNA_neworder_cutadapt_Wago_target.csv --sample_name PRG-1 --output Output/PRG-1_wago_site_level_$ts/

wait_group