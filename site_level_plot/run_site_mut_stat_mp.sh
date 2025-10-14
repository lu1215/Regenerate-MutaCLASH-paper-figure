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

run PRG-1_site_stat python SiteMutStat.py --csv_path Output/PRG-1_site_level_2025_09_26_17_41_22/PRG-1_site_level_result.csv --experiment_name PRG-1 --output_dir Output/PRG-1_site_level_2025_09_26_17_41_22/ 
run ALG-1_site_stat python SiteMutStat.py --csv_path Output/ALG-1_site_level_2025_09_26_17_41_22/ALG-1_site_level_result.csv --experiment_name ALG-1 --output_dir Output/ALG-1_site_level_2025_09_26_17_41_22/
run Homo_site_stat python SiteMutStat.py --csv_path Output/Homo_site_level_2025_09_26_17_41_22/Homo_site_level_result.csv --experiment_name Homo --output_dir Output/Homo_site_level_2025_09_26_17_41_22/
run Fly_site_stat python SiteMutStat.py --csv_path Output/Fly_site_level_2025_09_27_04_37_36/Fly_site_level_result.csv --experiment_name Fly --output_dir Output/Fly_site_level_2025_09_27_04_37_36

# svg version

run PRG-1_site_stat python SiteMutStat.py --csv_path Output/PRG-1_site_level_2025_09_26_17_41_22/PRG-1_site_level_result.csv --experiment_name PRG-1 --output_dir Output/PRG-1_site_level_2025_09_26_17_41_22/ --image_format svg
run ALG-1_site_stat python SiteMutStat.py --csv_path Output/ALG-1_site_level_2025_09_26_17_41_22/ALG-1_site_level_result.csv --experiment_name ALG-1 --output_dir Output/ALG-1_site_level_2025_09_26_17_41_22/ --image_format svg
run Homo_site_stat python SiteMutStat.py --csv_path Output/Homo_site_level_2025_09_26_17_41_22/Homo_site_level_result.csv --experiment_name Homo --output_dir Output/Homo_site_level_2025_09_26_17_41_22/ --image_format svg
run Fly_site_stat python SiteMutStat.py --csv_path Output/Fly_site_level_2025_09_27_04_37_36/Fly_site_level_result.csv --experiment_name Fly --output_dir Output/Fly_site_level_2025_09_27_04_37_36 --image_format svg

wait_group