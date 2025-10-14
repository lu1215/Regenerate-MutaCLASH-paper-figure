# -*- coding: utf-8 -*-
"""
22GAbuAnalysis_sig.py
- 在原始 22GAbuAnalysis 基礎上加入「顯著性分群」：
  依 single_read_mut_stat_significance_D / _M 切成 D=True、M=True、兩者皆 False 三組，
  與原本的 with/without mutation 流程互斥（啟用 --sig_grouping 時即走顯著性分群）。
- 相容 Python 3.6。

依賴：
  numpy, pandas, scipy, matplotlib, seaborn, (xlsx 需 openpyxl)

使用範例（啟用顯著性分群）：
  python 22GAbuAnalysis_sig.py --csv Input/PRG-1_rep1.csv \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode score \
    --out_dir Output/PRG-1_22G_test_score/ \
    --sig_grouping

  python 22GAbuAnalysis_sig.py --csv Input/ALG-1_rep3.csv \
    --target Input/add_two_HCLee.RNAseq.master.xlsx \
    --id_map Input/mRNA_WS275_IDtoName.csv \
    --mode readcount \
    --out_dir Output/ALG-1_abu_test_readcount/ \
    --boundaries 10000/4/1/0 \
    --RNA miRNA \
    --sig_grouping --sig_overlap M_priority
"""

import math
import gc
import os
from typing import Dict, Tuple, List
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")  # 伺服器/批次環境安全輸出
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------------------------------------------------------
#  Global constants & mappings
# -----------------------------------------------------------------------------
TITLE_MAP_GENE = {
    0: "all mRNAs",
    1: "CSR-1 target",
    2: "WAGO-1 target",
    8: "Germline target",
}

GROUP_LIST = {
    # "piRNA": [1, 2, 8],
    "piRNA": [2],
    "miRNA": [0]
}

MUT_ABBR = {"D": "del", "M": "subs", "A": "mut"}
pair_MUT = {"D": "deletion_readcount", "M": "mismatch_readcount", "A": "mutation_readcount"}
DN_LIST = ["all data", "high group", "middle group", "low group", "mid-low group"]

abu = {
    "piRNA": "22G",
    "miRNA": "mRNA"
}

score_col_dict = {
    ("miRNA", "true"): "max_mir_score",
    ("miRNA", "f"): "mir_score",
    ("piRNA", "true"): "max_pir_score",
    ("piRNA", "f"): "targeting_score",
}

rc_dict = {
    "true": "sum_readcount",
    "f": "read_count"
}

# -----------------------------------------------------------------------------
#  Small statistical wrappers
# -----------------------------------------------------------------------------
def KS_test(x, y):
    if len(x) > 0 and len(y) > 0:
        less = stats.mstats.ks_2samp(x, y, alternative='less')[1]
        greater = stats.mstats.ks_2samp(x, y, alternative='greater')[1]
        two_sided = stats.mstats.ks_2samp(x, y, alternative='two-sided')[1]
    else:
        less = greater = two_sided = 0
    return [two_sided, less, greater]

def T_test(x, y):
    if len(x) > 0 and len(y) > 0:
        d, two_sided = stats.ttest_ind(x, y, equal_var=False)
        if d < 0:
            greater = 1 - two_sided/2
            less = two_sided/2
        elif d > 0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
    return [two_sided, greater, less]

def U_test(x, y):
    if len(x) > 0 and len(y) > 0:
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided/2
            less = two_sided/2
        elif d > 0:
            greater = two_sided/2
            less = 1 - two_sided/2
        else:
            greater = less = two_sided/2
    else:
        less = greater = two_sided = 0
    return [two_sided, greater, less]

def permutation_test(data1, data2, num_permutations=1000, alternative='two-sided'):
    data1 = np.asarray(list(data1))
    data2 = np.asarray(list(data2))
    observed_diff = np.mean(data1) - np.mean(data2)
    combined = np.concatenate([data1, data2])
    count = 0
    n1 = len(data1)

    for _ in range(int(num_permutations)):
        np.random.shuffle(combined)
        new_data1 = combined[:n1]
        new_data2 = combined[n1:]
        new_diff = np.mean(new_data1) - np.mean(new_data2)

        if alternative == 'two-sided':
            if abs(new_diff) >= abs(observed_diff):
                count += 1
        elif alternative == 'greater':
            if new_diff >= observed_diff:
                count += 1
        elif alternative == 'less':
            if new_diff <= observed_diff:
                count += 1

    p_value = count / float(num_permutations)
    return observed_diff, p_value

# -----------------------------------------------------------------------------
#  IO helpers
# -----------------------------------------------------------------------------
def load_data(input_csv):
    return pd.read_csv(input_csv)

def read_target_excel(path):
    # 盡量用 openpyxl（適用於 .xlsx）；若環境未裝則退回預設引擎
    try:
        import openpyxl  # noqa: F401
        return pd.read_excel(path, engine="openpyxl")
    except Exception:
        return pd.read_excel(path)

def add_gene_id_if_needed(df, id_map_csv):
    if "Gene ID" not in df.columns:
        id_map = pd.read_csv(id_map_csv).rename(columns={"Gene name": "transcript_name"})
        df = pd.merge(df, id_map, on="transcript_name", how="inner")
    return df

# -----------------------------------------------------------------------------
#  Fold change helpers
# -----------------------------------------------------------------------------
def compute_fold_change(df, mode="piRNA"):
    if mode == 'miRNA':
        mut_cols = ['alg1.r1.rpkm', 'alg1.r2.rpkm', 'alg1.r3.rpkm']
        wt_cols = ['wt.r1.rpkm', 'wt.r2.rpkm', 'wt.r3.rpkm']
        min_vals = df[mut_cols + wt_cols].replace(0, np.nan).min()
        alpha = float(min_vals.min())

        df['avg_mRNA_WT_sum'] = df[wt_cols].sum(axis=1)/3.0 + alpha
        df['avg_mRNA_MUT_sum'] = df[mut_cols].sum(axis=1)/3.0 + alpha
        df['fold_change'] = np.log2(df['avg_mRNA_MUT_sum'] / df['avg_mRNA_WT_sum'])

        def _fc_wo(row):
            wt = row['avg_mRNA_WT_sum'] - alpha
            mut = row['avg_mRNA_MUT_sum'] - alpha
            if wt > 0 and mut > 0:
                return math.log(mut / wt, 2)
            return np.nan
        df['fold_change_without0'] = df.apply(_fc_wo, axis=1)

    else:
        alpha1 = df["22G_rc_WT"].replace(0, np.nan).min()
        alpha2 = df["22G_rc_MUT"].replace(0, np.nan).min()
        alpha = float(min(alpha1, alpha2))

        df["22G_rc_WT_alpha"] = df["22G_rc_WT"] + alpha
        df["22G_rc_MUT_alpha"] = df["22G_rc_MUT"] + alpha
        df["fold_change"] = np.log2(df["22G_rc_MUT_alpha"] / df["22G_rc_WT_alpha"])

        def _fc_wo(row):
            if row["22G_rc_WT"] and row["22G_rc_MUT"]:
                return math.log(row["22G_rc_MUT"] / row["22G_rc_WT"], 2)
            return np.nan
        df["fold_change_without0"] = df.apply(_fc_wo, axis=1)

    return df

def compute_pairing_fold_change(df, mode="piRNA"):
    if mode == 'miRNA':
        mut_cols = ['weight_alg1.r1.rpkm', 'weight_alg1.r2.rpkm', 'weight_alg1.r3.rpkm']
        wt_cols = ['weight_wt.r1.rpkm', 'weight_wt.r2.rpkm', 'weight_wt.r3.rpkm']
        min_vals = df[mut_cols + wt_cols].replace(0, np.nan).min()
        alpha = float(min_vals.min())

        df['avg_mRNA_WT_sum'] = df[wt_cols].sum(axis=1)/3.0 + alpha
        df['avg_mRNA_MUT_sum'] = df[mut_cols].sum(axis=1)/3.0 + alpha
        df['fold_change'] = np.log2(df['avg_mRNA_MUT_sum'] / df['avg_mRNA_WT_sum'])

        def _fc_wo(row):
            wt = row['avg_mRNA_WT_sum'] - alpha
            mut = row['avg_mRNA_MUT_sum'] - alpha
            if wt > 0 and mut > 0:
                return math.log(mut / wt, 2)
            return np.nan
        df['fold_change_without0'] = df.apply(_fc_wo, axis=1)

    else:
        alpha1 = df["weight_avg_22g_wt"].replace(0, np.nan).min()
        alpha2 = df["weight_avg_22g_mut"].replace(0, np.nan).min()
        alpha = float(min(alpha1, alpha2))

        df["weight_avg_22g_wt_alpha"] = df["weight_avg_22g_wt"] + alpha
        df["weight_avg_22g_mut_alpha"] = df["weight_avg_22g_mut"] + alpha
        df["fold_change"] = np.log2(df["weight_avg_22g_mut_alpha"] / df["weight_avg_22g_wt_alpha"])

        def _fc_wo(row):
            if row["weight_avg_22g_wt"] and row["weight_avg_22g_mut"]:
                return math.log(row["weight_avg_22g_mut"] / row["weight_avg_22g_wt"], 2)
            return np.nan
        df["fold_change_without0"] = df.apply(_fc_wo, axis=1)

    return df

# -----------------------------------------------------------------------------
#  Grouping helpers
# -----------------------------------------------------------------------------
def normalize_readcount(df, col="22G_rc_WT", factor=None):
    if factor is None:
        factor = df[col].sum() / 1000000.0
    return np.log10((df[col] / factor).replace(0, np.nan))

def split_data(df, column, boundaries=None, mode='manual'):
    """
    boundaries: (top, two_third, one_third, bot) 由大到小
    若 boundaries=None 或 mode='quantile'，改用 1/3、2/3 分位數自動分群。
    """
    all_df = df.copy()
    boundary_info = {}

    if mode == 'quantile' or boundaries is None:
        bot = all_df[column].min() - 1
        one_third = all_df[column].quantile(1/3.0)
        two_third = all_df[column].quantile(2/3.0)
        top = all_df[column].max() + 1

        boundary_info = {
            'min': bot,
            '1/3 quantile': one_third,
            '2/3 quantile': two_third,
            'max': top,
        }
    elif mode == 'manual':
        top, two_third, one_third, bot = boundaries
    else:
        raise ValueError("mode must be either 'manual' or 'quantile'")

    low = all_df[(all_df[column] <= one_third)]
    middle = all_df[(all_df[column] > one_third) & (all_df[column] <= two_third)]
    high = all_df[(all_df[column] > two_third) & (all_df[column] <= top)]
    midlow = all_df[(all_df[column] <= two_third)]

    return {
        "all data": all_df,
        "high group": high,
        "middle group": middle,
        "low group": low,
        "mid-low group": midlow,
        "boundaries": boundary_info,
    }

def assign_group_labels(series, boundaries, labels=("low", "mid", "high")):
    top, two_third, one_third, bot = boundaries
    group_labels = []
    for val in series:
        if val <= one_third:
            group_labels.append(labels[0])
        elif one_third < val <= two_third:
            group_labels.append(labels[1])
        elif two_third < val <= top:
            group_labels.append(labels[2])
        else:
            group_labels.append(None)
    return pd.Series(group_labels, index=series.index)

def add_two_mRNA_list(new, target, gene):
    if gene == 0:
        return new
    elif gene == 1:  # CSR-1 target
        csr1 = target[target['CSR1IP.N2__N2'] == True].reset_index(drop=True)
    elif gene == 2:  # WAGO-1 target
        csr1 = target[target['WAGO1IP__WAGO1Input'] == True].reset_index(drop=True)
    elif gene == 8:  # Germline target
        csr1 = target[target['ce.germline.genes.Ortiz.G3_2014.type'].notnull()].reset_index(drop=True)
    else:
        return new

    tmp1 = new[new['Gene ID'].isin(list(csr1['row_names']))].reset_index(drop=True)
    return tmp1

# --- Significance grouping helpers ---
def _coerce_bool_col(df, col):
    if col not in df.columns:
        raise ValueError("Missing required column: {}".format(col))
    def _to_bool(v):
        s = str(v).strip().lower()
        return (s == "true") or (s == "1") or (v is True) or (v == 1)
    return df[col].apply(_to_bool)

def split_by_significance(ana,
                          D_col="single_read_mut_stat_significance_D",
                          M_col="single_read_mut_stat_significance_M",
                          overlap_strategy="D_priority"):
    """
    產生三個 dataframe：
      - D 組：D True（預設 D 優先，若 D、M 同時 True 放進 D 組）
      - M 組：M True 且 D False（避免重複）
      - None 組：D False 且 M False
    overlap_strategy：'D_priority'、'M_priority'、'drop'
    **目前將overlap的參數忽略**
    """
    a = ana.copy()
    D = _coerce_bool_col(a, D_col)
    M = _coerce_bool_col(a, M_col)

    # if overlap_strategy == "drop":
    #     mask_D = D & (~M)
    #     mask_M = M & (~D)
    # elif overlap_strategy == "M_priority":
    #     mask_M = M
    #     mask_D = D & (~M)
    # else:  # "D_priority"
    #     mask_D = D
    #     mask_M = M & (~D)
    mask_D = D
    mask_M = M
    mask_none = (~D) & (~M)
    return a[mask_D], a[mask_M], a[mask_none]

def _grp_mask(series, level):
    """支援 'mid-low' 聯集，其餘為單一標籤精確比對。"""
    if level == "mid-low":
        return series.isin(["low", "mid"])
    return series == level

# -----------------------------------------------------------------------------
#  Plotting utilities
# -----------------------------------------------------------------------------
## without drawing data distribution
# def _plot_box_with_stats(ax, x, y, mut_n, title, ylabel, col1_name="", col2_name=""):
#     tmp1 = pd.Series(list(x))
#     tmp2 = pd.Series(list(y))
#     mean1 = round(float(np.mean(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
#     mean2 = round(float(np.mean(tmp2)), 3) if len(tmp2) > 0 else 'N/A'
#     median1 = round(float(np.median(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
#     median2 = round(float(np.median(tmp2)), 3) if len(tmp2) > 0 else 'N/A'

#     if col1_name != "" and col2_name != "":
#         col1_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col1_name, len(tmp1), mean1, median1)
#         col2_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col2_name, len(tmp2), mean2, median2)
#         tmp_df = pd.DataFrame({col1_label: tmp1, col2_label: tmp2})
#         sns.boxplot(
#             data=tmp_df, ax=ax, showfliers=True, width=0.5, linewidth=1, showmeans=False,
#             medianprops={"color": "orange"},
#             order=[col1_label, col2_label],
#             boxprops={"facecolor": "pink", "alpha": 0.7},
#             whiskerprops={"linestyle": "--"},
#         )
#     else:
#         tmp_df = pd.DataFrame({
#             "with {}\nN={}\nmean: {}\nmedian: {}".format(mut_n, len(tmp1), mean1, median1): tmp1,
#             "without mut\nN={}\nmean: {}\nmedian: {}".format(len(tmp2), mean2, median2): tmp2
#         })
#         sns.boxplot(
#             data=tmp_df, ax=ax, showfliers=True, width=0.5, linewidth=1, showmeans=False,
#             medianprops={"color": "orange"}, order=tmp_df.columns[::-1],
#             boxprops={"facecolor": "pink", "alpha": 0.7}, whiskerprops={"linestyle": "--"}
#         )

#     out = U_test(list(x), list(y))
#     U_m = np.format_float_scientific(out[1], precision=1)
#     U_c = np.format_float_scientific(out[2], precision=1)
#     out = T_test(list(x), list(y))
#     T_m = np.format_float_scientific(out[1], precision=1)
#     T_c = np.format_float_scientific(out[2], precision=1)
#     out = KS_test(list(x), list(y))
#     KS_m = np.format_float_scientific(out[1], precision=1)
#     KS_c = np.format_float_scientific(out[2], precision=1)
#     _, P_m = permutation_test(list(x), list(y), num_permutations=10000, alternative="greater")
#     _, P_c = permutation_test(list(x), list(y), num_permutations=10000, alternative="less")

#     if col1_name != "" and col2_name != "":
#         text = (
#             "U test: {} > {}: {}\n".format(col1_name, col2_name, U_m) +
#             "U test: {} < {}: {}\n".format(col1_name, col2_name, U_c) +
#             "-----------------------------------------\n" +
#             "T test: {} > {}: {}\n".format(col1_name, col2_name, T_m) +
#             "T test: {} < {}: {}\n".format(col1_name, col2_name, T_c) +
#             "-----------------------------------------\n" +
#             "KS test: {} > {}: {}\n".format(col1_name, col2_name, KS_m) +
#             "KS test: {} < {}: {}\n".format(col1_name, col2_name, KS_c) +
#             "-----------------------------------------\n" +
#             "Permutation test: {} > {}: {}\n".format(col1_name, col2_name, P_m) +
#             "Permutation test: {} < {}: {}".format(col1_name, col2_name, P_c)
#         )
#     else:
#         text = (
#             "U test: {} > {}: {}\n".format(mut_n, 'com', U_m) +
#             "U test: {} < {}: {}\n".format(mut_n, 'com', U_c) +
#             "-----------------------------------------\n" +
#             "T test: {} > {}: {}\n".format(mut_n, 'com', T_m) +
#             "T test: {} < {}: {}\n".format(mut_n, 'com', T_c) +
#             "-----------------------------------------\n" +
#             "KS test: {} > {}: {}\n".format(mut_n, 'com', KS_m) +
#             "KS test: {} < {}: {}\n".format(mut_n, 'com', KS_c) +
#             "-----------------------------------------\n" +
#             "Permutation test: {} > {}: {}\n".format(mut_n, 'com', P_m) +
#             "Permutation test: {} < {}: {}".format(mut_n, 'com', P_c)
#         )

#     ax.set_ylabel(ylabel, fontsize=12)
#     ax.set_title(title, fontsize=13)
#     ax.text(1.6, 0, text, fontsize=10, va='top')

#     result = {
#         "mut_name": mut_n,
#         "title": title,
#         "col1_name": col1_name if col1_name != "" else "with " + str(mut_n),
#         "col2_name": col2_name if col2_name != "" else "without mut",
#         "N_col1": int(len(x)),
#         "N_col2": int(len(y)),
#         "mean_col1": float(np.mean(x)) if len(x) > 0 else None,
#         "mean_col2": float(np.mean(y)) if len(y) > 0 else None,
#         "median_col1": float(np.median(x)) if len(x) > 0 else None,
#         "median_col2": float(np.median(y)) if len(y) > 0 else None,
#         "U_m": U_m, "U_c": U_c, "T_m": T_m, "T_c": T_c, "KS_m": KS_m, "KS_c": KS_c,
#         "P_m": P_m, "P_c": P_c
#     }
#     return result

# def fold_change_plot_with_and_without_0(
#     x1, y1, x2, y2, mut_n, title1, title2, out_png,
#     stats_records=None, col1_name="", col2_name=""
# ):
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
#     stat1 = _plot_box_with_stats(ax1, x1, y1, mut_n, title1, "log2(mut+α/wt+α)",
#                                  col1_name=col1_name, col2_name=col2_name)
#     stat2 = _plot_box_with_stats(ax2, x2, y2, mut_n, title2, "log2(mut/wt)",
#                                  col1_name=col1_name, col2_name=col2_name)
#     plt.tight_layout()
#     plt.savefig(out_png)
#     plt.clf()
#     plt.close()
#     gc.collect()
#     if stats_records is not None:
#         stats_records.append(stat1)
#         stats_records.append(stat2)

# def plot_single_box(x, y, mut_n, title, ylabel, out_png, col1_name="", col2_name="", stats_records=None):
#     fig, ax = plt.subplots(figsize=(8, 10))
#     stat1 = _plot_box_with_stats(ax, x, y, mut_n, title, ylabel, col1_name=col1_name, col2_name=col2_name)
#     plt.tight_layout()
#     plt.savefig(out_png)
#     plt.clf()
#     plt.close()
#     gc.collect()
#     if stats_records is not None:
#         stats_records.append(stat1)

# def plot_grouped_boxplot(groups, title, ylabel, out_png):
#     labeled_data = {}
#     for label, values in groups.items():
#         values = pd.Series(values)
#         mean = round(float(values.mean()), 3) if len(values) > 0 else 'N/A'
#         median = round(float(values.median()), 3) if len(values) > 0 else 'N/A'
#         labeled = "{}\nN={}\nmean: {}\nmedian: {}".format(label, len(values), mean, median)
#         labeled_data[labeled] = values

#     tmp_df = pd.DataFrame(labeled_data)
#     custom_order = ['low group', 'middle group', 'mid-low group', 'high group']
#     def get_base_label(colname):
#         return colname.split("\n")[0]
#     tmp_df = tmp_df[sorted(tmp_df.columns, key=lambda col: custom_order.index(get_base_label(col)))]

#     plt.figure(figsize=(7, 6))
#     sns.boxplot(
#         data=tmp_df, showfliers=True, width=0.5, linewidth=1, showmeans=False,
#         medianprops={"color": "orange"}, order=tmp_df.columns[::],
#         boxprops={"facecolor": "pink", "alpha": 0.7}, whiskerprops={"linestyle": "--"}
#     )
#     plt.ylabel(ylabel, fontsize=14)
#     plt.title(title, fontsize=13)
#     plt.tight_layout()
#     plt.savefig(out_png)
#     plt.clf()
#     plt.close()
#     gc.collect()

## violin plot
# =========================
# Core: 畫一張 violin（含統計文字），可選擇疊淡色點
# =========================
def _plot_violin_with_stats(
    ax, x, y, mut_n, title, ylabel, col1_name="", col2_name="",
    show_points=False,                 # 預設不疊點
    box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=6, point_alpha=0.20, mark_outliers=False,
    points_under_violin=True,         # True: 先畫點，讓點在小提琴下方
    violin_alpha=0.60,                # 小提琴透明度
    violin_inner="box",               # 'box'/'quartile'/None
    violin_cut=0,                     # 不向外延伸
    violin_bw=None,                   # 平滑帶寬（None=預設）
    violin_scale="width"              # 每組等寬
):
    import seaborn as sns
    tmp1 = pd.Series(list(x))
    tmp2 = pd.Series(list(y))
    mean1 = round(float(np.mean(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
    mean2 = round(float(np.mean(tmp2)), 3) if len(tmp2) > 0 else 'N/A'
    median1 = round(float(np.median(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
    median2 = round(float(np.median(tmp2)), 3) if len(tmp2) > 0 else 'N/A'

    # 組寬表 & 標籤
    if col1_name != "" and col2_name != "":
        col1_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col1_name, len(tmp1), mean1, median1)
        col2_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col2_name, len(tmp2), mean2, median2)
        wide_df = pd.DataFrame({col1_label: tmp1, col2_label: tmp2})
        order = [col1_label, col2_label]
    else:
        c1 = "with {}\nN={}\nmean: {}\nmedian: {}".format(mut_n, len(tmp1), mean1, median1)
        c2 = "without mut\nN={}\nmean: {}\nmedian: {}".format(len(tmp2), mean2, median2)
        wide_df = pd.DataFrame({c1: tmp1, c2: tmp2})
        order = list(wide_df.columns[::-1])

    # 先畫點（在小提琴下方）
    if show_points and points_under_violin:
        _add_point_cloud(
            ax, order, wide_df,
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=1
        )

    # 畫 violin（使用單色系，之後手動設 alpha）
    # 註：舊版 seaborn 用 'bw'；新版本是 'bw' 或 'bw_adjust'。為了相容性只用 'bw'。
    v = sns.violinplot(
        data=wide_df, order=order, ax=ax, width=box_width,
        inner=violin_inner, cut=violin_cut, bw=violin_bw, scale=violin_scale,
        linewidth=1.2, saturation=1.0, color="pink"
    )
    # 統一調整透明度
    try:
        for coll in v.collections:
            coll.set_alpha(violin_alpha)
    except Exception:
        pass

    # 後畫點（想疊在上方時）
    if show_points and not points_under_violin:
        _add_point_cloud(
            ax, order, wide_df,
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=3
        )

    # ===== 統計檢定（沿用你的流程） =====
    out = U_test(list(x), list(y)); U_m = _fmt_sci(out[1], precision=1); U_c = _fmt_sci(out[2], precision=1)
    out = T_test(list(x), list(y)); T_m = _fmt_sci(out[1], precision=1); T_c = _fmt_sci(out[2], precision=1)
    out = KS_test(list(x), list(y)); KS_m = _fmt_sci(out[1], precision=1); KS_c = _fmt_sci(out[2], precision=1)
    _, P_m = permutation_test(list(x), list(y), num_permutations=10000, alternative="greater")
    _, P_c = permutation_test(list(x), list(y), num_permutations=10000, alternative="less")

    if col1_name != "" and col2_name != "":
        text = (
            "U test: {} > {}: {}\n".format(col1_name, col2_name, U_m) +
            "U test: {} < {}: {}\n".format(col1_name, col2_name, U_c) +
            "-----------------------------------------\n" +
            "T test: {} > {}: {}\n".format(col1_name, col2_name, T_m) +
            "T test: {} < {}: {}\n".format(col1_name, col2_name, T_c) +
            "-----------------------------------------\n" +
            "KS test: {} > {}: {}\n".format(col1_name, col2_name, KS_m) +
            "KS test: {} < {}: {}\n".format(col1_name, col2_name, KS_c) +
            "-----------------------------------------\n" +
            "Permutation test: {} > {}: {}\n".format(col1_name, col2_name, P_m) +
            "Permutation test: {} < {}: {}".format(col1_name, col2_name, P_c)
        )
    else:
        text = (
            "U test: {} > {}: {}\n".format(mut_n, 'com', U_m) +
            "U test: {} < {}: {}\n".format(mut_n, 'com', U_c) +
            "-----------------------------------------\n" +
            "T test: {} > {}: {}\n".format(mut_n, 'com', T_m) +
            "T test: {} < {}: {}\n".format(mut_n, 'com', T_c) +
            "-----------------------------------------\n" +
            "KS test: {} > {}: {}\n".format(mut_n, 'com', KS_m) +
            "KS test: {} < {}: {}\n".format(mut_n, 'com', KS_c) +
            "-----------------------------------------\n" +
            "Permutation test: {} > {}: {}\n".format(mut_n, 'com', P_m) +
            "Permutation test: {} < {}: {}".format(mut_n, 'com', P_c)
        )

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.text(1.6, 0, text, fontsize=10, va='top')

    return {
        "mut_name": mut_n, "title": title,
        "col1_name": col1_name if col1_name != "" else "with " + str(mut_n),
        "col2_name": col2_name if col2_name != "" else "without mut"
    }


# =========================
# Wrapper: 單張 violin
# =========================
def plot_single_violin(
    x, y, mut_n, title, ylabel, out_png, col1_name="", col2_name="", stats_records=None,
    show_points=False, box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=6, point_alpha=0.20, mark_outliers=False,
    points_under_violin=True, violin_alpha=0.60, violin_inner="box", violin_cut=0, violin_bw=None, violin_scale="width"
):
    fig, ax = plt.subplots(figsize=(8, 10))

    _ = _plot_violin_with_stats(
        ax, x, y, mut_n, title, ylabel, col1_name=col1_name, col2_name=col2_name,
        show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
        max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers,
        points_under_violin=points_under_violin, violin_alpha=violin_alpha,
        violin_inner=violin_inner, violin_cut=violin_cut, violin_bw=violin_bw, violin_scale=violin_scale
    )

    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf(); plt.close(); gc.collect()
    if stats_records is not None:
        stats_records.append(_)


# =========================
# Wrapper: 上下兩張 violin（含 with/without 0）
# =========================
def fold_change_violin_with_and_without_0(
    x1, y1, x2, y2, mut_n, title1, title2, out_png,
    stats_records=None, col1_name="", col2_name="",
    show_points=False, box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=6, point_alpha=0.20, mark_outliers=False,
    points_under_violin=True, violin_alpha=0.60, violin_inner="box", violin_cut=0, violin_bw=None, violin_scale="width"
):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    _plot_violin_with_stats(
        ax1, x1, y1, mut_n, title1, "log2(mut+α/wt+α)",
        col1_name=col1_name, col2_name=col2_name,
        show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
        max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers,
        points_under_violin=points_under_violin, violin_alpha=violin_alpha,
        violin_inner=violin_inner, violin_cut=violin_cut, violin_bw=violin_bw, violin_scale=violin_scale
    )
    _plot_violin_with_stats(
        ax2, x2, y2, mut_n, title2, "log2(mut/wt)",
        col1_name=col1_name, col2_name=col2_name,
        show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
        max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers,
        points_under_violin=points_under_violin, violin_alpha=violin_alpha,
        violin_inner=violin_inner, violin_cut=violin_cut, violin_bw=violin_bw, violin_scale=violin_scale
    )

    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf(); plt.close(); gc.collect()


# =========================
# Wrapper: 多群組 violin
# =========================
def plot_grouped_violin(
    groups, title, ylabel, out_png,
    show_points=False, box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=5, point_alpha=0.20, mark_outliers=False,
    points_under_violin=True, violin_alpha=0.60, violin_inner="box", violin_cut=0, violin_bw=None, violin_scale="width"
):
    # 寬表 + 注記
    labeled_data = {}
    for label, values in groups.items():
        values = pd.Series(values)
        mean = round(float(values.mean()), 3) if len(values) > 0 else 'N/A'
        median = round(float(values.median()), 3) if len(values) > 0 else 'N/A'
        labeled = "{}\nN={}\nmean: {}\nmedian: {}".format(label, len(values), mean, median)
        labeled_data[labeled] = values
    wide_df = pd.DataFrame(labeled_data)

    custom_order = ['low group', 'middle group', 'mid-low group', 'high group']
    def _base(col):
        return col.split("\n")[0]
    order = sorted(wide_df.columns, key=lambda c: custom_order.index(_base(c)))

    fig, ax = plt.subplots(figsize=(7, 6))

    if show_points and points_under_violin:
        _add_point_cloud(
            ax, order, wide_df[order],
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=1
        )

    v = sns.violinplot(
        data=wide_df[order], order=order, ax=ax, width=box_width,
        inner=violin_inner, cut=violin_cut, bw=violin_bw, scale=violin_scale,
        linewidth=1.2, saturation=1.0, color="pink"
    )
    try:
        for coll in v.collections:
            coll.set_alpha(violin_alpha)
    except Exception:
        pass

    if show_points and not points_under_violin:
        _add_point_cloud(
            ax, order, wide_df[order],
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=3
        )

    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=13)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf(); plt.close(); gc.collect()

## drawing data point beside boxplot v2

# # ---- Python 3.5 / numpy 後備：科學記號格式化 ----
# def _fmt_sci(x, precision=1):
#     try:
#         return np.format_float_scientific(x, precision=precision)
#     except Exception:
#         # 後備，用一般格式化
#         return ("{0:." + str(precision) + "e}").format(x)


# # =========================
# # Helper: 疊加點雲（資料分佈）
# # =========================
# def _add_point_cloud(
#     ax, order, wide_df,
#     box_width=0.5,          # 與 seaborn.boxplot 的 width 對應
#     side="left",            # "left" 或 "right"：點雲擺在箱體哪一側
#     sep=0.07,               # 與箱體的水平間距
#     jitter=0.03,            # 點的水平抖動幅度
#     max_points=2000,        # 每個群組最多繪製的點數（大樣本抽樣以避免過慢）
#     point_size=6,
#     point_alpha=0.65,
#     mark_outliers=False     # False：離群值與一般點同樣式；True：離群值用黑菱形
# ):
#     """
#     將每一欄（群組）的原始資料以點雲形式疊加在箱型圖旁邊，避免與箱體重疊。
#     wide_df：寬表格（每欄一個群組），order：要顯示的欄名順序。
#     """
#     x_positions = np.arange(len(order))         # seaborn box 的 x 位置（0,1,2,...）
#     half = box_width / 2.0                      # 箱體半寬，用來計算外側中心點

#     for i, col in enumerate(order):
#         vals_full = wide_df[col].dropna().values
#         if len(vals_full) == 0:
#             continue

#         # 大樣本抽樣以加速繪圖
#         if len(vals_full) > max_points:
#             idx = np.random.choice(len(vals_full), size=max_points, replace=False)
#             vals = vals_full[idx]
#         else:
#             vals = vals_full

#         # 點雲中心位置：放在箱體外側，並留出 sep 的空隙
#         if side == "left":
#             center_x = x_positions[i] - (half + sep)
#         elif side == "right":
#             center_x = x_positions[i] + (half + sep)
#         else:  # fallback：中央（通常不建議）
#             center_x = x_positions[i]

#         xj = np.random.normal(loc=center_x, scale=jitter, size=len(vals))

#         if mark_outliers:
#             # 只有需要特別標示時才計算離群值
#             q1, q3 = np.percentile(vals_full, [25, 75])
#             iqr = q3 - q1
#             lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
#             is_out = (vals < lo) | (vals > hi)
#             # 一般點
#             ax.plot(
#                 xj[~is_out], vals[~is_out],
#                 marker='o', linestyle='None', ms=point_size, alpha=point_alpha, zorder=3
#             )
#             # 離群值（黑菱形）
#             if np.any(is_out):
#                 ax.plot(
#                     xj[is_out], vals[is_out],
#                     marker='D', linestyle='None', ms=point_size+1, color='black', alpha=1.0, zorder=4
#                 )
#         else:
#             # 全部同樣式（含離群值）
#             ax.plot(
#                 xj, vals, marker='o', linestyle='None',
#                 ms=point_size, alpha=point_alpha, zorder=3
#             )

#     # 給一點水平邊界，避免最邊的點被裁切
#     ax.margins(x=0.10)


# # =========================
# # Core: 畫一張 boxplot（含統計文字），可疊點雲
# # =========================
# def _plot_box_with_stats(
#     ax, x, y, mut_n, title, ylabel, col1_name="", col2_name="",
#     # 點雲參數（預設啟用、放在左側且不重疊）
#     show_points=False,
#     box_width=0.5, side="left", sep=0.08, jitter=0.03,
#     max_points=2000, point_size=6, point_alpha=0.65, mark_outliers=False
# ):
#     tmp1 = pd.Series(list(x))
#     tmp2 = pd.Series(list(y))
#     mean1 = round(float(np.mean(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
#     mean2 = round(float(np.mean(tmp2)), 3) if len(tmp2) > 0 else 'N/A'
#     median1 = round(float(np.median(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
#     median2 = round(float(np.median(tmp2)), 3) if len(tmp2) > 0 else 'N/A'

#     # 建立寬表 + 標籤
#     if col1_name != "" and col2_name != "":
#         col1_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col1_name, len(tmp1), mean1, median1)
#         col2_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col2_name, len(tmp2), mean2, median2)
#         wide_df = pd.DataFrame({col1_label: tmp1, col2_label: tmp2})
#         order = [col1_label, col2_label]
#     else:
#         c1 = "with {}\nN={}\nmean: {}\nmedian: {}".format(mut_n, len(tmp1), mean1, median1)
#         c2 = "without mut\nN={}\nmean: {}\nmedian: {}".format(len(tmp2), mean2, median2)
#         wide_df = pd.DataFrame({c1: tmp1, c2: tmp2})
#         order = list(wide_df.columns[::-1])  # 保持你原本的顯示順序

#     # Boxplot（關鍵：showfliers=False，避免與點雲重複）
#     sns.boxplot(
#         data=wide_df, ax=ax, showfliers=False,
#         width=box_width, linewidth=1, showmeans=False,
#         medianprops={"color": "orange"}, order=order,
#         boxprops={"facecolor": "pink", "alpha": 0.7},
#         whiskerprops={"linestyle": "--"}
#     )

#     # 疊上點雲（在箱體外側）
#     if show_points:
#         _add_point_cloud(
#             ax, order, wide_df,
#             box_width=box_width, side=side, sep=sep, jitter=jitter,
#             max_points=max_points, point_size=point_size,
#             point_alpha=point_alpha, mark_outliers=mark_outliers
#         )

#     # ===== 統計檢定（沿用你的既有流程） =====
#     out = U_test(list(x), list(y))
#     U_m = _fmt_sci(out[1], precision=1)
#     U_c = _fmt_sci(out[2], precision=1)

#     out = T_test(list(x), list(y))
#     T_m = _fmt_sci(out[1], precision=1)
#     T_c = _fmt_sci(out[2], precision=1)

#     out = KS_test(list(x), list(y))
#     KS_m = _fmt_sci(out[1], precision=1)
#     KS_c = _fmt_sci(out[2], precision=1)

#     _, P_m = permutation_test(list(x), list(y), num_permutations=10000, alternative="greater")
#     _, P_c = permutation_test(list(x), list(y), num_permutations=10000, alternative="less")

#     if col1_name != "" and col2_name != "":
#         text = (
#             "U test: {} > {}: {}\n".format(col1_name, col2_name, U_m) +
#             "U test: {} < {}: {}\n".format(col1_name, col2_name, U_c) +
#             "-----------------------------------------\n" +
#             "T test: {} > {}: {}\n".format(col1_name, col2_name, T_m) +
#             "T test: {} < {}: {}\n".format(col1_name, col2_name, T_c) +
#             "-----------------------------------------\n" +
#             "KS test: {} > {}: {}\n".format(col1_name, col2_name, KS_m) +
#             "KS test: {} < {}: {}\n".format(col1_name, col2_name, KS_c) +
#             "-----------------------------------------\n" +
#             "Permutation test: {} > {}: {}\n".format(col1_name, col2_name, P_m) +
#             "Permutation test: {} < {}: {}".format(col1_name, col2_name, P_c)
#         )
#     else:
#         text = (
#             "U test: {} > {}: {}\n".format(mut_n, 'com', U_m) +
#             "U test: {} < {}: {}\n".format(mut_n, 'com', U_c) +
#             "-----------------------------------------\n" +
#             "T test: {} > {}: {}\n".format(mut_n, 'com', T_m) +
#             "T test: {} < {}: {}\n".format(mut_n, 'com', T_c) +
#             "-----------------------------------------\n" +
#             "KS test: {} > {}: {}\n".format(mut_n, 'com', KS_m) +
#             "KS test: {} < {}: {}\n".format(mut_n, 'com', KS_c) +
#             "-----------------------------------------\n" +
#             "Permutation test: {} > {}: {}\n".format(mut_n, 'com', P_m) +
#             "Permutation test: {} < {}: {}".format(mut_n, 'com', P_c)
#         )

#     ax.set_ylabel(ylabel, fontsize=12)
#     ax.set_title(title, fontsize=13)
#     # 右側註記（維持你原本大致位置；若需要可再調整 x 座標）
#     ax.text(1.6, 0, text, fontsize=10, va='top')

#     result = {
#         "mut_name": mut_n,
#         "title": title,
#         "col1_name": col1_name if col1_name != "" else "with " + str(mut_n),
#         "col2_name": col2_name if col2_name != "" else "without mut",
#         "N_col1": int(len(x)),
#         "N_col2": int(len(y)),
#         "mean_col1": float(np.mean(x)) if len(x) > 0 else None,
#         "mean_col2": float(np.mean(y)) if len(y) > 0 else None,
#         "median_col1": float(np.median(x)) if len(x) > 0 else None,
#         "median_col2": float(np.median(y)) if len(y) > 0 else None,
#         "U_m": U_m, "U_c": U_c, "T_m": T_m, "T_c": T_c, "KS_m": KS_m, "KS_c": KS_c,
#         "P_m": P_m, "P_c": P_c
#     }
#     return result


# # =========================
# # Wrapper: 上下兩張 boxplot
# # =========================
# def fold_change_plot_with_and_without_0(
#     x1, y1, x2, y2, mut_n, title1, title2, out_png,
#     stats_records=None, col1_name="", col2_name="",
#     # 點雲參數（可視需要調整）
#     show_points=False, box_width=0.5, side="left", sep=0.08, jitter=0.03,
#     max_points=2000, point_size=6, point_alpha=0.65, mark_outliers=False
# ):
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

#     stat1 = _plot_box_with_stats(
#         ax1, x1, y1, mut_n, title1, "log2(mut+α/wt+α)",
#         col1_name=col1_name, col2_name=col2_name,
#         show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
#         max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers
#     )
#     stat2 = _plot_box_with_stats(
#         ax2, x2, y2, mut_n, title2, "log2(mut/wt)",
#         col1_name=col1_name, col2_name=col2_name,
#         show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
#         max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers
#     )

#     plt.tight_layout()
#     plt.savefig(out_png)
#     plt.clf(); plt.close(); gc.collect()

#     if stats_records is not None:
#         stats_records.append(stat1)
#         stats_records.append(stat2)


# # =========================
# # Wrapper: 單張圖
# # =========================
# def plot_single_box(
#     x, y, mut_n, title, ylabel, out_png, col1_name="", col2_name="", stats_records=None,
#     show_points=False, box_width=0.5, side="left", sep=0.08, jitter=0.03,
#     max_points=2000, point_size=6, point_alpha=0.65, mark_outliers=False
# ):
#     fig, ax = plt.subplots(figsize=(8, 10))

#     stat1 = _plot_box_with_stats(
#         ax, x, y, mut_n, title, ylabel,
#         col1_name=col1_name, col2_name=col2_name,
#         show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
#         max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers
#     )

#     plt.tight_layout()
#     plt.savefig(out_png)
#     plt.clf(); plt.close(); gc.collect()

#     if stats_records is not None:
#         stats_records.append(stat1)


# # =========================
# # Wrapper: 多群組 boxplot
# # =========================
# def plot_grouped_boxplot(
#     groups, title, ylabel, out_png,
#     show_points=False, box_width=0.5, side="left", sep=0.08, jitter=0.03,
#     max_points=2000, point_size=5, point_alpha=0.65, mark_outliers=False
# ):
#     # 準備寬表 + 標籤（含 N/mean/median）
#     labeled_data = {}
#     for label, values in groups.items():
#         values = pd.Series(values)
#         mean = round(float(values.mean()), 3) if len(values) > 0 else 'N/A'
#         median = round(float(values.median()), 3) if len(values) > 0 else 'N/A'
#         labeled = "{}\nN={}\nmean: {}\nmedian: {}".format(label, len(values), mean, median)
#         labeled_data[labeled] = values

#     wide_df = pd.DataFrame(labeled_data)

#     # 自訂排序（依你的既有邏輯）
#     custom_order = ['low group', 'middle group', 'mid-low group', 'high group']
#     def _base(col):
#         return col.split("\n")[0]

#     order = sorted(wide_df.columns, key=lambda c: custom_order.index(_base(c)))

#     fig, ax = plt.subplots(figsize=(7, 6))

#     sns.boxplot(
#         data=wide_df[order], ax=ax, showfliers=False,
#         width=box_width, linewidth=1, showmeans=False,
#         medianprops={"color": "orange"}, order=order,
#         boxprops={"facecolor": "pink", "alpha": 0.7},
#         whiskerprops={"linestyle": "--"}
#     )

#     if show_points:
#         _add_point_cloud(
#             ax, order, wide_df[order],
#             box_width=box_width, side=side, sep=sep, jitter=jitter,
#             max_points=max_points, point_size=point_size,
#             point_alpha=point_alpha, mark_outliers=mark_outliers
#         )

#     ax.set_ylabel(ylabel, fontsize=14)
#     ax.set_title(title, fontsize=13)
#     plt.tight_layout()
#     plt.savefig(out_png)
#     plt.clf(); plt.close(); gc.collect()

# v3
# ---- Python 3.5 / numpy 後備：科學記號格式化 ----
def _fmt_sci(x, precision=1):
    try:
        return np.format_float_scientific(x, precision=precision)
    except Exception:
        # 後備，用一般格式化
        return ("{0:." + str(precision) + "e}").format(x)


# =========================
# Helper: 疊加點雲（資料分佈）
# =========================
def _add_point_cloud(
    ax, order, wide_df,
    box_width=0.5,          # 與 seaborn.boxplot 的 width 對應
    side="overlap",         # "overlap"（與箱體重疊）或 "left"/"right"
    sep=0.0,                # 重疊時不需間距
    jitter=0.06,            # 點的水平抖動幅度（稍大一點讓密度更可讀）
    max_points=2000,        # 每個群組最多繪製的點數（抽樣避免過慢）
    point_size=6,
    point_alpha=0.20,       # 關鍵：非常淡，靠疊加顯示密度
    point_color="black",    # 單色＋低透明度，疊色自然變深
    mark_outliers=False,    # False：離群值與一般點同樣式；True：離群值用黑菱形
    zorder_points=3,        # 預設畫在箱體上方；若想壓在箱體下方，傳 1
    rasterized=True         # 向量輸出時點陣化點雲以減少檔案大小
):
    """
    將每一欄（群組）的原始資料以點雲形式疊加在箱型圖。
    side='overlap' 時點雲與箱體重疊；'left'/'right' 則貼側邊。
    """
    x_positions = np.arange(len(order))         # seaborn box 的 x 位置（0,1,2,...）
    half = box_width / 2.0                      # 箱體半寬，用來計算外側中心點

    for i, col in enumerate(order):
        vals_full = wide_df[col].dropna().values
        if len(vals_full) == 0:
            continue

        # 大樣本抽樣以加速繪圖
        if len(vals_full) > max_points:
            idx = np.random.choice(len(vals_full), size=max_points, replace=False)
            vals = vals_full[idx]
        else:
            vals = vals_full

        # 點雲中心位置：overlap=箱體中央，或放在箱體外側
        if side == "left":
            center_x = x_positions[i] - (half + max(sep, 0.02))
        elif side == "right":
            center_x = x_positions[i] + (half + max(sep, 0.02))
        else:  # 'overlap' 或其他 -> 正中
            center_x = x_positions[i]

        xj = np.random.normal(loc=center_x, scale=jitter, size=len(vals))

        if mark_outliers:
            # 只有需要特別標示時才計算離群值（以整體分佈估 IQR）
            q1, q3 = np.percentile(vals_full, [25, 75])
            iqr = q3 - q1
            lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
            is_out = (vals < lo) | (vals > hi)

            # 一般點（很淡）
            ax.scatter(
                xj[~is_out], vals[~is_out],
                s=point_size, alpha=point_alpha,
                c=point_color, edgecolors='none',
                zorder=zorder_points, rasterized=rasterized
            )
            # 離群值（略醒目）
            if np.any(is_out):
                ax.scatter(
                    xj[is_out], vals[is_out],
                    s=point_size + 6, alpha=min(1.0, point_alpha * 4.0),
                    c="black", marker='D', linewidths=0,
                    zorder=zorder_points + 1, rasterized=rasterized
                )
        else:
            # 全部同樣式（含離群值）
            ax.scatter(
                xj, vals, s=point_size, alpha=point_alpha,
                c=point_color, edgecolors='none',
                zorder=zorder_points, rasterized=rasterized
            )

    # 給一點水平邊界即可（重疊時不需要太大）
    ax.margins(x=0.05)


# =========================
# Core: 畫一張 boxplot（含統計文字），可疊點雲
# =========================
def _plot_box_with_stats(
    ax, x, y, mut_n, title, ylabel, col1_name="", col2_name="",
    # 點雲參數（預設啟用、與箱體重疊且畫在箱體下方）
    show_points=False,
    box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=6, point_alpha=0.20, mark_outliers=False,
    points_under_box=True   # True -> 先畫點，讓點在箱體下方；False -> 後畫點（在上方）
):
    tmp1 = pd.Series(list(x))
    tmp2 = pd.Series(list(y))
    mean1 = round(float(np.mean(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
    mean2 = round(float(np.mean(tmp2)), 3) if len(tmp2) > 0 else 'N/A'
    median1 = round(float(np.median(tmp1)), 3) if len(tmp1) > 0 else 'N/A'
    median2 = round(float(np.median(tmp2)), 3) if len(tmp2) > 0 else 'N/A'

    # 建立寬表 + 標籤
    if col1_name != "" and col2_name != "":
        col1_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col1_name, len(tmp1), mean1, median1)
        col2_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col2_name, len(tmp2), mean2, median2)
        wide_df = pd.DataFrame({col1_label: tmp1, col2_label: tmp2})
        order = [col1_label, col2_label]
    else:
        c1 = "with {}\nN={}\nmean: {}\nmedian: {}".format(mut_n, len(tmp1), mean1, median1)
        c2 = "without mut\nN={}\nmean: {}\nmedian: {}".format(len(tmp2), mean2, median2)
        wide_df = pd.DataFrame({c1: tmp1, c2: tmp2})
        order = list(wide_df.columns[::-1])  # 保持原本顯示順序

    # 先畫點（若要把點壓在箱體下方）
    if show_points and points_under_box:
        _add_point_cloud(
            ax, order, wide_df,
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=1
        )

    # Boxplot（關鍵：showfliers=False，避免與點雲重複；箱體略透明）
    sns.boxplot(
        data=wide_df, ax=ax, showfliers=False,
        width=box_width, linewidth=1.4, showmeans=False,
        medianprops={"color": "orange", "linewidth": 2.0}, order=order,
        boxprops={"facecolor": "pink", "alpha": 0.60},
        whiskerprops={"linestyle": "--", "linewidth": 1.2}
    )

    # 再畫點（若要把點疊在箱體上方、但很淡）
    if show_points and not points_under_box:
        _add_point_cloud(
            ax, order, wide_df,
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=3
        )

    # ===== 統計檢定（沿用既有流程；以下四個函式需由專案提供） =====
    out = U_test(list(x), list(y))
    U_m = _fmt_sci(out[1], precision=1)
    U_c = _fmt_sci(out[2], precision=1)

    out = T_test(list(x), list(y))
    T_m = _fmt_sci(out[1], precision=1)
    T_c = _fmt_sci(out[2], precision=1)

    out = KS_test(list(x), list(y))
    KS_m = _fmt_sci(out[1], precision=1)
    KS_c = _fmt_sci(out[2], precision=1)

    _, P_m = permutation_test(list(x), list(y), num_permutations=10000, alternative="greater")
    _, P_c = permutation_test(list(x), list(y), num_permutations=10000, alternative="less")

    if col1_name != "" and col2_name != "":
        text = (
            "U test: {} > {}: {}\n".format(col1_name, col2_name, U_m) +
            "U test: {} < {}: {}\n".format(col1_name, col2_name, U_c) +
            "-----------------------------------------\n" +
            "T test: {} > {}: {}\n".format(col1_name, col2_name, T_m) +
            "T test: {} < {}: {}\n".format(col1_name, col2_name, T_c) +
            "-----------------------------------------\n" +
            "KS test: {} > {}: {}\n".format(col1_name, col2_name, KS_m) +
            "KS test: {} < {}: {}\n".format(col1_name, col2_name, KS_c) +
            "-----------------------------------------\n" +
            "Permutation test: {} > {}: {}\n".format(col1_name, col2_name, P_m) +
            "Permutation test: {} < {}: {}".format(col1_name, col2_name, P_c)
        )
    else:
        text = (
            "U test: {} > {}: {}\n".format(mut_n, 'com', U_m) +
            "U test: {} < {}: {}\n".format(mut_n, 'com', U_c) +
            "-----------------------------------------\n" +
            "T test: {} > {}: {}\n".format(mut_n, 'com', T_m) +
            "T test: {} < {}: {}\n".format(mut_n, 'com', T_c) +
            "-----------------------------------------\n" +
            "KS test: {} > {}: {}\n".format(mut_n, 'com', KS_m) +
            "KS test: {} < {}: {}\n".format(mut_n, 'com', KS_c) +
            "-----------------------------------------\n" +
            "Permutation test: {} > {}: {}\n".format(mut_n, 'com', P_m) +
            "Permutation test: {} < {}: {}".format(mut_n, 'com', P_c)
        )

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13)
    # 右側註記（維持你原本大致位置；若需要可再調整 x 座標）
    ax.text(1.6, 0, text, fontsize=10, va='top')

    result = {
        "mut_name": mut_n,
        "title": title,
        "col1_name": col1_name if col1_name != "" else "with " + str(mut_n),
        "col2_name": col2_name if col2_name != "" else "without mut",
        "N_col1": int(len(x)),
        "N_col2": int(len(y)),
        "mean_col1": float(np.mean(x)) if len(x) > 0 else None,
        "mean_col2": float(np.mean(y)) if len(y) > 0 else None,
        "median_col1": float(np.median(x)) if len(x) > 0 else None,
        "median_col2": float(np.median(y)) if len(y) > 0 else None,
        "U_m": U_m, "U_c": U_c, "T_m": T_m, "T_c": T_c, "KS_m": KS_m, "KS_c": KS_c,
        "P_m": P_m, "P_c": P_c
    }
    return result


# =========================
# Wrapper: 上下兩張 boxplot
# =========================
def fold_change_plot_with_and_without_0(
    x1, y1, x2, y2, mut_n, title1, title2, out_png,
    stats_records=None, col1_name="", col2_name="",
    # 預設重疊、點淡、點在箱體下方
    show_points=False, box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=6, point_alpha=0.20, mark_outliers=False
):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    stat1 = _plot_box_with_stats(
        ax1, x1, y1, mut_n, title1, "log2(mut+α/wt+α)",
        col1_name=col1_name, col2_name=col2_name,
        show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
        max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers,
        points_under_box=True
    )
    stat2 = _plot_box_with_stats(
        ax2, x2, y2, mut_n, title2, "log2(mut/wt)",
        col1_name=col1_name, col2_name=col2_name,
        show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
        max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers,
        points_under_box=True
    )

    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf(); plt.close(); gc.collect()

    if stats_records is not None:
        stats_records.append(stat1)
        stats_records.append(stat2)


# =========================
# Wrapper: 單張圖
# =========================
def plot_single_box(
    x, y, mut_n, title, ylabel, out_png, col1_name="", col2_name="", stats_records=None,
    show_points=False, box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=6, point_alpha=0.20, mark_outliers=False
):
    fig, ax = plt.subplots(figsize=(8, 10))

    stat1 = _plot_box_with_stats(
        ax, x, y, mut_n, title, ylabel,
        col1_name=col1_name, col2_name=col2_name,
        show_points=show_points, box_width=box_width, side=side, sep=sep, jitter=jitter,
        max_points=max_points, point_size=point_size, point_alpha=point_alpha, mark_outliers=mark_outliers,
        points_under_box=True
    )

    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf(); plt.close(); gc.collect()

    if stats_records is not None:
        stats_records.append(stat1)


# =========================
# Wrapper: 多群組 boxplot
# =========================
def plot_grouped_boxplot(
    groups, title, ylabel, out_png,
    show_points=False, box_width=0.5, side="overlap", sep=0.0, jitter=0.06,
    max_points=2000, point_size=5, point_alpha=0.20, mark_outliers=False
):
    # 準備寬表 + 標籤（含 N/mean/median）
    labeled_data = {}
    for label, values in groups.items():
        values = pd.Series(values)
        mean = round(float(values.mean()), 3) if len(values) > 0 else 'N/A'
        median = round(float(values.median()), 3) if len(values) > 0 else 'N/A'
        labeled = "{}\nN={}\nmean: {}\nmedian: {}".format(label, len(values), mean, median)
        labeled_data[labeled] = values

    wide_df = pd.DataFrame(labeled_data)

    # 自訂排序（依你的既有邏輯）
    custom_order = ['low group', 'middle group', 'mid-low group', 'high group']
    def _base(col):
        return col.split("\n")[0]

    order = sorted(wide_df.columns, key=lambda c: custom_order.index(_base(c)))

    fig, ax = plt.subplots(figsize=(7, 6))

    sns.boxplot(
        data=wide_df[order], ax=ax, showfliers=False,
        width=box_width, linewidth=1.4, showmeans=False,
        medianprops={"color": "orange", "linewidth": 2.0}, order=order,
        boxprops={"facecolor": "pink", "alpha": 0.60},
        whiskerprops={"linestyle": "--", "linewidth": 1.2}
    )

    if show_points:
        _add_point_cloud(
            ax, order, wide_df[order],
            box_width=box_width, side=side, sep=sep, jitter=jitter,
            max_points=max_points, point_size=point_size,
            point_alpha=point_alpha, mark_outliers=mark_outliers,
            zorder_points=1  # 壓在箱體下方
        )

    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=13)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.clf(); plt.close(); gc.collect()

# -----------------------------------------------------------------------------
#  Main pipeline
# -----------------------------------------------------------------------------
def _run_pipeline(df, target, boundaries, group_by, score_type, out_dir,
                  mode="", RNA="piRNA", pairing=False, norm_factor=811.03,
                  use_significance=False, overlap_strategy="D_priority"):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    title_start = abu[RNA]

    for group in GROUP_LIST[RNA]:
    # for group in [2]:
        ana = add_two_mRNA_list(df, target, group)

        # 記數（並在 quantile 模式記錄邊界）
        split_all = split_data(ana, group_by, boundaries)
        log_filename = os.path.join(out_dir, "data_cnt_group{}_by_{}.log".format(group, group_by))
        with open(log_filename, "w", encoding="utf-8") as f:
            for group_name, df_group in split_all.items():
                if group_name == "boundaries":
                    binfo = split_all["boundaries"]
                    if binfo:
                        f.write("boundaries: {}\n".format(binfo))
                    continue
                f.write("{}: {} pairs\n".format(group_name, len(df_group)))
        del split_all

        # === 顯著性分群模式 ===
        if use_significance:
            D_df, M_df, N_df = split_by_significance(
                ana,
                D_col="single_read_mut_stat_significance_D",
                M_col="single_read_mut_stat_significance_M",
                overlap_strategy=overlap_strategy
            )
            splits_D = split_data(D_df, group_by, boundaries)
            splits_M = split_data(M_df, group_by, boundaries)
            splits_N = split_data(N_df, group_by, boundaries)

            for dn_label in DN_LIST:
                xD  = splits_D[dn_label]["fold_change"]
                xM  = splits_M[dn_label]["fold_change"]
                xN  = splits_N[dn_label]["fold_change"]
                xD2 = splits_D[dn_label].dropna()["fold_change_without0"]
                xM2 = splits_M[dn_label].dropna()["fold_change_without0"]
                xN2 = splits_N[dn_label].dropna()["fold_change_without0"]

                # D vs M
                # out_png = os.path.join(out_dir, "{}_FOLD_sig_group{}_D_vs_M_{}_by_{}.svg".format(
                #     title_start, group, dn_label.replace(" ", "_"), group_by))
                # fold_change_plot_with_and_without_0(
                #     xD, xM, xD2, xM2, mut_n="",
                #     title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     out_png=out_png, col1_name="D True", col2_name="M True"
                # )
                # D vs None
                # out_png = os.path.join(out_dir, "{}_FOLD_sig_group{}_D_vs_None_{}_by_{}.svg".format(
                #     title_start, group, dn_label.replace(" ", "_"), group_by))
                # fold_change_plot_with_and_without_0(
                #     xD, xN, xD2, xN2, mut_n="",
                #     title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     out_png=out_png, col1_name="D True", col2_name="Both False"
                # )
                out_png = os.path.join(out_dir, "{}_FOLD_sig_group{}_D_vs_None_{}_by_{}.svg".format(
                    title_start, group, dn_label.replace(" ", "_"), group_by))
                fold_change_plot_with_and_without_0(
                    xN, xD, xN2, xD2, mut_n="",
                    title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    out_png=out_png, col1_name="Both False", col2_name="D True"
                )
                # out_png = os.path.join(out_dir, "{}_FOLD_volin_sig_group{}_D_vs_None_{}_by_{}.svg".format(
                #     title_start, group, dn_label.replace(" ", "_"), group_by))
                # fold_change_violin_with_and_without_0(
                #     xN, xD, xN2, xD2, mut_n="",
                #     title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     out_png=out_png, col1_name="Both False", col2_name="D True"
                # )
                # M vs None
                # out_png = os.path.join(out_dir, "{}_FOLD_sig_group{}_M_vs_None_{}_by_{}.svg".format(
                #     title_start, group, dn_label.replace(" ", "_"), group_by))
                # fold_change_plot_with_and_without_0(
                #     xM, xN, xM2, xN2, mut_n="",
                #     title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     out_png=out_png, col1_name="M True", col2_name="Both False"
                # )
                out_png = os.path.join(out_dir, "{}_FOLD_sig_group{}_M_vs_None_{}_by_{}.svg".format(
                    title_start, group, dn_label.replace(" ", "_"), group_by))
                fold_change_plot_with_and_without_0(
                    xN, xM, xN2, xM2, mut_n="",
                    title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    out_png=out_png, col1_name="Both False", col2_name="M True"
                )

                # out_png = os.path.join(out_dir, "{}_FOLD_volin_sig_group{}_M_vs_None_{}_by_{}.svg".format(
                #     title_start, group, dn_label.replace(" ", "_"), group_by))
                # fold_change_violin_with_and_without_0(
                #     xN, xM, xN2, xM2, mut_n="",
                #     title1="{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                #     out_png=out_png, col1_name="Both False", col2_name="M True"
                # )

                # 若是 piRNA 再做 22G 讀數
                if RNA == "piRNA":
                    if pairing:
                        G22_D = [n/float(norm_factor) for n in splits_D[dn_label]["weight_avg_22g_wt"]]
                        G22_M = [n/float(norm_factor) for n in splits_M[dn_label]["weight_avg_22g_wt"]]
                        G22_N = [n/float(norm_factor) for n in splits_N[dn_label]["weight_avg_22g_wt"]]
                    else:
                        G22_D = [n/float(norm_factor) for n in splits_D[dn_label]["22G_rc_WT"]]
                        G22_M = [n/float(norm_factor) for n in splits_M[dn_label]["22G_rc_WT"]]
                        G22_N = [n/float(norm_factor) for n in splits_N[dn_label]["22G_rc_WT"]]

                    # plot_single_box(
                    #     x=G22_D, y=G22_M, mut_n="", ylabel="22G readcount",
                    #     title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    #     out_png=os.path.join(out_dir, "22G_sig_group{}_D_vs_M_{}_by_{}.svg".format(
                    #         group, dn_label.replace(" ", "_"), group_by)),
                    #     col1_name="D True", col2_name="M True"
                    # )
                    # plot_single_box(
                    #     x=G22_D, y=G22_N, mut_n="", ylabel="22G readcount",
                    #     title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    #     out_png=os.path.join(out_dir, "22G_sig_group{}_D_vs_None_{}_by_{}.svg".format(
                    #         group, dn_label.replace(" ", "_"), group_by)),
                    #     col1_name="D True", col2_name="Both False"
                    # )
                    # plot_single_box(
                    #     x=G22_M, y=G22_N, mut_n="", ylabel="22G readcount",
                    #     title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    #     out_png=os.path.join(out_dir, "22G_sig_group{}_M_vs_None_{}_by_{}.svg".format(
                    #         group, dn_label.replace(" ", "_"), group_by)),
                    #     col1_name="M True", col2_name="Both False"
                    # )
                    plot_single_box(
                        x=G22_N, y=G22_D, mut_n="", ylabel="22G readcount",
                        title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                        out_png=os.path.join(out_dir, "22G_sig_group{}_D_vs_None_{}_by_{}.svg".format(
                            group, dn_label.replace(" ", "_"), group_by)),
                        col1_name="Both False", col2_name="D True"
                    )

                    plot_single_box(
                        x=G22_N, y=G22_M, mut_n="", ylabel="22G readcount",
                        title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                        out_png=os.path.join(out_dir, "22G_sig_group{}_M_vs_None_{}_by_{}.svg".format(
                            group, dn_label.replace(" ", "_"), group_by)),
                        col1_name="Both False", col2_name="M True"
                    )

                    # plot_single_violin(
                    #     x=G22_N, y=G22_D, mut_n="", ylabel="22G readcount",
                    #     title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    #     out_png=os.path.join(out_dir, "22G_violin_sig_group{}_D_vs_None_{}_by_{}.svg".format(
                    #         group, dn_label.replace(" ", "_"), group_by)),
                    #     col1_name="Both False", col2_name="D True"
                    # )

                    # plot_single_violin(
                    #     x=G22_N, y=G22_M, mut_n="", ylabel="22G readcount",
                    #     title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    #     out_png=os.path.join(out_dir, "22G_violin_sig_group{}_M_vs_None_{}_by_{}.svg".format(
                    #         group, dn_label.replace(" ", "_"), group_by)),
                    #     col1_name="Both False", col2_name="M True"
                    # )
            # 顯著性分群完成，進入下一個 gene group
            continue

        # === 備用：維持舊的 with/without mutation 流程（當未啟用 sig_grouping） ===
        for mut_key, mut_name in MUT_ABBR.items():
            if mode == "pairing":
                with_mut = ana[ana[pair_MUT[mut_key]].astype(int) > 0]
                without_mut = ana[ana["mutation_readcount"].astype(int) == 0]
            else:
                with_mut = ana[~ana[mut_key].astype(str).isin(["[]"])]
                without_mut = ana[ana["A"].astype(str).isin(["[]"])]

            splits_with = split_data(with_mut, group_by, boundaries)
            splits_without = split_data(without_mut, group_by, boundaries)
            G22_abu_col = "weight_avg_22g_wt" if mode == "pairing" else "22G_rc_WT"
            if RNA == "piRNA":
                splits_all = split_data(ana, group_by, boundaries)
                high_group = [n/float(norm_factor) for n in splits_all["high group"][G22_abu_col]]
                mid_group  = [n/float(norm_factor) for n in splits_all["middle group"][G22_abu_col]]
                low_group  = [n/float(norm_factor) for n in splits_all["low group"][G22_abu_col]]
                midlow_group  = [n/float(norm_factor) for n in splits_all["mid-low group"][G22_abu_col]]
                out_png = os.path.join(out_dir, "22G_compare_group{}_by_{}.svg".format(group, group_by))
                # plot_grouped_boxplot(
                #     groups={"high group": high_group, "middle group": mid_group, "low group": low_group},
                #     title="22G level across {} groups ({})".format(group_by, TITLE_MAP_GENE[group]),
                #     ylabel="22G readcount", out_png=out_png
                # )
                # 三組兩兩比較
                # plot_single_box(low_group, mid_group, "", "22G level: low vs mid ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G readcount",
                #                 os.path.join(out_dir, "22G_low_vs_mid_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="low group", col2_name="mid group")
                # plot_single_box(mid_group, high_group, "", "22G level: mid vs high ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G readcount",
                #                 os.path.join(out_dir, "22G_mid_vs_high_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="mid group", col2_name="high group")
                # plot_single_box(low_group, high_group, "", "22G level: low vs high ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G readcount",
                #                 os.path.join(out_dir, "22G_low_vs_high_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="low group", col2_name="high group")
                # plot_single_box(midlow_group, high_group, "", "22G level: mid-low vs high ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G readcount",
                #                 os.path.join(out_dir, "22G_mid-low_vs_high_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="mid-low group", col2_name="high group")

                # fold change 的群間比較
                splits_all = split_data(ana, group_by, boundaries)
                high_group = splits_all["high group"]["fold_change"]
                mid_group  = splits_all["middle group"]["fold_change"]
                low_group  = splits_all["low group"]["fold_change"]
                midlow_group = splits_all["mid-low group"]["fold_change"]
                out_png = os.path.join(out_dir, "22G_foldchange_compare_group{}_by_{}.svg".format(group, group_by))
                # plot_grouped_boxplot(
                #     groups={"high group": high_group, "middle group": mid_group, "low group": low_group},
                #     title="22G foldchange across {} groups ({})".format(group_by, TITLE_MAP_GENE[group]),
                #     ylabel="22G foldchange", out_png=out_png
                # )
                # plot_single_box(low_group, mid_group, "", "22G foldchange: low vs mid ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G foldchange",
                #                 os.path.join(out_dir, "22G_foldchange_low_vs_mid_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="low group", col2_name="mid group")
                # plot_single_box(mid_group, high_group, "", "22G foldchange: mid vs high ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G foldchange",
                #                 os.path.join(out_dir, "22G_foldchange_mid_vs_high_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="mid group", col2_name="high group")
                # plot_single_box(low_group, high_group, "", "22G foldchange: low vs high ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G foldchange",
                #                 os.path.join(out_dir, "22G_foldchange_low_vs_high_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="low group", col2_name="high group")
                # plot_single_box(midlow_group, high_group, "", "22G foldchange: mid-low vs high ({})".format(TITLE_MAP_GENE[group]),
                #                 "22G foldchange",
                #                 os.path.join(out_dir, "22G_foldchange_mid-low_vs_high_group{}_by_{}.svg".format(group, group_by)),
                #                 col1_name="mid-low group", col2_name="high group")

            # with vs without per DN_LIST
            for dn_label in DN_LIST:
                x = splits_with[dn_label]["fold_change"]
                y = splits_without[dn_label]["fold_change"]
                x2 = splits_with[dn_label].dropna()["fold_change_without0"]
                y2 = splits_without[dn_label].dropna()["fold_change_without0"]
                out_png = os.path.join(out_dir, "FOLD_group{}_{}_{}_by_{}.svg".format(group, mut_key, dn_label.replace(" ", "_"), group_by))
                fold_change_plot_with_and_without_0(
                    x, y, x2, y2, mut_name,
                    "{} FOLD Change ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    "{} FOLD Change without 0 ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                    out_png
                )
                if RNA == "piRNA":
                    G22_x = [n/float(norm_factor) for n in splits_with[dn_label][G22_abu_col]]
                    G22_y = [n/float(norm_factor) for n in splits_without[dn_label][G22_abu_col]]
                    out_png_22G = os.path.join(out_dir, "22G_group{}_{}_{}_by_{}.svg".format(group, mut_key, dn_label.replace(" ", "_"), group_by))
                    plot_single_box(x=G22_x, y=G22_y, mut_n=mut_name,
                                    title="{} level ({}) - {}".format(title_start, dn_label, TITLE_MAP_GENE[group]),
                                    ylabel="22G readcount", out_png=out_png_22G)

# -----------------------------------------------------------------------------
#  Public API
# -----------------------------------------------------------------------------
def run_abu_score(input_csv, target_file, abu_region="10/0/-15/-30", score_col="score",
                  id_map_csv="../../data/reference/mRNA_WS275_IDtoName.csv",
                  out_dir="figure/G22_output/score", RNA="piRNA", norm_factor=811.03,
                  use_significance=False, overlap_strategy="D_priority", **kwargs):
    boundaries = tuple(float(x) for x in str(abu_region).split("/"))
    df = load_data(input_csv)
    df = add_gene_id_if_needed(df, id_map_csv)
    df = compute_fold_change(df, RNA)
    target = read_target_excel(target_file)
    _run_pipeline(df, target, boundaries, group_by=score_col, score_type=score_col, out_dir=out_dir,
                  RNA=RNA, norm_factor=float(norm_factor),
                  use_significance=bool(use_significance), overlap_strategy=str(overlap_strategy))

def run_abu_readcount(input_csv, target_file, abu_region="10000/4/1/0", readcount_col="read_count",
                      id_map_csv="../../data/reference/mRNA_WS275_IDtoName.csv",
                      out_dir="figure/G22_output/readcount", RNA="piRNA", norm_factor=811.03,
                      use_significance=False, overlap_strategy="D_priority", **kwargs):
    boundaries = tuple(float(x) for x in str(abu_region).split("/"))
    df = load_data(input_csv)
    df = add_gene_id_if_needed(df, id_map_csv)
    df = compute_fold_change(df, RNA)
    target = read_target_excel(target_file)
    _run_pipeline(df, target, boundaries, group_by=readcount_col, score_type=readcount_col, out_dir=out_dir,
                  RNA=RNA, norm_factor=float(norm_factor),
                  use_significance=bool(use_significance), overlap_strategy=str(overlap_strategy))

def run_pairing_abu_score(input_csv, target_file, abu_region="10/0/-15/-30", score_col="score",
                          id_map_csv="../../data/reference/mRNA_WS275_IDtoName.csv",
                          out_dir="figure/G22_output/score", RNA="piRNA", norm_factor=811.03, use_significance=True, **kwargs):
    boundaries = tuple(float(x) for x in str(abu_region).split("/"))
    df = load_data(input_csv)
    df = add_gene_id_if_needed(df, id_map_csv)
    df = compute_pairing_fold_change(df, RNA)
    target = read_target_excel(target_file)
    _run_pipeline(df, target, boundaries, group_by=score_col, score_type=score_col, out_dir=out_dir,
                  mode="pairing", RNA=RNA, norm_factor=float(norm_factor), pairing=True, use_significance=bool(use_significance))

def run_pairing_abu_readcount(input_csv, target_file, abu_region="10000/4/1/0", readcount_col="read_count",
                              id_map_csv="../../data/reference/mRNA_WS275_IDtoName.csv",
                              out_dir="figure/G22_output/readcount", RNA="piRNA", norm_factor=811.03, use_significance=True, **kwargs):
    boundaries = tuple(float(x) for x in str(abu_region).split("/"))
    df = load_data(input_csv)
    df = add_gene_id_if_needed(df, id_map_csv)
    df = compute_pairing_fold_change(df, RNA)
    target = read_target_excel(target_file)
    _run_pipeline(df, target, boundaries, group_by=readcount_col, score_type=readcount_col, out_dir=out_dir,
                  mode="pairing", RNA=RNA, norm_factor=float(norm_factor), pairing=True, use_significance=bool(use_significance))

# def run_combined(
#     input_csv, target_file,
#     score_col, readcount_col,
#     score_boundaries, rc_boundaries,
#     id_map_csv, out_dir,
#     RNA="piRNA", norm_factor=811.03, pairing=False, **kwargs
# ):
#     """
#     同時以 score 與 readcount 分箱（各 3 等分 + 'mid-low' 聯集），
#     在 4×4 = 16 個組合（含 mid-low）中做 with/without 比較並輸出圖與子集 CSV。
#     pairing=True 時使用 compute_pairing_fold_change，否則走 compute_fold_change。
#     """
#     # 讀資料 + 前處理
#     stats_records = []
#     title_start = abu[RNA]

#     df = load_data(input_csv)
#     df = add_gene_id_if_needed(df, id_map_csv)
#     df = compute_pairing_fold_change(df, RNA) if pairing else compute_fold_change(df, RNA)
#     target = read_target_excel(target_file)

#     sb = tuple(float(x) for x in str(score_boundaries).split("/"))
#     rb = tuple(float(x) for x in str(rc_boundaries).split("/"))

#     # 先在整體 df 上標好兩種分組欄位（low/mid/high）
#     df["score_grp"] = assign_group_labels(df[score_col], sb, labels=("low", "mid", "high"))
#     df["readcount_grp"] = assign_group_labels(df[readcount_col], rb, labels=("low", "mid", "high"))

#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)

#     group_levels = ("low", "mid", "mid-low", "high")

#     for group in GROUP_LIST[RNA]:
#         ana = add_two_mRNA_list(df, target, group)

#         for mut_key, mut_name in MUT_ABBR.items():
#             # with / without mutation 的切法依 pairing 與否而異
#             if pairing:
#                 with_mut = ana[ana[pair_MUT[mut_key]].astype(int) > 0]
#                 without_mut = ana[ana["mutation_readcount"].astype(int) == 0]
#             else:
#                 with_mut = ana[~ana[mut_key].astype(str).isin(["[]"])]
#                 without_mut = ana[ana["A"].astype(str).isin(["[]"])]

#             # 走 score×rc 的 16 個 cell（含 mid-low 聯集）
#             for sg in group_levels:
#                 for rg in group_levels:
#                     key = "{}_score_{}_rc".format(sg, rg)

#                     # 子集：with / without
#                     mask_s = _grp_mask(with_mut["score_grp"], sg)
#                     mask_r = _grp_mask(with_mut["readcount_grp"], rg)
#                     df_with = with_mut[mask_s & mask_r]

#                     mask_s_wo = _grp_mask(without_mut["score_grp"], sg)
#                     mask_r_wo = _grp_mask(without_mut["readcount_grp"], rg)
#                     df_wo = without_mut[mask_s_wo & mask_r_wo]

#                     # 匯出此 cell 的完整資料表（含有無變異混在一起，便於回溯）
#                     mask_s_all = _grp_mask(ana["score_grp"], sg)
#                     mask_r_all = _grp_mask(ana["readcount_grp"], rg)
#                     ana[mask_s_all & mask_r_all].to_csv(
#                         os.path.join(
#                             out_dir,
#                             "{}_FOLDChange_({})_{}.csv".format(title_start, TITLE_MAP_GENE[group], key)
#                         ),
#                         index=False
#                     )

#                     # 兩張 fold-change 圖（含 α 與去 0）
#                     out_png = os.path.join(
#                         out_dir,
#                         "{}_FOLD Change ({})_{}_{}.svg".format(
#                             title_start, TITLE_MAP_GENE[group], mut_name, key
#                         )
#                     )
#                     fold_change_plot_with_and_without_0(
#                         df_with["fold_change"],
#                         df_wo["fold_change"],
#                         df_with["fold_change_without0"].dropna(),
#                         df_wo["fold_change_without0"].dropna(),
#                         mut_name,
#                         "{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                         "{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                         out_png,
#                         stats_records
#                     )

#                     # piRNA 額外：22G 讀數
#                     if RNA == "piRNA":
#                         if pairing:
#                             G22_x = [n / float(norm_factor) for n in df_with["weight_avg_22g_wt"]]
#                             G22_y = [n / float(norm_factor) for n in df_wo["weight_avg_22g_wt"]]
#                         else:
#                             G22_x = [n / float(norm_factor) for n in df_with["22G_rc_WT"]]
#                             G22_y = [n / float(norm_factor) for n in df_wo["22G_rc_WT"]]

#                         out_png_22G = os.path.join(
#                             out_dir,
#                             "22G_group{}_({})_{}_{}.svg".format(group, TITLE_MAP_GENE[group], key, mut_name)
#                         )
#                         plot_single_box(
#                             x=G22_x, y=G22_y, mut_n=mut_name,
#                             title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                             ylabel="22G readcount",
#                             out_png=out_png_22G,
#                             stats_records=stats_records
#                         )

#     # 彙總統計
#     stats_df = pd.DataFrame(stats_records)
#     stats_df.to_csv(os.path.join(out_dir, "all_stats_summary.csv"), index=False)

# def run_combined(
#     input_csv, target_file,
#     score_col, readcount_col,
#     score_boundaries, rc_boundaries,
#     id_map_csv, out_dir,
#     RNA="piRNA", norm_factor=811.03, pairing=False,
#     overlap_strategy="union",  # 新增：D/M 交疊的解決策略，會傳給 split_by_significance
#     **kwargs
# ):
#     """
#     以 score 與 readcount 分箱（各 3 等分 + 'mid-low' 聯集），
#     在 4×4 = 16 個組合（含 mid-low）中，改用 single_read_mut_stat_significance_* 做分組，
#     針對每個 cell 輸出：
#       1) Both False (N) vs D True 的 fold-change 圖（含去 0 版本）
#       2) Both False (N) vs M True 的 fold-change 圖（含去 0 版本）
#       3) 若 RNA == 'piRNA'，再輸出 22G readcount 的盒型圖（同上兩種比較）
#     並保留每個 cell 的完整子表 CSV（含所有列，便於回溯），最後輸出 all_stats_summary.csv。

#     備註：
#     - D/M/None 的分群來自欄位：
#         single_read_mut_stat_significance_D, single_read_mut_stat_significance_M
#       其交疊由 overlap_strategy 控制，交由 split_by_significance() 處理。
#     """

#     # 讀資料 + 前處理
#     stats_records = []
#     title_start = abu[RNA]

#     df = load_data(input_csv)
#     df = add_gene_id_if_needed(df, id_map_csv)
#     df = compute_pairing_fold_change(df, RNA) if pairing else compute_fold_change(df, RNA)
#     target = read_target_excel(target_file)

#     # 邊界
#     sb = tuple(float(x) for x in str(score_boundaries).split("/"))
#     rb = tuple(float(x) for x in str(rc_boundaries).split("/"))

#     # 在整體 df 上標記分組欄位（low/mid/high；另外支援 mid-low 聯集）
#     df["score_grp"] = assign_group_labels(df[score_col], sb, labels=("low", "mid", "high"))
#     df["readcount_grp"] = assign_group_labels(df[readcount_col], rb, labels=("low", "mid", "high"))

#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)

#     group_levels = ("low", "mid", "mid-low", "high")

#     for group in GROUP_LIST[RNA]:
#         # 依基因群組挑資料
#         ana = add_two_mRNA_list(df, target, group)

#         # 走 score×rc 的 16 個 cell（含 mid-low 聯集）
#         for sg in group_levels:
#             for rg in group_levels:
#                 key = "{}_score_{}_rc".format(sg, rg)

#                 # 匯出此 cell 的完整資料表（含 D/M/None 混在一起，便於回溯）
#                 mask_s_all = _grp_mask(ana["score_grp"], sg)
#                 mask_r_all = _grp_mask(ana["readcount_grp"], rg)
#                 cell_all = ana[mask_s_all & mask_r_all]

#                 cell_csv = os.path.join(
#                     out_dir,
#                     "{}_FOLDChange_({})_{}.csv".format(title_start, TITLE_MAP_GENE[group], key)
#                 )
#                 cell_all.to_csv(cell_csv, index=False)

#                 # 用 single_read_mut_stat_significance 分群：D/M/None
#                 D_df, M_df, N_df = split_by_significance(
#                     cell_all,
#                     D_col="single_read_mut_stat_significance_D",
#                     M_col="single_read_mut_stat_significance_M",
#                     overlap_strategy=overlap_strategy
#                 )

#                 # --- Fold-change：N vs D ---
#                 xN  = N_df["fold_change"]
#                 xD  = D_df["fold_change"]
#                 xN2 = N_df["fold_change_without0"].dropna()
#                 xD2 = D_df["fold_change_without0"].dropna()

#                 out_png_D = os.path.join(
#                     out_dir,
#                     "{}_FOLD Change ({})_D_vs_None_{}.svg".format(
#                         title_start, TITLE_MAP_GENE[group], key
#                     )
#                 )
#                 fold_change_plot_with_and_without_0(
#                     xN, xD, xN2, xD2, mut_n="del",
#                     title1="{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                     title2="{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                     out_png=out_png_D, stats_records=stats_records,
#                     col1_name="Both False", col2_name="D True"
#                 )

#                 # --- Fold-change：N vs M ---
#                 xM  = M_df["fold_change"]
#                 xM2 = M_df["fold_change_without0"].dropna()

#                 out_png_M = os.path.join(
#                     out_dir,
#                     "{}_FOLD Change ({})_M_vs_None_{}.svg".format(
#                         title_start, TITLE_MAP_GENE[group], key
#                     )
#                 )
#                 fold_change_plot_with_and_without_0(
#                     xN, xM, xN2, xM2, mut_n="subs",
#                     title1="{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                     title2="{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                     out_png=out_png_M, stats_records=stats_records,
#                     col1_name="Both False", col2_name="M True"
#                 )

#                 # --- piRNA 額外：22G 讀數盒型圖（N vs D、N vs M） ---
#                 if RNA == "piRNA":
#                     if pairing:
#                         G22_N = [n/float(norm_factor) for n in N_df["weight_avg_22g_wt"]]
#                         G22_D = [n/float(norm_factor) for n in D_df["weight_avg_22g_wt"]]
#                         G22_M = [n/float(norm_factor) for n in M_df["weight_avg_22g_wt"]]
#                     else:
#                         G22_N = [n/float(norm_factor) for n in N_df["22G_rc_WT"]]
#                         G22_D = [n/float(norm_factor) for n in D_df["22G_rc_WT"]]
#                         G22_M = [n/float(norm_factor) for n in M_df["22G_rc_WT"]]

#                     out_png_22G_D = os.path.join(
#                         out_dir,
#                         "22G_group{}_({})_{}_D_vs_None.svg".format(group, TITLE_MAP_GENE[group], key)
#                     )
#                     plot_single_box(
#                         x=G22_N, y=G22_D, mut_n="del",
#                         title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                         ylabel="22G readcount",
#                         out_png=out_png_22G_D, stats_records=stats_records,
#                         col1_name="Both False", col2_name="D True"
#                     )

#                     out_png_22G_M = os.path.join(
#                         out_dir,
#                         "22G_group{}_({})_{}_M_vs_None.svg".format(group, TITLE_MAP_GENE[group], key)
#                     )
#                     plot_single_box(
#                         x=G22_N, y=G22_M, mut_n="subs",
#                         title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
#                         ylabel="22G readcount",
#                         out_png=out_png_22G_M, stats_records=stats_records,
#                         col1_name="Both False", col2_name="M True"
#                     )

#     # 彙總統計
#     stats_df = pd.DataFrame(stats_records)
#     stats_df.to_csv(os.path.join(out_dir, "all_stats_summary.csv"), index=False)

def run_combined(
    input_csv, target_file,
    score_col, readcount_col,
    score_boundaries, rc_boundaries,
    id_map_csv, out_dir,
    RNA="piRNA", norm_factor=811.03, pairing=False,
    overlap_strategy="union",
    **kwargs
):
    """
    score×readcount 分箱（含 'mid-low'），以 single_read_mut_stat_significance_* 做 D/M/None 分群，
    對每個 cell 輸出：
      - Fold-change：N vs D 與 N vs M 的 boxplot（含點雲）＋ violin（可選疊點）
      - 若 RNA == 'piRNA'：22G readcount 的 boxplot（含點雲）＋ violin
      - 該 cell 完整子表 CSV
    並彙整 all_stats_summary.csv。
    """
    stats_records = []
    title_start = abu[RNA]

    df = load_data(input_csv)
    df = add_gene_id_if_needed(df, id_map_csv)
    df = compute_pairing_fold_change(df, RNA) if pairing else compute_fold_change(df, RNA)
    target = read_target_excel(target_file)

    sb = tuple(float(x) for x in str(score_boundaries).split("/"))
    rb = tuple(float(x) for x in str(rc_boundaries).split("/"))

    df["score_grp"] = assign_group_labels(df[score_col], sb, labels=("low", "mid", "high"))
    df["readcount_grp"] = assign_group_labels(df[readcount_col], rb, labels=("low", "mid", "high"))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    group_levels = ("low", "mid", "mid-low", "high")

    for group in GROUP_LIST[RNA]:
        ana = add_two_mRNA_list(df, target, group)

        for sg in group_levels:
            for rg in group_levels:
                key = "{}_score_{}_rc".format(sg, rg)

                # 匯出該 cell 原始表
                mask_s_all = _grp_mask(ana["score_grp"], sg)
                mask_r_all = _grp_mask(ana["readcount_grp"], rg)
                cell_all = ana[mask_s_all & mask_r_all]
                cell_all.to_csv(
                    os.path.join(out_dir, "{}_FOLDChange_({})_{}.csv".format(
                        title_start, TITLE_MAP_GENE[group], key
                    )),
                    index=False
                )

                # 以顯著性切成 D/M/None
                D_df, M_df, N_df = split_by_significance(
                    cell_all,
                    D_col="single_read_mut_stat_significance_D",
                    M_col="single_read_mut_stat_significance_M",
                    overlap_strategy=overlap_strategy
                )

                # ---------- Fold-change：N vs D ----------
                xN  = N_df["fold_change"]
                xD  = D_df["fold_change"]
                xN2 = N_df["fold_change_without0"].dropna()
                xD2 = D_df["fold_change_without0"].dropna()

                # boxplot（含點雲）
                out_png_D_box = os.path.join(
                    out_dir,
                    "{}_FOLD Change ({})_D_vs_None_{}.svg".format(
                        title_start, TITLE_MAP_GENE[group], key
                    )
                )
                fold_change_plot_with_and_without_0(
                    xN, xD, xN2, xD2, mut_n="del",
                    title1="{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                    title2="{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                    out_png=out_png_D_box, stats_records=stats_records,
                    col1_name="Both False", col2_name="D True"
                )

                # violin（可選疊點；這裡預設不疊，想疊點可把 show_points=False）
                out_png_D_vio = os.path.join(
                    out_dir,
                    "{}_FOLD Violin ({})_D_vs_None_{}.svg".format(
                        title_start, TITLE_MAP_GENE[group], key
                    )
                )
                # fold_change_violin_with_and_without_0(
                #     xN, xD, xN2, xD2, mut_n="del",
                #     title1="{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                #     out_png=out_png_D_vio,
                #     col1_name="Both False", col2_name="D True"
                # )

                # ---------- Fold-change：N vs M ----------
                xM  = M_df["fold_change"]
                xM2 = M_df["fold_change_without0"].dropna()

                out_png_M_box = os.path.join(
                    out_dir,
                    "{}_FOLD Change ({})_M_vs_None_{}.svg".format(
                        title_start, TITLE_MAP_GENE[group], key
                    )
                )
                fold_change_plot_with_and_without_0(
                    xN, xM, xN2, xM2, mut_n="subs",
                    title1="{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                    title2="{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                    out_png=out_png_M_box, stats_records=stats_records,
                    col1_name="Both False", col2_name="M True"
                )

                out_png_M_vio = os.path.join(
                    out_dir,
                    "{}_FOLD Violin ({})_M_vs_None_{}.svg".format(
                        title_start, TITLE_MAP_GENE[group], key
                    )
                )
                # fold_change_violin_with_and_without_0(
                #     xN, xM, xN2, xM2, mut_n="subs",
                #     title1="{} FOLD Change ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                #     title2="{} FOLD Change without 0 ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                #     out_png=out_png_M_vio,
                #     col1_name="Both False", col2_name="M True"
                # )

                # ---------- 22G readcount（piRNA 才畫） ----------
                if RNA == "piRNA":
                    if pairing:
                        G22_N = [n/float(norm_factor) for n in N_df["weight_avg_22g_wt"]]
                        G22_D = [n/float(norm_factor) for n in D_df["weight_avg_22g_wt"]]
                        G22_M = [n/float(norm_factor) for n in M_df["weight_avg_22g_wt"]]
                    else:
                        G22_N = [n/float(norm_factor) for n in N_df["22G_rc_WT"]]
                        G22_D = [n/float(norm_factor) for n in D_df["22G_rc_WT"]]
                        G22_M = [n/float(norm_factor) for n in M_df["22G_rc_WT"]]

                    # boxplot（含點雲）
                    out_png_22G_D_box = os.path.join(
                        out_dir,
                        "22G_group{}_({})_{}_D_vs_None.svg".format(
                            group, TITLE_MAP_GENE[group], key
                        )
                    )
                    plot_single_box(
                        x=G22_N, y=G22_D, mut_n="del",
                        title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                        ylabel="22G readcount",
                        out_png=out_png_22G_D_box, stats_records=stats_records,
                        col1_name="Both False", col2_name="D True"
                    )

                    out_png_22G_M_box = os.path.join(
                        out_dir,
                        "22G_group{}_({})_{}_M_vs_None.svg".format(
                            group, TITLE_MAP_GENE[group], key
                        )
                    )
                    plot_single_box(
                        x=G22_N, y=G22_M, mut_n="subs",
                        title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                        ylabel="22G readcount",
                        out_png=out_png_22G_M_box, stats_records=stats_records,
                        col1_name="Both False", col2_name="M True"
                    )

                    # violin（這裡預設把原始點疊在小提琴下方；若不要點，把 show_points=False）
                    out_png_22G_D_vio = os.path.join(
                        out_dir,
                        "22G_violin_group{}_({})_{}_D_vs_None.svg".format(
                            group, TITLE_MAP_GENE[group], key
                        )
                    )
                    # plot_single_violin(
                    #     x=G22_N, y=G22_D, mut_n="del",
                    #     title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                    #     ylabel="22G readcount",
                    #     out_png=out_png_22G_D_vio, stats_records=stats_records,
                    #     col1_name="Both False", col2_name="D True",
                    #     show_points=False, points_under_violin=True
                    # )

                    out_png_22G_M_vio = os.path.join(
                        out_dir,
                        "22G_violin_group{}_({})_{}_M_vs_None.svg".format(
                            group, TITLE_MAP_GENE[group], key
                        )
                    )
                    # plot_single_violin(
                    #     x=G22_N, y=G22_M, mut_n="subs",
                    #     title="{} level ({} sc, {} rc) - {}".format(title_start, sg, rg, TITLE_MAP_GENE[group]),
                    #     ylabel="22G readcount",
                    #     out_png=out_png_22G_M_vio, stats_records=stats_records,
                    #     col1_name="Both False", col2_name="M True",
                    #     show_points=False, points_under_violin=True
                    # )

    # 彙總統計
    stats_df = pd.DataFrame(stats_records)
    stats_df.to_csv(os.path.join(out_dir, "all_stats_summary.csv"), index=False)


# -----------------------------------------------------------------------------
#  CLI
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="22G/mRNA fold-change analysis with optional significance grouping")

    parser.add_argument("--csv", required=True, help="input csv produced by snakemake")
    parser.add_argument("--target", required=True, help="add_two_HCLee.RNAseq.master.xlsx")
    parser.add_argument("--id_map", default="../../data/reference/mRNA_WS275_IDtoName.csv",
                        help="Gene ID to name mapping CSV")

    parser.add_argument("--mode", choices=["score", "readcount", "combined"], default="score")
    parser.add_argument("--boundaries", default="10/0/-15/-30",
                        help="top/2_3/1_3/bot boundaries, separated by '/'")
    parser.add_argument("--pairing", default="f",
                        help="'true' to use pairing-level pipeline for legacy modes")
    parser.add_argument("--RNA", default="piRNA", choices=["piRNA", "miRNA"])
    parser.add_argument("--norm_factor", default=811.03, type=float)
    parser.add_argument("--out_dir", default="figure/G22_output", help="output folder")

    # 額外：顯著性分群開關與策略
    parser.add_argument("--sig_grouping", action="store_true",
                        help="使用 single_read_mut_stat_significance_D/M 分三組作圖（D True / M True / both False）")
    parser.add_argument("--sig_overlap", choices=["D_priority", "M_priority", "drop"], default="D_priority",
                        help="D 與 M 同時為 True 的處理策略；預設 D_priority")

    # combined 模式參數（保留，未改動）
    parser.add_argument("--score_boundaries", default="10/0/-15/-30")
    parser.add_argument("--rc_boundaries", default="10000/4/1/0")

    args = parser.parse_args()

    score_col = score_col_dict.get((args.RNA, args.pairing.lower()), "targeting_score")
    rc_col = rc_dict.get(args.pairing.lower(), "read_count")

    use_sig = True
    overlap = str(args.sig_overlap)

    if args.mode == "combined":
        # raise NotImplementedError("combined 模式未在此版本整合顯著性分群；請改用 --mode score 或 --mode readcount。")
        run_combined(
            args.csv, args.target,
            score_col=score_col,
            readcount_col=rc_col,
            score_boundaries=args.score_boundaries,
            rc_boundaries=args.rc_boundaries,
            id_map_csv=args.id_map,
            out_dir=args.out_dir,
            RNA=args.RNA,
            norm_factor=args.norm_factor,
            pairing=(args.pairing.lower() == "true"),
            use_significance=True,                    # 強制用顯著性分群
            overlap_strategy=str(args.sig_overlap)    # 傳遞策略（function 目前允許重疊）
        )

    elif args.mode == "score" and args.pairing.lower() == "true":
        run_pairing_abu_score(args.csv, args.target, args.boundaries,
                              score_col=score_col, id_map_csv=args.id_map,
                              out_dir=args.out_dir, RNA=args.RNA, norm_factor=args.norm_factor)

    elif args.mode == "score":
        run_abu_score(args.csv, args.target, args.boundaries,
                      score_col=score_col, id_map_csv=args.id_map,
                      out_dir=args.out_dir, RNA=args.RNA, norm_factor=args.norm_factor,
                      use_significance=use_sig, overlap_strategy=overlap)

    elif args.mode == "readcount" and args.pairing.lower() == "true":
        run_pairing_abu_readcount(args.csv, args.target, args.boundaries,
                                  id_map_csv=args.id_map, readcount_col="sum_readcount",
                                  out_dir=args.out_dir, RNA=args.RNA, norm_factor=args.norm_factor)

    else:
        run_abu_readcount(args.csv, args.target, args.boundaries,
                          id_map_csv=args.id_map, readcount_col=rc_col,
                          out_dir=args.out_dir, RNA=args.RNA, norm_factor=args.norm_factor,
                          use_significance=use_sig, overlap_strategy=overlap)
