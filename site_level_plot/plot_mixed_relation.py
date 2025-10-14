"""
plot_mutation_relation.py   (Python 3.6 compatible)

功能
  1. 根據 species + mode 選擇分數欄位
  2. 根據 mut 類型產生 mutation (0/1) 及讀數欄位
  3. 產生圖表：
     A. 散點圖：metric vs. 突變比例
     B. Box-plot（含若干統計檢定）
     C. 進階：
        - mode=score    ：以 readcount 分成高/中/低，畫出「有/無 CIMS」的 score 分布
        - mode=readcount：以 score 分成高/中/低，畫出「有/無 CIMS」的 readcount 分布
        - mode=CIMS     ：在整體資料上，畫出「有/無 CIMS」的 score 與 readcount 分布
"""

from __future__ import print_function
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import ast

from scipy import stats
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.multitest import multipletests

# -------------------------------
# Helpers: statistics wrappers
# -------------------------------

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
            greater = 1 - two_sided / 2.0
            less = two_sided / 2.0
        elif d > 0:
            greater = two_sided / 2.0
            less = 1 - two_sided / 2.0
        else:
            greater = less = two_sided / 2.0
    else:
        less = greater = two_sided = 0
    return [two_sided, greater, less]

def U_test(x, y):
    if len(x) > 0 and len(y) > 0:
        d, two_sided = stats.ranksums(x, y)
        if d < 0:
            greater = 1 - two_sided / 2.0
            less = two_sided / 2.0
        elif d > 0:
            greater = two_sided / 2.0
            less = 1 - two_sided / 2.0
        else:
            greater = less = two_sided / 2.0
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
    p_value = float(count) / float(num_permutations)
    return observed_diff, p_value

# ----------------------------------
# Plot: box with stats text
# ----------------------------------

def _plot_box_with_stats(ax, x, y, mut_n, title, ylabel, col1_name="", col2_name=""):
    tmp1 = pd.Series(list(x))
    tmp2 = pd.Series(list(y))
    tmp1 = pd.to_numeric(tmp1, errors='coerce').dropna()
    tmp2 = pd.to_numeric(tmp2, errors='coerce').dropna()
    mean1 = round(np.mean(tmp1), 3) if len(tmp1) > 0 else 'N/A'
    mean2 = round(np.mean(tmp2), 3) if len(tmp2) > 0 else 'N/A'
    median1 = round(np.median(tmp1), 3) if len(tmp1) > 0 else 'N/A'
    median2 = round(np.median(tmp2), 3) if len(tmp2) > 0 else 'N/A'

    if col1_name != "" and col2_name != "":
        col1_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col1_name, len(tmp1), mean1, median1)
        col2_label = "{}\nN={}\nmean: {}\nmedian: {}".format(col2_name, len(tmp2), mean2, median2)
        tmp_df = pd.DataFrame({col2_label: tmp2, col1_label: tmp1})
        sns.boxplot(
            data=tmp_df,
            ax=ax,
            showfliers=False,
            width=0.5,
            linewidth=1,
            showmeans=False,
            medianprops={"color": "orange"},
            order=[col2_label, col1_label],  # put 'with' on the right
            boxprops={"facecolor": "pink", "alpha": 0.7},
            whiskerprops={"linestyle": "--"},
        )
    else:
        tmp_df = pd.DataFrame({
            "with {}\nN={}\nmean: {}\nmedian: {}".format(mut_n, len(tmp1), mean1, median1): tmp1,
            "without mut\nN={}\nmean: {}\nmedian: {}".format(len(tmp2), mean2, median2): tmp2,
        })
        sns.boxplot(
            data=tmp_df, ax=ax, showfliers=False, width=0.5, linewidth=1, showmeans=False,
            medianprops={"color": "orange"}, order=tmp_df.columns[::-1],
            boxprops={"facecolor": "pink", "alpha": 0.7}, whiskerprops={"linestyle": "--"}
        )

    # tests
    out = U_test(list(x), list(y))
    U_m = np.format_float_scientific(out[1], precision=1)
    U_c = np.format_float_scientific(out[2], precision=1)
    out = T_test(list(x), list(y))
    T_m = np.format_float_scientific(out[1], precision=1)
    T_c = np.format_float_scientific(out[2], precision=1)
    out = KS_test(list(x), list(y))
    KS_m = np.format_float_scientific(out[1], precision=1)
    KS_c = np.format_float_scientific(out[2], precision=1)
    _, P_m = permutation_test(list(x), list(y), num_permutations=10000, alternative="greater")
    _, P_c = permutation_test(list(x), list(y), num_permutations=10000, alternative="less")

    if col1_name != "" and col2_name != "":
        text = (
            "U test: {} > {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), U_m) +
            "U test: {} < {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), U_c) +
            "-----------------------------------------\n" +
            "T test: {} > {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), T_m) +
            "T test: {} < {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), T_c) +
            "-----------------------------------------\n" +
            "KS test: {} > {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), KS_m) +
            "KS test: {} < {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), KS_c) +
            "-----------------------------------------\n" +
            "Permutation test: {} > {}: {}\n".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), P_m) +
            "Permutation test: {} < {}: {}".format(col1_name.replace(" group", ""), col2_name.replace(" group", ""), P_c)
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
    ax.text(1.6, 0, text, fontsize=12, va='top')

    result = {
        "mut_name": mut_n,
        "title": title,
        "col1_name": col1_name if col1_name != "" else "with " + str(mut_n),
        "col2_name": col2_name if col2_name != "" else "without mut",
        "N_col1": len(x),
        "N_col2": len(y),
        "mean_col1": np.mean(x) if len(x) > 0 else None,
        "mean_col2": np.mean(y) if len(y) > 0 else None,
        "median_col1": np.median(x) if len(x) > 0 else None,
        "median_col2": np.median(y) if len(y) > 0 else None,
        "U_m": U_m, "U_c": U_c, "T_m": T_m, "T_c": T_c, "KS_m": KS_m, "KS_c": KS_c,
        "P_m": P_m, "P_c": P_c
    }
    return result

# -------------------------------
# Binning helpers
# -------------------------------

def parse_listlike(v):
    """Coerce a CSV cell to a Python list if possible; else empty list for NaNs/empties."""
    if isinstance(v, (list, tuple, set)):
        return list(v)
    if pd.isnull(v):
        return []
    s = v
    try:
        # try JSON/python literal list
        out = ast.literal_eval(str(s))
        if isinstance(out, (list, tuple, set)):
            return list(out)
    except Exception:
        pass
    # try semicolon/comma separated
    text = str(s).strip()
    if text in ("", "none", "None", "NA", "N/A", "[]"):
        return []
    if ";" in text:
        return [t for t in text.split(";") if t != ""]
    if "," in text:
        return [t for t in text.split(",") if t != ""]
    # fallback: treat as single token
    return [text]

def score_thresholds(species, series=None):
    # 固定門檻（pair 級）
    if species == "ALG-1" or species == "Homo":
        return 140.0, 100.0
    elif species == "PRG-1":
        return 0.0, -15.0
    # 其他（例如 Fly/Homo at hybrid）：若提供 series，使用分位數作為門檻
    if series is not None:
        try:
            q66, q33 = np.percentile(pd.Series(series).dropna(), [66.7, 33.3])
            return float(q66), float(q33)
        except Exception:
            pass
    # 後備：回傳 None（呼叫端應處理）
    return None, None

def readcount_thresholds():
    return 4.0, 1.0

def cut_three_groups(series, high, mid):
    cond_high = series > high
    cond_mid = (series <= high) & (series > mid)
    cond_low = series <= mid
    return np.select([cond_high, cond_mid, cond_low], ["high", "mid", "low"])  # dtype <U4

# -------------------------------
# Main
# -------------------------------

p = argparse.ArgumentParser(description="Visualise mutation vs. score/readcount/CIMS")
p.add_argument("--csv", required=True)
p.add_argument("--species", required=True, choices=["ALG-1", "PRG-1", "Fly", "Homo"])
p.add_argument("--mode", required=True, choices=["score", "readcount", "CIMS"])  # 新增 CIMS
p.add_argument("--mut", required=True, choices=["del", "all", "mis"])           # 要觀察的 CIMS 類型
p.add_argument("--output", required=True)
# pair 旗標保留，但預設也提供同名欄位的對應
p.add_argument('--pair', '-p', action='store_true', help='啟用 pair 模式（有給就為 True）')
args = p.parse_args()

# I/O
try:
    df = pd.read_csv(args.csv)
except Exception as e:
    sys.exit("read CSV failed: {}".format(e))
try:
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
except Exception as e:
    sys.exit("mkdir failed: {}".format(e))

# mutation 欄位對應
# mutation 欄位對應（pair vs hybrid）
if args.pair:
    mut_map = {
        "del": ("deletion_readcount", "Deletion readcount / Sum readcount"),
        "all": ("mutation_readcount", "All-mutation readcount / Sum readcount"),
        "mis": ("mismatch_readcount", "Mismatch readcount / Sum readcount"),
    }
else:
    # hybrid：欄位為 list（位置集合），後續以是否為空作為『是否有突變』
    mut_map = {
        "del": ("D", "Deletion"),
        "all": ("A", "All-mutation"),
        "mis": ("M", "Mismatch"),
    }
mut_col, mut_label = mut_map[args.mut]
mut_tag = "mut" if args.mut == "all" else args.mut
# hybrid 下，散點圖 y 標籤改為 presence（0/1）
if not args.pair:
    mut_label = mut_label + " (presence)"
if mut_col not in df.columns:
    sys.exit("'{}' missing".format(mut_col))

# metric 選擇
# Resolve metric column for pair vs hybrid
if args.pair:
    if args.mode == "readcount":
        metric_col = "sum_readcount"
    else:
        if args.species == "Fly":
            metric_col = "max_up_score"
        elif args.species  == "PRG-1":
            metric_col = "max_pir_score"
        else: # ALG-1 or Homo
            metric_col = "max_mir_score"
else:  # hybrid mode
    if args.mode == "readcount":
        metric_col = "read_count"
    else:  # CIMS/score use species-specific score at hybrid level
        if args.species == "PRG-1":
            metric_col = "targeting_score"
        elif args.species == "Fly":
            metric_col = "RNAup_score"
        else:  # ALG-1 or Homo
            metric_col = "mir_score"

if metric_col not in df.columns:
    sys.exit("'{}' missing".format(metric_col))

# 供後續方便：metric 欄位改名為 metric
# 同時決定「總讀數欄位」：pair 用 sum_readcount；hybrid 用 read_count
total_rc_col = "sum_readcount" if args.pair else "read_count"
if args.mode == "readcount" and metric_col != total_rc_col:
    sys.exit("metric_col/total_rc_col mismatch")

if metric_col == total_rc_col and total_rc_col not in df.columns:
    sys.exit("'{}' missing".format(total_rc_col))

Denominator = "metric" if metric_col == total_rc_col else total_rc_col

if total_rc_col not in df.columns:
    sys.exit("'{}' missing".format(total_rc_col))

if metric_col not in df.columns:
    sys.exit("'{}' missing".format(metric_col))

df = df.rename(columns={metric_col: "metric"})
if Denominator not in df.columns:
    sys.exit("'{}' missing".format(Denominator))

# ---- Mutation presence flags (pair vs hybrid) ----
if args.pair:
    df["_has_sel_mut"] = (df[mut_col] > 0).astype(int)
    if "mutation_readcount" in df.columns:
        df["_without_any_mut"] = (df["mutation_readcount"] == 0).astype(int)
    else:
        # fallback: derive from components if available
        comp = []
        for c in ["deletion_readcount", "mismatch_readcount", "mutation_readcount"]:
            if c in df.columns:
                comp.append(df[c].fillna(0))
        if len(comp) == 0:
            df["_without_any_mut"] = 1  # assume none
        else:
            s = sum(comp)
            df["_without_any_mut"] = (s == 0).astype(int)
else:
    # ensure D/M/A are lists
    for c in ["D", "M", "A"]:
        if c in df.columns:
            df[c] = df[c].apply(parse_listlike)
    def _nonempty_list(v):
        return 1 if isinstance(v, (list, tuple, set)) and len(v) > 0 else 0
    df["_has_sel_mut"] = df[mut_col].apply(_nonempty_list)
    if "A" in df.columns:
        df["_without_any_mut"] = df["A"].apply(_nonempty_list).map(lambda z: 1 if z == 0 else 0)
    else:
        d = df["D"].apply(_nonempty_list) if "D" in df.columns else 0
        m = df["M"].apply(_nonempty_list) if "M" in df.columns else 0
        df["_without_any_mut"] = ((d == 0) & (m == 0)).astype(int)

# 比例與相關（保留散點圖）
# if args.pair:
#     with np.errstate(divide='ignore', invalid='ignore'):
#         denom = df[Denominator].replace(0, np.nan)
#         ratio = df[mut_col] / denom
#         ratio = ratio.fillna(0.0)
# else:
#     # hybrid: use presence (0/1)
#     ratio = df["_has_sel_mut"].astype(float)

# try:
#     r2, p2 = pearsonr(df["metric"].values, ratio.values)
# except Exception:
#     r2, p2 = (np.nan, np.nan)

# corr2 = "Pearson r = {:.3f}".format(r2) if not np.isnan(r2) else "Pearson r = N/A"

# fig, ax = plt.subplots()
# ax.scatter(df["metric"], ratio, alpha=0.3, s=10)
# ax.set_xlabel(metric_col)
# ax.set_ylabel(mut_label)
# ax.set_title("Scatter of {} vs. {}".format(metric_col, mut_label))
# fig.tight_layout(rect=[0, 0.15, 1, 1])
# fig.text(0.5, 0.03, corr2, ha="center", va="bottom", fontsize=10)
# scatter_path = os.path.join(args.output, "scatter_{}_{}_{}.svg".format(args.species, metric_col, args.mut))
# fig.savefig(scatter_path, dpi=300, bbox_inches="tight")
# plt.close(fig)
# print("Scatter plot         →", scatter_path)

# 以 metric 三分組，並在各組做有/無該 CIMS 的比例檢定（延續舊版）
if args.mode in ("score", "readcount"):
    if args.mode == "score":
        # 舊版：以 score 定義三分
        bh, bm = score_thresholds(args.species)
    else:
        # 舊版：以 readcount 定義三分
        bh, bm = readcount_thresholds()
    df["group"] = cut_three_groups(df["metric"], bh, bm)
    df["has_mut"] = df["_has_sel_mut"]
    summary = df.groupby("group")["has_mut"].agg(["sum", "count"])
    group_stats = summary.rename(columns={"sum": "mutations", "count": "total"}).reset_index()
    group_stats["prop"] = group_stats["mutations"] / group_stats["total"].replace(0, np.nan)

    for _, row in group_stats.iterrows():
        print("{:<4s}  total = {:4d},  mutations = {:4d}  ({:.1%})".format(
            row["group"], int(row["total"]), int(row["mutations"]), 0.0 if pd.isnull(row["prop"]) else float(row["prop"]))
        )

    out_g = os.path.join(args.output, "group_summary_{}_{}_{}.csv".format(args.species, metric_col, args.mut))
    group_stats.to_csv(out_g, index=False)
    print("Group summary table       →", out_g)

    pairs = [("high", "mid"), ("high", "low"), ("mid", "low")]
    results = []
    for g1, g2 in pairs:
        if g1 in summary.index and g2 in summary.index:
            count = summary.loc[[g1, g2], "sum"].values
            nobs = summary.loc[[g1, g2], "count"].values
            try:
                z, p_raw = proportions_ztest(count, nobs, alternative="larger")
            except Exception:
                z, p_raw = (np.nan, np.nan)
            prop1, prop2 = (count / nobs.astype(float)).tolist()
            results.append({"pair": "{} vs {}".format(g1, g2), "prop1": prop1, "prop2": prop2, "z": z, "p_raw": p_raw})
    res_df = pd.DataFrame(results)
    if len(res_df) > 0:
        reject, p_adj, _, _ = multipletests(res_df["p_raw"].fillna(1.0), method="bonferroni")
        res_df["p_adj"] = p_adj
        res_df["signif (α=0.05)"] = reject
        out_t = os.path.join(args.output, "prop_test_{}_{}_{}.csv".format(args.species, metric_col, args.mut))
        res_df.to_csv(out_t, index=False)
        print("Pairwise proportion test table →", out_t)

# --------------------------------------------
# 主題 boxplot（舊版：整體 with/without；修正 y 定義）
# --------------------------------------------
try:
    x = df.loc[df["_has_sel_mut"] == 1, "metric"]
    y = df.loc[df["_without_any_mut"] == 1, "metric"]  # 修正：依所選 mut 類型
    fig, ax = plt.subplots(figsize=(10, 10))
    mut_tag = "mut" if args.mut == "all" else args.mut
    _ = _plot_box_with_stats(
        ax, x, y,
        mut_n='with {}'.format(mut_tag),
        title='{} boxplot'.format("metric" if metric_col == "sum_readcount" else metric_col),
        ylabel=metric_col,
        col1_name='with {}'.format(mut_tag),
        col2_name='without mut'
    )
    out_box = os.path.join(args.output, "boxplot_comp_{}_{}_{}.svg".format(args.species, metric_col, args.mut))
    plt.tight_layout()
    fig.savefig(out_box, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Boxplot (overall with/without) →", out_box)
except Exception as e:
    # 取出 traceback 物件（相容 Python 3.5）
    tb = sys.exc_info()[2]
    extracted = traceback.extract_tb(tb)
    last = extracted[-1] if extracted else None

    if last:
        filename, lineno, funcname, code = last
        print("[warn] overall boxplot failed: {msg} @ {file}:{line} in {func}".format(
            msg=e, file=filename, line=lineno, func=funcname))
        if code:
            print("    Code: {code}".format(code=code))
    else:
        # 後備：沒有可用的堆疊就只印訊息
        print("[warn] overall boxplot failed:", e)


# ---------------------------------------------------------
# 進階圖：依使用者需求新增的三種情境
# ---------------------------------------------------------

# A) mode=score → 以 readcount 三分組，為每個組別各自畫一張 boxplot（with/without 所選 CIMS 的 score 分布）
if args.mode == "score":
    rh, rm = readcount_thresholds()
    if total_rc_col not in df.columns:
        print("[warn] missing '{}' for score-mode grouping; skip per-group boxplots".format(total_rc_col))
    else:
        df["rc_group"] = cut_three_groups(df[total_rc_col], rh, rm)
        df["has_mut"] = df["_has_sel_mut"]
        order = ["high", "mid", "low"]
        for grp in order:
            sub = df[df["rc_group"] == grp]
            xg = sub.loc[sub["_has_sel_mut"] == 1, "metric"]
            yg = sub.loc[sub["_without_any_mut"] == 1, "metric"]
            if len(xg) == 0 and len(yg) == 0:
                print("[warn] rc_group='{}' has no data; skip".format(grp))
                continue
            fig, ax = plt.subplots(figsize=(10, 8))
            _plot_box_with_stats(
                ax, xg, yg,
                mut_n='with {}'.format(mut_tag),
                title='Score boxplot — readcount group: {} (thr {})'.format(grp, (rh, rm)),
                ylabel=metric_col,
                col1_name='with {}'.format(mut_tag),
                col2_name='without mut'
            )
            out_grp = os.path.join(
                args.output,
                'box_score_with_without_rc-{}_{}_{}.svg'.format(grp, args.species, args.mut)
            )
            plt.tight_layout()
            fig.savefig(out_grp, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("Per-group score boxplot (rc={}) →".format(grp), out_grp)

# B) mode=readcount → 以 score 三分組，為每個組別各自畫一張 boxplot（with/without 所選 CIMS 的 readcount 分布）
elif args.mode == "readcount":
    # 以 score 三分組，使用對應資料層級的 score 欄位
    if args.pair:
        if args.species == "PRG-1":
            score_group_col = "max_pir_score"
        elif args.species ==  "Fly":
            score_group_col = "max_up_score"
        else:
            score_group_col = "max_mir_score"
        rc_col = "sum_readcount"
    else:
        if args.species == "PRG-1":
            score_group_col = "targeting_score"
        elif args.species == "Fly":
            score_group_col = "RNAup_score"
        else:
            score_group_col = "mir_score"
    if score_group_col not in df.columns:
        print("[warn] missing '{}' for readcount-mode grouping; skip per-group boxplots".format(score_group_col))
    else:
        sh, sm = score_thresholds(args.species, df[score_group_col])
        if sh is None or sm is None:
            print("[warn] score thresholds unavailable; skip per-group boxplots")
        else:
            df["score_group"] = cut_three_groups(df[score_group_col], sh, sm)
            df["has_mut"] = df["_has_sel_mut"]
            order = ["high", "mid", "low"]
            for grp in order:
                sub = df[df["score_group"] == grp]
                xg = sub.loc[sub["_has_sel_mut"] == 1, "metric"]   # readcount（rename 為 metric）
                yg = sub.loc[sub["_without_any_mut"] == 1, "metric"]
                if len(xg) == 0 and len(yg) == 0:
                    print("[warn] score_group='{}' has no data; skip".format(grp))
                    continue
                fig, ax = plt.subplots(figsize=(10, 8))
                _plot_box_with_stats(
                    ax, xg, yg,
                    mut_n='with {}'.format(mut_tag),
                    title='Readcount boxplot — score group: {} (thr {})'.format(grp, (sh, sm)),
                    ylabel=total_rc_col,
                    col1_name='with {}'.format(mut_tag),
                    col2_name='without mut'
                )
                out_grp = os.path.join(
                    args.output,
                    'box_readcount_with_without_score-{}_{}_{}.svg'.format(grp, args.species, args.mut)
                )
                plt.tight_layout()
                fig.savefig(out_grp, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print("Per-group readcount boxplot (score={}) →".format(grp), out_grp)

# C) mode=CIMS → 整體上畫出有/無該 CIMS 的 score 與 readcount 分布 → 整體上畫出有/無該 CIMS 的 score 與 readcount 分布
elif args.mode == "CIMS":
    # score/readcount 欄位依層級決定
    if args.pair:
        if args.species == "PRG-1":
            score_col = "max_pir_score"
        elif args.species ==  "Fly":
            score_col = "max_up_score"
        else:
            score_col = "max_mir_score"
        rc_col = "sum_readcount"
    else:
        if args.species == "PRG-1":
            score_col = "targeting_score"
        elif args.species == "Fly":
            score_col = "RNAup_score"
        else:
            score_col = "mir_score"
        rc_col = "read_count"
    if score_col not in df.columns or rc_col not in df.columns:
        print("[warn] missing score or readcount columns for CIMS mode")
    else:
        for value_col, title, fname in [
            (score_col, "Score distribution (with/without {})".format(mut_tag), "box_score_with_without_{}_{}.svg"),
            (rc_col, "Readcount distribution (with/without {})".format(mut_tag), "box_readcount_with_without_{}_{}.svg"),
        ]:
            fig, ax = plt.subplots(figsize=(9, 6))
            x1 = df.loc[df["_has_sel_mut"] == 1, value_col]
            y1 = df.loc[df["_without_any_mut"] == 1, value_col]
            _plot_box_with_stats(
                ax, x1, y1,
                mut_n='with {}'.format(mut_tag),
                title=title,
                ylabel=value_col,
                col1_name='with {}'.format(mut_tag),
                col2_name='without mut'
            )
            outp = os.path.join(args.output, fname.format(args.species, args.mut))
            plt.tight_layout()
            fig.savefig(outp, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("CIMS mode boxplot →", outp)

# print("CIMS mode boxplot →", outp)

print("Done.")
