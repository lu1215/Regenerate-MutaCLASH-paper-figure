import argparse
import re
import os
import pandas as pd

def extract_parenthesis(text):
    m = re.search(r'\(([^)]*)\)', str(text))
    return m.group(1) if m else ""

def main():
    parser = argparse.ArgumentParser(description="依照 target（可選）、RNA 種類等條件，拆分並處理輸入 CSV，最後輸出 8 個子檔案。")
    parser.add_argument("-i", "--input", required=True, help="輸入的 CSV 檔路徑")
    parser.add_argument("-p", "--prefix", required=True, help="輸出檔案的前綴名稱")
    parser.add_argument("-r", "--rna", required=True, choices=["miRNA", "piRNA"], help="指定 RNA 種類")
    parser.add_argument("-t", "--target", default="", help="title 欄位要包含的目標字串")
    parser.add_argument("-o", "--output_dir", default=".", help="輸出檔案存放的資料夾（預設為目前資料夾）")
    parser.add_argument(
        "--suff_rev",
        action="store_true",      # 有帶這個旗標就會是 True
        default=False,            # 預設 False（沒帶旗標）
        help="反轉 p-value 方向/欄位 suffix 的選擇（使用 _c 取代 _m），並在輸出欄名標示 (common > mut)。"
    )
    args = parser.parse_args()

    # 建立輸出資料夾（若不存在）
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    df = pd.read_csv(args.input, dtype=str)

    if args.target:
        mask_target = df["title"].str.contains(args.target, case=False, na=False)
        df = df[mask_target].copy()

    mask_fold = df["title"].str.contains("FOLD Change", na=False)
    mask_22g = df["title"].str.contains("22G level", na=False)

    df_fold = df[mask_fold].copy().reset_index(drop=True)
    df_22g = df[mask_22g].copy().reset_index(drop=True)

    # mut_list = ["del", "mis"]
    mut_list = ["del", "subs"]

    def split_by_mut_and_without(df_group):
        result = {}
        for mut in mut_list:
            df_mut = df_group[df_group["mut_name"] == mut].copy()
            mask_w0 = df_mut["title"].str.contains("without 0", case=False, na=False)
            df_w0 = df_mut[mask_w0].copy().reset_index(drop=True)
            df_no_w0 = df_mut[~mask_w0].copy().reset_index(drop=True)
            result[(mut, "without0")] = df_w0
            result[(mut, "with0")] = df_no_w0
        return result

    split_fold = split_by_mut_and_without(df_fold)
    split_22g = split_by_mut_and_without(df_22g)

    # def process_and_save_fold(df_subset, mut, wtag, rna, prefix, output_dir):
    #     df_copy = df_subset.copy()
    #     if not df_copy.empty:
    #         df_copy.insert(0, "Parenthesis", df_copy["title"].map(extract_parenthesis))
    #     else:
    #         df_copy = pd.DataFrame(columns=["Parenthesis"])

    #     suffix = "_m" if rna == "miRNA" else "_c"
    #     cols_suffix = [c for c in df_subset.columns if c.endswith(suffix)]
    #     final_cols = ["Parenthesis"] + cols_suffix
    #     final_df = df_copy.reindex(columns=final_cols)

    #     fname = "{}_fold_{}_{}_{}.csv".format(prefix, mut, wtag, rna)
    #     path = os.path.join(output_dir, fname)
    #     final_df.to_csv(path, index=False, encoding="utf-8-sig")
    #     print("已輸出：{}".format(path))

    def process_and_save_fold(df_subset, mut, wtag, rna, prefix, output_dir, suffix_rev=False):
        df_copy = df_subset.copy()
        if not df_copy.empty:
            df_copy.insert(0, "Parenthesis", df_copy["title"].map(extract_parenthesis))
        else:
            df_copy = pd.DataFrame(columns=["Parenthesis"])
        if suffix_rev:
            suffix = "_c" if rna == "miRNA" else "_m"
            comparison_text = "(mut > common)" if suffix == "_c" else "(mut < common)"
        else:
            suffix = "_m" if rna == "miRNA" else "_c"
            comparison_text = "(mut > common)" if suffix == "_m" else "(mut < common)"

        # 找出所有符合後綴的欄位
        cols_suffix = [c for c in df_subset.columns if c.endswith(suffix)]
        final_cols = ["Parenthesis"] + cols_suffix
        final_df = df_copy.reindex(columns=final_cols)

        # 重新命名符合的欄位
        rename_dict = {col: col.replace(suffix, " {}".format(comparison_text)) for col in cols_suffix}
        final_df = final_df.rename(columns=rename_dict)

        fname = "{}_fold_{}_{}_{}.csv".format(prefix, mut, wtag, rna)
        path = os.path.join(output_dir, fname)
        final_df.to_csv(path, index=False, encoding="utf-8-sig")
        print("已輸出：{}".format(path))


    # def process_and_save_22g(df_subset, mut, wtag, prefix, output_dir):
    #     suffix = "_m"
    #     cols_suffix = [c for c in df_subset.columns if c.endswith(suffix)]
    #     final_cols = cols_suffix
    #     final_df = df_subset.reindex(columns=final_cols)

    #     fname = "{}_22G_{}_{}.csv".format(prefix, mut, wtag)
    #     path = os.path.join(output_dir, fname)
    #     final_df.to_csv(path, index=False, encoding="utf-8-sig")
    #     print("已輸出：{}".format(path))

    # def process_and_save_22g(df_subset, mut, wtag, prefix, output_dir):
    #     df_copy = df_subset.copy()
    #     if not df_copy.empty:
    #         df_copy.insert(0, "Parenthesis", df_copy["title"].map(extract_parenthesis))
    #     else:
    #         df_copy = pd.DataFrame(columns=["Parenthesis"])

    #     suffix = "_m"
    #     cols_suffix = [c for c in df_subset.columns if c.endswith(suffix)]
    #     final_cols = ["Parenthesis"] + cols_suffix
    #     final_df = df_copy.reindex(columns=final_cols)

    #     fname = "{}_22G_{}_{}.csv".format(prefix, mut, wtag)
    #     path = os.path.join(output_dir, fname)
    #     final_df.to_csv(path, index=False, encoding="utf-8-sig")
    #     print("已輸出：{}".format(path))

    def process_and_save_22g(df_subset, mut, wtag, prefix, output_dir, suffix_rev=False):
        df_copy = df_subset.copy()
        if not df_copy.empty:
            df_copy.insert(0, "Parenthesis", df_copy["title"].map(extract_parenthesis))
        else:
            df_copy = pd.DataFrame(columns=["Parenthesis"])

        if suffix_rev:
            suffix = "_c"
        else:
            suffix = "_m"
        comparison_text = "(mut > common)"

        cols_suffix = [c for c in df_subset.columns if c.endswith(suffix)]
        final_cols = ["Parenthesis"] + cols_suffix
        final_df = df_copy.reindex(columns=final_cols)

        rename_dict = {col: col.replace(suffix, " {}".format(comparison_text)) for col in cols_suffix}
        final_df = final_df.rename(columns=rename_dict)

        fname = "{}_22G_{}_{}.csv".format(prefix, mut, wtag)
        path = os.path.join(output_dir, fname)
        final_df.to_csv(path, index=False, encoding="utf-8-sig")
        print("已輸出：{}".format(path))


    for key, subdf in split_fold.items():
        process_and_save_fold(subdf, key[0], key[1], args.rna, args.prefix, args.output_dir, args.suff_rev)

    for key, subdf in split_22g.items():
        process_and_save_22g(subdf, key[0], key[1], args.prefix, args.output_dir, args.suff_rev)

if __name__ == "__main__":
    main()
