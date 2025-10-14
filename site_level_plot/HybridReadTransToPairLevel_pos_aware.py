# -*- coding: utf-8 -*-
# Filename: HybridReadTransToPairLevel_pos_aware.py
# Python: 3.5+

import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt


class SampleDataProcessor(object):
    """
    依據 sample_name 動態決定分組欄位：
      - 一般情況：以 ['transcript_name', 'regulator_name'] 分組
      - 若 sample_name == 'ALG-1'：
          以 ['transcript_name', 'regulator_name', 'mir_init_pos', 'mir_end_pos'] 分組
          並在輸出結果中包含 'mir_init_pos', 'mir_end_pos'

    功能：
      1. 讀取 CSV
      2. 依分組鍵彙總 read_count 與 row_count
      3. 統計 mutation/deletion/mismatch 讀數
      4. 加工 M 與 D（延用既有邏輯以 "[]" 判斷）
      5. 計算加權平均（22G_rc_WT/MUT；rpkm 欄位沿用原先取群組第一列的值）
      6. 取各類 score 的最大值（mir_score / targeting_score / RNAup_score）
      7. 輸出 DataFrame 或 CSV，必要時加上 mir_init_pos / mir_end_pos
    """

    def __init__(self, csv_path, sample_name):
        self.csv_path = csv_path
        self.sample_name = sample_name
        # 若 sample_name 正好是 'ALG-1'，啟用位置敏感分組
        # self.use_positions = (str(sample_name) == "ALG-1")
        self.df = None
        self.result_df = None

    def _group_keys(self):
        base = ['transcript_name', 'regulator_name']
        if str(self.sample_name) == "ALG-1" or str(self.sample_name) == "Homo":
            return base + ['mir_target_pos']
        elif str(self.sample_name) == "PRG-1":
            return base + ['pirscan_target_endpos']
        else: # Fly (RNAup)
            return base + ['RNAup_target_pos']
        return base

    def process(self):
        df = pd.read_csv(self.csv_path)
        self.df = df

        # 依照 "[]" 字串判斷是否有突變/缺失/錯配（沿用原本邏輯）
        mask_mut = (df['D'] != "[]") | (df['M'] != "[]")
        mask_del = df['D'] != "[]"
        mask_mis = df['M'] != "[]"

        group_keys = self._group_keys()

        # ---- 分組彙總 read_count ----
        total = df.groupby(group_keys)['read_count'].sum().reset_index().rename(
            columns={'read_count': 'sum_readcount'}
        )
        mutated = df[mask_mut].groupby(group_keys)['read_count'].sum().reset_index().rename(
            columns={'read_count': 'mutation_readcount'}
        )
        deleted = df[mask_del].groupby(group_keys)['read_count'].sum().reset_index().rename(
            columns={'read_count': 'deletion_readcount'}
        )
        mismatched = df[mask_mis].groupby(group_keys)['read_count'].sum().reset_index().rename(
            columns={'read_count': 'mismatch_readcount'}
        )

        # ---- 合併統計表 ----
        summary = total.merge(mutated, how='left', on=group_keys) \
                       .merge(deleted,  how='left', on=group_keys) \
                       .merge(mismatched, how='left', on=group_keys) \
                       .fillna(0)

        # 供後續以 tuple key 查詢
        summary.set_index(group_keys, inplace=True)

        rc_cols = ['22G_rc_WT', '22G_rc_MUT']
        rpkm_cols = ['wt.r1.rpkm', 'wt.r2.rpkm', 'wt.r3.rpkm',
                     'alg1.r1.rpkm', 'alg1.r2.rpkm', 'alg1.r3.rpkm']
        fly_rpm_cols = [
            'KGFP_1','KGFP_2','KGFP_3','KGFP_4','KGFP_5','KGFP_0',
            'PIWI_1','PIWI_2','PIWI_3','PIWI_4','PIWI_5','PIWI_0'
        ]

        fold_change_cols = ['foldchange_without0', 'foldchange_with0']

        def agg_group(g):
            # g.name 會是分組鍵的 tuple（或單值），統一轉成 tuple
            if isinstance(g.name, tuple):
                key_tuple = g.name
            else:
                key_tuple = (g.name,)

            # 用與 summary.set_index 相同的順序組成 key
            # 注意：group_keys 的順序需與 set_index 時一致
            # 這裡假設 pandas 保持 groupby keys 的順序（與定義一致）
            sum_rc = summary.at[key_tuple, 'sum_readcount']
            cnt = len(g)

            res = {}
            # 基本欄位
            res['row_count'] = cnt
            res['sum_readcount'] = sum_rc
            res['mutation_readcount'] = summary.at[key_tuple, 'mutation_readcount']
            res['deletion_readcount'] = summary.at[key_tuple, 'deletion_readcount']
            res['mismatch_readcount'] = summary.at[key_tuple, 'mismatch_readcount']

            # 將分組鍵值寫回輸出（含 ALG-1 的 mir_init_pos / mir_end_pos）
            # 以 group_keys 與 key_tuple 一一對應
            for idx, k in enumerate(group_keys):
                res[k] = key_tuple[idx]

            # 加權平均 22G 欄位
            for col in rc_cols:
                if col in g.columns:
                    out_name = 'weight_avg_' + col.lower().replace('22g_rc_', '22g_')
                    if sum_rc:
                        res[out_name] = (g[col] * g['read_count']).sum() / float(sum_rc)
                    else:
                        res[out_name] = np.nan

            # rpkm 欄位（沿用你原本做法：取群組第一列）
            for col in rpkm_cols:
                if col in g.columns:
                    out_name = 'weight_' + col
                    res[out_name] = g[col].iloc[0]

            # for col in fly_rpm_cols:
            #     if col in g.columns:
            #         out_name = 'weight_avg_' + col.lower()
            #         if sum_rc:
            #             res[out_name] = (g[col] * g['read_count']).sum() / float(sum_rc)
            #         else:
            #             res[out_name] = np.nan

            # 各 score 取最大值
            if 'mir_score' in g.columns:
                res['max_mir_score'] = g['mir_score'].max()
            if 'targeting_score' in g.columns:
                res['max_pir_score'] = g['targeting_score'].max()
            if 'RNAup_score' in g.columns:
                res['max_up_score'] = g['RNAup_score'].min()
            # save fly foldchange cols
            if self.sample_name == 'Fly':
                res['foldchange_without0'] = g['foldchange_without0'].iloc[0]
                res['foldchange_with0'] = g['foldchange_with0'].iloc[0]


            return pd.Series(res)

        grouped = df.groupby(group_keys)
        result = grouped.apply(agg_group).reset_index(drop=True)

        # 欄位順序：分組鍵 -> 統計欄位 -> 其他衍生欄位
        base_cols = list(group_keys) + [
            'row_count', 'sum_readcount',
            'mutation_readcount', 'deletion_readcount', 'mismatch_readcount'
        ]
        cols = list(base_cols)
        for col in result.columns:
            if col not in cols:
                cols.append(col)
        result = result[cols]

        self.result_df = result
        return result

    def save(self, output_folder=None):
        if self.result_df is None:
            raise RuntimeError("Please run process() before saving.")
        if output_folder is None:
            output_path = "./{}_site_level_result.csv".format(self.sample_name)
        else:
            output_path = "{}/{}_site_level_result.csv".format(output_folder, self.sample_name)
        self.result_df.to_csv(output_path, index=False)
        print("Saved to:", output_path)

    def analyze_score_distribution(self, output_folder=None):
        if self.result_df is None:
            raise RuntimeError("Please run process() before analyzing.")

        df = self.result_df.copy()

        # 選擇要分析的 score 欄位（沿用你的邏輯）
        if  (self.sample_name == 'ALG-1' or self.sample_name == 'Homo') and ('max_mir_score' in df.columns):
            score_col = 'max_mir_score'
            label = 'miRanda'
        elif self.sample_name == 'PRG-1' and 'max_pir_score' in df.columns:
            score_col = 'max_pir_score'
            label = 'pirScan'
        elif 'max_up_score' in df.columns:
            score_col = 'max_up_score'
            label = 'RNAup'
        else:
            print("No score column found to analyze.")
            return

        total = len(df)
        stat = df[score_col].value_counts().sort_index().reset_index()
        stat.columns = ['score', 'count']
        stat['proportion'] = (stat['count'] / float(total)).round(3)
        stat['cumulative'] = stat['proportion'].cumsum().round(3)

        # 存統計 CSV（首行加總筆數註解）
        stat_lines = ['# total count = {}'.format(total)]
        stat_lines += stat.to_csv(index=False).splitlines()

        if output_folder is None:
            csv_path = '{}_score_distribution.csv'.format(self.sample_name)
        else:
            csv_path = os.path.join(output_folder, '{}_score_distribution.csv'.format(self.sample_name))
        with open(csv_path, 'w') as f:
            for line in stat_lines:
                f.write(line + '\n')
        print("Score distribution CSV saved to:", csv_path)

        # 畫分數分布（折線圖）
        plt.figure(figsize=(8, 5))
        # 若未來要相容你之前的色彩設定可保留，也可移除 color 以採用預設色
        plt.plot(stat['score'], stat['count'], linestyle='-', linewidth=2)
        plt.xlabel('{} Score'.format(label))
        plt.ylabel('Count')
        plt.title('{} Score Distribution'.format(label))
        plt.tight_layout()
        plt.grid(True)
        plt.show()

        if output_folder is None:
            plot_path = '{}_score_distribution.svg'.format(self.sample_name)
        else:
            plot_path = os.path.join(output_folder, '{}_score_distribution.svg'.format(self.sample_name))
        plt.savefig(plot_path)
        plt.close()
        print("Score distribution plot saved to:", plot_path)


def main():
    parser = argparse.ArgumentParser(description="Process RNA CSV data by transcript/regulator group, optionally position-aware for ALG-1.")
    parser.add_argument("--csv_file", required=True, help="Input CSV file path")
    parser.add_argument("--sample_name", required=True, help="Sample name for output naming")
    parser.add_argument("--output", help="Optional output folder path")

    args = parser.parse_args()

    processor = SampleDataProcessor(args.csv_file, args.sample_name)
    processor.process()
    processor.save(args.output)
    processor.analyze_score_distribution(args.output)


if __name__ == "__main__":
    main()
