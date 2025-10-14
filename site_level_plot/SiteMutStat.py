import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
from scipy.stats import pearsonr
import matplotlib.ticker as ticker


class MutaCLASHAnalyzer:
    """
    讀取「已整理好的 aggregated CSV」，建立繪圖所需的資料視圖，
    並產生散點圖與分佈檔。
    """
    def __init__(self, csv_path, experiment_name, output_dir=None, image_format='png'):
        """
        :param csv_path: 路徑：已整理好的 aggregated CSV（直接包含 *_readcount 與 *_nor_readcount 欄位）
        :param experiment_name: 實驗名稱（用於輸出檔名/圖標題）
        :param output_dir: 輸出資料夾（預設為 csv 所在資料夾）
        """
        self.csv_path = csv_path
        self.experiment_name = experiment_name
        self.image_format = image_format.lower()
        # 直接把輸入 CSV 當作已聚合好的表
        self.df = pd.read_csv(csv_path)

        # 決定輸出資料夾
        if output_dir:
            self.output_dir = output_dir
        else:
            self.output_dir = os.path.dirname(os.path.abspath(csv_path))
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.agg_df = None
        self.grouped_dfs = {}

    def process_data(self):
        """
        不再做 groupby 聚合。
        直接把輸入 CSV 視為已聚合的結果，並建好 grouped_dfs 供繪圖使用。
        若缺少預期欄位，補 0 以避免繪圖失敗。
        """
        df = self.df.copy()

        # 預期會用到的欄位（若缺失就補 0）
        needed_cols = [
            'sum_readcount', 'read_count', 'nor_readcount',
            'mutation_readcount', 'deletion_readcount', 'mismatch_readcount',
            'mutation_nor_readcount', 'deletion_nor_readcount', 'mismatch_nor_readcount'
        ]
        for col in needed_cols:
            if col not in df.columns:
                # 盡量不終止流程：缺的欄位補 0
                df[col] = 0

        self.agg_df = df.reset_index(drop=True)
        self.grouped_dfs = {
            'mutation': self.agg_df[['sum_readcount', 'mutation_readcount', 'nor_readcount', 'mutation_nor_readcount']].copy(),
            'deletion': self.agg_df[['sum_readcount', 'deletion_readcount', 'nor_readcount', 'deletion_nor_readcount']].copy(),
            'mismatch': self.agg_df[['sum_readcount', 'mismatch_readcount', 'nor_readcount', 'mismatch_nor_readcount']].copy()
        }
        return self.agg_df

    # 🔻 原本用來輸出 _aggregated.csv 的方法已移除
    # def save_aggregated_csv(...): pass

    def generate_plots(self):
        """
        產生四張散點圖（使用已聚合好的欄位）：
          1. sum_readcount vs. mutation_readcount（linear）
          2. sum_readcount vs. mutation_readcount（log-log）
          3. nor_readcount vs. mutation_nor_readcount（linear）
          4. nor_readcount vs. mutation_nor_readcount（log-log）
        """
        if self.agg_df is None:
            raise ValueError("Data not processed yet. Call process_data() first.")

        df = self.agg_df
        plots = [
            ('sum_readcount', 'mutation_readcount', False),
            ('sum_readcount', 'mutation_readcount', True),
            # ('nor_readcount', 'mutation_nor_readcount', False),
            # ('nor_readcount', 'mutation_nor_readcount', True)
        ]

        epsilon = 1  # 避免 log(0)

        for x_col, y_col, log_scale in plots:
            plot_df = df.copy()

            if log_scale:
                plot_df[x_col] = plot_df[x_col] + epsilon
                plot_df[y_col] = plot_df[y_col] + epsilon

            fig, ax = plt.subplots()
            ax.scatter(plot_df[x_col], plot_df[y_col])
            ax.set_xlabel(x_col)
            ax.set_ylabel(y_col)
            title = "{}: {} vs. {}".format(self.experiment_name, y_col, x_col)
            if log_scale:
                ax.set_xscale('log')
                ax.set_yscale('log')
                title += ' (log-log)'
            ax.set_title(title)

            fname = "{}_{}_vs_{}".format(self.experiment_name, x_col, y_col)
            if log_scale:
                fname += '_loglog'
            fname += '.png'

            out_path = os.path.join(self.output_dir, fname)
            fig.savefig(out_path, dpi=300, bbox_inches='tight')
            plt.close(fig)

    def generate_distribution_files(self):
        """
        針對 sum_readcount / nor_readcount / mutation_readcount / deletion_readcount / mismatch_readcount
        產生分佈 CSV 與兩張散點圖（count 與 percent）。
        """
        if self.agg_df is None:
            raise ValueError("Data not processed yet. Call process_data() first.")

        for column in ['sum_readcount', 'nor_readcount', "mutation_readcount", "deletion_readcount", "mismatch_readcount"]:
            if column not in self.agg_df.columns:
                # 若缺少欄位則略過（已在 process_data 補 0，理論上不會進來）
                continue

            dist = self.agg_df[column].value_counts().sort_index().reset_index()
            dist.columns = [column, 'count_of_row']
            total = dist['count_of_row'].sum()
            dist['percent'] = (dist['count_of_row'] / total * 100).round(3)
            dist['accumulated_percent'] = dist['percent'].cumsum().round(3)

            header_line = "# count of data = {}".format(total)
            output_file = os.path.join(self.output_dir, "{}_{}_distribution.csv".format(self.experiment_name, column))

            with open(output_file, 'w') as f:
                f.write(header_line + "\n")
                dist.to_csv(f, index=False)

            # count_of_row 散點圖
            fig1, ax1 = plt.subplots()
            ax1.scatter(dist[column], dist['count_of_row'])
            ax1.set_xlabel(column)
            ax1.set_ylabel('count of row')
            ax1.set_title('{} distribution of {} (count)'.format(self.experiment_name, column))
            plot_path1 = os.path.join(self.output_dir, '{}_{}_distribution_count_plot.{}'.format(self.experiment_name, column, self.image_format))
            fig1.savefig(plot_path1, dpi=300, bbox_inches='tight')
            plt.close(fig1)

            # percent 散點圖
            fig2, ax2 = plt.subplots()
            ax2.scatter(dist[column], dist['percent'])
            ax2.set_xlabel(column)
            ax2.set_ylabel('percent')
            ax2.set_title('{} distribution of {} (percent)'.format(self.experiment_name, column))
            plot_path2 = os.path.join(self.output_dir, '{}_{}_distribution_percent_plot.{}'.format(self.experiment_name, column, self.image_format))
            fig2.savefig(plot_path2, dpi=300, bbox_inches='tight')
            plt.close(fig2)

    def _plot_group_comparison(self, df, x_col, y_col, title_suffix, log_scale=False):
        if log_scale:
            df[y_col] = df[y_col] + 0.1

        r_val, _ = pearsonr(df[x_col], df[y_col]) if len(df) > 1 else (np.nan, None)

        fig, ax = plt.subplots()
        ax.scatter(df[x_col], df[y_col])
        ax.set_xlabel(x_col)
        ax.set_ylabel('{} + 0.1'.format(y_col) if log_scale else y_col)

        title = '{} {} vs {}{}'.format(self.experiment_name, y_col, x_col, ' (log-log)' if log_scale else '')
        title += ' (Pearson r = {:.3f})'.format(r_val) if not np.isnan(r_val) else ' (r = N/A)'
        ax.set_title(title)

        if log_scale:
            ax.set_xscale('log')
            ax.set_yscale('log')

        fname = '{}_{}_vs_{}{}_{}.{}'.format(
            self.experiment_name, x_col, y_col,
            '_loglog' if log_scale else '', title_suffix, self.image_format
        )
        fig.savefig(os.path.join(self.output_dir, fname), dpi=300, bbox_inches='tight')
        plt.close(fig)

    def generate_grouped_plots(self):
        if self.agg_df is None:
            raise ValueError("Data not processed yet. Call process_data() first.")

        for group, df in self.grouped_dfs.items():
            if '{}_readcount'.format(group) not in df.columns:
                continue

            self._plot_group_comparison(df, 'sum_readcount', '{}_readcount'.format(group), group, log_scale=False)
            self._plot_group_comparison(df, 'sum_readcount', '{}_readcount'.format(group), group, log_scale=True)
            # self._plot_group_comparison(df, 'nor_readcount', '{}_nor_readcount'.format(group), group, log_scale=False)
            # self._plot_group_comparison(df, 'nor_readcount', '{}_nor_readcount'.format(group), group, log_scale=True)

    def run(self):
        """
        新流程：僅「載入/準備」→ 分佈檔 → 分組散點圖
        不再輸出 _aggregated.csv
        """
        self.process_data()
        # self.generate_plots()  # 如需也可開啟
        self.generate_distribution_files()
        self.generate_grouped_plots()
        return self.output_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='使用已整理好的 aggregated CSV 繪製統計圖與分佈。'
    )
    parser.add_argument('--csv_path', help='已整理好的 aggregated CSV 路徑')
    parser.add_argument('--experiment_name', help='實驗名稱')
    parser.add_argument('--output_dir', help='輸出資料夾', default=None)
    parser.add_argument('--image_format', help='輸出圖檔格式 (png 或 svg)', default='png')

    args = parser.parse_args()

    analyzer = MutaCLASHAnalyzer(args.csv_path, args.experiment_name, args.output_dir, args.image_format)
    analyzer.run()

# 使用方式範例（已聚合好的 CSV）
# python PairingMutationStat.py --csv_path Aggregated/PRG-1_aggregated.csv --experiment_name PRG-1 --output_dir Output/PRG-1_Pairing_analysis
# python PairingMutationStat.py --csv_path Aggregated/ALG-1_aggregated.csv --experiment_name ALG-1 --output_dir Output/ALG-1_Pairing_analysis
