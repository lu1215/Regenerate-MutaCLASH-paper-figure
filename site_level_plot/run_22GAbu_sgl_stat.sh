python 22GAbuGroupbysglstatsig.py  --csv Input/PRG-1_rep1_with_sgl_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/PRG-1_abu_score_sgl_stat/ 
python 22GAbuGroupbysglstatsig.py --csv Input/PRG-1_rep1_with_sgl_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/PRG-1_abu_readcount_sgl_stat --boundaries 10000/4/1/0
# python 22GAbuGroupbysglstatsig.py --csv Input/ALG-1_rep3_with_sgl_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode score --out_dir Output/ALG-1_abu_score_sgl_stat/ --boundaries 200/140/100/60 --RNA miRNA 
# python 22GAbuGroupbysglstatsig.py --csv Input/ALG-1_rep3_with_sgl_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode readcount --out_dir Output/ALG-1_abu_readcount_sgl_stat --boundaries 10000/4/1/0 --RNA miRNA

python 22GAbuGroupbysglstatsig.py --csv Input/PRG-1_rep1_with_sgl_stat.csv  --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/PRG-1_abu_sgl_stat_combined 
# python 22GAbuGroupbysglstatsig.py --csv Input/ALG-1_rep3_with_sgl_stat.csv --target Input/add_two_HCLee.RNAseq.master.xlsx --id_map Input/mRNA_WS275_IDtoName.csv --mode combined --out_dir Output/ALG-1_abu_sgl_stat_combined --RNA miRNA --score_boundaries 200/140/100/60 --rc_boundaries 10000/4/1/0

python process_stat_csv.py -i Output/PRG-1_abu_sgl_stat_combined/all_stats_summary.csv -p PRG-1_wago -r piRNA -o Output/PRG-1_abu_sgl_stat_combined -t WAGO-1 --suff_rev
# python process_stat_csv.py -i Output/ALG-1_abu_sgl_stat_combined/all_stats_summary.csv -p ALG-1 -r miRNA -o Output/ALG-1_abu_sgl_stat_combined --suff_rev

