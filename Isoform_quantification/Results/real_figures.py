import pandas as pd
import os
import generate_result_stat as result_process
import generate_stat as stats
from util import ExperimentFileProcessor
import matplotlib.pyplot as plt
from gen_stats import get_stats
import numpy as np


def plot_bar_graph(stats, models, fig_name, title):
    # Data preparation
    
    metrics = ['sc', 'pc', 'im']

    # Organize data into arrays for each metric
    # sc_values = [stats['kl_stats']['sc'], stats['sal_stats']['sc'], stats['e1_sstats']['sc'], stats['e4_sstats']['sc']]
    # pc_values = [stats['kl_stats']['pc'], stats['sal_stats']['pc'], stats['e1_sstats']['pc'], stats['e4_sstats']['pc']]
    # im_values = [stats['kl_stats']['im'], stats['sal_stats']['im'], stats['e1_sstats']['im'], stats['e4_sstats']['im']]

    sc_values = [stats[k]['sc'] for k in stats.keys()]
    pc_values = [stats[k]['pc'] for k in stats.keys()]
    im_values = [stats[k]['im'] for k in stats.keys()]

    # Set width of bars and positions of the bars
    bar_width = 0.25
    r1 = np.arange(len(models))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]

    # Create the plot
    plt.figure(figsize=(10, 6))

    plt.rcParams.update({'font.size': 14})  # Sets a base font size
    plt.rcParams['axes.titlesize'] = 16  # Sets the title font size
    plt.rcParams['axes.labelsize'] = 14  # Sets the axes labels font size


    # Create bars
    plt.bar(r1, sc_values, width=bar_width, label='SC', color='skyblue')
    plt.bar(r2, pc_values, width=bar_width, label='PC', color='lightgreen')
    plt.bar(r3, im_values, width=bar_width, label='IM', color='salmon')

    # Customize the plot
    plt.xlabel('Models')
    plt.ylabel('Correlation of Abundance')
    plt.title(title)
    plt.xticks([r + bar_width for r in range(len(models))], models)
    plt.legend()

    # Ensure y-axis starts from 0 and ends at 1.2
    plt.ylim(0, 1.2)

    # Add grid for better readability
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.savefig(fig_name)  # Save as a file
    print('saved')

def main():
    # Day 5 is best performing

    # Short Read Files
    kl_r1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_51/abundance.tsv' # kallisto
    kl_r2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_52/abundance.tsv'
    sal_r1 = '/gpfs/commons/home/sraghavendra/SOTA/salmon/salmon_quant51/quant.sf' # salmon
    sal_r2 = '/gpfs/commons/home/sraghavendra/SOTA/salmon/salmon_quant52/quant.sf'
    e1_sr1 = '../exprmnt_2025_01_21__17_21_14/files/output_files/output_PacIllu_VIGD_token_470401_sample1_file1_ds100num1aln51short_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_06_40_34.tsv' # exp 1
    e1_sr2 = '../exprmnt_2025_01_21__17_21_14/files/output_files/output_PacIllu_VIGD_token_5426678_sample1_file1_ds100num1aln52short_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_11_33_28.tsv'
    e4_sr1 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_24785461_sample2_file1_ds10num1aln01long_file2_ds100num1aln51short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_14_10_05.tsv' # exp 4
    e4_sr2 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_1196578_sample2_file1_ds10num1aln02long_file2_ds100num1aln52short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_14_38_40.tsv'

    # Long Read Files
    lrkl_r1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_51/abundance.tsv' #l r-kallisto
    lrkl_r2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_52/abundance.tsv'
    nc_r1 = '/gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_51.tsv' # NanoCount
    nc_r2 = '/gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_52.tsv'
    oar_r1 = '/gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1/.quant' # Oarfish
    oar_r2 = '/gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2/.quant'
    # mp_r1 = '' # MPAQT
    # mp_r2 = ''
    e1_lr1 = '../exprmnt_2025_01_21__22_53_07/files/output_files/output_PacIllu_VIGD_token_10807412_sample1_file1_ds10num1aln51long_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_23_27_48.tsv' #exp 1
    e1_lr2 = '../exprmnt_2025_01_21__22_53_07/files/output_files/output_PacIllu_VIGD_token_32391888_sample1_file1_ds10num1aln52long_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_23_44_03.tsv'
    e4_lr1 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_19932137_sample1_file1_ds10num1aln51long_file2_ds100num1aln01short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_08_40_46.tsv' #exp 4
    e4_lr2 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_16541151_sample1_file1_ds10num1aln52long_file2_ds100num1aln02short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_09_00_04.tsv'

    # Simulation LR Files
    gt1 = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/sample1_gt.tsv'
    gt2 = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/sample2_gt.tsv'
    kl1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_sim1/abundance.tsv'
    kl2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_sim2/abundance.tsv'

    # # SR stats
    # kl_stats = get_stats(kl_r1, kl_r2, 'target_id', 'target_id', 'tpm', 'tpm')
    # # print(kl_stats)

    # sal_stats = get_stats(sal_r1, sal_r2, 'Name', 'Name', 'TPM', 'TPM')
    # # print(sal_stats)

    # e1_sstats = get_stats(e1_sr1, e1_sr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # # print(e1_sstats)

    # e4_sstats = get_stats(e4_sr1, e4_sr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # # print(e4_sstats)

    # # LR stats
    # lrkl_stats = get_stats(lrkl_r1, lrkl_r2, 'target_id', 'target_id', 'tpm', 'tpm')
    # # print(lrkl_stats)

    # nc_stats = get_stats(nc_r1, nc_r2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # # print(nc_stats)

    # oar_stats = get_stats(oar_r1, oar_r2, 'tname', 'tname', 'num_reads', 'num_reads')
    # # print(oar_stats)

    # e1_lstats = get_stats(e1_lr1, e1_lr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # # print(e1_lstats)

    # e4_lstats = get_stats(e4_lr1, e4_lr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # # print(e4_lstats)

    # Sim LR Stats
    kl1_stats = get_stats(kl2, gt2, 'target_id', 'transcript_name', 'tpm', 'tpm')
    print(kl1_stats)

    # # SR Bar Plot
    # models = ['JOLI\nSingle Sample', 'JOLI\nMulti-Sample', 'kallisto', 'salmon']
    # stats = {'e1_sstats': e1_sstats, 'e4_sstats': e4_sstats, 'kl_stats': kl_stats, 'sal_stats': sal_stats}
    # fig_name = 'sr_real_res.png'
    # plot_bar_graph(stats, models, fig_name, 'Illumina SR')

    # # LR Bar Plot
    # models = ['JOLI\nSingle Sample', 'JOLI\nMulti-Sample', 'lr-kallisto', 'NanoCount', 'Oarfish']
    # stats = {'e1_lstats': e1_lstats, 'e4_lstats': e4_lstats, 'lrkl_stats': lrkl_stats, 'nc_stats': nc_stats, 'oar_stats': oar_stats}
    # fig_name = 'lr_real_res.png'
    # plot_bar_graph(stats, models, fig_name, 'PacBio LR')

    #Sim Results Bar Plot

if __name__ == "__main__":
    main()