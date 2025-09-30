#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import re
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from tqdm import tqdm
import matplotlib.ticker as ticker

# ---------- Helper Functions ----------

def get_identity_distribution_and_avg_mapq(bam_file):
    total_mapq = []
    total_nb_reads = 0

    sum_identity = 0
    sum_edit_distance = 0
    sum_read_length = 0
    count_identity = 0

    identity_values_hist = []
    edit_distances_hist = []
    read_lengths_hist = []

    if not os.path.exists(bam_file + ".bai"):
        pysam.index(bam_file)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.query_sequence is None:
                continue
            edit_distance = read.get_tag("NM")
            read_length = len(read.query_sequence)
            identity = (read_length - edit_distance) / read_length * 100

            identity_values_hist.append(identity)
            edit_distances_hist.append(edit_distance)
            read_lengths_hist.append(read_length)

            sum_identity += identity
            sum_edit_distance += edit_distance
            sum_read_length += read_length
            count_identity += 1

            total_mapq.append(read.mapping_quality)
            total_nb_reads += 1

    avg_mapq = np.mean(total_mapq) if total_mapq else 0
    avg_identity = sum_identity / count_identity if count_identity else 0

    return identity_values_hist, edit_distances_hist, read_lengths_hist, total_nb_reads, avg_mapq, avg_identity

def evenness_data(bam_file, N_tiles=500):
    refs = []
    positions = []
    depths = []

    for line in pysam.depth("-a", bam_file, split_lines=True):
        ref, pos, depth = line.strip().split("\t")
        refs.append(ref)
        positions.append(int(pos))
        depths.append(int(depth))

    df = pd.DataFrame({'Reference': refs, 'Position': positions, 'Depth': depths})
    mean_coverage = round(df['Depth'].mean(), 2)
    percent_bases_covered = round((df['Depth'] > 0).sum() / len(df['Depth']) * 100, 2)

    genome_length = len(df)
    tile_size = genome_length / N_tiles
    boc = np.zeros(genome_length, dtype=float)
    tile_covered_count = 0

    for i in range(N_tiles):
        start = int(i * tile_size)
        end = int(min((i + 1) * tile_size, genome_length))
        df_tile = df.iloc[start:end]
        if len(df_tile) == 0:
            continue
        frac_covered = (df_tile['Depth'] > 0).sum() / len(df_tile)
        boc[start:end] = frac_covered
        if frac_covered > 0:
            tile_covered_count += 1

    percent_tiles_covered = round(tile_covered_count / N_tiles * 100, 2)
    return df, boc, mean_coverage, percent_bases_covered, percent_tiles_covered

def reverse_complement(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(complement.get(b,b) for b in reversed(seq))

def calculate_dna_damage(bam_file, max_pos=31):
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    dna_damage = {"CpG": [0]*max_pos, "CT": [0]*max_pos, "other": [0]*max_pos}
    cpg_sites = [0]*max_pos
    c_sites = [0]*max_pos
    other_sites = [0]*max_pos

    for read in bamfile.fetch():
        if read.mapping_quality == 0 or read.query_sequence is None:
            continue
        seq = read.query_sequence
        if not read.has_tag("MD") or not read.has_tag("NM"):
            continue
        MD = read.get_tag("MD")
        reverse = read.is_reverse
        cigar = read.cigarstring
        if "N" in seq or any(x in cigar for x in "IDNSHP"):
            continue

        # rebuild reference sequence
        refseq, newread_seq = "", ""
        MDlist = re.findall(r'(\d+|\D+)', MD)
        counter = 0
        for m in MDlist:
            if m.isdigit():
                refseq += seq[counter:counter+int(m)]
                newread_seq += seq[counter:counter+int(m)]
                counter += int(m)
            elif '^' in m:
                continue
            else:
                refseq += m
                newread_seq += seq[counter]
                counter += len(m)

        if reverse:
            refseq = reverse_complement(refseq)
            newread_seq = reverse_complement(newread_seq)

        for i in range(min(max_pos, len(newread_seq)-1)):
            if "N" in refseq[i:i+2] or "N" in newread_seq[i:i+2]:
                continue

            # CpG and CT mismatches
            if refseq[i] == "C" and refseq[i+1] != "G":
                c_sites[i] += 1
                if newread_seq[i] == "T":
                    dna_damage["CT"][i] += 1
            elif refseq[i] == "C" and refseq[i+1] == "G":
                cpg_sites[i] += 1
                if newread_seq[i] == "T":
                    dna_damage["CpG"][i] += 1
            else:
                # "other" mismatches (exclude C/T & GA)
                if refseq[i] not in "CT" and newread_seq[i] not in "CT" and refseq[i]+newread_seq[i] != "GA":
                    other_sites[i] += 1
                    if refseq[i] != newread_seq[i]:
                        dna_damage["other"][i] += 1

    CT = [dna_damage["CT"][i]/c_sites[i] if c_sites[i] > 0 else 0 for i in range(max_pos)]
    CpG = [dna_damage["CpG"][i]/cpg_sites[i] if cpg_sites[i] > 0 else 0 for i in range(max_pos)]
    other = [dna_damage["other"][i]/other_sites[i] if other_sites[i] > 0 else 0 for i in range(max_pos)]
    return CT, CpG, other

# ---------- Plot a single BAM ----------

def process_bam_file(bam_file, pdf_pages, stats_row=None):
    # --- Compute read metrics ---
    identity_values, edit_distances, read_lengths, total_nb_reads, avg_mapq, avg_identity = get_identity_distribution_and_avg_mapq(bam_file)
    
    # --- Compute evenness metrics ---
    df_even, boc, mean_coverage, percent_bases_covered, percent_tiles_covered = evenness_data(bam_file)
    
    # --- Compute DNA damage ---
    CT, CpG, other = calculate_dna_damage(bam_file)
    
    fig = plt.figure(figsize=(22,20))
    gs = gridspec.GridSpec(3,2, height_ratios=[1,1,1.5])
    
    if stats_row is not None:
        plt.suptitle(f"{stats_row['reference']}", fontsize=22, weight='bold')
    
    # --- Edit Distance Histogram ---
    ax0 = plt.subplot(gs[0,0])
    ax0.hist(edit_distances, bins=50, color='#1f77b4', alpha=0.7, edgecolor='black')
    ax0.set_title("Edit Distance", fontsize=14, weight='bold')
    ax0.set_xlabel("Edit Distance (NM)", fontsize=12)
    ax0.set_ylabel("Count", fontsize=12)
    ax0.grid(True, linestyle='--', alpha=0.5)
    
    # --- Identity Histogram ---
    ax1 = plt.subplot(gs[0,1])
    ax1.hist(identity_values, bins=50, color='#2ca02c', alpha=0.7, edgecolor='black')
    ax1.set_title("% Identity", fontsize=14, weight='bold')
    ax1.set_xlabel("% Identity", fontsize=12)
    ax1.set_ylabel("Count", fontsize=12)
    ax1.set_xlim(85, 100)
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # --- Read Length Distribution (number of reads) ---
    ax2 = plt.subplot(gs[1,0])
    readlen_dict = [0]*301
    for rl in read_lengths:
        readlen_dict[min(300, rl)] += 1
    x_axis_readlen = [i for i in range(301) if readlen_dict[i] > 0]
    y_axis_readlen = [readlen_dict[i] for i in x_axis_readlen]  # absolute number of reads

    ax2.fill_between(x_axis_readlen, 0, y_axis_readlen, alpha=0.25, color='#ff7f0e')
    ax2.plot(x_axis_readlen, y_axis_readlen, color='#ff7f0e', linewidth=2, linestyle='-', marker='o', markersize=3)
    ax2.set_title("Read Length Distribution", fontsize=14, weight='bold')
    ax2.set_xlabel("Read length (bp)", fontsize=12)
    ax2.set_ylabel("Number of reads", fontsize=12)
    xmax = max(100, max(read_lengths))
    ax2.set_xlim(0, xmax)
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    # --- Evenness plot ---
    ax3 = plt.subplot(gs[1,1])
    ax3.plot(df_even['Position'], boc, color='black', linewidth=1)
    ax3.set_xlabel("Reference position", fontsize=12)
    ax3.set_ylabel("Fraction covered", fontsize=12)
    ax3.set_title(f"Evenness\nMean cov: {mean_coverage}, {percent_bases_covered}% genome, {percent_tiles_covered}% tiles", fontsize=14, weight='bold')
    ax3.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,_:'{:0.0e}'.format(x)))
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    # --- DNA damage plot ---
    ax_damage = plt.subplot(gs[2,0])
    x = list(range(31))
    ax_damage.plot(x, CT, label="C to T", color="#d62728", linewidth=3)
    ax_damage.plot(x, CpG, label="CpG to TpG", color="#1f77b4", linestyle='--', linewidth=2)
    ax_damage.plot(x, other, label="Other", color="#2ca02c", linestyle=':', linewidth=2)
    ax_damage.set_xlabel("Distance from read end (bp)", fontsize=12)
    ax_damage.set_ylabel("Mismatch frequency", fontsize=12)
    ax_damage.set_xticks(range(0,31,2))
    ax_damage.set_ylim(0,0.5)
    ax_damage.set_title("DNA damage by read position", fontsize=14, weight='bold')
    ax_damage.grid(True, linestyle='--', alpha=0.5)
    ax_damage.legend(loc="upper right", fontsize=12)
    
    # --- Metrics Table ---
    ax_table = plt.subplot(gs[2,1])
    ax_table.axis('off')
    if stats_row is not None:
        metrics = ['n_reads','mapping_quality','edit_distances','read_ani_mean','coverage_mean', 'breadth_exp_ratio','site_density','norm_entropy','norm_gini']
        values = [f"{stats_row[m]:.2f}" for m in metrics]
        table_data = list(zip([m.replace("_"," ") for m in metrics],values))
        table = ax_table.table(cellText=table_data, colLabels=['Metric','Value'], loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(16)
        table.scale(1,3)
        for key, cell in table.get_celld().items():
            if key[0]==0:
                cell.set_text_props(weight='bold', fontsize=16)
            cell.set_facecolor('#f0f0f0')

    plt.tight_layout(rect=[0,0.03,1,0.95])
    pdf_pages.savefig(fig)
    plt.close()

# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description="BAM plotting with table and DNA damage")
    parser.add_argument("--bam_file", type=str, required=True)
    parser.add_argument("--filterBAM_table", type=str, required=True)
    parser.add_argument("--out", type=str, required=True)
    parser.add_argument("--min_reads", type=int, default=100)
    parser.add_argument("--min_breadth", type=float, default=0.60)
    parser.add_argument("--min_norm_entropy", type=float, default=0.8)
    parser.add_argument("--min_norm_gini", type=float, default=0.3)
    args = parser.parse_args()

    tab = pd.read_csv(args.filterBAM_table, sep="\t")
    tab = tab[(tab['n_reads'] >= args.min_reads) &
              (tab['breadth_exp_ratio'] >= args.min_breadth) &
              (tab['norm_entropy'] >= args.min_norm_entropy) &
              (tab['norm_gini'] <= args.min_norm_gini)]

    with PdfPages(args.out) as pdf_pages:
        if tab.empty:
            fig = plt.figure(figsize=(8.5, 11))
            plt.axis("off")
            plt.text(0.5, 0.7, "No candidate bacteria", ha="center", va="center", fontsize=24, weight="bold")

            settings_text = (
                f"Current filtering settings:\n\n"
                f"Minimum reads: {args.min_reads}\n"
                f"Minimum breadth: {args.min_breadth}\n"
                f"Minimum normalized entropy: {args.min_norm_entropy}\n"
                f"Maximum normalized Gini: {args.min_norm_gini}"
            )
            plt.text(0.5, 0.4, settings_text, ha="center", va="center", fontsize=14)

            pdf_pages.savefig(fig)
            plt.close()
            return

        with pysam.AlignmentFile(args.bam_file,"rb") as bam_in:
            for idx,row in tqdm(tab.iterrows(),total=len(tab),desc="Processing references"):
                reference = row['reference']
                sub_bam_file = f"{os.path.splitext(args.bam_file)[0]}_{reference}.bam"

                with pysam.AlignmentFile(sub_bam_file,"wb",template=bam_in) as bam_out:
                    for read in bam_in.fetch(reference=reference):
                        bam_out.write(read)
                if not os.path.exists(sub_bam_file+".bai"):
                    pysam.index(sub_bam_file)
                process_bam_file(sub_bam_file,pdf_pages,stats_row=row)

if __name__=="__main__":
    main()

