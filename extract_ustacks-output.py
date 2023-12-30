#!/usr/bin/env python3

# Raphael 21/10/14; for Stacks 2.59; scRAD_GonyPop

# Extract the following parameters and results from ustacks output file and write to *.csv:
# read_file = Input file name (forward read file)
# sample_id = Unique sample ID
# m_par = Min depth of coverage to create a stack (m)
# M_par = Max distance allowed between stacks (M)
# N_par = Max distance allowed to align secondary reads (N)
# loaded_reads = Number of reads loaded into ustacks
# prim_stacks = Number of stacks formed by primary reads
# prim_reads = Number of primary reads
# rep_reads = Number of reads after repeat removal
# dist_reads = Number of reads after assembling stacks
# merge_reads = Number of reads after merging secondary stacks
# mean_cov = Mean final coverage
# stdev_cov = Standard deviation final coverage
# max_cov = Maximum final coverage
# final_reads = Final number of reads utilized in ustacks
#   (after merging secondary stacks, allowing for gaps if gapped assembly enabled)


import glob
import csv
import re

for file in glob.glob("*ustacks*.out"):  # ustacks output file path
    ustacks_sum = []
    filename = file.split(sep=".")[0]
    with open(file, "rt") as infile:
        file_check = False
        for line in infile:
            line = line.strip()  # remove new-line character from end of line
            if "ustacks parameters selected:" in line:
                file_check = True  # check if file is from ustacks
                ustacks_out = {}
            elif "Input file:" in line:
                read_file = line.split(sep="/")[-1][:-1]
            elif "Sample ID:" in line:
                sample_id = line.split(sep=" ")[-1]
            elif "Min depth of coverage to create a stack (m):" in line:
                m_par = line.split(sep=" ")[-1]
            elif "Max distance allowed between stacks (M):" in line:
                M_par = line.split(sep=" ")[-1]
            elif "Max distance allowed to align secondary reads:" in line:
                N_par = line.split(sep=" ")[-1]
            elif "Loaded" in line:
                loaded_reads = line.split(sep=" ")[1]
            elif "primary reads" in line:
                prim_stacks = line.split(sep=" ")[0]
                prim_reads = line.split(sep=" ")[3]
            elif "Coverage after repeat removal:" in line:
                rep_reads = line.split(sep="=")[-1].split(sep="(")[0]
            elif "Coverage after assembling stacks:" in line:
                dist_reads = line.split(sep="=")[-1].split(sep="(")[0]
            elif "Coverage after merging secondary stacks:" in line:
                merge_reads = line.split(sep="=")[-1].split(sep="(")[0]
            elif "stacks into" in line and "blacklisted" not in line:
                final_stacks = line.split(sep=" ")[-2]
                print(final_stacks)
            elif "Final coverage" in line:
                mean_cov = line.split(sep=";")[0].split(sep="=")[-1]
                stdev_cov = line.split(sep=";")[1].split(sep="=")[-1]
                max_cov = line.split(sep=";")[2].split(sep="=")[-1]
                final_reads = line.split(sep="=")[-1].split(sep="(")[0]
            elif "ustacks is done" in line:
                ustacks_out = {"read_file": read_file, "sample_id": sample_id, "m_par": m_par, "M_par": M_par,
                               "N_par": N_par, "loaded_reads": loaded_reads, "prim_stacks": prim_stacks,
                               "prim_reads": prim_reads, "rep_reads": rep_reads, "dist_reads": dist_reads,
                               "merge_reads": merge_reads, "final_reads": final_reads, "final_stacks": final_stacks,
                               "mean_cov": mean_cov, "stdev_cov": stdev_cov, "max_cov": max_cov}
                ustacks_sum.append(ustacks_out)
    if file_check is True:
        with open(".".join((filename, "csv")), "wt", newline="") as outfile:
            keys = ustacks_sum[0].keys()
            dict_writer = csv.DictWriter(outfile, keys)
            dict_writer.writeheader()
            dict_writer.writerows(ustacks_sum)
        print("%s: done" % file)
        print("Number of samples: %i" % len(ustacks_sum))
    else:
        print("%s: check file format!" % file)
