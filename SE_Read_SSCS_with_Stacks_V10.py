#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 1 11:32 2017 by Daniel E. Deatherage

consensus reads from single end amplicon sequencing

Suggested invocation:
SE_Read_SSCS_with_Stacks_V7.py -f <fastq_name.fastq> -p sample_name

steps:
    1. read fastq file
    2. identify reads with correct expected sequences on 5` and 3`
        command line allows for indels in barcode region, but requires perfect match on expected sequences
            perfect matches required in order for biological functionality
    3. output fasta file with the following format:
        >UniqueNumber
        Barcode region
    4. run stacks with following considerations
        1. allow 1 MI to start a stack // -m 1
        2. allow X mismatchs between stacks // -M X (command line option)
        3. first consideration eliminates the possible existence of secondary reads, therefore disable their mapping // -N 0
        4. happlotypes do not exist within MI, therefore disable // -H
        5. this assumed to be run on stampede, use maximal processors. currently set to 16 // -p 16
        6. do not eliminate reads which are highly represented // --keep_high_cov
        7. add required unique integer ID // -i 1 
    5. read in .tags.tsv file
        1. primary sequence names contain BC
        2. write consensus BC out to .fa file
    6. generate histogram of number of occurrences of each BC from bc.tags.tsv file
    7. write output stats to log file

Version 1 Adapted from sscs_with_stacks_v3.py
Version 2 removes consensus read MI steps and just looks at straight amplicon sequencing histograms
Version 3 looks for 5` and 3` matching sequences, and only creates histogram of sequences with perfect matches
Version 4 clean up for publication
Version 5 update to ustacks V1.48 (from 1.13)
Version 6 update to include stutter start of some number of bp
    New options:
        -rl --read_length: length of reads in fastq file
        -ps --padding_sequence: sequence 5' of 5' constant region
Version 7 update rework to simpler match of smaller region of interest
Version 8 update to look for correct read index in 3' primer
    New options:
        "-is", "--index_sequence": expected read index sequence
Version 9 updated to include lev distance histogram file
Version 10 added default options used in Al'Khafaji, Deatherage and Brock "Control of lineage-specific gene expression by functionalized gRNA barcodes"  ACS Synthetic Biology

@author: ded
"""

import re
import argparse
import os
import subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(description="Generate histogram of prevalance of all consensus BC detected")
parser.add_argument("-cr5", "--constant_region_5prime", default="AACACCG", help="expected 5' constant region to search for. Will be removed before anlyzing barcode region")
parser.add_argument("-cr3", "--constant_region_3prime", default="GTTTTA", help="expected 3' constant region to search for. Will be removed before anlyzing barcode region")
parser.add_argument("-i", "--indel_length", default=0, type=int, help='length of indels allowed in barcode region.')
parser.add_argument("-e", "--expected", default="NSNWNSNWNSNWNSNWNSN", help="expected degenerate sequence of barcode region. Used to determine presence of barcode")
parser.add_argument("-t", "--enrichment_target", default="GACATGGATCGCTAGAACCG", help="specific sequence enriched for")
parser.add_argument("-m", "--mutations", default=1, type=str, help="number of mutations allowed between barcode regions to group together")
parser.add_argument("-is", "--index_sequence", type=str, help="expected adapter index sequence to identify near 3' end of read that was used to assign read to sample", default='GCGGAC', required=True)
parser.add_argument("-f", "--fastq", help="input fastq file", required=True)
parser.add_argument("-p", "--prefix", help="prefix for output files", required=True)
parser.add_argument("--testing", action='store_true', default=False, help="break fastq read in after 50k reads to increase speed")
args = parser.parse_args()

main_dict =  {"raw_reads": 0, "barcodeless": 0, "incorrect_assignment": 0,  "dual_constant": 0, "Unique BC Count": 0 }
levenshtein_distance_dict = defaultdict(int)
assert re.match("^[A-Za-z0-9_]+$", args.prefix), "\nNon-alphanumeric characters found:\n%s\nPlease restrict to letters, numbers, or underscore characters in prefix name" % [x for x in args.prefix if not re.match("[A-Za-z0-9_]", x)]
assert args.fastq.endswith(".fastq"), "\nSequence file specified (%s) does not appear to be a fastq file. fastq files should end with '.fastq'" % args.fastq

def histogram_generation():
    """Generate histogram tsv file based on ustacks output"""
    reads = []
    hist_dict = {}
    total_reads = 0.0
    with open(args.prefix + ".tags.tsv") as bc_tags_in:
        for line in bc_tags_in:
            line = line.rstrip().split("\t")
            if re.match("^#", line[0]):
                continue  # ustacks V 1.48 has header line at top
            if re.match("model", line[6]):
                continue  # not doing anything with model lines
            elif re.match("consensus", line[6]):
                if len(reads) > 0:  # consesnsus sequence is first thing encountered, first time through will do nothing, this is expected, and hence the call after completion of the loop
                    assert consensus not in hist_dict, "duplicated barcode found: %s" % consensus
                    hist_dict[consensus] = len(reads)
                    total_reads += len(reads)

                reads = []  # reset sequences to empty lists
                consensus = line[9]  # consensus sequence built from all raw sequences
                main_dict["Unique BC Count"] += 1
            elif re.match("primary", line[6]):
                reads.append(line[9])
                if args.testing and main_dict["raw_reads"] % 200000 == 0:
                    break  # for testing subset rather than full read list
            elif re.match("secondary", line[6]):
                assert False, "Secondary read alignment found, ustacks, not behaving as intended.\n%s" % line
            else:
                assert False, "Unknown sequence type identified, don't know how to handle this.\n%s\n%s" % (line[6], line)

        #final assignment
        assert consensus not in hist_dict, "duplicated barcode found: %s" % consensus
        hist_dict[consensus] = len(reads) # must run on last consensus sequence, will obviously have reads assigned
        total_reads += len(reads)
        assert total_reads == main_dict["dual_constant"],  "different number of total reads counted in tags.tsv (%s) and dual constant region (%s)" % (total_reads, main_dict["dual_constant"])

        #print to histogram file
        print>> histo_out, "\t".join(map(str, ["con_read", "count", "frequency", "lev distance from: " + args.enrichment_target ]))
        for con_read in hist_dict:
            print>>histo_out, "\t".join(map(str, [con_read, hist_dict[con_read], hist_dict[con_read] / total_reads, levenshtein(args.enrichment_target, con_read)]))


def fastq_read_in(fastq):
    """read in fastq file while checking for correct formatting, and write fasta file containing only reads matching expectations"""
    line_count = 0
    with open(fastq, "r") as fastq, open(args.prefix + ".log.txt", "w", 0) as log:
        for line1 in fastq:
            line1 = line1.rstrip()
            line_count += 1
            if line_count % 4 == 1:  # line = header
                header1 = line1
                assert len(header1.split(" ")) == 2, "Header line does not agree with expectation that only a single space in header line.\nHeader:\n%s\nAt line:\n%s" % (header1, line_count)  # this may be inaccurate for different Illumina outputs
            elif line_count % 4 == 2:  # line = read
                read1 = line1
                assert re.match('^[ACTGN]*$', read1) , "Non-sequence characters found on the following sequence:\n%s\nAt line:\n%s" % (line1, line_count)
            elif line_count % 4 == 3:  # line = optional sequence descriptor. assume line to use 'standard' "+" notation rather than actual descriptor.
                assert line1 == "+", "Line does not display expected '+' sign on line3. Instead has:\nline3: %s\nat line: %s" % (line1, line_count)
            elif line_count % 4 == 0:  # line = quality, last line of individual read
                main_dict["raw_reads"] += 1

                bc = re.findall(search_string, read1)
                if bc:  # if BC not found, empty list is false
                    assert len(bc) == 1, "Multiple search string matches on single read: %s" % read1

                    index_search = re.findall(args.index_sequence, read1.split(bc[0])[-1])
                    assert len(bc) <=1, "Index sequence (%s) found multiple times 3' of barcode target for read: %s\nVerify index sequence and if correct, strengthen index search string." % (args.index_sequence, read1)
                    if len(index_search) == 0:
                        main_dict["incorrect_assignment"] += 1
                    else:  # assert above limits options to 0 or 1
                        main_dict["dual_constant"] += 1
                        print>> fasta_out, "\n".join(map(str, [">" + str(main_dict["raw_reads"]), bc[0]]))
                else:  # BC not found
                    main_dict["barcodeless"] +=1

            if line_count % 200000 == 0:
                print>>log, line_count / 4, "reads processed"
                if args.testing:
                    break  # for testing subset of reads rather than full read list
        print>>log, "Fastq read-in complete. %d total reads processed"  % (line_count / 4)

def levenshtein(s1, s2):
    """taken from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python"""
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1  # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    levenshtein_distance_dict[int(previous_row[-1])] += 1
    return previous_row[-1]


# Generate bc search string based on command line options
search_string = args.constant_region_5prime + "([ACTGN]{%d,%d})" % (len(args.expected) - args.indel_length, len(args.expected) + args.indel_length) + args.constant_region_3prime
with open(args.prefix + ".log.txt", "a") as log:
    print>>log, "fastq file %s will be evaluated for reads containing the bc string: %s" % (args.fastq, search_string)

# Read in fastq file, prepare fasta file for MI consensus read generation
assert not os.path.isfile(args.prefix + ".fasta"), "fasta file already found"  # NOT a starting point as file exists even if script didn't finish correctly, or was run in testing mode
with open(args.prefix + ".fasta", "w") as fasta_out:
    fastq_read_in(args.fastq)

# Run stacks to generate *.tags.tsv file containing reads with same barcode
stacks_command = ["ustacks", "-t", "fasta", "-f", str(args.prefix) + ".fasta", "-m", "1", "-M", str(args.mutations), "-N", "0", "-H", "-p", "16", "--keep_high_cov", "-i", "1"]
external_stacks_call = subprocess.Popen(stacks_command, stdout=open(args.prefix + ".log.txt", "a"), stderr=open(args.prefix + ".log.txt", "a"))
external_stacks_call.wait()

# Histogram generation
assert not os.path.isfile(args.prefix + ".histogram.tsv"), "Histogram file already exists, please remove or rename existing file: %s.histogram.tsv" % args.prefix
with open(args.prefix + ".histogram.tsv", "w") as histo_out:
    histogram_generation()

# Levenshtein distance histogram generation
assert not os.path.isfile(args.prefix + ".levenshtein_distance.histogram.tsv"), "levenshtein_distance histogram file already exists, please remove or rename existing file: %s.levenshtein_distance.histogram.tsv" % args.prefix
with open(args.prefix + ".levenshtein_distance.histogram.tsv", "w") as lev_out:
    for _ in sorted(levenshtein_distance_dict):
        print>>lev_out, "\t".join(map(str, [_, levenshtein_distance_dict[_]]))

# Print final data
with open(args.prefix + ".log.txt", "a") as log:
    for entry in ["raw_reads", "barcodeless", "incorrect_assignment", "dual_constant", "Unique BC Count"]:
        print>>log, "\t".join(map(str, [entry, main_dict[entry]]))
    print>>log, "Successfully complete"
