#!/usr/bin/env python
from __future__ import print_function  # support python3-style print
import sys
import os
import re
import argparse

# 2022 Alan McCulloch
# MIT License
# Copyright (c) 2022 Alan McCulloch
# Version: 1.0
# Email: alan.mcculloch@agresearch.co.nz

# ref /dataset/gseq_processing/active/bin/batch_qiime_prism/add_sample_name.py
# which does a similar thing - i.e. sample name is encoded in the
# filename, and you want to merge that with the contents of the file.

def splice_sample_names(args):
    """
    for each file, parse the filename to get sampleID and locusID,
    then parse the file contents and splice these in and list, print
    updated summary table to stdout.
    """
    for path in args["input_files"]:
        # parse the filename to get sampleID and locusID parts
        # e.g. from HERB01-A07-J11-65-14_S27_EPSPS.summary.txt we want to be
        # able to extract HERB01-A07-J11-65-14  EPSPS try each regular
        # expression (not case-sensitive) and use the first match
        # (usually there will be only one regexp and usually the default)
        sample_name = None
        locus = None

        for regexp in args["regexps"]:
            match = re.search(regexp, os.path.basename(path), re.IGNORECASE)
            if match is not None:
                if len(match.groups()) == 2:
                    (sample_name, locus) = match.groups()

        if sample_name is None or locus is None:
            raise Exception(
                "Sorry could not parse sampleID and locusID from %s using any of : %s" % (
                    path, str(args["regexps"])))

        # now open file, parse contents and splice in sample name and locus
        # file contents are like :
        # more /dataset/gseq_processing/active/bin/gtseq_prism/unit_test/HERB01-A07-J11-65-14_S27_ACCase1999.summary.txt
        # 0.0000  ATTCCCATGAGCGGTCTGTTCCTCGTGCTGGGCAAGTCTGGTTTCCAGATTCTGCTACCAAGACAGCGCAGGCAATGTTGGACTTCAA;size=2072      *       *       *       *       *       *       *       *       0       0   00       0       0       *       N
        # 0.0000  ATTCCCATGAGCGGTCTGTTTCTCGTGCTGGGCAAGTCTGGTTTCCAGATTCTGCTACCAAGACAGCGCAGGCAATGTTGGACTTCAA;size=42        *       *       *       *       *       *       *       *       0       0   00       0       0       *       N
        with open(path, "r") as savs:
            record_count = 0
            field_count = 0
            for record in savs:
                fields = re.split("\t", record.strip())
                if record_count == 0:
                    field_count = len(fields)
                record_count += 1
                # sanity check file - records should all be the same length
                if len(fields) != field_count:
                    raise Exception(
                        "oops : in %s, first rec had %d fields, but record %d has %d fields" % (
                            path, field_count, len(fields)))
                allele_and_size = re.split(";", fields[1])
                if len(allele_and_size) != 2:
                    raise Exception(
                        "oops : malformed sampleID field (%s) in record %d" %
                        fields[1], record_count)
                (allele, size) = allele_and_size
                # output record , e.g. like CCCGATTGAGAAGGATGCCAAGGAGGAAGTAAAGCTCTTCTTGGGCAACGCTGGAACTGCAATGCGGCCATTGACGGCAGCTGTAGTAGCTGCTGGTGGA  size=1125  HERB01-A07-J11-65-14  EPSPS
                print(allele, size, sample_name, locus, sep="\t")

def get_options():
    description = """
    for each file, parse the filename to get sampleID and locusID, 
    then parse the file contents and splice these in and list, print
    updated summary table to stdout.
    """
    long_description = """
    examples :
    python splice_sample_names.py /dataset/gseq_processing/active/bin/gtseq_prism/unit_test/HERB*.txt > updated.summary.txt
    """
    parser = argparse.ArgumentParser(description=description,
                                     epilog=long_description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_files',
                        type=str,
                        nargs="*",
                        help='input summary files from vASV workflow; vsearch_uchime3 summary output files.')
    parser.add_argument('-t', '--task', dest='task', required=False,
                        default="splice_sample_names", type=str,
                        choices=["splice_sample_names"],
                        help='task to execute; default: splice_sample_names')
    parser.add_argument('-x', '--regexps',
                        dest='regexps',
                        type=str,
                        default='^([^_]+)_.*_(\w+)\.summary\.txt$,^([^_]+)_.*_([\w\-]+)\.summary\.txt$',
                        help="""
                        default: '^([^_]+)_.*_(\w+)\.summary\.txt$,^([^_]+)_.*_([\w\-]+)\.summary\.txt$'; 
                        Comma-separated list of regexps to parse sample name; First default for 
                        e.g. HERB01-A07-J11-65-14_S27_EPSPS.summary.txt --> HERB01-A07-J11-65-14  EPSPS; 
                        Second default for loci with hyphens in name 
                        e.g. 'HERB01-A07-J11-65-14_S27_ACCase2027-2041.summary.txt'.""")
    parser.add_argument('-O', '--output_dir',
                        dest='output_dir',
                        type=str,
                        default=None,
                        help='output folder (not currently used - output written to stdout).')
    parser.add_argument('-G', '--genotyper',
                        dest='genotyper',
                        type=str,
                        default='vASV',
                        choices=['vASV'],
                        help='genotyper (not currently used).')
    args = vars(parser.parse_args())

    # check files exist
    for path in args["input_files"]:
        if not os.path.isfile(path) and not os.path.islink(path):
            raise Exception(
                "%s not found, or is not a file or link,  giving up" % path)

    # clean up regular expressions - remove leading and trailing quote characters that may have been inserted
    # (as an alternative to escaping shell metachars , when invoking this on the commandline - e.g. we may call it
    # like  -x "'-([^-_]+)_.*.fastq.gz'"
    args["regexps"] = re.split(",", args["regexps"])

    for i in range(0, len(args["regexps"])):
        args["regexps"][i] = re.sub("^[\"\']+", "", re.sub("[\"\']+$", "", args["regexps"][i]))
    return args

def main():
    opts = get_options()
    if opts["task"] == "splice_sample_names":
        splice_sample_names(opts)
    else:
        raise Exception(
            "Oops... --task %(splice_sample_names)s is not supported yet !" % opts)
    return 0

if __name__ == "__main__":
    main()
