#!/usr/bin/env python

# Copyright (C) 2022.  All rights reserved.
#
# See the file LICENSE for redistribution information.
#
# Author: Tobias Neumann
# Email: tobias.neumann.at<at>gmail.com

# Date located in: -
from __future__ import print_function
import sys, os, re

import subprocess

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from os.path import basename

from Bio import SeqIO
from Bio.Seq import Seq

import pysam

usage = "Extract potential translocations from Nanopore data"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--bam", type=str, required=True, dest="bamFile", help="Bam file")
parser.add_argument("-r", "--roi", type=str, required=False, default = "chr15:61983341-61992361", dest="roi", help="Region of interest (chr15:61983341-61992361)")
parser.add_argument("-s", "--split", action='store_true', dest="split", help="Split translocations by 100bp bait threshold.")

args = parser.parse_args()

###################
# Gather bait reads
###################

#bait = dict()

fields = args.roi.split(':')
roichr = fields[0]
roistart, roiend = fields[1].split('-')
roistart = int(roistart)
roiend = int(roiend)

offTarget = dict()

bamFile = pysam.AlignmentFile(args.bamFile, "rb")

readCount = 0

if args.split:
    fupper = open('translocations_greater100bp.bed', 'w')
    flower = open('translocations_smaller100bp.bed', 'w')

for read in bamFile.fetch(reference=roichr, start=roistart, end=roiend):

    if read.has_tag("SA"):
        #print(read)
        #print(bamFile.getrname(read.reference_id), end="\t")
        #print(read.reference_start, end="\t")
        #print(read.cigarstring)

        suppAlns = read.get_tag("SA").split(";")

        for aln in suppAlns:
            if aln != "":

                chr, start, strand, cigar, mapQ, NM = aln.split(",")

                if args.split:
                    if read.query_alignment_length >= 100:
                        print(chr + "\t" + start + "\t" + str(int(start) + 1) + "\t" + read.query_name + "\t0\t+", file = fupper)
                    else :
                        print(chr + "\t" + start + "\t" + str(int(start) + 1) + "\t" + read.query_name + "\t0\t+", file = flower)
                else :
                    print(chr + "\t" + start + "\t" + str(int(start) + 1) + "\t" + read.query_name + "\t0\t+")
