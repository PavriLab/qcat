from __future__ import print_function

import logging
import time
import sys
import os
import six
import re
import pysam

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from argparse import ArgumentParser, RawDescriptionHelpFormatter, ArgumentTypeError

from qcat.scanner import factory
from qcat import __version__

class SubRead:

    def __init__(self, sequence=None, quality=None, level=None, clade=None):
        self.sequence = sequence
        self.quality = quality
        self.level = level
        self.clade = clade

    def __str__(self):
     return str(self.level) + "-" + self.clade + "|" + self.sequence + "|" + self.quality

def get_output_file(output_files, out_folder, barcode, fastq):

    if barcode not in output_files.keys():
        if fastq:
            output_files[barcode] = open(
                os.path.join(out_folder, barcode + ".fastq"), "w")
        else:
            output_files[barcode] = open(
                os.path.join(out_folder, barcode + ".fasta"), "w")
    return output_files[barcode]

def partitionRead(read, partitions, subReads, linkerMargin, linkerSize, level = 0, clade = "C", quality = None):

    res = detector.scan(read, None, None, None)

    if res['barcode'] :

        if not res['barcode'].name in partitions:
            partitions[res['barcode'].name] = 0
        if not re.sub("_.*","",res['barcode'].name) in partitions:
            partitions[re.sub("_.*","",res['barcode'].name)] = 0

        partitions[res['barcode'].name] += 1
        partitions[re.sub("_.*","",res['barcode'].name)] += 1

        bcEnd = res['adapter_end'] + 1
        bcStart = res['adapter_end'] - len(res['barcode'].sequence) + 1

        if linkerMargin:
            bcEnd += linkerSize
            bcStart -= linkerSize

        left = read[:bcStart]
        right = read[bcEnd:]

        leftQual = None
        rightQual = None

        if quality:
            leftQual = quality[:bcStart]
            rightQual = quality[bcEnd:]

        partitionRead(left, partitions, subReads, linkerMargin, linkerSize, level + 1, "L", leftQual)
        partitionRead(right, partitions, subReads, linkerMargin, linkerSize, level + 1, "R", rightQual)

    else :
        if len(read) > 0:
                subReads.append(SubRead(read, quality, level, clade))


usage = "Python command-line tool for splitting ligated Oxford Nanopore reads from FASTQ files"

parser = ArgumentParser(description=usage,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-V", '--version',
                    action='version',
                    version='%(prog)s ' + __version__)
parser.add_argument("--stats",
                    dest="STATS",
                    action='store_true',
                    help="Only run stats")
parser.add_argument("-f", "--fastq",
                    type=str,
                    dest="fastq",
                    required=True,
                    help="Barcoded read file")
parser.add_argument('-b', "--barcodes",
                    dest="barcodes",
                    type=str,
                    required=True,
                    help="qcaget will use barcodes from this file")
parser.add_argument('-o', "--output",
                    dest="output",
                    type=str,
                    required=True,
                    help="File qcaget stats will be writen to ")

args = parser.parse_args()


detector = factory(mode="simple",
                   min_quality=60,
                   #kit="/groups/pavri/bioinfo/mihaela/TCSeq/nanopore/multifail/barcodes/barcodes_simple.fa")
                   kit=args.barcodes)

statsFile = open(args.output,"w")

linkerMargin = False
linkerSize = 12
junkFilterSize = 40

output_files = {}

fastq = True

bcDetectionStats = dict()
subreadStats = dict()

#with pysam.FastxFile('/groups/pavri/bioinfo/mihaela/TCSeq/nanopore/ligation_2/guppy/nanopore_ligation_2.fq.gz') as handle:
with pysam.FastxFile(args.fastq) as handle:
    for read in handle:
        partitions = {}
        subReads = []
        partitionRead(read.sequence, partitions, subReads, linkerMargin, linkerSize, 0, "C", read.quality)

        detectedBc = "none"

        if partitions:
            detectedBc = sorted(partitions.items(), key=lambda kv: kv[1], reverse=True)[0][0]
            detectedBc = re.sub("_.*","",detectedBc)

        if detectedBc not in bcDetectionStats :
            bcDetectionStats[detectedBc] = 0
        bcDetectionStats[detectedBc] += 1

        out_file = get_output_file(output_files, '.', detectedBc, fastq)

        subReadCount = 0

        for subRead in subReads:

            if fastq:
                print("@" + read.name + " " + read.comment + " level=" + str(subRead.level) + " clade=" + subRead.clade, subRead.sequence, "+", subRead.quality, sep="\n",
                      file=out_file)
            else:
                print(">" + name + " " + comment, subRead, sep="\n",
                      file=out_file)

            subReadCount += 1

        if subReadCount not in subreadStats:
            subreadStats[subReadCount] = 0
        subreadStats[subReadCount] += 1

print("Barcode\tCount",file=statsFile)
for bc in bcDetectionStats:
    print(bc + "\t" + str(bcDetectionStats[bc]),file=statsFile)
print(file=statsFile)
print("Detected subreads\tCount",file=statsFile)
for count in sorted(subreadStats):
    print(str(count) + "\t" + str(subreadStats[count]),file=statsFile)

statsFile.close()
