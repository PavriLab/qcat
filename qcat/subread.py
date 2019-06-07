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


detector = factory(mode="simple",
                   min_quality=60,
                   kit="/groups/pavri/bioinfo/mihaela/TCSeq/nanopore/multifail/barcodes/barcodes_simple.fa")


read = "TCATGATATGCTTCGTTCCGAATTTACGTATTGCTTGGGCTAGATATCCTCTACGGGCGACTACAAGCAGGAGCATCGACACTCTGAAGCCAAGGCCGATGGCGATTCACGGAGCTCTGCAGGGCTAAGTCCTGCTCGAAGGAGGCGGGGACCACGGGCGGCTGCTGGTCAGGCAGGCGTCACTGATAGTAGGGAGTCCACGCGTACCCTATAGTCTGACGACTACAAACGGAATCGAGAGGATATCTGAACTCCAATTGGAGAAGTAGAACGAACGACTACAAACCTGGAGATCGACCTCTGAAGCCAAGGAGCCGATGGCGATTCCTGGGCGTCTGAGGCTAAGTCCCTCTCCCGGGGAGAGCAGACTCCGGAGCAGCTGCTGGTCAGCGAGCGTCACTGACTATAGGGGAGTCCACGCGTGCCCTATAGTCACAGACGACTACAAACGGAATCGATCTATACTCCAATTGGTAGTTTGTATGAGAGACGAACGACTACAAACAGGATCAAGCCTCTGACCAAGGCCGGTGGCGATTCCTGGGCATCTGCTGAGGCTAGAGTCCTGCTCATGAAGGAGGCGGGGACTCGGAGCAGCTGCTGGTCCGACAGAGCGTCACCACGCGTGCCCTATAATTACGAACGGCTACAAACAGAATCGACTCCTTACACAAACTACCCAGAACTACTGCACCAACGAACTGACTACAAAACGGAATCGACCTCTGAAGCCAAGGCCGATGGCGATTCCTGGGCGTCTGCAGGGCTAAGTCCCTGCTCGAAGAGCTGAGGACTCGGAGCAGCTCACTAATTCCGACGAGCGTCACCCACGCGTGCCCTATAGTCACGGGCGACTACAAACGGAATCGATATACGGCGTTCTGGCAGATTTGTAGAGACGAACTGACTACAAACGGGAATCGACCTCGGCAAGGCCGATGGCGATTCCTGGGCGTCTGCGGGCTAAGTCCCTGCTCGAAGGAGAGCGGGGACTCCGGGACAGCTGCTAGTCCGACGAGCGTCACCCACGCATTGCCCTAGGAAGAGAGATGCAGACGACTACAAGCGGAATCGACTCCTTACAAACTACCAGAACTTTACTGCATAACAGACGACTACAAGCAGGACGGCCTCTGAAGCCAAGGCCGATGGCGATTTCACAGGCGTCTGCAGAATCCCTGCTCGAAGGAAGCGGGGACTCGGAGCAGCTGCTAGTCCGACGAGCGTCACTGATAGTAGGAAGTAAAAAGAGTGCGTCTCTCCCCCAACCCGCACACACACACACACACACACACACACACACACACACACACAGAGCCCCGAGGCAGATGTAAACCGTAGAGCTAGCTGAGATGCCAAACCGTAGACCTCTGGCAGTGCTGATAAACCTAAGCTACCTCTGGATTGGCTGTGGAGTCTCCCATCCTGGGTGATACCCAGAAGGCATGCGTGCATGTCCCGCGTGCCCTATGGTCACACAGGAACTGACTACAAGCGAATCGATATGCACAGTGAGAAGTTCTGGGTTTGTGTAAGGGAGTCGATTCCGTTTGTGATCGTCTGTGACTATAGGGCACGCGTGGAGGCTATACAGCGAACACCCTGTATCTCAAACAAACAGTTGACCCCTTTGCATGCATGAGCCCCATTGTGGGTTCTCAGGAAGTATCGACTGGACATGGATTACAGGTCCCCCTTGCTTCTCATTCTGCGGGCCAAGTGTCTTGGAAGTTCGGCTCACACAGGCACTTCGGCCTGAAATGAGCTCAGGGTCCATGTCAGCTACGAACACAATGAAACAAGGAGCTCAGCTGTTTGTTTGAAGGCGGGAATGATACTTCTGTGCCCTATAGTACGACGACTACAAACGGAATCGACTCCTTACACAAACTCCCAGAACTTCTGCATATCGATTCCGTTTGTAGTCGTCTCTTCTAGGGCACGCGTGGGTATCCTCATCTTTAAACTCCAGTTAACTTCTAAATACTTCCAGCAGAGCCCTGGGTGTTGGGTAGAGGTGCACACGTGTTCTCATGTCTCAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGTTGGGGGGGAGGCATGCGCTCTTTCCCTACTATGTTGACGCTGGCGTCGGACAACAGCTGCTCGGAATTTCAGCGACGAGGACAGCCCCTGCGCCAGGAATCGCCTCGGCCTTGGCTTCCAGAGGTCGATTCCGTTTTGTAGTCTAAAAATCTGTTATGCAGTAGTTCTGGTAGTTTGTGTAAGGGAGTCGATTCCGTTTGTAGTCGTCCTGTGACTATGAGGGCACGCGTGGGTGACGCTCGTCGGACTAGCAGCTGCTCCAGTCCCCCGCCTCCTTCGAGCAGGGACTTAGCCCTGCGACGCCCGAATCGCCATCGGCCTTGAAGGTCGATTCCGTTTGTAGTCGTCTGTTACTCCTTAAAGATCTGCCAGAACTACTGCATATCGATTCCGTTTGTAGTCGTCTGTGACTATAAGGCTGAGCCTGGCGCAGCGCGACACAATCACAGACCACTCCCCTTCAGCATGTACTTCCAAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGTTGGGGGAAGGCACTGCGCTCTTTTGCTCCTACCTCGGTGACGCTCGTCGGACTAGCAGCTGCTTGAGTCCGCCTCGCGAGCAGGGACTTGGCCCTGCGACGCCCAGGAATCGCCGTGGCCTTGGCTTCCAGAGGTCGATTCCGTTTGTAATGTTTATGCAGTAGTTCTGGAGTTTGTGTAAGGAGTCGATTCCGTTTGTAGTCGTCTGTGACTATAGGGCACGCGTGGTGACAGCTCGTCGGACTAGCAGCTGCTCCGAGTCCCCGCCTCCTTCGAGCGAGGACTTAGCCCTGCAGACCCTTTAGAATCGCATCGTGGCTTCAGAGATTCGATTCCGTTTGTAGTCGTCTGTCTCCTTACACAAACTACCAATTGGAGAAATGAACTCGATTCCGTTTGTAGTCGTCGGTGACTATGGGCGCGTGAATTTGACGCTCGTCGGACTAGCAGCTGCTCGAGTCCCCGCCTCCTTCGAGCAGGGACTTAGCCCTGCAGACTGCCAGGAATCGCCATCGGCTGATTCCGTTTGTAGTCGTCTGTTCTACTTTCTCAATTGGAGTTCAGATGCTCTCTTCGATTCCGTTTGTAGTCGTCTGTGACTATAGGGCACGCGGTGACACCGTCGGACTAGCAGCTGCTCCAGTCCCCGCACCTCCGCTCGAGCAGGGACAGCCCTGCGAACGCCCAGGAATCGCCATCGGCCTTGGCCCAAGAGGTCGATTCCGTTTGTAGTCGTCTGTAGAGGATGCTGAACTAGCAATACGTGTATG"

linkerMargin = False
linkerSize = 12
junkFilterSize = 40

output_files = {}

#res = detector.scan(read, None, None, None)

#partitions = {}
#subReads = []
#partitionRead(read, partitions, subReads, linkerMargin, linkerSize)

#detectedBc = sorted(partitions.items(), key=lambda kv: kv[1], reverse=True)[0][0]



#print(partitions)
#print(subReads)

#name = "test"
#comment = "test"
#quality = "&%&%&%"

# for subRead in subReads:
#
#     if fastq:
#         print("@" + name + " " + comment, subRead, "+", quality, sep="\n",
#               file=out_file)
#     else:
#         print(">" + name + " " + comment, subRead, sep="\n",
#               file=out_file)

fastq = True

with pysam.FastxFile('/groups/pavri/bioinfo/mihaela/TCSeq/nanopore/multifail/guppy/nanopore_multifail.fq.gz') as handle:
    for read in handle:
        partitions = {}
        subReads = []
        partitionRead(read.sequence, partitions, subReads, linkerMargin, linkerSize, 0, "C", read.quality)

        detectedBc = "none"

        if partitions:
            detectedBc = sorted(partitions.items(), key=lambda kv: kv[1], reverse=True)[0][0]
            detectedBc = re.sub("_.*","",detectedBc)

        out_file = get_output_file(output_files, '.', detectedBc, fastq)

        for subRead in subReads:

            if fastq:
                print("@" + read.name + " " + read.comment + " level=" + str(subRead.level) + " clade=" + subRead.clade, subRead.sequence, "+", subRead.quality, sep="\n",
                      file=out_file)
            else:
                print(">" + name + " " + comment, subRead, sep="\n",
                      file=out_file)
