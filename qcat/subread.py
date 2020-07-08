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


def scanRead(read, partitions, quality = None):

    res = detector.scan(read, None, None, None)

    if res['barcode'] :

        partitions.append(res['barcode'].name)

        bcEnd = res['adapter_end'] + 1
        bcStart = res['adapter_end'] - len(res['barcode'].sequence) + 1

        left = read[:bcStart]
        right = read[bcEnd:]

        leftQual = None
        rightQual = None

        if quality:
            leftQual = quality[:bcStart]
            rightQual = quality[bcEnd:]

        scanRead(left, partitions, leftQual)
        scanRead(right, partitions, rightQual)


def bipartRead(read, bcOrder):

    res = detector.scan(read, None, None, None)

    if res['barcode'] :

        bcEnd = res['adapter_end'] + 1
        bcStart = res['adapter_end'] - len(res['barcode'].sequence) + 1

        left = read[:bcStart]
        right = read[bcEnd:]

        bcOrder[bcStart] = dict()
        bcOrder[bcStart]['barcode'] = re.sub("_.*","",res['barcode'].name)
        bcOrder[bcStart]['end'] = res['adapter_end'] + 1

        bipartRead(left, bcOrder)
        bipartRead(right, bcOrder)


def runLigationQcat():
    statsFile = open(args.output,"w")

    linkerMargin = False
    linkerSize = 12
    junkFilterSize = 40

    output_files = {}

    fastq = True

    bcDetectionStats = dict()
    subreadStats = dict()

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

def runStats():
    statsFile = open(args.output,"w")

    readCount = 0
    linkerStats = dict()

    with pysam.FastxFile(args.fastq) as handle:
        for read in handle:

            res = detector.scan(read.sequence, None, None, None)

            if res['barcode'] :
                if not res['barcode'].name in linkerStats:
                    linkerStats[res['barcode'].name] = 0
                linkerStats[res['barcode'].name] += 1

            readCount += 1

    print("Linker\tCount",file=statsFile)
    print("Total\t" + str(readCount),file=statsFile)
    for stat in linkerStats:
        print(stat + "\t" + str(linkerStats[stat]),file=statsFile)
    statsFile.close()

def runScan():
    statsFile = open(args.output,"w")

    readCount = 0
    linkerStats = dict()
    linkerStats["EMPTY"] = 0
    linkerStats["TOTAL"] = 0

    with pysam.FastxFile(args.fastq) as handle:
        for read in handle:

            barcodes = []

            scanRead(read.sequence, barcodes)

            if len(barcodes) > 0:
                for barcode in barcodes:
                    if not barcode in linkerStats:
                        linkerStats[barcode] = 0
                    linkerStats[barcode] += 1
            else :
                linkerStats["EMPTY"] += 1

            linkerStats["TOTAL"] += 1

            if linkerStats["TOTAL"] % 100 == 0:
                progressString = "Processed " + str(linkerStats["TOTAL"]) + " reads..."
                print(f'{progressString}\r', end="", file=sys.stderr)

    print(file=sys.stderr)
    print("Linker\tCount",file=statsFile)
    for stat in linkerStats:
        print(stat + "\t" + str(linkerStats[stat]),file=statsFile)
    statsFile.close()

def splitByBarcodeChunk(read, bcOrder, quality = None) :
    curBc = ""

    startList = list()
    readList = list()

    for curStart in sorted(bcOrder.keys()):
        if (curBc == ""):
            curBc = bcOrder[curStart]['barcode']
        if curBc != bcOrder[curStart]['barcode'] :
            start = curStart
            end = bcOrder[curStart]['end']
            startList.append((start, end))
            curBc = bcOrder[curStart]['barcode']

    leftQual = None
    rightQual = None

    if len(startList) > 0:
        prevEnd = 0
        for start, end in startList:
            left = read[prevEnd:start]
            if quality:
                leftQual = quality[prevEnd:start]
            prevEnd = start
            readList.append((left, leftQual))
        right = read[prevEnd:]
        rightQual = quality[prevEnd:]
        readList.append((right, rightQual))
    else :
        readList.append((read, quality))

    return readList
    # if start >= 0:
    #     left = read[:start]
    #     right = read[start:]
    #     if quality:
    #         leftQual = quality[:start]
    #         rightQual = quality[start:]
    #     return [(left, leftQual), (right, rightQual)]
    # else :
    #     return [(read, quality)]


def runBipart():

    statsFile = open(args.output,"w")

    linkerMargin = False
    linkerSize = 12
    junkFilterSize = 40

    output_files = {}

    fastq = True

    bcDetectionStats = dict()
    subreadStats = dict()

    with pysam.FastxFile(args.fastq) as handle:
        for read in handle:
            bcOrder = dict()
            bipartRead(read.sequence, bcOrder)
            biparts = splitByBarcodeChunk(read.sequence, bcOrder,  read.quality)

            for bipart in biparts:
                partitions = {}
                subReads = []
                partitionRead(bipart[0], partitions, subReads, linkerMargin, linkerSize, 0, "C", bipart[1])

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
parser.add_argument("--scan",
                    dest="SCAN",
                    action='store_true',
                    help="Scan full reads")
parser.add_argument("--bipart",
                    dest="BIPART",
                    action='store_true',
                    help="Bipart reads before ligation")
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
                    help="File qcaget stats will be written to ")

args = parser.parse_args()


detector = factory(mode="simple",
                   min_quality=60,
                   kit=args.barcodes)

sys.setrecursionlimit(10000)

if args.STATS:
    runStats()
elif args.SCAN:
    runScan()
elif args.BIPART:
    runBipart()
else :
    runLigationQcat()
