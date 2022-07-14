subread [![DOI](https://zenodo.org/badge/190543385.svg)](https://zenodo.org/badge/latestdoi/190543385)
======================================================================================================

subread is an extension of the original [qcat](https://github.com/nanoporetech/qcat) Python command-line tool for demultiplexing Oxford Nanopore reads from FASTQ files. We create this separate branch for demultiplexing concatenated subreads separated by barcode-linker pairs mainly utilizing qcat's ***BarcodeScanner*** class to detect ONT simple barcodes - for a full list check [this file](https://raw.githubusercontent.com/nanoporetech/qcat/master/qcat/resources/kits/simple_standard.yml) provided by ONT.

In addition, it provides a simple utility script ```call_lam_htgts_translocations.py``` to extract translocations produced by the LAM-HTGTS assay from a bam file and outputs all detected translocations in a bed file.

Requirements
------------
* Linux or MacOS

Quick start
-----------
To demultiplex a multiplexed and subread-containing datasets, use the ```subread.py``` script contained in the ```qcat``` subfolder:
```bash
$ python qcat/subread.py -f <fastq_file> -b <barcodes_fasta_file> -o <output stats file>
```
After subread finished, please check the stats file to verify that a barcode was assigned to most of the reads.

Installation 
------------
To install qcaget manually, please make sure you have python3 and git available, and install as follows:
```bash
$ git clone https://github.com/t-neumann/qcat.git
$ cd qcat
$ git checkout subread
$ python3 setup.py install
```

What is the barcode file format?
--------------------------------

The barcode file format is a simple fasta file with both forward and reverse barcodes listed as demonstrated below:

```
>bc1_fwd
AAGAAAGTTGTCGGTGTCTTTGTG
>bc1_rev
CACAAAGACACCGACAACTTTCTT
>bc2_fwd
TCGATTCCGTTTGTAGTCGTCTGT
>bc2_rev
ACAGACGACTACAAACGGAATCGA
>bc3_fwd
GAGTCTTGTGTCCCAGTTACCAGG
>bc3_rev
CCTGGTAACTGGGACACAAGACTC
```


How do I call LAM-HTGTS translocations from a bam file?
-------------------------------------------------------

To extract LAM-HTGTS translocations from an aligned bam file, please use the provided ```call_lam_htgts_translocations.py``` utility script as described below:
```bash
$ python lam_htgts_util/call_lam_htgts_translocations.py -i <bam_file> -r <bait region in chr:start-end format> | sort -k1,1 -k2,2n > <translocations_bed_file>
```
