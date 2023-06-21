# C3POa

[![Github release](https://img.shields.io/github/tag/christopher-vollmers/C3POa.svg?label=Version)](https://github.com/christopher-vollmers/C3POa/tags)
[![Published in PNAS](https://img.shields.io/badge/Published%20in-PNAS-blue.svg)](https://doi.org/10.1073/pnas.1806447115)
[![GPL license](https://img.shields.io/badge/License-GPL-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

C3POa (**C**oncatemeric **C**onsensus **C**aller with **P**artial **O**rder **a**lignments) is a computational pipeline for calling consensi on R2C2 nanopore data.

This version of C3POa changes a lot. The old gonk branch can be found [here](https://github.com/rvolden/C3POa/tree/gonk). The even older version that uses water can be found [here](https://github.com/rvolden/C3POa/tree/water).

## Dependencies

- [Python 3](https://www.python.org/downloads/)
- [NumPy](https://pypi.org/project/numpy/)
- [SciPy](https://pypi.org/project/scipy/)
- [mappy](https://pypi.org/project/mappy/)
- [Cython](https://pypi.org/project/Cython/)
- [conk](https://github.com/rvolden/conk)
- [racon](https://github.com/isovic/racon)
- [editdistance](https://github.com/roy-ht/editdistance)
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)

To fetch and build dependencies, use setup.sh.

setup.sh will download and make the packages that you need to run C3POa (except for blat).

The setup script will install racon and conk. It will also download a "most linux"  blat executable into your C3POa directory.

```bash
chmod +x setup.sh
./setup.sh
```

If any of this doesn't work, you can download and install racon, conk, and blat manually and change where C3POa looks for them at the top of the C3POa.py and C3POa_postprocessing.py scripts.
Blat can built from [source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or you can get an [executable](http://hgdownload.soe.ucsc.edu/admin/exe/).
Please follow the documentation in the blat readme for make instructions.

--------------------------------------------------------------------------------

## Usage

After resolving all of the dependencies, you can run C3POa with python.

```bash
python3 C3POa.py -r path/to/reads.fastq -o path/to/where/C3POa/outputs/data/ -s splint.fasta -n 32
python3 C3POa_postprocessing.py -i path/to/where/C3POa/outputs/data/ -a adapter.fasta -x sampleSheet
```
Note that C3POa_postprocessing.py takes the output folder of C3POa.py as input

## C3POa.py

C3POa can be run on a nanopore 1D.fastq file or a directory containing multiple fastq files.

C3POa now has a --resume option that will look for a c3poa.log file in your output direcotry. 
It will look into the c3poa.log file in the output directory to determine if it was run previously on same of the fastqs in the input directory and then skip those.

If called on a directory, C3POa will now check the directory for new files coming in.
In essence, C3POa now does live consensus calling. 

You can start it at the same time as nanopore basecalling and just point it at the pass/ directory.
It will exit if no new fastq was created in the input directory for 30 minutes. 

Instead of pyabpoa, v3 of C3POa now uses the abpoa directly.
This allows us to use minimizer based msa for extra long inserts which reduces RAM requirements and prevents rare stalling events. 


```bash
python3 C3POa.py -r input/path 
                 -o output/path 
                 -s splint.fasta 
                 -n 32 

```

Arguments:
```
  --reads -r
                        FASTQ file that contains the long R2C2 reads or a directory containing multiple of these FASTQ files.
  --splint_file  -s
                        Path to the splint FASTA file.
  --out_path -o 
                        Directory where all the files will end up. Defaults to your current directory.
optional

  --lencutoff -l
                        Sets the length cutoff for your raw sequences. Anything shorter than the cutoff will be excluded. Defaults to 1000.
  --mdistcutoff -d
                        Sets the median distance cutoff for consensus sequences. Anything shorter will be excluded. Defaults to 500.
  --zero, -z            Use to exclude zero repeat reads. Defaults to True (includes zero repeats).
  --numThreads -n
                        Number of threads to use during multiprocessing. Defaults to 8.
  --groupSize -g
                        Number of reads processed by each thread in each iteration. Defaults to 100000.
  --blatThreads, -b     Use to chunk blat across the number of threads instead of by groupSize (faster).
  --compress_output, -co
                        Use to compress (gzip) both the consensus fasta and subread fastq output files.

  --resume, -u          If set, C3POa will look for c3poa.log file in output directory. If c3poa.log exists, files marked as processed in the file will be skipped.
                        Output will be appended to existing output files.
  --version, -v         Prints the C3POa version.


```

Example output read (readName_averageQuality_originalReadLength_numberOfRepeats_subreadLength):

```
>efbfbf09-7e2b-48e6-8e57-b3d36886739c_46.53_5798_2_1844
ACAGTCGATCATAGCTTAGCATGCATCGACGATCGATCGATCGA...
```

Example output directory tree:
```
output_dir
├── c3poa.log
├── tmp
│   └── splint_to_read_alignments.psl
├── Splint_1
│   ├── R2C2_Consensus.fasta
│   └── R2C2_Subreads.fastq
└── Splint_2
    ├── R2C2_Consensus.fasta
    └── R2C2_Subreads.fastq
```

--------------------------------------------------------------------------------

## C3POa_postprocessing.py

Trims and reorients consensus sequences generated by C3POa.py to 5'->3' direction.
If given a fasta of oligo dT indexes, it will also demux the reads by index.

```
-i  input_dir AND output_dir (has to be the output_dir used by C3POa.py)

-a  sequence of cDNA adapter sequences in fasta format. Sequence names must be
    '3Prime_adapter' and '5Prime_adapter' (or 'Adapter' if using -u flag)

optional:

-c  config file containing path to BLAT binary

-x  samplesheet

-n  number of threads to use

-u  use to ignore read directionality

-g  group size (number of reads given to each thread, default 1000)

-M  editdistance between read and best matching index in sample sheet has to be smaller than this number to return a match, default 2
    
-m  editdistance difference between read and best matching index and read and second best matching index has to be bigger than this number to return a match, default 1
    if editdistance between read and best matching index is 0, this cutoff does not apply. 

-bt split input by number of threads instead of groupSize

-co compress the output fasta/q files (gzip)

-v  print the C3POa version and exit
```

```bash
python3 C3POa_postprocessing.py -i /path/to/C3POa_directory  
                                -a /path/to/adapter.fasta
```

# sample sheet format

The sample sheet is a tab delimited file with one line per Splint/Index combination and what sample they are associated with.
(Multiple lines with different Splint/Index combinations can point to the same sample.)

The header line tells the script how to demultiplex so formatting that right is important. 
There are two types of indexes:

1) If your reads are directional, meaning you had a 5prime and 3prime adapter you can use 5prime and 3prime indexes which the program expects to find at the respective ends of the read. 
You can have an arbitrary number of those but we have only tested one 5 and one 3 index.
This is how your header line and first sample lines would look:

```bash
Name	Splint	5[0:16]	3[0:16]
Sample1	UMI_Splint1	ATATATAT	CGCGCGCG
Sample1	UMI_Splint3	CGATAGTG	ATCGAGTA
Sample2	UMI_Splint1	TGGTGGAT	TAGGACTA	
```

C3POa_postprocessing will look at the first 16 bases of the reads 5prime end for the first index and the first 16 bases of the 3prime end for the second index.
Include some extra bases around your index to have a little wiggle for the adapter matching.

2) If your reads are undirectional and you used the -u flag (for example regular Smartseq2 cDNA) but you used an index in either the TSO or oligodT primer, you can use an 'E' (for either) index.
You can only use one E index and they are not compatible with 5prime and 3prime indexes.
This is how your header line and first sample lines would look:

```bash
Name    Splint  E5[0:16]
Sample1 UMI_Splint1     ATATATAT
Sample1 UMI_Splint3     CGATAGTG
Sample2 UMI_Splint1     TGGTGGAT
```

The script will check both sides of the read for the index. 
The number (5 or 3) following the 'E' tells the script where the index should be (5' or 3') and the script will reorient the read based on where it finds the index so the read is 5'->3' in the output fasta.

Generally, the indexes you supply in the samplesheet should be in the orientation as they would appear in the PCR primers you amplified your library with.
If you run into any trouble here, just include all orientations of the indexes and see which ones get all the reads. 


Example output directory tree with demuxing:
```
output_dir
├── demultiplexed
│   ├── Sample1.fasta
│   ├── Sample2.fasta
│   └── Undetermined.fasta
├── Splint1
│   ├── R2C2_full_length_consensus_reads.fasta
│   ├── R2C2_full_length_consensus_reads_left_splint.fasta
│   └── R2C2_full_length_consensus_reads_right_splint.fasta
├── Splint2
│   ├── R2C2_full_length_consensus_reads.fasta
│   ├── R2C2_full_length_consensus_reads_left_splint.fasta
│   └── R2C2_full_length_consensus_reads_right_splint.fasta
└── no_index_found
    ├── R2C2_full_length_consensus_reads.fasta
    ├── R2C2_full_length_consensus_reads_left_splint.fasta
    └── R2C2_full_length_consensus_reads_right_splint.fasta
```
