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
- [pyabpoa](https://pypi.org/project/pyabpoa/)
- [mappy](https://pypi.org/project/mappy/)
- [tqdm](https://pypi.org/project/tqdm/)
- [Cython](https://pypi.org/project/Cython/)
- [conk](https://github.com/rvolden/conk)
- [racon](https://github.com/isovic/racon)
- [editdistance](https://github.com/roy-ht/editdistance)
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)

To fetch and build dependencies, use setup.sh.
setup.sh will download and make the packages that you need to run C3POa (except for blat).

The setup script **does not** install programs so you still need to put blat, and racon into your PATH.

```bash
chmod +x setup.sh
./setup.sh
```

Alternatively, you can grab all of the pip installable packages:
```bash
python3 -m pip install --user --upgrade scipy numpy pyabpoa==1.4.0 mappy Cython tqdm setuptools wheel
```
and then build conk and racon manually.

Blat can built from [source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or you can get an [executable](http://hgdownload.soe.ucsc.edu/admin/exe/).
Please follow the documentation in the blat readme for make instructions.

--------------------------------------------------------------------------------

## Usage

After resolving all of the dependencies, you can run C3POa with python.

```bash
python3 C3POa.py -r reads.fastq -o path/to/where/C3POa/outputs/data/ -s splint.fasta -n 32
python3 C3POa_postprocessing.py -i path/to/where/C3POa/outputs/data/ -a adapter.fasta -x sampleSheet
```

Note that C3POa_postprocessing.py takes the output folder of C3POa.py as input

## C3POa.py

Preprocessing is now built in.
Preprocessing takes raw 1D nanopore R2C2 reads in fastq (can be zipped) format, removes low quality and short reads and then finds splint sequences in those reads using BLAT.
Preprocessing will also demultiplex reads based on splints that are put into the splint fasta file.
The preprocessor will also look for the alignment psl file in case it was done before.
C3POa won't do the alignment if it finds `output_dir/tmp/splint_to_read_alignments.psl`.
By default, the input file will be chunked into fasta files of `len = group size`.
Use the `-b` option to have chunks be `len = number of reads / number of threads`.

The main algorithmic difference is we only align the splint to the read.
Command line tools have been replaced with their python APIs (except blat and racon).

```bash
python3 C3POa.py -r reads.fastq 
                 -o output/path 
                 -s splint.fasta 
                 -n 32 

```

Arguments:
```
  --reads -r
                        FASTQ file that contains the long R2C2 reads or a folder containing multiple of these FASTQ files.
  --splint_file  -s
                        Path to the splint FASTA file.
  --out_path -o 
                        Directory where all the files will end up. Defaults to your current directory.
optional

  --config
                        If you want to use a config file to specify paths to programs, specify them here. Use for racon and blat if they are not in your
                        path.
  --lencutoff -l
                        Sets the length cutoff for your raw sequences. Anything shorter than the cutoff will be excluded. Defaults to 1000.
  --mdistcutoff -d
                        Sets the median distance cutoff for consensus sequences. Anything shorter will be excluded. Defaults to 500.
  --zero, -z            Use to exclude zero repeat reads. Defaults to True (includes zero repeats).
  --numThreads -n
                        Number of threads to use during multiprocessing. Defaults to 1.
  --groupSize -g
                        Number of reads processed by each thread in each iteration. Defaults to 100000.
  --blatThreads, -b     Use to chunk blat across the number of threads instead of by groupSize (faster).
  --compress_output, -co
                        Use to compress (gzip) both the consensus fasta and subread fastq output files.
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
