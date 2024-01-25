import os
import sys
import numpy as np
import argparse
import multiprocessing as mp
import mappy as mm
from conk import conk
import gc
import gzip
import time

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

C3POaPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
abpoa=C3POaPath+'/abPOA-v1.4.1/bin/abpoa'
racon=C3POaPath+'/racon/build/bin/racon'
blat=C3POaPath+'/blat/blat'

from preprocess import preprocess
from call_peaks import call_peaks
from determine_consensus import determine_consensus

def rounding(x, base):
    '''Rounds to the nearest base, we use 50'''
    return int(base * round(float(x) / base))


def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Makes consensus sequences from R2C2 reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--reads', '-r', type=str, action='store',
                          help='FASTQ file that contains the long R2C2 reads or a folder containing multiple of these FASTQ files.')
    parser.add_argument('--splint_file', '-s', type=str, action='store',
                          help='Path to the splint FASTA file.')
    parser.add_argument('--out_path', '-o', type=str, action='store', default=os.getcwd(),
                        help='''Directory where all the files will end up.
                                Defaults to your current directory.''')
    parser.add_argument('--lencutoff', '-l', type=int, action='store', default=1000,
                        help='''Sets the length cutoff for your raw sequences. Anything
                                shorter than the cutoff will be excluded. Defaults to 1000.''')
    parser.add_argument('--mdistcutoff', '-d', type=int, action='store', default=500,
                        help='''Sets the median distance cutoff for consensus sequences.
                                Anything shorter will be excluded. Defaults to 500.''')
    parser.add_argument('--nosplint', '-ns', action='store_true',
                        help='''When set the first 300 bases of each read are used as splint''')
    parser.add_argument('--zero', '-z', action='store_false', default=True,
                        help='Use to exclude zero repeat reads. Defaults to True (includes zero repeats).')
    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')
    parser.add_argument('--peakFinderSettings', '-p', action='store', default='20,3,41,2')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def analyze_reads(args, read, splint, read_adapter, racon, tmp_dir,abpoa):

    peakFinderSettings=args.peakFinderSettings.split(',')
    penalty, iters, window, order = int(peakFinderSettings[0]),int(peakFinderSettings[1]),int(peakFinderSettings[2]),int(peakFinderSettings[3])

    name, seq, qual = read[0], read[1], read[2]
    seq_len = len(seq)
    use=False
    read_consensus = ''
    subs=[]
    scores = conk.conk(splint, seq, penalty)
    start=time.time()
    peaks = call_peaks(scores, args.mdistcutoff, iters, window, order)
    if  list(peaks):
        peaks = list(peaks + len(splint) // 2)
        for i in range(len(peaks) - 1, -1, -1):
            if peaks[i] >= seq_len:
                del peaks[i]
        if peaks:

            if args.nosplint:
                if len(peaks)>1:
                    subread_lens = np.diff(peaks)
                    relative_std = np.std(subread_lens)/np.mean(subread_lens)
                    if relative_std < 0.1:
                        use=True
                else:
                    use=True
            else:
                use=True

        # check for outliers in subread length
    if use:
        subreads, qual_subreads, dangling_subreads, qual_dangling_subreads = [], [], [], []
        if len(peaks) > 1:
            subread_lens = np.diff(peaks)
            subread_lens = [rounding(x, 50) for x in subread_lens]
            median_subread_len = np.median(subread_lens)
            for i in range(len(subread_lens)):
                bounds = [peaks[i], peaks[i+1]]
                if median_subread_len*0.8 <= subread_lens[i] <= median_subread_len*1.2:
                    subreads.append(seq[bounds[0]:bounds[1]])
                    qual_subreads.append(qual[bounds[0]:bounds[1]])
            if peaks[0] > 100:
                dangling_subreads.append(seq[:peaks[0]])
                qual_dangling_subreads.append(qual[:peaks[0]])
            if seq_len - peaks[-1] > 100:
                dangling_subreads.append(seq[peaks[-1]:])
                qual_dangling_subreads.append(qual[peaks[-1]:])
        else:
            dangling_subreads.append(seq[:peaks[0]])
            qual_dangling_subreads.append(qual[:peaks[0]])
            dangling_subreads.append(seq[peaks[0]:])
            qual_dangling_subreads.append(qual[peaks[0]:])


        consensus, repeats, subs = determine_consensus(
                args, read, subreads, qual_subreads, dangling_subreads, qual_dangling_subreads,
                racon, tmp_dir,abpoa
            )
        if consensus:
            numbers=[]
            for Q in qual:
                numbers.append(ord(Q)-33)
            avg_qual=str(round(np.average(numbers),1))
            cons_len = len(consensus)
            read_consensus ='>'+name+'_'+str(avg_qual)+'_'+str(seq_len)+'_'+str(repeats)+'_'+str(cons_len)+'\n'+consensus+'\n'
    return read_consensus,read_adapter,subs,peaks



def create_files(adapter,args,outDict,outSubDict,outCountDict):
    outCountDict[adapter]=0
    if not os.path.exists(args.out_path + adapter):
        os.mkdir(args.out_path + adapter)
    outCountDict[adapter]=0
    if args.compress_output:
        writeMode='ab+'
        outDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Consensus.fasta.gz',writeMode)
        outSubDict[adapter]=gzip.open(args.out_path + adapter +'/R2C2_Subreads.fastq.gz',writeMode)
    else:
        writeMode='a'
        outDict[adapter]=open(args.out_path + adapter +'/R2C2_Consensus.fasta',writeMode)
        outSubDict[adapter]=open(args.out_path + adapter +'/R2C2_Subreads.fastq',writeMode)
    return outDict,outSubDict,outCountDict


def main(args):
    log_file = open(args.out_path + 'c3poa.log', 'a+')
    pool = mp.Pool(args.numThreads)
    if not args.nosplint:
        splint_dict = {}
        for splint in mm.fastx_read(args.splint_file, read_comment=False):
            splint_dict[splint[0]] = [splint[1]]
            splint_dict[splint[0]].append(mm.revcomp(splint[1]))

    tmp_dir = args.out_path + 'tmp/'
    reads=args.reads
    fileStart=time.time()

    outDict={}
    outSubDict={}
    outCountDict={}
    previous=set()

    if not args.nosplint:
        adapter_dict, adapter_set, no_splint = preprocess(blat, args, tmp_dir, reads)
        for adapter in adapter_set:
            outDict,outSubDict,outCountDict=create_files(adapter,args,outDict,outSubDict,outCountDict)
    else:
        adapter_dict={}
        adapter_set=set()
        adapter_set.add('noSplint')
        for name,seq,q in mm.fastx_read(reads):
            adapter_dict[name]=seq[200:400]
        outDict,outSubDict,outCountDict=create_files('noSplint',args,outDict,outSubDict,outCountDict)
    total_reads=0
    short_reads=0
    no_splint_reads=0
    consNumber=0
    results={}
    processed_file=open(f'{reads}_processed','w')
    for name,seq,q in mm.fastx_read(reads, read_comment=False):
        total_reads+=1
        if len(seq) < args.lencutoff:
             short_reads+=1
             processed_file.write(f'{name}\n')
        elif name not in adapter_dict:
             no_splint_reads+=1
             processed_file.write(f'{name}\n')
        else:
            if not args.nosplint:
                read_adapter_info=adapter_dict[name]
                strand = read_adapter_info[1]
                if strand == '-':
                    # use reverse complement of the splint
                    splint = splint_dict[read_adapter_info[0]][1]
                else:
                    splint = splint_dict[read_adapter_info[0]][0]
                results[name]=pool.apply_async(analyze_reads,[args, [name,seq,q], splint, read_adapter_info[0], racon, tmp_dir,abpoa])
            else:
                 splint=adapter_dict[name]
                 results[name]=pool.apply_async(analyze_reads,[args, [name,seq,q], splint, 'noSplint', racon, tmp_dir,abpoa])

    pool.close()
    pool.join()
    gc.collect()


    adapters_with_reads=set()
    for index,result in results.items():
        consensus,adapter,subs,peaks = result.get()
        if consensus:
            if args.compress_output:
                consensus=consensus.encode()
            outDict[adapter].write(consensus)
            outCountDict[adapter]+=1
            consNumber+=1
            for subname,subseq,subq in subs:
                entry=f'@{subname}\n{subseq}\n+\n{subq}\n'
                if args.compress_output:
                    entry=entry.encode()
                outSubDict[adapter].write(entry)
        processed_file.write(f'{index}\n')
    processed_file.close()
    log_file.write(f'Too short reads: {short_reads}'+' ({:.2f}%)'.format((short_reads/total_reads)*100)+'\n')
    log_file.write(f'No splint reads: {no_splint_reads}'+' ({:.2f}%)'.format((no_splint_reads/total_reads)*100)+'\n')
    log_file.write(f'Successful consensus reads: {consNumber}'+' ({:.2f}%)'.format((consNumber/total_reads)*100)+'\n')

    fileEnd=time.time()
    duration=fileEnd-fileStart
    print(f'\tFinished generating {"{:,}".format(consNumber)} consensus sequences from {"{:,}".format(total_reads)} raw reads ({round((consNumber/total_reads)*100)}%) in {round(duration/60,1)} minutes.')
    for adapter in adapter_set:
        outDict[adapter].close()
        outSubDict[adapter].close()
        log_file.write(f'\t{outCountDict[adapter]} consensus reads generated for {adapter}\n')
    log_file.flush()
    log_file.close()


if __name__ == '__main__':
    args = parse_args()
    mp.set_start_method("spawn")
    main(args)

