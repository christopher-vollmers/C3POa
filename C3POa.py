#!/usr/bin/env python3
# Roger Volden

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


VERSION = "v3.2 - Bombad Consensus"

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
                        help='''When set the first bases 200-400 of each read are used as splint''')
    parser.add_argument('--zero', '-z', action='store_false', default=True,
                        help='Use to exclude zero repeat reads. Defaults to True (includes zero repeats).')
    parser.add_argument('--numThreads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')
    parser.add_argument('--peakFinderSettings', '-p', action='store', default='20,3,41,2',
                        help='Only set this if you have a really short splint (<50nt) and all your reads are discarded. Defaults to "20,3,41,2". Try "30,3,15,2" for a short splint. No promises though. We only tested C3POa for splints >100nt')
    parser.add_argument('--resume', '-u', action='store_true', default=False,
                        help='''If set, C3POa will look for processed.log file in output directory. 
                                If processed.log exists, reads marked as processed in the input will be skipped. 
                                Output will be appended to existing output files.''')
    parser.add_argument('--version', '-v', action='version', version=VERSION, help='Prints the C3POa version.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def getFileList(query_path,done):
    file_list=[]
    '''Takes a path and returns a list of fast5 files excluding the most recent fast5.'''
    for file in os.listdir(query_path):
        if 'fastq' in file and 'temp' not in file:
            if os.path.abspath(query_path+file) not in done:
                file_list.append(os.path.abspath(query_path+file))
    if file_list:
        exclude = max(file_list, key=os.path.getctime)
        if time.time()-os.path.getctime(exclude)<600:
            file_list.remove(exclude)
        file_list.sort(key=lambda x:os.path.getctime(x))
    return file_list





def main(args):
    argString=''
    argString+=f'-s {args.splint_file} '
    argString+=f'-o {args.out_path} '
    argString+=f'-l {args.lencutoff} '
    argString+=f'-d {args.mdistcutoff} '
    if args.nosplint:
        argString+=f'-ns '
    if args.zero:
        argString+=f'-z '
    argString+=f'-n {args.numThreads} '
    if args.compress_output:
        argString+=f'-co '
    argString+=f'-p {args.peakFinderSettings} '
    print(argString)

    print(f'C3POa {VERSION} \nGenerating consensus sequences from R2C2 read data')
    if not args.resume:
        print('\nRemoving old results\n')
        if not args.nosplint:
            splint_dict = {}
            for adapter,seq,q in mm.fastx_read(args.splint_file, read_comment=False):
                if os.path.isdir(f'{args.out_path}/{adapter}'):
                    os.system(f'rm -r {args.out_path}/{adapter}')
        else:
            if os.path.isdir(f'{args.out_path}/noSplint'):
                os.system(f'rm -r {args.out_path}/noSplint')

    if not args.out_path.endswith('/'):
        args.out_path += '/'
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)


    resume=args.resume
    done=set()
    processed_reads=set()
    if resume:
        print(f'\n--resume option is True: Looking for existing log file in {args.out_path}')
        if os.path.isfile(args.out_path + 'processed.log'):
            print('log file found')
            for line in open(args.out_path + 'processed.log'):
                processed_reads.add(line.strip())
        print(f'{len(processed_reads)} processed reads found in log file. They will be skipped\n')

    log_file = open(args.out_path + 'c3poa.log', 'a+')
    log_file.write(f'C3POa version: {VERSION}\n')
    iterate=True
    timeAtLastFile=time.time()
    timeSinceLastFile=0
    previous=set()
    consNumberTotal=0
    while iterate:
        fileTimes=[]
        fileStart=time.time()
        log_file.write('new iteration\n')
        tmp_dir = args.out_path + 'tmp/'
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        else:
            os.system(f'rm -r {tmp_dir}')
            os.mkdir(tmp_dir)

        print('Starting consensus calling iteration - if input is directory it will check for files that are new since last iteration')


        input_path=args.reads
        if os.path.isdir(input_path):
             print('\tRead input is a folder')
             if not input_path.endswith('/'):
                 input_path+='/'
             file_list=getFileList(input_path,done)

        elif os.path.isfile(input_path):
            print('\tRead input is a file')
            file_list=[]
            file_list.append(os.path.abspath(input_path))
            iterate=False
        else:
            print('\tno file provided')
            iterate=False
        print(f'\t{len(file_list)} file(s) provided')

        if len(file_list)==0:
            timeSinceLastFile=time.time()-timeAtLastFile
            print(f'\t{round(timeSinceLastFile/60,2)} minutes since last file was provided. Will terminate now if input was fastq file or after more than 30 minutes if input was directory')
            time.sleep(30)
            if timeSinceLastFile>1800:
                iterate=False
            continue
        else:
            timeAtLastFile=time.time()
            log_file.write(f'Total files to process: {len(file_list)}\n')
            tmp_file=open(f'{tmp_dir}/tmp_file','w')
            if args.resume:
                writeMode='a'
            else:
                writeMode='w'
            processed_file=open(f'{args.out_path}/processed.log',writeMode)
            for reads in file_list:
                total_reads=0
                print(f'\tProcessing file {reads}')
                log_file.write(f'Processing file {reads}\n')
                for name,seq,q in mm.fastx_read(reads, read_comment=False):
                    if name not in processed_reads:
                        total_reads+=1
                        tmp_file.write(f'@{name}\n{seq}\n+\n{q}\n')
                        if total_reads%10000==0:
                            if os.path.getsize(f'{tmp_dir}/tmp_file')>0:
                                os.system(f'python3 {C3POaPath}/generateConsensus.py -r {tmp_dir}/tmp_file {argString}')
                                for line in open(f'{tmp_dir}/tmp_file_processed'):
                                    processed_file.write(line)
                                tmp_file=open(f'{tmp_dir}/tmp_file','w')
                if os.path.getsize(f'{tmp_dir}/tmp_file')>0:
                    os.system(f'python3 {C3POaPath}/generateConsensus.py -r {tmp_dir}/tmp_file {argString}')
                    for line in open(f'{tmp_dir}/tmp_file_processed'):
                        processed_file.write(line)
                    tmp_file=open(f'{tmp_dir}/tmp_file','w')
                done.add(reads)

            print('\n')
    log_file.close()
    processed_file.close()
    print('\n')

if __name__ == '__main__':
    args = parse_args()
    mp.set_start_method("spawn")
    main(args)

