#!/usr/bin/env python3
# Roger Volden and Chris Vollmers

import sys
import os
import argparse
import mappy as mm
import multiprocessing as mp
import editdistance as ld
from glob import glob
import gzip
import gc


VERSION = "v3  - It's over Anakin. I have the consensus!"

C3POaPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
blat=C3POaPath+'/blat/blat'

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Reorients/demuxes/trims consensus reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--input_folder', '-i', type=str, action='store',
                        help='input_dir AND output_dir (has to be the output_dir used by C3POa.py)')
    parser.add_argument('--adapter_file', '-a', type=str, action='store',
                        help='Fasta file with adapter (3 and 5 prime) sequences')
    parser.add_argument('--samplesheet', '-x', type=str, action='store',
                        help='samplesheet with header line indicating where to find indexes')
    parser.add_argument('--undirectional', '-u', action='store_true',
                        help='''By default, your cDNA molecules are assumed to be
                                directional with two sequences named "3Prime_adapter"
                                and "5Prime_adapter" expected in your adapter_file in
                                fasta format. If you add this flag your cDNA molecules
                                are expected to be undirectional and only one sequence
                                named "Adapter" should be in your adapter_file in fasta
                                format''')

    parser.add_argument('--threads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--groupSize', '-g', type=int, default=10000,
                        help='Number of reads processed by each thread in each iteration. Defaults to 10000.')
    parser.add_argument('--maxDist', '-M', type=int, default=2,
                        help='editdistance between read and best matching index in sample sheet has to be smaller than this number to return a match')
    parser.add_argument('--minDist', '-m', type=int, default=1,
                        help='editdistance difference between read and best matching index and read and second best matching index has to be bigger than this number to return a match')
    parser.add_argument('--blatThreads', '-bt', action='store_true', default=False,
                        help='''Use to chunk blat across the number of threads instead of by groupSize (faster).''')
    parser.add_argument('--skip_trimming', '-s', action='store_true', default=False,
                        help='''Use to demultiplex an already trimmed dataset''')
    parser.add_argument('--compress_output', '-co', action='store_true', default=False,
                        help='Use to compress (gzip) both the consensus fasta and subread fastq output files.')
    parser.add_argument('--version', '-v', action='version', version=VERSION, help='Prints the C3POa version.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def get_file_len(inFile):
    '''Figure out how many reads for best chunk size for parallelization'''
    count = 0
    for _ in mm.fastx_read(inFile, read_comment=False):
        count += 1
    return count

def remove_files(path, pattern):
    fileList=glob(path + pattern)
    length=len(fileList)
    fileCounter=0
    for d in fileList:
        fileCounter+=1
        print(f'\tRemoving tmp file {fileCounter} of {length} ({round((fileCounter/length)*100,2)}%)',' '*60,end='\r')
        os.system(f'rm -r {d}')
    print('\tFinished removing tmp files',' '*60)

def cat_files(path, pattern, output,compress,pos):
    '''Use glob to get around bash argument list limitations'''
    if compress:
        output += '.gz'
        final_fh = gzip.open(output, 'wb+')
    else:
        final_fh = open(output, 'w+')

    fileList=glob(path + pattern)
    length=len(fileList)
    fileCounter=0
    for f in fileList:
        fileCounter+=1
        print(f'\tCombining tmp file {fileCounter} of {length} ({round((fileCounter/length)*100,2)}%)',' '*60,end='\r')
        with open(f) as fh:
            for line in fh:
                if compress:
                    line=line.encode()
                final_fh.write(line)
    final_fh.close()
    print(f'\tFinished combining tmp files {pos}/4',' '*60,end='\r')

def process(args, reads, blat, iteration,subfolder):
    tmp_dir = subfolder + 'post_tmp_' + str(iteration) + '/'
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    tmp_fa = tmp_dir + 'tmp_for_blat.fasta'
    tmp_fa_fh = open(tmp_fa, 'w+')
    for header, seq in reads.items():
        print('>' + header, file=tmp_fa_fh)
        print(seq, file=tmp_fa_fh)
    tmp_fa_fh.close()

    run_blat(tmp_dir, tmp_fa, args.adapter_file, blat)
    os.remove(tmp_fa)
    adapter_dict = parse_blat(tmp_dir, reads)
    write_fasta_file(args, tmp_dir, adapter_dict, reads)

def chunk_process(input_fasta,subfolder, args, blat):
    '''Split the input fasta into chunks and process'''
    num_reads=get_file_len(input_fasta)

    if args.blatThreads:
        chunk_size = (num_reads // args.threads) + 1
    else:
        chunk_size = args.groupSize

    total=num_reads // chunk_size + 1
    if chunk_size > num_reads:
        chunk_size = num_reads
        total=1


    pool = mp.Pool(args.threads)
    print('\tAligning with adapters to read with BLAT and processing alignments')
    iteration, current_num, tmp_reads, target = 1, 0, {}, chunk_size
    for read in mm.fastx_read(input_fasta, read_comment=False):
        tmp_reads[read[0]] = read[1]
        current_num += 1
        if current_num == target:
            pool.apply_async(
                process,
                args=(args, tmp_reads, blat, iteration, subfolder),
                callback=print(f'\tfinished BLAT process {iteration} of {total}',' '*60,end='\r')
            )
            iteration += 1
            target = chunk_size * iteration
            if target >= num_reads:
                target = num_reads
            tmp_reads = {}
    print('')
    pool.close()
    pool.join()

    flc = 'R2C2_full_length_consensus_reads.fasta'
    flc_a = 'R2C2_full_length_consensus_reads.adapters.fasta'
    flc_left = 'R2C2_full_length_consensus_reads_left_splint.fasta'
    flc_right = 'R2C2_full_length_consensus_reads_right_splint.fasta'
    pattern = 'post_tmp*/'
    cat_files(subfolder,pattern + flc_left,subfolder + '/' + flc_left,args.compress_output,1)
    cat_files(subfolder,pattern + flc_right,subfolder + '/' + flc_right,args.compress_output,2)
    cat_files(subfolder,pattern + flc_a,subfolder + '/' + flc_a,args.compress_output,3)
    cat_files(subfolder,pattern + flc,subfolder + '/' + flc,args.compress_output,4)
    print('')
    remove_files(subfolder, 'post_tmp*')

def read_fasta(inFile, indexes):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict, index_dict = {}, {}
    for read in mm.fastx_read(inFile, read_comment=False):
        readDict[read[0]] = read[1]
        if indexes:
            index_dict[read[1]] = read[0]
    if indexes:
        return readDict, index_dict
    return readDict

def run_blat(path, infile, adapter_fasta, blat):
    align_psl = path + 'adapter_to_consensus_alignment.psl'
#    if not os.path.exists(align_psl) or os.stat(align_psl).st_size == 0:
    os.system('{blat} -noHead -stepSize=1 -tileSize=6 -t=DNA -q=DNA -minScore=10 \
                  -minIdentity=10 -minMatch=1 -oneOff=1 {adapters} {reads} {psl} >{blat_msgs}'
                  .format(blat=blat, adapters=adapter_fasta, reads=infile, psl=align_psl, blat_msgs=path + 'blat_msgs.log'))
#    else:
#        print('Reading existing psl file', file=sys.stderr)

def parse_blat(path, reads):
    adapter_dict, iterator = {}, 0

    for name, sequence in reads.items():
        adapter_dict[name] = {}
        adapter_dict[name]['+'] = []
        adapter_dict[name]['-'] = []
        adapter_dict[name]['+'].append(('-', 1, 0))
        adapter_dict[name]['-'].append(('-', 1, len(sequence)))

    with open(path + 'adapter_to_consensus_alignment.psl') as f:
        for line in f:
            a = line.strip().split('\t')
            read_name, adapter, strand = a[9], a[13], a[8]
            if int(a[5]) < 50 and float(a[0]) > 10:
                if strand == '+':
                    start = int(a[11]) - int(a[15])
                    end = int(a[12]) + (int(a[14]) - int(a[16]))
                    position = end
                if strand == '-':
                    start = int(a[11]) - (int(a[14]) - int(a[16]))
                    end = int(a[12]) + int(a[15])
                    position = start
                adapter_dict[read_name][strand].append((adapter,
                                                        float(a[0]),
                                                        position))
    return adapter_dict

def match_index(seq, seq_to_idx,minDist,maxDist):
    dist_dict, dist_list = {}, []
    # there needs to be a better/more efficient way to do this.
    for idx_seq, idx in seq_to_idx.items():
        idx=tuple(sorted(list(idx)))
        if idx not in dist_dict:
            dist_dict[idx] = []
        query = seq
        dist = ld.eval(query, idx_seq)
        dist_dict[idx].append(dist)
    for idx, distances in dist_dict.items():
        dist_list.append((idx, min(distances)))
    dist_list = sorted(dist_list, key=lambda x: x[1])
    match = tuple()
    if dist_list:
        if dist_list[0][1]==0 and dist_list[1][1]!=0:
            match = dist_list[0][0]
        elif dist_list[0][1] < maxDist:
            if len(dist_list)>1:
                if dist_list[1][1] - dist_list[0][1] > minDist:
                    match = dist_list[0][0]
            else:
                match = dist_list[0][0]

    return match

def write_fasta_file(args, path, adapter_dict, reads):
    undirectional = args.undirectional

    out = open(path + 'R2C2_full_length_consensus_reads.fasta', 'w')
    outa = open(path + 'R2C2_full_length_consensus_reads.adapters.fasta', 'w')
    out3 = open(path + 'R2C2_full_length_consensus_reads_left_splint.fasta', 'w')
    out5 = open(path + 'R2C2_full_length_consensus_reads_right_splint.fasta', 'w')

    length=len(reads)

    count=0
    for name, sequence in reads.items():
        count+=1
        adapter_plus = sorted(adapter_dict[name]['+'],
                              key=lambda x: x[2], reverse=False)
        adapter_minus = sorted(adapter_dict[name]['-'],
                              key=lambda x: x[2], reverse=False)
        plus_list_name, plus_positions = [], []
        minus_list_name, minus_positions = [], []

        for adapter in adapter_plus:
            if adapter[0] != '-':
                plus_list_name.append(adapter[0])
                plus_positions.append(adapter[2])
        for adapter in adapter_minus:
            if adapter[0] != '-':
                minus_list_name.append(adapter[0])
                minus_positions.append(adapter[2])

        if len(plus_list_name) != 1 or len(minus_list_name) != 1:
            continue
        if minus_positions[0] <= plus_positions[0]:
            continue

        if undirectional:
            direction = '+'
        elif plus_list_name[0] != minus_list_name[0]:
            if plus_list_name[0] == '5Prime_adapter':
                direction = '+'
            else:
                direction = '-'
        else:
            continue


        seq = sequence[plus_positions[0]:minus_positions[0]]

        ada1 = sequence[max(0,plus_positions[0]-40):plus_positions[0]]
        ada2 = sequence[minus_positions[0]:minus_positions[0]+40]

        name += '_' + str(len(seq))
        if direction == '+':
            out.write('>%s\n%s\n' %(name, seq))
            outa.write('>%s\t%s\t%s\n' %(name, ada1, ada2))
            out5.write('>%s\n%s\n' %(name, mm.revcomp(sequence[:plus_positions[0]])))
            out3.write('>%s\n%s\n' %(name, sequence[minus_positions[0]:]))

        elif direction == '-':
            out.write('>%s\n%s\n' %(name, mm.revcomp(seq)))
            outa.write('>%s\t%s\t%s\n' %(name, mm.revcomp(ada2),mm.revcomp(ada1)))
            out3.write('>%s\n%s\n' %(name, mm.revcomp(sequence[:plus_positions[0]+40])))
            out5.write('>%s\n%s\n' %(name, sequence[minus_positions[0]:]))


    out.close()
    outa.close()
    out3.close()
    out5.close()


def readSamplesheet(readFolder,sampleSheet):

    countDict={}
    countDict['All']=0
    countDict['Undetermined']=0

    if os.path.exists(readFolder+'/demultiplexed'):
        os.system('rm -r %s/demultiplexed' %(readFolder))
    os.system('mkdir %s/demultiplexed' %(readFolder))
    indexDict={}
    lineCounter=0
    SplintOnly=False
    outDict={}
    outDict['Undetermined']=readFolder+'/demultiplexed/Undetermined.fasta'
    for line in open(sampleSheet):
        lineCounter+=1
        if lineCounter==1:
            categories=line.strip().split('\t')
            if categories[0]=='Name' and categories[1]=='Splint':
                indexes=''
                if len(categories)>2:
                    indexes=categories[2:]
                    for index in indexes:
                        if 'E' in index and len(indexes)>1:
                            print('\t\tsamplesheet is not formatted properly: if using "E" index, it has to be the only index')
            else:
                print('\t\tsamplesheet is not formatted properly: Needs columns named "Name" and "Splint"')
                sys.exit(1)

        else:
            a=line.strip().split('\t')
            libraryName=a[0]
            outDict[libraryName]=readFolder+'/demultiplexed/'+libraryName+'.fasta'
            countDict[libraryName]=0
            Splint=a[1]
            if Splint not in indexDict:
                indexDict[Splint]={}
            if indexes:
                sequences=('_').join(a[2:])
                indexDict[Splint][sequences]=libraryName
            else:
                SplintOnly=True
                indexDict[Splint]=libraryName

    for libraryName,filePath in outDict.items():
        outTemp=open(filePath,'w')
        outTemp.close()
    return indexDict, indexes, SplintOnly, countDict,outDict

# print(indexDict)


def findIndexSequence(sequence,pattern):

    IUPACdict={}
    IUPACdict['A']=set(['A'])
    IUPACdict['T']=set(['T'])
    IUPACdict['G']=set(['G'])
    IUPACdict['C']=set(['C'])
    IUPACdict['R']=set(['A','G'])
    IUPACdict['Y']=set(['C','T'])
    IUPACdict['S']=set(['G','C'])
    IUPACdict['W']=set(['A','T'])
    IUPACdict['K']=set(['G','T'])
    IUPACdict['M']=set(['A','C'])
    IUPACdict['B']=set(['C','G','T'])
    IUPACdict['D']=set(['A','G','T'])
    IUPACdict['H']=set(['A','C','T'])
    IUPACdict['V']=set(['A','C','G'])
    IUPACdict['N']=set(['A','T','G','C'])
    UMI=''
    valid=True
    direction,range1,left,variable,right = pattern.split('.')
    if direction not in ['5','3','E']:
        print('invalid pattern, direction has to be 5 or 3')
        valid=False
        reason='invalid direction'
    if direction=='3':
        sequence=mm.revcomp(sequence)


    start,end = int(range1.split(':')[0]),int(range1.split(':')[1])
    left_variable_start=''
    right_start=''
    if len(sequence)<100:
        valid=False
        reason='sequence shorter than 100nt'

    if valid:
        UMIpattern=left+variable+right
        valid=False
        reason='no UMI pattern match'



        for pos in range(start,end,1):
            matches=0
            match=sequence[pos:pos+len(UMIpattern)].upper()
            for index in range(0,len(UMIpattern),1):
                v_base=UMIpattern[index]
                s_base=match[index]
                if s_base in IUPACdict[v_base]:
                    matches+=1
            if len(UMIpattern)==matches:
                if right:
                    UMI=match[len(left):-len(right)]
                else:
                    UMI=match[len(left):]
                valid=True
                break

    if valid:
        return UMI,'UMI found that matched pattern '+pattern
    else:
        return '', reason





counter=0
def demultiplex(seq,seq_to_idx,minDist,maxDist,libraryName,number,total):
    if not libraryName:
        matchSet=[]
        for index,entries in seq_to_idx.items():
            if index[0]=='E':
                readSeq, reason = findIndexSequence(seq,'5'+index[2:])
                if readSeq:
                    matchSet.append(match_index(readSeq,entries,minDist,maxDist))
                else:
                    readSeq, reason = findIndexSequence(seq,'3'+index[2:])
                    if readSeq:
                        matchSet.append(match_index(readSeq,entries,minDist,maxDist))
                        if index[1]=='5':
                            seq=mm.revcomp(seq)
                    else:
                        matchSet.append('Undetermined')
            else:
                readSeq, reason = findIndexSequence(seq,index)
                matchSet.append(match_index(readSeq,entries,minDist,maxDist))

        if 'Undetermined' in matchSet:
            libraryName = 'Undetermined'
        elif len(matchSet)==1:
            if len(matchSet[0])==1:
                libraryName = matchSet[0][0]
            else:
                libraryName = 'Undetermined'
        else:
            root=set(matchSet[0])
            for matches in matchSet[1:]:
                root=root & set(matches)
            if len(root)==1:
                libraryName = list(root)[0]
            else:
                libraryName = 'Undetermined'
    print(f'\tfinished read {number} of {total} reads total ~{str(round((number/total)*100,2))}%',' '*20, end='\r')
    return libraryName,seq


def main(args):

    blat = 'blat'

    print(f'C3POa postprocessing {VERSION}\nFinding and trimming adapters; optional samplesheet based demultiplexing')


    input_folder=args.input_folder

    skip_trimming=args.skip_trimming
    minDist=args.minDist
    maxDist=args.maxDist


    if not skip_trimming:
        print(f'\n\nFinding and Trimming adapters in directory {input_folder}\n\n')

        for folder in os.listdir(input_folder):
            subfolder=input_folder+'/'+folder
            if os.path.isdir(subfolder):
                for file in os.listdir(subfolder):
                    if 'R2C2_Consensus.fasta' in file:
                        print('Finding and trimming adapters in file', subfolder, file)
                        input_fasta=subfolder+'/'+file
                        chunk_process(input_fasta,subfolder, args, blat)
    else:
        print(f'\n\nskipping the trimming step for {input_folder}\n\n')


    if args.samplesheet:
        print('\n\nStarting to demultiplex \n\n')
        print('Reading sample sheet')
        indexDict, indexes, SplintOnly, countDict,outDict = readSamplesheet(input_folder,args.samplesheet)
        for folder in os.listdir(input_folder):
            counter=0
            subfolder=input_folder+'/'+folder
            if os.path.isdir(subfolder):
                if folder in indexDict:
                    seq_to_idx={}
                    for sequence,name in indexDict[folder].items():
                        sequences = sequence.split('_')
                        for i in range(0,len(sequences),1):
                            if indexes[i] not in seq_to_idx:
                                seq_to_idx[indexes[i]]={}
                            if sequences[i] not in seq_to_idx[indexes[i]]:
                                seq_to_idx[indexes[i]][sequences[i]]=set()
                            seq_to_idx[indexes[i]][sequences[i]].add(name)

                    readFile=subfolder+'/R2C2_full_length_consensus_reads.fasta'
                    print('Demultiplexing file '+readFile)
                    if args.compress_output:
                        readFile+='.gz'

                    print('\tdetermining number of reads')
                    total_reads=0
                    for read in mm.fastx_read(readFile, read_comment=False):
                        total_reads+=1
                    print(f'\t{total_reads} reads to demultiplex')

                    demuxGroupsize=100000
                    target=min(demuxGroupsize,total_reads)
                    current_num=0
                    threads=args.threads
                    results={}
                    tmp_reads=[]
                    iteration=1
                    for read in mm.fastx_read(readFile, read_comment=False):
                        tmp_reads.append(read)
                        current_num += 1
                        if current_num == target:
                            pool = mp.Pool(args.threads)
                            length_tmp_reads=len(tmp_reads)
                            for index in range(length_tmp_reads):
                                tmp_read=tmp_reads[index]
                                name,seq = tmp_read[0],tmp_read[1]
                                libraryName=''
                                if len(seq)<30:
                                    libraryName = 'Undetermined'
                                elif SplintOnly:
                                    libraryName = indexDict[folder]
                                results[name]=pool.apply_async(demultiplex,[seq,seq_to_idx,minDist,maxDist,libraryName,index+current_num-length_tmp_reads,total_reads])
                            pool.close()
                            pool.join()
                            gc.collect()
                            for name in results:
                                libraryName,seq=results[name].get()
                                fh=open(outDict[libraryName],'a')
                                fh.write('>%s\n%s\n' %(name,seq))
                                fh.close()
                            results={}
                            iteration += 1
                            target = demuxGroupsize * iteration
                            if target >= total_reads:
                                target = total_reads
                            tmp_reads = []



        print('\nFinished demultiplexing')


def batchProcess(reads,seq_to_idx,minDist,maxDist,SplintOnly,batch):
    resultList=[]
    for name,seq in tqdm(reads,position=batch,ncols=100):
          if len(seq)<30:
               libraryName = 'Undetermined'
          elif SplintOnly:
               libraryName = indexDict[folder]
          else:
               libraryName, seq = demultiplex(seq,seq_to_idx,minDist,maxDist)
          resultList.append((name,seq,libraryName))
    return (resultList)


if __name__ == '__main__':
    args = parse_args()
    if not args.input_folder or not args.adapter_file:
        print('Reads (--input_folder/-i) and adapter (--adapter_file/-a) are required', file=sys.stderr)
        sys.exit(1)
    main(args)
