import sys
import os
import argparse
import mappy as mm
import multiprocessing as mp
import editdistance as ld
from glob import glob
import gzip
import gc
import pickle


C3POaPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
blat=C3POaPath+'/blat/blat'

def parse_args():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Reorients/demuxes/trims consensus reads.',
                                     add_help=True,
                                     prefix_chars='-')
    parser.add_argument('--results_pickle', '-r', type=str, action='store',
                        help='results will be dumped here')
    parser.add_argument('--tmp_reads_pickle', '-s', type=str, action='store',
                        help='location of the pickled rreads to be demuxed')
    parser.add_argument('--seq_to_idx_pickle', '-i', type=str, action='store',
                        help='location of the pickled index dictionary')
    parser.add_argument('--threads', '-n', type=int, default=1,
                        help='Number of threads to use during multiprocessing. Defaults to 1.')
    parser.add_argument('--maxDist', '-M', type=int, default=2,
                        help='editdistance between read and best matching index in sample sheet has to be smaller than this number to return a match')
    parser.add_argument('--minDist', '-m', type=int, default=1,
                        help='editdistance difference between read and best matching index and read and second best matching index has to be bigger than this number to return a match')
    parser.add_argument('--total_reads', '-t', action='store', 
                        help='total number of reads to demultiplex')
    parser.add_argument('--current_num', '-c', action='store',
                        help='running number')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()



def demultiplex(seq,seq_to_idx,minDist,maxDist,libraryName,number,total):
    if not libraryName:
        matchSet=[]
        for index,entries in seq_to_idx.items():
            if index[0]=='E':
                readSeq, reason = findIndexSequence(seq,'5'+index[2:])
                if readSeq:
                    matchSet.append(match_index(readSeq,entries,minDist,maxDist))
                    if index[1]=='3':
                        seq=mm.revcomp(seq)

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


def main():
    args=parse_args()
    seq_to_idx=pickle.load(open(args.seq_to_idx_pickle,'rb'))
    tmp_reads=pickle.load(open(args.tmp_reads_pickle,'rb'))
    pool = mp.Pool(args.threads)
    length_tmp_reads=len(tmp_reads)
    results={}
    for index in range(length_tmp_reads):
        tmp_read=tmp_reads[index]
        name,seq = tmp_read[0],tmp_read[1]
        libraryName=''
        if len(seq)<30:
             libraryName = 'Undetermined'
        results[name]=pool.apply_async(demultiplex,[seq,seq_to_idx,int(args.minDist),int(args.maxDist),libraryName,index+int(args.current_num)-length_tmp_reads,int(args.total_reads)])
    pool.close()
    pool.join()
    gc.collect()
    new_results={}
    for name in results:
        libraryName,seq=results[name].get()
        new_results[name]=[libraryName,seq]
    pickle.dump(new_results,open(args.results_pickle,'wb'))


main()



