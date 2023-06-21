#!/usr/bin/env python3
# Roger Volden

import editdistance
import pyabpoa as poa
import mappy as mm
import os
import subprocess
from consensus import pairwise_consensus
import time
import numpy as np

def determine_consensus(args, read, subreads, sub_qual, dangling_subreads, qual_dangling_subreads, racon, tmp_dir,abpoa):
    start=time.time()
    name, seq, qual = read[0], read[1], read[2]
    repeats = len(subreads)
    subs=[]
    final_cons,subs,abpoa_cons = '', [],''
    poa_aligner = poa.msa_aligner(match=5)
    if repeats == 0:
        final_cons,subs = '', []
        if args.zero and len(dangling_subreads) == 2:
            tmp_final_cons, tmp_subs = zero_repeats(name, seq, qual, dangling_subreads, qual_dangling_subreads,subs,tmp_dir,abpoa)
            if tmp_final_cons and len(tmp_final_cons) >= args.mdistcutoff:
                final_cons = tmp_final_cons
                subs = tmp_subs
    if not final_cons:
        # overlap file is where the mappy alignment will go (req. by racon)
        overlap_file = tmp_dir + '{name}_overlaps.paf'.format(name=name)
        overlap_fh = open(overlap_file, 'w+')
        # temporary subreads specific for the current read (req. by racon)
        tmp_subread_file = tmp_dir + f'{name}_subreads.fastq'
        tmp_subread_fh = open(tmp_subread_file, 'w+')
        abpoa_msa=f'{tmp_dir}{name}.msa'
        abpoa_fasta=f'{tmp_dir}{name}.fasta'
        # align subreads together using abPOA
        tmp_subread_fl_file = tmp_dir + f'{name}_fl_subreads.fastq'
        tmp_subread_fl_fh = open(tmp_subread_fl_file, 'w+')
        insert_lengths=[]
        for i in range(len(subreads)):
            tmp_subread_fl_fh.write(f'@{i+1}\n{subreads[i]}\n+\n{sub_qual[i]}\n')
            insert_lengths.append(len(subreads[i]))

        tmp_subread_fl_fh.close()
        if repeats == 1:
            abpoa_cons = subreads[0]
        elif repeats == 2:
            insert_length=np.median(insert_lengths)
            if insert_length<8000:
                os.system(f'{abpoa} -M 5 -r 1 {tmp_subread_fl_file} > {abpoa_msa} 2> abpoa.messages')
            else:
                os.system(f'{abpoa} -M 5 -r 1 -S {tmp_subread_fl_file} > {abpoa_msa} 2> abpoa.messages')
            if not os.path.getsize(abpoa_msa)==0:
                msaSequences=[]
                for msaName,msaSeq,msaQ in mm.fastx_read(abpoa_msa):
                    msaSequences.append(msaSeq)
                abpoa_cons = pairwise_consensus(msaSequences, subreads, sub_qual,name,'pair')
        elif repeats > 2:
            insert_length=np.median(insert_lengths)
            if insert_length<8000:
                os.system(f'{abpoa} -M 5 -r 0 {tmp_subread_fl_file} > {abpoa_fasta} 2> abpoa.messages')
            else:
                os.system(f'{abpoa} -M 5 -r 0 -S {tmp_subread_fl_file} > {abpoa_fasta} 2> abpoa.messages')

#            res = poa_aligner.msa(subreads, out_cons=True, out_msa=True)
#            abpoa_cons1 = res.cons_seq[0]


            cons_seq=''
            for consName,consSeq,consQ in mm.fastx_read(f'{abpoa_fasta}'):
                cons_seq=consSeq
            if cons_seq:
                abpoa_cons = cons_seq

#            print(len(abpoa_cons),len(abpoa_cons1),editdistance.eval(abpoa_cons,abpoa_cons1))


    if abpoa_cons:

        # have to write out the consensus seq because it's going to get polished by racon
        abpoa_out = tmp_dir + f'{name}_abpoa.fasta'
        abpoa_out_fh = open(abpoa_out, 'w+')
        abpoa_out_fh.write(f'>{name}\n{abpoa_cons}\n')
        abpoa_out_fh.close()
        # map each of the subreads to the poa consensus
        mm_align = mm.Aligner(seq=abpoa_cons, preset='map-ont')
        for i in range(repeats):
            subread = subreads[i]
            q = sub_qual[i]
            qname = name + '_' + str(i+1)
            subs.append((qname,subread,q))
            tmp_subread_fh.write(f'@{qname}\n{subread}\n+\n{q}\n')
            for hit in mm_align.map(subread):
                overlap_fh.write(f'{qname}\t{len(subread)}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\t{name}\t{hit.ctg_len}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.mapq}\n')

        for j in range(len(dangling_subreads)):
            subread = dangling_subreads[j]
            q = qual_dangling_subreads[j]
            if j == 0:
                qname = name + '_' + str(j)
            else:
                qname = name + '_' + str(i+2)
            subs.append((qname,subread,q))
            tmp_subread_fh.write(f'@{qname}\n{subread}\n+\n{q}\n')
            for hit in mm_align.map(subread):
                overlap_fh.write(f'{qname}\t{len(subread)}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\t{name}\t{hit.ctg_len}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.mapq}\n')
        overlap_fh.close()
        tmp_subread_fh.close()

        racon_cons_file = tmp_dir + '{name}_racon_cons.fasta'.format(name=name)
        racon_cons_fh = open(racon_cons_file, 'w+')
        racon_msgs_fh = open(tmp_dir + 'racon_messages.log', 'w+')

    # polish poa cons with the subreads
        subprocess.run([racon, tmp_subread_file, overlap_file, abpoa_out, '-q', '5', '-t', '1'],
                       stdout=racon_cons_fh, stderr=racon_msgs_fh)
        racon_cons_fh.close()
        racon_msgs_fh.close()

        final_cons = ''
        for read in mm.fastx_read(racon_cons_file, read_comment=False):
            final_cons = read[1]

        os.system(f'rm {tmp_dir}{name}*')

    return final_cons, repeats, subs

def zero_repeats(name, seq, qual, subreads, sub_qual,subs,tmp_dir,abpoa):
    # subread is the master subread fastq for this group
    for i in range(len(subreads)):
        subs.append((name + '_' + str(i),subreads[i],sub_qual[i]))

    mappy_res = []
    mm_align = mm.Aligner(seq=subreads[0], preset='map-ont', scoring=(20, 7, 10, 5))
    for hit in mm_align.map(subreads[1]):
        mappy_res = [hit.r_st, hit.r_en, hit.q_st, hit.q_en]
    if not mappy_res:
        return '',subs

    left = subreads[1][:mappy_res[2]]
    right = subreads[0][mappy_res[1]:]
    overlap_seq1 = subreads[0][mappy_res[0]:mappy_res[1]]
    overlap_qual1 = sub_qual[0][mappy_res[0]:mappy_res[1]]
    overlap_seq2 = subreads[1][mappy_res[2]:mappy_res[3]]
    overlap_qual2 = sub_qual[1][mappy_res[2]:mappy_res[3]]

    tmp_file=tmp_dir + f'{name}_overlaps.fasta'
    tmp_file_fh=open(tmp_file,'w')
    tmp_file_fh.write(f'>1\n{overlap_seq1}\n>2\n{overlap_seq2}\n')
    tmp_file_fh.close()
    if len(overlap_seq1)<8000:
        os.system(f'{abpoa} -M 5 -r 1 {tmp_file} > {tmp_dir}{name}.pw.msa 2> apboa.messages')
    else:
        os.system(f'{abpoa} -M 5 -r 1 -S {tmp_file} > {tmp_dir}{name}.pw.msa 2> apboa.messages')
    if os.path.getsize(f'{tmp_dir}{name}.pw.msa')==0:
        return '',subs
    else:
        msaSequences=[]
        for msaName,msaSeq,msaQ in mm.fastx_read(f'{tmp_dir}{name}.pw.msa'):
            msaSequences.append(msaSeq)
        abpoa_cons = pairwise_consensus(msaSequences, [overlap_seq1, overlap_seq2], [overlap_qual1, overlap_qual2],name,'zero')
        corrected_cons = left + abpoa_cons + right
        return corrected_cons,subs
