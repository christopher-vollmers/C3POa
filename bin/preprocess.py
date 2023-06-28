#!/usr/bin/env python3
# Roger Volden

import os
import sys
import mappy

def preprocess(blat, args, tmp_dir,fastq_file):
    tmp_fasta = tmp_dir + 'R2C2_temp_for_BLAT.fasta'
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'

    tmp_adapter_dict={}
#    print('\tAligning splints to reads with blat', file=sys.stderr)

    process(args, fastq_file, blat, tmp_dir)

    adapter_set = set()
    with open(align_psl) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            line = line.split('\t')
            read_name, adapter, strand = line[9], line[13], line[8]
            gaps, score = float(line[5]), float(line[0])
            if gaps < 50 and score > 50:
                if read_name not in tmp_adapter_dict:
                    tmp_adapter_dict[read_name]=[[None, 1, None]]
                tmp_adapter_dict[read_name].append([adapter, float(line[0]), strand])
                adapter_set.add(adapter)

    adapter_dict = {} # read_id: [adapter, strand]
    no_splint_reads = 0
    for name, alignments in tmp_adapter_dict.items():
        best = sorted(alignments, key=lambda x: x[1], reverse=True)[0]
        if not best[0]:
            no_splint_reads += 1
            continue
        adapter_set.add(best[0])
        adapter_dict[name] = [best[0], best[2]]
    return adapter_dict, adapter_set, no_splint_reads

def process(args, fastq_file, blat, tmp_dir):
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)
    tmp_fa = tmp_dir + 'tmp_for_blat.fasta'
    tmp_fa_fh = open(tmp_fa, 'w+')
    for header, seq,q in mappy.fastx_read(fastq_file):
        tmp_fa_fh.write(f'>{header}\n{seq}\n')
    tmp_fa_fh.close()
    align_psl = tmp_dir + 'splint_to_read_alignments.psl'
    b_msgs = tmp_dir + 'blat_messages.log'

    os.system('{blat} -noHead -stepSize=1 -t=DNA -q=DNA -minScore=15 \
              -minIdentity=10 {splint} {reads} {psl} >{blat_msgs}'
              .format(blat=blat, splint=args.splint_file, reads=tmp_fa, psl=align_psl, blat_msgs=b_msgs))
    os.remove(tmp_fa)

